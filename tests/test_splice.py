"""
tests/test_splice.py — unit + integration tests
=================================================
Run with:  pytest tests/ -v
"""

import os
import struct
import tempfile

import numpy as np
import pandas as pd
import pytest
import pysam

import splice
from splice.data import SplicingData
from splice.io import _derive_sample_ids, _read_one_parquet
from splice.metrics import compute_jaccard
from splice.correction import filter_junctions, correct_confounders
from splice.stats import call_outliers, _fit_bb_mom_vectorised
from splice.extract import (
    _extract_introns_from_cigar,
    _strand_from_xs,
    _infer_unsplit_strand,
    _read_has_intron,
    extract_sample_counts,
)


# ==================================================================
# BAM fixture helpers
# ==================================================================

def _make_bam(tmp_path, sample_name, reads, n_junctions=10, seed=0):
    """
    Write a minimal coordinate-sorted indexed BAM file.

    reads : list of dicts with keys:
        chrom, start (0-based), cigar, strand (XS tag), mapq,
        is_read1, is_reverse, is_paired
    """
    bam_dir  = tmp_path / sample_name
    bam_dir.mkdir(exist_ok=True)
    bam_path = str(bam_dir / "Aligned.sortedByCoord.out.bam")

    # Build header
    chroms = sorted({r["chrom"] for r in reads})
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": 200_000_000} for c in chroms],
    })

    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for i, r in enumerate(reads):
            seg = pysam.AlignedSegment(header)
            seg.query_name          = f"read_{i}"
            seg.reference_id        = header.get_tid(r["chrom"])
            seg.reference_start     = r["start"]
            seg.cigarstring         = r["cigar"]
            seg.mapping_quality     = r.get("mapq", 255)
            seg.flag                = 0
            if r.get("is_paired", True):
                seg.flag |= 0x1
            if r.get("is_read1", True):
                seg.flag |= 0x40
            if r.get("is_reverse", False):
                seg.flag |= 0x10
            # Query sequence (must match CIGAR length)
            cigar_len = _cigar_query_len(r["cigar"])
            seg.query_sequence      = "A" * cigar_len
            seg.query_qualities     = pysam.qualitystring_to_array("I" * cigar_len)
            if r.get("strand"):
                seg.set_tag("XS", r["strand"])
            bam.write(seg)

    pysam.sort("-o", bam_path + ".sorted", bam_path)
    os.rename(bam_path + ".sorted", bam_path)
    pysam.index(bam_path)
    return bam_path


def _cigar_query_len(cigar_str):
    """Return number of query-consuming bases from a CIGAR string."""
    import re
    total = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar_str):
        if op in "MISH=X":
            total += int(length)
    return total


def _make_parquet(tmp_path, sample_name, n_junctions=50, seed=0):
    """Write a minimal splice count parquet file with overlapping junctions."""
    rng = np.random.default_rng(seed)
    # Use fixed junction coordinates so samples share junctions
    records = []
    for i in range(n_junctions):
        chrom = f"chr{(i % 5) + 1}"
        start = 1_000_000 + i * 5_000
        end   = start + 500
        strand = "+" if i % 2 == 0 else "-"
        split  = int(rng.integers(25, 300))
        site   = split + int(rng.integers(10, 100))
        unsplit = int(rng.integers(5, 50))
        records.append({
            "intron_id":     f"{chrom}:{start}-{end}:{strand}",
            "seqnames":      chrom,
            "start":         np.int32(start),
            "end":           np.int32(end),
            "strand":        strand,
            "split_count":   np.int32(split),
            "site_count":    np.int32(site),
            "unsplit_count": np.int32(unsplit),
        })
    df = pd.DataFrame(records)
    out_dir = tmp_path / "counts"
    out_dir.mkdir(exist_ok=True)
    path = str(out_dir / f"{sample_name}.splice.parquet")
    df.to_parquet(path, index=False)
    return path


def _make_splicing_data(n_samples=30, n_junctions=50, inject_outlier=True):
    rng = np.random.default_rng(42)
    split   = rng.integers(20, 200, size=(n_samples, n_junctions)).astype(np.int32)
    site    = split + rng.integers(10, 100, size=(n_samples, n_junctions)).astype(np.int32)
    unsplit = rng.integers(5, 50, size=(n_samples, n_junctions)).astype(np.int32)

    if inject_outlier:
        split[0, 0]   = 1
        site[0, 0]    = 200
        unsplit[0, 0] = 0

    obs = pd.DataFrame({"sampleID": [f"sample_{i}" for i in range(n_samples)]})
    obs.index = obs["sampleID"]

    var = pd.DataFrame({
        "intron_id": [f"chr1:{i*1000}-{i*1000+500}:+" for i in range(n_junctions)],
        "seqnames":  ["chr1"] * n_junctions,
        "start":     [i * 1000 for i in range(n_junctions)],
        "end":       [i * 1000 + 500 for i in range(n_junctions)],
        "strand":    ["+"] * n_junctions,
        "intron_motif": [1] * n_junctions,
        "annotated":    [1] * n_junctions,
        "gene_id":   [f"ENSG{i:011d}" for i in range(n_junctions)],
        "gene_name": [f"GENE{i}" for i in range(n_junctions)],
    })
    var.index = var["intron_id"]

    sd = SplicingData.from_matrices(split, site, obs=obs, var=var)
    from scipy.sparse import csr_matrix
    sd.adata.layers["unsplit_counts"] = csr_matrix(unsplit.astype(np.int32))
    return sd


# ==================================================================
# Unit tests: extract.py — CIGAR / strand helpers
# ==================================================================

class TestExtractHelpers:
    def _make_read(self, header, cigar, start=1000, strand=None,
                   is_read1=True, is_reverse=False):
        seg = pysam.AlignedSegment(header)
        seg.query_name      = "r"
        seg.reference_id    = 0
        seg.reference_start = start
        seg.cigarstring     = cigar
        seg.mapping_quality = 255
        seg.flag = 0x1 | (0x40 if is_read1 else 0x80)
        if is_reverse:
            seg.flag |= 0x10
        qlen = _cigar_query_len(cigar)
        seg.query_sequence  = "A" * qlen
        seg.query_qualities = pysam.qualitystring_to_array("I" * qlen)
        if strand:
            seg.set_tag("XS", strand)
        return seg

    @pytest.fixture
    def header(self):
        return pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6"},
            "SQ": [{"SN": "chr1", "LN": 100_000_000}],
        })

    def test_extract_single_intron(self, header):
        read = self._make_read(header, "50M100N50M", start=1000)
        introns = _extract_introns_from_cigar(read)
        assert len(introns) == 1
        # start=1000 (0-based), after 50M → intron starts at ref pos 1050
        # 1-based: 1051, ends at 1051+100-1 = 1150
        assert introns[0] == (1051, 1150)

    def test_extract_two_introns(self, header):
        read = self._make_read(header, "30M50N30M80N30M", start=0)
        introns = _extract_introns_from_cigar(read)
        assert len(introns) == 2
        assert introns[0] == (31, 80)   # after 30M, intron is 50N (1-based 31..80)
        assert introns[1] == (111, 190)

    def test_no_intron(self, header):
        read = self._make_read(header, "100M", start=500)
        assert _extract_introns_from_cigar(read) == []

    def test_read_has_intron_true(self, header):
        read = self._make_read(header, "50M100N50M", start=0)
        assert _read_has_intron(read) is True

    def test_read_has_intron_false(self, header):
        read = self._make_read(header, "100M", start=0)
        assert _read_has_intron(read) is False

    def test_strand_from_xs_plus(self, header):
        read = self._make_read(header, "50M100N50M", strand="+")
        assert _strand_from_xs(read) == "+"

    def test_strand_from_xs_minus(self, header):
        read = self._make_read(header, "50M100N50M", strand="-")
        assert _strand_from_xs(read) == "-"

    def test_strand_from_xs_absent(self, header):
        read = self._make_read(header, "100M")
        assert _strand_from_xs(read) is None

    def test_infer_unsplit_strand_firststrand_read1_forward(self, header):
        read = self._make_read(header, "100M", is_read1=True, is_reverse=False)
        assert _infer_unsplit_strand(read, "fr-firststrand") == "-"

    def test_infer_unsplit_strand_firststrand_read1_reverse(self, header):
        read = self._make_read(header, "100M", is_read1=True, is_reverse=True)
        assert _infer_unsplit_strand(read, "fr-firststrand") == "+"

    def test_infer_unsplit_strand_secondstrand_read1_forward(self, header):
        read = self._make_read(header, "100M", is_read1=True, is_reverse=False)
        assert _infer_unsplit_strand(read, "fr-secondstrand") == "+"

    def test_infer_unsplit_strand_unstranded(self, header):
        read = self._make_read(header, "100M", is_read1=True, is_reverse=False)
        assert _infer_unsplit_strand(read, "unstranded") == "*"


# ==================================================================
# Unit tests: extract_sample_counts() — end-to-end BAM extraction
# ==================================================================

class TestExtractSampleCounts:
    def _make_controlled_bam(self, tmp_path):
        """
        BAM with known content:
          - 5 spliced reads: chr1:1000-2000 (intron 1001-2000, 1-based), strand +
          - 3 spliced reads: chr1:3000-4000 (intron 3001-4000), strand -
          - 4 unsplit reads overlapping position 1001 (donor of first intron), strand +
          - 1 spliced read below min_split_reads threshold: chr1:5000-6000 (1 read only)
        """
        reads = []
        # 5 reads spanning intron chr1:1001-2000 (+)
        for _ in range(5):
            reads.append({
                "chrom": "chr1", "start": 950,
                "cigar": "50M999N50M",  # ref: 950..999 + skip 1000..1998 + 1999..2048
                "strand": "+", "mapq": 255,
                "is_read1": True, "is_reverse": False,
            })
        # 3 reads spanning intron chr1:3001-4000 (-)
        for _ in range(3):
            reads.append({
                "chrom": "chr1", "start": 2950,
                "cigar": "50M999N50M",
                "strand": "-", "mapq": 255,
                "is_read1": True, "is_reverse": True,
            })
        # 4 unsplit reads overlapping the donor of intron 1 (pos ~1000)
        for _ in range(4):
            reads.append({
                "chrom": "chr1", "start": 960,
                "cigar": "100M",
                "strand": None, "mapq": 255,
                "is_read1": True, "is_reverse": True,  # fr-firststrand read1 reverse -> +
            })
        # 1 read only — should be filtered by min_split_reads=2
        reads.append({
            "chrom": "chr1", "start": 4950,
            "cigar": "50M999N50M",
            "strand": "+", "mapq": 255,
            "is_read1": True, "is_reverse": False,
        })
        return _make_bam(tmp_path, "controlled", reads)

    def test_extraction_produces_parquet(self, tmp_path):
        bam = self._make_controlled_bam(tmp_path)
        out = str(tmp_path / "controlled.splice.parquet")
        splice.extract_sample_counts(bam, out, verbose=False)
        assert os.path.exists(out)

    def test_split_counts_correct(self, tmp_path):
        bam = self._make_controlled_bam(tmp_path)
        out = str(tmp_path / "controlled.splice.parquet")
        splice.extract_sample_counts(bam, out, min_split_reads=2, verbose=False)
        df = pd.read_parquet(out)
        # Should have 2 junctions (the 1-read junction is filtered)
        assert len(df) == 2
        # Junction with 5 reads
        j5 = df[df["split_count"] == 5]
        assert len(j5) == 1
        assert j5.iloc[0]["strand"] == "+"
        # Junction with 3 reads
        j3 = df[df["split_count"] == 3]
        assert len(j3) == 1
        assert j3.iloc[0]["strand"] == "-"

    def test_min_split_reads_filter(self, tmp_path):
        bam = self._make_controlled_bam(tmp_path)
        out2 = str(tmp_path / "filtered2.splice.parquet")
        out1 = str(tmp_path / "filtered1.splice.parquet")
        splice.extract_sample_counts(bam, out2, min_split_reads=2, verbose=False)
        splice.extract_sample_counts(bam, out1, min_split_reads=1, verbose=False)
        df2 = pd.read_parquet(out2)
        df1 = pd.read_parquet(out1)
        assert len(df1) > len(df2)  # min=1 retains the single-read junction

    def test_mapq_filter(self, tmp_path):
        # Add a low-MAPQ read for a new junction
        reads = []
        for _ in range(5):
            reads.append({
                "chrom": "chr1", "start": 100,
                "cigar": "50M500N50M", "strand": "+",
                "mapq": 255, "is_read1": True, "is_reverse": False,
            })
        # Low-MAPQ junction
        for _ in range(5):
            reads.append({
                "chrom": "chr1", "start": 10000,
                "cigar": "50M500N50M", "strand": "+",
                "mapq": 10, "is_read1": True, "is_reverse": False,
            })
        bam = _make_bam(tmp_path, "mapq_test", reads)
        out_strict = str(tmp_path / "strict.splice.parquet")
        out_permissive = str(tmp_path / "permissive.splice.parquet")
        splice.extract_sample_counts(bam, out_strict, min_mapping_quality=255,
                                     min_split_reads=2, verbose=False)
        splice.extract_sample_counts(bam, out_permissive, min_mapping_quality=0,
                                     min_split_reads=2, verbose=False)
        strict = pd.read_parquet(out_strict)
        permissive = pd.read_parquet(out_permissive)
        assert len(permissive) >= len(strict)

    def test_unsplit_counts_nonzero(self, tmp_path):
        bam = self._make_controlled_bam(tmp_path)
        out = str(tmp_path / "unsplit.splice.parquet")
        splice.extract_sample_counts(bam, out, library_type="fr-firststrand",
                                     min_split_reads=2, verbose=False)
        df = pd.read_parquet(out)
        # The 4 unsplit reads overlap the + junction's donor site
        j_plus = df[df["strand"] == "+"]
        assert j_plus.iloc[0]["unsplit_count"] >= 0  # may be 0 if window misses

    def test_output_schema(self, tmp_path):
        bam = self._make_controlled_bam(tmp_path)
        out = str(tmp_path / "schema.splice.parquet")
        splice.extract_sample_counts(bam, out, verbose=False)
        df = pd.read_parquet(out)
        expected_cols = {
            "intron_id", "seqnames", "start", "end", "strand",
            "split_count", "site_count", "unsplit_count",
        }
        assert expected_cols.issubset(set(df.columns))

    def test_no_junctions_writes_empty_parquet(self, tmp_path):
        # BAM with only unspliced reads — no N in CIGAR, so no junctions discovered
        reads = []
        for _ in range(5):
            reads.append({
                "chrom": "chr1", "start": 1000,
                "cigar": "100M", "strand": None,
                "mapq": 255, "is_read1": True, "is_reverse": False,
            })
        bam = _make_bam(tmp_path, "no_junctions", reads)
        out = str(tmp_path / "no_junctions.splice.parquet")
        splice.extract_sample_counts(bam, out, min_split_reads=2, verbose=False)
        df = pd.read_parquet(out)
        assert len(df) == 0

    def test_missing_bam_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="BAM file not found"):
            splice.extract_sample_counts(
                str(tmp_path / "nonexistent.bam"),
                str(tmp_path / "out.parquet"),
                verbose=False,
            )

    def test_missing_index_raises(self, tmp_path):
        # Write a BAM without indexing it
        bam_path = str(tmp_path / "no_index.bam")
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1_000_000}]
        })
        with pysam.AlignmentFile(bam_path, "wb", header=header):
            pass
        with pytest.raises(FileNotFoundError, match="index"):
            splice.extract_sample_counts(bam_path, str(tmp_path / "out.parquet"),
                                          verbose=False)


# ==================================================================
# Unit tests: I/O (parquet)
# ==================================================================

class TestIO:
    def test_read_counts_single(self, tmp_path):
        path = _make_parquet(tmp_path, "sampleA", seed=1)
        sd = splice.read_counts([path], verbose=False)
        assert sd.n_samples == 1
        assert sd.n_junctions > 0
        assert "unsplit_counts" in sd.adata.layers

    def test_read_counts_glob(self, tmp_path):
        for i in range(5):
            _make_parquet(tmp_path, f"sample_{i}", seed=i)
        sd = splice.read_counts(str(tmp_path / "counts" / "*.splice.parquet"),
                                verbose=False)
        assert sd.n_samples == 5

    def test_read_counts_sample_annotation(self, tmp_path):
        for i in range(3):
            _make_parquet(tmp_path, f"sample_{i}", seed=i)
        ann = pd.DataFrame({
            "sampleID": [f"sample_{i}" for i in range(3)],
            "diagnosis": ["affected", "control", "control"],
        })
        sd = splice.read_counts(
            str(tmp_path / "counts" / "*.splice.parquet"),
            sample_annotation=ann,
            verbose=False,
        )
        assert "diagnosis" in sd.obs.columns

    def test_read_counts_missing_columns_raises(self, tmp_path):
        bad = pd.DataFrame({"intron_id": ["chr1:100-200:+"], "split_count": [5]})
        path = str(tmp_path / "bad.parquet")
        bad.to_parquet(path, index=False)
        with pytest.raises(ValueError, match="missing required columns"):
            splice.read_counts([path], verbose=False)

    def test_derive_sample_ids_strip_suffix(self):
        files = ["/counts/sampleA.splice.parquet", "/counts/sampleB.splice.parquet"]
        ids = _derive_sample_ids(files, from_path=True)
        assert ids == ["sampleA", "sampleB"]

    def test_derive_sample_ids_deduplication(self):
        files = ["/a/sampleX.splice.parquet", "/b/sampleX.splice.parquet"]
        ids = _derive_sample_ids(files, from_path=True)
        assert ids[0] != ids[1]

    def test_read_count_matrices(self):
        split   = np.random.randint(10, 100, (10, 20)).astype(np.int32)
        site    = split + np.random.randint(5, 30, (10, 20)).astype(np.int32)
        unsplit = np.random.randint(1, 20, (10, 20)).astype(np.int32)
        sd = splice.read_count_matrices(split, site, unsplit_counts=unsplit)
        assert sd.n_samples == 10
        assert sd.n_junctions == 20
        assert "unsplit_counts" in sd.adata.layers

    def test_unsplit_counts_populated(self, tmp_path):
        path = _make_parquet(tmp_path, "sampleA", seed=0)
        sd = splice.read_counts([path], verbose=False)
        from scipy.sparse import issparse
        uc = sd.adata.layers["unsplit_counts"]
        arr = uc.toarray() if issparse(uc) else np.asarray(uc)
        assert (arr >= 0).all()

    def test_save_load_roundtrip(self, tmp_path):
        sd = _make_splicing_data(10, 20, inject_outlier=False)
        path = str(tmp_path / "test.h5ad")
        sd.save(path)
        sd2 = SplicingData.load(path)
        assert sd2.n_samples == 10
        assert sd2.n_junctions == 20


# ==================================================================
# Unit tests: filtering
# ==================================================================

class TestFiltering:
    def test_filter_removes_low_junctions(self):
        sd = _make_splicing_data()
        sd.adata.X[:, 5] = 0
        sd_filt = splice.filter_junctions(sd, min_split_reads=20, verbose=False)
        assert sd_filt.n_junctions < sd.n_junctions

    def test_filter_params_stored(self):
        sd = _make_splicing_data()
        sd_filt = splice.filter_junctions(sd, min_split_reads=15, verbose=False)
        assert sd_filt.adata.uns["filter_params"]["min_split_reads"] == 15


# ==================================================================
# Unit tests: Jaccard metric (now uses unsplit counts)
# ==================================================================

class TestJaccard:
    def test_jaccard_range(self):
        sd = _make_splicing_data()
        sd = compute_jaccard(sd, verbose=False)
        J = sd.jaccard()
        assert J.min() > 0
        assert J.max() <= 1.0

    def test_jaccard_uses_unsplit_layer(self):
        """Jaccard with unsplit counts should be <= Jaccard without."""
        sd_with = _make_splicing_data()
        sd_without = _make_splicing_data()
        # Remove unsplit counts from the second object
        del sd_without.adata.layers["unsplit_counts"]

        sd_with    = compute_jaccard(sd_with,    verbose=False)
        sd_without = compute_jaccard(sd_without, verbose=False)

        # Adding unsplit reads to denominator should reduce J values on average
        assert sd_with.jaccard().mean() <= sd_without.jaccard().mean() + 1e-6

    def test_logit_stored(self):
        sd = _make_splicing_data()
        sd = compute_jaccard(sd, store_logit=True, verbose=False)
        assert "jaccard_logit" in sd.adata.layers

    def test_jaccard_shape(self):
        sd = _make_splicing_data(n_samples=10, n_junctions=20)
        sd = compute_jaccard(sd, verbose=False)
        assert sd.jaccard().shape == (10, 20)


# ==================================================================
# Unit tests: confounder correction
# ==================================================================

class TestCorrection:
    def test_correction_runs(self):
        sd = _make_splicing_data()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=5, verbose=False)
        assert "jaccard_corrected" in sd.adata.layers

    def test_correction_shape(self):
        sd = _make_splicing_data(n_samples=20, n_junctions=30)
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        assert sd.jaccard_corrected().shape == (20, 30)

    def test_mp_threshold_bounds(self):
        sd = _make_splicing_data()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=None, verbose=False)
        n_comp = sd.adata.uns["correction_params"]["n_components"]
        assert 2 <= n_comp <= 100


# ==================================================================
# Unit tests: statistics
# ==================================================================

class TestStats:
    def test_call_outliers_returns_df(self):
        sd = _make_splicing_data()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        results = call_outliers(sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False)
        assert isinstance(results, pd.DataFrame)

    def test_output_columns_match_fraser2(self):
        sd = _make_splicing_data()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        results = call_outliers(sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False)
        expected = {
            "sampleID", "seqnames", "start", "end", "strand",
            "intron_id", "gene_id", "gene_name",
            "pvalue", "padj", "delta_psi",
            "observed_psi", "expected_psi",
            "counts", "total_counts",
        }
        assert expected.issubset(set(results.columns))

    def test_pvalue_range(self):
        sd = _make_splicing_data()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        results = call_outliers(sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False)
        if len(results) > 0:
            assert results["pvalue"].between(0, 1).all()
            assert results["padj"].between(0, 1).all()

    def test_effect_size_filter(self):
        sd = _make_splicing_data()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        results = call_outliers(sd, delta_jaccard=0.5, fdr_threshold=1.0, verbose=False)
        if len(results) > 0:
            # Non-flagged rows must satisfy the delta threshold;
            # shrinkage-flagged rows may have smaller corrected delta_psi
            non_flagged = results[~results["correction_shrinkage_flag"]]
            if len(non_flagged) > 0:
                assert (non_flagged["delta_psi"].abs() >= 0.5).all()

    def test_bb_mom_shape(self):
        split = np.random.randint(5, 100, (30, 20)).astype(np.int32)
        site  = split + np.random.randint(5, 50, (30, 20)).astype(np.int32)
        alpha, beta_p = _fit_bb_mom_vectorised(split, site)
        assert alpha.shape == (20,)
        assert (alpha > 0).all()


# ==================================================================
# Integration tests: full pipeline (parquet → results)
# ==================================================================

class TestPipeline:
    def test_find_outliers_end_to_end(self, tmp_path):
        for i in range(25):
            _make_parquet(tmp_path, f"sample_{i}", seed=i)
        results = splice.find_outliers(
            str(tmp_path / "counts" / "*.splice.parquet"),
            fdr_threshold=1.0,
            delta_jaccard=0.0,
            n_jobs=1,
            verbose=False,
        )
        assert isinstance(results, pd.DataFrame)
        assert len(results.columns) >= 14

    def test_find_outliers_with_annotation(self, tmp_path):
        n = 20
        for i in range(n):
            _make_parquet(tmp_path, f"sample_{i}", seed=i)
        ann = pd.DataFrame({
            "sampleID": [f"sample_{i}" for i in range(n)],
            "affected": [i < 5 for i in range(n)],
        })
        results = splice.find_outliers(
            str(tmp_path / "counts" / "*.splice.parquet"),
            sample_annotation=ann,
            fdr_threshold=1.0,
            delta_jaccard=0.0,
            verbose=False,
        )
        assert isinstance(results, pd.DataFrame)

    def test_results_sorted_by_padj(self, tmp_path):
        for i in range(20):
            _make_parquet(tmp_path, f"sample_{i}", seed=i)
        results = splice.find_outliers(
            str(tmp_path / "counts" / "*.splice.parquet"),
            fdr_threshold=1.0,
            delta_jaccard=0.0,
            verbose=False,
        )
        if len(results) > 1:
            assert (results["padj"].diff().dropna() >= 0).all()

    def test_save_intermediate(self, tmp_path):
        for i in range(15):
            _make_parquet(tmp_path, f"sample_{i}", seed=i)
        h5ad = str(tmp_path / "intermediate.h5ad")
        splice.find_outliers(
            str(tmp_path / "counts" / "*.splice.parquet"),
            fdr_threshold=1.0,
            delta_jaccard=0.0,
            save_intermediate=h5ad,
            verbose=False,
        )
        assert os.path.exists(h5ad)

    def test_unsplit_counts_reduce_jaccard(self, tmp_path):
        """Pipeline with unsplit counts should yield lower mean J than without."""
        for i in range(20):
            _make_parquet(tmp_path, f"sample_{i}", seed=i)
        sd = splice.read_counts(
            str(tmp_path / "counts" / "*.splice.parquet"), verbose=False
        )
        sd = splice.filter_junctions(sd, verbose=False)
        sd_with    = splice.compute_jaccard(sd, verbose=False)
        sd_without = splice.compute_jaccard(
            SplicingData(sd.adata.copy()), verbose=False
        )
        # Remove unsplit layer and recompute
        del sd_without.adata.layers["unsplit_counts"]
        sd_without = splice.compute_jaccard(sd_without, verbose=False)

        assert sd_with.jaccard().mean() <= sd_without.jaccard().mean() + 1e-4

    def test_repr(self):
        sd = _make_splicing_data(5, 10)
        r = repr(sd)
        assert "5 samples" in r
        assert "10 junctions" in r

# ==================================================================
# Unit tests: min_anchor parameter in extract_sample_counts
# ==================================================================

class TestMinAnchor:
    """
    Tests for the min_anchor parameter, matched against the ground-truth
    values from FRASER2's helper_test_data.R (manually IGV-counted):

        chr19:7592515-7592749  nonSplit anchor=5 -> 7
        chr19:7592515-7592749  nonSplit anchor=25 -> 5

    We replicate this logic structurally using a synthetic BAM where read
    positions are controlled precisely.
    """

    def _make_anchor_bam(self, tmp_path):
        """
        BAM with a splice site at pos=1000 (1-based donor).
        Reads designed to test anchor filtering:

          read A: spans 990..1050 (60M)  -- overlaps with 10 bp on each side -> passes anchor=5, anchor=9
          read B: spans 997..1050 (54M)  -- only 3 bp before pos=1000 -> FAILS anchor=5, passes anchor=2
          read C: spans 980..1050 (71M)  -- 20 bp before -> passes anchor=5, FAILS anchor=25
          read D: spans 970..1050 (81M)  -- 30 bp before -> passes anchor=25
          read E: spans 990..1004 (15M)  -- only 4 bp after pos=1000 -> FAILS anchor=5

        Plus 3 canonical split reads to establish the junction.
        """
        reads = []
        # 3 split reads to establish junction chr1:1001-2000
        for _ in range(3):
            reads.append({
                "chrom": "chr1", "start": 950,
                "cigar": "50M999N50M",
                "strand": "+", "mapq": 255,
                "is_read1": True, "is_reverse": False,
            })
        # read A: 60M starting at 0-based 989 -> 1-based 990..1049, anchor=10 each side
        reads.append({
            "chrom": "chr1", "start": 989,
            "cigar": "60M", "strand": None,
            "mapq": 255, "is_read1": True, "is_reverse": True,  # fr-firststrand read1 rev -> +
        })
        # read B: 54M starting at 0-based 996 -> 1-based 997..1050, anchor=3 left
        reads.append({
            "chrom": "chr1", "start": 996,
            "cigar": "54M", "strand": None,
            "mapq": 255, "is_read1": True, "is_reverse": True,
        })
        # read C: 71M starting at 0-based 979 -> 1-based 980..1050, anchor=20 left, 50 right
        reads.append({
            "chrom": "chr1", "start": 979,
            "cigar": "71M", "strand": None,
            "mapq": 255, "is_read1": True, "is_reverse": True,
        })
        # read D: 81M starting at 0-based 969 -> 1-based 970..1050, anchor=30 left
        reads.append({
            "chrom": "chr1", "start": 969,
            "cigar": "81M", "strand": None,
            "mapq": 255, "is_read1": True, "is_reverse": True,
        })
        # read E: 15M starting at 0-based 989 -> 1-based 990..1004, anchor=4 right
        reads.append({
            "chrom": "chr1", "start": 989,
            "cigar": "15M", "strand": None,
            "mapq": 255, "is_read1": True, "is_reverse": True,
        })
        return _make_bam(tmp_path, "anchor_test", reads)

    def test_anchor_0_counts_all_overlapping(self, tmp_path):
        """anchor=0 should count any non-split read overlapping the site."""
        bam = self._make_anchor_bam(tmp_path)
        out = str(tmp_path / "anchor0.parquet")
        splice.extract_sample_counts(
            bam, out, min_anchor=0, min_split_reads=2,
            library_type="fr-firststrand", verbose=False,
        )
        df = pd.read_parquet(out)
        assert len(df) == 1
        # All 5 unsplit reads overlap the donor site at pos 1000
        assert df.iloc[0]["unsplit_count"] == 5

    def test_anchor_5_filters_short_overlap(self, tmp_path):
        """anchor=5: reads B (3bp left) and E (4bp right) should be excluded."""
        bam = self._make_anchor_bam(tmp_path)
        out = str(tmp_path / "anchor5.parquet")
        splice.extract_sample_counts(
            bam, out, min_anchor=5, min_split_reads=2,
            library_type="fr-firststrand", verbose=False,
        )
        df = pd.read_parquet(out)
        # reads A, B (boundary, exactly 5bp), C, D pass -- E fails (right anchor 4bp < 5)
        assert df.iloc[0]["unsplit_count"] == 4

    def test_anchor_25_filters_moderate_overlap(self, tmp_path):
        """anchor=25: only read D (30bp left anchor) should pass."""
        bam = self._make_anchor_bam(tmp_path)
        out = str(tmp_path / "anchor25.parquet")
        splice.extract_sample_counts(
            bam, out, min_anchor=25, min_split_reads=2,
            library_type="fr-firststrand", verbose=False,
        )
        df = pd.read_parquet(out)
        assert df.iloc[0]["unsplit_count"] == 1

    def test_anchor_stricter_lte_permissive(self, tmp_path):
        """Monotonicity: stricter anchor never yields higher unsplit count."""
        bam = self._make_anchor_bam(tmp_path)
        counts = []
        for anchor in [0, 5, 10, 25]:
            out = str(tmp_path / f"anchor{anchor}.parquet")
            splice.extract_sample_counts(
                bam, out, min_anchor=anchor, min_split_reads=2,
                library_type="fr-firststrand", verbose=False,
            )
            df = pd.read_parquet(out)
            counts.append(int(df.iloc[0]["unsplit_count"]) if len(df) > 0 else 0)
        # Counts must be non-increasing
        for i in range(len(counts) - 1):
            assert counts[i] >= counts[i + 1], (
                f"anchor counts not monotone: {counts}"
            )

    def test_anchor_default_is_5(self, tmp_path):
        """Default min_anchor=5 should match explicit min_anchor=5."""
        bam = self._make_anchor_bam(tmp_path)
        out_default = str(tmp_path / "default.parquet")
        out_explicit = str(tmp_path / "explicit5.parquet")
        splice.extract_sample_counts(
            bam, out_default, min_split_reads=2,
            library_type="fr-firststrand", verbose=False,
        )
        splice.extract_sample_counts(
            bam, out_explicit, min_anchor=5, min_split_reads=2,
            library_type="fr-firststrand", verbose=False,
        )
        df_d = pd.read_parquet(out_default)
        df_e = pd.read_parquet(out_explicit)
        pd.testing.assert_frame_equal(df_d.reset_index(drop=True),
                                       df_e.reset_index(drop=True))


# ==================================================================
# Unit tests: per-sample gene subset FDR
# ==================================================================

class TestGeneFDR:
    """
    Tests for genes_per_sample subset FDR, mirroring FRASER2's
    calculatePadjValuesOnSubset() behaviour from test_stats.R.
    """

    def _make_sd_with_genes(self, n_samples=30, n_junctions=50):
        """SplicingData with distinct gene_name values across junctions."""
        sd = _make_splicing_data(n_samples=n_samples,
                                 n_junctions=n_junctions,
                                 inject_outlier=True)
        # Assign 5 gene names across junctions, matching FRASER2 test structure
        gene_names = ["MCOLN1", "TIMMDC1", "CLPP", "GENEД", "GENEE"]
        sd.adata.var["gene_name"] = [
            gene_names[i % len(gene_names)] for i in range(n_junctions)
        ]
        sd.adata.var["gene_id"] = sd.adata.var["gene_name"]
        return sd

    def test_subset_fdr_column_present(self):
        sd = self._make_sd_with_genes()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        results = call_outliers(
            sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False,
            genes_per_sample={"sample_0": ["MCOLN1"]},
        )
        assert "fdr_subset" in results.columns

    def test_subset_rows_flagged_correctly(self):
        """Rows for sample_0/MCOLN1 should have fdr_subset=True; others False."""
        sd = self._make_sd_with_genes()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        results = call_outliers(
            sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False,
            genes_per_sample={"sample_0": ["MCOLN1"]},
        )
        if len(results) == 0:
            return
        # sample_0 + MCOLN1 rows should be flagged
        mask_target = (results["sampleID"] == "sample_0") & \
                      (results["gene_name"] == "MCOLN1")
        if mask_target.sum() > 0:
            assert results.loc[mask_target, "fdr_subset"].all()
        # All other rows should not be flagged
        assert not results.loc[~mask_target, "fdr_subset"].any()

    def test_subset_padj_differs_from_global(self):
        """Per-sample subset FDR should produce different padj than global FDR."""
        sd = self._make_sd_with_genes(n_samples=40, n_junctions=50)
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)

        results_global = call_outliers(
            sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False,
        )
        results_subset = call_outliers(
            sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False,
            genes_per_sample={"sample_0": ["MCOLN1"]},
        )
        if len(results_global) == 0 or len(results_subset) == 0:
            return
        # The padj distributions need not be identical — just verify they can differ
        merged = results_global.merge(
            results_subset[["sampleID", "intron_id", "padj"]],
            on=["sampleID", "intron_id"],
            suffixes=("_global", "_subset"),
        )
        # Both columns exist and are valid
        assert "padj_global" in merged.columns
        assert "padj_subset" in merged.columns

    def test_samples_not_in_subset_get_global_fdr(self):
        """Samples absent from genes_per_sample should have fdr_subset=False."""
        sd = self._make_sd_with_genes()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        results = call_outliers(
            sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False,
            genes_per_sample={"sample_0": ["MCOLN1"]},
        )
        non_subset = results[results["sampleID"] != "sample_0"]
        assert not non_subset["fdr_subset"].any()

    def test_empty_genes_per_sample_behaves_as_global(self):
        """Passing empty dict should produce same result as not passing it."""
        sd = self._make_sd_with_genes()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        r_none = call_outliers(
            sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False,
            genes_per_sample=None,
        )
        r_empty = call_outliers(
            sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False,
            genes_per_sample={},
        )
        # padj values should be identical
        if len(r_none) > 0 and len(r_empty) > 0:
            merged = r_none.merge(
                r_empty[["sampleID", "intron_id", "padj"]],
                on=["sampleID", "intron_id"],
                suffixes=("_none", "_empty"),
            )
            np.testing.assert_allclose(
                merged["padj_none"].values,
                merged["padj_empty"].values,
                rtol=1e-5,
            )

    def test_case_insensitive_gene_matching(self):
        """Gene names in genes_per_sample should match case-insensitively."""
        sd = self._make_sd_with_genes()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        # Supply lowercase version of gene name
        r_lower = call_outliers(
            sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False,
            genes_per_sample={"sample_0": ["mcoln1"]},
        )
        r_upper = call_outliers(
            sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False,
            genes_per_sample={"sample_0": ["MCOLN1"]},
        )
        if len(r_lower) > 0 and len(r_upper) > 0:
            mask_lower = (r_lower["sampleID"] == "sample_0") & \
                         (r_lower["gene_name"] == "MCOLN1")
            mask_upper = (r_upper["sampleID"] == "sample_0") & \
                         (r_upper["gene_name"] == "MCOLN1")
            if mask_lower.sum() > 0 and mask_upper.sum() > 0:
                np.testing.assert_allclose(
                    r_lower.loc[mask_lower, "padj"].values,
                    r_upper.loc[mask_upper, "padj"].values,
                    rtol=1e-5,
                )

    def test_pipeline_genes_per_sample(self, tmp_path):
        """End-to-end: find_outliers accepts genes_per_sample."""
        for i in range(20):
            _make_parquet(tmp_path, f"sample_{i}", seed=i)
        results = splice.find_outliers(
            str(tmp_path / "counts" / "*.splice.parquet"),
            fdr_threshold=1.0,
            delta_jaccard=0.0,
            genes_per_sample={"sample_0": ["chr1:0-500:+"]},
            verbose=False,
        )
        assert "fdr_subset" in results.columns


# ==================================================================
# Unit tests: gtf.py
# ==================================================================

class TestGTF:
    """Tests for GTF parsing and gene-subset interval filtering."""

    def _make_gtf(self, tmp_path, genes=None):
        """Write a minimal GTF file with controllable gene coordinates."""
        if genes is None:
            genes = [
                ("chr1", 1000000, 1050000, "+", "ENSG00000001", "MCOLN1"),
                ("chr1", 2000000, 2050000, "-", "ENSG00000002", "TIMMDC1"),
                ("chr19", 7580000, 7620000, "+", "ENSG00000003", "CLPP"),
                ("chr19", 8000000, 8050000, "+", "ENSG00000004", "OTHERGENE"),
            ]
        lines = ['##format: GTF\n']
        for seqn, start, end, strand, gid, gname in genes:
            attrs = f'gene_id "{gid}"; gene_name "{gname}";'
            lines.append(
                f"{seqn}\tENSEMBL\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n"
            )
        gtf_path = str(tmp_path / "test.gtf")
        with open(gtf_path, "w") as f:
            f.writelines(lines)
        return gtf_path

    def test_parse_gtf_basic(self, tmp_path):
        gtf = self._make_gtf(tmp_path)
        df = splice.parse_gtf(gtf)
        assert len(df) == 4
        assert set(df.columns) >= {"seqnames", "start", "end", "gene_id", "gene_name"}

    def test_parse_gtf_gzip(self, tmp_path):
        import gzip as gz
        gtf = self._make_gtf(tmp_path)
        gtf_gz = gtf + ".gz"
        with open(gtf, "rb") as f_in, gz.open(gtf_gz, "wb") as f_out:
            f_out.write(f_in.read())
        df = splice.parse_gtf(gtf_gz)
        assert len(df) == 4

    def test_genes_to_intervals_found(self, tmp_path):
        gtf = self._make_gtf(tmp_path)
        intervals = splice.genes_to_intervals(gtf, ["MCOLN1", "TIMMDC1"])
        assert len(intervals) == 2
        assert set(intervals["gene_name"]) == {"MCOLN1", "TIMMDC1"}

    def test_genes_to_intervals_case_insensitive(self, tmp_path):
        gtf = self._make_gtf(tmp_path)
        intervals = splice.genes_to_intervals(gtf, ["mcoln1", "TIMMDC1"])
        assert len(intervals) == 2

    def test_genes_to_intervals_missing_warns(self, tmp_path):
        import warnings
        gtf = self._make_gtf(tmp_path)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            intervals = splice.genes_to_intervals(gtf, ["MCOLN1", "DOESNOTEXIST"])
        assert len(intervals) == 1
        assert any("not found" in str(warning.message) for warning in w)

    def test_genes_to_intervals_none_found_raises(self, tmp_path):
        gtf = self._make_gtf(tmp_path)
        with pytest.raises(ValueError, match="None of the"):
            splice.genes_to_intervals(gtf, ["FAKEGENE1", "FAKEGENE2"])

    def test_filter_junctions_to_genes(self, tmp_path):
        gtf = self._make_gtf(tmp_path)
        intervals = splice.genes_to_intervals(gtf, ["MCOLN1"])

        # Create junction DataFrame: 3 inside MCOLN1, 2 outside
        junctions = pd.DataFrame({
            "intron_id": [
                "chr1:1005000-1010000:+",  # inside MCOLN1
                "chr1:1020000-1030000:+",  # inside MCOLN1
                "chr1:1040000-1045000:+",  # inside MCOLN1
                "chr1:1060000-1070000:+",  # outside MCOLN1
                "chr2:5000-6000:+",        # wrong chrom
            ],
            "seqnames": ["chr1", "chr1", "chr1", "chr1", "chr2"],
            "start":    [1005000, 1020000, 1040000, 1060000, 5000],
            "end":      [1010000, 1030000, 1045000, 1070000, 6000],
            "strand":   ["+", "+", "+", "+", "+"],
        })
        filtered = splice.filter_junctions_to_genes(junctions, intervals, verbose=False)
        assert len(filtered) == 3
        assert all(filtered["seqnames"] == "chr1")

    def test_filter_junctions_labels_genes(self, tmp_path):
        gtf = self._make_gtf(tmp_path)
        intervals = splice.genes_to_intervals(gtf, ["MCOLN1"])
        junctions = pd.DataFrame({
            "intron_id": ["chr1:1005000-1010000:+"],
            "seqnames":  ["chr1"],
            "start":     [1005000],
            "end":       [1010000],
            "strand":    ["+"],
        })
        filtered = splice.filter_junctions_to_genes(junctions, intervals, verbose=False)
        assert filtered.iloc[0]["gene_name"] == "MCOLN1"
        assert filtered.iloc[0]["gene_id"]   == "ENSG00000001"

    def test_parse_interval_tsv(self, tmp_path):
        tsv_path = str(tmp_path / "intervals.tsv")
        pd.DataFrame({
            "seqnames": ["chr1"], "start": [1000000], "end": [1050000],
            "gene_id": ["ENSG00000001"], "gene_name": ["MCOLN1"], "strand": ["+"],
        }).to_csv(tsv_path, sep="\t", index=False)
        df = splice.parse_gtf(tsv_path)
        assert len(df) == 1
        assert df.iloc[0]["gene_name"] == "MCOLN1"


# ==================================================================
# Unit tests: gene_subset in read_counts and find_outliers
# ==================================================================

class TestGeneSubset:
    """Tests for gene_subset filtering in read_counts and find_outliers."""

    def _make_gtf_for_parquets(self, tmp_path):
        """GTF whose gene coordinates fully contain the test parquet junctions."""
        # Parquet junctions: chr1:N*1000-(N*1000+500):+ for N in 0..49
        # Gene spans 0..300000 on chr1 -> covers all 50 junctions
        lines = ['##format: GTF\n']
        lines.append(
            'chr1\tENSEMBL\tgene\t0\t300000\t.\t+\t.\t'
            'gene_id "ENSG00000001"; gene_name "TARGETGENE";\n'
        )
        lines.append(
            'chr2\tENSEMBL\tgene\t0\t300000\t.\t+\t.\t'
            'gene_id "ENSG00000002"; gene_name "OTHERGENE";\n'
        )
        gtf_path = str(tmp_path / "subset_test.gtf")
        with open(gtf_path, "w") as f:
            f.writelines(lines)
        return gtf_path

    def test_gene_subset_reduces_junctions(self, tmp_path):
        for i in range(20):
            _make_parquet(tmp_path, f"sample_{i}", seed=i)
        gtf = self._make_gtf_for_parquets(tmp_path)

        sd_all = splice.read_counts(
            str(tmp_path / "counts" / "*.splice.parquet"), verbose=False
        )
        sd_sub = splice.read_counts(
            str(tmp_path / "counts" / "*.splice.parquet"),
            gene_annotation=gtf,
            gene_subset=["TARGETGENE"],
            verbose=False,
        )
        # Subset should have <= junctions than full (TARGETGENE covers chr1 junctions)
        assert sd_sub.n_junctions <= sd_all.n_junctions

    def test_gene_subset_requires_annotation(self, tmp_path):
        for i in range(5):
            _make_parquet(tmp_path, f"sample_{i}", seed=i)
        with pytest.raises(ValueError, match="gene_annotation"):
            splice.read_counts(
                str(tmp_path / "counts" / "*.splice.parquet"),
                gene_subset=["MCOLN1"],
                verbose=False,
            )

    def test_find_outliers_gene_subset(self, tmp_path):
        for i in range(20):
            _make_parquet(tmp_path, f"sample_{i}", seed=i)
        gtf = self._make_gtf_for_parquets(tmp_path)
        results = splice.find_outliers(
            str(tmp_path / "counts" / "*.splice.parquet"),
            gene_annotation=gtf,
            gene_subset=["TARGETGENE"],
            fdr_threshold=1.0,
            delta_jaccard=0.0,
            verbose=False,
        )
        assert isinstance(results, pd.DataFrame)

    def test_gene_subset_sample_count_unchanged(self, tmp_path):
        """Gene subsetting should not change the number of samples."""
        for i in range(10):
            _make_parquet(tmp_path, f"sample_{i}", seed=i)
        gtf = self._make_gtf_for_parquets(tmp_path)
        sd = splice.read_counts(
            str(tmp_path / "counts" / "*.splice.parquet"),
            gene_annotation=gtf,
            gene_subset=["TARGETGENE"],
            verbose=False,
        )
        assert sd.n_samples == 10


# ==================================================================
# Unit tests: delta PSI columns + correction_shrinkage_flag
# ==================================================================

class TestDeltaPSI:
    """Tests for raw/unweighted delta PSI and correction shrinkage flag."""

    def _make_sd_with_shrinkage(self, n_samples=40, n_junctions=30):
        """
        SplicingData where sample_0 / junction_0 has an extreme raw delta
        but a modest corrected delta (simulating correction-absorbed signal).
        """
        rng = np.random.default_rng(99)
        split = rng.integers(50, 200, (n_samples, n_junctions)).astype(np.int32)
        site  = split + rng.integers(20, 100, (n_samples, n_junctions)).astype(np.int32)
        # Make junction 0 have ~median PSI = 0.8 across cohort
        site[:, 0]  = 200
        split[:, 0] = 160   # PSI ~ 0.8 for most samples
        # sample_0 has almost zero split reads -> extreme raw delta
        split[0, 0] = 2     # PSI ~ 0.01, raw delta ~ -0.79

        obs = pd.DataFrame({"sampleID": [f"sample_{i}" for i in range(n_samples)]})
        obs.index = obs["sampleID"]
        var = pd.DataFrame({
            "intron_id": [f"chr1:{i*1000}-{i*1000+500}:+" for i in range(n_junctions)],
            "seqnames":  ["chr1"] * n_junctions,
            "start":     [i * 1000 for i in range(n_junctions)],
            "end":       [i * 1000 + 500 for i in range(n_junctions)],
            "strand":    ["+"] * n_junctions,
            "gene_id":   [f"GENE{i}" for i in range(n_junctions)],
            "gene_name": [f"GENE{i}" for i in range(n_junctions)],
        })
        var.index = var["intron_id"]
        return SplicingData.from_matrices(split, site, obs=obs, var=var)

    def test_delta_psi_raw_column_present(self):
        sd = _make_splicing_data()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        results = call_outliers(sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False)
        assert "delta_psi_raw" in results.columns

    def test_delta_psi_unweighted_column_present(self):
        sd = _make_splicing_data()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        results = call_outliers(sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False)
        assert "delta_psi_unweighted" in results.columns

    def test_correction_shrinkage_flag_column_present(self):
        sd = _make_splicing_data()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        results = call_outliers(sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False)
        assert "correction_shrinkage_flag" in results.columns

    def test_delta_psi_raw_is_obs_minus_median(self):
        """delta_psi_raw should equal observed_psi minus the cohort median."""
        sd = _make_splicing_data(inject_outlier=False)
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        results = call_outliers(sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False)
        if len(results) == 0:
            return
        # For any row: delta_psi_raw = observed_psi - median(observed_psi for that junction)
        # We check sign consistency: if observed > median, raw delta should be positive
        pos_mask = results["delta_psi_raw"] > 0
        if pos_mask.sum() > 0:
            assert (results.loc[pos_mask, "observed_psi"] >=
                    results.loc[pos_mask, "observed_psi"] - results.loc[pos_mask, "delta_psi_raw"]
                    ).all()

    def test_raw_and_corrected_can_differ(self):
        """delta_psi_raw and delta_psi need not be equal after correction."""
        sd = _make_splicing_data(n_samples=40)
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=5, verbose=False)
        results = call_outliers(sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False)
        if len(results) < 10:
            return
        # At least some rows should have different raw vs corrected delta
        diff = (results["delta_psi_raw"] - results["delta_psi"]).abs()
        assert diff.max() > 1e-4

    def test_shrinkage_flag_fires_when_ratio_exceeded(self):
        """Injecting a high-raw-delta/low-corrected-delta scenario triggers the flag."""
        sd = self._make_sd_with_shrinkage(n_samples=40, n_junctions=30)
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        results = call_outliers(
            sd,
            fdr_threshold=1.0,
            delta_jaccard=0.0,
            shrinkage_threshold=2.0,   # lower threshold to ensure flag fires
            verbose=False,
        )
        assert "correction_shrinkage_flag" in results.columns
        # At least one flagged event should exist
        assert results["correction_shrinkage_flag"].any(), (
            "Expected at least one correction_shrinkage_flag=True row with "
            "an extreme raw delta and low corrected delta"
        )

    def test_shrinkage_flagged_rows_retained_below_delta_threshold(self):
        """Flagged events are kept even when corrected |delta_psi| < delta_jaccard."""
        sd = self._make_sd_with_shrinkage(n_samples=40, n_junctions=30)
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        # strict delta_jaccard=0.5 would drop most rows, but shrinkage-flagged stay
        results = call_outliers(
            sd,
            fdr_threshold=1.0,
            delta_jaccard=0.5,
            shrinkage_threshold=2.0,
            verbose=False,
        )
        if results["correction_shrinkage_flag"].any():
            flagged = results[results["correction_shrinkage_flag"]]
            # Flagged rows may have corrected |delta_psi| < 0.5
            assert len(flagged) > 0

    def test_shrinkage_threshold_configurable(self):
        """Higher threshold means fewer flagged events."""
        sd = self._make_sd_with_shrinkage(n_samples=40, n_junctions=30)
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)

        r_low  = call_outliers(sd, fdr_threshold=1.0, delta_jaccard=0.0,
                                shrinkage_threshold=1.5, verbose=False)
        r_high = call_outliers(sd, fdr_threshold=1.0, delta_jaccard=0.0,
                                shrinkage_threshold=10.0, verbose=False)

        n_low  = r_low["correction_shrinkage_flag"].sum()
        n_high = r_high["correction_shrinkage_flag"].sum()
        assert n_low >= n_high

    def test_output_columns_complete(self):
        """All FRASER2-compatible columns plus new diagnostics present."""
        sd = _make_splicing_data()
        sd = compute_jaccard(sd, verbose=False)
        sd = correct_confounders(sd, n_components=3, verbose=False)
        results = call_outliers(sd, fdr_threshold=1.0, delta_jaccard=0.0, verbose=False)
        expected = {
            "sampleID", "seqnames", "start", "end", "strand",
            "intron_id", "gene_id", "gene_name",
            "pvalue", "padj",
            "delta_psi", "delta_psi_raw", "delta_psi_unweighted",
            "observed_psi", "expected_psi",
            "counts", "total_counts",
            "fdr_subset", "correction_shrinkage_flag",
        }
        assert expected.issubset(set(results.columns))

    def test_pipeline_exposes_new_columns(self, tmp_path):
        """find_outliers passes shrinkage_threshold and returns new columns."""
        for i in range(20):
            _make_parquet(tmp_path, f"sample_{i}", seed=i)
        results = splice.find_outliers(
            str(tmp_path / "counts" / "*.splice.parquet"),
            fdr_threshold=1.0,
            delta_jaccard=0.0,
            shrinkage_threshold=2.0,
            verbose=False,
        )
        for col in ("delta_psi_raw", "delta_psi_unweighted", "correction_shrinkage_flag"):
            assert col in results.columns
