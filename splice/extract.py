"""
extract.py -- BAM -> per-sample splice count parquet
=====================================================
Two-pass BAM scan producing a compact parquet file per sample.

Pass 1: discover all junctions from N-CIGAR reads; count split reads.
Pass 2: count unsplit (non-split) reads at each splice site with
        configurable anchor and strand filters.

When gene_subset + gene_annotation are provided, pass 1 scans the full
BAM but only retains junctions that fall within target gene coordinates.
Pass 2 then fetches reads ONLY within those gene windows, which is the
primary speed benefit for targeted panels.
"""

from __future__ import annotations

import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd

try:
    import pysam
except ImportError as e:
    raise ImportError(
        "pysam is required for BAM extraction. Install it with:\n"
        "  pip install pysam"
    ) from e

LibraryType = Literal["fr-firststrand", "fr-secondstrand", "unstranded"]

_CIGAR_N = 3


# ==================================================================
# Public API
# ==================================================================

def extract_sample_counts(
    bam_path: str,
    output_path: str,
    library_type: LibraryType = "fr-firststrand",
    min_mapping_quality: int = 255,
    min_split_reads: int = 2,
    min_anchor: int = 5,
    gene_subset: Optional[List[str]] = None,
    gene_annotation: Optional[Union[str, pd.DataFrame]] = None,
    verbose: bool = True,
) -> str:
    """
    Extract splice junction and site counts from a BAM file.

    Parameters
    ----------
    bam_path : path to coordinate-sorted, indexed BAM file (.bai required).
    output_path : path for the output parquet file.
        Convention: "<sample_id>.splice.parquet"
    library_type : strand protocol for unsplit read inference.
        'fr-firststrand' -- dUTP/reverse-stranded (Illumina TruSeq default).
        'fr-secondstrand' -- ligation/forward-stranded.
        'unstranded' -- strand not preserved.
    min_mapping_quality : minimum MAPQ (255 = STAR unique, 60 = HISAT2 unique).
    min_split_reads : minimum split reads to retain a junction (default 2).
    min_anchor : minimum bases a non-split read must extend past the splice
        site on each side (default 5, matches FRASER2 minAnchor).
    gene_subset : optional list of gene names or IDs to restrict output to.
        Junctions outside these genes are discarded. Dramatically reduces
        runtime for targeted panels when combined with gene_annotation.
        Example: ["MCOLN1", "TIMMDC1", "CLPP"]
    gene_annotation : GTF path (.gtf / .gtf.gz), pre-built interval TSV, or
        DataFrame with columns [seqnames, start, end, gene_id, gene_name].
        Required when gene_subset is provided.
    verbose : print progress messages.

    Returns
    -------
    output_path (str)

    Examples
    --------
    Full transcriptome:
    >>> splice.extract_sample_counts(
    ...     "sample.bam", "sample.splice.parquet",
    ...     library_type="fr-firststrand", min_anchor=5,
    ... )

    Targeted panel (much faster):
    >>> splice.extract_sample_counts(
    ...     "sample.bam", "sample.splice.parquet",
    ...     gene_subset=["MCOLN1", "TIMMDC1", "CLPP"],
    ...     gene_annotation="hg38.gtf.gz",
    ... )
    """
    bam_path    = str(bam_path)
    output_path = str(output_path)

    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    bai_path = bam_path + ".bai"
    bai_alt  = os.path.splitext(bam_path)[0] + ".bai"
    if not os.path.exists(bai_path) and not os.path.exists(bai_alt):
        raise FileNotFoundError(
            f"BAM index not found. Run:\n  samtools index {bam_path}"
        )

    # -- Resolve gene intervals if a subset was requested ----------
    target_intervals: Optional[pd.DataFrame] = None
    if gene_subset is not None:
        if gene_annotation is None:
            raise ValueError(
                "gene_annotation (GTF path or DataFrame) is required when "
                "gene_subset is provided."
            )
        from .gtf import parse_gtf, genes_to_intervals
        ann_df = parse_gtf(gene_annotation) if isinstance(gene_annotation, str) \
                 else gene_annotation
        target_intervals = genes_to_intervals(ann_df, gene_subset)
        if verbose:
            print(f"[splice] extract_sample_counts: targeting "
                  f"{len(target_intervals):,} gene intervals "
                  f"({len(gene_subset)} requested genes)")

    if verbose:
        sample_id = _sample_id_from_path(output_path)
        print(f"[splice] extract_sample_counts: {sample_id}  ({bam_path})")

    # -- Pass 1: discover junctions --------------------------------
    junctions, donor_totals, acceptor_totals = _pass1_split_reads(
        bam_path, min_mapping_quality
    )

    # Apply min_split_reads filter
    junctions = {k: v for k, v in junctions.items() if v >= min_split_reads}

    # Apply gene-subset filter: keep only junctions within target genes
    if target_intervals is not None and len(junctions) > 0:
        junctions = _filter_junctions_to_intervals(junctions, target_intervals)
        # Propagate gene annotation into donor/acceptor totals filtering
        donor_totals    = {k: v for k, v in donor_totals.items()
                           if any(k[0] == (seqn, pos, strand)[0]
                                  for seqn, pos, strand in [k])}
        # Simpler: rebuild donor/acceptor totals from retained junctions only
        donor_totals    = {(s, st, strand): donor_totals.get((s, st, strand), 0)
                           for s, st, en, strand in junctions}
        acceptor_totals = {(s, en, strand): acceptor_totals.get((s, en, strand), 0)
                           for s, st, en, strand in junctions}

    if not junctions:
        if verbose:
            print("[splice]   no junctions retained after filtering")
        _write_empty_parquet(output_path)
        return output_path

    if verbose:
        print(f"[splice]   {len(junctions):,} junctions retained")

    # -- Pass 2: unsplit reads (restricted to target regions) ------
    unsplit_counts = _pass2_unsplit_reads(
        bam_path, junctions, library_type, min_mapping_quality,
        min_anchor, target_intervals, verbose,
    )

    # -- Build result DataFrame ------------------------------------
    # Look up gene assignment for each junction
    gene_map: Dict[Tuple, Tuple[str, str]] = {}
    if target_intervals is not None:
        gene_map = _build_gene_map(junctions, target_intervals)

    records = []
    for (seqnames, start, end, strand), split_count in junctions.items():
        site_count = max(
            donor_totals.get(   (seqnames, start, strand), split_count),
            acceptor_totals.get((seqnames, end,   strand), split_count),
        )
        unsplit   = unsplit_counts.get((seqnames, start, end, strand), 0)
        gene_id, gene_name = gene_map.get(
            (seqnames, start, end, strand), ("", "")
        )
        records.append({
            "intron_id":     f"{seqnames}:{start}-{end}:{strand}",
            "seqnames":      seqnames,
            "start":         np.int32(start),
            "end":           np.int32(end),
            "strand":        strand,
            "split_count":   np.int32(split_count),
            "site_count":    np.int32(site_count),
            "unsplit_count": np.int32(unsplit),
            "gene_id":       gene_id,
            "gene_name":     gene_name,
        })

    df = pd.DataFrame(records)
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    df.to_parquet(output_path, index=False, compression="snappy")

    if verbose:
        size_mb = os.path.getsize(output_path) / 1e6
        print(f"[splice]   written -> {output_path}  ({size_mb:.2f} MB)")

    return output_path


# ==================================================================
# Pass 1: split reads
# ==================================================================

def _pass1_split_reads(bam_path, min_mapq):
    junctions:        Dict[Tuple, int] = defaultdict(int)
    donor_totals:     Dict[Tuple, int] = defaultdict(int)
    acceptor_totals:  Dict[Tuple, int] = defaultdict(int)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if _skip_read(read, min_mapq):
                continue
            strand = _strand_from_xs(read)
            if strand is None:
                continue
            seqnames = read.reference_name
            for intron_start, intron_end in _extract_introns_from_cigar(read):
                key = (seqnames, intron_start, intron_end, strand)
                junctions[key] += 1
                donor_totals[   (seqnames, intron_start, strand)] += 1
                acceptor_totals[(seqnames, intron_end,   strand)] += 1

    return dict(junctions), dict(donor_totals), dict(acceptor_totals)


# ==================================================================
# Pass 2: unsplit reads
# ==================================================================

def _pass2_unsplit_reads(
    bam_path, junctions, library_type, min_mapq,
    min_anchor, target_intervals, verbose,
):
    donor_sites:    Dict[Tuple, List] = defaultdict(list)
    acceptor_sites: Dict[Tuple, List] = defaultdict(list)

    for key in junctions:
        seqnames, start, end, strand = key
        donor_sites[   (seqnames, start, strand)].append(key)
        acceptor_sites[(seqnames, end,   strand)].append(key)

    site_unsplit: Dict[Tuple, int] = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        all_sites = set(donor_sites) | set(acceptor_sites)
        for seqnames, pos, strand in all_sites:
            count = _count_unsplit_at_site(
                bam, seqnames, pos, strand, library_type, min_mapq, min_anchor
            )
            site_unsplit[(seqnames, pos, strand)] = count

    unsplit_counts: Dict[Tuple, int] = {}
    for key in junctions:
        seqnames, start, end, strand = key
        unsplit_counts[key] = max(
            site_unsplit.get((seqnames, start, strand), 0),
            site_unsplit.get((seqnames, end,   strand), 0),
        )

    return unsplit_counts


def _count_unsplit_at_site(bam, seqnames, pos, strand, library_type, min_mapq, min_anchor):
    count = 0
    if min_anchor > 0:
        fetch_start = max(0, pos - min_anchor - 1)
        fetch_end   = pos + min_anchor
    else:
        fetch_start = max(0, pos - 2)
        fetch_end   = pos + 1

    try:
        for read in bam.fetch(seqnames, fetch_start, fetch_end):
            if _skip_read(read, min_mapq):
                continue
            if _read_has_intron(read):
                continue
            read_start_1based = read.reference_start + 1
            read_end_1based   = read.reference_end
            if min_anchor > 0:
                if read_start_1based > pos - min_anchor + 1:
                    continue
                if read_end_1based < pos + min_anchor:
                    continue
            else:
                if not (read_start_1based <= pos <= read_end_1based):
                    continue
            read_strand = _infer_unsplit_strand(read, library_type)
            if strand != "*" and read_strand != "*" and read_strand != strand:
                continue
            count += 1
    except (ValueError, KeyError):
        pass

    return count


# ==================================================================
# Gene subset helpers
# ==================================================================

def _filter_junctions_to_intervals(
    junctions: Dict[Tuple, int],
    intervals: pd.DataFrame,
) -> Dict[Tuple, int]:
    """Keep only junctions whose donor and acceptor both fall within a target gene."""
    # Build (seqnames -> list of (start, end)) lookup
    chrom_lookup: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for _, row in intervals.iterrows():
        chrom_lookup[str(row["seqnames"])].append(
            (int(row["start"]), int(row["end"]))
        )

    kept = {}
    for key, count in junctions.items():
        seqnames, j_start, j_end, strand = key
        gene_ranges = chrom_lookup.get(seqnames, [])
        for g_start, g_end in gene_ranges:
            if (g_start <= j_start <= g_end) and (g_start <= j_end <= g_end):
                kept[key] = count
                break
    return kept


def _build_gene_map(
    junctions: Dict[Tuple, int],
    intervals: pd.DataFrame,
) -> Dict[Tuple, Tuple[str, str]]:
    """Map each junction tuple to (gene_id, gene_name)."""
    chrom_lookup: Dict[str, List[Tuple]] = defaultdict(list)
    for _, row in intervals.iterrows():
        chrom_lookup[str(row["seqnames"])].append((
            int(row["start"]), int(row["end"]),
            str(row.get("gene_id", "")), str(row.get("gene_name", "")),
        ))

    result = {}
    for key in junctions:
        seqnames, j_start, j_end, strand = key
        for g_start, g_end, g_id, g_name in chrom_lookup.get(seqnames, []):
            if (g_start <= j_start <= g_end) and (g_start <= j_end <= g_end):
                result[key] = (g_id, g_name)
                break
        if key not in result:
            result[key] = ("", "")
    return result


# ==================================================================
# CIGAR helpers
# ==================================================================

def _extract_introns_from_cigar(read):
    introns = []
    ref_pos = read.reference_start
    for op, length in read.cigartuples:
        if op == _CIGAR_N:
            introns.append((ref_pos + 1, ref_pos + length))
            ref_pos += length
        elif op in (0, 2, 3, 7, 8):
            ref_pos += length
    return introns


def _read_has_intron(read):
    if read.cigartuples is None:
        return False
    return any(op == _CIGAR_N for op, _ in read.cigartuples)


def _read_overlaps_pos(read, pos):
    return (read.reference_start + 1) <= pos <= read.reference_end


# ==================================================================
# Strand helpers
# ==================================================================

def _strand_from_xs(read) -> Optional[str]:
    try:
        xs = read.get_tag("XS")
        if xs in ("+", "-"):
            return xs
    except KeyError:
        pass
    return None


def _infer_unsplit_strand(read, library_type: LibraryType) -> str:
    if library_type == "unstranded":
        return "*"
    is_read1   = read.is_read1
    is_reverse = read.is_reverse
    if library_type == "fr-firststrand":
        return ("-" if not is_reverse else "+") if is_read1 \
               else ("+" if not is_reverse else "-")
    if library_type == "fr-secondstrand":
        return ("+" if not is_reverse else "-") if is_read1 \
               else ("-" if not is_reverse else "+")
    return "*"


# ==================================================================
# Read filtering + utilities
# ==================================================================

def _skip_read(read, min_mapq: int) -> bool:
    return (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_duplicate
        or read.mapping_quality < min_mapq
        or read.cigartuples is None
    )


def _sample_id_from_path(path: str) -> str:
    p = Path(path)
    name = p.name
    for suffix in (".splice.parquet", ".parquet"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return p.stem


def _write_empty_parquet(output_path: str) -> None:
    pd.DataFrame(columns=[
        "intron_id", "seqnames", "start", "end", "strand",
        "split_count", "site_count", "unsplit_count", "gene_id", "gene_name",
    ]).to_parquet(output_path, index=False)
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
