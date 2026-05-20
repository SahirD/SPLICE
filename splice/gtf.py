"""
gtf.py — GTF parsing and gene-subset interval filtering
=========================================================
Provides utilities to:
  1. Parse a standard Ensembl/GENCODE GTF file into a gene interval DataFrame.
  2. Accept a pre-built TSV (seqnames, start, end, gene_id, gene_name) as a
     faster alternative to re-parsing a full GTF on every run.
  3. Given a list of gene names/IDs, return the corresponding genomic intervals.
  4. Filter a junction DataFrame (or SplicingData) to only junctions whose
     donor AND acceptor fall within at least one target gene interval.

GTF format notes
----------------
Standard GTF (GFF2) column layout (0-indexed, tab-separated):
  0  seqnames   (chromosome)
  1  source
  2  feature    (gene, transcript, exon, CDS, ...)
  3  start      (1-based, inclusive)
  4  end        (1-based, inclusive)
  5  score
  6  strand
  7  frame
  8  attributes (key "value"; key "value"; ...)

We extract 'gene' feature rows only and parse gene_id and gene_name from
column 8. This is fast even for large GTFs because we skip all non-gene rows.

Performance
-----------
For large GTFs (GENCODE hg38 ~ 60 MB uncompressed) parsing takes ~2-3 s.
The result can be saved as a TSV with to_csv() and reloaded via the TSV
path on subsequent runs — much faster than re-parsing the GTF each time.
Gzipped GTFs (.gtf.gz) are handled transparently via gzip.open().
"""

from __future__ import annotations

import gzip
import os
import re
from typing import List, Optional, Union

import numpy as np
import pandas as pd


# ==================================================================
# Public API
# ==================================================================

def parse_gtf(
    gtf_path: str,
    feature: str = "gene",
) -> pd.DataFrame:
    """
    Parse a GTF file and return a gene interval DataFrame.

    Accepts both plain-text (.gtf) and gzip-compressed (.gtf.gz) files.
    Also accepts a pre-built TSV/CSV if that is faster for repeated runs.

    Parameters
    ----------
    gtf_path : path to a GTF (.gtf / .gtf.gz) or pre-built interval TSV/CSV.
        TSV/CSV must contain columns: seqnames, start, end, gene_id, gene_name.
    feature : GTF feature type to extract (default 'gene').
        Use 'transcript' to get per-transcript intervals instead.

    Returns
    -------
    DataFrame with columns:
        seqnames, start (int32), end (int32), strand, gene_id, gene_name
    """
    path = str(gtf_path)

    # Pre-built TSV/CSV — skip GTF parsing entirely
    if path.endswith((".tsv", ".csv", ".tsv.gz", ".csv.gz")):
        return _load_interval_tsv(path)

    # GTF / GTF.GZ
    return _parse_gtf_file(path, feature)


def genes_to_intervals(
    gtf_or_df: Union[str, pd.DataFrame],
    gene_subset: List[str],
) -> pd.DataFrame:
    """
    Return genomic intervals for a list of gene names or IDs.

    Parameters
    ----------
    gtf_or_df : GTF file path (str) or pre-parsed interval DataFrame.
    gene_subset : list of gene_name or gene_id values (case-insensitive).

    Returns
    -------
    DataFrame with columns: seqnames, start, end, strand, gene_id, gene_name
    Rows are deduplicated by gene_id.

    Raises
    ------
    ValueError if none of the requested genes are found in the GTF.
    """
    if isinstance(gtf_or_df, pd.DataFrame):
        intervals = gtf_or_df.copy()
    else:
        intervals = parse_gtf(gtf_or_df)

    # Normalise to lower-case for matching
    query = frozenset(g.lower() for g in gene_subset)
    mask = (
        intervals["gene_name"].str.lower().isin(query) |
        intervals["gene_id"].str.lower().isin(query)
    )
    result = intervals[mask].drop_duplicates("gene_id").reset_index(drop=True)

    if len(result) == 0:
        raise ValueError(
            f"None of the {len(gene_subset)} requested genes were found in the "
            f"annotation. Check gene name spelling and that the correct GTF is "
            f"being used.\n"
            f"  First 5 requested: {gene_subset[:5]}\n"
            f"  First 5 in GTF:    {sorted(intervals['gene_name'].dropna().unique())[:5]}"
        )

    n_found = len(result)
    n_missing = len(gene_subset) - n_found
    if n_missing > 0:
        found_names = frozenset(result["gene_name"].str.lower()) | \
                      frozenset(result["gene_id"].str.lower())
        missing = [g for g in gene_subset if g.lower() not in found_names]
        import warnings
        warnings.warn(
            f"{n_missing} gene(s) not found in annotation and will be skipped: "
            f"{missing[:10]}{'...' if len(missing) > 10 else ''}",
            UserWarning,
            stacklevel=3,
        )

    return result


def filter_junctions_to_genes(
    junction_df: pd.DataFrame,
    intervals: pd.DataFrame,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Keep only junctions whose donor AND acceptor both fall within
    at least one target gene interval.

    Parameters
    ----------
    junction_df : DataFrame with columns seqnames, start, end (1-based intron coords).
    intervals : gene interval DataFrame from genes_to_intervals().
    verbose : print summary.

    Returns
    -------
    Filtered junction_df with a 'gene_id' and 'gene_name' column filled in
    from the overlapping interval.
    """
    if len(intervals) == 0:
        return junction_df.iloc[0:0].copy()

    n_before = len(junction_df)

    # Build a per-chromosome interval lookup for fast overlap
    # {seqnames: [(gene_start, gene_end, gene_id, gene_name), ...]}
    chrom_intervals: dict = {}
    for _, row in intervals.iterrows():
        chrom = str(row["seqnames"])
        if chrom not in chrom_intervals:
            chrom_intervals[chrom] = []
        chrom_intervals[chrom].append((
            int(row["start"]), int(row["end"]),
            str(row.get("gene_id", "")),
            str(row.get("gene_name", "")),
        ))

    # Vectorised overlap check per chromosome
    keep_mask = np.zeros(len(junction_df), dtype=bool)
    gene_id_col   = [""] * len(junction_df)
    gene_name_col = [""] * len(junction_df)

    jdf = junction_df.reset_index(drop=True)

    for chrom, civals in chrom_intervals.items():
        chr_mask = jdf["seqnames"].astype(str) == chrom
        chr_idx  = np.where(chr_mask)[0]
        if len(chr_idx) == 0:
            continue

        j_starts = jdf.loc[chr_idx, "start"].values.astype(np.int64)
        j_ends   = jdf.loc[chr_idx, "end"].values.astype(np.int64)

        for g_start, g_end, g_id, g_name in civals:
            # Junction overlaps gene if donor (j_start) and acceptor (j_end)
            # both fall within [g_start, g_end]
            donor_in    = (j_starts >= g_start) & (j_starts <= g_end)
            acceptor_in = (j_ends   >= g_start) & (j_ends   <= g_end)
            in_gene = donor_in & acceptor_in

            for local_i in np.where(in_gene)[0]:
                global_i = chr_idx[local_i]
                keep_mask[global_i]    = True
                gene_id_col[global_i]   = g_id
                gene_name_col[global_i] = g_name

    result = jdf[keep_mask].copy()
    result["gene_id"]   = [gene_id_col[i]   for i in np.where(keep_mask)[0]]
    result["gene_name"] = [gene_name_col[i] for i in np.where(keep_mask)[0]]

    if verbose:
        print(
            f"[splice] gene subset filter: {n_before:,} → {len(result):,} junctions "
            f"({n_before - len(result):,} removed, "
            f"{len(intervals):,} target genes)"
        )

    return result.reset_index(drop=True)


# ==================================================================
# Internal: GTF parsing
# ==================================================================

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]+)"')


def _parse_gtf_file(path: str, feature: str) -> pd.DataFrame:
    """Parse a GTF/GTF.GZ file, returning gene-level intervals."""
    open_fn = gzip.open if path.endswith(".gz") else open
    mode    = "rt"

    records = []
    with open_fn(path, mode, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != feature:
                continue

            seqnames = parts[0]
            start    = int(parts[3])   # 1-based inclusive
            end      = int(parts[4])
            strand   = parts[6]
            attrs    = parts[8]

            attr_dict = dict(_ATTR_RE.findall(attrs))
            gene_id   = attr_dict.get("gene_id",   "")
            gene_name = attr_dict.get("gene_name", gene_id)

            records.append((seqnames, start, end, strand, gene_id, gene_name))

    if not records:
        raise ValueError(
            f"No '{feature}' features found in {path}. "
            f"Check the file is a valid GTF and the feature type is correct."
        )

    df = pd.DataFrame(
        records,
        columns=["seqnames", "start", "end", "strand", "gene_id", "gene_name"],
    )
    df["start"] = df["start"].astype(np.int32)
    df["end"]   = df["end"].astype(np.int32)
    return df


def _load_interval_tsv(path: str) -> pd.DataFrame:
    """Load a pre-built interval TSV/CSV."""
    sep = "\t" if path.replace(".gz", "").endswith(".tsv") else ","
    df  = pd.read_csv(path, sep=sep, compression="infer")

    required = {"seqnames", "start", "end", "gene_id"}
    missing  = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Interval file {path} is missing required columns: {missing}. "
            f"Required: seqnames, start, end, gene_id, gene_name."
        )
    if "gene_name" not in df.columns:
        df["gene_name"] = df["gene_id"]
    if "strand" not in df.columns:
        df["strand"] = "*"

    df["start"] = df["start"].astype(np.int32)
    df["end"]   = df["end"].astype(np.int32)
    return df
