"""
io.py -- read per-sample parquet files and build SplicingData
=============================================================
read_counts()         -- primary entry point; reads *.splice.parquet files
read_count_matrices() -- kept for pre-built count matrices
"""

from __future__ import annotations

import glob
import os
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from tqdm import tqdm

from .data import SplicingData


# ==================================================================
# Public API
# ==================================================================

def read_counts(
    pattern: Union[str, List[str]],
    sample_annotation: Optional[Union[str, pd.DataFrame]] = None,
    sample_id_from_path: bool = True,
    gene_annotation: Optional[Union[str, pd.DataFrame]] = None,
    gene_subset: Optional[List[str]] = None,
    n_jobs: int = 1,
    verbose: bool = True,
) -> SplicingData:
    """
    Read per-sample splice count parquet files and return a SplicingData object.

    Parameters
    ----------
    pattern : glob pattern or list of parquet paths.
        e.g. "counts/*.splice.parquet"
    sample_annotation : TSV/CSV path or DataFrame with 'sampleID' column.
        Extra columns forwarded to SplicingData.obs.
    sample_id_from_path : derive sampleID from filename stem (default True).
    gene_annotation : GTF path (.gtf / .gtf.gz), pre-built interval TSV, or
        DataFrame with [seqnames, start, end, gene_id, gene_name].
        Used both for gene labelling and (when gene_subset is given) for
        interval-based junction filtering.
    gene_subset : list of gene names or IDs to restrict the junction universe
        to. Any junction whose donor AND acceptor do not fall within one of
        these genes is discarded before building count matrices.
        Requires gene_annotation.
        Example: ["MCOLN1", "TIMMDC1", "CLPP"]
    n_jobs : parallel workers for file reading.
    verbose : progress messages.

    Returns
    -------
    SplicingData with:
        X                        = split counts  (samples x junctions)
        layers["site_counts"]    = site denominator counts
        layers["unsplit_counts"] = unsplit read counts
    """
    files = _resolve_files(pattern)
    if not files:
        raise FileNotFoundError(f"No parquet files matched: {pattern}")
    if verbose:
        print(f"[splice] Found {len(files)} count files")

    # Parse gene annotation once upfront if needed
    parsed_ann: Optional[pd.DataFrame] = None
    target_intervals: Optional[pd.DataFrame] = None

    if gene_annotation is not None or gene_subset is not None:
        if gene_subset is not None and gene_annotation is None:
            raise ValueError(
                "gene_annotation (GTF path or DataFrame) is required when "
                "gene_subset is provided."
            )
        if gene_annotation is not None:
            from .gtf import parse_gtf, genes_to_intervals
            parsed_ann = parse_gtf(gene_annotation) \
                if isinstance(gene_annotation, str) else gene_annotation.copy()
            if gene_subset is not None:
                target_intervals = genes_to_intervals(parsed_ann, gene_subset)
                if verbose:
                    print(f"[splice] Gene subset: {len(target_intervals):,} intervals "
                          f"from {len(gene_subset)} requested genes")

    sample_ids = _derive_sample_ids(files, sample_id_from_path)
    records    = _read_files_parallel(files, sample_ids, n_jobs=n_jobs, verbose=verbose)

    split_mat, site_mat, unsplit_mat, junction_df = _build_matrices(
        records, verbose=verbose
    )

    # Apply gene-subset filter BEFORE building obs/var
    if target_intervals is not None:
        from .gtf import filter_junctions_to_genes
        junction_df_filt = filter_junctions_to_genes(
            junction_df, target_intervals, verbose=verbose
        )
        keep_ids  = set(junction_df_filt["intron_id"])
        keep_mask = np.array([
            jid in keep_ids for jid in junction_df["intron_id"]
        ])
        split_mat   = split_mat[:,   keep_mask]
        site_mat    = site_mat[:,    keep_mask]
        unsplit_mat = unsplit_mat[:, keep_mask]
        junction_df = junction_df_filt
    elif parsed_ann is not None:
        # Gene annotation provided without subsetting — just label junctions
        junction_df = _annotate_genes(junction_df, parsed_ann)
    else:
        junction_df["gene_id"]   = junction_df["intron_id"]
        junction_df["gene_name"] = junction_df["intron_id"]

    # If gene subset was used, gene labels are already filled by filter_junctions_to_genes
    if target_intervals is not None and "gene_id" not in junction_df.columns:
        junction_df["gene_id"]   = junction_df["intron_id"]
        junction_df["gene_name"] = junction_df["intron_id"]

    obs = _build_obs(sample_ids, sample_annotation)

    var = junction_df.copy()
    var.index = var["intron_id"]
    var.index.name = None

    sd = SplicingData.from_matrices(split_mat, site_mat, obs=obs, var=var)
    sd.adata.layers["unsplit_counts"] = csr_matrix(unsplit_mat.astype(np.int32))

    return sd


def read_count_matrices(
    split_counts: Union[str, pd.DataFrame, np.ndarray],
    site_counts:  Union[str, pd.DataFrame, np.ndarray],
    unsplit_counts: Optional[Union[str, pd.DataFrame, np.ndarray]] = None,
    sample_annotation:   Optional[Union[str, pd.DataFrame]] = None,
    junction_annotation: Optional[Union[str, pd.DataFrame]] = None,
) -> SplicingData:
    """Build SplicingData from pre-built count matrices."""
    sc = _load_matrix(split_counts)
    dc = _load_matrix(site_counts)
    uc = _load_matrix(unsplit_counts) if unsplit_counts is not None else None

    if sc.shape != dc.shape:
        raise ValueError(f"split_counts {sc.shape} != site_counts {dc.shape}")
    if uc is not None and uc.shape != sc.shape:
        raise ValueError(f"unsplit_counts {uc.shape} != split_counts {sc.shape}")

    obs = _build_obs([f"sample_{i}" for i in range(sc.shape[0])], sample_annotation)

    if junction_annotation is not None:
        var = _load_df(junction_annotation)
        if "intron_id" not in var.columns:
            raise ValueError("junction_annotation must contain 'intron_id' column")
        var = var.set_index("intron_id")
    else:
        var = pd.DataFrame(
            {"intron_id": [f"intron_{j}" for j in range(sc.shape[1])]}
        ).set_index("intron_id")

    sd = SplicingData.from_matrices(sc, dc, obs=obs, var=var)
    if uc is not None:
        sd.adata.layers["unsplit_counts"] = csr_matrix(uc.astype(np.int32))
    return sd


# ==================================================================
# Internal helpers
# ==================================================================

def _resolve_files(pattern):
    if isinstance(pattern, (list, tuple)):
        return [str(p) for p in pattern]
    return sorted(glob.glob(pattern, recursive=True))


def _derive_sample_ids(files, from_path):
    ids = []
    for f in files:
        p = Path(f)
        if not from_path:
            ids.append(p.stem)
            continue
        name = p.name
        for suffix in (".splice.parquet", ".parquet"):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        ids.append(name or p.stem)

    seen: Dict[str, int] = {}
    out = []
    for sid in ids:
        if sid in seen:
            seen[sid] += 1
            out.append(f"{sid}_{seen[sid]}")
        else:
            seen[sid] = 0
            out.append(sid)
    return out


def _read_one_parquet(filepath, sample_id):
    try:
        df = pd.read_parquet(filepath)
    except Exception as exc:
        raise RuntimeError(f"Failed to read {filepath}: {exc}") from exc

    required = {"intron_id", "split_count", "site_count"}
    missing  = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{filepath} is missing required columns: {missing}. "
            "Was this file produced by splice.extract_sample_counts()?"
        )

    if "unsplit_count" not in df.columns:
        df["unsplit_count"] = np.int32(0)

    for col in ("seqnames", "start", "end", "strand"):
        if col not in df.columns:
            if col == "seqnames":
                df["seqnames"] = df["intron_id"].str.split(":").str[0]
            elif col == "strand":
                df["strand"] = df["intron_id"].str.split(":").str[-1]
            else:
                coords = df["intron_id"].str.split(":").str[1].str.split("-")
                df["start"] = coords.str[0].astype(np.int32)
                df["end"]   = coords.str[1].astype(np.int32)

    return sample_id, df


def _read_files_parallel(files, sample_ids, n_jobs, verbose):
    from joblib import Parallel, delayed
    pairs = list(zip(files, sample_ids))
    if n_jobs == 1:
        results = []
        for fp, sid in tqdm(pairs, desc="Reading count files", disable=not verbose):
            results.append(_read_one_parquet(fp, sid))
        return results
    return Parallel(n_jobs=n_jobs)(
        delayed(_read_one_parquet)(fp, sid)
        for fp, sid in tqdm(pairs, desc="Reading count files", disable=not verbose)
    )


def _build_matrices(records, verbose):
    if verbose:
        print("[splice] Building junction universe ...")

    meta_frames = []
    for _, df in records:
        cols = ["intron_id", "seqnames", "start", "end", "strand"]
        meta_frames.append(df[cols].drop_duplicates("intron_id"))

    junction_df = (
        pd.concat(meta_frames, ignore_index=True)
        .drop_duplicates("intron_id")
        .reset_index(drop=True)
    )
    junction_df["gene_id"]   = ""
    junction_df["gene_name"] = ""

    n_samples   = len(records)
    n_junctions = len(junction_df)
    junc_index  = {jid: i for i, jid in enumerate(junction_df["intron_id"])}

    if verbose:
        print(f"[splice] {n_junctions:,} junctions x {n_samples:,} samples")

    split_mat   = np.zeros((n_samples, n_junctions), dtype=np.int32)
    site_mat    = np.zeros((n_samples, n_junctions), dtype=np.int32)
    unsplit_mat = np.zeros((n_samples, n_junctions), dtype=np.int32)

    for s_idx, (_, df) in enumerate(
        tqdm(records, desc="Filling matrices", disable=not verbose)
    ):
        j_idxs = df["intron_id"].map(junc_index).dropna().astype(int)
        valid  = df.loc[j_idxs.index]
        jv     = j_idxs.values
        split_mat[s_idx, jv]   = valid["split_count"].values
        site_mat[s_idx, jv]    = valid["site_count"].values
        unsplit_mat[s_idx, jv] = valid["unsplit_count"].values

    return split_mat, site_mat, unsplit_mat, junction_df


def _build_obs(sample_ids, annotation):
    obs = pd.DataFrame({"sampleID": sample_ids})
    obs.index = obs["sampleID"]
    obs.index.name = None
    if annotation is not None:
        ann = _load_df(annotation)
        if "sampleID" not in ann.columns:
            raise ValueError("sample_annotation must contain a 'sampleID' column")
        obs = obs.merge(ann, on="sampleID", how="left")
        obs.index = obs["sampleID"]
        obs.index.name = None
    return obs


def _annotate_genes(junction_df, gene_annotation):
    from .gtf import parse_gtf
    ann = parse_gtf(gene_annotation) if isinstance(gene_annotation, str) \
          else gene_annotation.copy()
    if "gene_name" not in ann.columns:
        ann["gene_name"] = ann["gene_id"]

    merged = junction_df.merge(
        ann[["seqnames", "start", "end", "gene_id", "gene_name"]].rename(
            columns={"start": "g_start", "end": "g_end"}
        ),
        on="seqnames", how="left",
    )
    in_gene = (
        (merged["start"] >= merged["g_start"]) &
        (merged["end"]   <= merged["g_end"])
    )
    merged = merged[in_gene].drop_duplicates("intron_id")
    junction_df = junction_df.merge(
        merged[["intron_id", "gene_id", "gene_name"]], on="intron_id", how="left",
    )
    junction_df["gene_id"]   = junction_df["gene_id"].fillna(junction_df["intron_id"])
    junction_df["gene_name"] = junction_df["gene_name"].fillna(junction_df["intron_id"])
    return junction_df


def _load_matrix(src):
    if isinstance(src, np.ndarray):
        return src
    if isinstance(src, pd.DataFrame):
        return src.values
    return pd.read_csv(src, sep=None, engine="python", index_col=0).values.astype(np.float32)


def _load_df(src):
    if isinstance(src, pd.DataFrame):
        return src.copy()
    sep = "\t" if str(src).endswith((".tsv", ".tab")) else ","
    return pd.read_csv(src, sep=sep)
