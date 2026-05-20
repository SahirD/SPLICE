"""
pipeline.py -- one-call interface
"""

from __future__ import annotations

from typing import Dict, List, Optional, Union

import pandas as pd

from .io import read_counts
from .correction import filter_junctions, correct_confounders
from .metrics import compute_jaccard
from .stats import call_outliers
from .data import SplicingData


def find_outliers(
    count_files: Union[str, List[str]],
    sample_annotation: Optional[Union[str, pd.DataFrame]] = None,
    gene_annotation: Optional[Union[str, pd.DataFrame]] = None,
    gene_subset: Optional[List[str]] = None,
    # Filtering
    min_split_reads: int = 20,
    min_site_reads: int = 10,
    site_count_percentile: float = 25.0,
    # Correction
    n_components: Optional[int] = None,
    # Calling
    fdr_threshold: float = 0.05,
    delta_jaccard: float = 0.1,
    min_coverage: int = 10,
    genes_per_sample: Optional[Dict[str, List[str]]] = None,
    shrinkage_threshold: float = 3.0,
    # I/O
    n_jobs: int = 1,
    save_intermediate: Optional[str] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Detect aberrant splicing outliers from per-sample splice count parquet files.

    Parameters
    ----------
    count_files : glob pattern or list of *.splice.parquet paths.
    sample_annotation : TSV/CSV or DataFrame with 'sampleID' column.
    gene_annotation : GTF path (.gtf/.gtf.gz), pre-built interval TSV, or
        DataFrame with [seqnames, start, end, gene_id, gene_name].
        Used for junction gene labelling and (with gene_subset) interval filtering.
    gene_subset : list of gene names or IDs to restrict the entire analysis to.
        Junctions outside these genes are dropped before any modelling.
        Dramatically reduces runtime for targeted panels.
        Requires gene_annotation.
        Example: ["MCOLN1", "TIMMDC1", "CLPP"]
    min_split_reads : junction filter k (FRASER2 default 20).
    min_site_reads : junction filter n at 25th percentile (FRASER2 default 10).
    site_count_percentile : percentile for min_site_reads (default 25).
    n_components : PCs for confounder correction. None = auto (Marchenko-Pastur).
    fdr_threshold : BH FDR cutoff (default 0.05).
    delta_jaccard : minimum |delta_psi| effect size (default 0.1).
        Also keeps events with correction_shrinkage_flag=True regardless of
        corrected delta_psi, since those may be real events the model absorbed.
    min_coverage : minimum site read coverage to report a call.
    genes_per_sample : dict sampleID -> [gene_names] for per-sample gene subset
        FDR correction (mirrors FRASER2 calculatePadjValuesOnSubset).
    shrinkage_threshold : ratio |delta_psi_raw|/|delta_psi| above which
        correction_shrinkage_flag is set (default 3.0). Flags events where the
        confounder correction has substantially reduced apparent effect size
        vs what a Sashimi plot would show.
    n_jobs : parallel workers.
    save_intermediate : .h5ad path to persist SplicingData after correction.
    verbose : progress messages.

    Returns
    -------
    pd.DataFrame with FRASER2-compatible columns plus diagnostics:

        sampleID, seqnames, start, end, strand,
        intron_id, gene_id, gene_name,
        pvalue, padj,
        delta_psi              -- observed_J minus model-expected (corrected)
        delta_psi_raw          -- observed_J minus cohort median  (Sashimi-concordant)
        delta_psi_unweighted   -- observed_J minus cohort mean
        observed_psi, expected_psi, counts, total_counts,
        fdr_subset,
        correction_shrinkage_flag

    Examples
    --------
    Full transcriptome:
    >>> results = splice.find_outliers("counts/*.splice.parquet")

    Targeted panel (faster, focused):
    >>> results = splice.find_outliers(
    ...     "counts/*.splice.parquet",
    ...     gene_annotation="hg38.gtf.gz",
    ...     gene_subset=["MCOLN1", "TIMMDC1", "CLPP"],
    ...     genes_per_sample={"patient_01": ["MCOLN1"]},
    ... )
    >>> # Events where Sashimi would show a bigger difference than padj suggests
    >>> flagged = results[results["correction_shrinkage_flag"]]
    """
    sd = read_counts(
        count_files,
        sample_annotation=sample_annotation,
        gene_annotation=gene_annotation,
        gene_subset=gene_subset,
        n_jobs=n_jobs,
        verbose=verbose,
    )

    sd = filter_junctions(
        sd,
        min_split_reads=min_split_reads,
        min_site_reads=min_site_reads,
        site_count_percentile=site_count_percentile,
        verbose=verbose,
    )

    sd = compute_jaccard(sd, verbose=verbose)
    sd = correct_confounders(sd, n_components=n_components, verbose=verbose)

    if save_intermediate:
        if verbose:
            print(f"[splice] Saving SplicingData -> {save_intermediate}")
        sd.save(save_intermediate)

    return call_outliers(
        sd,
        fdr_threshold=fdr_threshold,
        delta_jaccard=delta_jaccard,
        min_coverage=min_coverage,
        genes_per_sample=genes_per_sample,
        shrinkage_threshold=shrinkage_threshold,
        n_jobs=n_jobs,
        verbose=verbose,
    )
