"""
stats.py -- beta-binomial outlier scoring + delta PSI diagnostics
==================================================================
Adds three Sashimi-concordant delta columns alongside the existing
corrected delta_psi:

  delta_psi            = observed_J - expected_J_corrected  (existing, model-based)
  delta_psi_raw        = observed_J - cohort_median_J        (Sashimi-concordant)
  delta_psi_unweighted = observed_J - cohort_mean_J          (mean-based, stable at low N)
  correction_shrinkage_flag = True when |delta_psi_raw| / |delta_psi| > shrinkage_threshold
                              AND |delta_psi_raw| > 0.1
                              Signals the confounder correction has substantially
                              reduced the apparent effect size -- inspect Sashimi plot.

Per-sample gene subset FDR mirrors FRASER2 calculatePadjValuesOnSubset().
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.stats import betabinom
from tqdm import tqdm

from .data import SplicingData


# ==================================================================
# Public API
# ==================================================================

def call_outliers(
    sd: SplicingData,
    fdr_threshold: float = 0.05,
    delta_jaccard: float = 0.1,
    min_coverage: int = 10,
    use_corrected: bool = True,
    genes_per_sample: Optional[Dict[str, List[str]]] = None,
    shrinkage_threshold: float = 3.0,
    n_jobs: int = 1,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Identify aberrant splicing outliers via beta-binomial testing.

    Parameters
    ----------
    sd : SplicingData after compute_jaccard() and correct_confounders().
    fdr_threshold : BH FDR cutoff (default 0.05).
    delta_jaccard : minimum |delta_psi| effect size filter (default 0.1).
        Applied to the corrected delta_psi. Events with |delta_psi_raw| >= 0.1
        but small corrected delta_psi will still appear if correction_shrinkage_flag
        is True -- see shrinkage_threshold.
    min_coverage : minimum site read coverage to report a call.
    use_corrected : use confounder-corrected expected J (default True).
    genes_per_sample : dict mapping sampleID -> list of gene names/IDs for
        per-sample gene subset FDR (mirrors FRASER2 calculatePadjValuesOnSubset).
    shrinkage_threshold : ratio |delta_psi_raw| / |delta_psi| above which
        correction_shrinkage_flag is set to True (default 3.0). Also requires
        |delta_psi_raw| > 0.1 to avoid flagging noise.
    n_jobs : parallel workers.
    verbose : progress messages.

    Returns
    -------
    pd.DataFrame with FRASER2-compatible columns plus diagnostics:

        sampleID, seqnames, start, end, strand,
        intron_id, gene_id, gene_name,
        pvalue, padj,
        delta_psi              -- observed_J minus model-expected_J (corrected)
        delta_psi_raw          -- observed_J minus cohort median_J (Sashimi-concordant)
        delta_psi_unweighted   -- observed_J minus cohort mean_J
        observed_psi, expected_psi,
        counts, total_counts,
        fdr_subset,
        correction_shrinkage_flag -- True: correction may be masking a real event;
                                     inspect Sashimi plot regardless of padj
    """
    if "jaccard" not in sd.adata.layers:
        raise RuntimeError("Run splice.compute_jaccard(sd) before call_outliers().")

    J_raw = sd.jaccard()
    split = sd.split_counts()
    site  = sd.site_counts()

    # Corrected expected J
    if use_corrected and "jaccard_corrected" in sd.adata.layers:
        J_corrected = sd.adata.layers["jaccard_corrected"]
        J_logit = sd.adata.layers.get("jaccard_logit")
        if J_logit is None:
            eps = 1e-6
            J_clipped = np.clip(J_raw, eps, 1 - eps)
            J_logit = np.log(J_clipped / (1 - J_clipped))
        J_expected = _sigmoid(J_logit - J_corrected)
    else:
        J_expected = None

    # Cohort median and mean per junction (S,J) -> (J,) for raw deltas
    # Computed only over samples with non-zero site coverage
    safe_site  = np.where(site > 0, site, np.nan)
    J_for_ref  = np.where(site > 0, J_raw, np.nan)
    J_median   = np.nanmedian(J_for_ref, axis=0).astype(np.float32)   # (J,)
    J_mean     = np.nanmean(  J_for_ref, axis=0).astype(np.float32)   # (J,)

    n_samples, n_junctions = J_raw.shape
    sample_ids = list(sd.samples)
    var = sd.var.copy()

    if verbose:
        print(f"[splice] Fitting beta-binomial for {n_junctions:,} junctions ...")

    alpha, beta_param = _fit_bb_mom_vectorised(split, site)

    rows = _score_all_junctions(
        split, site, alpha, beta_param,
        J_raw, J_expected, J_median, J_mean,
        sample_ids, var, min_coverage, n_jobs, verbose,
    )

    if not rows:
        return _empty_result_df()

    results = pd.DataFrame(rows)
    results["fdr_subset"] = False

    # FDR correction
    if genes_per_sample:
        results = _apply_fdr_with_subsets(results, fdr_threshold, genes_per_sample)
    else:
        results = _apply_fdr_global(results, fdr_threshold)

    # correction_shrinkage_flag: model underestimates effect vs raw signal
    results["correction_shrinkage_flag"] = (
        (results["delta_psi_raw"].abs() > 0.1) &
        (results["delta_psi"].abs() > 0) &
        (results["delta_psi_raw"].abs() / results["delta_psi"].abs().clip(lower=1e-6)
         > shrinkage_threshold)
    )

    # Effect size filter on corrected delta — but keep shrinkage-flagged rows
    # even if their corrected delta_psi is small, since they may be real events
    # that the correction absorbed
    pass_effect = results["delta_psi"].abs() >= delta_jaccard
    pass_shrink = results["correction_shrinkage_flag"]
    results = results[pass_effect | pass_shrink].copy()

    # Sort: padj asc, |delta_psi_raw| desc (raw is Sashimi-concordant)
    results = results.sort_values(
        ["padj", "delta_psi_raw"],
        ascending=[True, False],
        key=lambda c: c.abs() if c.name == "delta_psi_raw" else c,
    ).reset_index(drop=True)

    if verbose:
        n_sig     = (results["padj"] < fdr_threshold).sum()
        n_shrink  = results["correction_shrinkage_flag"].sum()
        n_subset  = results["fdr_subset"].sum()
        print(
            f"[splice] call_outliers: {len(results):,} events "
            f"({n_sig:,} padj < {fdr_threshold}, "
            f"{n_shrink:,} correction_shrinkage_flag, "
            f"{n_subset:,} per-sample subset FDR)"
        )

    return results


# ==================================================================
# FDR correction
# ==================================================================

def _apply_fdr_global(results, fdr_threshold):
    from statsmodels.stats.multitest import multipletests
    results = results.copy()
    results["padj"] = np.nan
    for gene, grp in results.groupby("gene_id"):
        pvals = grp["pvalue"].values
        if len(pvals) == 0:
            continue
        try:
            _, padj, _, _ = multipletests(pvals, method="fdr_bh", alpha=fdr_threshold)
        except Exception:
            padj = pvals
        results.loc[grp.index, "padj"] = padj
    return results


def _apply_fdr_with_subsets(results, fdr_threshold, genes_per_sample):
    from statsmodels.stats.multitest import multipletests
    results = results.copy()
    results["padj"] = np.nan

    subset_lookup = {
        sid: frozenset(g.lower() for g in genes)
        for sid, genes in genes_per_sample.items()
    }
    subset_sample_mask = results["sampleID"].isin(subset_lookup)

    # Per-sample subset FDR
    subset_rows = results[subset_sample_mask].copy()
    if len(subset_rows) > 0:
        def _in_subset(row):
            gene_set = subset_lookup.get(row["sampleID"], frozenset())
            return (
                str(row.get("gene_name", "")).lower() in gene_set or
                str(row.get("gene_id",   "")).lower() in gene_set
            )
        in_subset_mask = subset_rows.apply(_in_subset, axis=1)
        results.loc[subset_rows.index[in_subset_mask], "fdr_subset"] = True
        testable = subset_rows[in_subset_mask].copy()
        for (sid, gene), grp in testable.groupby(["sampleID", "gene_id"]):
            pvals = grp["pvalue"].values
            if len(pvals) == 0:
                continue
            try:
                _, padj, _, _ = multipletests(pvals, method="fdr_bh", alpha=fdr_threshold)
            except Exception:
                padj = pvals
            results.loc[grp.index, "padj"] = padj

    # Standard global FDR for non-subset samples
    for gene, grp in results[~subset_sample_mask].groupby("gene_id"):
        pvals = grp["pvalue"].values
        if len(pvals) == 0:
            continue
        try:
            _, padj, _, _ = multipletests(pvals, method="fdr_bh", alpha=fdr_threshold)
        except Exception:
            padj = pvals
        results.loc[grp.index, "padj"] = padj

    return results


# ==================================================================
# BB fitting (method of moments, vectorised)
# ==================================================================

def _fit_bb_mom_vectorised(split, site):
    safe_site = np.where(site <= 0, 1, site)
    psi = split / safe_site
    mu    = psi.mean(axis=0)
    var   = psi.var(axis=0)
    n_bar = safe_site.mean(axis=0)
    eps      = 1e-9
    mu_var   = np.clip(mu * (1 - mu), eps, None)
    bin_var  = mu_var / np.maximum(n_bar, 1)
    excess   = var - bin_var
    phi      = np.where(excess > 0, excess / mu_var, 1e-4)
    phi      = np.clip(phi, 1e-6, 1 - 1e-6)
    conc     = (1 / phi) - 1
    alpha    = np.clip(mu * conc,       1e-4, 1e4).astype(np.float32)
    beta_p   = np.clip((1 - mu) * conc, 1e-4, 1e4).astype(np.float32)
    return alpha, beta_p


# ==================================================================
# p-value computation
# ==================================================================

def _score_all_junctions(
    split, site, alpha, beta_param,
    J_raw, J_expected, J_median, J_mean,
    sample_ids, var, min_coverage, n_jobs, verbose,
):
    n_samples, n_junctions = split.shape

    def _score_junction(j):
        a = float(alpha[j])
        b = float(beta_param[j])
        n_j   = site[:, j]
        k_j   = split[:, j]
        J_j   = J_raw[:, j]
        j_med = float(J_median[j])
        j_mn  = float(J_mean[j])

        if J_expected is not None:
            J_exp_j = J_expected[:, j]
        else:
            J_exp_j = np.full(n_samples, a / (a + b), dtype=np.float32)

        junc_rows = []
        for s in range(n_samples):
            n_s = int(n_j[s])
            if n_s < min_coverage:
                continue
            k_s = int(k_j[s])
            try:
                p_lower = betabinom.cdf(k_s, n_s, a, b)
                p_upper = 1.0 - betabinom.cdf(k_s - 1, n_s, a, b)
                pval = float(min(2 * min(p_lower, p_upper), 1.0))
            except Exception:
                pval = 1.0
            if pval >= 1.0:
                continue

            obs_j = float(J_j[s])
            exp_j = float(J_exp_j[s])

            junc_rows.append({
                "sampleID":            sample_ids[s],
                "intron_id":           var.index[j] if hasattr(var, "index") else f"intron_{j}",
                "pvalue":              pval,
                "delta_psi":           obs_j - exp_j,
                "delta_psi_raw":       obs_j - j_med,
                "delta_psi_unweighted": obs_j - j_mn,
                "observed_psi":        obs_j,
                "expected_psi":        exp_j,
                "counts":              k_s,
                "total_counts":        n_s,
            })
        return junc_rows

    if n_jobs == 1:
        all_rows = []
        for j in tqdm(range(n_junctions), desc="Scoring junctions", disable=not verbose):
            all_rows.extend(_score_junction(j))
    else:
        from joblib import Parallel, delayed
        nested = Parallel(n_jobs=n_jobs)(
            delayed(_score_junction)(j)
            for j in tqdm(range(n_junctions), desc="Scoring junctions", disable=not verbose)
        )
        all_rows = [row for chunk in nested for row in chunk]

    if not all_rows:
        return []

    df_rows   = pd.DataFrame(all_rows)
    var_reset = (
        var.reset_index(drop=True) if "intron_id" in var.columns
        else var.reset_index().rename(columns={"index": "intron_id"})
    )
    df_rows = df_rows.merge(
        var_reset[["intron_id", "seqnames", "start", "end", "strand",
                   "gene_id", "gene_name"]],
        on="intron_id", how="left",
    )
    return df_rows.to_dict("records")


def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def _empty_result_df():
    return pd.DataFrame(columns=[
        "sampleID", "seqnames", "start", "end", "strand",
        "intron_id", "gene_id", "gene_name",
        "pvalue", "padj",
        "delta_psi", "delta_psi_raw", "delta_psi_unweighted",
        "observed_psi", "expected_psi",
        "counts", "total_counts",
        "fdr_subset", "correction_shrinkage_flag",
    ])
