from __future__ import annotations
from typing import Optional
import numpy as np
from .data import SplicingData

def filter_junctions(
    sd: SplicingData,
    min_split_reads: int = 20,
    min_site_reads: int = 10,
    site_count_percentile: float = 25.0,
    verbose: bool = True,
) -> SplicingData:
    """
    Filters out low-coverage junctions across the cohort mirroring FRASER2 defaults.
    Returns a new SplicingData instance containing the filtered subset.
    """
    n_before = sd.n_junctions
    split = sd.split_counts()
    site = sd.site_counts()

    max_split = np.max(split, axis=0)
    pass_split = max_split >= min_split_reads

    perc_site = np.percentile(site, site_count_percentile, axis=0)
    pass_site = perc_site >= min_site_reads

    keep_mask = pass_split & pass_site
    
    new_adata = sd.adata[:, keep_mask].copy()
    
    new_adata.uns["filter_params"] = {
        "min_split_reads": min_split_reads,
        "min_site_reads": min_site_reads,
        "site_count_percentile": site_count_percentile,
    }

    if verbose:
        n_after = new_adata.shape[1]
        print(f"[splice] Filtered junctions: {n_before:,} -> {n_after:,} ({n_before - n_after:,} removed)")

    return SplicingData(new_adata)

def correct_confounders(
    sd: SplicingData,
    n_components: Optional[int] = None,
    verbose: bool = True,
) -> SplicingData:
    """
    Removes latent experimental confounders using randomized truncated SVD.
    Bypasses correction safely if the matrix has insufficient dimensions.
    """
    from sklearn.utils.extmath import randomized_svd

    if "jaccard_logit" not in sd.adata.layers:
        raise RuntimeError("Run compute_jaccard(sd, store_logit=True) before correcting confounders.")

    logit_mat = sd.adata.layers["jaccard_logit"]
    n_samples, n_features = logit_mat.shape

    # SAFETY GUARDRAIL: LAPACK SGETRF crashes if dimensions are 0.
    # SVD also requires at least 2 dimensions to project meaningful components.
    if n_samples < 2 or n_features < 2:
        if verbose:
            print(f"[splice] Skipping SVD correction: insufficient matrix dimensions ({n_samples}x{n_features}).")
        # Populate the corrected layer with uncorrected residuals to allow downstream execution
        sd.adata.layers["jaccard_corrected"] = np.zeros_like(logit_mat)
        sd.adata.uns["correction_params"] = {"n_components": 0}
        return sd

    feature_means = np.mean(logit_mat, axis=0)
    centered_mat = logit_mat - feature_means

    if n_components is None:
        beta = n_samples / n_features if n_samples < n_features else n_features / n_samples
        omega = 0.56 * beta**3 - 0.95 * beta**2 + 1.82 * beta + 1.43
        y = randomized_svd(centered_mat, n_components=min(n_samples, 100), n_iter=3, random_state=42)[1]
        cutoff = omega * np.median(y)
        n_components = int(np.sum(y > cutoff))
        n_components = np.clip(n_components, 2, 100)
        # Ensure we don't ask for more components than the matrix rank supports
        max_rank = min(n_samples, n_features)
        n_components = min(n_components, max_rank)
        if verbose:
            print(f"[splice] Marchenko-Pastur auto-selected components: {n_components}")

    if verbose:
        print(f"[splice] Correcting confounders via randomized SVD (q={n_components})...")

    U, Sigma, VT = randomized_svd(centered_mat, n_components=n_components, n_iter=5, random_state=42)
    
    reconstructed = np.dot(U, np.dot(np.diag(Sigma), VT))
    corrected_residuals = centered_mat - reconstructed

    sd.adata.layers["jaccard_corrected"] = corrected_residuals
    sd.adata.uns["correction_params"] = {"n_components": n_components}

    return sd
