from __future__ import annotations
import numpy as np
from .data import SplicingData

def compute_jaccard(
    sd: SplicingData,
    store_logit: bool = True,
    verbose: bool = True,
) -> SplicingData:
    """
    Computes the intron Jaccard index matrix across the cohort.
    J = split_reads / (site_totals + unsplit_reads)
    """
    if verbose:
        print("[splice] Computing intron Jaccard index matrix...")

    split = sd.split_counts()
    site = sd.site_counts()
    
    if "unsplit_counts" in sd.adata.layers:
        from scipy.sparse import issparse
        uc = sd.adata.layers["unsplit_counts"]
        unsplit = uc.toarray() if issparse(uc) else np.asarray(uc)
    else:
        unsplit = np.zeros_like(split)

    denom = site + unsplit
    safe_denom = np.where(denom <= 0, 1.0, denom)
    jaccard = np.clip(split / safe_denom, 0.0, 1.0).astype(np.float32)
    
    sd.adata.layers["jaccard"] = jaccard

    if store_logit:
        eps = 1e-6
        j_clipped = np.clip(jaccard, eps, 1.0 - eps)
        sd.adata.layers["jaccard_logit"] = np.log(j_clipped / (1.0 - j_clipped)).astype(np.float32)

    return sd
