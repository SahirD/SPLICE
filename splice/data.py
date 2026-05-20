from __future__ import annotations
from typing import Optional, Union
import numpy as np
import pandas as pd

class SplicingData:
    """
    Core data container wrapping an AnnData object for splicing cohorts.
    Stores split read counts (X), site totals, and unsplit read layers.
    """
    def __init__(self, adata):
        self.adata = adata

    @classmethod
    def from_matrices(
        cls,
        split_counts: np.ndarray,
        site_counts: np.ndarray,
        obs: pd.DataFrame,
        var: pd.DataFrame,
    ) -> SplicingData:
        X = split_counts.astype(np.float32)
        layers = {"site_counts": site_counts.astype(np.float32)}
        
        import anndata as ad
        adata = ad.AnnData(X=X, obs=obs, var=var, layers=layers)
        return cls(adata)

    @property
    def n_samples(self) -> int:
        return self.adata.shape[0]

    @property
    def n_junctions(self) -> int:
        return self.adata.shape[1]

    @property
    def samples(self) -> pd.Index:
        return self.adata.obs.index

    @property
    def obs(self) -> pd.DataFrame:
        return self.adata.obs

    @property
    def var(self) -> pd.DataFrame:
        return self.adata.var

    def split_counts(self) -> np.ndarray:
        return self.adata.X

    def site_counts(self) -> np.ndarray:
        return self.adata.layers["site_counts"]

    def jaccard(self) -> np.ndarray:
        if "jaccard" not in self.adata.layers:
            raise KeyError("Jaccard matrix not found. Run compute_jaccard() first.")
        return self.adata.layers["jaccard"]

    def jaccard_corrected(self) -> np.ndarray:
        if "jaccard_corrected" not in self.adata.layers:
            raise KeyError("Corrected matrix not found. Run correct_confounders() first.")
        return self.adata.layers["jaccard_corrected"]

    def save(self, filename: str) -> None:
        # Ensure all metadata columns are strings to prevent H5PY export errors
        for col in self.adata.obs.columns:
            self.adata.obs[col] = self.adata.obs[col].astype(str)
        for col in self.adata.var.columns:
            self.adata.var[col] = self.adata.var[col].astype(str)
            
        self.adata.write_h5ad(filename)

    @classmethod
    def load(cls, filename: str) -> SplicingData:
        import anndata as ad
        return cls(ad.read_h5ad(filename))

    def __repr__(self) -> str:
        return f"SplicingData object with {self.n_samples} samples and {self.n_junctions} junctions."
