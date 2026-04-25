# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [0.2.0] — 2025

### Added
- `extract_sample_counts()` — BAM → per-sample parquet extraction (replaces STAR SJ.out.tab input).
  Two-pass scan: split reads from N-CIGAR operations, unsplit reads at splice sites.
- `min_anchor` parameter in `extract_sample_counts()` — matches FRASER2's `minAnchor` filter.
  Default 5 bp. Unsplit reads must extend at least `min_anchor` bases past the splice site
  on both sides to be counted.
- `genes_per_sample` parameter in `call_outliers()` and `find_outliers()` — per-sample gene
  subset FDR correction, mirroring FRASER2's `calculatePadjValuesOnSubset()`. Clinically
  useful for focusing statistical power on candidate genes from exome/genome sequencing.
  Gene matching is case-insensitive. Results include a `fdr_subset` column.
- Full intron Jaccard index computation using BAM-extracted unsplit read counts.
  Falls back to PSI approximation when unsplit counts are unavailable.
- `SplicingData.save()` / `SplicingData.load()` for persisting the intermediate object
  as `.h5ad` (AnnData HDF5 format).

### Changed
- Primary input changed from STAR SJ.out.tab files to BAM files via the new
  `extract_sample_counts()` step.
- `read_star()` renamed to `read_counts()` — reads `.splice.parquet` files.
- `find_outliers()` parameter renamed from `star_files` to `count_files`.
- `setup.py` replaced by `pyproject.toml`.

### Removed
- Direct STAR SJ.out.tab reading (superseded by BAM extraction).

## [0.1.0] — 2025

### Added
- Initial release.
- `read_star()` — reads STAR SJ.out.tab files.
- `filter_junctions()` — FRASER2-compatible k/n/q filters.
- `compute_jaccard()` — intron Jaccard index.
- `correct_confounders()` — randomised SVD confounder correction with
  Marchenko-Pastur automatic component selection.
- `call_outliers()` — vectorised beta-binomial MLE + BH FDR.
- `find_outliers()` — single-function pipeline.
- `SplicingData` — AnnData wrapper with FRASER2-compatible metadata.
- FRASER2-compatible output column names.
