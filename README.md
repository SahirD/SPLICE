# NOTE: NOT ready for use

# splice-outlier

Fast, scalable detection of aberrant splicing events in RNA-seq cohorts.
A drop-in replacement for [FRASER2](https://github.com/gagneurlab/FRASER) designed for 10,000+ samples.

## Key differences from FRASER2

| | FRASER2 | splice-outlier |
|---|---|---|
| Language | R / Bioconductor | Python |
| Confounder correction | Denoising autoencoder (O(S²)) | Randomised truncated SVD (O(S·J·q)) |
| Input | BAM + STAR SJ files | BAM only (parquet intermediary) |
| Scale tested | ~500 samples | 10,000+ samples |
| Output format | FraserDataSet (R S4) | pandas DataFrame (FRASER2-compatible columns) |
| Per-sample gene subset FDR | ✓ | ✓ |

## Installation

```bash
# Core package (no BAM support)
pip install splice-outlier

# With BAM extraction support (requires pysam)
pip install "splice-outlier[bam]"

# Development install
git clone https://github.com/your-org/splice-outlier
cd splice-outlier
pip install -e ".[dev]"
```

## Quick start

### Step 1 — extract per-sample counts from BAM files

Run once per BAM, ideally on the HPC node where the BAM lives.
The BAM is not needed again after this step.

```python
import splice

splice.extract_sample_counts(
    "data/sample1/Aligned.sortedByCoord.out.bam",
    "counts/sample1.splice.parquet",
    library_type="fr-firststrand",   # or fr-secondstrand / unstranded
    min_anchor=5,                    # matches FRASER2 default
    min_mapping_quality=255,         # STAR unique mappers (use 60 for HISAT2)
)
```

Parallelise across a cohort:

```python
from joblib import Parallel, delayed

bam_files = [...]   # list of BAM paths

Parallel(n_jobs=32)(
    delayed(splice.extract_sample_counts)(
        bam,
        bam.replace(".bam", ".splice.parquet"),
        library_type="fr-firststrand",
    )
    for bam in bam_files
)
```

### Step 2 — detect outliers across the cohort

```python
results = splice.find_outliers(
    "counts/*.splice.parquet",
    sample_annotation="samples.tsv",   # TSV/CSV with 'sampleID' column
    gene_annotation="hg38_genes.tsv",  # optional: seqnames, start, end, gene_id, gene_name
    fdr_threshold=0.05,
    delta_jaccard=0.1,                 # FRASER2 default
    n_jobs=8,
)

results.to_csv("splicing_outliers.csv", index=False)
```

### Clinical use — per-sample gene subset FDR

Focus statistical power on candidate genes from exome/genome sequencing,
mirroring FRASER2's `calculatePadjValuesOnSubset()`:

```python
results = splice.find_outliers(
    "counts/*.splice.parquet",
    sample_annotation="samples.tsv",
    gene_annotation="hg38_genes.tsv",
    genes_per_sample={
        "patient_01": ["MCOLN1", "TIMMDC1"],
        "patient_02": ["CLPP"],
        "patient_03": ["TTN", "RYR1"],
    },
    fdr_threshold=0.05,
)

# Rows scored with per-sample gene subset FDR are flagged
subset_hits = results[results["fdr_subset"] & (results["padj"] < 0.05)]
```

## Step-by-step API

For more control, run each stage individually:

```python
import splice

# Read parquet files
sd = splice.read_counts(
    "counts/*.splice.parquet",
    sample_annotation="samples.tsv",
    gene_annotation="hg38_genes.tsv",
)

# Filter low-coverage junctions (FRASER2 defaults)
sd = splice.filter_junctions(sd, min_split_reads=20, min_site_reads=10)

# Compute intron Jaccard index
sd = splice.compute_jaccard(sd)

# Remove latent confounders (randomised SVD)
sd = splice.correct_confounders(sd, n_components=None)  # auto-select via Marchenko-Pastur

# Save intermediate object (skip expensive correction on re-runs)
sd.save("cohort_corrected.h5ad")

# Call outliers
results = splice.call_outliers(
    sd,
    fdr_threshold=0.05,
    delta_jaccard=0.1,
    genes_per_sample={"patient_01": ["MCOLN1"]},
)
```

## Output format

The results DataFrame mirrors FRASER2's output columns exactly so existing
downstream scripts, clinical report templates, and DROP pipeline steps
require no modification:

| Column | Description |
|---|---|
| `sampleID` | Sample identifier |
| `seqnames` | Chromosome |
| `start` | Intron start (1-based) |
| `end` | Intron end (1-based) |
| `strand` | `+`, `-`, or `*` |
| `intron_id` | `seqnames:start-end:strand` |
| `gene_id` | Ensembl gene ID |
| `gene_name` | HGNC symbol |
| `pvalue` | Two-sided beta-binomial p-value |
| `padj` | Benjamini-Hochberg FDR within gene |
| `delta_psi` | Observed minus expected Jaccard index |
| `observed_psi` | Observed Jaccard index |
| `expected_psi` | Model-expected Jaccard index |
| `counts` | Split reads for this junction |
| `total_counts` | Total reads at the splice site |
| `fdr_subset` | `True` if per-sample gene subset FDR was applied |

## Parameters

### `extract_sample_counts`

| Parameter | Default | Description |
|---|---|---|
| `library_type` | `fr-firststrand` | Strand protocol for unsplit read inference |
| `min_mapping_quality` | `255` | Minimum MAPQ (255 = STAR unique, 60 = HISAT2 unique) |
| `min_split_reads` | `2` | Minimum split reads to retain a junction |
| `min_anchor` | `5` | Minimum bases a non-split read must extend past each side of a splice site (matches FRASER2 `minAnchor`) |

### `find_outliers`

| Parameter | Default | Description |
|---|---|---|
| `min_split_reads` | `20` | Junction filter: minimum split reads in ≥1 sample |
| `min_site_reads` | `10` | Junction filter: minimum site count at 25th percentile |
| `n_components` | `None` | PCs to remove; `None` = auto via Marchenko-Pastur |
| `fdr_threshold` | `0.05` | BH FDR cutoff |
| `delta_jaccard` | `0.1` | Minimum \|ΔJ\| effect size |
| `genes_per_sample` | `None` | Dict of `sampleID → [gene_names]` for clinical subset FDR |
| `n_jobs` | `1` | Parallel workers |

## How it works

1. **BAM extraction** — two-pass scan: pass 1 discovers junctions from N-CIGAR reads and counts split reads; pass 2 counts non-split reads at each splice site with configurable anchor filtering.

2. **Jaccard index** — computes `J = s_canon / (s_canon + s_alt + u_unsplit)` per (sample, junction). With unsplit counts from the BAM this is the full intron Jaccard index as defined in FRASER2. Without them it reduces to a PSI-like approximation.

3. **Confounder correction** — replaces FRASER2's denoising autoencoder with randomised truncated SVD. The top-q principal components of the logit-Jaccard matrix are projected out, where q is auto-selected via the Marchenko-Pastur noise threshold. This is 10–100× faster at large cohort sizes.

4. **Outlier scoring** — beta-binomial parameters are fit per junction via vectorised method of moments. Two-sided p-values are computed and BH FDR is applied within each gene (or within per-sample gene subsets).

## Development

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

## Citation

If you use splice-outlier in published work, please also cite FRASER2 which this tool is based on:

> Scheller et al. (2023) Improved detection of aberrant splicing with FRASER 2.0 and the intron Jaccard index. *Am J Hum Genet* 110(12):2056–2067. https://doi.org/10.1016/j.ajhg.2023.10.014

