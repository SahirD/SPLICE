# Contributing

## Development setup

```bash
git clone https://github.com/your-org/splice-outlier
cd splice-outlier
pip install -e ".[dev]"
```

## Running tests

```bash
pytest tests/ -v
pytest tests/ -v --cov=splice --cov-report=term-missing
```

## Code structure

| Module | Responsibility |
|---|---|
| `extract.py` | BAM → parquet extraction (`extract_sample_counts`) |
| `io.py` | Parquet → `SplicingData` (`read_counts`) |
| `data.py` | `SplicingData` container (AnnData wrapper) |
| `metrics.py` | Intron Jaccard index computation |
| `correction.py` | Junction filtering + randomised SVD confounder correction |
| `stats.py` | Beta-binomial scoring + FDR correction |
| `pipeline.py` | `find_outliers()` one-call interface |

## Submitting changes

1. Fork the repository and create a feature branch.
2. Add tests for any new functionality — aim for all edge cases covered.
3. Ensure `pytest tests/ -v` passes with no failures.
4. Open a pull request with a clear description of the change and its motivation.

## Reporting issues

Please include:
- Python version and OS
- Package version (`import splice; print(splice.__version__)`)
- Minimal reproducible example
- Full traceback
