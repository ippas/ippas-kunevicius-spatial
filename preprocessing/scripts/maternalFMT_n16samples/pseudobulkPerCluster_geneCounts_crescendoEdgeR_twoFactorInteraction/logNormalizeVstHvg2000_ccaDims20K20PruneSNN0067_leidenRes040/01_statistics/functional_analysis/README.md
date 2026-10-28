# Functional analysis – maternal FMT model (Visium, offspring brain)

## Run
```bash
pip install -r requirements.txt
python run_functional_analysis.py                      # input_DE_genes.tsv -> results/
python run_functional_analysis.py my_genes.tsv my_out  # own gene list / output folder
```
The first run downloads the MSigDB v2026.1 gene sets (~13 MB, from GitHub) to `msigdb_v2026.1_sysdata.rda`.
Runtime: about 1–2 minutes.

## Input
Tab-separated file with at least the columns `gene` (mouse symbol) and `log2FC`.

## Output (`results/`)
| File | Content |
|---|---|
| `functional_analysis_results.xlsx` | `summary` (one row per database), one sheet with all tested terms per database (lists all / up / down), `input_genes`, `column_description` |
| `summary.tsv` | the summary table as text |
| `overview_databases.png/.pdf` | −log10(FDR) of the top term per database |
| `GO_CC_significant_terms.png/.pdf` | significant GO Cellular Component terms |
| `run_info.txt` | date, parameters, gene-set file checksum, software versions |

Methods: see `METHODS.md`. The `results/` folder in this package is the run used for the presentation.
