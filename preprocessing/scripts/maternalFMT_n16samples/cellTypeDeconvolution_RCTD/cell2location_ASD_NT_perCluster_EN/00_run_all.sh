#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_ROOT="${PROJECT_ROOT:-/home/mateusz/mnt/ippas-kunevicius-spatial}"
if [[ ! -d "$PROJECT_ROOT" ]]; then echo "ERROR: Missing project: $PROJECT_ROOT" >&2; exit 1; fi
VENV="$PROJECT_ROOT/.venv-cell2location"
if [[ ! -x "$VENV/bin/python" ]]; then python3 -m venv "$VENV"; fi
"$VENV/bin/python" -m pip install --upgrade pip
"$VENV/bin/python" -m pip install 'cell2location>=0.1.4' 'anndata' 'scanpy' 'pandas' 'numpy' 'scipy' 'statsmodels' 'torch'
DATA="$PROJECT_ROOT/results/maternalFMT_n16samples/cellTypeDeconvolution_cell2location/mouseBrain59/inputs"
Rscript -e 'for(p in c("Seurat","SeuratObject","Matrix","data.table","R.utils")){if(!requireNamespace(p,quietly=TRUE))stop("Missing R package: ",p)}'
Rscript "$SCRIPT_DIR/01_export_visium.R" "$PROJECT_ROOT" "$DATA"
"$VENV/bin/python" "$SCRIPT_DIR/02_cell2location_full.py"
"$VENV/bin/python" "$SCRIPT_DIR/03_compare_ASD_vs_Neurotypical.py"
