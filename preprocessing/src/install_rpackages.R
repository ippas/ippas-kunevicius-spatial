# ==============================================================================
# Install required R packages for spatial transcriptomics analysis
# Run once before the main analysis script
# ==============================================================================

cran_packages <- c(
  # Core data manipulation / visualization
  "tidyverse",
  "dplyr",
  "tidyr",
  "purrr",
  "tibble",
  "readr",
  "stringr",
  "ggplot2",

  # Seurat / spatial analysis
  "Seurat",
  "SeuratObject",
  "sctransform",

  # Parallelization and performance
  "future",
  "future.apply",
  "parallelly",
  "matrixStats",
  "data.table",

  # Statistics
  "moments",
  "broom",

  # Plotting helpers
  "patchwork",
  "cowplot",
  "ggrepel",
  "scales",
  "viridis",
  "pheatmap",

  # Miscellaneous commonly used utilities
  "RColorBrewer",
  "reshape2",
  "igraph"
)

installed <- rownames(installed.packages())

to_install <- setdiff(cran_packages, installed)

if (length(to_install) > 0) {
  install.packages(
    to_install,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
}

message("CRAN packages installed.")

purrr::walk(cran_packages, ~ library(.x, character.only = TRUE))

rm(cran_packages, installed, to_install)
