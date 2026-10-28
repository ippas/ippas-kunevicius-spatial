#!/usr/bin/env Rscript

# ==============================================================================
# maternalFMT ASD - Spatial Transcriptomics Shiny app | v9.14
# Dataset/clustering-specific application
# ==============================================================================

required_packages <- c(
  "shiny",
  "Seurat",
  "SeuratObject",
  "Matrix",
  "dplyr",
  "tibble",
  "ggplot2",
  "scales",
  "ggnewscale",
  "patchwork",
  "DT",
  "edgeR",
  "ComplexHeatmap",
  "circlize"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(ggnewscale)
  library(patchwork)
  library(DT)
})

options(stringsAsFactors = FALSE)

# ==============================================================================
# Application directory
# ==============================================================================

# Run the application from the project root.
app_dir <- "preprocessing/src/shinyApp/maternalFMT_shiny_spatial_v9_14"

if (!dir.exists(app_dir)) {
  stop(
    "Application directory does not exist: ",
    app_dir,
    "\nRun app.R from the project root: /home/mateusz/projects/ippas-kunevicius-spatial",
    call. = FALSE
  )
}

# ==============================================================================
# Configuration
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"
clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

cluster_column <- "leiden_res0.40"
sample_id_column <- "sample_ID"
group_column <- "fmt_donor_group"
sex_column <- "sex"

expression_assay <- "RNA"
expression_layer_prefix <- "data."
counts_layer_prefix <- "counts."
image_scale <- "lowres"

input_seurat_rdata_file <- file.path(
  project_root,
  "results",
  dataset_name,
  "seurat_clustering_analysis",
  "logNormalize_vst",
  "hvg2000",
  "02_cca",
  "dims20",
  "k20",
  "prune0067",
  "RData",
  paste0(
    "02_",
    dataset_name,
    "_logNormalizeVst_hvg2000_ccaDims20K20PruneSNN0067_",
    "res010to100by010_multiClusteringAndUmap.RData"
  )
)

analysis_prefix <- paste0(
  dataset_name,
  "_",
  clustering_name,
  "_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction"
)

edger_results_rdata_file <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics",
  "03_edgeRResults",
  paste0(
    analysis_prefix,
    "_edgeRResults.RData"
  )
)

edger_model_objects_rdata_file <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics",
  "05_objects",
  paste0(
    analysis_prefix,
    "_edgeRModelObjects.RData"
  )
)


raw_pseudobulk_counts_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics",
  "01_pseudobulkInputs",
  "perCluster_rawCounts"
)

log2_qn_expression_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics",
  "01_pseudobulkInputs",
  "perCluster_rawCounts_log2plus1Transformation_quantileNormalize"
)

# Anatomical labels from the selected-clustering analysis.
cluster_names <- c(
  "1"  = "posterior & sensory relay thalamic nuclei",
  "2"  = "isocortex, layers 4 & 5",
  "3"  = "cortical layers 1 & hippocampal neuropil",
  "4"  = "Fiber tracts",
  "5"  = "cortical subplate & deep olfactory areas",
  "6"  = "hypothalamus",
  "7"  = "isocortex, layer 6",
  "8"  = "isocortex, layer 2/3",
  "9"  = "reticular, ventral geniculate & habenular region",
  "10" = "striatum-like amygdala nuclei",
  "11" = "medial thalamic nuclei",
  "12" = "caudoputamen",
  "13" = "hippocampal CA fields, pyramidal layer",
  "14" = "meninges & vasculature",
  "15" = "ventricles",
  "16" = "dentate gyrus"
)

cluster_colors <- c(
  "1"  = "#e66063",
  "2"  = "#407ba7",
  "3"  = "#31cb00",
  "4"  = "#adb5bd",
  "5"  = "#9e0059",
  "6"  = "#b5838d",
  "7"  = "#002962",
  "8"  = "#004e89",
  "9"  = "#d02224",
  "10" = "#ff9505",
  "11" = "#9c191b",
  "12" = "#e85d04",
  "13" = "#119822",
  "14" = "#6B1E2D",
  "15" = "#ffd100",
  "16" = "#1e441e"
)

red_expression_palette <- c(
  "#D9D9D9",
  "#FEE5D9",
  "#FCAE91",
  "#FB6A4A",
  "#DE2D26",
  "#A50F15",
  "#67000D"
)

green_expression_palette <- c(
  "#D9D9D9",
  "#E5F5E0",
  "#A1D99B",
  "#74C476",
  "#31A354",
  "#006D2C",
  "#00441B"
)

# ==============================================================================
# Local application functions
# ==============================================================================

source(
  file.path(app_dir, "src", "load_data.R"),
  local = TRUE
)
source(
  file.path(app_dir, "src", "spatial_visualization.R"),
  local = TRUE
)
source(
  file.path(app_dir, "src", "statistics.R"),
  local = TRUE
)
source(
  file.path(app_dir, "src", "statistics_browser.R"),
  local = TRUE
)
source(
  file.path(app_dir, "src", "heatmap_browser.R"),
  local = TRUE
)

# ==============================================================================
# Load data once when the application starts
# ==============================================================================

message("Loading Shiny application data...")

app_data <- load_shiny_app_data(
  input_seurat_rdata_file = input_seurat_rdata_file,
  edger_results_rdata_file = edger_results_rdata_file,
  edger_model_objects_rdata_file = edger_model_objects_rdata_file,
  raw_pseudobulk_counts_dir = raw_pseudobulk_counts_dir,
  log2_qn_expression_dir = log2_qn_expression_dir,
  expression_assay = expression_assay
)

message("Seurat object: ", app_data$seurat_object_name)
message("Genes: ", format(length(app_data$genes), big.mark = ","))
message("edgeR result tables: ", length(app_data$edger_results))
message("edgeR test definitions: ", nrow(app_data$edger_test_definitions))
message("edgeR model clusters: ", length(app_data$edger_model_objects))
message("Raw pseudobulk matrix clusters: ", length(app_data$per_cluster_expression$raw_counts))
message("log2+1 QN matrix clusters: ", length(app_data$per_cluster_expression$log2_qn))

# ==============================================================================
# UI and server
# ==============================================================================

source(
  file.path(app_dir, "ui.R"),
  local = TRUE
)
source(
  file.path(app_dir, "server.R"),
  local = TRUE
)

shinyApp(
  ui = ui,
  server = server
)
