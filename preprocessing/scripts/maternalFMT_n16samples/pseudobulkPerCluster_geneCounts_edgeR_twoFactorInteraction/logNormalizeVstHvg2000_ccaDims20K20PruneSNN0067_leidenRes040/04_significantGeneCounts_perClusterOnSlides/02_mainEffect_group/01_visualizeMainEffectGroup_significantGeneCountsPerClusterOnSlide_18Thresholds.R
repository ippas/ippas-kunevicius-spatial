#!/usr/bin/env Rscript

# ==============================================================================
# 01_visualizeMainEffectGroup_significantGeneCountsPerClusterOnSlide_18Thresholds_fixed4.R
#
# Version: fixed4_nestedPatchworkGrob_2026-08-03
#
# Simple runner for spatial visualization of the number of statistically
# significant genes per cluster for the edgeR main donor-group effect:
#
#   ASD vs Neurotypical averaged equally across Male and Female
#
# The selected Visium section is displayed once in each of three columns:
#   1. total significant genes per cluster (UP + DOWN), green scale;
#   2. UP genes per cluster, red scale;
#   3. DOWN genes per cluster, blue scale.
#
# Eighteen threshold variants are produced:
#   FDR < 0.10, 0.05 or 0.01
#   combined with
#   |log2FC| > 0.5, 0.6, 0.7, 0.8, 0.9 or 1.0
#
# All plotting and export logic is stored in:
#   preprocessing/src/pseudobulkPerCluster_geneCounts_edgeR/
#   functions_significantGeneCountsPerClusterOnSlide_pseudobulkPerClusterEdgeR_fixed11.R
# ============================================================================== 


# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "scales",
  "grid"
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
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

options(stringsAsFactors = FALSE)


# ==============================================================================
# 2. Main configuration
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"

dataset_name <- "maternalFMT_n16samples"

clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

cluster_column <- "leiden_res0.40"

# ------------------------------------------------------------------------------
# Select the one Visium section to display in all three panels.
# Change only this value to generate the plots for another section.
# ------------------------------------------------------------------------------

sample_id_to_plot <- "15_1F"

clustering_parameters_label <- paste0(
  "LogNormalize + VST | HVG 2000 | ",
  "CCA dims 20 | k.param 20 | ",
  "prune.SNN 0.0667 | Leiden resolution 0.40 | 16 clusters"
)

# Fixed anatomical colours used in the selected-clustering plots.
# In each legend entry:
#   first dot  = this anatomical cluster colour;
#   second dot = colour intensity representing the significant-gene count.
cluster_base_colours <- c(
  # Cortex
  "2"  = "#407ba7",
  "8"  = "#004e89",
  "7"  = "#002962",

  # Hippocampus
  "3"  = "#31cb00",
  "13" = "#119822",
  "16" = "#1e441e",

  # Thalamus
  "1"  = "#e66063",
  "9"  = "#d02224",
  "11" = "#9c191b",

  # Hypothalamus
  "6"  = "#b5838d",

  # Other anatomical regions
  "14" = "#6B1E2D",
  "4"  = "#adb5bd",
  "10" = "#ff9505",
  "12" = "#e85d04",
  "5"  = "#9e0059",
  "15" = "#ffd100"
)

cluster_names <- c(
  "1"  = "Thalamus 1",
  "2"  = "Cortex 1",
  "3"  = "Hippocampus 1",
  "4"  = "Fiber tracts",
  "5"  = "Striatum-like amygdalar nuclei",
  "6"  = "Hypothalamus",
  "7"  = "Cortex 3",
  "8"  = "Cortex 2",
  "9"  = "Thalamus 2",
  "10" = "Cortical subplate",
  "11" = "Thalamus 3",
  "12" = "Caudoputamen",
  "13" = "Hippocampus 2",
  "14" = "Vessels",
  "15" = "Ventricles",
  "16" = "Hippocampus 3"
)


# ==============================================================================
# 3. Eighteen statistical threshold variants
# ==============================================================================

fdr_thresholds <- c(
  0.10,
  0.05,
  0.01
)

abs_log2fc_thresholds <- c(
  0.5,
  0.6,
  0.7,
  0.8,
  0.9,
  1.0
)


# ==============================================================================
# 4. Spatial plotting settings
# ==============================================================================

image_scale <- "lowres"

show_histology_image <- FALSE

point_size_no_image <- 3.00
point_size_with_image <- 2.10

# Each panel has a two-dot legend:
# anatomical cluster colour + significant-gene-count colour.
legend_ncol <- 2L

# Vertical distance between consecutive legend rows.
# 1.00 = compact/original spacing; 1.20 = moderately increased spacing.
legend_row_spacing <- 1.40

# Distance between the bottom of the legend and the panel title, in points.
legend_bottom_margin <- 8

pdf_width <- 27
pdf_height <- 13.5

# TRUE: do not regenerate complete variants.
# FALSE: overwrite all 18 variants during this run.
skip_completed_variants <- FALSE

continue_after_variant_error <- TRUE


# ==============================================================================
# 5. Colour scales
# ==============================================================================

# Total significant results: green scale.
green_palette_colors <- c(
  "#D9D9D9",
  "#E5F5E0",
  "#A1D99B",
  "#74C476",
  "#31A354",
  "#006D2C",
  "#00441B"
)

# UP in ASD: red scale.
red_palette_colors <- c(
  "#D9D9D9",
  "#FEE5D9",
  "#FCAE91",
  "#FB6A4A",
  "#DE2D26",
  "#A50F15",
  "#67000D"
)

# DOWN in ASD: blue scale.
blue_palette_colors <- c(
  "#D9D9D9",
  "#EFF3FF",
  "#BDD7E7",
  "#6BAED6",
  "#3182BD",
  "#08519C",
  "#08306B"
)


# ==============================================================================
# 6. Input and output paths
# ==============================================================================

functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "pseudobulkPerCluster_geneCounts_edgeR",
  paste0(
    "functions_significantGeneCountsPerClusterOnSlide_",
    "pseudobulkPerClusterEdgeR.R"
  )
)

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
    "_logNormalizeVst_hvg2000_",
    "ccaDims20K20PruneSNN0067_",
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

output_root_directory <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "04_significantGeneCounts_perClusterOnSlides",
  "02_mainEffect_group"
)

if (!file.exists(functions_file)) {
  stop(
    "Functions file does not exist:\n",
    functions_file,
    call. = FALSE
  )
}

if (!file.exists(input_seurat_rdata_file)) {
  stop(
    "Input Seurat RData file does not exist:\n",
    input_seurat_rdata_file,
    call. = FALSE
  )
}

if (!file.exists(edger_results_rdata_file)) {
  stop(
    "edgeR results RData file does not exist:\n",
    edger_results_rdata_file,
    call. = FALSE
  )
}

dir.create(
  output_root_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

source(functions_file)


# ==============================================================================
# 7. Run the 18-variant workflow
# ==============================================================================

visualization_output <-
  run_main_group_significant_gene_count_spatial_workflow(
    input_seurat_rdata_file = input_seurat_rdata_file,
    edger_results_rdata_file = edger_results_rdata_file,
    output_root_directory = output_root_directory,
    sample_id_to_plot = sample_id_to_plot,
    cluster_column = cluster_column,
    clustering_parameters_label = clustering_parameters_label,
    sample_id_column = "sample_ID",
    seurat_object_name = NULL,
    image_scale = image_scale,
    fdr_thresholds = fdr_thresholds,
    abs_log2fc_thresholds = abs_log2fc_thresholds,
    green_palette_colors = green_palette_colors,
    red_palette_colors = red_palette_colors,
    blue_palette_colors = blue_palette_colors,
    cluster_base_colours = cluster_base_colours,
    cluster_names = cluster_names,
    show_histology_image = show_histology_image,
    point_size_no_image = point_size_no_image,
    point_size_with_image = point_size_with_image,
    legend_ncol = legend_ncol,
    legend_row_spacing = legend_row_spacing,
    legend_bottom_margin = legend_bottom_margin,
    pdf_width = pdf_width,
    pdf_height = pdf_height,
    skip_completed_variants = skip_completed_variants,
    continue_after_variant_error = continue_after_variant_error
  )

message("\nRun-status file:")
message(visualization_output$run_status_file)

message("\nThreshold-overview file:")
message(visualization_output$threshold_overview_file)

message("\nVisualization output directory:")
message(visualization_output$output_root_directory)

# ==============================================================================
# End
# ==============================================================================
