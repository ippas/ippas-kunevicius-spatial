#!/usr/bin/env Rscript

# ==============================================================================
# 01_visualizeAllEdgeREffects_significantGeneCountsPerClusterOnSlide_3x3.R
# Version: V16_NONOVERLAPPING_DOTS_N_LABEL_2026-08-14
#
# Runner for one 3 x 3 spatial overview of significant cluster-gene counts from:
#   row 1: Group x Sex interaction
#   row 2: main FMT donor-group effect
#   row 3: main sex effect
#
# Columns:
#   1. all significant results
#   2. positive-logFC results
#   3. negative-logFC results
#
# Two PDFs are generated automatically:
#   1. independent colour scale for every panel;
#   2. shared scales across effects separately for All / Positive / Negative.
#
#
# Legend layout:
#   - one shared legend per biological effect (3 legends total);
#   - 4 columns x 4 rows = 16 clusters;
#   - four dots per cluster: anatomical colour, All, UP/positive, DOWN/negative;
#   - entry text: cluster ID, n/UP/DOWN counts, then anatomical name.
#
# Default statistical selection:
#   FDR < 0.05
#   abs(log2FC) >= 0.5
#   no max-mean-percent filter
#
# To enable the optional percent-positive filter, change only:
#   min_max_mean_percent <- NULL
# to, for example:
#   min_max_mean_percent <- 25
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
# Visium section displayed in all nine panels.
# ------------------------------------------------------------------------------

sample_id_to_plot <- "15_1F"

# ------------------------------------------------------------------------------
# Statistical thresholds.
# ------------------------------------------------------------------------------

fdr_threshold <- 0.05
abs_log2fc_threshold <- 0.5

# Three Max mean % variants are generated automatically.
#
# Max mean % is calculated separately for every edgeR effect:
#
#   Interaction:
#     max across Male Neurotypical, Male ASD,
#                Female Neurotypical, Female ASD
#
#   FMT donor group:
#     max across Neurotypical and ASD
#
#   Sex:
#     max across Male and Female
#
# The two percentage intervals do not overlap:
#   0 <= Max mean % < 25
#   25 <= Max mean % <= 100
#
# A result with exactly 25% is therefore included only in the 25-100% variant.
percent_filter_variants <- list(
  all_percentages = list(
    min_max_mean_percent = NULL,
    max_max_mean_percent = NULL
  ),
  max_mean_percent_0_to_lt25 = list(
    min_max_mean_percent = 0,
    max_max_mean_percent = 25
  ),
  max_mean_percent_25_to_100 = list(
    min_max_mean_percent = 25,
    max_max_mean_percent = 100
  )
)

# Exact edgeR tests used in the 3 x 3 figure.
interaction_test_id <- "Interaction"
group_test_id <- "Overall_Group_ASD_vs_Neurotypical"
sex_test_id <- "Overall_Sex_Female_vs_Male"

clustering_parameters_label <- paste0(
  "LogNormalize + VST | HVG 2000 | ",
  "CCA dims 20 | k.param 20 | ",
  "prune.SNN 0.0667 | Leiden resolution 0.40 | 16 clusters"
)


# ==============================================================================
# 3. Cluster colours and names
# ==============================================================================

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

# Exact anatomical labels used in the current heatmap workflow.
# Keep these labels synchronized with the heatmap scripts.
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

expected_cluster_ids <- as.character(seq_len(16L))

if (!identical(names(cluster_names), expected_cluster_ids)) {
  stop(
    "`cluster_names` must contain exactly cluster IDs 1-16 in order.",
    call. = FALSE
  )
}


# ==============================================================================
# 4. Spatial plotting settings
# ==============================================================================

image_scale <- "lowres"
show_histology_image <- FALSE

point_size_no_image <- 3.20
point_size_with_image <- 2.10

# One shared legend is drawn above each biological-effect row.
# There are therefore only three legends in the complete 3 x 3 figure.
# Each legend contains all 16 clusters using the layout defined below.
# ------------------------------------------------------------------------------
# Shared legend layout
# ------------------------------------------------------------------------------

effect_legend_ncol <- 2L

# Cluster entries in the shared legends.
# Displayed as: Ck (n = ..., UP = ..., DOWN = ...) anatomical region
effect_legend_text_size <- 20

# The functions now space the four dots automatically according to this size,
# so large values no longer make the dots overlap.
effect_legend_point_size <- 12
effect_legend_row_spacing <- 1.60

# Shared legend title and its explanatory subtitle.
effect_legend_title_size <- 28
effect_legend_subtitle_size <- 22

effect_legend_bottom_margin <- 8
effect_legend_relative_height <- 0.43

# Fine horizontal correction AFTER left/center/right alignment.
# Keep 0 for true centering.
effect_legend_x_offset <- 0

# Physical blank gap between the complete first and second legend columns,
# in millimetres.
#
# Examples:
#   10 = compact
#   20 = current setting
#   25 = wider
#   30 = clearly wider
effect_legend_column_spacing <- 30

# Fixed horizontal coordinate range used by the legend.
# Keep this unchanged while tuning `effect_legend_column_spacing`.
# Increase only if the right-hand legend column becomes clipped.
effect_legend_canvas_width <- 100

# Vertical space, in points, between:
#   effect title (e.g. "FMT donor group")
#   and the explanatory subtitle beginning with "Dots:".
effect_legend_title_subtitle_spacing <- 6

# Vertical space, in points, between:
#   the explanatory "Dots:" subtitle
#   and the first row of cluster legend entries.
effect_legend_subtitle_entries_spacing <- 18

# Alignment of the COMPLETE legend block.
# Allowed values:
#   "left"   = move the complete legend block to the left
#   "center" = center the complete legend block
#   "right"  = move the complete legend block to the right
#
# This does NOT change alignment of individual legend labels:
# each cluster label remains left-aligned relative to its four dots.
effect_legend_block_alignment <- "center"

# ------------------------------------------------------------------------------
# Individual spatial-panel text
# ------------------------------------------------------------------------------

panel_title_size <- 22
panel_subtitle_size <- 22


# ------------------------------------------------------------------------------
# Main figure text
# ------------------------------------------------------------------------------

main_title_size <- 28
main_subtitle_size <- 22

# Distance below the main title before the main subtitle.
main_title_subtitle_spacing <- 5

# Distance below the complete main subtitle before the first effect legend.
main_subtitle_bottom_spacing <- 12

# Wider page and slightly taller canvas improve readability of the shared
# 3-column legends with larger one-line labels.
pdf_width <- 28
pdf_height <- 50


# ==============================================================================
# 5. Colour scales
# ==============================================================================

# All significant results.
green_palette_colors <- c(
  "#D9D9D9",
  "#E5F5E0",
  "#A1D99B",
  "#74C476",
  "#31A354",
  "#006D2C",
  "#00441B"
)

# Positive logFC.
red_palette_colors <- c(
  "#D9D9D9",
  "#FEE5D9",
  "#FCAE91",
  "#FB6A4A",
  "#DE2D26",
  "#A50F15",
  "#67000D"
)

# Negative logFC.
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

# IMPORTANT: keep the existing server-side function filename.
# The downloadable V5 file should be pasted into this existing file; the runner
# intentionally continues to source the old/stable server filename.
functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "pseudobulkPerCluster_geneCounts_edgeR",
  "functions_significantGeneCountsPerClusterOnSlide_allEdgeREffects.R"
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
  "04_significantGeneCounts_perClusterOnSlides"
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
# 7. Run all three Max mean % variants
# ==============================================================================

visualization_outputs <- list()
run_overview_rows <- list()

number_of_variants <- length(
  percent_filter_variants
)

for (variant_index in seq_along(percent_filter_variants)) {

  variant_name <- names(
    percent_filter_variants
  )[[variant_index]]

  variant <- percent_filter_variants[[variant_index]]

  message("\n")
  message("============================================================")
  message(
    "Running Max mean % variant ",
    variant_index,
    "/",
    number_of_variants,
    ": ",
    variant_name
  )
  message("============================================================")

  visualization_output <-
    run_all_effects_significant_gene_count_spatial_workflow(
      input_seurat_rdata_file = input_seurat_rdata_file,
      edger_results_rdata_file = edger_results_rdata_file,
      output_root_directory = output_root_directory,
      sample_id_to_plot = sample_id_to_plot,
      cluster_column = cluster_column,
      clustering_parameters_label = clustering_parameters_label,
      sample_id_column = "sample_ID",
      seurat_object_name = NULL,
      image_scale = image_scale,
      interaction_test_id = interaction_test_id,
      group_test_id = group_test_id,
      sex_test_id = sex_test_id,
      fdr_threshold = fdr_threshold,
      abs_log2fc_threshold = abs_log2fc_threshold,
      min_max_mean_percent = variant$min_max_mean_percent,
      max_max_mean_percent = variant$max_max_mean_percent,
      green_palette_colors = green_palette_colors,
      red_palette_colors = red_palette_colors,
      blue_palette_colors = blue_palette_colors,
      cluster_base_colours = cluster_base_colours,
      cluster_names = cluster_names,
      show_histology_image = show_histology_image,
      point_size_no_image = point_size_no_image,
      point_size_with_image = point_size_with_image,
      effect_legend_ncol = effect_legend_ncol,
      effect_legend_text_size = effect_legend_text_size,
      effect_legend_point_size = effect_legend_point_size,
      effect_legend_row_spacing = effect_legend_row_spacing,
      effect_legend_title_size = effect_legend_title_size,
      effect_legend_subtitle_size = effect_legend_subtitle_size,
      effect_legend_bottom_margin = effect_legend_bottom_margin,
      effect_legend_relative_height = effect_legend_relative_height,
      effect_legend_x_offset = effect_legend_x_offset,
      effect_legend_column_spacing = effect_legend_column_spacing,
      effect_legend_canvas_width = effect_legend_canvas_width,
      effect_legend_title_subtitle_spacing = effect_legend_title_subtitle_spacing,
      effect_legend_subtitle_entries_spacing = effect_legend_subtitle_entries_spacing,
      effect_legend_block_alignment = effect_legend_block_alignment,
      panel_title_size = panel_title_size,
      panel_subtitle_size = panel_subtitle_size,
      main_title_size = main_title_size,
      main_subtitle_size = main_subtitle_size,
      main_title_subtitle_spacing = main_title_subtitle_spacing,
      main_subtitle_bottom_spacing = main_subtitle_bottom_spacing,
      pdf_width = pdf_width,
      pdf_height = pdf_height
    )

  visualization_outputs[[variant_name]] <- visualization_output

  run_overview_rows[[variant_index]] <- tibble::tibble(
    variant_order = variant_index,
    variant_name = variant_name,

    min_max_mean_percent = if (
      is.null(variant$min_max_mean_percent)
    ) {
      NA_real_
    } else {
      as.numeric(variant$min_max_mean_percent)
    },

    max_max_mean_percent = if (
      is.null(variant$max_max_mean_percent)
    ) {
      NA_real_
    } else {
      as.numeric(variant$max_max_mean_percent)
    },

    output_directory =
      visualization_output$output_directory,

    independent_scale_pdf =
      visualization_output$output_files$independent_scale_pdf,

    shared_scale_pdf =
      visualization_output$output_files$shared_scale_pdf,

    counts_tsv =
      visualization_output$output_files$counts_tsv,

    selected_results_tsv =
      visualization_output$output_files$selected_results_tsv,

    scale_limits_tsv =
      visualization_output$output_files$scale_limits_tsv,

    summary_tsv =
      visualization_output$output_files$summary_tsv
  )
}


# ==============================================================================
# 8. Save and print run overview
# ==============================================================================

run_overview <- dplyr::bind_rows(
  run_overview_rows
)

run_overview_file <- file.path(
  output_root_directory,
  paste0(
    "00_",
    dataset_name,
    "_sample_",
    sample_id_to_plot,
    "_significantGeneCounts_percentFilterVariants_runOverview.tsv"
  )
)

utils::write.table(
  x = run_overview,
  file = run_overview_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

message("\n")
message("============================================================")
message("All Max mean % variants completed successfully.")
message("============================================================")

message("\nVariants:")

print(
  run_overview |>
    dplyr::select(
      "variant_order",
      "variant_name",
      "min_max_mean_percent",
      "max_max_mean_percent",
      "output_directory"
    ),
  n = Inf
)

message("\nRun-overview TSV:")
message(
  normalizePath(
    run_overview_file,
    mustWork = TRUE
  )
)

message("\nOutput root directory:")
message(
  normalizePath(
    output_root_directory,
    mustWork = TRUE
  )
)

# ==============================================================================
# End
# ==============================================================================
