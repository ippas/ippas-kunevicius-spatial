#!/usr/bin/env Rscript

# ==============================================================================
# 01_visualizeMainEffectGroup_FDR005_absLog2FC07.R
#
# Simple runner for visualization of genes significant in the cluster-specific
# edgeR main donor-group effect:
#
#   ASD vs Neurotypical averaged equally across Male and Female
#
# Gene-selection rule:
#   FDR < 0.05 and |log2FC| > 0.7
#
# Test run:
#   the first three genes ordered by minimum FDR and maximum |log2FC|.
#   Set `maximum_number_of_genes <- NULL` to process every selected gene.
#
# All plotting and export logic is stored in:
#   preprocessing/src/pseudobulkPerCluster_geneCounts_edgeR/
#   functions_geneExpressionOnSlidesAndBarplots_pseudobulkPerClusterEdgeR.R
# ============================================================================== 


# ==============================================================================
# 1. Configuration
# ==============================================================================

options(stringsAsFactors = FALSE)

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"

dataset_name <- "maternalFMT_n16samples"

clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

cluster_column <- "leiden_res0.40"
assay_name <- "RNA"

fdr_threshold <- 0.05
abs_log2fc_threshold <- 0.7

# Test run: process only three selected genes.
# Change to NULL after checking the output.
maximum_number_of_genes <- NULL

# Optional explicit selection, for example:
# manual_target_genes <- c("Nrp2", "Nppc", "Gpc3")
# NULL means: select automatically from the edgeR results.
manual_target_genes <- NULL

# edgeR used only sample × cluster pseudobulks containing at least 20 spots.
# The barplots apply the same threshold.
min_spots_per_sample_cluster <- 20L

normalization_scale_factor <- 10000
upper_colour_quantile <- 0.99

show_histology_image <- FALSE
skip_completed_genes <- TRUE
continue_after_gene_error <- TRUE


# ==============================================================================
# 2. Samples retained after QC
# ==============================================================================

excluded_samples <- c(
  "20_1F",
  "12_3F",
  "15_1M",
  "20_3M"
)

metadata_file <- file.path(
  project_root,
  "data",
  "metadata_autismFMT.tsv"
)

if (!file.exists(metadata_file)) {
  stop("Metadata file does not exist: ", metadata_file, call. = FALSE)
}

metadata_all <- read.delim(
  metadata_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"sample_ID" %in% colnames(metadata_all)) {
  stop("Metadata does not contain sample_ID.", call. = FALSE)
}

included_sample_ids <- as.character(
  metadata_all$sample_ID[
    !metadata_all$sample_ID %in% excluded_samples
  ]
)

if (length(included_sample_ids) != 16L) {
  stop(
    "Expected 16 retained samples, but found ",
    length(included_sample_ids),
    ".",
    call. = FALSE
  )
}


# ==============================================================================
# 3. Colour settings
# ==============================================================================

red_palette_colors <- c(
  "#D9D9D9",
  "#FEE5D9",
  "#FCAE91",
  "#FB6A4A",
  "#DE2D26",
  "#A50F15",
  "#67000D"
)

green_palette_colors <- c(
  "#D9D9D9",
  "#E5F5E0",
  "#A1D99B",
  "#74C476",
  "#31A354",
  "#006D2C",
  "#00441B"
)

palette_values <- c(
  0.000,
  0.006,
  0.060,
  0.180,
  0.400,
  0.700,
  1.000
)

# Four barplot groups in the required order:
# Male Neurotypical, Male ASD, Female Neurotypical, Female ASD.
group_colors <- c(
  "male_neurotypical" = "#4C78A8",
  "male_asd" = "#C44E52",
  "female_neurotypical" = "#8AB8DB",
  "female_asd" = "#E78A8A"
)


# ==============================================================================
# 4. Plot sizes
# ==============================================================================

point_size_no_image <- 1.1
point_size_with_image <- 0.80

spatial_pdf_width <- 18
spatial_pdf_height <- 23

barplot_pdf_width <- 11
barplot_pdf_height <- 8


# ==============================================================================
# 5. Input and output paths
# ==============================================================================

functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "pseudobulkPerCluster_geneCounts_edgeR",
  "functions_geneExpressionOnSlidesAndBarplots_pseudobulkPerClusterEdgeR.R"
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
  "03_geneExpression_onSlidesAndBarplots",
  "02_mainEffect_group",
  "FDR005_absLog2FC07"
)

if (!file.exists(functions_file)) {
  stop("Functions file does not exist: ", functions_file, call. = FALSE)
}

if (!file.exists(input_seurat_rdata_file)) {
  stop(
    "Input Seurat RData file does not exist: ",
    input_seurat_rdata_file,
    call. = FALSE
  )
}

if (!file.exists(edger_results_rdata_file)) {
  stop(
    "edgeR results RData file does not exist: ",
    edger_results_rdata_file,
    call. = FALSE
  )
}

source(functions_file)


# ==============================================================================
# 6. Run visualization workflow
# ==============================================================================

visualization_output <- run_main_effect_group_gene_visualizations(
  input_seurat_rdata_file = input_seurat_rdata_file,
  edger_results_rdata_file = edger_results_rdata_file,
  output_root_directory = output_root_directory,
  cluster_column = cluster_column,
  assay_name = assay_name,
  sample_id_column = "sample_ID",
  group_column = "fmt_donor_group",
  sex_column = "sex",
  included_sample_ids = included_sample_ids,
  seurat_object_name = NULL,
  fdr_threshold = fdr_threshold,
  abs_log2fc_threshold = abs_log2fc_threshold,
  maximum_number_of_genes = maximum_number_of_genes,
  manual_target_genes = manual_target_genes,
  min_spots_per_sample_cluster = min_spots_per_sample_cluster,
  normalization_scale_factor = normalization_scale_factor,
  upper_colour_quantile = upper_colour_quantile,
  red_palette_colors = red_palette_colors,
  green_palette_colors = green_palette_colors,
  palette_values = palette_values,
  group_colors = group_colors,
  show_histology_image = show_histology_image,
  point_size_no_image = point_size_no_image,
  point_size_with_image = point_size_with_image,
  spatial_pdf_width = spatial_pdf_width,
  spatial_pdf_height = spatial_pdf_height,
  barplot_pdf_width = barplot_pdf_width,
  barplot_pdf_height = barplot_pdf_height,
  skip_completed_genes = skip_completed_genes,
  continue_after_gene_error = continue_after_gene_error
)

message("\nSelected genes file:")
message(visualization_output$selected_genes_file)

message("\nRun-status file:")
message(visualization_output$run_status_file)

message("\nVisualization output directory:")
message(visualization_output$output_root_directory)

# ==============================================================================
# End
# ==============================================================================
