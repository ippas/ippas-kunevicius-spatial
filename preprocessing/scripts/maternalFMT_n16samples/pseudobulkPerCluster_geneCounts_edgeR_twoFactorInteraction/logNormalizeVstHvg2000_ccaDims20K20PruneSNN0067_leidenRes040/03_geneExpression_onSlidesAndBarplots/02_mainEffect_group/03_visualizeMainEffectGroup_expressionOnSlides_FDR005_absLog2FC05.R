#!/usr/bin/env Rscript
# IMPORTANT: NO sample x cluster minimum-spot filter is applied in these visualization summaries.
# Every available sample x cluster entry is retained.
# MN/MA/FN/FA percentages shown in the statistics header are read directly
# from the corresponding sex-specific edgeR result tables and are NOT
# recalculated from the Seurat object.
#
# ==============================================================================
# 03_visualizeMainEffectGroup_expressionOnSlides_FDR005_absLog2FC05_fixedPercentages.R
#
# Spatial-only visualization of genes significant in the cluster-specific
# edgeR main donor-group effect:
#
#   ASD vs Neurotypical averaged equally across Male and Female
#
# Main-effect gene-selection rule:
#   FDR < 0.05 and |log2FC| > 0.5
#
# IMPORTANT POST-HOC BH RULE:
#
# The sex-specific BH correction is calculated HERE in the runner, before
# plotting.
#
# For each cluster independently:
#   1. select ALL genes passing the main-effect rule;
#   2. take Male ASD vs Neurotypical nominal P-values for these genes;
#   3. apply BH across this selected set -> male_FDR_selected;
#   4. do the same independently for Female -> female_FDR_selected.
#
# The existing edgeR FDR from the complete sex-specific contrast is retained as:
#   male_FDR_global / female_FDR_global.
#
# `maximum_number_of_genes` limits plotting only. It does NOT limit the genes
# entering the new BH correction.
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
abs_log2fc_threshold <- 0.5
# Percent-positive summaries use every available sample × cluster value; no spot-count threshold is applied.
# Test run only. This does NOT affect the post-hoc BH universe.
maximum_number_of_genes <- NULL
# Optional explicit selection, for example:
# manual_target_genes <- c("Arhgef40", "Ndufb1")
manual_target_genes <- NULL
normalization_scale_factor <- 10000
upper_colour_quantile <- 0.99
show_histology_image <- FALSE
point_size_no_image <- 1.25
point_size_with_image <- 0.80
# Statistics-line colours.
# Main-effect statistics remain black.
# Male statistics + MN/MA percentages use blue.
# Female statistics + FN/FA percentages use red.
male_text_colour <- "#4C78A8"
female_text_colour <- "#C44E52"
# Header text can be tuned independently from the spatial panels.
plot_title_size <- 20
plot_subtitle_size <- 8.2
legend_height_ratio <- 0.075
# Slightly taller than before because the header now contains one statistics
# line per significant cluster.
spatial_pdf_width <- 18
spatial_pdf_height <- 24.5
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
  stop(
    "Metadata file does not exist: ",
    metadata_file,
call. = FALSE
  )
}
metadata_all <- read.delim(
  metadata_file,
sep = "\t",
header = TRUE,
stringsAsFactors = FALSE,
check.names = FALSE
)
if (!"sample_ID" %in% colnames(metadata_all)) {
  stop(
    "Metadata does not contain sample_ID.",
call. = FALSE
  )
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
# 3. Spatial colour settings
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
# ==============================================================================
# 4. Input and output paths
# ==============================================================================
functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "pseudobulkPerCluster_geneCounts_edgeR",
  "functions_geneExpressionOnSlides_pseudobulkPerClusterEdgeR.R"
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
  "FDR005_absLog2FC05",
  "03_expressionOnSlidesOnly"
)
sex_posthoc_statistics_file <- file.path(
  output_root_directory,
  "00_sexPostHocStatistics_BHwithinMainEffectSelectedGenes.tsv"
)
required_files <- c(
  functions_file,
  input_seurat_rdata_file,
  edger_results_rdata_file
)
missing_required_files <- required_files[
  !file.exists(required_files)
]
if (length(missing_required_files) > 0L) {
  stop(
    paste0(
      "Required file(s) do not exist:\n",
      paste(
        missing_required_files,
collapse = "\n"
      )
    ),
call. = FALSE
  )
}
source(functions_file)
dir.create(
  output_root_directory,
recursive = TRUE,
showWarnings = FALSE
)
# ==============================================================================
# 5. Prepare sex-specific post-hoc BH statistics
#
# IMPORTANT:
# This is calculated from ALL genes passing the main-effect filter.
# `maximum_number_of_genes` is intentionally not used here.
# ==============================================================================
message("Preparing sex-specific post-hoc BH statistics...")
edger_objects_for_posthoc <- load_edger_results_for_group_visualization(
  edger_results_rdata_file
)
main_effect_selected <- edger_objects_for_posthoc$combined_results[[
  "Overall_Group_ASD_vs_Neurotypical"
]] |>
  dplyr::filter(
    is.finite(.data$FDR),
    is.finite(.data$logFC),
    .data$FDR < fdr_threshold,
    abs(.data$logFC) > abs_log2fc_threshold
  ) |>
  dplyr::transmute(
cluster_id = as.character(.data$cluster_id),
ensembl_gene_id = as.character(.data$ensembl_gene_id),
gene = dplyr::if_else(
      is.na(.data$gene) | .data$gene == "",
      as.character(.data$ensembl_gene_id),
      as.character(.data$gene)
    ),
main_log2FC = .data$logFC,
main_PValue = .data$PValue,
main_FDR = .data$FDR
  )
male_percentage_columns <- c(
  "mean_percent_positive_spots_Neurotypical_Male",
  "mean_percent_positive_spots_ASD_Male"
)
female_percentage_columns <- c(
  "mean_percent_positive_spots_Neurotypical_Female",
  "mean_percent_positive_spots_ASD_Female"
)
missing_male_percentage_columns <- setdiff(
  male_percentage_columns,
  colnames(
    edger_objects_for_posthoc$combined_results[[
      "ASD_Male_vs_Neurotypical_Male"
    ]]
  )
)
missing_female_percentage_columns <- setdiff(
  female_percentage_columns,
  colnames(
    edger_objects_for_posthoc$combined_results[[
      "ASD_Female_vs_Neurotypical_Female"
    ]]
  )
)
if (
  length(missing_male_percentage_columns) > 0L ||
  length(missing_female_percentage_columns) > 0L
) {
  stop(
    "Required percent-positive columns are missing from the sex-specific edgeR result tables. Male missing: ",
    paste(missing_male_percentage_columns, collapse = ", "),
    "; Female missing: ",
    paste(missing_female_percentage_columns, collapse = ", "),
    call. = FALSE
  )
}
male_statistics <- edger_objects_for_posthoc$combined_results[[
  "ASD_Male_vs_Neurotypical_Male"
]] |>
  dplyr::transmute(
cluster_id = as.character(.data$cluster_id),
ensembl_gene_id = as.character(.data$ensembl_gene_id),
male_log2FC = .data$logFC,
male_PValue = .data$PValue,
male_FDR_global = .data$FDR,
mean_percent_positive_spots_Neurotypical_Male =
      .data$mean_percent_positive_spots_Neurotypical_Male,
mean_percent_positive_spots_ASD_Male =
      .data$mean_percent_positive_spots_ASD_Male
  )
female_statistics <- edger_objects_for_posthoc$combined_results[[
  "ASD_Female_vs_Neurotypical_Female"
]] |>
  dplyr::transmute(
cluster_id = as.character(.data$cluster_id),
ensembl_gene_id = as.character(.data$ensembl_gene_id),
female_log2FC = .data$logFC,
female_PValue = .data$PValue,
female_FDR_global = .data$FDR,
mean_percent_positive_spots_Neurotypical_Female =
      .data$mean_percent_positive_spots_Neurotypical_Female,
mean_percent_positive_spots_ASD_Female =
      .data$mean_percent_positive_spots_ASD_Female
  )
sex_posthoc_statistics_table <- main_effect_selected |>
  dplyr::left_join(
    male_statistics,
by = c(
      "cluster_id",
      "ensembl_gene_id"
    )
  ) |>
  dplyr::left_join(
    female_statistics,
by = c(
      "cluster_id",
      "ensembl_gene_id"
    )
  ) |>
  dplyr::group_by(.data$cluster_id) |>
  dplyr::mutate(
n_main_effect_selected_genes_in_cluster = dplyr::n(),
male_FDR_selected = stats::p.adjust(
      .data$male_PValue,
method = "BH"
    ),
female_FDR_selected = stats::p.adjust(
      .data$female_PValue,
method = "BH"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
main_effect_fdr_threshold = fdr_threshold,
main_effect_abs_log2FC_threshold = abs_log2fc_threshold,
FDR_global_definition = paste0(
      "Original edgeR FDR from the complete sex-specific contrast"
    ),
FDR_selected_definition = paste0(
      "BH within genes passing the main-effect filter in the same cluster"
    )
  ) |>
  dplyr::arrange(
    suppressWarnings(as.numeric(.data$cluster_id)),
    .data$cluster_id,
    .data$main_FDR,
    dplyr::desc(abs(.data$main_log2FC)),
    .data$gene
  )
write_edger_group_visualization_tsv(
  sex_posthoc_statistics_table,
  sex_posthoc_statistics_file
)
message(
  "Main-effect-selected gene x cluster results entering post-hoc BH: ",
  nrow(sex_posthoc_statistics_table)
)
message(
  "Clusters represented: ",
  dplyr::n_distinct(
    sex_posthoc_statistics_table$cluster_id
  )
)
message("Saved sex-specific post-hoc statistics:")
message(sex_posthoc_statistics_file)
rm(
  edger_objects_for_posthoc,
  main_effect_selected,
  male_statistics,
  female_statistics
)
invisible(gc(verbose = FALSE))
# ==============================================================================
# 6. Run spatial-only visualization workflow
#
# MN/MA/FN/FA values in the statistics header are taken directly from the
# sex-specific edgeR result tables. They are not recalculated from Seurat.
# Sample-panel fractions remain sample-specific descriptive labels.
# ==============================================================================
visualization_output <- run_main_effect_group_gene_expression_on_slides(
input_seurat_rdata_file = input_seurat_rdata_file,
edger_results_rdata_file = edger_results_rdata_file,
sex_posthoc_statistics_table = sex_posthoc_statistics_table,
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
normalization_scale_factor = normalization_scale_factor,
upper_colour_quantile = upper_colour_quantile,
red_palette_colors = red_palette_colors,
green_palette_colors = green_palette_colors,
palette_values = palette_values,
show_histology_image = show_histology_image,
point_size_no_image = point_size_no_image,
point_size_with_image = point_size_with_image,
male_text_colour = male_text_colour,
female_text_colour = female_text_colour,
plot_title_size = plot_title_size,
plot_subtitle_size = plot_subtitle_size,
legend_height_ratio = legend_height_ratio,
pdf_width = spatial_pdf_width,
pdf_height = spatial_pdf_height,
continue_after_gene_error = continue_after_gene_error
)
message("\nSex-specific post-hoc statistics file:")
message(sex_posthoc_statistics_file)
message("\nSelected genes file:")
message(visualization_output$selected_genes_file)
message("\nSelected gene-cluster statistics file:")
message(
  visualization_output$selected_gene_cluster_statistics_file
)
message("\nRun-status file:")
message(visualization_output$run_status_file)
message("\nPer-gene PDF directory:")
message(
  visualization_output$individual_gene_pdf_directory
)
message("\nCombined PDF:")
message(visualization_output$combined_pdf_file)
message("\nCombined PDF created:")
message(visualization_output$combined_pdf_created)
message("\nVisualization output directory:")
message(visualization_output$output_root_directory)
# ==============================================================================
# End
# ==============================================================================
