#!/usr/bin/env Rscript

# ==============================================================================
# 01_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction_fixed.R
#
# Runner script for the cluster-specific pseudobulk edgeR analysis.
#
# The edgeR model and statistical tests are unchanged:
#
#   expression ~ fmt_donor_group * sex
#
# Added descriptive outputs:
#   - mean and sample SD of percent-positive spots for each biological level;
#   - mean and sample SD of TMM-normalized CPM for each biological level;
#   - sample-level TMM-CPM matrices per cluster;
#   - sample-level percent-positive matrices per cluster;
#   - sample/cluster library-size and group metadata;
#   - XLSX formatting with gene before chromosome, frozen first row and first
#     two columns, and at most five displayed decimal places.
#
# IMPORTANT:
# The functions file keeps its original name:
#
#   functions_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction.R
#
# Replace the existing functions file with the updated file before running this
# runner.
# ==============================================================================


# ==============================================================================
# 1. Project and analysis configuration
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"

dataset_name <- "maternalFMT_n16samples"

clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"


# ==============================================================================
# 2. Input files
# ==============================================================================

input_rdata_file <- file.path(
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

metadata_file <- file.path(
  project_root,
  "data",
  "metadata_autismFMT.tsv"
)

annotation_rds_file <- file.path(
  project_root,
  "data",
  "ensembl115",
  "ensembl115_mouse_proteinCodingGenes.rds"
)

annotation_metadata_file <- file.path(
  project_root,
  "data",
  "ensembl115",
  "ensembl115_mouse_proteinCodingGenes_metadata.tsv"
)

functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "pseudobulkPerCluster_geneCounts_edgeR",
  "functions_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction.R"
)


# ==============================================================================
# 3. Output directory
# ==============================================================================

output_statistics_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)


# ==============================================================================
# 4. Validate inputs and load functions
# ==============================================================================

required_input_files <- c(
  input_rdata_file,
  metadata_file,
  annotation_rds_file,
  annotation_metadata_file,
  functions_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0L) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_input_files, collapse = "\n"),
    call. = FALSE
  )
}

source(functions_file)


# ==============================================================================
# 5. Explicit analysis settings
# ==============================================================================

seurat_object_name <- NULL

# The selected clustering column is set explicitly to avoid accidental
# automatic selection of another resolution.
cluster_column <- "leiden_res0.40"

# A spot is positive for a gene when the raw count is at least this value.
# This affects descriptive percent-positive metrics only.
positive_spot_min_count <- 1L

# Sample-cluster pseudobulks below this number of spots are excluded.
min_spots_per_sample_cluster <- 0L

# XLSX displays at most this many decimal places. Full numeric precision remains
# present in R objects and TSV files.
xlsx_max_decimal_places <- 5L


# ==============================================================================
# 6. Run analysis
# ==============================================================================

analysis_output <- run_pseudobulk_per_cluster_edger_two_factor_interaction(
  project_root = project_root,
  dataset_name = dataset_name,
  clustering_name = clustering_name,
  input_rdata_file = input_rdata_file,
  output_statistics_dir = output_statistics_dir,
  metadata_file = metadata_file,
  annotation_rds_file = annotation_rds_file,
  annotation_metadata_file = annotation_metadata_file,
  seurat_object_name = seurat_object_name,
  assay_name = "RNA",
  sample_column = "sample_ID",
  cluster_column = cluster_column,
  excluded_samples = c(
    "20_1F",
    "12_3F",
    "15_1M",
    "20_3M"
  ),
  expected_number_of_samples = 16L,
  expected_number_of_clusters = 16L,
  min_spots_per_sample_cluster = min_spots_per_sample_cluster,
  min_samples_per_factor_cell = 2L,
  positive_spot_min_count = positive_spot_min_count,
  xlsx_max_decimal_places = xlsx_max_decimal_places,
  filter_by_expr_parameters = list(
    min.count = 10,
    min.total.count = 15,
    large.n = 10,
    min.prop = 0.7
  ),
  stop_on_cluster_error = TRUE
)


# ==============================================================================
# 7. Main output paths
# ==============================================================================

message("\nMain combined workbook:")
message(analysis_output$combined_xlsx_file)

message("\nRData containing edgeR result tables and descriptive QC:")
message(analysis_output$result_rdata_file)

message("\nCombined abundance and detection QC table:")
message(analysis_output$combined_gene_abundance_qc_file)

message("\nPer-cluster sample-level CPM and detection metrics:")
message(analysis_output$per_cluster_sample_metrics_dir)

message("\nDefinitions of descriptive metrics:")
message(analysis_output$descriptive_metric_definitions_file)

# ==============================================================================
# End
# ==============================================================================
