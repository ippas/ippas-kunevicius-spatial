#!/usr/bin/env Rscript

# ==============================================================================
# 01_heatmap_interaction_FDR005_perCluster_oneHeatmap_FDRlegendOnly_log2QN.R
#
# One heatmap for significant Group x Sex interaction results from the
# cluster-specific pseudobulk edgeR analysis.
#
# ROWS
#   Significant cluster-gene interaction results at FDR < 0.05.
#   The same gene can occur more than once when significant in several clusters.
#
# COLUMNS
#   Samples ordered as:
#     Male Neurotypical -> Male ASD -> Female Neurotypical -> Female ASD
#
# HEATMAP VALUES
#   The heatmap DOES NOT use TMM-normalized CPM.
#
#   For each cluster it reads the already prepared matrix:
#
#     raw pseudobulk counts
#       -> log2(raw count + 1)
#       -> quantile normalization across the 16 samples WITHIN that cluster
#
#   This heatmap script then performs only:
#
#     precomputed log2 + quantile-normalized values
#       -> row-wise Z-score
#       -> clipping to [-2, 2]
#
#   No second log2 transformation and no second quantile normalization are
#   performed here.
#
# INPUT EXPRESSION DIRECTORY
#   01_statistics/
#   01_pseudobulkInputs/
#   perCluster_rawCounts_log2plus1Transformation_quantileNormalize/
#
# OUTPUT DIRECTORY
#   The heatmap remains in the same location:
#
#   02_heatmaps/
#   01_interaction/
#   FDR005/
#
# IMPORTANT
#   edgeR significance is still read from *_edgeRResults.RData.
#   The log2 + QN matrices are used only for heatmap expression values.
# ==============================================================================

options(stringsAsFactors = FALSE)

# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tibble",
  "readr",
  "ComplexHeatmap",
  "circlize",
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
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ==============================================================================
# 2. Analysis configuration
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"

dataset_name <- "maternalFMT_n16samples"

clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

fdr_threshold <- 0.05

zscore_lower_limit <- -2
zscore_upper_limit <- 2

expected_number_of_samples <- 16L

# Row-clustering parameters.
# Base distance: 2 - Pearson correlation.
# Distances for rows from the same source cluster are multiplied by 0.80,
# whereas distances for rows from different clusters are multiplied by 1.20.
row_clustering_base_distance <- "2 - r"
row_clustering_r_scale <- 1
row_clustering_same_cluster_bonus <- 0.20
row_clustering_linkage_method <- "complete"

clustering_title <- paste(
  "LogNormalize + VST | HVG = 2000 | CCA dims = 20 |",
  "k = 20 | prune.SNN = 0.0667 | Leiden resolution = 0.40"
)

# ==============================================================================

# 3. Cluster colours and anatomical labels

# ==============================================================================

custom_cluster_colors <- c(

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

custom_cluster_labels <- c(

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

# ==============================================================================
# 4. Input and output paths
# ==============================================================================

statistics_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)

analysis_prefix <- paste0(
  dataset_name,
  "_",
  clustering_name,
  "_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction"
)

edgeR_results_rdata_file <- file.path(
  statistics_dir,
  "03_edgeRResults",
  paste0(
    analysis_prefix,
    "_edgeRResults.RData"
  )
)

metadata_file <- file.path(
  project_root,
  "data",
  "metadata_autismFMT.tsv"
)

log2_qn_input_dir <- file.path(
  statistics_dir,
  "01_pseudobulkInputs",
  "perCluster_rawCounts_log2plus1Transformation_quantileNormalize"
)

heatmap_output_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "02_heatmaps",
  "01_interaction",
  "FDR005"
)

dir.create(
  heatmap_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

required_input_files <- c(
  edgeR_results_rdata_file,
  metadata_file
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

if (!dir.exists(log2_qn_input_dir)) {
  stop(
    "Missing log2 + quantile-normalized input directory:\n",
    log2_qn_input_dir,
    call. = FALSE
  )
}

# ==============================================================================
# 5. Output files
# ==============================================================================

heatmap_pdf <- file.path(
  heatmap_output_dir,
  paste0(
    "01_",
    dataset_name,
    "_",
    clustering_name,
    "_interaction_FDRlt005_rawCountsLog2plus1QN_rowZscore_heatmap.pdf"
  )
)

selected_results_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "02_",
    dataset_name,
    "_interaction_FDRlt005_selectedClusterGeneResults.tsv"
  )
)

row_annotations_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "03_",
    dataset_name,
    "_interaction_FDRlt005_rowAnnotations.tsv"
  )
)

log2_qn_selected_matrix_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "04_",
    dataset_name,
    "_interaction_FDRlt005_log2plus1QN_selectedMatrix.tsv"
  )
)

zscore_matrix_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "05_",
    dataset_name,
    "_interaction_FDRlt005_log2plus1QN_rowZscoreMatrix.tsv"
  )
)

sample_order_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "06_",
    dataset_name,
    "_interaction_FDRlt005_sampleOrder.tsv"
  )
)

summary_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "07_",
    dataset_name,
    "_interaction_FDRlt005_heatmapSummary.tsv"
  )
)

# ==============================================================================
# 6. Helper functions
# ==============================================================================

find_interaction_test_id <- function(test_definitions) {

  required_columns <- c(

    "test_id",

    "sheet_name",

    "comparison"

  )

  missing_columns <- setdiff(

    required_columns,

    colnames(test_definitions)

  )

  if (length(missing_columns) > 0L) {

    stop(

      "Test-definition table is missing: ",

      paste(missing_columns, collapse = ", "),

call. = FALSE

    )

  }

  candidate_rows <- test_definitions |>

    dplyr::filter(

      grepl("interaction", .data$test_id, ignore.case = TRUE) |

        grepl("interaction", .data$sheet_name, ignore.case = TRUE) |

        grepl("interaction", .data$comparison, ignore.case = TRUE)

    )

  if (nrow(candidate_rows) != 1L) {

    stop(

      "Could not uniquely identify the interaction test. Candidates found: ",

      nrow(candidate_rows),

call. = FALSE

    )

  }

  as.character(candidate_rows$test_id[[1]])

}

calculate_row_zscores <- function(

expression_matrix,

lower_limit = -2,

upper_limit = 2

) {

  output_matrix <- matrix(

    NA_real_,

nrow = nrow(expression_matrix),

ncol = ncol(expression_matrix),

dimnames = dimnames(expression_matrix)

  )

  for (row_index in seq_len(nrow(expression_matrix))) {

    values <- as.numeric(expression_matrix[row_index, ])

    finite_values <- is.finite(values)

    if (sum(finite_values) < 2L) {

      output_matrix[row_index, finite_values] <- 0

      next

    }

    row_mean <- mean(values[finite_values])

    row_sd <- stats::sd(values[finite_values])

    if (!is.finite(row_sd) || row_sd == 0) {

      output_matrix[row_index, finite_values] <- 0

      next

    }

    output_matrix[row_index, finite_values] <-

      (values[finite_values] - row_mean) / row_sd

  }

  output_matrix[

    is.finite(output_matrix) & output_matrix < lower_limit

  ] <- lower_limit

  output_matrix[

    is.finite(output_matrix) & output_matrix > upper_limit

  ] <- upper_limit

  output_matrix

}

extract_cluster_id_from_log2_qn_filename <- function(file_name) {
  cluster_id <- sub(
    "^.*_cluster_([^_]+)_rawCounts_log2plus1Transformation_quantileNormalize_proteinCodingAutosomalAndX\\.tsv\\.gz$",
    "\\1",
    basename(file_name)
  )

  if (identical(cluster_id, basename(file_name))) {
    stop(
      "Could not extract cluster ID from log2 + QN file: ",
      basename(file_name),
      call. = FALSE
    )
  }

  cluster_id
}

read_log2_qn_cluster_matrix <- function(
    input_file,
    expected_sample_ids,
    cluster_id
) {
  input_df <- readr::read_tsv(
    input_file,
    show_col_types = FALSE,
    progress = FALSE
  )

  required_columns <- c(
    "ensembl_gene_id",
    "gene",
    "chromosome"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(input_df)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Cluster ",
      cluster_id,
      " log2 + QN file is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyDuplicated(input_df$ensembl_gene_id)) {
    stop(
      "Duplicated Ensembl gene IDs in log2 + QN file for cluster ",
      cluster_id,
      ".",
      call. = FALSE
    )
  }

  sample_columns <- setdiff(
    colnames(input_df),
    required_columns
  )

  missing_samples <- setdiff(
    expected_sample_ids,
    sample_columns
  )

  extra_samples <- setdiff(
    sample_columns,
    expected_sample_ids
  )

  if (length(missing_samples) > 0L || length(extra_samples) > 0L) {
    stop(
      "Unexpected sample columns in log2 + QN file for cluster ",
      cluster_id,
      ". Missing: ",
      paste(missing_samples, collapse = ", "),
      "; extra: ",
      paste(extra_samples, collapse = ", "),
      call. = FALSE
    )
  }

  input_df <- input_df |>
    dplyr::select(
      dplyr::all_of(required_columns),
      dplyr::all_of(expected_sample_ids)
    )

  expression_matrix <- as.matrix(
    input_df[
      ,
      expected_sample_ids,
      drop = FALSE
    ]
  )

  storage.mode(expression_matrix) <- "double"

  rownames(expression_matrix) <- as.character(
    input_df$ensembl_gene_id
  )

  if (anyNA(expression_matrix)) {
    stop(
      "NA values found in log2 + QN matrix for cluster ",
      cluster_id,
      ".",
      call. = FALSE
    )
  }

  if (any(!is.finite(expression_matrix))) {
    stop(
      "Non-finite values found in log2 + QN matrix for cluster ",
      cluster_id,
      ".",
      call. = FALSE
    )
  }

  list(
    matrix = expression_matrix,
    annotation = input_df |>
      dplyr::select(
        .data$ensembl_gene_id,
        .data$gene,
        .data$chromosome
      )
  )
}

modified_correlation_distance <- function(

expression_matrix,

cluster_labels,

base_distance = "2 - r",

r_scale = 1,

same_cluster_bonus = 0.25

) {

  if (!is.matrix(expression_matrix)) {

    expression_matrix <- as.matrix(expression_matrix)

  }

  if (nrow(expression_matrix) <= 1L) {

    return(stats::as.dist(matrix(0, nrow = 1, ncol = 1)))

  }

  cluster_labels <- as.character(cluster_labels)

  if (length(cluster_labels) != nrow(expression_matrix)) {

    stop(

      "The number of cluster labels must equal the number of heatmap rows.",

call. = FALSE

    )

  }

  if (!is.numeric(r_scale) || length(r_scale) != 1L || !is.finite(r_scale)) {

    stop("`r_scale` must be one finite numeric value.", call. = FALSE)

  }

  if (

    !is.numeric(same_cluster_bonus) ||

      length(same_cluster_bonus) != 1L ||

      !is.finite(same_cluster_bonus) ||

      same_cluster_bonus < 0 ||

      same_cluster_bonus >= 1

  ) {

    stop(

      "`same_cluster_bonus` must be a finite number in [0, 1).",

call. = FALSE

    )

  }

  cor_mat <- suppressWarnings(

    stats::cor(

      t(expression_matrix),

method = "pearson",

use = "pairwise.complete.obs"

    )

  )

  cor_mat[!is.finite(cor_mat)] <- 0

  diag(cor_mat) <- 1

  base_dist <- if (identical(base_distance, "2 - r")) {

    2 - (cor_mat * r_scale)

  } else if (identical(base_distance, "1 - r")) {

    1 - (cor_mat * r_scale)

  } else {

    stop(

      "`base_distance` must be either '2 - r' or '1 - r'.",

call. = FALSE

    )

  }

  same_cluster_matrix <- outer(

    cluster_labels,

    cluster_labels,

FUN = "=="

  )

  dist_mat <- base_dist

  dist_mat[same_cluster_matrix] <-

    base_dist[same_cluster_matrix] * (1 - same_cluster_bonus)

  dist_mat[!same_cluster_matrix] <-

    base_dist[!same_cluster_matrix] * (1 + same_cluster_bonus)

  diag(dist_mat) <- 0

  dist_mat[dist_mat < 0] <- 0

  stats::as.dist(dist_mat)

}

create_modified_row_hclust <- function(

expression_matrix,

annotation_table,

base_distance = "2 - r",

r_scale = 1,

same_cluster_bonus = 0.25,

linkage_method = "complete"

) {

  matched_annotation <- annotation_table[

    match(rownames(expression_matrix), annotation_table$row_id),

    ,

drop = FALSE

  ]

  if (

    nrow(matched_annotation) != nrow(expression_matrix) ||

      anyNA(matched_annotation$row_id) ||

      anyNA(matched_annotation$cluster_id)

  ) {

    stop(

      "Failed to match source-cluster annotations to all heatmap rows.",

call. = FALSE

    )

  }

  modified_distance <- modified_correlation_distance(

expression_matrix = expression_matrix,

cluster_labels = matched_annotation$cluster_id,

base_distance = base_distance,

r_scale = r_scale,

same_cluster_bonus = same_cluster_bonus

  )

  stats::hclust(

    modified_distance,

method = linkage_method

  )

}

format_cluster_label <- function(cluster_id) {

  paste0(

    "C",

    cluster_id,

    " | ",

    unname(custom_cluster_labels[cluster_id])

  )

}

make_fdr_marker_legend <- function(significant_results) {

  required_columns <- c(

    "FDR"

  )

  missing_columns <- setdiff(

    required_columns,

    colnames(significant_results)

  )

  if (length(missing_columns) > 0L) {

    stop(

      "FDR legend input is missing required column(s): ",

      paste(missing_columns, collapse = ", "),

call. = FALSE

    )

  }

  n_fdr_005 <- sum(

    !is.na(significant_results$FDR) &

      significant_results$FDR < 0.05

  )

  n_fdr_001 <- sum(

    !is.na(significant_results$FDR) &

      significant_results$FDR < 0.01

  )

  n_fdr_0001 <- sum(

    !is.na(significant_results$FDR) &

      significant_results$FDR < 0.001

  )

  ComplexHeatmap::Legend(

title = "FDR marker",

labels = c(

      paste0(

        "FDR < 0.05 (n = ",

        n_fdr_005,

        ")"

      ),

      paste0(

        "FDR < 0.01 (n = ",

        n_fdr_001,

        ")"

      ),

      paste0(

        "FDR < 0.001 (n = ",

        n_fdr_0001,

        ")"

      )

    ),

graphics = list(

      function(x, y, w, h) {

        grid::grid.text(

          "*",

x = x,

y = y,

gp = grid::gpar(

fontsize = 10,

fontface = "bold"

          )

        )

      },

      function(x, y, w, h) {

        grid::grid.text(

          "**",

x = x,

y = y,

gp = grid::gpar(

fontsize = 10,

fontface = "bold"

          )

        )

      },

      function(x, y, w, h) {

        grid::grid.text(

          "***",

x = x,

y = y,

gp = grid::gpar(

fontsize = 10,

fontface = "bold"

          )

        )

      }

    ),

grid_width = grid::unit(

      7,

      "mm"

    ),

grid_height = grid::unit(

      4,

      "mm"

    ),

title_gp = grid::gpar(

fontsize = 9,

fontface = "bold"

    ),

labels_gp = grid::gpar(

fontsize = 8

    )

  )

}

# ==============================================================================
# 7. Load edgeR result tables
# ==============================================================================

message("Loading edgeR results RData:")
message(edgeR_results_rdata_file)

load(edgeR_results_rdata_file)

required_objects <- c(
  "edgeR_perCluster_combinedResults",
  "edgeR_perCluster_testDefinitions"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1)
  )
]

if (length(missing_objects) > 0L) {
  stop(
    "Required object(s) missing after loading edgeR results RData: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

# ==============================================================================
# 8. Select significant interaction results
# ==============================================================================

interaction_test_id <- find_interaction_test_id(
  edgeR_perCluster_testDefinitions
)

message("\nSelected interaction test_id: ", interaction_test_id)

if (!interaction_test_id %in% names(edgeR_perCluster_combinedResults)) {
  stop(
    "Interaction test is absent from edgeR_perCluster_combinedResults: ",
    interaction_test_id,
    call. = FALSE
  )
}

interaction_results <-
  edgeR_perCluster_combinedResults[[interaction_test_id]] |>
  tibble::as_tibble()

required_interaction_columns <- c(
  "cluster_id",
  "ensembl_gene_id",
  "gene",
  "logFC",
  "PValue",
  "FDR",
  "mean_percent_positive_spots_Male_Neurotypical",
  "mean_percent_positive_spots_Male_ASD",
  "mean_percent_positive_spots_Female_Neurotypical",
  "mean_percent_positive_spots_Female_ASD"
)

missing_interaction_columns <- setdiff(
  required_interaction_columns,
  colnames(interaction_results)
)

if (length(missing_interaction_columns) > 0L) {
  stop(
    "Interaction results are missing required columns: ",
    paste(missing_interaction_columns, collapse = ", "),
    call. = FALSE
  )
}

significant_interaction <- interaction_results |>
  dplyr::filter(
    !is.na(.data$FDR),
    .data$FDR < fdr_threshold
  ) |>
  dplyr::mutate(
    cluster_id = as.character(.data$cluster_id),
    interaction_direction = dplyr::case_when(
      .data$logFC > 0 ~ "Positive logFC",
      .data$logFC < 0 ~ "Negative logFC",
      TRUE ~ "Zero logFC"
    ),
    max_mean_percent = pmax(
      .data$mean_percent_positive_spots_Male_Neurotypical,
      .data$mean_percent_positive_spots_Male_ASD,
      .data$mean_percent_positive_spots_Female_Neurotypical,
      .data$mean_percent_positive_spots_Female_ASD,
      na.rm = TRUE
    )
  ) |>
  dplyr::arrange(
    as.integer(.data$cluster_id),
    factor(
      .data$interaction_direction,
      levels = c(
        "Positive logFC",
        "Negative logFC",
        "Zero logFC"
      )
    ),
    .data$FDR,
    dplyr::desc(abs(.data$logFC))
  )

if (nrow(significant_interaction) == 0L) {
  stop(
    "No interaction result satisfies FDR < ",
    fdr_threshold,
    call. = FALSE
  )
}

significant_clusters <- sort(
  unique(
    as.integer(
      significant_interaction$cluster_id
    )
  )
)

significant_cluster_ids <- as.character(
  significant_clusters
)

message(
  "Significant cluster-gene results: ",
  nrow(significant_interaction)
)

message(
  "Unique genes: ",
  length(unique(significant_interaction$gene))
)

message(
  "Significant clusters: ",
  paste(significant_cluster_ids, collapse = ", ")
)

# ==============================================================================
# 9. Locate log2 + QN matrices and prepare sample metadata
# ==============================================================================

log2_qn_files <- list.files(
  log2_qn_input_dir,
  pattern = "_rawCounts_log2plus1Transformation_quantileNormalize_proteinCodingAutosomalAndX\\.tsv\\.gz$",
  full.names = TRUE
)

if (length(log2_qn_files) == 0L) {
  stop(
    "No log2 + QN per-cluster matrices found in:\n",
    log2_qn_input_dir,
    call. = FALSE
  )
}

log2_qn_cluster_ids <- vapply(
  log2_qn_files,
  extract_cluster_id_from_log2_qn_filename,
  character(1)
)

if (anyDuplicated(log2_qn_cluster_ids)) {
  stop(
    "More than one log2 + QN matrix was found for at least one cluster.",
    call. = FALSE
  )
}

names(log2_qn_files) <- log2_qn_cluster_ids

missing_qn_clusters <- setdiff(
  significant_cluster_ids,
  names(log2_qn_files)
)

if (length(missing_qn_clusters) > 0L) {
  stop(
    "Missing log2 + QN matrix for significant cluster(s): ",
    paste(missing_qn_clusters, collapse = ", "),
    call. = FALSE
  )
}

first_qn_file <- log2_qn_files[
  significant_cluster_ids[[1]]
]

first_qn_header <- readr::read_tsv(
  first_qn_file,
  n_max = 0,
  show_col_types = FALSE,
  progress = FALSE
)

expression_metadata_columns <- c(
  "ensembl_gene_id",
  "gene",
  "chromosome"
)

sample_ids_in_expression_files <- setdiff(
  colnames(first_qn_header),
  expression_metadata_columns
)

if (length(sample_ids_in_expression_files) != expected_number_of_samples) {
  stop(
    "Expected ",
    expected_number_of_samples,
    " sample columns in log2 + QN matrices but found ",
    length(sample_ids_in_expression_files),
    ".",
    call. = FALSE
  )
}

metadata_all <- readr::read_tsv(
  metadata_file,
  show_col_types = FALSE,
  progress = FALSE
)

required_metadata_columns <- c(
  "sample_ID",
  "fmt_donor_group",
  "sex"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata_all)
)

if (length(missing_metadata_columns) > 0L) {
  stop(
    "Metadata file is missing required column(s): ",
    paste(missing_metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

sample_metadata <- metadata_all |>
  dplyr::transmute(
    sample_ID = trimws(as.character(.data$sample_ID)),
    fmt_donor_group = trimws(as.character(.data$fmt_donor_group)),
    sex = trimws(as.character(.data$sex))
  ) |>
  dplyr::filter(
    .data$sample_ID %in% sample_ids_in_expression_files
  ) |>
  dplyr::distinct()

if (anyDuplicated(sample_metadata$sample_ID)) {
  stop(
    "Duplicated sample_ID values remain in sample metadata.",
    call. = FALSE
  )
}

missing_metadata_samples <- setdiff(
  sample_ids_in_expression_files,
  sample_metadata$sample_ID
)

if (length(missing_metadata_samples) > 0L) {
  stop(
    "Missing metadata for sample(s): ",
    paste(missing_metadata_samples, collapse = ", "),
    call. = FALSE
  )
}

sample_metadata <- sample_metadata |>
  dplyr::mutate(
    sex = factor(
      .data$sex,
      levels = c("Male", "Female")
    ),
    fmt_donor_group = factor(
      .data$fmt_donor_group,
      levels = c("Neurotypical", "ASD")
    )
  )

if (anyNA(sample_metadata$sex)) {
  stop(
    "Unexpected sex value in sample metadata.",
    call. = FALSE
  )
}

if (anyNA(sample_metadata$fmt_donor_group)) {
  stop(
    "Unexpected fmt_donor_group value in sample metadata.",
    call. = FALSE
  )
}

sample_metadata <- sample_metadata |>
  dplyr::arrange(
    .data$sex,
    .data$fmt_donor_group,
    .data$sample_ID
  ) |>
  dplyr::mutate(
    sample_order = dplyr::row_number(),
    combined_group = factor(
      paste(
        as.character(.data$sex),
        as.character(.data$fmt_donor_group),
        sep = "_"
      ),
      levels = c(
        "Male_Neurotypical",
        "Male_ASD",
        "Female_Neurotypical",
        "Female_ASD"
      )
    )
  )

if (nrow(sample_metadata) != expected_number_of_samples) {
  stop(
    "Expected ",
    expected_number_of_samples,
    " samples after metadata matching but found ",
    nrow(sample_metadata),
    ".",
    call. = FALSE
  )
}

sample_order <- as.character(
  sample_metadata$sample_ID
)

message("\nSample order used in the heatmap:")

print(
  sample_metadata |>
    dplyr::select(
      .data$sample_order,
      .data$sample_ID,
      .data$sex,
      .data$fmt_donor_group,
      .data$combined_group
    ),
  n = Inf
)

# ==============================================================================
# 10. Read precomputed log2(count + 1) + quantile-normalized matrices
# ==============================================================================

cluster_log2_qn_objects <- lapply(
  significant_cluster_ids,
  function(cluster_id_current) {
    input_file <- log2_qn_files[
      cluster_id_current
    ]

    message(
      "Reading log2 + QN matrix for cluster ",
      cluster_id_current,
      ": ",
      basename(input_file)
    )

    read_log2_qn_cluster_matrix(
      input_file = input_file,
      expected_sample_ids = sample_order,
      cluster_id = cluster_id_current
    )
  }
)

names(cluster_log2_qn_objects) <- significant_cluster_ids

cluster_log2_qn_matrices <- lapply(
  cluster_log2_qn_objects,
  function(x) {
    x$matrix
  }
)

# ==============================================================================
# 11. Build one cluster-gene x sample log2 + QN matrix
# ==============================================================================

sample_level_log2_qn_matrix <- matrix(
  NA_real_,
  nrow = nrow(significant_interaction),
  ncol = length(sample_order),
  dimnames = list(
    NULL,
    sample_order
  )
)

row_ids <- paste0(
  ifelse(
    is.na(significant_interaction$gene) |
      significant_interaction$gene == "",
    significant_interaction$ensembl_gene_id,
    significant_interaction$gene
  ),
  " [C",
  significant_interaction$cluster_id,
  "]"
)

# The cluster tag makes repeated gene symbols unambiguous without changing
# the official gene symbol itself, e.g. Ttr [C14] versus Ttr [C15].
if (anyDuplicated(row_ids) > 0L) {
  duplicated_labels <- unique(
    row_ids[
      duplicated(row_ids)
    ]
  )

  stop(
    "Duplicated gene-cluster row labels detected: ",
    paste(duplicated_labels, collapse = ", "),
    call. = FALSE
  )
}

rownames(sample_level_log2_qn_matrix) <- row_ids

for (row_index in seq_len(nrow(significant_interaction))) {
  current_cluster <-
    significant_interaction$cluster_id[[row_index]]

  current_gene_id <-
    significant_interaction$ensembl_gene_id[[row_index]]

  current_log2_qn_matrix <-
    cluster_log2_qn_matrices[[current_cluster]]

  gene_match <- match(
    current_gene_id,
    rownames(current_log2_qn_matrix)
  )

  if (is.na(gene_match)) {
    stop(
      "Gene ",
      current_gene_id,
      " is missing from the log2 + QN matrix for cluster ",
      current_cluster,
      ".",
      call. = FALSE
    )
  }

  sample_level_log2_qn_matrix[
    row_index,
    sample_order
  ] <- current_log2_qn_matrix[
    gene_match,
    sample_order,
    drop = TRUE
  ]
}

if (anyNA(sample_level_log2_qn_matrix)) {
  stop(
    "The assembled log2 + QN heatmap matrix contains NA values.",
    call. = FALSE
  )
}

if (any(!is.finite(sample_level_log2_qn_matrix))) {
  stop(
    "The assembled log2 + QN heatmap matrix contains non-finite values.",
    call. = FALSE
  )
}

# ==============================================================================
# 12. Row-wise Z-score for visualization
# ==============================================================================

# IMPORTANT:
# sample_level_log2_qn_matrix is already:
#
#   raw count -> log2(count + 1) -> quantile normalization
#
# Therefore no additional log transformation or quantile normalization is
# performed here.

zscore_matrix <- calculate_row_zscores(
  expression_matrix = sample_level_log2_qn_matrix,
  lower_limit = zscore_lower_limit,
  upper_limit = zscore_upper_limit
)

# ==============================================================================

# 13. Row annotations

# ==============================================================================

# Repetition is purely an annotation:

# number of significant interaction clusters containing the same gene.

gene_repetition_table <- significant_interaction |>

  dplyr::count(

    .data$ensembl_gene_id,

name = "n_repetition"

  )

row_annotation_table <- significant_interaction |>

  dplyr::transmute(

row_id = row_ids,

ensembl_gene_id = .data$ensembl_gene_id,

cluster_id = .data$cluster_id,

cluster_label = vapply(

      .data$cluster_id,

      format_cluster_label,

FUN.VALUE = character(1)

    ),

interaction_direction = .data$interaction_direction,

logFC = .data$logFC,

abs_logFC = abs(.data$logFC),

PValue = .data$PValue,

FDR = .data$FDR,

max_mean_percent = .data$max_mean_percent,

fdr_symbol = dplyr::case_when(

      .data$FDR < 0.001 ~ "***",

      .data$FDR < 0.01 ~ "**",

      .data$FDR < 0.05 ~ "*",

      TRUE ~ ""

    )

  ) |>

  dplyr::left_join(

    gene_repetition_table,

by = "ensembl_gene_id"

  )

cluster_factor_levels <- as.character(

  significant_clusters

)

cluster_annotation_colors <- custom_cluster_colors[

  cluster_factor_levels

]

names(

  cluster_annotation_colors

) <- cluster_factor_levels

# ------------------------------------------------------------------------------

# 13a. Source-cluster legend: number of significant genes per cluster

# ------------------------------------------------------------------------------

cluster_result_counts <- table(

  factor(

    row_annotation_table$cluster_id,

levels = cluster_factor_levels

  )

)

cluster_legend_labels <- paste0(

  "C",

  cluster_factor_levels,

  " (n = ",

  as.integer(

    cluster_result_counts

  ),

  "): ",

  unname(

    custom_cluster_labels[

      cluster_factor_levels

    ]

  )

)

# ------------------------------------------------------------------------------

# 13b. Max mean % annotation with values printed inside purple cells

# ------------------------------------------------------------------------------

max_mean_percent_color_function <- circlize::colorRamp2(

  c(

    0,

    25,

    50,

    75,

    100

  ),

  c(

    "#F3E8FF",

    "#DDD6FE",

    "#C084FC",

    "#9333EA",

    "#4C1D95"

  )

)

annotation_row_font_size <- if (

  nrow(

    row_annotation_table

  ) <= 80L

) {

  7.2

} else if (

  nrow(

    row_annotation_table

  ) <= 180L

) {

  5.8

} else {

  4.8

}

percent_annotation_labels <- ifelse(

  is.na(

    row_annotation_table$max_mean_percent

  ),

  "NA",

  paste0(

    round(

      row_annotation_table$max_mean_percent

    ),

    "%"

  )

)

percent_annotation_text_colour <- ifelse(

  !is.na(

    row_annotation_table$max_mean_percent

  ) &

    row_annotation_table$max_mean_percent > 70,

  "white",

  "black"

)

max_mean_percent_annotation <- ComplexHeatmap::AnnotationFunction(

which = "row",

width = grid::unit(

    5.5,

    "mm"

  ),

var_import = list(

values = row_annotation_table$max_mean_percent,

labels = percent_annotation_labels,

text_colours = percent_annotation_text_colour,

colour_function = max_mean_percent_color_function,

font_size = annotation_row_font_size

  ),

fun = function(index) {

    n_current <- length(

      index

    )

    y_positions <- 1 - (

      seq_len(

        n_current

      ) - 0.5

    ) / n_current

    current_values <- values[

      index

    ]

    current_labels <- labels[

      index

    ]

    current_text_colours <- text_colours[

      index

    ]

    current_fill <- rep(

      "grey90",

      n_current

    )

    finite_values <- is.finite(

      current_values

    )

    current_fill[

      finite_values

    ] <- colour_function(

      current_values[

        finite_values

      ]

    )

    grid::grid.rect(

x = grid::unit(

        0.5,

        "npc"

      ),

y = grid::unit(

        y_positions,

        "npc"

      ),

width = grid::unit(

        1,

        "npc"

      ),

height = grid::unit(

        1 / n_current,

        "npc"

      ),

gp = grid::gpar(

fill = current_fill,

col = NA

      )

    )

    grid::grid.text(

label = current_labels,

x = grid::unit(

        0.5,

        "npc"

      ),

y = grid::unit(

        y_positions,

        "npc"

      ),

gp = grid::gpar(

col = current_text_colours,

fontsize = font_size,

fontface = "bold"

      )

    )

  }

)

# ------------------------------------------------------------------------------

# 13c. Complete right-side annotation block

# ------------------------------------------------------------------------------

gene_symbol_annotation_width <-

  ComplexHeatmap::max_text_width(

    row_ids,

gp = grid::gpar(

fontsize = annotation_row_font_size

    )

  ) +

  grid::unit(

    1.5,

    "mm"

  )

right_annotation <- ComplexHeatmap::rowAnnotation(

Cluster = row_annotation_table$cluster_id,

  `Max mean %` = max_mean_percent_annotation,

  `FDR marker` = ComplexHeatmap::anno_text(

    row_annotation_table$fdr_symbol,

location = grid::unit(

      0.5,

      "npc"

    ),

just = "center",

gp = grid::gpar(

col = "black",

fontsize = annotation_row_font_size,

fontface = "bold"

    ),

width = grid::unit(

      3.5,

      "mm"

    )

  ),

  `n repetition` = ComplexHeatmap::anno_text(

    as.character(

      row_annotation_table$n_repetition

    ),

location = grid::unit(

      0.5,

      "npc"

    ),

just = "center",

gp = grid::gpar(

col = "black",

fontsize = annotation_row_font_size,

fontface = "bold"

    ),

width = grid::unit(

      4.0,

      "mm"

    )

  ),

  `Gene symbol` = ComplexHeatmap::anno_text(

    row_ids,

location = grid::unit(

      0,

      "npc"

    ),

just = "left",

gp = grid::gpar(

col = "black",

fontsize = annotation_row_font_size

    ),

width = gene_symbol_annotation_width

  ),

col = list(

Cluster = cluster_annotation_colors

  ),

simple_anno_size = grid::unit(

    4.2,

    "mm"

  ),

gap = grid::unit(

    0.60,

    "mm"

  ),

show_annotation_name = c(

Cluster = TRUE,

    `Max mean %` = TRUE,

    `FDR marker` = TRUE,

    `n repetition` = TRUE,

    `Gene symbol` = TRUE

  ),

annotation_label = c(

Cluster = "Cluster",

    `Max mean %` = "Max mean %",

    `FDR marker` = "FDR marker",

    `n repetition` = "n repetition",

    `Gene symbol` = "Gene symbol"

  ),

annotation_name_side = "top",

annotation_name_rot = c(

Cluster = 45,

    `Max mean %` = 45,

    `FDR marker` = 45,

    `n repetition` = 45,

    `Gene symbol` = 45

  ),

annotation_name_gp = grid::gpar(

fontsize = 8,

fontface = "bold"

  ),

annotation_legend_param = list(

Cluster = list(

title = "Source cluster",

labels = cluster_legend_labels,

at = cluster_factor_levels,

ncol = 1

    )

  )

)

# One global row clustering and one dendrogram.

# The source cluster modifies pairwise distances but does not impose a hard split.

hc_rows <- create_modified_row_hclust(

expression_matrix = zscore_matrix,

annotation_table = row_annotation_table,

base_distance = row_clustering_base_distance,

r_scale = row_clustering_r_scale,

same_cluster_bonus = row_clustering_same_cluster_bonus,

linkage_method = row_clustering_linkage_method

)

# ==============================================================================

# 14. Column annotation

# ==============================================================================

sex_colors <- c(

  "Male" = "#3C78D8",

  "Female" = "#C050A0"

)

group_colors <- c(

  "Neurotypical" = "grey70",

  "ASD" = "black"

)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(

Group = sample_metadata$fmt_donor_group,

Sex = sample_metadata$sex,

col = list(

Group = group_colors,

Sex = sex_colors

  ),

simple_anno_size = grid::unit(4.2, "mm"),

annotation_name_side = "left",

annotation_name_gp = grid::gpar(

fontsize = 9,

fontface = "bold"

  ),

annotation_legend_param = list(

Group = list(

title = "FMT donor group",

ncol = 1

    ),

Sex = list(

title = "Sex",

ncol = 1

    )

  )

)

# ==============================================================================
# 15. Heatmap title and dimensions
# ==============================================================================

number_of_results <- nrow(significant_interaction)
number_of_unique_genes <- length(unique(significant_interaction$gene))
number_of_clusters <- length(significant_cluster_ids)
number_of_samples <- ncol(zscore_matrix)

heatmap_title_lines <- c(
  "Group x Sex interaction: significant genes per cluster",
  paste0(
    "Selection: FDR < ",
    fdr_threshold,
    " | ",
    number_of_results,
    " cluster-gene results | ",
    number_of_unique_genes,
    " unique genes | ",
    number_of_clusters,
    " clusters"
  ),
  "Input: log2(count + 1) + quantile normalization",
  "Heatmap transformation: row Z-score of precomputed log2 + QN values, clipped to [-2, 2]",
  paste0(
    "Columns: Male Neurotypical -> Male ASD -> Female Neurotypical -> Female ASD",
    " | ",
    number_of_samples,
    " samples"
  ),
  paste0(
    "Row clustering: Pearson correlation on row Z-scores",
    " | base distance d = ",
    row_clustering_base_distance,
    " | r_scale = ",
    row_clustering_r_scale
  ),
  paste0(
    "Same-cluster bonus = ",
    row_clustering_same_cluster_bonus,
    " -> same cluster: d x ",
    formatC(
      1 - row_clustering_same_cluster_bonus,
      format = "f",
      digits = 2
    ),
    ", different clusters: d x ",
    formatC(
      1 + row_clustering_same_cluster_bonus,
      format = "f",
      digits = 2
    ),
    " | linkage = ",
    row_clustering_linkage_method
  ),
  clustering_title
)

heatmap_title <- paste(
  heatmap_title_lines,
  collapse = "\n"
)

row_font_size <- if (number_of_results <= 80L) {
  7.2
} else if (number_of_results <= 180L) {
  5.8
} else {
  4.8
}

row_height_mm <- if (number_of_results <= 80L) {
  3.1
} else if (number_of_results <= 180L) {
  2.35
} else {
  1.9
}

heatmap_body_height_mm <- number_of_results * row_height_mm
heatmap_body_width_mm <- number_of_samples * 6.2

pdf_width_inches <- 15

pdf_height_inches <- max(
  14,
  min(
    32,
    (heatmap_body_height_mm + 95) / 25.4
  )
)

# ==============================================================================

# 16. Create heatmap

# ==============================================================================

heatmap_color_function <- circlize::colorRamp2(

  c(-2, 0, 2),

  c("#2166AC", "white", "#B2182B")

)

heatmap_object <- ComplexHeatmap::Heatmap(

  zscore_matrix,

name = "Row Z-score",

col = heatmap_color_function,

top_annotation = top_annotation,

right_annotation = right_annotation,

cluster_columns = FALSE,

column_split = sample_metadata$combined_group,

cluster_column_slices = FALSE,

column_gap = grid::unit(1.2, "mm"),

cluster_rows = hc_rows,

row_dend_reorder = row_annotation_table$logFC,

row_gap = grid::unit(0, "mm"),

row_dend_width = grid::unit(18, "mm"),

show_column_names = TRUE,

column_names_side = "top",

column_names_rot = 45,

column_names_centered = TRUE,

column_names_gp = grid::gpar(

fontsize = 9,

fontface = "bold"

  ),

column_names_max_height = grid::unit(18, "mm"),

show_row_names = FALSE,

row_title_gp = grid::gpar(

fontsize = 8,

fontface = "bold"

  ),

na_col = "grey90",

rect_gp = grid::gpar(col = NA),

border = FALSE,

use_raster = FALSE,

width = grid::unit(heatmap_body_width_mm, "mm"),

height = grid::unit(heatmap_body_height_mm, "mm"),

column_title = heatmap_title,

column_title_side = "top",

column_title_gp = grid::gpar(

fontsize = 12,

fontface = "bold",

lineheight = 1.12

  ),

heatmap_legend_param = list(

title = "Row Z-score",

at = c(-2, -1, 0, 1, 2),

legend_height = grid::unit(46, "mm")

  )

)

# ==============================================================================
# 17. Export TSV files
# ==============================================================================

readr::write_tsv(
  significant_interaction,
  selected_results_tsv,
  na = "NA"
)

readr::write_tsv(
  row_annotation_table,
  row_annotations_tsv,
  na = "NA"
)

log2_qn_output <- tibble::tibble(
  row_id = rownames(sample_level_log2_qn_matrix)
) |>
  dplyr::bind_cols(
    tibble::as_tibble(
      as.data.frame(
        sample_level_log2_qn_matrix,
        check.names = FALSE
      ),
      .name_repair = "minimal"
    )
  )

readr::write_tsv(
  log2_qn_output,
  log2_qn_selected_matrix_tsv,
  na = "NA"
)

zscore_output <- tibble::tibble(
  row_id = rownames(zscore_matrix)
) |>
  dplyr::bind_cols(
    tibble::as_tibble(
      as.data.frame(
        zscore_matrix,
        check.names = FALSE
      ),
      .name_repair = "minimal"
    )
  )

readr::write_tsv(
  zscore_output,
  zscore_matrix_tsv,
  na = "NA"
)

sample_order_output <- sample_metadata |>
  dplyr::transmute(
    sample_order = .data$sample_order,
    sample_ID = .data$sample_ID,
    sex = as.character(.data$sex),
    fmt_donor_group = as.character(.data$fmt_donor_group),
    combined_group = as.character(.data$combined_group)
  )

readr::write_tsv(
  sample_order_output,
  sample_order_tsv,
  na = "NA"
)

summary_table <- tibble::tibble(
  dataset_name = dataset_name,
  clustering_name = clustering_name,
  interaction_test_id = interaction_test_id,
  fdr_threshold = fdr_threshold,
  significant_cluster_gene_results = number_of_results,
  unique_genes = number_of_unique_genes,
  significant_clusters = paste(
    significant_cluster_ids,
    collapse = ", "
  ),
  number_of_significant_clusters = number_of_clusters,
  samples_displayed = number_of_samples,
  matrix_input = paste0(
    "per-cluster raw pseudobulk counts transformed as ",
    "log2(count + 1) and quantile normalized within cluster"
  ),
  expression_input_directory = log2_qn_input_dir,
  additional_log2_in_heatmap_script = FALSE,
  additional_quantile_normalization_in_heatmap_script = FALSE,
  transformation_for_heatmap = "row Z-score -> clip [-2, 2]",
  row_clustering_base_distance = row_clustering_base_distance,
  row_clustering_r_scale = row_clustering_r_scale,
  row_clustering_same_cluster_bonus = row_clustering_same_cluster_bonus,
  row_clustering_linkage_method = row_clustering_linkage_method,
  output_pdf = heatmap_pdf
)

readr::write_tsv(
  summary_table,
  summary_tsv,
  na = "NA"
)

# ==============================================================================

# 18. Save PDF

# ==============================================================================

fdr_marker_legend <- make_fdr_marker_legend(

  significant_interaction

)

message("\nLegend threshold counts:")

message(

  "FDR < 0.05: n = ",

  sum(

    !is.na(significant_interaction$FDR) &

      significant_interaction$FDR < 0.05

  )

)

message(

  "FDR < 0.01: n = ",

  sum(

    !is.na(significant_interaction$FDR) &

      significant_interaction$FDR < 0.01

  )

)

message(

  "FDR < 0.001: n = ",

  sum(

    !is.na(significant_interaction$FDR) &

      significant_interaction$FDR < 0.001

  )

)

message("\nWriting heatmap PDF:")

message(heatmap_pdf)

message(

  "PDF dimensions: ",

  round(pdf_width_inches, 2),

  " x ",

  round(pdf_height_inches, 2),

  " inches"

)

grDevices::pdf(

file = heatmap_pdf,

width = pdf_width_inches,

height = pdf_height_inches,

onefile = TRUE,

useDingbats = FALSE

)

ComplexHeatmap::draw(

  heatmap_object,

heatmap_legend_side = "right",

annotation_legend_side = "right",

merge_legends = TRUE,

legend_grouping = "original",

annotation_legend_list = list(

    fdr_marker_legend

  ),

padding = grid::unit(

    c(3, 3, 3, 3),

    "mm"

  )

)

grDevices::dev.off()

# ==============================================================================

# 19. Validate outputs and report

# ==============================================================================

expected_output_files <- c(
  heatmap_pdf,
  selected_results_tsv,
  row_annotations_tsv,
  log2_qn_selected_matrix_tsv,
  zscore_matrix_tsv,
  sample_order_tsv,
  summary_tsv
)

missing_output_files <- expected_output_files[

  !file.exists(expected_output_files)

]

if (length(missing_output_files) > 0L) {

  stop(

    "Missing output file(s):\n",

    paste(missing_output_files, collapse = "\n"),

call. = FALSE

  )

}

empty_output_files <- expected_output_files[

  file.info(expected_output_files)$size == 0

]

if (length(empty_output_files) > 0L) {

  stop(

    "Empty output file(s):\n",

    paste(empty_output_files, collapse = "\n"),

call. = FALSE

  )

}

message("\nHeatmap completed successfully.")
message("Expression input directory:")
message(normalizePath(log2_qn_input_dir, mustWork = TRUE))
message("\nOutput directory:")
message(normalizePath(heatmap_output_dir, mustWork = TRUE))
message("\nHeatmap PDF:")
message(normalizePath(heatmap_pdf, mustWork = TRUE))

# ==============================================================================

# End

# ==============================================================================
