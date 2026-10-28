#!/usr/bin/env Rscript

# ==============================================================================
# 03_heatmap_effectDonorGroup_FDR005_crescendoCorrected_log2QN_v2.R
#
# Exploratory long heatmap for the main FMT donor-group effect after Crescendo
# batch correction.
#
# STATISTICAL SELECTION
#   Main donor-group test:
#     Overall_Group_ASD_vs_Neurotypical
#
#   Selection:
#     FDR < 0.05
#
#   IMPORTANT:
#     There is NO |log2FC| threshold in this heatmap.
#
# EXPRESSION VALUES
#   Crescendo-corrected pseudobulk counts, prepared by the standalone edgeR
#   runner:
#
#     corrected pseudobulk counts
#       -> log2(count + 1)
#       -> quantile normalization across the 16 biological samples
#          separately within each cluster
#       -> row-wise Z-score
#       -> clip to [-2, 2]
#
#   Therefore the heatmap expression values are based on the data AFTER
#   Crescendo correction and do not use the old uncorrected pseudobulk counts.
#
# ROWS
#   One row = one significant cluster-gene result.
#
#   The same gene may therefore occur more than once when significant in
#   several anatomical clusters, e.g.:
#     GeneA [C3]
#     GeneA [C13]
#
# COLUMNS
#   Samples are ordered:
#     Male Neurotypical
#       -> Male ASD
#       -> Female Neurotypical
#       -> Female ASD
#
# RIGHT ANNOTATIONS
#   Cluster
#   FDR marker:
#     *   FDR < 0.05
#     **  FDR < 0.01
#     *** FDR < 0.001
#   log2FC
#   Direction
#   n clusters
#   Gene symbol
#
# ROW CLUSTERING
#   Base distance:
#     2 - Pearson correlation
#
#   Rows from the same anatomical cluster are made slightly closer and rows
#   from different clusters slightly more distant.
#
# OUTPUT
#   1. One long PDF heatmap.
#   2. One TSV containing all displayed cluster-gene results.
#
# NOTE
#   This is an exploratory visualization. edgeR inference itself was performed
#   on Crescendo-corrected count data and is not performed in this script.
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

donor_group_test_id <- "Overall_Group_ASD_vs_Neurotypical"

fdr_threshold <- 0.05

zscore_lower_limit <- -2
zscore_upper_limit <- 2

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
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics",
  "02_edgeR"
)

analysis_prefix <- paste0(
  dataset_name,
  "_",
  clustering_name,
  "_crescendoCorrectedPseudobulk_edgeR_twoFactorInteraction"
)

edgeR_results_rdata_file <- file.path(
  statistics_dir,
  "03_edgeRResults",
  paste0(
    analysis_prefix,
    "_edgeRResults.RData"
  )
)

edgeR_model_objects_rdata_file <- file.path(
  statistics_dir,
  "05_objects",
  paste0(
    analysis_prefix,
    "_edgeRModelObjects.RData"
  )
)

log2_qn_input_dir <- file.path(
  statistics_dir,
  "01_pseudobulkInputs",
  "perCluster_crescendoCorrectedCounts_log2plus1_QN"
)

heatmap_output_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,
  "02_heatmaps",
  "02_mainEffect_donorGroup",
  "01_FDR005_crescendoCorrected"
)

dir.create(
  heatmap_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

heatmap_pdf_file <- file.path(
  heatmap_output_dir,
  paste0(
    "01_",
    dataset_name,
    "_effectDonorGroup_FDR005_",
    "crescendoCorrected_log2QN_heatmapLong.pdf"
  )
)

displayed_results_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "02_",
    dataset_name,
    "_effectDonorGroup_FDR005_",
    "crescendoCorrected_displayedResults.tsv"
  )
)

required_input_files <- c(
  edgeR_results_rdata_file,
  edgeR_model_objects_rdata_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0L) {
  stop(
    "Missing required RData file(s):\n",
    paste(missing_input_files, collapse = "\n"),
    call. = FALSE
  )
}

if (!dir.exists(log2_qn_input_dir)) {
  stop(
    "Missing Crescendo-corrected log2+QN input directory:\n",
    log2_qn_input_dir,
    call. = FALSE
  )
}

# ==============================================================================
# 5. Helper functions
# ==============================================================================

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
    is.finite(output_matrix) &
      output_matrix < lower_limit
  ] <- lower_limit

  output_matrix[
    is.finite(output_matrix) &
      output_matrix > upper_limit
  ] <- upper_limit

  output_matrix
}

modified_correlation_distance <- function(
    expression_matrix,
    cluster_labels,
    base_distance = "2 - r",
    r_scale = 1,
    same_cluster_bonus = 0.20
) {
  if (!is.matrix(expression_matrix)) {
    expression_matrix <- as.matrix(expression_matrix)
  }

  cluster_labels <- as.character(cluster_labels)

  if (length(cluster_labels) != nrow(expression_matrix)) {
    stop(
      "Number of cluster labels must equal number of heatmap rows.",
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
      "`base_distance` must be '2 - r' or '1 - r'.",
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
    base_dist[same_cluster_matrix] *
      (1 - same_cluster_bonus)

  dist_mat[!same_cluster_matrix] <-
    base_dist[!same_cluster_matrix] *
      (1 + same_cluster_bonus)

  diag(dist_mat) <- 0
  dist_mat[dist_mat < 0] <- 0

  stats::as.dist(dist_mat)
}

create_modified_row_hclust <- function(
    expression_matrix,
    annotation_table
) {
  matched_annotation <- annotation_table[
    match(
      rownames(expression_matrix),
      annotation_table$row_id
    ),
    ,
    drop = FALSE
  ]

  if (
    nrow(matched_annotation) != nrow(expression_matrix) ||
      anyNA(matched_annotation$row_id) ||
      anyNA(matched_annotation$cluster_id)
  ) {
    stop(
      "Failed to match row annotations to heatmap rows.",
      call. = FALSE
    )
  }

  modified_distance <- modified_correlation_distance(
    expression_matrix = expression_matrix,
    cluster_labels = matched_annotation$cluster_id,
    base_distance = row_clustering_base_distance,
    r_scale = row_clustering_r_scale,
    same_cluster_bonus = row_clustering_same_cluster_bonus
  )

  stats::hclust(
    modified_distance,
    method = row_clustering_linkage_method
  )
}

find_log2_qn_file_for_cluster <- function(
    cluster_id_current
) {
  candidate_file <- file.path(
    log2_qn_input_dir,
    paste0(
      analysis_prefix,
      "_cluster_",
      cluster_id_current,
      "_crescendoCorrectedCounts_log2plus1_QN.tsv.gz"
    )
  )

  if (!file.exists(candidate_file)) {
    stop(
      "Missing Crescendo-corrected log2+QN matrix for cluster ",
      cluster_id_current,
      ":\n",
      candidate_file,
      call. = FALSE
    )
  }

  candidate_file
}

# ==============================================================================
# 6. Load edgeR objects
# ==============================================================================

message("Loading edgeR results:")
message(edgeR_results_rdata_file)

load(edgeR_results_rdata_file)

message("\nLoading edgeR model objects:")
message(edgeR_model_objects_rdata_file)

load(edgeR_model_objects_rdata_file)

required_objects <- c(
  "edgeR_perCluster_combinedResults",
  "edgeR_perCluster_testDefinitions",
  "edgeR_perCluster_modelObjects"
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
    "Required object(s) missing after loading RData: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

# ==============================================================================
# 7. Validate and select donor-group main-effect results
# ==============================================================================

if (
  !donor_group_test_id %in%
    edgeR_perCluster_testDefinitions$test_id
) {
  stop(
    "Main donor-group test not found in test definitions: ",
    donor_group_test_id,
    call. = FALSE
  )
}

if (
  !donor_group_test_id %in%
    names(edgeR_perCluster_combinedResults)
) {
  stop(
    "Main donor-group test not found in combined edgeR results: ",
    donor_group_test_id,
    call. = FALSE
  )
}

donor_group_definition <- edgeR_perCluster_testDefinitions |>
  dplyr::filter(
    .data$test_id == donor_group_test_id
  )

selected_test_comparison <- as.character(
  donor_group_definition$comparison[[1]]
)

donor_group_results <- edgeR_perCluster_combinedResults[[donor_group_test_id]] |>
  tibble::as_tibble()

required_result_columns <- c(
  "cluster_id",
  "ensembl_gene_id",
  "gene",
  "logFC",
  "PValue",
  "FDR"
)

missing_result_columns <- setdiff(
  required_result_columns,
  colnames(donor_group_results)
)

if (length(missing_result_columns) > 0L) {
  stop(
    "Selected edgeR result table is missing required columns: ",
    paste(missing_result_columns, collapse = ", "),
    call. = FALSE
  )
}

significant_results <- donor_group_results |>
  dplyr::filter(
    !is.na(.data$FDR),
    !is.na(.data$logFC),
    .data$FDR < fdr_threshold
  ) |>
  dplyr::mutate(
    cluster_id = as.character(.data$cluster_id),
    abs_logFC = abs(.data$logFC),
    effect_direction = dplyr::case_when(
      .data$logFC > 0 ~ "UP",
      .data$logFC < 0 ~ "DOWN",
      TRUE ~ "ZERO"
    ),
    fdr_symbol = dplyr::case_when(
      .data$FDR < 0.001 ~ "***",
      .data$FDR < 0.01 ~ "**",
      .data$FDR < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) |>
  dplyr::group_by(
    .data$ensembl_gene_id
  ) |>
  dplyr::mutate(
    n_clusters_for_gene =
      dplyr::n_distinct(.data$cluster_id)
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(
    as.integer(.data$cluster_id),
    .data$FDR,
    dplyr::desc(.data$abs_logFC)
  )

if (nrow(significant_results) == 0L) {
  stop(
    "No donor-group result satisfies FDR < ",
    fdr_threshold,
    ".",
    call. = FALSE
  )
}

significant_results <- significant_results |>
  dplyr::mutate(
    row_id = paste0(
      dplyr::if_else(
        is.na(.data$gene) |
          .data$gene == "",
        .data$ensembl_gene_id,
        .data$gene
      ),
      " [C",
      .data$cluster_id,
      "]"
    )
  )

if (anyDuplicated(significant_results$row_id)) {
  duplicated_rows <- unique(
    significant_results$row_id[
      duplicated(significant_results$row_id)
    ]
  )

  stop(
    "Duplicated gene-cluster row IDs: ",
    paste(duplicated_rows, collapse = ", "),
    call. = FALSE
  )
}

message("\nSelected donor-group test:")
message(donor_group_test_id)

message("Comparison:")
message(selected_test_comparison)

message("\nFDR < 0.05 cluster-gene results: ", nrow(significant_results))

message(
  "Unique genes: ",
  length(unique(significant_results$ensembl_gene_id))
)

message(
  "Clusters represented: ",
  length(unique(significant_results$cluster_id))
)

# ==============================================================================
# 8. Extract common sample metadata and define column order
# ==============================================================================

significant_cluster_ids <- sort(
  unique(
    as.integer(significant_results$cluster_id)
  )
)

significant_cluster_ids <- as.character(
  significant_cluster_ids
)

missing_model_clusters <- setdiff(
  significant_cluster_ids,
  names(edgeR_perCluster_modelObjects)
)

if (length(missing_model_clusters) > 0L) {
  stop(
    "Model RData is missing significant cluster(s): ",
    paste(missing_model_clusters, collapse = ", "),
    call. = FALSE
  )
}

sample_metadata_list <- lapply(
  significant_cluster_ids,
  function(cluster_id_current) {
    model_object <- edgeR_perCluster_modelObjects[[cluster_id_current]]

    if (is.null(model_object$dge)) {
      stop(
        "DGEList missing for cluster ",
        cluster_id_current,
        ".",
        call. = FALSE
      )
    }

    model_object$dge$samples |>
      tibble::rownames_to_column(
        "dge_rowname"
      ) |>
      tibble::as_tibble() |>
      dplyr::transmute(
        sample_ID = as.character(.data$sample_ID),
        fmt_donor_group =
          as.character(.data$fmt_donor_group),
        sex = as.character(.data$sex)
      )
  }
)

sample_metadata <- sample_metadata_list |>
  dplyr::bind_rows() |>
  dplyr::distinct(
    .data$sample_ID,
    .data$fmt_donor_group,
    .data$sex
  )

metadata_conflicts <- sample_metadata |>
  dplyr::count(
    .data$sample_ID,
    name = "n_rows"
  ) |>
  dplyr::filter(
    .data$n_rows > 1L
  )

if (nrow(metadata_conflicts) > 0L) {
  stop(
    "Conflicting metadata found for sample(s): ",
    paste(
      metadata_conflicts$sample_ID,
      collapse = ", "
    ),
    call. = FALSE
  )
}

sample_metadata <- sample_metadata |>
  dplyr::mutate(
    sex = factor(
      .data$sex,
      levels = c(
        "Male",
        "Female"
      )
    ),
    fmt_donor_group = factor(
      .data$fmt_donor_group,
      levels = c(
        "Neurotypical",
        "ASD"
      )
    )
  ) |>
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

sample_order <- as.character(
  sample_metadata$sample_ID
)

if (length(sample_order) != 16L) {
  stop(
    "Expected 16 biological samples but found ",
    length(sample_order),
    ".",
    call. = FALSE
  )
}

message("\nSample order:")

print(
  sample_metadata |>
    dplyr::select(
      "sample_order",
      "sample_ID",
      "sex",
      "fmt_donor_group"
    ),
  n = Inf
)

# ==============================================================================
# 9. Load Crescendo-corrected log2+QN matrices
# ==============================================================================

cluster_log2_qn_matrices <- lapply(
  significant_cluster_ids,
  function(cluster_id_current) {
    input_file <- find_log2_qn_file_for_cluster(
      cluster_id_current
    )

    message(
      "Loading corrected log2+QN matrix for cluster ",
      cluster_id_current,
      ": ",
      basename(input_file)
    )

    input_table <- readr::read_tsv(
      input_file,
      show_col_types = FALSE,
      progress = FALSE
    ) |>
      tibble::as_tibble()

    required_annotation_columns <- c(
      "ensembl_gene_id",
      "gene",
      "chromosome"
    )

    missing_annotation_columns <- setdiff(
      required_annotation_columns,
      colnames(input_table)
    )

    if (length(missing_annotation_columns) > 0L) {
      stop(
        "Cluster ",
        cluster_id_current,
        " matrix is missing annotation column(s): ",
        paste(missing_annotation_columns, collapse = ", "),
        call. = FALSE
      )
    }

    missing_sample_columns <- setdiff(
      sample_order,
      colnames(input_table)
    )

    if (length(missing_sample_columns) > 0L) {
      stop(
        "Cluster ",
        cluster_id_current,
        " matrix is missing sample column(s): ",
        paste(missing_sample_columns, collapse = ", "),
        call. = FALSE
      )
    }

    if (anyDuplicated(input_table$ensembl_gene_id)) {
      stop(
        "Duplicated Ensembl IDs in corrected log2+QN matrix for cluster ",
        cluster_id_current,
        ".",
        call. = FALSE
      )
    }

    expression_matrix <- as.matrix(
      input_table[
        ,
        sample_order,
        drop = FALSE
      ]
    )

    storage.mode(expression_matrix) <- "double"

    rownames(expression_matrix) <-
      as.character(
        input_table$ensembl_gene_id
      )

    if (anyNA(expression_matrix)) {
      stop(
        "NA values found in corrected log2+QN matrix for cluster ",
        cluster_id_current,
        ".",
        call. = FALSE
      )
    }

    if (any(!is.finite(expression_matrix))) {
      stop(
        "Non-finite values found in corrected log2+QN matrix for cluster ",
        cluster_id_current,
        ".",
        call. = FALSE
      )
    }

    expression_matrix
  }
)

names(cluster_log2_qn_matrices) <-
  significant_cluster_ids

# ==============================================================================
# 10. Build cluster-gene x sample expression matrix
# ==============================================================================

sample_level_log2_qn_matrix <- matrix(
  NA_real_,
  nrow = nrow(significant_results),
  ncol = length(sample_order),
  dimnames = list(
    significant_results$row_id,
    sample_order
  )
)

for (row_index in seq_len(nrow(significant_results))) {
  current_cluster <-
    significant_results$cluster_id[[row_index]]

  current_gene_id <-
    significant_results$ensembl_gene_id[[row_index]]

  current_matrix <-
    cluster_log2_qn_matrices[[current_cluster]]

  gene_match <- match(
    current_gene_id,
    rownames(current_matrix)
  )

  if (is.na(gene_match)) {
    stop(
      "Gene ",
      current_gene_id,
      " missing from corrected log2+QN matrix for cluster ",
      current_cluster,
      ".",
      call. = FALSE
    )
  }

  sample_level_log2_qn_matrix[
    row_index,

  ] <-
    current_matrix[
      gene_match,
      sample_order,
      drop = TRUE
    ]
}

if (anyNA(sample_level_log2_qn_matrix)) {
  stop(
    "NA values found in assembled cluster-gene expression matrix.",
    call. = FALSE
  )
}

# ==============================================================================
# 11. Row Z-score
# ==============================================================================

zscore_matrix <- calculate_row_zscores(
  expression_matrix =
    sample_level_log2_qn_matrix,
  lower_limit =
    zscore_lower_limit,
  upper_limit =
    zscore_upper_limit
)

# ==============================================================================
# 12. Row clustering
# ==============================================================================

row_annotation_table <- significant_results |>
  dplyr::transmute(
    row_id = .data$row_id,
    cluster_id = .data$cluster_id,
    logFC = .data$logFC,
    FDR = .data$FDR,
    fdr_symbol = .data$fdr_symbol,
    effect_direction = .data$effect_direction,
    n_clusters_for_gene =
      .data$n_clusters_for_gene
  )

if (nrow(zscore_matrix) > 1L) {
  hc_rows <- create_modified_row_hclust(
    expression_matrix = zscore_matrix,
    annotation_table = row_annotation_table
  )

  cluster_rows_current <- hc_rows

  row_dend_reorder_current <-
    row_annotation_table$logFC
} else {
  cluster_rows_current <- FALSE
  row_dend_reorder_current <- FALSE
}

# ==============================================================================
# 13. Colours
# ==============================================================================

heatmap_color_function <- circlize::colorRamp2(
  c(
    -2,
    0,
    2
  ),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)

sex_colors <- c(
  "Male" = "#3C78D8",
  "Female" = "#C050A0"
)

group_colors <- c(
  "Neurotypical" = "grey70",
  "ASD" = "black"
)

direction_colors <- c(
  "UP" = "#B2182B",
  "DOWN" = "#2166AC",
  "ZERO" = "grey70"
)

# ==============================================================================
# 14. Top annotation
# ==============================================================================

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Group = sample_metadata$fmt_donor_group,
  Sex = sample_metadata$sex,
  col = list(
    Group = group_colors,
    Sex = sex_colors
  ),
  simple_anno_size = grid::unit(
    4.2,
    "mm"
  ),
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
# 15. Right annotations and legends
# ==============================================================================

cluster_factor_levels <- as.character(
  sort(
    unique(
      as.integer(
        significant_results$cluster_id
      )
    )
  )
)

cluster_annotation_colors <-
  custom_cluster_colors[
    cluster_factor_levels
  ]

names(cluster_annotation_colors) <-
  cluster_factor_levels

cluster_counts <- table(
  factor(
    significant_results$cluster_id,
    levels = cluster_factor_levels
  )
)

cluster_legend_labels <- paste0(
  "C",
  cluster_factor_levels,
  " (n = ",
  as.integer(cluster_counts),
  "): ",
  unname(
    custom_cluster_labels[
      cluster_factor_levels
    ]
  )
)

number_of_results <- nrow(
  significant_results
)

number_of_unique_genes <- length(
  unique(
    significant_results$ensembl_gene_id
  )
)

annotation_row_font_size <- if (
  number_of_results <= 100L
) {
  7.0
} else if (
  number_of_results <= 300L
) {
  5.6
} else if (
  number_of_results <= 700L
) {
  4.4
} else {
  3.6
}

logfc_labels <- sprintf(
  "%.2f",
  significant_results$logFC
)

gene_symbol_annotation_width <-
  ComplexHeatmap::max_text_width(
    significant_results$row_id,
    gp = grid::gpar(
      fontsize =
        annotation_row_font_size
    )
  ) +
  grid::unit(
    1.5,
    "mm"
  )

right_annotation <- ComplexHeatmap::rowAnnotation(
  Cluster =
    significant_results$cluster_id,

  `FDR marker` =
    ComplexHeatmap::anno_text(
      significant_results$fdr_symbol,
      location =
        grid::unit(
          0.5,
          "npc"
        ),
      just = "center",
      gp = grid::gpar(
        col = "black",
        fontsize =
          annotation_row_font_size,
        fontface = "bold"
      ),
      width =
        grid::unit(
          4.0,
          "mm"
        )
    ),

  log2FC =
    ComplexHeatmap::anno_text(
      logfc_labels,
      location =
        grid::unit(
          0.5,
          "npc"
        ),
      just = "center",
      gp = grid::gpar(
        col = "black",
        fontsize =
          annotation_row_font_size
      ),
      width =
        grid::unit(
          7.0,
          "mm"
        )
    ),

  Direction =
    significant_results$effect_direction,

  `n clusters` =
    ComplexHeatmap::anno_text(
      as.character(
        significant_results$n_clusters_for_gene
      ),
      location =
        grid::unit(
          0.5,
          "npc"
        ),
      just = "center",
      gp = grid::gpar(
        col = "black",
        fontsize =
          annotation_row_font_size,
        fontface = "bold"
      ),
      width =
        grid::unit(
          5.0,
          "mm"
        )
    ),

  `Gene symbol` =
    ComplexHeatmap::anno_text(
      significant_results$row_id,
      location =
        grid::unit(
          0,
          "npc"
        ),
      just = "left",
      gp = grid::gpar(
        col = "black",
        fontsize =
          annotation_row_font_size
      ),
      width =
        gene_symbol_annotation_width
    ),

  col = list(
    Cluster =
      cluster_annotation_colors,
    Direction =
      direction_colors
  ),

  simple_anno_size =
    grid::unit(
      4.2,
      "mm"
    ),

  gap =
    grid::unit(
      0.60,
      "mm"
    ),

  annotation_name_side = "top",

  annotation_name_rot = c(
    Cluster = 45,
    `FDR marker` = 45,
    log2FC = 45,
    Direction = 45,
    `n clusters` = 45,
    `Gene symbol` = 45
  ),

  annotation_name_gp =
    grid::gpar(
      fontsize = 8,
      fontface = "bold"
    ),

  annotation_legend_param = list(
    Cluster = list(
      title = "Source cluster",
      labels =
        cluster_legend_labels,
      at =
        cluster_factor_levels,
      ncol = 1
    ),
    Direction = list(
      title = "Donor-group effect",
      at = c(
        "UP",
        "DOWN"
      ),
      labels = c(
        "UP in ASD",
        "DOWN in ASD"
      )
    )
  )
)

fdr_count_005 <- sum(
  significant_results$FDR < 0.05,
  na.rm = TRUE
)

fdr_count_001 <- sum(
  significant_results$FDR < 0.01,
  na.rm = TRUE
)

fdr_count_0001 <- sum(
  significant_results$FDR < 0.001,
  na.rm = TRUE
)

fdr_marker_legend <- ComplexHeatmap::Legend(
  title = "FDR marker",
  labels = c(
    paste0(
      "FDR < 0.05 (n = ",
      fdr_count_005,
      ")"
    ),
    paste0(
      "FDR < 0.01 (n = ",
      fdr_count_001,
      ")"
    ),
    paste0(
      "FDR < 0.001 (n = ",
      fdr_count_0001,
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

# ==============================================================================
# 16. Heatmap title and dimensions
# ==============================================================================

number_of_clusters <- length(
  unique(
    significant_results$cluster_id
  )
)

number_of_samples <- ncol(
  zscore_matrix
)

number_up <- sum(
  significant_results$logFC > 0,
  na.rm = TRUE
)

number_down <- sum(
  significant_results$logFC < 0,
  na.rm = TRUE
)

heatmap_title <- paste(
  c(
    "FMT donor-group main effect after Crescendo batch correction",
    paste0(
      "Test: ",
      selected_test_comparison,
      " | test_id = ",
      donor_group_test_id
    ),
    paste0(
      "Selection: FDR < ",
      fdr_threshold,
      " | no log2FC threshold",
      " | ",
      number_of_results,
      " cluster-gene results",
      " | ",
      number_of_unique_genes,
      " unique genes",
      " | UP: ",
      number_up,
      " | DOWN: ",
      number_down
    ),
    paste0(
      "Expression input: Crescendo-corrected pseudobulk counts",
      " -> log2(count + 1)",
      " -> quantile normalization within each cluster"
    ),
    paste0(
      "Heatmap transformation: row Z-score clipped to [",
      zscore_lower_limit,
      ", ",
      zscore_upper_limit,
      "]"
    ),
    paste0(
      "Columns: Male Neurotypical -> Male ASD -> ",
      "Female Neurotypical -> Female ASD",
      " | ",
      number_of_samples,
      " samples"
    ),
    paste0(
      "Rows: significant cluster-gene results",
      " | ",
      number_of_clusters,
      " anatomical clusters represented"
    ),
    paste0(
      "Row clustering: Pearson correlation",
      " | base distance = ",
      row_clustering_base_distance,
      " | same-cluster bonus = ",
      row_clustering_same_cluster_bonus,
      " | linkage = ",
      row_clustering_linkage_method
    ),
    clustering_title
  ),
  collapse = "\n"
)

row_height_mm <- if (
  number_of_results <= 100L
) {
  3.1
} else if (
  number_of_results <= 300L
) {
  2.35
} else if (
  number_of_results <= 700L
) {
  1.90
} else {
  1.60
}

heatmap_body_height_mm <-
  number_of_results *
  row_height_mm

heatmap_body_width_mm <-
  number_of_samples *
  6.2

pdf_width_inches <- 14

pdf_height_inches <- max(
  14,
  (
    heatmap_body_height_mm +
      110
  ) /
    25.4
)

message(
  "\nPDF dimensions: ",
  round(
    pdf_width_inches,
    2
  ),
  " x ",
  round(
    pdf_height_inches,
    2
  ),
  " inches"
)

# ==============================================================================
# 17. Create heatmap
# ==============================================================================

heatmap_object <- ComplexHeatmap::Heatmap(
  zscore_matrix,

  name = "Row Z-score",

  col =
    heatmap_color_function,

  top_annotation =
    top_annotation,

  right_annotation =
    right_annotation,

  cluster_columns = FALSE,

  column_split =
    sample_metadata$combined_group,

  cluster_column_slices = FALSE,

  column_gap =
    grid::unit(
      1.2,
      "mm"
    ),

  cluster_rows =
    cluster_rows_current,

  row_dend_reorder =
    row_dend_reorder_current,

  row_gap =
    grid::unit(
      0,
      "mm"
    ),

  row_dend_width =
    grid::unit(
      18,
      "mm"
    ),

  show_column_names = TRUE,

  column_names_side = "top",

  column_names_rot = 45,

  column_names_centered = TRUE,

  column_names_gp =
    grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),

  column_names_max_height =
    grid::unit(
      18,
      "mm"
    ),

  show_row_names = FALSE,

  na_col = "grey90",

  rect_gp =
    grid::gpar(
      col = NA
    ),

  border = FALSE,

  use_raster = FALSE,

  width =
    grid::unit(
      heatmap_body_width_mm,
      "mm"
    ),

  height =
    grid::unit(
      heatmap_body_height_mm,
      "mm"
    ),

  column_title =
    heatmap_title,

  column_title_side = "top",

  column_title_gp =
    grid::gpar(
      fontsize = 12,
      fontface = "bold",
      lineheight = 1.12
    ),

  heatmap_legend_param = list(
    title = "Row Z-score",
    at = c(
      -2,
      -1,
      0,
      1,
      2
    ),
    legend_height =
      grid::unit(
        46,
        "mm"
      )
  )
)

# ==============================================================================
# 18. Write PDF
# ==============================================================================

message("\nWriting heatmap PDF:")
message(heatmap_pdf_file)

grDevices::pdf(
  file = heatmap_pdf_file,
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
  padding =
    grid::unit(
      c(
        3,
        3,
        3,
        3
      ),
      "mm"
    )
)

grDevices::dev.off()

if (
  !file.exists(heatmap_pdf_file) ||
    file.info(heatmap_pdf_file)$size == 0
) {
  stop(
    "Heatmap PDF was not created correctly: ",
    heatmap_pdf_file,
    call. = FALSE
  )
}

# ==============================================================================
# 19. Write displayed results TSV
# ==============================================================================

displayed_results <- significant_results |>
  dplyr::select(
    .data$cluster_id,
    .data$ensembl_gene_id,
    .data$gene,
    .data$logFC,
    .data$PValue,
    .data$FDR,
    .data$effect_direction,
    .data$n_clusters_for_gene,
    .data$row_id
  )

readr::write_tsv(
  displayed_results,
  displayed_results_tsv,
  na = "NA"
)

# ==============================================================================
# 20. Final report
# ==============================================================================

message("\n============================================================")
message("CRESCENDO-CORRECTED DONOR-GROUP HEATMAP FINISHED")
message("============================================================")

message(
  "Cluster-gene results displayed: ",
  number_of_results
)

message(
  "Unique genes displayed: ",
  number_of_unique_genes
)

message(
  "UP cluster-gene results: ",
  number_up
)

message(
  "DOWN cluster-gene results: ",
  number_down
)

message("\nPDF:")
message(
  normalizePath(
    heatmap_pdf_file,
    mustWork = TRUE
  )
)

message("\nDisplayed results TSV:")
message(
  normalizePath(
    displayed_results_tsv,
    mustWork = TRUE
  )
)

message("\nExpression source:")
message(
  normalizePath(
    log2_qn_input_dir,
    mustWork = TRUE
  )
)

# ==============================================================================
# End
# ==============================================================================


# ============================================================
# 21. Sex-specific expression differences: ASD vs NT
# ============================================================

# Use log2+QN expression BEFORE row-wise Z-scoring.

sex_group <- as.character(sample_metadata$sex)
donor_group <- as.character(sample_metadata$fmt_donor_group)

get_group_mean <- function(sex_value, group_value) {
  selected <- sex_group == sex_value &
    donor_group == group_value

  if (sum(selected) < 2L) {
    stop("Too few samples: ", sex_value, " ", group_value)
  }

  rowMeans(
    sample_level_log2_qn_matrix[, selected, drop = FALSE]
  )
}

# ASD - Neurotypical, separately for each sex

delta_male <-
  get_group_mean("Male", "ASD") -
  get_group_mean("Male", "Neurotypical")

delta_female <-
  get_group_mean("Female", "ASD") -
  get_group_mean("Female", "Neurotypical")

# One row per significant gene-cluster result

sex_effects <- significant_results |>
  dplyr::mutate(
    delta_male = delta_male,
    delta_female = delta_female,
    abs_delta_male = abs(delta_male),
    abs_delta_female = abs(delta_female)
  )

# Long format

sex_effects_long <- dplyr::bind_rows(
  sex_effects |>
    dplyr::transmute(
      row_id,
      direction = effect_direction,
      sex = "Male",
      delta = delta_male,
      abs_delta = abs_delta_male
    ),

  sex_effects |>
    dplyr::transmute(
      row_id,
      direction = effect_direction,
      sex = "Female",
      delta = delta_female,
      abs_delta = abs_delta_female
    )
)

# UP / DOWN / ALL

sex_effects_all <- dplyr::bind_rows(
  sex_effects_long |>
    dplyr::mutate(category = direction),

  sex_effects_long |>
    dplyr::mutate(category = "ALL")
) |>
  dplyr::filter(
    category %in% c("UP", "DOWN", "ALL")
  )

# Mean +/- SD across gene-cluster results

summary_table <- sex_effects_all |>
  dplyr::group_by(category, sex) |>
  dplyr::summarise(
    n = dplyr::n(),

    mean_delta = mean(delta),
    sd_delta = sd(delta),

    mean_abs_delta = mean(abs_delta),
    sd_abs_delta = sd(abs_delta),

    .groups = "drop"
  ) |>
  dplyr::mutate(
    category = factor(
      category,
      levels = c("UP", "DOWN", "ALL")
    ),
    sex = factor(
      sex,
      levels = c("Male", "Female")
    ),
    signed_mean_SD = sprintf(
      "%.3f +/- %.3f",
      mean_delta,
      sd_delta
    ),
    absolute_mean_SD = sprintf(
      "%.3f +/- %.3f",
      mean_abs_delta,
      sd_abs_delta
    )
  ) |>
  dplyr::arrange(category, sex)

cat("\n============================================\n")
cat("SEX-SPECIFIC ASD VS NEUROTYPICAL DIFFERENCES\n")
cat("============================================\n\n")

print(
  summary_table |>
    dplyr::select(
      category,
      sex,
      n,
      signed_mean_SD,
      absolute_mean_SD
    ),
  n = Inf
)

# Save results

readr::write_tsv(
  summary_table,
  file.path(
    heatmap_output_dir,
    "03_sexSpecific_meanSD_summary.tsv"
  )
)

readr::write_tsv(
  sex_effects,
  file.path(
    heatmap_output_dir,
    "04_sexSpecific_geneCluster_differences.tsv"
  )
)
