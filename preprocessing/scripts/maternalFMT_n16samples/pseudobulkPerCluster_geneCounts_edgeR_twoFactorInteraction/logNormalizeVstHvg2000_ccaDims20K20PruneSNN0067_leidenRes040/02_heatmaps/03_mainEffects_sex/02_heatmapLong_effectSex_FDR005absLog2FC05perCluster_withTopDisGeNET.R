#!/usr/bin/env Rscript

# ==============================================================================
# 01_heatmapLong_effectSex_FemaleVsMale_FDR005absLog2FC05perCluster_withDisgenetAnxiety.R
#
# Long heatmap for the main sex effect from the cluster-specific pseudobulk
# edgeR two-factor analysis, with an additional binary DisGeNET Anxiety
# overlap column placed between the expression heatmap and row annotations.
#
# TEST
#   Overall_Sex_Female_vs_Male
#
# INTERPRETATION
#   Positive logFC = higher expression in Female
#   Negative logFC = higher expression in Male
#
# SELECTION
#   FDR < 0.05
#   abs(log2FC) >= 0.5
#
# ROWS
#   Significant cluster-gene results.
#   The same gene can occur more than once if significant in several clusters.
#
# COLUMNS
#   Samples ordered as:
#     Male Neurotypical -> Male ASD -> Female Neurotypical -> Female ASD
#
# HEATMAP VALUES
#   TMM-normalized CPM per sample and cluster
#   -> log2(CPM + 1)
#   -> row-wise Z-score
#   -> clipped to [-2, 2]
#
# ROW CLUSTERING
#   One global hierarchical clustering.
#   Base distance:
#     2 - Pearson correlation
#
#   Rows derived from the same anatomical cluster are made slightly closer,
#   while rows derived from different clusters are made slightly more distant.
#
#
# DISGENET OVERLAP HEATMAP
#   One DisGeNET term is shown:
#     Anxiety
#
#   Binary coding:
#     1 = gene belongs to the DisGeNET Anxiety gene set
#     0 = gene does not belong to the DisGeNET Anxiety gene set
#
#   Colours:
#     1 = green
#     0 = light grey
#
#   Mouse/DisGeNET symbols are matched case-insensitively by gene symbol.
#   No orthology conversion is performed.
#
# RIGHT ANNOTATIONS
#   Cluster
#     Source cluster.
#
#   Max mean %
#     Maximum mean percentage of positive spots across:
#       Male
#       Female
#
#     For the sex main effect, these percentages are aggregated across
#     donor groups in the edgeR result table.
#
#   FDR marker
#       *   FDR < 0.05
#       **  FDR < 0.01
#       *** FDR < 0.001
#
#   |log2FC| marker
#       #   abs(log2FC) >= 0.5
#       ##  abs(log2FC) >= 0.8
#       ### abs(log2FC) >= 1.0
#
#   n clusters
#     Number of significant source clusters in which the same gene occurs
#     on this heatmap.
#
#   Gene symbol
#     Gene symbol with source-cluster tag, e.g. G3bp1 [C16].
#
# LEGEND COUNTS
#   Cluster legend:
#     n, UP and DOWN counts are numbers of displayed cluster-gene results.
#
#   FDR legend:
#     For every cumulative FDR threshold, the legend reports:
#       n = number of displayed cluster-gene results satisfying the threshold
#       UP = number with logFC > 0
#       DOWN = number with logFC < 0
#
#   |log2FC| legend:
#     For every cumulative absolute log2FC threshold, the legend reports:
#       n = number of displayed cluster-gene results satisfying the threshold
#       UP = number with logFC > 0
#       DOWN = number with logFC < 0
#
# OUTPUT
#   1. One long PDF heatmap.
#   2. TSV with displayed results per cluster.
#   3. TSV with FDR and |log2FC| threshold counts.
# ==============================================================================

options(stringsAsFactors = FALSE)


# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "GSA",
  "dplyr",
  "tibble",
  "tidyr",
  "readr",
  "edgeR",
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
  library(GSA)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(edgeR)
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

sex_test_id <- "Overall_Sex_Female_vs_Male"

fdr_threshold <- 0.05
abs_log2fc_threshold <- 0.5

zscore_lower_limit <- -2
zscore_upper_limit <- 2


parameter_block <- "FDR005absLog2FC05perCluster"

disgenet_library_url <-
  "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=DisGeNET"

disgenet_diseases <- c(
  "Anxiety"
)

disgenet_absent_colour <- "#F2F2F2"
disgenet_present_colour <- "#2EAD55"

# Missing sample-cluster CPM values are retained as NA.
# No imputation is performed.
impute_missing <- FALSE

# Row-clustering parameters.
row_clustering_base_distance <- "2 - r"
row_clustering_r_scale <- 1
row_clustering_same_cluster_bonus <- 0.15
row_clustering_linkage_method <- "complete"

# Exact clustering information shown in the technical heatmap title.
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

edgeR_model_objects_rdata_file <- file.path(
  statistics_dir,
  "05_objects",
  paste0(
    analysis_prefix,
    "_edgeRModelObjects.RData"
  )
)

heatmap_output_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "02_heatmaps",
  "03_mainEffects_sex",
  "01_long"
)

dir.create(
  heatmap_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
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


# ==============================================================================
# 5. Output files
# ==============================================================================

heatmap_pdf_file <- file.path(
  heatmap_output_dir,
  paste0(
    "01_",
    dataset_name,
    "_effectSex_FemaleVsMale_",
    parameter_block,
    "_withDisgenetAnxiety_heatmapLong.pdf"
  )
)

results_per_cluster_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "02_",
    dataset_name,
    "_effectSex_FemaleVsMale_",
    parameter_block,
    "_resultsPerCluster.tsv"
  )
)

threshold_counts_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "03_",
    dataset_name,
    "_effectSex_FemaleVsMale_",
    parameter_block,
    "_thresholdCounts.tsv"
  )
)


# ==============================================================================
# 6. Helper functions
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
  same_cluster_bonus = 0.15
) {

  if (!is.matrix(expression_matrix)) {
    expression_matrix <- as.matrix(expression_matrix)
  }

  if (nrow(expression_matrix) <= 1L) {
    return(
      stats::as.dist(
        matrix(
          0,
          nrow = 1,
          ncol = 1
        )
      )
    )
  }

  cluster_labels <- as.character(cluster_labels)

  if (length(cluster_labels) != nrow(expression_matrix)) {
    stop(
      "The number of cluster labels must equal the number of heatmap rows.",
      call. = FALSE
    )
  }

  if (!is.numeric(r_scale) || length(r_scale) != 1L || !is.finite(r_scale)) {
    stop(
      "`r_scale` must be one finite numeric value.",
      call. = FALSE
    )
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
  same_cluster_bonus = 0.15,
  linkage_method = "complete"
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


count_direction_for_subset <- function(data_subset) {

  tibble::tibble(
    n_results = nrow(data_subset),
    UP = sum(data_subset$logFC > 0, na.rm = TRUE),
    DOWN = sum(data_subset$logFC < 0, na.rm = TRUE)
  )
}


format_threshold_count_label <- function(
  threshold_label,
  count_row
) {

  paste0(
    threshold_label,
    " (n = ",
    count_row$n_results[[1]],
    ", UP: ",
    count_row$UP[[1]],
    ", DOWN: ",
    count_row$DOWN[[1]],
    ")"
  )
}


# ==============================================================================
# 7. Load edgeR objects
# ==============================================================================

message("Loading edgeR results RData:")
message(edgeR_results_rdata_file)
load(edgeR_results_rdata_file)

message("\nLoading edgeR model-object RData:")
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
    "Required object(s) missing after loading RData files: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}


# ==============================================================================
# 8. Select main sex-effect test
# ==============================================================================

if (!sex_test_id %in% names(edgeR_perCluster_combinedResults)) {

  message("\nAvailable tests:")
  message(
    paste(
      names(edgeR_perCluster_combinedResults),
      collapse = "\n"
    )
  )

  stop(
    "Sex main-effect test not found: ",
    sex_test_id,
    call. = FALSE
  )
}

sex_test_definition <- edgeR_perCluster_testDefinitions |>
  dplyr::filter(
    .data$test_id == sex_test_id
  )

if (nrow(sex_test_definition) != 1L) {
  stop(
    "Sex test_id does not map to exactly one test definition: ",
    sex_test_id,
    call. = FALSE
  )
}

selected_test_comparison <- as.character(
  sex_test_definition$comparison[[1]]
)

selected_test_sheet <- as.character(
  sex_test_definition$sheet_name[[1]]
)

message("\nSelected sex test_id: ", sex_test_id)
message("Selected comparison: ", selected_test_comparison)
message("Selected sheet: ", selected_test_sheet)


# ==============================================================================
# 9. Select significant sex-effect results
# ==============================================================================

sex_results <- edgeR_perCluster_combinedResults[[sex_test_id]] |>
  tibble::as_tibble()

required_result_columns_base <- c(
  "cluster_id",
  "ensembl_gene_id",
  "gene",
  "logFC",
  "PValue",
  "FDR"
)

missing_result_columns_base <- setdiff(
  required_result_columns_base,
  colnames(sex_results)
)

if (length(missing_result_columns_base) > 0L) {
  stop(
    "Sex-effect results are missing required columns: ",
    paste(missing_result_columns_base, collapse = ", "),
    call. = FALSE
  )
}

# Exact expected columns for the main sex effect.
male_percent_column <- "mean_percent_positive_spots_Male"
female_percent_column <- "mean_percent_positive_spots_Female"

if (
  !male_percent_column %in% colnames(sex_results) ||
  !female_percent_column %in% colnames(sex_results)
) {

  percent_related_columns <- grep(
    "percent|positive|spots|male|female",
    colnames(sex_results),
    ignore.case = TRUE,
    value = TRUE
  )

  message(
    "\nExpected sex-specific percent-positive columns were not found."
  )

  message(
    "Expected:"
  )

  message(
    paste(
      c(
        male_percent_column,
        female_percent_column
      ),
      collapse = "\n"
    )
  )

  message(
    "\nColumns potentially related to percent-positive spots:"
  )

  if (length(percent_related_columns) > 0L) {
    message(
      paste(
        percent_related_columns,
        collapse = "\n"
      )
    )
  } else {
    message("<none>")
  }

  stop(
    "Cannot calculate Max mean % for the sex-effect heatmap.",
    call. = FALSE
  )
}

significant_sex <- sex_results |>
  dplyr::filter(
    !is.na(.data$FDR),
    !is.na(.data$logFC),
    .data$FDR < fdr_threshold,
    abs(.data$logFC) >= abs_log2fc_threshold
  ) |>
  dplyr::mutate(
    cluster_id = as.character(.data$cluster_id),

    abs_logFC = abs(.data$logFC),

    effect_direction = dplyr::case_when(
      .data$logFC > 0 ~ "UP",
      .data$logFC < 0 ~ "DOWN",
      TRUE ~ "ZERO"
    ),

    max_mean_percent = pmax(
      .data[[male_percent_column]],
      .data[[female_percent_column]],
      na.rm = TRUE
    ),

    max_mean_percent = dplyr::if_else(
      is.infinite(.data$max_mean_percent),
      NA_real_,
      .data$max_mean_percent
    ),

    fdr_symbol = dplyr::case_when(
      .data$FDR < 0.001 ~ "***",
      .data$FDR < 0.01 ~ "**",
      .data$FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),

    log2fc_symbol = dplyr::case_when(
      .data$abs_logFC >= 1.0 ~ "###",
      .data$abs_logFC >= 0.8 ~ "##",
      .data$abs_logFC >= 0.5 ~ "#",
      TRUE ~ ""
    )
  ) |>
  dplyr::group_by(.data$ensembl_gene_id) |>
  dplyr::mutate(
    n_clusters_for_gene = dplyr::n_distinct(.data$cluster_id)
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(
    as.integer(.data$cluster_id),
    .data$FDR,
    dplyr::desc(.data$abs_logFC)
  )

if (nrow(significant_sex) == 0L) {
  stop(
    "No sex-effect result satisfies FDR < ",
    fdr_threshold,
    " and abs(log2FC) >= ",
    abs_log2fc_threshold,
    ".",
    call. = FALSE
  )
}

number_of_results <- nrow(significant_sex)

number_of_unique_genes <- length(
  unique(significant_sex$ensembl_gene_id)
)

significant_clusters <- sort(
  unique(
    as.integer(
      significant_sex$cluster_id
    )
  )
)

significant_cluster_ids <- as.character(
  significant_clusters
)

message(
  "\nSignificant cluster-gene results: ",
  number_of_results
)

message(
  "Unique significant genes: ",
  number_of_unique_genes
)

message(
  "Significant clusters: ",
  paste(significant_cluster_ids, collapse = ", ")
)


# ==============================================================================
# 10. Calculate threshold counts for legends
# ==============================================================================

fdr_count_005 <- count_direction_for_subset(
  significant_sex |>
    dplyr::filter(.data$FDR < 0.05)
)

fdr_count_001 <- count_direction_for_subset(
  significant_sex |>
    dplyr::filter(.data$FDR < 0.01)
)

fdr_count_0001 <- count_direction_for_subset(
  significant_sex |>
    dplyr::filter(.data$FDR < 0.001)
)

log2fc_count_05 <- count_direction_for_subset(
  significant_sex |>
    dplyr::filter(.data$abs_logFC >= 0.5)
)

log2fc_count_08 <- count_direction_for_subset(
  significant_sex |>
    dplyr::filter(.data$abs_logFC >= 0.8)
)

log2fc_count_10 <- count_direction_for_subset(
  significant_sex |>
    dplyr::filter(.data$abs_logFC >= 1.0)
)

fdr_legend_labels <- c(
  format_threshold_count_label(
    "FDR < 0.05",
    fdr_count_005
  ),
  format_threshold_count_label(
    "FDR < 0.01",
    fdr_count_001
  ),
  format_threshold_count_label(
    "FDR < 0.001",
    fdr_count_0001
  )
)

log2fc_legend_labels <- c(
  format_threshold_count_label(
    "|log2FC| >= 0.5",
    log2fc_count_05
  ),
  format_threshold_count_label(
    "|log2FC| >= 0.8",
    log2fc_count_08
  ),
  format_threshold_count_label(
    "|log2FC| >= 1.0",
    log2fc_count_10
  )
)

threshold_counts <- dplyr::bind_rows(
  tibble::tibble(
    statistic = "FDR",
    threshold = "< 0.05",
    fdr_count_005
  ),
  tibble::tibble(
    statistic = "FDR",
    threshold = "< 0.01",
    fdr_count_001
  ),
  tibble::tibble(
    statistic = "FDR",
    threshold = "< 0.001",
    fdr_count_0001
  ),
  tibble::tibble(
    statistic = "|log2FC|",
    threshold = ">= 0.5",
    log2fc_count_05
  ),
  tibble::tibble(
    statistic = "|log2FC|",
    threshold = ">= 0.8",
    log2fc_count_08
  ),
  tibble::tibble(
    statistic = "|log2FC|",
    threshold = ">= 1.0",
    log2fc_count_10
  )
)

readr::write_tsv(
  threshold_counts,
  threshold_counts_tsv,
  na = "NA"
)

message("\nThreshold counts used in legends:")

print(
  threshold_counts,
  n = Inf
)


# ==============================================================================
# 11. Check availability of cluster model objects
# ==============================================================================

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


# ==============================================================================
# 12. Extract sample metadata from cluster DGEList objects
# ==============================================================================

sample_metadata_list <- lapply(
  significant_cluster_ids,
  function(cluster_id_current) {

    model_object <- edgeR_perCluster_modelObjects[[cluster_id_current]]

    if (is.null(model_object$dge)) {
      stop(
        "DGEList is missing for cluster ",
        cluster_id_current,
        call. = FALSE
      )
    }

    dge_samples <- model_object$dge$samples |>
      tibble::rownames_to_column("dge_rowname") |>
      tibble::as_tibble()

    required_dge_sample_columns <- c(
      "sample_ID",
      "fmt_donor_group",
      "sex"
    )

    missing_dge_sample_columns <- setdiff(
      required_dge_sample_columns,
      colnames(dge_samples)
    )

    if (length(missing_dge_sample_columns) > 0L) {
      stop(
        "DGEList sample metadata for cluster ",
        cluster_id_current,
        " is missing: ",
        paste(missing_dge_sample_columns, collapse = ", "),
        call. = FALSE
      )
    }

    dge_samples |>
      dplyr::transmute(
        sample_ID = as.character(.data$sample_ID),
        fmt_donor_group = as.character(.data$fmt_donor_group),
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
    name = "n_metadata_rows"
  ) |>
  dplyr::filter(
    .data$n_metadata_rows > 1L
  )

if (nrow(metadata_conflicts) > 0L) {
  stop(
    "Conflicting Sex/Group metadata found for sample(s): ",
    paste(metadata_conflicts$sample_ID, collapse = ", "),
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
  )

if (anyNA(sample_metadata$sex)) {
  stop(
    "Unexpected sex value in DGEList metadata.",
    call. = FALSE
  )
}

if (anyNA(sample_metadata$fmt_donor_group)) {
  stop(
    "Unexpected fmt_donor_group value in DGEList metadata.",
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
# 13. Calculate sample-level TMM-normalized CPM
# ==============================================================================

cluster_cpm_matrices <- lapply(
  significant_cluster_ids,
  function(cluster_id_current) {

    dge <- edgeR_perCluster_modelObjects[[cluster_id_current]]$dge

    cpm_matrix <- edgeR::cpm(
      dge,
      normalized.lib.sizes = TRUE,
      log = FALSE
    )

    if (is.null(rownames(cpm_matrix))) {
      stop(
        "TMM-CPM matrix has no gene IDs for cluster ",
        cluster_id_current,
        call. = FALSE
      )
    }

    if (is.null(colnames(cpm_matrix))) {
      stop(
        "TMM-CPM matrix has no sample IDs for cluster ",
        cluster_id_current,
        call. = FALSE
      )
    }

    cpm_matrix
  }
)

names(cluster_cpm_matrices) <- significant_cluster_ids


# ==============================================================================
# 14. Build cluster-gene x sample CPM matrix
# ==============================================================================

row_ids <- paste0(
  ifelse(
    is.na(significant_sex$gene) |
      significant_sex$gene == "",
    significant_sex$ensembl_gene_id,
    significant_sex$gene
  ),
  " [C",
  significant_sex$cluster_id,
  "]"
)

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

significant_sex <- significant_sex |>
  dplyr::mutate(
    row_id = row_ids
  )


# ==============================================================================
# 14A. Download DisGeNET library and build binary overlap matrix
# ==============================================================================

message("\nDownloading DisGeNET gene-set library:")
message(disgenet_library_url)

disgenet <- GSA::GSA.read.gmt(
  url(disgenet_library_url)
)

disgenet_lists <- setNames(
  disgenet$genesets,
  disgenet$geneset.names
)

wanted_idx <- match(
  tolower(disgenet_diseases),
  tolower(names(disgenet_lists))
)

disease_summary <- tibble::tibble(
  requested = disgenet_diseases,

  matched = ifelse(
    is.na(wanted_idx),
    NA_character_,
    names(disgenet_lists)[wanted_idx]
  ),

  n_genes_in_disgenet = vapply(
    wanted_idx,
    function(i) {
      if (is.na(i)) {
        return(NA_integer_)
      }

      length(
        unique(
          disgenet_lists[[i]]
        )
      )
    },
    FUN.VALUE = integer(1)
  )
)

message("\nDisGeNET term matching:")

print(
  disease_summary,
  n = Inf
)

if (anyNA(wanted_idx)) {

  missing_disgenet_terms <- disgenet_diseases[
    is.na(wanted_idx)
  ]

  stop(
    "The following requested DisGeNET terms were not found: ",
    paste(
      missing_disgenet_terms,
      collapse = ", "
    ),
    call. = FALSE
  )
}

selected_disgenet <- disgenet_lists[
  wanted_idx
]

names(selected_disgenet) <- disgenet_diseases

selected_disgenet_upper <- lapply(
  selected_disgenet,
  function(gene_vector) {
    unique(
      toupper(
        trimws(
          as.character(
            gene_vector
          )
        )
      )
    )
  }
)

plot_gene_symbols <- ifelse(
  is.na(significant_sex$gene) |
    significant_sex$gene == "",
  NA_character_,
  as.character(
    significant_sex$gene
  )
)

plot_gene_symbols_upper <- toupper(
  trimws(
    plot_gene_symbols
  )
)

disgenet_binary_matrix <- vapply(
  selected_disgenet_upper,
  function(disgenet_gene_set) {

    as.integer(
      !is.na(plot_gene_symbols_upper) &
        plot_gene_symbols_upper %in% disgenet_gene_set
    )
  },
  FUN.VALUE = integer(
    nrow(
      significant_sex
    )
  )
)

disgenet_binary_matrix <- as.matrix(
  disgenet_binary_matrix
)

rownames(
  disgenet_binary_matrix
) <- row_ids

colnames(
  disgenet_binary_matrix
) <- disgenet_diseases

storage.mode(
  disgenet_binary_matrix
) <- "numeric"

if (
  !identical(
    rownames(disgenet_binary_matrix),
    row_ids
  )
) {
  stop(
    "DisGeNET matrix row order does not match expression heatmap rows.",
    call. = FALSE
  )
}

disgenet_overlap_counts <- tibble::tibble(
  disease = disgenet_diseases,
  n_cluster_gene_rows = as.integer(
    colSums(
      disgenet_binary_matrix,
      na.rm = TRUE
    )
  )
)

message("\nDisGeNET overlap counts among displayed cluster-gene rows:")

print(
  disgenet_overlap_counts,
  n = Inf
)


# ==============================================================================
# 14B. Build cluster-gene x sample CPM matrix
# ==============================================================================

sample_level_cpm_matrix <- matrix(
  NA_real_,
  nrow = nrow(significant_sex),
  ncol = length(sample_order),
  dimnames = list(
    row_ids,
    sample_order
  )
)

for (row_index in seq_len(nrow(significant_sex))) {

  current_cluster <- significant_sex$cluster_id[[row_index]]
  current_gene_id <- significant_sex$ensembl_gene_id[[row_index]]
  current_cpm_matrix <- cluster_cpm_matrices[[current_cluster]]

  gene_match <- match(
    current_gene_id,
    rownames(current_cpm_matrix)
  )

  if (is.na(gene_match)) {
    stop(
      "Gene ",
      current_gene_id,
      " is missing from the tested DGEList genes for cluster ",
      current_cluster,
      call. = FALSE
    )
  }

  available_samples <- intersect(
    sample_order,
    colnames(current_cpm_matrix)
  )

  sample_level_cpm_matrix[
    row_index,
    available_samples
  ] <- current_cpm_matrix[
    gene_match,
    available_samples,
    drop = TRUE
  ]
}

if (all(is.na(sample_level_cpm_matrix))) {
  stop(
    "The complete cluster-gene CPM matrix contains only NA values.",
    call. = FALSE
  )
}

number_missing_values <- sum(
  is.na(sample_level_cpm_matrix)
)

message(
  "\nMissing sample-cluster CPM values retained as NA: ",
  number_missing_values
)


# ==============================================================================
# 15. Transform CPM -> log2(CPM + 1) -> row Z-score
# ==============================================================================

log2_cpm_matrix <- log2(
  sample_level_cpm_matrix + 1
)

zscore_matrix <- calculate_row_zscores(
  expression_matrix = log2_cpm_matrix,
  lower_limit = zscore_lower_limit,
  upper_limit = zscore_upper_limit
)


# ==============================================================================
# 16. Build row annotation table
# ==============================================================================

row_annotation_table <- significant_sex |>
  dplyr::transmute(
    row_id = .data$row_id,
    cluster_id = .data$cluster_id,

    cluster_label = vapply(
      .data$cluster_id,
      format_cluster_label,
      FUN.VALUE = character(1)
    ),

    logFC = .data$logFC,
    abs_logFC = .data$abs_logFC,
    effect_direction = .data$effect_direction,
    PValue = .data$PValue,
    FDR = .data$FDR,
    max_mean_percent = .data$max_mean_percent,
    fdr_symbol = .data$fdr_symbol,
    log2fc_symbol = .data$log2fc_symbol,
    n_clusters_for_gene = .data$n_clusters_for_gene
  )


# ==============================================================================
# 17. Cluster heatmap rows
# ==============================================================================

if (nrow(zscore_matrix) > 1L) {

  hc_rows <- create_modified_row_hclust(
    expression_matrix = zscore_matrix,
    annotation_table = row_annotation_table,
    base_distance = row_clustering_base_distance,
    r_scale = row_clustering_r_scale,
    same_cluster_bonus = row_clustering_same_cluster_bonus,
    linkage_method = row_clustering_linkage_method
  )

  cluster_rows <- hc_rows
  row_dend_reorder <- row_annotation_table$logFC

} else {

  cluster_rows <- FALSE
  row_dend_reorder <- FALSE
}


# ==============================================================================
# 18. Common colours
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

sex_colors <- c(
  "Male" = "#3C78D8",
  "Female" = "#C050A0"
)

group_colors <- c(
  "Neurotypical" = "grey70",
  "ASD" = "black"
)

disgenet_color_function <- circlize::colorRamp2(
  c(
    0,
    1
  ),
  c(
    disgenet_absent_colour,
    disgenet_present_colour
  )
)


# ==============================================================================
# 19. Cluster legend
# ==============================================================================

cluster_factor_levels <- as.character(
  sort(
    unique(
      as.integer(
        row_annotation_table$cluster_id
      )
    )
  )
)

cluster_annotation_colors <- custom_cluster_colors[
  cluster_factor_levels
]

names(cluster_annotation_colors) <- cluster_factor_levels

cluster_result_counts <- table(
  factor(
    row_annotation_table$cluster_id,
    levels = cluster_factor_levels
  )
)

cluster_direction_counts <- row_annotation_table |>
  dplyr::mutate(
    direction = dplyr::case_when(
      .data$logFC > 0 ~ "UP",
      .data$logFC < 0 ~ "DOWN",
      TRUE ~ "ZERO"
    )
  ) |>
  dplyr::filter(
    .data$direction %in% c(
      "UP",
      "DOWN"
    )
  ) |>
  dplyr::count(
    .data$cluster_id,
    .data$direction,
    name = "n_direction"
  ) |>
  tidyr::pivot_wider(
    names_from = .data$direction,
    values_from = .data$n_direction,
    values_fill = 0L
  ) |>
  dplyr::mutate(
    cluster_id = as.character(.data$cluster_id)
  )

cluster_direction_counts_full <- tibble::tibble(
  cluster_id = cluster_factor_levels
) |>
  dplyr::left_join(
    cluster_direction_counts,
    by = "cluster_id"
  ) |>
  dplyr::mutate(
    UP = dplyr::coalesce(.data$UP, 0L),
    DOWN = dplyr::coalesce(.data$DOWN, 0L)
  )

cluster_legend_labels <- paste0(
  "C",
  cluster_factor_levels,
  " (n = ",
  as.integer(cluster_result_counts),
  ", UP: ",
  cluster_direction_counts_full$UP,
  ", DOWN: ",
  cluster_direction_counts_full$DOWN,
  "): ",
  unname(
    custom_cluster_labels[
      cluster_factor_levels
    ]
  )
)


# ==============================================================================
# 20. Top annotation
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
# 21. Annotation font size
# ==============================================================================

annotation_row_font_size <- if (number_of_results <= 100L) {
  7.0
} else if (number_of_results <= 300L) {
  5.6
} else if (number_of_results <= 700L) {
  4.4
} else {
  3.6
}


# ==============================================================================
# 22. Max mean % row annotation
# ==============================================================================

percent_annotation_labels <- ifelse(
  is.na(row_annotation_table$max_mean_percent),
  "NA",
  paste0(
    round(
      row_annotation_table$max_mean_percent
    ),
    "%"
  )
)

percent_annotation_text_colour <- ifelse(
  !is.na(row_annotation_table$max_mean_percent) &
    row_annotation_table$max_mean_percent > 70,
  "white",
  "black"
)

max_mean_percent_annotation <- ComplexHeatmap::AnnotationFunction(
  which = "row",

  # Slightly wider than the donor-group version so 3-digit percentages fit.
  width = grid::unit(
    7.0,
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

    n_current <- length(index)

    y_positions <- 1 - (
      seq_len(n_current) - 0.5
    ) / n_current

    current_values <- values[index]
    current_labels <- labels[index]
    current_text_colours <- text_colours[index]

    current_fill <- rep(
      "grey90",
      n_current
    )

    finite_values <- is.finite(
      current_values
    )

    current_fill[finite_values] <- colour_function(
      current_values[finite_values]
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


# ==============================================================================
# 23. FDR marker legend with threshold counts
# ==============================================================================

fdr_marker_legend <- ComplexHeatmap::Legend(
  title = "FDR marker",

  labels = fdr_legend_labels,

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
# 24. |log2FC| marker legend with threshold counts
# ==============================================================================

log2fc_marker_legend <- ComplexHeatmap::Legend(
  title = "|log2FC| marker",

  labels = log2fc_legend_labels,

  graphics = list(
    function(x, y, w, h) {
      grid::grid.text(
        "#",
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
        "##",
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
        "###",
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
    9,
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
# 25. Right-side annotations
# ==============================================================================

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

  `|log2FC| marker` = ComplexHeatmap::anno_text(
    row_annotation_table$log2fc_symbol,

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
      4.5,
      "mm"
    )
  ),

  `n clusters` = ComplexHeatmap::anno_text(
    as.character(
      row_annotation_table$n_clusters_for_gene
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
    `|log2FC| marker` = TRUE,
    `n clusters` = TRUE,
    `Gene symbol` = TRUE
  ),

  annotation_label = c(
    Cluster = "Cluster",
    `Max mean %` = "Max mean %",
    `FDR marker` = "FDR marker",
    `|log2FC| marker` = "|log2FC| marker",
    `n clusters` = "n clusters",
    `Gene symbol` = "Gene symbol"
  ),

  annotation_name_side = "top",

  annotation_name_rot = c(
    Cluster = 45,
    `Max mean %` = 45,
    `FDR marker` = 45,
    `|log2FC| marker` = 45,
    `n clusters` = 45,
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


# ==============================================================================
# 26. Technical heatmap title
# ==============================================================================

number_of_clusters <- length(
  cluster_factor_levels
)

number_of_samples <- ncol(
  zscore_matrix
)

heatmap_title_lines <- c(
  "Sex main effect: significant genes per cluster",

  paste0(
    "Test: ",
    selected_test_comparison,
    " | test_id = ",
    sex_test_id
  ),

  paste0(
    "Contrast direction: UP = Female > Male",
    " | DOWN = Male > Female"
  ),

  paste0(
    "Selection: FDR < ",
    fdr_threshold,
    " and |log2FC| >= ",
    abs_log2fc_threshold,
    " | ",
    number_of_results,
    " cluster-gene results | ",
    number_of_unique_genes,
    " unique genes | ",
    number_of_clusters,
    " clusters"
  ),

  "Input: sample-level TMM-normalized CPM from each cluster DGEList",

  paste0(
    "Transformation: log2(CPM + 1) -> row Z-score clipped to [",
    zscore_lower_limit,
    ", ",
    zscore_upper_limit,
    "]"
  ),

  paste0(
    "Missing sample-cluster CPM retained as NA (shown in grey)",
    " | missing values: ",
    number_missing_values
  ),

  paste0(
    "Columns: Male Neurotypical -> Male ASD -> Female Neurotypical -> Female ASD",
    " | ",
    number_of_samples,
    " samples"
  ),

  paste0(
    "Rows: expression heatmap -> DisGeNET Anxiety overlap column -> row annotations"
  ),

  paste0(
    "DisGeNET Anxiety: green = gene present in the gene set, light grey = absent"
  ),

  paste0(
    "Row annotations: Cluster = source cluster",
    " | Max mean % = maximum mean positive spots across Male and Female"
  ),

  paste0(
    "FDR marker and |log2FC| marker encode statistical significance and effect-size thresholds",
    " | n clusters = number of significant clusters containing the same gene"
  ),

  paste0(
    "Legend threshold counts are cumulative",
    " | n, UP and DOWN refer to displayed cluster-gene results"
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


# ==============================================================================
# 27. PDF dimensions
# ==============================================================================

row_height_mm <- if (number_of_results <= 100L) {
  3.1
} else if (number_of_results <= 300L) {
  2.35
} else if (number_of_results <= 700L) {
  1.90
} else {
  1.60
}

heatmap_body_height_mm <- number_of_results * row_height_mm
heatmap_body_width_mm <- number_of_samples * 6.2

# Slightly wider than the donor-group version because the FDR and log2FC
# legends now include n / UP / DOWN counts.
pdf_width_inches <- 15

pdf_height_inches <- max(
  14,
  (
    heatmap_body_height_mm +
      115
  ) /
    25.4
)

message(
  "\nHeatmap rows: ",
  number_of_results
)

message(
  "Unique genes: ",
  number_of_unique_genes
)

message(
  "PDF dimensions: ",
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
# 28. Create expression heatmap + DisGeNET overlap heatmap + row annotations
# ==============================================================================

expression_heatmap <- ComplexHeatmap::Heatmap(
  zscore_matrix,

  name = "Row Z-score",
  col = heatmap_color_function,

  top_annotation = top_annotation,

  cluster_columns = FALSE,

  column_split = sample_metadata$combined_group,
  cluster_column_slices = FALSE,

  column_gap = grid::unit(
    1.2,
    "mm"
  ),

  cluster_rows = cluster_rows,
  row_dend_reorder = row_dend_reorder,

  row_gap = grid::unit(
    0,
    "mm"
  ),

  row_dend_width = grid::unit(
    18,
    "mm"
  ),

  show_column_names = TRUE,
  column_names_side = "top",
  column_names_rot = 45,
  column_names_centered = TRUE,

  column_names_gp = grid::gpar(
    fontsize = 9,
    fontface = "bold"
  ),

  column_names_max_height = grid::unit(
    18,
    "mm"
  ),

  show_row_names = FALSE,

  na_col = "grey90",

  rect_gp = grid::gpar(
    col = NA
  ),

  border = FALSE,
  use_raster = FALSE,

  width = grid::unit(
    heatmap_body_width_mm,
    "mm"
  ),

  height = grid::unit(
    heatmap_body_height_mm,
    "mm"
  ),

  column_title = heatmap_title,
  column_title_side = "top",

  column_title_gp = grid::gpar(
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

    legend_height = grid::unit(
      46,
      "mm"
    )
  )
)

disgenet_heatmap <- ComplexHeatmap::Heatmap(
  disgenet_binary_matrix,

  name = "DisGeNET overlap",
  col = disgenet_color_function,

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  show_row_dend = FALSE,
  show_column_dend = FALSE,

  show_row_names = FALSE,
  show_column_names = TRUE,

  column_names_side = "top",
  column_names_rot = 45,
  column_names_centered = TRUE,

  column_names_gp = grid::gpar(
    fontsize = 6.5,
    fontface = "bold"
  ),

  column_names_max_height = grid::unit(
    34,
    "mm"
  ),

  rect_gp = grid::gpar(
    col = "white",
    lwd = 0.35
  ),

  border = FALSE,
  use_raster = FALSE,

  width = grid::unit(
    8.0,
    "mm"
  ),

  height = grid::unit(
    heatmap_body_height_mm,
    "mm"
  ),

  column_title = "DisGeNET",
  column_title_side = "top",

  column_title_gp = grid::gpar(
    fontsize = 9,
    fontface = "bold"
  ),

  heatmap_legend_param = list(
    title = "DisGeNET overlap",

    at = c(
      0,
      1
    ),

    labels = c(
      "No overlap",
      "Overlap"
    ),

    title_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),

    labels_gp = grid::gpar(
      fontsize = 8
    )
  )
)

heatmap_object <-
  expression_heatmap +
  disgenet_heatmap +
  right_annotation


# ==============================================================================
# 29. Write PDF
# ==============================================================================

message("\nWriting heatmap:")
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
    fdr_marker_legend,
    log2fc_marker_legend
  ),

  ht_gap = grid::unit(
    1.2,
    "mm"
  ),

  padding = grid::unit(
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


# ==============================================================================
# 30. Results-per-cluster table
# ==============================================================================

all_cluster_ids <- names(
  custom_cluster_labels
)

results_per_cluster <- tibble::tibble(
  cluster_id = all_cluster_ids,

  anatomical_region = unname(
    custom_cluster_labels[
      all_cluster_ids
    ]
  )
) |>
  dplyr::left_join(
    significant_sex |>
      dplyr::count(
        .data$cluster_id,
        name = "n_results"
      ),
    by = "cluster_id"
  ) |>
  dplyr::left_join(
    significant_sex |>
      dplyr::mutate(
        direction = dplyr::case_when(
          .data$logFC > 0 ~ "UP",
          .data$logFC < 0 ~ "DOWN",
          TRUE ~ "ZERO"
        )
      ) |>
      dplyr::filter(
        .data$direction %in% c(
          "UP",
          "DOWN"
        )
      ) |>
      dplyr::count(
        .data$cluster_id,
        .data$direction,
        name = "n_direction"
      ) |>
      tidyr::pivot_wider(
        names_from = .data$direction,
        values_from = .data$n_direction,
        values_fill = 0L
      ),
    by = "cluster_id"
  ) |>
  dplyr::mutate(
    n_results = dplyr::coalesce(
      .data$n_results,
      0L
    ),

    UP = dplyr::coalesce(
      .data$UP,
      0L
    ),

    DOWN = dplyr::coalesce(
      .data$DOWN,
      0L
    ),

    cluster_id = as.integer(
      .data$cluster_id
    )
  ) |>
  dplyr::arrange(
    .data$cluster_id
  )

readr::write_tsv(
  results_per_cluster,
  results_per_cluster_tsv,
  na = "NA"
)

message("\nResults per cluster:")

print(
  results_per_cluster,
  n = Inf
)


# ==============================================================================
# 31. Validate outputs
# ==============================================================================

expected_output_files <- c(
  heatmap_pdf_file,
  results_per_cluster_tsv,
  threshold_counts_tsv
)

missing_output_files <- expected_output_files[
  !file.exists(
    expected_output_files
  )
]

if (length(missing_output_files) > 0L) {
  stop(
    "Missing output file(s):\n",
    paste(
      missing_output_files,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

empty_output_files <- expected_output_files[
  file.info(
    expected_output_files
  )$size == 0
]

if (length(empty_output_files) > 0L) {
  stop(
    "Empty output file(s):\n",
    paste(
      empty_output_files,
      collapse = "\n"
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 32. Final report
# ==============================================================================

message("\n")

message(
  paste(
    rep(
      "=",
      78
    ),
    collapse = ""
  )
)

message(
  "Sex main-effect long heatmap with DisGeNET Anxiety overlap completed successfully."
)

message(
  "\nTest:"
)

message(
  sex_test_id
)

message(
  "\nSignificant cluster-gene results:"
)

message(
  number_of_results
)

message(
  "\nUnique significant genes:"
)

message(
  number_of_unique_genes
)

message(
  "\nOutput directory:"
)

message(
  normalizePath(
    heatmap_output_dir,
    mustWork = TRUE
  )
)

message(
  "\nHeatmap:"
)

message(
  normalizePath(
    heatmap_pdf_file,
    mustWork = TRUE
  )
)

message(
  "\nResults-per-cluster table:"
)

message(
  normalizePath(
    results_per_cluster_tsv,
    mustWork = TRUE
  )
)

message(
  "\nThreshold-count table:"
)

message(
  normalizePath(
    threshold_counts_tsv,
    mustWork = TRUE
  )
)


# ==============================================================================
# End
# ==============================================================================
