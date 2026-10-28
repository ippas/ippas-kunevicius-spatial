#!/usr/bin/env Rscript

# ==============================================================================
# 04_PCA_crescendoCorrected_log2QN_topVariableFeatures_slideID.R
#
# Exploratory PCA for residual batch / slide effect after Crescendo correction.
#
# INPUT
#   Per-cluster matrices created by the Crescendo + edgeR standalone runner:
#
#     Crescendo-corrected pseudobulk counts
#       -> log2(count + 1)
#       -> quantile normalization across the 16 samples within each cluster
#
# PCA FEATURES
#   Each gene in each anatomical cluster is treated as a separate feature:
#
#     GeneA_C1
#     GeneA_C2
#     ...
#
#   This avoids averaging expression across anatomically distinct clusters.
#
#   Features with zero variance are removed.
#   The TOP 2000 most variable gene-cluster features are used for PCA.
#
# PCA
#   Samples = 16 biological samples
#   Features = top variable gene-cluster features
#   prcomp(center = TRUE, scale. = TRUE)
#
# PLOT
#   colour = slide_ID
#   shape  = donor group
#   label  = sample_ID
#
#   Sample 13_1F is additionally marked with a black outline because it is a
#   useful internal sanity check for residual slide effect.
#
# OUTPUT
#   1. PDF with PC1 vs PC2
#   2. TSV with PCA scores and sample metadata
#   3. TSV with the top variable gene-cluster features used for PCA
#
# IMPORTANT
#   This PCA uses expression values AFTER Crescendo batch correction.
#   It does NOT use the old raw pseudobulk expression matrices.
# ==============================================================================

options(stringsAsFactors = FALSE)

# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tibble",
  "readr",
  "ggplot2",
  "ggrepel"
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
  library(ggplot2)
  library(ggrepel)
})

# ==============================================================================
# 2. Configuration
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"

dataset_name <- "maternalFMT_n16samples"

clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

number_of_top_variable_features <- 2000L

highlight_sample <- "13_1F"

pca_center <- TRUE
pca_scale <- TRUE

# ==============================================================================
# 3. Paths
# ==============================================================================

statistics_root <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)

edgeR_statistics_dir <- file.path(
  statistics_root,
  "02_edgeR"
)

analysis_prefix <- paste0(
  dataset_name,
  "_",
  clustering_name,
  "_crescendoCorrectedPseudobulk_edgeR_twoFactorInteraction"
)

log2_qn_input_dir <- file.path(
  edgeR_statistics_dir,
  "01_pseudobulkInputs",
  "perCluster_crescendoCorrectedCounts_log2plus1_QN"
)

crescendo_input_rds <- file.path(
  statistics_root,
  "01_crescendoBatchCorrection",
  "07_objects",
  paste0(
    dataset_name,
    "_",
    clustering_name,
    "_crescendoBatchCorrection_pseudobulkCounts_beforeEdgeR.rds"
  )
)

pca_output_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,
  "02_heatmaps",
  "00_QC_residualBatchEffect",
  "01_PCA_crescendoCorrected"
)

dir.create(
  pca_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

pca_pdf_file <- file.path(
  pca_output_dir,
  paste0(
    "01_",
    dataset_name,
    "_crescendoCorrected_log2QN_top",
    number_of_top_variable_features,
    "_PCA_coloredBySlideID.pdf"
  )
)

pca_scores_tsv <- file.path(
  pca_output_dir,
  paste0(
    "02_",
    dataset_name,
    "_crescendoCorrected_log2QN_top",
    number_of_top_variable_features,
    "_PCA_scores.tsv"
  )
)

pca_features_tsv <- file.path(
  pca_output_dir,
  paste0(
    "03_",
    dataset_name,
    "_crescendoCorrected_log2QN_top",
    number_of_top_variable_features,
    "_PCA_features.tsv"
  )
)

# ==============================================================================
# 4. Validate inputs
# ==============================================================================

message("\n============================================================")
message("STAGE 1: VALIDATING INPUTS")
message("============================================================")

if (!dir.exists(log2_qn_input_dir)) {
  stop(
    "Missing Crescendo-corrected log2+QN directory:\n",
    log2_qn_input_dir,
    call. = FALSE
  )
}

if (!file.exists(crescendo_input_rds)) {
  stop(
    "Missing Crescendo RDS:\n",
    crescendo_input_rds,
    call. = FALSE
  )
}

# ==============================================================================
# 5. Load sample metadata from the Crescendo object
# ==============================================================================

message("\n============================================================")
message("STAGE 2: LOADING SAMPLE METADATA")
message("============================================================")

crescendo_object <- readRDS(
  crescendo_input_rds
)

if (!"sample_metadata" %in% names(crescendo_object)) {
  stop(
    "Crescendo RDS does not contain sample_metadata.",
    call. = FALSE
  )
}

sample_metadata_raw <- tibble::as_tibble(
  crescendo_object$sample_metadata
)

required_sample_columns <- c(
  "sample_ID",
  "fmt_donor_group",
  "sex"
)

missing_sample_columns <- setdiff(
  required_sample_columns,
  colnames(sample_metadata_raw)
)

if (length(missing_sample_columns) > 0L) {
  stop(
    "sample_metadata is missing required columns: ",
    paste(missing_sample_columns, collapse = ", "),
    call. = FALSE
  )
}

slide_column_candidates <- c(
  "slide_ID",
  "slide_id",
  "slide",
  "Slide_ID",
  "SlideID",
  "Slide"
)

slide_column <- slide_column_candidates[
  slide_column_candidates %in%
    colnames(sample_metadata_raw)
]

if (length(slide_column) == 0L) {
  message("\nAvailable sample_metadata columns:")
  message(
    paste(
      colnames(sample_metadata_raw),
      collapse = "\n"
    )
  )

  stop(
    "Could not identify the slide-ID column in sample_metadata.",
    call. = FALSE
  )
}

slide_column <- slide_column[[1L]]

message("Using slide column: ", slide_column)

sample_metadata <- sample_metadata_raw |>
  dplyr::transmute(
    sample_ID =
      as.character(.data$sample_ID),

    fmt_donor_group =
      as.character(.data$fmt_donor_group),

    sex =
      as.character(.data$sex),

    slide_ID =
      as.character(.data[[slide_column]])
  ) |>
  dplyr::distinct(
    .data$sample_ID,
    .keep_all = TRUE
  )

if (nrow(sample_metadata) != 16L) {
  stop(
    "Expected 16 samples but found ",
    nrow(sample_metadata),
    ".",
    call. = FALSE
  )
}

if (anyNA(sample_metadata$slide_ID)) {
  stop(
    "NA slide_ID values found.",
    call. = FALSE
  )
}

message("\nSample metadata:")

print(
  sample_metadata,
  n = Inf
)

# ==============================================================================
# 6. Find all per-cluster Crescendo-corrected log2+QN matrices
# ==============================================================================

message("\n============================================================")
message("STAGE 3: FINDING PER-CLUSTER CORRECTED MATRICES")
message("============================================================")

input_files <- list.files(
  log2_qn_input_dir,
  pattern =
    "_cluster_[0-9]+_crescendoCorrectedCounts_log2plus1_QN\\.tsv\\.gz$",
  full.names = TRUE
)

if (length(input_files) == 0L) {
  stop(
    "No Crescendo-corrected log2+QN matrices found in:\n",
    log2_qn_input_dir,
    call. = FALSE
  )
}

extract_cluster_id <- function(file_name) {
  sub(
    ".*_cluster_([0-9]+)_crescendoCorrectedCounts_log2plus1_QN\\.tsv\\.gz$",
    "\\1",
    basename(file_name)
  )
}

cluster_ids <- vapply(
  input_files,
  extract_cluster_id,
  FUN.VALUE = character(1)
)

cluster_order <- order(
  as.integer(cluster_ids)
)

input_files <- input_files[
  cluster_order
]

cluster_ids <- cluster_ids[
  cluster_order
]

if (length(cluster_ids) != 16L) {
  stop(
    "Expected 16 cluster matrices but found ",
    length(cluster_ids),
    ". Clusters found: ",
    paste(cluster_ids, collapse = ", "),
    call. = FALSE
  )
}

message(
  "Clusters found: ",
  paste(cluster_ids, collapse = ", ")
)

# ==============================================================================
# 7. Load all cluster matrices and build feature x sample matrix
# ==============================================================================

message("\n============================================================")
message("STAGE 4: BUILDING GENE-CLUSTER FEATURE MATRIX")
message("============================================================")

sample_order <- sample_metadata$sample_ID

cluster_expression_list <- vector(
  "list",
  length(input_files)
)

feature_metadata_list <- vector(
  "list",
  length(input_files)
)

for (file_index in seq_along(input_files)) {
  cluster_id_current <- cluster_ids[[file_index]]
  input_file_current <- input_files[[file_index]]

  message(
    "[",
    file_index,
    "/",
    length(input_files),
    "] Cluster ",
    cluster_id_current,
    ": ",
    basename(input_file_current)
  )

  input_table <- readr::read_tsv(
    input_file_current,
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
      " matrix is missing annotation columns: ",
      paste(
        missing_annotation_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  missing_samples <- setdiff(
    sample_order,
    colnames(input_table)
  )

  if (length(missing_samples) > 0L) {
    stop(
      "Cluster ",
      cluster_id_current,
      " is missing samples: ",
      paste(
        missing_samples,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  expression_matrix_current <- as.matrix(
    input_table[
      ,
      sample_order,
      drop = FALSE
    ]
  )

  storage.mode(
    expression_matrix_current
  ) <- "double"

  if (
    anyNA(expression_matrix_current) ||
      any(!is.finite(expression_matrix_current))
  ) {
    stop(
      "NA or non-finite expression values in cluster ",
      cluster_id_current,
      ".",
      call. = FALSE
    )
  }

  gene_display <- ifelse(
    is.na(input_table$gene) |
      input_table$gene == "",
    input_table$ensembl_gene_id,
    input_table$gene
  )

  feature_id <- paste0(
    input_table$ensembl_gene_id,
    "__C",
    cluster_id_current
  )

  rownames(
    expression_matrix_current
  ) <- feature_id

  cluster_expression_list[[file_index]] <- expression_matrix_current

  feature_metadata_list[[file_index]] <- tibble::tibble(
    feature_id = feature_id,
    cluster_id =
      as.character(cluster_id_current),
    ensembl_gene_id =
      as.character(input_table$ensembl_gene_id),
    gene =
      as.character(gene_display),
    chromosome =
      as.character(input_table$chromosome)
  )
}

expression_matrix_all <- do.call(
  rbind,
  cluster_expression_list
)

feature_metadata_all <- dplyr::bind_rows(
  feature_metadata_list
)

if (anyDuplicated(
  rownames(expression_matrix_all)
)) {
  stop(
    "Duplicated gene-cluster feature IDs found.",
    call. = FALSE
  )
}

message(
  "\nCombined matrix: ",
  format(
    nrow(expression_matrix_all),
    big.mark = ","
  ),
  " gene-cluster features x ",
  ncol(expression_matrix_all),
  " samples"
)

# ==============================================================================
# 8. Calculate feature variance and select top variable features
# ==============================================================================

message("\n============================================================")
message("STAGE 5: SELECTING TOP VARIABLE FEATURES")
message("============================================================")

feature_variance <- apply(
  expression_matrix_all,
  1L,
  stats::var
)

valid_features <- is.finite(
  feature_variance
) &
  feature_variance > 0

message(
  "Features with non-zero finite variance: ",
  format(
    sum(valid_features),
    big.mark = ","
  )
)

if (sum(valid_features) < 2L) {
  stop(
    "Too few variable features for PCA.",
    call. = FALSE
  )
}

number_of_features_to_use <- min(
  number_of_top_variable_features,
  sum(valid_features)
)

selected_feature_indices <- which(
  valid_features
)[
  order(
    feature_variance[valid_features],
    decreasing = TRUE
  )[
    seq_len(
      number_of_features_to_use
    )
  ]
]

selected_feature_ids <- rownames(
  expression_matrix_all
)[
  selected_feature_indices
]

expression_matrix_pca <- expression_matrix_all[
  selected_feature_indices,
  ,
  drop = FALSE
]

feature_metadata_selected <- feature_metadata_all |>
  dplyr::filter(
    .data$feature_id %in%
      selected_feature_ids
  ) |>
  dplyr::mutate(
    variance =
      feature_variance[
        match(
          .data$feature_id,
          rownames(expression_matrix_all)
        )
      ],
    variance_rank =
      match(
        .data$feature_id,
        selected_feature_ids
      )
  ) |>
  dplyr::arrange(
    .data$variance_rank
  )

message(
  "Features selected for PCA: ",
  nrow(expression_matrix_pca)
)

# ==============================================================================
# 9. PCA
# ==============================================================================

message("\n============================================================")
message("STAGE 6: RUNNING PCA")
message("============================================================")

pca_object <- stats::prcomp(
  t(expression_matrix_pca),
  center = pca_center,
  scale. = pca_scale
)

variance_explained <- (
  pca_object$sdev^2 /
    sum(pca_object$sdev^2)
) * 100

pca_scores <- as.data.frame(
  pca_object$x,
  stringsAsFactors = FALSE
) |>
  tibble::rownames_to_column(
    "sample_ID"
  ) |>
  tibble::as_tibble() |>
  dplyr::left_join(
    sample_metadata,
    by = "sample_ID"
  )

if (anyNA(pca_scores$slide_ID)) {
  stop(
    "Failed to attach sample metadata to PCA scores.",
    call. = FALSE
  )
}

message(
  "PC1 variance explained: ",
  round(
    variance_explained[[1]],
    2
  ),
  "%"
)

message(
  "PC2 variance explained: ",
  round(
    variance_explained[[2]],
    2
  ),
  "%"
)

# ==============================================================================
# 10. PCA plot
# ==============================================================================

message("\n============================================================")
message("STAGE 7: DRAWING PCA")
message("============================================================")

pca_plot <- ggplot2::ggplot(
  pca_scores,
  ggplot2::aes(
    x = .data$PC1,
    y = .data$PC2,
    colour = .data$slide_ID,
    shape = .data$fmt_donor_group
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.25,
    colour = "grey80"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.25,
    colour = "grey80"
  ) +
  ggplot2::geom_point(
    size = 4,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(
      label = .data$sample_ID
    ),
    size = 3.4,
    max.overlaps = Inf,
    box.padding = 0.45,
    point.padding = 0.30,
    min.segment.length = 0
  ) +
  ggplot2::labs(
    title =
      "PCA after Crescendo batch correction",
    subtitle = paste0(
      "Top ",
      number_of_features_to_use,
      " variable gene-cluster features | colour = slide ID | shape = donor group"
    ),
    x = paste0(
      "PC1 (",
      round(
        variance_explained[[1]],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(
        variance_explained[[2]],
        1
      ),
      "%)"
    ),
    colour = "Slide ID",
    shape = "FMT donor group"
  ) +
  ggplot2::theme_classic(
    base_size = 12
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    legend.title =
      ggplot2::element_text(
        face = "bold"
      )
  )

if (highlight_sample %in% pca_scores$sample_ID) {
  highlight_data <- pca_scores |>
    dplyr::filter(
      .data$sample_ID ==
        .env$highlight_sample
    )

  pca_plot <- pca_plot +
    ggplot2::geom_point(
      data = highlight_data,
      mapping = ggplot2::aes(
        x = .data$PC1,
        y = .data$PC2
      ),
      inherit.aes = FALSE,
      shape = 21,
      fill = NA,
      colour = "black",
      size = 6,
      stroke = 1.2
    )
}

# ==============================================================================
# 11. Save outputs
# ==============================================================================

message("\n============================================================")
message("STAGE 8: SAVING OUTPUTS")
message("============================================================")

grDevices::pdf(
  pca_pdf_file,
  width = 9.5,
  height = 7.5,
  onefile = TRUE,
  useDingbats = FALSE
)

print(
  pca_plot
)

grDevices::dev.off()

readr::write_tsv(
  pca_scores,
  pca_scores_tsv,
  na = "NA"
)

readr::write_tsv(
  feature_metadata_selected,
  pca_features_tsv,
  na = "NA"
)

# ==============================================================================
# 12. Final report
# ==============================================================================

message("\n============================================================")
message("PCA FINISHED")
message("============================================================")

message(
  "Samples: ",
  nrow(pca_scores)
)

message(
  "Gene-cluster features available: ",
  format(
    nrow(expression_matrix_all),
    big.mark = ","
  )
)

message(
  "Top variable features used: ",
  number_of_features_to_use
)

message(
  "PC1: ",
  round(
    variance_explained[[1]],
    2
  ),
  "%"
)

message(
  "PC2: ",
  round(
    variance_explained[[2]],
    2
  ),
  "%"
)

message("\nPCA PDF:")
message(
  normalizePath(
    pca_pdf_file,
    mustWork = TRUE
  )
)

message("\nPCA scores:")
message(
  normalizePath(
    pca_scores_tsv,
    mustWork = TRUE
  )
)

message("\nPCA features:")
message(
  normalizePath(
    pca_features_tsv,
    mustWork = TRUE
  )
)

# ==============================================================================
# End
# ==============================================================================
