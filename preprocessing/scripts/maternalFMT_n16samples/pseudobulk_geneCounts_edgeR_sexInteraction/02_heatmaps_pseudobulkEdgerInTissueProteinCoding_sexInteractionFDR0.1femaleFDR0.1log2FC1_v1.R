#!/usr/bin/env Rscript

# ==============================================================================
# 03_pseudobulk_edgeR_sexInteraction_twoHeatmaps_all16samples_sexGrouped_titleFIXED.R
#
# PURPOSE
#   Generate two independent heatmaps from the sex-interaction edgeR analysis:
#
#   1. INTERACTION HEATMAP
#      Selected genes:
#        FDR < 0.1
#
#      Row handling:
#        all selected genes are clustered together;
#        no up/down classification is used because interaction logFC is a
#        difference of differences.
#
#   2. FEMALE HEATMAP
#      Selected genes:
#        FDR < 0.1
#        abs(logFC) > 1
#
#      Row handling:
#        genes are divided internally according to the sign of Female logFC;
#        the two blocks are clustered separately;
#        no direction labels are displayed.
#
#   BOTH HEATMAPS
#      - show all 16 retained samples;
#      - therefore Male expression is shown in both heatmaps;
#      - order columns with the same sex next to each other:
#          Neurotypical Male
#          ASD Male
#          Neurotypical Female
#          ASD Female
#      - show donor group and sex as top annotations;
#      - do not cluster columns;
#      - use:
#          raw pseudobulk counts
#          -> log2(count + 1)
#          -> quantile normalization across all 16 samples
#          -> row-wise Z-score
#          -> clipping to [-2, 2].
#
# INPUT DIRECTORY
#   /home/mateusz/projects/ippas-kunevicius-spatial/results/
#   maternalFMT_n16samples/pseudobulk_geneCounts_edgeR_sexInteraction
#
# MAIN INPUT FILES
#   sexInteraction_interaction_fullResults.tsv
#   sexInteraction_ASD_vs_Neurotypical_in_Female_fullResults.tsv
#   sexInteraction_rawCounts_testedGenes.tsv
#   sexInteraction_sampleMetadata.tsv
#
# OUTPUT
#   Two PDF heatmaps plus selected-gene, sample-order and Z-score TSV files.
# ==============================================================================


# ==============================================================================
# 1. Check and load packages
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tibble",
  "limma",
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
    "Missing required packages: ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\n\nInstall them with:\n",
    "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) ",
    "install.packages(\"BiocManager\")\n",
    "BiocManager::install(c(\"limma\", \"ComplexHeatmap\"))\n",
    "install.packages(c(\"dplyr\", \"tibble\", \"circlize\"))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(limma)
  library(ComplexHeatmap)
  library(circlize)
})


# ==============================================================================
# 2. Define project, analysis and output directories
# ==============================================================================

project_dir <- "/home/mateusz/projects/ippas-kunevicius-spatial"

analysis_dir <- file.path(
  project_dir,
  "results",
  "maternalFMT_n16samples",
  "pseudobulk_geneCounts_edgeR_sexInteraction"
)

heatmap_output_dir <- file.path(
  analysis_dir,
  "edgeR_selectedGenes_twoHeatmaps_all16samples"
)

dir.create(
  heatmap_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(analysis_dir)) {
  stop(
    "Interaction-analysis directory does not exist: ",
    analysis_dir
  )
}

if (!dir.exists(heatmap_output_dir)) {
  stop(
    "Could not create heatmap output directory: ",
    heatmap_output_dir
  )
}


# ==============================================================================
# 3. Define and validate input files
# ==============================================================================

interaction_results_file <- file.path(
  analysis_dir,
  "sexInteraction_interaction_fullResults.tsv"
)

female_results_file <- file.path(
  analysis_dir,
  paste0(
    "sexInteraction_ASD_vs_Neurotypical_",
    "in_Female_fullResults.tsv"
  )
)

raw_counts_file <- file.path(
  analysis_dir,
  "sexInteraction_rawCounts_testedGenes.tsv"
)

sample_metadata_file <- file.path(
  analysis_dir,
  "sexInteraction_sampleMetadata.tsv"
)

required_input_files <- c(
  interaction_results_file,
  female_results_file,
  raw_counts_file,
  sample_metadata_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0L) {
  stop(
    "Missing required input files:\n",
    paste(
      missing_input_files,
      collapse = "\n"
    ),
    "\n\nRun the final sex-interaction edgeR script first."
  )
}

message("Interaction results: ", interaction_results_file)
message("Female results: ", female_results_file)
message("Raw counts: ", raw_counts_file)
message("Sample metadata: ", sample_metadata_file)
message("Heatmap output directory: ", heatmap_output_dir)


# ==============================================================================
# 4. Define output files
# ==============================================================================

interaction_heatmap_pdf <- file.path(
  heatmap_output_dir,
  paste0(
    "sexInteraction_interaction_",
    "FDR_lt0.1_all16samples_",
    "rawLog2Quantile_rowZscore_heatmap.pdf"
  )
)

female_heatmap_pdf <- file.path(
  heatmap_output_dir,
  paste0(
    "sexInteraction_female_ASD_vs_Neurotypical_",
    "FDR_lt0.1_absLog2FC_gt1_all16samples_",
    "rawLog2Quantile_rowZscore_heatmap.pdf"
  )
)

interaction_selected_genes_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "sexInteraction_interaction_",
    "FDR_lt0.1_selectedGenes.tsv"
  )
)

female_selected_genes_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "sexInteraction_female_ASD_vs_Neurotypical_",
    "FDR_lt0.1_absLog2FC_gt1_selectedGenes.tsv"
  )
)

interaction_zscore_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "sexInteraction_interaction_",
    "FDR_lt0.1_rowZscoreMatrix.tsv"
  )
)

female_zscore_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "sexInteraction_female_ASD_vs_Neurotypical_",
    "FDR_lt0.1_absLog2FC_gt1_rowZscoreMatrix.tsv"
  )
)

sample_order_tsv <- file.path(
  heatmap_output_dir,
  "sexInteraction_twoHeatmaps_sampleOrder.tsv"
)

heatmap_summary_tsv <- file.path(
  heatmap_output_dir,
  "sexInteraction_twoHeatmaps_summary.tsv"
)


# ==============================================================================
# 5. Read and validate interaction results
# ==============================================================================

interaction_results <- read.delim(
  interaction_results_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  tibble::as_tibble()

required_interaction_columns <- c(
  "ensembl_gene_id",
  "gene",
  "logFC",
  "logCPM",
  "F",
  "PValue",
  "FDR"
)

missing_interaction_columns <- setdiff(
  required_interaction_columns,
  colnames(interaction_results)
)

if (length(missing_interaction_columns) > 0L) {
  stop(
    "Missing columns in interaction results: ",
    paste(
      missing_interaction_columns,
      collapse = ", "
    )
  )
}

if (anyDuplicated(
  interaction_results$ensembl_gene_id
)) {
  stop(
    "Duplicated Ensembl IDs in interaction results."
  )
}

interaction_selected_results <- interaction_results |>
  dplyr::filter(
    !is.na(FDR),
    FDR < 0.1
  ) |>
  dplyr::arrange(
    FDR,
    PValue,
    dplyr::desc(abs(logFC))
  )

if (nrow(interaction_selected_results) == 0L) {
  stop(
    "No interaction genes meet FDR < 0.1."
  )
}

message(
  "\nInteraction genes selected: ",
  nrow(interaction_selected_results),
  " | FDR < 0.1"
)


# ==============================================================================
# 6. Read and validate Female ASD-vs-Neurotypical results
# ==============================================================================

female_results <- read.delim(
  female_results_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  tibble::as_tibble()

required_female_columns <- c(
  "ensembl_gene_id",
  "gene",
  "logFC",
  "logCPM",
  "F",
  "PValue",
  "FDR"
)

missing_female_columns <- setdiff(
  required_female_columns,
  colnames(female_results)
)

if (length(missing_female_columns) > 0L) {
  stop(
    "Missing columns in Female results: ",
    paste(
      missing_female_columns,
      collapse = ", "
    )
  )
}

if (anyDuplicated(
  female_results$ensembl_gene_id
)) {
  stop(
    "Duplicated Ensembl IDs in Female results."
  )
}

female_selected_results <- female_results |>
  dplyr::filter(
    !is.na(FDR),
    !is.na(logFC),
    FDR < 0.1,
    abs(logFC) > 1
  ) |>
  dplyr::mutate(
    logFC_block = dplyr::if_else(
      logFC > 0,
      "positive_logFC",
      "negative_logFC"
    ),
    logFC_block = factor(
      logFC_block,
      levels = c(
        "positive_logFC",
        "negative_logFC"
      )
    )
  ) |>
  dplyr::arrange(
    logFC_block,
    FDR,
    PValue,
    dplyr::desc(abs(logFC))
  )

if (nrow(female_selected_results) == 0L) {
  stop(
    paste0(
      "No Female ASD-vs-Neurotypical genes meet ",
      "FDR < 0.1 and abs(logFC) > 1."
    )
  )
}

message(
  "Female genes selected: ",
  nrow(female_selected_results),
  " | FDR < 0.1 and abs(logFC) > 1"
)

message(
  "  Positive Female logFC genes: ",
  sum(
    female_selected_results$logFC_block == "positive_logFC"
  )
)

message(
  "  Negative Female logFC genes: ",
  sum(
    female_selected_results$logFC_block == "negative_logFC"
  )
)


# ==============================================================================
# 7. Read, order and validate sample metadata
#
# All 16 samples are retained in BOTH heatmaps.
# ==============================================================================

sample_metadata <- read.delim(
  sample_metadata_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  tibble::as_tibble()

required_metadata_columns <- c(
  "sample_ID",
  "fmt_donor_group",
  "sex"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(sample_metadata)
)

if (length(missing_metadata_columns) > 0L) {
  stop(
    "Missing columns in sample metadata: ",
    paste(
      missing_metadata_columns,
      collapse = ", "
    )
  )
}

if (anyDuplicated(
  sample_metadata$sample_ID
)) {
  stop(
    "Duplicated sample_ID values in metadata."
  )
}

sample_metadata <- sample_metadata |>
  dplyr::mutate(
    sample_ID = trimws(
      as.character(sample_ID)
    ),
    fmt_donor_group = factor(
      trimws(
        as.character(fmt_donor_group)
      ),
      levels = c(
        "Neurotypical",
        "ASD"
      )
    ),
    sex = factor(
      trimws(
        as.character(sex)
      ),
      levels = c(
        "Male",
        "Female"
      )
    )
  )

if (anyNA(
  sample_metadata$fmt_donor_group
)) {
  stop(
    "Unexpected fmt_donor_group values in metadata."
  )
}

if (anyNA(
  sample_metadata$sex
)) {
  stop(
    "Unexpected sex values in metadata."
  )
}

sample_metadata <- sample_metadata |>
  dplyr::arrange(
    sex,
    fmt_donor_group,
    sample_ID
  ) |>
  dplyr::mutate(
    Group_Sex = factor(
      paste(
        fmt_donor_group,
        sex,
        sep = "_"
      ),
      levels = c(
        "Neurotypical_Male",
        "ASD_Male",
        "Neurotypical_Female",
        "ASD_Female"
      )
    )
  )

if (nrow(sample_metadata) != 16L) {
  stop(
    "Expected 16 retained samples, but metadata contains ",
    nrow(sample_metadata),
    "."
  )
}

sample_order <- as.character(
  sample_metadata$sample_ID
)

sample_order_output <- sample_metadata |>
  dplyr::transmute(
    sample_order = dplyr::row_number(),
    sample_ID,
    fmt_donor_group =
      as.character(
        fmt_donor_group
      ),
    sex =
      as.character(
        sex
      ),
    combined_group =
      as.character(
        Group_Sex
      )
  )

message("\nColumn order used in both heatmaps:")
print(
  sample_order_output,
  n = Inf
)


# ==============================================================================
# 8. Read and validate raw pseudobulk counts
# ==============================================================================

raw_counts_table <- read.delim(
  raw_counts_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  tibble::as_tibble()

required_count_identifier_columns <- c(
  "ensembl_gene_id",
  "gene"
)

missing_count_identifier_columns <- setdiff(
  required_count_identifier_columns,
  colnames(raw_counts_table)
)

if (
  length(
    missing_count_identifier_columns
  ) > 0L
) {
  stop(
    "Missing identifier columns in raw counts: ",
    paste(
      missing_count_identifier_columns,
      collapse = ", "
    )
  )
}

if (anyDuplicated(
  raw_counts_table$ensembl_gene_id
)) {
  stop(
    "Duplicated Ensembl IDs in raw counts."
  )
}

raw_count_sample_columns <- setdiff(
  colnames(raw_counts_table),
  required_count_identifier_columns
)

missing_samples_in_counts <- setdiff(
  sample_order,
  raw_count_sample_columns
)

if (length(missing_samples_in_counts) > 0L) {
  stop(
    "Samples missing from raw counts: ",
    paste(
      missing_samples_in_counts,
      collapse = ", "
    )
  )
}

unexpected_samples_in_counts <- setdiff(
  raw_count_sample_columns,
  sample_order
)

if (length(unexpected_samples_in_counts) > 0L) {
  stop(
    "Unexpected sample columns in raw counts: ",
    paste(
      unexpected_samples_in_counts,
      collapse = ", "
    )
  )
}

raw_count_matrix <- as.matrix(
  raw_counts_table[
    ,
    sample_order,
    drop = FALSE
  ]
)

storage.mode(
  raw_count_matrix
) <- "numeric"

rownames(
  raw_count_matrix
) <- raw_counts_table$ensembl_gene_id

if (anyNA(raw_count_matrix)) {
  stop(
    "Raw count matrix contains NA values."
  )
}

if (any(!is.finite(
  raw_count_matrix
))) {
  stop(
    "Raw count matrix contains non-finite values."
  )
}

if (any(raw_count_matrix < 0)) {
  stop(
    "Raw count matrix contains negative values."
  )
}

if (!all(
  raw_count_matrix ==
    floor(raw_count_matrix)
)) {
  stop(
    "Raw count matrix contains non-integer values."
  )
}


# ==============================================================================
# 9. Transform the complete matrix for visualization
#
# raw counts -> log2(count + 1) -> quantile normalization across all 16 samples
# ==============================================================================

log2_raw_count_matrix <- log2(
  raw_count_matrix + 1
)

log2_quantile_normalized_matrix <-
  limma::normalizeBetweenArrays(
    log2_raw_count_matrix,
    method = "quantile"
  )

rownames(
  log2_quantile_normalized_matrix
) <- rownames(raw_count_matrix)

colnames(
  log2_quantile_normalized_matrix
) <- colnames(raw_count_matrix)

if (anyNA(
  log2_quantile_normalized_matrix
)) {
  stop(
    "Log2 quantile-normalized matrix contains NA values."
  )
}

if (any(!is.finite(
  log2_quantile_normalized_matrix
))) {
  stop(
    paste0(
      "Log2 quantile-normalized matrix contains ",
      "non-finite values."
    )
  )
}


# ==============================================================================
# 10. Helper: prepare selected-gene heatmap matrix
# ==============================================================================

prepare_heatmap_matrix <- function(
    selected_results,
    normalized_matrix,
    sample_order,
    use_regulation_split
) {

  selected_gene_ids <- unique(
    as.character(
      selected_results$ensembl_gene_id
    )
  )

  missing_selected_genes <- setdiff(
    selected_gene_ids,
    rownames(normalized_matrix)
  )

  if (length(missing_selected_genes) > 0L) {
    stop(
      "Selected genes missing from expression matrix: ",
      paste(
        head(
          missing_selected_genes,
          20
        ),
        collapse = ", "
      ),
      if (
        length(missing_selected_genes) > 20L
      ) {
        " ..."
      } else {
        ""
      }
    )
  }

  selected_expression <- normalized_matrix[
    selected_gene_ids,
    sample_order,
    drop = FALSE
  ]

  zscore_matrix <- t(
    scale(
      t(selected_expression),
      center = TRUE,
      scale = TRUE
    )
  )

  zero_variance_gene_ids <- rownames(
    zscore_matrix
  )[
    apply(
      zscore_matrix,
      1,
      function(values) {
        any(!is.finite(values))
      }
    )
  ]

  if (length(zero_variance_gene_ids) > 0L) {
    warning(
      length(zero_variance_gene_ids),
      " selected genes had zero variance; ",
      "their Z-scores were set to 0."
    )
  }

  zscore_matrix[
    !is.finite(zscore_matrix)
  ] <- 0

  zscore_matrix_clipped <- pmax(
    pmin(
      zscore_matrix,
      2
    ),
    -2
  )

  selected_annotation <- selected_results[
    match(
      rownames(zscore_matrix_clipped),
      selected_results$ensembl_gene_id
    ),
    ,
    drop = FALSE
  ]

  if (anyNA(
    selected_annotation$ensembl_gene_id
  )) {
    stop(
      "Selected annotation could not be matched to matrix rows."
    )
  }

  row_labels <- ifelse(
    is.na(selected_annotation$gene) |
      selected_annotation$gene == "",
    selected_annotation$ensembl_gene_id,
    selected_annotation$gene
  )

  rownames(
    zscore_matrix_clipped
  ) <- make.unique(
    row_labels
  )

  row_split <- NULL

  if (isTRUE(use_regulation_split)) {

    if (!"logFC_block" %in%
        colnames(selected_annotation)) {
      stop(
        "A row split was requested, ",
        "but logFC_block is missing."
      )
    }

    row_split <- factor(
      as.character(
        selected_annotation$logFC_block
      ),
      levels = c(
        "positive_logFC",
        "negative_logFC"
      )
    )
  }

  list(
    selected_results =
      selected_results,
    selected_annotation =
      selected_annotation,
    selected_expression =
      selected_expression,
    zscore_matrix =
      zscore_matrix,
    zscore_matrix_clipped =
      zscore_matrix_clipped,
    row_split =
      row_split,
    zero_variance_gene_ids =
      zero_variance_gene_ids
  )
}


# ==============================================================================
# 11. Prepare interaction and Female heatmap matrices
# ==============================================================================

interaction_heatmap_data <- prepare_heatmap_matrix(
  selected_results =
    interaction_selected_results,
  normalized_matrix =
    log2_quantile_normalized_matrix,
  sample_order =
    sample_order,
  use_regulation_split =
    FALSE
)

female_heatmap_data <- prepare_heatmap_matrix(
  selected_results =
    female_selected_results,
  normalized_matrix =
    log2_quantile_normalized_matrix,
  sample_order =
    sample_order,
  use_regulation_split =
    TRUE
)


# ==============================================================================
# 12. Shared column annotations
# ==============================================================================

column_annotation_data <- data.frame(
  Group =
    sample_metadata$fmt_donor_group,
  Sex =
    sample_metadata$sex,
  Group_Sex =
    sample_metadata$Group_Sex,
  row.names =
    sample_order,
  check.names = FALSE
)

annotation_colors <- list(
  Group = c(
    Neurotypical = "grey70",
    ASD = "black"
  ),
  Sex = c(
    Male = "#3C78D8",
    Female = "#C050A0"
  )
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Group =
    column_annotation_data$Group,
  Sex =
    column_annotation_data$Sex,
  col =
    annotation_colors,
  simple_anno_size =
    grid::unit(
      4.2,
      "mm"
    ),
  annotation_name_side =
    "left",
  annotation_name_gp =
    grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),
  annotation_legend_param = list(
    Group = list(
      title = "FMT donor group",
      direction = "vertical",
      ncol = 1
    ),
    Sex = list(
      title = "Sex",
      direction = "vertical",
      ncol = 1
    )
  )
)


# ==============================================================================
# 13. Shared colors and clustering distance
# ==============================================================================

heatmap_color_function <- circlize::colorRamp2(
  breaks = c(
    -2,
    0,
    2
  ),
  colors = c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)

row_correlation_distance <- function(
    expression_matrix
) {

  correlation_matrix <- stats::cor(
    t(expression_matrix),
    method = "pearson",
    use = "pairwise.complete.obs"
  )

  correlation_matrix[
    !is.finite(correlation_matrix)
  ] <- 0

  diag(correlation_matrix) <- 1

  stats::as.dist(
    1 - correlation_matrix
  )
}


# ==============================================================================
# 14. Helper: calculate compact PDF dimensions
#
# The width calculation intentionally reserves extra space for long multi-line
# titles so that the left edge of the title is not clipped in the exported PDF.
# ==============================================================================

calculate_pdf_dimensions <- function(
    heatmap_matrix,
    row_font_size,
    row_height_mm,
    title_lines,
    heatmap_body_width_mm
) {

  temporary_device <- FALSE

  if (grDevices::dev.cur() == 1L) {
    grDevices::pdf(NULL)
    temporary_device <- TRUE
  }

  row_names_width_mm <- grid::convertWidth(
    ComplexHeatmap::max_text_width(
      rownames(heatmap_matrix),
      gp = grid::gpar(
        fontsize = row_font_size
      )
    ),
    unitTo = "mm",
    valueOnly = TRUE
  )

  title_width_mm <- grid::convertWidth(
    ComplexHeatmap::max_text_width(
      title_lines,
      gp = grid::gpar(
        fontsize = 12,
        fontface = "bold"
      )
    ),
    unitTo = "mm",
    valueOnly = TRUE
  )

  if (temporary_device) {
    grDevices::dev.off()
  }

  heatmap_body_height_mm <-
    nrow(heatmap_matrix) *
    row_height_mm

  row_dendrogram_width_mm <- 20
  legend_column_width_mm <- 28
  internal_spacing_width_mm <- 12

  plot_content_width_mm <-
    row_dendrogram_width_mm +
    heatmap_body_width_mm +
    row_names_width_mm +
    legend_column_width_mm +
    internal_spacing_width_mm

  # Add a larger width reserve for long multi-line titles and extra
  # horizontal breathing room so the first characters are never clipped.
  pdf_width_mm <- max(
    title_width_mm + 24,
    plot_content_width_mm + 12
  )

  # Slightly taller title block because titles are now explicitly split
  # across five lines.
  title_block_height_mm <- 34
  column_names_height_mm <- 17
  top_annotation_height_mm <- 10
  outer_vertical_padding_mm <- 6

  pdf_height_mm <-
    title_block_height_mm +
    column_names_height_mm +
    top_annotation_height_mm +
    heatmap_body_height_mm +
    outer_vertical_padding_mm

  list(
    pdf_width_inches =
      pdf_width_mm / 25.4,
    pdf_height_inches =
      pdf_height_mm / 25.4,
    heatmap_body_height_mm =
      heatmap_body_height_mm,
    row_dendrogram_width_mm =
      row_dendrogram_width_mm
  )
}


# ==============================================================================
# 15. Helper: create and save one heatmap
# ==============================================================================

save_one_heatmap <- function(
    heatmap_data,
    output_pdf,
    complete_title,
    title_lines_for_measurement,
    split_rows,
    rasterize_large_heatmap = FALSE
) {

  heatmap_matrix <-
    heatmap_data$zscore_matrix_clipped

  number_of_genes <-
    nrow(heatmap_matrix)

  number_of_samples <-
    ncol(heatmap_matrix)

  if (
    number_of_genes == 0L ||
    number_of_samples == 0L
  ) {
    stop(
      "Cannot create an empty heatmap."
    )
  }

  row_font_size <- if (
    number_of_genes <= 50L
  ) {
    8
  } else if (
    number_of_genes <= 120L
  ) {
    6.7
  } else if (
    number_of_genes <= 220L
  ) {
    5.3
  } else {
    4.2
  }

  row_height_mm <- if (
    number_of_genes <= 50L
  ) {
    3.8
  } else if (
    number_of_genes <= 120L
  ) {
    2.8
  } else if (
    number_of_genes <= 220L
  ) {
    2.2
  } else {
    1.75
  }

  heatmap_body_width_mm <-
    number_of_samples *
    6.2

  dimensions <- calculate_pdf_dimensions(
    heatmap_matrix =
      heatmap_matrix,
    row_font_size =
      row_font_size,
    row_height_mm =
      row_height_mm,
    title_lines =
      title_lines_for_measurement,
    heatmap_body_width_mm =
      heatmap_body_width_mm
  )

  heatmap_arguments <- list(
    matrix =
      heatmap_matrix,
    name =
      "Row Z-score",
    col =
      heatmap_color_function,

    cluster_columns =
      FALSE,
    column_split =
      column_annotation_data$Group_Sex,
    cluster_column_slices =
      FALSE,
    column_gap =
      grid::unit(
        1.2,
        "mm"
      ),

    cluster_rows =
      TRUE,
    clustering_distance_rows =
      row_correlation_distance,
    clustering_method_rows =
      "complete",
    row_dend_width =
      grid::unit(
        dimensions$row_dendrogram_width_mm,
        "mm"
      ),

    show_column_names =
      TRUE,
    column_names_side =
      "top",
    column_names_rot =
      45,
    column_names_centered =
      TRUE,
    column_names_gp =
      grid::gpar(
        fontsize = 9.5,
        fontface = "bold"
      ),
    column_names_max_height =
      grid::unit(
        18,
        "mm"
      ),

    show_row_names =
      TRUE,
    row_names_side =
      "right",
    row_names_gp =
      grid::gpar(
        fontsize = row_font_size
      ),
    row_names_max_width =
      ComplexHeatmap::max_text_width(
        rownames(heatmap_matrix),
        gp = grid::gpar(
          fontsize = row_font_size
        )
      ) +
      grid::unit(
        3,
        "mm"
      ),

    top_annotation =
      top_annotation,

    width =
      grid::unit(
        heatmap_body_width_mm,
        "mm"
      ),
    height =
      grid::unit(
        dimensions$heatmap_body_height_mm,
        "mm"
      ),

    rect_gp =
      grid::gpar(
        col = NA
      ),
    border =
      FALSE,
    use_raster =
      isTRUE(
        rasterize_large_heatmap
      ),

    column_title =
      complete_title,
    column_title_side =
      "top",
    column_title_gp =
      grid::gpar(
        fontsize = 12,
        fontface = "bold",
        lineheight = 1.12
      ),

    heatmap_legend_param = list(
      title =
        "Row Z-score",
      direction =
        "vertical",
      at = c(
        -2,
        -1,
        0,
        1,
        2
      ),
      labels = c(
        "-2",
        "-1",
        "0",
        "1",
        "2"
      ),
      title_gp =
        grid::gpar(
          fontsize = 10,
          fontface = "bold"
        ),
      labels_gp =
        grid::gpar(
          fontsize = 9
        ),
      legend_height =
        grid::unit(
          46,
          "mm"
        )
    )
  )

  if (isTRUE(split_rows)) {

    heatmap_arguments$row_split <-
      heatmap_data$row_split

    heatmap_arguments$cluster_row_slices <-
      FALSE

    heatmap_arguments$row_gap <-
      grid::unit(
        1.5,
        "mm"
      )

    # Preserve the two blocks, but hide the labels.
    heatmap_arguments$row_title <-
      NULL
  }

  heatmap_object <- do.call(
    ComplexHeatmap::Heatmap,
    heatmap_arguments
  )

  grDevices::pdf(
    file =
      output_pdf,
    width =
      dimensions$pdf_width_inches,
    height =
      dimensions$pdf_height_inches,
    onefile =
      TRUE,
    useDingbats =
      FALSE
  )

  ComplexHeatmap::draw(
    heatmap_object,
    heatmap_legend_side =
      "right",
    annotation_legend_side =
      "right",
    merge_legends =
      TRUE,
    legend_grouping =
      "original",
    padding =
      grid::unit(
        c(
          2,
          4,
          2,
          4
        ),
        "mm"
      )
  )

  grDevices::dev.off()

  if (!file.exists(output_pdf)) {
    stop(
      "Heatmap PDF was not created: ",
      output_pdf
    )
  }

  if (file.info(output_pdf)$size == 0) {
    stop(
      "Heatmap PDF is empty: ",
      output_pdf
    )
  }

  list(
    heatmap_object =
      heatmap_object,
    output_pdf =
      output_pdf,
    genes =
      number_of_genes,
    samples =
      number_of_samples,
    row_font_size =
      row_font_size,
    pdf_width_inches =
      dimensions$pdf_width_inches,
    pdf_height_inches =
      dimensions$pdf_height_inches
  )
}


# ==============================================================================
# 16. Save interaction heatmap
#
# The interaction genes are freely clustered with no up/down split.
# ==============================================================================

interaction_title_lines <- c(
  "Sex-by-group interaction",
  "Selected genes: interaction FDR < 0.1",
  "Input: raw pseudobulk counts",
  "-> log2(count + 1) -> quantile normalization",
  paste0(
    "All 16 samples shown | Row Z-score clipped to [-2, 2] | ",
    nrow(
      interaction_heatmap_data$zscore_matrix_clipped
    ),
    " genes"
  )
)

interaction_title <- paste(
  interaction_title_lines,
  collapse = "\n"
)

interaction_heatmap_plot <- save_one_heatmap(
  heatmap_data =
    interaction_heatmap_data,
  output_pdf =
    interaction_heatmap_pdf,
  complete_title =
    interaction_title,
  title_lines_for_measurement =
    interaction_title_lines,
  split_rows =
    FALSE,
  rasterize_large_heatmap =
    FALSE
)


# ==============================================================================
# 17. Save Female ASD-vs-Neurotypical heatmap
#
# Female-selected genes are split internally by logFC sign without labels.
# Both Female and Male sample expression is shown.
# ==============================================================================

female_title_lines <- c(
  "Female ASD vs Neurotypical selected genes",
  "Selected genes: Female FDR < 0.1 and |log2FC| > 1",
  "Input: raw pseudobulk counts",
  "-> log2(count + 1) -> quantile normalization",
  paste0(
    "All 16 samples shown | Row Z-score clipped to [-2, 2] | ",
    nrow(
      female_heatmap_data$zscore_matrix_clipped
    ),
    " genes"
  )
)

female_title <- paste(
  female_title_lines,
  collapse = "\n"
)

female_heatmap_plot <- save_one_heatmap(
  heatmap_data =
    female_heatmap_data,
  output_pdf =
    female_heatmap_pdf,
  complete_title =
    female_title,
  title_lines_for_measurement =
    female_title_lines,
  split_rows =
    TRUE,
  rasterize_large_heatmap =
    FALSE
)


# ==============================================================================
# 18. Save selected genes and row Z-score matrices
# ==============================================================================

interaction_selected_output <-
  interaction_selected_results |>
  dplyr::select(
    ensembl_gene_id,
    gene,
    logFC,
    logCPM,
    F,
    PValue,
    FDR
  )

female_selected_output <-
  female_selected_results |>
  dplyr::select(
    ensembl_gene_id,
    gene,
    logFC,
    logCPM,
    F,
    PValue,
    FDR
  )

interaction_zscore_output <- tibble::tibble(
  gene =
    rownames(
      interaction_heatmap_data$zscore_matrix_clipped
    )
) |>
  dplyr::bind_cols(
    tibble::as_tibble(
      as.data.frame(
        interaction_heatmap_data$zscore_matrix_clipped,
        check.names = FALSE
      ),
      .name_repair = "minimal"
    )
  )

female_zscore_output <- tibble::tibble(
  gene =
    rownames(
      female_heatmap_data$zscore_matrix_clipped
    )
) |>
  dplyr::bind_cols(
    tibble::as_tibble(
      as.data.frame(
        female_heatmap_data$zscore_matrix_clipped,
        check.names = FALSE
      ),
      .name_repair = "minimal"
    )
  )


# ==============================================================================
# 19. Save summary and all TSV outputs
# ==============================================================================

heatmap_summary <- tibble::tibble(
  heatmap = c(
    "Interaction",
    "Female_ASD_vs_Neurotypical"
  ),
  selection = c(
    "FDR < 0.1",
    "FDR < 0.1 and abs(logFC) > 1"
  ),
  selected_genes = c(
    nrow(interaction_selected_output),
    nrow(female_selected_output)
  ),
  positive_logFC_genes = c(
    NA_integer_,
    sum(
      female_selected_results$logFC_block == "positive_logFC"
    )
  ),
  negative_logFC_genes = c(
    NA_integer_,
    sum(
      female_selected_results$logFC_block == "negative_logFC"
    )
  ),
  samples_displayed = c(
    nrow(sample_metadata),
    nrow(sample_metadata)
  ),
  Male_samples_displayed = c(
    sum(sample_metadata$sex == "Male"),
    sum(sample_metadata$sex == "Male")
  ),
  Female_samples_displayed = c(
    sum(sample_metadata$sex == "Female"),
    sum(sample_metadata$sex == "Female")
  ),
  row_order = c(
    "free hierarchical clustering",
    "two unlabeled logFC-sign blocks; clustered within each block"
  ),
  column_order = paste(
    "Neurotypical Male -> ASD Male",
    "-> Neurotypical Female -> ASD Female"
  ),
  preprocessing = paste(
    "raw pseudobulk counts -> log2(count + 1)",
    "-> quantile normalization -> row Z-score"
  ),
  output_pdf = c(
    interaction_heatmap_pdf,
    female_heatmap_pdf
  )
)

write_tsv <- function(
    object,
    output_file
) {
  write.table(
    object,
    file =
      output_file,
    sep =
      "\t",
    quote =
      FALSE,
    row.names =
      FALSE,
    col.names =
      TRUE,
    na =
      "NA"
  )
}

write_tsv(
  interaction_selected_output,
  interaction_selected_genes_tsv
)

write_tsv(
  female_selected_output,
  female_selected_genes_tsv
)

write_tsv(
  interaction_zscore_output,
  interaction_zscore_tsv
)

write_tsv(
  female_zscore_output,
  female_zscore_tsv
)

write_tsv(
  sample_order_output,
  sample_order_tsv
)

write_tsv(
  heatmap_summary,
  heatmap_summary_tsv
)


# ==============================================================================
# 20. Validate saved files
# ==============================================================================

expected_output_files <- c(
  interaction_heatmap_pdf,
  female_heatmap_pdf,
  interaction_selected_genes_tsv,
  female_selected_genes_tsv,
  interaction_zscore_tsv,
  female_zscore_tsv,
  sample_order_tsv,
  heatmap_summary_tsv
)

missing_output_files <- expected_output_files[
  !file.exists(
    expected_output_files
  )
]

if (length(missing_output_files) > 0L) {
  stop(
    "Missing output files:\n",
    paste(
      missing_output_files,
      collapse = "\n"
    )
  )
}

empty_output_files <- expected_output_files[
  file.info(
    expected_output_files
  )$size == 0
]

if (length(empty_output_files) > 0L) {
  stop(
    "Empty output files:\n",
    paste(
      empty_output_files,
      collapse = "\n"
    )
  )
}


# ==============================================================================
# 21. Final report
# ==============================================================================

message(
  "\n",
  paste(rep("=", 80), collapse = "")
)

message(
  "TWO SEX-INTERACTION HEATMAPS COMPLETED SUCCESSFULLY"
)

message(
  paste(rep("=", 80), collapse = "")
)

message(
  "\nInteraction heatmap:"
)

message(
  "  selected genes: ",
  nrow(interaction_selected_output)
)

message(
  "  samples shown: ",
  nrow(sample_metadata),
  " including ",
  sum(sample_metadata$sex == "Male"),
  " Male and ",
  sum(sample_metadata$sex == "Female"),
  " Female"
)

message(
  "  PDF: ",
  normalizePath(
    interaction_heatmap_pdf,
    mustWork = TRUE
  )
)

message(
  "\nFemale-selected heatmap:"
)

message(
  "  selected genes: ",
  nrow(female_selected_output)
)

message(
  "  positive Female logFC genes: ",
  sum(
    female_selected_results$logFC_block == "positive_logFC"
  )
)

message(
  "  negative Female logFC genes: ",
  sum(
    female_selected_results$logFC_block == "negative_logFC"
  )
)

message(
  "  samples shown: ",
  nrow(sample_metadata),
  " including ",
  sum(sample_metadata$sex == "Male"),
  " Male and ",
  sum(sample_metadata$sex == "Female"),
  " Female"
)

message(
  "  PDF: ",
  normalizePath(
    female_heatmap_pdf,
    mustWork = TRUE
  )
)

message(
  "\nOutput directory: ",
  normalizePath(
    heatmap_output_dir,
    mustWork = TRUE
  )
)

# ==============================================================================
# End
# ==============================================================================
