#!/usr/bin/env Rscript

# ==============================================================================
# 03_pseudobulk_edgeR_allHeatmaps_all16samples_proteinCodingGenes_autosomalAndXChromosome.R
#
# PURPOSE
#   Generate all pseudobulk edgeR heatmaps using one common script:
#
#   1. Sex-by-donor-group interaction:
#        FDR < 0.1
#
#   2. Overall donor-group effect, ASD vs Neurotypical:
#        FDR < 0.1 and abs(logFC) > 0.5
#
#   3. Overall sex effect, Female vs Male:
#        FDR < 0.1 and abs(logFC) > 0.5
#
#   4. ASD vs Neurotypical within Female samples:
#        FDR < 0.1 and abs(logFC) > 1
#
#   The plotting logic and appearance are retained from the previous heatmap
#   scripts:
#      - all 16 retained samples are shown;
#      - columns are ordered:
#          Neurotypical Male
#          ASD Male
#          Neurotypical Female
#          ASD Female
#      - columns are not clustered;
#      - genes are clustered using Pearson correlation distance and complete
#        linkage;
#      - main-effect and Female-specific genes are split into two unlabeled
#        logFC-sign blocks and clustered separately within each block;
#      - interaction genes are clustered together without an up/down split;
#      - raw counts -> log2(count + 1) -> quantile normalization;
#      - row-wise Z-scores are clipped to [-2, 2];
#      - PDF dimensions are adapted to the number and length of gene labels.
#
# ANALYSIS INPUT DIRECTORY
#   /home/mateusz/projects/ippas-kunevicius-spatial/results/
#   maternalFMT_n16samples/pseudobulk_geneCounts_edgeR_sexInteraction/
#   proteinCodingGenes_autosomalAndXChromosome
#
# OUTPUT DIRECTORY
#   /home/mateusz/projects/ippas-kunevicius-spatial/results/
#   maternalFMT_n16samples/pseudobulk_geneCounts_edgeR_sexInteraction/
#   proteinCodingGenes_autosomalAndXChromosome/heatmaps
#
# This script generates the heatmaps only for protein-coding genes located on
# mouse autosomes 1-19 and chromosome X. All plotting, filtering, clustering,
# transformation and aesthetic settings are unchanged.
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
    paste(missing_packages, collapse = ", "),
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
# 2. Define project and analysis directories
# ==============================================================================

project_dir <- "/home/mateusz/projects/ippas-kunevicius-spatial"

analysis_parent_dir <- file.path(
  project_dir,
  "results",
  "maternalFMT_n16samples",
  "pseudobulk_geneCounts_edgeR_sexInteraction"
)

analysis_subfolders <- c(
  "proteinCodingGenes_autosomalAndXChromosome"
)

if (!dir.exists(project_dir)) {
  stop("Project directory does not exist: ", project_dir)
}

if (!dir.exists(analysis_parent_dir)) {
  stop(
    "Parent edgeR analysis directory does not exist: ",
    analysis_parent_dir
  )
}

message("Project directory: ", project_dir)
message("Parent edgeR directory: ", analysis_parent_dir)
message(
  "Requested analysis subfolders: ",
  paste(analysis_subfolders, collapse = ", ")
)


# ==============================================================================
# 3. Shared constants
# ==============================================================================

expected_number_of_samples <- 16L

expected_sample_groups <- c(
  "Neurotypical_Male",
  "ASD_Male",
  "Neurotypical_Female",
  "ASD_Female"
)

heatmap_selection_settings <- list(
  interaction = list(
    name = "Interaction",
    results_filename = "sexInteraction_interaction_fullResults.tsv",
    output_prefix = "sexInteraction_interaction_FDR_lt0.1",
    selection_description = "FDR < 0.1",
    fdr_threshold = 0.1,
    absolute_logfc_threshold = NA_real_,
    use_logfc_split = FALSE,
    title_line_1 = "Sex-by-group interaction",
    title_line_2 = "Selected genes: interaction FDR < 0.1"
  ),
  overall_group = list(
    name = "Overall_Group_ASD_vs_Neurotypical",
    results_filename =
      "sexInteraction_overall_ASD_vs_Neurotypical_fullResults.tsv",
    output_prefix =
      "sexInteraction_overall_ASD_vs_Neurotypical_FDR_lt0.1_absLog2FC_gt0.5",
    selection_description = "FDR < 0.1 and abs(logFC) > 0.5",
    fdr_threshold = 0.1,
    absolute_logfc_threshold = 0.5,
    use_logfc_split = TRUE,
    title_line_1 = "Overall donor-group effect: ASD vs Neurotypical",
    title_line_2 = "Selected genes: FDR < 0.1 and |log2FC| > 0.5"
  ),
  overall_sex = list(
    name = "Overall_Sex_Female_vs_Male",
    results_filename =
      "sexInteraction_overall_Female_vs_Male_fullResults.tsv",
    output_prefix =
      "sexInteraction_overall_Female_vs_Male_FDR_lt0.1_absLog2FC_gt0.5",
    selection_description = "FDR < 0.1 and abs(logFC) > 0.5",
    fdr_threshold = 0.1,
    absolute_logfc_threshold = 0.5,
    use_logfc_split = TRUE,
    title_line_1 = "Overall sex effect: Female vs Male",
    title_line_2 = "Selected genes: FDR < 0.1 and |log2FC| > 0.5"
  ),
  female_group = list(
    name = "Female_ASD_vs_Neurotypical",
    results_filename =
      "sexInteraction_ASD_vs_Neurotypical_in_Female_fullResults.tsv",
    output_prefix =
      "sexInteraction_female_ASD_vs_Neurotypical_FDR_lt0.1_absLog2FC_gt1",
    selection_description = "FDR < 0.1 and abs(logFC) > 1",
    fdr_threshold = 0.1,
    absolute_logfc_threshold = 1,
    use_logfc_split = TRUE,
    title_line_1 = "Female ASD vs Neurotypical selected genes",
    title_line_2 =
      "Selected genes: Female FDR < 0.1 and |log2FC| > 1"
  )
)


# ==============================================================================
# 4. General input and output helpers
# ==============================================================================

read_tsv_as_tibble <- function(
    input_file
) {

  read.delim(
    input_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    tibble::as_tibble()
}


write_tsv <- function(
    object,
    output_file
) {

  write.table(
    object,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}


validate_nonempty_files <- function(
    expected_files
) {

  missing_files <- expected_files[
    !file.exists(expected_files)
  ]

  if (length(missing_files) > 0L) {
    stop(
      "Missing output files:\n",
      paste(missing_files, collapse = "\n")
    )
  }

  empty_files <- expected_files[
    file.info(expected_files)$size == 0
  ]

  if (length(empty_files) > 0L) {
    stop(
      "Empty output files:\n",
      paste(empty_files, collapse = "\n")
    )
  }

  invisible(TRUE)
}


# ==============================================================================
# 5. Read and validate one edgeR result table
# ==============================================================================

read_and_select_result <- function(
    results_file,
    heatmap_setting
) {

  results <- read_tsv_as_tibble(
    results_file
  )

  required_result_columns <- c(
    "ensembl_gene_id",
    "gene",
    "logFC",
    "logCPM",
    "F",
    "PValue",
    "FDR"
  )

  missing_result_columns <- setdiff(
    required_result_columns,
    colnames(results)
  )

  if (length(missing_result_columns) > 0L) {
    stop(
      "Missing columns in result file ",
      results_file,
      ": ",
      paste(missing_result_columns, collapse = ", ")
    )
  }

  if (anyDuplicated(results$ensembl_gene_id)) {
    stop(
      "Duplicated Ensembl IDs in result file: ",
      results_file
    )
  }

  selected_results <- results |>
    dplyr::filter(
      !is.na(FDR),
      !is.na(logFC),
      FDR < heatmap_setting$fdr_threshold
    )

  if (!is.na(heatmap_setting$absolute_logfc_threshold)) {
    selected_results <- selected_results |>
      dplyr::filter(
        abs(logFC) >
          heatmap_setting$absolute_logfc_threshold
      )
  }

  if (isTRUE(heatmap_setting$use_logfc_split)) {
    selected_results <- selected_results |>
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
  } else {
    selected_results <- selected_results |>
      dplyr::arrange(
        FDR,
        PValue,
        dplyr::desc(abs(logFC))
      )
  }

  if (nrow(selected_results) == 0L) {
    stop(
      "No genes meet ",
      heatmap_setting$selection_description,
      " for ",
      heatmap_setting$name,
      "."
    )
  }

  message(
    heatmap_setting$name,
    " selected genes: ",
    nrow(selected_results),
    " | ",
    heatmap_setting$selection_description
  )

  selected_results
}


# ==============================================================================
# 6. Read, order and validate sample metadata
# ==============================================================================

read_and_order_sample_metadata <- function(
    sample_metadata_file
) {

  sample_metadata <- read_tsv_as_tibble(
    sample_metadata_file
  )

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
      paste(missing_metadata_columns, collapse = ", ")
    )
  }

  if (anyDuplicated(sample_metadata$sample_ID)) {
    stop("Duplicated sample_ID values in metadata.")
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

  if (anyNA(sample_metadata$sample_ID) ||
      any(sample_metadata$sample_ID == "")) {
    stop("Missing or empty sample_ID values in metadata.")
  }

  if (anyNA(sample_metadata$fmt_donor_group)) {
    stop("Unexpected fmt_donor_group values in metadata.")
  }

  if (anyNA(sample_metadata$sex)) {
    stop("Unexpected sex values in metadata.")
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
        levels = expected_sample_groups
      )
    )

  if (nrow(sample_metadata) != expected_number_of_samples) {
    stop(
      "Expected ",
      expected_number_of_samples,
      " retained samples, but metadata contains ",
      nrow(sample_metadata),
      "."
    )
  }

  if (anyNA(sample_metadata$Group_Sex)) {
    stop("Could not assign one or more samples to Group_Sex.")
  }

  sample_order_output <- sample_metadata |>
    dplyr::transmute(
      sample_order = dplyr::row_number(),
      sample_ID,
      fmt_donor_group =
        as.character(fmt_donor_group),
      sex =
        as.character(sex),
      combined_group =
        as.character(Group_Sex)
    )

  message("\nColumn order used in all heatmaps:")
  print(
    sample_order_output,
    n = Inf
  )

  list(
    metadata = sample_metadata,
    sample_order = as.character(
      sample_metadata$sample_ID
    ),
    sample_order_output =
      sample_order_output
  )
}


# ==============================================================================
# 7. Read and validate raw pseudobulk counts
# ==============================================================================

read_raw_count_matrix <- function(
    raw_counts_file,
    sample_order
) {

  raw_counts_table <- read_tsv_as_tibble(
    raw_counts_file
  )

  required_identifier_columns <- c(
    "ensembl_gene_id",
    "gene"
  )

  missing_identifier_columns <- setdiff(
    required_identifier_columns,
    colnames(raw_counts_table)
  )

  if (length(missing_identifier_columns) > 0L) {
    stop(
      "Missing identifier columns in raw counts: ",
      paste(missing_identifier_columns, collapse = ", ")
    )
  }

  if (anyDuplicated(raw_counts_table$ensembl_gene_id)) {
    stop("Duplicated Ensembl IDs in raw counts.")
  }

  raw_count_sample_columns <- setdiff(
    colnames(raw_counts_table),
    required_identifier_columns
  )

  missing_samples_in_counts <- setdiff(
    sample_order,
    raw_count_sample_columns
  )

  if (length(missing_samples_in_counts) > 0L) {
    stop(
      "Samples missing from raw counts: ",
      paste(missing_samples_in_counts, collapse = ", ")
    )
  }

  unexpected_samples_in_counts <- setdiff(
    raw_count_sample_columns,
    sample_order
  )

  if (length(unexpected_samples_in_counts) > 0L) {
    stop(
      "Unexpected sample columns in raw counts: ",
      paste(unexpected_samples_in_counts, collapse = ", ")
    )
  }

  raw_count_matrix <- as.matrix(
    raw_counts_table[
      ,
      sample_order,
      drop = FALSE
    ]
  )

  storage.mode(raw_count_matrix) <- "numeric"

  rownames(raw_count_matrix) <-
    raw_counts_table$ensembl_gene_id

  if (anyNA(raw_count_matrix)) {
    stop("Raw count matrix contains NA values.")
  }

  if (any(!is.finite(raw_count_matrix))) {
    stop("Raw count matrix contains non-finite values.")
  }

  if (any(raw_count_matrix < 0)) {
    stop("Raw count matrix contains negative values.")
  }

  if (!all(
    raw_count_matrix ==
      floor(raw_count_matrix)
  )) {
    stop("Raw count matrix contains non-integer values.")
  }

  raw_count_matrix
}


# ==============================================================================
# 8. Transform the complete matrix for visualization
#
# raw counts -> log2(count + 1) -> quantile normalization across all 16 samples
# ==============================================================================

create_log2_quantile_normalized_matrix <- function(
    raw_count_matrix
) {

  log2_raw_count_matrix <- log2(
    raw_count_matrix + 1
  )

  normalized_matrix <-
    limma::normalizeBetweenArrays(
      log2_raw_count_matrix,
      method = "quantile"
    )

  rownames(normalized_matrix) <-
    rownames(raw_count_matrix)

  colnames(normalized_matrix) <-
    colnames(raw_count_matrix)

  if (anyNA(normalized_matrix)) {
    stop(
      "Log2 quantile-normalized matrix contains NA values."
    )
  }

  if (any(!is.finite(normalized_matrix))) {
    stop(
      paste0(
        "Log2 quantile-normalized matrix contains ",
        "non-finite values."
      )
    )
  }

  normalized_matrix
}


# ==============================================================================
# 9. Prepare one selected-gene heatmap matrix
# ==============================================================================

prepare_heatmap_matrix <- function(
    selected_results,
    normalized_matrix,
    sample_order,
    use_logfc_split
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
        head(missing_selected_genes, 20),
        collapse = ", "
      ),
      if (length(missing_selected_genes) > 20L) {
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

  if (anyNA(selected_annotation$ensembl_gene_id)) {
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

  rownames(zscore_matrix_clipped) <- make.unique(
    row_labels
  )

  row_split <- NULL

  if (isTRUE(use_logfc_split)) {

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
# 10. Shared colors and clustering
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
# 11. Calculate compact PDF dimensions
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

  pdf_width_mm <- max(
    title_width_mm + 24,
    plot_content_width_mm + 12
  )

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
# 12. Create and save one heatmap
# ==============================================================================

save_one_heatmap <- function(
    heatmap_data,
    output_pdf,
    complete_title,
    title_lines_for_measurement,
    split_rows,
    column_annotation_data,
    top_annotation,
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
    stop("Cannot create an empty heatmap.")
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

    if (is.null(heatmap_data$row_split)) {
      stop(
        "split_rows is TRUE, but row_split is missing."
      )
    }

    heatmap_arguments$row_split <-
      heatmap_data$row_split

    heatmap_arguments$cluster_row_slices <-
      FALSE

    heatmap_arguments$row_gap <-
      grid::unit(
        1.5,
        "mm"
      )

    number_of_row_blocks <- nlevels(
      droplevels(
        heatmap_data$row_split
      )
    )

    heatmap_arguments$row_title <-
      rep(
        "",
        number_of_row_blocks
      )

    heatmap_arguments$row_title_gp <-
      grid::gpar(
        fontsize = 0
      )
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

  draw_succeeded <- FALSE

  tryCatch(
    {
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

      draw_succeeded <- TRUE
    },
    error = function(error_condition) {
      stop(
        "Heatmap drawing failed for ",
        output_pdf,
        ": ",
        conditionMessage(error_condition)
      )
    },
    finally = {
      grDevices::dev.off()
    }
  )

  if (!draw_succeeded) {
    stop(
      "Heatmap was not drawn successfully: ",
      output_pdf
    )
  }

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
# 13. Convert one heatmap matrix to a TSV-ready table
# ==============================================================================

create_zscore_output <- function(
    heatmap_data
) {

  tibble::tibble(
    gene =
      rownames(
        heatmap_data$zscore_matrix_clipped
      )
  ) |>
    dplyr::bind_cols(
      tibble::as_tibble(
        as.data.frame(
          heatmap_data$zscore_matrix_clipped,
          check.names = FALSE
        ),
        .name_repair = "minimal"
      )
    )
}


# ==============================================================================
# 14. Run all four heatmaps for one edgeR analysis directory
# ==============================================================================

run_heatmaps_for_analysis <- function(
    analysis_dir,
    analysis_label
) {

  message(
    "\n",
    paste(rep("=", 80), collapse = "")
  )

  message(
    "PROCESSING ANALYSIS: ",
    analysis_label
  )

  message(
    paste(rep("=", 80), collapse = "")
  )

  message("Analysis directory: ", analysis_dir)

  heatmap_output_dir <- file.path(
    analysis_dir,
    "heatmaps"
  )

  dir.create(
    heatmap_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!dir.exists(heatmap_output_dir)) {
    stop(
      "Could not create heatmap output directory: ",
      heatmap_output_dir
    )
  }

  raw_counts_file <- file.path(
    analysis_dir,
    "sexInteraction_rawCounts_testedGenes.tsv"
  )

  sample_metadata_file <- file.path(
    analysis_dir,
    "sexInteraction_sampleMetadata.tsv"
  )

  result_files <- vapply(
    heatmap_selection_settings,
    function(setting) {
      file.path(
        analysis_dir,
        setting$results_filename
      )
    },
    FUN.VALUE = character(1)
  )

  required_input_files <- c(
    result_files,
    raw_counts_file,
    sample_metadata_file
  )

  missing_input_files <- required_input_files[
    !file.exists(required_input_files)
  ]

  if (length(missing_input_files) > 0L) {
    stop(
      "Missing required input files for ",
      analysis_label,
      ":\n",
      paste(missing_input_files, collapse = "\n"),
      "\n\nRun the corresponding edgeR script first."
    )
  }

  message("\nInput files:")
  message("  Raw counts: ", raw_counts_file)
  message("  Sample metadata: ", sample_metadata_file)

  for (result_file in result_files) {
    message("  Result table: ", result_file)
  }

  message("  Heatmap output: ", heatmap_output_dir)

  metadata_data <- read_and_order_sample_metadata(
    sample_metadata_file
  )

  sample_metadata <- metadata_data$metadata
  sample_order <- metadata_data$sample_order
  sample_order_output <-
    metadata_data$sample_order_output

  raw_count_matrix <- read_raw_count_matrix(
    raw_counts_file =
      raw_counts_file,
    sample_order =
      sample_order
  )

  normalized_matrix <-
    create_log2_quantile_normalized_matrix(
      raw_count_matrix
    )

  selected_results <- list()
  heatmap_data <- list()
  heatmap_pdf_files <- character(0)
  selected_gene_files <- character(0)
  zscore_files <- character(0)
  saved_heatmap_info <- list()

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

  for (heatmap_id in names(
    heatmap_selection_settings
  )) {

    setting <-
      heatmap_selection_settings[[heatmap_id]]

    results_file <- file.path(
      analysis_dir,
      setting$results_filename
    )

    selected_results[[heatmap_id]] <-
      read_and_select_result(
        results_file =
          results_file,
        heatmap_setting =
          setting
      )

    heatmap_data[[heatmap_id]] <-
      prepare_heatmap_matrix(
        selected_results =
          selected_results[[heatmap_id]],
        normalized_matrix =
          normalized_matrix,
        sample_order =
          sample_order,
        use_logfc_split =
          setting$use_logfc_split
      )

    output_pdf <- file.path(
      heatmap_output_dir,
      paste0(
        setting$output_prefix,
        "_all16samples_rawLog2Quantile_",
        "rowZscore_heatmap.pdf"
      )
    )

    selected_genes_tsv <- file.path(
      heatmap_output_dir,
      paste0(
        setting$output_prefix,
        "_selectedGenes.tsv"
      )
    )

    zscore_tsv <- file.path(
      heatmap_output_dir,
      paste0(
        setting$output_prefix,
        "_rowZscoreMatrix.tsv"
      )
    )

    title_lines <- c(
      setting$title_line_1,
      setting$title_line_2,
      "Input: raw pseudobulk counts",
      "-> log2(count + 1) -> quantile normalization",
      paste0(
        "All 16 samples shown | Row Z-score clipped to [-2, 2] | ",
        nrow(
          heatmap_data[[heatmap_id]]$zscore_matrix_clipped
        ),
        " genes"
      )
    )

    complete_title <- paste(
      title_lines,
      collapse = "\n"
    )

    saved_heatmap_info[[heatmap_id]] <-
      save_one_heatmap(
        heatmap_data =
          heatmap_data[[heatmap_id]],
        output_pdf =
          output_pdf,
        complete_title =
          complete_title,
        title_lines_for_measurement =
          title_lines,
        split_rows =
          setting$use_logfc_split,
        column_annotation_data =
          column_annotation_data,
        top_annotation =
          top_annotation,
        rasterize_large_heatmap =
          FALSE
      )

    selected_output <-
      selected_results[[heatmap_id]] |>
      dplyr::select(
        ensembl_gene_id,
        gene,
        logFC,
        logCPM,
        F,
        PValue,
        FDR
      )

    zscore_output <- create_zscore_output(
      heatmap_data[[heatmap_id]]
    )

    write_tsv(
      selected_output,
      selected_genes_tsv
    )

    write_tsv(
      zscore_output,
      zscore_tsv
    )

    heatmap_pdf_files[[heatmap_id]] <-
      output_pdf

    selected_gene_files[[heatmap_id]] <-
      selected_genes_tsv

    zscore_files[[heatmap_id]] <-
      zscore_tsv
  }

  sample_order_tsv <- file.path(
    heatmap_output_dir,
    "sexInteraction_allHeatmaps_sampleOrder.tsv"
  )

  heatmap_summary_tsv <- file.path(
    heatmap_output_dir,
    "sexInteraction_allHeatmaps_summary.tsv"
  )

  write_tsv(
    sample_order_output,
    sample_order_tsv
  )

  heatmap_summary <- dplyr::bind_rows(
    lapply(
      names(heatmap_selection_settings),
      function(heatmap_id) {

        setting <-
          heatmap_selection_settings[[heatmap_id]]

        selected <-
          selected_results[[heatmap_id]]

        tibble::tibble(
          analysis =
            analysis_label,
          heatmap =
            setting$name,
          selection =
            setting$selection_description,
          selected_genes =
            nrow(selected),
          positive_logFC_genes =
            if (isTRUE(setting$use_logfc_split)) {
              sum(
                selected$logFC_block ==
                  "positive_logFC"
              )
            } else {
              NA_integer_
            },
          negative_logFC_genes =
            if (isTRUE(setting$use_logfc_split)) {
              sum(
                selected$logFC_block ==
                  "negative_logFC"
              )
            } else {
              NA_integer_
            },
          samples_displayed =
            nrow(sample_metadata),
          Male_samples_displayed =
            sum(sample_metadata$sex == "Male"),
          Female_samples_displayed =
            sum(sample_metadata$sex == "Female"),
          row_order =
            if (isTRUE(setting$use_logfc_split)) {
              paste(
                "two unlabeled logFC-sign blocks;",
                "clustered within each block"
              )
            } else {
              "free hierarchical clustering"
            },
          column_order = paste(
            "Neurotypical Male -> ASD Male",
            "-> Neurotypical Female -> ASD Female"
          ),
          preprocessing = paste(
            "raw pseudobulk counts -> log2(count + 1)",
            "-> quantile normalization -> row Z-score"
          ),
          output_pdf =
            heatmap_pdf_files[[heatmap_id]]
        )
      }
    )
  )

  write_tsv(
    heatmap_summary,
    heatmap_summary_tsv
  )

  expected_output_files <- c(
    unname(heatmap_pdf_files),
    unname(selected_gene_files),
    unname(zscore_files),
    sample_order_tsv,
    heatmap_summary_tsv
  )

  validate_nonempty_files(
    expected_output_files
  )

  message(
    "\n",
    paste(rep("-", 80), collapse = "")
  )

  message(
    "HEATMAPS COMPLETED FOR: ",
    analysis_label
  )

  message(
    paste(rep("-", 80), collapse = "")
  )

  for (heatmap_id in names(
    heatmap_selection_settings
  )) {

    setting <-
      heatmap_selection_settings[[heatmap_id]]

    message(
      "\n",
      setting$name,
      ":"
    )

    message(
      "  selected genes: ",
      nrow(
        selected_results[[heatmap_id]]
      )
    )

    if (isTRUE(setting$use_logfc_split)) {
      message(
        "  positive logFC genes: ",
        sum(
          selected_results[[heatmap_id]]$logFC_block ==
            "positive_logFC"
        )
      )

      message(
        "  negative logFC genes: ",
        sum(
          selected_results[[heatmap_id]]$logFC_block ==
            "negative_logFC"
        )
      )
    }

    message(
      "  PDF: ",
      normalizePath(
        heatmap_pdf_files[[heatmap_id]],
        mustWork = TRUE
      )
    )
  }

  message(
    "\nOutput directory: ",
    normalizePath(
      heatmap_output_dir,
      mustWork = TRUE
    )
  )

  invisible(
    list(
      analysis =
        analysis_label,
      analysis_dir =
        analysis_dir,
      heatmap_output_dir =
        heatmap_output_dir,
      heatmap_summary =
        heatmap_summary,
      expected_output_files =
        expected_output_files
    )
  )
}


# ==============================================================================
# 15. Run the common heatmap workflow for the autosomal-and-X analysis
# ==============================================================================

available_analysis_directories <- file.path(
  analysis_parent_dir,
  analysis_subfolders
)

analysis_directory_exists <- dir.exists(
  available_analysis_directories
)

if (!all(analysis_directory_exists)) {

  missing_analysis_directories <-
    available_analysis_directories[
      !analysis_directory_exists
    ]

  warning(
    "The following analysis directories do not exist and will be skipped:\n",
    paste(
      missing_analysis_directories,
      collapse = "\n"
    )
  )
}

available_analysis_directories <-
  available_analysis_directories[
    analysis_directory_exists
  ]

available_analysis_labels <-
  analysis_subfolders[
    analysis_directory_exists
  ]

if (length(available_analysis_directories) == 0L) {
  stop(
    "None of the requested edgeR analysis directories exists."
  )
}

all_heatmap_runs <- vector(
  mode = "list",
  length = length(
    available_analysis_directories
  )
)

names(all_heatmap_runs) <-
  available_analysis_labels

for (analysis_index in seq_along(
  available_analysis_directories
)) {

  all_heatmap_runs[[analysis_index]] <-
    run_heatmaps_for_analysis(
      analysis_dir =
        available_analysis_directories[
          analysis_index
        ],
      analysis_label =
        available_analysis_labels[
          analysis_index
        ]
    )
}


# ==============================================================================
# 16. Final report
# ==============================================================================

message(
  "\n",
  paste(rep("=", 80), collapse = "")
)

message(
  "AUTOSOMAL-AND-X HEATMAP WORKFLOW COMPLETED SUCCESSFULLY"
)

message(
  paste(rep("=", 80), collapse = "")
)

message(
  "\nProcessed analyses: ",
  paste(
    names(all_heatmap_runs),
    collapse = ", "
  )
)

for (analysis_name in names(
  all_heatmap_runs
)) {
  message(
    "  ",
    analysis_name,
    ": ",
    normalizePath(
      all_heatmap_runs[[analysis_name]]$heatmap_output_dir,
      mustWork = TRUE
    )
  )
}

# ==============================================================================
# End
# ==============================================================================
