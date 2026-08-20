#!/usr/bin/env Rscript

# ==============================================================================
# Cluster-to-cluster correlation heatmaps for the main sex effect
#
# TEST
#   Overall_Sex_Female_vs_Male
#
# INTERPRETATION
#   Positive logFC = higher expression in Female
#   Negative logFC = higher expression in Male
#
# Four variants are generated from cluster-specific edgeR log2FC values:
#   1. ALL tested genes      + Pearson
#   2. ALL tested genes      + Spearman
#   3. SIGNIFICANT gene union + Pearson
#   4. SIGNIFICANT gene union + Spearman
#
# SIGNIFICANT gene union:
#   unique genes significant in at least one cluster at:
#     FDR < 0.05 and |log2FC| >= 0.5
#   Once this union is defined, the log2FC value is taken from every cluster in
#   which that gene was tested, regardless of significance in that cluster.
#
# Missing cluster-gene log2FC values:
#   retained as NA; no imputation; pairwise complete observations are used.
#
# Each variant writes:
#   - one PDF heatmap
#   - one XLSX workbook with:
#       correlation_matrix
#       pairwise_statistics
#       analysis_info
# ==============================================================================

options(stringsAsFactors = FALSE)


# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tibble",
  "tidyr",
  "ComplexHeatmap",
  "circlize",
  "grid",
  "openxlsx"
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
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(openxlsx)
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
p_adjust_method <- "BH"
linkage_method <- "complete"

clustering_title <- paste(
  "LogNormalize + VST | HVG = 2000 | CCA dims = 20 |",
  "k = 20 | prune.SNN = 0.0667 | Leiden resolution = 0.40"
)

correlation_variants <- tibble::tribble(
  ~file_index, ~gene_set,      ~method,     ~variant_id,
  1L,          "all",          "pearson",  "allTestedGenes_pearson",
  2L,          "all",          "spearman", "allTestedGenes_spearman",
  3L,          "significant",  "pearson",  "significantUnion_pearson",
  4L,          "significant",  "spearman", "significantUnion_spearman"
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

correlation_output_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "05_correlationHeatmap",
  "03_mainEffects_sex",
  "01_log2FCcorrelation_allVsSignificant"
)

dir.create(
  correlation_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(edgeR_results_rdata_file)) {
  stop(
    "Missing edgeR results RData file:\n",
    edgeR_results_rdata_file,
    call. = FALSE
  )
}


# ==============================================================================
# 5. Helper functions
# ==============================================================================

sort_cluster_ids <- function(cluster_ids) {

  cluster_ids <- unique(as.character(cluster_ids))

  cluster_numbers <- suppressWarnings(
    as.integer(cluster_ids)
  )

  if (anyNA(cluster_numbers)) {
    stop(
      "Cluster IDs must be integer-like values. Found: ",
      paste(cluster_ids[is.na(cluster_numbers)], collapse = ", "),
      call. = FALSE
    )
  }

  cluster_ids[order(cluster_numbers)]
}


build_logfc_matrix <- function(results_table, selected_gene_ids, cluster_ids) {

  selected_table <- results_table |>
    dplyr::filter(
      .data$ensembl_gene_id %in% selected_gene_ids,
      .data$cluster_id %in% cluster_ids
    ) |>
    dplyr::select(
      .data$ensembl_gene_id,
      .data$cluster_id,
      .data$logFC
    )

  duplicate_rows <- selected_table |>
    dplyr::count(
      .data$cluster_id,
      .data$ensembl_gene_id,
      name = "n_rows"
    ) |>
    dplyr::filter(.data$n_rows > 1L)

  if (nrow(duplicate_rows) > 0L) {
    stop(
      "Duplicated cluster-gene logFC rows detected. Example: cluster ",
      duplicate_rows$cluster_id[[1]],
      ", gene ",
      duplicate_rows$ensembl_gene_id[[1]],
      ".",
      call. = FALSE
    )
  }

  wide_table <- selected_table |>
    tidyr::pivot_wider(
      names_from = .data$cluster_id,
      values_from = .data$logFC
    )

  missing_cluster_columns <- setdiff(
    cluster_ids,
    colnames(wide_table)
  )

  if (length(missing_cluster_columns) > 0L) {
    for (cluster_id_current in missing_cluster_columns) {
      wide_table[[cluster_id_current]] <- NA_real_
    }
  }

  wide_table <- wide_table |>
    dplyr::arrange(.data$ensembl_gene_id)

  logfc_matrix <- as.matrix(
    wide_table[, cluster_ids, drop = FALSE]
  )

  storage.mode(logfc_matrix) <- "double"
  rownames(logfc_matrix) <- wide_table$ensembl_gene_id

  logfc_matrix
}


calculate_pairwise_correlations <- function(
  logfc_matrix,
  cluster_ids,
  correlation_method
) {

  if (!correlation_method %in% c("pearson", "spearman")) {
    stop(
      "Unknown correlation method: ",
      correlation_method,
      call. = FALSE
    )
  }

  n_clusters <- length(cluster_ids)
  display_cluster_ids <- paste0("C", cluster_ids)

  correlation_matrix <- matrix(
    NA_real_,
    nrow = n_clusters,
    ncol = n_clusters,
    dimnames = list(
      display_cluster_ids,
      display_cluster_ids
    )
  )

  diag(correlation_matrix) <- 1

  pairwise_results <- vector(
    mode = "list",
    length = choose(n_clusters, 2)
  )

  result_index <- 0L

  for (cluster_A_index in seq_len(n_clusters - 1L)) {
    for (cluster_B_index in seq.int(cluster_A_index + 1L, n_clusters)) {

      result_index <- result_index + 1L

      cluster_A_id_raw <- cluster_ids[[cluster_A_index]]
      cluster_B_id_raw <- cluster_ids[[cluster_B_index]]

      x <- logfc_matrix[, cluster_A_id_raw]
      y <- logfc_matrix[, cluster_B_id_raw]

      complete_index <-
        is.finite(x) &
        is.finite(y)

      x_complete <- x[complete_index]
      y_complete <- y[complete_index]
      n_genes <- length(x_complete)

      correlation_value <- NA_real_
      p_value <- NA_real_

      sufficient_data <-
        n_genes >= 3L &&
        length(unique(x_complete)) >= 2L &&
        length(unique(y_complete)) >= 2L

      if (sufficient_data) {

        correlation_test <- if (identical(correlation_method, "spearman")) {
          suppressWarnings(
            stats::cor.test(
              x_complete,
              y_complete,
              method = "spearman",
              exact = FALSE
            )
          )
        } else {
          stats::cor.test(
            x_complete,
            y_complete,
            method = "pearson"
          )
        }

        correlation_value <- unname(
          correlation_test$estimate[[1]]
        )

        p_value <- correlation_test$p.value
      }

      correlation_matrix[
        cluster_A_index,
        cluster_B_index
      ] <- correlation_value

      correlation_matrix[
        cluster_B_index,
        cluster_A_index
      ] <- correlation_value

      pairwise_results[[result_index]] <- tibble::tibble(
        cluster_A_id = paste0("C", cluster_A_id_raw),
        cluster_A_region = unname(custom_cluster_labels[cluster_A_id_raw]),
        cluster_B_id = paste0("C", cluster_B_id_raw),
        cluster_B_region = unname(custom_cluster_labels[cluster_B_id_raw]),
        correlation = correlation_value,
        n_genes = n_genes,
        p_value = p_value
      )
    }
  }

  pairwise_statistics <- dplyr::bind_rows(
    pairwise_results
  ) |>
    dplyr::mutate(
      FDR = stats::p.adjust(
        .data$p_value,
        method = p_adjust_method
      )
    )

  list(
    correlation_matrix = correlation_matrix,
    pairwise_statistics = pairwise_statistics
  )
}


write_correlation_xlsx <- function(
  xlsx_file,
  correlation_matrix,
  pairwise_statistics,
  analysis_info,
  cluster_ids
) {

  workbook <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(
    workbook,
    "correlation_matrix"
  )

  openxlsx::addWorksheet(
    workbook,
    "pairwise_statistics"
  )

  openxlsx::addWorksheet(
    workbook,
    "analysis_info"
  )

  matrix_table <- data.frame(
    cluster_id = rownames(correlation_matrix),
    anatomical_region = unname(custom_cluster_labels[cluster_ids]),
    correlation_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  openxlsx::writeData(
    workbook,
    sheet = "correlation_matrix",
    x = matrix_table,
    withFilter = FALSE
  )

  openxlsx::writeData(
    workbook,
    sheet = "pairwise_statistics",
    x = pairwise_statistics,
    withFilter = TRUE
  )

  openxlsx::writeData(
    workbook,
    sheet = "analysis_info",
    x = analysis_info,
    withFilter = FALSE
  )

  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    border = "Bottom"
  )

  openxlsx::addStyle(
    workbook,
    sheet = "correlation_matrix",
    style = header_style,
    rows = 1,
    cols = seq_len(ncol(matrix_table)),
    gridExpand = TRUE
  )

  openxlsx::addStyle(
    workbook,
    sheet = "pairwise_statistics",
    style = header_style,
    rows = 1,
    cols = seq_len(ncol(pairwise_statistics)),
    gridExpand = TRUE
  )

  openxlsx::addStyle(
    workbook,
    sheet = "analysis_info",
    style = header_style,
    rows = 1,
    cols = seq_len(ncol(analysis_info)),
    gridExpand = TRUE
  )

  correlation_style <- openxlsx::createStyle(
    numFmt = "0.0000"
  )

  pvalue_style <- openxlsx::createStyle(
    numFmt = "0.00E+00"
  )

  openxlsx::addStyle(
    workbook,
    sheet = "correlation_matrix",
    style = correlation_style,
    rows = 2:(nrow(matrix_table) + 1L),
    cols = 3:ncol(matrix_table),
    gridExpand = TRUE,
    stack = TRUE
  )

  correlation_col <- match(
    "correlation",
    colnames(pairwise_statistics)
  )

  p_value_col <- match(
    "p_value",
    colnames(pairwise_statistics)
  )

  fdr_col <- match(
    "FDR",
    colnames(pairwise_statistics)
  )

  openxlsx::addStyle(
    workbook,
    sheet = "pairwise_statistics",
    style = correlation_style,
    rows = 2:(nrow(pairwise_statistics) + 1L),
    cols = correlation_col,
    gridExpand = TRUE,
    stack = TRUE
  )

  openxlsx::addStyle(
    workbook,
    sheet = "pairwise_statistics",
    style = pvalue_style,
    rows = 2:(nrow(pairwise_statistics) + 1L),
    cols = c(p_value_col, fdr_col),
    gridExpand = TRUE,
    stack = TRUE
  )

  openxlsx::freezePane(
    workbook,
    sheet = "correlation_matrix",
    firstRow = TRUE,
    firstCol = TRUE
  )

  openxlsx::freezePane(
    workbook,
    sheet = "pairwise_statistics",
    firstRow = TRUE
  )

  openxlsx::setColWidths(
    workbook,
    sheet = "correlation_matrix",
    cols = 1,
    widths = 12
  )

  openxlsx::setColWidths(
    workbook,
    sheet = "correlation_matrix",
    cols = 2,
    widths = 45
  )

  openxlsx::setColWidths(
    workbook,
    sheet = "correlation_matrix",
    cols = 3:ncol(matrix_table),
    widths = 11
  )

  openxlsx::setColWidths(
    workbook,
    sheet = "pairwise_statistics",
    cols = 1:ncol(pairwise_statistics),
    widths = "auto"
  )

  openxlsx::setColWidths(
    workbook,
    sheet = "analysis_info",
    cols = 1,
    widths = 34
  )

  openxlsx::setColWidths(
    workbook,
    sheet = "analysis_info",
    cols = 2,
    widths = 95
  )

  openxlsx::saveWorkbook(
    workbook,
    xlsx_file,
    overwrite = TRUE
  )
}


make_cluster_legend <- function(cluster_ids) {

  display_cluster_ids <- paste0("C", cluster_ids)

  ComplexHeatmap::Legend(
    title = "Cluster / anatomical region",
    labels = paste0(
      display_cluster_ids,
      ": ",
      unname(custom_cluster_labels[cluster_ids])
    ),
    legend_gp = grid::gpar(
      fill = unname(custom_cluster_colors[cluster_ids]),
      col = NA
    ),
    ncol = 1,
    title_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    labels_gp = grid::gpar(
      fontsize = 8
    ),
    grid_width = grid::unit(4.2, "mm"),
    grid_height = grid::unit(4.2, "mm")
  )
}


write_correlation_heatmap <- function(
  pdf_file,
  correlation_matrix,
  cluster_ids,
  correlation_method,
  heatmap_title
) {

  display_cluster_ids <- paste0("C", cluster_ids)

  cluster_colors_display <- stats::setNames(
    unname(custom_cluster_colors[cluster_ids]),
    display_cluster_ids
  )

  correlation_color_function <- circlize::colorRamp2(
    c(-1, 0, 1),
    c(
      "#2166AC",
      "white",
      "#B2182B"
    )
  )

  cluster_factor <- factor(
    display_cluster_ids,
    levels = display_cluster_ids
  )

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Cluster = cluster_factor,
    col = list(
      Cluster = cluster_colors_display
    ),
    simple_anno_size = grid::unit(4.2, "mm"),
    show_legend = FALSE,
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    )
  )

  left_annotation <- ComplexHeatmap::rowAnnotation(
    Cluster = cluster_factor,
    col = list(
      Cluster = cluster_colors_display
    ),
    simple_anno_size = grid::unit(4.2, "mm"),
    show_legend = FALSE,
    annotation_name_side = "top",
    annotation_name_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    )
  )

  finite_off_diagonal <- correlation_matrix[
    row(correlation_matrix) != col(correlation_matrix)
  ]

  can_cluster <- all(
    is.finite(finite_off_diagonal)
  )

  if (can_cluster) {

    distance_matrix <- 1 - correlation_matrix
    diag(distance_matrix) <- 0

    cluster_hclust <- stats::hclust(
      stats::as.dist(distance_matrix),
      method = linkage_method
    )

    cluster_rows_current <- cluster_hclust
    cluster_columns_current <- cluster_hclust

  } else {

    warning(
      "At least one pairwise correlation is NA/non-finite. ",
      "Heatmap clustering is disabled for this variant; original cluster order is used."
    )

    cluster_rows_current <- FALSE
    cluster_columns_current <- FALSE
  }

  legend_title <- if (identical(correlation_method, "pearson")) {
    "Pearson r"
  } else {
    "Spearman rho"
  }

  heatmap_object <- ComplexHeatmap::Heatmap(
    correlation_matrix,
    name = legend_title,
    col = correlation_color_function,
    left_annotation = left_annotation,
    top_annotation = top_annotation,
    cluster_rows = cluster_rows_current,
    cluster_columns = cluster_columns_current,
    row_dend_reorder = FALSE,
    column_dend_reorder = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "left",
    column_names_side = "top",
    column_names_rot = 45,
    row_names_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    column_names_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    na_col = "grey90",
    rect_gp = grid::gpar(
      col = NA
    ),
    border = FALSE,
    use_raster = FALSE,
    width = grid::unit(115, "mm"),
    height = grid::unit(115, "mm"),
    column_title = heatmap_title,
    column_title_side = "top",
    column_title_gp = grid::gpar(
      fontsize = 12,
      fontface = "bold",
      lineheight = 1.12
    ),
    heatmap_legend_param = list(
      title = legend_title,
      at = c(-1, -0.5, 0, 0.5, 1),
      legend_height = grid::unit(46, "mm")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {

      value <- correlation_matrix[i, j]

      if (is.finite(value)) {

        text_color <- if (abs(value) >= 0.65) {
          "white"
        } else {
          "black"
        }

        grid::grid.text(
          label = formatC(
            value,
            format = "f",
            digits = 2
          ),
          x = x,
          y = y,
          gp = grid::gpar(
            col = text_color,
            fontsize = 7.5,
            fontface = "bold"
          )
        )
      }
    }
  )

  cluster_legend <- make_cluster_legend(
    cluster_ids
  )

  grDevices::pdf(
    file = pdf_file,
    width = 14,
    height = 12.5,
    onefile = TRUE,
    useDingbats = FALSE
  )

  ComplexHeatmap::draw(
    heatmap_object,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    annotation_legend_list = list(
      cluster_legend
    ),
    merge_legends = TRUE,
    legend_grouping = "original",
    padding = grid::unit(
      c(3, 3, 3, 3),
      "mm"
    )
  )

  grDevices::dev.off()
}


# ==============================================================================
# 6. Load edgeR results
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
    "Required object(s) missing after loading RData: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}


# ==============================================================================
# 7. Select main sex-effect test
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
  dplyr::filter(.data$test_id == sex_test_id)

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
# 8. Prepare cluster-specific sex-effect log2FC results and gene sets
# ==============================================================================

sex_results <- edgeR_perCluster_combinedResults[[sex_test_id]] |>
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
  colnames(sex_results)
)

if (length(missing_result_columns) > 0L) {
  stop(
    "Sex-effect results are missing required columns: ",
    paste(missing_result_columns, collapse = ", "),
    call. = FALSE
  )
}

sex_results <- sex_results |>
  dplyr::mutate(
    cluster_id = as.character(.data$cluster_id),
    ensembl_gene_id = as.character(.data$ensembl_gene_id),
    gene = as.character(.data$gene),
    logFC = as.numeric(.data$logFC),
    PValue = as.numeric(.data$PValue),
    FDR = as.numeric(.data$FDR)
  ) |>
  dplyr::filter(
    !is.na(.data$cluster_id),
    .data$cluster_id != "",
    !is.na(.data$ensembl_gene_id),
    .data$ensembl_gene_id != ""
  )

cluster_ids <- sort_cluster_ids(
  sex_results$cluster_id
)

unknown_cluster_ids <- setdiff(
  cluster_ids,
  names(custom_cluster_labels)
)

if (length(unknown_cluster_ids) > 0L) {
  stop(
    "Missing anatomical annotation for cluster(s): ",
    paste(unknown_cluster_ids, collapse = ", "),
    call. = FALSE
  )
}

all_gene_ids <- sex_results |>
  dplyr::filter(
    is.finite(.data$logFC)
  ) |>
  dplyr::distinct(.data$ensembl_gene_id) |>
  dplyr::pull(.data$ensembl_gene_id)

significant_cluster_gene_results <- sex_results |>
  dplyr::filter(
    !is.na(.data$FDR),
    is.finite(.data$logFC),
    .data$FDR < fdr_threshold,
    abs(.data$logFC) >= abs_log2fc_threshold
  )

significant_gene_ids <- significant_cluster_gene_results |>
  dplyr::distinct(.data$ensembl_gene_id) |>
  dplyr::pull(.data$ensembl_gene_id)

if (length(all_gene_ids) == 0L) {
  stop(
    "No finite sex-effect logFC values were found.",
    call. = FALSE
  )
}

if (length(significant_gene_ids) == 0L) {
  stop(
    "No genes satisfy FDR < ",
    fdr_threshold,
    " and |log2FC| >= ",
    abs_log2fc_threshold,
    " in any cluster.",
    call. = FALSE
  )
}

message(
  "\nClusters: ",
  paste(
    paste0("C", cluster_ids),
    collapse = ", "
  )
)

message(
  "ALL gene set: ",
  length(all_gene_ids),
  " unique genes with finite logFC in at least one cluster"
)

message(
  "SIGNIFICANT gene union: ",
  length(significant_gene_ids),
  " unique genes from ",
  nrow(significant_cluster_gene_results),
  " significant cluster-gene results"
)


# ==============================================================================
# 9. Generate four correlation heatmaps and XLSX files
# ==============================================================================

output_summary <- vector(
  mode = "list",
  length = nrow(correlation_variants)
)

for (variant_index in seq_len(nrow(correlation_variants))) {

  variant_current <- correlation_variants[
    variant_index,
    ,
    drop = FALSE
  ]

  file_index <- variant_current$file_index[[1]]
  gene_set <- variant_current$gene_set[[1]]
  correlation_method <- variant_current$method[[1]]
  variant_id <- variant_current$variant_id[[1]]

  message("\n")
  message(
    paste(
      rep("=", 78),
      collapse = ""
    )
  )
  message(
    "Generating: ",
    variant_id
  )

  if (identical(gene_set, "all")) {

    selected_gene_ids <- all_gene_ids

    gene_set_label <- "ALL tested genes"

    gene_set_definition <- paste0(
      "All genes with a finite cluster-specific edgeR log2FC in at least one cluster; ",
      "pairwise complete observations are used for each cluster pair"
    )

    output_gene_set_id <- "allTestedGenes"

    selection_note <- paste0(
      "No DE-significance filter is used for gene selection"
    )

  } else if (identical(gene_set, "significant")) {

    selected_gene_ids <- significant_gene_ids

    gene_set_label <- "SIGNIFICANT gene union"

    gene_set_definition <- paste0(
      "Union of unique genes significant in at least one cluster at FDR < ",
      fdr_threshold,
      " and |log2FC| >= ",
      abs_log2fc_threshold,
      "; for these genes, log2FC is then taken from every cluster in which the gene was tested"
    )

    output_gene_set_id <- paste0(
      "significantUnion_FDR005absLog2FC05"
    )

    selection_note <- paste0(
      "Gene selection is global across clusters and is not redefined separately for each cluster pair"
    )

  } else {
    stop(
      "Unknown gene set: ",
      gene_set,
      call. = FALSE
    )
  }

  logfc_matrix <- build_logfc_matrix(
    results_table = sex_results,
    selected_gene_ids = selected_gene_ids,
    cluster_ids = cluster_ids
  )

  correlation_results <- calculate_pairwise_correlations(
    logfc_matrix = logfc_matrix,
    cluster_ids = cluster_ids,
    correlation_method = correlation_method
  )

  correlation_matrix <- correlation_results$correlation_matrix
  pairwise_statistics <- correlation_results$pairwise_statistics

  number_of_clusters <- length(cluster_ids)
  number_of_cluster_pairs <- nrow(pairwise_statistics)
  number_of_selected_genes <- nrow(logfc_matrix)

  min_pairwise_genes <- min(
    pairwise_statistics$n_genes,
    na.rm = TRUE
  )

  max_pairwise_genes <- max(
    pairwise_statistics$n_genes,
    na.rm = TRUE
  )

  median_pairwise_genes <- stats::median(
    pairwise_statistics$n_genes,
    na.rm = TRUE
  )

  n_missing_correlations <- sum(
    !is.finite(pairwise_statistics$correlation)
  )

  correlation_method_label <- if (identical(correlation_method, "pearson")) {
    "Pearson"
  } else {
    "Spearman"
  }

  correlation_symbol_label <- if (identical(correlation_method, "pearson")) {
    "r"
  } else {
    "rho"
  }

  heatmap_title_lines <- c(
    paste0(
      "Sex main effect (Female vs Male): cluster-to-cluster correlation of log2FC profiles (",
      correlation_method_label,
      ")"
    ),

    paste0(
      "Test: ",
      selected_test_comparison,
      " | test_id = ",
      sex_test_id
    ),

    "Direction: positive log2FC = higher expression in Female | negative log2FC = higher expression in Male",

    paste0(
      "Gene set: ",
      gene_set_label,
      " | ",
      number_of_selected_genes,
      " unique genes in the union matrix"
    ),

    paste0(
      "Correlation input: cluster-specific edgeR log2FC",
      " | missing cluster-gene values retained as NA",
      " | no imputation"
    ),

    paste0(
      "Pairwise handling: complete observations separately for each cluster pair",
      " | n genes per pair = ",
      min_pairwise_genes,
      "-",
      max_pairwise_genes,
      " (median ",
      formatC(median_pairwise_genes, format = "f", digits = 0),
      ")"
    ),

    paste0(
      "Pairwise statistics: ",
      number_of_cluster_pairs,
      " unique cluster pairs",
      " | H0: correlation = 0",
      " | p-values adjusted by Benjamini-Hochberg within this heatmap"
    ),

    paste0(
      "Heatmap cells: ",
      correlation_symbol_label,
      " from -1 to +1",
      " | values printed to 2 decimals",
      " | exact correlation, n genes, p-value and FDR are provided in XLSX"
    ),

    paste0(
      "Cluster ordering: ",
      linkage_method,
      "-linkage hierarchical clustering using distance 1 - correlation"
    ),

    clustering_title
  )

  heatmap_title <- paste(
    heatmap_title_lines,
    collapse = "\n"
  )

  file_prefix <- paste0(
    sprintf("%02d", file_index),
    "_",
    dataset_name,
    "_effectSex_FemaleVsMale_",
    output_gene_set_id,
    "_",
    correlation_method,
    "_clusterCorrelation"
  )

  pdf_file <- file.path(
    correlation_output_dir,
    paste0(
      file_prefix,
      ".pdf"
    )
  )

  xlsx_file <- file.path(
    correlation_output_dir,
    paste0(
      file_prefix,
      ".xlsx"
    )
  )

  interpretation_note <- if (identical(gene_set, "significant")) {
    paste0(
      "Exploratory correlation inference: the SIGNIFICANT gene set was selected from the same DE results; ",
      "therefore correlation p-values/FDR should be treated descriptively rather than as independent confirmatory tests."
    )
  } else {
    paste0(
      "Exploratory correlation inference: genes are biologically dependent observations, so correlation p-values/FDR are reported ",
      "for completeness and should not be interpreted as fully independent confirmatory evidence."
    )
  }

  analysis_info <- tibble::tibble(
    parameter = c(
      "dataset",
      "effect",
      "test_id",
      "comparison",
      "log2FC_direction",
      "test_sheet",
      "gene_set",
      "gene_set_definition",
      "selection_note",
      "FDR_threshold_for_SIGNIFICANT",
      "abs_log2FC_threshold_for_SIGNIFICANT",
      "n_significant_cluster_gene_results",
      "n_significant_unique_genes",
      "n_selected_unique_genes",
      "n_clusters",
      "n_unique_cluster_pairs",
      "correlation_method",
      "missing_cluster_gene_values",
      "imputation",
      "pairwise_gene_handling",
      "min_n_genes_per_pair",
      "median_n_genes_per_pair",
      "max_n_genes_per_pair",
      "n_pairs_with_undefined_correlation",
      "correlation_test_null",
      "p_value_adjustment",
      "p_value_adjustment_scope",
      "heatmap_cluster_distance",
      "heatmap_linkage",
      "interpretation_note",
      "clustering_solution"
    ),
    value = c(
      dataset_name,
      "Sex main effect: Female vs Male",
      sex_test_id,
      selected_test_comparison,
      "positive = higher expression in Female; negative = higher expression in Male",
      selected_test_sheet,
      gene_set_label,
      gene_set_definition,
      selection_note,
      as.character(fdr_threshold),
      as.character(abs_log2fc_threshold),
      as.character(nrow(significant_cluster_gene_results)),
      as.character(length(significant_gene_ids)),
      as.character(number_of_selected_genes),
      as.character(number_of_clusters),
      as.character(number_of_cluster_pairs),
      correlation_method_label,
      "retained as NA",
      "none",
      "pairwise complete observations",
      as.character(min_pairwise_genes),
      as.character(median_pairwise_genes),
      as.character(max_pairwise_genes),
      as.character(n_missing_correlations),
      "correlation = 0",
      "Benjamini-Hochberg",
      paste0(
        "all finite p-values among the ",
        number_of_cluster_pairs,
        " unique cluster pairs for this heatmap"
      ),
      "1 - correlation",
      linkage_method,
      interpretation_note,
      clustering_title
    )
  )

  message(
    "Selected unique genes: ",
    number_of_selected_genes
  )

  message(
    "Pairwise n genes: min = ",
    min_pairwise_genes,
    ", median = ",
    median_pairwise_genes,
    ", max = ",
    max_pairwise_genes
  )

  message(
    "Writing PDF: ",
    pdf_file
  )

  write_correlation_heatmap(
    pdf_file = pdf_file,
    correlation_matrix = correlation_matrix,
    cluster_ids = cluster_ids,
    correlation_method = correlation_method,
    heatmap_title = heatmap_title
  )

  message(
    "Writing XLSX: ",
    xlsx_file
  )

  write_correlation_xlsx(
    xlsx_file = xlsx_file,
    correlation_matrix = correlation_matrix,
    pairwise_statistics = pairwise_statistics,
    analysis_info = analysis_info,
    cluster_ids = cluster_ids
  )

  output_summary[[variant_index]] <- tibble::tibble(
    file_index = file_index,
    variant_id = variant_id,
    gene_set = gene_set_label,
    correlation_method = correlation_method_label,
    n_selected_unique_genes = number_of_selected_genes,
    min_n_genes_per_pair = min_pairwise_genes,
    median_n_genes_per_pair = median_pairwise_genes,
    max_n_genes_per_pair = max_pairwise_genes,
    n_pairs = number_of_cluster_pairs,
    n_pairs_with_undefined_correlation = n_missing_correlations,
    pdf = pdf_file,
    xlsx = xlsx_file
  )
}


# ==============================================================================
# 10. Validate outputs and final report
# ==============================================================================

output_summary <- dplyr::bind_rows(
  output_summary
) |>
  dplyr::arrange(.data$file_index)

expected_output_files <- c(
  output_summary$pdf,
  output_summary$xlsx
)

missing_output_files <- expected_output_files[
  !file.exists(expected_output_files)
]

if (length(missing_output_files) > 0L) {
  stop(
    "Missing expected output file(s):\n",
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

message("\n")
message(
  paste(
    rep("=", 78),
    collapse = ""
  )
)
message("Cluster-correlation analysis completed successfully.")
message("\nOutput directory:")
message(normalizePath(correlation_output_dir, mustWork = TRUE))
message("\nGenerated outputs:")
print(output_summary, n = Inf)

# ==============================================================================
# End
# ==============================================================================
