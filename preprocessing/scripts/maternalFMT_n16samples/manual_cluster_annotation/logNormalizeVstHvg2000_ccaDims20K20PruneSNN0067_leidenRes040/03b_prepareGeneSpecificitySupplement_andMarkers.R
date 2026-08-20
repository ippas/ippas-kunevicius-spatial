#!/usr/bin/env Rscript

# ==============================================================================
# 03b_prepareGeneSpecificitySupplement.R
#
# Purpose:
# Convert the output of 03_calculateGeneSpecificity_onFindMarkers.R into a
# readable Excel supplement containing:
# 1. a summary of marker counts,
# 2. one sheet per cluster,
# 3. a notes sheet describing the reported variables.
# ==============================================================================

# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
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
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(openxlsx)
})

# ==============================================================================
# 2. Configuration
# ==============================================================================

project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

dataset_name <- "maternalFMT_n16samples"

configuration_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

expected_number_of_clusters <- 16L

bonferroni_adjusted_p_value_threshold <- 0.05

tau_thresholds <- seq(
  from = 0.1,
  to = 0.9,
  by = 0.1
)

number_of_decimal_places <- 5L

# ==============================================================================
# 3. Paths
# ==============================================================================

input_directory <- file.path(
  project_root,
  "results",
  dataset_name,
  "manual_cluster_annotation",
  configuration_name,
  "02_geneSpecificity_onFindMarkers"
)

input_file <- file.path(
  input_directory,
  paste0(
    "01_",
    dataset_name,
    "_leidenRes040_",
    "findAllMarkers_withGeneSpecificity.tsv"
  )
)

output_file <- file.path(
  input_directory,
  paste0(
    "02_",
    dataset_name,
    "_leidenRes040_",
    "findAllMarkers_geneSpecificitySupplement.xlsx"
  )
)

dir.create(
  input_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==============================================================================
# 4. Helper functions
# ==============================================================================

sort_cluster_ids <- function(cluster_ids) {

  cluster_ids <- unique(
    as.character(cluster_ids)
  )

  numeric_ids <- suppressWarnings(
    as.numeric(cluster_ids)
  )

  if (!anyNA(numeric_ids)) {
    return(
      cluster_ids[
        order(numeric_ids)
      ]
    )
  }

  sort(cluster_ids)
}

make_decimal_format <- function(number_of_digits) {

  paste0(
    "0.",
    paste0(
      rep("0", number_of_digits),
      collapse = ""
    )
  )
}

# ==============================================================================
# 5. Read input
# ==============================================================================

if (!file.exists(input_file)) {
  stop(
    "Input file does not exist:\n",
    input_file,
    call. = FALSE
  )
}

markers <- read.delim(
  file = input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "p_val",
  "avg_log2FC",
  "pct.1",
  "pct.2",
  "p_val_adj",
  "cluster",
  "gene",
  "best_other_cluster",
  "mean_expression_best_other",
  "expression_specificity",
  "expression_ratio_vs_best_other",
  "is_best_cluster",
  "tau",
  "gini",
  "shannon_specificity"
)

missing_columns <- setdiff(
  required_columns,
  colnames(markers)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing required column(s):\n",
    paste(missing_columns, collapse = "\n"),
    call. = FALSE
  )
}

if (nrow(markers) == 0L) {
  stop(
    "The input marker table is empty.",
    call. = FALSE
  )
}

number_of_rows_before_filtering <- nrow(markers)

# ==============================================================================
# 6. Standardize column types
# ==============================================================================

numeric_columns <- c(
  "p_val",
  "avg_log2FC",
  "pct.1",
  "pct.2",
  "p_val_adj",
  "mean_expression_best_other",
  "expression_specificity",
  "expression_ratio_vs_best_other",
  "tau",
  "gini",
  "shannon_specificity"
)

for (column_name in numeric_columns) {
  markers[[column_name]] <- suppressWarnings(
    as.numeric(markers[[column_name]])
  )
}

markers$cluster <- as.character(
  markers$cluster
)

markers$gene <- as.character(
  markers$gene
)

markers$best_other_cluster <- as.character(
  markers$best_other_cluster
)

if (!is.logical(markers$is_best_cluster)) {
  markers$is_best_cluster <- toupper(
    as.character(markers$is_best_cluster)
  ) == "TRUE"
}

# ==============================================================================
# 7. Apply Bonferroni-adjusted P-value threshold
# ==============================================================================

markers <- markers[
  !is.na(markers$p_val_adj) &
    markers$p_val_adj <
      bonferroni_adjusted_p_value_threshold,
  ,
  drop = FALSE
]

if (nrow(markers) == 0L) {
  stop(
    "No markers remained after applying p_val_adj < ",
    bonferroni_adjusted_p_value_threshold,
    ".",
    call. = FALSE
  )
}

# ==============================================================================
# 8. Define unique markers
# ==============================================================================

marker_clusters_per_gene <- split(
  markers$cluster,
  markers$gene
)

number_of_marker_clusters_per_gene <- vapply(
  marker_clusters_per_gene,
  function(cluster_ids) {
    length(
      unique(cluster_ids)
    )
  },
  FUN.VALUE = integer(1)
)

markers$is_unique_marker <-
  number_of_marker_clusters_per_gene[markers$gene] == 1L

# ==============================================================================
# 9. Rename columns
# ==============================================================================

colnames(markers)[
  colnames(markers) == "pct.1"
] <- "pct_in_cluster"

colnames(markers)[
  colnames(markers) == "pct.2"
] <- "pct_outside_cluster"

colnames(markers)[
  colnames(markers) == "p_val"
] <- "p_value"

colnames(markers)[
  colnames(markers) == "p_val_adj"
] <- "p_value_adj_bonferroni"

# ==============================================================================
# 10. Define final marker-table columns
# ==============================================================================

final_marker_column_order <- c(
  "gene",
  "pct_in_cluster",
  "pct_outside_cluster",
  "avg_log2FC",
  "p_value",
  "p_value_adj_bonferroni",
  "tau",
  "gini",
  "shannon_specificity",
  "is_unique_marker",
  "best_other_cluster",
  "mean_expression_best_other",
  "expression_specificity",
  "expression_ratio_vs_best_other",
  "is_best_cluster"
)

missing_export_columns <- setdiff(
  final_marker_column_order,
  colnames(markers)
)

if (length(missing_export_columns) > 0L) {
  stop(
    "Missing export column(s):\n",
    paste(missing_export_columns, collapse = "\n"),
    call. = FALSE
  )
}

# ==============================================================================
# 11. Determine cluster order
# ==============================================================================

cluster_ids <- sort_cluster_ids(
  markers$cluster
)

if (
  !is.null(expected_number_of_clusters) &&
    length(cluster_ids) != expected_number_of_clusters
) {
  stop(
    "Expected ",
    expected_number_of_clusters,
    " clusters but found ",
    length(cluster_ids),
    ": ",
    paste(cluster_ids, collapse = ", "),
    call. = FALSE
  )
}

# ==============================================================================
# 12. Prepare cluster summary
# ==============================================================================

summary_rows <- vector(
  mode = "list",
  length = length(cluster_ids)
)

for (cluster_index in seq_along(cluster_ids)) {

  current_cluster <- cluster_ids[cluster_index]

  current_markers <- markers[
    markers$cluster == current_cluster,
    ,
    drop = FALSE
  ]

  current_summary <- data.frame(
    cluster = current_cluster,
    n_markers = nrow(current_markers),
    n_unique_markers = sum(
      current_markers$is_unique_marker,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  for (tau_threshold in tau_thresholds) {

    threshold_label <- sprintf(
      "%.1f",
      tau_threshold
    )

    all_column_name <- paste0(
      "n_markers_tau_ge_",
      threshold_label
    )

    unique_column_name <- paste0(
      "n_unique_markers_tau_ge_",
      threshold_label
    )

    passes_tau_threshold <-
      !is.na(current_markers$tau) &
      current_markers$tau >= tau_threshold

    current_summary[[all_column_name]] <- sum(
      passes_tau_threshold,
      na.rm = TRUE
    )

    current_summary[[unique_column_name]] <- sum(
      passes_tau_threshold &
        current_markers$is_unique_marker,
      na.rm = TRUE
    )
  }

  summary_rows[[cluster_index]] <- current_summary
}

cluster_summary <- do.call(
  what = rbind,
  args = summary_rows
)

rownames(cluster_summary) <- NULL

# ==============================================================================
# 13. Add Total genes, Mean per cluster, and SD per cluster
# ==============================================================================

numeric_summary_columns <- setdiff(
  colnames(cluster_summary),
  "cluster"
)

# The Total row must contain numbers of distinct genes across all clusters.
# It must not be calculated by summing cluster-level marker assignments,
# because the same gene can be a marker of more than one cluster.

total_row <- data.frame(
  cluster = "Total",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

total_row$n_markers <- length(
  unique(markers$gene)
)

total_row$n_unique_markers <- length(
  unique(
    markers$gene[
      markers$is_unique_marker
    ]
  )
)

for (tau_threshold in tau_thresholds) {

  threshold_label <- sprintf(
    "%.1f",
    tau_threshold
  )

  all_column_name <- paste0(
    "n_markers_tau_ge_",
    threshold_label
  )

  unique_column_name <- paste0(
    "n_unique_markers_tau_ge_",
    threshold_label
  )

  passes_tau_threshold <-
    !is.na(markers$tau) &
    markers$tau >= tau_threshold

  total_row[[all_column_name]] <- length(
    unique(
      markers$gene[
        passes_tau_threshold
      ]
    )
  )

  total_row[[unique_column_name]] <- length(
    unique(
      markers$gene[
        passes_tau_threshold &
          markers$is_unique_marker
      ]
    )
  )
}

# Mean and SD are calculated from the 16 cluster-level counts.

mean_row <- data.frame(
  cluster = "Mean_per_cluster",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sd_row <- data.frame(
  cluster = "SD_per_cluster",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

for (column_name in numeric_summary_columns) {

  current_values <- cluster_summary[[column_name]]

  mean_row[[column_name]] <- mean(
    current_values,
    na.rm = TRUE
  )

  sd_row[[column_name]] <- sd(
    current_values,
    na.rm = TRUE
  )
}

summary_statistics <- rbind(
  total_row,
  mean_row,
  sd_row
)

summary_table <- rbind(
  cluster_summary,
  summary_statistics
)

rownames(summary_table) <- NULL

# ==============================================================================
# 14. Create workbook
# ==============================================================================

workbook <- createWorkbook(
  creator = "Mateusz Zięba"
)

decimal_format <- make_decimal_format(
  number_of_decimal_places
)

title_style <- createStyle(
  fontSize = 14,
  textDecoration = "bold",
  halign = "left",
  valign = "center"
)

group_header_style <- createStyle(
  fontSize = 10,
  textDecoration = "bold",
  fontColour = "#FFFFFF",
  fgFill = "#335C67",
  halign = "center",
  valign = "center",
  border = "TopBottomLeftRight",
  borderColour = "#FFFFFF"
)

subheader_style <- createStyle(
  fontSize = 10,
  textDecoration = "bold",
  fontColour = "#FFFFFF",
  fgFill = "#527985",
  halign = "center",
  valign = "center",
  border = "TopBottomLeftRight",
  borderColour = "#FFFFFF"
)

column_header_style <- createStyle(
  fontSize = 10,
  textDecoration = "bold",
  fontColour = "#FFFFFF",
  fgFill = "#335C67",
  halign = "center",
  valign = "center",
  wrapText = TRUE,
  border = "TopBottomLeftRight",
  borderColour = "#FFFFFF"
)

alternating_row_style <- createStyle(
  fgFill = "#F7F7F7"
)

total_row_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#E09F3E",
  border = "TopBottom"
)

mean_sd_row_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#F6E4C8"
)

integer_style <- createStyle(
  numFmt = "0"
)

decimal_style <- createStyle(
  numFmt = decimal_format
)

p_value_style <- createStyle(
  numFmt = "0.000E+00"
)

proportion_style <- createStyle(
  numFmt = "0.000"
)

wrap_text_style <- createStyle(
  wrapText = TRUE,
  valign = "top"
)

# ==============================================================================
# 15. Create summary sheet
# ==============================================================================

summary_sheet <- "00_Summary"

addWorksheet(
  workbook,
  summary_sheet,
  gridLines = FALSE
)

writeData(
  workbook,
  summary_sheet,
  x = "Marker counts by cluster and tau specificity threshold",
  startRow = 1,
  startCol = 1,
  colNames = FALSE
)

mergeCells(
  workbook,
  summary_sheet,
  cols = 1:ncol(cluster_summary),
  rows = 1
)

addStyle(
  workbook,
  summary_sheet,
  style = title_style,
  rows = 1,
  cols = 1,
  gridExpand = FALSE
)

writeData(
  workbook,
  summary_sheet,
  x = paste0(
    "Markers retained at Bonferroni-adjusted P < ",
    bonferroni_adjusted_p_value_threshold,
    "."
  ),
  startRow = 2,
  startCol = 1,
  colNames = FALSE
)

mergeCells(
  workbook,
  summary_sheet,
  cols = 1:ncol(cluster_summary),
  rows = 2
)

summary_group_header_row <- 4L
summary_subheader_row <- 5L
summary_data_start_row <- 6L

mergeCells(
  workbook,
  summary_sheet,
  cols = 1,
  rows = summary_group_header_row:summary_subheader_row
)

writeData(
  workbook,
  summary_sheet,
  x = "Cluster",
  startRow = summary_group_header_row,
  startCol = 1,
  colNames = FALSE
)

mergeCells(
  workbook,
  summary_sheet,
  cols = 2:3,
  rows = summary_group_header_row
)

writeData(
  workbook,
  summary_sheet,
  x = "All markers",
  startRow = summary_group_header_row,
  startCol = 2,
  colNames = FALSE
)

writeData(
  workbook,
  summary_sheet,
  x = matrix(
    c("All", "Unique"),
    nrow = 1L
  ),
  startRow = summary_subheader_row,
  startCol = 2,
  colNames = FALSE,
  rowNames = FALSE
)

current_column <- 4L

for (tau_threshold in tau_thresholds) {

  mergeCells(
    workbook,
    summary_sheet,
    cols = current_column:(current_column + 1L),
    rows = summary_group_header_row
  )

  writeData(
    workbook,
    summary_sheet,
    x = paste0(
      "Tau >= ",
      sprintf("%.1f", tau_threshold)
    ),
    startRow = summary_group_header_row,
    startCol = current_column,
    colNames = FALSE
  )

  writeData(
    workbook,
    summary_sheet,
    x = matrix(
      c("All", "Unique"),
      nrow = 1L
    ),
    startRow = summary_subheader_row,
    startCol = current_column,
    colNames = FALSE,
    rowNames = FALSE
  )

  current_column <- current_column + 2L
}

addStyle(
  workbook,
  summary_sheet,
  style = group_header_style,
  rows = summary_group_header_row,
  cols = 1:ncol(cluster_summary),
  gridExpand = TRUE
)

addStyle(
  workbook,
  summary_sheet,
  style = subheader_style,
  rows = summary_subheader_row,
  cols = 2:ncol(cluster_summary),
  gridExpand = TRUE
)

writeData(
  workbook,
  summary_sheet,
  x = summary_table,
  startRow = summary_data_start_row,
  startCol = 1,
  colNames = FALSE,
  rowNames = FALSE
)

cluster_data_end_row <-
  summary_data_start_row +
  nrow(cluster_summary) -
  1L

total_row_number <- cluster_data_end_row + 1L
mean_row_number <- total_row_number + 1L
sd_row_number <- mean_row_number + 1L

alternating_rows <- seq(
  from = summary_data_start_row,
  to = cluster_data_end_row,
  by = 2L
)

addStyle(
  workbook,
  summary_sheet,
  style = alternating_row_style,
  rows = alternating_rows,
  cols = 1:ncol(cluster_summary),
  gridExpand = TRUE
)

addStyle(
  workbook,
  summary_sheet,
  style = total_row_style,
  rows = total_row_number,
  cols = 1:ncol(cluster_summary),
  gridExpand = TRUE
)

addStyle(
  workbook,
  summary_sheet,
  style = mean_sd_row_style,
  rows = mean_row_number:sd_row_number,
  cols = 1:ncol(cluster_summary),
  gridExpand = TRUE
)

addStyle(
  workbook,
  summary_sheet,
  style = integer_style,
  rows = summary_data_start_row:total_row_number,
  cols = 2:ncol(cluster_summary),
  gridExpand = TRUE,
  stack = TRUE
)

addStyle(
  workbook,
  summary_sheet,
  style = decimal_style,
  rows = mean_row_number:sd_row_number,
  cols = 2:ncol(cluster_summary),
  gridExpand = TRUE,
  stack = TRUE
)

setColWidths(
  workbook,
  summary_sheet,
  cols = 1,
  widths = 20
)

setColWidths(
  workbook,
  summary_sheet,
  cols = 2:ncol(cluster_summary),
  widths = 11
)

setRowHeights(
  workbook,
  summary_sheet,
  rows = summary_group_header_row:summary_subheader_row,
  heights = 28
)

freezePane(
  workbook,
  summary_sheet,
  firstActiveRow = summary_data_start_row,
  firstActiveCol = 2
)

# ==============================================================================
# 16. Create one marker sheet per cluster
# ==============================================================================

for (cluster_id in cluster_ids) {

  sheet_name <- paste0(
    "Cluster_",
    cluster_id
  )

  cluster_markers <- markers[
    markers$cluster == cluster_id,
    ,
    drop = FALSE
  ]

  cluster_markers <- cluster_markers[
    order(
      -cluster_markers$tau,
      -cluster_markers$avg_log2FC,
      cluster_markers$p_value_adj_bonferroni,
      cluster_markers$gene
    ),
    ,
    drop = FALSE
  ]

  cluster_markers <- cluster_markers[
    ,
    final_marker_column_order,
    drop = FALSE
  ]

  addWorksheet(
    workbook,
    sheet_name,
    gridLines = FALSE
  )

  writeData(
    workbook,
    sheet_name,
    x = cluster_markers,
    startRow = 1,
    startCol = 1,
    colNames = TRUE,
    rowNames = FALSE,
    withFilter = TRUE
  )

  addStyle(
    workbook,
    sheet_name,
    style = column_header_style,
    rows = 1,
    cols = 1:ncol(cluster_markers),
    gridExpand = TRUE
  )

  if (nrow(cluster_markers) > 0L) {

    data_rows <- 2:(nrow(cluster_markers) + 1L)

    proportion_columns <- match(
      c(
        "pct_in_cluster",
        "pct_outside_cluster"
      ),
      colnames(cluster_markers)
    )

    addStyle(
      workbook,
      sheet_name,
      style = proportion_style,
      rows = data_rows,
      cols = proportion_columns,
      gridExpand = TRUE
    )

    decimal_columns <- match(
      c(
        "avg_log2FC",
        "tau",
        "gini",
        "shannon_specificity",
        "mean_expression_best_other",
        "expression_specificity",
        "expression_ratio_vs_best_other"
      ),
      colnames(cluster_markers)
    )

    addStyle(
      workbook,
      sheet_name,
      style = decimal_style,
      rows = data_rows,
      cols = decimal_columns,
      gridExpand = TRUE
    )

    p_value_columns <- match(
      c(
        "p_value",
        "p_value_adj_bonferroni"
      ),
      colnames(cluster_markers)
    )

    addStyle(
      workbook,
      sheet_name,
      style = p_value_style,
      rows = data_rows,
      cols = p_value_columns,
      gridExpand = TRUE
    )

  }

  setColWidths(
    workbook,
    sheet_name,
    cols = 1,
    widths = 18
  )

  setColWidths(
    workbook,
    sheet_name,
    cols = 2:ncol(cluster_markers),
    widths = 16
  )

  wide_columns <- which(
    colnames(cluster_markers) %in%
      c(
        "p_value_adj_bonferroni",
        "shannon_specificity",
        "is_unique_marker",
        "best_other_cluster",
        "mean_expression_best_other",
        "expression_specificity",
        "expression_ratio_vs_best_other",
        "is_best_cluster"
      )
  )

  setColWidths(
    workbook,
    sheet_name,
    cols = wide_columns,
    widths = 22
  )

  setRowHeights(
    workbook,
    sheet_name,
    rows = 1,
    heights = 42
  )

  freezePane(
    workbook,
    sheet_name,
    firstActiveRow = 2,
    firstActiveCol = 2
  )
}

# ==============================================================================
# 17. Create notes sheet
# ==============================================================================

notes_sheet <- "99_Notes"

addWorksheet(
  workbook,
  notes_sheet,
  gridLines = FALSE
)

notes <- data.frame(
  item = c(
    "Input",
    "Marker test",
    "Marker direction",
    "Minimum detection fraction",
    "Minimum average log2 fold change",
    "Statistical threshold",
    "Unique marker",
    "Specificity calculation",
    "Tau",
    "Normalized Gini coefficient",
    "Shannon specificity",
    "best_other_cluster",
    "mean_expression_best_other",
    "expression_specificity",
    "expression_ratio_vs_best_other",
    "is_best_cluster",
    "Summary: All",
    "Summary: Unique",
    "Total",
    "Mean_per_cluster",
    "SD_per_cluster"
  ),
  description = c(
    paste0(
      "Output of 03_calculateGeneSpecificity_onFindMarkers.R. ",
      "Rows not satisfying the Bonferroni-adjusted P-value threshold ",
      "were removed before preparing the workbook."
    ),
    paste0(
      "Wilcoxon rank-sum test comparing spots assigned to the target ",
      "cluster with all remaining spots."
    ),
    "Only positively enriched marker genes were considered.",
    paste0(
      "Genes were tested when detected in at least 25% of spots ",
      "in either comparison group."
    ),
    "Average log2 fold change of at least 0.25.",
    paste0(
      "Bonferroni-adjusted P value below ",
      bonferroni_adjusted_p_value_threshold,
      "."
    ),
    paste0(
      "A gene identified as a positive marker for only one cluster ",
      "among all clusters after statistical filtering."
    ),
    paste0(
      "Specificity was evaluated from mean expression calculated ",
      "separately for every cluster."
    ),
    "Higher values indicate greater cluster specificity.",
    "Higher values indicate greater cluster specificity.",
    "Higher values indicate greater cluster specificity.",
    "The other cluster with the highest mean expression of the gene.",
    "Mean expression of the gene in best_other_cluster.",
    paste0(
      "Fraction of the total mean expression across all clusters ",
      "assigned to the target cluster."
    ),
    paste0(
      "Mean expression in the target cluster divided by mean expression ",
      "in best_other_cluster."
    ),
    paste0(
      "TRUE when the target cluster has the highest mean expression ",
      "of the gene among all clusters."
    ),
    paste0(
      "Number of marker assignments satisfying the indicated tau ",
      "threshold. A gene may be counted for more than one cluster."
    ),
    paste0(
      "Number of markers satisfying the indicated tau threshold that ",
      "were identified as markers of only one cluster."
    ),
    paste0(
      "Number of distinct genes across all clusters. ",
      "The same gene is counted once even when it is a marker of multiple clusters."
    ),
    "Mean marker count across the 16 clusters.",
    "Standard deviation of marker counts across the 16 clusters."
  ),
  stringsAsFactors = FALSE
)

writeData(
  workbook,
  notes_sheet,
  x = "Marker identification and cluster-specificity supplement",
  startRow = 1,
  startCol = 1,
  colNames = FALSE
)

mergeCells(
  workbook,
  notes_sheet,
  cols = 1:2,
  rows = 1
)

addStyle(
  workbook,
  notes_sheet,
  style = title_style,
  rows = 1,
  cols = 1,
  gridExpand = FALSE
)

writeData(
  workbook,
  notes_sheet,
  x = notes,
  startRow = 3,
  startCol = 1,
  colNames = TRUE,
  rowNames = FALSE,
  withFilter = TRUE
)

addStyle(
  workbook,
  notes_sheet,
  style = column_header_style,
  rows = 3,
  cols = 1:2,
  gridExpand = TRUE
)

addStyle(
  workbook,
  notes_sheet,
  style = wrap_text_style,
  rows = 4:(nrow(notes) + 3L),
  cols = 1:2,
  gridExpand = TRUE
)

setColWidths(
  workbook,
  notes_sheet,
  cols = 1,
  widths = 34
)

setColWidths(
  workbook,
  notes_sheet,
  cols = 2,
  widths = 100
)

freezePane(
  workbook,
  notes_sheet,
  firstActiveRow = 4
)

# ==============================================================================
# 18. Save workbook
# ==============================================================================

saveWorkbook(
  workbook,
  output_file,
  overwrite = TRUE
)

# ==============================================================================
# 19. Console summary
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "Gene-specificity supplement completed\n"
)

cat(
  "============================================================\n"
)

cat(
  "Input rows before filtering: ",
  number_of_rows_before_filtering,
  "\n",
  sep = ""
)

cat(
  "Rows after Bonferroni-adjusted P < ",
  bonferroni_adjusted_p_value_threshold,
  ": ",
  nrow(markers),
  "\n",
  sep = ""
)

cat(
  "Clusters: ",
  length(cluster_ids),
  "\n",
  sep = ""
)

cat(
  "Unique marker genes: ",
  length(
    unique(
      markers$gene[
        markers$is_unique_marker
      ]
    )
  ),
  "\n",
  sep = ""
)

cat(
  "\nWorkbook saved to:\n",
  output_file,
  "\n",
  sep = ""
)
