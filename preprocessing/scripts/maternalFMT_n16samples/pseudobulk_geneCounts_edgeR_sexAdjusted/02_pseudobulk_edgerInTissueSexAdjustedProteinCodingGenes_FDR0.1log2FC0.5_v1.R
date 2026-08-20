#!/usr/bin/env Rscript

# ==============================================================================
# 03_pseudobulk_edgeR_sexAdjusted_heatmap_FDR0.1_absLogFC0.5_rawLog2Quantile_FIXED.R
#
# PURPOSE
#   Create one heatmap for genes selected from the sex-adjusted edgeR analysis:
#
#     FDR < 0.1
#     abs(logFC) > 0.5
#
# INPUT
#   Files written by:
#     02_pseudobulk_inTissue_edgeR_sexAdjusted_proteinCodingGenes_Ensembl115.R
#
# EXPRESSION VALUES USED FOR VISUALIZATION
#   1. Read the saved raw pseudobulk counts for genes tested by edgeR.
#   2. Transform the complete matrix as log2(raw count + 1).
#   3. Apply quantile normalization across all 16 samples.
#   4. Select the genes meeting the edgeR criteria.
#   5. Calculate a row-wise Z-score and clip it to [-2, 2].
#
# HEATMAP LAYOUT
#   - all 16 retained samples are shown;
#   - columns are ordered by donor group and then sex:
#       Neurotypical Male
#       Neurotypical Female
#       ASD Male
#       ASD Female
#   - donor group and sex are shown as separate top annotations;
#   - genes upregulated in ASD are shown first;
#   - genes downregulated in ASD are shown second;
#   - genes are clustered separately within the two blocks;
#   - the up/down block labels are intentionally hidden;
#   - columns are not clustered, preserving the requested sample grouping;
#   - row-wise Z-scores are clipped to [-2, 2].
#
# OUTPUT
#   PDF heatmap plus selected-gene, sample-order and Z-score TSV files.
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

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\n\nInstall them with:\n",
    "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) ",
    "install.packages(\"BiocManager\")\n",
    "BiocManager::install(c(\"ComplexHeatmap\", \"limma\"))\n",
    "install.packages(c(\"circlize\", \"dplyr\", \"tibble\"))"
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
# 2. Define input and output paths
# ==============================================================================

project_dir <- "/home/mateusz/projects/ippas-kunevicius-spatial"

analysis_dir <- file.path(
  project_dir,
  "results",
  "maternalFMT_n16samples",
  "pseudobulk_geneCounts_edgeR_sexAdjusted"
)

full_results_file <- file.path(
  analysis_dir,
  "sexAdjusted_ASD_vs_Neurotypical_fullResults.tsv"
)

raw_counts_file <- file.path(
  analysis_dir,
  "sexAdjusted_rawCounts_testedGenes.tsv"
)

sample_metadata_file <- file.path(
  analysis_dir,
  "sexAdjusted_sampleMetadata.tsv"
)

heatmap_output_dir <- file.path(
  analysis_dir,
  "edgeR_selectedGenes_heatmap_FDR0.1_absLogFC0.5"
)

dir.create(
  heatmap_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

heatmap_pdf <- file.path(
  heatmap_output_dir,
  paste0(
    "sexAdjusted_ASD_vs_Neurotypical_",
    "FDR_lt0.1_absLog2FC_gt0.5_",
    "rawCounts_log2QuantileNormalized_",
    "rowZscore_heatmap.pdf"
  )
)

selected_genes_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "sexAdjusted_ASD_vs_Neurotypical_",
    "FDR_lt0.1_absLog2FC_gt0.5_",
    "selectedGenes.tsv"
  )
)

sample_order_tsv <- file.path(
  heatmap_output_dir,
  "sexAdjusted_heatmap_sampleOrder.tsv"
)

zscore_matrix_tsv <- file.path(
  heatmap_output_dir,
  paste0(
    "sexAdjusted_ASD_vs_Neurotypical_",
    "FDR_lt0.1_absLog2FC_gt0.5_",
    "rowZscoreMatrix.tsv"
  )
)


# ==============================================================================
# 3. Check input paths
# ==============================================================================

required_input_files <- c(
  full_results_file,
  raw_counts_file,
  sample_metadata_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    "Missing required input files:\n",
    paste(
      missing_input_files,
      collapse = "\n"
    ),
    "\n\nRun the sex-adjusted edgeR script first."
  )
}

if (!dir.exists(heatmap_output_dir)) {
  stop(
    "Could not create heatmap output directory: ",
    heatmap_output_dir
  )
}

message("Analysis directory: ", analysis_dir)
message("Full edgeR results: ", full_results_file)
message("Raw pseudobulk counts: ", raw_counts_file)
message("Sample metadata: ", sample_metadata_file)
message("Heatmap output directory: ", heatmap_output_dir)


# ==============================================================================
# 4. Read and validate edgeR results
# ==============================================================================

sexAdjusted_fullResults <- read.delim(
  full_results_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  tibble::as_tibble()

required_result_columns <- c(
  "ensembl_gene_id",
  "gene",
  "logFC",
  "logCPM",
  "F",
  "PValue",
  "FDR",
  "regulation"
)

missing_result_columns <- setdiff(
  required_result_columns,
  colnames(sexAdjusted_fullResults)
)

if (length(missing_result_columns) > 0) {
  stop(
    "Missing columns in the full edgeR results: ",
    paste(
      missing_result_columns,
      collapse = ", "
    )
  )
}

if (anyDuplicated(
  sexAdjusted_fullResults$ensembl_gene_id
)) {
  stop(
    "Duplicated Ensembl gene IDs were found ",
    "in the full edgeR results."
  )
}

selected_results <- sexAdjusted_fullResults |>
  filter(
    !is.na(FDR),
    !is.na(logFC),
    FDR < 0.1,
    abs(logFC) > 0.5
  ) |>
  mutate(
    regulation = if_else(
      logFC > 0,
      "up",
      "down"
    ),
    regulation = factor(
      regulation,
      levels = c("up", "down")
    )
  ) |>
  arrange(
    regulation,
    FDR,
    PValue,
    desc(abs(logFC))
  )

if (nrow(selected_results) == 0) {
  stop(
    "No genes meet FDR < 0.1 and abs(logFC) > 0.5."
  )
}

message(
  "\nSelected genes: ",
  nrow(selected_results)
)

message(
  "  Up in ASD: ",
  sum(selected_results$regulation == "up")
)

message(
  "  Down in ASD: ",
  sum(selected_results$regulation == "down")
)


# ==============================================================================
# 5. Read and validate sample metadata
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

if (length(missing_metadata_columns) > 0) {
  stop(
    "Missing columns in sample metadata: ",
    paste(
      missing_metadata_columns,
      collapse = ", "
    )
  )
}

if (anyDuplicated(sample_metadata$sample_ID)) {
  stop(
    "Duplicated sample_ID values were found in metadata."
  )
}

sample_metadata <- sample_metadata |>
  mutate(
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

if (anyNA(sample_metadata$fmt_donor_group)) {
  stop(
    "Unexpected or missing fmt_donor_group values."
  )
}

if (anyNA(sample_metadata$sex)) {
  stop(
    "Unexpected or missing sex values."
  )
}

sample_metadata <- sample_metadata |>
  arrange(
    fmt_donor_group,
    sex,
    sample_ID
  ) |>
  mutate(
    Group_Sex = factor(
      paste(
        fmt_donor_group,
        sex,
        sep = "_"
      ),
      levels = c(
        "Neurotypical_Male",
        "Neurotypical_Female",
        "ASD_Male",
        "ASD_Female"
      )
    )
  )

sample_order <- as.character(
  sample_metadata$sample_ID
)

message("\nHeatmap sample order:")
print(
  sample_metadata |>
    dplyr::select(
      sample_ID,
      fmt_donor_group,
      sex,
      Group_Sex
    ),
  n = Inf
)


# ==============================================================================
# 6. Read raw pseudobulk counts and preprocess them for visualization
#
# Requested visualization processing:
#   raw pseudobulk counts -> log2(count + 1) -> quantile normalization
#
# This matrix is used only for the heatmap. The edgeR differential-expression
# test itself remains based on the original raw counts with TMM normalization.
# ==============================================================================

raw_counts_table <- read.delim(
  raw_counts_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  tibble::as_tibble()

required_expression_columns <- c(
  "ensembl_gene_id",
  "gene"
)

missing_expression_columns <- setdiff(
  required_expression_columns,
  colnames(raw_counts_table)
)

if (length(missing_expression_columns) > 0) {
  stop(
    "Missing identifier columns in raw counts: ",
    paste(
      missing_expression_columns,
      collapse = ", "
    )
  )
}

if (anyDuplicated(
  raw_counts_table$ensembl_gene_id
)) {
  stop(
    "Duplicated Ensembl gene IDs were found ",
    "in the raw-count table."
  )
}

expression_sample_columns <- setdiff(
  colnames(raw_counts_table),
  required_expression_columns
)

missing_expression_samples <- setdiff(
  sample_order,
  expression_sample_columns
)

if (length(missing_expression_samples) > 0) {
  stop(
    "Metadata samples missing from raw counts: ",
    paste(
      missing_expression_samples,
      collapse = ", "
    )
  )
}

unexpected_expression_samples <- setdiff(
  expression_sample_columns,
  sample_order
)

if (length(unexpected_expression_samples) > 0) {
  stop(
    "Unexpected sample columns in raw counts: ",
    paste(
      unexpected_expression_samples,
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
    "Raw pseudobulk count matrix contains NA values."
  )
}

if (any(!is.finite(
  raw_count_matrix
))) {
  stop(
    "Raw pseudobulk count matrix contains ",
    "non-finite values."
  )
}

if (any(raw_count_matrix < 0)) {
  stop(
    "Raw pseudobulk count matrix contains negative values."
  )
}

if (!all(
  raw_count_matrix == floor(raw_count_matrix)
)) {
  stop(
    "Raw pseudobulk count matrix contains non-integer values."
  )
}

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
    "The log2 quantile-normalized matrix contains NA values."
  )
}

if (any(!is.finite(
  log2_quantile_normalized_matrix
))) {
  stop(
    "The log2 quantile-normalized matrix contains ",
    "non-finite values."
  )
}


# ==============================================================================
# 7. Select genes and calculate row-wise Z-scores
# ==============================================================================

selected_gene_ids <- as.character(
  selected_results$ensembl_gene_id
)

missing_selected_genes <- setdiff(
  selected_gene_ids,
  rownames(log2_quantile_normalized_matrix)
)

if (length(missing_selected_genes) > 0) {
  stop(
    "Selected genes missing from the log2 quantile-normalized matrix: ",
    paste(
      head(
        missing_selected_genes,
        20
      ),
      collapse = ", "
    ),
    if (
      length(missing_selected_genes) > 20
    ) {
      " ..."
    } else {
      ""
    }
  )
}

selected_expression_matrix <-
  log2_quantile_normalized_matrix[
    selected_gene_ids,
    sample_order,
    drop = FALSE
  ]

zscore_matrix <- t(
  scale(
    t(selected_expression_matrix),
    center = TRUE,
    scale = TRUE
  )
)

zero_variance_genes <- rownames(
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

if (length(zero_variance_genes) > 0) {
  warning(
    length(zero_variance_genes),
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


# ==============================================================================
# 8. Prepare row labels and up/down split
# ==============================================================================

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
    "Selected-gene annotation order could not be aligned ",
    "to the heatmap matrix."
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

row_direction <- factor(
  ifelse(
    selected_annotation$regulation == "up",
    "Up in ASD",
    "Down in ASD"
  ),
  levels = c(
    "Up in ASD",
    "Down in ASD"
  )
)


# ==============================================================================
# 9. Prepare column annotations
# ==============================================================================

column_annotation_data <- data.frame(
  Group = sample_metadata$fmt_donor_group,
  Sex = sample_metadata$sex,
  Group_Sex = sample_metadata$Group_Sex,
  row.names = sample_order,
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
  Group = column_annotation_data$Group,
  Sex = column_annotation_data$Sex,
  col = annotation_colors,
  simple_anno_size = grid::unit(
    4.5,
    "mm"
  ),
  annotation_name_side = "left",
  annotation_name_gp = grid::gpar(
    fontsize = 10,
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
# 10. Correlation distance for clustering genes within up/down blocks
# ==============================================================================

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
# 11. Define heatmap colors and dimensions
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

number_of_genes <- nrow(
  zscore_matrix_clipped
)

number_of_samples <- ncol(
  zscore_matrix_clipped
)

row_font_size <- if (
  number_of_genes <= 80
) {
  8
} else if (
  number_of_genes <= 130
) {
  6.8
} else if (
  number_of_genes <= 180
) {
  5.8
} else {
  5
}

row_height_mm <- if (
  number_of_genes <= 80
) {
  3.6
} else if (
  number_of_genes <= 130
) {
  2.9
} else if (
  number_of_genes <= 180
) {
  2.2
} else {
  2.1
}

heatmap_body_height_mm <-
  number_of_genes *
  row_height_mm

heatmap_body_width_mm <-
  number_of_samples *
  6.2

temporary_measurement_device <- FALSE

if (grDevices::dev.cur() == 1L) {
  grDevices::pdf(NULL)
  temporary_measurement_device <- TRUE
}

row_names_width_mm <- grid::convertWidth(
  ComplexHeatmap::max_text_width(
    rownames(
      zscore_matrix_clipped
    ),
    gp = grid::gpar(
      fontsize = row_font_size
    )
  ),
  unitTo = "mm",
  valueOnly = TRUE
)

title_lines_for_measurement <- c(
  "Sex-adjusted ASD vs Neurotypical",
  "Selected genes: FDR < 0.1 and |log2FC| > 0.5",
  "Input: raw pseudobulk counts -> log2(count + 1) -> quantile normalization",
  paste0(
    "Row Z-score clipped to [-2, 2] | ",
    number_of_genes,
    " genes | ",
    number_of_samples,
    " samples"
  )
)

title_width_mm <- grid::convertWidth(
  ComplexHeatmap::max_text_width(
    title_lines_for_measurement,
    gp = grid::gpar(
      fontsize = 12,
      fontface = "bold"
    )
  ),
  unitTo = "mm",
  valueOnly = TRUE
)

if (temporary_measurement_device) {
  grDevices::dev.off()
}

row_dendrogram_width_mm <- 19
legend_column_width_mm <- 28
internal_spacing_width_mm <- 12

plot_content_width_mm <-
  row_dendrogram_width_mm +
  heatmap_body_width_mm +
  row_names_width_mm +
  legend_column_width_mm +
  internal_spacing_width_mm

pdf_width_mm <- max(
  title_width_mm + 6,
  plot_content_width_mm
)

title_block_height_mm <- 24
column_names_height_mm <- 16
top_annotation_height_mm <- 10
outer_vertical_padding_mm <- 4

pdf_height_mm <-
  title_block_height_mm +
  column_names_height_mm +
  top_annotation_height_mm +
  heatmap_body_height_mm +
  outer_vertical_padding_mm

pdf_width_inches <- pdf_width_mm / 25.4
pdf_height_inches <- pdf_height_mm / 25.4


# ==============================================================================
# 12. Create the heatmap
# ==============================================================================

complete_title <- paste0(
  "Sex-adjusted ASD vs Neurotypical",
  "\nSelected genes: FDR < 0.1 and |log2FC| > 0.5",
  "\nInput: raw pseudobulk counts -> log2(count + 1) -> quantile normalization",
  "\nRow Z-score clipped to [-2, 2] | ",
  number_of_genes,
  " genes | ",
  number_of_samples,
  " samples"
)

heatmap_object <- ComplexHeatmap::Heatmap(
  matrix = zscore_matrix_clipped,
  name = "Row Z-score",
  col = heatmap_color_function,

  # Preserve the explicitly requested sample grouping.
  cluster_columns = FALSE,
  column_split = column_annotation_data$Group_Sex,
  cluster_column_slices = FALSE,
  column_gap = grid::unit(
    1.2,
    "mm"
  ),

  # Upregulated genes first, downregulated genes second.
  row_split = row_direction,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  clustering_distance_rows =
    row_correlation_distance,
  clustering_method_rows = "complete",
  row_dend_width = grid::unit(
    row_dendrogram_width_mm,
    "mm"
  ),
  row_gap = grid::unit(
    1.5,
    "mm"
  ),
  # Preserve the up-first/down-second slices without printing their labels.
  row_title = NULL,

  show_column_names = TRUE,
  column_names_side = "top",
  column_names_rot = 45,
  column_names_centered = TRUE,
  column_names_gp = grid::gpar(
    fontsize = 10,
    fontface = "bold"
  ),
  column_names_max_height =
    grid::unit(
      18,
      "mm"
    ),

  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(
    fontsize = row_font_size
  ),
  row_names_max_width =
    ComplexHeatmap::max_text_width(
      rownames(
        zscore_matrix_clipped
      ),
      gp = grid::gpar(
        fontsize = row_font_size
      )
    ) +
    grid::unit(
      3,
      "mm"
    ),

  top_annotation = top_annotation,

  width = grid::unit(
    heatmap_body_width_mm,
    "mm"
  ),
  height = grid::unit(
    heatmap_body_height_mm,
    "mm"
  ),

  rect_gp = grid::gpar(
    col = NA
  ),
  border = FALSE,
  use_raster = FALSE,

  column_title = complete_title,
  column_title_side = "top",
  column_title_gp = grid::gpar(
    fontsize = 12,
    fontface = "bold",
    lineheight = 1.15
  ),

  heatmap_legend_param = list(
    title = "Row Z-score",
    direction = "vertical",
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
    title_gp = grid::gpar(
      fontsize = 10,
      fontface = "bold"
    ),
    labels_gp = grid::gpar(
      fontsize = 9
    ),
    legend_height = grid::unit(
      48,
      "mm"
    )
  )
)


# ==============================================================================
# 13. Save PDF
# ==============================================================================

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
  padding = grid::unit(
    c(
      1,
      1,
      1,
      1
    ),
    "mm"
  )
)

grDevices::dev.off()

if (!file.exists(heatmap_pdf)) {
  stop(
    "Heatmap PDF was not created: ",
    heatmap_pdf
  )
}

if (file.info(heatmap_pdf)$size == 0) {
  stop(
    "Heatmap PDF is empty: ",
    heatmap_pdf
  )
}


# ==============================================================================
# 14. Save selected genes, sample order and Z-score matrix
# ==============================================================================

selected_results_output <- selected_results |>
  mutate(
    regulation = as.character(
      regulation
    )
  ) |>
  dplyr::select(
    ensembl_gene_id,
    gene,
    logFC,
    logCPM,
    F,
    PValue,
    FDR,
    regulation
  )

sample_order_output <- sample_metadata |>
  transmute(
    sample_order = row_number(),
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

zscore_matrix_output <- tibble::tibble(
  gene = rownames(
    zscore_matrix_clipped
  ),
  regulation = as.character(
    row_direction
  )
) |>
  bind_cols(
    tibble::as_tibble(
      as.data.frame(
        zscore_matrix_clipped,
        check.names = FALSE
      ),
      .name_repair = "minimal"
    )
  )

write.table(
  selected_results_output,
  file = selected_genes_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

write.table(
  sample_order_output,
  file = sample_order_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

write.table(
  zscore_matrix_output,
  file = zscore_matrix_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)


# ==============================================================================
# 15. Final summary
# ==============================================================================

heatmap_summary <- tibble::tibble(
  comparison =
    "ASD_vs_Neurotypical_adjusted_for_sex",
  selection =
    "FDR < 0.1 and abs(logFC) > 0.5",
  selected_genes =
    number_of_genes,
  up_in_ASD =
    sum(
      selected_results$regulation == "up"
    ),
  down_in_ASD =
    sum(
      selected_results$regulation == "down"
    ),
  samples =
    number_of_samples,
  sample_order =
    paste(
      levels(
        sample_metadata$Group_Sex
      ),
      collapse = " -> "
    ),
  expression =
    paste(
      "raw pseudobulk counts -> log2(count + 1)",
      "-> quantile normalization"
    ),
  visualization =
    "row Z-score clipped to [-2, 2]",
  output_pdf =
    normalizePath(
      heatmap_pdf,
      mustWork = TRUE
    )
)

print(
  heatmap_summary,
  n = Inf,
  width = Inf
)

message(
  "\nHeatmap PDF saved: ",
  normalizePath(
    heatmap_pdf,
    mustWork = TRUE
  )
)

message(
  "Selected genes saved: ",
  normalizePath(
    selected_genes_tsv,
    mustWork = TRUE
  )
)

message(
  "Sample order saved: ",
  normalizePath(
    sample_order_tsv,
    mustWork = TRUE
  )
)

message(
  "Z-score matrix saved: ",
  normalizePath(
    zscore_matrix_tsv,
    mustWork = TRUE
  )
)

message(
  "\nHeatmap analysis completed successfully."
)

# ==============================================================================
# End
# ==============================================================================
