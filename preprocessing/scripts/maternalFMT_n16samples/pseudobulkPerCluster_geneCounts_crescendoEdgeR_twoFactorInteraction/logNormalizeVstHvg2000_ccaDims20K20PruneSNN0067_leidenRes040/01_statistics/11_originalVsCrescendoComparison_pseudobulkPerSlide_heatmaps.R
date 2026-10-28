#!/usr/bin/env Rscript
# ==============================================================================
# 04_originalVsCrescendo_pseudobulkPerTissue_highDetection_logFC_heatmap.R
#
# INPUT
#   Outputs from:
#     03_originalVsCrescendoComparison_pseudobulkPerTissue.R
#
# SELECTION
#   1. Mean spot detection >= detection_threshold_percent
#      Default selection source = ORIGINAL.
#
#   2. Absolute donor-group log2FC > abs_log2fc_threshold
#      Default = gene passes in ORIGINAL OR CRESCENDO.
#
#   No p-value / FDR threshold is used.
#
# HEATMAP
#   ORIGINAL + CRESCENDO log2(count+1)+QN matrices are joined first,
#   then ONE row-wise Z-score is calculated across all 32 columns.
#
#   Heatmap colours:
#     blue -> white -> red
#
# RIGHT SIDE
#   Original mean detection: purple intensity + numeric %
#   Crescendo mean detection: purple intensity + numeric %
#   Gene symbol
#
# OUTPUT
#   One PDF heatmap.
# ==============================================================================

options(stringsAsFactors = FALSE)

SCRIPT_BUILD <- "20260903_V11_SIMPLE_STATISTICS_SUMMARY"

# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tibble",
  "readr",
  "openxlsx",
  "msigdbr",
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
  library(openxlsx)
  library(msigdbr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ==============================================================================
# 2. User settings
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"

dataset_name <- "maternalFMT_n16samples"

clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

# Main test used for gene selection.
edgeR_sheet <- "Overall_Group"

# Mean percentage of spots with count > 0.
detection_threshold_percent <- 50

# Which data version must satisfy the detection threshold?
#
#   "original"  = ORIGINAL only
#   "crescendo" = CRESCENDO only
#   "either"    = ORIGINAL OR CRESCENDO
#   "both"      = ORIGINAL AND CRESCENDO
#
# Recommended default for this exploratory comparison:
#   "either"
detection_selection_source <- "either"

# Effect-size threshold.
abs_log2fc_threshold <- 0.2

# Which data version must satisfy |log2FC| > threshold?
#
#   "original"  = ORIGINAL only
#   "crescendo" = CRESCENDO only
#   "either"    = ORIGINAL OR CRESCENDO
#   "both"      = ORIGINAL AND CRESCENDO
#
# Recommended default for this exploratory comparison:
#   "either"
log2fc_selection_source <- "either"

# Restrict displayed genes to the union of these two MSigDB Mouse sets.
housekeeping_msigdb_id <- "MM1338"
housekeeping_msigdb_name <- "HOUNKPE_HOUSEKEEPING_GENES"

ribosomal_msigdb_id <- "MM15911"
ribosomal_msigdb_name <- "WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS"

valid_selection_sources <- c(
  "original",
  "crescendo",
  "either",
  "both"
)

if (!detection_selection_source %in% valid_selection_sources) {
  stop(
    "`detection_selection_source` must be one of: ",
    paste(valid_selection_sources, collapse = ", "),
    call. = FALSE
  )
}

if (!log2fc_selection_source %in% valid_selection_sources) {
  stop(
    "`log2fc_selection_source` must be one of: ",
    paste(valid_selection_sources, collapse = ", "),
    call. = FALSE
  )
}

# One joint row-wise Z-score after joining ORIGINAL + CRESCENDO.
#
# IMPORTANT:
# This is the ONLY parameter controlling the Z-score display range.
#
# Example:
#   zscore_range <- 2  -> range -2 to +2
#   zscore_range <- 3  -> range -3 to +3
#   zscore_range <- 4  -> range -4 to +4
#
# The same value controls:
#   1. clipping of calculated row Z-scores,
#   2. heatmap colour limits,
#   3. heatmap legend limits and ticks.
zscore_range <- 4

if (
  length(zscore_range) != 1L ||
  !is.numeric(zscore_range) ||
  !is.finite(zscore_range) ||
  zscore_range <= 0
) {
  stop(
    "`zscore_range` must be one positive finite number.",
    call. = FALSE
  )
}

zscore_lower_limit <- -zscore_range
zscore_upper_limit <- zscore_range

zscore_legend_ticks <- seq(
  -zscore_range,
  zscore_range,
  length.out = 5
)

# Row clustering.
row_clustering_distance <- "2 - Pearson r"
row_clustering_linkage <- "complete"

# ==============================================================================
# 3. Paths
# ==============================================================================

comparison_root <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics",
  "03_originalVsCrescendoComparison_pseudobulkPerTissue"
)

original_root <- file.path(
  comparison_root,
  "01_original"
)

crescendo_root <- file.path(
  comparison_root,
  "02_crescendo"
)

original_qn_file <- file.path(
  original_root,
  "03_log2plus1_QN",
  paste0(
    dataset_name,
    "_original_pseudobulkPerTissue_log2plus1_QN.tsv.gz"
  )
)

crescendo_qn_file <- file.path(
  crescendo_root,
  "03_log2plus1_QN",
  paste0(
    dataset_name,
    "_crescendo_pseudobulkPerTissue_log2plus1_QN.tsv.gz"
  )
)

original_detection_file <- file.path(
  original_root,
  "04_detectionPercent",
  paste0(
    dataset_name,
    "_original_pseudobulkPerTissue_",
    "spotDetectionPercent_meanSD.tsv.gz"
  )
)

crescendo_detection_file <- file.path(
  crescendo_root,
  "04_detectionPercent",
  paste0(
    dataset_name,
    "_crescendo_pseudobulkPerTissue_",
    "spotDetectionPercent_meanSD.tsv.gz"
  )
)

original_edger_file <- file.path(
  original_root,
  "05_edgeR",
  paste0(
    dataset_name,
    "_original_pseudobulkPerTissue_",
    "edgeR_allGenes_nineTests.xlsx"
  )
)

crescendo_edger_file <- file.path(
  crescendo_root,
  "05_edgeR",
  paste0(
    dataset_name,
    "_crescendo_pseudobulkPerTissue_",
    "edgeR_allGenes_nineTests.xlsx"
  )
)

metadata_file <- file.path(
  project_root,
  "data",
  "metadata_autismFMT.tsv"
)

output_dir <- file.path(
  comparison_root,
  "03_heatmaps"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

heatmap_pdf <- file.path(
  output_dir,
  paste0(
    "01_",
    dataset_name,
    "_originalVsCrescendo_",
    "meanDetection",
    detection_threshold_percent,
    "_absLog2FC",
    abs_log2fc_threshold,
    "_jointZscore_heatmap.pdf"
  )
)

selected_genes_tsv <- file.path(
  output_dir,
  paste0(
    "02_",
    dataset_name,
    "_originalVsCrescendo_",
    "meanDetection",
    detection_threshold_percent,
    "_absLog2FC",
    abs_log2fc_threshold,
    "_selectedGenes.tsv"
  )
)

nominal_pvalue_summary_tsv <- file.path(
  output_dir,
  paste0(
    "03_",
    dataset_name,
    "_originalVsCrescendo_",
    "selectedGenes_simpleStatisticsSummary.tsv"
  )
)

# ==============================================================================
# 4. Validate inputs
# ==============================================================================

required_files <- c(
  original_qn_file,
  crescendo_qn_file,
  original_detection_file,
  crescendo_detection_file,
  original_edger_file,
  crescendo_edger_file,
  metadata_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

message("SCRIPT BUILD: ", SCRIPT_BUILD)

# ==============================================================================
# 5. Helper functions
# ==============================================================================

calculate_row_zscores <- function(
    expression_matrix,
    lower_limit = -2,
    upper_limit = 2
) {
  output_matrix <- matrix(
    0,
    nrow = nrow(expression_matrix),
    ncol = ncol(expression_matrix),
    dimnames = dimnames(expression_matrix)
  )

  for (row_index in seq_len(nrow(expression_matrix))) {
    values <- as.numeric(
      expression_matrix[row_index, ]
    )

    row_sd <- stats::sd(values)

    if (is.finite(row_sd) && row_sd > 0) {
      output_matrix[row_index, ] <-
        (values - mean(values)) / row_sd
    }
  }

  output_matrix[
    output_matrix < lower_limit
  ] <- lower_limit

  output_matrix[
    output_matrix > upper_limit
  ] <- upper_limit

  output_matrix
}

read_expression_matrix <- function(
    input_file,
    sample_order
) {
  x <- readr::read_tsv(
    input_file,
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    tibble::as_tibble()

  required_columns <- c(
    "ensembl_gene_id",
    "gene",
    "chromosome"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(x)
  )

  if (length(missing_columns) > 0L) {
    stop(
      basename(input_file),
      " is missing column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  missing_samples <- setdiff(
    sample_order,
    colnames(x)
  )

  if (length(missing_samples) > 0L) {
    stop(
      basename(input_file),
      " is missing sample(s): ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  expression_matrix <- as.matrix(
    x[
      ,
      sample_order,
      drop = FALSE
    ]
  )

  storage.mode(expression_matrix) <- "double"

  rownames(expression_matrix) <-
    as.character(
      x$ensembl_gene_id
    )

  list(
    table = x,
    matrix = expression_matrix
  )
}

make_selection_flag <- function(
    original_value,
    crescendo_value,
    threshold,
    source,
    absolute = FALSE
) {
  original_test <- if (absolute) {
    abs(original_value) > threshold
  } else {
    original_value >= threshold
  }

  crescendo_test <- if (absolute) {
    abs(crescendo_value) > threshold
  } else {
    crescendo_value >= threshold
  }

  original_test[
    is.na(original_test)
  ] <- FALSE

  crescendo_test[
    is.na(crescendo_test)
  ] <- FALSE

  if (source == "original") {
    return(original_test)
  }

  if (source == "crescendo") {
    return(crescendo_test)
  }

  if (source == "either") {
    return(
      original_test |
        crescendo_test
    )
  }

  if (source == "both") {
    return(
      original_test &
        crescendo_test
    )
  }

  stop(
    "Unknown selection source: ",
    source,
    call. = FALSE
  )
}

# ==============================================================================
# 6. Sample metadata and column order
# ==============================================================================

# IMPORTANT:
# metadata_autismFMT.tsv contains metadata for 20 samples, including samples that
# were excluded from the final n=16 analysis. Therefore the heatmap MUST NOT use
# all metadata rows. The analysis sample set is taken directly from the prepared
# ORIGINAL and CRESCENDO pseudobulk/QN inputs.

original_qn_header <- readr::read_tsv(
  original_qn_file,
  n_max = 0,
  show_col_types = FALSE,
  progress = FALSE
)

crescendo_qn_header <- readr::read_tsv(
  crescendo_qn_file,
  n_max = 0,
  show_col_types = FALSE,
  progress = FALSE
)

annotation_columns <- c(
  "ensembl_gene_id",
  "gene",
  "chromosome"
)

original_analysis_samples <- setdiff(
  colnames(original_qn_header),
  annotation_columns
)

crescendo_analysis_samples <- setdiff(
  colnames(crescendo_qn_header),
  annotation_columns
)

if (!setequal(
  original_analysis_samples,
  crescendo_analysis_samples
)) {
  stop(
    "ORIGINAL and CRESCENDO QN matrices contain different sample sets.",
    call. = FALSE
  )
}

if (length(original_analysis_samples) != 16L) {
  stop(
    "Expected 16 samples in prepared pseudobulk/QN inputs but found ",
    length(original_analysis_samples),
    ".",
    call. = FALSE
  )
}

metadata <- read.delim(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
) |>
  tibble::as_tibble() |>
  dplyr::filter(
    .data$sample_ID %in% original_analysis_samples
  ) |>
  dplyr::transmute(
    sample_ID =
      as.character(.data$sample_ID),
    slide_ID =
      as.character(.data$slide_ID),
    fmt_donor_group = factor(
      as.character(
        .data$fmt_donor_group
      ),
      levels = c(
        "Neurotypical",
        "ASD"
      )
    ),
    sex = factor(
      as.character(
        .data$sex
      ),
      levels = c(
        "Male",
        "Female"
      )
    )
  ) |>
  dplyr::distinct(
    .data$sample_ID,
    .keep_all = TRUE
  ) |>
  dplyr::arrange(
    .data$sex,
    .data$fmt_donor_group,
    .data$sample_ID
  ) |>
  dplyr::mutate(
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

sample_order <- metadata$sample_ID

if (length(sample_order) != 16L) {
  stop(
    "After filtering metadata to the prepared analysis inputs, expected 16 samples but found ",
    length(sample_order),
    ".",
    call. = FALSE
  )
}

if (!setequal(
  sample_order,
  original_analysis_samples
)) {
  stop(
    "Filtered metadata and pseudobulk/QN sample sets do not match.",
    call. = FALSE
  )
}

message(
  "Analysis samples taken from prepared pseudobulk/QN inputs: ",
  length(sample_order)
)

# ==============================================================================
# 7. Read detection statistics
# ==============================================================================

message("\nLoading detection statistics...")

original_detection <- readr::read_tsv(
  original_detection_file,
  show_col_types = FALSE,
  progress = FALSE
) |>
  tibble::as_tibble() |>
  dplyr::select(
    "ensembl_gene_id",
    "gene",
    "mean_percent_positive",
    "sd_percent_positive"
  ) |>
  dplyr::rename(
    original_mean_detection =
      "mean_percent_positive",
    original_sd_detection =
      "sd_percent_positive"
  )

crescendo_detection <- readr::read_tsv(
  crescendo_detection_file,
  show_col_types = FALSE,
  progress = FALSE
) |>
  tibble::as_tibble() |>
  dplyr::select(
    "ensembl_gene_id",
    "mean_percent_positive",
    "sd_percent_positive"
  ) |>
  dplyr::rename(
    crescendo_mean_detection =
      "mean_percent_positive",
    crescendo_sd_detection =
      "sd_percent_positive"
  )

detection_stats <- original_detection |>
  dplyr::left_join(
    crescendo_detection,
    by = "ensembl_gene_id"
  )

# ==============================================================================
# 8. Read edgeR results
# ==============================================================================

message("Loading edgeR results...")

original_edger <- openxlsx::read.xlsx(
  original_edger_file,
  sheet = edgeR_sheet
) |>
  tibble::as_tibble() |>
  dplyr::transmute(
    ensembl_gene_id =
      as.character(.data$ensembl_gene_id),
    original_log2FC =
      as.numeric(.data$logFC),
    original_PValue =
      as.numeric(.data$PValue),
    original_FDR =
      as.numeric(.data$FDR)
  )

crescendo_edger <- openxlsx::read.xlsx(
  crescendo_edger_file,
  sheet = edgeR_sheet
) |>
  tibble::as_tibble() |>
  dplyr::transmute(
    ensembl_gene_id =
      as.character(.data$ensembl_gene_id),
    crescendo_log2FC =
      as.numeric(.data$logFC),
    crescendo_PValue =
      as.numeric(.data$PValue),
    crescendo_FDR =
      as.numeric(.data$FDR)
  )

# ==============================================================================
# 9. Load MSigDB Mouse housekeeping + ribosomal gene sets
# ==============================================================================

message("Loading MSigDB Mouse gene sets...")

msigdb_mouse <- msigdbr::msigdbr(
  species = "Mus musculus"
) |>
  tibble::as_tibble()

required_msigdb_columns <- c(
  "gs_id",
  "gs_name",
  "gene_symbol"
)

missing_msigdb_columns <- setdiff(
  required_msigdb_columns,
  colnames(msigdb_mouse)
)

if (length(missing_msigdb_columns) > 0L) {
  stop(
    "msigdbr output is missing column(s): ",
    paste(missing_msigdb_columns, collapse = ", "),
    call. = FALSE
  )
}

housekeeping_genes <- msigdb_mouse |>
  dplyr::filter(
    .data$gs_id == housekeeping_msigdb_id |
      .data$gs_name == housekeeping_msigdb_name
  ) |>
  dplyr::pull(.data$gene_symbol) |>
  unique()

ribosomal_genes <- msigdb_mouse |>
  dplyr::filter(
    .data$gs_id == ribosomal_msigdb_id |
      .data$gs_name == ribosomal_msigdb_name
  ) |>
  dplyr::pull(.data$gene_symbol) |>
  unique()

if (length(housekeeping_genes) == 0L) {
  stop(
    "Housekeeping MSigDB Mouse set not found: ",
    housekeeping_msigdb_id,
    " / ",
    housekeeping_msigdb_name,
    call. = FALSE
  )
}

if (length(ribosomal_genes) == 0L) {
  stop(
    "Ribosomal MSigDB Mouse set not found: ",
    ribosomal_msigdb_id,
    " / ",
    ribosomal_msigdb_name,
    call. = FALSE
  )
}

message(
  "Housekeeping set ",
  housekeeping_msigdb_id,
  ": ",
  length(housekeeping_genes),
  " genes"
)

message(
  "Ribosomal set ",
  ribosomal_msigdb_id,
  ": ",
  length(ribosomal_genes),
  " genes"
)

# ==============================================================================
# 10. Select genes
# ==============================================================================

selection_table <- detection_stats |>
  dplyr::left_join(
    original_edger,
    by = "ensembl_gene_id"
  ) |>
  dplyr::left_join(
    crescendo_edger,
    by = "ensembl_gene_id"
  ) |>
  dplyr::mutate(
    is_housekeeping =
      .data$gene %in% housekeeping_genes,
    is_ribosomal =
      .data$gene %in% ribosomal_genes,
    gene_set = dplyr::case_when(
      .data$is_housekeeping & .data$is_ribosomal ~ "Housekeeping + Ribosomal",
      .data$is_housekeeping ~ "Housekeeping",
      .data$is_ribosomal ~ "Ribosomal",
      TRUE ~ "Other"
    ),
    pass_detection = make_selection_flag(
      original_value =
        .data$original_mean_detection,
      crescendo_value =
        .data$crescendo_mean_detection,
      threshold =
        detection_threshold_percent,
      source =
        detection_selection_source,
      absolute = FALSE
    ),
    pass_log2fc = make_selection_flag(
      original_value =
        .data$original_log2FC,
      crescendo_value =
        .data$crescendo_log2FC,
      threshold =
        abs_log2fc_threshold,
      source =
        log2fc_selection_source,
      absolute = TRUE
    )
  ) |>
  dplyr::filter(
    .data$pass_detection,
    .data$pass_log2fc,
    .data$is_housekeeping | .data$is_ribosomal
  ) |>
  dplyr::arrange(
    dplyr::desc(
      pmax(
        abs(.data$original_log2FC),
        abs(.data$crescendo_log2FC),
        na.rm = TRUE
      )
    ),
    dplyr::desc(
      .data$original_mean_detection
    )
  )

if (nrow(selection_table) == 0L) {
  stop(
    "No genes passed selection. ",
    "Detection threshold = ",
    detection_threshold_percent,
    "%; |log2FC| > ",
    abs_log2fc_threshold,
    ".",
    call. = FALSE
  )
}

message(
  "Selected genes: ",
  nrow(selection_table)
)

# ==============================================================================
# 11. Load ORIGINAL and CRESCENDO QN matrices
# ==============================================================================

message("Loading log2+1+QN expression matrices...")

original_expression <- read_expression_matrix(
  original_qn_file,
  sample_order
)

crescendo_expression <- read_expression_matrix(
  crescendo_qn_file,
  sample_order
)

selected_gene_ids <- selection_table$ensembl_gene_id

missing_original <- setdiff(
  selected_gene_ids,
  rownames(
    original_expression$matrix
  )
)

missing_crescendo <- setdiff(
  selected_gene_ids,
  rownames(
    crescendo_expression$matrix
  )
)

if (
  length(missing_original) > 0L ||
  length(missing_crescendo) > 0L
) {
  stop(
    "Selected genes are missing from one or both QN matrices.",
    call. = FALSE
  )
}

original_matrix <- original_expression$matrix[
  selected_gene_ids,
  sample_order,
  drop = FALSE
]

crescendo_matrix <- crescendo_expression$matrix[
  selected_gene_ids,
  sample_order,
  drop = FALSE
]

gene_labels <- make.unique(
  ifelse(
    is.na(selection_table$gene) |
      selection_table$gene == "",
    selection_table$ensembl_gene_id,
    selection_table$gene
  )
)

rownames(original_matrix) <- gene_labels
rownames(crescendo_matrix) <- gene_labels

selection_table$row_label <-
  gene_labels

# ==============================================================================
# 12. Combine first, then ONE row Z-score
# ==============================================================================

combined_matrix <- cbind(
  original_matrix,
  crescendo_matrix
)

colnames(combined_matrix) <- c(
  paste0(
    "Original__",
    sample_order
  ),
  paste0(
    "Crescendo__",
    sample_order
  )
)

zscore_matrix <- calculate_row_zscores(
  combined_matrix,
  lower_limit =
    zscore_lower_limit,
  upper_limit =
    zscore_upper_limit
)

# ==============================================================================
# 13. Row clustering
# ==============================================================================

if (nrow(zscore_matrix) > 1L) {
  correlation_matrix <- suppressWarnings(
    stats::cor(
      t(zscore_matrix),
      method = "pearson",
      use = "pairwise.complete.obs"
    )
  )

  correlation_matrix[
    !is.finite(correlation_matrix)
  ] <- 0

  diag(correlation_matrix) <- 1

  row_hclust <- stats::hclust(
    stats::as.dist(
      2 - correlation_matrix
    ),
    method =
      row_clustering_linkage
  )
} else {
  row_hclust <- FALSE
}

# ==============================================================================
# 14. Top annotations
# ==============================================================================

combined_column_metadata <- dplyr::bind_rows(
  metadata |>
    dplyr::mutate(
      Data = "Original"
    ),
  metadata |>
    dplyr::mutate(
      Data = "Crescendo"
    )
) |>
  dplyr::mutate(
    Data = factor(
      .data$Data,
      levels = c(
        "Original",
        "Crescendo"
      )
    ),
    dataset_group = factor(
      paste(
        as.character(.data$Data),
        as.character(.data$combined_group),
        sep = " | "
      ),
      levels = c(
        "Original | Male_Neurotypical",
        "Original | Male_ASD",
        "Original | Female_Neurotypical",
        "Original | Female_ASD",
        "Crescendo | Male_Neurotypical",
        "Crescendo | Male_ASD",
        "Crescendo | Female_Neurotypical",
        "Crescendo | Female_ASD"
      )
    )
  )

data_colors <- c(
  "Original" = "#7F7F7F",
  "Crescendo" = "#E69F00"
)

sex_colors <- c(
  "Male" = "#3C78D8",
  "Female" = "#C050A0"
)

group_colors <- c(
  "Neurotypical" = "grey70",
  "ASD" = "black"
)

slide_levels <- unique(
  as.character(
    metadata$slide_ID
  )
)

slide_palette <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#CC79A7",
  "#0072B2",
  "#D55E00",
  "#F0E442",
  "#999999"
)

slide_colors <- setNames(
  slide_palette[
    seq_along(
      slide_levels
    )
  ],
  slide_levels
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Data =
    combined_column_metadata$Data,
  Slide =
    combined_column_metadata$slide_ID,
  Group =
    combined_column_metadata$fmt_donor_group,
  Sex =
    combined_column_metadata$sex,
  col = list(
    Data = data_colors,
    Slide = slide_colors,
    Group = group_colors,
    Sex = sex_colors
  ),
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
  show_legend = FALSE
)

# ==============================================================================
# 15. Colours
# ==============================================================================

heatmap_color_function <- circlize::colorRamp2(
  c(
    -zscore_range,
    0,
    zscore_range
  ),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)

detection_color_function <- circlize::colorRamp2(
  c(
    0,
    50,
    100
  ),
  c(
    "#FFFFFF",
    "#C5B0D5",
    "#54278F"
  )
)

gene_set_colors <- c(
  "Housekeeping" = "#4D9221",
  "Ribosomal" = "#E6AB02",
  "Housekeeping + Ribosomal" = "#666666"
)

# ==============================================================================
# 16. Detection side heatmaps
# ==============================================================================

original_detection_matrix <- matrix(
  selection_table$original_mean_detection,
  ncol = 1,
  dimnames = list(
    gene_labels,
    "Original detection"
  )
)

crescendo_detection_matrix <- matrix(
  selection_table$crescendo_mean_detection,
  ncol = 1,
  dimnames = list(
    gene_labels,
    "Crescendo detection"
  )
)

gene_set_matrix <- matrix(
  selection_table$gene_set,
  ncol = 1,
  dimnames = list(
    gene_labels,
    "Gene set"
  )
)

# ==============================================================================
# 17. Main heatmap
# ==============================================================================

n_genes <- nrow(
  zscore_matrix
)

row_font_size <- if (
  n_genes <= 30L
) {
  9
} else if (
  n_genes <= 75L
) {
  7
} else if (
  n_genes <= 150L
) {
  5.5
} else {
  4.2
}

row_height_mm <- if (
  n_genes <= 30L
) {
  5.0
} else if (
  n_genes <= 75L
) {
  3.5
} else if (
  n_genes <= 150L
) {
  2.6
} else {
  2.0
}

heatmap_body_height_mm <-
  n_genes *
  row_height_mm

pdf_height_inches <- max(
  10,
  (
    heatmap_body_height_mm +
      105
  ) /
    25.4
)

pdf_width_inches <- 16

heatmap_title <- paste0(
  "Whole-tissue pseudobulk: ORIGINAL vs CRESCENDO\n",
  "Selection: mean spot detection >= ",
  detection_threshold_percent,
  "% [source = ",
  detection_selection_source,
  "]",
  " and |log2FC| > ",
  abs_log2fc_threshold,
  " [source = ",
  log2fc_selection_source,
  "]",
  " | edgeR test = ",
  edgeR_sheet,
  " | no p-value/FDR threshold",
  "\nGene-set filter: ",
  housekeeping_msigdb_id,
  " housekeeping OR ",
  ribosomal_msigdb_id,
  " ribosomal",
  "\nExpression: log2(count + 1) + QN | ",
  "ONE joint row Z-score across ORIGINAL + CRESCENDO | ",
  "display range = [-",
  zscore_range,
  ", +",
  zscore_range,
  "] | ",
  "n = ",
  n_genes,
  " genes"
)

main_heatmap <- ComplexHeatmap::Heatmap(
  zscore_matrix,
  name = "Joint row Z-score",
  col =
    heatmap_color_function,
  top_annotation =
    top_annotation,
  cluster_columns = FALSE,
  column_split =
    combined_column_metadata$dataset_group,
  cluster_column_slices =
    FALSE,
  column_gap =
    grid::unit(
      1.2,
      "mm"
    ),
  column_labels =
    combined_column_metadata$sample_ID,
  cluster_rows =
    row_hclust,
  row_dend_reorder =
    FALSE,
  row_dend_width =
    grid::unit(
      16,
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
      fontsize = 7,
      fontface = "bold"
    ),
  show_row_names =
    FALSE,
  rect_gp =
    grid::gpar(
      col = NA
    ),
  border =
    FALSE,
  use_raster =
    FALSE,
  height =
    grid::unit(
      heatmap_body_height_mm,
      "mm"
    ),
  width =
    grid::unit(
      ncol(zscore_matrix) * 3.8,
      "mm"
    ),
  column_title =
    heatmap_title,
  column_title_gp =
    grid::gpar(
      fontsize = 11,
      fontface = "bold",
      lineheight = 1.12
    ),
  show_heatmap_legend = FALSE
)

gene_set_heatmap <- ComplexHeatmap::Heatmap(
  gene_set_matrix,
  name = "Gene set",
  col = gene_set_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  column_names_side = "top",
  column_names_rot = 45,
  column_names_gp = grid::gpar(
    fontsize = 7,
    fontface = "bold"
  ),
  show_row_names = FALSE,
  width = grid::unit(
    4.5,
    "mm"
  ),
  rect_gp = grid::gpar(
    col = "white",
    lwd = 0.5
  ),
  show_heatmap_legend = FALSE
)

original_detection_heatmap <- ComplexHeatmap::Heatmap(
  original_detection_matrix,
  name = "Original detection %",
  col =
    detection_color_function,
  cluster_rows =
    FALSE,
  cluster_columns =
    FALSE,
  show_column_names =
    TRUE,
  column_names_side =
    "top",
  column_names_rot =
    45,
  column_names_gp =
    grid::gpar(
      fontsize = 7,
      fontface = "bold"
    ),
  show_row_names =
    FALSE,
  width =
    grid::unit(
      5.5,
      "mm"
    ),
  rect_gp =
    grid::gpar(
      col = "white",
      lwd = 0.5
    ),
  cell_fun = function(
      j,
      i,
      x,
      y,
      width,
      height,
      fill
  ) {
    grid::grid.text(
      sprintf(
        "%.1f",
        original_detection_matrix[
          i,
          j
        ]
      ),
      x = x,
      y = y,
      gp = grid::gpar(
        fontsize =
          row_font_size,
        col = if (
          original_detection_matrix[
            i,
            j
          ] >= 70
        ) {
          "white"
        } else {
          "black"
        }
      )
    )
  },
  show_heatmap_legend = FALSE
)

crescendo_detection_heatmap <- ComplexHeatmap::Heatmap(
  crescendo_detection_matrix,
  name = "Crescendo detection %",
  col =
    detection_color_function,
  cluster_rows =
    FALSE,
  cluster_columns =
    FALSE,
  show_column_names =
    TRUE,
  column_names_side =
    "top",
  column_names_rot =
    45,
  column_names_gp =
    grid::gpar(
      fontsize = 7,
      fontface = "bold"
    ),
  show_row_names =
    TRUE,
  row_labels =
    gene_labels,
  row_names_side =
    "right",
  row_names_gp =
    grid::gpar(
      fontsize =
        row_font_size
    ),
  row_names_max_width =
    grid::unit(
      35,
      "mm"
    ),
  width =
    grid::unit(
      5.5,
      "mm"
    ),
  rect_gp =
    grid::gpar(
      col = "white",
      lwd = 0.5
    ),
  cell_fun = function(
      j,
      i,
      x,
      y,
      width,
      height,
      fill
  ) {
    grid::grid.text(
      sprintf(
        "%.1f",
        crescendo_detection_matrix[
          i,
          j
        ]
      ),
      x = x,
      y = y,
      gp = grid::gpar(
        fontsize =
          row_font_size,
        col = if (
          crescendo_detection_matrix[
            i,
            j
          ] >= 70
        ) {
          "white"
        } else {
          "black"
        }
      )
    )
  },
  show_heatmap_legend = FALSE
)

heatmap_list <-
  main_heatmap +
  gene_set_heatmap +
  original_detection_heatmap +
  crescendo_detection_heatmap

# ==============================================================================
# 18. Save selected genes
# ==============================================================================

readr::write_tsv(
  selection_table |>
    dplyr::select(
      "ensembl_gene_id",
      "gene",
      "gene_set",
      "is_housekeeping",
      "is_ribosomal",
      "original_mean_detection",
      "crescendo_mean_detection",
      "original_sd_detection",
      "crescendo_sd_detection",
      "original_log2FC",
      "crescendo_log2FC",
      "original_PValue",
      "crescendo_PValue",
      "original_FDR",
      "crescendo_FDR"
    ),
  selected_genes_tsv,
  na = "NA"
)

# Very simple statistics summary for genes displayed in the heatmap.
#
# For each criterion report:
#   - how many genes pass in ORIGINAL,
#   - how many genes pass in CRESCENDO,
#   - how many pass in BOTH datasets.

original_pass_log2fc <-
  !is.na(selection_table$original_log2FC) &
  abs(selection_table$original_log2FC) > abs_log2fc_threshold

crescendo_pass_log2fc <-
  !is.na(selection_table$crescendo_log2FC) &
  abs(selection_table$crescendo_log2FC) > abs_log2fc_threshold

original_pass_p <-
  !is.na(selection_table$original_PValue) &
  selection_table$original_PValue < 0.05

crescendo_pass_p <-
  !is.na(selection_table$crescendo_PValue) &
  selection_table$crescendo_PValue < 0.05

original_pass_p_and_log2fc <-
  original_pass_p &
  original_pass_log2fc

crescendo_pass_p_and_log2fc <-
  crescendo_pass_p &
  crescendo_pass_log2fc

nominal_pvalue_summary <- tibble::tibble(
  criterion = c(
    paste0(
      "|log2FC| > ",
      abs_log2fc_threshold
    ),
    "Nominal p < 0.05",
    paste0(
      "Nominal p < 0.05 & |log2FC| > ",
      abs_log2fc_threshold
    )
  ),
  n_original = c(
    sum(
      original_pass_log2fc,
      na.rm = TRUE
    ),
    sum(
      original_pass_p,
      na.rm = TRUE
    ),
    sum(
      original_pass_p_and_log2fc,
      na.rm = TRUE
    )
  ),
  n_crescendo = c(
    sum(
      crescendo_pass_log2fc,
      na.rm = TRUE
    ),
    sum(
      crescendo_pass_p,
      na.rm = TRUE
    ),
    sum(
      crescendo_pass_p_and_log2fc,
      na.rm = TRUE
    )
  ),
  n_both_original_and_crescendo = c(
    sum(
      original_pass_log2fc &
        crescendo_pass_log2fc,
      na.rm = TRUE
    ),
    sum(
      original_pass_p &
        crescendo_pass_p,
      na.rm = TRUE
    ),
    sum(
      original_pass_p_and_log2fc &
        crescendo_pass_p_and_log2fc,
      na.rm = TRUE
    )
  )
)

readr::write_tsv(
  nominal_pvalue_summary,
  nominal_pvalue_summary_tsv,
  na = "NA"
)

# ==============================================================================
# 19. One compact vertical legend column
# ==============================================================================

legend_zscore <- ComplexHeatmap::Legend(
  title = "Joint row Z-score",
  col_fun = heatmap_color_function,
  at = zscore_legend_ticks,
  labels = format(
    zscore_legend_ticks,
    trim = TRUE,
    scientific = FALSE
  ),
  legend_height = grid::unit(30, "mm"),
  title_gp = grid::gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  labels_gp = grid::gpar(
    fontsize = 7
  )
)

legend_gene_set <- ComplexHeatmap::Legend(
  title = "Gene set",
  at = c(
    "Housekeeping",
    "Ribosomal",
    "Housekeeping + Ribosomal"
  ),
  labels = c(
    paste0("Housekeeping (", housekeeping_msigdb_id, ")"),
    paste0("Ribosomal (", ribosomal_msigdb_id, ")"),
    paste0(
      "Housekeeping + Ribosomal (",
      housekeeping_msigdb_id,
      " + ",
      ribosomal_msigdb_id,
      ")"
    )
  ),
  legend_gp = grid::gpar(
    fill = unname(
      gene_set_colors[
        c(
          "Housekeeping",
          "Ribosomal",
          "Housekeeping + Ribosomal"
        )
      ]
    )
  ),
  title_gp = grid::gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  labels_gp = grid::gpar(
    fontsize = 7
  )
)

legend_detection <- ComplexHeatmap::Legend(
  title = "Mean detection %",
  col_fun = detection_color_function,
  at = c(
    0,
    25,
    50,
    75,
    100
  ),
  legend_height = grid::unit(
    26,
    "mm"
  ),
  title_gp = grid::gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  labels_gp = grid::gpar(
    fontsize = 7
  )
)

legend_data <- ComplexHeatmap::Legend(
  title = "Data",
  at = names(data_colors),
  legend_gp = grid::gpar(
    fill = unname(data_colors)
  ),
  title_gp = grid::gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  labels_gp = grid::gpar(
    fontsize = 7
  )
)

legend_slide <- ComplexHeatmap::Legend(
  title = "Slide",
  at = names(slide_colors),
  legend_gp = grid::gpar(
    fill = unname(slide_colors)
  ),
  title_gp = grid::gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  labels_gp = grid::gpar(
    fontsize = 7
  )
)

legend_group <- ComplexHeatmap::Legend(
  title = "Group",
  at = names(group_colors),
  legend_gp = grid::gpar(
    fill = unname(group_colors)
  ),
  title_gp = grid::gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  labels_gp = grid::gpar(
    fontsize = 7
  )
)

legend_sex <- ComplexHeatmap::Legend(
  title = "Sex",
  at = names(sex_colors),
  legend_gp = grid::gpar(
    fill = unname(sex_colors)
  ),
  title_gp = grid::gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  labels_gp = grid::gpar(
    fontsize = 7
  )
)

packed_legend <- ComplexHeatmap::packLegend(
  legend_zscore,
  legend_gene_set,
  legend_detection,
  legend_data,
  legend_slide,
  legend_group,
  legend_sex,
  direction = "vertical",
  gap = grid::unit(
    2.0,
    "mm"
  )
)

# ==============================================================================
# 20. Write PDF
# ==============================================================================

message("\nWriting heatmap:")
message(heatmap_pdf)

grDevices::pdf(
  heatmap_pdf,
  width =
    pdf_width_inches,
  height =
    pdf_height_inches,
  onefile = TRUE,
  useDingbats = FALSE
)

ComplexHeatmap::draw(
  heatmap_list,
  heatmap_legend_side =
    "right",
  annotation_legend_side =
    "right",
  merge_legends =
    FALSE,
  heatmap_legend_list = list(
    packed_legend
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

# ==============================================================================
# 21. Final report
# ==============================================================================

message("")
message("============================================================")
message("FINISHED")
message("============================================================")

message(
  "Detection threshold: ",
  detection_threshold_percent,
  "%"
)

message(
  "Detection selection source: ",
  detection_selection_source
)

message(
  "|log2FC| threshold: ",
  abs_log2fc_threshold
)

message(
  "log2FC selection source: ",
  log2fc_selection_source
)

message(
  "Selected genes: ",
  n_genes
)

message(
  "Z-score display range: -",
  zscore_range,
  " to +",
  zscore_range
)

message("")
message("PDF:")
message(heatmap_pdf)

message("")
message("Selected genes TSV:")
message(selected_genes_tsv)

message("")
message("Simple statistics summary:")
print(nominal_pvalue_summary)

message("")
message("Simple statistics summary TSV:")
message(nominal_pvalue_summary_tsv)
