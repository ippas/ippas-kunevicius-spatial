#!/usr/bin/env Rscript
# ==============================================================================
# 06_highDetectionGenes_unionFilters_originalPlusCrescendo_longHeatmap.R
#
# Requires object `gene_detection_summary` already present in the R session.
#
# ROW SELECTION = UNION of two criteria, evaluated per cluster-gene result:
#
#   mean detection > 80%
#
#   AND
#
#   FILTER A:
#     (|original log2FC| > 0.2 & original p-value < 0.05)
#     OR
#     (|Crescendo log2FC| > 0.2 & Crescendo p-value < 0.05)
#
#   OR FILTER B:
#     |original log2FC| > 0.5
#     OR
#     |Crescendo log2FC| > 0.5
#
# FILTER B has NO p-value requirement.
#
# HEATMAP:
#   - all selected cluster-gene rows from C1-C16 in ONE long heatmap;
#   - ORIGINAL and CRESCENDO samples are shown together;
#   - log2(count + 1) + QN values are combined first;
#   - ONE row-wise Z-score is then calculated across all 32 columns;
#   - Z-score is clipped to [-2, 2];
#   - blue -> white -> red;
#   - top annotations: Data, Slide, Group, Sex;
#   - row annotation: Cluster only;
#   - row clustering uses 2 - Pearson r and promotes rows from the same cluster.
# ==============================================================================

options(stringsAsFactors = FALSE)

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
    "Missing packages: ",
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
# 1. Settings
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"
clustering_name <- "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

test_id <- "Overall_Group_ASD_vs_Neurotypical"
cluster_ids <- as.character(1:16)

mean_detection_threshold <- 70
log2fc_threshold_filterA <- 0.3
pvalue_threshold_filterA <- 1
log2fc_threshold_filterB <- 0.3
pvalue_threshold_filterB <- 1

zscore_lower_limit <- -4
zscore_upper_limit <- 4

row_clustering_base_distance <- "2 - r"
row_clustering_same_cluster_bonus <- 0.10
row_clustering_linkage_method <- "complete"

if (!exists("gene_detection_summary")) {
  stop(
    "Object `gene_detection_summary` is not present. Run the detection script first.",
    call. = FALSE
  )
}

# ==============================================================================
# 2. Cluster colours and labels
# ==============================================================================

cluster_colors <- c(
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

cluster_labels <- c(
  "1"  = "posterior & sensory relay thalamic nuclei",
  "2"  = "isocortex, layers 4 & 5",
  "3"  = "cortical layers 1 & hippocampal neuropil",
  "4"  = "fiber tracts",
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
# 3. Input/output paths
# ==============================================================================

old_statistics_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)

crescendo_statistics_root <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)

crescendo_edger_dir <- file.path(
  crescendo_statistics_root,
  "02_edgeR"
)

old_edger_rdata <- file.path(
  old_statistics_dir,
  "03_edgeRResults",
  paste0(
    dataset_name,
    "_",
    clustering_name,
    "_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction_edgeRResults.RData"
  )
)

crescendo_edger_rdata <- file.path(
  crescendo_edger_dir,
  "03_edgeRResults",
  paste0(
    dataset_name,
    "_",
    clustering_name,
    "_crescendoCorrectedPseudobulk_edgeR_twoFactorInteraction_edgeRResults.RData"
  )
)

crescendo_input_rds <- file.path(
  crescendo_statistics_root,
  "01_crescendoBatchCorrection",
  "07_objects",
  paste0(
    dataset_name,
    "_",
    clustering_name,
    "_crescendoBatchCorrection_pseudobulkCounts_beforeEdgeR.rds"
  )
)

metadata_file <- file.path(
  project_root,
  "data",
  "metadata_autismFMT.tsv"
)

raw_qn_dir <- file.path(
  old_statistics_dir,
  "01_pseudobulkInputs",
  "perCluster_rawCounts_log2plus1Transformation_quantileNormalize"
)

crescendo_qn_dir <- file.path(
  crescendo_edger_dir,
  "01_pseudobulkInputs",
  "perCluster_crescendoCorrectedCounts_log2plus1_QN"
)

output_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,
  "02_heatmaps",
  "00_QC_residualBatchEffect",
  "03_highDetectionGenes_unionFilters_longHeatmap"
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
    "_highDetectionGenes_unionFilters_",
    "ORIGINALplusCRESCENDO_jointZscore_longHeatmap.pdf"
  )
)

# ==============================================================================
# 4. Helpers
# ==============================================================================

calculate_row_zscores <- function(
    expression_matrix,
    lower_limit = -2,
    upper_limit = 2
) {
  out <- matrix(
    0,
    nrow = nrow(expression_matrix),
    ncol = ncol(expression_matrix),
    dimnames = dimnames(expression_matrix)
  )

  for (i in seq_len(nrow(expression_matrix))) {
    x <- as.numeric(expression_matrix[i, ])
    s <- stats::sd(x)

    if (is.finite(s) && s > 0) {
      out[i, ] <- (x - mean(x)) / s
    }
  }

  out[out < lower_limit] <- lower_limit
  out[out > upper_limit] <- upper_limit

  out
}

modified_correlation_distance <- function(
    expression_matrix,
    cluster_labels_current,
    same_cluster_bonus = 0.20
) {
  cor_mat <- suppressWarnings(
    stats::cor(
      t(expression_matrix),
      method = "pearson",
      use = "pairwise.complete.obs"
    )
  )

  cor_mat[!is.finite(cor_mat)] <- 0
  diag(cor_mat) <- 1

  base_dist <- 2 - cor_mat

  same_cluster_matrix <- outer(
    cluster_labels_current,
    cluster_labels_current,
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

find_cluster_file <- function(
    input_dir,
    cluster_id,
    type = c("raw", "crescendo")
) {
  type <- match.arg(type)

  pattern <- if (type == "raw") {
    paste0(
      "_cluster_",
      cluster_id,
      "_rawCounts_log2plus1Transformation_quantileNormalize_",
      "proteinCodingAutosomalAndX\\.tsv\\.gz$"
    )
  } else {
    paste0(
      "_cluster_",
      cluster_id,
      "_crescendoCorrectedCounts_log2plus1_QN\\.tsv\\.gz$"
    )
  }

  files <- list.files(
    input_dir,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(files) != 1L) {
    stop(
      "Expected exactly one ",
      type,
      " file for cluster ",
      cluster_id,
      "; found ",
      length(files),
      ".",
      call. = FALSE
    )
  }

  files[[1]]
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

  matrix_out <- as.matrix(
    x[
      ,
      sample_order,
      drop = FALSE
    ]
  )

  storage.mode(matrix_out) <- "double"
  rownames(matrix_out) <- as.character(x$ensembl_gene_id)

  matrix_out
}

# ==============================================================================
# 5. Load edgeR results
# ==============================================================================

old_env <- new.env(parent = emptyenv())
crescendo_env <- new.env(parent = emptyenv())

load(old_edger_rdata, envir = old_env)
load(crescendo_edger_rdata, envir = crescendo_env)

old_results <- get(
  "edgeR_perCluster_combinedResults",
  envir = old_env
)[[test_id]] |>
  tibble::as_tibble() |>
  mutate(cluster_id = as.character(cluster_id))

crescendo_results <- get(
  "edgeR_perCluster_combinedResults",
  envir = crescendo_env
)[[test_id]] |>
  tibble::as_tibble() |>
  mutate(cluster_id = as.character(cluster_id))

# ==============================================================================
# 6. Sample metadata + order
# ==============================================================================

crescendo_object <- readRDS(crescendo_input_rds)

sample_metadata <- crescendo_object$sample_metadata |>
  tibble::as_tibble() |>
  transmute(
    sample_ID = as.character(sample_ID),
    fmt_donor_group = factor(
      as.character(fmt_donor_group),
      levels = c("Neurotypical", "ASD")
    ),
    sex = factor(
      as.character(sex),
      levels = c("Male", "Female")
    )
  ) |>
  distinct(sample_ID, .keep_all = TRUE)

slide_metadata <- read.delim(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
) |>
  tibble::as_tibble() |>
  transmute(
    sample_ID = as.character(sample_ID),
    slide_ID = as.character(slide_ID)
  ) |>
  distinct(sample_ID, .keep_all = TRUE)

sample_metadata <- sample_metadata |>
  left_join(
    slide_metadata,
    by = "sample_ID"
  ) |>
  arrange(
    sex,
    fmt_donor_group,
    sample_ID
  ) |>
  mutate(
    combined_group = factor(
      paste(
        as.character(sex),
        as.character(fmt_donor_group),
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

sample_order <- sample_metadata$sample_ID

if (length(sample_order) != 16L) {
  stop(
    "Expected 16 samples but found ",
    length(sample_order),
    ".",
    call. = FALSE
  )
}

# ==============================================================================
# 7. Select cluster-gene rows using UNION of filters
# ==============================================================================

detection_table <- gene_detection_summary |>
  select(
    ensembl_gene_id,
    gene,
    mean
  ) |>
  distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )

selected_rows <- vector(
  "list",
  length(cluster_ids)
)

names(selected_rows) <- cluster_ids

for (cluster_id_current in cluster_ids) {
  old_current <- old_results |>
    filter(
      cluster_id == cluster_id_current
    ) |>
    transmute(
      ensembl_gene_id,
      original_log2FC = logFC,
      original_pvalue = PValue,
      original_FDR = FDR
    )

  crescendo_current <- crescendo_results |>
    filter(
      cluster_id == cluster_id_current
    ) |>
    transmute(
      ensembl_gene_id,
      crescendo_log2FC = logFC,
      crescendo_pvalue = PValue,
      crescendo_FDR = FDR
    )

  selected_rows[[cluster_id_current]] <- detection_table |>
    left_join(
      old_current,
      by = "ensembl_gene_id"
    ) |>
    left_join(
      crescendo_current,
      by = "ensembl_gene_id"
    ) |>
    mutate(
      cluster_id = cluster_id_current,

      pass_filterA =
        (
          !is.na(original_log2FC) &
          !is.na(original_pvalue) &
          abs(original_log2FC) > log2fc_threshold_filterA &
          original_pvalue < pvalue_threshold_filterA
        ) |
        (
          !is.na(crescendo_log2FC) &
          !is.na(crescendo_pvalue) &
          abs(crescendo_log2FC) > log2fc_threshold_filterA &
          crescendo_pvalue < pvalue_threshold_filterA
        ),

      pass_filterB =
        (
          !is.na(original_log2FC) &
          abs(original_log2FC) > log2fc_threshold_filterB
        ) |
        (
          !is.na(crescendo_log2FC) &
          abs(crescendo_log2FC) > log2fc_threshold_filterB
        )
    ) |>
    filter(
      mean > mean_detection_threshold,
      pass_filterA | pass_filterB
    ) |>
    distinct(
      ensembl_gene_id,
      .keep_all = TRUE
    )
}

selected_results <- bind_rows(
  selected_rows
) |>
  mutate(
    cluster_id = as.character(cluster_id),
    row_id = paste0(
      ifelse(
        is.na(gene) | gene == "",
        ensembl_gene_id,
        gene
      ),
      " [C",
      cluster_id,
      "]"
    )
  ) |>
  arrange(
    as.integer(cluster_id),
    desc(mean),
    ensembl_gene_id
  )

if (nrow(selected_results) == 0L) {
  stop(
    "No cluster-gene rows passed the union filter.",
    call. = FALSE
  )
}

message("Selected cluster-gene rows: ", nrow(selected_results))
message(
  "Unique genes: ",
  dplyr::n_distinct(selected_results$ensembl_gene_id)
)

# ==============================================================================
# 8. Assemble ORIGINAL and CRESCENDO expression matrices
# ==============================================================================

original_matrix <- matrix(
  NA_real_,
  nrow = nrow(selected_results),
  ncol = length(sample_order),
  dimnames = list(
    selected_results$row_id,
    sample_order
  )
)

crescendo_matrix <- matrix(
  NA_real_,
  nrow = nrow(selected_results),
  ncol = length(sample_order),
  dimnames = list(
    selected_results$row_id,
    sample_order
  )
)

for (cluster_id_current in cluster_ids) {
  rows_current <- which(
    selected_results$cluster_id == cluster_id_current
  )

  if (length(rows_current) == 0L) {
    next
  }

  raw_file <- find_cluster_file(
    raw_qn_dir,
    cluster_id_current,
    type = "raw"
  )

  crescendo_file <- find_cluster_file(
    crescendo_qn_dir,
    cluster_id_current,
    type = "crescendo"
  )

  raw_all <- read_expression_matrix(
    raw_file,
    sample_order
  )

  crescendo_all <- read_expression_matrix(
    crescendo_file,
    sample_order
  )

  gene_ids_current <- selected_results$ensembl_gene_id[
    rows_current
  ]

  missing_raw <- setdiff(
    gene_ids_current,
    rownames(raw_all)
  )

  missing_crescendo <- setdiff(
    gene_ids_current,
    rownames(crescendo_all)
  )

  if (
    length(missing_raw) > 0L ||
    length(missing_crescendo) > 0L
  ) {
    stop(
      "Selected gene missing from expression matrix in cluster ",
      cluster_id_current,
      ".",
      call. = FALSE
    )
  }

  original_matrix[
    rows_current,
  ] <- raw_all[
    gene_ids_current,
    sample_order,
    drop = FALSE
  ]

  crescendo_matrix[
    rows_current,
  ] <- crescendo_all[
    gene_ids_current,
    sample_order,
    drop = FALSE
  ]
}

# ==============================================================================
# 9. Combine ORIGINAL + CRESCENDO, then ONE Z-score
# ==============================================================================

combined_matrix <- cbind(
  original_matrix,
  crescendo_matrix
)

colnames(combined_matrix) <- c(
  paste0("Original__", sample_order),
  paste0("Crescendo__", sample_order)
)

zscore_matrix <- calculate_row_zscores(
  combined_matrix,
  lower_limit = zscore_lower_limit,
  upper_limit = zscore_upper_limit
)

# ==============================================================================
# 10. Modified row clustering: promote same cluster
# ==============================================================================

if (nrow(zscore_matrix) > 1L) {
  modified_distance <- modified_correlation_distance(
    expression_matrix = zscore_matrix,
    cluster_labels_current = selected_results$cluster_id,
    same_cluster_bonus = row_clustering_same_cluster_bonus
  )

  row_hclust <- stats::hclust(
    modified_distance,
    method = row_clustering_linkage_method
  )
} else {
  row_hclust <- FALSE
}

# ==============================================================================
# 11. Top annotations
# ==============================================================================

combined_column_metadata <- bind_rows(
  sample_metadata |>
    mutate(Data = "Original"),
  sample_metadata |>
    mutate(Data = "Crescendo")
) |>
  mutate(
    Data = factor(
      Data,
      levels = c("Original", "Crescendo")
    ),
    dataset_group = factor(
      paste(
        as.character(Data),
        as.character(combined_group),
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

sex_colors <- c(
  "Male" = "#3C78D8",
  "Female" = "#C050A0"
)

group_colors <- c(
  "Neurotypical" = "grey70",
  "ASD" = "black"
)

data_colors <- c(
  "Original" = "#7F7F7F",
  "Crescendo" = "#E69F00"
)

slide_levels <- unique(
  as.character(sample_metadata$slide_ID)
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
  slide_palette[seq_along(slide_levels)],
  slide_levels
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Data = combined_column_metadata$Data,
  Slide = combined_column_metadata$slide_ID,
  Group = combined_column_metadata$fmt_donor_group,
  Sex = combined_column_metadata$sex,
  col = list(
    Data = data_colors,
    Slide = slide_colors,
    Group = group_colors,
    Sex = sex_colors
  ),
  simple_anno_size = grid::unit(4.2, "mm"),
  annotation_name_side = "left",
  annotation_name_gp = grid::gpar(
    fontsize = 9,
    fontface = "bold"
  )
)

# ==============================================================================
# 12. Cluster marker only on rows
# ==============================================================================

cluster_factor_levels <- as.character(
  sort(
    unique(
      as.integer(selected_results$cluster_id)
    )
  )
)

cluster_annotation_colors <- cluster_colors[
  cluster_factor_levels
]

names(cluster_annotation_colors) <- cluster_factor_levels

cluster_counts <- table(
  factor(
    selected_results$cluster_id,
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
    cluster_labels[
      cluster_factor_levels
    ]
  )
)

right_annotation <- ComplexHeatmap::rowAnnotation(
  Cluster = selected_results$cluster_id,
  col = list(
    Cluster = cluster_annotation_colors
  ),
  simple_anno_size = grid::unit(
    4.2,
    "mm"
  ),
  annotation_name_side = "top",
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
# 13. Heatmap colours
# ==============================================================================

heatmap_color_function <- circlize::colorRamp2(
  c(-2, 0, 2),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)

# ==============================================================================
# 14. Heatmap dimensions
# ==============================================================================

n_results <- nrow(selected_results)

row_font_size <- if (n_results <= 100L) {
  7
} else if (n_results <= 300L) {
  5.5
} else if (n_results <= 700L) {
  4.3
} else {
  3.5
}

row_height_mm <- if (n_results <= 100L) {
  3.1
} else if (n_results <= 300L) {
  2.35
} else if (n_results <= 700L) {
  1.90
} else {
  1.60
}

heatmap_body_height_mm <- n_results * row_height_mm

pdf_width_inches <- 20
pdf_height_inches <- max(
  14,
  (
    heatmap_body_height_mm + 110
  ) / 25.4
)

# ==============================================================================
# 15. Build heatmap
# ==============================================================================

heatmap_title <- paste0(
  "High-detection genes: ORIGINAL + CRESCENDO\n",
  "Selection = mean detection > ",
  mean_detection_threshold,
  "% AND [",
  "(|log2FC| > ",
  log2fc_threshold_filterA,
  " & p < ",
  pvalue_threshold_filterA,
  ") OR |log2FC| > ",
  log2fc_threshold_filterB,
  "] in ORIGINAL or CRESCENDO",
  "\nExpression = log2(count + 1) + QN | ONE joint row Z-score across 32 columns | clipped to [-2, 2]",
  "\nRow clustering = ",
  row_clustering_base_distance,
  " | same-cluster bonus = ",
  row_clustering_same_cluster_bonus,
  " | linkage = ",
  row_clustering_linkage_method,
  " | n = ",
  n_results,
  " cluster-gene rows"
)

heatmap_object <- ComplexHeatmap::Heatmap(
  zscore_matrix,
  name = "Joint row Z-score",
  col = heatmap_color_function,
  top_annotation = top_annotation,
  right_annotation = right_annotation,
  cluster_columns = FALSE,
  column_split = combined_column_metadata$dataset_group,
  cluster_column_slices = FALSE,
  column_gap = grid::unit(
    1.2,
    "mm"
  ),
  column_labels = combined_column_metadata$sample_ID,
  cluster_rows = row_hclust,
  row_dend_reorder = FALSE,
  row_dend_width = grid::unit(
    18,
    "mm"
  ),
  show_column_names = TRUE,
  column_names_side = "top",
  column_names_rot = 45,
  column_names_centered = TRUE,
  column_names_gp = grid::gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  show_row_names = TRUE,
  row_labels = ifelse(
    is.na(selected_results$gene) |
      selected_results$gene == "",
    selected_results$ensembl_gene_id,
    selected_results$gene
  ),
  row_names_side = "right",
  row_names_gp = grid::gpar(
    fontsize = row_font_size
  ),
  row_names_max_width = grid::unit(
    32,
    "mm"
  ),
  rect_gp = grid::gpar(
    col = NA
  ),
  border = FALSE,
  use_raster = FALSE,
  height = grid::unit(
    heatmap_body_height_mm,
    "mm"
  ),
  column_title = heatmap_title,
  column_title_side = "top",
  column_title_gp = grid::gpar(
    fontsize = 11,
    fontface = "bold",
    lineheight = 1.12
  ),
  heatmap_legend_param = list(
    title = "Joint row Z-score",
    at = c(-2, -1, 0, 1, 2),
    legend_height = grid::unit(
      45,
      "mm"
    )
  )
)

# ==============================================================================
# 16. Write PDF
# ==============================================================================

message("\nWriting PDF:")
message(heatmap_pdf)

grDevices::pdf(
  heatmap_pdf,
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
  padding = grid::unit(
    c(3, 3, 3, 3),
    "mm"
  )
)

grDevices::dev.off()

message("\n============================================================")
message("FINISHED")
message("============================================================")
message("Selected cluster-gene rows: ", n_results)
message("Unique genes: ", dplyr::n_distinct(selected_results$ensembl_gene_id))
message("PDF:")
message(heatmap_pdf)
