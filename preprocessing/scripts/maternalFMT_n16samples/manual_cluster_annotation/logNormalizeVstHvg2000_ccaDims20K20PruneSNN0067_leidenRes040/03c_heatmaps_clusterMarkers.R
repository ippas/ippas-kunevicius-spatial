#!/usr/bin/env Rscript

# ==============================================================================
# 03c_heatmap_findMarkers_tau080_sampleByCluster.R
#
# Purpose:
# Create one heatmap for marker genes with tau >= 0.80.
#
# ROWS
#   Unique marker genes retained at:
#   - Wilcoxon rank-sum test,
#   - only positive markers,
#   - min.pct = 0.25,
#   - avg_log2FC >= 0.25,
#   - Bonferroni-adjusted P < 0.05,
#   - tau >= 0.80.
#
#   When a gene was identified as a marker of more than one cluster, one
#   representative marker cluster is assigned for row ordering. Preference is
#   given to a cluster in which the gene has the highest mean expression,
#   followed by expression specificity, average log2 fold change, and adjusted
#   P value.
#
# COLUMNS
#   Every biological sample is repeated once for each of the 16 clusters:
#
#     Cluster 1:  sample 1 ... sample 16
#     Cluster 2:  sample 1 ... sample 16
#     ...
#     Cluster 16: sample 1 ... sample 16
#
#   Sample IDs are not displayed because they are repeated 16 times. Columns
#   are annotated by source cluster, FMT donor group, and sex.
#
# HEATMAP VALUES
#   RNA data layer
#   -> back-transformation from log1p using expm1()
#   -> mean normalized expression for each sample x cluster combination
#   -> log2(mean expression + 1)
#   -> row-wise Z-score
#   -> values clipped to [-2, 2].
#
# CLUSTERING PARAMETERS
#   LogNormalize + VST
#   HVG = 2000
#   CCA integration
#   dimensions = 20
#   k = 20
#   prune.SNN = 0.0667
#   Leiden resolution = 0.40
#
# IMPORTANT
#   No custom functions are defined in this script. All steps are written
#   directly in the analysis workflow.
# ==============================================================================

options(stringsAsFactors = FALSE)

# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
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
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ==============================================================================
# 2. Analysis configuration
# ==============================================================================

project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

dataset_name <- "maternalFMT_n16samples"

configuration_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

cluster_column <- "leiden_res0.40"
sample_column <- "sample_ID"
sex_column <- "sex"
group_column <- "fmt_donor_group"
assay_name <- "RNA"

expected_number_of_samples <- 16L
expected_number_of_clusters <- 16L

bonferroni_adjusted_p_value_threshold <- 0.05
minimum_detection_fraction <- 0.25
minimum_average_log2_fold_change <- 0.25
tau_threshold <- 0.80

zscore_lower_limit <- -2
zscore_upper_limit <- 2

sample_order <- c(
  "1_1F",
  "1_1Fd",
  "1_1M",
  "2_1M",
  "2_1Md",
  "3_1F",
  "3_1M",
  "5_1M",
  "5_3F",
  "12_1M",
  "13_1F",
  "13_1M",
  "15_1F",
  "18_1F",
  "18_1M",
  "23_1F"
)

clustering_description <- paste(
  "LogNormalize + VST | HVG = 2000 | CCA dims = 20 |",
  "k = 20 | prune.SNN = 0.0667 | Leiden resolution = 0.40"
)

marker_description <- paste(
  "Positive markers: Wilcoxon rank-sum test | min.pct = 0.25 |",
  "avg_log2FC >= 0.25 | Bonferroni-adjusted P < 0.05 | tau >= 0.80"
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
  "1"  = "Thalamus 1",
  "2"  = "Cortex 1",
  "3"  = "Hippocampus 1",
  "4"  = "Fiber tracts",
  "5"  = "Striatum-like amygdalar nuclei",
  "6"  = "Hypothalamus",
  "7"  = "Cortex 3",
  "8"  = "Cortex 2",
  "9"  = "Thalamus 3",
  "10" = "Cortical subplate",
  "11" = "Thalamus 2",
  "12" = "Caudoputamen",
  "13" = "Hippocampus 2",
  "14" = "Vessels",
  "15" = "Ventricles",
  "16" = "Hippocampus 3"
)

sex_colors <- c(
  "Male" = "#3C78D8",
  "Female" = "#C050A0"
)

group_colors <- c(
  "Neurotypical" = "grey70",
  "ASD" = "black"
)

# ==============================================================================
# 4. Input and output paths
# ==============================================================================

input_rdata_file <- file.path(
  project_root,
  "results",
  dataset_name,
  "seurat_clustering_analysis",
  "logNormalize_vst",
  "hvg2000",
  "02_cca",
  "dims20",
  "k20",
  "prune0067",
  "RData",
  paste0(
    "02_",
    dataset_name,
    "_logNormalizeVst_hvg2000_",
    "ccaDims20K20PruneSNN0067_",
    "res010to100by010_multiClusteringAndUmap.RData"
  )
)

input_markers_file <- file.path(
  project_root,
  "results",
  dataset_name,
  "manual_cluster_annotation",
  configuration_name,
  "02_geneSpecificity_onFindMarkers",
  paste0(
    "01_",
    dataset_name,
    "_leidenRes040_",
    "findAllMarkers_withGeneSpecificity.tsv"
  )
)

output_directory <- file.path(
  project_root,
  "results",
  dataset_name,
  "manual_cluster_annotation",
  configuration_name,
  "03_geneSpecificityHeatmap_tau080"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

heatmap_pdf <- file.path(
  output_directory,
  paste0(
    "01_",
    dataset_name,
    "_leidenRes040_",
    "findAllMarkers_Bonferroni005_tauGE080_",
    "sampleByClusterMeanExpression_rowZscore_heatmap.pdf"
  )
)

selected_genes_tsv <- file.path(
  output_directory,
  paste0(
    "02_",
    dataset_name,
    "_leidenRes040_",
    "findAllMarkers_Bonferroni005_tauGE080_selectedGenes.tsv"
  )
)

mean_expression_matrix_tsv <- file.path(
  output_directory,
  paste0(
    "03_",
    dataset_name,
    "_leidenRes040_",
    "sampleByCluster_meanExpressionMatrix.tsv"
  )
)

log2_expression_matrix_tsv <- file.path(
  output_directory,
  paste0(
    "04_",
    dataset_name,
    "_leidenRes040_",
    "sampleByCluster_log2MeanExpressionMatrix.tsv"
  )
)

zscore_matrix_tsv <- file.path(
  output_directory,
  paste0(
    "05_",
    dataset_name,
    "_leidenRes040_",
    "sampleByCluster_rowZscoreMatrix.tsv"
  )
)

column_metadata_tsv <- file.path(
  output_directory,
  paste0(
    "06_",
    dataset_name,
    "_leidenRes040_",
    "sampleByCluster_columnMetadata.tsv"
  )
)

summary_tsv <- file.path(
  output_directory,
  paste0(
    "07_",
    dataset_name,
    "_leidenRes040_",
    "heatmapSummary.tsv"
  )
)

required_input_files <- c(
  input_rdata_file,
  input_markers_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0L) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_input_files, collapse = "\n"),
    call. = FALSE
  )
}

# ==============================================================================
# 5. Load Seurat object
# ==============================================================================

load_environment <- new.env()

loaded_object_names <- load(
  input_rdata_file,
  envir = load_environment
)

seurat_object_names <- loaded_object_names[
  vapply(
    loaded_object_names,
    function(object_name) {
      inherits(
        get(object_name, envir = load_environment),
        "Seurat"
      )
    },
    FUN.VALUE = logical(1)
  )
]

if (length(seurat_object_names) != 1L) {
  stop(
    "Expected exactly one Seurat object in the RData file.",
    call. = FALSE
  )
}

seurat_object <- get(
  seurat_object_names[[1]],
  envir = load_environment
)

rm(load_environment)

DefaultAssay(seurat_object) <- assay_name

data_layers <- Layers(
  seurat_object[[assay_name]],
  search = "^data"
)

if (length(data_layers) == 0L) {
  stop(
    "No RNA data layer was found.",
    call. = FALSE
  )
}

if (length(data_layers) > 1L) {
  seurat_object <- JoinLayers(
    object = seurat_object,
    assay = assay_name
  )
}

# ==============================================================================
# 6. Validate spot metadata
# ==============================================================================

spot_metadata <- seurat_object[[]]

required_metadata_columns <- c(
  sample_column,
  cluster_column,
  sex_column,
  group_column
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(spot_metadata)
)

if (length(missing_metadata_columns) > 0L) {
  stop(
    "Seurat metadata is missing required column(s): ",
    paste(missing_metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

spot_metadata$sample_ID_heatmap <- as.character(
  spot_metadata[[sample_column]]
)

spot_metadata$cluster_ID_heatmap <- as.character(
  spot_metadata[[cluster_column]]
)

spot_metadata$sex_heatmap <- as.character(
  spot_metadata[[sex_column]]
)

spot_metadata$group_heatmap <- as.character(
  spot_metadata[[group_column]]
)

observed_samples <- unique(
  spot_metadata$sample_ID_heatmap
)

missing_samples <- setdiff(
  sample_order,
  observed_samples
)

unexpected_samples <- setdiff(
  observed_samples,
  sample_order
)

if (length(missing_samples) > 0L) {
  stop(
    "The following expected samples are missing from the Seurat object: ",
    paste(missing_samples, collapse = ", "),
    call. = FALSE
  )
}

if (length(unexpected_samples) > 0L) {
  stop(
    "The Seurat object contains unexpected samples: ",
    paste(unexpected_samples, collapse = ", "),
    call. = FALSE
  )
}

if (length(sample_order) != expected_number_of_samples) {
  stop(
    "Expected ",
    expected_number_of_samples,
    " sample IDs in sample_order but found ",
    length(sample_order),
    ".",
    call. = FALSE
  )
}

cluster_ids <- unique(
  spot_metadata$cluster_ID_heatmap
)

numeric_cluster_ids <- suppressWarnings(
  as.numeric(cluster_ids)
)

if (!anyNA(numeric_cluster_ids)) {
  cluster_ids <- cluster_ids[
    order(numeric_cluster_ids)
  ]
} else {
  cluster_ids <- sort(cluster_ids)
}

if (length(cluster_ids) != expected_number_of_clusters) {
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

missing_cluster_colors <- setdiff(
  cluster_ids,
  names(custom_cluster_colors)
)

missing_cluster_labels <- setdiff(
  cluster_ids,
  names(custom_cluster_labels)
)

if (length(missing_cluster_colors) > 0L) {
  stop(
    "Missing colours for cluster(s): ",
    paste(missing_cluster_colors, collapse = ", "),
    call. = FALSE
  )
}

if (length(missing_cluster_labels) > 0L) {
  stop(
    "Missing anatomical labels for cluster(s): ",
    paste(missing_cluster_labels, collapse = ", "),
    call. = FALSE
  )
}

# ==============================================================================
# 7. Read and filter marker results
# ==============================================================================

markers <- read.delim(
  input_markers_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_marker_columns <- c(
  "p_val",
  "avg_log2FC",
  "pct.1",
  "pct.2",
  "p_val_adj",
  "cluster",
  "gene",
  "expression_specificity",
  "is_best_cluster",
  "tau",
  "gini",
  "shannon_specificity"
)

missing_marker_columns <- setdiff(
  required_marker_columns,
  colnames(markers)
)

if (length(missing_marker_columns) > 0L) {
  stop(
    "Marker table is missing required column(s): ",
    paste(missing_marker_columns, collapse = ", "),
    call. = FALSE
  )
}

markers$p_val <- suppressWarnings(
  as.numeric(markers$p_val)
)

markers$avg_log2FC <- suppressWarnings(
  as.numeric(markers$avg_log2FC)
)

markers$pct.1 <- suppressWarnings(
  as.numeric(markers$pct.1)
)

markers$pct.2 <- suppressWarnings(
  as.numeric(markers$pct.2)
)

markers$p_val_adj <- suppressWarnings(
  as.numeric(markers$p_val_adj)
)

markers$expression_specificity <- suppressWarnings(
  as.numeric(markers$expression_specificity)
)

markers$tau <- suppressWarnings(
  as.numeric(markers$tau)
)

markers$gini <- suppressWarnings(
  as.numeric(markers$gini)
)

markers$shannon_specificity <- suppressWarnings(
  as.numeric(markers$shannon_specificity)
)

markers$cluster <- as.character(
  markers$cluster
)

markers$gene <- as.character(
  markers$gene
)

if (!is.logical(markers$is_best_cluster)) {
  markers$is_best_cluster <- toupper(
    as.character(markers$is_best_cluster)
  ) == "TRUE"
}

selected_marker_rows <- markers[
  !is.na(markers$p_val_adj) &
    markers$p_val_adj <
      bonferroni_adjusted_p_value_threshold &
    !is.na(markers$tau) &
    markers$tau >= tau_threshold,
  ,
  drop = FALSE
]

if (nrow(selected_marker_rows) == 0L) {
  stop(
    "No markers satisfy Bonferroni-adjusted P < ",
    bonferroni_adjusted_p_value_threshold,
    " and tau >= ",
    tau_threshold,
    ".",
    call. = FALSE
  )
}

# ==============================================================================
# 8. Assign one representative marker cluster to each gene
# ==============================================================================

representative_order <- order(
  selected_marker_rows$gene,
  -as.integer(selected_marker_rows$is_best_cluster),
  -selected_marker_rows$expression_specificity,
  -selected_marker_rows$avg_log2FC,
  selected_marker_rows$p_val_adj,
  selected_marker_rows$cluster
)

selected_marker_rows <- selected_marker_rows[
  representative_order,
  ,
  drop = FALSE
]

selected_gene_table <- selected_marker_rows[
  !duplicated(selected_marker_rows$gene),
  ,
  drop = FALSE
]

selected_gene_table$marker_cluster <- as.character(
  selected_gene_table$cluster
)

selected_gene_table$marker_cluster_label <- unname(
  custom_cluster_labels[
    selected_gene_table$marker_cluster
  ]
)

selected_gene_table$marker_cluster_display <- paste0(
  "C",
  selected_gene_table$marker_cluster,
  " | ",
  selected_gene_table$marker_cluster_label
)

selected_gene_table <- selected_gene_table[
  order(
    match(selected_gene_table$marker_cluster, cluster_ids),
    -selected_gene_table$tau,
    -selected_gene_table$avg_log2FC,
    selected_gene_table$p_val_adj,
    selected_gene_table$gene
  ),
  ,
  drop = FALSE
]

selected_genes <- unique(
  selected_gene_table$gene
)

missing_genes <- setdiff(
  selected_genes,
  rownames(seurat_object[[assay_name]])
)

if (length(missing_genes) > 0L) {
  stop(
    "Selected marker genes are missing from the RNA assay: ",
    paste(missing_genes, collapse = ", "),
    call. = FALSE
  )
}

message(
  "\nSelected marker assignments before gene deduplication: ",
  nrow(selected_marker_rows)
)

message(
  "Selected unique genes: ",
  length(selected_genes)
)

# ==============================================================================
# 9. Extract normalized expression and back-transform from log1p
# ==============================================================================

expression_log <- LayerData(
  object = seurat_object,
  assay = assay_name,
  layer = "data",
  features = selected_genes
)

expression_log <- expression_log[
  selected_genes,
  ,
  drop = FALSE
]

spot_metadata <- spot_metadata[
  match(
    colnames(expression_log),
    rownames(spot_metadata)
  ),
  ,
  drop = FALSE
]

if (anyNA(rownames(spot_metadata))) {
  stop(
    "Failed to align Seurat metadata to expression-matrix columns.",
    call. = FALSE
  )
}

if (inherits(expression_log, "sparseMatrix")) {
  expression_linear <- expression_log
  expression_linear@x <- expm1(
    expression_linear@x
  )
} else {
  expression_linear <- expm1(
    expression_log
  )
}

# ==============================================================================
# 10. Prepare 16 samples repeated for each of the 16 clusters
# ==============================================================================

number_of_columns <-
  length(cluster_ids) *
  length(sample_order)

mean_expression_matrix <- matrix(
  NA_real_,
  nrow = length(selected_genes),
  ncol = number_of_columns,
  dimnames = list(
    selected_genes,
    NULL
  )
)

column_metadata <- data.frame(
  column_index = seq_len(number_of_columns),
  column_id = character(number_of_columns),
  cluster_id = character(number_of_columns),
  cluster_label = character(number_of_columns),
  sample_ID = character(number_of_columns),
  sex = character(number_of_columns),
  fmt_donor_group = character(number_of_columns),
  number_of_spots = integer(number_of_columns),
  stringsAsFactors = FALSE
)

column_index <- 0L

for (cluster_id in cluster_ids) {

  for (sample_id in sample_order) {

    column_index <- column_index + 1L

    selected_spots <- which(
      spot_metadata$cluster_ID_heatmap == cluster_id &
        spot_metadata$sample_ID_heatmap == sample_id
    )

    column_id <- paste0(
      "C",
      cluster_id,
      "__",
      sample_id
    )

    column_metadata$column_id[column_index] <- column_id
    column_metadata$cluster_id[column_index] <- cluster_id
    column_metadata$cluster_label[column_index] <-
      unname(custom_cluster_labels[cluster_id])
    column_metadata$sample_ID[column_index] <- sample_id
    column_metadata$number_of_spots[column_index] <-
      length(selected_spots)

    sample_metadata_rows <- which(
      spot_metadata$sample_ID_heatmap == sample_id
    )

    if (length(sample_metadata_rows) == 0L) {
      stop(
        "No metadata was found for sample ",
        sample_id,
        ".",
        call. = FALSE
      )
    }

    sample_sex_values <- unique(
      spot_metadata$sex_heatmap[
        sample_metadata_rows
      ]
    )

    sample_group_values <- unique(
      spot_metadata$group_heatmap[
        sample_metadata_rows
      ]
    )

    if (length(sample_sex_values) != 1L) {
      stop(
        "Sample ",
        sample_id,
        " has inconsistent sex metadata.",
        call. = FALSE
      )
    }

    if (length(sample_group_values) != 1L) {
      stop(
        "Sample ",
        sample_id,
        " has inconsistent FMT donor-group metadata.",
        call. = FALSE
      )
    }

    column_metadata$sex[column_index] <-
      sample_sex_values

    column_metadata$fmt_donor_group[column_index] <-
      sample_group_values

    if (length(selected_spots) > 0L) {
      mean_expression_matrix[
        ,
        column_index
      ] <- Matrix::rowMeans(
        expression_linear[
          ,
          selected_spots,
          drop = FALSE
        ]
      )
    }
  }
}

colnames(mean_expression_matrix) <-
  column_metadata$column_id

if (column_index != number_of_columns) {
  stop(
    "Unexpected number of sample-by-cluster columns.",
    call. = FALSE
  )
}

# ==============================================================================
# 11. Transform mean expression to row Z-scores
# ==============================================================================

log2_mean_expression_matrix <- log2(
  mean_expression_matrix + 1
)

zscore_matrix <- matrix(
  NA_real_,
  nrow = nrow(log2_mean_expression_matrix),
  ncol = ncol(log2_mean_expression_matrix),
  dimnames = dimnames(log2_mean_expression_matrix)
)

for (row_index in seq_len(nrow(log2_mean_expression_matrix))) {

  row_values <- as.numeric(
    log2_mean_expression_matrix[
      row_index,
      ,
      drop = TRUE
    ]
  )

  finite_values <- is.finite(
    row_values
  )

  if (sum(finite_values) < 2L) {
    zscore_matrix[
      row_index,
      finite_values
    ] <- 0
    next
  }

  row_mean <- mean(
    row_values[finite_values]
  )

  row_sd <- stats::sd(
    row_values[finite_values]
  )

  if (!is.finite(row_sd) || row_sd == 0) {
    zscore_matrix[
      row_index,
      finite_values
    ] <- 0
    next
  }

  zscore_matrix[
    row_index,
    finite_values
  ] <- (
    row_values[finite_values] -
      row_mean
  ) / row_sd
}

zscore_matrix[
  is.finite(zscore_matrix) &
    zscore_matrix < zscore_lower_limit
] <- zscore_lower_limit

zscore_matrix[
  is.finite(zscore_matrix) &
    zscore_matrix > zscore_upper_limit
] <- zscore_upper_limit

# ==============================================================================
# 12. Arrange rows by representative marker cluster
# ==============================================================================

selected_gene_table <- selected_gene_table[
  match(
    rownames(zscore_matrix),
    selected_gene_table$gene
  ),
  ,
  drop = FALSE
]

if (anyNA(selected_gene_table$gene)) {
  stop(
    "Failed to align selected-gene metadata to the heatmap matrix.",
    call. = FALSE
  )
}

row_split_levels <- paste0(
  "C",
  cluster_ids,
  " | ",
  unname(custom_cluster_labels[cluster_ids])
)

row_split <- factor(
  selected_gene_table$marker_cluster_display,
  levels = row_split_levels
)

column_split_levels <- paste0(
  "C",
  cluster_ids
)

column_split <- factor(
  paste0(
    "C",
    column_metadata$cluster_id
  ),
  levels = column_split_levels
)

# ==============================================================================
# 13. Column annotations
# ==============================================================================

cluster_annotation <- factor(
  column_metadata$cluster_id,
  levels = cluster_ids
)

sex_annotation <- factor(
  column_metadata$sex,
  levels = c(
    "Male",
    "Female"
  )
)

group_annotation <- factor(
  column_metadata$fmt_donor_group,
  levels = c(
    "Neurotypical",
    "ASD"
  )
)

if (anyNA(sex_annotation)) {
  stop(
    "Unexpected sex value in column metadata.",
    call. = FALSE
  )
}

if (anyNA(group_annotation)) {
  stop(
    "Unexpected FMT donor-group value in column metadata.",
    call. = FALSE
  )
}

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Cluster = cluster_annotation,
  Group = group_annotation,
  Sex = sex_annotation,
  col = list(
    Cluster = custom_cluster_colors[cluster_ids],
    Group = group_colors,
    Sex = sex_colors
  ),
  simple_anno_size = grid::unit(
    3.2,
    "mm"
  ),
  annotation_name_side = "left",
  annotation_name_gp = grid::gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  annotation_legend_param = list(
    Cluster = list(
      title = "Cluster",
      at = cluster_ids,
      labels = paste0(
        "C",
        cluster_ids,
        " | ",
        unname(custom_cluster_labels[cluster_ids])
      ),
      ncol = 2
    ),
    Group = list(
      title = "FMT donor group"
    ),
    Sex = list(
      title = "Sex"
    )
  )
)

# ==============================================================================
# 14. Heatmap title and dimensions
# ==============================================================================

number_of_selected_marker_assignments <-
  nrow(selected_marker_rows)

number_of_selected_unique_genes <-
  length(selected_genes)

number_of_sample_cluster_columns <-
  ncol(zscore_matrix)

number_of_missing_sample_cluster_profiles <-
  sum(column_metadata$number_of_spots == 0L)

heatmap_title <- paste(
  c(
    "Cluster-specific marker genes across samples and anatomical clusters",
    paste0(
      number_of_selected_unique_genes,
      " unique genes | ",
      number_of_selected_marker_assignments,
      " marker assignments | ",
      number_of_sample_cluster_columns,
      " sample-by-cluster columns"
    ),
    marker_description,
    clustering_description,
    paste0(
      "Values: mean normalized expression per sample x cluster",
      " -> log2(mean + 1) -> row Z-score clipped to [-2, 2]"
    ),
    paste0(
      "Samples are repeated within every cluster; sample IDs are not displayed",
      " | empty sample-cluster combinations: ",
      number_of_missing_sample_cluster_profiles
    )
  ),
  collapse = "\n"
)

if (number_of_selected_unique_genes <= 150L) {
  row_font_size <- 6.0
  row_height_mm <- 2.4
} else if (number_of_selected_unique_genes <= 350L) {
  row_font_size <- 4.8
  row_height_mm <- 1.8
} else if (number_of_selected_unique_genes <= 700L) {
  row_font_size <- 3.8
  row_height_mm <- 1.45
} else {
  row_font_size <- 3.2
  row_height_mm <- 1.20
}

column_width_mm <- 1.15

heatmap_body_width_mm <-
  number_of_sample_cluster_columns *
  column_width_mm

heatmap_body_height_mm <-
  number_of_selected_unique_genes *
  row_height_mm

pdf_width_inches <- max(
  16,
  (
    heatmap_body_width_mm +
      115
  ) / 25.4
)

pdf_height_inches <- max(
  14,
  (
    heatmap_body_height_mm +
      115
  ) / 25.4
)

message(
  "\nHeatmap dimensions: ",
  round(pdf_width_inches, 2),
  " x ",
  round(pdf_height_inches, 2),
  " inches"
)

# ==============================================================================
# 15. Create heatmap
# ==============================================================================

heatmap_color_function <- circlize::colorRamp2(
  c(
    zscore_lower_limit,
    0,
    zscore_upper_limit
  ),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)

heatmap_object <- ComplexHeatmap::Heatmap(
  zscore_matrix,
  name = "Row Z-score",
  col = heatmap_color_function,
  top_annotation = top_annotation,

  cluster_columns = FALSE,
  column_split = column_split,
  cluster_column_slices = FALSE,
  column_gap = grid::unit(
    1.0,
    "mm"
  ),

  cluster_rows = TRUE,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "complete",
  row_split = row_split,
  cluster_row_slices = FALSE,
  row_gap = grid::unit(
    0.8,
    "mm"
  ),
  row_dend_width = grid::unit(
    12,
    "mm"
  ),

  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(
    fontsize = row_font_size
  ),
  row_names_max_width =
    ComplexHeatmap::max_text_width(
      rownames(zscore_matrix),
      gp = grid::gpar(
        fontsize = row_font_size
      )
    ) +
    grid::unit(
      1.5,
      "mm"
    ),

  row_title_gp = grid::gpar(
    fontsize = 6.5,
    fontface = "bold"
  ),
  row_title_rot = 0,

  column_title = column_split_levels,
  column_title_gp = grid::gpar(
    fontsize = 8,
    fontface = "bold"
  ),

  na_col = "grey90",
  rect_gp = grid::gpar(
    col = NA
  ),
  border = FALSE,
  use_raster = TRUE,
  raster_quality = 2,

  width = grid::unit(
    heatmap_body_width_mm,
    "mm"
  ),
  height = grid::unit(
    heatmap_body_height_mm,
    "mm"
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
      45,
      "mm"
    )
  )
)

# ==============================================================================
# 16. Export selected genes and matrices
# ==============================================================================

selected_gene_output <- selected_gene_table[
  ,
  c(
    "gene",
    "marker_cluster",
    "marker_cluster_label",
    "p_val",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj",
    "expression_specificity",
    "is_best_cluster",
    "tau",
    "gini",
    "shannon_specificity"
  ),
  drop = FALSE
]

colnames(selected_gene_output) <- c(
  "gene",
  "representative_marker_cluster",
  "representative_marker_cluster_label",
  "p_value",
  "avg_log2FC",
  "pct_in_cluster",
  "pct_outside_cluster",
  "p_value_adj_bonferroni",
  "expression_specificity",
  "is_best_cluster",
  "tau",
  "gini",
  "shannon_specificity"
)

write.table(
  selected_gene_output,
  file = selected_genes_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

mean_expression_output <- data.frame(
  gene = rownames(mean_expression_matrix),
  mean_expression_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.table(
  mean_expression_output,
  file = mean_expression_matrix_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

log2_expression_output <- data.frame(
  gene = rownames(log2_mean_expression_matrix),
  log2_mean_expression_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.table(
  log2_expression_output,
  file = log2_expression_matrix_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

zscore_output <- data.frame(
  gene = rownames(zscore_matrix),
  zscore_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.table(
  zscore_output,
  file = zscore_matrix_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

write.table(
  column_metadata,
  file = column_metadata_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

summary_table <- data.frame(
  dataset_name = dataset_name,
  configuration_name = configuration_name,
  assay_name = assay_name,
  cluster_column = cluster_column,
  sample_column = sample_column,
  marker_test = "Wilcoxon rank-sum test",
  only_positive_markers = TRUE,
  minimum_detection_fraction = minimum_detection_fraction,
  minimum_average_log2_fold_change =
    minimum_average_log2_fold_change,
  bonferroni_adjusted_p_value_threshold =
    bonferroni_adjusted_p_value_threshold,
  tau_threshold = tau_threshold,
  selected_marker_assignments =
    number_of_selected_marker_assignments,
  selected_unique_genes =
    number_of_selected_unique_genes,
  samples = length(sample_order),
  clusters = length(cluster_ids),
  sample_cluster_columns =
    number_of_sample_cluster_columns,
  empty_sample_cluster_combinations =
    number_of_missing_sample_cluster_profiles,
  expression_input =
    "Seurat RNA data layer back-transformed using expm1",
  aggregation =
    "mean expression within each sample x cluster combination",
  transformation =
    "log2(mean expression + 1) -> row Z-score -> clip [-2, 2]",
  row_clustering =
    "Pearson correlation distance with complete linkage within representative marker-cluster splits",
  column_clustering = FALSE,
  sample_ID_labels_displayed = FALSE,
  output_pdf = heatmap_pdf,
  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = summary_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

# ==============================================================================
# 17. Save PDF
# ==============================================================================

message("\nWriting heatmap PDF:")
message(heatmap_pdf)

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
  column_title = heatmap_title,
  column_title_gp = grid::gpar(
    fontsize = 11,
    fontface = "bold",
    lineheight = 1.10
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
# 18. Validate outputs
# ==============================================================================

expected_output_files <- c(
  heatmap_pdf,
  selected_genes_tsv,
  mean_expression_matrix_tsv,
  log2_expression_matrix_tsv,
  zscore_matrix_tsv,
  column_metadata_tsv,
  summary_tsv
)

missing_output_files <- expected_output_files[
  !file.exists(expected_output_files)
]

if (length(missing_output_files) > 0L) {
  stop(
    "Missing output file(s):\n",
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

# ==============================================================================
# 19. Console summary
# ==============================================================================

message("\nHeatmap completed successfully.")
message("Selected unique genes: ", number_of_selected_unique_genes)
message(
  "Selected marker assignments: ",
  number_of_selected_marker_assignments
)
message(
  "Sample-by-cluster columns: ",
  number_of_sample_cluster_columns
)
message(
  "Empty sample-by-cluster combinations: ",
  number_of_missing_sample_cluster_profiles
)
message("\nOutput directory:")
message(
  normalizePath(
    output_directory,
    mustWork = TRUE
  )
)
message("\nHeatmap PDF:")
message(
  normalizePath(
    heatmap_pdf,
    mustWork = TRUE
  )
)

# ==============================================================================
# End
# ==============================================================================
