#!/usr/bin/env Rscript

# ==============================================================================
# 05_PCA_perCluster_crescendoCorrected_log2QN_topVariableGenes_slideID.R
#
# Exploratory PCA performed separately for each of the 16 anatomical clusters
# after Crescendo batch correction.
#
# INPUT
#   Per-cluster matrices:
#
#     Crescendo-corrected pseudobulk counts
#       -> log2(count + 1)
#       -> quantile normalization across 16 samples
#
# PCA
#   One PCA independently for each cluster.
#
#   Within each cluster:
#     - remove genes with zero / non-finite variance;
#     - select TOP 2000 most variable genes;
#     - prcomp(center = TRUE, scale. = TRUE);
#     - samples = 16 biological samples.
#
# PLOT
#   colour = slide_ID
#   shape  = FMT donor group
#   label  = sample_ID
#
#   Sample 13_1F is highlighted with a black outline.
#
# OUTPUT
#   1. One multi-page PDF containing all 16 cluster-specific PCA plots.
#   2. Sixteen individual cluster-specific PCA PDFs.
#   3. One TSV with PCA scores from all clusters.
#   4. One TSV with PCA variance explained per cluster.
#   5. One TSV with top variable genes used for PCA in every cluster.
#
# IMPORTANT
#   These PCA plots use expression data AFTER Crescendo correction.
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

number_of_top_variable_genes <- 2000L

highlight_sample <- "13_1F"

pca_center <- TRUE
pca_scale <- TRUE

# ==============================================================================
# 3. Cluster labels
# ==============================================================================

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
# 4. Input paths
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

# ==============================================================================
# 5. Output paths
# ==============================================================================

output_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,
  "02_heatmaps",
  "00_QC_residualBatchEffect",
  "02_PCA_perCluster_crescendoCorrected"
)

individual_pdf_dir <- file.path(
  output_dir,
  "01_perCluster_PDF"
)

dir.create(
  individual_pdf_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

combined_pdf_file <- file.path(
  output_dir,
  paste0(
    "01_",
    dataset_name,
    "_crescendoCorrected_log2QN_PCA_perCluster_top",
    number_of_top_variable_genes,
    ".pdf"
  )
)

all_scores_tsv <- file.path(
  output_dir,
  paste0(
    "02_",
    dataset_name,
    "_crescendoCorrected_log2QN_PCA_perCluster_scores.tsv"
  )
)

variance_summary_tsv <- file.path(
  output_dir,
  paste0(
    "03_",
    dataset_name,
    "_crescendoCorrected_log2QN_PCA_perCluster_varianceExplained.tsv"
  )
)

selected_features_tsv <- file.path(
  output_dir,
  paste0(
    "04_",
    dataset_name,
    "_crescendoCorrected_log2QN_PCA_perCluster_topVariableGenes.tsv"
  )
)

# ==============================================================================
# 6. Validate inputs
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
# 7. Load sample metadata
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
  stop(
    "Could not identify slide-ID column in sample_metadata. Available columns: ",
    paste(
      colnames(sample_metadata_raw),
      collapse = ", "
    ),
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

sample_order <- sample_metadata$sample_ID

print(
  sample_metadata,
  n = Inf
)

# ==============================================================================
# 8. Find cluster matrices
# ==============================================================================

message("\n============================================================")
message("STAGE 3: FINDING PER-CLUSTER MATRICES")
message("============================================================")

input_files <- list.files(
  log2_qn_input_dir,
  pattern =
    "_cluster_[0-9]+_crescendoCorrectedCounts_log2plus1_QN\\.tsv\\.gz$",
  full.names = TRUE
)

if (length(input_files) == 0L) {
  stop(
    "No corrected per-cluster log2+QN matrices found.",
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
    "Expected 16 clusters but found ",
    length(cluster_ids),
    ". Clusters: ",
    paste(
      cluster_ids,
      collapse = ", "
    ),
    call. = FALSE
  )
}

message(
  "Clusters found: ",
  paste(
    cluster_ids,
    collapse = ", "
  )
)

# ==============================================================================
# 9. Helper: run PCA for one cluster
# ==============================================================================

run_cluster_pca <- function(
    cluster_id_current,
    input_file_current
) {
  message("\n------------------------------------------------------------")
  message(
    "Cluster ",
    cluster_id_current,
    ": ",
    custom_cluster_labels[[cluster_id_current]]
  )
  message("------------------------------------------------------------")

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
      " is missing annotation columns: ",
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

  expression_matrix <- as.matrix(
    input_table[
      ,
      sample_order,
      drop = FALSE
    ]
  )

  storage.mode(expression_matrix) <- "double"

  rownames(expression_matrix) <-
    as.character(
      input_table$ensembl_gene_id
    )

  if (
    anyNA(expression_matrix) ||
      any(!is.finite(expression_matrix))
  ) {
    stop(
      "NA or non-finite expression values in cluster ",
      cluster_id_current,
      ".",
      call. = FALSE
    )
  }

  gene_variance <- apply(
    expression_matrix,
    1L,
    stats::var
  )

  valid_genes <- is.finite(
    gene_variance
  ) &
    gene_variance > 0

  number_valid_genes <- sum(
    valid_genes
  )

  if (number_valid_genes < 2L) {
    stop(
      "Too few variable genes for PCA in cluster ",
      cluster_id_current,
      ".",
      call. = FALSE
    )
  }

  number_genes_to_use <- min(
    number_of_top_variable_genes,
    number_valid_genes
  )

  valid_gene_indices <- which(
    valid_genes
  )

  selected_relative_order <- order(
    gene_variance[
      valid_gene_indices
    ],
    decreasing = TRUE
  )[
    seq_len(
      number_genes_to_use
    )
  ]

  selected_gene_indices <-
    valid_gene_indices[
      selected_relative_order
    ]

  selected_gene_ids <- rownames(
    expression_matrix
  )[
    selected_gene_indices
  ]

  expression_matrix_pca <- expression_matrix[
    selected_gene_indices,
    ,
    drop = FALSE
  ]

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
    ) |>
    dplyr::mutate(
      cluster_id =
        as.character(
          cluster_id_current
        ),
      anatomical_region =
        unname(
          custom_cluster_labels[
            cluster_id_current
          ]
        )
    )

  gene_display <- ifelse(
    is.na(input_table$gene) |
      input_table$gene == "",
    input_table$ensembl_gene_id,
    input_table$gene
  )

  feature_table <- tibble::tibble(
    cluster_id =
      as.character(
        cluster_id_current
      ),
    anatomical_region =
      unname(
        custom_cluster_labels[
          cluster_id_current
        ]
      ),
    variance_rank =
      seq_len(
        number_genes_to_use
      ),
    ensembl_gene_id =
      as.character(
        input_table$ensembl_gene_id[
          selected_gene_indices
        ]
      ),
    gene =
      as.character(
        gene_display[
          selected_gene_indices
        ]
      ),
    chromosome =
      as.character(
        input_table$chromosome[
          selected_gene_indices
        ]
      ),
    variance =
      as.numeric(
        gene_variance[
          selected_gene_indices
        ]
      )
  )

  variance_table <- tibble::tibble(
    cluster_id =
      as.character(
        cluster_id_current
      ),
    anatomical_region =
      unname(
        custom_cluster_labels[
          cluster_id_current
        ]
      ),
    genes_available =
      nrow(expression_matrix),
    genes_with_nonzero_variance =
      number_valid_genes,
    genes_used_for_PCA =
      number_genes_to_use,
    PC1_percent =
      variance_explained[[1]],
    PC2_percent =
      variance_explained[[2]],
    PC1_plus_PC2_percent =
      variance_explained[[1]] +
      variance_explained[[2]]
  )

  pca_plot <- ggplot2::ggplot(
    pca_scores,
    ggplot2::aes(
      x = .data$PC1,
      y = .data$PC2,
      colour = .data$slide_ID,
      shape =
        .data$fmt_donor_group
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
        label =
          .data$sample_ID
      ),
      size = 3.3,
      max.overlaps = Inf,
      box.padding = 0.45,
      point.padding = 0.30,
      min.segment.length = 0
    ) +
    ggplot2::labs(
      title = paste0(
        "C",
        cluster_id_current,
        ": ",
        custom_cluster_labels[[cluster_id_current]]
      ),
      subtitle = paste0(
        "PCA after Crescendo correction | top ",
        number_genes_to_use,
        " variable genes | colour = slide ID | shape = donor group"
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

  if (
    highlight_sample %in%
      pca_scores$sample_ID
  ) {
    highlight_data <- pca_scores |>
      dplyr::filter(
        .data$sample_ID ==
          .env$highlight_sample
      )

    pca_plot <- pca_plot +
      ggplot2::geom_point(
        data =
          highlight_data,
        mapping =
          ggplot2::aes(
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

  message(
    "Genes available: ",
    nrow(expression_matrix)
  )

  message(
    "Variable genes: ",
    number_valid_genes
  )

  message(
    "Genes used: ",
    number_genes_to_use
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

  list(
    plot = pca_plot,
    scores = pca_scores,
    variance_summary =
      variance_table,
    selected_features =
      feature_table
  )
}

# ==============================================================================
# 10. Run all 16 cluster-specific PCAs
# ==============================================================================

message("\n============================================================")
message("STAGE 4: RUNNING PCA FOR ALL CLUSTERS")
message("============================================================")

pca_results <- vector(
  "list",
  length(cluster_ids)
)

names(pca_results) <-
  cluster_ids

for (cluster_index in seq_along(cluster_ids)) {
  cluster_id_current <-
    cluster_ids[[cluster_index]]

  pca_results[[cluster_id_current]] <- run_cluster_pca(
    cluster_id_current =
      cluster_id_current,
    input_file_current =
      input_files[[cluster_index]]
  )
}

# ==============================================================================
# 11. Save individual PDFs
# ==============================================================================

message("\n============================================================")
message("STAGE 5: SAVING INDIVIDUAL PCA PDFs")
message("============================================================")

individual_pdf_files <- character(
  length(cluster_ids)
)

names(individual_pdf_files) <-
  cluster_ids

for (cluster_id_current in cluster_ids) {
  individual_pdf_file <- file.path(
    individual_pdf_dir,
    paste0(
      "C",
      sprintf(
        "%02d",
        as.integer(
          cluster_id_current
        )
      ),
      "_",
      dataset_name,
      "_crescendoCorrected_PCA.pdf"
    )
  )

  ggplot2::ggsave(
    filename =
      individual_pdf_file,
    plot =
      pca_results[[cluster_id_current]]$plot,
    width = 9.5,
    height = 7.5,
    units = "in",
    device = grDevices::cairo_pdf
  )

  individual_pdf_files[[cluster_id_current]] <- individual_pdf_file

  message(
    "Saved C",
    cluster_id_current,
    ": ",
    individual_pdf_file
  )
}

# ==============================================================================
# 12. Save one multi-page PDF
# ==============================================================================

message("\n============================================================")
message("STAGE 6: SAVING COMBINED MULTI-PAGE PDF")
message("============================================================")

grDevices::pdf(
  combined_pdf_file,
  width = 9.5,
  height = 7.5,
  onefile = TRUE,
  useDingbats = FALSE
)

for (cluster_id_current in cluster_ids) {
  print(
    pca_results[[cluster_id_current]]$plot
  )
}

grDevices::dev.off()

# ==============================================================================
# 13. Combine and save TSV outputs
# ==============================================================================

message("\n============================================================")
message("STAGE 7: SAVING PCA TABLES")
message("============================================================")

all_scores <- dplyr::bind_rows(
  lapply(
    cluster_ids,
    function(cluster_id_current) {
      pca_results[[cluster_id_current]]$scores
    }
  )
) |>
  dplyr::mutate(
    cluster_id = factor(
      .data$cluster_id,
      levels = cluster_ids
    )
  ) |>
  dplyr::arrange(
    .data$cluster_id,
    .data$sample_ID
  ) |>
  dplyr::mutate(
    cluster_id =
      as.character(
        .data$cluster_id
      )
  )

variance_summary <- dplyr::bind_rows(
  lapply(
    cluster_ids,
    function(cluster_id_current) {
      pca_results[[cluster_id_current]]$variance_summary
    }
  )
)

selected_features <- dplyr::bind_rows(
  lapply(
    cluster_ids,
    function(cluster_id_current) {
      pca_results[[cluster_id_current]]$selected_features
    }
  )
)

readr::write_tsv(
  all_scores,
  all_scores_tsv,
  na = "NA"
)

readr::write_tsv(
  variance_summary,
  variance_summary_tsv,
  na = "NA"
)

readr::write_tsv(
  selected_features,
  selected_features_tsv,
  na = "NA"
)

# ==============================================================================
# 14. Final report
# ==============================================================================

message("\n============================================================")
message("PER-CLUSTER PCA FINISHED")
message("============================================================")

message(
  "Clusters analysed: ",
  length(cluster_ids)
)

message(
  "\nCombined multi-page PDF:\n",
  normalizePath(
    combined_pdf_file,
    mustWork = TRUE
  )
)

message(
  "\nIndividual PDFs:\n",
  normalizePath(
    individual_pdf_dir,
    mustWork = TRUE
  )
)

message(
  "\nPCA scores:\n",
  normalizePath(
    all_scores_tsv,
    mustWork = TRUE
  )
)

message(
  "\nVariance explained summary:\n",
  normalizePath(
    variance_summary_tsv,
    mustWork = TRUE
  )
)

message(
  "\nTop variable genes:\n",
  normalizePath(
    selected_features_tsv,
    mustWork = TRUE
  )
)

# ==============================================================================
# End
# ==============================================================================
