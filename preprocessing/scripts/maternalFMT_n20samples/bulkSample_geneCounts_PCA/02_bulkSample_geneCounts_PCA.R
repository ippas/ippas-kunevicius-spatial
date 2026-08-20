# ==============================================================================
# 01b_bulkSample_allBarcodes_geneCounts_PCA.R
#
# Purpose:
# - read Space Ranger raw_feature_bc_matrix.h5 files
# - include all detected barcodes:
#     1. tissue-associated barcodes
#     2. background / off-tissue barcodes
# - sum counts across all detected barcodes for every gene within each sample
# - create one bulk-like gene-expression profile per sample
# - perform PCA on:
#     1. raw summed gene counts
#     2. log2(count + 1)-transformed and quantile-normalized counts
# - highlight selected samples
# - save PCA plots, PCA scores and input summaries
#
# This script does not use filtered_feature_bc_matrix.
# This script does not use read_spatial_samples().
#
# No clustering.
# No integration.
# No low-expression filtering.
# ==============================================================================


# ==============================================================================
# 1. Check and load packages
# ==============================================================================

required_packages <- c(
  "Seurat",
  "Matrix",
  "ggplot2",
  "ggrepel",
  "hdf5r"
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
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}


suppressPackageStartupMessages({

  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(ggrepel)

})


message(
  "Seurat version: ",
  as.character(packageVersion("Seurat"))
)


# ==============================================================================
# 2. Define input paths
# ==============================================================================

path_to_data <- file.path(
  "data",
  "spacerangerCount_withoutJSON_GRCm39-2024-A_2026-03-28"
)


metadata_file <- file.path(
  "data",
  "metadata_autismFMT.tsv"
)


# ==============================================================================
# 3. Define output paths
# ==============================================================================

output_dir <- file.path(
  "results",
  "maternalFMT_n20samples",
  "bulkSample_allDetectedBarcodes_geneCounts_PCA"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


combined_pdf_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_allDetectedBarcodes_",
    "geneCounts_PCA_raw_vs_log2Quantile.pdf"
  )
)


raw_pdf_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_allDetectedBarcodes_",
    "geneCounts_PCA_raw.pdf"
  )
)


log2_qn_pdf_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_allDetectedBarcodes_",
    "geneCounts_PCA_log2Quantile.pdf"
  )
)


pca_scores_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_allDetectedBarcodes_",
    "geneCounts_PCA_scores.tsv"
  )
)


input_summary_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_allDetectedBarcodes_",
    "input_summary.tsv"
  )
)


# ==============================================================================
# 4. Check input paths
# ==============================================================================

if (!dir.exists(path_to_data)) {

  stop(
    "Space Ranger data directory does not exist: ",
    path_to_data
  )
}


if (!file.exists(metadata_file)) {

  stop(
    "Metadata file does not exist: ",
    metadata_file
  )
}


# ==============================================================================
# 5. Read metadata
# ==============================================================================

metadata_autismFMT <- read.delim(
  metadata_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


if (!"sample_ID" %in% colnames(metadata_autismFMT)) {

  stop(
    "Column 'sample_ID' is missing from metadata."
  )
}


metadata_autismFMT$sample_ID <- as.character(
  metadata_autismFMT$sample_ID
)


if (anyNA(metadata_autismFMT$sample_ID)) {

  stop(
    "Metadata contains missing sample_ID values."
  )
}


if (any(metadata_autismFMT$sample_ID == "")) {

  stop(
    "Metadata contains empty sample_ID values."
  )
}


if (anyDuplicated(metadata_autismFMT$sample_ID)) {

  duplicated_samples <- unique(
    metadata_autismFMT$sample_ID[
      duplicated(metadata_autismFMT$sample_ID)
    ]
  )

  stop(
    "Duplicated sample_ID values: ",
    paste(
      duplicated_samples,
      collapse = ", "
    )
  )
}


expected_samples <- metadata_autismFMT$sample_ID


message(
  "Samples in metadata: ",
  length(expected_samples)
)


# ==============================================================================
# 6. Define samples to highlight
# ==============================================================================

# Four tissue sections with nearly identical expression profiles
# that did not correspond to their H&E images

red_samples <- c(
  "20_1F",
  "12_3F",
  "15_1M",
  "20_3M"
)


# Tissue section with a similar expression profile,
# but correctly corresponding to its H&E image

purple_samples <- "2_1M"


sample_colors <- c(
  suspicious = "red",
  reference = "purple",
  remaining = "black"
)


# ==============================================================================
# 7. Find the Space Ranger directory for one sample
# ==============================================================================

find_sample_directory <- function(
    root_directory,
    sample_id
) {

  exact_directory <- file.path(
    root_directory,
    sample_id
  )


  if (dir.exists(exact_directory)) {

    return(
      normalizePath(
        exact_directory,
        mustWork = TRUE
      )
    )
  }


  possible_directories <- list.dirs(
    root_directory,
    recursive = TRUE,
    full.names = TRUE
  )


  matching_directories <- possible_directories[
    basename(possible_directories) == sample_id
  ]


  if (length(matching_directories) == 0) {

    stop(
      "Could not find a directory for sample: ",
      sample_id
    )
  }


  if (length(matching_directories) > 1) {

    stop(
      "More than one directory was found for sample ",
      sample_id,
      ":\n",
      paste(
        matching_directories,
        collapse = "\n"
      )
    )
  }


  normalizePath(
    matching_directories,
    mustWork = TRUE
  )
}


# ==============================================================================
# 8. Find raw_feature_bc_matrix.h5 for one sample
# ==============================================================================

find_raw_h5_file <- function(
    sample_directory,
    sample_id
) {

  standard_candidates <- c(
    file.path(
      sample_directory,
      "outs",
      "raw_feature_bc_matrix.h5"
    ),
    file.path(
      sample_directory,
      "raw_feature_bc_matrix.h5"
    )
  )


  existing_standard_candidates <- standard_candidates[
    file.exists(standard_candidates)
  ]


  if (length(existing_standard_candidates) == 1) {

    return(
      normalizePath(
        existing_standard_candidates,
        mustWork = TRUE
      )
    )
  }


  if (length(existing_standard_candidates) > 1) {

    stop(
      "More than one standard raw H5 file was found for sample ",
      sample_id,
      ":\n",
      paste(
        existing_standard_candidates,
        collapse = "\n"
      )
    )
  }


  recursive_candidates <- list.files(
    sample_directory,
    pattern = "^raw_feature_bc_matrix\\.h5$",
    recursive = TRUE,
    full.names = TRUE
  )


  recursive_candidates <- recursive_candidates[
    file.exists(recursive_candidates)
  ]


  if (length(recursive_candidates) == 0) {

    stop(
      "raw_feature_bc_matrix.h5 was not found for sample ",
      sample_id,
      " in directory:\n",
      sample_directory
    )
  }


  if (length(recursive_candidates) > 1) {

    stop(
      "More than one raw_feature_bc_matrix.h5 file was found for sample ",
      sample_id,
      ":\n",
      paste(
        recursive_candidates,
        collapse = "\n"
      )
    )
  }


  normalizePath(
    recursive_candidates,
    mustWork = TRUE
  )
}


# ==============================================================================
# 9. Extract Gene Expression matrix from Read10X_h5 output
# ==============================================================================

extract_gene_expression_matrix <- function(
    read10x_object,
    sample_id
) {

  if (inherits(read10x_object, "Matrix")) {

    return(
      read10x_object
    )
  }


  if (is.matrix(read10x_object)) {

    return(
      Matrix::Matrix(
        read10x_object,
        sparse = TRUE
      )
    )
  }


  if (is.list(read10x_object)) {

    if ("Gene Expression" %in% names(read10x_object)) {

      return(
        read10x_object[["Gene Expression"]]
      )
    }


    if (length(read10x_object) == 1) {

      return(
        read10x_object[[1]]
      )
    }


    stop(
      "Multiple feature types were found for sample ",
      sample_id,
      ", but no 'Gene Expression' matrix was present. Available types: ",
      paste(
        names(read10x_object),
        collapse = ", "
      )
    )
  }


  stop(
    "Unsupported object returned by Read10X_h5 for sample: ",
    sample_id
  )
}


# ==============================================================================
# 10. Read raw matrices and sum all detected barcodes
# ==============================================================================

sample_gene_counts <- vector(
  mode = "list",
  length = length(expected_samples)
)


names(sample_gene_counts) <- expected_samples


input_summary <- vector(
  mode = "list",
  length = length(expected_samples)
)


names(input_summary) <- expected_samples


for (sample_id in expected_samples) {

  message(
    "\nProcessing sample: ",
    sample_id
  )


  sample_directory <- find_sample_directory(
    root_directory = path_to_data,
    sample_id = sample_id
  )


  raw_h5_file <- find_raw_h5_file(
    sample_directory = sample_directory,
    sample_id = sample_id
  )


  message(
    "  Raw matrix: ",
    raw_h5_file
  )


  raw_10x_object <- Seurat::Read10X_h5(
    filename = raw_h5_file,
    use.names = TRUE,
    unique.features = TRUE
  )


  counts_matrix <- extract_gene_expression_matrix(
    read10x_object = raw_10x_object,
    sample_id = sample_id
  )


  if (!inherits(counts_matrix, "Matrix")) {

    counts_matrix <- Matrix::Matrix(
      counts_matrix,
      sparse = TRUE
    )
  }


  if (nrow(counts_matrix) == 0) {

    stop(
      "No genes were found for sample: ",
      sample_id
    )
  }


  if (ncol(counts_matrix) == 0) {

    stop(
      "No barcodes were found for sample: ",
      sample_id
    )
  }


  if (is.null(rownames(counts_matrix))) {

    stop(
      "Gene names are missing for sample: ",
      sample_id
    )
  }


  if (is.null(colnames(counts_matrix))) {

    stop(
      "Barcode names are missing for sample: ",
      sample_id
    )
  }


  if (anyDuplicated(rownames(counts_matrix))) {

    stop(
      "Duplicated gene names remain after Read10X_h5 for sample: ",
      sample_id
    )
  }


  counts_per_barcode <- Matrix::colSums(
    counts_matrix
  )


  number_of_all_detected_barcodes <- ncol(
    counts_matrix
  )


  number_of_nonzero_umi_barcodes <- sum(
    counts_per_barcode > 0
  )


  # Sum tissue-associated and background barcode counts together

  gene_counts <- Matrix::rowSums(
    counts_matrix
  )


  gene_counts <- setNames(
    as.numeric(gene_counts),
    rownames(counts_matrix)
  )


  sample_gene_counts[[sample_id]] <- gene_counts


  input_summary[[sample_id]] <- data.frame(
    sample_ID = sample_id,
    raw_matrix_file = raw_h5_file,
    genes = nrow(counts_matrix),
    all_detected_barcodes = number_of_all_detected_barcodes,
    barcodes_with_nonzero_UMI = number_of_nonzero_umi_barcodes,
    total_summed_UMIs = sum(gene_counts),
    stringsAsFactors = FALSE
  )


  message(
    "  genes: ",
    format(
      nrow(counts_matrix),
      big.mark = ","
    )
  )


  message(
    "  all detected barcodes: ",
    format(
      number_of_all_detected_barcodes,
      big.mark = ","
    )
  )


  message(
    "  barcodes with nonzero UMI count: ",
    format(
      number_of_nonzero_umi_barcodes,
      big.mark = ","
    )
  )


  message(
    "  total summed UMIs: ",
    format(
      sum(gene_counts),
      big.mark = ",",
      scientific = FALSE
    )
  )


  rm(
    raw_10x_object,
    counts_matrix,
    counts_per_barcode
  )


  invisible(
    gc()
  )
}


input_summary <- do.call(
  rbind,
  input_summary
)


rownames(input_summary) <- NULL


# ==============================================================================
# 11. Check and align genes between samples
# ==============================================================================

reference_genes <- names(
  sample_gene_counts[[1]]
)


identical_gene_sets <- vapply(
  sample_gene_counts,
  function(gene_counts) {

    setequal(
      names(gene_counts),
      reference_genes
    )
  },
  FUN.VALUE = logical(1)
)


if (!all(identical_gene_sets)) {

  problematic_samples <- names(
    identical_gene_sets
  )[!identical_gene_sets]


  stop(
    "Gene sets differ between samples: ",
    paste(
      problematic_samples,
      collapse = ", "
    )
  )
}


# Reorder every sample according to the first sample

sample_gene_counts <- lapply(
  sample_gene_counts,
  function(gene_counts) {

    gene_counts[
      reference_genes
    ]
  }
)


# ==============================================================================
# 12. Create gene x sample raw-count matrix
# ==============================================================================

gene_counts_per_sample_raw <- do.call(
  cbind,
  sample_gene_counts
)


rownames(
  gene_counts_per_sample_raw
) <- reference_genes


colnames(
  gene_counts_per_sample_raw
) <- expected_samples


storage.mode(
  gene_counts_per_sample_raw
) <- "numeric"


if (anyNA(gene_counts_per_sample_raw)) {

  stop(
    "Bulk all-barcode count matrix contains missing values."
  )
}


if (any(!is.finite(gene_counts_per_sample_raw))) {

  stop(
    "Bulk all-barcode count matrix contains non-finite values."
  )
}


if (any(gene_counts_per_sample_raw < 0)) {

  stop(
    "Bulk all-barcode count matrix contains negative values."
  )
}


message(
  "\nBulk all-barcode matrix dimensions: ",
  nrow(gene_counts_per_sample_raw),
  " genes x ",
  ncol(gene_counts_per_sample_raw),
  " samples"
)


# ==============================================================================
# 13. Define quantile-normalization function
# ==============================================================================

quantile_normalize <- function(expression_matrix) {

  expression_matrix <- as.matrix(
    expression_matrix
  )


  if (anyNA(expression_matrix)) {

    stop(
      "Expression matrix contains missing values."
    )
  }


  if (any(!is.finite(expression_matrix))) {

    stop(
      "Expression matrix contains non-finite values."
    )
  }


  if (ncol(expression_matrix) < 2) {

    stop(
      "Quantile normalization requires at least two samples."
    )
  }


  sorted_matrix <- apply(
    expression_matrix,
    MARGIN = 2,
    FUN = sort
  )


  target_distribution <- rowMeans(
    sorted_matrix
  )


  normalized_matrix <- matrix(
    NA_real_,
    nrow = nrow(expression_matrix),
    ncol = ncol(expression_matrix),
    dimnames = dimnames(expression_matrix)
  )


  for (sample_index in seq_len(ncol(expression_matrix))) {

    sample_values <- expression_matrix[
      ,
      sample_index
    ]


    sample_order <- order(
      sample_values
    )


    sorted_sample_values <- sample_values[
      sample_order
    ]


    normalized_sorted_values <- numeric(
      length(sorted_sample_values)
    )


    ties <- rle(
      sorted_sample_values
    )


    tie_end <- cumsum(
      ties$lengths
    )


    tie_start <- tie_end -
      ties$lengths +
      1


    for (tie_index in seq_along(ties$values)) {

      tie_positions <- tie_start[tie_index]:
        tie_end[tie_index]


      normalized_sorted_values[
        tie_positions
      ] <- mean(
        target_distribution[
          tie_positions
        ]
      )
    }


    normalized_matrix[
      sample_order,
      sample_index
    ] <- normalized_sorted_values
  }


  normalized_matrix
}


# ==============================================================================
# 14. Define PCA function
# ==============================================================================

run_sample_pca <- function(
    expression_matrix,
    analysis_id
) {

  expression_matrix <- as.matrix(
    expression_matrix
  )


  if (anyNA(expression_matrix)) {

    stop(
      analysis_id,
      ": expression matrix contains missing values."
    )
  }


  if (any(!is.finite(expression_matrix))) {

    stop(
      analysis_id,
      ": expression matrix contains non-finite values."
    )
  }


  if (nrow(expression_matrix) < 2) {

    stop(
      analysis_id,
      ": fewer than two genes are available."
    )
  }


  if (ncol(expression_matrix) < 3) {

    stop(
      analysis_id,
      ": fewer than three samples are available."
    )
  }


  # prcomp requires:
  # rows    = samples
  # columns = genes

  pca_result <- prcomp(
    x = t(expression_matrix),
    center = TRUE,
    scale. = FALSE
  )


  explained_variance <- (
    pca_result$sdev^2 /
      sum(pca_result$sdev^2)
  ) * 100


  pca_df <- data.frame(
    sample_ID = rownames(pca_result$x),
    PC1 = pca_result$x[, "PC1"],
    PC2 = pca_result$x[, "PC2"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )


  pca_df$sample_group <- ifelse(
    pca_df$sample_ID %in% red_samples,
    "suspicious",
    ifelse(
      pca_df$sample_ID %in% purple_samples,
      "reference",
      "remaining"
    )
  )


  list(
    analysis_id = analysis_id,
    input_matrix = expression_matrix,
    pca_result = pca_result,
    explained_variance = explained_variance,
    pca_df = pca_df
  )
}


# ==============================================================================
# 15. Define PCA plotting function
# ==============================================================================

make_pca_plot <- function(
    pca_object,
    data_processing_label
) {

  number_of_genes <- nrow(
    pca_object$input_matrix
  )


  number_of_samples <- ncol(
    pca_object$input_matrix
  )


  complete_plot_title <- paste0(
    "Bulk PCA: all detected barcodes summed per gene within each sample",
    "\n",
    "Input: raw_feature_bc_matrix (tissue-associated + background barcodes)",
    "\n",
    data_processing_label,
    "\n",
    number_of_genes,
    " genes | ",
    number_of_samples,
    " samples"
  )


  ggplot(
    data = pca_object$pca_df,
    mapping = aes(
      x = PC1,
      y = PC2
    )
  ) +

    geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      linetype = "dashed"
    ) +

    geom_vline(
      xintercept = 0,
      linewidth = 0.3,
      linetype = "dashed"
    ) +

    geom_point(
      aes(
        color = sample_group
      ),
      size = 3
    ) +

    ggrepel::geom_text_repel(
      aes(
        label = sample_ID,
        color = sample_group
      ),
      size = 3.2,
      seed = 7,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.35,
      min.segment.length = 0,
      segment.size = 0.25,
      show.legend = FALSE
    ) +

    scale_color_manual(
      values = sample_colors
    ) +

    labs(
      title = complete_plot_title,
      x = paste0(
        "PC1 (",
        round(
          pca_object$explained_variance[1],
          1
        ),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(
          pca_object$explained_variance[2],
          1
        ),
        "%)"
      )
    ) +

    coord_cartesian(
      clip = "off"
    ) +

    theme_classic(
      base_size = 11
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 10.5,
        lineheight = 1.08
      ),
      plot.margin = margin(
        t = 12,
        r = 25,
        b = 10,
        l = 10
      ),
      legend.position = "none"
    )
}


# ==============================================================================
# 16. PCA on raw counts summed across all detected barcodes
# ==============================================================================

pca_raw <- run_sample_pca(
  expression_matrix = gene_counts_per_sample_raw,
  analysis_id = "raw_all_detected_barcodes_gene_counts"
)


# ==============================================================================
# 17. Log2 transformation
# ==============================================================================

gene_counts_per_sample_log2 <- log2(
  gene_counts_per_sample_raw + 1
)


# ==============================================================================
# 18. Quantile normalization
# ==============================================================================

gene_counts_per_sample_log2_qn <- quantile_normalize(
  gene_counts_per_sample_log2
)


stopifnot(
  identical(
    dim(gene_counts_per_sample_raw),
    dim(gene_counts_per_sample_log2_qn)
  )
)


stopifnot(
  identical(
    rownames(gene_counts_per_sample_raw),
    rownames(gene_counts_per_sample_log2_qn)
  )
)


stopifnot(
  identical(
    colnames(gene_counts_per_sample_raw),
    colnames(gene_counts_per_sample_log2_qn)
  )
)


# ==============================================================================
# 19. PCA on log2 + quantile-normalized counts
# ==============================================================================

pca_log2_qn <- run_sample_pca(
  expression_matrix = gene_counts_per_sample_log2_qn,
  analysis_id = paste0(
    "log2_quantile_normalized_",
    "all_detected_barcodes_gene_counts"
  )
)


# ==============================================================================
# 20. Create PCA plots
# ==============================================================================

pca_plot_raw <- make_pca_plot(
  pca_object = pca_raw,
  data_processing_label = paste0(
    "Raw summed gene counts; ",
    "no transformation or normalization"
  )
)


pca_plot_log2_qn <- make_pca_plot(
  pca_object = pca_log2_qn,
  data_processing_label = paste0(
    "Summed gene counts after log2(count + 1) transformation ",
    "and quantile normalization"
  )
)


# ==============================================================================
# 21. Draw combined PCA grid
# ==============================================================================

draw_pca_grid <- function() {

  grid::grid.newpage()


  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(
        nrow = 1,
        ncol = 2
      )
    )
  )


  print(
    pca_plot_raw,
    vp = grid::viewport(
      layout.pos.row = 1,
      layout.pos.col = 1
    )
  )


  print(
    pca_plot_log2_qn,
    vp = grid::viewport(
      layout.pos.row = 1,
      layout.pos.col = 2
    )
  )
}


# ==============================================================================
# 22. Prepare PCA-score table
# ==============================================================================

prepare_pca_scores <- function(
    pca_object,
    analysis_id
) {

  pca_scores <- as.data.frame(
    pca_object$pca_result$x,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


  pca_scores$sample_ID <- rownames(
    pca_scores
  )


  rownames(pca_scores) <- NULL


  pca_scores$sample_group <- ifelse(
    pca_scores$sample_ID %in% red_samples,
    "suspicious",
    ifelse(
      pca_scores$sample_ID %in% purple_samples,
      "reference",
      "remaining"
    )
  )


  pca_scores$analysis <- analysis_id


  pca_scores$genes_used <- nrow(
    pca_object$input_matrix
  )


  pca_scores$samples_used <- ncol(
    pca_object$input_matrix
  )


  pca_scores$barcode_set <- paste0(
    "raw_feature_bc_matrix:",
    "tissue_associated_plus_background"
  )


  pca_scores$PC1_variance_percent <-
    pca_object$explained_variance[1]


  pca_scores$PC2_variance_percent <-
    pca_object$explained_variance[2]


  pc_columns <- grep(
    "^PC[0-9]+$",
    colnames(pca_scores),
    value = TRUE
  )


  pca_scores[
    ,
    c(
      "analysis",
      "sample_ID",
      "sample_group",
      "barcode_set",
      "genes_used",
      "samples_used",
      "PC1_variance_percent",
      "PC2_variance_percent",
      pc_columns
    ),
    drop = FALSE
  ]
}


pca_scores_raw <- prepare_pca_scores(
  pca_object = pca_raw,
  analysis_id = "raw_all_detected_barcodes_gene_counts"
)


pca_scores_log2_qn <- prepare_pca_scores(
  pca_object = pca_log2_qn,
  analysis_id = paste0(
    "log2_quantile_normalized_",
    "all_detected_barcodes_gene_counts"
  )
)


pca_scores_combined <- rbind(
  pca_scores_raw,
  pca_scores_log2_qn
)


# ==============================================================================
# 23. Save tables
# ==============================================================================

write.table(
  pca_scores_combined,
  file = pca_scores_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  input_summary,
  file = input_summary_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ==============================================================================
# 24. Save separate PCA PDFs
# ==============================================================================

ggsave(
  filename = raw_pdf_file,
  plot = pca_plot_raw,
  device = "pdf",
  width = 9,
  height = 8,
  units = "in",
  limitsize = FALSE
)


ggsave(
  filename = log2_qn_pdf_file,
  plot = pca_plot_log2_qn,
  device = "pdf",
  width = 9,
  height = 8,
  units = "in",
  limitsize = FALSE
)


# ==============================================================================
# 25. Save combined PCA PDF
# ==============================================================================

grDevices::pdf(
  file = combined_pdf_file,
  width = 17,
  height = 8,
  onefile = TRUE
)


draw_pca_grid()


grDevices::dev.off()


# ==============================================================================
# 26. Display combined PCA figure
# ==============================================================================

draw_pca_grid()


# ==============================================================================
# 27. Final messages
# ==============================================================================

message(
  "\nCombined PCA figure saved to: ",
  normalizePath(
    combined_pdf_file,
    mustWork = FALSE
  )
)


message(
  "Raw PCA figure saved to: ",
  normalizePath(
    raw_pdf_file,
    mustWork = FALSE
  )
)


message(
  "Log2 + quantile-normalized PCA figure saved to: ",
  normalizePath(
    log2_qn_pdf_file,
    mustWork = FALSE
  )
)


message(
  "PCA scores saved to: ",
  normalizePath(
    pca_scores_file,
    mustWork = FALSE
  )
)


message(
  "Input summary saved to: ",
  normalizePath(
    input_summary_file,
    mustWork = FALSE
  )
)


message(
  "All-detected-barcode bulk PCA completed."
)