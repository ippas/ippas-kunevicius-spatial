# ==============================================================================
# 02_bulkSample_allBarcodes_allGenes_directUMAP.R
#
# Purpose:
# - read raw_feature_bc_matrix.h5 for all Maternal FMT samples
# - include all barcodes present in the raw Space Ranger matrix:
#     1. in-tissue barcodes
#     2. off-tissue/background barcodes
# - sum counts across all barcodes for every gene within every sample
# - retain all genes without expression or variability filtering
# - calculate exactly two direct UMAP embeddings:
#     1. raw summed gene counts
#     2. log2(count + 1) followed by quantile normalization
# - highlight selected tissue sections
# - save separate and combined UMAP figures
#
# All barcodes.
# All genes.
# No PCA.
# No clustering.
# No integration.
# No gene filtering.
# ==============================================================================


# ==============================================================================
# 1. Check required packages
# ==============================================================================

required_packages <- c(
  "Seurat",
  "Matrix",
  "ggplot2",
  "ggrepel",
  "uwot",
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


# ==============================================================================
# 2. Load packages
# ==============================================================================

suppressPackageStartupMessages({

  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(ggrepel)
  library(uwot)

})


message(
  "Seurat version: ",
  as.character(packageVersion("Seurat"))
)


message(
  "uwot version: ",
  as.character(packageVersion("uwot"))
)


# ==============================================================================
# 3. Define input paths
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
# 4. Define output paths
# ==============================================================================

output_dir <- file.path(
  "results",
  "maternalFMT_n20samples",
  "bulkSample_allBarcodes_allGenes_directUMAP"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


combined_pdf_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_allBarcodes_allGenes_",
    "directUMAP_raw_vs_log2Quantile.pdf"
  )
)


raw_pdf_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_allBarcodes_allGenes_",
    "directUMAP_raw.pdf"
  )
)


log2_qn_pdf_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_allBarcodes_allGenes_",
    "directUMAP_log2Quantile.pdf"
  )
)


umap_coordinates_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_allBarcodes_allGenes_",
    "directUMAP_coordinates.tsv"
  )
)


input_summary_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_allBarcodes_allGenes_",
    "input_summary.tsv"
  )
)


# ==============================================================================
# 5. Define UMAP settings
# ==============================================================================

umap_n_neighbors <- 5

umap_min_dist <- 0.3

umap_metric <- "euclidean"

random_seed <- 7


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
# 7. Check input paths
# ==============================================================================

if (!dir.exists(path_to_data)) {

  stop(
    "Space Ranger output directory does not exist: ",
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
# 8. Read metadata
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


missing_highlighted_samples <- setdiff(
  c(
    red_samples,
    purple_samples
  ),
  expected_samples
)


if (length(missing_highlighted_samples) > 0) {

  warning(
    "Highlighted samples absent from metadata: ",
    paste(
      missing_highlighted_samples,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 9. Find the Space Ranger directory for one sample
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
      "Could not find a Space Ranger directory for sample: ",
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
# 10. Find raw_feature_bc_matrix.h5 for one sample
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
      "More than one raw_feature_bc_matrix.h5 was found for sample ",
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
# 11. Extract Gene Expression matrix from Read10X_h5 output
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
      ", but no 'Gene Expression' matrix was available. Types: ",
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
# 12. Read raw matrices and sum counts across all barcodes
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
    "  raw matrix: ",
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


  if (any(counts_matrix@x < 0)) {

    stop(
      "Negative counts were found for sample: ",
      sample_id
    )
  }


  counts_per_barcode <- Matrix::colSums(
    counts_matrix
  )


  # Sum all raw-matrix barcodes for every gene

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
    all_raw_matrix_barcodes = ncol(counts_matrix),
    barcodes_with_nonzero_UMI = sum(
      counts_per_barcode > 0
    ),
    barcodes_with_zero_UMI = sum(
      counts_per_barcode == 0
    ),
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
    "  all raw-matrix barcodes: ",
    format(
      ncol(counts_matrix),
      big.mark = ","
    )
  )


  message(
    "  barcodes with nonzero UMI: ",
    format(
      sum(counts_per_barcode > 0),
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
    counts_per_barcode,
    gene_counts
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
# 13. Check and align genes between samples
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


# Put genes in exactly the same order in every sample

sample_gene_counts <- lapply(
  sample_gene_counts,
  function(gene_counts) {

    gene_counts[
      reference_genes
    ]
  }
)


# ==============================================================================
# 14. Create gene x sample raw-count matrix
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
) <- "double"


if (anyNA(gene_counts_per_sample_raw)) {

  stop(
    "Raw all-barcode count matrix contains missing values."
  )
}


if (any(!is.finite(gene_counts_per_sample_raw))) {

  stop(
    "Raw all-barcode count matrix contains non-finite values."
  )
}


if (any(gene_counts_per_sample_raw < 0)) {

  stop(
    "Raw all-barcode count matrix contains negative values."
  )
}


message(
  "\nAll-barcode bulk matrix: ",
  nrow(gene_counts_per_sample_raw),
  " genes x ",
  ncol(gene_counts_per_sample_raw),
  " samples"
)


message(
  "All genes retained for both UMAP analyses: ",
  nrow(gene_counts_per_sample_raw)
)


# ==============================================================================
# 15. Define quantile-normalization function
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


  # Sort expression values within each sample

  sorted_matrix <- apply(
    expression_matrix,
    MARGIN = 2,
    FUN = sort
  )


  # Calculate the common target distribution

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


    # Give tied original values the same normalized value

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
# 16. Create log2(count + 1) matrix
# ==============================================================================

gene_counts_per_sample_log2 <- log2(
  gene_counts_per_sample_raw + 1
)


if (anyNA(gene_counts_per_sample_log2)) {

  stop(
    "Log2-transformed matrix contains missing values."
  )
}


if (any(!is.finite(gene_counts_per_sample_log2))) {

  stop(
    "Log2-transformed matrix contains non-finite values."
  )
}


# ==============================================================================
# 17. Quantile-normalize log2(count + 1) matrix
# ==============================================================================

gene_counts_per_sample_log2_qn <- quantile_normalize(
  gene_counts_per_sample_log2
)


# Confirm that both UMAP inputs contain exactly
# the same genes and samples

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


# Intermediate log2 matrix is no longer needed

rm(
  gene_counts_per_sample_log2
)


invisible(
  gc()
)


# ==============================================================================
# 18. Define direct UMAP function
# ==============================================================================

run_direct_umap <- function(
    expression_matrix,
    analysis_id,
    n_neighbors,
    min_dist,
    metric,
    seed
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


  number_of_genes <- nrow(
    expression_matrix
  )


  number_of_samples <- ncol(
    expression_matrix
  )


  if (number_of_genes < 2) {

    stop(
      analysis_id,
      ": fewer than two genes are available."
    )
  }


  if (number_of_samples < 4) {

    stop(
      analysis_id,
      ": fewer than four samples are available."
    )
  }


  # Current matrix:
  # rows    = genes
  # columns = samples
  #
  # UMAP input:
  # rows    = samples
  # columns = genes

  sample_by_gene_matrix <- t(
    expression_matrix
  )


  n_neighbors_used <- min(
    n_neighbors,
    number_of_samples - 1
  )


  if (n_neighbors_used < 2) {

    stop(
      analysis_id,
      ": n_neighbors must be at least 2."
    )
  }


  message(
    "\n",
    analysis_id,
    ": direct UMAP using ",
    number_of_samples,
    " samples and ",
    number_of_genes,
    " genes"
  )


  message(
    "  n_neighbors = ",
    n_neighbors_used,
    "; min_dist = ",
    min_dist,
    "; metric = ",
    metric,
    "; PCA = disabled"
  )


  umap_matrix <- uwot::umap(
    X = sample_by_gene_matrix,

    n_neighbors = n_neighbors_used,
    n_components = 2,

    metric = metric,
    min_dist = min_dist,

    scale = FALSE,

    # No PCA preprocessing

    pca = NULL,

    init = "spectral",

    seed = seed,

    # Single-threaded optimization for reproducibility

    n_threads = 1,
    n_sgd_threads = 1,

    fast_sgd = FALSE,

    ret_model = FALSE,
    verbose = TRUE
  )


  rownames(
    umap_matrix
  ) <- rownames(
    sample_by_gene_matrix
  )


  colnames(
    umap_matrix
  ) <- c(
    "UMAP1",
    "UMAP2"
  )


  umap_df <- data.frame(
    sample_ID = rownames(umap_matrix),
    UMAP1 = umap_matrix[, "UMAP1"],
    UMAP2 = umap_matrix[, "UMAP2"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )


  umap_df$sample_group <- ifelse(
    umap_df$sample_ID %in% red_samples,
    "suspicious",
    ifelse(
      umap_df$sample_ID %in% purple_samples,
      "reference",
      "remaining"
    )
  )


  list(
    analysis_id = analysis_id,
    umap_matrix = umap_matrix,
    umap_df = umap_df,
    genes_used = number_of_genes,
    samples_used = number_of_samples,
    n_neighbors = n_neighbors_used,
    min_dist = min_dist,
    metric = metric,
    seed = seed
  )
}


# ==============================================================================
# 19. UMAP 1: raw counts, all barcodes and all genes
# ==============================================================================

umap_raw <- run_direct_umap(
  expression_matrix = gene_counts_per_sample_raw,
  analysis_id = paste0(
    "raw_counts_",
    "all_barcodes_",
    "all_genes_",
    "direct_gene_space"
  ),
  n_neighbors = umap_n_neighbors,
  min_dist = umap_min_dist,
  metric = umap_metric,
  seed = random_seed
)


# ==============================================================================
# 20. UMAP 2: log2(count + 1) + quantile normalization,
#     all barcodes and all genes
# ==============================================================================

umap_log2_qn <- run_direct_umap(
  expression_matrix = gene_counts_per_sample_log2_qn,
  analysis_id = paste0(
    "log2_quantile_normalized_",
    "all_barcodes_",
    "all_genes_",
    "direct_gene_space"
  ),
  n_neighbors = umap_n_neighbors,
  min_dist = umap_min_dist,
  metric = umap_metric,
  seed = random_seed
)


# ==============================================================================
# 21. Define UMAP plotting function
# ==============================================================================

make_direct_umap_plot <- function(
    umap_object,
    processing_label
) {

  complete_plot_title <- paste0(
    "Direct bulk UMAP: all raw-matrix barcodes summed per gene within each sample",
    "\n",
    processing_label,
    "\n",
    umap_object$genes_used,
    " genes | ",
    umap_object$samples_used,
    " samples | all genes; no PCA"
  )


  ggplot(
    data = umap_object$umap_df,
    mapping = aes(
      x = UMAP1,
      y = UMAP2
    )
  ) +

    geom_point(
      aes(
        color = sample_group
      ),
      size = 3.5
    ) +

    ggrepel::geom_text_repel(
      aes(
        label = sample_ID,
        color = sample_group
      ),
      size = 3.4,
      seed = random_seed,
      max.overlaps = Inf,
      box.padding = 0.55,
      point.padding = 0.4,
      min.segment.length = 0,
      segment.size = 0.25,
      show.legend = FALSE
    ) +

    scale_color_manual(
      values = sample_colors
    ) +

    labs(
      title = complete_plot_title,
      x = "UMAP1",
      y = "UMAP2"
    ) +

    coord_equal(
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
        r = 28,
        b = 12,
        l = 12
      ),
      legend.position = "none"
    )
}


# ==============================================================================
# 22. Create exactly two UMAP plots
# ==============================================================================

umap_plot_raw <- make_direct_umap_plot(
  umap_object = umap_raw,
  processing_label = paste0(
    "Raw counts summed across all barcodes; ",
    "no transformation or normalization"
  )
)


umap_plot_log2_qn <- make_direct_umap_plot(
  umap_object = umap_log2_qn,
  processing_label = paste0(
    "Counts summed across all barcodes; ",
    "log2(count + 1) transformation and quantile normalization"
  )
)


# ==============================================================================
# 23. Define function for drawing the two-panel UMAP figure
# ==============================================================================

draw_umap_grid <- function() {

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
    umap_plot_raw,
    vp = grid::viewport(
      layout.pos.row = 1,
      layout.pos.col = 1
    )
  )


  print(
    umap_plot_log2_qn,
    vp = grid::viewport(
      layout.pos.row = 1,
      layout.pos.col = 2
    )
  )
}


# ==============================================================================
# 24. Prepare UMAP-coordinate tables
# ==============================================================================

prepare_umap_coordinates <- function(
    umap_object,
    processing_label
) {

  umap_coordinates <- umap_object$umap_df


  summary_order <- match(
    umap_coordinates$sample_ID,
    input_summary$sample_ID
  )


  if (anyNA(summary_order)) {

    stop(
      "Could not match all UMAP samples to input summary."
    )
  }


  matching_summary <- input_summary[
    summary_order,
    ,
    drop = FALSE
  ]


  umap_coordinates$analysis <-
    umap_object$analysis_id


  umap_coordinates$processing <-
    processing_label


  umap_coordinates$barcode_input <-
    "raw_feature_bc_matrix_all_barcodes"


  umap_coordinates$input_space <-
    "all_genes_direct_no_PCA"


  umap_coordinates$genes_used <-
    umap_object$genes_used


  umap_coordinates$samples_used <-
    umap_object$samples_used


  umap_coordinates$all_raw_matrix_barcodes <-
    matching_summary$all_raw_matrix_barcodes


  umap_coordinates$barcodes_with_nonzero_UMI <-
    matching_summary$barcodes_with_nonzero_UMI


  umap_coordinates$total_summed_UMIs <-
    matching_summary$total_summed_UMIs


  umap_coordinates$n_neighbors <-
    umap_object$n_neighbors


  umap_coordinates$min_dist <-
    umap_object$min_dist


  umap_coordinates$metric <-
    umap_object$metric


  umap_coordinates$seed <-
    umap_object$seed


  umap_coordinates[
    ,
    c(
      "analysis",
      "processing",
      "barcode_input",
      "input_space",
      "sample_ID",
      "sample_group",
      "genes_used",
      "samples_used",
      "all_raw_matrix_barcodes",
      "barcodes_with_nonzero_UMI",
      "total_summed_UMIs",
      "n_neighbors",
      "min_dist",
      "metric",
      "seed",
      "UMAP1",
      "UMAP2"
    ),
    drop = FALSE
  ]
}


umap_coordinates_raw <- prepare_umap_coordinates(
  umap_object = umap_raw,
  processing_label = "raw_gene_counts"
)


umap_coordinates_log2_qn <- prepare_umap_coordinates(
  umap_object = umap_log2_qn,
  processing_label = paste0(
    "log2_count_plus_1_",
    "quantile_normalized"
  )
)


umap_coordinates_combined <- rbind(
  umap_coordinates_raw,
  umap_coordinates_log2_qn
)


# ==============================================================================
# 25. Save UMAP coordinates
# ==============================================================================

write.table(
  umap_coordinates_combined,
  file = umap_coordinates_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ==============================================================================
# 26. Save input summary
# ==============================================================================

write.table(
  input_summary,
  file = input_summary_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ==============================================================================
# 27. Save raw-count UMAP
# ==============================================================================

ggsave(
  filename = raw_pdf_file,
  plot = umap_plot_raw,
  device = "pdf",
  width = 9,
  height = 8,
  units = "in",
  limitsize = FALSE
)


# ==============================================================================
# 28. Save log2 + quantile-normalized UMAP
# ==============================================================================

ggsave(
  filename = log2_qn_pdf_file,
  plot = umap_plot_log2_qn,
  device = "pdf",
  width = 9,
  height = 8,
  units = "in",
  limitsize = FALSE
)


# ==============================================================================
# 29. Save combined two-panel UMAP PDF
# ==============================================================================

grDevices::pdf(
  file = combined_pdf_file,
  width = 17,
  height = 8,
  onefile = TRUE
)


draw_umap_grid()


grDevices::dev.off()


# ==============================================================================
# 30. Display combined UMAP figure
# ==============================================================================

draw_umap_grid()


# ==============================================================================
# 31. Final messages
# ==============================================================================

message(
  "\nCombined two-panel UMAP figure saved to: ",
  normalizePath(
    combined_pdf_file,
    mustWork = FALSE
  )
)


message(
  "Raw UMAP figure saved to: ",
  normalizePath(
    raw_pdf_file,
    mustWork = FALSE
  )
)


message(
  "Log2 + quantile-normalized UMAP figure saved to: ",
  normalizePath(
    log2_qn_pdf_file,
    mustWork = FALSE
  )
)


message(
  "UMAP coordinates saved to: ",
  normalizePath(
    umap_coordinates_file,
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
  "Completed exactly two direct UMAP analyses:"
)


message(
  "  1. raw counts"
)


message(
  "  2. log2(count + 1) + quantile normalization"
)