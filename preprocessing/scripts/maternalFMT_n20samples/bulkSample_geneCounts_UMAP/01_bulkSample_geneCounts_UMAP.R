# ==============================================================================
# 02_bulkSample_geneCounts_directUMAP.R
#
# Purpose:
# - read raw Spatial count matrices for all Maternal FMT samples
# - create bulk sample-level expression profiles by summing all barcode-level
#   counts per gene within each sample
# - run UMAP directly in the complete gene-expression space, without PCA, for:
#     1. raw summed gene counts
#     2. log2(count + 1)-transformed and quantile-normalized summed gene counts
# - highlight selected tissue sections
# - save:
#     1. combined direct UMAP PDF
#     2. raw-count direct UMAP PDF
#     3. log2 + quantile-normalized direct UMAP PDF
#     4. combined direct UMAP-coordinate table
#
# No PCA.
# No clustering.
# No integration.
# No low-expression filtering.
# ==============================================================================


# ==============================================================================
# 1. Check required packages
# ==============================================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "ggplot2",
  "ggrepel",
  "uwot"
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
  library(SeuratObject)

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
  "SeuratObject version: ",
  as.character(packageVersion("SeuratObject"))
)

message(
  "uwot version: ",
  as.character(packageVersion("uwot"))
)


if (packageVersion("Seurat") < "5.0.0") {

  stop(
    "This workflow requires Seurat v5 or newer. ",
    "Current version: ",
    packageVersion("Seurat")
  )
}


# ==============================================================================
# 3. Load custom functions
# ==============================================================================

source(
  "preprocessing/src/functions_prepare_seurat_data.R"
)


# ==============================================================================
# 4. Define input paths
# ==============================================================================

path_to_data <- paste0(
  "data/",
  "spacerangerCount_withoutJSON_GRCm39-2024-A_2026-03-28/"
)

metadata_file <- "data/metadata_autismFMT.tsv"


# ==============================================================================
# 5. Define output paths
# ==============================================================================

output_dir <- paste0(
  "results/",
  "maternalFMT_n20samples/",
  "bulkSample_geneCounts_directUMAP"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


combined_pdf_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_geneCounts_directUMAP_",
    "raw_vs_log2Quantile.pdf"
  )
)


raw_pdf_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_geneCounts_",
    "directUMAP_raw.pdf"
  )
)


log2_qn_pdf_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_geneCounts_",
    "directUMAP_log2Quantile.pdf"
  )
)


umap_coordinates_file <- file.path(
  output_dir,
  paste0(
    "maternalFMT_bulkSample_geneCounts_",
    "directUMAP_coordinates.tsv"
  )
)


# ==============================================================================
# 6. Define direct UMAP settings
# ==============================================================================

umap_n_neighbors <- 5

umap_min_dist <- 0.3

umap_metric <- "euclidean"

random_seed <- 7


# ==============================================================================
# 7. Define samples to highlight
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


if (anyNA(metadata_autismFMT$sample_ID)) {

  stop(
    "Metadata contains missing sample_ID values."
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


message(
  "Samples in metadata: ",
  nrow(metadata_autismFMT)
)


# ==============================================================================
# 9. Read raw Spatial samples
# ==============================================================================

samples_list <- read_spatial_samples(
  path_to_data = path_to_data,
  metadata = metadata_autismFMT,
  sample_id_col = "sample_ID",
  min.cells = 0,
  min.features = 0,
  verbose = TRUE
)


message(
  "Loaded samples: ",
  length(samples_list)
)


# ==============================================================================
# 10. Check loaded samples
# ==============================================================================

expected_samples <- metadata_autismFMT$sample_ID

loaded_samples <- names(samples_list)


missing_samples <- setdiff(
  expected_samples,
  loaded_samples
)


unexpected_samples <- setdiff(
  loaded_samples,
  expected_samples
)


if (length(missing_samples) > 0) {

  stop(
    "Missing samples: ",
    paste(
      missing_samples,
      collapse = ", "
    )
  )
}


if (length(unexpected_samples) > 0) {

  stop(
    "Unexpected samples: ",
    paste(
      unexpected_samples,
      collapse = ", "
    )
  )
}


# Arrange samples according to metadata order

samples_list <- samples_list[
  expected_samples
]


# ==============================================================================
# 11. Sum all barcode-level counts per gene within each sample
# ==============================================================================

sample_gene_counts <- lapply(
  names(samples_list),
  function(sample_id) {

    message(
      "Processing sample: ",
      sample_id
    )


    seurat_object <- samples_list[[sample_id]]


    if (!"RNA" %in% Assays(seurat_object)) {

      stop(
        "RNA assay is missing in sample: ",
        sample_id
      )
    }


    # Raw count matrix:
    # rows    = genes
    # columns = barcodes/spots

    counts_matrix <- GetAssayData(
      object = seurat_object,
      assay = "RNA",
      layer = "counts"
    )


    if (nrow(counts_matrix) == 0) {

      stop(
        "No genes found in sample: ",
        sample_id
      )
    }


    if (ncol(counts_matrix) == 0) {

      stop(
        "No barcodes found in sample: ",
        sample_id
      )
    }


    # Bulk sample-level gene expression:
    # sum all barcode-level counts for every gene

    gene_counts <- Matrix::rowSums(
      counts_matrix
    )


    gene_counts <- setNames(
      as.numeric(gene_counts),
      rownames(counts_matrix)
    )


    message(
      "  genes: ",
      nrow(counts_matrix),
      "; barcodes: ",
      ncol(counts_matrix),
      "; total summed counts: ",
      format(
        sum(gene_counts),
        big.mark = ",",
        scientific = FALSE
      )
    )


    gene_counts
  }
)


names(sample_gene_counts) <- names(samples_list)


# ==============================================================================
# 12. Check gene names and gene order
# ==============================================================================

reference_genes <- names(
  sample_gene_counts[[1]]
)


identical_gene_order <- vapply(
  sample_gene_counts,
  function(x) {

    identical(
      names(x),
      reference_genes
    )
  },
  logical(1)
)


if (!all(identical_gene_order)) {

  problematic_samples <- names(
    identical_gene_order
  )[!identical_gene_order]

  stop(
    "Gene names or gene order differ in samples: ",
    paste(
      problematic_samples,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 13. Create bulk gene x sample raw-count matrix
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
) <- names(sample_gene_counts)


storage.mode(
  gene_counts_per_sample_raw
) <- "double"


if (anyNA(gene_counts_per_sample_raw)) {

  stop(
    "Raw bulk count matrix contains missing values."
  )
}


if (any(!is.finite(gene_counts_per_sample_raw))) {

  stop(
    "Raw bulk count matrix contains non-finite values."
  )
}


if (any(gene_counts_per_sample_raw < 0)) {

  stop(
    "Raw bulk count matrix contains negative values."
  )
}


message(
  "Bulk raw-count matrix dimensions: ",
  nrow(gene_counts_per_sample_raw),
  " genes x ",
  ncol(gene_counts_per_sample_raw),
  " samples"
)


# ==============================================================================
# 14. Define quantile-normalization function
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


  # Sort values independently within every sample

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


    # Preserve ties by assigning identical normalized values
    # to identical original values

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
# 15. Log2 transformation
# ==============================================================================

gene_counts_per_sample_log2 <- log2(
  gene_counts_per_sample_raw + 1
)


# ==============================================================================
# 16. Quantile normalization
# ==============================================================================

gene_counts_per_sample_log2_qn <- quantile_normalize(
  gene_counts_per_sample_log2
)


# Confirm identical matrix structure

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
# 17. Define direct UMAP function
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
  # UMAP requires:
  # rows    = observations/samples
  # columns = features/genes

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
    analysis_id,
    ": direct UMAP using ",
    number_of_samples,
    " samples and ",
    number_of_genes,
    " genes"
  )


  message(
    analysis_id,
    ": n_neighbors = ",
    n_neighbors_used,
    "; min_dist = ",
    min_dist,
    "; metric = ",
    metric,
    "; PCA = disabled"
  )


  set.seed(
    seed
  )


  umap_matrix <- uwot::umap(
    X = sample_by_gene_matrix,

    n_neighbors = n_neighbors_used,
    n_components = 2,

    metric = metric,
    min_dist = min_dist,

    # No feature scaling inside UMAP

    scale = FALSE,

    # No PCA preprocessing

    pca = NULL,

    # Random initialization prevents PCA-based initialization

    init = "random",

    # Deterministic single-thread settings

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
# 18. Direct UMAP on raw bulk gene counts
# ==============================================================================

umap_raw <- run_direct_umap(
  expression_matrix = gene_counts_per_sample_raw,
  analysis_id = "raw_bulk_gene_counts_direct_gene_space",
  n_neighbors = umap_n_neighbors,
  min_dist = umap_min_dist,
  metric = umap_metric,
  seed = random_seed
)


# ==============================================================================
# 19. Direct UMAP on log2 + quantile-normalized bulk gene counts
# ==============================================================================

umap_log2_qn <- run_direct_umap(
  expression_matrix = gene_counts_per_sample_log2_qn,
  analysis_id = paste0(
    "log2_quantile_normalized_",
    "bulk_gene_counts_direct_gene_space"
  ),
  n_neighbors = umap_n_neighbors,
  min_dist = umap_min_dist,
  metric = umap_metric,
  seed = random_seed
)


# ==============================================================================
# 20. Define direct UMAP plotting function
# ==============================================================================

make_direct_umap_plot <- function(
    umap_object,
    processing_label
) {

  complete_plot_title <- paste0(
    "Direct bulk UMAP: all barcode-level counts summed ",
    "per gene within each sample",
    "\n",
    processing_label,
    "\n",
    umap_object$genes_used,
    " genes | ",
    umap_object$samples_used,
    " samples | direct gene space; no PCA"
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
        size = 11.5,
        lineheight = 1.1
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
# 21. Create direct UMAP plots
# ==============================================================================

umap_plot_raw <- make_direct_umap_plot(
  umap_object = umap_raw,
  processing_label = paste0(
    "Raw summed gene counts; ",
    "no transformation or normalization"
  )
)


umap_plot_log2_qn <- make_direct_umap_plot(
  umap_object = umap_log2_qn,
  processing_label = paste0(
    "Summed gene counts after log2(count + 1) transformation ",
    "and quantile normalization"
  )
)


# ==============================================================================
# 22. Define function for drawing combined direct UMAP grid
# ==============================================================================

draw_direct_umap_grid <- function() {

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
# 23. Prepare direct UMAP-coordinate table
# ==============================================================================

prepare_umap_coordinates <- function(
    umap_object,
    processing_label
) {

  umap_coordinates <- umap_object$umap_df


  umap_coordinates$analysis <- umap_object$analysis_id

  umap_coordinates$processing <- processing_label

  umap_coordinates$input_space <- "all_genes_direct_no_PCA"

  umap_coordinates$genes_used <- umap_object$genes_used

  umap_coordinates$samples_used <- umap_object$samples_used

  umap_coordinates$n_neighbors <- umap_object$n_neighbors

  umap_coordinates$min_dist <- umap_object$min_dist

  umap_coordinates$metric <- umap_object$metric

  umap_coordinates$seed <- umap_object$seed


  umap_coordinates <- umap_coordinates[
    ,
    c(
      "analysis",
      "processing",
      "input_space",
      "sample_ID",
      "sample_group",
      "genes_used",
      "samples_used",
      "n_neighbors",
      "min_dist",
      "metric",
      "seed",
      "UMAP1",
      "UMAP2"
    ),
    drop = FALSE
  ]


  umap_coordinates
}


umap_coordinates_raw <- prepare_umap_coordinates(
  umap_object = umap_raw,
  processing_label = "raw_bulk_gene_counts"
)


umap_coordinates_log2_qn <- prepare_umap_coordinates(
  umap_object = umap_log2_qn,
  processing_label = paste0(
    "log2_count_plus_1_",
    "quantile_normalized_bulk_gene_counts"
  )
)


umap_coordinates_combined <- rbind(
  umap_coordinates_raw,
  umap_coordinates_log2_qn
)


# ==============================================================================
# 24. Save direct UMAP coordinates
# ==============================================================================

write.table(
  umap_coordinates_combined,
  file = umap_coordinates_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ==============================================================================
# 25. Save raw direct UMAP as separate PDF
# ==============================================================================

ggsave(
  filename = raw_pdf_file,
  plot = umap_plot_raw,
  device = "pdf",
  width = 9,
  height = 7.5,
  units = "in",
  limitsize = FALSE
)


# ==============================================================================
# 26. Save log2 + quantile-normalized direct UMAP as separate PDF
# ==============================================================================

ggsave(
  filename = log2_qn_pdf_file,
  plot = umap_plot_log2_qn,
  device = "pdf",
  width = 9,
  height = 7.5,
  units = "in",
  limitsize = FALSE
)


# ==============================================================================
# 27. Save combined direct UMAP PDF
# ==============================================================================

grDevices::pdf(
  file = combined_pdf_file,
  width = 16,
  height = 7.5,
  onefile = TRUE
)


draw_direct_umap_grid()


grDevices::dev.off()


# ==============================================================================
# 28. Display combined direct UMAP figure
# ==============================================================================

draw_direct_umap_grid()


# ==============================================================================
# 29. Final messages
# ==============================================================================

message(
  "Combined direct UMAP figure saved to: ",
  normalizePath(
    combined_pdf_file,
    mustWork = FALSE
  )
)


message(
  "Raw direct UMAP figure saved to: ",
  normalizePath(
    raw_pdf_file,
    mustWork = FALSE
  )
)


message(
  "Log2 + quantile-normalized direct UMAP figure saved to: ",
  normalizePath(
    log2_qn_pdf_file,
    mustWork = FALSE
  )
)


message(
  "Direct UMAP coordinates saved to: ",
  normalizePath(
    umap_coordinates_file,
    mustWork = FALSE
  )
)


message(
  "Direct bulk gene-space UMAP completed without PCA."
)