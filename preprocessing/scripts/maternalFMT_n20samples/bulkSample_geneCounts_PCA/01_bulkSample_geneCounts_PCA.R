# ==============================================================================
# 01_bulkSample_geneCounts_PCA.R
#
# Purpose:
# - read raw Spatial count matrices for all Maternal FMT samples
# - create bulk sample-level profiles by summing all barcode-level counts
#   for each gene within each sample
# - perform PCA on:
#     1. raw summed gene counts per sample
#     2. log2(count + 1)-transformed and quantile-normalized
#        summed gene counts per sample
# - highlight selected samples
# - save:
#     1. one combined PCA PDF
#     2. one raw-count PCA PDF
#     3. one log2 + quantile-normalized PCA PDF
#     4. one combined PCA-score table
#
# No clustering.
# No integration.
# No low-expression filtering.
# ==============================================================================


# ==============================================================================
# 1. Load packages
# ==============================================================================

suppressPackageStartupMessages({

  library(Seurat)
  library(SeuratObject)

  library(Matrix)

  library(ggplot2)
  library(ggrepel)

})


# ==============================================================================
# 2. Check Seurat version
# ==============================================================================

message(
  "Seurat version: ",
  as.character(packageVersion("Seurat"))
)

message(
  "SeuratObject version: ",
  as.character(packageVersion("SeuratObject"))
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
# 4. Define input and output paths
# ==============================================================================

path_to_data <- paste0(
  "data/",
  "spacerangerCount_withoutJSON_GRCm39-2024-A_2026-03-28/"
)

metadata_file <- "data/metadata_autismFMT.tsv"

output_dir <- paste0(
  "results/",
  "maternalFMT_n20samples/",
  "bulkSample_geneCounts_PCA"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


combined_pdf_file <- file.path(
  output_dir,
  "maternalFMT_bulkSample_geneCounts_PCA_raw_vs_log2Quantile.pdf"
)

raw_pdf_file <- file.path(
  output_dir,
  "maternalFMT_bulkSample_geneCounts_PCA_raw.pdf"
)

log2_qn_pdf_file <- file.path(
  output_dir,
  "maternalFMT_bulkSample_geneCounts_PCA_log2Quantile.pdf"
)

pca_scores_file <- file.path(
  output_dir,
  "maternalFMT_bulkSample_geneCounts_PCA_scores.tsv"
)


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
# 6. Read raw Spatial samples
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
# 7. Check loaded samples
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
# 8. Sum all barcode-level counts per gene within each sample
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


    # Raw matrix:
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


    # Bulk sample-level profile:
    # sum all barcode-level counts for each gene

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
# 9. Check gene order between samples
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
# 10. Create bulk gene x sample raw-count matrix
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
) <- "numeric"


if (anyNA(gene_counts_per_sample_raw)) {

  stop(
    "Raw bulk count matrix contains missing values."
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
# 11. Define quantile-normalization function
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


  if (ncol(expression_matrix) < 2) {

    stop(
      "Quantile normalization requires at least two samples."
    )
  }


  # Sort values independently within each sample

  sorted_matrix <- apply(
    expression_matrix,
    MARGIN = 2,
    FUN = sort
  )


  # Common target distribution

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


    # Handle tied values

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
# 12. Define samples to highlight
# ==============================================================================

red_samples <- c(
  "20_1F",
  "12_3F",
  "15_1M",
  "20_3M"
)

purple_samples <- "2_1M"


sample_colors <- c(
  suspicious = "red",
  reference = "purple",
  remaining = "black"
)


# ==============================================================================
# 13. Define PCA function
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


  # prcomp input:
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
# 14. Define PCA plotting function
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
    "Bulk PCA: all barcode-level counts summed per gene within each sample",
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
        size = 11.5,
        lineheight = 1.1
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
# 15. PCA 1: raw bulk gene counts
# ==============================================================================

pca_raw <- run_sample_pca(
  expression_matrix = gene_counts_per_sample_raw,
  analysis_id = "raw_bulk_gene_counts"
)


# ==============================================================================
# 16. Log2 transformation
# ==============================================================================

gene_counts_per_sample_log2 <- log2(
  gene_counts_per_sample_raw + 1
)


# ==============================================================================
# 17. Quantile normalization
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
# 18. PCA 2: log2-transformed and quantile-normalized bulk gene counts
# ==============================================================================

pca_log2_qn <- run_sample_pca(
  expression_matrix = gene_counts_per_sample_log2_qn,
  analysis_id = "log2_quantile_normalized_bulk_gene_counts"
)


# ==============================================================================
# 19. Create PCA plots
# ==============================================================================

pca_plot_raw <- make_pca_plot(
  pca_object = pca_raw,
  data_processing_label = paste0(
    "Raw bulk gene counts; no transformation or normalization"
  )
)


pca_plot_log2_qn <- make_pca_plot(
  pca_object = pca_log2_qn,
  data_processing_label = paste0(
    "Bulk gene counts after log2(count + 1) transformation ",
    "and quantile normalization"
  )
)


# ==============================================================================
# 20. Define function for drawing combined PCA grid
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
# 21. Prepare PCA-score table
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


  rownames(
    pca_scores
  ) <- NULL


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

  pca_scores$PC1_variance_percent <- pca_object$explained_variance[1]

  pca_scores$PC2_variance_percent <- pca_object$explained_variance[2]


  pc_columns <- grep(
    "^PC[0-9]+$",
    colnames(pca_scores),
    value = TRUE
  )


  pca_scores <- pca_scores[
    ,
    c(
      "analysis",
      "sample_ID",
      "sample_group",
      "genes_used",
      "samples_used",
      "PC1_variance_percent",
      "PC2_variance_percent",
      pc_columns
    ),
    drop = FALSE
  ]


  pca_scores
}


pca_scores_raw <- prepare_pca_scores(
  pca_object = pca_raw,
  analysis_id = "raw_bulk_gene_counts"
)


pca_scores_log2_qn <- prepare_pca_scores(
  pca_object = pca_log2_qn,
  analysis_id = "log2_quantile_normalized_bulk_gene_counts"
)


pca_scores_combined <- rbind(
  pca_scores_raw,
  pca_scores_log2_qn
)


# ==============================================================================
# 22. Save combined PCA-score table
# ==============================================================================

write.table(
  pca_scores_combined,
  file = pca_scores_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ==============================================================================
# 23. Save raw bulk-count PCA as separate PDF
# ==============================================================================

ggsave(
  filename = raw_pdf_file,
  plot = pca_plot_raw,
  device = "pdf",
  width = 9,
  height = 7.5,
  units = "in",
  limitsize = FALSE
)


# ==============================================================================
# 24. Save log2 + quantile-normalized bulk PCA as separate PDF
# ==============================================================================

ggsave(
  filename = log2_qn_pdf_file,
  plot = pca_plot_log2_qn,
  device = "pdf",
  width = 9,
  height = 7.5,
  units = "in",
  limitsize = FALSE
)


# ==============================================================================
# 25. Save both PCA plots in one combined PDF
# ==============================================================================

grDevices::pdf(
  file = combined_pdf_file,
  width = 16,
  height = 7.5,
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
  "Combined PCA figure saved to: ",
  normalizePath(
    combined_pdf_file,
    mustWork = FALSE
  )
)


message(
  "Raw bulk-count PCA figure saved to: ",
  normalizePath(
    raw_pdf_file,
    mustWork = FALSE
  )
)


message(
  "Log2 + quantile-normalized bulk PCA figure saved to: ",
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
  "Bulk sample-level gene-count PCA analysis completed."
)