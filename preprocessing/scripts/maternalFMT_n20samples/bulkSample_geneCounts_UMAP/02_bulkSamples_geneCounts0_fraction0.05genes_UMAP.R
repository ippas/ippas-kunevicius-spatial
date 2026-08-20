# ==============================================================================
# Filter genes, transform data and calculate direct UMAP
# ==============================================================================


# ==============================================================================
# 1. Select the top 50% of genes by mean CPM
# ==============================================================================

library_sizes <- colSums(
  gene_counts_per_sample_raw
)

if (any(!is.finite(library_sizes))) {
  stop("Library sizes contain non-finite values.")
}

if (any(library_sizes <= 0)) {
  stop("At least one sample has zero total counts.")
}


# CPM is calculated only to rank genes by their overall expression

gene_counts_per_sample_cpm <- sweep(
  gene_counts_per_sample_raw,
  MARGIN = 2,
  STATS = library_sizes / 1e6,
  FUN = "/"
)


mean_cpm_per_gene <- rowMeans(
  gene_counts_per_sample_cpm
)


if (anyNA(mean_cpm_per_gene)) {
  stop("Mean CPM values contain missing values.")
}

if (any(!is.finite(mean_cpm_per_gene))) {
  stop("Mean CPM values contain non-finite values.")
}


gene_fraction_to_keep <- 0.05

number_of_genes_before_filtering <- nrow(
  gene_counts_per_sample_raw
)

number_of_genes_to_keep <- ceiling(
  number_of_genes_before_filtering *
    gene_fraction_to_keep
)


gene_ranking <- order(
  mean_cpm_per_gene,
  decreasing = TRUE
)

genes_to_keep <- rownames(
  gene_counts_per_sample_raw
)[
  gene_ranking[
    seq_len(number_of_genes_to_keep)
  ]
]


message(
  "Genes before filtering: ",
  number_of_genes_before_filtering
)

message(
  "Genes retained for UMAP: ",
  length(genes_to_keep)
)

message(
  "Genes removed due to low expression: ",
  number_of_genes_before_filtering -
    length(genes_to_keep)
)

message(
  "Filtering method: top ",
  gene_fraction_to_keep * 100,
  "% of genes ranked by mean CPM across samples"
)


# ==============================================================================
# 2. Create filtered raw bulk gene-count matrix
# ==============================================================================

gene_counts_per_sample_raw_filtered <- gene_counts_per_sample_raw[
  genes_to_keep,
  ,
  drop = FALSE
]


if (
  nrow(gene_counts_per_sample_raw_filtered) !=
    number_of_genes_to_keep
) {
  stop("Incorrect number of genes after filtering.")
}


# ==============================================================================
# 3. Log2 transformation of the same filtered gene subset
# ==============================================================================

gene_counts_per_sample_log2_filtered <- log2(
  gene_counts_per_sample_raw_filtered + 1
)


# ==============================================================================
# 4. Quantile normalization
# ==============================================================================

gene_counts_per_sample_log2_qn_filtered <- quantile_normalize(
  gene_counts_per_sample_log2_filtered
)


# Confirm that both UMAP variants use exactly the same genes and samples

stopifnot(
  identical(
    dim(gene_counts_per_sample_raw_filtered),
    dim(gene_counts_per_sample_log2_qn_filtered)
  )
)

stopifnot(
  identical(
    rownames(gene_counts_per_sample_raw_filtered),
    rownames(gene_counts_per_sample_log2_qn_filtered)
  )
)

stopifnot(
  identical(
    colnames(gene_counts_per_sample_raw_filtered),
    colnames(gene_counts_per_sample_log2_qn_filtered)
  )
)


# ==============================================================================
# 5. Direct UMAP on filtered raw bulk gene counts
# ==============================================================================

umap_raw <- run_direct_umap(
  expression_matrix = gene_counts_per_sample_raw_filtered,
  analysis_id = paste0(
    "raw_bulk_gene_counts_",
    "top50percent_by_mean_CPM_",
    "direct_gene_space"
  ),
  n_neighbors = umap_n_neighbors,
  min_dist = umap_min_dist,
  metric = umap_metric,
  seed = random_seed
)


# ==============================================================================
# 6. Direct UMAP on filtered log2 + quantile-normalized counts
# ==============================================================================

umap_log2_qn <- run_direct_umap(
  expression_matrix = gene_counts_per_sample_log2_qn_filtered,
  analysis_id = paste0(
    "log2_quantile_normalized_",
    "bulk_gene_counts_",
    "top50percent_by_mean_CPM_",
    "direct_gene_space"
  ),
  n_neighbors = umap_n_neighbors,
  min_dist = umap_min_dist,
  metric = umap_metric,
  seed = random_seed
)


# ==============================================================================
# 7. Create UMAP plots
# ==============================================================================

umap_plot_raw <- make_direct_umap_plot(
  umap_object = umap_raw,
  processing_label = paste0(
    "Raw summed gene counts; top 50% of genes by mean CPM; ",
    "no transformation or normalization"
  )
)


umap_plot_log2_qn <- make_direct_umap_plot(
  umap_object = umap_log2_qn,
  processing_label = paste0(
    "Top 50% of genes by mean CPM; ",
    "log2(count + 1) transformation and quantile normalization"
  )
)


# ==============================================================================
# 8. Display both UMAP plots
# ==============================================================================

draw_direct_umap_grid()