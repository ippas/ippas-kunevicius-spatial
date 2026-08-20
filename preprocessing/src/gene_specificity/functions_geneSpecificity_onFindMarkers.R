# ==============================================================================
# functions_geneSpecificity_onFindMarkers.R
#
# Functions for calculating cluster specificity for genes identified by
# Seurat FindAllMarkers().
# ==============================================================================


sort_cluster_ids <- function(cluster_ids) {

  cluster_ids <- unique(
    as.character(cluster_ids)
  )

  numeric_ids <- suppressWarnings(
    as.numeric(cluster_ids)
  )

  if (!anyNA(numeric_ids)) {
    return(
      cluster_ids[order(numeric_ids)]
    )
  }

  sort(cluster_ids)
}


convert_log1p_to_linear <- function(expression_matrix) {

  if (inherits(expression_matrix, "sparseMatrix")) {
    expression_matrix@x <- expm1(expression_matrix@x)
    return(expression_matrix)
  }

  expm1(expression_matrix)
}


calculate_mean_expression_by_cluster <- function(
    expression_matrix,
    cluster_ids
) {

  if (ncol(expression_matrix) != length(cluster_ids)) {
    stop(
      "Number of expression-matrix columns and cluster IDs differs.",
      call. = FALSE
    )
  }

  cluster_levels <- sort_cluster_ids(cluster_ids)

  cluster_factor <- factor(
    as.character(cluster_ids),
    levels = cluster_levels
  )

  cluster_design <- Matrix::sparseMatrix(
    i = seq_along(cluster_factor),
    j = as.integer(cluster_factor),
    x = 1,
    dims = c(
      length(cluster_factor),
      length(cluster_levels)
    ),
    dimnames = list(
      colnames(expression_matrix),
      cluster_levels
    )
  )

  cluster_sizes <- Matrix::colSums(
    cluster_design
  )

  mean_expression <- expression_matrix %*%
    cluster_design

  mean_expression <- sweep(
    as.matrix(mean_expression),
    MARGIN = 2,
    STATS = cluster_sizes,
    FUN = "/"
  )

  mean_expression
}


calculate_tau <- function(expression_values) {

  number_of_clusters <- length(expression_values)
  maximum_expression <- max(expression_values)

  if (number_of_clusters < 2L ||
      maximum_expression <= 0) {
    return(NA_real_)
  }

  sum(
    1 - expression_values / maximum_expression
  ) / (
    number_of_clusters - 1
  )
}


calculate_gini <- function(expression_values) {

  number_of_clusters <- length(expression_values)
  total_expression <- sum(expression_values)

  if (number_of_clusters < 2L ||
      total_expression <= 0) {
    return(NA_real_)
  }

  sorted_expression <- sort(
    expression_values
  )

  gini <- sum(
    (
      2 * seq_len(number_of_clusters) -
        number_of_clusters -
        1
    ) *
      sorted_expression
  ) / (
    number_of_clusters *
      total_expression
  )

  # Normalize the maximum possible value to 1.
  gini * number_of_clusters /
    (number_of_clusters - 1)
}


calculate_shannon_specificity <- function(
    expression_values
) {

  number_of_clusters <- length(expression_values)
  total_expression <- sum(expression_values)

  if (number_of_clusters < 2L ||
      total_expression <= 0) {
    return(NA_real_)
  }

  probabilities <- expression_values /
    total_expression

  probabilities <- probabilities[
    probabilities > 0
  ]

  shannon_entropy <- -sum(
    probabilities * log(probabilities)
  )

  1 - shannon_entropy /
    log(number_of_clusters)
}


add_specificity_metrics_to_findmarkers <- function(
    marker_table,
    mean_expression_matrix,
    gene_column = "gene",
    cluster_column = "cluster"
) {

  marker_table[[gene_column]] <- as.character(
    marker_table[[gene_column]]
  )

  marker_table[[cluster_column]] <- as.character(
    marker_table[[cluster_column]]
  )

  specificity_rows <- lapply(
    seq_len(nrow(marker_table)),
    function(row_index) {

      gene <- marker_table[
        row_index,
        gene_column
      ]

      target_cluster <- marker_table[
        row_index,
        cluster_column
      ]

      expression_values <- as.numeric(
        mean_expression_matrix[
          gene,
          ,
          drop = TRUE
        ]
      )

      target_index <- match(
        target_cluster,
        colnames(mean_expression_matrix)
      )

      other_indices <- setdiff(
        seq_along(expression_values),
        target_index
      )

      best_other_index <- other_indices[
        which.max(
          expression_values[other_indices]
        )
      ]

      target_expression <- expression_values[
        target_index
      ]

      best_other_expression <- expression_values[
        best_other_index
      ]

      total_expression <- sum(
        expression_values
      )

      expression_ratio <- if (
        best_other_expression > 0
      ) {
        target_expression /
          best_other_expression
      } else if (
        target_expression > 0
      ) {
        Inf
      } else {
        NA_real_
      }

      data.frame(
        best_other_cluster =
          colnames(mean_expression_matrix)[
            best_other_index
          ],

        mean_expression_best_other =
          best_other_expression,

        expression_specificity =
          if (total_expression > 0) {
            target_expression /
              total_expression
          } else {
            NA_real_
          },

        expression_ratio_vs_best_other =
          expression_ratio,

        is_best_cluster =
          target_expression ==
          max(expression_values),

        tau =
          calculate_tau(
            expression_values
          ),

        gini =
          calculate_gini(
            expression_values
          ),

        shannon_specificity =
          calculate_shannon_specificity(
            expression_values
          )
      )
    }
  )

  specificity_table <- do.call(
    rbind,
    specificity_rows
  )

  rownames(specificity_table) <- NULL

  cbind(
    marker_table,
    specificity_table
  )
}


calculate_gene_specificity_on_findmarkers <- function(
    seurat_object,
    marker_table,
    cluster_column,
    assay_name = "RNA"
) {

  required_columns <- c(
    "gene",
    "cluster"
  )

  if (!all(
    required_columns %in% colnames(marker_table)
  )) {
    stop(
      "Marker table must contain gene and cluster columns.",
      call. = FALSE
    )
  }

  DefaultAssay(seurat_object) <- assay_name

  data_layers <- Layers(
    seurat_object[[assay_name]],
    search = "^data"
  )

  if (length(data_layers) > 1L) {
    seurat_object <- JoinLayers(
      object = seurat_object,
      assay = assay_name
    )
  }

  marker_genes <- unique(
    as.character(marker_table$gene)
  )

  missing_genes <- setdiff(
    marker_genes,
    rownames(seurat_object[[assay_name]])
  )

  if (length(missing_genes) > 0L) {
    stop(
      "Marker genes missing from RNA assay: ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  expression_log <- LayerData(
    object = seurat_object,
    assay = assay_name,
    layer = "data",
    features = marker_genes
  )

  cluster_ids <- as.character(
    seurat_object[[cluster_column]][
      colnames(expression_log),
      1
    ]
  )

  expression_linear <- convert_log1p_to_linear(
    expression_log
  )

  mean_expression_matrix <-
    calculate_mean_expression_by_cluster(
      expression_matrix = expression_linear,
      cluster_ids = cluster_ids
    )

  add_specificity_metrics_to_findmarkers(
    marker_table = marker_table,
    mean_expression_matrix =
      mean_expression_matrix
  )
}