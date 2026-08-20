# ==============================================================================
# Add Seurat clusters across multiple resolutions
#
# Behaviour:
# - reuses an existing SNN graph unless force_neighbors = TRUE
# - skips existing cluster columns unless force_clusters = TRUE
# - stores every resolution in metadata:
#   clusters_res0.01, clusters_res0.02, ..., clusters_res2
# - orders all clusters_res* metadata columns from lowest to highest resolution
# ==============================================================================

add_clusters_multi_resolution <- function(
    seurat_object,
    reduction = "integrated.rpca",
    dims = 1:30,
    resolution_start = 0.05,
    resolution_end = 2,
    resolution_step = 0.05,
    graph_names = c(
      "integrated_rpca_nn",
      "integrated_rpca_snn"
    ),
    cluster_prefix = "clusters_res",
    k.param = 20,
    algorithm = 1,
    random.seed = 7,
    force_neighbors = FALSE,
    force_clusters = FALSE,
    set_active_resolution = NULL,
    verbose = TRUE
) {

  # ============================================================================
  # 1. Validate input
  # ============================================================================

  if (!inherits(seurat_object, "Seurat")) {
    stop("`seurat_object` must be a Seurat object.")
  }

  if (!reduction %in% SeuratObject::Reductions(seurat_object)) {
    stop(
      "Reduction `", reduction, "` is absent from the object.\n",
      "Available reductions: ",
      paste(SeuratObject::Reductions(seurat_object), collapse = ", ")
    )
  }

  if (!is.numeric(dims) || length(dims) == 0) {
    stop("`dims` must be a non-empty numeric vector.")
  }

  available_dims <- ncol(
    SeuratObject::Embeddings(seurat_object[[reduction]])
  )

  if (max(dims) > available_dims) {
    stop(
      "Requested dimension ", max(dims),
      " exceeds available dimensions: ", available_dims
    )
  }

  if (length(graph_names) != 2 || any(graph_names == "")) {
    stop(
      "`graph_names` must contain exactly two names:\n",
      "1. nearest-neighbour graph\n",
      "2. SNN graph"
    )
  }

  if (anyDuplicated(graph_names) > 0) {
    stop("`graph_names` must contain two unique names.")
  }

  if (resolution_start <= 0 ||
      resolution_end <= 0 ||
      resolution_step <= 0 ||
      resolution_end < resolution_start) {
    stop("Invalid resolution range.")
  }

  if (!is.logical(force_neighbors) || length(force_neighbors) != 1) {
    stop("`force_neighbors` must be TRUE or FALSE.")
  }

  if (!is.logical(force_clusters) || length(force_clusters) != 1) {
    stop("`force_clusters` must be TRUE or FALSE.")
  }

  # ============================================================================
  # 2. Helper functions
  # ============================================================================

  format_resolution <- function(x) {
    x <- sprintf("%.2f", x)
    x <- sub("0+$", "", x)
    x <- sub("\\.$", "", x)
    return(x)
  }

  reorder_cluster_columns <- function(object, prefix) {

    metadata_columns <- colnames(object@meta.data)

    candidate_columns <- metadata_columns[
      startsWith(metadata_columns, prefix)
    ]

    if (length(candidate_columns) == 0) {
      return(object)
    }

    resolution_text <- substring(
      candidate_columns,
      first = nchar(prefix) + 1
    )

    resolution_numeric <- suppressWarnings(
      as.numeric(resolution_text)
    )

    sortable_columns <- candidate_columns[
      !is.na(resolution_numeric)
    ]

    sortable_resolutions <- resolution_numeric[
      !is.na(resolution_numeric)
    ]

    ordered_cluster_columns <- sortable_columns[
      order(sortable_resolutions, sortable_columns)
    ]

    non_cluster_columns <- metadata_columns[
      !metadata_columns %in% sortable_columns
    ]

    object@meta.data <- object@meta.data[
      ,
      c(non_cluster_columns, ordered_cluster_columns),
      drop = FALSE
    ]

    return(object)
  }

  # ============================================================================
  # 3. Define requested clustering resolutions
  # ============================================================================

  resolutions <- seq(
    from = resolution_start,
    to = resolution_end,
    by = resolution_step
  )

  resolutions <- round(resolutions, 8)

  cluster_columns <- paste0(
    cluster_prefix,
    vapply(
      resolutions,
      format_resolution,
      character(1)
    )
  )

  metadata_columns <- colnames(seurat_object[[]])

  cluster_exists <- cluster_columns %in% metadata_columns

  if (force_clusters) {
    run_clusters <- rep(TRUE, length(cluster_columns))
  } else {
    run_clusters <- !cluster_exists
  }

  resolutions_to_run <- resolutions[run_clusters]
  columns_to_run <- cluster_columns[run_clusters]

  if (verbose && any(cluster_exists) && !force_clusters) {
    message(
      "Keeping existing clustering column(s):\n",
      paste(cluster_columns[cluster_exists], collapse = ", ")
    )
  }

  # ============================================================================
  # 4. Preserve previous identities and generic cluster column
  # ============================================================================

  original_idents <- SeuratObject::Idents(seurat_object)

  had_seurat_clusters <- "seurat_clusters" %in% metadata_columns

  if (had_seurat_clusters) {
    original_seurat_clusters <- seurat_object$seurat_clusters
  }

  # ============================================================================
  # 5. Build/reuse the SNN graph
  # ============================================================================

  if (length(resolutions_to_run) > 0) {

    existing_graphs <- SeuratObject::Graphs(seurat_object)

    snn_graph_exists <- graph_names[2] %in% existing_graphs

    build_neighbors <- force_neighbors || !snn_graph_exists

    if (build_neighbors) {

      if (verbose) {
        message(
          "Building nearest-neighbour and SNN graphs from `",
          reduction, "`..."
        )
      }

      seurat_object <- Seurat::FindNeighbors(
        object = seurat_object,
        reduction = reduction,
        dims = dims,
        k.param = k.param,
        graph.name = graph_names,
        verbose = verbose
      )

    } else if (verbose) {

      message(
        "Reusing existing SNN graph: ",
        graph_names[2]
      )
    }
  }

  # ============================================================================
  # 6. Run requested resolutions
  # ============================================================================

  if (length(resolutions_to_run) > 0) {

    for (i in seq_along(resolutions_to_run)) {

      current_resolution <- resolutions_to_run[i]
      current_column <- columns_to_run[i]

      if (verbose) {
        message(
          "Clustering resolution ",
          format_resolution(current_resolution),
          " -> ",
          current_column
        )
      }

      seurat_object <- Seurat::FindClusters(
        object = seurat_object,
        graph.name = graph_names[2],
        resolution = current_resolution,
        algorithm = algorithm,
        random.seed = random.seed,
        cluster.name = current_column,
        verbose = verbose
      )
    }

  } else if (verbose) {

    message(
      "All requested clustering columns already exist. ",
      "No clustering was rerun."
    )
  }

  # ============================================================================
  # 7. Restore previous identities
  # ============================================================================

  if (had_seurat_clusters) {
    seurat_object$seurat_clusters <- original_seurat_clusters
  } else if ("seurat_clusters" %in% colnames(seurat_object[[]])) {
    seurat_object$seurat_clusters <- NULL
  }

  if (is.null(set_active_resolution)) {

    SeuratObject::Idents(seurat_object) <- original_idents

  } else {

    active_column <- paste0(
      cluster_prefix,
      format_resolution(set_active_resolution)
    )

    if (!active_column %in% colnames(seurat_object[[]])) {
      stop(
        "Requested active clustering column does not exist: ",
        active_column
      )
    }

    SeuratObject::Idents(seurat_object) <- active_column
  }

  # ============================================================================
  # 8. Reorder cluster metadata columns
  # ============================================================================

  seurat_object <- reorder_cluster_columns(
    object = seurat_object,
    prefix = cluster_prefix
  )

  # ============================================================================
  # 9. Store a compact summary in misc
  # ============================================================================

  final_metadata <- seurat_object[[]]

  all_cluster_columns <- colnames(final_metadata)[
    startsWith(colnames(final_metadata), cluster_prefix)
  ]

  all_cluster_resolutions <- suppressWarnings(
    as.numeric(
      substring(
        all_cluster_columns,
        first = nchar(cluster_prefix) + 1
      )
    )
  )

  valid_cluster_columns <- all_cluster_columns[
    !is.na(all_cluster_resolutions)
  ]

  valid_cluster_resolutions <- all_cluster_resolutions[
    !is.na(all_cluster_resolutions)
  ]

  cluster_summary <- data.frame(
    resolution = valid_cluster_resolutions,
    cluster_column = valid_cluster_columns,
    n_clusters = vapply(
      valid_cluster_columns,
      function(column_name) {
        length(
          unique(
            stats::na.omit(
              final_metadata[[column_name]]
            )
          )
        )
      },
      integer(1)
    ),
    stringsAsFactors = FALSE
  )

  cluster_summary <- cluster_summary[
    order(cluster_summary$resolution),
    ,
    drop = FALSE
  ]

  seurat_object@misc$cluster_resolution_summary <- cluster_summary

  return(seurat_object)
}

