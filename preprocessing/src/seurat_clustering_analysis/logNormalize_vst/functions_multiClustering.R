# ==============================================================================
# functions_multiClustering.R
#
# Purpose:
# Run multiple Seurat clustering algorithms across multiple resolutions on an
# existing SNN graph.
#
# Important behaviour:
# - clustering is performed directly on the SNN Graph object;
# - the input Seurat object is not modified during individual clustering runs;
# - no future parallelization is used;
# - all requested runs must finish successfully;
# - cluster assignments are validated against Seurat spot names;
# - all cluster columns are added to the Seurat object once, at the very end;
# - the generic `seurat_clusters` column is not created or overwritten;
# - existing active identities are preserved unless explicitly changed.
# ==============================================================================


format_clustering_resolution <- function(
    resolution,
    digits = 2L
) {

  formatC(
    as.numeric(resolution),
    format = "f",
    digits = as.integer(digits)
  )
}


format_clustering_elapsed_time <- function(seconds) {

  seconds <- max(
    0,
    as.numeric(seconds)
  )

  hours <- floor(seconds / 3600)
  minutes <- floor((seconds %% 3600) / 60)
  remaining_seconds <- floor(seconds %% 60)

  sprintf(
    "%02d:%02d:%02d",
    as.integer(hours),
    as.integer(minutes),
    as.integer(remaining_seconds)
  )
}


multi_clustering_message <- function(...) {

  cat(
    "[",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    "] ",
    paste0(..., collapse = ""),
    "\n",
    sep = ""
  )

  flush.console()
  invisible(NULL)
}


extract_single_clustering_result <- function(
    clustering_result,
    expected_spots
) {

  if (
    is.data.frame(clustering_result) ||
      is.matrix(clustering_result)
  ) {

    if (ncol(clustering_result) != 1L) {
      stop(
        "FindClusters returned ",
        ncol(clustering_result),
        " columns; exactly one column was expected.",
        call. = FALSE
      )
    }

    clusters <- clustering_result[[1]]
    names(clusters) <- rownames(clustering_result)

  } else if (
    is.factor(clustering_result) ||
      is.atomic(clustering_result)
  ) {

    clusters <- clustering_result

  } else {

    stop(
      "Unsupported FindClusters result class: ",
      paste(
        class(clustering_result),
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  if (is.null(names(clusters))) {

    if (length(clusters) != length(expected_spots)) {
      stop(
        "The clustering result has no spot names and its length does not ",
        "match the number of spots in the Seurat object.",
        call. = FALSE
      )
    }

    names(clusters) <- expected_spots
  }

  missing_spots <- setdiff(
    expected_spots,
    names(clusters)
  )

  extra_spots <- setdiff(
    names(clusters),
    expected_spots
  )

  if (
    length(missing_spots) > 0L ||
      length(extra_spots) > 0L
  ) {
    stop(
      "Spot names returned by FindClusters do not match the Seurat object. ",
      "Missing spots: ",
      length(missing_spots),
      "; extra spots: ",
      length(extra_spots),
      ".",
      call. = FALSE
    )
  }

  clusters <- clusters[expected_spots]

  if (anyNA(clusters)) {
    stop(
      "NA values were found in cluster assignments.",
      call. = FALSE
    )
  }

  cluster_labels <- as.character(clusters)

  if (all(grepl("^-?[0-9]+$", cluster_labels))) {

    numeric_levels <- sort(
      unique(
        as.integer(cluster_labels)
      )
    )

    return(
      factor(
        as.integer(cluster_labels),
        levels = numeric_levels
      )
    )
  }

  factor(
    cluster_labels,
    levels = sort(unique(cluster_labels))
  )
}


run_multi_clustering <- function(
    seurat_object,
    graph_name,
    algorithms = c(
      "louvain",
      "louvainRefined",
      "slm",
      "leiden"
    ),
    resolutions = c(
      0.20,
      0.30,
      0.40
    ),
    resolution_digits = 2L,
    modularity_fxn = 1L,
    n_start = 10L,
    n_iter = 10L,
    random_seed = 7L,
    group_singletons = TRUE,
    leiden_method = "leidenbase",
    leiden_objective_function = "modularity",
    overwrite_existing = FALSE,
    set_active_column = NULL,
    verbose = TRUE,
    verbose_findclusters = TRUE
) {

  # ============================================================================
  # 1. Validate the Seurat object and SNN graph
  # ============================================================================

  if (!inherits(seurat_object, "Seurat")) {
    stop(
      "`seurat_object` must be a Seurat object.",
      call. = FALSE
    )
  }

  if (
    !is.character(graph_name) ||
      length(graph_name) != 1L ||
      is.na(graph_name) ||
      graph_name == ""
  ) {
    stop(
      "`graph_name` must be one non-empty character value.",
      call. = FALSE
    )
  }

  existing_graphs <- SeuratObject::Graphs(
    seurat_object
  )

  if (!graph_name %in% existing_graphs) {
    stop(
      "Graph `",
      graph_name,
      "` is absent from the Seurat object.\nAvailable graphs: ",
      paste(existing_graphs, collapse = ", "),
      call. = FALSE
    )
  }

  graph_object <- seurat_object[[graph_name]]

  if (!inherits(graph_object, "Graph")) {
    stop(
      "Object `",
      graph_name,
      "` is not a Seurat Graph.",
      call. = FALSE
    )
  }

  spot_names <- colnames(seurat_object)

  if (length(spot_names) == 0L) {
    stop(
      "The Seurat object has no spot names.",
      call. = FALSE
    )
  }

  if (anyDuplicated(spot_names) > 0L) {
    stop(
      "Duplicated spot names were found in the Seurat object.",
      call. = FALSE
    )
  }

  if (
    nrow(graph_object) != length(spot_names) ||
      ncol(graph_object) != length(spot_names)
  ) {
    stop(
      "Graph dimensions do not match the number of Seurat spots.",
      call. = FALSE
    )
  }

  if (
    !identical(rownames(graph_object), spot_names) ||
      !identical(colnames(graph_object), spot_names)
  ) {
    stop(
      "Graph row/column names are not identically ordered as Seurat spots.",
      call. = FALSE
    )
  }


  # ============================================================================
  # 2. Validate algorithms and resolutions
  # ============================================================================

  algorithm_numbers <- c(
    louvain = 1L,
    louvainRefined = 2L,
    slm = 3L,
    leiden = 4L
  )

  if (
    !is.character(algorithms) ||
      length(algorithms) == 0L ||
      anyNA(algorithms)
  ) {
    stop(
      "`algorithms` must be a non-empty character vector.",
      call. = FALSE
    )
  }

  unknown_algorithms <- setdiff(
    algorithms,
    names(algorithm_numbers)
  )

  if (length(unknown_algorithms) > 0L) {
    stop(
      "Unsupported clustering algorithm(s): ",
      paste(
        unknown_algorithms,
        collapse = ", "
      ),
      ". Supported algorithms: ",
      paste(
        names(algorithm_numbers),
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  if (anyDuplicated(algorithms) > 0L) {
    stop(
      "`algorithms` contains duplicated values.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(resolutions) ||
      length(resolutions) == 0L ||
      anyNA(resolutions) ||
      any(!is.finite(resolutions)) ||
      any(resolutions <= 0)
  ) {
    stop(
      "`resolutions` must be a non-empty numeric vector containing ",
      "finite values greater than zero.",
      call. = FALSE
    )
  }

  resolutions <- sort(
    as.numeric(resolutions)
  )

  if (anyDuplicated(resolutions) > 0L) {
    stop(
      "`resolutions` contains duplicated values.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(resolution_digits) ||
      length(resolution_digits) != 1L ||
      resolution_digits < 1L
  ) {
    stop(
      "`resolution_digits` must be one positive integer.",
      call. = FALSE
    )
  }

  resolution_labels <- vapply(
    resolutions,
    format_clustering_resolution,
    character(1),
    digits = as.integer(resolution_digits)
  )

  cluster_columns <- unlist(
    lapply(
      algorithms,
      function(algorithm_name) {
        paste0(
          algorithm_name,
          "_res",
          resolution_labels
        )
      }
    ),
    use.names = FALSE
  )

  if (anyDuplicated(cluster_columns) > 0L) {
    stop(
      "Generated cluster column names are not unique.",
      call. = FALSE
    )
  }

  existing_cluster_columns <- intersect(
    cluster_columns,
    colnames(seurat_object[[]])
  )

  if (
    length(existing_cluster_columns) > 0L &&
      !isTRUE(overwrite_existing)
  ) {
    stop(
      "The following cluster columns already exist:\n",
      paste(
        existing_cluster_columns,
        collapse = ", "
      ),
      "\nSet `overwrite_existing = TRUE` to replace them.",
      call. = FALSE
    )
  }

  if ("leiden" %in% algorithms) {

    if (
      !leiden_method %in% c(
        "leidenbase",
        "igraph"
      )
    ) {
      stop(
        "`leiden_method` must be `leidenbase` or `igraph`.",
        call. = FALSE
      )
    }

    if (
      leiden_method == "leidenbase" &&
        !requireNamespace(
          "leidenbase",
          quietly = TRUE
        )
    ) {
      stop(
        "Package `leidenbase` is required for Leiden clustering.",
        call. = FALSE
      )
    }

    if (
      leiden_method == "igraph" &&
        !requireNamespace(
          "igraph",
          quietly = TRUE
        )
    ) {
      stop(
        "Package `igraph` is required for igraph-based Leiden clustering.",
        call. = FALSE
      )
    }
  }


  # ============================================================================
  # 3. Run all algorithm-resolution combinations directly on the graph
  #
  # The Seurat object is not modified inside these loops.
  # ============================================================================

  total_runs <- length(algorithms) * length(resolutions)
  current_run <- 0L
  workflow_started <- Sys.time()

  cluster_assignments <- data.frame(
    row.names = spot_names,
    check.names = FALSE
  )

  summary_rows <- vector(
    mode = "list",
    length = total_runs
  )

  if (verbose) {
    multi_clustering_message(
      "Starting clustering on graph `",
      graph_name,
      "`: ",
      length(algorithms),
      " algorithms x ",
      length(resolutions),
      " resolutions = ",
      total_runs,
      " runs."
    )
  }

  for (algorithm_name in algorithms) {

    algorithm_number <- unname(
      algorithm_numbers[[algorithm_name]]
    )

    if (verbose) {
      multi_clustering_message(
        "Starting algorithm `",
        algorithm_name,
        "` (Seurat algorithm = ",
        algorithm_number,
        ")."
      )
    }

    for (i in seq_along(resolutions)) {

      current_run <- current_run + 1L

      current_resolution <- resolutions[[i]]
      current_resolution_label <- resolution_labels[[i]]

      current_column <- paste0(
        algorithm_name,
        "_res",
        current_resolution_label
      )

      run_started <- Sys.time()

      if (verbose) {
        multi_clustering_message(
          sprintf(
            "[%02d/%02d] ",
            current_run,
            total_runs
          ),
          "START | algorithm=",
          algorithm_name,
          " | resolution=",
          current_resolution_label,
          " | column=",
          current_column
        )
      }

      findclusters_arguments <- list(
        object = graph_object,
        modularity.fxn = as.integer(modularity_fxn),
        resolution = current_resolution,
        algorithm = as.integer(algorithm_number),
        n.start = as.integer(n_start),
        n.iter = as.integer(n_iter),
        random.seed = as.integer(random_seed),
        group.singletons = group_singletons,
        verbose = verbose_findclusters
      )

      if (algorithm_number == 4L) {

        findclusters_arguments$leiden_method <-
          leiden_method

        findclusters_arguments$leiden_objective_function <-
          leiden_objective_function
      }

      clustering_result <- do.call(
        what = Seurat::FindClusters,
        args = findclusters_arguments
      )

      cluster_vector <- extract_single_clustering_result(
        clustering_result = clustering_result,
        expected_spots = spot_names
      )

      cluster_assignments[[current_column]] <-
        cluster_vector

      run_finished <- Sys.time()

      run_elapsed_seconds <- as.numeric(
        difftime(
          run_finished,
          run_started,
          units = "secs"
        )
      )

      total_elapsed_seconds <- as.numeric(
        difftime(
          run_finished,
          workflow_started,
          units = "secs"
        )
      )

      n_clusters <- nlevels(
        cluster_vector
      )

      summary_rows[[current_run]] <- data.frame(
        run = current_run,
        totalRuns = total_runs,
        algorithm = algorithm_name,
        algorithmNumber = algorithm_number,
        graphName = graph_name,
        resolution = current_resolution,
        clusterColumn = current_column,
        nSpots = length(cluster_vector),
        nClusters = n_clusters,
        modularityFunction = as.integer(modularity_fxn),
        nStart = as.integer(n_start),
        nIter = as.integer(n_iter),
        randomSeed = as.integer(random_seed),
        groupSingletons = group_singletons,
        leidenMethod = if (
          algorithm_number == 4L
        ) {
          leiden_method
        } else {
          NA_character_
        },
        leidenObjectiveFunction = if (
          algorithm_number == 4L
        ) {
          leiden_objective_function
        } else {
          NA_character_
        },
        elapsedSeconds = run_elapsed_seconds,
        elapsedFormatted =
          format_clustering_elapsed_time(
            run_elapsed_seconds
          ),
        completedAt = format(
          run_finished,
          "%Y-%m-%d %H:%M:%S %z"
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      if (verbose) {
        multi_clustering_message(
          sprintf(
            "[%02d/%02d] ",
            current_run,
            total_runs
          ),
          "DONE  | algorithm=",
          algorithm_name,
          " | resolution=",
          current_resolution_label,
          " | clusters=",
          n_clusters,
          " | runTime=",
          format_clustering_elapsed_time(
            run_elapsed_seconds
          ),
          " | totalTime=",
          format_clustering_elapsed_time(
            total_elapsed_seconds
          )
        )
      }

      rm(
        clustering_result,
        cluster_vector
      )

      invisible(
        gc(verbose = FALSE)
      )
    }
  }


  # ============================================================================
  # 4. Validate the complete assignment table
  #
  # No metadata has been added to the Seurat object yet.
  # ============================================================================

  if (
    nrow(cluster_assignments) != ncol(seurat_object) ||
      ncol(cluster_assignments) != length(cluster_columns)
  ) {
    stop(
      "Final cluster-assignment table has unexpected dimensions.",
      call. = FALSE
    )
  }

  if (!identical(
    rownames(cluster_assignments),
    spot_names
  )) {
    stop(
      "Cluster-assignment rows are not identically ordered as Seurat spots.",
      call. = FALSE
    )
  }

  if (!identical(
    colnames(cluster_assignments),
    cluster_columns
  )) {
    stop(
      "Cluster-assignment columns differ from the requested columns.",
      call. = FALSE
    )
  }

  if (anyNA(cluster_assignments)) {
    stop(
      "NA values were found in the completed cluster-assignment table.",
      call. = FALSE
    )
  }

  clustering_summary <- do.call(
    rbind,
    summary_rows
  )

  rownames(clustering_summary) <- NULL


  # ============================================================================
  # 5. Add all clustering columns to the Seurat object once
  #
  # This is the first and only metadata modification of the Seurat object.
  # ============================================================================

  if (verbose) {
    multi_clustering_message(
      "All ",
      total_runs,
      " clustering runs completed and validated. Adding ",
      length(cluster_columns),
      " metadata columns to the Seurat object in one operation."
    )
  }

  original_idents <- SeuratObject::Idents(
    seurat_object
  )

  updated_seurat_object <- SeuratObject::AddMetaData(
    object = seurat_object,
    metadata = cluster_assignments
  )

  if (is.null(set_active_column)) {

    SeuratObject::Idents(
      updated_seurat_object
    ) <- original_idents

  } else {

    if (
      !is.character(set_active_column) ||
        length(set_active_column) != 1L ||
        !set_active_column %in% cluster_columns
    ) {
      stop(
        "`set_active_column` must be NULL or one of the newly created ",
        "cluster columns.",
        call. = FALSE
      )
    }

    SeuratObject::Idents(
      updated_seurat_object
    ) <- set_active_column
  }


  # ============================================================================
  # 6. Store a compact record in object@misc
  # ============================================================================

  workflow_finished <- Sys.time()

  workflow_elapsed_seconds <- as.numeric(
    difftime(
      workflow_finished,
      workflow_started,
      units = "secs"
    )
  )

  run_id <- paste0(
    graph_name,
    "__multiClustering__",
    format(
      workflow_finished,
      "%Y%m%d_%H%M%S"
    )
  )

  if (
    is.null(
      updated_seurat_object@misc$multiClustering
    )
  ) {
    updated_seurat_object@misc$multiClustering <-
      list()
  }

  updated_seurat_object@misc$multiClustering[[run_id]] <-
    list(
      graphName = graph_name,
      algorithms = algorithms,
      algorithmNumbers =
        unname(algorithm_numbers[algorithms]),
      resolutions = resolutions,
      clusterColumns = cluster_columns,
      modularityFunction = as.integer(modularity_fxn),
      nStart = as.integer(n_start),
      nIter = as.integer(n_iter),
      randomSeed = as.integer(random_seed),
      groupSingletons = group_singletons,
      leidenMethod = if (
        "leiden" %in% algorithms
      ) {
        leiden_method
      } else {
        NA_character_
      },
      leidenObjectiveFunction = if (
        "leiden" %in% algorithms
      ) {
        leiden_objective_function
      } else {
        NA_character_
      },
      workflowElapsedSeconds =
        workflow_elapsed_seconds,
      completedAt = format(
        workflow_finished,
        "%Y-%m-%d %H:%M:%S %z"
      ),
      clusteringSummary = clustering_summary
    )

  cluster_assignments_table <- data.frame(
    spot = rownames(cluster_assignments),
    cluster_assignments,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    multi_clustering_message(
      "Finished multi-clustering. Total elapsed time: ",
      format_clustering_elapsed_time(
        workflow_elapsed_seconds
      ),
      "."
    )
  }

  list(
    seurat_object = updated_seurat_object,
    cluster_assignments =
      cluster_assignments_table,
    clustering_summary =
      clustering_summary,
    run_id = run_id
  )
}
