# ==============================================================================
# functions_multiClustering_externalTable.R
#
# Purpose:
# Run Seurat clustering directly on an existing SNN graph and store cluster
# assignments in an external table.
#
# Important:
# - The input Seurat object is never modified.
# - No Seurat RData file is written or overwritten.
# - Clustering is sequential: one algorithm and one resolution at a time.
# - No future parallelization is used.
# ==============================================================================


load_single_seurat_object <- function(rdata_file) {

  if (!file.exists(rdata_file)) {
    stop("Input RData file does not exist: ", rdata_file, call. = FALSE)
  }

  input_environment <- new.env(parent = emptyenv())
  loaded_names <- load(rdata_file, envir = input_environment)

  seurat_names <- loaded_names[
    vapply(
      loaded_names,
      function(object_name) {
        inherits(get(object_name, envir = input_environment), "Seurat")
      },
      logical(1)
    )
  ]

  if (length(seurat_names) != 1L) {
    stop(
      "Expected exactly one Seurat object in `",
      rdata_file,
      "`, but found ",
      length(seurat_names),
      ".",
      call. = FALSE
    )
  }

  list(
    object = get(seurat_names, envir = input_environment),
    object_name = seurat_names[[1]]
  )
}


format_resolution <- function(resolution, digits = 2L) {

  formatC(
    as.numeric(resolution),
    format = "f",
    digits = as.integer(digits)
  )
}


format_elapsed_time <- function(seconds) {

  seconds <- max(0, as.numeric(seconds))

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


progress_message <- function(...) {

  cat(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste0(..., collapse = ""),
    "\n",
    sep = ""
  )

  flush.console()
  invisible(NULL)
}


extract_single_clustering <- function(
    clustering_result,
    expected_spots
) {

  if (is.data.frame(clustering_result) || is.matrix(clustering_result)) {

    if (ncol(clustering_result) != 1L) {
      stop(
        "FindClusters returned ",
        ncol(clustering_result),
        " columns; exactly one was expected.",
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
      paste(class(clustering_result), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (is.null(names(clusters))) {

    if (length(clusters) != length(expected_spots)) {
      stop(
        "Clustering result has no spot names and its length does not match ",
        "the number of expected spots.",
        call. = FALSE
      )
    }

    names(clusters) <- expected_spots
  }

  missing_spots <- setdiff(expected_spots, names(clusters))
  extra_spots <- setdiff(names(clusters), expected_spots)

  if (length(missing_spots) > 0L || length(extra_spots) > 0L) {
    stop(
      "Spot names returned by FindClusters do not match the input graph.",
      call. = FALSE
    )
  }

  clusters <- clusters[expected_spots]

  cluster_character <- as.character(clusters)

  if (all(grepl("^-?[0-9]+$", cluster_character))) {
    return(as.integer(cluster_character))
  }

  cluster_character
}


run_multi_clustering_to_external_table <- function(
    graph_object,
    spot_names,
    algorithms,
    resolutions,
    modularity_fxn = 1L,
    n_start = 10L,
    n_iter = 10L,
    random_seed = 7L,
    group_singletons = TRUE,
    leiden_method = "leidenbase",
    leiden_objective_function = "modularity",
    verbose_findclusters = TRUE
) {

  if (!inherits(graph_object, "Graph")) {
    stop("`graph_object` must be a Seurat Graph object.", call. = FALSE)
  }

  if (nrow(graph_object) != length(spot_names)) {
    stop(
      "Graph vertex count does not match the number of spot names.",
      call. = FALSE
    )
  }

  algorithm_numbers <- c(
    louvain = 1L,
    louvainRefined = 2L,
    slm = 3L,
    leiden = 4L
  )

  unknown_algorithms <- setdiff(algorithms, names(algorithm_numbers))

  if (length(unknown_algorithms) > 0L) {
    stop(
      "Unsupported clustering algorithm(s): ",
      paste(unknown_algorithms, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  total_runs <- length(algorithms) * length(resolutions)
  current_run <- 0L
  workflow_started <- Sys.time()

  cluster_assignments <- data.frame(
    spot = spot_names,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  summary_rows <- vector("list", total_runs)

  progress_message(
    "Starting clustering directly on the SNN graph: ",
    length(algorithms),
    " algorithms x ",
    length(resolutions),
    " resolutions = ",
    total_runs,
    " runs."
  )

  for (algorithm_label in algorithms) {

    algorithm_number <- unname(
      algorithm_numbers[[algorithm_label]]
    )

    progress_message(
      "Starting algorithm `",
      algorithm_label,
      "` (Seurat algorithm = ",
      algorithm_number,
      ")."
    )

    for (resolution in resolutions) {

      current_run <- current_run + 1L
      resolution_label <- format_resolution(resolution, digits = 2L)

      output_column <- paste0(
        algorithm_label,
        "_res",
        resolution_label
      )

      run_started <- Sys.time()

      progress_message(
        sprintf("[%03d/%03d] ", current_run, total_runs),
        "START | algorithm=",
        algorithm_label,
        " | resolution=",
        resolution_label,
        " | column=",
        output_column
      )

      findclusters_arguments <- list(
        object = graph_object,
        modularity.fxn = as.integer(modularity_fxn),
        resolution = as.numeric(resolution),
        algorithm = as.integer(algorithm_number),
        n.start = as.integer(n_start),
        n.iter = as.integer(n_iter),
        random.seed = as.integer(random_seed),
        group.singletons = group_singletons,
        verbose = verbose_findclusters
      )

      if (algorithm_number == 4L) {
        findclusters_arguments$leiden_method <- leiden_method
        findclusters_arguments$leiden_objective_function <-
          leiden_objective_function
      }

      clustering_result <- do.call(
        Seurat::FindClusters,
        findclusters_arguments
      )

      cluster_vector <- extract_single_clustering(
        clustering_result = clustering_result,
        expected_spots = spot_names
      )

      cluster_assignments[[output_column]] <- cluster_vector

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

      n_clusters <- length(unique(cluster_vector))

      summary_rows[[current_run]] <- data.frame(
        run = current_run,
        totalRuns = total_runs,
        algorithm = algorithm_label,
        algorithmNumber = algorithm_number,
        resolution = as.numeric(resolution),
        metadataColumn = output_column,
        nClusters = n_clusters,
        elapsedSeconds = run_elapsed_seconds,
        elapsedFormatted = format_elapsed_time(run_elapsed_seconds),
        completedAt = format(
          run_finished,
          "%Y-%m-%d %H:%M:%S %z"
        ),
        stringsAsFactors = FALSE
      )

      remaining_runs <- total_runs - current_run

      mean_seconds_per_run <- total_elapsed_seconds / current_run
      estimated_remaining_seconds <-
        mean_seconds_per_run * remaining_runs

      estimated_finish <- Sys.time() + estimated_remaining_seconds

      progress_message(
        sprintf("[%03d/%03d] ", current_run, total_runs),
        "DONE  | algorithm=",
        algorithm_label,
        " | resolution=",
        resolution_label,
        " | clusters=",
        n_clusters,
        " | runTime=",
        format_elapsed_time(run_elapsed_seconds),
        " | totalTime=",
        format_elapsed_time(total_elapsed_seconds),
        " | remaining~",
        format_elapsed_time(estimated_remaining_seconds),
        " | ETA=",
        format(estimated_finish, "%Y-%m-%d %H:%M:%S")
      )

      rm(clustering_result, cluster_vector)
      invisible(gc(verbose = FALSE))
    }
  }

  clustering_summary <- do.call(
    rbind,
    summary_rows
  )

  rownames(clustering_summary) <- NULL

  list(
    cluster_assignments = cluster_assignments,
    clustering_summary = clustering_summary
  )
}
