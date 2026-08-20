# ==============================================================================
# functions_multiClustering.R
#
# Purpose:
# Reusable Seurat v5 functions for running multiple graph-based clustering
# algorithms over multiple resolution values on an existing SNN graph.
#
# Supported Seurat FindClusters algorithms:
#   1. Louvain
#   2. Louvain with multilevel refinement
#   3. SLM
#   4. Leiden
#
# This file does not perform integration, FindNeighbors(), UMAP, plotting,
# or table export.
# ==============================================================================


# ==============================================================================
# 1. General validation and I/O
# ==============================================================================

check_multi_clustering_packages <- function(packages) {
  if (!is.character(packages) || length(packages) == 0L) {
    stop("`packages` must be a non-empty character vector.", call. = FALSE)
  }

  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0L) {
    stop(
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


assert_seurat_v5_for_clustering <- function() {
  check_multi_clustering_packages(c("Seurat", "SeuratObject"))

  seurat_version <- utils::packageVersion("Seurat")
  seurat_object_version <- utils::packageVersion("SeuratObject")

  message("Seurat version: ", as.character(seurat_version))
  message("SeuratObject version: ", as.character(seurat_object_version))

  if (seurat_version < "5.0.0") {
    stop(
      "This workflow requires Seurat v5 or newer. Current version: ",
      as.character(seurat_version),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


assert_seurat_object_for_clustering <- function(object) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }

  invisible(TRUE)
}


load_single_seurat_for_clustering <- function(rdata_file) {
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
      "Expected exactly one Seurat object in `", rdata_file,
      "`, but found ", length(seurat_names), ".",
      call. = FALSE
    )
  }

  message("Loaded Seurat object: ", seurat_names)
  get(seurat_names, envir = input_environment)
}


save_single_seurat_after_clustering <- function(
    object,
    object_name,
    output_file,
    compress = TRUE
) {
  assert_seurat_object_for_clustering(object)

  if (
    length(object_name) != 1L || is.na(object_name) ||
      !nzchar(object_name) || make.names(object_name) != object_name
  ) {
    stop(
      "`object_name` must be one valid, non-empty R object name.",
      call. = FALSE
    )
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  output_environment <- new.env(parent = emptyenv())
  assign(object_name, object, envir = output_environment)

  save(
    list = object_name,
    file = output_file,
    envir = output_environment,
    compress = compress
  )

  invisible(output_file)
}


# ==============================================================================
# 2. Formatting helpers
# ==============================================================================

emit_clustering_progress <- function(..., log_connection = NULL) {
  progress_text <- paste0(..., collapse = "")
  progress_line <- paste0(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
    progress_text
  )

  cat(progress_line, "\n", sep = "")
  flush.console()

  if (
    !is.null(log_connection) &&
      inherits(log_connection, "connection") &&
      isOpen(log_connection)
  ) {
    flush(log_connection)
  }

  invisible(progress_line)
}


format_elapsed_seconds <- function(seconds) {
  if (length(seconds) != 1L || is.na(seconds) || !is.finite(seconds)) {
    return(NA_character_)
  }

  seconds <- max(0, as.numeric(seconds))

  hours <- floor(seconds / 3600)
  minutes <- floor((seconds %% 3600) / 60)
  remaining_seconds <- seconds %% 60

  sprintf(
    "%02d:%02d:%05.2f",
    as.integer(hours),
    as.integer(minutes),
    remaining_seconds
  )
}


validate_resolutions <- function(resolutions) {
  if (
    !is.numeric(resolutions) || length(resolutions) == 0L ||
      anyNA(resolutions) || any(!is.finite(resolutions)) ||
      any(resolutions <= 0)
  ) {
    stop(
      "`resolutions` must be a non-empty numeric vector containing values > 0.",
      call. = FALSE
    )
  }

  rounded_resolutions <- round(as.numeric(resolutions), digits = 10)

  if (anyDuplicated(rounded_resolutions) > 0L) {
    stop("`resolutions` contains duplicated values.", call. = FALSE)
  }

  invisible(TRUE)
}


format_resolution_value <- function(resolution, digits = 2L) {
  if (
    length(resolution) != 1L || is.na(resolution) ||
      !is.finite(resolution) || resolution <= 0
  ) {
    stop("`resolution` must be one finite number greater than 0.", call. = FALSE)
  }

  formatC(
    as.numeric(resolution),
    format = "f",
    digits = as.integer(digits)
  )
}


normalize_clustering_algorithms <- function(algorithms) {
  supported_algorithms <- c(
    louvain = 1L,
    louvainRefined = 2L,
    slm = 3L,
    leiden = 4L
  )

  if (is.character(algorithms)) {
    normalized_names <- gsub("[^a-zA-Z0-9]", "", algorithms)
    normalized_names <- tolower(normalized_names)

    aliases <- c(
      louvain = "louvain",
      originalLouvain = "louvain",
      louvainrefined = "louvainRefined",
      refinedlouvain = "louvainRefined",
      louvainmultilevelrefinement = "louvainRefined",
      slm = "slm",
      smartlocalmoving = "slm",
      leiden = "leiden"
    )
    names(aliases) <- tolower(names(aliases))

    unsupported <- setdiff(normalized_names, names(aliases))

    if (length(unsupported) > 0L) {
      stop(
        "Unsupported clustering algorithm name(s): ",
        paste(algorithms[normalized_names %in% unsupported], collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    canonical_names <- unname(aliases[normalized_names])
    result <- supported_algorithms[canonical_names]
    names(result) <- canonical_names
  } else if (is.numeric(algorithms)) {
    algorithms <- as.integer(algorithms)

    if (anyNA(algorithms) || any(!algorithms %in% supported_algorithms)) {
      stop(
        "Numeric `algorithms` must contain only 1, 2, 3, or 4.",
        call. = FALSE
      )
    }

    canonical_names <- names(supported_algorithms)[
      match(algorithms, supported_algorithms)
    ]
    result <- algorithms
    names(result) <- canonical_names
  } else {
    stop(
      "`algorithms` must be a character vector or a numeric vector.",
      call. = FALSE
    )
  }

  if (anyDuplicated(names(result)) > 0L) {
    stop("`algorithms` contains duplicated methods.", call. = FALSE)
  }

  result
}


make_clustering_column_name <- function(
    algorithm_label,
    resolution,
    resolution_digits = 2L
) {
  if (
    length(algorithm_label) != 1L || is.na(algorithm_label) ||
      !nzchar(algorithm_label)
  ) {
    stop("`algorithm_label` must be one non-empty value.", call. = FALSE)
  }

  paste0(
    algorithm_label,
    "_res",
    format_resolution_value(resolution, digits = resolution_digits)
  )
}


# ==============================================================================
# 3. Graph and dependency validation
# ==============================================================================

validate_snn_graph_for_clustering <- function(object, graph_name) {
  assert_seurat_object_for_clustering(object)

  if (
    length(graph_name) != 1L || is.na(graph_name) || !nzchar(graph_name)
  ) {
    stop("`graph_name` must be one non-empty value.", call. = FALSE)
  }

  if (!graph_name %in% names(object@graphs)) {
    stop(
      "The Seurat object does not contain graph `", graph_name,
      "`. Available graphs: ",
      paste(names(object@graphs), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  graph_object <- object[[graph_name]]

  if (nrow(graph_object) != ncol(object)) {
    stop(
      "Graph `", graph_name, "` contains ", nrow(graph_object),
      " vertices, but the object contains ", ncol(object), " spots.",
      call. = FALSE
    )
  }

  message(
    "Validated graph `", graph_name, "`: ",
    nrow(graph_object), " vertices."
  )

  invisible(TRUE)
}


validate_leiden_dependency <- function(leiden_method) {
  leiden_method <- match.arg(leiden_method, c("leidenbase", "igraph"))

  if (leiden_method == "leidenbase") {
    check_multi_clustering_packages("leidenbase")
  }

  if (leiden_method == "igraph") {
    check_multi_clustering_packages("igraph")
  }

  invisible(leiden_method)
}


# ==============================================================================
# 4. Universal multi-algorithm, multi-resolution clustering
# ==============================================================================

run_seurat_multi_clustering <- function(
    object,
    graph_name,
    resolutions = round(seq(0.05, 2.00, by = 0.05), 2),
    algorithms = c("louvain", "louvainRefined", "slm", "leiden"),
    modularity_fxn = 1L,
    n_start = 10L,
    n_iter = 10L,
    random_seed = 7L,
    group_singletons = TRUE,
    leiden_method = "leidenbase",
    leiden_objective_function = "modularity",
    resolution_digits = 2L,
    overwrite_existing = FALSE,
    restore_original_idents = TRUE,
    remove_seurat_clusters_column = TRUE,
    verbose = FALSE,
    log_connection = NULL
) {
  assert_seurat_v5_for_clustering()
  validate_snn_graph_for_clustering(object, graph_name)
  validate_resolutions(resolutions)

  algorithm_map <- normalize_clustering_algorithms(algorithms)
  resolutions <- round(as.numeric(resolutions), digits = 10)

  if (
    length(modularity_fxn) != 1L || is.na(modularity_fxn) ||
      !modularity_fxn %in% c(1L, 2L)
  ) {
    stop("`modularity_fxn` must be 1 or 2.", call. = FALSE)
  }

  for (parameter in c("n_start", "n_iter", "random_seed")) {
    value <- get(parameter)

    if (
      length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 1 || value != as.integer(value)
    ) {
      stop("`", parameter, "` must be one positive integer.", call. = FALSE)
    }
  }

  if (!is.logical(group_singletons) || length(group_singletons) != 1L) {
    stop("`group_singletons` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(overwrite_existing) || length(overwrite_existing) != 1L) {
    stop("`overwrite_existing` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(restore_original_idents) || length(restore_original_idents) != 1L) {
    stop("`restore_original_idents` must be TRUE or FALSE.", call. = FALSE)
  }

  if (
    !is.logical(remove_seurat_clusters_column) ||
      length(remove_seurat_clusters_column) != 1L
  ) {
    stop(
      "`remove_seurat_clusters_column` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  leiden_method <- match.arg(leiden_method, c("leidenbase", "igraph"))
  leiden_objective_function <- match.arg(
    leiden_objective_function,
    c("modularity", "CPM")
  )

  if (4L %in% unname(algorithm_map)) {
    validate_leiden_dependency(leiden_method)
  }

  original_idents <- SeuratObject::Idents(object)
  had_seurat_clusters <- "seurat_clusters" %in% colnames(object[[]])

  if (had_seurat_clusters) {
    original_seurat_clusters <- object[["seurat_clusters"]][, 1]
    names(original_seurat_clusters) <- colnames(object)
  } else {
    original_seurat_clusters <- NULL
  }

  expected_columns <- unlist(
    lapply(
      names(algorithm_map),
      function(algorithm_label) {
        vapply(
          resolutions,
          function(resolution) {
            make_clustering_column_name(
              algorithm_label = algorithm_label,
              resolution = resolution,
              resolution_digits = resolution_digits
            )
          },
          character(1)
        )
      }
    ),
    use.names = FALSE
  )

  existing_columns <- intersect(expected_columns, colnames(object[[]]))

  if (length(existing_columns) > 0L && !isTRUE(overwrite_existing)) {
    stop(
      "The object already contains requested clustering column(s): ",
      paste(existing_columns, collapse = ", "),
      ". Set `overwrite_existing = TRUE` to replace them.",
      call. = FALSE
    )
  }

  run_started <- Sys.time()
  total_runs <- length(algorithm_map) * length(resolutions)
  summary_rows <- vector(
    mode = "list",
    length = total_runs
  )
  summary_index <- 0L
  run_index <- 0L

  emit_clustering_progress(
    "Starting multi-clustering on graph `", graph_name, "` with ",
    length(algorithm_map), " algorithms and ",
    length(resolutions), " resolution values (",
    total_runs, " clustering runs in total).",
    log_connection = log_connection
  )

  for (algorithm_label in names(algorithm_map)) {
    algorithm_number <- unname(algorithm_map[[algorithm_label]])
    algorithm_started <- Sys.time()

    emit_clustering_progress(
      "Starting algorithm `", algorithm_label,
      "` (Seurat algorithm = ", algorithm_number, ").",
      log_connection = log_connection
    )

    for (resolution in resolutions) {
      clustering_column <- make_clustering_column_name(
        algorithm_label = algorithm_label,
        resolution = resolution,
        resolution_digits = resolution_digits
      )

      clustering_started <- Sys.time()
      run_index <- run_index + 1L

      emit_clustering_progress(
        sprintf("[%03d/%03d] START | ", run_index, total_runs),
        "algorithm=", algorithm_label,
        " | resolution=", format_resolution_value(resolution, resolution_digits),
        " | column=", clustering_column,
        log_connection = log_connection
      )

      clustering_arguments <- list(
        object = object,
        graph.name = graph_name,
        cluster.name = clustering_column,
        modularity.fxn = as.integer(modularity_fxn),
        resolution = resolution,
        algorithm = as.integer(algorithm_number),
        n.start = as.integer(n_start),
        n.iter = as.integer(n_iter),
        random.seed = as.integer(random_seed),
        group.singletons = group_singletons,
        verbose = verbose
      )

      if (algorithm_number == 4L) {
        clustering_arguments$leiden_method <- leiden_method
        clustering_arguments$leiden_objective_function <-
          leiden_objective_function
      }

      object <- do.call(
        what = Seurat::FindClusters,
        args = clustering_arguments
      )

      if (!clustering_column %in% colnames(object[[]])) {
        stop(
          "FindClusters() did not create expected metadata column `",
          clustering_column,
          "`.",
          call. = FALSE
        )
      }

      cluster_values <- object[[clustering_column]][, 1]
      n_clusters <- length(unique(as.character(cluster_values)))

      clustering_finished <- Sys.time()
      elapsed_seconds <- as.numeric(
        difftime(
          clustering_finished,
          clustering_started,
          units = "secs"
        )
      )

      summary_index <- summary_index + 1L
      summary_rows[[summary_index]] <- data.frame(
        algorithm = algorithm_label,
        algorithmNumber = as.integer(algorithm_number),
        resolution = as.numeric(resolution),
        metadataColumn = clustering_column,
        nClusters = as.integer(n_clusters),
        elapsedSeconds = elapsed_seconds,
        elapsedTime = format_elapsed_seconds(elapsed_seconds),
        stringsAsFactors = FALSE
      )

      elapsed_total_seconds <- as.numeric(
        difftime(clustering_finished, run_started, units = "secs")
      )
      mean_seconds_per_run <- elapsed_total_seconds / run_index
      estimated_remaining_seconds <- mean_seconds_per_run * (total_runs - run_index)
      estimated_finish_time <- Sys.time() + estimated_remaining_seconds

      emit_clustering_progress(
        sprintf("[%03d/%03d] DONE  | ", run_index, total_runs),
        "algorithm=", algorithm_label,
        " | resolution=", format_resolution_value(resolution, resolution_digits),
        " | clusters=", n_clusters,
        " | runTime=", format_elapsed_seconds(elapsed_seconds),
        " | totalTime=", format_elapsed_seconds(elapsed_total_seconds),
        " | remaining~", format_elapsed_seconds(estimated_remaining_seconds),
        " | ETA=", format(estimated_finish_time, "%Y-%m-%d %H:%M:%S"),
        log_connection = log_connection
      )
    }

    algorithm_elapsed <- as.numeric(
      difftime(Sys.time(), algorithm_started, units = "secs")
    )

    emit_clustering_progress(
      "Completed algorithm `", algorithm_label,
      "`; elapsed ", format_elapsed_seconds(algorithm_elapsed), ".",
      log_connection = log_connection
    )
  }

  clustering_summary <- do.call(rbind, summary_rows)
  rownames(clustering_summary) <- NULL

  if (isTRUE(remove_seurat_clusters_column)) {
    if (had_seurat_clusters) {
      object[["seurat_clusters"]] <- original_seurat_clusters
    } else if ("seurat_clusters" %in% colnames(object[[]])) {
      object[["seurat_clusters"]] <- NULL
    }
  }

  if (isTRUE(restore_original_idents)) {
    SeuratObject::Idents(object) <- original_idents
  }

  total_elapsed_seconds <- as.numeric(
    difftime(Sys.time(), run_started, units = "secs")
  )

  clustering_info <- list(
    graphName = graph_name,
    algorithms = algorithm_map,
    resolutions = resolutions,
    modularityFunction = as.integer(modularity_fxn),
    nStart = as.integer(n_start),
    nIter = as.integer(n_iter),
    randomSeed = as.integer(random_seed),
    groupSingletons = group_singletons,
    leidenMethod = leiden_method,
    leidenObjectiveFunction = leiden_objective_function,
    resolutionDigits = as.integer(resolution_digits),
    clusteringColumns = expected_columns,
    nClusteringRuns = length(expected_columns),
    totalElapsedSeconds = total_elapsed_seconds,
    totalElapsedTime = format_elapsed_seconds(total_elapsed_seconds),
    summary = clustering_summary,
    completedAt = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  )

  if (is.null(object@misc$seuratClusteringAnalysis)) {
    object@misc$seuratClusteringAnalysis <- list()
  }

  object@misc$seuratClusteringAnalysis$multiClustering <- clustering_info

  emit_clustering_progress(
    "Multi-clustering completed: ", length(expected_columns),
    " metadata columns; total elapsed ",
    format_elapsed_seconds(total_elapsed_seconds), ".",
    log_connection = log_connection
  )

  list(
    object = object,
    summary = clustering_summary,
    clustering_info = clustering_info
  )
}
