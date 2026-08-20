# ==============================================================================
# 02_run_rpcaDims20K20PruneSNN0067_multiClustering_externalTable.R
#
# Purpose:
# Run Louvain, Louvain refined, SLM and Leiden clustering directly on the
# existing RPCA SNN graph for resolutions 0.05-2.00.
#
# Output:
# - external TSV.GZ table with spot-level cluster assignments;
# - external TSV table with run times and numbers of clusters;
# - log file.
#
# The input Seurat object and all RData files remain unchanged.
# ==============================================================================


main <- function() {

  log_message <- function(...) {
    cat(
      paste0(..., collapse = ""),
      "\n",
      sep = ""
    )
    flush.console()
    invisible(NULL)
  }

  # =============================================================================
  # 1. Configuration
  # =============================================================================

  analysis_started <- Sys.time()

  dataset_name <- "maternalFMT_n16samples"
  normalization_label <- "logNormalizeVst"
  nfeatures <- 2000L

  graph_configuration_tag <- "rpcaDims20K20PruneSNN0067"
  snn_graph_name <- paste0(
    graph_configuration_tag,
    "_snn"
  )

  clustering_algorithms <- c(
    "louvain",
    "louvainRefined",
    "slm",
    "leiden"
  )

  resolutions <- round(
    seq(
      from = 0.05,
      to = 2.00,
      by = 0.05
    ),
    digits = 2
  )

  modularity_fxn <- 1L
  n_start <- 10L
  n_iter <- 10L
  seed <- 7L
  group_singletons <- TRUE

  leiden_method <- "leidenbase"
  leiden_objective_function <- "modularity"


  # =============================================================================
  # 2. Paths
  # =============================================================================

  project_root <- normalizePath(
    getwd(),
    mustWork = TRUE
  )

  functions_file <- file.path(
    project_root,
    "preprocessing",
    "src",
    "seurat_clustering_analysis",
    "logNormalize_vst",
    "functions_multiClustering_externalTable.R"
  )

  analysis_root <- file.path(
    project_root,
    "results",
    dataset_name,
    "seurat_clustering_analysis",
    "logNormalize_vst",
    paste0("hvg", nfeatures),
    "02_rpca",
    "dims20",
    "k20",
    "prune0067"
  )

  input_rdata <- file.path(
    analysis_root,
    "RData",
    paste0(
      "01_",
      dataset_name,
      "_",
      normalization_label,
      "_hvg",
      nfeatures,
      "_",
      graph_configuration_tag,
      "_neighborsGraph.RData"
    )
  )

  tables_directory <- file.path(
    analysis_root,
    "tables"
  )

  logs_directory <- file.path(
    analysis_root,
    "logs"
  )

  dir.create(
    tables_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    logs_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  output_prefix <- paste0(
    "02_",
    dataset_name,
    "_",
    normalization_label,
    "_hvg",
    nfeatures,
    "_",
    graph_configuration_tag,
    "_multiClustering"
  )

  cluster_assignments_file <- file.path(
    tables_directory,
    paste0(
      output_prefix,
      "_clusterAssignments.tsv.gz"
    )
  )

  clustering_summary_file <- file.path(
    tables_directory,
    paste0(
      output_prefix,
      "_clusteringSummary.tsv"
    )
  )

  log_file <- file.path(
    logs_directory,
    paste0(
      output_prefix,
      "_externalTableRun.log"
    )
  )

  if (!file.exists(functions_file)) {
    stop(
      "Functions file does not exist: ",
      functions_file,
      call. = FALSE
    )
  }

  if (!file.exists(input_rdata)) {
    stop(
      "Input RData file does not exist: ",
      input_rdata,
      call. = FALSE
    )
  }


  # =============================================================================
  # 3. Logging
  # =============================================================================

  log_connection <- file(
    log_file,
    open = "wt"
  )

  sink(
    log_connection,
    type = "output",
    split = TRUE
  )

  on.exit(
    {
      while (sink.number(type = "output") > 0L) {
        sink(type = "output")
      }

      close(log_connection)
    },
    add = TRUE
  )

  log_message(
    "Analysis started: ",
    format(
      analysis_started,
      "%Y-%m-%d %H:%M:%S %z"
    )
  )

  log_message("Project root: ", project_root)
  log_message("Input RData: ", input_rdata)
  log_message("SNN graph: ", snn_graph_name)
  log_message(
    "Cluster assignments table: ",
    cluster_assignments_file
  )
  log_message(
    "Clustering summary table: ",
    clustering_summary_file
  )
  log_message("Log file: ", log_file)
  log_message(
    "Important: no RData file will be created or overwritten."
  )


  # =============================================================================
  # 4. Packages and functions
  # =============================================================================

  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
  })

  source(functions_file)

  if (
    utils::packageVersion("Seurat") < "5.0.0"
  ) {
    stop(
      "Seurat v5 or newer is required.",
      call. = FALSE
    )
  }

  if (
    "leiden" %in% clustering_algorithms &&
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

  set.seed(seed)


  # =============================================================================
  # 5. Load the Seurat object and extract only the SNN graph
  # =============================================================================

  loaded <- load_single_seurat_object(
    input_rdata
  )

  seurat_object <- loaded$object

  log_message(
    "Loaded Seurat object: ",
    loaded$object_name
  )

  log_message(
    "Input object: ",
    nrow(seurat_object),
    " genes x ",
    ncol(seurat_object),
    " spots."
  )

  log_message(
    "Stored graphs: ",
    paste(
      names(seurat_object@graphs),
      collapse = ", "
    )
  )

  if (
    !snn_graph_name %in%
      names(seurat_object@graphs)
  ) {
    stop(
      "SNN graph `",
      snn_graph_name,
      "` is absent from the Seurat object.",
      call. = FALSE
    )
  }

  graph_object <- seurat_object[[snn_graph_name]]
  spot_names <- colnames(seurat_object)

  log_message(
    "Extracted graph `",
    snn_graph_name,
    "` with ",
    nrow(graph_object),
    " vertices."
  )

  # The complete Seurat object is no longer needed for clustering.
  rm(seurat_object, loaded)
  invisible(gc(verbose = FALSE))


  # =============================================================================
  # 6. Run clustering directly on the graph
  # =============================================================================

  clustering_results <-
    run_multi_clustering_to_external_table(
      graph_object = graph_object,
      spot_names = spot_names,
      algorithms = clustering_algorithms,
      resolutions = resolutions,
      modularity_fxn = modularity_fxn,
      n_start = n_start,
      n_iter = n_iter,
      random_seed = seed,
      group_singletons = group_singletons,
      leiden_method = leiden_method,
      leiden_objective_function =
        leiden_objective_function,
      verbose_findclusters = TRUE
    )


  # =============================================================================
  # 7. Save only external tables
  # =============================================================================

  cluster_assignments_connection <- gzfile(
    cluster_assignments_file,
    open = "wt"
  )

  write.table(
    clustering_results$cluster_assignments,
    file = cluster_assignments_connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )

  close(cluster_assignments_connection)

  write.table(
    clustering_results$clustering_summary,
    file = clustering_summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )

  analysis_finished <- Sys.time()

  total_elapsed_seconds <- as.numeric(
    difftime(
      analysis_finished,
      analysis_started,
      units = "secs"
    )
  )

  log_message(
    "Saved cluster assignments: ",
    cluster_assignments_file
  )

  log_message(
    "Saved clustering summary: ",
    clustering_summary_file
  )

  log_message(
    "Number of clustering columns: ",
    ncol(
      clustering_results$cluster_assignments
    ) - 1L
  )

  log_message(
    "Total elapsed time: ",
    format_elapsed_time(
      total_elapsed_seconds
    )
  )

  log_message(
    "Analysis completed: ",
    format(
      analysis_finished,
      "%Y-%m-%d %H:%M:%S %z"
    )
  )

  log_message(
    "No RData file was created or overwritten."
  )

  invisible(TRUE)
}


tryCatch(
  main(),
  error = function(error_condition) {

    cat(
      "ERROR: ",
      conditionMessage(error_condition),
      "\n",
      sep = ""
    )
    flush.console()

    quit(
      save = "no",
      status = 1L,
      runLast = FALSE
    )
  }
)
