# ==============================================================================
# 01_prepare_rpcaDims20K20PruneSNN0067_neighborsGraph.R
#
# Purpose:
# Create one pre-clustering Seurat object for maternalFMT_n16samples:
#   - input: LogNormalize + VST, 2000 HVGs, PCA50 parent object,
#   - integration: RPCA,
#   - dimensions: 1:20,
#   - k.param: 20,
#   - prune.SNN: 1/15 (~0.0667),
#   - output: integrated reduction plus kNN/SNN graphs.
#
# This script performs no FindClusters() and no UMAP calculation.
# It writes one RData object and one log file.
# ============================================================================== 


main <- function() {

  # =============================================================================
  # 1. Analysis configuration
  # =============================================================================

  dataset_name <- "maternalFMT_n16samples"
  normalization_label <- "logNormalizeVst"
  nfeatures <- 2000L
  pca_components_in_parent <- 50L

  integration_method <- "rpca"
  dims <- 1:20
  k_param <- 20L
  prune_snn <- 1 / 15
  seed <- 7L

  expected_n_samples <- 16L
  sample_id_col <- "sample_ID"
  assay <- "RNA"
  original_reduction <- "pca"
  integrated_reduction <- "integrated.rpca"
  normalization_method <- "LogNormalize"

  # Keep method-specific integration parameters at Seurat defaults for now.
  integration_method_args <- list()


  # =============================================================================
  # 2. Project paths
  # =============================================================================

  project_root <- Sys.getenv(
    "IPPAS_SPATIAL_PROJECT_ROOT",
    unset = getwd()
  )
  project_root <- normalizePath(project_root, mustWork = TRUE)

  functions_file <- file.path(
    project_root,
    "preprocessing",
    "src",
    "seurat_clustering_analysis",
    "logNormalize_vst",
    "functions_integration_neighborsGraph.R"
  )

  input_rdata <- file.path(
    project_root,
    "results",
    dataset_name,
    "seurat_clustering_analysis",
    "logNormalize_vst",
    paste0("hvg", nfeatures),
    "00_preprocessing",
    "RData",
    paste0(
      "01_", dataset_name, "_", normalization_label,
      "_hvg", nfeatures,
      "_pca", pca_components_in_parent,
      "_parentObject.RData"
    )
  )

  if (!file.exists(functions_file)) {
    stop("Functions file does not exist: ", functions_file, call. = FALSE)
  }

  source(functions_file)

  configuration_tag <- make_integration_graph_tag(
    integration_method = integration_method,
    dims = dims,
    k_param = k_param,
    prune_snn = prune_snn
  )

  graph_names <- make_graph_names(configuration_tag)

  dims_folder <- format_dims_label(dims)
  k_folder <- paste0("k", k_param)
  prune_folder <- paste0("prune", format_prune_snn_code(prune_snn))

  output_root <- file.path(
    project_root,
    "results",
    dataset_name,
    "seurat_clustering_analysis",
    "logNormalize_vst",
    paste0("hvg", nfeatures),
    "02_rpca",
    dims_folder,
    k_folder,
    prune_folder
  )

  rdata_directory <- file.path(output_root, "RData")

  # Keep execution logs outside results/. The log tree mirrors the analysis
  # structure used for scripts and results.
  logs_directory <- file.path(
    project_root,
    "preprocessing",
    "logs",
    dataset_name,
    "seurat_clustering_analysis",
    "logNormalize_vst",
    paste0("hvg", nfeatures),
    "02_rpca",
    dims_folder,
    k_folder,
    prune_folder
  )

  dir.create(rdata_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_directory, recursive = TRUE, showWarnings = FALSE)

  file_prefix <- paste0(
    dataset_name,
    "_", normalization_label,
    "_hvg", nfeatures,
    "_", configuration_tag
  )

  output_rdata <- file.path(
    rdata_directory,
    paste0("01_", file_prefix, "_neighborsGraph.RData")
  )

  output_object_name <- paste0(file_prefix, "_neighborsGraph")

  log_file <- file.path(
    logs_directory,
    paste0(file_prefix, "_integrationNeighborsRun.log")
  )


  # =============================================================================
  # 3. Start logging
  # =============================================================================

  log_connection <- file(log_file, open = "wt")
  sink(log_connection, type = "output", split = TRUE)
  sink(log_connection, type = "message")

  on.exit(
    {
      while (sink.number(type = "message") > 2L) {
        sink(type = "message")
      }
      while (sink.number(type = "output") > 0L) {
        sink(type = "output")
      }
      close(log_connection)
    },
    add = TRUE
  )

  message("Analysis started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
  message("Project root: ", project_root)
  message("Input RData: ", input_rdata)
  message("Output root: ", output_root)
  message("Output RData: ", output_rdata)
  message("Log file: ", log_file)
  message("Integration method: ", integration_method)
  message("Configuration tag: ", configuration_tag)
  message("Graph names: ", paste(graph_names, collapse = ", "))


  # =============================================================================
  # 4. Packages and version checks
  # =============================================================================

  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
  })

  assert_seurat_v5()
  set.seed(seed)


  # =============================================================================
  # 5. Load parent object
  # =============================================================================

  parent_object <- load_single_seurat_object(input_rdata)

  message(
    "Parent object dimensions: ", nrow(parent_object), " genes x ",
    ncol(parent_object), " spots."
  )


  # =============================================================================
  # 6. RPCA integration and kNN/SNN graph construction
  # =============================================================================

  result <- run_integration_and_neighbors_graph(
    object = parent_object,
    integration_method = integration_method,
    dims = dims,
    k_param = k_param,
    prune_snn = prune_snn,
    assay = assay,
    original_reduction = original_reduction,
    integrated_reduction = integrated_reduction,
    features = SeuratObject::VariableFeatures(parent_object),
    sample_id_col = sample_id_col,
    expected_n_samples = expected_n_samples,
    normalization_method = normalization_method,
    integration_method_args = integration_method_args,
    graph_names = graph_names,
    nn_method = "annoy",
    n_trees = 50L,
    distance_metric = "euclidean",
    l2_norm = FALSE,
    seed = seed,
    overwrite_graphs = FALSE,
    verbose = TRUE
  )

  neighbors_object <- result$object

  rm(parent_object, result)
  invisible(gc())


  # =============================================================================
  # 7. Confirm that this remains a pre-clustering object
  # =============================================================================

  if (!integrated_reduction %in% names(neighbors_object@reductions)) {
    stop(
      "Final object is missing reduction `", integrated_reduction, "`.",
      call. = FALSE
    )
  }

  missing_graphs <- setdiff(graph_names, names(neighbors_object@graphs))

  if (length(missing_graphs) > 0L) {
    stop(
      "Final object is missing graph(s): ",
      paste(missing_graphs, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  clustering_columns <- grep(
    pattern = "(^seurat_clusters$|(^|_)res[._]?[0-9]+)",
    x = colnames(neighbors_object[[]]),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(clustering_columns) > 0L) {
    stop(
      "Unexpected clustering columns found: ",
      paste(clustering_columns, collapse = ", "),
      ". This stage must remain pre-clustering.",
      call. = FALSE
    )
  }


  # =============================================================================
  # 8. Save the integrated pre-clustering object
  # =============================================================================

  save_single_seurat_object(
    object = neighbors_object,
    object_name = output_object_name,
    output_file = output_rdata,
    compress = TRUE
  )

  message("Saved RPCA + neighbors graph object: ", output_rdata)
  message("Saved R object name: ", output_object_name)
  message(
    "Final object: ", nrow(neighbors_object), " genes x ",
    ncol(neighbors_object), " spots."
  )
  message(
    "Integrated reduction dimensions: ",
    ncol(Seurat::Embeddings(neighbors_object, reduction = integrated_reduction))
  )
  message("Stored graphs: ", paste(graph_names, collapse = ", "))
  message("No clustering or UMAP was performed.")
  message("Analysis completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))

  invisible(output_rdata)
}


tryCatch(
  main(),
  error = function(error_condition) {
    message("ERROR: ", conditionMessage(error_condition))
    quit(status = 1L, save = "no")
  }
)
