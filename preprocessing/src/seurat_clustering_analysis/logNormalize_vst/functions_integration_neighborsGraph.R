# ==============================================================================
# functions_integration_neighborsGraph.R
#
# Reusable helpers for the pre-clustering stage:
# - load/save one Seurat object from/to RData;
# - create stable analysis tags and graph names;
# - validate the PCA50 parent object;
# - run RPCA integration once for a selected dimensions setting;
# - build kNN/SNN graphs for a selected k.param and prune.SNN;
# - optionally run both operations in one call.
#
# This file performs no clustering, UMAP calculation, plotting or direct file
# writing other than through the explicit save helper.
# ==============================================================================


assert_seurat_v5 <- function() {

  seurat_version <- utils::packageVersion("Seurat")
  seurat_object_version <- utils::packageVersion("SeuratObject")

  message("Seurat version: ", as.character(seurat_version))
  message("SeuratObject version: ", as.character(seurat_object_version))

  if (seurat_version < "5.0.0") {
    stop(
      "Seurat v5 or newer is required. Current version: ",
      as.character(seurat_version),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


load_single_seurat_object <- function(rdata_file) {

  if (!file.exists(rdata_file)) {
    stop(
      "Input RData file does not exist: ",
      rdata_file,
      call. = FALSE
    )
  }

  input_environment <- new.env(parent = emptyenv())
  loaded_names <- load(rdata_file, envir = input_environment)

  seurat_names <- loaded_names[
    vapply(
      loaded_names,
      function(object_name) {
        inherits(
          get(object_name, envir = input_environment, inherits = FALSE),
          "Seurat"
        )
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
      ". Objects loaded: ",
      paste(loaded_names, collapse = ", "),
      call. = FALSE
    )
  }

  message("Loaded Seurat object: ", seurat_names[[1]])

  get(
    seurat_names[[1]],
    envir = input_environment,
    inherits = FALSE
  )
}


save_single_seurat_object <- function(
    object,
    object_name,
    output_file,
    compress = TRUE
) {

  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }

  if (
    !is.character(object_name) ||
      length(object_name) != 1L ||
      is.na(object_name) ||
      object_name == ""
  ) {
    stop("`object_name` must be one non-empty string.", call. = FALSE)
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  save_environment <- new.env(parent = emptyenv())
  assign(object_name, object, envir = save_environment)

  save(
    list = object_name,
    file = output_file,
    envir = save_environment,
    compress = compress
  )

  if (
    !file.exists(output_file) ||
      is.na(file.info(output_file)$size) ||
      file.info(output_file)$size <= 0L
  ) {
    stop(
      "RData file was not saved correctly: ",
      output_file,
      call. = FALSE
    )
  }

  invisible(output_file)
}


format_prune_snn_code <- function(prune_snn) {

  if (
    length(prune_snn) != 1L ||
      is.na(prune_snn) ||
      prune_snn < 0 ||
      prune_snn > 1
  ) {
    stop(
      "`prune_snn` must be one numeric value between 0 and 1.",
      call. = FALSE
    )
  }

  sprintf("%04d", round(as.numeric(prune_snn) * 1000))
}


format_dims_label <- function(dims) {

  if (length(dims) == 0L || anyNA(dims) || any(dims < 1L)) {
    stop("`dims` must contain positive integer dimensions.", call. = FALSE)
  }

  paste0("dims", max(as.integer(dims)))
}


make_integration_graph_tag <- function(
    integration_method,
    dims,
    k_param,
    prune_snn
) {

  integration_method <- tolower(as.character(integration_method))

  method_tag <- switch(
    integration_method,
    rpca = "rpca",
    stop(
      "Unsupported integration method: ",
      integration_method,
      ". This workflow currently supports `rpca`.",
      call. = FALSE
    )
  )

  paste0(
    method_tag,
    "Dims",
    max(as.integer(dims)),
    "K",
    as.integer(k_param),
    "PruneSNN",
    format_prune_snn_code(prune_snn)
  )
}


make_graph_names <- function(configuration_tag) {

  c(
    paste0(configuration_tag, "_nn"),
    paste0(configuration_tag, "_snn")
  )
}


validate_parent_object_for_integration <- function(
    object,
    dims,
    assay = "RNA",
    original_reduction = "pca",
    sample_id_col = "sample_ID",
    expected_n_samples = 16L
) {

  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }

  if (!assay %in% names(object@assays)) {
    stop(
      "The input object does not contain assay `",
      assay,
      "`.",
      call. = FALSE
    )
  }

  if (!sample_id_col %in% colnames(object[[]])) {
    stop(
      "The input object does not contain metadata column `",
      sample_id_col,
      "`.",
      call. = FALSE
    )
  }

  sample_ids <- unique(as.character(object[[sample_id_col]][, 1]))

  if (length(sample_ids) != expected_n_samples) {
    stop(
      "Expected ",
      expected_n_samples,
      " samples, but found ",
      length(sample_ids),
      ".",
      call. = FALSE
    )
  }

  if (!original_reduction %in% SeuratObject::Reductions(object)) {
    stop(
      "The input object does not contain reduction `",
      original_reduction,
      "`.",
      call. = FALSE
    )
  }

  available_dimensions <- ncol(
    SeuratObject::Embeddings(object[[original_reduction]])
  )

  if (max(dims) > available_dimensions) {
    stop(
      "Requested dimensions extend to ",
      max(dims),
      ", but reduction `",
      original_reduction,
      "` contains only ",
      available_dimensions,
      " dimensions.",
      call. = FALSE
    )
  }

  count_layers <- SeuratObject::Layers(
    object[[assay]],
    search = "^counts"
  )

  data_layers <- SeuratObject::Layers(
    object[[assay]],
    search = "^data"
  )

  if (length(count_layers) < 2L || length(data_layers) < 2L) {
    stop(
      "RPCA integration requires sample-specific assay layers. Found ",
      length(count_layers),
      " count layer(s) and ",
      length(data_layers),
      " normalized data layer(s).",
      call. = FALSE
    )
  }

  message(
    "Validated parent object: ",
    nrow(object),
    " genes x ",
    ncol(object),
    " spots; ",
    length(sample_ids),
    " samples; ",
    available_dimensions,
    " dimensions in `",
    original_reduction,
    "`."
  )

  invisible(TRUE)
}


validate_integrated_object <- function(
    object,
    dims,
    integrated_reduction = "integrated.rpca"
) {

  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }

  if (!integrated_reduction %in% SeuratObject::Reductions(object)) {
    stop(
      "Integrated reduction `",
      integrated_reduction,
      "` is absent from the object.",
      call. = FALSE
    )
  }

  integrated_embeddings <- SeuratObject::Embeddings(
    object[[integrated_reduction]]
  )

  if (ncol(integrated_embeddings) < max(dims)) {
    stop(
      "Integrated reduction contains ",
      ncol(integrated_embeddings),
      " dimensions; at least ",
      max(dims),
      " are required.",
      call. = FALSE
    )
  }

  if (
    nrow(integrated_embeddings) != ncol(object) ||
      !identical(rownames(integrated_embeddings), colnames(object))
  ) {
    stop(
      "Integrated embeddings are not aligned with the Seurat object.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


initialise_integration_provenance <- function(object) {

  if (is.null(object@misc$integrationNeighborsGraph)) {
    object@misc$integrationNeighborsGraph <- list()
  }

  object
}


run_integration_only <- function(
    object,
    integration_method = "rpca",
    dims = 1:20,
    assay = "RNA",
    original_reduction = "pca",
    integrated_reduction = "integrated.rpca",
    features = NULL,
    sample_id_col = "sample_ID",
    expected_n_samples = 16L,
    normalization_method = "LogNormalize",
    integration_method_args = list(),
    seed = 7L,
    verbose = TRUE
) {

  integration_method <- tolower(as.character(integration_method))
  dims <- as.integer(dims)

  if (!identical(integration_method, "rpca")) {
    stop(
      "Only `rpca` integration is supported by this workflow.",
      call. = FALSE
    )
  }

  validate_parent_object_for_integration(
    object = object,
    dims = dims,
    assay = assay,
    original_reduction = original_reduction,
    sample_id_col = sample_id_col,
    expected_n_samples = expected_n_samples
  )

  if (is.null(features)) {
    features <- SeuratObject::VariableFeatures(object[[assay]])
  }

  if (length(features) == 0L) {
    stop("No variable features were supplied for integration.", call. = FALSE)
  }

  SeuratObject::DefaultAssay(object) <- assay

  if (integrated_reduction %in% SeuratObject::Reductions(object)) {
    message(
      "Removing pre-existing reduction `",
      integrated_reduction,
      "` before integration."
    )
    object[[integrated_reduction]] <- NULL
  }

  set.seed(as.integer(seed))
  integration_started <- Sys.time()

  message(
    "Running RPCA integration using dimensions 1:",
    max(dims),
    "."
  )

  integration_arguments <- c(
    list(
      object = object,
      method = Seurat::RPCAIntegration,
      orig.reduction = original_reduction,
      new.reduction = integrated_reduction,
      features = features,
      dims = dims,
      normalization.method = normalization_method,
      verbose = verbose
    ),
    integration_method_args
  )

  object <- do.call(
    Seurat::IntegrateLayers,
    integration_arguments
  )

  integration_finished <- Sys.time()

  validate_integrated_object(
    object = object,
    dims = dims,
    integrated_reduction = integrated_reduction
  )

  object <- initialise_integration_provenance(object)

  provenance_name <- paste0(
    "rpcaDims",
    max(dims),
    "Integration"
  )

  object@misc$integrationNeighborsGraph[[provenance_name]] <- list(
    stage = "integration",
    integrationMethod = integration_method,
    assay = assay,
    originalReduction = original_reduction,
    integratedReduction = integrated_reduction,
    dims = dims,
    nFeatures = length(features),
    normalizationMethod = normalization_method,
    sampleIdColumn = sample_id_col,
    expectedNSamples = as.integer(expected_n_samples),
    seed = as.integer(seed),
    elapsedSeconds = as.numeric(
      difftime(integration_finished, integration_started, units = "secs")
    ),
    completedAt = format(
      integration_finished,
      "%Y-%m-%d %H:%M:%S %z"
    )
  )

  message(
    "RPCA integration completed: ",
    ncol(SeuratObject::Embeddings(object[[integrated_reduction]])),
    " dimensions stored in `",
    integrated_reduction,
    "`."
  )

  object
}


build_neighbors_graph_from_integrated <- function(
    object,
    dims = 1:20,
    k_param = 20L,
    prune_snn = 1 / 15,
    integrated_reduction = "integrated.rpca",
    graph_names,
    nn_method = "annoy",
    n_trees = 50L,
    distance_metric = "euclidean",
    l2_norm = FALSE,
    seed = 7L,
    overwrite_graphs = FALSE,
    verbose = TRUE
) {

  dims <- as.integer(dims)

  if (
    length(k_param) != 1L ||
      is.na(k_param) ||
      k_param < 2L
  ) {
    stop("`k_param` must be one integer greater than 1.", call. = FALSE)
  }

  if (
    length(prune_snn) != 1L ||
      is.na(prune_snn) ||
      prune_snn < 0 ||
      prune_snn > 1
  ) {
    stop(
      "`prune_snn` must be one numeric value between 0 and 1.",
      call. = FALSE
    )
  }

  if (length(graph_names) != 2L || anyNA(graph_names) || any(graph_names == "")) {
    stop(
      "`graph_names` must contain exactly two non-empty graph names.",
      call. = FALSE
    )
  }

  validate_integrated_object(
    object = object,
    dims = dims,
    integrated_reduction = integrated_reduction
  )

  existing_graphs <- intersect(graph_names, SeuratObject::Graphs(object))

  if (length(existing_graphs) > 0L && !isTRUE(overwrite_graphs)) {
    stop(
      "Expected graph name(s) already exist: ",
      paste(existing_graphs, collapse = ", "),
      ". Set `overwrite_graphs = TRUE` to replace them.",
      call. = FALSE
    )
  }

  if (length(existing_graphs) > 0L) {
    for (graph_name in existing_graphs) {
      object[[graph_name]] <- NULL
    }
  }

  set.seed(as.integer(seed))
  graph_started <- Sys.time()

  message(
    "Building kNN/SNN graphs: dims=1:",
    max(dims),
    "; k.param=",
    as.integer(k_param),
    "; prune.SNN=",
    formatC(prune_snn, format = "f", digits = 6L),
    "."
  )

  object <- Seurat::FindNeighbors(
    object = object,
    reduction = integrated_reduction,
    dims = dims,
    k.param = as.integer(k_param),
    compute.SNN = TRUE,
    prune.SNN = as.numeric(prune_snn),
    nn.method = nn_method,
    n.trees = as.integer(n_trees),
    annoy.metric = distance_metric,
    l2.norm = l2_norm,
    graph.name = graph_names,
    verbose = verbose
  )

  graph_finished <- Sys.time()

  missing_graphs <- setdiff(graph_names, SeuratObject::Graphs(object))

  if (length(missing_graphs) > 0L) {
    stop(
      "FindNeighbors did not create expected graph(s): ",
      paste(missing_graphs, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  object <- initialise_integration_provenance(object)

  graph_provenance_name <- paste0(
    "dims",
    max(dims),
    "K",
    as.integer(k_param),
    "Prune",
    format_prune_snn_code(prune_snn),
    "Graph"
  )

  object@misc$integrationNeighborsGraph[[graph_provenance_name]] <- list(
    stage = "neighborsGraph",
    integratedReduction = integrated_reduction,
    dims = dims,
    kParam = as.integer(k_param),
    pruneSNN = as.numeric(prune_snn),
    graphNames = graph_names,
    nnMethod = nn_method,
    nTrees = as.integer(n_trees),
    distanceMetric = distance_metric,
    l2Norm = l2_norm,
    seed = as.integer(seed),
    elapsedSeconds = as.numeric(
      difftime(graph_finished, graph_started, units = "secs")
    ),
    completedAt = format(
      graph_finished,
      "%Y-%m-%d %H:%M:%S %z"
    )
  )

  message(
    "Created graphs: ",
    paste(graph_names, collapse = ", "),
    "."
  )

  object
}


run_integration_and_neighbors_graph <- function(
    object,
    integration_method = "rpca",
    dims = 1:20,
    k_param = 20L,
    prune_snn = 1 / 15,
    assay = "RNA",
    original_reduction = "pca",
    integrated_reduction = "integrated.rpca",
    features = NULL,
    sample_id_col = "sample_ID",
    expected_n_samples = 16L,
    normalization_method = "LogNormalize",
    integration_method_args = list(),
    graph_names,
    nn_method = "annoy",
    n_trees = 50L,
    distance_metric = "euclidean",
    l2_norm = FALSE,
    seed = 7L,
    overwrite_graphs = FALSE,
    verbose = TRUE
) {

  integrated_object <- run_integration_only(
    object = object,
    integration_method = integration_method,
    dims = dims,
    assay = assay,
    original_reduction = original_reduction,
    integrated_reduction = integrated_reduction,
    features = features,
    sample_id_col = sample_id_col,
    expected_n_samples = expected_n_samples,
    normalization_method = normalization_method,
    integration_method_args = integration_method_args,
    seed = seed,
    verbose = verbose
  )

  neighbors_object <- build_neighbors_graph_from_integrated(
    object = integrated_object,
    dims = dims,
    k_param = k_param,
    prune_snn = prune_snn,
    integrated_reduction = integrated_reduction,
    graph_names = graph_names,
    nn_method = nn_method,
    n_trees = n_trees,
    distance_metric = distance_metric,
    l2_norm = l2_norm,
    seed = seed,
    overwrite_graphs = overwrite_graphs,
    verbose = verbose
  )

  list(
    object = neighbors_object,
    integratedReduction = integrated_reduction,
    graphNames = graph_names
  )
}
