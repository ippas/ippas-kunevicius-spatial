# ==============================================================================
# functions_integration_neighborsGraph.R
#
# Purpose:
# Reusable Seurat v5 functions for:
#   1. loading and saving a single Seurat object,
#   2. running one of the supported integration branches,
#   3. constructing a kNN/SNN graph from a selected dimensional reduction.
#
# Supported integration branches:
#   - noIntegration
#   - RPCA
#   - CCA
#   - Harmony
#   - FastMNN
#   - scVI
#
# This file performs no clustering and no UMAP calculation.
# ============================================================================== 


# ==============================================================================
# 1. General validation and I/O
# ==============================================================================

check_required_r_packages <- function(packages) {
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


assert_seurat_v5 <- function() {
  check_required_r_packages(c("Seurat", "SeuratObject"))

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


assert_seurat_object <- function(object) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }

  invisible(TRUE)
}


validate_positive_integer_vector <- function(x, argument_name) {
  if (
    length(x) == 0L || anyNA(x) || any(!is.finite(x)) ||
      any(x < 1) || any(x != as.integer(x))
  ) {
    stop(
      "`", argument_name, "` must contain positive integers.",
      call. = FALSE
    )
  }

  if (anyDuplicated(x) > 0L) {
    stop("`", argument_name, "` contains duplicated values.", call. = FALSE)
  }

  invisible(TRUE)
}


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
      "Expected exactly one Seurat object in `", rdata_file,
      "`, but found ", length(seurat_names), ".",
      call. = FALSE
    )
  }

  message("Loaded Seurat object: ", seurat_names)
  get(seurat_names, envir = input_environment)
}


save_single_seurat_object <- function(
    object,
    object_name,
    output_file,
    compress = TRUE
) {
  assert_seurat_object(object)

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
# 2. Naming helpers
# ==============================================================================

normalize_integration_method <- function(method) {
  if (length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("`method` must be one non-empty value.", call. = FALSE)
  }

  normalized <- tolower(gsub("[^a-zA-Z0-9]", "", method))

  aliases <- c(
    nointegration = "noIntegration",
    none = "noIntegration",
    unintegrated = "noIntegration",
    pca = "noIntegration",
    rpca = "rpca",
    cca = "cca",
    harmony = "harmony",
    fastmnn = "fastMNN",
    mnn = "fastMNN",
    scvi = "scVI"
  )

  if (!normalized %in% names(aliases)) {
    stop(
      "Unsupported integration method: `", method, "`. Supported values: ",
      "noIntegration, RPCA, CCA, Harmony, FastMNN, scVI.",
      call. = FALSE
    )
  }

  unname(aliases[[normalized]])
}


get_default_integrated_reduction <- function(method) {
  method <- normalize_integration_method(method)

  reductions <- c(
    noIntegration = "pca",
    rpca = "integrated.rpca",
    cca = "integrated.cca",
    harmony = "harmony",
    fastMNN = "integrated.mnn",
    scVI = "integrated.scvi"
  )

  unname(reductions[[method]])
}


format_dims_label <- function(dims) {
  validate_positive_integer_vector(dims, "dims")
  dims <- as.integer(dims)

  if (identical(dims, seq_len(max(dims)))) {
    return(paste0("dims", max(dims)))
  }

  paste0("dims", min(dims), "To", max(dims))
}


format_prune_snn_code <- function(prune_snn, multiplier = 1000L) {
  if (
    length(prune_snn) != 1L || is.na(prune_snn) ||
      !is.finite(prune_snn) || prune_snn < 0 || prune_snn > 1
  ) {
    stop("`prune_snn` must be one number between 0 and 1.", call. = FALSE)
  }

  code <- as.integer(round(prune_snn * multiplier))
  sprintf("%04d", code)
}


make_integration_graph_tag <- function(
    integration_method,
    dims,
    k_param,
    prune_snn
) {
  method <- normalize_integration_method(integration_method)
  validate_positive_integer_vector(dims, "dims")
  validate_positive_integer_vector(k_param, "k_param")

  if (length(k_param) != 1L) {
    stop("`k_param` must contain exactly one integer.", call. = FALSE)
  }

  paste0(
    method,
    "Dims", max(dims),
    "K", as.integer(k_param),
    "PruneSNN", format_prune_snn_code(prune_snn)
  )
}


make_graph_names <- function(configuration_tag) {
  if (
    length(configuration_tag) != 1L || is.na(configuration_tag) ||
      !nzchar(configuration_tag)
  ) {
    stop("`configuration_tag` must be one non-empty value.", call. = FALSE)
  }

  c(
    paste0(configuration_tag, "_nn"),
    paste0(configuration_tag, "_snn")
  )
}


# ==============================================================================
# 3. Integration method discovery and validation
# ==============================================================================

find_exported_function <- function(function_name, packages) {
  for (package_name in packages) {
    if (
      requireNamespace(package_name, quietly = TRUE) &&
        function_name %in% getNamespaceExports(package_name)
    ) {
      return(getExportedValue(package_name, function_name))
    }
  }

  stop(
    "Could not find exported function `", function_name,
    "` in package(s): ", paste(packages, collapse = ", "), ".",
    call. = FALSE
  )
}


get_integration_method_function <- function(method) {
  method <- normalize_integration_method(method)

  if (method == "noIntegration") {
    return(NULL)
  }

  if (method == "rpca") {
    return(find_exported_function("RPCAIntegration", "Seurat"))
  }

  if (method == "cca") {
    return(find_exported_function("CCAIntegration", "Seurat"))
  }

  if (method == "harmony") {
    check_required_r_packages("harmony")
    return(find_exported_function("HarmonyIntegration", "Seurat"))
  }

  if (method == "fastMNN") {
    check_required_r_packages(c("SeuratWrappers", "batchelor"))
    return(
      find_exported_function(
        "FastMNNIntegration",
        c("SeuratWrappers", "Seurat")
      )
    )
  }

  if (method == "scVI") {
    check_required_r_packages(c("SeuratWrappers", "reticulate"))
    return(
      find_exported_function(
        "scVIIntegration",
        c("SeuratWrappers", "Seurat")
      )
    )
  }

  stop("Internal error: unsupported integration method.", call. = FALSE)
}


validate_reduction_dimensions <- function(object, reduction, dims) {
  assert_seurat_object(object)
  validate_positive_integer_vector(dims, "dims")

  if (!reduction %in% names(object@reductions)) {
    stop(
      "The Seurat object does not contain reduction `",
      reduction,
      "`.",
      call. = FALSE
    )
  }

  available_dimensions <- ncol(
    Seurat::Embeddings(object, reduction = reduction)
  )

  if (max(dims) > available_dimensions) {
    stop(
      "Requested dimensions extend to ", max(dims),
      ", but reduction `", reduction,
      "` contains only ", available_dimensions, " dimensions.",
      call. = FALSE
    )
  }

  invisible(available_dimensions)
}


validate_parent_object_for_integration <- function(
    object,
    dims,
    assay = "RNA",
    original_reduction = "pca",
    sample_id_col = "sample_ID",
    expected_n_samples = NULL,
    require_sample_specific_layers = TRUE
) {
  assert_seurat_object(object)
  validate_positive_integer_vector(dims, "dims")

  if (!assay %in% names(object@assays)) {
    stop("The input object does not contain assay `", assay, "`.", call. = FALSE)
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
  sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]

  if (!is.null(expected_n_samples)) {
    validate_positive_integer_vector(expected_n_samples, "expected_n_samples")

    if (length(expected_n_samples) != 1L) {
      stop("`expected_n_samples` must contain one integer.", call. = FALSE)
    }

    if (length(sample_ids) != expected_n_samples) {
      stop(
        "Expected ", expected_n_samples, " samples, but found ",
        length(sample_ids), ".",
        call. = FALSE
      )
    }
  }

  validate_reduction_dimensions(object, original_reduction, dims)

  count_layers <- SeuratObject::Layers(object[[assay]], search = "^counts")
  data_layers <- SeuratObject::Layers(object[[assay]], search = "^data")

  if (
    isTRUE(require_sample_specific_layers) &&
      (length(count_layers) < 2L || length(data_layers) < 2L)
  ) {
    stop(
      "Integration requires sample-specific assay layers. Found ",
      length(count_layers), " count layer(s) and ",
      length(data_layers), " normalized data layer(s) in assay `",
      assay, "`.",
      call. = FALSE
    )
  }

  message(
    "Validated parent object: ", nrow(object), " genes x ",
    ncol(object), " spots; ", length(sample_ids), " samples."
  )
  message(
    "Assay layers: ", length(count_layers), " count layer(s), ",
    length(data_layers), " normalized data layer(s)."
  )

  invisible(TRUE)
}


# ==============================================================================
# 4. Universal integration function
# ==============================================================================

run_seurat_integration <- function(
    object,
    method,
    dims,
    assay = "RNA",
    original_reduction = "pca",
    new_reduction = NULL,
    features = NULL,
    sample_id_col = "sample_ID",
    expected_n_samples = NULL,
    normalization_method = "LogNormalize",
    method_args = list(),
    verbose = TRUE
) {
  assert_seurat_v5()

  method <- normalize_integration_method(method)
  validate_positive_integer_vector(dims, "dims")

  if (!is.list(method_args)) {
    stop("`method_args` must be a list.", call. = FALSE)
  }

  if (is.null(new_reduction)) {
    new_reduction <- get_default_integrated_reduction(method)
  }

  require_layers <- method != "noIntegration"

  validate_parent_object_for_integration(
    object = object,
    dims = dims,
    assay = assay,
    original_reduction = original_reduction,
    sample_id_col = sample_id_col,
    expected_n_samples = expected_n_samples,
    require_sample_specific_layers = require_layers
  )

  if (method == "noIntegration") {
    message(
      "No integration requested. Using reduction `",
      original_reduction,
      "` directly."
    )

    integration_info <- list(
      method = method,
      input_reduction = original_reduction,
      output_reduction = original_reduction,
      dims = as.integer(dims),
      assay = assay,
      normalization_method = normalization_method,
      method_args = method_args,
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
    )

    object@misc$seuratClusteringAnalysis$integration <- integration_info

    return(
      list(
        object = object,
        reduction = original_reduction,
        method = method,
        integration_info = integration_info
      )
    )
  }

  integration_function <- get_integration_method_function(method)

  if (is.null(features)) {
    features <- SeuratObject::VariableFeatures(object)
  }

  if (!is.character(features) || length(features) == 0L) {
    stop(
      "No integration features were supplied or found in VariableFeatures().",
      call. = FALSE
    )
  }

  protected_arguments <- c(
    "object", "method", "orig.reduction", "new.reduction", "assay", "features"
  )

  conflicting_arguments <- intersect(names(method_args), protected_arguments)

  if (length(conflicting_arguments) > 0L) {
    stop(
      "`method_args` must not redefine protected argument(s): ",
      paste(conflicting_arguments, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  integration_arguments <- list(
    object = object,
    method = integration_function,
    orig.reduction = original_reduction,
    new.reduction = new_reduction,
    assay = assay,
    features = features,
    verbose = verbose
  )

  if (method %in% c("rpca", "cca")) {
    integration_arguments$dims <- dims
    integration_arguments$normalization.method <- normalization_method
  }

  integration_arguments <- c(integration_arguments, method_args)

  message(
    "Running integration method `", method,
    "`; output reduction: `", new_reduction, "`..."
  )

  integrated_object <- do.call(
    what = Seurat::IntegrateLayers,
    args = integration_arguments
  )

  available_dimensions <- validate_reduction_dimensions(
    integrated_object,
    new_reduction,
    dims
  )

  integration_info <- list(
    method = method,
    input_reduction = original_reduction,
    output_reduction = new_reduction,
    requested_dims = as.integer(dims),
    available_output_dimensions = available_dimensions,
    assay = assay,
    normalization_method = normalization_method,
    nfeatures = length(features),
    method_args = method_args,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  )

  integrated_object@misc$seuratClusteringAnalysis$integration <- integration_info

  list(
    object = integrated_object,
    reduction = new_reduction,
    method = method,
    integration_info = integration_info
  )
}


# ==============================================================================
# 5. Universal kNN/SNN graph function
# ==============================================================================

build_seurat_neighbors_graph <- function(
    object,
    reduction,
    dims,
    k_param = 20L,
    prune_snn = 1 / 15,
    integration_method = "noIntegration",
    graph_names = NULL,
    nn_method = "annoy",
    n_trees = 50L,
    distance_metric = "euclidean",
    l2_norm = FALSE,
    seed = 7L,
    overwrite = FALSE,
    verbose = TRUE
) {
  assert_seurat_v5()
  assert_seurat_object(object)
  validate_positive_integer_vector(dims, "dims")
  validate_positive_integer_vector(k_param, "k_param")
  validate_positive_integer_vector(n_trees, "n_trees")
  validate_positive_integer_vector(seed, "seed")

  if (length(k_param) != 1L) {
    stop("`k_param` must contain exactly one integer.", call. = FALSE)
  }

  if (length(n_trees) != 1L) {
    stop("`n_trees` must contain exactly one integer.", call. = FALSE)
  }

  if (length(seed) != 1L) {
    stop("`seed` must contain exactly one integer.", call. = FALSE)
  }

  if (
    length(prune_snn) != 1L || is.na(prune_snn) ||
      !is.finite(prune_snn) || prune_snn < 0 || prune_snn > 1
  ) {
    stop("`prune_snn` must be one number between 0 and 1.", call. = FALSE)
  }

  if (k_param >= ncol(object)) {
    stop(
      "`k_param` must be smaller than the number of spots. Found k_param = ",
      k_param, " and ", ncol(object), " spots.",
      call. = FALSE
    )
  }

  validate_reduction_dimensions(object, reduction, dims)

  configuration_tag <- make_integration_graph_tag(
    integration_method = integration_method,
    dims = dims,
    k_param = k_param,
    prune_snn = prune_snn
  )

  if (is.null(graph_names)) {
    graph_names <- make_graph_names(configuration_tag)
  }

  if (
    !is.character(graph_names) || length(graph_names) != 2L ||
      anyNA(graph_names) || any(!nzchar(graph_names)) ||
      anyDuplicated(graph_names) > 0L
  ) {
    stop(
      "`graph_names` must contain exactly two unique, non-empty names: ",
      "first kNN, then SNN.",
      call. = FALSE
    )
  }

  existing_graphs <- intersect(graph_names, names(object@graphs))

  if (length(existing_graphs) > 0L && !isTRUE(overwrite)) {
    stop(
      "Graph(s) already exist and `overwrite = FALSE`: ",
      paste(existing_graphs, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (length(existing_graphs) > 0L && isTRUE(overwrite)) {
    for (graph_name in existing_graphs) {
      object[[graph_name]] <- NULL
    }
  }

  set.seed(seed)

  message(
    "Building kNN/SNN graphs from reduction `", reduction,
    "` using dims 1:", max(dims),
    ", k.param = ", k_param,
    ", prune.SNN = ", format(prune_snn, digits = 8), "..."
  )

  object <- Seurat::FindNeighbors(
    object = object,
    reduction = reduction,
    dims = dims,
    k.param = as.integer(k_param),
    compute.SNN = TRUE,
    prune.SNN = prune_snn,
    nn.method = nn_method,
    n.trees = as.integer(n_trees),
    annoy.metric = distance_metric,
    graph.name = graph_names,
    l2.norm = l2_norm,
    verbose = verbose
  )

  missing_graphs <- setdiff(graph_names, names(object@graphs))

  if (length(missing_graphs) > 0L) {
    stop(
      "FindNeighbors did not create expected graph(s): ",
      paste(missing_graphs, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  graph_info <- list(
    configuration_tag = configuration_tag,
    reduction = reduction,
    dims = as.integer(dims),
    k_param = as.integer(k_param),
    prune_snn = prune_snn,
    graph_names = graph_names,
    nn_method = nn_method,
    n_trees = as.integer(n_trees),
    distance_metric = distance_metric,
    l2_norm = l2_norm,
    seed = as.integer(seed),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  )

  object@misc$seuratClusteringAnalysis$neighborsGraph <- graph_info

  list(
    object = object,
    configuration_tag = configuration_tag,
    graph_names = graph_names,
    graph_info = graph_info
  )
}


# ==============================================================================
# 6. Convenience wrapper: integration followed by neighbor graph construction
# ==============================================================================

run_integration_and_neighbors_graph <- function(
    object,
    integration_method,
    dims,
    k_param = 20L,
    prune_snn = 1 / 15,
    assay = "RNA",
    original_reduction = "pca",
    integrated_reduction = NULL,
    features = NULL,
    sample_id_col = "sample_ID",
    expected_n_samples = NULL,
    normalization_method = "LogNormalize",
    integration_method_args = list(),
    graph_names = NULL,
    nn_method = "annoy",
    n_trees = 50L,
    distance_metric = "euclidean",
    l2_norm = FALSE,
    seed = 7L,
    overwrite_graphs = FALSE,
    verbose = TRUE
) {
  integration_result <- run_seurat_integration(
    object = object,
    method = integration_method,
    dims = dims,
    assay = assay,
    original_reduction = original_reduction,
    new_reduction = integrated_reduction,
    features = features,
    sample_id_col = sample_id_col,
    expected_n_samples = expected_n_samples,
    normalization_method = normalization_method,
    method_args = integration_method_args,
    verbose = verbose
  )

  graph_result <- build_seurat_neighbors_graph(
    object = integration_result$object,
    reduction = integration_result$reduction,
    dims = dims,
    k_param = k_param,
    prune_snn = prune_snn,
    integration_method = integration_result$method,
    graph_names = graph_names,
    nn_method = nn_method,
    n_trees = n_trees,
    distance_metric = distance_metric,
    l2_norm = l2_norm,
    seed = seed,
    overwrite = overwrite_graphs,
    verbose = verbose
  )

  list(
    object = graph_result$object,
    integration_method = integration_result$method,
    reduction = integration_result$reduction,
    configuration_tag = graph_result$configuration_tag,
    graph_names = graph_result$graph_names,
    integration_info = integration_result$integration_info,
    graph_info = graph_result$graph_info
  )
}
