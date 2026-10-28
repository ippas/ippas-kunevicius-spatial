# ==============================================================================
# Data loading
# ==============================================================================

load_single_seurat_object <- function(input_seurat_rdata_file) {
  if (!file.exists(input_seurat_rdata_file)) {
    stop(
      "Input Seurat RData file does not exist:\n",
      input_seurat_rdata_file,
      call. = FALSE
    )
  }

  load_environment <- new.env(parent = globalenv())

  loaded_names <- load(
    input_seurat_rdata_file,
    envir = load_environment
  )

  seurat_object_names <- loaded_names[
    vapply(
      loaded_names,
      function(object_name) {
        inherits(
          get(
            object_name,
            envir = load_environment,
            inherits = FALSE
          ),
          "Seurat"
        )
      },
      FUN.VALUE = logical(1)
    )
  ]

  if (length(seurat_object_names) != 1L) {
    stop(
      "Expected exactly one Seurat object in the RData file. Found: ",
      paste(seurat_object_names, collapse = ", "),
      call. = FALSE
    )
  }

  selected_name <- seurat_object_names[[1]]

  list(
    object = get(
      selected_name,
      envir = load_environment,
      inherits = FALSE
    ),
    object_name = selected_name
  )
}

load_edger_results_and_test_definitions <- function(
  edger_results_rdata_file
) {
  if (!file.exists(edger_results_rdata_file)) {
    stop(
      "edgeR results RData file does not exist:\n",
      edger_results_rdata_file,
      call. = FALSE
    )
  }

  result_environment <- new.env(parent = globalenv())

  loaded_names <- load(
    edger_results_rdata_file,
    envir = result_environment
  )

  required_objects <- c(
    "edgeR_perCluster_combinedResults",
    "edgeR_perCluster_testDefinitions"
  )

  missing_objects <- setdiff(
    required_objects,
    loaded_names
  )

  if (length(missing_objects) > 0L) {
    stop(
      "The edgeR results RData file is missing required object(s): ",
      paste(missing_objects, collapse = ", "),
      ". Loaded objects: ",
      paste(loaded_names, collapse = ", "),
      call. = FALSE
    )
  }

  combined_results <- get(
    "edgeR_perCluster_combinedResults",
    envir = result_environment,
    inherits = FALSE
  )

  test_definitions <- get(
    "edgeR_perCluster_testDefinitions",
    envir = result_environment,
    inherits = FALSE
  )

  if (!is.list(combined_results) || length(combined_results) == 0L) {
    stop(
      "edgeR_perCluster_combinedResults is empty or is not a list.",
      call. = FALSE
    )
  }

  if (!is.data.frame(test_definitions)) {
    stop(
      "edgeR_perCluster_testDefinitions is not a data.frame.",
      call. = FALSE
    )
  }

  required_definition_columns <- c(
    "test_id",
    "sheet_name",
    "comparison"
  )

  missing_definition_columns <- setdiff(
    required_definition_columns,
    colnames(test_definitions)
  )

  if (length(missing_definition_columns) > 0L) {
    stop(
      "edgeR_perCluster_testDefinitions is missing column(s): ",
      paste(missing_definition_columns, collapse = ", "),
      call. = FALSE
    )
  }

  list(
    combined_results = combined_results,
    test_definitions = test_definitions,
    loaded_names = loaded_names
  )
}


load_edger_model_objects <- function(
  edger_model_objects_rdata_file
) {
  if (!file.exists(edger_model_objects_rdata_file)) {
    stop(
      "edgeR model-object RData file does not exist:\n",
      edger_model_objects_rdata_file,
      call. = FALSE
    )
  }

  model_environment <- new.env(parent = globalenv())

  loaded_names <- load(
    edger_model_objects_rdata_file,
    envir = model_environment
  )

  required_object <- "edgeR_perCluster_modelObjects"

  if (!required_object %in% loaded_names) {
    stop(
      "The edgeR model-object RData file is missing `",
      required_object,
      "`. Loaded objects: ",
      paste(loaded_names, collapse = ", "),
      call. = FALSE
    )
  }

  model_objects <- get(
    required_object,
    envir = model_environment,
    inherits = FALSE
  )

  if (!is.list(model_objects) || length(model_objects) == 0L) {
    stop(
      "edgeR_perCluster_modelObjects is empty or is not a list.",
      call. = FALSE
    )
  }

  list(
    model_objects = model_objects,
    loaded_names = loaded_names
  )
}

get_expression_gene_names <- function(
  seurat_object,
  expression_assay
) {
  if (!expression_assay %in% SeuratObject::Assays(seurat_object)) {
    stop(
      "Expression assay is absent from the Seurat object: ",
      expression_assay,
      call. = FALSE
    )
  }

  genes <- rownames(
    seurat_object[[expression_assay]]
  )

  if (is.null(genes) || length(genes) == 0L) {
    stop(
      "No gene names were found in the selected expression assay.",
      call. = FALSE
    )
  }

  sort(
    unique(
      as.character(genes)
    )
  )
}

load_shiny_app_data <- function(
  input_seurat_rdata_file,
  edger_results_rdata_file,
  edger_model_objects_rdata_file,
  expression_assay
) {
  seurat_loaded <- load_single_seurat_object(
    input_seurat_rdata_file
  )

  edger_loaded <- load_edger_results_and_test_definitions(
    edger_results_rdata_file
  )

  edger_models_loaded <- load_edger_model_objects(
    edger_model_objects_rdata_file
  )

  genes <- get_expression_gene_names(
    seurat_object = seurat_loaded$object,
    expression_assay = expression_assay
  )

  list(
    seurat_object = seurat_loaded$object,
    seurat_object_name = seurat_loaded$object_name,
    edger_results = edger_loaded$combined_results,
    edger_test_definitions = edger_loaded$test_definitions,
    edger_model_objects = edger_models_loaded$model_objects,
    genes = genes
  )
}
