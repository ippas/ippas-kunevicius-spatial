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


canonical_expression_cluster_id <- function(x) {
  sub("^C", "", as.character(x), ignore.case = TRUE)
}

extract_expression_cluster_id <- function(file_name) {
  file_name <- basename(file_name)
  matched <- regexec(
    "_cluster_([^_]+)_rawCounts",
    file_name,
    perl = TRUE
  )
  pieces <- regmatches(file_name, matched)[[1]]

  if (length(pieces) != 2L) {
    stop(
      "Could not extract cluster ID from expression-matrix file: ",
      file_name,
      call. = FALSE
    )
  }

  canonical_expression_cluster_id(pieces[[2]])
}

sort_expression_cluster_ids <- function(cluster_ids) {
  cluster_ids <- as.character(cluster_ids)
  numeric_ids <- suppressWarnings(as.integer(cluster_ids))

  if (all(!is.na(numeric_ids))) {
    return(cluster_ids[order(numeric_ids)])
  }

  sort(cluster_ids)
}

read_per_cluster_expression_matrix <- function(file_name) {
  table <- utils::read.delim(
    file_name,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  metadata_columns <- c(
    "ensembl_gene_id",
    "gene",
    "chromosome"
  )

  missing_metadata_columns <- setdiff(
    metadata_columns,
    colnames(table)
  )

  if (length(missing_metadata_columns) > 0L) {
    stop(
      "Expression matrix is missing required metadata column(s): ",
      paste(missing_metadata_columns, collapse = ", "),
      "\nFile: ",
      file_name,
      call. = FALSE
    )
  }

  sample_columns <- setdiff(
    colnames(table),
    metadata_columns
  )

  if (length(sample_columns) == 0L) {
    stop(
      "No sample columns were found in expression matrix:\n",
      file_name,
      call. = FALSE
    )
  }

  matrix_values <- as.matrix(
    table[
      ,
      sample_columns,
      drop = FALSE
    ]
  )
  storage.mode(matrix_values) <- "double"

  if (anyNA(matrix_values) || any(!is.finite(matrix_values))) {
    stop(
      "Expression matrix contains NA or non-finite values:\n",
      file_name,
      call. = FALSE
    )
  }

  gene_ids <- as.character(table$ensembl_gene_id)

  if (
    anyNA(gene_ids) ||
      any(!nzchar(gene_ids)) ||
      anyDuplicated(gene_ids)
  ) {
    stop(
      "Expression matrix contains missing or duplicated Ensembl gene IDs:\n",
      file_name,
      call. = FALSE
    )
  }

  rownames(matrix_values) <- gene_ids

  list(
    matrix = matrix_values,
    gene_metadata = table[
      ,
      metadata_columns,
      drop = FALSE
    ],
    sample_ids = sample_columns
  )
}

load_per_cluster_expression_matrices <- function(
  raw_pseudobulk_counts_dir,
  log2_qn_expression_dir
) {
  if (!dir.exists(raw_pseudobulk_counts_dir)) {
    stop(
      "Raw pseudobulk-count directory does not exist:\n",
      raw_pseudobulk_counts_dir,
      call. = FALSE
    )
  }

  if (!dir.exists(log2_qn_expression_dir)) {
    stop(
      "log2(count + 1) + quantile-normalized directory does not exist:\n",
      log2_qn_expression_dir,
      call. = FALSE
    )
  }

  raw_files <- list.files(
    raw_pseudobulk_counts_dir,
    pattern = "_rawCounts_proteinCodingAutosomalAndX\\.tsv\\.gz$",
    full.names = TRUE
  )

  qn_files <- list.files(
    log2_qn_expression_dir,
    pattern = "_rawCounts_log2plus1Transformation_quantileNormalize_proteinCodingAutosomalAndX\\.tsv\\.gz$",
    full.names = TRUE
  )

  if (length(raw_files) == 0L) {
    stop(
      "No per-cluster raw pseudobulk-count files were found in:\n",
      raw_pseudobulk_counts_dir,
      call. = FALSE
    )
  }

  if (length(qn_files) == 0L) {
    stop(
      "No per-cluster log2+1 QN files were found in:\n",
      log2_qn_expression_dir,
      call. = FALSE
    )
  }

  raw_cluster_ids <- vapply(
    raw_files,
    extract_expression_cluster_id,
    character(1)
  )
  qn_cluster_ids <- vapply(
    qn_files,
    extract_expression_cluster_id,
    character(1)
  )

  if (anyDuplicated(raw_cluster_ids)) {
    stop(
      "Duplicated raw pseudobulk-count files were found for cluster(s): ",
      paste(unique(raw_cluster_ids[duplicated(raw_cluster_ids)]), collapse = ", "),
      call. = FALSE
    )
  }

  if (anyDuplicated(qn_cluster_ids)) {
    stop(
      "Duplicated log2+1 QN files were found for cluster(s): ",
      paste(unique(qn_cluster_ids[duplicated(qn_cluster_ids)]), collapse = ", "),
      call. = FALSE
    )
  }

  if (!setequal(raw_cluster_ids, qn_cluster_ids)) {
    stop(
      "Raw pseudobulk and log2+1 QN cluster sets do not match.\n",
      "Raw: ",
      paste(sort_expression_cluster_ids(raw_cluster_ids), collapse = ", "),
      "\nlog2+1 QN: ",
      paste(sort_expression_cluster_ids(qn_cluster_ids), collapse = ", "),
      call. = FALSE
    )
  }

  cluster_ids <- sort_expression_cluster_ids(raw_cluster_ids)

  raw_files <- raw_files[
    match(cluster_ids, raw_cluster_ids)
  ]
  qn_files <- qn_files[
    match(cluster_ids, qn_cluster_ids)
  ]

  raw_counts <- vector("list", length(cluster_ids))
  log2_qn <- vector("list", length(cluster_ids))
  gene_metadata <- vector("list", length(cluster_ids))

  names(raw_counts) <- cluster_ids
  names(log2_qn) <- cluster_ids
  names(gene_metadata) <- cluster_ids

  for (i in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[[i]]

    message(
      "Loading raw/log2+1 QN matrices for cluster C",
      cluster_id,
      " (",
      i,
      "/",
      length(cluster_ids),
      ")"
    )

    raw_loaded <- read_per_cluster_expression_matrix(
      raw_files[[i]]
    )
    qn_loaded <- read_per_cluster_expression_matrix(
      qn_files[[i]]
    )

    if (
      !identical(
        rownames(raw_loaded$matrix),
        rownames(qn_loaded$matrix)
      )
    ) {
      stop(
        "Raw and log2+1 QN gene order differs for cluster C",
        cluster_id,
        ".",
        call. = FALSE
      )
    }

    if (
      !identical(
        colnames(raw_loaded$matrix),
        colnames(qn_loaded$matrix)
      )
    ) {
      stop(
        "Raw and log2+1 QN sample order differs for cluster C",
        cluster_id,
        ".",
        call. = FALSE
      )
    }

    raw_gene_metadata <- raw_loaded$gene_metadata
    qn_gene_metadata <- qn_loaded$gene_metadata

    if (
      !identical(
        as.character(raw_gene_metadata$ensembl_gene_id),
        as.character(qn_gene_metadata$ensembl_gene_id)
      ) ||
        !identical(
          as.character(raw_gene_metadata$gene),
          as.character(qn_gene_metadata$gene)
        ) ||
        !identical(
          as.character(raw_gene_metadata$chromosome),
          as.character(qn_gene_metadata$chromosome)
        )
    ) {
      stop(
        "Raw and log2+1 QN gene metadata differs for cluster C",
        cluster_id,
        ".",
        call. = FALSE
      )
    }

    raw_counts[[cluster_id]] <- raw_loaded$matrix
    log2_qn[[cluster_id]] <- qn_loaded$matrix
    gene_metadata[[cluster_id]] <- raw_gene_metadata
  }

  list(
    raw_counts = raw_counts,
    log2_qn = log2_qn,
    gene_metadata = gene_metadata,
    cluster_ids = cluster_ids,
    raw_dir = raw_pseudobulk_counts_dir,
    log2_qn_dir = log2_qn_expression_dir
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
  raw_pseudobulk_counts_dir,
  log2_qn_expression_dir,
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

  per_cluster_expression <- load_per_cluster_expression_matrices(
    raw_pseudobulk_counts_dir = raw_pseudobulk_counts_dir,
    log2_qn_expression_dir = log2_qn_expression_dir
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
    per_cluster_expression = per_cluster_expression,
    genes = genes
  )
}
