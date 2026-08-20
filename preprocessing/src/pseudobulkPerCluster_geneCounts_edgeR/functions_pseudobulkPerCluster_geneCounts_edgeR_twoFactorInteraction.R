#!/usr/bin/env Rscript

# ==============================================================================
# functions_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction_FIXED_20260803.R
#
# Dedicated functions for one analysis workflow:
#   pseudobulk raw gene counts per sample and cluster followed by edgeR
#   quasi-likelihood modelling with:
#
#     expression ~ fmt_donor_group * sex
#
# The code is intentionally dedicated to this workflow. It does not modify or
# reuse the existing whole-tissue pseudobulk functions.
#
# Tested result sets per cluster:
#   1. Interaction
#   2. Overall sex effect
#   3. Overall donor-group effect
#   4-9. All six pairwise comparisons between the four sex-by-group cells
#
# Gene universe:
#   Ensembl 115 Mus musculus protein-coding genes on chromosomes 1-19 and X.
#   Chromosome Y, mitochondrial and non-canonical genes are excluded before
#   filterByExpr and edgeR modelling.
#
# Additional descriptive metrics added to every full-result table:
#   - mean and sample SD of the sample-level percent of positive spots;
#   - mean and sample SD of sample-level TMM-normalized CPM.
#
# Most tests receive eight descriptive columns:
#   two compared biological levels x four metrics.
# The interaction receives sixteen descriptive columns:
#   four sex-by-group cells x four metrics.
#
# These metrics are descriptive only. They do not alter filterByExpr, TMM
# normalization, dispersion estimation, model fitting, contrasts, P-values or
# FDR values.
# ==============================================================================

PSEUDOBULK_EDGER_FUNCTIONS_BUILD <- "FIXED_20260803_V1"

check_required_packages_pseudobulk_per_cluster <- function() {
  required_packages <- c(
    "Seurat",
    "SeuratObject",
    "Matrix",
    "edgeR",
    "limma",
    "dplyr",
    "tibble",
    "openxlsx"
  )

  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing_packages) > 0L) {
    stop(
      "Missing required packages: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

write_tsv_pseudobulk_per_cluster <- function(object, output_file) {
  output_connection <- if (grepl("\\.gz$", output_file, ignore.case = TRUE)) {
    gzfile(output_file, open = "wt")
  } else {
    file(output_file, open = "wt")
  }

  on.exit(close(output_connection), add = TRUE)

  write.table(
    object,
    file = output_connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )

  invisible(output_file)
}

sanitize_file_component_pseudobulk_per_cluster <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)

  ifelse(x == "", "unnamed", x)
}

natural_order_cluster_ids_pseudobulk_per_cluster <- function(cluster_ids) {
  cluster_ids <- unique(as.character(cluster_ids))
  numeric_ids <- suppressWarnings(as.numeric(cluster_ids))

  if (all(!is.na(numeric_ids))) {
    return(cluster_ids[order(numeric_ids)])
  }

  cluster_ids[order(cluster_ids)]
}

load_single_seurat_object_from_rdata <- function(
    input_rdata_file,
    requested_object_name = NULL
) {
  if (!file.exists(input_rdata_file)) {
    stop("Input RData file does not exist: ", input_rdata_file, call. = FALSE)
  }

  load_environment <- new.env(parent = emptyenv())
  loaded_object_names <- load(input_rdata_file, envir = load_environment)

  if (length(loaded_object_names) == 0L) {
    stop("No objects were loaded from: ", input_rdata_file, call. = FALSE)
  }

  if (!is.null(requested_object_name)) {
    if (!requested_object_name %in% loaded_object_names) {
      stop(
        "Requested Seurat object '", requested_object_name,
        "' was not found in the RData file. Available objects: ",
        paste(loaded_object_names, collapse = ", "),
        call. = FALSE
      )
    }

    requested_object <- get(requested_object_name, envir = load_environment)

    if (!inherits(requested_object, "Seurat")) {
      stop(
        "Requested object '", requested_object_name,
        "' is not a Seurat object.",
        call. = FALSE
      )
    }

    return(
      list(
        object = requested_object,
        object_name = requested_object_name,
        loaded_object_names = loaded_object_names
      )
    )
  }

  seurat_object_names <- loaded_object_names[
    vapply(
      loaded_object_names,
      function(object_name) {
        inherits(get(object_name, envir = load_environment), "Seurat")
      },
      FUN.VALUE = logical(1)
    )
  ]

  if (length(seurat_object_names) == 0L) {
    stop(
      "No Seurat object was found in the RData file. Loaded objects: ",
      paste(loaded_object_names, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(seurat_object_names) > 1L) {
    stop(
      "More than one Seurat object was found in the RData file: ",
      paste(seurat_object_names, collapse = ", "),
      ". Set seurat_object_name explicitly in the runner script.",
      call. = FALSE
    )
  }

  selected_name <- seurat_object_names[[1]]

  list(
    object = get(selected_name, envir = load_environment),
    object_name = selected_name,
    loaded_object_names = loaded_object_names
  )
}

resolve_cluster_column_pseudobulk_per_cluster <- function(
    seurat_metadata,
    requested_cluster_column = NULL,
    expected_number_of_clusters = NULL
) {
  if (!is.null(requested_cluster_column)) {
    if (!requested_cluster_column %in% colnames(seurat_metadata)) {
      stop(
        "Requested cluster column does not exist: ",
        requested_cluster_column,
        "\nAvailable metadata columns:\n",
        paste(colnames(seurat_metadata), collapse = "\n"),
        call. = FALSE
      )
    }

    return(requested_cluster_column)
  }

  candidate_names <- colnames(seurat_metadata)[
    grepl("cluster|leiden", colnames(seurat_metadata), ignore.case = TRUE)
  ]

  if (length(candidate_names) == 0L) {
    stop(
      "No cluster-like metadata column was found. Set cluster_column explicitly ",
      "in the runner script.",
      call. = FALSE
    )
  }

  candidate_summary <- data.frame(
    column = candidate_names,
    n_clusters = vapply(
      candidate_names,
      function(column_name) {
        length(unique(stats::na.omit(seurat_metadata[[column_name]])))
      },
      FUN.VALUE = integer(1)
    ),
    stringsAsFactors = FALSE
  )

  if (!is.null(expected_number_of_clusters)) {
    candidate_summary <- candidate_summary[
      candidate_summary$n_clusters == expected_number_of_clusters,
      ,
      drop = FALSE
    ]
  }

  if (nrow(candidate_summary) == 0L) {
    stop(
      "No cluster-like metadata column has the expected number of clusters. ",
      "Set cluster_column explicitly in the runner script.",
      call. = FALSE
    )
  }

  preferred_resolution <- grepl(
    "res[^0-9]*0?[._-]?4|040|0[._-]?4",
    candidate_summary$column,
    ignore.case = TRUE
  )

  preferred_leiden <- grepl(
    "leiden",
    candidate_summary$column,
    ignore.case = TRUE
  )

  preferred_candidates <- candidate_summary[
    preferred_resolution & preferred_leiden,
    ,
    drop = FALSE
  ]

  if (nrow(preferred_candidates) == 1L) {
    return(preferred_candidates$column[[1]])
  }

  if (nrow(candidate_summary) == 1L) {
    return(candidate_summary$column[[1]])
  }

  stop(
    "Cluster column could not be selected unambiguously. Set cluster_column ",
    "explicitly in the runner script. Candidates:\n",
    paste(
      paste0(
        candidate_summary$column,
        " (", candidate_summary$n_clusters, " clusters)"
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

extract_raw_counts_from_seurat_layers <- function(
    seurat_object,
    assay_name
) {
  if (!assay_name %in% names(seurat_object@assays)) {
    stop(
      "Assay '", assay_name, "' is not present in the Seurat object. Available assays: ",
      paste(names(seurat_object@assays), collapse = ", "),
      call. = FALSE
    )
  }

  assay_object <- seurat_object[[assay_name]]

  count_layer_names <- tryCatch(
    SeuratObject::Layers(
      object = assay_object,
      search = "^counts"
    ),
    error = function(e) character(0)
  )

  if (length(count_layer_names) == 0L) {
    count_matrix <- tryCatch(
      Seurat::GetAssayData(
        object = seurat_object,
        assay = assay_name,
        slot = "counts"
      ),
      error = function(e) NULL
    )

    if (is.null(count_matrix)) {
      stop(
        "No raw count layer could be extracted from assay '",
        assay_name,
        "'.",
        call. = FALSE
      )
    }

    if (!inherits(count_matrix, "Matrix")) {
      count_matrix <- Matrix::Matrix(count_matrix, sparse = TRUE)
    }

    return(
      list(
        counts = count_matrix,
        count_layers = "counts"
      )
    )
  }

  count_layer_matrices <- lapply(
    count_layer_names,
    function(layer_name) {
      layer_matrix <- SeuratObject::LayerData(
        object = seurat_object,
        assay = assay_name,
        layer = layer_name
      )

      if (!inherits(layer_matrix, "Matrix")) {
        layer_matrix <- Matrix::Matrix(layer_matrix, sparse = TRUE)
      }

      layer_matrix
    }
  )

  reference_features <- rownames(count_layer_matrices[[1]])

  if (is.null(reference_features)) {
    stop("Feature names are missing from the raw count layers.", call. = FALSE)
  }

  for (layer_index in seq_along(count_layer_matrices)) {
    current_features <- rownames(count_layer_matrices[[layer_index]])

    if (!setequal(current_features, reference_features)) {
      stop(
        "Raw count layers do not contain identical feature sets. Problematic layer: ",
        count_layer_names[[layer_index]],
        call. = FALSE
      )
    }

    count_layer_matrices[[layer_index]] <- count_layer_matrices[[layer_index]][
      reference_features,
      ,
      drop = FALSE
    ]
  }

  count_matrix <- Reduce(
    function(left_matrix, right_matrix) {
      cbind(left_matrix, right_matrix)
    },
    count_layer_matrices
  )

  if (anyDuplicated(colnames(count_matrix))) {
    duplicated_cells <- unique(
      colnames(count_matrix)[duplicated(colnames(count_matrix))]
    )

    stop(
      "Duplicated spot names were found after joining raw count layers: ",
      paste(head(duplicated_cells, 20L), collapse = ", "),
      call. = FALSE
    )
  }

  list(
    counts = count_matrix,
    count_layers = count_layer_names
  )
}

read_and_validate_ensembl115_annotation <- function(
    annotation_rds_file,
    annotation_metadata_file,
    ensembl_release = 115L,
    genome_assembly = "GRCm39",
    analysis_chromosomes = c(as.character(1:19), "X")
) {
  if (!file.exists(annotation_rds_file)) {
    stop("Annotation RDS file does not exist: ", annotation_rds_file, call. = FALSE)
  }

  if (!file.exists(annotation_metadata_file)) {
    stop(
      "Annotation metadata file does not exist: ",
      annotation_metadata_file,
      call. = FALSE
    )
  }

  annotation_metadata <- read.delim(
    annotation_metadata_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_metadata_columns <- c(
    "genome_assembly",
    "ensembl_release",
    "retained_feature_type",
    "retained_gene_biotype",
    "number_of_protein_coding_genes"
  )

  missing_metadata_columns <- setdiff(
    required_metadata_columns,
    colnames(annotation_metadata)
  )

  if (length(missing_metadata_columns) > 0L) {
    stop(
      "Missing annotation metadata columns: ",
      paste(missing_metadata_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(annotation_metadata) != 1L) {
    stop("Annotation metadata must contain exactly one row.", call. = FALSE)
  }

  if (as.integer(annotation_metadata$ensembl_release[[1]]) != ensembl_release) {
    stop(
      "Unexpected Ensembl release in annotation metadata: ",
      annotation_metadata$ensembl_release[[1]],
      ". Expected: ", ensembl_release,
      call. = FALSE
    )
  }

  if (!identical(
    as.character(annotation_metadata$genome_assembly[[1]]),
    genome_assembly
  )) {
    stop(
      "Unexpected genome assembly in annotation metadata: ",
      annotation_metadata$genome_assembly[[1]],
      ". Expected: ", genome_assembly,
      call. = FALSE
    )
  }

  annotation_all <- readRDS(annotation_rds_file) |>
    tibble::as_tibble()

  required_annotation_columns <- c(
    "chromosome",
    "ensembl_gene_id",
    "gene",
    "gene_biotype"
  )

  missing_annotation_columns <- setdiff(
    required_annotation_columns,
    colnames(annotation_all)
  )

  if (length(missing_annotation_columns) > 0L) {
    stop(
      "Missing columns in annotation RDS: ",
      paste(missing_annotation_columns, collapse = ", "),
      call. = FALSE
    )
  }

  annotation_all <- annotation_all |>
    dplyr::transmute(
      chromosome = trimws(as.character(chromosome)),
      ensembl_gene_id = sub(
        "\\.[0-9]+$",
        "",
        trimws(as.character(ensembl_gene_id))
      ),
      gene = trimws(as.character(gene)),
      gene_biotype = trimws(as.character(gene_biotype))
    ) |>
    dplyr::filter(
      !is.na(chromosome),
      chromosome != "",
      !is.na(ensembl_gene_id),
      ensembl_gene_id != "",
      gene_biotype == "protein_coding"
    ) |>
    dplyr::mutate(
      gene = dplyr::if_else(
        is.na(gene) | gene == "",
        ensembl_gene_id,
        gene
      )
    ) |>
    dplyr::arrange(chromosome, ensembl_gene_id, gene) |>
    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

  if (nrow(annotation_all) == 0L) {
    stop("The annotation contains no protein-coding genes.", call. = FALSE)
  }

  expected_number_of_genes <- as.integer(
    annotation_metadata$number_of_protein_coding_genes[[1]]
  )

  if (nrow(annotation_all) != expected_number_of_genes) {
    stop(
      "Annotation RDS contains ", nrow(annotation_all),
      " protein-coding genes, while metadata reports ",
      expected_number_of_genes, ".",
      call. = FALSE
    )
  }

  annotation_analysis <- annotation_all |>
    dplyr::filter(chromosome %in% analysis_chromosomes) |>
    dplyr::mutate(
      chromosome = factor(chromosome, levels = analysis_chromosomes)
    ) |>
    dplyr::arrange(chromosome, ensembl_gene_id) |>
    dplyr::mutate(chromosome = as.character(chromosome))

  if (nrow(annotation_analysis) == 0L) {
    stop(
      "No protein-coding genes on chromosomes 1-19 or X remain.",
      call. = FALSE
    )
  }

  list(
    annotation_all = annotation_all,
    annotation_analysis = annotation_analysis,
    annotation_metadata = annotation_metadata
  )
}

map_seurat_features_to_analysis_annotation <- function(
    count_matrix,
    annotation_all,
    annotation_analysis
) {
  feature_names <- rownames(count_matrix)

  if (is.null(feature_names)) {
    stop("Raw count matrix has no feature names.", call. = FALSE)
  }

  if (anyDuplicated(feature_names)) {
    stop("Raw count matrix contains duplicated feature names.", call. = FALSE)
  }

  normalized_feature_names <- sub("\\.[0-9]+$", "", feature_names)
  ensembl_fraction <- mean(grepl("^ENSMUSG[0-9]+$", normalized_feature_names))
  feature_identifier_type <- if (ensembl_fraction >= 0.5) {
    "ensembl_gene_id"
  } else {
    "gene_symbol"
  }

  if (feature_identifier_type == "ensembl_gene_id") {
    feature_mapping <- tibble::tibble(
      feature_id = feature_names,
      normalized_feature_id = normalized_feature_names
    ) |>
      dplyr::left_join(
        annotation_all,
        by = c("normalized_feature_id" = "ensembl_gene_id")
      ) |>
      dplyr::rename(ensembl_gene_id = normalized_feature_id)
  } else {
    unique_symbol_annotation <- annotation_all |>
      dplyr::add_count(gene, name = "symbol_annotation_count") |>
      dplyr::filter(symbol_annotation_count == 1L) |>
      dplyr::select(-symbol_annotation_count)

    feature_mapping <- tibble::tibble(
      feature_id = feature_names,
      feature_gene_symbol = feature_names
    ) |>
      dplyr::left_join(
        unique_symbol_annotation,
        by = c("feature_gene_symbol" = "gene")
      ) |>
      dplyr::mutate(
        gene = feature_gene_symbol
      ) |>
      dplyr::select(
        feature_id,
        chromosome,
        ensembl_gene_id,
        gene,
        gene_biotype
      )
  }

  analysis_gene_ids <- annotation_analysis$ensembl_gene_id

  feature_mapping <- feature_mapping |>
    dplyr::mutate(
      found_in_ensembl115_protein_coding_annotation = !is.na(ensembl_gene_id),
      located_on_autosomal_or_X_chromosome =
        !is.na(ensembl_gene_id) & ensembl_gene_id %in% analysis_gene_ids,
      passed_protein_coding_autosomal_and_X_filter =
        found_in_ensembl115_protein_coding_annotation &
        located_on_autosomal_or_X_chromosome,
      mapping_status = dplyr::case_when(
        passed_protein_coding_autosomal_and_X_filter ~
          "passed_protein_coding_autosomal_and_X_filter",
        found_in_ensembl115_protein_coding_annotation ~
          "removed_protein_coding_gene_outside_autosomal_and_X",
        TRUE ~
          "removed_not_in_unambiguous_Ensembl115_protein_coding_annotation"
      )
    )

  retained_mapping <- feature_mapping |>
    dplyr::filter(passed_protein_coding_autosomal_and_X_filter) |>
    dplyr::left_join(
      annotation_analysis |>
        dplyr::select(
          ensembl_gene_id,
          chromosome_analysis = chromosome,
          gene_analysis = gene,
          gene_biotype_analysis = gene_biotype
        ),
      by = "ensembl_gene_id"
    ) |>
    dplyr::mutate(
      chromosome = chromosome_analysis,
      gene = gene_analysis,
      gene_biotype = gene_biotype_analysis
    ) |>
    dplyr::select(
      feature_id,
      chromosome,
      ensembl_gene_id,
      gene,
      gene_biotype
    )

  if (nrow(retained_mapping) == 0L) {
    stop(
      "No Seurat features matched Ensembl 115 protein-coding genes on ",
      "chromosomes 1-19 or X.",
      call. = FALSE
    )
  }

  if (anyDuplicated(retained_mapping$ensembl_gene_id)) {
    duplicated_ensembl_ids <- unique(
      retained_mapping$ensembl_gene_id[
        duplicated(retained_mapping$ensembl_gene_id)
      ]
    )

    stop(
      "Multiple Seurat features mapped to the same Ensembl gene ID: ",
      paste(head(duplicated_ensembl_ids, 20L), collapse = ", "),
      call. = FALSE
    )
  }

  retained_counts <- count_matrix[
    retained_mapping$feature_id,
    ,
    drop = FALSE
  ]

  rownames(retained_counts) <- retained_mapping$ensembl_gene_id

  retained_annotation <- retained_mapping |>
    dplyr::select(
      chromosome,
      ensembl_gene_id,
      gene,
      gene_biotype
    )

  if (!identical(
    rownames(retained_counts),
    retained_annotation$ensembl_gene_id
  )) {
    stop("Retained count matrix and annotation order do not match.", call. = FALSE)
  }

  list(
    counts = retained_counts,
    annotation = retained_annotation,
    mapping_status = feature_mapping,
    feature_identifier_type = feature_identifier_type
  )
}

aggregate_counts_per_sample_and_cluster <- function(
    count_matrix,
    seurat_metadata,
    sample_column,
    cluster_column,
    retained_sample_ids
) {
  if (!sample_column %in% colnames(seurat_metadata)) {
    stop(
      "Sample column does not exist in Seurat metadata: ",
      sample_column,
      call. = FALSE
    )
  }

  if (!cluster_column %in% colnames(seurat_metadata)) {
    stop(
      "Cluster column does not exist in Seurat metadata: ",
      cluster_column,
      call. = FALSE
    )
  }

  if (is.null(colnames(count_matrix))) {
    stop("Raw count matrix has no spot names.", call. = FALSE)
  }

  missing_metadata_spots <- setdiff(
    colnames(count_matrix),
    rownames(seurat_metadata)
  )

  if (length(missing_metadata_spots) > 0L) {
    stop(
      length(missing_metadata_spots),
      " count-matrix spots are absent from Seurat metadata. Examples: ",
      paste(head(missing_metadata_spots, 20L), collapse = ", "),
      call. = FALSE
    )
  }

  spot_metadata <- seurat_metadata[
    colnames(count_matrix),
    ,
    drop = FALSE
  ]

  spot_metadata$sample_ID_analysis <- trimws(
    as.character(spot_metadata[[sample_column]])
  )

  spot_metadata$cluster_id_analysis <- trimws(
    as.character(spot_metadata[[cluster_column]])
  )

  valid_spots <-
    !is.na(spot_metadata$sample_ID_analysis) &
    spot_metadata$sample_ID_analysis != "" &
    !is.na(spot_metadata$cluster_id_analysis) &
    spot_metadata$cluster_id_analysis != "" &
    spot_metadata$sample_ID_analysis %in% retained_sample_ids

  if (!any(valid_spots)) {
    stop(
      "No spots remain after selecting retained samples and non-missing clusters.",
      call. = FALSE
    )
  }

  count_matrix <- count_matrix[, valid_spots, drop = FALSE]
  spot_metadata <- spot_metadata[valid_spots, , drop = FALSE]

  group_separator <- "|||"
  pseudobulk_group_id <- paste(
    spot_metadata$sample_ID_analysis,
    spot_metadata$cluster_id_analysis,
    sep = group_separator
  )

  group_levels <- unique(pseudobulk_group_id)
  group_factor <- factor(pseudobulk_group_id, levels = group_levels)

  aggregation_matrix <- Matrix::sparseMatrix(
    i = seq_along(group_factor),
    j = as.integer(group_factor),
    x = 1,
    dims = c(length(group_factor), length(group_levels)),
    dimnames = list(colnames(count_matrix), group_levels)
  )

  pseudobulk_counts <- count_matrix %*% aggregation_matrix

  split_group_levels <- strsplit(
    group_levels,
    split = group_separator,
    fixed = TRUE
  )

  group_metadata <- data.frame(
    pseudobulk_id = group_levels,
    sample_ID = vapply(split_group_levels, `[[`, character(1), 1L),
    cluster_id = vapply(split_group_levels, `[[`, character(1), 2L),
    number_of_spots = as.integer(tabulate(as.integer(group_factor))),
    stringsAsFactors = FALSE
  )

  group_metadata$total_raw_counts <- as.numeric(Matrix::colSums(pseudobulk_counts))
  group_metadata$detected_genes <- as.integer(
    Matrix::colSums(pseudobulk_counts > 0)
  )

  if (!identical(colnames(pseudobulk_counts), group_metadata$pseudobulk_id)) {
    stop("Pseudobulk count columns and group metadata do not match.", call. = FALSE)
  }

  list(
    counts = pseudobulk_counts,
    group_metadata = group_metadata,
    spot_metadata = spot_metadata
  )
}


# ==============================================================================
# Descriptive expression and detection metrics
# ==============================================================================

compute_spot_detection_percentages_per_sample <- function(
    spot_count_matrix,
    spot_metadata,
    included_metadata,
    cluster_id,
    tested_gene_ids,
    positive_spot_min_count = 1L
) {
  positive_spot_min_count <- as.integer(positive_spot_min_count)

  if (
    length(positive_spot_min_count) != 1L ||
      is.na(positive_spot_min_count) ||
      positive_spot_min_count < 1L
  ) {
    stop(
      "`positive_spot_min_count` must be one positive integer.",
      call. = FALSE
    )
  }

  required_spot_metadata_columns <- c(
    "sample_ID_analysis",
    "cluster_id_analysis"
  )

  missing_spot_metadata_columns <- setdiff(
    required_spot_metadata_columns,
    colnames(spot_metadata)
  )

  if (length(missing_spot_metadata_columns) > 0L) {
    stop(
      "Missing spot-metadata columns required for detection percentages: ",
      paste(missing_spot_metadata_columns, collapse = ", "),
      call. = FALSE
    )
  }

  missing_tested_genes <- setdiff(
    tested_gene_ids,
    rownames(spot_count_matrix)
  )

  if (length(missing_tested_genes) > 0L) {
    stop(
      "Tested genes are absent from the spot-level count matrix. Examples: ",
      paste(head(missing_tested_genes, 20L), collapse = ", "),
      call. = FALSE
    )
  }

  included_sample_ids <- as.character(included_metadata$sample_ID)

  selected_spot_rows <-
    as.character(spot_metadata$cluster_id_analysis) == as.character(cluster_id) &
    as.character(spot_metadata$sample_ID_analysis) %in% included_sample_ids

  selected_spot_metadata <- spot_metadata[
    selected_spot_rows,
    ,
    drop = FALSE
  ]

  if (nrow(selected_spot_metadata) == 0L) {
    stop(
      "No spot-level observations were found for cluster ",
      cluster_id,
      " and the included samples.",
      call. = FALSE
    )
  }

  selected_spot_ids <- rownames(selected_spot_metadata)

  missing_selected_spots <- setdiff(
    selected_spot_ids,
    colnames(spot_count_matrix)
  )

  if (length(missing_selected_spots) > 0L) {
    stop(
      "Selected spots are absent from the spot-level count matrix. Examples: ",
      paste(head(missing_selected_spots, 20L), collapse = ", "),
      call. = FALSE
    )
  }

  detection_percentage_matrix <- matrix(
    NA_real_,
    nrow = length(tested_gene_ids),
    ncol = length(included_sample_ids),
    dimnames = list(
      tested_gene_ids,
      included_sample_ids
    )
  )

  for (sample_id in included_sample_ids) {
    sample_spot_ids <- rownames(selected_spot_metadata)[
      as.character(selected_spot_metadata$sample_ID_analysis) == sample_id
    ]

    if (length(sample_spot_ids) == 0L) {
      stop(
        "No spots were found for included sample ",
        sample_id,
        " in cluster ",
        cluster_id,
        ".",
        call. = FALSE
      )
    }

    sample_gene_counts <- spot_count_matrix[
      tested_gene_ids,
      sample_spot_ids,
      drop = FALSE
    ]

    positive_spot_counts <- Matrix::rowSums(
      sample_gene_counts >= positive_spot_min_count
    )

    detection_percentage_matrix[, sample_id] <-
      100 * as.numeric(positive_spot_counts) / length(sample_spot_ids)
  }

  detection_percentage_matrix
}



compute_four_cell_gene_abundance_qc <- function(
    dge,
    model_metadata,
    detection_percentage_matrix,
    sample_spot_counts,
    gene_annotation,
    cluster_id
) {
  sample_ids <- as.character(model_metadata$sample_ID)
  tested_gene_ids <- rownames(dge$counts)

  sample_spot_counts <- as.numeric(sample_spot_counts[sample_ids])
  names(sample_spot_counts) <- sample_ids

  if (
    anyNA(sample_spot_counts) ||
      any(!is.finite(sample_spot_counts)) ||
      any(sample_spot_counts <= 0)
  ) {
    stop(
      "Missing or invalid per-sample spot counts while computing descriptive metrics.",
      call. = FALSE
    )
  }

  if (!identical(colnames(dge$counts), sample_ids)) {
    stop(
      "DGEList columns and model metadata order do not match while computing descriptive metrics.",
      call. = FALSE
    )
  }

  if (!identical(colnames(detection_percentage_matrix), sample_ids)) {
    stop(
      "Detection-percentage columns and model metadata order do not match.",
      call. = FALSE
    )
  }

  if (!identical(rownames(detection_percentage_matrix), tested_gene_ids)) {
    stop(
      "Detection-percentage genes and tested DGEList genes do not match.",
      call. = FALSE
    )
  }

  tmm_cpm_matrix <- edgeR::cpm(
    dge,
    normalized.lib.sizes = TRUE,
    log = FALSE
  )

  if (!identical(dimnames(tmm_cpm_matrix), dimnames(dge$counts))) {
    stop(
      "TMM-normalized CPM matrix does not match the DGEList dimensions.",
      call. = FALSE
    )
  }

  model_metadata <- model_metadata |>
    dplyr::mutate(
      sample_ID = as.character(sample_ID),
      sex = as.character(sex),
      fmt_donor_group = as.character(fmt_donor_group),
      factor_cell = paste(sex, fmt_donor_group, sep = "_")
    )

  factor_cell_levels <- c(
    "Male_Neurotypical",
    "Male_ASD",
    "Female_Neurotypical",
    "Female_ASD"
  )

  unexpected_factor_cells <- setdiff(
    unique(model_metadata$factor_cell),
    factor_cell_levels
  )

  if (length(unexpected_factor_cells) > 0L) {
    stop(
      "Unexpected sex-by-group factor cells: ",
      paste(unexpected_factor_cells, collapse = ", "),
      call. = FALSE
    )
  }

  safe_row_sd <- function(x) {
    if (ncol(x) <= 1L) {
      return(rep(NA_real_, nrow(x)))
    }
    apply(x, 1L, stats::sd)
  }

  cell_metric_tables <- lapply(
    factor_cell_levels,
    function(factor_cell) {
      cell_sample_ids <- model_metadata$sample_ID[
        model_metadata$factor_cell == factor_cell
      ]

      if (length(cell_sample_ids) == 0L) {
        stop(
          "No samples were found for factor cell: ",
          factor_cell,
          call. = FALSE
        )
      }

      cell_detection <- detection_percentage_matrix[
        ,
        cell_sample_ids,
        drop = FALSE
      ]

      cell_cpm <- tmm_cpm_matrix[
        ,
        cell_sample_ids,
        drop = FALSE
      ]

      cell_spot_counts <- sample_spot_counts[cell_sample_ids]

      cell_positive_spot_counts <- sweep(
        cell_detection,
        MARGIN = 2L,
        STATS = cell_spot_counts / 100,
        FUN = "*"
      )

      cell_positive_spot_counts <- round(cell_positive_spot_counts)

      tibble::tibble(
        ensembl_gene_id = tested_gene_ids,

        !!paste0("mean_percent_positive_spots_", factor_cell) :=
          rowMeans(cell_detection),

        !!paste0("sd_percent_positive_spots_", factor_cell) :=
          safe_row_sd(cell_detection),

        !!paste0("median_percent_positive_spots_", factor_cell) :=
          apply(cell_detection, 1L, stats::median),

        !!paste0("samples_with_positive_spots_", factor_cell) :=
          as.integer(rowSums(cell_detection > 0)),

        !!paste0("samples_total_", factor_cell) :=
          as.integer(length(cell_sample_ids)),

        !!paste0("total_positive_spots_", factor_cell) :=
          as.integer(rowSums(cell_positive_spot_counts)),

        !!paste0("mean_positive_spots_per_sample_", factor_cell) :=
          rowMeans(cell_positive_spot_counts),

        !!paste0("mean_TMM_CPM_", factor_cell) :=
          rowMeans(cell_cpm),

        !!paste0("sd_TMM_CPM_", factor_cell) :=
          safe_row_sd(cell_cpm),

        !!paste0("median_TMM_CPM_", factor_cell) :=
          apply(cell_cpm, 1L, stats::median)
      )
    }
  )

  full_qc_table <- Reduce(
    function(x, y) {
      dplyr::left_join(x, y, by = "ensembl_gene_id")
    },
    cell_metric_tables
  ) |>
    dplyr::left_join(
      gene_annotation,
      by = "ensembl_gene_id"
    ) |>
    dplyr::mutate(
      cluster_id = as.character(cluster_id),
      .before = 1
    ) |>
    dplyr::select(
      cluster_id,
      ensembl_gene_id,
      gene,
      chromosome,
      dplyr::everything()
    )

  list(
    tmm_cpm_matrix = tmm_cpm_matrix,
    full_qc_table = full_qc_table
  )
}



make_descriptive_result_tables <- function(
    tmm_cpm_matrix,
    detection_percentage_matrix,
    model_metadata,
    test_definitions
) {
  if (!identical(dimnames(tmm_cpm_matrix), dimnames(detection_percentage_matrix))) {
    stop(
      "TMM-CPM and detection-percentage matrices must have identical dimensions and dimnames.",
      call. = FALSE
    )
  }

  model_metadata <- model_metadata |>
    dplyr::mutate(
      sample_ID = as.character(sample_ID),
      fmt_donor_group = as.character(fmt_donor_group),
      sex = as.character(sex)
    )

  if (!identical(colnames(tmm_cpm_matrix), model_metadata$sample_ID)) {
    stop(
      "Sample order in descriptive matrices does not match model metadata.",
      call. = FALSE
    )
  }

  tested_gene_ids <- rownames(tmm_cpm_matrix)

  safe_row_sd <- function(x) {
    if (ncol(x) <= 1L) {
      return(rep(NA_real_, nrow(x)))
    }
    apply(x, 1L, stats::sd)
  }

  summarize_samples <- function(sample_ids, displayed_name) {
    sample_ids <- as.character(sample_ids)

    if (length(sample_ids) == 0L) {
      stop(
        "No samples were selected for descriptive level: ",
        displayed_name,
        call. = FALSE
      )
    }

    missing_samples <- setdiff(sample_ids, colnames(tmm_cpm_matrix))

    if (length(missing_samples) > 0L) {
      stop(
        "Missing samples in descriptive matrices: ",
        paste(missing_samples, collapse = ", "),
        call. = FALSE
      )
    }

    detection_subset <- detection_percentage_matrix[
      ,
      sample_ids,
      drop = FALSE
    ]

    cpm_subset <- tmm_cpm_matrix[
      ,
      sample_ids,
      drop = FALSE
    ]

    tibble::tibble(
      ensembl_gene_id = tested_gene_ids,

      !!paste0("mean_percent_positive_spots_", displayed_name) :=
        rowMeans(detection_subset),

      !!paste0("sd_percent_positive_spots_", displayed_name) :=
        safe_row_sd(detection_subset),

      !!paste0("mean_TMM_CPM_", displayed_name) :=
        rowMeans(cpm_subset),

      !!paste0("sd_TMM_CPM_", displayed_name) :=
        safe_row_sd(cpm_subset)
    )
  }

  select_samples <- function(
      sex_value = NULL,
      group_value = NULL
  ) {
    keep <- rep(TRUE, nrow(model_metadata))

    if (!is.null(sex_value)) {
      keep <- keep & model_metadata$sex == sex_value
    }

    if (!is.null(group_value)) {
      keep <- keep & model_metadata$fmt_donor_group == group_value
    }

    model_metadata$sample_ID[keep]
  }

  join_metric_tables <- function(...) {
    metric_tables <- list(...)

    Reduce(
      function(x, y) {
        dplyr::left_join(x, y, by = "ensembl_gene_id")
      },
      metric_tables
    )
  }

  descriptive_tables <- list(
    Interaction = join_metric_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Male_Neurotypical"
      ),
      summarize_samples(
        select_samples("Male", "ASD"),
        "Male_ASD"
      ),
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Female_Neurotypical"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "Female_ASD"
      )
    ),

    Overall_Sex_Female_vs_Male = join_metric_tables(
      summarize_samples(
        select_samples(sex_value = "Male"),
        "Male"
      ),
      summarize_samples(
        select_samples(sex_value = "Female"),
        "Female"
      )
    ),

    Overall_Group_ASD_vs_Neurotypical = join_metric_tables(
      summarize_samples(
        select_samples(group_value = "Neurotypical"),
        "Neurotypical"
      ),
      summarize_samples(
        select_samples(group_value = "ASD"),
        "ASD"
      )
    ),

    ASD_Male_vs_Neurotypical_Male = join_metric_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Neurotypical_Male"
      ),
      summarize_samples(
        select_samples("Male", "ASD"),
        "ASD_Male"
      )
    ),

    ASD_Female_vs_Neurotypical_Female = join_metric_tables(
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Neurotypical_Female"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "ASD_Female"
      )
    ),

    Neurotypical_Female_vs_Neurotypical_Male = join_metric_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Neurotypical_Male"
      ),
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Neurotypical_Female"
      )
    ),

    ASD_Female_vs_ASD_Male = join_metric_tables(
      summarize_samples(
        select_samples("Male", "ASD"),
        "ASD_Male"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "ASD_Female"
      )
    ),

    ASD_Female_vs_Neurotypical_Male = join_metric_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Neurotypical_Male"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "ASD_Female"
      )
    ),

    ASD_Male_vs_Neurotypical_Female = join_metric_tables(
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Neurotypical_Female"
      ),
      summarize_samples(
        select_samples("Male", "ASD"),
        "ASD_Male"
      )
    )
  )

  expected_test_ids <- as.character(test_definitions$test_id)

  if (!setequal(names(descriptive_tables), expected_test_ids)) {
    stop(
      "Descriptive table names do not match the edgeR test definitions.",
      call. = FALSE
    )
  }

  descriptive_tables[expected_test_ids]
}


create_two_factor_design_and_contrasts <- function(cluster_sample_metadata) {
  model_metadata <- cluster_sample_metadata |>
    dplyr::transmute(
      sample_ID = as.character(sample_ID),
      fmt_donor_group = factor(
        as.character(fmt_donor_group),
        levels = c("Neurotypical", "ASD")
      ),
      sex = factor(
        as.character(sex),
        levels = c("Male", "Female")
      )
    )

  if (anyNA(model_metadata$fmt_donor_group) || anyNA(model_metadata$sex)) {
    stop("Missing or unexpected sex/group values in cluster metadata.", call. = FALSE)
  }

  design <- model.matrix(
    ~ fmt_donor_group * sex,
    data = model_metadata
  )

  rownames(design) <- model_metadata$sample_ID

  expected_original_columns <- c(
    "(Intercept)",
    "fmt_donor_groupASD",
    "sexFemale",
    "fmt_donor_groupASD:sexFemale"
  )

  if (!identical(colnames(design), expected_original_columns)) {
    stop(
      "Unexpected design columns: ",
      paste(colnames(design), collapse = ", "),
      call. = FALSE
    )
  }

  colnames(design) <- c(
    "Intercept",
    "groupASD",
    "sexFemale",
    "groupASD_sexFemale"
  )

  if (qr(design)$rank != ncol(design)) {
    stop("The cluster-specific design matrix is not full rank.", call. = FALSE)
  }

  contrasts <- limma::makeContrasts(
    Interaction = groupASD_sexFemale,

    Overall_Sex_Female_vs_Male =
      sexFemale + 0.5 * groupASD_sexFemale,

    Overall_Group_ASD_vs_Neurotypical =
      groupASD + 0.5 * groupASD_sexFemale,

    ASD_Male_vs_Neurotypical_Male =
      groupASD,

    ASD_Female_vs_Neurotypical_Female =
      groupASD + groupASD_sexFemale,

    Neurotypical_Female_vs_Neurotypical_Male =
      sexFemale,

    ASD_Female_vs_ASD_Male =
      sexFemale + groupASD_sexFemale,

    ASD_Female_vs_Neurotypical_Male =
      groupASD + sexFemale + groupASD_sexFemale,

    ASD_Male_vs_Neurotypical_Female =
      groupASD - sexFemale,

    levels = design
  )

  test_definitions <- tibble::tribble(
    ~test_id, ~sheet_name, ~comparison, ~positive_logFC_meaning, ~negative_logFC_meaning, ~add_regulation,

    "Interaction",
    "Interaction",
    "(ASD - Neurotypical)_Female - (ASD - Neurotypical)_Male",
    "ASD effect is more positive in Female than Male",
    "ASD effect is more negative in Female than Male",
    FALSE,

    "Overall_Sex_Female_vs_Male",
    "Overall_Sex",
    "Female vs Male averaged equally across donor groups",
    "higher expression in Female",
    "higher expression in Male",
    TRUE,

    "Overall_Group_ASD_vs_Neurotypical",
    "Overall_Group",
    "ASD vs Neurotypical averaged equally across sexes",
    "higher expression in ASD",
    "higher expression in Neurotypical",
    TRUE,

    "ASD_Male_vs_Neurotypical_Male",
    "Group_in_Male",
    "ASD Male vs Neurotypical Male",
    "higher expression in ASD Male",
    "higher expression in Neurotypical Male",
    TRUE,

    "ASD_Female_vs_Neurotypical_Female",
    "Group_in_Female",
    "ASD Female vs Neurotypical Female",
    "higher expression in ASD Female",
    "higher expression in Neurotypical Female",
    TRUE,

    "Neurotypical_Female_vs_Neurotypical_Male",
    "Sex_in_Neurotypical",
    "Neurotypical Female vs Neurotypical Male",
    "higher expression in Neurotypical Female",
    "higher expression in Neurotypical Male",
    TRUE,

    "ASD_Female_vs_ASD_Male",
    "Sex_in_ASD",
    "ASD Female vs ASD Male",
    "higher expression in ASD Female",
    "higher expression in ASD Male",
    TRUE,

    "ASD_Female_vs_Neurotypical_Male",
    "ASD_F_vs_NT_M",
    "ASD Female vs Neurotypical Male",
    "higher expression in ASD Female",
    "higher expression in Neurotypical Male",
    TRUE,

    "ASD_Male_vs_Neurotypical_Female",
    "ASD_M_vs_NT_F",
    "ASD Male vs Neurotypical Female",
    "higher expression in ASD Male",
    "higher expression in Neurotypical Female",
    TRUE
  )

  if (!identical(colnames(contrasts), test_definitions$test_id)) {
    stop("Contrast order does not match test definitions.", call. = FALSE)
  }

  list(
    metadata = model_metadata,
    design = design,
    contrasts = contrasts,
    test_definitions = test_definitions
  )
}

create_edger_result_table_per_cluster <- function(
    qlf_test,
    gene_annotation,
    add_regulation,
    descriptive_metrics_table
) {
  if (is.null(descriptive_metrics_table)) {
    stop(
      "`descriptive_metrics_table` cannot be NULL.",
      call. = FALSE
    )
  }

  if (!"ensembl_gene_id" %in% colnames(descriptive_metrics_table)) {
    stop(
      "The descriptive metric table must contain `ensembl_gene_id`.",
      call. = FALSE
    )
  }

  result_table <- edgeR::topTags(
    qlf_test,
    n = Inf,
    sort.by = "PValue"
  )$table |>
    as.data.frame(stringsAsFactors = FALSE) |>
    tibble::rownames_to_column(var = "ensembl_gene_id") |>
    dplyr::left_join(
      gene_annotation,
      by = "ensembl_gene_id"
    ) |>
    dplyr::mutate(
      gene = dplyr::if_else(
        is.na(gene) | gene == "",
        ensembl_gene_id,
        gene
      )
    )

  if (add_regulation) {
    result_table <- result_table |>
      dplyr::mutate(
        regulation = dplyr::if_else(logFC >= 0, "up", "down")
      ) |>
      dplyr::select(
        ensembl_gene_id,
        gene,
        chromosome,
        logFC,
        logCPM,
        F,
        PValue,
        FDR,
        regulation
      )
  } else {
    result_table <- result_table |>
      dplyr::select(
        ensembl_gene_id,
        gene,
        chromosome,
        logFC,
        logCPM,
        F,
        PValue,
        FDR
      )
  }

  result_table <- result_table |>
    dplyr::left_join(
      descriptive_metrics_table,
      by = "ensembl_gene_id"
    )

  descriptive_columns <- setdiff(
    colnames(descriptive_metrics_table),
    "ensembl_gene_id"
  )

  if (anyNA(result_table[, descriptive_columns, drop = FALSE])) {
    stop(
      "Missing descriptive metrics after joining them to an edgeR result table.",
      call. = FALSE
    )
  }

  result_table |>
    dplyr::arrange(PValue, dplyr::desc(abs(logFC)))
}


write_result_workbook_per_cluster <- function(
    result_list,
    test_definitions,
    output_file,
    include_cluster_id = FALSE,
    xlsx_max_decimal_places = 5L
) {
  xlsx_max_decimal_places <- as.integer(xlsx_max_decimal_places)

  if (
    length(xlsx_max_decimal_places) != 1L ||
      is.na(xlsx_max_decimal_places) ||
      xlsx_max_decimal_places < 0L ||
      xlsx_max_decimal_places > 15L
  ) {
    stop(
      "`xlsx_max_decimal_places` must be one integer between 0 and 15.",
      call. = FALSE
    )
  }

  workbook <- openxlsx::createWorkbook()

  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#D9EAF7",
    border = "Bottom"
  )

  decimal_pattern <- if (xlsx_max_decimal_places == 0L) {
    "#,##0"
  } else {
    paste0(
      "#,##0.",
      paste(rep("#", xlsx_max_decimal_places), collapse = "")
    )
  }

  scientific_pattern <- if (xlsx_max_decimal_places == 0L) {
    "0E+00"
  } else {
    paste0(
      "0.",
      paste(rep("#", xlsx_max_decimal_places), collapse = ""),
      "E+00"
    )
  }

  decimal_style <- openxlsx::createStyle(numFmt = decimal_pattern)
  scientific_style <- openxlsx::createStyle(numFmt = scientific_pattern)
  integer_style <- openxlsx::createStyle(numFmt = "#,##0")

  for (test_index in seq_len(nrow(test_definitions))) {
    test_id <- test_definitions$test_id[[test_index]]
    sheet_name <- test_definitions$sheet_name[[test_index]]
    result_table <- result_list[[test_id]]

    if (is.null(result_table)) {
      stop("Missing result table for workbook test: ", test_id, call. = FALSE)
    }

    if (nrow(result_table) > 1048575L) {
      stop(
        "Result table exceeds the Excel row limit for sheet: ",
        sheet_name,
        call. = FALSE
      )
    }

    expected_first_columns <- if (include_cluster_id) {
      c("cluster_id", "ensembl_gene_id", "gene", "chromosome")
    } else {
      c("ensembl_gene_id", "gene", "chromosome")
    }

    if (!identical(
      colnames(result_table)[seq_along(expected_first_columns)],
      expected_first_columns
    )) {
      stop(
        "Unexpected leading-column order in sheet ",
        sheet_name,
        ". Expected: ",
        paste(expected_first_columns, collapse = ", "),
        "; observed: ",
        paste(
          colnames(result_table)[seq_along(expected_first_columns)],
          collapse = ", "
        ),
        call. = FALSE
      )
    }

    descriptive_columns <- grep(
      "^(mean|sd)_(percent_positive_spots|TMM_CPM)_",
      colnames(result_table),
      value = TRUE
    )

    expected_descriptive_count <- if (identical(test_id, "Interaction")) {
      16L
    } else {
      8L
    }

    if (length(descriptive_columns) != expected_descriptive_count) {
      stop(
        "Sheet ", sheet_name,
        " has ", length(descriptive_columns),
        " descriptive columns, but expected ",
        expected_descriptive_count,
        ". The XLSX will not be written with incomplete descriptive metrics.",
        call. = FALSE
      )
    }

    openxlsx::addWorksheet(workbook, sheetName = sheet_name)

    openxlsx::writeData(
      workbook,
      sheet = sheet_name,
      x = result_table,
      startRow = 1,
      startCol = 1,
      headerStyle = header_style,
      withFilter = TRUE
    )

    openxlsx::freezePane(
      workbook,
      sheet = sheet_name,
      firstActiveRow = 2,
      firstActiveCol = 3
    )

    numeric_columns <- which(
      vapply(result_table, is.numeric, FUN.VALUE = logical(1))
    )

    scientific_columns <- which(
      colnames(result_table) %in% c("PValue", "FDR")
    )

    integer_columns <- which(
      grepl(
        "^(cluster_id|samples_total_|samples_with_positive_spots_|total_positive_spots_)",
        colnames(result_table)
      )
    )

    regular_numeric_columns <- setdiff(
      numeric_columns,
      union(scientific_columns, integer_columns)
    )

    if (length(regular_numeric_columns) > 0L && nrow(result_table) > 0L) {
      openxlsx::addStyle(
        workbook,
        sheet = sheet_name,
        style = decimal_style,
        rows = 2:(nrow(result_table) + 1L),
        cols = regular_numeric_columns,
        gridExpand = TRUE,
        stack = TRUE
      )
    }

    if (length(scientific_columns) > 0L && nrow(result_table) > 0L) {
      openxlsx::addStyle(
        workbook,
        sheet = sheet_name,
        style = scientific_style,
        rows = 2:(nrow(result_table) + 1L),
        cols = scientific_columns,
        gridExpand = TRUE,
        stack = TRUE
      )
    }

    if (length(integer_columns) > 0L && nrow(result_table) > 0L) {
      openxlsx::addStyle(
        workbook,
        sheet = sheet_name,
        style = integer_style,
        rows = 2:(nrow(result_table) + 1L),
        cols = integer_columns,
        gridExpand = TRUE,
        stack = TRUE
      )
    }

    column_widths <- rep(13, ncol(result_table))
    names(column_widths) <- colnames(result_table)

    if ("cluster_id" %in% names(column_widths)) {
      column_widths[["cluster_id"]] <- 12
    }
    if ("ensembl_gene_id" %in% names(column_widths)) {
      column_widths[["ensembl_gene_id"]] <- 22
    }
    if ("gene" %in% names(column_widths)) {
      column_widths[["gene"]] <- 18
    }
    if ("chromosome" %in% names(column_widths)) {
      column_widths[["chromosome"]] <- 11
    }
    if ("regulation" %in% names(column_widths)) {
      column_widths[["regulation"]] <- 12
    }

    descriptive_metric_columns <- grepl(
      "^(mean|sd)_(percent_positive_spots|TMM_CPM)_",
      names(column_widths)
    )

    column_widths[descriptive_metric_columns] <- 29

    openxlsx::setColWidths(
      workbook,
      sheet = sheet_name,
      cols = seq_len(ncol(result_table)),
      widths = unname(column_widths)
    )
  }

  openxlsx::saveWorkbook(
    workbook,
    file = output_file,
    overwrite = TRUE
  )

  invisible(output_file)
}


run_one_cluster_edger_two_factor_interaction <- function(
    cluster_id,
    pseudobulk_counts_all,
    pseudobulk_group_metadata,
    spot_count_matrix_all,
    spot_metadata_all,
    sample_metadata,
    gene_annotation,
    min_spots_per_sample_cluster,
    min_samples_per_factor_cell,
    positive_spot_min_count,
    filter_by_expr_parameters
) {
  cluster_group_metadata <- pseudobulk_group_metadata |>
    dplyr::filter(.data$cluster_id == .env$cluster_id) |>
    dplyr::left_join(
      sample_metadata |>
        dplyr::select(sample_ID, fmt_donor_group, sex),
      by = "sample_ID"
    ) |>
    dplyr::mutate(
      included_by_spot_threshold =
        number_of_spots >= min_spots_per_sample_cluster,
      inclusion_status = dplyr::if_else(
        included_by_spot_threshold,
        "included",
        "excluded_below_minimum_spots"
      )
    ) |>
    dplyr::arrange(match(sample_ID, sample_metadata$sample_ID))

  included_metadata <- cluster_group_metadata |>
    dplyr::filter(included_by_spot_threshold)

  if (nrow(included_metadata) == 0L) {
    return(
      list(
        status = "skipped",
        reason = "no_sample_passed_minimum_spot_threshold",
        cluster_group_metadata = cluster_group_metadata
      )
    )
  }

  if (anyNA(included_metadata$fmt_donor_group) || anyNA(included_metadata$sex)) {
    stop(
      "Missing sample metadata after joining cluster pseudobulk metadata.",
      call. = FALSE
    )
  }

  factor_cell_counts <- included_metadata |>
    dplyr::count(sex, fmt_donor_group, name = "n_samples")

  complete_factor_grid <- expand.grid(
    sex = c("Male", "Female"),
    fmt_donor_group = c("Neurotypical", "ASD"),
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(
      factor_cell_counts,
      by = c("sex", "fmt_donor_group")
    ) |>
    dplyr::mutate(n_samples = dplyr::coalesce(n_samples, 0L))

  if (any(complete_factor_grid$n_samples < min_samples_per_factor_cell)) {
    return(
      list(
        status = "skipped",
        reason = paste0(
          "at_least_one_factor_cell_has_fewer_than_",
          min_samples_per_factor_cell,
          "_samples"
        ),
        cluster_group_metadata = cluster_group_metadata,
        factor_cell_counts = complete_factor_grid
      )
    )
  }

  selected_pseudobulk_ids <- included_metadata$pseudobulk_id

  cluster_counts <- pseudobulk_counts_all[
    ,
    selected_pseudobulk_ids,
    drop = FALSE
  ]

  colnames(cluster_counts) <- included_metadata$sample_ID

  if (any(Matrix::colSums(cluster_counts) <= 0)) {
    return(
      list(
        status = "skipped",
        reason = "at_least_one_included_sample_has_zero_library_size",
        cluster_group_metadata = cluster_group_metadata,
        factor_cell_counts = complete_factor_grid
      )
    )
  }

  design_objects <- create_two_factor_design_and_contrasts(
    included_metadata
  )

  if (!identical(
    colnames(cluster_counts),
    design_objects$metadata$sample_ID
  )) {
    stop("Count columns and design metadata order do not match.", call. = FALSE)
  }

  dge_unfiltered <- edgeR::DGEList(
    counts = cluster_counts,
    samples = data.frame(
      sample_ID = design_objects$metadata$sample_ID,
      fmt_donor_group = as.character(
        design_objects$metadata$fmt_donor_group
      ),
      sex = as.character(design_objects$metadata$sex),
      row.names = design_objects$metadata$sample_ID,
      stringsAsFactors = FALSE
    )
  )

  keep_genes <- edgeR::filterByExpr(
    dge_unfiltered,
    design = design_objects$design,
    min.count = filter_by_expr_parameters$min.count,
    min.total.count = filter_by_expr_parameters$min.total.count,
    large.n = filter_by_expr_parameters$large.n,
    min.prop = filter_by_expr_parameters$min.prop
  )

  if (!any(keep_genes)) {
    return(
      list(
        status = "skipped",
        reason = "no_gene_passed_filterByExpr",
        cluster_group_metadata = cluster_group_metadata,
        factor_cell_counts = complete_factor_grid,
        keep_genes = keep_genes
      )
    )
  }

  filter_status <- tibble::tibble(
    ensembl_gene_id = rownames(dge_unfiltered$counts),
    total_count = rowSums(dge_unfiltered$counts),
    mean_count = rowMeans(dge_unfiltered$counts),
    samples_with_nonzero_count = rowSums(dge_unfiltered$counts > 0),
    passed_filterByExpr = keep_genes,
    filterByExpr_status = dplyr::if_else(
      keep_genes,
      "passed_filterByExpr",
      "filtered_out_by_filterByExpr"
    )
  ) |>
    dplyr::left_join(gene_annotation, by = "ensembl_gene_id") |>
    dplyr::select(
      ensembl_gene_id,
      gene,
      chromosome,
      total_count,
      mean_count,
      samples_with_nonzero_count,
      passed_filterByExpr,
      filterByExpr_status
    )

  dge <- dge_unfiltered[
    keep_genes,
    ,
    keep.lib.sizes = FALSE
  ]

  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  detection_percentage_matrix <-
    compute_spot_detection_percentages_per_sample(
      spot_count_matrix = spot_count_matrix_all,
      spot_metadata = spot_metadata_all,
      included_metadata = included_metadata,
      cluster_id = cluster_id,
      tested_gene_ids = rownames(dge$counts),
      positive_spot_min_count = positive_spot_min_count
    )

  gene_abundance_objects <- compute_four_cell_gene_abundance_qc(
    dge = dge,
    model_metadata = design_objects$metadata,
    detection_percentage_matrix = detection_percentage_matrix,
    sample_spot_counts = stats::setNames(
      included_metadata$number_of_spots,
      included_metadata$sample_ID
    ),
    gene_annotation = gene_annotation,
    cluster_id = cluster_id
  )

  descriptive_result_tables <- make_descriptive_result_tables(
    tmm_cpm_matrix = gene_abundance_objects$tmm_cpm_matrix,
    detection_percentage_matrix = detection_percentage_matrix,
    model_metadata = design_objects$metadata,
    test_definitions = design_objects$test_definitions
  )

  dge <- edgeR::estimateDisp(
    dge,
    design = design_objects$design,
    robust = TRUE
  )

  fit <- edgeR::glmQLFit(
    dge,
    design = design_objects$design,
    robust = TRUE
  )

  test_objects <- lapply(
    design_objects$test_definitions$test_id,
    function(test_id) {
      edgeR::glmQLFTest(
        fit,
        contrast = design_objects$contrasts[, test_id]
      )
    }
  )

  names(test_objects) <- design_objects$test_definitions$test_id

  result_tables <- lapply(
    seq_len(nrow(design_objects$test_definitions)),
    function(test_index) {
      test_id <- design_objects$test_definitions$test_id[[test_index]]
      add_regulation <- design_objects$test_definitions$add_regulation[[test_index]]

      create_edger_result_table_per_cluster(
        qlf_test = test_objects[[test_id]],
        gene_annotation = gene_annotation,
        add_regulation = add_regulation,
        descriptive_metrics_table =
          descriptive_result_tables[[test_id]]
      )
    }
  )

  names(result_tables) <- design_objects$test_definitions$test_id

  analysis_summary <- dplyr::bind_rows(
    lapply(
      design_objects$test_definitions$test_id,
      function(test_id) {
        result_table <- result_tables[[test_id]]
        test_definition <- design_objects$test_definitions |>
          dplyr::filter(.data$test_id == .env$test_id)

        tibble::tibble(
          cluster_id = as.character(cluster_id),
          test_id = test_id,
          sheet_name = test_definition$sheet_name,
          comparison = test_definition$comparison,
          tested_genes = nrow(result_table),
          significant_P_0.05 = sum(
            result_table$PValue < 0.05,
            na.rm = TRUE
          ),
          significant_FDR_0.05 = sum(
            result_table$FDR < 0.05,
            na.rm = TRUE
          ),
          significant_FDR_0.05_abs_logFC_ge_0.5 = sum(
            result_table$FDR < 0.05 &
              abs(result_table$logFC) >= 0.5,
            na.rm = TRUE
          ),
          significant_FDR_0.10 = sum(
            result_table$FDR < 0.10,
            na.rm = TRUE
          ),
          significant_FDR_0.10_abs_logFC_ge_0.5 = sum(
            result_table$FDR < 0.10 &
              abs(result_table$logFC) >= 0.5,
            na.rm = TRUE
          )
        )
      }
    )
  )

  design_output <- data.frame(
    sample_ID = rownames(design_objects$design),
    design_objects$design,
    row.names = NULL,
    check.names = FALSE
  )

  contrast_output <- data.frame(
    coefficient = rownames(design_objects$contrasts),
    design_objects$contrasts,
    row.names = NULL,
    check.names = FALSE
  )

  library_summary <- tibble::tibble(
    cluster_id = as.character(cluster_id),
    sample_ID = colnames(dge$counts),
    fmt_donor_group = as.character(design_objects$metadata$fmt_donor_group),
    sex = as.character(design_objects$metadata$sex),
    number_of_spots = as.integer(
      included_metadata$number_of_spots[
        match(colnames(dge$counts), included_metadata$sample_ID)
      ]
    ),
    raw_library_size = as.numeric(dge$samples$lib.size),
    normalization_factor = as.numeric(dge$samples$norm.factors),
    effective_library_size =
      as.numeric(dge$samples$lib.size * dge$samples$norm.factors)
  )

  list(
    status = "completed",
    reason = NA_character_,
    cluster_group_metadata = cluster_group_metadata,
    factor_cell_counts = complete_factor_grid,
    counts_before_filterByExpr = cluster_counts,
    counts_tested_genes = dge$counts,
    keep_genes = keep_genes,
    filter_status = filter_status,
    design = design_objects$design,
    design_output = design_output,
    contrasts = design_objects$contrasts,
    contrast_output = contrast_output,
    test_definitions = design_objects$test_definitions,
    dge = dge,
    fit = fit,
    test_objects = test_objects,
    results = result_tables,
    descriptive_result_tables = descriptive_result_tables,
    gene_abundance_qc = gene_abundance_objects$full_qc_table,
    detection_percentage_matrix = detection_percentage_matrix,
    tmm_cpm_matrix = gene_abundance_objects$tmm_cpm_matrix,
    analysis_summary = analysis_summary,
    library_summary = library_summary
  )
}

run_pseudobulk_per_cluster_edger_two_factor_interaction <- function(
    project_root,
    dataset_name,
    clustering_name,
    input_rdata_file,
    output_statistics_dir,
    metadata_file,
    annotation_rds_file,
    annotation_metadata_file,
    seurat_object_name = NULL,
    assay_name = "RNA",
    sample_column = "sample_ID",
    cluster_column = NULL,
    excluded_samples = c("20_1F", "12_3F", "15_1M", "20_3M"),
    expected_number_of_samples = 16L,
    expected_number_of_clusters = 16L,
    min_spots_per_sample_cluster = 20L,
    min_samples_per_factor_cell = 2L,
    positive_spot_min_count = 1L,
    xlsx_max_decimal_places = 5L,
    filter_by_expr_parameters = list(
      min.count = 10,
      min.total.count = 15,
      large.n = 10,
      min.prop = 0.7
    ),
    stop_on_cluster_error = TRUE
) {
  check_required_packages_pseudobulk_per_cluster()

  if (!dir.exists(project_root)) {
    stop("Project root does not exist: ", project_root, call. = FALSE)
  }

  if (!file.exists(metadata_file)) {
    stop("Metadata file does not exist: ", metadata_file, call. = FALSE)
  }

  positive_spot_min_count <- as.integer(positive_spot_min_count)

  if (
    length(positive_spot_min_count) != 1L ||
      is.na(positive_spot_min_count) ||
      positive_spot_min_count < 1L
  ) {
    stop(
      "`positive_spot_min_count` must be one positive integer.",
      call. = FALSE
    )
  }

  xlsx_max_decimal_places <- as.integer(xlsx_max_decimal_places)

  if (
    length(xlsx_max_decimal_places) != 1L ||
      is.na(xlsx_max_decimal_places) ||
      xlsx_max_decimal_places < 0L ||
      xlsx_max_decimal_places > 15L
  ) {
    stop(
      "`xlsx_max_decimal_places` must be one integer between 0 and 15.",
      call. = FALSE
    )
  }

  required_filter_parameters <- c(
    "min.count",
    "min.total.count",
    "large.n",
    "min.prop"
  )

  missing_filter_parameters <- setdiff(
    required_filter_parameters,
    names(filter_by_expr_parameters)
  )

  if (length(missing_filter_parameters) > 0L) {
    stop(
      "Missing filterByExpr parameters: ",
      paste(missing_filter_parameters, collapse = ", "),
      call. = FALSE
    )
  }

  output_subdirectories <- list(
    pseudobulk_inputs = file.path(
      output_statistics_dir,
      "01_pseudobulkInputs"
    ),
    design_metadata = file.path(
      output_statistics_dir,
      "02_designAndSampleMetadata"
    ),
    edger_results = file.path(
      output_statistics_dir,
      "03_edgeRResults"
    ),
    diagnostics = file.path(
      output_statistics_dir,
      "04_modelDiagnostics"
    ),
    objects = file.path(
      output_statistics_dir,
      "05_objects"
    ),
    logs = file.path(
      output_statistics_dir,
      "06_logs"
    )
  )

  invisible(
    lapply(
      c(output_statistics_dir, unname(output_subdirectories)),
      dir.create,
      recursive = TRUE,
      showWarnings = FALSE
    )
  )

  per_cluster_xlsx_dir <- file.path(
    output_subdirectories$edger_results,
    "01_perCluster_XLSX"
  )

  combined_results_dir <- file.path(
    output_subdirectories$edger_results,
    "02_combinedResults"
  )

  per_cluster_counts_dir <- file.path(
    output_subdirectories$pseudobulk_inputs,
    "perCluster_rawCounts"
  )

  per_cluster_sample_metrics_dir <- file.path(
    output_subdirectories$pseudobulk_inputs,
    "perCluster_sampleLevelMetrics"
  )

  per_cluster_design_dir <- file.path(
    output_subdirectories$design_metadata,
    "perCluster"
  )

  per_cluster_diagnostics_dir <- file.path(
    output_subdirectories$diagnostics,
    "perCluster"
  )

  invisible(
    lapply(
      c(
        per_cluster_xlsx_dir,
        combined_results_dir,
        per_cluster_counts_dir,
        per_cluster_sample_metrics_dir,
        per_cluster_design_dir,
        per_cluster_diagnostics_dir
      ),
      dir.create,
      recursive = TRUE,
      showWarnings = FALSE
    )
  )

  analysis_prefix <- paste0(
    dataset_name,
    "_",
    clustering_name,
    "_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction"
  )

  message("\nLoading input RData: ", input_rdata_file)

  loaded_seurat <- load_single_seurat_object_from_rdata(
    input_rdata_file = input_rdata_file,
    requested_object_name = seurat_object_name
  )

  seurat_object <- loaded_seurat$object
  seurat_object_name_resolved <- loaded_seurat$object_name
  seurat_metadata <- seurat_object[[]]

  cluster_column_resolved <- resolve_cluster_column_pseudobulk_per_cluster(
    seurat_metadata = seurat_metadata,
    requested_cluster_column = cluster_column,
    expected_number_of_clusters = expected_number_of_clusters
  )

  message("Selected Seurat object: ", seurat_object_name_resolved)
  message("Selected assay: ", assay_name)
  message("Selected sample column: ", sample_column)
  message("Selected cluster column: ", cluster_column_resolved)

  metadata_all <- read.delim(
    metadata_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    tibble::as_tibble()

  required_metadata_columns <- c(
    "sample_ID",
    "fmt_donor_group",
    "sex"
  )

  missing_metadata_columns <- setdiff(
    required_metadata_columns,
    colnames(metadata_all)
  )

  if (length(missing_metadata_columns) > 0L) {
    stop(
      "Missing metadata columns: ",
      paste(missing_metadata_columns, collapse = ", "),
      call. = FALSE
    )
  }

  metadata_all <- metadata_all |>
    dplyr::mutate(
      sample_ID = trimws(as.character(sample_ID)),
      fmt_donor_group = trimws(as.character(fmt_donor_group)),
      sex = trimws(as.character(sex))
    )

  if (anyDuplicated(metadata_all$sample_ID)) {
    stop("Metadata contains duplicated sample_ID values.", call. = FALSE)
  }

  sample_metadata <- metadata_all |>
    dplyr::filter(!sample_ID %in% excluded_samples) |>
    dplyr::mutate(
      fmt_donor_group = factor(
        fmt_donor_group,
        levels = c("Neurotypical", "ASD")
      ),
      sex = factor(
        sex,
        levels = c("Male", "Female")
      )
    )

  if (nrow(sample_metadata) != expected_number_of_samples) {
    stop(
      "Expected ", expected_number_of_samples,
      " retained samples after QC exclusion, but found ",
      nrow(sample_metadata), ".",
      call. = FALSE
    )
  }

  if (anyNA(sample_metadata$fmt_donor_group) || anyNA(sample_metadata$sex)) {
    stop("Unexpected group or sex values in metadata.", call. = FALSE)
  }

  count_extraction <- extract_raw_counts_from_seurat_layers(
    seurat_object = seurat_object,
    assay_name = assay_name
  )

  raw_count_matrix <- count_extraction$counts

  annotation_objects <- read_and_validate_ensembl115_annotation(
    annotation_rds_file = annotation_rds_file,
    annotation_metadata_file = annotation_metadata_file
  )

  mapped_counts <- map_seurat_features_to_analysis_annotation(
    count_matrix = raw_count_matrix,
    annotation_all = annotation_objects$annotation_all,
    annotation_analysis = annotation_objects$annotation_analysis
  )

  rm(raw_count_matrix)
  invisible(gc())

  pseudobulk_objects <- aggregate_counts_per_sample_and_cluster(
    count_matrix = mapped_counts$counts,
    seurat_metadata = seurat_metadata,
    sample_column = sample_column,
    cluster_column = cluster_column_resolved,
    retained_sample_ids = sample_metadata$sample_ID
  )

  cluster_ids <- natural_order_cluster_ids_pseudobulk_per_cluster(
    pseudobulk_objects$group_metadata$cluster_id
  )

  if (length(cluster_ids) != expected_number_of_clusters) {
    stop(
      "Expected ", expected_number_of_clusters,
      " clusters, but found ", length(cluster_ids),
      ": ", paste(cluster_ids, collapse = ", "),
      call. = FALSE
    )
  }

  pseudobulk_input_summary <- pseudobulk_objects$group_metadata |>
    dplyr::left_join(
      sample_metadata |>
        dplyr::mutate(
          fmt_donor_group = as.character(fmt_donor_group),
          sex = as.character(sex)
        ) |>
        dplyr::select(sample_ID, fmt_donor_group, sex),
      by = "sample_ID"
    ) |>
    dplyr::mutate(
      min_spots_required = min_spots_per_sample_cluster,
      included_by_spot_threshold =
        number_of_spots >= min_spots_per_sample_cluster
    ) |>
    dplyr::arrange(
      match(cluster_id, cluster_ids),
      match(sample_ID, sample_metadata$sample_ID)
    )

  write_tsv_pseudobulk_per_cluster(
    pseudobulk_input_summary,
    file.path(
      output_subdirectories$pseudobulk_inputs,
      paste0(analysis_prefix, "_pseudobulkInputSummary_allClusters.tsv")
    )
  )

  write_tsv_pseudobulk_per_cluster(
    mapped_counts$mapping_status,
    file.path(
      output_subdirectories$pseudobulk_inputs,
      paste0(analysis_prefix, "_featureMappingStatus.tsv.gz")
    )
  )

  cluster_run_objects <- vector("list", length(cluster_ids))
  names(cluster_run_objects) <- cluster_ids

  cluster_status_rows <- vector("list", length(cluster_ids))
  names(cluster_status_rows) <- cluster_ids

  first_completed_test_definitions <- NULL

  for (cluster_id_current in cluster_ids) {
    message("\n", paste(rep("=", 80), collapse = ""))
    message("Processing cluster: ", cluster_id_current)
    message(paste(rep("=", 80), collapse = ""))

    cluster_slug <- sanitize_file_component_pseudobulk_per_cluster(
      cluster_id_current
    )

    cluster_result <- tryCatch(
      run_one_cluster_edger_two_factor_interaction(
        cluster_id = cluster_id_current,
        pseudobulk_counts_all = pseudobulk_objects$counts,
        pseudobulk_group_metadata = pseudobulk_objects$group_metadata,
        spot_count_matrix_all = mapped_counts$counts,
        spot_metadata_all = pseudobulk_objects$spot_metadata,
        sample_metadata = sample_metadata,
        gene_annotation = mapped_counts$annotation,
        min_spots_per_sample_cluster = min_spots_per_sample_cluster,
        min_samples_per_factor_cell = min_samples_per_factor_cell,
        positive_spot_min_count = positive_spot_min_count,
        filter_by_expr_parameters = filter_by_expr_parameters
      ),
      error = function(e) {
        list(
          status = "error",
          reason = conditionMessage(e)
        )
      }
    )

    cluster_run_objects[[cluster_id_current]] <- cluster_result

    cluster_status_rows[[cluster_id_current]] <- tibble::tibble(
      cluster_id = cluster_id_current,
      status = cluster_result$status,
      reason = cluster_result$reason
    )

    if (!identical(cluster_result$status, "completed")) {
      message(
        "Cluster ", cluster_id_current,
        " was not completed. Status: ", cluster_result$status,
        "; reason: ", cluster_result$reason
      )
      next
    }

    if (is.null(first_completed_test_definitions)) {
      first_completed_test_definitions <- cluster_result$test_definitions
    }

    raw_count_output <- tibble::tibble(
      ensembl_gene_id = rownames(cluster_result$counts_before_filterByExpr)
    ) |>
      dplyr::left_join(
        mapped_counts$annotation |>
          dplyr::select(ensembl_gene_id, gene, chromosome),
        by = "ensembl_gene_id"
      ) |>
      dplyr::bind_cols(
        tibble::as_tibble(
          as.data.frame(
            cluster_result$counts_before_filterByExpr,
            check.names = FALSE
          ),
          .name_repair = "minimal"
        )
      )

    write_tsv_pseudobulk_per_cluster(
      raw_count_output,
      file.path(
        per_cluster_counts_dir,
        paste0(
          analysis_prefix,
          "_cluster_", cluster_slug,
          "_rawCounts_proteinCodingAutosomalAndX.tsv.gz"
        )
      )
    )

    sample_level_cpm_output <- tibble::tibble(
      ensembl_gene_id = rownames(cluster_result$tmm_cpm_matrix)
    ) |>
      dplyr::left_join(
        mapped_counts$annotation |>
          dplyr::select(ensembl_gene_id, gene, chromosome),
        by = "ensembl_gene_id"
      ) |>
      dplyr::select(
        ensembl_gene_id,
        gene,
        chromosome
      ) |>
      dplyr::bind_cols(
        tibble::as_tibble(
          as.data.frame(
            cluster_result$tmm_cpm_matrix,
            check.names = FALSE
          ),
          .name_repair = "minimal"
        )
      )

    sample_level_detection_output <- tibble::tibble(
      ensembl_gene_id = rownames(cluster_result$detection_percentage_matrix)
    ) |>
      dplyr::left_join(
        mapped_counts$annotation |>
          dplyr::select(ensembl_gene_id, gene, chromosome),
        by = "ensembl_gene_id"
      ) |>
      dplyr::select(
        ensembl_gene_id,
        gene,
        chromosome
      ) |>
      dplyr::bind_cols(
        tibble::as_tibble(
          as.data.frame(
            cluster_result$detection_percentage_matrix,
            check.names = FALSE
          ),
          .name_repair = "minimal"
        )
      )

    write_tsv_pseudobulk_per_cluster(
      sample_level_cpm_output,
      file.path(
        per_cluster_sample_metrics_dir,
        paste0(
          analysis_prefix,
          "_cluster_", cluster_slug,
          "_sampleLevel_TMMnormalizedCPM_testedGenes.tsv.gz"
        )
      )
    )

    write_tsv_pseudobulk_per_cluster(
      sample_level_detection_output,
      file.path(
        per_cluster_sample_metrics_dir,
        paste0(
          analysis_prefix,
          "_cluster_", cluster_slug,
          "_sampleLevel_percentPositiveSpots_testedGenes.tsv.gz"
        )
      )
    )

    write_tsv_pseudobulk_per_cluster(
      cluster_result$library_summary,
      file.path(
        per_cluster_sample_metrics_dir,
        paste0(
          analysis_prefix,
          "_cluster_", cluster_slug,
          "_sampleLevel_libraryAndGroupMetadata.tsv"
        )
      )
    )

    write_tsv_pseudobulk_per_cluster(
      cluster_result$cluster_group_metadata,
      file.path(
        per_cluster_design_dir,
        paste0(
          analysis_prefix,
          "_cluster_", cluster_slug,
          "_sampleInclusionAndPseudobulkSummary.tsv"
        )
      )
    )

    write_tsv_pseudobulk_per_cluster(
      cluster_result$design_output,
      file.path(
        per_cluster_design_dir,
        paste0(
          analysis_prefix,
          "_cluster_", cluster_slug,
          "_designMatrix.tsv"
        )
      )
    )

    write_tsv_pseudobulk_per_cluster(
      cluster_result$contrast_output,
      file.path(
        per_cluster_design_dir,
        paste0(
          analysis_prefix,
          "_cluster_", cluster_slug,
          "_contrastMatrix.tsv"
        )
      )
    )

    write_tsv_pseudobulk_per_cluster(
      cluster_result$filter_status,
      file.path(
        per_cluster_diagnostics_dir,
        paste0(
          analysis_prefix,
          "_cluster_", cluster_slug,
          "_filterByExprStatus.tsv.gz"
        )
      )
    )

    write_tsv_pseudobulk_per_cluster(
      cluster_result$library_summary,
      file.path(
        per_cluster_diagnostics_dir,
        paste0(
          analysis_prefix,
          "_cluster_", cluster_slug,
          "_librarySizeSummary.tsv"
        )
      )
    )

    write_tsv_pseudobulk_per_cluster(
      cluster_result$gene_abundance_qc,
      file.path(
        per_cluster_diagnostics_dir,
        paste0(
          analysis_prefix,
          "_cluster_", cluster_slug,
          "_geneAbundanceAndDetectionQC.tsv.gz"
        )
      )
    )

    write_tsv_pseudobulk_per_cluster(
      cluster_result$analysis_summary,
      file.path(
        per_cluster_diagnostics_dir,
        paste0(
          analysis_prefix,
          "_cluster_", cluster_slug,
          "_edgeRAnalysisSummary.tsv"
        )
      )
    )

    cluster_xlsx_file <- file.path(
      per_cluster_xlsx_dir,
      paste0(
        analysis_prefix,
        "_cluster_", cluster_slug,
        "_fullResults_nineSheets.xlsx"
      )
    )

    write_result_workbook_per_cluster(
      result_list = cluster_result$results,
      test_definitions = cluster_result$test_definitions,
      output_file = cluster_xlsx_file,
      include_cluster_id = FALSE,
      xlsx_max_decimal_places = xlsx_max_decimal_places
    )

    message("Completed cluster: ", cluster_id_current)
    message("Cluster XLSX: ", cluster_xlsx_file)
  }

  cluster_status <- dplyr::bind_rows(cluster_status_rows)

  write_tsv_pseudobulk_per_cluster(
    cluster_status,
    file.path(
      output_subdirectories$diagnostics,
      paste0(analysis_prefix, "_clusterAnalysisStatus.tsv")
    )
  )

  completed_cluster_ids <- cluster_status |>
    dplyr::filter(status == "completed") |>
    dplyr::pull(cluster_id)

  if (length(completed_cluster_ids) == 0L) {
    stop("No cluster completed the edgeR analysis.", call. = FALSE)
  }

  combined_results <- lapply(
    first_completed_test_definitions$test_id,
    function(test_id) {
      combined_test_result <- dplyr::bind_rows(
        lapply(
          completed_cluster_ids,
          function(cluster_id_current) {
            cluster_run_objects[[cluster_id_current]]$results[[test_id]] |>
              dplyr::mutate(
                cluster_id = as.character(cluster_id_current),
                .before = 1
              )
          }
        )
      )

      combined_test_result |>
        dplyr::mutate(
          cluster_id = factor(cluster_id, levels = cluster_ids)
        ) |>
        dplyr::arrange(cluster_id, PValue, dplyr::desc(abs(logFC))) |>
        dplyr::mutate(cluster_id = as.character(cluster_id))
    }
  )

  names(combined_results) <- first_completed_test_definitions$test_id

  combined_xlsx_file <- file.path(
    combined_results_dir,
    paste0(
      analysis_prefix,
      "_allClusters_fullResults_nineSheets.xlsx"
    )
  )

  write_result_workbook_per_cluster(
    result_list = combined_results,
    test_definitions = first_completed_test_definitions,
    output_file = combined_xlsx_file,
    include_cluster_id = TRUE,
    xlsx_max_decimal_places = xlsx_max_decimal_places
  )

  for (test_id in first_completed_test_definitions$test_id) {
    write_tsv_pseudobulk_per_cluster(
      combined_results[[test_id]],
      file.path(
        combined_results_dir,
        paste0(
          analysis_prefix,
          "_allClusters_", test_id,
          "_fullResults.tsv.gz"
        )
      )
    )
  }

  combined_gene_abundance_qc <- dplyr::bind_rows(
    lapply(
      completed_cluster_ids,
      function(cluster_id_current) {
        cluster_run_objects[[cluster_id_current]]$gene_abundance_qc
      }
    )
  ) |>
    dplyr::mutate(
      cluster_id = factor(cluster_id, levels = cluster_ids)
    ) |>
    dplyr::arrange(
      .data$cluster_id,
      .data$ensembl_gene_id
    ) |>
    dplyr::mutate(cluster_id = as.character(cluster_id))

  combined_gene_abundance_qc_file <- file.path(
    output_subdirectories$diagnostics,
    paste0(
      analysis_prefix,
      "_geneAbundanceAndDetectionQC_allClusters.tsv.gz"
    )
  )

  write_tsv_pseudobulk_per_cluster(
    combined_gene_abundance_qc,
    combined_gene_abundance_qc_file
  )

  combined_analysis_summary <- dplyr::bind_rows(
    lapply(
      completed_cluster_ids,
      function(cluster_id_current) {
        cluster_run_objects[[cluster_id_current]]$analysis_summary
      }
    )
  ) |>
    dplyr::mutate(
      cluster_id = factor(cluster_id, levels = cluster_ids)
    ) |>
    dplyr::arrange(
      .data$cluster_id,
      match(.data$test_id, first_completed_test_definitions$test_id)
    ) |>
    dplyr::mutate(cluster_id = as.character(cluster_id))

  write_tsv_pseudobulk_per_cluster(
    combined_analysis_summary,
    file.path(
      output_subdirectories$diagnostics,
      paste0(analysis_prefix, "_edgeRAnalysisSummary_allClusters.tsv")
    )
  )

  write_tsv_pseudobulk_per_cluster(
    first_completed_test_definitions,
    file.path(
      output_subdirectories$design_metadata,
      paste0(analysis_prefix, "_testDefinitions.tsv")
    )
  )

  descriptive_metric_definitions <- tibble::tribble(
    ~metric, ~definition, ~weighting_and_scope,
    "mean_percent_positive_spots",
    paste0(
      "For each sample and cluster: 100 * number of spots with raw count >= ",
      positive_spot_min_count,
      " divided by the number of included spots. The displayed value is the mean of sample-level percentages."
    ),
    paste0(
      "Each sample has equal weight. Only samples with at least ",
      min_spots_per_sample_cluster,
      " spots in the cluster are included. Overall effects use equal-weight marginal means across the other factor."
    ),
    "sd_percent_positive_spots",
    "Sample standard deviation of the sample-level percent-positive values within the displayed biological level.",
    "Calculated across included biological samples; descriptive only.",
    "mean_TMM_CPM",
    paste0(
      "Arithmetic mean of sample-level CPM values calculated by edgeR::cpm from the TMM-normalized DGEList with normalized.lib.sizes=TRUE."
    ),
    paste0(
      "Each sample has equal weight. Overall effects use equal-weight marginal means across the other factor. This is descriptive and does not change the edgeR model or filtering."
    ),
    "sd_TMM_CPM",
    "Sample standard deviation of sample-level TMM-normalized CPM values within the displayed biological level.",
    "Calculated across included biological samples; descriptive only.",
    "samples_with_positive_spots",
    paste0(
      "Number of included samples in the sex-by-group cell with at least one spot having raw count >= ",
      positive_spot_min_count,
      "."
    ),
    "Saved in the detailed geneAbundanceAndDetectionQC tables, not added to the nine main result sheets.",
    "total_positive_spots",
    paste0(
      "Total number of positive spots across included samples in the sex-by-group cell, using raw count >= ",
      positive_spot_min_count,
      "."
    ),
    "Saved in the detailed geneAbundanceAndDetectionQC tables, not added to the nine main result sheets.",
    "median_percent_positive_spots_and_median_TMM_CPM",
    "Median sample-level detection percentage and median sample-level TMM-normalized CPM within each sex-by-group cell.",
    "Saved in the detailed geneAbundanceAndDetectionQC tables to reveal one-sample-driven or strongly skewed patterns.",
    "logCPM",
    "The standard edgeR abundance column returned with the QL test result; it is not a group-specific mean.",
    "Retained unchanged from edgeR::topTags."
  )

  descriptive_metric_definitions_file <- file.path(
    output_subdirectories$design_metadata,
    paste0(analysis_prefix, "_descriptiveMetricDefinitions.tsv")
  )

  write_tsv_pseudobulk_per_cluster(
    descriptive_metric_definitions,
    descriptive_metric_definitions_file
  )

  edgeR_perCluster_results <- lapply(
    completed_cluster_ids,
    function(cluster_id_current) {
      cluster_run_objects[[cluster_id_current]]$results
    }
  )
  names(edgeR_perCluster_results) <- completed_cluster_ids

  edgeR_perCluster_combinedResults <- combined_results
  edgeR_perCluster_analysisSummary <- combined_analysis_summary
  edgeR_perCluster_clusterStatus <- cluster_status
  edgeR_perCluster_testDefinitions <- first_completed_test_definitions

  edgeR_perCluster_geneAbundanceQC <- lapply(
    completed_cluster_ids,
    function(cluster_id_current) {
      cluster_run_objects[[cluster_id_current]]$gene_abundance_qc
    }
  )
  names(edgeR_perCluster_geneAbundanceQC) <- completed_cluster_ids

  edgeR_perCluster_combinedGeneAbundanceQC <- combined_gene_abundance_qc
  edgeR_perCluster_descriptiveMetricDefinitions <-
    descriptive_metric_definitions

  result_rdata_file <- file.path(
    output_subdirectories$edger_results,
    paste0(analysis_prefix, "_edgeRResults.RData")
  )

  save(
    edgeR_perCluster_results,
    edgeR_perCluster_combinedResults,
    edgeR_perCluster_analysisSummary,
    edgeR_perCluster_clusterStatus,
    edgeR_perCluster_testDefinitions,
    edgeR_perCluster_geneAbundanceQC,
    edgeR_perCluster_combinedGeneAbundanceQC,
    edgeR_perCluster_descriptiveMetricDefinitions,
    file = result_rdata_file,
    compress = "xz"
  )

  edgeR_perCluster_modelObjects <- lapply(
    completed_cluster_ids,
    function(cluster_id_current) {
      cluster_object <- cluster_run_objects[[cluster_id_current]]

      list(
        dge = cluster_object$dge,
        fit = cluster_object$fit,
        test_objects = cluster_object$test_objects,
        design = cluster_object$design,
        contrasts = cluster_object$contrasts,
        keep_genes = cluster_object$keep_genes
      )
    }
  )
  names(edgeR_perCluster_modelObjects) <- completed_cluster_ids

  model_rdata_file <- file.path(
    output_subdirectories$objects,
    paste0(analysis_prefix, "_edgeRModelObjects.RData")
  )

  save(
    edgeR_perCluster_modelObjects,
    file = model_rdata_file,
    compress = "xz"
  )

  analysis_parameters <- c(
    paste0("generated=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("project_root=", project_root),
    paste0("dataset_name=", dataset_name),
    paste0("clustering_name=", clustering_name),
    paste0("input_rdata_file=", input_rdata_file),
    paste0("seurat_object_name=", seurat_object_name_resolved),
    paste0("assay_name=", assay_name),
    paste0("count_layers=", paste(count_extraction$count_layers, collapse = ",")),
    paste0("sample_column=", sample_column),
    paste0("cluster_column=", cluster_column_resolved),
    paste0("expected_number_of_samples=", expected_number_of_samples),
    paste0("expected_number_of_clusters=", expected_number_of_clusters),
    paste0("completed_clusters=", paste(completed_cluster_ids, collapse = ",")),
    paste0("min_spots_per_sample_cluster=", min_spots_per_sample_cluster),
    paste0("min_samples_per_factor_cell=", min_samples_per_factor_cell),
    paste0("positive_spot_min_count=", positive_spot_min_count),
    paste0("xlsx_max_decimal_places=", xlsx_max_decimal_places),
    "descriptive_detection_metric=mean and sample SD of sample-level percent-positive spots",
    "descriptive_abundance_metric=mean and sample SD of sample-level TMM-normalized CPM",
    "descriptive_metrics_used_for_filtering_or_testing=FALSE",
    "overall_descriptive_means=equal-weight marginal means across the other factor",
    paste0("filterByExpr_min.count=", filter_by_expr_parameters$min.count),
    paste0("filterByExpr_min.total.count=", filter_by_expr_parameters$min.total.count),
    paste0("filterByExpr_large.n=", filter_by_expr_parameters$large.n),
    paste0("filterByExpr_min.prop=", filter_by_expr_parameters$min.prop),
    paste0("feature_identifier_type=", mapped_counts$feature_identifier_type),
    "retained_gene_biotype=protein_coding",
    "retained_chromosomes=1-19,X",
    "excluded_chromosomes=Y,MT,non-canonical",
    "model=expression ~ fmt_donor_group * sex",
    "number_of_tests_per_cluster=9",
    paste0("combined_xlsx=", combined_xlsx_file),
    paste0("result_rdata=", result_rdata_file),
    paste0("model_rdata=", model_rdata_file)
  )

  writeLines(
    analysis_parameters,
    con = file.path(
      output_subdirectories$logs,
      paste0(analysis_prefix, "_analysisParameters.txt")
    )
  )

  writeLines(
    capture.output(utils::sessionInfo()),
    con = file.path(
      output_subdirectories$logs,
      paste0(analysis_prefix, "_sessionInfo.txt")
    )
  )

  expected_main_files <- c(
    combined_xlsx_file,
    result_rdata_file,
    model_rdata_file,
    file.path(
      output_subdirectories$diagnostics,
      paste0(analysis_prefix, "_clusterAnalysisStatus.tsv")
    ),
    combined_gene_abundance_qc_file,
    descriptive_metric_definitions_file
  )

  missing_main_files <- expected_main_files[!file.exists(expected_main_files)]

  if (length(missing_main_files) > 0L) {
    stop(
      "Missing main output files:\n",
      paste(missing_main_files, collapse = "\n"),
      call. = FALSE
    )
  }

  error_clusters <- cluster_status |>
    dplyr::filter(status == "error")

  message("\n", paste(rep("=", 80), collapse = ""))
  message("PSEUDOBULK PER-CLUSTER edgeR ANALYSIS FINISHED")
  message(paste(rep("=", 80), collapse = ""))
  message("Clusters discovered: ", length(cluster_ids))
  message("Clusters completed: ", length(completed_cluster_ids))
  message("Combined XLSX: ", combined_xlsx_file)
  message("Result RData: ", result_rdata_file)
  message("Model-object RData: ", model_rdata_file)
  message(
    "Combined abundance/detection QC: ",
    combined_gene_abundance_qc_file
  )
  message(
    "Descriptive metric definitions: ",
    descriptive_metric_definitions_file
  )
  message(
    "Per-cluster sample-level CPM and detection metrics: ",
    per_cluster_sample_metrics_dir
  )

  if (nrow(error_clusters) > 0L && isTRUE(stop_on_cluster_error)) {
    stop(
      "At least one cluster failed with an error. See clusterAnalysisStatus.tsv.",
      call. = FALSE
    )
  }

  invisible(
    list(
      output_statistics_dir = output_statistics_dir,
      cluster_ids = cluster_ids,
      completed_cluster_ids = completed_cluster_ids,
      cluster_status = cluster_status,
      combined_xlsx_file = combined_xlsx_file,
      result_rdata_file = result_rdata_file,
      model_rdata_file = model_rdata_file,
      combined_gene_abundance_qc_file =
        combined_gene_abundance_qc_file,
      descriptive_metric_definitions_file =
        descriptive_metric_definitions_file,
      per_cluster_sample_metrics_dir =
        per_cluster_sample_metrics_dir,
      results = edgeR_perCluster_results,
      combined_results = edgeR_perCluster_combinedResults,
      analysis_summary = edgeR_perCluster_analysisSummary,
      gene_abundance_qc = edgeR_perCluster_geneAbundanceQC,
      combined_gene_abundance_qc =
        edgeR_perCluster_combinedGeneAbundanceQC,
      descriptive_metric_definitions =
        edgeR_perCluster_descriptiveMetricDefinitions
    )
  )
}

# ==============================================================================
# End of functions file
# ==============================================================================
