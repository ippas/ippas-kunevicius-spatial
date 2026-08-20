# ==============================================================================
# functions_logNormalizeVst_preprocessing.R
#
# Purpose:
# Reusable functions for preparing a Seurat v5 parent object for the
# LogNormalize + VST branch of the maternalFMT spatial clustering analysis.
#
# Main workflow:
#   read metadata -> exclude samples -> read Space Ranger counts/coordinates
#   -> merge samples -> LogNormalize -> VST HVG selection -> ScaleData -> PCA
#
# Notes:
# - Raw RNA counts remain available in sample-specific Seurat v5 layers.
# - No integration and no clustering are performed in this file.
# - Spatial images are loaded only when the standard Space Ranger files needed
#   by Load10X_Spatial() are present. Otherwise, counts and coordinates are
#   loaded without stopping the analysis.
# ==============================================================================  


# ==============================================================================
# 1. General validation helpers
# ==============================================================================

check_required_packages <- function(
    packages = c("Seurat", "SeuratObject", "Matrix", "ggplot2")
) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


assert_seurat_v5 <- function() {
  seurat_version <- utils::packageVersion("Seurat")
  seurat_object_version <- utils::packageVersion("SeuratObject")

  message("Seurat version: ", as.character(seurat_version))
  message("SeuratObject version: ", as.character(seurat_object_version))

  if (seurat_version < "5.0.0") {
    stop(
      "This workflow requires Seurat v5 or newer. Current version: ",
      as.character(seurat_version),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


create_output_directories <- function(output_root) {
  subdirectories <- c("RData", "tables", "plots", "config", "logs")

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  output_paths <- stats::setNames(
    file.path(output_root, subdirectories),
    subdirectories
  )

  invisible(
    lapply(
      output_paths,
      dir.create,
      recursive = TRUE,
      showWarnings = FALSE
    )
  )

  output_paths
}


write_tsv <- function(data, file) {
  utils::write.table(
    data,
    file = file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA"
  )

  invisible(file)
}


escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\?.])", "\\\\\\1", x)
}


# ==============================================================================
# 2. Metadata handling
# ==============================================================================

read_and_filter_sample_metadata <- function(
    metadata_file,
    excluded_sample_ids = character(0),
    sample_id_col = "sample_ID",
    expected_n_samples = NULL
) {
  if (!file.exists(metadata_file)) {
    stop("Metadata file does not exist: ", metadata_file, call. = FALSE)
  }

  metadata <- utils::read.delim(
    metadata_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!sample_id_col %in% colnames(metadata)) {
    stop(
      "Metadata does not contain column `", sample_id_col, "`.",
      call. = FALSE
    )
  }

  sample_ids <- metadata[[sample_id_col]]

  if (anyNA(sample_ids) || any(sample_ids == "")) {
    stop("Metadata contains missing or empty sample IDs.", call. = FALSE)
  }

  if (anyDuplicated(sample_ids) > 0) {
    duplicated_ids <- unique(sample_ids[duplicated(sample_ids)])
    stop(
      "Duplicated sample IDs in metadata: ",
      paste(duplicated_ids, collapse = ", "),
      call. = FALSE
    )
  }

  missing_excluded <- setdiff(excluded_sample_ids, sample_ids)

  if (length(missing_excluded) > 0) {
    stop(
      "The following requested exclusions are absent from metadata: ",
      paste(missing_excluded, collapse = ", "),
      call. = FALSE
    )
  }

  filtered_metadata <- metadata[
    !metadata[[sample_id_col]] %in% excluded_sample_ids,
    ,
    drop = FALSE
  ]

  if (!is.null(expected_n_samples) && nrow(filtered_metadata) != expected_n_samples) {
    stop(
      "Expected ", expected_n_samples, " samples after filtering, but found ",
      nrow(filtered_metadata), ".",
      call. = FALSE
    )
  }

  rownames(filtered_metadata) <- filtered_metadata[[sample_id_col]]

  message(
    "Metadata: retained ", nrow(filtered_metadata), " samples; excluded ",
    length(excluded_sample_ids), "."
  )

  filtered_metadata
}


# ==============================================================================
# 3. Space Ranger path detection
# ==============================================================================

find_sample_directory <- function(path_to_data, sample_id) {
  if (!dir.exists(path_to_data)) {
    stop("Data root does not exist: ", path_to_data, call. = FALSE)
  }

  direct_candidate <- file.path(path_to_data, sample_id)

  if (dir.exists(direct_candidate)) {
    return(normalizePath(direct_candidate, mustWork = TRUE))
  }

  all_directories <- list.dirs(
    path_to_data,
    recursive = TRUE,
    full.names = TRUE
  )

  exact_matches <- all_directories[basename(all_directories) == sample_id]

  if (length(exact_matches) == 1) {
    return(normalizePath(exact_matches, mustWork = TRUE))
  }

  if (length(exact_matches) > 1) {
    stop(
      "Multiple exact directories found for sample `", sample_id, "`: ",
      paste(exact_matches, collapse = "; "),
      call. = FALSE
    )
  }

  sample_pattern <- paste0("^", escape_regex(sample_id), "($|[-_.])")
  prefix_matches <- all_directories[
    grepl(sample_pattern, basename(all_directories))
  ]

  if (length(prefix_matches) == 1) {
    warning(
      "Using non-exact directory match for sample `", sample_id, "`: ",
      prefix_matches,
      call. = FALSE
    )
    return(normalizePath(prefix_matches, mustWork = TRUE))
  }

  stop(
    "Could not uniquely identify a data directory for sample `",
    sample_id,
    "` under: ",
    path_to_data,
    call. = FALSE
  )
}


locate_counts_input <- function(sample_directory) {
  direct_h5_candidates <- c(
    file.path(sample_directory, "outs", "filtered_feature_bc_matrix.h5"),
    file.path(sample_directory, "filtered_feature_bc_matrix.h5")
  )

  existing_h5 <- direct_h5_candidates[file.exists(direct_h5_candidates)]

  if (length(existing_h5) == 0) {
    recursive_h5 <- list.files(
      sample_directory,
      pattern = "^filtered_feature_bc_matrix\\.h5$",
      recursive = TRUE,
      full.names = TRUE
    )
    existing_h5 <- recursive_h5[file.exists(recursive_h5)]
  }

  existing_h5 <- unique(normalizePath(existing_h5, mustWork = FALSE))

  if (length(existing_h5) == 1) {
    return(list(type = "h5", path = existing_h5))
  }

  if (length(existing_h5) > 1) {
    stop(
      "Multiple filtered_feature_bc_matrix.h5 files found under: ",
      sample_directory,
      call. = FALSE
    )
  }

  direct_directory_candidates <- c(
    file.path(sample_directory, "outs", "filtered_feature_bc_matrix"),
    file.path(sample_directory, "filtered_feature_bc_matrix")
  )

  existing_directories <- direct_directory_candidates[
    dir.exists(direct_directory_candidates)
  ]

  if (length(existing_directories) == 0) {
    recursive_directories <- list.dirs(
      sample_directory,
      recursive = TRUE,
      full.names = TRUE
    )
    existing_directories <- recursive_directories[
      basename(recursive_directories) == "filtered_feature_bc_matrix"
    ]
  }

  existing_directories <- unique(
    normalizePath(existing_directories, mustWork = FALSE)
  )

  if (length(existing_directories) == 1) {
    return(list(type = "directory", path = existing_directories))
  }

  if (length(existing_directories) > 1) {
    stop(
      "Multiple filtered_feature_bc_matrix directories found under: ",
      sample_directory,
      call. = FALSE
    )
  }

  stop(
    "No filtered Space Ranger count matrix found under: ",
    sample_directory,
    call. = FALSE
  )
}


locate_spatial_directory <- function(sample_directory, counts_input_path = NULL) {
  candidates <- c(
    file.path(sample_directory, "outs", "spatial"),
    file.path(sample_directory, "spatial")
  )

  if (!is.null(counts_input_path)) {
    candidates <- c(candidates, file.path(dirname(counts_input_path), "spatial"))
  }

  existing <- unique(candidates[dir.exists(candidates)])

  if (length(existing) == 1) {
    return(normalizePath(existing, mustWork = TRUE))
  }

  if (length(existing) > 1) {
    position_score <- vapply(
      existing,
      function(path) {
        any(file.exists(file.path(
          path,
          c("tissue_positions.csv", "tissue_positions_list.csv")
        )))
      },
      logical(1)
    )

    scored <- existing[position_score]

    if (length(scored) == 1) {
      return(normalizePath(scored, mustWork = TRUE))
    }
  }

  recursive_directories <- list.dirs(
    sample_directory,
    recursive = TRUE,
    full.names = TRUE
  )
  recursive_spatial <- recursive_directories[
    basename(recursive_directories) == "spatial"
  ]

  valid_recursive <- recursive_spatial[
    vapply(
      recursive_spatial,
      function(path) {
        any(file.exists(file.path(
          path,
          c("tissue_positions.csv", "tissue_positions_list.csv")
        )))
      },
      logical(1)
    )
  ]

  if (length(valid_recursive) == 1) {
    return(normalizePath(valid_recursive, mustWork = TRUE))
  }

  if (length(valid_recursive) > 1) {
    stop(
      "Multiple spatial directories with tissue positions found under: ",
      sample_directory,
      call. = FALSE
    )
  }

  NULL
}


locate_tissue_positions_file <- function(spatial_directory) {
  if (is.null(spatial_directory)) {
    return(NULL)
  }

  candidates <- file.path(
    spatial_directory,
    c("tissue_positions.csv", "tissue_positions_list.csv")
  )

  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0) {
    return(NULL)
  }

  if (length(existing) > 1) {
    existing <- existing[basename(existing) == "tissue_positions.csv"]
  }

  normalizePath(existing[[1]], mustWork = TRUE)
}


# ==============================================================================
# 4. Counts and coordinates
# ==============================================================================

extract_gene_expression_matrix <- function(read10x_result) {
  if (inherits(read10x_result, "Matrix")) {
    return(read10x_result)
  }

  if (!is.list(read10x_result) || length(read10x_result) == 0) {
    stop("Read10X returned an unsupported object.", call. = FALSE)
  }

  if ("Gene Expression" %in% names(read10x_result)) {
    return(read10x_result[["Gene Expression"]])
  }

  matrix_elements <- vapply(
    read10x_result,
    inherits,
    logical(1),
    what = "Matrix"
  )

  if (sum(matrix_elements) == 1) {
    return(read10x_result[[which(matrix_elements)]])
  }

  stop(
    "Could not uniquely identify the Gene Expression matrix in Read10X output.",
    call. = FALSE
  )
}


read_counts_matrix <- function(counts_input) {
  if (identical(counts_input$type, "h5")) {
    result <- Seurat::Read10X_h5(
      filename = counts_input$path,
      use.names = TRUE,
      unique.features = TRUE
    )
  } else if (identical(counts_input$type, "directory")) {
    result <- Seurat::Read10X(
      data.dir = counts_input$path,
      gene.column = 2,
      unique.features = TRUE,
      strip.suffix = FALSE
    )
  } else {
    stop("Unsupported counts input type: ", counts_input$type, call. = FALSE)
  }

  counts <- extract_gene_expression_matrix(result)

  if (nrow(counts) == 0 || ncol(counts) == 0) {
    stop("The count matrix is empty: ", counts_input$path, call. = FALSE)
  }

  counts
}


read_tissue_positions <- function(tissue_positions_file) {
  if (is.null(tissue_positions_file)) {
    return(NULL)
  }

  if (basename(tissue_positions_file) == "tissue_positions_list.csv") {
    positions <- utils::read.csv(
      tissue_positions_file,
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if (ncol(positions) < 6) {
      stop(
        "Unexpected tissue_positions_list.csv format: ",
        tissue_positions_file,
        call. = FALSE
      )
    }

    positions <- positions[, seq_len(6), drop = FALSE]
    colnames(positions) <- c(
      "barcode",
      "in_tissue",
      "array_row",
      "array_col",
      "pxl_row_in_fullres",
      "pxl_col_in_fullres"
    )
  } else {
    positions <- utils::read.csv(
      tissue_positions_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    expected_columns <- c(
      "barcode",
      "in_tissue",
      "array_row",
      "array_col",
      "pxl_row_in_fullres",
      "pxl_col_in_fullres"
    )

    missing_columns <- setdiff(expected_columns, colnames(positions))

    if (length(missing_columns) > 0) {
      stop(
        "Missing columns in tissue positions file: ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }

    positions <- positions[, expected_columns, drop = FALSE]
  }

  if (anyDuplicated(positions$barcode) > 0) {
    stop(
      "Duplicated barcodes in tissue positions file: ",
      tissue_positions_file,
      call. = FALSE
    )
  }

  rownames(positions) <- positions$barcode
  positions
}


can_load_standard_spatial_image <- function(counts_input, spatial_directory) {
  if (!identical(counts_input$type, "h5") || is.null(spatial_directory)) {
    return(FALSE)
  }

  required_files <- c(
    file.path(spatial_directory, "scalefactors_json.json"),
    file.path(spatial_directory, "tissue_lowres_image.png")
  )

  all(file.exists(required_files)) &&
    identical(
      normalizePath(dirname(spatial_directory), mustWork = FALSE),
      normalizePath(dirname(counts_input$path), mustWork = FALSE)
    )
}


add_metadata_row_to_object <- function(object, metadata_row, sample_id_col) {
  for (column_name in colnames(metadata_row)) {
    object[[column_name]] <- rep(
      metadata_row[[column_name]][[1]],
      ncol(object)
    )
  }

  if (!sample_id_col %in% colnames(object[[]])) {
    stop("Failed to add sample metadata to Seurat object.", call. = FALSE)
  }

  object
}


add_coordinates_to_object <- function(object, positions) {
  coordinate_columns <- c(
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"
  )

  if (is.null(positions)) {
    for (column_name in coordinate_columns) {
      object[[column_name]] <- NA_real_
    }
    return(object)
  }

  matched_positions <- positions[
    match(colnames(object), positions$barcode),
    ,
    drop = FALSE
  ]

  for (column_name in coordinate_columns) {
    object[[column_name]] <- matched_positions[[column_name]]
  }

  n_missing <- sum(is.na(matched_positions$barcode))

  if (n_missing > 0) {
    warning(
      n_missing,
      " barcodes in the count matrix were absent from the tissue positions file.",
      call. = FALSE
    )
  }

  object
}


read_one_spatial_sample <- function(
    sample_id,
    metadata_row,
    path_to_data,
    sample_id_col = "sample_ID",
    mitochondrial_pattern = "^mt-",
    min_cells = 0,
    min_features = 0,
    load_images_when_possible = TRUE,
    verbose = TRUE
) {
  if (verbose) {
    message("Reading sample: ", sample_id)
  }

  sample_directory <- find_sample_directory(path_to_data, sample_id)
  counts_input <- locate_counts_input(sample_directory)
  spatial_directory <- locate_spatial_directory(
    sample_directory,
    counts_input_path = counts_input$path
  )
  positions_file <- locate_tissue_positions_file(spatial_directory)

  image_loaded <- FALSE

  if (
    load_images_when_possible &&
      can_load_standard_spatial_image(counts_input, spatial_directory)
  ) {
    object <- tryCatch(
      Seurat::Load10X_Spatial(
        data.dir = dirname(counts_input$path),
        filename = basename(counts_input$path),
        assay = "RNA",
        slice = sample_id,
        filter.matrix = TRUE,
        to.upper = FALSE
      ),
      error = function(error) {
        warning(
          "Load10X_Spatial failed for sample `", sample_id,
          "`; falling back to counts + coordinate metadata. Error: ",
          conditionMessage(error),
          call. = FALSE
        )
        NULL
      }
    )

    image_loaded <- !is.null(object)
  } else {
    object <- NULL
  }

  if (is.null(object)) {
    counts <- read_counts_matrix(counts_input)

    object <- Seurat::CreateSeuratObject(
      counts = counts,
      project = sample_id,
      assay = "RNA",
      min.cells = min_cells,
      min.features = min_features
    )
  }

  # Store standard Visium coordinates in meta.data regardless of whether a
  # Seurat spatial image object could be loaded.
  positions <- read_tissue_positions(positions_file)
  object <- add_coordinates_to_object(object, positions)

  object <- add_metadata_row_to_object(
    object = object,
    metadata_row = metadata_row,
    sample_id_col = sample_id_col
  )

  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    object,
    pattern = mitochondrial_pattern,
    assay = "RNA"
  )

  object <- SeuratObject::RenameCells(
    object,
    add.cell.id = sample_id
  )

  SeuratObject::DefaultAssay(object) <- "RNA"
  SeuratObject::Project(object) <- sample_id

  manifest_row <- data.frame(
    sample_ID = sample_id,
    sample_directory = sample_directory,
    counts_input_type = counts_input$type,
    counts_input_path = counts_input$path,
    spatial_directory = if (is.null(spatial_directory)) NA_character_ else spatial_directory,
    tissue_positions_file = if (is.null(positions_file)) NA_character_ else positions_file,
    spatial_image_loaded = image_loaded,
    n_features = nrow(object),
    n_spots = ncol(object),
    stringsAsFactors = FALSE
  )

  list(object = object, manifest = manifest_row)
}


read_spatial_samples_from_metadata <- function(
    metadata,
    path_to_data,
    sample_id_col = "sample_ID",
    mitochondrial_pattern = "^mt-",
    min_cells = 0,
    min_features = 0,
    load_images_when_possible = TRUE,
    verbose = TRUE
) {
  sample_ids <- metadata[[sample_id_col]]

  results <- lapply(
    sample_ids,
    function(sample_id) {
      read_one_spatial_sample(
        sample_id = sample_id,
        metadata_row = metadata[sample_id, , drop = FALSE],
        path_to_data = path_to_data,
        sample_id_col = sample_id_col,
        mitochondrial_pattern = mitochondrial_pattern,
        min_cells = min_cells,
        min_features = min_features,
        load_images_when_possible = load_images_when_possible,
        verbose = verbose
      )
    }
  )

  samples_list <- stats::setNames(
    lapply(results, function(x) x$object),
    sample_ids
  )

  input_manifest <- do.call(
    rbind,
    lapply(results, function(x) x$manifest)
  )

  rownames(input_manifest) <- NULL

  list(
    samples_list = samples_list,
    input_manifest = input_manifest
  )
}


# ==============================================================================
# 5. QC summaries
# ==============================================================================

safe_numeric_summary <- function(x, function_name) {
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  switch(
    function_name,
    min = min(x),
    median = stats::median(x),
    mean = mean(x),
    max = max(x),
    stop("Unsupported summary function: ", function_name, call. = FALSE)
  )
}


summarize_qc_by_sample <- function(samples_list) {
  summaries <- lapply(
    names(samples_list),
    function(sample_id) {
      object <- samples_list[[sample_id]]
      metadata <- object[[]]

      data.frame(
        sample_ID = sample_id,
        n_spots = ncol(object),
        n_features_total = nrow(object),
        min_nCount_RNA = safe_numeric_summary(metadata$nCount_RNA, "min"),
        median_nCount_RNA = safe_numeric_summary(metadata$nCount_RNA, "median"),
        mean_nCount_RNA = safe_numeric_summary(metadata$nCount_RNA, "mean"),
        max_nCount_RNA = safe_numeric_summary(metadata$nCount_RNA, "max"),
        min_nFeature_RNA = safe_numeric_summary(metadata$nFeature_RNA, "min"),
        median_nFeature_RNA = safe_numeric_summary(metadata$nFeature_RNA, "median"),
        mean_nFeature_RNA = safe_numeric_summary(metadata$nFeature_RNA, "mean"),
        max_nFeature_RNA = safe_numeric_summary(metadata$nFeature_RNA, "max"),
        median_percent_mt = safe_numeric_summary(metadata$percent.mt, "median"),
        mean_percent_mt = safe_numeric_summary(metadata$percent.mt, "mean"),
        max_percent_mt = safe_numeric_summary(metadata$percent.mt, "max"),
        n_in_tissue = if ("in_tissue" %in% colnames(metadata)) {
          sum(metadata$in_tissue == 1, na.rm = TRUE)
        } else {
          NA_integer_
        },
        stringsAsFactors = FALSE
      )
    }
  )

  result <- do.call(rbind, summaries)
  rownames(result) <- NULL
  result
}


# ==============================================================================
# 6. Merge and LogNormalize + VST preprocessing
# ==============================================================================

merge_spatial_samples_v5 <- function(
    samples_list,
    project = "MaternalFMT_n16samples",
    sample_id_col = "sample_ID",
    verbose = TRUE
) {
  if (!is.list(samples_list) || length(samples_list) < 2) {
    stop("`samples_list` must contain at least two Seurat objects.", call. = FALSE)
  }

  sample_ids <- names(samples_list)

  if (is.null(sample_ids) || any(sample_ids == "")) {
    stop("`samples_list` must be named using sample IDs.", call. = FALSE)
  }

  if (verbose) {
    message("Merging ", length(samples_list), " samples...")
  }

  merged_object <- merge(
    x = samples_list[[1]],
    y = samples_list[-1],
    project = project,
    merge.data = FALSE,
    collapse = FALSE
  )

  SeuratObject::DefaultAssay(merged_object) <- "RNA"

  if (!sample_id_col %in% colnames(merged_object[[]])) {
    stop(
      "Merged object does not contain metadata column `",
      sample_id_col,
      "`.",
      call. = FALSE
    )
  }

  count_layers <- SeuratObject::Layers(
    merged_object[["RNA"]],
    search = "^counts"
  )

  if (length(count_layers) == 1) {
    if (verbose) {
      message("Splitting the RNA assay into sample-specific layers...")
    }

    merged_object[["RNA"]] <- split(
      merged_object[["RNA"]],
      f = merged_object[[sample_id_col, drop = TRUE]]
    )
  }

  count_layers <- SeuratObject::Layers(
    merged_object[["RNA"]],
    search = "^counts"
  )

  if (length(count_layers) < 2) {
    stop(
      "The merged RNA assay does not contain sample-specific count layers.",
      call. = FALSE
    )
  }

  if (verbose) {
    message("Merged object contains ", length(count_layers), " RNA count layers.")
  }

  merged_object
}


prepare_logNormalize_vst_parent_object <- function(
    samples_list,
    nfeatures = 2000,
    scale_factor = 10000,
    npcs = 50,
    project = "MaternalFMT_n16samples",
    sample_id_col = "sample_ID",
    seed = 7,
    verbose = TRUE
) {
  if (nfeatures < 1 || npcs < 2 || scale_factor <= 0) {
    stop("Invalid preprocessing parameters.", call. = FALSE)
  }

  object <- merge_spatial_samples_v5(
    samples_list = samples_list,
    project = project,
    sample_id_col = sample_id_col,
    verbose = verbose
  )

  if (verbose) {
    message("Running LogNormalize with scale.factor = ", scale_factor, "...")
  }

  object <- Seurat::NormalizeData(
    object = object,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = scale_factor,
    verbose = verbose
  )

  if (verbose) {
    message("Selecting ", nfeatures, " HVGs using VST...")
  }

  object <- Seurat::FindVariableFeatures(
    object = object,
    assay = "RNA",
    selection.method = "vst",
    nfeatures = nfeatures,
    verbose = verbose
  )

  variable_features <- SeuratObject::VariableFeatures(object[["RNA"]])

  if (length(variable_features) == 0) {
    stop("No variable features were identified.", call. = FALSE)
  }

  if (verbose) {
    message("Scaling selected HVGs without regressing covariates...")
  }

  object <- Seurat::ScaleData(
    object = object,
    assay = "RNA",
    features = variable_features,
    vars.to.regress = NULL,
    do.center = TRUE,
    do.scale = TRUE,
    verbose = verbose
  )

  set.seed(seed)

  if (verbose) {
    message("Running PCA with ", npcs, " components...")
  }

  object <- Seurat::RunPCA(
    object = object,
    assay = "RNA",
    features = variable_features,
    npcs = npcs,
    reduction.name = "pca",
    reduction.key = "PC_",
    seed.use = seed,
    verbose = verbose
  )

  object@misc$logNormalizeVst_preprocessing <- list(
    normalization_method = "LogNormalize",
    scale_factor = scale_factor,
    variable_feature_method = "vst",
    nfeatures_requested = nfeatures,
    nfeatures_identified = length(variable_features),
    vars_to_regress = NULL,
    npcs = npcs,
    seed = seed,
    project = project,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  )

  list(
    object = object,
    variable_features = variable_features
  )
}


# ==============================================================================
# 7. Output tables and plots
# ==============================================================================

create_variable_features_table <- function(object) {
  variable_features <- SeuratObject::VariableFeatures(object[["RNA"]])

  data.frame(
    rank = seq_along(variable_features),
    gene = variable_features,
    stringsAsFactors = FALSE
  )
}


create_pca_standard_deviation_table <- function(object, reduction = "pca") {
  standard_deviations <- SeuratObject::Stdev(object, reduction = reduction)
  variance <- standard_deviations^2
  variance_fraction <- variance / sum(variance)

  data.frame(
    PC = seq_along(standard_deviations),
    standard_deviation = standard_deviations,
    variance = variance,
    variance_fraction = variance_fraction,
    cumulative_variance_fraction = cumsum(variance_fraction),
    stringsAsFactors = FALSE
  )
}


create_preprocessing_summary <- function(
    object,
    n_samples,
    nfeatures_requested,
    scale_factor,
    npcs,
    seed
) {
  count_layers <- SeuratObject::Layers(object[["RNA"]], search = "^counts")
  data_layers <- SeuratObject::Layers(object[["RNA"]], search = "^data")

  data.frame(
    dataset = "maternalFMT_n16samples",
    n_samples = n_samples,
    n_spots = ncol(object),
    n_genes = nrow(object),
    normalization_method = "LogNormalize",
    scale_factor = scale_factor,
    variable_feature_method = "vst",
    nfeatures_requested = nfeatures_requested,
    nfeatures_identified = length(SeuratObject::VariableFeatures(object[["RNA"]])),
    vars_to_regress = "none",
    npcs = npcs,
    n_count_layers = length(count_layers),
    n_data_layers = length(data_layers),
    seed = seed,
    stringsAsFactors = FALSE
  )
}


create_run_parameters_table <- function(
    path_to_data,
    metadata_file,
    excluded_sample_ids,
    output_root,
    nfeatures,
    scale_factor,
    npcs,
    seed,
    mitochondrial_pattern,
    min_cells,
    min_features,
    load_images_when_possible
) {
  values <- list(
    path_to_data = path_to_data,
    metadata_file = metadata_file,
    excluded_sample_ids = paste(excluded_sample_ids, collapse = ","),
    output_root = output_root,
    normalization_method = "LogNormalize",
    scale_factor = scale_factor,
    variable_feature_method = "vst",
    nfeatures = nfeatures,
    vars_to_regress = "none",
    npcs = npcs,
    seed = seed,
    mitochondrial_pattern = mitochondrial_pattern,
    min_cells = min_cells,
    min_features = min_features,
    load_images_when_possible = load_images_when_possible
  )

  data.frame(
    parameter = names(values),
    value = vapply(values, as.character, character(1)),
    stringsAsFactors = FALSE
  )
}


save_preprocessing_plots <- function(
    object,
    output_directory,
    file_prefix,
    npcs = 50
) {
  variable_feature_plot <- Seurat::VariableFeaturePlot(object, assay = "RNA")
  ggplot2::ggsave(
    filename = file.path(
      output_directory,
      paste0(file_prefix, "_variableFeatures.pdf")
    ),
    plot = variable_feature_plot,
    width = 9,
    height = 7,
    units = "in"
  )

  elbow_plot <- Seurat::ElbowPlot(
    object,
    reduction = "pca",
    ndims = npcs
  )
  ggplot2::ggsave(
    filename = file.path(
      output_directory,
      paste0(file_prefix, "_pcaElbow.pdf")
    ),
    plot = elbow_plot,
    width = 9,
    height = 6,
    units = "in"
  )

  pca_sample_plot <- Seurat::DimPlot(
    object,
    reduction = "pca",
    group.by = "sample_ID",
    pt.size = 0.05,
    raster = TRUE
  )
  ggplot2::ggsave(
    filename = file.path(
      output_directory,
      paste0(file_prefix, "_pcaBySample.pdf")
    ),
    plot = pca_sample_plot,
    width = 11,
    height = 8,
    units = "in"
  )

  invisible(TRUE)
}
