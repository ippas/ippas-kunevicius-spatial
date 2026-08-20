# ==============================================================================
# Read individual Space Ranger count matrices into Seurat objects
# ==============================================================================

read_spatial_samples <- function(
    path_to_data,
    metadata,
    sample_id_col = "sample_ID",
    min.cells = 0,
    min.features = 0,
    verbose = TRUE
) {

  if (!is.data.frame(metadata)) {
    stop("`metadata` must be a data.frame.")
  }

  if (!sample_id_col %in% colnames(metadata)) {
    stop("Metadata does not contain column: ", sample_id_col)
  }

  sample_ids <- as.character(metadata[[sample_id_col]])

  if (anyDuplicated(sample_ids) > 0) {
    stop(
      "Duplicated sample IDs in metadata: ",
      paste(unique(sample_ids[duplicated(sample_ids)]), collapse = ", ")
    )
  }

  read_one_sample <- function(sample_id) {

    matrix_dir <- file.path(
      path_to_data,
      sample_id,
      "outs",
      "filtered_feature_bc_matrix"
    )

    if (!dir.exists(matrix_dir)) {
      stop(
        "Count-matrix folder not found for sample ",
        sample_id,
        ": ",
        matrix_dir
      )
    }

    if (verbose) {
      message("Reading sample: ", sample_id)
    }

    counts <- Seurat::Read10X(
      data.dir = matrix_dir,
      gene.column = 2,
      unique.features = TRUE
    )

    # Protection in case Read10X returns several feature types
    if (is.list(counts)) {

      if ("Gene Expression" %in% names(counts)) {
        counts <- counts[["Gene Expression"]]
      } else {
        stop(
          "Read10X returned multiple feature types for sample ",
          sample_id,
          ", but no 'Gene Expression' matrix was found."
        )
      }
    }

    sample_object <- Seurat::CreateSeuratObject(
      counts = counts,
      project = sample_id,
      min.cells = min.cells,
      min.features = min.features
    )

    sample_metadata <- metadata[
      metadata[[sample_id_col]] == sample_id,
      ,
      drop = FALSE
    ]

    cell_metadata <- sample_metadata[
      rep(1, ncol(sample_object)),
      ,
      drop = FALSE
    ]

    rownames(cell_metadata) <- colnames(sample_object)

    sample_object <- Seurat::AddMetaData(
      object = sample_object,
      metadata = cell_metadata
    )

    return(sample_object)
  }

  samples_list <- lapply(sample_ids, read_one_sample)

  names(samples_list) <- sample_ids

  return(samples_list)
}


# ==============================================================================
# Add Space Ranger spatial coordinates to individual Seurat objects
#
# Reads tissue_positions.csv for each sample and adds:
# - in_tissue
# - array_row
# - array_col
# - pxl_row_in_fullres
# - pxl_col_in_fullres
#
# No image files or scalefactors_json.json are required.
# ==============================================================================

add_spatial_coordinates_to_samples <- function(
    samples_list,
    path_to_data,
    verbose = TRUE
) {

  if (!is.list(samples_list) || length(samples_list) == 0) {
    stop("`samples_list` must be a non-empty named list of Seurat objects.")
  }

  sample_ids <- names(samples_list)

  if (is.null(sample_ids) || any(sample_ids == "")) {
    stop("`samples_list` must have sample IDs as list names.")
  }

  required_columns <- c(
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"
  )

  for (sample_id in sample_ids) {

    if (verbose) {
      message("Adding spatial coordinates: ", sample_id)
    }

    spatial_dir <- file.path(
      path_to_data,
      sample_id,
      "outs",
      "spatial"
    )

    positions_file_new <- file.path(
      spatial_dir,
      "tissue_positions.csv"
    )

    positions_file_old <- file.path(
      spatial_dir,
      "tissue_positions_list.csv"
    )

    if (file.exists(positions_file_new)) {

      positions <- utils::read.csv(
        positions_file_new,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      colnames(positions) <- tolower(colnames(positions))

    } else if (file.exists(positions_file_old)) {

      positions <- utils::read.csv(
        positions_file_old,
        header = FALSE,
        stringsAsFactors = FALSE,
        col.names = required_columns
      )

    } else {
      stop(
        "No tissue-position file found for sample `", sample_id, "`.\n",
        "Expected one of:\n",
        positions_file_new, "\n",
        positions_file_old
      )
    }

    missing_columns <- setdiff(required_columns, colnames(positions))

    if (length(missing_columns) > 0) {
      stop(
        "Spatial-position file for sample `", sample_id,
        "` is missing column(s): ",
        paste(missing_columns, collapse = ", ")
      )
    }

    positions <- positions[, required_columns, drop = FALSE]

    positions$barcode <- as.character(positions$barcode)

    if (anyDuplicated(positions$barcode) > 0) {
      stop(
        "Duplicated barcodes in tissue-position file for sample `",
        sample_id, "`."
      )
    }

    sample_object <- samples_list[[sample_id]]

    object_barcodes <- colnames(sample_object)

    barcode_match <- match(
      object_barcodes,
      positions$barcode
    )

    unmatched_barcodes <- object_barcodes[is.na(barcode_match)]

    if (length(unmatched_barcodes) > 0) {
      stop(
        "Spatial coordinates are missing for ",
        length(unmatched_barcodes),
        " barcode(s) in sample `", sample_id, "`.\n",
        "Examples: ",
        paste(head(unmatched_barcodes, 10), collapse = ", ")
      )
    }

    coordinates_for_object <- positions[
      barcode_match,
      ,
      drop = FALSE
    ]

    rownames(coordinates_for_object) <- object_barcodes

    spatial_metadata <- data.frame(
      in_tissue = as.integer(coordinates_for_object$in_tissue),
      array_row = as.integer(coordinates_for_object$array_row),
      array_col = as.integer(coordinates_for_object$array_col),
      pxl_row_in_fullres = as.numeric(
        coordinates_for_object$pxl_row_in_fullres
      ),
      pxl_col_in_fullres = as.numeric(
        coordinates_for_object$pxl_col_in_fullres
      ),
      row.names = object_barcodes,
      check.names = FALSE
    )

    samples_list[[sample_id]] <- SeuratObject::AddMetaData(
      object = sample_object,
      metadata = spatial_metadata
    )
  }

  return(samples_list)
}

# ==============================================================================
# Add QC metrics to individual Seurat objects and summarise samples
# ==============================================================================

add_qc_metrics_to_samples <- function(
    samples_list,
    mitochondrial_pattern = "^mt-",
    verbose = TRUE
) {

  if (!is.list(samples_list) || length(samples_list) == 0) {
    stop("`samples_list` must be a non-empty named list of Seurat objects.")
  }

  sample_ids <- names(samples_list)

  if (is.null(sample_ids) || any(sample_ids == "")) {
    stop("`samples_list` must be a named list of Seurat objects.")
  }

  for (sample_id in sample_ids) {

    if (verbose) {
      message("Calculating QC metrics: ", sample_id)
    }

    sample_object <- samples_list[[sample_id]]

    sample_object[["percent_mt"]] <- Seurat::PercentageFeatureSet(
      object = sample_object,
      pattern = mitochondrial_pattern
    )

    samples_list[[sample_id]] <- sample_object
  }

  qc_summary <- purrr::imap_dfr(
    samples_list,
    function(sample_object, sample_id) {

      metadata <- sample_object[[]]

      tibble::tibble(
        sample_ID = sample_id,
        n_spots = ncol(sample_object),
        n_genes = nrow(sample_object),

        median_UMI_per_spot = median(metadata$nCount_RNA),
        median_genes_per_spot = median(metadata$nFeature_RNA),
        median_percent_mt = median(metadata$percent_mt),

        mean_UMI_per_spot = mean(metadata$nCount_RNA),
        mean_genes_per_spot = mean(metadata$nFeature_RNA),
        mean_percent_mt = mean(metadata$percent_mt),

        min_UMI_per_spot = min(metadata$nCount_RNA),
        max_UMI_per_spot = max(metadata$nCount_RNA),

        min_genes_per_spot = min(metadata$nFeature_RNA),
        max_genes_per_spot = max(metadata$nFeature_RNA)
      )
    }
  )

  return(
    list(
      samples_list = samples_list,
      qc_summary = qc_summary
    )
  )
}



