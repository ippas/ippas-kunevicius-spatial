# ==============================================================================
# Add reconstructed Visium spatial images to individual Seurat objects
#
# Input layout:
# raw/<sample_ID>/<sample_ID>.png
# <path_to_data>/<sample_ID>/outs/spatial/tissue_positions.csv
#
# Original scalefactors_json.json is not required.
# A compatibility cache is created in results/spatial_image_cache/.
# ==============================================================================

add_spatial_images_to_samples <- function(
    samples_list,
    path_to_data,
    raw_image_root = "raw",
    cache_dir = "results/spatial_image_cache",
    assay = "RNA",
    spot_diameter_to_pitch = 0.55,
    verbose = TRUE
) {

  required_packages <- c(
    "png",
    "jsonlite",
    "RANN",
    "tibble",
    "dplyr"
  )

  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_packages) > 0) {
    stop(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", ")
    )
  }

  if (!is.list(samples_list) || length(samples_list) == 0) {
    stop("`samples_list` must be a non-empty named list of Seurat objects.")
  }

  sample_ids <- names(samples_list)

  if (is.null(sample_ids) || any(sample_ids == "")) {
    stop("`samples_list` must have sample IDs as list names.")
  }

  if (!all(vapply(samples_list, inherits, logical(1), what = "Seurat"))) {
    stop("Every element of `samples_list` must be a Seurat object.")
  }

  dir.create(
    cache_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  required_position_columns <- c(
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"
  )

  # ---------------------------------------------------------------------------
  # Create symlink where possible; otherwise copy the file.
  # ---------------------------------------------------------------------------
  link_or_copy_file <- function(source_file, target_file) {

    if (file.exists(target_file) || nzchar(Sys.readlink(target_file))) {
      unlink(target_file)
    }

    linked <- file.symlink(
      from = normalizePath(source_file, mustWork = TRUE),
      to = target_file
    )

    if (!linked) {
      copied <- file.copy(
        from = source_file,
        to = target_file,
        overwrite = TRUE
      )

      if (!copied) {
        stop(
          "Could not link or copy file:\n",
          source_file,
          "\n->\n",
          target_file
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Read current or legacy Space Ranger tissue-position file.
  # ---------------------------------------------------------------------------
  read_positions_file <- function(spatial_dir, sample_id) {

    current_file <- file.path(
      spatial_dir,
      "tissue_positions.csv"
    )

    legacy_file <- file.path(
      spatial_dir,
      "tissue_positions_list.csv"
    )

    if (file.exists(current_file)) {

      positions <- utils::read.csv(
        current_file,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      colnames(positions) <- tolower(colnames(positions))

    } else if (file.exists(legacy_file)) {

      positions <- utils::read.csv(
        legacy_file,
        header = FALSE,
        stringsAsFactors = FALSE,
        col.names = required_position_columns
      )

    } else {
      stop(
        "No tissue-position file found for sample `", sample_id, "`.\n",
        "Expected:\n",
        current_file,
        "\nor\n",
        legacy_file
      )
    }

    missing_columns <- setdiff(
      required_position_columns,
      colnames(positions)
    )

    if (length(missing_columns) > 0) {
      stop(
        "Tissue-position file for sample `", sample_id,
        "` is missing column(s): ",
        paste(missing_columns, collapse = ", ")
      )
    }

    positions <- positions[
      ,
      required_position_columns,
      drop = FALSE
    ]

    positions$barcode <- as.character(positions$barcode)

    if (anyDuplicated(positions$barcode) > 0) {
      stop(
        "Duplicated barcodes in tissue-position file for sample `",
        sample_id, "`."
      )
    }

    return(positions)
  }

  image_summary_list <- vector(
    mode = "list",
    length = length(sample_ids)
  )

  names(image_summary_list) <- sample_ids

  # ===========================================================================
  # Process one sample at a time
  # ===========================================================================

  for (sample_id in sample_ids) {

    if (verbose) {
      message("Adding spatial image: ", sample_id)
    }

    sample_object <- samples_list[[sample_id]]

    spatial_dir <- file.path(
      path_to_data,
      sample_id,
      "outs",
      "spatial"
    )

    if (!dir.exists(spatial_dir)) {
      stop(
        "Space Ranger spatial directory not found for sample `",
        sample_id, "`:\n",
        spatial_dir
      )
    }

    raw_png_file <- file.path(
      raw_image_root,
      sample_id,
      paste0(sample_id, ".png")
    )

    if (!file.exists(raw_png_file)) {
      stop(
        "Original PNG image not found for sample `",
        sample_id, "`:\n",
        raw_png_file
      )
    }

    # Prefer the hires Space Ranger PNG; use lowres only if hires is absent.
    spatial_image_candidates <- file.path(
      spatial_dir,
      c(
        "tissue_hires_image.png",
        "tissue_lowres_image.png"
      )
    )

    spatial_image_candidates <- spatial_image_candidates[
      file.exists(spatial_image_candidates)
    ]

    if (length(spatial_image_candidates) == 0) {
      stop(
        "No Space Ranger tissue PNG found for sample `",
        sample_id, "`.\nExpected tissue_hires_image.png or tissue_lowres_image.png."
      )
    }

    spatial_image_file <- spatial_image_candidates[1]

    positions <- read_positions_file(
      spatial_dir = spatial_dir,
      sample_id = sample_id
    )

    object_barcodes <- SeuratObject::Cells(sample_object)

    barcode_match <- match(
      object_barcodes,
      positions$barcode
    )

    if (anyNA(barcode_match)) {
      missing_barcodes <- object_barcodes[is.na(barcode_match)]

      stop(
        "Missing spatial coordinates for ",
        length(missing_barcodes),
        " barcode(s) in sample `", sample_id, "`.\nExamples: ",
        paste(head(missing_barcodes, 10), collapse = ", ")
      )
    }

    positions_for_object <- positions[
      barcode_match,
      ,
      drop = FALSE
    ]

    # -------------------------------------------------------------------------
    # Calculate scale from the original full-resolution PNG to Space Ranger PNG.
    # -------------------------------------------------------------------------
    raw_image <- png::readPNG(raw_png_file)
    spatial_image <- png::readPNG(spatial_image_file)

    raw_height <- nrow(raw_image)
    raw_width <- ncol(raw_image)

    spatial_height <- nrow(spatial_image)
    spatial_width <- ncol(spatial_image)

    scale_x <- spatial_width / raw_width
    scale_y <- spatial_height / raw_height

    relative_scale_difference <- abs(scale_x - scale_y) /
      mean(c(scale_x, scale_y))

    if (relative_scale_difference > 0.01) {
      stop(
        "Image aspect-ratio mismatch for sample `", sample_id, "`.\n",
        "Raw PNG: ", raw_width, " x ", raw_height, "\n",
        "Space Ranger PNG: ", spatial_width, " x ", spatial_height, "\n",
        "These files probably do not originate from the same source image."
      )
    }

    image_scale_factor <- mean(c(scale_x, scale_y))

    # -------------------------------------------------------------------------
    # Estimate Visium spot diameter from nearest-neighbour spot distance.
    # 55 um spot diameter / 100 um center-to-center pitch = 0.55.
    # -------------------------------------------------------------------------
    coordinate_matrix <- as.matrix(
      positions_for_object[
        ,
        c("pxl_col_in_fullres", "pxl_row_in_fullres"),
        drop = FALSE
      ]
    )

    nearest_neighbor_distances <- RANN::nn2(
      data = coordinate_matrix,
      query = coordinate_matrix,
      k = 2
    )$nn.dists[, 2]

    nearest_neighbor_distances <- nearest_neighbor_distances[
      is.finite(nearest_neighbor_distances) &
        nearest_neighbor_distances > 0
    ]

    if (length(nearest_neighbor_distances) == 0) {
      stop(
        "Could not estimate spot pitch for sample `",
        sample_id, "`."
      )
    }

    median_spot_pitch_fullres <- median(nearest_neighbor_distances)

    spot_diameter_fullres <- median_spot_pitch_fullres *
      spot_diameter_to_pitch

    # -------------------------------------------------------------------------
    # Create a Seurat-compatible cache folder.
    # -------------------------------------------------------------------------
    sample_cache_dir <- file.path(
      cache_dir,
      sample_id
    )

    dir.create(
      sample_cache_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    # Both names point to the selected Space Ranger PNG.
    # This makes Seurat's default lowres plotting work consistently.
    link_or_copy_file(
      source_file = spatial_image_file,
      target_file = file.path(
        sample_cache_dir,
        "tissue_hires_image.png"
      )
    )

    link_or_copy_file(
      source_file = spatial_image_file,
      target_file = file.path(
        sample_cache_dir,
        "tissue_lowres_image.png"
      )
    )

    # Read10X_Image expects this legacy no-header format.
    utils::write.table(
      positions,
      file = file.path(
        sample_cache_dir,
        "tissue_positions_list.csv"
      ),
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )

    reconstructed_scale_factors <- list(
      fiducial_diameter_fullres = spot_diameter_fullres,
      spot_diameter_fullres = spot_diameter_fullres,
      tissue_hires_scalef = image_scale_factor,
      tissue_lowres_scalef = image_scale_factor
    )

    jsonlite::write_json(
      x = reconstructed_scale_factors,
      path = file.path(
        sample_cache_dir,
        "scalefactors_json.json"
      ),
      auto_unbox = TRUE,
      pretty = TRUE
    )

    # -------------------------------------------------------------------------
    # Read reconstructed image and attach it to the Seurat object.
    # -------------------------------------------------------------------------
    image_name <- paste0("image_", sample_id)

    spatial_image_object <- Seurat::Read10X_Image(
      image.dir = sample_cache_dir,
      image.name = "tissue_lowres_image.png",
      assay = assay,
      slice = image_name,
      filter.matrix = FALSE,
      image.type = "VisiumV2"
    )

    spatial_image_object <- spatial_image_object[
      SeuratObject::Cells(sample_object)
    ]

    sample_object[[image_name]] <- spatial_image_object

    samples_list[[sample_id]] <- sample_object

    image_summary_list[[sample_id]] <- tibble::tibble(
      sample_ID = sample_id,
      raw_png_file = raw_png_file,
      spatial_image_file = spatial_image_file,
      raw_width = raw_width,
      raw_height = raw_height,
      spatial_width = spatial_width,
      spatial_height = spatial_height,
      image_scale_factor = image_scale_factor,
      median_spot_pitch_fullres = median_spot_pitch_fullres,
      spot_diameter_fullres = spot_diameter_fullres,
      image_name = image_name
    )
  }

  image_summary <- dplyr::bind_rows(image_summary_list)

  return(
    list(
      samples_list = samples_list,
      image_summary = image_summary
    )
  )
}