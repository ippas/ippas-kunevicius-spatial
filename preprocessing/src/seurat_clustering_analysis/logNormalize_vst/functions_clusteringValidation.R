# ==============================================================================
# functions_clusteringValidation.R
#
# Purpose:
# Validate multiple Seurat clustering solutions generated from the same graph.
#
# Metrics:
# 1. Transcriptomic silhouette in an integrated embedding.
# 2. Spatial silhouette calculated separately for each tissue section.
# 3. CHAOS calculated separately for each tissue section.
# 4. PAS calculated separately for each tissue section.
# 5. ARI stability between clustering solutions.
# 6. Cluster-size diagnostics.
#
# Main output:
# - one validation-summary TSV for every clustering algorithm;
# - one spatial-by-sample TSV for every clustering algorithm;
# - one ARI TSV for every clustering algorithm;
# - one combined validation-summary TSV;
# - one combined spatial-by-sample TSV;
# - one pairwise ARI table for all clustering solutions;
# - one parameter table;
# - PNG and PDF metric plots for every algorithm;
# - one combined PNG and PDF metric plot;
# - error bars only for spatial ASW, CHAOS and PAS on separate
#   algorithm plots; no error bars on the combined all-methods plot.
#
# Important:
# - transcriptomic silhouette is calculated on one reproducible common sample
#   of spots shared by all clustering solutions;
# - spatial silhouette is calculated per tissue section;
# - CHAOS and PAS are calculated on all available spots in each section;
# - physical coordinates are normalized by the median nearest-neighbour spot
#   distance within each section before CHAOS is calculated;
# - lower CHAOS and PAS are better;
# - higher transcriptomic ASW, spatial ASW and ARI are better.
# ==============================================================================


# ==============================================================================
# 1. General utilities
# ==============================================================================

clustering_validation_message <- function(...) {

  cat(
    "[",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    "] ",
    paste0(..., collapse = ""),
    "\n",
    sep = ""
  )

  flush.console()
  invisible(NULL)
}


format_validation_elapsed_time <- function(seconds) {

  seconds <- max(
    0,
    as.numeric(seconds)
  )

  hours <- floor(
    seconds / 3600
  )

  minutes <- floor(
    (seconds %% 3600) / 60
  )

  remaining_seconds <- floor(
    seconds %% 60
  )

  sprintf(
    "%02d:%02d:%02d",
    as.integer(hours),
    as.integer(minutes),
    as.integer(remaining_seconds)
  )
}


write_validation_tsv <- function(
    data,
    filename
) {

  dir.create(
    dirname(filename),
    recursive = TRUE,
    showWarnings = FALSE
  )

  utils::write.table(
    x = data,
    file = filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )

  invisible(filename)
}


save_validation_plot <- function(
    plot_object,
    output_prefix,
    save_png = TRUE,
    save_pdf = TRUE,
    png_width_in = 14,
    png_height_in = 12,
    pdf_width_in = 14,
    pdf_height_in = 12,
    dpi = 300
) {

  output_files <- list()

  dir.create(
    dirname(output_prefix),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (isTRUE(save_png)) {

    png_file <- paste0(
      output_prefix,
      ".png"
    )

    ggplot2::ggsave(
      filename = png_file,
      plot = plot_object,
      width = png_width_in,
      height = png_height_in,
      units = "in",
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )

    output_files$png <- png_file
  }

  if (isTRUE(save_pdf)) {

    pdf_file <- paste0(
      output_prefix,
      ".pdf"
    )

    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot_object,
      width = pdf_width_in,
      height = pdf_height_in,
      units = "in",
      bg = "white",
      limitsize = FALSE
    )

    output_files$pdf <- pdf_file
  }

  output_files
}


parse_clustering_validation_column <- function(
    cluster_column
) {

  if (
    !is.character(cluster_column) ||
      length(cluster_column) != 1L ||
      is.na(cluster_column) ||
      cluster_column == ""
  ) {
    stop(
      "`cluster_column` must be one non-empty character value.",
      call. = FALSE
    )
  }

  parsed <- regexec(
    pattern = "^(.*)_res([0-9]+(?:\\.[0-9]+)?)$",
    text = cluster_column
  )

  parsed_match <- regmatches(
    cluster_column,
    parsed
  )[[1]]

  if (length(parsed_match) != 3L) {
    stop(
      "Could not parse clustering column `",
      cluster_column,
      "`. Expected a name similar to `leiden_res0.20`.",
      call. = FALSE
    )
  }

  list(
    algorithm = parsed_match[[2]],
    resolution = as.numeric(
      parsed_match[[3]]
    )
  )
}


build_clustering_solution_table <- function(
    cluster_columns
) {

  if (
    !is.character(cluster_columns) ||
      length(cluster_columns) == 0L ||
      anyNA(cluster_columns) ||
      any(cluster_columns == "")
  ) {
    stop(
      "`cluster_columns` must be a non-empty character vector.",
      call. = FALSE
    )
  }

  if (anyDuplicated(cluster_columns) > 0L) {
    stop(
      "`cluster_columns` contains duplicated names.",
      call. = FALSE
    )
  }

  parsed_columns <- lapply(
    cluster_columns,
    parse_clustering_validation_column
  )

  solution_table <- data.frame(
    clusterColumn = cluster_columns,
    algorithm = vapply(
      parsed_columns,
      function(x) x$algorithm,
      character(1)
    ),
    resolution = vapply(
      parsed_columns,
      function(x) x$resolution,
      numeric(1)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  solution_table <- solution_table[
    order(
      solution_table$algorithm,
      solution_table$resolution
    ),
    ,
    drop = FALSE
  ]

  rownames(solution_table) <- NULL

  solution_table
}


build_validation_parameter_subtitle <- function(
    integration_method,
    normalization_label,
    n_hvg,
    reduction,
    dims,
    k_param,
    prune_snn,
    transcriptomic_sample_n,
    spatial_asw_max_spots_per_sample,
    pas_k,
    pas_min_different,
    ari_comparison_mode,
    seed
) {

  line_1 <- paste0(
    "Integration: ",
    integration_method,
    " | normalization: ",
    normalization_label,
    " | HVG: ",
    n_hvg,
    " | reduction: ",
    reduction,
    " | dims: ",
    min(dims),
    "–",
    max(dims),
    " | k.param: ",
    k_param,
    " | prune.SNN: ",
    formatC(
      prune_snn,
      format = "f",
      digits = 4L
    )
  )

  line_2 <- paste0(
    "Transcriptomic ASW common sample: n=",
    format(
      transcriptomic_sample_n,
      big.mark = " ",
      scientific = FALSE
    ),
    " | spatial ASW max spots/sample: ",
    spatial_asw_max_spots_per_sample,
    " | PAS: k=",
    pas_k,
    ", threshold=",
    pas_min_different,
    " | ARI: ",
    ari_comparison_mode,
    " | seed=",
    seed
  )

  paste(
    line_1,
    line_2,
    sep = "\n"
  )
}


# ==============================================================================
# 2. Input validation
# ==============================================================================

validate_clustering_validation_inputs <- function(
    seurat_object,
    cluster_columns,
    reduction,
    dims,
    sample_order
) {

  if (!inherits(seurat_object, "Seurat")) {
    stop(
      "`seurat_object` must be a Seurat object.",
      call. = FALSE
    )
  }

  metadata_table <- seurat_object[[]]

  missing_cluster_columns <- setdiff(
    cluster_columns,
    colnames(metadata_table)
  )

  if (length(missing_cluster_columns) > 0L) {
    stop(
      "The following clustering columns are absent from the Seurat object:\n",
      paste(
        missing_cluster_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  for (cluster_column in cluster_columns) {

    cluster_values <- metadata_table[[cluster_column]]

    if (anyNA(cluster_values)) {
      stop(
        "NA values were found in clustering column `",
        cluster_column,
        "`.",
        call. = FALSE
      )
    }

    if (length(unique(cluster_values)) < 2L) {
      stop(
        "Clustering column `",
        cluster_column,
        "` contains fewer than two clusters.",
        call. = FALSE
      )
    }
  }

  if (!"sample_ID" %in% colnames(metadata_table)) {
    stop(
      "The Seurat object does not contain `sample_ID` metadata.",
      call. = FALSE
    )
  }

  if (
    !is.character(sample_order) ||
      length(sample_order) == 0L ||
      anyNA(sample_order) ||
      any(sample_order == "")
  ) {
    stop(
      "`sample_order` must be a non-empty character vector.",
      call. = FALSE
    )
  }

  if (anyDuplicated(sample_order) > 0L) {
    stop(
      "`sample_order` contains duplicated sample IDs.",
      call. = FALSE
    )
  }

  missing_samples <- setdiff(
    sample_order,
    unique(
      as.character(metadata_table$sample_ID)
    )
  )

  if (length(missing_samples) > 0L) {
    stop(
      "The following requested sample IDs are absent from metadata: ",
      paste(
        missing_samples,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  available_reductions <- SeuratObject::Reductions(
    seurat_object
  )

  if (!reduction %in% available_reductions) {
    stop(
      "Reduction `",
      reduction,
      "` is absent from the Seurat object.\nAvailable reductions: ",
      paste(
        available_reductions,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  if (
    !is.numeric(dims) ||
      length(dims) == 0L ||
      anyNA(dims) ||
      any(dims < 1)
  ) {
    stop(
      "`dims` must be a non-empty vector of positive dimensions.",
      call. = FALSE
    )
  }

  available_dims <- ncol(
    SeuratObject::Embeddings(
      seurat_object[[reduction]]
    )
  )

  if (max(dims) > available_dims) {
    stop(
      "Requested dimension ",
      max(dims),
      " exceeds the ",
      available_dims,
      " dimensions available in reduction `",
      reduction,
      "`.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# ==============================================================================
# 3. Common spot sampling for transcriptomic silhouette
# ==============================================================================

select_common_transcriptomic_validation_spots <- function(
    metadata_table,
    cluster_columns,
    max_spots = 4000L,
    minimum_spots_per_cluster = 20L,
    seed = 7L,
    verbose = TRUE
) {

  if (is.null(rownames(metadata_table))) {
    stop(
      "`metadata_table` must have spot names as row names.",
      call. = FALSE
    )
  }

  all_spots <- rownames(
    metadata_table
  )

  if (
    !is.numeric(max_spots) ||
      length(max_spots) != 1L ||
      is.na(max_spots) ||
      max_spots < 2L
  ) {
    stop(
      "`max_spots` must be one integer greater than or equal to 2.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(minimum_spots_per_cluster) ||
      length(minimum_spots_per_cluster) != 1L ||
      is.na(minimum_spots_per_cluster) ||
      minimum_spots_per_cluster < 1L
  ) {
    stop(
      "`minimum_spots_per_cluster` must be one positive integer.",
      call. = FALSE
    )
  }

  set.seed(
    as.integer(seed)
  )

  mandatory_spots <- character(0)

  for (cluster_column in cluster_columns) {

    cluster_values <- as.character(
      metadata_table[[cluster_column]]
    )

    cluster_levels <- unique(
      cluster_values
    )

    for (cluster_level in cluster_levels) {

      available_spots <- all_spots[
        cluster_values == cluster_level
      ]

      number_to_select <- min(
        length(available_spots),
        as.integer(minimum_spots_per_cluster)
      )

      mandatory_spots <- union(
        mandatory_spots,
        sample(
          available_spots,
          size = number_to_select,
          replace = FALSE
        )
      )
    }
  }

  remaining_spots <- setdiff(
    all_spots,
    mandatory_spots
  )

  target_number <- min(
    length(all_spots),
    as.integer(max_spots)
  )

  if (length(mandatory_spots) < target_number) {

    number_to_add <- min(
      target_number - length(mandatory_spots),
      length(remaining_spots)
    )

    additional_spots <- if (number_to_add > 0L) {
      sample(
        remaining_spots,
        size = number_to_add,
        replace = FALSE
      )
    } else {
      character(0)
    }

    selected_spots <- c(
      mandatory_spots,
      additional_spots
    )

  } else {

    # With many clustering solutions, the union of per-cluster mandatory spots
    # can become very large. Never exceed `max_spots`, because an exact
    # silhouette uses a quadratic distance matrix.
    selected_spots <- sample(
      mandatory_spots,
      size = target_number,
      replace = FALSE
    )

    if (
      verbose &&
        length(mandatory_spots) > target_number
    ) {
      clustering_validation_message(
        "The mandatory silhouette candidate pool contains ",
        length(mandatory_spots),
        " spots. It was reproducibly capped at `max_spots = ",
        target_number,
        "`."
      )
    }
  }

  selected_spots <- unique(
    selected_spots
  )

  if (verbose) {
    clustering_validation_message(
      "Selected common transcriptomic silhouette sample: ",
      length(selected_spots),
      " of ",
      length(all_spots),
      " spots."
    )
  }

  selected_spots
}


# ==============================================================================
# 4. Spatial coordinates
# ==============================================================================

create_validation_image_map <- function(
    seurat_object,
    sample_order
) {

  image_names <- SeuratObject::Images(
    seurat_object
  )

  if (length(image_names) == 0L) {
    stop(
      "No spatial images were found in the Seurat object.",
      call. = FALSE
    )
  }

  metadata_table <- seurat_object[[]]

  mapping_rows <- lapply(
    image_names,
    function(image_name) {

      image_cells <- SeuratObject::Cells(
        seurat_object[[image_name]]
      )

      if (length(image_cells) == 0L) {
        stop(
          "Spatial image `",
          image_name,
          "` contains zero spots.",
          call. = FALSE
        )
      }

      sample_ids <- unique(
        as.character(
          metadata_table[
            image_cells,
            "sample_ID",
            drop = TRUE
          ]
        )
      )

      if (length(sample_ids) != 1L) {
        stop(
          "Spatial image `",
          image_name,
          "` maps to ",
          length(sample_ids),
          " sample IDs. Exactly one was expected.",
          call. = FALSE
        )
      }

      data.frame(
        sample_ID = sample_ids[[1]],
        imageName = image_name,
        nSpots = length(image_cells),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )

  image_map <- do.call(
    rbind,
    mapping_rows
  )

  rownames(image_map) <- NULL

  if (anyDuplicated(image_map$sample_ID) > 0L) {
    stop(
      "More than one image was mapped to the same sample ID.",
      call. = FALSE
    )
  }

  missing_images <- setdiff(
    sample_order,
    image_map$sample_ID
  )

  if (length(missing_images) > 0L) {
    stop(
      "The following requested sample IDs do not have spatial images: ",
      paste(
        missing_images,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  image_map <- image_map[
    match(
      sample_order,
      image_map$sample_ID
    ),
    ,
    drop = FALSE
  ]

  rownames(image_map) <- NULL

  image_map
}


extract_spatial_coordinates_from_image <- function(
    seurat_object,
    image_name,
    sample_id
) {

  image_object <- seurat_object[[image_name]]

  image_cells <- SeuratObject::Cells(
    image_object
  )

  coordinate_table <- tryCatch(
    {
      SeuratObject::GetTissueCoordinates(
        object = image_object,
        scale = NULL
      )
    },
    error = function(first_error) {
      SeuratObject::GetTissueCoordinates(
        object = image_object
      )
    }
  )

  coordinate_table <- as.data.frame(
    coordinate_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (nrow(coordinate_table) == 0L) {
    stop(
      "No tissue coordinates were returned for image `",
      image_name,
      "`.",
      call. = FALSE
    )
  }

  coordinate_column_names <- colnames(
    coordinate_table
  )

  spot_column_candidates <- c(
    "cell",
    "Cell",
    "barcode",
    "Barcode",
    "spot",
    "Spot",
    "ID"
  )

  spot_column <- spot_column_candidates[
    spot_column_candidates %in%
      coordinate_column_names
  ]

  if (length(spot_column) > 0L) {

    spot_names <- as.character(
      coordinate_table[[spot_column[[1]]]]
    )

  } else if (
    !is.null(rownames(coordinate_table)) &&
      !identical(
        rownames(coordinate_table),
        as.character(
          seq_len(nrow(coordinate_table))
        )
      )
  ) {

    spot_names <- rownames(
      coordinate_table
    )

  } else if (
    nrow(coordinate_table) ==
      length(image_cells)
  ) {

    spot_names <- image_cells

  } else {

    stop(
      "Could not identify spot names in coordinates for image `",
      image_name,
      "`.",
      call. = FALSE
    )
  }

  coordinate_pairs <- list(
    c("x", "y"),
    c("imagecol", "imagerow"),
    c(
      "pxl_col_in_fullres",
      "pxl_row_in_fullres"
    )
  )

  selected_coordinate_pair <- NULL

  for (coordinate_pair in coordinate_pairs) {

    if (
      all(
        coordinate_pair %in%
          coordinate_column_names
      )
    ) {
      selected_coordinate_pair <-
        coordinate_pair

      break
    }
  }

  if (is.null(selected_coordinate_pair)) {
    stop(
      "Could not identify x/y coordinate columns for image `",
      image_name,
      "`. Available columns: ",
      paste(
        coordinate_column_names,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  normalized_coordinates <- data.frame(
    spot = spot_names,
    sample_ID = sample_id,
    x = as.numeric(
      coordinate_table[[selected_coordinate_pair[[1]]]]
    ),
    y = as.numeric(
      coordinate_table[[selected_coordinate_pair[[2]]]]
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (
    anyNA(normalized_coordinates$x) ||
      anyNA(normalized_coordinates$y)
  ) {
    stop(
      "NA values were found in spatial coordinates for sample `",
      sample_id,
      "`.",
      call. = FALSE
    )
  }

  if (anyDuplicated(normalized_coordinates$spot) > 0L) {
    stop(
      "Duplicated spot names were found in spatial coordinates for sample `",
      sample_id,
      "`.",
      call. = FALSE
    )
  }

  missing_image_cells <- setdiff(
    image_cells,
    normalized_coordinates$spot
  )

  if (length(missing_image_cells) > 0L) {
    stop(
      "Spatial coordinates are missing for ",
      length(missing_image_cells),
      " spots in sample `",
      sample_id,
      "`.",
      call. = FALSE
    )
  }

  normalized_coordinates <- normalized_coordinates[
    match(
      image_cells,
      normalized_coordinates$spot
    ),
    ,
    drop = FALSE
  ]

  rownames(normalized_coordinates) <-
    normalized_coordinates$spot

  normalized_coordinates
}


extract_all_spatial_coordinates <- function(
    seurat_object,
    sample_order,
    verbose = TRUE
) {

  image_map <- create_validation_image_map(
    seurat_object = seurat_object,
    sample_order = sample_order
  )

  coordinate_list <- vector(
    mode = "list",
    length = nrow(image_map)
  )

  names(coordinate_list) <-
    image_map$sample_ID

  for (sample_index in seq_len(nrow(image_map))) {

    sample_id <- image_map$sample_ID[[sample_index]]

    image_name <- image_map$imageName[[sample_index]]

    if (verbose) {
      clustering_validation_message(
        "Extracting spatial coordinates: ",
        sample_id,
        " (",
        image_name,
        ")."
      )
    }

    coordinate_list[[sample_id]] <-
      extract_spatial_coordinates_from_image(
        seurat_object = seurat_object,
        image_name = image_name,
        sample_id = sample_id
      )
  }

  list(
    imageMap = image_map,
    coordinatesBySample = coordinate_list
  )
}


estimate_spatial_spot_pitch <- function(
    coordinate_matrix
) {

  if (nrow(coordinate_matrix) < 2L) {
    stop(
      "At least two spots are required to estimate spatial spot pitch.",
      call. = FALSE
    )
  }

  nearest_neighbour_result <- RANN::nn2(
    data = coordinate_matrix,
    query = coordinate_matrix,
    k = 2L
  )

  nearest_neighbour_distances <-
    nearest_neighbour_result$nn.dists[
      ,
      2L
    ]

  nearest_neighbour_distances <-
    nearest_neighbour_distances[
      is.finite(
        nearest_neighbour_distances
      ) &
        nearest_neighbour_distances > 0
    ]

  if (
    length(nearest_neighbour_distances) ==
      0L
  ) {
    stop(
      "Could not estimate nearest-neighbour spot distance.",
      call. = FALSE
    )
  }

  stats::median(
    nearest_neighbour_distances
  )
}


# ==============================================================================
# 5. Silhouette metrics
# ==============================================================================

calculate_silhouette_summary_from_distance <- function(
    distance_object,
    cluster_labels
) {

  cluster_labels <- as.character(
    cluster_labels
  )

  n_observations <- attr(
    distance_object,
    "Size"
  )

  if (is.null(n_observations)) {
    stop(
      "`distance_object` must be a valid `dist` object.",
      call. = FALSE
    )
  }

  if (
    n_observations !=
      length(cluster_labels)
  ) {
    stop(
      "The size of `distance_object` does not match cluster labels.",
      call. = FALSE
    )
  }

  if (n_observations < 3L) {
    return(
      data.frame(
        mean = NA_real_,
        median = NA_real_,
        sd = NA_real_,
        fractionNegative = NA_real_,
        nEvaluated = n_observations,
        nClustersEvaluated =
          length(unique(cluster_labels)),
        stringsAsFactors = FALSE
      )
    )
  }

  factor_labels <- factor(
    cluster_labels
  )

  if (nlevels(factor_labels) < 2L) {
    return(
      data.frame(
        mean = NA_real_,
        median = NA_real_,
        sd = NA_real_,
        fractionNegative = NA_real_,
        nEvaluated = n_observations,
        nClustersEvaluated =
          nlevels(factor_labels),
        stringsAsFactors = FALSE
      )
    )
  }

  silhouette_object <- cluster::silhouette(
    x = as.integer(
      factor_labels
    ),
    dist = distance_object
  )

  silhouette_values <- as.numeric(
    silhouette_object[
      ,
      "sil_width"
    ]
  )

  data.frame(
    mean = mean(
      silhouette_values,
      na.rm = TRUE
    ),
    median = stats::median(
      silhouette_values,
      na.rm = TRUE
    ),
    sd = stats::sd(
      silhouette_values,
      na.rm = TRUE
    ),
    fractionNegative = mean(
      silhouette_values < 0,
      na.rm = TRUE
    ),
    nEvaluated = length(
      silhouette_values
    ),
    nClustersEvaluated = nlevels(
      factor_labels
    ),
    stringsAsFactors = FALSE
  )
}


calculate_silhouette_summary <- function(
    data_matrix,
    cluster_labels
) {

  data_matrix <- as.matrix(
    data_matrix
  )

  cluster_labels <- as.character(
    cluster_labels
  )

  if (
    nrow(data_matrix) !=
      length(cluster_labels)
  ) {
    stop(
      "The number of rows in `data_matrix` does not match cluster labels.",
      call. = FALSE
    )
  }

  distance_object <- stats::dist(
    data_matrix,
    method = "euclidean"
  )

  calculate_silhouette_summary_from_distance(
    distance_object = distance_object,
    cluster_labels = cluster_labels
  )
}


calculate_transcriptomic_silhouette <- function(
    seurat_object,
    cluster_column,
    reduction,
    dims,
    selected_spots,
    distance_object = NULL
) {

  cluster_labels <- seurat_object[[]][
    selected_spots,
    cluster_column,
    drop = TRUE
  ]

  if (is.null(distance_object)) {

    embedding_matrix <- SeuratObject::Embeddings(
      seurat_object[[reduction]]
    )

    embedding_subset <- embedding_matrix[
      selected_spots,
      dims,
      drop = FALSE
    ]

    distance_object <- stats::dist(
      embedding_subset,
      method = "euclidean"
    )
  }

  calculate_silhouette_summary_from_distance(
    distance_object = distance_object,
    cluster_labels = cluster_labels
  )
}


select_common_spatial_silhouette_spots <- function(
    sample_spots,
    metadata_table,
    cluster_columns,
    max_spots,
    minimum_spots_per_cluster = 10L,
    seed = 7L
) {

  if (length(sample_spots) <= max_spots) {
    return(sample_spots)
  }

  sample_metadata <- metadata_table[
    sample_spots,
    ,
    drop = FALSE
  ]

  select_common_transcriptomic_validation_spots(
    metadata_table = sample_metadata,
    cluster_columns = cluster_columns,
    max_spots = max_spots,
    minimum_spots_per_cluster =
      minimum_spots_per_cluster,
    seed = seed,
    verbose = FALSE
  )
}


calculate_spatial_silhouette <- function(
    coordinate_table,
    cluster_labels,
    selected_spots,
    distance_object = NULL
) {

  selected_cluster_labels <- cluster_labels[
    selected_spots
  ]

  if (is.null(distance_object)) {

    selected_coordinates <- coordinate_table[
      selected_spots,
      c("x", "y"),
      drop = FALSE
    ]

    distance_object <- stats::dist(
      selected_coordinates,
      method = "euclidean"
    )
  }

  calculate_silhouette_summary_from_distance(
    distance_object = distance_object,
    cluster_labels = selected_cluster_labels
  )
}


# ==============================================================================
# 6. CHAOS and PAS
# ==============================================================================

calculate_chaos_score <- function(
    coordinate_matrix,
    cluster_labels,
    normalize_by_spot_pitch = TRUE
) {

  coordinate_matrix <- as.matrix(
    coordinate_matrix
  )

  cluster_labels <- as.character(
    cluster_labels
  )

  if (
    nrow(coordinate_matrix) !=
      length(cluster_labels)
  ) {
    stop(
      "Coordinate and cluster-label lengths differ in CHAOS calculation.",
      call. = FALSE
    )
  }

  spot_pitch <- estimate_spatial_spot_pitch(
    coordinate_matrix
  )

  if (isTRUE(normalize_by_spot_pitch)) {
    coordinate_matrix <-
      coordinate_matrix / spot_pitch
  }

  cluster_levels <- unique(
    cluster_labels
  )

  nearest_same_cluster_distances <-
    rep(
      NA_real_,
      nrow(coordinate_matrix)
    )

  for (cluster_level in cluster_levels) {

    cluster_indices <- which(
      cluster_labels == cluster_level
    )

    if (length(cluster_indices) < 2L) {
      next
    }

    cluster_coordinates <-
      coordinate_matrix[
        cluster_indices,
        ,
        drop = FALSE
      ]

    nearest_result <- RANN::nn2(
      data = cluster_coordinates,
      query = cluster_coordinates,
      k = 2L
    )

    nearest_same_cluster_distances[
      cluster_indices
    ] <- nearest_result$nn.dists[
      ,
      2L
    ]
  }

  valid_distances <- is.finite(
    nearest_same_cluster_distances
  )

  data.frame(
    CHAOS = if (any(valid_distances)) {
      mean(
        nearest_same_cluster_distances[
          valid_distances
        ]
      )
    } else {
      NA_real_
    },
    spotPitch = spot_pitch,
    nValidSpots = sum(
      valid_distances
    ),
    nExcludedSingletonClusterSpots =
      sum(!valid_distances),
    normalizedBySpotPitch =
      isTRUE(normalize_by_spot_pitch),
    stringsAsFactors = FALSE
  )
}


calculate_pas_score <- function(
    coordinate_matrix,
    cluster_labels,
    k = 10L,
    minimum_different_neighbours = 6L
) {

  coordinate_matrix <- as.matrix(
    coordinate_matrix
  )

  cluster_labels <- as.character(
    cluster_labels
  )

  if (
    nrow(coordinate_matrix) !=
      length(cluster_labels)
  ) {
    stop(
      "Coordinate and cluster-label lengths differ in PAS calculation.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(k) ||
      length(k) != 1L ||
      is.na(k) ||
      k < 1L
  ) {
    stop(
      "`k` must be one positive integer.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(minimum_different_neighbours) ||
      length(minimum_different_neighbours) != 1L ||
      is.na(minimum_different_neighbours) ||
      minimum_different_neighbours < 1L
  ) {
    stop(
      "`minimum_different_neighbours` must be one positive integer.",
      call. = FALSE
    )
  }

  if (minimum_different_neighbours > k) {
    stop(
      "`minimum_different_neighbours` cannot exceed `k`.",
      call. = FALSE
    )
  }

  if (nrow(coordinate_matrix) <= k) {
    stop(
      "PAS requires more than ",
      k,
      " spots.",
      call. = FALSE
    )
  }

  nearest_result <- RANN::nn2(
    data = coordinate_matrix,
    query = coordinate_matrix,
    k = as.integer(k) + 1L
  )

  neighbour_indices <- nearest_result$nn.idx[
    ,
    -1L,
    drop = FALSE
  ]

  different_neighbour_counts <- vapply(
    seq_len(nrow(neighbour_indices)),
    function(spot_index) {

      sum(
        cluster_labels[
          neighbour_indices[
            spot_index,
            ,
            drop = TRUE
          ]
        ] !=
          cluster_labels[[spot_index]]
      )
    },
    integer(1)
  )

  abnormal_spots <-
    different_neighbour_counts >=
    as.integer(
      minimum_different_neighbours
    )

  pas_fraction <- mean(
    abnormal_spots
  )

  data.frame(
    PAS_fraction = pas_fraction,
    PAS_percent = 100 * pas_fraction,
    nAbnormalSpots = sum(abnormal_spots),
    nSpots = length(abnormal_spots),
    k = as.integer(k),
    minimumDifferentNeighbours =
      as.integer(
        minimum_different_neighbours
      ),
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 7. ARI stability
# ==============================================================================

calculate_adjusted_rand_index <- function(
    labels_1,
    labels_2
) {

  labels_1 <- as.character(
    labels_1
  )

  labels_2 <- as.character(
    labels_2
  )

  if (length(labels_1) != length(labels_2)) {
    stop(
      "ARI label vectors must have equal lengths.",
      call. = FALSE
    )
  }

  contingency_table <- table(
    labels_1,
    labels_2
  )

  choose_two <- function(x) {
    x * (x - 1) / 2
  }

  observed_index <- sum(
    choose_two(
      contingency_table
    )
  )

  row_index <- sum(
    choose_two(
      rowSums(contingency_table)
    )
  )

  column_index <- sum(
    choose_two(
      colSums(contingency_table)
    )
  )

  total_pairs <- choose_two(
    sum(contingency_table)
  )

  if (total_pairs == 0) {
    return(NA_real_)
  }

  expected_index <-
    row_index * column_index /
    total_pairs

  maximum_index <-
    0.5 * (
      row_index +
        column_index
    )

  denominator <-
    maximum_index -
    expected_index

  if (
    abs(denominator) <
      .Machine$double.eps
  ) {

    partitions_identical <-
      all(
        rowSums(contingency_table > 0) == 1L
      ) &&
      all(
        colSums(contingency_table > 0) == 1L
      )

    return(
      if (partitions_identical) {
        1
      } else {
        0
      }
    )
  }

  (
    observed_index -
      expected_index
  ) / denominator
}


calculate_pairwise_ari_table <- function(
    metadata_table,
    solution_table,
    comparison_mode = c(
      "all",
      "withinAlgorithmAndSameResolution",
      "adjacentAndSameResolution"
    )
) {

  comparison_mode <- match.arg(
    comparison_mode
  )

  cluster_columns <-
    solution_table$clusterColumn

  if (length(cluster_columns) < 2L) {
    return(
      data.frame(
        clusterColumn1 = character(0),
        algorithm1 = character(0),
        resolution1 = numeric(0),
        clusterColumn2 = character(0),
        algorithm2 = character(0),
        resolution2 = numeric(0),
        ARI = numeric(0),
        sameAlgorithm = logical(0),
        sameResolution = logical(0),
        stringsAsFactors = FALSE
      )
    )
  }

  pairs <- utils::combn(
    cluster_columns,
    m = 2L,
    simplify = FALSE
  )

  keep_pair <- vapply(
    pairs,
    function(current_pair) {

      solution_1 <- solution_table[
        solution_table$clusterColumn ==
          current_pair[[1]],
        ,
        drop = FALSE
      ]

      solution_2 <- solution_table[
        solution_table$clusterColumn ==
          current_pair[[2]],
        ,
        drop = FALSE
      ]

      same_algorithm <-
        solution_1$algorithm ==
        solution_2$algorithm

      same_resolution <-
        isTRUE(
          all.equal(
            solution_1$resolution,
            solution_2$resolution
          )
        )

      if (comparison_mode == "all") {
        return(TRUE)
      }

      if (
        comparison_mode ==
          "withinAlgorithmAndSameResolution"
      ) {
        return(
          same_algorithm ||
            same_resolution
        )
      }

      resolution_difference <- abs(
        solution_1$resolution -
          solution_2$resolution
      )

      algorithm_resolutions <- sort(
        unique(
          solution_table$resolution[
            solution_table$algorithm ==
              solution_1$algorithm
          ]
        )
      )

      adjacent_step <- if (
        length(algorithm_resolutions) > 1L
      ) {
        min(
          diff(algorithm_resolutions)
        )
      } else {
        Inf
      }

      adjacent_same_algorithm <-
        same_algorithm &&
        isTRUE(
          all.equal(
            resolution_difference,
            adjacent_step,
            tolerance = 1e-10
          )
        )

      adjacent_same_algorithm ||
        same_resolution
    },
    logical(1)
  )

  pairs <- pairs[
    keep_pair
  ]

  ari_rows <- lapply(
    pairs,
    function(current_pair) {

      column_1 <- current_pair[[1]]
      column_2 <- current_pair[[2]]

      solution_1 <- solution_table[
        solution_table$clusterColumn ==
          column_1,
        ,
        drop = FALSE
      ]

      solution_2 <- solution_table[
        solution_table$clusterColumn ==
          column_2,
        ,
        drop = FALSE
      ]

      data.frame(
        clusterColumn1 = column_1,
        algorithm1 = solution_1$algorithm,
        resolution1 = solution_1$resolution,
        clusterColumn2 = column_2,
        algorithm2 = solution_2$algorithm,
        resolution2 = solution_2$resolution,
        ARI = calculate_adjusted_rand_index(
          labels_1 =
            metadata_table[[column_1]],
          labels_2 =
            metadata_table[[column_2]]
        ),
        sameAlgorithm =
          solution_1$algorithm ==
          solution_2$algorithm,
        sameResolution =
          isTRUE(
            all.equal(
              solution_1$resolution,
              solution_2$resolution
            )
          ),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )

  if (length(ari_rows) == 0L) {
    return(
      data.frame(
        clusterColumn1 = character(0),
        algorithm1 = character(0),
        resolution1 = numeric(0),
        clusterColumn2 = character(0),
        algorithm2 = character(0),
        resolution2 = numeric(0),
        ARI = numeric(0),
        sameAlgorithm = logical(0),
        sameResolution = logical(0),
        stringsAsFactors = FALSE
      )
    )
  }

  ari_table <- do.call(
    rbind,
    ari_rows
  )

  rownames(ari_table) <- NULL

  ari_table
}


add_adjacent_resolution_ari <- function(
    validation_summary,
    pairwise_ari
) {

  validation_summary$ARI_vsPreviousResolution <-
    NA_real_

  validation_summary$ARI_vsNextResolution <-
    NA_real_

  algorithms <- unique(
    validation_summary$algorithm
  )

  for (algorithm_name in algorithms) {

    algorithm_indices <- which(
      validation_summary$algorithm ==
        algorithm_name
    )

    algorithm_indices <- algorithm_indices[
      order(
        validation_summary$resolution[
          algorithm_indices
        ]
      )
    ]

    if (length(algorithm_indices) < 2L) {
      next
    }

    for (
      local_index in
      seq_along(algorithm_indices)
    ) {

      current_index <-
        algorithm_indices[[local_index]]

      current_column <-
        validation_summary$clusterColumn[[current_index]]

      if (local_index > 1L) {

        previous_index <-
          algorithm_indices[[local_index - 1L]]

        previous_column <-
          validation_summary$clusterColumn[[previous_index]]

        ari_row <- pairwise_ari[
          (
            pairwise_ari$clusterColumn1 ==
              previous_column &
              pairwise_ari$clusterColumn2 ==
              current_column
          ) |
            (
              pairwise_ari$clusterColumn1 ==
                current_column &
                pairwise_ari$clusterColumn2 ==
                previous_column
            ),
          ,
          drop = FALSE
        ]

        if (nrow(ari_row) == 1L) {
          validation_summary$ARI_vsPreviousResolution[[current_index]] <- ari_row$ARI[[1]]
        }
      }

      if (
        local_index <
          length(algorithm_indices)
      ) {

        next_index <-
          algorithm_indices[[local_index + 1L]]

        next_column <-
          validation_summary$clusterColumn[[next_index]]

        ari_row <- pairwise_ari[
          (
            pairwise_ari$clusterColumn1 ==
              current_column &
              pairwise_ari$clusterColumn2 ==
              next_column
          ) |
            (
              pairwise_ari$clusterColumn1 ==
                next_column &
                pairwise_ari$clusterColumn2 ==
                current_column
            ),
          ,
          drop = FALSE
        ]

        if (nrow(ari_row) == 1L) {
          validation_summary$ARI_vsNextResolution[[current_index]] <- ari_row$ARI[[1]]
        }
      }
    }
  }

  validation_summary
}


# ==============================================================================
# 8. Cluster-size diagnostics
# ==============================================================================

calculate_cluster_size_statistics <- function(
    metadata_table,
    cluster_column,
    sample_order
) {

  cluster_values <- as.character(
    metadata_table[[cluster_column]]
  )

  cluster_counts <- table(
    cluster_values
  )

  cluster_percentages <-
    100 * cluster_counts /
    sum(cluster_counts)

  cluster_by_sample <- table(
    factor(
      metadata_table$sample_ID,
      levels = sample_order
    ),
    cluster_values
  )

  clusters_present_all_samples <- sum(
    colSums(cluster_by_sample > 0) ==
      length(sample_order)
  )

  n_clusters <- length(
    cluster_counts
  )

  data.frame(
    nClusters = n_clusters,
    minClusterN = min(cluster_counts),
    medianClusterN = stats::median(cluster_counts),
    meanClusterN = mean(cluster_counts),
    sdClusterN = stats::sd(cluster_counts),
    maxClusterN = max(cluster_counts),
    smallestClusterPercent =
      min(cluster_percentages),
    largestClusterPercent =
      max(cluster_percentages),
    clusterSizeCV = if (
      mean(cluster_counts) > 0
    ) {
      stats::sd(cluster_counts) /
        mean(cluster_counts)
    } else {
      NA_real_
    },
    nClustersPresentAllSamples =
      clusters_present_all_samples,
    nClustersAbsentAtLeastOneSample =
      n_clusters -
      clusters_present_all_samples,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


# ==============================================================================
# 9. Spatial validation for one clustering solution
# ==============================================================================

calculate_spatial_validation_for_clustering <- function(
    seurat_object,
    cluster_column,
    coordinates_by_sample,
    spatial_silhouette_spots_by_sample,
    spatial_distance_by_sample = NULL,
    pas_k = 10L,
    pas_minimum_different_neighbours = 6L,
    normalize_chaos_by_spot_pitch = TRUE
) {

  metadata_table <- seurat_object[[]]

  sample_ids <- names(
    coordinates_by_sample
  )

  sample_rows <- vector(
    mode = "list",
    length = length(sample_ids)
  )

  for (sample_index in seq_along(sample_ids)) {

    sample_id <- sample_ids[[sample_index]]

    coordinate_table <-
      coordinates_by_sample[[sample_id]]

    sample_spots <- rownames(
      coordinate_table
    )

    cluster_labels <- as.character(
      metadata_table[
        sample_spots,
        cluster_column,
        drop = TRUE
      ]
    )

    names(cluster_labels) <-
      sample_spots

    selected_spatial_silhouette_spots <-
      spatial_silhouette_spots_by_sample[[sample_id]]

    current_spatial_distance <- if (
      !is.null(spatial_distance_by_sample)
    ) {
      spatial_distance_by_sample[[sample_id]]
    } else {
      NULL
    }

    spatial_silhouette <-
      calculate_spatial_silhouette(
        coordinate_table = coordinate_table,
        cluster_labels = cluster_labels,
        selected_spots =
          selected_spatial_silhouette_spots,
        distance_object =
          current_spatial_distance
      )

    chaos_result <- calculate_chaos_score(
      coordinate_matrix = coordinate_table[
        ,
        c("x", "y"),
        drop = FALSE
      ],
      cluster_labels = cluster_labels,
      normalize_by_spot_pitch =
        normalize_chaos_by_spot_pitch
    )

    pas_result <- calculate_pas_score(
      coordinate_matrix = coordinate_table[
        ,
        c("x", "y"),
        drop = FALSE
      ],
      cluster_labels = cluster_labels,
      k = pas_k,
      minimum_different_neighbours =
        pas_minimum_different_neighbours
    )

    sample_rows[[sample_index]] <- data.frame(
      sample_ID = sample_id,
      clusterColumn = cluster_column,
      nSpots = length(sample_spots),
      nClustersPresent =
        length(unique(cluster_labels)),
      spatialASW_mean =
        spatial_silhouette$mean,
      spatialASW_median =
        spatial_silhouette$median,
      spatialASW_sd =
        spatial_silhouette$sd,
      spatialASW_fractionNegative =
        spatial_silhouette$fractionNegative,
      spatialASW_nEvaluated =
        spatial_silhouette$nEvaluated,
      CHAOS = chaos_result$CHAOS,
      spotPitch =
        chaos_result$spotPitch,
      CHAOS_nValidSpots =
        chaos_result$nValidSpots,
      CHAOS_nExcludedSingletonClusterSpots =
        chaos_result$nExcludedSingletonClusterSpots,
      PAS_fraction =
        pas_result$PAS_fraction,
      PAS_percent =
        pas_result$PAS_percent,
      PAS_nAbnormalSpots =
        pas_result$nAbnormalSpots,
      PAS_k =
        pas_result$k,
      PAS_minimumDifferentNeighbours =
        pas_result$minimumDifferentNeighbours,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  spatial_table <- do.call(
    rbind,
    sample_rows
  )

  rownames(spatial_table) <- NULL

  spatial_table
}


aggregate_spatial_validation <- function(
    spatial_validation_table
) {

  data.frame(
    spatialASW_meanAcrossSamples = mean(
      spatial_validation_table$spatialASW_mean,
      na.rm = TRUE
    ),
    spatialASW_sdAcrossSamples = stats::sd(
      spatial_validation_table$spatialASW_mean,
      na.rm = TRUE
    ),
    spatialASW_medianAcrossSamples = stats::median(
      spatial_validation_table$spatialASW_mean,
      na.rm = TRUE
    ),
    spatialASW_meanFractionNegativeAcrossSamples =
      mean(
        spatial_validation_table$spatialASW_fractionNegative,
        na.rm = TRUE
      ),
    CHAOS_meanAcrossSamples = mean(
      spatial_validation_table$CHAOS,
      na.rm = TRUE
    ),
    CHAOS_sdAcrossSamples = stats::sd(
      spatial_validation_table$CHAOS,
      na.rm = TRUE
    ),
    CHAOS_medianAcrossSamples = stats::median(
      spatial_validation_table$CHAOS,
      na.rm = TRUE
    ),
    PAS_meanAcrossSamples = mean(
      spatial_validation_table$PAS_fraction,
      na.rm = TRUE
    ),
    PAS_sdAcrossSamples = stats::sd(
      spatial_validation_table$PAS_fraction,
      na.rm = TRUE
    ),
    PAS_medianAcrossSamples = stats::median(
      spatial_validation_table$PAS_fraction,
      na.rm = TRUE
    ),
    PAS_meanPercentAcrossSamples = mean(
      spatial_validation_table$PAS_percent,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


# ==============================================================================
# 10. Plot-table preparation
# ==============================================================================

prepare_validation_metric_plot_table <- function(
    validation_summary
) {

  metric_specification <- data.frame(
    valueColumn = c(
      "transcriptomicASW_mean",
      "spatialASW_meanAcrossSamples",
      "CHAOS_meanAcrossSamples",
      "PAS_meanAcrossSamples",
      "ARI_vsPreviousResolution",
      "nClusters",
      "smallestClusterPercent",
      "largestClusterPercent"
    ),
    sdColumn = c(
      "transcriptomicASW_sd",
      "spatialASW_sdAcrossSamples",
      "CHAOS_sdAcrossSamples",
      "PAS_sdAcrossSamples",
      NA,
      NA,
      NA,
      NA
    ),
    metricLabel = c(
      "Transcriptomic ASW (higher = better)",
      "Spatial ASW, mean across samples (higher = better)",
      "CHAOS, mean across samples (lower = better)",
      "PAS, mean across samples (lower = better)",
      "ARI versus previous resolution (higher = more stable)",
      "Number of clusters",
      "Smallest cluster (% of spots)",
      "Largest cluster (% of spots)"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  plot_rows <- lapply(
    seq_len(nrow(metric_specification)),
    function(metric_index) {

      value_column <-
        metric_specification$valueColumn[[metric_index]]

      sd_column <-
        metric_specification$sdColumn[[metric_index]]

      data.frame(
        algorithm =
          validation_summary$algorithm,
        resolution =
          validation_summary$resolution,
        clusterColumn =
          validation_summary$clusterColumn,
        metric =
          metric_specification$metricLabel[[metric_index]],
        value =
          validation_summary[[value_column]],
        valueSD = if (
          !is.na(sd_column) &&
            sd_column %in%
            colnames(validation_summary)
        ) {
          validation_summary[[sd_column]]
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )

  plot_table <- do.call(
    rbind,
    plot_rows
  )

  rownames(plot_table) <- NULL

  plot_table$metric <- factor(
    plot_table$metric,
    levels =
      metric_specification$metricLabel
  )

  plot_table
}


create_validation_metrics_plot <- function(
    validation_summary,
    algorithm = NULL,
    plot_title,
    plot_subtitle,
    combined_methods = FALSE
) {

  plot_table <- prepare_validation_metric_plot_table(
    validation_summary
  )

  if (!is.null(algorithm)) {
    plot_table <- plot_table[
      plot_table$algorithm ==
        algorithm,
      ,
      drop = FALSE
    ]
  }

  if (nrow(plot_table) == 0L) {
    stop(
      "No rows are available for validation plotting.",
      call. = FALSE
    )
  }

  # Error bars are intentionally restricted:
  #
  # - combined all-methods plot:
  #     no error bars, because four overlapping algorithms make the figure
  #     difficult to read;
  #
  # - separate algorithm plots:
  #     mean ± SD is shown only for metrics aggregated across tissue sections:
  #       Spatial ASW, CHAOS and PAS;
  #
  # - transcriptomic ASW:
  #     its SD describes variation among individual spots and is retained in
  #     TSV tables, but is not shown as an uncertainty interval on the plot.
  error_bar_metrics <- c(
    "Spatial ASW, mean across samples (higher = better)",
    "CHAOS, mean across samples (lower = better)",
    "PAS, mean across samples (lower = better)"
  )

  plot_table$showErrorBar <-
    !isTRUE(combined_methods) &
    as.character(plot_table$metric) %in%
      error_bar_metrics &
    is.finite(plot_table$valueSD)

  plot_table$ymin <- ifelse(
    plot_table$showErrorBar,
    plot_table$value -
      plot_table$valueSD,
    NA_real_
  )

  plot_table$ymax <- ifelse(
    plot_table$showErrorBar,
    plot_table$value +
      plot_table$valueSD,
    NA_real_
  )

  error_bar_table <- plot_table[
    plot_table$showErrorBar,
    ,
    drop = FALSE
  ]

  if (isTRUE(combined_methods)) {

    output_plot <- ggplot2::ggplot(
      plot_table,
      ggplot2::aes(
        x = resolution,
        y = value,
        group = algorithm,
        colour = algorithm
      )
    )

  } else {

    output_plot <- ggplot2::ggplot(
      plot_table,
      ggplot2::aes(
        x = resolution,
        y = value,
        group = 1
      )
    )
  }

  output_plot <- output_plot +
    ggplot2::geom_line(
      linewidth = 0.7,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      size = 3,
      na.rm = TRUE
    )

  if (nrow(error_bar_table) > 0L) {

    if (isTRUE(combined_methods)) {

      output_plot <- output_plot +
        ggplot2::geom_errorbar(
          data = error_bar_table,
          ggplot2::aes(
            x = resolution,
            y = value,
            ymin = ymin,
            ymax = ymax,
            group = algorithm,
            colour = algorithm
          ),
          width = 0.015,
          linewidth = 0.6,
          inherit.aes = FALSE,
          na.rm = TRUE
        )

    } else {

      output_plot <- output_plot +
        ggplot2::geom_errorbar(
          data = error_bar_table,
          ggplot2::aes(
            x = resolution,
            y = value,
            ymin = ymin,
            ymax = ymax
          ),
          width = 0.015,
          linewidth = 0.6,
          inherit.aes = FALSE,
          na.rm = TRUE
        )
    }
  }

  output_plot +
    ggplot2::facet_wrap(
      facets = ggplot2::vars(metric),
      scales = "free_y",
      ncol = 2
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(
        unique(
          plot_table$resolution
        )
      )
    ) +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Clustering resolution",
      y = "Metric value",
      colour = "Clustering algorithm"
    ) +
    ggplot2::theme_classic(
      base_size = 11
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 17,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        size = 9,
        hjust = 0.5,
        lineheight = 1.08,
        margin = ggplot2::margin(
          b = 10
        )
      ),
      strip.background = ggplot2::element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.6
      ),
      strip.text = ggplot2::element_text(
        size = 9,
        face = "bold"
      ),
      panel.border = ggplot2::element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.5
      ),
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = if (
        isTRUE(combined_methods)
      ) {
        "top"
      } else {
        "none"
      },
      plot.margin = ggplot2::margin(
        t = 10,
        r = 10,
        b = 10,
        l = 10
      )
    )
}


# ==============================================================================
# 11. Main validation workflow
# ==============================================================================

run_clustering_validation <- function(
    seurat_object,
    cluster_columns,
    reduction = "integrated.rpca",
    dims = 1:20,
    sample_order,
    output_dir,
    analysis_prefix =
      "logNormalizeVst_hvg2000_rpcaDims20K20PruneSNN0067",
    integration_method = "RPCA",
    normalization_label = "LogNormalize + VST",
    n_hvg = 2000L,
    k_param = 20L,
    prune_snn = 1 / 15,
    transcriptomic_silhouette_max_spots = 4000L,
    transcriptomic_silhouette_minimum_spots_per_cluster = 20L,
    spatial_silhouette_max_spots_per_sample = 3000L,
    spatial_silhouette_minimum_spots_per_cluster = 10L,
    pas_k = 10L,
    pas_minimum_different_neighbours = 6L,
    normalize_chaos_by_spot_pitch = TRUE,
    ari_comparison_mode = c(
      "all",
      "withinAlgorithmAndSameResolution",
      "adjacentAndSameResolution"
    ),
    seed = 7L,
    save_png = TRUE,
    save_pdf = TRUE,
    png_width_in = 14,
    png_height_in = 13,
    pdf_width_in = 14,
    pdf_height_in = 13,
    dpi = 300,
    verbose = TRUE
) {

  # ============================================================================
  # 11.1 Required packages and input checks
  # ============================================================================

  ari_comparison_mode <- match.arg(
    ari_comparison_mode
  )

  required_packages <- c(
    "Seurat",
    "SeuratObject",
    "cluster",
    "RANN",
    "ggplot2"
  )

  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_packages) > 0L) {
    stop(
      "Missing required R package(s): ",
      paste(
        missing_packages,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  workflow_started <- Sys.time()

  solution_table <- build_clustering_solution_table(
    cluster_columns
  )

  cluster_columns <- solution_table$clusterColumn

  validate_clustering_validation_inputs(
    seurat_object = seurat_object,
    cluster_columns = cluster_columns,
    reduction = reduction,
    dims = dims,
    sample_order = sample_order
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  tables_directory <- file.path(
    output_dir,
    "tables"
  )

  figures_directory <- file.path(
    output_dir,
    "figures"
  )

  dir.create(
    tables_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    figures_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  metadata_table <- seurat_object[[]]

  if (verbose) {
    clustering_validation_message(
      "Starting clustering validation for ",
      length(cluster_columns),
      " solutions."
    )
  }


  # ============================================================================
  # 11.2 Select common transcriptomic and spatial silhouette samples
  # ============================================================================

  transcriptomic_spots <-
    select_common_transcriptomic_validation_spots(
      metadata_table = metadata_table,
      cluster_columns = cluster_columns,
      max_spots =
        transcriptomic_silhouette_max_spots,
      minimum_spots_per_cluster =
        transcriptomic_silhouette_minimum_spots_per_cluster,
      seed = seed,
      verbose = verbose
    )

  spatial_data <- extract_all_spatial_coordinates(
    seurat_object = seurat_object,
    sample_order = sample_order,
    verbose = verbose
  )

  coordinates_by_sample <-
    spatial_data$coordinatesBySample

  image_map <- spatial_data$imageMap

  spatial_silhouette_spots_by_sample <-
    vector(
      mode = "list",
      length = length(sample_order)
    )

  names(
    spatial_silhouette_spots_by_sample
  ) <- sample_order

  for (sample_index in seq_along(sample_order)) {

    sample_id <- sample_order[[sample_index]]

    sample_spots <- rownames(
      coordinates_by_sample[[sample_id]]
    )

    spatial_silhouette_spots_by_sample[[sample_id]] <- select_common_spatial_silhouette_spots(
      sample_spots = sample_spots,
      metadata_table = metadata_table,
      cluster_columns = cluster_columns,
      max_spots =
        spatial_silhouette_max_spots_per_sample,
      minimum_spots_per_cluster =
        spatial_silhouette_minimum_spots_per_cluster,
      seed =
        as.integer(seed) +
        sample_index
    )
  }

  if (verbose) {
    clustering_validation_message(
      "Precomputing one transcriptomic distance matrix shared by all ",
      length(cluster_columns),
      " clustering solutions."
    )
  }

  transcriptomic_embedding <- SeuratObject::Embeddings(
    seurat_object[[reduction]]
  )[
    transcriptomic_spots,
    dims,
    drop = FALSE
  ]

  transcriptomic_distance <- stats::dist(
    transcriptomic_embedding,
    method = "euclidean"
  )

  rm(transcriptomic_embedding)

  if (verbose) {
    clustering_validation_message(
      "Precomputing one spatial distance matrix per sample."
    )
  }

  spatial_distance_by_sample <- lapply(
    sample_order,
    function(sample_id) {

      selected_spots <-
        spatial_silhouette_spots_by_sample[[sample_id]]

      selected_coordinates <-
        coordinates_by_sample[[sample_id]][
          selected_spots,
          c("x", "y"),
          drop = FALSE
        ]

      stats::dist(
        selected_coordinates,
        method = "euclidean"
      )
    }
  )

  names(spatial_distance_by_sample) <-
    sample_order

  invisible(
    gc(verbose = FALSE)
  )


  # ============================================================================
  # 11.3 Calculate validation metrics for every clustering solution
  # ============================================================================

  summary_rows <- vector(
    mode = "list",
    length = nrow(solution_table)
  )

  spatial_rows <- vector(
    mode = "list",
    length = nrow(solution_table)
  )

  for (
    solution_index in
    seq_len(nrow(solution_table))
  ) {

    cluster_column <-
      solution_table$clusterColumn[[solution_index]]

    algorithm_name <-
      solution_table$algorithm[[solution_index]]

    resolution_value <-
      solution_table$resolution[[solution_index]]

    solution_started <- Sys.time()

    if (verbose) {
      clustering_validation_message(
        "[",
        solution_index,
        "/",
        nrow(solution_table),
        "] START | ",
        cluster_column
      )
    }

    transcriptomic_silhouette <-
      calculate_transcriptomic_silhouette(
        seurat_object = seurat_object,
        cluster_column = cluster_column,
        reduction = reduction,
        dims = dims,
        selected_spots =
          transcriptomic_spots,
        distance_object =
          transcriptomic_distance
      )

    spatial_validation <-
      calculate_spatial_validation_for_clustering(
        seurat_object = seurat_object,
        cluster_column = cluster_column,
        coordinates_by_sample =
          coordinates_by_sample,
        spatial_silhouette_spots_by_sample =
          spatial_silhouette_spots_by_sample,
        spatial_distance_by_sample =
          spatial_distance_by_sample,
        pas_k = pas_k,
        pas_minimum_different_neighbours =
          pas_minimum_different_neighbours,
        normalize_chaos_by_spot_pitch =
          normalize_chaos_by_spot_pitch
      )

    spatial_validation$algorithm <-
      algorithm_name

    spatial_validation$resolution <-
      resolution_value

    spatial_validation <- spatial_validation[
      ,
      c(
        "algorithm",
        "resolution",
        setdiff(
          colnames(spatial_validation),
          c(
            "algorithm",
            "resolution"
          )
        )
      ),
      drop = FALSE
    ]

    spatial_aggregate <-
      aggregate_spatial_validation(
        spatial_validation
      )

    cluster_size_statistics <-
      calculate_cluster_size_statistics(
        metadata_table = metadata_table,
        cluster_column = cluster_column,
        sample_order = sample_order
      )

    solution_finished <- Sys.time()

    elapsed_seconds <- as.numeric(
      difftime(
        solution_finished,
        solution_started,
        units = "secs"
      )
    )

    summary_rows[[solution_index]] <-
      data.frame(
        algorithm = algorithm_name,
        resolution = resolution_value,
        clusterColumn = cluster_column,
        nSpots = ncol(seurat_object),
        transcriptomicASW_mean =
          transcriptomic_silhouette$mean,
        transcriptomicASW_median =
          transcriptomic_silhouette$median,
        transcriptomicASW_sd =
          transcriptomic_silhouette$sd,
        transcriptomicASW_fractionNegative =
          transcriptomic_silhouette$fractionNegative,
        transcriptomicASW_nEvaluated =
          transcriptomic_silhouette$nEvaluated,
        spatial_aggregate,
        cluster_size_statistics,
        validationElapsedSeconds =
          elapsed_seconds,
        validationElapsedFormatted =
          format_validation_elapsed_time(
            elapsed_seconds
          ),
        integrationMethod =
          integration_method,
        normalization =
          normalization_label,
        nHVG = as.integer(n_hvg),
        reduction = reduction,
        dimsStart = min(dims),
        dimsEnd = max(dims),
        kParam = as.integer(k_param),
        pruneSNN = as.numeric(prune_snn),
        PAS_k = as.integer(pas_k),
        PAS_minimumDifferentNeighbours =
          as.integer(
            pas_minimum_different_neighbours
          ),
        ARIcomparisonMode =
          ari_comparison_mode,
        seed = as.integer(seed),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

    spatial_rows[[solution_index]] <-
      spatial_validation

    if (verbose) {
      clustering_validation_message(
        "[",
        solution_index,
        "/",
        nrow(solution_table),
        "] DONE  | ",
        cluster_column,
        " | elapsed=",
        format_validation_elapsed_time(
          elapsed_seconds
        )
      )
    }

    invisible(
      gc(verbose = FALSE)
    )
  }

  validation_summary <- do.call(
    rbind,
    summary_rows
  )

  rownames(validation_summary) <- NULL

  spatial_validation_by_sample <- do.call(
    rbind,
    spatial_rows
  )

  rownames(spatial_validation_by_sample) <-
    NULL


  # ============================================================================
  # 11.4 ARI stability
  # ============================================================================

  pairwise_ari <- calculate_pairwise_ari_table(
    metadata_table = metadata_table,
    solution_table = solution_table,
    comparison_mode =
      ari_comparison_mode
  )

  validation_summary <-
    add_adjacent_resolution_ari(
      validation_summary =
        validation_summary,
      pairwise_ari = pairwise_ari
    )

  validation_summary <- validation_summary[
    order(
      validation_summary$algorithm,
      validation_summary$resolution
    ),
    ,
    drop = FALSE
  ]

  rownames(validation_summary) <- NULL


  # ============================================================================
  # 11.5 Parameter and combined tables
  # ============================================================================

  parameter_table <- data.frame(
    parameter = c(
      "integration_method",
      "normalization",
      "n_hvg",
      "reduction",
      "dims",
      "k_param",
      "prune_snn",
      "transcriptomic_silhouette_common_sample_n",
      "transcriptomic_silhouette_max_spots",
      "transcriptomic_silhouette_minimum_spots_per_cluster",
      "spatial_silhouette_max_spots_per_sample",
      "spatial_silhouette_minimum_spots_per_cluster",
      "pas_k",
      "pas_minimum_different_neighbours",
      "normalize_chaos_by_spot_pitch",
      "ari_comparison_mode",
      "seed"
    ),
    value = c(
      integration_method,
      normalization_label,
      as.character(n_hvg),
      reduction,
      paste0(
        min(dims),
        "-",
        max(dims)
      ),
      as.character(k_param),
      formatC(
        prune_snn,
        format = "f",
        digits = 6L
      ),
      as.character(
        length(transcriptomic_spots)
      ),
      as.character(
        transcriptomic_silhouette_max_spots
      ),
      as.character(
        transcriptomic_silhouette_minimum_spots_per_cluster
      ),
      as.character(
        spatial_silhouette_max_spots_per_sample
      ),
      as.character(
        spatial_silhouette_minimum_spots_per_cluster
      ),
      as.character(pas_k),
      as.character(
        pas_minimum_different_neighbours
      ),
      as.character(
        normalize_chaos_by_spot_pitch
      ),
      ari_comparison_mode,
      as.character(seed)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  output_files <- list()

  output_files$combinedSummary <-
    write_validation_tsv(
      data = validation_summary,
      filename = file.path(
        tables_directory,
        paste0(
          analysis_prefix,
          "_allMethods_clusteringValidationSummary.tsv"
        )
      )
    )

  output_files$combinedSpatialBySample <-
    write_validation_tsv(
      data = spatial_validation_by_sample,
      filename = file.path(
        tables_directory,
        paste0(
          analysis_prefix,
          "_allMethods_spatialValidationBySample.tsv"
        )
      )
    )

  output_files$pairwiseARI <-
    write_validation_tsv(
      data = pairwise_ari,
      filename = file.path(
        tables_directory,
        paste0(
          analysis_prefix,
          "_clusteringStabilityARI.tsv"
        )
      )
    )

  output_files$parameters <-
    write_validation_tsv(
      data = parameter_table,
      filename = file.path(
        tables_directory,
        paste0(
          analysis_prefix,
          "_clusteringValidationParameters.tsv"
        )
      )
    )

  output_files$imageMap <-
    write_validation_tsv(
      data = image_map,
      filename = file.path(
        tables_directory,
        paste0(
          analysis_prefix,
          "_clusteringValidationImageMap.tsv"
        )
      )
    )


  # ============================================================================
  # 11.6 Separate tables and plots for every clustering algorithm
  # ============================================================================

  parameter_subtitle <-
    build_validation_parameter_subtitle(
      integration_method =
        integration_method,
      normalization_label =
        normalization_label,
      n_hvg = n_hvg,
      reduction = reduction,
      dims = dims,
      k_param = k_param,
      prune_snn = prune_snn,
      transcriptomic_sample_n =
        length(transcriptomic_spots),
      spatial_asw_max_spots_per_sample =
        spatial_silhouette_max_spots_per_sample,
      pas_k = pas_k,
      pas_min_different =
        pas_minimum_different_neighbours,
      ari_comparison_mode =
        ari_comparison_mode,
      seed = seed
    )

  algorithms <- unique(
    validation_summary$algorithm
  )

  output_files$byAlgorithm <- list()

  for (algorithm_name in algorithms) {

    algorithm_summary <-
      validation_summary[
        validation_summary$algorithm ==
          algorithm_name,
        ,
        drop = FALSE
      ]

    algorithm_spatial <-
      spatial_validation_by_sample[
        spatial_validation_by_sample$algorithm ==
          algorithm_name,
        ,
        drop = FALSE
      ]

    algorithm_ari <- pairwise_ari[
      pairwise_ari$algorithm1 ==
        algorithm_name &
        pairwise_ari$algorithm2 ==
        algorithm_name,
      ,
      drop = FALSE
    ]

    algorithm_table_prefix <- file.path(
      tables_directory,
      paste0(
        analysis_prefix,
        "_",
        algorithm_name
      )
    )

    algorithm_figure_prefix <- file.path(
      figures_directory,
      paste0(
        analysis_prefix,
        "_",
        algorithm_name,
        "_clusteringValidationMetrics"
      )
    )

    algorithm_outputs <- list()

    algorithm_outputs$summary <-
      write_validation_tsv(
        data = algorithm_summary,
        filename = paste0(
          algorithm_table_prefix,
          "_clusteringValidationSummary.tsv"
        )
      )

    algorithm_outputs$spatialBySample <-
      write_validation_tsv(
        data = algorithm_spatial,
        filename = paste0(
          algorithm_table_prefix,
          "_spatialValidationBySample.tsv"
        )
      )

    algorithm_outputs$ARI <-
      write_validation_tsv(
        data = algorithm_ari,
        filename = paste0(
          algorithm_table_prefix,
          "_clusteringStabilityARI.tsv"
        )
      )

    algorithm_plot <-
      create_validation_metrics_plot(
        validation_summary =
          algorithm_summary,
        algorithm = algorithm_name,
        plot_title = paste0(
          "Clustering validation: ",
          algorithm_name,
          " | ",
          nrow(algorithm_summary),
          " resolutions"
        ),
        plot_subtitle =
          parameter_subtitle,
        combined_methods = FALSE
      )

    algorithm_outputs$plot <-
      save_validation_plot(
        plot_object = algorithm_plot,
        output_prefix =
          algorithm_figure_prefix,
        save_png = save_png,
        save_pdf = save_pdf,
        png_width_in =
          png_width_in,
        png_height_in =
          png_height_in,
        pdf_width_in =
          pdf_width_in,
        pdf_height_in =
          pdf_height_in,
        dpi = dpi
      )

    output_files$byAlgorithm[[algorithm_name]] <- algorithm_outputs
  }


  # ============================================================================
  # 11.7 Combined all-methods plot
  # ============================================================================

  combined_plot <-
    create_validation_metrics_plot(
      validation_summary =
        validation_summary,
      algorithm = NULL,
      plot_title = paste0(
        "Clustering validation across ",
        length(algorithms),
        " algorithms and ",
        length(unique(validation_summary$resolution)),
        " resolutions"
      ),
      plot_subtitle =
        parameter_subtitle,
      combined_methods = TRUE
    )

  output_files$combinedPlot <-
    save_validation_plot(
      plot_object = combined_plot,
      output_prefix = file.path(
        figures_directory,
        paste0(
          analysis_prefix,
          "_allMethods_clusteringValidationMetrics"
        )
      ),
      save_png = save_png,
      save_pdf = save_pdf,
      png_width_in =
        png_width_in,
      png_height_in =
        png_height_in,
      pdf_width_in =
        pdf_width_in,
      pdf_height_in =
        pdf_height_in,
      dpi = dpi
    )


  # ============================================================================
  # 11.8 Final report
  # ============================================================================

  workflow_finished <- Sys.time()

  workflow_elapsed_seconds <- as.numeric(
    difftime(
      workflow_finished,
      workflow_started,
      units = "secs"
    )
  )

  if (verbose) {
    clustering_validation_message(
      "Clustering validation completed in ",
      format_validation_elapsed_time(
        workflow_elapsed_seconds
      ),
      "."
    )

    clustering_validation_message(
      "Output directory: ",
      normalizePath(
        output_dir,
        mustWork = TRUE
      )
    )
  }

  list(
    validationSummary =
      validation_summary,
    spatialValidationBySample =
      spatial_validation_by_sample,
    pairwiseARI =
      pairwise_ari,
    parameterTable =
      parameter_table,
    imageMap =
      image_map,
    transcriptomicSilhouetteSpots =
      transcriptomic_spots,
    spatialSilhouetteSpotsBySample =
      spatial_silhouette_spots_by_sample,
    outputFiles =
      output_files,
    workflowElapsedSeconds =
      workflow_elapsed_seconds
  )
}
