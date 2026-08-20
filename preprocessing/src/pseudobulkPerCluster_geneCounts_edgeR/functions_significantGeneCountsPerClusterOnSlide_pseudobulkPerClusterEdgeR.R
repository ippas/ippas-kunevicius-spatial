# ==============================================================================
# functions_significantGeneCountsPerClusterOnSlide_pseudobulkPerClusterEdgeR.R
#
# Version: fixed4_nestedPatchworkGrob_2026-08-03
#
# Dedicated functions for visualizing the number of statistically significant
# genes per spatial cluster for the cluster-specific edgeR analysis.
#
# Current workflow:
#   main donor-group effect
#   ASD vs Neurotypical averaged equally across Male and Female
#
# For one selected Visium section, the same section is displayed three times:
#   1. total significant results per cluster (UP + DOWN), green scale;
#   2. UP results per cluster, red scale;
#   3. DOWN results per cluster, blue scale.
#
# Colour intensity represents the number of significant cluster-gene results.
# Each panel contains a clustering-style legend with one entry per cluster and
# the corresponding number of significant genes.
#
# The calling script supplies:
#   - one sample_ID to plot;
#   - FDR thresholds;
#   - absolute log2FC thresholds;
#   - paths to the Seurat and edgeR RData files.
# ============================================================================== 


# ==============================================================================
# 1. General helpers
# ==============================================================================

sort_spatial_edger_cluster_ids <- function(cluster_ids) {

  cluster_ids <- unique(as.character(cluster_ids))
  numeric_ids <- suppressWarnings(as.numeric(cluster_ids))

  if (!anyNA(numeric_ids)) {
    return(cluster_ids[order(numeric_ids)])
  }

  sort(cluster_ids)
}


write_spatial_edger_tsv <- function(data, filename) {

  utils::write.table(
    x = data,
    file = filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}


format_spatial_edger_number <- function(x, digits = 2L) {

  formatC(
    as.numeric(x),
    format = "f",
    digits = as.integer(digits)
  )
}


format_spatial_edger_integer <- function(x) {

  format(
    as.integer(x),
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


hex_to_rgb01_spatial_edger <- function(hex_colour) {

  rgb_matrix <- grDevices::col2rgb(hex_colour) / 255
  as.numeric(rgb_matrix[, 1])
}


mix_spatial_edger_colour_with_white <- function(base_colour, intensity) {

  intensity <- max(min(as.numeric(intensity), 1), 0)
  base_rgb <- hex_to_rgb01_spatial_edger(base_colour)
  mixed_rgb <- 1 - (1 - base_rgb) * intensity
  grDevices::rgb(mixed_rgb[[1]], mixed_rgb[[2]], mixed_rgb[[3]])
}


build_spatial_edger_threshold_label <- function(
    fdr_threshold,
    abs_log2fc_threshold
) {

  fdr_code <- sprintf(
    "%03d",
    as.integer(round(as.numeric(fdr_threshold) * 100))
  )

  logfc_code <- sprintf(
    "%02d",
    as.integer(round(as.numeric(abs_log2fc_threshold) * 10))
  )

  paste0(
    "FDR",
    fdr_code,
    "_absLog2FC",
    logfc_code
  )
}


compute_spatial_edger_square_limits <- function(
    xmin,
    xmax,
    ymin,
    ymax,
    padding_fraction = 0.03
) {

  x_range <- xmax - xmin
  y_range <- ymax - ymin
  maximum_range <- max(x_range, y_range)

  if (!is.finite(maximum_range) || maximum_range <= 0) {
    maximum_range <- 1
  }

  half_side <- maximum_range / 2 * (1 + padding_fraction)
  x_center <- (xmin + xmax) / 2
  y_center <- (ymin + ymax) / 2

  list(
    x_limits = c(
      x_center - half_side,
      x_center + half_side
    ),
    y_limits = c(
      y_center - half_side,
      y_center + half_side
    )
  )
}


standardize_spatial_edger_coordinates <- function(coordinate_table) {

  coordinate_table <- as.data.frame(
    coordinate_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"cell" %in% colnames(coordinate_table)) {
    coordinate_table$cell <- rownames(coordinate_table)
  }

  x_candidates <- c(
    "imagecol",
    "x",
    "col",
    "pxl_col_in_fullres"
  )

  y_candidates <- c(
    "imagerow",
    "y",
    "row",
    "pxl_row_in_fullres"
  )

  x_column <- x_candidates[x_candidates %in% colnames(coordinate_table)][1]
  y_column <- y_candidates[y_candidates %in% colnames(coordinate_table)][1]

  if (is.na(x_column) || is.na(y_column)) {
    stop(
      "Could not identify x/y columns in tissue coordinates. Available: ",
      paste(colnames(coordinate_table), collapse = ", "),
      call. = FALSE
    )
  }

  data.frame(
    cell = as.character(coordinate_table$cell),
    x_plot = as.numeric(coordinate_table[[x_column]]),
    y_plot = as.numeric(coordinate_table[[y_column]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


# ==============================================================================
# 2. Load the Seurat object and edgeR results
# ==============================================================================

load_single_seurat_for_significant_gene_count_plot <- function(
    input_seurat_rdata_file,
    seurat_object_name = NULL
) {

  if (!file.exists(input_seurat_rdata_file)) {
    stop(
      "Input Seurat RData file does not exist:\n",
      input_seurat_rdata_file,
      call. = FALSE
    )
  }

  load_environment <- new.env(parent = globalenv())

  loaded_object_names <- load(
    input_seurat_rdata_file,
    envir = load_environment
  )

  if (!is.null(seurat_object_name)) {

    if (!seurat_object_name %in% loaded_object_names) {
      stop(
        "Requested Seurat object '",
        seurat_object_name,
        "' is absent from the RData file. Available objects: ",
        paste(loaded_object_names, collapse = ", "),
        call. = FALSE
      )
    }

    selected_object <- get(
      seurat_object_name,
      envir = load_environment,
      inherits = FALSE
    )

    if (!inherits(selected_object, "Seurat")) {
      stop(
        "Requested object is not a Seurat object: ",
        seurat_object_name,
        call. = FALSE
      )
    }

    return(
      list(
        object = selected_object,
        object_name = seurat_object_name
      )
    )
  }

  seurat_object_names <- loaded_object_names[
    vapply(
      loaded_object_names,
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
      logical(1)
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


load_main_group_edger_results_for_significant_gene_count_plot <- function(
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

  required_object <- "edgeR_perCluster_combinedResults"

  if (!required_object %in% loaded_names) {
    stop(
      "The edgeR RData file does not contain ",
      required_object,
      ". Loaded objects: ",
      paste(loaded_names, collapse = ", "),
      call. = FALSE
    )
  }

  combined_results <- get(
    required_object,
    envir = result_environment,
    inherits = FALSE
  )

  main_group_test_id <- "Overall_Group_ASD_vs_Neurotypical"

  if (!main_group_test_id %in% names(combined_results)) {
    stop(
      "Main group result table is absent: ",
      main_group_test_id,
      "\nAvailable result tables: ",
      paste(names(combined_results), collapse = ", "),
      call. = FALSE
    )
  }

  main_group_results <- combined_results[[main_group_test_id]]

  required_columns <- c(
    "cluster_id",
    "ensembl_gene_id",
    "gene",
    "logFC",
    "PValue",
    "FDR"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(main_group_results)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Missing main-group result column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  main_group_results$cluster_id <- as.character(
    main_group_results$cluster_id
  )

  list(
    main_group_results = main_group_results,
    main_group_test_id = main_group_test_id,
    loaded_names = loaded_names
  )
}


# ==============================================================================
# 3. Prepare one selected Visium section
# ==============================================================================

prepare_selected_section_for_significant_gene_count_plot <- function(
    seurat_object,
    sample_id_to_plot,
    cluster_column,
    sample_id_column = "sample_ID",
    image_scale = "lowres",
    panel_padding_fraction = 0.03
) {

  if (!inherits(seurat_object, "Seurat")) {
    stop("`seurat_object` must be a Seurat object.", call. = FALSE)
  }

  metadata_table <- seurat_object[[]]

  required_columns <- c(
    sample_id_column,
    cluster_column
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(metadata_table)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Missing Seurat metadata column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  sample_id_to_plot <- as.character(sample_id_to_plot)[[1]]

  available_sample_ids <- unique(
    as.character(metadata_table[[sample_id_column]])
  )

  if (!sample_id_to_plot %in% available_sample_ids) {
    stop(
      "Requested sample_ID is absent from the Seurat object: ",
      sample_id_to_plot,
      "\nAvailable sample IDs: ",
      paste(available_sample_ids, collapse = ", "),
      call. = FALSE
    )
  }

  image_names <- SeuratObject::Images(seurat_object)

  if (length(image_names) == 0L) {
    stop("No spatial images are present in the Seurat object.", call. = FALSE)
  }

  image_mapping_rows <- lapply(
    image_names,
    function(image_name) {

      image_cells <- SeuratObject::Cells(
        seurat_object[[image_name]]
      )

      shared_cells <- intersect(
        image_cells,
        rownames(metadata_table)
      )

      image_sample_ids <- unique(
        as.character(
          metadata_table[
            shared_cells,
            sample_id_column,
            drop = TRUE
          ]
        )
      )

      if (length(image_sample_ids) != 1L) {
        stop(
          "Spatial image ",
          image_name,
          " does not map to exactly one sample_ID.",
          call. = FALSE
        )
      }

      data.frame(
        sample_ID = image_sample_ids[[1]],
        image_name = image_name,
        n_spots = length(shared_cells),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )

  image_mapping <- do.call(
    rbind,
    image_mapping_rows
  )

  matching_images <- image_mapping$image_name[
    image_mapping$sample_ID == sample_id_to_plot
  ]

  if (length(matching_images) != 1L) {
    stop(
      "Expected exactly one spatial image for sample ",
      sample_id_to_plot,
      ", but found ",
      length(matching_images),
      ".",
      call. = FALSE
    )
  }

  image_name <- matching_images[[1]]

  image_cells <- SeuratObject::Cells(
    seurat_object[[image_name]]
  )

  coordinate_table <- SeuratObject::GetTissueCoordinates(
    object = seurat_object,
    image = image_name,
    scale = image_scale
  )

  coordinate_table <- standardize_spatial_edger_coordinates(
    coordinate_table
  )

  common_cells <- image_cells[
    image_cells %in% coordinate_table$cell &
      image_cells %in% rownames(metadata_table)
  ]

  if (length(common_cells) == 0L) {
    stop(
      "No shared cells among image, coordinates and metadata for sample ",
      sample_id_to_plot,
      ".",
      call. = FALSE
    )
  }

  coordinate_table <- coordinate_table[
    match(common_cells, coordinate_table$cell),
    ,
    drop = FALSE
  ]

  spot_metadata <- metadata_table[
    common_cells,
    ,
    drop = FALSE
  ]

  coordinate_table$sample_ID <- sample_id_to_plot
  coordinate_table$cluster_id <- as.character(
    spot_metadata[[cluster_column]]
  )

  complete_cluster_levels <- sort_spatial_edger_cluster_ids(
    metadata_table[[cluster_column]]
  )

  coordinate_table$cluster_id <- factor(
    coordinate_table$cluster_id,
    levels = complete_cluster_levels
  )

  raw_image <- tryCatch(
    SeuratObject::GetImage(
      object = seurat_object,
      image = image_name,
      mode = "raw"
    ),
    error = function(error_condition) NULL
  )

  raw_dimensions <- dim(raw_image)

  if (!is.null(raw_dimensions) && length(raw_dimensions) >= 2L) {
    image_height <- as.numeric(raw_dimensions[[1]])
    image_width <- as.numeric(raw_dimensions[[2]])
  } else {
    image_height <- max(coordinate_table$y_plot, na.rm = TRUE)
    image_width <- max(coordinate_table$x_plot, na.rm = TRUE)
  }

  frame_limits <- compute_spatial_edger_square_limits(
    xmin = 0,
    xmax = image_width,
    ymin = 0,
    ymax = image_height,
    padding_fraction = panel_padding_fraction
  )

  list(
    sample_ID = sample_id_to_plot,
    image_name = image_name,
    image_array = raw_image,
    image_height = image_height,
    image_width = image_width,
    plot_data = coordinate_table,
    cluster_levels = complete_cluster_levels,
    frame_limits = frame_limits,
    n_spots = nrow(coordinate_table)
  )
}


# ==============================================================================
# 4. Count significant genes per cluster
# ==============================================================================

summarize_main_group_significant_genes_per_cluster <- function(
    main_group_results,
    cluster_levels,
    fdr_threshold,
    abs_log2fc_threshold
) {

  cluster_levels <- sort_spatial_edger_cluster_ids(cluster_levels)

  significant_results <- main_group_results |>
    dplyr::filter(
      is.finite(.data$FDR),
      is.finite(.data$logFC),
      .data$FDR < fdr_threshold,
      abs(.data$logFC) > abs_log2fc_threshold
    ) |>
    dplyr::mutate(
      direction = dplyr::if_else(
        .data$logFC > 0,
        "UP",
        "DOWN"
      )
    )

  all_cluster_table <- tibble::tibble(
    cluster_id = cluster_levels
  )

  total_counts <- significant_results |>
    dplyr::count(
      .data$cluster_id,
      name = "n_total"
    )

  up_counts <- significant_results |>
    dplyr::filter(.data$direction == "UP") |>
    dplyr::count(
      .data$cluster_id,
      name = "n_up"
    )

  down_counts <- significant_results |>
    dplyr::filter(.data$direction == "DOWN") |>
    dplyr::count(
      .data$cluster_id,
      name = "n_down"
    )

  cluster_count_table <- all_cluster_table |>
    dplyr::left_join(
      total_counts,
      by = "cluster_id"
    ) |>
    dplyr::left_join(
      up_counts,
      by = "cluster_id"
    ) |>
    dplyr::left_join(
      down_counts,
      by = "cluster_id"
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(
          c(
            "n_total",
            "n_up",
            "n_down"
          )
        ),
        ~ dplyr::coalesce(as.integer(.x), 0L)
      )
    )

  if (!all(
    cluster_count_table$n_total ==
      cluster_count_table$n_up + cluster_count_table$n_down
  )) {
    stop(
      "Internal error: total count is not equal to UP + DOWN.",
      call. = FALSE
    )
  }

  summary_table <- tibble::tibble(
    comparison = "ASD vs Neurotypical averaged equally across sexes",
    fdr_threshold = as.numeric(fdr_threshold),
    abs_log2fc_threshold = as.numeric(abs_log2fc_threshold),
    total_significant_results = nrow(significant_results),
    unique_significant_genes = dplyr::n_distinct(
      significant_results$ensembl_gene_id
    ),
    up_significant_results = sum(
      significant_results$direction == "UP"
    ),
    unique_up_genes = dplyr::n_distinct(
      significant_results$ensembl_gene_id[
        significant_results$direction == "UP"
      ]
    ),
    down_significant_results = sum(
      significant_results$direction == "DOWN"
    ),
    unique_down_genes = dplyr::n_distinct(
      significant_results$ensembl_gene_id[
        significant_results$direction == "DOWN"
      ]
    ),
    number_of_clusters = length(cluster_levels)
  )

  list(
    significant_results = significant_results,
    cluster_count_table = cluster_count_table,
    summary_table = summary_table
  )
}


# ==============================================================================
# 5. Colour mapping, dual-dot legend and panel creation
# ==============================================================================

map_spatial_edger_cluster_counts_to_colours <- function(
    cluster_count_table,
    count_column,
    palette_colors,
    zero_colour = "#D9D9D9"
) {

  count_values <- as.numeric(
    cluster_count_table[[count_column]]
  )

  maximum_count <- max(
    count_values,
    na.rm = TRUE
  )

  if (!is.finite(maximum_count) || maximum_count <= 0) {
    mapped_colours <- rep(
      zero_colour,
      length(count_values)
    )
  } else {

    colour_function <- scales::col_numeric(
      palette = palette_colors,
      domain = c(0, maximum_count),
      na.color = zero_colour
    )

    mapped_colours <- colour_function(
      count_values
    )

    mapped_colours[count_values == 0] <- zero_colour
  }

  names(mapped_colours) <- as.character(
    cluster_count_table$cluster_id
  )

  list(
    colours = mapped_colours,
    maximum_count = maximum_count
  )
}


build_spatial_edger_cluster_legend_labels <- function(
    cluster_count_table,
    count_column,
    cluster_names = NULL
) {

  labels <- vapply(
    seq_len(nrow(cluster_count_table)),
    function(row_index) {

      cluster_id <- as.character(
        cluster_count_table$cluster_id[[row_index]]
      )

      count_value <- as.integer(
        cluster_count_table[[count_column]][[row_index]]
      )

      anatomical_label <- ""

      if (
        !is.null(cluster_names) &&
          cluster_id %in% names(cluster_names) &&
          !is.na(cluster_names[[cluster_id]]) &&
          nzchar(cluster_names[[cluster_id]])
      ) {
        anatomical_label <- paste0(
          " — ",
          cluster_names[[cluster_id]]
        )
      }

      paste0(
        "C",
        cluster_id,
        anatomical_label,
        ": ",
        count_value,
        if (count_value == 1L) " gene" else " genes"
      )
    },
    character(1)
  )

  names(labels) <- as.character(
    cluster_count_table$cluster_id
  )

  labels
}


create_spatial_edger_dual_dot_legend <- function(
    cluster_count_table,
    count_column,
    result_colours,
    cluster_base_colours,
    cluster_names = NULL,
    legend_ncol = 2L,
    legend_text_size = 14.4,
    legend_point_size = 5.8,
    legend_row_spacing = 1.20,
    legend_bottom_margin = 8
) {

  cluster_ids <- as.character(
    cluster_count_table$cluster_id
  )

  legend_row_spacing <- as.numeric(legend_row_spacing)
  legend_bottom_margin <- as.numeric(legend_bottom_margin)

  if (
    length(legend_row_spacing) != 1L ||
      !is.finite(legend_row_spacing) ||
      legend_row_spacing <= 0
  ) {
    stop(
      "`legend_row_spacing` must be one positive finite number.",
      call. = FALSE
    )
  }

  if (
    length(legend_bottom_margin) != 1L ||
      !is.finite(legend_bottom_margin)
  ) {
    stop(
      "`legend_bottom_margin` must be one finite number.",
      call. = FALSE
    )
  }

  if (is.null(cluster_base_colours)) {
    stop(
      "`cluster_base_colours` must be supplied for the dual-dot legend.",
      call. = FALSE
    )
  }

  missing_base_colours <- setdiff(
    cluster_ids,
    names(cluster_base_colours)
  )

  if (length(missing_base_colours) > 0L) {
    stop(
      "Missing anatomical colour(s) for cluster ID(s): ",
      paste(missing_base_colours, collapse = ", "),
      call. = FALSE
    )
  }

  missing_result_colours <- setdiff(
    cluster_ids,
    names(result_colours)
  )

  if (length(missing_result_colours) > 0L) {
    stop(
      "Missing result colour(s) for cluster ID(s): ",
      paste(missing_result_colours, collapse = ", "),
      call. = FALSE
    )
  }

  legend_ncol <- max(
    1L,
    min(
      as.integer(legend_ncol),
      length(cluster_ids)
    )
  )

  legend_nrow <- ceiling(
    length(cluster_ids) / legend_ncol
  )

  entry_index <- seq_along(cluster_ids)

  legend_column <- ceiling(
    entry_index / legend_nrow
  )

  legend_row <- (
    (entry_index - 1L) %% legend_nrow
  ) + 1L

  x_spacing <- 12
  x_base <- (
    legend_column - 1L
  ) * x_spacing

  y_position <- (
    legend_nrow - legend_row + 1L
  ) * legend_row_spacing

  legend_labels <- build_spatial_edger_cluster_legend_labels(
    cluster_count_table = cluster_count_table,
    count_column = count_column,
    cluster_names = cluster_names
  )

  anatomical_dot_data <- data.frame(
    x = x_base + 0.25,
    y = y_position,
    colour = unname(
      cluster_base_colours[cluster_ids]
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  result_dot_data <- data.frame(
    x = x_base + 0.85,
    y = y_position,
    colour = unname(
      result_colours[cluster_ids]
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  text_data <- data.frame(
    x = x_base + 1.45,
    y = y_position,
    label = unname(
      legend_labels[cluster_ids]
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = anatomical_dot_data,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        colour = .data$colour
      ),
      shape = 16,
      size = legend_point_size,
      stroke = 0
    ) +
    ggplot2::geom_point(
      data = result_dot_data,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        colour = .data$colour
      ),
      shape = 16,
      size = legend_point_size,
      stroke = 0
    ) +
    ggplot2::geom_text(
      data = text_data,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        label = .data$label
      ),
      family = "DejaVu Sans",
      size = legend_text_size / 3.88,
      hjust = 0,
      vjust = 0.5
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::coord_cartesian(
      xlim = c(
        -0.20,
        legend_ncol * x_spacing - 0.20
      ),
      ylim = c(
        0.40,
        legend_nrow * legend_row_spacing + 0.60
      ),
      clip = "off"
    ) +
    ggplot2::theme_void(
      base_family = "DejaVu Sans"
    ) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      plot.margin = ggplot2::margin(
        t = 0,
        r = 6,
        b = legend_bottom_margin,
        l = 6
      )
    )
}


create_spatial_edger_panel_title_block <- function(
    panel_title,
    title_size = 28
) {

  ggplot2::ggplot() +
    ggplot2::labs(
      title = panel_title
    ) +
    ggplot2::theme_void(
      base_family = "DejaVu Sans"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = title_size,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(
          t = 0,
          b = -4
        )
      ),
      plot.margin = ggplot2::margin(
        t = 0,
        r = 3,
        b = -5,
        l = 3
      )
    )
}


create_spatial_edger_panel_subtitle_block <- function(
    panel_subtitle,
    subtitle_size = 18
) {

  ggplot2::ggplot() +
    ggplot2::labs(
      title = panel_subtitle
    ) +
    ggplot2::theme_void(
      base_family = "DejaVu Sans"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = subtitle_size,
        hjust = 0.5,
        lineheight = 1.05,
        margin = ggplot2::margin(
          t = 0,
          b = -6
        )
      ),
      plot.margin = ggplot2::margin(
        t = -2,
        r = 3,
        b = -8,
        l = 3
      )
    )
}


create_significant_gene_count_spatial_panel <- function(
    section_context,
    cluster_count_table,
    count_column,
    panel_title,
    panel_subtitle,
    palette_colors,
    cluster_base_colours,
    cluster_names = NULL,
    show_histology_image = FALSE,
    point_size_no_image = 3.00,
    point_size_with_image = 2.10,
    panel_border_linewidth = 1.5,
    legend_ncol = 2L,
    legend_text_size = 14.4,
    legend_point_size = 5.8,
    legend_row_spacing = 1.20,
    legend_bottom_margin = 8
) {

  colour_info <- map_spatial_edger_cluster_counts_to_colours(
    cluster_count_table = cluster_count_table,
    count_column = count_column,
    palette_colors = palette_colors
  )

  plot_data <- section_context$plot_data

  cluster_levels <- as.character(
    cluster_count_table$cluster_id
  )

  plot_data$cluster_id <- factor(
    as.character(plot_data$cluster_id),
    levels = cluster_levels
  )

  if (
    isTRUE(show_histology_image) &&
      !is.null(section_context$image_array)
  ) {

    spatial_plot <- ggplot2::ggplot() +
      ggplot2::annotation_raster(
        raster = section_context$image_array,
        xmin = 0,
        xmax = section_context$image_width,
        ymin = section_context$image_height,
        ymax = 0
      ) +
      ggplot2::geom_point(
        data = plot_data,
        ggplot2::aes(
          x = .data$x_plot,
          y = .data$y_plot,
          colour = .data$cluster_id
        ),
        size = point_size_with_image,
        shape = 16,
        stroke = 0
      )

  } else {

    spatial_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data$x_plot,
        y = .data$y_plot,
        colour = .data$cluster_id
      )
    ) +
      ggplot2::geom_point(
        size = point_size_no_image,
        shape = 16,
        stroke = 0
      )
  }

  spatial_plot <- spatial_plot +
    ggplot2::scale_colour_manual(
      values = colour_info$colours,
      limits = cluster_levels,
      breaks = cluster_levels,
      drop = FALSE,
      na.value = "#D9D9D9",
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      limits = section_context$frame_limits$x_limits,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_reverse(
      limits = rev(section_context$frame_limits$y_limits),
      expand = c(0, 0)
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = panel_title,
      subtitle = panel_subtitle
    ) +
    ggplot2::theme_void(
      base_family = "DejaVu Sans"
    ) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = panel_border_linewidth
      ),
      plot.title = ggplot2::element_text(
        size = 28,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(t = 0, b = 1)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 18,
        hjust = 0.5,
        lineheight = 1.0,
        margin = ggplot2::margin(t = 0, b = 2)
      ),
      plot.margin = ggplot2::margin(
        t = 0,
        r = 8,
        b = 8,
        l = 8
      )
    )

  dual_dot_legend <- create_spatial_edger_dual_dot_legend(
    cluster_count_table = cluster_count_table,
    count_column = count_column,
    result_colours = colour_info$colours,
    cluster_base_colours = cluster_base_colours,
    cluster_names = cluster_names,
    legend_ncol = legend_ncol,
    legend_text_size = legend_text_size,
    legend_point_size = legend_point_size,
    legend_row_spacing = legend_row_spacing,
    legend_bottom_margin = legend_bottom_margin
  )

  panel_patchwork <- (
    dual_dot_legend /
      spatial_plot
  ) +
    patchwork::plot_layout(
      ncol = 1L,
      heights = c(
        0.17,
        1
      )
    )

  # Convert the nested four-row patchwork into one grob before combining the
  # three biological panels. Without this step, patchwork flattens the nested
  # layouts and tries to place their internal plots in the outer 1 x 3 grid.
  panel_grob <- patchwork::patchworkGrob(
    panel_patchwork
  )

  patchwork::wrap_elements(
    full = panel_grob,
    clip = FALSE
  )
}


# ==============================================================================
# 6. Create and save one threshold variant
# ==============================================================================

create_main_group_significant_gene_count_three_panel_plot <- function(
    section_context,
    threshold_summary,
    clustering_parameters_label,
    green_palette_colors,
    red_palette_colors,
    blue_palette_colors,
    cluster_base_colours = NULL,
    cluster_names = NULL,
    show_histology_image = FALSE,
    point_size_no_image = 3.00,
    point_size_with_image = 2.10,
    legend_ncol = 2L,
    legend_row_spacing = 1.20,
    legend_bottom_margin = 8
) {

  cluster_count_table <- threshold_summary$cluster_count_table
  summary_table <- threshold_summary$summary_table

  fdr_threshold <- summary_table$fdr_threshold[[1]]
  abs_log2fc_threshold <- summary_table$abs_log2fc_threshold[[1]]

  total_results <- summary_table$total_significant_results[[1]]
  unique_genes <- summary_table$unique_significant_genes[[1]]
  up_results <- summary_table$up_significant_results[[1]]
  unique_up_genes <- summary_table$unique_up_genes[[1]]
  down_results <- summary_table$down_significant_results[[1]]
  unique_down_genes <- summary_table$unique_down_genes[[1]]

  total_panel <- create_significant_gene_count_spatial_panel(
    section_context = section_context,
    cluster_count_table = cluster_count_table,
    count_column = "n_total",
    panel_title = "Total significant genes",
    panel_subtitle = paste0(
      "UP + DOWN | ",
      format_spatial_edger_integer(total_results),
      " cluster-gene results | ",
      format_spatial_edger_integer(unique_genes),
      " unique genes"
    ),
    palette_colors = green_palette_colors,
    cluster_base_colours = cluster_base_colours,
    cluster_names = cluster_names,
    show_histology_image = show_histology_image,
    point_size_no_image = point_size_no_image,
    point_size_with_image = point_size_with_image,
    legend_ncol = legend_ncol,
    legend_row_spacing = legend_row_spacing,
    legend_bottom_margin = legend_bottom_margin
  )

  up_panel <- create_significant_gene_count_spatial_panel(
    section_context = section_context,
    cluster_count_table = cluster_count_table,
    count_column = "n_up",
    panel_title = "UP in ASD",
    panel_subtitle = paste0(
      format_spatial_edger_integer(up_results),
      " cluster-gene results | ",
      format_spatial_edger_integer(unique_up_genes),
      " unique genes"
    ),
    palette_colors = red_palette_colors,
    cluster_base_colours = cluster_base_colours,
    cluster_names = cluster_names,
    show_histology_image = show_histology_image,
    point_size_no_image = point_size_no_image,
    point_size_with_image = point_size_with_image,
    legend_ncol = legend_ncol,
    legend_row_spacing = legend_row_spacing,
    legend_bottom_margin = legend_bottom_margin
  )

  down_panel <- create_significant_gene_count_spatial_panel(
    section_context = section_context,
    cluster_count_table = cluster_count_table,
    count_column = "n_down",
    panel_title = "DOWN in ASD",
    panel_subtitle = paste0(
      format_spatial_edger_integer(down_results),
      " cluster-gene results | ",
      format_spatial_edger_integer(unique_down_genes),
      " unique genes"
    ),
    palette_colors = blue_palette_colors,
    cluster_base_colours = cluster_base_colours,
    cluster_names = cluster_names,
    show_histology_image = show_histology_image,
    point_size_no_image = point_size_no_image,
    point_size_with_image = point_size_with_image,
    legend_ncol = legend_ncol,
    legend_row_spacing = legend_row_spacing,
    legend_bottom_margin = legend_bottom_margin
  )

  plot_title <- paste0(
    "Main group effect (ASD vs Neurotypical): significant genes per spatial cluster"
  )

  plot_subtitle <- paste0(
    "Criteria: FDR < ",
    format_spatial_edger_number(fdr_threshold, digits = 2L),
    " and |log2FC| > ",
    format_spatial_edger_number(abs_log2fc_threshold, digits = 2L),
    " | all significant cluster-gene results: ",
    format_spatial_edger_integer(total_results),
    " | unique genes: ",
    format_spatial_edger_integer(unique_genes),
    "\nSample: ",
    section_context$sample_ID,
    " | tissue spots: ",
    format_spatial_edger_integer(section_context$n_spots),
    "\n",
    clustering_parameters_label
  )

  three_panel_grid <- patchwork::wrap_plots(
    plots = list(
      total_panel,
      up_panel,
      down_panel
    ),
    nrow = 1L,
    ncol = 3L
  )

  three_panel_grid +
    patchwork::plot_annotation(
      title = plot_title,
      subtitle = plot_subtitle,
      theme = ggplot2::theme(
        text = ggplot2::element_text(
          family = "DejaVu Sans"
        ),
        plot.title = ggplot2::element_text(
          size = 36,
          face = "bold",
          hjust = 0.5,
          margin = ggplot2::margin(b = 5)
        ),
        plot.subtitle = ggplot2::element_text(
          size = 20,
          hjust = 0.5,
          lineheight = 1.08,
          margin = ggplot2::margin(b = 8)
        ),
        plot.margin = ggplot2::margin(
          t = 8,
          r = 10,
          b = 8,
          l = 10
        )
      )
    )
}


get_main_group_significant_gene_count_output_files <- function(
    output_root_directory,
    sample_id_to_plot,
    fdr_threshold,
    abs_log2fc_threshold
) {

  threshold_label <- build_spatial_edger_threshold_label(
    fdr_threshold = fdr_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold
  )

  threshold_output_directory <- file.path(
    output_root_directory,
    threshold_label
  )

  safe_sample_id <- gsub(
    "[^A-Za-z0-9._-]+",
    "-",
    as.character(sample_id_to_plot)
  )

  file_prefix <- paste0(
    "mainEffectGroup_",
    threshold_label,
    "_sample_",
    safe_sample_id
  )

  list(
    threshold_label = threshold_label,
    threshold_output_directory = threshold_output_directory,
    plot_pdf = file.path(
      threshold_output_directory,
      paste0(
        "01_",
        file_prefix,
        "_significantGeneCountsPerCluster_threePanels.pdf"
      )
    ),
    cluster_count_tsv = file.path(
      threshold_output_directory,
      paste0(
        "02_",
        file_prefix,
        "_significantGeneCountsPerCluster.tsv"
      )
    ),
    threshold_summary_tsv = file.path(
      threshold_output_directory,
      paste0(
        "03_",
        file_prefix,
        "_thresholdSummary.tsv"
      )
    ),
    significant_results_tsv = file.path(
      threshold_output_directory,
      paste0(
        "04_",
        file_prefix,
        "_significantClusterGeneResults.tsv"
      )
    )
  )
}


save_main_group_significant_gene_count_variant <- function(
    section_context,
    main_group_results,
    output_root_directory,
    fdr_threshold,
    abs_log2fc_threshold,
    clustering_parameters_label,
    green_palette_colors,
    red_palette_colors,
    blue_palette_colors,
    cluster_base_colours = NULL,
    cluster_names = NULL,
    show_histology_image = FALSE,
    point_size_no_image = 3.00,
    point_size_with_image = 2.10,
    legend_ncol = 2L,
    legend_row_spacing = 1.20,
    legend_bottom_margin = 8,
    pdf_width = 27,
    pdf_height = 12.0,
    skip_completed = TRUE
) {

  output_files <- get_main_group_significant_gene_count_output_files(
    output_root_directory = output_root_directory,
    sample_id_to_plot = section_context$sample_ID,
    fdr_threshold = fdr_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold
  )

  expected_files <- unlist(
    output_files[
      c(
        "plot_pdf",
        "cluster_count_tsv",
        "threshold_summary_tsv",
        "significant_results_tsv"
      )
    ],
    use.names = FALSE
  )

  if (
    isTRUE(skip_completed) &&
      all(file.exists(expected_files)) &&
      all(file.info(expected_files)$size > 0)
  ) {
    message(
      "Skipping completed threshold variant: ",
      output_files$threshold_label
    )

    existing_summary <- read.delim(
      output_files$threshold_summary_tsv,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    return(
      invisible(
        list(
          output_files = output_files,
          summary_table = existing_summary,
          status = "skipped_existing"
        )
      )
    )
  }

  dir.create(
    output_files$threshold_output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  threshold_summary <- summarize_main_group_significant_genes_per_cluster(
    main_group_results = main_group_results,
    cluster_levels = section_context$cluster_levels,
    fdr_threshold = fdr_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold
  )

  threshold_summary$cluster_count_table <-
    threshold_summary$cluster_count_table |>
    dplyr::mutate(
      sample_ID = section_context$sample_ID,
      fdr_threshold = as.numeric(fdr_threshold),
      abs_log2fc_threshold = as.numeric(abs_log2fc_threshold),
      .before = 1
    )

  threshold_summary$summary_table <-
    threshold_summary$summary_table |>
    dplyr::mutate(
      sample_ID = section_context$sample_ID,
      threshold_label = output_files$threshold_label,
      .before = 1
    )

  three_panel_plot <- create_main_group_significant_gene_count_three_panel_plot(
    section_context = section_context,
    threshold_summary = threshold_summary,
    clustering_parameters_label = clustering_parameters_label,
    green_palette_colors = green_palette_colors,
    red_palette_colors = red_palette_colors,
    blue_palette_colors = blue_palette_colors,
    cluster_base_colours = cluster_base_colours,
    cluster_names = cluster_names,
    show_histology_image = show_histology_image,
    point_size_no_image = point_size_no_image,
    point_size_with_image = point_size_with_image,
    legend_ncol = legend_ncol,
    legend_row_spacing = legend_row_spacing,
    legend_bottom_margin = legend_bottom_margin
  )

  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  ggplot2::ggsave(
    filename = output_files$plot_pdf,
    plot = three_panel_plot,
    device = pdf_device,
    width = pdf_width,
    height = pdf_height,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )

  write_spatial_edger_tsv(
    threshold_summary$cluster_count_table,
    output_files$cluster_count_tsv
  )

  write_spatial_edger_tsv(
    threshold_summary$summary_table,
    output_files$threshold_summary_tsv
  )

  write_spatial_edger_tsv(
    threshold_summary$significant_results,
    output_files$significant_results_tsv
  )

  missing_files <- expected_files[
    !file.exists(expected_files)
  ]

  if (length(missing_files) > 0L) {
    stop(
      "Missing output file(s) for ",
      output_files$threshold_label,
      ":\n",
      paste(missing_files, collapse = "\n"),
      call. = FALSE
    )
  }

  empty_files <- expected_files[
    file.info(expected_files)$size == 0
  ]

  if (length(empty_files) > 0L) {
    stop(
      "Empty output file(s) for ",
      output_files$threshold_label,
      ":\n",
      paste(empty_files, collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(
    list(
      output_files = output_files,
      summary_table = threshold_summary$summary_table,
      status = "completed"
    )
  )
}


# ==============================================================================
# 7. Run all requested threshold variants
# ==============================================================================

run_main_group_significant_gene_count_spatial_workflow <- function(
    input_seurat_rdata_file,
    edger_results_rdata_file,
    output_root_directory,
    sample_id_to_plot,
    cluster_column,
    clustering_parameters_label,
    sample_id_column = "sample_ID",
    seurat_object_name = NULL,
    image_scale = "lowres",
    fdr_thresholds = c(0.10, 0.05, 0.01),
    abs_log2fc_thresholds = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
    green_palette_colors = c(
      "#D9D9D9",
      "#E5F5E0",
      "#A1D99B",
      "#74C476",
      "#31A354",
      "#006D2C",
      "#00441B"
    ),
    red_palette_colors = c(
      "#D9D9D9",
      "#FEE5D9",
      "#FCAE91",
      "#FB6A4A",
      "#DE2D26",
      "#A50F15",
      "#67000D"
    ),
    blue_palette_colors = c(
      "#D9D9D9",
      "#EFF3FF",
      "#BDD7E7",
      "#6BAED6",
      "#3182BD",
      "#08519C",
      "#08306B"
    ),
    cluster_base_colours = NULL,
    cluster_names = NULL,
    show_histology_image = FALSE,
    point_size_no_image = 3.00,
    point_size_with_image = 2.10,
    legend_ncol = 2L,
    legend_row_spacing = 1.20,
    legend_bottom_margin = 8,
    pdf_width = 27,
    pdf_height = 12.0,
    skip_completed_variants = TRUE,
    continue_after_variant_error = TRUE
) {

  dir.create(
    output_root_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  message("Loading Seurat object: ", input_seurat_rdata_file)

  seurat_loaded <- load_single_seurat_for_significant_gene_count_plot(
    input_seurat_rdata_file = input_seurat_rdata_file,
    seurat_object_name = seurat_object_name
  )

  message("Selected Seurat object: ", seurat_loaded$object_name)

  section_context <- prepare_selected_section_for_significant_gene_count_plot(
    seurat_object = seurat_loaded$object,
    sample_id_to_plot = sample_id_to_plot,
    cluster_column = cluster_column,
    sample_id_column = sample_id_column,
    image_scale = image_scale
  )

  message(
    "Selected section: ",
    section_context$sample_ID,
    " | image: ",
    section_context$image_name,
    " | spots: ",
    section_context$n_spots,
    " | clusters: ",
    length(section_context$cluster_levels)
  )

  rm(seurat_loaded)
  invisible(gc(verbose = FALSE))

  message("Loading edgeR results: ", edger_results_rdata_file)

  edger_loaded <-
    load_main_group_edger_results_for_significant_gene_count_plot(
      edger_results_rdata_file = edger_results_rdata_file
    )

  main_group_results <- edger_loaded$main_group_results

  expected_number_of_variants <-
    length(fdr_thresholds) * length(abs_log2fc_thresholds)

  message(
    "Threshold variants to process: ",
    expected_number_of_variants
  )

  run_status_rows <- list()
  overview_rows <- list()
  variant_index <- 0L

  run_status_file <- file.path(
    output_root_directory,
    paste0(
      "00_mainEffectGroup_significantGeneCountsPerCluster_",
      "runStatus.tsv"
    )
  )

  threshold_overview_file <- file.path(
    output_root_directory,
    paste0(
      "00_mainEffectGroup_significantGeneCountsPerCluster_",
      "thresholdOverview.tsv"
    )
  )

  for (fdr_threshold in fdr_thresholds) {
    for (abs_log2fc_threshold in abs_log2fc_thresholds) {

      variant_index <- variant_index + 1L

      threshold_label <- build_spatial_edger_threshold_label(
        fdr_threshold = fdr_threshold,
        abs_log2fc_threshold = abs_log2fc_threshold
      )

      variant_started_at <- Sys.time()

      message(
        "\n[",
        variant_index,
        "/",
        expected_number_of_variants,
        "] ",
        threshold_label,
        " | FDR < ",
        fdr_threshold,
        " | |log2FC| > ",
        abs_log2fc_threshold
      )

      variant_result <- tryCatch(
        {
          save_main_group_significant_gene_count_variant(
            section_context = section_context,
            main_group_results = main_group_results,
            output_root_directory = output_root_directory,
            fdr_threshold = fdr_threshold,
            abs_log2fc_threshold = abs_log2fc_threshold,
            clustering_parameters_label = clustering_parameters_label,
            green_palette_colors = green_palette_colors,
            red_palette_colors = red_palette_colors,
            blue_palette_colors = blue_palette_colors,
            cluster_base_colours = cluster_base_colours,
            cluster_names = cluster_names,
            show_histology_image = show_histology_image,
            point_size_no_image = point_size_no_image,
            point_size_with_image = point_size_with_image,
            legend_ncol = legend_ncol,
            legend_row_spacing = legend_row_spacing,
            legend_bottom_margin = legend_bottom_margin,
            pdf_width = pdf_width,
            pdf_height = pdf_height,
            skip_completed = skip_completed_variants
          )
        },
        error = function(error_condition) error_condition
      )

      variant_finished_at <- Sys.time()

      elapsed_seconds <- as.numeric(
        difftime(
          variant_finished_at,
          variant_started_at,
          units = "secs"
        )
      )

      if (inherits(variant_result, "error")) {

        current_status <- "failed"
        current_message <- conditionMessage(variant_result)
        current_plot_file <- NA_character_

        message(
          "FAILED: ",
          threshold_label,
          "\n",
          current_message
        )

      } else {

        current_status <- variant_result$status
        current_message <- ""
        current_plot_file <- variant_result$output_files$plot_pdf

        overview_rows[[length(overview_rows) + 1L]] <-
          variant_result$summary_table

        message(
          "Completed: ",
          threshold_label,
          " | plot: ",
          current_plot_file
        )
      }

      run_status_rows[[length(run_status_rows) + 1L]] <- data.frame(
        order = variant_index,
        total_variants = expected_number_of_variants,
        sample_ID = section_context$sample_ID,
        threshold_label = threshold_label,
        fdr_threshold = as.numeric(fdr_threshold),
        abs_log2fc_threshold = as.numeric(abs_log2fc_threshold),
        status = current_status,
        message = current_message,
        plot_file = current_plot_file,
        elapsed_seconds = round(elapsed_seconds, 2),
        started_at = format(
          variant_started_at,
          "%Y-%m-%d %H:%M:%S"
        ),
        finished_at = format(
          variant_finished_at,
          "%Y-%m-%d %H:%M:%S"
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      write_spatial_edger_tsv(
        do.call(rbind, run_status_rows),
        run_status_file
      )

      if (
        current_status == "failed" &&
          !isTRUE(continue_after_variant_error)
      ) {
        stop(current_message, call. = FALSE)
      }

      invisible(gc(verbose = FALSE))
    }
  }

  final_status_table <- do.call(
    rbind,
    run_status_rows
  )

  write_spatial_edger_tsv(
    final_status_table,
    run_status_file
  )

  if (length(overview_rows) > 0L) {

    threshold_overview <- dplyr::bind_rows(
      overview_rows
    ) |>
      dplyr::arrange(
        dplyr::desc(.data$fdr_threshold),
        .data$abs_log2fc_threshold
      )

    write_spatial_edger_tsv(
      threshold_overview,
      threshold_overview_file
    )
  } else {
    threshold_overview <- tibble::tibble()
  }

  message("\n============================================================")
  message("Finished significant-gene-count spatial visualization.")
  message("Sample: ", section_context$sample_ID)
  message("Completed/skipped: ", sum(final_status_table$status != "failed"))
  message("Failed: ", sum(final_status_table$status == "failed"))
  message("Output directory: ", output_root_directory)
  message("Run status: ", run_status_file)
  message("Threshold overview: ", threshold_overview_file)
  message("============================================================")

  invisible(
    list(
      section_context = section_context,
      run_status = final_status_table,
      threshold_overview = threshold_overview,
      run_status_file = run_status_file,
      threshold_overview_file = threshold_overview_file,
      output_root_directory = output_root_directory
    )
  )
}

# ==============================================================================
# End
# ==============================================================================
