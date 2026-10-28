# ==============================================================================
# Spatial gene-expression visualization
# ==============================================================================

sort_cluster_ids <- function(cluster_ids) {
  cluster_ids <- unique(as.character(cluster_ids))
  numeric_ids <- suppressWarnings(as.numeric(cluster_ids))
  if (!anyNA(numeric_ids)) {
    return(cluster_ids[order(numeric_ids)])
  }
  sort(cluster_ids)
}

standardize_spatial_coordinates <- function(coordinate_table) {
  coordinate_table <- as.data.frame(
    coordinate_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!"cell" %in% colnames(coordinate_table)) {
    coordinate_table$cell <- rownames(coordinate_table)
  }
  x_candidates <- c("imagecol", "x", "col", "pxl_col_in_fullres")
  y_candidates <- c("imagerow", "y", "row", "pxl_row_in_fullres")
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

compute_square_field_limits <- function(
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

get_available_cluster_ids <- function(
  seurat_object,
  cluster_column
) {
  metadata_table <- seurat_object[[]]
  if (!cluster_column %in% colnames(metadata_table)) {
    stop(
      "Missing Seurat metadata column: ",
      cluster_column,
      call. = FALSE
    )
  }
  sort_cluster_ids(metadata_table[[cluster_column]])
}

get_sample_layout_metadata <- function(
  seurat_object,
  sample_id_column,
  group_column = "fmt_donor_group",
  sex_column = "sex"
) {
  metadata_table <- seurat_object[[]]
  required_columns <- c(
    sample_id_column,
    group_column,
    sex_column
  )
  missing_columns <- setdiff(
    required_columns,
    colnames(metadata_table)
  )
  if (length(missing_columns) > 0L) {
    stop(
      "Missing Seurat metadata column(s) required for the four-column spatial layout: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  sample_metadata <- metadata_table |>
    tibble::rownames_to_column("cell") |>
    dplyr::transmute(
      sample_ID = as.character(.data[[sample_id_column]]),
      fmt_donor_group = as.character(.data[[group_column]]),
      sex = as.character(.data[[sex_column]])
    ) |>
    dplyr::distinct()
  conflicts <- sample_metadata |>
    dplyr::count(
      .data$sample_ID,
      name = "n_metadata_rows"
    ) |>
    dplyr::filter(.data$n_metadata_rows > 1L)
  if (nrow(conflicts) > 0L) {
    stop(
      "Conflicting Sex/Group metadata found for sample(s): ",
      paste(conflicts$sample_ID, collapse = ", "),
      call. = FALSE
    )
  }
  sample_metadata <- sample_metadata |>
    dplyr::mutate(
      sex = factor(
        .data$sex,
        levels = c("Male", "Female")
      ),
      fmt_donor_group = factor(
        .data$fmt_donor_group,
        levels = c("Neurotypical", "ASD")
      )
    )
  if (anyNA(sample_metadata$sex)) {
    stop(
      "Unexpected sex value in Seurat metadata. Expected Male/Female.",
      call. = FALSE
    )
  }
  if (anyNA(sample_metadata$fmt_donor_group)) {
    stop(
      "Unexpected FMT donor-group value in Seurat metadata. Expected Neurotypical/ASD.",
      call. = FALSE
    )
  }
  group_levels <- c(
    "Male Neurotypical",
    "Male ASD",
    "Female Neurotypical",
    "Female ASD"
  )
  sample_metadata |>
    dplyr::mutate(
      display_group = factor(
        paste(
          as.character(.data$sex),
          as.character(.data$fmt_donor_group)
        ),
        levels = group_levels
      )
    ) |>
    dplyr::arrange(
      .data$display_group,
      .data$sample_ID
    ) |>
    dplyr::mutate(
      sample_order = dplyr::row_number()
    )
}

get_ordered_sample_ids <- function(
  seurat_object,
  sample_id_column,
  group_column = "fmt_donor_group",
  sex_column = "sex"
) {
  get_sample_layout_metadata(
    seurat_object = seurat_object,
    sample_id_column = sample_id_column,
    group_column = group_column,
    sex_column = sex_column
  ) |>
    dplyr::pull(.data$sample_ID)
}

build_image_to_sample_mapping <- function(
  seurat_object,
  sample_id_column
) {
  metadata_table <- seurat_object[[]]
  image_names <- SeuratObject::Images(seurat_object)
  if (length(image_names) == 0L) {
    stop(
      "No spatial images are present in the Seurat object.",
      call. = FALSE
    )
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
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )
  dplyr::bind_rows(image_mapping_rows)
}

prepare_spatial_section <- function(
  seurat_object,
  sample_id,
  sample_id_column,
  cluster_column,
  image_mapping,
  image_scale = "lowres",
  panel_padding_fraction = 0.03
) {
  metadata_table <- seurat_object[[]]
  sample_id <- as.character(sample_id)[[1]]
  matching_images <- image_mapping$image_name[
    image_mapping$sample_ID == sample_id
  ]
  if (length(matching_images) != 1L) {
    stop(
      "Expected exactly one spatial image for sample ",
      sample_id,
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
  coordinate_table <- standardize_spatial_coordinates(
    coordinate_table
  )
  common_cells <- image_cells[
    image_cells %in% coordinate_table$cell &
      image_cells %in% rownames(metadata_table)
  ]
  if (length(common_cells) == 0L) {
    stop(
      "No shared cells among image, coordinates and metadata for sample ",
      sample_id,
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
  coordinate_table$sample_ID <- sample_id
  coordinate_table$cluster_id <- as.character(
    spot_metadata[[cluster_column]]
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
    image_height <- max(
      coordinate_table$y_plot,
      na.rm = TRUE
    )
    image_width <- max(
      coordinate_table$x_plot,
      na.rm = TRUE
    )
  }
  frame_limits <- compute_square_field_limits(
    xmin = 0,
    xmax = image_width,
    ymin = 0,
    ymax = image_height,
    padding_fraction = panel_padding_fraction
  )
  list(
    sample_ID = sample_id,
    image_name = image_name,
    image_array = raw_image,
    image_height = image_height,
    image_width = image_width,
    frame_limits = frame_limits,
    cells = common_cells,
    plot_data = coordinate_table
  )
}

prepare_all_spatial_sections <- function(
  seurat_object,
  sample_ids,
  sample_id_column,
  cluster_column,
  image_scale = "lowres"
) {
  image_mapping <- build_image_to_sample_mapping(
    seurat_object = seurat_object,
    sample_id_column = sample_id_column
  )
  section_contexts <- lapply(
    sample_ids,
    function(sample_id) {
      prepare_spatial_section(
        seurat_object = seurat_object,
        sample_id = sample_id,
        sample_id_column = sample_id_column,
        cluster_column = cluster_column,
        image_mapping = image_mapping,
        image_scale = image_scale
      )
    }
  )
  names(section_contexts) <- sample_ids
  section_contexts
}

add_gene_expression_to_section <- function(
  section_context,
  seurat_object,
  gene,
  expression_assay,
  expression_layer_prefix = "data.",
  counts_layer_prefix = "counts."
) {
  expression_layer <- paste0(
    expression_layer_prefix,
    section_context$sample_ID
  )
  available_layers <- SeuratObject::Layers(
    seurat_object[[expression_assay]]
  )
  if (!expression_layer %in% available_layers) {
    stop(
      "Expression layer is absent for sample ",
      section_context$sample_ID,
      ": ",
      expression_layer,
      call. = FALSE
    )
  }
  expression_matrix <- SeuratObject::LayerData(
    object = seurat_object[[expression_assay]],
    layer = expression_layer,
    features = gene,
    cells = section_context$cells,
    fast = FALSE
  )
  if (nrow(expression_matrix) != 1L) {
    stop(
      "Expected exactly one expression row for gene ",
      gene,
      " in layer ",
      expression_layer,
      ".",
      call. = FALSE
    )
  }
  missing_cells <- setdiff(
    section_context$cells,
    colnames(expression_matrix)
  )
  if (length(missing_cells) > 0L) {
    stop(
      "Some tissue spots are absent from expression layer ",
      expression_layer,
      ": ",
      length(missing_cells),
      call. = FALSE
    )
  }

  counts_layer <- paste0(
    counts_layer_prefix,
    section_context$sample_ID
  )

  if (!counts_layer %in% available_layers) {
    stop(
      "Counts layer is absent for sample ",
      section_context$sample_ID,
      ": ",
      counts_layer,
      call. = FALSE
    )
  }

  counts_matrix <- SeuratObject::LayerData(
    object = seurat_object[[expression_assay]],
    layer = counts_layer,
    features = gene,
    cells = section_context$cells,
    fast = FALSE
  )

  if (nrow(counts_matrix) != 1L) {
    stop(
      "Expected exactly one counts row for gene ",
      gene,
      " in layer ",
      counts_layer,
      ".",
      call. = FALSE
    )
  }

  counts_values <- as.numeric(
    counts_matrix[
      1,
      section_context$cells,
      drop = TRUE
    ]
  )

  percent_positive <- 100 * mean(
    counts_values > 0,
    na.rm = TRUE
  )

  plot_data <- section_context$plot_data
  plot_data$expression <- as.numeric(
    expression_matrix[
      1,
      section_context$cells,
      drop = TRUE
    ]
  )
  plot_data$percent_positive <- percent_positive
  plot_data
}

add_gene_expression_to_all_sections <- function(
  section_contexts,
  seurat_object,
  gene,
  expression_assay,
  expression_layer_prefix = "data.",
  counts_layer_prefix = "counts.",
  sample_order = names(section_contexts)
) {
  plot_data_list <- lapply(
    section_contexts,
    function(section_context) {
      add_gene_expression_to_section(
        section_context = section_context,
        seurat_object = seurat_object,
        gene = gene,
        expression_assay = expression_assay,
        expression_layer_prefix = expression_layer_prefix,
        counts_layer_prefix = counts_layer_prefix
      )
    }
  )
  plot_data <- dplyr::bind_rows(plot_data_list)
  plot_data$sample_ID <- factor(
    plot_data$sample_ID,
    levels = sample_order
  )
  plot_data
}

create_group_header_plot <- function(group_label) {
  ggplot2::ggplot() +
    ggplot2::labs(title = group_label) +
    ggplot2::theme_void(base_family = "DejaVu Sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(t = 7, b = 2)
      ),
      plot.margin = ggplot2::margin(4, 2, 0, 2)
    )
}

create_single_gene_spatial_panel <- function(
  section_context,
  section_plot_data,
  indicated_clusters,
  red_palette,
  green_palette,
  expression_limits,
  point_size,
  gene
) {
  section_plot_data$indicated <-
    as.character(section_plot_data$cluster_id) %in%
    as.character(indicated_clusters)

  indicated_data <- section_plot_data[
    section_plot_data$indicated,
    ,
    drop = FALSE
  ]

  remaining_data <- section_plot_data[
    !section_plot_data$indicated,
    ,
    drop = FALSE
  ]

  p <- ggplot2::ggplot()

  # Remaining, non-indicated expression is shown in green.
  if (nrow(remaining_data) > 0L) {
    p <- p +
      ggplot2::geom_point(
        data = remaining_data,
        ggplot2::aes(
          x = .data$x_plot,
          y = .data$y_plot,
          colour = .data$expression
        ),
        size = point_size,
        shape = 16,
        stroke = 0
      ) +
      ggplot2::scale_colour_gradientn(
        colours = green_palette,
        limits = expression_limits,
        oob = scales::squish,
        name = paste0(gene, " — Remaining expression"),
        guide = ggplot2::guide_colorbar(
          direction = "horizontal",
          title.position = "top",
          title.hjust = 0.5,
          label.position = "bottom",
          barwidth = grid::unit(54, "mm"),
          barheight = grid::unit(5.6, "mm")
        )
      )
  }

  if (nrow(remaining_data) > 0L && nrow(indicated_data) > 0L) {
    p <- p + ggnewscale::new_scale_colour()
  }

  # Indicated clusters are shown in red.
  if (nrow(indicated_data) > 0L) {
    p <- p +
      ggplot2::geom_point(
        data = indicated_data,
        ggplot2::aes(
          x = .data$x_plot,
          y = .data$y_plot,
          colour = .data$expression
        ),
        size = point_size,
        shape = 16,
        stroke = 0
      ) +
      ggplot2::scale_colour_gradientn(
        colours = red_palette,
        limits = expression_limits,
        oob = scales::squish,
        name = paste0(gene, " — Indicated clusters"),
        guide = ggplot2::guide_colorbar(
          direction = "horizontal",
          title.position = "top",
          title.hjust = 0.5,
          label.position = "bottom",
          barwidth = grid::unit(54, "mm"),
          barheight = grid::unit(5.6, "mm")
        )
      )
  }

  percent_positive_values <- unique(
    section_plot_data$percent_positive[
      is.finite(section_plot_data$percent_positive)
    ]
  )

  percent_positive <- if (length(percent_positive_values) > 0L) {
    percent_positive_values[[1]]
  } else {
    NA_real_
  }

  sample_title <- if (is.finite(percent_positive)) {
    paste0(
      section_context$sample_ID,
      " (",
      formatC(
        percent_positive,
        format = "f",
        digits = 1
      ),
      "%)"
    )
  } else {
    section_context$sample_ID
  }

  p +
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
      title = sample_title
    ) +
    ggplot2::theme_void(base_family = "DejaVu Sans") +
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
        linewidth = 0.65
      ),
      plot.title = ggplot2::element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 7)
      ),
      plot.margin = ggplot2::margin(
        t = 11,
        r = 1,
        b = 11,
        l = 1
      ),
      legend.title = ggplot2::element_text(
        size = 14,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 14
      )
    )
}

create_gene_spatial_plot_all_sections <- function(
  plot_data,
  section_contexts,
  sample_layout,
  gene,
  indicated_clusters = character(0),
  red_palette,
  green_palette,
  point_size = 1.4
) {
  indicated_clusters <- as.character(indicated_clusters)

  finite_expression <- plot_data$expression[
    is.finite(plot_data$expression)
  ]

  if (length(finite_expression) == 0L) {
    expression_limits <- c(0, 1)
  } else {
    expression_max <- max(
      finite_expression,
      na.rm = TRUE
    )

    if (!is.finite(expression_max) || expression_max <= 0) {
      expression_max <- 1
    }

    expression_limits <- c(
      0,
      expression_max
    )
  }

  group_levels <- levels(
    sample_layout$display_group
  )

  samples_by_group <- lapply(
    group_levels,
    function(group_label) {
      as.character(
        sample_layout$sample_ID[
          sample_layout$display_group == group_label
        ]
      )
    }
  )

  names(samples_by_group) <- group_levels

  max_samples_per_group <- max(
    vapply(
      samples_by_group,
      length,
      FUN.VALUE = integer(1)
    )
  )

  # The patchwork design has a dedicated full-width guide row above the
  # four group headings. This keeps the legend from consuming panel width.
  plot_symbols <- c(
    LETTERS,
    letters,
    as.character(0:9)
  )

  number_of_plots <-
    1L +
    4L +
    4L * max_samples_per_group

  if (number_of_plots > length(plot_symbols)) {
    stop(
      "Too many spatial panels for the patchwork design.",
      call. = FALSE
    )
  }

  symbols <- plot_symbols[
    seq_len(number_of_plots)
  ]

  legend_symbol <- symbols[[1]]
  group_symbols <- symbols[2:5]
  sample_symbols <- symbols[
    6:number_of_plots
  ]

  design_rows <- c(
    paste0(
      rep(legend_symbol, 4L),
      collapse = ""
    ),
    paste0(
      group_symbols,
      collapse = ""
    )
  )

  for (row_index in seq_len(max_samples_per_group)) {
    current_symbols <- sample_symbols[
      ((row_index - 1L) * 4L + 1L):
        (row_index * 4L)
    ]

    design_rows <- c(
      design_rows,
      paste0(
        current_symbols,
        collapse = ""
      )
    )
  }

  design <- paste(
    design_rows,
    collapse = "\n"
  )

  layout_plots <- list()

  layout_plots[[legend_symbol]] <- patchwork::guide_area()

  for (group_index in seq_along(group_levels)) {
    layout_plots[[group_symbols[[group_index]]]] <-
      create_group_header_plot(
        group_levels[[group_index]]
      )
  }

  sample_plot_index <- 0L

  for (row_index in seq_len(max_samples_per_group)) {
    for (group_label in group_levels) {
      sample_plot_index <- sample_plot_index + 1L
      current_symbol <- sample_symbols[[sample_plot_index]]
      group_samples <- samples_by_group[[group_label]]

      if (row_index <= length(group_samples)) {
        sample_id <- group_samples[[row_index]]

        section_plot_data <- plot_data[
          as.character(plot_data$sample_ID) == sample_id,
          ,
          drop = FALSE
        ]

        layout_plots[[current_symbol]] <-
          create_single_gene_spatial_panel(
            section_context = section_contexts[[sample_id]],
            section_plot_data = section_plot_data,
            indicated_clusters = indicated_clusters,
            red_palette = red_palette,
            green_palette = green_palette,
            expression_limits = expression_limits,
            point_size = point_size,
            gene = gene
          )
      } else {
        layout_plots[[current_symbol]] <-
          patchwork::plot_spacer()
      }
    }
  }

  patchwork::wrap_plots(
    layout_plots,
    design = design,
    guides = "collect",
    widths = rep(1, 4),
    heights = c(
      0.24,
      0.10,
      rep(1, max_samples_per_group)
    )
  ) &
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.margin = ggplot2::margin(0, 0, 10, 0),
      legend.box.spacing = grid::unit(2, "mm"),
      legend.title = ggplot2::element_text(
        size = 14,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 14
      )
    )
}
