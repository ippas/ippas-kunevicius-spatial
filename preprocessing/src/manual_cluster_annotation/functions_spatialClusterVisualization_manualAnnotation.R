# ==============================================================================
# functions_spatialClusterVisualization_manualAnnotation.R
#
# Purpose:
# Sandbox copy used during manual cluster annotation.
# Create multi-panel spatial cluster plots from a Seurat object that already
# contains Visium images and clustering metadata.
#
# Main features:
# - validates clustering columns and spatial images;
# - builds a sample_ID <-> image_name mapping automatically;
# - creates a global cluster summary for the selected clustering column;
# - adds summary information directly into legend labels:
#     cluster_id (n=1234; 3.1%; mean±SD=77.1±12.4)
# - generates one multi-panel plot (for example, 16 samples in a 4x4 layout);
# - automatically adds the number of clusters to the main plot title;
# - draws a black border around every sample panel;
# - places the shared legend in a dedicated row under the main title and
#   above all sample titles;
# - uses a one-line legend title: <cluster_column> (global n, %, mean±SD);
# - enlarges legend points;
# - provides two neutral working palettes supporting up to 256 clusters;
#   the first 30 working colours are unchanged from the development workflow:
#     working30 and dark3;
# - saves PNG, PDF and summary tables.
# ==============================================================================


format_cluster_percentage <- function(x, digits = 1L) {

  formatC(
    as.numeric(x),
    format = "f",
    digits = as.integer(digits)
  )
}


format_cluster_decimal <- function(x, digits = 1L) {

  formatC(
    as.numeric(x),
    format = "f",
    digits = as.integer(digits)
  )
}


parse_cluster_column_name <- function(cluster_column) {

  if (
    !is.character(cluster_column) ||
      length(cluster_column) != 1L ||
      is.na(cluster_column) ||
      cluster_column == ""
  ) {
    stop("`cluster_column` must be one non-empty character value.", call. = FALSE)
  }

  parsed <- regexec(
    pattern = "^(.*)_res([0-9]+(?:\\.[0-9]+)?)$",
    text = cluster_column
  )

  parsed_match <- regmatches(cluster_column, parsed)[[1]]

  if (length(parsed_match) != 3L) {
    stop(
      "Could not parse clustering column name: ",
      cluster_column,
      ". Expected format similar to `leiden_res0.20`.",
      call. = FALSE
    )
  }

  list(
    algorithm = parsed_match[[2]],
    resolution = parsed_match[[3]]
  )
}


build_spatial_cluster_output_prefix <- function(
    analysis_prefix,
    cluster_column
) {

  parsed <- parse_cluster_column_name(cluster_column)

  paste0(
    analysis_prefix,
    "_",
    parsed$algorithm,
    "Res",
    parsed$resolution
  )
}


create_spatial_image_map <- function(seurat_object) {

  if (!inherits(seurat_object, "Seurat")) {
    stop("`seurat_object` must be a Seurat object.", call. = FALSE)
  }

  if (!"sample_ID" %in% colnames(seurat_object[[]])) {
    stop(
      "The Seurat object does not contain `sample_ID` metadata.",
      call. = FALSE
    )
  }

  image_names <- SeuratObject::Images(seurat_object)

  if (length(image_names) == 0L) {
    stop(
      "No spatial images were found in the Seurat object.",
      call. = FALSE
    )
  }

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
          "` contains zero cells.",
          call. = FALSE
        )
      }

      image_sample_ids <- unique(
        seurat_object[[]][image_cells, "sample_ID", drop = TRUE]
      )

      if (length(image_sample_ids) != 1L) {
        stop(
          "Spatial image `",
          image_name,
          "` is associated with ",
          length(image_sample_ids),
          " sample IDs. Exactly one was expected.",
          call. = FALSE
        )
      }

      data.frame(
        sample_ID = image_sample_ids[[1]],
        image_name = image_name,
        n_spots = length(image_cells),
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
    duplicated_samples <- unique(
      image_map$sample_ID[duplicated(image_map$sample_ID)]
    )

    stop(
      "Multiple spatial images were mapped to the same sample_ID: ",
      paste(duplicated_samples, collapse = ", "),
      call. = FALSE
    )
  }

  image_map
}


build_cluster_summary_table <- function(
    seurat_object,
    cluster_column,
    sample_order = NULL
) {

  if (!cluster_column %in% colnames(seurat_object[[]])) {
    stop(
      "Clustering column is absent from the Seurat object: ",
      cluster_column,
      call. = FALSE
    )
  }

  if (!"sample_ID" %in% colnames(seurat_object[[]])) {
    stop(
      "The Seurat object does not contain `sample_ID` metadata.",
      call. = FALSE
    )
  }

  metadata_table <- seurat_object[[]]

  if (is.null(sample_order)) {
    sample_order <- unique(as.character(metadata_table$sample_ID))
  }

  metadata_table <- metadata_table[
    metadata_table$sample_ID %in% sample_order,
    ,
    drop = FALSE
  ]

  cluster_values <- metadata_table[[cluster_column]]

  if (anyNA(cluster_values)) {
    stop(
      "NA values were found in clustering column `",
      cluster_column,
      "`.",
      call. = FALSE
    )
  }

  cluster_values <- as.character(cluster_values)
  sample_values <- as.character(metadata_table$sample_ID)

  cluster_counts <- sort(
    table(cluster_values),
    decreasing = FALSE
  )

  summary_table <- data.frame(
    cluster_id = names(cluster_counts),
    n_spots = as.integer(cluster_counts),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  all_numeric_clusters <- all(
    grepl("^-?[0-9]+$", summary_table$cluster_id)
  )

  if (all_numeric_clusters) {
    summary_table <- summary_table[
      order(as.integer(summary_table$cluster_id)),
      ,
      drop = FALSE
    ]
  } else {
    summary_table <- summary_table[
      order(summary_table$cluster_id),
      ,
      drop = FALSE
    ]
  }

  summary_table$spot_fraction <- summary_table$n_spots / sum(summary_table$n_spots)
  summary_table$spot_percentage <- 100 * summary_table$spot_fraction

  sample_levels <- sample_order
  cluster_levels <- summary_table$cluster_id

  count_table_by_sample <- table(
    factor(sample_values, levels = sample_levels),
    factor(cluster_values, levels = cluster_levels)
  )

  summary_table$mean_spots_per_sample <- vapply(
    seq_len(ncol(count_table_by_sample)),
    function(i) mean(count_table_by_sample[, i]),
    numeric(1)
  )

  summary_table$sd_spots_per_sample <- vapply(
    seq_len(ncol(count_table_by_sample)),
    function(i) stats::sd(count_table_by_sample[, i]),
    numeric(1)
  )

  summary_table$legend_label <- paste0(
    summary_table$cluster_id,
    " (n=",
    format(summary_table$n_spots, big.mark = " ", scientific = FALSE),
    "; ",
    format_cluster_percentage(summary_table$spot_percentage, digits = 1L),
    "%; mean\u00b1SD=",
    format_cluster_decimal(summary_table$mean_spots_per_sample, digits = 1L),
    "\u00b1",
    format_cluster_decimal(summary_table$sd_spots_per_sample, digits = 1L),
    ")"
  )

  rownames(summary_table) <- NULL

  cluster_counts_by_sample <- as.data.frame.matrix(count_table_by_sample)
  cluster_counts_by_sample$sample_ID <- rownames(cluster_counts_by_sample)
  rownames(cluster_counts_by_sample) <- NULL

  cluster_counts_by_sample <- cluster_counts_by_sample[
    ,
    c("sample_ID", cluster_levels),
    drop = FALSE
  ]

  list(
    cluster_summary = summary_table,
    cluster_counts_by_sample = cluster_counts_by_sample
  )
}


add_display_cluster_column <- function(
    seurat_object,
    cluster_column,
    sample_order = NULL,
    display_column = NULL
) {

  if (is.null(display_column)) {
    display_column <- paste0(
      cluster_column,
      "__display"
    )
  }

  summary_result <- build_cluster_summary_table(
    seurat_object = seurat_object,
    cluster_column = cluster_column,
    sample_order = sample_order
  )

  summary_table <- summary_result$cluster_summary

  label_map <- setNames(
    object = summary_table$legend_label,
    nm = summary_table$cluster_id
  )

  cluster_values <- as.character(
    seurat_object[[]][[cluster_column]]
  )

  display_values <- unname(
    label_map[cluster_values]
  )

  if (anyNA(display_values)) {
    stop(
      "Could not map all cluster values to legend labels for column `",
      cluster_column,
      "`.",
      call. = FALSE
    )
  }

  seurat_object[[display_column]] <- factor(
    display_values,
    levels = summary_table$legend_label
  )

  list(
    seurat_object = seurat_object,
    display_column = display_column,
    cluster_summary = summary_table,
    cluster_counts_by_sample = summary_result$cluster_counts_by_sample
  )
}


validate_spatial_cluster_plot_inputs <- function(
    seurat_object,
    cluster_column,
    sample_order = NULL
) {

  if (!inherits(seurat_object, "Seurat")) {
    stop("`seurat_object` must be a Seurat object.", call. = FALSE)
  }

  if (!cluster_column %in% colnames(seurat_object[[]])) {
    stop(
      "Missing clustering column: ",
      cluster_column,
      call. = FALSE
    )
  }

  if (!"sample_ID" %in% colnames(seurat_object[[]])) {
    stop(
      "The Seurat object does not contain `sample_ID` metadata.",
      call. = FALSE
    )
  }

  image_map <- create_spatial_image_map(
    seurat_object = seurat_object
  )

  if (is.null(sample_order)) {
    sample_order <- unique(
      as.character(seurat_object[[]]$sample_ID)
    )
  }

  missing_samples <- setdiff(
    sample_order,
    image_map$sample_ID
  )

  if (length(missing_samples) > 0L) {
    stop(
      "The following requested sample(s) do not have spatial images: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  validation_summary <- data.frame(
    n_samples_requested = length(sample_order),
    n_images_available = nrow(image_map),
    n_clusters = length(unique(seurat_object[[]][[cluster_column]])),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  list(
    image_map = image_map,
    sample_order = sample_order,
    validation_summary = validation_summary
  )
}


generate_cluster_palette <- function(
    n_clusters,
    palette_name = c(
      "working30",
      "dark3"
    )
) {

  palette_name <- match.arg(
    palette_name
  )

  if (
    !is.numeric(n_clusters) ||
      length(n_clusters) != 1L ||
      is.na(n_clusters) ||
      n_clusters < 1L
  ) {
    stop(
      "`n_clusters` must be one positive integer.",
      call. = FALSE
    )
  }

  n_clusters <- as.integer(
    n_clusters
  )

  # The first 30 colours remain exactly the same as in the development
  # workflow. Additional deterministic HCL colours are appended only when a
  # high-resolution clustering contains more than 30 clusters.
  maximum_supported_clusters <- 256L

  if (n_clusters > maximum_supported_clusters) {
    stop(
      "The current working palettes support at most ",
      maximum_supported_clusters,
      " clusters. Requested: ",
      n_clusters,
      ".",
      call. = FALSE
    )
  }

  working30_base <- c(
    "#E41A1C",
    "#377EB8",
    "#4DAF4A",
    "#984EA3",
    "#FF7F00",
    "#A65628",
    "#F781BF",
    "#00A6D6",
    "#7A8B00",
    "#D7308F",
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#66A61E",
    "#E6AB02",
    "#A6761D",
    "#1F78B4",
    "#33A02C",
    "#6A3D9A",
    "#B15928",
    "#17BECF",
    "#BCBD22",
    "#E377C2",
    "#8C564B",
    "#2CA02C",
    "#FF9896",
    "#9467BD",
    "#C49C94",
    "#DBDB8D",
    "#7F7F7F"
  )

  number_of_extra_colours <-
    maximum_supported_clusters -
    length(working30_base)

  golden_angle <- 137.507764

  extra_colour_indices <- seq_len(
    number_of_extra_colours
  )

  extra_colours <- grDevices::hcl(
    h = (
      extra_colour_indices *
        golden_angle
    ) %% 360,
    c = rep(
      c(
        85,
        65,
        75,
        55
      ),
      length.out =
        number_of_extra_colours
    ),
    l = rep(
      c(
        52,
        68,
        42,
        78
      ),
      length.out =
        number_of_extra_colours
    ),
    fixup = TRUE
  )

  working_extended <- c(
    working30_base,
    extra_colours
  )

  # Generate the complete fixed palette first and subset it afterwards. This
  # keeps cluster colours stable between solutions containing different numbers
  # of clusters.
  dark3_extended <- grDevices::hcl.colors(
    n = maximum_supported_clusters,
    palette = "Dark 3"
  )

  complete_palette <- switch(
    palette_name,
    working30 = working_extended,
    dark3 = dark3_extended
  )

  complete_palette[
    seq_len(n_clusters)
  ]
}


extract_shared_plot_legend <- function(plot_object) {

  plot_gtable <- ggplot2::ggplotGrob(plot_object)

  guide_indices <- which(
    grepl(
      "^guide-box",
      plot_gtable$layout$name
    )
  )

  if (length(guide_indices) == 0L) {
    stop(
      "Could not find a legend in the supplied ggplot object.",
      call. = FALSE
    )
  }

  non_empty_guide_indices <- guide_indices[
    !vapply(
      plot_gtable$grobs[guide_indices],
      inherits,
      logical(1),
      what = "zeroGrob"
    )
  ]

  if (length(non_empty_guide_indices) == 0L) {
    stop(
      "All detected guide boxes were empty.",
      call. = FALSE
    )
  }

  plot_gtable$grobs[[non_empty_guide_indices[[1]]]]
}


plot_spatial_clusters_all_samples <- function(
    seurat_object,
    cluster_column,
    sample_order = NULL,
    ncol = 4,
    plot_title = NULL,
    plot_subtitle = NULL,
    include_cluster_count_in_title = TRUE,
    show_image = FALSE,
    image_alpha = 1,
    crop = FALSE,
    pt.size.factor = 1.8,
    legend_position = "top",
    legend_ncol = 4,
    legend_point_size = 6,
    legend_height_ratio = 0.16,
    palette_name = c(
      "working30",
      "dark3"
    ),
    custom_cluster_colors = NULL,
    output_dir = NULL,
    output_prefix = NULL,
    save_png = TRUE,
    save_pdf = TRUE,
    png_width_in = 18,
    png_height_in = 22,
    pdf_width_in = 18,
    pdf_height_in = 22,
    dpi = 300,
    verbose = TRUE
) {

  palette_name <- match.arg(
    palette_name
  )

  required_packages <- c(
    "Seurat",
    "SeuratObject",
    "ggplot2",
    "patchwork",
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

  if (length(missing_packages) > 0L) {
    stop(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  input_check <- validate_spatial_cluster_plot_inputs(
    seurat_object = seurat_object,
    cluster_column = cluster_column,
    sample_order = sample_order
  )

  image_map <- input_check$image_map
  sample_order <- input_check$sample_order

  display_result <- add_display_cluster_column(
    seurat_object = seurat_object,
    cluster_column = cluster_column,
    sample_order = sample_order
  )

  plotting_object <- display_result$seurat_object
  display_column <- display_result$display_column
  cluster_summary <- display_result$cluster_summary
  cluster_counts_by_sample <- display_result$cluster_counts_by_sample

  cluster_levels <- cluster_summary$legend_label
  cluster_ids <- as.character(cluster_summary$cluster_id)

  default_cluster_colors <- generate_cluster_palette(
    n_clusters = length(cluster_levels),
    palette_name = palette_name
  )

  if (is.null(custom_cluster_colors)) {

    cluster_colors_by_id <- setNames(
      object = default_cluster_colors,
      nm = cluster_ids
    )

  } else {

    if (
      !is.character(custom_cluster_colors) ||
        is.null(names(custom_cluster_colors)) ||
        anyNA(custom_cluster_colors) ||
        anyNA(names(custom_cluster_colors)) ||
        any(!nzchar(custom_cluster_colors)) ||
        any(!nzchar(names(custom_cluster_colors)))
    ) {
      stop(
        "`custom_cluster_colors` must be a named character vector, ",
        "for example `c(\"1\" = \"#CBC9E2\")`.",
        call. = FALSE
      )
    }

    if (anyDuplicated(names(custom_cluster_colors)) > 0L) {
      stop(
        "Duplicated cluster IDs were found in `custom_cluster_colors`: ",
        paste(
          unique(
            names(custom_cluster_colors)[
              duplicated(names(custom_cluster_colors))
            ]
          ),
          collapse = ", "
        ),
        call. = FALSE
      )
    }

    invalid_colours <- custom_cluster_colors[
      !vapply(
        custom_cluster_colors,
        function(colour_value) {
          tryCatch(
            {
              grDevices::col2rgb(colour_value)
              TRUE
            },
            error = function(error_condition) FALSE
          )
        },
        logical(1)
      )
    ]

    if (length(invalid_colours) > 0L) {
      stop(
        "Invalid colour value(s) in `custom_cluster_colors`: ",
        paste(unique(invalid_colours), collapse = ", "),
        call. = FALSE
      )
    }

    cluster_colors_by_id <- setNames(
      object = default_cluster_colors,
      nm = cluster_ids
    )

    matching_cluster_ids <- intersect(
      cluster_ids,
      names(custom_cluster_colors)
    )

    cluster_colors_by_id[matching_cluster_ids] <-
      custom_cluster_colors[matching_cluster_ids]

    missing_custom_cluster_ids <- setdiff(
      cluster_ids,
      names(custom_cluster_colors)
    )

    if (
      length(missing_custom_cluster_ids) > 0L &&
        isTRUE(verbose)
    ) {
      warning(
        "No custom colour was supplied for cluster(s): ",
        paste(missing_custom_cluster_ids, collapse = ", "),
        ". Default palette colours will be used for these clusters.",
        call. = FALSE
      )
    }

    unused_custom_cluster_ids <- setdiff(
      names(custom_cluster_colors),
      cluster_ids
    )

    if (
      length(unused_custom_cluster_ids) > 0L &&
        isTRUE(verbose)
    ) {
      warning(
        "Custom colours were supplied for cluster ID(s) absent from `",
        cluster_column,
        "`: ",
        paste(unused_custom_cluster_ids, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  cluster_palette <- setNames(
    object = unname(cluster_colors_by_id[cluster_ids]),
    nm = cluster_levels
  )

  if (is.null(plot_title)) {
    plot_title <- paste0(
      "Spatial clustering: ",
      cluster_column
    )
  }

  if (isTRUE(include_cluster_count_in_title)) {
    plot_title <- paste0(
      plot_title,
      " | ",
      nrow(cluster_summary),
      " clusters"
    )
  }

  if (is.null(plot_subtitle)) {
    plot_subtitle <- paste0(
      "Panels: ",
      length(sample_order),
      " samples | legend shows global n, %, mean\u00b1SD across samples"
    )
  }

  if (verbose) {
    message(
      "Creating spatial cluster plots for ",
      length(sample_order),
      " sample(s) using clustering column `",
      cluster_column,
      "`."
    )
  }

  sample_plots <- lapply(
    sample_order,
    function(sample_id) {

      current_image_name <- image_map$image_name[
        match(sample_id, image_map$sample_ID)
      ]

      sample_cells <- colnames(plotting_object)[
        plotting_object[[]]$sample_ID == sample_id
      ]

      sample_object <- subset(
        x = plotting_object,
        cells = sample_cells
      )

      sample_group <- unique(
        sample_object[[]]$fmt_donor_group
      )

      sample_sex <- unique(
        sample_object[[]]$sex
      )

      sample_n_spots <- ncol(sample_object)

      sample_subtitle <- paste0(
        "Group: ",
        sample_group[[1]],
        " | Sex: ",
        sample_sex[[1]],
        "\nSpots: ",
        format(sample_n_spots, big.mark = " ", scientific = FALSE)
      )

      spatial_plot <- Seurat::SpatialDimPlot(
        object = sample_object,
        group.by = display_column,
        images = current_image_name,
        cols = cluster_palette,
        crop = crop,
        image.alpha = if (show_image) image_alpha else 0,
        pt.size.factor = pt.size.factor,
        combine = TRUE,
        stroke = NA
      ) +
        ggplot2::labs(
          title = sample_id,
          subtitle = sample_subtitle,
          color = paste0(
            cluster_column,
            " (global n, %, mean\u00b1SD)"
          ),
          fill = paste0(
            cluster_column,
            " (global n, %, mean\u00b1SD)"
          )
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 12,
            face = "bold",
            hjust = 0.5
          ),
          plot.subtitle = ggplot2::element_text(
            size = 8,
            hjust = 0.5,
            lineheight = 1.0
          ),
          legend.position = legend_position,
          legend.title = ggplot2::element_text(
            size = 9,
            face = "bold"
          ),
          legend.text = ggplot2::element_text(
            size = 8
          ),
          panel.border = ggplot2::element_rect(
            colour = "black",
            fill = NA,
            linewidth = 0.6
          ),
          plot.margin = ggplot2::margin(
            t = 6,
            r = 6,
            b = 6,
            l = 6
          )
        ) +
        ggplot2::guides(
          colour = ggplot2::guide_legend(
            ncol = legend_ncol,
            override.aes = list(size = legend_point_size, alpha = 1),
            title.position = "top"
          ),
          fill = ggplot2::guide_legend(
            ncol = legend_ncol,
            override.aes = list(size = legend_point_size, alpha = 1),
            title.position = "top"
          )
        )

      spatial_plot
    }
  )

  names(sample_plots) <- sample_order

  # ============================================================================
  # Build the shared legend as a separate row
  #
  # This guarantees the following vertical order:
  #   main title
  #   main subtitle
  #   shared legend
  #   sample titles and subtitles
  #   spatial panels
  #
  # The legend therefore cannot overlap the first row of sample titles.
  # ============================================================================

  legend_source_plot <- sample_plots[[1]] +
    ggplot2::theme(
      legend.position = "top",
      legend.box = "vertical",
      legend.justification = "center",
      legend.title = ggplot2::element_text(
        size = 10,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 9
      ),
      legend.margin = ggplot2::margin(
        t = 2,
        r = 4,
        b = 8,
        l = 4
      )
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = legend_ncol,
        override.aes = list(
          size = legend_point_size,
          alpha = 1
        ),
        title.position = "top"
      ),
      fill = ggplot2::guide_legend(
        ncol = legend_ncol,
        override.aes = list(
          size = legend_point_size,
          alpha = 1
        ),
        title.position = "top"
      )
    )

  shared_legend_grob <- extract_shared_plot_legend(
    legend_source_plot
  )

  legend_panel <- patchwork::wrap_elements(
    full = shared_legend_grob
  )

  sample_plots_without_legends <- lapply(
    sample_plots,
    function(sample_plot) {
      sample_plot +
        ggplot2::theme(
          legend.position = "none"
        )
    }
  )

  sample_grid <- patchwork::wrap_plots(
    sample_plots_without_legends,
    ncol = ncol
  )

  combined_plot <- (
    legend_panel /
      sample_grid
  ) +
    patchwork::plot_layout(
      heights = c(
        legend_height_ratio,
        1
      )
    ) +
    patchwork::plot_annotation(
      title = plot_title,
      subtitle = plot_subtitle,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 18,
          face = "bold",
          hjust = 0.5
        ),
        plot.subtitle = ggplot2::element_text(
          size = 10,
          hjust = 0.5,
          lineheight = 1.08,
          margin = ggplot2::margin(
            t = 2,
            b = 6
          )
        )
      )
    )

  output_files <- list()

  if (!is.null(output_dir)) {

    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    if (is.null(output_prefix)) {
      output_prefix <- cluster_column
    }

    summary_file <- file.path(
      output_dir,
      paste0(
        output_prefix,
        "_clusterSummary.tsv"
      )
    )

    image_map_file <- file.path(
      output_dir,
      paste0(
        output_prefix,
        "_imageMap.tsv"
      )
    )

    counts_by_sample_file <- file.path(
      output_dir,
      paste0(
        output_prefix,
        "_clusterCountsBySample.tsv"
      )
    )

    utils::write.table(
      x = cluster_summary,
      file = summary_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )

    utils::write.table(
      x = image_map,
      file = image_map_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )

    utils::write.table(
      x = cluster_counts_by_sample,
      file = counts_by_sample_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )

    output_files$clusterSummary <- summary_file
    output_files$imageMap <- image_map_file
    output_files$clusterCountsBySample <- counts_by_sample_file

    if (isTRUE(save_png)) {

      png_file <- file.path(
        output_dir,
        paste0(
          output_prefix,
          "_spatialClusters.png"
        )
      )

      ggplot2::ggsave(
        filename = png_file,
        plot = combined_plot,
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

      pdf_file <- file.path(
        output_dir,
        paste0(
          output_prefix,
          "_spatialClusters.pdf"
        )
      )

      ggplot2::ggsave(
        filename = pdf_file,
        plot = combined_plot,
        width = pdf_width_in,
        height = pdf_height_in,
        units = "in",
        bg = "white",
        limitsize = FALSE
      )

      output_files$pdf <- pdf_file
    }
  }

  list(
    plot = combined_plot,
    cluster_summary = cluster_summary,
    cluster_counts_by_sample = cluster_counts_by_sample,
    image_map = image_map,
    output_files = output_files,
    display_column = display_column,
    palette_name = palette_name,
    cluster_colors_by_id = cluster_colors_by_id,
    cluster_palette = cluster_palette
  )
}


# ==============================================================================
# Save one selected clustering as a PDF
#
# This wrapper is intentionally small. It uses the plotting function above,
# but writes only the requested PDF. It does not modify the Seurat object and
# does not require a manual-annotation column yet.
# ==============================================================================

save_selected_spatial_clustering_pdf <- function(
    seurat_object,
    cluster_column,
    output_pdf,
    sample_order = NULL,
    ncol = 4L,
    plot_title = NULL,
    plot_subtitle = NULL,
    show_image = FALSE,
    image_alpha = 1,
    crop = FALSE,
    pt.size.factor = 1.8,
    legend_ncol = 4L,
    legend_point_size = 6,
    legend_height_ratio = 0.16,
    palette_name = "working30",
    custom_cluster_colors = NULL,
    pdf_width_in = 18,
    pdf_height_in = 22,
    verbose = TRUE
) {

  if (!inherits(seurat_object, "Seurat")) {
    stop(
      "`seurat_object` must be a Seurat object.",
      call. = FALSE
    )
  }

  if (
    !is.character(cluster_column) ||
      length(cluster_column) != 1L ||
      is.na(cluster_column) ||
      !nzchar(cluster_column)
  ) {
    stop(
      "`cluster_column` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (!cluster_column %in% colnames(seurat_object[[]])) {
    stop(
      "Clustering column is absent from the Seurat object: ",
      cluster_column,
      call. = FALSE
    )
  }

  if (
    !is.character(output_pdf) ||
      length(output_pdf) != 1L ||
      is.na(output_pdf) ||
      !nzchar(output_pdf)
  ) {
    stop(
      "`output_pdf` must be one non-empty file path.",
      call. = FALSE
    )
  }

  if (!grepl("\\.pdf$", output_pdf, ignore.case = TRUE)) {
    stop(
      "`output_pdf` must end with `.pdf`.",
      call. = FALSE
    )
  }

  dir.create(
    dirname(output_pdf),
    recursive = TRUE,
    showWarnings = FALSE
  )

  plot_result <- plot_spatial_clusters_all_samples(
    seurat_object = seurat_object,
    cluster_column = cluster_column,
    sample_order = sample_order,
    ncol = ncol,
    plot_title = plot_title,
    plot_subtitle = plot_subtitle,
    include_cluster_count_in_title = TRUE,
    show_image = show_image,
    image_alpha = image_alpha,
    crop = crop,
    pt.size.factor = pt.size.factor,
    legend_position = "top",
    legend_ncol = legend_ncol,
    legend_point_size = legend_point_size,
    legend_height_ratio = legend_height_ratio,
    palette_name = palette_name,
    custom_cluster_colors = custom_cluster_colors,
    output_dir = NULL,
    output_prefix = NULL,
    save_png = FALSE,
    save_pdf = FALSE,
    pdf_width_in = pdf_width_in,
    pdf_height_in = pdf_height_in,
    verbose = verbose
  )

  ggplot2::ggsave(
    filename = output_pdf,
    plot = plot_result$plot,
    width = pdf_width_in,
    height = pdf_height_in,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )

  if (
    !file.exists(output_pdf) ||
      is.na(file.info(output_pdf)$size) ||
      file.info(output_pdf)$size <= 0L
  ) {
    stop(
      "PDF was not saved correctly: ",
      output_pdf,
      call. = FALSE
    )
  }

  if (isTRUE(verbose)) {
    message(
      "Saved selected spatial clustering PDF: ",
      output_pdf
    )
  }

  invisible(
    list(
      output_pdf = output_pdf,
      cluster_column = cluster_column,
      plot = plot_result$plot,
      cluster_summary = plot_result$cluster_summary,
      cluster_counts_by_sample = plot_result$cluster_counts_by_sample,
      image_map = plot_result$image_map,
      cluster_colors_by_id = plot_result$cluster_colors_by_id,
      cluster_palette = plot_result$cluster_palette
    )
  )
}
