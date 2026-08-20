# ==============================================================================
# plot_spatial_clusters_all_samples.R
#
# Plot spatial clusters across all samples from a Seurat object.
#
# Seurat v5 / SeuratObject v5 compatible.
#
# Main features:
# - clusters are shown with one fixed colour mapping across all sections
# - complete common legend, including clusters absent from a particular section
# - legend can be placed at the top and its title is above the legend symbols
# - full spatial image can be shown or hidden
# - black frame around every section can be enabled independently of histology
# - optional black border around spots; disabled by default
# - automatic title contains clustering resolution and number of clusters
# - legend can show mean +/- SD, range, and presence across sections
# - optional PNG / SVG / TSV export
# ==============================================================================


plot_spatial_clusters_all_samples <- function(
    seurat_object,
    cluster_column,

    # Metadata and naming
    project_label = "MaternalFMT",
    cluster_prefix = "clusters_res",
    sample_id_col = "sample_ID",
    fmt_donor_group_col = "fmt_donor_group",
    sex_col = "sex",
    in_tissue_col = "in_tissue",

    # Which sections and how they are arranged
    images = NULL,
    sample_order = NULL,
    ncol = 4,
    panel_size_in = 6,
    crop = FALSE,

    # Cluster palette
    palette = "glasbey",

    # Spatial spots
    pt.size.factor = 1.2,
    spot_alpha = 1,
    show_spot_border = FALSE,
    spot_border_colour = "black",
    spot_border_width = 0.30,
    spot_border_alpha = 1,

    # Histology image and frame around each panel
    show_image = TRUE,
    image_alpha = 1,
    panel_border = TRUE,
    panel_border_colour = "black",
    panel_border_width = 0.75,

    # Optional labels placed directly on spatial sections
    label = FALSE,
    label.size = 7,

    # Panel title styling
    sample_title_size = 18,

    # Figure / legend title
    outer_title = NULL,
    legend_title = NULL,
    legend_title_as_plot_title = TRUE,
    show_outer_title = FALSE,
    outer_title_size = 24,

    # Legend styling
    legend_position = "top",
    legend_ncol = 4,
    legend_stats = "mean_sd",
    legend_title_size = 22,
    legend_text_size = 14,
    legend_point_size = 10,
    legend_key_width_cm = 0.75,
    legend_key_height_cm = 0.55,

    # Saving
    save_png = FALSE,
    save_svg = FALSE,
    save_cluster_summary_tsv = FALSE,
    output_dir = NULL,
    file_prefix = NULL,
    output_width = NULL,
    output_height = NULL,
    png_dpi = 300,

    verbose = TRUE
) {

  # ============================================================================
  # 1. Validate input
  # ============================================================================

  if (!inherits(seurat_object, "Seurat")) {
    stop("`seurat_object` must be a Seurat object.")
  }

  required_packages <- c(
    "Seurat",
    "SeuratObject",
    "ggplot2",
    "patchwork",
    "cowplot"
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

  if (save_svg && !requireNamespace("svglite", quietly = TRUE)) {
    stop(
      "SVG export requires package `svglite`.\n",
      "Install it once with: install.packages('svglite')"
    )
  }

  if ((save_png || save_svg || save_cluster_summary_tsv) &&
      is.null(output_dir)) {
    stop(
      "Set `output_dir` when saving PNG, SVG, or TSV output."
    )
  }

  if (!is.numeric(ncol) || length(ncol) != 1 || ncol < 1) {
    stop("`ncol` must be one positive integer.")
  }

  if (!is.numeric(panel_size_in) ||
      length(panel_size_in) != 1 ||
      panel_size_in <= 0) {
    stop("`panel_size_in` must be one positive numeric value.")
  }

  if (!is.numeric(pt.size.factor) ||
      length(pt.size.factor) != 1 ||
      pt.size.factor <= 0) {
    stop("`pt.size.factor` must be one positive numeric value.")
  }

  if (!is.numeric(spot_alpha) ||
      length(spot_alpha) != 1 ||
      spot_alpha < 0 ||
      spot_alpha > 1) {
    stop("`spot_alpha` must be a number between 0 and 1.")
  }

  if (!is.numeric(image_alpha) ||
      length(image_alpha) != 1 ||
      image_alpha < 0 ||
      image_alpha > 1) {
    stop("`image_alpha` must be a number between 0 and 1.")
  }

  legend_position <- match.arg(
    legend_position,
    choices = c("top", "bottom", "right")
  )

  legend_stats <- match.arg(
    legend_stats,
    choices = c("none", "mean_sd", "mean_sd_range")
  )

  metadata <- seurat_object[[]]

  required_metadata_columns <- c(
    cluster_column,
    sample_id_col,
    fmt_donor_group_col,
    sex_col
  )

  missing_metadata_columns <- setdiff(
    required_metadata_columns,
    colnames(metadata)
  )

  if (length(missing_metadata_columns) > 0) {
    stop(
      "Missing metadata column(s): ",
      paste(missing_metadata_columns, collapse = ", ")
    )
  }

  # ============================================================================
  # 2. Helper functions
  # ============================================================================

  format_resolution <- function(x) {
    x <- as.character(x)
    x <- sub("0+$", "", x)
    x <- sub("\\.$", "", x)
    x
  }

  get_first_non_missing <- function(x) {
    x <- unique(as.character(x[!is.na(x)]))

    if (length(x) == 0) {
      return("NA")
    }

    x[1]
  }

  order_cluster_levels <- function(x) {
    x <- unique(as.character(x[!is.na(x)]))
    x_numeric <- suppressWarnings(as.numeric(x))

    if (length(x) > 0 && all(!is.na(x_numeric))) {
      return(x[order(x_numeric)])
    }

    sort(x)
  }

  get_cluster_palette <- function(n, palette_name) {

    # High-contrast Glasbey-like colours.
    # These are intentionally hard-coded instead of calling
    # grDevices::hcl.colors("Glasbey"), because base R does not provide
    # a palette named "Glasbey" on this system.
    glasbey_like <- c(
      "#E41A1C", "#20D737", "#377EFB", "#FF9D00", "#B20D68", "#00C8B4",
      "#8FAA00", "#00BCEB", "#C77CFF", "#A50FDE", "#FF7F7F", "#386CB0",
      "#C65A17", "#0C7C59", "#8C6D00", "#AE017E", "#E6D300", "#00A65A",
      "#F0027F", "#F768A1", "#C9A227", "#9E77B5", "#6A3D9A", "#66A61E",
      "#6BAED6", "#A65628", "#7FC97F", "#1B9E77", "#1F78B4", "#E7298A",
      "#7B3294", "#2C7FB8", "#54278F", "#A6CEE3", "#FB9A99", "#666666"
    )

    dark <- c(
      "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
      "#A6761D", "#666666", "#1F78B4", "#E31A1C", "#6A3D9A", "#B15928",
      "#33A02C", "#FF7F00", "#A6CEE3", "#FB9A99", "#CAB2D6", "#B2DF8A"
    )

    pastel <- c(
      "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
      "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
      "#A6CEE3", "#FDBF6F", "#CAB2D6", "#B2DF8A", "#FB9A99", "#E5C494"
    )

    palette_name <- tolower(palette_name)

    palette_values <- switch(
      palette_name,
      "glasbey" = glasbey_like,
      "vivid" = glasbey_like,
      "dark" = dark,
      "pastel" = pastel,
      "hue" = character(0),
      stop(
        "`palette` must be one of: glasbey, vivid, dark, pastel, hue"
      )
    )

    # If there are more clusters than the selected fixed palette supplies,
    # generate a long HCL sequence without relying on named base-R palettes.
    if (palette_name == "hue" || n > length(palette_values)) {
      return(
        grDevices::hcl(
          h = seq(15, 375, length.out = n + 1)[-1],
          c = 100,
          l = 65,
          fixup = TRUE
        )
      )
    }

    palette_values[seq_len(n)]
  }

  # ============================================================================
  # 3. Select and order spatial images
  # ============================================================================

  available_images <- SeuratObject::Images(seurat_object)

  if (length(available_images) == 0) {
    stop("No spatial images are stored in this Seurat object.")
  }

  if (is.null(images)) {
    images <- available_images
  }

  missing_images <- setdiff(images, available_images)

  if (length(missing_images) > 0) {
    stop(
      "Requested image(s) are missing:\n",
      paste(missing_images, collapse = "\n")
    )
  }

  image_sample_ids <- sub("^image_", "", images)

  if (!is.null(sample_order)) {
    sample_order <- as.character(sample_order)

    missing_from_order <- setdiff(image_sample_ids, sample_order)

    if (length(missing_from_order) > 0) {
      stop(
        "`sample_order` does not contain sample(s): ",
        paste(missing_from_order, collapse = ", ")
      )
    }

    images <- images[
      order(match(image_sample_ids, sample_order))
    ]

    image_sample_ids <- sub("^image_", "", images)
  }

  # ============================================================================
  # 4. Global cluster levels and a fixed global colour mapping
  # ============================================================================

  cluster_values <- as.character(metadata[[cluster_column]])
  cluster_levels <- order_cluster_levels(cluster_values)
  n_clusters <- length(cluster_levels)

  if (n_clusters == 0) {
    stop(
      "Cluster column `", cluster_column,
      "` contains no non-missing values."
    )
  }

  cluster_colors <- get_cluster_palette(
    n = n_clusters,
    palette_name = palette
  )

  names(cluster_colors) <- cluster_levels

  # Use a local copy so factor levels are fixed for all images, including
  # images that do not contain a particular cluster.
  plot_object <- seurat_object

  plot_object[[cluster_column]] <- factor(
    cluster_values,
    levels = cluster_levels
  )

  # ============================================================================
  # 5. Cluster-size statistics across the displayed sections
  # ============================================================================

  summary_metadata <- metadata[
    as.character(metadata[[sample_id_col]]) %in% image_sample_ids,
    ,
    drop = FALSE
  ]

  if (in_tissue_col %in% colnames(summary_metadata)) {
    summary_metadata <- summary_metadata[
      summary_metadata[[in_tissue_col]] == 1,
      ,
      drop = FALSE
    ]
  }

  # table() retains zero values for a cluster absent from a given section.
  count_matrix <- table(
    factor(
      as.character(summary_metadata[[sample_id_col]]),
      levels = image_sample_ids
    ),
    factor(
      as.character(summary_metadata[[cluster_column]]),
      levels = cluster_levels
    )
  )

  cluster_size_summary <- data.frame(
    cluster = cluster_levels,
    mean_spots_per_section = colMeans(count_matrix),
    sd_spots_per_section = apply(count_matrix, 2, stats::sd),
    min_spots_per_section = apply(count_matrix, 2, min),
    max_spots_per_section = apply(count_matrix, 2, max),
    n_sections_present = colSums(count_matrix > 0),
    n_sections_total = nrow(count_matrix),
    stringsAsFactors = FALSE
  )

  cluster_size_summary$sd_spots_per_section[
    is.na(cluster_size_summary$sd_spots_per_section)
  ] <- 0

  make_legend_label <- function(i) {
    cluster_id <- cluster_size_summary$cluster[i]

    if (legend_stats == "none") {
      return(cluster_id)
    }

    mean_value <- round(cluster_size_summary$mean_spots_per_section[i])
    sd_value <- round(cluster_size_summary$sd_spots_per_section[i])

    label_text <- paste0(
      cluster_id,
      " | ",
      mean_value,
      " +/- ",
      sd_value,
      " spots/section",
      " | present: ",
      cluster_size_summary$n_sections_present[i],
      "/",
      cluster_size_summary$n_sections_total[i]
    )

    if (legend_stats == "mean_sd") {
      return(label_text)
    }

    paste0(
      label_text,
      " | range: ",
      cluster_size_summary$min_spots_per_section[i],
      "-",
      cluster_size_summary$max_spots_per_section[i]
    )
  }

  legend_labels <- vapply(
    seq_len(n_clusters),
    make_legend_label,
    character(1)
  )

  names(legend_labels) <- cluster_levels

  # ============================================================================
  # 6. Automatic figure title and legend title
  # ============================================================================

  if (is.null(outer_title)) {
    if (startsWith(cluster_column, cluster_prefix)) {

      resolution_value <- substring(
        cluster_column,
        first = nchar(cluster_prefix) + 1
      )

      outer_title <- paste0(
        project_label,
        " — clustering resolution ",
        format_resolution(resolution_value),
        " (",
        n_clusters,
        " clusters)"
      )

    } else {

      outer_title <- paste0(
        project_label,
        " — ",
        cluster_column,
        " (",
        n_clusters,
        " clusters)"
      )
    }
  }

  if (is.null(legend_title)) {
    if (legend_title_as_plot_title) {
      legend_title <- outer_title
    } else {
      legend_title <- paste0("Clusters (n = ", n_clusters, ")")
    }
  }

  # ============================================================================
  # 7. Section titles
  # ============================================================================

  create_sample_title <- function(sample_id) {

    sample_metadata <- metadata[
      as.character(metadata[[sample_id_col]]) == sample_id,
      ,
      drop = FALSE
    ]

    fmt_group <- get_first_non_missing(
      sample_metadata[[fmt_donor_group_col]]
    )

    sex_value <- get_first_non_missing(
      sample_metadata[[sex_col]]
    )

    paste0(
      sample_id,
      "\nFMT: ", fmt_group,
      "\nSex: ", sex_value
    )
  }

  # ============================================================================
  # 8. Complete legend with title positioned above all legend symbols
  # ============================================================================

  if (is.null(legend_ncol)) {
    legend_ncol <- if (legend_position == "right") 1 else min(4, n_clusters)
  }

  legend_data <- data.frame(
    cluster = factor(cluster_levels, levels = cluster_levels),
    x = seq_along(cluster_levels),
    y = 1
  )

  legend_plot <- ggplot2::ggplot(
    legend_data,
    ggplot2::aes(x = x, y = y, fill = cluster)
  ) +
    ggplot2::geom_point(
      shape = 21,
      size = legend_point_size,
      colour = "black",
      stroke = 0.35
    ) +
    ggplot2::scale_fill_manual(
      values = cluster_colors,
      limits = cluster_levels,
      breaks = cluster_levels,
      labels = legend_labels,
      drop = FALSE,
      name = legend_title
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        ncol = legend_ncol,
        byrow = TRUE,
        title.position = "top",
        title.hjust = 0.5,
        label.position = "right",
        override.aes = list(
          shape = 21,
          size = legend_point_size,
          alpha = 1
        )
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = legend_position,
      legend.justification = "center",
      legend.title = ggplot2::element_text(
        size = legend_title_size,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 10)
      ),
      legend.text = ggplot2::element_text(
        size = legend_text_size,
        lineheight = 0.92
      ),
      legend.key.width = grid::unit(legend_key_width_cm, "cm"),
      legend.key.height = grid::unit(legend_key_height_cm, "cm"),
      legend.spacing.x = grid::unit(0.25, "cm"),
      legend.spacing.y = grid::unit(0.22, "cm")
    )

  common_legend <- cowplot::get_legend(legend_plot)

  # ============================================================================
  # 9. Plot all spatial sections
  # ============================================================================

  effective_image_alpha <- if (show_image) image_alpha else 0

  effective_stroke <- if (show_spot_border) {
    spot_border_width
  } else {
    0
  }

  effective_stroke_alpha <- if (show_spot_border) {
    spot_border_alpha
  } else {
    0
  }

  panel_border_element <- if (panel_border) {
    ggplot2::element_rect(
      colour = panel_border_colour,
      fill = NA,
      linewidth = panel_border_width
    )
  } else {
    ggplot2::element_blank()
  }

  plot_list <- lapply(seq_along(images), function(i) {

    image_name <- images[i]
    sample_id <- image_sample_ids[i]

    if (verbose) {
      message("Plotting ", cluster_column, ": ", sample_id)
    }

    current_plot <- Seurat::SpatialDimPlot(
      object = plot_object,
      group.by = cluster_column,
      images = image_name,
      cols = cluster_colors,
      crop = crop,
      pt.size.factor = pt.size.factor,
      alpha = spot_alpha,
      image.alpha = effective_image_alpha,
      shape = 21,
      stroke = effective_stroke,
      stroke.alpha = effective_stroke_alpha,
      label = label,
      label.size = label.size,
      combine = FALSE
    )

    if (is.list(current_plot) && !inherits(current_plot, "ggplot")) {
      current_plot <- current_plot[[1]]
    }

    current_plot +
      ggplot2::ggtitle(create_sample_title(sample_id)) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = sample_title_size,
          lineheight = 0.86
        ),
        legend.position = "none",
        panel.background = ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
        panel.border = panel_border_element,
        plot.background = ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
        plot.margin = ggplot2::margin(
          t = 18,
          r = 5,
          b = 5,
          l = 5
        )
      )
  })

  # ============================================================================
  # 10. Combine square section panels with the full legend
  # ============================================================================

  panel_plot <- patchwork::wrap_plots(
    plot_list,
    ncol = ncol
  )

  legend_element <- patchwork::wrap_elements(full = common_legend)

  if (legend_position == "top") {

    final_plot <- (
      legend_element / panel_plot
    ) +
      patchwork::plot_layout(heights = c(0.28, 1))

  } else if (legend_position == "bottom") {

    final_plot <- (
      panel_plot / legend_element
    ) +
      patchwork::plot_layout(heights = c(1, 0.28))

  } else {

    final_plot <- (
      panel_plot | legend_element
    ) +
      patchwork::plot_layout(widths = c(1, 0.42))
  }

  if (show_outer_title) {
    final_plot <- final_plot +
      patchwork::plot_annotation(
        title = outer_title,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = outer_title_size,
            face = "bold",
            hjust = 0.5
          )
        )
      )
  }

  # ============================================================================
  # 11. Save PNG, SVG, and/or TSV
  # ============================================================================

  output_files <- character(0)

  if (save_png || save_svg || save_cluster_summary_tsv) {

    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    if (is.null(file_prefix)) {
      file_prefix <- paste0(
        tolower(project_label),
        "_",
        cluster_column,
        "_spatial_clusters"
      )
    }

    file_prefix <- gsub(
      pattern = "[^A-Za-z0-9._-]",
      replacement = "_",
      x = file_prefix
    )

    n_rows <- ceiling(length(images) / ncol)

    # Width and height are derived so each section gets approximately a
    # panel_size_in x panel_size_in square canvas.
    if (is.null(output_width)) {
      output_width <- ncol * panel_size_in

      if (legend_position == "right") {
        output_width <- output_width + 5.5
      }
    }

    if (is.null(output_height)) {
      output_height <- n_rows * panel_size_in

      if (legend_position %in% c("top", "bottom")) {
        output_height <- output_height + max(
          2.5,
          ceiling(n_clusters / legend_ncol) * 0.72 + 1
        )
      }

      if (show_outer_title) {
        output_height <- output_height + 0.65
      }
    }

    if (save_png) {
      png_file <- file.path(
        output_dir,
        paste0(file_prefix, ".png")
      )

      ggplot2::ggsave(
        filename = png_file,
        plot = final_plot,
        width = output_width,
        height = output_height,
        units = "in",
        dpi = png_dpi,
        bg = "white",
        limitsize = FALSE
      )

      output_files <- c(output_files, png_file)

      if (verbose) {
        message("Saved PNG: ", png_file)
      }
    }

    if (save_svg) {
      svg_file <- file.path(
        output_dir,
        paste0(file_prefix, ".svg")
      )

      ggplot2::ggsave(
        filename = svg_file,
        plot = final_plot,
        width = output_width,
        height = output_height,
        units = "in",
        device = svglite::svglite,
        bg = "white",
        limitsize = FALSE
      )

      output_files <- c(output_files, svg_file)

      if (verbose) {
        message("Saved SVG: ", svg_file)
      }
    }

    if (save_cluster_summary_tsv) {
      summary_file <- file.path(
        output_dir,
        paste0(file_prefix, "_cluster_size_summary.tsv")
      )

      utils::write.table(
        cluster_size_summary,
        file = summary_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )

      output_files <- c(output_files, summary_file)

      if (verbose) {
        message("Saved cluster summary TSV: ", summary_file)
      }
    }
  }

  attr(final_plot, "cluster_size_summary") <- cluster_size_summary
  attr(final_plot, "cluster_colors") <- cluster_colors
  attr(final_plot, "output_files") <- output_files

  return(final_plot)
}


# # ==============================================================================
# # Example
# # ==============================================================================
# plot_spatial_clusters_all_samples(
#   seurat_object = maternalFMT_integrated_n20samples,
#   cluster_column = "clusters_res0.2",
#   project_label = "MaternalFMT",
#   sample_order = metadata_autismFMT$sample_ID,
#   ncol = 7,
#   panel_size_in = 6,
#   palette = "glasbey",

#   # Full image or no histology:
#   show_image = TRUE,
#   image_alpha = 1,
#   crop = FALSE,

#   # Frame around every panel. Independent of show_image:
#   panel_border = TRUE,
#   panel_border_colour = "black",
#   panel_border_width = 0.8,

#   # Spot borders are OFF by default:
#   show_spot_border = FALSE,
#   pt.size.factor = 1.2,

#   # Legend above all panels:
#   legend_position = "top",
#   legend_ncol = 5,
#   legend_stats = "mean_sd",
#   legend_title_as_plot_title = TRUE,
#   show_outer_title = FALSE,
#   legend_title_size = 22,
#   legend_text_size = 14,
#   legend_point_size = 10,

#   save_png = TRUE,
#   save_svg = FALSE,
#   save_cluster_summary_tsv = FALSE,
#   output_dir = "results/maternalFMT_n20samples/spatial_clusters_withImage",
#   png_dpi = 300,
#   verbose = TRUE
# )


