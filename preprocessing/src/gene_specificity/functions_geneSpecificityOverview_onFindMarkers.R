# ==============================================================================
# functions_geneSpecificityOverview_onFindMarkers.R
#
# Functions for creating per-cluster overview plots of gene-specificity metrics
# calculated for genes identified by Seurat FindAllMarkers().
# ==============================================================================


sort_cluster_ids_for_plot <- function(cluster_ids) {

  cluster_ids <- unique(
    as.character(cluster_ids)
  )

  numeric_ids <- suppressWarnings(
    as.numeric(cluster_ids)
  )

  if (!anyNA(numeric_ids)) {
    return(
      cluster_ids[order(numeric_ids)]
    )
  }

  sort(cluster_ids)
}


validate_named_colours <- function(
    colours,
    required_names,
    argument_name
) {

  if (
    !is.character(colours) ||
      is.null(names(colours)) ||
      anyNA(colours) ||
      anyNA(names(colours))
  ) {
    stop(
      "`",
      argument_name,
      "` must be a named character vector.",
      call. = FALSE
    )
  }

  missing_names <- setdiff(
    required_names,
    names(colours)
  )

  if (length(missing_names) > 0L) {
    stop(
      "Missing colours in `",
      argument_name,
      "` for: ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }

  invalid_colours <- colours[
    !vapply(
      colours,
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
      "Invalid colour value(s) in `",
      argument_name,
      "`: ",
      paste(unique(invalid_colours), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


get_gene_specificity_plot_theme <- function(
    theme_style = c(
      "bw",
      "classic",
      "minimal",
      "light",
      "grey"
    ),
    base_size = 9
) {

  theme_style <- match.arg(
    theme_style
  )

  switch(
    theme_style,
    bw = ggplot2::theme_bw(
      base_size = base_size
    ),
    classic = ggplot2::theme_classic(
      base_size = base_size
    ),
    minimal = ggplot2::theme_minimal(
      base_size = base_size
    ),
    light = ggplot2::theme_light(
      base_size = base_size
    ),
    grey = ggplot2::theme_grey(
      base_size = base_size
    )
  )
}


extract_shared_plot_legend <- function(plot_object) {

  plot_gtable <- ggplot2::ggplotGrob(
    plot_object
  )

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


calculate_specificity_threshold_summary <- function(
    marker_table,
    metric_columns = c(
      "tau",
      "gini",
      "shannon_specificity"
    ),
    thresholds = seq(
      0,
      1,
      by = 0.1
    )
) {

  required_columns <- c(
    "cluster",
    "gene",
    metric_columns
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(marker_table)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Missing column(s) in marker table: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  marker_table$cluster <- as.character(
    marker_table$cluster
  )

  cluster_ids <- sort_cluster_ids_for_plot(
    marker_table$cluster
  )

  summary_rows <- vector(
    mode = "list",
    length = length(cluster_ids) *
      length(metric_columns)
  )

  row_index <- 1L

  for (cluster_id in cluster_ids) {

    cluster_table <- marker_table[
      marker_table$cluster == cluster_id,
      ,
      drop = FALSE
    ]

    number_of_markers <- nrow(
      cluster_table
    )

    for (metric_column in metric_columns) {

      metric_values <- as.numeric(
        cluster_table[[metric_column]]
      )

      number_of_genes <- vapply(
        thresholds,
        function(threshold) {
          sum(
            metric_values >= threshold,
            na.rm = TRUE
          )
        },
        integer(1)
      )

      summary_rows[[row_index]] <- data.frame(
        cluster = cluster_id,
        metric = metric_column,
        threshold = thresholds,
        number_of_genes = number_of_genes,
        percentage_of_findmarkers_genes =
          100 * number_of_genes /
          number_of_markers,
        number_of_findmarkers_genes =
          number_of_markers,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      row_index <- row_index + 1L
    }
  }

  summary_table <- do.call(
    rbind,
    summary_rows
  )

  rownames(summary_table) <- NULL

  summary_table$cluster <- factor(
    summary_table$cluster,
    levels = cluster_ids
  )

  summary_table$metric <- factor(
    summary_table$metric,
    levels = metric_columns
  )

  summary_table
}


build_specificity_panel_title <- function(
    cluster_id,
    cluster_summary,
    number_of_markers,
    cluster_colour,
    cluster_names = NULL,
    cluster_symbol = "■",
    cluster_symbol_size_pt = 76
) {

  cluster_id <- as.character(
    cluster_id
  )

  summary_row <- cluster_summary[
    as.character(cluster_summary$cluster_id) ==
      cluster_id,
    ,
    drop = FALSE
  ]

  if (nrow(summary_row) != 1L) {
    stop(
      "Could not identify one cluster-summary row for cluster ",
      cluster_id,
      ".",
      call. = FALSE
    )
  }

  displayed_cluster_name <- ""

  if (
    !is.null(cluster_names) &&
      cluster_id %in% names(cluster_names) &&
      !is.na(cluster_names[[cluster_id]]) &&
      nzchar(cluster_names[[cluster_id]])
  ) {
    displayed_cluster_name <- paste0(
      " — ",
      cluster_names[[cluster_id]]
    )
  }

  paste0(
    "<span style='color:",
    cluster_colour,
    ";font-size:",
    cluster_symbol_size_pt,
    "pt;font-weight:bold;vertical-align:middle;line-height:0.8'>",
    cluster_symbol,
    "</span> ",
    "Cluster ",
    cluster_id,
    displayed_cluster_name,
    "<br>",
    "<span style='font-size:8pt;font-weight:normal'>",
    "spots: n=",
    format(
      summary_row$n_spots,
      big.mark = " ",
      scientific = FALSE
    ),
    "; ",
    formatC(
      summary_row$spot_percentage,
      format = "f",
      digits = 1
    ),
    "%; mean±SD=",
    formatC(
      summary_row$mean_spots_per_sample,
      format = "f",
      digits = 1
    ),
    "±",
    formatC(
      summary_row$sd_spots_per_sample,
      format = "f",
      digits = 1
    ),
    "; markers=",
    format(
      number_of_markers,
      big.mark = " ",
      scientific = FALSE
    ),
    "</span>"
  )
}


plot_gene_specificity_threshold_grid <- function(
    marker_table,
    cluster_summary,
    custom_cluster_colors,
    cluster_names = NULL,
    metric_columns = c(
      "tau",
      "gini",
      "shannon_specificity"
    ),
    metric_labels = c(
      "tau" = "Tau",
      "gini" = "Gini",
      "shannon_specificity" =
        "Shannon specificity"
    ),
    metric_colors = c(
      "tau" = "#1B9E77",
      "gini" = "#D95F02",
      "shannon_specificity" =
        "#7570B3"
    ),
    metric_linetypes = c(
      "tau" = "solid",
      "gini" = "dashed",
      "shannon_specificity" =
        "dotdash"
    ),
    metric_shapes = c(
      "tau" = 16,
      "gini" = 17,
      "shannon_specificity" = 15
    ),
    thresholds = seq(
      0,
      1,
      by = 0.1
    ),
    y_mode = c(
      "count",
      "percentage"
    ),
    ncol = 4L,
    plot_title = NULL,
    plot_subtitle = NULL,
    line_width = 0.8,
    point_size = 1.8,
    theme_style = c(
      "bw",
      "classic",
      "minimal",
      "light",
      "grey"
    ),
    legend_position = "top",
    legend_title_size = 14,
    legend_text_size = 13,
    legend_point_size = 5.5,
    legend_line_width = 2.2,
    legend_key_width_cm = 2.0,
    legend_key_height_cm = 0.8,
    legend_height_ratio = 0.05,
    cluster_symbol = "■",
    cluster_symbol_size_pt = 76
) {

  y_mode <- match.arg(
    y_mode
  )

  theme_style <- match.arg(
    theme_style
  )

  cluster_ids <- sort_cluster_ids_for_plot(
    marker_table$cluster
  )

  validate_named_colours(
    colours = custom_cluster_colors,
    required_names = cluster_ids,
    argument_name = "custom_cluster_colors"
  )

  validate_named_colours(
    colours = metric_colors,
    required_names = metric_columns,
    argument_name = "metric_colors"
  )

  if (!all(
    metric_columns %in% names(metric_labels)
  )) {
    stop(
      "`metric_labels` must contain names for all metric columns.",
      call. = FALSE
    )
  }

  if (!all(
    metric_columns %in% names(metric_linetypes)
  )) {
    stop(
      "`metric_linetypes` must contain names for all metric columns.",
      call. = FALSE
    )
  }

  if (!all(
    metric_columns %in% names(metric_shapes)
  )) {
    stop(
      "`metric_shapes` must contain names for all metric columns.",
      call. = FALSE
    )
  }

  plot_data <- calculate_specificity_threshold_summary(
    marker_table = marker_table,
    metric_columns = metric_columns,
    thresholds = thresholds
  )

  if (y_mode == "count") {

    y_column <- "number_of_genes"
    y_label <- "Number of genes"

    y_limits <- c(
      0,
      max(
        plot_data$number_of_findmarkers_genes
      )
    )

  } else {

    y_column <-
      "percentage_of_findmarkers_genes"

    y_label <-
      "FindAllMarkers genes (%)"

    y_limits <- c(
      0,
      100
    )
  }

  base_theme <- get_gene_specificity_plot_theme(
    theme_style = theme_style,
    base_size = 9
  )

  cluster_plots <- lapply(
    cluster_ids,
    function(cluster_id) {

      cluster_data <- plot_data[
        as.character(plot_data$cluster) ==
          cluster_id,
        ,
        drop = FALSE
      ]

      number_of_markers <- unique(
        cluster_data$number_of_findmarkers_genes
      )

      panel_title <- build_specificity_panel_title(
        cluster_id = cluster_id,
        cluster_summary = cluster_summary,
        number_of_markers =
          number_of_markers,
        cluster_colour =
          custom_cluster_colors[[cluster_id]],
        cluster_names = cluster_names,
        cluster_symbol = cluster_symbol,
        cluster_symbol_size_pt =
          cluster_symbol_size_pt
      )

      ggplot2::ggplot(
        cluster_data,
        ggplot2::aes(
          x = threshold,
          y = .data[[y_column]],
          colour = metric,
          linetype = metric,
          shape = metric,
          group = metric
        )
      ) +
        ggplot2::geom_line(
          linewidth = line_width
        ) +
        ggplot2::geom_point(
          size = point_size
        ) +
        ggplot2::scale_colour_manual(
          values =
            metric_colors[metric_columns],
          breaks = metric_columns,
          labels =
            metric_labels[metric_columns],
          drop = FALSE
        ) +
        ggplot2::scale_linetype_manual(
          values =
            metric_linetypes[metric_columns],
          breaks = metric_columns,
          labels =
            metric_labels[metric_columns],
          drop = FALSE
        ) +
        ggplot2::scale_shape_manual(
          values =
            metric_shapes[metric_columns],
          breaks = metric_columns,
          labels =
            metric_labels[metric_columns],
          drop = FALSE
        ) +
        ggplot2::scale_x_continuous(
          limits = c(
            0,
            1
          ),
          breaks = seq(
            0,
            1,
            by = 0.1
          ),
          labels = sprintf(
            "%.1f",
            seq(
              0,
              1,
              by = 0.1
            )
          ),
          expand = ggplot2::expansion(
            mult = c(
              0,
              0
            )
          )
        ) +
        ggplot2::scale_y_continuous(
          limits = y_limits,
          expand = ggplot2::expansion(
            mult = c(
              0,
              0.05
            )
          )
        ) +
        ggplot2::labs(
          title = panel_title,
          x = "Specificity threshold",
          y = y_label,
          colour = "Metric",
          linetype = "Metric",
          shape = "Metric"
        ) +
        base_theme +
        ggplot2::theme(
          plot.title = ggtext::element_markdown(
            size = 9.5,
            face = "bold",
            hjust = 0,
            lineheight = 0.95,
            margin = ggplot2::margin(
              b = 5
            )
          ),
          axis.title = ggplot2::element_text(
            size = 8.5
          ),
          axis.text = ggplot2::element_text(
            size = 7.5
          ),
          axis.text.x = ggplot2::element_text(
            angle = 45,
            hjust = 1
          ),
          panel.grid.minor =
            ggplot2::element_blank(),
          panel.border =
            ggplot2::element_rect(
              colour = "black",
              fill = NA,
              linewidth = 0.5
            ),
          legend.position = "none",
          plot.margin = ggplot2::margin(
            t = 6,
            r = 6,
            b = 6,
            l = 6
          )
        )
    }
  )

  legend_source_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = threshold,
      y = .data[[y_column]],
      colour = metric,
      linetype = metric,
      shape = metric,
      group = metric
    )
  ) +
    ggplot2::geom_line(
      linewidth = line_width
    ) +
    ggplot2::geom_point(
      size = point_size
    ) +
    ggplot2::scale_colour_manual(
      values =
        metric_colors[metric_columns],
      breaks = metric_columns,
      labels =
        metric_labels[metric_columns],
      drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      values =
        metric_linetypes[metric_columns],
      breaks = metric_columns,
      labels =
        metric_labels[metric_columns],
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values =
        metric_shapes[metric_columns],
      breaks = metric_columns,
      labels =
        metric_labels[metric_columns],
      drop = FALSE
    ) +
    ggplot2::labs(
      colour = "Metric",
      linetype = "Metric",
      shape = "Metric"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = legend_position,
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.title = ggplot2::element_text(
        size = legend_title_size,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = legend_text_size
      ),
      legend.key.width = grid::unit(
        legend_key_width_cm,
        "cm"
      ),
      legend.key.height = grid::unit(
        legend_key_height_cm,
        "cm"
      ),
      legend.margin = ggplot2::margin(
        t = 1,
        r = 4,
        b = 4,
        l = 4
      )
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = legend_line_width,
          size = legend_point_size
        )
      ),
      linetype = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = legend_line_width,
          size = legend_point_size
        )
      ),
      shape = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          size = legend_point_size
        )
      )
    )

  shared_legend_grob <- extract_shared_plot_legend(
    legend_source_plot
  )

  legend_panel <- patchwork::wrap_elements(
    full = shared_legend_grob
  )

  sample_grid <- patchwork::wrap_plots(
    cluster_plots,
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
          size = 17,
          face = "bold",
          hjust = 0.5
        ),
        plot.subtitle =
          ggplot2::element_text(
            size = 10,
            hjust = 0.5,
            lineheight = 1.05,
            margin = ggplot2::margin(
              b = 5
            )
          )
      )
    )

  list(
    plot = combined_plot,
    plot_data = plot_data
  )
}
