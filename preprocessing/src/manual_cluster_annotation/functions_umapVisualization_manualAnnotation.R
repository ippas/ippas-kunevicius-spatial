# ==============================================================================
# functions_umapVisualization_manualAnnotation.R
#
# Purpose:
# Sandbox functions used during manual cluster annotation.
# Calculate or reuse one UMAP representation and save a selected clustering as
# a publication-quality UMAP plot.
#
# Main features:
# - does not modify the production `functions_umapVisualization.R` file;
# - calculates UMAP only when the requested reduction is absent;
# - maps colours directly to the original cluster IDs;
# - accepts optional human-readable cluster names without changing metadata;
# - keeps the legend statistics used by the spatial-clustering plot:
#     cluster_id — cluster_name (n=...; ...%; mean±SD=...±...)
# - places the shared legend in a dedicated row above the UMAP;
# - saves one selected UMAP clustering as a PDF.
#
# Dependency:
# Source `functions_spatialClusterVisualization_manualAnnotation.R` before this
# file. The UMAP functions reuse validated summary, palette and legend helpers
# defined there.
# ==============================================================================


validate_manual_annotation_umap_helper_functions <- function() {

  required_functions <- c(
    "build_cluster_summary_table",
    "format_cluster_percentage",
    "format_cluster_decimal",
    "generate_cluster_palette",
    "extract_shared_plot_legend"
  )

  missing_functions <- required_functions[
    !vapply(
      required_functions,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE
    )
  ]

  if (length(missing_functions) > 0L) {
    stop(
      "Required helper function(s) are unavailable: ",
      paste(missing_functions, collapse = ", "),
      ". Source `functions_spatialClusterVisualization_manualAnnotation.R` ",
      "before `functions_umapVisualization_manualAnnotation.R`.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


validate_named_cluster_character_vector <- function(
    x,
    argument_name,
    validate_colours = FALSE
) {

  if (is.null(x)) {
    return(invisible(TRUE))
  }

  if (
    !is.character(x) ||
      is.null(names(x)) ||
      length(x) == 0L ||
      anyNA(x) ||
      anyNA(names(x)) ||
      any(!nzchar(x)) ||
      any(!nzchar(names(x)))
  ) {
    stop(
      "`",
      argument_name,
      "` must be a non-empty named character vector, for example ",
      "`c(\"1\" = \"Cortex\")`.",
      call. = FALSE
    )
  }

  if (anyDuplicated(names(x)) > 0L) {
    duplicated_ids <- unique(
      names(x)[duplicated(names(x))]
    )

    stop(
      "Duplicated cluster IDs were found in `",
      argument_name,
      "`: ",
      paste(duplicated_ids, collapse = ", "),
      call. = FALSE
    )
  }

  if (isTRUE(validate_colours)) {
    invalid_values <- x[
      !vapply(
        x,
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

    if (length(invalid_values) > 0L) {
      stop(
        "Invalid colour value(s) in `",
        argument_name,
        "`: ",
        paste(unique(invalid_values), collapse = ", "),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


build_manual_annotation_umap_display_column <- function(
    seurat_object,
    cluster_column,
    sample_order = NULL,
    custom_cluster_labels = NULL,
    include_cluster_id_in_label = TRUE,
    display_column = NULL,
    verbose = TRUE
) {

  validate_manual_annotation_umap_helper_functions()

  if (!inherits(seurat_object, "Seurat")) {
    stop(
      "`seurat_object` must be a Seurat object.",
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

  validate_named_cluster_character_vector(
    x = custom_cluster_labels,
    argument_name = "custom_cluster_labels",
    validate_colours = FALSE
  )

  if (is.null(display_column)) {
    display_column <- paste0(
      cluster_column,
      "__manualUmapDisplay"
    )
  }

  summary_result <- build_cluster_summary_table(
    seurat_object = seurat_object,
    cluster_column = cluster_column,
    sample_order = sample_order
  )

  cluster_summary <- summary_result$cluster_summary
  cluster_ids <- as.character(cluster_summary$cluster_id)

  cluster_names <- setNames(
    object = cluster_ids,
    nm = cluster_ids
  )

  custom_label_was_supplied <- setNames(
    object = rep(FALSE, length(cluster_ids)),
    nm = cluster_ids
  )

  if (!is.null(custom_cluster_labels)) {
    matching_cluster_ids <- intersect(
      cluster_ids,
      names(custom_cluster_labels)
    )

    cluster_names[matching_cluster_ids] <-
      custom_cluster_labels[matching_cluster_ids]

    custom_label_was_supplied[matching_cluster_ids] <- TRUE

    missing_custom_cluster_ids <- setdiff(
      cluster_ids,
      names(custom_cluster_labels)
    )

    if (
      length(missing_custom_cluster_ids) > 0L &&
        isTRUE(verbose)
    ) {
      warning(
        "No custom name was supplied for cluster(s): ",
        paste(missing_custom_cluster_ids, collapse = ", "),
        ". Original cluster IDs will be used for these clusters.",
        call. = FALSE
      )
    }

    unused_custom_cluster_ids <- setdiff(
      names(custom_cluster_labels),
      cluster_ids
    )

    if (
      length(unused_custom_cluster_ids) > 0L &&
        isTRUE(verbose)
    ) {
      warning(
        "Custom names were supplied for cluster ID(s) absent from `",
        cluster_column,
        "`: ",
        paste(unused_custom_cluster_ids, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  display_names <- cluster_ids

  for (cluster_id in cluster_ids) {
    if (isTRUE(custom_label_was_supplied[[cluster_id]])) {
      if (isTRUE(include_cluster_id_in_label)) {
        display_names[cluster_id == cluster_ids] <- paste0(
          cluster_id,
          " \u2014 ",
          cluster_names[[cluster_id]]
        )
      } else {
        display_names[cluster_id == cluster_ids] <-
          cluster_names[[cluster_id]]
      }
    }
  }

  cluster_summary$cluster_name <- unname(
    cluster_names[cluster_ids]
  )

  cluster_summary$display_name <- display_names

  if (anyDuplicated(cluster_summary$display_name) > 0L) {
    duplicated_display_names <- unique(
      cluster_summary$display_name[
        duplicated(cluster_summary$display_name)
      ]
    )

    stop(
      "Cluster display names must be unique. Duplicated value(s): ",
      paste(duplicated_display_names, collapse = ", "),
      ". Keep `include_cluster_id_in_label = TRUE` or provide unique names.",
      call. = FALSE
    )
  }

  cluster_summary$legend_label <- paste0(
    cluster_summary$display_name,
    " (n=",
    format(
      cluster_summary$n_spots,
      big.mark = " ",
      scientific = FALSE
    ),
    "; ",
    format_cluster_percentage(
      cluster_summary$spot_percentage,
      digits = 1L
    ),
    "%; mean\u00b1SD=",
    format_cluster_decimal(
      cluster_summary$mean_spots_per_sample,
      digits = 1L
    ),
    "\u00b1",
    format_cluster_decimal(
      cluster_summary$sd_spots_per_sample,
      digits = 1L
    ),
    ")"
  )

  display_map <- setNames(
    object = cluster_summary$legend_label,
    nm = cluster_ids
  )

  cluster_values <- as.character(
    seurat_object[[]][[cluster_column]]
  )

  display_values <- unname(
    display_map[cluster_values]
  )

  if (anyNA(display_values)) {
    stop(
      "Could not map all cluster values to UMAP legend labels for column `",
      cluster_column,
      "`.",
      call. = FALSE
    )
  }

  seurat_object[[display_column]] <- factor(
    display_values,
    levels = cluster_summary$legend_label
  )

  list(
    seurat_object = seurat_object,
    display_column = display_column,
    cluster_summary = cluster_summary,
    cluster_counts_by_sample = summary_result$cluster_counts_by_sample,
    cluster_ids = cluster_ids,
    cluster_names_by_id = cluster_names,
    display_names_by_id = setNames(
      object = cluster_summary$display_name,
      nm = cluster_ids
    ),
    legend_labels_by_id = setNames(
      object = cluster_summary$legend_label,
      nm = cluster_ids
    )
  )
}


resolve_manual_annotation_umap_cluster_colours <- function(
    cluster_ids,
    palette_name = c(
      "working30",
      "dark3"
    ),
    custom_cluster_colors = NULL,
    verbose = TRUE
) {

  validate_manual_annotation_umap_helper_functions()

  palette_name <- match.arg(
    palette_name
  )

  cluster_ids <- as.character(cluster_ids)

  if (
    length(cluster_ids) == 0L ||
      anyNA(cluster_ids) ||
      any(!nzchar(cluster_ids)) ||
      anyDuplicated(cluster_ids) > 0L
  ) {
    stop(
      "`cluster_ids` must contain unique, non-empty cluster IDs.",
      call. = FALSE
    )
  }

  validate_named_cluster_character_vector(
    x = custom_cluster_colors,
    argument_name = "custom_cluster_colors",
    validate_colours = TRUE
  )

  default_cluster_colors <- generate_cluster_palette(
    n_clusters = length(cluster_ids),
    palette_name = palette_name
  )

  cluster_colors_by_id <- setNames(
    object = default_cluster_colors,
    nm = cluster_ids
  )

  if (!is.null(custom_cluster_colors)) {
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
        "Custom colours were supplied for cluster ID(s) absent from the ",
        "selected clustering: ",
        paste(unused_custom_cluster_ids, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  cluster_colors_by_id
}


run_manual_annotation_umap_if_missing <- function(
    seurat_object,
    input_reduction = "integrated.cca",
    dims = 1:20,
    umap_reduction_name = "umap.ccaDims20",
    umap_reduction_key = "UMAPCCA_",
    n_neighbors = 30L,
    min_dist = 0.30,
    spread = 1,
    metric = "cosine",
    seed_use = 7L,
    force_umap = FALSE,
    verbose = TRUE
) {

  required_packages <- c(
    "Seurat",
    "SeuratObject"
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

  if (!inherits(seurat_object, "Seurat")) {
    stop(
      "`seurat_object` must be a Seurat object.",
      call. = FALSE
    )
  }

  available_reductions <- SeuratObject::Reductions(
    seurat_object
  )

  if (!input_reduction %in% available_reductions) {
    stop(
      "Input reduction `",
      input_reduction,
      "` is absent from the Seurat object.\nAvailable reductions: ",
      paste(available_reductions, collapse = ", "),
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

  dims <- unique(as.integer(dims))

  available_dims <- ncol(
    SeuratObject::Embeddings(
      seurat_object[[input_reduction]]
    )
  )

  if (max(dims) > available_dims) {
    stop(
      "Requested UMAP dimension ",
      max(dims),
      " exceeds the ",
      available_dims,
      " dimensions available in reduction `",
      input_reduction,
      "`.",
      call. = FALSE
    )
  }

  if (
    !is.character(umap_reduction_name) ||
      length(umap_reduction_name) != 1L ||
      is.na(umap_reduction_name) ||
      !nzchar(umap_reduction_name)
  ) {
    stop(
      "`umap_reduction_name` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (
    !is.character(umap_reduction_key) ||
      length(umap_reduction_key) != 1L ||
      is.na(umap_reduction_key) ||
      !grepl("_$", umap_reduction_key)
  ) {
    stop(
      "`umap_reduction_key` must be one character value ending in `_`.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(n_neighbors) ||
      length(n_neighbors) != 1L ||
      is.na(n_neighbors) ||
      n_neighbors < 2
  ) {
    stop(
      "`n_neighbors` must be one integer greater than or equal to 2.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(min_dist) ||
      length(min_dist) != 1L ||
      is.na(min_dist) ||
      min_dist < 0
  ) {
    stop(
      "`min_dist` must be one non-negative numeric value.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(spread) ||
      length(spread) != 1L ||
      is.na(spread) ||
      spread <= 0
  ) {
    stop(
      "`spread` must be one positive numeric value.",
      call. = FALSE
    )
  }

  umap_exists <- umap_reduction_name %in%
    SeuratObject::Reductions(seurat_object)

  if (umap_exists && !isTRUE(force_umap)) {
    umap_embeddings <- SeuratObject::Embeddings(
      seurat_object[[umap_reduction_name]]
    )

    if (
      nrow(umap_embeddings) != ncol(seurat_object) ||
        ncol(umap_embeddings) != 2L ||
        !identical(
          rownames(umap_embeddings),
          colnames(seurat_object)
        )
    ) {
      stop(
        "Existing UMAP reduction `",
        umap_reduction_name,
        "` is not aligned with the Seurat object.",
        call. = FALSE
      )
    }

    if (isTRUE(verbose)) {
      message(
        "Reusing existing UMAP reduction: ",
        umap_reduction_name
      )
    }

    return(
      list(
        seurat_object = seurat_object,
        umap_reduction_name = umap_reduction_name,
        umap_was_calculated = FALSE
      )
    )
  }

  if (isTRUE(verbose)) {
    message("")
    message("============================================================")
    message("Calculating UMAP for manual cluster annotation")
    message("============================================================")
    message("Input reduction: ", input_reduction)
    message(
      "Dimensions: ",
      min(dims),
      "-",
      max(dims)
    )
    message("Output reduction: ", umap_reduction_name)
    message("n.neighbors: ", n_neighbors)
    message("min.dist: ", min_dist)
    message("spread: ", spread)
    message("metric: ", metric)
    message("seed.use: ", seed_use)
  }

  umap_started <- Sys.time()

  updated_seurat_object <- Seurat::RunUMAP(
    object = seurat_object,
    reduction = input_reduction,
    dims = dims,
    reduction.name = umap_reduction_name,
    reduction.key = umap_reduction_key,
    umap.method = "uwot",
    n.neighbors = as.integer(n_neighbors),
    min.dist = as.numeric(min_dist),
    spread = as.numeric(spread),
    metric = metric,
    seed.use = as.integer(seed_use),
    return.model = FALSE,
    verbose = verbose
  )

  umap_finished <- Sys.time()

  umap_elapsed_seconds <- as.numeric(
    difftime(
      umap_finished,
      umap_started,
      units = "secs"
    )
  )

  umap_embeddings <- SeuratObject::Embeddings(
    updated_seurat_object[[umap_reduction_name]]
  )

  if (
    nrow(umap_embeddings) != ncol(updated_seurat_object) ||
      ncol(umap_embeddings) != 2L ||
      !identical(
        rownames(umap_embeddings),
        colnames(updated_seurat_object)
      )
  ) {
    stop(
      "Calculated UMAP reduction has invalid dimensions or spot ordering.",
      call. = FALSE
    )
  }

  if (is.null(updated_seurat_object@misc$manualAnnotationUmap)) {
    updated_seurat_object@misc$manualAnnotationUmap <- list()
  }

  updated_seurat_object@misc$manualAnnotationUmap[[umap_reduction_name]] <-
    list(
      inputReduction = input_reduction,
      dims = dims,
      nNeighbors = as.integer(n_neighbors),
      minDist = as.numeric(min_dist),
      spread = as.numeric(spread),
      metric = metric,
      seedUse = as.integer(seed_use),
      umapMethod = "uwot",
      elapsedSeconds = umap_elapsed_seconds,
      completedAt = format(
        umap_finished,
        "%Y-%m-%d %H:%M:%S %z"
      )
    )

  if (isTRUE(verbose)) {
    message(
      "UMAP completed in ",
      round(umap_elapsed_seconds, 2),
      " seconds."
    )
  }

  list(
    seurat_object = updated_seurat_object,
    umap_reduction_name = umap_reduction_name,
    umap_was_calculated = TRUE
  )
}


plot_manual_annotation_umap_clusters <- function(
    seurat_object,
    cluster_column,
    umap_reduction_name,
    sample_order = NULL,
    plot_title = NULL,
    plot_subtitle = NULL,
    include_cluster_count_in_title = TRUE,
    palette_name = c(
      "working30",
      "dark3"
    ),
    custom_cluster_colors = NULL,
    custom_cluster_labels = NULL,
    include_cluster_id_in_label = TRUE,
    pt_size = 0.45,
    point_alpha = 0.90,
    shuffle = TRUE,
    shuffle_seed = 7L,
    raster = TRUE,
    raster_dpi = c(
      512,
      512
    ),
    legend_ncol = 3L,
    legend_point_size = 6,
    legend_height_ratio = 0.32,
    output_dir = NULL,
    output_prefix = NULL,
    save_png = TRUE,
    save_pdf = TRUE,
    png_width_in = 12,
    png_height_in = 12,
    pdf_width_in = 12,
    pdf_height_in = 12,
    dpi = 300,
    verbose = TRUE
) {

  validate_manual_annotation_umap_helper_functions()

  palette_name <- match.arg(
    palette_name
  )

  required_packages <- c(
    "Seurat",
    "SeuratObject",
    "ggplot2",
    "patchwork"
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

  if (!inherits(seurat_object, "Seurat")) {
    stop(
      "`seurat_object` must be a Seurat object.",
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
    !umap_reduction_name %in%
      SeuratObject::Reductions(seurat_object)
  ) {
    stop(
      "UMAP reduction is absent from the Seurat object: ",
      umap_reduction_name,
      call. = FALSE
    )
  }

  umap_embeddings <- SeuratObject::Embeddings(
    seurat_object[[umap_reduction_name]]
  )

  if (
    nrow(umap_embeddings) != ncol(seurat_object) ||
      ncol(umap_embeddings) != 2L ||
      !identical(
        rownames(umap_embeddings),
        colnames(seurat_object)
      )
  ) {
    stop(
      "UMAP reduction is not aligned with the Seurat object.",
      call. = FALSE
    )
  }

  if (!"sample_ID" %in% colnames(seurat_object[[]])) {
    stop(
      "The Seurat object does not contain `sample_ID` metadata.",
      call. = FALSE
    )
  }

  if (is.null(sample_order)) {
    sample_order <- unique(
      as.character(
        seurat_object[[]]$sample_ID
      )
    )
  }

  missing_samples <- setdiff(
    sample_order,
    unique(as.character(seurat_object[[]]$sample_ID))
  )

  if (length(missing_samples) > 0L) {
    stop(
      "The following requested sample(s) are absent from the Seurat object: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  display_result <- build_manual_annotation_umap_display_column(
    seurat_object = seurat_object,
    cluster_column = cluster_column,
    sample_order = sample_order,
    custom_cluster_labels = custom_cluster_labels,
    include_cluster_id_in_label = include_cluster_id_in_label,
    verbose = verbose
  )

  plotting_object <- display_result$seurat_object
  display_column <- display_result$display_column
  cluster_summary <- display_result$cluster_summary
  cluster_ids <- display_result$cluster_ids
  cluster_levels <- cluster_summary$legend_label

  cluster_colors_by_id <-
    resolve_manual_annotation_umap_cluster_colours(
      cluster_ids = cluster_ids,
      palette_name = palette_name,
      custom_cluster_colors = custom_cluster_colors,
      verbose = verbose
    )

  cluster_palette <- setNames(
    object = unname(cluster_colors_by_id[cluster_ids]),
    nm = cluster_levels
  )

  if (is.null(plot_title)) {
    plot_title <- paste0(
      "UMAP clustering: ",
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
      "Reduction: ",
      umap_reduction_name,
      " | spots: ",
      format(
        ncol(plotting_object),
        big.mark = " ",
        scientific = FALSE
      ),
      " | legend: global n, %, mean\u00b1SD across samples"
    )
  }

  if (isTRUE(verbose)) {
    message(
      "Creating manual-annotation UMAP plot for `",
      cluster_column,
      "` using reduction `",
      umap_reduction_name,
      "`."
    )
  }

  set.seed(
    as.integer(shuffle_seed)
  )

  umap_plot_with_legend <- Seurat::DimPlot(
    object = plotting_object,
    reduction = umap_reduction_name,
    group.by = display_column,
    cols = cluster_palette,
    pt.size = pt_size,
    alpha = point_alpha,
    order = NULL,
    shuffle = shuffle,
    seed = as.integer(shuffle_seed),
    raster = raster,
    raster.dpi = raster_dpi,
    combine = TRUE
  ) +
    ggplot2::labs(
      x = "UMAP_1",
      y = "UMAP_2",
      colour = paste0(
        cluster_column,
        " (global n, %, mean\u00b1SD)"
      )
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(
      base_size = 12
    ) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.7
      ),
      axis.title = ggplot2::element_text(
        size = 12,
        face = "bold"
      ),
      axis.text = ggplot2::element_text(
        size = 9
      ),
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
      ),
      plot.margin = ggplot2::margin(
        t = 8,
        r = 8,
        b = 8,
        l = 8
      )
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = as.integer(legend_ncol),
        override.aes = list(
          size = legend_point_size,
          alpha = 1
        ),
        title.position = "top"
      )
    )

  shared_legend_grob <- extract_shared_plot_legend(
    umap_plot_with_legend
  )

  legend_panel <- patchwork::wrap_elements(
    full = shared_legend_grob
  )

  umap_plot_without_legend <-
    umap_plot_with_legend +
    ggplot2::theme(
      legend.position = "none"
    )

  combined_plot <- (
    legend_panel /
      umap_plot_without_legend
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

    if (isTRUE(save_png)) {
      png_file <- file.path(
        output_dir,
        paste0(
          output_prefix,
          "_umapClusters.png"
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
          "_umapClusters.pdf"
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
    cluster_counts_by_sample = display_result$cluster_counts_by_sample,
    output_files = output_files,
    display_column = display_column,
    palette_name = palette_name,
    cluster_colors_by_id = cluster_colors_by_id,
    cluster_palette = cluster_palette,
    cluster_names_by_id = display_result$cluster_names_by_id,
    display_names_by_id = display_result$display_names_by_id,
    legend_labels_by_id = display_result$legend_labels_by_id,
    umap_reduction_name = umap_reduction_name
  )
}


# ==============================================================================
# Save one selected clustering as a UMAP PDF
#
# The wrapper calculates the requested UMAP only when it is missing, then saves
# one PDF. It returns the updated Seurat object because a newly calculated UMAP
# exists only in memory unless the caller explicitly saves the object.
# ==============================================================================

save_selected_umap_clustering_pdf <- function(
    seurat_object,
    cluster_column,
    output_pdf,
    input_reduction = "integrated.cca",
    dims = 1:20,
    umap_reduction_name = "umap.ccaDims20",
    umap_reduction_key = "UMAPCCA_",
    n_neighbors = 30L,
    min_dist = 0.30,
    spread = 1,
    metric = "cosine",
    seed_use = 7L,
    force_umap = FALSE,
    sample_order = NULL,
    plot_title = NULL,
    plot_subtitle = NULL,
    palette_name = "working30",
    custom_cluster_colors = NULL,
    custom_cluster_labels = NULL,
    include_cluster_id_in_label = TRUE,
    pt_size = 0.45,
    point_alpha = 0.90,
    shuffle = TRUE,
    shuffle_seed = 7L,
    raster = TRUE,
    raster_dpi = c(
      512,
      512
    ),
    legend_ncol = 3L,
    legend_point_size = 6,
    legend_height_ratio = 0.32,
    pdf_width_in = 12,
    pdf_height_in = 12,
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

  umap_result <- run_manual_annotation_umap_if_missing(
    seurat_object = seurat_object,
    input_reduction = input_reduction,
    dims = dims,
    umap_reduction_name = umap_reduction_name,
    umap_reduction_key = umap_reduction_key,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    spread = spread,
    metric = metric,
    seed_use = seed_use,
    force_umap = force_umap,
    verbose = verbose
  )

  plotting_result <- plot_manual_annotation_umap_clusters(
    seurat_object = umap_result$seurat_object,
    cluster_column = cluster_column,
    umap_reduction_name = umap_result$umap_reduction_name,
    sample_order = sample_order,
    plot_title = plot_title,
    plot_subtitle = plot_subtitle,
    include_cluster_count_in_title = TRUE,
    palette_name = palette_name,
    custom_cluster_colors = custom_cluster_colors,
    custom_cluster_labels = custom_cluster_labels,
    include_cluster_id_in_label = include_cluster_id_in_label,
    pt_size = pt_size,
    point_alpha = point_alpha,
    shuffle = shuffle,
    shuffle_seed = shuffle_seed,
    raster = raster,
    raster_dpi = raster_dpi,
    legend_ncol = legend_ncol,
    legend_point_size = legend_point_size,
    legend_height_ratio = legend_height_ratio,
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
    plot = plotting_result$plot,
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
      "Saved selected UMAP clustering PDF: ",
      output_pdf
    )
  }

  invisible(
    list(
      seurat_object = umap_result$seurat_object,
      output_pdf = output_pdf,
      cluster_column = cluster_column,
      umap_reduction_name = umap_result$umap_reduction_name,
      umap_was_calculated = umap_result$umap_was_calculated,
      plot = plotting_result$plot,
      cluster_summary = plotting_result$cluster_summary,
      cluster_counts_by_sample = plotting_result$cluster_counts_by_sample,
      cluster_colors_by_id = plotting_result$cluster_colors_by_id,
      cluster_palette = plotting_result$cluster_palette,
      cluster_names_by_id = plotting_result$cluster_names_by_id,
      display_names_by_id = plotting_result$display_names_by_id,
      legend_labels_by_id = plotting_result$legend_labels_by_id
    )
  )
}
