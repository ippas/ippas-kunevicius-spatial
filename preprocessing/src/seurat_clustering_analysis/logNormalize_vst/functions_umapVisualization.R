# ==============================================================================
# functions_umapVisualization.R
#
# Purpose:
# Calculate one UMAP representation from an existing integrated reduction and
# generate UMAP cluster plots for multiple clustering metadata columns.
#
# Important:
# - UMAP is calculated once and reused for every clustering column;
# - cluster colours are identical to the spatial-cluster plots;
# - legend labels contain:
#     cluster_id (n=...; ...%; mean±SD=...±...)
# - the legend is placed in a dedicated row above the UMAP;
# - the legend uses fewer columns so long labels do not artificially widen
#   the main UMAP panel;
# - PNG and PDF outputs are supported.
#
# Dependency:
# Source `functions_spatialClusterVisualization.R` before this file.
# The UMAP functions reuse the validated cluster-summary and palette helpers
# defined there.
# ==============================================================================


validate_umap_helper_functions <- function() {

  required_functions <- c(
    "add_display_cluster_column",
    "generate_cluster_palette",
    "extract_shared_plot_legend",
    "parse_cluster_column_name"
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
      ". Source `functions_spatialClusterVisualization.R` before ",
      "`functions_umapVisualization.R`.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


build_umap_cluster_output_prefix <- function(
    analysis_prefix,
    cluster_column
) {

  validate_umap_helper_functions()

  parsed <- parse_cluster_column_name(
    cluster_column
  )

  paste0(
    analysis_prefix,
    "_",
    parsed$algorithm,
    "Res",
    parsed$resolution
  )
}


run_umap_if_missing <- function(
    seurat_object,
    input_reduction = "integrated.rpca",
    dims = 1:20,
    umap_reduction_name = "umap.rpcaDims20",
    umap_reduction_key = "UMAPRPCA_",
    n_neighbors = 30L,
    min_dist = 0.30,
    spread = 1,
    metric = "cosine",
    seed_use = 7L,
    force_umap = FALSE,
    verbose = TRUE
) {

  # ============================================================================
  # 1. Validate input object and reduction
  # ============================================================================

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

  dims <- as.integer(dims)

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
      umap_reduction_name == ""
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


  # ============================================================================
  # 2. Reuse a valid existing UMAP unless force_umap = TRUE
  # ============================================================================

  umap_exists <- umap_reduction_name %in%
    SeuratObject::Reductions(seurat_object)

  if (umap_exists && !isTRUE(force_umap)) {

    umap_embeddings <- SeuratObject::Embeddings(
      seurat_object[[umap_reduction_name]]
    )

    if (
      nrow(umap_embeddings) != ncol(seurat_object) ||
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

    if (verbose) {
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


  # ============================================================================
  # 3. Calculate one UMAP from the integrated reduction
  # ============================================================================

  if (verbose) {
    message("")
    message("============================================================")
    message("Calculating UMAP")
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

  if (
    is.null(
      updated_seurat_object@misc$umapVisualization
    )
  ) {
    updated_seurat_object@misc$umapVisualization <-
      list()
  }

  updated_seurat_object@misc$umapVisualization[[umap_reduction_name]] <-
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

  if (verbose) {
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


plot_umap_clusters <- function(
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

  validate_umap_helper_functions()

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


  # ============================================================================
  # 1. Validate UMAP and clustering inputs
  # ============================================================================

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

  if (is.null(sample_order)) {
    sample_order <- unique(
      as.character(
        seurat_object[[]]$sample_ID
      )
    )
  }


  # ============================================================================
  # 2. Create display labels and fixed categorical colours
  # ============================================================================

  display_result <- add_display_cluster_column(
    seurat_object = seurat_object,
    cluster_column = cluster_column,
    sample_order = sample_order
  )

  plotting_object <- display_result$seurat_object
  display_column <- display_result$display_column
  cluster_summary <- display_result$cluster_summary

  cluster_levels <- cluster_summary$legend_label

  cluster_palette <- setNames(
    object = generate_cluster_palette(
      n_clusters = length(cluster_levels),
      palette_name = palette_name
    ),
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
      "Spots: ",
      format(
        ncol(plotting_object),
        big.mark = " ",
        scientific = FALSE
      ),
      " | legend: global n, %, mean±SD across samples"
    )
  }

  if (verbose) {
    message(
      "Creating UMAP cluster plot for `",
      cluster_column,
      "` using reduction `",
      umap_reduction_name,
      "`."
    )
  }


  # ============================================================================
  # 3. Create UMAP and extract the shared legend
  # ============================================================================

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
        " (global n, %, mean±SD)"
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


  # ============================================================================
  # 4. Save PNG and PDF
  # ============================================================================

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
    output_files = output_files,
    display_column = display_column,
    palette_name = palette_name,
    cluster_palette = cluster_palette,
    umap_reduction_name = umap_reduction_name
  )
}
