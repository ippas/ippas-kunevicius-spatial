# ==============================================================================
# functions_markerGeneVisualization.R
#
# Functions for visualizing individual marker genes on Visium sections.
#
# The functions are independent of marker-selection criteria. The calling
# script supplies:
# - target_gene,
# - one or more marker_clusters.
#
# Five visualizations are supported:
# 1. per-spot expression using one red scale;
# 2. sample-by-cluster aggregated expression using one red scale;
# 3. per-spot expression with marker clusters in red and other clusters in green;
# 4. sample-by-cluster aggregated expression with marker clusters in red and
#    other clusters in green;
# 5. cluster-level barplot with sample points and pooled percent-positive spots.
# ==============================================================================


# ==============================================================================
# 1. General helpers
# ==============================================================================

sort_marker_cluster_ids <- function(cluster_ids) {

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


format_marker_number <- function(
    x,
    digits = 3L
) {

  format(
    round(
      as.numeric(x),
      digits = as.integer(digits)
    ),
    nsmall = as.integer(digits),
    trim = TRUE,
    scientific = FALSE
  )
}


write_marker_tsv <- function(
    data,
    filename
) {

  utils::write.table(
    x = data,
    file = filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}


standardize_marker_sex <- function(x) {

  x_lower <- tolower(
    trimws(
      as.character(x)
    )
  )

  output <- rep(
    NA_character_,
    length(x_lower)
  )

  output[
    x_lower %in% c(
      "m",
      "male"
    )
  ] <- "Male"

  output[
    x_lower %in% c(
      "f",
      "female"
    )
  ] <- "Female"

  unsupported <- unique(
    x[
      is.na(output)
    ]
  )

  if (length(unsupported) > 0L) {
    stop(
      "Unsupported sex label(s): ",
      paste(unsupported, collapse = ", "),
      call. = FALSE
    )
  }

  output
}


standardize_marker_group <- function(x) {

  x_character <- trimws(
    as.character(x)
  )

  x_lower <- tolower(
    x_character
  )

  output <- rep(
    NA_character_,
    length(x_lower)
  )

  output[
    grepl(
      "asd|autism",
      x_lower
    )
  ] <- "ASD"

  output[
    grepl(
      "neurotypical|control|typical|(^|[^a-z])nt([^a-z]|$)",
      x_lower
    )
  ] <- "Neurotypical"

  unsupported <- unique(
    x_character[
      is.na(output)
    ]
  )

  if (length(unsupported) > 0L) {
    stop(
      "Unsupported donor-group label(s): ",
      paste(unsupported, collapse = ", "),
      call. = FALSE
    )
  }

  output
}


make_marker_group_key <- function(
    sex_std,
    group_std
) {

  paste0(
    tolower(sex_std),
    "_",
    tolower(group_std)
  )
}


compute_marker_square_limits <- function(
    xmin,
    xmax,
    ymin,
    ymax,
    padding_fraction = 0.03
) {

  x_range <- xmax - xmin
  y_range <- ymax - ymin

  maximum_range <- max(
    x_range,
    y_range
  )

  if (
    !is.finite(maximum_range) ||
      maximum_range <= 0
  ) {
    maximum_range <- 1
  }

  half_side <- (
    maximum_range / 2
  ) * (
    1 + padding_fraction
  )

  x_center <- (
    xmin + xmax
  ) / 2

  y_center <- (
    ymin + ymax
  ) / 2

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


extract_marker_plot_legend <- function(plot_object) {

  plot_gtable <- ggplot2::ggplotGrob(
    plot_object
  )

  guide_indices <- which(
    grepl(
      "^guide-box",
      plot_gtable$layout$name
    )
  )

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
      "Could not extract a non-empty legend.",
      call. = FALSE
    )
  }

  plot_gtable$grobs[[
    non_empty_guide_indices[[1]]
  ]]
}


map_marker_values_to_palette <- function(
    values,
    palette_colors,
    palette_values,
    limits,
    na_color = "#D9D9D9"
) {

  if (
    length(palette_colors) !=
      length(palette_values)
  ) {
    stop(
      "`palette_colors` and `palette_values` must have equal lengths.",
      call. = FALSE
    )
  }

  limits <- as.numeric(
    limits
  )

  if (
    length(limits) != 2L ||
      anyNA(limits) ||
      limits[[2]] <= limits[[1]]
  ) {
    stop(
      "`limits` must contain two increasing numeric values.",
      call. = FALSE
    )
  }

  normalized_values <- (
    values - limits[[1]]
  ) / (
    limits[[2]] - limits[[1]]
  )

  normalized_values <- pmin(
    pmax(
      normalized_values,
      0
    ),
    1
  )

  palette_function <- scales::gradient_n_pal(
    colours = palette_colors,
    values = palette_values
  )

  output_colors <- palette_function(
    normalized_values
  )

  output_colors[
    is.na(values) |
      !is.finite(values)
  ] <- na_color

  output_colors
}


create_marker_colourbar_legend <- function(
    palette_colors,
    palette_values,
    colour_max,
    legend_title,
    legend_title_size = 9,
    legend_text_size = 8,
    legend_bar_width_mm = 42,
    legend_bar_height_mm = 4.5
) {

  if (
    !is.finite(colour_max) ||
      colour_max <= 0
  ) {
    colour_max <- 1
  }

  legend_data <- data.frame(
    x = seq(
      0,
      colour_max,
      length.out = 100
    ),
    y = 1
  )

  legend_source <- ggplot2::ggplot(
    legend_data,
    ggplot2::aes(
      x = x,
      y = y,
      colour = x
    )
  ) +
    ggplot2::geom_point(
      size = 2
    ) +
    ggplot2::scale_colour_gradientn(
      colours = palette_colors,
      values = palette_values,
      limits = c(
        0,
        colour_max
      ),
      oob = scales::squish
    ) +
    ggplot2::labs(
      colour = legend_title
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.title = ggplot2::element_text(
        size = legend_title_size,
        face = "bold",
        hjust = 0.5
      ),
      legend.text = ggplot2::element_text(
        size = legend_text_size
      ),
      legend.margin = ggplot2::margin(
        t = 1,
        r = 4,
        b = 2,
        l = 4
      )
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_colourbar(
        title.position = "top",
        direction = "horizontal",
        barwidth = grid::unit(
          legend_bar_width_mm,
          "mm"
        ),
        barheight = grid::unit(
          legend_bar_height_mm,
          "mm"
        )
      )
    )

  patchwork::wrap_elements(
    full = extract_marker_plot_legend(
      legend_source
    )
  )
}


format_marker_specificity_summary <- function(
    marker_metrics_table,
    digits = 3L
) {

  if (
    is.null(marker_metrics_table) ||
      nrow(marker_metrics_table) == 0L
  ) {
    return(
      "Specificity metrics: unavailable"
    )
  }

  marker_metrics_table <- as.data.frame(
    marker_metrics_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  format_metric <- function(
      column_name,
      displayed_name
  ) {

    if (
      !column_name %in%
        colnames(marker_metrics_table)
    ) {
      return(NULL)
    }

    metric_values <- unique(
      marker_metrics_table[[column_name]]
    )

    metric_values <- metric_values[
      !is.na(metric_values)
    ]

    if (length(metric_values) == 0L) {
      return(NULL)
    }

    paste0(
      displayed_name,
      "=",
      formatC(
        as.numeric(metric_values[[1]]),
        format = "f",
        digits = digits
      )
    )
  }

  global_metrics <- Filter(
    Negate(is.null),
    list(
      format_metric(
        "tau",
        "Tau"
      ),
      format_metric(
        "gini",
        "Gini"
      ),
      format_metric(
        "shannon_specificity",
        "Shannon"
      )
    )
  )

  global_label <- if (
    length(global_metrics) > 0L
  ) {
    paste(
      global_metrics,
      collapse = " | "
    )
  } else {
    "global metrics unavailable"
  }

  cluster_labels <- vapply(
    seq_len(
      nrow(marker_metrics_table)
    ),
    function(row_index) {

      current_row <- marker_metrics_table[
        row_index,
        ,
        drop = FALSE
      ]

      cluster_id <- if (
        "cluster" %in%
          colnames(current_row)
      ) {
        as.character(
          current_row$cluster[[1]]
        )
      } else {
        "NA"
      }

      cluster_metrics <- character(0)

      if (
        "expression_specificity" %in%
          colnames(current_row) &&
          is.finite(
            current_row$expression_specificity[[1]]
          )
      ) {
        cluster_metrics <- c(
          cluster_metrics,
          paste0(
            "ExprSpec=",
            formatC(
              current_row$expression_specificity[[1]],
              format = "f",
              digits = digits
            )
          )
        )
      }

      if (
        "expression_ratio_vs_best_other" %in%
          colnames(current_row) &&
          is.finite(
            current_row$expression_ratio_vs_best_other[[1]]
          )
      ) {
        cluster_metrics <- c(
          cluster_metrics,
          paste0(
            "ratio=",
            formatC(
              current_row$expression_ratio_vs_best_other[[1]],
              format = "f",
              digits = digits
            )
          )
        )
      }

      if (
        "best_other_cluster" %in%
          colnames(current_row) &&
          !is.na(
            current_row$best_other_cluster[[1]]
          )
      ) {
        cluster_metrics <- c(
          cluster_metrics,
          paste0(
            "bestOther=C",
            as.character(
              current_row$best_other_cluster[[1]]
            )
          )
        )
      }

      if (
        "is_best_cluster" %in%
          colnames(current_row) &&
          !is.na(
            current_row$is_best_cluster[[1]]
          )
      ) {
        cluster_metrics <- c(
          cluster_metrics,
          paste0(
            "isBest=",
            as.character(
              current_row$is_best_cluster[[1]]
            )
          )
        )
      }

      if (length(cluster_metrics) == 0L) {
        paste0(
          "C",
          cluster_id
        )
      } else {
        paste0(
          "C",
          cluster_id,
          " [",
          paste(
            cluster_metrics,
            collapse = ", "
          ),
          "]"
        )
      }
    },
    character(1)
  )

  paste0(
    "Specificity: ",
    global_label,
    " | marker cluster(s): ",
    paste(
      cluster_labels,
      collapse = "; "
    )
  )
}


# ==============================================================================
# 2. Sample layout and spatial-image mapping
# ==============================================================================

build_marker_sample_layout <- function(
    seurat_object,
    sample_id_column = "sample_ID",
    group_column = "fmt_donor_group",
    sex_column = "sex",
    column_order = c(
      "male_neurotypical",
      "male_asd",
      "female_neurotypical",
      "female_asd"
    ),
    column_titles = c(
      "male_neurotypical" = "Male neurotypical",
      "male_asd" = "Male ASD",
      "female_neurotypical" = "Female neurotypical",
      "female_asd" = "Female ASD"
    )
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
      "Missing metadata column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  sample_order <- unique(
    as.character(
      metadata_table[[sample_id_column]]
    )
  )

  sample_rows <- lapply(
    sample_order,
    function(sample_id) {

      sample_metadata <- metadata_table[
        as.character(
          metadata_table[[sample_id_column]]
        ) == sample_id,
        ,
        drop = FALSE
      ]

      group_values <- unique(
        as.character(
          sample_metadata[[group_column]]
        )
      )

      sex_values <- unique(
        as.character(
          sample_metadata[[sex_column]]
        )
      )

      if (
        length(group_values) != 1L ||
          length(sex_values) != 1L
      ) {
        stop(
          "Sample ",
          sample_id,
          " does not have one unique group and sex value.",
          call. = FALSE
        )
      }

      sex_std <- standardize_marker_sex(
        sex_values
      )

      group_std <- standardize_marker_group(
        group_values
      )

      data.frame(
        sample_ID = sample_id,
        group_raw = group_values,
        sex_raw = sex_values,
        group_std = group_std,
        sex_std = sex_std,
        column_key = make_marker_group_key(
          sex_std,
          group_std
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )

  sample_table <- do.call(
    rbind,
    sample_rows
  )

  rownames(sample_table) <- NULL

  unsupported_keys <- setdiff(
    unique(sample_table$column_key),
    column_order
  )

  if (length(unsupported_keys) > 0L) {
    stop(
      "Unsupported group-column key(s): ",
      paste(unsupported_keys, collapse = ", "),
      call. = FALSE
    )
  }

  samples_by_column <- lapply(
    column_order,
    function(column_key) {
      sample_table$sample_ID[
        sample_table$column_key ==
          column_key
      ]
    }
  )

  names(samples_by_column) <-
    column_order

  list(
    sample_table = sample_table,
    sample_order = sample_order,
    samples_by_column = samples_by_column,
    max_rows = max(
      lengths(samples_by_column)
    ),
    column_order = column_order,
    column_titles = column_titles
  )
}


create_marker_image_map <- function(
    seurat_object,
    sample_id_column = "sample_ID"
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

      sample_ids <- unique(
        as.character(
          metadata_table[
            image_cells,
            sample_id_column,
            drop = TRUE
          ]
        )
      )

      if (length(sample_ids) != 1L) {
        stop(
          "Spatial image ",
          image_name,
          " maps to ",
          length(sample_ids),
          " sample IDs.",
          call. = FALSE
        )
      }

      data.frame(
        sample_ID = sample_ids[[1]],
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
    stop(
      "More than one image maps to the same sample ID.",
      call. = FALSE
    )
  }

  image_map
}


standardize_marker_coordinates <- function(
    coordinate_table
) {

  coordinate_table <- as.data.frame(
    coordinate_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (
    !"cell" %in% colnames(coordinate_table)
  ) {
    coordinate_table$cell <- rownames(
      coordinate_table
    )
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

  x_column <- x_candidates[
    x_candidates %in%
      colnames(coordinate_table)
  ][1]

  y_column <- y_candidates[
    y_candidates %in%
      colnames(coordinate_table)
  ][1]

  if (
    is.na(x_column) ||
      is.na(y_column)
  ) {
    stop(
      "Could not identify x/y columns in tissue coordinates. Available: ",
      paste(
        colnames(coordinate_table),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  data.frame(
    cell = as.character(
      coordinate_table$cell
    ),
    x_plot = as.numeric(
      coordinate_table[[x_column]]
    ),
    y_plot = as.numeric(
      coordinate_table[[y_column]]
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


get_marker_image_dimensions <- function(
    seurat_object,
    image_name,
    coordinate_table
) {

  raw_image <- tryCatch(
    SeuratObject::GetImage(
      object = seurat_object,
      image = image_name,
      mode = "raw"
    ),
    error = function(error_condition) NULL
  )

  raw_dimensions <- dim(
    raw_image
  )

  if (
    !is.null(raw_dimensions) &&
      length(raw_dimensions) >= 2L
  ) {
    return(
      list(
        image_array = raw_image,
        image_height =
          as.numeric(raw_dimensions[[1]]),
        image_width =
          as.numeric(raw_dimensions[[2]])
      )
    )
  }

  list(
    image_array = NULL,
    image_height = max(
      coordinate_table$y_plot,
      na.rm = TRUE
    ),
    image_width = max(
      coordinate_table$x_plot,
      na.rm = TRUE
    )
  )
}


# ==============================================================================
# 3. Shared context and gene-level data preparation
# ==============================================================================

prepare_marker_gene_visualization_context <- function(
    seurat_object,
    cluster_column,
    assay_name = "RNA",
    sample_id_column = "sample_ID",
    group_column = "fmt_donor_group",
    sex_column = "sex",
    image_scale = "lowres",
    panel_padding_fraction = 0.03
) {

  if (!inherits(seurat_object, "Seurat")) {
    stop(
      "`seurat_object` must be a Seurat object.",
      call. = FALSE
    )
  }

  metadata_table <- seurat_object[[]]

  if (
    !cluster_column %in%
      colnames(metadata_table)
  ) {
    stop(
      "Missing clustering column: ",
      cluster_column,
      call. = FALSE
    )
  }

  SeuratObject::DefaultAssay(
    seurat_object
  ) <- assay_name

  count_layers <- SeuratObject::Layers(
    seurat_object[[assay_name]],
    search = "^counts"
  )

  if (length(count_layers) > 1L) {
    seurat_object <- SeuratObject::JoinLayers(
      object = seurat_object,
      assay = assay_name
    )
  }

  counts_matrix <- SeuratObject::LayerData(
    object = seurat_object,
    assay = assay_name,
    layer = "counts"
  )

  total_umi_per_spot <- Matrix::colSums(
    counts_matrix
  )

  sample_layout <- build_marker_sample_layout(
    seurat_object = seurat_object,
    sample_id_column = sample_id_column,
    group_column = group_column,
    sex_column = sex_column
  )

  image_map <- create_marker_image_map(
    seurat_object = seurat_object,
    sample_id_column = sample_id_column
  )

  missing_images <- setdiff(
    sample_layout$sample_order,
    image_map$sample_ID
  )

  if (length(missing_images) > 0L) {
    stop(
      "No spatial image for sample(s): ",
      paste(missing_images, collapse = ", "),
      call. = FALSE
    )
  }

  sample_cache <- lapply(
    sample_layout$sample_order,
    function(sample_id) {

      image_name <- image_map$image_name[
        match(
          sample_id,
          image_map$sample_ID
        )
      ]

      image_cells <- SeuratObject::Cells(
        seurat_object[[image_name]]
      )

      coordinate_table <-
        SeuratObject::GetTissueCoordinates(
          object = seurat_object,
          image = image_name,
          scale = image_scale
        )

      coordinate_table <-
        standardize_marker_coordinates(
          coordinate_table
        )

      common_cells <- image_cells[
        image_cells %in%
          coordinate_table$cell &
          image_cells %in%
          colnames(counts_matrix)
      ]

      if (length(common_cells) == 0L) {
        stop(
          "No shared cells among image, coordinates and counts for sample ",
          sample_id,
          ".",
          call. = FALSE
        )
      }

      coordinate_table <- coordinate_table[
        match(
          common_cells,
          coordinate_table$cell
        ),
        ,
        drop = FALSE
      ]

      image_info <- get_marker_image_dimensions(
        seurat_object = seurat_object,
        image_name = image_name,
        coordinate_table = coordinate_table
      )

      frame_limits <- compute_marker_square_limits(
        xmin = 0,
        xmax = image_info$image_width,
        ymin = 0,
        ymax = image_info$image_height,
        padding_fraction =
          panel_padding_fraction
      )

      sample_metadata <- sample_layout$sample_table[
        sample_layout$sample_table$sample_ID ==
          sample_id,
        ,
        drop = FALSE
      ]

      spot_metadata <- metadata_table[
        common_cells,
        ,
        drop = FALSE
      ]

      coordinate_table$sample_ID <- sample_id
      coordinate_table$group_std <-
        sample_metadata$group_std[[1]]
      coordinate_table$sex_std <-
        sample_metadata$sex_std[[1]]
      coordinate_table$column_key <-
        sample_metadata$column_key[[1]]
      coordinate_table$cluster <- as.character(
        spot_metadata[[cluster_column]]
      )

      list(
        sample_ID = sample_id,
        image_name = image_name,
        cells = common_cells,
        coordinates = coordinate_table,
        image_array = image_info$image_array,
        image_width = image_info$image_width,
        image_height = image_info$image_height,
        frame_limits = frame_limits
      )
    }
  )

  names(sample_cache) <-
    sample_layout$sample_order

  list(
    seurat_object = seurat_object,
    counts_matrix = counts_matrix,
    total_umi_per_spot =
      total_umi_per_spot,
    cluster_column = cluster_column,
    assay_name = assay_name,
    sample_layout = sample_layout,
    image_map = image_map,
    sample_cache = sample_cache,
    cluster_levels =
      sort_marker_cluster_ids(
        metadata_table[[cluster_column]]
      )
  )
}


match_marker_gene_feature <- function(
    feature_names,
    target_gene
) {

  feature_names <- as.character(
    feature_names
  )

  target_gene <- as.character(
    target_gene
  )[[1]]

  canonicalize <- function(x) {
    x <- trimws(
      as.character(x)
    )
    x <- sub(
      "\\.[0-9]+$",
      "",
      x
    )
    x <- sub(
      "-[0-9]+$",
      "",
      x
    )
    x
  }

  candidate_sets <- list(
    exact = which(
      feature_names == target_gene
    ),
    case_insensitive = which(
      tolower(feature_names) ==
        tolower(target_gene)
    ),
    canonical_exact = which(
      canonicalize(feature_names) ==
        canonicalize(target_gene)
    ),
    canonical_case_insensitive = which(
      tolower(
        canonicalize(feature_names)
      ) ==
        tolower(
          canonicalize(target_gene)
        )
    )
  )

  for (match_type in names(candidate_sets)) {

    matching_indices <- unique(
      candidate_sets[[match_type]]
    )

    if (length(matching_indices) == 1L) {
      return(
        list(
          matched_index =
            matching_indices[[1]],
          matched_name =
            feature_names[
              matching_indices[[1]]
            ],
          match_type = match_type
        )
      )
    }

    if (length(matching_indices) > 1L) {
      stop(
        "Ambiguous feature match for gene ",
        target_gene,
        ".",
        call. = FALSE
      )
    }
  }

  stop(
    "Gene was not found in the counts matrix: ",
    target_gene,
    call. = FALSE
  )
}


prepare_marker_gene_data <- function(
    context,
    target_gene,
    marker_clusters,
    marker_metrics_table = NULL,
    clustering_parameters_label = NULL,
    normalization_scale_factor = 10000,
    upper_colour_quantile = 0.99
) {

  marker_clusters <- sort_marker_cluster_ids(
    marker_clusters
  )

  specificity_summary_label <-
    format_marker_specificity_summary(
      marker_metrics_table =
        marker_metrics_table
    )

  number_of_clusters <- length(
    context$cluster_levels
  )

  clustering_summary_label <- paste0(
    if (
      !is.null(
        clustering_parameters_label
      ) &&
        nzchar(
          clustering_parameters_label
        )
    ) {
      clustering_parameters_label
    } else {
      "Clustering parameters not supplied"
    },
    " | ",
    number_of_clusters,
    " clusters"
  )

  missing_clusters <- setdiff(
    marker_clusters,
    context$cluster_levels
  )

  if (length(missing_clusters) > 0L) {
    stop(
      "Marker cluster(s) absent from the Seurat object: ",
      paste(missing_clusters, collapse = ", "),
      call. = FALSE
    )
  }

  gene_match <- match_marker_gene_feature(
    feature_names = rownames(
      context$counts_matrix
    ),
    target_gene = target_gene
  )

  gene_counts <- as.numeric(
    context$counts_matrix[
      gene_match$matched_index,
      ,
      drop = TRUE
    ]
  )

  names(gene_counts) <- colnames(
    context$counts_matrix
  )

  sample_plot_data <- vector(
    mode = "list",
    length = length(
      context$sample_layout$sample_order
    )
  )

  names(sample_plot_data) <-
    context$sample_layout$sample_order

  sample_summary_rows <- list()
  sample_cluster_rows <- list()

  for (
    sample_id in
    context$sample_layout$sample_order
  ) {

    sample_cache <-
      context$sample_cache[[sample_id]]

    cells <- sample_cache$cells

    raw_counts <- as.numeric(
      gene_counts[cells]
    )

    total_umi <- as.numeric(
      context$total_umi_per_spot[cells]
    )

    log_normalized_expression <- numeric(
      length(cells)
    )

    valid_spots <- total_umi > 0

    log_normalized_expression[
      valid_spots
    ] <- log1p(
      raw_counts[valid_spots] /
        total_umi[valid_spots] *
        normalization_scale_factor
    )

    plot_data <- sample_cache$coordinates

    plot_data$target_gene <-
      target_gene
    plot_data$matched_gene <-
      gene_match$matched_name
    plot_data$target_raw_count <-
      raw_counts
    plot_data$total_UMI <-
      total_umi
    plot_data$logNormalized_expression <-
      log_normalized_expression
    plot_data$is_marker_cluster <-
      plot_data$cluster %in%
      marker_clusters

    cluster_ids_in_sample <-
      sort_marker_cluster_ids(
        plot_data$cluster
      )

    cluster_summary_rows <- lapply(
      cluster_ids_in_sample,
      function(cluster_id) {

        cluster_indices <- which(
          plot_data$cluster ==
            cluster_id
        )

        cluster_raw_count <- sum(
          plot_data$target_raw_count[
            cluster_indices
          ]
        )

        cluster_total_umi <- sum(
          plot_data$total_UMI[
            cluster_indices
          ]
        )

        cluster_expression <- if (
          cluster_total_umi > 0
        ) {
          log1p(
            cluster_raw_count /
              cluster_total_umi *
              normalization_scale_factor
          )
        } else {
          0
        }

        data.frame(
          sample_ID = sample_id,
          group_std =
            plot_data$group_std[[1]],
          sex_std =
            plot_data$sex_std[[1]],
          column_key =
            plot_data$column_key[[1]],
          cluster = cluster_id,
          is_marker_cluster =
            cluster_id %in%
            marker_clusters,
          n_spots =
            length(cluster_indices),
          n_positive_spots =
            sum(
              plot_data$target_raw_count[
                cluster_indices
              ] > 0
            ),
          percent_positive_spots =
            100 * mean(
              plot_data$target_raw_count[
                cluster_indices
              ] > 0
            ),
          total_target_raw_count =
            cluster_raw_count,
          total_UMI =
            cluster_total_umi,
          cluster_logNormalized_expression =
            cluster_expression,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    )

    cluster_summary_table <- do.call(
      rbind,
      cluster_summary_rows
    )

    plot_data <-
      merge(
        x = plot_data,
        y = cluster_summary_table[
          ,
          c(
            "cluster",
            "cluster_logNormalized_expression"
          ),
          drop = FALSE
        ],
        by = "cluster",
        all.x = TRUE,
        sort = FALSE
      )

    plot_data <- plot_data[
      match(
        cells,
        plot_data$cell
      ),
      ,
      drop = FALSE
    ]

    marker_cluster_indices <- which(
      plot_data$is_marker_cluster
    )

    sample_summary_rows[[
      length(sample_summary_rows) + 1L
    ]] <- data.frame(
      sample_ID = sample_id,
      group_std =
        plot_data$group_std[[1]],
      sex_std =
        plot_data$sex_std[[1]],
      column_key =
        plot_data$column_key[[1]],
      target_gene = target_gene,
      marker_clusters = paste(
        marker_clusters,
        collapse = ","
      ),
      n_tissue_spots =
        nrow(plot_data),
      n_positive_tissue_spots =
        sum(
          plot_data$target_raw_count >
            0
        ),
      percent_positive_tissue_spots =
        100 * mean(
          plot_data$target_raw_count >
            0
        ),
      n_marker_cluster_spots =
        length(marker_cluster_indices),
      n_positive_marker_cluster_spots =
        sum(
          plot_data$target_raw_count[
            marker_cluster_indices
          ] > 0
        ),
      percent_positive_marker_cluster_spots =
        if (
          length(marker_cluster_indices) > 0
        ) {
          100 * mean(
            plot_data$target_raw_count[
              marker_cluster_indices
            ] > 0
          )
        } else {
          NA_real_
        },
      mean_logNormalized_all_spots =
        mean(
          plot_data$logNormalized_expression
        ),
      sd_logNormalized_all_spots =
        stats::sd(
          plot_data$logNormalized_expression
        ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    sample_cluster_rows[[
      length(sample_cluster_rows) + 1L
    ]] <- cluster_summary_table

    sample_plot_data[[sample_id]] <-
      plot_data
  }

  sample_summary <- do.call(
    rbind,
    sample_summary_rows
  )

  sample_cluster_summary <- do.call(
    rbind,
    sample_cluster_rows
  )

  rownames(sample_summary) <- NULL
  rownames(sample_cluster_summary) <- NULL

  all_plot_data <- do.call(
    rbind,
    sample_plot_data
  )

  positive_spot_values <-
    all_plot_data$logNormalized_expression[
      is.finite(
        all_plot_data$logNormalized_expression
      ) &
        all_plot_data$logNormalized_expression >
        0
    ]

  positive_cluster_values <-
    sample_cluster_summary$cluster_logNormalized_expression[
      is.finite(
        sample_cluster_summary$cluster_logNormalized_expression
      ) &
        sample_cluster_summary$cluster_logNormalized_expression >
        0
    ]

  calculate_colour_max <- function(
      positive_values
  ) {

    if (length(positive_values) == 0L) {
      return(1)
    }

    colour_max <- as.numeric(
      stats::quantile(
        positive_values,
        probs = upper_colour_quantile,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
    )

    if (
      !is.finite(colour_max) ||
        colour_max <= 0
    ) {
      colour_max <- max(
        positive_values,
        na.rm = TRUE
      )
    }

    if (
      !is.finite(colour_max) ||
        colour_max <= 0
    ) {
      colour_max <- 1
    }

    colour_max
  }

  list(
    target_gene = target_gene,
    matched_gene = gene_match$matched_name,
    marker_clusters = marker_clusters,
    sample_plot_data = sample_plot_data,
    sample_summary = sample_summary,
    sample_cluster_summary =
      sample_cluster_summary,
    spot_colour_max =
      calculate_colour_max(
        positive_spot_values
      ),
    cluster_colour_max =
      calculate_colour_max(
        positive_cluster_values
      ),
    normalization_scale_factor =
      normalization_scale_factor,
    upper_colour_quantile =
      upper_colour_quantile,
    marker_metrics_table =
      marker_metrics_table,
    specificity_summary_label =
      specificity_summary_label,
    clustering_summary_label =
      clustering_summary_label,
    number_of_clusters =
      number_of_clusters
  )
}


# ==============================================================================
# 4. Panel titles and four-column layout
# ==============================================================================

make_marker_spatial_subtitle <- function(
    plot_data,
    target_gene,
    marker_clusters,
    emphasize_marker_clusters = FALSE,
    aggregated = FALSE
) {

  positive_all <- sum(
    plot_data$target_raw_count >
      0
  )

  percent_positive_all <- 100 * mean(
    plot_data$target_raw_count >
      0
  )

  subtitle_lines <- c(
    paste0(
      "Group: ",
      plot_data$group_std[[1]],
      " | Sex: ",
      plot_data$sex_std[[1]]
    ),
    paste0(
      "Tissue spots: ",
      format(
        nrow(plot_data),
        big.mark = " ",
        scientific = FALSE
      ),
      " | ",
      target_gene,
      "+: ",
      format(
        positive_all,
        big.mark = " ",
        scientific = FALSE
      ),
      " (",
      formatC(
        percent_positive_all,
        format = "f",
        digits = 1
      ),
      "%)"
    )
  )

  if (isTRUE(emphasize_marker_clusters)) {

    marker_data <- plot_data[
      plot_data$is_marker_cluster,
      ,
      drop = FALSE
    ]

    marker_positive <- sum(
      marker_data$target_raw_count >
        0
    )

    marker_percent <- if (
      nrow(marker_data) > 0
    ) {
      100 * mean(
        marker_data$target_raw_count >
          0
      )
    } else {
      NA_real_
    }

    subtitle_lines <- c(
      subtitle_lines,
      paste0(
        "Marker cluster(s) ",
        paste(
          marker_clusters,
          collapse = ", "
        ),
        ": ",
        format(
          nrow(marker_data),
          big.mark = " ",
          scientific = FALSE
        ),
        " spots | ",
        target_gene,
        "+: ",
        format(
          marker_positive,
          big.mark = " ",
          scientific = FALSE
        ),
        " (",
        if (
          is.finite(marker_percent)
        ) {
          formatC(
            marker_percent,
            format = "f",
            digits = 1
          )
        } else {
          "NA"
        },
        "%)"
      )
    )
  }

  paste(
    subtitle_lines,
    collapse = "\n"
  )
}


create_marker_spatial_panel <- function(
    sample_cache,
    plot_data,
    value_column,
    colour_column,
    target_gene,
    marker_clusters,
    show_histology_image = FALSE,
    emphasize_marker_clusters = FALSE,
    aggregated = FALSE,
    point_size_no_image = 1.1,
    point_size_with_image = 0.80,
    panel_border_linewidth = 0.8
) {

  title_text <- paste0(
    "Sample: ",
    sample_cache$sample_ID
  )

  subtitle_text <- make_marker_spatial_subtitle(
    plot_data = plot_data,
    target_gene = target_gene,
    marker_clusters = marker_clusters,
    emphasize_marker_clusters =
      emphasize_marker_clusters,
    aggregated = aggregated
  )

  if (
    isTRUE(show_histology_image) &&
      !is.null(sample_cache$image_array)
  ) {

    output_plot <- ggplot2::ggplot() +
      ggplot2::annotation_raster(
        raster = sample_cache$image_array,
        xmin = 0,
        xmax = sample_cache$image_width,
        ymin = sample_cache$image_height,
        ymax = 0
      ) +
      ggplot2::geom_point(
        data = plot_data,
        ggplot2::aes(
          x = x_plot,
          y = y_plot,
          colour = .data[[colour_column]]
        ),
        size = point_size_with_image,
        shape = 16,
        stroke = 0
      )

  } else {

    output_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = x_plot,
        y = y_plot,
        colour = .data[[colour_column]]
      )
    ) +
      ggplot2::geom_point(
        size = point_size_no_image,
        shape = 16,
        stroke = 0
      )
  }

  output_plot +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(
      limits =
        sample_cache$frame_limits$x_limits,
      expand = c(
        0,
        0
      )
    ) +
    ggplot2::scale_y_reverse(
      limits = rev(
        sample_cache$frame_limits$y_limits
      ),
      expand = c(
        0,
        0
      )
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = title_text,
      subtitle = subtitle_text
    ) +
    ggplot2::theme_void(
      base_family = "DejaVu Sans"
    ) +
    ggplot2::theme(
      plot.background =
        ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
      panel.background =
        ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
      panel.border =
        ggplot2::element_rect(
          colour = "black",
          fill = NA,
          linewidth =
            panel_border_linewidth
        ),
      plot.title =
        ggplot2::element_text(
          size = 12.5,
          face = "bold",
          hjust = 0.5,
          margin = ggplot2::margin(
            b = 4
          )
        ),
      plot.subtitle =
        ggplot2::element_text(
          size = 7.2,
          hjust = 0.5,
          lineheight = 1.05,
          margin = ggplot2::margin(
            b = 6
          )
        ),
      plot.margin =
        ggplot2::margin(
          t = 8,
          r = 8,
          b = 8,
          l = 8
        )
    )
}


arrange_marker_plots_four_columns <- function(
    plot_list_named,
    sample_layout,
    legend_panel,
    plot_title,
    plot_subtitle,
    legend_height_ratio = 0.075
) {

  plots_by_column <- lapply(
    sample_layout$column_order,
    function(column_key) {

      sample_ids <-
        sample_layout$samples_by_column[[
          column_key
        ]]

      current_plots <- lapply(
        sample_ids,
        function(sample_id) {
          plot_list_named[[sample_id]]
        }
      )

      if (
        length(current_plots) <
          sample_layout$max_rows
      ) {
        current_plots <- c(
          current_plots,
          rep(
            list(
              patchwork::plot_spacer()
            ),
            sample_layout$max_rows -
              length(current_plots)
          )
        )
      }

      current_plots
    }
  )

  names(plots_by_column) <-
    sample_layout$column_order

  interleaved_plots <- list()

  for (
    row_index in
    seq_len(sample_layout$max_rows)
  ) {
    for (
      column_key in
      sample_layout$column_order
    ) {
      interleaved_plots[[
        length(interleaved_plots) + 1L
      ]] <- plots_by_column[[
        column_key
      ]][[row_index]]
    }
  }

  spatial_grid <- patchwork::wrap_plots(
    interleaved_plots,
    ncol = 4L
  )

  (
    legend_panel /
      spatial_grid
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
        text = ggplot2::element_text(
          family = "DejaVu Sans"
        ),
        plot.title =
          ggplot2::element_text(
            size = 20,
            face = "bold",
            hjust = 0.5
          ),
        plot.subtitle =
          ggplot2::element_text(
            size = 9.5,
            hjust = 0.5,
            lineheight = 1.04,
            margin = ggplot2::margin(
              b = 4
            )
          ),
        plot.margin =
          ggplot2::margin(
            t = 6,
            r = 8,
            b = 8,
            l = 8
          )
      )
    )
}


# ==============================================================================
# 5. Five requested plotting functions
# ==============================================================================

plot_marker_gene_per_spot_all_clusters <- function(
    context,
    gene_data,
    palette_colors,
    palette_values,
    show_histology_image = FALSE,
    point_size_no_image = 1.1,
    point_size_with_image = 0.80
) {

  colour_limits <- c(
    0,
    gene_data$spot_colour_max
  )

  panels <- lapply(
    context$sample_layout$sample_order,
    function(sample_id) {

      plot_data <-
        gene_data$sample_plot_data[[
          sample_id
        ]]

      plot_data$plot_colour <-
        map_marker_values_to_palette(
          values =
            plot_data$logNormalized_expression,
          palette_colors =
            palette_colors,
          palette_values =
            palette_values,
          limits =
            colour_limits
        )

      create_marker_spatial_panel(
        sample_cache =
          context$sample_cache[[sample_id]],
        plot_data = plot_data,
        value_column =
          "logNormalized_expression",
        colour_column =
          "plot_colour",
        target_gene =
          gene_data$target_gene,
        marker_clusters =
          gene_data$marker_clusters,
        show_histology_image =
          show_histology_image,
        emphasize_marker_clusters =
          FALSE,
        aggregated =
          FALSE,
        point_size_no_image =
          point_size_no_image,
        point_size_with_image =
          point_size_with_image
      )
    }
  )

  names(panels) <-
    context$sample_layout$sample_order

  legend_panel <- create_marker_colourbar_legend(
    palette_colors = palette_colors,
    palette_values = palette_values,
    colour_max =
      gene_data$spot_colour_max,
    legend_title = paste0(
      gene_data$target_gene,
      " per-spot expression\n",
      "log1p(count / total UMI × ",
      format(
        gene_data$normalization_scale_factor,
        big.mark = ",",
        scientific = FALSE
      ),
      ")"
    )
  )

  arrange_marker_plots_four_columns(
    plot_list_named = panels,
    sample_layout =
      context$sample_layout,
    legend_panel = legend_panel,
    plot_title = paste0(
      gene_data$target_gene,
      " per-spot spatial expression"
    ),
    plot_subtitle = paste0(
      "Columns: Male neurotypical | Male ASD | Female neurotypical | Female ASD",
      " | marker cluster(s): ",
      paste(
        gene_data$marker_clusters,
        collapse = ", "
      ),
      " | common colour range: 0–",
      format_marker_number(
        gene_data$spot_colour_max
      ),
      " | upper limit = ",
      gene_data$upper_colour_quantile *
        100,
      "th percentile of positive values",
      "
",
      gene_data$clustering_summary_label,
      "
",
      gene_data$specificity_summary_label
    )
  )
}


plot_marker_gene_cluster_aggregated_all_clusters <- function(
    context,
    gene_data,
    palette_colors,
    palette_values,
    show_histology_image = FALSE,
    point_size_no_image = 1.1,
    point_size_with_image = 0.80
) {

  colour_limits <- c(
    0,
    gene_data$cluster_colour_max
  )

  panels <- lapply(
    context$sample_layout$sample_order,
    function(sample_id) {

      plot_data <-
        gene_data$sample_plot_data[[
          sample_id
        ]]

      plot_data$plot_colour <-
        map_marker_values_to_palette(
          values =
            plot_data$cluster_logNormalized_expression,
          palette_colors =
            palette_colors,
          palette_values =
            palette_values,
          limits =
            colour_limits
        )

      create_marker_spatial_panel(
        sample_cache =
          context$sample_cache[[sample_id]],
        plot_data = plot_data,
        value_column =
          "cluster_logNormalized_expression",
        colour_column =
          "plot_colour",
        target_gene =
          gene_data$target_gene,
        marker_clusters =
          gene_data$marker_clusters,
        show_histology_image =
          show_histology_image,
        emphasize_marker_clusters =
          FALSE,
        aggregated =
          TRUE,
        point_size_no_image =
          point_size_no_image,
        point_size_with_image =
          point_size_with_image
      )
    }
  )

  names(panels) <-
    context$sample_layout$sample_order

  legend_panel <- create_marker_colourbar_legend(
    palette_colors = palette_colors,
    palette_values = palette_values,
    colour_max =
      gene_data$cluster_colour_max,
    legend_title = paste0(
      gene_data$target_gene,
      " sample × cluster expression\n",
      "log1p(sum count / sum UMI × ",
      format(
        gene_data$normalization_scale_factor,
        big.mark = ",",
        scientific = FALSE
      ),
      ")"
    )
  )

  arrange_marker_plots_four_columns(
    plot_list_named = panels,
    sample_layout =
      context$sample_layout,
    legend_panel = legend_panel,
    plot_title = paste0(
      gene_data$target_gene,
      " cluster-aggregated spatial expression"
    ),
    plot_subtitle = paste0(
      "One normalized value per sample × cluster; all spots in a cluster share that value",
      " | marker cluster(s): ",
      paste(
        gene_data$marker_clusters,
        collapse = ", "
      ),
      " | common colour range: 0–",
      format_marker_number(
        gene_data$cluster_colour_max
      ),
      "
",
      gene_data$clustering_summary_label,
      "
",
      gene_data$specificity_summary_label
    )
  )
}


plot_marker_gene_per_spot_target_vs_other <- function(
    context,
    gene_data,
    target_palette_colors,
    other_palette_colors,
    palette_values,
    show_histology_image = FALSE,
    point_size_no_image = 1.1,
    point_size_with_image = 0.80
) {

  colour_limits <- c(
    0,
    gene_data$spot_colour_max
  )

  panels <- lapply(
    context$sample_layout$sample_order,
    function(sample_id) {

      plot_data <-
        gene_data$sample_plot_data[[
          sample_id
        ]]

      target_colours <-
        map_marker_values_to_palette(
          values =
            plot_data$logNormalized_expression,
          palette_colors =
            target_palette_colors,
          palette_values =
            palette_values,
          limits =
            colour_limits
        )

      other_colours <-
        map_marker_values_to_palette(
          values =
            plot_data$logNormalized_expression,
          palette_colors =
            other_palette_colors,
          palette_values =
            palette_values,
          limits =
            colour_limits
        )

      plot_data$plot_colour <- ifelse(
        plot_data$is_marker_cluster,
        target_colours,
        other_colours
      )

      create_marker_spatial_panel(
        sample_cache =
          context$sample_cache[[sample_id]],
        plot_data = plot_data,
        value_column =
          "logNormalized_expression",
        colour_column =
          "plot_colour",
        target_gene =
          gene_data$target_gene,
        marker_clusters =
          gene_data$marker_clusters,
        show_histology_image =
          show_histology_image,
        emphasize_marker_clusters =
          TRUE,
        aggregated =
          FALSE,
        point_size_no_image =
          point_size_no_image,
        point_size_with_image =
          point_size_with_image
      )
    }
  )

  names(panels) <-
    context$sample_layout$sample_order

  target_legend <-
    create_marker_colourbar_legend(
      palette_colors =
        target_palette_colors,
      palette_values =
        palette_values,
      colour_max =
        gene_data$spot_colour_max,
      legend_title = paste0(
        "Marker cluster(s) ",
        paste(
          gene_data$marker_clusters,
          collapse = ", "
        ),
        " — ",
        gene_data$target_gene,
        " per-spot expression"
      ),
      legend_bar_width_mm = 34
    )

  other_legend <-
    create_marker_colourbar_legend(
      palette_colors =
        other_palette_colors,
      palette_values =
        palette_values,
      colour_max =
        gene_data$spot_colour_max,
      legend_title = paste0(
        "Other clusters — ",
        gene_data$target_gene,
        " per-spot expression"
      ),
      legend_bar_width_mm = 34
    )

  legend_panel <- target_legend +
    other_legend +
    patchwork::plot_layout(
      ncol = 2
    )

  arrange_marker_plots_four_columns(
    plot_list_named = panels,
    sample_layout =
      context$sample_layout,
    legend_panel = legend_panel,
    plot_title = paste0(
      gene_data$target_gene,
      " per-spot expression: marker cluster(s) versus other clusters"
    ),
    plot_subtitle = paste0(
      "Marker cluster(s) shown with the red palette: ",
      paste(
        gene_data$marker_clusters,
        collapse = ", "
      ),
      " | all other clusters shown with the green palette",
      " | both palettes use the same expression limits: 0–",
      format_marker_number(
        gene_data$spot_colour_max
      ),
      "
",
      gene_data$clustering_summary_label,
      "
",
      gene_data$specificity_summary_label
    )
  )
}


plot_marker_gene_cluster_aggregated_target_vs_other <- function(
    context,
    gene_data,
    target_palette_colors,
    other_palette_colors,
    palette_values,
    show_histology_image = FALSE,
    point_size_no_image = 1.1,
    point_size_with_image = 0.80
) {

  colour_limits <- c(
    0,
    gene_data$cluster_colour_max
  )

  panels <- lapply(
    context$sample_layout$sample_order,
    function(sample_id) {

      plot_data <-
        gene_data$sample_plot_data[[
          sample_id
        ]]

      target_colours <-
        map_marker_values_to_palette(
          values =
            plot_data$cluster_logNormalized_expression,
          palette_colors =
            target_palette_colors,
          palette_values =
            palette_values,
          limits =
            colour_limits
        )

      other_colours <-
        map_marker_values_to_palette(
          values =
            plot_data$cluster_logNormalized_expression,
          palette_colors =
            other_palette_colors,
          palette_values =
            palette_values,
          limits =
            colour_limits
        )

      plot_data$plot_colour <- ifelse(
        plot_data$is_marker_cluster,
        target_colours,
        other_colours
      )

      create_marker_spatial_panel(
        sample_cache =
          context$sample_cache[[sample_id]],
        plot_data = plot_data,
        value_column =
          "cluster_logNormalized_expression",
        colour_column =
          "plot_colour",
        target_gene =
          gene_data$target_gene,
        marker_clusters =
          gene_data$marker_clusters,
        show_histology_image =
          show_histology_image,
        emphasize_marker_clusters =
          TRUE,
        aggregated =
          TRUE,
        point_size_no_image =
          point_size_no_image,
        point_size_with_image =
          point_size_with_image
      )
    }
  )

  names(panels) <-
    context$sample_layout$sample_order

  target_legend <-
    create_marker_colourbar_legend(
      palette_colors =
        target_palette_colors,
      palette_values =
        palette_values,
      colour_max =
        gene_data$cluster_colour_max,
      legend_title = paste0(
        "Marker cluster(s) ",
        paste(
          gene_data$marker_clusters,
          collapse = ", "
        ),
        " — sample × cluster expression"
      ),
      legend_bar_width_mm = 34
    )

  other_legend <-
    create_marker_colourbar_legend(
      palette_colors =
        other_palette_colors,
      palette_values =
        palette_values,
      colour_max =
        gene_data$cluster_colour_max,
      legend_title =
        "Other clusters — sample × cluster expression",
      legend_bar_width_mm = 34
    )

  legend_panel <- target_legend +
    other_legend +
    patchwork::plot_layout(
      ncol = 2
    )

  arrange_marker_plots_four_columns(
    plot_list_named = panels,
    sample_layout =
      context$sample_layout,
    legend_panel = legend_panel,
    plot_title = paste0(
      gene_data$target_gene,
      " cluster-aggregated expression: marker cluster(s) versus other clusters"
    ),
    plot_subtitle = paste0(
      "One normalized value per sample × cluster",
      " | marker cluster(s) shown in red: ",
      paste(
        gene_data$marker_clusters,
        collapse = ", "
      ),
      " | other clusters shown in green",
      " | common limits: 0–",
      format_marker_number(
        gene_data$cluster_colour_max
      ),
      "
",
      gene_data$clustering_summary_label,
      "
",
      gene_data$specificity_summary_label
    )
  )
}


plot_marker_gene_barplot_per_cluster <- function(
    context,
    gene_data,
    target_cluster_color = "#B30000",
    other_cluster_color = "#1B9E77",
    cluster_names = NULL,
    bar_width = 0.68,
    point_size = 3.0,
    jitter_width = 0.08,
    errorbar_width = 0.18
) {

  cluster_levels <-
    context$cluster_levels

  sample_cluster_data <-
    gene_data$sample_cluster_summary

  complete_grid <- expand.grid(
    sample_ID =
      context$sample_layout$sample_order,
    cluster =
      cluster_levels,
    stringsAsFactors = FALSE
  )

  sample_cluster_data <- merge(
    x = complete_grid,
    y = sample_cluster_data,
    by = c(
      "sample_ID",
      "cluster"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  sample_cluster_data$cluster <- factor(
    sample_cluster_data$cluster,
    levels = cluster_levels
  )

  sample_cluster_data$is_marker_cluster <-
    as.character(
      sample_cluster_data$cluster
    ) %in%
    gene_data$marker_clusters

  sample_cluster_data$cluster_role <- ifelse(
    sample_cluster_data$is_marker_cluster,
    "Marker cluster",
    "Other cluster"
  )

  cluster_summary <- dplyr::as_tibble(
    sample_cluster_data
  ) |>
    dplyr::group_by(
      cluster
    ) |>
    dplyr::summarise(
      n_samples_with_cluster =
        sum(
          !is.na(
            cluster_logNormalized_expression
          )
        ),
      total_spots =
        sum(
          n_spots,
          na.rm = TRUE
        ),
      positive_spots =
        sum(
          n_positive_spots,
          na.rm = TRUE
        ),
      pooled_percent_positive =
        if (
          sum(
            n_spots,
            na.rm = TRUE
          ) > 0
        ) {
          100 *
            sum(
              n_positive_spots,
              na.rm = TRUE
            ) /
            sum(
              n_spots,
              na.rm = TRUE
            )
        } else {
          NA_real_
        },
      mean_expression =
        mean(
          cluster_logNormalized_expression,
          na.rm = TRUE
        ),
      sd_expression =
        stats::sd(
          cluster_logNormalized_expression,
          na.rm = TRUE
        ),
      se_expression =
        sd_expression /
        sqrt(
          n_samples_with_cluster
        ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      cluster_id =
        as.character(cluster),
      is_marker_cluster =
        cluster_id %in%
        gene_data$marker_clusters,
      cluster_role =
        ifelse(
          is_marker_cluster,
          "Marker cluster",
          "Other cluster"
        ),
      ymin =
        pmax(
          0,
          mean_expression -
            se_expression
        ),
      ymax =
        mean_expression +
        se_expression
    )

  role_colors <- c(
    "Marker cluster" =
      target_cluster_color,
    "Other cluster" =
      other_cluster_color
  )

  cluster_labels <- vapply(
    cluster_levels,
    function(cluster_id) {

      cluster_row <- cluster_summary[
        cluster_summary$cluster_id ==
          cluster_id,
        ,
        drop = FALSE
      ]

      displayed_name <- ""

      if (
        !is.null(cluster_names) &&
          cluster_id %in%
          names(cluster_names) &&
          !is.na(
            cluster_names[[cluster_id]]
          ) &&
          nzchar(
            cluster_names[[cluster_id]]
          )
      ) {
        displayed_name <- paste0(
          "<br>",
          cluster_names[[cluster_id]]
        )
      }

      label_color <- if (
        cluster_id %in%
          gene_data$marker_clusters
      ) {
        target_cluster_color
      } else {
        other_cluster_color
      }

      paste0(
        "<span style='color:",
        label_color,
        ";font-weight:bold'>Cluster ",
        cluster_id,
        displayed_name,
        "</span>",
        "<br>",
        if (
          is.finite(
            cluster_row$pooled_percent_positive
          )
        ) {
          formatC(
            cluster_row$pooled_percent_positive,
            format = "f",
            digits = 1
          )
        } else {
          "NA"
        },
        "% ",
        gene_data$target_gene,
        "+ spots"
      )
    },
    character(1)
  )

  names(cluster_labels) <-
    cluster_levels

  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = cluster_summary,
      ggplot2::aes(
        x = cluster,
        y = mean_expression,
        colour = cluster_role
      ),
      width = bar_width,
      fill = NA,
      linewidth = 0.9
    ) +
    ggplot2::geom_errorbar(
      data = cluster_summary,
      ggplot2::aes(
        x = cluster,
        ymin = ymin,
        ymax = ymax,
        colour = cluster_role
      ),
      width = errorbar_width,
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = sample_cluster_data[
        !is.na(
          sample_cluster_data$cluster_logNormalized_expression
        ),
        ,
        drop = FALSE
      ],
      ggplot2::aes(
        x = cluster,
        y = cluster_logNormalized_expression,
        fill = cluster_role
      ),
      position = ggplot2::position_jitter(
        width = jitter_width,
        height = 0,
        seed = 123
      ),
      shape = 21,
      size = point_size,
      stroke = 0.7,
      colour = "black"
    ) +
    ggplot2::scale_colour_manual(
      values = role_colors,
      breaks = names(role_colors)
    ) +
    ggplot2::scale_fill_manual(
      values = role_colors,
      breaks = names(role_colors)
    ) +
    ggplot2::scale_x_discrete(
      limits = cluster_levels,
      labels = cluster_labels,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.10
        )
      )
    ) +
    ggplot2::labs(
      title = paste0(
        gene_data$target_gene,
        " expression across spatial clusters"
      ),
      subtitle = paste0(
        "Empty bars: mean sample × cluster expression",
        " | dots: sample × cluster values",
        " | whiskers: mean ± SE",
        "\nMarker cluster(s): ",
        paste(
          gene_data$marker_clusters,
          collapse = ", "
        ),
        " | below each cluster: pooled % of spots with detectable ",
        gene_data$target_gene,
        " expression",
        "
",
        gene_data$clustering_summary_label,
        "
",
        gene_data$specificity_summary_label
      ),
      x = NULL,
      y = paste0(
        gene_data$target_gene,
        " cluster-aggregated expression",
        "\nlog1p(sum count / sum UMI × ",
        format(
          gene_data$normalization_scale_factor,
          big.mark = ",",
          scientific = FALSE
        ),
        ")"
      ),
      colour = NULL,
      fill = NULL
    ) +
    ggplot2::theme_classic(
      base_family = "DejaVu Sans"
    ) +
    ggplot2::theme(
      plot.title =
        ggplot2::element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
      plot.subtitle =
        ggplot2::element_text(
          size = 9.5,
          hjust = 0.5,
          lineheight = 1.08,
          margin = ggplot2::margin(
            b = 12
          )
        ),
      axis.title.y =
        ggplot2::element_text(
          size = 11,
          face = "bold",
          margin = ggplot2::margin(
            r = 10
          )
        ),
      axis.text.x =
        ggtext::element_markdown(
          size = 8,
          lineheight = 0.95,
          angle = 45,
          hjust = 1,
          vjust = 1,
          margin = ggplot2::margin(
            t = 6
          )
        ),
      axis.text.y =
        ggplot2::element_text(
          size = 9
        ),
      axis.line =
        ggplot2::element_line(
          linewidth = 0.7
        ),
      axis.ticks =
        ggplot2::element_line(
          linewidth = 0.7
        ),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text =
        ggplot2::element_text(
          size = 10
        ),
      plot.margin =
        ggplot2::margin(
          t = 10,
          r = 16,
          b = 12,
          l = 16
        )
    ) +
    ggplot2::guides(
      colour = "none",
      fill = ggplot2::guide_legend(
        override.aes = list(
          shape = 21,
          size = 4,
          colour = "black"
        )
      )
    )
}


build_marker_cluster_file_label <- function(
    marker_clusters
) {

  marker_clusters <- sort_marker_cluster_ids(
    marker_clusters
  )

  safe_cluster_ids <- gsub(
    pattern = "[^A-Za-z0-9._-]+",
    replacement = "-",
    x = marker_clusters
  )

  paste0(
    "markerClusters_",
    paste(
      safe_cluster_ids,
      collapse = "-"
    )
  )
}


# ==============================================================================
# 6. Output orchestration
# ==============================================================================

get_marker_gene_output_files <- function(
    output_root_directory,
    target_gene,
    marker_clusters
) {

  marker_cluster_label <-
    build_marker_cluster_file_label(
      marker_clusters
    )

  gene_output_directory <- file.path(
    output_root_directory,
    paste0(
      target_gene,
      "_",
      marker_cluster_label
    )
  )

  file_prefix <- paste0(
    target_gene,
    "_",
    marker_cluster_label
  )

  list(
    gene_output_directory =
      gene_output_directory,

    per_spot_all_clusters_pdf =
      file.path(
        gene_output_directory,
        paste0(
          "01_",
          file_prefix,
          "_perSpotExpression_allClusters.pdf"
        )
      ),

    cluster_aggregated_all_clusters_pdf =
      file.path(
        gene_output_directory,
        paste0(
          "02_",
          file_prefix,
          "_clusterAggregatedExpression_allClusters.pdf"
        )
      ),

    per_spot_target_vs_other_pdf =
      file.path(
        gene_output_directory,
        paste0(
          "03_",
          file_prefix,
          "_perSpotExpression_markerClustersRed_otherClustersGreen.pdf"
        )
      ),

    cluster_aggregated_target_vs_other_pdf =
      file.path(
        gene_output_directory,
        paste0(
          "04_",
          file_prefix,
          "_clusterAggregatedExpression_markerClustersRed_otherClustersGreen.pdf"
        )
      ),

    cluster_barplot_pdf =
      file.path(
        gene_output_directory,
        paste0(
          "05_",
          file_prefix,
          "_clusterExpression_barplot.pdf"
        )
      ),

    sample_summary_tsv =
      file.path(
        gene_output_directory,
        paste0(
          "06_",
          file_prefix,
          "_sampleSummary.tsv"
        )
      ),

    sample_cluster_summary_tsv =
      file.path(
        gene_output_directory,
        paste0(
          "07_",
          file_prefix,
          "_sampleClusterSummary.tsv"
        )
      )
  )
}

marker_gene_outputs_complete <- function(
    output_files
) {

  expected_files <- unlist(
    output_files[
      setdiff(
        names(output_files),
        "gene_output_directory"
      )
    ],
    use.names = FALSE
  )

  all(
    file.exists(
      expected_files
    )
  )
}


save_marker_gene_visualizations <- function(
    context,
    target_gene,
    marker_clusters,
    marker_metrics_table = NULL,
    clustering_parameters_label = NULL,
    output_root_directory,
    normalization_scale_factor = 10000,
    upper_colour_quantile = 0.99,
    red_palette_colors = c(
      "#D9D9D9",
      "#FEE5D9",
      "#FCAE91",
      "#FB6A4A",
      "#DE2D26",
      "#A50F15",
      "#67000D"
    ),
    green_palette_colors = c(
      "#D9D9D9",
      "#E5F5E0",
      "#A1D99B",
      "#74C476",
      "#31A354",
      "#006D2C",
      "#00441B"
    ),
    palette_values = c(
      0.000,
      0.006,
      0.060,
      0.180,
      0.400,
      0.700,
      1.000
    ),
    target_cluster_color = "#B30000",
    other_cluster_color = "#1B9E77",
    cluster_names = NULL,
    show_histology_image = FALSE,
    point_size_no_image = 1.1,
    point_size_with_image = 0.80,
    spatial_pdf_width = 18,
    spatial_pdf_height = 23,
    barplot_pdf_width = 15,
    barplot_pdf_height = 8,
    skip_completed = TRUE
) {

  output_files <- get_marker_gene_output_files(
    output_root_directory =
      output_root_directory,
    target_gene = target_gene,
    marker_clusters =
      marker_clusters
  )

  if (
    isTRUE(skip_completed) &&
      marker_gene_outputs_complete(
        output_files
      )
  ) {
    message(
      "Skipping completed gene: ",
      target_gene
    )

    return(
      invisible(output_files)
    )
  }

  dir.create(
    output_files$gene_output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  gene_data <- prepare_marker_gene_data(
    context = context,
    target_gene = target_gene,
    marker_clusters = marker_clusters,
    marker_metrics_table =
      marker_metrics_table,
    clustering_parameters_label =
      clustering_parameters_label,
    normalization_scale_factor =
      normalization_scale_factor,
    upper_colour_quantile =
      upper_colour_quantile
  )

  write_marker_tsv(
    gene_data$sample_summary,
    output_files$sample_summary_tsv
  )

  write_marker_tsv(
    gene_data$sample_cluster_summary,
    output_files$sample_cluster_summary_tsv
  )

  plot_1 <-
    plot_marker_gene_per_spot_all_clusters(
      context = context,
      gene_data = gene_data,
      palette_colors =
        red_palette_colors,
      palette_values =
        palette_values,
      show_histology_image =
        show_histology_image,
      point_size_no_image =
        point_size_no_image,
      point_size_with_image =
        point_size_with_image
    )

  plot_2 <-
    plot_marker_gene_cluster_aggregated_all_clusters(
      context = context,
      gene_data = gene_data,
      palette_colors =
        red_palette_colors,
      palette_values =
        palette_values,
      show_histology_image =
        show_histology_image,
      point_size_no_image =
        point_size_no_image,
      point_size_with_image =
        point_size_with_image
    )

  plot_3 <-
    plot_marker_gene_per_spot_target_vs_other(
      context = context,
      gene_data = gene_data,
      target_palette_colors =
        red_palette_colors,
      other_palette_colors =
        green_palette_colors,
      palette_values =
        palette_values,
      show_histology_image =
        show_histology_image,
      point_size_no_image =
        point_size_no_image,
      point_size_with_image =
        point_size_with_image
    )

  plot_4 <-
    plot_marker_gene_cluster_aggregated_target_vs_other(
      context = context,
      gene_data = gene_data,
      target_palette_colors =
        red_palette_colors,
      other_palette_colors =
        green_palette_colors,
      palette_values =
        palette_values,
      show_histology_image =
        show_histology_image,
      point_size_no_image =
        point_size_no_image,
      point_size_with_image =
        point_size_with_image
    )

  plot_5 <-
    plot_marker_gene_barplot_per_cluster(
      context = context,
      gene_data = gene_data,
      target_cluster_color =
        target_cluster_color,
      other_cluster_color =
        other_cluster_color,
      cluster_names =
        cluster_names
    )

  pdf_device <- if (
    capabilities("cairo")
  ) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  ggplot2::ggsave(
    filename =
      output_files$per_spot_all_clusters_pdf,
    plot = plot_1,
    device = pdf_device,
    width = spatial_pdf_width,
    height = spatial_pdf_height,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )

  ggplot2::ggsave(
    filename =
      output_files$cluster_aggregated_all_clusters_pdf,
    plot = plot_2,
    device = pdf_device,
    width = spatial_pdf_width,
    height = spatial_pdf_height,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )

  ggplot2::ggsave(
    filename =
      output_files$per_spot_target_vs_other_pdf,
    plot = plot_3,
    device = pdf_device,
    width = spatial_pdf_width,
    height = spatial_pdf_height,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )

  ggplot2::ggsave(
    filename =
      output_files$cluster_aggregated_target_vs_other_pdf,
    plot = plot_4,
    device = pdf_device,
    width = spatial_pdf_width,
    height = spatial_pdf_height,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )

  ggplot2::ggsave(
    filename =
      output_files$cluster_barplot_pdf,
    plot = plot_5,
    device = pdf_device,
    width = barplot_pdf_width,
    height = barplot_pdf_height,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )

  if (
    !marker_gene_outputs_complete(
      output_files
    )
  ) {
    stop(
      "Not all expected outputs were created for gene: ",
      target_gene,
      call. = FALSE
    )
  }

  invisible(
    output_files
  )
}
