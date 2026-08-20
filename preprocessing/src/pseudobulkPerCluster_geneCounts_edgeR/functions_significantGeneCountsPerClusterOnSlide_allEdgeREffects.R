# ==============================================================================
# functions_significantGeneCountsPerClusterOnSlide_allEdgeREffects.R
#
# Functions for one 3 x 3 spatial overview of significant cluster-gene counts
# from the cluster-specific pseudobulk edgeR Group x Sex model.
#
# Rows:
#   1. Group x Sex interaction
#   2. Main FMT donor-group effect
#   3. Main sex effect
#
# Columns:
#   1. all significant results
#   2. positive-logFC results
#   3. negative-logFC results
#
# Two scale modes are exported by the workflow:
#   - independent: every panel uses its own maximum count;
#   - shared_by_direction: one maximum for all three ALL panels, one maximum
#     for all three POSITIVE panels and one maximum for all three NEGATIVE panels.
#
# Filtering is shared across all three edgeR effects:
#   FDR < fdr_threshold
#   abs(logFC) >= abs_log2fc_threshold
#   optional Max mean %-positive interval filter:
#     lower bound: max_mean_percent >= min_max_mean_percent
#     upper bound: max_mean_percent < max_max_mean_percent
#     except upper bound 100, which is treated as <= 100.
#
# When both bounds are NULL, the percent-positive filter is disabled.
# Standard variants:
#   all percentages: lower = NULL, upper = NULL
#   0 to <25%:      lower = 0,    upper = 25
#   25 to 100%:     lower = 25,   upper = 100
#
# Legend design in this version:
#   - three legends total: one for Interaction, one for FMT donor group,
#     one for Sex;
#   - every legend is shared within one biological-effect row and is
#     displayed in a 3-column layout with larger one-line labels;
#   - four dots per cluster: anatomical colour, All, UP/positive, DOWN/negative;
#   - legend entry reports All / UP / DOWN cluster-gene counts and the current
#     anatomical cluster label used in the heatmap workflow.
#
# Intended server-side filename remains:
#   functions_significantGeneCountsPerClusterOnSlide_allEdgeREffects.R
# ==============================================================================


# ==============================================================================
# V17 diagnostic marker
# ==============================================================================

SPATIAL_EDGER_ALL_EFFECTS_FUNCTIONS_VERSION <- "V17_HALF_DOT_TEXT_GAP_2026-08-14"
message(
  "Loaded functions_significantGeneCountsPerClusterOnSlide_allEdgeREffects: ",
  SPATIAL_EDGER_ALL_EFFECTS_FUNCTIONS_VERSION
)

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
  formatC(as.numeric(x), format = "f", digits = as.integer(digits))
}

format_spatial_edger_integer <- function(x) {
  format(
    as.integer(x),
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

sanitize_spatial_edger_filename_component <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "-", as.character(x))
}

format_spatial_edger_percent_code <- function(x) {

  value <- format(
    as.numeric(x),
    scientific = FALSE,
    trim = TRUE,
    digits = 6
  )

  gsub("\\.", "p", value)
}


build_spatial_edger_parameter_label <- function(
  fdr_threshold,
  abs_log2fc_threshold,
  min_max_mean_percent = NULL,
  max_max_mean_percent = NULL
) {

  fdr_code <- sprintf(
    "%03d",
    as.integer(round(as.numeric(fdr_threshold) * 100))
  )

  logfc_code <- sprintf(
    "%02d",
    as.integer(round(as.numeric(abs_log2fc_threshold) * 10))
  )

  percent_code <- if (
    is.null(min_max_mean_percent) &&
      is.null(max_max_mean_percent)
  ) {
    "All"
  } else if (
    !is.null(min_max_mean_percent) &&
      !is.null(max_max_mean_percent)
  ) {

    lower_code <- format_spatial_edger_percent_code(
      min_max_mean_percent
    )

    upper_code <- format_spatial_edger_percent_code(
      max_max_mean_percent
    )

    if (as.numeric(max_max_mean_percent) >= 100) {
      paste0(
        lower_code,
        "to",
        upper_code
      )
    } else {
      paste0(
        lower_code,
        "toLt",
        upper_code
      )
    }

  } else if (!is.null(min_max_mean_percent)) {

    paste0(
      "Ge",
      format_spatial_edger_percent_code(
        min_max_mean_percent
      )
    )

  } else {

    upper_code <- format_spatial_edger_percent_code(
      max_max_mean_percent
    )

    if (as.numeric(max_max_mean_percent) >= 100) {
      paste0(
        "Le",
        upper_code
      )
    } else {
      paste0(
        "Lt",
        upper_code
      )
    }
  }

  paste0(
    "FDR",
    fdr_code,
    "_absLog2FC",
    logfc_code,
    "_maxMeanPct",
    percent_code
  )
}


format_spatial_edger_filter_label <- function(
  fdr_threshold,
  abs_log2fc_threshold,
  min_max_mean_percent = NULL,
  max_max_mean_percent = NULL
) {

  percent_text <- if (
    is.null(min_max_mean_percent) &&
      is.null(max_max_mean_percent)
  ) {

    "Max mean % filter: none"

  } else if (
    !is.null(min_max_mean_percent) &&
      !is.null(max_max_mean_percent)
  ) {

    upper_operator <- if (
      as.numeric(max_max_mean_percent) >= 100
    ) {
      "<="
    } else {
      "<"
    }

    paste0(
      format(as.numeric(min_max_mean_percent), trim = TRUE),
      "% <= Max mean % ",
      upper_operator,
      " ",
      format(as.numeric(max_max_mean_percent), trim = TRUE),
      "%"
    )

  } else if (!is.null(min_max_mean_percent)) {

    paste0(
      "Max mean % >= ",
      format(as.numeric(min_max_mean_percent), trim = TRUE),
      "%"
    )

  } else {

    upper_operator <- if (
      as.numeric(max_max_mean_percent) >= 100
    ) {
      "<="
    } else {
      "<"
    }

    paste0(
      "Max mean % ",
      upper_operator,
      " ",
      format(as.numeric(max_max_mean_percent), trim = TRUE),
      "%"
    )
  }

  paste0(
    "FDR < ",
    format(as.numeric(fdr_threshold), trim = TRUE),
    " | |log2FC| >= ",
    format(as.numeric(abs_log2fc_threshold), trim = TRUE),
    " | ",
    percent_text
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
    x_limits = c(x_center - half_side, x_center + half_side),
    y_limits = c(y_center - half_side, y_center + half_side)
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


# ==============================================================================
# 2. Load input objects
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
          get(object_name, envir = load_environment, inherits = FALSE),
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
    object = get(selected_name, envir = load_environment, inherits = FALSE),
    object_name = selected_name
  )
}

load_all_edger_effects_for_significant_gene_count_plot <- function(
  edger_results_rdata_file,
  interaction_test_id = "Interaction",
  group_test_id = "Overall_Group_ASD_vs_Neurotypical",
  sex_test_id = "Overall_Sex_Female_vs_Male"
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
  required_objects <- c(
    "edgeR_perCluster_combinedResults",
    "edgeR_perCluster_testDefinitions"
  )
  missing_objects <- setdiff(required_objects, loaded_names)
  if (length(missing_objects) > 0L) {
    stop(
      "The edgeR RData file is missing required object(s): ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }
  combined_results <- get(
    "edgeR_perCluster_combinedResults",
    envir = result_environment,
    inherits = FALSE
  )
  test_definitions <- get(
    "edgeR_perCluster_testDefinitions",
    envir = result_environment,
    inherits = FALSE
  )
  requested_test_ids <- c(
    interaction = interaction_test_id,
    group = group_test_id,
    sex = sex_test_id
  )
  missing_test_ids <- requested_test_ids[
    !requested_test_ids %in% names(combined_results)
  ]
  if (length(missing_test_ids) > 0L) {
    stop(
      "Requested edgeR test(s) are absent: ",
      paste(missing_test_ids, collapse = ", "),
      "\nAvailable tests: ",
      paste(names(combined_results), collapse = ", "),
      call. = FALSE
    )
  }
  test_definition_subset <- test_definitions |>
    dplyr::filter(.data$test_id %in% unname(requested_test_ids))
  list(
    results = list(
      interaction = tibble::as_tibble(combined_results[[interaction_test_id]]),
      group = tibble::as_tibble(combined_results[[group_test_id]]),
      sex = tibble::as_tibble(combined_results[[sex_test_id]])
    ),
    test_ids = requested_test_ids,
    test_definitions = test_definition_subset,
    loaded_names = loaded_names
  )
}


# ==============================================================================
# 3. Prepare selected Visium section
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
  required_columns <- c(sample_id_column, cluster_column)
  missing_columns <- setdiff(required_columns, colnames(metadata_table))
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
      image_cells <- SeuratObject::Cells(seurat_object[[image_name]])
      shared_cells <- intersect(image_cells, rownames(metadata_table))
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
  image_mapping <- do.call(rbind, image_mapping_rows)
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
  image_cells <- SeuratObject::Cells(seurat_object[[image_name]])
  coordinate_table <- SeuratObject::GetTissueCoordinates(
    object = seurat_object,
    image = image_name,
    scale = image_scale
  )
  coordinate_table <- standardize_spatial_edger_coordinates(coordinate_table)
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
  spot_metadata <- metadata_table[common_cells, , drop = FALSE]
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
# 4. Effect definitions and max-mean-percent calculation
# ==============================================================================

get_spatial_edger_effect_definitions <- function(
  interaction_test_id = "Interaction",
  group_test_id = "Overall_Group_ASD_vs_Neurotypical",
  sex_test_id = "Overall_Sex_Female_vs_Male"
) {
  tibble::tibble(
    effect_key = c("interaction", "group", "sex"),
    effect_order = c(1L, 2L, 3L),
    test_id = c(interaction_test_id, group_test_id, sex_test_id),
    row_label = c(
      "Group x Sex interaction",
      "FMT donor group",
      "Sex"
    ),
    positive_label = c(
      "ASD effect higher in Female",
      "UP in ASD",
      "UP in Female"
    ),
    negative_label = c(
      "ASD effect higher in Male",
      "DOWN in ASD",
      "DOWN in Female"
    ),
    positive_interpretation = c(
      "positive interaction logFC",
      "ASD > Neurotypical",
      "Female > Male"
    ),
    negative_interpretation = c(
      "negative interaction logFC",
      "ASD < Neurotypical",
      "Female < Male"
    )
  )
}

resolve_spatial_edger_percent_columns <- function(result_table, effect_key) {

  exact_columns <- switch(
    effect_key,
    interaction = c(
      "mean_percent_positive_spots_Male_Neurotypical",
      "mean_percent_positive_spots_Male_ASD",
      "mean_percent_positive_spots_Female_Neurotypical",
      "mean_percent_positive_spots_Female_ASD"
    ),
    group = c(
      "mean_percent_positive_spots_Neurotypical",
      "mean_percent_positive_spots_ASD"
    ),
    sex = c(
      "mean_percent_positive_spots_Male",
      "mean_percent_positive_spots_Female"
    ),
    stop("Unknown effect_key: ", effect_key, call. = FALSE)
  )

  if (all(exact_columns %in% colnames(result_table))) {
    return(exact_columns)
  }

  available_mean_percent_columns <- grep(
    "^mean_percent_positive_spots_",
    colnames(result_table),
    value = TRUE
  )

  stop(
    "Could not resolve the exact percent-positive columns required for effect '",
    effect_key,
    "'. Expected: ",
    paste(exact_columns, collapse = ", "),
    "\nAvailable mean percent-positive columns: ",
    if (length(available_mean_percent_columns) == 0L) {
      "<none>"
    } else {
      paste(available_mean_percent_columns, collapse = ", ")
    },
    "\nThe Max mean % filter is defined from the effect-specific groups and will not use marginal fallback columns.",
    call. = FALSE
  )
}

calculate_spatial_edger_max_mean_percent <- function(
  result_table,
  percent_columns
) {
  percent_matrix <- as.matrix(
    result_table[, percent_columns, drop = FALSE]
  )
  storage.mode(percent_matrix) <- "numeric"
  max_mean_percent <- apply(
    percent_matrix,
    1L,
    function(values) {
      finite_values <- values[is.finite(values)]
      if (length(finite_values) == 0L) {
        return(NA_real_)
      }
      max(finite_values)
    }
  )
  as.numeric(max_mean_percent)
}

validate_spatial_edger_result_table <- function(result_table, effect_key) {
  required_columns <- c(
    "cluster_id",
    "ensembl_gene_id",
    "gene",
    "logFC",
    "PValue",
    "FDR"
  )
  missing_columns <- setdiff(required_columns, colnames(result_table))
  if (length(missing_columns) > 0L) {
    stop(
      "edgeR result table for '",
      effect_key,
      "' is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


# ==============================================================================
# 5. Select and summarize significant results
# ==============================================================================

prepare_spatial_edger_effect_results <- function(
  result_table,
  effect_definition,
  cluster_levels,
  fdr_threshold,
  abs_log2fc_threshold,
  min_max_mean_percent = NULL,
  max_max_mean_percent = NULL,
  cluster_names = NULL
) {
  effect_key <- effect_definition$effect_key[[1]]
  validate_spatial_edger_result_table(result_table, effect_key)

  # Max mean % is an OPTIONAL interval criterion.
  # Percent-positive columns are inspected only when at least one interval
  # boundary is supplied. When both boundaries are NULL, Max mean % is not
  # calculated and percent-positive columns are not required.
  result_table <- result_table |>
    tibble::as_tibble() |>
    dplyr::mutate(
      cluster_id = as.character(.data$cluster_id)
    )

  percent_filter_enabled <- !(
    is.null(min_max_mean_percent) &&
      is.null(max_max_mean_percent)
  )

  if (!isTRUE(percent_filter_enabled)) {
    # Percent-positive columns are deliberately not inspected when the interval
    # filter is disabled.
    percent_columns <- character(0)
    result_table <- result_table |>
      dplyr::mutate(
        max_mean_percent = NA_real_
      )
  } else {
    percent_columns <- resolve_spatial_edger_percent_columns(
      result_table = result_table,
      effect_key = effect_key
    )
    result_table <- result_table |>
      dplyr::mutate(
        max_mean_percent = calculate_spatial_edger_max_mean_percent(
          result_table = result_table,
          percent_columns = percent_columns
        )
      )
  }
  significant_results <- result_table |>
    dplyr::filter(
      is.finite(.data$FDR),
      is.finite(.data$logFC),
      .data$FDR < .env$fdr_threshold,
      abs(.data$logFC) >= .env$abs_log2fc_threshold
    )
  if (isTRUE(percent_filter_enabled)) {

    significant_results <- significant_results |>
      dplyr::filter(
        is.finite(.data$max_mean_percent)
      )

    if (!is.null(min_max_mean_percent)) {
      significant_results <- significant_results |>
        dplyr::filter(
          .data$max_mean_percent >= .env$min_max_mean_percent
        )
    }

    if (!is.null(max_max_mean_percent)) {

      if (as.numeric(max_max_mean_percent) >= 100) {
        significant_results <- significant_results |>
          dplyr::filter(
            .data$max_mean_percent <= .env$max_max_mean_percent
          )
      } else {
        significant_results <- significant_results |>
          dplyr::filter(
            .data$max_mean_percent < .env$max_max_mean_percent
          )
      }
    }
  }

  significant_results <- significant_results |>
    dplyr::mutate(
      effect_key = .env$effect_key,
      effect_order = effect_definition$effect_order[[1]],
      test_id = effect_definition$test_id[[1]],
      effect_label = effect_definition$row_label[[1]],
      direction = dplyr::case_when(
        .data$logFC > 0 ~ "positive",
        .data$logFC < 0 ~ "negative",
        TRUE ~ "zero"
      ),
      direction_label = dplyr::case_when(
        .data$logFC > 0 ~ effect_definition$positive_label[[1]],
        .data$logFC < 0 ~ effect_definition$negative_label[[1]],
        TRUE ~ "Zero logFC"
      ),
      fdr_threshold = as.numeric(fdr_threshold),
      abs_log2fc_threshold = as.numeric(abs_log2fc_threshold),
      min_max_mean_percent = if (is.null(min_max_mean_percent)) {
        NA_real_
      } else {
        as.numeric(min_max_mean_percent)
      },
      max_max_mean_percent = if (is.null(max_max_mean_percent)) {
        NA_real_
      } else {
        as.numeric(max_max_mean_percent)
      }
    ) |>
    dplyr::arrange(
      as.integer(.data$cluster_id),
      .data$FDR,
      dplyr::desc(abs(.data$logFC))
    )
  cluster_levels <- sort_spatial_edger_cluster_ids(cluster_levels)
  all_cluster_table <- tibble::tibble(cluster_id = cluster_levels)
  if (!is.null(cluster_names)) {
    missing_cluster_names <- setdiff(cluster_levels, names(cluster_names))
    if (length(missing_cluster_names) > 0L) {
      stop(
        "Missing cluster name(s) for cluster ID(s): ",
        paste(missing_cluster_names, collapse = ", "),
        call. = FALSE
      )
    }
    all_cluster_table <- all_cluster_table |>
      dplyr::mutate(
        cluster_name = unname(cluster_names[.data$cluster_id])
      )
  } else {
    all_cluster_table <- all_cluster_table |>
      dplyr::mutate(cluster_name = NA_character_)
  }
  total_counts <- significant_results |>
    dplyr::count(.data$cluster_id, name = "n_all")
  positive_counts <- significant_results |>
    dplyr::filter(.data$direction == "positive") |>
    dplyr::count(.data$cluster_id, name = "n_positive")
  negative_counts <- significant_results |>
    dplyr::filter(.data$direction == "negative") |>
    dplyr::count(.data$cluster_id, name = "n_negative")
  cluster_count_table <- all_cluster_table |>
    dplyr::left_join(total_counts, by = "cluster_id") |>
    dplyr::left_join(positive_counts, by = "cluster_id") |>
    dplyr::left_join(negative_counts, by = "cluster_id") |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("n_all", "n_positive", "n_negative")),
        ~ dplyr::coalesce(as.integer(.x), 0L)
      ),
      effect_key = .env$effect_key,
      effect_order = effect_definition$effect_order[[1]],
      test_id = effect_definition$test_id[[1]],
      effect_label = effect_definition$row_label[[1]],
      fdr_threshold = as.numeric(fdr_threshold),
      abs_log2fc_threshold = as.numeric(abs_log2fc_threshold),
      min_max_mean_percent = if (is.null(min_max_mean_percent)) {
        NA_real_
      } else {
        as.numeric(min_max_mean_percent)
      },
      max_max_mean_percent = if (is.null(max_max_mean_percent)) {
        NA_real_
      } else {
        as.numeric(max_max_mean_percent)
      },
      .before = 1
    )
  if (!all(
    cluster_count_table$n_all ==
      cluster_count_table$n_positive + cluster_count_table$n_negative
  )) {
    stop(
      "Internal count error for effect '",
      effect_key,
      "': n_all != n_positive + n_negative.",
      call. = FALSE
    )
  }
  summary_table <- tibble::tibble(
    effect_key = effect_key,
    effect_order = effect_definition$effect_order[[1]],
    test_id = effect_definition$test_id[[1]],
    effect_label = effect_definition$row_label[[1]],
    positive_label = effect_definition$positive_label[[1]],
    negative_label = effect_definition$negative_label[[1]],
    fdr_threshold = as.numeric(fdr_threshold),
    abs_log2fc_threshold = as.numeric(abs_log2fc_threshold),
    min_max_mean_percent = if (is.null(min_max_mean_percent)) {
      NA_real_
    } else {
      as.numeric(min_max_mean_percent)
    },
    max_max_mean_percent = if (is.null(max_max_mean_percent)) {
      NA_real_
    } else {
      as.numeric(max_max_mean_percent)
    },
    n_all = nrow(significant_results),
    unique_genes_all = dplyr::n_distinct(significant_results$ensembl_gene_id),
    n_positive = sum(significant_results$direction == "positive"),
    unique_genes_positive = dplyr::n_distinct(
      significant_results$ensembl_gene_id[
        significant_results$direction == "positive"
      ]
    ),
    n_negative = sum(significant_results$direction == "negative"),
    unique_genes_negative = dplyr::n_distinct(
      significant_results$ensembl_gene_id[
        significant_results$direction == "negative"
      ]
    ),
    percent_columns = if (length(percent_columns) == 0L) {
      "not_used"
    } else {
      paste(percent_columns, collapse = ";")
    }
  )
  list(
    significant_results = significant_results,
    cluster_count_table = cluster_count_table,
    summary_table = summary_table,
    effect_definition = effect_definition,
    percent_columns = percent_columns
  )
}

prepare_all_spatial_edger_effect_summaries <- function(
  edger_results,
  effect_definitions,
  cluster_levels,
  fdr_threshold,
  abs_log2fc_threshold,
  min_max_mean_percent = NULL,
  max_max_mean_percent = NULL,
  cluster_names = NULL
) {
  output <- lapply(
    effect_definitions$effect_key,
    function(effect_key) {
      effect_definition <- effect_definitions |>
        dplyr::filter(.data$effect_key == .env$effect_key)
      prepare_spatial_edger_effect_results(
        result_table = edger_results[[effect_key]],
        effect_definition = effect_definition,
        cluster_levels = cluster_levels,
        fdr_threshold = fdr_threshold,
        abs_log2fc_threshold = abs_log2fc_threshold,
        min_max_mean_percent = min_max_mean_percent,
        max_max_mean_percent = max_max_mean_percent,
        cluster_names = cluster_names
      )
    }
  )
  names(output) <- effect_definitions$effect_key
  output
}


# ==============================================================================
# 6. Convert effect summaries to one 9-panel count table
# ==============================================================================

build_spatial_edger_nine_panel_count_table <- function(effect_summaries) {
  panel_rows <- lapply(
    names(effect_summaries),
    function(effect_key) {
      effect_summary <- effect_summaries[[effect_key]]
      count_table <- effect_summary$cluster_count_table
      effect_definition <- effect_summary$effect_definition
      all_panel <- count_table |>
        dplyr::transmute(
          effect_key = .data$effect_key,
          effect_order = .data$effect_order,
          test_id = .data$test_id,
          effect_label = .data$effect_label,
          panel_key = "all",
          panel_order = 1L,
          scale_group = "all",
          panel_label = "All significant",
          panel_interpretation = "positive + negative logFC",
          cluster_id = .data$cluster_id,
          cluster_name = .data$cluster_name,
          n_results = .data$n_all,
          fdr_threshold = .data$fdr_threshold,
          abs_log2fc_threshold = .data$abs_log2fc_threshold,
          min_max_mean_percent = .data$min_max_mean_percent,
          max_max_mean_percent = .data$max_max_mean_percent
        )
      positive_panel <- count_table |>
        dplyr::transmute(
          effect_key = .data$effect_key,
          effect_order = .data$effect_order,
          test_id = .data$test_id,
          effect_label = .data$effect_label,
          panel_key = "positive",
          panel_order = 2L,
          scale_group = "positive",
          panel_label = effect_definition$positive_label[[1]],
          panel_interpretation = effect_definition$positive_interpretation[[1]],
          cluster_id = .data$cluster_id,
          cluster_name = .data$cluster_name,
          n_results = .data$n_positive,
          fdr_threshold = .data$fdr_threshold,
          abs_log2fc_threshold = .data$abs_log2fc_threshold,
          min_max_mean_percent = .data$min_max_mean_percent,
          max_max_mean_percent = .data$max_max_mean_percent
        )
      negative_panel <- count_table |>
        dplyr::transmute(
          effect_key = .data$effect_key,
          effect_order = .data$effect_order,
          test_id = .data$test_id,
          effect_label = .data$effect_label,
          panel_key = "negative",
          panel_order = 3L,
          scale_group = "negative",
          panel_label = effect_definition$negative_label[[1]],
          panel_interpretation = effect_definition$negative_interpretation[[1]],
          cluster_id = .data$cluster_id,
          cluster_name = .data$cluster_name,
          n_results = .data$n_negative,
          fdr_threshold = .data$fdr_threshold,
          abs_log2fc_threshold = .data$abs_log2fc_threshold,
          min_max_mean_percent = .data$min_max_mean_percent,
          max_max_mean_percent = .data$max_max_mean_percent
        )
      dplyr::bind_rows(all_panel, positive_panel, negative_panel)
    }
  )
  dplyr::bind_rows(panel_rows) |>
    dplyr::arrange(
      .data$effect_order,
      .data$panel_order,
      as.integer(.data$cluster_id)
    )
}

build_spatial_edger_panel_summary_table <- function(
  nine_panel_count_table,
  effect_summaries
) {
  result_rows <- lapply(
    names(effect_summaries),
    function(effect_key) {
      significant_results <- effect_summaries[[effect_key]]$significant_results
      effect_definition <- effect_summaries[[effect_key]]$effect_definition
      tibble::tibble(
        effect_key = effect_key,
        effect_order = effect_definition$effect_order[[1]],
        panel_key = c("all", "positive", "negative"),
        panel_order = c(1L, 2L, 3L),
        panel_label = c(
          "All significant",
          effect_definition$positive_label[[1]],
          effect_definition$negative_label[[1]]
        ),
        n_cluster_gene_results = c(
          nrow(significant_results),
          sum(significant_results$direction == "positive"),
          sum(significant_results$direction == "negative")
        ),
        n_unique_genes = c(
          dplyr::n_distinct(significant_results$ensembl_gene_id),
          dplyr::n_distinct(
            significant_results$ensembl_gene_id[
              significant_results$direction == "positive"
            ]
          ),
          dplyr::n_distinct(
            significant_results$ensembl_gene_id[
              significant_results$direction == "negative"
            ]
          )
        )
      )
    }
  )
  dplyr::bind_rows(result_rows) |>
    dplyr::arrange(.data$effect_order, .data$panel_order)
}


# ==============================================================================
# 7. Scale definitions
# ==============================================================================

calculate_spatial_edger_scale_limits <- function(
  nine_panel_count_table,
  scale_mode = c("independent", "shared_by_direction")
) {
  scale_mode <- match.arg(scale_mode)
  if (identical(scale_mode, "independent")) {
    return(
      nine_panel_count_table |>
        dplyr::group_by(
          .data$effect_key,
          .data$effect_order,
          .data$panel_key,
          .data$panel_order,
          .data$scale_group
        ) |>
        dplyr::summarise(
          scale_max = max(.data$n_results, na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::mutate(scale_mode = scale_mode, .before = 1) |>
        dplyr::arrange(.data$effect_order, .data$panel_order)
    )
  }
  nine_panel_count_table |>
    dplyr::group_by(.data$scale_group) |>
    dplyr::summarise(
      scale_max = max(.data$n_results, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(scale_mode = scale_mode, .before = 1) |>
    dplyr::arrange(
      factor(.data$scale_group, levels = c("all", "positive", "negative"))
    )
}

get_spatial_edger_panel_scale_max <- function(
  scale_limits,
  effect_key,
  panel_key,
  scale_group,
  scale_mode
) {
  if (identical(scale_mode, "independent")) {
    matched <- scale_limits |>
      dplyr::filter(
        .data$effect_key == .env$effect_key,
        .data$panel_key == .env$panel_key
      )
  } else {
    matched <- scale_limits |>
      dplyr::filter(.data$scale_group == .env$scale_group)
  }
  if (nrow(matched) != 1L) {
    stop(
      "Could not uniquely resolve scale maximum for ",
      effect_key,
      " / ",
      panel_key,
      " (scale mode: ", scale_mode, "). Matched rows: ", nrow(matched),
      ".",
      call. = FALSE
    )
  }
  as.numeric(matched$scale_max[[1]])
}


# ==============================================================================
# 8. Colour mapping and one shared four-dot legend per biological effect
# ==============================================================================

map_spatial_edger_cluster_counts_to_colours <- function(
  cluster_count_table,
  count_column,
  palette_colors,
  maximum_count_override = NULL,
  zero_colour = "#D9D9D9"
) {
  count_values <- as.numeric(cluster_count_table[[count_column]])

  maximum_count <- if (is.null(maximum_count_override)) {
    max(count_values, na.rm = TRUE)
  } else {
    as.numeric(maximum_count_override)
  }

  if (
    length(maximum_count) != 1L ||
      !is.finite(maximum_count) ||
      maximum_count < 0
  ) {
    stop("Invalid colour-scale maximum.", call. = FALSE)
  }

  if (maximum_count <= 0) {
    mapped_colours <- rep(zero_colour, length(count_values))
  } else {
    colour_function <- scales::col_numeric(
      palette = palette_colors,
      domain = c(0, maximum_count),
      na.color = zero_colour
    )

    clipped_values <- pmin(
      pmax(count_values, 0),
      maximum_count
    )

    mapped_colours <- colour_function(clipped_values)
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


build_spatial_edger_effect_legend_table <- function(
  nine_panel_count_table,
  effect_key,
  cluster_names
) {
  effect_table <- nine_panel_count_table |>
    dplyr::filter(
      .data$effect_key == .env$effect_key
    )

  if (nrow(effect_table) == 0L) {
    stop(
      "No rows found for legend effect: ",
      effect_key,
      call. = FALSE
    )
  }

  effect_labels <- unique(
    as.character(effect_table$effect_label)
  )

  if (length(effect_labels) != 1L) {
    stop(
      "Could not uniquely resolve effect label for legend: ",
      effect_key,
      call. = FALSE
    )
  }

  cluster_table <- effect_table |>
    dplyr::distinct(
      .data$cluster_id,
      .data$cluster_name
    ) |>
    dplyr::arrange(
      as.integer(.data$cluster_id)
    )

  all_counts <- effect_table |>
    dplyr::filter(.data$panel_key == "all") |>
    dplyr::transmute(
      cluster_id = .data$cluster_id,
      n_all = .data$n_results
    )

  positive_counts <- effect_table |>
    dplyr::filter(.data$panel_key == "positive") |>
    dplyr::transmute(
      cluster_id = .data$cluster_id,
      n_positive = .data$n_results
    )

  negative_counts <- effect_table |>
    dplyr::filter(.data$panel_key == "negative") |>
    dplyr::transmute(
      cluster_id = .data$cluster_id,
      n_negative = .data$n_results
    )

  legend_table <- cluster_table |>
    dplyr::left_join(
      all_counts,
      by = "cluster_id"
    ) |>
    dplyr::left_join(
      positive_counts,
      by = "cluster_id"
    ) |>
    dplyr::left_join(
      negative_counts,
      by = "cluster_id"
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(
          c(
            "n_all",
            "n_positive",
            "n_negative"
          )
        ),
        ~ dplyr::coalesce(as.integer(.x), 0L)
      )
    )

  if (!is.null(cluster_names)) {
    cluster_ids <- as.character(
      legend_table$cluster_id
    )

    missing_names <- setdiff(
      cluster_ids,
      names(cluster_names)
    )

    if (length(missing_names) > 0L) {
      stop(
        "Missing cluster name(s) for legend: ",
        paste(missing_names, collapse = ", "),
        call. = FALSE
      )
    }

    legend_table$cluster_name <- unname(
      cluster_names[cluster_ids]
    )
  }

  list(
    effect_label = effect_labels[[1]],
    legend_table = legend_table
  )
}


get_spatial_edger_effect_legend_subtitle <- function(effect_key) {
  switch(
    effect_key,
    interaction = paste0(
      "Dots: anatomical cluster | all | positive interaction ",
      "(ASD effect higher in Female) | negative interaction ",
      "(ASD effect higher in Male)"
    ),
    group = "Dots: anatomical cluster | all | UP in ASD | DOWN in ASD",
    sex = "Dots: anatomical cluster | all | UP in Female | DOWN in Female",
    stop(
      "Unknown effect_key for legend subtitle: ",
      effect_key,
      call. = FALSE
    )
  )
}


create_spatial_edger_effect_four_dot_legend <- function(
  nine_panel_count_table,
  effect_key,
  scale_limits,
  scale_mode,
  green_palette_colors,
  red_palette_colors,
  blue_palette_colors,
  cluster_base_colours,
  cluster_names,
  legend_ncol = 2L,
  legend_text_size = 17,
  legend_point_size = 5.2,
  legend_row_spacing = 1.60,
  legend_title_size = 22,
  legend_subtitle_size = 22,
  legend_bottom_margin = 4,
  legend_x_offset = 0,
  legend_column_spacing = 50,
  legend_canvas_width = 100,
  legend_title_subtitle_spacing = 8,
  legend_subtitle_entries_spacing = 14,
  legend_block_alignment = "center"
) {

  legend_block_alignment <- match.arg(
    as.character(legend_block_alignment),
    choices = c(
      "left",
      "center",
      "right"
    )
  )

  numeric_positive_parameters <- c(
    legend_text_size,
    legend_point_size,
    legend_row_spacing,
    legend_title_size,
    legend_subtitle_size,
    legend_column_spacing,
    legend_canvas_width
  )

  if (
    any(!is.finite(numeric_positive_parameters)) ||
      any(numeric_positive_parameters <= 0)
  ) {
    stop(
      "Legend size/spacing parameters must be positive finite numbers.",
      call. = FALSE
    )
  }

  numeric_nonnegative_parameters <- c(
    legend_bottom_margin,
    legend_title_subtitle_spacing,
    legend_subtitle_entries_spacing
  )

  if (
    any(!is.finite(numeric_nonnegative_parameters)) ||
      any(numeric_nonnegative_parameters < 0)
  ) {
    stop(
      "Legend margin parameters must be finite non-negative numbers.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(legend_x_offset) ||
      length(legend_x_offset) != 1L ||
      !is.finite(legend_x_offset)
  ) {
    stop(
      "`legend_x_offset` must be one finite numeric value.",
      call. = FALSE
    )
  }

  legend_input <- build_spatial_edger_effect_legend_table(
    nine_panel_count_table = nine_panel_count_table,
    effect_key = effect_key,
    cluster_names = cluster_names
  )

  legend_table <- legend_input$legend_table

  cluster_ids <- as.character(
    legend_table$cluster_id
  )

  if (is.null(cluster_base_colours)) {
    stop(
      "`cluster_base_colours` must be supplied for the four-dot legend.",
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

  all_maximum <- get_spatial_edger_panel_scale_max(
    scale_limits = scale_limits,
    effect_key = effect_key,
    panel_key = "all",
    scale_group = "all",
    scale_mode = scale_mode
  )

  positive_maximum <- get_spatial_edger_panel_scale_max(
    scale_limits = scale_limits,
    effect_key = effect_key,
    panel_key = "positive",
    scale_group = "positive",
    scale_mode = scale_mode
  )

  negative_maximum <- get_spatial_edger_panel_scale_max(
    scale_limits = scale_limits,
    effect_key = effect_key,
    panel_key = "negative",
    scale_group = "negative",
    scale_mode = scale_mode
  )

  all_colour_info <- map_spatial_edger_cluster_counts_to_colours(
    cluster_count_table = legend_table,
    count_column = "n_all",
    palette_colors = green_palette_colors,
    maximum_count_override = all_maximum
  )

  positive_colour_info <- map_spatial_edger_cluster_counts_to_colours(
    cluster_count_table = legend_table,
    count_column = "n_positive",
    palette_colors = red_palette_colors,
    maximum_count_override = positive_maximum
  )

  negative_colour_info <- map_spatial_edger_cluster_counts_to_colours(
    cluster_count_table = legend_table,
    count_column = "n_negative",
    palette_colors = blue_palette_colors,
    maximum_count_override = negative_maximum
  )

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

  legend_column <- (
    (entry_index - 1L) %% legend_ncol
  ) + 1L

  legend_row <- ceiling(
    entry_index / legend_ncol
  )

  legend_labels <- paste0(
    "C",
    cluster_ids,
    " (n = ",
    legend_table$n_all,
    ", UP = ",
    legend_table$n_positive,
    ", DOWN = ",
    legend_table$n_negative,
    ") ",
    legend_table$cluster_name
  )

  # ---------------------------------------------------------------------------
  # Build each legend entry as a grid grob.
  #
  # The four dots and the text remain left-aligned inside each entry.
  # The COMPLETE set of entries is then placed as one block using a viewport
  # with just = left / center / right. This fixes the old pseudo-centering
  # problem, where only column anchor positions were centered.
  # ---------------------------------------------------------------------------

  # Dot spacing is derived from the actual point size.
  # This prevents the four legend dots from overlapping when large points
  # are requested in the runner (e.g. effect_legend_point_size = 12).
  point_radius_lines <- max(
    0.11,
    as.numeric(legend_point_size) / 28
  )

  point_diameter_lines <- 2 * point_radius_lines

  # Physical blank space between neighbouring dots, expressed in "lines".
  # The gap scales gently with point size but never becomes too small.
  point_gap_lines <- max(
    0.16,
    point_radius_lines * 0.40
  )

  point_center_spacing_lines <- point_diameter_lines +
    point_gap_lines

  dot_x_lines <- 0.30 +
    (
      0:3
    ) * point_center_spacing_lines

  # Leave a clear gap between the fourth dot and the text label.
  text_x_lines <- max(dot_x_lines) +
    point_radius_lines +
    max(
      0.225,
      point_radius_lines * 0.40
    )

  label_grobs <- lapply(
    seq_along(cluster_ids),
    function(i) {
      grid::textGrob(
        label = legend_labels[[i]],
        x = grid::unit(
          text_x_lines,
          "lines"
        ),
        y = grid::unit(
          0.5,
          "npc"
        ),
        just = c(
          "left",
          "center"
        ),
        gp = grid::gpar(
          fontfamily = "DejaVu Sans",
          fontsize = legend_text_size
        )
      )
    }
  )

  entry_grobs <- lapply(
    seq_along(cluster_ids),
    function(i) {

      dot_colours <- c(
        unname(
          cluster_base_colours[
            cluster_ids[[i]]
          ]
        ),
        unname(
          all_colour_info$colours[
            cluster_ids[[i]]
          ]
        ),
        unname(
          positive_colour_info$colours[
            cluster_ids[[i]]
          ]
        ),
        unname(
          negative_colour_info$colours[
            cluster_ids[[i]]
          ]
        )
      )

      dot_grobs <- lapply(
        seq_along(dot_x_lines),
        function(dot_index) {
          grid::circleGrob(
            x = grid::unit(
              dot_x_lines[[dot_index]],
              "lines"
            ),
            y = grid::unit(
              0.5,
              "npc"
            ),
            r = grid::unit(
              point_radius_lines,
              "lines"
            ),
            gp = grid::gpar(
              fill = dot_colours[[dot_index]],
              col = NA
            )
          )
        }
      )

      grid::gTree(
        children = do.call(
          grid::gList,
          c(
            dot_grobs,
            list(
              label_grobs[[i]]
            )
          )
        )
      )
    }
  )

  # Determine a practical width for each legend column from the longest label
  # in that column. `grobWidth()` gives a physical grid unit, so centering is
  # based on the visible legend content rather than on arbitrary x coordinates.
  column_widths <- lapply(
    seq_len(legend_ncol),
    function(column_index) {

      current_indices <- which(
        legend_column == column_index
      )

      longest_index <- current_indices[
        which.max(
          nchar(
            legend_labels[current_indices],
            type = "width"
          )
        )
      ]

      grid::unit(
        text_x_lines + 0.15,
        "lines"
      ) +
        grid::grobWidth(
          label_grobs[[longest_index]]
        )
    }
  )

  row_height <- grid::unit(
    max(
      legend_text_size * legend_row_spacing,
      legend_point_size * 2.2
    ),
    "pt"
  )

  entry_column_positions <- seq(
    1L,
    by = 2L,
    length.out = legend_ncol
  )

  entry_layout_widths <- list()

  for (column_index in seq_len(legend_ncol)) {

    entry_layout_widths[[length(entry_layout_widths) + 1L]] <-
      column_widths[[column_index]]

    if (column_index < legend_ncol) {
      entry_layout_widths[[length(entry_layout_widths) + 1L]] <-
        grid::unit(
          legend_column_spacing,
          "mm"
        )
    }
  }

  entry_layout_widths <- do.call(
    grid::unit.c,
    entry_layout_widths
  )

  entries_layout <- grid::grid.layout(
    nrow = legend_nrow,
    ncol = length(entry_layout_widths),
    widths = entry_layout_widths,
    heights = rep(
      row_height,
      legend_nrow
    ),
    respect = FALSE
  )

  entries_children <- lapply(
    seq_along(entry_grobs),
    function(i) {

      current_grob <- entry_grobs[[i]]

      current_grob$vp <- grid::viewport(
        layout.pos.row = legend_row[[i]],
        layout.pos.col = entry_column_positions[
          legend_column[[i]]
        ]
      )

      current_grob
    }
  )

  entries_block_grob <- grid::gTree(
    children = do.call(
      grid::gList,
      entries_children
    ),
    vp = grid::viewport(
      layout = entries_layout
    )
  )

  entries_block_width <- sum(
    entry_layout_widths
  )

  entries_block_height <- sum(
    rep(
      row_height,
      legend_nrow
    )
  )

  # `legend_canvas_width` is now interpreted as a percentage of the full legend
  # width. 100 = use the full available row width.
  legend_canvas_fraction <- as.numeric(
    legend_canvas_width
  ) / 100

  legend_canvas_fraction <- max(
    legend_canvas_fraction,
    0.01
  )

  block_x <- switch(
    legend_block_alignment,
    "left" = grid::unit(0, "npc") +
      grid::unit(legend_x_offset, "pt"),
    "center" = grid::unit(0.5, "npc") +
      grid::unit(legend_x_offset, "pt"),
    "right" = grid::unit(1, "npc") +
      grid::unit(legend_x_offset, "pt")
  )

  block_just <- switch(
    legend_block_alignment,
    "left" = "left",
    "center" = "center",
    "right" = "right"
  )

  entries_canvas_grob <- grid::gTree(
    children = grid::gList(
      entries_block_grob
    ),
    vp = grid::viewport(
      x = block_x,
      y = grid::unit(
        0.5,
        "npc"
      ),
      width = grid::unit(
        legend_canvas_fraction,
        "npc"
      ),
      height = entries_block_height,
      just = c(
        block_just,
        "center"
      ),
      clip = "off"
    )
  )

  # The entry block itself must keep its natural physical width. It is anchored
  # left/center/right inside the canvas viewport above.
  entries_canvas_grob$children[[1]]$vp <- grid::viewport(
    x = switch(
      legend_block_alignment,
      "left" = grid::unit(0, "npc"),
      "center" = grid::unit(0.5, "npc"),
      "right" = grid::unit(1, "npc")
    ),
    y = grid::unit(
      0.5,
      "npc"
    ),
    width = entries_block_width,
    height = entries_block_height,
    just = c(
      block_just,
      "center"
    ),
    layout = entries_layout,
    clip = "off"
  )

  title_grob <- grid::textGrob(
    label = legend_input$effect_label,
    x = grid::unit(
      0.5,
      "npc"
    ),
    just = "center",
    gp = grid::gpar(
      fontfamily = "DejaVu Sans",
      fontsize = legend_title_size,
      fontface = "bold"
    )
  )

  subtitle_grob <- grid::textGrob(
    label = get_spatial_edger_effect_legend_subtitle(
      effect_key
    ),
    x = grid::unit(
      0.5,
      "npc"
    ),
    just = "center",
    gp = grid::gpar(
      fontfamily = "DejaVu Sans",
      fontsize = legend_subtitle_size
    )
  )

  title_height <- grid::unit(
    legend_title_size * 1.25,
    "pt"
  )

  subtitle_height <- grid::unit(
    legend_subtitle_size * 1.25,
    "pt"
  )

  full_layout <- grid::grid.layout(
    nrow = 6L,
    ncol = 1L,
    heights = grid::unit.c(
      title_height,
      grid::unit(
        legend_title_subtitle_spacing,
        "pt"
      ),
      subtitle_height,
      grid::unit(
        legend_subtitle_entries_spacing,
        "pt"
      ),
      entries_block_height,
      grid::unit(
        legend_bottom_margin,
        "pt"
      )
    )
  )

  title_grob$vp <- grid::viewport(
    layout.pos.row = 1L
  )

  subtitle_grob$vp <- grid::viewport(
    layout.pos.row = 3L
  )

  entries_canvas_grob$vp <- grid::viewport(
    layout.pos.row = 5L,
    clip = "off"
  )

  complete_legend_grob <- grid::gTree(
    children = grid::gList(
      title_grob,
      subtitle_grob,
      entries_canvas_grob
    ),
    vp = grid::viewport(
      layout = full_layout,
      clip = "off"
    )
  )

  patchwork::wrap_elements(
    full = complete_legend_grob,
    clip = FALSE
  )
}

# ==============================================================================
# 9. Create one spatial panel without an internal legend
# ==============================================================================

create_significant_gene_count_spatial_panel <- function(
  section_context,
  panel_count_table,
  panel_title,
  panel_subtitle,
  palette_colors,
  maximum_count,
  show_histology_image = FALSE,
  point_size_no_image = 3.00,
  point_size_with_image = 2.10,
  panel_border_linewidth = 1.2,
  panel_title_size = 21,
  panel_subtitle_size = 16
) {
  colour_info <- map_spatial_edger_cluster_counts_to_colours(
    cluster_count_table = panel_count_table,
    count_column = "n_results",
    palette_colors = palette_colors,
    maximum_count_override = maximum_count
  )

  plot_data <- section_context$plot_data

  cluster_levels <- as.character(
    panel_count_table$cluster_id
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

  spatial_plot +
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
      limits = rev(
        section_context$frame_limits$y_limits
      ),
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
        size = panel_title_size,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(
          t = 0,
          b = 2
        )
      ),
      plot.subtitle = ggplot2::element_text(
        size = panel_subtitle_size,
        face = "plain",
        hjust = 0.5,
        lineheight = 1.02,
        margin = ggplot2::margin(
          t = 0,
          b = 5
        )
      ),
      plot.margin = ggplot2::margin(
        t = 0,
        r = 8,
        b = 8,
        l = 8
      )
    )
}


# ==============================================================================
# 10. Create one 3 x 3 figure with three shared legends
# ==============================================================================

get_spatial_edger_palette_for_panel <- function(
  panel_key,
  green_palette_colors,
  red_palette_colors,
  blue_palette_colors
) {
  switch(
    panel_key,
    all = green_palette_colors,
    positive = red_palette_colors,
    negative = blue_palette_colors,
    stop(
      "Unknown panel_key: ",
      panel_key,
      call. = FALSE
    )
  )
}


create_spatial_edger_nine_panel_figure <- function(
  section_context,
  nine_panel_count_table,
  panel_summary_table,
  scale_limits,
  scale_mode,
  clustering_parameters_label,
  fdr_threshold,
  abs_log2fc_threshold,
  min_max_mean_percent,
  max_max_mean_percent,
  green_palette_colors,
  red_palette_colors,
  blue_palette_colors,
  cluster_base_colours,
  cluster_names,
  show_histology_image = FALSE,
  point_size_no_image = 3.00,
  point_size_with_image = 2.10,
  effect_legend_ncol = 3L,
  effect_legend_text_size = 13.5,
  effect_legend_point_size = 5.2,
  effect_legend_row_spacing = 2.10,
  effect_legend_title_size = 15,
  effect_legend_subtitle_size = 11.8,
  effect_legend_bottom_margin = 4,
  effect_legend_relative_height = 0.43,
  panel_title_size = 21,
  panel_subtitle_size = 16,
  main_title_size = 32,
  main_subtitle_size = 17,
  main_title_subtitle_spacing = 5,
  main_subtitle_bottom_spacing = 8,
  effect_legend_x_offset = 0,
  effect_legend_column_spacing = 50,
  effect_legend_canvas_width = 100,
  effect_legend_title_subtitle_spacing = 8,
  effect_legend_subtitle_entries_spacing = 12,
  effect_legend_block_alignment = "center"
) {
  effect_keys <- nine_panel_count_table |>
    dplyr::distinct(
      .data$effect_key,
      .data$effect_order
    ) |>
    dplyr::arrange(
      .data$effect_order
    ) |>
    dplyr::pull(
      "effect_key"
    )

  if (!identical(
    as.character(effect_keys),
    c("interaction", "group", "sex")
  )) {
    stop(
      "Expected effect order: interaction, group, sex. Observed: ",
      paste(effect_keys, collapse = ", "),
      call. = FALSE
    )
  }

  row_blocks <- lapply(
    effect_keys,
    function(effect_key_current) {
      effect_panel_info <- nine_panel_count_table |>
        dplyr::filter(
          .data$effect_key == .env$effect_key_current
        ) |>
        dplyr::distinct(
          .data$effect_key,
          .data$effect_label,
          .data$panel_key,
          .data$panel_order,
          .data$scale_group,
          .data$panel_label,
          .data$panel_interpretation
        ) |>
        dplyr::arrange(
          .data$panel_order
        )

      if (!identical(
        as.character(effect_panel_info$panel_key),
        c("all", "positive", "negative")
      )) {
        stop(
          "Expected panel order all/positive/negative for effect: ",
          effect_key_current,
          call. = FALSE
        )
      }

      effect_panel_plots <- lapply(
        seq_len(nrow(effect_panel_info)),
        function(panel_index) {
          panel_info <- effect_panel_info[
            panel_index,
            ,
            drop = FALSE
          ]

          panel_count_table <- nine_panel_count_table |>
            dplyr::filter(
              .data$effect_key == .env$effect_key_current,
              .data$panel_key == panel_info$panel_key[[1]]
            ) |>
            dplyr::arrange(
              as.integer(.data$cluster_id)
            )

          panel_summary <- panel_summary_table |>
            dplyr::filter(
              .data$effect_key == .env$effect_key_current,
              .data$panel_key == panel_info$panel_key[[1]]
            )

          if (nrow(panel_summary) != 1L) {
            stop(
              "Could not uniquely resolve panel summary for ",
              effect_key_current,
              " / ",
              panel_info$panel_key[[1]],
              ".",
              call. = FALSE
            )
          }

          maximum_count <- get_spatial_edger_panel_scale_max(
            scale_limits = scale_limits,
            effect_key = effect_key_current,
            panel_key = panel_info$panel_key[[1]],
            scale_group = panel_info$scale_group[[1]],
            scale_mode = scale_mode
          )

          palette_colors <- get_spatial_edger_palette_for_panel(
            panel_key = panel_info$panel_key[[1]],
            green_palette_colors = green_palette_colors,
            red_palette_colors = red_palette_colors,
            blue_palette_colors = blue_palette_colors
          )

          panel_title <- paste0(
            panel_info$effect_label[[1]],
            " — ",
            panel_info$panel_label[[1]]
          )

          panel_subtitle <- paste0(
            format_spatial_edger_integer(
              panel_summary$n_cluster_gene_results[[1]]
            ),
            " cluster-gene results | ",
            format_spatial_edger_integer(
              panel_summary$n_unique_genes[[1]]
            ),
            " unique genes | scale max = ",
            format_spatial_edger_integer(
              maximum_count
            )
          )

          create_significant_gene_count_spatial_panel(
            section_context = section_context,
            panel_count_table = panel_count_table,
            panel_title = panel_title,
            panel_subtitle = panel_subtitle,
            palette_colors = palette_colors,
            maximum_count = maximum_count,
            show_histology_image = show_histology_image,
            point_size_no_image = point_size_no_image,
            point_size_with_image = point_size_with_image,
            panel_title_size = panel_title_size,
            panel_subtitle_size = panel_subtitle_size
          )
        }
      )

      shared_effect_legend <- create_spatial_edger_effect_four_dot_legend(
        nine_panel_count_table = nine_panel_count_table,
        effect_key = effect_key_current,
        scale_limits = scale_limits,
        scale_mode = scale_mode,
        green_palette_colors = green_palette_colors,
        red_palette_colors = red_palette_colors,
        blue_palette_colors = blue_palette_colors,
        cluster_base_colours = cluster_base_colours,
        cluster_names = cluster_names,
        legend_ncol = effect_legend_ncol,
        legend_text_size = effect_legend_text_size,
        legend_point_size = effect_legend_point_size,
        legend_row_spacing = effect_legend_row_spacing,
        legend_title_size = effect_legend_title_size,
        legend_subtitle_size = effect_legend_subtitle_size,
        legend_bottom_margin = effect_legend_bottom_margin,
        legend_x_offset = effect_legend_x_offset,
        legend_column_spacing = effect_legend_column_spacing,
        legend_canvas_width = effect_legend_canvas_width,
        legend_title_subtitle_spacing = effect_legend_title_subtitle_spacing,
        legend_subtitle_entries_spacing = effect_legend_subtitle_entries_spacing,
        legend_block_alignment = effect_legend_block_alignment
      )

      three_panel_row <- patchwork::wrap_plots(
        plots = effect_panel_plots,
        nrow = 1L,
        ncol = 3L
      )

      row_patchwork <- (
        shared_effect_legend /
          three_panel_row
      ) +
        patchwork::plot_layout(
          ncol = 1L,
          heights = c(
            effect_legend_relative_height,
            1
          )
        )

      row_grob <- patchwork::patchworkGrob(
        row_patchwork
      )

      patchwork::wrap_elements(
        full = row_grob,
        clip = FALSE
      )
    }
  )

  scale_mode_label <- if (
    identical(scale_mode, "independent")
  ) {
    "Scale mode: independent scale for each of the 9 panels"
  } else {
    paste0(
      "Scale mode: shared across effects separately for All, Positive and ",
      "Negative panels"
    )
  }

  plot_title <- paste0(
    "Significant genes per spatial cluster across edgeR effects"
  )

  plot_subtitle <- paste0(
    "Selection: ",
    format_spatial_edger_filter_label(
      fdr_threshold = fdr_threshold,
      abs_log2fc_threshold = abs_log2fc_threshold,
      min_max_mean_percent = min_max_mean_percent,
      max_max_mean_percent = max_max_mean_percent
    ),
    "\nSample: ",
    section_context$sample_ID,
    " | tissue spots: ",
    format_spatial_edger_integer(
      section_context$n_spots
    ),
    " | ",
    scale_mode_label,
    "\n",
    clustering_parameters_label,
    "\nInteraction sign: positive = ASD effect higher in Female; negative = ASD effect higher in Male"
  )

  three_row_grid <- patchwork::wrap_plots(
    plots = row_blocks,
    nrow = 3L,
    ncol = 1L
  )

  three_row_grid +
    patchwork::plot_annotation(
      title = plot_title,
      subtitle = plot_subtitle,
      theme = ggplot2::theme(
        text = ggplot2::element_text(
          family = "DejaVu Sans"
        ),
        plot.title = ggplot2::element_text(
          size = main_title_size,
          face = "bold",
          hjust = 0.5,
          margin = ggplot2::margin(
            b = main_title_subtitle_spacing
          )
        ),
        plot.subtitle = ggplot2::element_text(
          size = main_subtitle_size,
          hjust = 0.5,
          lineheight = 1.08,
          margin = ggplot2::margin(
            b = main_subtitle_bottom_spacing
          )
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


# ==============================================================================
# 11. Output file definitions
# ==============================================================================

get_all_effect_significant_gene_count_output_files <- function(
  output_root_directory,
  sample_id_to_plot,
  fdr_threshold,
  abs_log2fc_threshold,
  min_max_mean_percent = NULL,
  max_max_mean_percent = NULL
) {
  parameter_label <- build_spatial_edger_parameter_label(
    fdr_threshold = fdr_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold,
    min_max_mean_percent = min_max_mean_percent,
    max_max_mean_percent = max_max_mean_percent
  )
  parameter_output_directory <- file.path(
    output_root_directory,
    parameter_label
  )
  safe_sample_id <- sanitize_spatial_edger_filename_component(
    sample_id_to_plot
  )
  file_prefix <- paste0(
    "allEffects_",
    parameter_label,
    "_sample_",
    safe_sample_id
  )
  list(
    parameter_label = parameter_label,
    parameter_output_directory = parameter_output_directory,
    independent_scale_pdf = file.path(
      parameter_output_directory,
      paste0(
        "01_",
        file_prefix,
        "_significantGeneCounts_3x3_independentScales.pdf"
      )
    ),
    shared_scale_pdf = file.path(
      parameter_output_directory,
      paste0(
        "02_",
        file_prefix,
        "_significantGeneCounts_3x3_sharedAllPositiveNegativeScales.pdf"
      )
    ),
    counts_tsv = file.path(
      parameter_output_directory,
      paste0(
        "03_",
        file_prefix,
        "_significantGeneCountsPerCluster.tsv"
      )
    ),
    selected_results_tsv = file.path(
      parameter_output_directory,
      paste0(
        "04_",
        file_prefix,
        "_selectedClusterGeneResults.tsv"
      )
    ),
    scale_limits_tsv = file.path(
      parameter_output_directory,
      paste0(
        "05_",
        file_prefix,
        "_scaleLimits.tsv"
      )
    ),
    summary_tsv = file.path(
      parameter_output_directory,
      paste0(
        "06_",
        file_prefix,
        "_analysisSummary.tsv"
      )
    )
  )
}


# ==============================================================================
# 12. Main workflow
# ==============================================================================

run_all_effects_significant_gene_count_spatial_workflow <- function(
  input_seurat_rdata_file,
  edger_results_rdata_file,
  output_root_directory,
  sample_id_to_plot,
  cluster_column,
  clustering_parameters_label,
  sample_id_column = "sample_ID",
  seurat_object_name = NULL,
  image_scale = "lowres",
  interaction_test_id = "Interaction",
  group_test_id = "Overall_Group_ASD_vs_Neurotypical",
  sex_test_id = "Overall_Sex_Female_vs_Male",
  fdr_threshold = 0.05,
  abs_log2fc_threshold = 0.5,
  min_max_mean_percent = NULL,
  max_max_mean_percent = NULL,
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
  cluster_base_colours,
  cluster_names,
  show_histology_image = FALSE,
  point_size_no_image = 3.00,
  point_size_with_image = 2.10,
  effect_legend_ncol = 3L,
  effect_legend_text_size = 13.5,
  effect_legend_point_size = 5.2,
  effect_legend_row_spacing = 2.10,
  effect_legend_title_size = 15,
  effect_legend_subtitle_size = 11.8,
  effect_legend_bottom_margin = 4,
  effect_legend_relative_height = 0.43,
  effect_legend_x_offset = 0,
  effect_legend_column_spacing = 50,
  effect_legend_canvas_width = 100,
  effect_legend_title_subtitle_spacing = 8,
  effect_legend_subtitle_entries_spacing = 12,
  effect_legend_block_alignment = "center",
  panel_title_size = 21,
  panel_subtitle_size = 16,
  main_title_size = 32,
  main_subtitle_size = 17,
  main_title_subtitle_spacing = 5,
  main_subtitle_bottom_spacing = 8,
  pdf_width = 36,
  pdf_height = 34
) {
  if (
    !is.numeric(fdr_threshold) ||
      length(fdr_threshold) != 1L ||
      !is.finite(fdr_threshold) ||
      fdr_threshold <= 0 ||
      fdr_threshold >= 1
  ) {
    stop("`fdr_threshold` must be one finite value between 0 and 1.", call. = FALSE)
  }
  if (
    !is.numeric(effect_legend_column_spacing) ||
      length(effect_legend_column_spacing) != 1L ||
      !is.finite(effect_legend_column_spacing) ||
      effect_legend_column_spacing <= 0
  ) {
    stop(
      "`effect_legend_column_spacing` must be one positive finite numeric value.",
      call. = FALSE
    )
  }

  main_heading_spacing_values <- c(
    main_title_subtitle_spacing,
    main_subtitle_bottom_spacing
  )

  if (
    any(!is.finite(main_heading_spacing_values)) ||
      any(main_heading_spacing_values < 0)
  ) {
    stop(
      "Main title/subtitle spacing values must be finite non-negative numbers.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(effect_legend_canvas_width) ||
      length(effect_legend_canvas_width) != 1L ||
      !is.finite(effect_legend_canvas_width) ||
      effect_legend_canvas_width <= 0
  ) {
    stop(
      "`effect_legend_canvas_width` must be one positive finite numeric value.",
      call. = FALSE
    )
  }

  effect_legend_vertical_spacing <- c(
    effect_legend_title_subtitle_spacing,
    effect_legend_subtitle_entries_spacing
  )

  if (
    any(!is.finite(effect_legend_vertical_spacing)) ||
      any(effect_legend_vertical_spacing < 0)
  ) {
    stop(
      "Legend vertical-spacing settings must be finite non-negative values.",
      call. = FALSE
    )
  }

  effect_legend_block_alignment <- match.arg(
    as.character(effect_legend_block_alignment),
    choices = c(
      "left",
      "center",
      "right"
    )
  )
  if (
    !is.numeric(abs_log2fc_threshold) ||
      length(abs_log2fc_threshold) != 1L ||
      !is.finite(abs_log2fc_threshold) ||
      abs_log2fc_threshold < 0
  ) {
    stop("`abs_log2fc_threshold` must be one non-negative finite value.", call. = FALSE)
  }
  if (!is.null(min_max_mean_percent)) {
    if (
      !is.numeric(min_max_mean_percent) ||
        length(min_max_mean_percent) != 1L ||
        !is.finite(min_max_mean_percent) ||
        min_max_mean_percent < 0 ||
        min_max_mean_percent > 100
    ) {
      stop(
        "`min_max_mean_percent` must be NULL or one finite value in [0, 100].",
        call. = FALSE
      )
    }
  }

  if (!is.null(max_max_mean_percent)) {
    if (
      !is.numeric(max_max_mean_percent) ||
        length(max_max_mean_percent) != 1L ||
        !is.finite(max_max_mean_percent) ||
        max_max_mean_percent < 0 ||
        max_max_mean_percent > 100
    ) {
      stop(
        "`max_max_mean_percent` must be NULL or one finite value in [0, 100].",
        call. = FALSE
      )
    }
  }

  if (
    !is.null(min_max_mean_percent) &&
      !is.null(max_max_mean_percent) &&
      min_max_mean_percent >= max_max_mean_percent
  ) {
    stop(
      "`min_max_mean_percent` must be smaller than `max_max_mean_percent`.",
      call. = FALSE
    )
  }
  if (is.null(cluster_base_colours) || is.null(cluster_names)) {
    stop(
      "Both `cluster_base_colours` and `cluster_names` must be supplied.",
      call. = FALSE
    )
  }
  dir.create(
    output_root_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  output_files <- get_all_effect_significant_gene_count_output_files(
    output_root_directory = output_root_directory,
    sample_id_to_plot = sample_id_to_plot,
    fdr_threshold = fdr_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold,
    min_max_mean_percent = min_max_mean_percent,
    max_max_mean_percent = max_max_mean_percent
  )
  dir.create(
    output_files$parameter_output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  message("============================================================")
  message("Spatial significant-gene-count overview: Interaction + Group + Sex")
  message("============================================================")
  message("Parameters: ", output_files$parameter_label)
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
  edger_loaded <- load_all_edger_effects_for_significant_gene_count_plot(
    edger_results_rdata_file = edger_results_rdata_file,
    interaction_test_id = interaction_test_id,
    group_test_id = group_test_id,
    sex_test_id = sex_test_id
  )
  effect_definitions <- get_spatial_edger_effect_definitions(
    interaction_test_id = interaction_test_id,
    group_test_id = group_test_id,
    sex_test_id = sex_test_id
  )

  message("\nResolved edgeR tests and available mean-percent columns:")
  for (diagnostic_effect_key in c("interaction", "group", "sex")) {
    diagnostic_definition <- effect_definitions |>
      dplyr::filter(.data$effect_key == .env$diagnostic_effect_key)
    diagnostic_percent_columns <- grep(
      "^mean_percent_positive_spots_",
      colnames(edger_loaded$results[[diagnostic_effect_key]]),
      value = TRUE
    )
    message(
      "  ", diagnostic_effect_key,
      " | test_id = ", diagnostic_definition$test_id[[1]],
      " | percent columns = ",
      if (length(diagnostic_percent_columns) == 0L) {
        "<none>"
      } else {
        paste(diagnostic_percent_columns, collapse = ", ")
      }
    )
  }

  message(
    "Percent-positive filter: ",
    if (
      is.null(min_max_mean_percent) &&
        is.null(max_max_mean_percent)
    ) {
      "DISABLED (percent-positive columns will not be inspected)"
    } else {
      format_spatial_edger_filter_label(
        fdr_threshold = fdr_threshold,
        abs_log2fc_threshold = abs_log2fc_threshold,
        min_max_mean_percent = min_max_mean_percent,
        max_max_mean_percent = max_max_mean_percent
      )
    }
  )

  effect_summaries <- prepare_all_spatial_edger_effect_summaries(
    edger_results = edger_loaded$results,
    effect_definitions = effect_definitions,
    cluster_levels = section_context$cluster_levels,
    fdr_threshold = fdr_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold,
    min_max_mean_percent = min_max_mean_percent,
    max_max_mean_percent = max_max_mean_percent,
    cluster_names = cluster_names
  )
  summary_table <- dplyr::bind_rows(
    lapply(effect_summaries, function(x) x$summary_table)
  ) |>
    dplyr::arrange(.data$effect_order)
  expected_effect_keys <- c("interaction", "group", "sex")
  expected_effect_labels <- c(
    "Group x Sex interaction",
    "FMT donor group",
    "Sex"
  )
  if (
    nrow(summary_table) != 3L ||
    !identical(as.character(summary_table$effect_key), expected_effect_keys) ||
    !identical(as.character(summary_table$effect_label), expected_effect_labels)
  ) {
    stop(
      "Internal effect-mapping error. Expected interaction/group/sex exactly once.\n",
      "Observed effect_key: ", paste(summary_table$effect_key, collapse = ", "), "\n",
      "Observed effect_label: ", paste(summary_table$effect_label, collapse = ", "),
      call. = FALSE
    )
  }
  message("\nSelected results:")
  print(
    summary_table |>
      dplyr::select(
        "effect_label",
        "test_id",
        "n_all",
        "n_positive",
        "n_negative",
        "unique_genes_all"
      ),
    n = Inf
  )
  nine_panel_count_table <- build_spatial_edger_nine_panel_count_table(
    effect_summaries
  )
  panel_summary_table <- build_spatial_edger_panel_summary_table(
    nine_panel_count_table = nine_panel_count_table,
    effect_summaries = effect_summaries
  )
  independent_scale_limits <- calculate_spatial_edger_scale_limits(
    nine_panel_count_table = nine_panel_count_table,
    scale_mode = "independent"
  )
  shared_scale_limits <- calculate_spatial_edger_scale_limits(
    nine_panel_count_table = nine_panel_count_table,
    scale_mode = "shared_by_direction"
  )
  all_scale_limits <- dplyr::bind_rows(
    independent_scale_limits,
    shared_scale_limits
  )
  selected_results_table <- dplyr::bind_rows(
    lapply(effect_summaries, function(x) x$significant_results)
  ) |>
    dplyr::arrange(
      .data$effect_order,
      as.integer(.data$cluster_id),
      .data$FDR,
      dplyr::desc(abs(.data$logFC))
    )
  independent_figure <- create_spatial_edger_nine_panel_figure(
    section_context = section_context,
    nine_panel_count_table = nine_panel_count_table,
    panel_summary_table = panel_summary_table,
    scale_limits = independent_scale_limits,
    scale_mode = "independent",
    clustering_parameters_label = clustering_parameters_label,
    fdr_threshold = fdr_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold,
    min_max_mean_percent = min_max_mean_percent,
    max_max_mean_percent = max_max_mean_percent,
    green_palette_colors = green_palette_colors,
    red_palette_colors = red_palette_colors,
    blue_palette_colors = blue_palette_colors,
    cluster_base_colours = cluster_base_colours,
    cluster_names = cluster_names,
    show_histology_image = show_histology_image,
    point_size_no_image = point_size_no_image,
    point_size_with_image = point_size_with_image,
    effect_legend_ncol = effect_legend_ncol,
    effect_legend_text_size = effect_legend_text_size,
    effect_legend_point_size = effect_legend_point_size,
    effect_legend_row_spacing = effect_legend_row_spacing,
    effect_legend_title_size = effect_legend_title_size,
    effect_legend_subtitle_size = effect_legend_subtitle_size,
    effect_legend_bottom_margin = effect_legend_bottom_margin,
    effect_legend_relative_height = effect_legend_relative_height,
    effect_legend_x_offset = effect_legend_x_offset,
    effect_legend_column_spacing = effect_legend_column_spacing,
    effect_legend_canvas_width = effect_legend_canvas_width,
    effect_legend_title_subtitle_spacing = effect_legend_title_subtitle_spacing,
    effect_legend_subtitle_entries_spacing = effect_legend_subtitle_entries_spacing,
    effect_legend_block_alignment = effect_legend_block_alignment,
    panel_title_size = panel_title_size,
    panel_subtitle_size = panel_subtitle_size,
    main_title_size = main_title_size,
    main_subtitle_size = main_subtitle_size,
    main_title_subtitle_spacing = main_title_subtitle_spacing,
    main_subtitle_bottom_spacing = main_subtitle_bottom_spacing
  )
  shared_figure <- create_spatial_edger_nine_panel_figure(
    section_context = section_context,
    nine_panel_count_table = nine_panel_count_table,
    panel_summary_table = panel_summary_table,
    scale_limits = shared_scale_limits,
    scale_mode = "shared_by_direction",
    clustering_parameters_label = clustering_parameters_label,
    fdr_threshold = fdr_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold,
    min_max_mean_percent = min_max_mean_percent,
    max_max_mean_percent = max_max_mean_percent,
    green_palette_colors = green_palette_colors,
    red_palette_colors = red_palette_colors,
    blue_palette_colors = blue_palette_colors,
    cluster_base_colours = cluster_base_colours,
    cluster_names = cluster_names,
    show_histology_image = show_histology_image,
    point_size_no_image = point_size_no_image,
    point_size_with_image = point_size_with_image,
    effect_legend_ncol = effect_legend_ncol,
    effect_legend_text_size = effect_legend_text_size,
    effect_legend_point_size = effect_legend_point_size,
    effect_legend_row_spacing = effect_legend_row_spacing,
    effect_legend_title_size = effect_legend_title_size,
    effect_legend_subtitle_size = effect_legend_subtitle_size,
    effect_legend_bottom_margin = effect_legend_bottom_margin,
    effect_legend_relative_height = effect_legend_relative_height,
    effect_legend_x_offset = effect_legend_x_offset,
    effect_legend_column_spacing = effect_legend_column_spacing,
    effect_legend_canvas_width = effect_legend_canvas_width,
    effect_legend_title_subtitle_spacing = effect_legend_title_subtitle_spacing,
    effect_legend_subtitle_entries_spacing = effect_legend_subtitle_entries_spacing,
    effect_legend_block_alignment = effect_legend_block_alignment,
    panel_title_size = panel_title_size,
    panel_subtitle_size = panel_subtitle_size,
    main_title_size = main_title_size,
    main_subtitle_size = main_subtitle_size,
    main_title_subtitle_spacing = main_title_subtitle_spacing,
    main_subtitle_bottom_spacing = main_subtitle_bottom_spacing
  )
  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }
  message("\nWriting independent-scale PDF:")
  message(output_files$independent_scale_pdf)
  ggplot2::ggsave(
    filename = output_files$independent_scale_pdf,
    plot = independent_figure,
    device = pdf_device,
    width = pdf_width,
    height = pdf_height,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )
  message("\nWriting shared-scale PDF:")
  message(output_files$shared_scale_pdf)
  ggplot2::ggsave(
    filename = output_files$shared_scale_pdf,
    plot = shared_figure,
    device = pdf_device,
    width = pdf_width,
    height = pdf_height,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )
  write_spatial_edger_tsv(
    nine_panel_count_table,
    output_files$counts_tsv
  )
  write_spatial_edger_tsv(
    selected_results_table,
    output_files$selected_results_tsv
  )
  write_spatial_edger_tsv(
    all_scale_limits,
    output_files$scale_limits_tsv
  )
  analysis_summary <- summary_table |>
    dplyr::mutate(
      sample_ID = section_context$sample_ID,
      n_spots = section_context$n_spots,
      parameter_label = output_files$parameter_label,
      clustering_parameters = clustering_parameters_label,
      output_independent_scale_pdf = output_files$independent_scale_pdf,
      output_shared_scale_pdf = output_files$shared_scale_pdf,
      .before = 1
    )
  write_spatial_edger_tsv(
    analysis_summary,
    output_files$summary_tsv
  )
  expected_output_files <- c(
    output_files$independent_scale_pdf,
    output_files$shared_scale_pdf,
    output_files$counts_tsv,
    output_files$selected_results_tsv,
    output_files$scale_limits_tsv,
    output_files$summary_tsv
  )
  missing_output_files <- expected_output_files[
    !file.exists(expected_output_files)
  ]
  if (length(missing_output_files) > 0L) {
    stop(
      "Missing output file(s):\n",
      paste(missing_output_files, collapse = "\n"),
      call. = FALSE
    )
  }
  empty_output_files <- expected_output_files[
    file.info(expected_output_files)$size == 0
  ]
  if (length(empty_output_files) > 0L) {
    stop(
      "Empty output file(s):\n",
      paste(empty_output_files, collapse = "\n"),
      call. = FALSE
    )
  }
  message("\n============================================================")
  message("Completed spatial significant-gene-count overview.")
  message("Output directory:")
  message(
    normalizePath(
      output_files$parameter_output_directory,
      mustWork = TRUE
    )
  )
  message("============================================================")
  invisible(
    list(
      section_context = section_context,
      effect_summaries = effect_summaries,
      summary_table = summary_table,
      panel_summary_table = panel_summary_table,
      nine_panel_count_table = nine_panel_count_table,
      scale_limits = all_scale_limits,
      selected_results = selected_results_table,
      output_files = output_files,
      output_directory = output_files$parameter_output_directory
    )
  )
}

# ==============================================================================
# End
# ==============================================================================
