#!/usr/bin/env Rscript
# ==============================================================================
# 02_volcanoPlotPerCluster_allEdgeREffects_grid_twoVersions.R
#
# Per-cluster volcano plots shown as one faceted grid per biological effect.
#
# Effects:
#   1. Group x Sex interaction
#   2. Main FMT donor-group effect
#   3. Main Sex effect
#
# Outputs:
#   exactly three PDFs, one PDF per effect
#   each PDF contains two pages:
#     page 1 = 4 x 4 grid without gene labels
#     page 2 = 4 x 4 grid with labels for top significant genes
#
# Volcano axes:
#   x = log2FC
#   y = -log10(FDR)
#
# Threshold lines:
#   horizontal: FDR < 0.05
#   vertical:   |log2FC| >= 0.5
#   dashed black
#
# Point colours:
#   significant UP   = red
#   significant DOWN = blue
#   all other genes  = light grey
#
# Important requested visual settings:
#   - no point outlines
#   - lighter grey than before for non-significant points
#
# Current Max mean % filter:
#   0 to 100%
# ==============================================================================

options(stringsAsFactors = FALSE)

# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tibble",
  "ggplot2",
  "ggrepel"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})

# ==============================================================================
# 2. Analysis configuration
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"

clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

fdr_threshold <- 0.05
abs_log2fc_threshold <- 0.5

min_max_mean_percent <- 0
max_max_mean_percent <- 100

top_n_labels_up <- 10L
top_n_labels_down <- 10L

facet_ncol <- 4

colour_up <- "#B2182B"
colour_down <- "#2166AC"
colour_not_significant <- "grey85"

point_size <- 0.95
point_alpha <- 0.80
point_shape <- 16
point_stroke <- 0

threshold_linewidth <- 0.40
threshold_linetype <- "dashed"

label_text_size <- 2.45
label_box_padding <- 0.22
label_point_padding <- 0.12
label_segment_colour <- "grey45"
label_segment_linewidth <- 0.22
label_seed <- 1234

pdf_width <- 16
pdf_height <- 13

# ==============================================================================
# 3. Cluster labels
# ==============================================================================

custom_cluster_labels <- c(
  "1"  = "posterior & sensory relay thalamic nuclei",
  "2"  = "isocortex, layers 4 & 5",
  "3"  = "cortical layers 1 & hippocampal neuropil",
  "4"  = "Fiber tracts",
  "5"  = "cortical subplate & deep olfactory areas",
  "6"  = "hypothalamus",
  "7"  = "isocortex, layer 6",
  "8"  = "isocortex, layer 2/3",
  "9"  = "reticular, ventral geniculate & habenular region",
  "10" = "striatum-like amygdala nuclei",
  "11" = "medial thalamic nuclei",
  "12" = "caudoputamen",
  "13" = "hippocampal CA fields, pyramidal layer",
  "14" = "meninges & vasculature",
  "15" = "ventricles",
  "16" = "dentate gyrus"
)

# ==============================================================================
# 4. edgeR effect definitions
# ==============================================================================

effect_definitions <- tibble::tribble(
  ~effect_key, ~effect_order, ~test_id, ~file_stub, ~effect_title,
  "interaction", 1L, "Interaction", "interaction", "Group x Sex interaction",
  "group", 2L, "Overall_Group_ASD_vs_Neurotypical", "effectDonorGroup", "FMT donor-group main effect",
  "sex", 3L, "Overall_Sex_Female_vs_Male", "effectSex", "Sex main effect"
)

# ==============================================================================
# 5. Input and output paths
# ==============================================================================

statistics_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)

analysis_prefix <- paste0(
  dataset_name,
  "_",
  clustering_name,
  "_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction"
)

edgeR_results_rdata_file <- file.path(
  statistics_dir,
  "03_edgeRResults",
  paste0(
    analysis_prefix,
    "_edgeRResults.RData"
  )
)

volcano_output_root <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "06_volcanoPlot"
)

if (!file.exists(edgeR_results_rdata_file)) {
  stop(
    "edgeR results RData file does not exist:\n",
    edgeR_results_rdata_file,
    call. = FALSE
  )
}

# ==============================================================================
# 6. Helpers
# ==============================================================================

sort_cluster_ids <- function(cluster_ids) {
  cluster_ids <- unique(as.character(cluster_ids))
  numeric_ids <- suppressWarnings(as.numeric(cluster_ids))
  if (!anyNA(numeric_ids)) {
    return(cluster_ids[order(numeric_ids)])
  }
  sort(cluster_ids)
}

format_number_for_filename <- function(x) {
  output <- format(
    as.numeric(x),
    scientific = FALSE,
    trim = TRUE,
    digits = 8
  )
  gsub("\\.", "p", output)
}

build_parameter_label <- function(
  fdr_threshold,
  abs_log2fc_threshold,
  min_max_mean_percent = NULL,
  max_max_mean_percent = NULL
) {
  fdr_code <- gsub(
    "\\.",
    "",
    formatC(as.numeric(fdr_threshold), format = "f", digits = 2)
  )

  logfc_code <- gsub(
    "\\.",
    "",
    formatC(as.numeric(abs_log2fc_threshold), format = "f", digits = 1)
  )

  percent_code <- if (
    !is.null(min_max_mean_percent) &&
      !is.null(max_max_mean_percent)
  ) {
    paste0(
      format_number_for_filename(min_max_mean_percent),
      "to",
      format_number_for_filename(max_max_mean_percent)
    )
  } else {
    "All"
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

format_percent_filter_label <- function(
  min_max_mean_percent = NULL,
  max_max_mean_percent = NULL
) {
  if (
    is.null(min_max_mean_percent) &&
      is.null(max_max_mean_percent)
  ) {
    return("Max mean % filter: none")
  }

  upper_operator <- if (
    !is.null(max_max_mean_percent) &&
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
}

make_cluster_facet_label <- function(cluster_id) {
  cluster_id <- as.character(cluster_id)
  anatomical_label <- if (cluster_id %in% names(custom_cluster_labels)) {
    unname(custom_cluster_labels[[cluster_id]])
  } else {
    "anatomical label unavailable"
  }

  paste0(
    "C",
    cluster_id,
    "\n",
    anatomical_label
  )
}

# ==============================================================================
# 7. Resolve percent-positive columns
# ==============================================================================

resolve_percent_columns <- function(result_table, effect_key) {

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

  available_columns <- grep(
    "^mean_percent_positive_spots_",
    colnames(result_table),
    value = TRUE
  )

  stop(
    "Could not resolve exact percent-positive columns for effect '",
    effect_key,
    "'.\nExpected: ",
    paste(exact_columns, collapse = ", "),
    "\nAvailable: ",
    if (length(available_columns) == 0L) {
      "<none>"
    } else {
      paste(available_columns, collapse = ", ")
    },
    call. = FALSE
  )
}

calculate_max_mean_percent <- function(result_table, percent_columns) {

  percent_matrix <- as.matrix(
    result_table[, percent_columns, drop = FALSE]
  )

  storage.mode(percent_matrix) <- "numeric"

  apply(
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
}

# ==============================================================================
# 8. Load edgeR objects
# ==============================================================================

message("Loading edgeR results:")
message(edgeR_results_rdata_file)

result_environment <- new.env(parent = globalenv())

loaded_objects <- load(
  edgeR_results_rdata_file,
  envir = result_environment
)

required_objects <- c(
  "edgeR_perCluster_combinedResults",
  "edgeR_perCluster_testDefinitions"
)

missing_objects <- setdiff(required_objects, loaded_objects)

if (length(missing_objects) > 0L) {
  stop(
    "Required object(s) missing from edgeR RData: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

edgeR_perCluster_combinedResults <- get(
  "edgeR_perCluster_combinedResults",
  envir = result_environment,
  inherits = FALSE
)

edgeR_perCluster_testDefinitions <- get(
  "edgeR_perCluster_testDefinitions",
  envir = result_environment,
  inherits = FALSE
)

missing_test_ids <- effect_definitions$test_id[
  !effect_definitions$test_id %in% names(edgeR_perCluster_combinedResults)
]

if (length(missing_test_ids) > 0L) {
  stop(
    "Requested edgeR test(s) missing: ",
    paste(missing_test_ids, collapse = ", "),
    call. = FALSE
  )
}

# ==============================================================================
# 9. Output directory
# ==============================================================================

parameter_label <- build_parameter_label(
  fdr_threshold = fdr_threshold,
  abs_log2fc_threshold = abs_log2fc_threshold,
  min_max_mean_percent = min_max_mean_percent,
  max_max_mean_percent = max_max_mean_percent
)

parameter_output_dir <- file.path(
  volcano_output_root,
  parameter_label
)

dir.create(
  parameter_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("\nOutput directory:")
message(parameter_output_dir)

# ==============================================================================
# 10. Prepare effect table
# ==============================================================================

prepare_effect_table <- function(
  result_table,
  effect_key,
  fdr_threshold,
  abs_log2fc_threshold,
  min_max_mean_percent = NULL,
  max_max_mean_percent = NULL
) {

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
    colnames(result_table)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "edgeR result table for effect '",
      effect_key,
      "' is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  result_table <- result_table |>
    tibble::as_tibble() |>
    dplyr::mutate(
      cluster_id = as.character(.data$cluster_id),
      gene_label = dplyr::if_else(
        is.na(.data$gene) | .data$gene == "",
        as.character(.data$ensembl_gene_id),
        as.character(.data$gene)
      )
    )

  percent_columns <- resolve_percent_columns(
    result_table = result_table,
    effect_key = effect_key
  )

  result_table <- result_table |>
    dplyr::mutate(
      max_mean_percent = calculate_max_mean_percent(
        result_table = result_table,
        percent_columns = percent_columns
      )
    ) |>
    dplyr::filter(
      is.finite(.data$logFC),
      is.finite(.data$FDR),
      .data$FDR >= 0,
      is.finite(.data$max_mean_percent)
    )

  if (!is.null(min_max_mean_percent)) {
    result_table <- result_table |>
      dplyr::filter(
        .data$max_mean_percent >= .env$min_max_mean_percent
      )
  }

  if (!is.null(max_max_mean_percent)) {

    if (as.numeric(max_max_mean_percent) >= 100) {
      result_table <- result_table |>
        dplyr::filter(
          .data$max_mean_percent <= .env$max_max_mean_percent
        )
    } else {
      result_table <- result_table |>
        dplyr::filter(
          .data$max_mean_percent < .env$max_max_mean_percent
        )
    }
  }

  result_table |>
    dplyr::mutate(
      FDR_for_plot = pmax(
        .data$FDR,
        .Machine$double.xmin
      ),
      minus_log10_FDR = -log10(.data$FDR_for_plot),
      significance = dplyr::case_when(
        .data$FDR < .env$fdr_threshold &
          .data$logFC >= .env$abs_log2fc_threshold ~ "UP",
        .data$FDR < .env$fdr_threshold &
          .data$logFC <= -.env$abs_log2fc_threshold ~ "DOWN",
        TRUE ~ "Not significant"
      ),
      significance = factor(
        .data$significance,
        levels = c("DOWN", "Not significant", "UP")
      ),
      cluster_facet = vapply(
        .data$cluster_id,
        make_cluster_facet_label,
        FUN.VALUE = character(1)
      )
    )
}

# ==============================================================================
# 11. Label selection
# ==============================================================================

select_label_table <- function(
  prepared_table,
  top_n_labels_up,
  top_n_labels_down
) {

  top_up <- prepared_table |>
    dplyr::filter(.data$significance == "UP") |>
    dplyr::group_by(.data$cluster_id) |>
    dplyr::arrange(
      .data$FDR,
      dplyr::desc(abs(.data$logFC)),
      dplyr::desc(.data$logFC),
      .data$gene_label,
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = top_n_labels_up) |>
    dplyr::ungroup()

  top_down <- prepared_table |>
    dplyr::filter(.data$significance == "DOWN") |>
    dplyr::group_by(.data$cluster_id) |>
    dplyr::arrange(
      .data$FDR,
      dplyr::desc(abs(.data$logFC)),
      .data$logFC,
      .data$gene_label,
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = top_n_labels_down) |>
    dplyr::ungroup()

  dplyr::bind_rows(top_up, top_down)
}

# ==============================================================================
# 12. Build one faceted grid plot
# ==============================================================================

build_effect_grid_plot <- function(
  prepared_table,
  effect_definition,
  comparison_text,
  with_labels,
  fdr_threshold,
  abs_log2fc_threshold,
  min_max_mean_percent,
  max_max_mean_percent,
  top_n_labels_up,
  top_n_labels_down,
  facet_ncol
) {

  cluster_ids_sorted <- sort_cluster_ids(
    prepared_table$cluster_id
  )

  facet_levels <- vapply(
    cluster_ids_sorted,
    make_cluster_facet_label,
    FUN.VALUE = character(1)
  )

  prepared_table <- prepared_table |>
    dplyr::mutate(
      cluster_facet = factor(
        .data$cluster_facet,
        levels = facet_levels
      )
    )

  label_table <- if (isTRUE(with_labels)) {
    select_label_table(
      prepared_table = prepared_table,
      top_n_labels_up = top_n_labels_up,
      top_n_labels_down = top_n_labels_down
    ) |>
      dplyr::mutate(
        cluster_facet = factor(
          .data$cluster_facet,
          levels = facet_levels
        )
      )
  } else {
    prepared_table[0, , drop = FALSE]
  }

  n_up <- sum(prepared_table$significance == "UP", na.rm = TRUE)
  n_down <- sum(prepared_table$significance == "DOWN", na.rm = TRUE)

  threshold_y <- -log10(fdr_threshold)

  max_abs_logfc <- max(
    c(abs(prepared_table$logFC), abs_log2fc_threshold),
    na.rm = TRUE
  )

  x_limit <- max(
    abs_log2fc_threshold * 1.25,
    max_abs_logfc * 1.03
  )

  max_y <- max(
    c(prepared_table$minus_log10_FDR, threshold_y),
    na.rm = TRUE
  )

  y_limit <- max(
    threshold_y * 1.25,
    max_y * 1.03
  )

  version_label <- if (isTRUE(with_labels)) {
    paste0(
      "grid with labels: top ",
      top_n_labels_up,
      " UP + top ",
      top_n_labels_down,
      " DOWN per cluster"
    )
  } else {
    "grid without labels"
  }

  percent_filter_label <- format_percent_filter_label(
    min_max_mean_percent = min_max_mean_percent,
    max_max_mean_percent = max_max_mean_percent
  )

  plot_object <- ggplot2::ggplot(
    prepared_table,
    ggplot2::aes(
      x = .data$logFC,
      y = .data$minus_log10_FDR
    )
  ) +
    ggplot2::geom_hline(
      yintercept = threshold_y,
      linetype = threshold_linetype,
      colour = "black",
      linewidth = threshold_linewidth
    ) +
    ggplot2::geom_vline(
      xintercept = c(
        -abs_log2fc_threshold,
        abs_log2fc_threshold
      ),
      linetype = threshold_linetype,
      colour = "black",
      linewidth = threshold_linewidth
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$significance),
      size = point_size,
      alpha = point_alpha,
      shape = point_shape,
      stroke = point_stroke
    ) +
    ggplot2::facet_wrap(
      ~ cluster_facet,
      ncol = facet_ncol
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "DOWN" = colour_down,
        "Not significant" = colour_not_significant,
        "UP" = colour_up
      ),
      breaks = c("UP", "DOWN", "Not significant"),
      labels = c("UP", "DOWN", "Not significant"),
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-x_limit, x_limit),
      ylim = c(0, y_limit),
      clip = "off"
    ) +
    ggplot2::labs(
      title = effect_definition$effect_title[[1]],
      subtitle = paste0(
        comparison_text,
        "\n",
        "FDR < ",
        format(fdr_threshold, trim = TRUE),
        " | |log2FC| >= ",
        format(abs_log2fc_threshold, trim = TRUE),
        " | ",
        percent_filter_label,
        "\n",
        "all genes after percent filtering: ",
        format(nrow(prepared_table), big.mark = ",", scientific = FALSE, trim = TRUE),
        " | significant UP: ",
        n_up,
        " | significant DOWN: ",
        n_down,
        " | ",
        version_label
      ),
      x = "log2FC",
      y = expression(-log[10](FDR))
    ) +
    ggplot2::theme_classic(
      base_family = "DejaVu Sans",
      base_size = 11
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 5)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10.5,
        hjust = 0.5,
        lineheight = 1.06,
        margin = ggplot2::margin(b = 8)
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        size = 8.8,
        face = "bold",
        lineheight = 0.95
      ),
      axis.title = ggplot2::element_text(
        size = 12,
        face = "bold"
      ),
      axis.text = ggplot2::element_text(
        size = 8.7
      ),
      axis.line = ggplot2::element_line(
        linewidth = 0.35
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.35
      ),
      legend.position = "bottom",
      legend.justification = "center",
      legend.text = ggplot2::element_text(
        size = 10
      ),
      panel.spacing = grid::unit(0.85, "lines"),
      plot.margin = ggplot2::margin(
        t = 10,
        r = 8,
        b = 10,
        l = 8
      )
    )

  if (isTRUE(with_labels) && nrow(label_table) > 0L) {
    plot_object <- plot_object +
      ggrepel::geom_text_repel(
        data = label_table,
        ggplot2::aes(label = .data$gene_label),
        size = label_text_size,
        family = "DejaVu Sans",
        box.padding = label_box_padding,
        point.padding = label_point_padding,
        min.segment.length = 0,
        segment.colour = label_segment_colour,
        segment.linewidth = label_segment_linewidth,
        max.overlaps = Inf,
        seed = label_seed,
        show.legend = FALSE
      )
  }

  plot_object
}

# ==============================================================================
# 13. Write one PDF per effect, with 2 pages
# ==============================================================================

write_effect_pdf <- function(
  effect_definition,
  result_table,
  test_definitions,
  output_file
) {

  effect_key <- effect_definition$effect_key[[1]]
  test_id <- effect_definition$test_id[[1]]

  prepared_table <- prepare_effect_table(
    result_table = result_table,
    effect_key = effect_key,
    fdr_threshold = fdr_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold,
    min_max_mean_percent = min_max_mean_percent,
    max_max_mean_percent = max_max_mean_percent
  )

  if (nrow(prepared_table) == 0L) {
    stop(
      "Effect '",
      effect_key,
      "' contains zero rows after filtering.",
      call. = FALSE
    )
  }

  matching_definition <- test_definitions |>
    dplyr::filter(
      .data$test_id == .env$test_id
    )

  comparison_text <- if (nrow(matching_definition) == 1L) {
    as.character(matching_definition$comparison[[1]])
  } else {
    test_id
  }

  message("\n============================================================")
  message("Effect: ", effect_definition$effect_title[[1]])
  message("Writing PDF: ", output_file)

  cluster_summary <- prepared_table |>
    dplyr::group_by(.data$cluster_id) |>
    dplyr::summarise(
      n_genes = dplyr::n(),
      n_up = sum(.data$significance == "UP"),
      n_down = sum(.data$significance == "DOWN"),
      .groups = "drop"
    ) |>
    dplyr::arrange(as.integer(.data$cluster_id))

  print(cluster_summary, n = Inf)

  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  pdf_device(
    file = output_file,
    width = pdf_width,
    height = pdf_height,
    onefile = TRUE
  )

  device_was_closed <- FALSE

  on.exit(
    {
      if (!device_was_closed) {
        grDevices::dev.off()
      }
    },
    add = TRUE
  )

  plot_without_labels <- build_effect_grid_plot(
    prepared_table = prepared_table,
    effect_definition = effect_definition,
    comparison_text = comparison_text,
    with_labels = FALSE,
    fdr_threshold = fdr_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold,
    min_max_mean_percent = min_max_mean_percent,
    max_max_mean_percent = max_max_mean_percent,
    top_n_labels_up = top_n_labels_up,
    top_n_labels_down = top_n_labels_down,
    facet_ncol = facet_ncol
  )

  print(plot_without_labels)

  plot_with_labels <- build_effect_grid_plot(
    prepared_table = prepared_table,
    effect_definition = effect_definition,
    comparison_text = comparison_text,
    with_labels = TRUE,
    fdr_threshold = fdr_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold,
    min_max_mean_percent = min_max_mean_percent,
    max_max_mean_percent = max_max_mean_percent,
    top_n_labels_up = top_n_labels_up,
    top_n_labels_down = top_n_labels_down,
    facet_ncol = facet_ncol
  )

  print(plot_with_labels)

  grDevices::dev.off()
  device_was_closed <- TRUE

  if (!file.exists(output_file) || file.info(output_file)$size == 0) {
    stop(
      "PDF was not created correctly: ",
      output_file,
      call. = FALSE
    )
  }

  invisible(
    list(
      output_file = output_file,
      cluster_summary = cluster_summary
    )
  )
}

# ==============================================================================
# 14. Define output files
# ==============================================================================

effect_output_files <- setNames(
  file.path(
    parameter_output_dir,
    paste0(
      sprintf("%02d", effect_definitions$effect_order),
      "_",
      effect_definitions$file_stub,
      "_volcanoPerCluster_grid_twoVersions.pdf"
    )
  ),
  effect_definitions$effect_key
)

# ==============================================================================
# 15. Run all three effects
# ==============================================================================

effect_outputs <- vector(
  mode = "list",
  length = nrow(effect_definitions)
)

names(effect_outputs) <- effect_definitions$effect_key

for (effect_index in seq_len(nrow(effect_definitions))) {

  effect_definition_current <- effect_definitions[
    effect_index,
    ,
    drop = FALSE
  ]

  effect_key_current <- effect_definition_current$effect_key[[1]]
  test_id_current <- effect_definition_current$test_id[[1]]

  effect_outputs[[effect_key_current]] <- write_effect_pdf(
    effect_definition = effect_definition_current,
    result_table = edgeR_perCluster_combinedResults[[test_id_current]],
    test_definitions = edgeR_perCluster_testDefinitions,
    output_file = effect_output_files[[effect_key_current]]
  )
}

# ==============================================================================
# 16. Final validation
# ==============================================================================

expected_output_files <- unname(effect_output_files)

missing_output_files <- expected_output_files[
  !file.exists(expected_output_files)
]

if (length(missing_output_files) > 0L) {
  stop(
    "Missing output PDF(s):\n",
    paste(missing_output_files, collapse = "\n"),
    call. = FALSE
  )
}

empty_output_files <- expected_output_files[
  file.info(expected_output_files)$size == 0
]

if (length(empty_output_files) > 0L) {
  stop(
    "Empty output PDF(s):\n",
    paste(empty_output_files, collapse = "\n"),
    call. = FALSE
  )
}

message("\n============================================================")
message("Completed volcano grid workflow.")
message("Output directory:")
message(normalizePath(parameter_output_dir, mustWork = TRUE))
message("\nGenerated PDFs:")

for (output_file in expected_output_files) {
  message(normalizePath(output_file, mustWork = TRUE))
}

message("\nEach PDF has exactly 2 pages:")
message("  page 1 = grid without labels")
message("  page 2 = grid with labels")
message("============================================================")

# ==============================================================================
# End
# ==============================================================================
