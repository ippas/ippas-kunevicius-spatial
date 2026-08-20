#!/usr/bin/env Rscript

# ==============================================================================
# 05_plotMarkerGeneVisualizationOnSlide_highTau.R
#
# Purpose:
# Select FindAllMarkers genes with Tau above the configured threshold and
# generate five marker-gene visualizations for each selected gene.
#
# The function module does not apply a Tau threshold. Marker selection is
# performed only in this calling script.
#
# If one gene is a FindAllMarkers marker for several clusters, it is processed
# once and all corresponding cluster IDs are supplied together.
# ==============================================================================


# ==============================================================================
# 1. Packages
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
  library(scales)
  library(grid)
})


# ==============================================================================
# 2. Configuration
# ==============================================================================

options(
  stringsAsFactors = FALSE
)

project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

dataset_name <- "maternalFMT_n16samples"

configuration_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

cluster_column <- "leiden_res0.40"

assay_name <- "RNA"

clustering_parameters_label <- paste0(
  "LogNormalize + VST | HVG 2000 | ",
  "CCA dims 20 | k.param 20 | ",
  "prune.SNN 0.0667 | Leiden resolution 0.40"
)

tau_threshold <- 0.60

# The request was "above 0.8", therefore the filter is Tau > 0.8.
tau_comparison_label <- paste0(
  "tauAbove",
  gsub(
    "\\.",
    "",
    formatC(
      tau_threshold,
      format = "f",
      digits = 2
    )
  )
)


# For the first test run only two genes are processed.
# Set to NULL to process every selected high-Tau gene.
maximum_number_of_genes <- NULL


# Optional manual gene selection.
#
# NULL:
#   take the highest-Tau genes automatically.
#
# Character vector:
#   process only the requested genes, for example:
#   manual_target_genes <- c("Ace", "Hcrt")
manual_target_genes <- NULL


# Optional anatomical cluster names.
# Original numerical cluster IDs are always retained.
cluster_names <- NULL

# Example:
# cluster_names <- c(
#   "1" = "Thalamus 1",
#   "2" = "Cortex 1"
# )


normalization_scale_factor <- 10000

upper_colour_quantile <- 0.99

show_histology_image <- FALSE

skip_completed_genes <- TRUE

continue_after_gene_error <- TRUE


# ==============================================================================
# 3. Colour settings
# ==============================================================================

red_palette_colors <- c(
  "#D9D9D9",
  "#FEE5D9",
  "#FCAE91",
  "#FB6A4A",
  "#DE2D26",
  "#A50F15",
  "#67000D"
)

green_palette_colors <- c(
  "#D9D9D9",
  "#E5F5E0",
  "#A1D99B",
  "#74C476",
  "#31A354",
  "#006D2C",
  "#00441B"
)

palette_values <- c(
  0.000,
  0.006,
  0.060,
  0.180,
  0.400,
  0.700,
  1.000
)

target_cluster_color <- "#B30000"

other_cluster_color <- "#1B9E77"


# ==============================================================================
# 4. Plot sizes
# ==============================================================================

point_size_no_image <- 1.1

point_size_with_image <- 0.80

spatial_pdf_width <- 18

spatial_pdf_height <- 23

barplot_pdf_width <- 15

barplot_pdf_height <- 8


# ==============================================================================
# 5. Paths
# ==============================================================================

functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "marker_gene_visualization",
  "functions_markerGeneVisualization.R"
)

input_rdata_file <- file.path(
  project_root,
  "results",
  dataset_name,
  "seurat_clustering_analysis",
  "logNormalize_vst",
  "hvg2000",
  "02_cca",
  "dims20",
  "k20",
  "prune0067",
  "RData",
  paste0(
    "02_",
    dataset_name,
    "_logNormalizeVst_hvg2000_",
    "ccaDims20K20PruneSNN0067_",
    "res010to100by010_multiClusteringAndUmap.RData"
  )
)

input_specificity_file <- file.path(
  project_root,
  "results",
  dataset_name,
  "manual_cluster_annotation",
  configuration_name,
  "02_geneSpecificity_onFindMarkers",
  paste0(
    "01_",
    dataset_name,
    "_leidenRes040_",
    "findAllMarkers_withGeneSpecificity.tsv"
  )
)

output_root_directory <- file.path(
  project_root,
  "results",
  dataset_name,
  "manual_cluster_annotation",
  configuration_name,
  "04_markerGeneVisualizationOnSlide_highTau"
)

selected_genes_file <- file.path(
  output_root_directory,
  paste0(
    "00_",
    dataset_name,
    "_leidenRes040_",
    "selectedMarkerGenes_",
    tau_comparison_label,
    ".tsv"
  )
)

run_status_file <- file.path(
  output_root_directory,
  paste0(
    "00_",
    dataset_name,
    "_leidenRes040_",
    "markerGeneVisualization_runStatus.tsv"
  )
)

if (!file.exists(functions_file)) {
  stop(
    "Functions file does not exist:\n",
    functions_file,
    call. = FALSE
  )
}

if (!file.exists(input_rdata_file)) {
  stop(
    "Input RData file does not exist:\n",
    input_rdata_file,
    call. = FALSE
  )
}

if (!file.exists(input_specificity_file)) {
  stop(
    "Input specificity file does not exist:\n",
    input_specificity_file,
    call. = FALSE
  )
}

dir.create(
  output_root_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

source(
  functions_file
)


# ==============================================================================
# 6. Load the Seurat object
# ==============================================================================

load_environment <- new.env(
  parent = globalenv()
)

loaded_object_names <- load(
  input_rdata_file,
  envir = load_environment
)

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
    "Expected exactly one Seurat object in the RData file.\n",
    "Seurat objects found: ",
    paste(
      seurat_object_names,
      collapse = ", "
    ),
    call. = FALSE
  )
}

seurat_object <- get(
  seurat_object_names[[1]],
  envir = load_environment,
  inherits = FALSE
)

rm(
  load_environment
)

message(
  "Loaded Seurat object: ",
  seurat_object_names[[1]]
)

message(
  "Object dimensions: ",
  nrow(seurat_object),
  " genes x ",
  ncol(seurat_object),
  " spots"
)


# ==============================================================================
# 7. Select high-Tau marker genes
# ==============================================================================

specificity_table <- read.delim(
  input_specificity_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_specificity_columns <- c(
  "gene",
  "cluster",
  "tau"
)

missing_specificity_columns <- setdiff(
  required_specificity_columns,
  colnames(specificity_table)
)

if (length(missing_specificity_columns) > 0L) {
  stop(
    "Missing specificity column(s): ",
    paste(
      missing_specificity_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

specificity_table$gene <- as.character(
  specificity_table$gene
)

specificity_table$cluster <- as.character(
  specificity_table$cluster
)

high_tau_marker_table <- specificity_table |>
  dplyr::filter(
    is.finite(tau),
    tau > tau_threshold
  )

if (nrow(high_tau_marker_table) == 0L) {
  stop(
    "No FindAllMarkers genes passed Tau > ",
    tau_threshold,
    ".",
    call. = FALSE
  )
}

specificity_metric_columns <- intersect(
  c(
    "gene",
    "cluster",
    "tau",
    "gini",
    "shannon_specificity",
    "expression_specificity",
    "expression_ratio_vs_best_other",
    "best_other_cluster",
    "mean_expression_best_other",
    "is_best_cluster"
  ),
  colnames(
    high_tau_marker_table
  )
)

gene_selection_table <- high_tau_marker_table |>
  dplyr::group_by(
    gene
  ) |>
  dplyr::summarise(
    tau = max(
      tau,
      na.rm = TRUE
    ),
    marker_clusters = list(
      sort_marker_cluster_ids(
        cluster
      )
    ),
    marker_metrics = list(
      dplyr::pick(
        dplyr::all_of(
          setdiff(
            specificity_metric_columns,
            "gene"
          )
        )
      ) |>
        dplyr::distinct(
          cluster,
          .keep_all = TRUE
        ) |>
        dplyr::arrange(
          suppressWarnings(
            as.numeric(
              cluster
            )
          )
        )
    ),
    number_of_marker_clusters =
      dplyr::n_distinct(
        cluster
      ),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(tau),
    gene
  )

if (!is.null(manual_target_genes)) {

  manual_target_genes <- unique(
    as.character(
      manual_target_genes
    )
  )

  missing_manual_genes <- setdiff(
    manual_target_genes,
    gene_selection_table$gene
  )

  if (length(missing_manual_genes) > 0L) {
    stop(
      "Manual gene(s) absent from the Tau-selected marker table: ",
      paste(
        missing_manual_genes,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  gene_selection_table <-
    gene_selection_table[
      match(
        manual_target_genes,
        gene_selection_table$gene
      ),
      ,
      drop = FALSE
    ]
}

if (!is.null(maximum_number_of_genes)) {

  maximum_number_of_genes <- as.integer(
    maximum_number_of_genes
  )

  if (maximum_number_of_genes < 1L) {
    stop(
      "`maximum_number_of_genes` must be NULL or a positive integer.",
      call. = FALSE
    )
  }

  gene_selection_table <- utils::head(
    gene_selection_table,
    n = maximum_number_of_genes
  )
}

if (nrow(gene_selection_table) == 0L) {
  stop(
    "No genes remain after applying the selection settings.",
    call. = FALSE
  )
}

selected_genes_output <- gene_selection_table |>
  dplyr::mutate(
    marker_clusters = vapply(
      marker_clusters,
      paste,
      character(1),
      collapse = ","
    ),
    specificity_summary = vapply(
      marker_metrics,
      format_marker_specificity_summary,
      character(1)
    )
  ) |>
  dplyr::select(
    -marker_metrics
  )

write_marker_tsv(
  selected_genes_output,
  selected_genes_file
)

message(
  "Genes selected for this run: ",
  nrow(gene_selection_table)
)

message(
  paste0(
    "  ",
    gene_selection_table$gene,
    " | Tau=",
    formatC(
      gene_selection_table$tau,
      format = "f",
      digits = 3
    ),
    " | marker cluster(s)=",
    vapply(
      gene_selection_table$marker_clusters,
      paste,
      character(1),
      collapse = ","
    ),
    collapse = "\n"
  )
)


# ==============================================================================
# 8. Prepare shared spatial and count context once
# ==============================================================================

message(
  "\nPreparing shared Visium context..."
)

visualization_context <-
  prepare_marker_gene_visualization_context(
    seurat_object = seurat_object,
    cluster_column = cluster_column,
    assay_name = assay_name,
    sample_id_column = "sample_ID",
    group_column = "fmt_donor_group",
    sex_column = "sex",
    image_scale = "lowres",
    panel_padding_fraction = 0.03
  )

rm(
  seurat_object
)

invisible(
  gc(
    verbose = FALSE
  )
)

message(
  "Shared context prepared for ",
  length(
    visualization_context$sample_layout$sample_order
  ),
  " samples and ",
  length(
    visualization_context$cluster_levels
  ),
  " clusters."
)


# ==============================================================================
# 9. Process selected genes
# ==============================================================================

run_started_at <- Sys.time()

status_rows <- list()

for (
  gene_index in
  seq_len(
    nrow(gene_selection_table)
  )
) {

  target_gene <-
    gene_selection_table$gene[[
      gene_index
    ]]

  marker_clusters <-
    gene_selection_table$marker_clusters[[
      gene_index
    ]]

  marker_metrics_table <-
    gene_selection_table$marker_metrics[[
      gene_index
    ]]

  gene_started_at <- Sys.time()

  message(
    "\n[",
    gene_index,
    "/",
    nrow(gene_selection_table),
    "] Processing ",
    target_gene,
    " | marker cluster(s): ",
    paste(
      marker_clusters,
      collapse = ", "
    ),
    "
",
    format_marker_specificity_summary(
      marker_metrics_table
    ),
    "
Output folder: ",
    paste0(
      target_gene,
      "_",
      build_marker_cluster_file_label(
        marker_clusters
      )
    )
  )

  result <- tryCatch(
    {
      save_marker_gene_visualizations(
        context =
          visualization_context,
        target_gene =
          target_gene,
        marker_clusters =
          marker_clusters,
        marker_metrics_table =
          marker_metrics_table,
        clustering_parameters_label =
          clustering_parameters_label,
        output_root_directory =
          output_root_directory,
        normalization_scale_factor =
          normalization_scale_factor,
        upper_colour_quantile =
          upper_colour_quantile,
        red_palette_colors =
          red_palette_colors,
        green_palette_colors =
          green_palette_colors,
        palette_values =
          palette_values,
        target_cluster_color =
          target_cluster_color,
        other_cluster_color =
          other_cluster_color,
        cluster_names =
          cluster_names,
        show_histology_image =
          show_histology_image,
        point_size_no_image =
          point_size_no_image,
        point_size_with_image =
          point_size_with_image,
        spatial_pdf_width =
          spatial_pdf_width,
        spatial_pdf_height =
          spatial_pdf_height,
        barplot_pdf_width =
          barplot_pdf_width,
        barplot_pdf_height =
          barplot_pdf_height,
        skip_completed =
          skip_completed_genes
      )
    },
    error = function(error_condition) {
      error_condition
    }
  )

  gene_finished_at <- Sys.time()

  elapsed_seconds <- as.numeric(
    difftime(
      gene_finished_at,
      gene_started_at,
      units = "secs"
    )
  )

  if (inherits(result, "error")) {

    current_status <- "failed"
    current_message <- conditionMessage(
      result
    )

    message(
      "FAILED: ",
      target_gene,
      "\n",
      current_message
    )

  } else {

    current_status <- "completed"
    current_message <- ""

    message(
      "Completed: ",
      target_gene
    )
  }

  status_rows[[
    length(status_rows) + 1L
  ]] <- data.frame(
    order = gene_index,
    total_genes =
      nrow(gene_selection_table),
    target_gene =
      target_gene,
    tau =
      gene_selection_table$tau[[
        gene_index
      ]],
    marker_clusters =
      paste(
        marker_clusters,
        collapse = ","
      ),
    output_subdirectory =
      paste0(
        target_gene,
        "_",
        build_marker_cluster_file_label(
          marker_clusters
        )
      ),
    status =
      current_status,
    message =
      current_message,
    elapsed_seconds =
      round(
        elapsed_seconds,
        2
      ),
    started_at =
      format(
        gene_started_at,
        "%Y-%m-%d %H:%M:%S"
      ),
    finished_at =
      format(
        gene_finished_at,
        "%Y-%m-%d %H:%M:%S"
      ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  current_status_table <- do.call(
    rbind,
    status_rows
  )

  write_marker_tsv(
    current_status_table,
    run_status_file
  )

  if (
    current_status == "failed" &&
      !isTRUE(
        continue_after_gene_error
      )
  ) {
    stop(
      current_message,
      call. = FALSE
    )
  }

  invisible(
    gc(
      verbose = FALSE
    )
  )
}


# ==============================================================================
# 10. Final summary
# ==============================================================================

run_finished_at <- Sys.time()

final_status_table <- do.call(
  rbind,
  status_rows
)

write_marker_tsv(
  final_status_table,
  run_status_file
)

message(
  "\n============================================================"
)

message(
  "Finished marker-gene visualization run."
)

message(
  "Completed: ",
  sum(
    final_status_table$status ==
      "completed"
  )
)

message(
  "Failed: ",
  sum(
    final_status_table$status ==
      "failed"
  )
)

message(
  "Elapsed time: ",
  round(
    as.numeric(
      difftime(
        run_finished_at,
        run_started_at,
        units = "mins"
      )
    ),
    2
  ),
  " min"
)

message(
  "Output directory:\n",
  output_root_directory
)

message(
  "============================================================"
)
