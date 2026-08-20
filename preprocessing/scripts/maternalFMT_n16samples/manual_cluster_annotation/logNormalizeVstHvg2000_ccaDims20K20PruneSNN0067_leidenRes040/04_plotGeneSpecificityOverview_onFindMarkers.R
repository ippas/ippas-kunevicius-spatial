# ==============================================================================
# 04_plotGeneSpecificityOverview_onFindMarkers.R
#
# Purpose:
# Create two 16-panel overview plots for genes identified by FindAllMarkers():
# 1. number of genes passing consecutive specificity thresholds;
# 2. percentage of FindAllMarkers genes passing consecutive thresholds.
#
# Three specificity metrics are shown in every panel:
# - tau,
# - gini,
# - shannon_specificity.
# ==============================================================================


# ==============================================================================
# 1. Configuration
# ==============================================================================

options(
  stringsAsFactors = FALSE
)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
})

project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

dataset_name <- "maternalFMT_n16samples"

configuration_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

cluster_column <- "leiden_res0.40"

specificity_metrics <- c(
  "tau",
  "gini",
  "shannon_specificity"
)

specificity_thresholds <- seq(
  0,
  1,
  by = 0.1
)


# Optional anatomical names.
# Keep NULL to show only the original cluster numbers.
cluster_names <- NULL


# Colours of squares shown before cluster names in panel titles.
custom_cluster_colors <- c(
  # Cortex
  "2"  = "#407ba7",
  "8"  = "#004e89",
  "7"  = "#002962",

  # Hippocampus
  "3"  = "#31cb00",
  "13" = "#119822",
  "16" = "#1e441e",

  # Thalamus
  "1"  = "#e66063",
  "9"  = "#d02224",
  "11" = "#9c191b",

  "6"  = "#b5838d",
  "14" = "#6B1E2D",
  "4"  = "#adb5bd",
  "10" = "#ff9505",
  "12" = "#e85d04",
  "5"  = "#9e0059",
  "15" = "#ffd100"
)


# Colours of the three specificity curves.
specificity_metric_colors <- c(
  "tau" = "#1B9E77",
  "gini" = "#D95F02",
  "shannon_specificity" = "#7570B3"
)

specificity_metric_labels <- c(
  "tau" = "Tau",
  "gini" = "Gini",
  "shannon_specificity" =
    "Shannon specificity"
)


# Plot appearance.
# Available themes:
# "bw", "classic", "minimal", "light", "grey"
plot_theme_style <- "bw"

legend_position <- "top"
legend_title_size <- 14
legend_text_size <- 13
legend_point_size <- 5.5
legend_line_width <- 2.2
legend_key_width_cm <- 2.0
legend_key_height_cm <- 0.8
legend_height_ratio <- 0.05

cluster_title_symbol <- "■"

# Increased strongly on purpose.
cluster_title_symbol_size_pt <- 76


# ==============================================================================
# 2. Paths
# ==============================================================================

spatial_functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "manual_cluster_annotation",
  "functions_spatialClusterVisualization_manualAnnotation.R"
)

specificity_plot_functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "gene_specificity",
  "functions_geneSpecificityOverview_onFindMarkers.R"
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

output_directory <- file.path(
  project_root,
  "results",
  dataset_name,
  "manual_cluster_annotation",
  configuration_name,
  "03_geneSpecificityOverview_onFindMarkers"
)

output_count_pdf <- file.path(
  output_directory,
  paste0(
    "01_",
    dataset_name,
    "_leidenRes040_",
    "tauGiniShannonThreshold_",
    "numberOfGenes.pdf"
  )
)

output_percentage_pdf <- file.path(
  output_directory,
  paste0(
    "02_",
    dataset_name,
    "_leidenRes040_",
    "tauGiniShannonThreshold_",
    "percentageOfFindMarkersGenes.pdf"
  )
)

output_threshold_summary <- file.path(
  output_directory,
  paste0(
    "03_",
    dataset_name,
    "_leidenRes040_",
    "tauGiniShannonThreshold_",
    "summary.tsv"
  )
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

source(
  spatial_functions_file
)

source(
  specificity_plot_functions_file
)


# ==============================================================================
# 3. Load Seurat object
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
    "Expected exactly one Seurat object in the RData file.",
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


# ==============================================================================
# 4. Load gene-specificity results
# ==============================================================================

markers_with_specificity <- read.delim(
  input_specificity_file,
  check.names = FALSE
)

markers_with_specificity$cluster <- as.character(
  markers_with_specificity$cluster
)


# ==============================================================================
# 5. Build cluster summary using all samples in the Seurat object
# ==============================================================================

cluster_summary <- build_cluster_summary_table(
  seurat_object = seurat_object,
  cluster_column = cluster_column,
  sample_order = NULL
)$cluster_summary


# ==============================================================================
# 6. Plot titles
# ==============================================================================

plot_subtitle <- paste0(
  "LogNormalize + VST | HVG 2000 | CCA dims 20 | ",
  "k.param 20 | prune.SNN 0.0667 | Leiden resolution 0.40 | ",
  "FindAllMarkers: min.pct 0.25, log2FC 0.25, adjusted P ≤ 0.05 | ",
  "specificity threshold: 0.0–1.0, step 0.1"
)


# ==============================================================================
# 7. Number of genes
# ==============================================================================

count_result <- plot_gene_specificity_threshold_grid(
  marker_table =
    markers_with_specificity,
  cluster_summary =
    cluster_summary,
  custom_cluster_colors =
    custom_cluster_colors,
  cluster_names =
    cluster_names,
  metric_columns =
    specificity_metrics,
  metric_labels =
    specificity_metric_labels,
  metric_colors =
    specificity_metric_colors,
  thresholds =
    specificity_thresholds,
  y_mode =
    "count",
  ncol =
    4L,
  plot_title =
    paste0(
      "Gene specificity among FindAllMarkers genes: ",
      "number of genes"
    ),
  plot_subtitle =
    plot_subtitle,
  theme_style =
    plot_theme_style,
  legend_position =
    legend_position,
  legend_title_size =
    legend_title_size,
  legend_text_size =
    legend_text_size,
  legend_point_size =
    legend_point_size,
  legend_line_width =
    legend_line_width,
  legend_key_width_cm =
    legend_key_width_cm,
  legend_key_height_cm =
    legend_key_height_cm,
  legend_height_ratio =
    legend_height_ratio,
  cluster_symbol =
    cluster_title_symbol,
  cluster_symbol_size_pt =
    cluster_title_symbol_size_pt
)

ggplot2::ggsave(
  filename = output_count_pdf,
  plot = count_result$plot,
  width = 18,
  height = 16,
  units = "in",
  bg = "white",
  limitsize = FALSE
)


# ==============================================================================
# 8. Percentage of FindAllMarkers genes
# ==============================================================================

percentage_result <- plot_gene_specificity_threshold_grid(
  marker_table =
    markers_with_specificity,
  cluster_summary =
    cluster_summary,
  custom_cluster_colors =
    custom_cluster_colors,
  cluster_names =
    cluster_names,
  metric_columns =
    specificity_metrics,
  metric_labels =
    specificity_metric_labels,
  metric_colors =
    specificity_metric_colors,
  thresholds =
    specificity_thresholds,
  y_mode =
    "percentage",
  ncol =
    4L,
  plot_title =
    paste0(
      "Gene specificity among FindAllMarkers genes: ",
      "percentage of genes"
    ),
  plot_subtitle =
    plot_subtitle,
  theme_style =
    plot_theme_style,
  legend_position =
    legend_position,
  legend_title_size =
    legend_title_size,
  legend_text_size =
    legend_text_size,
  legend_point_size =
    legend_point_size,
  legend_line_width =
    legend_line_width,
  legend_key_width_cm =
    legend_key_width_cm,
  legend_key_height_cm =
    legend_key_height_cm,
  legend_height_ratio =
    legend_height_ratio,
  cluster_symbol =
    cluster_title_symbol,
  cluster_symbol_size_pt =
    cluster_title_symbol_size_pt
)

ggplot2::ggsave(
  filename = output_percentage_pdf,
  plot = percentage_result$plot,
  width = 18,
  height = 16,
  units = "in",
  bg = "white",
  limitsize = FALSE
)


# ==============================================================================
# 9. Save values used in plots
# ==============================================================================

write.table(
  count_result$plot_data,
  file = output_threshold_summary,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# ==============================================================================
# 10. Summary
# ==============================================================================

message(
  "Saved number-of-genes plot:\n",
  output_count_pdf
)

message(
  "Saved percentage plot:\n",
  output_percentage_pdf
)

message(
  "Saved threshold summary:\n",
  output_threshold_summary
)
