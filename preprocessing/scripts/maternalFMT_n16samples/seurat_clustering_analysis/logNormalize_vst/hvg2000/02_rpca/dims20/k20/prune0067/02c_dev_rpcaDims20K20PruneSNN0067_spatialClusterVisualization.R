# ==============================================================================
# 02c_dev_rpcaDims20K20PruneSNN0067_spatialClusterVisualization.R
#
# Combined development visualization workflow:
#
# 1. Load the Seurat object produced by 02b_dev.
# 2. Calculate one UMAP from integrated.rpca dimensions 1:20.
# 3. Generate spatial cluster plots for all 12 clustering columns.
# 4. Generate UMAP cluster plots for all 12 clustering columns.
# 5. Save every plot as PNG and PDF.
# 6. Optionally save a new RData containing the calculated UMAP reduction.
#
# Clustering columns:
# - Louvain, Louvain refined, SLM and Leiden;
# - resolutions 0.20, 0.30 and 0.40.
#
# Both plot types use:
# - one fixed neutral working palette;
# - legend labels:
#     cluster_id (n=...; ...%; mean±SD=...±...)
# - a one-line legend title:
#     <cluster_column> (global n, %, mean±SD)
# ==============================================================================


# ==============================================================================
# 0. Clean session
# ==============================================================================

rm(list = ls())

options(
  stringsAsFactors = FALSE
)


# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "patchwork",
  "dplyr",
  "uwot"
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
    "Missing required R package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})


# ==============================================================================
# 2. Analysis configuration
# ==============================================================================

dataset_name <- "maternalFMT_n16samples"

analysis_prefix <-
  "logNormalizeVst_hvg2000_rpcaDims20K20PruneSNN0067"

graph_configuration_tag <-
  "rpcaDims20K20PruneSNN0067"

normalization_label <- "LogNormalize + VST"
n_hvg <- 2000L

integration_method <- "RPCA"
integration_reduction <- "integrated.rpca"
integration_dims <- 1:20

k_param <- 20L
prune_snn <- 1 / 15

clustering_algorithms <- c(
  "louvain",
  "louvainRefined",
  "slm",
  "leiden"
)

clustering_resolutions <- c(
  0.20,
  0.30,
  0.40
)

selected_cluster_columns <- unlist(
  lapply(
    clustering_algorithms,
    function(algorithm_name) {

      paste0(
        algorithm_name,
        "_res",
        formatC(
          clustering_resolutions,
          format = "f",
          digits = 2L
        )
      )
    }
  ),
  use.names = FALSE
)

sample_order <- c(
  "1_1F",
  "1_1Fd",
  "1_1M",
  "2_1M",
  "2_1Md",
  "3_1F",
  "3_1M",
  "5_1M",
  "5_3F",
  "12_1M",
  "13_1F",
  "13_1M",
  "15_1F",
  "18_1F",
  "18_1M",
  "23_1F"
)


# ==============================================================================
# 3. UMAP configuration
# ==============================================================================

umap_reduction_name <- "umap.rpcaDims20"
umap_reduction_key <- "UMAPRPCA_"

umap_n_neighbors <- 30L
umap_min_dist <- 0.30
umap_spread <- 1
umap_metric <- "cosine"
umap_seed <- 7L

force_umap <- FALSE

umap_pt_size <- 0.05
umap_point_alpha <- 1
umap_raster <- FALSE
umap_raster_dpi <- c(
  512,
  512
)

# Fewer legend columns produce more rows and reduce artificial plot width.
umap_legend_ncol <- 3L
umap_legend_point_size <- 6
umap_legend_height_ratio <- 0.32

umap_png_width_in <- 12
umap_png_height_in <- 12

umap_pdf_width_in <- 12
umap_pdf_height_in <- 12


# ==============================================================================
# 4. Spatial-plot configuration
# ==============================================================================

show_image <- FALSE
image_alpha <- 1
crop <- FALSE
pt.size.factor <- 1.8

spatial_figure_ncol <- 4

spatial_legend_position <- "top"
spatial_legend_ncol <- 4
spatial_legend_point_size <- 6
spatial_legend_height_ratio <- 0.16

spatial_png_width_in <- 18
spatial_png_height_in <- 22

spatial_pdf_width_in <- 18
spatial_pdf_height_in <- 22


# ==============================================================================
# 5. Shared output configuration
# ==============================================================================

include_cluster_count_in_title <- TRUE

# Neutral working palette:
# - "working30": fixed 30-colour working palette;
# - "dark3": previous neutral alternative.
#
# Colours do not encode anatomical or biological meaning.
cluster_palette_name <- "working30"

save_png <- TRUE
save_pdf <- TRUE
dpi <- 300

# Save a new object containing the calculated UMAP.
# The input 02b RData is never overwritten.
save_updated_rdata <- TRUE


# ==============================================================================
# 6. Paths
# ==============================================================================

project_root <- normalizePath(
  "/home/mateusz/projects/ippas-kunevicius-spatial",
  winslash = "/",
  mustWork = TRUE
)

spatial_functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "seurat_clustering_analysis",
  "logNormalize_vst",
  "functions_spatialClusterVisualization.R"
)

umap_functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "seurat_clustering_analysis",
  "logNormalize_vst",
  "functions_umapVisualization.R"
)

analysis_root <- file.path(
  project_root,
  "results",
  dataset_name,
  "seurat_clustering_analysis",
  "logNormalize_vst",
  "hvg2000",
  "02_rpca",
  "dims20",
  "k20",
  "prune0067"
)

rdata_directory <- file.path(
  analysis_root,
  "RData"
)

input_object_name <- paste0(
  dataset_name,
  "_logNormalizeVst_hvg2000_",
  graph_configuration_tag,
  "_multiClustering"
)

input_rdata_file <- file.path(
  rdata_directory,
  paste0(
    "02b_dev_",
    dataset_name,
    "_logNormalizeVst_hvg2000_",
    graph_configuration_tag,
    "_multiClustering.RData"
  )
)

output_object_name <- paste0(
  dataset_name,
  "_logNormalizeVst_hvg2000_",
  graph_configuration_tag,
  "_multiClusteringAndUmap"
)

output_rdata_file <- file.path(
  rdata_directory,
  paste0(
    "02c_dev_",
    dataset_name,
    "_logNormalizeVst_hvg2000_",
    graph_configuration_tag,
    "_multiClusteringAndUmap.RData"
  )
)

spatial_figures_root <- file.path(
  analysis_root,
  "figures",
  "spatial_clusters"
)

umap_figures_root <- file.path(
  analysis_root,
  "figures",
  "umap_clusters"
)

dir.create(
  spatial_figures_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  umap_figures_root,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(spatial_functions_file)) {
  stop(
    "Spatial visualization functions file does not exist:\n",
    spatial_functions_file
  )
}

if (!file.exists(umap_functions_file)) {
  stop(
    "UMAP visualization functions file does not exist:\n",
    umap_functions_file
  )
}

if (!file.exists(input_rdata_file)) {
  stop(
    "Input RData file does not exist:\n",
    input_rdata_file
  )
}

source(
  spatial_functions_file
)

source(
  umap_functions_file
)


# ==============================================================================
# 7. Load and validate Seurat object
# ==============================================================================

loaded_object_names <- load(
  input_rdata_file
)

if (!input_object_name %in% loaded_object_names) {
  stop(
    "Expected Seurat object `",
    input_object_name,
    "` was not found in:\n",
    input_rdata_file,
    "\nObjects found: ",
    paste(
      loaded_object_names,
      collapse = ", "
    )
  )
}

seurat_object <- get(
  input_object_name,
  inherits = FALSE
)

if (!inherits(seurat_object, "Seurat")) {
  stop(
    "Loaded object `",
    input_object_name,
    "` is not a Seurat object."
  )
}

missing_cluster_columns <- setdiff(
  selected_cluster_columns,
  colnames(seurat_object[[]])
)

if (length(missing_cluster_columns) > 0L) {
  stop(
    "The following clustering columns are missing in the Seurat object:\n",
    paste(
      missing_cluster_columns,
      collapse = ", "
    )
  )
}

message("")
message("============================================================")
message("Loaded input Seurat object")
message("============================================================")
message("Object: ", input_object_name)
message(
  "Dimensions: ",
  nrow(seurat_object),
  " genes x ",
  ncol(seurat_object),
  " spots"
)
message("Input RData: ", input_rdata_file)
message(
  "Spatial images detected: ",
  length(
    SeuratObject::Images(seurat_object)
  )
)
message(
  "Clustering columns selected: ",
  length(selected_cluster_columns)
)
message(
  "Cluster palette: ",
  cluster_palette_name
)


# ==============================================================================
# 8. Calculate or reuse one UMAP
# ==============================================================================

umap_result <- run_umap_if_missing(
  seurat_object = seurat_object,
  input_reduction = integration_reduction,
  dims = integration_dims,
  umap_reduction_name = umap_reduction_name,
  umap_reduction_key = umap_reduction_key,
  n_neighbors = umap_n_neighbors,
  min_dist = umap_min_dist,
  spread = umap_spread,
  metric = umap_metric,
  seed_use = umap_seed,
  force_umap = force_umap,
  verbose = TRUE
)

seurat_object <- umap_result$seurat_object

rm(umap_result)

invisible(
  gc(verbose = FALSE)
)


# ==============================================================================
# 9. Shared header text
# ==============================================================================

integration_parameter_line <- paste0(
  "Integration: ",
  integration_method,
  " | normalization: ",
  normalization_label,
  " | HVG: ",
  n_hvg,
  " | dims: ",
  min(integration_dims),
  "–",
  max(integration_dims),
  " | k.param: ",
  k_param,
  " | prune.SNN: ",
  formatC(
    prune_snn,
    format = "f",
    digits = 4L
  )
)

umap_parameter_line <- paste0(
  "Integration: ",
  integration_method,
  " | normalization: ",
  normalization_label,
  " | HVG: ",
  n_hvg,
  " | dims: ",
  min(integration_dims),
  "–",
  max(integration_dims),
  " | UMAP n.neighbors: ",
  umap_n_neighbors,
  " | min.dist: ",
  formatC(
    umap_min_dist,
    format = "f",
    digits = 2L
  ),
  " | metric: ",
  umap_metric,
  " | seed: ",
  umap_seed
)


# ==============================================================================
# 10. Generate spatial and UMAP plots for all 12 clusterings
# ==============================================================================

for (
  clustering_index in
  seq_along(selected_cluster_columns)
) {

  cluster_column <-
    selected_cluster_columns[[clustering_index]]

  message("")
  message("============================================================")
  message(
    "[",
    clustering_index,
    "/",
    length(selected_cluster_columns),
    "] Processing: ",
    cluster_column
  )
  message("============================================================")


  # ---------------------------------------------------------------------------
  # 10a. Spatial plot
  # ---------------------------------------------------------------------------

  spatial_output_dir <- file.path(
    spatial_figures_root,
    cluster_column
  )

  spatial_output_prefix <-
    build_spatial_cluster_output_prefix(
      analysis_prefix = analysis_prefix,
      cluster_column = cluster_column
    )

  spatial_plot_title <- paste0(
    "Spatial clustering across 16 maternal FMT samples: ",
    cluster_column
  )

  spatial_information_line <- paste0(
    "4 × 4 panel layout | legend: global n, %, mean±SD across samples | ",
    "palette: ",
    cluster_palette_name,
    " | show_image = ",
    if (show_image) "TRUE" else "FALSE"
  )

  spatial_plot_subtitle <- paste(
    integration_parameter_line,
    spatial_information_line,
    sep = "\n"
  )

  spatial_plot_result <-
    plot_spatial_clusters_all_samples(
      seurat_object = seurat_object,
      cluster_column = cluster_column,
      sample_order = sample_order,
      ncol = spatial_figure_ncol,
      plot_title = spatial_plot_title,
      plot_subtitle = spatial_plot_subtitle,
      include_cluster_count_in_title =
        include_cluster_count_in_title,
      show_image = show_image,
      image_alpha = image_alpha,
      crop = crop,
      pt.size.factor = pt.size.factor,
      legend_position =
        spatial_legend_position,
      legend_ncol =
        spatial_legend_ncol,
      legend_point_size =
        spatial_legend_point_size,
      legend_height_ratio =
        spatial_legend_height_ratio,
      palette_name =
        cluster_palette_name,
      output_dir =
        spatial_output_dir,
      output_prefix =
        spatial_output_prefix,
      save_png = save_png,
      save_pdf = save_pdf,
      png_width_in =
        spatial_png_width_in,
      png_height_in =
        spatial_png_height_in,
      pdf_width_in =
        spatial_pdf_width_in,
      pdf_height_in =
        spatial_pdf_height_in,
      dpi = dpi,
      verbose = TRUE
    )

  message("Spatial outputs:")

  print(
    spatial_plot_result$output_files
  )

  rm(spatial_plot_result)

  invisible(
    gc(verbose = FALSE)
  )


  # ---------------------------------------------------------------------------
  # 10b. UMAP plot
  # ---------------------------------------------------------------------------

  umap_output_dir <- file.path(
    umap_figures_root,
    cluster_column
  )

  umap_output_prefix <-
    build_umap_cluster_output_prefix(
      analysis_prefix = analysis_prefix,
      cluster_column = cluster_column
    )

  umap_plot_title <- paste0(
    "UMAP clustering across 16 maternal FMT samples: ",
    cluster_column
  )

  umap_information_line <- paste0(
    "Spots: ",
    format(
      ncol(seurat_object),
      big.mark = " ",
      scientific = FALSE
    ),
    " | legend: global n, %, mean±SD across samples | ",
    "palette: ",
    cluster_palette_name
  )

  umap_plot_subtitle <- paste(
    umap_parameter_line,
    umap_information_line,
    sep = "\n"
  )

  umap_plot_result <- plot_umap_clusters(
    seurat_object = seurat_object,
    cluster_column = cluster_column,
    umap_reduction_name =
      umap_reduction_name,
    sample_order = sample_order,
    plot_title = umap_plot_title,
    plot_subtitle = umap_plot_subtitle,
    include_cluster_count_in_title =
      include_cluster_count_in_title,
    palette_name =
      cluster_palette_name,
    pt_size = umap_pt_size,
    point_alpha =
      umap_point_alpha,
    shuffle = TRUE,
    shuffle_seed = umap_seed,
    raster = umap_raster,
    raster_dpi =
      umap_raster_dpi,
    legend_ncol =
      umap_legend_ncol,
    legend_point_size =
      umap_legend_point_size,
    legend_height_ratio =
      umap_legend_height_ratio,
    output_dir =
      umap_output_dir,
    output_prefix =
      umap_output_prefix,
    save_png = save_png,
    save_pdf = save_pdf,
    png_width_in =
      umap_png_width_in,
    png_height_in =
      umap_png_height_in,
    pdf_width_in =
      umap_pdf_width_in,
    pdf_height_in =
      umap_pdf_height_in,
    dpi = dpi,
    verbose = TRUE
  )

  message("UMAP outputs:")

  print(
    umap_plot_result$output_files
  )

  rm(umap_plot_result)

  invisible(
    gc(verbose = FALSE)
  )
}


# ==============================================================================
# 11. Save a new RData containing the UMAP reduction
# ==============================================================================

if (isTRUE(save_updated_rdata)) {

  assign(
    output_object_name,
    seurat_object
  )

  save(
    list = output_object_name,
    file = output_rdata_file,
    compress = TRUE
  )

  if (
    !file.exists(output_rdata_file) ||
      is.na(
        file.info(output_rdata_file)$size
      ) ||
      file.info(output_rdata_file)$size <= 0L
  ) {
    stop(
      "The output RData file was not created correctly:\n",
      output_rdata_file
    )
  }
}


# ==============================================================================
# 12. Final report
# ==============================================================================

message("")
message("============================================================")
message("SPATIAL AND UMAP VISUALIZATION COMPLETED")
message("============================================================")
message(
  "Processed clustering columns: ",
  length(selected_cluster_columns)
)
message("Spatial figures:")
message(spatial_figures_root)
message("UMAP figures:")
message(umap_figures_root)

if (isTRUE(save_updated_rdata)) {
  message("Output RData with UMAP:")
  message(output_rdata_file)
}
