# ==============================================================================
# 01_plotSelectedClustering_fourOutputs_v3.R
#
# FOUR-OUTPUT VERSION
#
# Outputs:
# 1. Spatial clustering: cluster IDs only
# 2. UMAP clustering: cluster IDs only
# 3. Spatial clustering: cluster IDs + anatomical names
# 4. UMAP clustering: cluster IDs + anatomical names
#
# UMAP points:
# - point_alpha = 1
# - raster = FALSE
# - standard circular ggplot points
# ==============================================================================


# ==============================================================================
# 1. Configuration
# ==============================================================================

options(
  stringsAsFactors = FALSE
)

script_version <- "2026-08-02_fourOutputs_v4"

message("")
message("============================================================")
message("RUNNING: ", script_version)
message("THIS SCRIPT MUST GENERATE FOUR PDF FILES")
message("============================================================")
message("")

project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

dataset_name <- "maternalFMT_n16samples"

configuration_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

cluster_column <- "leiden_res0.40"

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
# 2. Cluster colours and names
# ==============================================================================

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

  # Other anatomical regions
  "14" = "#6B1E2D",
  "4"  = "#adb5bd",
  "10" = "#ff9505",
  "12" = "#e85d04",
  "5"  = "#9e0059",
  "15" = "#ffd100"
)

custom_cluster_labels <- c(
  "1"  = "Thalamus 1",
  "2"  = "Cortex 1",
  "3"  = "Hippocampus 1",
  "4"  = "Fiber tracts",
  "5"  = "Striatum-like amygdalar nuclei",
  "6"  = "Hypothalamus",
  "7"  = "Cortex 3",
  "8"  = "Cortex 2",
  "9"  = "Thalamus 2",
  "10" = "Cortical subplate",
  "11" = "Thalamus 3",
  "12" = "Caudoputamen",
  "13" = "Hippocampus 2",
  "14" = "Vessels",
  "15" = "Ventricles",
  "16" = "Hippocampus 3"
)


# ==============================================================================
# 3. UMAP configuration
# ==============================================================================

umap_input_reduction <- "integrated.cca"
umap_dims <- 1:20
umap_reduction_name <- "umap.ccaDims20"
umap_reduction_key <- "UMAPCCA_"
umap_n_neighbors <- 30L
umap_min_dist <- 0.30
umap_spread <- 1
umap_metric <- "cosine"
umap_seed <- 7L

# Existing UMAP is reused.
force_umap <- FALSE


# ==============================================================================
# 4. Paths
# ==============================================================================

spatial_functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "manual_cluster_annotation",
  "functions_spatialClusterVisualization_manualAnnotation.R"
)

umap_functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "manual_cluster_annotation",
  "functions_umapVisualization_manualAnnotation.R"
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

output_directory <- file.path(
  project_root,
  "results",
  dataset_name,
  "manual_cluster_annotation",
  configuration_name,
  "00_selectedClustering"
)

output_spatial_ids_pdf <- file.path(
  output_directory,
  paste0(
    "01_",
    dataset_name,
    "_",
    configuration_name,
    "_spatialClusters.pdf"
  )
)

output_umap_ids_pdf <- file.path(
  output_directory,
  paste0(
    "02_",
    dataset_name,
    "_",
    configuration_name,
    "_umapClusters.pdf"
  )
)

output_spatial_named_pdf <- file.path(
  output_directory,
  paste0(
    "03_",
    dataset_name,
    "_",
    configuration_name,
    "_spatialClusters_named.pdf"
  )
)

output_umap_named_pdf <- file.path(
  output_directory,
  paste0(
    "04_",
    dataset_name,
    "_",
    configuration_name,
    "_umapClusters_named.pdf"
  )
)


# ==============================================================================
# 5. Validate paths
# ==============================================================================

if (!dir.exists(file.path(project_root, "preprocessing")) ||
    !dir.exists(file.path(project_root, "results"))) {
  stop(
    "Run this script from the project root:\n",
    "/home/mateusz/projects/ippas-kunevicius-spatial",
    call. = FALSE
  )
}

required_input_files <- c(
  spatial_functions_file,
  umap_functions_file,
  input_rdata_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0L) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_input_files, collapse = "\n"),
    call. = FALSE
  )
}

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 6. Load packages and functions
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

source(spatial_functions_file)
source(umap_functions_file)


# ==============================================================================
# 7. Load one Seurat object
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
    "Objects loaded: ",
    paste(loaded_object_names, collapse = ", "),
    "\nSeurat objects found: ",
    paste(seurat_object_names, collapse = ", "),
    call. = FALSE
  )
}

seurat_object <- get(
  seurat_object_names[[1]],
  envir = load_environment,
  inherits = FALSE
)

message("Loaded Seurat object: ", seurat_object_names[[1]])
message(
  "Object dimensions: ",
  nrow(seurat_object),
  " genes x ",
  ncol(seurat_object),
  " spots"
)
message(
  "Available reductions: ",
  paste(
    SeuratObject::Reductions(seurat_object),
    collapse = ", "
  )
)


# ==============================================================================
# 8. Validate cluster IDs, colours and labels
# ==============================================================================

if (!cluster_column %in% colnames(seurat_object[[]])) {
  stop(
    "Missing clustering column: ",
    cluster_column,
    call. = FALSE
  )
}

cluster_ids <- unique(
  as.character(seurat_object[[]][[cluster_column]])
)

if (all(grepl("^[0-9]+$", cluster_ids))) {
  cluster_ids <- cluster_ids[
    order(as.integer(cluster_ids))
  ]
} else {
  cluster_ids <- sort(cluster_ids)
}

missing_colour_ids <- setdiff(
  cluster_ids,
  names(custom_cluster_colors)
)

if (length(missing_colour_ids) > 0L) {
  stop(
    "Missing colours for cluster ID(s): ",
    paste(missing_colour_ids, collapse = ", "),
    call. = FALSE
  )
}

missing_label_ids <- setdiff(
  cluster_ids,
  names(custom_cluster_labels)
)

if (length(missing_label_ids) > 0L) {
  stop(
    "Missing names for cluster ID(s): ",
    paste(missing_label_ids, collapse = ", "),
    call. = FALSE
  )
}


# ==============================================================================
# 9. OUTPUT 1/4: spatial plot with cluster IDs
# ==============================================================================

message("")
message("============================================================")
message("OUTPUT 1/4: SPATIAL — CLUSTER IDs")
message("============================================================")

save_selected_spatial_clustering_pdf(
  seurat_object = seurat_object,
  cluster_column = cluster_column,
  output_pdf = output_spatial_ids_pdf,
  sample_order = sample_order,
  ncol = 4L,
  plot_title = paste0(
    "Selected clustering for manual annotation: ",
    "CCA, Leiden resolution 0.40"
  ),
  plot_subtitle = paste0(
    "LogNormalize + VST | HVG 2000 | ",
    "dims 20 | k.param 20 | prune.SNN 0.0667"
  ),
  show_image = FALSE,
  image_alpha = 1,
  crop = FALSE,
  pt.size.factor = 2.1,
  legend_ncol = 4L,
  legend_point_size = 6,
  legend_height_ratio = 0.16,
  palette_name = "working30",
  custom_cluster_colors = custom_cluster_colors,
  pdf_width_in = 18,
  pdf_height_in = 22,
  verbose = TRUE
)


# ==============================================================================
# 10. OUTPUT 2/4: UMAP with cluster IDs
# ==============================================================================

message("")
message("============================================================")
message("OUTPUT 2/4: UMAP — CLUSTER IDs")
message("============================================================")

umap_ids_result <- save_selected_umap_clustering_pdf(
  seurat_object = seurat_object,
  cluster_column = cluster_column,
  output_pdf = output_umap_ids_pdf,

  input_reduction = umap_input_reduction,
  dims = umap_dims,
  umap_reduction_name = umap_reduction_name,
  umap_reduction_key = umap_reduction_key,
  n_neighbors = umap_n_neighbors,
  min_dist = umap_min_dist,
  spread = umap_spread,
  metric = umap_metric,
  seed_use = umap_seed,
  force_umap = force_umap,

  sample_order = sample_order,
  plot_title = paste0(
    "Selected clustering for manual annotation:\n",
    "CCA UMAP, Leiden resolution 0.40"
  ),
  plot_subtitle = paste0(
    "LogNormalize + VST | HVG 2000 | CCA dims 20 | ",
    "UMAP n.neighbors 30 | min.dist 0.30 | cosine"
  ),

  palette_name = "working30",
  custom_cluster_colors = custom_cluster_colors,
  custom_cluster_labels = NULL,
  include_cluster_id_in_label = TRUE,

  pt_size = 0.05,
  point_alpha = 1,
  shuffle = TRUE,
  shuffle_seed = umap_seed,
  raster = FALSE,
  raster_dpi = c(512, 512),

  legend_ncol = 3L,
  legend_point_size = 6,
  legend_height_ratio = 0.28,

  pdf_width_in = 12,
  pdf_height_in = 12,
  verbose = TRUE
)

seurat_object <- umap_ids_result$seurat_object


# ==============================================================================
# 11. Build temporary named cluster column for spatial plot
# ==============================================================================

named_spatial_column <- paste0(
  cluster_column,
  "__namedSpatial"
)

named_cluster_values_by_id <- setNames(
  object = paste0(
    cluster_ids,
    " — ",
    unname(custom_cluster_labels[cluster_ids])
  ),
  nm = cluster_ids
)

original_cluster_values <- as.character(
  seurat_object[[]][[cluster_column]]
)

named_cluster_values <- unname(
  named_cluster_values_by_id[original_cluster_values]
)

if (anyNA(named_cluster_values)) {
  stop(
    "Could not create the named spatial cluster column.",
    call. = FALSE
  )
}

seurat_object[[named_spatial_column]] <- factor(
  named_cluster_values,
  levels = unname(named_cluster_values_by_id[cluster_ids])
)

named_spatial_colors <- setNames(
  object = unname(custom_cluster_colors[cluster_ids]),
  nm = unname(named_cluster_values_by_id[cluster_ids])
)


# ==============================================================================
# 12. OUTPUT 3/4: spatial plot with cluster IDs + names
# ==============================================================================

message("")
message("============================================================")
message("OUTPUT 3/4: SPATIAL — CLUSTER IDs + NAMES")
message("============================================================")

save_selected_spatial_clustering_pdf(
  seurat_object = seurat_object,
  cluster_column = named_spatial_column,
  output_pdf = output_spatial_named_pdf,
  sample_order = sample_order,
  ncol = 4L,
  plot_title = paste0(
    "Selected clustering for manual annotation: ",
    "CCA, Leiden resolution 0.40 | named clusters"
  ),
  plot_subtitle = paste0(
    "LogNormalize + VST | HVG 2000 | ",
    "dims 20 | k.param 20 | prune.SNN 0.0667"
  ),
  show_image = FALSE,
  image_alpha = 1,
  crop = FALSE,
  pt.size.factor = 2.1,

  # Three columns for the spatial named plot.
  legend_ncol = 3L,
  legend_point_size = 6,
  legend_height_ratio = 0.28,

  palette_name = "working30",
  custom_cluster_colors = named_spatial_colors,

  pdf_width_in = 18,
  pdf_height_in = 22,
  verbose = TRUE
)


# ==============================================================================
# 13. OUTPUT 4/4: UMAP with cluster IDs + names
# ==============================================================================

message("")
message("============================================================")
message("OUTPUT 4/4: UMAP — CLUSTER IDs + NAMES")
message("============================================================")

umap_named_result <- save_selected_umap_clustering_pdf(
  seurat_object = seurat_object,
  cluster_column = cluster_column,
  output_pdf = output_umap_named_pdf,

  input_reduction = umap_input_reduction,
  dims = umap_dims,
  umap_reduction_name = umap_reduction_name,
  umap_reduction_key = umap_reduction_key,
  n_neighbors = umap_n_neighbors,
  min_dist = umap_min_dist,
  spread = umap_spread,
  metric = umap_metric,
  seed_use = umap_seed,

  # UMAP was already calculated or reused above.
  force_umap = FALSE,

  sample_order = sample_order,
  plot_title = paste0(
    "Selected clustering for manual annotation:\n",
    "CCA UMAP, Leiden resolution 0.40\n",
    "named clusters"
  ),
  plot_subtitle = paste0(
    "LogNormalize + VST | HVG 2000 | CCA dims 20 | ",
    "UMAP n.neighbors 30 | min.dist 0.30 | cosine"
  ),

  palette_name = "working30",
  custom_cluster_colors = custom_cluster_colors,
  custom_cluster_labels = custom_cluster_labels,
  include_cluster_id_in_label = TRUE,

  pt_size = 0.1125,
  point_alpha = 1,
  shuffle = TRUE,
  shuffle_seed = umap_seed,
  raster = FALSE,
  raster_dpi = c(512, 512),

  # Exactly two legend columns for the named UMAP.
  legend_ncol = 2L,
  legend_point_size = 6,
  legend_height_ratio = 0.42,

  pdf_width_in = 12,
  pdf_height_in = 12,
  verbose = TRUE
)


# ==============================================================================
# 14. Mandatory validation: all four files must exist
# ==============================================================================

expected_output_files <- c(
  spatial_ids = output_spatial_ids_pdf,
  umap_ids = output_umap_ids_pdf,
  spatial_named = output_spatial_named_pdf,
  umap_named = output_umap_named_pdf
)

output_validation <- vapply(
  expected_output_files,
  function(output_file) {
    file.exists(output_file) &&
      !is.na(file.info(output_file)$size) &&
      file.info(output_file)$size > 0L
  },
  logical(1)
)

if (!all(output_validation)) {
  failed_outputs <- expected_output_files[
    !output_validation
  ]

  stop(
    "The script did not generate all four required PDF files.\n",
    "Missing or empty output(s):\n",
    paste(failed_outputs, collapse = "\n"),
    call. = FALSE
  )
}

message("")
message("============================================================")
message("SUCCESS: ALL FOUR PDF FILES WERE GENERATED")
message("============================================================")

for (output_name in names(expected_output_files)) {
  message(
    output_name,
    ":\n",
    expected_output_files[[output_name]]
  )
}

message("")
message("Finished script version: ", script_version)
