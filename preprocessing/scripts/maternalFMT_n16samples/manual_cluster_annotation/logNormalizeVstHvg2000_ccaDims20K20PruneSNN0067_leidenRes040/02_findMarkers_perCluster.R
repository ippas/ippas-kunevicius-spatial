# ==============================================================================
# 02_findMarkers_perCluster.R
#
# Purpose:
# Find positive marker genes for each Leiden cluster and save:
# 1. marker results,
# 2. number of markers per cluster.
# ==============================================================================


# ==============================================================================
# 1. Configuration
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
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

min_pct <- 0.25
logfc_threshold <- 0.25
adjusted_p_value_threshold <- 0.05


# ==============================================================================
# 2. Paths
# ==============================================================================

input_file <- file.path(
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
  "01_findMarkers_perCluster"
)

markers_output_file <- file.path(
  output_directory,
  paste0(
    "01_",
    dataset_name,
    "_leidenRes040_",
    "findAllMarkers_minPct025_log2FC025_padj005.tsv"
  )
)

summary_output_file <- file.path(
  output_directory,
  paste0(
    "02_",
    dataset_name,
    "_leidenRes040_",
    "numberOfMarkersPerCluster.tsv"
  )
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 3. Load Seurat object
# ==============================================================================

load_environment <- new.env()

loaded_object_names <- load(
  input_file,
  envir = load_environment
)

seurat_object_names <- loaded_object_names[
  vapply(
    loaded_object_names,
    function(object_name) {
      inherits(
        get(object_name, envir = load_environment),
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
  envir = load_environment
)

rm(load_environment)


# ==============================================================================
# 4. Prepare object
# ==============================================================================

DefaultAssay(seurat_object) <- "RNA"

rna_data_layers <- grep(
  pattern = "^data",
  x = Layers(seurat_object[["RNA"]]),
  value = TRUE
)

if (length(rna_data_layers) > 1L) {
  seurat_object <- JoinLayers(
    object = seurat_object,
    assay = "RNA"
  )
}

Idents(seurat_object) <- cluster_column

cluster_levels <- as.character(
  sort(
    unique(
      as.numeric(
        as.character(Idents(seurat_object))
      )
    )
  )
)


# ==============================================================================
# 5. Find markers
# ==============================================================================

markers <- FindAllMarkers(
  object = seurat_object,
  assay = "RNA",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold,
  return.thresh = adjusted_p_value_threshold,
  verbose = TRUE
)

markers$cluster <- as.character(
  markers$cluster
)

markers <- markers[
  order(
    as.numeric(markers$cluster),
    -markers$avg_log2FC
  ),
]


# ==============================================================================
# 6. Save marker results
# ==============================================================================

write.table(
  markers,
  file = markers_output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# ==============================================================================
# 7. Count and save markers per cluster
# ==============================================================================

marker_counts <- table(
  factor(
    markers$cluster,
    levels = cluster_levels
  )
)

marker_summary <- data.frame(
  cluster = names(marker_counts),
  number_of_markers = as.integer(marker_counts)
)

write.table(
  marker_summary,
  file = summary_output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# ==============================================================================
# 8. Summary
# ==============================================================================

cat("\nNumber of markers per cluster:\n\n")

print(
  marker_summary,
  row.names = FALSE
)

cat(
  "\nTotal number of markers: ",
  nrow(markers),
  "\n",
  sep = ""
)

cat(
  "\nMarker results:\n",
  markers_output_file,
  "\n",
  sep = ""
)

cat(
  "\nMarker count summary:\n",
  summary_output_file,
  "\n",
  sep = ""
)