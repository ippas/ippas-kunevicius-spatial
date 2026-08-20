# ==============================================================================
# 03_calculateGeneSpecificity_onFindMarkers.R
#
# Purpose:
# Calculate cluster-specificity metrics for genes identified by FindAllMarkers().
# ==============================================================================


# ==============================================================================
# 1. Configuration
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
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


# ==============================================================================
# 2. Paths
# ==============================================================================

functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "gene_specificity",
  "functions_geneSpecificity_onFindMarkers.R"
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

input_markers_file <- file.path(
  project_root,
  "results",
  dataset_name,
  "manual_cluster_annotation",
  configuration_name,
  "01_findMarkers_perCluster",
  paste0(
    "01_",
    dataset_name,
    "_leidenRes040_",
    "findAllMarkers_minPct025_log2FC025_padj005.tsv"
  )
)

output_directory <- file.path(
  project_root,
  "results",
  dataset_name,
  "manual_cluster_annotation",
  configuration_name,
  "02_geneSpecificity_onFindMarkers"
)

output_file <- file.path(
  output_directory,
  paste0(
    "01_",
    dataset_name,
    "_leidenRes040_",
    "findAllMarkers_withGeneSpecificity.tsv"
  )
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

source(functions_file)


# ==============================================================================
# 3. Load Seurat object
# ==============================================================================

load_environment <- new.env()

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
          envir = load_environment
        ),
        "Seurat"
      )
    },
    logical(1)
  )
]

if (length(seurat_object_names) != 1L) {
  stop(
    "Expected exactly one Seurat object.",
    call. = FALSE
  )
}

seurat_object <- get(
  seurat_object_names[[1]],
  envir = load_environment
)

rm(load_environment)


# ==============================================================================
# 4. Load FindAllMarkers results
# ==============================================================================

markers <- read.delim(
  input_markers_file,
  check.names = FALSE
)


# ==============================================================================
# 5. Calculate gene specificity
# ==============================================================================

markers_with_specificity <-
  calculate_gene_specificity_on_findmarkers(
    seurat_object = seurat_object,
    marker_table = markers,
    cluster_column = cluster_column,
    assay_name = "RNA"
  )


# ==============================================================================
# 6. Save results
# ==============================================================================

write.table(
  markers_with_specificity,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# ==============================================================================
# 7. Summary
# ==============================================================================

new_columns <- setdiff(
  colnames(markers_with_specificity),
  colnames(markers)
)

cat(
  "\nMarkers analysed: ",
  nrow(markers_with_specificity),
  "\n",
  sep = ""
)

cat(
  "New columns: ",
  length(new_columns),
  "\n",
  sep = ""
)

cat(
  paste0(
    "  - ",
    new_columns,
    collapse = "\n"
  ),
  "\n"
)

cat(
  "\nSaved to:\n",
  output_file,
  "\n",
  sep = ""
)