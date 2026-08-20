# ==============================================================================
# 02d_dev_rpcaDims20K20PruneSNN0067_clusteringValidation.R
#
# Validation of 12 clustering solutions:
# 4 algorithms × 3 resolutions.
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})


# ==============================================================================
# 1. Paths
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"

analysis_root <- file.path(
  project_root,
  "results",
  "maternalFMT_n16samples",
  "seurat_clustering_analysis",
  "logNormalize_vst",
  "hvg2000",
  "02_rpca",
  "dims20",
  "k20",
  "prune0067"
)

source(
  file.path(
    project_root,
    "preprocessing",
    "src",
    "seurat_clustering_analysis",
    "logNormalize_vst",
    "functions_clusteringValidation.R"
  )
)


# ==============================================================================
# 2. Load the object containing all 12 clustering solutions
# ==============================================================================

input_object_name <- paste0(
  "maternalFMT_n16samples_",
  "logNormalizeVst_hvg2000_",
  "rpcaDims20K20PruneSNN0067_",
  "multiClustering"
)

load(
  file.path(
    analysis_root,
    "RData",
    paste0(
      "02b_dev_",
      input_object_name,
      ".RData"
    )
  )
)

seurat_object <- get(
  input_object_name
)


# ==============================================================================
# 3. Clustering solutions
# ==============================================================================

algorithms <- c(
  "louvain",
  "louvainRefined",
  "slm",
  "leiden"
)

resolutions <- c(
  0.20,
  0.30,
  0.40
)

cluster_columns <- unlist(
  lapply(
    algorithms,
    function(algorithm) {
      paste0(
        algorithm,
        "_res",
        formatC(
          resolutions,
          format = "f",
          digits = 2
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
# 4. Run validation
# ==============================================================================

validation_results <- run_clustering_validation(
  seurat_object = seurat_object,
  cluster_columns = cluster_columns,
  reduction = "integrated.rpca",
  dims = 1:20,
  sample_order = sample_order,

  output_dir = file.path(
    analysis_root,
    "clustering_validation"
  ),

  analysis_prefix =
    "logNormalizeVst_hvg2000_rpcaDims20K20PruneSNN0067",

  integration_method = "RPCA",
  normalization_label = "LogNormalize + VST",
  n_hvg = 2000L,
  k_param = 20L,
  prune_snn = 1 / 15,

  transcriptomic_silhouette_max_spots = 4000L,
  spatial_silhouette_max_spots_per_sample = 3000L,

  pas_k = 10L,
  pas_minimum_different_neighbours = 6L,

  seed = 7L,
  save_png = TRUE,
  save_pdf = TRUE,
  verbose = TRUE
)
