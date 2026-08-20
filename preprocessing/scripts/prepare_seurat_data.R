# ==============================================================================
# 01_read_risperidone_data.R
#
# Purpose:
# Read risperidone 3q29 metadata and verify Spatial / Space Ranger input folders.
# Compatible with Seurat v5.
# ==============================================================================


# ==============================================================================
# 1. Load packages
# ==============================================================================

suppressPackageStartupMessages({

  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(stringr)

  library(ggplot2)

  library(Seurat)
  library(SeuratObject)

  library(future)
  library(future.apply)

  library(data.table)
  library(matrixStats)
  library(moments)

})


# ==============================================================================
# 2. Check Seurat version
# ==============================================================================

message("Seurat version: ", as.character(packageVersion("Seurat")))
message("SeuratObject version: ", as.character(packageVersion("SeuratObject")))

if (packageVersion("Seurat") < "5.0.0") {
  stop(
    "This workflow should be run with Seurat v5 or newer. ",
    "Current version: ", packageVersion("Seurat")
  )
}


# ==============================================================================
# 3. Load custom functions
# ==============================================================================
source("preprocessing/src/add_clusters_multi_resolution.R")
source("preprocessing/src/add_spatial_images_to_samples.R")
source("preprocessing/src/functions_prepare_seurat_data.R")
source("preprocessing/src/integrate_seurat_samples_rpca.R")


# ==============================================================================
# 4. Global settings
# ==============================================================================

path_to_data <- "data/spacerangerCount_withoutJSON_GRCm39-2024-A_2026-03-28/"

metadata_file <- "data/metadata_autismFMT.tsv"

annotation_file <- paste0(
  "data/risperidone-3q29/",
  "gene-annotation/peaks-annotate-reduction.tsv"
)



# ==============================================================================
# 5. Read and inspect metadata
# ==============================================================================

metadata_autismFMT <- read.delim(
  metadata_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

samples_list <- read_spatial_samples(
  path_to_data = path_to_data,
  metadata = metadata_autismFMT,
  sample_id_col = "sample_ID",
  min.cells = 0,
  min.features = 0,
  verbose = TRUE
)

samples_list$`23_1M` %>% str

samples_list <- add_spatial_coordinates_to_samples(
  samples_list = samples_list,
  path_to_data = path_to_data,
  verbose = TRUE
)

samples_list[["23_1M"]][[]] |>
  dplyr::select(
    sample_ID,
    in_tissue,
    array_row,
    array_col,
    pxl_row_in_fullres,
    pxl_col_in_fullres
  ) |>
  head()

qc_results <- add_qc_metrics_to_samples(
  samples_list = samples_list,
  mitochondrial_pattern = "^mt-",
  verbose = TRUE
)

samples_list <- qc_results$samples_list
qc_summary <- qc_results$qc_summary

print(qc_summary)

write.table(
  qc_summary,
  file = file.path(
    "results/maternalFMT_n20samples",
    "maternalFMT_n20samples_QCsummaryBySample.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

image_results <- add_spatial_images_to_samples(
  samples_list = samples_list,
  path_to_data = path_to_data,
  raw_image_root = "raw",
  cache_dir = "results/spatial_image_cache",
  assay = "RNA",
  verbose = TRUE
)

samples_list <- image_results$samples_list
spatial_image_summary <- image_results$image_summary


maternalFMT_integrated_n20samples <- integrate_seurat_samples_rpca(
  samples_list = samples_list,
  nfeatures = 2000,
  dims = 1:30,
  project = "MaternalFMT_n20samples",
  seed = 7,
  verbose = TRUE
)

save(
  maternalFMT_integrated_n20samples,
  file = "results/maternalFMT_n20samples/maternalFMT_integrated_n20samples.RData"
)

# load("results/maternalFMT_n20samples/maternalFMT_integrated_n20samples.RData")

maternalFMT_integrated_n20samples <- add_clusters_multi_resolution(
  seurat_object = maternalFMT_integrated_n20samples,
  reduction = "integrated.rpca",
  dims = 1:30,
  resolution_start = 0.05,
  resolution_end = 2,
  resolution_step = 0.05,
  graph_names = c(
    "integrated_rpca_nn",
    "integrated_rpca_snn"
  ),
  cluster_prefix = "clusters_res",
  k.param = 20,
  algorithm = 1,
  random.seed = 7,
  verbose = TRUE
)

save(
  maternalFMT_integrated_n20samples,
  file = "results/maternalFMT_n20samples/maternalFMT_integrated_n20samples_withClusters.RData"
)

load("results/maternalFMT_n20samples/maternalFMT_integrated_n20samples_withClusters.RData")

maternalFMT_integrated_n20samples <- add_clusters_multi_resolution(
  seurat_object = maternalFMT_integrated_n20samples,
  reduction = "integrated.rpca",
  dims = 1:30,
  resolution_start = 0.2,
  resolution_end = 0.30,
  resolution_step = 0.01,
  graph_names = c(
    "integrated_rpca_nn",
    "integrated_rpca_snn"
  ),
  cluster_prefix = "clusters_res",
  k.param = 20,
  algorithm = 1,
  random.seed = 7,
  force_neighbors = FALSE,
  force_clusters = FALSE,
  verbose = TRUE
)

maternalFMT_integrated_n20samples[[]]$clusters_res0.1 %>% 
  as.character() %>% unique() %>% length()
