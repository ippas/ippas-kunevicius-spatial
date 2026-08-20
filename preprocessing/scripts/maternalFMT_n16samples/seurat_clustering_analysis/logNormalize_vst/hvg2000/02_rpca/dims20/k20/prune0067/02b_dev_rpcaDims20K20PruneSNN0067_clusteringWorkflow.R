# ==============================================================================
# 02b_dev_rpcaDims20K20PruneSNN0067_clusteringWorkflow.R
#
# Development workflow:
# - input: existing RPCA neighbours/SNN graph Seurat object;
# - algorithms: Louvain, Louvain refined, SLM and Leiden;
# - resolutions: 0.20, 0.30 and 0.40 for every algorithm;
# - total clustering runs: 4 algorithms x 3 resolutions = 12;
# - every run is performed directly on the existing SNN graph;
# - the Seurat object is not modified during individual clustering runs;
# - all 12 cluster columns are added once, only after every run finishes and
#   the complete assignment table passes validation;
# - the input 01_*_neighborsGraph.RData file is never overwritten;
# - the clustered Seurat object is saved as a new 02b_dev_* RData file.
#
# Intended execution:
# Run the complete script from Positron.
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
  "leidenbase"
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

graph_configuration_tag <-
  "rpcaDims20K20PruneSNN0067"

snn_graph_name <- paste0(
  graph_configuration_tag,
  "_snn"
)

clustering_algorithms <- c(
  "louvain",
  "louvainRefined",
  "slm",
  "leiden"
)

# The same resolution values are used for every clustering algorithm.
clustering_resolutions <- c(
  0.20,
  0.30,
  0.40
)

modularity_fxn <- 1L
n_start <- 10L
n_iter <- 10L
random_seed <- 7L
group_singletons <- TRUE

leiden_method <- "leidenbase"
leiden_objective_function <- "modularity"

set_active_column <- NULL

overwrite_existing_cluster_columns <- FALSE


# ==============================================================================
# 3. Paths
# ==============================================================================

project_root <- normalizePath(
  "/home/mateusz/projects/ippas-kunevicius-spatial",
  winslash = "/",
  mustWork = TRUE
)

functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "seurat_clustering_analysis",
  "logNormalize_vst",
  "functions_multiClustering.R"
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

tables_directory <- file.path(
  analysis_root,
  "tables"
)

dir.create(
  rdata_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

input_object_name <- paste0(
  dataset_name,
  "_logNormalizeVst_hvg2000_",
  graph_configuration_tag,
  "_neighborsGraph"
)

output_object_name <- paste0(
  dataset_name,
  "_logNormalizeVst_hvg2000_",
  graph_configuration_tag,
  "_multiClustering"
)

input_rdata_file <- file.path(
  rdata_directory,
  paste0(
    "01_",
    input_object_name,
    ".RData"
  )
)

output_prefix <- paste0(
  "02b_dev_",
  dataset_name,
  "_logNormalizeVst_hvg2000_",
  graph_configuration_tag,
  "_multiClustering"
)

output_rdata_file <- file.path(
  rdata_directory,
  paste0(
    output_prefix,
    ".RData"
  )
)

cluster_assignments_file <- file.path(
  tables_directory,
  paste0(
    output_prefix,
    "_clusterAssignments.tsv.gz"
  )
)

clustering_summary_file <- file.path(
  tables_directory,
  paste0(
    output_prefix,
    "_clusteringSummary.tsv"
  )
)

if (!file.exists(functions_file)) {
  stop(
    "Functions file does not exist:\n",
    functions_file
  )
}

if (!file.exists(input_rdata_file)) {
  stop(
    "Input RData file does not exist:\n",
    input_rdata_file
  )
}

source(functions_file)


# ==============================================================================
# 4. Load exactly one expected Seurat object
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
message("SNN graph: ", snn_graph_name)
message(
  "Algorithms: ",
  paste(
    clustering_algorithms,
    collapse = ", "
  )
)
message(
  "Resolutions for every algorithm: ",
  paste(
    formatC(
      clustering_resolutions,
      format = "f",
      digits = 2L
    ),
    collapse = ", "
  )
)
message(
  "Total clustering runs: ",
  length(clustering_algorithms) *
    length(clustering_resolutions)
)


# ==============================================================================
# 5. Run all algorithms at resolutions 0.20, 0.30 and 0.40
#
# The Seurat object is modified only once, after all 12 clustering runs have
# finished and the complete assignment table has passed validation.
# ==============================================================================

clustering_result <- run_multi_clustering(
  seurat_object = seurat_object,
  graph_name = snn_graph_name,
  algorithms = clustering_algorithms,
  resolutions = clustering_resolutions,
  resolution_digits = 2L,
  modularity_fxn = modularity_fxn,
  n_start = n_start,
  n_iter = n_iter,
  random_seed = random_seed,
  group_singletons = group_singletons,
  leiden_method = leiden_method,
  leiden_objective_function =
    leiden_objective_function,
  overwrite_existing =
    overwrite_existing_cluster_columns,
  set_active_column =
    set_active_column,
  verbose = TRUE,
  verbose_findclusters = TRUE
)

clustered_seurat_object <-
  clustering_result$seurat_object


# ==============================================================================
# 6. Verify all newly added metadata columns
# ==============================================================================

resolution_labels <- vapply(
  clustering_resolutions,
  format_clustering_resolution,
  character(1),
  digits = 2L
)

expected_cluster_columns <- unlist(
  lapply(
    clustering_algorithms,
    function(algorithm_name) {
      paste0(
        algorithm_name,
        "_res",
        resolution_labels
      )
    }
  ),
  use.names = FALSE
)

missing_cluster_columns <- setdiff(
  expected_cluster_columns,
  colnames(
    clustered_seurat_object[[]]
  )
)

if (length(missing_cluster_columns) > 0L) {
  stop(
    "Expected cluster columns were not added:\n",
    paste(
      missing_cluster_columns,
      collapse = ", "
    )
  )
}

for (cluster_column in expected_cluster_columns) {

  cluster_values <-
    clustered_seurat_object[[]][
      ,
      cluster_column
    ]

  if (
    length(cluster_values) !=
      ncol(clustered_seurat_object)
  ) {
    stop(
      "Cluster column has an invalid length: ",
      cluster_column
    )
  }

  if (anyNA(cluster_values)) {
    stop(
      "Cluster column contains NA values: ",
      cluster_column
    )
  }
}


# ==============================================================================
# 7. Save external clustering tables
# ==============================================================================

cluster_assignments_connection <- gzfile(
  cluster_assignments_file,
  open = "wt"
)

utils::write.table(
  x = clustering_result$cluster_assignments,
  file = cluster_assignments_connection,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

close(
  cluster_assignments_connection
)

utils::write.table(
  x = clustering_result$clustering_summary,
  file = clustering_summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)


# ==============================================================================
# 8. Save the clustered Seurat object as a new RData file
#
# The original 01_*_neighborsGraph.RData file is never overwritten.
# ==============================================================================

assign(
  output_object_name,
  clustered_seurat_object
)

save(
  list = output_object_name,
  file = output_rdata_file,
  compress = TRUE
)

if (
  !file.exists(output_rdata_file) ||
    is.na(file.info(output_rdata_file)$size) ||
    file.info(output_rdata_file)$size <= 0L
) {
  stop(
    "The output RData file was not created correctly:\n",
    output_rdata_file
  )
}


# ==============================================================================
# 9. Final report
# ==============================================================================

message("")
message("============================================================")
message("DEVELOPMENT MULTI-CLUSTERING WORKFLOW COMPLETED")
message("============================================================")

message("")
message("Input RData was not modified:")
message(input_rdata_file)

message("")
message("Output Seurat object:")
message(output_object_name)

message("")
message("Output RData:")
message(output_rdata_file)

message("")
message("Cluster assignments:")
message(cluster_assignments_file)

message("")
message("Clustering summary:")
message(clustering_summary_file)

message("")
message("Added metadata columns:")
message(
  paste(
    expected_cluster_columns,
    collapse = ", "
  )
)

message("")
message("Clustering summary:")

print(
  clustering_result$clustering_summary[
    ,
    c(
      "algorithm",
      "resolution",
      "clusterColumn",
      "nSpots",
      "nClusters",
      "elapsedFormatted"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)
