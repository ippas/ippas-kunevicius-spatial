# ==============================================================================
# 00_prepare_logNormalizeVst_hvg2000.R
#
# Purpose:
# Prepare the maternalFMT_n16samples Seurat v5 parent object for downstream
# integration and clustering analyses using:
#   LogNormalize -> VST (2000 HVGs) -> ScaleData -> PCA (50 PCs)
#
# Samples excluded from the original 20-sample metadata:
#   20_1F, 12_3F, 15_1M, 20_3M
#
# Run from the project root, for example:
#   Rscript preprocessing/scripts/maternalFMT_n16samples/\
#     seurat_clustering_analysis/logNormalize_vst/hvg2000/00_preprocessing/\
#     00_prepare_logNormalizeVst_hvg2000.R
#
# Alternatively, set the project root explicitly:
#   IPPAS_SPATIAL_PROJECT_ROOT=/path/to/ippas-kunevicius-spatial Rscript ...
# ==============================================================================  


main <- function() {

  # ============================================================================
  # 1. Project paths and fixed analysis parameters
  # ============================================================================

  project_root <- Sys.getenv(
    "IPPAS_SPATIAL_PROJECT_ROOT",
    unset = getwd()
  )
  project_root <- normalizePath(project_root, mustWork = TRUE)

  path_to_data <- file.path(
    project_root,
    "data",
    "spacerangerCount_withoutJSON_GRCm39-2024-A_2026-03-28"
  )

  metadata_file <- file.path(
    project_root,
    "data",
    "metadata_autismFMT.tsv"
  )

  functions_file <- file.path(
    project_root,
    "preprocessing",
    "src",
    "seurat_clustering_analysis",
    "logNormalize_vst",
    "functions_logNormalizeVst_preprocessing.R"
  )

  output_root <- file.path(
    project_root,
    "results",
    "maternalFMT_n16samples",
    "seurat_clustering_analysis",
    "logNormalize_vst",
    "hvg2000",
    "00_preprocessing"
  )

  dataset_name <- "maternalFMT_n16samples"
  project_name <- "MaternalFMT_n16samples"

  excluded_sample_ids <- c(
    "20_1F",
    "12_3F",
    "15_1M",
    "20_3M"
  )

  expected_n_samples <- 16L
  sample_id_col <- "sample_ID"

  nfeatures <- 2000L
  scale_factor <- 10000
  npcs <- 50L
  seed <- 7L

  mitochondrial_pattern <- "^mt-"
  min_cells <- 0L
  min_features <- 0L
  load_images_when_possible <- TRUE

  file_prefix <- paste0(
    dataset_name,
    "_logNormalizeVst_hvg", nfeatures,
    "_pca", npcs
  )


  # ============================================================================
  # 2. Create output folders and start logging
  # ============================================================================

  output_paths <- stats::setNames(
    file.path(
      output_root,
      c("RData", "tables", "plots", "config", "logs")
    ),
    c("RData", "tables", "plots", "config", "logs")
  )

  invisible(
    lapply(
      output_paths,
      dir.create,
      recursive = TRUE,
      showWarnings = FALSE
    )
  )

  log_file <- file.path(
    output_paths[["logs"]],
    paste0(file_prefix, "_preprocessingRun.log")
  )

  log_connection <- file(log_file, open = "wt")
  sink(log_connection, type = "output", split = TRUE)
  sink(log_connection, type = "message")

  on.exit(
    {
      while (sink.number(type = "message") > 2) {
        sink(type = "message")
      }
      while (sink.number(type = "output") > 0) {
        sink(type = "output")
      }
      close(log_connection)
    },
    add = TRUE
  )

  message("Analysis started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
  message("Project root: ", project_root)
  message("Output root: ", output_root)


  # ============================================================================
  # 3. Load packages and custom functions
  # ============================================================================

  required_packages <- c(
    "Seurat",
    "SeuratObject",
    "Matrix",
    "ggplot2"
  )

  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
    library(ggplot2)
  })

  if (!file.exists(functions_file)) {
    stop("Functions file does not exist: ", functions_file, call. = FALSE)
  }

  source(functions_file)

  check_required_packages(required_packages)
  assert_seurat_v5()
  output_paths <- create_output_directories(output_root)

  set.seed(seed)


  # ============================================================================
  # 4. Read and filter metadata
  # ============================================================================

  metadata_n16 <- read_and_filter_sample_metadata(
    metadata_file = metadata_file,
    excluded_sample_ids = excluded_sample_ids,
    sample_id_col = sample_id_col,
    expected_n_samples = expected_n_samples
  )

  message(
    "Samples retained: ",
    paste(metadata_n16[[sample_id_col]], collapse = ", ")
  )

  write_tsv(
    metadata_n16,
    file.path(
      output_paths[["tables"]],
      paste0(file_prefix, "_sampleMetadata.tsv")
    )
  )


  # ============================================================================
  # 5. Read Space Ranger samples, coordinates and QC metrics
  # ============================================================================

  spatial_input <- read_spatial_samples_from_metadata(
    metadata = metadata_n16,
    path_to_data = path_to_data,
    sample_id_col = sample_id_col,
    mitochondrial_pattern = mitochondrial_pattern,
    min_cells = min_cells,
    min_features = min_features,
    load_images_when_possible = load_images_when_possible,
    verbose = TRUE
  )

  samples_list <- spatial_input$samples_list
  input_manifest <- spatial_input$input_manifest

  if (length(samples_list) != expected_n_samples) {
    stop(
      "Expected ", expected_n_samples, " loaded samples, but found ",
      length(samples_list), ".",
      call. = FALSE
    )
  }

  qc_summary <- summarize_qc_by_sample(samples_list)

  write_tsv(
    input_manifest,
    file.path(
      output_paths[["tables"]],
      paste0(file_prefix, "_inputManifest.tsv")
    )
  )

  write_tsv(
    qc_summary,
    file.path(
      output_paths[["tables"]],
      paste0(file_prefix, "_qcSummaryBySample.tsv")
    )
  )


  # ============================================================================
  # 6. Prepare the LogNormalize + VST parent object
  # ============================================================================

  preprocessing_result <- prepare_logNormalize_vst_parent_object(
    samples_list = samples_list,
    nfeatures = nfeatures,
    scale_factor = scale_factor,
    npcs = npcs,
    project = project_name,
    sample_id_col = sample_id_col,
    seed = seed,
    verbose = TRUE
  )

  maternalFMT_n16samples_logNormalizeVst_hvg2000_pca50_parentObject <-
    preprocessing_result$object


  # ============================================================================
  # 7. Save the parent object
  # ============================================================================

  rdata_file <- file.path(
    output_paths[["RData"]],
    paste0("01_", file_prefix, "_parentObject.RData")
  )

  save(
    maternalFMT_n16samples_logNormalizeVst_hvg2000_pca50_parentObject,
    file = rdata_file,
    compress = TRUE
  )

  message("Saved parent object: ", rdata_file)


  # ============================================================================
  # 8. Save tables and configuration
  # ============================================================================

  variable_features_table <- create_variable_features_table(
    maternalFMT_n16samples_logNormalizeVst_hvg2000_pca50_parentObject
  )

  pca_standard_deviations <- create_pca_standard_deviation_table(
    maternalFMT_n16samples_logNormalizeVst_hvg2000_pca50_parentObject,
    reduction = "pca"
  )

  preprocessing_summary <- create_preprocessing_summary(
    object = maternalFMT_n16samples_logNormalizeVst_hvg2000_pca50_parentObject,
    n_samples = expected_n_samples,
    nfeatures_requested = nfeatures,
    scale_factor = scale_factor,
    npcs = npcs,
    seed = seed
  )

  run_parameters <- create_run_parameters_table(
    path_to_data = path_to_data,
    metadata_file = metadata_file,
    excluded_sample_ids = excluded_sample_ids,
    output_root = output_root,
    nfeatures = nfeatures,
    scale_factor = scale_factor,
    npcs = npcs,
    seed = seed,
    mitochondrial_pattern = mitochondrial_pattern,
    min_cells = min_cells,
    min_features = min_features,
    load_images_when_possible = load_images_when_possible
  )

  write_tsv(
    variable_features_table,
    file.path(
      output_paths[["tables"]],
      paste0(file_prefix, "_variableFeatures.tsv")
    )
  )

  write_tsv(
    pca_standard_deviations,
    file.path(
      output_paths[["tables"]],
      paste0(file_prefix, "_pcaStandardDeviations.tsv")
    )
  )

  write_tsv(
    preprocessing_summary,
    file.path(
      output_paths[["tables"]],
      paste0(file_prefix, "_preprocessingSummary.tsv")
    )
  )

  write_tsv(
    run_parameters,
    file.path(
      output_paths[["config"]],
      paste0(file_prefix, "_runParameters.tsv")
    )
  )

  session_info_file <- file.path(
    output_paths[["config"]],
    paste0(file_prefix, "_sessionInfo.txt")
  )

  utils::capture.output(
    utils::sessionInfo(),
    file = session_info_file
  )


  # ============================================================================
  # 9. Save diagnostic plots
  # ============================================================================

  save_preprocessing_plots(
    object = maternalFMT_n16samples_logNormalizeVst_hvg2000_pca50_parentObject,
    output_directory = output_paths[["plots"]],
    file_prefix = file_prefix,
    npcs = npcs
  )


  # ============================================================================
  # 10. Final checks
  # ============================================================================

  final_object <- maternalFMT_n16samples_logNormalizeVst_hvg2000_pca50_parentObject

  if (length(unique(final_object$sample_ID)) != expected_n_samples) {
    stop("Final object does not contain exactly 16 sample IDs.", call. = FALSE)
  }

  if (length(SeuratObject::VariableFeatures(final_object[["RNA"]])) == 0) {
    stop("Final object has no variable features.", call. = FALSE)
  }

  if (!"pca" %in% names(final_object@reductions)) {
    stop("Final object does not contain the PCA reduction.", call. = FALSE)
  }

  message("Final object: ", nrow(final_object), " genes x ", ncol(final_object), " spots.")
  message(
    "Variable features: ",
    length(SeuratObject::VariableFeatures(final_object[["RNA"]]))
  )
  message("PCA components: ", ncol(SeuratObject::Embeddings(final_object, "pca")))
  message("Analysis completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))

  invisible(final_object)
}


main()
