# ==============================================================================
# 02_run_rpcaPruneSNN0067_dimsKGrid_completeWorkflow.R
#
# Full terminal workflow for the maternalFMT_n16samples RPCA grid.
#
# Grid:
# - RPCA dimensions: 20, 30, 40 and 50;
# - k.param: 20, 30, 40 and 50;
# - prune.SNN: 0.0667 (1/15);
# - clustering algorithms:
#     Louvain, Louvain refined, SLM and Leiden;
# - clustering resolutions:
#     0.05-2.00 with a step of 0.05.
#
# Complete workflow:
# 1. Load the LogNormalize + VST / HVG2000 / PCA50 parent object.
# 2. Run RPCA integration once for each dims value and save a reusable cache.
# 3. Build and save every 01_*_neighborsGraph.RData object for k=20/30/40/50.
# 4. Run 160 clustering solutions per graph: 4 algorithms x 40 resolutions.
# 5. Calculate UMAP and save one 02_*_multiClusteringAndUmap.RData object.
# 6. Save external clustering-assignment and clustering-summary tables.
# 7. Generate spatial and UMAP plots for all clustering solutions.
# 8. Run transcriptomic ASW, spatial ASW, CHAOS, PAS, ARI stability and
#    cluster-size validation.
# 9. Save grid-wide summary tables across all 16 graph configurations.
#
# Three stage logs are written in every configuration-specific `logs/` folder:
# - stage 01: RPCA cache, neighbors graph, clustering and UMAP;
# - stage 02: spatial and UMAP visualizations;
# - stage 03: clustering validation.
#
# Recommended location:
# preprocessing/scripts/maternalFMT_n16samples/seurat_clustering_analysis/
# logNormalize_vst/hvg2000/02_rpca/
#
# Intended terminal execution from the project root:
# Rscript preprocessing/scripts/maternalFMT_n16samples/\
# seurat_clustering_analysis/logNormalize_vst/hvg2000/02_rpca/\
# 02_run_rpcaPruneSNN0067_dimsKGrid_completeWorkflow.R
# ==============================================================================


# ==============================================================================
# 0. Utilities
# ==============================================================================

options(
  stringsAsFactors = FALSE
)


workflow_message <- function(...) {

  cat(
    "[",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    "] ",
    paste0(..., collapse = ""),
    "\n",
    sep = ""
  )

  flush.console()
  invisible(NULL)
}


format_elapsed_time <- function(seconds) {

  seconds <- max(
    0,
    as.numeric(seconds)
  )

  sprintf(
    "%02d:%02d:%02d",
    as.integer(
      floor(seconds / 3600)
    ),
    as.integer(
      floor(
        (seconds %% 3600) / 60
      )
    ),
    as.integer(
      floor(seconds %% 60)
    )
  )
}


write_tsv <- function(
    data,
    filename
) {

  dir.create(
    dirname(filename),
    recursive = TRUE,
    showWarnings = FALSE
  )

  utils::write.table(
    x = data,
    file = filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )

  invisible(filename)
}


write_gzipped_tsv <- function(
    data,
    filename
) {

  dir.create(
    dirname(filename),
    recursive = TRUE,
    showWarnings = FALSE
  )

  output_connection <- gzfile(
    filename,
    open = "wt"
  )

  on.exit(
    close(output_connection),
    add = TRUE
  )

  utils::write.table(
    x = data,
    file = output_connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )

  invisible(filename)
}


run_logged_stage <- function(
    log_file,
    stage_label,
    code
) {

  dir.create(
    dirname(log_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  log_connection <- file(
    log_file,
    open = "wt"
  )

  sink(
    log_connection,
    type = "output",
    split = TRUE
  )

  on.exit(
    {
      while (
        sink.number(type = "output") >
          0L
      ) {
        sink(
          type = "output"
        )
      }

      close(log_connection)
    },
    add = TRUE
  )

  stage_started <- Sys.time()

  workflow_message(
    "============================================================"
  )

  workflow_message(
    "START STAGE: ",
    stage_label
  )

  workflow_message(
    "Log file: ",
    log_file
  )

  stage_result <- withCallingHandlers(
    force(code),
    message = function(message_condition) {

      workflow_message(
        "MESSAGE: ",
        conditionMessage(
          message_condition
        )
      )

      invokeRestart(
        "muffleMessage"
      )
    },
    warning = function(warning_condition) {

      workflow_message(
        "WARNING: ",
        conditionMessage(
          warning_condition
        )
      )

      invokeRestart(
        "muffleWarning"
      )
    }
  )

  stage_finished <- Sys.time()

  workflow_message(
    "FINISH STAGE: ",
    stage_label,
    " | elapsed=",
    format_elapsed_time(
      as.numeric(
        difftime(
          stage_finished,
          stage_started,
          units = "secs"
        )
      )
    )
  )

  workflow_message(
    "============================================================"
  )

  stage_result
}


write_skipped_stage_log <- function(
    log_file,
    stage_label,
    reason
) {

  dir.create(
    dirname(log_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  writeLines(
    c(
      paste0(
        "Stage: ",
        stage_label
      ),
      paste0(
        "Status: SKIPPED"
      ),
      paste0(
        "Reason: ",
        reason
      ),
      paste0(
        "Time: ",
        format(
          Sys.time(),
          "%Y-%m-%d %H:%M:%S %z"
        )
      )
    ),
    con = log_file
  )

  invisible(log_file)
}


load_expected_seurat_object <- function(
    rdata_file,
    expected_object_name
) {

  load_environment <- new.env(
    parent = globalenv()
  )

  loaded_object_names <- load(
    rdata_file,
    envir = load_environment
  )

  if (
    !expected_object_name %in%
      loaded_object_names
  ) {
    stop(
      "Expected object `",
      expected_object_name,
      "` was not found in:\n",
      rdata_file,
      "\nObjects found: ",
      paste(
        loaded_object_names,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  output_object <- get(
    expected_object_name,
    envir = load_environment,
    inherits = FALSE
  )

  if (!inherits(output_object, "Seurat")) {
    stop(
      "Object `",
      expected_object_name,
      "` is not a Seurat object.",
      call. = FALSE
    )
  }

  output_object
}


save_named_object <- function(
    object,
    object_name,
    rdata_file
) {

  dir.create(
    dirname(rdata_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  save_environment <- new.env(
    parent = emptyenv()
  )

  assign(
    object_name,
    object,
    envir = save_environment
  )

  save(
    list = object_name,
    file = rdata_file,
    envir = save_environment,
    compress = TRUE
  )

  if (
    !file.exists(rdata_file) ||
      is.na(file.info(rdata_file)$size) ||
      file.info(rdata_file)$size <= 0L
  ) {
    stop(
      "RData file was not saved correctly:\n",
      rdata_file,
      call. = FALSE
    )
  }

  invisible(rdata_file)
}


all_required_files_exist <- function(
    files
) {

  files <- files[
    !is.na(files) &
      files != ""
  ]

  length(files) > 0L &&
    all(
      file.exists(files)
    )
}


format_resolution <- function(
    resolution
) {

  formatC(
    resolution,
    format = "f",
    digits = 2L
  )
}


build_cluster_columns <- function(
    algorithms,
    resolutions
) {

  resolution_labels <- vapply(
    resolutions,
    format_resolution,
    character(1)
  )

  unlist(
    lapply(
      algorithms,
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
}


validate_cluster_columns <- function(
    seurat_object,
    expected_cluster_columns
) {

  metadata_table <- seurat_object[[]]

  missing_cluster_columns <- setdiff(
    expected_cluster_columns,
    colnames(metadata_table)
  )

  if (length(missing_cluster_columns) > 0L) {
    stop(
      "Expected clustering columns are missing:\n",
      paste(
        missing_cluster_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  for (
    cluster_column in
    expected_cluster_columns
  ) {

    cluster_values <-
      metadata_table[[cluster_column]]

    if (
      length(cluster_values) !=
        ncol(seurat_object)
    ) {
      stop(
        "Invalid clustering-column length: ",
        cluster_column,
        call. = FALSE
      )
    }

    if (anyNA(cluster_values)) {
      stop(
        "NA values found in clustering column: ",
        cluster_column,
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


get_dynamic_legend_layout <- function(
    n_clusters
) {

  spatial_ncol <- if (
    n_clusters <= 24L
  ) {
    4L
  } else if (
    n_clusters <= 40L
  ) {
    5L
  } else {
    6L
  }

  umap_ncol <- if (
    n_clusters <= 18L
  ) {
    3L
  } else if (
    n_clusters <= 36L
  ) {
    4L
  } else {
    5L
  }

  spatial_rows <- ceiling(
    n_clusters /
      spatial_ncol
  )

  umap_rows <- ceiling(
    n_clusters /
      umap_ncol
  )

  list(
    spatialNcol = spatial_ncol,
    spatialHeightRatio = max(
      0.16,
      min(
        0.65,
        0.045 * spatial_rows
      )
    ),
    spatialPlotHeight = 22 +
      max(
        0,
        spatial_rows - 5
      ) * 0.55,
    umapNcol = umap_ncol,
    umapHeightRatio = max(
      0.32,
      min(
        0.85,
        0.065 * umap_rows
      )
    ),
    umapPlotHeight = 12 +
      max(
        0,
        umap_rows - 5
      ) * 0.50
  )
}


prepend_configuration_columns <- function(
    data,
    dims_value,
    k_value,
    graph_configuration_tag
) {

  data.frame(
    rpcaDims = as.integer(
      dims_value
    ),
    graphK = as.integer(
      k_value
    ),
    graphConfiguration =
      graph_configuration_tag,
    data,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


# ==============================================================================
# 1. Main workflow
# ==============================================================================

main <- function() {

  analysis_started <- Sys.time()


  # ============================================================================
  # 1.1 Configuration
  # ============================================================================

  dataset_name <- "maternalFMT_n16samples"

  normalization_tag <- "logNormalizeVst"
  normalization_label <-
    "LogNormalize + VST"

  n_hvg <- 2000L

  integration_method <- "RPCA"
  integration_reduction <-
    "integrated.rpca"

  rpca_dims_values <- c(
    20L,
    30L,
    40L,
    50L
  )

  graph_k_values <- c(
    20L,
    30L,
    40L,
    50L
  )

  prune_snn <- 1 / 15
  prune_directory <- "prune0067"
  prune_tag <- "PruneSNN0067"

  clustering_algorithms <- c(
    "louvain",
    "louvainRefined",
    "slm",
    "leiden"
  )

  clustering_resolutions <- round(
    seq(
      from = 0.05,
      to = 2.00,
      by = 0.05
    ),
    digits = 2L
  )

  modularity_fxn <- 1L
  n_start <- 10L
  n_iter <- 10L
  random_seed <- 7L
  group_singletons <- TRUE

  leiden_method <- "leidenbase"
  leiden_objective_function <-
    "modularity"

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


  # ============================================================================
  # 1.2 UMAP configuration
  # ============================================================================

  umap_n_neighbors <- 30L
  umap_min_dist <- 0.30
  umap_spread <- 1
  umap_metric <- "cosine"
  umap_seed <- 7L

  # Recalculate only when the expected reduction is absent.
  force_umap <- FALSE


  # ============================================================================
  # 1.3 Visualization configuration
  # ============================================================================

  cluster_palette_name <- "working30"

  show_image <- FALSE
  image_alpha <- 1
  crop <- FALSE
  spatial_point_size_factor <- 1.8
  spatial_panel_columns <- 4L

  umap_point_size <- 0.05
  umap_point_alpha <- 1
  umap_raster <- FALSE
  umap_raster_dpi <- c(
    512,
    512
  )

  save_spatial_png <- TRUE
  save_spatial_pdf <- TRUE
  save_umap_png <- TRUE
  save_umap_pdf <- TRUE
  save_validation_png <- TRUE
  save_validation_pdf <- TRUE

  figure_dpi <- 300L


  # ============================================================================
  # 1.4 Validation configuration
  # ============================================================================

  transcriptomic_silhouette_max_spots <-
    4000L

  transcriptomic_silhouette_minimum_spots_per_cluster <-
    20L

  spatial_silhouette_max_spots_per_sample <-
    3000L

  spatial_silhouette_minimum_spots_per_cluster <-
    10L

  pas_k <- 10L

  pas_minimum_different_neighbours <-
    6L

  normalize_chaos_by_spot_pitch <-
    TRUE

  # For 160 solutions, this avoids the unnecessary all-by-all comparison of
  # different algorithms at different resolutions. It retains:
  # - all resolution comparisons within each algorithm;
  # - comparisons between algorithms at the same resolution.
  ari_comparison_mode <-
    "withinAlgorithmAndSameResolution"


  # ============================================================================
  # 1.5 Resume and error handling
  # ============================================================================

  resume_completed_stages <- TRUE
  overwrite_existing_figures <- FALSE
  overwrite_existing_validation <- FALSE

  continue_after_configuration_error <-
    TRUE

  fail_at_end_if_any_stage_failed <-
    TRUE


  # ============================================================================
  # 2. Project paths and function files
  # ============================================================================

  project_root_candidate <- Sys.getenv(
    "IPPAS_SPATIAL_PROJECT_ROOT",
    unset = Sys.getenv(
      "SPATIAL_ASD_PROJECT_ROOT",
      unset = getwd()
    )
  )

  project_root <- normalizePath(
    project_root_candidate,
    winslash = "/",
    mustWork = TRUE
  )

  expected_project_directories <- c(
    file.path(
      project_root,
      "preprocessing"
    ),
    file.path(
      project_root,
      "results"
    )
  )

  if (
    !all(
      dir.exists(
        expected_project_directories
      )
    )
  ) {
    stop(
      "Run the script from the project root or define ",
      "`SPATIAL_ASD_PROJECT_ROOT`.\nResolved project root: ",
      project_root,
      call. = FALSE
    )
  }

  functions_directory <- file.path(
    project_root,
    "preprocessing",
    "src",
    "seurat_clustering_analysis",
    "logNormalize_vst"
  )

  functions_files <- c(
    file.path(
      functions_directory,
      "functions_integration_neighborsGraph.R"
    ),
    file.path(
      functions_directory,
      "functions_multiClustering.R"
    ),
    file.path(
      functions_directory,
      "functions_spatialClusterVisualization.R"
    ),
    file.path(
      functions_directory,
      "functions_umapVisualization.R"
    ),
    file.path(
      functions_directory,
      "functions_clusteringValidation.R"
    )
  )

  missing_functions_files <-
    functions_files[
      !file.exists(
        functions_files
      )
    ]

  if (
    length(missing_functions_files) >
      0L
  ) {
    stop(
      "Missing function file(s):\n",
      paste(
        missing_functions_files,
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  parent_rdata_file <- file.path(
    project_root,
    "results",
    dataset_name,
    "seurat_clustering_analysis",
    "logNormalize_vst",
    paste0(
      "hvg",
      n_hvg
    ),
    "00_preprocessing",
    "RData",
    paste0(
      "01_",
      dataset_name,
      "_",
      normalization_tag,
      "_hvg",
      n_hvg,
      "_pca50_parentObject.RData"
    )
  )

  if (!file.exists(parent_rdata_file)) {
    stop(
      "PCA50 parent RData does not exist:\n",
      parent_rdata_file,
      call. = FALSE
    )
  }

  results_rpca_root <- file.path(
    project_root,
    "results",
    dataset_name,
    "seurat_clustering_analysis",
    "logNormalize_vst",
    paste0(
      "hvg",
      n_hvg
    ),
    "02_rpca"
  )

  grid_summary_root <- file.path(
    results_rpca_root,
    "grid_summary",
    "prune0067"
  )

  grid_tables_directory <- file.path(
    grid_summary_root,
    "tables"
  )

  dir.create(
    grid_tables_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  grid_prefix <-
    "02_rpcaPruneSNN0067_dims20to50_k20to50"


  # ============================================================================
  # 3. Packages and functions
  # ============================================================================

  required_packages <- c(
    "Seurat",
    "SeuratObject",
    "leidenbase",
    "ggplot2",
    "patchwork",
    "dplyr",
    "uwot",
    "cluster",
    "RANN"
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
      ".",
      call. = FALSE
    )
  }

  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
  })

  for (functions_file in functions_files) {
    source(functions_file)
  }

  if (
    utils::packageVersion("Seurat") <
      "5.0.0"
  ) {
    stop(
      "Seurat v5 or newer is required.",
      call. = FALSE
    )
  }

  set.seed(
    random_seed
  )


  # ============================================================================
  # 4. Configuration grid and manifest
  # ============================================================================

  configuration_grid <- expand.grid(
    rpcaDims = rpca_dims_values,
    graphK = graph_k_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  configuration_grid <- configuration_grid[
    order(
      configuration_grid$rpcaDims,
      configuration_grid$graphK
    ),
    ,
    drop = FALSE
  ]

  rownames(configuration_grid) <- NULL

  configuration_grid$graphConfiguration <-
    paste0(
      "rpcaDims",
      configuration_grid$rpcaDims,
      "K",
      configuration_grid$graphK,
      prune_tag
    )

  configuration_table_file <- file.path(
    grid_tables_directory,
    paste0(
      grid_prefix,
      "_configurationGrid.tsv"
    )
  )

  write_tsv(
    configuration_grid,
    configuration_table_file
  )

  manifest_file <- file.path(
    grid_tables_directory,
    paste0(
      grid_prefix,
      "_runManifest.tsv"
    )
  )

  manifest_rows <- list()
  combined_clustering_summaries <- list()
  combined_validation_summaries <- list()

  expected_cluster_columns <-
    build_cluster_columns(
      algorithms =
        clustering_algorithms,
      resolutions =
        clustering_resolutions
    )

  workflow_message(
    "Project root: ",
    project_root
  )

  workflow_message(
    "PCA50 parent RData: ",
    parent_rdata_file
  )

  workflow_message(
    "Configurations: ",
    nrow(configuration_grid)
  )

  workflow_message(
    "Clustering solutions per configuration: ",
    length(expected_cluster_columns)
  )

  workflow_message(
    "Total clustering runs: ",
    nrow(configuration_grid) *
      length(expected_cluster_columns)
  )

  workflow_message(
    "Grid configuration table: ",
    configuration_table_file
  )


  # ============================================================================
  # 5. Process every dims x k graph configuration
  # ============================================================================

  for (
    configuration_index in
    seq_len(nrow(configuration_grid))
  ) {

    configuration_started <- Sys.time()

    dims_value <-
      configuration_grid$rpcaDims[[configuration_index]]

    k_value <-
      configuration_grid$graphK[[configuration_index]]

    graph_configuration_tag <-
      configuration_grid$graphConfiguration[[configuration_index]]

    snn_graph_name <- paste0(
      graph_configuration_tag,
      "_snn"
    )

    analysis_prefix <- paste0(
      normalization_tag,
      "_hvg",
      n_hvg,
      "_",
      graph_configuration_tag
    )

    analysis_root <- file.path(
      results_rpca_root,
      paste0(
        "dims",
        dims_value
      ),
      paste0(
        "k",
        k_value
      ),
      prune_directory
    )

    rdata_directory <- file.path(
      analysis_root,
      "RData"
    )

    tables_directory <- file.path(
      analysis_root,
      "tables"
    )

    logs_directory <- file.path(
      project_root,
      "preprocessing",
      "logs",
      dataset_name,
      "seurat_clustering_analysis",
      "logNormalize_vst",
      paste0(
        "hvg",
        n_hvg
      ),
      "02_rpca",
      paste0(
        "dims",
        dims_value
      ),
      paste0(
        "k",
        k_value
      ),
      prune_directory
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

    validation_directory <- file.path(
      analysis_root,
      "clustering_validation"
    )

    for (
      current_directory in
      c(
        rdata_directory,
        tables_directory,
        logs_directory,
        spatial_figures_root,
        umap_figures_root,
        validation_directory
      )
    ) {
      dir.create(
        current_directory,
        recursive = TRUE,
        showWarnings = FALSE
      )
    }

    integrated_cache_directory <- file.path(
      results_rpca_root,
      paste0(
        "dims",
        dims_value
      ),
      "RData"
    )

    dir.create(
      integrated_cache_directory,
      recursive = TRUE,
      showWarnings = FALSE
    )

    integrated_cache_object_name <- paste0(
      dataset_name,
      "_",
      normalization_tag,
      "_hvg",
      n_hvg,
      "_rpcaDims",
      dims_value,
      "_integrated"
    )

    integrated_cache_rdata_file <- file.path(
      integrated_cache_directory,
      paste0(
        "01_",
        integrated_cache_object_name,
        ".RData"
      )
    )

    graph_names <- make_graph_names(
      graph_configuration_tag
    )

    input_object_name <- paste0(
      dataset_name,
      "_",
      normalization_tag,
      "_hvg",
      n_hvg,
      "_",
      graph_configuration_tag,
      "_neighborsGraph"
    )

    output_object_name <- paste0(
      dataset_name,
      "_",
      normalization_tag,
      "_hvg",
      n_hvg,
      "_",
      graph_configuration_tag,
      "_multiClusteringAndUmap"
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
      "02_",
      output_object_name
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

    umap_reduction_name <- paste0(
      "umap.rpcaDims",
      dims_value
    )

    umap_reduction_key <- "UMAP_"

    stage_01_log <- file.path(
      logs_directory,
      paste0(
        output_prefix,
        "_stage01_integrationNeighborsClusteringAndUmap.log"
      )
    )

    stage_02_log <- file.path(
      logs_directory,
      paste0(
        output_prefix,
        "_stage02_spatialAndUmapVisualization.log"
      )
    )

    stage_03_log <- file.path(
      logs_directory,
      paste0(
        output_prefix,
        "_stage03_clusteringValidation.log"
      )
    )

    workflow_message("")
    workflow_message(
      "############################################################"
    )

    workflow_message(
      "CONFIGURATION ",
      configuration_index,
      "/",
      nrow(configuration_grid),
      ": ",
      graph_configuration_tag
    )

    workflow_message(
      "Analysis root: ",
      analysis_root
    )

    workflow_message(
      "############################################################"
    )

    stage_01_status <- "not_run"
    stage_02_status <- "not_run"
    stage_03_status <- "not_run"

    stage_01_elapsed <- NA_real_
    stage_02_elapsed <- NA_real_
    stage_03_elapsed <- NA_real_

    configuration_error <- NA_character_
    seurat_object <- NULL


    # ==========================================================================
    # STAGE 01: integration cache + neighbors graph + clustering + UMAP
    # ==========================================================================


      stage_01_started <- Sys.time()

      stage_01_result <- tryCatch(
        {
          run_logged_stage(
            log_file = stage_01_log,
            stage_label =
              "01 integration, neighbors graph, clustering and UMAP",
            code = {

              workflow_message(
                "Configuration: ",
                graph_configuration_tag
              )

              workflow_message(
                "Input RData: ",
                input_rdata_file
              )

              workflow_message(
                "Output RData: ",
                output_rdata_file
              )

              workflow_message(
                "SNN graph: ",
                snn_graph_name
              )

              workflow_message(
                "Algorithms: ",
                paste(
                  clustering_algorithms,
                  collapse = ", "
                )
              )

              workflow_message(
                "Resolutions: 0.05-2.00 by 0.05"
              )

              workflow_message(
                "Expected clustering columns: ",
                length(
                  expected_cluster_columns
                )
              )

              existing_stage_outputs_complete <-
                all_required_files_exist(
                  c(
                    input_rdata_file,
                    output_rdata_file,
                    cluster_assignments_file,
                    clustering_summary_file
                  )
                )

              if (
                isTRUE(
                  resume_completed_stages
                ) &&
                  existing_stage_outputs_complete
              ) {

                workflow_message(
                  "Existing stage-01 outputs detected."
                )

                existing_object <-
                  load_expected_seurat_object(
                    rdata_file =
                      output_rdata_file,
                    expected_object_name =
                      output_object_name
                  )

                validate_cluster_columns(
                  seurat_object =
                    existing_object,
                  expected_cluster_columns =
                    expected_cluster_columns
                )

                umap_present <-
                  umap_reduction_name %in%
                  SeuratObject::Reductions(
                    existing_object
                  )

                if (!umap_present) {

                  workflow_message(
                    "Expected UMAP is absent. Calculating it and updating ",
                    "the existing 02 RData file."
                  )

                  umap_result <- run_umap_if_missing(
                    seurat_object =
                      existing_object,
                    input_reduction =
                      integration_reduction,
                    dims =
                      seq_len(dims_value),
                    umap_reduction_name =
                      umap_reduction_name,
                    umap_reduction_key =
                      umap_reduction_key,
                    n_neighbors =
                      umap_n_neighbors,
                    min_dist =
                      umap_min_dist,
                    spread =
                      umap_spread,
                    metric =
                      umap_metric,
                    seed_use =
                      umap_seed,
                    force_umap =
                      FALSE,
                    verbose = TRUE
                  )

                  existing_object <-
                    umap_result$seurat_object

                  rm(umap_result)

                  save_named_object(
                    object =
                      existing_object,
                    object_name =
                      output_object_name,
                    rdata_file =
                      output_rdata_file
                  )

                  stage_status <-
                    "updated_missing_umap"

                } else {

                  workflow_message(
                    "Expected UMAP is already present: ",
                    umap_reduction_name
                  )

                  stage_status <-
                    "skipped_complete"
                }

                list(
                  seuratObject =
                    existing_object,
                  stageStatus =
                    stage_status
                )

              } else {

                if (!file.exists(input_rdata_file)) {

                  workflow_message(
                    "The configuration-specific neighbors graph is absent and will be created."
                  )

                  if (file.exists(integrated_cache_rdata_file)) {

                    workflow_message(
                      "Loading reusable RPCA integration cache: ",
                      integrated_cache_rdata_file
                    )

                    integrated_object <-
                      load_expected_seurat_object(
                        rdata_file =
                          integrated_cache_rdata_file,
                        expected_object_name =
                          integrated_cache_object_name
                      )

                    validate_integrated_object(
                      object = integrated_object,
                      dims = seq_len(dims_value),
                      integrated_reduction =
                        integration_reduction
                    )

                  } else {

                    workflow_message(
                      "RPCA integration cache is absent. Loading PCA50 parent object: ",
                      parent_rdata_file
                    )

                    parent_object <-
                      load_single_seurat_object(
                        parent_rdata_file
                      )

                    workflow_message(
                      "Running reusable RPCA integration for dimensions 1:",
                      dims_value,
                      "."
                    )

                    integrated_object <-
                      run_integration_only(
                        object = parent_object,
                        integration_method = "rpca",
                        dims = seq_len(dims_value),
                        assay = "RNA",
                        original_reduction = "pca",
                        integrated_reduction =
                          integration_reduction,
                        features =
                          SeuratObject::VariableFeatures(
                            parent_object
                          ),
                        sample_id_col = "sample_ID",
                        expected_n_samples = 16L,
                        normalization_method =
                          "LogNormalize",
                        integration_method_args = list(),
                        seed = random_seed,
                        verbose = TRUE
                      )

                    save_named_object(
                      object = integrated_object,
                      object_name =
                        integrated_cache_object_name,
                      rdata_file =
                        integrated_cache_rdata_file
                    )

                    workflow_message(
                      "Saved reusable RPCA integration cache: ",
                      integrated_cache_rdata_file
                    )

                    rm(parent_object)
                    invisible(gc(verbose = FALSE))
                  }

                  workflow_message(
                    "Building configuration-specific kNN/SNN graphs: ",
                    paste(graph_names, collapse = ", ")
                  )

                  neighbors_object <-
                    build_neighbors_graph_from_integrated(
                      object = integrated_object,
                      dims = seq_len(dims_value),
                      k_param = k_value,
                      prune_snn = prune_snn,
                      integrated_reduction =
                        integration_reduction,
                      graph_names = graph_names,
                      nn_method = "annoy",
                      n_trees = 50L,
                      distance_metric = "euclidean",
                      l2_norm = FALSE,
                      seed = random_seed,
                      overwrite_graphs = FALSE,
                      verbose = TRUE
                    )

                  save_named_object(
                    object = neighbors_object,
                    object_name = input_object_name,
                    rdata_file = input_rdata_file
                  )

                  workflow_message(
                    "Saved 01 neighbors-graph RData: ",
                    input_rdata_file
                  )

                  rm(integrated_object, neighbors_object)
                  invisible(gc(verbose = FALSE))

                } else {

                  workflow_message(
                    "Existing 01 neighbors-graph RData will be reused: ",
                    input_rdata_file
                  )
                }

                workflow_message(
                  "Loading input neighbors/SNN object."
                )

                input_object <-
                  load_expected_seurat_object(
                    rdata_file =
                      input_rdata_file,
                    expected_object_name =
                      input_object_name
                  )

                workflow_message(
                  "Input object: ",
                  nrow(input_object),
                  " genes x ",
                  ncol(input_object),
                  " spots."
                )

                if (
                  !snn_graph_name %in%
                    SeuratObject::Graphs(
                      input_object
                    )
                ) {
                  stop(
                    "SNN graph `",
                    snn_graph_name,
                    "` is absent from the input object.",
                    call. = FALSE
                  )
                }

                clustering_result <-
                  run_multi_clustering(
                    seurat_object =
                      input_object,
                    graph_name =
                      snn_graph_name,
                    algorithms =
                      clustering_algorithms,
                    resolutions =
                      clustering_resolutions,
                    resolution_digits = 2L,
                    modularity_fxn =
                      modularity_fxn,
                    n_start =
                      n_start,
                    n_iter =
                      n_iter,
                    random_seed =
                      random_seed,
                    group_singletons =
                      group_singletons,
                    leiden_method =
                      leiden_method,
                    leiden_objective_function =
                      leiden_objective_function,
                    overwrite_existing =
                      FALSE,
                    set_active_column =
                      NULL,
                    verbose = TRUE,
                    verbose_findclusters =
                      TRUE
                  )

                clustered_object <-
                  clustering_result$seurat_object

                validate_cluster_columns(
                  seurat_object =
                    clustered_object,
                  expected_cluster_columns =
                    expected_cluster_columns
                )

                workflow_message(
                  "Saving clustering tables."
                )

                write_gzipped_tsv(
                  data =
                    clustering_result$cluster_assignments,
                  filename =
                    cluster_assignments_file
                )

                write_tsv(
                  data =
                    clustering_result$clustering_summary,
                  filename =
                    clustering_summary_file
                )

                workflow_message(
                  "Calculating UMAP if absent: ",
                  umap_reduction_name
                )

                umap_result <- run_umap_if_missing(
                  seurat_object =
                    clustered_object,
                  input_reduction =
                    integration_reduction,
                  dims =
                    seq_len(dims_value),
                  umap_reduction_name =
                    umap_reduction_name,
                  umap_reduction_key =
                    umap_reduction_key,
                  n_neighbors =
                    umap_n_neighbors,
                  min_dist =
                    umap_min_dist,
                  spread =
                    umap_spread,
                  metric =
                    umap_metric,
                  seed_use =
                    umap_seed,
                  force_umap =
                    force_umap,
                  verbose = TRUE
                )

                clustered_object <-
                  umap_result$seurat_object

                rm(
                  input_object,
                  umap_result
                )

                workflow_message(
                  "Saving 02 RData containing clustering metadata and UMAP."
                )

                save_named_object(
                  object =
                    clustered_object,
                  object_name =
                    output_object_name,
                  rdata_file =
                    output_rdata_file
                )

                rm(
                  clustering_result
                )

                invisible(
                  gc(verbose = FALSE)
                )

                list(
                  seuratObject =
                    clustered_object,
                  stageStatus =
                    "completed"
                )
              }
            }
          )
        },
        error = function(error_condition) {

          list(
            error =
              conditionMessage(
                error_condition
              )
          )
        }
      )

      stage_01_elapsed <- as.numeric(
        difftime(
          Sys.time(),
          stage_01_started,
          units = "secs"
        )
      )

      if (
        !is.null(stage_01_result$error)
      ) {

        stage_01_status <- paste0(
          "error: ",
          stage_01_result$error
        )

        configuration_error <-
          stage_01_result$error

        workflow_message(
          "STAGE 01 ERROR: ",
          stage_01_result$error
        )

      } else {

        stage_01_status <-
          stage_01_result$stageStatus

        seurat_object <-
          stage_01_result$seuratObject
      }


      # ========================================================================
      # STAGE 02: spatial + UMAP plots
      # ========================================================================

      if (is.null(seurat_object)) {

        stage_02_status <-
          "skipped_stage01_failed"

        write_skipped_stage_log(
          stage_02_log,
          "02 spatial and UMAP visualization",
          "Stage 01 did not produce a valid Seurat object."
        )

      } else {

        stage_02_started <- Sys.time()

        stage_02_result <- tryCatch(
          {
            run_logged_stage(
              log_file = stage_02_log,
              stage_label =
                "02 spatial and UMAP visualization",
              code = {

                integration_parameter_line <- paste0(
                  "Integration: ",
                  integration_method,
                  " | normalization: ",
                  normalization_label,
                  " | HVG: ",
                  n_hvg,
                  " | dims: 1–",
                  dims_value,
                  " | k.param: ",
                  k_value,
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
                  " | dims: 1–",
                  dims_value,
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

                plots_created <- 0L
                plots_skipped <- 0L

                for (
                  clustering_index in
                  seq_along(
                    expected_cluster_columns
                  )
                ) {

                  cluster_column <-
                    expected_cluster_columns[[clustering_index]]

                  n_clusters <- length(
                    unique(
                      seurat_object[[]][[cluster_column]]
                    )
                  )

                  legend_layout <-
                    get_dynamic_legend_layout(
                      n_clusters
                    )

                  workflow_message(
                    "[",
                    clustering_index,
                    "/",
                    length(
                      expected_cluster_columns
                    ),
                    "] ",
                    cluster_column,
                    " | clusters=",
                    n_clusters
                  )

                  spatial_output_dir <- file.path(
                    spatial_figures_root,
                    cluster_column
                  )

                  umap_output_dir <- file.path(
                    umap_figures_root,
                    cluster_column
                  )

                  spatial_output_prefix <-
                    build_spatial_cluster_output_prefix(
                      analysis_prefix =
                        analysis_prefix,
                      cluster_column =
                        cluster_column
                    )

                  umap_output_prefix <-
                    build_umap_cluster_output_prefix(
                      analysis_prefix =
                        analysis_prefix,
                      cluster_column =
                        cluster_column
                    )

                  spatial_expected_files <- c(
                    file.path(
                      spatial_output_dir,
                      paste0(
                        spatial_output_prefix,
                        "_clusterSummary.tsv"
                      )
                    ),
                    file.path(
                      spatial_output_dir,
                      paste0(
                        spatial_output_prefix,
                        "_imageMap.tsv"
                      )
                    ),
                    file.path(
                      spatial_output_dir,
                      paste0(
                        spatial_output_prefix,
                        "_clusterCountsBySample.tsv"
                      )
                    ),
                    if (
                      isTRUE(
                        save_spatial_png
                      )
                    ) {
                      file.path(
                        spatial_output_dir,
                        paste0(
                          spatial_output_prefix,
                          "_spatialClusters.png"
                        )
                      )
                    } else {
                      NA_character_
                    },
                    if (
                      isTRUE(
                        save_spatial_pdf
                      )
                    ) {
                      file.path(
                        spatial_output_dir,
                        paste0(
                          spatial_output_prefix,
                          "_spatialClusters.pdf"
                        )
                      )
                    } else {
                      NA_character_
                    }
                  )

                  umap_expected_files <- c(
                    if (
                      isTRUE(save_umap_png)
                    ) {
                      file.path(
                        umap_output_dir,
                        paste0(
                          umap_output_prefix,
                          "_umapClusters.png"
                        )
                      )
                    } else {
                      NA_character_
                    },
                    if (
                      isTRUE(save_umap_pdf)
                    ) {
                      file.path(
                        umap_output_dir,
                        paste0(
                          umap_output_prefix,
                          "_umapClusters.pdf"
                        )
                      )
                    } else {
                      NA_character_
                    }
                  )

                  spatial_complete <-
                    all_required_files_exist(
                      spatial_expected_files
                    )

                  umap_complete <-
                    all_required_files_exist(
                      umap_expected_files
                    )

                  if (
                    isTRUE(
                      resume_completed_stages
                    ) &&
                      !isTRUE(
                        overwrite_existing_figures
                      ) &&
                      spatial_complete &&
                      umap_complete
                  ) {

                    workflow_message(
                      "SKIP existing spatial and UMAP outputs: ",
                      cluster_column
                    )

                    plots_skipped <-
                      plots_skipped + 1L

                    next
                  }

                  if (
                    !spatial_complete ||
                      isTRUE(
                        overwrite_existing_figures
                      )
                  ) {

                    spatial_plot_title <- paste0(
                      "Spatial clustering across 16 maternal FMT samples: ",
                      cluster_column
                    )

                    spatial_plot_subtitle <- paste(
                      integration_parameter_line,
                      paste0(
                        "4 × 4 panel layout | legend: global n, %, mean±SD ",
                        "across samples | palette: ",
                        cluster_palette_name,
                        " | show_image = ",
                        if (show_image) {
                          "TRUE"
                        } else {
                          "FALSE"
                        }
                      ),
                      sep = "\n"
                    )

                    spatial_plot_result <-
                      plot_spatial_clusters_all_samples(
                        seurat_object =
                          seurat_object,
                        cluster_column =
                          cluster_column,
                        sample_order =
                          sample_order,
                        ncol =
                          spatial_panel_columns,
                        plot_title =
                          spatial_plot_title,
                        plot_subtitle =
                          spatial_plot_subtitle,
                        include_cluster_count_in_title =
                          TRUE,
                        show_image =
                          show_image,
                        image_alpha =
                          image_alpha,
                        crop =
                          crop,
                        pt.size.factor =
                          spatial_point_size_factor,
                        legend_position =
                          "top",
                        legend_ncol =
                          legend_layout$spatialNcol,
                        legend_point_size =
                          6,
                        legend_height_ratio =
                          legend_layout$spatialHeightRatio,
                        palette_name =
                          cluster_palette_name,
                        output_dir =
                          spatial_output_dir,
                        output_prefix =
                          spatial_output_prefix,
                        save_png =
                          save_spatial_png,
                        save_pdf =
                          save_spatial_pdf,
                        png_width_in =
                          18,
                        png_height_in =
                          legend_layout$spatialPlotHeight,
                        pdf_width_in =
                          18,
                        pdf_height_in =
                          legend_layout$spatialPlotHeight,
                        dpi =
                          figure_dpi,
                        verbose =
                          TRUE
                      )

                    rm(
                      spatial_plot_result
                    )

                    invisible(
                      gc(verbose = FALSE)
                    )
                  }

                  if (
                    !umap_complete ||
                      isTRUE(
                        overwrite_existing_figures
                      )
                  ) {

                    umap_plot_title <- paste0(
                      "UMAP clustering across 16 maternal FMT samples: ",
                      cluster_column
                    )

                    umap_plot_subtitle <- paste(
                      umap_parameter_line,
                      paste0(
                        "Spots: ",
                        format(
                          ncol(seurat_object),
                          big.mark = " ",
                          scientific = FALSE
                        ),
                        " | legend: global n, %, mean±SD across samples | ",
                        "palette: ",
                        cluster_palette_name
                      ),
                      sep = "\n"
                    )

                    umap_plot_result <-
                      plot_umap_clusters(
                        seurat_object =
                          seurat_object,
                        cluster_column =
                          cluster_column,
                        umap_reduction_name =
                          umap_reduction_name,
                        sample_order =
                          sample_order,
                        plot_title =
                          umap_plot_title,
                        plot_subtitle =
                          umap_plot_subtitle,
                        include_cluster_count_in_title =
                          TRUE,
                        palette_name =
                          cluster_palette_name,
                        pt_size =
                          umap_point_size,
                        point_alpha =
                          umap_point_alpha,
                        shuffle =
                          TRUE,
                        shuffle_seed =
                          umap_seed,
                        raster =
                          umap_raster,
                        raster_dpi =
                          umap_raster_dpi,
                        legend_ncol =
                          legend_layout$umapNcol,
                        legend_point_size =
                          6,
                        legend_height_ratio =
                          legend_layout$umapHeightRatio,
                        output_dir =
                          umap_output_dir,
                        output_prefix =
                          umap_output_prefix,
                        save_png =
                          save_umap_png,
                        save_pdf =
                          save_umap_pdf,
                        png_width_in =
                          12,
                        png_height_in =
                          legend_layout$umapPlotHeight,
                        pdf_width_in =
                          12,
                        pdf_height_in =
                          legend_layout$umapPlotHeight,
                        dpi =
                          figure_dpi,
                        verbose =
                          TRUE
                      )

                    rm(
                      umap_plot_result
                    )

                    invisible(
                      gc(verbose = FALSE)
                    )
                  }

                  plots_created <-
                    plots_created + 1L
                }

                list(
                  stageStatus =
                    "completed",
                  plotsCreated =
                    plots_created,
                  plotsSkipped =
                    plots_skipped
                )
              }
            )
          },
          error = function(error_condition) {

            list(
              error =
                conditionMessage(
                  error_condition
                )
            )
          }
        )

        stage_02_elapsed <- as.numeric(
          difftime(
            Sys.time(),
            stage_02_started,
            units = "secs"
          )
        )

        if (
          !is.null(stage_02_result$error)
        ) {

          stage_02_status <- paste0(
            "error: ",
            stage_02_result$error
          )

          configuration_error <- paste(
            na.omit(
              c(
                configuration_error,
                stage_02_result$error
              )
            ),
            collapse = " | "
          )

          workflow_message(
            "STAGE 02 ERROR: ",
            stage_02_result$error
          )

        } else {

          stage_02_status <-
            stage_02_result$stageStatus
        }
      }


      # ========================================================================
      # STAGE 03: clustering validation
      # ========================================================================

      if (is.null(seurat_object)) {

        stage_03_status <-
          "skipped_stage01_failed"

        write_skipped_stage_log(
          stage_03_log,
          "03 clustering validation",
          "Stage 01 did not produce a valid Seurat object."
        )

      } else {

        validation_summary_file <- file.path(
          validation_directory,
          "tables",
          paste0(
            analysis_prefix,
            "_allMethods_clusteringValidationSummary.tsv"
          )
        )

        validation_plot_png <- file.path(
          validation_directory,
          "figures",
          paste0(
            analysis_prefix,
            "_allMethods_clusteringValidationMetrics.png"
          )
        )

        validation_plot_pdf <- file.path(
          validation_directory,
          "figures",
          paste0(
            analysis_prefix,
            "_allMethods_clusteringValidationMetrics.pdf"
          )
        )

        stage_03_started <- Sys.time()

        stage_03_result <- tryCatch(
          {
            run_logged_stage(
              log_file = stage_03_log,
              stage_label =
                "03 clustering validation",
              code = {

                validation_outputs_complete <-
                  all_required_files_exist(
                    c(
                      validation_summary_file,
                      if (save_validation_png) {
                        validation_plot_png
                      } else {
                        NA_character_
                      },
                      if (save_validation_pdf) {
                        validation_plot_pdf
                      } else {
                        NA_character_
                      }
                    )
                  )

                if (
                  isTRUE(
                    resume_completed_stages
                  ) &&
                    !isTRUE(
                      overwrite_existing_validation
                    ) &&
                    validation_outputs_complete
                ) {

                  workflow_message(
                    "Existing validation outputs detected. Skipping stage."
                  )

                  list(
                    stageStatus =
                      "skipped_complete"
                  )

                } else {

                  validation_result <-
                    run_clustering_validation(
                      seurat_object =
                        seurat_object,
                      cluster_columns =
                        expected_cluster_columns,
                      reduction =
                        integration_reduction,
                      dims =
                        seq_len(dims_value),
                      sample_order =
                        sample_order,
                      output_dir =
                        validation_directory,
                      analysis_prefix =
                        analysis_prefix,
                      integration_method =
                        integration_method,
                      normalization_label =
                        normalization_label,
                      n_hvg =
                        n_hvg,
                      k_param =
                        k_value,
                      prune_snn =
                        prune_snn,
                      transcriptomic_silhouette_max_spots =
                        transcriptomic_silhouette_max_spots,
                      transcriptomic_silhouette_minimum_spots_per_cluster =
                        transcriptomic_silhouette_minimum_spots_per_cluster,
                      spatial_silhouette_max_spots_per_sample =
                        spatial_silhouette_max_spots_per_sample,
                      spatial_silhouette_minimum_spots_per_cluster =
                        spatial_silhouette_minimum_spots_per_cluster,
                      pas_k =
                        pas_k,
                      pas_minimum_different_neighbours =
                        pas_minimum_different_neighbours,
                      normalize_chaos_by_spot_pitch =
                        normalize_chaos_by_spot_pitch,
                      ari_comparison_mode =
                        ari_comparison_mode,
                      seed =
                        random_seed,
                      save_png =
                        save_validation_png,
                      save_pdf =
                        save_validation_pdf,
                      dpi =
                        figure_dpi,
                      verbose =
                        TRUE
                    )

                  rm(
                    validation_result
                  )

                  invisible(
                    gc(verbose = FALSE)
                  )

                  list(
                    stageStatus =
                      "completed"
                  )
                }
              }
            )
          },
          error = function(error_condition) {

            list(
              error =
                conditionMessage(
                  error_condition
                )
            )
          }
        )

        stage_03_elapsed <- as.numeric(
          difftime(
            Sys.time(),
            stage_03_started,
            units = "secs"
          )
        )

        if (
          !is.null(stage_03_result$error)
        ) {

          stage_03_status <- paste0(
            "error: ",
            stage_03_result$error
          )

          configuration_error <- paste(
            na.omit(
              c(
                configuration_error,
                stage_03_result$error
              )
            ),
            collapse = " | "
          )

          workflow_message(
            "STAGE 03 ERROR: ",
            stage_03_result$error
          )

        } else {

          stage_03_status <-
            stage_03_result$stageStatus
        }
      }


    # --------------------------------------------------------------------------
    # Collect configuration-level summaries
    # --------------------------------------------------------------------------

    if (file.exists(clustering_summary_file)) {

      clustering_summary <- utils::read.delim(
        clustering_summary_file,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      combined_clustering_summaries[[graph_configuration_tag]] <- prepend_configuration_columns(
        data =
          clustering_summary,
        dims_value =
          dims_value,
        k_value =
          k_value,
        graph_configuration_tag =
          graph_configuration_tag
      )
    }

    validation_summary_file <- file.path(
      validation_directory,
      "tables",
      paste0(
        analysis_prefix,
        "_allMethods_clusteringValidationSummary.tsv"
      )
    )

    if (file.exists(validation_summary_file)) {

      validation_summary <- utils::read.delim(
        validation_summary_file,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      combined_validation_summaries[[graph_configuration_tag]] <- prepend_configuration_columns(
        data =
          validation_summary,
        dims_value =
          dims_value,
        k_value =
          k_value,
        graph_configuration_tag =
          graph_configuration_tag
      )
    }

    configuration_finished <- Sys.time()

    manifest_rows[[graph_configuration_tag]] <- data.frame(
      configurationIndex =
        configuration_index,
      rpcaDims =
        dims_value,
      graphK =
        k_value,
      graphConfiguration =
        graph_configuration_tag,
      parentRData =
        parent_rdata_file,
      integratedCacheRData =
        integrated_cache_rdata_file,
      inputRData =
        input_rdata_file,
      outputRData =
        output_rdata_file,
      stage01Status =
        stage_01_status,
      stage01ElapsedSeconds =
        stage_01_elapsed,
      stage01ElapsedFormatted =
        if (is.na(stage_01_elapsed)) {
          NA_character_
        } else {
          format_elapsed_time(
            stage_01_elapsed
          )
        },
      stage02Status =
        stage_02_status,
      stage02ElapsedSeconds =
        stage_02_elapsed,
      stage02ElapsedFormatted =
        if (is.na(stage_02_elapsed)) {
          NA_character_
        } else {
          format_elapsed_time(
            stage_02_elapsed
          )
        },
      stage03Status =
        stage_03_status,
      stage03ElapsedSeconds =
        stage_03_elapsed,
      stage03ElapsedFormatted =
        if (is.na(stage_03_elapsed)) {
          NA_character_
        } else {
          format_elapsed_time(
            stage_03_elapsed
          )
        },
      configurationElapsedSeconds =
        as.numeric(
          difftime(
            configuration_finished,
            configuration_started,
            units = "secs"
          )
        ),
      configurationElapsedFormatted =
        format_elapsed_time(
          as.numeric(
            difftime(
              configuration_finished,
              configuration_started,
              units = "secs"
            )
          )
        ),
      errorMessage =
        configuration_error,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    current_manifest <- do.call(
      rbind,
      manifest_rows
    )

    rownames(current_manifest) <- NULL

    write_tsv(
      current_manifest,
      manifest_file
    )

    rm(
      seurat_object
    )

    invisible(
      gc(verbose = FALSE)
    )

    if (
      !is.na(configuration_error) &&
        !isTRUE(
          continue_after_configuration_error
        )
    ) {
      stop(
        "Configuration failed and `continue_after_configuration_error = FALSE`: ",
        graph_configuration_tag,
        "\n",
        configuration_error,
        call. = FALSE
      )
    }
  }


  # ============================================================================
  # 6. Grid-wide summaries
  # ============================================================================

  if (
    length(
      combined_clustering_summaries
    ) > 0L
  ) {

    combined_clustering_summary <- do.call(
      rbind,
      combined_clustering_summaries
    )

    rownames(
      combined_clustering_summary
    ) <- NULL

    write_tsv(
      combined_clustering_summary,
      file.path(
        grid_tables_directory,
        paste0(
          grid_prefix,
          "_allConfigurations_clusteringSummary.tsv"
        )
      )
    )
  }

  if (
    length(
      combined_validation_summaries
    ) > 0L
  ) {

    combined_validation_summary <- do.call(
      rbind,
      combined_validation_summaries
    )

    rownames(
      combined_validation_summary
    ) <- NULL

    write_tsv(
      combined_validation_summary,
      file.path(
        grid_tables_directory,
        paste0(
          grid_prefix,
          "_allConfigurations_clusteringValidationSummary.tsv"
        )
      )
    )
  }


  # ============================================================================
  # 7. Final report
  # ============================================================================

  final_manifest <- if (
    length(manifest_rows) > 0L
  ) {
    do.call(
      rbind,
      manifest_rows
    )
  } else {
    data.frame()
  }

  failed_stage_pattern <- "^error:|^skipped_missing_input"

  failed_configurations <- if (
    nrow(final_manifest) > 0L
  ) {
    final_manifest[
      grepl(
        failed_stage_pattern,
        final_manifest$stage01Status
      ) |
        grepl(
          "^error:",
          final_manifest$stage02Status
        ) |
        grepl(
          "^error:",
          final_manifest$stage03Status
        ),
      ,
      drop = FALSE
    ]
  } else {
    final_manifest
  }

  analysis_finished <- Sys.time()

  workflow_message("")
  workflow_message(
    "============================================================"
  )

  workflow_message(
    "FULL RPCA DIMS x K GRID WORKFLOW FINISHED"
  )

  workflow_message(
    "Configurations processed: ",
    nrow(configuration_grid)
  )

  workflow_message(
    "Failed configurations: ",
    nrow(failed_configurations)
  )

  workflow_message(
    "Grid manifest: ",
    manifest_file
  )

  workflow_message(
    "Grid summary directory: ",
    grid_summary_root
  )

  workflow_message(
    "Total elapsed time: ",
    format_elapsed_time(
      as.numeric(
        difftime(
          analysis_finished,
          analysis_started,
          units = "secs"
        )
      )
    )
  )

  workflow_message(
    "============================================================"
  )

  if (
    nrow(failed_configurations) > 0L &&
      isTRUE(
        fail_at_end_if_any_stage_failed
      )
  ) {
    stop(
      "One or more graph configurations failed. Inspect the run manifest and ",
      "the three stage logs in each configuration-specific `logs/` folder.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


tryCatch(
  main(),
  error = function(error_condition) {

    cat(
      "ERROR: ",
      conditionMessage(
        error_condition
      ),
      "\n",
      sep = ""
    )

    flush.console()

    quit(
      save = "no",
      status = 1L,
      runLast = FALSE
    )
  }
)
