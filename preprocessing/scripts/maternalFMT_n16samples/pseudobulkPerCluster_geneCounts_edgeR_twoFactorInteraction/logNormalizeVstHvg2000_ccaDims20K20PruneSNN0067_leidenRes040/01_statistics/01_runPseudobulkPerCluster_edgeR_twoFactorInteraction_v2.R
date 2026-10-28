#!/usr/bin/env Rscript

# ==============================================================================
# 01_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction_fixed_withLog2QN.R
#
# PURPOSE
#   1. Run the existing cluster-specific pseudobulk edgeR analysis unchanged:
#
#        expression ~ fmt_donor_group * sex
#
#   2. After the edgeR workflow finishes, take the per-cluster RAW pseudobulk
#      count matrices already written by the analysis and, for EACH CLUSTER
#      SEPARATELY:
#
#        raw pseudobulk counts
#          -> log2(count + 1)
#          -> quantile normalization across biological samples
#          -> save the final log2 + quantile-normalized matrix
#
#   3. Add mean and sample SD of log2(count + 1) + QN values to every edgeR
#      full-result table using the biological levels relevant to each test.
#
# IMPORTANT
#   - edgeR itself is NOT run on log2-transformed or quantile-normalized data.
#   - edgeR continues to use raw counts, filterByExpr, TMM and the QL GLM.
#   - The additional log2 + QN matrices are used for visualization,
#     heatmaps, PCA and descriptive group-level mean/SD summaries.
#   - Quantile normalization is performed independently for every cluster.
#
# INPUT RAW-COUNT DIRECTORY CREATED BY THE edgeR WORKFLOW
#   01_statistics/
#   01_pseudobulkInputs/
#   perCluster_rawCounts/
#
# NEW OUTPUT DIRECTORY
#   01_statistics/
#   01_pseudobulkInputs/
#   perCluster_rawCounts_log2plus1Transformation_quantileNormalize/
# ==============================================================================


# ==============================================================================
# 1. Project and analysis configuration
# ==============================================================================

options(stringsAsFactors = FALSE)

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"

dataset_name <- "maternalFMT_n16samples"

clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"


# ==============================================================================
# 2. Input files
# ==============================================================================

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

metadata_file <- file.path(
  project_root,
  "data",
  "metadata_autismFMT.tsv"
)

annotation_rds_file <- file.path(
  project_root,
  "data",
  "ensembl115",
  "ensembl115_mouse_proteinCodingGenes.rds"
)

annotation_metadata_file <- file.path(
  project_root,
  "data",
  "ensembl115",
  "ensembl115_mouse_proteinCodingGenes_metadata.tsv"
)

functions_file <- file.path(
  project_root,
  "preprocessing",
  "src",
  "pseudobulkPerCluster_geneCounts_edgeR",
  "functions_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction.R"
)


# ==============================================================================
# 3. Output directory
# ==============================================================================

output_statistics_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)


# ==============================================================================
# 4. Validate inputs and load functions
# ==============================================================================

required_input_files <- c(
  input_rdata_file,
  metadata_file,
  annotation_rds_file,
  annotation_metadata_file,
  functions_file
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

source(functions_file)

if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  stop(
    paste0(
      "Package 'preprocessCore' is required for quantile normalization.\n",
      "Install it with:\n",
      "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) ",
      "install.packages(\"BiocManager\")\n",
      "BiocManager::install(\"preprocessCore\")"
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 5. Explicit analysis settings
# ==============================================================================

seurat_object_name <- NULL

# Explicit clustering column.
cluster_column <- "leiden_res0.40"

# A spot is positive for a gene when the raw count is at least this value.
# This affects descriptive percent-positive metrics only.
positive_spot_min_count <- 1L

# No sample-cluster pseudobulk is removed on the basis of spot count.
min_spots_per_sample_cluster <- 0L

# XLSX display precision only.
xlsx_max_decimal_places <- 5L

# Expected dimensions.
expected_number_of_samples <- 16L
expected_number_of_clusters <- 16L


# ==============================================================================
# 6. Helper functions for log2(count + 1) + quantile normalization
# ==============================================================================

write_tsv_gz <- function(x, output_file) {

  output_connection <- gzfile(
    output_file,
    open = "wt",
    encoding = "UTF-8"
  )

  on.exit(
    close(output_connection),
    add = TRUE
  )

  utils::write.table(
    x,
    file = output_connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )

  invisible(output_file)
}


extract_cluster_id_from_raw_count_filename <- function(file_name) {

  cluster_id <- sub(
    "^.*_cluster_([^_]+)_rawCounts_proteinCodingAutosomalAndX\\.tsv\\.gz$",
    "\\1",
    basename(file_name)
  )

  if (identical(cluster_id, basename(file_name))) {
    stop(
      "Could not extract cluster ID from file name: ",
      basename(file_name),
      call. = FALSE
    )
  }

  cluster_id
}


natural_cluster_order <- function(cluster_ids) {

  numeric_part <- suppressWarnings(
    as.integer(
      sub("^.*?([0-9]+)$", "\\1", cluster_ids)
    )
  )

  if (all(!is.na(numeric_part))) {
    return(cluster_ids[order(numeric_part)])
  }

  sort(cluster_ids)
}


prepare_log2plus1_quantile_normalized_raw_counts <- function(
    output_statistics_dir,
    expected_number_of_samples = 16L,
    expected_number_of_clusters = 16L
) {

  raw_count_dir <- file.path(
    output_statistics_dir,
    "01_pseudobulkInputs",
    "perCluster_rawCounts"
  )

  output_dir <- file.path(
    output_statistics_dir,
    "01_pseudobulkInputs",
    "perCluster_rawCounts_log2plus1Transformation_quantileNormalize"
  )

  if (!dir.exists(raw_count_dir)) {
    stop(
      "Raw pseudobulk-count directory does not exist:\n",
      raw_count_dir,
      call. = FALSE
    )
  }

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  raw_count_files <- list.files(
    raw_count_dir,
    pattern = "_rawCounts_proteinCodingAutosomalAndX\\.tsv\\.gz$",
    full.names = TRUE
  )

  if (length(raw_count_files) == 0L) {
    stop(
      "No per-cluster raw-count files were found in:\n",
      raw_count_dir,
      call. = FALSE
    )
  }

  cluster_ids <- vapply(
    raw_count_files,
    extract_cluster_id_from_raw_count_filename,
    character(1)
  )

  if (anyDuplicated(cluster_ids)) {
    duplicated_clusters <- unique(
      cluster_ids[duplicated(cluster_ids)]
    )

    stop(
      "Duplicated raw-count files were found for cluster(s): ",
      paste(duplicated_clusters, collapse = ", "),
      call. = FALSE
    )
  }

  ordered_cluster_ids <- natural_cluster_order(cluster_ids)

  raw_count_files <- raw_count_files[
    match(ordered_cluster_ids, cluster_ids)
  ]

  cluster_ids <- ordered_cluster_ids

  if (length(raw_count_files) != expected_number_of_clusters) {
    stop(
      "Expected ",
      expected_number_of_clusters,
      " per-cluster raw-count files but found ",
      length(raw_count_files),
      ".\nFiles:\n",
      paste(basename(raw_count_files), collapse = "\n"),
      call. = FALSE
    )
  }

  message("\n============================================================")
  message("LOG2(COUNT + 1) + QUANTILE NORMALIZATION")
  message("============================================================")
  message("Input directory:")
  message(raw_count_dir)
  message("\nOutput directory:")
  message(output_dir)
  message("\nClusters: ", paste(cluster_ids, collapse = ", "))

  summary_rows <- vector(
    mode = "list",
    length = length(raw_count_files)
  )

  names(summary_rows) <- cluster_ids

  metadata_columns <- c(
    "ensembl_gene_id",
    "gene",
    "chromosome"
  )

  for (i in seq_along(raw_count_files)) {

    input_file <- raw_count_files[[i]]
    cluster_id <- cluster_ids[[i]]

    message("\n------------------------------------------------------------")
    message(
      "[",
      i,
      "/",
      length(raw_count_files),
      "] Cluster ",
      cluster_id
    )
    message("Reading: ", basename(input_file))

    raw_df <- utils::read.delim(
      input_file,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    missing_metadata_columns <- setdiff(
      metadata_columns,
      colnames(raw_df)
    )

    if (length(missing_metadata_columns) > 0L) {
      stop(
        "Cluster ",
        cluster_id,
        " is missing required metadata column(s): ",
        paste(missing_metadata_columns, collapse = ", "),
        call. = FALSE
      )
    }

    sample_columns <- setdiff(
      colnames(raw_df),
      metadata_columns
    )

    if (length(sample_columns) != expected_number_of_samples) {
      stop(
        "Cluster ",
        cluster_id,
        " has ",
        length(sample_columns),
        " sample columns; expected ",
        expected_number_of_samples,
        ".\nDetected sample columns:\n",
        paste(sample_columns, collapse = "\n"),
        call. = FALSE
      )
    }

    non_numeric_sample_columns <- sample_columns[
      !vapply(
        raw_df[sample_columns],
        is.numeric,
        logical(1)
      )
    ]

    if (length(non_numeric_sample_columns) > 0L) {

      raw_df[non_numeric_sample_columns] <- lapply(
        raw_df[non_numeric_sample_columns],
        function(x) {
          suppressWarnings(
            as.numeric(x)
          )
        }
      )
    }

    raw_count_matrix <- as.matrix(
      raw_df[
        ,
        sample_columns,
        drop = FALSE
      ]
    )

    storage.mode(raw_count_matrix) <- "double"

    if (anyNA(raw_count_matrix)) {
      stop(
        "NA values were found after numeric conversion in cluster ",
        cluster_id,
        ".",
        call. = FALSE
      )
    }

    if (any(!is.finite(raw_count_matrix))) {
      stop(
        "Non-finite raw counts were found in cluster ",
        cluster_id,
        ".",
        call. = FALSE
      )
    }

    if (any(raw_count_matrix < 0)) {
      stop(
        "Negative raw counts were found in cluster ",
        cluster_id,
        ".",
        call. = FALSE
      )
    }

    # --------------------------------------------------------------------------
    # Step 1: log2(count + 1)
    # --------------------------------------------------------------------------

    log2plus1_matrix <- log2(
      raw_count_matrix + 1
    )

    # --------------------------------------------------------------------------
    # Step 2: quantile normalization across samples WITHIN this cluster
    # --------------------------------------------------------------------------

    quantile_normalized_matrix <-
      preprocessCore::normalize.quantiles(
        log2plus1_matrix,
        copy = TRUE
      )

    rownames(quantile_normalized_matrix) <-
      rownames(log2plus1_matrix)

    colnames(quantile_normalized_matrix) <-
      colnames(log2plus1_matrix)

    if (
      !identical(
        dim(quantile_normalized_matrix),
        dim(raw_count_matrix)
      )
    ) {
      stop(
        "Dimension mismatch after quantile normalization for cluster ",
        cluster_id,
        ".",
        call. = FALSE
      )
    }

    if (any(!is.finite(quantile_normalized_matrix))) {
      stop(
        "Non-finite values were generated after quantile normalization for cluster ",
        cluster_id,
        ".",
        call. = FALSE
      )
    }

    output_df <- data.frame(
      raw_df[
        ,
        metadata_columns,
        drop = FALSE
      ],
      quantile_normalized_matrix,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    output_file <- file.path(
      output_dir,
      sub(
        "_rawCounts_proteinCodingAutosomalAndX\\.tsv\\.gz$",
        paste0(
          "_rawCounts_log2plus1Transformation_",
          "quantileNormalize_",
          "proteinCodingAutosomalAndX.tsv.gz"
        ),
        basename(input_file)
      )
    )

    write_tsv_gz(
      output_df,
      output_file
    )

    summary_rows[[i]] <- data.frame(
      cluster_id = cluster_id,
      input_file = basename(input_file),
      output_file = basename(output_file),
      number_of_genes = nrow(raw_count_matrix),
      number_of_samples = ncol(raw_count_matrix),
      raw_count_min = min(raw_count_matrix),
      raw_count_max = max(raw_count_matrix),
      log2plus1_min = min(log2plus1_matrix),
      log2plus1_max = max(log2plus1_matrix),
      quantile_normalized_min = min(quantile_normalized_matrix),
      quantile_normalized_max = max(quantile_normalized_matrix),
      stringsAsFactors = FALSE
    )

    message(
      "Saved: ",
      basename(output_file)
    )
    message(
      "Genes: ",
      format(nrow(raw_count_matrix), big.mark = ","),
      " | samples: ",
      ncol(raw_count_matrix)
    )

    rm(
      raw_df,
      raw_count_matrix,
      log2plus1_matrix,
      quantile_normalized_matrix,
      output_df
    )

    invisible(gc())
  }

  normalization_summary <- do.call(
    rbind,
    summary_rows
  )

  rownames(normalization_summary) <- NULL

  summary_file <- file.path(
    output_dir,
    paste0(
      dataset_name,
      "_",
      clustering_name,
      "_perCluster_rawCounts_",
      "log2plus1Transformation_quantileNormalize_",
      "summary.tsv"
    )
  )

  utils::write.table(
    normalization_summary,
    file = summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )

  definitions_file <- file.path(
    output_dir,
    paste0(
      dataset_name,
      "_",
      clustering_name,
      "_perCluster_rawCounts_",
      "log2plus1Transformation_quantileNormalize_",
      "definitions.txt"
    )
  )

  writeLines(
    c(
      "Input: per-cluster raw pseudobulk gene-count matrices.",
      "Genes: protein-coding genes on chromosomes 1-19 and X.",
      "Input raw-count files are taken before filterByExpr.",
      "Transformation step: log2(raw_count + 1).",
      paste0(
        "Normalization step: preprocessCore::normalize.quantiles() ",
        "across the ",
        expected_number_of_samples,
        " biological samples."
      ),
      "Quantile normalization is performed independently for each cluster.",
      "No row-wise Z-score is applied in this script.",
      "These transformed matrices are descriptive/visualization outputs only.",
      "They are not used as input to edgeR.",
      "edgeR remains based on raw counts -> filterByExpr -> TMM -> QL GLM."
    ),
    con = definitions_file
  )

  message("\n============================================================")
  message("LOG2 + QUANTILE NORMALIZATION FINISHED")
  message("============================================================")
  message("Output directory:")
  message(output_dir)
  message("\nSummary:")
  message(summary_file)
  message("\nDefinitions:")
  message(definitions_file)

  list(
    input_dir = raw_count_dir,
    output_dir = output_dir,
    output_files = file.path(
      output_dir,
      normalization_summary$output_file
    ),
    summary_file = summary_file,
    definitions_file = definitions_file,
    summary = normalization_summary
  )
}



# ==============================================================================
# 7. Add group-wise mean and SD of log2(count + 1) + QN to edgeR result tables
# ==============================================================================

canonical_cluster_id <- function(x) {
  sub("^C", "", as.character(x), ignore.case = TRUE)
}

safe_row_sd_log2qn <- function(x) {
  if (ncol(x) <= 1L) {
    return(rep(NA_real_, nrow(x)))
  }
  apply(x, 1L, stats::sd)
}

make_log2qn_descriptive_tables <- function(
    log2_qn_matrix,
    model_metadata,
    test_definitions
) {
  model_metadata <- model_metadata |>
    dplyr::mutate(
      sample_ID = as.character(.data$sample_ID),
      fmt_donor_group = as.character(.data$fmt_donor_group),
      sex = as.character(.data$sex)
    )

  if (!all(model_metadata$sample_ID %in% colnames(log2_qn_matrix))) {
    missing_samples <- setdiff(
      model_metadata$sample_ID,
      colnames(log2_qn_matrix)
    )
    stop(
      "Missing samples in log2+QN matrix: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  tested_gene_ids <- rownames(log2_qn_matrix)

  summarize_samples <- function(sample_ids, displayed_name) {
    sample_ids <- as.character(sample_ids)

    if (length(sample_ids) == 0L) {
      stop(
        "No samples were selected for log2+QN descriptive level: ",
        displayed_name,
        call. = FALSE
      )
    }

    values_subset <- log2_qn_matrix[
      ,
      sample_ids,
      drop = FALSE
    ]

    tibble::tibble(
      ensembl_gene_id = tested_gene_ids,
      !!paste0("mean_log2plus1_QN_", displayed_name) :=
        rowMeans(values_subset),
      !!paste0("sd_log2plus1_QN_", displayed_name) :=
        safe_row_sd_log2qn(values_subset)
    )
  }

  select_samples <- function(
      sex_value = NULL,
      group_value = NULL
  ) {
    keep <- rep(TRUE, nrow(model_metadata))

    if (!is.null(sex_value)) {
      keep <- keep & model_metadata$sex == sex_value
    }

    if (!is.null(group_value)) {
      keep <- keep & model_metadata$fmt_donor_group == group_value
    }

    model_metadata$sample_ID[keep]
  }

  join_metric_tables <- function(...) {
    metric_tables <- list(...)
    Reduce(
      function(x, y) {
        dplyr::left_join(x, y, by = "ensembl_gene_id")
      },
      metric_tables
    )
  }

  descriptive_tables <- list(
    Interaction = join_metric_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Male_Neurotypical"
      ),
      summarize_samples(
        select_samples("Male", "ASD"),
        "Male_ASD"
      ),
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Female_Neurotypical"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "Female_ASD"
      )
    ),

    Overall_Sex_Female_vs_Male = join_metric_tables(
      summarize_samples(
        select_samples(sex_value = "Male"),
        "Male"
      ),
      summarize_samples(
        select_samples(sex_value = "Female"),
        "Female"
      )
    ),

    Overall_Group_ASD_vs_Neurotypical = join_metric_tables(
      summarize_samples(
        select_samples(group_value = "Neurotypical"),
        "Neurotypical"
      ),
      summarize_samples(
        select_samples(group_value = "ASD"),
        "ASD"
      )
    ),

    ASD_Male_vs_Neurotypical_Male = join_metric_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Neurotypical_Male"
      ),
      summarize_samples(
        select_samples("Male", "ASD"),
        "ASD_Male"
      )
    ),

    ASD_Female_vs_Neurotypical_Female = join_metric_tables(
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Neurotypical_Female"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "ASD_Female"
      )
    ),

    Neurotypical_Female_vs_Neurotypical_Male = join_metric_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Neurotypical_Male"
      ),
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Neurotypical_Female"
      )
    ),

    ASD_Female_vs_ASD_Male = join_metric_tables(
      summarize_samples(
        select_samples("Male", "ASD"),
        "ASD_Male"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "ASD_Female"
      )
    ),

    ASD_Female_vs_Neurotypical_Male = join_metric_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Neurotypical_Male"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "ASD_Female"
      )
    ),

    ASD_Male_vs_Neurotypical_Female = join_metric_tables(
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Neurotypical_Female"
      ),
      summarize_samples(
        select_samples("Male", "ASD"),
        "ASD_Male"
      )
    )
  )

  expected_test_ids <- as.character(test_definitions$test_id)

  if (!setequal(names(descriptive_tables), expected_test_ids)) {
    stop(
      "log2+QN descriptive table names do not match edgeR test definitions.",
      call. = FALSE
    )
  }

  descriptive_tables[expected_test_ids]
}

add_log2qn_descriptive_metrics_to_edger_outputs <- function(
    output_statistics_dir,
    log2_qn_output,
    xlsx_max_decimal_places = 5L
) {
  analysis_prefix <- paste0(
    dataset_name,
    "_",
    clustering_name,
    "_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction"
  )

  result_rdata_file <- file.path(
    output_statistics_dir,
    "03_edgeRResults",
    paste0(
      analysis_prefix,
      "_edgeRResults.RData"
    )
  )

  model_rdata_file <- file.path(
    output_statistics_dir,
    "05_objects",
    paste0(
      analysis_prefix,
      "_edgeRModelObjects.RData"
    )
  )

  if (!file.exists(result_rdata_file)) {
    stop(
      "Missing edgeR results RData:\n",
      result_rdata_file,
      call. = FALSE
    )
  }

  if (!file.exists(model_rdata_file)) {
    stop(
      "Missing edgeR model-object RData:\n",
      model_rdata_file,
      call. = FALSE
    )
  }

  result_env <- new.env(parent = emptyenv())
  model_env <- new.env(parent = emptyenv())

  load(
    result_rdata_file,
    envir = result_env
  )

  load(
    model_rdata_file,
    envir = model_env
  )

  required_result_objects <- c(
    "edgeR_perCluster_results",
    "edgeR_perCluster_combinedResults",
    "edgeR_perCluster_analysisSummary",
    "edgeR_perCluster_clusterStatus",
    "edgeR_perCluster_testDefinitions",
    "edgeR_perCluster_geneAbundanceQC",
    "edgeR_perCluster_combinedGeneAbundanceQC",
    "edgeR_perCluster_descriptiveMetricDefinitions"
  )

  missing_result_objects <- required_result_objects[
    !vapply(
      required_result_objects,
      exists,
      logical(1),
      envir = result_env,
      inherits = FALSE
    )
  ]

  if (length(missing_result_objects) > 0L) {
    stop(
      "Missing object(s) in edgeR results RData: ",
      paste(missing_result_objects, collapse = ", "),
      call. = FALSE
    )
  }

  if (
    !exists(
      "edgeR_perCluster_modelObjects",
      envir = model_env,
      inherits = FALSE
    )
  ) {
    stop(
      "Object `edgeR_perCluster_modelObjects` is missing from model RData.",
      call. = FALSE
    )
  }

  edgeR_perCluster_results <- get(
    "edgeR_perCluster_results",
    envir = result_env
  )

  edgeR_perCluster_combinedResults <- get(
    "edgeR_perCluster_combinedResults",
    envir = result_env
  )

  edgeR_perCluster_analysisSummary <- get(
    "edgeR_perCluster_analysisSummary",
    envir = result_env
  )

  edgeR_perCluster_clusterStatus <- get(
    "edgeR_perCluster_clusterStatus",
    envir = result_env
  )

  edgeR_perCluster_testDefinitions <- get(
    "edgeR_perCluster_testDefinitions",
    envir = result_env
  )

  edgeR_perCluster_geneAbundanceQC <- get(
    "edgeR_perCluster_geneAbundanceQC",
    envir = result_env
  )

  edgeR_perCluster_combinedGeneAbundanceQC <- get(
    "edgeR_perCluster_combinedGeneAbundanceQC",
    envir = result_env
  )

  edgeR_perCluster_descriptiveMetricDefinitions <- get(
    "edgeR_perCluster_descriptiveMetricDefinitions",
    envir = result_env
  )

  edgeR_perCluster_modelObjects <- get(
    "edgeR_perCluster_modelObjects",
    envir = model_env
  )

  completed_cluster_ids <- edgeR_perCluster_clusterStatus |>
    dplyr::filter(.data$status == "completed") |>
    dplyr::pull(.data$cluster_id) |>
    as.character()

  if (length(completed_cluster_ids) == 0L) {
    stop(
      "No completed clusters were found while adding log2+QN metrics.",
      call. = FALSE
    )
  }

  qn_summary <- log2_qn_output$summary

  if (
    is.null(qn_summary) ||
      !all(c("cluster_id", "output_file") %in% colnames(qn_summary))
  ) {
    stop(
      "log2_qn_output$summary does not contain cluster_id/output_file.",
      call. = FALSE
    )
  }

  message("\n============================================================")
  message("ADDING log2(count + 1) + QN GROUP MEAN/SD TO edgeR TABLES")
  message("============================================================")

  for (cluster_id_current in completed_cluster_ids) {
    canonical_current <- canonical_cluster_id(cluster_id_current)

    qn_matches <- which(
      canonical_cluster_id(qn_summary$cluster_id) == canonical_current
    )

    if (length(qn_matches) != 1L) {
      stop(
        "Expected exactly one log2+QN file for cluster ",
        cluster_id_current,
        " but found ",
        length(qn_matches),
        ".",
        call. = FALSE
      )
    }

    qn_file <- file.path(
      log2_qn_output$output_dir,
      qn_summary$output_file[[qn_matches]]
    )

    if (!file.exists(qn_file)) {
      stop(
        "Missing log2+QN file for cluster ",
        cluster_id_current,
        ":\n",
        qn_file,
        call. = FALSE
      )
    }

    message(
      "Cluster ",
      cluster_id_current,
      ": ",
      basename(qn_file)
    )

    qn_df <- utils::read.delim(
      qn_file,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    required_gene_columns <- c(
      "ensembl_gene_id",
      "gene",
      "chromosome"
    )

    missing_gene_columns <- setdiff(
      required_gene_columns,
      colnames(qn_df)
    )

    if (length(missing_gene_columns) > 0L) {
      stop(
        "log2+QN file for cluster ",
        cluster_id_current,
        " is missing: ",
        paste(missing_gene_columns, collapse = ", "),
        call. = FALSE
      )
    }

    sample_columns <- setdiff(
      colnames(qn_df),
      required_gene_columns
    )

    qn_matrix <- as.matrix(
      qn_df[
        ,
        sample_columns,
        drop = FALSE
      ]
    )

    storage.mode(qn_matrix) <- "double"

    rownames(qn_matrix) <- qn_df$ensembl_gene_id

    if (
      anyNA(qn_matrix) ||
        any(!is.finite(qn_matrix))
    ) {
      stop(
        "NA or non-finite values found in log2+QN matrix for cluster ",
        cluster_id_current,
        ".",
        call. = FALSE
      )
    }

    model_object <- edgeR_perCluster_modelObjects[[cluster_id_current]]

    if (is.null(model_object)) {
      model_match <- names(edgeR_perCluster_modelObjects)[
        canonical_cluster_id(
          names(edgeR_perCluster_modelObjects)
        ) == canonical_current
      ]

      if (length(model_match) == 1L) {
        model_object <- edgeR_perCluster_modelObjects[[model_match]]
      }
    }

    if (is.null(model_object) || is.null(model_object$dge)) {
      stop(
        "DGEList model object not found for cluster ",
        cluster_id_current,
        ".",
        call. = FALSE
      )
    }

    model_metadata <- model_object$dge$samples |>
      tibble::rownames_to_column("dge_rowname") |>
      tibble::as_tibble()

    required_metadata_columns <- c(
      "sample_ID",
      "fmt_donor_group",
      "sex"
    )

    missing_metadata_columns <- setdiff(
      required_metadata_columns,
      colnames(model_metadata)
    )

    if (length(missing_metadata_columns) > 0L) {
      stop(
        "DGEList metadata for cluster ",
        cluster_id_current,
        " is missing: ",
        paste(missing_metadata_columns, collapse = ", "),
        call. = FALSE
      )
    }

    tested_gene_ids <- rownames(model_object$dge$counts)

    missing_tested_genes <- setdiff(
      tested_gene_ids,
      rownames(qn_matrix)
    )

    if (length(missing_tested_genes) > 0L) {
      stop(
        "log2+QN file for cluster ",
        cluster_id_current,
        " is missing ",
        length(missing_tested_genes),
        " tested edgeR gene(s).",
        call. = FALSE
      )
    }

    qn_matrix_tested <- qn_matrix[
      tested_gene_ids,
      ,
      drop = FALSE
    ]

    qn_descriptive_tables <- make_log2qn_descriptive_tables(
      log2_qn_matrix = qn_matrix_tested,
      model_metadata = model_metadata,
      test_definitions = edgeR_perCluster_testDefinitions
    )

    result_cluster_name <- names(edgeR_perCluster_results)[
      canonical_cluster_id(
        names(edgeR_perCluster_results)
      ) == canonical_current
    ]

    if (length(result_cluster_name) != 1L) {
      stop(
        "Could not uniquely match edgeR result list for cluster ",
        cluster_id_current,
        ".",
        call. = FALSE
      )
    }

    for (test_id in edgeR_perCluster_testDefinitions$test_id) {
      result_table <- edgeR_perCluster_results[[result_cluster_name]][[test_id]]

      if (is.null(result_table)) {
        stop(
          "Missing result table for cluster ",
          cluster_id_current,
          ", test ",
          test_id,
          ".",
          call. = FALSE
        )
      }

      old_qn_columns <- grep(
        "^(mean|sd)_log2plus1_QN_",
        colnames(result_table),
        value = TRUE
      )

      if (length(old_qn_columns) > 0L) {
        result_table <- result_table |>
          dplyr::select(
            -dplyr::all_of(old_qn_columns)
          )
      }

      result_table <- result_table |>
        dplyr::left_join(
          qn_descriptive_tables[[test_id]],
          by = "ensembl_gene_id"
        )

      new_qn_columns <- grep(
        "^(mean|sd)_log2plus1_QN_",
        colnames(result_table),
        value = TRUE
      )

      if (
        length(new_qn_columns) == 0L ||
          anyNA(
            result_table[
              ,
              new_qn_columns,
              drop = FALSE
            ]
          )
      ) {
        stop(
          "Missing log2+QN descriptive metrics after join for cluster ",
          cluster_id_current,
          ", test ",
          test_id,
          ".",
          call. = FALSE
        )
      }

      expected_qn_columns <- if (
        identical(test_id, "Interaction")
      ) {
        8L
      } else {
        4L
      }

      if (length(new_qn_columns) != expected_qn_columns) {
        stop(
          "Unexpected number of log2+QN columns for cluster ",
          cluster_id_current,
          ", test ",
          test_id,
          ": found ",
          length(new_qn_columns),
          ", expected ",
          expected_qn_columns,
          ".",
          call. = FALSE
        )
      }

      edgeR_perCluster_results[[result_cluster_name]][[test_id]] <- result_table
    }

    rm(
      qn_df,
      qn_matrix,
      qn_matrix_tested,
      qn_descriptive_tables
    )

    invisible(gc())
  }

  cluster_order <- as.character(
    edgeR_perCluster_clusterStatus$cluster_id
  )

  edgeR_perCluster_combinedResults <- lapply(
    edgeR_perCluster_testDefinitions$test_id,
    function(test_id) {
      combined_test_result <- dplyr::bind_rows(
        lapply(
          completed_cluster_ids,
          function(cluster_id_current) {
            canonical_current <- canonical_cluster_id(
              cluster_id_current
            )

            result_cluster_name <- names(
              edgeR_perCluster_results
            )[
              canonical_cluster_id(
                names(edgeR_perCluster_results)
              ) == canonical_current
            ]

            edgeR_perCluster_results[[result_cluster_name]][[test_id]] |>
              dplyr::mutate(
                cluster_id = as.character(
                  cluster_id_current
                ),
                .before = 1
              )
          }
        )
      )

      combined_test_result |>
        dplyr::mutate(
          cluster_id = factor(
            .data$cluster_id,
            levels = cluster_order
          )
        ) |>
        dplyr::arrange(
          .data$cluster_id,
          .data$PValue,
          dplyr::desc(abs(.data$logFC))
        ) |>
        dplyr::mutate(
          cluster_id = as.character(
            .data$cluster_id
          )
        )
    }
  )

  names(edgeR_perCluster_combinedResults) <-
    edgeR_perCluster_testDefinitions$test_id

  qn_definition_rows <- tibble::tribble(
    ~metric, ~definition, ~weighting_and_scope,
    "mean_log2plus1_QN",
    paste0(
      "Arithmetic mean of sample-level values after ",
      "log2(raw pseudobulk count + 1) transformation and ",
      "quantile normalization across all biological samples ",
      "within each cluster."
    ),
    paste0(
      "Calculated within the biological level displayed for each test. ",
      "Interaction uses the four Sex x donor-group cells. ",
      "Main effects use the two compared levels. ",
      "Post hoc tests use the two specific Sex x donor-group cells. ",
      "Descriptive only; not used by edgeR."
    ),
    "sd_log2plus1_QN",
    paste0(
      "Sample standard deviation of sample-level values after ",
      "log2(raw pseudobulk count + 1) transformation and ",
      "quantile normalization."
    ),
    paste0(
      "Calculated across biological samples within the displayed level. ",
      "Descriptive only; not used by edgeR."
    )
  )

  edgeR_perCluster_descriptiveMetricDefinitions <-
    edgeR_perCluster_descriptiveMetricDefinitions |>
    dplyr::filter(
      !.data$metric %in% c(
        "mean_log2plus1_QN",
        "sd_log2plus1_QN"
      )
    ) |>
    dplyr::bind_rows(
      qn_definition_rows
    )

  per_cluster_xlsx_dir <- file.path(
    output_statistics_dir,
    "03_edgeRResults",
    "01_perCluster_XLSX"
  )

  combined_results_dir <- file.path(
    output_statistics_dir,
    "03_edgeRResults",
    "02_combinedResults"
  )

  for (cluster_id_current in completed_cluster_ids) {
    canonical_current <- canonical_cluster_id(
      cluster_id_current
    )

    result_cluster_name <- names(
      edgeR_perCluster_results
    )[
      canonical_cluster_id(
        names(edgeR_perCluster_results)
      ) == canonical_current
    ]

    cluster_slug <-
      sanitize_file_component_pseudobulk_per_cluster(
        cluster_id_current
      )

    cluster_xlsx_file <- file.path(
      per_cluster_xlsx_dir,
      paste0(
        analysis_prefix,
        "_cluster_",
        cluster_slug,
        "_fullResults_nineSheets.xlsx"
      )
    )

    write_result_workbook_per_cluster(
      result_list =
        edgeR_perCluster_results[[result_cluster_name]],
      test_definitions =
        edgeR_perCluster_testDefinitions,
      output_file = cluster_xlsx_file,
      include_cluster_id = FALSE,
      xlsx_max_decimal_places =
        xlsx_max_decimal_places
    )
  }

  combined_xlsx_file <- file.path(
    combined_results_dir,
    paste0(
      analysis_prefix,
      "_allClusters_fullResults_nineSheets.xlsx"
    )
  )

  write_result_workbook_per_cluster(
    result_list =
      edgeR_perCluster_combinedResults,
    test_definitions =
      edgeR_perCluster_testDefinitions,
    output_file = combined_xlsx_file,
    include_cluster_id = TRUE,
    xlsx_max_decimal_places =
      xlsx_max_decimal_places
  )

  for (test_id in edgeR_perCluster_testDefinitions$test_id) {
    write_tsv_pseudobulk_per_cluster(
      edgeR_perCluster_combinedResults[[test_id]],
      file.path(
        combined_results_dir,
        paste0(
          analysis_prefix,
          "_allClusters_",
          test_id,
          "_fullResults.tsv.gz"
        )
      )
    )
  }

  descriptive_metric_definitions_file <- file.path(
    output_statistics_dir,
    "02_designAndSampleMetadata",
    paste0(
      analysis_prefix,
      "_descriptiveMetricDefinitions.tsv"
    )
  )

  write_tsv_pseudobulk_per_cluster(
    edgeR_perCluster_descriptiveMetricDefinitions,
    descriptive_metric_definitions_file
  )

  save(
    edgeR_perCluster_results,
    edgeR_perCluster_combinedResults,
    edgeR_perCluster_analysisSummary,
    edgeR_perCluster_clusterStatus,
    edgeR_perCluster_testDefinitions,
    edgeR_perCluster_geneAbundanceQC,
    edgeR_perCluster_combinedGeneAbundanceQC,
    edgeR_perCluster_descriptiveMetricDefinitions,
    file = result_rdata_file,
    compress = "xz"
  )

  summary_rows <- dplyr::bind_rows(
    lapply(
      edgeR_perCluster_testDefinitions$test_id,
      function(test_id) {
        qn_columns <- grep(
          "^(mean|sd)_log2plus1_QN_",
          colnames(
            edgeR_perCluster_combinedResults[[test_id]]
          ),
          value = TRUE
        )

        tibble::tibble(
          test_id = test_id,
          n_log2plus1_QN_columns = length(qn_columns),
          log2plus1_QN_columns = paste(
            qn_columns,
            collapse = " | "
          )
        )
      }
    )
  )

  summary_file <- file.path(
    output_statistics_dir,
    "02_designAndSampleMetadata",
    paste0(
      analysis_prefix,
      "_log2plus1QN_groupMeanSD_summary.tsv"
    )
  )

  write_tsv_pseudobulk_per_cluster(
    summary_rows,
    summary_file
  )

  message("\nlog2+QN descriptive metrics added successfully.")
  message("Updated combined XLSX:")
  message(combined_xlsx_file)
  message("Updated edgeR result RData:")
  message(result_rdata_file)
  message("Metric summary:")
  message(summary_file)

  list(
    combined_xlsx_file = combined_xlsx_file,
    result_rdata_file = result_rdata_file,
    descriptive_metric_definitions_file =
      descriptive_metric_definitions_file,
    summary_file = summary_file
  )
}

# ==============================================================================
# 8. Run edgeR analysis
# ==============================================================================

analysis_output <- run_pseudobulk_per_cluster_edger_two_factor_interaction(
  project_root = project_root,
  dataset_name = dataset_name,
  clustering_name = clustering_name,
  input_rdata_file = input_rdata_file,
  output_statistics_dir = output_statistics_dir,
  metadata_file = metadata_file,
  annotation_rds_file = annotation_rds_file,
  annotation_metadata_file = annotation_metadata_file,
  seurat_object_name = seurat_object_name,
  assay_name = "RNA",
  sample_column = "sample_ID",
  cluster_column = cluster_column,
  excluded_samples = c(
    "20_1F",
    "12_3F",
    "15_1M",
    "20_3M"
  ),
  expected_number_of_samples = expected_number_of_samples,
  expected_number_of_clusters = expected_number_of_clusters,
  min_spots_per_sample_cluster = min_spots_per_sample_cluster,
  min_samples_per_factor_cell = 2L,
  positive_spot_min_count = positive_spot_min_count,
  xlsx_max_decimal_places = xlsx_max_decimal_places,
  filter_by_expr_parameters = list(
    min.count = 10,
    min.total.count = 15,
    large.n = 10,
    min.prop = 0.7
  ),
  stop_on_cluster_error = TRUE
)


# ==============================================================================
# 9. Create per-cluster log2(count + 1) + quantile-normalized matrices
# ==============================================================================

log2_qn_output <- prepare_log2plus1_quantile_normalized_raw_counts(
  output_statistics_dir = output_statistics_dir,
  expected_number_of_samples = expected_number_of_samples,
  expected_number_of_clusters = expected_number_of_clusters
)


# ==============================================================================
# 10. Add log2(count + 1) + QN group mean/SD to all edgeR result tables
# ==============================================================================

log2_qn_statistics_output <- add_log2qn_descriptive_metrics_to_edger_outputs(
  output_statistics_dir = output_statistics_dir,
  log2_qn_output = log2_qn_output,
  xlsx_max_decimal_places = xlsx_max_decimal_places
)


# ==============================================================================
# 11. Main output paths
# ==============================================================================

message("\n============================================================")
message("MAIN OUTPUT PATHS")
message("============================================================")

message("\nMain combined edgeR workbook:")
message(analysis_output$combined_xlsx_file)

message("\nRData containing edgeR result tables and descriptive QC:")
message(analysis_output$result_rdata_file)

message("\nCombined abundance and detection QC table:")
message(analysis_output$combined_gene_abundance_qc_file)

message("\nPer-cluster sample-level TMM-CPM and detection metrics:")
message(analysis_output$per_cluster_sample_metrics_dir)

message("\nDefinitions of descriptive edgeR metrics:")
message(analysis_output$descriptive_metric_definitions_file)

message("\nPer-cluster log2(count + 1) + quantile-normalized matrices:")
message(log2_qn_output$output_dir)

message("\nLog2 + QN transformation summary:")
message(log2_qn_output$summary_file)

message("\nLog2 + QN definitions:")
message(log2_qn_output$definitions_file)

message("\nUpdated combined edgeR workbook with log2+QN mean/SD:")
message(log2_qn_statistics_output$combined_xlsx_file)

message("\nUpdated descriptive metric definitions:")
message(log2_qn_statistics_output$descriptive_metric_definitions_file)

message("\nSummary of log2+QN columns added per edgeR test:")
message(log2_qn_statistics_output$summary_file)


# ==============================================================================
# End
# ==============================================================================
