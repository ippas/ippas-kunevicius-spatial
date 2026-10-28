#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

SCRIPT_BUILD <- "20260831_BEFORE_AFTER_V2"

message("\n============================================================")
message("SCRIPT BUILD: ", SCRIPT_BUILD)
message("THIS VERSION CALCULATES BOTH BEFORE AND AFTER CRESCENDO")
message("============================================================")

# ==============================================================================
# 06_PERMANOVA_slideID_crescendo.R
#
# Global PERMANOVA QC for slide-associated variation:
#
#   BEFORE Crescendo
#   AFTER Crescendo
#
# plus sensitivity analysis after removing sample 13_1F.
#
# AFTER INPUT
#   Crescendo-corrected pseudobulk counts
#     -> log2(count + 1)
#     -> quantile normalization within each cluster
#
# BEFORE INPUT
#   Original raw pseudobulk counts
#     -> log2(count + 1)
#     -> quantile normalization within each cluster
#
# GLOBAL FEATURES
#   Each gene in each anatomical cluster is treated as a separate feature:
#
#     ENSMUSG...__C1
#     ENSMUSG...__C2
#     ...
#
# FEATURE SETS
#   ALL variable features
#   TOP 10000
#   TOP 5000
#   TOP 2000
#   TOP 1000
#   TOP 500
#
# IMPORTANT
#   Feature variance and feature ranking are calculated ONCE from the
#   AFTER-Crescendo n=16 matrix.
#
#   Exactly the SAME ranked gene-cluster features are then used:
#     - AFTER Crescendo, n=16
#     - AFTER Crescendo, n=15 without 13_1F
#     - BEFORE Crescendo, n=16
#     - BEFORE Crescendo, n=15 without 13_1F
#
#   This gives a strict 1:1 BEFORE vs AFTER comparison.
#
# DISTANCE
#   Euclidean distance between biological samples.
#
# PERMANOVA
#   vegan::adonis2(distance ~ slide_ID, by = "terms")
#
# OUTPUTS
#   1. AFTER Crescendo PERMANOVA table
#   2. BEFORE Crescendo PERMANOVA table
#   3. BEFORE vs AFTER comparison table
#   4. Fixed feature ranking
#   5. Slide counts by scenario
#
# INTERPRETATION
#   slide_R2:
#     fraction of multivariate expression variation associated with slide_ID.
#
#   slide_p_value:
#     permutation p-value for association with slide_ID.
#
#   R2_reduction_percent:
#     positive value = smaller slide-associated component AFTER Crescendo
#     zero           = no change
#     negative value = larger slide-associated component AFTER Crescendo
#
# IMPORTANT CAVEAT
#   slide_ID is partially associated with biological variables in this design.
#   Therefore this is a test of slide-associated structure, not proof that
#   every slide-associated component is purely technical batch.
# ==============================================================================

# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tibble",
  "readr",
  "vegan"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(vegan)
})

# ==============================================================================
# 2. Configuration
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"

dataset_name <- "maternalFMT_n16samples"

clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

excluded_sample <- "13_1F"

feature_thresholds <- c(
  Inf,
  10000L,
  5000L,
  2000L,
  1000L,
  500L
)

number_of_permutations <- 9999L

random_seed <- 20260831L

distance_method <- "euclidean"

# ==============================================================================
# 3. AFTER-Crescendo paths
# ==============================================================================

crescendo_statistics_root <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)

edgeR_statistics_dir_after <- file.path(
  crescendo_statistics_root,
  "02_edgeR"
)

after_analysis_prefix <- paste0(
  dataset_name,
  "_",
  clustering_name,
  "_crescendoCorrectedPseudobulk_edgeR_twoFactorInteraction"
)

after_log2_qn_input_dir <- file.path(
  edgeR_statistics_dir_after,
  "01_pseudobulkInputs",
  "perCluster_crescendoCorrectedCounts_log2plus1_QN"
)

crescendo_input_rds <- file.path(
  crescendo_statistics_root,
  "01_crescendoBatchCorrection",
  "07_objects",
  paste0(
    dataset_name,
    "_",
    clustering_name,
    "_crescendoBatchCorrection_pseudobulkCounts_beforeEdgeR.rds"
  )
)

# ==============================================================================
# 4. BEFORE-Crescendo paths
# ==============================================================================

before_statistics_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)

before_analysis_prefix <- paste0(
  dataset_name,
  "_",
  clustering_name,
  "_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction"
)

before_log2_qn_input_dir <- file.path(
  before_statistics_dir,
  "01_pseudobulkInputs",
  "perCluster_rawCounts_log2plus1Transformation_quantileNormalize"
)

# ==============================================================================
# 5. Output paths
# ==============================================================================

output_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,
  "02_heatmaps",
  "00_QC_residualBatchEffect",
  "03_PERMANOVA_slideID"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

after_results_tsv <- file.path(
  output_dir,
  paste0(
    "01_",
    dataset_name,
    "_AFTER_crescendo_PERMANOVA_slideID.tsv"
  )
)

before_results_tsv <- file.path(
  output_dir,
  paste0(
    "02_",
    dataset_name,
    "_BEFORE_crescendo_PERMANOVA_slideID.tsv"
  )
)

comparison_results_tsv <- file.path(
  output_dir,
  paste0(
    "03_",
    dataset_name,
    "_PERMANOVA_slideID_BEFORE_vs_AFTER_crescendo.tsv"
  )
)

feature_ranking_tsv <- file.path(
  output_dir,
  paste0(
    "04_",
    dataset_name,
    "_crescendo_geneClusterFeatureVarianceRanking.tsv.gz"
  )
)

slide_counts_tsv <- file.path(
  output_dir,
  paste0(
    "05_",
    dataset_name,
    "_PERMANOVA_slideCountsByScenario.tsv"
  )
)

# ==============================================================================
# 6. Validate inputs
# ==============================================================================

message("\n============================================================")
message("STAGE 1: VALIDATING INPUTS")
message("============================================================")

if (!file.exists(crescendo_input_rds)) {
  stop(
    "Missing Crescendo RDS:\n",
    crescendo_input_rds,
    call. = FALSE
  )
}

if (!dir.exists(after_log2_qn_input_dir)) {
  stop(
    "Missing AFTER-Crescendo log2+QN directory:\n",
    after_log2_qn_input_dir,
    call. = FALSE
  )
}

if (!dir.exists(before_log2_qn_input_dir)) {
  stop(
    "Missing BEFORE-Crescendo log2+QN directory:\n",
    before_log2_qn_input_dir,
    call. = FALSE
  )
}

message("AFTER directory:")
message(after_log2_qn_input_dir)

message("\nBEFORE directory:")
message(before_log2_qn_input_dir)

# ==============================================================================
# 7. Load sample metadata
# ==============================================================================

message("\n============================================================")
message("STAGE 2: LOADING SAMPLE METADATA")
message("============================================================")

crescendo_object <- readRDS(
  crescendo_input_rds
)

if (!"sample_metadata" %in% names(crescendo_object)) {
  stop(
    "Crescendo RDS does not contain sample_metadata.",
    call. = FALSE
  )
}

sample_metadata_raw <- tibble::as_tibble(
  crescendo_object$sample_metadata
)

required_sample_columns <- c(
  "sample_ID",
  "fmt_donor_group",
  "sex"
)

missing_sample_columns <- setdiff(
  required_sample_columns,
  colnames(sample_metadata_raw)
)

if (length(missing_sample_columns) > 0L) {
  stop(
    "sample_metadata is missing columns: ",
    paste(missing_sample_columns, collapse = ", "),
    call. = FALSE
  )
}

slide_candidates <- c(
  "slide_ID",
  "slide_id",
  "slide",
  "Slide_ID",
  "SlideID",
  "Slide"
)

slide_column <- slide_candidates[
  slide_candidates %in%
    colnames(sample_metadata_raw)
]

if (length(slide_column) == 0L) {
  stop(
    "Could not identify slide-ID column. Available columns: ",
    paste(colnames(sample_metadata_raw), collapse = ", "),
    call. = FALSE
  )
}

slide_column <- slide_column[[1L]]

sample_metadata <- sample_metadata_raw |>
  dplyr::transmute(
    sample_ID =
      as.character(.data$sample_ID),
    fmt_donor_group =
      as.character(.data$fmt_donor_group),
    sex =
      as.character(.data$sex),
    slide_ID =
      as.character(.data[[slide_column]])
  ) |>
  dplyr::distinct(
    .data$sample_ID,
    .keep_all = TRUE
  )

if (nrow(sample_metadata) != 16L) {
  stop(
    "Expected 16 samples; found ",
    nrow(sample_metadata),
    ".",
    call. = FALSE
  )
}

if (!excluded_sample %in% sample_metadata$sample_ID) {
  stop(
    "Excluded sample not found: ",
    excluded_sample,
    call. = FALSE
  )
}

if (anyNA(sample_metadata$slide_ID)) {
  stop(
    "NA slide_ID values found.",
    call. = FALSE
  )
}

sample_order <- sample_metadata$sample_ID

message("Using slide column: ", slide_column)

print(
  sample_metadata,
  n = Inf
)

# ==============================================================================
# 8. Find AFTER per-cluster matrices
# ==============================================================================

message("\n============================================================")
message("STAGE 3: FINDING AFTER-CRESCENDO MATRICES")
message("============================================================")

after_input_files <- list.files(
  after_log2_qn_input_dir,
  pattern =
    "_cluster_[0-9]+_crescendoCorrectedCounts_log2plus1_QN\\.tsv\\.gz$",
  full.names = TRUE
)

extract_after_cluster_id <- function(
    file_name
) {
  sub(
    ".*_cluster_([0-9]+)_crescendoCorrectedCounts_log2plus1_QN\\.tsv\\.gz$",
    "\\1",
    basename(file_name)
  )
}

after_cluster_ids <- vapply(
  after_input_files,
  extract_after_cluster_id,
  FUN.VALUE = character(1)
)

after_order <- order(
  as.integer(after_cluster_ids)
)

after_input_files <-
  after_input_files[
    after_order
  ]

after_cluster_ids <-
  after_cluster_ids[
    after_order
  ]

if (length(after_cluster_ids) != 16L) {
  stop(
    "Expected 16 AFTER matrices; found ",
    length(after_cluster_ids),
    ".",
    call. = FALSE
  )
}

message(
  "AFTER clusters: ",
  paste(after_cluster_ids, collapse = ", ")
)

# ==============================================================================
# 9. Build AFTER global gene-cluster matrix
# ==============================================================================

message("\n============================================================")
message("STAGE 4: BUILDING AFTER GLOBAL FEATURE MATRIX")
message("============================================================")

after_expression_list <- vector(
  "list",
  length(after_input_files)
)

after_feature_metadata_list <- vector(
  "list",
  length(after_input_files)
)

for (
  i in
    seq_along(after_input_files)
) {
  cluster_id_current <-
    after_cluster_ids[[i]]

  message(
    "[",
    i,
    "/16] Cluster ",
    cluster_id_current
  )

  x <- readr::read_tsv(
    after_input_files[[i]],
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    tibble::as_tibble()

  required_columns <- c(
    "ensembl_gene_id",
    "gene",
    "chromosome"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(x)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "AFTER cluster ",
      cluster_id_current,
      " missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  missing_samples <- setdiff(
    sample_order,
    colnames(x)
  )

  if (length(missing_samples) > 0L) {
    stop(
      "AFTER cluster ",
      cluster_id_current,
      " missing samples: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  m <- as.matrix(
    x[
      ,
      sample_order,
      drop = FALSE
    ]
  )

  storage.mode(m) <-
    "double"

  if (
    anyNA(m) ||
      any(!is.finite(m))
  ) {
    stop(
      "NA/non-finite AFTER values in cluster ",
      cluster_id_current,
      ".",
      call. = FALSE
    )
  }

  feature_id <- paste0(
    x$ensembl_gene_id,
    "__C",
    cluster_id_current
  )

  rownames(m) <-
    feature_id

  after_expression_list[[i]] <- m

  after_feature_metadata_list[[i]] <- tibble::tibble(
    feature_id =
      feature_id,
    cluster_id =
      as.character(cluster_id_current),
    ensembl_gene_id =
      as.character(x$ensembl_gene_id),
    gene =
      as.character(x$gene),
    chromosome =
      as.character(x$chromosome)
  )
}

after_expression_all <- do.call(
  rbind,
  after_expression_list
)

after_feature_metadata <- dplyr::bind_rows(
  after_feature_metadata_list
)

if (anyDuplicated(
  rownames(after_expression_all)
)) {
  stop(
    "Duplicated AFTER gene-cluster feature IDs.",
    call. = FALSE
  )
}

message(
  "AFTER matrix: ",
  format(
    nrow(after_expression_all),
    big.mark = ","
  ),
  " features x ",
  ncol(after_expression_all),
  " samples"
)

# ==============================================================================
# 10. Rank features ONCE using AFTER n=16 matrix
# ==============================================================================

message("\n============================================================")
message("STAGE 5: RANKING FEATURES ON AFTER n=16 MATRIX")
message("============================================================")

after_feature_variance <- apply(
  after_expression_all,
  1L,
  stats::var
)

valid_after_features <-
  is.finite(
    after_feature_variance
  ) &
  after_feature_variance > 0

ranked_feature_ids <- rownames(
  after_expression_all
)[
  valid_after_features
][
  order(
    after_feature_variance[
      valid_after_features
    ],
    decreasing = TRUE
  )
]

ranked_feature_variance <- after_feature_variance[
  match(
    ranked_feature_ids,
    rownames(after_expression_all)
  )
]

message(
  "Variable AFTER features: ",
  format(
    length(ranked_feature_ids),
    big.mark = ","
  )
)

feature_ranking <- tibble::tibble(
  feature_rank =
    seq_along(ranked_feature_ids),
  feature_id =
    ranked_feature_ids,
  variance_after_n16 =
    as.numeric(
      ranked_feature_variance
    )
) |>
  dplyr::left_join(
    after_feature_metadata,
    by = "feature_id"
  ) |>
  dplyr::select(
    "feature_rank",
    "feature_id",
    "cluster_id",
    "ensembl_gene_id",
    "gene",
    "chromosome",
    "variance_after_n16"
  )

readr::write_tsv(
  feature_ranking,
  feature_ranking_tsv,
  na = "NA"
)

# ==============================================================================
# 11. Find BEFORE per-cluster matrices
# ==============================================================================

message("\n============================================================")
message("STAGE 6: FINDING BEFORE-CRESCENDO MATRICES")
message("============================================================")

find_before_file_for_cluster <- function(
    cluster_id_current
) {
  candidate_files <- c(
    file.path(
      before_log2_qn_input_dir,
      paste0(
        before_analysis_prefix,
        "_cluster_",
        cluster_id_current,
        "_rawCounts_log2plus1Transformation_quantileNormalize_",
        "proteinCodingAutosomalAndX.tsv.gz"
      )
    ),
    file.path(
      before_log2_qn_input_dir,
      paste0(
        before_analysis_prefix,
        "_cluster_C",
        cluster_id_current,
        "_rawCounts_log2plus1Transformation_quantileNormalize_",
        "proteinCodingAutosomalAndX.tsv.gz"
      )
    )
  )

  existing_files <- candidate_files[
    file.exists(candidate_files)
  ]

  if (length(existing_files) != 1L) {
    stop(
      "Expected exactly one BEFORE file for cluster ",
      cluster_id_current,
      "; found ",
      length(existing_files),
      ". Checked:\n",
      paste(candidate_files, collapse = "\n"),
      call. = FALSE
    )
  }

  existing_files[[1L]]
}

before_input_files <- vapply(
  as.character(1:16),
  find_before_file_for_cluster,
  FUN.VALUE = character(1)
)

names(before_input_files) <-
  as.character(1:16)

message(
  "BEFORE clusters: ",
  paste(
    names(before_input_files),
    collapse = ", "
  )
)

# ==============================================================================
# 12. Build BEFORE global gene-cluster matrix
# ==============================================================================

message("\n============================================================")
message("STAGE 7: BUILDING BEFORE GLOBAL FEATURE MATRIX")
message("============================================================")

before_expression_list <- vector(
  "list",
  16L
)

names(before_expression_list) <-
  as.character(1:16)

for (
  cluster_id_current in
    as.character(1:16)
) {
  message(
    "[",
    cluster_id_current,
    "/16] Cluster ",
    cluster_id_current
  )

  x <- readr::read_tsv(
    before_input_files[[cluster_id_current]],
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    tibble::as_tibble()

  required_columns <- c(
    "ensembl_gene_id",
    "gene",
    "chromosome"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(x)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "BEFORE cluster ",
      cluster_id_current,
      " missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  missing_samples <- setdiff(
    sample_order,
    colnames(x)
  )

  if (length(missing_samples) > 0L) {
    stop(
      "BEFORE cluster ",
      cluster_id_current,
      " missing samples: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyDuplicated(
    x$ensembl_gene_id
  )) {
    stop(
      "Duplicated Ensembl IDs in BEFORE cluster ",
      cluster_id_current,
      ".",
      call. = FALSE
    )
  }

  m <- as.matrix(
    x[
      ,
      sample_order,
      drop = FALSE
    ]
  )

  storage.mode(m) <-
    "double"

  if (
    anyNA(m) ||
      any(!is.finite(m))
  ) {
    stop(
      "NA/non-finite BEFORE values in cluster ",
      cluster_id_current,
      ".",
      call. = FALSE
    )
  }

  feature_id <- paste0(
    x$ensembl_gene_id,
    "__C",
    cluster_id_current
  )

  rownames(m) <-
    feature_id

  before_expression_list[[cluster_id_current]] <- m
}

before_expression_all <- do.call(
  rbind,
  before_expression_list
)

if (anyDuplicated(
  rownames(before_expression_all)
)) {
  stop(
    "Duplicated BEFORE gene-cluster feature IDs.",
    call. = FALSE
  )
}

message(
  "BEFORE matrix: ",
  format(
    nrow(before_expression_all),
    big.mark = ","
  ),
  " features x ",
  ncol(before_expression_all),
  " samples"
)

# ==============================================================================
# 13. Force strict common AFTER-ranked feature universe
# ==============================================================================

message("\n============================================================")
message("STAGE 8: VALIDATING STRICT COMMON FEATURE SET")
message("============================================================")

missing_before_features <- setdiff(
  ranked_feature_ids,
  rownames(before_expression_all)
)

if (length(missing_before_features) > 0L) {
  stop(
    "Strict 1:1 comparison impossible: ",
    length(missing_before_features),
    " AFTER-ranked features are absent BEFORE. First examples: ",
    paste(
      head(
        missing_before_features,
        20
      ),
      collapse = ", "
    ),
    call. = FALSE
  )
}

after_expression_fixed <- after_expression_all[
  ranked_feature_ids,
  sample_order,
  drop = FALSE
]

before_expression_fixed <- before_expression_all[
  ranked_feature_ids,
  sample_order,
  drop = FALSE
]

if (!identical(
  rownames(after_expression_fixed),
  rownames(before_expression_fixed)
)) {
  stop(
    "BEFORE and AFTER feature order differs.",
    call. = FALSE
  )
}

if (!identical(
  colnames(after_expression_fixed),
  colnames(before_expression_fixed)
)) {
  stop(
    "BEFORE and AFTER sample order differs.",
    call. = FALSE
  )
}

message(
  "Strict common feature set: ",
  format(
    nrow(after_expression_fixed),
    big.mark = ","
  ),
  " features"
)

# ==============================================================================
# 14. Define fixed feature sets
# ==============================================================================

make_feature_set <- function(
    threshold
) {
  if (is.infinite(threshold)) {
    return(
      list(
        feature_set =
          "ALL",
        requested_top_n =
          NA_integer_,
        feature_ids =
          ranked_feature_ids
      )
    )
  }

  n_use <- min(
    as.integer(threshold),
    length(ranked_feature_ids)
  )

  list(
    feature_set =
      paste0(
        "TOP_",
        as.integer(threshold)
      ),
    requested_top_n =
      as.integer(threshold),
    feature_ids =
      ranked_feature_ids[
        seq_len(n_use)
      ]
  )
}

feature_sets <- lapply(
  feature_thresholds,
  make_feature_set
)

names(feature_sets) <- vapply(
  feature_sets,
  function(x) {
    x$feature_set
  },
  FUN.VALUE = character(1)
)

# ==============================================================================
# 15. Define sample scenarios
# ==============================================================================

sample_scenarios <- list(
  allSamples_n16 = list(
    scenario =
      "allSamples_n16",
    excluded_sample =
      NA_character_,
    sample_ids =
      sample_order
  ),

  without_13_1F_n15 = list(
    scenario =
      "without_13_1F_n15",
    excluded_sample =
      excluded_sample,
    sample_ids =
      sample_order[
        sample_order !=
          excluded_sample
      ]
  )
)

slide_counts <- dplyr::bind_rows(
  lapply(
    sample_scenarios,
    function(scenario_object) {
      sample_metadata |>
        dplyr::filter(
          .data$sample_ID %in%
            scenario_object$sample_ids
        ) |>
        dplyr::count(
          .data$slide_ID,
          name = "n_samples"
        ) |>
        dplyr::mutate(
          scenario =
            scenario_object$scenario,
          excluded_sample =
            scenario_object$excluded_sample,
          .before = 1
        )
    }
  )
)

readr::write_tsv(
  slide_counts,
  slide_counts_tsv,
  na = "NA"
)

# ==============================================================================
# 16. PERMANOVA function
# ==============================================================================

run_permanova <- function(
    expression_matrix,
    data_stage,
    scenario_object,
    feature_set_object,
    seed_offset
) {
  sample_ids_current <-
    scenario_object$sample_ids

  feature_ids_current <-
    feature_set_object$feature_ids

  expr <- expression_matrix[
    feature_ids_current,
    sample_ids_current,
    drop = FALSE
  ]

  metadata <- sample_metadata[
    match(
      sample_ids_current,
      sample_metadata$sample_ID
    ),
    ,
    drop = FALSE
  ]

  if (anyNA(
    metadata$sample_ID
  )) {
    stop(
      "Failed to match sample metadata.",
      call. = FALSE
    )
  }

  metadata$slide_ID <- droplevels(
    factor(
      metadata$slide_ID
    )
  )

  distance_object <- stats::dist(
    t(expr),
    method = distance_method
  )

  set.seed(
    random_seed +
      seed_offset
  )

  fit <- vegan::adonis2(
    distance_object ~ slide_ID,
    data = metadata,
    permutations =
      number_of_permutations,
    by = "terms"
  )

  fit_df <- as.data.frame(
    fit
  )

  fit_df$term <- rownames(
    fit_df
  )

  slide_row <- fit_df[
    fit_df$term ==
      "slide_ID",
    ,
    drop = FALSE
  ]

  residual_row <- fit_df[
    fit_df$term ==
      "Residual",
    ,
    drop = FALSE
  ]

  if (nrow(slide_row) != 1L) {
    stop(
      "Could not identify slide_ID row in adonis2 output. Terms: ",
      paste(fit_df$term, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(residual_row) != 1L) {
    stop(
      "Could not identify Residual row in adonis2 output. Terms: ",
      paste(fit_df$term, collapse = ", "),
      call. = FALSE
    )
  }

  tibble::tibble(
    data_stage =
      data_stage,

    scenario =
      scenario_object$scenario,

    excluded_sample =
      scenario_object$excluded_sample,

    n_samples =
      length(sample_ids_current),

    n_slides =
      nlevels(
        metadata$slide_ID
      ),

    feature_set =
      feature_set_object$feature_set,

    requested_top_n =
      feature_set_object$requested_top_n,

    n_features_used =
      length(feature_ids_current),

    distance =
      distance_method,

    slide_Df =
      as.numeric(
        slide_row$Df[[1]]
      ),

    slide_SumOfSqs =
      as.numeric(
        slide_row$SumOfSqs[[1]]
      ),

    slide_R2 =
      as.numeric(
        slide_row$R2[[1]]
      ),

    slide_F =
      as.numeric(
        slide_row$F[[1]]
      ),

    slide_p_value =
      as.numeric(
        slide_row$`Pr(>F)`[[1]]
      ),

    residual_Df =
      as.numeric(
        residual_row$Df[[1]]
      ),

    residual_SumOfSqs =
      as.numeric(
        residual_row$SumOfSqs[[1]]
      ),

    residual_R2 =
      as.numeric(
        residual_row$R2[[1]]
      ),

    n_permutations =
      number_of_permutations,

    feature_ranking_reference =
      "AFTER_Crescendo_variance_ranking_n16"
  )
}

# ==============================================================================
# 17. Run AFTER and BEFORE PERMANOVA
# ==============================================================================

message("\n============================================================")
message("STAGE 9: RUNNING AFTER AND BEFORE PERMANOVA")
message("============================================================")

after_results_list <- list()
before_results_list <- list()

counter <- 0L

for (
  scenario_name in
    names(sample_scenarios)
) {
  scenario_object <-
    sample_scenarios[[scenario_name]]

  message("\n------------------------------------------------------------")
  message(
    "SCENARIO: ",
    scenario_object$scenario
  )
  message("------------------------------------------------------------")

  for (
    feature_set_name in
      names(feature_sets)
  ) {
    counter <-
      counter + 1L

    feature_set_object <-
      feature_sets[[feature_set_name]]

    message(
      "[",
      counter,
      "/12] ",
      feature_set_object$feature_set,
      " | n features = ",
      format(
        length(
          feature_set_object$feature_ids
        ),
        big.mark = ","
      )
    )

    after_current <- run_permanova(
      expression_matrix =
        after_expression_fixed,
      data_stage =
        "AFTER_Crescendo",
      scenario_object =
        scenario_object,
      feature_set_object =
        feature_set_object,
      seed_offset =
        counter
    )

    before_current <- run_permanova(
      expression_matrix =
        before_expression_fixed,
      data_stage =
        "BEFORE_Crescendo",
      scenario_object =
        scenario_object,
      feature_set_object =
        feature_set_object,
      seed_offset =
        counter
    )

    after_results_list[[counter]] <- after_current


    before_results_list[[counter]] <- before_current

    message(
      "  BEFORE: R2 = ",
      formatC(
        before_current$slide_R2,
        format = "f",
        digits = 4
      ),
      " | F = ",
      formatC(
        before_current$slide_F,
        format = "f",
        digits = 4
      ),
      " | p = ",
      formatC(
        before_current$slide_p_value,
        format = "g",
        digits = 4
      )
    )

    message(
      "  AFTER : R2 = ",
      formatC(
        after_current$slide_R2,
        format = "f",
        digits = 4
      ),
      " | F = ",
      formatC(
        after_current$slide_F,
        format = "f",
        digits = 4
      ),
      " | p = ",
      formatC(
        after_current$slide_p_value,
        format = "g",
        digits = 4
      )
    )
  }
}

after_results <- dplyr::bind_rows(
  after_results_list
)

before_results <- dplyr::bind_rows(
  before_results_list
)

# ==============================================================================
# 18. Build direct BEFORE vs AFTER comparison
# ==============================================================================

message("\n============================================================")
message("STAGE 10: BUILDING BEFORE VS AFTER TABLE")
message("============================================================")

comparison_results <- before_results |>
  dplyr::select(
    "scenario",
    "excluded_sample",
    "n_samples",
    "n_slides",
    "feature_set",
    "requested_top_n",
    "n_features_used",
    slide_R2_before =
      "slide_R2",
    slide_F_before =
      "slide_F",
    slide_p_before =
      "slide_p_value"
  ) |>
  dplyr::left_join(
    after_results |>
      dplyr::select(
        "scenario",
        "feature_set",
        slide_R2_after =
          "slide_R2",
        slide_F_after =
          "slide_F",
        slide_p_after =
          "slide_p_value"
      ),
    by = c(
      "scenario",
      "feature_set"
    )
  ) |>
  dplyr::mutate(
    delta_R2_after_minus_before =
      .data$slide_R2_after -
      .data$slide_R2_before,

    R2_reduction_percent =
      100 *
      (
        .data$slide_R2_before -
          .data$slide_R2_after
      ) /
      .data$slide_R2_before,

    R2_ratio_after_to_before =
      .data$slide_R2_after /
      .data$slide_R2_before
  )

if (
  anyNA(
    comparison_results$slide_R2_after
  )
) {
  stop(
    "Some BEFORE rows failed to match AFTER rows.",
    call. = FALSE
  )
}

# ==============================================================================
# 19. Save outputs
# ==============================================================================

message("\n============================================================")
message("STAGE 11: SAVING OUTPUTS")
message("============================================================")

readr::write_tsv(
  after_results,
  after_results_tsv,
  na = "NA"
)

readr::write_tsv(
  before_results,
  before_results_tsv,
  na = "NA"
)

readr::write_tsv(
  comparison_results,
  comparison_results_tsv,
  na = "NA"
)

# ==============================================================================
# 20. Console result
# ==============================================================================

message("\n============================================================")
message("MAIN BEFORE VS AFTER RESULT")
message("============================================================")

print(
  comparison_results |>
    dplyr::select(
      "scenario",
      "feature_set",
      "n_features_used",
      "slide_R2_before",
      "slide_R2_after",
      "R2_reduction_percent",
      "slide_p_before",
      "slide_p_after"
    ),
  n = Inf
)

message("\n============================================================")
message("PERMANOVA BEFORE VS AFTER CRESCENDO FINISHED")
message("============================================================")

message("\nAFTER results:")
message(
  normalizePath(
    after_results_tsv,
    mustWork = TRUE
  )
)

message("\nBEFORE results:")
message(
  normalizePath(
    before_results_tsv,
    mustWork = TRUE
  )
)

message("\nBEFORE vs AFTER:")
message(
  normalizePath(
    comparison_results_tsv,
    mustWork = TRUE
  )
)

message("\nFeature ranking:")
message(
  normalizePath(
    feature_ranking_tsv,
    mustWork = TRUE
  )
)

message("\nSlide counts:")
message(
  normalizePath(
    slide_counts_tsv,
    mustWork = TRUE
  )
)

# ==============================================================================
# End
# ==============================================================================
