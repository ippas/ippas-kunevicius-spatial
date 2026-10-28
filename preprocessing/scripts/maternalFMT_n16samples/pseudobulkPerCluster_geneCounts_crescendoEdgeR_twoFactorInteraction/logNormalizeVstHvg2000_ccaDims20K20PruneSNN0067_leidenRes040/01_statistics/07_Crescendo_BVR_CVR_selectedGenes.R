#!/usr/bin/env Rscript

# ==============================================================================
# 07_BVR_CVR_crescendo.R
#
# Crescendo method-specific QC:
#   1. Calculate BVR and CVR only for selected genes:
#        - all OLD edgeR genes with FDR < 0.05 (no log2FC threshold)
#        - selected ribosomal genes
#        - Actb
#   2. Use raw spot counts BEFORE and Crescendo-corrected spot counts AFTER.
#   3. Use Leiden C1-C16 as the user-defined biological-cluster proxy.
#   4. Use slide_ID as batch.
#   5. Use observed raw nUMI as offset in BOTH models.
#   6. Extract the original edgeR donor-group DEGs:
#        FDR < 0.05 & |log2FC| >= 0.5
#   7. Compare the same old cluster-gene signals with edgeR AFTER Crescendo.
#   8. Summarize whether old DEGs preferentially show reduced batch variance.
#
# Crescendo paper metric:
#
#   count ~ 1 +
#     (1 | biological_cluster) +
#     (1 | batch) +
#     (1 | biological_cluster:batch) +
#     offset(log(observed_nUMI))
#
#   batch variance =
#     Var(batch) + Var(biological_cluster:batch)
#
#   BVR = batch_variance_AFTER / batch_variance_BEFORE
#   CVR = biological_cluster_variance_AFTER /
#         biological_cluster_variance_BEFORE
#
# Target zone defined in the Crescendo paper:
#
#   BVR < 1 AND CVR >= 0.5
#
# IMPORTANT FOR THIS VISIUM DATASET
#   Leiden clusters are used as the discrete biological-cluster proxy.
#   They are anatomical/spatial clusters, not purified cell types.
#
# IMPORTANT FOR OLD-vs-AFTER edgeR
#   The script separately reports old DEGs that were NOT present in the
#   Crescendo-corrected gene universe. Those genes must not be described as
#   "removed by correction", because they were excluded before correction /
#   downstream testing.
#
# Paper:
#   Batch correcting single-cell spatial transcriptomics count data with
#   Crescendo improves visualization and detection of spatial gene patterns
#   Genome Biology (2025)
# ==============================================================================

options(stringsAsFactors = FALSE)

SCRIPT_BUILD <- "20260831_BVR_CVR_SELECTED_ONE_GENE_V3"

message("\n============================================================")
message("SCRIPT BUILD: ", SCRIPT_BUILD)
message("CRESCENDO BVR/CVR — SELECTED GENES, ONE GENE AT A TIME")
message("============================================================")

set.seed(20260831)

# ==============================================================================
# 1. Configuration
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"

dataset_name <- "maternalFMT_n16samples"

clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

cluster_column <- "leiden_res0.40"
sample_column <- "sample_ID"
assay_name <- "RNA"

main_test_id <- "Overall_Group_ASD_vs_Neurotypical"

de_fdr_threshold <- 0.05

# Strict threshold retained for the separate OLD-vs-AFTER DEG survival summary.
de_abs_log2fc_threshold <- 0.5

# BVR/CVR pilot selection:
# all OLD genes significant at FDR < 0.05, WITHOUT a log2FC threshold.
control_genes <- c(
  "Rpl39",
  "Rpl36",
  "Rpl38",
  "Rpl37",
  "Rpl29",
  "Rpl13",
  "Rps27",
  "Rps28",
  "Rps29",
  "Rps3a1",
  "Actb"
)

# Process exactly one selected gene at a time.
# BEFORE and AFTER models are fitted sequentially for that gene.
process_one_gene_at_a_time <- TRUE

# The Crescendo paper notes that very small pre-correction batch variance makes
# BVR less informative. In one real-data benchmark, genes with batch variance
# > 0.001 were emphasized.
minimum_pre_batch_variance_for_evaluation <- 0.001

# Full run by default.
test_mode <- FALSE
test_number_of_genes <- 20L

# Resume finished individual-gene results.
resume_existing_genes <- TRUE

# One-gene-at-a-time mode is intentionally sequential so the console always
# shows exactly which gene is being calculated and each completed gene can be
# saved immediately.

# lme4::glmer settings.
# nAGQ = 1 is the standard Laplace approximation.
glmer_nAGQ <- 1L
glmer_maxfun <- 100000L

# ==============================================================================
# 2. Input paths
# ==============================================================================

input_seurat_rdata_file <- file.path(
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

crescendo_statistics_root <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)

crescendo_root <- file.path(
  crescendo_statistics_root,
  "01_crescendoBatchCorrection"
)

crescendo_input_rds <- file.path(
  crescendo_root,
  "07_objects",
  paste0(
    dataset_name,
    "_",
    clustering_name,
    "_crescendoBatchCorrection_pseudobulkCounts_beforeEdgeR.rds"
  )
)

corrected_spot_chunks_dir <- file.path(
  crescendo_root,
  "03_correctedSpotCounts_chunks_v13"
)

old_edger_rdata_file <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics",
  "03_edgeRResults",
  paste0(
    dataset_name,
    "_",
    clustering_name,
    "_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction_",
    "edgeRResults.RData"
  )
)

after_edger_results_dir <- file.path(
  crescendo_statistics_root,
  "02_edgeR",
  "03_edgeRResults"
)

# ==============================================================================
# 3. Output paths
# ==============================================================================

output_root <- file.path(
  crescendo_root,
  "08_BVR_CVR"
)

output_chunks_dir <- file.path(
  output_root,
  "01_geneResults"
)

output_tables_dir <- file.path(
  output_root,
  "02_tables"
)

output_plots_dir <- file.path(
  output_root,
  "03_plots"
)

output_logs_dir <- file.path(
  output_root,
  "04_logs"
)

for (dir_current in c(
  output_root,
  output_chunks_dir,
  output_tables_dir,
  output_plots_dir,
  output_logs_dir
)) {
  dir.create(
    dir_current,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

all_gene_metrics_tsv <- file.path(
  output_tables_dir,
  "01_selectedGenes_BVR_CVR.tsv.gz"
)

selected_gene_list_tsv <- file.path(
  output_tables_dir,
  "00_selectedGenes_for_BVR_CVR.tsv"
)

global_summary_tsv <- file.path(
  output_tables_dir,
  "02_BVR_CVR_summary.tsv"
)

old_deg_genes_tsv <- file.path(
  output_tables_dir,
  "03_oldEdgeR_DEG_FDR005_absLog2FC05_uniqueGenes.tsv"
)

old_deg_cluster_gene_tsv <- file.path(
  output_tables_dir,
  "04_oldEdgeR_DEG_FDR005_absLog2FC05_clusterGene.tsv.gz"
)

after_deg_genes_tsv <- file.path(
  output_tables_dir,
  "05_afterCrescendoEdgeR_DEG_FDR005_absLog2FC05_uniqueGenes.tsv"
)

after_deg_cluster_gene_tsv <- file.path(
  output_tables_dir,
  "06_afterCrescendoEdgeR_DEG_FDR005_absLog2FC05_clusterGene.tsv.gz"
)

old_vs_after_cluster_gene_tsv <- file.path(
  output_tables_dir,
  "07_oldDEG_clusterGene_BEFORE_vs_AFTER.tsv.gz"
)

old_gene_survival_tsv <- file.path(
  output_tables_dir,
  "08_oldDEG_uniqueGenes_survivalAndBVR_CVR.tsv"
)

deg_count_summary_tsv <- file.path(
  output_tables_dir,
  "09_DEG_countSummary_BEFORE_vs_AFTER.tsv"
)

deg_category_bvr_summary_tsv <- file.path(
  output_tables_dir,
  "10_BVR_CVR_by_DEGcategory.tsv"
)

model_status_tsv <- file.path(
  output_tables_dir,
  "11_BVR_CVR_modelStatusSummary.tsv"
)

all_genes_pdf <- file.path(
  output_plots_dir,
  "01_BVR_CVR_allGenes.pdf"
)

old_deg_pdf <- file.path(
  output_plots_dir,
  "02_BVR_CVR_oldDEGs.pdf"
)

# ==============================================================================
# 4. Packages
# ==============================================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "lme4",
  "dplyr",
  "tibble",
  "readr",
  "ggplot2"
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
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(lme4)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

# ==============================================================================
# 5. Helpers
# ==============================================================================

stage_message <- function(stage_number, text) {
  message("\n============================================================")
  message("STAGE ", stage_number, ": ", text)
  message("============================================================")
}

natural_order <- function(x) {
  x <- unique(as.character(x))
  numeric_x <- suppressWarnings(as.numeric(x))

  if (!anyNA(numeric_x)) {
    return(x[order(numeric_x)])
  }

  sort(x)
}

strip_ensembl_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

safe_character <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

collapse_unique_natural <- function(x) {
  paste(
    natural_order(unique(as.character(x))),
    collapse = ","
  )
}

find_single_seurat_object <- function(
    rdata_file,
    required_cluster_column
) {
  env <- new.env(parent = emptyenv())

  loaded_names <- load(
    rdata_file,
    envir = env
  )

  seurat_names <- loaded_names[
    vapply(
      loaded_names,
      function(object_name) {
        inherits(
          get(
            object_name,
            envir = env,
            inherits = FALSE
          ),
          "Seurat"
        )
      },
      FUN.VALUE = logical(1)
    )
  ]

  if (length(seurat_names) == 0L) {
    stop(
      "No Seurat object found in:\n",
      rdata_file,
      call. = FALSE
    )
  }

  with_cluster <- seurat_names[
    vapply(
      seurat_names,
      function(object_name) {
        object_current <- get(
          object_name,
          envir = env,
          inherits = FALSE
        )

        required_cluster_column %in%
          colnames(object_current@meta.data)
      },
      FUN.VALUE = logical(1)
    )
  ]

  if (length(with_cluster) == 1L) {
    chosen <- with_cluster[[1L]]
  } else if (length(with_cluster) > 1L) {
    object_sizes <- vapply(
      with_cluster,
      function(object_name) {
        ncol(
          get(
            object_name,
            envir = env,
            inherits = FALSE
          )
        )
      },
      FUN.VALUE = integer(1)
    )

    chosen <- with_cluster[[which.max(object_sizes)]]
  } else if (length(seurat_names) == 1L) {
    chosen <- seurat_names[[1L]]
  } else {
    stop(
      "Could not uniquely select a Seurat object.",
      call. = FALSE
    )
  }

  list(
    object_name = chosen,
    object = get(
      chosen,
      envir = env,
      inherits = FALSE
    )
  )
}

get_raw_counts <- function(
    seurat_object,
    assay_name
) {
  if (!assay_name %in% names(seurat_object@assays)) {
    stop(
      "Assay missing: ",
      assay_name,
      call. = FALSE
    )
  }

  assay_object <- seurat_object[[assay_name]]

  if (inherits(assay_object, "Assay5")) {
    all_layers <- SeuratObject::Layers(
      assay_object
    )

    count_layers <- all_layers[
      all_layers == "counts" |
        startsWith(all_layers, "counts.")
    ]

    if (length(count_layers) == 0L) {
      stop(
        "No raw counts layer found.",
        call. = FALSE
      )
    }

    if (length(count_layers) > 1L) {
      message(
        "Joining ",
        length(count_layers),
        " split counts.* layers..."
      )

      assay_object <- SeuratObject::JoinLayers(
        assay_object
      )
    }

    counts <- SeuratObject::LayerData(
      assay_object,
      layer = "counts"
    )

    return(
      methods::as(
        counts,
        "dgCMatrix"
      )
    )
  }

  counts <- Seurat::GetAssayData(
    seurat_object,
    assay = assay_name,
    slot = "counts"
  )

  methods::as(
    counts,
    "dgCMatrix"
  )
}

find_after_edger_rdata <- function(
    directory
) {
  if (!dir.exists(directory)) {
    stop(
      "AFTER edgeR directory missing:\n",
      directory,
      call. = FALSE
    )
  }

  candidates <- list.files(
    directory,
    pattern = "_edgeRResults\\.RData$",
    full.names = TRUE
  )

  if (length(candidates) == 0L) {
    stop(
      "No AFTER edgeRResults RData found in:\n",
      directory,
      call. = FALSE
    )
  }

  preferred <- candidates[
    grepl(
      "crescendoCorrected",
      basename(candidates),
      ignore.case = TRUE
    )
  ]

  if (length(preferred) == 1L) {
    return(preferred[[1L]])
  }

  if (length(candidates) == 1L) {
    return(candidates[[1L]])
  }

  stop(
    "Multiple AFTER edgeRResults RData files found:\n",
    paste(candidates, collapse = "\n"),
    call. = FALSE
  )
}

load_main_edger_result <- function(
    rdata_file,
    test_id
) {
  env <- new.env(parent = emptyenv())

  loaded_names <- load(
    rdata_file,
    envir = env
  )

  if (
    "edgeR_perCluster_combinedResults" %in%
      loaded_names
  ) {
    combined <- get(
      "edgeR_perCluster_combinedResults",
      envir = env,
      inherits = FALSE
    )
  } else {
    candidate_names <- loaded_names[
      vapply(
        loaded_names,
        function(object_name) {
          object_current <- get(
            object_name,
            envir = env,
            inherits = FALSE
          )

          is.list(object_current) &&
            test_id %in% names(object_current)
        },
        FUN.VALUE = logical(1)
      )
    ]

    if (length(candidate_names) != 1L) {
      stop(
        "Could not identify combined edgeR result object in:\n",
        rdata_file,
        call. = FALSE
      )
    }

    combined <- get(
      candidate_names[[1L]],
      envir = env,
      inherits = FALSE
    )
  }

  if (!test_id %in% names(combined)) {
    stop(
      "Test not found: ",
      test_id,
      call. = FALSE
    )
  }

  result <- tibble::as_tibble(
    combined[[test_id]]
  )

  required_columns <- c(
    "cluster_id",
    "ensembl_gene_id",
    "gene",
    "logFC",
    "PValue",
    "FDR"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(result)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "edgeR table missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  result |>
    dplyr::mutate(
      cluster_id = as.character(.data$cluster_id),
      ensembl_gene_id =
        strip_ensembl_version(.data$ensembl_gene_id),
      gene = as.character(.data$gene)
    )
}

select_deg_cluster_gene <- function(
    result_table
) {
  result_table |>
    dplyr::filter(
      is.finite(.data$FDR),
      is.finite(.data$logFC),
      .data$FDR < de_fdr_threshold,
      abs(.data$logFC) >=
        de_abs_log2fc_threshold
    ) |>
    dplyr::arrange(
      .data$FDR,
      dplyr::desc(abs(.data$logFC))
    )
}

select_fdr_only_cluster_gene <- function(
    result_table
) {
  result_table |>
    dplyr::filter(
      is.finite(.data$FDR),
      .data$FDR < de_fdr_threshold
    ) |>
    dplyr::arrange(
      .data$FDR,
      dplyr::desc(abs(.data$logFC))
    )
}

summarize_deg_genes <- function(
    significant_cluster_gene
) {
  significant_cluster_gene |>
    dplyr::group_by(
      .data$ensembl_gene_id,
      .data$gene
    ) |>
    dplyr::summarise(
      n_significant_clusters =
        dplyr::n_distinct(.data$cluster_id),
      significant_clusters =
        collapse_unique_natural(.data$cluster_id),
      minimum_FDR =
        min(.data$FDR, na.rm = TRUE),
      maximum_abs_log2FC =
        max(abs(.data$logFC), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$minimum_FDR,
      dplyr::desc(.data$maximum_abs_log2FC)
    )
}

make_feature_annotation <- function(
    feature_ids,
    old_results,
    after_results
) {
  lookup <- dplyr::bind_rows(
    old_results |>
      dplyr::select(
        "ensembl_gene_id",
        "gene"
      ),
    after_results |>
      dplyr::select(
        "ensembl_gene_id",
        "gene"
      )
  ) |>
    dplyr::mutate(
      ensembl_gene_id =
        strip_ensembl_version(.data$ensembl_gene_id),
      gene = as.character(.data$gene),
      gene_lower =
        tolower(trimws(.data$gene))
    ) |>
    dplyr::filter(
      !is.na(.data$ensembl_gene_id),
      .data$ensembl_gene_id != ""
    ) |>
    dplyr::distinct(
      .data$ensembl_gene_id,
      .data$gene,
      .data$gene_lower
    )

  ensembl_lookup <- lookup |>
    dplyr::group_by(
      .data$ensembl_gene_id
    ) |>
    dplyr::slice(1L) |>
    dplyr::ungroup()

  symbol_lookup <- lookup |>
    dplyr::filter(
      !is.na(.data$gene_lower),
      .data$gene_lower != ""
    ) |>
    dplyr::add_count(
      .data$gene_lower,
      name = "symbol_n"
    ) |>
    dplyr::filter(
      .data$symbol_n == 1L
    ) |>
    dplyr::select(
      -"symbol_n"
    )

  output <- tibble::tibble(
    feature_id =
      as.character(feature_ids),
    feature_canonical =
      strip_ensembl_version(feature_ids),
    feature_lower =
      tolower(trimws(as.character(feature_ids)))
  )

  ensembl_index <- match(
    output$feature_canonical,
    ensembl_lookup$ensembl_gene_id
  )

  symbol_index <- match(
    output$feature_lower,
    symbol_lookup$gene_lower
  )

  output$ensembl_gene_id <- NA_character_
  output$gene <- NA_character_
  output$edgeR_mapping_type <- NA_character_

  exact_ensembl <- !is.na(
    ensembl_index
  )

  output$ensembl_gene_id[
    exact_ensembl
  ] <- ensembl_lookup$ensembl_gene_id[
    ensembl_index[
      exact_ensembl
    ]
  ]

  output$gene[
    exact_ensembl
  ] <- ensembl_lookup$gene[
    ensembl_index[
      exact_ensembl
    ]
  ]

  output$edgeR_mapping_type[
    exact_ensembl
  ] <- "ensembl"

  symbol_only <- !exact_ensembl &
    !is.na(symbol_index)

  output$ensembl_gene_id[
    symbol_only
  ] <- symbol_lookup$ensembl_gene_id[
    symbol_index[
      symbol_only
    ]
  ]

  output$gene[
    symbol_only
  ] <- symbol_lookup$gene[
    symbol_index[
      symbol_only
    ]
  ]

  output$edgeR_mapping_type[
    symbol_only
  ] <- "gene_symbol_unique"

  output
}

is_gene_in_summary <- function(
    feature_annotation,
    gene_summary
) {
  summary_ensembl <- unique(
    strip_ensembl_version(
      gene_summary$ensembl_gene_id
    )
  )

  summary_symbols <- unique(
    tolower(
      trimws(
        safe_character(
          gene_summary$gene
        )
      )
    )
  )

  (
    !is.na(
      feature_annotation$ensembl_gene_id
    ) &
      feature_annotation$ensembl_gene_id %in%
        summary_ensembl
  ) |
    (
      !is.na(
        feature_annotation$gene
      ) &
        tolower(
          trimws(
            feature_annotation$gene
          )
        ) %in%
          summary_symbols
    )
}

annotate_deg_summary_with_crescendo_universe <- function(
    deg_summary,
    feature_annotation,
    after_deg_summary = NULL
) {
  crescendo_ensembl <- unique(
    feature_annotation$ensembl_gene_id[
      !is.na(
        feature_annotation$ensembl_gene_id
      )
    ]
  )

  crescendo_symbols <- unique(
    tolower(
      trimws(
        safe_character(
          feature_annotation$gene
        )
      )
    )
  )

  output <- deg_summary |>
    dplyr::mutate(
      in_crescendo_corrected_gene_universe =
        strip_ensembl_version(
          .data$ensembl_gene_id
        ) %in%
          crescendo_ensembl |
          tolower(
            trimws(
              safe_character(
                .data$gene
              )
            )
          ) %in%
            crescendo_symbols
    )

  if (!is.null(after_deg_summary)) {
    after_ensembl <- unique(
      strip_ensembl_version(
        after_deg_summary$ensembl_gene_id
      )
    )

    after_symbols <- unique(
      tolower(
        trimws(
          safe_character(
            after_deg_summary$gene
          )
        )
      )
    )

    output <- output |>
      dplyr::mutate(
        significant_after_crescendo =
          strip_ensembl_version(
            .data$ensembl_gene_id
          ) %in%
            after_ensembl |
            tolower(
              trimws(
                safe_character(
                  .data$gene
                )
              )
            ) %in%
              after_symbols,
        DEG_survival_status =
          dplyr::case_when(
            !.data$in_crescendo_corrected_gene_universe ~
              "not_in_crescendo_gene_universe",
            .data$significant_after_crescendo ~
              "retained_after_crescendo",
            TRUE ~
              "lost_after_crescendo"
          )
      )
  }

  output
}

extract_variance_component <- function(
    variance_table,
    group_name
) {
  values <- variance_table$vcov[
    variance_table$grp ==
      group_name
  ]

  if (length(values) == 0L) {
    return(NA_real_)
  }

  sum(
    as.numeric(values),
    na.rm = TRUE
  )
}

fit_variance_model <- function(
    counts_vector,
    model_metadata
) {
  model_data <- model_metadata

  model_data$count <- as.numeric(
    counts_vector
  )

  warning_messages <- character(0)

  fit_once <- function(optimizer_name) {
    withCallingHandlers(
      lme4::glmer(
        count ~
          1 +
          offset(offset_log_nUMI) +
          (1 | biological_cluster) +
          (1 | slide_ID) +
          (1 | cluster_slide),
        data = model_data,
        family = stats::poisson(
          link = "log"
        ),
        nAGQ = glmer_nAGQ,
        control = lme4::glmerControl(
          optimizer = optimizer_name,
          calc.derivs = FALSE,
          optCtrl = list(
            maxfun = glmer_maxfun
          )
        )
      ),
      warning = function(w) {
        warning_messages <<- c(
          warning_messages,
          conditionMessage(w)
        )

        invokeRestart(
          "muffleWarning"
        )
      }
    )
  }

  fit <- tryCatch(
    fit_once("bobyqa"),
    error = function(e) {
      NULL
    }
  )

  optimizer_used <- "bobyqa"

  if (is.null(fit)) {
    warning_messages <- c(
      warning_messages,
      "bobyqa failed; retried with nloptwrap"
    )

    fit <- tryCatch(
      fit_once("nloptwrap"),
      error = function(e) {
        structure(
          list(
            message =
              conditionMessage(e)
          ),
          class = "captured_fit_error"
        )
      }
    )

    optimizer_used <- "nloptwrap"
  }

  if (
    inherits(
      fit,
      "captured_fit_error"
    )
  ) {
    return(
      list(
        fit_ok = FALSE,
        optimizer = optimizer_used,
        singular = NA,
        biological_cluster_variance = NA_real_,
        slide_variance = NA_real_,
        cluster_slide_variance = NA_real_,
        batch_variance = NA_real_,
        warning_message =
          paste(
            unique(warning_messages),
            collapse = " | "
          ),
        error_message =
          fit$message
      )
    )
  }

  variance_table <- as.data.frame(
    lme4::VarCorr(fit)
  )

  biological_cluster_variance <-
    extract_variance_component(
      variance_table,
      "biological_cluster"
    )

  slide_variance <-
    extract_variance_component(
      variance_table,
      "slide_ID"
    )

  cluster_slide_variance <-
    extract_variance_component(
      variance_table,
      "cluster_slide"
    )

  batch_variance <-
    slide_variance +
      cluster_slide_variance

  singular_status <- tryCatch(
    lme4::isSingular(
      fit,
      tol = 1e-4
    ),
    error = function(e) {
      NA
    }
  )

  list(
    fit_ok =
      all(
        is.finite(
          c(
            biological_cluster_variance,
            slide_variance,
            cluster_slide_variance,
            batch_variance
          )
        )
      ),
    optimizer =
      optimizer_used,
    singular =
      singular_status,
    biological_cluster_variance =
      biological_cluster_variance,
    slide_variance =
      slide_variance,
    cluster_slide_variance =
      cluster_slide_variance,
    batch_variance =
      batch_variance,
    warning_message =
      paste(
        unique(warning_messages),
        collapse = " | "
      ),
    error_message =
      ""
  )
}

fit_one_gene_bvr_cvr <- function(
    gene_index,
    raw_chunk,
    corrected_chunk,
    model_metadata,
    verbose = TRUE
) {
  feature_id <- rownames(
    corrected_chunk
  )[[gene_index]]

  raw_counts_gene <- as.numeric(
    raw_chunk[
      gene_index,
      ,
      drop = TRUE
    ]
  )

  corrected_counts_gene <- as.numeric(
    corrected_chunk[
      gene_index,
      ,
      drop = TRUE
    ]
  )

  if (isTRUE(verbose)) {
    message("  BEFORE model...")
  }

  before_fit <- fit_variance_model(
    counts_vector =
      raw_counts_gene,
    model_metadata =
      model_metadata
  )

  if (isTRUE(verbose)) {
    message("  AFTER model...")
  }

  after_fit <- fit_variance_model(
    counts_vector =
      corrected_counts_gene,
    model_metadata =
      model_metadata
  )

  bvr <- if (
    isTRUE(before_fit$fit_ok) &&
      isTRUE(after_fit$fit_ok) &&
      is.finite(
        before_fit$batch_variance
      ) &&
      before_fit$batch_variance > 0
  ) {
    after_fit$batch_variance /
      before_fit$batch_variance
  } else {
    NA_real_
  }

  cvr <- if (
    isTRUE(before_fit$fit_ok) &&
      isTRUE(after_fit$fit_ok) &&
      is.finite(
        before_fit$biological_cluster_variance
      ) &&
      before_fit$biological_cluster_variance > 0
  ) {
    after_fit$biological_cluster_variance /
      before_fit$biological_cluster_variance
  } else {
    NA_real_
  }

  tibble::tibble(
    feature_id =
      feature_id,

    total_raw_counts =
      sum(raw_counts_gene),

    total_corrected_counts =
      sum(corrected_counts_gene),

    detected_spots_before =
      sum(raw_counts_gene > 0),

    detected_spots_after =
      sum(corrected_counts_gene > 0),

    fit_ok_before =
      before_fit$fit_ok,

    fit_ok_after =
      after_fit$fit_ok,

    singular_before =
      before_fit$singular,

    singular_after =
      after_fit$singular,

    optimizer_before =
      before_fit$optimizer,

    optimizer_after =
      after_fit$optimizer,

    biological_cluster_variance_before =
      before_fit$biological_cluster_variance,

    biological_cluster_variance_after =
      after_fit$biological_cluster_variance,

    slide_variance_before =
      before_fit$slide_variance,

    slide_variance_after =
      after_fit$slide_variance,

    cluster_slide_variance_before =
      before_fit$cluster_slide_variance,

    cluster_slide_variance_after =
      after_fit$cluster_slide_variance,

    batch_variance_before =
      before_fit$batch_variance,

    batch_variance_after =
      after_fit$batch_variance,

    BVR =
      bvr,

    CVR =
      cvr,

    metric_valid =
      is.finite(bvr) &
        is.finite(cvr),

    pre_batch_variance_ge_0p001 =
      is.finite(
        before_fit$batch_variance
      ) &
        before_fit$batch_variance >=
          minimum_pre_batch_variance_for_evaluation,

    BVR_lt_1 =
      is.finite(bvr) &
        bvr < 1,

    CVR_ge_0p5 =
      is.finite(cvr) &
        cvr >= 0.5,

    CVR_ge_1 =
      is.finite(cvr) &
        cvr >= 1,

    target_zone =
      is.finite(bvr) &
        is.finite(cvr) &
        bvr < 1 &
        cvr >= 0.5,

    warning_before =
      before_fit$warning_message,

    warning_after =
      after_fit$warning_message,

    error_before =
      before_fit$error_message,

    error_after =
      after_fit$error_message
  )
}

summarize_metric_subset <- function(
    data,
    set_name,
    evaluation_scope
) {
  if (
    identical(
      evaluation_scope,
      "meaningful_pre_batch_variance"
    )
  ) {
    data <- data |>
      dplyr::filter(
        .data$pre_batch_variance_ge_0p001
      )
  }

  valid <- data |>
    dplyr::filter(
      .data$metric_valid
    )

  tibble::tibble(
    gene_set =
      set_name,
    evaluation_scope =
      evaluation_scope,
    n_genes_total =
      nrow(data),
    n_genes_valid_metrics =
      nrow(valid),
    median_batch_variance_before =
      if (nrow(valid) > 0L) {
        stats::median(
          valid$batch_variance_before,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    median_batch_variance_after =
      if (nrow(valid) > 0L) {
        stats::median(
          valid$batch_variance_after,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    median_BVR =
      if (nrow(valid) > 0L) {
        stats::median(
          valid$BVR,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    median_CVR =
      if (nrow(valid) > 0L) {
        stats::median(
          valid$CVR,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    percent_BVR_lt_1 =
      if (nrow(valid) > 0L) {
        100 *
          mean(
            valid$BVR < 1,
            na.rm = TRUE
          )
      } else {
        NA_real_
      },
    percent_CVR_ge_0p5 =
      if (nrow(valid) > 0L) {
        100 *
          mean(
            valid$CVR >= 0.5,
            na.rm = TRUE
          )
      } else {
        NA_real_
      },
    percent_CVR_ge_1 =
      if (nrow(valid) > 0L) {
        100 *
          mean(
            valid$CVR >= 1,
            na.rm = TRUE
          )
      } else {
        NA_real_
      },
    percent_target_zone =
      if (nrow(valid) > 0L) {
        100 *
          mean(
            valid$BVR < 1 &
              valid$CVR >= 0.5,
            na.rm = TRUE
          )
      } else {
        NA_real_
      }
  )
}

# ==============================================================================
# 6. Validate files
# ==============================================================================

stage_message(
  1,
  "VALIDATING INPUTS"
)

required_files <- c(
  input_seurat_rdata_file,
  crescendo_input_rds,
  old_edger_rdata_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Missing required files:\n",
    paste(
      missing_files,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

if (!dir.exists(corrected_spot_chunks_dir)) {
  stop(
    "Missing corrected spot chunk directory:\n",
    corrected_spot_chunks_dir,
    call. = FALSE
  )
}

after_edger_rdata_file <-
  find_after_edger_rdata(
    after_edger_results_dir
  )

message("Old edgeR:")
message(old_edger_rdata_file)

message("\nAfter-Crescendo edgeR:")
message(after_edger_rdata_file)

message("\nCorrected spot chunks:")
message(corrected_spot_chunks_dir)

message("\nProcessing mode: one selected gene at a time")
message("glmer nAGQ: ", glmer_nAGQ)

# ==============================================================================
# 7. Load Crescendo object and edgeR results
# ==============================================================================

stage_message(
  2,
  "LOADING CRESCENDO OBJECT AND edgeR RESULTS"
)

crescendo_object <- readRDS(
  crescendo_input_rds
)

required_crescendo_elements <- c(
  "sample_metadata",
  "gene_metadata"
)

missing_crescendo_elements <- setdiff(
  required_crescendo_elements,
  names(crescendo_object)
)

if (length(missing_crescendo_elements) > 0L) {
  stop(
    "Crescendo RDS missing elements: ",
    paste(
      missing_crescendo_elements,
      collapse = ", "
    ),
    call. = FALSE
  )
}

sample_metadata <- tibble::as_tibble(
  crescendo_object$sample_metadata
)

gene_metadata <- tibble::as_tibble(
  crescendo_object$gene_metadata
)

if (!"feature_id" %in% colnames(gene_metadata)) {
  stop(
    "Crescendo gene_metadata lacks feature_id.",
    call. = FALSE
  )
}

crescendo_gene_ids <- as.character(
  gene_metadata$feature_id
)

if (anyDuplicated(crescendo_gene_ids)) {
  stop(
    "Duplicated Crescendo feature IDs.",
    call. = FALSE
  )
}

message(
  "Crescendo-corrected genes: ",
  format(
    length(crescendo_gene_ids),
    big.mark = ","
  )
)

old_main_results <- load_main_edger_result(
  old_edger_rdata_file,
  main_test_id
)

after_main_results <- load_main_edger_result(
  after_edger_rdata_file,
  main_test_id
)

old_sig_cluster_gene <- select_deg_cluster_gene(
  old_main_results
)

after_sig_cluster_gene <- select_deg_cluster_gene(
  after_main_results
)

old_deg_genes <- summarize_deg_genes(
  old_sig_cluster_gene
)

after_deg_genes <- summarize_deg_genes(
  after_sig_cluster_gene
)

# Genes selected for the BVR/CVR pilot:
# OLD edgeR FDR < 0.05, regardless of log2FC.
old_fdr005_cluster_gene <- select_fdr_only_cluster_gene(
  old_main_results
)

old_fdr005_genes <- old_fdr005_cluster_gene |>
  dplyr::group_by(
    .data$ensembl_gene_id,
    .data$gene
  ) |>
  dplyr::summarise(
    minimum_FDR =
      min(.data$FDR, na.rm = TRUE),
    maximum_abs_log2FC =
      max(abs(.data$logFC), na.rm = TRUE),
    n_FDR005_clusters =
      dplyr::n_distinct(.data$cluster_id),
    FDR005_clusters =
      collapse_unique_natural(.data$cluster_id),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    .data$minimum_FDR,
    dplyr::desc(.data$maximum_abs_log2FC)
  )

message(
  "\nOLD edgeR significant cluster-gene results: ",
  format(
    nrow(old_sig_cluster_gene),
    big.mark = ","
  )
)

message(
  "OLD edgeR unique genes: ",
  format(
    nrow(old_deg_genes),
    big.mark = ","
  )
)

message(
  "\nAFTER Crescendo significant cluster-gene results: ",
  format(
    nrow(after_sig_cluster_gene),
    big.mark = ","
  )
)

message(
  "AFTER Crescendo unique genes: ",
  format(
    nrow(after_deg_genes),
    big.mark = ","
  )
)

message(
  "\nOLD FDR < 0.05 unique genes selected for BVR/CVR: ",
  format(
    nrow(old_fdr005_genes),
    big.mark = ","
  )
)

message(
  "Additional control genes requested: ",
  paste(
    control_genes,
    collapse = ", "
  )
)

# ==============================================================================
# 8. Map Crescendo features to edgeR identifiers
# ==============================================================================

stage_message(
  3,
  "MAPPING CRESCENDO FEATURES TO edgeR GENE IDS"
)

feature_annotation <- make_feature_annotation(
  feature_ids =
    crescendo_gene_ids,
  old_results =
    old_main_results,
  after_results =
    after_main_results
)

feature_annotation$old_DEG <-
  is_gene_in_summary(
    feature_annotation,
    old_deg_genes
  )

feature_annotation$after_DEG <-
  is_gene_in_summary(
    feature_annotation,
    after_deg_genes
  )

feature_annotation$DEG_category <-
  dplyr::case_when(
    feature_annotation$old_DEG &
      feature_annotation$after_DEG ~
      "old_retained_after",
    feature_annotation$old_DEG &
      !feature_annotation$after_DEG ~
      "old_lost_after",
    !feature_annotation$old_DEG &
      feature_annotation$after_DEG ~
      "new_after_only",
    TRUE ~
      "never_DEG"
  )

old_fdr005_ensembl <- unique(
  strip_ensembl_version(
    old_fdr005_genes$ensembl_gene_id
  )
)

old_fdr005_symbols <- unique(
  tolower(
    trimws(
      safe_character(
        old_fdr005_genes$gene
      )
    )
  )
)

control_symbols_lower <- tolower(
  control_genes
)

feature_annotation$selected_old_FDR005 <-
  (
    !is.na(feature_annotation$ensembl_gene_id) &
      feature_annotation$ensembl_gene_id %in%
        old_fdr005_ensembl
  ) |
    (
      !is.na(feature_annotation$gene) &
        tolower(feature_annotation$gene) %in%
          old_fdr005_symbols
    )

feature_annotation$selected_control_gene <-
  (
    !is.na(feature_annotation$gene) &
      tolower(feature_annotation$gene) %in%
        control_symbols_lower
  ) |
    (
      tolower(feature_annotation$feature_id) %in%
        control_symbols_lower
    )

feature_annotation$selected_for_BVR_CVR <-
  feature_annotation$selected_old_FDR005 |
    feature_annotation$selected_control_gene

feature_annotation$selection_reason <-
  dplyr::case_when(
    feature_annotation$selected_old_FDR005 &
      feature_annotation$selected_control_gene ~
      "old_FDR005_and_control",
    feature_annotation$selected_old_FDR005 ~
      "old_FDR005",
    feature_annotation$selected_control_gene ~
      "control_gene",
    TRUE ~
      "not_selected"
  )

selected_feature_annotation <- feature_annotation |>
  dplyr::filter(
    .data$selected_for_BVR_CVR
  )

selected_crescendo_gene_ids <-
  selected_feature_annotation$feature_id

if (length(selected_crescendo_gene_ids) == 0L) {
  stop(
    "No selected genes were found in the Crescendo-corrected gene universe.",
    call. = FALSE
  )
}

message(
  "Features mapped to edgeR IDs: ",
  sum(
    !is.na(
      feature_annotation$ensembl_gene_id
    )
  ),
  "/",
  nrow(feature_annotation)
)

message(
  "\nSelected genes available in Crescendo universe: ",
  length(selected_crescendo_gene_ids)
)

selected_control_found <- selected_feature_annotation |>
  dplyr::filter(
    .data$selected_control_gene
  ) |>
  dplyr::pull(
    .data$gene
  ) |>
  unique()

selected_control_found <- selected_control_found[
  !is.na(selected_control_found)
]

selected_control_missing <- control_genes[
  !tolower(control_genes) %in%
    tolower(selected_control_found)
]

message(
  "Control genes found: ",
  if (length(selected_control_found) > 0L) {
    paste(selected_control_found, collapse = ", ")
  } else {
    "none"
  }
)

if (length(selected_control_missing) > 0L) {
  message(
    "Control genes NOT found in Crescendo universe: ",
    paste(selected_control_missing, collapse = ", ")
  )
}

selected_gene_table <- selected_feature_annotation |>
  dplyr::left_join(
    old_fdr005_genes |>
      dplyr::select(
        "ensembl_gene_id",
        old_minimum_FDR =
          "minimum_FDR",
        old_maximum_abs_log2FC =
          "maximum_abs_log2FC",
        "n_FDR005_clusters",
        "FDR005_clusters"
      ),
    by = "ensembl_gene_id"
  ) |>
  dplyr::arrange(
    dplyr::desc(.data$selected_control_gene),
    .data$old_minimum_FDR,
    .data$gene
  )

readr::write_tsv(
  selected_gene_table,
  selected_gene_list_tsv
)

# ==============================================================================
# 9. Annotate old DEG universe before metric calculation
# ==============================================================================

stage_message(
  4,
  "PREPARING OLD DEG SURVIVAL TABLES"
)

old_deg_genes_annotated <-
  annotate_deg_summary_with_crescendo_universe(
    deg_summary =
      old_deg_genes,
    feature_annotation =
      feature_annotation,
    after_deg_summary =
      after_deg_genes
  )

after_deg_genes_annotated <-
  annotate_deg_summary_with_crescendo_universe(
    deg_summary =
      after_deg_genes,
    feature_annotation =
      feature_annotation
  )

readr::write_tsv(
  old_sig_cluster_gene,
  old_deg_cluster_gene_tsv
)

readr::write_tsv(
  after_sig_cluster_gene,
  after_deg_cluster_gene_tsv
)

readr::write_tsv(
  old_deg_genes_annotated,
  old_deg_genes_tsv
)

readr::write_tsv(
  after_deg_genes_annotated,
  after_deg_genes_tsv
)

# Strict old significant cluster-gene rows compared with the exact same
# cluster-gene entry after Crescendo.
after_same_cluster_gene <- after_main_results |>
  dplyr::select(
    cluster_id =
      "cluster_id",
    ensembl_gene_id =
      "ensembl_gene_id",
    after_gene =
      "gene",
    after_logFC =
      "logFC",
    after_PValue =
      "PValue",
    after_FDR =
      "FDR"
  )

old_vs_after_cluster_gene <-
  old_sig_cluster_gene |>
    dplyr::transmute(
      cluster_id =
        .data$cluster_id,
      ensembl_gene_id =
        .data$ensembl_gene_id,
      gene =
        .data$gene,
      old_logFC =
        .data$logFC,
      old_PValue =
        .data$PValue,
      old_FDR =
        .data$FDR
    ) |>
    dplyr::left_join(
      after_same_cluster_gene,
      by = c(
        "cluster_id",
        "ensembl_gene_id"
      )
    ) |>
    dplyr::mutate(
      tested_after_crescendo =
        is.finite(.data$after_FDR) &
          is.finite(.data$after_logFC),

      significant_after_same_cluster =
        .data$tested_after_crescendo &
          .data$after_FDR <
            de_fdr_threshold &
          abs(.data$after_logFC) >=
            de_abs_log2fc_threshold,

      abs_log2FC_change_after_minus_before =
        abs(.data$after_logFC) -
          abs(.data$old_logFC),

      abs_log2FC_ratio_after_to_before =
        dplyr::if_else(
          is.finite(.data$after_logFC) &
            abs(.data$old_logFC) > 0,
          abs(.data$after_logFC) /
            abs(.data$old_logFC),
          NA_real_
        ),

      cluster_gene_status =
        dplyr::case_when(
          !.data$tested_after_crescendo ~
            "not_tested_after",
          .data$significant_after_same_cluster ~
            "retained_significant_same_cluster",
          TRUE ~
            "lost_significance_same_cluster"
        )
    )

readr::write_tsv(
  old_vs_after_cluster_gene,
  old_vs_after_cluster_gene_tsv
)

# ==============================================================================
# 10. Load Seurat raw counts and construct model metadata
# ==============================================================================

stage_message(
  5,
  "LOADING RAW SPOT COUNTS"
)

seurat_loaded <- find_single_seurat_object(
  rdata_file =
    input_seurat_rdata_file,
  required_cluster_column =
    cluster_column
)

seurat_object <- seurat_loaded$object

message(
  "Selected Seurat object: ",
  seurat_loaded$object_name
)

if (!sample_column %in%
    colnames(
      seurat_object@meta.data
    )) {
  stop(
    "Seurat metadata lacks sample_ID.",
    call. = FALSE
  )
}

sample_ids_use <- unique(
  as.character(
    sample_metadata$sample_ID
  )
)

keep_spots <- as.character(
  seurat_object@meta.data[[sample_column]]
) %in%
  sample_ids_use

seurat_object <- subset(
  seurat_object,
  cells =
    colnames(seurat_object)[
      keep_spots
    ]
)

if (
  length(
    unique(
      as.character(
        seurat_object@meta.data[[sample_column]]
      )
    )
  ) !=
    length(sample_ids_use)
) {
  stop(
    "Retained Seurat samples do not match Crescendo sample metadata.",
    call. = FALSE
  )
}

spot_metadata <- seurat_object@meta.data[
  ,
  c(
    sample_column,
    cluster_column
  ),
  drop = FALSE
]

spot_metadata$sample_ID <- as.character(
  spot_metadata[[sample_column]]
)

spot_metadata$biological_cluster <- as.character(
  spot_metadata[[cluster_column]]
)

sample_metadata_for_join <- sample_metadata |>
  dplyr::mutate(
    sample_ID =
      as.character(.data$sample_ID)
  ) |>
  dplyr::select(
    "sample_ID",
    "slide_ID",
    dplyr::any_of(
      c(
        "fmt_donor_group",
        "sex"
      )
    )
  ) |>
  dplyr::distinct(
    .data$sample_ID,
    .keep_all = TRUE
  )

spot_metadata <- tibble::rownames_to_column(
  as.data.frame(spot_metadata),
  var = "spot_ID"
) |>
  dplyr::left_join(
    sample_metadata_for_join,
    by = "sample_ID"
  )

if (anyNA(spot_metadata$slide_ID)) {
  stop(
    "Some spots are missing slide_ID.",
    call. = FALSE
  )
}

rownames(spot_metadata) <-
  spot_metadata$spot_ID

raw_counts_all <- get_raw_counts(
  seurat_object,
  assay_name
)

missing_spots <- setdiff(
  rownames(spot_metadata),
  colnames(raw_counts_all)
)

if (length(missing_spots) > 0L) {
  stop(
    "Raw count matrix is missing spots: ",
    length(missing_spots),
    call. = FALSE
  )
}

raw_counts_all <- raw_counts_all[
  ,
  rownames(spot_metadata),
  drop = FALSE
]

raw_nUMI <- Matrix::colSums(
  raw_counts_all
)

if (any(raw_nUMI <= 0)) {
  stop(
    "At least one spot has raw nUMI <= 0.",
    call. = FALSE
  )
}

missing_selected_genes <- setdiff(
  selected_crescendo_gene_ids,
  rownames(raw_counts_all)
)

if (length(missing_selected_genes) > 0L) {
  stop(
    "Raw counts missing ",
    length(missing_selected_genes),
    " selected Crescendo genes.",
    call. = FALSE
  )
}

raw_counts <- raw_counts_all[
  selected_crescendo_gene_ids,
  ,
  drop = FALSE
]

rm(
  raw_counts_all,
  seurat_object
)

invisible(
  gc()
)

model_metadata <- data.frame(
  biological_cluster =
    factor(
      spot_metadata$biological_cluster
    ),
  slide_ID =
    factor(
      spot_metadata$slide_ID
    ),
  offset_log_nUMI =
    log(
      as.numeric(
        raw_nUMI[
          rownames(spot_metadata)
        ]
      )
    ),
  row.names =
    rownames(spot_metadata),
  stringsAsFactors = FALSE
)

model_metadata$cluster_slide <-
  interaction(
    model_metadata$biological_cluster,
    model_metadata$slide_ID,
    drop = TRUE,
    sep = "__"
  )

message(
  "Raw matrix for SELECTED-GENE BVR/CVR: ",
  format(
    nrow(raw_counts),
    big.mark = ","
  ),
  " genes x ",
  format(
    ncol(raw_counts),
    big.mark = ","
  ),
  " spots"
)

message(
  "Biological clusters: ",
  nlevels(
    model_metadata$biological_cluster
  )
)

message(
  "Slides: ",
  nlevels(
    model_metadata$slide_ID
  )
)

# ==============================================================================
# 11. Find and validate corrected spot chunks
# ==============================================================================

stage_message(
  6,
  "FINDING CRESCENDO CORRECTED SPOT CHUNKS"
)

corrected_chunk_files <- list.files(
  corrected_spot_chunks_dir,
  pattern =
    "_crescendoCorrected_spotCounts\\.rds$",
  full.names = TRUE
)

if (length(corrected_chunk_files) == 0L) {
  stop(
    "No corrected spot chunk RDS files found.",
    call. = FALSE
  )
}

extract_chunk_number <- function(file_name) {
  as.integer(
    sub(
      ".*_chunk_([0-9]+)_genes_.*",
      "\\1",
      basename(file_name)
    )
  )
}

chunk_numbers <- vapply(
  corrected_chunk_files,
  extract_chunk_number,
  FUN.VALUE = integer(1)
)

if (anyNA(chunk_numbers)) {
  stop(
    "Could not parse chunk number from all files.",
    call. = FALSE
  )
}

chunk_order <- order(
  chunk_numbers
)

corrected_chunk_files <-
  corrected_chunk_files[
    chunk_order
  ]

chunk_numbers <-
  chunk_numbers[
    chunk_order
  ]

message(
  "Corrected chunks found: ",
  length(corrected_chunk_files)
)

# Validate all gene IDs without keeping all corrected matrices in RAM.
chunk_gene_ids <- character(0)

for (file_index in seq_along(corrected_chunk_files)) {
  chunk_object <- readRDS(
    corrected_chunk_files[[file_index]]
  )

  if (
    !"corrected_counts_gene_by_spot" %in%
      names(chunk_object)
  ) {
    stop(
      "Chunk missing corrected_counts_gene_by_spot:\n",
      corrected_chunk_files[[file_index]],
      call. = FALSE
    )
  }

  corrected_current <-
    chunk_object$corrected_counts_gene_by_spot

  chunk_gene_ids <- c(
    chunk_gene_ids,
    rownames(corrected_current)
  )

  rm(
    chunk_object,
    corrected_current
  )
}

if (anyDuplicated(chunk_gene_ids)) {
  stop(
    "Duplicated genes across corrected chunks.",
    call. = FALSE
  )
}

if (!setequal(
  chunk_gene_ids,
  crescendo_gene_ids
)) {
  stop(
    "Corrected chunk genes do not equal Crescendo gene_metadata universe.",
    call. = FALSE
  )
}

message(
  "Validated corrected gene universe: ",
  format(
    length(chunk_gene_ids),
    big.mark = ","
  )
)

# ==============================================================================
# 12. Calculate BVR/CVR one selected gene at a time
# ==============================================================================

stage_message(
  7,
  "CALCULATING BVR/CVR — ONE SELECTED GENE AT A TIME"
)

analysis_start_time <- Sys.time()

# Each completed gene receives its own result file. This makes progress fully
# visible and allows exact resume after interruption.
safe_gene_file_id <- function(feature_id) {
  gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    as.character(feature_id)
  )
}

gene_result_file <- function(feature_id) {
  file.path(
    output_chunks_dir,
    paste0(
      "BVR_CVR_gene_",
      safe_gene_file_id(feature_id),
      ".tsv.gz"
    )
  )
}

# Map every selected feature to the corrected Crescendo source chunk.
selected_gene_to_chunk <- vector(
  "list",
  length(selected_crescendo_gene_ids)
)

names(selected_gene_to_chunk) <-
  selected_crescendo_gene_ids

for (
  chunk_index in
    seq_along(corrected_chunk_files)
) {
  chunk_object <- readRDS(
    corrected_chunk_files[[chunk_index]]
  )

  corrected_ids <- rownames(
    chunk_object$corrected_counts_gene_by_spot
  )

  selected_here <- intersect(
    corrected_ids,
    selected_crescendo_gene_ids
  )

  if (length(selected_here) > 0L) {
    for (feature_id_current in selected_here) {
      selected_gene_to_chunk[[feature_id_current]] <- corrected_chunk_files[[chunk_index]]
    }
  }

  rm(chunk_object)
}

missing_chunk_mapping <- names(
  selected_gene_to_chunk
)[
  vapply(
    selected_gene_to_chunk,
    is.null,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_chunk_mapping) > 0L) {
  stop(
    "Could not locate corrected chunk for selected genes: ",
    paste(
      head(missing_chunk_mapping, 20),
      collapse = ", "
    ),
    call. = FALSE
  )
}

# Process genes in an informative order:
# controls first, then OLD FDR<0.05 genes ordered by old minimum FDR.
selected_gene_progress <- selected_gene_table |>
  dplyr::mutate(
    progress_order_control =
      dplyr::if_else(
        .data$selected_control_gene,
        0L,
        1L
      )
  ) |>
  dplyr::arrange(
    .data$progress_order_control,
    .data$old_minimum_FDR,
    .data$gene,
    .data$feature_id
  )

if (isTRUE(test_mode)) {
  selected_gene_progress <- head(
    selected_gene_progress,
    test_number_of_genes
  )
}

n_selected_to_run <- nrow(
  selected_gene_progress
)

message(
  "Selected genes to evaluate: ",
  n_selected_to_run
)

message(
  "Processing mode: ONE GENE AT A TIME"
)

message(
  "Each completed gene is saved immediately in:"
)

message(output_chunks_dir)

current_loaded_chunk_file <- NA_character_
current_corrected_chunk <- NULL

for (
  gene_counter in
    seq_len(n_selected_to_run)
) {
  gene_info <- selected_gene_progress[
    gene_counter,
    ,
    drop = FALSE
  ]

  feature_id_current <-
    gene_info$feature_id[[1L]]

  gene_symbol_current <-
    gene_info$gene[[1L]]

  if (
    is.na(gene_symbol_current) ||
      gene_symbol_current == ""
  ) {
    gene_symbol_current <-
      feature_id_current
  }

  selection_reason_current <-
    gene_info$selection_reason[[1L]]

  old_fdr_current <-
    gene_info$old_minimum_FDR[[1L]]

  result_file_current <-
    gene_result_file(
      feature_id_current
    )

  message("\n------------------------------------------------------------")
  message(
    "[",
    gene_counter,
    "/",
    n_selected_to_run,
    "] ",
    gene_symbol_current,
    " | ",
    feature_id_current
  )

  message(
    "Reason: ",
    selection_reason_current,
    " | OLD min FDR: ",
    if (is.finite(old_fdr_current)) {
      format(
        old_fdr_current,
        scientific = TRUE,
        digits = 3
      )
    } else {
      "control only"
    }
  )

  if (
    resume_existing_genes &&
      file.exists(result_file_current) &&
      file.info(result_file_current)$size > 0
  ) {
    message("Status: already finished -> SKIP")
    next
  }

  gene_start_time <- Sys.time()

  corrected_file_current <-
    selected_gene_to_chunk[[feature_id_current]]

  # Avoid repeatedly loading the same original Crescendo chunk when consecutive
  # selected genes happen to reside in it.
  if (
    is.na(current_loaded_chunk_file) ||
      !identical(
        current_loaded_chunk_file,
        corrected_file_current
      )
  ) {
    chunk_object <- readRDS(
      corrected_file_current
    )

    current_corrected_chunk <- methods::as(
      chunk_object$corrected_counts_gene_by_spot,
      "dgCMatrix"
    )

    current_loaded_chunk_file <-
      corrected_file_current

    rm(chunk_object)
  }

  if (!feature_id_current %in%
      rownames(current_corrected_chunk)) {
    stop(
      "Selected feature not present in loaded corrected chunk: ",
      feature_id_current,
      call. = FALSE
    )
  }

  missing_chunk_spots <- setdiff(
    rownames(model_metadata),
    colnames(current_corrected_chunk)
  )

  if (length(missing_chunk_spots) > 0L) {
    stop(
      "Corrected chunk missing ",
      length(missing_chunk_spots),
      " spots.",
      call. = FALSE
    )
  }

  corrected_gene <- current_corrected_chunk[
    feature_id_current,
    rownames(model_metadata),
    drop = FALSE
  ]

  raw_gene <- raw_counts[
    feature_id_current,
    rownames(model_metadata),
    drop = FALSE
  ]

  result_current <- tryCatch(
    fit_one_gene_bvr_cvr(
      gene_index = 1L,
      raw_chunk = raw_gene,
      corrected_chunk = corrected_gene,
      model_metadata = model_metadata,
      verbose = TRUE
    ),
    error = function(e) {
      tibble::tibble(
        feature_id =
          feature_id_current,
        total_raw_counts =
          NA_real_,
        total_corrected_counts =
          NA_real_,
        detected_spots_before =
          NA_integer_,
        detected_spots_after =
          NA_integer_,
        fit_ok_before =
          FALSE,
        fit_ok_after =
          FALSE,
        singular_before =
          NA,
        singular_after =
          NA,
        optimizer_before =
          NA_character_,
        optimizer_after =
          NA_character_,
        biological_cluster_variance_before =
          NA_real_,
        biological_cluster_variance_after =
          NA_real_,
        slide_variance_before =
          NA_real_,
        slide_variance_after =
          NA_real_,
        cluster_slide_variance_before =
          NA_real_,
        cluster_slide_variance_after =
          NA_real_,
        batch_variance_before =
          NA_real_,
        batch_variance_after =
          NA_real_,
        BVR =
          NA_real_,
        CVR =
          NA_real_,
        metric_valid =
          FALSE,
        pre_batch_variance_ge_0p001 =
          FALSE,
        BVR_lt_1 =
          FALSE,
        CVR_ge_0p5 =
          FALSE,
        CVR_ge_1 =
          FALSE,
        target_zone =
          FALSE,
        warning_before =
          "",
        warning_after =
          "",
        error_before =
          conditionMessage(e),
        error_after =
          conditionMessage(e)
      )
    }
  )

  result_current <- result_current |>
    dplyr::mutate(
      gene =
        gene_symbol_current,
      selection_reason =
        selection_reason_current,
      old_minimum_FDR =
        old_fdr_current,
      .after =
        "feature_id"
    )

  readr::write_tsv(
    result_current,
    result_file_current
  )

  gene_elapsed_minutes <- as.numeric(
    difftime(
      Sys.time(),
      gene_start_time,
      units = "mins"
    )
  )

  completed_result_files <- vapply(
    selected_gene_progress$feature_id,
    function(feature_id_check) {
      file_check <-
        gene_result_file(
          feature_id_check
        )

      file.exists(file_check) &&
        file.info(file_check)$size > 0
    },
    FUN.VALUE = logical(1)
  )

  n_completed_now <- sum(
    completed_result_files
  )

  total_elapsed_minutes <- as.numeric(
    difftime(
      Sys.time(),
      analysis_start_time,
      units = "mins"
    )
  )

  mean_minutes_per_completed_gene <- if (
    n_completed_now > 0L
  ) {
    total_elapsed_minutes /
      n_completed_now
  } else {
    NA_real_
  }

  remaining_genes <-
    n_selected_to_run -
      n_completed_now

  estimated_remaining_hours <- if (
    is.finite(
      mean_minutes_per_completed_gene
    )
  ) {
    mean_minutes_per_completed_gene *
      remaining_genes /
      60
  } else {
    NA_real_
  }

  if (
    isTRUE(
      result_current$metric_valid[[1L]]
    )
  ) {
    message(
      "Result: BVR = ",
      formatC(
        result_current$BVR[[1L]],
        format = "f",
        digits = 3
      ),
      " | CVR = ",
      formatC(
        result_current$CVR[[1L]],
        format = "f",
        digits = 3
      ),
      " | target = ",
      if (
        isTRUE(
          result_current$target_zone[[1L]]
        )
      ) {
        "YES"
      } else {
        "NO"
      }
    )
  } else {
    message(
      "Result: metric invalid; inspect warning/error columns."
    )
  }

  message(
    "Saved: ",
    basename(result_file_current)
  )

  message(
    "Gene time: ",
    sprintf(
      "%.2f",
      gene_elapsed_minutes
    ),
    " min | completed: ",
    n_completed_now,
    "/",
    n_selected_to_run,
    " | estimated remaining: ",
    if (is.finite(estimated_remaining_hours)) {
      paste0(
        sprintf(
          "%.2f",
          estimated_remaining_hours
        ),
        " h"
      )
    } else {
      "NA"
    }
  )
}

# Collect exact per-gene result files corresponding to this selected run.
metric_chunk_files_existing <- vapply(
  selected_gene_progress$feature_id,
  gene_result_file,
  FUN.VALUE = character(1)
)

metric_chunk_files_existing <-
  metric_chunk_files_existing[
    file.exists(
      metric_chunk_files_existing
    )
  ]

# ==============================================================================
# 13. Combine metric chunks
# ==============================================================================

stage_message(
  8,
  "COMBINING BVR/CVR RESULTS"
)

if (length(metric_chunk_files_existing) == 0L) {
  stop(
    "No metric chunk files were generated.",
    call. = FALSE
  )
}

all_metrics <- dplyr::bind_rows(
  lapply(
    metric_chunk_files_existing,
    function(file_current) {
      readr::read_tsv(
        file_current,
        show_col_types = FALSE,
        progress = FALSE
      )
    }
  )
)

if (!isTRUE(test_mode)) {
  if (nrow(all_metrics) !=
      nrow(selected_gene_progress)) {
    stop(
      "Expected ",
      nrow(selected_gene_progress),
      " selected-gene BVR/CVR rows; found ",
      nrow(all_metrics),
      ".",
      call. = FALSE
    )
  }

  if (anyDuplicated(
    all_metrics$feature_id
  )) {
    stop(
      "Duplicated feature IDs in combined BVR/CVR results.",
      call. = FALSE
    )
  }
}

all_metrics_annotated <-
  all_metrics |>
    dplyr::left_join(
      feature_annotation |>
        dplyr::select(
          "feature_id",
          "ensembl_gene_id",
          "gene",
          "edgeR_mapping_type",
          "old_DEG",
          "after_DEG",
          "DEG_category",
          "selected_old_FDR005",
          "selected_control_gene",
          "selected_for_BVR_CVR"
        ),
      by = "feature_id"
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$old_DEG),
      .data$BVR,
      dplyr::desc(.data$CVR)
    )

readr::write_tsv(
  all_metrics_annotated,
  all_gene_metrics_tsv
)

# ==============================================================================
# 14. Join BVR/CVR to the old ~1200-gene list
# ==============================================================================

stage_message(
  9,
  "ANNOTATING OLD DEGs WITH BVR/CVR"
)

metrics_gene_lookup <-
  all_metrics_annotated |>
    dplyr::filter(
      !is.na(.data$ensembl_gene_id)
    ) |>
    dplyr::select(
      "ensembl_gene_id",
      "feature_id",
      "BVR",
      "CVR",
      "batch_variance_before",
      "batch_variance_after",
      "biological_cluster_variance_before",
      "biological_cluster_variance_after",
      "pre_batch_variance_ge_0p001",
      "target_zone",
      "metric_valid"
    ) |>
    dplyr::distinct(
      .data$ensembl_gene_id,
      .keep_all = TRUE
    )

old_gene_survival <-
  old_deg_genes_annotated |>
    dplyr::left_join(
      metrics_gene_lookup,
      by = "ensembl_gene_id"
    ) |>
    dplyr::arrange(
      .data$DEG_survival_status,
      .data$BVR,
      .data$minimum_FDR
    )

readr::write_tsv(
  old_gene_survival,
  old_gene_survival_tsv
)

# ==============================================================================
# 15. DEG count summary
# ==============================================================================

stage_message(
  10,
  "SUMMARIZING DEG LOSS / RETENTION"
)

n_old_unique <- nrow(
  old_deg_genes_annotated
)

n_old_in_crescendo <- sum(
  old_deg_genes_annotated$
    in_crescendo_corrected_gene_universe,
  na.rm = TRUE
)

n_old_not_in_crescendo <-
  n_old_unique -
    n_old_in_crescendo

n_after_unique <- nrow(
  after_deg_genes_annotated
)

n_old_retained <- sum(
  old_deg_genes_annotated$
    significant_after_crescendo &
    old_deg_genes_annotated$
      in_crescendo_corrected_gene_universe,
  na.rm = TRUE
)

n_old_lost <- sum(
  !old_deg_genes_annotated$
    significant_after_crescendo &
    old_deg_genes_annotated$
      in_crescendo_corrected_gene_universe,
  na.rm = TRUE
)

n_after_new <- sum(
  !after_deg_genes_annotated$
    ensembl_gene_id %in%
      old_deg_genes_annotated$
        ensembl_gene_id,
  na.rm = TRUE
)

deg_count_summary <- tibble::tibble(
  metric = c(
    "old_unique_DEGs",
    "old_DEGs_in_Crescendo_gene_universe",
    "old_DEGs_not_in_Crescendo_gene_universe",
    "after_Crescendo_unique_DEGs",
    "old_DEGs_retained_after_Crescendo",
    "old_DEGs_lost_after_Crescendo_among_corrected",
    "new_after_Crescendo_DEGs",
    "percent_old_corrected_DEGs_lost",
    "percent_old_corrected_DEGs_retained"
  ),
  value = c(
    n_old_unique,
    n_old_in_crescendo,
    n_old_not_in_crescendo,
    n_after_unique,
    n_old_retained,
    n_old_lost,
    n_after_new,
    if (n_old_in_crescendo > 0L) {
      100 * n_old_lost /
        n_old_in_crescendo
    } else {
      NA_real_
    },
    if (n_old_in_crescendo > 0L) {
      100 * n_old_retained /
        n_old_in_crescendo
    } else {
      NA_real_
    }
  )
)

readr::write_tsv(
  deg_count_summary,
  deg_count_summary_tsv
)

# ==============================================================================
# 16. BVR/CVR summaries
# ==============================================================================

stage_message(
  11,
  "SUMMARIZING BVR/CVR"
)

metric_sets <- list(
  all_selected_genes =
    all_metrics_annotated,

  old_FDR005_selected =
    all_metrics_annotated |>
      dplyr::filter(
        .data$selected_old_FDR005
      ),

  control_genes =
    all_metrics_annotated |>
      dplyr::filter(
        .data$selected_control_gene
      ),

  old_DEG =
    all_metrics_annotated |>
      dplyr::filter(
        .data$old_DEG
      ),

  old_DEG_lost_after =
    all_metrics_annotated |>
      dplyr::filter(
        .data$DEG_category ==
          "old_lost_after"
      ),

  old_DEG_retained_after =
    all_metrics_annotated |>
      dplyr::filter(
        .data$DEG_category ==
          "old_retained_after"
      ),

  after_DEG =
    all_metrics_annotated |>
      dplyr::filter(
        .data$after_DEG
      ),

  never_DEG =
    all_metrics_annotated |>
      dplyr::filter(
        .data$DEG_category ==
          "never_DEG"
      )
)

summary_rows <- list()
summary_counter <- 0L

for (set_name in names(metric_sets)) {
  for (scope_current in c(
    "all_valid",
    "meaningful_pre_batch_variance"
  )) {
    summary_counter <-
      summary_counter + 1L

    summary_rows[[summary_counter]] <- summarize_metric_subset(
      data =
        metric_sets[[set_name]],
      set_name =
        set_name,
      evaluation_scope =
        scope_current
    )
  }
}

metric_summary <- dplyr::bind_rows(
  summary_rows
)

readr::write_tsv(
  metric_summary,
  global_summary_tsv
)

readr::write_tsv(
  metric_summary |>
    dplyr::filter(
      .data$gene_set %in%
        c(
          "old_DEG",
          "old_DEG_lost_after",
          "old_DEG_retained_after",
          "after_DEG",
          "never_DEG"
        )
    ),
  deg_category_bvr_summary_tsv
)

model_status_summary <- tibble::tibble(
  n_genes =
    nrow(all_metrics_annotated),
  n_fit_ok_before =
    sum(
      all_metrics_annotated$fit_ok_before,
      na.rm = TRUE
    ),
  n_fit_ok_after =
    sum(
      all_metrics_annotated$fit_ok_after,
      na.rm = TRUE
    ),
  n_valid_BVR_CVR =
    sum(
      all_metrics_annotated$metric_valid,
      na.rm = TRUE
    ),
  n_singular_before =
    sum(
      all_metrics_annotated$singular_before,
      na.rm = TRUE
    ),
  n_singular_after =
    sum(
      all_metrics_annotated$singular_after,
      na.rm = TRUE
    ),
  n_pre_batch_variance_ge_0p001 =
    sum(
      all_metrics_annotated$
        pre_batch_variance_ge_0p001,
      na.rm = TRUE
    )
)

readr::write_tsv(
  model_status_summary,
  model_status_tsv
)

# ==============================================================================
# 17. BVR/CVR plots
# ==============================================================================

stage_message(
  12,
  "CREATING BVR/CVR PLOTS"
)

plot_data_all <-
  all_metrics_annotated |>
    dplyr::filter(
      .data$metric_valid,
      .data$pre_batch_variance_ge_0p001,
      is.finite(.data$BVR),
      is.finite(.data$CVR),
      .data$BVR > 0,
      .data$CVR > 0
    )

if (nrow(plot_data_all) > 0L) {
  p_all <- ggplot2::ggplot(
    plot_data_all,
    ggplot2::aes(
      x = .data$CVR,
      y = .data$BVR
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        shape = .data$old_DEG
      ),
      alpha = 0.45,
      size = 1.4
    ) +
    ggplot2::geom_vline(
      xintercept = 0.5,
      linetype = 2,
      linewidth = 0.6
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = 2,
      linewidth = 0.6
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title =
        "Crescendo BVR/CVR — selected OLD FDR<0.05 genes + controls",
      subtitle =
        paste0(
          "Leiden C1-C16 used as biological-cluster proxy | ",
          "slide_ID = batch | pre batch variance >= ",
          minimum_pre_batch_variance_for_evaluation,
          "\nTarget zone: BVR < 1 and CVR >= 0.5"
        ),
      x =
        "CVR: biological-cluster variance AFTER / BEFORE",
      y =
        "BVR: batch variance AFTER / BEFORE",
      shape =
        "Old DEG"
    ) +
    ggplot2::theme_bw(
      base_size = 11
    )

  ggplot2::ggsave(
    filename =
      all_genes_pdf,
    plot =
      p_all,
    width =
      9,
    height =
      7,
    units =
      "in",
    limitsize =
      FALSE
  )
}

plot_data_old <-
  plot_data_all |>
    dplyr::filter(
      .data$old_DEG
    )

if (nrow(plot_data_old) > 0L) {
  p_old <- ggplot2::ggplot(
    plot_data_old,
    ggplot2::aes(
      x = .data$CVR,
      y = .data$BVR
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        shape = .data$after_DEG
      ),
      alpha = 0.65,
      size = 1.8
    ) +
    ggplot2::geom_vline(
      xintercept = 0.5,
      linetype = 2,
      linewidth = 0.6
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = 2,
      linewidth = 0.6
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title =
        "Crescendo BVR/CVR — genes significant before correction",
      subtitle =
        paste0(
          "Old selection: FDR < ",
          de_fdr_threshold,
          " & |log2FC| >= ",
          de_abs_log2fc_threshold,
          " | pre batch variance >= ",
          minimum_pre_batch_variance_for_evaluation,
          "\nTarget zone: BVR < 1 and CVR >= 0.5"
        ),
      x =
        "CVR: biological-cluster variance AFTER / BEFORE",
      y =
        "BVR: batch variance AFTER / BEFORE",
      shape =
        "Still DEG after Crescendo"
    ) +
    ggplot2::theme_bw(
      base_size = 11
    )

  ggplot2::ggsave(
    filename =
      old_deg_pdf,
    plot =
      p_old,
    width =
      9,
    height =
      7,
    units =
      "in",
    limitsize =
      FALSE
  )
}

# ==============================================================================
# 18. Console interpretation table
# ==============================================================================

stage_message(
  13,
  "FINAL SUMMARY"
)

message("\nDEG COUNTS:")
print(
  deg_count_summary,
  n = Inf
)

message(
  "\nBVR/CVR SUMMARY — meaningful pre-batch variance:"
)

print(
  metric_summary |>
    dplyr::filter(
      .data$evaluation_scope ==
        "meaningful_pre_batch_variance"
    ) |>
    dplyr::select(
      "gene_set",
      "n_genes_total",
      "n_genes_valid_metrics",
      "median_BVR",
      "median_CVR",
      "percent_BVR_lt_1",
      "percent_CVR_ge_0p5",
      "percent_target_zone"
    ),
  n = Inf
)

message("\n============================================================")
message("HOW TO READ THE RESULT")
message("============================================================")
message("BVR < 1        = slide-associated variance decreased")
message("CVR >= 0.5     = biological-cluster variance reasonably preserved")
message("Target zone    = BVR < 1 AND CVR >= 0.5")
message("")
message(
  "Most informative comparison for the old DE signal:"
)
message(
  "  old_DEG_lost_after vs old_DEG_retained_after"
)
message("")
message(
  "If old_DEG_lost_after has clearly lower BVR and good CVR, "
)
message(
  "the disappearance of many old DEGs is consistent with removal of "
)
message(
  "slide-associated technical variation rather than simple destruction "
)
message(
  "of biological cluster structure."
)
message("")
message(
  "Genes not present in the Crescendo-corrected universe are reported "
)
message(
  "separately and must NOT be counted as genes removed by correction."
)

message("\n============================================================")
message("BVR/CVR ANALYSIS FINISHED")
message("============================================================")

message("\nSelected gene list:")
message(
  normalizePath(
    selected_gene_list_tsv,
    mustWork = TRUE
  )
)

message("\nSelected-gene BVR/CVR:")
message(
  normalizePath(
    all_gene_metrics_tsv,
    mustWork = TRUE
  )
)

message("\nBVR/CVR summary:")
message(
  normalizePath(
    global_summary_tsv,
    mustWork = TRUE
  )
)

message("\nOld DEG list:")
message(
  normalizePath(
    old_deg_genes_tsv,
    mustWork = TRUE
  )
)

message("\nOld cluster-gene BEFORE vs AFTER:")
message(
  normalizePath(
    old_vs_after_cluster_gene_tsv,
    mustWork = TRUE
  )
)

message("\nOld DEG survival + BVR/CVR:")
message(
  normalizePath(
    old_gene_survival_tsv,
    mustWork = TRUE
  )
)

message("\nDEG count summary:")
message(
  normalizePath(
    deg_count_summary_tsv,
    mustWork = TRUE
  )
)

# ==============================================================================
# End
# ==============================================================================
