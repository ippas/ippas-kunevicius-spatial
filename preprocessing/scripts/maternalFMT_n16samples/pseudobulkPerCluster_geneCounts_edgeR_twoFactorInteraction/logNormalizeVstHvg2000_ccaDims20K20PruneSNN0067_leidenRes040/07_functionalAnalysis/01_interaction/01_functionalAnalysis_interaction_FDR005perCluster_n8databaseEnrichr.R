#!/usr/bin/env Rscript

# ==============================================================================
# 01_functionalAnalysis_interaction_multiEnrichrLibraries_v01_fullRun.R
#
# Functional enrichment analysis for the Group x Sex interaction from the
# cluster-specific pseudobulk edgeR two-factor model.
#
# ENRICHR LIBRARIES
#   ChEA_2022
#   Reactome_Pathways_2024
#   WikiPathways_2024_Mouse
#   GO_Biological_Process_2026
#   GO_Cellular_Component_2026
#   GO_Molecular_Function_2026
#   MGI_Mammalian_Phenotype_Level_4_2024
#   DisGeNET
#
# INTERACTION GENE SELECTION
#   FDR < 0.05
#   No |logFC| threshold.
#
# INTERACTION COEFFICIENT VARIANTS
#   ALL
#   PositiveLogFC
#   NegativeLogFC
#
# IMPORTANT
#   PositiveLogFC / NegativeLogFC describe only the sign of the interaction
#   coefficient. They must NOT be interpreted or labelled as UP / DOWN.
#
# ENRICHMENT
#   Standard Enrichr through enrichR.
#   No custom background / Speedrichr.
#   No explicit ortholog conversion.
#   Unique edgeR gene symbols are submitted directly to Enrichr.
#
# PER-DATABASE WORK
#   3 interaction-coefficient variants x 17 scopes = 51 Enrichr lists.
#   One progress bar is shown independently for each database.
#
# PER-DATABASE OUTPUT
#   <functional-analysis-dir>/<library_id>/xlsx/ : exactly 3 XLSX files
#
# XLSX
#   Each workbook contains 17 sheets:
#     genes_from_allClusters
#     C1-C16 with short anatomical names
#
#   Visual highlight only:
#     Enrichr FDR < 0.10 and overlap >= 2 genes
#     fill = #DEDCE6
#
# LONG-RUN SAFETY
#   Completed databases can be skipped automatically.
#   If one database fails, later databases can continue.
#   Failures are summarized at the end.
#
# DEVELOPMENT MODE
#   test_first_n_lists <- NULL  -> full run
#   test_first_n_lists <- 3L    -> first 3 lists per database, etc.
# ==============================================================================

options(stringsAsFactors = FALSE)

# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tibble",
  "openxlsx",
  "enrichR"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them, for example with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(openxlsx)
  library(enrichR)
})

# ==============================================================================
# 2. Analysis configuration
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"
clustering_name <- "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

# Exact interaction selection used by the interaction heatmap.
interaction_fdr_threshold <- 0.05

# Enrichr result criteria.
enrichment_fdr_threshold <- 0.05
min_overlap_genes <- 3L

# XLSX visual highlighting only.
xlsx_highlight_fdr_threshold <- 0.10
xlsx_highlight_min_overlap_genes <- 2L

enrichr_libraries <- data.frame(
  library_id = c(
    "ChEA_2022",
    "Reactome_Pathways_2024",
    "WikiPathways_2024_Mouse",
    "GO_Biological_Process_2026",
    "GO_Cellular_Component_2026",
    "GO_Molecular_Function_2026",
    "MGI_Mammalian_Phenotype_Level_4_2024",
    "DisGeNET"
  ),
  display_name = c(
    "ChEA 2022",
    "Reactome Pathways 2024",
    "WikiPathways 2024 Mouse",
    "GO Biological Process 2026",
    "GO Cellular Component 2026",
    "GO Molecular Function 2026",
    "MGI Mammalian Phenotype Level 4 2024",
    "DisGeNET"
  ),
  stringsAsFactors = FALSE
)

enrichr_libraries <- data.frame(
  library_id = c(
    "GWAS_Catalog_2025",
    "GeDiPNet_2023"
  ),
  display_name = c(
    "GWAS Catalog 2025",
    "GeDiPNet 2023"
  ),
  stringsAsFactors = FALSE
)

# NULL = all 51 lists per database.
# Integer = first N lists per database.
test_first_n_lists <- NULL

# Skip a database when all expected 3 XLSX files already exist and are non-empty.
skip_completed_databases <- TRUE

# Continue with later databases if one database fails.
continue_after_database_error <- TRUE

# Enrichr retry settings.
enrichr_sleep_seconds <- 1
enrichr_max_attempts <- 3L
enrichr_retry_wait_seconds <- 3

parameter_block <- "FDR005perCluster"

interaction_subsets <- c(
  "ALL",
  "PositiveLogFC",
  "NegativeLogFC"
)

analysis_variants <- data.frame(
  xlsx_index = 1:3,
  subset_id = interaction_subsets,
  subset_label = c(
    "all significant interaction genes",
    "positive interaction logFC",
    "negative interaction logFC"
  ),
  file_suffix = c(
    "ALL_enrichrResults.xlsx",
    "PositiveLogFC_enrichrResults.xlsx",
    "NegativeLogFC_enrichrResults.xlsx"
  ),
  stringsAsFactors = FALSE
)

planned_total_lists <- nrow(analysis_variants) * 17L

if (!is.null(test_first_n_lists)) {
  if (
    length(test_first_n_lists) != 1L ||
      is.na(test_first_n_lists) ||
      !is.finite(test_first_n_lists) ||
      test_first_n_lists < 1 ||
      test_first_n_lists != as.integer(test_first_n_lists)
  ) {
    stop(
      "`test_first_n_lists` must be NULL or one positive integer.",
      call. = FALSE
    )
  }

  test_first_n_lists <- as.integer(test_first_n_lists)
}

test_mode <- !is.null(test_first_n_lists)

progress_total_lists <- if (test_mode) {
  min(test_first_n_lists, planned_total_lists)
} else {
  planned_total_lists
}

# ==============================================================================
# 3. Cluster anatomical labels
# ==============================================================================

custom_cluster_labels <- c(
  "1"  = "posterior & sensory relay thalamic nuclei",
  "2"  = "isocortex, layers 4 & 5",
  "3"  = "cortical layers 1 & hippocampal neuropil",
  "4"  = "Fiber tracts",
  "5"  = "cortical subplate & deep olfactory areas",
  "6"  = "hypothalamus",
  "7"  = "isocortex, layer 6",
  "8"  = "isocortex, layer 2/3",
  "9"  = "reticular, ventral geniculate & habenular region",
  "10" = "striatum-like amygdala nuclei",
  "11" = "medial thalamic nuclei",
  "12" = "caudoputamen",
  "13" = "hippocampal CA fields, pyramidal layer",
  "14" = "meninges & vasculature",
  "15" = "ventricles",
  "16" = "dentate gyrus"
)

expected_cluster_ids <- names(custom_cluster_labels)

# Short labels are used only for XLSX worksheet names.
# Excel/openxlsx worksheet names are limited to 31 characters.
xlsx_cluster_sheet_labels <- c(
  "1"  = "C1_posteriorSensoryThalamus",
  "2"  = "C2_isocortex_L4_L5",
  "3"  = "C3_corticalL1_hippNeuropil",
  "4"  = "C4_fiberTracts",
  "5"  = "C5_cortSubplate_deepOlfactory",
  "6"  = "C6_hypothalamus",
  "7"  = "C7_isocortex_L6",
  "8"  = "C8_isocortex_L2_3",
  "9"  = "C9_reticGeniculateHabenula",
  "10" = "C10_striatumLikeAmygdala",
  "11" = "C11_medialThalamus",
  "12" = "C12_caudoputamen",
  "13" = "C13_hippCA_pyramidal",
  "14" = "C14_meninges_vasculature",
  "15" = "C15_ventricles",
  "16" = "C16_dentateGyrus"
)

if (!identical(names(xlsx_cluster_sheet_labels), expected_cluster_ids)) {
  stop(
    "XLSX cluster sheet labels are not defined for exactly C1-C16 in order.",
    call. = FALSE
  )
}

if (any(nchar(xlsx_cluster_sheet_labels) > 31L)) {
  stop(
    "At least one XLSX worksheet name exceeds the Excel 31-character limit: ",
    paste(
      xlsx_cluster_sheet_labels[nchar(xlsx_cluster_sheet_labels) > 31L],
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (anyDuplicated(tolower(xlsx_cluster_sheet_labels)) > 0L) {
  stop(
    "Duplicated XLSX worksheet names detected (case-insensitive).",
    call. = FALSE
  )
}

# ==============================================================================
# 4. Input and output paths
# ==============================================================================

statistics_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)

analysis_prefix <- paste0(
  dataset_name,
  "_",
  clustering_name,
  "_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction"
)

edgeR_results_rdata_file <- file.path(
  statistics_dir,
  "03_edgeRResults",
  paste0(analysis_prefix, "_edgeRResults.RData")
)

functional_analysis_root <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "07_functionalAnalysis",
  "01_interaction"
)

dir.create(
  functional_analysis_root,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(edgeR_results_rdata_file)) {
  stop(
    "Missing edgeR results RData file:\n",
    edgeR_results_rdata_file,
    call. = FALSE
  )
}

build_database_output_paths <- function(database) {
  library_output_dir <- file.path(
    functional_analysis_root,
    database
  )

  if (test_mode) {
    library_output_dir <- file.path(
      library_output_dir,
      paste0("TEST_first", progress_total_lists)
    )
  }

  xlsx_output_dir <- file.path(
    library_output_dir,
    "xlsx"
  )

  dir.create(
    xlsx_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  output_prefix <- paste0(
    dataset_name,
    "_interaction_",
    parameter_block,
    "_",
    database
  )

  xlsx_files <- setNames(
    vapply(
      seq_len(nrow(analysis_variants)),
      function(i) {
        file.path(
          xlsx_output_dir,
          paste0(
            sprintf("%02d", analysis_variants$xlsx_index[[i]]),
            "_",
            output_prefix,
            "_",
            analysis_variants$file_suffix[[i]]
          )
        )
      },
      FUN.VALUE = character(1)
    ),
    analysis_variants$subset_id
  )

  list(
    library_output_dir = library_output_dir,
    xlsx_output_dir = xlsx_output_dir,
    output_prefix = output_prefix,
    xlsx_files = xlsx_files
  )
}

database_outputs_complete <- function(paths) {
  if (test_mode) {
    return(FALSE)
  }

  expected_files <- unname(paths$xlsx_files)

  if (!all(file.exists(expected_files))) {
    return(FALSE)
  }

  file_sizes <- file.info(expected_files)$size

  all(
    !is.na(file_sizes) &
      file_sizes > 0
  )
}

# ==============================================================================
# 5. Helper: identify the Group x Sex interaction test
# ==============================================================================

find_interaction_test_id <- function(test_definitions) {
  required_columns <- c(
    "test_id",
    "sheet_name",
    "comparison"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(test_definitions)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Test-definition table is missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  candidate_rows <- test_definitions |>
    dplyr::mutate(
      test_id_chr = as.character(.data$test_id),
      sheet_name_chr = as.character(.data$sheet_name),
      comparison_chr = as.character(.data$comparison)
    ) |>
    dplyr::filter(
      grepl("interaction", .data$test_id_chr, ignore.case = TRUE) |
        grepl("interaction", .data$sheet_name_chr, ignore.case = TRUE) |
        grepl("interaction", .data$comparison_chr, ignore.case = TRUE)
    )

  if (nrow(candidate_rows) != 1L) {
    message("\nInteraction-test candidates:")
    print(
      candidate_rows |>
        dplyr::select(
          .data$test_id,
          .data$sheet_name,
          .data$comparison
        ),
      n = Inf
    )

    stop(
      "Could not uniquely identify the interaction test. Candidates found: ",
      nrow(candidate_rows),
      call. = FALSE
    )
  }

  as.character(candidate_rows$test_id[[1]])
}

# ==============================================================================
# 6. Helpers: gene cleaning and empty-result tables
# ==============================================================================

clean_gene_vector <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x)]
  sort(unique(x))
}

empty_enrichment_table <- function(
  analysis_status,
  interaction_subset,
  cluster_id,
  cluster_name,
  n_query_genes
) {
  tibble::tibble(
    row_id = NA_integer_,
    term = NA_character_,
    p_value = NA_real_,
    FDR = NA_real_,
    odds_ratio = NA_real_,
    combined_score = NA_real_,
    overlap = NA_character_,
    overlap_genes = NA_character_,
    n_overlap_genes = NA_integer_,
    term_gene_set_size = NA_integer_,
    n_query_genes = as.integer(n_query_genes),
    significant = FALSE,
    interaction_subset = interaction_subset,
    cluster_id = cluster_id,
    cluster_name = cluster_name,
    analysis_status = analysis_status
  )
}

format_odds_ratio_for_excel <- function(x) {
  suppressWarnings(as.numeric(x))
}

# ==============================================================================
# 7. Helper: one database-wide progress bar
# ==============================================================================

format_elapsed_time <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) {
    return("--:--:--")
  }

  seconds <- as.integer(round(seconds))
  hours <- seconds %/% 3600L
  minutes <- (seconds %% 3600L) %/% 60L
  secs <- seconds %% 60L

  sprintf(
    "%02d:%02d:%02d",
    hours,
    minutes,
    secs
  )
}

create_database_progress <- function(database, total_lists, width = 32L) {
  progress <- new.env(parent = emptyenv())

  progress$database <- database
  progress$total <- as.integer(total_lists)
  progress$completed <- 0L
  progress$width <- as.integer(width)
  progress$start_time <- Sys.time()
  progress$last_label <- "starting"

  progress
}

draw_database_progress <- function(
  progress,
  current_label = NULL,
  finish = FALSE
) {
  if (!is.null(current_label)) {
    progress$last_label <- current_label
  }

  completed <- progress$completed
  total <- progress$total

  fraction <- if (total > 0L) {
    completed / total
  } else {
    1
  }

  fraction <- min(
    max(fraction, 0),
    1
  )

  filled <- floor(
    progress$width * fraction
  )

  bar <- paste0(
    "[",
    paste(rep("=", filled), collapse = ""),
    paste(rep(" ", progress$width - filled), collapse = ""),
    "]"
  )

  elapsed <- as.numeric(
    difftime(
      Sys.time(),
      progress$start_time,
      units = "secs"
    )
  )

  mean_per_list <- if (completed > 0L) {
    elapsed / completed
  } else {
    NA_real_
  }

  eta <- if (
    completed > 0L &&
      completed < total
  ) {
    mean_per_list * (total - completed)
  } else if (completed >= total) {
    0
  } else {
    NA_real_
  }

  avg_label <- if (is.finite(mean_per_list)) {
    sprintf(
      "%.1fs/list",
      mean_per_list
    )
  } else {
    "--.-s/list"
  }

  line <- paste0(
    "\r",
    progress$database,
    " ",
    bar,
    " ",
    completed,
    "/",
    total,
    " (",
    sprintf("%5.1f", 100 * fraction),
    "%)",
    " | elapsed ",
    format_elapsed_time(elapsed),
    " | avg ",
    avg_label,
    " | ETA ",
    format_elapsed_time(eta),
    " | ",
    progress$last_label
  )

  # Padding ensures that a shorter new label overwrites the previous label.
  cat(
    sprintf(
      "%-220s",
      line
    )
  )

  flush.console()

  if (finish) {
    cat("\n")
  }

  invisible(progress)
}

advance_database_progress <- function(
  progress,
  current_label
) {
  progress$completed <- progress$completed + 1L

  draw_database_progress(
    progress = progress,
    current_label = current_label,
    finish = progress$completed >= progress$total
  )
}

# ==============================================================================
# 8. Helper: standard Enrichr request with retries
# ==============================================================================

run_enrichr_with_retry <- function(
  query_genes,
  database
) {
  last_error <- NULL

  for (attempt in seq_len(enrichr_max_attempts)) {
    result <- NULL

    # enrichR prints status lines for each request. Capture them so the terminal
    # shows one persistent database-wide progress bar instead.
    invisible(
      capture.output(
        result <- tryCatch(
          suppressMessages(
            suppressWarnings(
              enrichR::enrichr(
                genes = query_genes,
                databases = database,
                sleepTime = enrichr_sleep_seconds
              )
            )
          ),
          error = function(e) {
            last_error <<- e
            NULL
          }
        ),
        type = "output"
      )
    )

    if (
      !is.null(result) &&
        database %in% names(result)
    ) {
      database_result <- result[[database]]

      if (
        is.null(database_result) ||
          nrow(database_result) == 0L
      ) {
        return(data.frame())
      }

      return(database_result)
    }

    # Some enrichR versions can throw a column-name assignment error when the
    # server returns zero rows. Treat that parser failure as a valid empty result.
    if (
      !is.null(last_error) &&
        grepl(
          "names.*attribute.*same length.*vector.*\\[0\\]",
          conditionMessage(last_error),
          ignore.case = TRUE
        )
    ) {
      return(data.frame())
    }

    if (attempt < enrichr_max_attempts) {
      Sys.sleep(
        enrichr_retry_wait_seconds * attempt
      )
    }
  }

  stop(
    "Enrichr failed after ",
    enrichr_max_attempts,
    " attempts. Last error: ",
    if (is.null(last_error)) {
      "unknown"
    } else {
      conditionMessage(last_error)
    },
    call. = FALSE
  )
}

# ==============================================================================
# 9. Helper: standardize raw Enrichr output
# ==============================================================================

standardize_enrichr_results <- function(
  raw_results,
  interaction_subset,
  cluster_id,
  cluster_name,
  n_query_genes
) {
  if (
    is.null(raw_results) ||
      nrow(raw_results) == 0L
  ) {
    return(
      empty_enrichment_table(
        analysis_status = "NO_ENRICHMENT_RESULTS",
        interaction_subset = interaction_subset,
        cluster_id = cluster_id,
        cluster_name = cluster_name,
        n_query_genes = n_query_genes
      )
    )
  }

  required_raw_columns <- c(
    "Term",
    "Overlap",
    "P.value",
    "Adjusted.P.value",
    "Odds.Ratio",
    "Combined.Score",
    "Genes"
  )

  missing_raw_columns <- setdiff(
    required_raw_columns,
    colnames(raw_results)
  )

  if (length(missing_raw_columns) > 0L) {
    stop(
      "Enrichr output is missing required column(s): ",
      paste(missing_raw_columns, collapse = ", "),
      call. = FALSE
    )
  }

  overlap_gene_lists <- strsplit(
    as.character(raw_results$Genes),
    ";",
    fixed = TRUE
  )

  overlap_gene_lists <- lapply(
    overlap_gene_lists,
    clean_gene_vector
  )

  term <- as.character(
    raw_results$Term
  )

  overlap <- as.character(
    raw_results$Overlap
  )

  n_overlap_genes <- suppressWarnings(
    as.integer(
      sub(
        "/.*$",
        "",
        overlap
      )
    )
  )

  term_gene_set_size <- suppressWarnings(
    as.integer(
      sub(
        "^.*/",
        "",
        overlap
      )
    )
  )

  missing_overlap_n <- is.na(
    n_overlap_genes
  )

  if (any(missing_overlap_n)) {
    n_overlap_genes[missing_overlap_n] <- lengths(
      overlap_gene_lists
    )[missing_overlap_n]
  }

  output <- tibble::tibble(
    term = term,
    p_value = as.numeric(raw_results$P.value),
    FDR = as.numeric(raw_results$Adjusted.P.value),
    odds_ratio = format_odds_ratio_for_excel(raw_results$Odds.Ratio),
    combined_score = as.numeric(raw_results$Combined.Score),
    overlap = overlap,
    overlap_genes = vapply(
      overlap_gene_lists,
      paste,
      collapse = "|",
      FUN.VALUE = character(1)
    ),
    n_overlap_genes = as.integer(n_overlap_genes),
    term_gene_set_size = as.integer(term_gene_set_size),
    n_query_genes = as.integer(n_query_genes),
    interaction_subset = interaction_subset,
    cluster_id = cluster_id,
    cluster_name = cluster_name,
    analysis_status = "OK"
  ) |>
    dplyr::mutate(
      significant =
        !is.na(.data$FDR) &
        .data$FDR < enrichment_fdr_threshold &
        !is.na(.data$n_overlap_genes) &
        .data$n_overlap_genes >= min_overlap_genes
    ) |>
    dplyr::arrange(
      .data$p_value,
      .data$FDR,
      dplyr::desc(.data$combined_score),
      .data$term
    ) |>
    dplyr::mutate(
      row_id = dplyr::row_number()
    ) |>
    dplyr::select(
      .data$row_id,
      .data$term,
      .data$p_value,
      .data$FDR,
      .data$odds_ratio,
      .data$combined_score,
      .data$overlap,
      .data$overlap_genes,
      .data$n_overlap_genes,
      .data$term_gene_set_size,
      .data$n_query_genes,
      .data$significant,
      .data$interaction_subset,
      .data$cluster_id,
      .data$cluster_name,
      .data$analysis_status
    )

  output
}

# ==============================================================================
# 10. Helper: run one enrichment list
# ==============================================================================

run_one_enrichment <- function(
  query_genes,
  interaction_subset,
  cluster_id,
  cluster_name,
  database
) {
  query_genes <- clean_gene_vector(
    query_genes
  )

  n_query_genes <- length(
    query_genes
  )

  if (n_query_genes == 0L) {
    return(
      empty_enrichment_table(
        analysis_status = "NO_QUERY_GENES",
        interaction_subset = interaction_subset,
        cluster_id = cluster_id,
        cluster_name = cluster_name,
        n_query_genes = n_query_genes
      )
    )
  }

  raw_results <- run_enrichr_with_retry(
    query_genes = query_genes,
    database = database
  )

  standardize_enrichr_results(
    raw_results = raw_results,
    interaction_subset = interaction_subset,
    cluster_id = cluster_id,
    cluster_name = cluster_name,
    n_query_genes = n_query_genes
  )
}

# ==============================================================================
# 11. Helper: write one XLSX workbook
# ==============================================================================

write_interaction_workbook <- function(
  interaction_subset,
  result_list,
  output_file
) {
  workbook <- openxlsx::createWorkbook()

  header_style <- openxlsx::createStyle(
    fontColour = "white",
    fgFill = "#1F4E78",
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    border = "TopBottomLeftRight",
    borderColour = "#D9D9D9",
    borderStyle = "thin"
  )

  # p-values/FDR: up to four decimal places, scientific notation below 0.0001.
  pvalue_style <- openxlsx::createStyle(
    numFmt = "[<0.0001]0.0000E+00;0.####"
  )

  two_decimal_style <- openxlsx::createStyle(
    numFmt = "0.00"
  )

  # Same pale-lavender visual highlight as in the main-effect workflow.
  highlight_style <- openxlsx::createStyle(
    fgFill = "#DEDCE6"
  )

  # Explicit borders keep highlighted cells visually separated.
  body_border_style <- openxlsx::createStyle(
    border = "TopBottomLeftRight",
    borderColour = "#D9D9D9",
    borderStyle = "thin"
  )

  sheet_keys <- c(
    "allGenes",
    expected_cluster_ids
  )

  sheet_names <- c(
    "genes_from_allClusters",
    unname(
      xlsx_cluster_sheet_labels[
        expected_cluster_ids
      ]
    )
  )

  export_columns <- c(
    "row_id",
    "term",
    "p_value",
    "FDR",
    "odds_ratio",
    "combined_score",
    "overlap",
    "overlap_genes",
    "n_overlap_genes",
    "term_gene_set_size",
    "n_query_genes"
  )

  for (sheet_index in seq_along(sheet_keys)) {
    key <- sheet_keys[[sheet_index]]
    sheet_name <- sheet_names[[sheet_index]]
    data_current <- result_list[[key]]

    missing_export_columns <- setdiff(
      export_columns,
      colnames(data_current)
    )

    if (length(missing_export_columns) > 0L) {
      stop(
        "XLSX export table is missing column(s): ",
        paste(missing_export_columns, collapse = ", "),
        call. = FALSE
      )
    }

    data_export <- data_current |>
      dplyr::select(
        dplyr::all_of(
          export_columns
        )
      )

    # In test mode, non-executed scopes are exported as header-only sheets.
    if (
      "analysis_status" %in% colnames(data_current) &&
        nrow(data_current) > 0L &&
        all(data_current$analysis_status == "NOT_RUN_TEST_MODE")
    ) {
      data_export <- data_export[
        0,
        ,
        drop = FALSE
      ]
    }

    data_export$odds_ratio <- suppressWarnings(
      as.numeric(
        data_export$odds_ratio
      )
    )

    # Excel cannot store +/-Inf as numeric cells.
    infinite_odds_rows <- which(
      is.infinite(
        data_export$odds_ratio
      )
    )

    if (length(infinite_odds_rows) > 0L) {
      data_export$odds_ratio[infinite_odds_rows] <- NA_real_
    }

    openxlsx::addWorksheet(
      workbook,
      sheetName = sheet_name,
      gridLines = TRUE
    )

    openxlsx::writeData(
      workbook,
      sheet = sheet_name,
      x = data_export,
      startRow = 1,
      startCol = 1,
      headerStyle = header_style,
      withFilter = TRUE
    )

    if (length(infinite_odds_rows) > 0L) {
      odds_ratio_col <- match(
        "odds_ratio",
        colnames(data_export)
      )

      for (row_current in infinite_odds_rows) {
        openxlsx::writeData(
          workbook,
          sheet = sheet_name,
          x = "Inf",
          startRow = row_current + 1L,
          startCol = odds_ratio_col,
          colNames = FALSE,
          rowNames = FALSE
        )
      }
    }

    openxlsx::freezePane(
      workbook,
      sheet = sheet_name,
      firstRow = TRUE
    )

    widths <- c(
      9,   # row_id
      68,  # term
      16,  # p_value
      16,  # FDR
      14,  # odds_ratio
      18,  # combined_score
      12,  # overlap
      70,  # overlap_genes
      18,  # n_overlap_genes
      20,  # term_gene_set_size
      16   # n_query_genes
    )

    openxlsx::setColWidths(
      workbook,
      sheet = sheet_name,
      cols = seq_along(widths),
      widths = widths
    )

    openxlsx::setRowHeights(
      workbook,
      sheet = sheet_name,
      rows = 1,
      heights = 28
    )

    if (nrow(data_export) > 0L) {
      pvalue_cols <- match(
        c(
          "p_value",
          "FDR"
        ),
        colnames(data_export)
      )

      two_decimal_cols <- match(
        c(
          "odds_ratio",
          "combined_score"
        ),
        colnames(data_export)
      )

      openxlsx::addStyle(
        workbook,
        sheet = sheet_name,
        style = body_border_style,
        rows = 2:(nrow(data_export) + 1L),
        cols = seq_along(export_columns),
        gridExpand = TRUE,
        stack = TRUE
      )

      openxlsx::addStyle(
        workbook,
        sheet = sheet_name,
        style = pvalue_style,
        rows = 2:(nrow(data_export) + 1L),
        cols = pvalue_cols,
        gridExpand = TRUE,
        stack = TRUE
      )

      openxlsx::addStyle(
        workbook,
        sheet = sheet_name,
        style = two_decimal_style,
        rows = 2:(nrow(data_export) + 1L),
        cols = two_decimal_cols,
        gridExpand = TRUE,
        stack = TRUE
      )

      # Visual aid only; does not filter the workbook.
      highlight_rows <- which(
        !is.na(data_current$FDR) &
          data_current$FDR < xlsx_highlight_fdr_threshold &
          !is.na(data_current$n_overlap_genes) &
          data_current$n_overlap_genes >= xlsx_highlight_min_overlap_genes
      )

      if (length(highlight_rows) > 0L) {
        openxlsx::addStyle(
          workbook,
          sheet = sheet_name,
          style = highlight_style,
          rows = highlight_rows + 1L,
          cols = seq_along(export_columns),
          gridExpand = TRUE,
          stack = TRUE
        )
      }
    }
  }

  openxlsx::saveWorkbook(
    workbook,
    output_file,
    overwrite = TRUE
  )

  if (
    !file.exists(output_file) ||
      file.info(output_file)$size == 0
  ) {
    stop(
      "XLSX file was not created correctly: ",
      output_file,
      call. = FALSE
    )
  }

  invisible(
    output_file
  )
}

# ==============================================================================
# 12. Load edgeR results
# ==============================================================================

message(
  "Loading edgeR results RData:"
)

message(
  edgeR_results_rdata_file
)

load(
  edgeR_results_rdata_file
)

required_objects <- c(
  "edgeR_perCluster_combinedResults",
  "edgeR_perCluster_testDefinitions"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1)
  )
]

if (length(missing_objects) > 0L) {
  stop(
    "Required object(s) missing after loading RData: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

# ==============================================================================
# 13. Select the interaction test and significant cluster-gene results
# ==============================================================================

interaction_test_id <- find_interaction_test_id(
  edgeR_perCluster_testDefinitions
)

message(
  "\nSelected interaction test_id: ",
  interaction_test_id
)

if (!interaction_test_id %in% names(edgeR_perCluster_combinedResults)) {
  stop(
    "Interaction test is absent from edgeR_perCluster_combinedResults: ",
    interaction_test_id,
    call. = FALSE
  )
}

interaction_test_definition <- edgeR_perCluster_testDefinitions |>
  dplyr::filter(
    .data$test_id == interaction_test_id
  )

if (nrow(interaction_test_definition) != 1L) {
  stop(
    "Selected interaction test_id does not map to exactly one test definition: ",
    interaction_test_id,
    call. = FALSE
  )
}

selected_test_comparison <- as.character(
  interaction_test_definition$comparison[[1]]
)

message(
  "Selected comparison: ",
  selected_test_comparison
)

interaction_results <- edgeR_perCluster_combinedResults[[interaction_test_id]] |>
  tibble::as_tibble()

required_interaction_columns <- c(
  "cluster_id",
  "gene",
  "logFC",
  "FDR"
)

missing_interaction_columns <- setdiff(
  required_interaction_columns,
  colnames(interaction_results)
)

if (length(missing_interaction_columns) > 0L) {
  stop(
    "Interaction result table is missing required column(s): ",
    paste(missing_interaction_columns, collapse = ", "),
    call. = FALSE
  )
}

interaction_results <- interaction_results |>
  dplyr::mutate(
    cluster_id = as.character(.data$cluster_id),
    gene_symbol = trimws(
      as.character(
        .data$gene
      )
    ),
    gene_symbol = dplyr::na_if(
      .data$gene_symbol,
      ""
    )
  )

unexpected_cluster_ids <- setdiff(
  unique(interaction_results$cluster_id),
  expected_cluster_ids
)

if (length(unexpected_cluster_ids) > 0L) {
  warning(
    "Unexpected cluster IDs in interaction results: ",
    paste(
      sort(unexpected_cluster_ids),
      collapse = ", "
    )
  )
}

significant_interaction_results <- interaction_results |>
  dplyr::filter(
    !is.na(.data$FDR),
    .data$FDR < interaction_fdr_threshold
  ) |>
  dplyr::mutate(
    interaction_direction = dplyr::case_when(
      !is.na(.data$logFC) & .data$logFC > 0 ~ "PositiveLogFC",
      !is.na(.data$logFC) & .data$logFC < 0 ~ "NegativeLogFC",
      !is.na(.data$logFC) & .data$logFC == 0 ~ "ZeroLogFC",
      TRUE ~ "MissingLogFC"
    )
  )

if (nrow(significant_interaction_results) == 0L) {
  stop(
    "No interaction result satisfies FDR < ",
    interaction_fdr_threshold,
    call. = FALSE
  )
}

message(
  "\nInteraction gene selection:"
)

message(
  "  FDR < ",
  interaction_fdr_threshold
)

message(
  "  |logFC| threshold: NONE"
)

message(
  "  significant cluster-gene results = ",
  nrow(significant_interaction_results)
)

message(
  "  unique non-missing gene symbols = ",
  significant_interaction_results |>
    dplyr::filter(
      !is.na(.data$gene_symbol)
    ) |>
    dplyr::summarise(
      n = dplyr::n_distinct(
        .data$gene_symbol
      )
    ) |>
    dplyr::pull(
      .data$n
    )
)

message(
  "  positive-logFC cluster-gene results = ",
  sum(
    significant_interaction_results$interaction_direction == "PositiveLogFC"
  )
)

message(
  "  negative-logFC cluster-gene results = ",
  sum(
    significant_interaction_results$interaction_direction == "NegativeLogFC"
  )
)

message(
  "  zero-logFC cluster-gene results = ",
  sum(
    significant_interaction_results$interaction_direction == "ZeroLogFC"
  )
)

message(
  "  missing-logFC cluster-gene results = ",
  sum(
    significant_interaction_results$interaction_direction == "MissingLogFC"
  )
)

message(
  "  significant rows without a gene symbol = ",
  sum(
    is.na(
      significant_interaction_results$gene_symbol
    )
  )
)

# ==============================================================================
# 14. Validate requested Enrichr libraries
# ==============================================================================

message(
  "\nChecking requested Enrichr libraries..."
)

enrichR::setEnrichrSite(
  "Enrichr"
)

available_enrichr_databases <- enrichR::listEnrichrDbs()

if (!"libraryName" %in% colnames(available_enrichr_databases)) {
  stop(
    "Unexpected output from enrichR::listEnrichrDbs(): libraryName column missing.",
    call. = FALSE
  )
}

missing_enrichr_libraries <- setdiff(
  enrichr_libraries$library_id,
  available_enrichr_databases$libraryName
)

if (length(missing_enrichr_libraries) > 0L) {
  message(
    "\nUnavailable requested Enrichr libraries:"
  )

  for (missing_library in missing_enrichr_libraries) {
    message(
      "  ",
      missing_library
    )

    distances <- as.numeric(
      utils::adist(
        missing_library,
        available_enrichr_databases$libraryName
      )
    )

    closest <- available_enrichr_databases$libraryName[
      order(distances)[
        seq_len(
          min(
            5L,
            length(distances)
          )
        )
      ]
    ]

    message(
      "    closest current library names: ",
      paste(
        closest,
        collapse = ", "
      )
    )
  }

  stop(
    "At least one requested Enrichr library is not currently available.",
    call. = FALSE
  )
}

message(
  "All requested Enrichr libraries are available:"
)

for (i in seq_len(nrow(enrichr_libraries))) {
  message(
    "  ",
    i,
    "/",
    nrow(enrichr_libraries),
    " | ",
    enrichr_libraries$library_id[[i]],
    " | ",
    enrichr_libraries$display_name[[i]]
  )
}

# ==============================================================================
# 15. Build the 3 sets of Enrichr query lists
# ==============================================================================

apply_interaction_subset_filter <- function(
  data,
  subset_id
) {
  if (subset_id == "ALL") {
    return(
      data
    )
  }

  if (subset_id == "PositiveLogFC") {
    return(
      data |>
        dplyr::filter(
          !is.na(.data$logFC),
          .data$logFC > 0
        )
    )
  }

  if (subset_id == "NegativeLogFC") {
    return(
      data |>
        dplyr::filter(
          !is.na(.data$logFC),
          .data$logFC < 0
        )
    )
  }

  stop(
    "Unknown interaction subset: ",
    subset_id,
    call. = FALSE
  )
}

query_by_variant <- setNames(
  vector(
    "list",
    nrow(analysis_variants)
  ),
  analysis_variants$subset_id
)

for (variant_index in seq_len(nrow(analysis_variants))) {
  subset_id <- analysis_variants$subset_id[[variant_index]]

  subset_data <- apply_interaction_subset_filter(
    significant_interaction_results,
    subset_id
  ) |>
    dplyr::filter(
      !is.na(.data$gene_symbol)
    )

  query_by_variant[[subset_id]] <- list()

  # Global list = union of unique gene symbols across all significant clusters
  # after applying the requested interaction-coefficient sign filter.
  #
  # The same gene can legitimately occur in both global PositiveLogFC and
  # NegativeLogFC workbooks if its interaction coefficient has opposite signs
  # in different clusters.
  query_by_variant[[subset_id]][["allGenes"]] <- clean_gene_vector(
    subset_data$gene_symbol
  )

  for (cluster_id_current in expected_cluster_ids) {
    cluster_data <- subset_data |>
      dplyr::filter(
        .data$cluster_id == cluster_id_current
      )

    query_by_variant[[subset_id]][[cluster_id_current]] <- clean_gene_vector(
      cluster_data$gene_symbol
    )
  }
}

message(
  "\nQuery-list sizes submitted directly to standard Enrichr:"
)

for (variant_index in seq_len(nrow(analysis_variants))) {
  subset_id <- analysis_variants$subset_id[[variant_index]]

  message(
    "  ",
    subset_id,
    " | global = ",
    length(
      query_by_variant[[subset_id]][["allGenes"]]
    )
  )
}

message(
  "  NOTE: PositiveLogFC / NegativeLogFC refer to the sign of the interaction coefficient, not UP / DOWN expression."
)

message(
  "  NOTE: standard Enrichr mode is used with no custom background and no ortholog conversion."
)

# ==============================================================================
# 16. Run one complete Enrichr database
# ==============================================================================

run_complete_database <- function(
  database,
  database_display_name,
  database_index,
  n_databases
) {
  paths <- build_database_output_paths(
    database
  )

  message(
    "\n",
    paste(
      rep("=", 88),
      collapse = ""
    )
  )

  message(
    "DATABASE ",
    database_index,
    "/",
    n_databases,
    ": ",
    database_display_name,
    " [",
    database,
    "]"
  )

  message(
    paste(
      rep("=", 88),
      collapse = ""
    )
  )

  if (
    !test_mode &&
      isTRUE(skip_completed_databases) &&
      database_outputs_complete(paths)
  ) {
    message(
      "Complete output already exists; skipping this database."
    )

    return(
      list(
        status = "SKIPPED_COMPLETE",
        database = database,
        display_name = database_display_name,
        output_dir = paths$library_output_dir,
        executed_lists = 0L
      )
    )
  }

  analysis_results <- setNames(
    vector(
      "list",
      nrow(analysis_variants)
    ),
    analysis_variants$subset_id
  )

  # Pre-create all scopes. In test mode this allows unfinished scopes to remain
  # as header-only worksheets in any workbook that was partially executed.
  for (variant_index in seq_len(nrow(analysis_variants))) {
    subset_id <- analysis_variants$subset_id[[variant_index]]

    analysis_results[[subset_id]] <- list()

    analysis_results[[subset_id]][["allGenes"]] <- empty_enrichment_table(
      analysis_status = "NOT_RUN_TEST_MODE",
      interaction_subset = subset_id,
      cluster_id = "ALL",
      cluster_name = "all clusters",
      n_query_genes = length(query_by_variant[[subset_id]][["allGenes"]])
    )

    for (cluster_id_current in expected_cluster_ids) {
      analysis_results[[subset_id]][[cluster_id_current]] <- empty_enrichment_table(
        analysis_status = "NOT_RUN_TEST_MODE",
        interaction_subset = subset_id,
        cluster_id = cluster_id_current,
        cluster_name = unname(custom_cluster_labels[[cluster_id_current]]),
        n_query_genes = length(query_by_variant[[subset_id]][[cluster_id_current]])
      )
    }
  }

  if (test_mode) {
    message(
      "TEST MODE: executing first ",
      progress_total_lists,
      "/",
      planned_total_lists,
      " lists for this database."
    )
  } else {
    message(
      "FULL MODE: executing all ",
      planned_total_lists,
      " lists for this database."
    )
  }

  database_progress <- create_database_progress(
    database = database_display_name,
    total_lists = progress_total_lists
  )

  draw_database_progress(
    database_progress,
    current_label = "starting"
  )

  executed_variant_ids <- character(0)
  executed_list_count <- 0L
  stop_running <- FALSE

  for (variant_index in seq_len(nrow(analysis_variants))) {
    if (stop_running) {
      break
    }

    subset_id <- analysis_variants$subset_id[[variant_index]]

    for (scope_key in c("allGenes", expected_cluster_ids)) {
      if (executed_list_count >= progress_total_lists) {
        stop_running <- TRUE
        break
      }

      if (identical(scope_key, "allGenes")) {
        cluster_id_current <- "ALL"
        cluster_name_current <- "all clusters"
        scope_label <- "genes_from_allClusters"
      } else {
        cluster_id_current <- scope_key
        cluster_name_current <- unname(custom_cluster_labels[[scope_key]])
        scope_label <- paste0("C", scope_key)
      }

      analysis_results[[subset_id]][[scope_key]] <- run_one_enrichment(
        query_genes = query_by_variant[[subset_id]][[scope_key]],
        interaction_subset = subset_id,
        cluster_id = cluster_id_current,
        cluster_name = cluster_name_current,
        database = database
      )

      executed_list_count <- executed_list_count + 1L

      executed_variant_ids <- unique(
        c(
          executed_variant_ids,
          subset_id
        )
      )

      advance_database_progress(
        database_progress,
        paste0(
          subset_id,
          " | ",
          scope_label
        )
      )
    }
  }

  variants_to_write <- if (test_mode) {
    analysis_variants$subset_id[
      analysis_variants$subset_id %in% executed_variant_ids
    ]
  } else {
    analysis_variants$subset_id
  }

  message(
    "\nWriting ",
    length(variants_to_write),
    if (test_mode) {
      " test XLSX workbook(s)..."
    } else {
      " XLSX workbooks..."
    }
  )

  for (subset_id in variants_to_write) {
    write_interaction_workbook(
      interaction_subset = subset_id,
      result_list = analysis_results[[subset_id]],
      output_file = paths$xlsx_files[[subset_id]]
    )
  }

  message(
    "XLSX workbook writing completed."
  )

  expected_output_files <- unname(
    paths$xlsx_files[
      variants_to_write
    ]
  )

  missing_output_files <- expected_output_files[
    !file.exists(expected_output_files)
  ]

  if (length(missing_output_files) > 0L) {
    stop(
      "Missing output file(s) for ",
      database,
      ":\n",
      paste(
        missing_output_files,
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  empty_output_files <- expected_output_files[
    file.info(expected_output_files)$size == 0
  ]

  if (length(empty_output_files) > 0L) {
    stop(
      "Empty output file(s) for ",
      database,
      ":\n",
      paste(
        empty_output_files,
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  message(
    "Database completed successfully: ",
    database_display_name
  )

  message(
    "Output directory: ",
    normalizePath(
      paths$library_output_dir,
      mustWork = TRUE
    )
  )

  list(
    status = if (test_mode) {
      "TEST_COMPLETED"
    } else {
      "COMPLETED"
    },
    database = database,
    display_name = database_display_name,
    output_dir = paths$library_output_dir,
    executed_lists = executed_list_count
  )
}

# ==============================================================================
# 17. Run all selected databases
# ==============================================================================

message(
  "\n",
  paste(
    rep("=", 88),
    collapse = ""
  )
)

message(
  "MULTI-DATABASE ENRICHR RUN: GROUP x SEX INTERACTION"
)

message(
  paste(
    rep("=", 88),
    collapse = ""
  )
)

message(
  "Databases selected: ",
  nrow(enrichr_libraries)
)

message(
  "Planned lists per database: ",
  progress_total_lists,
  if (test_mode) {
    paste0(
      " (TEST subset of ",
      planned_total_lists,
      ")"
    )
  } else {
    ""
  }
)

message(
  "Maximum list analyses in this run: ",
  progress_total_lists * nrow(enrichr_libraries)
)

message(
  "One independent progress bar will be shown for each database."
)

database_run_summary <- vector(
  "list",
  nrow(enrichr_libraries)
)

database_failures <- list()

multi_database_start_time <- Sys.time()

for (database_index in seq_len(nrow(enrichr_libraries))) {
  database <- enrichr_libraries$library_id[[database_index]]
  database_display_name <- enrichr_libraries$display_name[[database_index]]

  result_current <- tryCatch(
    run_complete_database(
      database = database,
      database_display_name = database_display_name,
      database_index = database_index,
      n_databases = nrow(enrichr_libraries)
    ),
    error = function(e) {
      database_failures[[database]] <<- conditionMessage(e)

      message(
        "\nERROR in database ",
        database_display_name,
        " [",
        database,
        "]:"
      )

      message(
        conditionMessage(e)
      )

      if (!isTRUE(continue_after_database_error)) {
        stop(e)
      }

      list(
        status = "FAILED",
        database = database,
        display_name = database_display_name,
        output_dir = build_database_output_paths(database)$library_output_dir,
        executed_lists = NA_integer_
      )
    }
  )

  database_run_summary[[database_index]] <- result_current

  invisible(
    gc()
  )
}

multi_database_elapsed <- as.numeric(
  difftime(
    Sys.time(),
    multi_database_start_time,
    units = "secs"
  )
)

database_run_summary_df <- dplyr::bind_rows(
  lapply(
    database_run_summary,
    tibble::as_tibble
  )
)

# ==============================================================================
# 18. Final report
# ==============================================================================

message(
  "\n",
  paste(
    rep("=", 88),
    collapse = ""
  )
)

message(
  "MULTI-DATABASE INTERACTION FUNCTIONAL ANALYSIS FINISHED"
)

message(
  paste(
    rep("=", 88),
    collapse = ""
  )
)

message(
  "Total elapsed time: ",
  format_elapsed_time(
    multi_database_elapsed
  )
)

message(
  "\nDatabase status:"
)

for (i in seq_len(nrow(database_run_summary_df))) {
  message(
    "  ",
    i,
    "/",
    nrow(database_run_summary_df),
    " | ",
    database_run_summary_df$display_name[[i]],
    " | ",
    database_run_summary_df$status[[i]]
  )
}

message(
  "\nInteraction selection: FDR < ",
  interaction_fdr_threshold,
  " | no |logFC| threshold"
)

message(
  "Interaction subsets: ALL | PositiveLogFC | NegativeLogFC"
)

message(
  "PositiveLogFC / NegativeLogFC describe the interaction-coefficient sign and are not UP / DOWN labels."
)

message(
  "Enrichr significant-term definition retained internally: FDR < ",
  enrichment_fdr_threshold,
  " & overlap >= ",
  min_overlap_genes,
  " genes"
)

message(
  "XLSX highlight: FDR < ",
  xlsx_highlight_fdr_threshold,
  " & overlap >= ",
  xlsx_highlight_min_overlap_genes,
  " genes | #DEDCE6"
)

message(
  "Comparison: ",
  selected_test_comparison
)

message(
  "Test ID: ",
  interaction_test_id
)

message(
  "Per completed database: 3 XLSX workbooks x 17 sheets."
)

if (length(database_failures) > 0L) {
  message(
    "\nFAILED DATABASES:"
  )

  for (database in names(database_failures)) {
    message(
      "  ",
      database,
      " -> ",
      database_failures[[database]]
    )
  }

  stop(
    length(database_failures),
    " Enrichr database(s) failed. Other databases were allowed to finish; see summary above.",
    call. = FALSE
  )
}

message(
  "\nAll requested databases completed or were already complete."
)

# ==============================================================================
# End
# ==============================================================================
