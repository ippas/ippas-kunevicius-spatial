#!/usr/bin/env Rscript

# ==============================================================================
# 01_functionalAnalysis_effectDonorGroup_multiEnrichrLibraries_v12_fullRun.R
#
# Functional enrichment analysis for the main FMT donor-group effect from the
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
# NOTE
#   Enrichr provides separate WikiPathways 2024 Human and Mouse collections.
#   This workflow uses WikiPathways_2024_Mouse because the spatial dataset is
#   mouse. Replace it with WikiPathways_2024_Human in the configuration block
#   if the human pathway collection is desired instead.
#
# DE GENE SELECTION
#   FDR < 0.05
#   abs(log2FC) >= 0.5
#
# DIRECTION VARIANTS
#   UPplusDOWN
#   UP
#   DOWN
#
# MAX MEAN % VARIANTS
#   all
#   maxPctGE25: Max mean % >= 25
#   maxPctLT25: Max mean % < 25
#
# ENRICHMENT
#   Standard Enrichr through enrichR.
#   No custom background / Speedrichr.
#   No explicit ortholog conversion.
#   Unique edgeR gene symbols are submitted directly to Enrichr.
#
# PER-DATABASE WORK
#   3 Max mean % variants x 3 directions x 17 scopes = 153 Enrichr lists.
#   One progress bar is shown independently for each database.
#
# PER-DATABASE OUTPUT
#   <functional-analysis-dir>/<library_id>/xlsx/ : 9 XLSX files
#   <functional-analysis-dir>/<library_id>/pdf/  : 4 PDF files
#
# XLSX
#   17 sheets: genes_from_allClusters and C1-C16 with short anatomical names.
#   Visual highlight only: FDR < 0.10 and overlap >= 2 genes; fill = #DEDCE6.
#
# BARPLOTS
#   Top 10 terms by smallest nominal Enrichr P-value, without FDR filtering.
#   Significant-term count in titles: FDR < 0.05 and overlap >= 3 genes.
#
# LONG-RUN SAFETY
#   Completed databases can be skipped automatically. If one database fails,
#   the remaining databases continue and failures are summarized at the end.
# ==============================================================================

options(stringsAsFactors = FALSE)

# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tibble",
  "stringr",
  "ggplot2",
  "patchwork",
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
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
  library(enrichR)
})

# ==============================================================================
# 2. Analysis configuration
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"
clustering_name <- "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

de_fdr_threshold <- 0.05
de_abs_log2fc_threshold <- 0.5
max_mean_percent_threshold <- 25
enrichment_fdr_threshold <- 0.05
min_overlap_genes <- 3L
xlsx_highlight_fdr_threshold <- 0.10
xlsx_highlight_min_overlap_genes <- 2L
global_top_n <- 10L
cluster_top_n <- 10L

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

# NULL = all 153 lists per database. Integer = first N lists per database.
test_first_n_lists <- NULL

# Skip a database when all expected 9 XLSX + 4 PDF files already exist and
# are non-empty. This is useful when rerunning a long multi-database job.
skip_completed_databases <- TRUE

# Continue with later databases if one database fails. At the very end the
# script stops with a failure summary if any database failed.
continue_after_database_error <- TRUE

enrichr_sleep_seconds <- 1
enrichr_max_attempts <- 3L
enrichr_retry_wait_seconds <- 3

parameter_block <- "FDR005absLog2FC05perCluster"
analysis_directions <- c("UPplusDOWN", "UP", "DOWN")
analysis_subsets <- c("all", "maxPctGE25", "maxPctLT25")

analysis_variants <- data.frame(
  xlsx_index = 1:9,
  subset_id = rep(analysis_subsets, each = length(analysis_directions)),
  direction = rep(analysis_directions, times = length(analysis_subsets)),
  stringsAsFactors = FALSE
)

analysis_variants$variant_id <- paste(analysis_variants$subset_id, analysis_variants$direction, sep = "__")
analysis_variants$file_suffix <- ifelse(
  analysis_variants$subset_id == "all",
  paste0(analysis_variants$direction, "_enrichrResults.xlsx"),
  paste0(analysis_variants$subset_id, "_", analysis_variants$direction, "_enrichrResults.xlsx")
)
analysis_variants$subset_label <- c(
  rep("no Max mean % filter", 3),
  rep("Max mean % >= 25%", 3),
  rep("Max mean % < 25%", 3)
)

planned_total_lists <- nrow(analysis_variants) * 17L

if (!is.null(test_first_n_lists)) {
  if (length(test_first_n_lists) != 1L || is.na(test_first_n_lists) || !is.finite(test_first_n_lists) || test_first_n_lists < 1 || test_first_n_lists != as.integer(test_first_n_lists)) {
    stop("`test_first_n_lists` must be NULL or one positive integer.", call. = FALSE)
  }
  test_first_n_lists <- as.integer(test_first_n_lists)
}

test_mode <- !is.null(test_first_n_lists)
progress_total_lists <- if (test_mode) min(test_first_n_lists, planned_total_lists) else planned_total_lists

# ==============================================================================
# 3. Cluster colours and anatomical labels
# ==============================================================================

custom_cluster_colors <- c(
  "1"  = "#e66063",
  "2"  = "#407ba7",
  "3"  = "#31cb00",
  "4"  = "#adb5bd",
  "5"  = "#9e0059",
  "6"  = "#b5838d",
  "7"  = "#002962",
  "8"  = "#004e89",
  "9"  = "#d02224",
  "10" = "#ff9505",
  "11" = "#9c191b",
  "12" = "#e85d04",
  "13" = "#119822",
  "14" = "#6B1E2D",
  "15" = "#ffd100",
  "16" = "#1e441e"
)

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

# Short anatomical labels used only in XLSX worksheet names.
# Microsoft Excel/openxlsx limits worksheet names to 31 characters.
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

# Global plots do not represent one anatomical cluster, so direction-specific
# colours are used only in the global overview.
global_direction_colors <- c(
  "UPplusDOWN" = "grey35",
  "UP" = "#B2182B",
  "DOWN" = "#2166AC"
)

# ==============================================================================
# 4. Input and common output paths
# ==============================================================================

statistics_dir <- file.path(
  project_root, "results", dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name, "01_statistics"
)

analysis_prefix <- paste0(
  dataset_name, "_", clustering_name,
  "_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction"
)

edgeR_results_rdata_file <- file.path(
  statistics_dir, "03_edgeRResults", paste0(analysis_prefix, "_edgeRResults.RData")
)

functional_analysis_root <- file.path(
  project_root, "results", dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name, "07_functionalAnalysis", "02_mainEffects_donorGroup"
)

dir.create(functional_analysis_root, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(edgeR_results_rdata_file)) {
  stop("Missing edgeR results RData file:\n", edgeR_results_rdata_file, call. = FALSE)
}

build_database_output_paths <- function(database) {
  library_output_dir <- file.path(functional_analysis_root, database)
  if (test_mode) library_output_dir <- file.path(library_output_dir, paste0("TEST_first", progress_total_lists))

  xlsx_output_dir <- file.path(library_output_dir, "xlsx")
  pdf_output_dir <- file.path(library_output_dir, "pdf")
  dir.create(xlsx_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pdf_output_dir, recursive = TRUE, showWarnings = FALSE)

  output_prefix <- paste0(dataset_name, "_effectDonorGroup_", parameter_block, "_", database)

  xlsx_files <- setNames(
    vapply(seq_len(nrow(analysis_variants)), function(i) {
      file.path(
        xlsx_output_dir,
        paste0(sprintf("%02d", analysis_variants$xlsx_index[[i]]), "_", output_prefix, "_", analysis_variants$file_suffix[[i]])
      )
    }, FUN.VALUE = character(1)),
    analysis_variants$variant_id
  )

  pdf_files <- c(
    global = file.path(pdf_output_dir, paste0("01_", output_prefix, "_global_UPplusDOWN_UP_DOWN_top10.pdf")),
    UPplusDOWN = file.path(pdf_output_dir, paste0("02_", output_prefix, "_clusters_UPplusDOWN_top10_grid4x4.pdf")),
    UP = file.path(pdf_output_dir, paste0("03_", output_prefix, "_clusters_UP_top10_grid4x4.pdf")),
    DOWN = file.path(pdf_output_dir, paste0("04_", output_prefix, "_clusters_DOWN_top10_grid4x4.pdf"))
  )

  list(
    library_output_dir = library_output_dir,
    xlsx_output_dir = xlsx_output_dir,
    pdf_output_dir = pdf_output_dir,
    output_prefix = output_prefix,
    xlsx_files = xlsx_files,
    pdf_files = pdf_files
  )
}

database_outputs_complete <- function(paths) {
  if (test_mode) return(FALSE)
  expected <- c(unname(paths$xlsx_files), unname(paths$pdf_files))
  if (!all(file.exists(expected))) return(FALSE)
  sizes <- file.info(expected)$size
  all(!is.na(sizes) & sizes > 0)
}

# ==============================================================================
# 5. Helper: identify the donor-group main-effect test
# ==============================================================================

find_donor_group_test_id <- function(test_definitions) {
  required_columns <- c("test_id", "sheet_name", "comparison")
  missing_columns <- setdiff(required_columns, colnames(test_definitions))

  if (length(missing_columns) > 0L) {
    stop(
      "Test-definition table is missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  test_table <- test_definitions |>
    dplyr::mutate(
      test_id_chr = as.character(.data$test_id),
      sheet_name_chr = as.character(.data$sheet_name),
      comparison_chr = as.character(.data$comparison),
      test_id_lower = tolower(.data$test_id_chr),
      sheet_name_lower = tolower(.data$sheet_name_chr),
      comparison_lower = tolower(.data$comparison_chr)
    )

  eligible_rows <- test_table |>
    dplyr::filter(
      !grepl("interaction", .data$test_id_lower),
      !grepl("interaction", .data$sheet_name_lower),
      !grepl("sex", .data$test_id_lower),
      !grepl("sex", .data$sheet_name_lower),
      !grepl("male", .data$test_id_lower),
      !grepl("male", .data$sheet_name_lower),
      !grepl("female", .data$test_id_lower),
      !grepl("female", .data$sheet_name_lower),
      !grepl("simple", .data$test_id_lower),
      !grepl("simple", .data$sheet_name_lower)
    )

  scored_rows <- eligible_rows |>
    dplyr::mutate(
      donor_test_score =
        5L * (grepl("donor", .data$test_id_lower) | grepl("fmt", .data$test_id_lower)) +
        4L * (grepl("donor", .data$sheet_name_lower) | grepl("fmt", .data$sheet_name_lower)) +
        3L * grepl("group", .data$test_id_lower) +
        2L * grepl("group", .data$sheet_name_lower) +
        3L * (grepl("asd", .data$comparison_lower) & grepl("neurotypical", .data$comparison_lower)) +
        1L * (grepl("asd", .data$test_id_lower) | grepl("neurotypical", .data$test_id_lower)) +
        5L * (grepl("overall", .data$test_id_lower) | grepl("overall", .data$sheet_name_lower))
    ) |>
    dplyr::filter(.data$donor_test_score > 0L)

  if (nrow(scored_rows) == 0L) {
    message("\nAll available edgeR test definitions:")
    print(
      test_definitions |>
        dplyr::select(.data$test_id, .data$sheet_name, .data$comparison),
      n = Inf
    )
    stop("Could not identify any donor-group main-effect candidate.", call. = FALSE)
  }

  best_score <- max(scored_rows$donor_test_score)
  best_rows <- scored_rows |>
    dplyr::filter(.data$donor_test_score == best_score)

  if (nrow(best_rows) != 1L) {
    message("\nBest donor-group test candidates:")
    print(
      best_rows |>
        dplyr::select(
          .data$test_id,
          .data$sheet_name,
          .data$comparison,
          .data$donor_test_score
        ),
      n = Inf
    )
    stop(
      "Could not uniquely identify the donor-group main-effect test. Best-score candidates: ",
      nrow(best_rows),
      call. = FALSE
    )
  }

  selected_test_id <- as.character(best_rows$test_id[[1]])
  selected_test_id_lower <- tolower(selected_test_id)
  selected_sheet_lower <- tolower(as.character(best_rows$sheet_name[[1]]))

  if (
    grepl("sex", selected_test_id_lower) ||
    grepl("sex", selected_sheet_lower) ||
    grepl("interaction", selected_test_id_lower) ||
    grepl("interaction", selected_sheet_lower)
  ) {
    stop(
      "Safety check failed: selected donor-group test appears to be a sex or interaction test: ",
      selected_test_id,
      call. = FALSE
    )
  }

  selected_test_id
}

# ==============================================================================
# 6. Helper: general formatting and empty-result handling
# ==============================================================================

clean_gene_vector <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x)]
  sort(unique(x))
}

empty_enrichment_table <- function(
  analysis_status,
  direction,
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
    direction = direction,
    cluster_id = cluster_id,
    cluster_name = cluster_name,
    analysis_status = analysis_status
  )
}

format_odds_ratio_for_excel <- function(x) {
  suppressWarnings(as.numeric(x))
}

# ==============================================================================
# 7. Helper: one database-wide progress bar and standard Enrichr request
# ============================================================================== 

format_elapsed_time <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) {
    return("--:--:--")
  }

  seconds <- as.integer(round(seconds))
  hours <- seconds %/% 3600L
  minutes <- (seconds %% 3600L) %/% 60L
  secs <- seconds %% 60L

  sprintf("%02d:%02d:%02d", hours, minutes, secs)
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

draw_database_progress <- function(progress, current_label = NULL, finish = FALSE) {
  if (!is.null(current_label)) {
    progress$last_label <- current_label
  }

  completed <- progress$completed
  total <- progress$total
  fraction <- if (total > 0L) completed / total else 1
  fraction <- min(max(fraction, 0), 1)

  filled <- floor(progress$width * fraction)
  bar <- paste0(
    "[",
    paste(rep("=", filled), collapse = ""),
    paste(rep(" ", progress$width - filled), collapse = ""),
    "]"
  )

  elapsed <- as.numeric(difftime(Sys.time(), progress$start_time, units = "secs"))
  mean_per_list <- if (completed > 0L) elapsed / completed else NA_real_
  eta <- if (completed > 0L && completed < total) {
    mean_per_list * (total - completed)
  } else if (completed >= total) {
    0
  } else {
    NA_real_
  }

  avg_label <- if (is.finite(mean_per_list)) {
    sprintf("%.1fs/list", mean_per_list)
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

  # Pad the line so a shorter subsequent label fully overwrites the previous one.
  cat(sprintf("%-220s", line))
  flush.console()

  if (finish) {
    cat("\n")
  }

  invisible(progress)
}

advance_database_progress <- function(progress, current_label) {
  progress$completed <- progress$completed + 1L
  draw_database_progress(
    progress = progress,
    current_label = current_label,
    finish = progress$completed >= progress$total
  )
}

run_enrichr_with_retry <- function(query_genes, database) {
  last_error <- NULL

  for (attempt in seq_len(enrichr_max_attempts)) {
    result <- NULL

    # enrichR prints several status lines for every request. Capture them so the
    # terminal shows one persistent database-wide progress bar instead.
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

    if (!is.null(result) && database %in% names(result)) {
      database_result <- result[[database]]

      if (is.null(database_result) || nrow(database_result) == 0L) {
        return(data.frame())
      }

      return(database_result)
    }

    # Some enrichR versions can throw a column-name assignment error when the
    # server returns zero rows. Treat that specific parser failure as a valid
    # empty enrichment result rather than retrying the same empty list 3 times.
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
      Sys.sleep(enrichr_retry_wait_seconds * attempt)
    }
  }

  stop(
    "Enrichr failed after ",
    enrichr_max_attempts,
    " attempts. Last error: ",
    if (is.null(last_error)) "unknown" else conditionMessage(last_error),
    call. = FALSE
  )
}

# ==============================================================================
# 8. Helper: convert raw Enrichr output to the standardized XLSX table
# ==============================================================================

standardize_enrichr_results <- function(
  raw_results,
  direction,
  cluster_id,
  cluster_name,
  n_query_genes
) {
  if (is.null(raw_results) || nrow(raw_results) == 0L) {
    return(
      empty_enrichment_table(
        analysis_status = "NO_ENRICHMENT_RESULTS",
        direction = direction,
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

  missing_raw_columns <- setdiff(required_raw_columns, colnames(raw_results))
  if (length(missing_raw_columns) > 0L) {
    stop(
      "Enrichr output is missing required column(s): ",
      paste(missing_raw_columns, collapse = ", "),
      call. = FALSE
    )
  }

  overlap_gene_lists <- strsplit(as.character(raw_results$Genes), ";", fixed = TRUE)
  overlap_gene_lists <- lapply(overlap_gene_lists, clean_gene_vector)

  # Keep the Enrichr term exactly as returned by the selected library.
  # No GO-ID or other database-specific identifier extraction is performed.
  term <- as.character(raw_results$Term)

  overlap <- as.character(raw_results$Overlap)
  n_overlap_genes <- suppressWarnings(as.integer(sub("/.*$", "", overlap)))
  term_gene_set_size <- suppressWarnings(as.integer(sub("^.*/", "", overlap)))

  # Fallback only if Enrichr returns an unexpected Overlap representation.
  missing_overlap_n <- is.na(n_overlap_genes)
  if (any(missing_overlap_n)) {
    n_overlap_genes[missing_overlap_n] <- lengths(overlap_gene_lists)[missing_overlap_n]
  }

  output <- tibble::tibble(
    term = term,
    p_value = as.numeric(raw_results$P.value),
    FDR = as.numeric(raw_results$Adjusted.P.value),
    odds_ratio = format_odds_ratio_for_excel(raw_results$Odds.Ratio),
    combined_score = as.numeric(raw_results$Combined.Score),
    overlap = overlap,
    overlap_genes = vapply(overlap_gene_lists, paste, collapse = "|", FUN.VALUE = character(1)),
    n_overlap_genes = as.integer(n_overlap_genes),
    term_gene_set_size = as.integer(term_gene_set_size),
    n_query_genes = as.integer(n_query_genes),
    direction = direction,
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
    dplyr::mutate(row_id = dplyr::row_number()) |>
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
      .data$direction,
      .data$cluster_id,
      .data$cluster_name,
      .data$analysis_status
    )

  output
}

# ==============================================================================
# 9. Helper: build one analysis unit
# ==============================================================================

run_one_enrichment <- function(query_genes, direction, cluster_id, cluster_name, database) {
  query_genes <- clean_gene_vector(query_genes)
  n_query_genes <- length(query_genes)

  if (n_query_genes == 0L) {
    return(empty_enrichment_table(
      analysis_status = "NO_QUERY_GENES",
      direction = direction,
      cluster_id = cluster_id,
      cluster_name = cluster_name,
      n_query_genes = n_query_genes
    ))
  }

  raw_results <- run_enrichr_with_retry(query_genes = query_genes, database = database)

  standardize_enrichr_results(
    raw_results = raw_results,
    direction = direction,
    cluster_id = cluster_id,
    cluster_name = cluster_name,
    n_query_genes = n_query_genes
  )
}

# ==============================================================================
# 10. Helper: XLSX writer
# ==============================================================================

write_direction_workbook <- function(direction, result_list, output_file) {
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

  # Display p-values/FDR with up to four decimal places.
  # Values below 0.0001 are displayed in scientific notation with four
  # decimal places in the mantissa while the underlying cells remain numeric.
  pvalue_style <- openxlsx::createStyle(
    numFmt = "[<0.0001]0.0000E+00;0.####"
  )

  two_decimal_style <- openxlsx::createStyle(
    numFmt = "0.00"
  )

  # Very pale lavender fill for potentially interesting rows. This is visual only.
  highlight_style <- openxlsx::createStyle(
    fgFill = "#DEDCE6"
  )

  # Explicit body borders are required because spreadsheet gridlines may become
  # visually hidden when a fill colour is applied. Light-grey thin borders keep
  # all rows and columns clearly separated, including highlighted rows.
  body_border_style <- openxlsx::createStyle(
    border = "TopBottomLeftRight",
    borderColour = "#D9D9D9",
    borderStyle = "thin"
  )

  sheet_keys <- c("allGenes", expected_cluster_ids)
  sheet_names <- c(
    "genes_from_allClusters",
    unname(xlsx_cluster_sheet_labels[expected_cluster_ids])
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

    missing_export_columns <- setdiff(export_columns, colnames(data_current))
    if (length(missing_export_columns) > 0L) {
      stop(
        "XLSX export table is missing column(s): ",
        paste(missing_export_columns, collapse = ", "),
        call. = FALSE
      )
    }

    data_export <- data_current |>
      dplyr::select(dplyr::all_of(export_columns))

    # In development test mode, scopes that were not executed are retained in
    # analysis_results only as internal placeholders. Export them as header-only
    # worksheets rather than as a visually confusing row of NA values.
    if (
      "analysis_status" %in% colnames(data_current) &&
        nrow(data_current) > 0L &&
        all(data_current$analysis_status == "NOT_RUN_TEST_MODE")
    ) {
      data_export <- data_export[0, , drop = FALSE]
    }

    # Coerce Odds Ratio back to numeric here as well so this XLSX-writing
    # section can be rerun on analysis_results created by an older script
    # version where odds_ratio was stored as character text.
    data_export$odds_ratio <- suppressWarnings(
      as.numeric(data_export$odds_ratio)
    )

    # Keep Odds Ratio numeric for sorting/filtering. Excel/LibreOffice cannot
    # store +/-Inf as a normal numeric cell, so only those rare cells are
    # temporarily written as NA and then replaced by the literal text "Inf".
    infinite_odds_rows <- which(is.infinite(data_export$odds_ratio))
    if (length(infinite_odds_rows) > 0L) {
      data_export$odds_ratio[infinite_odds_rows] <- NA_real_
    }

    openxlsx::addWorksheet(workbook, sheetName = sheet_name, gridLines = TRUE)
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
      odds_ratio_col <- match("odds_ratio", colnames(data_export))

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

    openxlsx::freezePane(workbook, sheet = sheet_name, firstRow = TRUE)

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
      pvalue_cols <- match(c("p_value", "FDR"), colnames(data_export))
      two_decimal_cols <- match(
        c("odds_ratio", "combined_score"),
        colnames(data_export)
      )

      # Apply explicit borders to the complete worksheet body before any
      # number formats or highlighting. stack = TRUE below ensures that later
      # styles preserve these borders.
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

      # Highlight complete result rows when FDR < 0.10 and at least two genes
      # overlap with the Enrichr term. The criterion is intentionally only a
      # spreadsheet visual aid and is not used for filtering or plotting.
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

  openxlsx::saveWorkbook(workbook, output_file, overwrite = TRUE)

  if (!file.exists(output_file) || file.info(output_file)$size == 0) {
    stop("XLSX file was not created correctly: ", output_file, call. = FALSE)
  }
}

# ==============================================================================
# 11. Helper: barplots
# ==============================================================================

prepare_plot_data <- function(enrichment_table, top_n, wrap_width = 44L) {
  if (is.null(enrichment_table) || nrow(enrichment_table) == 0L) {
    return(tibble::tibble())
  }

  # IMPORTANT:
  # The displayed Top N is selected from ALL Enrichr terms with a valid
  # nominal P-value. There is deliberately NO FDR filter here and NO minimum
  # overlap filter here. FDR < 0.05 and overlap >= 3 genes are used only to
  # calculate the number of statistically significant terms reported in titles.
  enrichment_table |>
    dplyr::filter(
      .data$analysis_status == "OK",
      !is.na(.data$p_value),
      is.finite(.data$p_value),
      .data$p_value >= 0
    ) |>
    dplyr::arrange(
      .data$p_value,
      .data$FDR,
      dplyr::desc(.data$combined_score),
      .data$term
    ) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(
      minus_log10_p = -log10(pmax(.data$p_value, .Machine$double.xmin)),
      term_plot = paste0(.data$term, " (n = ", .data$n_overlap_genes, ")"),
      term_plot = stringr::str_wrap(.data$term_plot, width = wrap_width),
      term_plot = factor(.data$term_plot, levels = rev(.data$term_plot))
    )
}

get_result_counts <- function(enrichment_table) {
  if (is.null(enrichment_table) || nrow(enrichment_table) == 0L) {
    return(list(
      status = "NO_RESULTS",
      n_query = 0L,
      n_all_terms = 0L,
      n_sig_terms = 0L
    ))
  }

  status <- as.character(enrichment_table$analysis_status[[1]])
  n_query <- suppressWarnings(as.integer(enrichment_table$n_query_genes[[1]]))
  n_all_terms <- sum(
    enrichment_table$analysis_status == "OK" &
      !is.na(enrichment_table$p_value),
    na.rm = TRUE
  )
  n_sig_terms <- sum(
    enrichment_table$analysis_status == "OK" & enrichment_table$significant,
    na.rm = TRUE
  )

  list(
    status = status,
    n_query = ifelse(is.na(n_query), 0L, n_query),
    n_all_terms = n_all_terms,
    n_sig_terms = n_sig_terms
  )
}

make_empty_plot <- function(title, message_text, title_size = 11, message_size = 4.2) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = message_text,
      size = message_size,
      lineheight = 1.15,
      fontface = "plain"
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = title_size,
        hjust = 0,
        lineheight = 1.05
      ),
      plot.margin = ggplot2::margin(7, 7, 7, 7)
    )
}

make_global_barplot <- function(enrichment_table, direction) {
  counts <- get_result_counts(enrichment_table)
  plot_data <- prepare_plot_data(
    enrichment_table = enrichment_table,
    top_n = global_top_n,
    wrap_width = 58L
  )

  direction_label <- switch(
    direction,
    "UPplusDOWN" = "UP + DOWN",
    "UP" = "UP",
    "DOWN" = "DOWN",
    direction
  )

  title_text <- paste0(
    direction_label,
    " | all clusters | ",
    counts$n_query,
    " query genes | ",
    counts$n_sig_terms,
    " significant terms"
  )

  if (nrow(plot_data) == 0L) {
    message_text <- if (counts$status != "OK") {
      paste0("Analysis status: ", counts$status)
    } else {
      "No enrichment terms with a valid nominal P-value were returned."
    }

    return(
      make_empty_plot(
        title = title_text,
        message_text = message_text,
        title_size = 13,
        message_size = 4.5
      )
    )
  }

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$minus_log10_p, y = .data$term_plot)
  ) +
    ggplot2::geom_col(
      fill = unname(global_direction_colors[[direction]]),
      width = 0.72
    ) +
    ggplot2::labs(
      title = title_text,
      x = expression(-log[10](italic(P) * "-value")),
      y = NULL
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 13,
        lineheight = 1.05
      ),
      axis.text.y = ggplot2::element_text(
        size = 10.5,
        lineheight = 0.95,
        colour = "black"
      ),
      axis.text.x = ggplot2::element_text(size = 9.5, colour = "black"),
      axis.title.x = ggplot2::element_text(size = 11, face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(9, 10, 9, 10)
    )
}

make_cluster_barplot <- function(enrichment_table, cluster_id, direction) {
  counts <- get_result_counts(enrichment_table)
  plot_data <- prepare_plot_data(
    enrichment_table = enrichment_table,
    top_n = cluster_top_n,
    wrap_width = 42L
  )

  cluster_name <- unname(custom_cluster_labels[[cluster_id]])
  title_text <- paste0(
    "C",
    cluster_id,
    " | ",
    cluster_name,
    "\n",
    counts$n_query,
    " query genes | ",
    counts$n_sig_terms,
    " significant terms"
  )

  if (nrow(plot_data) == 0L) {
    message_text <- if (counts$status != "OK") {
      paste0("Status: ", counts$status)
    } else {
      "No enrichment terms with a valid nominal P-value were returned."
    }

    return(
      make_empty_plot(
        title = title_text,
        message_text = message_text,
        title_size = 9.4,
        message_size = 3.8
      )
    )
  }

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$minus_log10_p, y = .data$term_plot)
  ) +
    ggplot2::geom_col(
      fill = unname(custom_cluster_colors[[cluster_id]]),
      width = 0.72
    ) +
    ggplot2::labs(
      title = title_text,
      x = NULL,
      y = NULL
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 9.4,
        lineheight = 1.05,
        colour = "black"
      ),
      axis.text.y = ggplot2::element_text(
        size = 7.6,
        lineheight = 0.92,
        colour = "black"
      ),
      axis.text.x = ggplot2::element_text(
        size = 7.5,
        colour = "black"
      ),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(6, 7, 6, 7)
    )
}

# ==============================================================================
# 12. Load edgeR results
# ==============================================================================

message("Loading edgeR results RData:")
message(edgeR_results_rdata_file)
load(edgeR_results_rdata_file)

required_objects <- c(
  "edgeR_perCluster_combinedResults",
  "edgeR_perCluster_testDefinitions"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1))
]

if (length(missing_objects) > 0L) {
  stop(
    "Required object(s) missing after loading RData: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

# ==============================================================================
# 13. Select donor-group main effect and prepare all tested genes
# ==============================================================================

donor_group_test_id <- find_donor_group_test_id(edgeR_perCluster_testDefinitions)

message("\nSelected donor-group test_id: ", donor_group_test_id)

if (!donor_group_test_id %in% names(edgeR_perCluster_combinedResults)) {
  stop(
    "Selected donor-group test is absent from edgeR_perCluster_combinedResults: ",
    donor_group_test_id,
    call. = FALSE
  )
}

donor_group_test_definition <- edgeR_perCluster_testDefinitions |>
  dplyr::filter(.data$test_id == donor_group_test_id)

if (nrow(donor_group_test_definition) != 1L) {
  stop(
    "Selected donor-group test_id does not map to exactly one test definition: ",
    donor_group_test_id,
    call. = FALSE
  )
}

selected_test_comparison <- as.character(donor_group_test_definition$comparison[[1]])
message("Selected comparison: ", selected_test_comparison)

donor_group_results <- edgeR_perCluster_combinedResults[[donor_group_test_id]] |>
  tibble::as_tibble()

required_result_columns <- c(
  "cluster_id",
  "ensembl_gene_id",
  "gene",
  "logFC",
  "PValue",
  "FDR",
  "mean_percent_positive_spots_Neurotypical",
  "mean_percent_positive_spots_ASD"
)

missing_result_columns <- setdiff(required_result_columns, colnames(donor_group_results))

if (length(missing_result_columns) > 0L) {
  stop(
    "Donor-group result table is missing required column(s): ",
    paste(missing_result_columns, collapse = ", "),
    call. = FALSE
  )
}

donor_group_results <- donor_group_results |>
  dplyr::mutate(
    cluster_id = as.character(.data$cluster_id),
    gene_symbol = trimws(as.character(.data$gene)),
    gene_symbol = dplyr::na_if(.data$gene_symbol, ""),
    max_mean_percent = pmax(
      .data$mean_percent_positive_spots_Neurotypical,
      .data$mean_percent_positive_spots_ASD,
      na.rm = TRUE
    ),
    max_mean_percent = dplyr::if_else(
      is.infinite(.data$max_mean_percent),
      NA_real_,
      .data$max_mean_percent
    )
  )

unexpected_cluster_ids <- setdiff(unique(donor_group_results$cluster_id), expected_cluster_ids)
if (length(unexpected_cluster_ids) > 0L) {
  warning(
    "Unexpected cluster IDs in donor-group results: ",
    paste(sort(unexpected_cluster_ids), collapse = ", ")
  )
}

significant_donor_group_results <- donor_group_results |>
  dplyr::filter(
    !is.na(.data$FDR),
    !is.na(.data$logFC),
    .data$FDR < de_fdr_threshold,
    abs(.data$logFC) >= de_abs_log2fc_threshold
  )

message("\nDE selection:")
message("  FDR < ", de_fdr_threshold)
message("  |log2FC| >= ", de_abs_log2fc_threshold)
message("  cluster-gene results = ", nrow(significant_donor_group_results))
message(
  "  unique non-missing gene symbols = ",
  significant_donor_group_results |>
    dplyr::filter(!is.na(.data$gene_symbol)) |>
    dplyr::summarise(n = dplyr::n_distinct(.data$gene_symbol)) |>
    dplyr::pull(.data$n)
)
message(
  "  significant rows without a gene symbol = ",
  sum(is.na(significant_donor_group_results$gene_symbol))
)

# ==============================================================================
# 14. Validate all requested Enrichr libraries
# ==============================================================================

message("\nChecking requested Enrichr libraries...")
enrichR::setEnrichrSite("Enrichr")
available_enrichr_databases <- enrichR::listEnrichrDbs()

if (!"libraryName" %in% colnames(available_enrichr_databases)) {
  stop("Unexpected output from enrichR::listEnrichrDbs(): libraryName column missing.", call. = FALSE)
}

missing_enrichr_libraries <- setdiff(enrichr_libraries$library_id, available_enrichr_databases$libraryName)

if (length(missing_enrichr_libraries) > 0L) {
  message("\nUnavailable requested Enrichr libraries:")
  for (missing_library in missing_enrichr_libraries) {
    message("  ", missing_library)
    distances <- as.numeric(utils::adist(missing_library, available_enrichr_databases$libraryName))
    closest <- available_enrichr_databases$libraryName[order(distances)[seq_len(min(5L, length(distances)))]]
    message("    closest current library names: ", paste(closest, collapse = ", "))
  }
  stop("At least one requested Enrichr library is not currently available.", call. = FALSE)
}

message("All requested Enrichr libraries are available:")
for (i in seq_len(nrow(enrichr_libraries))) {
  message("  ", i, "/", nrow(enrichr_libraries), " | ", enrichr_libraries$library_id[[i]], " | ", enrichr_libraries$display_name[[i]])
}

# ==============================================================================
# 15. Build the 9 sets of Enrichr query lists
# ============================================================================== 

apply_direction_filter <- function(data, direction) {
  if (direction == "UP") {
    return(data |> dplyr::filter(.data$logFC > 0))
  }

  if (direction == "DOWN") {
    return(data |> dplyr::filter(.data$logFC < 0))
  }

  if (direction == "UPplusDOWN") {
    return(data)
  }

  stop("Unknown direction: ", direction, call. = FALSE)
}

apply_cluster_subset_filter <- function(data, subset_id) {
  if (subset_id == "all") {
    return(data)
  }

  if (subset_id == "maxPctGE25") {
    return(
      data |>
        dplyr::filter(
          !is.na(.data$max_mean_percent),
          .data$max_mean_percent >= max_mean_percent_threshold
        )
    )
  }

  if (subset_id == "maxPctLT25") {
    return(
      data |>
        dplyr::filter(
          !is.na(.data$max_mean_percent),
          .data$max_mean_percent < max_mean_percent_threshold
        )
    )
  }

  stop("Unknown subset_id: ", subset_id, call. = FALSE)
}

build_global_query <- function(data_direction, subset_id) {
  data_direction <- data_direction |>
    dplyr::filter(!is.na(.data$gene_symbol))

  if (subset_id == "all") {
    return(clean_gene_vector(data_direction$gene_symbol))
  }

  # One gene can be significant in several clusters. For a global Max mean %
  # split, first collapse to one row per symbol using the maximum Max mean %
  # observed among its significant cluster-gene occurrences. This guarantees
  # complementary global >=25 and <25 lists within the same direction.
  gene_level <- data_direction |>
    dplyr::group_by(.data$gene_symbol) |>
    dplyr::summarise(
      global_max_mean_percent = if (all(is.na(.data$max_mean_percent))) {
        NA_real_
      } else {
        max(.data$max_mean_percent, na.rm = TRUE)
      },
      .groups = "drop"
    )

  if (subset_id == "maxPctGE25") {
    gene_level <- gene_level |>
      dplyr::filter(
        !is.na(.data$global_max_mean_percent),
        .data$global_max_mean_percent >= max_mean_percent_threshold
      )
  } else if (subset_id == "maxPctLT25") {
    gene_level <- gene_level |>
      dplyr::filter(
        !is.na(.data$global_max_mean_percent),
        .data$global_max_mean_percent < max_mean_percent_threshold
      )
  } else {
    stop("Unknown subset_id: ", subset_id, call. = FALSE)
  }

  clean_gene_vector(gene_level$gene_symbol)
}

query_by_variant <- setNames(
  vector("list", nrow(analysis_variants)),
  analysis_variants$variant_id
)

for (variant_index in seq_len(nrow(analysis_variants))) {
  variant <- analysis_variants[variant_index, , drop = FALSE]
  variant_id <- variant$variant_id[[1]]
  subset_id <- variant$subset_id[[1]]
  direction <- variant$direction[[1]]

  direction_data <- apply_direction_filter(
    significant_donor_group_results,
    direction
  )

  query_by_variant[[variant_id]] <- list()
  query_by_variant[[variant_id]][["allGenes"]] <- build_global_query(
    data_direction = direction_data,
    subset_id = subset_id
  )

  for (cluster_id_current in expected_cluster_ids) {
    cluster_data <- direction_data |>
      dplyr::filter(.data$cluster_id == cluster_id_current)

    cluster_data <- apply_cluster_subset_filter(
      data = cluster_data,
      subset_id = subset_id
    )

    query_by_variant[[variant_id]][[cluster_id_current]] <- clean_gene_vector(
      cluster_data$gene_symbol
    )
  }
}

message("\nQuery-list sizes submitted directly to standard Enrichr:")
for (variant_index in seq_len(nrow(analysis_variants))) {
  variant <- analysis_variants[variant_index, , drop = FALSE]
  variant_id <- variant$variant_id[[1]]

  message(
    "  ",
    variant$subset_id[[1]],
    " | ",
    variant$direction[[1]],
    " | global = ",
    length(query_by_variant[[variant_id]][["allGenes"]])
  )
}

message(
  "  NOTE: standard Enrichr mode is used with no custom background and no ",
  "ortholog conversion."
)

# ==============================================================================
# 16. Run one complete Enrichr database
# ==============================================================================

unfiltered_variant_id <- function(direction) paste("all", direction, sep = "__")

run_complete_database <- function(database, database_display_name, database_index, n_databases) {
  paths <- build_database_output_paths(database)

  message("\n", paste(rep("=", 88), collapse = ""))
  message("DATABASE ", database_index, "/", n_databases, ": ", database_display_name, " [", database, "]")
  message(paste(rep("=", 88), collapse = ""))

  if (!test_mode && isTRUE(skip_completed_databases) && database_outputs_complete(paths)) {
    message("Complete output already exists; skipping this database.")
    return(list(status = "SKIPPED_COMPLETE", database = database, display_name = database_display_name, output_dir = paths$library_output_dir, executed_lists = 0L))
  }

  legacy_unnumbered_xlsx_files <- file.path(paths$xlsx_output_dir, paste0(paths$output_prefix, "_", analysis_variants$file_suffix))
  legacy_unnumbered_xlsx_files <- legacy_unnumbered_xlsx_files[file.exists(legacy_unnumbered_xlsx_files)]
  if (!test_mode && length(legacy_unnumbered_xlsx_files) > 0L) unlink(legacy_unnumbered_xlsx_files)

  analysis_results <- setNames(vector("list", nrow(analysis_variants)), analysis_variants$variant_id)

  for (variant_index in seq_len(nrow(analysis_variants))) {
    variant <- analysis_variants[variant_index, , drop = FALSE]
    variant_id <- variant$variant_id[[1]]
    direction <- variant$direction[[1]]
    analysis_results[[variant_id]] <- list()

    analysis_results[[variant_id]][["allGenes"]] <- empty_enrichment_table(
      analysis_status = "NOT_RUN_TEST_MODE", direction = direction, cluster_id = "ALL", cluster_name = "all clusters",
      n_query_genes = length(query_by_variant[[variant_id]][["allGenes"]])
    )

    for (cluster_id_current in expected_cluster_ids) {
      analysis_results[[variant_id]][[cluster_id_current]] <- empty_enrichment_table(
        analysis_status = "NOT_RUN_TEST_MODE", direction = direction, cluster_id = cluster_id_current,
        cluster_name = unname(custom_cluster_labels[[cluster_id_current]]),
        n_query_genes = length(query_by_variant[[variant_id]][[cluster_id_current]])
      )
    }
  }

  if (test_mode) {
    message("TEST MODE: executing first ", progress_total_lists, "/", planned_total_lists, " lists for this database.")
  } else {
    message("FULL MODE: executing all ", planned_total_lists, " lists for this database.")
  }

  database_progress <- create_database_progress(database = database_display_name, total_lists = progress_total_lists)
  draw_database_progress(database_progress, current_label = "starting")

  executed_variant_ids <- character(0)
  executed_list_count <- 0L
  stop_running <- FALSE

  for (variant_index in seq_len(nrow(analysis_variants))) {
    if (stop_running) break
    variant <- analysis_variants[variant_index, , drop = FALSE]
    variant_id <- variant$variant_id[[1]]
    subset_id <- variant$subset_id[[1]]
    direction <- variant$direction[[1]]

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

      analysis_results[[variant_id]][[scope_key]] <- run_one_enrichment(
        query_genes = query_by_variant[[variant_id]][[scope_key]],
        direction = direction,
        cluster_id = cluster_id_current,
        cluster_name = cluster_name_current,
        database = database
      )

      executed_list_count <- executed_list_count + 1L
      executed_variant_ids <- unique(c(executed_variant_ids, variant_id))
      advance_database_progress(database_progress, paste0(subset_id, " | ", direction, " | ", scope_label))
    }
  }

  variants_to_write <- if (test_mode) {
    analysis_variants$variant_id[analysis_variants$variant_id %in% executed_variant_ids]
  } else {
    analysis_variants$variant_id
  }

  message("\nWriting ", length(variants_to_write), if (test_mode) " test XLSX workbook(s)..." else " XLSX workbooks...")
  for (variant_id in variants_to_write) {
    variant_index <- match(variant_id, analysis_variants$variant_id)
    variant <- analysis_variants[variant_index, , drop = FALSE]
    write_direction_workbook(
      direction = variant$direction[[1]],
      result_list = analysis_results[[variant_id]],
      output_file = paths$xlsx_files[[variant_id]]
    )
  }
  message("XLSX workbook writing completed.")

  if (test_mode) {
    message("TEST MODE: PDFs skipped for ", database_display_name, ".")
    expected_output_files <- unname(paths$xlsx_files[variants_to_write])
  } else {
    message("Creating global PDF...")
    global_plots <- lapply(analysis_directions, function(direction) {
      make_global_barplot(
        enrichment_table = analysis_results[[unfiltered_variant_id(direction)]][["allGenes"]],
        direction = direction
      )
    })
    names(global_plots) <- analysis_directions

    global_combined_plot <- global_plots[["UPplusDOWN"]] / global_plots[["UP"]] / global_plots[["DOWN"]] +
      patchwork::plot_annotation(
        title = paste0("FMT donor-group main effect | ", database_display_name),
        subtitle = paste0(
          "DE: FDR < ", de_fdr_threshold, " & |log2FC| >= ", de_abs_log2fc_threshold,
          " | Significant-term count: FDR < ", enrichment_fdr_threshold, " & overlap >= ", min_overlap_genes,
          " genes | Bars: Top ", global_top_n, " by nominal P-value, without an FDR filter"
        ),
        caption = "Standard Enrichr analysis without a custom background; query = unique edgeR gene symbols submitted directly to Enrichr.",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 17),
          plot.subtitle = ggplot2::element_text(size = 11.5),
          plot.caption = ggplot2::element_text(size = 9.5, hjust = 0)
        )
      )

    grDevices::pdf(file = paths$pdf_files[["global"]], width = 16, height = 21, onefile = TRUE, useDingbats = FALSE)
    print(global_combined_plot)
    grDevices::dev.off()

    for (direction in analysis_directions) {
      message("Creating cluster-grid PDF: ", direction)
      cluster_plots <- lapply(expected_cluster_ids, function(cluster_id_current) {
        make_cluster_barplot(
          enrichment_table = analysis_results[[unfiltered_variant_id(direction)]][[cluster_id_current]],
          cluster_id = cluster_id_current,
          direction = direction
        )
      })

      direction_label <- switch(direction, UPplusDOWN = "UP + DOWN genes", UP = "UP genes", DOWN = "DOWN genes", direction)

      cluster_grid <- patchwork::wrap_plots(cluster_plots, ncol = 4, nrow = 4, byrow = TRUE) +
        patchwork::plot_annotation(
          title = paste0("FMT donor-group main effect | ", direction_label, " | ", database_display_name),
          subtitle = paste0(
            "DE: FDR < ", de_fdr_threshold, " & |log2FC| >= ", de_abs_log2fc_threshold,
            " | Significant-term count: FDR < ", enrichment_fdr_threshold, " & overlap >= ", min_overlap_genes,
            " genes | Bars: Top ", cluster_top_n, " by nominal P-value without an FDR filter | Bars = -log10(P-value)"
          ),
          caption = "Standard Enrichr analysis without a custom background. Top 10 bars are selected only by nominal P-value; FDR is used only for the significant-term count. X-axis scales are independent between panels; gene counts are submitted query genes.",
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 18),
            plot.subtitle = ggplot2::element_text(size = 11),
            plot.caption = ggplot2::element_text(size = 9.5, hjust = 0)
          )
        )

      grDevices::pdf(file = paths$pdf_files[[direction]], width = 36, height = 28, onefile = TRUE, useDingbats = FALSE)
      print(cluster_grid)
      grDevices::dev.off()
    }

    expected_output_files <- c(unname(paths$xlsx_files), unname(paths$pdf_files))
  }

  missing_output_files <- expected_output_files[!file.exists(expected_output_files)]
  if (length(missing_output_files) > 0L) {
    stop("Missing output file(s) for ", database, ":\n", paste(missing_output_files, collapse = "\n"), call. = FALSE)
  }

  empty_output_files <- expected_output_files[file.info(expected_output_files)$size == 0]
  if (length(empty_output_files) > 0L) {
    stop("Empty output file(s) for ", database, ":\n", paste(empty_output_files, collapse = "\n"), call. = FALSE)
  }

  message("Database completed successfully: ", database_display_name)
  message("Output directory: ", normalizePath(paths$library_output_dir, mustWork = TRUE))

  list(
    status = if (test_mode) "TEST_COMPLETED" else "COMPLETED",
    database = database,
    display_name = database_display_name,
    output_dir = paths$library_output_dir,
    executed_lists = executed_list_count
  )
}

# ==============================================================================
# 17. Run all selected databases
# ==============================================================================

message("\n", paste(rep("=", 88), collapse = ""))
message("MULTI-DATABASE ENRICHR RUN")
message(paste(rep("=", 88), collapse = ""))
message("Databases selected: ", nrow(enrichr_libraries))
message("Planned lists per database: ", progress_total_lists, if (test_mode) paste0(" (TEST subset of ", planned_total_lists, ")") else "")
message("Maximum list analyses in this run: ", progress_total_lists * nrow(enrichr_libraries))
message("One independent progress bar will be shown for each database.")

database_run_summary <- vector("list", nrow(enrichr_libraries))
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
      message("\nERROR in database ", database_display_name, " [", database, "]:")
      message(conditionMessage(e))
      if (!isTRUE(continue_after_database_error)) stop(e)
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
  invisible(gc())
}

multi_database_elapsed <- as.numeric(difftime(Sys.time(), multi_database_start_time, units = "secs"))
database_run_summary_df <- dplyr::bind_rows(lapply(database_run_summary, tibble::as_tibble))

# ==============================================================================
# 18. Final report
# ==============================================================================

message("\n", paste(rep("=", 88), collapse = ""))
message("MULTI-DATABASE FUNCTIONAL ANALYSIS FINISHED")
message(paste(rep("=", 88), collapse = ""))
message("Total elapsed time: ", format_elapsed_time(multi_database_elapsed))
message("\nDatabase status:")
for (i in seq_len(nrow(database_run_summary_df))) {
  message("  ", i, "/", nrow(database_run_summary_df), " | ", database_run_summary_df$display_name[[i]], " | ", database_run_summary_df$status[[i]])
}

message("\nCommon DE selection: FDR < ", de_fdr_threshold, " & |log2FC| >= ", de_abs_log2fc_threshold)
message("Barplots: Top ", global_top_n, " terms by smallest nominal Enrichr P-value; no FDR filter")
message("Significant-term count: FDR < ", enrichment_fdr_threshold, " & overlap >= ", min_overlap_genes, " genes")
message("XLSX highlight: FDR < ", xlsx_highlight_fdr_threshold, " & overlap >= ", xlsx_highlight_min_overlap_genes, " genes | #DEDCE6")
message("Comparison: ", selected_test_comparison)
message("Test ID: ", donor_group_test_id)

if (length(database_failures) > 0L) {
  message("\nFAILED DATABASES:")
  for (database in names(database_failures)) message("  ", database, " -> ", database_failures[[database]])
  stop(length(database_failures), " Enrichr database(s) failed. Other databases were allowed to finish; see summary above.", call. = FALSE)
}

message("\nAll requested databases completed or were already complete.")

# ==============================================================================
# End
# ==============================================================================
