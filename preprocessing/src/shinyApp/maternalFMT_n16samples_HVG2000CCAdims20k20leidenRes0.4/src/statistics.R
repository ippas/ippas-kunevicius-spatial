# ==============================================================================
# Selected-gene edgeR statistics
# ==============================================================================

prepare_test_definition_table <- function(
  test_definitions
) {
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

  test_definitions |>
    dplyr::mutate(
      test_id_chr = as.character(.data$test_id),
      sheet_name_chr = as.character(.data$sheet_name),
      comparison_chr = as.character(.data$comparison),
      test_id_lower = tolower(.data$test_id_chr),
      sheet_name_lower = tolower(.data$sheet_name_chr),
      comparison_lower = tolower(.data$comparison_chr)
    )
}

select_unique_best_test <- function(
  scored_rows,
  score_column,
  effect_label,
  test_definitions
) {
  if (nrow(scored_rows) == 0L) {
    stop(
      "Could not identify any edgeR candidate for ",
      effect_label,
      ". Available test_id values: ",
      paste(test_definitions$test_id, collapse = ", "),
      call. = FALSE
    )
  }

  score_values <- scored_rows[[score_column]]
  best_score <- max(score_values)

  best_rows <- scored_rows[
    score_values == best_score,
    ,
    drop = FALSE
  ]

  if (nrow(best_rows) != 1L) {
    stop(
      "Could not uniquely identify the edgeR test for ",
      effect_label,
      ". Best candidates: ",
      paste(best_rows$test_id_chr, collapse = ", "),
      ". Available test_id values: ",
      paste(test_definitions$test_id, collapse = ", "),
      call. = FALSE
    )
  }

  as.character(
    best_rows$test_id_chr[[1]]
  )
}

find_donor_group_test_id <- function(
  test_definitions
) {
  test_table <- prepare_test_definition_table(
    test_definitions
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
      donor_score =
        5L * (
          grepl("donor", .data$test_id_lower) |
            grepl("fmt", .data$test_id_lower)
        ) +
        4L * (
          grepl("donor", .data$sheet_name_lower) |
            grepl("fmt", .data$sheet_name_lower)
        ) +
        3L * grepl("group", .data$test_id_lower) +
        2L * grepl("group", .data$sheet_name_lower) +
        3L * (
          grepl("asd", .data$comparison_lower) &
            grepl("neurotypical", .data$comparison_lower)
        ) +
        1L * (
          grepl("asd", .data$test_id_lower) |
            grepl("neurotypical", .data$test_id_lower)
        ) +
        5L * (
          grepl("overall", .data$test_id_lower) |
            grepl("overall", .data$sheet_name_lower)
        )
    ) |>
    dplyr::filter(
      .data$donor_score > 0L
    )

  select_unique_best_test(
    scored_rows = scored_rows,
    score_column = "donor_score",
    effect_label = "donor-group main effect",
    test_definitions = test_definitions
  )
}

find_sex_test_id <- function(
  test_definitions
) {
  test_table <- prepare_test_definition_table(
    test_definitions
  )

  eligible_rows <- test_table |>
    dplyr::filter(
      !grepl("interaction", .data$test_id_lower),
      !grepl("interaction", .data$sheet_name_lower),
      !grepl("simple", .data$test_id_lower),
      !grepl("simple", .data$sheet_name_lower)
    )

  scored_rows <- eligible_rows |>
    dplyr::mutate(
      sex_score =
        6L * grepl("sex", .data$test_id_lower) +
        5L * grepl("sex", .data$sheet_name_lower) +
        4L * (
          grepl("female", .data$comparison_lower) &
            grepl("male", .data$comparison_lower)
        ) +
        2L * (
          grepl("female", .data$test_id_lower) |
            grepl("male", .data$test_id_lower)
        ) +
        5L * (
          grepl("overall", .data$test_id_lower) |
            grepl("overall", .data$sheet_name_lower)
        )
    ) |>
    dplyr::filter(
      .data$sex_score > 0L
    )

  select_unique_best_test(
    scored_rows = scored_rows,
    score_column = "sex_score",
    effect_label = "sex main effect",
    test_definitions = test_definitions
  )
}

find_interaction_test_id <- function(
  test_definitions
) {
  test_table <- prepare_test_definition_table(
    test_definitions
  )

  eligible_rows <- test_table |>
    dplyr::filter(
      grepl("interaction", .data$test_id_lower) |
        grepl("interaction", .data$sheet_name_lower) |
        grepl("interaction", .data$comparison_lower)
    )

  scored_rows <- eligible_rows |>
    dplyr::mutate(
      interaction_score =
        8L * grepl("interaction", .data$test_id_lower) +
        7L * grepl("interaction", .data$sheet_name_lower) +
        4L * grepl("interaction", .data$comparison_lower) +
        3L * (
          grepl("group", .data$test_id_lower) |
            grepl("donor", .data$test_id_lower) |
            grepl("fmt", .data$test_id_lower)
        ) +
        3L * grepl("sex", .data$test_id_lower) +
        4L * (
          grepl("overall", .data$test_id_lower) |
            grepl("overall", .data$sheet_name_lower)
        )
    ) |>
    dplyr::filter(
      .data$interaction_score > 0L
    )

  select_unique_best_test(
    scored_rows = scored_rows,
    score_column = "interaction_score",
    effect_label = "Group × Sex interaction",
    test_definitions = test_definitions
  )
}

resolve_shiny_effect_test_ids <- function(
  test_definitions,
  edger_results
) {
  selected_ids <- c(
    group = find_donor_group_test_id(
      test_definitions
    ),
    sex = find_sex_test_id(
      test_definitions
    ),
    interaction = find_interaction_test_id(
      test_definitions
    )
  )

  missing_result_tables <- setdiff(
    unname(selected_ids),
    names(edger_results)
  )

  if (length(missing_result_tables) > 0L) {
    stop(
      "Selected edgeR test(s) are absent from edgeR_perCluster_combinedResults: ",
      paste(missing_result_tables, collapse = ", "),
      call. = FALSE
    )
  }

  selected_ids
}

get_edger_effect_table <- function(
  edger_results,
  effect_test_ids,
  effect
) {
  if (!effect %in% names(effect_test_ids)) {
    stop(
      "Unknown edgeR effect: ",
      effect,
      call. = FALSE
    )
  }

  result_id <- unname(
    effect_test_ids[[effect]]
  )

  result_table <- edger_results[[result_id]]

  required_columns <- c(
    "cluster_id",
    "gene",
    "logFC",
    "PValue",
    "FDR"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(result_table)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Missing edgeR result column(s) in ",
      result_id,
      ": ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  result_table$cluster_id <- as.character(
    result_table$cluster_id
  )

  list(
    result_id = result_id,
    table = result_table
  )
}

get_test_definition <- function(
  edger_test_definitions,
  result_id
) {
  definition <- edger_test_definitions |>
    dplyr::filter(
      as.character(.data$test_id) == result_id
    )

  if (nrow(definition) != 1L) {
    return(
      list(
        comparison = result_id,
        sheet_name = NA_character_
      )
    )
  }

  list(
    comparison = as.character(
      definition$comparison[[1]]
    ),
    sheet_name = as.character(
      definition$sheet_name[[1]]
    )
  )
}

get_gene_edger_statistics <- function(
  edger_results,
  edger_test_definitions,
  effect_test_ids,
  effect,
  gene_symbol,
  cluster_names = NULL
) {
  effect_data <- get_edger_effect_table(
    edger_results = edger_results,
    effect_test_ids = effect_test_ids,
    effect = effect
  )

  gene_results <- effect_data$table |>
    dplyr::filter(
      .data$gene == .env$gene_symbol
    ) |>
    dplyr::select(
      .data$cluster_id,
      .data$logFC,
      .data$PValue,
      .data$FDR
    )

  if (anyDuplicated(gene_results$cluster_id) > 0L) {
    duplicated_clusters <- unique(
      gene_results$cluster_id[duplicated(gene_results$cluster_id)]
    )
    stop(
      "Selected gene has more than one edgeR row in cluster(s): ",
      paste(paste0("C", duplicated_clusters), collapse = ", "),
      " | gene: ",
      gene_symbol,
      " | test_id: ",
      effect_data$result_id,
      call. = FALSE
    )
  }

  definition <- get_test_definition(
    edger_test_definitions = edger_test_definitions,
    result_id = effect_data$result_id
  )

  if (nrow(gene_results) == 0L) {
    return(
      list(
        result_id = effect_data$result_id,
        comparison = definition$comparison,
        raw = gene_results,
        display = data.frame()
      )
    )
  }

  cluster_order <- sort_cluster_ids(
    gene_results$cluster_id
  )

  gene_results <- gene_results[
    match(
      cluster_order,
      gene_results$cluster_id
    ),
    ,
    drop = FALSE
  ]

  cluster_labels <- paste0(
    "C",
    gene_results$cluster_id
  )

  anatomical_region <- rep(
    "",
    nrow(gene_results)
  )

  if (!is.null(cluster_names)) {
    has_name <- gene_results$cluster_id %in%
      names(cluster_names)

    anatomical_region[has_name] <- unname(
      cluster_names[
        gene_results$cluster_id[has_name]
      ]
    )
  }

  display_table <- data.frame(
    Cluster = cluster_labels,
    Region = anatomical_region,
    log2FC = formatC(
      gene_results$logFC,
      format = "f",
      digits = 3
    ),
    PValue = format(
      gene_results$PValue,
      scientific = TRUE,
      digits = 3
    ),
    FDR = format(
      gene_results$FDR,
      scientific = TRUE,
      digits = 3
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  list(
    result_id = effect_data$result_id,
    comparison = definition$comparison,
    raw = gene_results,
    display = display_table
  )
}

get_significant_clusters_for_gene <- function(
  edger_results,
  effect_test_ids,
  effect,
  gene_symbol,
  fdr_threshold = 0.05,
  abs_log2fc_threshold = 0.5
) {
  effect_data <- get_edger_effect_table(
    edger_results = edger_results,
    effect_test_ids = effect_test_ids,
    effect = effect
  )

  significant_rows <- effect_data$table |>
    dplyr::filter(
      .data$gene == .env$gene_symbol,
      is.finite(.data$FDR),
      is.finite(.data$logFC),
      .data$FDR < fdr_threshold,
      abs(.data$logFC) >= abs_log2fc_threshold
    )

  sort_cluster_ids(
    significant_rows$cluster_id
  )
}
