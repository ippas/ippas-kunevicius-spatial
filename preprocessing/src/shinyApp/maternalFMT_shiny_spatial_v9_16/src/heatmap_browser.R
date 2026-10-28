# ==============================================================================
# Manual gene-expression heatmap browser
# ============================================================================== 

manual_heatmap_metadata_test_id <- function(edger_results) {
  preferred_id <- "Overall_Group_ASD_vs_Neurotypical"

  if (preferred_id %in% names(edger_results)) {
    return(preferred_id)
  }

  required_columns <- c(
    "cluster_id",
    "ensembl_gene_id",
    "gene",
    "mean_percent_positive_spots_Neurotypical",
    "mean_percent_positive_spots_ASD"
  )

  candidates <- names(edger_results)[
    vapply(
      edger_results,
      function(x) {
        is.data.frame(x) && all(required_columns %in% colnames(x))
      },
      logical(1)
    )
  ]

  if (length(candidates) == 0L) {
    stop(
      "Could not find an edgeR result table containing donor-group percent-positive columns required for the heatmap annotation.",
      call. = FALSE
    )
  }

  candidates[[1]]
}

parse_manual_heatmap_genes <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(trimws(x))) {
    return(character(0))
  }

  genes <- unlist(
    strsplit(
      x,
      split = "[,;[:space:]]+",
      perl = TRUE
    ),
    use.names = FALSE
  )

  genes <- trimws(genes)
  genes <- genes[nzchar(genes)]

  if (length(genes) == 0L) {
    return(character(0))
  }

  genes[!duplicated(tolower(genes))]
}

manual_heatmap_max_mean_percent <- function(neurotypical, asd) {
  neurotypical <- suppressWarnings(as.numeric(neurotypical))
  asd <- suppressWarnings(as.numeric(asd))

  output <- pmax(neurotypical, asd, na.rm = TRUE)
  output[is.infinite(output)] <- NA_real_
  output
}

manual_heatmap_metadata_table <- function(
  edger_results,
  selected_clusters,
  per_cluster_expression
) {
  selected_clusters <- sort_cluster_ids(
    intersect(
      as.character(selected_clusters),
      names(per_cluster_expression$gene_metadata)
    )
  )

  if (length(selected_clusters) == 0L) {
    return(
      tibble::tibble(
        cluster_id = character(0),
        ensembl_gene_id = character(0),
        gene = character(0),
        chromosome = character(0),
        max_mean_percent = numeric(0)
      )
    )
  }

  expression_metadata <- dplyr::bind_rows(
    lapply(
      selected_clusters,
      function(cluster_id) {
        current <- per_cluster_expression$gene_metadata[[cluster_id]]

        tibble::tibble(
          cluster_id = cluster_id,
          ensembl_gene_id = as.character(current$ensembl_gene_id),
          gene = as.character(current$gene),
          chromosome = as.character(current$chromosome)
        )
      }
    )
  ) |>
    dplyr::mutate(
      gene = dplyr::if_else(
        is.na(.data$gene) | .data$gene == "",
        .data$ensembl_gene_id,
        .data$gene
      )
    )

  test_id <- manual_heatmap_metadata_test_id(edger_results)
  percent_table <- tibble::as_tibble(edger_results[[test_id]])

  required_percent_columns <- c(
    "cluster_id",
    "ensembl_gene_id",
    "mean_percent_positive_spots_Neurotypical",
    "mean_percent_positive_spots_ASD"
  )

  missing_percent_columns <- setdiff(
    required_percent_columns,
    colnames(percent_table)
  )

  if (length(missing_percent_columns) > 0L) {
    stop(
      "Heatmap percent-positive metadata is missing required column(s): ",
      paste(missing_percent_columns, collapse = ", "),
      call. = FALSE
    )
  }

  percent_table <- percent_table |>
    dplyr::transmute(
      cluster_id = as.character(.data$cluster_id),
      ensembl_gene_id = as.character(.data$ensembl_gene_id),
      max_mean_percent = manual_heatmap_max_mean_percent(
        .data$mean_percent_positive_spots_Neurotypical,
        .data$mean_percent_positive_spots_ASD
      )
    ) |>
    dplyr::filter(.data$cluster_id %in% selected_clusters) |>
    dplyr::distinct(
      .data$cluster_id,
      .data$ensembl_gene_id,
      .keep_all = TRUE
    )

  expression_metadata |>
    dplyr::left_join(
      percent_table,
      by = c("cluster_id", "ensembl_gene_id")
    ) |>
    dplyr::distinct()
}

manual_heatmap_gene_availability <- function(
  gene_text,
  selected_clusters,
  edger_results,
  per_cluster_expression
) {
  requested <- parse_manual_heatmap_genes(gene_text)

  if (length(requested) == 0L || length(selected_clusters) == 0L) {
    return(
      list(
        requested = requested,
        found = character(0),
        missing = requested
      )
    )
  }

  metadata <- manual_heatmap_metadata_table(
    edger_results = edger_results,
    selected_clusters = selected_clusters,
    per_cluster_expression = per_cluster_expression
  )

  gene_lower <- tolower(metadata$gene)
  ensembl_lower <- tolower(metadata$ensembl_gene_id)

  found <- requested[
    vapply(
      requested,
      function(requested_gene) {
        token <- tolower(requested_gene)
        any(gene_lower == token | ensembl_lower == token)
      },
      logical(1)
    )
  ]

  list(
    requested = requested,
    found = found,
    missing = requested[!tolower(requested) %in% tolower(found)]
  )
}

manual_heatmap_sample_metadata <- function(
  edger_model_objects,
  selected_clusters
) {
  selected_clusters <- intersect(
    as.character(selected_clusters),
    names(edger_model_objects)
  )

  metadata_list <- lapply(
    selected_clusters,
    function(cluster_id) {
      model_object <- edger_model_objects[[cluster_id]]

      if (is.null(model_object$dge) || is.null(model_object$dge$samples)) {
        return(NULL)
      }

      samples <- model_object$dge$samples |>
        tibble::rownames_to_column("dge_rowname") |>
        tibble::as_tibble()

      required_columns <- c(
        "sample_ID",
        "fmt_donor_group",
        "sex"
      )

      if (!all(required_columns %in% colnames(samples))) {
        stop(
          "DGEList sample metadata for cluster ",
          cluster_id,
          " is missing: ",
          paste(setdiff(required_columns, colnames(samples)), collapse = ", "),
          call. = FALSE
        )
      }

      samples |>
        dplyr::transmute(
          sample_ID = as.character(.data$sample_ID),
          fmt_donor_group = as.character(.data$fmt_donor_group),
          sex = as.character(.data$sex)
        )
    }
  )

  metadata <- dplyr::bind_rows(metadata_list) |>
    dplyr::distinct(
      .data$sample_ID,
      .data$fmt_donor_group,
      .data$sex
    )

  if (nrow(metadata) == 0L) {
    stop(
      "No sample metadata could be recovered from the selected cluster DGEList objects.",
      call. = FALSE
    )
  }

  conflicts <- metadata |>
    dplyr::count(.data$sample_ID, name = "n") |>
    dplyr::filter(.data$n > 1L)

  if (nrow(conflicts) > 0L) {
    stop(
      "Conflicting group/sex metadata were found for sample(s): ",
      paste(conflicts$sample_ID, collapse = ", "),
      call. = FALSE
    )
  }

  metadata |>
    dplyr::mutate(
      sex = factor(
        .data$sex,
        levels = c("Male", "Female")
      ),
      fmt_donor_group = factor(
        .data$fmt_donor_group,
        levels = c("Neurotypical", "ASD")
      )
    ) |>
    dplyr::arrange(
      .data$sex,
      .data$fmt_donor_group,
      .data$sample_ID
    ) |>
    dplyr::mutate(
      combined_group = factor(
        paste(
          as.character(.data$sex),
          as.character(.data$fmt_donor_group),
          sep = "_"
        ),
        levels = c(
          "Male_Neurotypical",
          "Male_ASD",
          "Female_Neurotypical",
          "Female_ASD"
        )
      )
    )
}

calculate_manual_heatmap_row_zscores <- function(
  expression_matrix,
  lower_limit = -2,
  upper_limit = 2
) {
  output <- matrix(
    NA_real_,
    nrow = nrow(expression_matrix),
    ncol = ncol(expression_matrix),
    dimnames = dimnames(expression_matrix)
  )

  for (row_index in seq_len(nrow(expression_matrix))) {
    values <- as.numeric(expression_matrix[row_index, ])
    finite_values <- is.finite(values)

    if (sum(finite_values) < 2L) {
      output[row_index, finite_values] <- 0
      next
    }

    row_sd <- stats::sd(values[finite_values])

    if (!is.finite(row_sd) || row_sd == 0) {
      output[row_index, finite_values] <- 0
      next
    }

    output[row_index, finite_values] <-
      (values[finite_values] - mean(values[finite_values])) / row_sd
  }

  output[is.finite(output) & output < lower_limit] <- lower_limit
  output[is.finite(output) & output > upper_limit] <- upper_limit
  output
}

manual_heatmap_modified_distance <- function(
  expression_matrix,
  cluster_labels,
  same_cluster_bonus = 0
) {
  if (!is.matrix(expression_matrix)) {
    expression_matrix <- as.matrix(expression_matrix)
  }

  if (nrow(expression_matrix) <= 1L) {
    return(stats::as.dist(matrix(0, nrow = 1L, ncol = 1L)))
  }

  if (
    !is.numeric(same_cluster_bonus) ||
      length(same_cluster_bonus) != 1L ||
      !is.finite(same_cluster_bonus) ||
      same_cluster_bonus < 0 ||
      same_cluster_bonus >= 1
  ) {
    stop(
      "Same-cluster clustering bonus must be a finite number in [0, 1).",
      call. = FALSE
    )
  }

  correlation_matrix <- suppressWarnings(
    stats::cor(
      t(expression_matrix),
      method = "pearson",
      use = "pairwise.complete.obs"
    )
  )

  correlation_matrix[!is.finite(correlation_matrix)] <- 0
  diag(correlation_matrix) <- 1

  base_distance <- 2 - correlation_matrix

  same_cluster <- outer(
    as.character(cluster_labels),
    as.character(cluster_labels),
    FUN = "=="
  )

  distance_matrix <- base_distance
  distance_matrix[same_cluster] <-
    base_distance[same_cluster] * (1 - same_cluster_bonus)
  distance_matrix[!same_cluster] <-
    base_distance[!same_cluster] * (1 + same_cluster_bonus)

  diag(distance_matrix) <- 0
  distance_matrix[distance_matrix < 0] <- 0

  stats::as.dist(distance_matrix)
}

manual_heatmap_row_hclust <- function(
  expression_matrix,
  cluster_labels,
  same_cluster_bonus = 0
) {
  if (nrow(expression_matrix) <= 1L) {
    return(FALSE)
  }

  stats::hclust(
    manual_heatmap_modified_distance(
      expression_matrix = expression_matrix,
      cluster_labels = cluster_labels,
      same_cluster_bonus = same_cluster_bonus
    ),
    method = "complete"
  )
}

prepare_manual_gene_heatmap <- function(
  gene_text,
  selected_clusters,
  same_cluster_bonus,
  zscore_limit,
  edger_results,
  edger_model_objects,
  per_cluster_expression
) {
  requested_genes <- parse_manual_heatmap_genes(gene_text)

  if (length(requested_genes) == 0L) {
    return(
      list(
        valid = FALSE,
        message = "Enter at least one gene symbol.",
        requested = character(0),
        found = character(0),
        missing = character(0)
      )
    )
  }

  zscore_limit <- suppressWarnings(as.numeric(zscore_limit))
  if (
    length(zscore_limit) != 1L ||
      !is.finite(zscore_limit) ||
      zscore_limit <= 0
  ) {
    return(
      list(
        valid = FALSE,
        message = "Row Z-score limit must be one positive number.",
        requested = requested_genes,
        found = character(0),
        missing = requested_genes
      )
    )
  }

  available_expression_clusters <- intersect(
    names(per_cluster_expression$raw_counts),
    names(per_cluster_expression$log2_qn)
  )

  selected_clusters <- sort_cluster_ids(
    intersect(
      as.character(selected_clusters),
      available_expression_clusters
    )
  )

  if (length(selected_clusters) == 0L) {
    return(
      list(
        valid = FALSE,
        message = "Select at least one cluster.",
        requested = requested_genes,
        found = character(0),
        missing = requested_genes
      )
    )
  }

  metadata <- manual_heatmap_metadata_table(
    edger_results = edger_results,
    selected_clusters = selected_clusters,
    per_cluster_expression = per_cluster_expression
  )

  requested_rows <- lapply(
    seq_along(requested_genes),
    function(request_index) {
      requested_gene <- requested_genes[[request_index]]
      token <- tolower(requested_gene)

      hits <- metadata[
        tolower(metadata$gene) == token |
          tolower(metadata$ensembl_gene_id) == token,
        ,
        drop = FALSE
      ]

      if (nrow(hits) == 0L) {
        return(NULL)
      }

      hits |>
        dplyr::mutate(
          requested_gene = requested_gene,
          requested_order = request_index
        )
    }
  )

  row_table <- dplyr::bind_rows(requested_rows)

  if (nrow(row_table) == 0L) {
    return(
      list(
        valid = FALSE,
        message = "None of the requested genes were found in the selected per-cluster expression matrices.",
        requested = requested_genes,
        found = character(0),
        missing = requested_genes
      )
    )
  }

  row_table <- row_table |>
    dplyr::mutate(
      cluster_order = match(.data$cluster_id, selected_clusters)
    ) |>
    dplyr::arrange(
      .data$requested_order,
      .data$cluster_order,
      .data$ensembl_gene_id
    ) |>
    dplyr::distinct(
      .data$requested_gene,
      .data$cluster_id,
      .data$ensembl_gene_id,
      .keep_all = TRUE
    )

  found_genes <- requested_genes[
    tolower(requested_genes) %in%
      unique(tolower(row_table$requested_gene))
  ]
  missing_genes <- requested_genes[
    !tolower(requested_genes) %in% tolower(found_genes)
  ]

  sample_metadata <- manual_heatmap_sample_metadata(
    edger_model_objects = edger_model_objects,
    selected_clusters = selected_clusters
  )
  sample_order <- as.character(sample_metadata$sample_ID)

  row_labels <- paste0(
    row_table$gene,
    " [C",
    row_table$cluster_id,
    "]"
  )
  row_labels <- make.unique(row_labels, sep = " #")
  row_table$row_id <- row_labels

  raw_counts_matrix <- matrix(
    NA_real_,
    nrow = nrow(row_table),
    ncol = length(sample_order),
    dimnames = list(row_labels, sample_order)
  )

  log2_qn_matrix <- matrix(
    NA_real_,
    nrow = nrow(row_table),
    ncol = length(sample_order),
    dimnames = list(row_labels, sample_order)
  )

  for (cluster_id in unique(row_table$cluster_id)) {
    row_indexes <- which(row_table$cluster_id == cluster_id)

    raw_current <- per_cluster_expression$raw_counts[[cluster_id]]
    qn_current <- per_cluster_expression$log2_qn[[cluster_id]]

    if (is.null(raw_current) || is.null(qn_current)) {
      stop(
        "Precomputed raw/log2+1 QN matrix is missing for cluster C",
        cluster_id,
        ".",
        call. = FALSE
      )
    }

    gene_indexes_raw <- match(
      row_table$ensembl_gene_id[row_indexes],
      rownames(raw_current)
    )
    gene_indexes_qn <- match(
      row_table$ensembl_gene_id[row_indexes],
      rownames(qn_current)
    )

    if (anyNA(gene_indexes_raw) || anyNA(gene_indexes_qn)) {
      stop(
        "Requested heatmap gene is missing from the precomputed expression matrix for cluster C",
        cluster_id,
        ".",
        call. = FALSE
      )
    }

    available_samples <- Reduce(
      intersect,
      list(
        sample_order,
        colnames(raw_current),
        colnames(qn_current)
      )
    )

    if (length(available_samples) == 0L) {
      stop(
        "No sample columns overlap for cluster C",
        cluster_id,
        ".",
        call. = FALSE
      )
    }

    raw_counts_matrix[row_indexes, available_samples] <- raw_current[
      gene_indexes_raw,
      available_samples,
      drop = FALSE
    ]

    log2_qn_matrix[row_indexes, available_samples] <- qn_current[
      gene_indexes_qn,
      available_samples,
      drop = FALSE
    ]
  }

  zscore_matrix <- calculate_manual_heatmap_row_zscores(
    expression_matrix = log2_qn_matrix,
    lower_limit = -zscore_limit,
    upper_limit = zscore_limit
  )

  same_cluster_bonus <- suppressWarnings(
    as.numeric(same_cluster_bonus)
  )
  if (
    length(same_cluster_bonus) != 1L ||
      !is.finite(same_cluster_bonus) ||
      same_cluster_bonus < 0 ||
      same_cluster_bonus >= 1
  ) {
    return(
      list(
        valid = FALSE,
        message = "Same-cluster bonus must be between 0 and less than 1.",
        requested = requested_genes,
        found = found_genes,
        missing = missing_genes
      )
    )
  }

  row_hclust <- manual_heatmap_row_hclust(
    expression_matrix = zscore_matrix,
    cluster_labels = row_table$cluster_id,
    same_cluster_bonus = same_cluster_bonus
  )

  row_height_mm <- 4.2
  heatmap_body_height_mm <- nrow(zscore_matrix) * row_height_mm
  heatmap_body_width_mm <- ncol(zscore_matrix) * 6.2

  browser_height_px <- max(
    900,
    round(
      (heatmap_body_height_mm + 125) * 3.78
    )
  )

  pdf_height_inches <- max(
    8,
    (heatmap_body_height_mm + 90) / 25.4
  )

  list(
    valid = TRUE,
    message = NULL,
    requested = requested_genes,
    found = found_genes,
    missing = missing_genes,
    row_table = row_table,
    sample_metadata = sample_metadata,
    raw_counts_matrix = raw_counts_matrix,
    log2_qn_matrix = log2_qn_matrix,
    zscore_matrix = zscore_matrix,
    row_hclust = row_hclust,
    row_height_mm = row_height_mm,
    heatmap_body_height_mm = heatmap_body_height_mm,
    heatmap_body_width_mm = heatmap_body_width_mm,
    browser_height_px = browser_height_px,
    pdf_width_inches = 14,
    pdf_height_inches = pdf_height_inches,
    same_cluster_bonus = same_cluster_bonus,
    zscore_limit = zscore_limit
  )
}

manual_heatmap_percent_annotation <- function(
  values,
  font_size,
  colour_function
) {
  labels <- ifelse(
    is.finite(values),
    paste0(round(values), "%"),
    "NA"
  )

  text_colours <- ifelse(
    is.finite(values) & values > 70,
    "white",
    "black"
  )

  ComplexHeatmap::AnnotationFunction(
    which = "row",
    width = grid::unit(6.0, "mm"),
    var_import = list(
      values = values,
      labels = labels,
      text_colours = text_colours,
      colour_function = colour_function,
      font_size = font_size
    ),
    fun = function(index) {
      n_current <- length(index)
      y_positions <- 1 - (seq_len(n_current) - 0.5) / n_current
      current_values <- values[index]
      current_labels <- labels[index]
      current_text_colours <- text_colours[index]
      current_fill <- rep("grey90", n_current)
      finite_values <- is.finite(current_values)

      current_fill[finite_values] <- colour_function(
        current_values[finite_values]
      )

      grid::grid.rect(
        x = grid::unit(0.5, "npc"),
        y = grid::unit(y_positions, "npc"),
        width = grid::unit(1, "npc"),
        height = grid::unit(1 / n_current, "npc"),
        gp = grid::gpar(fill = current_fill, col = NA)
      )

      grid::grid.text(
        label = current_labels,
        x = grid::unit(0.5, "npc"),
        y = grid::unit(y_positions, "npc"),
        gp = grid::gpar(
          col = current_text_colours,
          fontsize = font_size,
          fontface = "bold"
        )
      )
    }
  )
}


manual_heatmap_logfc_annotation <- function(
  values,
  font_size
) {
  labels <- vapply(
    values,
    function(value) {
      if (!is.finite(value)) {
        return("NA")
      }

      label <- formatC(
        round(value, 2),
        format = "f",
        digits = 2
      )

      sub("\\.?0+$", "", label)
    },
    character(1)
  )

  finite_values <- values[is.finite(values)]
  colour_limit <- if (length(finite_values) == 0L) {
    1
  } else {
    max(abs(finite_values))
  }

  if (!is.finite(colour_limit) || colour_limit <= 0) {
    colour_limit <- 1
  }

  colour_function <- circlize::colorRamp2(
    c(-colour_limit, 0, colour_limit),
    c("#2166AC", "white", "#B2182B")
  )

  text_colours <- ifelse(
    is.finite(values) & abs(values) >= colour_limit * 0.65,
    "white",
    "black"
  )

  ComplexHeatmap::AnnotationFunction(
    which = "row",
    width = grid::unit(7.0, "mm"),
    var_import = list(
      values = values,
      labels = labels,
      text_colours = text_colours,
      colour_function = colour_function,
      font_size = font_size
    ),
    fun = function(index) {
      n_current <- length(index)
      y_positions <- 1 - (seq_len(n_current) - 0.5) / n_current
      current_values <- values[index]
      current_labels <- labels[index]
      current_text_colours <- text_colours[index]
      current_fill <- rep("grey90", n_current)
      finite_current <- is.finite(current_values)

      current_fill[finite_current] <- colour_function(
        current_values[finite_current]
      )

      grid::grid.rect(
        x = grid::unit(0.5, "npc"),
        y = grid::unit(y_positions, "npc"),
        width = grid::unit(1, "npc"),
        height = grid::unit(1 / n_current, "npc"),
        gp = grid::gpar(fill = current_fill, col = NA)
      )

      grid::grid.text(
        label = current_labels,
        x = grid::unit(0.5, "npc"),
        y = grid::unit(y_positions, "npc"),
        gp = grid::gpar(
          col = current_text_colours,
          fontsize = font_size,
          fontface = "bold"
        )
      )
    }
  )
}

build_manual_gene_heatmap <- function(
  heatmap_state,
  cluster_names,
  cluster_colors
) {
  if (!isTRUE(heatmap_state$valid)) {
    stop("Cannot build an invalid heatmap state.", call. = FALSE)
  }

  row_table <- heatmap_state$row_table
  sample_metadata <- heatmap_state$sample_metadata
  zscore_matrix <- heatmap_state$zscore_matrix

  displayed_clusters <- sort_cluster_ids(
    unique(row_table$cluster_id)
  )

  cluster_annotation_colors <- cluster_colors[displayed_clusters]
  names(cluster_annotation_colors) <- displayed_clusters

  cluster_legend_labels <- vapply(
    displayed_clusters,
    function(cluster_id) {
      region <- cluster_names[[cluster_id]]
      if (is.null(region) || is.na(region) || !nzchar(region)) {
        paste0("C", cluster_id)
      } else {
        paste0("C", cluster_id, " — ", region)
      }
    },
    character(1)
  )

  zscore_limit <- heatmap_state$zscore_limit

  heatmap_colour_function <- circlize::colorRamp2(
    c(-zscore_limit, 0, zscore_limit),
    c("#2166AC", "white", "#B2182B")
  )

  percent_colour_function <- circlize::colorRamp2(
    c(0, 25, 50, 75, 100),
    c(
      "#F3E8FF",
      "#DDD6FE",
      "#C084FC",
      "#9333EA",
      "#4C1D95"
    )
  )

  group_colors <- c(
    "Neurotypical" = "grey70",
    "ASD" = "black"
  )

  sex_colors <- c(
    "Male" = "#3C78D8",
    "Female" = "#C050A0"
  )

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Group = sample_metadata$fmt_donor_group,
    Sex = sample_metadata$sex,
    col = list(
      Group = group_colors,
      Sex = sex_colors
    ),
    simple_anno_size = grid::unit(4.2, "mm"),
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    annotation_legend_param = list(
      Group = list(title = "FMT donor group", ncol = 1),
      Sex = list(title = "Sex", ncol = 1)
    )
  )

  annotation_font_size <- 7.5

  gene_annotation_width <- ComplexHeatmap::max_text_width(
    row_table$row_id,
    gp = grid::gpar(fontsize = annotation_font_size)
  ) + grid::unit(1.5, "mm")

  has_significance <-
    !is.null(heatmap_state$significance_test_id) &&
      nzchar(heatmap_state$significance_test_id)

  if (has_significance) {
    right_annotation <- ComplexHeatmap::rowAnnotation(
      Cluster = row_table$cluster_id,
      `Max mean %` = manual_heatmap_percent_annotation(
        values = row_table$max_mean_percent,
        font_size = annotation_font_size,
        colour_function = percent_colour_function
      ),
      log2FC = manual_heatmap_logfc_annotation(
        values = row_table$significance_logfc,
        font_size = annotation_font_size
      ),
      Significance = ComplexHeatmap::anno_text(
        row_table$significance_symbol,
        location = grid::unit(0.5, "npc"),
        just = "center",
        gp = grid::gpar(
          col = "black",
          fontsize = heatmap_state$significance_star_size,
          fontface = "bold"
        ),
        width = grid::unit(7.0, "mm")
      ),
      Gene = ComplexHeatmap::anno_text(
        row_table$row_id,
        location = grid::unit(0, "npc"),
        just = "left",
        gp = grid::gpar(
          col = "black",
          fontsize = annotation_font_size
        ),
        width = gene_annotation_width
      ),
      col = list(
        Cluster = cluster_annotation_colors
      ),
      simple_anno_size = grid::unit(4.2, "mm"),
      gap = grid::unit(0.6, "mm"),
      show_annotation_name = c(
        Cluster = TRUE,
        `Max mean %` = TRUE,
        log2FC = TRUE,
        Significance = TRUE,
        Gene = TRUE
      ),
      annotation_label = c(
        Cluster = "Cluster",
        `Max mean %` = "Max mean %",
        log2FC = "log2FC",
        Significance = heatmap_state$significance_metric,
        Gene = "Gene"
      ),
      annotation_name_side = "top",
      annotation_name_rot = c(
        Cluster = 45,
        `Max mean %` = 45,
        log2FC = 45,
        Significance = 45,
        Gene = 45
      ),
      annotation_name_gp = grid::gpar(
        fontsize = 8,
        fontface = "bold"
      ),
      annotation_legend_param = list(
        Cluster = list(
          title = "Source cluster",
          labels = cluster_legend_labels,
          at = displayed_clusters,
          ncol = 1
        )
      )
    )
  } else {
    right_annotation <- ComplexHeatmap::rowAnnotation(
      Cluster = row_table$cluster_id,
      `Max mean %` = manual_heatmap_percent_annotation(
        values = row_table$max_mean_percent,
        font_size = annotation_font_size,
        colour_function = percent_colour_function
      ),
      Gene = ComplexHeatmap::anno_text(
        row_table$row_id,
        location = grid::unit(0, "npc"),
        just = "left",
        gp = grid::gpar(
          col = "black",
          fontsize = annotation_font_size
        ),
        width = gene_annotation_width
      ),
      col = list(
        Cluster = cluster_annotation_colors
      ),
      simple_anno_size = grid::unit(4.2, "mm"),
      gap = grid::unit(0.6, "mm"),
      show_annotation_name = c(
        Cluster = TRUE,
        `Max mean %` = TRUE,
        Gene = TRUE
      ),
      annotation_label = c(
        Cluster = "Cluster",
        `Max mean %` = "Max mean %",
        Gene = "Gene"
      ),
      annotation_name_side = "top",
      annotation_name_rot = c(
        Cluster = 45,
        `Max mean %` = 45,
        Gene = 45
      ),
      annotation_name_gp = grid::gpar(
        fontsize = 8,
        fontface = "bold"
      ),
      annotation_legend_param = list(
        Cluster = list(
          title = "Source cluster",
          labels = cluster_legend_labels,
          at = displayed_clusters,
          ncol = 1
        )
      )
    )
  }

  title <- paste0(
    "Manual gene-expression heatmap · ",
    length(heatmap_state$found),
    "/",
    length(heatmap_state$requested),
    " requested genes found"
  )

  heatmap <- ComplexHeatmap::Heatmap(
    zscore_matrix,
    name = "Row Z-score",
    col = heatmap_colour_function,
    top_annotation = top_annotation,
    right_annotation = right_annotation,
    cluster_columns = FALSE,
    column_split = sample_metadata$combined_group,
    cluster_column_slices = FALSE,
    column_gap = grid::unit(1.2, "mm"),
    cluster_rows = heatmap_state$row_hclust,
    row_dend_reorder = FALSE,
    row_gap = grid::unit(0, "mm"),
    row_dend_width = grid::unit(18, "mm"),
    show_column_names = TRUE,
    column_names_side = "top",
    column_names_rot = 45,
    column_names_centered = TRUE,
    column_names_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    column_names_max_height = grid::unit(18, "mm"),
    show_row_names = FALSE,
    na_col = "grey90",
    rect_gp = grid::gpar(col = NA),
    border = FALSE,
    use_raster = FALSE,
    width = grid::unit(
      heatmap_state$heatmap_body_width_mm,
      "mm"
    ),
    height = grid::unit(
      heatmap_state$heatmap_body_height_mm,
      "mm"
    ),
    column_title = title,
    column_title_side = "top",
    column_title_gp = grid::gpar(
      fontsize = 12,
      fontface = "bold"
    ),
    heatmap_legend_param = list(
      title = "Row Z-score",
      at = seq(-zscore_limit, zscore_limit, length.out = 5),
      legend_height = grid::unit(42, "mm")
    )
  )

  percent_legend <- ComplexHeatmap::Legend(
    title = "Max mean positive spots (%)",
    col_fun = percent_colour_function,
    at = c(0, 25, 50, 75, 100)
  )

  significance_legend <- if (has_significance) {
    manual_heatmap_significance_legend(
      label = heatmap_state$significance_label,
      metric = heatmap_state$significance_metric,
      threshold_one = heatmap_state$significance_threshold_one,
      threshold_two = heatmap_state$significance_threshold_two,
      threshold_three = heatmap_state$significance_threshold_three,
      star_size = heatmap_state$significance_star_size
    )
  } else {
    NULL
  }

  list(
    heatmap = heatmap,
    percent_legend = percent_legend,
    significance_legend = significance_legend
  )
}

draw_manual_gene_heatmap <- function(
  heatmap_state,
  cluster_names,
  cluster_colors
) {
  objects <- build_manual_gene_heatmap(
    heatmap_state = heatmap_state,
    cluster_names = cluster_names,
    cluster_colors = cluster_colors
  )

  extra_legends <- list(objects$percent_legend)
  if (!is.null(objects$significance_legend)) {
    extra_legends <- c(
      extra_legends,
      list(objects$significance_legend)
    )
  }

  ComplexHeatmap::draw(
    objects$heatmap,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE,
    legend_grouping = "original",
    annotation_legend_list = extra_legends,
    padding = grid::unit(c(3, 3, 3, 3), "mm")
  )

  invisible(TRUE)
}

save_manual_gene_heatmap_pdf <- function(
  file,
  heatmap_state,
  cluster_names,
  cluster_colors
) {
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(
      filename = file,
      width = heatmap_state$pdf_width_inches,
      height = heatmap_state$pdf_height_inches,
      onefile = TRUE
    )
  } else {
    grDevices::pdf(
      file = file,
      width = heatmap_state$pdf_width_inches,
      height = heatmap_state$pdf_height_inches,
      onefile = TRUE,
      useDingbats = FALSE
    )
  }

  on.exit(grDevices::dev.off(), add = TRUE)

  draw_manual_gene_heatmap(
    heatmap_state = heatmap_state,
    cluster_names = cluster_names,
    cluster_colors = cluster_colors
  )

  invisible(file)
}


manual_heatmap_matrix_display_table <- function(
  heatmap_state,
  matrix_name
) {
  if (!isTRUE(heatmap_state$valid)) {
    return(tibble::tibble())
  }

  matrix_object <- heatmap_state[[matrix_name]]
  if (is.null(matrix_object)) {
    return(tibble::tibble())
  }

  matrix_object <- as.matrix(matrix_object)
  row_table <- heatmap_state$row_table

  if (nrow(matrix_object) != nrow(row_table)) {
    stop(
      "Heatmap matrix and row metadata have different numbers of rows.",
      call. = FALSE
    )
  }

  output <- data.frame(
    cluster_id = paste0("C", row_table$cluster_id),
    gene = row_table$gene,
    matrix_object,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  tibble::as_tibble(output)
}

manual_heatmap_significance_choices <- function(edger_results) {
  catalog <- list(
    "Donor group main effect" = c(
      "Overall_Group_ASD_vs_Neurotypical"
    ),
    "Sex main effect" = c(
      "Overall_Sex_Female_vs_Male"
    ),
    "Group × Sex interaction" = c(
      "Interaction"
    ),
    "ASD vs Neurotypical — Male" = c(
      "ASD_vs_Neurotypical_in_Male",
      "ASD_Male_vs_Neurotypical_Male"
    ),
    "ASD vs Neurotypical — Female" = c(
      "ASD_vs_Neurotypical_in_Female",
      "ASD_Female_vs_Neurotypical_Female"
    ),
    "Female vs Male — Neurotypical" = c(
      "Neurotypical_Female_vs_Neurotypical_Male",
      "Female_vs_Male_in_Neurotypical"
    ),
    "Female vs Male — ASD" = c(
      "ASD_Female_vs_ASD_Male",
      "Female_vs_Male_in_ASD"
    )
  )

  available_ids <- names(edger_results)
  choices <- c("None" = "none")

  for (label in names(catalog)) {
    matches <- catalog[[label]][catalog[[label]] %in% available_ids]
    if (length(matches) > 0L) {
      choices[[label]] <- matches[[1]]
    }
  }

  choices
}

manual_heatmap_significance_label <- function(
  test_id,
  edger_results
) {
  test_id <- if (is.null(test_id) || length(test_id) == 0L) {
    ""
  } else {
    as.character(test_id[[1]])
  }
  if (tolower(test_id) == "none") {
    test_id <- ""
  }

  if (!nzchar(test_id)) {
    return("None")
  }

  choices <- manual_heatmap_significance_choices(edger_results)
  matched_label <- names(choices)[unname(choices) == test_id]

  if (length(matched_label) == 0L) {
    return(test_id)
  }

  matched_label[[1]]
}

manual_heatmap_add_significance <- function(
  heatmap_state,
  test_id,
  edger_results,
  metric = "FDR",
  threshold_one = 0.05,
  threshold_two = 0.01,
  threshold_three = 0.001,
  star_size = 10
) {
  if (!isTRUE(heatmap_state$valid)) {
    return(heatmap_state)
  }

  test_id <- if (is.null(test_id) || length(test_id) == 0L) {
    ""
  } else {
    as.character(test_id[[1]])
  }
  if (tolower(test_id) == "none") {
    test_id <- ""
  }

  metric <- as.character(metric[[1]])
  if (!metric %in% c("FDR", "PValue")) {
    stop(
      "`metric` must be either 'FDR' or 'PValue'.",
      call. = FALSE
    )
  }

  thresholds <- c(
    threshold_one,
    threshold_two,
    threshold_three
  )
  thresholds <- suppressWarnings(as.numeric(thresholds))

  if (
    length(thresholds) != 3L ||
      any(!is.finite(thresholds)) ||
      any(thresholds < 0) ||
      any(thresholds > 1)
  ) {
    stop(
      "Significance-star thresholds must be finite values between 0 and 1.",
      call. = FALSE
    )
  }

  if (!(threshold_three <= threshold_two && threshold_two <= threshold_one)) {
    stop(
      "Significance-star thresholds must satisfy: *** <= ** <= *.",
      call. = FALSE
    )
  }

  star_size <- suppressWarnings(as.numeric(star_size[[1]]))
  if (!is.finite(star_size) || star_size <= 0) {
    stop(
      "Significance-star font size must be a positive number.",
      call. = FALSE
    )
  }

  state <- heatmap_state
  state$significance_test_id <- test_id
  state$significance_label <- manual_heatmap_significance_label(
    test_id = test_id,
    edger_results = edger_results
  )
  state$significance_metric <- metric
  state$significance_threshold_one <- threshold_one
  state$significance_threshold_two <- threshold_two
  state$significance_threshold_three <- threshold_three
  state$significance_star_size <- star_size
  state$row_table$significance_value <- NA_real_
  state$row_table$significance_logfc <- NA_real_
  state$row_table$significance_symbol <- ""

  if (!nzchar(test_id)) {
    return(state)
  }

  if (!test_id %in% names(edger_results)) {
    stop(
      "Selected significance table is unavailable: ",
      test_id,
      call. = FALSE
    )
  }

  result_table <- tibble::as_tibble(edger_results[[test_id]])
  required_columns <- c(
    "cluster_id",
    "ensembl_gene_id",
    "logFC",
    metric
  )
  missing_columns <- setdiff(
    required_columns,
    colnames(result_table)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Selected significance table is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  significance_table <- result_table |>
    dplyr::transmute(
      cluster_id = as.character(.data$cluster_id),
      ensembl_gene_id = as.character(.data$ensembl_gene_id),
      significance_value = suppressWarnings(
        as.numeric(.data[[metric]])
      ),
      significance_logfc = suppressWarnings(
        as.numeric(.data$logFC)
      )
    ) |>
    dplyr::distinct(
      .data$cluster_id,
      .data$ensembl_gene_id,
      .keep_all = TRUE
    )

  state$row_table <- state$row_table |>
    dplyr::select(
      -dplyr::any_of(c(
        "significance_value",
        "significance_logfc",
        "significance_symbol"
      ))
    ) |>
    dplyr::left_join(
      significance_table,
      by = c("cluster_id", "ensembl_gene_id")
    ) |>
    dplyr::mutate(
      significance_symbol = dplyr::case_when(
        is.finite(.data$significance_value) &
          .data$significance_value < threshold_three ~ "***",
        is.finite(.data$significance_value) &
          .data$significance_value < threshold_two ~ "**",
        is.finite(.data$significance_value) &
          .data$significance_value < threshold_one ~ "*",
        TRUE ~ ""
      )
    )

  state
}


manual_heatmap_apply_significance_row_filter <- function(
  heatmap_state,
  filter_significant_rows = FALSE,
  minimum_abs_logfc = 0,
  minimum_positive_percent = 0,
  maximum_positive_percent = 100
) {
  if (!isTRUE(heatmap_state$valid)) {
    return(heatmap_state)
  }

  minimum_abs_logfc <- suppressWarnings(
    as.numeric(minimum_abs_logfc[[1]])
  )
  if (!is.finite(minimum_abs_logfc) || minimum_abs_logfc < 0) {
    stop(
      "Minimum |log2FC| for significance filtering must be a finite non-negative number.",
      call. = FALSE
    )
  }

  minimum_positive_percent <- suppressWarnings(
    as.numeric(minimum_positive_percent[[1]])
  )
  maximum_positive_percent <- suppressWarnings(
    as.numeric(maximum_positive_percent[[1]])
  )

  if (
    !is.finite(minimum_positive_percent) ||
      !is.finite(maximum_positive_percent) ||
      minimum_positive_percent < 0 ||
      maximum_positive_percent > 100 ||
      minimum_positive_percent > maximum_positive_percent
  ) {
    stop(
      "Positive-spot thresholds must satisfy 0 <= minimum <= maximum <= 100.",
      call. = FALSE
    )
  }

  state <- heatmap_state
  state$significance_filter_enabled <- isTRUE(filter_significant_rows)
  state$significance_min_abs_logfc <- minimum_abs_logfc
  state$significance_min_positive_percent <- minimum_positive_percent
  state$significance_max_positive_percent <- maximum_positive_percent
  state$significance_total_rows <- nrow(state$row_table)

  gene_key <- if ("requested_gene" %in% colnames(state$row_table)) {
    as.character(state$row_table$requested_gene)
  } else {
    as.character(state$row_table$gene)
  }
  state$significance_total_genes <- length(unique(tolower(gene_key)))

  has_significance <-
    !is.null(state$significance_test_id) &&
      nzchar(state$significance_test_id)

  if (!has_significance) {
    state$significance_matching_rows <- 0L
    state$significance_matching_genes <- 0L

    if (isTRUE(filter_significant_rows)) {
      state$valid <- FALSE
      state$message <- "Select a Statistics metric before enabling significant-row filtering."
    }

    return(state)
  }

  matching_rows <-
    is.finite(state$row_table$significance_value) &
    state$row_table$significance_value < state$significance_threshold_one &
    is.finite(state$row_table$significance_logfc) &
    abs(state$row_table$significance_logfc) >= minimum_abs_logfc &
    is.finite(state$row_table$max_mean_percent) &
    state$row_table$max_mean_percent >= minimum_positive_percent &
    state$row_table$max_mean_percent <= maximum_positive_percent

  state$significance_matching_rows <- sum(matching_rows)
  state$significance_matching_genes <- length(
    unique(tolower(gene_key[matching_rows]))
  )

  if (!isTRUE(filter_significant_rows)) {
    return(state)
  }

  if (!any(matching_rows)) {
    state$valid <- FALSE
    state$message <- paste0(
      "No gene × cluster rows satisfy ",
      state$significance_metric,
      " < ",
      manual_heatmap_format_threshold(state$significance_threshold_one),
      " and |log2FC| ≥ ",
      format(minimum_abs_logfc, trim = TRUE),
      " and Max mean positive spots (%) is between ",
      format(minimum_positive_percent, trim = TRUE),
      " and ",
      format(maximum_positive_percent, trim = TRUE),
      "."
    )
    return(state)
  }

  state$row_table <- state$row_table[
    matching_rows,
    ,
    drop = FALSE
  ]

  matrix_names <- c(
    "raw_counts_matrix",
    "log2_qn_matrix",
    "zscore_matrix"
  )
  for (matrix_name in matrix_names) {
    state[[matrix_name]] <- state[[matrix_name]][
      matching_rows,
      ,
      drop = FALSE
    ]
  }

  state$row_hclust <- manual_heatmap_row_hclust(
    expression_matrix = state$zscore_matrix,
    cluster_labels = state$row_table$cluster_id,
    same_cluster_bonus = state$same_cluster_bonus
  )

  state$heatmap_body_height_mm <-
    nrow(state$zscore_matrix) * state$row_height_mm
  state$heatmap_body_width_mm <-
    ncol(state$zscore_matrix) * 6.2
  state$browser_height_px <- max(
    900,
    round((state$heatmap_body_height_mm + 125) * 3.78)
  )
  state$pdf_height_inches <- max(
    8,
    (state$heatmap_body_height_mm + 90) / 25.4
  )

  state
}

manual_heatmap_format_threshold <- function(value) {
  value <- suppressWarnings(as.numeric(value[[1]]))

  if (!is.finite(value)) {
    return("NA")
  }

  format(
    value,
    scientific = FALSE,
    trim = TRUE,
    digits = 8
  )
}

manual_heatmap_significance_legend <- function(
  label,
  metric,
  threshold_one,
  threshold_two,
  threshold_three,
  star_size
) {
  ComplexHeatmap::Legend(
    title = paste0(metric, " significance — ", label),
    labels = c(
      paste0(
        metric,
        " < ",
        manual_heatmap_format_threshold(threshold_one)
      ),
      paste0(
        metric,
        " < ",
        manual_heatmap_format_threshold(threshold_two)
      ),
      paste0(
        metric,
        " < ",
        manual_heatmap_format_threshold(threshold_three)
      )
    ),
    graphics = list(
      function(x, y, w, h) {
        grid::grid.text(
          "*",
          x = x,
          y = y,
          gp = grid::gpar(
            fontsize = star_size,
            fontface = "bold"
          )
        )
      },
      function(x, y, w, h) {
        grid::grid.text(
          "**",
          x = x,
          y = y,
          gp = grid::gpar(
            fontsize = star_size,
            fontface = "bold"
          )
        )
      },
      function(x, y, w, h) {
        grid::grid.text(
          "***",
          x = x,
          y = y,
          gp = grid::gpar(
            fontsize = star_size,
            fontface = "bold"
          )
        )
      }
    ),
    grid_width = grid::unit(9, "mm"),
    grid_height = grid::unit(5, "mm"),
    title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 8)
  )
}

manual_heatmap_statistics_test_id <- function(
  effect,
  edger_results
) {
  candidates <- switch(
    as.character(effect),
    donor = c(
      "Overall_Group_ASD_vs_Neurotypical"
    ),
    sex = c(
      "Overall_Sex_Female_vs_Male"
    ),
    interaction = c(
      "Interaction"
    ),
    character(0)
  )

  matched <- candidates[candidates %in% names(edger_results)]

  if (length(matched) == 0L) {
    stop(
      "Could not find the requested edgeR statistics table. Available tables: ",
      paste(names(edger_results), collapse = ", "),
      call. = FALSE
    )
  }

  matched[[1]]
}

manual_heatmap_statistics_display_table <- function(
  heatmap_state,
  effect,
  edger_results
) {
  if (!isTRUE(heatmap_state$valid)) {
    return(tibble::tibble())
  }

  test_id <- manual_heatmap_statistics_test_id(
    effect = effect,
    edger_results = edger_results
  )

  result_table <- tibble::as_tibble(edger_results[[test_id]]) |>
    dplyr::mutate(
      cluster_id = as.character(.data$cluster_id),
      ensembl_gene_id = as.character(.data$ensembl_gene_id),
      gene = dplyr::if_else(
        is.na(.data$gene) | .data$gene == "",
        .data$ensembl_gene_id,
        as.character(.data$gene)
      )
    )

  keys <- heatmap_state$row_table |>
    dplyr::select(
      .data$cluster_id,
      .data$ensembl_gene_id
    ) |>
    dplyr::distinct()

  result_table <- result_table |>
    dplyr::inner_join(
      keys,
      by = c("cluster_id", "ensembl_gene_id")
    )

  if (nrow(result_table) == 0L) {
    return(tibble::tibble())
  }

  preferred_columns <- c(
    "cluster_id",
    "gene",
    "logFC",
    "F",
    "PValue",
    "FDR",
    "regulation"
  )

  displayed_columns <- intersect(
    preferred_columns,
    colnames(result_table)
  )

  output <- result_table |>
    dplyr::select(dplyr::all_of(displayed_columns)) |>
    dplyr::mutate(
      cluster_id = paste0("C", .data$cluster_id)
    ) |>
    dplyr::arrange(
      as.integer(sub("^C", "", .data$cluster_id)),
      .data$gene
    )

  names(output)[names(output) == "cluster_id"] <- "Cluster"
  names(output)[names(output) == "gene"] <- "Gene"
  names(output)[names(output) == "logFC"] <- "log2FC"
  names(output)[names(output) == "regulation"] <- "Regulation"

  attr(output, "test_id") <- test_id
  output
}

heatmap_browser_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::fluidRow(
    shiny::column(
      width = 3,
      shiny::div(
        class = "heatmap-browser-sidebar",

        shiny::h4("Heatmap settings"),

        shiny::h4("Clusters"),
        shiny::div(
          class = "heatmap-cluster-checkboxes",
          shiny::checkboxGroupInput(
            inputId = ns("clusters"),
            label = NULL,
            choices = NULL,
            selected = NULL
          )
        ),

        shiny::fluidRow(
          shiny::column(
            width = 6,
            shiny::actionButton(
              inputId = ns("select_all_clusters"),
              label = "Select all",
              width = "100%"
            )
          ),
          shiny::column(
            width = 6,
            shiny::actionButton(
              inputId = ns("clear_clusters"),
              label = "Clear",
              width = "100%"
            )
          )
        ),
        shiny::div(
          class = "parameter-note",
          "Only selected source clusters are included."
        ),

        shiny::hr(),

        shiny::numericInput(
          inputId = ns("same_cluster_bonus"),
          label = "Same-cluster clustering bonus",
          value = 0,
          min = 0,
          max = 0.9,
          step = 0.05,
          width = "100%"
        ),
        shiny::div(
          class = "parameter-note",
          "0 = no preference. Higher values promote neighboring rows from the same cluster."
        ),

        shiny::fluidRow(
          shiny::column(
            width = 6,
            shiny::numericInput(
              inputId = ns("zscore_limit"),
              label = "Row Z-score scale",
              value = 2,
              min = 0.5,
              max = 6,
              step = 0.25,
              width = "100%"
            ),
            shiny::div(
              class = "parameter-note",
              "Symmetric colour/clipping limit."
            )
          ),
          shiny::column(
            width = 6,
            shiny::numericInput(
              inputId = ns("star_font_size"),
              label = "Star size",
              value = 10,
              min = 5,
              max = 24,
              step = 1,
              width = "100%"
            ),
            shiny::div(
              class = "parameter-note",
              "Significance-star font size."
            )
          )
        ),

        shiny::selectInput(
          inputId = ns("significance_test"),
          label = "Statistics metric",
          choices = c("None" = "none"),
          selected = "none",
          width = "100%"
        ),
        shiny::div(
          class = "parameter-note",
          "The selected test supplies per-row log2FC and optional significance stars."
        ),

        shiny::selectInput(
          inputId = ns("significance_metric"),
          label = "Significance based on",
          choices = c(
            "FDR" = "FDR",
            "PValue" = "PValue"
          ),
          selected = "FDR",
          width = "100%"
        ),

        shiny::div(
          class = "heatmap-star-threshold-row",
          shiny::div(
            class = "heatmap-star-threshold-pair",
            shiny::span(class = "heatmap-star-label", "*"),
            shiny::numericInput(
              inputId = ns("star_threshold_one"),
              label = NULL,
              value = 0.05,
              min = 0,
              max = 1,
              step = 0.005,
              width = "100%"
            )
          ),
          shiny::div(
            class = "heatmap-star-threshold-pair",
            shiny::span(class = "heatmap-star-label", "**"),
            shiny::numericInput(
              inputId = ns("star_threshold_two"),
              label = NULL,
              value = 0.01,
              min = 0,
              max = 1,
              step = 0.001,
              width = "100%"
            )
          ),
          shiny::div(
            class = "heatmap-star-threshold-pair",
            shiny::span(class = "heatmap-star-label", "***"),
            shiny::numericInput(
              inputId = ns("star_threshold_three"),
              label = NULL,
              value = 0.001,
              min = 0,
              max = 1,
              step = 0.001,
              width = "100%"
            )
          )
        ),
        shiny::div(
          class = "parameter-note",
          "Star thresholds: *** ≤ ** ≤ *."
        ),

        shiny::div(
          class = "heatmap-threshold-filter-toggle",
          shiny::checkboxInput(
            inputId = ns("show_only_significant"),
            label = "Show only genes by thresholds",
            value = FALSE
          )
        ),
        shiny::div(
          class = "parameter-note heatmap-threshold-filter-toggle-note",
          "Switch on to display only gene × cluster rows meeting all thresholds below."
        ),

        shiny::fluidRow(
          shiny::column(
            width = 4,
            shiny::numericInput(
              inputId = ns("significance_min_positive_percent"),
              label = "Min positive spots (%)",
              value = 0,
              min = 0,
              max = 100,
              step = 1,
              width = "100%"
            )
          ),
          shiny::column(
            width = 4,
            shiny::numericInput(
              inputId = ns("significance_max_positive_percent"),
              label = "Max positive spots (%)",
              value = 100,
              min = 0,
              max = 100,
              step = 1,
              width = "100%"
            )
          ),
          shiny::column(
            width = 4,
            shiny::numericInput(
              inputId = ns("significance_logfc_threshold"),
              label = "Minimum |log2FC|",
              value = 0,
              min = 0,
              step = 0.1,
              width = "100%"
            )
          )
        ),
        shiny::div(
          class = "parameter-note",
          "Thresholds use the selected FDR/PValue, |log2FC|, and the Max mean positive spots (%) annotation."
        ),
        shiny::uiOutput(ns("significance_filter_summary")),

        shiny::textAreaInput(
          inputId = ns("genes"),
          label = "Genes",
          value = "",
          placeholder = "Drd2\nCnr1\nGrin2a",
          rows = 8,
          resize = "vertical",
          width = "100%"
        ),
        shiny::div(
          class = "parameter-note",
          "Enter gene symbols separated by new lines, spaces, commas or semicolons."
        ),

        shiny::uiOutput(ns("gene_status")),

        shiny::br(),

        shiny::div(
          class = "heatmap-action-grid",
          shiny::actionButton(
            inputId = ns("generate"),
            label = "Generate heatmap",
            class = "btn-primary",
            width = "100%"
          ),
          shiny::downloadButton(
            outputId = ns("download_heatmap"),
            label = "Download heatmap",
            class = "btn-default",
            width = "100%"
          )
        )
      )
    ),

    shiny::column(
      width = 9,
      shiny::div(
        class = "heatmap-browser-main",
        shiny::tabsetPanel(
          id = ns("heatmap_view"),
          type = "tabs",

          shiny::tabPanel(
            title = "Heatmap",
            shiny::div(
              class = "heatmap-inner-tab",
              shiny::uiOutput(ns("summary")),
              shiny::div(
                class = "heatmap-browser-description",
                "Rows are gene × source-cluster pseudobulks. Columns are fixed and are not clustered: Male Neurotypical → Male ASD → Female Neurotypical → Female ASD. ",
                "Max mean % is the larger of the Neurotypical and ASD mean percentages of spots with detected expression in that source cluster."
              ),
              shiny::plotOutput(
                outputId = ns("heatmap_plot"),
                width = "100%"
              )
            )
          ),

          shiny::tabPanel(
            title = "Statistics",
            shiny::div(
              class = "heatmap-inner-tab",
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  shiny::selectInput(
                    inputId = ns("statistics_effect"),
                    label = "Statistics",
                    choices = c(
                      "Donor group main effect" = "donor",
                      "Sex main effect" = "sex",
                      "Group × Sex interaction" = "interaction"
                    ),
                    selected = "donor",
                    width = "100%"
                  )
                )
              ),
              shiny::uiOutput(ns("statistics_summary")),
              DT::DTOutput(ns("statistics_table"))
            )
          ),

          shiny::tabPanel(
            title = "Raw counts",
            shiny::div(
              class = "heatmap-inner-tab",
              shiny::div(
                class = "heatmap-data-note",
                "Raw pseudobulk UMI counts. Each row is a gene × source cluster; sample columns contain counts summed across spots from that sample and cluster."
              ),
              DT::DTOutput(ns("raw_counts_table"))
            )
          ),

          shiny::tabPanel(
            title = "Z-score input",
            shiny::div(
              class = "heatmap-inner-tab",
              shiny::div(
                class = "heatmap-data-note",
                "Values used directly to calculate the row Z-score are loaded from the precomputed per-cluster matrices: raw pseudobulk counts → log2(count + 1) → quantile normalization across all 16 biological samples, independently within each cluster."
              ),
              DT::DTOutput(ns("zscore_input_table"))
            )
          )
        )
      )
    )
  )
}

heatmap_browser_server <- function(
  id,
  edger_results,
  edger_model_objects,
  per_cluster_expression,
  cluster_ids,
  cluster_names,
  cluster_colors
) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      available_clusters <- sort_cluster_ids(
        Reduce(
          intersect,
          list(
            as.character(cluster_ids),
            names(edger_model_objects),
            names(per_cluster_expression$raw_counts),
            names(per_cluster_expression$log2_qn)
          )
        )
      )

      cluster_labels <- vapply(
        available_clusters,
        function(cluster_id) {
          region <- cluster_names[[cluster_id]]
          if (is.null(region) || is.na(region) || !nzchar(region)) {
            paste0("C", cluster_id)
          } else {
            paste0("C", cluster_id, " — ", region)
          }
        },
        character(1)
      )

      cluster_choices <- stats::setNames(
        available_clusters,
        cluster_labels
      )

      shiny::updateCheckboxGroupInput(
        session = session,
        inputId = "clusters",
        choices = cluster_choices,
        selected = available_clusters
      )

      significance_choices <- manual_heatmap_significance_choices(
        edger_results
      )

      shiny::updateSelectInput(
        session = session,
        inputId = "significance_test",
        choices = significance_choices,
        selected = "none"
      )

      shiny::observeEvent(
        input$select_all_clusters,
        {
          shiny::updateCheckboxGroupInput(
            session = session,
            inputId = "clusters",
            choices = cluster_choices,
            selected = available_clusters
          )
        }
      )

      shiny::observeEvent(
        input$clear_clusters,
        {
          shiny::updateCheckboxGroupInput(
            session = session,
            inputId = "clusters",
            choices = cluster_choices,
            selected = character(0)
          )
        }
      )

      live_gene_status <- shiny::reactive({
        manual_heatmap_gene_availability(
          gene_text = input$genes,
          selected_clusters = input$clusters,
          edger_results = edger_results,
          per_cluster_expression = per_cluster_expression
        )
      })

      output$gene_status <- shiny::renderUI({
        status <- live_gene_status()

        if (length(status$requested) == 0L) {
          return(
            shiny::div(
              class = "heatmap-gene-status heatmap-gene-status-empty",
              "Paste genes to check availability."
            )
          )
        }

        shiny::tagList(
          shiny::div(
            class = "heatmap-gene-status heatmap-gene-found",
            shiny::span(class = "heatmap-status-icon", "✓"),
            paste0(
              " Found (",
              length(status$found),
              "/",
              length(status$requested),
              "): ",
              if (length(status$found) > 0L) {
                paste(status$found, collapse = ", ")
              } else {
                "none"
              }
            )
          ),
          if (length(status$missing) > 0L) {
            shiny::div(
              class = "heatmap-gene-status heatmap-gene-missing",
              shiny::span(class = "heatmap-status-icon", "✕"),
              paste0(
                " Not found: ",
                paste(status$missing, collapse = ", ")
              )
            )
          }
        )
      })

      generated_heatmap <- shiny::eventReactive(
        input$generate,
        {
          shiny::req(input$same_cluster_bonus)
          shiny::req(input$zscore_limit)

          prepare_manual_gene_heatmap(
            gene_text = input$genes,
            selected_clusters = input$clusters,
            same_cluster_bonus = input$same_cluster_bonus,
            zscore_limit = input$zscore_limit,
            edger_results = edger_results,
            edger_model_objects = edger_model_objects,
            per_cluster_expression = per_cluster_expression
          )
        },
        ignoreInit = TRUE
      )

      display_heatmap_state <- shiny::reactive({
        state <- generated_heatmap()

        if (!isTRUE(state$valid)) {
          return(state)
        }

        state <- manual_heatmap_add_significance(
          heatmap_state = state,
          test_id = input$significance_test,
          edger_results = edger_results,
          metric = input$significance_metric,
          threshold_one = input$star_threshold_one,
          threshold_two = input$star_threshold_two,
          threshold_three = input$star_threshold_three,
          star_size = input$star_font_size
        )

        manual_heatmap_apply_significance_row_filter(
          heatmap_state = state,
          filter_significant_rows = input$show_only_significant,
          minimum_abs_logfc = input$significance_logfc_threshold,
          minimum_positive_percent = input$significance_min_positive_percent,
          maximum_positive_percent = input$significance_max_positive_percent
        )
      })

      output$significance_filter_summary <- shiny::renderUI({
        if (is.null(input$generate) || input$generate == 0L) {
          return(
            shiny::div(
              class = "heatmap-significance-filter-summary heatmap-summary-muted",
              "Generate the heatmap to evaluate significant rows."
            )
          )
        }

        state <- display_heatmap_state()

        has_test <-
          !is.null(state$significance_test_id) &&
            nzchar(state$significance_test_id)

        if (!has_test) {
          return(
            shiny::div(
              class = "heatmap-significance-filter-summary heatmap-summary-muted",
              "Statistics metric: None."
            )
          )
        }

        shiny::div(
          class = "heatmap-significance-filter-summary",
          shiny::strong(
            paste0(
              "Significant genes: ",
              state$significance_matching_genes,
              "/",
              state$significance_total_genes
            )
          ),
          paste0(
            " · significant gene × cluster rows: ",
            state$significance_matching_rows,
            "/",
            state$significance_total_rows
          ),
          shiny::br(),
          shiny::span(
            class = "heatmap-significance-filter-criterion",
            paste0(
              state$significance_metric,
              " < ",
              manual_heatmap_format_threshold(
                state$significance_threshold_one
              ),
              " · |log2FC| ≥ ",
              format(
                state$significance_min_abs_logfc,
                trim = TRUE
              ),
              " · positive spots ",
              format(
                state$significance_min_positive_percent,
                trim = TRUE
              ),
              "–",
              format(
                state$significance_max_positive_percent,
                trim = TRUE
              ),
              "%"
            )
          )
        )
      })

      output$summary <- shiny::renderUI({
        if (input$generate == 0L) {
          return(
            shiny::div(
              class = "heatmap-summary heatmap-summary-muted",
              "Enter genes and click Generate heatmap."
            )
          )
        }

        state <- display_heatmap_state()

        if (!isTRUE(state$valid)) {
          return(
            shiny::div(
              class = "heatmap-summary heatmap-summary-error",
              state$message
            )
          )
        }

        shiny::div(
          class = "heatmap-summary",
          shiny::strong(
            paste0(
              length(state$found),
              "/",
              length(state$requested),
              " requested genes found"
            )
          ),
          paste0(
            " · ",
            nrow(state$row_table),
            " gene × cluster rows · ",
            length(unique(state$row_table$cluster_id)),
            " source clusters · row Z-score ±",
            format(state$zscore_limit, trim = TRUE),
            if (nzchar(state$significance_test_id)) {
              paste0(
                " · ",
                state$significance_metric,
                " stars: ",
                state$significance_label,
                " · thresholds ",
                manual_heatmap_format_threshold(
                  state$significance_threshold_one
                ),
                "/",
                manual_heatmap_format_threshold(
                  state$significance_threshold_two
                ),
                "/",
                manual_heatmap_format_threshold(
                  state$significance_threshold_three
                )
              )
            } else {
              ""
            },
            if (isTRUE(state$significance_filter_enabled)) {
              paste0(
                " · displayed significant rows ",
                state$significance_matching_rows,
                "/",
                state$significance_total_rows,
                " at |log2FC| ≥ ",
                format(state$significance_min_abs_logfc, trim = TRUE),
                " and positive spots ",
                format(state$significance_min_positive_percent, trim = TRUE),
                "–",
                format(state$significance_max_positive_percent, trim = TRUE),
                "%"
              )
            } else {
              ""
            }
          )
        )
      })

      output$heatmap_plot <- shiny::renderPlot(
        {
          shiny::validate(
            shiny::need(
              input$generate > 0L,
              "Enter genes and click Generate heatmap."
            )
          )

          state <- display_heatmap_state()

          shiny::validate(
            shiny::need(
              isTRUE(state$valid),
              state$message
            )
          )

          draw_manual_gene_heatmap(
            heatmap_state = state,
            cluster_names = cluster_names,
            cluster_colors = cluster_colors
          )
        },
        height = function() {
          if (is.null(input$generate) || input$generate == 0L) {
            return(850)
          }

          state <- display_heatmap_state()
          if (!isTRUE(state$valid)) {
            return(850)
          }

          state$browser_height_px
        },
        res = 96
      )

      selected_statistics <- shiny::reactive({
        shiny::validate(
          shiny::need(
            input$generate > 0L,
            "Generate the gene set first."
          )
        )

        state <- generated_heatmap()
        shiny::validate(
          shiny::need(isTRUE(state$valid), state$message)
        )

        manual_heatmap_statistics_display_table(
          heatmap_state = state,
          effect = input$statistics_effect,
          edger_results = edger_results
        )
      })

      output$statistics_summary <- shiny::renderUI({
        if (input$generate == 0L) {
          return(
            shiny::div(
              class = "heatmap-summary heatmap-summary-muted",
              "Generate the gene set first."
            )
          )
        }

        state <- generated_heatmap()
        if (!isTRUE(state$valid)) {
          return(
            shiny::div(
              class = "heatmap-summary heatmap-summary-error",
              state$message
            )
          )
        }

        test_id <- tryCatch(
          manual_heatmap_statistics_test_id(
            effect = input$statistics_effect,
            edger_results = edger_results
          ),
          error = function(e) NA_character_
        )

        shiny::div(
          class = "heatmap-summary",
          if (is.na(test_id)) {
            "Selected statistics table is unavailable."
          } else {
            paste0("edgeR test_id: ", test_id)
          }
        )
      })

      output$statistics_table <- DT::renderDT(
        {
          table <- selected_statistics()

          shiny::validate(
            shiny::need(
              nrow(table) > 0L,
              "No statistics were found for the generated gene × cluster rows."
            )
          )

          dt <- DT::datatable(
            table,
            rownames = FALSE,
            filter = "top",
            class = "stripe hover compact",
            options = list(
              pageLength = 25,
              lengthChange = TRUE,
              scrollX = TRUE,
              autoWidth = TRUE,
              order = list(),
              dom = "lftip"
            )
          )

          if ("log2FC" %in% colnames(table)) {
            dt <- DT::formatRound(dt, columns = "log2FC", digits = 3)
          }
          if ("F" %in% colnames(table)) {
            dt <- DT::formatRound(dt, columns = "F", digits = 3)
          }
          p_columns <- intersect(c("PValue", "FDR"), colnames(table))
          if (length(p_columns) > 0L) {
            dt <- DT::formatSignif(dt, columns = p_columns, digits = 3)
          }

          dt
        },
        server = FALSE
      )

      output$raw_counts_table <- DT::renderDT(
        {
          shiny::validate(
            shiny::need(
              input$generate > 0L,
              "Generate the gene set first."
            )
          )

          state <- generated_heatmap()
          shiny::validate(
            shiny::need(isTRUE(state$valid), state$message)
          )

          table <- manual_heatmap_matrix_display_table(
            heatmap_state = state,
            matrix_name = "raw_counts_matrix"
          )

          DT::datatable(
            table,
            rownames = FALSE,
            class = "stripe hover compact",
            options = list(
              pageLength = 25,
              lengthChange = TRUE,
              scrollX = TRUE,
              autoWidth = TRUE,
              order = list(),
              dom = "lftip"
            )
          )
        },
        server = FALSE
      )

      output$zscore_input_table <- DT::renderDT(
        {
          shiny::validate(
            shiny::need(
              input$generate > 0L,
              "Generate the gene set first."
            )
          )

          state <- generated_heatmap()
          shiny::validate(
            shiny::need(isTRUE(state$valid), state$message)
          )

          table <- manual_heatmap_matrix_display_table(
            heatmap_state = state,
            matrix_name = "log2_qn_matrix"
          )

          numeric_columns <- names(table)[
            vapply(table, is.numeric, logical(1))
          ]

          dt <- DT::datatable(
            table,
            rownames = FALSE,
            class = "stripe hover compact",
            options = list(
              pageLength = 25,
              lengthChange = TRUE,
              scrollX = TRUE,
              autoWidth = TRUE,
              order = list(),
              dom = "lftip"
            )
          )

          if (length(numeric_columns) > 0L) {
            dt <- DT::formatRound(
              dt,
              columns = numeric_columns,
              digits = 4
            )
          }

          dt
        },
        server = FALSE
      )

      output$download_heatmap <- shiny::downloadHandler(
        filename = function() {
          state <- display_heatmap_state()
          n_genes <- if (isTRUE(state$valid)) length(state$found) else 0L
          paste0(
            "manual_gene_expression_heatmap_",
            n_genes,
            "genes.pdf"
          )
        },
        content = function(file) {
          state <- display_heatmap_state()

          if (!isTRUE(state$valid)) {
            stop(
              "Generate a valid heatmap before downloading.",
              call. = FALSE
            )
          }

          save_manual_gene_heatmap_pdf(
            file = file,
            heatmap_state = state,
            cluster_names = cluster_names,
            cluster_colors = cluster_colors
          )
        },
        contentType = "application/pdf"
      )
    }
  )
}
