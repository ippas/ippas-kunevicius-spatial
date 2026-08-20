#!/usr/bin/env Rscript

# ==============================================================================
# functions_geneExpressionOnSlidesAndBarplots_pseudobulkPerClusterEdgeR.R
#
# FIXED VERSION:
# - barplots are strictly monochrome;
# - bars have white fill and a thick black border;
# - sample points have white fill and a black border;
# - a style marker forces regeneration of older coloured barplots.
#
# Dedicated visualization functions for the cluster-specific pseudobulk edgeR
# workflow:
#
#   expression ~ fmt_donor_group * sex
#
# This module is intentionally separate from marker-gene visualization code.
# It visualizes genes selected from the edgeR main donor-group effect and adds
# the corresponding cluster-specific edgeR statistics to the barplots.
#
# For each selected gene, the module creates:
#   1. per-spot expression using one red scale;
#   2. sample-by-cluster aggregated expression using one red scale;
#   3. per-spot expression with significant clusters in red and all remaining
#      clusters in green;
#   4. sample-by-cluster aggregated expression with significant clusters in red
#      and all remaining clusters in green;
#   5. one four-group barplot for every significant cluster.
#
# Barplot group order:
#   Male Neurotypical, Male ASD, Female Neurotypical, Female ASD.
#
# Barplot statistics:
#   - title: main donor-group effect log2FC, nominal P-value and FDR;
#   - brackets: ASD vs Neurotypical within Male and Female separately, showing
#     nominal P-value and FDR;
#   - x-axis: pooled percentage of spots with raw count > 0 under every bar.
# ============================================================================== 


# ==============================================================================
# 1. Package checks and general helpers
# ==============================================================================

check_required_packages_edger_group_visualization <- function() {

  required_packages <- c(
    "Seurat",
    "SeuratObject",
    "Matrix",
    "dplyr",
    "tibble",
    "ggplot2",
    "patchwork",
    "scales",
    "grid"
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
      "Missing required packages: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


write_edger_group_visualization_tsv <- function(data, filename) {

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


sanitize_edger_group_file_component <- function(x) {

  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)

  ifelse(x == "", "unnamed", x)
}


sort_edger_group_cluster_ids <- function(cluster_ids) {

  cluster_ids <- unique(as.character(cluster_ids))
  numeric_ids <- suppressWarnings(as.numeric(cluster_ids))

  if (length(cluster_ids) == 0L) {
    return(character(0))
  }

  if (!anyNA(numeric_ids)) {
    return(cluster_ids[order(numeric_ids)])
  }

  sort(cluster_ids)
}


format_edger_group_number <- function(x, digits = 3L) {

  if (length(x) == 0L || is.na(x[[1]]) || !is.finite(x[[1]])) {
    return("NA")
  }

  formatC(
    as.numeric(x[[1]]),
    format = "f",
    digits = as.integer(digits)
  )
}


format_edger_group_probability <- function(x) {

  if (length(x) == 0L || is.na(x[[1]]) || !is.finite(x[[1]])) {
    return("NA")
  }

  x <- as.numeric(x[[1]])

  if (x < 0.001) {
    return(formatC(x, format = "e", digits = 2))
  }

  formatC(x, format = "f", digits = 3)
}


standardize_edger_group_sex <- function(x) {

  x_character <- trimws(as.character(x))
  x_lower <- tolower(x_character)

  output <- rep(NA_character_, length(x_lower))
  output[x_lower %in% c("m", "male")] <- "Male"
  output[x_lower %in% c("f", "female")] <- "Female"

  unsupported <- unique(x_character[is.na(output)])

  if (length(unsupported) > 0L) {
    stop(
      "Unsupported sex label(s): ",
      paste(unsupported, collapse = ", "),
      call. = FALSE
    )
  }

  output
}


standardize_edger_group_donor_group <- function(x) {

  x_character <- trimws(as.character(x))
  x_lower <- tolower(x_character)

  output <- rep(NA_character_, length(x_lower))

  output[grepl("asd|autism", x_lower)] <- "ASD"
  output[
    grepl(
      "neurotypical|control|typical|(^|[^a-z])nt([^a-z]|$)",
      x_lower
    )
  ] <- "Neurotypical"

  unsupported <- unique(x_character[is.na(output)])

  if (length(unsupported) > 0L) {
    stop(
      "Unsupported donor-group label(s): ",
      paste(unsupported, collapse = ", "),
      call. = FALSE
    )
  }

  output
}


make_edger_group_key <- function(sex_std, group_std) {

  paste0(
    tolower(sex_std),
    "_",
    tolower(group_std)
  )
}


build_significant_cluster_label <- function(cluster_ids) {

  cluster_ids <- sort_edger_group_cluster_ids(cluster_ids)

  paste0(
    "significantClusters_",
    paste0("C", sanitize_edger_group_file_component(cluster_ids), collapse = "-")
  )
}


# ==============================================================================
# 2. Input loading
# ==============================================================================

load_single_seurat_for_edger_group_visualization <- function(
    input_rdata_file,
    requested_object_name = NULL
) {

  if (!file.exists(input_rdata_file)) {
    stop(
      "Input Seurat RData file does not exist:\n",
      input_rdata_file,
      call. = FALSE
    )
  }

  load_environment <- new.env(parent = globalenv())
  loaded_object_names <- load(input_rdata_file, envir = load_environment)

  if (!is.null(requested_object_name)) {

    if (!requested_object_name %in% loaded_object_names) {
      stop(
        "Requested Seurat object was not found: ",
        requested_object_name,
        call. = FALSE
      )
    }

    selected_object <- get(
      requested_object_name,
      envir = load_environment,
      inherits = FALSE
    )

    if (!inherits(selected_object, "Seurat")) {
      stop(
        "Requested object is not a Seurat object: ",
        requested_object_name,
        call. = FALSE
      )
    }

    return(
      list(
        object = selected_object,
        object_name = requested_object_name
      )
    )
  }

  seurat_object_names <- loaded_object_names[
    vapply(
      loaded_object_names,
      function(object_name) {
        inherits(
          get(
            object_name,
            envir = load_environment,
            inherits = FALSE
          ),
          "Seurat"
        )
      },
      FUN.VALUE = logical(1)
    )
  ]

  if (length(seurat_object_names) != 1L) {
    stop(
      "Expected exactly one Seurat object in the input RData file. Found: ",
      paste(seurat_object_names, collapse = ", "),
      call. = FALSE
    )
  }

  selected_name <- seurat_object_names[[1]]

  list(
    object = get(
      selected_name,
      envir = load_environment,
      inherits = FALSE
    ),
    object_name = selected_name
  )
}


load_edger_results_for_group_visualization <- function(edger_results_rdata_file) {

  if (!file.exists(edger_results_rdata_file)) {
    stop(
      "edgeR results RData file does not exist:\n",
      edger_results_rdata_file,
      call. = FALSE
    )
  }

  result_environment <- new.env(parent = globalenv())
  loaded_names <- load(edger_results_rdata_file, envir = result_environment)

  required_objects <- c(
    "edgeR_perCluster_combinedResults",
    "edgeR_perCluster_testDefinitions"
  )

  missing_objects <- setdiff(required_objects, loaded_names)

  if (length(missing_objects) > 0L) {
    stop(
      "Missing required edgeR object(s) in RData: ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  combined_results <- get(
    "edgeR_perCluster_combinedResults",
    envir = result_environment,
    inherits = FALSE
  )

  test_definitions <- get(
    "edgeR_perCluster_testDefinitions",
    envir = result_environment,
    inherits = FALSE
  )

  required_test_ids <- c(
    "Overall_Group_ASD_vs_Neurotypical",
    "ASD_Male_vs_Neurotypical_Male",
    "ASD_Female_vs_Neurotypical_Female"
  )

  missing_test_ids <- setdiff(
    required_test_ids,
    names(combined_results)
  )

  if (length(missing_test_ids) > 0L) {
    stop(
      "Missing required edgeR result table(s): ",
      paste(missing_test_ids, collapse = ", "),
      "\nAvailable result tables: ",
      paste(names(combined_results), collapse = ", "),
      call. = FALSE
    )
  }

  list(
    combined_results = combined_results,
    test_definitions = test_definitions,
    loaded_names = loaded_names
  )
}


# ==============================================================================
# 3. Gene selection from the main donor-group effect
# ==============================================================================

select_main_group_genes_for_visualization <- function(
    combined_results,
    fdr_threshold = 0.05,
    abs_log2fc_threshold = 0.7,
    manual_target_genes = NULL,
    maximum_number_of_genes = 3L
) {

  main_group_results <- combined_results[[
    "Overall_Group_ASD_vs_Neurotypical"
  ]]

  required_result_columns <- c(
    "cluster_id",
    "ensembl_gene_id",
    "gene",
    "logFC",
    "PValue",
    "FDR"
  )

  missing_result_columns <- setdiff(
    required_result_columns,
    colnames(main_group_results)
  )

  if (length(missing_result_columns) > 0L) {
    stop(
      "Missing main-group result column(s): ",
      paste(missing_result_columns, collapse = ", "),
      call. = FALSE
    )
  }

  significant_results <- main_group_results |>
    dplyr::filter(
      is.finite(.data$FDR),
      is.finite(.data$logFC),
      .data$FDR < fdr_threshold,
      abs(.data$logFC) > abs_log2fc_threshold
    ) |>
    dplyr::mutate(
      cluster_id = as.character(.data$cluster_id),
      gene = dplyr::if_else(
        is.na(.data$gene) | .data$gene == "",
        .data$ensembl_gene_id,
        as.character(.data$gene)
      )
    )

  if (nrow(significant_results) == 0L) {
    stop(
      "No genes passed FDR < ",
      fdr_threshold,
      " and |log2FC| > ",
      abs_log2fc_threshold,
      ".",
      call. = FALSE
    )
  }

  gene_selection_table <- significant_results |>
    dplyr::group_by(
      .data$ensembl_gene_id,
      .data$gene
    ) |>
    dplyr::summarise(
      significant_clusters = list(
        sort_edger_group_cluster_ids(.data$cluster_id)
      ),
      number_of_significant_clusters = dplyr::n_distinct(.data$cluster_id),
      minimum_FDR = min(.data$FDR, na.rm = TRUE),
      minimum_PValue = min(.data$PValue, na.rm = TRUE),
      maximum_abs_log2FC = max(abs(.data$logFC), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$minimum_FDR,
      dplyr::desc(.data$maximum_abs_log2FC),
      .data$gene
    )

  if (!is.null(manual_target_genes)) {

    manual_target_genes <- unique(as.character(manual_target_genes))

    gene_selection_table <- gene_selection_table |>
      dplyr::mutate(
        manual_match_key = dplyr::case_when(
          .data$gene %in% manual_target_genes ~ .data$gene,
          .data$ensembl_gene_id %in% manual_target_genes ~ .data$ensembl_gene_id,
          TRUE ~ NA_character_
        )
      ) |>
      dplyr::filter(!is.na(.data$manual_match_key)) |>
      dplyr::mutate(
        manual_order = match(.data$manual_match_key, manual_target_genes)
      ) |>
      dplyr::arrange(.data$manual_order) |>
      dplyr::select(-manual_match_key, -manual_order)

    selected_keys <- c(
      gene_selection_table$gene,
      gene_selection_table$ensembl_gene_id
    )

    missing_manual_genes <- manual_target_genes[
      !manual_target_genes %in% selected_keys
    ]

    if (length(missing_manual_genes) > 0L) {
      stop(
        "Manual gene(s) did not pass the configured thresholds: ",
        paste(missing_manual_genes, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (!is.null(maximum_number_of_genes)) {

    maximum_number_of_genes <- as.integer(maximum_number_of_genes)

    if (
      length(maximum_number_of_genes) != 1L ||
        is.na(maximum_number_of_genes) ||
        maximum_number_of_genes < 1L
    ) {
      stop(
        "`maximum_number_of_genes` must be NULL or one positive integer.",
        call. = FALSE
      )
    }

    gene_selection_table <- utils::head(
      gene_selection_table,
      n = maximum_number_of_genes
    )
  }

  if (nrow(gene_selection_table) == 0L) {
    stop("No genes remain after gene selection.", call. = FALSE)
  }

  selected_gene_ids <- gene_selection_table$ensembl_gene_id

  selected_significant_results <- significant_results |>
    dplyr::filter(.data$ensembl_gene_id %in% selected_gene_ids) |>
    dplyr::mutate(
      cluster_numeric = suppressWarnings(as.numeric(.data$cluster_id))
    ) |>
    dplyr::arrange(
      .data$gene,
      .data$cluster_numeric,
      .data$cluster_id
    ) |>
    dplyr::select(-cluster_numeric)

  list(
    gene_selection_table = gene_selection_table,
    selected_significant_results = selected_significant_results,
    all_significant_results = significant_results,
    main_group_results = main_group_results
  )
}


flatten_selected_gene_table <- function(
    gene_selection_table,
    fdr_threshold,
    abs_log2fc_threshold
) {

  gene_selection_table |>
    dplyr::mutate(
      significant_clusters = vapply(
        .data$significant_clusters,
        paste,
        character(1),
        collapse = ","
      ),
      fdr_threshold = fdr_threshold,
      abs_log2FC_threshold = abs_log2fc_threshold,
      selection_rule = paste0(
        "FDR < ",
        fdr_threshold,
        " and |log2FC| > ",
        abs_log2fc_threshold
      )
    )
}


# ==============================================================================
# 4. Sample layout and spatial-image mapping
# ==============================================================================

build_edger_group_sample_layout <- function(
    seurat_object,
    sample_id_column = "sample_ID",
    group_column = "fmt_donor_group",
    sex_column = "sex",
    included_sample_ids = NULL,
    column_order = c(
      "male_neurotypical",
      "male_asd",
      "female_neurotypical",
      "female_asd"
    ),
    column_titles = c(
      "male_neurotypical" = "Male neurotypical",
      "male_asd" = "Male ASD",
      "female_neurotypical" = "Female neurotypical",
      "female_asd" = "Female ASD"
    )
) {

  metadata_table <- seurat_object[[]]

  required_columns <- c(
    sample_id_column,
    group_column,
    sex_column
  )

  missing_columns <- setdiff(required_columns, colnames(metadata_table))

  if (length(missing_columns) > 0L) {
    stop(
      "Missing metadata column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  sample_order <- unique(as.character(metadata_table[[sample_id_column]]))

  if (!is.null(included_sample_ids)) {
    included_sample_ids <- as.character(included_sample_ids)
    sample_order <- sample_order[sample_order %in% included_sample_ids]

    missing_in_object <- setdiff(included_sample_ids, sample_order)

    if (length(missing_in_object) > 0L) {
      stop(
        "Included sample(s) absent from Seurat object: ",
        paste(missing_in_object, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (length(sample_order) == 0L) {
    stop("No samples remain in the visualization layout.", call. = FALSE)
  }

  sample_rows <- lapply(
    sample_order,
    function(sample_id) {

      sample_metadata <- metadata_table[
        as.character(metadata_table[[sample_id_column]]) == sample_id,
        ,
        drop = FALSE
      ]

      group_values <- unique(as.character(sample_metadata[[group_column]]))
      sex_values <- unique(as.character(sample_metadata[[sex_column]]))

      if (length(group_values) != 1L || length(sex_values) != 1L) {
        stop(
          "Sample ",
          sample_id,
          " does not have one unique donor group and sex.",
          call. = FALSE
        )
      }

      sex_std <- standardize_edger_group_sex(sex_values)
      group_std <- standardize_edger_group_donor_group(group_values)

      data.frame(
        sample_ID = sample_id,
        group_std = group_std,
        sex_std = sex_std,
        group_key = make_edger_group_key(sex_std, group_std),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )

  sample_table <- do.call(rbind, sample_rows)
  rownames(sample_table) <- NULL

  unsupported_keys <- setdiff(unique(sample_table$group_key), column_order)

  if (length(unsupported_keys) > 0L) {
    stop(
      "Unsupported sex-by-group key(s): ",
      paste(unsupported_keys, collapse = ", "),
      call. = FALSE
    )
  }

  samples_by_column <- lapply(
    column_order,
    function(group_key) {
      sample_table$sample_ID[sample_table$group_key == group_key]
    }
  )

  names(samples_by_column) <- column_order

  list(
    sample_table = sample_table,
    sample_order = sample_order,
    samples_by_column = samples_by_column,
    max_rows = max(lengths(samples_by_column)),
    column_order = column_order,
    column_titles = column_titles
  )
}


create_edger_group_image_map <- function(
    seurat_object,
    sample_id_column = "sample_ID"
) {

  image_names <- SeuratObject::Images(seurat_object)

  if (length(image_names) == 0L) {
    stop("No spatial images were found in the Seurat object.", call. = FALSE)
  }

  metadata_table <- seurat_object[[]]

  mapping_rows <- lapply(
    image_names,
    function(image_name) {

      image_cells <- SeuratObject::Cells(seurat_object[[image_name]])

      sample_ids <- unique(
        as.character(
          metadata_table[
            image_cells,
            sample_id_column,
            drop = TRUE
          ]
        )
      )

      sample_ids <- sample_ids[!is.na(sample_ids) & sample_ids != ""]

      if (length(sample_ids) != 1L) {
        stop(
          "Spatial image ",
          image_name,
          " maps to ",
          length(sample_ids),
          " sample IDs.",
          call. = FALSE
        )
      }

      data.frame(
        sample_ID = sample_ids[[1]],
        image_name = image_name,
        n_spots = length(image_cells),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )

  image_map <- do.call(rbind, mapping_rows)
  rownames(image_map) <- NULL

  duplicated_samples <- unique(
    image_map$sample_ID[duplicated(image_map$sample_ID)]
  )

  if (length(duplicated_samples) > 0L) {
    stop(
      "More than one spatial image maps to sample(s): ",
      paste(duplicated_samples, collapse = ", "),
      call. = FALSE
    )
  }

  image_map
}


standardize_edger_group_coordinates <- function(coordinate_table) {

  coordinate_table <- as.data.frame(
    coordinate_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"cell" %in% colnames(coordinate_table)) {
    coordinate_table$cell <- rownames(coordinate_table)
  }

  x_candidates <- c(
    "imagecol",
    "x",
    "col",
    "pxl_col_in_fullres"
  )

  y_candidates <- c(
    "imagerow",
    "y",
    "row",
    "pxl_row_in_fullres"
  )

  x_column <- x_candidates[x_candidates %in% colnames(coordinate_table)][1]
  y_column <- y_candidates[y_candidates %in% colnames(coordinate_table)][1]

  if (is.na(x_column) || is.na(y_column)) {
    stop(
      "Could not identify x/y columns in tissue coordinates. Available: ",
      paste(colnames(coordinate_table), collapse = ", "),
      call. = FALSE
    )
  }

  data.frame(
    cell = as.character(coordinate_table$cell),
    x_plot = as.numeric(coordinate_table[[x_column]]),
    y_plot = as.numeric(coordinate_table[[y_column]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


compute_edger_group_square_limits <- function(
    xmin,
    xmax,
    ymin,
    ymax,
    padding_fraction = 0.03
) {

  x_range <- xmax - xmin
  y_range <- ymax - ymin
  maximum_range <- max(x_range, y_range)

  if (!is.finite(maximum_range) || maximum_range <= 0) {
    maximum_range <- 1
  }

  half_side <- maximum_range / 2 * (1 + padding_fraction)
  x_center <- (xmin + xmax) / 2
  y_center <- (ymin + ymax) / 2

  list(
    x_limits = c(x_center - half_side, x_center + half_side),
    y_limits = c(y_center - half_side, y_center + half_side)
  )
}


get_edger_group_image_dimensions <- function(
    seurat_object,
    image_name,
    coordinate_table
) {

  raw_image <- tryCatch(
    SeuratObject::GetImage(
      object = seurat_object,
      image = image_name,
      mode = "raw"
    ),
    error = function(error_condition) NULL
  )

  raw_dimensions <- dim(raw_image)

  if (!is.null(raw_dimensions) && length(raw_dimensions) >= 2L) {
    return(
      list(
        image_array = raw_image,
        image_height = as.numeric(raw_dimensions[[1]]),
        image_width = as.numeric(raw_dimensions[[2]])
      )
    )
  }

  list(
    image_array = NULL,
    image_height = max(coordinate_table$y_plot, na.rm = TRUE),
    image_width = max(coordinate_table$x_plot, na.rm = TRUE)
  )
}


# ==============================================================================
# 5. Shared spatial/count context
# ==============================================================================

prepare_edger_group_visualization_context <- function(
    seurat_object,
    cluster_column,
    assay_name = "RNA",
    sample_id_column = "sample_ID",
    group_column = "fmt_donor_group",
    sex_column = "sex",
    included_sample_ids = NULL,
    image_scale = "lowres",
    panel_padding_fraction = 0.03
) {

  if (!inherits(seurat_object, "Seurat")) {
    stop("`seurat_object` must be a Seurat object.", call. = FALSE)
  }

  available_assays <- SeuratObject::Assays(seurat_object)

  if (!assay_name %in% available_assays) {
    stop(
      "Assay '",
      assay_name,
      "' is not present. Available assays: ",
      paste(available_assays, collapse = ", "),
      call. = FALSE
    )
  }

  metadata_table <- seurat_object[[]]

  if (!cluster_column %in% colnames(metadata_table)) {
    stop(
      "Missing clustering column: ",
      cluster_column,
      call. = FALSE
    )
  }

  SeuratObject::DefaultAssay(seurat_object) <- assay_name

  count_layers <- SeuratObject::Layers(
    seurat_object[[assay_name]],
    search = "^counts"
  )

  if (length(count_layers) == 0L) {
    stop(
      "No counts layer was found in assay: ",
      assay_name,
      call. = FALSE
    )
  }

  if (length(count_layers) > 1L) {
    seurat_object <- SeuratObject::JoinLayers(
      object = seurat_object,
      assay = assay_name
    )
  }

  counts_matrix <- SeuratObject::LayerData(
    object = seurat_object,
    assay = assay_name,
    layer = "counts"
  )

  if (!inherits(counts_matrix, "Matrix")) {
    counts_matrix <- Matrix::Matrix(counts_matrix, sparse = TRUE)
  }

  total_umi_per_spot <- Matrix::colSums(counts_matrix)

  sample_layout <- build_edger_group_sample_layout(
    seurat_object = seurat_object,
    sample_id_column = sample_id_column,
    group_column = group_column,
    sex_column = sex_column,
    included_sample_ids = included_sample_ids
  )

  image_map <- create_edger_group_image_map(
    seurat_object = seurat_object,
    sample_id_column = sample_id_column
  )

  missing_images <- setdiff(
    sample_layout$sample_order,
    image_map$sample_ID
  )

  if (length(missing_images) > 0L) {
    stop(
      "No spatial image for sample(s): ",
      paste(missing_images, collapse = ", "),
      call. = FALSE
    )
  }

  sample_cache <- lapply(
    sample_layout$sample_order,
    function(sample_id) {

      image_name <- image_map$image_name[
        match(sample_id, image_map$sample_ID)
      ]

      image_cells <- SeuratObject::Cells(seurat_object[[image_name]])

      coordinate_table <- SeuratObject::GetTissueCoordinates(
        object = seurat_object,
        image = image_name,
        scale = image_scale
      )

      coordinate_table <- standardize_edger_group_coordinates(
        coordinate_table
      )

      common_cells <- image_cells[
        image_cells %in% coordinate_table$cell &
          image_cells %in% colnames(counts_matrix)
      ]

      if (length(common_cells) == 0L) {
        stop(
          "No shared cells among image, coordinates and counts for sample ",
          sample_id,
          ".",
          call. = FALSE
        )
      }

      coordinate_table <- coordinate_table[
        match(common_cells, coordinate_table$cell),
        ,
        drop = FALSE
      ]

      image_info <- get_edger_group_image_dimensions(
        seurat_object = seurat_object,
        image_name = image_name,
        coordinate_table = coordinate_table
      )

      frame_limits <- compute_edger_group_square_limits(
        xmin = 0,
        xmax = image_info$image_width,
        ymin = 0,
        ymax = image_info$image_height,
        padding_fraction = panel_padding_fraction
      )

      sample_metadata <- sample_layout$sample_table[
        sample_layout$sample_table$sample_ID == sample_id,
        ,
        drop = FALSE
      ]

      spot_metadata <- metadata_table[
        common_cells,
        ,
        drop = FALSE
      ]

      coordinate_table$sample_ID <- sample_id
      coordinate_table$group_std <- sample_metadata$group_std[[1]]
      coordinate_table$sex_std <- sample_metadata$sex_std[[1]]
      coordinate_table$group_key <- sample_metadata$group_key[[1]]
      coordinate_table$cluster_id <- as.character(
        spot_metadata[[cluster_column]]
      )

      list(
        sample_ID = sample_id,
        image_name = image_name,
        cells = common_cells,
        coordinates = coordinate_table,
        image_array = image_info$image_array,
        image_width = image_info$image_width,
        image_height = image_info$image_height,
        frame_limits = frame_limits
      )
    }
  )

  names(sample_cache) <- sample_layout$sample_order

  list(
    seurat_object = seurat_object,
    counts_matrix = counts_matrix,
    total_umi_per_spot = total_umi_per_spot,
    cluster_column = cluster_column,
    assay_name = assay_name,
    sample_layout = sample_layout,
    image_map = image_map,
    sample_cache = sample_cache,
    cluster_levels = sort_edger_group_cluster_ids(
      metadata_table[[cluster_column]]
    )
  )
}


match_edger_group_gene_feature <- function(
    feature_names,
    target_gene,
    target_ensembl_gene_id
) {

  feature_names <- as.character(feature_names)
  target_gene <- as.character(target_gene)[[1]]
  target_ensembl_gene_id <- as.character(target_ensembl_gene_id)[[1]]

  canonicalize <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("\\.[0-9]+$", "", x)
    x <- sub("-[0-9]+$", "", x)
    x
  }

  candidate_sets <- list(
    ensembl_exact = which(
      canonicalize(feature_names) == canonicalize(target_ensembl_gene_id)
    ),
    symbol_exact = which(feature_names == target_gene),
    symbol_case_insensitive = which(
      tolower(feature_names) == tolower(target_gene)
    ),
    symbol_canonical = which(
      tolower(canonicalize(feature_names)) ==
        tolower(canonicalize(target_gene))
    )
  )

  for (match_type in names(candidate_sets)) {

    matching_indices <- unique(candidate_sets[[match_type]])

    if (length(matching_indices) == 1L) {
      return(
        list(
          matched_index = matching_indices[[1]],
          matched_name = feature_names[[matching_indices[[1]]]],
          match_type = match_type
        )
      )
    }

    if (length(matching_indices) > 1L) {
      stop(
        "Ambiguous feature match for gene ",
        target_gene,
        " (",
        target_ensembl_gene_id,
        ").",
        call. = FALSE
      )
    }
  }

  stop(
    "Gene was not found in the Seurat counts matrix: ",
    target_gene,
    " (",
    target_ensembl_gene_id,
    ").",
    call. = FALSE
  )
}


calculate_edger_group_colour_max <- function(
    positive_values,
    upper_colour_quantile
) {

  positive_values <- positive_values[
    is.finite(positive_values) & positive_values > 0
  ]

  if (length(positive_values) == 0L) {
    return(1)
  }

  colour_max <- as.numeric(
    stats::quantile(
      positive_values,
      probs = upper_colour_quantile,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )

  if (!is.finite(colour_max) || colour_max <= 0) {
    colour_max <- max(positive_values, na.rm = TRUE)
  }

  if (!is.finite(colour_max) || colour_max <= 0) {
    colour_max <- 1
  }

  colour_max
}


prepare_edger_group_gene_data <- function(
    context,
    target_gene,
    target_ensembl_gene_id,
    significant_clusters,
    normalization_scale_factor = 10000,
    upper_colour_quantile = 0.99
) {

  significant_clusters <- sort_edger_group_cluster_ids(
    significant_clusters
  )

  missing_clusters <- setdiff(
    significant_clusters,
    context$cluster_levels
  )

  if (length(missing_clusters) > 0L) {
    stop(
      "Significant cluster(s) absent from Seurat object: ",
      paste(missing_clusters, collapse = ", "),
      call. = FALSE
    )
  }

  gene_match <- match_edger_group_gene_feature(
    feature_names = rownames(context$counts_matrix),
    target_gene = target_gene,
    target_ensembl_gene_id = target_ensembl_gene_id
  )

  gene_counts <- as.numeric(
    context$counts_matrix[
      gene_match$matched_index,
      ,
      drop = TRUE
    ]
  )

  names(gene_counts) <- colnames(context$counts_matrix)

  sample_plot_data <- vector(
    mode = "list",
    length = length(context$sample_layout$sample_order)
  )

  names(sample_plot_data) <- context$sample_layout$sample_order

  sample_cluster_rows <- list()
  sample_summary_rows <- list()

  for (sample_id in context$sample_layout$sample_order) {

    sample_cache <- context$sample_cache[[sample_id]]
    cells <- sample_cache$cells

    raw_counts <- as.numeric(gene_counts[cells])
    total_umi <- as.numeric(context$total_umi_per_spot[cells])

    log_normalized_expression <- numeric(length(cells))
    valid_spots <- total_umi > 0

    log_normalized_expression[valid_spots] <- log1p(
      raw_counts[valid_spots] /
        total_umi[valid_spots] *
        normalization_scale_factor
    )

    plot_data <- sample_cache$coordinates
    plot_data$target_gene <- target_gene
    plot_data$ensembl_gene_id <- target_ensembl_gene_id
    plot_data$target_raw_count <- raw_counts
    plot_data$total_UMI <- total_umi
    plot_data$logNormalized_expression <- log_normalized_expression
    plot_data$is_significant_cluster <-
      plot_data$cluster_id %in% significant_clusters

    cluster_ids_in_sample <- sort_edger_group_cluster_ids(
      plot_data$cluster_id
    )

    cluster_summary_rows <- lapply(
      cluster_ids_in_sample,
      function(cluster_id_current) {

        cluster_indices <- which(
          plot_data$cluster_id == cluster_id_current
        )

        cluster_raw_count <- sum(
          plot_data$target_raw_count[cluster_indices]
        )

        cluster_total_umi <- sum(
          plot_data$total_UMI[cluster_indices]
        )

        cluster_expression <- if (cluster_total_umi > 0) {
          log1p(
            cluster_raw_count /
              cluster_total_umi *
              normalization_scale_factor
          )
        } else {
          0
        }

        data.frame(
          sample_ID = sample_id,
          group_std = plot_data$group_std[[1]],
          sex_std = plot_data$sex_std[[1]],
          group_key = plot_data$group_key[[1]],
          cluster_id = cluster_id_current,
          is_significant_cluster =
            cluster_id_current %in% significant_clusters,
          number_of_spots = length(cluster_indices),
          positive_spots = sum(
            plot_data$target_raw_count[cluster_indices] > 0
          ),
          percent_positive_spots = 100 * mean(
            plot_data$target_raw_count[cluster_indices] > 0
          ),
          total_target_raw_count = cluster_raw_count,
          total_UMI = cluster_total_umi,
          cluster_logNormalized_expression = cluster_expression,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    )

    cluster_summary_table <- do.call(rbind, cluster_summary_rows)

    plot_data <- merge(
      x = plot_data,
      y = cluster_summary_table[
        ,
        c(
          "cluster_id",
          "cluster_logNormalized_expression"
        ),
        drop = FALSE
      ],
      by = "cluster_id",
      all.x = TRUE,
      sort = FALSE
    )

    plot_data <- plot_data[
      match(cells, plot_data$cell),
      ,
      drop = FALSE
    ]

    sample_plot_data[[sample_id]] <- plot_data

    sample_cluster_rows[[length(sample_cluster_rows) + 1L]] <-
      cluster_summary_table

    sample_summary_rows[[length(sample_summary_rows) + 1L]] <- data.frame(
      sample_ID = sample_id,
      group_std = plot_data$group_std[[1]],
      sex_std = plot_data$sex_std[[1]],
      group_key = plot_data$group_key[[1]],
      target_gene = target_gene,
      ensembl_gene_id = target_ensembl_gene_id,
      number_of_tissue_spots = nrow(plot_data),
      positive_tissue_spots = sum(plot_data$target_raw_count > 0),
      percent_positive_tissue_spots = 100 * mean(
        plot_data$target_raw_count > 0
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  sample_cluster_summary <- do.call(rbind, sample_cluster_rows)
  sample_summary <- do.call(rbind, sample_summary_rows)

  rownames(sample_cluster_summary) <- NULL
  rownames(sample_summary) <- NULL

  all_plot_data <- do.call(rbind, sample_plot_data)

  spot_colour_max <- calculate_edger_group_colour_max(
    all_plot_data$logNormalized_expression,
    upper_colour_quantile = upper_colour_quantile
  )

  cluster_colour_max <- calculate_edger_group_colour_max(
    sample_cluster_summary$cluster_logNormalized_expression,
    upper_colour_quantile = upper_colour_quantile
  )

  list(
    target_gene = target_gene,
    target_ensembl_gene_id = target_ensembl_gene_id,
    matched_feature = gene_match$matched_name,
    match_type = gene_match$match_type,
    significant_clusters = significant_clusters,
    sample_plot_data = sample_plot_data,
    sample_cluster_summary = sample_cluster_summary,
    sample_summary = sample_summary,
    spot_colour_max = spot_colour_max,
    cluster_colour_max = cluster_colour_max,
    normalization_scale_factor = normalization_scale_factor,
    upper_colour_quantile = upper_colour_quantile
  )
}


# ==============================================================================
# 6. Colour mapping, legends and four-column spatial layout
# ==============================================================================

map_edger_group_values_to_palette <- function(
    values,
    palette_colors,
    palette_values,
    limits,
    na_color = "#D9D9D9"
) {

  if (length(palette_colors) != length(palette_values)) {
    stop(
      "`palette_colors` and `palette_values` must have equal lengths.",
      call. = FALSE
    )
  }

  limits <- as.numeric(limits)

  if (
    length(limits) != 2L ||
      anyNA(limits) ||
      limits[[2]] <= limits[[1]]
  ) {
    stop(
      "`limits` must contain two increasing numeric values.",
      call. = FALSE
    )
  }

  normalized_values <- (values - limits[[1]]) / (limits[[2]] - limits[[1]])
  normalized_values <- pmin(pmax(normalized_values, 0), 1)

  palette_function <- scales::gradient_n_pal(
    colours = palette_colors,
    values = palette_values
  )

  output_colors <- palette_function(normalized_values)
  output_colors[is.na(values) | !is.finite(values)] <- na_color

  output_colors
}


extract_edger_group_plot_legend <- function(plot_object) {

  plot_gtable <- ggplot2::ggplotGrob(plot_object)

  guide_indices <- which(
    grepl("^guide-box", plot_gtable$layout$name)
  )

  non_empty_guide_indices <- guide_indices[
    !vapply(
      plot_gtable$grobs[guide_indices],
      inherits,
      logical(1),
      what = "zeroGrob"
    )
  ]

  if (length(non_empty_guide_indices) == 0L) {
    stop("Could not extract a non-empty legend.", call. = FALSE)
  }

  plot_gtable$grobs[[non_empty_guide_indices[[1]]]]
}


create_edger_group_colourbar_legend <- function(
    palette_colors,
    palette_values,
    colour_max,
    legend_title,
    legend_bar_width_mm = 40,
    legend_bar_height_mm = 4.5
) {

  if (!is.finite(colour_max) || colour_max <= 0) {
    colour_max <- 1
  }

  legend_data <- data.frame(
    x = seq(0, colour_max, length.out = 100),
    y = 1
  )

  legend_source <- ggplot2::ggplot(
    legend_data,
    ggplot2::aes(x = .data$x, y = .data$y, colour = .data$x)
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_gradientn(
      colours = palette_colors,
      values = palette_values,
      limits = c(0, colour_max),
      oob = scales::squish
    ) +
    ggplot2::labs(colour = legend_title) +
    ggplot2::theme_void(base_family = "DejaVu Sans") +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.title = ggplot2::element_text(
        size = 9,
        face = "bold",
        hjust = 0.5
      ),
      legend.text = ggplot2::element_text(size = 8),
      legend.margin = ggplot2::margin(t = 1, r = 4, b = 2, l = 4)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_colourbar(
        title.position = "top",
        direction = "horizontal",
        barwidth = grid::unit(legend_bar_width_mm, "mm"),
        barheight = grid::unit(legend_bar_height_mm, "mm")
      )
    )

  patchwork::wrap_elements(
    full = extract_edger_group_plot_legend(legend_source)
  )
}


create_edger_group_spatial_panel <- function(
    sample_cache,
    plot_data,
    colour_column,
    target_gene,
    show_histology_image = FALSE,
    point_size_no_image = 1.1,
    point_size_with_image = 0.80,
    panel_border_linewidth = 0.8
) {

  title_text <- paste0("Sample: ", sample_cache$sample_ID)

  subtitle_text <- paste0(
    plot_data$sex_std[[1]],
    " | ",
    plot_data$group_std[[1]],
    "\n",
    target_gene,
    "+ spots: ",
    sum(plot_data$target_raw_count > 0),
    "/",
    nrow(plot_data),
    " (",
    formatC(
      100 * mean(plot_data$target_raw_count > 0),
      format = "f",
      digits = 1
    ),
    "%)"
  )

  if (
    isTRUE(show_histology_image) &&
      !is.null(sample_cache$image_array)
  ) {

    output_plot <- ggplot2::ggplot() +
      ggplot2::annotation_raster(
        raster = sample_cache$image_array,
        xmin = 0,
        xmax = sample_cache$image_width,
        ymin = sample_cache$image_height,
        ymax = 0
      ) +
      ggplot2::geom_point(
        data = plot_data,
        ggplot2::aes(
          x = .data$x_plot,
          y = .data$y_plot,
          colour = .data[[colour_column]]
        ),
        size = point_size_with_image,
        shape = 16,
        stroke = 0
      )

  } else {

    output_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data$x_plot,
        y = .data$y_plot,
        colour = .data[[colour_column]]
      )
    ) +
      ggplot2::geom_point(
        size = point_size_no_image,
        shape = 16,
        stroke = 0
      )
  }

  output_plot +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(
      limits = sample_cache$frame_limits$x_limits,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_reverse(
      limits = rev(sample_cache$frame_limits$y_limits),
      expand = c(0, 0)
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = title_text,
      subtitle = subtitle_text
    ) +
    ggplot2::theme_void(base_family = "DejaVu Sans") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = panel_border_linewidth
      ),
      plot.title = ggplot2::element_text(
        size = 12.5,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 4)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 7.5,
        hjust = 0.5,
        lineheight = 1.05,
        margin = ggplot2::margin(b = 6)
      ),
      plot.margin = ggplot2::margin(t = 8, r = 8, b = 8, l = 8)
    )
}


arrange_edger_group_spatial_plots_four_columns <- function(
    plot_list_named,
    sample_layout,
    legend_panel,
    plot_title,
    plot_subtitle,
    legend_height_ratio = 0.075
) {

  plots_by_column <- lapply(
    sample_layout$column_order,
    function(column_key) {

      sample_ids <- sample_layout$samples_by_column[[column_key]]

      current_plots <- lapply(
        sample_ids,
        function(sample_id) {
          plot_list_named[[sample_id]]
        }
      )

      if (length(current_plots) < sample_layout$max_rows) {
        current_plots <- c(
          current_plots,
          rep(
            list(patchwork::plot_spacer()),
            sample_layout$max_rows - length(current_plots)
          )
        )
      }

      current_plots
    }
  )

  names(plots_by_column) <- sample_layout$column_order

  interleaved_plots <- list()

  for (row_index in seq_len(sample_layout$max_rows)) {
    for (column_key in sample_layout$column_order) {
      interleaved_plots[[length(interleaved_plots) + 1L]] <-
        plots_by_column[[column_key]][[row_index]]
    }
  }

  spatial_grid <- patchwork::wrap_plots(
    interleaved_plots,
    ncol = 4L
  )

  (
    legend_panel /
      spatial_grid
  ) +
    patchwork::plot_layout(
      heights = c(legend_height_ratio, 1)
    ) +
    patchwork::plot_annotation(
      title = plot_title,
      subtitle = plot_subtitle,
      theme = ggplot2::theme(
        text = ggplot2::element_text(family = "DejaVu Sans"),
        plot.title = ggplot2::element_text(
          size = 20,
          face = "bold",
          hjust = 0.5
        ),
        plot.subtitle = ggplot2::element_text(
          size = 9.5,
          hjust = 0.5,
          lineheight = 1.05,
          margin = ggplot2::margin(b = 4)
        ),
        plot.margin = ggplot2::margin(t = 6, r = 8, b = 8, l = 8)
      )
    )
}


plot_edger_group_gene_spatial <- function(
    context,
    gene_data,
    expression_level = c("spot", "cluster"),
    highlight_significant_clusters = FALSE,
    red_palette_colors,
    green_palette_colors,
    palette_values,
    selection_label,
    show_histology_image = FALSE,
    point_size_no_image = 1.1,
    point_size_with_image = 0.80
) {

  expression_level <- match.arg(expression_level)

  if (expression_level == "spot") {
    value_column <- "logNormalized_expression"
    colour_max <- gene_data$spot_colour_max
    expression_description <- "per-spot expression"
    aggregation_description <- paste0(
      "log1p(count / total UMI × ",
      format(
        gene_data$normalization_scale_factor,
        big.mark = ",",
        scientific = FALSE
      ),
      ")"
    )
  } else {
    value_column <- "cluster_logNormalized_expression"
    colour_max <- gene_data$cluster_colour_max
    expression_description <- "sample × cluster aggregated expression"
    aggregation_description <- paste0(
      "log1p(sum count / sum UMI × ",
      format(
        gene_data$normalization_scale_factor,
        big.mark = ",",
        scientific = FALSE
      ),
      ")"
    )
  }

  colour_limits <- c(0, colour_max)

  panels <- lapply(
    context$sample_layout$sample_order,
    function(sample_id) {

      plot_data <- gene_data$sample_plot_data[[sample_id]]
      expression_values <- plot_data[[value_column]]

      red_colours <- map_edger_group_values_to_palette(
        values = expression_values,
        palette_colors = red_palette_colors,
        palette_values = palette_values,
        limits = colour_limits
      )

      if (isTRUE(highlight_significant_clusters)) {

        green_colours <- map_edger_group_values_to_palette(
          values = expression_values,
          palette_colors = green_palette_colors,
          palette_values = palette_values,
          limits = colour_limits
        )

        plot_data$plot_colour <- ifelse(
          plot_data$is_significant_cluster,
          red_colours,
          green_colours
        )

      } else {
        plot_data$plot_colour <- red_colours
      }

      create_edger_group_spatial_panel(
        sample_cache = context$sample_cache[[sample_id]],
        plot_data = plot_data,
        colour_column = "plot_colour",
        target_gene = gene_data$target_gene,
        show_histology_image = show_histology_image,
        point_size_no_image = point_size_no_image,
        point_size_with_image = point_size_with_image
      )
    }
  )

  names(panels) <- context$sample_layout$sample_order

  significant_cluster_text <- paste(
    paste0("C", gene_data$significant_clusters),
    collapse = ", "
  )

  if (isTRUE(highlight_significant_clusters)) {

    red_legend <- create_edger_group_colourbar_legend(
      palette_colors = red_palette_colors,
      palette_values = palette_values,
      colour_max = colour_max,
      legend_title = paste0(
        "Statistically significant cluster(s) ",
        significant_cluster_text,
        " — ",
        expression_description
      ),
      legend_bar_width_mm = 34
    )

    green_legend <- create_edger_group_colourbar_legend(
      palette_colors = green_palette_colors,
      palette_values = palette_values,
      colour_max = colour_max,
      legend_title = paste0(
        "Other clusters — ",
        expression_description
      ),
      legend_bar_width_mm = 34
    )

    legend_panel <- red_legend +
      green_legend +
      patchwork::plot_layout(ncol = 2)

    plot_title <- paste0(
      gene_data$target_gene,
      " — ",
      expression_description,
      ": significant clusters red, other clusters green"
    )

  } else {

    legend_panel <- create_edger_group_colourbar_legend(
      palette_colors = red_palette_colors,
      palette_values = palette_values,
      colour_max = colour_max,
      legend_title = paste0(
        gene_data$target_gene,
        " — ",
        expression_description,
        "\n",
        aggregation_description
      )
    )

    plot_title <- paste0(
      gene_data$target_gene,
      " — ",
      expression_description,
      " using one red scale"
    )
  }

  plot_subtitle <- paste0(
    "Statistically significant main-group cluster(s): ",
    significant_cluster_text,
    " | selection: ",
    selection_label,
    "\nColumns: Male neurotypical | Male ASD | Female neurotypical | Female ASD",
    " | common colour range: 0–",
    format_edger_group_number(colour_max),
    " | upper limit = ",
    gene_data$upper_colour_quantile * 100,
    "th percentile of positive values"
  )

  arrange_edger_group_spatial_plots_four_columns(
    plot_list_named = panels,
    sample_layout = context$sample_layout,
    legend_panel = legend_panel,
    plot_title = plot_title,
    plot_subtitle = plot_subtitle
  )
}


# ==============================================================================
# 7. edgeR statistics for one gene and its significant clusters
# ==============================================================================

extract_gene_cluster_edger_statistics <- function(
    combined_results,
    target_ensembl_gene_id,
    significant_clusters
) {

  significant_clusters <- as.character(significant_clusters)

  main_group_table <- combined_results[[
    "Overall_Group_ASD_vs_Neurotypical"
  ]] |>
    dplyr::filter(
      .data$ensembl_gene_id == target_ensembl_gene_id,
      as.character(.data$cluster_id) %in% significant_clusters
    ) |>
    dplyr::transmute(
      cluster_id = as.character(.data$cluster_id),
      ensembl_gene_id = .data$ensembl_gene_id,
      gene = .data$gene,
      main_group_log2FC = .data$logFC,
      main_group_PValue = .data$PValue,
      main_group_FDR = .data$FDR
    )

  male_group_table <- combined_results[[
    "ASD_Male_vs_Neurotypical_Male"
  ]] |>
    dplyr::filter(
      .data$ensembl_gene_id == target_ensembl_gene_id,
      as.character(.data$cluster_id) %in% significant_clusters
    ) |>
    dplyr::transmute(
      cluster_id = as.character(.data$cluster_id),
      male_group_log2FC = .data$logFC,
      male_group_PValue = .data$PValue,
      male_group_FDR = .data$FDR
    )

  female_group_table <- combined_results[[
    "ASD_Female_vs_Neurotypical_Female"
  ]] |>
    dplyr::filter(
      .data$ensembl_gene_id == target_ensembl_gene_id,
      as.character(.data$cluster_id) %in% significant_clusters
    ) |>
    dplyr::transmute(
      cluster_id = as.character(.data$cluster_id),
      female_group_log2FC = .data$logFC,
      female_group_PValue = .data$PValue,
      female_group_FDR = .data$FDR
    )

  statistics_table <- main_group_table |>
    dplyr::left_join(male_group_table, by = "cluster_id") |>
    dplyr::left_join(female_group_table, by = "cluster_id") |>
    dplyr::mutate(
      cluster_numeric = suppressWarnings(as.numeric(.data$cluster_id))
    ) |>
    dplyr::arrange(.data$cluster_numeric, .data$cluster_id) |>
    dplyr::select(-cluster_numeric)

  missing_clusters <- setdiff(
    significant_clusters,
    statistics_table$cluster_id
  )

  if (length(missing_clusters) > 0L) {
    stop(
      "Missing edgeR statistics for cluster(s): ",
      paste(missing_clusters, collapse = ", "),
      call. = FALSE
    )
  }

  statistics_table
}


# ==============================================================================
# 8. Four-group barplot for one statistically significant cluster
# ==============================================================================

plot_edger_group_barplot_for_cluster <- function(
    gene_data,
    statistics_row,
    min_spots_per_sample_cluster = 20L,
    group_colors = c(
      "male_neurotypical" = "#4C78A8",
      "male_asd" = "#C44E52",
      "female_neurotypical" = "#8AB8DB",
      "female_asd" = "#E78A8A"
    ),
    point_size = 3.0,
    jitter_width = 0.07,
    bar_width = 0.70,
    errorbar_width = 0.17
) {

  cluster_id_current <- as.character(statistics_row$cluster_id[[1]])

  group_order <- c(
    "male_neurotypical",
    "male_asd",
    "female_neurotypical",
    "female_asd"
  )

  cluster_data <- gene_data$sample_cluster_summary |>
    dplyr::filter(
      as.character(.data$cluster_id) == cluster_id_current,
      .data$number_of_spots >= min_spots_per_sample_cluster
    ) |>
    dplyr::mutate(
      group_key = factor(.data$group_key, levels = group_order)
    )

  if (nrow(cluster_data) == 0L) {
    stop(
      "No sample × cluster values passed the spot threshold for cluster ",
      cluster_id_current,
      ".",
      call. = FALSE
    )
  }

  missing_groups <- setdiff(
    group_order,
    unique(as.character(cluster_data$group_key))
  )

  if (length(missing_groups) > 0L) {
    stop(
      "Cluster ",
      cluster_id_current,
      " is missing barplot group(s) after the spot threshold: ",
      paste(missing_groups, collapse = ", "),
      call. = FALSE
    )
  }

  group_summary <- cluster_data |>
    dplyr::group_by(.data$group_key) |>
    dplyr::summarise(
      number_of_samples = dplyr::n(),
      total_spots = sum(.data$number_of_spots, na.rm = TRUE),
      positive_spots = sum(.data$positive_spots, na.rm = TRUE),
      pooled_percent_positive = 100 *
        sum(.data$positive_spots, na.rm = TRUE) /
        sum(.data$number_of_spots, na.rm = TRUE),
      mean_expression = mean(
        .data$cluster_logNormalized_expression,
        na.rm = TRUE
      ),
      sd_expression = stats::sd(
        .data$cluster_logNormalized_expression,
        na.rm = TRUE
      ),
      se_expression = stats::sd(
        .data$cluster_logNormalized_expression,
        na.rm = TRUE
      ) / sqrt(dplyr::n()),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      group_key = factor(.data$group_key, levels = group_order),
      ymin = pmax(0, .data$mean_expression - .data$se_expression),
      ymax = .data$mean_expression + .data$se_expression
    ) |>
    dplyr::arrange(.data$group_key)

  if (any(!is.finite(group_summary$se_expression))) {
    group_summary$se_expression[
      !is.finite(group_summary$se_expression)
    ] <- 0
    group_summary$ymin <- pmax(
      0,
      group_summary$mean_expression - group_summary$se_expression
    )
    group_summary$ymax <-
      group_summary$mean_expression + group_summary$se_expression
  }

  group_label_names <- c(
    "male_neurotypical" = "Male\nNeurotypical",
    "male_asd" = "Male\nASD",
    "female_neurotypical" = "Female\nNeurotypical",
    "female_asd" = "Female\nASD"
  )

  percent_lookup <- setNames(
    group_summary$pooled_percent_positive,
    as.character(group_summary$group_key)
  )

  x_axis_labels <- vapply(
    group_order,
    function(group_key_current) {
      paste0(
        group_label_names[[group_key_current]],
        "\n",
        formatC(
          percent_lookup[[group_key_current]],
          format = "f",
          digits = 1
        ),
        "% spots+"
      )
    },
    character(1)
  )

  names(x_axis_labels) <- group_order

  observed_max <- max(
    c(
      cluster_data$cluster_logNormalized_expression,
      group_summary$ymax
    ),
    na.rm = TRUE
  )

  if (!is.finite(observed_max) || observed_max <= 0) {
    observed_max <- 1
  }

  bracket_height <- observed_max * 0.06
  male_bracket_y <- observed_max * 1.16
  female_bracket_y <- observed_max * 1.34
  plot_upper_limit <- observed_max * 1.55

  male_label <- paste0(
    "P=",
    format_edger_group_probability(statistics_row$male_group_PValue),
    " | FDR=",
    format_edger_group_probability(statistics_row$male_group_FDR)
  )

  female_label <- paste0(
    "P=",
    format_edger_group_probability(statistics_row$female_group_PValue),
    " | FDR=",
    format_edger_group_probability(statistics_row$female_group_FDR)
  )

  main_title <- paste0(
    gene_data$target_gene,
    " | cluster C",
    cluster_id_current,
    " | main group effect (ASD vs Neurotypical)"
  )

  main_subtitle <- paste0(
    "log2FC=",
    format_edger_group_number(statistics_row$main_group_log2FC),
    " | nominal P=",
    format_edger_group_probability(statistics_row$main_group_PValue),
    " | FDR=",
    format_edger_group_probability(statistics_row$main_group_FDR),
    "\nBars: mean | whiskers: mean ± SE | dots: individual sample × cluster values",
    " | included samples have ≥",
    min_spots_per_sample_cluster,
    " spots in this cluster"
  )

  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = group_summary,
      ggplot2::aes(
        x = .data$group_key,
        y = .data$mean_expression
      ),
      width = bar_width,
      fill = "white",
      colour = "black",
      linewidth = 1.50
    ) +
    ggplot2::geom_errorbar(
      data = group_summary,
      ggplot2::aes(
        x = .data$group_key,
        ymin = .data$ymin,
        ymax = .data$ymax
      ),
      width = errorbar_width,
      colour = "black",
      linewidth = 0.90
    ) +
    ggplot2::geom_point(
      data = cluster_data,
      ggplot2::aes(
        x = .data$group_key,
        y = .data$cluster_logNormalized_expression
      ),
      position = ggplot2::position_jitter(
        width = jitter_width,
        height = 0,
        seed = 123
      ),
      shape = 21,
      size = point_size,
      stroke = 0.95,
      fill = "white",
      colour = "black"
    ) +
    ggplot2::annotate(
      "segment",
      x = 1,
      xend = 2,
      y = male_bracket_y,
      yend = male_bracket_y,
      linewidth = 0.8
    ) +
    ggplot2::annotate(
      "segment",
      x = 1,
      xend = 1,
      y = male_bracket_y - bracket_height,
      yend = male_bracket_y,
      linewidth = 0.8
    ) +
    ggplot2::annotate(
      "segment",
      x = 2,
      xend = 2,
      y = male_bracket_y - bracket_height,
      yend = male_bracket_y,
      linewidth = 0.8
    ) +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = male_bracket_y + bracket_height * 0.35,
      label = male_label,
      size = 3.5,
      family = "DejaVu Sans"
    ) +
    ggplot2::annotate(
      "segment",
      x = 3,
      xend = 4,
      y = female_bracket_y,
      yend = female_bracket_y,
      linewidth = 0.8
    ) +
    ggplot2::annotate(
      "segment",
      x = 3,
      xend = 3,
      y = female_bracket_y - bracket_height,
      yend = female_bracket_y,
      linewidth = 0.8
    ) +
    ggplot2::annotate(
      "segment",
      x = 4,
      xend = 4,
      y = female_bracket_y - bracket_height,
      yend = female_bracket_y,
      linewidth = 0.8
    ) +
    ggplot2::annotate(
      "text",
      x = 3.5,
      y = female_bracket_y + bracket_height * 0.35,
      label = female_label,
      size = 3.5,
      family = "DejaVu Sans"
    ) +
    ggplot2::scale_x_discrete(
      limits = group_order,
      labels = x_axis_labels,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, plot_upper_limit),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::labs(
      title = main_title,
      subtitle = main_subtitle,
      x = NULL,
      y = paste0(
        gene_data$target_gene,
        " sample × cluster expression\n",
        "log1p(sum count / sum UMI × ",
        format(
          gene_data$normalization_scale_factor,
          big.mark = ",",
          scientific = FALSE
        ),
        ")"
      )
    ) +
    ggplot2::theme_classic(base_family = "DejaVu Sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        size = 9.5,
        hjust = 0.5,
        lineheight = 1.08,
        margin = ggplot2::margin(b = 12)
      ),
      axis.title.y = ggplot2::element_text(
        size = 11,
        face = "bold",
        margin = ggplot2::margin(r = 10)
      ),
      axis.text.x = ggplot2::element_text(
        size = 9,
        face = "bold",
        lineheight = 0.95,
        margin = ggplot2::margin(t = 8)
      ),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.line = ggplot2::element_line(linewidth = 0.7),
      axis.ticks = ggplot2::element_line(linewidth = 0.7),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 10, r = 18, b = 12, l = 18)
    )
}


# ==============================================================================
# 9. Output paths for one gene
# ==============================================================================

get_edger_group_gene_output_files <- function(
    output_root_directory,
    target_gene,
    significant_clusters
) {

  significant_cluster_label <- build_significant_cluster_label(
    significant_clusters
  )

  gene_file_label <- sanitize_edger_group_file_component(target_gene)

  gene_output_directory <- file.path(
    output_root_directory,
    paste0(
      gene_file_label,
      "_",
      significant_cluster_label
    )
  )

  file_prefix <- paste0(
    gene_file_label,
    "_",
    significant_cluster_label
  )

  barplot_directory <- file.path(
    gene_output_directory,
    "05_barplots_perSignificantCluster"
  )

  list(
    gene_output_directory = gene_output_directory,
    barplot_directory = barplot_directory,
    per_spot_red_pdf = file.path(
      gene_output_directory,
      paste0("01_", file_prefix, "_perSpotExpression_redScale.pdf")
    ),
    cluster_aggregated_red_pdf = file.path(
      gene_output_directory,
      paste0(
        "02_",
        file_prefix,
        "_clusterAggregatedExpression_redScale.pdf"
      )
    ),
    per_spot_red_green_pdf = file.path(
      gene_output_directory,
      paste0(
        "03_",
        file_prefix,
        "_perSpotExpression_significantClustersRed_otherClustersGreen.pdf"
      )
    ),
    cluster_aggregated_red_green_pdf = file.path(
      gene_output_directory,
      paste0(
        "04_",
        file_prefix,
        "_clusterAggregatedExpression_significantClustersRed_otherClustersGreen.pdf"
      )
    ),
    edger_statistics_tsv = file.path(
      gene_output_directory,
      paste0(
        "06_",
        file_prefix,
        "_edgeRStatistics_significantClusters.tsv"
      )
    ),
    sample_cluster_summary_tsv = file.path(
      gene_output_directory,
      paste0(
        "07_",
        file_prefix,
        "_sampleClusterExpressionSummary.tsv"
      )
    ),
    sample_summary_tsv = file.path(
      gene_output_directory,
      paste0(
        "08_",
        file_prefix,
        "_samplePercentPositiveSpots.tsv"
      )
    ),
    barplot_style_marker = file.path(
      gene_output_directory,
      "09_barplotStyle_whiteFill_thickBlackOutline_v1.txt"
    )
  )
}


edger_group_gene_outputs_complete <- function(
    output_files,
    significant_clusters,
    target_gene
) {

  fixed_expected_files <- unlist(
    output_files[
      setdiff(
        names(output_files),
        c("gene_output_directory", "barplot_directory")
      )
    ],
    use.names = FALSE
  )

  barplot_files <- file.path(
    output_files$barplot_directory,
    paste0(
      "05_",
      sanitize_edger_group_file_component(target_gene),
      "_cluster_C",
      sanitize_edger_group_file_component(significant_clusters),
      "_fourGroups_withSexSpecificEdgeRStatistics.pdf"
    )
  )

  all(file.exists(c(fixed_expected_files, barplot_files)))
}


# ==============================================================================
# 10. Save all outputs for one gene
# ==============================================================================

save_edger_group_gene_visualizations <- function(
    context,
    combined_results,
    target_gene,
    target_ensembl_gene_id,
    significant_clusters,
    output_root_directory,
    selection_label,
    min_spots_per_sample_cluster = 20L,
    normalization_scale_factor = 10000,
    upper_colour_quantile = 0.99,
    red_palette_colors = c(
      "#D9D9D9",
      "#FEE5D9",
      "#FCAE91",
      "#FB6A4A",
      "#DE2D26",
      "#A50F15",
      "#67000D"
    ),
    green_palette_colors = c(
      "#D9D9D9",
      "#E5F5E0",
      "#A1D99B",
      "#74C476",
      "#31A354",
      "#006D2C",
      "#00441B"
    ),
    palette_values = c(
      0.000,
      0.006,
      0.060,
      0.180,
      0.400,
      0.700,
      1.000
    ),
    group_colors = c(
      "male_neurotypical" = "#4C78A8",
      "male_asd" = "#C44E52",
      "female_neurotypical" = "#8AB8DB",
      "female_asd" = "#E78A8A"
    ),
    show_histology_image = FALSE,
    point_size_no_image = 1.1,
    point_size_with_image = 0.80,
    spatial_pdf_width = 18,
    spatial_pdf_height = 23,
    barplot_pdf_width = 11,
    barplot_pdf_height = 8,
    skip_completed = TRUE
) {

  significant_clusters <- sort_edger_group_cluster_ids(
    significant_clusters
  )

  output_files <- get_edger_group_gene_output_files(
    output_root_directory = output_root_directory,
    target_gene = target_gene,
    significant_clusters = significant_clusters
  )

  if (
    isTRUE(skip_completed) &&
      edger_group_gene_outputs_complete(
        output_files = output_files,
        significant_clusters = significant_clusters,
        target_gene = target_gene
      )
  ) {
    message("Skipping completed gene: ", target_gene)
    return(invisible(output_files))
  }

  dir.create(
    output_files$gene_output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    output_files$barplot_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  gene_data <- prepare_edger_group_gene_data(
    context = context,
    target_gene = target_gene,
    target_ensembl_gene_id = target_ensembl_gene_id,
    significant_clusters = significant_clusters,
    normalization_scale_factor = normalization_scale_factor,
    upper_colour_quantile = upper_colour_quantile
  )

  statistics_table <- extract_gene_cluster_edger_statistics(
    combined_results = combined_results,
    target_ensembl_gene_id = target_ensembl_gene_id,
    significant_clusters = significant_clusters
  )

  write_edger_group_visualization_tsv(
    statistics_table,
    output_files$edger_statistics_tsv
  )

  write_edger_group_visualization_tsv(
    gene_data$sample_cluster_summary,
    output_files$sample_cluster_summary_tsv
  )

  write_edger_group_visualization_tsv(
    gene_data$sample_summary,
    output_files$sample_summary_tsv
  )

  plot_1 <- plot_edger_group_gene_spatial(
    context = context,
    gene_data = gene_data,
    expression_level = "spot",
    highlight_significant_clusters = FALSE,
    red_palette_colors = red_palette_colors,
    green_palette_colors = green_palette_colors,
    palette_values = palette_values,
    selection_label = selection_label,
    show_histology_image = show_histology_image,
    point_size_no_image = point_size_no_image,
    point_size_with_image = point_size_with_image
  )

  plot_2 <- plot_edger_group_gene_spatial(
    context = context,
    gene_data = gene_data,
    expression_level = "cluster",
    highlight_significant_clusters = FALSE,
    red_palette_colors = red_palette_colors,
    green_palette_colors = green_palette_colors,
    palette_values = palette_values,
    selection_label = selection_label,
    show_histology_image = show_histology_image,
    point_size_no_image = point_size_no_image,
    point_size_with_image = point_size_with_image
  )

  plot_3 <- plot_edger_group_gene_spatial(
    context = context,
    gene_data = gene_data,
    expression_level = "spot",
    highlight_significant_clusters = TRUE,
    red_palette_colors = red_palette_colors,
    green_palette_colors = green_palette_colors,
    palette_values = palette_values,
    selection_label = selection_label,
    show_histology_image = show_histology_image,
    point_size_no_image = point_size_no_image,
    point_size_with_image = point_size_with_image
  )

  plot_4 <- plot_edger_group_gene_spatial(
    context = context,
    gene_data = gene_data,
    expression_level = "cluster",
    highlight_significant_clusters = TRUE,
    red_palette_colors = red_palette_colors,
    green_palette_colors = green_palette_colors,
    palette_values = palette_values,
    selection_label = selection_label,
    show_histology_image = show_histology_image,
    point_size_no_image = point_size_no_image,
    point_size_with_image = point_size_with_image
  )

  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  ggplot2::ggsave(
    filename = output_files$per_spot_red_pdf,
    plot = plot_1,
    device = pdf_device,
    width = spatial_pdf_width,
    height = spatial_pdf_height,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = output_files$cluster_aggregated_red_pdf,
    plot = plot_2,
    device = pdf_device,
    width = spatial_pdf_width,
    height = spatial_pdf_height,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = output_files$per_spot_red_green_pdf,
    plot = plot_3,
    device = pdf_device,
    width = spatial_pdf_width,
    height = spatial_pdf_height,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = output_files$cluster_aggregated_red_green_pdf,
    plot = plot_4,
    device = pdf_device,
    width = spatial_pdf_width,
    height = spatial_pdf_height,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )

  for (cluster_id_current in significant_clusters) {

    statistics_row <- statistics_table |>
      dplyr::filter(.data$cluster_id == cluster_id_current)

    if (nrow(statistics_row) != 1L) {
      stop(
        "Expected exactly one statistics row for cluster ",
        cluster_id_current,
        ".",
        call. = FALSE
      )
    }

    barplot <- plot_edger_group_barplot_for_cluster(
      gene_data = gene_data,
      statistics_row = statistics_row,
      min_spots_per_sample_cluster = min_spots_per_sample_cluster,
      group_colors = group_colors
    )

    barplot_file <- file.path(
      output_files$barplot_directory,
      paste0(
        "05_",
        sanitize_edger_group_file_component(target_gene),
        "_cluster_C",
        sanitize_edger_group_file_component(cluster_id_current),
        "_fourGroups_withSexSpecificEdgeRStatistics.pdf"
      )
    )

    ggplot2::ggsave(
      filename = barplot_file,
      plot = barplot,
      device = pdf_device,
      width = barplot_pdf_width,
      height = barplot_pdf_height,
      units = "in",
      limitsize = FALSE,
      bg = "white"
    )
  }

  writeLines(
    text = c(
      "barplot_style=white_fill_thick_black_outline",
      "bar_fill=white",
      "bar_border=black",
      "bar_border_linewidth=1.50",
      "point_fill=white",
      "point_border=black"
    ),
    con = output_files$barplot_style_marker,
    useBytes = TRUE
  )

  if (
    !edger_group_gene_outputs_complete(
      output_files = output_files,
      significant_clusters = significant_clusters,
      target_gene = target_gene
    )
  ) {
    stop(
      "Not all expected outputs were created for gene: ",
      target_gene,
      call. = FALSE
    )
  }

  invisible(output_files)
}


# ==============================================================================
# 11. Complete main-effect-group visualization workflow
# ==============================================================================

run_main_effect_group_gene_visualizations <- function(
    input_seurat_rdata_file,
    edger_results_rdata_file,
    output_root_directory,
    cluster_column = "leiden_res0.40",
    assay_name = "RNA",
    sample_id_column = "sample_ID",
    group_column = "fmt_donor_group",
    sex_column = "sex",
    included_sample_ids = NULL,
    seurat_object_name = NULL,
    fdr_threshold = 0.05,
    abs_log2fc_threshold = 0.7,
    maximum_number_of_genes = 3L,
    manual_target_genes = NULL,
    min_spots_per_sample_cluster = 20L,
    normalization_scale_factor = 10000,
    upper_colour_quantile = 0.99,
    red_palette_colors = c(
      "#D9D9D9",
      "#FEE5D9",
      "#FCAE91",
      "#FB6A4A",
      "#DE2D26",
      "#A50F15",
      "#67000D"
    ),
    green_palette_colors = c(
      "#D9D9D9",
      "#E5F5E0",
      "#A1D99B",
      "#74C476",
      "#31A354",
      "#006D2C",
      "#00441B"
    ),
    palette_values = c(
      0.000,
      0.006,
      0.060,
      0.180,
      0.400,
      0.700,
      1.000
    ),
    group_colors = c(
      "male_neurotypical" = "#4C78A8",
      "male_asd" = "#C44E52",
      "female_neurotypical" = "#8AB8DB",
      "female_asd" = "#E78A8A"
    ),
    show_histology_image = FALSE,
    point_size_no_image = 1.1,
    point_size_with_image = 0.80,
    spatial_pdf_width = 18,
    spatial_pdf_height = 23,
    barplot_pdf_width = 11,
    barplot_pdf_height = 8,
    skip_completed_genes = TRUE,
    continue_after_gene_error = TRUE
) {

  check_required_packages_edger_group_visualization()

  dir.create(
    output_root_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  selection_label <- paste0(
    "FDR < ",
    fdr_threshold,
    " and |log2FC| > ",
    abs_log2fc_threshold
  )

  selected_genes_file <- file.path(
    output_root_directory,
    "00_selectedGenes_mainEffectGroup.tsv"
  )

  selected_gene_cluster_statistics_file <- file.path(
    output_root_directory,
    "00_selectedGeneClusterStatistics_mainEffectGroup.tsv"
  )

  run_status_file <- file.path(
    output_root_directory,
    "00_runStatus_mainEffectGroup.tsv"
  )

  message("Loading edgeR results: ", edger_results_rdata_file)

  edger_objects <- load_edger_results_for_group_visualization(
    edger_results_rdata_file
  )

  selection_objects <- select_main_group_genes_for_visualization(
    combined_results = edger_objects$combined_results,
    fdr_threshold = fdr_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold,
    manual_target_genes = manual_target_genes,
    maximum_number_of_genes = maximum_number_of_genes
  )

  selected_genes_output <- flatten_selected_gene_table(
    gene_selection_table = selection_objects$gene_selection_table,
    fdr_threshold = fdr_threshold,
    abs_log2fc_threshold = abs_log2fc_threshold
  )

  write_edger_group_visualization_tsv(
    selected_genes_output,
    selected_genes_file
  )

  write_edger_group_visualization_tsv(
    selection_objects$selected_significant_results,
    selected_gene_cluster_statistics_file
  )

  message("Selected genes: ", nrow(selection_objects$gene_selection_table))
  message(
    paste0(
      "  ",
      selection_objects$gene_selection_table$gene,
      " | significant cluster(s): ",
      vapply(
        selection_objects$gene_selection_table$significant_clusters,
        paste,
        character(1),
        collapse = ","
      ),
      collapse = "\n"
    )
  )

  message("Loading Seurat object: ", input_seurat_rdata_file)

  loaded_seurat <- load_single_seurat_for_edger_group_visualization(
    input_rdata_file = input_seurat_rdata_file,
    requested_object_name = seurat_object_name
  )

  message("Loaded Seurat object: ", loaded_seurat$object_name)
  message(
    "Object dimensions: ",
    nrow(loaded_seurat$object),
    " genes x ",
    ncol(loaded_seurat$object),
    " spots"
  )

  message("Preparing shared spatial/count context...")

  visualization_context <- prepare_edger_group_visualization_context(
    seurat_object = loaded_seurat$object,
    cluster_column = cluster_column,
    assay_name = assay_name,
    sample_id_column = sample_id_column,
    group_column = group_column,
    sex_column = sex_column,
    included_sample_ids = included_sample_ids,
    image_scale = "lowres",
    panel_padding_fraction = 0.03
  )

  rm(loaded_seurat)
  invisible(gc(verbose = FALSE))

  status_rows <- list()
  run_started_at <- Sys.time()

  for (gene_index in seq_len(nrow(selection_objects$gene_selection_table))) {

    target_gene <- selection_objects$gene_selection_table$gene[[gene_index]]
    target_ensembl_gene_id <-
      selection_objects$gene_selection_table$ensembl_gene_id[[gene_index]]
    significant_clusters <-
      selection_objects$gene_selection_table$significant_clusters[[gene_index]]

    gene_started_at <- Sys.time()

    message(
      "\n[",
      gene_index,
      "/",
      nrow(selection_objects$gene_selection_table),
      "] Processing ",
      target_gene,
      " | significant cluster(s): ",
      paste(significant_clusters, collapse = ", ")
    )

    result <- tryCatch(
      {
        save_edger_group_gene_visualizations(
          context = visualization_context,
          combined_results = edger_objects$combined_results,
          target_gene = target_gene,
          target_ensembl_gene_id = target_ensembl_gene_id,
          significant_clusters = significant_clusters,
          output_root_directory = output_root_directory,
          selection_label = selection_label,
          min_spots_per_sample_cluster = min_spots_per_sample_cluster,
          normalization_scale_factor = normalization_scale_factor,
          upper_colour_quantile = upper_colour_quantile,
          red_palette_colors = red_palette_colors,
          green_palette_colors = green_palette_colors,
          palette_values = palette_values,
          group_colors = group_colors,
          show_histology_image = show_histology_image,
          point_size_no_image = point_size_no_image,
          point_size_with_image = point_size_with_image,
          spatial_pdf_width = spatial_pdf_width,
          spatial_pdf_height = spatial_pdf_height,
          barplot_pdf_width = barplot_pdf_width,
          barplot_pdf_height = barplot_pdf_height,
          skip_completed = skip_completed_genes
        )
      },
      error = function(error_condition) {
        error_condition
      }
    )

    gene_finished_at <- Sys.time()

    elapsed_seconds <- as.numeric(
      difftime(
        gene_finished_at,
        gene_started_at,
        units = "secs"
      )
    )

    if (inherits(result, "error")) {
      current_status <- "failed"
      current_message <- conditionMessage(result)
      message("FAILED: ", target_gene, "\n", current_message)
    } else {
      current_status <- "completed"
      current_message <- ""
      message("Completed: ", target_gene)
    }

    status_rows[[length(status_rows) + 1L]] <- data.frame(
      order = gene_index,
      total_genes = nrow(selection_objects$gene_selection_table),
      target_gene = target_gene,
      ensembl_gene_id = target_ensembl_gene_id,
      significant_clusters = paste(significant_clusters, collapse = ","),
      status = current_status,
      message = current_message,
      elapsed_seconds = round(elapsed_seconds, 2),
      started_at = format(gene_started_at, "%Y-%m-%d %H:%M:%S"),
      finished_at = format(gene_finished_at, "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    current_status_table <- do.call(rbind, status_rows)

    write_edger_group_visualization_tsv(
      current_status_table,
      run_status_file
    )

    if (
      current_status == "failed" &&
        !isTRUE(continue_after_gene_error)
    ) {
      stop(current_message, call. = FALSE)
    }

    invisible(gc(verbose = FALSE))
  }

  run_finished_at <- Sys.time()
  final_status_table <- do.call(rbind, status_rows)

  write_edger_group_visualization_tsv(
    final_status_table,
    run_status_file
  )

  message("\n============================================================")
  message("Finished main donor-group-effect visualization run.")
  message(
    "Completed: ",
    sum(final_status_table$status == "completed")
  )
  message(
    "Failed: ",
    sum(final_status_table$status == "failed")
  )
  message(
    "Elapsed time: ",
    round(
      as.numeric(
        difftime(run_finished_at, run_started_at, units = "mins")
      ),
      2
    ),
    " min"
  )
  message("Output directory:\n", output_root_directory)
  message("============================================================")

  invisible(
    list(
      selected_genes = selection_objects$gene_selection_table,
      selected_gene_cluster_statistics =
        selection_objects$selected_significant_results,
      run_status = final_status_table,
      output_root_directory = output_root_directory,
      selected_genes_file = selected_genes_file,
      run_status_file = run_status_file
    )
  )
}

# ==============================================================================
# End
# ==============================================================================
