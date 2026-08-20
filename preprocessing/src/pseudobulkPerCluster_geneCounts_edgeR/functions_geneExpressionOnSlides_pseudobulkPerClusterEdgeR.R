#!/usr/bin/env Rscript

# ==============================================================================
# functions_geneExpressionOnSlides_pseudobulkPerClusterEdgeR.R
#
# ONE unified functions file for the spatial visualization of the cluster-
# specific pseudobulk edgeR main donor-group effect.
#
# Key rules:
#   - no sample x cluster minimum-spot filter is applied in visualization;
#   - MN/MA/FN/FA percentages shown in statistics headers are read directly
#     from the corresponding edgeR pairwise-comparison result tables;
#   - these group percentages are NOT recomputed from the Seurat object;
#   - Seurat raw counts are still used for the spatial expression maps and
#     sample-level cluster labels shown under individual sections.
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
    panel_border_linewidth = 0.8,
    subtitle_text_override = NULL
) {

  title_text <- paste0("Sample: ", sample_cache$sample_ID)

  if (is.null(subtitle_text_override)) {
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
  } else {
    subtitle_text <- as.character(subtitle_text_override)[[1]]
  }

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
      plot.title = ggplot2::element_text(
        size = 11,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 2)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 7.3,
        hjust = 0.5,
        lineheight = 1.05,
        margin = ggplot2::margin(b = 3)
      ),
      plot.margin = ggplot2::margin(t = 3, r = 3, b = 3, l = 3),
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = panel_border_linewidth
      )
    )
}


# ==============================================================================
# Spatial-only statistics/header/output helpers
# ==============================================================================

check_gene_expression_on_slides_dependencies <- function() {

  required_functions <- c(

    "check_required_packages_edger_group_visualization",

    "write_edger_group_visualization_tsv",

    "sanitize_edger_group_file_component",

    "sort_edger_group_cluster_ids",

    "format_edger_group_number",

    "format_edger_group_probability",

    "build_significant_cluster_label",

    "load_single_seurat_for_edger_group_visualization",

    "load_edger_results_for_group_visualization",

    "select_main_group_genes_for_visualization",

    "flatten_selected_gene_table",

    "prepare_edger_group_visualization_context",

    "prepare_edger_group_gene_data",

    "map_edger_group_values_to_palette",

    "create_edger_group_spatial_panel",

    "create_edger_group_colourbar_legend"

  )

  missing_functions <- required_functions[

    !vapply(

      required_functions,

      exists,

      logical(1),

mode = "function",

inherits = TRUE

    )

  ]

  if (length(missing_functions) > 0L) {

    stop(

      paste0(

        "Missing helper function(s): ",

        paste(missing_functions, collapse = ", "),

        "\nSource functions_geneExpressionOnSlidesAndBarplots_",

        "pseudobulkPerClusterEdgeR.R before this module."

      ),

call. = FALSE

    )

  }

  check_required_packages_edger_group_visualization()

  if (!requireNamespace("ggtext", quietly = TRUE)) {

    stop(

      "Missing required package: ggtext",

call. = FALSE

    )

  }

  invisible(TRUE)

}



# ==============================================================================

# 2. Statistics-table validation and formatting

# ==============================================================================

validate_sex_posthoc_statistics_table <- function(

sex_posthoc_statistics_table

) {

  required_columns <- c(

    "cluster_id",

    "ensembl_gene_id",

    "gene",

    "main_log2FC",

    "main_PValue",

    "main_FDR",

    "male_log2FC",

    "male_PValue",

    "male_FDR_global",

    "male_FDR_selected",

    "female_log2FC",

    "female_PValue",

    "female_FDR_global",

    "female_FDR_selected",

    "mean_percent_positive_spots_Neurotypical_Male",

    "mean_percent_positive_spots_ASD_Male",

    "mean_percent_positive_spots_Neurotypical_Female",

    "mean_percent_positive_spots_ASD_Female",

    "n_main_effect_selected_genes_in_cluster"

  )

  missing_columns <- setdiff(

    required_columns,

    colnames(sex_posthoc_statistics_table)

  )

  if (length(missing_columns) > 0L) {

    stop(

      "Missing sex-posthoc statistics column(s): ",

      paste(missing_columns, collapse = ", "),

call. = FALSE

    )

  }

  invisible(TRUE)

}



format_gene_expression_percent_positive <- function(x) {

  if (

    length(x) == 0L ||

      is.na(x[[1]]) ||

      !is.finite(x[[1]])

  ) {

    return("NA")

  }

  paste0(

    formatC(

      as.numeric(x[[1]]),

format = "f",

digits = 1

    ),

    "%"

  )

}



build_gene_expression_sample_cluster_subtitle <- function(
    gene_data,
    sample_id,
    clusters_per_line = 2L
) {

  sample_id <- as.character(sample_id)[[1]]
  clusters_per_line <- max(1L, as.integer(clusters_per_line)[[1]])

  plot_data <- gene_data$sample_plot_data[[sample_id]]

  if (is.null(plot_data) || nrow(plot_data) == 0L) {
    stop("No plot data found for sample: ", sample_id, call. = FALSE)
  }

  sample_cluster_table <- gene_data$sample_cluster_summary |>
    dplyr::filter(
      .data$sample_ID == .env$sample_id,
      as.character(.data$cluster_id) %in% gene_data$significant_clusters
    ) |>
    dplyr::mutate(
      cluster_id = as.character(.data$cluster_id)
    )

  cluster_labels <- vapply(
    gene_data$significant_clusters,
    function(target_cluster_id) {

      current_row <- sample_cluster_table |>
        dplyr::filter(.data$cluster_id == .env$target_cluster_id)

      if (nrow(current_row) == 0L) {
        return(paste0("C", target_cluster_id, ": absent"))
      }

      current_row <- current_row[1, , drop = FALSE]

      paste0(
        "C",
        target_cluster_id,
        ": ",
        current_row$positive_spots[[1]],
        "/",
        current_row$number_of_spots[[1]],
        " (",
        formatC(
          current_row$percent_positive_spots[[1]],
          format = "f",
          digits = 1
        ),
        "%)"
      )
    },
    character(1)
  )

  label_groups <- split(
    cluster_labels,
    ceiling(seq_along(cluster_labels) / clusters_per_line)
  )

  cluster_text <- paste(
    vapply(
      label_groups,
      paste,
      character(1),
      collapse = " | "
    ),
    collapse = "\n"
  )

  paste0(
    plot_data$sex_std[[1]],
    " | ",
    plot_data$group_std[[1]],
    "\n",
    cluster_text
  )
}


format_gene_expression_statistics_line <- function(

statistics_row,

male_text_colour = "#4C78A8",

female_text_colour = "#C44E52"

) {

  if (nrow(statistics_row) != 1L) {

    stop(

      "Expected exactly one statistics row per gene x cluster.",

call. = FALSE

    )

  }

  main_text <- paste0(

    "C",

    statistics_row$cluster_id[[1]],

    " | Main: log2FC=",

    format_edger_group_number(

      statistics_row$main_log2FC,

digits = 3L

    ),

    ", P=",

    format_edger_group_probability(

      statistics_row$main_PValue

    ),

    ", FDR=",

    format_edger_group_probability(

      statistics_row$main_FDR

    )

  )

  male_text <- paste0(

    "Male: log2FC=",

    format_edger_group_number(

      statistics_row$male_log2FC,

digits = 3L

    ),

    ", P=",

    format_edger_group_probability(

      statistics_row$male_PValue

    ),

    ", FDR(global)=",

    format_edger_group_probability(

      statistics_row$male_FDR_global

    ),

    ", FDR(selected)=",

    format_edger_group_probability(

      statistics_row$male_FDR_selected

    ),

    " | MN:",

    format_gene_expression_percent_positive(

      statistics_row$mean_percent_positive_spots_Neurotypical_Male

    ),

    " | MA:",

    format_gene_expression_percent_positive(

      statistics_row$mean_percent_positive_spots_ASD_Male

    )

  )

  female_text <- paste0(

    "Female: log2FC=",

    format_edger_group_number(

      statistics_row$female_log2FC,

digits = 3L

    ),

    ", P=",

    format_edger_group_probability(

      statistics_row$female_PValue

    ),

    ", FDR(global)=",

    format_edger_group_probability(

      statistics_row$female_FDR_global

    ),

    ", FDR(selected)=",

    format_edger_group_probability(

      statistics_row$female_FDR_selected

    ),

    " | FN:",

    format_gene_expression_percent_positive(

      statistics_row$mean_percent_positive_spots_Neurotypical_Female

    ),

    " | FA:",

    format_gene_expression_percent_positive(

      statistics_row$mean_percent_positive_spots_ASD_Female

    )

  )

  paste0(

    main_text,

    " | ",

    "<span style='color:",

    male_text_colour,

    ";'><b>",

    male_text,

    "</b></span>",

    " | ",

    "<span style='color:",

    female_text_colour,

    ";'><b>",

    female_text,

    "</b></span>"

  )

}



build_gene_expression_statistics_subtitle <- function(

sex_posthoc_statistics_table,

target_ensembl_gene_id,

significant_clusters,

selection_label,

colour_max,

upper_colour_quantile,


male_text_colour = "#4C78A8",

female_text_colour = "#C44E52"

) {

  significant_clusters <- sort_edger_group_cluster_ids(

    significant_clusters

  )

  gene_statistics <- sex_posthoc_statistics_table |>

    dplyr::filter(

      .data$ensembl_gene_id == target_ensembl_gene_id,

      as.character(.data$cluster_id) %in% significant_clusters

    ) |>

    dplyr::mutate(

cluster_numeric = suppressWarnings(

        as.numeric(.data$cluster_id)

      )

    ) |>

    dplyr::arrange(

      .data$cluster_numeric,

      .data$cluster_id

    ) |>

    dplyr::select(-.data$cluster_numeric)

  missing_clusters <- setdiff(

    significant_clusters,

    gene_statistics$cluster_id

  )

  if (length(missing_clusters) > 0L) {

    stop(

      "Missing post-hoc statistics for cluster(s): ",

      paste(missing_clusters, collapse = ", "),

call. = FALSE

    )

  }

  statistics_lines <- vapply(

    seq_len(nrow(gene_statistics)),

    function(row_index) {

      statistics_row <- gene_statistics[

        row_index,

        ,

drop = FALSE

      ]

      format_gene_expression_statistics_line(

statistics_row = statistics_row,

male_text_colour = male_text_colour,

female_text_colour = female_text_colour

      )

    },

    character(1)

  )

  explanatory_line <- paste0(

    "FDR(global) = original edgeR FDR; ",

    "FDR(selected) = BH within genes passing the main-effect filter ",

    "in the same cluster. MN/MA/FN/FA percentages are read directly from ",

    "the corresponding sex-specific edgeR result tables; they are not recalculated for the header."

  )

  layout_line <- paste0(

    "Selection: ",

    selection_label,

    " | Columns: Male neurotypical | Male ASD | Female neurotypical | Female ASD",

    " | common colour range: 0–",

    format_edger_group_number(

      colour_max,

digits = 3L

    ),

    " | upper limit = ",

    upper_colour_quantile * 100,

    "th percentile of positive values"

  )

  paste(

    c(

      statistics_lines,

      explanatory_line,

      layout_line

    ),

collapse = "<br>"

  )

}



# ==============================================================================

# 3. Output paths

# ==============================================================================

get_gene_expression_on_slides_output_paths <- function(

output_root_directory,

target_gene,

significant_clusters

) {

  individual_gene_pdf_directory <- file.path(

    output_root_directory,

    "01_individualGenePdfs"

  )

  dir.create(

    individual_gene_pdf_directory,

recursive = TRUE,

showWarnings = FALSE

  )

  gene_label <- sanitize_edger_group_file_component(

    target_gene

  )

  significant_cluster_label <- build_significant_cluster_label(

    significant_clusters

  )

  gene_pdf_file <- file.path(

    individual_gene_pdf_directory,

    paste0(

      gene_label,

      "_",

      significant_cluster_label,

      "_perSpotExpression_significantClustersRed_otherClustersGreen.pdf"

    )

  )

  list(

individual_gene_pdf_directory = individual_gene_pdf_directory,

gene_pdf_file = gene_pdf_file

  )

}



# ==============================================================================

# 4. Combine per-gene PDFs

# ==============================================================================

combine_gene_expression_on_slides_pdfs <- function(

input_pdf_files,

output_pdf_file

) {

  input_pdf_files <- input_pdf_files[

    file.exists(input_pdf_files)

  ]

  if (length(input_pdf_files) == 0L) {

    return(

      list(

success = FALSE,

method = NA_character_,

message = "No successfully created per-gene PDFs were available."

      )

    )

  }

  if (file.exists(output_pdf_file)) {

    unlink(output_pdf_file)

  }

  if (length(input_pdf_files) == 1L) {

    copied <- file.copy(

from = input_pdf_files[[1]],

to = output_pdf_file,

overwrite = TRUE

    )

    return(

      list(

success = isTRUE(copied) && file.exists(output_pdf_file),

method = "file.copy",

message = "Only one gene PDF existed; copied it as the combined PDF."

      )

    )

  }

  if (requireNamespace("qpdf", quietly = TRUE)) {

    qpdf_ok <- tryCatch(

      {

        qpdf::pdf_combine(

input = input_pdf_files,

output = output_pdf_file

        )

        TRUE

      },

error = function(error_condition) {

        FALSE

      }

    )

    if (

      isTRUE(qpdf_ok) &&

        file.exists(output_pdf_file)

    ) {

      return(

        list(

success = TRUE,

method = "qpdf::pdf_combine",

message = "Combined using qpdf::pdf_combine."

        )

      )

    }

  }

  pdfunite_path <- Sys.which("pdfunite")

  if (nzchar(pdfunite_path)) {

    pdfunite_status <- tryCatch(

      {

        system2(

command = pdfunite_path,

args = c(

            input_pdf_files,

            output_pdf_file

          )

        )

      },

error = function(error_condition) {

        1L

      }

    )

    if (

      identical(

        as.integer(pdfunite_status),

        0L

      ) &&

        file.exists(output_pdf_file)

    ) {

      return(

        list(

success = TRUE,

method = "pdfunite",

message = "Combined using pdfunite."

        )

      )

    }

  }

  list(

success = FALSE,

method = NA_character_,

message = paste0(

      "Per-gene PDFs were created, but automatic merging failed. ",

      "Install the R package qpdf or make pdfunite available in PATH."

    )

  )

}



# ==============================================================================

# 5. Spatial plot with cluster-wise statistics in the subtitle

# ==============================================================================

plot_gene_expression_on_slides_with_statistics <- function(

context,

gene_data,

sex_posthoc_statistics_table,

red_palette_colors,

green_palette_colors,

palette_values,

selection_label,


show_histology_image = FALSE,

point_size_no_image = 1.1,

point_size_with_image = 0.80,

male_text_colour = "#4C78A8",

female_text_colour = "#C44E52",

plot_title_size = 20,

plot_subtitle_size = 8.2,

legend_height_ratio = 0.075

) {

  colour_max <- gene_data$spot_colour_max

  colour_limits <- c(0, colour_max)

  panels <- lapply(

    context$sample_layout$sample_order,

    function(sample_id) {

      plot_data <- gene_data$sample_plot_data[[sample_id]]

      expression_values <- plot_data$logNormalized_expression

      red_colours <- map_edger_group_values_to_palette(

values = expression_values,

palette_colors = red_palette_colors,

palette_values = palette_values,

limits = colour_limits

      )

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

      panel_subtitle <- build_gene_expression_sample_cluster_subtitle(

gene_data = gene_data,

sample_id = sample_id
      )

      create_edger_group_spatial_panel(

sample_cache = context$sample_cache[[sample_id]],

plot_data = plot_data,

colour_column = "plot_colour",

target_gene = gene_data$target_gene,

show_histology_image = show_histology_image,

point_size_no_image = point_size_no_image,

point_size_with_image = point_size_with_image,

subtitle_text_override = panel_subtitle

      )

    }

  )

  names(panels) <- context$sample_layout$sample_order

  significant_cluster_text <- paste(

    paste0(

      "C",

      gene_data$significant_clusters

    ),

collapse = ", "

  )

  red_legend <- create_edger_group_colourbar_legend(

palette_colors = red_palette_colors,

palette_values = palette_values,

colour_max = colour_max,

legend_title = paste0(

      "Statistically significant cluster(s) ",

      significant_cluster_text,

      " — per-spot expression"

    ),

legend_bar_width_mm = 34

  )

  green_legend <- create_edger_group_colourbar_legend(

palette_colors = green_palette_colors,

palette_values = palette_values,

colour_max = colour_max,

legend_title = "Other clusters — per-spot expression",

legend_bar_width_mm = 34

  )

  legend_panel <- red_legend +

    green_legend +

    patchwork::plot_layout(ncol = 2)

  plots_by_column <- lapply(

    context$sample_layout$column_order,

    function(column_key) {

      sample_ids <- context$sample_layout$samples_by_column[[column_key]]

      current_plots <- lapply(

        sample_ids,

        function(sample_id) {

          panels[[sample_id]]

        }

      )

      if (length(current_plots) < context$sample_layout$max_rows) {

        current_plots <- c(

          current_plots,

          rep(

            list(patchwork::plot_spacer()),

            context$sample_layout$max_rows - length(current_plots)

          )

        )

      }

      current_plots

    }

  )

  names(plots_by_column) <- context$sample_layout$column_order

  interleaved_plots <- list()

  for (row_index in seq_len(context$sample_layout$max_rows)) {

    for (column_key in context$sample_layout$column_order) {

      interleaved_plots[[length(interleaved_plots) + 1L]] <-

        plots_by_column[[column_key]][[row_index]]

    }

  }

  spatial_grid <- patchwork::wrap_plots(

    interleaved_plots,

ncol = 4L

  )

  plot_title <- paste0(

    gene_data$target_gene,

    " — per-spot expression: significant clusters red, other clusters green"

  )

  plot_subtitle <- build_gene_expression_statistics_subtitle(

sex_posthoc_statistics_table = sex_posthoc_statistics_table,

target_ensembl_gene_id = gene_data$target_ensembl_gene_id,

significant_clusters = gene_data$significant_clusters,

selection_label = selection_label,

colour_max = colour_max,

upper_colour_quantile = gene_data$upper_colour_quantile,


male_text_colour = male_text_colour,

female_text_colour = female_text_colour

  )

  (

    legend_panel /

      spatial_grid

  ) +

    patchwork::plot_layout(

heights = c(

        legend_height_ratio,

        1

      )

    ) +

    patchwork::plot_annotation(

title = plot_title,

subtitle = plot_subtitle,

theme = ggplot2::theme(

text = ggplot2::element_text(

family = "DejaVu Sans"

        ),

plot.title = ggplot2::element_text(

size = plot_title_size,

face = "bold",

hjust = 0.5

        ),

plot.subtitle = ggtext::element_markdown(

size = plot_subtitle_size,

hjust = 0.5,

lineheight = 1.08,

margin = ggplot2::margin(b = 5)

        ),

plot.margin = ggplot2::margin(

t = 6,

r = 8,

b = 8,

l = 8

        )

      )

    )

}



# ==============================================================================

# 6. Save one gene

# ==============================================================================

save_gene_expression_on_slides_pdf <- function(

context,

sex_posthoc_statistics_table,

target_gene,

target_ensembl_gene_id,

significant_clusters,

output_root_directory,

selection_label,


normalization_scale_factor = 10000,

upper_colour_quantile = 0.99,

red_palette_colors,

green_palette_colors,

palette_values,

show_histology_image = FALSE,

point_size_no_image = 1.1,

point_size_with_image = 0.80,

male_text_colour = "#4C78A8",

female_text_colour = "#C44E52",

plot_title_size = 20,

plot_subtitle_size = 8.2,

legend_height_ratio = 0.075,

pdf_width = 18,

pdf_height = 24.5

) {

  significant_clusters <- sort_edger_group_cluster_ids(

    significant_clusters

  )

  output_paths <- get_gene_expression_on_slides_output_paths(

output_root_directory = output_root_directory,

target_gene = target_gene,

significant_clusters = significant_clusters

  )

  gene_data <- prepare_edger_group_gene_data(

context = context,

target_gene = target_gene,

target_ensembl_gene_id = target_ensembl_gene_id,

significant_clusters = significant_clusters,

normalization_scale_factor = normalization_scale_factor,

upper_colour_quantile = upper_colour_quantile

  )

  spatial_plot <- plot_gene_expression_on_slides_with_statistics(

context = context,

gene_data = gene_data,

sex_posthoc_statistics_table = sex_posthoc_statistics_table,

red_palette_colors = red_palette_colors,

green_palette_colors = green_palette_colors,

palette_values = palette_values,

selection_label = selection_label,


show_histology_image = show_histology_image,

point_size_no_image = point_size_no_image,

point_size_with_image = point_size_with_image,

male_text_colour = male_text_colour,

female_text_colour = female_text_colour,

plot_title_size = plot_title_size,

plot_subtitle_size = plot_subtitle_size,

legend_height_ratio = legend_height_ratio

  )

  pdf_device <- if (capabilities("cairo")) {

    grDevices::cairo_pdf

  } else {

    grDevices::pdf

  }

  ggplot2::ggsave(

filename = output_paths$gene_pdf_file,

plot = spatial_plot,

device = pdf_device,

width = pdf_width,

height = pdf_height,

units = "in",

limitsize = FALSE,

bg = "white"

  )

  if (!file.exists(output_paths$gene_pdf_file)) {

    stop(

      "Per-gene spatial PDF was not created: ",

      output_paths$gene_pdf_file,

call. = FALSE

    )

  }

  invisible(

    list(

gene_pdf_file = output_paths$gene_pdf_file,

gene_data = gene_data

    )

  )

}



# ==============================================================================

# 7. Complete workflow

# ==============================================================================

run_main_effect_group_gene_expression_on_slides <- function(

input_seurat_rdata_file,

edger_results_rdata_file,

sex_posthoc_statistics_table,

output_root_directory,

cluster_column = "leiden_res0.40",

assay_name = "RNA",

sample_id_column = "sample_ID",

group_column = "fmt_donor_group",

sex_column = "sex",

included_sample_ids = NULL,

seurat_object_name = NULL,

fdr_threshold = 0.05,

abs_log2fc_threshold = 0.5,


maximum_number_of_genes = 2L,

manual_target_genes = NULL,

normalization_scale_factor = 10000,

upper_colour_quantile = 0.99,

red_palette_colors,

green_palette_colors,

palette_values,

show_histology_image = FALSE,

point_size_no_image = 1.1,

point_size_with_image = 0.80,

male_text_colour = "#4C78A8",

female_text_colour = "#C44E52",

plot_title_size = 20,

plot_subtitle_size = 8.2,

legend_height_ratio = 0.075,

pdf_width = 18,

pdf_height = 24.5,

continue_after_gene_error = TRUE

) {

  check_gene_expression_on_slides_dependencies()

  validate_sex_posthoc_statistics_table(

    sex_posthoc_statistics_table

  )

  dir.create(

    output_root_directory,

recursive = TRUE,

showWarnings = FALSE

  )

  individual_gene_pdf_directory <- file.path(

    output_root_directory,

    "01_individualGenePdfs"

  )

  dir.create(

    individual_gene_pdf_directory,

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

    "00_runStatus_mainEffectGroup_expressionOnSlides.tsv"

  )

  combined_pdf_file <- file.path(

    output_root_directory,

    "02_mainEffectGroup_expressionOnSlides_FDR005_absLog2FC05.pdf"

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

  message(

    "Selected genes: ",

    nrow(selection_objects$gene_selection_table)

  )

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

  created_gene_pdfs <- character(0)

  run_started_at <- Sys.time()

  for (gene_index in seq_len(nrow(selection_objects$gene_selection_table))) {

    target_gene <- selection_objects$gene_selection_table$gene[[gene_index]]

    target_ensembl_gene_id <- selection_objects$gene_selection_table$ensembl_gene_id[[gene_index]]

    significant_clusters <- selection_objects$gene_selection_table$significant_clusters[[gene_index]]

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

        save_gene_expression_on_slides_pdf(

context = visualization_context,

sex_posthoc_statistics_table = sex_posthoc_statistics_table,

target_gene = target_gene,

target_ensembl_gene_id = target_ensembl_gene_id,

significant_clusters = significant_clusters,

output_root_directory = output_root_directory,

selection_label = selection_label,


normalization_scale_factor = normalization_scale_factor,

upper_colour_quantile = upper_colour_quantile,

red_palette_colors = red_palette_colors,

green_palette_colors = green_palette_colors,

palette_values = palette_values,

show_histology_image = show_histology_image,

point_size_no_image = point_size_no_image,

point_size_with_image = point_size_with_image,

male_text_colour = male_text_colour,

female_text_colour = female_text_colour,

plot_title_size = plot_title_size,

plot_subtitle_size = plot_subtitle_size,

legend_height_ratio = legend_height_ratio,

pdf_width = pdf_width,

pdf_height = pdf_height

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

      current_gene_pdf_file <- ""

      message("FAILED: ", target_gene, "\n", current_message)

    } else {

      current_status <- "completed"

      current_message <- ""

      current_gene_pdf_file <- result$gene_pdf_file

      created_gene_pdfs <- c(

        created_gene_pdfs,

        current_gene_pdf_file

      )

      message("Completed: ", target_gene)

      message("Saved gene PDF: ", current_gene_pdf_file)

    }

    status_rows[[length(status_rows) + 1L]] <- data.frame(

order = gene_index,

total_genes = nrow(selection_objects$gene_selection_table),

target_gene = target_gene,

ensembl_gene_id = target_ensembl_gene_id,

significant_clusters = paste(significant_clusters, collapse = ","),

status = current_status,

message = current_message,

gene_pdf_file = current_gene_pdf_file,

elapsed_seconds = round(elapsed_seconds, 2),

started_at = format(

        gene_started_at,

        "%Y-%m-%d %H:%M:%S"

      ),

finished_at = format(

        gene_finished_at,

        "%Y-%m-%d %H:%M:%S"

      ),

stringsAsFactors = FALSE,

check.names = FALSE

    )

    write_edger_group_visualization_tsv(

      do.call(rbind, status_rows),

      run_status_file

    )

    if (

      current_status == "failed" &&

        !isTRUE(continue_after_gene_error)

    ) {

      stop(

        current_message,

call. = FALSE

      )

    }

    invisible(gc(verbose = FALSE))

  }

  combine_status <- combine_gene_expression_on_slides_pdfs(

input_pdf_files = created_gene_pdfs,

output_pdf_file = combined_pdf_file

  )

  run_finished_at <- Sys.time()

  final_status_table <- do.call(

    rbind,

    status_rows

  )

  write_edger_group_visualization_tsv(

    final_status_table,

    run_status_file

  )

  message("\n============================================================")

  message("Finished spatial-only gene-expression visualization run.")

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

        difftime(

          run_finished_at,

          run_started_at,

units = "mins"

        )

      ),

      2

    ),

    " min"

  )

  message("Per-gene PDF directory:")

  message(individual_gene_pdf_directory)

  message("Combined PDF:")

  message(combined_pdf_file)

  message("Combine status:")

  message(combine_status$message)

  message("Output directory:")

  message(output_root_directory)

  message("============================================================")

  invisible(

    list(

selected_genes = selection_objects$gene_selection_table,

selected_gene_cluster_statistics =

        selection_objects$selected_significant_results,

run_status = final_status_table,

output_root_directory = output_root_directory,

individual_gene_pdf_directory = individual_gene_pdf_directory,

combined_pdf_file = combined_pdf_file,

combined_pdf_created = isTRUE(combine_status$success),

combine_status = combine_status,

selected_genes_file = selected_genes_file,

selected_gene_cluster_statistics_file =

        selected_gene_cluster_statistics_file,

run_status_file = run_status_file

    )

  )

}

# ==============================================================================

# End

# ==============================================================================