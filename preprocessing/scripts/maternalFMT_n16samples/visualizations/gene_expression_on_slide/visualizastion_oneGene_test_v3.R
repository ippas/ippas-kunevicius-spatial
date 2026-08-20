#!/usr/bin/env Rscript

# ==============================================================================
# 05_geneExpression_spatialMaps_groupBarplot_v4.R
#
# Purpose:
# - visualize one selected gene in 16 accepted maternal FMT Visium samples
# - use one common colour scale across all 16 samples for per-spot expression
# - arrange panels into 4 fixed columns:
#     1. Male Neurotypical
#     2. Male ASD
#     3. Female Neurotypical
#     4. Female ASD
# - save all outputs directly inside a gene-specific folder, without
#   additional pdf/ or summary/ subfolders
# - create three main outputs:
#     1. standard spatial-map PDF (true per-spot expression)
#     2. group barplot PDF (sample means by group)
#     3. mean-filled spatial-map PDF (all spots in a sample filled with the
#        sample-level mean expression)
# - save sample-level and group-level summaries in TSV format
#
# Notes:
# - in the no-image mode, point positions are based on scaled Space Ranger
#   pixel coordinates, but the histology image itself is not drawn.
# - Per-spot normalization is:
#     log1p(gene_count / total_UMI_in_spot * 10000)
# - total UMI is always calculated from ALL genes
# - group barplot uses sample means as observations (not spots)
# ==============================================================================


# ==============================================================================
# 1. Packages
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(png)
  library(jsonlite)
  library(grid)
})


# ==============================================================================
# 2. User settings
# ==============================================================================

project_dir <- "/home/mateusz/projects/ippas-kunevicius-spatial"

path_to_data <- file.path(
  project_dir,
  "data",
  "spacerangerCount_withoutJSON_GRCm39-2024-A_2026-03-28"
)

metadata_file <- file.path(
  project_dir,
  "data",
  "metadata_autismFMT.tsv"
)

output_root_dir <- file.path(
  project_dir,
  "results",
  "maternalFMT_n16samples",
  "geneExpression_visualtionSlide",
  "test_visulation_v3"
)

# --------------------------------------------------------------------------
# Main switches
# --------------------------------------------------------------------------

target_gene <- "Arc"

# FALSE = lightweight mode, no histology image, common slide-like frame
# TRUE  = include histology image in the background
show_histology_image <- FALSE

# FALSE = preserve original Space Ranger orientation
# TRUE  = rotate the no-image view by 90 degrees clockwise
rotate_no_image_90 <- FALSE

normalization_scale_factor <- 10000
upper_colour_quantile <- 0.99

# --------------------------------------------------------------------------
# Plot appearance
# --------------------------------------------------------------------------

number_of_columns <- 4
combined_pdf_width_inches <- 18
combined_pdf_height_inches <- 23
mean_fill_pdf_width_inches <- 18
mean_fill_pdf_height_inches <- 23

barplot_pdf_width_inches <- 8.7
barplot_pdf_height_inches <- 7.0
barplot_bar_width <- 0.68
barplot_point_size <- 3.2
barplot_jitter_width <- 0.08
barplot_errorbar_width <- 0.18

point_size_no_image <- 1.1
point_size_with_image <- 0.80

panel_padding_fraction_no_image <- 0.03
panel_padding_fraction_with_image <- 0.03

panel_border_linewidth <- 0.8

# Strong red palette for spatial maps
palette_colors <- c(
  "#D9D9D9",
  "#FEE5D9",
  "#FCAE91",
  "#FB6A4A",
  "#DE2D26",
  "#A50F15",
  "#67000D"
)

palette_values <- c(
  0.000,
  0.006,
  0.060,
  0.180,
  0.400,
  0.700,
  1.000
)

# --------------------------------------------------------------------------
# Samples
# --------------------------------------------------------------------------

selected_samples <- c(
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

excluded_samples <- c(
  "12_3F",
  "15_1M",
  "20_1F",
  "20_3M"
)

expected_counts_by_column <- c(
  male_neurotypical = 3,
  male_asd = 5,
  female_neurotypical = 3,
  female_asd = 5
)


# ==============================================================================
# 3. Check paths and create output
# ==============================================================================

if (!dir.exists(project_dir)) {
  stop("Project directory does not exist: ", project_dir)
}

if (!dir.exists(path_to_data)) {
  stop("Space Ranger data directory does not exist: ", path_to_data)
}

if (!file.exists(metadata_file)) {
  stop("Metadata file does not exist: ", metadata_file)
}

gene_output_dir <- file.path(output_root_dir, target_gene)

for (directory in c(output_root_dir, gene_output_dir)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

mode_label <- if (show_histology_image) {
  "withHistology"
} else if (rotate_no_image_90) {
  "noHistology_pixelFrame_rotated90"
} else {
  "noHistology_pixelFrame"
}

output_spatial_pdf <- file.path(
  gene_output_dir,
  paste0(
    "05_",
    target_gene,
    "_logNormalized_spatialMaps_n16_",
    mode_label,
    "_fourColumnsSexGroup.pdf"
  )
)

output_barplot_pdf <- file.path(
  gene_output_dir,
  paste0(
    "05_",
    target_gene,
    "_sampleMean_groupBarplot_meanSE_withPercentPositive.pdf"
  )
)

output_mean_fill_pdf <- file.path(
  gene_output_dir,
  paste0(
    "05_",
    target_gene,
    "_sampleMeanFilled_spatialMaps_n16_",
    mode_label,
    "_fourColumnsSexGroup.pdf"
  )
)

output_sample_summary_tsv <- file.path(
  gene_output_dir,
  paste0(
    "05_",
    target_gene,
    "_logNormalized_spatialMaps_n16_",
    mode_label,
    "_summary.tsv"
  )
)

output_group_summary_tsv <- file.path(
  gene_output_dir,
  paste0(
    "05_",
    target_gene,
    "_sampleMean_group_summary.tsv"
  )
)

message("Target gene: ", target_gene)
message("Show histology image: ", show_histology_image)
message("Rotate no-image mode by 90 degrees: ", rotate_no_image_90)
message("Output root: ", output_root_dir)
message("Gene output directory: ", gene_output_dir)
message("Selected samples: ", paste(selected_samples, collapse = ", "))
message("Excluded samples: ", paste(excluded_samples, collapse = ", "))


# ==============================================================================
# 4. Helper functions
# ==============================================================================

choose_matrix_file <- function(outs_dir) {

  candidate_files <- c(
    file.path(outs_dir, "filtered_feature_bc_matrix.h5"),
    file.path(outs_dir, "raw_feature_bc_matrix.h5")
  )

  existing_files <- candidate_files[file.exists(candidate_files)]

  if (length(existing_files) == 0) {
    stop(
      "Neither filtered_feature_bc_matrix.h5 nor raw_feature_bc_matrix.h5 exists in: ",
      outs_dir
    )
  }

  existing_files[[1]]
}


choose_tissue_positions_file <- function(spatial_dir) {

  candidate_files <- c(
    file.path(spatial_dir, "tissue_positions.csv"),
    file.path(spatial_dir, "tissue_positions_list.csv")
  )

  existing_files <- candidate_files[file.exists(candidate_files)]

  if (length(existing_files) == 0) {
    stop("No tissue positions file found in: ", spatial_dir)
  }

  existing_files[[1]]
}


read_tissue_positions <- function(filename) {

  if (basename(filename) == "tissue_positions.csv") {
    coordinates <- read.csv(
      filename,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    coordinates <- read.csv(
      filename,
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      col.names = c(
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres"
      )
    )
  }

  required_columns <- c(
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"
  )

  missing_columns <- setdiff(required_columns, colnames(coordinates))

  if (length(missing_columns) > 0) {
    stop(
      "Missing columns in ",
      filename,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  coordinates |>
    transmute(
      barcode = trimws(as.character(barcode)),
      in_tissue = as.integer(in_tissue),
      array_row = as.numeric(array_row),
      array_col = as.numeric(array_col),
      pxl_row_in_fullres = as.numeric(pxl_row_in_fullres),
      pxl_col_in_fullres = as.numeric(pxl_col_in_fullres)
    )
}


extract_gene_expression_matrix <- function(read10x_result) {

  if (inherits(read10x_result, "Matrix") || is.matrix(read10x_result)) {
    return(read10x_result)
  }

  if (!is.list(read10x_result)) {
    stop("Read10X_h5 returned an unsupported object type.")
  }

  if ("Gene Expression" %in% names(read10x_result)) {
    return(read10x_result[["Gene Expression"]])
  }

  if (length(read10x_result) == 1) {
    return(read10x_result[[1]])
  }

  stop(
    "Read10X_h5 returned multiple matrices, but 'Gene Expression' could not be identified."
  )
}


find_target_gene <- function(feature_names, requested_gene) {

  exact_index <- which(feature_names == requested_gene)
  if (length(exact_index) == 1) {
    return(exact_index)
  }

  ci_index <- which(tolower(feature_names) == tolower(requested_gene))
  if (length(ci_index) == 1) {
    warning(
      "Gene ",
      requested_gene,
      " matched case-insensitively as ",
      feature_names[ci_index]
    )
    return(ci_index)
  }

  similar_features <- grep(
    requested_gene,
    feature_names,
    ignore.case = TRUE,
    value = TRUE
  )

  stop(
    "Target gene was not found: ",
    requested_gene,
    if (length(similar_features) > 0) {
      paste0(
        "\nSimilar feature names: ",
        paste(head(similar_features, 20), collapse = ", ")
      )
    } else {
      ""
    }
  )
}


format_metadata_value <- function(x) {
  if (length(x) == 0 || is.na(x) || trimws(as.character(x)) == "") {
    return("NA")
  }
  as.character(x)
}


standardize_sex <- function(x) {

  x_lower <- tolower(trimws(as.character(x)))

  if (x_lower %in% c("m", "male")) {
    return("Male")
  }

  if (x_lower %in% c("f", "female")) {
    return("Female")
  }

  stop("Unsupported sex label: ", x)
}


standardize_group <- function(x) {

  x_lower <- tolower(trimws(as.character(x)))

  if (grepl("asd|autism", x_lower)) {
    return("ASD")
  }

  if (grepl("neurotypical|control|typical|nt", x_lower)) {
    return("Neurotypical")
  }

  stop(
    "Unsupported donor-group label: ",
    x,
    "\nExpected something matching ASD/autism or neurotypical/control."
  )
}


compute_square_limits <- function(xmin, xmax, ymin, ymax, padding_fraction = 0.03) {

  x_range <- xmax - xmin
  y_range <- ymax - ymin
  max_range <- max(x_range, y_range)

  if (!is.finite(max_range) || max_range <= 0) {
    max_range <- 1
  }

  half_side <- (max_range / 2) * (1 + padding_fraction)
  x_center <- (xmin + xmax) / 2
  y_center <- (ymin + ymax) / 2

  list(
    x_limits = c(x_center - half_side, x_center + half_side),
    y_limits = c(y_center - half_side, y_center + half_side)
  )
}


choose_image_file_and_scale <- function(spatial_dir) {

  image_candidates <- c(
    file.path(spatial_dir, "tissue_lowres_image.png"),
    file.path(spatial_dir, "tissue_hires_image.png")
  )

  image_file <- image_candidates[file.exists(image_candidates)][1]

  if (is.na(image_file) || length(image_file) == 0) {
    stop("No tissue image PNG found in: ", spatial_dir)
  }

  scalefactors_file <- file.path(spatial_dir, "scalefactors_json.json")
  if (!file.exists(scalefactors_file)) {
    stop("scalefactors_json.json not found in: ", spatial_dir)
  }

  scalefactors <- jsonlite::fromJSON(scalefactors_file)

  if (basename(image_file) == "tissue_lowres_image.png") {
    scale_factor <- scalefactors$tissue_lowres_scalef
  } else {
    scale_factor <- scalefactors$tissue_hires_scalef
  }

  if (!is.numeric(scale_factor) || length(scale_factor) != 1 || is.na(scale_factor)) {
    stop("Could not determine image scale factor for: ", image_file)
  }

  image_array <- png::readPNG(image_file)
  image_height <- dim(image_array)[1]
  image_width <- dim(image_array)[2]
  rm(image_array)

  list(
    image_file = image_file,
    scale_factor = scale_factor,
    image_width = image_width,
    image_height = image_height
  )
}


make_column_key <- function(sex_std, group_std) {

  if (sex_std == "Male" && group_std == "Neurotypical") {
    return("male_neurotypical")
  }

  if (sex_std == "Male" && group_std == "ASD") {
    return("male_asd")
  }

  if (sex_std == "Female" && group_std == "Neurotypical") {
    return("female_neurotypical")
  }

  if (sex_std == "Female" && group_std == "ASD") {
    return("female_asd")
  }

  stop(
    "Unsupported sex/group combination: ",
    sex_std,
    " / ",
    group_std
  )
}


make_plot_title_line <- function(plot_data) {
  paste0(
    "Sample: ",
    plot_data$sample_ID[[1]]
  )
}


make_plot_subtitle_line <- function(plot_data) {

  n_positive <- sum(plot_data$target_raw_count > 0)
  percent_positive <- 100 * mean(plot_data$target_raw_count > 0)

  mean_expression_all_spots <- mean(
    plot_data$logNormalized_expression,
    na.rm = TRUE
  )

  sd_expression_all_spots <- sd(
    plot_data$logNormalized_expression,
    na.rm = TRUE
  )

  paste0(
    "Group: ",
    format_metadata_value(plot_data$group_std[[1]]),
    " | Sex: ",
    format_metadata_value(plot_data$sex_std[[1]]),
    "\nTissue spots: ",
    format(nrow(plot_data), big.mark = " ", scientific = FALSE),
    " | ",
    target_gene,
    "+: ",
    format(n_positive, big.mark = " ", scientific = FALSE),
    " (",
    format(round(percent_positive, 1), nsmall = 1),
    "%)",
    "\nMean ± SD (all spots): ",
    format(round(mean_expression_all_spots, 3), nsmall = 3, trim = TRUE),
    " ± ",
    format(round(sd_expression_all_spots, 3), nsmall = 3, trim = TRUE)
  )
}


make_mean_fill_subtitle_line <- function(plot_data) {

  n_positive <- sum(plot_data$target_raw_count > 0)
  percent_positive <- 100 * mean(plot_data$target_raw_count > 0)

  mean_expression_all_spots <- mean(
    plot_data$logNormalized_expression,
    na.rm = TRUE
  )

  paste0(
    "Group: ",
    format_metadata_value(plot_data$group_std[[1]]),
    " | Sex: ",
    format_metadata_value(plot_data$sex_std[[1]]),
    "\nTissue spots: ",
    format(nrow(plot_data), big.mark = " ", scientific = FALSE),
    " | ",
    target_gene,
    "+: ",
    format(n_positive, big.mark = " ", scientific = FALSE),
    " (",
    format(round(percent_positive, 1), nsmall = 1),
    "%)",
    "\nAll spots filled with sample mean: ",
    format(round(mean_expression_all_spots, 3), nsmall = 3, trim = TRUE)
  )
}


# ==============================================================================
# 5. Read metadata and define sample order
# ==============================================================================

metadata_autismFMT <- read.delim(
  metadata_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_metadata_columns <- c(
  "sample_ID",
  "fmt_donor_group",
  "sex"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata_autismFMT)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

sample_table <- metadata_autismFMT |>
  transmute(
    sample_ID = trimws(as.character(sample_ID)),
    donor_group_raw = trimws(as.character(fmt_donor_group)),
    sex_raw = trimws(as.character(sex))
  ) |>
  filter(
    !is.na(sample_ID),
    sample_ID != "",
    sample_ID %in% selected_samples
  ) |>
  distinct(sample_ID, .keep_all = TRUE)

missing_selected_samples <- setdiff(selected_samples, sample_table$sample_ID)

if (length(missing_selected_samples) > 0) {
  stop(
    "Selected samples are missing in metadata: ",
    paste(missing_selected_samples, collapse = ", ")
  )
}

sample_table <- tibble(
  sample_ID = selected_samples
) |>
  left_join(sample_table, by = "sample_ID") |>
  mutate(
    sex_std = vapply(sex_raw, standardize_sex, character(1)),
    group_std = vapply(donor_group_raw, standardize_group, character(1)),
    column_key = mapply(make_column_key, sex_std, group_std)
  )

if (anyDuplicated(sample_table$sample_ID) > 0) {
  stop("Duplicated sample IDs were detected in sample_table.")
}

message("")
message("Sample assignment:")
print(
  sample_table |>
    select(sample_ID, sex_raw, donor_group_raw, sex_std, group_std, column_key),
  n = Inf,
  width = Inf
)

observed_counts_by_column <- table(sample_table$column_key)

message("")
message("Observed counts by requested column layout:")
print(observed_counts_by_column)

for (key in names(expected_counts_by_column)) {
  observed_n <- if (key %in% names(observed_counts_by_column)) {
    as.integer(observed_counts_by_column[[key]])
  } else {
    0L
  }

  expected_n <- expected_counts_by_column[[key]]

  if (observed_n != expected_n) {
    warning(
      "Column '",
      key,
      "' has ",
      observed_n,
      " samples, but ",
      expected_n,
      " were expected."
    )
  }
}


# ==============================================================================
# 6. Process samples one at a time
# ==============================================================================

sample_plot_data <- vector("list", length(selected_samples))
names(sample_plot_data) <- selected_samples

sample_summary <- vector("list", length(selected_samples))
names(sample_summary) <- selected_samples

image_info_list <- vector("list", length(selected_samples))
names(image_info_list) <- selected_samples

for (sample_id in selected_samples) {

  message("")
  message("============================================================")
  message("Processing sample: ", sample_id)
  message("============================================================")

  sample_outs_dir <- file.path(path_to_data, sample_id, "outs")
  spatial_dir <- file.path(sample_outs_dir, "spatial")

  if (!dir.exists(sample_outs_dir)) {
    stop("Missing outs directory for sample ", sample_id, ": ", sample_outs_dir)
  }

  if (!dir.exists(spatial_dir)) {
    stop("Missing spatial directory for sample ", sample_id, ": ", spatial_dir)
  }

  matrix_file <- choose_matrix_file(sample_outs_dir)
  positions_file <- choose_tissue_positions_file(spatial_dir)

  coordinates <- read_tissue_positions(positions_file) |>
    filter(
      !is.na(barcode),
      barcode != "",
      in_tissue == 1
    )

  if (nrow(coordinates) == 0) {
    stop("No in-tissue spots found for sample: ", sample_id)
  }

  read10x_result <- Seurat::Read10X_h5(
    filename = matrix_file,
    use.names = TRUE,
    unique.features = TRUE
  )

  counts <- extract_gene_expression_matrix(read10x_result)
  rm(read10x_result)

  if (nrow(counts) == 0 || ncol(counts) == 0) {
    stop("Empty count matrix for sample: ", sample_id)
  }

  common_barcodes <- coordinates$barcode[
    coordinates$barcode %in% colnames(counts)
  ]

  if (length(common_barcodes) == 0) {
    stop("No shared barcodes between coordinates and matrix for sample: ", sample_id)
  }

  coordinates <- coordinates |>
    filter(barcode %in% common_barcodes) |>
    arrange(match(barcode, common_barcodes))

  counts <- counts[, coordinates$barcode, drop = FALSE]

  if (!identical(colnames(counts), coordinates$barcode)) {
    stop("Barcode order mismatch for sample: ", sample_id)
  }

  target_gene_index <- find_target_gene(
    feature_names = rownames(counts),
    requested_gene = target_gene
  )

  target_gene_name_in_matrix <- rownames(counts)[target_gene_index]

  total_umi_per_spot <- Matrix::colSums(counts)
  target_raw_counts <- as.numeric(counts[target_gene_index, , drop = TRUE])

  log_normalized_expression <- numeric(length(total_umi_per_spot))
  valid_spots <- total_umi_per_spot > 0

  log_normalized_expression[valid_spots] <- log1p(
    target_raw_counts[valid_spots] /
      total_umi_per_spot[valid_spots] *
      normalization_scale_factor
  )

  sample_metadata <- sample_table |>
    filter(sample_ID == sample_id)

  if (nrow(sample_metadata) != 1) {
    stop("Expected one metadata row for sample ", sample_id)
  }

  plot_data <- coordinates |>
    mutate(
      sample_ID = sample_id,
      donor_group_raw = sample_metadata$donor_group_raw[[1]],
      sex_raw = sample_metadata$sex_raw[[1]],
      group_std = sample_metadata$group_std[[1]],
      sex_std = sample_metadata$sex_std[[1]],
      column_key = sample_metadata$column_key[[1]],
      target_gene = target_gene_name_in_matrix,
      target_raw_count = target_raw_counts,
      total_UMI = as.numeric(total_umi_per_spot),
      logNormalized_expression = log_normalized_expression
    ) |>
    select(
      sample_ID,
      donor_group_raw,
      sex_raw,
      group_std,
      sex_std,
      column_key,
      barcode,
      in_tissue,
      array_row,
      array_col,
      pxl_row_in_fullres,
      pxl_col_in_fullres,
      target_gene,
      target_raw_count,
      total_UMI,
      logNormalized_expression
    )

  sample_plot_data[[sample_id]] <- plot_data
  image_info_list[[sample_id]] <- choose_image_file_and_scale(spatial_dir)

  sample_summary[[sample_id]] <- tibble(
    sample_ID = sample_id,
    donor_group_raw = sample_metadata$donor_group_raw[[1]],
    sex_raw = sample_metadata$sex_raw[[1]],
    group_std = sample_metadata$group_std[[1]],
    sex_std = sample_metadata$sex_std[[1]],
    column_key = sample_metadata$column_key[[1]],
    matrix_file = matrix_file,
    positions_file = positions_file,
    target_gene_requested = target_gene,
    target_gene_in_matrix = target_gene_name_in_matrix,
    n_tissue_spots = nrow(plot_data),
    total_UMI_all_tissue_spots = sum(plot_data$total_UMI),
    median_total_UMI_per_tissue_spot = median(plot_data$total_UMI),
    n_target_positive_spots = sum(plot_data$target_raw_count > 0),
    percent_target_positive_spots = 100 * mean(plot_data$target_raw_count > 0),
    total_target_raw_count = sum(plot_data$target_raw_count),
    mean_target_logNormalized_all_spots =
      mean(plot_data$logNormalized_expression),
    sd_target_logNormalized_all_spots =
      sd(plot_data$logNormalized_expression),
    median_target_logNormalized_all_spots =
      median(plot_data$logNormalized_expression),
    max_target_logNormalized_all_spots =
      max(plot_data$logNormalized_expression)
  )

  rm(
    counts,
    coordinates,
    total_umi_per_spot,
    target_raw_counts,
    log_normalized_expression,
    plot_data
  )

  invisible(gc(verbose = FALSE))
}


# ==============================================================================
# 7. Combine processed data and define common colour scales
# ==============================================================================

all_plot_data <- bind_rows(sample_plot_data)

summary_table <- bind_rows(sample_summary) |>
  mutate(
    sample_ID = factor(sample_ID, levels = selected_samples)
  ) |>
  arrange(sample_ID) |>
  mutate(
    sample_ID = as.character(sample_ID)
  )

write.table(
  summary_table,
  file = output_sample_summary_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

barplot_group_order <- c(
  "male_neurotypical",
  "male_asd",
  "female_neurotypical",
  "female_asd"
)

barplot_group_labels <- c(
  male_neurotypical = "Male neurotypical",
  male_asd = "Male ASD",
  female_neurotypical = "Female neurotypical",
  female_asd = "Female ASD"
)

barplot_axis_labels_base <- c(
  "Male neurotypical" = "Male\nneurotypical",
  "Male ASD" = "Male\nASD",
  "Female neurotypical" = "Female\nneurotypical",
  "Female ASD" = "Female\nASD"
)

barplot_sample_data <- summary_table |>
  mutate(
    group_key = factor(column_key, levels = barplot_group_order),
    group_label = factor(
      barplot_group_labels[as.character(group_key)],
      levels = unname(barplot_group_labels[barplot_group_order])
    )
  )

if (any(is.na(barplot_sample_data$group_key))) {
  stop("At least one sample could not be assigned to a barplot group.")
}

group_summary_table <- barplot_sample_data |>
  group_by(group_key, group_label) |>
  summarise(
    n_samples = n(),
    mean_sampleMean_logNormalized =
      mean(mean_target_logNormalized_all_spots, na.rm = TRUE),
    sd_sampleMean_logNormalized =
      sd(mean_target_logNormalized_all_spots, na.rm = TRUE),
    se_sampleMean_logNormalized =
      sd_sampleMean_logNormalized / sqrt(n_samples),
    mean_percent_target_positive_spots =
      mean(percent_target_positive_spots, na.rm = TRUE),
    sd_percent_target_positive_spots =
      sd(percent_target_positive_spots, na.rm = TRUE),
    se_percent_target_positive_spots =
      sd_percent_target_positive_spots / sqrt(n_samples),
    .groups = "drop"
  ) |>
  arrange(group_key) |>
  mutate(
    group_key = as.character(group_key),
    group_label = as.character(group_label)
  )

write.table(
  group_summary_table,
  file = output_group_summary_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

all_expression_values <- all_plot_data$logNormalized_expression
positive_expression_values <- all_expression_values[
  is.finite(all_expression_values) &
    all_expression_values > 0
]

common_colour_min <- 0

if (length(positive_expression_values) == 0) {
  warning("No positive ", target_gene, " expression values found. Using range 0-1.")
  common_colour_max <- 1
} else {
  common_colour_max <- as.numeric(
    quantile(
      positive_expression_values,
      probs = upper_colour_quantile,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )

  if (!is.finite(common_colour_max) || common_colour_max <= 0) {
    common_colour_max <- max(positive_expression_values, na.rm = TRUE)
  }
}

sample_mean_colour_min <- 0
sample_mean_values <- summary_table$mean_target_logNormalized_all_spots
sample_mean_positive_values <- sample_mean_values[
  is.finite(sample_mean_values) & sample_mean_values > 0
]

if (length(sample_mean_positive_values) == 0) {
  warning("No positive sample means found. Using range 0-1 for mean-fill map.")
  sample_mean_colour_max <- 1
} else {
  sample_mean_colour_max <- max(sample_mean_positive_values, na.rm = TRUE)

  if (!is.finite(sample_mean_colour_max) || sample_mean_colour_max <= 0) {
    sample_mean_colour_max <- 1
  }
}

message("")
message("Per-spot colour scale:")
message("Minimum: ", round(common_colour_min, 4))
message(
  "Maximum: ",
  round(common_colour_max, 4),
  " (",
  upper_colour_quantile * 100,
  "th percentile of positive per-spot values)"
)

message("")
message("Mean-fill colour scale:")
message("Minimum: ", round(sample_mean_colour_min, 4))
message("Maximum: ", round(sample_mean_colour_max, 4), " (max sample mean)")


# ==============================================================================
# 8. Define no-image geometry strategy
# ==============================================================================

if (!show_histology_image) {
  message("")
  message(
    "No-image mode uses scaled Space Ranger pixel coordinates, not array_row/array_col."
  )
  message(
    "This preserves x/y geometry and helps keep spots visually circular."
  )
  message(
    "Optional 90-degree rotation in no-image mode: ", rotate_no_image_90
  )
}


# ==============================================================================
# 9. Plot functions
# ==============================================================================

prepare_no_image_geometry <- function(sample_id) {

  plot_data <- sample_plot_data[[sample_id]]
  image_info <- image_info_list[[sample_id]]

  if (is.null(image_info)) {
    stop("Missing image information for sample: ", sample_id)
  }

  plot_data <- plot_data |>
    mutate(
      x_plot = pxl_col_in_fullres * image_info$scale_factor,
      y_plot = pxl_row_in_fullres * image_info$scale_factor
    )

  if (rotate_no_image_90) {

    plot_data <- plot_data |>
      mutate(
        x_use = y_plot,
        y_use = image_info$image_width - x_plot
      )

    frame_limits <- compute_square_limits(
      xmin = 0,
      xmax = image_info$image_height,
      ymin = 0,
      ymax = image_info$image_width,
      padding_fraction = panel_padding_fraction_no_image
    )

  } else {

    plot_data <- plot_data |>
      mutate(
        x_use = x_plot,
        y_use = y_plot
      )

    frame_limits <- compute_square_limits(
      xmin = 0,
      xmax = image_info$image_width,
      ymin = 0,
      ymax = image_info$image_height,
      padding_fraction = panel_padding_fraction_no_image
    )
  }

  list(
    plot_data = plot_data,
    frame_limits = frame_limits
  )
}


create_no_image_plot <- function(sample_id) {

  geometry <- prepare_no_image_geometry(sample_id)
  plot_data <- geometry$plot_data
  frame_limits <- geometry$frame_limits

  title_text <- make_plot_title_line(plot_data)
  subtitle_text <- make_plot_subtitle_line(plot_data)

  ggplot(
    plot_data,
    aes(
      x = x_use,
      y = y_use,
      colour = logNormalized_expression
    )
  ) +
    geom_point(
      size = point_size_no_image,
      shape = 16,
      stroke = 0
    ) +
    scale_x_continuous(
      limits = frame_limits$x_limits,
      expand = c(0, 0)
    ) +
    scale_y_reverse(
      limits = rev(frame_limits$y_limits),
      expand = c(0, 0)
    ) +
    coord_fixed() +
    scale_colour_gradientn(
      colours = palette_colors,
      values = palette_values,
      limits = c(common_colour_min, common_colour_max),
      oob = scales::squish,
      na.value = "#D9D9D9"
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      colour = paste0(
        target_gene,
        "\nlog1p(count/UMI × 10,000)"
      )
    ) +
    theme_void(base_family = "DejaVu Sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = panel_border_linewidth
      ),
      plot.title = element_text(
        size = 12.5,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 7.5,
        hjust = 0.5,
        lineheight = 1.08,
        margin = margin(b = 6)
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    ) +
    guides(
      colour = guide_colourbar(
        title.position = "top",
        direction = "horizontal",
        barwidth = grid::unit(42, "mm"),
        barheight = grid::unit(4.5, "mm")
      )
    )
}


create_with_image_plot <- function(sample_id) {

  plot_data <- sample_plot_data[[sample_id]]
  image_info <- image_info_list[[sample_id]]

  if (is.null(image_info)) {
    stop("Missing image information for sample: ", sample_id)
  }

  image_array <- png::readPNG(image_info$image_file)
  image_height <- dim(image_array)[1]
  image_width <- dim(image_array)[2]
  scale_factor <- image_info$scale_factor

  plot_data <- plot_data |>
    mutate(
      x_plot = pxl_col_in_fullres * scale_factor,
      y_plot = pxl_row_in_fullres * scale_factor
    )

  image_limits <- compute_square_limits(
    xmin = 0,
    xmax = image_width,
    ymin = 0,
    ymax = image_height,
    padding_fraction = panel_padding_fraction_with_image
  )

  title_text <- make_plot_title_line(plot_data)
  subtitle_text <- make_plot_subtitle_line(plot_data)

  ggplot() +
    annotation_raster(
      raster = image_array,
      xmin = 0,
      xmax = image_width,
      ymin = image_height,
      ymax = 0
    ) +
    geom_point(
      data = plot_data,
      aes(
        x = x_plot,
        y = y_plot,
        colour = logNormalized_expression
      ),
      size = point_size_with_image,
      shape = 16,
      stroke = 0
    ) +
    scale_x_continuous(
      limits = image_limits$x_limits,
      expand = c(0, 0)
    ) +
    scale_y_reverse(
      limits = rev(image_limits$y_limits),
      expand = c(0, 0)
    ) +
    coord_fixed() +
    scale_colour_gradientn(
      colours = palette_colors,
      values = palette_values,
      limits = c(common_colour_min, common_colour_max),
      oob = scales::squish,
      na.value = "#D9D9D9"
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      colour = paste0(
        target_gene,
        "\nlog1p(count/UMI × 10,000)"
      )
    ) +
    theme_void(base_family = "DejaVu Sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = panel_border_linewidth
      ),
      plot.title = element_text(
        size = 12.5,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 7.5,
        hjust = 0.5,
        lineheight = 1.08,
        margin = margin(b = 6)
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    ) +
    guides(
      colour = guide_colourbar(
        title.position = "top",
        direction = "horizontal",
        barwidth = grid::unit(42, "mm"),
        barheight = grid::unit(4.5, "mm")
      )
    )
}


create_gene_plot <- function(sample_id) {
  if (show_histology_image) {
    create_with_image_plot(sample_id)
  } else {
    create_no_image_plot(sample_id)
  }
}


create_no_image_mean_fill_plot <- function(sample_id) {

  geometry <- prepare_no_image_geometry(sample_id)
  plot_data <- geometry$plot_data
  frame_limits <- geometry$frame_limits

  plot_data <- plot_data |>
    mutate(
      sample_mean_fill_value = mean(logNormalized_expression, na.rm = TRUE)
    )

  title_text <- make_plot_title_line(plot_data)
  subtitle_text <- make_mean_fill_subtitle_line(plot_data)

  ggplot(
    plot_data,
    aes(
      x = x_use,
      y = y_use,
      colour = sample_mean_fill_value
    )
  ) +
    geom_point(
      size = point_size_no_image,
      shape = 16,
      stroke = 0
    ) +
    scale_x_continuous(
      limits = frame_limits$x_limits,
      expand = c(0, 0)
    ) +
    scale_y_reverse(
      limits = rev(frame_limits$y_limits),
      expand = c(0, 0)
    ) +
    coord_fixed() +
    scale_colour_gradientn(
      colours = palette_colors,
      values = palette_values,
      limits = c(sample_mean_colour_min, sample_mean_colour_max),
      oob = scales::squish,
      na.value = "#D9D9D9"
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      colour = paste0(
        target_gene,
        " sample mean\nlog1p(count/UMI × 10,000)"
      )
    ) +
    theme_void(base_family = "DejaVu Sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = panel_border_linewidth
      ),
      plot.title = element_text(
        size = 12.5,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 7.5,
        hjust = 0.5,
        lineheight = 1.08,
        margin = margin(b = 6)
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    ) +
    guides(
      colour = guide_colourbar(
        title.position = "top",
        direction = "horizontal",
        barwidth = grid::unit(42, "mm"),
        barheight = grid::unit(4.5, "mm")
      )
    )
}


create_with_image_mean_fill_plot <- function(sample_id) {

  plot_data <- sample_plot_data[[sample_id]]
  image_info <- image_info_list[[sample_id]]

  if (is.null(image_info)) {
    stop("Missing image information for sample: ", sample_id)
  }

  image_array <- png::readPNG(image_info$image_file)
  image_height <- dim(image_array)[1]
  image_width <- dim(image_array)[2]
  scale_factor <- image_info$scale_factor

  plot_data <- plot_data |>
    mutate(
      x_plot = pxl_col_in_fullres * scale_factor,
      y_plot = pxl_row_in_fullres * scale_factor,
      sample_mean_fill_value = mean(logNormalized_expression, na.rm = TRUE)
    )

  image_limits <- compute_square_limits(
    xmin = 0,
    xmax = image_width,
    ymin = 0,
    ymax = image_height,
    padding_fraction = panel_padding_fraction_with_image
  )

  title_text <- make_plot_title_line(plot_data)
  subtitle_text <- make_mean_fill_subtitle_line(plot_data)

  ggplot() +
    annotation_raster(
      raster = image_array,
      xmin = 0,
      xmax = image_width,
      ymin = image_height,
      ymax = 0
    ) +
    geom_point(
      data = plot_data,
      aes(
        x = x_plot,
        y = y_plot,
        colour = sample_mean_fill_value
      ),
      size = point_size_with_image,
      shape = 16,
      stroke = 0
    ) +
    scale_x_continuous(
      limits = image_limits$x_limits,
      expand = c(0, 0)
    ) +
    scale_y_reverse(
      limits = rev(image_limits$y_limits),
      expand = c(0, 0)
    ) +
    coord_fixed() +
    scale_colour_gradientn(
      colours = palette_colors,
      values = palette_values,
      limits = c(sample_mean_colour_min, sample_mean_colour_max),
      oob = scales::squish,
      na.value = "#D9D9D9"
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      colour = paste0(
        target_gene,
        " sample mean\nlog1p(count/UMI × 10,000)"
      )
    ) +
    theme_void(base_family = "DejaVu Sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = panel_border_linewidth
      ),
      plot.title = element_text(
        size = 12.5,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 7.5,
        hjust = 0.5,
        lineheight = 1.08,
        margin = margin(b = 6)
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    ) +
    guides(
      colour = guide_colourbar(
        title.position = "top",
        direction = "horizontal",
        barwidth = grid::unit(42, "mm"),
        barheight = grid::unit(4.5, "mm")
      )
    )
}


create_gene_mean_fill_plot <- function(sample_id) {
  if (show_histology_image) {
    create_with_image_mean_fill_plot(sample_id)
  } else {
    create_no_image_mean_fill_plot(sample_id)
  }
}


# ==============================================================================
# 10. Build spatial plots
# ==============================================================================

gene_plots <- lapply(selected_samples, create_gene_plot)
names(gene_plots) <- selected_samples

gene_mean_fill_plots <- lapply(selected_samples, create_gene_mean_fill_plot)
names(gene_mean_fill_plots) <- selected_samples


# ==============================================================================
# 11. Arrange plots into 4 requested columns
# ==============================================================================

column_order <- c(
  "male_neurotypical",
  "male_asd",
  "female_neurotypical",
  "female_asd"
)

samples_by_column <- lapply(
  column_order,
  function(column_key) {
    sample_table |>
      filter(column_key == !!column_key) |>
      pull(sample_ID)
  }
)

names(samples_by_column) <- column_order

max_rows <- max(lengths(samples_by_column))

message("")
message("Samples by output column:")
for (column_key in column_order) {
  message(
    column_key,
    ": ",
    paste(samples_by_column[[column_key]], collapse = ", ")
  )
}

pad_plots_by_column <- function(plot_list_named) {
  out <- lapply(
    column_order,
    function(column_key) {
      sample_ids <- samples_by_column[[column_key]]
      plot_list <- lapply(sample_ids, function(sid) plot_list_named[[sid]])

      if (length(plot_list) < max_rows) {
        plot_list <- c(
          plot_list,
          rep(list(patchwork::plot_spacer()), max_rows - length(plot_list))
        )
      }

      plot_list
    }
  )

  names(out) <- column_order
  out
}

interleave_padded_plots <- function(plots_by_column_padded) {
  interleaved_plots <- list()

  for (row_i in seq_len(max_rows)) {
    for (column_key in column_order) {
      interleaved_plots[[length(interleaved_plots) + 1]] <-
        plots_by_column_padded[[column_key]][[row_i]]
    }
  }

  interleaved_plots
}

plots_by_column_padded <- pad_plots_by_column(gene_plots)
mean_fill_plots_by_column_padded <- pad_plots_by_column(gene_mean_fill_plots)

interleaved_plots <- interleave_padded_plots(plots_by_column_padded)
interleaved_mean_fill_plots <- interleave_padded_plots(mean_fill_plots_by_column_padded)

layout_description <- paste(
  "Columns:",
  "1 = Male neurotypical",
  "2 = Male ASD",
  "3 = Female neurotypical",
  "4 = Female ASD",
  sep = " | "
)

subtitle_mode <- if (show_histology_image) {
  "Histology image shown"
} else {
  "No histology image; true Space Ranger pixel-coordinate frame preserved across samples"
}

combined_plot <- wrap_plots(
  interleaved_plots,
  ncol = number_of_columns,
  guides = "collect"
) +
  plot_annotation(
    title = paste0(
      target_gene,
      " spatial expression in 16 maternal FMT Visium samples"
    ),
    subtitle = paste0(
      subtitle_mode,
      " | ",
      layout_description,
      " | common colour range: 0–",
      round(common_colour_max, 3),
      " | upper limit = ",
      upper_colour_quantile * 100,
      "th percentile of positive values"
    ),
    theme = theme(
      text = element_text(family = "DejaVu Sans"),
      plot.title = element_text(
        size = 20,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 10,
        hjust = 0.5,
        lineheight = 1.04
      ),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center"
  )

mean_fill_combined_plot <- wrap_plots(
  interleaved_mean_fill_plots,
  ncol = number_of_columns,
  guides = "collect"
) +
  plot_annotation(
    title = paste0(
      target_gene,
      " sample-mean filled spatial maps in 16 maternal FMT Visium samples"
    ),
    subtitle = paste0(
      subtitle_mode,
      " | ",
      layout_description,
      " | every tissue spot in a sample is filled with that sample's mean expression",
      " | common colour range: 0–",
      round(sample_mean_colour_max, 3)
    ),
    theme = theme(
      text = element_text(family = "DejaVu Sans"),
      plot.title = element_text(
        size = 20,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 10,
        hjust = 0.5,
        lineheight = 1.04
      ),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center"
  )


# ==============================================================================
# 12. Create and save the four-group barplot and PDFs
# ==============================================================================

barplot_group_data <- group_summary_table |>
  mutate(
    group_label = factor(
      group_label,
      levels = unname(barplot_group_labels[barplot_group_order])
    ),
    ymin = pmax(
      0,
      mean_sampleMean_logNormalized - se_sampleMean_logNormalized
    ),
    ymax = mean_sampleMean_logNormalized + se_sampleMean_logNormalized
  )

barplot_sample_data <- barplot_sample_data |>
  mutate(
    group_label = factor(
      as.character(group_label),
      levels = unname(barplot_group_labels[barplot_group_order])
    )
  )

barplot_axis_labels <- setNames(
  nm = group_summary_table$group_label,
  object = paste0(
    unname(barplot_axis_labels_base[group_summary_table$group_label]),
    "\n",
    format(
      round(group_summary_table$mean_percent_target_positive_spots, 1),
      nsmall = 1
    ),
    "% ",
    target_gene,
    "+ spots"
  )
)

group_barplot <- ggplot() +
  geom_col(
    data = barplot_group_data,
    aes(
      x = group_label,
      y = mean_sampleMean_logNormalized
    ),
    width = barplot_bar_width,
    fill = NA,
    colour = "black",
    linewidth = 0.9
  ) +
  geom_errorbar(
    data = barplot_group_data,
    aes(
      x = group_label,
      ymin = ymin,
      ymax = ymax
    ),
    width = barplot_errorbar_width,
    linewidth = 0.8
  ) +
  geom_point(
    data = barplot_sample_data,
    aes(
      x = group_label,
      y = mean_target_logNormalized_all_spots
    ),
    position = position_jitter(
      width = barplot_jitter_width,
      height = 0,
      seed = 123
    ),
    shape = 21,
    size = barplot_point_size,
    stroke = 0.8,
    fill = "#B30000",
    colour = "black"
  ) +
  scale_x_discrete(
    labels = barplot_axis_labels
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    title = paste0(
      target_gene,
      " expression across maternal FMT groups"
    ),
    subtitle = paste0(
      "Empty bars: group mean | dots: sample means | whiskers: mean ± SE",
      "\nBelow each group: mean % of tissue spots with detectable ",
      target_gene,
      " expression"
    ),
    x = NULL,
    y = paste0(
      "Mean ",
      target_gene,
      " log-normalized expression per sample",
      "\nlog1p(count / total UMI × 10,000)"
    )
  ) +
  theme_classic(base_family = "DejaVu Sans") +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5,
      lineheight = 1.08,
      margin = margin(b = 12)
    ),
    axis.title.y = element_text(
      size = 11,
      face = "bold",
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 10,
      face = "bold",
      lineheight = 0.95,
      margin = margin(t = 6)
    ),
    axis.text.y = element_text(size = 9),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    plot.margin = margin(t = 12, r = 16, b = 12, l = 16)
  )

if (capabilities("cairo")) {
  pdf_device <- grDevices::cairo_pdf
} else {
  pdf_device <- grDevices::pdf
  warning("Cairo is unavailable; using standard pdf device.")
}

ggsave(
  filename = output_spatial_pdf,
  plot = combined_plot,
  device = pdf_device,
  width = combined_pdf_width_inches,
  height = combined_pdf_height_inches,
  units = "in",
  limitsize = FALSE,
  bg = "white"
)

ggsave(
  filename = output_barplot_pdf,
  plot = group_barplot,
  device = pdf_device,
  width = barplot_pdf_width_inches,
  height = barplot_pdf_height_inches,
  units = "in",
  limitsize = FALSE,
  bg = "white"
)

ggsave(
  filename = output_mean_fill_pdf,
  plot = mean_fill_combined_plot,
  device = pdf_device,
  width = mean_fill_pdf_width_inches,
  height = mean_fill_pdf_height_inches,
  units = "in",
  limitsize = FALSE,
  bg = "white"
)

required_output_files <- c(
  output_spatial_pdf,
  output_barplot_pdf,
  output_mean_fill_pdf,
  output_sample_summary_tsv,
  output_group_summary_tsv
)

missing_output_files <- required_output_files[
  !file.exists(required_output_files)
]

if (length(missing_output_files) > 0) {
  stop(
    "The following output files were not created:\n",
    paste(missing_output_files, collapse = "\n")
  )
}


# ==============================================================================
# 13. Final report
# ==============================================================================

spatial_pdf_size_mb <- file.info(output_spatial_pdf)$size / 1024^2
barplot_pdf_size_mb <- file.info(output_barplot_pdf)$size / 1024^2
mean_fill_pdf_size_mb <- file.info(output_mean_fill_pdf)$size / 1024^2

message("")
message("============================================================")
message("Analysis completed successfully.")
message("============================================================")
message("Target gene: ", target_gene)
message("Show histology image: ", show_histology_image)
message("Rotate no-image mode by 90 degrees: ", rotate_no_image_90)
message("Gene output directory: ", normalizePath(gene_output_dir, mustWork = TRUE))
message("Spatial-map PDF: ", normalizePath(output_spatial_pdf, mustWork = TRUE))
message("Spatial-map PDF size: ", round(spatial_pdf_size_mb, 2), " MB")
message("Group barplot PDF: ", normalizePath(output_barplot_pdf, mustWork = TRUE))
message("Group barplot PDF size: ", round(barplot_pdf_size_mb, 2), " MB")
message("Mean-fill spatial-map PDF: ", normalizePath(output_mean_fill_pdf, mustWork = TRUE))
message("Mean-fill spatial-map PDF size: ", round(mean_fill_pdf_size_mb, 2), " MB")
message(
  "Sample summary TSV: ",
  normalizePath(output_sample_summary_tsv, mustWork = TRUE)
)
message(
  "Group summary TSV: ",
  normalizePath(output_group_summary_tsv, mustWork = TRUE)
)
message("============================================================")
