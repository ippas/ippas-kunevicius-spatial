#!/usr/bin/env Rscript

# ==============================================================================
# 05_multigene_spatialMaps_groupBarplots_mainDonorGroup_FIXED_v2.R
#
# Definitive fixed multi-gene version. Prevents dplyr data-mask shadowing of
# target_gene and uses robust feature-name matching.
# Generates the same five outputs as the single-gene v4 script for all genes
# listed below. Each gene is written directly into:
#
#   <output_root_dir>/<GENE>/
#
# The 16 Space Ranger matrices are loaded only once. During the gene loop the
# console reports progress as [current/total], elapsed time, completed genes,
# skipped complete genes and failures.
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

script_version <- "2026-07-28_FIXED_v2"


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
  "pseudobulk_geneCounts_edgeR_sexInteraction",
  "proteinCodingGenes_autosomalAndXChromosome",
  "geneExpression_visualizationOnSlideAndBarplot",
  "mainDonorGroup"
)

# FALSE = no histology image
# TRUE  = histology image in the background
show_histology_image <- FALSE

# Used only when show_histology_image is FALSE
rotate_no_image_90 <- FALSE

normalization_scale_factor <- 10000
upper_colour_quantile <- 0.99

# A rerun skips a gene only when all five expected outputs already exist.
skip_completed_genes <- TRUE

# TRUE = report a failed gene and continue with the next gene.
continue_after_gene_error <- TRUE

run_status_file <- file.path(
  output_root_dir,
  "05_multigene_run_status.tsv"
)

gene_availability_file <- file.path(
  output_root_dir,
  "05_targetGene_availability.tsv"
)

target_genes <- c(
  "Nrp2", "Nppc", "Tmsb10b", "Gpc3", "Syndig1l", "Lrrc8d", "Rtn2", "Nectin1",
  "Lgals3bp", "Smox", "Ccdc152", "Mid1", "Phactr2", "Pdpn", "Arhgef40", "Sirt6",
  "Dbndd2", "Thbs1", "Tgfbr3", "Ksr2", "Shc4", "Trpm3", "Ier2", "Nr5a1",
  "Kcne5", "Hes5", "Sh3rf1", "Sap30l", "Ltbp1", "Wfikkn2", "Gpc4", "Ccdc85b",
  "Hdhd5", "Pcsk1n", "Crybb1", "Them6", "Vcpkmt", "Lpl", "Cd36", "Entpd7",
  "Lrrc10b", "Slc13a4", "Slc22a5", "Carmil2", "Prxl2b", "Hcrt", "Itga4", "Ogn",
  "Qdpr", "Ulk4", "Cryab", "Supt16", "Cd164l2", "Grid1", "Vps72", "Dab2",
  "Dusp1", "Tcap", "Fn3krp", "Pfdn2", "Acvr1c", "Cgnl1", "Actb", "Apoa1",
  "Tsku", "Efcab12", "Kctd12", "Pcolce", "Kics2", "Dcn", "Rraga", "Zar1",
  "D17H6S53E", "Gpc1", "Fezf1", "Tmem87a", "Vgll3", "Rnd2", "Calcrl", "Prrc2c",
  "Bag1", "Teddm2", "Pacs2", "Ppp1r14b", "Cavin1", "Cped1", "Mmrn2", "Rep15",
  "Gal3st1", "Adamtsl5", "Dhdh", "Thrsp", "Gab3", "St3gal1", "Hsf2", "Ccdc73",
  "Samhd1", "Inmt", "Dynlt1b", "Bgn", "Nme7", "Timmdc1", "Gal3st3", "Rpsa",
  "Scp2", "Pter", "Xk", "Arnt", "Nid2", "Aldh1a3", "Syngap1", "Il17rd",
  "Tuba1c", "Ankmy1", "Fam98b", "Adam23", "Loxl2", "Tgfbi", "Map3k15", "Cacna1s",
  "Slc38a6", "Pcdhb7", "Lum", "F13a1", "Col8a2", "Tsc22d4", "Tm4sf5", "Lsm14b",
  "Fth1", "Pdlim2", "Gm11837", "Enho", "Slc6a20a", "Npr3", "Krt73", "Foxk1",
  "Lipg", "Lmnb1", "Ftl1", "Myl6b", "Septin10", "Maff", "Car12", "Notch2",
  "Myadm", "Gbp2", "Cenpf", "Abhd2", "Myh3", "Mrc1", "Sertad3", "Pam",
  "Cyp4f17", "Cdh19", "Ccbe1", "Cx3cr1", "Tmie", "Cd248", "Maml3", "Ppp1r3b",
  "Rps4x", "Tssk6", "Jag1", "Rps2", "Prg4", "Ogdhl", "Rps3a1", "Frrs1",
  "Zim1", "Frmpd3", "Nid1", "Rgs9", "Snx22", "Slc16a8", "Aox3", "Zdhhc14",
  "Pltp", "Cdkn1b", "Hspg2", "Wdcp", "Nfkb2", "Nkx6-2", "Tmem275", "Drd2",
  "Bco2", "Hyou1", "Ptgis", "Slc16a12", "Gpbp1l1", "Fli1", "Hrct1", "Polh",
  "Myh9", "Depp1", "Itpr3", "Khsrp", "Ccn3", "Hmgn2", "Kremen2", "Pcdhb13",
  "Pik3c2b", "Sun2", "Nat2", "Plekhf1", "Lepr", "Myof", "Tbx1", "Penk",
  "Inhba", "Sptbn5", "Ppwd1", "Rbm45", "Gemin4", "Prlr", "Sp140", "Krt9",
  "Ggcx", "Zfp619", "Col8a1", "Klhl33", "H2-Ab1", "Itpripl1", "Col1a2", "Tfap4",
  "Mycbpap", "Pabpn1", "Aass", "Itgad", "Dact2", "1190005I06Rik", "Gcnt1", "Zar1l",
  "Lrmda", "Fos", "Steap1", "Crb2", "A2m", "Col9a3", "Ube2srt", "Aqp1",
  "Rnls", "Sema7a", "Cfap74", "Slc39a4", "Ceacam1", "Mpzl2", "Fhip1a", "Tbx3",
  "Hsd17b11", "Snai1", "Lbhd2", "Ntn4", "Ace", "Zfp691", "Wnt2", "Meis1",
  "Slc16a9", "Edn3", "Klf8", "Slc24a5", "Glis3", "Isl1", "Ppp1r1b", "Mapk4",
  "Slc45a3"
)

if (length(target_genes) == 0) {
  stop("target_genes is empty.")
}

if (anyDuplicated(target_genes) > 0) {
  duplicated_genes <- unique(target_genes[duplicated(target_genes)])
  stop("Duplicated target genes: ", paste(duplicated_genes, collapse = ", "))
}


# ==============================================================================
# 3. Plot settings
# ==============================================================================

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


# ==============================================================================
# 4. Samples
# ==============================================================================

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

column_order <- c(
  "male_neurotypical",
  "male_asd",
  "female_neurotypical",
  "female_asd"
)

column_titles <- c(
  male_neurotypical = "Male neurotypical",
  male_asd = "Male ASD",
  female_neurotypical = "Female neurotypical",
  female_asd = "Female ASD"
)

barplot_group_labels <- column_titles

barplot_axis_labels_base <- c(
  "Male neurotypical" = "Male\nneurotypical",
  "Male ASD" = "Male\nASD",
  "Female neurotypical" = "Female\nneurotypical",
  "Female ASD" = "Female\nASD"
)


# ==============================================================================
# 5. General helper functions
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


match_requested_genes <- function(feature_names, requested_genes) {

  canonicalize_feature_name <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("\\.[0-9]+$", "", x)
    x <- sub("-[0-9]+$", "", x)
    x
  }

  feature_names <- as.character(feature_names)
  feature_names_lower <- tolower(feature_names)
  feature_names_canonical <- canonicalize_feature_name(feature_names)
  feature_names_canonical_lower <- tolower(feature_names_canonical)

  matched_index <- rep(NA_integer_, length(requested_genes))
  matched_name <- rep(NA_character_, length(requested_genes))
  match_type <- rep("not_found", length(requested_genes))

  for (gene_i in seq_along(requested_genes)) {

    requested_gene <- trimws(as.character(requested_genes[[gene_i]]))
    requested_gene_lower <- tolower(requested_gene)
    requested_gene_canonical <- canonicalize_feature_name(requested_gene)
    requested_gene_canonical_lower <- tolower(requested_gene_canonical)

    candidate_sets <- list(
      exact = which(feature_names == requested_gene),
      case_insensitive = which(feature_names_lower == requested_gene_lower),
      canonical_exact = which(
        feature_names_canonical == requested_gene_canonical
      ),
      canonical_case_insensitive = which(
        feature_names_canonical_lower == requested_gene_canonical_lower
      )
    )

    matched <- FALSE

    for (current_match_type in names(candidate_sets)) {
      current_indices <- unique(candidate_sets[[current_match_type]])

      if (length(current_indices) == 1) {
        matched_index[[gene_i]] <- current_indices[[1]]
        matched_name[[gene_i]] <- feature_names[[current_indices[[1]]]]
        match_type[[gene_i]] <- current_match_type
        matched <- TRUE
        break
      }

      if (length(current_indices) > 1) {
        match_type[[gene_i]] <- "ambiguous"
        matched <- TRUE
        break
      }
    }

    if (!matched) {
      match_type[[gene_i]] <- "not_found"
    }
  }

  tibble(
    requested_gene = as.character(requested_genes),
    matched_index = matched_index,
    matched_name = matched_name,
    match_type = match_type
  )
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
    "\nExpected ASD/autism or neurotypical/control."
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


compute_square_limits <- function(
  xmin,
  xmax,
  ymin,
  ymax,
  padding_fraction = 0.03
) {

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

  if (
    !is.numeric(scale_factor) ||
      length(scale_factor) != 1 ||
      is.na(scale_factor)
  ) {
    stop("Could not determine image scale factor for: ", image_file)
  }

  image_array <- png::readPNG(image_file)
  image_height <- dim(image_array)[1]
  image_width <- dim(image_array)[2]

  list(
    image_file = image_file,
    image_array = if (show_histology_image) image_array else NULL,
    scale_factor = scale_factor,
    image_width = image_width,
    image_height = image_height
  )
}


format_number <- function(x, digits = 3) {
  format(
    round(x, digits),
    nsmall = digits,
    trim = TRUE,
    scientific = FALSE
  )
}


write_tsv <- function(data, filename) {
  write.table(
    data,
    file = filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}


# ==============================================================================
# 6. Validate input and read metadata
# ==============================================================================

message("============================================================")
message("SCRIPT VERSION: ", script_version)
message("IMPORTANT: legacy vectorIndex bug is removed in this version.")
message("============================================================")

if (!dir.exists(project_dir)) {
  stop("Project directory does not exist: ", project_dir)
}

if (!dir.exists(path_to_data)) {
  stop("Space Ranger data directory does not exist: ", path_to_data)
}

if (!file.exists(metadata_file)) {
  stop("Metadata file does not exist: ", metadata_file)
}

dir.create(
  output_root_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

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

missing_selected_samples <- setdiff(
  selected_samples,
  sample_table$sample_ID
)

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

observed_counts_by_column <- table(sample_table$column_key)

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

samples_by_column <- lapply(
  column_order,
  function(current_column_key) {
    sample_table |>
      filter(column_key == current_column_key) |>
      pull(sample_ID)
  }
)

names(samples_by_column) <- column_order
max_rows <- max(lengths(samples_by_column))

mode_label <- if (show_histology_image) {
  "withHistology"
} else if (rotate_no_image_90) {
  "noHistology_pixelFrame_rotated90"
} else {
  "noHistology_pixelFrame"
}

if (capabilities("cairo")) {
  pdf_device <- grDevices::cairo_pdf
} else {
  pdf_device <- grDevices::pdf
  warning("Cairo is unavailable; using standard pdf device.")
}


# ==============================================================================
# 7. Load each sample once and cache all requested genes
# ==============================================================================

message("")
message("============================================================")
message("Loading and caching 16 samples")
message("Requested genes: ", length(target_genes))
message("Output root: ", output_root_dir)
message("============================================================")

sample_cache <- vector("list", length(selected_samples))
names(sample_cache) <- selected_samples

availability_rows <- vector(
  "list",
  length(selected_samples)
)

for (sample_i in seq_along(selected_samples)) {

  sample_id <- selected_samples[[sample_i]]

  message(
    sprintf(
      "[sample %02d/%02d] Loading %s",
      sample_i,
      length(selected_samples),
      sample_id
    )
  )

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
    stop(
      "No shared barcodes between coordinates and matrix for sample: ",
      sample_id
    )
  }

  coordinates <- coordinates |>
    filter(barcode %in% common_barcodes) |>
    arrange(match(barcode, common_barcodes))

  counts <- counts[, coordinates$barcode, drop = FALSE]

  if (!identical(colnames(counts), coordinates$barcode)) {
    stop("Barcode order mismatch for sample: ", sample_id)
  }

  gene_matches <- match_requested_genes(
    feature_names = rownames(counts),
    requested_genes = target_genes
  )

  ambiguous_genes <- gene_matches |>
    filter(match_type == "ambiguous") |>
    pull(requested_gene)

  if (length(ambiguous_genes) > 0) {
    stop(
      "Ambiguous feature matches in sample ",
      sample_id,
      ": ",
      paste(ambiguous_genes, collapse = ", ")
    )
  }

  found_matches <- gene_matches |>
    filter(!is.na(matched_index))

  target_counts <- counts[
    found_matches$matched_index,
    ,
    drop = FALSE
  ]

  rownames(target_counts) <- found_matches$requested_gene

  matched_gene_names <- setNames(
    found_matches$matched_name,
    found_matches$requested_gene
  )

  total_umi_per_spot <- as.numeric(Matrix::colSums(counts))

  sample_metadata <- sample_table |>
    filter(sample_ID == sample_id)

  image_info <- choose_image_file_and_scale(spatial_dir)

  coordinates <- coordinates |>
    mutate(
      sample_ID = sample_id,
      donor_group_raw = sample_metadata$donor_group_raw[[1]],
      sex_raw = sample_metadata$sex_raw[[1]],
      group_std = sample_metadata$group_std[[1]],
      sex_std = sample_metadata$sex_std[[1]],
      column_key = sample_metadata$column_key[[1]],
      x_plot = pxl_col_in_fullres * image_info$scale_factor,
      y_plot = pxl_row_in_fullres * image_info$scale_factor
    )

  if (rotate_no_image_90) {

    coordinates <- coordinates |>
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

    coordinates <- coordinates |>
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

  image_limits <- compute_square_limits(
    xmin = 0,
    xmax = image_info$image_width,
    ymin = 0,
    ymax = image_info$image_height,
    padding_fraction = panel_padding_fraction_with_image
  )

  sample_cache[[sample_id]] <- list(
    coordinates = coordinates,
    target_counts = target_counts,
    matched_gene_names = matched_gene_names,
    total_umi_per_spot = total_umi_per_spot,
    matrix_file = matrix_file,
    positions_file = positions_file,
    image_info = image_info,
    frame_limits = frame_limits,
    image_limits = image_limits
  )

  availability_rows[[sample_i]] <- gene_matches |>
    transmute(
      sample_ID = sample_id,
      requested_gene = requested_gene,
      matched_gene = matched_name,
      match_type = match_type,
      found = !is.na(matched_index)
    )

  rm(
    counts,
    coordinates,
    target_counts,
    total_umi_per_spot,
    gene_matches,
    found_matches
  )

  invisible(gc(verbose = FALSE))

  message(
    sprintf(
      "[sample %02d/%02d] Cached %s",
      sample_i,
      length(selected_samples),
      sample_id
    )
  )
}

availability_long <- bind_rows(availability_rows)

availability_table <- availability_long |>
  group_by(requested_gene) |>
  summarise(
    n_samples_found = sum(found),
    n_samples_missing = sum(!found),
    missing_samples = paste(sample_ID[!found], collapse = ","),
    matched_feature_names = paste(
      sort(unique(na.omit(matched_gene))),
      collapse = ","
    ),
    .groups = "drop"
  ) |>
  mutate(
    requested_gene = factor(requested_gene, levels = target_genes)
  ) |>
  arrange(requested_gene) |>
  mutate(
    requested_gene = as.character(requested_gene)
  )

write_tsv(
  availability_table,
  gene_availability_file
)

message("")
message("Sample cache completed.")
message(
  "Genes found in all samples: ",
  sum(availability_table$n_samples_missing == 0),
  "/",
  nrow(availability_table)
)


# ==============================================================================
# 8. Gene-level data preparation
# ==============================================================================

build_gene_data <- function(target_gene) {

  availability_row <- availability_table |>
    filter(requested_gene == .env$target_gene)

  if (nrow(availability_row) != 1) {
    stop("Availability row missing for gene: ", target_gene)
  }

  if (availability_row$n_samples_missing[[1]] > 0) {
    stop(
      "Gene ",
      target_gene,
      " is missing in ",
      availability_row$n_samples_missing[[1]],
      " sample(s): ",
      availability_row$missing_samples[[1]]
    )
  }

  sample_plot_data <- vector(
    "list",
    length(selected_samples)
  )

  sample_summary <- vector(
    "list",
    length(selected_samples)
  )

  names(sample_plot_data) <- selected_samples
  names(sample_summary) <- selected_samples

  # Keep scalar gene names in variables that cannot be shadowed by
  # columns created later inside dplyr::mutate() / tibble().
  requested_gene_scalar <- as.character(target_gene)[[1]]

  for (sample_id in selected_samples) {

    current_cache <- sample_cache[[sample_id]]

    cached_gene_row <- match(
      requested_gene_scalar,
      rownames(current_cache$target_counts)
    )

    if (is.na(cached_gene_row)) {
      stop(
        "Cached target-count row is missing for ",
        requested_gene_scalar,
        " in sample ",
        sample_id
      )
    }

    matched_name_position <- match(
      requested_gene_scalar,
      names(current_cache$matched_gene_names)
    )

    if (is.na(matched_name_position)) {
      stop(
        "Cached matched-name entry is missing for ",
        requested_gene_scalar,
        " in sample ",
        sample_id
      )
    }

    matched_gene_scalar <- unname(
      current_cache$matched_gene_names[[matched_name_position]]
    )

    if (
      length(matched_gene_scalar) != 1 ||
        is.na(matched_gene_scalar) ||
        matched_gene_scalar == ""
    ) {
      stop(
        "Matched matrix gene name is invalid for ",
        requested_gene_scalar,
        " in sample ",
        sample_id
      )
    }

    target_raw_counts <- as.numeric(
      current_cache$target_counts[
        cached_gene_row,
        ,
        drop = TRUE
      ]
    )

    total_umi_per_spot <- current_cache$total_umi_per_spot

    log_normalized_expression <- numeric(
      length(total_umi_per_spot)
    )

    valid_spots <- total_umi_per_spot > 0

    log_normalized_expression[valid_spots] <- log1p(
      target_raw_counts[valid_spots] /
        total_umi_per_spot[valid_spots] *
        normalization_scale_factor
    )

    plot_data <- current_cache$coordinates |>
      mutate(
        target_gene = .env$requested_gene_scalar,
        target_gene_in_matrix = .env$matched_gene_scalar,
        target_raw_count = .env$target_raw_counts,
        total_UMI = .env$total_umi_per_spot,
        logNormalized_expression = .env$log_normalized_expression
      )

    sample_plot_data[[sample_id]] <- plot_data

    sample_summary[[sample_id]] <- tibble(
      sample_ID = sample_id,
      donor_group_raw = plot_data$donor_group_raw[[1]],
      sex_raw = plot_data$sex_raw[[1]],
      group_std = plot_data$group_std[[1]],
      sex_std = plot_data$sex_std[[1]],
      column_key = plot_data$column_key[[1]],
      matrix_file = current_cache$matrix_file,
      positions_file = current_cache$positions_file,
      target_gene_requested = requested_gene_scalar,
      target_gene_in_matrix = matched_gene_scalar,
      n_tissue_spots = nrow(plot_data),
      total_UMI_all_tissue_spots = sum(plot_data$total_UMI),
      median_total_UMI_per_tissue_spot = median(plot_data$total_UMI),
      n_target_positive_spots = sum(plot_data$target_raw_count > 0),
      percent_target_positive_spots =
        100 * mean(plot_data$target_raw_count > 0),
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
  }

  summary_table <- bind_rows(sample_summary) |>
    mutate(
      sample_ID = factor(sample_ID, levels = selected_samples)
    ) |>
    arrange(sample_ID) |>
    mutate(
      sample_ID = as.character(sample_ID)
    )

  barplot_sample_data <- summary_table |>
    mutate(
      group_key = factor(column_key, levels = column_order),
      group_label = factor(
        barplot_group_labels[as.character(group_key)],
        levels = unname(barplot_group_labels[column_order])
      )
    )

  group_summary_table <- barplot_sample_data |>
    group_by(group_key, group_label) |>
    summarise(
      n_samples = n(),
      total_tissue_spots = sum(n_tissue_spots),
      total_target_positive_spots = sum(n_target_positive_spots),
      pooled_percent_target_positive_spots =
        100 * sum(n_target_positive_spots) / sum(n_tissue_spots),
      mean_sampleMean_logNormalized =
        mean(mean_target_logNormalized_all_spots, na.rm = TRUE),
      sd_sampleMean_logNormalized =
        sd(mean_target_logNormalized_all_spots, na.rm = TRUE),
      se_sampleMean_logNormalized =
        sd(mean_target_logNormalized_all_spots, na.rm = TRUE) / sqrt(n()),
      mean_percent_target_positive_spots =
        mean(percent_target_positive_spots, na.rm = TRUE),
      sd_percent_target_positive_spots =
        sd(percent_target_positive_spots, na.rm = TRUE),
      se_percent_target_positive_spots =
        sd(percent_target_positive_spots, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) |>
    arrange(group_key) |>
    mutate(
      group_key = as.character(group_key),
      group_label = as.character(group_label)
    )

  all_plot_data <- bind_rows(sample_plot_data)

  positive_expression_values <- all_plot_data$logNormalized_expression[
    is.finite(all_plot_data$logNormalized_expression) &
      all_plot_data$logNormalized_expression > 0
  ]

  if (length(positive_expression_values) == 0) {
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
      common_colour_max <- max(
        positive_expression_values,
        na.rm = TRUE
      )
    }
  }

  sample_mean_values <-
    summary_table$mean_target_logNormalized_all_spots

  sample_mean_positive_values <- sample_mean_values[
    is.finite(sample_mean_values) &
      sample_mean_values > 0
  ]

  if (length(sample_mean_positive_values) == 0) {
    sample_mean_colour_max <- 1
  } else {
    sample_mean_colour_max <- max(
      sample_mean_positive_values,
      na.rm = TRUE
    )

    if (
      !is.finite(sample_mean_colour_max) ||
        sample_mean_colour_max <= 0
    ) {
      sample_mean_colour_max <- 1
    }
  }

  list(
    sample_plot_data = sample_plot_data,
    summary_table = summary_table,
    barplot_sample_data = barplot_sample_data,
    group_summary_table = group_summary_table,
    common_colour_max = common_colour_max,
    sample_mean_colour_max = sample_mean_colour_max
  )
}


# ==============================================================================
# 9. Plot functions
# ==============================================================================

make_standard_subtitle <- function(plot_data, target_gene) {

  n_positive <- sum(plot_data$target_raw_count > 0)
  percent_positive <- 100 * mean(plot_data$target_raw_count > 0)

  paste0(
    "Group: ",
    plot_data$group_std[[1]],
    " | Sex: ",
    plot_data$sex_std[[1]],
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
    format_number(
      mean(plot_data$logNormalized_expression, na.rm = TRUE)
    ),
    " ± ",
    format_number(
      sd(plot_data$logNormalized_expression, na.rm = TRUE)
    )
  )
}


make_mean_fill_subtitle <- function(plot_data, target_gene) {

  n_positive <- sum(plot_data$target_raw_count > 0)
  percent_positive <- 100 * mean(plot_data$target_raw_count > 0)
  sample_mean <- mean(
    plot_data$logNormalized_expression,
    na.rm = TRUE
  )

  paste0(
    "Group: ",
    plot_data$group_std[[1]],
    " | Sex: ",
    plot_data$sex_std[[1]],
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
    format_number(sample_mean)
  )
}


create_spatial_panel <- function(
  sample_id,
  plot_data,
  target_gene,
  fill_mode = c("per_spot", "sample_mean"),
  colour_max
) {

  fill_mode <- match.arg(fill_mode)
  current_cache <- sample_cache[[sample_id]]

  if (fill_mode == "per_spot") {
    plot_data$plot_value <- plot_data$logNormalized_expression
    subtitle_text <- make_standard_subtitle(
      plot_data,
      target_gene
    )
    legend_title <- paste0(
      target_gene,
      "\nlog1p(count/UMI × 10,000)"
    )
  } else {
    sample_mean <- mean(
      plot_data$logNormalized_expression,
      na.rm = TRUE
    )
    plot_data$plot_value <- sample_mean
    subtitle_text <- make_mean_fill_subtitle(
      plot_data,
      target_gene
    )
    legend_title <- paste0(
      target_gene,
      " sample mean\nlog1p(count/UMI × 10,000)"
    )
  }

  title_text <- paste0("Sample: ", sample_id)

  if (show_histology_image) {

    image_info <- current_cache$image_info
    image_limits <- current_cache$image_limits

    output_plot <- ggplot() +
      annotation_raster(
        raster = image_info$image_array,
        xmin = 0,
        xmax = image_info$image_width,
        ymin = image_info$image_height,
        ymax = 0
      ) +
      geom_point(
        data = plot_data,
        aes(
          x = x_plot,
          y = y_plot,
          colour = plot_value
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
      )

  } else {

    frame_limits <- current_cache$frame_limits

    output_plot <- ggplot(
      plot_data,
      aes(
        x = x_use,
        y = y_use,
        colour = plot_value
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
      )
  }

  output_plot +
    coord_fixed() +
    scale_colour_gradientn(
      colours = palette_colors,
      values = palette_values,
      limits = c(0, colour_max),
      oob = scales::squish,
      na.value = "#D9D9D9"
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      colour = legend_title
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


arrange_four_columns <- function(
  plot_list_named,
  plot_title,
  plot_subtitle
) {

  plots_by_column_padded <- lapply(
    column_order,
    function(current_column_key) {

      sample_ids <- samples_by_column[[current_column_key]]

      current_plots <- lapply(
        sample_ids,
        function(sample_id) {
          plot_list_named[[sample_id]]
        }
      )

      if (length(current_plots) < max_rows) {
        current_plots <- c(
          current_plots,
          rep(
            list(patchwork::plot_spacer()),
            max_rows - length(current_plots)
          )
        )
      }

      current_plots
    }
  )

  names(plots_by_column_padded) <- column_order

  interleaved_plots <- list()

  for (row_i in seq_len(max_rows)) {
    for (current_column_key in column_order) {
      interleaved_plots[[length(interleaved_plots) + 1]] <-
        plots_by_column_padded[[current_column_key]][[row_i]]
    }
  }

  wrap_plots(
    interleaved_plots,
    ncol = number_of_columns,
    guides = "collect"
  ) +
    plot_annotation(
      title = plot_title,
      subtitle = plot_subtitle,
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
}


create_group_barplot <- function(
  target_gene,
  barplot_sample_data,
  group_summary_table
) {

  barplot_group_data <- group_summary_table |>
    mutate(
      group_label = factor(
        group_label,
        levels = unname(barplot_group_labels[column_order])
      ),
      ymin = pmax(
        0,
        mean_sampleMean_logNormalized -
          se_sampleMean_logNormalized
      ),
      ymax = mean_sampleMean_logNormalized +
        se_sampleMean_logNormalized
    )

  barplot_sample_data <- barplot_sample_data |>
    mutate(
      group_label = factor(
        as.character(group_label),
        levels = unname(barplot_group_labels[column_order])
      )
    )

  barplot_axis_labels <- setNames(
    object = paste0(
      unname(
        barplot_axis_labels_base[
          group_summary_table$group_label
        ]
      ),
      "\n",
      format(
        round(
          group_summary_table$pooled_percent_target_positive_spots,
          1
        ),
        nsmall = 1
      ),
      "% ",
      target_gene,
      "+ spots"
    ),
    nm = group_summary_table$group_label
  )

  ggplot() +
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
        "Empty bars: group mean | dots: sample means | ",
        "whiskers: mean ± SE",
        "\nBelow each group: pooled % of tissue spots with ",
        "detectable ",
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
}


# ==============================================================================
# 10. One-gene output function
# ==============================================================================

get_gene_output_files <- function(target_gene) {

  gene_output_dir <- file.path(
    output_root_dir,
    target_gene
  )

  list(
    gene_output_dir = gene_output_dir,
    spatial_pdf = file.path(
      gene_output_dir,
      paste0(
        "05_",
        target_gene,
        "_logNormalized_spatialMaps_n16_",
        mode_label,
        "_fourColumnsSexGroup.pdf"
      )
    ),
    barplot_pdf = file.path(
      gene_output_dir,
      paste0(
        "05_",
        target_gene,
        "_sampleMean_groupBarplot_meanSE_",
        "withPercentPositive.pdf"
      )
    ),
    mean_fill_pdf = file.path(
      gene_output_dir,
      paste0(
        "05_",
        target_gene,
        "_sampleMeanFilled_spatialMaps_n16_",
        mode_label,
        "_fourColumnsSexGroup.pdf"
      )
    ),
    sample_summary_tsv = file.path(
      gene_output_dir,
      paste0(
        "05_",
        target_gene,
        "_logNormalized_spatialMaps_n16_",
        mode_label,
        "_summary.tsv"
      )
    ),
    group_summary_tsv = file.path(
      gene_output_dir,
      paste0(
        "05_",
        target_gene,
        "_sampleMean_group_summary.tsv"
      )
    )
  )
}


gene_outputs_complete <- function(output_files) {

  expected_files <- unlist(
    output_files[
      c(
        "spatial_pdf",
        "barplot_pdf",
        "mean_fill_pdf",
        "sample_summary_tsv",
        "group_summary_tsv"
      )
    ],
    use.names = FALSE
  )

  all(file.exists(expected_files))
}


process_one_gene <- function(target_gene) {

  output_files <- get_gene_output_files(target_gene)

  dir.create(
    output_files$gene_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  gene_data <- build_gene_data(target_gene)

  write_tsv(
    gene_data$summary_table,
    output_files$sample_summary_tsv
  )

  write_tsv(
    gene_data$group_summary_table,
    output_files$group_summary_tsv
  )

  standard_panels <- lapply(
    selected_samples,
    function(sample_id) {
      create_spatial_panel(
        sample_id = sample_id,
        plot_data =
          gene_data$sample_plot_data[[sample_id]],
        target_gene = target_gene,
        fill_mode = "per_spot",
        colour_max = gene_data$common_colour_max
      )
    }
  )

  names(standard_panels) <- selected_samples

  mean_fill_panels <- lapply(
    selected_samples,
    function(sample_id) {
      create_spatial_panel(
        sample_id = sample_id,
        plot_data =
          gene_data$sample_plot_data[[sample_id]],
        target_gene = target_gene,
        fill_mode = "sample_mean",
        colour_max = gene_data$sample_mean_colour_max
      )
    }
  )

  names(mean_fill_panels) <- selected_samples

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
    paste0(
      "No histology image; true Space Ranger ",
      "pixel-coordinate frame preserved across samples"
    )
  }

  standard_combined_plot <- arrange_four_columns(
    plot_list_named = standard_panels,
    plot_title = paste0(
      target_gene,
      " spatial expression in 16 maternal FMT Visium samples"
    ),
    plot_subtitle = paste0(
      subtitle_mode,
      " | ",
      layout_description,
      " | common colour range: 0–",
      round(gene_data$common_colour_max, 3),
      " | upper limit = ",
      upper_colour_quantile * 100,
      "th percentile of positive values"
    )
  )

  mean_fill_combined_plot <- arrange_four_columns(
    plot_list_named = mean_fill_panels,
    plot_title = paste0(
      target_gene,
      " sample-mean filled spatial maps in 16 maternal FMT Visium samples"
    ),
    plot_subtitle = paste0(
      subtitle_mode,
      " | ",
      layout_description,
      " | every tissue spot in a sample is filled with ",
      "that sample's mean expression",
      " | common colour range: 0–",
      round(gene_data$sample_mean_colour_max, 3)
    )
  )

  group_barplot <- create_group_barplot(
    target_gene = target_gene,
    barplot_sample_data = gene_data$barplot_sample_data,
    group_summary_table = gene_data$group_summary_table
  )

  ggsave(
    filename = output_files$spatial_pdf,
    plot = standard_combined_plot,
    device = pdf_device,
    width = combined_pdf_width_inches,
    height = combined_pdf_height_inches,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )

  ggsave(
    filename = output_files$barplot_pdf,
    plot = group_barplot,
    device = pdf_device,
    width = barplot_pdf_width_inches,
    height = barplot_pdf_height_inches,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )

  ggsave(
    filename = output_files$mean_fill_pdf,
    plot = mean_fill_combined_plot,
    device = pdf_device,
    width = mean_fill_pdf_width_inches,
    height = mean_fill_pdf_height_inches,
    units = "in",
    limitsize = FALSE,
    bg = "white"
  )

  if (!gene_outputs_complete(output_files)) {
    stop(
      "Not all five expected outputs were created for gene: ",
      target_gene
    )
  }

  rm(
    gene_data,
    standard_panels,
    mean_fill_panels,
    standard_combined_plot,
    mean_fill_combined_plot,
    group_barplot
  )

  invisible(gc(verbose = FALSE))

  output_files
}


# ==============================================================================
# 11. Main gene loop with progress reporting
# ==============================================================================

message("")
message("============================================================")
message("Starting multi-gene loop")
message("Genes to process: ", length(target_genes))
message("============================================================")

run_started_at <- Sys.time()
status_rows <- list()

for (gene_i in seq_along(target_genes)) {

  target_gene <- target_genes[[gene_i]]
  gene_started_at <- Sys.time()
  output_files <- get_gene_output_files(target_gene)

  message("")
  message(
    sprintf(
      "[%03d/%03d] START %s",
      gene_i,
      length(target_genes),
      target_gene
    )
  )

  current_status <- "completed"
  current_message <- ""
  output_directory <- output_files$gene_output_dir

  if (
    skip_completed_genes &&
      gene_outputs_complete(output_files)
  ) {

    current_status <- "skipped_complete"
    current_message <- "All five expected files already exist."

  } else {

    result <- tryCatch(
      {
        process_one_gene(target_gene)
      },
      error = function(error_condition) {
        error_condition
      }
    )

    if (inherits(result, "error")) {
      current_status <- "failed"
      current_message <- conditionMessage(result)
    }
  }

  gene_finished_at <- Sys.time()
  elapsed_seconds <- as.numeric(
    difftime(
      gene_finished_at,
      gene_started_at,
      units = "secs"
    )
  )

  status_rows[[length(status_rows) + 1]] <- tibble(
    order = gene_i,
    total_genes = length(target_genes),
    target_gene = target_gene,
    status = current_status,
    message = current_message,
    elapsed_seconds = round(elapsed_seconds, 2),
    output_directory = output_directory,
    started_at = format(
      gene_started_at,
      "%Y-%m-%d %H:%M:%S"
    ),
    finished_at = format(
      gene_finished_at,
      "%Y-%m-%d %H:%M:%S"
    )
  )

  current_status_table <- bind_rows(status_rows)

  write_tsv(
    current_status_table,
    run_status_file
  )

  n_completed <- sum(
    current_status_table$status == "completed"
  )

  n_skipped <- sum(
    current_status_table$status == "skipped_complete"
  )

  n_failed <- sum(
    current_status_table$status == "failed"
  )

  message(
    sprintf(
      paste0(
        "[%03d/%03d] END %s | status=%s | ",
        "%.1f s | completed=%d | skipped=%d | failed=%d"
      ),
      gene_i,
      length(target_genes),
      target_gene,
      current_status,
      elapsed_seconds,
      n_completed,
      n_skipped,
      n_failed
    )
  )

  if (current_status == "failed") {
    message("ERROR: ", current_message)

    if (!continue_after_gene_error) {
      stop(
        "Stopping after failure for gene ",
        target_gene,
        ": ",
        current_message
      )
    }
  }

  invisible(gc(verbose = FALSE))
}


# ==============================================================================
# 12. Final report
# ==============================================================================

final_status_table <- bind_rows(status_rows)

total_elapsed_minutes <- as.numeric(
  difftime(
    Sys.time(),
    run_started_at,
    units = "mins"
  )
)

n_completed <- sum(
  final_status_table$status == "completed"
)

n_skipped <- sum(
  final_status_table$status == "skipped_complete"
)

n_failed <- sum(
  final_status_table$status == "failed"
)

message("")
message("============================================================")
message("Multi-gene analysis finished")
message("============================================================")
message("Requested genes: ", length(target_genes))
message("Completed in this run: ", n_completed)
message("Skipped because complete: ", n_skipped)
message("Failed: ", n_failed)
message("Elapsed time: ", round(total_elapsed_minutes, 2), " min")
message("Output root: ", normalizePath(output_root_dir, mustWork = TRUE))
message("Run status: ", normalizePath(run_status_file, mustWork = TRUE))
message(
  "Gene availability: ",
  normalizePath(gene_availability_file, mustWork = TRUE)
)

if (n_failed > 0) {
  message(
    "Failed genes: ",
    paste(
      final_status_table$target_gene[
        final_status_table$status == "failed"
      ],
      collapse = ", "
    )
  )
}

message("============================================================")
