# ==============================================================================
# 01_pseudobulk_inTissue_edgeR_sexStratified_n16samples.R
#
# PURPOSE
#   1. Read raw Space Ranger gene-count matrices for maternal FMT Visium samples.
#   2. Retain only spatial barcodes marked as in_tissue == 1.
#   3. Sum counts across all in-tissue spots within each sample, creating one
#      pseudobulk raw-count profile per sample.
#   4. Exclude the four samples removed after QC:
#        20_1F, 12_3F, 15_1M, 20_3M
#   5. Run two independent edgeR analyses:
#        - Male:   ASD vs Neurotypical
#        - Female: ASD vs Neurotypical
#   6. Apply edgeR::filterByExpr separately within each sex-specific analysis.
#   7. Keep all outputs as R objects only. Nothing is written to disk.
#
# INTERPRETATION OF logFC
#   Positive logFC: up in ASD relative to Neurotypical.
#   Negative logFC: down in ASD relative to Neurotypical.
#
# MAIN OBJECTS CREATED AT THE END
#   gene_counts_per_sample_raw_in_tissue
#       Gene x sample pseudobulk raw-count matrix for the 16 retained samples.
#
#   sample_inclusion_status
#       Metadata for all 20 samples, including inclusion/exclusion status.
#
#   sample_metadata_analysis
#       Metadata for the 16 retained samples used in edgeR.
#
#   pseudobulk_input_summary
#       Per-sample summary of raw barcodes, in-tissue barcodes and summed UMIs.
#
#   edgeR_results_by_sex
#       List containing two complete result data frames:
#         edgeR_results_by_sex$male
#         edgeR_results_by_sex$female
#
#   filterByExpr_status_by_sex
#       List showing, for every gene, whether it passed filterByExpr:
#         filterByExpr_status_by_sex$male
#         filterByExpr_status_by_sex$female
#
#   filterByExpr_status_combined
#       Combined male/female filtering status for every gene.
#
#   edgeR_analysis_objects
#       Complete edgeR objects for inspection: design, DGEList, fit and QL test.
#
#   edgeR_analysis_summary
#       Summary of sample and gene numbers for both analyses.
# ==============================================================================


# ==============================================================================
# 1. Check and load packages
# ==============================================================================
required_packages <- c(
  "Seurat",
  "Matrix",
  "edgeR",
  "limma",
  "dplyr",
  "tibble",
  "hdf5r"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tibble)
})

message("Seurat version: ", as.character(packageVersion("Seurat")))
message("edgeR version: ", as.character(packageVersion("edgeR")))


# ==============================================================================
# 2. Define project and input paths
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


# ==============================================================================
# 3. Define samples excluded after QC
# ==============================================================================
excluded_samples <- c(
  "20_1F",
  "12_3F",
  "15_1M",
  "20_3M"
)

expected_number_of_all_samples <- 20L
expected_number_of_retained_samples <- 16L


# ==============================================================================
# 4. Define filterByExpr parameters explicitly
#
# These are the standard edgeR filterByExpr.default thresholds. They are written
# explicitly here so that the filtering settings are visible and reproducible.
# ==============================================================================
filterByExpr_parameters <- list(
  min.count = 10,
  min.total.count = 15,
  large.n = 10,
  min.prop = 0.7
)


# ==============================================================================
# 5. Check input paths
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

message("Project directory: ", project_dir)
message("Space Ranger directory: ", path_to_data)
message("Metadata file: ", metadata_file)


# ==============================================================================
# 6. Read and validate metadata
# ==============================================================================
metadata_autismFMT <- read.delim(
  metadata_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  tibble::as_tibble()

required_metadata_columns <- c(
  "sample_ID",
  "slide_ID",
  "slide_area",
  "experiment",
  "mouse_genotype",
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

metadata_autismFMT <- metadata_autismFMT |>
  mutate(
    across(
      all_of(required_metadata_columns),
      ~ trimws(as.character(.x))
    )
  )

if (anyNA(metadata_autismFMT$sample_ID)) {
  stop("Metadata contains missing sample_ID values.")
}

if (any(metadata_autismFMT$sample_ID == "")) {
  stop("Metadata contains empty sample_ID values.")
}

if (anyDuplicated(metadata_autismFMT$sample_ID)) {
  duplicated_samples <- unique(
    metadata_autismFMT$sample_ID[duplicated(metadata_autismFMT$sample_ID)]
  )

  stop(
    "Duplicated sample_ID values: ",
    paste(duplicated_samples, collapse = ", ")
  )
}

if (nrow(metadata_autismFMT) != expected_number_of_all_samples) {
  stop(
    "Expected ",
    expected_number_of_all_samples,
    " samples in metadata, but found ",
    nrow(metadata_autismFMT),
    "."
  )
}

missing_excluded_samples <- setdiff(
  excluded_samples,
  metadata_autismFMT$sample_ID
)

if (length(missing_excluded_samples) > 0) {
  stop(
    "Excluded samples missing from metadata: ",
    paste(missing_excluded_samples, collapse = ", ")
  )
}

unexpected_group_values <- setdiff(
  unique(metadata_autismFMT$fmt_donor_group),
  c("ASD", "Neurotypical")
)

if (length(unexpected_group_values) > 0) {
  stop(
    "Unexpected fmt_donor_group values: ",
    paste(unexpected_group_values, collapse = ", "),
    ". Expected only ASD and Neurotypical."
  )
}

unexpected_sex_values <- setdiff(
  unique(metadata_autismFMT$sex),
  c("Male", "Female")
)

if (length(unexpected_sex_values) > 0) {
  stop(
    "Unexpected sex values: ",
    paste(unexpected_sex_values, collapse = ", "),
    ". Expected only Male and Female."
  )
}


# ==============================================================================
# 7. Create sample inclusion/exclusion status
# ==============================================================================
sample_inclusion_status <- metadata_autismFMT |>
  mutate(
    included_in_analysis = !sample_ID %in% excluded_samples,
    analysis_status = if_else(
      included_in_analysis,
      "included",
      "excluded_after_QC"
    ),
    exclusion_reason = if_else(
      included_in_analysis,
      NA_character_,
      "excluded_spatial_QC"
    ),
    biological_group = if_else(
      fmt_donor_group == "ASD",
      "experimental_ASD",
      "control_Neurotypical"
    )
  )

sample_metadata_analysis <- sample_inclusion_status |>
  filter(included_in_analysis) |>
  mutate(
    fmt_donor_group = factor(
      fmt_donor_group,
      levels = c("Neurotypical", "ASD")
    ),
    sex = factor(
      sex,
      levels = c("Male", "Female")
    )
  )

if (nrow(sample_metadata_analysis) != expected_number_of_retained_samples) {
  stop(
    "Expected ",
    expected_number_of_retained_samples,
    " retained samples, but found ",
    nrow(sample_metadata_analysis),
    "."
  )
}

retained_samples <- as.character(sample_metadata_analysis$sample_ID)

message("\nSamples retained for analysis: ", length(retained_samples))
message(paste(retained_samples, collapse = ", "))
message("\nSamples excluded after QC: ", paste(excluded_samples, collapse = ", "))

message("\nSample numbers by sex and group:")
print(
  sample_metadata_analysis |>
    count(sex, fmt_donor_group, name = "n_samples"),
  n = Inf
)


# ==============================================================================
# 8. Find the Space Ranger directory for one sample
# ==============================================================================
find_sample_directory <- function(root_directory, sample_id) {

  exact_directory <- file.path(root_directory, sample_id)

  if (dir.exists(exact_directory)) {
    return(normalizePath(exact_directory, mustWork = TRUE))
  }

  possible_directories <- list.dirs(
    root_directory,
    recursive = TRUE,
    full.names = TRUE
  )

  matching_directories <- possible_directories[
    basename(possible_directories) == sample_id
  ]

  if (length(matching_directories) == 0) {
    stop("Could not find a directory for sample: ", sample_id)
  }

  if (length(matching_directories) > 1) {
    stop(
      "More than one directory was found for sample ",
      sample_id,
      ":\n",
      paste(matching_directories, collapse = "\n")
    )
  }

  normalizePath(matching_directories, mustWork = TRUE)
}


# ==============================================================================
# 9. Find raw_feature_bc_matrix.h5 for one sample
# ==============================================================================
find_raw_h5_file <- function(sample_directory, sample_id) {

  standard_candidates <- c(
    file.path(sample_directory, "outs", "raw_feature_bc_matrix.h5"),
    file.path(sample_directory, "raw_feature_bc_matrix.h5")
  )

  existing_standard_candidates <- standard_candidates[
    file.exists(standard_candidates)
  ]

  if (length(existing_standard_candidates) == 1) {
    return(normalizePath(existing_standard_candidates, mustWork = TRUE))
  }

  if (length(existing_standard_candidates) > 1) {
    stop(
      "More than one standard raw H5 file was found for sample ",
      sample_id,
      ":\n",
      paste(existing_standard_candidates, collapse = "\n")
    )
  }

  recursive_candidates <- list.files(
    sample_directory,
    pattern = "^raw_feature_bc_matrix\\.h5$",
    recursive = TRUE,
    full.names = TRUE
  )

  recursive_candidates <- recursive_candidates[
    file.exists(recursive_candidates)
  ]

  if (length(recursive_candidates) == 0) {
    stop(
      "raw_feature_bc_matrix.h5 was not found for sample ",
      sample_id,
      " in directory:\n",
      sample_directory
    )
  }

  if (length(recursive_candidates) > 1) {
    stop(
      "More than one raw_feature_bc_matrix.h5 file was found for sample ",
      sample_id,
      ":\n",
      paste(recursive_candidates, collapse = "\n")
    )
  }

  normalizePath(recursive_candidates, mustWork = TRUE)
}


# ==============================================================================
# 10. Find tissue-position file for one sample
#
# Space Ranger versions may use either:
#   tissue_positions.csv       - newer format, usually with a header
#   tissue_positions_list.csv  - older format, usually without a header
# ==============================================================================
find_tissue_positions_file <- function(sample_directory, sample_id) {

  preferred_candidates <- c(
    file.path(sample_directory, "outs", "spatial", "tissue_positions.csv"),
    file.path(sample_directory, "outs", "spatial", "tissue_positions_list.csv"),
    file.path(sample_directory, "spatial", "tissue_positions.csv"),
    file.path(sample_directory, "spatial", "tissue_positions_list.csv")
  )

  existing_preferred_candidates <- preferred_candidates[
    file.exists(preferred_candidates)
  ]

  if (length(existing_preferred_candidates) >= 1) {
    # Prefer tissue_positions.csv over tissue_positions_list.csv when both exist.
    preferred_order <- order(
      basename(existing_preferred_candidates) != "tissue_positions.csv"
    )

    return(
      normalizePath(
        existing_preferred_candidates[preferred_order][1],
        mustWork = TRUE
      )
    )
  }

  recursive_candidates <- list.files(
    sample_directory,
    pattern = "^tissue_positions(_list)?\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  recursive_candidates <- recursive_candidates[
    file.exists(recursive_candidates)
  ]

  if (length(recursive_candidates) == 0) {
    stop(
      "No tissue_positions.csv or tissue_positions_list.csv file was found for sample ",
      sample_id,
      "."
    )
  }

  preferred_new_format <- recursive_candidates[
    basename(recursive_candidates) == "tissue_positions.csv"
  ]

  if (length(preferred_new_format) == 1) {
    return(normalizePath(preferred_new_format, mustWork = TRUE))
  }

  if (length(recursive_candidates) > 1) {
    stop(
      "More than one tissue-position file was found for sample ",
      sample_id,
      ":\n",
      paste(recursive_candidates, collapse = "\n")
    )
  }

  normalizePath(recursive_candidates, mustWork = TRUE)
}


# ==============================================================================
# 11. Read tissue-position file robustly
# ==============================================================================
read_tissue_positions <- function(tissue_positions_file, sample_id) {

  first_line <- readLines(
    tissue_positions_file,
    n = 1,
    warn = FALSE
  )

  has_header <- grepl(
    "barcode",
    tolower(first_line),
    fixed = TRUE
  ) && grepl(
    "in_tissue",
    tolower(first_line),
    fixed = TRUE
  )

  if (has_header) {
    tissue_positions <- read.csv(
      tissue_positions_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    tissue_positions <- read.csv(
      tissue_positions_file,
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if (ncol(tissue_positions) < 6) {
      stop(
        "Headerless tissue-position file has fewer than six columns for sample: ",
        sample_id
      )
    }

    colnames(tissue_positions)[1:6] <- c(
      "barcode",
      "in_tissue",
      "array_row",
      "array_col",
      "pxl_row_in_fullres",
      "pxl_col_in_fullres"
    )
  }

  required_position_columns <- c("barcode", "in_tissue")

  missing_position_columns <- setdiff(
    required_position_columns,
    colnames(tissue_positions)
  )

  if (length(missing_position_columns) > 0) {
    stop(
      "Missing tissue-position columns for sample ",
      sample_id,
      ": ",
      paste(missing_position_columns, collapse = ", ")
    )
  }

  tissue_positions <- tissue_positions |>
    transmute(
      barcode = trimws(as.character(barcode)),
      in_tissue = suppressWarnings(as.integer(in_tissue))
    )

  if (anyNA(tissue_positions$barcode) || any(tissue_positions$barcode == "")) {
    stop("Missing or empty barcodes in tissue-position file for sample: ", sample_id)
  }

  if (anyNA(tissue_positions$in_tissue)) {
    stop("Non-numeric in_tissue values found for sample: ", sample_id)
  }

  unexpected_in_tissue_values <- setdiff(
    unique(tissue_positions$in_tissue),
    c(0L, 1L)
  )

  if (length(unexpected_in_tissue_values) > 0) {
    stop(
      "Unexpected in_tissue values for sample ",
      sample_id,
      ": ",
      paste(unexpected_in_tissue_values, collapse = ", ")
    )
  }

  if (anyDuplicated(tissue_positions$barcode)) {
    stop("Duplicated barcodes in tissue-position file for sample: ", sample_id)
  }

  tissue_positions
}


# ==============================================================================
# 12. Extract Gene Expression matrix from Read10X_h5 output
# ==============================================================================
extract_gene_expression_matrix <- function(read10x_object, sample_id) {

  if (inherits(read10x_object, "Matrix")) {
    return(read10x_object)
  }

  if (is.matrix(read10x_object)) {
    return(Matrix::Matrix(read10x_object, sparse = TRUE))
  }

  if (is.list(read10x_object)) {

    if ("Gene Expression" %in% names(read10x_object)) {
      return(read10x_object[["Gene Expression"]])
    }

    if (length(read10x_object) == 1) {
      return(read10x_object[[1]])
    }

    stop(
      "Multiple feature types were found for sample ",
      sample_id,
      ", but no 'Gene Expression' matrix was present. Available types: ",
      paste(names(read10x_object), collapse = ", ")
    )
  }

  stop("Unsupported object returned by Read10X_h5 for sample: ", sample_id)
}


# ==============================================================================
# 13. Read raw matrices and sum counts only across in-tissue spots
# ==============================================================================
sample_gene_counts <- vector(
  mode = "list",
  length = length(retained_samples)
)

names(sample_gene_counts) <- retained_samples

pseudobulk_input_summary_list <- vector(
  mode = "list",
  length = length(retained_samples)
)

names(pseudobulk_input_summary_list) <- retained_samples

for (sample_id in retained_samples) {

  message("\nProcessing sample: ", sample_id)

  sample_directory <- find_sample_directory(
    root_directory = path_to_data,
    sample_id = sample_id
  )

  raw_h5_file <- find_raw_h5_file(
    sample_directory = sample_directory,
    sample_id = sample_id
  )

  tissue_positions_file <- find_tissue_positions_file(
    sample_directory = sample_directory,
    sample_id = sample_id
  )

  message("  Raw matrix: ", raw_h5_file)
  message("  Tissue positions: ", tissue_positions_file)

  raw_10x_object <- Seurat::Read10X_h5(
    filename = raw_h5_file,
    use.names = TRUE,
    unique.features = TRUE
  )

  counts_matrix <- extract_gene_expression_matrix(
    read10x_object = raw_10x_object,
    sample_id = sample_id
  )

  if (!inherits(counts_matrix, "Matrix")) {
    counts_matrix <- Matrix::Matrix(counts_matrix, sparse = TRUE)
  }

  if (nrow(counts_matrix) == 0) {
    stop("No genes were found for sample: ", sample_id)
  }

  if (ncol(counts_matrix) == 0) {
    stop("No barcodes were found for sample: ", sample_id)
  }

  if (is.null(rownames(counts_matrix))) {
    stop("Gene names are missing for sample: ", sample_id)
  }

  if (is.null(colnames(counts_matrix))) {
    stop("Barcode names are missing for sample: ", sample_id)
  }

  if (anyDuplicated(rownames(counts_matrix))) {
    stop(
      "Duplicated gene names remain after Read10X_h5 for sample: ",
      sample_id
    )
  }

  tissue_positions <- read_tissue_positions(
    tissue_positions_file = tissue_positions_file,
    sample_id = sample_id
  )

  in_tissue_barcodes_from_positions <- tissue_positions |>
    filter(in_tissue == 1L) |>
    pull(barcode)

  if (length(in_tissue_barcodes_from_positions) == 0) {
    stop("No in-tissue barcodes were found for sample: ", sample_id)
  }

  matched_in_tissue_barcodes <- intersect(
    in_tissue_barcodes_from_positions,
    colnames(counts_matrix)
  )

  missing_in_tissue_barcodes <- setdiff(
    in_tissue_barcodes_from_positions,
    colnames(counts_matrix)
  )

  if (length(matched_in_tissue_barcodes) == 0) {
    stop(
      "None of the in-tissue barcodes matched the raw count matrix for sample: ",
      sample_id
    )
  }

  if (length(missing_in_tissue_barcodes) > 0) {
    warning(
      sample_id,
      ": ",
      length(missing_in_tissue_barcodes),
      " in-tissue barcodes were absent from raw_feature_bc_matrix.h5."
    )
  }

  in_tissue_counts_matrix <- counts_matrix[
    ,
    matched_in_tissue_barcodes,
    drop = FALSE
  ]

  gene_counts <- Matrix::rowSums(in_tissue_counts_matrix)

  gene_counts <- setNames(
    as.numeric(gene_counts),
    rownames(in_tissue_counts_matrix)
  )

  sample_gene_counts[[sample_id]] <- gene_counts

  sample_metadata_row <- sample_metadata_analysis |>
    filter(sample_ID == sample_id)

  pseudobulk_input_summary_list[[sample_id]] <- data.frame(
    sample_ID = sample_id,
    sex = as.character(sample_metadata_row$sex),
    fmt_donor_group = as.character(sample_metadata_row$fmt_donor_group),
    raw_matrix_file = raw_h5_file,
    tissue_positions_file = tissue_positions_file,
    genes = nrow(counts_matrix),
    raw_matrix_barcodes = ncol(counts_matrix),
    in_tissue_barcodes_in_positions = length(in_tissue_barcodes_from_positions),
    matched_in_tissue_barcodes = length(matched_in_tissue_barcodes),
    missing_in_tissue_barcodes = length(missing_in_tissue_barcodes),
    total_in_tissue_UMIs = sum(gene_counts),
    stringsAsFactors = FALSE
  )

  message("  genes: ", format(nrow(counts_matrix), big.mark = ","))
  message(
    "  raw-matrix barcodes: ",
    format(ncol(counts_matrix), big.mark = ",")
  )
  message(
    "  matched in-tissue barcodes: ",
    format(length(matched_in_tissue_barcodes), big.mark = ",")
  )
  message(
    "  total in-tissue UMIs: ",
    format(sum(gene_counts), big.mark = ",", scientific = FALSE)
  )

  rm(
    raw_10x_object,
    counts_matrix,
    tissue_positions,
    in_tissue_counts_matrix
  )

  invisible(gc())
}

pseudobulk_input_summary <- do.call(
  rbind,
  pseudobulk_input_summary_list
)

rownames(pseudobulk_input_summary) <- NULL


# ==============================================================================
# 14. Check and align genes between samples
# ==============================================================================
reference_genes <- names(sample_gene_counts[[1]])

identical_gene_sets <- vapply(
  sample_gene_counts,
  function(gene_counts) {
    setequal(names(gene_counts), reference_genes)
  },
  FUN.VALUE = logical(1)
)

if (!all(identical_gene_sets)) {
  problematic_samples <- names(identical_gene_sets)[!identical_gene_sets]

  stop(
    "Gene sets differ between samples: ",
    paste(problematic_samples, collapse = ", ")
  )
}

sample_gene_counts <- lapply(
  sample_gene_counts,
  function(gene_counts) {
    gene_counts[reference_genes]
  }
)


# ==============================================================================
# 15. Create gene x sample pseudobulk raw-count matrix
# ==============================================================================
gene_counts_per_sample_raw_in_tissue <- do.call(
  cbind,
  sample_gene_counts
)

rownames(gene_counts_per_sample_raw_in_tissue) <- reference_genes
colnames(gene_counts_per_sample_raw_in_tissue) <- retained_samples
storage.mode(gene_counts_per_sample_raw_in_tissue) <- "numeric"

if (anyNA(gene_counts_per_sample_raw_in_tissue)) {
  stop("Pseudobulk in-tissue count matrix contains missing values.")
}

if (any(!is.finite(gene_counts_per_sample_raw_in_tissue))) {
  stop("Pseudobulk in-tissue count matrix contains non-finite values.")
}

if (any(gene_counts_per_sample_raw_in_tissue < 0)) {
  stop("Pseudobulk in-tissue count matrix contains negative values.")
}

if (!all(gene_counts_per_sample_raw_in_tissue == floor(gene_counts_per_sample_raw_in_tissue))) {
  stop("Pseudobulk count matrix contains non-integer count values.")
}

if (!identical(
  colnames(gene_counts_per_sample_raw_in_tissue),
  as.character(sample_metadata_analysis$sample_ID)
)) {
  stop("Count-matrix columns do not match the retained metadata order.")
}

message(
  "\nPseudobulk in-tissue matrix dimensions: ",
  nrow(gene_counts_per_sample_raw_in_tissue),
  " genes x ",
  ncol(gene_counts_per_sample_raw_in_tissue),
  " samples"
)


# ==============================================================================
# 16. Function for one sex-specific edgeR analysis
# ==============================================================================
run_edgeR_sex_specific <- function(
    sex_label,
    raw_count_matrix,
    analysis_metadata,
    filter_parameters
) {

  message("\n", paste(rep("=", 80), collapse = ""))
  message("edgeR analysis: ", sex_label, " | ASD vs Neurotypical")
  message(paste(rep("=", 80), collapse = ""))

  sex_metadata <- analysis_metadata |>
    filter(as.character(sex) == sex_label)

  if (nrow(sex_metadata) == 0) {
    stop("No samples were found for sex: ", sex_label)
  }

  sex_sample_ids <- as.character(sex_metadata$sample_ID)

  missing_count_samples <- setdiff(
    sex_sample_ids,
    colnames(raw_count_matrix)
  )

  if (length(missing_count_samples) > 0) {
    stop(
      "Samples missing from count matrix for ",
      sex_label,
      ": ",
      paste(missing_count_samples, collapse = ", ")
    )
  }

  sex_counts <- raw_count_matrix[
    ,
    sex_sample_ids,
    drop = FALSE
  ]

  if (!identical(colnames(sex_counts), sex_sample_ids)) {
    stop("Sample order mismatch in the ", sex_label, " analysis.")
  }

  group <- factor(
    as.character(sex_metadata$fmt_donor_group),
    levels = c("Neurotypical", "ASD")
  )

  group_counts <- table(group)

  if (any(group_counts < 2)) {
    stop(
      sex_label,
      " analysis has fewer than two samples in at least one group: ",
      paste(names(group_counts), group_counts, sep = "=", collapse = ", ")
    )
  }

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  rownames(design) <- sex_sample_ids

  contrast_matrix <- limma::makeContrasts(
    ASD_vs_Neurotypical = ASD - Neurotypical,
    levels = design
  )

  dge_unfiltered <- edgeR::DGEList(
    counts = sex_counts,
    group = group,
    samples = data.frame(
      sample_ID = sex_sample_ids,
      sex = sex_label,
      fmt_donor_group = as.character(group),
      row.names = sex_sample_ids,
      stringsAsFactors = FALSE
    )
  )

  keep_genes <- edgeR::filterByExpr(
    dge_unfiltered,
    design = design,
    min.count = filter_parameters$min.count,
    min.total.count = filter_parameters$min.total.count,
    large.n = filter_parameters$large.n,
    min.prop = filter_parameters$min.prop
  )

  if (!is.logical(keep_genes) || length(keep_genes) != nrow(sex_counts)) {
    stop("Unexpected filterByExpr output for sex: ", sex_label)
  }

  cpm_unfiltered <- edgeR::cpm(
    dge_unfiltered,
    log = FALSE
  )

  filter_status <- tibble::tibble(
    gene = rownames(sex_counts),
    total_count = rowSums(sex_counts),
    mean_count = rowMeans(sex_counts),
    samples_with_nonzero_count = rowSums(sex_counts > 0),
    samples_with_count_ge_10 = rowSums(sex_counts >= 10),
    mean_CPM = rowMeans(cpm_unfiltered),
    max_CPM = apply(cpm_unfiltered, 1, max),
    passed_filterByExpr = keep_genes,
    filter_status = ifelse(
      keep_genes,
      "passed_filterByExpr",
      "filtered_out"
    )
  )

  dge_filtered <- dge_unfiltered[
    keep_genes,
    ,
    keep.lib.sizes = FALSE
  ]

  if (nrow(dge_filtered) == 0) {
    stop("No genes passed filterByExpr for sex: ", sex_label)
  }

  dge_filtered <- edgeR::calcNormFactors(
    dge_filtered,
    method = "TMM"
  )

  dge_filtered <- edgeR::estimateDisp(
    dge_filtered,
    design = design,
    robust = TRUE
  )

  ql_fit <- edgeR::glmQLFit(
    dge_filtered,
    design = design,
    robust = TRUE
  )

  ql_test <- edgeR::glmQLFTest(
    ql_fit,
    contrast = contrast_matrix[, "ASD_vs_Neurotypical"]
  )

  full_results <- edgeR::topTags(
    ql_test,
    n = Inf,
    sort.by = "PValue"
  )$table |>
    as.data.frame(stringsAsFactors = FALSE) |>
    tibble::rownames_to_column(var = "gene") |>
    mutate(
      sex = sex_label,
      comparison = "ASD_vs_Neurotypical",
      control_group = "Neurotypical",
      experimental_group = "ASD",
      significant_FDR_0.05 = FDR <= 0.05,
      regulation = case_when(
        logFC > 0 ~ "up",
        logFC < 0 ~ "down",
        TRUE ~ "unchanged"
      )
    ) |>
    select(
      gene,
      sex,
      comparison,
      control_group,
      experimental_group,
      logFC,
      logCPM,
      F,
      PValue,
      FDR,
      significant_FDR_0.05,
      regulation
    )

  analysis_summary <- tibble::tibble(
    sex = sex_label,
    comparison = "ASD_vs_Neurotypical",
    n_Neurotypical = unname(group_counts["Neurotypical"]),
    n_ASD = unname(group_counts["ASD"]),
    genes_before_filterByExpr = nrow(sex_counts),
    genes_passed_filterByExpr = sum(keep_genes),
    genes_filtered_out = sum(!keep_genes),
    significant_genes_FDR_0.05 = sum(full_results$FDR <= 0.05)
  )

  message("Samples: ", ncol(sex_counts))
  message(
    "Groups: Neurotypical=",
    unname(group_counts["Neurotypical"]),
    ", ASD=",
    unname(group_counts["ASD"])
  )
  message("Genes before filterByExpr: ", nrow(sex_counts))
  message("Genes passed filterByExpr: ", sum(keep_genes))
  message("Genes filtered out: ", sum(!keep_genes))
  message(
    "Genes significant at FDR <= 0.05: ",
    sum(full_results$FDR <= 0.05)
  )

  list(
    results = full_results,
    filter_status = filter_status,
    analysis_summary = analysis_summary,
    metadata = sex_metadata,
    counts_unfiltered = sex_counts,
    design = design,
    contrast_matrix = contrast_matrix,
    keep_genes = keep_genes,
    dge_unfiltered = dge_unfiltered,
    dge_filtered = dge_filtered,
    ql_fit = ql_fit,
    ql_test = ql_test
  )
}


# ==============================================================================
# 17. Run the two independent edgeR analyses
# ==============================================================================
edgeR_male_analysis <- run_edgeR_sex_specific(
  sex_label = "Male",
  raw_count_matrix = gene_counts_per_sample_raw_in_tissue,
  analysis_metadata = sample_metadata_analysis,
  filter_parameters = filterByExpr_parameters
)

edgeR_female_analysis <- run_edgeR_sex_specific(
  sex_label = "Female",
  raw_count_matrix = gene_counts_per_sample_raw_in_tissue,
  analysis_metadata = sample_metadata_analysis,
  filter_parameters = filterByExpr_parameters
)


# ==============================================================================
# 18. Main result list: two full data frames
# ==============================================================================
edgeR_results_by_sex <- list(
  male = edgeR_male_analysis$results,
  female = edgeR_female_analysis$results
)

# Optional convenient aliases.
male_results_df <- edgeR_results_by_sex$male
female_results_df <- edgeR_results_by_sex$female


# ==============================================================================
# 19. Filtering-status objects
#
# Each table includes all genes from the original raw-count matrix, including
# genes that did not pass filterByExpr and therefore were not tested by edgeR.
# ==============================================================================
filterByExpr_status_by_sex <- list(
  male = edgeR_male_analysis$filter_status,
  female = edgeR_female_analysis$filter_status
)

filterByExpr_status_combined <- filterByExpr_status_by_sex$male |>
  select(
    gene,
    male_total_count = total_count,
    male_mean_CPM = mean_CPM,
    male_passed_filterByExpr = passed_filterByExpr
  ) |>
  full_join(
    filterByExpr_status_by_sex$female |>
      select(
        gene,
        female_total_count = total_count,
        female_mean_CPM = mean_CPM,
        female_passed_filterByExpr = passed_filterByExpr
      ),
    by = "gene"
  ) |>
  mutate(
    filter_category = case_when(
      male_passed_filterByExpr & female_passed_filterByExpr ~
        "passed_in_both",
      male_passed_filterByExpr & !female_passed_filterByExpr ~
        "passed_male_only",
      !male_passed_filterByExpr & female_passed_filterByExpr ~
        "passed_female_only",
      TRUE ~
        "filtered_out_in_both"
    )
  )


# ==============================================================================
# 20. Keep complete edgeR analysis objects in one list
# ==============================================================================
edgeR_analysis_objects <- list(
  male = edgeR_male_analysis,
  female = edgeR_female_analysis
)

edgeR_analysis_summary <- bind_rows(
  edgeR_male_analysis$analysis_summary,
  edgeR_female_analysis$analysis_summary
)


# ==============================================================================
# 21. Basic final validation
# ==============================================================================
if (!all(c("male", "female") %in% names(edgeR_results_by_sex))) {
  stop("edgeR_results_by_sex does not contain both male and female results.")
}

if (nrow(edgeR_results_by_sex$male) != sum(edgeR_male_analysis$keep_genes)) {
  stop("Male result-row count does not match the number of tested genes.")
}

if (nrow(edgeR_results_by_sex$female) != sum(edgeR_female_analysis$keep_genes)) {
  stop("Female result-row count does not match the number of tested genes.")
}

if (nrow(filterByExpr_status_by_sex$male) != nrow(gene_counts_per_sample_raw_in_tissue)) {
  stop("Male filtering-status table does not contain all genes.")
}

if (nrow(filterByExpr_status_by_sex$female) != nrow(gene_counts_per_sample_raw_in_tissue)) {
  stop("Female filtering-status table does not contain all genes.")
}


# ==============================================================================
# 22. Display main outputs in the R console
# ==============================================================================
message("\n", paste(rep("=", 80), collapse = ""))
message("ANALYSIS COMPLETED SUCCESSFULLY")
message(paste(rep("=", 80), collapse = ""))

message("\nSample inclusion status:")
print(
  sample_inclusion_status |>
    select(
      sample_ID,
      fmt_donor_group,
      sex,
      included_in_analysis,
      analysis_status,
      exclusion_reason
    ),
  n = Inf
)

message("\nedgeR analysis summary:")
print(edgeR_analysis_summary, row.names = FALSE)

message("\nfilterByExpr categories across male and female analyses:")
print(
  filterByExpr_status_combined |>
    count(filter_category, name = "n_genes") |>
    arrange(desc(n_genes)),
  n = Inf
)

message("\nFirst rows of male ASD vs Neurotypical results:")
print(head(edgeR_results_by_sex$male, 10), row.names = FALSE)

message("\nFirst rows of female ASD vs Neurotypical results:")
print(head(edgeR_results_by_sex$female, 10), row.names = FALSE)

message("\nMain objects available in the R environment:")
message("  gene_counts_per_sample_raw_in_tissue")
message("  sample_inclusion_status")
message("  sample_metadata_analysis")
message("  pseudobulk_input_summary")
message("  filterByExpr_parameters")
message("  edgeR_results_by_sex$male")
message("  edgeR_results_by_sex$female")
message("  male_results_df")
message("  female_results_df")
message("  filterByExpr_status_by_sex$male")
message("  filterByExpr_status_by_sex$female")
message("  filterByExpr_status_combined")
message("  edgeR_analysis_objects")
message("  edgeR_analysis_summary")

# Examples for interactive inspection:
# View(edgeR_results_by_sex$male)
# View(edgeR_results_by_sex$female)
# View(filterByExpr_status_by_sex$male)
# View(filterByExpr_status_by_sex$female)
# View(filterByExpr_status_combined)

# Examples of later filtering, without changing the stored full results:
# male_FDR_005 <- edgeR_results_by_sex$male |>
#   filter(FDR <= 0.05)
#
# female_FDR_005 <- edgeR_results_by_sex$female |>
#   filter(FDR <= 0.05)

# ==============================================================================
# End
# ==============================================================================


# edgeR_results_by_sex$male %>% 
#   filter(PValue < 0.01) %>% 
#   filter(abs(logFC) > 0.5)

# edgeR_results_by_sex$female %>% 
#   filter(FDR < 0.05) %>% 
#   filter(abs(logFC) > 0.5) %>% 
#   .$gene %>% cat(sep = "\n")
