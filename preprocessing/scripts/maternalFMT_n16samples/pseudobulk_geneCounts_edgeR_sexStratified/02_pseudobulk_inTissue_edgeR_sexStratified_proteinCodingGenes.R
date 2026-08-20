# ==============================================================================
# 02_pseudobulk_inTissue_edgeR_sexStratified_proteinCodingGenes_Ensembl115.R
#
# PURPOSE
#   1. Read raw Space Ranger gene-count matrices for maternal FMT Visium samples.
#   2. Retain only spatial barcodes marked as in_tissue == 1.
#   3. Sum counts across all in-tissue spots within each sample, creating one
#      pseudobulk raw-count profile per sample.
#   4. Exclude the four samples removed after QC:
#        20_1F, 12_3F, 15_1M, 20_3M
#   5. Retrieve mouse gene annotation from the pinned Ensembl 115 BioMart
#      archive, selected as the newest completed fixed archive.
#   6. Retain only genes annotated by BioMart as gene_biotype == protein_coding.
#   7. Run two independent edgeR analyses:
#        - Male:   ASD vs Neurotypical
#        - Female: ASD vs Neurotypical
#   8. Apply edgeR::filterByExpr separately within each sex-specific analysis,
#      after the protein-coding filter.
#   9. Keep all outputs as R objects only. Nothing is written to disk.
#
# INTERPRETATION OF logFC
#   Positive logFC: up in ASD relative to Neurotypical.
#   Negative logFC: down in ASD relative to Neurotypical.
#
# MAIN OBJECTS CREATED AT THE END
#   gene_counts_per_sample_raw_in_tissue_allGenes
#       Gene x sample pseudobulk raw-count matrix before BioMart filtering.
#
#   gene_counts_per_sample_raw_in_tissue_proteinCoding
#       Protein-coding gene x sample matrix used for edgeR.
#
#   gene_counts_per_sample_raw_in_tissue
#       Alias of the protein-coding matrix used for edgeR.
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
#   biomart_version_info
#       R, Bioconductor, biomaRt package, Ensembl release, assembly and host.
#
#   biomart_mouse_gene_annotation
#       Mouse Ensembl 115 annotation downloaded from BioMart.
#
#   protein_coding_filter_status
#       Protein-coding-filter status for every gene in the raw matrix.
#
#   filterByExpr_status_combined
#       Male/female filterByExpr status for protein-coding genes.
#
#   all_gene_filter_status_combined
#       Complete status for every original gene: BioMart protein-coding filter
#       followed by the sex-specific filterByExpr results.
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
  "biomaRt",
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
    paste(missing_packages, collapse = ", "),
    "\n\nInstall missing Bioconductor packages with:\n",
    "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) ",
    "install.packages(\"BiocManager\")\n",
    "BiocManager::install(c(\"edgeR\", ",
    "\"limma\", \"biomaRt\"), ask = FALSE, update = FALSE)"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tibble)
})

message("R version: ", R.version.string)
message("Seurat version: ", as.character(packageVersion("Seurat")))
message("edgeR version: ", as.character(packageVersion("edgeR")))
message("limma version: ", as.character(packageVersion("limma")))
message("biomaRt package version: ", as.character(packageVersion("biomaRt")))


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

# Ensembl release is pinned for reproducibility.
# Ensembl 115 (September 2025) is the newest completed archive available
# at the time this script was prepared. It is intentionally not tied to
# the annotation release used in another project.
ensembl_release <- 115L
ensembl_biomart <- "genes"
ensembl_dataset <- "mmusculus_gene_ensembl"
reference_annotation_label <- paste0(
  "Protein-coding annotation from fixed Ensembl ",
  ensembl_release,
  " BioMart archive"
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
    use.names = FALSE,
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
    stop("Gene identifiers are missing for sample: ", sample_id)
  }

  # Remove an optional Ensembl version suffix, for example ENSMUSG... .12.
  normalized_ensembl_gene_ids <- sub(
    "\\.[0-9]+$",
    "",
    rownames(counts_matrix)
  )

  if (anyDuplicated(normalized_ensembl_gene_ids)) {
    duplicated_ids <- unique(
      normalized_ensembl_gene_ids[duplicated(normalized_ensembl_gene_ids)]
    )

    stop(
      "Duplicated Ensembl gene IDs after removing version suffixes for sample ",
      sample_id,
      ": ",
      paste(head(duplicated_ids, 20), collapse = ", ")
    )
  }

  rownames(counts_matrix) <- normalized_ensembl_gene_ids

  if (is.null(colnames(counts_matrix))) {
    stop("Barcode names are missing for sample: ", sample_id)
  }

  if (anyDuplicated(rownames(counts_matrix))) {
    stop(
      "Duplicated gene identifiers remain after Read10X_h5 for sample: ",
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
gene_counts_per_sample_raw_in_tissue_allGenes <- do.call(
  cbind,
  sample_gene_counts
)

rownames(gene_counts_per_sample_raw_in_tissue_allGenes) <- reference_genes
colnames(gene_counts_per_sample_raw_in_tissue_allGenes) <- retained_samples
storage.mode(gene_counts_per_sample_raw_in_tissue_allGenes) <- "numeric"

if (anyNA(gene_counts_per_sample_raw_in_tissue_allGenes)) {
  stop("Pseudobulk in-tissue count matrix contains missing values.")
}

if (any(!is.finite(gene_counts_per_sample_raw_in_tissue_allGenes))) {
  stop("Pseudobulk in-tissue count matrix contains non-finite values.")
}

if (any(gene_counts_per_sample_raw_in_tissue_allGenes < 0)) {
  stop("Pseudobulk in-tissue count matrix contains negative values.")
}

if (!all(
  gene_counts_per_sample_raw_in_tissue_allGenes ==
    floor(gene_counts_per_sample_raw_in_tissue_allGenes)
)) {
  stop("Pseudobulk count matrix contains non-integer count values.")
}

if (!identical(
  colnames(gene_counts_per_sample_raw_in_tissue_allGenes),
  as.character(sample_metadata_analysis$sample_ID)
)) {
  stop("Count-matrix columns do not match the retained metadata order.")
}

message(
  "\nPseudobulk in-tissue all-gene matrix dimensions: ",
  nrow(gene_counts_per_sample_raw_in_tissue_allGenes),
  " genes x ",
  ncol(gene_counts_per_sample_raw_in_tissue_allGenes),
  " samples"
)


# ==============================================================================
# 16. Connect to the pinned Ensembl 115 BioMart archive
#
# IMPORTANT
# - packageVersion("biomaRt") reports the installed R package version.
# - ensembl_release reports the annotation database release.
# These are different version numbers and both are recorded below.
# ==============================================================================
connect_mouse_ensembl_biomart <- function(
    ensembl_release,
    ensembl_biomart,
    ensembl_dataset
) {

  message(
    "Connecting to Ensembl BioMart release ",
    ensembl_release,
    " using dataset: ",
    ensembl_dataset
  )

  tryCatch(
    biomaRt::useEnsembl(
      biomart = ensembl_biomart,
      dataset = ensembl_dataset,
      version = ensembl_release
    ),
    error = function(error_condition) {
      stop(
        "Could not connect to Ensembl release ",
        ensembl_release,
        " BioMart archive. Original error:\n",
        conditionMessage(error_condition)
      )
    }
  )
}

mouse_ensembl_mart <- connect_mouse_ensembl_biomart(
  ensembl_release = ensembl_release,
  ensembl_biomart = ensembl_biomart,
  ensembl_dataset = ensembl_dataset
)

# Verify the requested Ensembl release and retrieve dataset/assembly metadata.
ensembl_mart_version_table <- biomaRt::listEnsembl(
  version = ensembl_release
) |>
  tibble::as_tibble()

ensembl_genes_version_row <- ensembl_mart_version_table |>
  dplyr::filter(biomart == ensembl_biomart)

if (nrow(ensembl_genes_version_row) != 1) {
  stop(
    "Could not uniquely verify Ensembl Genes release ",
    ensembl_release,
    "."
  )
}

if (!grepl(
  paste0("(^|[^0-9])", ensembl_release, "([^0-9]|$)"),
  ensembl_genes_version_row$version[[1]]
)) {
  stop(
    "BioMart release verification failed. Reported mart version: ",
    ensembl_genes_version_row$version[[1]]
  )
}

biomart_dataset_table <- biomaRt::listDatasets(
  mouse_ensembl_mart
) |>
  tibble::as_tibble()

biomart_mouse_dataset_info <- biomart_dataset_table |>
  dplyr::filter(dataset == ensembl_dataset)

if (nrow(biomart_mouse_dataset_info) != 1) {
  stop(
    "Could not uniquely identify dataset metadata for: ",
    ensembl_dataset
  )
}

bioconductor_version <- if (
  requireNamespace("BiocManager", quietly = TRUE)
) {
  as.character(BiocManager::version())
} else {
  NA_character_
}

biomart_version_info <- tibble::tibble(
  r_version = R.version.string,
  bioconductor_version = bioconductor_version,
  biomaRt_package_version = as.character(packageVersion("biomaRt")),
  ensembl_release = ensembl_release,
  ensembl_release_date = "September 2025",
  ensembl_selection_policy = "newest completed fixed archive at script creation",
  ensembl_mart_version = ensembl_genes_version_row$version[[1]],
  biomart_name = ensembl_biomart,
  dataset = ensembl_dataset,
  dataset_description = biomart_mouse_dataset_info$description[[1]],
  genome_assembly = biomart_mouse_dataset_info$version[[1]],
  reference_annotation = reference_annotation_label,
  biomart_host = methods::slot(mouse_ensembl_mart, "host"),
  query_date = as.character(Sys.Date())
)

message("\nBioMart reproducibility information:")
print(biomart_version_info, width = Inf)

# ==============================================================================
# 17. Download mouse gene annotation from BioMart
# ==============================================================================
required_biomart_attributes <- c(
  "ensembl_gene_id",
  "external_gene_name",
  "gene_biotype"
)

available_biomart_attributes <- biomaRt::listAttributes(
  mouse_ensembl_mart
)$name

missing_biomart_attributes <- setdiff(
  required_biomart_attributes,
  available_biomart_attributes
)

if (length(missing_biomart_attributes) > 0) {
  stop(
    "Required BioMart attributes are unavailable: ",
    paste(missing_biomart_attributes, collapse = ", ")
  )
}

biomart_mouse_gene_annotation <- biomaRt::getBM(
  attributes = required_biomart_attributes,
  mart = mouse_ensembl_mart,
  uniqueRows = TRUE
) |>
  tibble::as_tibble() |>
  transmute(
    ensembl_gene_id = sub(
      "\\.[0-9]+$",
      "",
      trimws(as.character(ensembl_gene_id))
    ),
    external_gene_name = trimws(as.character(external_gene_name)),
    gene_biotype = trimws(as.character(gene_biotype))
  ) |>
  filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != ""
  ) |>
  mutate(
    external_gene_name = if_else(
      is.na(external_gene_name),
      "",
      external_gene_name
    ),
    gene_biotype = if_else(
      is.na(gene_biotype),
      "",
      gene_biotype
    ),
    has_gene_symbol = external_gene_name != ""
  ) |>
  arrange(
    ensembl_gene_id,
    desc(has_gene_symbol),
    external_gene_name
  ) |>
  distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  ) |>
  mutate(
    gene = if_else(
      external_gene_name == "",
      ensembl_gene_id,
      external_gene_name
    )
  ) |>
  select(
    ensembl_gene_id,
    gene,
    external_gene_name,
    gene_biotype
  )

if (nrow(biomart_mouse_gene_annotation) == 0) {
  stop("BioMart returned no mouse gene annotation.")
}

if (anyDuplicated(biomart_mouse_gene_annotation$ensembl_gene_id)) {
  stop("Duplicated Ensembl gene IDs remain in BioMart annotation.")
}

message(
  "BioMart mouse genes retrieved: ",
  format(nrow(biomart_mouse_gene_annotation), big.mark = ",")
)


# ==============================================================================
# 18. Filter the pseudobulk matrix to protein-coding genes
#
# This filter is applied before filterByExpr and before edgeR modelling.
# ==============================================================================
protein_coding_filter_status <- tibble::tibble(
  ensembl_gene_id = rownames(
    gene_counts_per_sample_raw_in_tissue_allGenes
  )
) |>
  left_join(
    biomart_mouse_gene_annotation,
    by = "ensembl_gene_id"
  ) |>
  mutate(
    found_in_biomart = !is.na(gene_biotype),
    passed_protein_coding_filter =
      found_in_biomart & gene_biotype == "protein_coding",
    protein_coding_filter_status = case_when(
      !found_in_biomart ~ "filtered_out_not_found_in_biomart",
      gene_biotype == "protein_coding" ~ "passed_protein_coding_filter",
      TRUE ~ "filtered_out_non_protein_coding"
    ),
    gene = if_else(
      is.na(gene) | gene == "",
      ensembl_gene_id,
      gene
    )
  )

protein_coding_gene_ids <- protein_coding_filter_status |>
  filter(passed_protein_coding_filter) |>
  pull(ensembl_gene_id)

if (length(protein_coding_gene_ids) == 0) {
  stop("No protein-coding genes matched between BioMart and the count matrix.")
}

gene_counts_per_sample_raw_in_tissue_proteinCoding <-
  gene_counts_per_sample_raw_in_tissue_allGenes[
    protein_coding_gene_ids,
    ,
    drop = FALSE
  ]

# Convenient alias: this is the matrix used in the edgeR analyses below.
gene_counts_per_sample_raw_in_tissue <-
  gene_counts_per_sample_raw_in_tissue_proteinCoding

protein_coding_gene_annotation <- protein_coding_filter_status |>
  filter(passed_protein_coding_filter) |>
  select(
    ensembl_gene_id,
    gene,
    external_gene_name,
    gene_biotype
  )

if (!identical(
  rownames(gene_counts_per_sample_raw_in_tissue_proteinCoding),
  protein_coding_gene_annotation$ensembl_gene_id
)) {
  stop("Protein-coding count matrix and annotation order do not match.")
}

message(
  "Genes before BioMart protein-coding filter: ",
  format(nrow(gene_counts_per_sample_raw_in_tissue_allGenes), big.mark = ",")
)
message(
  "Protein-coding genes retained: ",
  format(nrow(gene_counts_per_sample_raw_in_tissue_proteinCoding), big.mark = ",")
)
message(
  "Genes removed by protein-coding filter: ",
  format(
    nrow(gene_counts_per_sample_raw_in_tissue_allGenes) -
      nrow(gene_counts_per_sample_raw_in_tissue_proteinCoding),
    big.mark = ","
  )
)


# ==============================================================================
# 19. Function for one sex-specific edgeR analysis
# ==============================================================================
run_edgeR_sex_specific <- function(
    sex_label,
    raw_count_matrix,
    analysis_metadata,
    filter_parameters,
    gene_annotation
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
    ensembl_gene_id = rownames(sex_counts),
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
  ) |>
    left_join(
      gene_annotation,
      by = "ensembl_gene_id"
    ) |>
    mutate(
      gene = if_else(
        is.na(gene) | gene == "",
        ensembl_gene_id,
        gene
      )
    ) |>
    select(
      ensembl_gene_id,
      gene,
      gene_biotype,
      everything(),
      -external_gene_name
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
    tibble::rownames_to_column(var = "ensembl_gene_id") |>
    left_join(
      gene_annotation,
      by = "ensembl_gene_id"
    ) |>
    mutate(
      gene = if_else(
        is.na(gene) | gene == "",
        ensembl_gene_id,
        gene
      ),
      sex = sex_label,
      comparison = "ASD_vs_Neurotypical",
      control_group = "Neurotypical",
      experimental_group = "ASD",
      significant_FDR_0.05 = FDR <= 0.05,
      regulation = if_else(logFC >= 0, "up", "down")
    ) |>
    select(
      ensembl_gene_id,
      gene,
      gene_biotype,
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
# 20. Run the two independent edgeR analyses
# ==============================================================================
edgeR_male_analysis <- run_edgeR_sex_specific(
  sex_label = "Male",
  raw_count_matrix = gene_counts_per_sample_raw_in_tissue,
  analysis_metadata = sample_metadata_analysis,
  filter_parameters = filterByExpr_parameters,
  gene_annotation = protein_coding_gene_annotation
)

edgeR_female_analysis <- run_edgeR_sex_specific(
  sex_label = "Female",
  raw_count_matrix = gene_counts_per_sample_raw_in_tissue,
  analysis_metadata = sample_metadata_analysis,
  filter_parameters = filterByExpr_parameters,
  gene_annotation = protein_coding_gene_annotation
)


# ==============================================================================
# 21. Main result list: two full data frames
# ==============================================================================
edgeR_results_by_sex <- list(
  male = edgeR_male_analysis$results,
  female = edgeR_female_analysis$results
)

# Optional convenient aliases.
male_results_df <- edgeR_results_by_sex$male
female_results_df <- edgeR_results_by_sex$female


# ==============================================================================
# 22. Filtering-status objects
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
    ensembl_gene_id,
    gene,
    gene_biotype,
    male_total_count = total_count,
    male_mean_CPM = mean_CPM,
    male_passed_filterByExpr = passed_filterByExpr
  ) |>
  full_join(
    filterByExpr_status_by_sex$female |>
      select(
        ensembl_gene_id,
        gene,
        gene_biotype,
        female_total_count = total_count,
        female_mean_CPM = mean_CPM,
        female_passed_filterByExpr = passed_filterByExpr
      ),
    by = c("ensembl_gene_id", "gene", "gene_biotype")
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

all_gene_filter_status_combined <- protein_coding_filter_status |>
  left_join(
    filterByExpr_status_combined |>
      select(
        ensembl_gene_id,
        male_passed_filterByExpr,
        female_passed_filterByExpr,
        filterByExpr_category = filter_category
      ),
    by = "ensembl_gene_id"
  ) |>
  mutate(
    final_filter_status = case_when(
      !passed_protein_coding_filter ~ protein_coding_filter_status,
      filterByExpr_category == "passed_in_both" ~
        "protein_coding_passed_filterByExpr_in_both",
      filterByExpr_category == "passed_male_only" ~
        "protein_coding_passed_filterByExpr_male_only",
      filterByExpr_category == "passed_female_only" ~
        "protein_coding_passed_filterByExpr_female_only",
      filterByExpr_category == "filtered_out_in_both" ~
        "protein_coding_filtered_out_by_filterByExpr_in_both",
      TRUE ~ "unexpected_filter_status"
    )
  )


# ==============================================================================
# 23. Keep complete edgeR analysis objects in one list
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
# 24. Basic final validation
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

if (
  nrow(filterByExpr_status_by_sex$male) !=
    nrow(gene_counts_per_sample_raw_in_tissue_proteinCoding)
) {
  stop("Male filtering-status table does not contain all protein-coding genes.")
}

if (
  nrow(filterByExpr_status_by_sex$female) !=
    nrow(gene_counts_per_sample_raw_in_tissue_proteinCoding)
) {
  stop("Female filtering-status table does not contain all protein-coding genes.")
}

if (
  nrow(protein_coding_filter_status) !=
    nrow(gene_counts_per_sample_raw_in_tissue_allGenes)
) {
  stop("Protein-coding status table does not contain every original gene.")
}

if (any(edgeR_results_by_sex$male$gene_biotype != "protein_coding")) {
  stop("Non-protein-coding genes are present in male edgeR results.")
}

if (any(edgeR_results_by_sex$female$gene_biotype != "protein_coding")) {
  stop("Non-protein-coding genes are present in female edgeR results.")
}


# ==============================================================================
# 25. Display main outputs in the R console
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

message("\nBioMart protein-coding filter status:")
print(
  protein_coding_filter_status |>
    count(protein_coding_filter_status, name = "n_genes") |>
    arrange(desc(n_genes)),
  n = Inf
)

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
message("  gene_counts_per_sample_raw_in_tissue_allGenes")
message("  gene_counts_per_sample_raw_in_tissue_proteinCoding")
message("  gene_counts_per_sample_raw_in_tissue")
message("  biomart_version_info")
message("  biomart_mouse_gene_annotation")
message("  protein_coding_gene_annotation")
message("  protein_coding_filter_status")
message("  all_gene_filter_status_combined")
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
# View(biomart_version_info)
# View(protein_coding_filter_status)
# View(all_gene_filter_status_combined)

# Examples of later filtering, without changing the stored full results:
# male_FDR_005 <- edgeR_results_by_sex$male |>
#   filter(FDR <= 0.05)
#
# female_FDR_005 <- edgeR_results_by_sex$female |>
#   filter(FDR <= 0.05)

# ==============================================================================
# End
# ==============================================================================
edgeR_results_by_sex$male %>% 
  filter(PValue < 0.01) %>% 
  filter(abs(logFC) > 0.5) %>% dim 


edgeR_results_by_sex$female %>% 
  # filter(PValue < 0.01) %>% 
  # filter(abs(logFC) > 0.5) %>% dim
  filter(FDR < 0.05) %>% 
  filter(abs(logFC) > 0.5) %>% dim 

