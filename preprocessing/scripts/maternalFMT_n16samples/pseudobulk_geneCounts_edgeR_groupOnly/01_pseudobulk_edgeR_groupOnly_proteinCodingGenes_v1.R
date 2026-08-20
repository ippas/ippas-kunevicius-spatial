#!/usr/bin/env Rscript

# ==============================================================================
# 02_pseudobulk_inTissue_edgeR_groupOnly_proteinCodingGenes_Ensembl115.R
#
# PURPOSE
#   1. Read raw Space Ranger gene-count matrices for maternal FMT Visium samples.
#   2. Retain only spatial barcodes marked as in_tissue == 1.
#   3. Sum counts across all in-tissue spots within each sample, creating one
#      pseudobulk raw-count profile per sample.
#   4. Exclude the four samples removed after QC:
#        20_1F, 12_3F, 15_1M, 20_3M
#   5. Read the locally saved Ensembl 115 Mus musculus protein-coding annotation.
#   6. Retain genes present in the local protein-coding annotation.
#   7. Apply edgeR::filterByExpr using the group-only design.
#   8. Run one edgeR quasi-likelihood analysis:
#        ASD vs Neurotypical
#      using all 16 retained samples and the model:
#        expression ~ fmt_donor_group
#   9. Save full edgeR results and complete filtering/QC summaries as TSV files.
#
# INTERPRETATION OF logFC
#   Positive logFC: higher expression in ASD than Neurotypical.
#   Negative logFC: lower expression in ASD than Neurotypical.
#
# MAIN RESULT COLUMNS
#   ensembl_gene_id
#   gene
#   logFC
#   logCPM
#   F
#   PValue
#   FDR
#   regulation
#
# OUTPUT DIRECTORY
#   /home/mateusz/projects/ippas-kunevicius-spatial/results/
#   maternalFMT_n16samples/pseudobulk_geneCounts_edgeR_groupOnly
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
    paste(missing_packages, collapse = ", "),
    "\n\nInstall missing Bioconductor packages with:\n",
    "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) ",
    "install.packages(\"BiocManager\")\n",
    "BiocManager::install(c(\"edgeR\", \"limma\"), ",
    "ask = FALSE, update = FALSE)"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tibble)
})

message("R version: ", R.version.string)
message("Seurat version: ", as.character(packageVersion("Seurat")))
message("edgeR version: ", as.character(packageVersion("edgeR")))
message("limma version: ", as.character(packageVersion("limma")))


# ==============================================================================
# 2. Define project, input and output paths
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

output_dir <- file.path(
  project_dir,
  "results",
  "maternalFMT_n16samples",
  "pseudobulk_geneCounts_edgeR_groupOnly"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Local Ensembl 115 annotation generated once from the official GTF file.
ensembl_release <- 115L
ensembl_genome_assembly <- "GRCm39"

ensembl115_directory <- file.path(
  project_dir,
  "data",
  "ensembl115"
)

ensembl115_protein_coding_rds <- file.path(
  ensembl115_directory,
  "ensembl115_mouse_proteinCodingGenes.rds"
)

ensembl115_annotation_metadata_file <- file.path(
  ensembl115_directory,
  "ensembl115_mouse_proteinCodingGenes_metadata.tsv"
)

reference_annotation_label <- paste0(
  "Mus musculus protein-coding genes from official Ensembl ",
  ensembl_release,
  " ",
  ensembl_genome_assembly,
  " GTF"
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
# ==============================================================================

filterByExpr_parameters <- list(
  min.count = 10,
  min.total.count = 15,
  large.n = 10,
  min.prop = 0.7
)


# ==============================================================================
# 5. Check input and output paths
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

if (!dir.exists(ensembl115_directory)) {
  stop(
    "Local Ensembl 115 directory does not exist: ",
    ensembl115_directory
  )
}

if (!file.exists(ensembl115_protein_coding_rds)) {
  stop(
    "Local Ensembl 115 protein-coding RDS file does not exist: ",
    ensembl115_protein_coding_rds,
    "\nRun download_Ensembl115_mouse_proteinCodingGenes_fromGTF.R first."
  )
}

if (!file.exists(ensembl115_annotation_metadata_file)) {
  stop(
    "Local Ensembl 115 annotation metadata file does not exist: ",
    ensembl115_annotation_metadata_file
  )
}

if (!dir.exists(output_dir)) {
  stop("Could not create output directory: ", output_dir)
}

message("Project directory: ", project_dir)
message("Space Ranger directory: ", path_to_data)
message("Metadata file: ", metadata_file)
message("Local Ensembl 115 annotation: ", ensembl115_protein_coding_rds)
message("Output directory: ", output_dir)




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
# 16. Read and validate the local Ensembl 115 protein-coding annotation
# ==============================================================================
ensembl115_annotation_metadata <- read.delim(
  ensembl115_annotation_metadata_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  tibble::as_tibble()

required_annotation_metadata_columns <- c(
  "species",
  "genome_assembly",
  "ensembl_release",
  "retained_feature_type",
  "retained_gene_biotype",
  "number_of_protein_coding_genes",
  "source_url"
)

missing_annotation_metadata_columns <- setdiff(
  required_annotation_metadata_columns,
  colnames(ensembl115_annotation_metadata)
)

if (length(missing_annotation_metadata_columns) > 0) {
  stop(
    "Missing columns in local Ensembl 115 metadata: ",
    paste(missing_annotation_metadata_columns, collapse = ", ")
  )
}

if (nrow(ensembl115_annotation_metadata) != 1L) {
  stop(
    "Expected exactly one row in local Ensembl 115 metadata, but found ",
    nrow(ensembl115_annotation_metadata),
    "."
  )
}

if (as.integer(ensembl115_annotation_metadata$ensembl_release[[1]]) !=
    ensembl_release) {
  stop(
    "Local annotation has Ensembl release ",
    ensembl115_annotation_metadata$ensembl_release[[1]],
    ", expected release ",
    ensembl_release,
    "."
  )
}

if (!identical(
  as.character(ensembl115_annotation_metadata$genome_assembly[[1]]),
  ensembl_genome_assembly
)) {
  stop(
    "Local annotation has genome assembly ",
    ensembl115_annotation_metadata$genome_assembly[[1]],
    ", expected ",
    ensembl_genome_assembly,
    "."
  )
}

if (!identical(
  as.character(ensembl115_annotation_metadata$retained_feature_type[[1]]),
  "gene"
)) {
  stop("Local Ensembl 115 annotation was not generated from gene records.")
}

if (!identical(
  as.character(ensembl115_annotation_metadata$retained_gene_biotype[[1]]),
  "protein_coding"
)) {
  stop("Local Ensembl 115 annotation is not restricted to protein-coding genes.")
}

ensembl115_protein_coding_annotation_all <- readRDS(
  ensembl115_protein_coding_rds
) |>
  tibble::as_tibble()

required_local_annotation_columns <- c(
  "ensembl_gene_id",
  "gene",
  "gene_biotype"
)

missing_local_annotation_columns <- setdiff(
  required_local_annotation_columns,
  colnames(ensembl115_protein_coding_annotation_all)
)

if (length(missing_local_annotation_columns) > 0) {
  stop(
    "Missing columns in local Ensembl 115 RDS annotation: ",
    paste(missing_local_annotation_columns, collapse = ", ")
  )
}

ensembl115_protein_coding_annotation_all <-
  ensembl115_protein_coding_annotation_all |>
  transmute(
    ensembl_gene_id = sub(
      "\\.[0-9]+$",
      "",
      trimws(as.character(ensembl_gene_id))
    ),
    gene = trimws(as.character(gene)),
    gene_biotype = trimws(as.character(gene_biotype))
  ) |>
  filter(
    !is.na(ensembl_gene_id),
    ensembl_gene_id != "",
    gene_biotype == "protein_coding"
  ) |>
  mutate(
    gene = if_else(
      is.na(gene) | gene == "",
      ensembl_gene_id,
      gene
    )
  ) |>
  arrange(
    ensembl_gene_id,
    desc(gene != ensembl_gene_id),
    gene
  ) |>
  distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )

if (nrow(ensembl115_protein_coding_annotation_all) == 0L) {
  stop("The local Ensembl 115 RDS contains no protein-coding genes.")
}

if (anyDuplicated(
  ensembl115_protein_coding_annotation_all$ensembl_gene_id
)) {
  stop("Duplicated Ensembl gene IDs remain in the local annotation.")
}

expected_local_protein_coding_genes <- as.integer(
  ensembl115_annotation_metadata$number_of_protein_coding_genes[[1]]
)

if (nrow(ensembl115_protein_coding_annotation_all) !=
    expected_local_protein_coding_genes) {
  stop(
    "The RDS contains ",
    nrow(ensembl115_protein_coding_annotation_all),
    " protein-coding genes, but its metadata reports ",
    expected_local_protein_coding_genes,
    "."
  )
}

bioconductor_version <- if (
  requireNamespace("BiocManager", quietly = TRUE)
) {
  as.character(BiocManager::version())
} else {
  NA_character_
}

reference_annotation_info <- tibble::tibble(
  r_version = R.version.string,
  bioconductor_version = bioconductor_version,
  edgeR_package_version = as.character(packageVersion("edgeR")),
  ensembl_release = ensembl_release,
  genome_assembly = ensembl_genome_assembly,
  annotation_source = "official Ensembl GTF saved locally",
  source_url = as.character(
    ensembl115_annotation_metadata$source_url[[1]]
  ),
  retained_feature_type = "gene",
  retained_gene_biotype = "protein_coding",
  genes_in_local_annotation =
    nrow(ensembl115_protein_coding_annotation_all),
  local_rds_file = normalizePath(
    ensembl115_protein_coding_rds,
    mustWork = TRUE
  ),
  local_metadata_file = normalizePath(
    ensembl115_annotation_metadata_file,
    mustWork = TRUE
  ),
  reference_annotation = reference_annotation_label
)

message("\nLocal Ensembl 115 annotation information:")
print(reference_annotation_info, width = Inf)


# ==============================================================================
# 17. Filter the pseudobulk matrix to local Ensembl 115 protein-coding genes
#
# This filter is applied before filterByExpr and before edgeR modelling.
# ==============================================================================
protein_coding_filter_status <- tibble::tibble(
  ensembl_gene_id = rownames(
    gene_counts_per_sample_raw_in_tissue_allGenes
  )
) |>
  left_join(
    ensembl115_protein_coding_annotation_all,
    by = "ensembl_gene_id"
  ) |>
  mutate(
    found_in_ensembl115_protein_coding_annotation =
      !is.na(gene_biotype),
    passed_protein_coding_filter =
      found_in_ensembl115_protein_coding_annotation &
      gene_biotype == "protein_coding",
    protein_coding_filter_status = case_when(
      passed_protein_coding_filter ~
        "passed_ensembl115_protein_coding_filter",
      TRUE ~
        "filtered_out_not_in_ensembl115_protein_coding_annotation"
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

if (length(protein_coding_gene_ids) == 0L) {
  stop(
    "No genes matched between the count matrix and the local ",
    "Ensembl 115 protein-coding annotation."
  )
}

gene_counts_per_sample_raw_in_tissue_proteinCoding <-
  gene_counts_per_sample_raw_in_tissue_allGenes[
    protein_coding_gene_ids,
    ,
    drop = FALSE
  ]

# Convenient alias: this is the protein-coding matrix used below.
gene_counts_per_sample_raw_in_tissue <-
  gene_counts_per_sample_raw_in_tissue_proteinCoding

protein_coding_gene_annotation <- protein_coding_filter_status |>
  filter(passed_protein_coding_filter) |>
  dplyr::select(
    ensembl_gene_id,
    gene,
    gene_biotype
  )

if (!identical(
  rownames(gene_counts_per_sample_raw_in_tissue_proteinCoding),
  protein_coding_gene_annotation$ensembl_gene_id
)) {
  stop("Protein-coding count matrix and annotation order do not match.")
}

message(
  "Genes before local Ensembl 115 protein-coding filter: ",
  format(
    nrow(gene_counts_per_sample_raw_in_tissue_allGenes),
    big.mark = ","
  )
)
message(
  "Protein-coding genes retained: ",
  format(
    nrow(gene_counts_per_sample_raw_in_tissue_proteinCoding),
    big.mark = ","
  )
)
message(
  "Genes removed by local protein-coding filter: ",
  format(
    nrow(gene_counts_per_sample_raw_in_tissue_allGenes) -
      nrow(gene_counts_per_sample_raw_in_tissue_proteinCoding),
    big.mark = ","
  )
)


# ==============================================================================

# ==============================================================================
# 18. Prepare the group-only edgeR design
# ==============================================================================

groupOnly_metadata <- sample_metadata_analysis |>
  dplyr::select(
    sample_ID,
    fmt_donor_group,
    sex
  ) |>
  dplyr::mutate(
    fmt_donor_group = factor(
      as.character(fmt_donor_group),
      levels = c("Neurotypical", "ASD")
    ),
    sex = factor(
      as.character(sex),
      levels = c("Male", "Female")
    )
  )

if (!identical(
  as.character(groupOnly_metadata$sample_ID),
  colnames(gene_counts_per_sample_raw_in_tissue)
)) {
  stop("Metadata order does not match the count-matrix column order.")
}

if (anyNA(groupOnly_metadata$fmt_donor_group)) {
  stop("Missing or unexpected fmt_donor_group values.")
}

group <- groupOnly_metadata$fmt_donor_group

group_counts <- table(group)

if (any(group_counts < 2L)) {
  stop(
    "At least one donor group has fewer than two samples: ",
    paste(
      names(group_counts),
      as.integer(group_counts),
      sep = "=",
      collapse = ", "
    )
  )
}

groupOnly_design <- model.matrix(
  ~ 0 + group
)

colnames(groupOnly_design) <- levels(group)
rownames(groupOnly_design) <- as.character(
  groupOnly_metadata$sample_ID
)

if (qr(groupOnly_design)$rank != ncol(groupOnly_design)) {
  stop("The group-only design matrix is not full rank.")
}

groupOnly_contrast <- limma::makeContrasts(
  ASD_vs_Neurotypical = ASD - Neurotypical,
  levels = groupOnly_design
)

message("\nGroup-only design matrix:")
print(groupOnly_design)

message("\nSample numbers by donor group:")
print(group_counts)


# ==============================================================================
# 19. Create DGEList and apply filterByExpr
# ==============================================================================

groupOnly_dge_unfiltered <- edgeR::DGEList(
  counts = gene_counts_per_sample_raw_in_tissue,
  group = group,
  samples = data.frame(
    sample_ID = as.character(groupOnly_metadata$sample_ID),
    fmt_donor_group = as.character(group),
    sex = as.character(groupOnly_metadata$sex),
    row.names = as.character(groupOnly_metadata$sample_ID),
    stringsAsFactors = FALSE
  )
)

groupOnly_cpm_unfiltered <- edgeR::cpm(
  groupOnly_dge_unfiltered,
  log = FALSE
)

groupOnly_keep_genes <- edgeR::filterByExpr(
  groupOnly_dge_unfiltered,
  design = groupOnly_design,
  min.count = filterByExpr_parameters$min.count,
  min.total.count = filterByExpr_parameters$min.total.count,
  large.n = filterByExpr_parameters$large.n,
  min.prop = filterByExpr_parameters$min.prop
)

if (
  !is.logical(groupOnly_keep_genes) ||
  length(groupOnly_keep_genes) !=
    nrow(groupOnly_dge_unfiltered$counts)
) {
  stop("Unexpected filterByExpr output.")
}

if (!any(groupOnly_keep_genes)) {
  stop("No genes passed filterByExpr.")
}

filterByExpr_status <- tibble::tibble(
  ensembl_gene_id = rownames(
    groupOnly_dge_unfiltered$counts
  ),
  total_count = rowSums(
    groupOnly_dge_unfiltered$counts
  ),
  mean_count = rowMeans(
    groupOnly_dge_unfiltered$counts
  ),
  samples_with_nonzero_count = rowSums(
    groupOnly_dge_unfiltered$counts > 0
  ),
  samples_with_count_ge_10 = rowSums(
    groupOnly_dge_unfiltered$counts >= 10
  ),
  mean_CPM = rowMeans(
    groupOnly_cpm_unfiltered
  ),
  max_CPM = apply(
    groupOnly_cpm_unfiltered,
    1,
    max
  ),
  passed_filterByExpr = groupOnly_keep_genes,
  filterByExpr_status = if_else(
    groupOnly_keep_genes,
    "passed_filterByExpr",
    "filtered_out_by_filterByExpr"
  )
) |>
  left_join(
    protein_coding_gene_annotation,
    by = "ensembl_gene_id"
  ) |>
  mutate(
    gene = if_else(
      is.na(gene) | gene == "",
      ensembl_gene_id,
      gene
    )
  ) |>
  dplyr::select(
    ensembl_gene_id,
    gene,
    gene_biotype,
    total_count,
    mean_count,
    samples_with_nonzero_count,
    samples_with_count_ge_10,
    mean_CPM,
    max_CPM,
    passed_filterByExpr,
    filterByExpr_status
  )

groupOnly_dge <- groupOnly_dge_unfiltered[
  groupOnly_keep_genes,
  ,
  keep.lib.sizes = FALSE
]

groupOnly_dge <- edgeR::calcNormFactors(
  groupOnly_dge,
  method = "TMM"
)

message(
  "\nGenes before filterByExpr: ",
  format(
    nrow(groupOnly_dge_unfiltered$counts),
    big.mark = ","
  )
)

message(
  "Genes retained after filterByExpr: ",
  format(
    nrow(groupOnly_dge$counts),
    big.mark = ","
  )
)

message(
  "Genes removed by filterByExpr: ",
  format(
    sum(!groupOnly_keep_genes),
    big.mark = ","
  )
)


# ==============================================================================
# 20. Estimate dispersions and run the edgeR quasi-likelihood test
# ==============================================================================

groupOnly_dge <- edgeR::estimateDisp(
  groupOnly_dge,
  design = groupOnly_design,
  robust = TRUE
)

groupOnly_fit <- edgeR::glmQLFit(
  groupOnly_dge,
  design = groupOnly_design,
  robust = TRUE
)

groupOnly_test <- edgeR::glmQLFTest(
  groupOnly_fit,
  contrast = groupOnly_contrast[
    ,
    "ASD_vs_Neurotypical"
  ]
)


# ==============================================================================
# 21. Create the simple full-results table
# ==============================================================================

groupOnly_fullResults <- edgeR::topTags(
  groupOnly_test,
  n = Inf,
  sort.by = "PValue"
)$table |>
  as.data.frame(
    stringsAsFactors = FALSE
  ) |>
  tibble::rownames_to_column(
    var = "ensembl_gene_id"
  ) |>
  left_join(
    protein_coding_gene_annotation,
    by = "ensembl_gene_id"
  ) |>
  mutate(
    gene = if_else(
      is.na(gene) | gene == "",
      ensembl_gene_id,
      gene
    ),
    regulation = if_else(
      logFC >= 0,
      "up",
      "down"
    )
  ) |>
  dplyr::select(
    ensembl_gene_id,
    gene,
    logFC,
    logCPM,
    F,
    PValue,
    FDR,
    regulation
  ) |>
  arrange(
    PValue,
    desc(abs(logFC))
  )

if (
  nrow(groupOnly_fullResults) !=
    nrow(groupOnly_dge$counts)
) {
  stop(
    "The number of result rows does not match ",
    "the number of tested genes."
  )
}

if (anyDuplicated(
  groupOnly_fullResults$ensembl_gene_id
)) {
  stop("Duplicated Ensembl gene IDs are present in full results.")
}


# ==============================================================================
# 22. Create complete mapping and filtering status tables
# ==============================================================================

protein_coding_mapping_status <- protein_coding_filter_status |>
  dplyr::select(
    ensembl_gene_id,
    gene,
    gene_biotype,
    found_in_ensembl115_protein_coding_annotation,
    passed_protein_coding_filter,
    protein_coding_filter_status
  )

all_gene_filter_status <- protein_coding_mapping_status |>
  left_join(
    filterByExpr_status |>
      dplyr::select(
        ensembl_gene_id,
        passed_filterByExpr,
        filterByExpr_status
      ),
    by = "ensembl_gene_id"
  ) |>
  mutate(
    passed_filterByExpr = if_else(
      is.na(passed_filterByExpr),
      FALSE,
      passed_filterByExpr
    ),
    final_filter_status = case_when(
      !passed_protein_coding_filter ~
        "removed_by_local_protein_coding_mapping",
      passed_protein_coding_filter &
        !passed_filterByExpr ~
        "removed_by_filterByExpr",
      passed_protein_coding_filter &
        passed_filterByExpr ~
        "tested_by_edgeR",
      TRUE ~
        "unexpected_filter_status"
    )
  ) |>
  dplyr::select(
    ensembl_gene_id,
    gene,
    gene_biotype,
    passed_protein_coding_filter,
    protein_coding_filter_status,
    passed_filterByExpr,
    filterByExpr_status,
    final_filter_status
  )


# ==============================================================================
# 23. Create filtering and analysis summaries
# ==============================================================================

genes_in_raw_matrix <- nrow(
  gene_counts_per_sample_raw_in_tissue_allGenes
)

genes_in_local_annotation <- nrow(
  ensembl115_protein_coding_annotation_all
)

genes_after_protein_coding_mapping <- nrow(
  gene_counts_per_sample_raw_in_tissue_proteinCoding
)

genes_removed_by_protein_coding_mapping <-
  genes_in_raw_matrix -
  genes_after_protein_coding_mapping

local_annotation_genes_absent_from_count_matrix <-
  length(
    setdiff(
      ensembl115_protein_coding_annotation_all$ensembl_gene_id,
      rownames(
        gene_counts_per_sample_raw_in_tissue_allGenes
      )
    )
  )

genes_before_filterByExpr <- nrow(
  groupOnly_dge_unfiltered$counts
)

genes_after_filterByExpr <- nrow(
  groupOnly_dge$counts
)

genes_removed_by_filterByExpr <-
  genes_before_filterByExpr -
  genes_after_filterByExpr

significant_FDR_0.05 <- sum(
  groupOnly_fullResults$FDR <= 0.05,
  na.rm = TRUE
)

significant_FDR_0.05_logFC_0.5 <- sum(
  groupOnly_fullResults$FDR <= 0.05 &
    abs(groupOnly_fullResults$logFC) > 0.5,
  na.rm = TRUE
)

filtering_summary <- tibble::tibble(
  metric = c(
    "samples_in_metadata",
    "samples_excluded_after_QC",
    "samples_retained_for_edgeR",
    "Neurotypical_samples",
    "ASD_samples",
    "genes_in_raw_count_matrix",
    "genes_in_local_Ensembl115_protein_coding_annotation",
    "genes_retained_after_protein_coding_mapping",
    "genes_removed_by_protein_coding_mapping",
    "local_annotation_genes_absent_from_count_matrix",
    "genes_before_filterByExpr",
    "genes_retained_after_filterByExpr",
    "genes_removed_by_filterByExpr",
    "genes_tested_by_edgeR",
    "significant_genes_FDR_0.05",
    "significant_genes_FDR_0.05_abs_logFC_gt_0.5",
    "filterByExpr_min_count",
    "filterByExpr_min_total_count",
    "filterByExpr_large_n",
    "filterByExpr_min_prop"
  ),
  value = c(
    nrow(metadata_autismFMT),
    length(excluded_samples),
    nrow(groupOnly_metadata),
    unname(group_counts["Neurotypical"]),
    unname(group_counts["ASD"]),
    genes_in_raw_matrix,
    genes_in_local_annotation,
    genes_after_protein_coding_mapping,
    genes_removed_by_protein_coding_mapping,
    local_annotation_genes_absent_from_count_matrix,
    genes_before_filterByExpr,
    genes_after_filterByExpr,
    genes_removed_by_filterByExpr,
    nrow(groupOnly_fullResults),
    significant_FDR_0.05,
    significant_FDR_0.05_logFC_0.5,
    filterByExpr_parameters$min.count,
    filterByExpr_parameters$min.total.count,
    filterByExpr_parameters$large.n,
    filterByExpr_parameters$min.prop
  )
)

analysis_summary <- tibble::tibble(
  model = "groupOnly",
  formula = "expression ~ fmt_donor_group",
  comparison = "ASD_vs_Neurotypical",
  control_group = "Neurotypical",
  experimental_group = "ASD",
  n_samples = nrow(groupOnly_metadata),
  n_Neurotypical = unname(
    group_counts["Neurotypical"]
  ),
  n_ASD = unname(
    group_counts["ASD"]
  ),
  genes_tested = nrow(groupOnly_fullResults),
  significant_FDR_0.05 = significant_FDR_0.05,
  significant_FDR_0.05_abs_logFC_gt_0.5 =
    significant_FDR_0.05_logFC_0.5,
  positive_logFC_interpretation =
    "higher expression in ASD than Neurotypical",
  negative_logFC_interpretation =
    "lower expression in ASD than Neurotypical",
  edgeR_version = as.character(
    packageVersion("edgeR")
  ),
  ensembl_release = ensembl_release,
  genome_assembly = ensembl_genome_assembly
)


# ==============================================================================
# 24. Prepare small supporting output tables
# ==============================================================================

groupOnly_design_output <- data.frame(
  sample_ID = rownames(groupOnly_design),
  groupOnly_design,
  row.names = NULL,
  check.names = FALSE
)

groupOnly_sample_metadata_output <- groupOnly_metadata |>
  mutate(
    fmt_donor_group = as.character(
      fmt_donor_group
    ),
    sex = as.character(
      sex
    )
  )


# ==============================================================================
# 25. Define output files
# ==============================================================================

output_files <- list(
  full_results = file.path(
    output_dir,
    "groupOnly_ASD_vs_Neurotypical_fullResults.tsv"
  ),
  analysis_summary = file.path(
    output_dir,
    "groupOnly_analysisSummary.tsv"
  ),
  filtering_summary = file.path(
    output_dir,
    "groupOnly_filteringSummary.tsv"
  ),
  all_gene_filter_status = file.path(
    output_dir,
    "groupOnly_allGeneFilterStatus.tsv"
  ),
  protein_coding_mapping_status = file.path(
    output_dir,
    "groupOnly_proteinCodingMappingStatus.tsv"
  ),
  filterByExpr_status = file.path(
    output_dir,
    "groupOnly_filterByExprStatus.tsv"
  ),
  sample_inclusion_status = file.path(
    output_dir,
    "groupOnly_sampleInclusionStatus.tsv"
  ),
  pseudobulk_input_summary = file.path(
    output_dir,
    "groupOnly_pseudobulkInputSummary.tsv"
  ),
  sample_metadata = file.path(
    output_dir,
    "groupOnly_sampleMetadata.tsv"
  ),
  design_matrix = file.path(
    output_dir,
    "groupOnly_designMatrix.tsv"
  ),
  reference_annotation_info = file.path(
    output_dir,
    "groupOnly_referenceAnnotationInfo.tsv"
  )
)


# ==============================================================================
# 26. Save all outputs as TSV files
# ==============================================================================

write_tsv_base <- function(
    object,
    output_file
) {
  write.table(
    object,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

write_tsv_base(
  groupOnly_fullResults,
  output_files$full_results
)

write_tsv_base(
  analysis_summary,
  output_files$analysis_summary
)

write_tsv_base(
  filtering_summary,
  output_files$filtering_summary
)

write_tsv_base(
  all_gene_filter_status,
  output_files$all_gene_filter_status
)

write_tsv_base(
  protein_coding_mapping_status,
  output_files$protein_coding_mapping_status
)

write_tsv_base(
  filterByExpr_status,
  output_files$filterByExpr_status
)

write_tsv_base(
  sample_inclusion_status,
  output_files$sample_inclusion_status
)

write_tsv_base(
  pseudobulk_input_summary,
  output_files$pseudobulk_input_summary
)

write_tsv_base(
  groupOnly_sample_metadata_output,
  output_files$sample_metadata
)

write_tsv_base(
  groupOnly_design_output,
  output_files$design_matrix
)

write_tsv_base(
  reference_annotation_info,
  output_files$reference_annotation_info
)


# ==============================================================================
# 27. Validate saved files
# ==============================================================================

missing_output_files <- unlist(
  output_files
)[
  !file.exists(
    unlist(output_files)
  )
]

if (length(missing_output_files) > 0L) {
  stop(
    "The following output files were not created:\n",
    paste(
      missing_output_files,
      collapse = "\n"
    )
  )
}

empty_output_files <- unlist(
  output_files
)[
  file.info(
    unlist(output_files)
  )$size == 0
]

if (length(empty_output_files) > 0L) {
  stop(
    "The following output files are empty:\n",
    paste(
      empty_output_files,
      collapse = "\n"
    )
  )
}


# ==============================================================================
# 28. Final report
# ==============================================================================

message(
  "\n",
  paste(rep("=", 80), collapse = "")
)

message("GROUP-ONLY edgeR ANALYSIS COMPLETED SUCCESSFULLY")

message(
  paste(rep("=", 80), collapse = "")
)

message("\nComparison: ASD vs Neurotypical")
message("Samples used: ", nrow(groupOnly_metadata))
message(
  "Genes in raw count matrix: ",
  format(genes_in_raw_matrix, big.mark = ",")
)
message(
  "Genes removed by protein-coding mapping: ",
  format(
    genes_removed_by_protein_coding_mapping,
    big.mark = ","
  )
)
message(
  "Genes removed by filterByExpr: ",
  format(
    genes_removed_by_filterByExpr,
    big.mark = ","
  )
)
message(
  "Genes tested by edgeR: ",
  format(
    nrow(groupOnly_fullResults),
    big.mark = ","
  )
)
message(
  "Significant genes at FDR <= 0.05: ",
  significant_FDR_0.05
)
message(
  "Significant genes at FDR <= 0.05 and |logFC| > 0.5: ",
  significant_FDR_0.05_logFC_0.5
)

message("\nFull edgeR result columns:")
message(
  paste(
    colnames(groupOnly_fullResults),
    collapse = ", "
  )
)

message("\nSaved files:")
for (output_name in names(output_files)) {
  message(
    "  ",
    output_name,
    ": ",
    output_files[[output_name]]
  )
}

message("\nMain result object in R:")
message("  groupOnly_fullResults")

message("\nTop 10 results:")
print(
  head(groupOnly_fullResults, 10),
  row.names = FALSE
)

# ==============================================================================
# End
# ==============================================================================
# groupOnly_fullResults %>% filter(FDR < 0.1) %>% 
#   filter(abs(logFC) > 0.5) %>%
#   .$regulation %>% table
