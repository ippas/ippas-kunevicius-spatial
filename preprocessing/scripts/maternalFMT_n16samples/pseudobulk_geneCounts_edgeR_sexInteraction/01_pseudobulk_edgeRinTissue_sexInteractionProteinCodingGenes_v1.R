#!/usr/bin/env Rscript

# ==============================================================================
# 02_pseudobulk_inTissue_edgeR_sexInteraction_proteinCodingGenes_Ensembl115.R
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
#   7. Apply edgeR::filterByExpr using the sex-interaction design.
#   8. Fit one edgeR quasi-likelihood model using all 16 retained samples:
#        expression ~ fmt_donor_group * sex
#   9. Test seven biologically useful contrasts:
#        - sex-by-group interaction;
#        - overall Female vs Male effect, averaged across donor groups;
#        - overall ASD vs Neurotypical effect, averaged across sexes;
#        - ASD vs Neurotypical among Male samples;
#        - ASD vs Neurotypical among Female samples;
#        - Female vs Male among Neurotypical samples;
#        - Female vs Male among ASD samples.
#  10. Save all seven complete edgeR result tables as separate TSV files.
#  11. Save six XLSX sheets:
#        1. Interaction
#        2. Overall_Sex
#        3. Overall_Group
#        4. Group_in_Male
#        5. Group_in_Female
#        6. Sex_within_Group, containing two full tables side by side.
#  12. Save raw pseudobulk counts for tested genes.
#  13. Save log2(raw count + 1) followed by quantile normalization.
#  14. Save mapping, filterByExpr and gene-filter-stage summaries.
#  15. Generate a human-readable README describing the model, all seven
#      comparisons, sample numbers, filtering counts and output files.
#
# INTERPRETATION OF logFC
#   Interaction:
#     positive logFC = the ASD-vs-Neurotypical effect is more positive in Female
#                      than in Male samples;
#     negative logFC = the ASD-vs-Neurotypical effect is more negative in Female
#                      than in Male samples.
#     The Interaction result intentionally has NO regulation column because
#     labels such as up/down would be biologically misleading for a difference
#     of differences.
#   Overall sex:
#     positive = higher expression in Female than Male, averaged across groups.
#   Overall group:
#     positive = higher expression in ASD than Neurotypical, averaged across sex.
#   Group-within-sex:
#     positive = higher expression in ASD than Neurotypical within that sex.
#   Sex-within-group:
#     positive = higher expression in Female than Male within that donor group.
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
#   maternalFMT_n16samples/pseudobulk_geneCounts_edgeR_sexInteraction/
#   proteinCodingGenes_allChromosomes
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
  "hdf5r",
  "openxlsx"
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
    "ask = FALSE, update = FALSE)\n",
    "install.packages(\"openxlsx\")"
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
  "pseudobulk_geneCounts_edgeR_sexInteraction",
  "proteinCodingGenes_allChromosomes"
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

# ==============================================================================

# ==============================================================================
# 18. Prepare interaction design and seven contrasts
#
# Model:
#   expression ~ fmt_donor_group * sex
#
# Reference levels:
#   fmt_donor_group = Neurotypical
#   sex = Male
#
# Renamed coefficients:
#   Intercept
#   groupASD
#   sexFemale
#   groupASD_sexFemale
#
# Cell means implied by the model:
#   Neurotypical Male = Intercept
#   ASD Male          = Intercept + groupASD
#   Neurotypical Female =
#       Intercept + sexFemale
#   ASD Female =
#       Intercept + groupASD + sexFemale + groupASD_sexFemale
#
# Overall effects use equal-weight marginal means:
#   Overall sex =
#       mean(Female cells) - mean(Male cells)
#     = sexFemale + 0.5 * groupASD_sexFemale
#
#   Overall group =
#       mean(ASD cells) - mean(Neurotypical cells)
#     = groupASD + 0.5 * groupASD_sexFemale
#
# In this dataset, sex and group proportions are balanced across the other
# factor, so these marginal contrasts also represent the requested comparisons
# of all Female vs all Male and all ASD vs all Neurotypical samples.
# ==============================================================================

sexInteraction_metadata <- sample_metadata_analysis |>
  dplyr::select(
    sample_ID,
    fmt_donor_group,
    sex
  ) |>
  dplyr::mutate(
    fmt_donor_group = factor(
      as.character(fmt_donor_group),
      levels = c(
        "Neurotypical",
        "ASD"
      )
    ),
    sex = factor(
      as.character(sex),
      levels = c(
        "Male",
        "Female"
      )
    )
  )

if (!identical(
  as.character(
    sexInteraction_metadata$sample_ID
  ),
  colnames(
    gene_counts_per_sample_raw_in_tissue
  )
)) {
  stop(
    "Metadata and count-matrix sample orders do not match."
  )
}

if (
  anyNA(
    sexInteraction_metadata$fmt_donor_group
  ) ||
  anyNA(
    sexInteraction_metadata$sex
  )
) {
  stop(
    "Missing or unexpected fmt_donor_group/sex values."
  )
}

sex_group_counts <- table(
  sexInteraction_metadata$sex,
  sexInteraction_metadata$fmt_donor_group
)

if (any(sex_group_counts < 2L)) {
  stop(
    "At least one sex-by-group cell has fewer than two samples."
  )
}

sexInteraction_design <- model.matrix(
  ~ fmt_donor_group * sex,
  data = sexInteraction_metadata
)

rownames(
  sexInteraction_design
) <- as.character(
  sexInteraction_metadata$sample_ID
)

expected_original_columns <- c(
  "(Intercept)",
  "fmt_donor_groupASD",
  "sexFemale",
  "fmt_donor_groupASD:sexFemale"
)

if (!identical(
  colnames(sexInteraction_design),
  expected_original_columns
)) {
  stop(
    "Unexpected design columns: ",
    paste(
      colnames(sexInteraction_design),
      collapse = ", "
    )
  )
}

colnames(
  sexInteraction_design
) <- c(
  "Intercept",
  "groupASD",
  "sexFemale",
  "groupASD_sexFemale"
)

if (
  qr(sexInteraction_design)$rank !=
    ncol(sexInteraction_design)
) {
  stop(
    "The interaction design matrix is not full rank."
  )
}

sexInteraction_contrasts <- limma::makeContrasts(

  # Difference between the ASD effect in Female and Male samples.
  Interaction =
    groupASD_sexFemale,

  # Female vs Male averaged equally across Neurotypical and ASD.
  Overall_Sex_Female_vs_Male =
    sexFemale +
    0.5 * groupASD_sexFemale,

  # ASD vs Neurotypical averaged equally across Male and Female.
  Overall_Group_ASD_vs_Neurotypical =
    groupASD +
    0.5 * groupASD_sexFemale,

  # Group effects within each sex.
  ASD_vs_Neurotypical_in_Male =
    groupASD,

  ASD_vs_Neurotypical_in_Female =
    groupASD +
    groupASD_sexFemale,

  # Sex effects within each donor group.
  Female_vs_Male_in_Neurotypical =
    sexFemale,

  Female_vs_Male_in_ASD =
    sexFemale +
    groupASD_sexFemale,

  levels = sexInteraction_design
)

expected_contrast_names <- c(
  "Interaction",
  "Overall_Sex_Female_vs_Male",
  "Overall_Group_ASD_vs_Neurotypical",
  "ASD_vs_Neurotypical_in_Male",
  "ASD_vs_Neurotypical_in_Female",
  "Female_vs_Male_in_Neurotypical",
  "Female_vs_Male_in_ASD"
)

if (!identical(
  colnames(sexInteraction_contrasts),
  expected_contrast_names
)) {
  stop(
    "Unexpected contrast names: ",
    paste(
      colnames(sexInteraction_contrasts),
      collapse = ", "
    )
  )
}

message("\nInteraction design matrix:")
print(sexInteraction_design)

message("\nSeven contrast vectors:")
print(sexInteraction_contrasts)

message("\nSample numbers by sex and donor group:")
print(sex_group_counts)


# ==============================================================================
# 19. Create DGEList and apply filterByExpr
# ==============================================================================

sexInteraction_dge_unfiltered <- edgeR::DGEList(
  counts =
    gene_counts_per_sample_raw_in_tissue,
  samples = data.frame(
    sample_ID = as.character(
      sexInteraction_metadata$sample_ID
    ),
    fmt_donor_group = as.character(
      sexInteraction_metadata$fmt_donor_group
    ),
    sex = as.character(
      sexInteraction_metadata$sex
    ),
    row.names = as.character(
      sexInteraction_metadata$sample_ID
    ),
    stringsAsFactors = FALSE
  )
)

unfiltered_cpm <- edgeR::cpm(
  sexInteraction_dge_unfiltered,
  log = FALSE
)

sexInteraction_keep_genes <- edgeR::filterByExpr(
  sexInteraction_dge_unfiltered,
  design = sexInteraction_design,
  min.count =
    filterByExpr_parameters$min.count,
  min.total.count =
    filterByExpr_parameters$min.total.count,
  large.n =
    filterByExpr_parameters$large.n,
  min.prop =
    filterByExpr_parameters$min.prop
)

if (
  !is.logical(sexInteraction_keep_genes) ||
  length(sexInteraction_keep_genes) !=
    nrow(
      sexInteraction_dge_unfiltered$counts
    )
) {
  stop("Unexpected filterByExpr output.")
}

if (!any(sexInteraction_keep_genes)) {
  stop("No genes passed filterByExpr.")
}

filterByExpr_status <- tibble::tibble(
  ensembl_gene_id = rownames(
    sexInteraction_dge_unfiltered$counts
  ),
  total_count = rowSums(
    sexInteraction_dge_unfiltered$counts
  ),
  mean_count = rowMeans(
    sexInteraction_dge_unfiltered$counts
  ),
  samples_with_nonzero_count = rowSums(
    sexInteraction_dge_unfiltered$counts > 0
  ),
  samples_with_count_ge_10 = rowSums(
    sexInteraction_dge_unfiltered$counts >= 10
  ),
  mean_CPM = rowMeans(
    unfiltered_cpm
  ),
  max_CPM = apply(
    unfiltered_cpm,
    1,
    max
  ),
  passed_filterByExpr =
    sexInteraction_keep_genes,
  filterByExpr_status = dplyr::if_else(
    sexInteraction_keep_genes,
    "passed_filterByExpr",
    "filtered_out_by_filterByExpr"
  )
) |>
  dplyr::left_join(
    protein_coding_gene_annotation,
    by = "ensembl_gene_id"
  ) |>
  dplyr::mutate(
    gene = dplyr::if_else(
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

sexInteraction_dge <-
  sexInteraction_dge_unfiltered[
    sexInteraction_keep_genes,
    ,
    keep.lib.sizes = FALSE
  ]

sexInteraction_dge <- edgeR::calcNormFactors(
  sexInteraction_dge,
  method = "TMM"
)

message(
  "\nGenes before filterByExpr: ",
  format(
    nrow(
      sexInteraction_dge_unfiltered$counts
    ),
    big.mark = ","
  )
)

message(
  "Genes retained after filterByExpr: ",
  format(
    nrow(
      sexInteraction_dge$counts
    ),
    big.mark = ","
  )
)

message(
  "Genes removed by filterByExpr: ",
  format(
    sum(!sexInteraction_keep_genes),
    big.mark = ","
  )
)


# ==============================================================================
# 20. Fit edgeR quasi-likelihood model
# ==============================================================================

sexInteraction_dge <- edgeR::estimateDisp(
  sexInteraction_dge,
  design = sexInteraction_design,
  robust = TRUE
)

sexInteraction_fit <- edgeR::glmQLFit(
  sexInteraction_dge,
  design = sexInteraction_design,
  robust = TRUE
)


# ==============================================================================
# 21. Run all seven tests
# ==============================================================================

run_contrast_test <- function(
    contrast_name
) {

  edgeR::glmQLFTest(
    sexInteraction_fit,
    contrast =
      sexInteraction_contrasts[
        ,
        contrast_name
      ]
  )
}

sexInteraction_test_objects <- lapply(
  expected_contrast_names,
  run_contrast_test
)

names(
  sexInteraction_test_objects
) <- expected_contrast_names


# ==============================================================================
# 22. Create simple full-result tables
# ==============================================================================

create_simple_result <- function(
    qlf_test
) {

  result <- edgeR::topTags(
    qlf_test,
    n = Inf,
    sort.by = "PValue"
  )$table |>
    as.data.frame(
      stringsAsFactors = FALSE
    ) |>
    tibble::rownames_to_column(
      var = "ensembl_gene_id"
    ) |>
    dplyr::left_join(
      protein_coding_gene_annotation,
      by = "ensembl_gene_id"
    ) |>
    dplyr::mutate(
      gene = dplyr::if_else(
        is.na(gene) | gene == "",
        ensembl_gene_id,
        gene
      ),
      regulation = dplyr::if_else(
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
    dplyr::arrange(
      PValue,
      dplyr::desc(abs(logFC))
    )

  expected_result_columns <- c(
    "ensembl_gene_id",
    "gene",
    "logFC",
    "logCPM",
    "F",
    "PValue",
    "FDR",
    "regulation"
  )

  if (!identical(
    colnames(result),
    expected_result_columns
  )) {
    stop(
      "Unexpected edgeR result columns."
    )
  }

  if (anyDuplicated(
    result$ensembl_gene_id
  )) {
    stop(
      "Duplicated Ensembl gene IDs in an edgeR result."
    )
  }

  result
}

sexInteraction_results <- lapply(
  sexInteraction_test_objects,
  create_simple_result
)

tested_genes_n <- nrow(
  sexInteraction_dge$counts
)

for (
  result_name in names(
    sexInteraction_results
  )
) {
  if (
    nrow(
      sexInteraction_results[[result_name]]
    ) != tested_genes_n
  ) {
    stop(
      "Unexpected row count for result: ",
      result_name
    )
  }
}

# The interaction is a difference of differences. Its logFC sign is meaningful,
# but calling it up/down would be misleading. Therefore, remove regulation only
# from the interaction table. The other six result tables retain regulation.
sexInteraction_results$Interaction <-
  sexInteraction_results$Interaction |>
  dplyr::select(
    -regulation
  )

# Convenient aliases.
sexInteraction_fullResults_interaction <-
  sexInteraction_results$Interaction

sexInteraction_fullResults_overallSex <-
  sexInteraction_results$Overall_Sex_Female_vs_Male

sexInteraction_fullResults_overallGroup <-
  sexInteraction_results$Overall_Group_ASD_vs_Neurotypical

sexInteraction_fullResults_groupMale <-
  sexInteraction_results$ASD_vs_Neurotypical_in_Male

sexInteraction_fullResults_groupFemale <-
  sexInteraction_results$ASD_vs_Neurotypical_in_Female

sexInteraction_fullResults_sexNeurotypical <-
  sexInteraction_results$Female_vs_Male_in_Neurotypical

sexInteraction_fullResults_sexASD <-
  sexInteraction_results$Female_vs_Male_in_ASD


# ==============================================================================
# 23. Create raw and log2 + quantile-normalized expression matrices
# ==============================================================================

sexInteraction_rawCounts_testedGenes <-
  sexInteraction_dge$counts

sexInteraction_log2RawCounts <- log2(
  sexInteraction_rawCounts_testedGenes + 1
)

sexInteraction_log2RawCounts_QN <-
  limma::normalizeBetweenArrays(
    sexInteraction_log2RawCounts,
    method = "quantile"
  )

rownames(
  sexInteraction_log2RawCounts_QN
) <- rownames(
  sexInteraction_rawCounts_testedGenes
)

colnames(
  sexInteraction_log2RawCounts_QN
) <- colnames(
  sexInteraction_rawCounts_testedGenes
)

matrix_to_gene_table <- function(
    expression_matrix
) {

  matched_annotation <-
    protein_coding_gene_annotation[
      match(
        rownames(expression_matrix),
        protein_coding_gene_annotation$ensembl_gene_id
      ),
      ,
      drop = FALSE
    ]

  if (anyNA(
    matched_annotation$ensembl_gene_id
  )) {
    stop(
      "Expression genes could not be matched to annotation."
    )
  }

  tibble::tibble(
    ensembl_gene_id =
      matched_annotation$ensembl_gene_id,
    gene =
      matched_annotation$gene
  ) |>
    dplyr::bind_cols(
      tibble::as_tibble(
        as.data.frame(
          expression_matrix,
          check.names = FALSE
        ),
        .name_repair = "minimal"
      )
    )
}

sexInteraction_rawCounts_output <-
  matrix_to_gene_table(
    sexInteraction_rawCounts_testedGenes
  )

sexInteraction_log2RawCounts_QN_output <-
  matrix_to_gene_table(
    sexInteraction_log2RawCounts_QN
  )


# ==============================================================================
# 24. Create filtering-status tables
# ==============================================================================

protein_coding_mapping_status <-
  protein_coding_filter_status |>
  dplyr::select(
    ensembl_gene_id,
    gene,
    gene_biotype,
    found_in_ensembl115_protein_coding_annotation,
    passed_protein_coding_filter,
    protein_coding_filter_status
  )

all_gene_filter_status <-
  protein_coding_mapping_status |>
  dplyr::left_join(
    filterByExpr_status |>
      dplyr::select(
        ensembl_gene_id,
        passed_filterByExpr,
        filterByExpr_status
      ),
    by = "ensembl_gene_id"
  ) |>
  dplyr::mutate(
    passed_filterByExpr =
      dplyr::if_else(
        is.na(passed_filterByExpr),
        FALSE,
        passed_filterByExpr
      ),
    final_filter_status =
      dplyr::case_when(
        !passed_protein_coding_filter ~
          paste0(
            "removed_not_in_local_Ensembl115_",
            "protein_coding_annotation"
          ),
        passed_protein_coding_filter &
          !passed_filterByExpr ~
          "removed_by_filterByExpr",
        passed_protein_coding_filter &
          passed_filterByExpr ~
          "tested_by_edgeR",
        TRUE ~
          "unexpected_filter_status"
      )
  )

gene_filtering_stage <-
  all_gene_filter_status |>
  dplyr::transmute(
    gene,
    filter_stage =
      final_filter_status
  ) |>
  dplyr::arrange(
    filter_stage,
    gene
  )


# ==============================================================================
# 25. Create test definitions and summaries
# ==============================================================================

test_definitions <- tibble::tribble(
  ~test, ~comparison, ~positive_logFC_meaning, ~negative_logFC_meaning,

  "Interaction",
  "(ASD - Neurotypical)_Female - (ASD - Neurotypical)_Male",
  "positive interaction logFC: ASD effect is more positive in Female than Male",
  "negative interaction logFC: ASD effect is more negative in Female than Male",

  "Overall_Sex_Female_vs_Male",
  "Female vs Male averaged equally across donor groups",
  "higher expression in Female than Male",
  "lower expression in Female than Male",

  "Overall_Group_ASD_vs_Neurotypical",
  "ASD vs Neurotypical averaged equally across sexes",
  "higher expression in ASD than Neurotypical",
  "lower expression in ASD than Neurotypical",

  "ASD_vs_Neurotypical_in_Male",
  "ASD vs Neurotypical within Male",
  "higher expression in ASD males",
  "lower expression in ASD males",

  "ASD_vs_Neurotypical_in_Female",
  "ASD vs Neurotypical within Female",
  "higher expression in ASD females",
  "lower expression in ASD females",

  "Female_vs_Male_in_Neurotypical",
  "Female vs Male within Neurotypical",
  "higher expression in Neurotypical females",
  "lower expression in Neurotypical females",

  "Female_vs_Male_in_ASD",
  "Female vs Male within ASD",
  "higher expression in ASD females",
  "lower expression in ASD females"
)

summarize_one_result <- function(
    test_name
) {

  result_table <-
    sexInteraction_results[[test_name]]

  test_definition <-
    test_definitions |>
    dplyr::filter(
      test == test_name
    )

  tibble::tibble(
    model =
      "sexInteraction",
    formula =
      "expression ~ fmt_donor_group * sex",
    test =
      test_name,
    comparison =
      test_definition$comparison,
    tested_genes =
      nrow(result_table),
    significant_P_0.05 =
      sum(
        result_table$PValue < 0.05,
        na.rm = TRUE
      ),
    significant_FDR_0.05 =
      sum(
        result_table$FDR < 0.05,
        na.rm = TRUE
      ),
    significant_FDR_0.05_abs_logFC_gt_0.5 =
      sum(
        result_table$FDR < 0.05 &
          abs(result_table$logFC) > 0.5,
        na.rm = TRUE
      ),
    significant_FDR_0.1 =
      sum(
        result_table$FDR < 0.1,
        na.rm = TRUE
      ),
    significant_FDR_0.1_abs_logFC_gt_0.5 =
      sum(
        result_table$FDR < 0.1 &
          abs(result_table$logFC) > 0.5,
        na.rm = TRUE
      ),
    positive_logFC_interpretation =
      test_definition$positive_logFC_meaning,
    negative_logFC_interpretation =
      test_definition$negative_logFC_meaning
  )
}

edgeR_analysis_summary <-
  dplyr::bind_rows(
    lapply(
      expected_contrast_names,
      summarize_one_result
    )
  )

genes_in_raw_matrix <- nrow(
  gene_counts_per_sample_raw_in_tissue_allGenes
)

genes_after_mapping <- nrow(
  gene_counts_per_sample_raw_in_tissue_proteinCoding
)

genes_removed_by_mapping <-
  genes_in_raw_matrix -
  genes_after_mapping

genes_before_filterByExpr <- nrow(
  sexInteraction_dge_unfiltered$counts
)

genes_after_filterByExpr <- nrow(
  sexInteraction_dge$counts
)

genes_removed_by_filterByExpr <-
  genes_before_filterByExpr -
  genes_after_filterByExpr

filtering_summary <- tibble::tibble(
  metric = c(
    "samples_in_metadata",
    "samples_excluded_after_QC",
    "samples_retained_for_edgeR",
    "Male_Neurotypical_samples",
    "Male_ASD_samples",
    "Female_Neurotypical_samples",
    "Female_ASD_samples",
    "genes_in_raw_count_matrix",
    paste0(
      "genes_in_local_Ensembl115_",
      "protein_coding_annotation"
    ),
    "genes_retained_after_protein_coding_mapping",
    "genes_removed_by_protein_coding_mapping",
    "genes_before_filterByExpr",
    "genes_retained_after_filterByExpr",
    "genes_removed_by_filterByExpr",
    "genes_tested_by_edgeR",
    "number_of_edgeR_tests",
    "number_of_XLSX_sheets",
    "filterByExpr_min_count",
    "filterByExpr_min_total_count",
    "filterByExpr_large_n",
    "filterByExpr_min_prop"
  ),
  value = c(
    nrow(metadata_autismFMT),
    length(excluded_samples),
    nrow(sexInteraction_metadata),
    unname(
      sex_group_counts[
        "Male",
        "Neurotypical"
      ]
    ),
    unname(
      sex_group_counts[
        "Male",
        "ASD"
      ]
    ),
    unname(
      sex_group_counts[
        "Female",
        "Neurotypical"
      ]
    ),
    unname(
      sex_group_counts[
        "Female",
        "ASD"
      ]
    ),
    genes_in_raw_matrix,
    nrow(
      ensembl115_protein_coding_annotation_all
    ),
    genes_after_mapping,
    genes_removed_by_mapping,
    genes_before_filterByExpr,
    genes_after_filterByExpr,
    genes_removed_by_filterByExpr,
    tested_genes_n,
    length(expected_contrast_names),
    6L,
    filterByExpr_parameters$min.count,
    filterByExpr_parameters$min.total.count,
    filterByExpr_parameters$large.n,
    filterByExpr_parameters$min.prop
  )
)


# ==============================================================================
# 26. Prepare supporting output tables
# ==============================================================================

design_output <- data.frame(
  sample_ID =
    rownames(sexInteraction_design),
  sexInteraction_design,
  row.names = NULL,
  check.names = FALSE
)

contrast_output <- data.frame(
  coefficient =
    rownames(sexInteraction_contrasts),
  sexInteraction_contrasts,
  row.names = NULL,
  check.names = FALSE
)

sample_metadata_output <-
  sexInteraction_metadata |>
  dplyr::mutate(
    fmt_donor_group =
      as.character(fmt_donor_group),
    sex =
      as.character(sex)
  )


# ==============================================================================
# 27. Define output files
# ==============================================================================

result_tsv_files <- c(
  Interaction = file.path(
    output_dir,
    "sexInteraction_interaction_fullResults.tsv"
  ),
  Overall_Sex_Female_vs_Male = file.path(
    output_dir,
    "sexInteraction_overall_Female_vs_Male_fullResults.tsv"
  ),
  Overall_Group_ASD_vs_Neurotypical = file.path(
    output_dir,
    "sexInteraction_overall_ASD_vs_Neurotypical_fullResults.tsv"
  ),
  ASD_vs_Neurotypical_in_Male = file.path(
    output_dir,
    "sexInteraction_ASD_vs_Neurotypical_in_Male_fullResults.tsv"
  ),
  ASD_vs_Neurotypical_in_Female = file.path(
    output_dir,
    "sexInteraction_ASD_vs_Neurotypical_in_Female_fullResults.tsv"
  ),
  Female_vs_Male_in_Neurotypical = file.path(
    output_dir,
    "sexInteraction_Female_vs_Male_in_Neurotypical_fullResults.tsv"
  ),
  Female_vs_Male_in_ASD = file.path(
    output_dir,
    "sexInteraction_Female_vs_Male_in_ASD_fullResults.tsv"
  )
)

output_files <- list(
  results_xlsx = file.path(
    output_dir,
    "sexInteraction_fullResults_sixSheets.xlsx"
  ),
  raw_counts = file.path(
    output_dir,
    "sexInteraction_rawCounts_testedGenes.tsv"
  ),
  log2_qn = file.path(
    output_dir,
    paste0(
      "sexInteraction_log2RawCounts_",
      "quantileNormalized_testedGenes.tsv"
    )
  ),
  filtering_stage = file.path(
    output_dir,
    "sexInteraction_geneFilteringStage.tsv"
  ),
  analysis_summary = file.path(
    output_dir,
    "sexInteraction_edgeRAnalysisSummary.tsv"
  ),
  test_definitions = file.path(
    output_dir,
    "sexInteraction_testDefinitions.tsv"
  ),
  filtering_summary = file.path(
    output_dir,
    "sexInteraction_filteringSummary.tsv"
  ),
  all_gene_filter_status = file.path(
    output_dir,
    "sexInteraction_allGeneFilterStatus.tsv"
  ),
  mapping_status = file.path(
    output_dir,
    "sexInteraction_proteinCodingMappingStatus.tsv"
  ),
  filterByExpr_status = file.path(
    output_dir,
    "sexInteraction_filterByExprStatus.tsv"
  ),
  sample_inclusion = file.path(
    output_dir,
    "sexInteraction_sampleInclusionStatus.tsv"
  ),
  pseudobulk_summary = file.path(
    output_dir,
    "sexInteraction_pseudobulkInputSummary.tsv"
  ),
  sample_metadata = file.path(
    output_dir,
    "sexInteraction_sampleMetadata.tsv"
  ),
  design = file.path(
    output_dir,
    "sexInteraction_designMatrix.tsv"
  ),
  contrasts = file.path(
    output_dir,
    "sexInteraction_contrastMatrix.tsv"
  ),
  reference_info = file.path(
    output_dir,
    "sexInteraction_referenceAnnotationInfo.tsv"
  ),
  readme = file.path(
    output_dir,
    "README_sexInteraction_analysis.txt"
  )
)


# ==============================================================================
# 28. Save all seven result TSV files and supporting TSV files
# ==============================================================================

write_tsv <- function(
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

for (
  test_name in names(result_tsv_files)
) {
  write_tsv(
    sexInteraction_results[[test_name]],
    result_tsv_files[[test_name]]
  )
}

write_tsv(
  sexInteraction_rawCounts_output,
  output_files$raw_counts
)

write_tsv(
  sexInteraction_log2RawCounts_QN_output,
  output_files$log2_qn
)

write_tsv(
  gene_filtering_stage,
  output_files$filtering_stage
)

write_tsv(
  edgeR_analysis_summary,
  output_files$analysis_summary
)

write_tsv(
  test_definitions,
  output_files$test_definitions
)

write_tsv(
  filtering_summary,
  output_files$filtering_summary
)

write_tsv(
  all_gene_filter_status,
  output_files$all_gene_filter_status
)

write_tsv(
  protein_coding_mapping_status,
  output_files$mapping_status
)

write_tsv(
  filterByExpr_status,
  output_files$filterByExpr_status
)

write_tsv(
  sample_inclusion_status,
  output_files$sample_inclusion
)

write_tsv(
  pseudobulk_input_summary,
  output_files$pseudobulk_summary
)

write_tsv(
  sample_metadata_output,
  output_files$sample_metadata
)

write_tsv(
  design_output,
  output_files$design
)

write_tsv(
  contrast_output,
  output_files$contrasts
)

write_tsv(
  reference_annotation_info,
  output_files$reference_info
)


# ==============================================================================
# 29. Create one XLSX workbook with exactly six sheets
#
# Sheets 1-5 each contain one complete result table.
#
# Sheet 6 contains two complete sex-within-group tables side by side:
#   columns A:H  = Female vs Male in Neurotypical
#   columns J:Q  = Female vs Male in ASD
# ==============================================================================

workbook <- openxlsx::createWorkbook()

workbook_sheet_names <- c(
  "Interaction",
  "Overall_Sex",
  "Overall_Group",
  "Group_in_Male",
  "Group_in_Female",
  "Sex_within_Group"
)

for (
  sheet_name in workbook_sheet_names
) {
  openxlsx::addWorksheet(
    workbook,
    sheetName = sheet_name
  )
}

single_table_sheets <- list(
  Interaction =
    sexInteraction_fullResults_interaction,
  Overall_Sex =
    sexInteraction_fullResults_overallSex,
  Overall_Group =
    sexInteraction_fullResults_overallGroup,
  Group_in_Male =
    sexInteraction_fullResults_groupMale,
  Group_in_Female =
    sexInteraction_fullResults_groupFemale
)

for (
  sheet_name in names(
    single_table_sheets
  )
) {
  openxlsx::writeDataTable(
    workbook,
    sheet = sheet_name,
    x = single_table_sheets[[sheet_name]],
    startRow = 1,
    startCol = 1,
    withFilter = TRUE,
    tableStyle = "TableStyleMedium2"
  )

  openxlsx::freezePane(
    workbook,
    sheet = sheet_name,
    firstRow = TRUE
  )

  openxlsx::setColWidths(
    workbook,
    sheet = sheet_name,
    cols = seq_len(
      ncol(
        single_table_sheets[[sheet_name]]
      )
    ),
    widths = "auto"
  )
}

# Sixth sheet: two complete tables side by side.
openxlsx::writeData(
  workbook,
  sheet = "Sex_within_Group",
  x = "Female vs Male in Neurotypical",
  startRow = 1,
  startCol = 1
)

openxlsx::writeDataTable(
  workbook,
  sheet = "Sex_within_Group",
  x =
    sexInteraction_fullResults_sexNeurotypical,
  startRow = 3,
  startCol = 1,
  withFilter = TRUE,
  tableStyle = "TableStyleMedium2",
  tableName = "SexNeurotypical"
)

openxlsx::writeData(
  workbook,
  sheet = "Sex_within_Group",
  x = "Female vs Male in ASD",
  startRow = 1,
  startCol = 10
)

openxlsx::writeDataTable(
  workbook,
  sheet = "Sex_within_Group",
  x =
    sexInteraction_fullResults_sexASD,
  startRow = 3,
  startCol = 10,
  withFilter = TRUE,
  tableStyle = "TableStyleMedium4",
  tableName = "SexASD"
)

openxlsx::freezePane(
  workbook,
  sheet = "Sex_within_Group",
  firstActiveRow = 3,
  firstActiveCol = 1
)

openxlsx::setColWidths(
  workbook,
  sheet = "Sex_within_Group",
  cols = c(
    1:8,
    10:17
  ),
  widths = "auto"
)

openxlsx::saveWorkbook(
  workbook,
  file =
    output_files$results_xlsx,
  overwrite = TRUE
)


# ==============================================================================
# 30. Generate human-readable README
#
# filterByExpr is performed exactly ONCE using the complete interaction design:
#   expression ~ fmt_donor_group * sex
#
# The same retained gene set is used for all seven contrasts. There is no
# contrast-specific filterByExpr step.
# ==============================================================================

format_integer_for_readme <- function(
    value
) {
  format(
    as.integer(value),
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

comparison_readme_table <- edgeR_analysis_summary |>
  dplyr::select(
    test,
    comparison,
    tested_genes,
    significant_P_0.05,
    significant_FDR_0.05,
    significant_FDR_0.05_abs_logFC_gt_0.5,
    significant_FDR_0.1,
    significant_FDR_0.1_abs_logFC_gt_0.5
  )

comparison_readme_header <- paste(
  colnames(comparison_readme_table),
  collapse = "\t"
)

comparison_readme_rows <- apply(
  comparison_readme_table,
  1,
  function(row_values) {
    paste(
      row_values,
      collapse = "\t"
    )
  }
)

readme_lines <- c(
  "Maternal FMT pseudobulk edgeR analysis: sex-by-group interaction",
  paste(rep("=", 72), collapse = ""),
  "",
  paste0(
    "Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  paste0(
    "Output directory: ",
    normalizePath(
      output_dir,
      mustWork = TRUE
    )
  ),
  "",
  "1. ANALYSIS MODEL",
  "-----------------",
  "Model:",
  "  expression ~ fmt_donor_group * sex",
  "",
  "Expanded model:",
  "  expression ~ fmt_donor_group + sex + fmt_donor_group:sex",
  "",
  "Reference levels:",
  "  fmt_donor_group = Neurotypical",
  "  sex = Male",
  "",
  "The model was fitted once using all 16 retained samples.",
  "",
  "2. SAMPLE NUMBERS",
  "-----------------",
  paste0(
    "Samples in metadata: ",
    format_integer_for_readme(
      nrow(metadata_autismFMT)
    )
  ),
  paste0(
    "Samples excluded after QC: ",
    format_integer_for_readme(
      length(excluded_samples)
    )
  ),
  paste0(
    "Samples retained for edgeR: ",
    format_integer_for_readme(
      nrow(sexInteraction_metadata)
    )
  ),
  "",
  paste0(
    "Male Neurotypical: ",
    format_integer_for_readme(
      sex_group_counts[
        "Male",
        "Neurotypical"
      ]
    )
  ),
  paste0(
    "Male ASD: ",
    format_integer_for_readme(
      sex_group_counts[
        "Male",
        "ASD"
      ]
    )
  ),
  paste0(
    "Female Neurotypical: ",
    format_integer_for_readme(
      sex_group_counts[
        "Female",
        "Neurotypical"
      ]
    )
  ),
  paste0(
    "Female ASD: ",
    format_integer_for_readme(
      sex_group_counts[
        "Female",
        "ASD"
      ]
    )
  ),
  "",
  "3. GENE FILTERING",
  "-----------------",
  paste0(
    "Genes in original raw count matrix: ",
    format_integer_for_readme(
      genes_in_raw_matrix
    )
  ),
  paste0(
    "Genes in local Ensembl 115 protein-coding annotation: ",
    format_integer_for_readme(
      nrow(
        ensembl115_protein_coding_annotation_all
      )
    )
  ),
  paste0(
    "Genes retained after protein-coding mapping: ",
    format_integer_for_readme(
      genes_after_mapping
    )
  ),
  paste0(
    "Genes removed during protein-coding mapping: ",
    format_integer_for_readme(
      genes_removed_by_mapping
    )
  ),
  "",
  paste0(
    "Genes entering filterByExpr: ",
    format_integer_for_readme(
      genes_before_filterByExpr
    )
  ),
  paste0(
    "Genes retained by filterByExpr: ",
    format_integer_for_readme(
      genes_after_filterByExpr
    )
  ),
  paste0(
    "Genes removed by filterByExpr: ",
    format_integer_for_readme(
      genes_removed_by_filterByExpr
    )
  ),
  "",
  "filterByExpr parameters:",
  paste0(
    "  min.count = ",
    filterByExpr_parameters$min.count
  ),
  paste0(
    "  min.total.count = ",
    filterByExpr_parameters$min.total.count
  ),
  paste0(
    "  large.n = ",
    filterByExpr_parameters$large.n
  ),
  paste0(
    "  min.prop = ",
    filterByExpr_parameters$min.prop
  ),
  "",
  "IMPORTANT:",
  "  filterByExpr was run exactly ONCE using the full interaction design.",
  paste0(
    "  The same ",
    format_integer_for_readme(
      tested_genes_n
    ),
    " genes were tested in every one of the seven comparisons."
  ),
  "  No separate filtering was performed for individual contrasts.",
  "",
  "4. TESTED COMPARISONS",
  "---------------------",
  "",
  "1) Interaction",
  "   (ASD - Neurotypical)_Female - (ASD - Neurotypical)_Male",
  "   Tests whether the ASD effect differs between Female and Male.",
  "",
  "2) Overall sex",
  "   Female vs Male, averaged equally across Neurotypical and ASD.",
  "",
  "3) Overall group",
  "   ASD vs Neurotypical, averaged equally across Male and Female.",
  "",
  "4) Group effect in Male",
  "   ASD vs Neurotypical among Male samples.",
  "",
  "5) Group effect in Female",
  "   ASD vs Neurotypical among Female samples.",
  "",
  "6) Sex effect in Neurotypical",
  "   Female vs Male among Neurotypical samples.",
  "",
  "7) Sex effect in ASD",
  "   Female vs Male among ASD samples.",
  "",
  "5. GENES TESTED AND SIGNIFICANT IN EACH COMPARISON",
  "---------------------------------------------------",
  comparison_readme_header,
  comparison_readme_rows,
  "",
  "6. XLSX WORKBOOK",
  "----------------",
  paste0(
    "File: ",
    basename(
      output_files$results_xlsx
    )
  ),
  "",
  "Sheets:",
  "  1. Interaction",
  "  2. Overall_Sex",
  "  3. Overall_Group",
  "  4. Group_in_Male",
  "  5. Group_in_Female",
  "  6. Sex_within_Group",
  "",
  "The Sex_within_Group sheet contains two complete tables side by side:",
  "  Female vs Male in Neurotypical",
  "  Female vs Male in ASD",
  "",
  "7. EXPRESSION MATRICES",
  "----------------------",
  paste0(
    "Raw counts: ",
    basename(
      output_files$raw_counts
    )
  ),
  "  Integer pseudobulk counts for genes retained by filterByExpr.",
  "",
  paste0(
    "Visualization matrix: ",
    basename(
      output_files$log2_qn
    )
  ),
  "  raw counts -> log2(count + 1) -> quantile normalization",
  "  This transformed matrix was not used for edgeR testing.",
  "",
  "8. COMPLETE RESULT COLUMNS",
  "--------------------------",
  "Interaction result:",
  "  ensembl_gene_id",
  "  gene",
  "  logFC",
  "  logCPM",
  "  F",
  "  PValue",
  "  FDR",
  "",
  "The Interaction result intentionally has no regulation column.",
  "Its logFC represents a difference of differences rather than a direct",
  "increase or decrease of expression in one biological group.",
  "",
  "All six non-interaction results additionally contain:",
  "  regulation",
  "",
  "For non-interaction results:",
  "  regulation = up   when logFC >= 0",
  "  regulation = down when logFC < 0",
  "",
  "The biological interpretation of positive and negative logFC is recorded in:",
  paste0(
    "  ",
    basename(
      output_files$test_definitions
    )
  ),
  "",
  "9. FILTERING FILES",
  "------------------",
  paste0(
    "Simple gene/filter-stage table: ",
    basename(
      output_files$filtering_stage
    )
  ),
  paste0(
    "Complete gene filtering status: ",
    basename(
      output_files$all_gene_filter_status
    )
  ),
  paste0(
    "Protein-coding mapping status: ",
    basename(
      output_files$mapping_status
    )
  ),
  paste0(
    "filterByExpr status: ",
    basename(
      output_files$filterByExpr_status
    )
  ),
  paste0(
    "Filtering summary: ",
    basename(
      output_files$filtering_summary
    )
  ),
  "",
  "End of README"
)

writeLines(
  text = readme_lines,
  con = output_files$readme,
  useBytes = TRUE
)


# ==============================================================================
# 31. Validate all outputs
# ==============================================================================

expected_interaction_columns <- c(
  "ensembl_gene_id",
  "gene",
  "logFC",
  "logCPM",
  "F",
  "PValue",
  "FDR"
)

expected_standard_result_columns <- c(
  "ensembl_gene_id",
  "gene",
  "logFC",
  "logCPM",
  "F",
  "PValue",
  "FDR",
  "regulation"
)

if (!identical(
  colnames(
    sexInteraction_fullResults_interaction
  ),
  expected_interaction_columns
)) {
  stop(
    "Unexpected columns in the Interaction result: ",
    paste(
      colnames(
        sexInteraction_fullResults_interaction
      ),
      collapse = ", "
    )
  )
}

for (
  standard_result_name in setdiff(
    names(sexInteraction_results),
    "Interaction"
  )
) {
  if (!identical(
    colnames(
      sexInteraction_results[[standard_result_name]]
    ),
    expected_standard_result_columns
  )) {
    stop(
      "Unexpected columns in result: ",
      standard_result_name
    )
  }
}

if (!identical(
  openxlsx::getSheetNames(
    output_files$results_xlsx
  ),
  workbook_sheet_names
)) {
  stop(
    "The XLSX workbook does not contain ",
    "the expected six sheets."
  )
}

all_expected_files <- c(
  unname(result_tsv_files),
  unlist(output_files)
)

missing_files <- all_expected_files[
  !file.exists(all_expected_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Missing output files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

empty_files <- all_expected_files[
  file.info(all_expected_files)$size == 0
]

if (length(empty_files) > 0L) {
  stop(
    "Empty output files:\n",
    paste(
      empty_files,
      collapse = "\n"
    )
  )
}

if (!identical(
  colnames(
    sexInteraction_rawCounts_output
  ),
  colnames(
    sexInteraction_log2RawCounts_QN_output
  )
)) {
  stop(
    "Raw and normalized matrix columns do not match."
  )
}

if (!identical(
  sexInteraction_rawCounts_output$ensembl_gene_id,
  sexInteraction_log2RawCounts_QN_output$ensembl_gene_id
)) {
  stop(
    "Raw and normalized matrix gene orders do not match."
  )
}


# ==============================================================================
# 32. Final report
# ==============================================================================

message(
  "\n",
  paste(rep("=", 80), collapse = "")
)

message(
  "SEX-INTERACTION edgeR ANALYSIS COMPLETED SUCCESSFULLY"
)

message(
  paste(rep("=", 80), collapse = "")
)

message(
  "\nModel: expression ~ fmt_donor_group * sex"
)

message(
  "Samples used: ",
  nrow(sexInteraction_metadata)
)

message(
  "Genes tested: ",
  tested_genes_n
)

message(
  "edgeR result tables: ",
  length(sexInteraction_results)
)

message(
  "XLSX sheets: ",
  length(workbook_sheet_names)
)

message("\nXLSX sheet order:")
message(
  paste(
    seq_along(workbook_sheet_names),
    workbook_sheet_names,
    sep = ". ",
    collapse = "\n"
  )
)

message("\nAnalysis summary:")
print(
  edgeR_analysis_summary,
  n = Inf,
  width = Inf
)

message("\nSeven result TSV files:")
for (
  test_name in names(
    result_tsv_files
  )
) {
  message(
    "  ",
    test_name,
    ": ",
    result_tsv_files[[test_name]]
  )
}

message("\nAdditional files:")
for (
  output_name in names(
    output_files
  )
) {
  message(
    "  ",
    output_name,
    ": ",
    output_files[[output_name]]
  )
}

message(
  "\nHuman-readable README: ",
  normalizePath(
    output_files$readme,
    mustWork = TRUE
  )
)

message(
  paste0(
    "filterByExpr was performed once; all seven contrasts use ",
    tested_genes_n,
    " tested genes."
  )
)

message(
  paste0(
    "\nInteraction result has no regulation column. ",
    "Interpret the sign using sexInteraction_testDefinitions.tsv."
  )
)

message("\nMain result objects in R:")
message(
  "  sexInteraction_fullResults_interaction"
)
message(
  "  sexInteraction_fullResults_overallSex"
)
message(
  "  sexInteraction_fullResults_overallGroup"
)
message(
  "  sexInteraction_fullResults_groupMale"
)
message(
  "  sexInteraction_fullResults_groupFemale"
)
message(
  "  sexInteraction_fullResults_sexNeurotypical"
)
message(
  "  sexInteraction_fullResults_sexASD"
)

# ==============================================================================
# End
# ==============================================================================
sexInteraction_fullResults_interaction %>% 
  filter(FDR < 0.1) %>% dim

sexInteraction_fullResults_groupFemale %>% 
  filter(FDR < 0.1) %>% 
  filter(abs(logFC) > 1) %>% dim

sexInteraction_fullResults_groupMale %>% 
  filter(PValue < 0.05) %>% 
  filter(abs(logFC) > 1)


sexInteraction_fullResults_overallSex  %>% 
  filter(FDR< 0.1) %>% 
  filter(abs(logFC) > 0.5) %>% dim

sexInteraction_fullResults_overallGroup  %>% 
  filter(FDR< 0.1) %>% 
  filter(abs(logFC) > 0.5) %>% dim

sexInteraction_fullResults_groupMale %>% 
  head
# ==============================================================================
# check data
# ==============================================================================
fmt <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

initial_n <- nrow(gene_counts_per_sample_raw_in_tissue_allGenes)
protein_coding_n <- nrow(gene_counts_per_sample_raw_in_tissue_proteinCoding)
low_abundance_removed_n <- sum(!sexInteraction_keep_genes)
tested_n <- nrow(sexInteraction_dge$counts)

interaction_n <- sum(
  sexInteraction_fullResults_interaction$FDR < 0.1,
  na.rm = TRUE
)

group_sig <- sexInteraction_fullResults_overallGroup |>
  dplyr::filter(
    FDR < 0.1,
    abs(logFC) > 0.5
  )

sex_sig <- sexInteraction_fullResults_overallSex |>
  dplyr::filter(
    FDR < 0.1,
    abs(logFC) > 0.5
  )

cat(
  sprintf(
    paste0(
      "The initial pseudobulk count matrix contained %s transcripts, ",
      "of which %s were annotated as protein-coding genes. ",
      "After filtering, %s low-abundance genes were removed, ",
      "leaving %s genes for differential expression analysis across the 16 samples.\n",
      "At FDR < 0.1, %s genes showed a significant sex-by-donor-group interaction. ",
      "For the main effects, donor group: %s genes (ASD, UP: %s & DOWN: %s) ",
      "and sex: %s genes (Female, UP: %s & DOWN: %s).\n"
    ),
    fmt(initial_n),
    fmt(protein_coding_n),
    fmt(low_abundance_removed_n),
    fmt(tested_n),
    fmt(interaction_n),
    fmt(nrow(group_sig)),
    fmt(sum(group_sig$logFC >= 0)),
    fmt(sum(group_sig$logFC < 0)),
    fmt(nrow(sex_sig)),
    fmt(sum(sex_sig$logFC >= 0)),
    fmt(sum(sex_sig$logFC < 0))
  )
)
