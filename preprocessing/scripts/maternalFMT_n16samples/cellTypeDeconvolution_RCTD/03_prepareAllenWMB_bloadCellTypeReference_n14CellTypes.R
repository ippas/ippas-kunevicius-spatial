#!/usr/bin/env Rscript

message(
  "VERSION: 02 prepare Allen WMB Yao2023 broad cell-type reference | ",
  "n14 cell types | max10000 cells/type | seed7 | manual H5AD dimnames"
)

# ==============================================================================
# 1. Required packages
# ==============================================================================

required_packages <- c(
  "data.table",
  "Matrix",
  "HDF5Array",
  "DelayedArray",
  "rhdf5"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {

  stop(
    "Missing required package(s):\n",
    paste(missing_packages, collapse = "\n"),
    "\n\nInstall with:\n\n",
    "install.packages('data.table')\n\n",
    "if (!requireNamespace('BiocManager', quietly = TRUE)) ",
    "install.packages('BiocManager')\n",
    "BiocManager::install(c(",
    "'HDF5Array', 'DelayedArray', 'rhdf5'",
    "))\n",
    call. = FALSE
  )
}


# ==============================================================================
# 2. Configuration
# ==============================================================================

project_root <-
  "/home/mateusz/projects/ippas-kunevicius-spatial"

reference_root <- file.path(
  project_root,
  "data",
  "reference",
  "allenWholeMouseBrain_Yao2023"
)

cell_metadata_file <- file.path(
  reference_root,
  "metadata",
  "WMB-10X",
  "20241115",
  "views",
  "cell_metadata_with_cluster_annotation.csv"
)

gene_metadata_file <- file.path(
  reference_root,
  "metadata",
  "WMB-10X",
  "20241115",
  "gene.csv"
)

taxonomy_file <- file.path(
  reference_root,
  "metadata",
  "WMB-taxonomy",
  "20231215",
  "views",
  "cluster_to_cluster_annotation_membership_pivoted.csv"
)

expression_dir <- file.path(
  reference_root,
  "expression_matrices",
  "WMB-10Xv3",
  "20230630"
)

output_dir <- file.path(
  reference_root,
  "processed",
  "broadCellTypeReference_n14CellTypes"
)

max_cells_per_type <- 10000L

min_cells_per_type <- 25L

random_seed <- 7L


# ==============================================================================
# 3. Broad cell-type definition
# ==============================================================================

broad_cell_type_levels <- c(
  "Excitatory_neurons",
  "Inhibitory_neurons",
  "Immature_neurons",
  "Other_neurons",
  "Astrocytes",
  "Ependymal_ventricular",
  "OPC",
  "Oligodendrocytes",
  "Endothelial",
  "Pericytes",
  "Other_vascular",
  "Microglia",
  "BAM",
  "Other_immune"
)


# ==============================================================================
# 4. Validate inputs
# ==============================================================================

required_files <- c(
  cell_metadata_file,
  gene_metadata_file,
  taxonomy_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {

  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

if (!dir.exists(expression_dir)) {

  stop(
    "Expression directory does not exist:\n",
    expression_dir,
    call. = FALSE
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(output_dir)) {

  stop(
    "Could not create output directory:\n",
    output_dir,
    call. = FALSE
  )
}

message("\nReference root:")
message(reference_root)

message("\nOutput directory:")
message(output_dir)


# ==============================================================================
# 5. Helper: extract numeric Allen annotation prefix
# ==============================================================================

extract_numeric_prefix <- function(x) {

  x <- trimws(x)

  result <- rep(
    NA_integer_,
    length(x)
  )

  has_prefix <- grepl(
    "^[0-9]+",
    x
  )

  result[has_prefix] <- as.integer(
    sub(
      "^([0-9]+).*",
      "\\1",
      x[has_prefix]
    )
  )

  result
}


# ==============================================================================
# 6. Helper: read obs/var index directly from H5AD
#
# Allen H5AD:
#
# obs index = cell_label
# var index = gene_identifier
#
# HDF5Array can read the matrix but does not always expose these as dimnames.
# ==============================================================================

read_h5ad_index <- function(
    filepath,
    group,
    fallback_names,
    expected_length
) {

  group_path <- paste0(
    "/",
    group
  )

  attributes <- tryCatch(
    rhdf5::h5readAttributes(
      filepath,
      group_path
    ),
    error = function(e) {
      list()
    }
  )

  index_attribute <- attributes[["_index"]]

  candidates <- unique(
    c(
      as.character(index_attribute),
      fallback_names,
      "_index",
      "index"
    )
  )

  candidates <- candidates[
    !is.na(candidates) &
      nzchar(candidates)
  ]

  for (candidate in candidates) {

    dataset_path <- paste0(
      "/",
      group,
      "/",
      candidate
    )

    values <- tryCatch(
      rhdf5::h5read(
        filepath,
        dataset_path
      ),
      error = function(e) {
        NULL
      }
    )

    if (is.null(values)) {
      next
    }

    values <- as.character(values)

    if (length(values) == expected_length) {

      message(
        "Read ",
        group,
        " index from: ",
        dataset_path
      )

      return(values)
    }
  }

  # --------------------------------------------------------------------------
  # Diagnostic fallback
  # --------------------------------------------------------------------------

  message(
    "\nCould not automatically identify ",
    group,
    " index."
  )

  message(
    "Available HDF5 objects under /",
    group,
    ":"
  )

  h5_structure <- rhdf5::h5ls(
    filepath,
    recursive = TRUE
  )

  print(
    h5_structure[
      grepl(
        paste0("^/", group),
        h5_structure[["group"]]
      ) |
        h5_structure[["group"]] == paste0("/", group),
      ,
      drop = FALSE
    ]
  )

  stop(
    "Could not find ",
    group,
    " index of expected length ",
    expected_length,
    " in H5AD file:\n",
    filepath,
    call. = FALSE
  )
}


# ==============================================================================
# 7. Read cell metadata
# ==============================================================================

message("\n============================================================")
message("Reading Allen cell metadata")
message("============================================================")

metadata_columns <- names(
  data.table::fread(
    cell_metadata_file,
    nrows = 0L,
    showProgress = FALSE
  )
)

required_metadata_columns <- c(
  "cell_label",
  "feature_matrix_label",
  "class",
  "subclass",
  "neurotransmitter"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  metadata_columns
)

if (length(missing_metadata_columns) > 0L) {

  stop(
    "Missing cell metadata column(s):\n",
    paste(
      missing_metadata_columns,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

cell_metadata <- data.table::fread(
  cell_metadata_file,
  select = required_metadata_columns,
  showProgress = TRUE
)

message(
  "All cells in metadata: ",
  format(
    nrow(cell_metadata),
    big.mark = ","
  )
)


# ==============================================================================
# 8. Keep WMB-10Xv3 cells
# ==============================================================================

cell_metadata <- cell_metadata[
  startsWith(
    feature_matrix_label,
    "WMB-10Xv3-"
  )
]

message(
  "WMB-10Xv3 cells: ",
  format(
    nrow(cell_metadata),
    big.mark = ","
  )
)

if (nrow(cell_metadata) == 0L) {

  stop(
    "No WMB-10Xv3 cells found.",
    call. = FALSE
  )
}

if (anyDuplicated(
  cell_metadata[["cell_label"]]
) > 0L) {

  stop(
    "Duplicated cell_label values detected.",
    call. = FALSE
  )
}


# ==============================================================================
# 9. Extract Allen class/subclass numeric codes
# ==============================================================================

cell_metadata[
  ,
  class_code := extract_numeric_prefix(
    class
  )
]

cell_metadata[
  ,
  subclass_code := extract_numeric_prefix(
    subclass
  )
]


# ==============================================================================
# 10. Assign 14 broad cell types
# ==============================================================================

message("\n============================================================")
message("Assigning 14 broad cell types")
message("============================================================")

cell_metadata[
  ,
  broad_cell_type := NA_character_
]


# ------------------------------------------------------------------------------
# Astrocytes
# ------------------------------------------------------------------------------

cell_metadata[
  subclass_code %in% 316:320,
  broad_cell_type := "Astrocytes"
]


# ------------------------------------------------------------------------------
# Ependymal / ventricular
# ------------------------------------------------------------------------------

cell_metadata[
  subclass_code %in% 321:325,
  broad_cell_type := "Ependymal_ventricular"
]


# ------------------------------------------------------------------------------
# OPC
# ------------------------------------------------------------------------------

cell_metadata[
  subclass_code == 326L,
  broad_cell_type := "OPC"
]


# ------------------------------------------------------------------------------
# Oligodendrocytes
# ------------------------------------------------------------------------------

cell_metadata[
  subclass_code == 327L,
  broad_cell_type := "Oligodendrocytes"
]


# ------------------------------------------------------------------------------
# Vascular
# ------------------------------------------------------------------------------

cell_metadata[
  subclass_code == 331L,
  broad_cell_type := "Pericytes"
]

cell_metadata[
  subclass_code == 333L,
  broad_cell_type := "Endothelial"
]

cell_metadata[
  subclass_code %in% c(
    329L,
    330L,
    332L
  ),
  broad_cell_type := "Other_vascular"
]


# ------------------------------------------------------------------------------
# Immune
# ------------------------------------------------------------------------------

cell_metadata[
  subclass_code == 334L,
  broad_cell_type := "Microglia"
]

cell_metadata[
  subclass_code == 335L,
  broad_cell_type := "BAM"
]

cell_metadata[
  subclass_code %in% c(
    336L,
    337L,
    338L
  ),
  broad_cell_type := "Other_immune"
]


# ------------------------------------------------------------------------------
# Immature neurons
# ------------------------------------------------------------------------------

cell_metadata[
  subclass_code %in% c(
    38L,
    45L
  ),
  broad_cell_type := "Immature_neurons"
]


# ------------------------------------------------------------------------------
# Mature neuronal classes = Allen classes 01-29
# ------------------------------------------------------------------------------

neuronal_cell <-
  !is.na(cell_metadata[["class_code"]]) &
  cell_metadata[["class_code"]] >= 1L &
  cell_metadata[["class_code"]] <= 29L


# ------------------------------------------------------------------------------
# Excitatory = Glut
# ------------------------------------------------------------------------------

cell_metadata[
  is.na(broad_cell_type) &
    neuronal_cell &
    neurotransmitter == "Glut",
  broad_cell_type := "Excitatory_neurons"
]


# ------------------------------------------------------------------------------
# Inhibitory = GABA
# ------------------------------------------------------------------------------

cell_metadata[
  is.na(broad_cell_type) &
    neuronal_cell &
    neurotransmitter == "GABA",
  broad_cell_type := "Inhibitory_neurons"
]


# ------------------------------------------------------------------------------
# Remaining neuronal neurotransmitter classes
# ------------------------------------------------------------------------------

cell_metadata[
  is.na(broad_cell_type) &
    neuronal_cell,
  broad_cell_type := "Other_neurons"
]


# ==============================================================================
# 11. Inspect unassigned cells
# ==============================================================================

unassigned_cells <- cell_metadata[
  is.na(broad_cell_type)
]

unassigned_summary <- unassigned_cells[
  ,
  .N,
  by = .(
    class,
    subclass,
    neurotransmitter
  )
][
  order(-N)
]

unassigned_summary_file <- file.path(
  output_dir,
  "allenWMB_unassignedCells.tsv"
)

data.table::fwrite(
  unassigned_summary,
  unassigned_summary_file,
  sep = "\t"
)

# Only class 32 OEC is intentionally excluded.

unexpected_unassigned <- unassigned_cells[
  is.na(class_code) |
    class_code != 32L
]

if (nrow(unexpected_unassigned) > 0L) {

  message("\nUnexpected unassigned cell types:")

  print(
    unexpected_unassigned[
      ,
      .N,
      by = .(
        class,
        subclass,
        neurotransmitter
      )
    ][
      order(-N)
    ]
  )

  stop(
    "Unexpected Allen annotations remain unassigned.\n",
    "Inspect:\n",
    unassigned_summary_file,
    call. = FALSE
  )
}

message(
  "Intentionally excluded OEC cells: ",
  format(
    nrow(unassigned_cells),
    big.mark = ","
  )
)


# ==============================================================================
# 12. Save taxonomy mapping
# ==============================================================================

mapping_table <- unique(
  cell_metadata[
    !is.na(broad_cell_type),
    .(
      class,
      subclass,
      neurotransmitter,
      broad_cell_type
    )
  ]
)

data.table::setorder(
  mapping_table,
  broad_cell_type,
  class,
  subclass,
  neurotransmitter
)

mapping_file <- file.path(
  output_dir,
  "allenWMB_toBroadCellType_n14_mapping.tsv"
)

data.table::fwrite(
  mapping_table,
  mapping_file,
  sep = "\t"
)


# ==============================================================================
# 13. Cell counts before sampling
# ==============================================================================

cell_counts_before <- cell_metadata[
  !is.na(broad_cell_type),
  .N,
  by = broad_cell_type
]

cell_counts_before[
  ,
  broad_cell_type := factor(
    broad_cell_type,
    levels = broad_cell_type_levels
  )
]

data.table::setorder(
  cell_counts_before,
  broad_cell_type
)

message("\nCell counts before sampling:")

print(
  cell_counts_before
)

missing_broad_types <- setdiff(
  broad_cell_type_levels,
  as.character(
    cell_counts_before[["broad_cell_type"]]
  )
)

if (length(missing_broad_types) > 0L) {

  stop(
    "Missing broad cell type(s):\n",
    paste(
      missing_broad_types,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

too_small_types <- cell_counts_before[
  N < min_cells_per_type
]

if (nrow(too_small_types) > 0L) {

  stop(
    "Broad cell type(s) have fewer than ",
    min_cells_per_type,
    " cells:\n",
    paste(
      paste0(
        too_small_types[["broad_cell_type"]],
        ": ",
        too_small_types[["N"]]
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

data.table::fwrite(
  cell_counts_before,
  file.path(
    output_dir,
    "cellCounts_beforeSampling.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 14. Downsample max 10,000 cells per broad type
# ==============================================================================

message(
  "\nSampling maximum ",
  format(
    max_cells_per_type,
    big.mark = ","
  ),
  " cells per broad cell type..."
)

set.seed(
  random_seed
)

mapped_cells <- cell_metadata[
  !is.na(broad_cell_type)
]

selected_cells <- mapped_cells[
  ,
  {

    n_select <- min(
      .N,
      max_cells_per_type
    )

    .SD[
      sample.int(
        .N,
        n_select,
        replace = FALSE
      )
    ]

  },
  by = broad_cell_type
]

selected_cells[
  ,
  broad_cell_type := factor(
    broad_cell_type,
    levels = broad_cell_type_levels
  )
]

data.table::setorder(
  selected_cells,
  broad_cell_type,
  feature_matrix_label,
  cell_label
)

if (anyDuplicated(
  selected_cells[["cell_label"]]
) > 0L) {

  stop(
    "Duplicated selected cell labels.",
    call. = FALSE
  )
}


# ==============================================================================
# 15. Cell counts after sampling
# ==============================================================================

cell_counts_after <- selected_cells[
  ,
  .N,
  by = broad_cell_type
]

data.table::setorder(
  cell_counts_after,
  broad_cell_type
)

message("\nCell counts after sampling:")

print(
  cell_counts_after
)

data.table::fwrite(
  cell_counts_after,
  file.path(
    output_dir,
    "cellCounts_afterSampling.tsv"
  ),
  sep = "\t"
)

data.table::fwrite(
  selected_cells,
  file.path(
    output_dir,
    "selectedCells_n14_max10000PerType_seed7.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 16. Selected cells by regional feature matrix
# ==============================================================================

selected_by_matrix <- selected_cells[
  ,
  .N,
  by = .(
    broad_cell_type,
    feature_matrix_label
  )
]

data.table::setorder(
  selected_by_matrix,
  broad_cell_type,
  feature_matrix_label
)

data.table::fwrite(
  selected_by_matrix,
  file.path(
    output_dir,
    "selectedCellCounts_byFeatureMatrix_seed7.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 17. Find raw WMB-10Xv3 H5AD files
# ==============================================================================

raw_h5ad_files <- list.files(
  expression_dir,
  pattern = "-raw\\.h5ad$",
  full.names = TRUE
)

if (length(raw_h5ad_files) == 0L) {

  stop(
    "No raw H5AD files found in:\n",
    expression_dir,
    call. = FALSE
  )
}

file_table <- data.table::data.table(
  feature_matrix_label = sub(
    "-raw\\.h5ad$",
    "",
    basename(raw_h5ad_files)
  ),
  h5ad_file = raw_h5ad_files
)

data.table::setorder(
  file_table,
  feature_matrix_label
)

required_feature_matrices <- sort(
  unique(
    as.character(
      selected_cells[["feature_matrix_label"]]
    )
  )
)

missing_feature_matrices <- setdiff(
  required_feature_matrices,
  file_table[["feature_matrix_label"]]
)

if (length(missing_feature_matrices) > 0L) {

  stop(
    "Missing raw H5AD matrix/matrices:\n",
    paste(
      missing_feature_matrices,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

message(
  "\nRaw H5AD files required: ",
  length(required_feature_matrices)
)


# ==============================================================================
# 18. Extract selected cells from each H5AD
# ==============================================================================

message("\n============================================================")
message("Extracting selected raw-count cells from Allen H5AD files")
message("============================================================")

reference_matrices <- vector(
  "list",
  length(required_feature_matrices)
)

names(reference_matrices) <-
  required_feature_matrices

reference_gene_ids <- NULL


for (i in seq_along(required_feature_matrices)) {

  matrix_label <-
    required_feature_matrices[i]

  h5ad_file <- file_table[
    feature_matrix_label == matrix_label,
    h5ad_file
  ][1]

  wanted_cells <- as.character(
    selected_cells[
      feature_matrix_label == matrix_label,
      cell_label
    ]
  )

  message("\n============================================================")

  message(
    "[",
    i,
    "/",
    length(required_feature_matrices),
    "] ",
    matrix_label
  )

  message(
    "File: ",
    basename(h5ad_file)
  )

  message(
    "Selected cells: ",
    format(
      length(wanted_cells),
      big.mark = ","
    )
  )


  # ---------------------------------------------------------------------------
  # Open central H5AD matrix on disk.
  #
  # The Allen file contains X in AnnData format.
  # HDF5Array presents it as genes x cells.
  # ---------------------------------------------------------------------------

  message(
    "Opening H5AD matrix on disk..."
  )

  h5_matrix <- suppressWarnings(
    HDF5Array::H5ADMatrix(
      h5ad_file
    )
  )

  message(
    "H5AD dimensions: ",
    format(
      nrow(h5_matrix),
      big.mark = ","
    ),
    " genes x ",
    format(
      ncol(h5_matrix),
      big.mark = ","
    ),
    " cells"
  )


  # ---------------------------------------------------------------------------
  # IMPORTANT:
  #
  # HDF5Array does not expose the Allen obs / var indices as dimnames.
  # Read them explicitly.
  # ---------------------------------------------------------------------------

  cell_labels <- read_h5ad_index(
    filepath = h5ad_file,
    group = "obs",
    fallback_names = c(
      "cell_label"
    ),
    expected_length = ncol(
      h5_matrix
    )
  )

  gene_ids <- read_h5ad_index(
    filepath = h5ad_file,
    group = "var",
    fallback_names = c(
      "gene_identifier"
    ),
    expected_length = nrow(
      h5_matrix
    )
  )

  message(
    "Cell labels read: ",
    format(
      length(cell_labels),
      big.mark = ","
    )
  )

  message(
    "Gene identifiers read: ",
    format(
      length(gene_ids),
      big.mark = ","
    )
  )


  # ---------------------------------------------------------------------------
  # Match selected Allen cells to H5AD column indices
  # ---------------------------------------------------------------------------

  cell_index <- match(
    wanted_cells,
    cell_labels
  )

  if (anyNA(cell_index)) {

    missing_cells <- wanted_cells[
      is.na(cell_index)
    ]

    stop(
      "Selected Allen cells were not found in H5AD matrix ",
      matrix_label,
      ":\n",
      paste(
        head(
          missing_cells,
          100L
        ),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Ensure identical gene universe/order across regional matrices
  # ---------------------------------------------------------------------------

  current_gene_ids <- gene_ids

  if (is.null(reference_gene_ids)) {

    reference_gene_ids <-
      current_gene_ids

  } else if (!identical(
    reference_gene_ids,
    current_gene_ids
  )) {

    stop(
      "Gene identifiers/order differs between Allen H5AD files.\n",
      "Problematic matrix: ",
      matrix_label,
      call. = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Extract ONLY selected cells
  # ---------------------------------------------------------------------------

  message(
    "Reading only selected cells from disk..."
  )

  selected_delayed <- h5_matrix[
    ,
    cell_index,
    drop = FALSE
  ]

  message(
    "Converting selected counts to sparse dgCMatrix..."
  )

  selected_sparse <- methods::as(
    selected_delayed,
    "dgCMatrix"
  )

  rownames(selected_sparse) <-
    reference_gene_ids

  colnames(selected_sparse) <-
    wanted_cells

  reference_matrices[[matrix_label]] <- selected_sparse

  message(
    "Matrix extracted: ",
    format(
      nrow(selected_sparse),
      big.mark = ","
    ),
    " genes x ",
    format(
      ncol(selected_sparse),
      big.mark = ","
    ),
    " cells"
  )

  message(
    "Non-zero counts: ",
    format(
      length(
        selected_sparse@x
      ),
      big.mark = ","
    )
  )

  rm(
    h5_matrix,
    selected_delayed,
    selected_sparse,
    cell_labels,
    gene_ids,
    cell_index
  )

  invisible(
    gc()
  )
}


# ==============================================================================
# 19. Combine regional matrices
# ==============================================================================

message("\n============================================================")
message("Combining regional Allen matrices")
message("============================================================")

reference_counts <- do.call(
  cbind,
  reference_matrices
)

rm(
  reference_matrices
)

invisible(
  gc()
)

if (anyDuplicated(
  colnames(reference_counts)
) > 0L) {

  stop(
    "Duplicated cell labels in combined reference.",
    call. = FALSE
  )
}

message(
  "Combined dimensions: ",
  format(
    nrow(reference_counts),
    big.mark = ","
  ),
  " genes x ",
  format(
    ncol(reference_counts),
    big.mark = ","
  ),
  " cells"
)


# ==============================================================================
# 20. Read Allen gene metadata
# ==============================================================================

message("\nReading Allen gene metadata...")

gene_metadata <- data.table::fread(
  gene_metadata_file,
  showProgress = FALSE
)

required_gene_columns <- c(
  "gene_identifier",
  "gene_symbol"
)

missing_gene_columns <- setdiff(
  required_gene_columns,
  names(gene_metadata)
)

if (length(missing_gene_columns) > 0L) {

  stop(
    "Missing gene metadata column(s):\n",
    paste(
      missing_gene_columns,
      collapse = "\n"
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 21. Convert Ensembl IDs to mouse gene symbols
# ==============================================================================

gene_match <- match(
  rownames(reference_counts),
  gene_metadata[["gene_identifier"]]
)

gene_symbols <- gene_metadata[["gene_symbol"]][gene_match]

valid_gene <-
  !is.na(gene_symbols) &
  nzchar(gene_symbols)

duplicated_gene_symbol <- duplicated(
  gene_symbols
)

keep_gene <-
  valid_gene &
  !duplicated_gene_symbol

message(
  "Genes before symbol filtering: ",
  format(
    nrow(reference_counts),
    big.mark = ","
  )
)

message(
  "Genes without valid symbol: ",
  sum(
    !valid_gene
  )
)

message(
  "Duplicated gene symbols removed: ",
  sum(
    duplicated_gene_symbol &
      valid_gene
  )
)

reference_counts <- reference_counts[
  keep_gene,
  ,
  drop = FALSE
]

rownames(reference_counts) <-
  gene_symbols[
    keep_gene
  ]

message(
  "Genes retained: ",
  format(
    nrow(reference_counts),
    big.mark = ","
  )
)


# ==============================================================================
# 22. Match selected metadata to count matrix
# ==============================================================================

selected_cells[
  ,
  broad_cell_type := as.character(
    broad_cell_type
  )
]

reference_cell_metadata <- selected_cells[
  match(
    colnames(reference_counts),
    cell_label
  )
]

if (anyNA(
  reference_cell_metadata[["cell_label"]]
)) {

  stop(
    "Failed to match reference columns to cell metadata.",
    call. = FALSE
  )
}

if (!identical(
  reference_cell_metadata[["cell_label"]],
  colnames(reference_counts)
)) {

  stop(
    "Reference metadata and count matrix cell order differ.",
    call. = FALSE
  )
}


# ==============================================================================
# 23. Construct RCTD cell-type vector
# ==============================================================================

cell_types <- factor(
  reference_cell_metadata[["broad_cell_type"]],
  levels = broad_cell_type_levels
)

names(cell_types) <-
  reference_cell_metadata[["cell_label"]]


# ==============================================================================
# 24. Calculate reference nUMI
# ==============================================================================

nUMI <- Matrix::colSums(
  reference_counts
)

names(nUMI) <-
  colnames(reference_counts)

if (any(
  nUMI <= 0
)) {

  stop(
    "Reference contains cells with zero total UMI counts.",
    call. = FALSE
  )
}


# ==============================================================================
# 25. Final QC
# ==============================================================================

final_cell_counts <- sort(
  table(
    cell_types
  ),
  decreasing = TRUE
)

message("\n============================================================")
message("FINAL REFERENCE")
message("============================================================")

message(
  "Genes: ",
  format(
    nrow(reference_counts),
    big.mark = ","
  )
)

message(
  "Cells: ",
  format(
    ncol(reference_counts),
    big.mark = ","
  )
)

message(
  "Broad cell types: ",
  nlevels(cell_types)
)

print(
  final_cell_counts
)

if (nlevels(cell_types) != 14L) {

  stop(
    "Reference does not contain exactly 14 broad cell types.",
    call. = FALSE
  )
}

if (any(
  final_cell_counts <
    min_cells_per_type
)) {

  stop(
    "At least one broad cell type contains fewer than ",
    min_cells_per_type,
    " cells.",
    call. = FALSE
  )
}


# ==============================================================================
# 26. Construct final reference object
# ==============================================================================

reference_object <- list(

  counts = reference_counts,

  cell_types = cell_types,

  nUMI = nUMI,

  cell_metadata = reference_cell_metadata,

  broad_cell_type_levels =
    broad_cell_type_levels,

  parameters = list(

    source =
      "Allen Whole Mouse Brain WMB-10Xv3 / Yao et al. 2023",

    taxonomy =
      "WMB-taxonomy 20231215",

    metadata =
      "WMB-10X 20241115",

    expression =
      "WMB-10Xv3 20230630 raw counts",

    n_broad_cell_types =
      14L,

    max_cells_per_type =
      max_cells_per_type,

    min_cells_per_type =
      min_cells_per_type,

    random_seed =
      random_seed,

    excluded_class =
      "32 OEC",

    count_scale =
      "raw counts"
  )
)


# ==============================================================================
# 27. Save final RDS
# ==============================================================================

reference_rds_file <- file.path(
  output_dir,
  paste0(
    "allenWMB_broadCellTypeReference_",
    "n14CellTypes_",
    "max10000PerType_",
    "seed7.rds"
  )
)

message("\nSaving reference RDS...")

saveRDS(
  reference_object,
  reference_rds_file,
  compress = FALSE
)


# ==============================================================================
# 28. Save reference cell metadata
# ==============================================================================

reference_metadata_file <- file.path(
  output_dir,
  "referenceCellMetadata_n14CellTypes_seed7.tsv"
)

data.table::fwrite(
  reference_cell_metadata,
  reference_metadata_file,
  sep = "\t"
)


# ==============================================================================
# 29. Save final cell counts
# ==============================================================================

final_cell_counts_table <-
  data.table::data.table(

    broad_cell_type =
      names(final_cell_counts),

    n_cells =
      as.integer(
        final_cell_counts
      )
  )

final_cell_counts_file <- file.path(
  output_dir,
  "finalReference_cellCounts_n14CellTypes_seed7.tsv"
)

data.table::fwrite(
  final_cell_counts_table,
  final_cell_counts_file,
  sep = "\t"
)


# ==============================================================================
# 30. Save parameters
# ==============================================================================

parameters_table <- data.table::data.table(

  parameter = c(
    "source",
    "taxonomy",
    "metadata",
    "expression",
    "n_broad_cell_types",
    "max_cells_per_type",
    "min_cells_per_type",
    "random_seed",
    "excluded_class",
    "count_scale"
  ),

  value = c(
    "Allen Whole Mouse Brain WMB-10Xv3 / Yao et al. 2023",
    "WMB-taxonomy 20231215",
    "WMB-10X 20241115",
    "WMB-10Xv3 20230630 raw counts",
    "14",
    as.character(
      max_cells_per_type
    ),
    as.character(
      min_cells_per_type
    ),
    as.character(
      random_seed
    ),
    "32 OEC",
    "raw counts"
  )
)

parameters_file <- file.path(
  output_dir,
  "referenceParameters_seed7.tsv"
)

data.table::fwrite(
  parameters_table,
  parameters_file,
  sep = "\t"
)


# ==============================================================================
# 31. Final summary
# ==============================================================================

reference_size_gib <- file.info(
  reference_rds_file
)[["size"]] / 1024^3

message("\n============================================================")
message("REFERENCE PREPARATION FINISHED")
message("============================================================")

message(
  "Reference RDS:"
)

message(
  reference_rds_file
)

message(
  "\nRDS size: ",
  round(
    reference_size_gib,
    2
  ),
  " GiB"
)

message(
  "\nFinal cell counts:"
)

message(
  final_cell_counts_file
)

message(
  "\nMapping:"
)

message(
  mapping_file
)

message(
  "\nParameters:"
)

message(
  parameters_file
)

message("============================================================")