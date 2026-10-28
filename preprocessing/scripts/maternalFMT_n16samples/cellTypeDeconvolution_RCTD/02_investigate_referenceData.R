#!/usr/bin/env Rscript
message("VERSION: 01 inspect Allen Whole Mouse Brain cell-type annotations")
# ==============================================================================
# 1. Paths
# ==============================================================================
project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
reference_dir <- file.path(
  project_root,
  "data",
  "reference",
  "allenWholeMouseBrain_Yao2023"
)
cell_metadata_file <- file.path(
  reference_dir,
  "metadata",
  "WMB-10X",
  "20241115",
  "views",
  "cell_metadata_with_cluster_annotation.csv"
)
taxonomy_file <- file.path(
  reference_dir,
  "metadata",
  "WMB-taxonomy",
  "20231215",
  "views",
  "cluster_to_cluster_annotation_membership_pivoted.csv"
)
gene_file <- file.path(
  reference_dir,
  "metadata",
  "WMB-10X",
  "20241115",
  "gene.csv"
)
required_files <- c(
  cell_metadata_file,
  taxonomy_file,
  gene_file
)
missing_files <- required_files[
  !file.exists(required_files)
]
if (length(missing_files) > 0L) {
  stop(
    "Missing files:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}
# ==============================================================================
# 2. Read taxonomy
# ==============================================================================
message("\nReading taxonomy...")
taxonomy <- read.csv(
  taxonomy_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
message("Taxonomy rows: ", nrow(taxonomy))
message("Taxonomy columns:")
print(colnames(taxonomy))
# ==============================================================================
# 3. Inspect taxonomy annotation levels
# ==============================================================================
message("\n============================================================")
message("TAXONOMY LEVELS")
message("============================================================")
possible_annotation_columns <- grep(
  "class|subclass|supertype|cluster|neurotrans",
  colnames(taxonomy),
  ignore.case = TRUE,
  value = TRUE
)
for (column in possible_annotation_columns) {
  values <- taxonomy[[column]]
  values <- values[
    !is.na(values) &
      values != ""
  ]
  message("\n--- ", column, " ---")
  message("Unique values: ", length(unique(values)))
  if (length(unique(values)) <= 100L) {
    print(sort(unique(values)))
  }
}
# ==============================================================================
# 4. Read cell metadata
# ==============================================================================
message("\nReading cell metadata...")
cell_metadata <- read.csv(
  cell_metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
message("Cells: ", nrow(cell_metadata))
message("Cell metadata columns:")
print(colnames(cell_metadata))
# ==============================================================================
# 5. Find broad annotation columns in cell metadata
# ==============================================================================
message("\n============================================================")
message("CELL METADATA ANNOTATION COLUMNS")
message("============================================================")
annotation_columns <- grep(
  "class|subclass|supertype|cluster|neurotrans",
  colnames(cell_metadata),
  ignore.case = TRUE,
  value = TRUE
)
print(annotation_columns)
# ==============================================================================
# 6. Counts per broad annotation
# ==============================================================================
for (column in annotation_columns) {
  values <- cell_metadata[[column]]
  values <- values[
    !is.na(values) &
      values != ""
  ]
  n_unique <- length(unique(values))
  message("\n============================================================")
  message(column)
  message("Unique values: ", n_unique)
  message("============================================================")
  counts <- sort(
    table(values),
    decreasing = TRUE
  )
  if (n_unique <= 100L) {
    print(counts)
  } else {
    print(head(counts, 50))
  }
}
# ==============================================================================
# 7. Read gene metadata
# ==============================================================================
message("\nReading gene metadata...")
gene_metadata <- read.csv(
  gene_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
message("Genes: ", nrow(gene_metadata))
message("Gene metadata columns:")
print(colnames(gene_metadata))
# ==============================================================================
# 8. Final summary
# ==============================================================================
message("\n============================================================")
message("DONE")
message("============================================================")
message("Cells in reference: ", nrow(cell_metadata))
message("Taxonomy rows: ", nrow(taxonomy))
message("Genes: ", nrow(gene_metadata))
message("============================================================")


taxonomy[
  taxonomy[["class"]] %in% c(
    "04 DG-IMN Glut",
    "05 OB-IMN GABA",
    "30 Astro-Epen",
    "31 OPC-Oligo",
    "33 Vascular",
    "34 Immune"
  ),
  c(
    "class",
    "subclass",
    "neurotransmitter"
  )
] |>
  unique() |>
  print(row.names = FALSE)