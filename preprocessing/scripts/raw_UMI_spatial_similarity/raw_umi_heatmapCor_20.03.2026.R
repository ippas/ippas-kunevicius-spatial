# ==============================================================================
# Raw UMI correlation analysis for standard Visium samples
#
# Input:
# - samples_list: named list of Seurat objects
# - output_dir: directory where results will be saved
#
# Analysis:
# - exclude samples containing "undetermined" in their names
# - extract raw UMI counts from nCount_RNA
# - match samples by spatial barcode
# - calculate Pearson and Spearman correlations
# - save correlation matrices
# - save heatmaps
#
# Heatmap color scale:
# - negative correlations: blue
# - zero correlation: white
# - positive correlations: red
# ==============================================================================


# ==============================================================================
# 1. Packages
# ==============================================================================

suppressPackageStartupMessages({

  library(Seurat)
  library(dplyr)
  library(purrr)
  library(pheatmap)

})


# ==============================================================================
# 2. Check required objects
# ==============================================================================

if (!exists("samples_list")) {

  stop(
    "Object samples_list was not found."
  )

}


if (!is.list(samples_list)) {

  stop(
    "samples_list must be a named list of Seurat objects."
  )

}


if (is.null(names(samples_list))) {

  stop(
    "samples_list must have sample names."
  )

}


if (!exists("output_dir")) {

  stop(
    "Object output_dir was not found."
  )

}


if (!dir.exists(output_dir)) {

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

}


# ==============================================================================
# 3. Select standard samples
# ==============================================================================

standard_sample_names <- names(samples_list)[
  !grepl(
    pattern = "undetermined",
    x = names(samples_list),
    ignore.case = TRUE
  )
]


if (length(standard_sample_names) < 2) {

  stop(
    "Fewer than two standard samples were found."
  )

}


standard_samples_list <- samples_list[
  standard_sample_names
]


message(
  "Number of standard samples: ",
  length(standard_samples_list)
)


message(
  "Standard samples: ",
  paste(
    names(standard_samples_list),
    collapse = ", "
  )
)


# ==============================================================================
# 4. Extract raw UMI values per spatial barcode
# ==============================================================================

umi_vectors <- imap(
  standard_samples_list,
  function(seurat_object, sample_id) {

    if (!inherits(seurat_object, "Seurat")) {

      stop(
        "Object for sample ",
        sample_id,
        " is not a Seurat object."
      )

    }


    metadata <- seurat_object[[]]


    if (!"nCount_RNA" %in% colnames(metadata)) {

      stop(
        "Column nCount_RNA was not found for sample: ",
        sample_id
      )

    }


    barcode_names <- colnames(
      seurat_object
    )


    raw_umi <- as.numeric(
      metadata[["nCount_RNA"]]
    )


    if (length(barcode_names) != length(raw_umi)) {

      stop(
        "Number of barcodes does not match the number of nCount_RNA values for sample: ",
        sample_id
      )

    }


    if (anyDuplicated(barcode_names) > 0) {

      stop(
        "Duplicated spatial barcodes were detected for sample: ",
        sample_id
      )

    }


    if (any(!is.finite(raw_umi))) {

      stop(
        "Non-finite nCount_RNA values were detected for sample: ",
        sample_id
      )

    }


    names(raw_umi) <- barcode_names

    raw_umi

  }
)


# ==============================================================================
# 5. Find spatial barcodes shared by all samples
# ==============================================================================

shared_barcodes <- reduce(
  map(
    umi_vectors,
    names
  ),
  intersect
)


if (length(shared_barcodes) == 0) {

  stop(
    "No spatial barcodes are shared by all standard samples."
  )

}


message(
  "Number of spatial barcodes shared by all samples: ",
  length(shared_barcodes)
)


# ==============================================================================
# 6. Create raw UMI matrix
#
# Rows:
# - spatial barcodes
#
# Columns:
# - samples
# ==============================================================================

raw_umi_matrix <- vapply(
  umi_vectors,
  function(umi_vector) {

    as.numeric(
      umi_vector[
        shared_barcodes
      ]
    )

  },
  numeric(
    length(shared_barcodes)
  )
)


raw_umi_matrix <- as.matrix(
  raw_umi_matrix
)


rownames(
  raw_umi_matrix
) <- shared_barcodes


colnames(
  raw_umi_matrix
) <- names(
  umi_vectors
)


storage.mode(
  raw_umi_matrix
) <- "numeric"


# ==============================================================================
# 7. Validate raw UMI matrix
# ==============================================================================

if (any(!is.finite(raw_umi_matrix))) {

  stop(
    "raw_umi_matrix contains NA, NaN or Inf values."
  )

}


if (nrow(raw_umi_matrix) < 2) {

  stop(
    "Fewer than two shared barcodes are available."
  )

}


sample_variances <- apply(
  raw_umi_matrix,
  2,
  var
)


if (any(!is.finite(sample_variances))) {

  stop(
    "Non-finite sample variances were detected."
  )

}


if (any(sample_variances == 0)) {

  zero_variance_samples <- names(
    sample_variances[
      sample_variances == 0
    ]
  )

  stop(
    "Zero variance in raw UMI counts for samples: ",
    paste(
      zero_variance_samples,
      collapse = ", "
    )
  )

}


message(
  "Raw UMI matrix dimensions: ",
  nrow(raw_umi_matrix),
  " barcodes x ",
  ncol(raw_umi_matrix),
  " samples"
)


# ==============================================================================
# 8. Calculate Pearson and Spearman correlations
# ==============================================================================

pearson_matrix <- cor(
  x = raw_umi_matrix,
  use = "everything",
  method = "pearson"
)


spearman_matrix <- cor(
  x = raw_umi_matrix,
  use = "everything",
  method = "spearman"
)


# ==============================================================================
# 9. Validate correlation matrices
# ==============================================================================

if (any(!is.finite(pearson_matrix))) {

  stop(
    "Pearson correlation matrix contains NA, NaN or Inf values."
  )

}


if (any(!is.finite(spearman_matrix))) {

  stop(
    "Spearman correlation matrix contains NA, NaN or Inf values."
  )

}


diag(
  pearson_matrix
) <- 1


diag(
  spearman_matrix
) <- 1


pearson_range <- range(
  pearson_matrix[
    upper.tri(pearson_matrix)
  ],
  na.rm = TRUE
)


spearman_range <- range(
  spearman_matrix[
    upper.tri(spearman_matrix)
  ],
  na.rm = TRUE
)


message(
  "Pearson range excluding diagonal: ",
  sprintf(
    "%.4f to %.4f",
    pearson_range[1],
    pearson_range[2]
  )
)


message(
  "Spearman range excluding diagonal: ",
  sprintf(
    "%.4f to %.4f",
    spearman_range[1],
    spearman_range[2]
  )
)


# ==============================================================================
# 10. Save raw UMI matrix
# ==============================================================================

raw_umi_matrix_file <- file.path(
  output_dir,
  "01_raw_UMI_per_barcode_standard_samples.tsv"
)


write.table(
  x = raw_umi_matrix,
  file = raw_umi_matrix_file,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)


message(
  "Raw UMI matrix saved: ",
  raw_umi_matrix_file
)


# ==============================================================================
# 11. Save correlation matrices
# ==============================================================================

pearson_matrix_file <- file.path(
  output_dir,
  "01_raw_UMI_per_barcode_Pearson_matrix_standard_samples.tsv"
)


spearman_matrix_file <- file.path(
  output_dir,
  "01_raw_UMI_per_barcode_Spearman_matrix_standard_samples.tsv"
)


write.table(
  x = pearson_matrix,
  file = pearson_matrix_file,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)


write.table(
  x = spearman_matrix,
  file = spearman_matrix_file,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)


message(
  "Pearson matrix saved: ",
  pearson_matrix_file
)


message(
  "Spearman matrix saved: ",
  spearman_matrix_file
)


# ==============================================================================
# 12. Prepare heatmap color scale
#
# Breaks:
# - -1 to 0: blue to white
# -  0 to 1: white to red
#
# The negative and positive parts are defined separately so that:
# - only values below zero are blue
# - zero is exactly white
# - all values above zero are on the white-red scale
# ==============================================================================

negative_breaks <- seq(
  from = -1,
  to = 0,
  length.out = 101
)


positive_breaks <- seq(
  from = 0,
  to = 1,
  length.out = 101
)[-1]


heatmap_breaks <- c(
  negative_breaks,
  positive_breaks
)


negative_colors <- colorRampPalette(
  c(
    "navy",
    "white"
  )
)(
  length(negative_breaks) - 1
)


positive_colors <- colorRampPalette(
  c(
    "white",
    "firebrick3"
  )
)(
  length(positive_breaks)
)


heatmap_colors <- c(
  negative_colors,
  positive_colors
)


if (
  length(heatmap_colors) !=
  length(heatmap_breaks) - 1
) {

  stop(
    "Number of heatmap colors does not match number of intervals."
  )

}


# ==============================================================================
# 13. Pearson heatmap
# ==============================================================================

pearson_heatmap_file <- file.path(
  output_dir,
  "01_raw_UMI_per_barcode_Pearson_heatmap_standard_samples.pdf"
)


pheatmap(
  mat = pearson_matrix,

  color = heatmap_colors,
  breaks = heatmap_breaks,

  cluster_rows = TRUE,
  cluster_cols = TRUE,

  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",

  clustering_method = "complete",

  display_numbers = TRUE,
  number_format = "%.3f",

  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 9,
  fontsize_number = 6,

  border_color = "grey80",

  main = paste0(
    "Pearson correlation of raw UMI",
    "\nper spatial barcode"
  ),

  filename = pearson_heatmap_file,

  width = 12,
  height = 11
)


message(
  "Pearson heatmap saved: ",
  pearson_heatmap_file
)


# ==============================================================================
# 14. Spearman heatmap
# ==============================================================================

spearman_heatmap_file <- file.path(
  output_dir,
  "01_raw_UMI_per_barcode_Spearman_heatmap_standard_samples.pdf"
)


pheatmap(
  mat = spearman_matrix,

  color = heatmap_colors,
  breaks = heatmap_breaks,

  cluster_rows = TRUE,
  cluster_cols = TRUE,

  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",

  clustering_method = "complete",

  display_numbers = TRUE,
  number_format = "%.3f",

  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 9,
  fontsize_number = 6,

  border_color = "grey80",

  main = paste0(
    "Spearman correlation of raw UMI",
    "\nper spatial barcode"
  ),

  filename = spearman_heatmap_file,

  width = 12,
  height = 11
)


message(
  "Spearman heatmap saved: ",
  spearman_heatmap_file
)


# ==============================================================================
# 15. Final summary
# ==============================================================================

message(
  "Correlation analysis completed successfully."
)


message(
  "Output directory: ",
  output_dir
)