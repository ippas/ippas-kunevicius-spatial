# ==============================================================================
# 03_pseudobulk_edgeR_sexStratified_proteinCoding_heatmaps_log2Quantile_zScore_FIXED_v2.R
#
# PURPOSE
#   Create two separate heatmaps for genes selected from the sex-stratified
#   protein-coding edgeR analyses.
#
#     Male:
#       PValue < 0.01 and abs(logFC) > 0.5
#
#     Female:
#       FDR < 0.05 and abs(logFC) > 0.5
#
# EXPRESSION VALUES USED FOR VISUALIZATION
#   1. Start from raw in-tissue pseudobulk counts.
#   2. Subset samples separately for Male and Female.
#   3. Transform the complete protein-coding matrix as log2(raw count + 1).
#   4. Apply quantile normalization across samples within each sex.
#   5. Select the genes meeting the edgeR criteria.
#   6. Calculate a row-wise Z-score for every selected gene.
#   7. Clip Z-scores to the range from -2 to 2.
#
# FIXED VISUALIZATION
#   - sample names are displayed ABOVE the heatmap matrix
#   - columns are narrower
#   - row dendrogram is wider
#   - gene labels are larger
#   - color range is fixed at [-2, 2]
#   - expression and normalization descriptions are split across two lines
#   - legends are stacked vertically in one column
#   - PDF dimensions are fitted closely to the title, heatmap and legends
#   - the existing PDF filenames are retained and overwritten
#
# OUTPUT
#   results/maternalFMT_n16samples/
#     pseudobulk_geneCounts_edgeR_sexStratified/
#       edgeR_selectedGenes_heatmaps_log2Quantile_zScore/
#
# REQUIRED OBJECTS IN THE ACTIVE R SESSION
#   gene_counts_per_sample_raw_in_tissue
#   sample_metadata_analysis
#   edgeR_results_by_sex
# ==============================================================================


# ==============================================================================
# 1. Required packages
# ==============================================================================
required_packages <- c(
  "dplyr",
  "limma",
  "ComplexHeatmap",
  "circlize"
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
    "\nInstall missing packages before running this script.",
    "\nFor ComplexHeatmap use: BiocManager::install('ComplexHeatmap')",
    "\nFor circlize use: install.packages('circlize')"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(limma)
  library(ComplexHeatmap)
  library(circlize)
})


# ==============================================================================
# 2. Check required analysis objects
# ==============================================================================
required_objects <- c(
  "gene_counts_per_sample_raw_in_tissue",
  "sample_metadata_analysis",
  "edgeR_results_by_sex"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    envir = .GlobalEnv,
    inherits = FALSE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Missing required R objects: ",
    paste(missing_objects, collapse = ", "),
    "\nRun the protein-coding pseudobulk edgeR script first."
  )
}

if (!all(c("male", "female") %in% names(edgeR_results_by_sex))) {
  stop("edgeR_results_by_sex must contain both $male and $female.")
}


# ==============================================================================
# 3. Define output directory and KEEP THE SAME PDF FILENAMES
# ==============================================================================
project_dir <- "/home/mateusz/projects/ippas-kunevicius-spatial"

heatmap_output_dir <- file.path(
  project_dir,
  "results",
  "maternalFMT_n16samples",
  "pseudobulk_geneCounts_edgeR_sexStratified",
  "edgeR_selectedGenes_heatmaps_log2Quantile_zScore"
)

dir.create(
  heatmap_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(heatmap_output_dir)) {
  stop("Could not create heatmap output directory: ", heatmap_output_dir)
}

# IMPORTANT: these names are unchanged, so the old PDFs are overwritten.
male_heatmap_pdf <- file.path(
  heatmap_output_dir,
  paste0(
    "male_ASD_vs_Neurotypical_",
    "PValue_lt0.01_absLog2FC_gt0.5_",
    "log2Quantile_rowZscore_heatmap.pdf"
  )
)

female_heatmap_pdf <- file.path(
  heatmap_output_dir,
  paste0(
    "female_ASD_vs_Neurotypical_",
    "FDR_lt0.05_absLog2FC_gt0.5_",
    "log2Quantile_rowZscore_heatmap.pdf"
  )
)


# ==============================================================================
# 4. Validate the raw pseudobulk matrix
# ==============================================================================
raw_count_matrix <- as.matrix(
  gene_counts_per_sample_raw_in_tissue
)

storage.mode(raw_count_matrix) <- "numeric"

if (nrow(raw_count_matrix) == 0 || ncol(raw_count_matrix) == 0) {
  stop("The raw pseudobulk count matrix is empty.")
}

if (is.null(rownames(raw_count_matrix))) {
  stop("The raw pseudobulk count matrix has no gene identifiers.")
}

if (is.null(colnames(raw_count_matrix))) {
  stop("The raw pseudobulk count matrix has no sample identifiers.")
}

if (anyNA(raw_count_matrix)) {
  stop("The raw pseudobulk count matrix contains NA values.")
}

if (any(!is.finite(raw_count_matrix))) {
  stop("The raw pseudobulk count matrix contains non-finite values.")
}

if (any(raw_count_matrix < 0)) {
  stop("The raw pseudobulk count matrix contains negative values.")
}

if (anyDuplicated(rownames(raw_count_matrix))) {
  stop("Duplicated gene identifiers were found in the raw count matrix.")
}

if (anyDuplicated(colnames(raw_count_matrix))) {
  stop("Duplicated sample identifiers were found in the raw count matrix.")
}


# ==============================================================================
# 5. Validate analysis metadata
# ==============================================================================
required_metadata_columns <- c(
  "sample_ID",
  "fmt_donor_group",
  "sex"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(sample_metadata_analysis)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    "Missing columns in sample_metadata_analysis: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

if (anyDuplicated(sample_metadata_analysis$sample_ID)) {
  stop("Duplicated sample_ID values were found in sample_metadata_analysis.")
}

missing_metadata_samples <- setdiff(
  colnames(raw_count_matrix),
  as.character(sample_metadata_analysis$sample_ID)
)

if (length(missing_metadata_samples) > 0) {
  stop(
    "Samples from the count matrix are missing from metadata: ",
    paste(missing_metadata_samples, collapse = ", ")
  )
}


# ==============================================================================
# 6. Select genes separately for Male and Female
# ==============================================================================
male_selected_results <- edgeR_results_by_sex$male |>
  filter(
    PValue < 0.01,
    abs(logFC) > 0.5
  ) |>
  arrange(PValue)

female_selected_results <- edgeR_results_by_sex$female |>
  filter(
    FDR < 0.05,
    abs(logFC) > 0.5
  ) |>
  arrange(FDR, PValue)

if (nrow(male_selected_results) == 0) {
  stop("No male genes meet PValue < 0.01 and abs(logFC) > 0.5.")
}

if (nrow(female_selected_results) == 0) {
  stop("No female genes meet FDR < 0.05 and abs(logFC) > 0.5.")
}

message(
  "Male genes selected for heatmap: ",
  nrow(male_selected_results),
  " | PValue < 0.01 and abs(logFC) > 0.5"
)

message(
  "Female genes selected for heatmap: ",
  nrow(female_selected_results),
  " | FDR < 0.05 and abs(logFC) > 0.5"
)


# ==============================================================================
# 7. Prepare one sex-specific heatmap matrix
# ==============================================================================
prepare_heatmap_data <- function(
    sex_label,
    selected_results,
    raw_count_matrix,
    sample_metadata
) {

  sex_metadata <- sample_metadata |>
    filter(as.character(sex) == sex_label) |>
    mutate(
      fmt_donor_group = factor(
        as.character(fmt_donor_group),
        levels = c("Neurotypical", "ASD")
      )
    ) |>
    arrange(fmt_donor_group, sample_ID)

  if (nrow(sex_metadata) == 0) {
    stop("No metadata rows found for sex: ", sex_label)
  }

  if (anyNA(sex_metadata$fmt_donor_group)) {
    stop(
      "Unexpected fmt_donor_group values for sex: ",
      sex_label,
      ". Expected only Neurotypical and ASD."
    )
  }

  sex_sample_ids <- as.character(sex_metadata$sample_ID)

  missing_samples <- setdiff(
    sex_sample_ids,
    colnames(raw_count_matrix)
  )

  if (length(missing_samples) > 0) {
    stop(
      "Samples missing from the raw count matrix for ",
      sex_label,
      ": ",
      paste(missing_samples, collapse = ", ")
    )
  }

  raw_counts_sex <- raw_count_matrix[
    ,
    sex_sample_ids,
    drop = FALSE
  ]

  # Visualization processing only:
  # raw counts -> log2(count + 1) -> quantile normalization.
  log2_counts_sex <- log2(raw_counts_sex + 1)

  log2_quantile_counts_sex <- limma::normalizeBetweenArrays(
    log2_counts_sex,
    method = "quantile"
  )

  rownames(log2_quantile_counts_sex) <- rownames(raw_counts_sex)
  colnames(log2_quantile_counts_sex) <- colnames(raw_counts_sex)

  selected_gene_ids <- unique(
    as.character(selected_results$ensembl_gene_id)
  )

  missing_selected_genes <- setdiff(
    selected_gene_ids,
    rownames(log2_quantile_counts_sex)
  )

  if (length(missing_selected_genes) > 0) {
    stop(
      "Selected genes missing from the expression matrix for ",
      sex_label,
      ": ",
      paste(head(missing_selected_genes, 20), collapse = ", "),
      if (length(missing_selected_genes) > 20) " ..." else ""
    )
  }

  selected_expression <- log2_quantile_counts_sex[
    selected_gene_ids,
    ,
    drop = FALSE
  ]

  # Row-wise Z-score across samples.
  z_score_matrix <- t(
    scale(
      t(selected_expression),
      center = TRUE,
      scale = TRUE
    )
  )

  zero_variance_genes <- rownames(z_score_matrix)[
    apply(z_score_matrix, 1, function(x) any(!is.finite(x)))
  ]

  if (length(zero_variance_genes) > 0) {
    warning(
      sex_label,
      ": ",
      length(zero_variance_genes),
      " selected genes had zero variance after preprocessing; ",
      "their Z-scores were set to 0."
    )
  }

  z_score_matrix[!is.finite(z_score_matrix)] <- 0

  # FIXED requested range: -2 to 2.
  z_score_matrix_clipped <- pmax(
    pmin(z_score_matrix, 2),
    -2
  )

  selected_annotation <- selected_results |>
    distinct(ensembl_gene_id, .keep_all = TRUE) |>
    select(
      ensembl_gene_id,
      gene,
      regulation,
      logFC,
      PValue,
      FDR
    )

  selected_annotation <- selected_annotation[
    match(
      rownames(z_score_matrix_clipped),
      selected_annotation$ensembl_gene_id
    ),
    ,
    drop = FALSE
  ]

  row_labels <- ifelse(
    is.na(selected_annotation$gene) |
      selected_annotation$gene == "",
    selected_annotation$ensembl_gene_id,
    selected_annotation$gene
  )

  rownames(z_score_matrix_clipped) <- make.unique(row_labels)

  column_annotation <- data.frame(
    Group = factor(
      as.character(sex_metadata$fmt_donor_group),
      levels = c("Neurotypical", "ASD")
    ),
    row.names = sex_sample_ids,
    check.names = FALSE
  )

  list(
    sex = sex_label,
    metadata = sex_metadata,
    selected_results = selected_results,
    selected_annotation = selected_annotation,
    raw_counts = raw_counts_sex,
    log2_counts = log2_counts_sex,
    log2_quantile_counts_all_genes = log2_quantile_counts_sex,
    log2_quantile_counts_selected_genes = selected_expression,
    z_score_matrix = z_score_matrix,
    z_score_matrix_clipped = z_score_matrix_clipped,
    column_annotation = column_annotation,
    zero_variance_genes = zero_variance_genes
  )
}


# ==============================================================================
# 8. Prepare Male and Female heatmap data
# ==============================================================================
male_heatmap_data <- prepare_heatmap_data(
  sex_label = "Male",
  selected_results = male_selected_results,
  raw_count_matrix = raw_count_matrix,
  sample_metadata = sample_metadata_analysis
)

female_heatmap_data <- prepare_heatmap_data(
  sex_label = "Female",
  selected_results = female_selected_results,
  raw_count_matrix = raw_count_matrix,
  sample_metadata = sample_metadata_analysis
)

heatmap_data_by_sex <- list(
  male = male_heatmap_data,
  female = female_heatmap_data
)

heatmap_selected_genes_by_sex <- list(
  male = male_selected_results,
  female = female_selected_results
)


# ==============================================================================
# 9. Fixed color scale and annotation colors
# ==============================================================================
heatmap_color_function <- circlize::colorRamp2(
  breaks = c(-2, 0, 2),
  colors = c("blue", "white", "red")
)

annotation_colors <- list(
  Group = c(
    Neurotypical = "grey70",
    ASD = "black"
  )
)


# ==============================================================================
# 10. Correlation distance for clustering genes
# ==============================================================================
row_correlation_distance <- function(expression_matrix) {

  correlation_matrix <- stats::cor(
    t(expression_matrix),
    method = "pearson",
    use = "pairwise.complete.obs"
  )

  correlation_matrix[!is.finite(correlation_matrix)] <- 0
  diag(correlation_matrix) <- 1

  stats::as.dist(1 - correlation_matrix)
}


# ==============================================================================
# 11. Save one ComplexHeatmap to PDF
# ==============================================================================
save_sex_heatmap <- function(
    heatmap_data,
    output_pdf,
    selection_label
) {

  heatmap_matrix <- heatmap_data$z_score_matrix_clipped

  number_of_genes <- nrow(heatmap_matrix)
  number_of_samples <- ncol(heatmap_matrix)

  if (number_of_genes == 0 || number_of_samples == 0) {
    stop("Heatmap matrix is empty for sex: ", heatmap_data$sex)
  }

  # Gene-label size is deliberately larger than in the previous version.
  row_font_size <- if (number_of_genes <= 70) {
    9
  } else if (number_of_genes <= 120) {
    7.5
  } else {
    6.5
  }

  # Row height is selected so the larger labels remain readable.
  row_height_mm <- if (number_of_genes <= 70) {
    3.8
  } else if (number_of_genes <= 120) {
    3.2
  } else {
    2.7
  }

  heatmap_body_height_mm <- number_of_genes * row_height_mm

  # Narrower columns: 6.5 mm per sample.
  heatmap_body_width_mm <- number_of_samples * 6.5

  # A substantially wider row dendrogram than before.
  row_dendrogram_width <- grid::unit(28, "mm")

  # --------------------------------------------------------------------------
  # Tight PDF dimensions
  # --------------------------------------------------------------------------
  # Estimate the space needed by gene labels, titles, dendrogram and legends.
  # This avoids the large white side margins produced by a fixed 9.5-inch page.

  temporary_measurement_device <- FALSE

  if (grDevices::dev.cur() == 1L) {
    grDevices::pdf(NULL)
    temporary_measurement_device <- TRUE
  }

  row_names_width_mm <- grid::convertWidth(
    ComplexHeatmap::max_text_width(
      rownames(heatmap_matrix),
      gp = grid::gpar(fontsize = row_font_size)
    ),
    unitTo = "mm",
    valueOnly = TRUE
  )

  title_lines_for_measurement <- c(
    paste0(heatmap_data$sex, ": ASD vs Neurotypical"),
    paste0("Selected genes: ", selection_label),
    "Expression: log2(raw in-tissue pseudobulk counts + 1)",
    "Normalization: quantile normalization within sex",
    paste0(
      "Row Z-score clipped to [-2, 2] | ",
      number_of_genes,
      " genes | ",
      number_of_samples,
      " samples"
    )
  )

  title_width_mm <- grid::convertWidth(
    ComplexHeatmap::max_text_width(
      title_lines_for_measurement,
      gp = grid::gpar(fontsize = 13, fontface = "bold")
    ),
    unitTo = "mm",
    valueOnly = TRUE
  )

  if (temporary_measurement_device) {
    grDevices::dev.off()
  }

  # Width components: dendrogram + matrix + row names + vertically stacked
  # legends + small internal gaps and outer padding.
  legend_column_width_mm <- 30
  internal_spacing_width_mm <- 18

  plot_content_width_mm <-
    28 +
    heatmap_body_width_mm +
    row_names_width_mm +
    legend_column_width_mm +
    internal_spacing_width_mm

  pdf_width_mm <- max(
    title_width_mm + 8,
    plot_content_width_mm
  )

  # Height components: title block, column labels, top annotation, heatmap body
  # and small outer margins. No large fixed minimum is used.
  title_block_height_mm <- 33
  column_names_height_mm <- 18
  top_annotation_height_mm <- 6
  outer_vertical_padding_mm <- 8

  pdf_height_mm <-
    title_block_height_mm +
    column_names_height_mm +
    top_annotation_height_mm +
    heatmap_body_height_mm +
    outer_vertical_padding_mm

  pdf_width_inches <- pdf_width_mm / 25.4
  pdf_height_inches <- pdf_height_mm / 25.4

  complete_title <- paste0(
    heatmap_data$sex,
    ": ASD vs Neurotypical",
    "\nSelected genes: ",
    selection_label,
    "\nExpression: log2(raw in-tissue pseudobulk counts + 1)",
    "\nNormalization: quantile normalization within sex",
    "\nRow Z-score clipped to [-2, 2] | ",
    number_of_genes,
    " genes | ",
    number_of_samples,
    " samples"
  )

  group_values <- heatmap_data$column_annotation$Group
  names(group_values) <- rownames(heatmap_data$column_annotation)

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Group = group_values,
    col = annotation_colors,
    simple_anno_size = grid::unit(4.5, "mm"),
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(
      fontsize = 10,
      fontface = "bold"
    ),
    annotation_legend_param = list(
      Group = list(
        title = "Group",
        direction = "vertical",
        ncol = 1,
        title_gp = grid::gpar(
          fontsize = 10,
          fontface = "bold"
        ),
        labels_gp = grid::gpar(fontsize = 9)
      )
    )
  )

  heatmap_object <- ComplexHeatmap::Heatmap(
    matrix = heatmap_matrix,
    name = "Row Z-score",
    col = heatmap_color_function,

    # Row clustering and larger dendrogram.
    cluster_rows = TRUE,
    clustering_distance_rows = row_correlation_distance,
    clustering_method_rows = "complete",
    row_dend_width = row_dendrogram_width,
    show_row_dend = TRUE,

    # Keep the samples in Neurotypical -> ASD order.
    cluster_columns = FALSE,

    # REQUIRED: sample names at the TOP of the heatmap matrix.
    show_column_names = TRUE,
    column_names_side = "top",
    column_names_rot = 45,
    column_names_centered = TRUE,
    column_names_gp = grid::gpar(
      fontsize = 11,
      fontface = "bold"
    ),
    column_names_max_height = grid::unit(18, "mm"),

    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = row_font_size),
    row_names_max_width = ComplexHeatmap::max_text_width(
      rownames(heatmap_matrix),
      gp = grid::gpar(fontsize = row_font_size)
    ) + grid::unit(3, "mm"),

    top_annotation = top_annotation,

    # Narrow heatmap body and explicit row height.
    width = grid::unit(heatmap_body_width_mm, "mm"),
    height = grid::unit(heatmap_body_height_mm, "mm"),

    rect_gp = grid::gpar(col = NA),
    border = FALSE,
    use_raster = FALSE,

    column_title = complete_title,
    column_title_side = "top",
    column_title_gp = grid::gpar(
      fontsize = 13,
      fontface = "bold",
      lineheight = 1.15
    ),

    heatmap_legend_param = list(
      title = "Row Z-score",
      direction = "vertical",
      at = c(-2, -1, 0, 1, 2),
      labels = c("-2", "-1", "0", "1", "2"),
      title_gp = grid::gpar(
        fontsize = 10,
        fontface = "bold"
      ),
      labels_gp = grid::gpar(fontsize = 9),
      legend_height = grid::unit(48, "mm")
    )
  )

  # pdf() overwrites an existing PDF with the same filename.
  grDevices::pdf(
    file = output_pdf,
    width = pdf_width_inches,
    height = pdf_height_inches,
    onefile = TRUE,
    useDingbats = FALSE
  )

  ComplexHeatmap::draw(
    heatmap_object,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",

    # Put both legends in one vertical column: Row Z-score above Group.
    merge_legends = TRUE,
    legend_grouping = "original",

    # Tight but safe margins around all heatmap components.
    padding = grid::unit(c(3, 3, 3, 3), "mm")
  )

  grDevices::dev.off()

  if (!file.exists(output_pdf)) {
    stop("Heatmap PDF was not created: ", output_pdf)
  }

  message(
    heatmap_data$sex,
    " heatmap saved and overwritten: ",
    normalizePath(output_pdf, mustWork = TRUE)
  )

  list(
    heatmap_object = heatmap_object,
    output_pdf = output_pdf,
    number_of_genes = number_of_genes,
    number_of_samples = number_of_samples,
    z_score_range = c(-2, 2),
    column_names_side = "top",
    heatmap_body_width_mm = heatmap_body_width_mm,
    row_dendrogram_width_mm = 28,
    pdf_width_inches = pdf_width_inches,
    pdf_height_inches = pdf_height_inches,
    legends_layout = "one vertical column"
  )
}


# ==============================================================================
# 12. Save the two independent PDF heatmaps
# ==============================================================================
male_heatmap_plot <- save_sex_heatmap(
  heatmap_data = male_heatmap_data,
  output_pdf = male_heatmap_pdf,
  selection_label = "PValue < 0.01 and |log2FC| > 0.5"
)

female_heatmap_plot <- save_sex_heatmap(
  heatmap_data = female_heatmap_data,
  output_pdf = female_heatmap_pdf,
  selection_label = "FDR < 0.05 and |log2FC| > 0.5"
)

heatmap_plots_by_sex <- list(
  male = male_heatmap_plot,
  female = female_heatmap_plot
)

heatmap_output_files <- c(
  male = male_heatmap_pdf,
  female = female_heatmap_pdf
)


# ==============================================================================
# 13. Final summary
# ==============================================================================
heatmap_summary <- tibble::tibble(
  sex = c("Male", "Female"),
  comparison = "ASD_vs_Neurotypical",
  selection_criterion = c(
    "PValue < 0.01 and abs(logFC) > 0.5",
    "FDR < 0.05 and abs(logFC) > 0.5"
  ),
  selected_genes = c(
    nrow(male_selected_results),
    nrow(female_selected_results)
  ),
  samples = c(
    ncol(male_heatmap_data$z_score_matrix_clipped),
    ncol(female_heatmap_data$z_score_matrix_clipped)
  ),
  transformation = "log2(raw pseudobulk counts + 1)",
  normalization = "quantile normalization within sex",
  visualization = "row Z-score clipped to [-2, 2]",
  column_names_side = "top",
  title_expression_layout = "expression and normalization on separate lines",
  legends_layout = "Row Z-score above Group in one vertical column",
  output_pdf = c(
    male_heatmap_pdf,
    female_heatmap_pdf
  )
)

print(heatmap_summary, n = Inf, width = Inf)

message("\nHeatmap output directory:")
message(normalizePath(heatmap_output_dir, mustWork = TRUE))
message("\nHeatmap analysis completed successfully.")


# ==============================================================================
# End
# ==============================================================================
