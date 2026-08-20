# ==============================================================================
# 06_log2_UMI_maps_20_samples_common_scale_bottomLegend_blueYellowRed.R
#
# Purpose:
# - read 20 standard Visium samples as separate Seurat objects
# - calculate raw UMI statistics separately for every sample
# - create spatial maps using log2(UMI + 1)
# - use one common color scale for all samples: 0 to global max across all samples
# - place one shared horizontal legend at the top of the figure
# - save all 20 plots into one PDF
#
# Notes:
# - only standard samples are included
# - no Undetermined FASTQ analyses are used here
# - raw UMI values are not modified; only the visualization uses log2(UMI + 1)
# ==============================================================================

# ==============================================================================
# 1. Packages
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(scales)
  library(viridisLite)
})

# ==============================================================================
# 2. Project paths
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
  "maternalFMT_n20samples",
  "umi_per_slide"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_pdf <- file.path(
  output_dir,
  "06_log2_UMI_maps_20_samples_common_scale_bottomLegend_blueYellowRed.pdf"
)

output_summary_tsv <- file.path(
  output_dir,
  "04_log2_UMI_summary_20_samples_common_scale_0_to_globalMax.tsv"
)

# ==============================================================================
# 3. Plot and PDF settings
# ==============================================================================
number_of_columns <- 4
pdf_width_inches <- 20
pdf_height_inches <- 29
spot_size_factor <- 1.6

# Common blue -> yellow -> red palette: low UMI = blue, high UMI = red.
palette_colors <- grDevices::colorRampPalette(
  c(
    "#08306B",
    "#2171B5",
    "#41B6C4",
    "#FFFFBF",
    "#FDAE61",
    "#D73027"
  )
)(256)

# ==============================================================================
# 4. Check input paths
# ==============================================================================
if (!dir.exists(project_dir)) {
  stop("Project directory does not exist: ", project_dir)
}

if (!dir.exists(path_to_data)) {
  stop("Standard Space Ranger data directory does not exist: ", path_to_data)
}

if (!file.exists(metadata_file)) {
  stop("Metadata file does not exist: ", metadata_file)
}

message("Project directory: ", project_dir)
message("Standard Space Ranger data directory: ", path_to_data)
message("Metadata file: ", metadata_file)
message("Output directory: ", output_dir)

# ==============================================================================
# 5. Read metadata
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

# ==============================================================================
# 6. Prepare metadata for the 20 standard samples
# ==============================================================================
metadata_samples <- metadata_autismFMT |>
  transmute(
    metadata_sample_ID = trimws(as.character(sample_ID)),
    group = trimws(as.character(fmt_donor_group)),
    sex = trimws(as.character(sex))
  ) |>
  filter(!is.na(metadata_sample_ID), metadata_sample_ID != "") |>
  distinct(metadata_sample_ID, .keep_all = TRUE)

if (nrow(metadata_samples) == 0) {
  stop("No samples were found in the metadata file.")
}

if (anyDuplicated(metadata_samples$metadata_sample_ID) > 0) {
  stop("Duplicated sample IDs were detected after metadata preparation.")
}

if (nrow(metadata_samples) != 20) {
  stop(
    "Expected exactly 20 standard samples in metadata, but found: ",
    nrow(metadata_samples)
  )
}

message("Number of standard samples in metadata: ", nrow(metadata_samples))

# ==============================================================================
# 7. Define Space Ranger sources for standard samples
# ==============================================================================
standard_sample_sources <- metadata_samples |>
  transmute(
    sample_ID = metadata_sample_ID,
    metadata_sample_ID = metadata_sample_ID,
    sample_type = "standard",
    outs_dir = file.path(path_to_data, metadata_sample_ID, "outs")
  )

sample_sources <- standard_sample_sources
selected_samples <- sample_sources$sample_ID

if (length(selected_samples) == 0) {
  stop("No samples were selected.")
}

if (anyDuplicated(selected_samples) > 0) {
  duplicated_sample_ids <- unique(selected_samples[duplicated(selected_samples)])
  stop(
    "Duplicated plotting sample IDs were detected: ",
    paste(duplicated_sample_ids, collapse = ", ")
  )
}

sample_sources <- sample_sources |>
  left_join(metadata_samples, by = "metadata_sample_ID")

missing_group_samples <- sample_sources$sample_ID[
  is.na(sample_sources$group) | trimws(sample_sources$group) == ""
]

if (length(missing_group_samples) > 0) {
  stop(
    "Missing FMT donor group for samples: ",
    paste(missing_group_samples, collapse = ", ")
  )
}

missing_sex_samples <- sample_sources$sample_ID[
  is.na(sample_sources$sex) | trimws(sample_sources$sex) == ""
]

if (length(missing_sex_samples) > 0) {
  stop(
    "Missing sex for samples: ",
    paste(missing_sex_samples, collapse = ", ")
  )
}

number_of_rows <- ceiling(length(selected_samples) / number_of_columns)

message("Total number of plotted samples: ", length(selected_samples))
message("PDF layout: ", number_of_columns, " columns x ", number_of_rows, " rows")
message("Samples: ", paste(selected_samples, collapse = ", "))

# ==============================================================================
# 8. Check Space Ranger matrices and spatial folders
# ==============================================================================
matrix_files <- file.path(sample_sources$outs_dir, "raw_feature_bc_matrix.h5")
missing_matrix_samples <- sample_sources$sample_ID[!file.exists(matrix_files)]

if (length(missing_matrix_samples) > 0) {
  missing_matrix_information <- sample_sources |>
    filter(sample_ID %in% missing_matrix_samples) |>
    transmute(
      information = paste0(sample_ID, ": ", file.path(outs_dir, "raw_feature_bc_matrix.h5"))
    ) |>
    pull(information)

  stop(
    "Missing raw_feature_bc_matrix.h5 for samples:\n",
    paste(missing_matrix_information, collapse = "\n")
  )
}

spatial_directories <- file.path(sample_sources$outs_dir, "spatial")
missing_spatial_samples <- sample_sources$sample_ID[!dir.exists(spatial_directories)]

if (length(missing_spatial_samples) > 0) {
  missing_spatial_information <- sample_sources |>
    filter(sample_ID %in% missing_spatial_samples) |>
    transmute(
      information = paste0(sample_ID, ": ", file.path(outs_dir, "spatial"))
    ) |>
    pull(information)

  stop(
    "Missing spatial directories for samples:\n",
    paste(missing_spatial_information, collapse = "\n")
  )
}

message("All required Space Ranger matrices and spatial folders were found.")

# ==============================================================================
# 9. Read samples as separate Seurat objects
# ==============================================================================
samples_list <- setNames(
  map(
    selected_samples,
    function(sample_id) {
      message("Reading sample: ", sample_id)

      sample_row <- sample_sources |>
        filter(sample_ID == sample_id)

      if (nrow(sample_row) != 1) {
        stop(
          "Expected exactly one input path for sample: ",
          sample_id,
          ". Found: ",
          nrow(sample_row)
        )
      }

      sample_outs_dir <- sample_row$outs_dir[[1]]

      seurat_object <- Load10X_Spatial(
        data.dir = sample_outs_dir,
        filename = "raw_feature_bc_matrix.h5",
        assay = "RNA",
        slice = sample_id,
        filter.matrix = FALSE
      )

      if (!"nCount_RNA" %in% colnames(seurat_object[[]])) {
        stop("nCount_RNA is missing in Seurat metadata for sample: ", sample_id)
      }

      if (length(Images(seurat_object)) == 0) {
        stop("No histological image was loaded for sample: ", sample_id)
      }

      if (ncol(seurat_object) == 0) {
        stop("The Seurat object contains no spots for sample: ", sample_id)
      }

      seurat_object$plot_sample_ID <- sample_id
      seurat_object$metadata_sample_ID <- sample_row$metadata_sample_ID[[1]]
      seurat_object$sample_type <- sample_row$sample_type[[1]]
      seurat_object$group <- sample_row$group[[1]]
      seurat_object$sex <- sample_row$sex[[1]]
      seurat_object$log2_nCount_RNA <- log2(seurat_object$nCount_RNA + 1)

      seurat_object
    }
  ),
  selected_samples
)

if (length(samples_list) != length(selected_samples)) {
  stop("The number of Seurat objects does not match the number of selected samples.")
}

message("All samples were loaded successfully.")

# ==============================================================================
# 10. Calculate raw UMI statistics and log2 scale summaries
# ==============================================================================
umi_summary <- imap_dfr(
  samples_list,
  function(seurat_object, sample_id) {
    umi <- as.numeric(seurat_object[[]][["nCount_RNA"]])
    log2_umi <- as.numeric(seurat_object[[]][["log2_nCount_RNA"]])

    if (length(umi) == 0) {
      stop("No UMI values found for sample: ", sample_id)
    }

    if (all(is.na(umi))) {
      stop("All UMI values are NA for sample: ", sample_id)
    }

    tibble(
      sample_ID = sample_id,
      n_spots = length(umi),
      n_spots_with_UMI = sum(umi > 0, na.rm = TRUE),
      percent_spots_with_UMI = 100 * sum(umi > 0, na.rm = TRUE) / length(umi),
      total_UMI = sum(umi, na.rm = TRUE),
      mean_UMI_per_spot = mean(umi, na.rm = TRUE),
      median_UMI_per_spot = median(umi, na.rm = TRUE),
      sd_UMI_per_spot = sd(umi, na.rm = TRUE),
      min_UMI_per_spot = min(umi, na.rm = TRUE),
      max_UMI_per_spot = max(umi, na.rm = TRUE),
      min_log2_UMI_per_spot = min(log2_umi, na.rm = TRUE),
      max_log2_UMI_per_spot = max(log2_umi, na.rm = TRUE)
    )
  }
) |>
  left_join(
    sample_sources |>
      select(sample_ID, metadata_sample_ID, sample_type, group, sex, outs_dir),
    by = "sample_ID"
  ) |>
  select(
    sample_ID,
    metadata_sample_ID,
    sample_type,
    group,
    sex,
    everything()
  )

if (nrow(umi_summary) != length(samples_list)) {
  stop("The number of UMI summary rows does not match the number of loaded Seurat objects.")
}

if (anyDuplicated(umi_summary$sample_ID) > 0) {
  stop("Duplicated sample IDs detected in umi_summary.")
}

message("UMI summary:")
print(umi_summary, n = Inf, width = Inf)

# ==============================================================================
# 11. Save summary table
# ==============================================================================
write.table(
  umi_summary,
  file = output_summary_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

if (!file.exists(output_summary_tsv)) {
  stop("The UMI summary file was not created: ", output_summary_tsv)
}

message("UMI summary saved: ")
message(normalizePath(output_summary_tsv, mustWork = TRUE))

# ==============================================================================
# 12. Global common color scale for all samples
# ==============================================================================
all_log2_values <- unlist(
  lapply(
    samples_list,
    function(seurat_object) {
      as.numeric(seurat_object[[]][["log2_nCount_RNA"]])
    }
  ),
  use.names = FALSE
)

if (length(all_log2_values) == 0 || all(is.na(all_log2_values))) {
  stop("Global log2(UMI + 1) values could not be calculated.")
}

global_min_log2 <- 0
global_max_log2 <- max(all_log2_values, na.rm = TRUE)

message("Common log2(UMI + 1) scale:")
message("Minimum: ", round(global_min_log2, 3))
message("Maximum: ", round(global_max_log2, 3))

# ==============================================================================
# 13. Helper function for displaying metadata
# ==============================================================================
format_metadata_value <- function(x) {
  if (length(x) == 0 || is.na(x) || trimws(as.character(x)) == "") {
    return("NA")
  }
  as.character(x)
}

# ==============================================================================
# 14. Create log2 UMI maps
# ==============================================================================
umi_plots <- imap(
  samples_list,
  function(seurat_object, sample_id) {
    message("Creating plot for sample: ", sample_id)

    summary_row_number <- match(sample_id, umi_summary$sample_ID)
    if (is.na(summary_row_number)) {
      stop("No UMI summary found for sample: ", sample_id)
    }

    stats_row <- umi_summary[summary_row_number, , drop = FALSE]

    image_names <- Images(seurat_object)
    if (length(image_names) == 0) {
      stop("No spatial image found for sample: ", sample_id)
    }

    image_name <- image_names[1]
    group_text <- format_metadata_value(stats_row$group)
    sex_text <- format_metadata_value(stats_row$sex)

    subtitle_text <- paste0(
      "Group: ", group_text,
      " | Sex: ", sex_text,
      "\nRaw total UMI: ",
      format(stats_row$total_UMI, big.mark = " ", scientific = FALSE, trim = TRUE),
      " | spots: ",
      format(stats_row$n_spots, big.mark = " ", scientific = FALSE, trim = TRUE),
      "\nSpots with UMI: ",
      format(stats_row$n_spots_with_UMI, big.mark = " ", scientific = FALSE, trim = TRUE),
      " (",
      format(round(stats_row$percent_spots_with_UMI, 1), nsmall = 1, trim = TRUE),
      "%)",
      "\nMean raw UMI: ",
      format(round(stats_row$mean_UMI_per_spot, 1), big.mark = " ", scientific = FALSE, trim = TRUE),
      " | median: ",
      format(round(stats_row$median_UMI_per_spot, 1), big.mark = " ", scientific = FALSE, trim = TRUE),
      " | max: ",
      format(round(stats_row$max_UMI_per_spot, 0), big.mark = " ", scientific = FALSE, trim = TRUE)
    )

    p <- SpatialFeaturePlot(
      object = seurat_object,
      features = "log2_nCount_RNA",
      images = image_name,
      image.alpha = 1,
      crop = FALSE,
      pt.size.factor = spot_size_factor,
      min.cutoff = global_min_log2,
      max.cutoff = global_max_log2
    )

    p +
      labs(
        title = paste0("Sample: ", sample_id),
        subtitle = subtitle_text,
        fill = "log2(UMI + 1)",
        colour = "log2(UMI + 1)"
      ) +
      scale_fill_gradientn(
        colours = palette_colors,
        limits = c(global_min_log2, global_max_log2),
        oob = scales::squish
      ) +
      scale_colour_gradientn(
        colours = palette_colors,
        limits = c(global_min_log2, global_max_log2),
        oob = scales::squish
      ) +
      guides(
        fill = guide_colorbar(
          title = "log2(UMI + 1)",
          title.position = "top",
          direction = "horizontal",
          barwidth = grid::unit(70, "mm"),
          barheight = grid::unit(6.5, "mm")
        ),
        colour = guide_colorbar(
          title = "log2(UMI + 1)",
          title.position = "top",
          direction = "horizontal",
          barwidth = grid::unit(70, "mm"),
          barheight = grid::unit(6.5, "mm")
        )
      ) +
      theme(
        text = element_text(family = "DejaVu Sans"),
        plot.title = element_text(size = 13.5, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 7.3, hjust = 0.5, lineheight = 1.06),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.margin = margin(t = 10, r = 6, b = 6, l = 6)
      )
  }
)

if (length(umi_plots) != length(selected_samples)) {
  stop("The number of plots does not match the number of selected samples.")
}

# ==============================================================================
# 15. Combine all maps with shared legend
# ==============================================================================
all_umi_plot <- wrap_plots(
  umi_plots,
  ncol = number_of_columns,
  guides = "collect"
) +
  plot_annotation(
    title = "log2(UMI + 1) maps for maternal FMT Visium samples",
    subtitle = paste0(
      "20 standard Visium samples",
      " | one common color scale: 0–",
      round(global_max_log2, 3),
      " | total: ",
      length(selected_samples),
      " samples",
      " | layout: ",
      number_of_columns,
      " columns × ",
      number_of_rows,
      " rows"
    ),
    theme = theme(
      text = element_text(family = "DejaVu Sans"),
      plot.title = element_text(size = 21, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      plot.margin = margin(t = 8, r = 8, b = 6, l = 8)
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.margin = margin(t = 10, r = 0, b = 2, l = 0)
  )

# ==============================================================================
# 16. Display combined plot
# ==============================================================================
print(all_umi_plot)

# ==============================================================================
# 17. Select PDF device
# ==============================================================================
if (capabilities("cairo")) {
  pdf_device <- grDevices::cairo_pdf
  message("Using cairo_pdf device.")
} else {
  pdf_device <- grDevices::pdf
  warning("Cairo PDF is not available. Using the standard PDF device.")
}

# ==============================================================================
# 18. Save combined plot as PDF
# ==============================================================================
message("Saving PDF...")

ggsave(
  filename = output_pdf,
  plot = all_umi_plot,
  device = pdf_device,
  width = pdf_width_inches,
  height = pdf_height_inches,
  units = "in",
  limitsize = FALSE
)

# ==============================================================================
# 19. Validate output
# ==============================================================================
if (!file.exists(output_pdf)) {
  stop("The PDF file was not created: ", output_pdf)
}

output_pdf_size_mb <- file.info(output_pdf)$size / 1024^2

message("PDF saved successfully:")
message(normalizePath(output_pdf, mustWork = TRUE))
message("PDF size: ", round(output_pdf_size_mb, 2), " MB")
message("UMI summary:")
message(normalizePath(output_summary_tsv, mustWork = TRUE))
message("Analysis completed successfully.")
# ==============================================================================
# End
# ==============================================================================
