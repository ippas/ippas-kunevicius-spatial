# ==============================================================================
# 01_raw_UMI_maps_all_samples_plus_undetermined.R
#
# Purpose:
# - read 20 standard Visium samples as separate Seurat objects
# - read 4 additional Space Ranger outputs generated from Undetermined FASTQ
# - calculate raw UMI statistics separately for every sample
# - plot raw UMI over histological images
# - add sample ID, FMT donor group, sex, input type and UMI summary
# - save all 24 plots into one PDF
#
# Output layout:
# - 4 columns
# - 6 rows for 24 samples
#
# No normalization, integration, clustering or correlation analysis.
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
  "01_raw_UMI_maps_20_samples_plus_4_undetermined_4_columns.pdf"
)


output_summary_tsv <- file.path(
  output_dir,
  "01_raw_UMI_summary_20_samples_plus_4_undetermined.tsv"
)


# ==============================================================================
# 3. Plot and PDF settings
# ==============================================================================

number_of_columns <- 4


# Large custom PDF page suitable for 24 samples.
pdf_width_inches <- 20
pdf_height_inches <- 33


# Spatial spot size.
spot_size_factor <- 1.6


# ==============================================================================
# 4. Check input paths
# ==============================================================================

if (!dir.exists(project_dir)) {

  stop(
    "Project directory does not exist: ",
    project_dir
  )

}


if (!dir.exists(path_to_data)) {

  stop(
    "Standard Space Ranger data directory does not exist: ",
    path_to_data
  )

}


if (!file.exists(metadata_file)) {

  stop(
    "Metadata file does not exist: ",
    metadata_file
  )

}


message(
  "Project directory: ",
  project_dir
)

message(
  "Standard Space Ranger data directory: ",
  path_to_data
)

message(
  "Metadata file: ",
  metadata_file
)

message(
  "Output directory: ",
  output_dir
)


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
    paste(
      missing_metadata_columns,
      collapse = ", "
    )
  )

}


# ==============================================================================
# 6. Prepare metadata for the original 20 samples
# ==============================================================================

metadata_samples <- metadata_autismFMT |>
  transmute(

    metadata_sample_ID = trimws(
      as.character(sample_ID)
    ),

    group = trimws(
      as.character(fmt_donor_group)
    ),

    sex = trimws(
      as.character(sex)
    )

  ) |>
  filter(
    !is.na(metadata_sample_ID),
    metadata_sample_ID != ""
  ) |>
  distinct(
    metadata_sample_ID,
    .keep_all = TRUE
  )


if (nrow(metadata_samples) == 0) {

  stop(
    "No samples were found in the metadata file."
  )

}


if (anyDuplicated(metadata_samples$metadata_sample_ID) > 0) {

  stop(
    "Duplicated sample IDs were detected after metadata preparation."
  )

}


message(
  "Number of original samples in metadata: ",
  nrow(metadata_samples)
)


# ==============================================================================
# 7. Define Space Ranger sources for standard samples
# ==============================================================================

standard_sample_sources <- metadata_samples |>
  transmute(

    sample_ID = metadata_sample_ID,

    metadata_sample_ID = metadata_sample_ID,

    sample_type = "standard",

    input_description = "Standard FASTQ",

    outs_dir = file.path(
      path_to_data,
      metadata_sample_ID,
      "outs"
    )

  )


# ==============================================================================
# 8. Define Space Ranger sources for Undetermined samples
# ==============================================================================

undetermined_sample_sources <- tribble(

  ~sample_ID,
  ~metadata_sample_ID,
  ~sample_type,
  ~input_description,
  ~outs_dir,

  "12_3F_undetermined",
  "12_3F",
  "undetermined",
  "Undetermined FASTQ, lane L002",
  file.path(
    project_dir,
    "data",
    "spacerangerCount_undetermined_12_3F_L002_V12N14_348_C1",
    "outs"
  ),

  "15_1M_undetermined",
  "15_1M",
  "undetermined",
  "Undetermined FASTQ, lane L001",
  file.path(
    project_dir,
    "data",
    "spacerangerCount_undetermined_15_1M_L001_V12N14_349_A1",
    "outs"
  ),

  "20_1F_undetermined",
  "20_1F",
  "undetermined",
  "Undetermined FASTQ, lane L001",
  file.path(
    project_dir,
    "data",
    "spacerangerCount_undetermined_20_1F_L001_V12N14_348_A1",
    "outs"
  ),

  "20_3M_undetermined",
  "20_3M",
  "undetermined",
  "Undetermined FASTQ, lane L002",
  file.path(
    project_dir,
    "data",
    "spacerangerCount_undetermined_20_3M_L002_V12N14_349_D1",
    "outs"
  )

)


# ==============================================================================
# 9. Combine sample sources with metadata
# ==============================================================================

sample_sources <- bind_rows(
  standard_sample_sources,
  undetermined_sample_sources
) |>
  left_join(
    metadata_samples,
    by = "metadata_sample_ID"
  )


selected_samples <- sample_sources$sample_ID


if (length(selected_samples) == 0) {

  stop(
    "No samples were selected."
  )

}


if (anyDuplicated(selected_samples) > 0) {

  duplicated_sample_ids <- unique(
    selected_samples[
      duplicated(selected_samples)
    ]
  )

  stop(
    "Duplicated plotting sample IDs were detected: ",
    paste(
      duplicated_sample_ids,
      collapse = ", "
    )
  )

}


missing_group_samples <- sample_sources$sample_ID[
  is.na(sample_sources$group) |
    trimws(sample_sources$group) == ""
]


if (length(missing_group_samples) > 0) {

  stop(
    "Missing FMT donor group for samples: ",
    paste(
      missing_group_samples,
      collapse = ", "
    )
  )

}


missing_sex_samples <- sample_sources$sample_ID[
  is.na(sample_sources$sex) |
    trimws(sample_sources$sex) == ""
]


if (length(missing_sex_samples) > 0) {

  stop(
    "Missing sex for samples: ",
    paste(
      missing_sex_samples,
      collapse = ", "
    )
  )

}


number_of_rows <- ceiling(
  length(selected_samples) / number_of_columns
)


message(
  "Total number of plotted samples: ",
  length(selected_samples)
)

message(
  "Number of standard samples: ",
  sum(sample_sources$sample_type == "standard")
)

message(
  "Number of Undetermined samples: ",
  sum(sample_sources$sample_type == "undetermined")
)

message(
  "PDF layout: ",
  number_of_columns,
  " columns x ",
  number_of_rows,
  " rows"
)

message(
  "Samples:"
)

message(
  paste(
    selected_samples,
    collapse = ", "
  )
)


# ==============================================================================
# 10. Check Space Ranger matrices and spatial folders
# ==============================================================================

matrix_files <- file.path(
  sample_sources$outs_dir,
  "raw_feature_bc_matrix.h5"
)


missing_matrix_samples <- sample_sources$sample_ID[
  !file.exists(matrix_files)
]


if (length(missing_matrix_samples) > 0) {

  missing_matrix_information <- sample_sources |>
    filter(
      sample_ID %in% missing_matrix_samples
    ) |>
    transmute(
      information = paste0(
        sample_ID,
        ": ",
        file.path(
          outs_dir,
          "raw_feature_bc_matrix.h5"
        )
      )
    ) |>
    pull(information)


  stop(
    "Missing raw_feature_bc_matrix.h5 for samples:\n",
    paste(
      missing_matrix_information,
      collapse = "\n"
    )
  )

}


spatial_directories <- file.path(
  sample_sources$outs_dir,
  "spatial"
)


missing_spatial_samples <- sample_sources$sample_ID[
  !dir.exists(spatial_directories)
]


if (length(missing_spatial_samples) > 0) {

  missing_spatial_information <- sample_sources |>
    filter(
      sample_ID %in% missing_spatial_samples
    ) |>
    transmute(
      information = paste0(
        sample_ID,
        ": ",
        file.path(
          outs_dir,
          "spatial"
        )
      )
    ) |>
    pull(information)


  stop(
    "Missing spatial directories for samples:\n",
    paste(
      missing_spatial_information,
      collapse = "\n"
    )
  )

}


message(
  "All required Space Ranger matrices and spatial folders were found."
)


# ==============================================================================
# 11. Read samples as separate Seurat objects
# ==============================================================================

samples_list <- setNames(
  map(
    selected_samples,
    function(sample_id) {

      message(
        "Reading sample: ",
        sample_id
      )


      sample_row <- sample_sources |>
        filter(
          sample_ID == sample_id
        )


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

        stop(
          "nCount_RNA is missing in Seurat metadata for sample: ",
          sample_id
        )

      }


      if (length(Images(seurat_object)) == 0) {

        stop(
          "No histological image was loaded for sample: ",
          sample_id
        )

      }


      if (ncol(seurat_object) == 0) {

        stop(
          "The Seurat object contains no spots for sample: ",
          sample_id
        )

      }


      seurat_object$plot_sample_ID <- sample_id
      seurat_object$metadata_sample_ID <- sample_row$metadata_sample_ID[[1]]
      seurat_object$sample_type <- sample_row$sample_type[[1]]
      seurat_object$group <- sample_row$group[[1]]
      seurat_object$sex <- sample_row$sex[[1]]


      seurat_object

    }
  ),
  selected_samples
)


if (length(samples_list) != length(selected_samples)) {

  stop(
    "The number of Seurat objects does not match ",
    "the number of selected samples."
  )

}


message(
  "All samples were loaded successfully."
)


# ==============================================================================
# 12. Calculate raw UMI statistics
# ==============================================================================

umi_summary <- imap_dfr(
  samples_list,
  function(seurat_object, sample_id) {

    umi <- as.numeric(
      seurat_object[[]][["nCount_RNA"]]
    )


    if (length(umi) == 0) {

      stop(
        "No UMI values found for sample: ",
        sample_id
      )

    }


    if (all(is.na(umi))) {

      stop(
        "All UMI values are NA for sample: ",
        sample_id
      )

    }


    tibble(

      sample_ID = sample_id,

      n_spots = length(umi),

      n_spots_with_UMI = sum(
        umi > 0,
        na.rm = TRUE
      ),

      percent_spots_with_UMI = 100 * sum(
        umi > 0,
        na.rm = TRUE
      ) / length(umi),

      total_UMI = sum(
        umi,
        na.rm = TRUE
      ),

      mean_UMI_per_spot = mean(
        umi,
        na.rm = TRUE
      ),

      median_UMI_per_spot = median(
        umi,
        na.rm = TRUE
      ),

      sd_UMI_per_spot = sd(
        umi,
        na.rm = TRUE
      ),

      min_UMI_per_spot = min(
        umi,
        na.rm = TRUE
      ),

      max_UMI_per_spot = max(
        umi,
        na.rm = TRUE
      )

    )

  }
) |>
  left_join(
    sample_sources |>
      select(
        sample_ID,
        metadata_sample_ID,
        sample_type,
        input_description,
        group,
        sex,
        outs_dir
      ),
    by = "sample_ID"
  ) |>
  select(
    sample_ID,
    metadata_sample_ID,
    sample_type,
    input_description,
    group,
    sex,
    everything()
  )


if (nrow(umi_summary) != length(samples_list)) {

  stop(
    "The number of UMI summary rows does not match ",
    "the number of loaded Seurat objects."
  )

}


if (anyDuplicated(umi_summary$sample_ID) > 0) {

  stop(
    "Duplicated sample IDs detected in umi_summary."
  )

}


message(
  "Raw UMI summary:"
)


print(
  umi_summary,
  n = Inf,
  width = Inf
)


# ==============================================================================
# 13. Save raw UMI summary
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

  stop(
    "The UMI summary file was not created: ",
    output_summary_tsv
  )

}


message(
  "UMI summary saved:"
)

message(
  normalizePath(
    output_summary_tsv,
    mustWork = TRUE
  )
)


# ==============================================================================
# 14. Helper function for displaying metadata
# ==============================================================================

format_metadata_value <- function(x) {

  if (
    length(x) == 0 ||
    is.na(x) ||
    trimws(as.character(x)) == ""
  ) {

    return("NA")

  }


  as.character(x)

}


# ==============================================================================
# 15. Create raw UMI maps
# ==============================================================================

umi_plots <- imap(
  samples_list,
  function(seurat_object, sample_id) {

    message(
      "Creating plot for sample: ",
      sample_id
    )


    summary_row_number <- match(
      sample_id,
      umi_summary$sample_ID
    )


    if (is.na(summary_row_number)) {

      stop(
        "No UMI summary found for sample: ",
        sample_id
      )

    }


    stats_row <- umi_summary[
      summary_row_number,
      ,
      drop = FALSE
    ]


    image_names <- Images(
      seurat_object
    )


    if (length(image_names) == 0) {

      stop(
        "No spatial image found for sample: ",
        sample_id
      )

    }


    image_name <- image_names[1]


    group_text <- format_metadata_value(
      stats_row$group
    )


    sex_text <- format_metadata_value(
      stats_row$sex
    )


    input_text <- format_metadata_value(
      stats_row$input_description
    )


    subtitle_text <- paste0(
      "Group: ",
      group_text,
      " | Sex: ",
      sex_text,
      "\nInput: ",
      input_text,
      "\nTotal UMI: ",
      format(
        stats_row$total_UMI,
        big.mark = " ",
        scientific = FALSE,
        trim = TRUE
      ),
      " | spots: ",
      format(
        stats_row$n_spots,
        big.mark = " ",
        scientific = FALSE,
        trim = TRUE
      ),
      "\nSpots with UMI: ",
      format(
        stats_row$n_spots_with_UMI,
        big.mark = " ",
        scientific = FALSE,
        trim = TRUE
      ),
      " (",
      format(
        round(
          stats_row$percent_spots_with_UMI,
          1
        ),
        nsmall = 1,
        trim = TRUE
      ),
      "%)",
      "\nMean: ",
      format(
        round(
          stats_row$mean_UMI_per_spot,
          1
        ),
        big.mark = " ",
        scientific = FALSE,
        trim = TRUE
      ),
      " | median: ",
      format(
        round(
          stats_row$median_UMI_per_spot,
          1
        ),
        big.mark = " ",
        scientific = FALSE,
        trim = TRUE
      ),
      " | max: ",
      format(
        round(
          stats_row$max_UMI_per_spot,
          0
        ),
        big.mark = " ",
        scientific = FALSE,
        trim = TRUE
      )
    )


    SpatialFeaturePlot(
      object = seurat_object,
      features = "nCount_RNA",
      images = image_name,
      image.alpha = 1,
      crop = FALSE,
      pt.size.factor = spot_size_factor
    ) +
      labs(
        title = paste0(
          "Sample: ",
          sample_id
        ),
        subtitle = subtitle_text,
        fill = "Raw UMI",
        colour = "Raw UMI"
      ) +
      guides(

        fill = guide_colorbar(
          title = "Raw UMI",
          title.position = "top",
          direction = "vertical",
          barheight = grid::unit(
            24,
            "mm"
          ),
          barwidth = grid::unit(
            3.5,
            "mm"
          )
        ),

        colour = guide_colorbar(
          title = "Raw UMI",
          title.position = "top",
          direction = "vertical",
          barheight = grid::unit(
            24,
            "mm"
          ),
          barwidth = grid::unit(
            3.5,
            "mm"
          )
        )

      ) +
      theme(

        text = element_text(
          family = "DejaVu Sans"
        ),

        plot.title = element_text(
          size = 14,
          face = "bold",
          hjust = 0.5
        ),

        plot.subtitle = element_text(
          size = 7.5,
          hjust = 0.5,
          lineheight = 1.08
        ),

        legend.position = "right",

        legend.direction = "vertical",

        legend.title = element_text(
          size = 8.5,
          face = "bold"
        ),

        legend.text = element_text(
          size = 7
        ),

        plot.margin = margin(
          t = 7,
          r = 7,
          b = 7,
          l = 7
        )

      )

  }
)


if (length(umi_plots) != length(selected_samples)) {

  stop(
    "The number of plots does not match ",
    "the number of selected samples."
  )

}


# ==============================================================================
# 16. Combine all maps
# ==============================================================================

all_umi_plot <- wrap_plots(
  umi_plots,
  ncol = number_of_columns,
  guides = "keep"
) +
  plot_annotation(

    title = paste0(
      "Raw UMI maps for maternal FMT Visium samples"
    ),

    subtitle = paste0(
      "20 standard samples + 4 Undetermined FASTQ analyses",
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

      text = element_text(
        family = "DejaVu Sans"
      ),

      plot.title = element_text(
        size = 22,
        face = "bold",
        hjust = 0.5
      ),

      plot.subtitle = element_text(
        size = 12,
        hjust = 0.5
      ),

      plot.margin = margin(
        t = 12,
        r = 10,
        b = 10,
        l = 10
      )

    )

  )


# ==============================================================================
# 17. Display combined plot
# ==============================================================================

print(
  all_umi_plot
)


# ==============================================================================
# 18. Select PDF device
# ==============================================================================

if (capabilities("cairo")) {

  pdf_device <- grDevices::cairo_pdf

  message(
    "Using cairo_pdf device."
  )

} else {

  pdf_device <- grDevices::pdf

  warning(
    "Cairo PDF is not available. ",
    "Using the standard PDF device."
  )

}


# ==============================================================================
# 19. Save combined plot as PDF
# ==============================================================================

message(
  "Saving PDF..."
)


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
# 20. Validate output
# ==============================================================================

if (!file.exists(output_pdf)) {

  stop(
    "The PDF file was not created: ",
    output_pdf
  )

}


output_pdf_size_mb <- file.info(output_pdf)$size / 1024^2


message(
  "PDF saved successfully:"
)

message(
  normalizePath(
    output_pdf,
    mustWork = TRUE
  )
)

message(
  "PDF size: ",
  round(
    output_pdf_size_mb,
    2
  ),
  " MB"
)


message(
  "UMI summary:"
)

message(
  normalizePath(
    output_summary_tsv,
    mustWork = TRUE
  )
)


message(
  "Analysis completed successfully."
)


# ==============================================================================
# End
# ==============================================================================

# ==============================================================================
# 12B. Compare spatial barcodes and raw UMI between samples
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)


# ------------------------------------------------------------------------------
# 1. Extract raw UMI for every spatial barcode
# ------------------------------------------------------------------------------

umi_by_barcode <- imap_dfr(
  samples_list,
  function(seurat_object, sample_id) {

    tibble(
      sample_ID = sample_id,
      barcode = colnames(seurat_object),
      raw_UMI = as.numeric(
        seurat_object[[]][["nCount_RNA"]]
      )
    )

  }
)


# ------------------------------------------------------------------------------
# 2. Create barcode sets
# ------------------------------------------------------------------------------

barcode_sets <- split(
  umi_by_barcode$barcode,
  umi_by_barcode$sample_ID
)


# ------------------------------------------------------------------------------
# 3. Generate unique sample pairs
# ------------------------------------------------------------------------------

sample_pairs <- combn(
  names(samples_list),
  2,
  simplify = FALSE
)


# ------------------------------------------------------------------------------
# 4. Compare barcode overlap and raw UMI correlation
# ------------------------------------------------------------------------------

barcode_comparison <- map_dfr(
  sample_pairs,
  function(pair) {

    sample_1 <- pair[[1]]
    sample_2 <- pair[[2]]

    barcodes_1 <- barcode_sets[[sample_1]]
    barcodes_2 <- barcode_sets[[sample_2]]

    shared_barcodes <- intersect(
      barcodes_1,
      barcodes_2
    )

    n_shared <- length(shared_barcodes)

    sample_1_data <- umi_by_barcode |>
      filter(
        sample_ID == sample_1,
        barcode %in% shared_barcodes
      ) |>
      select(
        barcode,
        raw_UMI_1 = raw_UMI
      )

    sample_2_data <- umi_by_barcode |>
      filter(
        sample_ID == sample_2,
        barcode %in% shared_barcodes
      ) |>
      select(
        barcode,
        raw_UMI_2 = raw_UMI
      )

    joined_data <- inner_join(
      sample_1_data,
      sample_2_data,
      by = "barcode"
    )

    # Correlation only when enough shared non-zero barcodes exist.
    correlation_data <- joined_data |>
      filter(
        raw_UMI_1 > 0 |
        raw_UMI_2 > 0
      )

    if (
      nrow(correlation_data) >= 100 &&
      sd(correlation_data$raw_UMI_1) > 0 &&
      sd(correlation_data$raw_UMI_2) > 0
    ) {

      pearson <- cor(
        correlation_data$raw_UMI_1,
        correlation_data$raw_UMI_2,
        method = "pearson",
        use = "complete.obs"
      )

      spearman <- cor(
        correlation_data$raw_UMI_1,
        correlation_data$raw_UMI_2,
        method = "spearman",
        use = "complete.obs"
      )

    } else {

      pearson <- NA_real_
      spearman <- NA_real_

    }

    tibble(
      sample_1 = sample_1,
      sample_2 = sample_2,

      n_barcodes_1 = length(barcodes_1),
      n_barcodes_2 = length(barcodes_2),

      n_shared_barcodes = n_shared,

      fraction_shared_sample_1 =
        n_shared / length(barcodes_1),

      fraction_shared_sample_2 =
        n_shared / length(barcodes_2),

      n_shared_nonzero_barcodes =
        nrow(correlation_data),

      pearson_raw_UMI = pearson,
      spearman_raw_UMI = spearman
    )

  }
) |>
  arrange(
    desc(spearman_raw_UMI),
    desc(n_shared_nonzero_barcodes)
  )


# ------------------------------------------------------------------------------
# 5. Display the most similar pairs
# ------------------------------------------------------------------------------
barcode_comparison %>% 
  filter(!grepl("undetermined", sample_1)) %>% 
  filter(!grepl("undetermined", sample_2)) %>% 
  select(!c(n_barcodes_1, n_barcodes_2, n_shared_nonzero_barcodes, fraction_shared_sample_1, fraction_shared_sample_2))  %>% 
  as.data.frame() 

# ==============================================================================
# Heatmaps of raw UMI correlations per spatial barcode
# ==============================================================================

