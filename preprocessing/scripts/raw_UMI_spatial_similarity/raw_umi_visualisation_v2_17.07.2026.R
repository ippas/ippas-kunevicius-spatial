# ==============================================================================
# 01_raw_UMI_maps_all_samples.R
#
# Purpose:
# - read all 20 Visium samples as separate Seurat objects
# - calculate raw UMI statistics separately for every sample
# - plot raw UMI over histological images
# - add sample ID, FMT donor group, sex and UMI summary
# - save all plots into one PDF
#
# Output layout:
# - 4 columns
# - 5 rows for 20 samples
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
  "01_raw_UMI_maps_all_samples_4_columns.pdf"
)


# ==============================================================================
# 3. Plot and PDF settings
# ==============================================================================

number_of_columns <- 4


# Large custom PDF page suitable for 20 samples.
pdf_width_inches <- 20
pdf_height_inches <- 28


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
    "Space Ranger data directory does not exist: ",
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
  "Space Ranger data directory: ",
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
# 6. Prepare sample metadata
# ==============================================================================

metadata_samples <- metadata_autismFMT |>
  transmute(
    sample_ID = trimws(
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
    !is.na(sample_ID),
    sample_ID != ""
  ) |>
  distinct(
    sample_ID,
    .keep_all = TRUE
  )


selected_samples <- metadata_samples$sample_ID


if (length(selected_samples) == 0) {

  stop(
    "No samples were found in the metadata file."
  )

}


if (anyDuplicated(selected_samples) > 0) {

  stop(
    "Duplicated sample IDs were detected after metadata preparation."
  )

}


message(
  "Number of samples: ",
  length(selected_samples)
)


number_of_rows <- ceiling(
  length(selected_samples) / number_of_columns
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
# 7. Check Space Ranger matrices and spatial folders
# ==============================================================================

matrix_files <- file.path(
  path_to_data,
  selected_samples,
  "outs",
  "raw_feature_bc_matrix.h5"
)


missing_matrix_samples <- selected_samples[
  !file.exists(matrix_files)
]


if (length(missing_matrix_samples) > 0) {

  stop(
    "Missing raw_feature_bc_matrix.h5 for samples: ",
    paste(
      missing_matrix_samples,
      collapse = ", "
    )
  )

}


spatial_directories <- file.path(
  path_to_data,
  selected_samples,
  "outs",
  "spatial"
)


missing_spatial_samples <- selected_samples[
  !dir.exists(spatial_directories)
]


if (length(missing_spatial_samples) > 0) {

  stop(
    "Missing spatial directory for samples: ",
    paste(
      missing_spatial_samples,
      collapse = ", "
    )
  )

}


message(
  "All required Space Ranger matrices and spatial folders were found."
)


# ==============================================================================
# 8. Read samples as separate Seurat objects
# ==============================================================================

samples_list <- setNames(
  map(
    selected_samples,
    function(sample_id) {

      message(
        "Reading sample: ",
        sample_id
      )


      sample_outs_dir <- file.path(
        path_to_data,
        sample_id,
        "outs"
      )


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
# 9. Calculate raw UMI statistics
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


    tibble(
      sample_ID = sample_id,

      n_spots = length(umi),

      n_spots_with_UMI = sum(
        umi > 0,
        na.rm = TRUE
      ),

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

      max_UMI_per_spot = max(
        umi,
        na.rm = TRUE
      )

    )

  }
) |>
  left_join(
    metadata_samples,
    by = "sample_ID"
  ) |>
  select(
    sample_ID,
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
  n = Inf
)


# ==============================================================================
# 10. Helper function for displaying metadata
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
# 11. Create raw UMI maps
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


    subtitle_text <- paste0(
      "Group: ",
      group_text,
      " | Sex: ",
      sex_text,
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
          size = 15,
          face = "bold",
          hjust = 0.5
        ),

        plot.subtitle = element_text(
          size = 8.5,
          hjust = 0.5,
          lineheight = 1.1
        ),

        legend.position = "right",
        legend.direction = "vertical",

        legend.title = element_text(
          size = 9,
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
# 12. Combine all maps
# ==============================================================================

all_umi_plot <- wrap_plots(
  umi_plots,
  ncol = number_of_columns,
  guides = "keep"
) +
  plot_annotation(
    title = "Raw UMI maps for maternal FMT Visium samples",

    subtitle = paste0(
      "Number of samples: ",
      length(selected_samples),
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
# 13. Display combined plot
# ==============================================================================

print(
  all_umi_plot
)


# ==============================================================================
# 14. Select PDF device
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
# 15. Save combined plot as PDF
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
# 16. Validate output
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


# ==============================================================================
# End
# ==============================================================================