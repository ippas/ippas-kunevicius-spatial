# ==============================================================================
# 01_raw_UMI_maps_and_correlations_all_samples.R
#
# Purpose:
# - read all 20 Visium samples as separate Seurat objects
# - calculate raw UMI statistics separately for every sample
# - plot raw UMI over histological images
# - add sample ID, FMT donor group, sex and UMI summary
# - calculate Pearson and Spearman correlations between samples
#
# No normalization, integration or clustering.
# ==============================================================================


# ==============================================================================
# 1. Packages
# ==============================================================================

suppressPackageStartupMessages({

  library(Seurat)

  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)

  library(ggplot2)
  library(patchwork)

})


# ==============================================================================
# 2. Paths
# ==============================================================================

path_to_data <- file.path(
  "data",
  "spacerangerCount_withoutJSON_GRCm39-2024-A_2026-03-28"
)

metadata_file <- file.path(
  "data",
  "metadata_autismFMT.tsv"
)


# ==============================================================================
# 3. Read metadata
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


metadata_samples <- metadata_autismFMT |>
  transmute(
    sample_ID = as.character(sample_ID),
    group = as.character(fmt_donor_group),
    sex = as.character(sex)
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

message(
  "Number of samples: ",
  length(selected_samples)
)


# ==============================================================================
# 4. Check Space Ranger files
# ==============================================================================

missing_matrix_files <- selected_samples[
  !file.exists(
    file.path(
      path_to_data,
      selected_samples,
      "outs",
      "raw_feature_bc_matrix.h5"
    )
  )
]

if (length(missing_matrix_files) > 0) {

  stop(
    "Missing raw matrices for samples: ",
    paste(
      missing_matrix_files,
      collapse = ", "
    )
  )

}


# ==============================================================================
# 5. Read all samples as separate Seurat objects
# ==============================================================================

samples_list <- setNames(
  map(
    selected_samples,
    function(sample_id) {

      message("Reading sample: ", sample_id)

      Load10X_Spatial(
        data.dir = file.path(
          path_to_data,
          sample_id,
          "outs"
        ),
        filename = "raw_feature_bc_matrix.h5",
        assay = "RNA",
        slice = sample_id,
        filter.matrix = FALSE
      )

    }
  ),
  selected_samples
)


# ==============================================================================
# 6. Calculate UMI statistics separately for every sample
# ==============================================================================

umi_summary <- imap_dfr(
  samples_list,
  function(seurat_object, sample_id) {

    umi <- as.numeric(
      seurat_object[[]][["nCount_RNA"]]
    )

    tibble(
      sample_ID = sample_id,
      n_spots = length(umi),
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
      )
    )

  }
) |>
  left_join(
    metadata_samples,
    by = "sample_ID"
  )


if (nrow(umi_summary) != length(samples_list)) {

  stop(
    "The number of UMI summary rows does not match ",
    "the number of Seurat objects."
  )

}

if (anyDuplicated(umi_summary$sample_ID) > 0) {

  stop(
    "Duplicated sample IDs detected in umi_summary."
  )

}


print(
  umi_summary,
  n = Inf
)


# ==============================================================================
# 7. Plot raw UMI over histological images
# ==============================================================================

umi_plots <- imap(
  samples_list,
  function(seurat_object, sample_id) {

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


    image_name <- Images(
      seurat_object
    )[1]


    subtitle_text <- paste0(
      "Group: ",
      stats_row$group,
      " | Sex: ",
      stats_row$sex,
      "\nUMI summary: total = ",
      format(
        stats_row$total_UMI,
        big.mark = " ",
        scientific = FALSE,
        trim = TRUE
      ),
      " | spots = ",
      stats_row$n_spots,
      "\nMean/spot = ",
      round(
        stats_row$mean_UMI_per_spot,
        1
      ),
      " | median/spot = ",
      round(
        stats_row$median_UMI_per_spot,
        1
      ),
      " | SD/spot = ",
      round(
        stats_row$sd_UMI_per_spot,
        1
      )
    )


    SpatialFeaturePlot(
      object = seurat_object,
      features = "nCount_RNA",
      images = image_name,
      image.alpha = 1,
      crop = FALSE,
      pt.size.factor = 1.6
    ) +
      labs(
        title = paste0(
          "Sample ID: ",
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
            32,
            "mm"
          ),
          barwidth = grid::unit(
            4,
            "mm"
          )
        ),
        colour = guide_colorbar(
          title = "Raw UMI",
          title.position = "top",
          direction = "vertical",
          barheight = grid::unit(
            32,
            "mm"
          ),
          barwidth = grid::unit(
            4,
            "mm"
          )
        )
      ) +
      theme(
        text = element_text(
          family = "DejaVu Sans"
        ),
        plot.title = element_text(
          size = 18,
          face = "bold",
          hjust = 0.5
        ),
        plot.subtitle = element_text(
          size = 11,
          hjust = 0.5,
          lineheight = 1.15
        ),
        legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_text(
          size = 11,
          face = "bold"
        ),
        legend.text = element_text(
          size = 9
        ),
        plot.margin = margin(
          t = 10,
          r = 10,
          b = 10,
          l = 10
        )
      )

  }
)


all_umi_plot <- wrap_plots(
  umi_plots,
  ncol = 6,
  guides = "keep"
)

all_umi_plot


# ==============================================================================
# 8. Extract raw UMI and array coordinates
# ==============================================================================

spot_data <- imap_dfr(
  samples_list,
  function(seurat_object, sample_id) {

    image_name <- Images(
      seurat_object
    )[1]


    coordinates <- slot(
      seurat_object[[image_name]],
      "coordinates"
    ) |>
      rownames_to_column(
        "barcode"
      )


    if (
      all(
        c(
          "row",
          "col"
        ) %in% colnames(coordinates)
      )
    ) {

      coordinates <- coordinates |>
        transmute(
          barcode = barcode,
          array_row = row,
          array_col = col
        )

    } else if (
      all(
        c(
          "array_row",
          "array_col"
        ) %in% colnames(coordinates)
      )
    ) {

      coordinates <- coordinates |>
        transmute(
          barcode = barcode,
          array_row = array_row,
          array_col = array_col
        )

    } else {

      stop(
        "Array coordinates not found for sample: ",
        sample_id
      )

    }


    umi_data <- tibble(
      barcode = colnames(
        seurat_object
      ),
      UMI = as.numeric(
        seurat_object[[]][["nCount_RNA"]]
      )
    )


    coordinates |>
      inner_join(
        umi_data,
        by = "barcode"
      ) |>
      mutate(
        sample_ID = sample_id,
        .before = 1
      )

  }
)


print(
  spot_data |>
    count(sample_ID),
  n = Inf
)


# ==============================================================================
# 9. Check duplicated array positions
# ==============================================================================

duplicate_coordinates <- spot_data |>
  count(
    sample_ID,
    array_row,
    array_col
  ) |>
  filter(
    n > 1
  )

if (nrow(duplicate_coordinates) > 0) {

  stop(
    "Duplicated array coordinates detected."
  )

}


# ==============================================================================
# 10. Prepare array-position by sample UMI matrix
# ==============================================================================

umi_wide <- spot_data |>
  select(
    sample_ID,
    array_row,
    array_col,
    UMI
  ) |>
  pivot_wider(
    id_cols = c(
      array_row,
      array_col
    ),
    names_from = sample_ID,
    values_from = UMI
  )


umi_matrix <- umi_wide |>
  select(
    all_of(
      selected_samples
    )
  ) |>
  as.matrix()


# ==============================================================================
# 11. Pearson and Spearman correlations
# ==============================================================================

umi_cor_pearson <- cor(
  umi_matrix,
  method = "pearson",
  use = "pairwise.complete.obs"
)

umi_cor_spearman <- cor(
  umi_matrix,
  method = "spearman",
  use = "pairwise.complete.obs"
)


message("Pearson correlation matrix:")

print(
  round(
    umi_cor_pearson,
    3
  )
)


message("Spearman correlation matrix:")

print(
  round(
    umi_cor_spearman,
    3
  )
)


# ==============================================================================
# 12. Pairwise correlation table
# ==============================================================================

pearson_long <- as.data.frame(
  as.table(
    umi_cor_pearson
  ),
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  transmute(
    sample_1 = as.character(
      Var1
    ),
    sample_2 = as.character(
      Var2
    ),
    pearson_r = as.numeric(
      Freq
    )
  )


spearman_long <- as.data.frame(
  as.table(
    umi_cor_spearman
  ),
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  transmute(
    sample_1 = as.character(
      Var1
    ),
    sample_2 = as.character(
      Var2
    ),
    spearman_rho = as.numeric(
      Freq
    )
  )


umi_correlations <- pearson_long |>
  left_join(
    spearman_long,
    by = c(
      "sample_1",
      "sample_2"
    )
  ) |>
  filter(
    sample_1 < sample_2
  ) |>
  arrange(
    desc(
      pearson_r
    )
  )


print(
  umi_correlations,
  n = Inf
)