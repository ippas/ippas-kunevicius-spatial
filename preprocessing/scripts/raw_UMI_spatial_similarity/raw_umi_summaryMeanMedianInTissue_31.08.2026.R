#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ==============================================================================
# 1. Input table
# ==============================================================================

umi_summary <- data.frame(
  Sample = c(
    "1_1F",
    "1_1Fd",
    "1_1M",
    "2_1M",
    "2_1Md",
    "3_1F",
    "3_1M",
    "5_1M",
    "5_3F",
    "12_1M",
    "13_1M",
    "13_1F",
    "15_1F",
    "18_1F",
    "18_1M",
    "23_1F"
  ),
  Slide_ID = c(
    "V12N14-292",
    "V12N14-292",
    "V12N14-298",
    "V12N14-298",
    "V12N14-297",
    "V12N14-292",
    "V12N14-298",
    "V12N14-298",
    "V12N14-324",
    "V12N14-297",
    "V12N14-297",
    "V12N14-292",
    "V12N14-324",
    "V12N14-324",
    "V12N14-297",
    "V12N14-324"
  ),
  Donor_group = c(
    "Neurotypical",
    "Neurotypical",
    "Neurotypical",
    "Neurotypical",
    "Neurotypical",
    "Neurotypical",
    "Neurotypical",
    "ASD",
    "ASD",
    "ASD",
    "ASD",
    "ASD",
    "ASD",
    "ASD",
    "ASD",
    "ASD"
  ),
  Sex = c(
    "Female",
    "Female",
    "Male",
    "Male",
    "Male",
    "Female",
    "Male",
    "Male",
    "Female",
    "Male",
    "Male",
    "Female",
    "Female",
    "Female",
    "Male",
    "Female"
  ),
  Mean_UMI_spot = c(
    4010,
    3498,
    5319,
    8356,
    3926,
    4660,
    6991,
    3903,
    10662,
    2813,
    3642,
    1809,
    9676,
    8621,
    2413,
    11156
  ),
  SD_UMI_spot = c(
    4742,
    3146,
    3183,
    4305,
    2042,
    3651,
    3413,
    2083,
    4136,
    1295,
    2395,
    1839,
    5436,
    3740,
    1252,
    4801
  ),
  Median_UMI_spot = c(
    2482,
    2425,
    5010,
    7820,
    3629,
    3559,
    6495,
    3779,
    10054,
    2592,
    2797,
    1182,
    8794,
    8074,
    2180,
    10699
  ),
  IQR_UMI_spot = c(
    3246,
    4240,
    4056,
    5862,
    2732,
    5293,
    4328,
    2678,
    5214,
    1628,
    3092,
    1706,
    6856,
    4930,
    1438,
    6537
  ),
  Q25_UMI_spot = c(
    1344,
    1015,
    3031,
    5136,
    2380,
    1762,
    4625,
    2479,
    7742,
    1883,
    1937,
    656,
    5702,
    5880,
    1549,
    7565
  ),
  Q75_UMI_spot = c(
    4591,
    5254,
    7087,
    10998,
    5112,
    7055,
    8952,
    5157,
    12956,
    3511,
    5029,
    2361,
    12558,
    10810,
    2987,
    14102
  ),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 2. Prepare groups and factor order
# ==============================================================================

umi_summary <- umi_summary |>
  mutate(
    Group = paste(
      Donor_group,
      Sex,
      sep = " | "
    ),
    Group = factor(
      Group,
      levels = c(
        "Neurotypical | Female",
        "Neurotypical | Male",
        "ASD | Female",
        "ASD | Male"
      )
    ),
    Slide_ID = factor(
      Slide_ID,
      levels = c(
        "V12N14-292",
        "V12N14-297",
        "V12N14-298",
        "V12N14-324"
      )
    )
  )

# ==============================================================================
# 3. Summary: Mean UMI/spot by donor group and sex
# ==============================================================================

mean_group_data <- umi_summary |>
  group_by(Group) |>
  summarise(
    bar_value = mean(
      Mean_UMI_spot,
      na.rm = TRUE
    ),
    sd_value = sd(
      Mean_UMI_spot,
      na.rm = TRUE
    ),
    n = n(),
    .groups = "drop"
  )

# ==============================================================================
# 4. Summary: Median UMI/spot by donor group and sex
# ==============================================================================

median_group_data <- umi_summary |>
  group_by(Group) |>
  summarise(
    bar_value = mean(
      Median_UMI_spot,
      na.rm = TRUE
    ),
    sd_value = sd(
      Median_UMI_spot,
      na.rm = TRUE
    ),
    n = n(),
    .groups = "drop"
  )

# ==============================================================================
# 5. Summary: Mean UMI/spot by Slide ID
# ==============================================================================

mean_slide_data <- umi_summary |>
  group_by(Slide_ID) |>
  summarise(
    bar_value = mean(
      Mean_UMI_spot,
      na.rm = TRUE
    ),
    sd_value = sd(
      Mean_UMI_spot,
      na.rm = TRUE
    ),
    n = n(),
    .groups = "drop"
  )

# ==============================================================================
# 6. Summary: Median UMI/spot by Slide ID
# ==============================================================================

median_slide_data <- umi_summary |>
  group_by(Slide_ID) |>
  summarise(
    bar_value = mean(
      Median_UMI_spot,
      na.rm = TRUE
    ),
    sd_value = sd(
      Median_UMI_spot,
      na.rm = TRUE
    ),
    n = n(),
    .groups = "drop"
  )

# ==============================================================================
# 7. Plot 1: Mean UMI/spot by donor group and sex
# ==============================================================================

plot_mean_group <- ggplot(
  mean_group_data,
  aes(
    x = Group,
    y = bar_value
  )
) +
  geom_col(
    width = 0.7,
    fill = "grey75",
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = bar_value - sd_value,
      ymax = bar_value + sd_value
    ),
    width = 0.2,
    linewidth = 0.5
  ) +
  geom_point(
    data = umi_summary,
    aes(
      x = Group,
      y = Mean_UMI_spot
    ),
    inherit.aes = FALSE,
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    size = 2
  ) +
  labs(
    title = "Mean UMI per spot",
    x = NULL,
    y = "Mean UMI/spot"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )

# ==============================================================================
# 8. Plot 2: Median UMI/spot by donor group and sex
# ==============================================================================

plot_median_group <- ggplot(
  median_group_data,
  aes(
    x = Group,
    y = bar_value
  )
) +
  geom_col(
    width = 0.7,
    fill = "grey75",
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = bar_value - sd_value,
      ymax = bar_value + sd_value
    ),
    width = 0.2,
    linewidth = 0.5
  ) +
  geom_point(
    data = umi_summary,
    aes(
      x = Group,
      y = Median_UMI_spot
    ),
    inherit.aes = FALSE,
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    size = 2
  ) +
  labs(
    title = "Median UMI per spot",
    x = NULL,
    y = "Median UMI/spot"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )

# ==============================================================================
# 9. Plot 3: Mean UMI/spot by Slide ID
# ==============================================================================

plot_mean_slide <- ggplot(
  mean_slide_data,
  aes(
    x = Slide_ID,
    y = bar_value
  )
) +
  geom_col(
    width = 0.7,
    fill = "grey75",
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = bar_value - sd_value,
      ymax = bar_value + sd_value
    ),
    width = 0.2,
    linewidth = 0.5
  ) +
  geom_point(
    data = umi_summary,
    aes(
      x = Slide_ID,
      y = Mean_UMI_spot
    ),
    inherit.aes = FALSE,
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    size = 2
  ) +
  labs(
    title = "Mean UMI per spot by Slide ID",
    x = "Slide ID",
    y = "Mean UMI/spot"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )

# ==============================================================================
# 10. Plot 4: Median UMI/spot by Slide ID
# ==============================================================================

plot_median_slide <- ggplot(
  median_slide_data,
  aes(
    x = Slide_ID,
    y = bar_value
  )
) +
  geom_col(
    width = 0.7,
    fill = "grey75",
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = bar_value - sd_value,
      ymax = bar_value + sd_value
    ),
    width = 0.2,
    linewidth = 0.5
  ) +
  geom_point(
    data = umi_summary,
    aes(
      x = Slide_ID,
      y = Median_UMI_spot
    ),
    inherit.aes = FALSE,
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    size = 2
  ) +
  labs(
    title = "Median UMI per spot by Slide ID",
    x = "Slide ID",
    y = "Median UMI/spot"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )

# ==============================================================================
# 11. Display plots
# ==============================================================================

print(plot_mean_group)

print(plot_median_group)

print(plot_mean_slide)

print(plot_median_slide)

# ==============================================================================
# 12. Print summary tables
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("Mean UMI/spot by donor group and sex\n")
cat("============================================================\n")

print(mean_group_data)

cat("\n")
cat("============================================================\n")
cat("Median UMI/spot by donor group and sex\n")
cat("============================================================\n")

print(median_group_data)

cat("\n")
cat("============================================================\n")
cat("Mean UMI/spot by Slide ID\n")
cat("============================================================\n")

print(mean_slide_data)

cat("\n")
cat("============================================================\n")
cat("Median UMI/spot by Slide ID\n")
cat("============================================================\n")

print(median_slide_data)