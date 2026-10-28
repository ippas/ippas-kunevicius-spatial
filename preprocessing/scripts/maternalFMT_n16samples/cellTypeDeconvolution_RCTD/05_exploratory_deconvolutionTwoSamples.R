library(data.table)
library(ggplot2)
library(writexl)

# ==============================================================================
# PATH
# ==============================================================================

dir <- "/home/mateusz/projects/ippas-kunevicius-spatial/results/maternalFMT_n16samples/cellTypeDeconvolution_RCTD/logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040/00_test_twoSamples"

files <- file.path(
  dir,
  "perSpot",
  c(
    "3_1F_RCTD_full_cellTypeProportions_perSpot.tsv",
    "15_1F_RCTD_full_cellTypeProportions_perSpot.tsv"
  )
)

plot_dir <- file.path(dir, "plots", "barplots_perCluster")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# READ
# ==============================================================================

dat <- rbindlist(lapply(files, fread))

metadata_columns <- c(
  "spot", "sample_ID", "image", "x", "y", "cluster"
)

cell_types <- setdiff(names(dat), metadata_columns)

# ==============================================================================
# NORMALIZE RCTD WEIGHTS PER SPOT
# ==============================================================================

total_weight <- rowSums(dat[, ..cell_types])

dat[
  ,
  (cell_types) := lapply(.SD, function(x) x / total_weight),
  .SDcols = cell_types
]

print(summary(rowSums(dat[, ..cell_types])))

# ==============================================================================
# LONG FORMAT
# ==============================================================================

long <- melt(
  dat,
  id.vars = metadata_columns,
  measure.vars = cell_types,
  variable.name = "cell_type",
  value.name = "proportion"
)

long[, proportion := proportion * 100]

# ==============================================================================
# WHOLE SECTION
# ==============================================================================

section <- long[
  ,
  .(
    mean = mean(proportion),
    SD = sd(proportion)
  ),
  by = .(cell_type, sample_ID)
]

section_mean <- dcast(
  section,
  cell_type ~ sample_ID,
  value.var = "mean"
)

section_sd <- dcast(
  section,
  cell_type ~ sample_ID,
  value.var = "SD"
)

setnames(
  section_mean,
  setdiff(names(section_mean), "cell_type"),
  paste0(setdiff(names(section_mean), "cell_type"), "_mean")
)

setnames(
  section_sd,
  setdiff(names(section_sd), "cell_type"),
  paste0(setdiff(names(section_sd), "cell_type"), "_SD")
)

section_table <- merge(
  section_mean,
  section_sd,
  by = "cell_type"
)

setcolorder(
  section_table,
  c(
    "cell_type",
    "3_1F_mean", "3_1F_SD",
    "15_1F_mean", "15_1F_SD"
  )
)

# ==============================================================================
# PER CLUSTER
# ==============================================================================

cluster_summary <- long[
  ,
  .(
    mean = mean(proportion),
    SD = sd(proportion)
  ),
  by = .(cluster, cell_type, sample_ID)
]

cluster_mean <- dcast(
  cluster_summary,
  cluster + cell_type ~ sample_ID,
  value.var = "mean"
)

cluster_sd <- dcast(
  cluster_summary,
  cluster + cell_type ~ sample_ID,
  value.var = "SD"
)

setnames(
  cluster_mean,
  setdiff(names(cluster_mean), c("cluster", "cell_type")),
  paste0(
    setdiff(names(cluster_mean), c("cluster", "cell_type")),
    "_mean"
  )
)

setnames(
  cluster_sd,
  setdiff(names(cluster_sd), c("cluster", "cell_type")),
  paste0(
    setdiff(names(cluster_sd), c("cluster", "cell_type")),
    "_SD"
  )
)

cluster_table <- merge(
  cluster_mean,
  cluster_sd,
  by = c("cluster", "cell_type")
)

setcolorder(
  cluster_table,
  c(
    "cluster", "cell_type",
    "3_1F_mean", "3_1F_SD",
    "15_1F_mean", "15_1F_SD"
  )
)

# ==============================================================================
# SAVE XLSX
# ==============================================================================

write_xlsx(
  list(
    whole_section = as.data.frame(section_table),
    per_cluster = as.data.frame(cluster_table)
  ),
  path = file.path(
    dir,
    "RCTD_cellTypeProportions_3_1F_vs_15_1F.xlsx"
  )
)

# ==============================================================================
# BARPLOTS PER CELL TYPE
# ==============================================================================

cluster_summary[, cluster := factor(
  cluster,
  levels = sort(unique(as.numeric(as.character(cluster))))
)]

for (ct in unique(cluster_summary$cell_type)) {

  plot_data <- cluster_summary[cell_type == ct]

  p <- ggplot(
    plot_data,
    aes(
      x = cluster,
      y = mean,
      fill = sample_ID
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7
    ) +
    geom_errorbar(
      aes(
        ymin = pmax(0, mean - SD),
        ymax = mean + SD
      ),
      position = position_dodge(width = 0.8),
      width = 0.2
    ) +
    labs(
      title = ct,
      x = "Cluster",
      y = "Estimated proportion (%)",
      fill = "Sample"
    ) +
    theme_classic(base_size = 12)

  ggsave(
    file.path(
      plot_dir,
      paste0(ct, "_proportion_perCluster.pdf")
    ),
    p,
    width = 9,
    height = 5
  )
}


# ==============================================================================
# BARPLOT - WHOLE SECTION
# ==============================================================================

p <- ggplot(
  section,
  aes(
    x = cell_type,
    y = mean,
    fill = sample_ID
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean - SD),
      ymax = mean + SD
    ),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  labs(
    x = "Cell type",
    y = "Estimated proportion (%)",
    fill = "Sample"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  file.path(
    dir,
    "plots",
    "RCTD_cellTypeProportions_wholeSection_3_1F_vs_15_1F.pdf"
  ),
  p,
  width = 11,
  height = 6
)
