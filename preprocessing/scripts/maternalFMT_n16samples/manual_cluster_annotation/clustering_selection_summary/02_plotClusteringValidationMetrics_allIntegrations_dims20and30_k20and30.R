#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
options(scipen = 999)
required_packages <- c(
  "data.table",
  "ggplot2"
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
    paste0(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them with:\ninstall.packages(c(",
      paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}
library(data.table)
library(ggplot2)
message_line <- function(text = "", char = "=") {
  message(strrep(char, 72))
  if (nzchar(text)) {
    message(text)
    message(strrep(char, 72))
  }
}
read_validation_tsv <- function(file) {
  header_line <- readLines(
    file,
    n = 1L,
    warn = FALSE
  )
  if (length(header_line) != 1L) {
    stop(
      "Could not read the header from:\n",
      file,
      call. = FALSE
    )
  }
  header <- strsplit(
    header_line,
    "\t",
    fixed = TRUE
  )[[1L]]
  header[1L] <- sub(
    "^\\ufeff",
    "",
    header[1L]
  )
  first_two_lines <- readLines(
    file,
    n = 2L,
    warn = FALSE
  )
  if (length(first_two_lines) >= 2L) {
    first_data_nfields <- length(
      strsplit(
        first_two_lines[2L],
        "\t",
        fixed = TRUE
      )[[1L]]
    )
    if (length(header) != first_data_nfields) {
      merged_name <- "clusterSizeCVnClustersPresentAllSamples"
      if (
        merged_name %in% header &&
        length(header) + 1L == first_data_nfields
      ) {
        merged_index <- match(
          merged_name,
          header
        )
        header_before <- if (merged_index > 1L) {
          header[seq_len(merged_index - 1L)]
        } else {
          character(0)
        }
        header_after <- if (merged_index < length(header)) {
          header[(merged_index + 1L):length(header)]
        } else {
          character(0)
        }
        header <- c(
          header_before,
          "clusterSizeCV",
          "nClustersPresentAllSamples",
          header_after
        )
      } else {
        stop(
          "Header/data column-count mismatch in file:\n",
          file,
          "\nHeader fields: ",
          length(header),
          "\nFirst data-row fields: ",
          first_data_nfields,
          call. = FALSE
        )
      }
    }
  }
  result <- fread(
    file = file,
    sep = "\t",
    header = FALSE,
    skip = 1L,
    col.names = header,
    na.strings = c("NA", ""),
    quote = "",
    data.table = TRUE,
    showProgress = FALSE
  )
  if (ncol(result) != length(header)) {
    stop(
      "Unexpected number of columns after reading:\n",
      file,
      call. = FALSE
    )
  }
  result
}
assert_columns <- function(
  data,
  required,
  file
) {
  missing_columns <- setdiff(
    required,
    names(data)
  )
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      "\nFile: ",
      file,
      call. = FALSE
    )
  }
}
format_axis_values <- function(x) {
  formatC(
    x,
    format = "f",
    digits = 3L
  )
}
# -----------------------------------------------------------------------------
# User-editable parameters
# -----------------------------------------------------------------------------
project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"
input_root <- file.path(
  project_root,
  "results",
  dataset_name,
  "seurat_clustering_analysis",
  "logNormalize_vst",
  "hvg2000"
)
output_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "manual_cluster_annotation",
  "clustering_selection_summary",
  "figures"
)
integration_methods <- c(
  "cca",
  "rpca",
  "harmony",
  "fastMNN",
  "scVI",
  "noIntegration"
)
dims_values <- c(
  20L,
  30L
)
k_values <- c(
  20L,
  30L
)
clustering_algorithms <- c(
  "leiden",
  "louvain",
  "louvainRefined",
  "slm"
)
resolution_values <- seq(
  0.1,
  1.0,
  by = 0.1
)
prune_directory <- "prune0067"
# Colours can be changed here.
algorithm_colors <- c(
  leiden = "#E09F3E",
  louvain = "#335C67",
  louvainRefined = "#9E2A2B",
  slm = "#540B0E"
)
line_width <- 0.65
point_size <- 1.75
plot_width_in <- 16
plot_height_in <- 19
png_dpi <- 300
save_pdf <- TRUE
save_png <- TRUE
integration_labels <- c(
  cca = "CCA",
  rpca = "RPCA",
  harmony = "Harmony",
  fastMNN = "FastMNN",
  scVI = "scVI",
  noIntegration = "No integration"
)
parameter_panel_levels <- c(
  "dims20_k20",
  "dims20_k30",
  "dims30_k20",
  "dims30_k30"
)
parameter_panel_labels <- c(
  dims20_k20 = "dims 20 | k = 20",
  dims20_k30 = "dims 20 | k = 30",
  dims30_k20 = "dims 30 | k = 20",
  dims30_k30 = "dims 30 | k = 30"
)
# -----------------------------------------------------------------------------
# ASW annotation parameters
# The annotation is added only to the ASW plot.
# -----------------------------------------------------------------------------
asw_annotation_integration <- "cca"
asw_annotation_dims <- 20L
asw_annotation_k <- 20L
asw_annotation_algorithm <- "leiden"
asw_annotation_resolution <- 0.40
asw_annotation_colour <- "#000000"
asw_annotation_linewidth <- 0.75
asw_annotation_text_size <- 3.2
asw_annotation_label <- paste0(
  "Highest  ASW for\n",
  "Leiden clustering"
)
# -----------------------------------------------------------------------------
# Find and read all validation-summary files
# -----------------------------------------------------------------------------
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
message_line(
  "Collecting clustering-validation summaries"
)
message(
  "Input root:  ",
  input_root
)
message(
  "Output dir:  ",
  output_dir
)
if (!dir.exists(input_root)) {
  stop(
    "Input directory does not exist:\n",
    input_root,
    call. = FALSE
  )
}
validation_files <- list.files(
  path = input_root,
  pattern = "_allMethods_clusteringValidationSummary\\.tsv$",
  full.names = TRUE,
  recursive = TRUE
)
validation_files <- normalizePath(
  validation_files,
  winslash = "/",
  mustWork = TRUE
)
path_pattern <- paste0(
  "/02_",
  "(",
  paste(
    integration_methods,
    collapse = "|"
  ),
  ")",
  "/dims(20|30)",
  "/k(20|30)",
  "/",
  prune_directory,
  "/clustering_validation/tables/"
)
validation_files <- validation_files[
  grepl(
    path_pattern,
    validation_files,
    perl = TRUE
  )
]
if (length(validation_files) == 0L) {
  stop(
    paste0(
      "No allMethods_clusteringValidationSummary.tsv files ",
      "were found under:\n",
      input_root
    ),
    call. = FALSE
  )
}
file_map <- rbindlist(
  lapply(
    validation_files,
    function(file) {
      integration <- sub(
        paste0(
          "^.*?/02_([^/]+)/dims[0-9]+/k[0-9]+/",
          prune_directory,
          "/.*$"
        ),
        "\\1",
        file,
        perl = TRUE
      )
      dims <- as.integer(
        sub(
          "^.*?/dims([0-9]+)/k[0-9]+/.*$",
          "\\1",
          file,
          perl = TRUE
        )
      )
      k <- as.integer(
        sub(
          "^.*?/k([0-9]+)/.*$",
          "\\1",
          file,
          perl = TRUE
        )
      )
      data.table(
        integrationMethod = integration,
        dims = dims,
        k = k,
        file = file
      )
    }
  )
)
file_map <- file_map[
  integrationMethod %in% integration_methods &
    dims %in% dims_values &
    k %in% k_values
]
file_map[
  ,
  integrationOrder := match(
    integrationMethod,
    integration_methods
  )
]
setorder(
  file_map,
  integrationOrder,
  dims,
  k
)
file_map[
  ,
  integrationOrder := NULL
]
duplicate_groups <- file_map[
  ,
  .N,
  by = .(
    integrationMethod,
    dims,
    k
  )
][N != 1L]
if (nrow(duplicate_groups) > 0L) {
  stop(
    paste0(
      "Expected exactly one validation-summary file ",
      "per integration/dims/k group.\n"
    ),
    paste(
      capture.output(
        print(duplicate_groups)
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}
expected_groups <- CJ(
  integrationMethod = integration_methods,
  dims = dims_values,
  k = k_values,
  unique = TRUE
)
missing_groups <- expected_groups[
  !file_map,
  on = .(
    integrationMethod,
    dims,
    k
  )
]
if (nrow(missing_groups) > 0L) {
  stop(
    "Missing validation-summary files for these groups:\n",
    paste(
      capture.output(
        print(missing_groups)
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}
message(
  "Validation-summary files found: ",
  nrow(file_map),
  " / 24"
)
all_validation_data <- vector(
  "list",
  nrow(file_map)
)
required_columns <- c(
  "algorithm",
  "resolution",
  "nClusters",
  "transcriptomicASW_mean",
  "CHAOS_meanAcrossSamples",
  "PAS_meanAcrossSamples"
)
for (file_index in seq_len(nrow(file_map))) {
  current <- file_map[file_index]
  message(
    sprintf(
      "[%02d/%02d] %s, dims%d, k%d",
      file_index,
      nrow(file_map),
      current$integrationMethod,
      current$dims,
      current$k
    )
  )
  current_data <- read_validation_tsv(
    current$file
  )
  assert_columns(
    current_data,
    required_columns,
    current$file
  )
  current_data <- current_data[
    algorithm %in% clustering_algorithms &
      round(
        as.numeric(resolution),
        10L
      ) %in% round(
        resolution_values,
        10L
      )
  ]
  current_data[
    ,
    `:=`(
      integrationMethod = current$integrationMethod,
      dims = current$dims,
      k = current$k,
      resolution = round(
        as.numeric(resolution),
        1L
      ),
      nClusters = as.integer(nClusters),
      transcriptomicASW_mean = as.numeric(
        transcriptomicASW_mean
      ),
      CHAOS_meanAcrossSamples = as.numeric(
        CHAOS_meanAcrossSamples
      ),
      PAS_meanAcrossSamples = as.numeric(
        PAS_meanAcrossSamples
      )
    )
  ]
  current_data <- current_data[
    ,
    .(
      integrationMethod,
      dims,
      k,
      algorithm,
      resolution,
      nClusters,
      transcriptomicASW_mean,
      CHAOS_meanAcrossSamples,
      PAS_meanAcrossSamples
    )
  ]
  all_validation_data[[file_index]] <- current_data
}
plot_data <- rbindlist(
  all_validation_data,
  use.names = TRUE,
  fill = FALSE
)
duplicate_configurations <- plot_data[
  ,
  .N,
  by = .(
    integrationMethod,
    dims,
    k,
    algorithm,
    resolution
  )
][N != 1L]
if (nrow(duplicate_configurations) > 0L) {
  stop(
    "Duplicated clustering configurations were found:\n",
    paste(
      capture.output(
        print(duplicate_configurations)
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}
expected_configurations <- CJ(
  integrationMethod = integration_methods,
  dims = dims_values,
  k = k_values,
  algorithm = clustering_algorithms,
  resolution = round(
    resolution_values,
    1L
  ),
  unique = TRUE
)
missing_configurations <- expected_configurations[
  !plot_data,
  on = .(
    integrationMethod,
    dims,
    k,
    algorithm,
    resolution
  )
]
if (nrow(missing_configurations) > 0L) {
  stop(
    paste0(
      "Missing clustering configurations in ",
      "validation-summary files:\n"
    ),
    paste(
      capture.output(
        print(missing_configurations)
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}
if (nrow(plot_data) != 960L) {
  stop(
    "Expected 960 clustering configurations, but collected ",
    nrow(plot_data),
    ".",
    call. = FALSE
  )
}
message(
  "Clustering configurations collected: ",
  nrow(plot_data),
  " / 960"
)
# -----------------------------------------------------------------------------
# Prepare factor order and panel labels
# -----------------------------------------------------------------------------
plot_data[
  ,
  integrationMethod := factor(
    integrationMethod,
    levels = integration_methods
  )
]
plot_data[
  ,
  integrationLabel := factor(
    integration_labels[
      as.character(integrationMethod)
    ],
    levels = integration_labels[
      integration_methods
    ]
  )
]
plot_data[
  ,
  parameterPanel := paste0(
    "dims",
    dims,
    "_k",
    k
  )
]
plot_data[
  ,
  parameterPanel := factor(
    parameterPanel,
    levels = parameter_panel_levels
  )
]
plot_data[
  ,
  algorithm := factor(
    algorithm,
    levels = clustering_algorithms
  )
]
setorder(
  plot_data,
  integrationMethod,
  dims,
  k,
  algorithm,
  resolution
)
# -----------------------------------------------------------------------------
# Prepare the ASW annotation
# -----------------------------------------------------------------------------
selected_asw_row <- plot_data[
  as.character(integrationMethod) == asw_annotation_integration &
    dims == asw_annotation_dims &
    k == asw_annotation_k &
    as.character(algorithm) == asw_annotation_algorithm &
    abs(
      resolution - asw_annotation_resolution
    ) < 1e-8
]
if (nrow(selected_asw_row) != 1L) {
  stop(
    paste0(
      "Expected exactly one row for the ASW annotation, found ",
      nrow(selected_asw_row),
      "."
    ),
    call. = FALSE
  )
}
selected_asw_value <- selected_asw_row[["transcriptomicASW_mean"]][1L]
global_asw_maximum <- max(
  plot_data$transcriptomicASW_mean,
  na.rm = TRUE
)
if (
  abs(
    selected_asw_value - global_asw_maximum
  ) > 1e-8
) {
  stop(
    paste0(
      "The selected CCA/dims20/k20/Leiden/resolution0.40 ",
      "configuration is not the configuration with the highest ",
      "transcriptomic ASW.\n",
      "Selected ASW: ",
      selected_asw_value,
      "\nGlobal maximum ASW: ",
      global_asw_maximum
    ),
    call. = FALSE
  )
}
asw_annotation_line <- selected_asw_row[
  ,
  .(
    integrationLabel,
    parameterPanel,
    xintercept = asw_annotation_resolution
  )
]
asw_annotation_text <- selected_asw_row[
  ,
  .(
    integrationLabel,
    parameterPanel,
    x = asw_annotation_resolution + 0.025,
    y = selected_asw_value + 0.012,
    label = asw_annotation_label
  )
]
message(
  "ASW annotation value: ",
  formatC(
    selected_asw_value,
    format = "f",
    digits = 5L
  )
)
# -----------------------------------------------------------------------------
# Plotting function used only within this script
# -----------------------------------------------------------------------------
make_metric_plot <- function(
  data,
  metric_column,
  plot_title,
  plot_subtitle,
  y_axis_title,
  y_expand_top = 0.08
) {
  current_data <- copy(data)
  current_data[
    ,
    metricValue := get(metric_column)
  ]
  ggplot(
    current_data,
    aes(
      x = resolution,
      y = metricValue,
      colour = algorithm,
      group = algorithm
    )
  ) +
    geom_line(
      linewidth = line_width,
      na.rm = TRUE
    ) +
    geom_point(
      size = point_size,
      na.rm = TRUE
    ) +
    facet_grid(
      rows = vars(integrationLabel),
      cols = vars(parameterPanel),
      labeller = labeller(
        parameterPanel = parameter_panel_labels
      ),
      drop = FALSE,
      switch = "y"
    ) +
    scale_colour_manual(
      values = algorithm_colors,
      breaks = clustering_algorithms,
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = resolution_values,
      limits = c(0.1, 1.0),
      expand = expansion(
        mult = c(0.02, 0.02)
      )
    ) +
    scale_y_continuous(
      labels = format_axis_values,
      expand = expansion(
        mult = c(
          0.05,
          y_expand_top
        )
      )
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Clustering resolution",
      y = y_axis_title,
      colour = "Clustering algorithm"
    ) +
    theme_bw(
      base_size = 10.5
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 17,
        hjust = 0.5,
        margin = margin(
          b = 5
        )
      ),
      plot.subtitle = element_text(
        size = 10.5,
        hjust = 0.5,
        margin = margin(
          b = 9
        )
      ),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_text(
        face = "bold"
      ),
      legend.key.width = grid::unit(
        1.25,
        "lines"
      ),
      strip.background = element_rect(
        fill = "grey96",
        colour = "black",
        linewidth = 0.45
      ),
      strip.text.x = element_text(
        face = "bold",
        size = 10
      ),
      strip.text.y.left = element_text(
        face = "bold",
        size = 10,
        angle = 90
      ),
      strip.placement = "outside",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.spacing.x = grid::unit(
        0.75,
        "lines"
      ),
      panel.spacing.y = grid::unit(
        0.8,
        "lines"
      ),
      axis.title = element_text(
        face = "bold"
      ),
      axis.text.x = element_text(
        size = 7.5,
        angle = 0,
        vjust = 0.5
      ),
      axis.text.y = element_text(
        size = 8
      ),
      plot.margin = margin(
        10,
        12,
        10,
        12
      )
    )
}
save_metric_plot <- function(
  plot,
  file_stem
) {
  saved_files <- character(0)
  if (isTRUE(save_pdf)) {
    pdf_file <- file.path(
      output_dir,
      paste0(
        file_stem,
        ".pdf"
      )
    )
    ggsave(
      filename = pdf_file,
      plot = plot,
      width = plot_width_in,
      height = plot_height_in,
      units = "in",
      limitsize = FALSE
    )
    saved_files <- c(
      saved_files,
      pdf_file
    )
  }
  if (isTRUE(save_png)) {
    png_file <- file.path(
      output_dir,
      paste0(
        file_stem,
        ".png"
      )
    )
    ggsave(
      filename = png_file,
      plot = plot,
      width = plot_width_in,
      height = plot_height_in,
      units = "in",
      dpi = png_dpi,
      limitsize = FALSE
    )
    saved_files <- c(
      saved_files,
      png_file
    )
  }
  saved_files
}
# -----------------------------------------------------------------------------
# Generate the ASW plot
# -----------------------------------------------------------------------------
message_line(
  "Generating plots"
)
asw_plot <- make_metric_plot(
  data = plot_data,
  metric_column = "transcriptomicASW_mean",
  plot_title = "Transcriptomic ASW across clustering configurations",
  plot_subtitle = paste0(
    "Rows: integration method | columns: dimensions and k | ",
    "colours: clustering algorithm | ",
    "higher values indicate better separation"
  ),
  y_axis_title = "Transcriptomic ASW",
  y_expand_top = 0.18
)
# Add the vertical line and text only to the ASW plot.
asw_plot <- asw_plot +
  geom_vline(
    data = asw_annotation_line,
    mapping = aes(
      xintercept = xintercept
    ),
    inherit.aes = FALSE,
    colour = asw_annotation_colour,
    linetype = "dashed",
    linewidth = asw_annotation_linewidth
  ) +
  geom_text(
    data = asw_annotation_text,
    mapping = aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    colour = asw_annotation_colour,
    size = asw_annotation_text_size,
    fontface = "bold",
    hjust = 0,
    vjust = 0.5,
    lineheight = 0.95
  )
# -----------------------------------------------------------------------------
# Generate CHAOS and PAS plots without annotations
# -----------------------------------------------------------------------------
chaos_plot <- make_metric_plot(
  data = plot_data,
  metric_column = "CHAOS_meanAcrossSamples",
  plot_title = "CHAOS across clustering configurations",
  plot_subtitle = paste0(
    "Rows: integration method | columns: dimensions and k | ",
    "colours: clustering algorithm | ",
    "lower values indicate better spatial coherence"
  ),
  y_axis_title = "CHAOS, mean across samples"
)
pas_plot <- make_metric_plot(
  data = plot_data,
  metric_column = "PAS_meanAcrossSamples",
  plot_title = "PAS across clustering configurations",
  plot_subtitle = paste0(
    "Rows: integration method | columns: dimensions and k | ",
    "colours: clustering algorithm | ",
    "lower values indicate better spatial coherence"
  ),
  y_axis_title = "PAS, mean across samples"
)
# -----------------------------------------------------------------------------
# Save plots
# -----------------------------------------------------------------------------
saved_files <- c(
  save_metric_plot(
    asw_plot,
    paste0(
      "01_",
      dataset_name,
      "_transcriptomicASW_allIntegrations_dims20and30_k20and30"
    )
  ),
  save_metric_plot(
    chaos_plot,
    paste0(
      "02_",
      dataset_name,
      "_CHAOS_allIntegrations_dims20and30_k20and30"
    )
  ),
  save_metric_plot(
    pas_plot,
    paste0(
      "03_",
      dataset_name,
      "_PAS_allIntegrations_dims20and30_k20and30"
    )
  )
)
message_line(
  "Completed"
)
message(
  "Files saved:"
)
for (file in saved_files) {
  message(
    "  ",
    file
  )
}