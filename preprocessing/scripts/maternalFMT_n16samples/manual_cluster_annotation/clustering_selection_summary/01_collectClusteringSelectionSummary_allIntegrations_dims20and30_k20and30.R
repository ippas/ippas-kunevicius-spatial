#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(scipen = 999)
options(openxlsx.maxWidth = 40)

required_packages <- c("data.table", "openxlsx")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
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
library(openxlsx)

message_line <- function(text = "", char = "=") {
  message(strrep(char, 72))
  if (nzchar(text)) {
    message(text)
    message(strrep(char, 72))
  }
}

find_single_file <- function(directory, pattern, exclude_pattern = NULL) {
  files <- list.files(
    path = directory,
    pattern = pattern,
    full.names = TRUE,
    recursive = FALSE
  )
  if (!is.null(exclude_pattern)) {
    files <- files[!grepl(exclude_pattern, basename(files), perl = TRUE)]
  }
  files <- sort(files)
  list(
    path = if (length(files) == 1L) files else NA_character_,
    match_count = length(files),
    matches = files,
    status = if (length(files) == 1L) {
      "OK"
    } else if (length(files) == 0L) {
      "MISSING"
    } else {
      "MULTIPLE_MATCHES"
    }
  )
}

read_tsv_strict <- function(file) {
  header_line <- readLines(file, n = 1L, warn = FALSE)
  if (length(header_line) != 1L) {
    stop("Could not read the header from: ", file, call. = FALSE)
  }
  header <- strsplit(header_line, "\t", fixed = TRUE)[[1L]]
  header[1L] <- sub("^\\ufeff", "", header[1L])

  first_two_lines <- readLines(file, n = 2L, warn = FALSE)
  if (length(first_two_lines) >= 2L) {
    first_data_nfields <- length(
      strsplit(first_two_lines[2L], "\t", fixed = TRUE)[[1L]]
    )

    if (length(header) != first_data_nfields) {
      merged_name <- "clusterSizeCVnClustersPresentAllSamples"
      if (
        merged_name %in% header &&
        length(header) + 1L == first_data_nfields
      ) {
        merged_index <- match(merged_name, header)
        header <- c(
          header[seq_len(merged_index - 1L)],
                "nClustersPresentAllSamples",
          header[(merged_index + 1L):length(header)]
        )
        warning(
          "Repaired a merged validation-summary header in: ",
          basename(file),
          call. = FALSE
        )
      } else {
        stop(
          "Header/data column-count mismatch in file:\n",
          file,
          "\nHeader fields: ", length(header),
          "\nFirst data-row fields: ", first_data_nfields,
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

assert_columns <- function(data, required, label, file) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing, collapse = ", "),
      "\nFile: ", file,
      call. = FALSE
    )
  }
}

format_resolution <- function(x) sprintf("%.2f", as.numeric(x))

format_up_to_5_decimals <- function(x) {
  formatted <- formatC(
    round(as.numeric(x), 5L),
    format = "f",
    digits = 5L
  )
  sub("\\.?0+$", "", formatted)
}

make_configuration_id <- function(integration, dims, k, algorithm, resolution) {
  paste0(
    integration,
    "_dims", dims,
    "_k", k,
    "_", algorithm,
    "_res", format_resolution(resolution)
  )
}

safe_row_mean <- function(x, y) {
  result <- rowMeans(cbind(x, y), na.rm = TRUE)
  result[is.nan(result)] <- NA_real_
  result
}

safe_row_sd <- function(x, y) {
  vapply(
    seq_along(x),
    function(index) {
      values <- c(x[index], y[index])
      values <- values[!is.na(values)]
      if (length(values) >= 2L) stats::sd(values) else NA_real_
    },
    numeric(1)
  )
}

extract_sample_id <- function(spot) {
  sample_id <- sub("_[ACGTN]+-[0-9]+$", "", spot, perl = TRUE)
  not_parsed <- which(sample_id == spot)
  if (length(not_parsed) > 0L) {
    examples <- paste(head(spot[not_parsed], 5L), collapse = ", ")
    stop(
      "Could not extract sample_ID from some spot names. Examples: ",
      examples,
      call. = FALSE
    )
  }
  sample_id
}

write_table_sheet <- function(
  workbook,
  sheet_name,
  data,
  table_name,
  freeze_first_column = TRUE
) {
  addWorksheet(workbook, sheet_name, gridLines = FALSE)
  writeDataTable(
    workbook,
    sheet = sheet_name,
    x = data,
    startRow = 1L,
    startCol = 1L,
    tableName = table_name,
    tableStyle = "TableStyleMedium2",
    withFilter = TRUE
  )
  freezePane(
    workbook,
    sheet = sheet_name,
    firstRow = TRUE,
    firstCol = freeze_first_column
  )
  setColWidths(
    workbook,
    sheet = sheet_name,
    cols = seq_len(ncol(data)),
    widths = "auto"
  )
  setRowHeights(workbook, sheet = sheet_name, rows = 1L, heights = 32)
}

apply_numeric_format <- function(workbook, sheet, data, columns, num_fmt) {
  columns <- intersect(columns, names(data))
  if (length(columns) == 0L || nrow(data) == 0L) {
    return(invisible(NULL))
  }
  style <- createStyle(numFmt = num_fmt)
  addStyle(
    workbook,
    sheet = sheet,
    style = style,
    rows = 2L:(nrow(data) + 1L),
    cols = match(columns, names(data)),
    gridExpand = TRUE,
    stack = TRUE
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
  "clustering_selection_summary"
)

output_xlsx <- file.path(
  output_dir,
  paste0(
    "01_", dataset_name,
    "_clusteringSelectionSummary_allIntegrations_dims20and30_k20and30.xlsx"
  )
)

integration_methods <- c(
  "cca",
  "rpca",
  "harmony",
  "fastMNN",
  "scVI",
  "noIntegration"
)
dims_values <- c(20L, 30L)
k_values <- c(20L, 30L)
clustering_algorithms <- c(
  "louvain",
  "louvainRefined",
  "slm",
  "leiden"
)
resolution_values <- seq(0.1, 1.0, by = 0.1)
prune_directory <- "prune0067"

selected_integration <- "cca"
selected_dims <- 20L
selected_k <- 20L
selected_algorithm <- "leiden"
selected_resolution <- 0.40

sample_order <- c(
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
  "13_1F",
  "13_1M",
  "15_1F",
  "18_1F",
  "18_1M",
  "23_1F"
)

expected_configuration_groups <-
  length(integration_methods) * length(dims_values) * length(k_values)
expected_clustering_rows <-
  expected_configuration_groups *
  length(clustering_algorithms) *
  length(resolution_values)

# -----------------------------------------------------------------------------
# Collect input tables
# -----------------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message_line("Collecting clustering-selection data")
message("Input root:  ", input_root)
message("Output file: ", output_xlsx)

all_configuration_list <- list()
all_sample_metrics_list <- list()
all_cluster_sizes_list <- list()
all_ari_list <- list()
manifest_list <- list()

configuration_group_index <- 0L

for (integration in integration_methods) {
  for (dims in dims_values) {
    for (k in k_values) {
      configuration_group_index <- configuration_group_index + 1L
      group_id <- paste0(integration, "_dims", dims, "_k", k)
      message(
        sprintf(
          "[%d/%d] %s",
          configuration_group_index,
          expected_configuration_groups,
          group_id
        )
      )

      configuration_dir <- file.path(
        input_root,
        paste0("02_", integration),
        paste0("dims", dims),
        paste0("k", k),
        prune_directory
      )
      tables_dir <- file.path(configuration_dir, "tables")
      validation_tables_dir <- file.path(
        configuration_dir,
        "clustering_validation",
        "tables"
      )

      file_results <- list(
        clusteringSummary = find_single_file(
          tables_dir,
          "_res010to100by010_multiClusteringAndUmap_clusteringSummary\\.tsv$",
          exclude_pattern = "testNoIntegration|_dev_"
        ),
        validationSummary = find_single_file(
          validation_tables_dir,
          "_allMethods_clusteringValidationSummary\\.tsv$"
        ),
        spatialValidationBySample = find_single_file(
          validation_tables_dir,
          "_allMethods_spatialValidationBySample\\.tsv$"
        ),
        clusteringStabilityARI = find_single_file(
          validation_tables_dir,
          "_clusteringStabilityARI\\.tsv$",
          exclude_pattern = paste0(
            "_(",
            paste(clustering_algorithms, collapse = "|"),
            ")_clusteringStabilityARI\\.tsv$"
          )
        ),
        clusterAssignments = find_single_file(
          tables_dir,
          "_res010to100by010_multiClusteringAndUmap_clusterAssignments\\.tsv\\.gz$",
          exclude_pattern = "testNoIntegration|_dev_"
        )
      )

      group_manifest <- rbindlist(
        lapply(
          names(file_results),
          function(input_type) {
            result <- file_results[[input_type]]
            data.table(
              configurationGroupID = group_id,
              integrationMethod = integration,
              dims = dims,
              k = k,
              pruneDirectory = prune_directory,
              inputType = input_type,
              matchedFile = result$path,
              fileExists = !is.na(result$path) && file.exists(result$path),
              matchCount = result$match_count,
              status = result$status,
              allMatches = paste(result$matches, collapse = " | ")
            )
          }
        ),
        use.names = TRUE,
        fill = TRUE
      )
      manifest_list[[length(manifest_list) + 1L]] <- group_manifest

      if (any(group_manifest$status != "OK")) {
        warning(
          "Skipping incomplete configuration group: ",
          group_id,
          call. = FALSE
        )
        next
      }

      clustering_summary_file <- file_results$clusteringSummary$path
      validation_summary_file <- file_results$validationSummary$path
      sample_metrics_file <- file_results$spatialValidationBySample$path
      ari_file <- file_results$clusteringStabilityARI$path
      assignments_file <- file_results$clusterAssignments$path

      clustering_summary <- read_tsv_strict(clustering_summary_file)
      validation_summary <- read_tsv_strict(validation_summary_file)
      sample_metrics <- read_tsv_strict(sample_metrics_file)
      ari_details <- read_tsv_strict(ari_file)
      assignments <- fread(
        assignments_file,
        sep = "\t",
        header = TRUE,
        na.strings = c("NA", ""),
        quote = "",
        data.table = TRUE,
        showProgress = FALSE
      )

      assert_columns(
        clustering_summary,
        c(
          "algorithm",
          "resolution",
          "clusterColumn",
          "nSpots",
          "nClusters"
        ),
        "Clustering summary",
        clustering_summary_file
      )
      assert_columns(
        validation_summary,
        c(
          "algorithm",
          "resolution",
          "clusterColumn",
          "nSpots",
          "transcriptomicASW_mean",
          "transcriptomicASW_sd",
          "transcriptomicASW_fractionNegative",
          "CHAOS_meanAcrossSamples",
          "CHAOS_sdAcrossSamples",
          "PAS_meanAcrossSamples",
          "PAS_sdAcrossSamples",
          "PAS_meanPercentAcrossSamples",
          "nClusters",
          "normalization",
          "nHVG",
          "reduction",
          "dimsEnd",
          "kParam",
          "pruneSNN",
          "ARI_vsPreviousResolution",
          "ARI_vsNextResolution"
        ),
        "Clustering validation summary",
        validation_summary_file
      )
      assert_columns(
        sample_metrics,
        c(
          "algorithm",
          "resolution",
          "sample_ID",
          "clusterColumn",
          "nSpots",
          "nClustersPresent",
          "CHAOS",
          "CHAOS_nValidSpots",
          "CHAOS_nExcludedSingletonClusterSpots",
          "PAS_percent",
          "PAS_nAbnormalSpots",
          "PAS_k"
        ),
        "Spatial validation by sample",
        sample_metrics_file
      )
      assert_columns(
        ari_details,
        c(
          "algorithm1",
          "resolution1",
          "algorithm2",
          "resolution2",
          "ARI"
        ),
        "ARI details",
        ari_file
      )
      assert_columns(
        assignments,
        c("spot"),
        "Cluster assignments",
        assignments_file
      )

      expected_rows_per_group <-
        length(clustering_algorithms) * length(resolution_values)
      if (nrow(validation_summary) != expected_rows_per_group) {
        warning(
          group_id,
          ": expected ", expected_rows_per_group,
          " validation-summary rows, found ", nrow(validation_summary),
          call. = FALSE
        )
      }

      clustering_cluster_counts <- clustering_summary[
        ,
        .(
          algorithm,
          resolution = as.numeric(resolution),
          clusterColumn,
          nClustersFromClusteringSummary = as.integer(nClusters)
        )
      ]
      validation_summary[, resolution := as.numeric(resolution)]

      configuration_summary <- merge(
        validation_summary,
        clustering_cluster_counts,
        by = c("algorithm", "resolution", "clusterColumn"),
        all.x = TRUE,
        sort = FALSE
      )

      mismatch <- configuration_summary[
        !is.na(nClustersFromClusteringSummary) &
          as.integer(nClusters) != nClustersFromClusteringSummary
      ]
      if (nrow(mismatch) > 0L) {
        stop(
          "Number-of-clusters mismatch between source tables for ",
          group_id,
          call. = FALSE
        )
      }

      configuration_summary[
        ,
        `:=`(
          configurationID = make_configuration_id(
            integration,
            dims,
            k,
            algorithm,
            resolution
          ),
          integrationMethod = integration,
          dims = dims,
          k = k,
          ARI_adjacentMean = safe_row_mean(
            as.numeric(ARI_vsPreviousResolution),
            as.numeric(ARI_vsNextResolution)
          ),
          ARI_adjacentSD = safe_row_sd(
            as.numeric(ARI_vsPreviousResolution),
            as.numeric(ARI_vsNextResolution)
          )
        )
      ]

      wanted_configuration_columns <- c(
        "configurationID",
        "integrationMethod",
        "normalization",
        "nHVG",
        "reduction",
        "dims",
        "k",
        "pruneSNN",
        "algorithm",
        "resolution",
        "clusterColumn",
        "nSpots",
        "nClusters",
        "transcriptomicASW_mean",
        "transcriptomicASW_sd",
        "transcriptomicASW_fractionNegative",
        "CHAOS_meanAcrossSamples",
        "CHAOS_sdAcrossSamples",
        "PAS_meanAcrossSamples",
        "PAS_sdAcrossSamples",
        "PAS_meanPercentAcrossSamples",
        "ARI_vsPreviousResolution",
        "ARI_vsNextResolution",
        "ARI_adjacentMean",
        "ARI_adjacentSD"
      )
      configuration_summary <- configuration_summary[
        ,
        ..wanted_configuration_columns
      ]
      all_configuration_list[[length(all_configuration_list) + 1L]] <-
        configuration_summary

      sample_metrics[, resolution := as.numeric(resolution)]
      sample_metrics[
        ,
        `:=`(
          configurationID = make_configuration_id(
            integration,
            dims,
            k,
            algorithm,
            resolution
          ),
          integrationMethod = integration,
          dims = dims,
          k = k
        )
      ]

      wanted_sample_columns <- c(
        "configurationID",
        "integrationMethod",
        "dims",
        "k",
        "algorithm",
        "resolution",
        "clusterColumn",
        "sample_ID",
        "nSpots",
        "nClustersPresent",
        "CHAOS",
        "CHAOS_nValidSpots",
        "CHAOS_nExcludedSingletonClusterSpots",
        "PAS_percent",
        "PAS_nAbnormalSpots",
        "PAS_k"
      )
      sample_metrics <- sample_metrics[, ..wanted_sample_columns]
      all_sample_metrics_list[[length(all_sample_metrics_list) + 1L]] <-
        sample_metrics

      ari_details[
        ,
        `:=`(
          integrationMethod = integration,
          dims = dims,
          k = k,
          configurationID1 = make_configuration_id(
            integration,
            dims,
            k,
            algorithm1,
            resolution1
          ),
          configurationID2 = make_configuration_id(
            integration,
            dims,
            k,
            algorithm2,
            resolution2
          )
        )
      ]
      ari_details <- ari_details[
        ,
        .(
          integrationMethod,
          dims,
          k,
          configurationID1,
          algorithm1,
          resolution1,
          configurationID2,
          algorithm2,
          resolution2,
          ARI
        )
      ]
      all_ari_list[[length(all_ari_list) + 1L]] <- ari_details

      assignments[, sample_ID := extract_sample_id(spot)]

      observed_sample_spots <- assignments[, .(assignmentNSpots = .N), by = sample_ID]
      source_sample_spots <- unique(sample_metrics[, .(sample_ID, nSpots)])
      if (source_sample_spots[, anyDuplicated(sample_ID)] > 0L) {
        stop(
          "Per-sample validation table contains inconsistent nSpots values for ",
          group_id,
          call. = FALSE
        )
      }
      sample_spot_check <- merge(
        source_sample_spots,
        observed_sample_spots,
        by = "sample_ID",
        all = TRUE
      )
      if (
        anyNA(sample_spot_check$nSpots) ||
        anyNA(sample_spot_check$assignmentNSpots) ||
        any(sample_spot_check$nSpots != sample_spot_check$assignmentNSpots)
      ) {
        stop(
          "Per-sample spot counts differ between validation and assignment tables for ",
          group_id,
          call. = FALSE
        )
      }
      if (uniqueN(assignments$sample_ID) != length(sample_order)) {
        warning(
          group_id,
          ": expected ", length(sample_order),
          " samples, found ", uniqueN(assignments$sample_ID),
          call. = FALSE
        )
      }

      expected_cluster_columns <- unique(validation_summary$clusterColumn)
      missing_assignment_columns <- setdiff(
        expected_cluster_columns,
        names(assignments)
      )
      if (length(missing_assignment_columns) > 0L) {
        stop(
          "Cluster-assignment table is missing columns for ",
          group_id,
          ": ",
          paste(missing_assignment_columns, collapse = ", "),
          call. = FALSE
        )
      }

      group_samples <- unique(assignments$sample_ID)
      group_cluster_sizes <- vector(
        mode = "list",
        length = length(expected_cluster_columns)
      )

      for (column_index in seq_along(expected_cluster_columns)) {
        cluster_column <- expected_cluster_columns[column_index]
        validation_row <- validation_summary[clusterColumn == cluster_column]
        if (nrow(validation_row) != 1L) {
          stop(
            "Expected one validation row for cluster column: ",
            cluster_column,
            call. = FALSE
          )
        }

        algorithm <- validation_row$algorithm[1L]
        resolution <- as.numeric(validation_row$resolution[1L])
        cluster_values <- as.integer(assignments[[cluster_column]])
        if (anyNA(cluster_values)) {
          stop(
            "Missing cluster assignments detected in column: ",
            cluster_column,
            call. = FALSE
          )
        }
        cluster_ids <- sort(unique(cluster_values))
        expected_n_clusters <- as.integer(validation_row$nClusters[1L])
        if (length(cluster_ids) != expected_n_clusters) {
          stop(
            "Number of unique cluster IDs differs from nClusters for ",
            group_id,
            ", column ", cluster_column,
            ". Unique IDs: ", length(cluster_ids),
            "; nClusters: ", expected_n_clusters,
            call. = FALSE
          )
        }

        observed_counts <- assignments[
          ,
          .(nSpotsInCluster = .N),
          by = .(
            sample_ID,
            clusterID = as.integer(get(cluster_column))
          )
        ]
        complete_grid <- CJ(
          sample_ID = group_samples,
          clusterID = cluster_ids,
          unique = TRUE
        )
        cluster_sizes <- merge(
          complete_grid,
          observed_counts,
          by = c("sample_ID", "clusterID"),
          all.x = TRUE,
          sort = FALSE
        )
        cluster_sizes[is.na(nSpotsInCluster), nSpotsInCluster := 0L]
        cluster_sizes[
          ,
          `:=`(
            configurationID = make_configuration_id(
              integration,
              dims,
              k,
              algorithm,
              resolution
            ),
            integrationMethod = integration,
            dims = dims,
            k = k,
            algorithm = algorithm,
            resolution = resolution,
            clusterColumn = cluster_column
          )
        ]
        setcolorder(
          cluster_sizes,
          c(
            "configurationID",
            "integrationMethod",
            "dims",
            "k",
            "algorithm",
            "resolution",
            "clusterColumn",
            "sample_ID",
            "clusterID",
            "nSpotsInCluster"
          )
        )
        group_cluster_sizes[[column_index]] <- cluster_sizes
      }

      all_cluster_sizes_list[[length(all_cluster_sizes_list) + 1L]] <-
        rbindlist(group_cluster_sizes, use.names = TRUE, fill = TRUE)

      rm(
        clustering_summary,
        validation_summary,
        sample_metrics,
        ari_details,
        assignments,
        group_cluster_sizes
      )
      invisible(gc(verbose = FALSE))
    }
  }
}

input_status <- rbindlist(manifest_list, use.names = TRUE, fill = TRUE)
all_configurations <- rbindlist(
  all_configuration_list,
  use.names = TRUE,
  fill = TRUE
)
all_metrics_by_sample <- rbindlist(
  all_sample_metrics_list,
  use.names = TRUE,
  fill = TRUE
)
all_cluster_sizes <- rbindlist(
  all_cluster_sizes_list,
  use.names = TRUE,
  fill = TRUE
)
all_ari_details <- rbindlist(all_ari_list, use.names = TRUE, fill = TRUE)

if (nrow(all_configurations) == 0L) {
  stop(
    "No complete clustering configuration groups were found. Check the terminal warnings for missing or ambiguous input files.",
    call. = FALSE
  )
}

all_configurations[
  ,
  integrationMethod := factor(
    integrationMethod,
    levels = integration_methods,
    ordered = TRUE
  )
]
all_configurations[
  ,
  algorithm := factor(
    algorithm,
    levels = clustering_algorithms,
    ordered = TRUE
  )
]
setorder(
  all_configurations,
  integrationMethod,
  dims,
  k,
  algorithm,
  resolution
)
all_configurations[, integrationMethod := as.character(integrationMethod)]
all_configurations[, algorithm := as.character(algorithm)]

all_metrics_by_sample[
  ,
  integrationMethod := factor(
    integrationMethod,
    levels = integration_methods,
    ordered = TRUE
  )
]
all_metrics_by_sample[
  ,
  algorithm := factor(
    algorithm,
    levels = clustering_algorithms,
    ordered = TRUE
  )
]
all_metrics_by_sample[
  ,
  sampleOrder := match(sample_ID, sample_order)
]
setorder(
  all_metrics_by_sample,
  integrationMethod,
  dims,
  k,
  algorithm,
  resolution,
  sampleOrder,
  sample_ID
)
all_metrics_by_sample[, sampleOrder := NULL]
all_metrics_by_sample[, integrationMethod := as.character(integrationMethod)]
all_metrics_by_sample[, algorithm := as.character(algorithm)]

all_cluster_sizes[
  ,
  integrationMethod := factor(
    integrationMethod,
    levels = integration_methods,
    ordered = TRUE
  )
]
all_cluster_sizes[
  ,
  algorithm := factor(
    algorithm,
    levels = clustering_algorithms,
    ordered = TRUE
  )
]
all_cluster_sizes[, sampleOrder := match(sample_ID, sample_order)]
setorder(
  all_cluster_sizes,
  integrationMethod,
  dims,
  k,
  algorithm,
  resolution,
  sampleOrder,
  sample_ID,
  clusterID
)
all_cluster_sizes[, sampleOrder := NULL]
all_cluster_sizes[, integrationMethod := as.character(integrationMethod)]
all_cluster_sizes[, algorithm := as.character(algorithm)]

setorder(
  all_ari_details,
  integrationMethod,
  dims,
  k,
  algorithm1,
  resolution1,
  algorithm2,
  resolution2
)
setorder(input_status, integrationMethod, dims, k, inputType)

selected_summary <- all_configurations[
  integrationMethod == selected_integration &
    dims == selected_dims &
    k == selected_k &
    algorithm == selected_algorithm &
    abs(resolution - selected_resolution) < 1e-8
]
selected_metrics <- all_metrics_by_sample[
  integrationMethod == selected_integration &
    dims == selected_dims &
    k == selected_k &
    algorithm == selected_algorithm &
    abs(resolution - selected_resolution) < 1e-8
]
selected_cluster_sizes <- all_cluster_sizes[
  integrationMethod == selected_integration &
    dims == selected_dims &
    k == selected_k &
    algorithm == selected_algorithm &
    abs(resolution - selected_resolution) < 1e-8
]

if (nrow(selected_summary) != 1L) {
  stop(
    "Expected exactly one selected clustering row, found ",
    nrow(selected_summary),
    call. = FALSE
  )
}
if (nrow(selected_metrics) == 0L) {
  stop("No per-sample metrics found for the selected clustering.", call. = FALSE)
}
if (nrow(selected_cluster_sizes) == 0L) {
  stop("No per-sample cluster sizes found for the selected clustering.", call. = FALSE)
}

if (nrow(all_configurations) != expected_clustering_rows) {
  warning(
    "Expected ", expected_clustering_rows,
    " rows in AllConfigurations, collected ", nrow(all_configurations),
    ". Check the terminal warnings for missing or ambiguous input files.",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# Prepare selected-clustering tables
# -----------------------------------------------------------------------------

selected_parameters <- data.table(
  parameter = c(
    "dataset",
    "integration method",
    "normalization",
    "number of HVGs",
    "reduction",
    "dims",
    "k",
    "prune SNN",
    "clustering algorithm",
    "resolution",
    "cluster column"
  ),
  value = c(
    dataset_name,
    selected_summary$integrationMethod,
    selected_summary$normalization,
    selected_summary$nHVG,
    selected_summary$reduction,
    selected_summary$dims,
    selected_summary$k,
    format_up_to_5_decimals(selected_summary$pruneSNN),
    selected_summary$algorithm,
    format_resolution(selected_summary$resolution),
    selected_summary$clusterColumn
  )
)

selected_global_metrics <- data.table(
  metric = c(
    "number of spots",
    "number of clusters",
    "transcriptomic ASW mean",
    "transcriptomic ASW SD",
    "fraction of negative transcriptomic ASW",
    "CHAOS mean across samples",
    "CHAOS SD across samples",
    "PAS mean across samples",
    "PAS SD across samples",
    "PAS mean percent across samples",
    "ARI vs previous resolution",
    "ARI vs next resolution",
    "ARI adjacent mean",
    "ARI adjacent SD"
  ),
  value = unlist(
    selected_summary[
      ,
      .(
        nSpots,
        nClusters,
        transcriptomicASW_mean,
        transcriptomicASW_sd,
        transcriptomicASW_fractionNegative,
        CHAOS_meanAcrossSamples,
        CHAOS_sdAcrossSamples,
        PAS_meanAcrossSamples,
        PAS_sdAcrossSamples,
        PAS_meanPercentAcrossSamples,
        ARI_vsPreviousResolution,
        ARI_vsNextResolution,
        ARI_adjacentMean,
        ARI_adjacentSD
      )
    ],
    use.names = FALSE
  )
)

selected_global_metrics_split <- ceiling(nrow(selected_global_metrics) / 2L)
selected_global_metrics_left <- selected_global_metrics[
  seq_len(selected_global_metrics_split)
]
selected_global_metrics_right <- selected_global_metrics[
  (selected_global_metrics_split + 1L):nrow(selected_global_metrics)
]

selected_metrics_export <- copy(selected_metrics)
selected_cluster_sizes_export <- copy(selected_cluster_sizes)
all_configurations_export <- copy(all_configurations)
all_metrics_export <- copy(all_metrics_by_sample)
all_cluster_sizes_export <- copy(all_cluster_sizes)

# -----------------------------------------------------------------------------
# Write Excel workbook
# -----------------------------------------------------------------------------

message_line("Writing Excel workbook")

workbook <- createWorkbook(creator = "Mateusz Zięba")

write_table_sheet(
  workbook,
  "01_AllConfigurations",
  all_configurations_export,
  "AllConfigurationsTable",
  freeze_first_column = TRUE
)

addWorksheet(workbook, "02_SelectedClustering", gridLines = FALSE)
mergeCells(workbook, "02_SelectedClustering", cols = 1:7, rows = 1)
writeData(
  workbook,
  "02_SelectedClustering",
  paste0(
    "Selected clustering: CCA, dims20, k20, Leiden, resolution 0.40"
  ),
  startRow = 1L,
  startCol = 1L
)

title_style <- createStyle(
  fontSize = 16,
  textDecoration = "bold",
  fontColour = "#FFFFFF",
  fgFill = "#1F4E78",
  halign = "left",
  valign = "center"
)
section_style <- createStyle(
  fontSize = 12,
  textDecoration = "bold",
  fontColour = "#FFFFFF",
  fgFill = "#5B9BD5",
  halign = "left",
  valign = "center"
)

addStyle(
  workbook,
  "02_SelectedClustering",
  title_style,
  rows = 1L,
  cols = 1:7,
  gridExpand = TRUE,
  stack = TRUE
)
setRowHeights(workbook, "02_SelectedClustering", rows = 1L, heights = 28)

writeData(
  workbook,
  "02_SelectedClustering",
  "Parameters",
  startRow = 3L,
  startCol = 1L
)
addStyle(
  workbook,
  "02_SelectedClustering",
  section_style,
  rows = 3L,
  cols = 1:2,
  gridExpand = TRUE,
  stack = TRUE
)
writeDataTable(
  workbook,
  "02_SelectedClustering",
  selected_parameters,
  startRow = 4L,
  startCol = 1L,
  tableName = "SelectedParametersTable",
  tableStyle = "TableStyleMedium2",
  withFilter = FALSE
)

mergeCells(
  workbook,
  "02_SelectedClustering",
  cols = 4:7,
  rows = 3L
)
writeData(
  workbook,
  "02_SelectedClustering",
  "Global validation metrics",
  startRow = 3L,
  startCol = 4L
)
addStyle(
  workbook,
  "02_SelectedClustering",
  section_style,
  rows = 3L,
  cols = 4:7,
  gridExpand = TRUE,
  stack = TRUE
)
writeDataTable(
  workbook,
  "02_SelectedClustering",
  selected_global_metrics_left,
  startRow = 4L,
  startCol = 4L,
  tableName = "SelectedGlobalMetricsLeftTable",
  tableStyle = "TableStyleMedium2",
  withFilter = FALSE
)
writeDataTable(
  workbook,
  "02_SelectedClustering",
  selected_global_metrics_right,
  startRow = 4L,
  startCol = 6L,
  tableName = "SelectedGlobalMetricsRightTable",
  tableStyle = "TableStyleMedium2",
  withFilter = FALSE
)

sample_metrics_start_row <- 20L
writeData(
  workbook,
  "02_SelectedClustering",
  "Validation metrics per sample",
  startRow = sample_metrics_start_row,
  startCol = 1L
)
addStyle(
  workbook,
  "02_SelectedClustering",
  section_style,
  rows = sample_metrics_start_row,
  cols = 1:ncol(selected_metrics_export),
  gridExpand = TRUE,
  stack = TRUE
)
writeDataTable(
  workbook,
  "02_SelectedClustering",
  selected_metrics_export,
  startRow = sample_metrics_start_row + 1L,
  startCol = 1L,
  tableName = "SelectedMetricsBySampleTable",
  tableStyle = "TableStyleMedium2",
  withFilter = TRUE
)

cluster_sizes_start_row <-
  sample_metrics_start_row + nrow(selected_metrics_export) + 5L
writeData(
  workbook,
  "02_SelectedClustering",
  "Cluster sizes per sample",
  startRow = cluster_sizes_start_row,
  startCol = 1L
)
addStyle(
  workbook,
  "02_SelectedClustering",
  section_style,
  rows = cluster_sizes_start_row,
  cols = 1:ncol(selected_cluster_sizes_export),
  gridExpand = TRUE,
  stack = TRUE
)
writeDataTable(
  workbook,
  "02_SelectedClustering",
  selected_cluster_sizes_export,
  startRow = cluster_sizes_start_row + 1L,
  startCol = 1L,
  tableName = "SelectedClusterSizesTable",
  tableStyle = "TableStyleMedium2",
  withFilter = TRUE
)
freezePane(
  workbook,
  "02_SelectedClustering",
  firstActiveRow = sample_metrics_start_row + 2L,
  firstActiveCol = 2L
)
setColWidths(
  workbook,
  "02_SelectedClustering",
  cols = 1:max(
    ncol(selected_metrics_export),
    ncol(selected_cluster_sizes_export)
  ),
  widths = "auto"
)

write_table_sheet(
  workbook,
  "03_MetricsPerSample",
  all_metrics_export,
  "AllMetricsBySampleTable",
  freeze_first_column = TRUE
)
write_table_sheet(
  workbook,
  "04_ClusterSizePerSample",
  all_cluster_sizes_export,
  "AllClusterSizesTable",
  freeze_first_column = TRUE
)
write_table_sheet(
  workbook,
  "05_ARI_Details",
  all_ari_details,
  "ARIDetailsTable",
  freeze_first_column = TRUE
)
addWorksheet(workbook, "06_Notes", gridLines = FALSE)
notes <- data.table(
  item = c(
    "Scope",
    "Selected clustering",
    "Spatial ASW",
    "Transcriptomic ASW SD",
    "CHAOS and PAS SD",
    "ARI adjacent mean",
    "ARI adjacent SD",
    "Cluster sizes",
    "Cluster identifiers"
  ),
  description = c(
    paste0(
      "All integration methods; dims 20 and 30; k 20 and 30; ",
      "Louvain, Louvain refined, SLM and Leiden; resolutions 0.10-1.00."
    ),
    "CCA, dims20, k20, Leiden, resolution 0.40.",
    "Excluded from the workbook as requested.",
    "Copied directly from transcriptomicASW_sd in the validation summary.",
    "Copied directly from SD across samples in the validation summary.",
    paste0(
      "Derived as the mean of ARI versus the previous and next resolution. ",
      "At resolution boundaries it equals the single available adjacent ARI."
    ),
    paste0(
      "Derived as the sample SD of ARI versus the previous and next ",
      "resolution; NA when only one adjacent ARI is available."
    ),
    paste0(
      "Calculated only by counting spot assignments for each ",
      "sample_ID x cluster combination."
    ),
    "clusterID preserves the original zero-based cluster ID from the source assignment table."
  )
)
writeDataTable(
  workbook,
  "06_Notes",
  notes,
  startRow = 1L,
  startCol = 1L,
  tableName = "NotesTable",
  tableStyle = "TableStyleMedium2",
  withFilter = FALSE
)
setColWidths(workbook, "06_Notes", cols = 1L, widths = 34)
setColWidths(workbook, "06_Notes", cols = 2L, widths = 90)
setRowHeights(workbook, "06_Notes", rows = 2L:(nrow(notes) + 1L), heights = 36)
addStyle(
  workbook,
  "06_Notes",
  createStyle(wrapText = TRUE, valign = "top"),
  rows = 2L:(nrow(notes) + 1L),
  cols = 1:2,
  gridExpand = TRUE,
  stack = TRUE
)
freezePane(workbook, "06_Notes", firstRow = TRUE)

apply_numeric_format(
  workbook,
  "01_AllConfigurations",
  all_configurations_export,
  c("resolution"),
  "0.00"
)
apply_numeric_format(
  workbook,
  "01_AllConfigurations",
  all_configurations_export,
  c(
    "transcriptomicASW_mean",
    "transcriptomicASW_sd",
    "transcriptomicASW_fractionNegative",
    "CHAOS_meanAcrossSamples",
    "CHAOS_sdAcrossSamples",
    "PAS_meanAcrossSamples",
    "PAS_sdAcrossSamples",
    "PAS_meanPercentAcrossSamples",
    "ARI_vsPreviousResolution",
    "ARI_vsNextResolution",
    "ARI_adjacentMean",
    "ARI_adjacentSD"
  ),
  "0.#####"
)
apply_numeric_format(
  workbook,
  "03_MetricsPerSample",
  all_metrics_export,
  c("resolution"),
  "0.00"
)
apply_numeric_format(
  workbook,
  "03_MetricsPerSample",
  all_metrics_export,
  c("CHAOS", "PAS_percent"),
  "0.#####"
)
apply_numeric_format(
  workbook,
  "04_ClusterSizePerSample",
  all_cluster_sizes_export,
  c("resolution"),
  "0.00"
)
apply_numeric_format(
  workbook,
  "05_ARI_Details",
  all_ari_details,
  c("resolution1", "resolution2"),
  "0.00"
)
apply_numeric_format(
  workbook,
  "05_ARI_Details",
  all_ari_details,
  c("ARI"),
  "0.#####"
)

selected_global_value_style <- createStyle(numFmt = "0.#####")
addStyle(
  workbook,
  "02_SelectedClustering",
  selected_global_value_style,
  rows = 5L:(4L + nrow(selected_global_metrics_left)),
  cols = 5L,
  gridExpand = TRUE,
  stack = TRUE
)
addStyle(
  workbook,
  "02_SelectedClustering",
  selected_global_value_style,
  rows = 5L:(4L + nrow(selected_global_metrics_right)),
  cols = 7L,
  gridExpand = TRUE,
  stack = TRUE
)
# The helper above assumes tables begin in row 1, so format the selected tables directly.
selected_metrics_numeric_style <- createStyle(numFmt = "0.#####")
selected_metrics_numeric_columns <- intersect(
  c("CHAOS", "PAS_percent"),
  names(selected_metrics_export)
)
if (length(selected_metrics_numeric_columns) > 0L) {
  addStyle(
    workbook,
    "02_SelectedClustering",
    selected_metrics_numeric_style,
    rows = (sample_metrics_start_row + 2L):(
      sample_metrics_start_row + 1L + nrow(selected_metrics_export)
    ),
    cols = match(
      selected_metrics_numeric_columns,
      names(selected_metrics_export)
    ),
    gridExpand = TRUE,
    stack = TRUE
  )
}
selected_resolution_style <- createStyle(numFmt = "0.00")
addStyle(
  workbook,
  "02_SelectedClustering",
  selected_resolution_style,
  rows = (sample_metrics_start_row + 2L):(
    sample_metrics_start_row + 1L + nrow(selected_metrics_export)
  ),
  cols = match("resolution", names(selected_metrics_export)),
  gridExpand = TRUE,
  stack = TRUE
)
addStyle(
  workbook,
  "02_SelectedClustering",
  selected_resolution_style,
  rows = (cluster_sizes_start_row + 2L):(
    cluster_sizes_start_row + 1L + nrow(selected_cluster_sizes_export)
  ),
  cols = match("resolution", names(selected_cluster_sizes_export)),
  gridExpand = TRUE,
  stack = TRUE
)

saveWorkbook(workbook, output_xlsx, overwrite = TRUE)

message_line("Completed")
message("Workbook saved to:")
message(output_xlsx)
message("AllConfigurations rows: ", nrow(all_configurations))
message("AllMetricsBySample rows: ", nrow(all_metrics_by_sample))
message("AllClusterSizes rows: ", nrow(all_cluster_sizes))
message("ARI detail rows: ", nrow(all_ari_details))
message("Selected per-sample rows: ", nrow(selected_metrics))
message("Selected cluster-size rows: ", nrow(selected_cluster_sizes))
message("Missing or ambiguous input files: ", sum(input_status$status != "OK"))
