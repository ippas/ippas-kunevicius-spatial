#!/usr/bin/env Rscript

message("VERSION: 02 RCTD deconvolution | all n16 samples test | full mode | Allen WMB n14 reference | per-spot normalized proportions")

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(spacexr)
})

# ==============================================================================
# 1. PARAMETERS
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"
clustering_name <- "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

sample_column <- "sample_ID"
cluster_column <- "leiden_res0.40"
assay_name <- "RNA"
samples_to_run <- c(
  "1_1F", "1_1Fd", "1_1M", "2_1M", "2_1Md",
  "3_1F", "3_1M", "5_1M", "5_3F", "12_1M",
  "13_1F", "13_1M", "15_1F", "18_1F", "18_1M", "23_1F"
)

reference_file <- file.path(
  project_root,
  "data", "reference", "allenWholeMouseBrain_Yao2023", "processed",
  "broadCellTypeReference_n14CellTypes",
  "allenWMB_broadCellTypeReference_n14CellTypes_max10000PerType_seed7.rds"
)

input_rdata_file <- file.path(
  project_root,
  "results", dataset_name, "seurat_clustering_analysis",
  "logNormalize_vst", "hvg2000", "02_cca", "dims20", "k20", "prune0067", "RData",
  paste0(
    "02_", dataset_name, "_logNormalizeVst_hvg2000_",
    "ccaDims20K20PruneSNN0067_",
    "res010to100by010_multiClusteringAndUmap.RData"
  )
)

output_root <- file.path(
  project_root,
  "results", dataset_name, "cellTypeDeconvolution_RCTD",
  clustering_name, "01_test_allSamples_n16"
)

output_dirs <- list(
  perSpot = file.path(output_root, "perSpot"),
  perCluster = file.path(output_root, "perCluster"),
  RDS = file.path(output_root, "RDS"),
  plots = file.path(output_root, "plots"),
  logs = file.path(output_root, "logs")
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

physical_cores <- parallel::detectCores(logical = FALSE)
if (is.na(physical_cores) || physical_cores < 1L) physical_cores <- 1L
max_cores <- max(1L, floor(physical_cores / 2L))

message("Reference: ", reference_file)
message("Seurat RData: ", input_rdata_file)
message("Samples: ", paste(samples_to_run, collapse = ", "))
message("Physical cores detected: ", physical_cores)
message("RCTD max_cores: ", max_cores)
message("Output: ", output_root)

# ==============================================================================
# 2. HELPERS
# ==============================================================================

as_dgCMatrix <- function(x) {
  if (inherits(x, "dgCMatrix")) return(x)
  methods::as(x, "dgCMatrix")
}

get_raw_counts_for_cells <- function(object, assay, cells) {
  assay_object <- object[[assay]]

  # Seurat v5: collect all count layers containing the requested spots.
  if ("Layers" %in% getNamespaceExports("SeuratObject") &&
      "LayerData" %in% getNamespaceExports("SeuratObject")) {

    count_layers <- tryCatch(
      SeuratObject::Layers(assay_object, search = "^counts"),
      error = function(e) character()
    )

    if (length(count_layers) > 0L) {
      matrices <- lapply(count_layers, function(layer_name) {
        x <- SeuratObject::LayerData(assay_object, layer = layer_name)
        selected <- intersect(cells, colnames(x))
        if (length(selected) == 0L) return(NULL)
        x[, selected, drop = FALSE]
      })
      matrices <- Filter(Negate(is.null), matrices)

      if (length(matrices) > 0L) {
        counts <- if (length(matrices) == 1L) matrices[[1]] else do.call(cbind, matrices)
        if (anyDuplicated(colnames(counts))) {
          stop("Duplicated spot barcodes were found across RNA count layers.")
        }
        missing_cells <- setdiff(cells, colnames(counts))
        if (length(missing_cells) > 0L) {
          stop(
            "Raw counts were not found for ", length(missing_cells),
            " requested spots. Example: ", paste(head(missing_cells, 5L), collapse = ", ")
          )
        }
        return(as_dgCMatrix(counts[, cells, drop = FALSE]))
      }
    }
  }

  # Seurat v4 / single-layer fallback.
  counts <- tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = "counts"),
    error = function(e1) Seurat::GetAssayData(object, assay = assay, slot = "counts")
  )

  missing_cells <- setdiff(cells, colnames(counts))
  if (length(missing_cells) > 0L) {
    stop(
      "Raw counts were not found for ", length(missing_cells),
      " requested spots. Example: ", paste(head(missing_cells, 5L), collapse = ", ")
    )
  }

  as_dgCMatrix(counts[, cells, drop = FALSE])
}

find_image_for_cells <- function(object, sample_cells) {
  image_names <- Seurat::Images(object)
  if (length(image_names) == 0L) stop("No spatial images were found in the Seurat object.")

  overlap_n <- vapply(image_names, function(image_name) {
    image_cells <- tryCatch(
      SeuratObject::Cells(object[[image_name]]),
      error = function(e) character()
    )
    length(intersect(sample_cells, image_cells))
  }, integer(1))

  best_index <- which.max(overlap_n)
  if (length(best_index) == 0L || overlap_n[best_index] == 0L) {
    stop("Could not match sample spots to any Seurat spatial image.")
  }

  best_name <- image_names[best_index]
  message(
    "Matched spatial image: ", best_name,
    " (", overlap_n[best_index], "/", length(sample_cells), " sample spots)"
  )
  best_name
}

get_visium_coordinates <- function(object, image_name, sample_cells) {
  image_object <- object[[image_name]]

  coords <- tryCatch(
    SeuratObject::GetTissueCoordinates(image_object, scale = NULL),
    error = function(e) Seurat::GetTissueCoordinates(image_object, scale = NULL)
  )
  coords <- as.data.frame(coords)

  if (length(intersect(rownames(coords), sample_cells)) == 0L) {
    id_candidates <- c("cell", "barcode", "barcodes", "ID")
    id_candidates <- id_candidates[id_candidates %in% colnames(coords)]
    if (length(id_candidates) > 0L) {
      rownames(coords) <- as.character(coords[[id_candidates[1]]])
    }
  }

  missing_cells <- setdiff(sample_cells, rownames(coords))
  if (length(missing_cells) > 0L) {
    stop(
      "Coordinates are missing for ", length(missing_cells),
      " sample spots. Example: ", paste(head(missing_cells, 5L), collapse = ", ")
    )
  }

  coords <- coords[sample_cells, , drop = FALSE]

  if (all(c("x", "y") %in% colnames(coords))) {
    xy <- coords[, c("x", "y"), drop = FALSE]
  } else if (all(c("imagecol", "imagerow") %in% colnames(coords))) {
    xy <- data.frame(
      x = coords[["imagecol"]], y = coords[["imagerow"]], row.names = rownames(coords)
    )
  } else if (all(c("pxl_col_in_fullres", "pxl_row_in_fullres") %in% colnames(coords))) {
    xy <- data.frame(
      x = coords[["pxl_col_in_fullres"]],
      y = coords[["pxl_row_in_fullres"]],
      row.names = rownames(coords)
    )
  } else {
    numeric_columns <- names(coords)[vapply(coords, is.numeric, logical(1))]
    if (length(numeric_columns) < 2L) {
      stop("Could not identify x/y coordinate columns.")
    }
    message("Using coordinate columns: ", paste(numeric_columns[1:2], collapse = ", "))
    xy <- coords[, numeric_columns[1:2], drop = FALSE]
    colnames(xy) <- c("x", "y")
  }

  xy[["x"]] <- as.numeric(xy[["x"]])
  xy[["y"]] <- as.numeric(xy[["y"]])
  if (any(!is.finite(xy[["x"]])) || any(!is.finite(xy[["y"]]))) {
    stop("Non-finite spatial coordinates were detected.")
  }
  xy
}

load_seurat_object <- function(rdata_file) {
  if (!file.exists(rdata_file)) stop("Seurat RData file does not exist: ", rdata_file)

  e <- new.env(parent = emptyenv())
  loaded_names <- load(rdata_file, envir = e)
  candidates <- loaded_names[vapply(loaded_names, function(x) inherits(e[[x]], "Seurat"), logical(1))]

  if (length(candidates) == 0L) {
    stop("No Seurat object was found in RData. Loaded objects: ", paste(loaded_names, collapse = ", "))
  }

  if (length(candidates) > 1L) {
    n_cells <- vapply(candidates, function(x) ncol(e[[x]]), integer(1))
    selected <- candidates[which.max(n_cells)]
    message("Multiple Seurat objects found. Using the largest: ", selected)
  } else {
    selected <- candidates
    message("Seurat object: ", selected)
  }

  e[[selected]]
}

prepare_reference_vectors <- function(reference_data) {
  required_names <- c("counts", "cell_types", "nUMI")
  missing_names <- setdiff(required_names, names(reference_data))
  if (length(missing_names) > 0L) {
    stop("Reference RDS is missing: ", paste(missing_names, collapse = ", "))
  }

  counts <- as_dgCMatrix(reference_data[["counts"]])
  cell_types <- reference_data[["cell_types"]]
  nUMI <- reference_data[["nUMI"]]

  if (is.null(colnames(counts)) || is.null(rownames(counts))) {
    stop("Reference counts must have gene and cell names.")
  }

  if (is.null(names(cell_types))) names(cell_types) <- colnames(counts)
  cell_types <- cell_types[colnames(counts)]
  if (anyNA(cell_types)) stop("Reference cell_types could not be aligned to count-matrix columns.")
  cell_types <- factor(as.character(cell_types))
  names(cell_types) <- colnames(counts)

  if (is.null(names(nUMI))) names(nUMI) <- colnames(counts)
  nUMI <- as.numeric(nUMI[colnames(counts)])
  names(nUMI) <- colnames(counts)
  if (anyNA(nUMI)) stop("Reference nUMI could not be aligned to count-matrix columns.")

  list(counts = counts, cell_types = cell_types, nUMI = nUMI)
}

plot_spatial_weights <- function(per_spot_long, sample_id, output_file) {
  p <- ggplot(per_spot_long, aes(x = x, y = y, color = proportion)) +
    geom_point(size = 0.65, alpha = 0.95) +
    scale_color_viridis_c(limits = c(0, 1), oob = scales::squish) +
    scale_y_reverse() +
    coord_fixed() +
    facet_wrap(~ broad_cell_type, ncol = 4) +
    labs(
      title = paste0(sample_id, " | RCTD full-mode proportions"),
      x = NULL, y = NULL, color = "Proportion"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(size = 8)
    )

  ggsave(output_file, p, width = 12, height = 10, units = "in")
}

# ==============================================================================
# 3. CHECK RCTD API
# ==============================================================================

if (!file.exists(reference_file)) stop("Reference RDS does not exist: ", reference_file)

spacexr_namespace <- asNamespace("spacexr")
has_bioconductor_api <-
  exists("createRctd", envir = spacexr_namespace, mode = "function", inherits = FALSE) &&
  exists("runRctd", envir = spacexr_namespace, mode = "function", inherits = FALSE)
has_legacy_api <-
  exists("create.RCTD", envir = spacexr_namespace, mode = "function", inherits = FALSE) &&
  exists("run.RCTD", envir = spacexr_namespace, mode = "function", inherits = FALSE)

if (!has_bioconductor_api && !has_legacy_api) {
  stop("Installed spacexr exposes neither createRctd/runRctd nor create.RCTD/run.RCTD.")
}

rctd_api <- if (has_bioconductor_api) "Bioconductor" else "legacy"
message("spacexr version: ", as.character(utils::packageVersion("spacexr")))
message("RCTD API selected: ", rctd_api)

# ==============================================================================
# 4. LOAD REFERENCE
# ==============================================================================

message("\n============================================================")
message("Loading Allen WMB reference")
message("============================================================")

reference_data <- readRDS(reference_file)
reference_vectors <- prepare_reference_vectors(reference_data)
reference_counts <- reference_vectors[["counts"]]
reference_cell_types <- reference_vectors[["cell_types"]]
reference_nUMI <- reference_vectors[["nUMI"]]

message(
  "Reference dimensions: ", format(nrow(reference_counts), big.mark = ","),
  " genes x ", format(ncol(reference_counts), big.mark = ","), " cells"
)
print(table(reference_cell_types))

if (rctd_api == "Bioconductor") {
  required_bioc <- c("SummarizedExperiment", "SpatialExperiment", "S4Vectors")
  missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc) > 0L) {
    stop("Missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  }

  reference_coldata <- S4Vectors::DataFrame(
    cell_type = reference_cell_types,
    nUMI = reference_nUMI,
    row.names = colnames(reference_counts)
  )

  reference_object <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = reference_counts),
    colData = reference_coldata
  )
} else {
  reference_object <- spacexr::Reference(
    counts = reference_counts,
    cell_types = reference_cell_types,
    nUMI = reference_nUMI
  )
}

rm(reference_data, reference_vectors, reference_counts, reference_cell_types, reference_nUMI)
gc()

# ==============================================================================
# 5. LOAD SEURAT OBJECT
# ==============================================================================

message("\n============================================================")
message("Loading spatial Seurat object")
message("============================================================")

seurat_object <- load_seurat_object(input_rdata_file)
metadata <- seurat_object[[]]

if (!sample_column %in% colnames(metadata)) stop("Missing metadata column: ", sample_column)
if (!cluster_column %in% colnames(metadata)) stop("Missing metadata column: ", cluster_column)

missing_samples <- setdiff(samples_to_run, unique(as.character(metadata[[sample_column]])))
if (length(missing_samples) > 0L) {
  stop("Requested samples not found: ", paste(missing_samples, collapse = ", "))
}

message("Total spots in Seurat object: ", format(ncol(seurat_object), big.mark = ","))

# ==============================================================================
# 6. RUN RCTD FOR EACH SAMPLE
# ==============================================================================

combined_per_spot <- list()
combined_per_cluster <- list()
combined_sample_summary <- list()

for (sample_i in seq_along(samples_to_run)) {

  sample_id <- samples_to_run[sample_i]

  message("\n============================================================")
  message("RCTD SAMPLE [", sample_i, "/", length(samples_to_run), "]: ", sample_id)
  message("============================================================")

  sample_cells <- rownames(metadata)[as.character(metadata[[sample_column]]) == sample_id]
  if (length(sample_cells) == 0L) stop("No spots found for sample: ", sample_id)
  message("Sample spots before RCTD filtering: ", format(length(sample_cells), big.mark = ","))

  message("Extracting raw RNA counts...")
  spatial_counts <- get_raw_counts_for_cells(seurat_object, assay_name, sample_cells)
  if (any(spatial_counts@x < 0)) stop("Negative values detected in spatial raw counts.")
  message(
    "Spatial count matrix: ", format(nrow(spatial_counts), big.mark = ","),
    " genes x ", format(ncol(spatial_counts), big.mark = ","), " spots"
  )

  message("Matching Visium image and coordinates...")
  image_name <- find_image_for_cells(seurat_object, sample_cells)
  coords <- get_visium_coordinates(seurat_object, image_name, sample_cells)

  spatial_counts <- spatial_counts[, rownames(coords), drop = FALSE]
  sample_cells <- colnames(spatial_counts)
  nUMI_spatial <- Matrix::colSums(spatial_counts)
  names(nUMI_spatial) <- sample_cells

  message(
    "Spatial UMI range: ", format(min(nUMI_spatial), big.mark = ","),
    " - ", format(max(nUMI_spatial), big.mark = ",")
  )

  if (rctd_api == "Bioconductor") {

    spatial_coldata <- S4Vectors::DataFrame(
      sample_ID = rep(sample_id, length(sample_cells)),
      nUMI = as.numeric(nUMI_spatial),
      row.names = sample_cells
    )

    spatial_object <- SpatialExperiment::SpatialExperiment(
      assays = list(counts = spatial_counts),
      colData = spatial_coldata,
      spatialCoords = as.matrix(coords[, c("x", "y"), drop = FALSE])
    )

    message("Running createRctd()...")
    rctd_data <- spacexr::createRctd(
      spatial_experiment = spatial_object,
      reference_experiment = reference_object,
      cell_type_col = "cell_type",
      ref_n_cells_min = 25,
      ref_n_cells_max = 10000
    )

    message("Running runRctd(..., rctd_mode = 'full')...")
    rctd_result <- spacexr::runRctd(
      rctd_data = rctd_data,
      rctd_mode = "full",
      max_cores = max_cores
    )

    weights <- as.matrix(SummarizedExperiment::assay(rctd_result, "weights"))
    proportions <- t(weights)

  } else {

    spatial_object <- spacexr::SpatialRNA(
      coords = coords,
      counts = spatial_counts,
      nUMI = nUMI_spatial
    )

    message("Running create.RCTD(test_mode = FALSE)...")
    rctd_data <- spacexr::create.RCTD(
      spatialRNA = spatial_object,
      reference = reference_object,
      max_cores = max_cores,
      test_mode = FALSE,
      CELL_MIN_INSTANCE = 25
    )

    message("Running run.RCTD(..., doublet_mode = 'full')...")
    rctd_result <- spacexr::run.RCTD(rctd_data, doublet_mode = "full")
    raw_weights <- as.matrix(rctd_result@results$weights)
    proportions <- as.matrix(spacexr::normalize_weights(raw_weights))
  }

  if (nrow(proportions) == 0L || ncol(proportions) == 0L) {
    stop("RCTD returned an empty proportion matrix for sample ", sample_id)
  }

  # Ensure rows = spots and columns = cell populations.
  if (length(intersect(rownames(proportions), sample_cells)) == 0L &&
      length(intersect(colnames(proportions), sample_cells)) > 0L) {
    proportions <- t(proportions)
  }

  retained_spots <- intersect(sample_cells, rownames(proportions))
  if (length(retained_spots) == 0L) {
    stop("RCTD results could not be matched back to sample spot barcodes.")
  }
  proportions <- proportions[retained_spots, , drop = FALSE]

  # Convert RCTD weights to relative proportions that sum to 1 within each spot.
  raw_row_sums <- rowSums(proportions)
  if (any(!is.finite(raw_row_sums)) || any(raw_row_sums <= 0)) {
    stop("Non-finite or non-positive RCTD weight sum detected for sample ", sample_id)
  }
  proportions <- proportions / raw_row_sums
  row_sums <- rowSums(proportions)

  message(
    "RCTD retained spots: ", format(nrow(proportions), big.mark = ","),
    " / ", format(length(sample_cells), big.mark = ",")
  )
  message(
    "Raw RCTD weight row-sum range: ", signif(min(raw_row_sums), 6),
    " - ", signif(max(raw_row_sums), 6)
  )
  message(
    "Normalized proportion row-sum range: ", signif(min(row_sums), 6),
    " - ", signif(max(row_sums), 6)
  )

  # Save full RCTD result object.
  rctd_rds_file <- file.path(output_dirs$RDS, paste0(sample_id, "_RCTD_full.rds"))
  saveRDS(rctd_result, rctd_rds_file, compress = FALSE)
  message("Saved RCTD object: ", rctd_rds_file)

  # Per-spot wide table.
  per_spot <- data.table(
    spot = rownames(proportions),
    sample_ID = sample_id,
    image = image_name,
    x = coords[rownames(proportions), "x"],
    y = coords[rownames(proportions), "y"],
    cluster = as.character(metadata[rownames(proportions), cluster_column, drop = TRUE])
  )
  proportions_dt <- as.data.table(proportions)
  setnames(proportions_dt, colnames(proportions))
  per_spot <- cbind(per_spot, proportions_dt)

  fwrite(
    per_spot,
    file.path(output_dirs$perSpot, paste0(sample_id, "_RCTD_full_cellTypeProportions_perSpot.tsv")),
    sep = "\t"
  )

  # Long table.
  cell_type_columns <- colnames(proportions)
  per_spot_long <- melt(
    per_spot,
    id.vars = c("spot", "sample_ID", "image", "x", "y", "cluster"),
    measure.vars = cell_type_columns,
    variable.name = "broad_cell_type",
    value.name = "proportion"
  )

  fwrite(
    per_spot_long,
    file.path(output_dirs$perSpot, paste0(sample_id, "_RCTD_full_cellTypeProportions_perSpot_long.tsv")),
    sep = "\t"
  )

  # Descriptive summary per Seurat cluster.
  per_cluster <- per_spot_long[
    ,
    .(
      n_spots = uniqueN(spot),
      mean_proportion = mean(proportion, na.rm = TRUE),
      sd_proportion = sd(proportion, na.rm = TRUE),
      median_proportion = median(proportion, na.rm = TRUE)
    ),
    by = .(sample_ID, cluster, broad_cell_type)
  ]

  fwrite(
    per_cluster,
    file.path(output_dirs$perCluster, paste0(sample_id, "_RCTD_full_cellTypeProportions_perCluster.tsv")),
    sep = "\t"
  )

  # Whole-section descriptive summary.
  sample_summary <- per_spot_long[
    ,
    .(
      n_spots = uniqueN(spot),
      mean_proportion = mean(proportion, na.rm = TRUE),
      sd_proportion = sd(proportion, na.rm = TRUE),
      median_proportion = median(proportion, na.rm = TRUE)
    ),
    by = .(sample_ID, broad_cell_type)
  ]

  fwrite(
    sample_summary,
    file.path(output_dirs$perSpot, paste0(sample_id, "_RCTD_full_cellTypeProportions_sectionSummary.tsv")),
    sep = "\t"
  )

  # Spatial QC plot: all 14 populations.
  plot_file <- file.path(
    output_dirs$plots,
    paste0(sample_id, "_RCTD_full_cellTypeProportions_spatial.pdf")
  )
  plot_spatial_weights(per_spot_long, sample_id, plot_file)
  message("Saved spatial QC plot: ", plot_file)

  combined_per_spot[[sample_id]] <- per_spot_long
  combined_per_cluster[[sample_id]] <- per_cluster
  combined_sample_summary[[sample_id]] <- sample_summary

  rm(
    spatial_counts, coords, nUMI_spatial, spatial_object, rctd_data,
    rctd_result, proportions, per_spot, proportions_dt, per_spot_long,
    per_cluster, sample_summary
  )
  gc()
  message("Finished sample: ", sample_id)
}

# ==============================================================================
# 7. COMBINE ALL-SAMPLE OUTPUTS
# ==============================================================================

message("\n============================================================")
message("Combining all-sample test outputs")
message("============================================================")

combined_per_spot_dt <- rbindlist(combined_per_spot, use.names = TRUE, fill = TRUE)
combined_per_cluster_dt <- rbindlist(combined_per_cluster, use.names = TRUE, fill = TRUE)
combined_sample_summary_dt <- rbindlist(combined_sample_summary, use.names = TRUE, fill = TRUE)

fwrite(
  combined_per_spot_dt,
  file.path(output_dirs$perSpot, "allSamples_n16_RCTD_full_cellTypeProportions_perSpot_long.tsv"),
  sep = "\t"
)
fwrite(
  combined_per_cluster_dt,
  file.path(output_dirs$perCluster, "allSamples_n16_RCTD_full_cellTypeProportions_perCluster.tsv"),
  sep = "\t"
)
fwrite(
  combined_sample_summary_dt,
  file.path(output_dirs$perSpot, "allSamples_n16_RCTD_full_cellTypeProportions_sectionSummary.tsv"),
  sep = "\t"
)

run_parameters <- data.table(
  parameter = c(
    "dataset_name", "clustering_name", "sample_column", "cluster_column",
    "assay_name", "n_samples", "samples", "reference_file", "rctd_mode", "normalization",
    "spacexr_version", "rctd_api", "physical_cores", "max_cores"
  ),
  value = c(
    dataset_name, clustering_name, sample_column, cluster_column,
    assay_name, length(samples_to_run), paste(samples_to_run, collapse = ","), reference_file, "full", "per_spot_sum1",
    as.character(utils::packageVersion("spacexr")), rctd_api,
    physical_cores, max_cores
  )
)

fwrite(run_parameters, file.path(output_dirs$logs, "runParameters.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(output_dirs$logs, "sessionInfo.txt"))

message("\n============================================================")
message("RCTD ALL-SAMPLE TEST FINISHED")
message("============================================================")
message("Samples: ", paste(samples_to_run, collapse = ", "))
message("Mode: full")
message("Reference: 14 prepared Allen WMB cell populations")
message("Output: ", output_root)
message("============================================================")
