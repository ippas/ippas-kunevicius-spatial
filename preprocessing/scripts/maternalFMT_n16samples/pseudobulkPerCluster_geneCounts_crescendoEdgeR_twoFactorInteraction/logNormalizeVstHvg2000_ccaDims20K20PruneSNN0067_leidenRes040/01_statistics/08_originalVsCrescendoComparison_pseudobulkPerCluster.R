#!/usr/bin/env Rscript
# ============================================================================== 
# 03_originalVsCrescendoComparison_pseudobulkPerCluster.R
#
# Prepare one clean ORIGINAL vs CRESCENDO comparison dataset for later heatmaps.
#
# COMMON GENE UNIVERSE
#   Genes retained by Crescendo AND protein-coding on chromosomes 1-19 or X.
#   The same gene universe and same 16 samples are used for both variants.
#
# OUTPUT PER VARIANT
#   01_pseudobulkCounts/   per-cluster sample pseudobulk counts
#   02_log2plus1/          log2(count + 1)
#   03_log2plus1_QN/       quantile normalization across 16 samples per cluster
#   04_detectionPercent/   % spots count > 0 per sample + mean + SD
#   05_edgeR/              one XLSX, 9 tests, all tested genes, all clusters
#
# edgeR model:
#   expression ~ fmt_donor_group * sex
# ============================================================================== 

options(stringsAsFactors = FALSE)

SCRIPT_BUILD <- "20260903_V2_SYMBOL_MAPPING"

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "tibble",
  "edgeR", "limma", "preprocessCore", "openxlsx"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(edgeR)
  library(limma)
  library(preprocessCore)
  library(openxlsx)
})

# ============================================================================== 
# 1. Configuration
# ============================================================================== 

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"
clustering_name <- "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

sample_column <- "sample_ID"
assay_name <- "RNA"

expected_number_of_samples <- 16L
expected_number_of_clusters <- 16L
analysis_chromosomes <- c(as.character(1:19), "X")

filter_by_expr_parameters <- list(
  min.count = 10,
  min.total.count = 15,
  large.n = 10,
  min.prop = 0.7
)

# ============================================================================== 
# 2. Input paths
# ============================================================================== 

input_rdata_file <- file.path(
  project_root,
  "results",
  dataset_name,
  "seurat_clustering_analysis",
  "logNormalize_vst",
  "hvg2000",
  "02_cca",
  "dims20",
  "k20",
  "prune0067",
  "RData",
  paste0(
    "02_", dataset_name,
    "_logNormalizeVst_hvg2000_ccaDims20K20PruneSNN0067_",
    "res010to100by010_multiClusteringAndUmap.RData"
  )
)

annotation_rds_file <- file.path(
  project_root,
  "data",
  "ensembl115",
  "ensembl115_mouse_proteinCodingGenes.rds"
)

crescendo_statistics_root <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)

crescendo_root <- file.path(
  crescendo_statistics_root,
  "01_crescendoBatchCorrection"
)

crescendo_input_rds <- file.path(
  crescendo_root,
  "07_objects",
  paste0(
    dataset_name,
    "_",
    clustering_name,
    "_crescendoBatchCorrection_pseudobulkCounts_beforeEdgeR.rds"
  )
)

corrected_spot_chunks_dir <- file.path(
  crescendo_root,
  "03_correctedSpotCounts_chunks_v13"
)

# ============================================================================== 
# 3. Output structure
# ============================================================================== 

comparison_root <- file.path(
  crescendo_statistics_root,
  "03_originalVsCrescendoComparison_pseudobulkPerCluster"
)

original_root <- file.path(comparison_root, "01_original")
crescendo_comparison_root <- file.path(comparison_root, "02_crescendo")

make_variant_dirs <- function(root_dir) {
  dirs <- list(
    pseudobulk = file.path(root_dir, "01_pseudobulkCounts"),
    log2plus1 = file.path(root_dir, "02_log2plus1"),
    log2plus1_qn = file.path(root_dir, "03_log2plus1_QN"),
    detection = file.path(root_dir, "04_detectionPercent"),
    edger = file.path(root_dir, "05_edgeR")
  )

  invisible(lapply(c(root_dir, unname(dirs)), dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

original_dirs <- make_variant_dirs(original_root)
crescendo_dirs <- make_variant_dirs(crescendo_comparison_root)

# ============================================================================== 
# 4. Helpers
# ============================================================================== 

stage_message <- function(stage_number, text) {
  message("\n============================================================")
  message("STAGE ", stage_number, ": ", text)
  message("============================================================")
}

write_tsv_gz <- function(x, output_file) {
  con <- gzfile(output_file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  utils::write.table(
    x,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  invisible(output_file)
}

natural_cluster_order <- function(x) {
  x <- unique(as.character(x))
  x_num <- suppressWarnings(as.integer(x))
  if (all(!is.na(x_num))) return(x[order(x_num)])
  sort(x)
}

safe_row_sd <- function(x) {
  if (ncol(x) <= 1L) return(rep(NA_real_, nrow(x)))
  apply(x, 1L, stats::sd)
}

load_seurat_object <- function(input_file) {
  env <- new.env(parent = emptyenv())
  loaded <- load(input_file, envir = env)
  seurat_names <- loaded[
    vapply(loaded, function(x) inherits(get(x, envir = env), "Seurat"), logical(1))
  ]
  if (length(seurat_names) != 1L) {
    stop("Expected exactly one Seurat object. Found: ", paste(seurat_names, collapse = ", "), call. = FALSE)
  }
  get(seurat_names[[1L]], envir = env)
}

get_raw_counts <- function(seurat_object, assay_name) {
  assay_object <- seurat_object[[assay_name]]

  if (inherits(assay_object, "Assay5")) {
    count_layers <- SeuratObject::Layers(assay_object)
    count_layers <- count_layers[count_layers == "counts" | startsWith(count_layers, "counts.")]

    if (length(count_layers) == 0L) {
      stop("No raw count layer found in assay ", assay_name, ".", call. = FALSE)
    }

    if (length(count_layers) > 1L) {
      assay_object <- SeuratObject::JoinLayers(assay_object)
    }

    counts <- SeuratObject::LayerData(assay_object, layer = "counts")
  } else {
    counts <- Seurat::GetAssayData(seurat_object, assay = assay_name, slot = "counts")
  }

  methods::as(counts, "dgCMatrix")
}

matrix_to_output_table <- function(matrix_object, gene_annotation) {
  ids <- rownames(matrix_object)
  annotation_current <- gene_annotation[
    match(ids, gene_annotation$ensembl_gene_id),
    ,
    drop = FALSE
  ]

  if (anyNA(annotation_current$ensembl_gene_id)) {
    stop("Failed to annotate one or more matrix genes.", call. = FALSE)
  }

  dplyr::bind_cols(
    annotation_current |>
      dplyr::select(.data$ensembl_gene_id, .data$gene, .data$chromosome),
    tibble::as_tibble(
      as.data.frame(matrix_object, check.names = FALSE),
      .name_repair = "minimal"
    )
  )
}

create_log2plus1 <- function(count_matrix) {
  x <- as.matrix(count_matrix)
  storage.mode(x) <- "double"
  log2(x + 1)
}

create_quantile_normalized <- function(log2_matrix) {
  x <- as.matrix(log2_matrix)
  storage.mode(x) <- "double"
  qn <- preprocessCore::normalize.quantiles(x, copy = TRUE)
  rownames(qn) <- rownames(x)
  colnames(qn) <- colnames(x)
  qn
}

create_design_and_contrasts <- function(sample_metadata) {
  model_metadata <- sample_metadata |>
    dplyr::transmute(
      sample_ID = as.character(.data$sample_ID),
      fmt_donor_group = factor(as.character(.data$fmt_donor_group), levels = c("Neurotypical", "ASD")),
      sex = factor(as.character(.data$sex), levels = c("Male", "Female"))
    )

  design <- stats::model.matrix(~ fmt_donor_group * sex, data = model_metadata)
  rownames(design) <- model_metadata$sample_ID

  expected_columns <- c(
    "(Intercept)",
    "fmt_donor_groupASD",
    "sexFemale",
    "fmt_donor_groupASD:sexFemale"
  )

  if (!identical(colnames(design), expected_columns)) {
    stop("Unexpected design columns: ", paste(colnames(design), collapse = ", "), call. = FALSE)
  }

  colnames(design) <- c("Intercept", "groupASD", "sexFemale", "groupASD_sexFemale")

  contrasts <- limma::makeContrasts(
    Interaction = groupASD_sexFemale,
    Overall_Sex_Female_vs_Male = sexFemale + 0.5 * groupASD_sexFemale,
    Overall_Group_ASD_vs_Neurotypical = groupASD + 0.5 * groupASD_sexFemale,
    ASD_Male_vs_Neurotypical_Male = groupASD,
    ASD_Female_vs_Neurotypical_Female = groupASD + groupASD_sexFemale,
    Neurotypical_Female_vs_Neurotypical_Male = sexFemale,
    ASD_Female_vs_ASD_Male = sexFemale + groupASD_sexFemale,
    ASD_Female_vs_Neurotypical_Male = groupASD + sexFemale + groupASD_sexFemale,
    ASD_Male_vs_Neurotypical_Female = groupASD - sexFemale,
    levels = design
  )

  test_definitions <- tibble::tribble(
    ~test_id, ~sheet_name, ~comparison,
    "Interaction", "Interaction", "(ASD - Neurotypical)_Female - (ASD - Neurotypical)_Male",
    "Overall_Sex_Female_vs_Male", "Overall_Sex", "Female vs Male averaged equally across donor groups",
    "Overall_Group_ASD_vs_Neurotypical", "Overall_Group", "ASD vs Neurotypical averaged equally across sexes",
    "ASD_Male_vs_Neurotypical_Male", "Group_in_Male", "ASD Male vs Neurotypical Male",
    "ASD_Female_vs_Neurotypical_Female", "Group_in_Female", "ASD Female vs Neurotypical Female",
    "Neurotypical_Female_vs_Neurotypical_Male", "Sex_in_Neurotypical", "Neurotypical Female vs Neurotypical Male",
    "ASD_Female_vs_ASD_Male", "Sex_in_ASD", "ASD Female vs ASD Male",
    "ASD_Female_vs_Neurotypical_Male", "ASD_F_vs_NT_M", "ASD Female vs Neurotypical Male",
    "ASD_Male_vs_Neurotypical_Female", "ASD_M_vs_NT_F", "ASD Male vs Neurotypical Female"
  )

  list(
    metadata = model_metadata,
    design = design,
    contrasts = contrasts,
    test_definitions = test_definitions
  )
}

write_edger_workbook <- function(result_list, test_definitions, output_file) {
  wb <- openxlsx::createWorkbook()
  header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9EAF7", border = "Bottom")
  decimal_style <- openxlsx::createStyle(numFmt = "0.00000")
  scientific_style <- openxlsx::createStyle(numFmt = "0.00000E+00")

  for (i in seq_len(nrow(test_definitions))) {
    test_id <- test_definitions$test_id[[i]]
    sheet_name <- test_definitions$sheet_name[[i]]
    x <- result_list[[test_id]]

    openxlsx::addWorksheet(wb, sheetName = sheet_name)
    openxlsx::writeData(wb, sheet = sheet_name, x = x, headerStyle = header_style, withFilter = TRUE)
    openxlsx::freezePane(wb, sheet = sheet_name, firstActiveRow = 2, firstActiveCol = 5)

    scientific_columns <- which(colnames(x) %in% c("PValue", "FDR"))
    numeric_columns <- which(vapply(x, is.numeric, logical(1)))
    decimal_columns <- setdiff(numeric_columns, scientific_columns)

    if (nrow(x) > 0L && length(decimal_columns) > 0L) {
      openxlsx::addStyle(
        wb, sheet = sheet_name, style = decimal_style,
        rows = 2:(nrow(x) + 1L), cols = decimal_columns,
        gridExpand = TRUE, stack = TRUE
      )
    }

    if (nrow(x) > 0L && length(scientific_columns) > 0L) {
      openxlsx::addStyle(
        wb, sheet = sheet_name, style = scientific_style,
        rows = 2:(nrow(x) + 1L), cols = scientific_columns,
        gridExpand = TRUE, stack = TRUE
      )
    }

    openxlsx::setColWidths(
      wb,
      sheet = sheet_name,
      cols = seq_len(ncol(x)),
      widths = c(10, 20, 18, 12, 12, 12, 12, 14, 14)
    )
  }

  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  invisible(output_file)
}

run_edger_all_clusters <- function(
    counts_by_cluster,
    sample_metadata,
    gene_annotation,
    output_file,
    variant_label
) {
  design_objects <- create_design_and_contrasts(sample_metadata)
  cluster_ids <- natural_cluster_order(names(counts_by_cluster))

  per_test_results <- lapply(design_objects$test_definitions$test_id, function(x) list())
  names(per_test_results) <- design_objects$test_definitions$test_id

  for (cluster_index in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[[cluster_index]]

    message(
      "edgeR ", variant_label, " [", cluster_index, "/", length(cluster_ids), "] cluster ", cluster_id
    )

    counts_current <- counts_by_cluster[[cluster_id]][
      ,
      design_objects$metadata$sample_ID,
      drop = FALSE
    ]

    if (!identical(colnames(counts_current), rownames(design_objects$design))) {
      stop("Count/design order mismatch in cluster ", cluster_id, ".", call. = FALSE)
    }

    dge <- edgeR::DGEList(counts = as.matrix(counts_current))

    keep <- edgeR::filterByExpr(
      dge,
      design = design_objects$design,
      min.count = filter_by_expr_parameters$min.count,
      min.total.count = filter_by_expr_parameters$min.total.count,
      large.n = filter_by_expr_parameters$large.n,
      min.prop = filter_by_expr_parameters$min.prop
    )

    if (!any(keep)) {
      stop("No genes passed filterByExpr in ", variant_label, ", cluster ", cluster_id, ".", call. = FALSE)
    }

    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    dge <- edgeR::estimateDisp(dge, design = design_objects$design, robust = TRUE)
    fit <- edgeR::glmQLFit(dge, design = design_objects$design, robust = TRUE)

    for (test_id in design_objects$test_definitions$test_id) {
      qlf <- edgeR::glmQLFTest(
        fit,
        contrast = design_objects$contrasts[, test_id]
      )

      result <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table |>
        as.data.frame(stringsAsFactors = FALSE) |>
        tibble::rownames_to_column("ensembl_gene_id") |>
        tibble::as_tibble() |>
        dplyr::left_join(gene_annotation, by = "ensembl_gene_id") |>
        dplyr::transmute(
          cluster_id = as.character(cluster_id),
          ensembl_gene_id = .data$ensembl_gene_id,
          gene = .data$gene,
          chromosome = .data$chromosome,
          logFC = .data$logFC,
          logCPM = .data$logCPM,
          F = .data$F,
          PValue = .data$PValue,
          FDR = .data$FDR
        )

      per_test_results[[test_id]][[cluster_id]] <- result
    }
  }

  combined_results <- lapply(per_test_results, dplyr::bind_rows)

  for (test_id in names(combined_results)) {
    combined_results[[test_id]] <- combined_results[[test_id]] |>
      dplyr::mutate(cluster_id = factor(.data$cluster_id, levels = cluster_ids)) |>
      dplyr::arrange(.data$cluster_id, .data$PValue, dplyr::desc(abs(.data$logFC))) |>
      dplyr::mutate(cluster_id = as.character(.data$cluster_id))
  }

  write_edger_workbook(
    result_list = combined_results,
    test_definitions = design_objects$test_definitions,
    output_file = output_file
  )

  combined_results
}

# ============================================================================== 
# 5. Validate inputs
# ============================================================================== 

stage_message(1, "VALIDATING INPUTS")

required_files <- c(input_rdata_file, annotation_rds_file, crescendo_input_rds)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"), call. = FALSE)
}

if (!dir.exists(corrected_spot_chunks_dir)) {
  stop("Missing corrected spot-count chunk directory:\n", corrected_spot_chunks_dir, call. = FALSE)
}

message("SCRIPT BUILD: ", SCRIPT_BUILD)
message("Output root:")
message(comparison_root)

# ============================================================================== 
# 6. Load Crescendo object and define common gene universe
# ============================================================================== 

stage_message(2, "LOADING CRESCENDO OBJECT AND DEFINING COMMON GENE UNIVERSE")

crescendo_object <- readRDS(crescendo_input_rds)

required_crescendo_elements <- c(
  "sample_metadata",
  "gene_metadata",
  "raw_counts_by_cluster",
  "crescendo_corrected_counts_by_cluster"
)

missing_crescendo_elements <- setdiff(required_crescendo_elements, names(crescendo_object))

if (length(missing_crescendo_elements) > 0L) {
  stop(
    "Crescendo RDS missing element(s): ",
    paste(missing_crescendo_elements, collapse = ", "),
    call. = FALSE
  )
}

sample_metadata <- crescendo_object$sample_metadata |>
  tibble::as_tibble() |>
  dplyr::transmute(
    sample_ID = as.character(.data$sample_ID),
    fmt_donor_group = factor(as.character(.data$fmt_donor_group), levels = c("Neurotypical", "ASD")),
    sex = factor(as.character(.data$sex), levels = c("Male", "Female"))
  ) |>
  dplyr::distinct(.data$sample_ID, .keep_all = TRUE) |>
  dplyr::arrange(.data$sex, .data$fmt_donor_group, .data$sample_ID)

if (nrow(sample_metadata) != expected_number_of_samples) {
  stop(
    "Expected ", expected_number_of_samples, " samples but found ", nrow(sample_metadata), ".",
    call. = FALSE
  )
}

sample_order <- sample_metadata$sample_ID

# Crescendo gene_metadata contains gene symbols, not Ensembl IDs.
# Map the Crescendo symbols to the Ensembl v115 annotation.
gene_metadata_crescendo <- crescendo_object$gene_metadata |>
  tibble::as_tibble()

required_gene_metadata_columns <- c(
  "feature_id",
  "gene"
)

missing_gene_metadata_columns <- setdiff(
  required_gene_metadata_columns,
  colnames(gene_metadata_crescendo)
)

if (length(missing_gene_metadata_columns) > 0L) {
  stop(
    "Crescendo gene_metadata is missing column(s): ",
    paste(missing_gene_metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

gene_metadata_crescendo <- gene_metadata_crescendo |>
  dplyr::transmute(
    feature_id = as.character(.data$feature_id),
    gene = as.character(.data$gene)
  ) |>
  dplyr::filter(
    !is.na(.data$feature_id),
    .data$feature_id != "",
    !is.na(.data$gene),
    .data$gene != ""
  ) |>
  dplyr::distinct(
    .data$feature_id,
    .keep_all = TRUE
  )

annotation <- readRDS(annotation_rds_file) |>
  tibble::as_tibble()

required_annotation_columns <- c(
  "chromosome",
  "ensembl_gene_id",
  "gene",
  "gene_biotype"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  colnames(annotation)
)

if (length(missing_annotation_columns) > 0L) {
  stop(
    "Ensembl v115 annotation is missing column(s): ",
    paste(missing_annotation_columns, collapse = ", "),
    call. = FALSE
  )
}

annotation <- annotation |>
  dplyr::transmute(
    chromosome = as.character(.data$chromosome),
    ensembl_gene_id = sub(
      "\\.[0-9]+$",
      "",
      as.character(.data$ensembl_gene_id)
    ),
    gene = as.character(.data$gene),
    gene_biotype = as.character(.data$gene_biotype)
  ) |>
  dplyr::filter(
    .data$gene_biotype == "protein_coding",
    .data$chromosome %in% analysis_chromosomes,
    !is.na(.data$gene),
    .data$gene != "",
    !is.na(.data$ensembl_gene_id),
    .data$ensembl_gene_id != ""
  ) |>
  dplyr::distinct(
    .data$ensembl_gene_id,
    .keep_all = TRUE
  )

# Keep only unambiguous gene-symbol mappings in Ensembl v115.
annotation_unique_symbols <- annotation |>
  dplyr::add_count(
    .data$gene,
    name = "n_symbol_matches"
  ) |>
  dplyr::filter(
    .data$n_symbol_matches == 1L
  ) |>
  dplyr::select(
    -.data$n_symbol_matches
  )

gene_universe <- gene_metadata_crescendo |>
  dplyr::inner_join(
    annotation_unique_symbols |>
      dplyr::select(
        .data$ensembl_gene_id,
        .data$gene,
        .data$chromosome
      ),
    by = "gene"
  ) |>
  dplyr::transmute(
    feature_id = .data$feature_id,
    ensembl_gene_id = .data$ensembl_gene_id,
    gene = .data$gene,
    chromosome = .data$chromosome
  ) |>
  dplyr::distinct(
    .data$feature_id,
    .keep_all = TRUE
  )

if (nrow(gene_universe) == 0L) {
  stop(
    "No Crescendo gene symbols matched protein-coding Ensembl v115 genes on chr1-19/X.",
    call. = FALSE
  )
}

message(
  "Common Crescendo / protein-coding chr1-19,X gene universe: ",
  format(nrow(gene_universe), big.mark = ",")
)

comparison_feature_ids <- gene_universe$feature_id
comparison_ensembl_ids <- gene_universe$ensembl_gene_id

# ============================================================================== 
# 7. Prepare matched pseudobulk matrices
# ============================================================================== 

stage_message(3, "PREPARING MATCHED PSEUDOBULK MATRICES")

raw_source_by_cluster <- crescendo_object$raw_counts_by_cluster
crescendo_source_by_cluster <- crescendo_object$crescendo_corrected_counts_by_cluster

cluster_ids <- natural_cluster_order(
  intersect(names(raw_source_by_cluster), names(crescendo_source_by_cluster))
)

if (length(cluster_ids) != expected_number_of_clusters) {
  stop(
    "Expected ", expected_number_of_clusters, " common clusters but found ", length(cluster_ids), ".",
    call. = FALSE
  )
}

original_counts_by_cluster <- vector("list", length(cluster_ids))
crescendo_counts_by_cluster <- vector("list", length(cluster_ids))
names(original_counts_by_cluster) <- cluster_ids
names(crescendo_counts_by_cluster) <- cluster_ids

for (cluster_index in seq_along(cluster_ids)) {
  cluster_id <- cluster_ids[[cluster_index]]

  message("[", cluster_index, "/", length(cluster_ids), "] cluster ", cluster_id)

  raw_current <- raw_source_by_cluster[[cluster_id]]
  crescendo_current <- crescendo_source_by_cluster[[cluster_id]]

  missing_raw_features <- setdiff(comparison_feature_ids, rownames(raw_current))
  missing_crescendo_features <- setdiff(comparison_feature_ids, rownames(crescendo_current))

  if (length(missing_raw_features) > 0L || length(missing_crescendo_features) > 0L) {
    stop(
      "Common gene universe is not present in both pseudobulk matrices for cluster ",
      cluster_id,
      ".",
      call. = FALSE
    )
  }

  raw_current <- raw_current[comparison_feature_ids, sample_order, drop = FALSE]
  crescendo_current <- crescendo_current[comparison_feature_ids, sample_order, drop = FALSE]

  rownames(raw_current) <- comparison_ensembl_ids
  rownames(crescendo_current) <- comparison_ensembl_ids

  original_counts_by_cluster[[cluster_id]] <- raw_current
  crescendo_counts_by_cluster[[cluster_id]] <- crescendo_current
}

rm(raw_source_by_cluster, crescendo_source_by_cluster)

# ============================================================================== 
# 8. Save counts, log2+1, QN
# ============================================================================== 

stage_message(4, "WRITING PSEUDOBULK, LOG2+1 AND QN MATRICES")

write_variant_matrices <- function(counts_by_cluster, dirs, variant_name) {
  for (cluster_index in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[[cluster_index]]

    message(
      variant_name, " [", cluster_index, "/", length(cluster_ids), "] cluster ", cluster_id
    )

    counts_current <- counts_by_cluster[[cluster_id]]
    log2_current <- create_log2plus1(counts_current)
    qn_current <- create_quantile_normalized(log2_current)

    prefix <- paste0(
      dataset_name, "_", clustering_name, "_", variant_name, "_cluster_", cluster_id
    )

    write_tsv_gz(
      matrix_to_output_table(counts_current, gene_universe),
      file.path(dirs$pseudobulk, paste0(prefix, "_pseudobulkCounts.tsv.gz"))
    )

    write_tsv_gz(
      matrix_to_output_table(log2_current, gene_universe),
      file.path(dirs$log2plus1, paste0(prefix, "_log2plus1.tsv.gz"))
    )

    write_tsv_gz(
      matrix_to_output_table(qn_current, gene_universe),
      file.path(dirs$log2plus1_qn, paste0(prefix, "_log2plus1_QN.tsv.gz"))
    )

    rm(log2_current, qn_current)
    invisible(gc())
  }
}

write_variant_matrices(
  counts_by_cluster = original_counts_by_cluster,
  dirs = original_dirs,
  variant_name = "original"
)

write_variant_matrices(
  counts_by_cluster = crescendo_counts_by_cluster,
  dirs = crescendo_dirs,
  variant_name = "crescendo"
)

# ============================================================================== 
# 9. Raw spot counts and spot metadata
# ============================================================================== 

stage_message(5, "LOADING RAW SPOT COUNTS FOR DETECTION PERCENTAGES")

seurat_object <- load_seurat_object(input_rdata_file)
raw_spot_counts <- get_raw_counts(seurat_object, assay_name)
spot_metadata <- seurat_object[[]]

missing_spots_in_metadata <- setdiff(colnames(raw_spot_counts), rownames(spot_metadata))

if (length(missing_spots_in_metadata) > 0L) {
  stop("Some raw-count spots are missing from Seurat metadata.", call. = FALSE)
}

spot_metadata <- spot_metadata[colnames(raw_spot_counts), , drop = FALSE]
spot_metadata$sample_ID_analysis <- as.character(spot_metadata[[sample_column]])

keep_spots <- spot_metadata$sample_ID_analysis %in% sample_order
raw_spot_counts <- raw_spot_counts[, keep_spots, drop = FALSE]
spot_metadata <- spot_metadata[keep_spots, , drop = FALSE]

if (dplyr::n_distinct(spot_metadata$sample_ID_analysis) != expected_number_of_samples) {
  stop("Unexpected number of samples in retained spot metadata.", call. = FALSE)
}

missing_raw_features <- setdiff(comparison_feature_ids, rownames(raw_spot_counts))

if (length(missing_raw_features) > 0L) {
  stop(
    "Raw Seurat matrix is missing ", length(missing_raw_features), " common comparison genes.",
    call. = FALSE
  )
}

raw_spot_counts <- raw_spot_counts[comparison_feature_ids, , drop = FALSE]
rownames(raw_spot_counts) <- comparison_ensembl_ids

sample_spot_ids <- setNames(
  lapply(
    sample_order,
    function(sample_id) {
      rownames(spot_metadata)[spot_metadata$sample_ID_analysis == sample_id]
    }
  ),
  sample_order
)

# ============================================================================== 
# 10. ORIGINAL detection
# ============================================================================== 

stage_message(6, "CALCULATING ORIGINAL SPOT-LEVEL DETECTION")

original_detection <- matrix(
  NA_real_,
  nrow = length(comparison_ensembl_ids),
  ncol = length(sample_order),
  dimnames = list(comparison_ensembl_ids, sample_order)
)

for (sample_id in sample_order) {
  spots <- sample_spot_ids[[sample_id]]

  if (length(spots) == 0L) {
    stop("No spots found for sample ", sample_id, ".", call. = FALSE)
  }

  original_detection[, sample_id] <- 100 *
    as.numeric(Matrix::rowSums(raw_spot_counts[, spots, drop = FALSE] > 0)) /
    length(spots)
}

# ============================================================================== 
# 11. CRESCENDO detection from corrected spot chunks
# ============================================================================== 

stage_message(7, "CALCULATING CRESCENDO SPOT-LEVEL DETECTION")

corrected_chunk_files <- list.files(
  corrected_spot_chunks_dir,
  pattern = "_crescendoCorrected_spotCounts\\.rds$",
  full.names = TRUE
)

if (length(corrected_chunk_files) == 0L) {
  stop("No corrected spot-count chunks found.", call. = FALSE)
}

extract_chunk_number <- function(file_name) {
  as.integer(sub(".*_chunk_([0-9]+)_genes_.*", "\\1", basename(file_name)))
}

chunk_numbers <- vapply(corrected_chunk_files, extract_chunk_number, integer(1))
corrected_chunk_files <- corrected_chunk_files[order(chunk_numbers)]

crescendo_detection <- matrix(
  NA_real_,
  nrow = length(comparison_ensembl_ids),
  ncol = length(sample_order),
  dimnames = list(comparison_ensembl_ids, sample_order)
)

feature_to_ensembl <- setNames(gene_universe$ensembl_gene_id, gene_universe$feature_id)
all_retained_spot_ids <- rownames(spot_metadata)

for (chunk_index in seq_along(corrected_chunk_files)) {
  message("Corrected detection chunk ", chunk_index, "/", length(corrected_chunk_files))

  chunk_object <- readRDS(corrected_chunk_files[[chunk_index]])

  if (!"corrected_counts_gene_by_spot" %in% names(chunk_object)) {
    stop(
      "Chunk missing corrected_counts_gene_by_spot:\n",
      corrected_chunk_files[[chunk_index]],
      call. = FALSE
    )
  }

  corrected_chunk <- methods::as(
    chunk_object$corrected_counts_gene_by_spot,
    "dgCMatrix"
  )

  features_here <- intersect(rownames(corrected_chunk), comparison_feature_ids)

  if (length(features_here) == 0L) {
    rm(chunk_object, corrected_chunk)
    next
  }

  missing_spots <- setdiff(all_retained_spot_ids, colnames(corrected_chunk))

  if (length(missing_spots) > 0L) {
    stop(
      "Corrected chunk is missing ", length(missing_spots), " retained spots.",
      call. = FALSE
    )
  }

  corrected_chunk <- corrected_chunk[
    features_here,
    all_retained_spot_ids,
    drop = FALSE
  ]

  ensembl_here <- unname(feature_to_ensembl[features_here])

  for (sample_id in sample_order) {
    spots <- sample_spot_ids[[sample_id]]

    crescendo_detection[ensembl_here, sample_id] <- 100 *
      as.numeric(Matrix::rowSums(corrected_chunk[, spots, drop = FALSE] > 0)) /
      length(spots)
  }

  rm(chunk_object, corrected_chunk)
  invisible(gc())
}

if (anyNA(crescendo_detection)) {
  missing_genes <- rownames(crescendo_detection)[apply(is.na(crescendo_detection), 1, any)]
  stop(
    "Crescendo detection was not calculated for ", length(missing_genes), " genes.",
    call. = FALSE
  )
}

# ============================================================================== 
# 12. Write detection tables
# ============================================================================== 

stage_message(8, "WRITING DETECTION TABLES")

make_detection_table <- function(detection_matrix) {
  annotation_current <- gene_universe[
    match(rownames(detection_matrix), gene_universe$ensembl_gene_id),
    ,
    drop = FALSE
  ]

  dplyr::bind_cols(
    annotation_current |>
      dplyr::select(.data$ensembl_gene_id, .data$gene, .data$chromosome),
    tibble::tibble(
      mean_percent_positive = rowMeans(detection_matrix),
      sd_percent_positive = safe_row_sd(detection_matrix)
    ),
    tibble::as_tibble(
      as.data.frame(detection_matrix, check.names = FALSE),
      .name_repair = "minimal"
    )
  ) |>
    dplyr::arrange(dplyr::desc(.data$mean_percent_positive), .data$sd_percent_positive)
}

original_detection_table <- make_detection_table(original_detection)
crescendo_detection_table <- make_detection_table(crescendo_detection)

write_tsv_gz(
  original_detection_table,
  file.path(
    original_dirs$detection,
    paste0(dataset_name, "_original_spotDetectionPercent_all16samples_meanSD.tsv.gz")
  )
)

write_tsv_gz(
  crescendo_detection_table,
  file.path(
    crescendo_dirs$detection,
    paste0(dataset_name, "_crescendo_spotDetectionPercent_all16samples_meanSD.tsv.gz")
  )
)

# ============================================================================== 
# 13. edgeR ORIGINAL
# ============================================================================== 

stage_message(9, "RUNNING edgeR — ORIGINAL")

original_edger_xlsx <- file.path(
  original_dirs$edger,
  paste0(dataset_name, "_original_edgeR_allGenes_allClusters_nineTests.xlsx")
)

original_edger_results <- run_edger_all_clusters(
  counts_by_cluster = original_counts_by_cluster,
  sample_metadata = sample_metadata,
  gene_annotation = gene_universe |>
    dplyr::select(.data$ensembl_gene_id, .data$gene, .data$chromosome),
  output_file = original_edger_xlsx,
  variant_label = "ORIGINAL"
)

# ============================================================================== 
# 14. edgeR CRESCENDO
# ============================================================================== 

stage_message(10, "RUNNING edgeR — CRESCENDO")

crescendo_edger_xlsx <- file.path(
  crescendo_dirs$edger,
  paste0(dataset_name, "_crescendo_edgeR_allGenes_allClusters_nineTests.xlsx")
)

crescendo_edger_results <- run_edger_all_clusters(
  counts_by_cluster = crescendo_counts_by_cluster,
  sample_metadata = sample_metadata,
  gene_annotation = gene_universe |>
    dplyr::select(.data$ensembl_gene_id, .data$gene, .data$chromosome),
  output_file = crescendo_edger_xlsx,
  variant_label = "CRESCENDO"
)

# ============================================================================== 
# 15. Final report
# ============================================================================== 

stage_message(11, "FINISHED")

message("Samples: ", nrow(sample_metadata))
message("Clusters: ", length(cluster_ids))
message(
  "Common protein-coding chr1-19/X genes: ",
  format(nrow(gene_universe), big.mark = ",")
)

message("\nORIGINAL output:")
message(original_root)

message("\nCRESCENDO output:")
message(crescendo_comparison_root)

message("\nORIGINAL edgeR XLSX:")
message(original_edger_xlsx)

message("\nCRESCENDO edgeR XLSX:")
message(crescendo_edger_xlsx)

message("\nComparison root:")
message(comparison_root)
