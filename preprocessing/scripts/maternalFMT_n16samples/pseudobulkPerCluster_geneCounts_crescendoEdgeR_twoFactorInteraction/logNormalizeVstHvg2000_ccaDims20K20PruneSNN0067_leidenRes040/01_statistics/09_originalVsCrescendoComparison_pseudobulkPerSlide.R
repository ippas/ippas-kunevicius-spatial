#!/usr/bin/env Rscript
# ==============================================================================
# 03_originalVsCrescendoComparison_pseudobulkPerTissue.R
#
# PURPOSE
#   Prepare a clean whole-tissue / whole-section comparison of ORIGINAL raw
#   counts versus CRESCENDO-corrected counts.
#
#   There is NO Leiden-cluster split in this script.
#
# ANALYSIS UNIT
#   One biological sample = one whole tissue section.
#   All spots belonging to a sample are summed together.
#
# COMMON GENE UNIVERSE
#   Crescendo gene symbols
#     -> map to Ensembl v115
#     -> protein-coding only
#     -> chromosomes 1-19 and X
#     -> require presence in ORIGINAL raw spot counts
#     -> require presence in CRESCENDO-corrected spot-count chunks
#
# FOR EACH DATA VARIANT
#
#   01_original
#     ORIGINAL raw spot counts
#
#   02_crescendo
#     CRESCENDO-corrected spot counts
#
# PREPARED OUTPUTS
#
#   1. pseudobulk counts:
#        genes x 16 biological samples
#
#   2. log2(count + 1)
#
#   3. log2(count + 1) + quantile normalization across the 16 samples
#
#   4. spot-level detection percentage:
#        % spots with count > 0 for every gene and every biological sample
#        + mean across samples
#        + SD across samples
#
#   5. edgeR:
#        expression ~ fmt_donor_group * sex
#
#      Tests:
#        Interaction
#        Overall Sex
#        Overall Group
#        ASD Male vs Neurotypical Male
#        ASD Female vs Neurotypical Female
#        Neurotypical Female vs Neurotypical Male
#        ASD Female vs ASD Male
#        ASD Female vs Neurotypical Male
#        ASD Male vs Neurotypical Female
#
#      Two XLSX files are produced:
#        ORIGINAL
#        CRESCENDO
#
#      All genes retained by filterByExpr are written.
#      No result filtering by p-value, FDR or log2FC is performed.
#
# IMPORTANT
#   - all 16 biological samples are used
#   - edgeR uses count-scale pseudobulk values only
#   - log2+1 and QN are descriptive inputs for later heatmaps
#   - ORIGINAL and CRESCENDO use the same sample order and same gene universe
# ==============================================================================

options(stringsAsFactors = FALSE)

SCRIPT_BUILD <- "20260903_V1_PSEUDOBULK_PER_TISSUE"

# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "dplyr",
  "tibble",
  "edgeR",
  "limma",
  "preprocessCore",
  "openxlsx"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
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
# 2. Configuration
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"

dataset_name <- "maternalFMT_n16samples"

clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

assay_name <- "RNA"

sample_column <- "sample_ID"

expected_number_of_samples <- 16L

analysis_chromosomes <- c(
  as.character(1:19),
  "X"
)

filter_by_expr_parameters <- list(
  min.count = 10,
  min.total.count = 15,
  large.n = 10,
  min.prop = 0.7
)

# ==============================================================================
# 3. Input paths
# ==============================================================================

seurat_rdata_file <- file.path(
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
    "02_",
    dataset_name,
    "_logNormalizeVst_hvg2000_",
    "ccaDims20K20PruneSNN0067_",
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
# 4. Output structure
# ==============================================================================

comparison_root <- file.path(
  crescendo_statistics_root,
  "03_originalVsCrescendoComparison_pseudobulkPerTissue"
)

original_root <- file.path(
  comparison_root,
  "01_original"
)

crescendo_comparison_root <- file.path(
  comparison_root,
  "02_crescendo"
)

make_variant_dirs <- function(root_dir) {
  dirs <- list(
    pseudobulk = file.path(
      root_dir,
      "01_pseudobulkCounts"
    ),
    log2plus1 = file.path(
      root_dir,
      "02_log2plus1"
    ),
    log2plus1_qn = file.path(
      root_dir,
      "03_log2plus1_QN"
    ),
    detection = file.path(
      root_dir,
      "04_detectionPercent"
    ),
    edger = file.path(
      root_dir,
      "05_edgeR"
    )
  )

  for (dir_current in c(
    root_dir,
    unname(dirs)
  )) {
    dir.create(
      dir_current,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  dirs
}

original_dirs <- make_variant_dirs(
  original_root
)

crescendo_dirs <- make_variant_dirs(
  crescendo_comparison_root
)

# ==============================================================================
# 5. Helper functions
# ==============================================================================

stage_message <- function(
    stage_number,
    text
) {
  message("")
  message("============================================================")
  message("STAGE ", stage_number, ": ", text)
  message("============================================================")
}

write_tsv_gz <- function(
    x,
    output_file
) {
  connection <- gzfile(
    output_file,
    open = "wt",
    encoding = "UTF-8"
  )

  on.exit(
    close(connection),
    add = TRUE
  )

  utils::write.table(
    x,
    file = connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )

  invisible(output_file)
}

load_seurat_object <- function(
    input_file
) {
  env <- new.env(
    parent = emptyenv()
  )

  loaded_objects <- load(
    input_file,
    envir = env
  )

  seurat_object_names <- loaded_objects[
    vapply(
      loaded_objects,
      function(object_name) {
        inherits(
          get(
            object_name,
            envir = env
          ),
          "Seurat"
        )
      },
      logical(1)
    )
  ]

  if (length(seurat_object_names) != 1L) {
    stop(
      "Expected exactly one Seurat object in RData. Found: ",
      paste(
        seurat_object_names,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  get(
    seurat_object_names[[1]],
    envir = env
  )
}

get_raw_counts <- function(
    seurat_object,
    assay_name
) {
  assay_object <- seurat_object[[assay_name]]

  if (inherits(
    assay_object,
    "Assay5"
  )) {
    assay_layers <- SeuratObject::Layers(
      assay_object
    )

    count_layers <- assay_layers[
      assay_layers == "counts" |
        startsWith(
          assay_layers,
          "counts."
        )
    ]

    if (length(count_layers) == 0L) {
      stop(
        "No raw count layer found in assay ",
        assay_name,
        ".",
        call. = FALSE
      )
    }

    if (length(count_layers) > 1L) {
      assay_object <- SeuratObject::JoinLayers(
        assay_object
      )
    }

    counts <- SeuratObject::LayerData(
      assay_object,
      layer = "counts"
    )
  } else {
    counts <- Seurat::GetAssayData(
      seurat_object,
      assay = assay_name,
      slot = "counts"
    )
  }

  methods::as(
    counts,
    "dgCMatrix"
  )
}

calculate_row_sd <- function(
    x
) {
  if (ncol(x) <= 1L) {
    return(
      rep(
        NA_real_,
        nrow(x)
      )
    )
  }

  apply(
    x,
    1L,
    stats::sd
  )
}

create_log2plus1 <- function(
    count_matrix
) {
  x <- as.matrix(
    count_matrix
  )

  storage.mode(x) <- "double"

  log2(
    x + 1
  )
}

create_quantile_normalized <- function(
    log2_matrix
) {
  x <- as.matrix(
    log2_matrix
  )

  storage.mode(x) <- "double"

  qn <- preprocessCore::normalize.quantiles(
    x,
    copy = TRUE
  )

  rownames(qn) <- rownames(x)
  colnames(qn) <- colnames(x)

  qn
}

matrix_to_output_table <- function(
    matrix_object,
    gene_universe
) {
  gene_ids <- rownames(
    matrix_object
  )

  annotation_current <- gene_universe[
    match(
      gene_ids,
      gene_universe$ensembl_gene_id
    ),
    ,
    drop = FALSE
  ]

  if (anyNA(
    annotation_current$ensembl_gene_id
  )) {
    stop(
      "Failed to annotate one or more matrix genes.",
      call. = FALSE
    )
  }

  dplyr::bind_cols(
    annotation_current |>
      dplyr::select(
        "ensembl_gene_id",
        "gene",
        "chromosome"
      ),
    tibble::as_tibble(
      as.data.frame(
        matrix_object,
        check.names = FALSE
      ),
      .name_repair = "minimal"
    )
  )
}

create_design_and_contrasts <- function(
    sample_metadata
) {
  model_metadata <- sample_metadata |>
    dplyr::transmute(
      sample_ID =
        as.character(.data$sample_ID),
      fmt_donor_group = factor(
        as.character(
          .data$fmt_donor_group
        ),
        levels = c(
          "Neurotypical",
          "ASD"
        )
      ),
      sex = factor(
        as.character(
          .data$sex
        ),
        levels = c(
          "Male",
          "Female"
        )
      )
    )

  design <- stats::model.matrix(
    ~ fmt_donor_group * sex,
    data = model_metadata
  )

  rownames(design) <-
    model_metadata$sample_ID

  expected_design_columns <- c(
    "(Intercept)",
    "fmt_donor_groupASD",
    "sexFemale",
    "fmt_donor_groupASD:sexFemale"
  )

  if (!identical(
    colnames(design),
    expected_design_columns
  )) {
    stop(
      "Unexpected design columns: ",
      paste(
        colnames(design),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  colnames(design) <- c(
    "Intercept",
    "groupASD",
    "sexFemale",
    "groupASD_sexFemale"
  )

  contrasts <- limma::makeContrasts(
    Interaction =
      groupASD_sexFemale,

    Overall_Sex_Female_vs_Male =
      sexFemale +
      0.5 * groupASD_sexFemale,

    Overall_Group_ASD_vs_Neurotypical =
      groupASD +
      0.5 * groupASD_sexFemale,

    ASD_Male_vs_Neurotypical_Male =
      groupASD,

    ASD_Female_vs_Neurotypical_Female =
      groupASD +
      groupASD_sexFemale,

    Neurotypical_Female_vs_Neurotypical_Male =
      sexFemale,

    ASD_Female_vs_ASD_Male =
      sexFemale +
      groupASD_sexFemale,

    ASD_Female_vs_Neurotypical_Male =
      groupASD +
      sexFemale +
      groupASD_sexFemale,

    ASD_Male_vs_Neurotypical_Female =
      groupASD -
      sexFemale,

    levels = design
  )

  test_definitions <- tibble::tribble(
    ~test_id,
    ~sheet_name,
    ~comparison,

    "Interaction",
    "Interaction",
    "(ASD - Neurotypical)_Female - (ASD - Neurotypical)_Male",

    "Overall_Sex_Female_vs_Male",
    "Overall_Sex",
    "Female vs Male averaged equally across donor groups",

    "Overall_Group_ASD_vs_Neurotypical",
    "Overall_Group",
    "ASD vs Neurotypical averaged equally across sexes",

    "ASD_Male_vs_Neurotypical_Male",
    "Group_in_Male",
    "ASD Male vs Neurotypical Male",

    "ASD_Female_vs_Neurotypical_Female",
    "Group_in_Female",
    "ASD Female vs Neurotypical Female",

    "Neurotypical_Female_vs_Neurotypical_Male",
    "Sex_in_Neurotypical",
    "Neurotypical Female vs Neurotypical Male",

    "ASD_Female_vs_ASD_Male",
    "Sex_in_ASD",
    "ASD Female vs ASD Male",

    "ASD_Female_vs_Neurotypical_Male",
    "ASD_F_vs_NT_M",
    "ASD Female vs Neurotypical Male",

    "ASD_Male_vs_Neurotypical_Female",
    "ASD_M_vs_NT_F",
    "ASD Male vs Neurotypical Female"
  )

  list(
    metadata = model_metadata,
    design = design,
    contrasts = contrasts,
    test_definitions =
      test_definitions
  )
}

run_edger <- function(
    count_matrix,
    sample_metadata,
    gene_universe,
    output_file,
    variant_label
) {
  stage_message(
    "",
    paste0(
      "edgeR — ",
      variant_label
    )
  )

  design_objects <-
    create_design_and_contrasts(
      sample_metadata
    )

  count_matrix <- count_matrix[
    ,
    design_objects$metadata$sample_ID,
    drop = FALSE
  ]

  if (!identical(
    colnames(count_matrix),
    rownames(
      design_objects$design
    )
  )) {
    stop(
      "Count/design sample order mismatch for ",
      variant_label,
      ".",
      call. = FALSE
    )
  }

  dge <- edgeR::DGEList(
    counts = as.matrix(
      count_matrix
    )
  )

  keep <- edgeR::filterByExpr(
    dge,
    design =
      design_objects$design,
    min.count =
      filter_by_expr_parameters$min.count,
    min.total.count =
      filter_by_expr_parameters$min.total.count,
    large.n =
      filter_by_expr_parameters$large.n,
    min.prop =
      filter_by_expr_parameters$min.prop
  )

  message(
    variant_label,
    " genes before filterByExpr: ",
    nrow(dge)
  )

  message(
    variant_label,
    " genes after filterByExpr: ",
    sum(keep)
  )

  if (!any(keep)) {
    stop(
      "No genes passed filterByExpr for ",
      variant_label,
      ".",
      call. = FALSE
    )
  }

  dge <- dge[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]

  dge <- edgeR::calcNormFactors(
    dge,
    method = "TMM"
  )

  dge <- edgeR::estimateDisp(
    dge,
    design =
      design_objects$design,
    robust = TRUE
  )

  fit <- edgeR::glmQLFit(
    dge,
    design =
      design_objects$design,
    robust = TRUE
  )

  result_list <- vector(
    "list",
    nrow(
      design_objects$test_definitions
    )
  )

  names(result_list) <-
    design_objects$test_definitions$test_id

  for (
    test_index in seq_len(
      nrow(
        design_objects$test_definitions
      )
    )
  ) {
    test_id <-
      design_objects$test_definitions$test_id[[test_index]]

    message(
      variant_label,
      " edgeR test ",
      test_index,
      "/",
      nrow(
        design_objects$test_definitions
      ),
      ": ",
      test_id
    )

    qlf <- edgeR::glmQLFTest(
      fit,
      contrast =
        design_objects$contrasts[
          ,
          test_id
        ]
    )

    result <- edgeR::topTags(
      qlf,
      n = Inf,
      sort.by = "PValue"
    )$table |>
      as.data.frame(
        stringsAsFactors = FALSE
      ) |>
      tibble::rownames_to_column(
        "ensembl_gene_id"
      ) |>
      tibble::as_tibble() |>
      dplyr::left_join(
        gene_universe |>
          dplyr::select(
            "ensembl_gene_id",
            "gene",
            "chromosome"
          ),
        by = "ensembl_gene_id"
      ) |>
      dplyr::transmute(
        ensembl_gene_id =
          .data$ensembl_gene_id,
        gene =
          .data$gene,
        chromosome =
          .data$chromosome,
        logFC =
          .data$logFC,
        logCPM =
          .data$logCPM,
        F =
          .data$F,
        PValue =
          .data$PValue,
        FDR =
          .data$FDR
      )

    result_list[[test_id]] <-
      result
  }

  workbook <- openxlsx::createWorkbook()

  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#D9EAF7",
    border = "Bottom"
  )

  scientific_style <- openxlsx::createStyle(
    numFmt = "0.00000E+00"
  )

  decimal_style <- openxlsx::createStyle(
    numFmt = "0.00000"
  )

  for (
    test_index in seq_len(
      nrow(
        design_objects$test_definitions
      )
    )
  ) {
    test_id <-
      design_objects$test_definitions$test_id[[test_index]]

    sheet_name <-
      design_objects$test_definitions$sheet_name[[test_index]]

    result <- result_list[[test_id]]

    openxlsx::addWorksheet(
      workbook,
      sheetName = sheet_name
    )

    openxlsx::writeData(
      workbook,
      sheet = sheet_name,
      x = result,
      headerStyle = header_style,
      withFilter = TRUE
    )

    openxlsx::freezePane(
      workbook,
      sheet = sheet_name,
      firstActiveRow = 2,
      firstActiveCol = 4
    )

    numeric_columns <- which(
      vapply(
        result,
        is.numeric,
        logical(1)
      )
    )

    scientific_columns <- which(
      colnames(result) %in%
        c(
          "PValue",
          "FDR"
        )
    )

    decimal_columns <- setdiff(
      numeric_columns,
      scientific_columns
    )

    if (
      nrow(result) > 0L &&
      length(decimal_columns) > 0L
    ) {
      openxlsx::addStyle(
        workbook,
        sheet = sheet_name,
        style = decimal_style,
        rows = 2:(nrow(result) + 1L),
        cols = decimal_columns,
        gridExpand = TRUE,
        stack = TRUE
      )
    }

    if (
      nrow(result) > 0L &&
      length(scientific_columns) > 0L
    ) {
      openxlsx::addStyle(
        workbook,
        sheet = sheet_name,
        style = scientific_style,
        rows = 2:(nrow(result) + 1L),
        cols = scientific_columns,
        gridExpand = TRUE,
        stack = TRUE
      )
    }

    openxlsx::setColWidths(
      workbook,
      sheet = sheet_name,
      cols = seq_len(
        ncol(result)
      ),
      widths = "auto"
    )
  }

  openxlsx::saveWorkbook(
    workbook,
    output_file,
    overwrite = TRUE
  )

  list(
    results = result_list,
    keep = keep,
    dge = dge,
    fit = fit,
    design =
      design_objects$design,
    test_definitions =
      design_objects$test_definitions
  )
}

# ==============================================================================
# 6. Validate inputs
# ==============================================================================

stage_message(
  1,
  "VALIDATING INPUTS"
)

required_files <- c(
  seurat_rdata_file,
  annotation_rds_file,
  crescendo_input_rds
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (length(missing_files) > 0L) {
  stop(
    "Missing required input file(s):\n",
    paste(
      missing_files,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

if (!dir.exists(
  corrected_spot_chunks_dir
)) {
  stop(
    "Missing corrected spot-count chunk directory:\n",
    corrected_spot_chunks_dir,
    call. = FALSE
  )
}

message(
  "SCRIPT BUILD: ",
  SCRIPT_BUILD
)

message(
  "Output root:"
)

message(
  comparison_root
)

# ==============================================================================
# 7. Load Crescendo metadata and Ensembl annotation
# ==============================================================================

stage_message(
  2,
  "LOADING CRESCENDO METADATA AND ENSEMBL V115"
)

crescendo_object <- readRDS(
  crescendo_input_rds
)

if (!"sample_metadata" %in%
    names(crescendo_object)) {
  stop(
    "Crescendo RDS does not contain sample_metadata.",
    call. = FALSE
  )
}

if (!"gene_metadata" %in%
    names(crescendo_object)) {
  stop(
    "Crescendo RDS does not contain gene_metadata.",
    call. = FALSE
  )
}

sample_metadata <- crescendo_object$sample_metadata |>
  tibble::as_tibble() |>
  dplyr::transmute(
    sample_ID =
      as.character(.data$sample_ID),
    fmt_donor_group = factor(
      as.character(
        .data$fmt_donor_group
      ),
      levels = c(
        "Neurotypical",
        "ASD"
      )
    ),
    sex = factor(
      as.character(
        .data$sex
      ),
      levels = c(
        "Male",
        "Female"
      )
    )
  ) |>
  dplyr::distinct(
    .data$sample_ID,
    .keep_all = TRUE
  ) |>
  dplyr::arrange(
    .data$sex,
    .data$fmt_donor_group,
    .data$sample_ID
  )

if (
  nrow(sample_metadata) !=
    expected_number_of_samples
) {
  stop(
    "Expected ",
    expected_number_of_samples,
    " biological samples but found ",
    nrow(sample_metadata),
    ".",
    call. = FALSE
  )
}

sample_order <-
  sample_metadata$sample_ID

gene_metadata_crescendo <-
  crescendo_object$gene_metadata |>
  tibble::as_tibble()

required_gene_metadata_columns <- c(
  "feature_id",
  "gene"
)

missing_gene_metadata_columns <- setdiff(
  required_gene_metadata_columns,
  colnames(
    gene_metadata_crescendo
  )
)

if (
  length(
    missing_gene_metadata_columns
  ) > 0L
) {
  stop(
    "Crescendo gene_metadata is missing column(s): ",
    paste(
      missing_gene_metadata_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

gene_metadata_crescendo <-
  gene_metadata_crescendo |>
  dplyr::transmute(
    feature_id =
      as.character(.data$feature_id),
    gene =
      as.character(.data$gene)
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

annotation <- readRDS(
  annotation_rds_file
) |>
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

if (
  length(
    missing_annotation_columns
  ) > 0L
) {
  stop(
    "Ensembl v115 annotation is missing column(s): ",
    paste(
      missing_annotation_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

annotation <- annotation |>
  dplyr::transmute(
    chromosome =
      as.character(.data$chromosome),
    ensembl_gene_id = sub(
      "\\.[0-9]+$",
      "",
      as.character(
        .data$ensembl_gene_id
      )
    ),
    gene =
      as.character(.data$gene),
    gene_biotype =
      as.character(
        .data$gene_biotype
      )
  ) |>
  dplyr::filter(
    .data$gene_biotype ==
      "protein_coding",
    .data$chromosome %in%
      analysis_chromosomes,
    !is.na(.data$gene),
    .data$gene != "",
    !is.na(
      .data$ensembl_gene_id
    ),
    .data$ensembl_gene_id != ""
  ) |>
  dplyr::distinct(
    .data$ensembl_gene_id,
    .keep_all = TRUE
  )

annotation_unique_symbols <- annotation |>
  dplyr::add_count(
    .data$gene,
    name = "n_symbol_matches"
  ) |>
  dplyr::filter(
    .data$n_symbol_matches == 1L
  ) |>
  dplyr::select(
    -"n_symbol_matches"
  )

gene_universe <- gene_metadata_crescendo |>
  dplyr::inner_join(
    annotation_unique_symbols |>
      dplyr::select(
        "ensembl_gene_id",
        "gene",
        "chromosome"
      ),
    by = "gene"
  ) |>
  dplyr::transmute(
    feature_id =
      .data$feature_id,
    ensembl_gene_id =
      .data$ensembl_gene_id,
    gene =
      .data$gene,
    chromosome =
      .data$chromosome
  ) |>
  dplyr::distinct(
    .data$feature_id,
    .keep_all = TRUE
  )

message(
  "Mapped Crescendo genes -> protein-coding chr1-19/X: ",
  format(
    nrow(gene_universe),
    big.mark = ","
  )
)

# ==============================================================================
# 8. Load ORIGINAL raw spot counts
# ==============================================================================

stage_message(
  3,
  "LOADING ORIGINAL RAW SPOT COUNTS"
)

seurat_object <- load_seurat_object(
  seurat_rdata_file
)

raw_spot_counts <- get_raw_counts(
  seurat_object,
  assay_name
)

spot_metadata <- seurat_object[[]]

if (!sample_column %in%
    colnames(spot_metadata)) {
  stop(
    "Sample column `",
    sample_column,
    "` is missing from Seurat metadata.",
    call. = FALSE
  )
}

missing_spots_in_metadata <- setdiff(
  colnames(raw_spot_counts),
  rownames(spot_metadata)
)

if (
  length(
    missing_spots_in_metadata
  ) > 0L
) {
  stop(
    "Some raw-count spots are missing from Seurat metadata.",
    call. = FALSE
  )
}

spot_metadata <- spot_metadata[
  colnames(raw_spot_counts),
  ,
  drop = FALSE
]

spot_metadata$sample_ID_analysis <-
  as.character(
    spot_metadata[[sample_column]]
  )

keep_spots <-
  spot_metadata$sample_ID_analysis %in%
  sample_order

raw_spot_counts <- raw_spot_counts[
  ,
  keep_spots,
  drop = FALSE
]

spot_metadata <- spot_metadata[
  keep_spots,
  ,
  drop = FALSE
]

if (
  dplyr::n_distinct(
    spot_metadata$sample_ID_analysis
  ) != expected_number_of_samples
) {
  stop(
    "Retained spot metadata does not contain exactly 16 samples.",
    call. = FALSE
  )
}

# Require ORIGINAL presence.
gene_universe <- gene_universe |>
  dplyr::filter(
    .data$feature_id %in%
      rownames(raw_spot_counts)
  )

message(
  "After requiring presence in ORIGINAL raw counts: ",
  format(
    nrow(gene_universe),
    big.mark = ","
  )
)

if (nrow(gene_universe) == 0L) {
  stop(
    "No common genes remain.",
    call. = FALSE
  )
}

sample_spot_ids <- setNames(
  lapply(
    sample_order,
    function(sample_id) {
      rownames(spot_metadata)[
        spot_metadata$sample_ID_analysis ==
          sample_id
      ]
    }
  ),
  sample_order
)

sample_spot_counts <- vapply(
  sample_spot_ids,
  length,
  integer(1)
)

if (any(
  sample_spot_counts == 0L
)) {
  stop(
    "At least one biological sample contains zero spots.",
    call. = FALSE
  )
}

message(
  "Retained spots: ",
  format(
    ncol(raw_spot_counts),
    big.mark = ","
  )
)

message(
  "Spots per sample:"
)

print(
  sample_spot_counts
)

# ==============================================================================
# 9. ORIGINAL: pseudobulk and detection
# ==============================================================================

stage_message(
  4,
  "ORIGINAL — PSEUDOBULK AND DETECTION"
)

original_feature_ids <-
  gene_universe$feature_id

original_ensembl_ids <-
  gene_universe$ensembl_gene_id

raw_spot_counts_common <-
  raw_spot_counts[
    original_feature_ids,
    ,
    drop = FALSE
  ]

rownames(raw_spot_counts_common) <-
  original_ensembl_ids

original_pseudobulk <- matrix(
  0,
  nrow = nrow(
    gene_universe
  ),
  ncol = length(
    sample_order
  ),
  dimnames = list(
    original_ensembl_ids,
    sample_order
  )
)

original_detection <- matrix(
  NA_real_,
  nrow = nrow(
    gene_universe
  ),
  ncol = length(
    sample_order
  ),
  dimnames = list(
    original_ensembl_ids,
    sample_order
  )
)

for (
  sample_index in seq_along(
    sample_order
  )
) {
  sample_id <-
    sample_order[[sample_index]]

  message(
    "ORIGINAL [",
    sample_index,
    "/",
    length(sample_order),
    "] ",
    sample_id
  )

  spots <-
    sample_spot_ids[[sample_id]]

  counts_current <-
    raw_spot_counts_common[
      ,
      spots,
      drop = FALSE
    ]

  original_pseudobulk[
    ,
    sample_id
  ] <- as.numeric(
    Matrix::rowSums(
      counts_current
    )
  )

  original_detection[
    ,
    sample_id
  ] <- 100 *
    as.numeric(
      Matrix::rowSums(
        counts_current > 0
      )
    ) /
    length(spots)
}

# ==============================================================================
# 10. CRESCENDO: read corrected chunks, pseudobulk and detection
# ==============================================================================

stage_message(
  5,
  "CRESCENDO — PSEUDOBULK AND DETECTION FROM CORRECTED SPOT CHUNKS"
)

corrected_chunk_files <- list.files(
  corrected_spot_chunks_dir,
  pattern =
    "_crescendoCorrected_spotCounts\\.rds$",
  full.names = TRUE
)

if (
  length(
    corrected_chunk_files
  ) == 0L
) {
  stop(
    "No corrected spot-count chunks found in:\n",
    corrected_spot_chunks_dir,
    call. = FALSE
  )
}

extract_chunk_number <- function(
    file_name
) {
  x <- sub(
    ".*_chunk_([0-9]+)_genes_.*",
    "\\1",
    basename(file_name)
  )

  suppressWarnings(
    as.integer(x)
  )
}

chunk_numbers <- vapply(
  corrected_chunk_files,
  extract_chunk_number,
  integer(1)
)

if (anyNA(chunk_numbers)) {
  corrected_chunk_files <- sort(
    corrected_chunk_files
  )
} else {
  corrected_chunk_files <-
    corrected_chunk_files[
      order(chunk_numbers)
    ]
}

crescendo_pseudobulk <- matrix(
  NA_real_,
  nrow = nrow(
    gene_universe
  ),
  ncol = length(
    sample_order
  ),
  dimnames = list(
    gene_universe$ensembl_gene_id,
    sample_order
  )
)

crescendo_detection <- matrix(
  NA_real_,
  nrow = nrow(
    gene_universe
  ),
  ncol = length(
    sample_order
  ),
  dimnames = list(
    gene_universe$ensembl_gene_id,
    sample_order
  )
)

feature_to_ensembl <- setNames(
  gene_universe$ensembl_gene_id,
  gene_universe$feature_id
)

all_retained_spots <-
  rownames(spot_metadata)

seen_feature_ids <- character()

for (
  chunk_index in seq_along(
    corrected_chunk_files
  )
) {
  message(
    "CRESCENDO chunk ",
    chunk_index,
    "/",
    length(
      corrected_chunk_files
    ),
    ": ",
    basename(
      corrected_chunk_files[[chunk_index]]
    )
  )

  chunk_object <- readRDS(
    corrected_chunk_files[[chunk_index]]
  )

  if (
    !"corrected_counts_gene_by_spot" %in%
      names(chunk_object)
  ) {
    stop(
      "Corrected chunk does not contain `corrected_counts_gene_by_spot`:\n",
      corrected_chunk_files[[chunk_index]],
      call. = FALSE
    )
  }

  corrected_chunk <-
    chunk_object$corrected_counts_gene_by_spot

  corrected_chunk <- methods::as(
    corrected_chunk,
    "dgCMatrix"
  )

  features_here <- intersect(
    rownames(corrected_chunk),
    gene_universe$feature_id
  )

  if (
    length(
      features_here
    ) == 0L
  ) {
    rm(
      chunk_object,
      corrected_chunk
    )

    next
  }

  duplicated_seen <- intersect(
    features_here,
    seen_feature_ids
  )

  if (
    length(
      duplicated_seen
    ) > 0L
  ) {
    stop(
      "The same corrected feature occurred in multiple chunks. Example: ",
      duplicated_seen[[1]],
      call. = FALSE
    )
  }

  seen_feature_ids <- c(
    seen_feature_ids,
    features_here
  )

  missing_spots <- setdiff(
    all_retained_spots,
    colnames(
      corrected_chunk
    )
  )

  if (
    length(
      missing_spots
    ) > 0L
  ) {
    stop(
      "Corrected chunk is missing ",
      length(missing_spots),
      " retained spots.",
      call. = FALSE
    )
  }

  corrected_chunk <-
    corrected_chunk[
      features_here,
      all_retained_spots,
      drop = FALSE
    ]

  ensembl_here <- unname(
    feature_to_ensembl[
      features_here
    ]
  )

  for (
    sample_id in sample_order
  ) {
    spots <-
      sample_spot_ids[[sample_id]]

    counts_current <-
      corrected_chunk[
        ,
        spots,
        drop = FALSE
      ]

    crescendo_pseudobulk[
      ensembl_here,
      sample_id
    ] <- as.numeric(
      Matrix::rowSums(
        counts_current
      )
    )

    crescendo_detection[
      ensembl_here,
      sample_id
    ] <- 100 *
      as.numeric(
        Matrix::rowSums(
          counts_current > 0
        )
      ) /
      length(spots)
  }

  rm(
    chunk_object,
    corrected_chunk
  )

  invisible(
    gc()
  )
}

missing_corrected_features <- setdiff(
  gene_universe$feature_id,
  seen_feature_ids
)

if (
  length(
    missing_corrected_features
  ) > 0L
) {
  stop(
    "CRESCENDO chunks are missing ",
    length(missing_corrected_features),
    " genes from the common universe. Example: ",
    missing_corrected_features[[1]],
    call. = FALSE
  )
}

if (
  anyNA(
    crescendo_pseudobulk
  )
) {
  stop(
    "NA values remain in CRESCENDO pseudobulk matrix.",
    call. = FALSE
  )
}

if (
  anyNA(
    crescendo_detection
  )
) {
  stop(
    "NA values remain in CRESCENDO detection matrix.",
    call. = FALSE
  )
}

# ==============================================================================
# 11. Prepare log2+1 and QN matrices
# ==============================================================================

stage_message(
  6,
  "PREPARING LOG2+1 AND QUANTILE-NORMALIZED MATRICES"
)

original_log2plus1 <-
  create_log2plus1(
    original_pseudobulk
  )

original_log2plus1_qn <-
  create_quantile_normalized(
    original_log2plus1
  )

crescendo_log2plus1 <-
  create_log2plus1(
    crescendo_pseudobulk
  )

crescendo_log2plus1_qn <-
  create_quantile_normalized(
    crescendo_log2plus1
  )

# ==============================================================================
# 12. Prepare detection tables
# ==============================================================================

stage_message(
  7,
  "PREPARING DETECTION TABLES"
)

make_detection_table <- function(
    detection_matrix
) {
  annotation_current <- gene_universe[
    match(
      rownames(
        detection_matrix
      ),
      gene_universe$ensembl_gene_id
    ),
    ,
    drop = FALSE
  ]

  dplyr::bind_cols(
    annotation_current |>
      dplyr::select(
        "ensembl_gene_id",
        "gene",
        "chromosome"
      ),
    tibble::tibble(
      mean_percent_positive =
        rowMeans(
          detection_matrix
        ),
      sd_percent_positive =
        calculate_row_sd(
          detection_matrix
        )
    ),
    tibble::as_tibble(
      as.data.frame(
        detection_matrix,
        check.names = FALSE
      ),
      .name_repair = "minimal"
    )
  ) |>
    dplyr::arrange(
      dplyr::desc(
        .data$mean_percent_positive
      ),
      .data$sd_percent_positive
    )
}

original_detection_table <-
  make_detection_table(
    original_detection
  )

crescendo_detection_table <-
  make_detection_table(
    crescendo_detection
  )

# ==============================================================================
# 13. Write ORIGINAL inputs
# ==============================================================================

stage_message(
  8,
  "WRITING ORIGINAL INPUTS"
)

original_prefix <- paste0(
  dataset_name,
  "_original_pseudobulkPerTissue"
)

write_tsv_gz(
  matrix_to_output_table(
    original_pseudobulk,
    gene_universe
  ),
  file.path(
    original_dirs$pseudobulk,
    paste0(
      original_prefix,
      "_counts.tsv.gz"
    )
  )
)

write_tsv_gz(
  matrix_to_output_table(
    original_log2plus1,
    gene_universe
  ),
  file.path(
    original_dirs$log2plus1,
    paste0(
      original_prefix,
      "_log2plus1.tsv.gz"
    )
  )
)

write_tsv_gz(
  matrix_to_output_table(
    original_log2plus1_qn,
    gene_universe
  ),
  file.path(
    original_dirs$log2plus1_qn,
    paste0(
      original_prefix,
      "_log2plus1_QN.tsv.gz"
    )
  )
)

write_tsv_gz(
  original_detection_table,
  file.path(
    original_dirs$detection,
    paste0(
      original_prefix,
      "_spotDetectionPercent_meanSD.tsv.gz"
    )
  )
)

# ==============================================================================
# 14. Write CRESCENDO inputs
# ==============================================================================

stage_message(
  9,
  "WRITING CRESCENDO INPUTS"
)

crescendo_prefix <- paste0(
  dataset_name,
  "_crescendo_pseudobulkPerTissue"
)

write_tsv_gz(
  matrix_to_output_table(
    crescendo_pseudobulk,
    gene_universe
  ),
  file.path(
    crescendo_dirs$pseudobulk,
    paste0(
      crescendo_prefix,
      "_counts.tsv.gz"
    )
  )
)

write_tsv_gz(
  matrix_to_output_table(
    crescendo_log2plus1,
    gene_universe
  ),
  file.path(
    crescendo_dirs$log2plus1,
    paste0(
      crescendo_prefix,
      "_log2plus1.tsv.gz"
    )
  )
)

write_tsv_gz(
  matrix_to_output_table(
    crescendo_log2plus1_qn,
    gene_universe
  ),
  file.path(
    crescendo_dirs$log2plus1_qn,
    paste0(
      crescendo_prefix,
      "_log2plus1_QN.tsv.gz"
    )
  )
)

write_tsv_gz(
  crescendo_detection_table,
  file.path(
    crescendo_dirs$detection,
    paste0(
      crescendo_prefix,
      "_spotDetectionPercent_meanSD.tsv.gz"
    )
  )
)

# ==============================================================================
# 15. edgeR ORIGINAL
# ==============================================================================

stage_message(
  10,
  "RUNNING edgeR — ORIGINAL"
)

original_edger_xlsx <- file.path(
  original_dirs$edger,
  paste0(
    dataset_name,
    "_original_pseudobulkPerTissue_",
    "edgeR_allGenes_nineTests.xlsx"
  )
)

original_edger <- run_edger(
  count_matrix =
    original_pseudobulk,
  sample_metadata =
    sample_metadata,
  gene_universe =
    gene_universe,
  output_file =
    original_edger_xlsx,
  variant_label =
    "ORIGINAL"
)

# ==============================================================================
# 16. edgeR CRESCENDO
# ==============================================================================

stage_message(
  11,
  "RUNNING edgeR — CRESCENDO"
)

crescendo_edger_xlsx <- file.path(
  crescendo_dirs$edger,
  paste0(
    dataset_name,
    "_crescendo_pseudobulkPerTissue_",
    "edgeR_allGenes_nineTests.xlsx"
  )
)

crescendo_edger <- run_edger(
  count_matrix =
    crescendo_pseudobulk,
  sample_metadata =
    sample_metadata,
  gene_universe =
    gene_universe,
  output_file =
    crescendo_edger_xlsx,
  variant_label =
    "CRESCENDO"
)

# ==============================================================================
# 17. Final validation and report
# ==============================================================================

stage_message(
  12,
  "FINAL VALIDATION"
)

if (!identical(
  rownames(
    original_pseudobulk
  ),
  rownames(
    crescendo_pseudobulk
  )
)) {
  stop(
    "ORIGINAL and CRESCENDO gene order differs.",
    call. = FALSE
  )
}

if (!identical(
  colnames(
    original_pseudobulk
  ),
  colnames(
    crescendo_pseudobulk
  )
)) {
  stop(
    "ORIGINAL and CRESCENDO sample order differs.",
    call. = FALSE
  )
}

if (!file.exists(
  original_edger_xlsx
)) {
  stop(
    "ORIGINAL edgeR XLSX was not created.",
    call. = FALSE
  )
}

if (!file.exists(
  crescendo_edger_xlsx
)) {
  stop(
    "CRESCENDO edgeR XLSX was not created.",
    call. = FALSE
  )
}

stage_message(
  13,
  "FINISHED"
)

message(
  "Samples: ",
  ncol(
    original_pseudobulk
  )
)

message(
  "Common protein-coding chr1-19/X genes: ",
  format(
    nrow(
      original_pseudobulk
    ),
    big.mark = ","
  )
)

message("")
message(
  "ORIGINAL pseudobulk dimensions: ",
  nrow(
    original_pseudobulk
  ),
  " genes x ",
  ncol(
    original_pseudobulk
  ),
  " samples"
)

message(
  "CRESCENDO pseudobulk dimensions: ",
  nrow(
    crescendo_pseudobulk
  ),
  " genes x ",
  ncol(
    crescendo_pseudobulk
  ),
  " samples"
)

message("")
message("ORIGINAL output:")
message(original_root)

message("")
message("CRESCENDO output:")
message(crescendo_comparison_root)

message("")
message("ORIGINAL edgeR:")
message(original_edger_xlsx)

message("")
message("CRESCENDO edgeR:")
message(crescendo_edger_xlsx)

message("")
message("Comparison root:")
message(comparison_root)
