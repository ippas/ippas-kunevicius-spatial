#!/usr/bin/env Rscript
# ==============================================================================
# 02b_pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction_without13_1F.R
#
# STANDALONE SCRIPT
#
# Full per-cluster edgeR sensitivity analysis using Crescendo-corrected pseudobulk counts.
#
# SAMPLE SENSITIVITY ANALYSIS:
#   sample 13_1F is excluded from ALL edgeR statistical calculations.
#   Crescendo correction itself is NOT recalculated here; this script uses the
#   Outputs are written under:
#     02_edgeR/02b_without13_1F/
#
#   already generated Crescendo-corrected counts and removes 13_1F before edgeR.
# No external functions/src file is sourced.
#
# Count input:
#   crescendo_corrected_counts_by_cluster
#
# Model:
#   expression ~ fmt_donor_group * sex
#
# Tests per cluster:
#   1. Interaction
#   2. Overall sex effect
#   3. Overall donor-group effect
#   4-9. Six pairwise comparisons between the four Sex x donor-group cells
#
# Statistical workflow:
#   Crescendo-corrected pseudobulk counts
#     -> protein-coding chromosomes 1-19 + X
#     -> edgeR::filterByExpr()
#     -> TMM
#     -> estimateDisp(robust=TRUE)
#     -> glmQLFit(robust=TRUE)
#     -> glmQLFTest()
#     -> BH FDR
#
# Descriptive outputs:
#   - raw-spot percent-positive metrics from the original Seurat raw counts;
#   - TMM-normalized CPM;
#   - log2(Crescendo corrected pseudobulk count + 1) + quantile normalization.
#
# IMPORTANT:
#   - raw Seurat counts are NEVER used as edgeR count input in this script;
#   - log2+QN values are NEVER used as edgeR count input;
#   - the edgeR model/contrasts match the previous raw-pseudobulk workflow.
# ==============================================================================

options(stringsAsFactors = FALSE)

SCRIPT_BUILD <- "20260831_02b_WITHOUT13_1F_V1"

# ==============================================================================
# 1. PROJECT CONFIGURATION
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"
clustering_name <- "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

cluster_column <- "leiden_res0.40"
sample_column <- "sample_ID"
assay_name <- "RNA"

excluded_sample_ids <- c("13_1F")
expected_number_of_samples <- 15L
expected_number_of_clusters <- 16L
min_spots_per_sample_cluster <- 0L
min_samples_per_factor_cell <- 2L
positive_spot_min_count <- 1L
xlsx_max_decimal_places <- 5L

filter_by_expr_parameters <- list(
  min.count = 10,
  min.total.count = 15,
  large.n = 10,
  min.prop = 0.7
)

# ==============================================================================
# 2. INPUT FILES
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

annotation_metadata_file <- file.path(
  project_root,
  "data",
  "ensembl115",
  "ensembl115_mouse_proteinCodingGenes_metadata.tsv"
)

crescendo_statistics_root <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)

crescendo_input_rds <- file.path(
  crescendo_statistics_root,
  "01_crescendoBatchCorrection",
  "07_objects",
  paste0(
    dataset_name,
    "_",
    clustering_name,
    "_crescendoBatchCorrection_pseudobulkCounts_beforeEdgeR.rds"
  )
)

# ==============================================================================
# 3. OUTPUT DIRECTORY
# ==============================================================================

output_statistics_dir <- file.path(
  crescendo_statistics_root,
  "02_edgeR",
  "02b_without13_1F"
)

analysis_prefix <- paste0(
  dataset_name,
  "_",
  clustering_name,
  "_crescendoCorrectedPseudobulk_edgeR_twoFactorInteraction_without13_1F"
)

# ==============================================================================
# 4. BASIC HELPERS
# ==============================================================================

stage_message <- function(stage_number, stage_name) {
  message("\n============================================================")
  message("STAGE ", stage_number, ": ", stage_name)
  message("============================================================")
}

write_tsv <- function(object, output_file) {
  output_connection <- if (grepl("\\.gz$", output_file, ignore.case = TRUE)) {
    gzfile(output_file, open = "wt", encoding = "UTF-8")
  } else {
    file(output_file, open = "wt", encoding = "UTF-8")
  }
  on.exit(close(output_connection), add = TRUE)
  utils::write.table(
    object,
    file = output_connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  invisible(output_file)
}

sanitize_file_component <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(x == "", "unnamed", x)
}

natural_cluster_order <- function(cluster_ids) {
  cluster_ids <- unique(as.character(cluster_ids))
  numeric_ids <- suppressWarnings(as.numeric(cluster_ids))
  if (all(!is.na(numeric_ids))) {
    return(cluster_ids[order(numeric_ids)])
  }
  cluster_ids[order(cluster_ids)]
}

safe_row_sd <- function(x) {
  if (ncol(x) <= 1L) {
    return(rep(NA_real_, nrow(x)))
  }
  apply(x, 1L, stats::sd)
}

# ==============================================================================
# 5. PACKAGE CHECK
# ==============================================================================

check_required_packages <- function() {
  required_packages <- c(
    "Seurat",
    "SeuratObject",
    "Matrix",
    "edgeR",
    "limma",
    "dplyr",
    "tibble",
    "openxlsx",
    "preprocessCore"
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

  invisible(TRUE)
}

# ==============================================================================
# 6. SEURAT HELPERS
# ==============================================================================

load_single_seurat_object_from_rdata <- function(
    input_rdata_file,
    requested_object_name = NULL
) {
  load_environment <- new.env(parent = emptyenv())
  loaded_object_names <- load(input_rdata_file, envir = load_environment)

  if (length(loaded_object_names) == 0L) {
    stop("No objects loaded from RData.", call. = FALSE)
  }

  if (!is.null(requested_object_name)) {
    if (!requested_object_name %in% loaded_object_names) {
      stop(
        "Requested object not found: ",
        requested_object_name,
        call. = FALSE
      )
    }
    requested_object <- get(requested_object_name, envir = load_environment)
    if (!inherits(requested_object, "Seurat")) {
      stop("Requested object is not a Seurat object.", call. = FALSE)
    }
    return(
      list(
        object = requested_object,
        object_name = requested_object_name
      )
    )
  }

  seurat_object_names <- loaded_object_names[
    vapply(
      loaded_object_names,
      function(object_name) {
        inherits(get(object_name, envir = load_environment), "Seurat")
      },
      FUN.VALUE = logical(1)
    )
  ]

  if (length(seurat_object_names) != 1L) {
    stop(
      "Expected exactly one Seurat object; found: ",
      paste(seurat_object_names, collapse = ", "),
      call. = FALSE
    )
  }

  selected_name <- seurat_object_names[[1L]]

  list(
    object = get(selected_name, envir = load_environment),
    object_name = selected_name
  )
}

resolve_cluster_column <- function(
    seurat_metadata,
    requested_cluster_column,
    expected_number_of_clusters
) {
  if (!requested_cluster_column %in% colnames(seurat_metadata)) {
    stop(
      "Cluster column does not exist: ",
      requested_cluster_column,
      call. = FALSE
    )
  }

  n_clusters <- length(
    unique(
      stats::na.omit(
        seurat_metadata[[requested_cluster_column]]
      )
    )
  )

  if (n_clusters != expected_number_of_clusters) {
    stop(
      "Cluster column ",
      requested_cluster_column,
      " contains ",
      n_clusters,
      " clusters; expected ",
      expected_number_of_clusters,
      ".",
      call. = FALSE
    )
  }

  requested_cluster_column
}

extract_raw_counts_from_seurat_layers <- function(
    seurat_object,
    assay_name
) {
  if (!assay_name %in% names(seurat_object@assays)) {
    stop(
      "Assay not present: ",
      assay_name,
      call. = FALSE
    )
  }

  assay_object <- seurat_object[[assay_name]]

  count_layer_names <- tryCatch(
    SeuratObject::Layers(
      object = assay_object,
      search = "^counts"
    ),
    error = function(e) character(0)
  )

  if (length(count_layer_names) == 0L) {
    count_matrix <- Seurat::GetAssayData(
      object = seurat_object,
      assay = assay_name,
      slot = "counts"
    )

    if (!inherits(count_matrix, "Matrix")) {
      count_matrix <- Matrix::Matrix(
        count_matrix,
        sparse = TRUE
      )
    }

    return(
      list(
        counts = count_matrix,
        count_layers = "counts"
      )
    )
  }

  count_layer_matrices <- lapply(
    count_layer_names,
    function(layer_name) {
      layer_matrix <- SeuratObject::LayerData(
        object = assay_object,
        layer = layer_name
      )

      if (!inherits(layer_matrix, "Matrix")) {
        layer_matrix <- Matrix::Matrix(
          layer_matrix,
          sparse = TRUE
        )
      }

      layer_matrix
    }
  )

  reference_features <- rownames(
    count_layer_matrices[[1L]]
  )

  for (layer_index in seq_along(count_layer_matrices)) {
    current_features <- rownames(
      count_layer_matrices[[layer_index]]
    )

    if (!setequal(current_features, reference_features)) {
      stop(
        "Count layers do not contain identical feature sets. Problematic layer: ",
        count_layer_names[[layer_index]],
        call. = FALSE
      )
    }

    count_layer_matrices[[layer_index]] <-
      count_layer_matrices[[layer_index]][
        reference_features,
        ,
        drop = FALSE
      ]
  }

  count_matrix <- Reduce(
    cbind,
    count_layer_matrices
  )

  if (anyDuplicated(colnames(count_matrix))) {
    stop(
      "Duplicated spot names after combining count layers.",
      call. = FALSE
    )
  }

  list(
    counts = count_matrix,
    count_layers = count_layer_names
  )
}

# ==============================================================================
# 7. ANNOTATION HELPERS
# ==============================================================================

read_and_validate_annotation <- function(
    annotation_rds_file,
    annotation_metadata_file,
    ensembl_release = 115L,
    genome_assembly = "GRCm39",
    analysis_chromosomes = c(as.character(1:19), "X")
) {
  annotation_metadata <- utils::read.delim(
    annotation_metadata_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_metadata_columns <- c(
    "genome_assembly",
    "ensembl_release",
    "retained_feature_type",
    "retained_gene_biotype",
    "number_of_protein_coding_genes"
  )

  missing_metadata_columns <- setdiff(
    required_metadata_columns,
    colnames(annotation_metadata)
  )

  if (length(missing_metadata_columns) > 0L) {
    stop(
      "Missing annotation metadata columns: ",
      paste(missing_metadata_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(annotation_metadata) != 1L) {
    stop(
      "Annotation metadata must contain exactly one row.",
      call. = FALSE
    )
  }

  if (
    as.integer(annotation_metadata$ensembl_release[[1L]]) !=
      ensembl_release
  ) {
    stop(
      "Unexpected Ensembl release.",
      call. = FALSE
    )
  }

  if (
    as.character(annotation_metadata$genome_assembly[[1L]]) !=
      genome_assembly
  ) {
    stop(
      "Unexpected genome assembly.",
      call. = FALSE
    )
  }

  annotation_all <- readRDS(
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
    colnames(annotation_all)
  )

  if (length(missing_annotation_columns) > 0L) {
    stop(
      "Missing annotation RDS columns: ",
      paste(missing_annotation_columns, collapse = ", "),
      call. = FALSE
    )
  }

  annotation_all <- annotation_all |>
    dplyr::transmute(
      chromosome = trimws(as.character(.data$chromosome)),
      ensembl_gene_id = sub(
        "\\.[0-9]+$",
        "",
        trimws(as.character(.data$ensembl_gene_id))
      ),
      gene = trimws(as.character(.data$gene)),
      gene_biotype = trimws(as.character(.data$gene_biotype))
    ) |>
    dplyr::filter(
      !is.na(.data$chromosome),
      .data$chromosome != "",
      !is.na(.data$ensembl_gene_id),
      .data$ensembl_gene_id != "",
      .data$gene_biotype == "protein_coding"
    ) |>
    dplyr::mutate(
      gene = dplyr::if_else(
        is.na(.data$gene) | .data$gene == "",
        .data$ensembl_gene_id,
        .data$gene
      )
    ) |>
    dplyr::distinct(
      .data$ensembl_gene_id,
      .keep_all = TRUE
    )

  annotation_analysis <- annotation_all |>
    dplyr::filter(
      .data$chromosome %in% analysis_chromosomes
    ) |>
    dplyr::mutate(
      chromosome = factor(
        .data$chromosome,
        levels = analysis_chromosomes
      )
    ) |>
    dplyr::arrange(
      .data$chromosome,
      .data$ensembl_gene_id
    ) |>
    dplyr::mutate(
      chromosome = as.character(.data$chromosome)
    )

  list(
    annotation_all = annotation_all,
    annotation_analysis = annotation_analysis,
    annotation_metadata = annotation_metadata
  )
}

map_features_to_annotation <- function(
    feature_names,
    annotation_all,
    annotation_analysis
) {
  feature_names <- as.character(feature_names)

  if (anyDuplicated(feature_names)) {
    stop(
      "Feature names are duplicated.",
      call. = FALSE
    )
  }

  normalized_features <- sub(
    "\\.[0-9]+$",
    "",
    feature_names
  )

  ensembl_fraction <- mean(
    grepl(
      "^ENSMUSG[0-9]+$",
      normalized_features
    )
  )

  if (ensembl_fraction >= 0.5) {
    mapping <- tibble::tibble(
      feature_id = feature_names,
      normalized_feature_id = normalized_features
    ) |>
      dplyr::left_join(
        annotation_all,
        by = c(
          "normalized_feature_id" = "ensembl_gene_id"
        )
      ) |>
      dplyr::rename(
        ensembl_gene_id = .data$normalized_feature_id
      )
    identifier_type <- "ensembl_gene_id"
  } else {
    unique_symbol_annotation <- annotation_all |>
      dplyr::add_count(
        .data$gene,
        name = "symbol_annotation_count"
      ) |>
      dplyr::filter(
        .data$symbol_annotation_count == 1L
      ) |>
      dplyr::select(
        -.data$symbol_annotation_count
      )

    mapping <- tibble::tibble(
      feature_id = feature_names,
      gene_lookup = feature_names
    ) |>
      dplyr::left_join(
        unique_symbol_annotation,
        by = c("gene_lookup" = "gene")
      ) |>
      dplyr::mutate(
        gene = .data$gene_lookup
      ) |>
      dplyr::select(
        .data$feature_id,
        .data$chromosome,
        .data$ensembl_gene_id,
        .data$gene,
        .data$gene_biotype
      )
    identifier_type <- "gene_symbol"
  }

  analysis_gene_ids <- annotation_analysis$ensembl_gene_id

  mapping <- mapping |>
    dplyr::mutate(
      found_in_protein_coding_annotation =
        !is.na(.data$ensembl_gene_id),
      located_on_chr1to19_or_X =
        !is.na(.data$ensembl_gene_id) &
          .data$ensembl_gene_id %in% analysis_gene_ids,
      keep_for_edger =
        .data$found_in_protein_coding_annotation &
          .data$located_on_chr1to19_or_X
    )

  retained_mapping <- mapping |>
    dplyr::filter(
      .data$keep_for_edger
    ) |>
    dplyr::select(
      .data$feature_id,
      .data$ensembl_gene_id
    ) |>
    dplyr::distinct(
      .data$ensembl_gene_id,
      .keep_all = TRUE
    )

  if (nrow(retained_mapping) == 0L) {
    stop(
      "No features matched protein-coding genes on chromosomes 1-19/X.",
      call. = FALSE
    )
  }

  retained_annotation <- annotation_analysis[
    match(
      retained_mapping$ensembl_gene_id,
      annotation_analysis$ensembl_gene_id
    ),
    ,
    drop = FALSE
  ] |>
    dplyr::select(
      .data$chromosome,
      .data$ensembl_gene_id,
      .data$gene,
      .data$gene_biotype
    )

  if (anyNA(retained_annotation$ensembl_gene_id)) {
    stop(
      "Annotation matching produced NA Ensembl IDs.",
      call. = FALSE
    )
  }

  list(
    mapping_status = mapping,
    retained_mapping = retained_mapping,
    retained_annotation = retained_annotation,
    identifier_type = identifier_type
  )
}

# ==============================================================================
# 8. DESIGN AND CONTRASTS
# ==============================================================================

create_two_factor_design_and_contrasts <- function(
    cluster_sample_metadata
) {
  model_metadata <- cluster_sample_metadata |>
    dplyr::transmute(
      sample_ID = as.character(.data$sample_ID),
      fmt_donor_group = factor(
        as.character(.data$fmt_donor_group),
        levels = c("Neurotypical", "ASD")
      ),
      sex = factor(
        as.character(.data$sex),
        levels = c("Male", "Female")
      )
    )

  if (
    anyNA(model_metadata$fmt_donor_group) ||
      anyNA(model_metadata$sex)
  ) {
    stop(
      "Missing or unexpected sex/group values.",
      call. = FALSE
    )
  }

  design <- stats::model.matrix(
    ~ fmt_donor_group * sex,
    data = model_metadata
  )

  rownames(design) <- model_metadata$sample_ID

  expected_original_columns <- c(
    "(Intercept)",
    "fmt_donor_groupASD",
    "sexFemale",
    "fmt_donor_groupASD:sexFemale"
  )

  if (!identical(
    colnames(design),
    expected_original_columns
  )) {
    stop(
      "Unexpected design columns: ",
      paste(colnames(design), collapse = ", "),
      call. = FALSE
    )
  }

  colnames(design) <- c(
    "Intercept",
    "groupASD",
    "sexFemale",
    "groupASD_sexFemale"
  )

  if (qr(design)$rank != ncol(design)) {
    stop(
      "Design matrix is not full rank.",
      call. = FALSE
    )
  }

  contrasts <- limma::makeContrasts(
    Interaction = groupASD_sexFemale,

    Overall_Sex_Female_vs_Male =
      sexFemale + 0.5 * groupASD_sexFemale,

    Overall_Group_ASD_vs_Neurotypical =
      groupASD + 0.5 * groupASD_sexFemale,

    ASD_Male_vs_Neurotypical_Male =
      groupASD,

    ASD_Female_vs_Neurotypical_Female =
      groupASD + groupASD_sexFemale,

    Neurotypical_Female_vs_Neurotypical_Male =
      sexFemale,

    ASD_Female_vs_ASD_Male =
      sexFemale + groupASD_sexFemale,

    ASD_Female_vs_Neurotypical_Male =
      groupASD + sexFemale + groupASD_sexFemale,

    ASD_Male_vs_Neurotypical_Female =
      groupASD - sexFemale,

    levels = design
  )

  test_definitions <- tibble::tribble(
    ~test_id, ~sheet_name, ~comparison, ~add_regulation,

    "Interaction",
    "Interaction",
    "(ASD - Neurotypical)_Female - (ASD - Neurotypical)_Male",
    FALSE,

    "Overall_Sex_Female_vs_Male",
    "Overall_Sex",
    "Female vs Male averaged equally across donor groups",
    TRUE,

    "Overall_Group_ASD_vs_Neurotypical",
    "Overall_Group",
    "ASD vs Neurotypical averaged equally across sexes",
    TRUE,

    "ASD_Male_vs_Neurotypical_Male",
    "Group_in_Male",
    "ASD Male vs Neurotypical Male",
    TRUE,

    "ASD_Female_vs_Neurotypical_Female",
    "Group_in_Female",
    "ASD Female vs Neurotypical Female",
    TRUE,

    "Neurotypical_Female_vs_Neurotypical_Male",
    "Sex_in_Neurotypical",
    "Neurotypical Female vs Neurotypical Male",
    TRUE,

    "ASD_Female_vs_ASD_Male",
    "Sex_in_ASD",
    "ASD Female vs ASD Male",
    TRUE,

    "ASD_Female_vs_Neurotypical_Male",
    "ASD_F_vs_NT_M",
    "ASD Female vs Neurotypical Male",
    TRUE,

    "ASD_Male_vs_Neurotypical_Female",
    "ASD_M_vs_NT_F",
    "ASD Male vs Neurotypical Female",
    TRUE
  )

  if (!identical(
    colnames(contrasts),
    test_definitions$test_id
  )) {
    stop(
      "Contrast order does not match test definitions.",
      call. = FALSE
    )
  }

  list(
    metadata = model_metadata,
    design = design,
    contrasts = contrasts,
    test_definitions = test_definitions
  )
}

# ==============================================================================
# 9. DESCRIPTIVE METRICS
# ==============================================================================

compute_detection_percentages <- function(
    raw_spot_count_matrix,
    raw_spot_metadata,
    included_metadata,
    cluster_id,
    tested_gene_ids,
    positive_spot_min_count
) {
  sample_ids <- as.character(
    included_metadata$sample_ID
  )

  detection_matrix <- matrix(
    NA_real_,
    nrow = length(tested_gene_ids),
    ncol = length(sample_ids),
    dimnames = list(
      tested_gene_ids,
      sample_ids
    )
  )

  selected_cluster_spots <-
    as.character(
      raw_spot_metadata$cluster_id_analysis
    ) == as.character(cluster_id)

  for (sample_id in sample_ids) {
    sample_spot_ids <- rownames(
      raw_spot_metadata
    )[
      selected_cluster_spots &
        as.character(
          raw_spot_metadata$sample_ID_analysis
        ) == sample_id
    ]

    if (length(sample_spot_ids) == 0L) {
      stop(
        "No raw spots for sample ",
        sample_id,
        " in cluster ",
        cluster_id,
        ".",
        call. = FALSE
      )
    }

    sample_gene_counts <- raw_spot_count_matrix[
      tested_gene_ids,
      sample_spot_ids,
      drop = FALSE
    ]

    detection_matrix[, sample_id] <-
      100 *
      as.numeric(
        Matrix::rowSums(
          sample_gene_counts >=
            positive_spot_min_count
        )
      ) /
      length(sample_spot_ids)
  }

  detection_matrix
}

make_sample_level_descriptive_tables <- function(
    tmm_cpm_matrix,
    detection_percentage_matrix,
    model_metadata,
    test_definitions
) {
  model_metadata <- model_metadata |>
    dplyr::mutate(
      sample_ID = as.character(.data$sample_ID),
      fmt_donor_group =
        as.character(.data$fmt_donor_group),
      sex = as.character(.data$sex)
    )

  tested_gene_ids <- rownames(
    tmm_cpm_matrix
  )

  summarize_samples <- function(
      sample_ids,
      displayed_name
  ) {
    sample_ids <- as.character(sample_ids)

    detection_subset <-
      detection_percentage_matrix[
        ,
        sample_ids,
        drop = FALSE
      ]

    cpm_subset <- tmm_cpm_matrix[
      ,
      sample_ids,
      drop = FALSE
    ]

    tibble::tibble(
      ensembl_gene_id = tested_gene_ids,

      !!paste0(
        "mean_percent_positive_spots_",
        displayed_name
      ) :=
        rowMeans(detection_subset),

      !!paste0(
        "sd_percent_positive_spots_",
        displayed_name
      ) :=
        safe_row_sd(detection_subset),

      !!paste0(
        "mean_TMM_CPM_",
        displayed_name
      ) :=
        rowMeans(cpm_subset),

      !!paste0(
        "sd_TMM_CPM_",
        displayed_name
      ) :=
        safe_row_sd(cpm_subset)
    )
  }

  select_samples <- function(
      sex_value = NULL,
      group_value = NULL
  ) {
    keep <- rep(
      TRUE,
      nrow(model_metadata)
    )

    if (!is.null(sex_value)) {
      keep <- keep &
        model_metadata$sex == sex_value
    }

    if (!is.null(group_value)) {
      keep <- keep &
        model_metadata$fmt_donor_group ==
          group_value
    }

    model_metadata$sample_ID[keep]
  }

  join_tables <- function(...) {
    Reduce(
      function(x, y) {
        dplyr::left_join(
          x,
          y,
          by = "ensembl_gene_id"
        )
      },
      list(...)
    )
  }

  out <- list(
    Interaction = join_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Male_Neurotypical"
      ),
      summarize_samples(
        select_samples("Male", "ASD"),
        "Male_ASD"
      ),
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Female_Neurotypical"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "Female_ASD"
      )
    ),

    Overall_Sex_Female_vs_Male = join_tables(
      summarize_samples(
        select_samples(sex_value = "Male"),
        "Male"
      ),
      summarize_samples(
        select_samples(sex_value = "Female"),
        "Female"
      )
    ),

    Overall_Group_ASD_vs_Neurotypical = join_tables(
      summarize_samples(
        select_samples(
          group_value = "Neurotypical"
        ),
        "Neurotypical"
      ),
      summarize_samples(
        select_samples(
          group_value = "ASD"
        ),
        "ASD"
      )
    ),

    ASD_Male_vs_Neurotypical_Male = join_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Neurotypical_Male"
      ),
      summarize_samples(
        select_samples("Male", "ASD"),
        "ASD_Male"
      )
    ),

    ASD_Female_vs_Neurotypical_Female = join_tables(
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Neurotypical_Female"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "ASD_Female"
      )
    ),

    Neurotypical_Female_vs_Neurotypical_Male = join_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Neurotypical_Male"
      ),
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Neurotypical_Female"
      )
    ),

    ASD_Female_vs_ASD_Male = join_tables(
      summarize_samples(
        select_samples("Male", "ASD"),
        "ASD_Male"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "ASD_Female"
      )
    ),

    ASD_Female_vs_Neurotypical_Male = join_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Neurotypical_Male"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "ASD_Female"
      )
    ),

    ASD_Male_vs_Neurotypical_Female = join_tables(
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Neurotypical_Female"
      ),
      summarize_samples(
        select_samples("Male", "ASD"),
        "ASD_Male"
      )
    )
  )

  expected_ids <- test_definitions$test_id

  if (!setequal(names(out), expected_ids)) {
    stop(
      "Descriptive tables do not match tests.",
      call. = FALSE
    )
  }

  out[expected_ids]
}

create_log2qn_matrix <- function(
    corrected_count_matrix
) {
  dense_counts <- as.matrix(
    corrected_count_matrix
  )

  storage.mode(dense_counts) <- "double"

  log2plus1 <- log2(
    dense_counts + 1
  )

  qn <- preprocessCore::normalize.quantiles(
    log2plus1,
    copy = TRUE
  )

  rownames(qn) <- rownames(
    dense_counts
  )

  colnames(qn) <- colnames(
    dense_counts
  )

  qn
}

make_log2qn_descriptive_tables <- function(
    log2_qn_matrix,
    model_metadata,
    test_definitions
) {
  model_metadata <- model_metadata |>
    dplyr::mutate(
      sample_ID = as.character(.data$sample_ID),
      fmt_donor_group =
        as.character(.data$fmt_donor_group),
      sex = as.character(.data$sex)
    )

  summarize_samples <- function(
      sample_ids,
      displayed_name
  ) {
    values <- log2_qn_matrix[
      ,
      sample_ids,
      drop = FALSE
    ]

    tibble::tibble(
      ensembl_gene_id = rownames(
        log2_qn_matrix
      ),

      !!paste0(
        "mean_log2plus1_QN_",
        displayed_name
      ) :=
        rowMeans(values),

      !!paste0(
        "sd_log2plus1_QN_",
        displayed_name
      ) :=
        safe_row_sd(values)
    )
  }

  select_samples <- function(
      sex_value = NULL,
      group_value = NULL
  ) {
    keep <- rep(
      TRUE,
      nrow(model_metadata)
    )

    if (!is.null(sex_value)) {
      keep <- keep &
        model_metadata$sex == sex_value
    }

    if (!is.null(group_value)) {
      keep <- keep &
        model_metadata$fmt_donor_group ==
          group_value
    }

    model_metadata$sample_ID[keep]
  }

  join_tables <- function(...) {
    Reduce(
      function(x, y) {
        dplyr::left_join(
          x,
          y,
          by = "ensembl_gene_id"
        )
      },
      list(...)
    )
  }

  out <- list(
    Interaction = join_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Male_Neurotypical"
      ),
      summarize_samples(
        select_samples("Male", "ASD"),
        "Male_ASD"
      ),
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Female_Neurotypical"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "Female_ASD"
      )
    ),

    Overall_Sex_Female_vs_Male = join_tables(
      summarize_samples(
        select_samples(sex_value = "Male"),
        "Male"
      ),
      summarize_samples(
        select_samples(sex_value = "Female"),
        "Female"
      )
    ),

    Overall_Group_ASD_vs_Neurotypical = join_tables(
      summarize_samples(
        select_samples(
          group_value = "Neurotypical"
        ),
        "Neurotypical"
      ),
      summarize_samples(
        select_samples(
          group_value = "ASD"
        ),
        "ASD"
      )
    ),

    ASD_Male_vs_Neurotypical_Male = join_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Neurotypical_Male"
      ),
      summarize_samples(
        select_samples("Male", "ASD"),
        "ASD_Male"
      )
    ),

    ASD_Female_vs_Neurotypical_Female = join_tables(
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Neurotypical_Female"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "ASD_Female"
      )
    ),

    Neurotypical_Female_vs_Neurotypical_Male = join_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Neurotypical_Male"
      ),
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Neurotypical_Female"
      )
    ),

    ASD_Female_vs_ASD_Male = join_tables(
      summarize_samples(
        select_samples("Male", "ASD"),
        "ASD_Male"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "ASD_Female"
      )
    ),

    ASD_Female_vs_Neurotypical_Male = join_tables(
      summarize_samples(
        select_samples("Male", "Neurotypical"),
        "Neurotypical_Male"
      ),
      summarize_samples(
        select_samples("Female", "ASD"),
        "ASD_Female"
      )
    ),

    ASD_Male_vs_Neurotypical_Female = join_tables(
      summarize_samples(
        select_samples("Female", "Neurotypical"),
        "Neurotypical_Female"
      ),
      summarize_samples(
        select_samples("Male", "ASD"),
        "ASD_Male"
      )
    )
  )

  expected_ids <- test_definitions$test_id

  if (!setequal(names(out), expected_ids)) {
    stop(
      "log2+QN descriptive tables do not match tests.",
      call. = FALSE
    )
  }

  out[expected_ids]
}

# ==============================================================================
# 10. RESULT HELPERS
# ==============================================================================

create_result_table <- function(
    qlf_test,
    gene_annotation,
    add_regulation,
    descriptive_table,
    log2qn_table
) {
  result <- edgeR::topTags(
    qlf_test,
    n = Inf,
    sort.by = "PValue"
  )$table |>
    as.data.frame(
      stringsAsFactors = FALSE
    ) |>
    tibble::rownames_to_column(
      var = "ensembl_gene_id"
    ) |>
    dplyr::left_join(
      gene_annotation,
      by = "ensembl_gene_id"
    ) |>
    dplyr::mutate(
      gene = dplyr::if_else(
        is.na(.data$gene) |
          .data$gene == "",
        .data$ensembl_gene_id,
        .data$gene
      )
    )

  if (add_regulation) {
    result <- result |>
      dplyr::mutate(
        regulation = dplyr::if_else(
          .data$logFC >= 0,
          "up",
          "down"
        )
      ) |>
      dplyr::select(
        .data$ensembl_gene_id,
        .data$gene,
        .data$chromosome,
        .data$logFC,
        .data$logCPM,
        .data$F,
        .data$PValue,
        .data$FDR,
        .data$regulation
      )
  } else {
    result <- result |>
      dplyr::select(
        .data$ensembl_gene_id,
        .data$gene,
        .data$chromosome,
        .data$logFC,
        .data$logCPM,
        .data$F,
        .data$PValue,
        .data$FDR
      )
  }

  result |>
    dplyr::left_join(
      descriptive_table,
      by = "ensembl_gene_id"
    ) |>
    dplyr::left_join(
      log2qn_table,
      by = "ensembl_gene_id"
    ) |>
    dplyr::arrange(
      .data$PValue,
      dplyr::desc(abs(.data$logFC))
    )
}

write_result_workbook <- function(
    result_list,
    test_definitions,
    output_file,
    include_cluster_id = FALSE,
    xlsx_max_decimal_places = 5L
) {
  workbook <- openxlsx::createWorkbook()

  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#D9EAF7",
    border = "Bottom"
  )

  decimal_pattern <- paste0(
    "#,##0.",
    paste(
      rep("#", xlsx_max_decimal_places),
      collapse = ""
    )
  )

  scientific_pattern <- paste0(
    "0.",
    paste(
      rep("#", xlsx_max_decimal_places),
      collapse = ""
    ),
    "E+00"
  )

  decimal_style <- openxlsx::createStyle(
    numFmt = decimal_pattern
  )

  scientific_style <- openxlsx::createStyle(
    numFmt = scientific_pattern
  )

  for (test_index in seq_len(nrow(test_definitions))) {
    test_id <- test_definitions$test_id[[test_index]]
    sheet_name <- test_definitions$sheet_name[[test_index]]
    result_table <- result_list[[test_id]]

    openxlsx::addWorksheet(
      workbook,
      sheetName = sheet_name
    )

    openxlsx::writeData(
      workbook,
      sheet = sheet_name,
      x = result_table,
      startRow = 1,
      startCol = 1,
      headerStyle = header_style,
      withFilter = TRUE
    )

    openxlsx::freezePane(
      workbook,
      sheet = sheet_name,
      firstActiveRow = 2,
      firstActiveCol = if (include_cluster_id) 4 else 3
    )

    numeric_columns <- which(
      vapply(
        result_table,
        is.numeric,
        FUN.VALUE = logical(1)
      )
    )

    scientific_columns <- which(
      colnames(result_table) %in%
        c("PValue", "FDR")
    )

    regular_numeric_columns <- setdiff(
      numeric_columns,
      scientific_columns
    )

    if (
      length(regular_numeric_columns) > 0L &&
        nrow(result_table) > 0L
    ) {
      openxlsx::addStyle(
        workbook,
        sheet = sheet_name,
        style = decimal_style,
        rows = 2:(nrow(result_table) + 1L),
        cols = regular_numeric_columns,
        gridExpand = TRUE,
        stack = TRUE
      )
    }

    if (
      length(scientific_columns) > 0L &&
        nrow(result_table) > 0L
    ) {
      openxlsx::addStyle(
        workbook,
        sheet = sheet_name,
        style = scientific_style,
        rows = 2:(nrow(result_table) + 1L),
        cols = scientific_columns,
        gridExpand = TRUE,
        stack = TRUE
      )
    }

    openxlsx::setColWidths(
      workbook,
      sheet = sheet_name,
      cols = seq_len(ncol(result_table)),
      widths = "auto"
    )
  }

  openxlsx::saveWorkbook(
    workbook,
    file = output_file,
    overwrite = TRUE
  )

  invisible(output_file)
}

# ==============================================================================
# 11. LOAD AND PREPARE CRESCENDO INPUT
# ==============================================================================

prepare_crescendo_inputs <- function() {
  crescendo_object <- readRDS(
    crescendo_input_rds
  )

  required_objects <- c(
    "sample_metadata",
    "group_table",
    "gene_metadata",
    "crescendo_corrected_counts_by_cluster"
  )

  missing_objects <- setdiff(
    required_objects,
    names(crescendo_object)
  )

  if (length(missing_objects) > 0L) {
    stop(
      "Missing Crescendo RDS object(s): ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  sample_metadata_all <- tibble::as_tibble(
    crescendo_object$sample_metadata
  ) |>
    dplyr::mutate(
      sample_ID = as.character(.data$sample_ID),
      fmt_donor_group = factor(
        as.character(.data$fmt_donor_group),
        levels = c("Neurotypical", "ASD")
      ),
      sex = factor(
        as.character(.data$sex),
        levels = c("Male", "Female")
      )
    )

  missing_excluded_samples <- setdiff(
    excluded_sample_ids,
    sample_metadata_all$sample_ID
  )

  if (length(missing_excluded_samples) > 0L) {
    stop(
      "Requested excluded sample(s) not found in Crescendo metadata: ",
      paste(missing_excluded_samples, collapse = ", "),
      call. = FALSE
    )
  }

  sample_metadata <- sample_metadata_all |>
    dplyr::filter(
      !.data$sample_ID %in% excluded_sample_ids
    ) |>
    droplevels()

  message(
    "Excluded sample(s): ",
    paste(excluded_sample_ids, collapse = ", ")
  )

  message(
    "Samples retained for edgeR: ",
    nrow(sample_metadata)
  )

  print(
    sample_metadata |>
      dplyr::select(
        .data$sample_ID,
        .data$fmt_donor_group,
        .data$sex
      ),
    n = Inf
  )

  if (nrow(sample_metadata) != expected_number_of_samples) {
    stop(
      "Unexpected number of samples in Crescendo RDS: ",
      nrow(sample_metadata),
      call. = FALSE
    )
  }

  corrected_counts_by_cluster <-
    crescendo_object$crescendo_corrected_counts_by_cluster

  cluster_ids <- natural_cluster_order(
    names(corrected_counts_by_cluster)
  )

  if (length(cluster_ids) != expected_number_of_clusters) {
    stop(
      "Unexpected number of corrected clusters: ",
      length(cluster_ids),
      call. = FALSE
    )
  }

  reference_features <- rownames(
    corrected_counts_by_cluster[[cluster_ids[[1L]]]]
  )

  annotation_objects <- read_and_validate_annotation(
    annotation_rds_file = annotation_rds_file,
    annotation_metadata_file = annotation_metadata_file
  )

  feature_mapping <- map_features_to_annotation(
    feature_names = reference_features,
    annotation_all =
      annotation_objects$annotation_all,
    annotation_analysis =
      annotation_objects$annotation_analysis
  )

  retained_mapping <-
    feature_mapping$retained_mapping

  gene_annotation <-
    feature_mapping$retained_annotation

  message(
    "Crescendo corrected genes: ",
    format(
      length(reference_features),
      big.mark = ","
    )
  )

  message(
    "Protein-coding chr1-19/X genes retained for edgeR: ",
    format(
      nrow(retained_mapping),
      big.mark = ","
    )
  )

  group_table <- tibble::as_tibble(
    crescendo_object$group_table
  )

  required_group_columns <- c(
    "group_id",
    "sample_ID",
    "cluster_id",
    "n_spots"
  )

  missing_group_columns <- setdiff(
    required_group_columns,
    colnames(group_table)
  )

  if (length(missing_group_columns) > 0L) {
    stop(
      "Missing Crescendo group_table columns: ",
      paste(missing_group_columns, collapse = ", "),
      call. = FALSE
    )
  }

  group_metadata <- group_table |>
    dplyr::transmute(
      pseudobulk_id = as.character(.data$group_id),
      sample_ID = as.character(.data$sample_ID),
      cluster_id = as.character(.data$cluster_id),
      number_of_spots = as.integer(.data$n_spots)
    ) |>
    dplyr::filter(
      .data$sample_ID %in% sample_metadata$sample_ID
    ) |>
    dplyr::left_join(
      sample_metadata |>
        dplyr::mutate(
          fmt_donor_group =
            as.character(.data$fmt_donor_group),
          sex = as.character(.data$sex)
        ) |>
        dplyr::select(
          .data$sample_ID,
          .data$fmt_donor_group,
          .data$sex
        ),
      by = "sample_ID"
    )

  corrected_cluster_list <- vector(
    "list",
    length(cluster_ids)
  )

  names(corrected_cluster_list) <-
    cluster_ids

  corrected_by_cluster <- vector(
    "list",
    length(cluster_ids)
  )

  names(corrected_by_cluster) <-
    cluster_ids

  for (cluster_id_current in cluster_ids) {
    corrected_current <-
      corrected_counts_by_cluster[[cluster_id_current]]

    missing_retained_samples <- setdiff(
      sample_metadata$sample_ID,
      colnames(corrected_current)
    )

    if (length(missing_retained_samples) > 0L) {
      stop(
        "Corrected count matrix is missing retained sample(s) in cluster ",
        cluster_id_current,
        ": ",
        paste(missing_retained_samples, collapse = ", "),
        call. = FALSE
      )
    }

    if (any(excluded_sample_ids %in% sample_metadata$sample_ID)) {
      stop(
        "Excluded sample is still present in analysis metadata.",
        call. = FALSE
      )
    }

    corrected_current <- corrected_current[
      retained_mapping$feature_id,
      sample_metadata$sample_ID,
      drop = FALSE
    ]

    rownames(corrected_current) <-
      retained_mapping$ensembl_gene_id

    if (any(corrected_current < 0)) {
      stop(
        "Negative corrected counts in cluster ",
        cluster_id_current,
        ".",
        call. = FALSE
      )
    }

    group_current <- group_metadata |>
      dplyr::filter(
        .data$cluster_id ==
          .env$cluster_id_current
      )

    group_current <- group_current[
      match(
        sample_metadata$sample_ID,
        group_current$sample_ID
      ),
      ,
      drop = FALSE
    ]

    corrected_for_all <-
      corrected_current

    colnames(corrected_for_all) <-
      group_current$pseudobulk_id

    corrected_cluster_list[[cluster_id_current]] <- corrected_for_all

    corrected_by_cluster[[cluster_id_current]] <- corrected_current
  }

  pseudobulk_counts_all <- do.call(
    cbind,
    corrected_cluster_list
  )

  if (!inherits(
    pseudobulk_counts_all,
    "Matrix"
  )) {
    pseudobulk_counts_all <- Matrix::Matrix(
      pseudobulk_counts_all,
      sparse = TRUE
    )
  }

  group_metadata <- group_metadata[
    match(
      colnames(pseudobulk_counts_all),
      group_metadata$pseudobulk_id
    ),
    ,
    drop = FALSE
  ]

  list(
    sample_metadata = sample_metadata,
    cluster_ids = cluster_ids,
    pseudobulk_counts_all =
      pseudobulk_counts_all,
    corrected_by_cluster =
      corrected_by_cluster,
    group_metadata = group_metadata,
    gene_annotation = gene_annotation,
    feature_mapping =
      feature_mapping$mapping_status
  )
}

# ==============================================================================
# 12. LOAD RAW SPOT COUNTS FOR DESCRIPTIVE DETECTION METRICS
# ==============================================================================

prepare_raw_spot_detection_input <- function(
    target_gene_ids,
    sample_metadata,
    annotation_objects
) {
  loaded_seurat <- load_single_seurat_object_from_rdata(
    input_rdata_file = input_rdata_file,
    requested_object_name = NULL
  )

  seurat_object <- loaded_seurat$object
  seurat_metadata <- seurat_object[[]]

  cluster_column_resolved <- resolve_cluster_column(
    seurat_metadata = seurat_metadata,
    requested_cluster_column = cluster_column,
    expected_number_of_clusters =
      expected_number_of_clusters
  )

  count_extraction <- extract_raw_counts_from_seurat_layers(
    seurat_object = seurat_object,
    assay_name = assay_name
  )

  raw_counts <- count_extraction$counts

  raw_feature_mapping <- map_features_to_annotation(
    feature_names = rownames(raw_counts),
    annotation_all =
      annotation_objects$annotation_all,
    annotation_analysis =
      annotation_objects$annotation_analysis
  )

  retained_raw_mapping <-
    raw_feature_mapping$retained_mapping

  missing_target_ids <- setdiff(
    target_gene_ids,
    retained_raw_mapping$ensembl_gene_id
  )

  if (length(missing_target_ids) > 0L) {
    stop(
      "Raw Seurat counts are missing ",
      length(missing_target_ids),
      " corrected edgeR gene(s).",
      call. = FALSE
    )
  }

  raw_feature_ids_for_targets <-
    retained_raw_mapping$feature_id[
      match(
        target_gene_ids,
        retained_raw_mapping$ensembl_gene_id
      )
    ]

  raw_counts <- raw_counts[
    raw_feature_ids_for_targets,
    ,
    drop = FALSE
  ]

  rownames(raw_counts) <- target_gene_ids

  spot_metadata <- seurat_metadata[
    colnames(raw_counts),
    ,
    drop = FALSE
  ]

  spot_metadata$sample_ID_analysis <-
    as.character(
      spot_metadata[[sample_column]]
    )

  spot_metadata$cluster_id_analysis <-
    as.character(
      spot_metadata[[cluster_column_resolved]]
    )

  keep_spots <-
    spot_metadata$sample_ID_analysis %in%
      sample_metadata$sample_ID &
    !is.na(
      spot_metadata$cluster_id_analysis
    ) &
    spot_metadata$cluster_id_analysis != ""

  raw_counts <- raw_counts[
    ,
    keep_spots,
    drop = FALSE
  ]

  spot_metadata <- spot_metadata[
    keep_spots,
    ,
    drop = FALSE
  ]

  list(
    raw_counts = raw_counts,
    spot_metadata = spot_metadata,
    count_layers =
      count_extraction$count_layers,
    seurat_object_name =
      loaded_seurat$object_name
  )
}

# ==============================================================================
# 13. ONE-CLUSTER edgeR
# ==============================================================================

run_one_cluster <- function(
    cluster_id,
    prepared,
    raw_spot_input
) {
  cluster_group_metadata <-
    prepared$group_metadata |>
    dplyr::filter(
      .data$cluster_id ==
        .env$cluster_id
    ) |>
    dplyr::mutate(
      included_by_spot_threshold =
        .data$number_of_spots >=
          min_spots_per_sample_cluster
    ) |>
    dplyr::arrange(
      match(
        .data$sample_ID,
        prepared$sample_metadata$sample_ID
      )
    )

  included_metadata <-
    cluster_group_metadata |>
    dplyr::filter(
      .data$included_by_spot_threshold
    )

  factor_cell_counts <-
    included_metadata |>
    dplyr::count(
      .data$sex,
      .data$fmt_donor_group,
      name = "n_samples"
    )

  complete_factor_grid <- expand.grid(
    sex = c("Male", "Female"),
    fmt_donor_group =
      c("Neurotypical", "ASD"),
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(
      factor_cell_counts,
      by = c(
        "sex",
        "fmt_donor_group"
      )
    ) |>
    dplyr::mutate(
      n_samples = dplyr::coalesce(
        .data$n_samples,
        0L
      )
    )

  if (
    any(
      complete_factor_grid$n_samples <
        min_samples_per_factor_cell
    )
  ) {
    return(
      list(
        status = "skipped",
        reason = "insufficient_samples_in_factor_cell"
      )
    )
  }

  selected_group_ids <-
    included_metadata$pseudobulk_id

  cluster_counts <-
    prepared$pseudobulk_counts_all[
      ,
      selected_group_ids,
      drop = FALSE
    ]

  colnames(cluster_counts) <-
    included_metadata$sample_ID

  design_objects <-
    create_two_factor_design_and_contrasts(
      included_metadata
    )

  if (
    !identical(
      colnames(cluster_counts),
      design_objects$metadata$sample_ID
    )
  ) {
    stop(
      "Count/design sample order mismatch.",
      call. = FALSE
    )
  }

  dge_unfiltered <- edgeR::DGEList(
    counts = cluster_counts,
    samples = data.frame(
      sample_ID =
        design_objects$metadata$sample_ID,
      fmt_donor_group =
        as.character(
          design_objects$metadata$fmt_donor_group
        ),
      sex =
        as.character(
          design_objects$metadata$sex
        ),
      row.names =
        design_objects$metadata$sample_ID,
      stringsAsFactors = FALSE
    )
  )

  keep_genes <- edgeR::filterByExpr(
    dge_unfiltered,
    design = design_objects$design,
    min.count =
      filter_by_expr_parameters$min.count,
    min.total.count =
      filter_by_expr_parameters$min.total.count,
    large.n =
      filter_by_expr_parameters$large.n,
    min.prop =
      filter_by_expr_parameters$min.prop
  )

  if (!any(keep_genes)) {
    return(
      list(
        status = "skipped",
        reason = "no_gene_passed_filterByExpr"
      )
    )
  }

  filter_status <- tibble::tibble(
    ensembl_gene_id =
      rownames(dge_unfiltered$counts),
    total_count =
      rowSums(dge_unfiltered$counts),
    mean_count =
      rowMeans(dge_unfiltered$counts),
    samples_with_nonzero_count =
      rowSums(dge_unfiltered$counts > 0),
    passed_filterByExpr = keep_genes
  ) |>
    dplyr::left_join(
      prepared$gene_annotation,
      by = "ensembl_gene_id"
    ) |>
    dplyr::select(
      .data$ensembl_gene_id,
      .data$gene,
      .data$chromosome,
      .data$total_count,
      .data$mean_count,
      .data$samples_with_nonzero_count,
      .data$passed_filterByExpr
    )

  dge <- dge_unfiltered[
    keep_genes,
    ,
    keep.lib.sizes = FALSE
  ]

  dge <- edgeR::calcNormFactors(
    dge,
    method = "TMM"
  )

  tested_gene_ids <- rownames(
    dge$counts
  )

  detection_percentage_matrix <-
    compute_detection_percentages(
      raw_spot_count_matrix =
        raw_spot_input$raw_counts,
      raw_spot_metadata =
        raw_spot_input$spot_metadata,
      included_metadata =
        included_metadata,
      cluster_id = cluster_id,
      tested_gene_ids = tested_gene_ids,
      positive_spot_min_count =
        positive_spot_min_count
    )

  tmm_cpm_matrix <- edgeR::cpm(
    dge,
    normalized.lib.sizes = TRUE,
    log = FALSE
  )

  descriptive_tables <-
    make_sample_level_descriptive_tables(
      tmm_cpm_matrix =
        tmm_cpm_matrix,
      detection_percentage_matrix =
        detection_percentage_matrix,
      model_metadata =
        design_objects$metadata,
      test_definitions =
        design_objects$test_definitions
    )

  log2_qn_matrix <- create_log2qn_matrix(
    cluster_counts
  )

  log2_qn_tested <- log2_qn_matrix[
    tested_gene_ids,
    design_objects$metadata$sample_ID,
    drop = FALSE
  ]

  log2qn_tables <-
    make_log2qn_descriptive_tables(
      log2_qn_matrix =
        log2_qn_tested,
      model_metadata =
        design_objects$metadata,
      test_definitions =
        design_objects$test_definitions
    )

  dge <- edgeR::estimateDisp(
    dge,
    design = design_objects$design,
    robust = TRUE
  )

  fit <- edgeR::glmQLFit(
    dge,
    design = design_objects$design,
    robust = TRUE
  )

  test_objects <- lapply(
    design_objects$test_definitions$test_id,
    function(test_id) {
      edgeR::glmQLFTest(
        fit,
        contrast =
          design_objects$contrasts[
            ,
            test_id
          ]
      )
    }
  )

  names(test_objects) <-
    design_objects$test_definitions$test_id

  result_tables <- lapply(
    seq_len(
      nrow(
        design_objects$test_definitions
      )
    ),
    function(test_index) {
      test_id <-
        design_objects$test_definitions$test_id[[test_index]]

      create_result_table(
        qlf_test =
          test_objects[[test_id]],
        gene_annotation =
          prepared$gene_annotation,
        add_regulation =
          design_objects$test_definitions$add_regulation[[test_index]],
        descriptive_table =
          descriptive_tables[[test_id]],
        log2qn_table =
          log2qn_tables[[test_id]]
      )
    }
  )

  names(result_tables) <-
    design_objects$test_definitions$test_id

  analysis_summary <- dplyr::bind_rows(
    lapply(
      design_objects$test_definitions$test_id,
      function(test_id) {
        x <- result_tables[[test_id]]

        tibble::tibble(
          cluster_id =
            as.character(cluster_id),
          test_id = test_id,
          tested_genes = nrow(x),

          significant_P_0.01 =
            sum(
              x$PValue < 0.01,
              na.rm = TRUE
            ),

          significant_P_0.05 =
            sum(
              x$PValue < 0.05,
              na.rm = TRUE
            ),

          significant_FDR_0.01 =
            sum(
              x$FDR < 0.01,
              na.rm = TRUE
            ),

          significant_FDR_0.05 =
            sum(
              x$FDR < 0.05,
              na.rm = TRUE
            ),

          significant_FDR_0.10 =
            sum(
              x$FDR < 0.10,
              na.rm = TRUE
            ),

          significant_FDR_0.05_abs_logFC_ge_0.5 =
            sum(
              x$FDR < 0.05 &
                abs(x$logFC) >= 0.5,
              na.rm = TRUE
            ),

          significant_FDR_0.10_abs_logFC_ge_0.5 =
            sum(
              x$FDR < 0.10 &
                abs(x$logFC) >= 0.5,
              na.rm = TRUE
            )
        )
      }
    )
  )

  library_summary <- tibble::tibble(
    cluster_id =
      as.character(cluster_id),
    sample_ID =
      colnames(dge$counts),
    fmt_donor_group =
      as.character(
        design_objects$metadata$fmt_donor_group
      ),
    sex =
      as.character(
        design_objects$metadata$sex
      ),
    number_of_spots =
      included_metadata$number_of_spots[
        match(
          colnames(dge$counts),
          included_metadata$sample_ID
        )
      ],
    library_size =
      as.numeric(dge$samples$lib.size),
    normalization_factor =
      as.numeric(dge$samples$norm.factors),
    effective_library_size =
      as.numeric(
        dge$samples$lib.size *
          dge$samples$norm.factors
      )
  )

  list(
    status = "completed",
    reason = NA_character_,
    cluster_group_metadata =
      cluster_group_metadata,
    counts_before_filterByExpr =
      cluster_counts,
    log2_qn_matrix =
      log2_qn_matrix,
    keep_genes = keep_genes,
    filter_status = filter_status,
    design = design_objects$design,
    contrasts = design_objects$contrasts,
    test_definitions =
      design_objects$test_definitions,
    dge = dge,
    fit = fit,
    test_objects = test_objects,
    results = result_tables,
    detection_percentage_matrix =
      detection_percentage_matrix,
    tmm_cpm_matrix =
      tmm_cpm_matrix,
    analysis_summary =
      analysis_summary,
    library_summary =
      library_summary
  )
}

# ==============================================================================
# 14. VALIDATE INPUTS AND CREATE OUTPUT DIRECTORIES
# ==============================================================================

stage_message(
  1,
  "VALIDATING INPUTS"
)

check_required_packages()

required_files <- c(
  input_rdata_file,
  annotation_rds_file,
  annotation_metadata_file,
  crescendo_input_rds
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

output_subdirectories <- list(
  pseudobulk_inputs =
    file.path(
      output_statistics_dir,
      "01_pseudobulkInputs"
    ),
  design_metadata =
    file.path(
      output_statistics_dir,
      "02_designAndSampleMetadata"
    ),
  edger_results =
    file.path(
      output_statistics_dir,
      "03_edgeRResults"
    ),
  diagnostics =
    file.path(
      output_statistics_dir,
      "04_modelDiagnostics"
    ),
  objects =
    file.path(
      output_statistics_dir,
      "05_objects"
    ),
  logs =
    file.path(
      output_statistics_dir,
      "06_logs"
    )
)

per_cluster_counts_dir <- file.path(
  output_subdirectories$pseudobulk_inputs,
  "perCluster_crescendoCorrectedCounts"
)

per_cluster_qn_dir <- file.path(
  output_subdirectories$pseudobulk_inputs,
  "perCluster_crescendoCorrectedCounts_log2plus1_QN"
)

per_cluster_sample_metrics_dir <- file.path(
  output_subdirectories$pseudobulk_inputs,
  "perCluster_sampleLevelMetrics"
)

per_cluster_design_dir <- file.path(
  output_subdirectories$design_metadata,
  "perCluster"
)

per_cluster_xlsx_dir <- file.path(
  output_subdirectories$edger_results,
  "01_perCluster_XLSX"
)

combined_results_dir <- file.path(
  output_subdirectories$edger_results,
  "02_combinedResults"
)

per_cluster_diagnostics_dir <- file.path(
  output_subdirectories$diagnostics,
  "perCluster"
)

all_output_dirs <- c(
  output_statistics_dir,
  unname(output_subdirectories),
  per_cluster_counts_dir,
  per_cluster_qn_dir,
  per_cluster_sample_metrics_dir,
  per_cluster_design_dir,
  per_cluster_xlsx_dir,
  combined_results_dir,
  per_cluster_diagnostics_dir
)

invisible(
  lapply(
    all_output_dirs,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

message("SCRIPT BUILD: ", SCRIPT_BUILD)
message(
  "SENSITIVITY ANALYSIS — excluded sample(s): ",
  paste(excluded_sample_ids, collapse = ", ")
)
message("Expected samples in edgeR: ", expected_number_of_samples)
message("Crescendo input:")
message(crescendo_input_rds)
message("\nOutput:")
message(output_statistics_dir)

# ==============================================================================
# 15. PREPARE CRESCENDO PSEUDOBULK
# ==============================================================================

stage_message(
  2,
  "LOADING CRESCENDO-CORRECTED PSEUDOBULK COUNTS"
)

prepared <- prepare_crescendo_inputs()

annotation_objects <- read_and_validate_annotation(
  annotation_rds_file =
    annotation_rds_file,
  annotation_metadata_file =
    annotation_metadata_file
)

write_tsv(
  prepared$feature_mapping,
  file.path(
    output_subdirectories$pseudobulk_inputs,
    paste0(
      analysis_prefix,
      "_crescendoFeatureMappingStatus.tsv.gz"
    )
  )
)

# ==============================================================================
# 16. PREPARE RAW SPOT COUNTS FOR DESCRIPTIVE DETECTION METRICS
# ==============================================================================

stage_message(
  3,
  "LOADING RAW SPOT COUNTS FOR DESCRIPTIVE DETECTION METRICS"
)

raw_spot_input <- prepare_raw_spot_detection_input(
  target_gene_ids =
    prepared$gene_annotation$ensembl_gene_id,
  sample_metadata =
    prepared$sample_metadata,
  annotation_objects =
    annotation_objects
)

message(
  "Raw spot matrix used only for percent-positive metrics: ",
  nrow(raw_spot_input$raw_counts),
  " genes x ",
  ncol(raw_spot_input$raw_counts),
  " spots"
)

# ==============================================================================
# 17. RUN edgeR PER CLUSTER
# ==============================================================================

stage_message(
  4,
  "RUNNING edgeR PER CLUSTER"
)

cluster_ids <- prepared$cluster_ids

cluster_run_objects <- vector(
  "list",
  length(cluster_ids)
)

names(cluster_run_objects) <-
  cluster_ids

cluster_status_rows <- vector(
  "list",
  length(cluster_ids)
)

names(cluster_status_rows) <-
  cluster_ids

first_test_definitions <- NULL

for (cluster_index in seq_along(cluster_ids)) {
  cluster_id_current <-
    cluster_ids[[cluster_index]]

  message("\n------------------------------------------------------------")
  message(
    "[",
    cluster_index,
    "/",
    length(cluster_ids),
    "] Cluster ",
    cluster_id_current
  )
  message("------------------------------------------------------------")

  cluster_result <- tryCatch(
    run_one_cluster(
      cluster_id =
        cluster_id_current,
      prepared = prepared,
      raw_spot_input =
        raw_spot_input
    ),
    error = function(e) {
      list(
        status = "error",
        reason = conditionMessage(e)
      )
    }
  )

  cluster_run_objects[[cluster_id_current]] <- cluster_result

  cluster_status_rows[[cluster_id_current]] <- tibble::tibble(
    cluster_id =
      cluster_id_current,
    status =
      cluster_result$status,
    reason =
      cluster_result$reason
  )

  if (!identical(
    cluster_result$status,
    "completed"
  )) {
    message(
      "Cluster ",
      cluster_id_current,
      " not completed: ",
      cluster_result$status,
      " | ",
      cluster_result$reason
    )
    next
  }

  if (is.null(first_test_definitions)) {
    first_test_definitions <-
      cluster_result$test_definitions
  }

  cluster_slug <- sanitize_file_component(
    cluster_id_current
  )

  # Corrected pseudobulk counts before filterByExpr.
  corrected_count_output <- tibble::tibble(
    ensembl_gene_id =
      rownames(
        cluster_result$counts_before_filterByExpr
      )
  ) |>
    dplyr::left_join(
      prepared$gene_annotation |>
        dplyr::select(
          .data$ensembl_gene_id,
          .data$gene,
          .data$chromosome
        ),
      by = "ensembl_gene_id"
    ) |>
    dplyr::bind_cols(
      tibble::as_tibble(
        as.data.frame(
          cluster_result$counts_before_filterByExpr,
          check.names = FALSE
        ),
        .name_repair = "minimal"
      )
    )

  write_tsv(
    corrected_count_output,
    file.path(
      per_cluster_counts_dir,
      paste0(
        analysis_prefix,
        "_cluster_",
        cluster_slug,
        "_crescendoCorrectedCounts_proteinCodingAutosomalAndX.tsv.gz"
      )
    )
  )

  # log2+1 + quantile normalization.
  log2qn_output <- tibble::tibble(
    ensembl_gene_id =
      rownames(
        cluster_result$log2_qn_matrix
      )
  ) |>
    dplyr::left_join(
      prepared$gene_annotation |>
        dplyr::select(
          .data$ensembl_gene_id,
          .data$gene,
          .data$chromosome
        ),
      by = "ensembl_gene_id"
    ) |>
    dplyr::bind_cols(
      tibble::as_tibble(
        as.data.frame(
          cluster_result$log2_qn_matrix,
          check.names = FALSE
        ),
        .name_repair = "minimal"
      )
    )

  write_tsv(
    log2qn_output,
    file.path(
      per_cluster_qn_dir,
      paste0(
        analysis_prefix,
        "_cluster_",
        cluster_slug,
        "_crescendoCorrectedCounts_log2plus1_QN.tsv.gz"
      )
    )
  )

  # Sample-level TMM CPM.
  tmm_cpm_output <- tibble::tibble(
    ensembl_gene_id =
      rownames(
        cluster_result$tmm_cpm_matrix
      )
  ) |>
    dplyr::left_join(
      prepared$gene_annotation |>
        dplyr::select(
          .data$ensembl_gene_id,
          .data$gene,
          .data$chromosome
        ),
      by = "ensembl_gene_id"
    ) |>
    dplyr::bind_cols(
      tibble::as_tibble(
        as.data.frame(
          cluster_result$tmm_cpm_matrix,
          check.names = FALSE
        ),
        .name_repair = "minimal"
      )
    )

  write_tsv(
    tmm_cpm_output,
    file.path(
      per_cluster_sample_metrics_dir,
      paste0(
        analysis_prefix,
        "_cluster_",
        cluster_slug,
        "_sampleLevel_TMMnormalizedCPM.tsv.gz"
      )
    )
  )

  # Sample-level percent-positive.
  detection_output <- tibble::tibble(
    ensembl_gene_id =
      rownames(
        cluster_result$detection_percentage_matrix
      )
  ) |>
    dplyr::left_join(
      prepared$gene_annotation |>
        dplyr::select(
          .data$ensembl_gene_id,
          .data$gene,
          .data$chromosome
        ),
      by = "ensembl_gene_id"
    ) |>
    dplyr::bind_cols(
      tibble::as_tibble(
        as.data.frame(
          cluster_result$detection_percentage_matrix,
          check.names = FALSE
        ),
        .name_repair = "minimal"
      )
    )

  write_tsv(
    detection_output,
    file.path(
      per_cluster_sample_metrics_dir,
      paste0(
        analysis_prefix,
        "_cluster_",
        cluster_slug,
        "_sampleLevel_rawPercentPositiveSpots.tsv.gz"
      )
    )
  )

  # Design / contrasts / QC.
  design_output <- data.frame(
    sample_ID =
      rownames(
        cluster_result$design
      ),
    cluster_result$design,
    row.names = NULL,
    check.names = FALSE
  )

  contrast_output <- data.frame(
    coefficient =
      rownames(
        cluster_result$contrasts
      ),
    cluster_result$contrasts,
    row.names = NULL,
    check.names = FALSE
  )

  write_tsv(
    cluster_result$cluster_group_metadata,
    file.path(
      per_cluster_design_dir,
      paste0(
        analysis_prefix,
        "_cluster_",
        cluster_slug,
        "_sampleInclusion.tsv"
      )
    )
  )

  write_tsv(
    design_output,
    file.path(
      per_cluster_design_dir,
      paste0(
        analysis_prefix,
        "_cluster_",
        cluster_slug,
        "_designMatrix.tsv"
      )
    )
  )

  write_tsv(
    contrast_output,
    file.path(
      per_cluster_design_dir,
      paste0(
        analysis_prefix,
        "_cluster_",
        cluster_slug,
        "_contrastMatrix.tsv"
      )
    )
  )

  write_tsv(
    cluster_result$filter_status,
    file.path(
      per_cluster_diagnostics_dir,
      paste0(
        analysis_prefix,
        "_cluster_",
        cluster_slug,
        "_filterByExprStatus.tsv.gz"
      )
    )
  )

  write_tsv(
    cluster_result$library_summary,
    file.path(
      per_cluster_diagnostics_dir,
      paste0(
        analysis_prefix,
        "_cluster_",
        cluster_slug,
        "_librarySizeSummary.tsv"
      )
    )
  )

  write_tsv(
    cluster_result$analysis_summary,
    file.path(
      per_cluster_diagnostics_dir,
      paste0(
        analysis_prefix,
        "_cluster_",
        cluster_slug,
        "_edgeRAnalysisSummary.tsv"
      )
    )
  )

  cluster_xlsx_file <- file.path(
    per_cluster_xlsx_dir,
    paste0(
      analysis_prefix,
      "_cluster_",
      cluster_slug,
      "_fullResults_nineSheets.xlsx"
    )
  )

  write_result_workbook(
    result_list =
      cluster_result$results,
    test_definitions =
      cluster_result$test_definitions,
    output_file =
      cluster_xlsx_file,
    include_cluster_id = FALSE,
    xlsx_max_decimal_places =
      xlsx_max_decimal_places
  )

  message(
    "Completed cluster ",
    cluster_id_current,
    " | tested genes: ",
    nrow(cluster_result$dge$counts)
  )
}

# ==============================================================================
# 18. COMBINE RESULTS
# ==============================================================================

stage_message(
  5,
  "COMBINING edgeR RESULTS"
)

cluster_status <- dplyr::bind_rows(
  cluster_status_rows
)

cluster_status_file <- file.path(
  output_subdirectories$diagnostics,
  paste0(
    analysis_prefix,
    "_clusterAnalysisStatus.tsv"
  )
)

write_tsv(
  cluster_status,
  cluster_status_file
)

completed_cluster_ids <- cluster_status |>
  dplyr::filter(
    .data$status == "completed"
  ) |>
  dplyr::pull(
    .data$cluster_id
  ) |>
  as.character()

if (length(completed_cluster_ids) == 0L) {
  stop(
    "No cluster completed edgeR.",
    call. = FALSE
  )
}

combined_results <- lapply(
  first_test_definitions$test_id,
  function(test_id) {
    combined <- dplyr::bind_rows(
      lapply(
        completed_cluster_ids,
        function(cluster_id_current) {
          cluster_run_objects[[cluster_id_current]]$results[[test_id]] |>
            dplyr::mutate(
              cluster_id =
                cluster_id_current,
              .before = 1
            )
        }
      )
    )

    combined |>
      dplyr::mutate(
        cluster_id = factor(
          .data$cluster_id,
          levels = cluster_ids
        )
      ) |>
      dplyr::arrange(
        .data$cluster_id,
        .data$PValue,
        dplyr::desc(
          abs(.data$logFC)
        )
      ) |>
      dplyr::mutate(
        cluster_id =
          as.character(.data$cluster_id)
      )
  }
)

names(combined_results) <-
  first_test_definitions$test_id

combined_xlsx_file <- file.path(
  combined_results_dir,
  paste0(
    analysis_prefix,
    "_allClusters_fullResults_nineSheets.xlsx"
  )
)

write_result_workbook(
  result_list = combined_results,
  test_definitions =
    first_test_definitions,
  output_file =
    combined_xlsx_file,
  include_cluster_id = TRUE,
  xlsx_max_decimal_places =
    xlsx_max_decimal_places
)

for (test_id in first_test_definitions$test_id) {
  write_tsv(
    combined_results[[test_id]],
    file.path(
      combined_results_dir,
      paste0(
        analysis_prefix,
        "_allClusters_",
        test_id,
        "_fullResults.tsv.gz"
      )
    )
  )
}

combined_analysis_summary <- dplyr::bind_rows(
  lapply(
    completed_cluster_ids,
    function(cluster_id_current) {
      cluster_run_objects[[cluster_id_current]]$analysis_summary
    }
  )
)

combined_analysis_summary_file <- file.path(
  output_subdirectories$diagnostics,
  paste0(
    analysis_prefix,
    "_edgeRAnalysisSummary_allClusters.tsv"
  )
)

write_tsv(
  combined_analysis_summary,
  combined_analysis_summary_file
)

# ==============================================================================
# 19. UNIQUE SIGNIFICANT GENES ACROSS CLUSTERS
# ==============================================================================

stage_message(
  6,
  "SUMMARIZING UNIQUE SIGNIFICANT GENES"
)

unique_gene_summary <- dplyr::bind_rows(
  lapply(
    first_test_definitions$test_id,
    function(test_id) {
      x <- combined_results[[test_id]]

      tibble::tibble(
        test_id = test_id,

        unique_genes_tested =
          dplyr::n_distinct(
            x$ensembl_gene_id
          ),

        unique_P_lt_0.01 =
          dplyr::n_distinct(
            x$ensembl_gene_id[
              x$PValue < 0.01
            ]
          ),

        unique_P_lt_0.05 =
          dplyr::n_distinct(
            x$ensembl_gene_id[
              x$PValue < 0.05
            ]
          ),

        unique_FDR_lt_0.01 =
          dplyr::n_distinct(
            x$ensembl_gene_id[
              x$FDR < 0.01
            ]
          ),

        unique_FDR_lt_0.05 =
          dplyr::n_distinct(
            x$ensembl_gene_id[
              x$FDR < 0.05
            ]
          ),

        unique_FDR_lt_0.10 =
          dplyr::n_distinct(
            x$ensembl_gene_id[
              x$FDR < 0.10
            ]
          ),

        unique_FDR_lt_0.01_abs_logFC_ge_0.5 =
          dplyr::n_distinct(
            x$ensembl_gene_id[
              x$FDR < 0.01 &
                abs(x$logFC) >= 0.5
            ]
          ),

        unique_FDR_lt_0.05_abs_logFC_ge_0.5 =
          dplyr::n_distinct(
            x$ensembl_gene_id[
              x$FDR < 0.05 &
                abs(x$logFC) >= 0.5
            ]
          ),

        unique_FDR_lt_0.10_abs_logFC_ge_0.5 =
          dplyr::n_distinct(
            x$ensembl_gene_id[
              x$FDR < 0.10 &
                abs(x$logFC) >= 0.5
            ]
          )
      )
    }
  )
)

unique_gene_summary_file <- file.path(
  output_subdirectories$diagnostics,
  paste0(
    analysis_prefix,
    "_uniqueSignificantGenesAcrossClusters.tsv"
  )
)

write_tsv(
  unique_gene_summary,
  unique_gene_summary_file
)

print(
  unique_gene_summary,
  n = Inf,
  width = Inf
)

# ==============================================================================
# 20. SAVE R OBJECTS
# ==============================================================================

stage_message(
  7,
  "SAVING edgeR OBJECTS"
)

edgeR_perCluster_results <- lapply(
  completed_cluster_ids,
  function(cluster_id_current) {
    cluster_run_objects[[cluster_id_current]]$results
  }
)

names(edgeR_perCluster_results) <-
  completed_cluster_ids

edgeR_perCluster_combinedResults <-
  combined_results

edgeR_perCluster_analysisSummary <-
  combined_analysis_summary

edgeR_perCluster_clusterStatus <-
  cluster_status

edgeR_perCluster_testDefinitions <-
  first_test_definitions

result_rdata_file <- file.path(
  output_subdirectories$edger_results,
  paste0(
    analysis_prefix,
    "_edgeRResults.RData"
  )
)

save(
  edgeR_perCluster_results,
  edgeR_perCluster_combinedResults,
  edgeR_perCluster_analysisSummary,
  edgeR_perCluster_clusterStatus,
  edgeR_perCluster_testDefinitions,
  file = result_rdata_file,
  compress = "xz"
)

edgeR_perCluster_modelObjects <- lapply(
  completed_cluster_ids,
  function(cluster_id_current) {
    cluster_object <-
      cluster_run_objects[[cluster_id_current]]

    list(
      dge = cluster_object$dge,
      fit = cluster_object$fit,
      test_objects =
        cluster_object$test_objects,
      design =
        cluster_object$design,
      contrasts =
        cluster_object$contrasts,
      keep_genes =
        cluster_object$keep_genes
    )
  }
)

names(edgeR_perCluster_modelObjects) <-
  completed_cluster_ids

model_rdata_file <- file.path(
  output_subdirectories$objects,
  paste0(
    analysis_prefix,
    "_edgeRModelObjects.RData"
  )
)

save(
  edgeR_perCluster_modelObjects,
  file = model_rdata_file,
  compress = "xz"
)

# ==============================================================================
# 21. LOG ANALYSIS SETTINGS
# ==============================================================================

stage_message(
  8,
  "WRITING ANALYSIS LOGS"
)

write_tsv(
  first_test_definitions,
  file.path(
    output_subdirectories$design_metadata,
    paste0(
      analysis_prefix,
      "_testDefinitions.tsv"
    )
  )
)

analysis_parameters <- c(
  paste0(
    "generated=",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  paste0(
    "script_build=",
    SCRIPT_BUILD
  ),
  paste0(
    "excluded_sample_ids=",
    paste(excluded_sample_ids, collapse = ",")
  ),
  paste0(
    "number_of_samples_used=",
    expected_number_of_samples
  ),
  "crescendo_recalculated_after_sample_exclusion=FALSE",
  paste0(
    "crescendo_input_rds=",
    crescendo_input_rds
  ),
  "edgeR_count_input=crescendo_corrected_pseudobulk_counts",
  "raw_spot_counts_used_for_edgeR=FALSE",
  "raw_spot_counts_used_for_percent_positive_metrics=TRUE",
  "log2QN_used_for_edgeR=FALSE",
  "model=expression ~ fmt_donor_group * sex",
  "filterByExpr=TRUE",
  paste0(
    "filterByExpr_min.count=",
    filter_by_expr_parameters$min.count
  ),
  paste0(
    "filterByExpr_min.total.count=",
    filter_by_expr_parameters$min.total.count
  ),
  paste0(
    "filterByExpr_large.n=",
    filter_by_expr_parameters$large.n
  ),
  paste0(
    "filterByExpr_min.prop=",
    filter_by_expr_parameters$min.prop
  ),
  "TMM=TRUE",
  "estimateDisp_robust=TRUE",
  "glmQLFit_robust=TRUE",
  "number_of_tests_per_cluster=9",
  "gene_universe=Ensembl115 protein_coding chromosomes 1-19,X",
  paste0(
    "completed_clusters=",
    paste(
      completed_cluster_ids,
      collapse = ","
    )
  ),
  paste0(
    "combined_xlsx=",
    combined_xlsx_file
  ),
  paste0(
    "unique_gene_summary=",
    unique_gene_summary_file
  ),
  paste0(
    "result_rdata=",
    result_rdata_file
  ),
  paste0(
    "model_rdata=",
    model_rdata_file
  )
)

writeLines(
  analysis_parameters,
  con = file.path(
    output_subdirectories$logs,
    paste0(
      analysis_prefix,
      "_analysisParameters.txt"
    )
  )
)

writeLines(
  capture.output(
    utils::sessionInfo()
  ),
  con = file.path(
    output_subdirectories$logs,
    paste0(
      analysis_prefix,
      "_sessionInfo.txt"
    )
  )
)

# ==============================================================================
# 22. FINISH
# ==============================================================================

stage_message(
  9,
  "CRESCENDO edgeR ANALYSIS FINISHED"
)

message(
  "Clusters completed: ",
  length(completed_cluster_ids),
  "/",
  length(cluster_ids)
)

message(
  "Samples used in edgeR: ",
  expected_number_of_samples,
  " | excluded: ",
  paste(excluded_sample_ids, collapse = ", ")
)

message("\nCombined XLSX:")
message(combined_xlsx_file)

message("\nUnique significant genes summary:")
message(unique_gene_summary_file)

message("\nedgeR results RData:")
message(result_rdata_file)

message("\nedgeR model RData:")
message(model_rdata_file)

message("\nOutput directory:")
message(output_statistics_dir)

error_clusters <- cluster_status |>
  dplyr::filter(
    .data$status == "error"
  )

if (nrow(error_clusters) > 0L) {
  stop(
    "At least one cluster failed. See clusterAnalysisStatus.tsv.",
    call. = FALSE
  )
}
