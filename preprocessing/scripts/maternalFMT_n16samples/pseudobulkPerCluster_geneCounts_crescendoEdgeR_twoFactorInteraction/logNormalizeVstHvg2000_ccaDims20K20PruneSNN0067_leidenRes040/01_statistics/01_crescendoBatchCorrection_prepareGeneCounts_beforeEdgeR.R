#!/usr/bin/env Rscript
# ==============================================================================
# 01_crescendoBatchCorrection_prepareGeneCounts_beforeEdgeR.R
#
# PURPOSE
#   1. Load the same Seurat object and sample metadata used by the existing
#      pseudobulk-per-cluster edgeR workflow.
#   2. Keep the same 16 biological samples and the existing Leiden clusters.
#   3. Extract raw RNA spot counts.
#   4. Build Harmony soft-cluster assignments required by Crescendo.
#   5. Run Crescendo on raw spot-level counts using Slide ID as the batch.
#   6. Sum Crescendo-corrected spot counts to sample x cluster pseudobulk counts.
#   7. Save corrected pseudobulk matrices ready for a later edgeR analysis.
#
# IMPORTANT
#   - Crescendo is run on ALL retained spots together, not separately per cluster.
#   - Existing C1-C16 / leiden_res0.40 labels are NOT changed.
#   - Existing clusters are used only for the final pseudobulk SUM after Crescendo.
#   - Crescendo is run in gene chunks because its intermediate matrices are dense.
#   - Crescendo models each gene independently, so chunking does not change the
#     gene-wise model.
#   - Slide ID is partially confounded with biological group/sex in this dataset.
#     This workflow must therefore be interpreted as a conservative sensitivity /
#     rescue analysis for slide-associated technical variation.
# ==============================================================================

message("SCRIPT VERSION: v13 - fresh full run + dense-safe Matching -1 repair")

options(stringsAsFactors = FALSE)
set.seed(20260831)

# ==============================================================================
# 1. Project and analysis configuration
# ==============================================================================

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"

clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

analysis_name <-
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction"

cluster_column <- "leiden_res0.40"
sample_column <- "sample_ID"
assay_name <- "RNA"

excluded_samples <- c(
  "20_1F",
  "12_3F",
  "15_1M",
  "20_3M"
)

expected_number_of_samples <- 16L
expected_number_of_clusters <- 16L

# Crescendo / Harmony settings.
harmony_preferred_reductions <- c(
  "pca",
  "integrated.cca",
  "cca"
)
harmony_n_dims <- 20L
harmony_nclust <- 20L
harmony_theta <- 2
harmony_verbose <- TRUE

# Crescendo estimation is downsampled for computational stability.
# In the current Crescendo implementation any prop != 1 activates the package's
# batch/cluster-aware downsampling, whose internal default is prop = 0.1.
crescendo_prop <- 0.1
crescendo_min_cells <- 50L
crescendo_lambda <- NULL
crescendo_alpha <- 0
crescendo_seed <- 20260831L
crescendo_verbose <- TRUE

# Very sparse genes can make Poisson glmnet unstable.
# Keep a gene only if it is detected (>0 raw count) in at least this fraction
# of spots in AT LEAST ONE existing Leiden cluster.
# This preserves cluster-specific genes while removing extremely sparse genes.
gene_min_detected_fraction_within_any_cluster <- 0.10

detected_cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cores) || detected_cores < 2L) {
  detected_cores <- 2L
}
crescendo_mc_cores <- min(40L, max(1L, detected_cores - 2L))

# Gene chunk size.
gene_chunk_size <- 100L

# FALSE = full analysis. TRUE = small technical test run.
test_mode <- FALSE
test_number_of_genes <- 20L

# If TRUE, completed chunk files are reused.
resume_existing_chunks <- FALSE

# ==============================================================================
# 2. Input files
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

metadata_file <- file.path(
  project_root,
  "data",
  "metadata_autismFMT.tsv"
)

annotation_rds_file <- file.path(
  project_root,
  "data",
  "ensembl115",
  "ensembl115_mouse_proteinCodingGenes.rds"
)

# ==============================================================================
# 3. Output directories
# ==============================================================================

output_statistics_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  analysis_name,
  clustering_name,
  "01_statistics"
)

output_crescendo_dir <- file.path(
  output_statistics_dir,
  "01_crescendoBatchCorrection"
)

output_metadata_dir <- file.path(
  output_crescendo_dir,
  "01_metadata"
)

output_harmony_dir <- file.path(
  output_crescendo_dir,
  "02_harmony"
)

output_corrected_spot_chunks_dir <- file.path(
  output_crescendo_dir,
  "03_correctedSpotCounts_chunks_v13"
)

output_pseudobulk_chunks_dir <- file.path(
  output_crescendo_dir,
  "04_pseudobulkChunks_v13"
)

output_pseudobulk_dir <- file.path(
  output_crescendo_dir,
  "05_pseudobulkCorrectedCounts"
)

output_pseudobulk_per_cluster_dir <- file.path(
  output_pseudobulk_dir,
  "01_perCluster"
)

output_qc_dir <- file.path(
  output_crescendo_dir,
  "06_qc"
)

output_object_dir <- file.path(
  output_crescendo_dir,
  "07_objects"
)

all_output_dirs <- c(
  output_statistics_dir,
  output_crescendo_dir,
  output_metadata_dir,
  output_harmony_dir,
  output_corrected_spot_chunks_dir,
  output_pseudobulk_chunks_dir,
  output_pseudobulk_dir,
  output_pseudobulk_per_cluster_dir,
  output_qc_dir,
  output_object_dir
)

for (dir_current in all_output_dirs) {
  dir.create(
    dir_current,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

analysis_prefix <- paste0(
  dataset_name,
  "_",
  clustering_name,
  "_crescendoBatchCorrection"
)

# ==============================================================================
# 4. Required packages
# ==============================================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "harmony",
  "crescendo"
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
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      "\n\nInstallation examples:\n",
      "install.packages(c(\"Seurat\", \"harmony\"))\n",
      "install.packages(\"remotes\")\n",
      "remotes::install_github(\"immunogenomics/crescendo\")\n"
    ),
    call. = FALSE
  )
}

# ==============================================================================
# 5. Helper functions
# ==============================================================================

stage_message <- function(stage_number, text) {
  message("\n============================================================")
  message("STAGE ", stage_number, ": ", text)
  message("============================================================")
}

write_tsv_gz <- function(x, output_file) {
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

sanitize_file_component <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x
}

natural_order <- function(x) {
  x <- as.character(x)
  numeric_part <- suppressWarnings(
    as.integer(
      sub("^.*?([0-9]+)$", "\\1", x)
    )
  )
  if (all(!is.na(numeric_part))) {
    return(x[order(numeric_part)])
  }
  sort(x)
}

detect_column <- function(
  data,
  candidates,
  displayed_name,
  required = TRUE
) {
  exact_match <- candidates[candidates %in% colnames(data)]
  if (length(exact_match) > 0L) {
    return(exact_match[[1L]])
  }

  normalized_names <- tolower(
    gsub("[^a-z0-9]+", "", colnames(data))
  )
  normalized_candidates <- tolower(
    gsub("[^a-z0-9]+", "", candidates)
  )

  for (candidate_index in seq_along(normalized_candidates)) {
    match_index <- which(
      normalized_names == normalized_candidates[[candidate_index]]
    )
    if (length(match_index) > 0L) {
      return(colnames(data)[match_index[[1L]]])
    }
  }

  if (required) {
    stop(
      "Could not detect ",
      displayed_name,
      " column. Available columns:\n",
      paste(colnames(data), collapse = "\n"),
      call. = FALSE
    )
  }

  NULL
}

find_seurat_object <- function(
  input_rdata_file,
  cluster_column
) {
  load_environment <- new.env(parent = emptyenv())

  loaded_object_names <- load(
    input_rdata_file,
    envir = load_environment
  )

  seurat_object_names <- loaded_object_names[
    vapply(
      loaded_object_names,
      function(object_name) {
        inherits(
          get(
            object_name,
            envir = load_environment,
            inherits = FALSE
          ),
          "Seurat"
        )
      },
      logical(1)
    )
  ]

  if (length(seurat_object_names) == 0L) {
    stop(
      "No Seurat object was found in:\n",
      input_rdata_file,
      call. = FALSE
    )
  }

  objects_with_cluster <- seurat_object_names[
    vapply(
      seurat_object_names,
      function(object_name) {
        object_current <- get(
          object_name,
          envir = load_environment,
          inherits = FALSE
        )
        cluster_column %in% colnames(object_current@meta.data)
      },
      logical(1)
    )
  ]

  if (length(objects_with_cluster) == 1L) {
    selected_object_name <- objects_with_cluster[[1L]]
  } else if (length(objects_with_cluster) > 1L) {
    object_sizes <- vapply(
      objects_with_cluster,
      function(object_name) {
        ncol(
          get(
            object_name,
            envir = load_environment,
            inherits = FALSE
          )
        )
      },
      integer(1)
    )
    selected_object_name <- objects_with_cluster[
      which.max(object_sizes)
    ][[1L]]
    warning(
      "Multiple Seurat objects contain cluster column '",
      cluster_column,
      "'. Using the object with the largest number of spots: ",
      selected_object_name,
      call. = FALSE
    )
  } else if (length(seurat_object_names) == 1L) {
    selected_object_name <- seurat_object_names[[1L]]
  } else {
    stop(
      "Multiple Seurat objects were found, but none could be uniquely ",
      "selected using cluster column '",
      cluster_column,
      "'. Objects: ",
      paste(seurat_object_names, collapse = ", "),
      call. = FALSE
    )
  }

  list(
    object_name = selected_object_name,
    object = get(
      selected_object_name,
      envir = load_environment,
      inherits = FALSE
    )
  )
}

get_raw_count_matrix <- function(
  seurat_object,
  assay_name
) {
  if (!(assay_name %in% names(seurat_object@assays))) {
    stop(
      "Assay '",
      assay_name,
      "' is missing from the Seurat object.",
      call. = FALSE
    )
  }

  assay_object <- seurat_object[[assay_name]]

  # --------------------------------------------------------------------------
  # Seurat v5 Assay5:
  # Raw counts may be stored as one layer per sample:
  # counts.sample1, counts.sample2, ...
  #
  # Calling LayerData(..., layer = "counts") directly in that situation uses
  # only the first matching layer. Therefore the split count layers must be
  # joined before raw counts are extracted.
  # --------------------------------------------------------------------------
  if (inherits(assay_object, "Assay5")) {
    all_layers <- SeuratObject::Layers(assay_object)

    count_layers <- all_layers[
      all_layers == "counts" |
        startsWith(all_layers, "counts.")
    ]

    if (length(count_layers) == 0L) {
      stop(
        "No raw-count layer named 'counts' or starting with 'counts.' was found in assay '",
        assay_name,
        "'. Available layers:
",
        paste(all_layers, collapse = "
"),
        call. = FALSE
      )
    }

    message(
      "Raw-count layers detected in Seurat v5 assay: ",
      length(count_layers)
    )
    message(
      paste(count_layers, collapse = "
")
    )

    if (length(count_layers) > 1L) {
      message(
        "Joining all split counts.* layers into one 'counts' layer..."
      )

      assay_object <- SeuratObject::JoinLayers(
        object = assay_object
      )

      joined_layers <- SeuratObject::Layers(assay_object)

      message(
        "Layers after JoinLayers(): ",
        paste(joined_layers, collapse = ", ")
      )
    }

    count_matrix <- SeuratObject::LayerData(
      object = assay_object,
      layer = "counts"
    )

    count_matrix <- methods::as(
      count_matrix,
      "dgCMatrix"
    )

    message(
      "Joined raw-count matrix: ",
      nrow(count_matrix),
      " genes x ",
      ncol(count_matrix),
      " spots"
    )

    expected_spots <- colnames(seurat_object)

    missing_spots <- setdiff(
      expected_spots,
      colnames(count_matrix)
    )

    extra_spots <- setdiff(
      colnames(count_matrix),
      expected_spots
    )

    if (length(missing_spots) > 0L) {
      stop(
        "Joined raw-count matrix is missing ",
        length(missing_spots),
        " retained Seurat spots.",
        call. = FALSE
      )
    }

    if (length(extra_spots) > 0L) {
      message(
        "Removing ",
        length(extra_spots),
        " count-matrix spots not present in the retained Seurat object."
      )
    }

    count_matrix <- count_matrix[
      ,
      expected_spots,
      drop = FALSE
    ]

    if (ncol(count_matrix) != ncol(seurat_object)) {
      stop(
        "After reordering, raw-count matrix has ",
        ncol(count_matrix),
        " spots but Seurat has ",
        ncol(seurat_object),
        ".",
        call. = FALSE
      )
    }

    return(count_matrix)
  }

  # --------------------------------------------------------------------------
  # Seurat v3/v4 fallback.
  # --------------------------------------------------------------------------
  count_matrix <- tryCatch(
    {
      Seurat::GetAssayData(
        object = seurat_object,
        assay = assay_name,
        slot = "counts"
      )
    },
    error = function(e) {
      NULL
    }
  )

  if (!is.null(count_matrix)) {
    return(
      methods::as(
        count_matrix,
        "dgCMatrix"
      )
    )
  }

  stop(
    "Could not extract raw counts from assay '",
    assay_name,
    "'.",
    call. = FALSE
  )
}

select_harmony_reduction <- function(
  seurat_object,
  preferred_reductions
) {
  available_reductions <- Seurat::Reductions(seurat_object)

  selected <- preferred_reductions[
    preferred_reductions %in% available_reductions
  ]

  if (length(selected) == 0L) {
    stop(
      "No suitable reduction for Harmony was found.\n",
      "Preferred: ",
      paste(preferred_reductions, collapse = ", "),
      "\nAvailable: ",
      paste(available_reductions, collapse = ", "),
      call. = FALSE
    )
  }

  selected[[1L]]
}

extract_reduction_embeddings <- function(
  seurat_object,
  reduction_name,
  spot_ids,
  n_dims
) {
  embeddings <- Seurat::Embeddings(
    seurat_object,
    reduction = reduction_name
  )

  missing_spots <- setdiff(
    spot_ids,
    rownames(embeddings)
  )

  if (length(missing_spots) > 0L) {
    stop(
      "Reduction '",
      reduction_name,
      "' is missing ",
      length(missing_spots),
      " retained spot(s).",
      call. = FALSE
    )
  }

  embeddings <- embeddings[
    spot_ids,
    ,
    drop = FALSE
  ]

  n_dims_use <- min(
    n_dims,
    ncol(embeddings)
  )

  embeddings[
    ,
    seq_len(n_dims_use),
    drop = FALSE
  ]
}

strip_ensembl_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

load_and_match_protein_coding_genes <- function(
  annotation_rds_file,
  count_features
) {
  annotation_object <- readRDS(
    annotation_rds_file
  )

  message(
    "Protein-coding RDS class: ",
    paste(class(annotation_object), collapse = ", ")
  )

  candidate_vectors <- list()

  add_candidate <- function(name, values) {
    if (is.null(values)) {
      return(invisible(NULL))
    }

    if (is.factor(values)) {
      values <- as.character(values)
    }

    if (
      is.character(values) ||
        is.integer(values) ||
        is.numeric(values)
    ) {
      values <- as.character(values)
      values <- values[
        !is.na(values) &
          nzchar(values)
      ]

      if (length(values) > 0L) {
        candidate_vectors[[name]] <<- values
      }
    }

    invisible(NULL)
  }

  if (
    is.atomic(annotation_object) &&
      is.null(dim(annotation_object))
  ) {
    add_candidate(
      "atomic_vector",
      annotation_object
    )
  }

  if (
    is.data.frame(annotation_object) ||
      is.matrix(annotation_object)
  ) {
    annotation_df <- as.data.frame(
      annotation_object,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    for (column_name in colnames(annotation_df)) {
      add_candidate(
        paste0("column:", column_name),
        annotation_df[[column_name]]
      )
    }

    annotation_rownames <- rownames(annotation_df)

    if (
      !is.null(annotation_rownames) &&
        !identical(
          annotation_rownames,
          as.character(seq_len(nrow(annotation_df)))
        )
    ) {
      add_candidate(
        "rownames",
        annotation_rownames
      )
    }
  }

  if (is.list(annotation_object)) {
    list_names <- names(annotation_object)

    if (is.null(list_names)) {
      list_names <- paste0(
        "element",
        seq_along(annotation_object)
      )
    }

    for (i in seq_along(annotation_object)) {
      element <- annotation_object[[i]]
      element_name <- list_names[[i]]

      if (
        is.atomic(element) &&
          is.null(dim(element))
      ) {
        add_candidate(
          paste0("list:", element_name),
          element
        )
      } else if (
        is.data.frame(element) ||
          is.matrix(element)
      ) {
        element_df <- as.data.frame(
          element,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )

        for (column_name in colnames(element_df)) {
          add_candidate(
            paste0(
              "list:",
              element_name,
              "$",
              column_name
            ),
            element_df[[column_name]]
          )
        }

        element_rownames <- rownames(element_df)

        if (
          !is.null(element_rownames) &&
            !identical(
              element_rownames,
              as.character(seq_len(nrow(element_df)))
            )
        ) {
          add_candidate(
            paste0(
              "list:",
              element_name,
              "$rownames"
            ),
            element_rownames
          )
        }
      }
    }
  }

  if (length(candidate_vectors) == 0L) {
    stop(
      "Could not extract any candidate gene-ID vector from protein-coding RDS.",
      call. = FALSE
    )
  }

  count_features <- as.character(
    count_features
  )

  count_features_versionless <- strip_ensembl_version(
    count_features
  )

  candidate_scores <- lapply(
    names(candidate_vectors),
    function(candidate_name) {
      candidate_values <- unique(
        candidate_vectors[[candidate_name]]
      )

      candidate_values_versionless <- unique(
        strip_ensembl_version(
          candidate_values
        )
      )

      exact_matches <- sum(
        count_features %in% candidate_values
      )

      versionless_matches <- sum(
        count_features_versionless %in%
          candidate_values_versionless
      )

      data.frame(
        candidate = candidate_name,
        n_candidate_values = length(candidate_values),
        exact_matches = exact_matches,
        versionless_matches = versionless_matches,
        best_matches = max(
          exact_matches,
          versionless_matches
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  candidate_scores <- do.call(
    rbind,
    candidate_scores
  )

  candidate_scores <- candidate_scores[
    order(
      -candidate_scores$best_matches,
      -candidate_scores$exact_matches
    ),
    ,
    drop = FALSE
  ]

  best_candidate_name <- candidate_scores$candidate[[1L]]
  best_candidate_values <- unique(
    candidate_vectors[[best_candidate_name]]
  )

  use_versionless <- candidate_scores$versionless_matches[[1L]] >
    candidate_scores$exact_matches[[1L]]

  if (use_versionless) {
    keep <- count_features_versionless %in%
      unique(
        strip_ensembl_version(
          best_candidate_values
        )
      )
    matching_mode <- "version-stripped"
  } else {
    keep <- count_features %in%
      best_candidate_values
    matching_mode <- "exact"
  }

  matched_features <- count_features[
    keep
  ]

  if (length(matched_features) < 1000L) {
    stop(
      "Only ",
      length(matched_features),
      " raw-count features matched the protein-coding RDS. ",
      "Best candidate: ",
      best_candidate_name,
      " (",
      matching_mode,
      "). This is unexpectedly low.",
      call. = FALSE
    )
  }

  message(
    "Protein-coding matching source: ",
    best_candidate_name
  )
  message(
    "Protein-coding matching mode: ",
    matching_mode
  )
  message(
    "Protein-coding genes matched to raw counts: ",
    format(
      length(matched_features),
      big.mark = ","
    )
  )

  list(
    genes = matched_features,
    candidate_scores = candidate_scores,
    best_candidate = best_candidate_name,
    matching_mode = matching_mode
  )
}

select_protein_coding_autosomal_x_genes <- function(
  count_features,
  annotation_metadata_file
) {
  annotation <- utils::read.delim(
    annotation_metadata_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (nrow(annotation) == 0L) {
    stop(
      "Annotation metadata file is empty:\n",
      annotation_metadata_file,
      call. = FALSE
    )
  }

  chromosome_column <- detect_column(
    annotation,
    candidates = c(
      "chromosome",
      "chromosome_name",
      "chromosome/scaffold name",
      "chr"
    ),
    displayed_name = "chromosome"
  )

  ensembl_column <- detect_column(
    annotation,
    candidates = c(
      "ensembl_gene_id",
      "ensembl gene id",
      "gene_id",
      "ensembl_id"
    ),
    displayed_name = "Ensembl gene ID",
    required = FALSE
  )

  gene_symbol_column <- detect_column(
    annotation,
    candidates = c(
      "gene",
      "gene_symbol",
      "external_gene_name",
      "mgi_symbol",
      "symbol"
    ),
    displayed_name = "gene symbol",
    required = FALSE
  )

  chromosome_normalized <- sub(
    "^chr",
    "",
    as.character(annotation[[chromosome_column]]),
    ignore.case = TRUE
  )

  keep_chromosome <- chromosome_normalized %in% c(
    as.character(1:19),
    "X"
  )

  annotation <- annotation[
    keep_chromosome,
    ,
    drop = FALSE
  ]

  if (nrow(annotation) == 0L) {
    stop(
      "No autosomal/X genes were retained from annotation metadata.",
      call. = FALSE
    )
  }

  count_features_character <- as.character(count_features)

  matching_strategies <- list()

  if (!is.null(ensembl_column)) {
    annotation_ensembl <- as.character(annotation[[ensembl_column]])

    matching_strategies[["ensembl_exact"]] <- list(
      annotation_key = annotation_ensembl,
      count_key = count_features_character
    )

    matching_strategies[["ensembl_versionless"]] <- list(
      annotation_key = strip_ensembl_version(annotation_ensembl),
      count_key = strip_ensembl_version(count_features_character)
    )
  }

  if (!is.null(gene_symbol_column)) {
    matching_strategies[["gene_symbol_exact"]] <- list(
      annotation_key = as.character(annotation[[gene_symbol_column]]),
      count_key = count_features_character
    )
  }

  if (length(matching_strategies) == 0L) {
    stop(
      "Annotation metadata has neither an Ensembl-ID nor gene-symbol column.",
      call. = FALSE
    )
  }

  strategy_scores <- vapply(
    matching_strategies,
    function(strategy) {
      sum(
        strategy$count_key %in% strategy$annotation_key
      )
    },
    integer(1)
  )

  selected_strategy_name <- names(
    which.max(strategy_scores)
  )[[1L]]

  selected_strategy <- matching_strategies[[selected_strategy_name]]

  annotation_match_index <- match(
    selected_strategy$count_key,
    selected_strategy$annotation_key
  )

  keep_count_features <- !is.na(annotation_match_index)

  selected_features <- count_features_character[
    keep_count_features
  ]

  selected_annotation <- annotation[
    annotation_match_index[keep_count_features],
    ,
    drop = FALSE
  ]

  if (length(selected_features) < 1000L) {
    stop(
      "Only ",
      length(selected_features),
      " count features matched protein-coding autosomal/X annotation. ",
      "Selected matching strategy: ",
      selected_strategy_name,
      ". This is unexpectedly low.",
      call. = FALSE
    )
  }

  ensembl_values <- if (!is.null(ensembl_column)) {
    as.character(selected_annotation[[ensembl_column]])
  } else {
    rep(NA_character_, length(selected_features))
  }

  gene_values <- if (!is.null(gene_symbol_column)) {
    as.character(selected_annotation[[gene_symbol_column]])
  } else {
    selected_features
  }

  chromosome_values <- as.character(
    selected_annotation[[chromosome_column]]
  )

  gene_metadata <- data.frame(
    feature_id = selected_features,
    ensembl_gene_id = ensembl_values,
    gene = gene_values,
    chromosome = chromosome_values,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (anyDuplicated(gene_metadata$feature_id)) {
    stop(
      "Duplicated feature IDs were generated after annotation matching.",
      call. = FALSE
    )
  }

  list(
    genes = selected_features,
    gene_metadata = gene_metadata,
    matching_strategy = selected_strategy_name,
    strategy_scores = strategy_scores
  )
}

make_spot_to_group_matrix <- function(
  spot_metadata,
  sample_levels,
  cluster_levels
) {
  group_table <- expand.grid(
    sample_ID = sample_levels,
    cluster_id = cluster_levels,
    stringsAsFactors = FALSE
  )

  group_table$group_id <- paste(
    group_table$sample_ID,
    group_table$cluster_id,
    sep = "__CLUSTER__"
  )

  spot_group_id <- paste(
    spot_metadata$sample_ID,
    spot_metadata$cluster_id,
    sep = "__CLUSTER__"
  )

  group_index <- match(
    spot_group_id,
    group_table$group_id
  )

  if (anyNA(group_index)) {
    stop(
      "Some spots could not be mapped to sample x cluster groups.",
      call. = FALSE
    )
  }

  aggregation_matrix <- Matrix::sparseMatrix(
    i = seq_len(nrow(spot_metadata)),
    j = group_index,
    x = 1,
    dims = c(
      nrow(spot_metadata),
      nrow(group_table)
    ),
    dimnames = list(
      rownames(spot_metadata),
      group_table$group_id
    )
  )

  list(
    aggregation_matrix = aggregation_matrix,
    group_table = group_table
  )
}

validate_count_matrix <- function(
  x,
  displayed_name
) {
  if (any(!is.finite(x@x))) {
    stop(
      displayed_name,
      " contains non-finite values.",
      call. = FALSE
    )
  }

  if (any(x@x < 0)) {
    stop(
      displayed_name,
      " contains negative values.",
      call. = FALSE
    )
  }

  non_integer_distance <- abs(
    x@x - round(x@x)
  )

  if (
    length(non_integer_distance) > 0L &&
      max(non_integer_distance) > 1e-8
  ) {
    stop(
      displayed_name,
      " contains non-integer values.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

repair_crescendo_matching_sentinels <- function(
  corrected_cells_by_genes,
  raw_counts_chunk,
  crescendo_metadata,
  chunk_label,
  output_qc_dir
) {
  # Each chunk is only ~40k spots x 100 genes, so conversion to a standard
  # numeric matrix is modest in memory (~32 MB) and avoids class-specific
  # behaviour of Matrix sparse objects in comparisons/which().
  corrected_dense <- as.matrix(
    corrected_cells_by_genes
  )

  storage.mode(corrected_dense) <- "double"

  negative_positions <- which(
    corrected_dense < 0,
    arr.ind = TRUE
  )

  if (nrow(negative_positions) == 0L) {
    corrected_sparse <- Matrix::Matrix(
      corrected_dense,
      sparse = TRUE
    )

    return(
      list(
        corrected_cells_by_genes = corrected_sparse,
        n_repaired = 0L,
        n_genes_repaired = 0L
      )
    )
  }

  negative_values <- corrected_dense[
    negative_positions
  ]

  if (any(negative_values != -1)) {
    stop(
      "Crescendo returned negative values other than the internal ",
      "-1 Matching failure sentinel in ",
      chunk_label,
      ". Values observed: ",
      paste(
        sort(unique(negative_values)),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  failed_spot_indices <- negative_positions[, 1L]
  failed_gene_indices <- negative_positions[, 2L]

  fallback_raw_counts <- numeric(
    nrow(negative_positions)
  )

  for (k in seq_len(nrow(negative_positions))) {
    fallback_raw_counts[[k]] <- as.numeric(
      raw_counts_chunk[
        failed_gene_indices[[k]],
        failed_spot_indices[[k]]
      ]
    )
  }

  qc_table <- data.frame(
    chunk_label = chunk_label,
    spot_ID = rownames(corrected_dense)[
      failed_spot_indices
    ],
    gene = colnames(corrected_dense)[
      failed_gene_indices
    ],
    crescendo_returned_value = as.numeric(
      negative_values
    ),
    fallback_raw_count = fallback_raw_counts,
    sample_ID = crescendo_metadata$sample_ID[
      failed_spot_indices
    ],
    slide_ID = as.character(
      crescendo_metadata$slide_ID[
        failed_spot_indices
      ]
    ),
    original_cluster = crescendo_metadata$original_cluster[
      failed_spot_indices
    ],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  qc_file <- file.path(
    output_qc_dir,
    paste0(
      analysis_prefix,
      "_",
      chunk_label,
      "_crescendoMatching_minus1Sentinels.tsv.gz"
    )
  )

  write_tsv_gz(
    qc_table,
    qc_file
  )

  # Conservative fallback:
  # for entries where Crescendo Matching failed and returned -1,
  # retain the original raw count only for that exact spot x gene.
  corrected_dense[
    negative_positions
  ] <- fallback_raw_counts

  corrected_sparse <- Matrix::Matrix(
    corrected_dense,
    sparse = TRUE
  )

  message(
    "Crescendo Matching returned ",
    format(nrow(negative_positions), big.mark = ","),
    " '-1' sentinel value(s) across ",
    length(unique(failed_gene_indices)),
    " gene(s) in ",
    chunk_label,
    "."
  )
  message(
    "Conservative fallback applied: original raw count retained only for ",
    "those failed spot x gene entries."
  )
  message(
    "QC written to: ",
    qc_file
  )

  list(
    corrected_cells_by_genes = corrected_sparse,
    n_repaired = nrow(negative_positions),
    n_genes_repaired = length(unique(failed_gene_indices))
  )
}

# ==============================================================================
# 6. Validate input files
# ==============================================================================

stage_message(1, "VALIDATING INPUT FILES")

required_input_files <- c(
  input_rdata_file,
  metadata_file,
  annotation_rds_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0L) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_input_files, collapse = "\n"),
    call. = FALSE
  )
}

message("Input Seurat RData:")
message(input_rdata_file)
message("\nSample metadata:")
message(metadata_file)
message("\nProtein-coding gene annotation RDS:")
message(annotation_rds_file)
message("\nOutput:")
message(output_crescendo_dir)

# ==============================================================================
# 7. Load Seurat object and retain the same 16 samples
# ==============================================================================

stage_message(2, "LOADING SEURAT OBJECT")

seurat_input <- find_seurat_object(
  input_rdata_file = input_rdata_file,
  cluster_column = cluster_column
)

seurat_object_name <- seurat_input$object_name
seurat_object <- seurat_input$object

message("Selected Seurat object: ", seurat_object_name)
message("Spots before exclusions: ", format(ncol(seurat_object), big.mark = ","))

if (!(sample_column %in% colnames(seurat_object@meta.data))) {
  stop(
    "Seurat metadata is missing sample column '",
    sample_column,
    "'.",
    call. = FALSE
  )
}

if (!(cluster_column %in% colnames(seurat_object@meta.data))) {
  stop(
    "Seurat metadata is missing cluster column '",
    cluster_column,
    "'.",
    call. = FALSE
  )
}

seurat_sample_ids <- as.character(
  seurat_object@meta.data[[sample_column]]
)

keep_spots <- !seurat_sample_ids %in% excluded_samples

seurat_object <- subset(
  seurat_object,
  cells = colnames(seurat_object)[keep_spots]
)

retained_samples <- unique(
  as.character(
    seurat_object@meta.data[[sample_column]]
  )
)

if (length(retained_samples) != expected_number_of_samples) {
  stop(
    "Expected ",
    expected_number_of_samples,
    " retained samples but found ",
    length(retained_samples),
    ":\n",
    paste(sort(retained_samples), collapse = "\n"),
    call. = FALSE
  )
}

message("Spots after exclusions: ", format(ncol(seurat_object), big.mark = ","))
message("Retained samples: ", length(retained_samples))

# ==============================================================================
# 8. Load and standardize sample metadata
# ==============================================================================

stage_message(3, "PREPARING SAMPLE AND SPOT METADATA")

sample_metadata_raw <- utils::read.delim(
  metadata_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata_sample_column <- detect_column(
  sample_metadata_raw,
  candidates = c(
    "sample_ID",
    "Sample",
    "sample",
    "sample_id"
  ),
  displayed_name = "sample ID"
)

metadata_slide_column <- detect_column(
  sample_metadata_raw,
  candidates = c(
    "Slide ID",
    "slide_ID",
    "slide_id",
    "SlideID",
    "slide"
  ),
  displayed_name = "Slide ID"
)

metadata_group_column <- detect_column(
  sample_metadata_raw,
  candidates = c(
    "fmt_donor_group",
    "Donor group",
    "donor_group",
    "group"
  ),
  displayed_name = "donor group"
)

metadata_sex_column <- detect_column(
  sample_metadata_raw,
  candidates = c(
    "sex",
    "Sex"
  ),
  displayed_name = "sex"
)

sample_metadata <- data.frame(
  sample_ID = as.character(
    sample_metadata_raw[[metadata_sample_column]]
  ),
  slide_ID = as.character(
    sample_metadata_raw[[metadata_slide_column]]
  ),
  fmt_donor_group = as.character(
    sample_metadata_raw[[metadata_group_column]]
  ),
  sex = as.character(
    sample_metadata_raw[[metadata_sex_column]]
  ),
  stringsAsFactors = FALSE
)

sample_metadata <- sample_metadata[
  sample_metadata$sample_ID %in% retained_samples,
  ,
  drop = FALSE
]

sample_metadata <- unique(sample_metadata)

if (anyDuplicated(sample_metadata$sample_ID)) {
  duplicated_samples <- unique(
    sample_metadata$sample_ID[
      duplicated(sample_metadata$sample_ID)
    ]
  )
  stop(
    "Sample metadata has non-unique retained sample IDs: ",
    paste(duplicated_samples, collapse = ", "),
    call. = FALSE
  )
}

missing_metadata_samples <- setdiff(
  retained_samples,
  sample_metadata$sample_ID
)

if (length(missing_metadata_samples) > 0L) {
  stop(
    "Sample metadata is missing retained sample(s): ",
    paste(missing_metadata_samples, collapse = ", "),
    call. = FALSE
  )
}

# Preserve Seurat sample order.
sample_order <- unique(
  as.character(
    seurat_object@meta.data[[sample_column]]
  )
)

sample_metadata <- sample_metadata[
  match(
    sample_order,
    sample_metadata$sample_ID
  ),
  ,
  drop = FALSE
]

seurat_meta <- seurat_object@meta.data

spot_metadata <- data.frame(
  spot_id = rownames(seurat_meta),
  sample_ID = as.character(
    seurat_meta[[sample_column]]
  ),
  cluster_id = as.character(
    seurat_meta[[cluster_column]]
  ),
  stringsAsFactors = FALSE
)

sample_match <- match(
  spot_metadata$sample_ID,
  sample_metadata$sample_ID
)

spot_metadata$slide_ID <- sample_metadata$slide_ID[
  sample_match
]
spot_metadata$fmt_donor_group <- sample_metadata$fmt_donor_group[
  sample_match
]
spot_metadata$sex <- sample_metadata$sex[
  sample_match
]

if (anyNA(spot_metadata$slide_ID)) {
  stop(
    "Slide ID is missing for one or more retained spots.",
    call. = FALSE
  )
}

cluster_order <- natural_order(
  unique(spot_metadata$cluster_id)
)

if (length(cluster_order) != expected_number_of_clusters) {
  stop(
    "Expected ",
    expected_number_of_clusters,
    " clusters but found ",
    length(cluster_order),
    ": ",
    paste(cluster_order, collapse = ", "),
    call. = FALSE
  )
}

rownames(spot_metadata) <- spot_metadata$spot_id

sample_slide_summary <- unique(
  sample_metadata[
    ,
    c(
      "sample_ID",
      "slide_ID",
      "fmt_donor_group",
      "sex"
    ),
    drop = FALSE
  ]
)

write_tsv_gz(
  sample_slide_summary,
  file.path(
    output_metadata_dir,
    paste0(
      analysis_prefix,
      "_sampleSlideGroupSex_metadata.tsv.gz"
    )
  )
)

spot_count_summary <- as.data.frame(
  xtabs(
    ~ sample_ID + slide_ID + fmt_donor_group + sex + cluster_id,
    data = spot_metadata
  )
)

spot_count_summary <- spot_count_summary[
  spot_count_summary$Freq > 0,
  ,
  drop = FALSE
]

colnames(spot_count_summary)[
  colnames(spot_count_summary) == "Freq"
] <- "n_spots"

write_tsv_gz(
  spot_count_summary,
  file.path(
    output_metadata_dir,
    paste0(
      analysis_prefix,
      "_spotCounts_sampleSlideGroupSexCluster.tsv.gz"
    )
  )
)

message("Slides: ", paste(unique(sample_metadata$slide_ID), collapse = ", "))
message("Clusters: ", paste(cluster_order, collapse = ", "))

# ==============================================================================
# 9. Extract raw counts, retain protein-coding genes, then apply expression filter
# ==============================================================================

stage_message(4, "EXTRACTING RAW SPOT COUNTS AND FILTERING GENES")

raw_counts_all <- get_raw_count_matrix(
  seurat_object = seurat_object,
  assay_name = assay_name
)

if (!identical(
  colnames(raw_counts_all),
  rownames(spot_metadata)
)) {
  missing_count_spots <- setdiff(
    rownames(spot_metadata),
    colnames(raw_counts_all)
  )

  if (length(missing_count_spots) > 0L) {
    stop(
      "Raw count matrix is missing ",
      length(missing_count_spots),
      " retained spots.",
      call. = FALSE
    )
  }

  raw_counts_all <- raw_counts_all[
    ,
    rownames(spot_metadata),
    drop = FALSE
  ]
}

raw_nUMI <- Matrix::colSums(
  raw_counts_all
)

if (any(raw_nUMI <= 0)) {
  stop(
    "At least one retained spot has total raw UMI <= 0.",
    call. = FALSE
  )
}

number_of_raw_genes <- nrow(
  raw_counts_all
)

message(
  "Raw genes before protein-coding filter: ",
  format(
    number_of_raw_genes,
    big.mark = ","
  )
)

# --------------------------------------------------------------------------
# Step 1: keep only the already-prepared protein-coding gene set.
# --------------------------------------------------------------------------

protein_coding_match <- load_and_match_protein_coding_genes(
  annotation_rds_file = annotation_rds_file,
  count_features = rownames(raw_counts_all)
)

protein_coding_genes <- protein_coding_match$genes

raw_counts <- raw_counts_all[
  protein_coding_genes,
  ,
  drop = FALSE
]

rm(raw_counts_all)
invisible(gc())

number_of_protein_coding_genes <- nrow(
  raw_counts
)

write_tsv_gz(
  protein_coding_match$candidate_scores,
  file.path(
    output_metadata_dir,
    paste0(
      analysis_prefix,
      "_proteinCodingRDS_matchingCandidates.tsv.gz"
    )
  )
)

message(
  "Protein-coding genes before expression filter: ",
  format(
    number_of_protein_coding_genes,
    big.mark = ","
  )
)

# --------------------------------------------------------------------------
# Step 2: cluster-aware expression filter.
#
# Keep a protein-coding gene if count > 0 in >=10% of spots in at least one
# existing Leiden cluster. This removes extremely sparse genes that can make
# Poisson glmnet unstable while preserving genes specific to small regions.
# --------------------------------------------------------------------------

cluster_ids_for_filter <- as.character(
  spot_metadata$cluster_id
)

cluster_levels_for_filter <- natural_order(
  unique(cluster_ids_for_filter)
)

gene_total_counts <- Matrix::rowSums(
  raw_counts
)

max_detected_fraction_any_cluster <- rep(
  0,
  nrow(raw_counts)
)

names(max_detected_fraction_any_cluster) <- rownames(
  raw_counts
)

cluster_filter_summary_rows <- vector(
  mode = "list",
  length = length(cluster_levels_for_filter)
)

names(cluster_filter_summary_rows) <-
  cluster_levels_for_filter

for (cluster_filter_current in cluster_levels_for_filter) {
  cluster_spot_indices <- which(
    cluster_ids_for_filter == cluster_filter_current
  )

  if (length(cluster_spot_indices) == 0L) {
    next
  }

  detected_fraction_current <- Matrix::rowSums(
    raw_counts[
      ,
      cluster_spot_indices,
      drop = FALSE
    ] > 0
  ) / length(cluster_spot_indices)

  max_detected_fraction_any_cluster <- pmax(
    max_detected_fraction_any_cluster,
    as.numeric(detected_fraction_current)
  )

  cluster_filter_summary_rows[[cluster_filter_current]] <-
    data.frame(
      cluster_id = cluster_filter_current,
      n_spots = length(cluster_spot_indices),
      minimum_detected_spots_at_10pct = ceiling(
        length(cluster_spot_indices) *
          gene_min_detected_fraction_within_any_cluster
      ),
      stringsAsFactors = FALSE
    )
}

cluster_filter_summary <- do.call(
  rbind,
  cluster_filter_summary_rows
)

rownames(cluster_filter_summary) <- NULL

write_tsv_gz(
  cluster_filter_summary,
  file.path(
    output_metadata_dir,
    paste0(
      analysis_prefix,
      "_clusterAwareGeneFilter_clusterSummary.tsv.gz"
    )
  )
)

gene_filter_qc <- data.frame(
  feature_id = rownames(raw_counts),
  total_raw_counts = as.numeric(
    gene_total_counts
  ),
  max_detected_fraction_any_cluster = as.numeric(
    max_detected_fraction_any_cluster
  ),
  keep_for_crescendo = as.numeric(
    max_detected_fraction_any_cluster
  ) >= gene_min_detected_fraction_within_any_cluster,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write_tsv_gz(
  gene_filter_qc,
  file.path(
    output_metadata_dir,
    paste0(
      analysis_prefix,
      "_proteinCoding_geneFilter_clusterAware10pct_beforeCrescendo.tsv.gz"
    )
  )
)

genes_use <- gene_filter_qc$feature_id[
  gene_filter_qc$keep_for_crescendo
]

if (length(genes_use) == 0L) {
  stop(
    "No protein-coding genes passed the cluster-aware Crescendo expression filter.",
    call. = FALSE
  )
}

raw_counts <- raw_counts[
  genes_use,
  ,
  drop = FALSE
]

if (test_mode) {
  genes_use <- head(
    genes_use,
    test_number_of_genes
  )

  raw_counts <- raw_counts[
    genes_use,
    ,
    drop = FALSE
  ]
}

number_of_final_crescendo_genes <- length(
  genes_use
)

gene_metadata <- data.frame(
  feature_id = genes_use,
  gene = genes_use,
  total_raw_counts = as.numeric(
    gene_total_counts[genes_use]
  ),
  max_detected_fraction_any_cluster = as.numeric(
    max_detected_fraction_any_cluster[genes_use]
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_filter_stage_summary <- data.frame(
  stage = c(
    "raw_RNA_genes",
    "matched_protein_coding_genes",
    "protein_coding_genes_after_clusterAware10pct_filter"
  ),
  n_genes = c(
    number_of_raw_genes,
    number_of_protein_coding_genes,
    number_of_final_crescendo_genes
  ),
  stringsAsFactors = FALSE
)

write_tsv_gz(
  gene_filter_stage_summary,
  file.path(
    output_metadata_dir,
    paste0(
      analysis_prefix,
      "_geneFilter_stageSummary.tsv.gz"
    )
  )
)

write_tsv_gz(
  gene_metadata,
  file.path(
    output_metadata_dir,
    paste0(
      analysis_prefix,
      "_genesUsed_proteinCoding_clusterAware10pct.tsv.gz"
    )
  )
)

message(
  "Cluster-aware detection threshold: >= ",
  gene_min_detected_fraction_within_any_cluster * 100,
  "% of spots in at least one cluster"
)
message(
  "Final genes selected for Crescendo: ",
  format(
    number_of_final_crescendo_genes,
    big.mark = ","
  )
)
message(
  "Filtered raw count dimensions: ",
  nrow(raw_counts),
  " genes x ",
  ncol(raw_counts),
  " spots"
)

# ==============================================================================
# 10. Prepare Harmony input and obtain soft-cluster matrix R
# ==============================================================================

stage_message(5, "RUNNING HARMONY FOR CRESCENDO SOFT CLUSTERS")

selected_reduction <- select_harmony_reduction(
  seurat_object = seurat_object,
  preferred_reductions = harmony_preferred_reductions
)

harmony_embeddings <- extract_reduction_embeddings(
  seurat_object = seurat_object,
  reduction_name = selected_reduction,
  spot_ids = rownames(spot_metadata),
  n_dims = harmony_n_dims
)

crescendo_metadata <- data.frame(
  sample_ID = spot_metadata$sample_ID,
  slide_ID = factor(spot_metadata$slide_ID),
  fmt_donor_group = spot_metadata$fmt_donor_group,
  sex = spot_metadata$sex,
  original_cluster = spot_metadata$cluster_id,
  nUMI = as.numeric(
    raw_nUMI[
      rownames(spot_metadata)
    ]
  ),
  stringsAsFactors = FALSE,
  row.names = rownames(spot_metadata)
)

if (!identical(
  rownames(crescendo_metadata),
  rownames(harmony_embeddings)
)) {
  stop(
    "Harmony embeddings and Crescendo metadata are not in identical spot order.",
    call. = FALSE
  )
}

harmony_object_file <- file.path(
  output_harmony_dir,
  paste0(
    analysis_prefix,
    "_harmonyForCrescendo_v13_nclust",
    harmony_nclust,
    "_object.rds"
  )
)

harmony_r_file <- file.path(
  output_harmony_dir,
  paste0(
    analysis_prefix,
    "_harmonyForCrescendo_v13_nclust",
    harmony_nclust,
    "_R.rds"
  )
)

if (
  resume_existing_chunks &&
    file.exists(harmony_r_file)
) {
  message("Using existing Harmony R matrix:")
  message(harmony_r_file)

  harmony_R <- readRDS(
    harmony_r_file
  )

  if (
    ncol(harmony_R) != nrow(crescendo_metadata) ||
      nrow(harmony_R) != harmony_nclust
  ) {
    stop(
      "Existing Harmony R matrix dimensions are ",
      paste(dim(harmony_R), collapse = " x "),
      " but expected ",
      harmony_nclust,
      " soft clusters x ",
      nrow(crescendo_metadata),
      " spots. Delete the incompatible Harmony output and rerun.",
      call. = FALSE
    )
  }
} else {
  message("Input reduction: ", selected_reduction)
  message("Harmony dimensions: ", ncol(harmony_embeddings))
  message("Harmony soft clusters (nclust): ", harmony_nclust)
  message("Harmony batch variable: slide_ID")

  harmony_object <- harmony::RunHarmony(
    data_mat = harmony_embeddings,
    meta_data = crescendo_metadata,
    vars_use = "slide_ID",
    theta = harmony_theta,
    nclust = harmony_nclust,
    return_object = TRUE,
    verbose = harmony_verbose
  )

  harmony_R <- harmony_object$R

  if (is.null(harmony_R)) {
    stop(
      "Harmony object does not contain soft-cluster matrix R.",
      call. = FALSE
    )
  }

  saveRDS(
    harmony_object,
    harmony_object_file,
    compress = "xz"
  )

  saveRDS(
    harmony_R,
    harmony_r_file,
    compress = "xz"
  )

  rm(harmony_object)
  invisible(gc())
}

if (ncol(harmony_R) != nrow(crescendo_metadata)) {
  stop(
    "Harmony R must be cluster x spot with ",
    nrow(crescendo_metadata),
    " spots, but dimensions are ",
    paste(dim(harmony_R), collapse = " x "),
    ".",
    call. = FALSE
  )
}

message(
  "Harmony R dimensions: ",
  nrow(harmony_R),
  " soft clusters x ",
  ncol(harmony_R),
  " spots"
)

# ==============================================================================
# 11. Create sample x cluster aggregation matrix
# ==============================================================================

stage_message(6, "PREPARING PSEUDOBULK AGGREGATION MATRIX")

aggregation_info <- make_spot_to_group_matrix(
  spot_metadata = spot_metadata,
  sample_levels = sample_order,
  cluster_levels = cluster_order
)

aggregation_matrix <- aggregation_info$aggregation_matrix
group_table <- aggregation_info$group_table

group_table$n_spots <- as.numeric(
  Matrix::colSums(
    aggregation_matrix
  )
)

group_table <- merge(
  group_table,
  sample_metadata,
  by = "sample_ID",
  all.x = TRUE,
  sort = FALSE
)

group_table <- group_table[
  match(
    colnames(aggregation_matrix),
    group_table$group_id
  ),
  ,
  drop = FALSE
]

if (!identical(
  group_table$group_id,
  colnames(aggregation_matrix)
)) {
  stop(
    "Internal sample x cluster group order mismatch.",
    call. = FALSE
  )
}

write_tsv_gz(
  group_table,
  file.path(
    output_metadata_dir,
    paste0(
      analysis_prefix,
      "_pseudobulkGroupDefinitions.tsv.gz"
    )
  )
)

# ==============================================================================
# 12. Split genes into chunks
# ==============================================================================

stage_message(7, "PREPARING GENE CHUNKS")

gene_chunk_index <- ceiling(
  seq_along(genes_use) / gene_chunk_size
)

gene_chunks <- split(
  genes_use,
  gene_chunk_index
)

message("Gene chunk size: ", gene_chunk_size)
message("Number of chunks: ", length(gene_chunks))
message("Crescendo cores: ", crescendo_mc_cores)

corrected_pseudobulk_chunks <- vector(
  mode = "list",
  length = length(gene_chunks)
)

raw_pseudobulk_chunks <- vector(
  mode = "list",
  length = length(gene_chunks)
)

corrected_spot_library_sizes <- setNames(
  numeric(nrow(crescendo_metadata)),
  rownames(crescendo_metadata)
)

chunk_summary_rows <- vector(
  mode = "list",
  length = length(gene_chunks)
)

# ==============================================================================
# 13. Run Crescendo chunk-by-chunk and pseudobulk each corrected chunk
# ==============================================================================

stage_message(8, "RUNNING CRESCENDO BATCH CORRECTION")

analysis_start_time <- Sys.time()

for (chunk_number in seq_along(gene_chunks)) {
  chunk_start_time <- Sys.time()
  n_matching_sentinels_repaired <- 0L
  n_genes_with_matching_sentinels <- 0L
  genes_chunk <- gene_chunks[[chunk_number]]

  chunk_label <- sprintf(
    "chunk_%04d_genes_%05d-%05d",
    chunk_number,
    min(which(gene_chunk_index == chunk_number)),
    max(which(gene_chunk_index == chunk_number))
  )

  corrected_spot_chunk_file <- file.path(
    output_corrected_spot_chunks_dir,
    paste0(
      analysis_prefix,
      "_",
      chunk_label,
      "_crescendoCorrected_spotCounts.rds"
    )
  )

  pseudobulk_chunk_file <- file.path(
    output_pseudobulk_chunks_dir,
    paste0(
      analysis_prefix,
      "_",
      chunk_label,
      "_rawAndCorrected_pseudobulk.rds"
    )
  )

  message("\n------------------------------------------------------------")
  message(
    "[",
    chunk_number,
    "/",
    length(gene_chunks),
    "] ",
    chunk_label,
    " | genes: ",
    length(genes_chunk)
  )
  message("------------------------------------------------------------")

  if (
    resume_existing_chunks &&
      file.exists(corrected_spot_chunk_file) &&
      file.exists(pseudobulk_chunk_file)
  ) {
    message("Reusing completed Crescendo chunk.")

    corrected_chunk_saved <- readRDS(
      corrected_spot_chunk_file
    )

    pseudobulk_chunk_saved <- readRDS(
      pseudobulk_chunk_file
    )

    corrected_spot_counts_chunk <- corrected_chunk_saved$corrected_counts_gene_by_spot

    if (!identical(
      rownames(corrected_spot_counts_chunk),
      genes_chunk
    )) {
      stop(
        "Existing corrected chunk has unexpected genes: ",
        corrected_spot_chunk_file,
        call. = FALSE
      )
    }

    corrected_pseudobulk_chunks[[chunk_number]] <-
      pseudobulk_chunk_saved$corrected_pseudobulk

    raw_pseudobulk_chunks[[chunk_number]] <-
      pseudobulk_chunk_saved$raw_pseudobulk
  } else {
    set.seed(
      crescendo_seed + chunk_number
    )

    raw_counts_chunk <- raw_counts[
      genes_chunk,
      ,
      drop = FALSE
    ]

    message("Running Crescendo...")

    corrected_cells_by_genes <- crescendo::crescendo(
      Ycounts = raw_counts_chunk,
      meta_data = crescendo_metadata,
      R = harmony_R,
      batch_var = "slide_ID",
      genes_use = NULL,
      prop = crescendo_prop,
      min_cells = crescendo_min_cells,
      lambda = crescendo_lambda,
      alpha = crescendo_alpha,
      seed = crescendo_seed + chunk_number,
      mc.cores = crescendo_mc_cores,
      return_obj = FALSE,
      verbose = crescendo_verbose
    )

    corrected_cells_by_genes <- methods::as(
      corrected_cells_by_genes,
      "dgCMatrix"
    )

    if (!identical(
      rownames(corrected_cells_by_genes),
      rownames(crescendo_metadata)
    )) {
      stop(
        "Crescendo corrected matrix has unexpected spot order.",
        call. = FALSE
      )
    }

    if (!identical(
      colnames(corrected_cells_by_genes),
      genes_chunk
    )) {
      stop(
        "Crescendo corrected matrix has unexpected gene order.",
        call. = FALSE
      )
    }

    matching_repair <- repair_crescendo_matching_sentinels(
      corrected_cells_by_genes = corrected_cells_by_genes,
      raw_counts_chunk = raw_counts_chunk,
      crescendo_metadata = crescendo_metadata,
      chunk_label = chunk_label,
      output_qc_dir = output_qc_dir
    )

    corrected_cells_by_genes <-
      matching_repair$corrected_cells_by_genes

    n_matching_sentinels_repaired <-
      matching_repair$n_repaired

    n_genes_with_matching_sentinels <-
      matching_repair$n_genes_repaired

    rm(matching_repair)

    validate_count_matrix(
      corrected_cells_by_genes,
      displayed_name = paste0(
        "Crescendo corrected matrix ",
        chunk_label
      )
    )

    corrected_spot_counts_chunk <- Matrix::t(
      corrected_cells_by_genes
    )

    corrected_spot_counts_chunk <- Matrix::drop0(
      methods::as(
        corrected_spot_counts_chunk,
        "dgCMatrix"
      )
    )

    raw_pseudobulk_chunk <- raw_counts_chunk %*%
      aggregation_matrix

    corrected_pseudobulk_chunk <- corrected_spot_counts_chunk %*%
      aggregation_matrix

    raw_pseudobulk_chunk <- Matrix::drop0(
      methods::as(
        raw_pseudobulk_chunk,
        "dgCMatrix"
      )
    )

    corrected_pseudobulk_chunk <- Matrix::drop0(
      methods::as(
        corrected_pseudobulk_chunk,
        "dgCMatrix"
      )
    )

    validate_count_matrix(
      corrected_pseudobulk_chunk,
      displayed_name = paste0(
        "Crescendo corrected pseudobulk ",
        chunk_label
      )
    )

    saveRDS(
      list(
        chunk_number = chunk_number,
        chunk_label = chunk_label,
        genes = genes_chunk,
        n_matching_sentinels_repaired =
          n_matching_sentinels_repaired,
        n_genes_with_matching_sentinels =
          n_genes_with_matching_sentinels,
        corrected_counts_gene_by_spot = corrected_spot_counts_chunk
      ),
      corrected_spot_chunk_file,
      compress = "xz"
    )

    saveRDS(
      list(
        chunk_number = chunk_number,
        chunk_label = chunk_label,
        genes = genes_chunk,
        group_table = group_table,
        raw_pseudobulk = raw_pseudobulk_chunk,
        corrected_pseudobulk = corrected_pseudobulk_chunk
      ),
      pseudobulk_chunk_file,
      compress = "xz"
    )

    corrected_pseudobulk_chunks[[chunk_number]] <-
      corrected_pseudobulk_chunk

    raw_pseudobulk_chunks[[chunk_number]] <-
      raw_pseudobulk_chunk

    rm(
      raw_counts_chunk,
      corrected_cells_by_genes,
      raw_pseudobulk_chunk,
      corrected_pseudobulk_chunk
    )

    invisible(gc())
  }

  corrected_spot_library_sizes <- corrected_spot_library_sizes +
    as.numeric(
      Matrix::colSums(
        corrected_spot_counts_chunk
      )
    )

  chunk_summary_rows[[chunk_number]] <- data.frame(
    chunk_number = chunk_number,
    chunk_label = chunk_label,
    n_genes = length(genes_chunk),
    n_matching_sentinels_repaired =
      n_matching_sentinels_repaired,
    n_genes_with_matching_sentinels =
      n_genes_with_matching_sentinels,
    first_gene = genes_chunk[[1L]],
    last_gene = genes_chunk[[length(genes_chunk)]],
    corrected_spot_chunk_file = basename(
      corrected_spot_chunk_file
    ),
    pseudobulk_chunk_file = basename(
      pseudobulk_chunk_file
    ),
    stringsAsFactors = FALSE
  )

  rm(corrected_spot_counts_chunk)
  invisible(gc())

  chunk_elapsed_minutes <- as.numeric(
    difftime(
      Sys.time(),
      chunk_start_time,
      units = "mins"
    )
  )

  total_elapsed_minutes <- as.numeric(
    difftime(
      Sys.time(),
      analysis_start_time,
      units = "mins"
    )
  )

  mean_minutes_per_chunk <- total_elapsed_minutes / chunk_number
  remaining_chunks <- length(gene_chunks) - chunk_number
  estimated_remaining_hours <-
    (mean_minutes_per_chunk * remaining_chunks) / 60

  message(
    "Chunk time: ",
    sprintf("%.2f", chunk_elapsed_minutes),
    " min | mean: ",
    sprintf("%.2f", mean_minutes_per_chunk),
    " min/chunk | estimated remaining: ",
    sprintf("%.2f", estimated_remaining_hours),
    " h"
  )
}

# ==============================================================================
# 14. Combine pseudobulk chunks
# ==============================================================================

stage_message(9, "COMBINING PSEUDOBULK CHUNKS")

raw_pseudobulk_all <- do.call(
  rbind,
  raw_pseudobulk_chunks
)

corrected_pseudobulk_all <- do.call(
  rbind,
  corrected_pseudobulk_chunks
)

raw_pseudobulk_all <- methods::as(
  raw_pseudobulk_all,
  "dgCMatrix"
)

corrected_pseudobulk_all <- methods::as(
  corrected_pseudobulk_all,
  "dgCMatrix"
)

if (!identical(
  rownames(raw_pseudobulk_all),
  genes_use
)) {
  stop(
    "Combined raw pseudobulk gene order is incorrect.",
    call. = FALSE
  )
}

if (!identical(
  rownames(corrected_pseudobulk_all),
  genes_use
)) {
  stop(
    "Combined corrected pseudobulk gene order is incorrect.",
    call. = FALSE
  )
}

if (!identical(
  colnames(corrected_pseudobulk_all),
  group_table$group_id
)) {
  stop(
    "Combined corrected pseudobulk group order is incorrect.",
    call. = FALSE
  )
}

validate_count_matrix(
  corrected_pseudobulk_all,
  displayed_name = "Combined Crescendo corrected pseudobulk"
)

message(
  "Combined corrected pseudobulk dimensions: ",
  nrow(corrected_pseudobulk_all),
  " genes x ",
  ncol(corrected_pseudobulk_all),
  " sample-cluster groups"
)

# ==============================================================================
# 15. Save per-cluster pseudobulk matrices ready for edgeR
# ==============================================================================

stage_message(10, "WRITING PER-CLUSTER COUNT MATRICES")

corrected_counts_by_cluster <- vector(
  mode = "list",
  length = length(cluster_order)
)

raw_counts_by_cluster <- vector(
  mode = "list",
  length = length(cluster_order)
)

names(corrected_counts_by_cluster) <- cluster_order
names(raw_counts_by_cluster) <- cluster_order

for (cluster_id_current in cluster_order) {
  cluster_group_indices <- which(
    group_table$cluster_id == cluster_id_current
  )

  cluster_group_table <- group_table[
    cluster_group_indices,
    ,
    drop = FALSE
  ]

  sample_indices <- match(
    sample_order,
    cluster_group_table$sample_ID
  )

  if (anyNA(sample_indices)) {
    stop(
      "Cluster ",
      cluster_id_current,
      " is missing one or more sample definitions.",
      call. = FALSE
    )
  }

  cluster_group_indices <- cluster_group_indices[
    sample_indices
  ]

  raw_cluster_matrix <- raw_pseudobulk_all[
    ,
    cluster_group_indices,
    drop = FALSE
  ]

  corrected_cluster_matrix <- corrected_pseudobulk_all[
    ,
    cluster_group_indices,
    drop = FALSE
  ]

  colnames(raw_cluster_matrix) <- sample_order
  colnames(corrected_cluster_matrix) <- sample_order

  raw_counts_by_cluster[[cluster_id_current]] <-
    raw_cluster_matrix

  corrected_counts_by_cluster[[cluster_id_current]] <-
    corrected_cluster_matrix

  gene_metadata_cluster <- gene_metadata[
    match(
      rownames(corrected_cluster_matrix),
      gene_metadata$feature_id
    ),
    ,
    drop = FALSE
  ]

  output_corrected_df <- data.frame(
    gene_metadata_cluster[
      ,
      c(
        "feature_id",
        "gene"
      ),
      drop = FALSE
    ],
    as.matrix(corrected_cluster_matrix),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  output_raw_df <- data.frame(
    gene_metadata_cluster[
      ,
      c(
        "feature_id",
        "gene"
      ),
      drop = FALSE
    ],
    as.matrix(raw_cluster_matrix),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  cluster_slug <- sanitize_file_component(
    cluster_id_current
  )

  corrected_cluster_file <- file.path(
    output_pseudobulk_per_cluster_dir,
    paste0(
      analysis_prefix,
      "_cluster_",
      cluster_slug,
      "_crescendoCorrectedPseudobulkCounts_",
      "proteinCoding_clusterAware10pct.tsv.gz"
    )
  )

  raw_cluster_file <- file.path(
    output_pseudobulk_per_cluster_dir,
    paste0(
      analysis_prefix,
      "_cluster_",
      cluster_slug,
      "_rawPseudobulkCounts_sameSpots_",
      "proteinCoding_clusterAware10pct.tsv.gz"
    )
  )

  write_tsv_gz(
    output_corrected_df,
    corrected_cluster_file
  )

  write_tsv_gz(
    output_raw_df,
    raw_cluster_file
  )

  message(
    "Saved cluster ",
    cluster_id_current,
    ": ",
    basename(corrected_cluster_file)
  )

  rm(
    raw_cluster_matrix,
    corrected_cluster_matrix,
    output_corrected_df,
    output_raw_df
  )

  invisible(gc())
}

# ==============================================================================
# 16. QC summaries
# ==============================================================================

stage_message(11, "WRITING QC SUMMARIES")

chunk_summary <- do.call(
  rbind,
  chunk_summary_rows
)

rownames(chunk_summary) <- NULL

write_tsv_gz(
  chunk_summary,
  file.path(
    output_qc_dir,
    paste0(
      analysis_prefix,
      "_crescendoChunkSummary.tsv.gz"
    )
  )
)

spot_library_size_qc <- data.frame(
  spot_id = rownames(crescendo_metadata),
  sample_ID = crescendo_metadata$sample_ID,
  slide_ID = as.character(crescendo_metadata$slide_ID),
  fmt_donor_group = crescendo_metadata$fmt_donor_group,
  sex = crescendo_metadata$sex,
  cluster_id = crescendo_metadata$original_cluster,
  raw_total_UMI_allGenes = crescendo_metadata$nUMI,
  corrected_total_counts_selectedGenes = corrected_spot_library_sizes,
  stringsAsFactors = FALSE
)

write_tsv_gz(
  spot_library_size_qc,
  file.path(
    output_qc_dir,
    paste0(
      analysis_prefix,
      "_spotLibrarySize_beforeAfter.tsv.gz"
    )
  )
)

raw_group_library_sizes <- Matrix::colSums(
  raw_pseudobulk_all
)

corrected_group_library_sizes <- Matrix::colSums(
  corrected_pseudobulk_all
)

pseudobulk_library_size_qc <- group_table

pseudobulk_library_size_qc$raw_total_counts_selectedGenes <-
  as.numeric(
    raw_group_library_sizes[
      pseudobulk_library_size_qc$group_id
    ]
  )

pseudobulk_library_size_qc$crescendo_total_counts_selectedGenes <-
  as.numeric(
    corrected_group_library_sizes[
      pseudobulk_library_size_qc$group_id
    ]
  )

write_tsv_gz(
  pseudobulk_library_size_qc,
  file.path(
    output_qc_dir,
    paste0(
      analysis_prefix,
      "_pseudobulkLibrarySize_beforeAfter.tsv.gz"
    )
  )
)

slide_group_table <- as.data.frame(
  xtabs(
    ~ slide_ID + fmt_donor_group + sex,
    data = sample_metadata
  )
)

colnames(slide_group_table)[
  colnames(slide_group_table) == "Freq"
] <- "n_samples"

slide_group_table <- slide_group_table[
  slide_group_table$n_samples > 0,
  ,
  drop = FALSE
]

write_tsv_gz(
  slide_group_table,
  file.path(
    output_qc_dir,
    paste0(
      analysis_prefix,
      "_slideGroupSexConfoundingSummary.tsv.gz"
    )
  )
)

# ==============================================================================
# 17. Save compact R object for the next edgeR stage
# ==============================================================================

stage_message(12, "SAVING OBJECTS FOR THE NEXT edgeR STAGE")

settings <- list(
  project_root = project_root,
  dataset_name = dataset_name,
  clustering_name = clustering_name,
  analysis_name = analysis_name,
  input_rdata_file = input_rdata_file,
  metadata_file = metadata_file,
  annotation_rds_file = annotation_rds_file,
  seurat_object_name = seurat_object_name,
  assay_name = assay_name,
  sample_column = sample_column,
  cluster_column = cluster_column,
  excluded_samples = excluded_samples,
  sample_order = sample_order,
  cluster_order = cluster_order,
  selected_reduction_for_harmony = selected_reduction,
  harmony_n_dims = ncol(harmony_embeddings),
  harmony_nclust = harmony_nclust,
  harmony_theta = harmony_theta,
  crescendo_batch_var = "slide_ID",
  crescendo_prop = crescendo_prop,
  crescendo_min_cells = crescendo_min_cells,
  crescendo_lambda = crescendo_lambda,
  crescendo_alpha = crescendo_alpha,
  crescendo_seed = crescendo_seed,
  crescendo_mc_cores = crescendo_mc_cores,
  protein_coding_matching_source =
    protein_coding_match$best_candidate,
  protein_coding_matching_mode =
    protein_coding_match$matching_mode,
  gene_min_detected_fraction_within_any_cluster =
    gene_min_detected_fraction_within_any_cluster,
  gene_chunk_size = gene_chunk_size,
  test_mode = test_mode
)

edgeR_input_object <- list(
  settings = settings,
  sample_metadata = sample_metadata,
  group_table = group_table,
  gene_metadata = gene_metadata,
  raw_counts_by_cluster = raw_counts_by_cluster,
  crescendo_corrected_counts_by_cluster =
    corrected_counts_by_cluster
)

edgeR_input_object_file <- file.path(
  output_object_dir,
  paste0(
    analysis_prefix,
    "_pseudobulkCounts_beforeEdgeR.rds"
  )
)

saveRDS(
  edgeR_input_object,
  edgeR_input_object_file,
  compress = "xz"
)

settings_file <- file.path(
  output_object_dir,
  paste0(
    analysis_prefix,
    "_settings.rds"
  )
)

saveRDS(
  settings,
  settings_file,
  compress = "xz"
)

# ==============================================================================
# 18. Definitions / methods record
# ==============================================================================

definitions_file <- file.path(
  output_crescendo_dir,
  paste0(
    analysis_prefix,
    "_definitions.txt"
  )
)

writeLines(
  c(
    "Crescendo batch-correction preprocessing before edgeR.",
    "",
    paste0("Input Seurat RData: ", input_rdata_file),
    paste0("Input sample metadata: ", metadata_file),
    paste0("RNA assay: ", assay_name),
    paste0("Sample column: ", sample_column),
    paste0("Existing cluster column: ", cluster_column),
    paste0(
      "Excluded samples: ",
      paste(excluded_samples, collapse = ", ")
    ),
    paste0(
      "Retained biological samples: ",
      length(sample_order)
    ),
    paste0(
      "Retained existing clusters: ",
      length(cluster_order)
    ),
    "",
    "Crescendo input:",
    "  raw spot-level RNA counts",
    "  batch variable = Slide ID",
    paste0(
      "  Harmony soft clusters from reduction = ",
      selected_reduction
    ),
    paste0(
      "  Harmony dimensions = ",
      ncol(harmony_embeddings)
    ),
    paste0(
      "  Harmony soft clusters (nclust) = ",
      harmony_nclust
    ),
    "",
    "Gene scope:",
    paste0(
      "  protein-coding genes loaded from: ",
      annotation_rds_file
    ),
    paste0(
      "  protein-coding matching source = ",
      protein_coding_match$best_candidate
    ),
    paste0(
      "  protein-coding matching mode = ",
      protein_coding_match$matching_mode
    ),
    paste0(
      "  retained only if detected (>0 count) in at least ",
      gene_min_detected_fraction_within_any_cluster * 100,
      "% of spots in at least one existing Leiden cluster"
    ),
    paste0(
      "  number of corrected genes = ",
      length(genes_use)
    ),
    "",
    "Crescendo settings:",
    paste0("  prop = ", crescendo_prop),
    paste0("  min_cells = ", crescendo_min_cells),
    paste0("  alpha = ", crescendo_alpha),
    paste0("  gene chunk size = ", gene_chunk_size),
    paste0("  mc.cores = ", crescendo_mc_cores),
    "",
    "Post-Crescendo aggregation:",
    "  corrected spot counts are summed by sample x existing Leiden cluster",
    "  no re-clustering is performed",
    "  existing cluster labels are preserved",
    "",
    "Output for downstream differential expression:",
    paste0("  ", edgeR_input_object_file),
    "",
    "Interpretation note:",
    "  Slide ID is partially confounded with donor group and sex.",
    "  Therefore this Crescendo workflow is intended as a conservative",
    "  slide-effect sensitivity/rescue analysis rather than proof that all",
    "  corrected-away signal was purely technical."
  ),
  con = definitions_file
)

# ==============================================================================
# 19. Session information
# ==============================================================================

session_info_file <- file.path(
  output_crescendo_dir,
  paste0(
    analysis_prefix,
    "_sessionInfo.txt"
  )
)

capture.output(
  sessionInfo(),
  file = session_info_file
)

# ==============================================================================
# 20. Final summary
# ==============================================================================

stage_message(13, "CRESCENDO PREPROCESSING FINISHED")

message("Output directory:")
message(output_crescendo_dir)

message("\nCorrected pseudobulk object ready for edgeR:")
message(edgeR_input_object_file)

message("\nPer-cluster corrected count tables:")
message(output_pseudobulk_per_cluster_dir)

message("\nGenes corrected: ", format(length(genes_use), big.mark = ","))
message("Spots corrected: ", format(nrow(crescendo_metadata), big.mark = ","))
message("Samples: ", length(sample_order))
message("Clusters: ", length(cluster_order))
message("Slides: ", length(unique(sample_metadata$slide_ID)))

message("\nIMPORTANT:")
message(
  "Use 'crescendo_corrected_counts_by_cluster' from the saved RDS ",
  "as count input for the next edgeR stage."
)
