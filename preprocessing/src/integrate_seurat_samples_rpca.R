# ==============================================================================
# Merge and integrate Spatial Transcriptomics samples with Seurat v5 RPCA
#
# Integration is used only for UMAP and clustering.
# Raw RNA counts remain available for pseudobulk and differential expression.
# ==============================================================================

integrate_seurat_samples_rpca <- function(
    samples_list,
    nfeatures = 2000,
    dims = 1:30,
    project = "MaternalFMT",
    seed = 7,
    verbose = TRUE
) {

  if (!is.list(samples_list) || length(samples_list) < 2) {
    stop("`samples_list` must be a named list containing at least 2 Seurat objects.")
  }

  sample_ids <- names(samples_list)

  if (is.null(sample_ids) || any(sample_ids == "")) {
    stop("`samples_list` must have sample IDs as names.")
  }

  if (anyDuplicated(sample_ids) > 0) {
    stop("Duplicated sample IDs in `samples_list`.")
  }

  if (!all(vapply(samples_list, inherits, logical(1), what = "Seurat"))) {
    stop("Every element of `samples_list` must be a Seurat object.")
  }

  if (verbose) {
    message("Merging ", length(samples_list), " samples...")
  }

  integrated_object <- merge(
    x = samples_list[[1]],
    y = samples_list[-1],
    add.cell.ids = sample_ids,
    project = project,
    merge.data = FALSE,
    collapse = FALSE
  )

  SeuratObject::DefaultAssay(integrated_object) <- "RNA"

  if (!"sample_ID" %in% colnames(integrated_object[[]])) {
    stop("Merged object does not contain `sample_ID` metadata.")
  }

  count_layers <- SeuratObject::Layers(
    integrated_object[["RNA"]],
    search = "^counts"
  )

  # Safety fallback: normally merge() in Seurat v5 already creates one count
  # layer per sample. Split only if there is a single count layer.
  if (length(count_layers) == 1) {

    if (verbose) {
      message("Splitting RNA assay into sample-specific layers...")
    }

    integrated_object[["RNA"]] <- split(
      integrated_object[["RNA"]],
      f = integrated_object$sample_ID
    )
  }

  count_layers <- SeuratObject::Layers(
    integrated_object[["RNA"]],
    search = "^counts"
  )

  if (length(count_layers) < 2) {
    stop("Could not create sample-specific RNA count layers.")
  }

  if (verbose) {
    message("Normalizing RNA layers...")
  }

  integrated_object <- Seurat::NormalizeData(
    object = integrated_object,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = verbose
  )

  if (verbose) {
    message("Finding ", nfeatures, " variable features...")
  }

  integrated_object <- Seurat::FindVariableFeatures(
    object = integrated_object,
    assay = "RNA",
    selection.method = "vst",
    nfeatures = nfeatures,
    verbose = verbose
  )

  integration_features <- SeuratObject::VariableFeatures(
    integrated_object
  )

  if (length(integration_features) == 0) {
    stop("No variable features were identified.")
  }

  if (verbose) {
    message("Scaling data and running PCA...")
  }

  set.seed(seed)

  integrated_object <- Seurat::ScaleData(
    object = integrated_object,
    assay = "RNA",
    features = integration_features,
    verbose = verbose
  )

  integrated_object <- Seurat::RunPCA(
    object = integrated_object,
    assay = "RNA",
    features = integration_features,
    npcs = max(dims),
    reduction.name = "pca",
    verbose = verbose,
    seed.use = seed
  )

  if (verbose) {
    message("Running Seurat v5 RPCA integration...")
  }

  integrated_object <- Seurat::IntegrateLayers(
    object = integrated_object,
    method = Seurat::RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    dims = dims,
    verbose = verbose
  )

  # Integration is complete; join layers back into one RNA counts/data layer.
  integrated_object[["RNA"]] <- SeuratObject::JoinLayers(
    integrated_object[["RNA"]]
  )

  if (verbose) {
    message("Building integrated neighbor graph and UMAP...")
  }

  integrated_object <- Seurat::FindNeighbors(
    object = integrated_object,
    reduction = "integrated.rpca",
    dims = dims,
    verbose = verbose
  )

  integrated_object <- Seurat::RunUMAP(
    object = integrated_object,
    reduction = "integrated.rpca",
    dims = dims,
    reduction.name = "umap.integrated",
    reduction.key = "UMAPINT_",
    seed.use = seed,
    verbose = verbose
  )

  return(integrated_object)
}
