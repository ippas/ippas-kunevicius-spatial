library(dplyr)
library(tidyr)
library(magrittr)

load(
  "/home/mateusz/projects/ippas-kunevicius-spatial/results/"
  |> paste0(
    "maternalFMT_n16samples/",
    "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction/",
    "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040/",
    "01_statistics/03_edgeRResults/",
    "maternalFMT_n16samples_",
    "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040_",
    "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction_",
    "edgeRResults.RData"
  )
)

group_results <- edgeR_perCluster_combinedResults[["Overall_Group_ASD_vs_Neurotypical"]]

edgeR_perCluster_combinedResults[["Overall_Group_ASD_vs_Neurotypical"]] %>% 
  filter(grepl("Homer1", gene)) 

edgeR_perCluster_combinedResults[["Overall_Group_ASD_vs_Neurotypical"]] %>% 
  filter(grepl("Shank", gene)) 




edgeR_perCluster_combinedResults$ASD_Male_vs_Neurotypical_Male %>% 
  filter(gene == "Ndufb1") %>% filter(cluster_id == "16")

edgeR_perCluster_combinedResults$ASD_Female_vs_Neurotypical_Female %>% 
  filter(gene == "Ndufb1") %>% filter(cluster_id == "16")


input_seurat_rdata_file <- "/home/mateusz/projects/ippas-kunevicius-spatial/results/maternalFMT_n16samples/seurat_clustering_analysis/logNormalize_vst/hvg2000/02_cca/dims20/k20/prune0067/RData/02_maternalFMT_n16samples_logNormalizeVst_hvg2000_ccaDims20K20PruneSNN0067_res010to100by010_multiClusteringAndUmap.RData"
annotation_rds_file <- "/home/mateusz/projects/ippas-kunevicius-spatial/data/ensembl115/ensembl115_mouse_proteinCodingGenes.rds"

load(input_seurat_rdata_file)

seurat_object <- get(ls()[sapply(ls(), function(x) inherits(get(x), "Seurat"))][1])
protein_coding <- readRDS(annotation_rds_file)

seurat_genes <- rownames(seurat_object[["RNA"]])

protein_coding_genes <- if ("ensembl_gene_id" %in% colnames(protein_coding)) {
  protein_coding$ensembl_gene_id
} else if ("gene_id" %in% colnames(protein_coding)) {
  protein_coding$gene_id
} else {
  rownames(protein_coding)
}

cat("Geny w danych RNA:", length(unique(seurat_genes)), "\n")
cat("Geny protein-coding BioMart:", length(unique(protein_coding_genes)), "\n")
cat("Protein-coding obecne w moich danych:", length(intersect(seurat_genes, protein_coding_genes)), "\n")
# 11657
# group_results <- edgeR_perCluster_combinedResults[["Overall_Sex_Female_vs_Male"]]

group_significant <- group_results |>
  # filter(
  #   FDR < 0.05,
  #   abs(logFC) > 0.5
  # ) |>
  filter(
    FDR < 0.05
  ) |>
  mutate(
    direction = if_else(
      logFC > 0,
      "up",
      "down"
    )
  )

group_significant  %>% 
  # filter(direction == "up") %>% 
  group_by(cluster_id)  %>% 
  nest

group_significant  %>% 
  # filter(regulation == "down") %>% 
  group_by(cluster_id)  %>% 
  nest


crossing(
  FDR_threshold = c(0.05, 0.01),
  logFC_threshold = c(0, 0.5, 0.8, 1)
) %>%
  pmap_dfr(function(FDR_threshold, logFC_threshold) {

    tmp <- group_results %>%
      filter(
        FDR < FDR_threshold,
        abs(logFC) > logFC_threshold
      )

    tibble(
      FDR_threshold = FDR_threshold,
      logFC_threshold = logFC_threshold,
      n_results = nrow(tmp),
      n_results_up = sum(tmp$logFC > 0),
      n_results_down = sum(tmp$logFC < 0),
      n_genes_unique = n_distinct(tmp$gene)
    )
  }) %>%
  select(-c(n_results_up, n_results_down)) %>%
  arrange(match(FDR_threshold, c(0.05, 0.01)), logFC_threshold)

group_significant <- group_results |>
  # filter(
  #   FDR < 0.05,
  #   abs(logFC) > 0.5
  # ) |>
  mutate(
    direction = if_else(
      logFC > 0,
      "UP",
      "DOWN"
    )
  )
group_significant %>% 
    filter(
    FDR < 0.05,
    abs(logFC) > 0.5
  )  %>% 
  pull(gene) %>% as.character() %>% cat(sep = "\n")

group_significant$gene %>% 
  table %>% as.data.frame() %>% arrange(desc(Freq)) %>% head(20)
#!/usr/bin/env Rscript

# ==============================================================================
# 02_summarizeInteraction_FDR005_perCluster.R
#
# Simple summary of the edgeR interaction results per spatial cluster.
# Significance threshold: FDR < 0.05.
# ============================================================================== 

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"
clustering_name <-
  "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"

fdr_threshold <- 0.05

analysis_prefix <- paste0(
  dataset_name,
  "_",
  clustering_name,
  "_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction"
)

statistics_dir <- file.path(
  project_root,
  "results",
  dataset_name,
  "pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction",
  clustering_name,
  "01_statistics"
)

result_rdata_file <- file.path(
  statistics_dir,
  "03_edgeRResults",
  paste0(
    analysis_prefix,
    "_edgeRResults.RData"
  )
)


load(result_rdata_file)

if (!exists("edgeR_perCluster_combinedResults")) {
  stop(
    "Object `edgeR_perCluster_combinedResults` was not found in: ",
    result_rdata_file,
    call. = FALSE
  )
}

if (!"Interaction" %in% names(edgeR_perCluster_combinedResults)) {
  stop(
    "The `Interaction` result table was not found in ",
    "`edgeR_perCluster_combinedResults`.",
    call. = FALSE
  )
}

interaction_results <-
  edgeR_perCluster_combinedResults[["Interaction"]]

required_columns <- c(
  "cluster_id",
  "ensembl_gene_id",
  "gene",
  "logFC",
  "PValue",
  "FDR"
)

missing_columns <- setdiff(
  required_columns,
  colnames(interaction_results)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

interaction_results$cluster_id <-
  as.character(interaction_results$cluster_id)

cluster_ids <- unique(interaction_results$cluster_id)

cluster_ids <- cluster_ids[
  order(
    suppressWarnings(as.numeric(cluster_ids)),
    cluster_ids,
    na.last = TRUE
  )
]

significant_interaction <- interaction_results[
  !is.na(interaction_results$FDR) &
    interaction_results$FDR < fdr_threshold,
  ,
  drop = FALSE
]

summary_per_cluster <- do.call(
  rbind,
  lapply(
    cluster_ids,
    function(cluster_id_current) {
      cluster_all <- interaction_results[
        interaction_results$cluster_id == cluster_id_current,
        ,
        drop = FALSE
      ]

      cluster_significant <- significant_interaction[
        significant_interaction$cluster_id == cluster_id_current,
        ,
        drop = FALSE
      ]

      data.frame(
        cluster_id = cluster_id_current,
        tested_genes = nrow(cluster_all),
        significant_genes_FDR_lt_0.05 = nrow(cluster_significant),
        positive_logFC = sum(
          cluster_significant$logFC > 0,
          na.rm = TRUE
        ),
        negative_logFC = sum(
          cluster_significant$logFC < 0,
          na.rm = TRUE
        ),
        zero_logFC = sum(
          cluster_significant$logFC == 0,
          na.rm = TRUE
        ),
        minimum_FDR = if (nrow(cluster_significant) > 0L) {
          min(cluster_significant$FDR, na.rm = TRUE)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )
)

summary_per_cluster <- summary_per_cluster[
  order(
    suppressWarnings(as.numeric(summary_per_cluster$cluster_id)),
    summary_per_cluster$cluster_id,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

if (nrow(significant_interaction) > 0L) {
  significant_interaction <- significant_interaction[
    order(
      suppressWarnings(as.numeric(significant_interaction$cluster_id)),
      significant_interaction$cluster_id,
      significant_interaction$FDR,
      -abs(significant_interaction$logFC),
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
}

significant_interaction %>%  
  select(!c(chromosome, mean_TMM_CPM_Male_Neurotypical, sd_TMM_CPM_Male_Neurotypical,
  mean_TMM_CPM_Male_ASD, sd_TMM_CPM_Male_ASD, mean_TMM_CPM_Female_Neurotypical, sd_TMM_CPM_Female_Neurotypical, mean_TMM_CPM_Female_ASD, sd_TMM_CPM_Female_ASD)) %>% 
  select(!c(sd_percent_positive_spots_Male_Neurotypical, sd_percent_positive_spots_Male_ASD, sd_percent_positive_spots_Female_Neurotypical, sd_percent_positive_spots_Female_ASD)) %>%  
  mutate(max_mean_percent = pmax(mean_percent_positive_spots_Male_Neurotypical, mean_percent_positive_spots_Male_ASD, mean_percent_positive_spots_Female_Neurotypical, mean_percent_positive_spots_Female_ASD, na.rm = TRUE)) %>% 
  pull(max_mean_percent) %>% hist

project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"

input_rdata_file <- "/home/mateusz/projects/ippas-kunevicius-spatial/results/maternalFMT_n16samples/seurat_clustering_analysis/logNormalize_vst/hvg2000/02_cca/dims20/k20/prune0067/RData/02_maternalFMT_n16samples_logNormalizeVst_hvg2000_ccaDims20K20PruneSNN0067_res010to100by010_multiClusteringAndUmap.RData"

annotation_rds_file <- "/home/mateusz/projects/ippas-kunevicius-spatial/data/ensembl115/ensembl115_mouse_proteinCodingGenes.rds"

annotation_metadata_file <- "/home/mateusz/projects/ippas-kunevicius-spatial/data/ensembl115/ensembl115_mouse_proteinCodingGenes_metadata.tsv"

functions_file <- "/home/mateusz/projects/ippas-kunevicius-spatial/preprocessing/src/pseudobulkPerCluster_geneCounts_edgeR/functions_pseudobulkPerCluster_geneCounts_edgeR_twoFactorInteraction.R"

source(functions_file)

seurat_loaded <- load_single_seurat_object_from_rdata(
  input_rdata_file = input_rdata_file,
  requested_object_name = NULL
)

raw_counts <- extract_raw_counts_from_seurat_layers(
  seurat_object = seurat_loaded$object,
  assay_name = "RNA"
)$counts

annotation <- read_and_validate_ensembl115_annotation(
  annotation_rds_file = annotation_rds_file,
  annotation_metadata_file = annotation_metadata_file
)

mapped <- map_seurat_features_to_analysis_annotation(
  count_matrix = raw_counts,
  annotation_all = annotation$annotation_all,
  annotation_analysis = annotation$annotation_analysis
)

cat("Wszystkie geny/features w danych:", nrow(raw_counts), "\n")
cat("Protein-coding chr1-19,X w BioMart:", nrow(annotation$annotation_analysis), "\n")
cat("Protein-coding chr1-19,X OBECNE W DANYCH:", nrow(mapped$counts), "\n")
cat("Typ identyfikatorów:", mapped$feature_identifier_type, "\n")


global_results <- read.delim(
  "/home/mateusz/projects/ippas-kunevicius-spatial/results/maternalFMT_n16samples/pseudobulk_geneCounts_edgeR_sexInteraction/proteinCodingGenes_autosomalAndXChromosome/sexInteraction_overall_ASD_vs_Neurotypical_fullResults.tsv"
)

global_significant <- global_results %>%
  filter(
    FDR < 0.05,
    abs(logFC) > 0.5
  )

cat("Significant genes:", nrow(global_significant), "\n")
cat("UP:", sum(global_significant$logFC > 0), "\n")
cat("DOWN:", sum(global_significant$logFC < 0), "\n")

global_results %>% head


