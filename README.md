# ippas-kunevicius-spatial

#### Visium spatial transcriptomics · *Mus musculus* · coronal brain sections
Offspring brains after maternal faecal microbiota transplantation (FMT) from autism spectrum disorder (ASD) or neurotypical (NT) human donors. 20 sections were sequenced; 16 sections (8 NT, 8 ASD; both sexes) passed quality control and were analysed.

---

## Repository layout

- [`preprocessing/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing)
  - [`spaceranger_analysis/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/spaceranger_analysis) – Space Ranger alignment (shell)
  - [`scripts/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/scripts) – analysis scripts, one folder per analysis
  - [`src/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/src) – shared R functions sourced by the scripts
- `data/` – sample metadata, gene annotation (Ensembl 115); not tracked in git
- `results/` – outputs, mirroring the script folders; not tracked in git
- `raw/`, `tools/` – FASTQ, images, Space Ranger; not tracked in git

Result folders follow the same path as the script folders:
`preprocessing/scripts/maternalFMT_n16samples/<analysis>/…` → `results/maternalFMT_n16samples/<analysis>/…`.

Shared functions ([`preprocessing/src/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/src)):

| Folder | Used for |
|-------------------------------|----------|
| [`seurat_clustering_analysis/logNormalize_vst/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/src/seurat_clustering_analysis/logNormalize_vst) | preprocessing, integration, clustering, UMAP, validation metrics |
| [`manual_cluster_annotation/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/src/manual_cluster_annotation) | spatial and UMAP plots of the annotated clusters |
| [`gene_specificity/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/src/gene_specificity) | tau, Gini and Shannon specificity of marker genes |
| [`pseudobulkPerCluster_geneCounts_edgeR/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/src/pseudobulkPerCluster_geneCounts_edgeR) | edgeR helpers and spatial expression maps |

---

## Analysis steps

### 1. Read alignment
[`preprocessing/spaceranger_analysis/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/spaceranger_analysis)
- `spaceranger count` for each section, mouse reference `refdata-gex-GRCm39-2024-A`, automatic image alignment.
- Sample sheet: `data/metadata_autismFMT.tsv` (sample ID, slide ID, donor group, sex).

### 2. Section quality control
[`preprocessing/scripts/raw_UMI_spatial_similarity/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/scripts/raw_UMI_spatial_similarity)
- Raw UMI counts of all barcodes (including off-tissue) plotted on the tissue image, as raw values and as log2(UMI + 1).
- Sections whose UMI pattern did not match the tissue image were excluded: **20_1F, 12_3F, 15_1M, 20_3M**.

### 3. Normalisation and dimensionality reduction
[`preprocessing/scripts/maternalFMT_n16samples/seurat_clustering_analysis/logNormalize_vst/hvg2000/00_preprocessing/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/scripts/maternalFMT_n16samples/seurat_clustering_analysis/logNormalize_vst/hvg2000/00_preprocessing)
- In-tissue spots of the 16 sections merged into one Seurat v5 object.
- LogNormalize (scale factor 10,000), 2,000 highly variable genes (VST), PCA (50 components).
- Output: parent object used by all clustering configurations.

### 4. Integration and clustering grid
[`preprocessing/scripts/maternalFMT_n16samples/seurat_clustering_analysis/logNormalize_vst/hvg2000/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/scripts/maternalFMT_n16samples/seurat_clustering_analysis/logNormalize_vst/hvg2000)
- Integration: none (PCA), CCA, RPCA, Harmony, FastMNN, scVI.
- For each method: dims 20 or 30, SNN graph with k = 20 or 30 (prune.SNN = 1/15).
- Clustering: Louvain, Louvain refined, SLM, Leiden; resolutions 0.1–1.0 (step 0.1) → 960 solutions.
- For each solution: UMAP, spatial plots and validation metrics (transcriptomic silhouette, CHAOS, PAS, ARI between resolutions).
- Output per configuration: `results/…/hvg2000/02_<method>/dims<d>/k<k>/prune0067/`.

### 5. Selection of the clustering solution
[`preprocessing/scripts/maternalFMT_n16samples/manual_cluster_annotation/clustering_selection_summary/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/scripts/maternalFMT_n16samples/manual_cluster_annotation/clustering_selection_summary)
- Validation metrics of all 960 solutions collected into one table.
- Selected: highest mean transcriptomic silhouette – **CCA, dims 20, k 20, Leiden, resolution 0.4 (16 clusters)**.

### 6. Anatomical annotation and marker genes
[`preprocessing/scripts/maternalFMT_n16samples/manual_cluster_annotation/logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/scripts/maternalFMT_n16samples/manual_cluster_annotation/logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040)
- Spatial and UMAP plots of the selected clustering; clusters compared with coronal plate 100960224 of the Allen Mouse Brain Atlas.
- Marker genes: `FindAllMarkers` on the RNA assay (Wilcoxon, positive markers, min.pct 0.25, log2FC ≥ 0.25, Bonferroni-adjusted P < 0.05).
- Cluster specificity of each marker (tau, Gini, Shannon) calculated from mean expression per cluster.

### 7. Batch correction and differential expression
[`preprocessing/scripts/maternalFMT_n16samples/pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction/logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040/01_statistics/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/scripts/maternalFMT_n16samples/pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction/logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040/01_statistics)

*Crescendo batch correction*
- Raw spot counts of protein-coding genes (Ensembl 115) detected in ≥ 10% of spots in at least one cluster (10,203 genes).
- Batch variable: slide ID. Biological structure: Harmony soft clusters (20 PCs, nclust = 20).
- Crescendo v1.0 run on all spots together, in chunks of 100 genes (prop = 0.1, min_cells = 50, alpha = 0).
- Corrected spot counts summed per sample × cluster (pseudobulk).

*edgeR*
- Per cluster: genes on chr 1–19 and X (10,182), `filterByExpr` (min.count 10, min.total.count 15), TMM normalisation.
- Model `~ donor group * sex` (reference: Neurotypical, Male); `estimateDisp` and `glmQLFit` with robust = TRUE.
- Nine contrasts with `glmQLFTest`: interaction, main effects of sex and donor group, six pairwise comparisons.
- Benjamini–Hochberg FDR within each cluster and contrast; significant at FDR < 0.05.
- Main result, donor-group effect: **79 cluster–gene results, 73 unique genes, 6 clusters**.

*Visualisation (same folder)*
- Heatmap: log2(count + 1), quantile normalisation within cluster, row Z-score; row clustering with distance 0.8 × (2 − r) within and 1.2 × (2 − r) between clusters.
- Spatial maps: log2(1 + corrected count) per spot, clipped at the 99th percentile and scaled to 0–1.
- Bar plots: log2(TMM-CPM + 1) and TMM-CPM per sample, mean ± SEM, stars from edgeR FDR.

### 8. Functional enrichment
[`preprocessing/scripts/maternalFMT_n16samples/pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction/logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040/01_statistics/functional_analysis/`](https://github.com/ippas/ippas-kunevicius-spatial/tree/master/preprocessing/scripts/maternalFMT_n16samples/pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction/logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040/01_statistics/functional_analysis)
- Python script `run_functional_analysis.py`; input gene list `input_DE_genes.tsv`; MSigDB gene sets bundled as `msigdb_v2026.1_sysdata.rda`; dependencies in `requirements.txt`.
- Over-representation analysis of the 73 genes, one-sided hypergeometric test.
- Gene sets: MSigDB 2026.1 – mouse (GO BP, GO MF, GO CC, Reactome, WikiPathways, Hallmarks) and human (KEGG, HPO).
- Background: all genes of each database; sets ≥ 10 genes; ≥ 2 overlapping genes; BH-FDR per database, significant at FDR < 0.10.
- Result: only GO Cellular Component significant (10 synaptic and membrane terms).

---

## Software

R 4.5.2 · Seurat 5.5.1 · Harmony 2.0.5 · Crescendo 1.0 · edgeR 4.8.2 · Python 3.14.4 (SciPy 1.18.1) · Space Ranger 4.1.0
