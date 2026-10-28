# Functional analysis – methods

## Input
- Differentially expressed genes from the cluster-wise pseudobulk edgeR analysis
  (maternal FMT model, ASD vs neurotypical donor group, model donor group × sex, Crescendo-corrected data),
  FDR < 0.05, no log2FC threshold.
- 79 cluster–gene results = **73 unique genes** (37 up, 36 down in ASD-donor offspring).
  A gene significant in several clusters was counted once; direction = sign of log2FC.
- File: `input_DE_genes.tsv` (gene, cluster, log2FC). It was parsed from the heatmap PDF of the DE analysis;
  for a final version use the edgeR result tables directly (same columns `gene`, `log2FC`).
- Three gene lists were tested: all changed genes, up-regulated genes, down-regulated genes.

## Gene sets
- **MSigDB v2026.1**, taken from the data of the `msigdbdf` R package (release v26.1.0,
  github.com/igordot/msigdbdf, file `R/sysdata.rda`; the data source of the `msigdbr` package).
  The script downloads the file automatically and checks that it is identical to the one used for the presentation.
- Mouse gene sets (*Mus musculus*, MSigDB mouse collections):
  GO Biological Process, GO Molecular Function, GO Cellular Component (M5),
  Reactome, WikiPathways (M2 canonical pathways), Hallmarks (MH).
- Human gene sets (no mouse version in MSigDB), mouse symbols matched in upper case:
  KEGG (legacy + MEDICUS), Human Phenotype Ontology (HPO).

## Statistics
- Over-representation analysis (ORA): one-sided hypergeometric test, P(X ≥ k)
  (equivalent to a one-sided Fisher exact test).
- Background (universe): all genes annotated in the given database.
- Gene sets with ≥ 10 genes were tested, no upper size limit; terms with < 2 overlapping genes were assigned p = 1.
- Multiple testing: Benjamini–Hochberg FDR over **all tested terms of a database**, separately for each gene list.
- Significance: FDR < 0.10.
- Software: Python 3.13, pandas, numpy, scipy, pyreadr, matplotlib (exact versions in `results/run_info.txt`).

## Results (all 73 genes)
| Database | Terms tested | FDR < 0.10 | Top term (FDR) |
|---|---|---|---|
| GO Biological Process | 5,681 | 0 | Radial glial cell differentiation (0.12) |
| GO Molecular Function | 1,245 | 0 | Monoacylglycerol lipase activity (1) |
| **GO Cellular Component** | **758** | **10** | **Synaptic membrane, 11 genes (0.003)** |
| Reactome | 1,053 | 0 | Lysosome vesicle biogenesis (0.29) |
| KEGG | 589 | 0 | Endocytosis (1) |
| WikiPathways | 179 | 0 | Oxidative damage response (1) |
| MSigDB Hallmarks | 50 | 0 | UV response down (0.16) |
| Human Phenotype Ontology | 4,384 | 0 | Generalized myoclonic-atonic seizure (0.11) |

- GO Cellular Component, 10 significant terms: synaptic membrane, postsynaptic membrane, neuron-to-neuron synapse,
  glutamatergic synapse, postsynapse, plasma membrane region, extrinsic component of synaptic membrane,
  climbing fiber–Purkinje cell synapse, recycling endosome, cell trailing edge.
  Together they contain **21 unique genes** (29% of the changed genes).
- Down-regulated genes alone: GO CC neuron-to-neuron synapse (8 genes, FDR = 0.007);
  HPO generalized myoclonic-atonic seizure and myoclonic absence seizure (FDR = 0.012; Scn1a, Syngap1, Ogdhl / Kcnma1).
- Up-regulated genes alone: no significant term.

## Notes and limitations
- Compared with Enrichr, FDR values here are more conservative: Enrichr corrects only over the terms that share
  at least one gene with the input list, here the correction is over all tested terms. Gene-set versions also differ
  (Enrichr libraries are based on human annotations). The direction of the results is the same (synaptic terms).
- The background is the whole annotated genome of each database, not the genes detected in the Visium data;
  a background of detected genes would be more appropriate and is easy to add once the full edgeR tables are available.
- Genes from different clusters are pooled (too few genes per cluster for cluster-specific enrichment).

## Opis po polsku (do pracy / prezentacji)
Analizę nadreprezentacji (ORA) przeprowadzono dla 73 unikalnych genów o zmienionej ekspresji (FDR < 0,05),
osobno dla wszystkich genów oraz dla genów o podwyższonej i obniżonej ekspresji. Wykorzystano zbiory genów
MSigDB v2026.1 (mysie: GO BP, GO MF, GO CC, Reactome, WikiPathways, Hallmarks; ludzkie: KEGG, HPO).
Istotność oceniano jednostronnym testem hipergeometrycznym (tło: wszystkie geny z adnotacją w danej bazie,
zbiory ≥ 10 genów), a wartości p korygowano metodą Benjaminiego–Hochberga osobno dla każdej bazy;
za istotne uznano FDR < 0,10. Istotne terminy uzyskano wyłącznie dla GO Cellular Component
(10 terminów, 21 unikalnych genów; najsilniej: błona synaptyczna, 11 genów, FDR = 0,003).
