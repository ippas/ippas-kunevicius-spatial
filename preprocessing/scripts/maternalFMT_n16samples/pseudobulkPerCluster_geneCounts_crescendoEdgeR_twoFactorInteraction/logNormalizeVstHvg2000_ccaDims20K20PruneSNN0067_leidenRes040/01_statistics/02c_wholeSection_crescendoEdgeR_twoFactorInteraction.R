#!/usr/bin/env Rscript
# =============================================================================
# 02c_wholeSection_crescendoEdgeR_twoFactorInteraction.R
# One pseudobulk per biological sample/whole section (n=16).
# INPUT: Crescendo-corrected per-cluster pseudobulk counts (NOT raw Seurat counts).
# Sum nonnegative corrected counts across disjoint clusters, by sample & gene.
# MODEL: ~ fmt_donor_group * sex; same nine contrasts as per-cluster analysis.
# FDR: BH per test across tested genes (not across the nine tests).
# Heatmap: donor main effect FDR < 0.05, no logFC cutoff, log2+QN -> row Z.
# WARNING: Summing independently corrected cluster-level counts is NOT
# mathematically equivalent to correcting a whole-section pseudobulk directly.
# =============================================================================
options(stringsAsFactors = FALSE)
SCRIPT_BUILD <- "20261019_02c_wholeSection_v2_heatmapMatrixFix"

# ---- Paths and settings ------------------------------------------------------
project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"
clustering_name <- "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"
expected_samples <- 16L
expected_clusters <- 16L
minimum_cell_n <- 2L
fdr_heatmap <- 0.05
filter_params <- list(min.count=10, min.total.count=15, large.n=10, min.prop=0.7)
input_root <- file.path(project_root,"results",dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,"01_statistics")
input_rds <- file.path(input_root,"01_crescendoBatchCorrection","07_objects",
  paste0(dataset_name,"_",clustering_name,
         "_crescendoBatchCorrection_pseudobulkCounts_beforeEdgeR.rds"))
annotation_rds <- file.path(project_root,"data","ensembl115",
                            "ensembl115_mouse_proteinCodingGenes.rds")
annotation_meta <- file.path(project_root,"data","ensembl115",
                            "ensembl115_mouse_proteinCodingGenes_metadata.tsv")
output_dir <- file.path(input_root,"02c_wholeSection_edgeR")
prefix <- paste0("02c_",dataset_name,"_wholeSection_Crescendo_edgeR")
dirs <- setNames(file.path(output_dir,c("01_inputs","02_design","03_results",
                                    "04_summary","05_heatmap","06_objects","07_logs")),
                 c("inputs","design","results","summary","heatmap","objects","logs"))
for (d in c(output_dir,dirs)) dir.create(d,recursive=TRUE,showWarnings=FALSE)

# ---- Dependencies ------------------------------------------------------------
needed <- c("edgeR","limma","dplyr","tibble","openxlsx",
            "preprocessCore","ComplexHeatmap","circlize","grid")
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Install missing R packages: ",paste(missing,collapse=", "))
for (f in c(input_rds,annotation_rds,annotation_meta))
  if(!file.exists(f)) stop("Missing input: ",f)
write_tsv <- function(x,path) {
  con <- if(grepl("\\.gz$",path)) gzfile(path,"wt") else file(path,"wt")
  on.exit(close(con)); utils::write.table(x,con,sep="\t",quote=FALSE,
    row.names=FALSE,col.names=TRUE,na="NA")
}
stopif <- function(cond,...) if(isTRUE(cond)) stop(...,call.=FALSE)
message("Build: ",SCRIPT_BUILD,"\nInput: ",input_rds,"\nOutput: ",output_dir)

# ---- Load and validate sample / cluster data --------------------------------
x <- readRDS(input_rds)
stopif(!all(c("sample_metadata","group_table",
              "crescendo_corrected_counts_by_cluster") %in% names(x)),
       "Required components missing from Crescendo RDS")
meta <- as.data.frame(x$sample_metadata)
stopif(!all(c("sample_ID","fmt_donor_group","sex") %in% names(meta)),
       "Sample metadata lack required fields")
meta$sample_ID <- as.character(meta$sample_ID)
meta$fmt_donor_group <- factor(as.character(meta$fmt_donor_group),
                              levels=c("Neurotypical","ASD"))
meta$sex <- factor(as.character(meta$sex),levels=c("Male","Female"))
stopif(anyNA(meta[,c("sample_ID","fmt_donor_group","sex")]),
       "Missing/unexpected sample metadata")
stopif(anyDuplicated(meta$sample_ID)>0L,"Duplicated sample IDs")
stopif(nrow(meta)!=expected_samples,"Expected ",expected_samples," samples")
meta <- meta[order(meta$sex,meta$fmt_donor_group,meta$sample_ID),,drop=FALSE]
rownames(meta) <- meta$sample_ID
meta$cell <- factor(paste(meta$sex,meta$fmt_donor_group,sep="_"),
  levels=c("Male_Neurotypical","Male_ASD","Female_Neurotypical","Female_ASD"))
cell_n <- table(meta$cell)
stopif(any(cell_n < minimum_cell_n),"Insufficient samples in at least one cell")

clusters <- x$crescendo_corrected_counts_by_cluster
stopif(!is.list(clusters) || length(clusters)!=expected_clusters,
       "Expected ",expected_clusters," cluster count matrices")
stopif(is.null(names(clusters)) || anyDuplicated(names(clusters))>0L,
       "Cluster list must have unique names")
stopif(any(vapply(clusters,is.null,logical(1))),"NULL cluster count matrix")
cluster_ids <- names(clusters)
if(all(!is.na(suppressWarnings(as.integer(cluster_ids)))))
  cluster_ids <- cluster_ids[order(as.integer(cluster_ids))]
features <- rownames(clusters[[cluster_ids[1]]])
stopif(is.null(features) || anyDuplicated(features)>0L,"Invalid gene IDs")

# Validate that every cluster has the same genes and all 16 sample columns.
for (cid in cluster_ids) {
  m <- clusters[[cid]]
  stopif(is.null(rownames(m)) || !setequal(rownames(m),features),
         "Gene mismatch in cluster ",cid)
  stopif(is.null(colnames(m)) || anyDuplicated(colnames(m))>0L ||
         !setequal(colnames(m),meta$sample_ID),
         "Sample mismatch in cluster ",cid)
  stopif(any(!is.finite(m)),"Non-finite values in cluster ",cid)
  stopif(any(m < 0),"Negative values in cluster ",cid)
}
# Check group table for coverage of every sample x cluster and spot counts.
gt <- as.data.frame(x$group_table)
stopif(!all(c("sample_ID","cluster_id","n_spots") %in% names(gt)),
       "Invalid Crescendo group_table")
gt$sample_ID <- as.character(gt$sample_ID)
gt$cluster_id <- as.character(gt$cluster_id)
stopif(anyDuplicated(paste(gt$sample_ID,gt$cluster_id))>0L,
       "Duplicate sample/cluster in group_table")
expected_pairs <- as.vector(outer(meta$sample_ID,cluster_ids,paste,sep="::"))
provided_pairs <- paste(gt$sample_ID,gt$cluster_id,sep="::")
stopif(!setequal(expected_pairs,provided_pairs),
       "group_table does not cover exactly all sample x cluster combinations")
stopif(anyNA(gt$n_spots) || any(gt$n_spots < 0),"Invalid spot counts")
spot_count <- stats::aggregate(n_spots ~ sample_ID,data=gt,FUN=sum)
meta$n_spots_whole_section <- spot_count$n_spots[match(meta$sample_ID,spot_count$sample_ID)]
stopif(any(meta$n_spots_whole_section<=0),"Sample has zero spots")

# ---- Gene annotation ---------------------------------------------------------
am <- utils::read.delim(annotation_meta,check.names=FALSE)
stopif(nrow(am)!=1L || !all(c("ensembl_release","genome_assembly") %in% names(am)) ||
       as.integer(am$ensembl_release[1])!=115L ||
       as.character(am$genome_assembly[1])!="GRCm39",
       "Annotation metadata not Ensembl 115 / GRCm39")
ann <- as.data.frame(readRDS(annotation_rds))
stopif(!all(c("chromosome","ensembl_gene_id","gene","gene_biotype") %in% names(ann)),
       "Unexpected annotation columns")
ann$ensembl_gene_id <- sub("\\.[0-9]+$","",trimws(as.character(ann$ensembl_gene_id)))
ann$gene <- trimws(as.character(ann$gene))
ann$chromosome <- sub("^chr","",trimws(as.character(ann$chromosome)))
ann <- ann[!is.na(ann$ensembl_gene_id) & ann$ensembl_gene_id!="" &
             ann$gene_biotype=="protein_coding" &
             ann$chromosome %in% c(as.character(1:19),"X"),
           c("chromosome","ensembl_gene_id","gene","gene_biotype")]
ann <- ann[!duplicated(ann$ensembl_gene_id),,drop=FALSE]
normalized <- sub("\\.[0-9]+$","",features)
identifier_type <- if(mean(grepl("^ENSMUSG[0-9]+$",normalized))>=0.5)
  "ensembl_gene_id" else "gene_symbol"
if(identifier_type=="ensembl_gene_id") {
  ai <- match(normalized,ann$ensembl_gene_id)
} else {
  unique_symbols <- ann$gene[!duplicated(ann$gene) &
                               !duplicated(ann$gene,fromLast=TRUE)]
  ai <- match(features,ann$gene)
  ai[!features %in% unique_symbols] <- NA_integer_
}
map <- data.frame(feature_id=features,ensembl_gene_id=ann$ensembl_gene_id[ai],
                  gene=ann$gene[ai],chromosome=ann$chromosome[ai],
                  kept=!is.na(ai))
map$kept[duplicated(map$ensembl_gene_id) & !is.na(map$ensembl_gene_id)] <- FALSE
kept <- which(map$kept)
stopif(length(kept)==0L,"No protein-coding chr1-19/X genes mapped")
gene_ann <- map[kept,c("ensembl_gene_id","gene","chromosome")]
gene_ann$gene[is.na(gene_ann$gene) | gene_ann$gene==""] <-
  gene_ann$ensembl_gene_id[is.na(gene_ann$gene) | gene_ann$gene==""]
write_tsv(map,file.path(dirs["inputs"],paste0(prefix,"_geneMapping.tsv.gz")))

# ---- Whole-section pseudobulk: sum counts across clusters per sample ---------
# Summation takes place BEFORE edgeR filtering / TMM / transformation.
whole_counts <- matrix(0,nrow=length(kept),ncol=nrow(meta),
                       dimnames=list(gene_ann$ensembl_gene_id,meta$sample_ID))
for(cid in cluster_ids) {
  m <- clusters[[cid]][features[kept],meta$sample_ID,drop=FALSE]
  whole_counts <- whole_counts + as.matrix(m)
  message("Summed cluster ",cid)
}
rm(clusters,x)
gc(verbose=FALSE)
stopif(any(!is.finite(whole_counts)) || any(whole_counts<0),
       "Invalid whole-section pseudobulk counts")
stopif(any(colSums(whole_counts)<=0),"Sample has zero pseudobulk library size")
message("Whole-section counts: ",nrow(whole_counts)," genes x ",ncol(whole_counts)," samples")
write_tsv(cbind(gene_ann,as.data.frame(whole_counts,check.names=FALSE)),
 file.path(dirs["inputs"],paste0(prefix,"_correctedWholeSectionCounts.tsv.gz")))
write_tsv(meta[,c("sample_ID","fmt_donor_group","sex","cell","n_spots_whole_section")],
 file.path(dirs["design"],paste0(prefix,"_sampleMetadata.tsv")))

# ---- Design and nine contrasts ----------------------------------------------
design <- model.matrix(~fmt_donor_group*sex,data=meta)
stopif(qr(design)$rank!=ncol(design),"Design is not full rank")
stopif(!identical(colnames(design),c("(Intercept)","fmt_donor_groupASD",
                                     "sexFemale","fmt_donor_groupASD:sexFemale")),
       "Unexpected design columns")
colnames(design) <- c("Intercept","groupASD","sexFemale","groupASD_sexFemale")
rownames(design) <- meta$sample_ID
contrasts <- limma::makeContrasts(
 Interaction=groupASD_sexFemale,
 Overall_Sex_Female_vs_Male=sexFemale+0.5*groupASD_sexFemale,
 Overall_Group_ASD_vs_Neurotypical=groupASD+0.5*groupASD_sexFemale,
 ASD_Male_vs_Neurotypical_Male=groupASD,
 ASD_Female_vs_Neurotypical_Female=groupASD+groupASD_sexFemale,
 Neurotypical_Female_vs_Neurotypical_Male=sexFemale,
 ASD_Female_vs_ASD_Male=sexFemale+groupASD_sexFemale,
 ASD_Female_vs_Neurotypical_Male=groupASD+sexFemale+groupASD_sexFemale,
 ASD_Male_vs_Neurotypical_Female=groupASD-sexFemale,
 levels=design)
test_ids <- colnames(contrasts)
sheets <- c("Interaction","Overall_Sex","Overall_Group","Group_in_Male",
            "Group_in_Female","Sex_in_Neurotypical","Sex_in_ASD",
            "ASD_F_vs_NT_M","ASD_M_vs_NT_F")
test_defs <- data.frame(test_id=test_ids,sheet_name=sheets,
 comparison=c("(ASD-NT) Female minus (ASD-NT) Male",
              "Female minus Male averaged equally across donor groups",
              "ASD minus NT averaged equally across sexes",
              "ASD Male minus NT Male","ASD Female minus NT Female",
              "NT Female minus NT Male","ASD Female minus ASD Male",
              "ASD Female minus NT Male","ASD Male minus NT Female"))
write_tsv(data.frame(sample_ID=rownames(design),design,check.names=FALSE),
 file.path(dirs["design"],paste0(prefix,"_design.tsv")))
write_tsv(data.frame(coefficient=rownames(contrasts),contrasts,check.names=FALSE),
 file.path(dirs["design"],paste0(prefix,"_contrasts.tsv")))
write_tsv(test_defs,file.path(dirs["design"],paste0(prefix,"_tests.tsv")))

# ---- edgeR ------------------------------------------------------------------
dge0 <- edgeR::DGEList(counts=whole_counts)
keep_expr <- edgeR::filterByExpr(dge0,design=design,
 min.count=filter_params$min.count,min.total.count=filter_params$min.total.count,
 large.n=filter_params$large.n,min.prop=filter_params$min.prop)
write_tsv(data.frame(gene_ann,total_count=rowSums(whole_counts),
                      samples_nonzero=rowSums(whole_counts>0),
                      passed_filterByExpr=as.logical(keep_expr)),
 file.path(dirs["inputs"],paste0(prefix,"_filterByExpr.tsv.gz")))
stopif(!any(keep_expr),"No genes passed filterByExpr")
dge <- dge0[keep_expr,,keep.lib.sizes=FALSE]
dge <- edgeR::calcNormFactors(dge,method="TMM")
dge$samples$sample_ID <- meta$sample_ID
dge$samples$fmt_donor_group <- as.character(meta$fmt_donor_group)
dge$samples$sex <- as.character(meta$sex)
dge$samples$n_spots <- meta$n_spots_whole_section
sample_qc <- data.frame(sample_ID=meta$sample_ID,
 fmt_donor_group=as.character(meta$fmt_donor_group),sex=as.character(meta$sex),
 n_spots=meta$n_spots_whole_section,library_size=dge$samples$lib.size,
 TMM_factor=dge$samples$norm.factors,
 effective_library_size=dge$samples$lib.size*dge$samples$norm.factors)
write_tsv(sample_qc,file.path(dirs["summary"],paste0(prefix,"_sampleLibraryQC.tsv")))
tmm_cpm <- edgeR::cpm(dge,log=FALSE,normalized.lib.sizes=TRUE)
write_tsv(cbind(gene_ann[match(rownames(tmm_cpm),gene_ann$ensembl_gene_id),],
                as.data.frame(tmm_cpm,check.names=FALSE)),
 file.path(dirs["inputs"],paste0(prefix,"_TMM_CPM_testedGenes.tsv.gz")))
dge <- edgeR::estimateDisp(dge,design=design,robust=TRUE)
fit <- edgeR::glmQLFit(dge,design=design,robust=TRUE)

# Descriptive QN input is calculated on all mapped genes before selection,
# then subset for heatmap. Never fed to edgeR.
qn <- preprocessCore::normalize.quantiles(log2(whole_counts+1),copy=TRUE)
dimnames(qn) <- dimnames(whole_counts)
write_tsv(cbind(gene_ann,as.data.frame(qn,check.names=FALSE)),
 file.path(dirs["inputs"],paste0(prefix,"_log2plus1_QN_allMappedGenes.tsv.gz")))

results <- setNames(vector("list",length(test_ids)),test_ids)
qlf_tests <- setNames(vector("list",length(test_ids)),test_ids)
for (i in seq_along(test_ids)) {
  id <- test_ids[i]
  tt <- edgeR::glmQLFTest(fit,contrast=contrasts[,id])
  qlf_tests[[id]] <- tt
  tbl <- edgeR::topTags(tt,n=Inf,sort.by="PValue")$table
  tbl <- data.frame(ensembl_gene_id=rownames(tbl),tbl,check.names=FALSE)
  tbl <- merge(tbl,gene_ann,by="ensembl_gene_id",sort=FALSE)
  tbl <- tbl[order(tbl$PValue,-abs(tbl$logFC)),,drop=FALSE]
  tbl$regulation <- ifelse(tbl$logFC>0,"UP",ifelse(tbl$logFC<0,"DOWN","ZERO"))
  tbl <- tbl[,c("ensembl_gene_id","gene","chromosome","logFC","logCPM",
               "F","PValue","FDR","regulation")]
  rownames(tbl) <- NULL
  results[[id]] <- tbl
  write_tsv(tbl,file.path(dirs["results"],paste0(prefix,"_",id,"_full.tsv.gz")))
}

# ---- XLSX workbook -----------------------------------------------------------
wb <- openxlsx::createWorkbook()
header <- openxlsx::createStyle(textDecoration="bold",fgFill="#D9EAF7")
for(i in seq_along(test_ids)) {
  sh <- sheets[i]; df <- results[[test_ids[i]]]
  openxlsx::addWorksheet(wb,sh)
  openxlsx::writeData(wb,sh,df,headerStyle=header,withFilter=TRUE)
  openxlsx::freezePane(wb,sh,firstActiveRow=2,firstActiveCol=3)
  openxlsx::setColWidths(wb,sh,1:ncol(df),widths=c(20,18,12,12,12,12,16,16,12))
  for(nm in c("PValue","FDR")) {
    style <- openxlsx::createStyle(numFmt="0.00000E+00")
    openxlsx::addStyle(wb,sh,style,rows=2:(nrow(df)+1),
                      cols=match(nm,names(df)),gridExpand=TRUE,stack=TRUE)
  }
}
xlsx_path <- file.path(dirs["results"],paste0(prefix,"_nineTests.xlsx"))
openxlsx::saveWorkbook(wb,xlsx_path,overwrite=TRUE)

# ---- Summary FDR counts, UP / DOWN ------------------------------------------
thresholds <- c(0.05,0.01,0.001)
summary_list <- lapply(test_ids,function(id) {
  r <- results[[id]]
  do.call(rbind,lapply(thresholds,function(th) {
    sig <- !is.na(r$FDR) & r$FDR<th
    data.frame(test_id=id,FDR_threshold=th,genes_tested=nrow(r),
               genes_significant=sum(sig),UP=sum(sig & r$logFC>0),
               DOWN=sum(sig & r$logFC<0),
               significant_abs_logFC_ge_0.5=sum(sig & abs(r$logFC)>=0.5))
  }))
})
summary_fdr <- do.call(rbind,summary_list)
rownames(summary_fdr) <- NULL
summary_path <- file.path(dirs["summary"],paste0(prefix,"_FDRcounts.tsv"))
write_tsv(summary_fdr,summary_path)
message("\nFDR counts:\n")
print(summary_fdr,row.names=FALSE)

# ---- Simple donor main-effect heatmap ----------------------------------------
donor_id <- "Overall_Group_ASD_vs_Neurotypical"
donor <- results[[donor_id]]
sig <- donor[is.finite(donor$FDR) & donor$FDR<fdr_heatmap,,drop=FALSE]
sig <- sig[order(sig$logFC,decreasing=TRUE),,drop=FALSE]
heatmap_results_path <- file.path(dirs["heatmap"],
                                  paste0(prefix,"_donor_FDR005_genes.tsv"))
write_tsv(sig,heatmap_results_path)
heatmap_path <- file.path(dirs["heatmap"],paste0(prefix,"_donor_FDR005_heatmap.pdf"))
if (nrow(sig)==0L) {
  message("No donor-group genes with FDR < ",fdr_heatmap,"; heatmap not generated")
  writeLines("No significant donor main-effect genes at FDR < 0.05; no heatmap generated.",
             file.path(dirs["heatmap"],paste0(prefix,"_noSignificantGenes.txt")))
} else {
  expr <- qn[sig$ensembl_gene_id,meta$sample_ID,drop=FALSE]
  z <- t(scale(t(expr)))
  z[!is.finite(z)] <- 0
  # Clip values in place: pmin/pmax with a scalar first may drop matrix dimensions.
  z[z < -2] <- -2
  z[z > 2] <- 2
  stopifnot(is.matrix(z), identical(dim(z), dim(expr)))
  symbols <- ifelse(is.na(sig$gene)|sig$gene=="",sig$ensembl_gene_id,sig$gene)
  rownames(z) <- make.unique(symbols)
  # Save precisely the heatmap's values for reproducibility.
  write_tsv(data.frame(ensembl_gene_id=sig$ensembl_gene_id,gene=symbols,
                       as.data.frame(z,check.names=FALSE)),
   file.path(dirs["heatmap"],paste0(prefix,"_donor_FDR005_rowZ.tsv.gz")))
  cell <- factor(as.character(meta$cell),levels=levels(meta$cell))
  top <- ComplexHeatmap::HeatmapAnnotation(
    Group=as.character(meta$fmt_donor_group),Sex=as.character(meta$sex),
    col=list(Group=c(Neurotypical="grey70",ASD="black"),
             Sex=c(Male="#3C78D8",Female="#C050A0")))
  stars <- ifelse(sig$FDR<0.001,"***",ifelse(sig$FDR<0.01,"**","*"))
  ra <- ComplexHeatmap::rowAnnotation(
    `FDR`=ComplexHeatmap::anno_text(stars,
      gp=grid::gpar(fontsize=8,fontface="bold"),
      width=grid::unit(8,"mm")),
    `log2FC`=ComplexHeatmap::anno_text(sprintf("%.2f",sig$logFC),
      gp=grid::gpar(fontsize=7),width=grid::unit(13,"mm")),
    Direction=ifelse(sig$logFC>0,"UP","DOWN"),
    col=list(Direction=c(UP="#B2182B",DOWN="#2166AC")))
  # Pearson distance; no cluster-specific distance modifications.
  row_dist <- function(m) {
    cor_mat <- suppressWarnings(stats::cor(t(m),use="pairwise.complete.obs"))
    if(is.null(dim(cor_mat))) cor_mat <- matrix(1,nrow=nrow(m),ncol=nrow(m))
    cor_mat[!is.finite(cor_mat)] <- 0
    diag(cor_mat) <- 1
    distance_matrix <- 1 - cor_mat
    distance_matrix[distance_matrix < 0] <- 0
    diag(distance_matrix) <- 0
    stats::as.dist(distance_matrix)
  }
  hm <- ComplexHeatmap::Heatmap(z,name="Row Z-score",
    col=circlize::colorRamp2(c(-2,0,2),c("#2166AC","white","#B2182B")),
    top_annotation=top,right_annotation=ra,
    column_split=cell,cluster_columns=FALSE,cluster_column_slices=FALSE,
    cluster_rows=nrow(z)>1L,
    clustering_distance_rows=if(nrow(z)>1L) row_dist else "euclidean",
    clustering_method_rows="complete",
    show_row_names=TRUE,row_names_side="right",
    row_names_gp=grid::gpar(fontsize=if(nrow(z)>150) 5 else 8),
    column_names_rot=45,column_names_gp=grid::gpar(fontsize=8),
    column_title=paste0("Whole-section donor-group main effect | FDR < 0.05 | ",
                        nrow(sig)," genes | UP ",sum(sig$logFC>0),
                        " / DOWN ",sum(sig$logFC<0),
                        "\nCrescendo summed pseudobulk: log2(count+1) + QN + row Z-score"),
    column_title_gp=grid::gpar(fontsize=11,fontface="bold"),
    use_raster=nrow(z)>1000L)
  grDevices::pdf(heatmap_path,width=12,
                 height=max(8,min(100,(nrow(z)*3+90)/25.4)),useDingbats=FALSE)
  tryCatch(ComplexHeatmap::draw(hm,heatmap_legend_side="right",
                               annotation_legend_side="right",merge_legends=TRUE),
           finally=grDevices::dev.off())
  message("Heatmap: ",heatmap_path)
}

# ---- Save native R objects / log --------------------------------------------
saveRDS(list(script_build=SCRIPT_BUILD,metadata=meta,
             whole_section_corrected_counts=whole_counts,
             annotation=gene_ann,keep_genes=keep_expr,design=design,
             contrasts=contrasts,test_definitions=test_defs,
             dge=dge,fit=fit,tests=qlf_tests,results=results,
             tmm_cpm=tmm_cpm,log2_qn=qn,FDR_summary=summary_fdr),
        file.path(dirs["objects"],paste0(prefix,"_analysisObjects.rds")),
        compress="xz")
writeLines(c(paste("script_build",SCRIPT_BUILD),
             paste("input",normalizePath(input_rds)),
             paste("output",normalizePath(output_dir)),
             "input_type: sum of Crescendo-corrected per-cluster pseudobulk counts",
             "experimental_unit: biological sample (whole brain section)",
             "model: ~ fmt_donor_group * sex",
             "sex reference: Male; group reference: Neurotypical",
             "FDR: BH separately across tested genes within each contrast",
             "donor_main_effect: average of ASD-NT in males and females (equal weights)",
             "visualization_only: log2(count+1) -> quantile normalization -> gene Z-score",
             "caveat: per-cluster-corrected count sums not equivalent to whole-section correction"),
           file.path(dirs["logs"],paste0(prefix,"_README.txt")))
writeLines(capture.output(sessionInfo()),
 file.path(dirs["logs"],paste0(prefix,"_sessionInfo.txt")))
message("\nDONE\nResults: ",xlsx_path,"\nSummary: ",summary_path,
        "\nOutput: ",output_dir)