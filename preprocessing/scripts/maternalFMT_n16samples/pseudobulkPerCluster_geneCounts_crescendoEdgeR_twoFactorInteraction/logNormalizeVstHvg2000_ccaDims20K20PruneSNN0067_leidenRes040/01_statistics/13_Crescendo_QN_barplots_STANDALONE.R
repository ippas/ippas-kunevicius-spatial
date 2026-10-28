#!/usr/bin/env Rscript
# Standalone TWO-GROUP barplots (Neurotypical vs ASD), Crescendo pseudobulk,
# y = log1p(corrected_gene_sum / raw_total_UMI_allGenes_in_cluster * 10000).
# Bars = mean; whiskers = mean +/- SE; dots = biological samples.
# edgeR main-effect FDR gives stars; sex remains a covariate in the original model.
# No QN, no TMM, no Z-score; these are descriptive plots, not new tests.
options(stringsAsFactors = FALSE)

# ======================== SETTINGS ==========================================
project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"
clustering_name <- "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"
fdr_threshold <- 0.05
abs_log2fc_threshold <- 0  # all FDR < 0.05; choose 0.5 for effect-size filtering
maximum_number_of_genes <- NULL  # 5 for quick test; NULL = all
manual_target_genes <- NULL       # e.g. c("Neto2", "Aak1")
min_spots_per_sample_cluster <- 20L  # as in older barplot; change to 0 for all
normalization_scale_factor <- 10000
plot_width <- 7.5
plot_height <- 6.3
write_individual_pdfs <- TRUE
write_combined_pdf <- TRUE

# ======================== PATHS =============================================
stats_root <- file.path(project_root,"results",dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name,"01_statistics")
crescendo_root <- file.path(stats_root,"01_crescendoBatchCorrection")
input_rds <- file.path(crescendo_root,"07_objects",paste0(dataset_name,"_",clustering_name,
  "_crescendoBatchCorrection_pseudobulkCounts_beforeEdgeR.rds"))
qc_file <- file.path(crescendo_root,"06_qc",paste0(dataset_name,"_",clustering_name,
  "_crescendoBatchCorrection_spotLibrarySize_beforeAfter.tsv.gz"))
edgeR_file <- file.path(stats_root,"02_edgeR","03_edgeRResults",paste0(dataset_name,"_",
  clustering_name,"_crescendoCorrectedPseudobulk_edgeR_twoFactorInteraction_edgeRResults.RData"))
output_dir <- file.path(stats_root,"15_Crescendo_barplots_Neurotypical_vs_ASD")
plot_dir <- file.path(output_dir,"individual_pdf")

for (f in c(input_rds,qc_file,edgeR_file)) if (!file.exists(f)) stop("Missing input: ",f)
if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2: install.packages('ggplot2')")
dir.create(plot_dir,recursive=TRUE,showWarnings=FALSE)

read_tsv <- function(f) utils::read.delim(f,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
write_tsv <- function(x,f) utils::write.table(x,f,sep="\t",row.names=FALSE,quote=FALSE,na="NA")
sanitize <- function(x) gsub("[^A-Za-z0-9._-]+","_",as.character(x))
fmt_p <- function(x) if (!is.finite(x)) "NA" else if (x<0.001) formatC(x,digits=2,format="e") else sprintf("%.3f",x)
normalize_sex <- function(x) {z<-tolower(trimws(as.character(x)));ifelse(z %in% c("male","m"),"Male",ifelse(z %in% c("female","f"),"Female",NA_character_))}
normalize_group <- function(x) {z<-tolower(trimws(as.character(x)));ifelse(grepl("asd|autism",z),"ASD",ifelse(grepl("neurotypical|control|typical|^nt$",z),"Neurotypical",NA_character_))}

# ======================== edgeR DEG SELECTION ================================
e <- new.env(parent=globalenv())
load(edgeR_file,envir=e)
if (!exists("edgeR_perCluster_combinedResults",envir=e,inherits=FALSE)) stop("No edgeR_perCluster_combinedResults in edgeR file")
tests <- e$edgeR_perCluster_combinedResults
main_test <- "Overall_Group_ASD_vs_Neurotypical"
if (!main_test %in% names(tests)) stop("Missing main-effect comparison")
res <- as.data.frame(tests[[main_test]])
needed <- c("cluster_id","ensembl_gene_id","gene","logFC","PValue","FDR")
if (!all(needed %in% names(res))) stop("Missing edgeR columns: ",paste(setdiff(needed,names(res)),collapse=", "))
res$cluster_id <- as.character(res$cluster_id)
res$ensembl_gene_id <- as.character(res$ensembl_gene_id)
res$gene <- as.character(res$gene)
res$gene[is.na(res$gene)|res$gene==""] <- res$ensembl_gene_id[is.na(res$gene)|res$gene==""]
res <- res[is.finite(res$FDR)&is.finite(res$logFC)&res$FDR<fdr_threshold&abs(res$logFC)>abs_log2fc_threshold,,drop=FALSE]
if (!is.null(manual_target_genes)) res <- res[res$gene %in% manual_target_genes|res$ensembl_gene_id %in% manual_target_genes,,drop=FALSE]
res <- res[order(res$FDR,-abs(res$logFC)),,drop=FALSE]
if (!is.null(maximum_number_of_genes)) {
  ids <- head(unique(res$ensembl_gene_id),as.integer(maximum_number_of_genes))
  res <- res[res$ensembl_gene_id %in% ids,,drop=FALSE]
}
res <- res[!duplicated(paste(res$cluster_id,res$ensembl_gene_id)),,drop=FALSE]
if (!nrow(res)) stop("No DEG meets selection criteria")
message("Selected ",nrow(res)," cluster-gene results, ",length(unique(res$ensembl_gene_id))," unique genes")

# ======================== CRESCENDO RDS ======================================
message("Reading corrected pseudobulk matrices...")
obj <- readRDS(input_rds)
required_obj <- c("sample_metadata","gene_metadata","crescendo_corrected_counts_by_cluster")
if (!all(required_obj %in% names(obj))) stop("Missing Crescendo RDS components: ",paste(setdiff(required_obj,names(obj)),collapse=", "))
meta <- as.data.frame(obj$sample_metadata)
meta$sample_ID <- as.character(meta$sample_ID)
meta$sex_std <- normalize_sex(meta$sex)
meta$group_std <- normalize_group(meta$fmt_donor_group)
meta$group <- meta$group_std
if (anyNA(meta$sex_std)|anyNA(meta$group_std)|anyDuplicated(meta$sample_ID)) stop("Invalid biological sample metadata")
groups <- c("Neurotypical","ASD")
meta$group <- factor(meta$group,levels=groups)
if (nrow(meta)!=16L) warning("Expected 16 samples; got ",nrow(meta))

# ======================== RAW UMI DENOMINATORS ===============================
# QC is calculated from ALL raw RNA genes, NOT only the Crescendo gene subset.
# Sum raw spot UMI per sample x cluster; also count spots for the legacy >=20 filter.
message("Reading raw whole-transcriptome spot library sizes from Crescendo QC...")
qc <- read_tsv(qc_file)
qc_required <- c("spot_id","sample_ID","cluster_id","raw_total_UMI_allGenes")
if (!all(qc_required %in% colnames(qc))) stop("QC missing columns: ",paste(setdiff(qc_required,names(qc)),collapse=", "))
qc$sample_ID <- as.character(qc$sample_ID)
qc$cluster_id <- as.character(qc$cluster_id)
qc$raw_total_UMI_allGenes <- as.numeric(qc$raw_total_UMI_allGenes)
if (anyDuplicated(qc$spot_id)) stop("Duplicate spots in QC file")
if (any(!is.finite(qc$raw_total_UMI_allGenes))|any(qc$raw_total_UMI_allGenes<0)) stop("Invalid raw UMI values")
qc <- qc[qc$sample_ID %in% meta$sample_ID,,drop=FALSE]
key <- paste(qc$sample_ID,qc$cluster_id,sep="__C__")
raw_umi <- tapply(qc$raw_total_UMI_allGenes,key,sum)
n_spots <- table(key)
# Denominator represents exactly the raw UMI from spots belonging to the sample-cluster.

# ======================== PER-GENE EXTRACTION ================================
# Gene metadata relates Ensembl gene IDs to rownames(feature_id) and symbols.
gm <- as.data.frame(obj$gene_metadata)
message("Gene annotation fields: ",paste(names(gm),collapse=", "))
if (!"feature_id" %in% names(gm)) stop("Crescendo gene_metadata lacks feature_id")
if (!"gene" %in% names(gm)) stop("Crescendo gene_metadata lacks gene")
matrices <- obj$crescendo_corrected_counts_by_cluster
records <- list(); summaries <- list(); statuses <- list(); completed <- 0L
combined_file <- file.path(output_dir,"ALL_Crescendo_NT_vs_ASD_barplots.pdf")
if (write_combined_pdf) grDevices::pdf(combined_file,width=plot_width,height=plot_height,onefile=TRUE,useDingbats=FALSE)

# We preserve all individual sample values and never use spot-level observations
# as independent replicates for the barplot.
for (i in seq_len(nrow(res))) {
  r <- res[i,,drop=FALSE]; cl <- as.character(r$cluster_id)
  status <- tryCatch({
    mat <- matrices[[cl]]
    if (is.null(mat)) stop("corrected pseudobulk cluster not found")
    candidates <- unique(as.character(gm$feature_id[as.character(gm$gene)==r$gene]))
    ens_matches <- as.character(gm$feature_id[as.character(gm$ensembl_gene_id)==r$ensembl_gene_id])
    candidates <- unique(c(ens_matches,candidates,r$ensembl_gene_id,r$gene))
    candidates <- candidates[!is.na(candidates)&nzchar(candidates)]
    idx <- which(rownames(mat) %in% candidates)
    if (length(idx)!=1L) stop("Gene not uniquely mapped in corrected counts (",length(idx)," matches)")
    if (!all(meta$sample_ID %in% colnames(mat))) stop("Not all biological samples present in matrix")
    corrected <- as.numeric(mat[idx,meta$sample_ID,drop=TRUE])
    denom_key <- paste(meta$sample_ID,cl,sep="__C__")
    denominator <- as.numeric(raw_umi[denom_key]); spots <- as.integer(n_spots[denom_key]); spots[is.na(spots)] <- 0L
    if (any(!is.finite(corrected))|any(corrected<0)) stop("Nonfinite/negative corrected counts")
    if (any(is.na(denominator)&spots>0L)) stop("Missing raw UMI denominator for existing spots")
    keep <- spots>=min_spots_per_sample_cluster & is.finite(denominator)&denominator>0
    if (!any(keep)) stop("No samples pass min_spots_per_sample_cluster")
    d <- data.frame(gene=r$gene,ensembl_gene_id=r$ensembl_gene_id,cluster_id=cl,
      sample_ID=meta$sample_ID,group=meta$group,sex=meta$sex_std,
      donor_group=meta$group_std,n_spots=spots,corrected_gene_sum=corrected,
      raw_total_UMI_allGenes=denominator,included=keep,
      expression_log1p=ifelse(keep,log1p(corrected/pmax(denominator,1)*normalization_scale_factor),NA_real_),
      main_log2FC=r$logFC,main_PValue=r$PValue,main_FDR=r$FDR,stringsAsFactors=FALSE)
    records[[length(records)+1L]] <- d
    dd <- d[d$included,,drop=FALSE]
    s <- do.call(rbind,lapply(groups,function(gr){
      v <- dd$expression_log1p[dd$group==gr]
      data.frame(group=gr,n=length(v),mean=if(length(v))mean(v) else NA_real_,
        SD=if(length(v)>1L)sd(v) else NA_real_,
        SE=if(length(v)>1L)sd(v)/sqrt(length(v)) else NA_real_)
    }))
    s$group <- factor(s$group,levels=groups)
    s$edgeR_FDR <- r$FDR
    s$gene <- r$gene;s$ensembl_gene_id<-r$ensembl_gene_id;s$cluster_id<-cl
    summaries[[length(summaries)+1L]] <- s
    # Stars are based on BH-adjusted FDR from the original sex-adjusted edgeR main effect.
    star <- if (r$FDR < 1e-4) "****" else if (r$FDR < 1e-3) "***" else if (r$FDR < 0.01) "**" else if (r$FDR < 0.05) "*" else "ns"
    upper <- max(c(dd$expression_log1p, s$mean + ifelse(is.na(s$SE), 0, s$SE)), na.rm=TRUE)
    y_bracket <- upper + max(0.05, upper * 0.07)
    tick_height <- max(0.02, upper * 0.018)
    star_y <- y_bracket + max(0.045, upper * 0.035)
    top_y <- star_y + max(0.12, upper * 0.065)
    subtitle <- paste0("Sex-adjusted edgeR main group effect: log2FC=",sprintf("%.3f",r$logFC),
      " | nominal P=",fmt_p(r$PValue)," | FDR=",fmt_p(r$FDR),
      "\nBars: mean | whiskers: mean +/- SE | dots: biological samples | >=",
      min_spots_per_sample_cluster," spots/cluster")
    p <- ggplot2::ggplot() +
      ggplot2::geom_col(data=s,ggplot2::aes(x=group,y=mean),width=0.7,
        fill="white",color="black",linewidth=0.95,na.rm=TRUE)+
      ggplot2::geom_errorbar(data=s,ggplot2::aes(x=group,ymin=mean-SE,ymax=mean+SE),
        width=0.17,linewidth=0.7,na.rm=TRUE)+
      ggplot2::geom_point(data=dd,ggplot2::aes(x=group,y=expression_log1p),
        position=ggplot2::position_jitter(width=0.08,height=0,seed=42),
        shape=21,fill="white",color="black",size=2.6,stroke=0.65)+
      ggplot2::annotate("segment", x=1, xend=1, y=y_bracket-tick_height, yend=y_bracket, linewidth=0.65)+
      ggplot2::annotate("segment", x=1, xend=2, y=y_bracket, yend=y_bracket, linewidth=0.65)+
      ggplot2::annotate("segment", x=2, xend=2, y=y_bracket, yend=y_bracket-tick_height, linewidth=0.65)+
      ggplot2::annotate("text", x=1.5, y=star_y, label=star, size=5.5, fontface="bold")+
      ggplot2::scale_x_discrete(drop=FALSE)+
      ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,0)))+
      ggplot2::coord_cartesian(ylim=c(0,top_y),clip="off")+
      ggplot2::labs(title=paste0(r$gene," | cluster C",cl," | ASD vs Neurotypical"),
        subtitle=subtitle,x=NULL,
        y="log(1 + corrected gene sum / raw cluster UMI x 10,000)")+
      ggplot2::theme_classic(base_size=12)+
      ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=14,hjust=0.5),
        plot.subtitle=ggplot2::element_text(size=9.5,hjust=0.5),
        axis.text.x=ggplot2::element_text(size=10,color="black"),
        axis.title.y=ggplot2::element_text(size=11,face="bold"),
        axis.line=ggplot2::element_line(linewidth=0.7),
        plot.margin=ggplot2::margin(12,15,15,12))
    keyfile <- paste0(sanitize(r$gene),"_",sanitize(r$ensembl_gene_id),"_C",sanitize(cl),".pdf")
    if (write_individual_pdfs) ggplot2::ggsave(file.path(plot_dir,keyfile),plot=p,
      width=plot_width,height=plot_height,device="pdf",useDingbats=FALSE)
    if (write_combined_pdf) print(p)
    completed <- completed+1L
    message("[",i,"/",nrow(res),"] ",r$gene," C",cl," | ",nrow(dd)," samples")
    "OK"
  },error=function(err){warning(r$gene," C",cl,": ",conditionMessage(err)); paste0("ERROR: ",conditionMessage(err))})
  statuses[[length(statuses)+1L]] <- data.frame(gene=r$gene,cluster_id=cl,status=status)
}
if (write_combined_pdf) grDevices::dev.off()
if (length(records)) write_tsv(do.call(rbind,records),file.path(output_dir,"per_sample_expression_all_including_excluded.tsv"))
if (length(summaries)) write_tsv(do.call(rbind,summaries),file.path(output_dir,"group_mean_SE_SD.tsv"))
write_tsv(do.call(rbind,statuses),file.path(output_dir,"run_status.tsv"))
write_tsv(res,file.path(output_dir,"selected_edgeR_significant_genes.tsv"))
message("Done: ",completed," / ",nrow(res)," plots. Output: ",output_dir)
if (!completed) stop("No plots created; see run_status.tsv")