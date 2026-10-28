#!/usr/bin/env Rscript
# =============================================================================
# 13b_Crescendo_edgeR_GraphPadBarplots_STANDALONE.R
# GraphPad-style barplots for the edgeR main donor-group effect, using
# Crescendo-corrected pseudobulk expression and the ORIGINAL edgeR TMM factors.
#
#   - No new statistical tests: log2FC and FDR are read from the original edgeR.
#   - No rerun of Crescendo or edgeR. The saved DGEList contains TMM factors.
#   - 1 point = 1 biological sample in a specified anatomical cluster.
#   - TWO styles: combined samples (unmarked) and sex-marked dots.
#   - All output goes to a completely NEW, timestamped directory.
#   - Individual PDFs, multipage PDFs, per-sample TSV and Prism-ready CSVs.
#
# Default Y axis:
#   log2(1 + edgeR::cpm(original_TMM_DGEList, log = FALSE))
# Also exports conventional edgeR logCPM (prior.count=2), without using it
# for the default plot. Neither QN nor row Z-scores are used.
# =============================================================================
options(stringsAsFactors = FALSE)

# -------------------------- USER SETTINGS ------------------------------------
project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
dataset_name <- "maternalFMT_n16samples"
clustering_name <- "logNormalizeVstHvg2000_ccaDims20K20PruneSNN0067_leidenRes040"
fdr_threshold <- 0.05
min_abs_log2fc <- 0
manual_genes <- NULL          # e.g. c("Syngap1", "Neto2") or NULL for all
max_results <- NULL           # e.g. 5 for a quick test; NULL for all
# Choices: "log2CPMplus1" (default, nonnegative) or "edgeRlogCPM"
y_metric <- "log2CPMplus1"
plot_width <- 5.1
plot_height <- 5.7
include_stars <- TRUE
# Use a new directory for each run. NEVER write to original edgeR folders.
run_label <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ------------------------------ PATHS ----------------------------------------
stats_root <- file.path(project_root, "results", dataset_name,
  "pseudobulkPerCluster_geneCounts_crescendoEdgeR_twoFactorInteraction",
  clustering_name, "01_statistics")
analysis_prefix <- paste(dataset_name, clustering_name,
  "crescendoCorrectedPseudobulk_edgeR_twoFactorInteraction", sep = "_")
original_model_file <- file.path(stats_root, "02_edgeR", "05_objects",
  paste0(analysis_prefix, "_edgeRModelObjects.RData"))
original_results_file <- file.path(stats_root, "02_edgeR", "03_edgeRResults",
  paste0(analysis_prefix, "_edgeRResults.RData"))
output_root <- file.path(stats_root, "16_GraphPad_donorGroupBarplots", paste0("run_", run_label))
out_data <- file.path(output_root, "01_data")
out_prism <- file.path(output_root, "02_GraphPad_import_CSV")
out_plain <- file.path(output_root, "03_barplots_combined")
out_sex <- file.path(output_root, "04_barplots_sexMarked")
out_multipage <- file.path(output_root, "05_allBarplots_PDF")

# -------------------------- VALIDATION ---------------------------------------
for (f in c(original_model_file, original_results_file)) {
  if (!file.exists(f)) stop("Missing ORIGINAL file (no data changed): ", f)
}
for (pkg in c("edgeR", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing R package: ", pkg, ". Install it before running this script.")
  }
}
if (!y_metric %in% c("log2CPMplus1", "edgeRlogCPM")) stop("Invalid y_metric")
if (dir.exists(output_root)) stop("Output already exists; refusing to overwrite: ", output_root,
  "\nSet run_label to a unique value.")
for (d in c(out_data, out_prism, out_plain, out_sex, out_multipage)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

write_tsv <- function(x, p) utils::write.table(x, p, sep="\t", quote=FALSE,
  row.names=FALSE, col.names=TRUE, na="NA", fileEncoding="UTF-8")
clean_name <- function(x) gsub("[^a-zA-Z0-9._-]", "_", as.character(x))
fmt_p <- function(x) ifelse(is.finite(x), ifelse(x < .001,
  formatC(x, digits=2, format="e"), sprintf("%.3f", x)), "NA")
normalize_cluster <- function(x) sub("^[Cc]", "", as.character(x))
normalize_group <- function(x) {
  z <- tolower(trimws(as.character(x)))
  ifelse(grepl("asd|autism", z), "ASD",
    ifelse(grepl("neurotypical|control|typical|^nt$", z), "Neurotypical", NA_character_))
}
normalize_sex <- function(x) {
  z <- tolower(trimws(as.character(x)))
  ifelse(z %in% c("m", "male"), "Male",
    ifelse(z %in% c("f", "female"), "Female", NA_character_))
}

cluster_labels <- c(
  "1" = "Posterior and sensory relay thalamic nuclei",
  "2" = "Isocortex, layers 4/5",
  "3" = "Cortical layer 1 and hippocampal neuropil",
  "4" = "Fiber tracts",
  "5" = "Cortical subplate and deep olfactory areas",
  "6" = "Hypothalamus",
  "7" = "Isocortex, layer 6",
  "8" = "Isocortex, layers 2/3",
  "9" = "Reticular, ventral geniculate and habenular regions",
  "10" = "Striatum-like amygdalar nuclei",
  "11" = "Medial thalamic nuclei",
  "12" = "Caudoputamen",
  "13" = "Hippocampal CA pyramidal layers",
  "14" = "Meninges and vasculature",
  "15" = "Ventricles",
  "16" = "Dentate gyrus"
)

# -------------------------- LOAD ORIGINAL RESULTS ----------------------------
result_env <- new.env(parent=globalenv())
model_env <- new.env(parent=globalenv())
load(original_results_file, envir=result_env)
load(original_model_file, envir=model_env)
if (!exists("edgeR_perCluster_combinedResults", envir=result_env, inherits=FALSE))
  stop("No edgeR_perCluster_combinedResults in original results")
if (!exists("edgeR_perCluster_modelObjects", envir=model_env, inherits=FALSE))
  stop("No edgeR_perCluster_modelObjects in original model file")
results <- result_env$edgeR_perCluster_combinedResults
models <- model_env$edgeR_perCluster_modelObjects
main_name <- "Overall_Group_ASD_vs_Neurotypical"
if (!main_name %in% names(results)) stop("Missing main donor-group contrast")
main <- as.data.frame(results[[main_name]], stringsAsFactors=FALSE)
needed <- c("cluster_id", "ensembl_gene_id", "gene", "logFC", "FDR")
if (!all(needed %in% colnames(main))) stop("Missing columns: ",
  paste(setdiff(needed,colnames(main)),collapse=", "))
main$cluster_id <- normalize_cluster(main$cluster_id)
main$ensembl_gene_id <- as.character(main$ensembl_gene_id)
main$gene <- as.character(main$gene)
blank <- is.na(main$gene) | main$gene == ""
main$gene[blank] <- main$ensembl_gene_id[blank]
main <- main[is.finite(main$FDR) & main$FDR < fdr_threshold &
  is.finite(main$logFC) & abs(main$logFC) > min_abs_log2fc, , drop=FALSE]
if (!is.null(manual_genes)) main <- main[
  main$gene %in% manual_genes | main$ensembl_gene_id %in% manual_genes,
  , drop=FALSE]
main <- main[order(main$FDR, -abs(main$logFC)), , drop=FALSE]
main <- main[!duplicated(paste(main$cluster_id,main$ensembl_gene_id)),,drop=FALSE]
if (!is.null(max_results)) main <- head(main, as.integer(max_results))
if (!nrow(main)) stop("No gene-cluster results passed the filter")
message("Selected ", nrow(main), " gene-cluster results from original edgeR")
write_tsv(main, file.path(out_data, "01_selectedSignificantGeneClusterResults.tsv"))

# -------------------------- SAMPLE-LEVEL EXPRESSION --------------------------
per_sample <- list()
lib_info <- list()
status <- list()
i_out <- 0L
for (cl in unique(main$cluster_id)) {
  selected_cl <- main[main$cluster_id == cl, , drop=FALSE]
  if (!cl %in% names(models)) {
    for (j in seq_len(nrow(selected_cl))) status[[length(status)+1L]] <- data.frame(
      gene=selected_cl$gene[j], cluster_id=cl, status="missing_model")
    warning("Missing saved model for C", cl)
    next
  }
  dge <- models[[cl]]$dge
  if (is.null(dge) || is.null(dge$counts) || is.null(dge$samples))
    stop("Invalid DGEList for cluster ",cl)
  if (is.null(dge$samples$norm.factors) || is.null(dge$samples$lib.size))
    stop("Missing original TMM factors or library sizes for C",cl)

  raw_cpm <- edgeR::cpm(dge, log=FALSE, normalized.lib.sizes=TRUE)
  edgeR_log_cpm <- edgeR::cpm(dge, log=TRUE, prior.count=2,
    normalized.lib.sizes=TRUE)
  sample_ids <- colnames(dge$counts)
  md <- as.data.frame(dge$samples)
  if (!"fmt_donor_group" %in% names(md) || !"sex" %in% names(md))
    stop("Missing sample sex/donor-group metadata in model for C",cl)
  if (!is.null(md$sample_ID) && !identical(as.character(md$sample_ID),sample_ids))
    stop("Sample_ID order differs from DGEList counts for C",cl)
  md$sample_ID <- sample_ids
  md$group <- normalize_group(md$fmt_donor_group)
  md$sex_std <- normalize_sex(md$sex)
  if (anyNA(md$group) || anyNA(md$sex_std) || anyDuplicated(md$sample_ID))
    stop("Invalid sample metadata for C",cl)
  if (length(sample_ids) != 16L) warning("C",cl," has ",length(sample_ids)," samples (expected 16)")

  lib_info[[length(lib_info)+1L]] <- data.frame(
    cluster_id=cl, sample_ID=sample_ids,
    donor_group=md$group, sex=md$sex_std,
    original_library_size=as.numeric(md$lib.size),
    original_TMM_factor=as.numeric(md$norm.factors),
    effective_library_size=as.numeric(md$lib.size*md$norm.factors),
    stringsAsFactors=FALSE
  )

  for (j in seq_len(nrow(selected_cl))) {
    rr <- selected_cl[j,,drop=FALSE]
    ensembl <- as.character(rr$ensembl_gene_id)
    if (!ensembl %in% rownames(raw_cpm)) {
      status[[length(status)+1L]] <- data.frame(
        gene=rr$gene, cluster_id=cl, status="gene_not_in_filtered_model")
      next
    }
    current <- data.frame(
      gene=as.character(rr$gene),
      ensembl_gene_id=ensembl,
      cluster_id=cl,
      cluster_annotation=if (cl %in% names(cluster_labels)) unname(cluster_labels[[cl]]) else "",
      sample_ID=sample_ids,
      donor_group=md$group,
      sex=md$sex_std,
      TMM_CPM=as.numeric(raw_cpm[ensembl,sample_ids]),
      log2_CPM_plus1=log2(1+as.numeric(raw_cpm[ensembl,sample_ids])),
      edgeR_log2CPM=as.numeric(edgeR_log_cpm[ensembl,sample_ids]),
      main_log2FC=as.numeric(rr$logFC),
      main_FDR=as.numeric(rr$FDR),
      stringsAsFactors=FALSE
    )
    if (any(!is.finite(current$log2_CPM_plus1)) || any(!is.finite(current$edgeR_log2CPM)))
      stop("Non-finite expression for ",rr$gene," C",cl)
    i_out <- i_out + 1L
    per_sample[[i_out]] <- current
    status[[length(status)+1L]] <- data.frame(
      gene=rr$gene, cluster_id=cl, status="ok")
  }
}
if (!length(per_sample)) stop("No expressions could be extracted from original DGELists")
long <- do.call(rbind,per_sample)
rownames(long) <- NULL
long$plot_value <- if (y_metric=="log2CPMplus1") long$log2_CPM_plus1 else long$edgeR_log2CPM
write_tsv(long, file.path(out_data,"02_perSample_geneClusterExpression.tsv"))
write_tsv(do.call(rbind,lib_info), file.path(out_data,"03_originalTMM_libraryFactors.tsv"))
write_tsv(do.call(rbind,status), file.path(out_data,"04_extractionStatus.tsv"))

# -------------------------- SUMMARIES & GRAPHPAD TABLES ----------------------
make_group_summary <- function(d) {
  do.call(rbind, lapply(c("Neurotypical","ASD"), function(gr) {
    x <- d$plot_value[d$donor_group == gr]
    data.frame(donor_group=gr, n=length(x),
      mean=if(length(x)) mean(x) else NA_real_,
      SD=if(length(x)>1L) stats::sd(x) else NA_real_,
      SEM=if(length(x)>1L) stats::sd(x)/sqrt(length(x)) else NA_real_,
      stringsAsFactors=FALSE)
  }))
}
keys <- unique(paste(long$cluster_id,long$ensembl_gene_id,sep="__"))
plot_summaries <- list()
summary_rows <- list()
plot_records <- list()
for (key in keys) {
  split <- strsplit(key,"__",fixed=TRUE)[[1]]
  cl <- split[[1]]
  ens <- split[[2]]
  d <- long[long$cluster_id==cl & long$ensembl_gene_id==ens,,drop=FALSE]
  d <- d[order(match(d$donor_group,c("Neurotypical","ASD")),d$sex,d$sample_ID),,drop=FALSE]
  s <- make_group_summary(d)
  s$gene <- d$gene[[1]]
  s$ensembl_gene_id <- ens
  s$cluster_id <- cl
  s$main_log2FC <- d$main_log2FC[[1]]
  s$main_FDR <- d$main_FDR[[1]]
  summary_rows[[length(summary_rows)+1L]] <- s

  # GraphPad Prism column format: each column is a donor group; each row a replicate.
  nt <- d$plot_value[d$donor_group=="Neurotypical"]
  asd <- d$plot_value[d$donor_group=="ASD"]
  nmax <- max(length(nt),length(asd))
  pad <- function(x) c(x,rep(NA_real_,nmax-length(x)))
  stem <- paste0(clean_name(d$gene[[1]]),"_C",clean_name(cl))
  if (anyDuplicated(paste(long$gene,long$cluster_id,long$ensembl_gene_id))) {
    # Per-sample rows are expected to repeat. Real file-name conflicts are
    # handled by adding the gene's stable Ensembl identifier when necessary.
    duplicate_stem <- any(vapply(plot_records, function(z) identical(z$stem,stem), logical(1)))
    if (duplicate_stem) stem <- paste0(stem,"_",clean_name(ens))
  }
  utils::write.csv(data.frame(Neurotypical=pad(nt),ASD=pad(asd),check.names=FALSE),
    file=file.path(out_prism,paste0(stem,"_GraphPad.csv")),row.names=FALSE,na="")
  plot_records[[length(plot_records)+1L]] <- list(d=d,s=s,stem=stem)
}
write_tsv(do.call(rbind,summary_rows),file.path(out_data,"05_groupMean_SEM_SD.tsv"))

# -------------------------- GRAPH-PAD-LIKE PLOTTING --------------------------
# Bar: group mean; whisker: SEM; dot: one biological sample.
# The stars and FDR refer to the ORIGINAL main effect from the interaction model,
# not to an unadjusted t test on these plotted values.
make_plot <- function(d, s, sex_marked=FALSE) {
  labels <- c("Neurotypical","ASD")
  d$group_order <- match(d$donor_group,labels)
  d$sex <- factor(d$sex,levels=c("Male","Female"))
  s$group_order <- match(s$donor_group,labels)
  s$SEM[!is.finite(s$SEM)] <- 0
  s$lower <- s$mean - s$SEM
  s$upper <- s$mean + s$SEM
  # Fixed deterministic horizontal scatter; no statistical jitter in Y.
  d$xpos <- d$group_order
  for (gr in 1:2) {
    idx <- which(d$group_order == gr)
    if (length(idx)) d$xpos[idx] <- gr + seq(-0.11,0.11,length.out=length(idx))
  }
  ymx <- max(c(d$plot_value,s$upper,0),na.rm=TRUE)
  ymn <- min(c(d$plot_value,s$lower,0),na.rm=TRUE)
  yr <- max(0.4,ymx-ymn)
  gap <- 0.09 * yr
  bracket_y <- ymx + gap
  star_y <- bracket_y + 0.55*gap
  upper_lim <- star_y + 1.2*gap
  lower_lim <- if (ymn < 0) ymn-0.18*gap else 0
  fdr <- d$main_FDR[[1]]
  star <- if (fdr < 1e-4) "****" else if (fdr < 1e-3) "***" else
    if (fdr < 0.01) "**" else if (fdr < 0.05) "*" else "ns"
  cl <- d$cluster_id[[1]]
  region <- d$cluster_annotation[[1]]
  # Compact cluster label avoids oversized axes / titles.
  region <- if (nzchar(region)) strwrap(region, width=38) else ""
  title <- paste0(d$gene[[1]],"  |  C",cl)
  y_label <- if (y_metric=="log2CPMplus1")
    "log2(TMM-normalized CPM + 1)" else "edgeR log2 CPM (prior.count = 2)"
  sub <- paste0(paste(region,collapse="\n"),"\n",
    "edgeR main donor-group effect: log2FC = ", sprintf("%.3f",d$main_log2FC[[1]]),
    ", FDR = ",fmt_p(fdr))
  s$barfill <- ifelse(s$donor_group=="Neurotypical","#ECECEC","#A3A3A3")
  p <- ggplot2::ggplot() +
    ggplot2::geom_col(data=s, ggplot2::aes(x=group_order,y=mean),
      width=0.61,fill=s$barfill,color="black",linewidth=0.85) +
    ggplot2::geom_errorbar(data=s,ggplot2::aes(x=group_order,ymin=lower,ymax=upper),
      width=0.13,linewidth=0.85,color="black")
  if (sex_marked) {
    p <- p +
      ggplot2::geom_point(data=d,ggplot2::aes(x=xpos,y=plot_value,
        colour=sex,shape=sex),size=3.25,stroke=0.7) +
      ggplot2::scale_colour_manual(values=c(Male="#356EA6",Female="#C35683"),drop=FALSE) +
      ggplot2::scale_shape_manual(values=c(Male=16,Female=17),drop=FALSE) +
      ggplot2::guides(colour=ggplot2::guide_legend(title="Sample sex"),shape="none") +
      ggplot2::theme(legend.position="bottom",legend.text=ggplot2::element_text(size=10))
  } else {
    p <- p + ggplot2::geom_point(data=d,ggplot2::aes(x=xpos,y=plot_value),
      shape=21,fill="white",colour="black",size=3.2,stroke=0.75)
  }
  if (include_stars) {
    p <- p +
      ggplot2::annotate("segment",x=1,xend=2,y=bracket_y,yend=bracket_y,linewidth=0.8) +
      ggplot2::annotate("segment",x=1,xend=1,y=bracket_y,yend=bracket_y-0.26*gap,linewidth=0.8) +
      ggplot2::annotate("segment",x=2,xend=2,y=bracket_y,yend=bracket_y-0.26*gap,linewidth=0.8) +
      ggplot2::annotate("text",x=1.5,y=star_y,label=star,size=5.5,fontface="bold")
  }
  p + ggplot2::scale_x_continuous(breaks=1:2,labels=labels,limits=c(0.5,2.5),
      expand=ggplot2::expansion(mult=0)) +
    ggplot2::scale_y_continuous(limits=c(lower_lim,upper_lim),
      expand=ggplot2::expansion(mult=0)) +
    ggplot2::labs(title=title,subtitle=sub,x=NULL,y=y_label,
      caption="Bars: mean ± SEM  |  Each dot: one biological sample") +
    ggplot2::theme_classic(base_size=12) +
    ggplot2::theme(
      plot.title=ggplot2::element_text(face="bold",size=15,hjust=0.5,margin=ggplot2::margin(b=3)),
      plot.subtitle=ggplot2::element_text(size=10,hjust=0.5,lineheight=1.15,margin=ggplot2::margin(b=14)),
      plot.caption=ggplot2::element_text(size=9,hjust=0.5,margin=ggplot2::margin(t=12)),
      axis.title.y=ggplot2::element_text(size=11,face="bold",margin=ggplot2::margin(r=7)),
      axis.text.x=ggplot2::element_text(size=11,face="bold",colour="black"),
      axis.text.y=ggplot2::element_text(size=10,colour="black"),
      axis.line=ggplot2::element_line(linewidth=0.85,color="black"),
      axis.ticks=ggplot2::element_line(linewidth=0.75,color="black"),
      plot.margin=ggplot2::margin(t=15,r=18,b=12,l=12)
    )
}

# -------------------------- OUTPUT PLOTS -------------------------------------
combined_path <- file.path(out_multipage,"ALL_donorGroup_barplots_combined.pdf")
sex_path <- file.path(out_multipage,"ALL_donorGroup_barplots_sexMarked.pdf")
# Separate passes prevent changes to current PDF devices when saving individual plots.
for (record in plot_records) {
  ggplot2::ggsave(file.path(out_plain,paste0(record$stem,"_barplot.pdf")),
    plot=make_plot(record$d,record$s,sex_marked=FALSE),
    width=plot_width,height=plot_height,units="in",device="pdf",useDingbats=FALSE)
  ggplot2::ggsave(file.path(out_sex,paste0(record$stem,"_barplot_sexMarked.pdf")),
    plot=make_plot(record$d,record$s,sex_marked=TRUE),
    width=plot_width,height=plot_height,units="in",device="pdf",useDingbats=FALSE)
}
for (sex_marked in c(FALSE,TRUE)) {
  dest <- if (sex_marked) sex_path else combined_path
  grDevices::pdf(dest,width=plot_width,height=plot_height,onefile=TRUE,useDingbats=FALSE)
  for (record in plot_records) print(make_plot(record$d,record$s,sex_marked=sex_marked))
  grDevices::dev.off()
}
readme <- c(
  "GraphPad-style barplots, main donor-group effect (ASD vs Neurotypical)",
  paste("Run:",run_label),
  paste("Input edgeR statistics:",original_results_file),
  paste("Input original edgeR model:",original_model_file),
  "No original inputs are modified; no new differential-expression tests are run.",
  "Original edgeR test: Overall_Group_ASD_vs_Neurotypical, sex-averaged contrast from interaction model.",
  "TMM factors and library sizes come directly from the saved edgeR DGEList for each cluster.",
  paste("Plotted metric:", if(y_metric=="log2CPMplus1") "log2(1 + TMM-normalized CPM)" else "edgeR log2 CPM with prior.count=2"),
  "The long-format TSV exports both metrics, raw CPM and original sex/group metadata.",
  "CSV files in 02_GraphPad_import_CSV are wide, 2-column Prism-ready tables.",
  "Bars are unadjusted descriptive group means (mean +/- SEM); points are biological samples.",
  "Sex-marked version shows male blue circles and female pink triangles.",
  "Statistics displayed (log2FC, FDR and stars) come from original edgeR, NOT from a new t test on plot values.",
  "QN and per-gene Z-score are not applied to plotted expression.",
  paste("Selected gene-cluster rows:",length(plot_records))
)
writeLines(readme,file.path(output_root,"README.txt"))
message("SUCCESS: ",length(plot_records)," gene-cluster barplots in two styles")
message("NEW output only: ",output_root)