#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat);library(SeuratObject);library(Matrix);library(data.table)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2)stop('Usage: Rscript 01_export_visium.R <project_root> <output_folder>')
root <- normalizePath(args[1],mustWork=TRUE); out <- args[2]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
p <- file.path(root,'results','maternalFMT_n16samples','seurat_clustering_analysis','logNormalize_vst','hvg2000','02_cca','dims20','k20','prune0067','RData','02_maternalFMT_n16samples_logNormalizeVst_hvg2000_ccaDims20K20PruneSNN0067_res010to100by010_multiClusteringAndUmap.RData')
if(!file.exists(p))stop('Missing Seurat object: ',p)
e <- new.env(); nm <- load(p,envir=e); nm <- nm[vapply(nm,function(n)inherits(e[[n]],'Seurat'),logical(1))]; if(length(nm)<1) stop('No Seurat object in RData')
s <- e[[nm[which.max(vapply(nm,function(n)ncol(e[[n]]),numeric(1)))]]]; md <- s[[]]
stopifnot(all(c('sample_ID','leiden_res0.40') %in% colnames(md)))
samples <- c('1_1F','1_1Fd','1_1M','2_1M','2_1Md','3_1F','3_1M','5_1M','5_3F','12_1M','13_1F','13_1M','15_1F','18_1F','18_1M','23_1F')
for(id in samples){
 d <- file.path(out,id); dir.create(d,recursive=TRUE,showWarnings=FALSE)
 if(file.exists(file.path(d,'matrix.mtx.gz')) && file.exists(file.path(d,'spots.tsv'))) {message('SKIP ',id);next}
 cells <- rownames(md)[as.character(md$sample_ID)==id];if(!length(cells))stop('No spots for ',id)
 assay <- s[['RNA']]; lyr <- SeuratObject::Layers(assay,search='^counts'); if(!length(lyr))stop('Missing raw counts')
 matlist <- lapply(lyr,function(l){x<-SeuratObject::LayerData(assay,layer=l);k<-intersect(cells,colnames(x));if(length(k))x[,k,drop=FALSE] else NULL}); matlist<-Filter(Negate(is.null),matlist)
 if(!length(matlist))stop('No raw counts for ',id)
 if(length(matlist)>1 && !all(vapply(matlist,function(m)identical(rownames(m),rownames(matlist[[1]])),logical(1))))stop('Inconsistent gene indices across count layers')
 x<- if(length(matlist)==1)matlist[[1]] else do.call(cbind,matlist)
 if(anyDuplicated(colnames(x)))stop('Duplicated barcodes')
 if(!setequal(cells,colnames(x)))stop('Missing spots in raw counts for ',id)
 x <- methods::as(x[,cells,drop=FALSE],'dgCMatrix')
 if(any(x@x<0) || any(abs(x@x-round(x@x))>1e-6))stop('Expected nonnegative integer counts')
 Matrix::writeMM(x,file.path(d,'matrix.mtx'));R.utils::gzip(file.path(d,'matrix.mtx'),overwrite=TRUE)
 fwrite(data.table(gene=rownames(x)),file.path(d,'genes.tsv'),sep='\t')
 fwrite(data.table(spot=colnames(x),sample_ID=id,cluster=as.character(md[colnames(x),'leiden_res0.40'])),file.path(d,'spots.tsv'),sep='\t')
 message('EXPORTED ',id,' ',ncol(x),' spots')
}
