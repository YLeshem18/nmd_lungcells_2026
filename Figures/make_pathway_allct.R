#!/usr/bin/env Rscript
# Render the ER protein-processing KEGG pathway (hsa04141) colored by SR gene logFC,
# one diagram per cell type (AT2, LAE, FB, MV). Mirrors kegg_pathview.Rmd pathview() call.
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(org.Hs.eg.db); library(AnnotationDbi); library(pathview) })
DATA<-"nmd_fig_data"
WORK<-"fig_panels"   # cached hsa04141.xml/.png here if present; else pathview fetches from KEGG
setwd(WORK)
CT<-c("AT2","LAE","FB","MV"); CODE<-c(AT2="at2",LAE="dd",FB="fb",MV="mv"); KID<-"hsa04141"
sv<-function(x) sub("\\..*$","",x)

for(ct in CT){
  d<-fread(list.files(DATA,pattern=sprintf("nmd_mashr_dge_%s_.*csv$",CODE[[ct]]),full.names=TRUE)[1])
  d<-d[!is.na(d$logFC),]
  m<-suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=sv(unique(d$ensembl_gene_id_version)),
       columns="ENTREZID", keytype="ENSEMBL"))
  m<-m[!is.na(m$ENTREZID)&!duplicated(m$ENSEMBL),]
  d$ens<-sv(d$ensembl_gene_id_version); d$entrez<-m$ENTREZID[match(d$ens,m$ENSEMBL)]
  d<-d[!is.na(d$entrez),]; d<-d[order(-abs(d$logFC))]; d<-d[!duplicated(d$entrez)]
  gene_vec<-setNames(d$logFC, d$entrez)
  cat(sprintf("%s: %d entrez genes -> pathview %s\n", ct, length(gene_vec), KID))
  tryCatch(pathview(gene.data=gene_vec, pathway.id=KID, species="hsa", out.suffix=ct,
    limit=list(gene=2,cpd=1), low=list(gene="#2166AC"), mid=list(gene="#FFFFFF"),
    high=list(gene="#B2182B"), kegg.native=TRUE, same.layer=FALSE),
    error=function(e) cat("  ",ct,"error:",conditionMessage(e),"\n"))
}
cat("Done. hsa04141.<CT>.png written to fig_panels/.\n")
