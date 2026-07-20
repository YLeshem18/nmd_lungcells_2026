#!/usr/bin/env Rscript
# Panel A: % transcriptional output lost to NMD, isoform level (Fig 2A)
# Panel B: productive output, NMD vs non-NMD genes (custom log2FC)
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(scales); library(edgeR)
})
select<-dplyr::select; filter<-dplyr::filter; rename<-dplyr::rename
options(bitmapType="cairo")
DATA_DIR<-"nmd_fig_data"; OUT_DIR<-"fig_panels"; P<-function(f) file.path(DATA_DIR,f)
CELL_TYPES<-c("AT2","LAE","FB","MV"); LFSR_NMD_THR<-0.05
FLOOR_BOTH<-1; FLOOR_DMSO<-2; pass_floor<-function(td,ts)(td>=FLOOR_BOTH&ts>=FLOOR_BOTH)|(td>=FLOOR_DMSO)

THEME <- theme_classic(base_size=13) +
  theme(axis.title=element_text(size=13), axis.title.y=element_text(margin=margin(r=4)),
        axis.text=element_text(size=12),
        axis.line=element_line(linewidth=0.4),
        legend.text=element_text(size=11), legend.title=element_blank(),
        legend.key.size=unit(0.4,"cm"), legend.position="top",
        strip.background=element_blank(), strip.text=element_text(face="bold",size=13),
        panel.grid=element_blank(), plot.margin=margin(8,6,4,24))

# ============================ PANEL A ============================
cat("Panel A: % output lost (isoform level)...\n")
dge <- readRDS(P("dge_isoform_longread_2026.3.3.rds"))   # UNFILTERED isoforms
sa <- as.data.frame(dge$samples)
sa$ct <- as.character(sa$ct); sa$ct[sa$ct=="AT"]<-"AT2"; sa$ct[sa$ct=="DD"]<-"LAE"
sa$treatment<-as.character(sa$treatment); sa$donor<-as.character(sa$id); sa$bam<-rownames(sa)
sa <- sa[sa$ct %in% CELL_TYPES,]
cpm_mat <- cpm(dge[, sa$bam], normalized.lib.sizes=FALSE, log=FALSE)   # plain CPM, no TMM
glob <- lapply(CELL_TYPES, function(ct){
  ss <- sa[sa$ct==ct,]; don <- intersect(ss$donor[ss$treatment=="Smg1i"], ss$donor[ss$treatment=="DMSO"])
  sp<-0; tot<-0
  for(d in don){
    sid<-ss$bam[ss$donor==d & ss$treatment=="Smg1i"][1]; did<-ss$bam[ss$donor==d & ss$treatment=="DMSO"][1]
    delta<-cpm_mat[,sid]-cpm_mat[,did]; sp<-sp+sum(pmax(delta,0)); tot<-tot+sum(cpm_mat[,sid])
  }
  data.frame(ct=ct, pct=100*sp/tot, n_donor=length(don))
}) %>% bind_rows()
cat("  ", paste(sprintf("%s=%.2f%%",glob$ct,glob$pct), collapse="  "), "\n")
# horizontal bars, ordered high -> low (highest at top)
glob$ct <- factor(glob$ct, levels = glob$ct[order(glob$pct)])
pA <- ggplot(glob, aes(ct, pct, fill=ct)) +
  geom_col(width=0.72) +
  geom_text(aes(label=sprintf("%.2f%%", pct)), hjust=-0.15, size=4) +
  coord_flip(clip="off") +
  scale_fill_manual(values=c(AT2="#F28E84", FB="#2CC6C9", MV="#C8A2EA", LAE="#9DBE2E"), guide="none") +
  scale_y_continuous(limits=c(0,21), breaks=seq(0,20,5), expand=expansion(mult=c(0,0))) +
  labs(x="Cell type", y="% Output Lost To NMD") +
  theme_minimal(base_size=13) +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=12),
        panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major.x=element_line(colour="grey88"),
        plot.margin=margin(6,18,6,6), legend.position="none")
ggsave(file.path(OUT_DIR,"panelA_outputlost.png"), pA, width=6, height=4, dpi=300, bg="white")
cat("  wrote panelA_outputlost.png\n")

# ============================ PANEL B ============================
cat("Panel B: productive NMD vs non-NMD...\n")
pm_iso<-read.csv(P("mashr_isoform_posterior_means_2026.3.10.csv"),check.names=FALSE)
lfsr_iso<-read.csv(P("mashr_isoform_lfsr_2026.3.10.csv"),check.names=FALSE)
colnames(pm_iso)[1]<-"isoform_id"; colnames(lfsr_iso)[1]<-"isoform_id"
dgef<-readRDS(P("dge_isoform_longread_filtered_2026.3.3.rds"))
iso2gene<-data.frame(isoform_id=dgef$genes$txid, gene_id=dgef$genes$gene_id, stringsAsFactors=FALSE)
samples<-as.data.frame(dgef$samples)
samples$cell_type<-as.character(samples$ct); samples$cell_type[samples$cell_type=="AT"]<-"AT2"; samples$cell_type[samples$cell_type=="DD"]<-"LAE"
samples$treatment<-as.character(samples$treatment); samples$sample_id<-as.character(samples$id)
samples<-samples[samples$cell_type %in% CELL_TYPES,]
counts_iso<-dgef$counts[, rownames(samples)]
resolve<-function(df,ct){cn<-colnames(df); if(ct%in%cn) return(ct); h<-cn[grepl(paste0("_in_",ct,"$"),cn)]; if(length(h)==1) return(h); stop("no col ",ct)}
build_adj<-function(ct){
  l<-lfsr_iso[[resolve(lfsr_iso,ct)]]; p<-pm_iso[[resolve(pm_iso,ct)]]
  is_unprod<-!is.na(l)&!is.na(p)&l<LFSR_NMD_THR&p>0
  cls<-setNames(ifelse(is_unprod,"unprod","prod"),lfsr_iso$isoform_id)
  cs<-rownames(samples)[samples$cell_type==ct]; cnt<-counts_iso[,cs,drop=FALSE]
  trt<-setNames(samples[cs,"treatment"],cs); donor<-setNames(samples[cs,"sample_id"],cs)
  rawcpm<-sweep(cnt,2,colSums(cnt),"/")*1e6
  ids<-rownames(rawcpm); gene<-iso2gene$gene_id[match(ids,iso2gene$isoform_id)]; class_v<-cls[ids]
  anchor<-!is.na(gene)&!grepl("::",gene,fixed=TRUE)&!is.na(class_v)
  dmso<-cs[trt[cs]=="DMSO"]; smg<-cs[trt[cs]=="Smg1i"]; dmso_mean<-rowMeans(rawcpm[,dmso,drop=FALSE])
  prod_idx<-anchor&class_v=="prod"; Pdmso<-tapply(dmso_mean[prod_idx],gene[prod_idx],sum)
  iu<-which(anchor&class_v=="unprod"); pg<-Pdmso[gene[iu]]; hp<-!is.na(pg)&pg>0; share<-ifelse(hp,dmso_mean[iu]/pg,NA_real_)
  ia<-iu[hp]; ina<-iu[!hp]; adj<-rawcpm
  for(s in smg){col<-rawcpm[,s]; Pgs<-tapply(col[prod_idx],gene[prod_idx],sum); nc<-col
    nc[ia]<-Pgs[gene[ia]]*share[hp]; nc[ina]<-dmso_mean[ina]; nc[!is.finite(nc)]<-0; adj[,s]<-nc*(1e6/sum(nc))}
  list(rawcpm=rawcpm,adj=adj,gene=gene,class_v=class_v,anchor=anchor,trt=trt,donor=donor)
}
process_ct<-function(ct){
  o<-build_adj(ct); k<-o$anchor; fk<-paste(o$gene[k],o$class_v[k],sep="|||")
  rawf<-rowsum(o$rawcpm[k,,drop=FALSE],fk); adjf<-rowsum(o$adj[k,,drop=FALSE],fk)
  fm<-data.frame(feature=rownames(rawf)); fm$gene<-sub("\\|\\|\\|.*$","",fm$feature); fm$class<-sub("^.*\\|\\|\\|","",fm$feature)
  sm<-data.frame(sample=colnames(rawf),donor=o$donor[colnames(rawf)],trt=o$trt[colnames(rawf)])
  melt<-function(m,v){d<-as.data.frame(m);d$feature<-rownames(m);pivot_longer(d,-feature,names_to="sample",values_to=v)}
  long<-melt(adjf,"adj")%>%left_join(melt(rawf,"raw"),by=c("feature","sample"))%>%left_join(fm,by="feature")%>%left_join(sm,by="sample")
  dm<-long%>%filter(trt=="DMSO")%>%select(gene,class,donor,cpm=adj)%>%pivot_wider(names_from=class,values_from=cpm,values_fill=0)%>%rename(prodCPM_DMSO=prod,unprodCPM_DMSO=unprod)
  sg<-long%>%filter(trt=="Smg1i")%>%select(gene,class,donor,adj,raw)%>%pivot_wider(names_from=class,values_from=c(adj,raw),values_fill=0)%>%rename(prodCPM_SMG_adj=adj_prod,prodCPM_SMG_raw=raw_prod,unprodCPM_SMG=raw_unprod)
  nmd_g<-unique(o$gene[o$anchor&o$class_v=="unprod"])
  full_join(dm,sg,by=c("gene","donor"))%>%mutate(ct=ct,gene_category=ifelse(gene%in%nmd_g,"NMD","non-NMD"))
}
gs<-bind_rows(lapply(CELL_TYPES,process_ct))%>%rename(gene_id=gene)
base<-gs%>%filter(!is.na(prodCPM_DMSO),!is.na(prodCPM_SMG_raw))%>%
  mutate(across(c(prodCPM_DMSO,prodCPM_SMG_raw,prodCPM_SMG_adj,unprodCPM_DMSO,unprodCPM_SMG),~coalesce(.x,0)))%>%
  group_by(ct,gene_id,gene_category)%>%
  summarise(prodCPM_DMSO=mean(prodCPM_DMSO),prodCPM_SMG_adj=mean(prodCPM_SMG_adj),
            unprodCPM_DMSO=mean(unprodCPM_DMSO),unprodCPM_SMG=mean(unprodCPM_SMG),prodCPM_SMG_raw=mean(prodCPM_SMG_raw),.groups="drop")%>%
  mutate(totalCPM_DMSO=prodCPM_DMSO+unprodCPM_DMSO,totalCPM_SMG=prodCPM_SMG_raw+unprodCPM_SMG)%>%
  filter(pass_floor(totalCPM_DMSO,totalCPM_SMG))%>%
  mutate(custom=log2((prodCPM_SMG_adj+1)/(prodCPM_DMSO+1)))
base$ct<-factor(base$ct,levels=CELL_TYPES); base$gene_category<-factor(base$gene_category,levels=c("NMD","non-NMD"))
cat("  median custom per CT (NMD):", paste(sprintf("%s=%.2f",CELL_TYPES,
   sapply(CELL_TYPES,function(c) median(base$custom[base$ct==c&base$gene_category=="NMD"]))),collapse="  "),"\n")
pB <- ggplot(base, aes(gene_category, custom, fill=gene_category)) +
  geom_hline(yintercept=0, colour="grey60") +
  geom_boxplot(outlier.size=0.3, outlier.alpha=0.2, width=0.6, colour="grey25") +
  facet_wrap(~ct, nrow=1) + coord_cartesian(ylim=c(-1.2,1.2)) +
  scale_fill_manual(values=c("NMD"="#c0392b","non-NMD"="#7f8c8d")) +
  labs(x=NULL, y="Productive log2FC (Smg1i Vs DMSO)") +
  THEME + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(file.path(OUT_DIR,"panelB_nmd_vs_nonmd.png"), pB, width=6, height=4, dpi=300, bg="white")
cat("  wrote panelB_nmd_vs_nonmd.png\n")
