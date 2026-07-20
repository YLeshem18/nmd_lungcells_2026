#!/usr/bin/env Rscript
# Three example-gene rows: isoform STRUCTURE (exons colored by NMD class,
# intron chevrons) + donor-paired log2FC CI. SHMT2/MV, SRSF2/MV, PCNA/LAE.
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork); library(edgeR); library(Isopair)
})
select<-dplyr::select; filter<-dplyr::filter
options(bitmapType="cairo")
DATA<-"nmd_fig_data"; GEX<-"gene_ex"; OUT<-"fig_panels"; LFSR<-0.05
GTF<-file.path(GEX,"sub.isoforms.gtf")
PAL<-c(Unproductive="#c0392b", Productive="#16a085")
EX <- tibble::tribble(~G,~CT,~sub,
  "SHMT2","MV","ISR/ATF4 induction overrides NMD: unproductive isoform surges, productive output still rises",
  "SRSF2","MV","Splicing-factor autoregulation: unproductive isoform surges, productive output falls",
  "PCNA","LAE","Transcriptional repression: whole-gene output falls, composition preserved")

pm<-read.csv(file.path(DATA,"mashr_isoform_posterior_means_2026.3.10.csv"),check.names=FALSE)
lf<-read.csv(file.path(DATA,"mashr_isoform_lfsr_2026.3.10.csv"),check.names=FALSE)
colnames(pm)[1]<-"isoform_id"; colnames(lf)[1]<-"isoform_id"
dge<-readRDS(file.path(DATA,"dge_isoform_longread_filtered_2026.3.3.rds"))
i2g<-data.frame(t=dge$genes$txid, sym=dge$genes$hgnc_symbol, stringsAsFactors=FALSE)
sa<-as.data.frame(dge$samples); sa$ct<-as.character(sa$ct); sa$ct[sa$ct=="AT"]<-"AT2"; sa$ct[sa$ct=="DD"]<-"LAE"
sa$treatment<-as.character(sa$treatment); sa$donor<-as.character(sa$id); sa$bam<-rownames(sa)
resolve<-function(df,ct) colnames(df)[grepl(paste0("_in_",ct,"$"),colnames(df))][1]
structures<-Isopair::parseIsoformStructures(GTF, verbose=FALSE)

exons_of<-function(iso){ r<-structures[structures$isoform_id==iso,]
  data.frame(s=r$exon_starts[[1]], e=r$exon_ends[[1]]) }
strand_of<-function(iso) structures$strand[match(iso,structures$isoform_id)]

build_structure <- function(iso, cls, gene, strand){
  yc<-setNames(seq(length(iso),1L), iso)          # unproductive already ordered on top
  regs<-bind_rows(lapply(iso,function(t){ ex<-exons_of(t)
    data.frame(isoform=t, xmin=ex$s, xmax=ex$e, y=yc[t]) }))
  regs$cls<-cls[regs$isoform]; H<-0.30
  regs$ymin<-regs$y-H; regs$ymax<-regs$y+H
  intr<-bind_rows(lapply(iso,function(t){ ex<-exons_of(t); ex<-ex[order(ex$s),]
    if(nrow(ex)<2) return(NULL); data.frame(isoform=t, x=head(ex$e,-1), xend=tail(ex$s,-1), y=yc[t]) }))
  intr<-intr[intr$xend>intr$x,,drop=FALSE]
  chev<-bind_rows(lapply(seq_len(nrow(intr)), function(i){
    a<-intr$x[i]; b<-intr$xend[i]; len<-b-a; if(len<300) return(NULL)
    n<-max(1,floor(len/1500)); cx<-a+(len/(n+1))*seq_len(n)
    data.frame(x=cx, y=intr$y[i], lab=if(strand=="-") "‹" else "›") }))
  xr<-range(c(regs$xmin,regs$xmax)); pad<-diff(xr)*0.02
  labs_df<-data.frame(x=xr[1]-pad, y=yc[iso], isoform=iso, cls=cls[iso])
  ggplot()+
    { if(nrow(intr)) geom_segment(data=intr,aes(x=x,xend=xend,y=y,yend=y),colour="grey60",linewidth=0.35) }+
    { if(!is.null(chev)&&nrow(chev)) geom_text(data=chev,aes(x=x,y=y,label=lab),colour="grey55",size=3) }+
    geom_rect(data=regs,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=cls),colour="black",linewidth=0.25)+
    geom_text(data=labs_df,aes(x=x,y=y,label=isoform,colour=cls),hjust=1,size=2.6,fontface="bold")+
    scale_fill_manual(values=PAL,breaks=c("Unproductive","Productive"),name=NULL)+
    scale_colour_manual(values=PAL,guide="none")+
    scale_y_continuous(breaks=NULL,limits=c(0.3,length(iso)+0.7))+
    coord_cartesian(xlim=c(xr[1]-diff(xr)*0.35, xr[2]+pad),clip="off")+
    labs(x=sprintf("Genomic position (%s strand)",strand),y=NULL)+
    theme_classic(base_size=12)+
    theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),
          axis.line=element_line(linewidth=0.4),axis.title.x=element_text(size=12),
          legend.position="bottom",legend.text=element_text(size=10),
          plot.margin=margin(4,6,4,4))
}

build_ci <- function(iso, cls, CT){
  cs<-sa$bam[sa$ct==CT]; cpmm<-cpm(dge[,cs],log=FALSE); donors<-unique(sa$donor[sa$ct==CT])
  l<-lf[[resolve(lf,CT)]]; names(l)<-lf$isoform_id
  fc<-bind_rows(lapply(iso,function(t){
    v<-vapply(donors,function(d){ dc<-sa$bam[sa$ct==CT&sa$donor==d&sa$treatment=="DMSO"]; sc<-sa$bam[sa$ct==CT&sa$donor==d&sa$treatment=="Smg1i"]
      if(length(dc)!=1||length(sc)!=1) return(NA_real_); log2((cpmm[t,sc]+1)/(cpmm[t,dc]+1)) },numeric(1)); v<-v[is.finite(v)]
    data.frame(isoform=t,logFC=mean(v),se=ifelse(length(v)>1,sd(v)/sqrt(length(v)),NA)) }))
  fc$cls<-cls[fc$isoform]; fc$sig<-!is.na(l[fc$isoform])&l[fc$isoform]<LFSR
  fc$CI.L<-fc$logFC-1.96*fc$se; fc$CI.R<-fc$logFC+1.96*fc$se
  fc$isoform<-factor(fc$isoform,levels=rev(iso))   # rev -> unproductive (iso[1]) on TOP, matching structure
  ggplot(fc,aes(logFC,isoform))+
    geom_vline(xintercept=0,colour="grey60",linetype="dashed")+
    geom_pointrange(aes(xmin=CI.L,xmax=CI.R,colour=cls,shape=sig),linewidth=0.6,size=0.45)+
    scale_colour_manual(values=PAL,guide="none")+
    scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1), breaks=c(TRUE,FALSE),
                       labels=c("lfsr < 0.05","n.s."), name=NULL)+
    labs(x="Smg1i log2FC (95% CI)",y=NULL)+
    guides(shape=guide_legend(override.aes=list(colour="grey25", size=0.6)))+
    theme_classic(base_size=12)+
    theme(axis.line=element_line(linewidth=0.4),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
          legend.position="top", legend.text=element_text(size=9), legend.key.size=unit(0.35,"cm"),
          plot.margin=margin(4,8,4,2))
}

for(i in seq_len(nrow(EX))){
  G<-EX$G[i]; CT<-EX$CT[i]
  iso<-intersect(i2g$t[i2g$sym==G], structures$isoform_id)
  cs<-sa$bam[sa$ct==CT]; cpmm<-cpm(dge[,cs],log=FALSE)
  iso<-intersect(iso, rownames(cpmm)[rowMeans(cpmm)>=1])
  l<-lf[[resolve(lf,CT)]]; p<-pm[[resolve(pm,CT)]]; names(l)<-lf$isoform_id; names(p)<-pm$isoform_id
  cls<-ifelse(!is.na(l[iso])&!is.na(p[iso])&l[iso]<LFSR&p[iso]>0,"Unproductive","Productive"); names(cls)<-iso
  iso<-iso[order(factor(cls[iso],levels=c("Unproductive","Productive")))]  # unproductive first -> top
  strand<-strand_of(iso[1])
  pS<-build_structure(iso,cls,G,strand); pC<-build_ci(iso,cls,CT)
  row<-(pS | pC) + plot_layout(widths=c(2.4,1)) +
    plot_annotation(title=bquote(bold(italic(.(G))~"in"~.(CT))),
      theme=theme(plot.title=element_text(face="bold",size=14)))
  ggsave(file.path(OUT,sprintf("gene_row_%s.png",G)), row, width=12, height=3.4, dpi=300, bg="white")
  cat("wrote gene_row_",G,".png (",length(iso)," isoforms, ",sum(cls=="Unproductive")," unproductive)\n",sep="")
}
