#!/usr/bin/env Rscript
# Section 3 supplement: sensitivity of "% transcriptional output lost" to a
# per-feature CPM-difference floor, at isoform and gene level. A noise control:
# only count a feature's SMG1i-DMSO increase toward "loss" if it clears t CPM.
# If the estimate were noise (many tiny ups), it would collapse toward 0 as t rises.
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tibble); library(tidyr)
  library(ggplot2); library(edgeR); library(patchwork)
})
select<-dplyr::select; filter<-dplyr::filter
options(bitmapType="cairo")
DATA <- "nmd_fig_data"
HERE <- "supplement_figures"; dir.create(HERE, recursive=TRUE, showWarnings=FALSE)
P <- function(f) file.path(DATA,f)
CT<-c("AT2","LAE","FB","MV"); CODE<-c(AT2="at2",LAE="dd",FB="fb",MV="mv")
CTCOL <- c(AT2="#F28E84", FB="#2CC6C9", MV="#C8A2EA", LAE="#9DBE2E")
THR <- c(0,0.25,0.5,1,2,3,5,7.5,10)

dge<-readRDS(P("dge_isoform_longread_2026.3.3.rds"))
samp<-as.data.frame(dge$samples)%>%rownames_to_column("bam")%>%
  mutate(ct=as.character(ct),treatment=as.character(treatment),donor=as.character(id))
samp$ct[samp$ct=="AT"]<-"AT2"; samp$ct[samp$ct=="DD"]<-"LAE"; samp<-samp[samp$ct%in%CT,]
cpm_mat<-cpm(dge[,samp$bam],normalized.lib.sizes=FALSE,log=FALSE)
i2g<-dge$genes%>%select(any_of(c("txid","gene_id")))
paired<-samp%>%group_by(ct,donor)%>%summarise(s=any(treatment=="Smg1i"),d=any(treatment=="DMSO"),.groups="drop")%>%filter(s&d)

keep<-intersect(rownames(cpm_mat),i2g$txid[!is.na(i2g$gene_id)&i2g$gene_id!=""])
gid<-i2g$gene_id[match(keep,i2g$txid)]; cpm_gene<-rowsum(cpm_mat[keep,samp$bam],group=gid)

# unrestricted all-positive-delta metric (as reported in the manuscript), swept over t
sweep<-function(M){bind_rows(lapply(CT,function(ct){
  dons<-paired$donor[paired$ct==ct]; num<-setNames(numeric(length(THR)),as.character(THR)); den<-0
  for(d in dons){
    sid<-samp$bam[samp$ct==ct&samp$donor==d&samp$treatment=="Smg1i"][1]
    did<-samp$bam[samp$ct==ct&samp$donor==d&samp$treatment=="DMSO"][1]
    delta<-M[,sid]-M[,did]; den<-den+sum(M[,sid])
    for(i in seq_along(THR)) num[i]<-num[i]+sum(delta[delta>=THR[i]])
  }
  data.frame(ct=ct,threshold=THR,pct=100*num/den)
}))}
iso<-sweep(cpm_mat)   %>% mutate(level="Isoform level")
gen<-sweep(cpm_gene)  %>% mutate(level="Gene level")
dat<-bind_rows(iso,gen)
dat$level<-factor(dat$level,levels=c("Isoform level","Gene level"))
dat$ct<-factor(dat$ct,levels=CT)
fwrite(dat,file.path(HERE,"data_output_lost_threshold.csv"))

# ---- compositional deflation: among expressed (DMSO CPM>=1) non-NMD isoforms,
#      fraction with negative vs positive SMG1i-DMSO delta (donor-pooled) ----
nmd<-lapply(CT,function(ct){f<-list.files(DATA,pattern=sprintf("nmd_mashr_die_%s_.*csv",CODE[[ct]]),full.names=TRUE)[1]
  d<-fread(f); d$txid[d$nmd_responsive%in%c(TRUE,"TRUE")]}); names(nmd)<-CT
defl<-bind_rows(lapply(CT,function(ct){
  dons<-paired$donor[paired$ct==ct]; up<-0; dn<-0; tot<-0; nmdset<-nmd[[ct]]
  for(d in dons){
    sid<-samp$bam[samp$ct==ct&samp$donor==d&samp$treatment=="Smg1i"][1]
    did<-samp$bam[samp$ct==ct&samp$donor==d&samp$treatment=="DMSO"][1]
    delta<-cpm_mat[,sid]-cpm_mat[,did]; keepx<-cpm_mat[,did]>=1 & !(names(delta)%in%nmdset)
    up<-up+sum(delta[keepx]>0); dn<-dn+sum(delta[keepx]<0); tot<-tot+sum(keepx)
  }
  data.frame(ct=ct, Decreased=100*dn/tot, Increased=100*up/tot)
}))
fwrite(defl,file.path(HERE,"data_nonnmd_deflation.csv"))
defl_l<-pivot_longer(defl,c(Decreased,Increased),names_to="dir",values_to="pct")
defl_l$ct<-factor(defl_l$ct,levels=CT)
defl_l$dir<-factor(defl_l$dir,levels=c("Increased","Decreased"))

pA<-ggplot(dat,aes(threshold,pct,colour=ct))+
  geom_line(linewidth=0.9)+geom_point(size=1.8)+
  facet_wrap(~level,nrow=1)+
  scale_colour_manual(values=CTCOL,name="Cell type")+
  scale_x_continuous(breaks=c(0,1,2,3,5,7.5,10))+
  scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(0,0.05)))+
  labs(x="Per-feature CPM-difference floor (SMG1i - DMSO >= t)",
       y="% transcriptional output lost")+
  theme_bw(base_size=13)+
  theme(strip.background=element_rect(fill="grey92"),
        panel.grid.minor=element_blank(),legend.position="top",
        axis.title=element_text(size=13))

pB<-ggplot(defl_l,aes(ct,pct,fill=dir))+
  geom_col(width=0.7)+
  geom_hline(yintercept=50,linetype="dashed",colour="grey40")+
  scale_fill_manual(values=c(Increased="#C0392B",Decreased="#4C78B0"),name=NULL)+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25),expand=expansion(mult=c(0,0)))+
  labs(x="Cell type", y="Non-NMD isoforms (%)",
       subtitle="Expressed non-NMD isoforms: SMG1i vs DMSO direction (compositional deflation)")+
  theme_bw(base_size=13)+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),
        legend.position="top",plot.subtitle=element_text(size=10.5),
        axis.title=element_text(size=13))

comb<-(pA/pB)+patchwork::plot_layout(heights=c(1,0.9))+
  patchwork::plot_annotation(tag_levels="A")&
  ggplot2::theme(plot.tag=ggplot2::element_text(face="bold",size=15))
ggsave(file.path(HERE,"figureS3_output_lost_threshold.png"),comb,width=9,height=8,dpi=300,bg="white")
ggsave(file.path(HERE,"figureS3_output_lost_threshold.pdf"),comb,width=9,height=8,bg="white")
cat("Wrote figure. Threshold sweep:\n")
print(as.data.frame(pivot_wider(dat%>%mutate(pct=round(pct,1)),names_from=threshold,values_from=pct,names_prefix="t>=")))
cat("\nNon-NMD deflation split:\n"); print(defl)
