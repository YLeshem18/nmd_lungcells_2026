#!/usr/bin/env Rscript
# SFx - Normalization Sanity Check (productive analysis): non-NMD genes, DMSO productive CPM
# vs custom-adjusted SMG1i productive CPM, per cell type (controls land on the identity line;
# Spearman >= 0.96). Also reports rank preservation (custom vs uncorrected/Kitagawa logFC, Spearman >= 0.99).
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(tidyr); library(ggplot2); library(scales) })
select<-dplyr::select; filter<-dplyr::filter; rename<-dplyr::rename
options(bitmapType="cairo")
DATA<-"nmd_fig_data"
HERE<-"supplement_figures"
dir.create(HERE, recursive=TRUE, showWarnings=FALSE)
P<-function(f) file.path(DATA,f); CT<-c("AT2","LAE","FB","MV"); LF<-0.05; LF_UNS<-0.20
FLOOR_BOTH<-1; FLOOR_DMSO<-2; pass_floor<-function(td,ts)(td>=FLOOR_BOTH&ts>=FLOOR_BOTH)|(td>=FLOOR_DMSO)

dge<-readRDS(P("dge_isoform_longread_filtered_2026.3.3.rds"))
pm<-read.csv(P("mashr_isoform_posterior_means_2026.3.10.csv"),check.names=FALSE)
lf<-read.csv(P("mashr_isoform_lfsr_2026.3.10.csv"),check.names=FALSE)
colnames(pm)[1]<-"isoform_id"; colnames(lf)[1]<-"isoform_id"
i2g<-data.frame(isoform_id=dge$genes$txid, gene_id=dge$genes$gene_id)
samp<-as.data.frame(dge$samples); samp$ct<-as.character(samp$ct); samp$treatment<-as.character(samp$treatment)
samp$ct[samp$ct=="AT"]<-"AT2"; samp$ct[samp$ct=="DD"]<-"LAE"; samp$sid<-as.character(samp$id)
samp<-samp[samp$ct%in%CT,]; counts<-dge$counts[,rownames(samp)]
resolve<-function(df,ct) df[[ grep(paste0("_in_",ct,"$"),colnames(df))[1] ]]

build_adj<-function(ct){
  l<-resolve(lf,ct); p<-resolve(pm,ct)
  unp<-!is.na(l)&!is.na(p)&l<LF&p>0; uns<-!is.na(l)&!is.na(p)&l>=LF&l<LF_UNS&p>0; dwn<-!is.na(l)&!is.na(p)&l<LF&p<0
  cls<-setNames(ifelse(unp,"unprod","prod"),lf$isoform_id)
  cs<-rownames(samp)[samp$ct==ct]; cnt<-counts[,cs,drop=FALSE]
  trt<-setNames(samp[cs,"treatment"],cs); don<-setNames(samp[cs,"sid"],cs)
  raw<-sweep(cnt,2,colSums(cnt),"/")*1e6; ids<-rownames(raw); g<-i2g$gene_id[match(ids,i2g$isoform_id)]; cl<-cls[ids]
  unpg<-unique(g[!is.na(g)&cl=="unprod"]); unsg<-setdiff(unique(g[!is.na(g)&(uns[ids]|dwn[ids])]),unpg)
  anc<-!is.na(g)&!grepl("::",g,fixed=TRUE)&!is.na(cl)
  dmso<-cs[trt[cs]=="DMSO"]; smg<-cs[trt[cs]=="Smg1i"]; dm<-rowMeans(raw[,dmso,drop=FALSE])
  pidx<-anc&cl=="prod"; Pd<-tapply(dm[pidx],g[pidx],sum)
  iu<-which(anc&cl=="unprod"); pg<-Pd[g[iu]]; hp<-!is.na(pg)&pg>0; sh<-ifelse(hp,dm[iu]/pg,NA); ia<-iu[hp]; ina<-iu[!hp]
  adj<-raw
  for(s in smg){col<-raw[,s]; Pgs<-tapply(col[pidx],g[pidx],sum); nc<-col
    nc[ia]<-Pgs[g[ia]]*sh[hp]; nc[ina]<-dm[ina]; nc[!is.finite(nc)]<-0; adj[,s]<-nc*(1e6/sum(nc))}
  list(raw=raw,adj=adj,gene=g,cl=cl,anc=anc,trt=trt,don=don,unpg=unpg,unsg=unsg)
}
process_ct<-function(ct){
  o<-build_adj(ct); k<-o$anc; fk<-paste(o$gene[k],o$cl[k],sep="|||")
  rf<-rowsum(o$raw[k,,drop=FALSE],fk); af<-rowsum(o$adj[k,,drop=FALSE],fk)
  fm<-data.frame(feature=rownames(rf)); fm$gene<-sub("\\|\\|\\|.*","",fm$feature); fm$class<-sub(".*\\|\\|\\|","",fm$feature)
  sm<-data.frame(sample=colnames(rf),trt=o$trt[colnames(rf)])
  melt<-function(m,v){d<-as.data.frame(m);d$feature<-rownames(m);pivot_longer(d,-feature,names_to="sample",values_to=v)}
  long<-melt(af,"adj")%>%left_join(melt(rf,"raw"),by=c("feature","sample"))%>%left_join(fm,by="feature")%>%left_join(sm,by="sample")
  agg<-long%>%group_by(gene,class,trt)%>%summarise(adj=mean(adj),raw=mean(raw),.groups="drop")
  dmw<-agg%>%filter(trt=="DMSO")%>%select(gene,class,adj)%>%pivot_wider(names_from=class,values_from=adj,values_fill=0)%>%rename(prodCPM_DMSO=prod,unprodCPM_DMSO=unprod)
  sg<-agg%>%filter(trt=="Smg1i")%>%select(gene,class,adj,raw)%>%pivot_wider(names_from=class,values_from=c(adj,raw),values_fill=0)%>%
    rename(prodCPM_SMG_adj=adj_prod,prodCPM_SMG_raw=raw_prod,unprodCPM_SMG=raw_unprod)
  full_join(dmw,sg,by="gene")%>%mutate(ct=ct,
    gene_category=case_when(gene%in%o$unpg~"NMD",gene%in%o$unsg~"unsure",TRUE~"non-NMD"))
}
gs<-bind_rows(lapply(CT,process_ct))
base<-gs%>%filter(!is.na(prodCPM_DMSO),!is.na(prodCPM_SMG_raw))%>%
  mutate(across(c(prodCPM_DMSO,unprodCPM_DMSO,prodCPM_SMG_raw,unprodCPM_SMG,prodCPM_SMG_adj),~coalesce(.x,0)),
         totalCPM_DMSO=prodCPM_DMSO+unprodCPM_DMSO,totalCPM_SMG=prodCPM_SMG_raw+unprodCPM_SMG)%>%
  filter(pass_floor(totalCPM_DMSO,totalCPM_SMG))%>%
  mutate(custom=log2((prodCPM_SMG_adj+1)/(prodCPM_DMSO+1)),
         Kitagawa=log2((prodCPM_SMG_raw+1)/(prodCPM_DMSO+1)))
base$ct<-factor(base$ct,levels=CT)

nn<-base%>%filter(gene_category=="non-NMD")
labs<-nn%>%group_by(ct)%>%summarise(
  sp=round(cor(prodCPM_DMSO,prodCPM_SMG_adj,method="spearman"),2),
  slope=round(coef(lm(log2(prodCPM_SMG_adj+1)~log2(prodCPM_DMSO+1)))[2],2),.groups="drop")
rank_pres<-base%>%filter(gene_category%in%c("NMD","non-NMD"))%>%group_by(ct)%>%
  summarise(spearman_custom_vs_uncorrected=round(cor(custom,Kitagawa,method="spearman"),3),.groups="drop")
cat("Non-NMD DMSO vs adjusted-SMG1i (Spearman, slope):\n"); print(as.data.frame(labs))
cat("\nRank preservation custom vs uncorrected logFC (Spearman):\n"); print(as.data.frame(rank_pres))
fwrite(labs%>%left_join(rank_pres,by="ct"), file.path(HERE,"data_norm_sanity.csv"))

xmax<-max(log2(nn$prodCPM_DMSO+1),log2(nn$prodCPM_SMG_adj+1))
p<-ggplot(nn, aes(log2(prodCPM_DMSO+1), log2(prodCPM_SMG_adj+1)))+
  geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey55")+
  geom_point(alpha=0.10,size=0.5,colour="#4C78B0")+
  geom_smooth(method="lm",se=FALSE,colour="#C0392B",linewidth=0.7)+
  geom_text(data=labs,aes(x=0,y=xmax,label=sprintf("Spearman = %.2f\nslope = %.2f",sp,slope)),
            hjust=0,vjust=1,size=4,inherit.aes=FALSE)+
  facet_wrap(~ct,ncol=2)+
  labs(x="log2(DMSO productive CPM + 1)", y="log2(adjusted SMG1i productive CPM + 1)")+
  theme_minimal(base_size=13)+
  theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold",size=15),
        plot.margin=margin(8,10,8,8))
ggsave(file.path(HERE,"figureS3_norm_sanity.png"), p, width=10, height=8, dpi=300, bg="white")
ggsave(file.path(HERE,"figureS3_norm_sanity.pdf"), p, width=10, height=8, bg="white")
cat("\nWrote figureS3_norm_sanity.{png,pdf}\n")
