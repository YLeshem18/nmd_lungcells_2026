#!/usr/bin/env Rscript
# SFx - Reproducibility battery, restricted to the 4 tests retained in the manuscript:
#   1) Sign-flip (Rademacher) permutation test  -> permutation p
#   2) Cross-donor sign concordance vs exact binomial -> fold over expected (+ p)
#   3) Effect-size variance ratio vs sign-flip null
#   4) Storey pi0 -> fraction of genes carrying a real effect (1 - pi0)
# Computed on donor-paired adjusted (custom) productive log2FC, MAIN_CATS (NMD + non-NMD).
# Mirrors productive_response.Rmd noise tests. FB shown but flagged provisional.
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(limma)
})
select<-dplyr::select; filter<-dplyr::filter; rename<-dplyr::rename
options(bitmapType="cairo")
DATA<-"nmd_fig_data"
HERE<-"supplement_figures"
dir.create(HERE, recursive=TRUE, showWarnings=FALSE)
P<-function(f) file.path(DATA,f)
CT<-c("AT2","LAE","FB","MV"); LF<-0.05; LF_UNS<-0.20
rowSds<-function(X){n<-ncol(X); sqrt((rowMeans(X^2)-rowMeans(X)^2)*n/(n-1))}

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
  list(adj=adj,gene=g,cl=cl,anc=anc,trt=trt,don=don,unpg=unpg,unsg=unsg)
}

# per-donor adjusted productive log2FC matrix (gene x donor), MAIN_CATS genes
prod_logfc_mat<-function(ct){
  o<-build_adj(ct); k<-o$anc&o$cl=="prod"
  cs<-colnames(o$adj); don<-o$don[cs]; trt<-o$trt[cs]
  pe<-rowsum(o$adj[k,,drop=FALSE], o$gene[k])
  cat_of<-ifelse(rownames(pe)%in%o$unpg,"NMD",ifelse(rownames(pe)%in%o$unsg,"unsure","non-NMD"))
  pe<-pe[cat_of%in%c("NMD","non-NMD"),,drop=FALSE]
  dons<-unique(don)
  M<-sapply(dons,function(d){ dc<-cs[don==d&trt=="DMSO"]; sc<-cs[don==d&trt=="Smg1i"]
    if(length(dc)!=1||length(sc)!=1) return(rep(NA,nrow(pe)))
    log2((pe[,sc]+1)/(pe[,dc]+1)) })
  rownames(M)<-rownames(pe); M[rowSums(is.finite(M))==ncol(M),,drop=FALSE]
}

# limma-trend P-values (for Storey), same fit as fit_custom_limma
limmaP<-function(ct){
  o<-build_adj(ct); k<-o$anc&o$cl=="prod"; pe<-rowsum(o$adj[k,,drop=FALSE],o$gene[k])
  cat_of<-ifelse(rownames(pe)%in%o$unpg,"NMD",ifelse(rownames(pe)%in%o$unsg,"unsure","non-NMD"))
  cs<-colnames(pe); sp<-data.frame(row.names=cs,treatment=relevel(factor(o$trt[cs]),ref="DMSO"),sid=factor(o$don[cs]))
  keep<-rowSums(pe>=1)>=ceiling(ncol(pe)/2); E<-log2(pe[keep,,drop=FALSE]+1)
  des<-model.matrix(~treatment,data=sp); cf<-duplicateCorrelation(E,des,block=sp$sid)
  fit<-eBayes(lmFit(E,des,block=sp$sid,correlation=cf$consensus.correlation),trend=TRUE)
  ci<-grep("Smg1i",colnames(des)); tt<-topTable(fit,coef=ci,number=Inf,sort.by="none")
  data.frame(gene_id=rownames(tt), P=tt$P.Value, cat=cat_of[match(rownames(tt),rownames(pe))])
}
storey_pi0<-function(p,lambda=0.5){p<-p[is.finite(p)]; min(1, mean(p>lambda)/(1-lambda))}

res<-list()
for(ct in CT){
  M<-prod_logfc_mat(ct); n<-ncol(M); G<-nrow(M)
  # 1) sign-flip permutation
  tstat<-function(X) rowMeans(X)/(rowSds(X)/sqrt(n)); cut<-qt(0.975,n-1)
  obs<-sum(abs(tstat(M))>=cut,na.rm=TRUE); set.seed(1); B<-2000; null<-numeric(B)
  for(b in 1:B){S<-matrix(sample(c(-1,1),G*n,TRUE),G,n); null[b]<-sum(abs(tstat(M*S))>=cut,na.rm=TRUE)}
  perm_p<-(sum(null>=obs)+1)/(B+1)
  # 2) concordance
  allpos<-sum(rowSums(M>0)==n); allneg<-sum(rowSums(M<0)==n); conc<-allpos+allneg
  exp_conc<-G*2^(1-n); fold<-conc/exp_conc
  bt<-binom.test(conc,G,p=2^(1-n),alternative="greater")
  # 3) variance ratio
  set.seed(1); effnull<-unlist(lapply(1:200,function(b){S<-matrix(sample(c(-1,1),G*n,TRUE),G,n); rowMeans(M*S)}))
  vr<-var(rowMeans(M))/var(effnull)
  # 4) Storey pi0
  lp<-limmaP(ct); p<-lp$P[lp$cat%in%c("NMD","non-NMD")]; pi0<-storey_pi0(p); nonnull<-1-pi0
  res[[ct]]<-data.frame(ct=ct, n_genes=G, perm_p=perm_p, concordance_fold=round(fold,2),
                        variance_ratio=round(vr,2), nonnull_fraction=round(nonnull,3),
                        storey_pi0=round(pi0,3))
}
batt<-bind_rows(res); cat("Reproducibility battery (real values):\n"); print(as.data.frame(batt))
fwrite(batt, file.path(HERE,"data_reproducibility_battery.csv"))

# ---- figure: tiles, tests x cell type. Pass threshold baked into each row label. ----
long<-batt %>% transmute(ct,
  `Sign-flip permutation test\n(pass: permutation p < 0.05)`=perm_p,
  `Cross-donor sign concordance\n(pass: fold > 1× expected)`=concordance_fold,
  `Effect-size variance ratio\n(pass: obs / sign-flip null > 1)`=variance_ratio,
  `Storey π₀ non-null fraction\n(pass: 1−π₀ > 0.1)`=nonnull_fraction) %>%
  pivot_longer(-ct, names_to="test", values_to="stat")
rule<-function(test,stat) dplyr::case_when(
  grepl("permutation",test) ~ stat<0.05,
  grepl("concordance",test) ~ stat>1,
  grepl("variance",test)    ~ stat>1,
  grepl("Storey",test)      ~ stat>0.1)
long<-long %>% mutate(pass=rule(test,stat),
  label=ifelse(grepl("permutation",test),
               ifelse(stat<1e-3, formatC(stat,format="e",digits=1), sprintf("%.3f",stat)),
               sprintf("%.2f",stat)),
  ct=factor(ct,levels=CT),
  test=factor(test,levels=rev(unique(test))))
p<-ggplot(long,aes(ct,test,fill=pass))+
  geom_tile(color="white",linewidth=1.1)+
  geom_text(aes(label=label),fontface="bold",size=4.4,color=ifelse(long$pass,"grey10","white"))+
  scale_fill_manual(values=c(`TRUE`="#4393c3",`FALSE`="#d6604d"),
                    labels=c(`TRUE`="meets criterion",`FALSE`="fails criterion"),name=NULL)+
  scale_x_discrete(position="top")+
  labs(x=NULL,y=NULL)+
  theme_minimal(base_size=12)+
  theme(panel.grid=element_blank(), axis.text.x.top=element_text(face="bold",size=15,colour="grey10"),
        axis.text.y=element_text(size=13,colour="grey10",face="bold",lineheight=0.95),
        legend.position="top", legend.text=element_text(size=12),
        plot.margin=margin(8,10,8,8))
ggsave(file.path(HERE,"figureS3_reproducibility.png"),p,width=9,height=5,dpi=300,bg="white")
ggsave(file.path(HERE,"figureS3_reproducibility.pdf"),p,width=9,height=5,bg="white")
cat("Wrote figureS3_reproducibility.{png,pdf}\n")
