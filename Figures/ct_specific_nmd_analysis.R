#!/usr/bin/env Rscript
# Cell-type-specific NMD: isoform-level drivers + gene-vs-isoform sharing + short-read gene level.
# Writes a self-contained HTML report with embedded figures.
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tibble); library(edgeR)
  library(ggplot2); library(tidyr); library(jsonlite); library(patchwork)
})
options(bitmapType="cairo")
DATA<-"nmd_fig_data"; P<-function(f) file.path(DATA,f)
OUT <-"ct_specific_nmd_report"
dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
CT<-c("AT2","LAE","FB","MV"); CODE<-c(AT2="at2",LAE="dd",FB="fb",MV="mv"); EXPR<-1
CTCOL<-c(AT2="#e74c3c",LAE="#3498db",FB="#16a085",MV="#9b59b6")
imgs<-list(); tbls<-list()
emb<-function(path) sprintf('<img src="data:image/png;base64,%s" style="max-width:100%%;border:1px solid #ddd;border-radius:6px;">',
                            jsonlite::base64_enc(readBin(path,"raw",file.info(path)$size)))
saveimg<-function(p,name,w=8,h=4){f<-file.path(OUT,name);ggsave(f,p,width=w,height=h,dpi=140,bg="white");emb(f)}
tbl<-function(df){ df<-as.data.frame(df)
  paste0("<table><thead><tr>",paste0("<th>",names(df),"</th>",collapse=""),"</tr></thead><tbody>",
    paste0(apply(df,1,function(r) paste0("<tr>",paste0("<td>",r,"</td>",collapse=""),"</tr>")),collapse=""),
    "</tbody></table>") }

# ============================ ISOFORM LEVEL ============================
dge<-readRDS(P("dge_isoform_longread_2026.3.3.rds"))
s<-as.data.frame(dge$samples)%>%rownames_to_column("bam")%>%mutate(ct=as.character(ct),treatment=as.character(treatment))
s$ct[s$ct=="AT"]<-"AT2"; s$ct[s$ct=="DD"]<-"LAE"; s<-s[s$ct%in%CT,]
cpm<-cpm(dge[,s$bam],normalized.lib.sizes=FALSE,log=FALSE); tx<-rownames(cpm)
smg<-sapply(CT,function(ct) rowMeans(cpm[,s$bam[s$ct==ct&s$treatment=="Smg1i"],drop=FALSE]))
dm <-sapply(CT,function(ct) rowMeans(cpm[,s$bam[s$ct==ct&s$treatment=="DMSO"],drop=FALSE]))
gene<-dge$genes$gene_id; names(gene)<-tx
nmdmat<-sapply(CT,function(ct){f<-list.files(DATA,pattern=sprintf("nmd_mashr_die_%s_.*csv",CODE[[ct]]),full.names=TRUE)[1]
  d<-fread(f);setNames(d$nmd_responsive%in%c(TRUE,"TRUE"),d$txid)[tx]}); nmdmat[is.na(nmdmat)]<-FALSE; rownames(nmdmat)<-tx
n_nmd<-rowSums(nmdmat); spec<-which(n_nmd==1); shared<-which(n_nmd>=2)

iso_share<-data.frame(n_cell_types=1:4, n_isoforms=as.integer(table(factor(n_nmd[n_nmd>=1],1:4))))
tbls$iso_share<-iso_share

# specific: producibility (SMG1i) + host gene expressed (DMSO)
S<-smg[spec,,drop=FALSE]; N<-nmdmat[spec,,drop=FALSE]; col<-max.col(N,ties.method="first")
Soth<-S; Soth[cbind(1:nrow(S),col)]<-NA
prod_elsewhere<-apply(Soth>=EXPR,1,any,na.rm=TRUE)
ok<-!is.na(gene)&gene!=""; gdm<-rowsum(dm[ok,,drop=FALSE],gene[ok]); geneidx<-match(gene,rownames(gdm))
Gd<-gdm[geneidx[spec],,drop=FALSE]; Goth<-Gd; Goth[cbind(1:nrow(Gd),col)]<-NA
gene_elsewhere<-apply(Goth>=EXPR,1,any,na.rm=TRUE)
spec_drivers<-data.frame(
  metric=c("Isoform produced (SMG1i>=1) in >=1 other CT","Isoform NOT produced elsewhere",
           "Host gene expressed (DMSO>=1) in >=1 other CT"),
  pct=round(100*c(mean(prod_elsewhere),mean(!prod_elsewhere),mean(gene_elsewhere,na.rm=TRUE)),1))
tbls$spec_drivers<-spec_drivers

# DMSO NMD-CT vs non-NMD-CT (producible-elsewhere set)
D<-dm[spec,,drop=FALSE]; dmso_nmd<-D[cbind(1:nrow(D),col)]
Dp<-D; Dp[cbind(1:nrow(D),col)]<-NA; Dp[!(Soth>=EXPR)]<-NA
dmso_oth<-rowMeans(Dp,na.rm=TRUE)
pe<-prod_elsewhere & is.finite(dmso_oth)
dmso_dt<-data.table(dmso_nmd=dmso_nmd[pe],dmso_oth=dmso_oth[pe])
w<-wilcox.test(dmso_dt$dmso_nmd,dmso_dt$dmso_oth,paired=TRUE)
dmso_tab<-data.frame(
  quantity=c("median DMSO in NMD cell type","median DMSO in non-NMD cell type(s)",
             "% isoforms lower in NMD CT","paired Wilcoxon p"),
  value=c(round(median(dmso_dt$dmso_nmd),3),round(median(dmso_dt$dmso_oth),3),
          paste0(round(100*mean(dmso_dt$dmso_nmd<dmso_dt$dmso_oth),1),"%"),
          signif(w$p.value,3)))
tbls$dmso_tab<-dmso_tab
# figure C: paired DMSO
fc_dat<-melt(dmso_dt,measure.vars=c("dmso_nmd","dmso_oth"),variable.name="grp",value.name="cpm")
fc_dat$grp<-factor(ifelse(fc_dat$grp=="dmso_nmd","NMD cell type","non-NMD cell type"),
                   levels=c("NMD cell type","non-NMD cell type"))
figC<-ggplot(fc_dat,aes(grp,log2(cpm+1),fill=grp))+
  geom_boxplot(width=0.55,outlier.size=0.3,outlier.alpha=0.2,colour="grey25")+
  scale_fill_manual(values=c("NMD cell type"="#c0392b","non-NMD cell type"="#4C78B0"),guide="none")+
  labs(x=NULL,y="DMSO expression, log2(CPM+1)",
       title="Isoform DMSO expression: NMD vs non-NMD cell type (specific isoforms produced in both)")+
  theme_bw(base_size=12)
imgs$C<-saveimg(figC,"figC_dmso_nmd_vs_not.png",7,4)

# same gene, different isoform, other CT
long<-rbindlist(lapply(CT,function(c){t<-tx[nmdmat[,c]];data.table(tx=t,gene=gene[t],ct=c)}))[!is.na(gene)&gene!=""]
total_pairs<-long[,.(total=.N),by=gene]; in_ct<-long[,.(n_in_ct=.N),by=.(gene,ct)]
sc<-CT[max.col(nmdmat[tx[spec],,drop=FALSE],ties.method="first")]
sdt<-data.table(tx=tx[spec],gene=gene[spec],ct=sc)[!is.na(gene)&gene!=""]
sdt<-merge(merge(sdt,total_pairs,by="gene",all.x=TRUE),in_ct,by=c("gene","ct"),all.x=TRUE)
sdt[,other_iso_other_ct:=total-n_in_ct]
split_tab<-data.frame(
  category=c("Same gene NMD via a DIFFERENT isoform in a DIFFERENT cell type",
             "Gene NMD only in this one cell type"),
  n=c(sum(sdt$other_iso_other_ct>=1),sum(sdt$other_iso_other_ct==0)))
split_tab$pct<-round(100*split_tab$n/sum(split_tab$n),1)
tbls$split_tab<-split_tab
figA<-ggplot(split_tab,aes(x=reorder(category,-n),y=n,fill=category))+
  geom_col(width=0.6)+geom_text(aes(label=sprintf("%d (%.1f%%)",n,pct)),vjust=-0.4,size=4)+
  scale_fill_manual(values=c("#2166AC","#F4A582"),guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0,0.12)))+
  labs(x=NULL,y="Cell-type-specific NMD isoforms",
       title="Is the specific isoform's gene NMD-regulated elsewhere via a different isoform?")+
  theme_bw(base_size=12)+theme(axis.text.x=element_text(size=9))
imgs$A<-saveimg(figA,"figA_split.png",8,4.2)

# gene- vs isoform-level sharing
gene_nmd_ct<-unique(long[,.(gene,ct)])            # gene NMD in CT if >=1 NMD isoform
gshare<-gene_nmd_ct[,.(nct=.N),by=gene]
gene_share<-data.frame(n_cell_types=1:4, n_genes=as.integer(table(factor(gshare$nct,1:4))))
tbls$gene_share<-gene_share
shr<-rbind(
  data.frame(level="Isoform",n_cell_types=1:4,pct=100*iso_share$n_isoforms/sum(iso_share$n_isoforms)),
  data.frame(level="Gene",   n_cell_types=1:4,pct=100*gene_share$n_genes/sum(gene_share$n_genes)))
figB<-ggplot(shr,aes(factor(n_cell_types),pct,fill=level))+
  geom_col(position=position_dodge(0.75),width=0.7)+
  geom_text(aes(label=sprintf("%.0f%%",pct)),position=position_dodge(0.75),vjust=-0.3,size=3.5)+
  scale_fill_manual(values=c("Isoform"="#5B8FF9","Gene"="#F6BD16"),name=NULL)+
  scale_y_continuous(expand=expansion(mult=c(0,0.12)))+
  labs(x="Number of cell types NMD-susceptible",y="% of features",
       title="NMD sharing is higher at the gene level than the isoform level")+
  theme_bw(base_size=12)+theme(legend.position="top")
imgs$B<-saveimg(figB,"figB_gene_vs_iso_sharing.png",8,4.2)
# reconciliation stat
genes_with_spec<-unique(sdt$gene)
gshare_map<-setNames(gshare$nct,gshare$gene)
recon<-data.frame(
  statistic=c("Genes with >=1 cell-type-specific NMD isoform",
              "  of those, gene is NMD in >=2 cell types (shared gene)",
              "  of those, gene is NMD in only 1 cell type"),
  value=c(length(genes_with_spec),
          sprintf("%d (%.1f%%)",sum(gshare_map[genes_with_spec]>=2),100*mean(gshare_map[genes_with_spec]>=2)),
          sprintf("%d (%.1f%%)",sum(gshare_map[genes_with_spec]==1),100*mean(gshare_map[genes_with_spec]==1))))
tbls$recon<-recon

# ============================ SHORT-READ GENE LEVEL ============================
g<-readRDS(P("dge_gene_unfiltered_2026.1.2.rds"))
gs<-as.data.frame(g$samples); gs$ct<-as.character(gs$ct); gs$treatment<-as.character(gs$treatment)
gs$ct[gs$ct=="AT"]<-"AT2"; gs$ct[gs$ct=="DD"]<-"LAE"
keep<-gs$ct%in%CT; gs<-gs[keep,]; gc_mat<-cpm(g$counts,normalized.lib.sizes=FALSE,log=FALSE)[,rownames(g$samples)[keep]]
rownames(gc_mat)<-sub("\\..*","",rownames(gc_mat))
sr_dmso<-sapply(CT,function(ct) rowMeans(gc_mat[,rownames(gs)[gs$ct==ct&gs$treatment=="DMSO"],drop=FALSE]))
srlist<-lapply(CT,function(ct){f<-list.files(DATA,pattern=sprintf("nmd_mashr_dge_%s_.*csv",CODE[[ct]]),full.names=TRUE)[1]
  d<-fread(f);gid<-sub("\\..*","",d$ensembl_gene_id_version);setNames(d$nmd_responsive%in%c(TRUE,"TRUE"),gid)}); names(srlist)<-CT
allg<-Reduce(union,lapply(srlist,names))
srmat<-sapply(CT,function(ct){v<-srlist[[ct]][allg];v[is.na(v)]<-FALSE;v}); rownames(srmat)<-allg
srn<-rowSums(srmat); sr_spec<-which(srn==1); sr_shared<-which(srn>=2)
sr_share<-data.frame(n_cell_types=1:4,n_genes=as.integer(table(factor(srn[srn>=1],1:4))))
tbls$sr_share<-sr_share
# expression of SR-specific genes in other CTs
common<-intersect(rownames(sr_dmso),names(sr_spec))
SM<-srmat[common,,drop=FALSE]; E<-sr_dmso[common,,drop=FALSE]
scol<-max.col(SM,ties.method="first")
e_nmd<-E[cbind(1:nrow(E),scol)]; Eo<-E; Eo[cbind(1:nrow(E),scol)]<-NA
e_oth_mean<-rowMeans(Eo,na.rm=TRUE); e_oth_max<-apply(Eo,1,max,na.rm=TRUE)
sr_expr_tab<-data.frame(
  quantity=c("SR-specific NMD genes analyzed",
             "Gene expressed (DMSO>=1) in >=1 other CT",
             "median DMSO CPM in NMD cell type",
             "median DMSO CPM in other cell types (mean)",
             "% with DMSO lower in NMD cell type",
             "paired Wilcoxon p (NMD-CT vs other-CT mean)"),
  value=c(nrow(E),
          sprintf("%.1f%%",100*mean(e_oth_max>=EXPR)),
          round(median(e_nmd),2),
          round(median(e_oth_mean),2),
          sprintf("%.1f%%",100*mean(e_nmd<e_oth_mean)),
          signif(wilcox.test(e_nmd,e_oth_mean,paired=TRUE)$p.value,3)))
tbls$sr_expr_tab<-sr_expr_tab
srfig_dat<-melt(data.table(NMD=e_nmd,other=e_oth_mean),measure.vars=c("NMD","other"),
                variable.name="grp",value.name="cpm")
srfig_dat$grp<-factor(ifelse(srfig_dat$grp=="NMD","NMD cell type","other cell types"),
                      levels=c("NMD cell type","other cell types"))
figD<-ggplot(srfig_dat,aes(grp,log2(cpm+1),fill=grp))+
  geom_boxplot(width=0.55,outlier.size=0.3,outlier.alpha=0.2,colour="grey25")+
  scale_fill_manual(values=c("NMD cell type"="#c0392b","other cell types"="#4C78B0"),guide="none")+
  labs(x=NULL,y="Short-read gene DMSO expression, log2(CPM+1)",
       title="SR gene DMSO expression: NMD cell type vs other cell types (genes NMD in only 1 CT)")+
  theme_bw(base_size=12)
imgs$D<-saveimg(figD,"figD_sr_gene_expr.png",7,4)

# composite supplement figure (A-D)
comp<-(figA + labs(title="A. Specific isoform's gene NMD elsewhere via a different isoform")) /
      (figB + labs(title="B. NMD sharing: gene vs isoform level")) /
      ((figC + labs(title="C. Isoform DMSO: NMD vs non-NMD CT")) |
       (figD + labs(title="D. Short-read gene DMSO: single-CT-NMD genes")))
ggsave(file.path(OUT,"figS_ct_specific_nmd_composite.png"),comp,width=11,height=13,dpi=150,bg="white")
ggsave(file.path(OUT,"figS_ct_specific_nmd_composite.pdf"),comp,width=11,height=13,bg="white")

# ============================ HTML ============================
css<-"body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;max-width:1000px;margin:24px auto;padding:0 18px;color:#222;line-height:1.5}
h1{border-bottom:3px solid #2166AC;padding-bottom:6px} h2{color:#2166AC;margin-top:34px;border-bottom:1px solid #ddd}
table{border-collapse:collapse;margin:12px 0;font-size:14px} th,td{border:1px solid #ccc;padding:5px 10px;text-align:left}
th{background:#f0f4f8} .k{background:#eef7ee;padding:10px 14px;border-left:4px solid #16a085;border-radius:4px;margin:12px 0}
figure{margin:16px 0}"
H<-function(...) paste0(...)
html<-H("<!doctype html><html><head><meta charset='utf-8'><title>Cell-type-specific NMD analysis</title><style>",css,"</style></head><body>",
"<h1>What drives cell-type-specific NMD susceptibility?</h1>",
"<p>Analysis of long-read isoform-level NMD calls (mashr <code>nmd_responsive</code>, lfsr&lt;0.05 &amp; posterior mean&gt;0) across AT2, LAE, FB, MV, plus a short-read gene-level parallel. 'Expressed/produced in a cell type' = mean CPM &ge; 1 (DMSO for baseline, SMG1i for de-repressed production).</p>",
"<div class='k'><b>Summary.</b> Cell-type-specific NMD is driven by <b>which isoform each cell type makes, not by whether the gene is an NMD target or how much it is expressed.</b> (i) 80.5% of cell-type-specific NMD isoforms belong to genes that are NMD-regulated in other cell types too, via a different isoform &mdash; so NMD-target <i>gene</i> identity is largely conserved while the triggering isoform varies. (ii) At the short-read gene level, genes called NMD in a single cell type are expressed at essentially the <b>same</b> level in the other cell types (median 19.3 vs 18.6 CPM; 75.8% expressed elsewhere) &mdash; expression gating does <b>not</b> explain specificity. (iii) For the minority (~10%) of specific isoforms produced in multiple cell types, the isoform sits ~2.6&times; lower at DMSO where it is NMD-targeted (p&asymp;10<sup>-63</sup>), a footprint of genuine cell-type-selective decay. Together: genes are broadly and equally available across cell types; cell-type-specific alternative splicing determines which isoform is made and degraded, and gene-level totals barely move because the NMD isoform is a minor fraction.</div>",

"<h2>1. Isoform-level sharing</h2>",tbl(tbls$iso_share),

"<h2>2. Why are specific isoforms specific?</h2>",tbl(tbls$spec_drivers),
"<div class='k'>~90% of cell-type-specific NMD isoforms are simply <b>not produced</b> in the other cell types, yet the <b>host gene is expressed</b> in other cell types ~96% of the time &mdash; so specificity is upstream, in which isoform each cell type splices/makes, not in whether the gene is on.</div>",

"<h2>3. Differential decay: DMSO expression, NMD vs non-NMD cell type</h2>",
"<p>For the ~10% of specific isoforms that <i>are</i> produced in other cell types, is the isoform lower at baseline where it is NMD-targeted?</p>",
tbl(tbls$dmso_tab),imgs$C,
"<div class='k'>Where the isoform is NMD-susceptible it sits ~2.6&times; lower at DMSO than in cell types where it escapes &mdash; a footprint of active, cell-type-selective degradation visible even without SMG1i.</div>",

"<h2>4. Same gene, different isoform (the key result)</h2>",tbl(tbls$split_tab),imgs$A,
"<div class='k'>80.5% of 'cell-type-specific' NMD isoforms belong to genes that are NMD-regulated in other cell types too &mdash; via a <b>different isoform</b>. Specificity is at the isoform level; the gene remains a conserved NMD target.</div>",

"<h2>5. Gene- vs isoform-level sharing reconciliation</h2>",tbl(tbls$recon),imgs$B,
"<p>Gene-level distribution:</p>",tbl(tbls$gene_share),

"<h2>6. Short-read gene level: expression of single-CT NMD genes elsewhere</h2>",
"<p>Genes called NMD in only one cell type (short-read mashr) &mdash; are they expressed in the other cell types, and is their DMSO expression lower where they are NMD-targeted?</p>",
tbl(tbls$sr_share),tbl(tbls$sr_expr_tab),imgs$D,
"</body></html>")
writeLines(html, file.path(OUT,"ct_specific_nmd_report.html"))
cat("Wrote", file.path(OUT,"ct_specific_nmd_report.html"),"\n")
cat("\n== key numbers ==\n")
cat("Isoform specific/shared:",sum(n_nmd==1),"/",sum(n_nmd>=2),"\n")
print(split_tab); print(recon); print(sr_share); print(sr_expr_tab)
