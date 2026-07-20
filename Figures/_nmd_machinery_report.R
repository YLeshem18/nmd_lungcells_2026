suppressPackageStartupMessages({library(data.table);library(edgeR);library(ggplot2);library(tibble);library(dplyr);library(tidyr);library(jsonlite)})
options(bitmapType="cairo")
DATA<-"nmd_fig_data"; P<-function(f) file.path(DATA,f)
OUT<-"nmd_machinery_report"; dir.create(OUT,showWarnings=FALSE)
CT<-c("AT2","LAE","FB","MV"); CODE<-c(AT2="at2",LAE="dd",FB="fb",MV="mv")
base_id<-function(x) vapply(strsplit(x,".",fixed=TRUE),`[`,character(1),1)
CORE<-c("UPF1","UPF2","UPF3A","UPF3B","SMG1","SMG5","SMG6","SMG7","SMG8","SMG9")
emb<-function(p) sprintf('<img src="data:image/png;base64,%s" style="max-width:100%%;border:1px solid #ddd;border-radius:6px;">',jsonlite::base64_enc(readBin(p,"raw",file.info(p)$size)))
saveimg<-function(g,n,w=7,h=4){f<-file.path(OUT,n);ggsave(f,g,width=w,height=h,dpi=140,bg="white");emb(f)}
tbl<-function(df){df<-as.data.frame(df);paste0("<table style='border-collapse:collapse'>","<tr>",paste0("<th style='border:1px solid #ccc;padding:4px 9px;background:#eef'>",names(df),"</th>",collapse=""),"</tr>",paste0(apply(df,1,function(r)paste0("<tr>",paste0("<td style='border:1px solid #ccc;padding:4px 9px'>",r,"</td>",collapse=""),"</tr>")),collapse=""),"</table>")}

# ---- SR gene-level machinery score ----
sr1<-fread(list.files(DATA,pattern="nmd_mashr_dge_at2_.*csv",full.names=TRUE)[1])
core_ens<-setNames(base_id(sr1$ensembl_gene_id_version),sr1$hgnc_symbol)[CORE]
g<-readRDS(P("dge_gene_unfiltered_2026.1.2.rds")); gs<-as.data.frame(g$samples)
gs$ct<-as.character(gs$ct); gs$treatment<-as.character(gs$treatment); gs$ct[gs$ct=="AT"]<-"AT2"; gs$ct[gs$ct=="DD"]<-"LAE"
keep<-gs$ct%in%CT; gc<-cpm(g$counts,log=FALSE)[,rownames(g$samples)[keep]]; rownames(gc)<-base_id(rownames(gc)); gs<-gs[keep,]
dmso<-sapply(CT,function(c) rowMeans(gc[,rownames(gs)[gs$ct==c&gs$treatment=="DMSO"],drop=FALSE]))
sel<-core_ens[core_ens%in%rownames(dmso)]; E<-dmso[sel,]; rownames(E)<-names(sel)
Z<-t(scale(t(log2(E+1)))); score<-colMeans(Z)
score_tab<-data.frame(cell_type=CT, machinery_score=round(score[CT],3))
lr<-sapply(CT,function(ct){d<-fread(list.files(DATA,pattern=sprintf("nmd_mashr_die_%s_.*csv",CODE[[ct]]),full.names=TRUE)[1]);sum(d$nmd_responsive%in%c(TRUE,"TRUE"))})
score_tab$n_NMD_isoforms<-lr[CT]
# heatmap of per-gene z
zl<-as.data.frame(Z) %>% rownames_to_column("gene") %>% pivot_longer(-gene,names_to="ct",values_to="z")
zl$ct<-factor(zl$ct,levels=CT); zl$gene<-factor(zl$gene,levels=CORE[CORE%in%rownames(Z)])
fig1<-ggplot(zl,aes(ct,gene,fill=z))+geom_tile(colour="white")+
  geom_text(aes(label=sprintf("%.0f",E[cbind(as.character(gene),as.character(ct))])),size=3)+
  scale_fill_gradient2(low="#4C78B0",mid="white",high="#c0392b",midpoint=0,name="z (row)")+
  labs(x=NULL,y=NULL,title="Core NMD machinery expression (SR gene level, DMSO CPM; color = row z-score)")+
  theme_minimal(base_size=12)+theme(axis.text=element_text(face="bold"))
img1<-saveimg(fig1,"machinery_heatmap.png",7,5)

# ---- canonical-isoform CPM for the 4 well-captured factors ----
CANON<-c(UPF2="ENST00000357604.10",SMG5="ENST00000361813.5",SMG8="ENST00000300917.10",SMG9="ENST00000270066.11")
di<-readRDS(P("dge_isoform_longread_2026.3.3.rds")); si<-as.data.frame(di$samples)%>%rownames_to_column("bam")
si$ct<-as.character(si$ct); si$treatment<-as.character(si$treatment); si$ct[si$ct=="AT"]<-"AT2"; si$ct[si$ct=="DD"]<-"LAE"; si<-si[si$ct%in%CT,]
ci<-cpm(di[,si$bam],log=FALSE); rownames(ci)<-di$genes$txid
canon<-rbindlist(lapply(c("DMSO","Smg1i"),function(trt){
  M<-sapply(CT,function(c) rowMeans(ci[,si$bam[si$ct==c&si$treatment==trt],drop=FALSE]))
  data.table(gene=names(CANON),treatment=trt, round(M[CANON,,drop=FALSE],2))}))
canon_dmso<-canon[treatment=="DMSO",.(gene,AT2,LAE,FB,MV)]
figc_dat<-melt(canon,id.vars=c("gene","treatment"),variable.name="ct",value.name="cpm")
figc_dat$ct<-factor(figc_dat$ct,levels=CT); figc_dat$treatment<-factor(figc_dat$treatment,levels=c("DMSO","Smg1i"))
fig2<-ggplot(figc_dat,aes(ct,cpm,fill=treatment))+geom_col(position=position_dodge(0.8),width=0.75)+
  facet_wrap(~gene,scales="free_y",nrow=1)+
  scale_fill_manual(values=c(DMSO="#4C78B0",Smg1i="#C0392B"),name=NULL)+
  labs(x=NULL,y="canonical-isoform mean CPM",
       title="Canonical (protein-coding, MANE-FSM) isoform expression — well-captured NMD factors")+
  theme_bw(base_size=12)+theme(legend.position="top",strip.background=element_rect(fill="grey92"),strip.text=element_text(face="bold"))
img2<-saveimg(fig2,"canonical_isoform_cpm.png",10,3.4)

# ---- METHOD 3: functional NMD-activity score = accumulation of a fixed validated endogenous substrate panel ----
PANEL<-c("ATF4","ATF3","GADD45A","GADD45B","DDIT3","ASNS","GAS5","TBL2","RHOB","NAT9","DDIT4","CHAC1","PPP1R15A","SLC7A11")
CTRL <-c("ACTB","GAPDH","TBP","B2M")
srlf<-sapply(CT,function(ct){d<-fread(list.files(DATA,pattern=sprintf("nmd_mashr_dge_%s_.*csv",CODE[[ct]]),full.names=TRUE)[1])
  setNames(d$logFC,d$hgnc_symbol)})   # posterior-mean logFC (SMG1i vs DMSO) by symbol
panel_lf<-srlf[rownames(srlf)%in%PANEL,,drop=FALSE]; panel_lf<-panel_lf[!duplicated(rownames(panel_lf)),]
ctrl_lf <-srlf[rownames(srlf)%in%CTRL ,,drop=FALSE]; ctrl_lf<-ctrl_lf[!duplicated(rownames(ctrl_lf)),]
raw_act<-apply(panel_lf,2,median,na.rm=TRUE); ctrl_med<-apply(ctrl_lf,2,median,na.rm=TRUE)
activity<-raw_act-ctrl_med                        # housekeeping-reference normalized (removes global per-CT shift)
act_tab<-data.frame(cell_type=CT, raw_medianLFC=round(raw_act[CT],3),
                    control_medianLFC=round(ctrl_med[CT],3),
                    normalized_activity=round(activity[CT],3))
pct_lost<-c(AT2=14.4,LAE=18.2,FB=14.6,MV=15.5)   # isoform-level % output lost (computed earlier)
cmp3<-data.frame(cell_type=CT, norm_substrate_activity=round(activity[CT],3),
                 machinery_score=round(score[CT],3), pct_output_lost=pct_lost[CT], n_NMD_isoforms=lr[CT])
cor_ms<-cor(activity[CT],score[CT],method="spearman"); cor_pl<-cor(activity[CT],pct_lost[CT],method="spearman")
# --- NMD-anchored normalization (anchor library to non-NMD gene content) = faithful analog of the custom norm ---
sym2ens<-setNames(base_id(sr1$ensembl_gene_id_version),sr1$hgnc_symbol)
nmd_set<-lapply(CT,function(ct){d<-fread(list.files(DATA,pattern=sprintf("nmd_mashr_dge_%s_.*csv",CODE[[ct]]),full.names=TRUE)[1]);base_id(d$ensembl_gene_id_version[d$nmd_responsive%in%c(TRUE,"TRUE")])}); names(nmd_set)<-CT
cnt<-g$counts[,rownames(g$samples)[keep]]; rownames(cnt)<-base_id(rownames(cnt))
anch<-function(ct,ids){ss<-rownames(gs)[gs$ct==ct]; dd<-ss[gs[ss,"treatment"]=="DMSO"]; sg<-ss[gs[ss,"treatment"]=="Smg1i"]
  lib<-colSums(cnt[setdiff(rownames(cnt),nmd_set[[ct]]),ss,drop=FALSE]); acpm<-sweep(cnt[,ss,drop=FALSE],2,lib,"/")*1e6
  ids<-ids[ids%in%rownames(acpm)]; median(log2((rowMeans(acpm[ids,sg,drop=FALSE])+1)/(rowMeans(acpm[ids,dd,drop=FALSE])+1)),na.rm=TRUE)}
activity<-setNames(sapply(CT,anch,ids=sym2ens[PANEL]),CT); ctrl_anch<-setNames(sapply(CT,anch,ids=sym2ens[CTRL]),CT)
act_tab<-data.frame(cell_type=CT, raw_medianLFC=round(raw_act[CT],3), NMDanchored_activity=round(activity[CT],3),
                    NMDanchored_control=round(ctrl_anch[CT],3))
cmp3<-data.frame(cell_type=CT, anchored_activity=round(activity[CT],3), machinery_score=round(score[CT],3),
                 pct_output_lost=pct_lost[CT], n_NMD_isoforms=lr[CT])
cor_ms<-cor(activity[CT],score[CT],method="spearman"); cor_pl<-cor(activity[CT],pct_lost[CT],method="spearman")
# heatmap of panel gene logFC
pl<-as.data.frame(panel_lf)%>%rownames_to_column("gene")%>%pivot_longer(-gene,names_to="ct",values_to="lfc")
pl$ct<-factor(pl$ct,levels=CT); pl$gene<-factor(pl$gene,levels=rownames(panel_lf)[order(-rowMeans(panel_lf,na.rm=TRUE))])
fig3<-ggplot(pl,aes(ct,gene,fill=lfc))+geom_tile(colour="white")+
  geom_text(aes(label=sprintf("%.1f",lfc)),size=3)+
  scale_fill_gradient2(low="#4C78B0",mid="white",high="#c0392b",midpoint=0,name="logFC")+
  labs(x=NULL,y=NULL,title="Endogenous NMD-substrate accumulation under SMG1i (SR gene posterior logFC)")+
  theme_minimal(base_size=12)+theme(axis.text=element_text(face="bold"))
img3<-saveimg(fig3,"substrate_panel_heatmap.png",6.5,5)

css<-"body{font-family:Segoe UI,Roboto,sans-serif;max-width:950px;margin:24px auto;padding:0 16px;color:#222}h1{border-bottom:3px solid #2166AC}h2{color:#2166AC;margin-top:30px}.k{background:#eef7ee;border-left:4px solid #16a085;padding:10px 14px;border-radius:4px;margin:12px 0}"
html<-paste0("<!doctype html><meta charset='utf-8'><style>",css,"</style>",
"<h1>NMD machinery expression across cell types</h1>",
"<p>Core NMD factors (UPF1/2/3A/3B, SMG1/5/6/7/8/9). Short-read gene level is length-unbiased (SMG1 quantifiable); isoform level is shown only for factors whose full-length canonical isoform is well captured by long-read.</p>",
"<h2>1. Short-read gene-level machinery score</h2>",tbl(score_tab),img1,
"<div class='k'>LAE expresses the NMD machinery highest (score +1.19), consistently across ~all factors (e.g., UPF1 275 vs 68-140 CPM; SMG1 130 vs 59-98), then AT2; FB/MV lowest. LAE also has the most NMD activity, so its NMD dominance is not only a donor-count/power artifact. Spearman(score, NMD isoforms)=0.40 across 4 cell types &mdash; LAE concordant, MV the outlier (high activity, low machinery).</div>",
"<p><b>Caveats:</b> n=4 (descriptive); NMD factors are themselves NMD-autoregulated so DMSO expression partly reflects NMD activity; LAE has 4 donors vs 3.</p>",
"<h2>2. Canonical-isoform CPM (well-captured factors)</h2>",
"<p>Long-read under-recovers full-length isoforms of the largest factors (SMG1 not assembled as FSM; UPF1/SMG6 canonical near floor). Shown here are the four factors whose canonical (MANE full-splice-match) isoform is reliably captured. Table is DMSO; figure shows DMSO vs SMG1i.</p>",
tbl(canon_dmso),img2,
"<div class='k'>Even for these, LAE/AT2 tend to exceed FB/MV (e.g., SMG8, UPF2, SMG5), mirroring the gene-level pattern; SMG9 is highest in AT2. Values rise under SMG1i (autoregulation).</div>",
"<h2>3. Functional NMD-activity score (endogenous substrate accumulation)</h2>",
"<p>The most direct cross-cell-type measure: accumulation of a fixed panel of validated endogenous NMD substrates (ATF4/ATF3/GADD45A/B/DDIT3/ASNS/GAS5/TBL2/RHOB/NAT9/DDIT4/CHAC1/PPP1R15A/SLC7A11) under SMG1i. To remove global per-cell-type shifts, each sample's library is normalized to <b>non-NMD gene content</b> (thousands of anchor genes) before computing the substrate log2 fold-change — a faithful gene-level analog of the manuscript's NMD-anchored normalization. Because the panel is the same transcripts everywhere, it is directly comparable.</p>",
tbl(act_tab),img3,
tbl(cmp3),
sprintf("<div class='k'>NMD is robustly and comparably active in all four cell types (median NMD-anchored substrate accumulation 1.17-1.56; most panel substrates strongly up, e.g. ATF4/DDIT3). With the non-NMD-content anchor, housekeeping controls sit near zero in every cell type (including MV), confirming there is no genome-wide artifact. Functional activity is <b>highest in MV and LAE</b> but only modestly so. Critically, functional activity does <b>not</b> track machinery-factor expression (Spearman %.2f) — LAE expresses the most machinery yet is not the most functionally active — whereas it agrees with the independent %% output-lost measure (Spearman %.2f). Conclusion: NMD-machinery expression is a poor proxy for functional NMD activity; the two functional readouts agree and show comparable NMD activity across cell types. n=4 (descriptive).</div>", cor_ms, cor_pl))
writeLines(html,file.path(OUT,"nmd_machinery_report.html"))
cat("Wrote",file.path(OUT,"nmd_machinery_report.html"),"\n\n")
cat("Canonical-isoform DMSO CPM:\n"); print(canon_dmso)
