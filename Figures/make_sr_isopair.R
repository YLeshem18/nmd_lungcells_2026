#!/usr/bin/env Rscript
# SR-protein Isopair figures (rbp_sr.Rmd): per-gene non-NMD vs NMD isoform structure
# (5'UTR/CDS/3'UTR, events, ref-AUG PTC marker) + CPM side panel AS GROUPED BARS.
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(edgeR); library(Biostrings); library(Isopair)
})
select<-dplyr::select; filter<-dplyr::filter
options(bitmapType="cairo")
DATA<-"nmd_fig_data"
HERE<-"supplement_figures"
dir.create(HERE, recursive=TRUE, showWarnings=FALSE)
P<-function(f) file.path(DATA,f)
CELL_TYPES<-c("AT2","LAE","FB","MV"); LFSR_NMD_THR<-0.05; sr_family<-paste0("SRSF",1:12)
sr_num<-function(x) as.integer(sub("^SRSF","",x))
ct_disp<-c(AT2="AT2",LAE="LAE",FB="FB",MV="MV")

dge_iso<-readRDS(P("dge_isoform_longread_filtered_2026.3.3.rds"))
pm_iso<-read.csv(P("mashr_isoform_posterior_means_2026.3.10.csv"),check.names=FALSE)
lfsr_iso<-read.csv(P("mashr_isoform_lfsr_2026.3.10.csv"),check.names=FALSE)
colnames(pm_iso)[1]<-"isoform_id"; colnames(lfsr_iso)[1]<-"isoform_id"
iso2gene<-data.frame(isoform_id=dge_iso$genes$txid, gene_id=dge_iso$genes$gene_id, hgnc_symbol=dge_iso$genes$hgnc_symbol)
g2sym<-distinct(iso2gene,gene_id,hgnc_symbol); sym_of<-function(g) g2sym$hgnc_symbol[match(g,g2sym$gene_id)]
resolve_mash_col<-function(df,ct){cn<-colnames(df); if(ct%in%cn) return(ct); h<-cn[grepl(paste0("_in_",ct,"$"),cn)]; h[1]}

samp<-dge_iso$samples; trt_col<-intersect(c("treatment","trt"),colnames(samp))[1]
is_smg<-grepl("SMG|Smg1i",samp[[trt_col]],ignore.case=TRUE)
cpm_iso<-cpm(dge_iso); mean_dmso<-rowMeans(cpm_iso[,!is_smg,drop=FALSE]); mean_smg<-rowMeans(cpm_iso[,is_smg,drop=FALSE])
ct_col<-intersect(c("ct","cell_type"),colnames(samp))[1]
samp_ct<-as.character(samp[[ct_col]]); samp_ct[samp_ct=="AT"]<-"AT2"; samp_ct[samp_ct=="DD"]<-"LAE"
expr_ct<-sapply(CELL_TYPES, function(c) rowMeans(cpm_iso[, samp_ct==c & is_smg, drop=FALSE]))

nmd_iso_all<-unique(unlist(lapply(CELL_TYPES,function(ct){
  l<-lfsr_iso[[resolve_mash_col(lfsr_iso,ct)]]; p<-pm_iso[[resolve_mash_col(pm_iso,ct)]]
  lfsr_iso$isoform_id[!is.na(l)&!is.na(p)&l<LFSR_NMD_THR&p>0]})))
sr_gene_ids<-iso2gene %>% filter(hgnc_symbol %in% sr_family) %>% pull(gene_id) %>% unique()
pick_top<-function(ids,score){ids<-ids[ids%in%names(score)]; ids<-ids[!is.na(score[ids])]; if(!length(ids)) return(NA_character_); ids[which.max(score[ids])]}
names(mean_dmso)<-rownames(cpm_iso); names(mean_smg)<-rownames(cpm_iso)

sr_pairs<-bind_rows(lapply(sr_gene_ids,function(g){
  iso<-intersect(iso2gene$isoform_id[iso2gene$gene_id==g], lfsr_iso$isoform_id)
  nmd<-iso[iso%in%nmd_iso_all]; non<-iso[!(iso%in%nmd_iso_all)]
  if(!length(nmd)||!length(non)) return(NULL)
  ref<-pick_top(non,mean_dmso); cnmd<-pick_top(nmd,mean_smg); non2<-setdiff(non,ref)
  cctl<-if(length(non2)) pick_top(non2,mean_dmso) else NA_character_
  out<-data.frame(gene_id=g,reference_isoform_id=ref,comparator_isoform_id=cnmd,pair_class="NMD")
  if(!is.na(cctl)) out<-bind_rows(out,data.frame(gene_id=g,reference_isoform_id=ref,comparator_isoform_id=cctl,pair_class="Control"))
  out})) %>% mutate(sym=sym_of(gene_id))
wanted<-unique(c(sr_pairs$reference_isoform_id,sr_pairs$comparator_isoform_id))

structures<-Isopair::parseIsoformStructures(P("srsf.gtf"))
ue<-Isopair::buildUnionExons(structures)
sr_prof<-Isopair::buildProfiles(sr_pairs,structures,ue$union_exons,ue$isoform_union_mapping,verify=TRUE,verbose=FALSE) %>%
  left_join(sr_pairs %>% select(reference_isoform_id,comparator_isoform_id,sym,pair_class),by=c("reference_isoform_id","comparator_isoform_id"))
poison_comp<-sr_prof %>% select(detailed_events) %>% tidyr::unnest(detailed_events) %>%
  filter(event_type=="SE",direction=="GAIN") %>% distinct(comparator_isoform_id) %>% pull(comparator_isoform_id)
sr_prof<-sr_prof %>% mutate(has_poison=comparator_isoform_id %in% poison_comp)

cds_meta<-Isopair::extractCdsAnnotations(P("srsf.cds.gff3"), isoform_ids=wanted, verbose=FALSE)
fai<-fasta.index(P("srsf.fasta"),seqtype="DNA"); fai$id<-sub("\\s.*","",fai$desc)
seqs<-readDNAStringSet(fai[fai$id %in% wanted,]); sequences<-setNames(toupper(as.character(seqs)),sub("\\s.*","",names(seqs)))
orfs<-Isopair::enumerateOrfs(structures,cds_meta,sequences)
ref_by_gene<-sr_pairs %>% distinct(gene_id,reference_isoform_id)
prim<-bind_rows(lapply(seq_len(nrow(ref_by_gene)),function(i){
  g<-ref_by_gene$gene_id[i]; ref<-ref_by_gene$reference_isoform_id[i]
  ids<-sr_pairs %>% filter(gene_id==g) %>% { unique(c(.$reference_isoform_id,.$comparator_isoform_id)) }
  Isopair::selectPrimaryOrf(orfs %>% filter(isoform_id%in%ids),structures %>% filter(isoform_id%in%ids),
                            cds_meta %>% filter(isoform_id%in%ids),dominant_isoform_id=ref)}))
tx2genomic<-function(id,tx_pos){r<-structures[structures$isoform_id==id,]; if(nrow(r)==0||is.na(tx_pos)) return(NA_integer_)
  s<-r$exon_starts[[1]]; e<-r$exon_ends[[1]]; strand<-r$strand; o<-order(s); s<-s[o]; e<-e[o]; lens<-e-s+1
  if(strand=="-"){s<-rev(s); e<-rev(e); lens<-rev(lens)}; cum<-cumsum(lens); prev<-c(0,head(cum,-1))
  ex<-which(tx_pos<=cum)[1]; if(is.na(ex)) return(NA_integer_); into<-tx_pos-prev[ex]-1
  if(strand=="-") e[ex]-into else s[ex]+into}
ptc_x<-prim %>% left_join(orfs %>% select(isoform_id,atg_tx_pos,stop_tx_pos),by=c("isoform_id","primary_atg_tx_pos"="atg_tx_pos")) %>%
  rowwise() %>% mutate(stop_x=tx2genomic(isoform_id,stop_tx_pos),atg_x=tx2genomic(isoform_id,primary_atg_tx_pos),
                       has_ptc=primary_n_downstream_ejc>=1) %>% ungroup()
ref_cds<-cds_meta %>% filter(isoform_id %in% sr_pairs$reference_isoform_id) %>% transmute(isoform_id,cds_start,cds_end=cds_stop,strand)
comp_cds<-ptc_x %>% filter(isoform_id %in% sr_pairs$comparator_isoform_id,!is.na(atg_x),!is.na(stop_x)) %>%
  mutate(strand=structures$strand[match(isoform_id,structures$isoform_id)]) %>%
  transmute(isoform_id,cds_start=pmin(atg_x,stop_x),cds_end=pmax(atg_x,stop_x),strand)
cds_for_plot<-bind_rows(ref_cds,comp_cds)

# ---- plotting helpers (your code) ----
relabel_tracks<-function(p){for(i in seq_along(p$layers)){d<-p$layers[[i]]$data
  if(is.data.frame(d)&&"label"%in%names(d)){d$label<-sub("^Reference","non-NMD",sub("^Comparator","NMD",d$label)); p$layers[[i]]$data<-d}}; p}
shrink_utr<-function(p,utr_h=0.10,cds_h=0.25){
  k<-which(vapply(p$layers,function(L) is.data.frame(L$data)&&all(c("region","y_pos","exon_start","exon_end")%in%names(L$data)),logical(1)))
  if(!length(k)) return(p); k<-k[1]; d<-p$layers[[k]]$data; h<-ifelse(d$region%in%c("5UTR","3UTR"),utr_h,cds_h)
  d$.ymin<-d$y_pos-h; d$.ymax<-d$y_pos+h
  p$layers[[k]]<-ggplot2::geom_rect(data=d,inherit.aes=FALSE,colour=NA,
    mapping=ggplot2::aes(xmin=exon_start,xmax=exon_end,ymin=.data$.ymin,ymax=.data$.ymax,fill=region)); p}

# CPM side panel AS GROUPED BARS (replaces the heatmap)
TRACKCOL<-c("non-NMD"="#4C78B0","NMD"="#C0392B")
mk_side<-function(r){
  ids<-c("non-NMD"=r$reference_isoform_id,"NMD"=r$comparator_isoform_id)
  d<-do.call(rbind,lapply(names(ids),function(lab)
    data.frame(track=lab,ct=CELL_TYPES,cpm=expr_ct[match(ids[[lab]],rownames(expr_ct)),CELL_TYPES])))
  d$ctd<-factor(unname(ct_disp[d$ct]),levels=CELL_TYPES); d$track<-factor(d$track,levels=c("non-NMD","NMD"))
  ggplot(d,aes(ctd,cpm,fill=track))+
    geom_col(position=position_dodge(width=0.78),width=0.72)+
    geom_text(aes(label=round(cpm)),position=position_dodge(width=0.78),vjust=-0.25,size=3.1)+
    scale_fill_manual(values=TRACKCOL,name=NULL)+
    scale_y_continuous(expand=expansion(mult=c(0,0.15)))+
    labs(x=NULL,y="mean CPM",subtitle="mean CPM (SMG1i)")+
    theme_minimal(base_size=13)+
    theme(legend.position="top",panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=11,face="bold"),plot.subtitle=element_text(size=13,hjust=0.5,face="bold"))
}
get_exons<-function(id){r<-structures[structures$isoform_id==id,]; data.frame(exon_start=r$exon_starts[[1]],exon_end=r$exon_ends[[1]])}
get_strand<-function(id) structures$strand[match(id,structures$isoform_id)]
build_one<-function(r){
  refx<-get_exons(r$reference_isoform_id); compx<-get_exons(r$comparator_isoform_id); strand<-get_strand(r$reference_isoform_id)
  p<-Isopair::plotIsoformPair(reference_exons=refx,comparator_exons=compx,
    events=if("detailed_events"%in%colnames(sr_prof)) r$detailed_events[[1]] else NULL,
    gene_id=r$gene_id,reference_id=r$reference_isoform_id,comparator_id=r$comparator_isoform_id,
    strand=strand,cds_metadata=cds_for_plot,show_events=TRUE,pair_label=r$sym)
  p<-shrink_utr(relabel_tracks(p))
  # italicize the gene symbol in the figure title (keep gene id roman)
  p<-p+ggplot2::labs(title=bquote(italic(.(r$sym))*" "*.(r$gene_id)))
  px<-ptc_x %>% filter(isoform_id==r$comparator_isoform_id, has_ptc %in% TRUE)
  if(nrow(px)) p<-p+ggplot2::annotate("text",x=px$stop_x,y=1.45,label="*",size=7,colour="red",fontface="bold")
  patchwork::wrap_plots(p,mk_side(r),widths=c(5,1.6))
}

nmd_prof<-sr_prof %>% filter(pair_class=="NMD") %>% arrange(sr_num(sym))
cat("SR NMD pairs:\n"); print(nmd_prof %>% select(sym,comparator_isoform_id,n_se,has_poison) %>% as.data.frame())
fwrite(sr_prof %>% filter(pair_class=="NMD") %>%
  left_join(ptc_x %>% select(isoform_id,has_ptc),by=c("comparator_isoform_id"="isoform_id")) %>%
  select(sym,comparator_isoform_id,n_se,has_poison,has_ptc) %>% arrange(sr_num(sym)),
  file.path(HERE,"data_sr_events.csv"))
for(i in seq_len(nrow(nmd_prof))){
  r<-nmd_prof[i,]
  if(!(r$reference_isoform_id%in%structures$isoform_id)||!(r$comparator_isoform_id%in%structures$isoform_id)){
    cat(sprintf("%s skipped (not in GTF)\n",r$sym)); next}
  g<-build_one(r)
  ggsave(file.path(HERE,sprintf("SR_%s.png",r$sym)), g, width=12, height=3.2, dpi=300, bg="white")
  cat(sprintf("wrote SR_%s.png\n",r$sym))
}
cat("Done. SR Isopair figures in figureS2_sr_isopair/.\n")
