#!/usr/bin/env Rscript
# Panel B (option: two opposing productive programs) — productive log2FC for
# GO-defined ISR/UPR (up) vs proliferation/DNA-repair (down) gene sets, per CT.
# Uses the custom NMD-anchored normalization (build_adj/process_ct) from
# productive_response.Rmd, on the data in nmd_fig_data/.
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(edgeR); library(org.Hs.eg.db); library(AnnotationDbi)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
options(bitmapType="cairo")
DATA_DIR <- "nmd_fig_data"; OUT_DIR <- "fig_panels"
P <- function(f) file.path(DATA_DIR, f)
CELL_TYPES <- c("AT2","LAE","FB","MV"); LFSR_NMD_THR <- 0.05
FLOOR_BOTH<-1; FLOOR_DMSO<-2
pass_floor <- function(td,ts) (td>=FLOOR_BOTH & ts>=FLOOR_BOTH) | (td>=FLOOR_DMSO)
svf <- function(x) sub("\\..*$","",x)

# ---- load ----
dge_iso  <- readRDS(P("dge_isoform_longread_filtered_2026.3.3.rds"))
pm_iso   <- read.csv(P("mashr_isoform_posterior_means_2026.3.10.csv"), check.names=FALSE)
lfsr_iso <- read.csv(P("mashr_isoform_lfsr_2026.3.10.csv"), check.names=FALSE)
colnames(pm_iso)[1] <- "isoform_id"; colnames(lfsr_iso)[1] <- "isoform_id"

iso2gene <- data.frame(isoform_id=dge_iso$genes$txid, gene_id=dge_iso$genes$gene_id,
                       hgnc_symbol=dge_iso$genes$hgnc_symbol, stringsAsFactors=FALSE)
n_iso_tab <- iso2gene %>% count(gene_id, name="n_iso")

samples <- as.data.frame(dge_iso$samples)
samples$cell_type <- as.character(samples$ct)
samples$cell_type[samples$cell_type=="AT"] <- "AT2"
samples$cell_type[samples$cell_type=="DD"] <- "LAE"
samples$treatment <- as.character(samples$treatment)
samples$sample_id <- as.character(samples$id)
samples <- samples[samples$cell_type %in% CELL_TYPES, ]
counts_iso <- dge_iso$counts[, rownames(samples)]

resolve <- function(df, ct){ cn<-colnames(df); if(ct %in% cn) return(ct)
  h<-cn[grepl(paste0("_in_",ct,"$"),cn)]; if(length(h)==1) return(h); stop("no col ",ct) }

# ---- custom NMD-anchored normalization (verbatim logic) ----
build_adj <- function(ct){
  l <- lfsr_iso[[resolve(lfsr_iso,ct)]]; p <- pm_iso[[resolve(pm_iso,ct)]]
  is_unprod <- !is.na(l)&!is.na(p)&l<LFSR_NMD_THR&p>0
  cls <- setNames(ifelse(is_unprod,"unprod","prod"), lfsr_iso$isoform_id)
  cs <- rownames(samples)[samples$cell_type==ct]; cnt <- counts_iso[, cs, drop=FALSE]
  trt <- setNames(samples[cs,"treatment"], cs); donor <- setNames(samples[cs,"sample_id"], cs)
  rawcpm <- sweep(cnt,2,colSums(cnt),"/")*1e6
  ids <- rownames(rawcpm); gene <- iso2gene$gene_id[match(ids, iso2gene$isoform_id)]
  class_v <- cls[ids]
  anchor <- !is.na(gene) & !grepl("::",gene,fixed=TRUE) & !is.na(class_v)
  dmso <- cs[trt[cs]=="DMSO"]; smg <- cs[trt[cs]=="Smg1i"]
  dmso_mean <- rowMeans(rawcpm[,dmso,drop=FALSE])
  prod_idx <- anchor & class_v=="prod"
  Pdmso <- tapply(dmso_mean[prod_idx], gene[prod_idx], sum)
  iu <- which(anchor & class_v=="unprod"); pg <- Pdmso[gene[iu]]
  hp <- !is.na(pg)&pg>0; share <- ifelse(hp, dmso_mean[iu]/pg, NA_real_)
  ia <- iu[hp]; ina <- iu[!hp]
  adj <- rawcpm
  for (s in smg){ col<-rawcpm[,s]; Pgs<-tapply(col[prod_idx],gene[prod_idx],sum); nc<-col
    nc[ia]<-Pgs[gene[ia]]*share[hp]; nc[ina]<-dmso_mean[ina]; nc[!is.finite(nc)]<-0
    adj[,s] <- nc * (1e6/sum(nc)) }
  list(rawcpm=rawcpm, adj=adj, gene=gene, class_v=class_v, anchor=anchor, trt=trt, donor=donor)
}
process_ct <- function(ct){
  o<-build_adj(ct); k<-o$anchor; fk<-paste(o$gene[k],o$class_v[k],sep="|||")
  rawf<-rowsum(o$rawcpm[k,,drop=FALSE],fk); adjf<-rowsum(o$adj[k,,drop=FALSE],fk)
  fm<-data.frame(feature=rownames(rawf)); fm$gene<-sub("\\|\\|\\|.*$","",fm$feature); fm$class<-sub("^.*\\|\\|\\|","",fm$feature)
  sm<-data.frame(sample=colnames(rawf),donor=o$donor[colnames(rawf)],trt=o$trt[colnames(rawf)])
  melt<-function(m,v){d<-as.data.frame(m);d$feature<-rownames(m);pivot_longer(d,-feature,names_to="sample",values_to=v)}
  long<-melt(adjf,"adj") %>% left_join(melt(rawf,"raw"),by=c("feature","sample")) %>%
    left_join(fm,by="feature") %>% left_join(sm,by="sample")
  dm<-long %>% filter(trt=="DMSO") %>% select(gene,class,donor,cpm=adj) %>%
    pivot_wider(names_from=class,values_from=cpm,values_fill=0) %>%
    rename(prodCPM_DMSO=prod, unprodCPM_DMSO=unprod)
  sg<-long %>% filter(trt=="Smg1i") %>% select(gene,class,donor,adj,raw) %>%
    pivot_wider(names_from=class,values_from=c(adj,raw),values_fill=0) %>%
    rename(prodCPM_SMG_adj=adj_prod, prodCPM_SMG_raw=raw_prod, unprodCPM_SMG=raw_unprod)
  nmd_g<-unique(o$gene[o$anchor & o$class_v=="unprod"])
  full_join(dm,sg,by=c("gene","donor")) %>%
    mutate(ct=ct, gene_category=ifelse(gene %in% nmd_g,"NMD","non-NMD"))
}
gene_sample <- bind_rows(lapply(CELL_TYPES, process_ct)) %>% rename(gene_id=gene)
base <- gene_sample %>%
  filter(!is.na(prodCPM_DMSO), !is.na(prodCPM_SMG_raw)) %>%
  mutate(across(c(prodCPM_DMSO,prodCPM_SMG_raw,prodCPM_SMG_adj,unprodCPM_DMSO,unprodCPM_SMG), ~coalesce(.x,0))) %>%
  group_by(ct, gene_id, gene_category) %>%
  summarise(prodCPM_DMSO=mean(prodCPM_DMSO), prodCPM_SMG_adj=mean(prodCPM_SMG_adj),
            unprodCPM_DMSO=mean(unprodCPM_DMSO), unprodCPM_SMG=mean(unprodCPM_SMG),
            prodCPM_SMG_raw=mean(prodCPM_SMG_raw), .groups="drop") %>%
  mutate(totalCPM_DMSO=prodCPM_DMSO+unprodCPM_DMSO, totalCPM_SMG=prodCPM_SMG_raw+unprodCPM_SMG) %>%
  filter(pass_floor(totalCPM_DMSO,totalCPM_SMG)) %>%
  mutate(custom = log2((prodCPM_SMG_adj+1)/(prodCPM_DMSO+1)))

# ---- GO programs ----
go_genes <- function(go_ids){
  e <- tryCatch(unique(unlist(AnnotationDbi::mget(go_ids, org.Hs.egGO2ALLEGS))), error=function(...) character(0))
  unique(na.omit(AnnotationDbi::mapIds(org.Hs.eg.db, keys=e, column="ENSEMBL", keytype="ENTREZID", multiVals="first")))
}
up_go <- go_genes(c("GO:0006986","GO:0140467"))                 # ISR/UPR
dn_go <- go_genes(c("GO:0006281","GO:0006260","GO:0000278"))    # DNA repair + replication + mitosis
both <- intersect(up_go,dn_go); up_go<-setdiff(up_go,both); dn_go<-setdiff(dn_go,both)

pg <- base %>% mutate(g0=svf(gene_id),
        program=case_when(g0 %in% up_go ~ "ISR/UPR (up)",
                          g0 %in% dn_go ~ "Proliferation /\nDNA-repair (down)", TRUE ~ NA_character_)) %>%
  filter(!is.na(program))
pg$ct <- factor(pg$ct, levels=CELL_TYPES)

cat("Program gene-CT counts:\n"); print(pg %>% count(ct, program) %>% pivot_wider(names_from=program, values_from=n))

pg$program <- factor(pg$program, levels=c("ISR/UPR (up)","Proliferation /\nDNA-repair (down)"))
THEME <- theme_classic(base_size=13) +
  theme(axis.title=element_text(size=14), axis.text.y=element_text(size=11),
        axis.line=element_line(linewidth=0.4),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position="top", legend.text=element_text(size=11), legend.title=element_blank(),
        legend.key.size=unit(0.4,"cm"),
        strip.background=element_blank(),
        strip.text=element_text(face="bold", size=13), panel.grid=element_blank(),
        plot.margin=margin(4,6,4,6))

pB <- ggplot(pg, aes(program, custom, fill=program)) +
  geom_hline(yintercept=0, colour="grey60") +
  geom_boxplot(outlier.size=0.3, outlier.alpha=0.25, width=0.6, colour="grey25") +
  facet_wrap(~ct, nrow=1) + coord_cartesian(ylim=c(-1.2,1.2)) +
  scale_fill_manual(values=c("ISR/UPR (up)"="#16a085","Proliferation /\nDNA-repair (down)"="#c0392b"),
                    labels=c("ISR/UPR (up)","Proliferation / DNA-repair (down)")) +
  labs(x=NULL, y="Productive log2FC (Smg1i Vs DMSO)") +
  THEME

ggsave(file.path(OUT_DIR,"panelB_twoprogram.png"), pB, width=6, height=4, dpi=300, bg="white")
cat("wrote panelB_twoprogram.png\n")
