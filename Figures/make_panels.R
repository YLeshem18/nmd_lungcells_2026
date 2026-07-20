#!/usr/bin/env Rscript
# =====================================================================
# make_panels.R — manuscript multipanel figure A-F (NMD long-read paper)
# Single patchwork composite, figures/lib conventions:
#   - theme_classic, NO panel border / frame
#   - two-font hierarchy, compact legends, all data shown
#   - panels are equal-sized cells (patchwork), tags A-F one size
#   A SR vs LR gene correlation | B isoform landscape sharing
#   C DIE volcano              | D productive logFC SR vs LR
#   E NMD-susceptible sharing  | F KEGG ER pathway (LAE)
#
# Usage: Rscript make_panels.R [DATA_DIR] [OUT_DIR] [REPO_DIR]
# =====================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(scales); library(patchwork); library(magick)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
options(bitmapType = "cairo")

args     <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(args) >= 1) args[1] else "nmd_fig_data"
OUT_DIR  <- if (length(args) >= 2) args[2] else "fig_panels"
REPO_DIR <- if (length(args) >= 3) args[3] else "repo_clone_test/nmd_isopair_analysis"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
P <- function(f) file.path(DATA_DIR, f)

CT_LEVELS <- c("AT2", "LAE", "FB", "MV")
LFSR_THR  <- 0.05
code2disp <- c(at = "AT2", at2 = "AT2", dd = "LAE", lae = "LAE", fb = "FB", mv = "MV")

# ---- theme: clean, no frame, two fonts ----
THEME <- theme_classic(base_size = 13) +
  theme(axis.title   = element_text(size = 14),
        axis.text    = element_text(size = 11),
        axis.line    = element_line(linewidth = 0.4),
        legend.text  = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.key.size = unit(0.40, "cm"),
        strip.background = element_blank(),
        strip.text   = element_text(face = "bold", size = 13),
        plot.title    = element_blank(),
        plot.subtitle = element_blank(),
        plot.margin   = margin(4, 6, 4, 6))

run <- function(fn, nm) tryCatch(fn(), error = function(e){ cat("PANEL", nm, "FAILED:", conditionMessage(e), "\n"); NULL })

# =====================================================================
# A — SR vs LR gene correlation (hexbin); thin colorbar legend
# =====================================================================
panelA <- function(){
  suppressPackageStartupMessages(library(edgeR))
  CT_OK <- c("AT","DD","FB","MV")
  lr_key <- function(nm){ f <- strsplit(nm,"_")[[1]]
    di <- which(grepl("^[0-9]{3}[A-Z]$", f))[1]; if (is.na(di)) return(NA_character_)
    ct <- paste(f[2:(di-1)], collapse="_"); if (!(ct %in% CT_OK)) return(NA_character_)
    paste(ct, f[di], f[length(f)], sep="_") }
  agg <- function(mat, keys){ ok <- !is.na(keys); t(rowsum(t(mat[, ok, drop=FALSE]), keys[ok])) }

  g <- readRDS(P("dge_gene_unfiltered_2026.1.2.rds"))
  sct <- as.character(g$samples$ct); sct[sct=="AT2"] <- "AT"
  skey <- paste(sct, g$samples$sample_id, g$samples$treatment, sep="_"); skey[!(sct %in% CT_OK)] <- NA
  cm <- g$counts; rownames(cm) <- sub("\\..*","",rownames(cm))
  s_short <- agg(cm, skey)

  ci <- fread(P("nmd_lungcells_filtered.count_matrix.txt"), data.table=FALSE)
  di <- readRDS(P("dge_isoform_longread_2026.3.3.rds"))
  tx2gene <- unique(di$genes[!is.na(di$genes$txid), c("txid","gene_id")])
  rownames(ci) <- ci$id; ci$id <- NULL; cmat <- round(as.matrix(ci))
  ix <- match(rownames(cmat), tx2gene$txid); has <- !is.na(ix)
  lg <- rowsum(cmat[has,], group = sub("\\..*","",tx2gene$gene_id[ix[has]]))
  s_long <- agg(lg, vapply(colnames(lg), lr_key, character(1)))

  cg <- intersect(rownames(s_short), rownames(s_long)); cs <- intersect(colnames(s_short), colnames(s_long))
  s <- s_short[cg,cs]; l <- s_long[cg,cs]; kp <- rowSums(s)>0 | rowSums(l)>0; s<-s[kp,]; l<-l[kp,]
  cs1 <- cpm(calcNormFactors(DGEList(s),"TMM"), normalized.lib.sizes=TRUE, log=TRUE)
  cl1 <- cpm(calcNormFactors(DGEList(l),"TMM"), normalized.lib.sizes=TRUE, log=TRUE)
  rp <- mean(vapply(seq_len(ncol(cs1)), function(i) cor(cs1[,i],cl1[,i]), numeric(1)))
  rs <- mean(vapply(seq_len(ncol(cs1)), function(i) cor(cs1[,i],cl1[,i],method="spearman"), numeric(1)))
  pd <- data.frame(short=as.vector(cs1), long=as.vector(cl1))
  ggplot(pd, aes(short, long)) +
    geom_hex(bins = 55) +
    scale_fill_gradientn(colors=c("white","lightblue","blue","darkblue","red"), trans="log10", name="Count") +
    annotate("text", x=min(pd$short), y=max(pd$long), hjust=0, vjust=1, size=3.6,
             label=sprintf("Spearman = %.3f\nPearson = %.3f", rs, rp)) +
    labs(x="Short-Read Gene log2(CPM)", y="Long-Read Gene log2(CPM)") +
    guides(fill = guide_colourbar(barwidth=unit(0.35,"cm"), barheight=unit(2.6,"cm"))) +
    THEME + theme(legend.position="right", panel.grid=element_blank())
}

# =====================================================================
# B — isoform landscape: shared vs cell-type-specific (split y-axis)
# =====================================================================
panelB <- function(){
  dge <- readRDS(P("dge_isoform_longread_filtered_2026.3.3.rds"))
  ct <- as.character(dge$samples$ct); ct[ct=="AT"]<-"AT2"; ct[ct=="DD"]<-"LAE"; dge$samples$ct <- ct
  keep <- ct %in% CT_LEVELS & as.character(dge$samples$treatment)=="DMSO"
  dge <- dge[, keep, keep.lib.sizes=FALSE]; samp <- dge$samples
  cts <- sort(unique(as.character(samp$ct)))
  ct_iso <- lapply(cts, function(c){ idx<-which(as.character(samp$ct)==c)
    rownames(dge$counts)[rowSums(dge$counts[,idx,drop=FALSE])>0] }); names(ct_iso) <- cts
  il <- rbindlist(lapply(names(ct_iso), function(c) data.table(iso=ct_iso[[c]], ct=c)))
  spec <- il[, .(n=uniqueN(ct)), by=iso][n==1, iso]
  us <- data.frame(ct=names(ct_iso), n_total=sapply(ct_iso,length),
                   n_specific=sapply(ct_iso,function(x) sum(x %in% spec))) %>%
    mutate(n_shared=n_total-n_specific) %>%
    pivot_longer(c(n_shared,n_specific), names_to="category", values_to="n") %>%
    mutate(category=ifelse(category=="n_specific","Cell-Type-Specific","Shared (2+)"))
  ord <- us %>% group_by(ct) %>% summarise(t=sum(n),.groups="drop") %>% arrange(desc(t)) %>% pull(ct)
  us$ct <- factor(us$ct, levels=ord)
  fv <- c("Cell-Type-Specific"="#F4A582","Shared (2+)"="#2166AC")
  ymax <- max(tapply(us$n, us$ct, sum)); brk <- floor(min(tapply(us$n[us$category=="Shared (2+)"], us$ct[us$category=="Shared (2+)"], sum)))*0.985
  top <- ggplot(us, aes(ct,n,fill=category)) + geom_col(alpha=.85) +
    geom_text(aes(label=comma(n)), position=position_stack(vjust=.5), size=3) +
    scale_fill_manual(values=fv) + scale_y_continuous(labels=comma) +
    coord_cartesian(ylim=c(brk,NA)) + labs(y=NULL,x=NULL,fill=NULL) +
    THEME + theme(legend.position="top", axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                  axis.line.x=element_blank(), panel.grid=element_blank())
  bot <- ggplot(us, aes(ct,n,fill=category)) + geom_col(alpha=.85) +
    scale_fill_manual(values=fv) + scale_y_continuous(labels=comma) +
    coord_cartesian(ylim=c(0, brk*0.45)) + labs(x="Cell Type", y=NULL, fill=NULL) +
    THEME + theme(legend.position="none", panel.grid=element_blank())
  (top / bot) + plot_layout(heights=c(3,1)) &
    theme(plot.margin=margin(2,4,2,4))
}

# =====================================================================
# C — DIE volcano; compact one-row legend (no clipping)
# =====================================================================
load_die <- function(){
  files <- list.files(DATA_DIR, pattern="^nmd_mashr_die_.*\\.csv$", full.names=TRUE)
  rbindlist(lapply(files, function(f){ code<-sub("^nmd_mashr_die_(.+)_2026.*","\\1",basename(f))
    d<-fread(f); d$ct<-code2disp[[code]]; d }), fill=TRUE)
}
panelC <- function(){
  d <- load_die(); d$ct <- factor(d$ct, levels=CT_LEVELS)
  d[, status := fifelse(logFC>0 & adj.P.Val<LFSR_THR, "NMD susceptible",
                fifelse(logFC<0 & adj.P.Val<LFSR_THR, "Other Significant", "Not significant"))]
  d$status <- factor(d$status, levels=c("NMD susceptible","Other Significant","Not significant"))
  ggplot(d, aes(logFC, -log10(adj.P.Val+1e-300), color=status)) +
    geom_point(alpha=.5, size=.5) + facet_wrap(~ct, ncol=2) +
    scale_color_manual(values=c("NMD susceptible"="steelblue","Other Significant"="coral","Not significant"="gray70")) +
    coord_cartesian(ylim=c(0,50)) +
    geom_hline(yintercept=-log10(LFSR_THR), linetype=2, color="black") +
    geom_vline(xintercept=0, linetype=2, color="black") +
    labs(x="Posterior Mean log2FC (SMG1i Vs DMSO)", y=expression(-log[10]("lfsr")), color=NULL) +
    guides(color=guide_legend(override.aes=list(size=4, alpha=1), nrow=1)) +
    THEME + theme(legend.position="top", legend.text=element_text(size=11),
                  legend.spacing.x=unit(0.2,"cm"))
}

# =====================================================================
# D & E — SR gene vs LR isoform
# =====================================================================
load_sr <- function(){
  files <- list.files(DATA_DIR, pattern="^nmd_mashr_dge_.*\\.csv$", full.names=TRUE)
  rbindlist(lapply(files, function(f){ code<-sub("^nmd_mashr_dge_(.+)_2026.*","\\1",basename(f))
    d<-fread(f); d$ct<-code2disp[[code]]; setnames(d,"ensembl_gene_id_version","gene_id",skip_absent=TRUE)
    d[,.(ct,gene_id,logFC,adj.P.Val)] }), fill=TRUE)
}
panelDE <- function(){
  sr<-load_sr(); die<-load_die()
  sr[, gc:=sub("\\..*","",gene_id)]; die[, gc:=sub("\\..*","",gene_id)]
  EPS<-1e-3
  lrg <- die[, .(logFC=sum(pmax(-log10(adj.P.Val),EPS)*logFC)/sum(pmax(-log10(adj.P.Val),EPS)),
                 lfsr=min(adj.P.Val,na.rm=TRUE)), by=.(ct,gc)]
  comb <- rbind(sr[,.(analysis="Short-Read Gene",ct,gc,logFC,lfsr=adj.P.Val)],
                lrg[,.(analysis="Long-Read Isoform",ct,gc,logFC,lfsr)])
  comb[,`:=`(analysis=factor(analysis,c("Short-Read Gene","Long-Read Isoform")),
             ct=factor(ct,CT_LEVELS), is_nmd=lfsr<LFSR_THR & logFC>0)]
  nmd<-comb[is_nmd==TRUE]; om<-median(nmd$logFC,na.rm=TRUE)
  pD <- ggplot(nmd, aes(ct,logFC,fill=analysis)) +
    geom_hline(yintercept=0,color="grey85") +
    geom_boxplot(alpha=.85, outlier.size=.4, outlier.alpha=.4, position=position_dodge(.75), width=.7) +
    scale_fill_manual(values=c("Short-Read Gene"="#7B8DC9","Long-Read Isoform"="#7BB07B")) +
    labs(x="Cell Type", y="Posterior Mean logFC (SMG1i Vs DMSO)", fill=NULL) +
    coord_cartesian(ylim=c(0, quantile(nmd$logFC,.99,na.rm=TRUE))) +
    THEME + theme(legend.position="top", panel.grid=element_blank())

  srf <- sr[adj.P.Val<LFSR_THR & logFC>0, .(analysis="Short-Read Gene", ct, feature=gc)]
  lrf <- die[adj.P.Val<LFSR_THR & logFC>0, .(analysis="Long-Read Isoform", ct, feature=txid)]
  ft <- rbind(srf,lrf); ft[, analysis:=factor(analysis,c("Short-Read Gene","Long-Read Isoform"))]
  sh <- ft[, .(n_cts=uniqueN(ct)), by=.(analysis,feature)]
  shl <- merge(ft, sh, by=c("analysis","feature")); shl[, bin:=factor(as.character(n_cts), levels=as.character(1:4))]
  ss <- shl[, .(n=.N), by=.(analysis,ct,bin)]; ss[, ct:=factor(ct,CT_LEVELS)]
  tot <- ss[, .(n_total=sum(n), pct=round(100*sum(n[bin=="1"])/sum(n),1)), by=.(analysis,ct)]
  fp <- c("1"="#D6604D","2"="#92C5DE","3"="#4393C3","4"="#2166AC")
  pE <- ggplot(ss, aes(ct,n,fill=bin)) + geom_col(alpha=.9,width=.7) +
    geom_text(data=tot, aes(ct,n_total,label=paste0(comma(n_total),"\n(",pct,"%)"),fill=NULL),
              vjust=-.2, size=2.7, lineheight=.85) +
    facet_wrap(~analysis, ncol=2, scales="free_y") +
    scale_fill_manual(values=fp, name="Shared in # CTs", labels=c("1"="1 (unique)","2"="2","3"="3","4"="4 (all)")) +
    scale_y_continuous(labels=comma, expand=expansion(mult=c(0,.18))) +
    labs(x="Cell Type", y="Number Of Features") +
    guides(fill=guide_legend(nrow=1)) +
    THEME + theme(legend.position="top", panel.grid=element_blank())
  list(D=pD, E=pE)
}

# =====================================================================
# F — KEGG ER pathway (LAE): pad the frameless pathview PNG to the 6x4 cell
#     directly with magick (NO ggplot, so NO panel border around it).
# =====================================================================
writeF <- function(out){
  f <- file.path(OUT_DIR, "hsa04141.LAE.png")
  if (!file.exists(f)) { cat("  [F] hsa04141.LAE.png missing\n"); return(invisible()) }
  W <- as.integer(CELL_W*300); H <- as.integer(CELL_H*300)   # 6x4 in @300dpi
  img <- magick::image_read(f)
  inf <- magick::image_info(img); b <- 2L                    # KEGG bakes a 1px black frame
  img <- magick::image_crop(img, sprintf("%dx%d+%d+%d", inf$width-2*b, inf$height-2*b, b, b))
  # label the KEGG color key (top-right) as logFC
  img <- magick::image_annotate(img, "logFC", size=30, weight=700, color="black",
                                gravity="northwest", location="+900+12")
  img <- magick::image_resize(img, paste0(W, "x", H))                 # fit, keep aspect
  img <- magick::image_extent(img, paste0(W, "x", H), color="white", gravity="center")
  magick::image_write(img, out)
  cat("  wrote", basename(out), "(frameless)\n")
}

# ---- build all panels, save each at the multipanel system's uniform cell ----
# (figures/multipanel convention: CELL_W x CELL_H = 6.0 x 4.0 in, NO tight bbox,
#  so the composite GridSpec(3,2) has identical cells with zero distortion.)
CELL_W <- 6.0; CELL_H <- 4.0
pA<-run(panelA,"A"); pB<-run(panelB,"B"); pC<-run(panelC,"C")
de<-run(panelDE,"D/E"); pD<-de$D; pE<-de$E
blank <- function() ggplot()+theme_void()
panels <- list(A=pA, B=pB, C=pC, D=pD, E=pE)          # A-E are ggplots
files  <- c(A="panelA_correlation.png", B="panelB_landscape.png",
            C="panelC_volcano.png",    D="panelD_logfc_box.png",
            E="panelE_sharing.png")
for (k in names(panels)) {
  p <- panels[[k]]; if (is.null(p)) p <- blank()
  ggsave(file.path(OUT_DIR, files[k]), p, width = CELL_W, height = CELL_H, dpi = 300, bg = "white")
  cat("  wrote", files[k], "\n")
}
writeF(file.path(OUT_DIR, "panelF_pathway.png"))      # F via magick, no frame
cat("Panels saved at", CELL_W, "x", CELL_H, "in. Compose with figure_composite.py\n")
