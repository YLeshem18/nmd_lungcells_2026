#!/usr/bin/env Rscript
## =============================================================================
## build_manuscript_tables.R
##   Assemble the manuscript supplementary result workbooks (run on Channing).
##   Reads the CSV/TSV outputs already written by the analysis Rmds (see the
##   "made by" note on each input below) and packages them into six .xlsx files,
##   one worksheet per cell type where relevant.
##
##   Nothing here recomputes primary results — it only collects and formats what
##   the upstream Rmds produced. Each section fails soft: if an input is missing
##   the workbook is skipped with a message rather than erroring the whole run.
##
##   Output: BASE/tmp/manuscript_tables/*.xlsx
##   Cell types: AT2, LAE (files coded "dd"), FB, MV.
## =============================================================================
suppressPackageStartupMessages({ library(openxlsx); library(data.table) })

BASE    <- "/udd/reyle/nmd_lungcells_2026"
OUT_DIR <- file.path(BASE, "tmp", "manuscript_tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CT_CODE <- c(AT2 = "at2", LAE = "dd", FB = "fb", MV = "mv")   # LAE files coded "dd"

## ---- input locations (all pinned) ----
## Each path holds outputs written by an upstream Rmd in this repo:
DGE_DIR       <- file.path(BASE, "shortread_dge", "mashr")                              # gene mashr DGE CSVs      <- Gene-Level_DGE_Summary_mashR.Rmd
DIE_DIR       <- file.path(BASE, "long_read",     "mashr")                              # isoform mashr DIE CSVs   <- Isoform-Level_DIE_Summary_p1.Rmd
GSEA_ISO      <- file.path(BASE, "tmp", "mashr_die_enrichment")                         # isoform GSEA (fgsea)     <- die_mashr_enrichment_part2_2026.3.10.Rmd
RBP_DIR       <- file.path(BASE, "long_read", "nmd_rbp_enrichment")                     # RBP census / targets     <- rbp_sr.Rmd
PROD_DIR      <- file.path(BASE, "long_read", "productive_response", "basics_ordered")  # productive-response CSVs <- productive_response.Rmd
GSEA_GENE_SIG <- file.path(BASE, "shortread_dge", "gsea_classified", "all_significant_pathways.csv")  #  gene GSEA <- Gene-Level_DGE_Summary_mashR.Rmd
CENSUS        <- file.path(BASE, "New_NMD_Files", "gerstberger_2014_rbp_census.csv")    # static reference: Gerstberger et al. 2014 RBP census (New_NMD_Files/)

## ---- helpers ----
newest <- function(dir, pat){ f <- list.files(dir, pat, full.names = TRUE)
  if (!length(f)) return(NA_character_); f[which.max(file.mtime(f))] }
rd <- function(path) if (is.na(path) || !file.exists(path)) NULL else as.data.frame(fread(path))
find_either <- function(pat){ f <- c(newest(PROD_DIR, pat), newest(OUT_DIR, pat)); f <- f[!is.na(f)]
  if (length(f)) f[which.max(file.mtime(f))] else NA_character_ }
save_tabs <- function(tabs, file){
  tabs <- tabs[!vapply(tabs, is.null, logical(1))]
  if (!length(tabs)) { message("  [skip] ", file, " — no inputs found"); return(invisible()) }
  wb <- createWorkbook()
  for (nm in names(tabs)) { sn <- substr(gsub("[\\[\\]:*?/\\\\]", "_", nm), 1, 31)
    addWorksheet(wb, sn); writeData(wb, sn, tabs[[nm]]); freezePane(wb, sn, firstRow = TRUE) }
  saveWorkbook(wb, file.path(OUT_DIR, file), overwrite = TRUE)
  cat("  wrote", file, "(", length(tabs), "tabs )\n")
}
split_by_ct <- function(df, keep = c("pathway","collection","pval","padj","ES","NES","size","leadingEdge")){
  if (is.null(df)) return(NULL)
  ctcol <- intersect(c("cell_type","celltype","ct"), names(df))[1]
  if (is.na(ctcol)) return(list(all = df))
  cols <- intersect(keep, names(df))
  cts  <- unique(df[[ctcol]])
  setNames(lapply(cts, function(c){ s <- df[df[[ctcol]] == c, , drop = FALSE]
    s[order(s[["padj"]]), cols, drop = FALSE] }), cts)
}
sv <- function(x) sub("\\..*$", "", x)   # strip Ensembl version suffix

## =============================================================================
## 1. Gene-level pathways (significant only), tab per cell type
##    Source: Gene-Level_DGE_Summary_mashR.Rmd -> shortread_dge/gsea_classified/
##            all_significant_pathways.csv (fgsea on mashr gene effects, classified).
## =============================================================================
cat("[1] gene pathways\n")
save_tabs(split_by_ct(rd(GSEA_GENE_SIG)), "1_gene_pathways_by_celltype.xlsx")

## =============================================================================
## 2. Isoform-level pathways (significant only: NES>0 & padj<0.05), tab per cell type
##    Source: die_mashr_enrichment_part2_2026.3.10.Rmd -> tmp/mashr_die_enrichment/
##            mashr_gsea_largest_effect*.tsv (fgsea on the largest-effect isoform per gene).
## =============================================================================
cat("[2] isoform pathways\n")
ig <- rd(newest(GSEA_ISO, "mashr_gsea_largest_effect.*\\.tsv$"))
if (!is.null(ig)) ig <- ig[!is.na(ig$padj) & ig$padj < 0.05 & ig$NES > 0, ]
save_tabs(split_by_ct(ig), "2_isoform_pathways_by_celltype.xlsx")

## =============================================================================
## 3. RBP results (enrichment computed in-script) + DE-by-CT summary
##    Sources: nmd_target_rbps*.csv from rbp_sr.Rmd -> long_read/nmd_rbp_enrichment/;
##             DIE CSVs from Isoform-Level_DIE_Summary_p1.Rmd (long_read/mashr/);
##             DGE CSVs from Gene-Level_DGE_Summary_mashR.Rmd (shortread_dge/mashr/);
##             RBP census static ref (Gerstberger 2014, New_NMD_Files/).
##    The Fisher enrichment (NMD-target genes vs RBP census) is recomputed here
##    because rbp_sr.Rmd does not persist the enrichment table itself.
## =============================================================================
cat("[3] RBP + DE-by-CT\n")
census  <- rd(CENSUS)
rbp_ens <- if (!is.null(census)) unique(sv(census$ensembl_gene_id)) else character(0)
rbp_sym <- if (!is.null(census)) unique(census$gene_symbol)          else character(0)

die_list <- setNames(lapply(names(CT_CODE), function(ct)
  rd(newest(DIE_DIR, sprintf("nmd_mashr_die_%s_.*csv", CT_CODE[[ct]])))), names(CT_CODE))
die_list <- die_list[!vapply(die_list, is.null, logical(1))]

rbp_tabs <- list(nmd_target_rbps = rd(newest(RBP_DIR, "nmd_target_rbps.*csv")))
if (length(die_list) && !is.null(census)) {
  allrows <- rbindlist(lapply(die_list, function(d) data.table(
    gene_id = d$gene_id, sym = d$hgnc_symbol,
    nmd = d$nmd_responsive %in% c(TRUE, "TRUE"))), fill = TRUE)
  symmap <- allrows[, .(sym = sym[1]), by = gene_id]; symv <- setNames(symmap$sym, symmap$gene_id)
  bg    <- unique(allrows$gene_id)
  isrbp <- (sv(bg) %in% rbp_ens) | (symv[bg] %in% rbp_sym)
  fish  <- function(hit){ ft <- tryCatch(fisher.test(table(bg %in% hit, isrbp), alternative = "greater"),
                                          error = function(e) list(estimate = NA, p.value = NA))
    list(n = length(hit), n_rbp = sum(isrbp[bg %in% hit]),
         pct = round(100*mean(isrbp[bg %in% hit]), 1),
         OR = round(as.numeric(ft$estimate), 2), p = signif(ft$p.value, 3)) }
  o <- fish(unique(allrows$gene_id[allrows$nmd]))
  rbp_tabs$enrichment_overall <- data.frame(
    set = c("NMD-target", "background (DIE-tested)"),
    n_genes = c(o$n, length(bg)), n_RBP = c(o$n_rbp, sum(isrbp)),
    pct_RBP = c(o$pct, round(100*mean(isrbp), 1)), OR = c(o$OR, NA), p = c(o$p, NA))
  rbp_tabs$enrichment_by_ct <- rbindlist(lapply(names(die_list), function(ct){
    g <- unique(die_list[[ct]]$gene_id[die_list[[ct]]$nmd_responsive %in% c(TRUE, "TRUE")]); f <- fish(g)
    data.table(ct = ct, n_nmd = f$n, n_rbp = f$n_rbp, pct_rbp = f$pct, OR = f$OR, p = f$p) }))
} else message("  [note] census or DIE CSVs not found — enrichment tabs skipped")

de_by_ct <- rbindlist(lapply(names(CT_CODE), function(ct){
  d <- rd(newest(DGE_DIR, sprintf("nmd_mashr_dge_%s_.*csv", CT_CODE[[ct]]))); if (is.null(d)) return(NULL)
  data.table(cell_type = ct, n_tested = nrow(d),
             n_up_NMD = sum(d$adj.P.Val < 0.05 & d$logFC > 0, na.rm = TRUE),
             n_down   = sum(d$adj.P.Val < 0.05 & d$logFC < 0, na.rm = TRUE),
             n_nonsig = sum(d$adj.P.Val >= 0.05, na.rm = TRUE)) }), fill = TRUE)
rbp_tabs$DE_by_celltype <- as.data.frame(de_by_ct)
save_tabs(rbp_tabs, "3_rbp_and_DE_by_celltype.xlsx")

## =============================================================================
## 4. Productive results
##    Source: productive_response.Rmd -> long_read/productive_response/basics_ordered/
##            (mechanism-split and DOWN-gene topGO tables may land in tmp/manuscript_tables/;
##             find_either() searches both PROD_DIR and OUT_DIR).
## =============================================================================
cat("[4] productive results\n")
prod4 <- list()
# (a) mashr direction counts per CT (up/down/ns), by gene category
pdm <- rd(newest(PROD_DIR, "productive_direction_mashr_custom.*csv"))
if (!is.null(pdm)) { setDT(pdm); pdm <- pdm[gene_category %in% c("NMD","non-NMD")]
  prod4$mashr_direction_counts <- dcast(pdm, ct + gene_category ~ dir,
                                         fun.aggregate = length, value.var = "gene_id") }
# (b) Kitagawa level-vs-composition summary per CT (both % driven)
kd <- rd(newest(PROD_DIR, "kitagawa_decomposition.*csv"))
if (!is.null(kd)) { setDT(kd); kd <- kd[gene_category %in% c("NMD","non-NMD")]
  kd[, driver := ifelse(abs(composition) > abs(level), "composition", "level")]
  prod4$kitagawa_level_vs_comp <- kd[, .(median_level = round(median(level),2),
      median_composition = round(median(composition),2), median_dP = round(median(dP),2),
      pct_composition_driven = round(100*mean(driver=="composition"),1),
      pct_level_driven       = round(100*mean(driver=="level"),1)), by = .(ct, gene_category)] }
# (c) level-vs-composition by OUR metric (classify_mech: transcriptional/level vs composition), up+down
ms <- rbindlist(list(rd(find_either("productive_mech_split_down.*csv")),
                     rd(find_either("productive_mech_split_up.*csv"))), fill = TRUE)
if (nrow(ms)) prod4$mechanism_split <- ms
# (d) overall pathway enrichment (GO terms recurring in >=2 CTs)
prod4$pathways_overall <- rd(newest(PROD_DIR, "topgo_consensus.*csv"))
# (e) DOWN-gene GO-BP enrichment per mechanism (topGO / run_bp) — two separate tabs
prod4$down_pathways_transcriptional <- rd(find_either("down_pathways_transcriptional.*csv"))
prod4$down_pathways_composition     <- rd(find_either("down_pathways_composition.*csv"))
save_tabs(prod4, "4_productive_results.xlsx")

## =============================================================================
## 5. mashr DGE — gene level, tab per cell type
##    Source: Gene-Level_DGE_Summary_mashR.Rmd -> shortread_dge/mashr/nmd_mashr_dge_<ct>_*.csv
## =============================================================================
cat("[5] mashr DGE by CT\n")
save_tabs(setNames(lapply(names(CT_CODE), function(ct)
  rd(newest(DGE_DIR, sprintf("nmd_mashr_dge_%s_.*csv", CT_CODE[[ct]])))), names(CT_CODE)),
  "5_mashr_dge_by_celltype.xlsx")

## =============================================================================
## 6. mashr DIE — isoform level, tab per cell type
##    Source: Isoform-Level_DIE_Summary_p1.Rmd -> long_read/mashr/nmd_mashr_die_<ct>_*.csv
## =============================================================================
cat("[6] mashr DIE by CT\n")
save_tabs(setNames(lapply(names(CT_CODE), function(ct)
  rd(newest(DIE_DIR, sprintf("nmd_mashr_die_%s_.*csv", CT_CODE[[ct]])))), names(CT_CODE)),
  "6_mashr_die_by_celltype.xlsx")

cat("\nDone. Workbooks in:", OUT_DIR, "\n")
