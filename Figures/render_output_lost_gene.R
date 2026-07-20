# =============================================================================
# Supplemental figure: % transcriptional output lost to NMD — GENE level
# Source Rmd : transcriptional_output.Rmd
# Chunk      : gene-level-barplot  (plot p_gene_global)
# Upstream deps inlined: libraries, load, matched-pairs, normalize,
#                        gene-level-output, gene-level-global-metric
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(edgeR)
  library(tibble)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
count  <- dplyr::count;  slice  <- dplyr::slice;  first  <- dplyr::first

DATA <- "nmd_fig_data"
P    <- function(f) file.path(DATA, f)
HERE <- "supplement_figures"
dir.create(HERE, recursive = TRUE, showWarnings = FALSE)

DGELIST_PATH <- P("dge_isoform_longread_2026.3.3.rds")
CT_KEEP <- c("LAE", "AT2", "FB", "MV")

# ---- Load DGEList and parse metadata -------------------------------------
dge        <- readRDS(DGELIST_PATH)
samp <- dge$samples %>%
  rownames_to_column("bam_id") %>%
  mutate(ct = as.character(ct), treatment = as.character(treatment), donor_id = id)
samp$ct[samp$ct == "AT"] <- "AT2"
samp$ct[samp$ct == "DD"] <- "LAE"
samp <- samp %>% filter(ct %in% CT_KEEP)

iso_info <- dge$genes %>% select(any_of(c("txid", "gene_id", "hgnc_symbol")))

# ---- Matched donor pairs -------------------------------------------------
paired_donors <- samp %>%
  group_by(ct, donor_id) %>%
  summarise(has_smg1i = any(treatment == "Smg1i"),
            has_dmso  = any(treatment == "DMSO"), .groups = "drop") %>%
  filter(has_smg1i & has_dmso)

# ---- Plain CPM (no TMM) --------------------------------------------------
cpm_mat <- cpm(dge[, samp$bam_id], normalized.lib.sizes = FALSE, log = FALSE)

# ---- Gene-level aggregation + delta --------------------------------------
iso_to_gene <- iso_info %>%
  select(txid, gene_id, hgnc_symbol) %>%
  filter(!is.na(gene_id) & gene_id != "")

keep_iso <- intersect(rownames(cpm_mat), iso_to_gene$txid)
cpm_gene <- cpm_mat[keep_iso, samp$bam_id]
gene_ids <- iso_to_gene$gene_id[match(keep_iso, iso_to_gene$txid)]
stopifnot(length(gene_ids) == nrow(cpm_gene), !any(is.na(gene_ids)))
cpm_gene_agg <- rowsum(cpm_gene, group = gene_ids)

gene_anno <- iso_to_gene %>%
  group_by(gene_id) %>%
  summarise(hgnc_symbol = first(hgnc_symbol), n_isoforms = n(), .groups = "drop")

gene_results_list <- list()
for (ct_i in sort(unique(paired_donors$ct))) {
  donors_i <- paired_donors %>% filter(ct == ct_i) %>% pull(donor_id)
  for (d in donors_i) {
    smg1i_id <- samp %>% filter(ct == ct_i, donor_id == d, treatment == "Smg1i") %>% pull(bam_id)
    dmso_id  <- samp %>% filter(ct == ct_i, donor_id == d, treatment == "DMSO")  %>% pull(bam_id)
    stopifnot(length(smg1i_id) == 1, length(dmso_id) == 1)
    val_smg1i <- cpm_gene_agg[, smg1i_id]
    val_dmso  <- cpm_gene_agg[, dmso_id]
    delta     <- val_smg1i - val_dmso
    delta_pos <- pmax(delta, 0)
    frac_lost <- ifelse(val_smg1i > 0, delta_pos / val_smg1i, 0)
    gene_results_list[[paste(ct_i, d, sep = "_")]] <- data.frame(
      gene_id = rownames(cpm_gene_agg), ct = ct_i, donor_id = d,
      smg1i_val = val_smg1i, dmso_val = val_dmso, delta = delta,
      delta_pos = delta_pos, frac_lost = frac_lost, stringsAsFactors = FALSE)
  }
}
gene_results_all <- bind_rows(gene_results_list) %>% left_join(gene_anno, by = "gene_id")

# ---- Global gene-level metric --------------------------------------------
gene_global_metric <- gene_results_all %>%
  group_by(ct) %>%
  summarise(sum_delta_pos   = sum(delta_pos),
            total_smg1i     = sum(smg1i_val),
            pct_output_lost = 100 * sum_delta_pos / total_smg1i,
            n_gene_donor    = n(), .groups = "drop") %>%
  arrange(desc(pct_output_lost))

# ---- Figure --------------------------------------------------------------
p_gene_global <- ggplot(gene_global_metric,
                        aes(x = reorder(ct, pct_output_lost),
                            y = pct_output_lost, fill = ct)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0(round(pct_output_lost, 2), "%")),
            hjust = -0.1, size = 4) +
  coord_flip() +
  labs(title    = "% Transcriptional Output Lost To NMD By Cell Type",
       subtitle = "Gene-level",
       x = "Cell Type", y = "% Output Lost To NMD") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(color = "grey40")) +
  expand_limits(y = max(gene_global_metric$pct_output_lost) * 1.15)

ggsave(file.path(HERE, "pct_output_lost_gene_level_barplot.png"),
       p_gene_global, width = 10, height = 6, dpi = 300, bg = "white")
