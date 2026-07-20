# =============================================================================
# Supplemental figure: Per-isoform fraction of output lost to NMD (boxplot)
# Source Rmd : transcriptional_output.Rmd
# Chunk      : boxplot-isoform  (plot p_isoform_dist)
# Upstream deps inlined: libraries, load, matched-pairs, normalize,
#                        compute-frac-lost, + filtered DGEList for filt_isoforms
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

DGELIST_PATH      <- P("dge_isoform_longread_2026.3.3.rds")
DGELIST_FILT_PATH <- P("dge_isoform_longread_filtered_2026.3.3.rds")
CT_KEEP <- c("LAE", "AT2", "FB", "MV")

# ---- Load DGEList and parse metadata -------------------------------------
dge  <- readRDS(DGELIST_PATH)
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

# ---- Per-isoform fraction lost -------------------------------------------
results_list <- list()
for (ct_i in sort(unique(paired_donors$ct))) {
  donors_i <- paired_donors %>% filter(ct == ct_i) %>% pull(donor_id)
  for (d in donors_i) {
    smg1i_id <- samp %>% filter(ct == ct_i, donor_id == d, treatment == "Smg1i") %>% pull(bam_id)
    dmso_id  <- samp %>% filter(ct == ct_i, donor_id == d, treatment == "DMSO")  %>% pull(bam_id)
    stopifnot(length(smg1i_id) == 1, length(dmso_id) == 1)
    val_smg1i <- cpm_mat[, smg1i_id]
    val_dmso  <- cpm_mat[, dmso_id]
    delta     <- val_smg1i - val_dmso
    delta_pos <- pmax(delta, 0)
    frac_lost <- ifelse(val_smg1i > 0, delta_pos / val_smg1i, 0)
    results_list[[paste(ct_i, d, sep = "_")]] <- data.frame(
      txid = rownames(cpm_mat), ct = ct_i, donor_id = d,
      smg1i_val = val_smg1i, dmso_val = val_dmso, delta = delta,
      delta_pos = delta_pos, frac_lost = frac_lost, stringsAsFactors = FALSE)
  }
}
results_all <- bind_rows(results_list)
stopifnot(!anyDuplicated(iso_info$txid))
results_all <- results_all %>% left_join(iso_info, by = "txid")

# ---- Restrict to filtered isoforms (what limma-voom actually tested) -----
dge_filt      <- readRDS(DGELIST_FILT_PATH)
filt_isoforms <- rownames(dge_filt$counts)
plot_data     <- results_all %>% filter(frac_lost > 0, txid %in% filt_isoforms)

# ---- Figure --------------------------------------------------------------
p_isoform_dist <- ggplot(plot_data, aes(x = ct, y = frac_lost, fill = ct)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.3, outlier.alpha = 0.3) +
  labs(title    = "Per-Isoform Fraction Lost To NMD",
       subtitle = "(Smg1i − DMSO) / Smg1i CPM, positive deltas, filtered isoforms",
       x = "Cell Type", y = "Fraction Lost To NMD") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(color = "grey40"))

ggsave(file.path(HERE, "isoform_frac_lost_boxplot.png"),
       p_isoform_dist, width = 10, height = 6, dpi = 300, bg = "white")
