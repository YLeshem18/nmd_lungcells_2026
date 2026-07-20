# =============================================================================
# Supplemental figure: Pairwise expression similarity across cell types
# Source Rmd : Isoform_Landscape.Rmd
# Chunk      : corr-jaccard-heatmaps  (plot p_cor, Spearman heatmap)
# Upstream deps inlined: libraries, config, load-dgelist,
#                        build-cpm-matrices (cpm_ct_mat), corr-all (cor_spearman)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(edgeR)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
count  <- dplyr::count;  slice  <- dplyr::slice;  first  <- dplyr::first

DATA <- "nmd_fig_data"
P    <- function(f) file.path(DATA, f)
HERE <- "supplement_figures"
dir.create(HERE, recursive = TRUE, showWarnings = FALSE)

DGE_FILT_PATH <- P("dge_isoform_longread_filtered_2026.3.3.rds")
CT_KEEP <- c("LAE", "AT2", "FB", "MV")
relabel_ct <- function(ct_vec) {
  ct_vec <- as.character(ct_vec)
  ct_vec[ct_vec == "AT"] <- "AT2"
  ct_vec[ct_vec == "DD"] <- "LAE"
  ct_vec
}

# ---- Load filtered DGEList; keep 4 primary CTs, DMSO only ----------------
dge_filt <- readRDS(DGE_FILT_PATH)
dge_filt$samples$ct <- relabel_ct(dge_filt$samples$ct)
ct_char   <- as.character(dge_filt$samples$ct)
trt_char  <- as.character(dge_filt$samples$treatment)
keep_samp <- ct_char %in% CT_KEEP & trt_char == "DMSO"
dge_filt  <- dge_filt[, keep_samp, keep.lib.sizes = FALSE]

samp_filt <- dge_filt$samples
cts       <- sort(unique(as.character(samp_filt$ct)))

# ---- Mean log-CPM per cell type (all samples in the restricted set) ------
build_cpm_matrix <- function(dge, samp, cts, treatment = NULL) {
  out <- list()
  for (ct_i in cts) {
    idx <- if (is.null(treatment)) {
      which(as.character(samp$ct) == ct_i)
    } else {
      which(as.character(samp$ct) == ct_i & as.character(samp$treatment) == treatment)
    }
    if (length(idx) == 0) next
    ct_cpm <- cpm(dge[, idx], log = TRUE, prior.count = 1)
    out[[ct_i]] <- rowMeans(ct_cpm)
  }
  do.call(cbind, out)
}
cpm_ct_mat <- build_cpm_matrix(dge_filt, samp_filt, cts)

# ---- Spearman correlation + heatmap (the figure) -------------------------
cor_spearman <- cor(cpm_ct_mat, method = "spearman")

cor_df <- melt(cor_spearman)
colnames(cor_df) <- c("CT1", "CT2", "Correlation")

p_cor <- ggplot(cor_df, aes(x = CT1, y = CT2, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(Correlation, 2)),
            size = 4, fontface = "bold",
            color = ifelse(cor_df$Correlation > 0.85, "white", "grey20")) +
  scale_fill_gradient2(low = "#D73027", mid = "white", high = "#2166AC",
                       midpoint = min(cor_spearman),
                       limits   = c(min(cor_spearman), 1),
                       name = "Spearman r") +
  labs(title    = "Pairwise Expression Similarity Across Cell Types",
       subtitle = "Spearman correlation of mean log-CPM per isoform",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
        axis.text.y   = element_text(face = "bold", size = 14),
        panel.grid    = element_blank(),
        plot.title    = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40", size = 13))

ggsave(file.path(HERE, "pairwise_expression_similarity_spearman.png"),
       p_cor, width = 10, height = 8, dpi = 300, bg = "white")
