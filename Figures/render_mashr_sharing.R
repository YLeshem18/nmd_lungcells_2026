# =============================================================================
# Supplemental figure: mashr pairwise sharing (LR isoform)
# Source Rmd : Isoform_Level_Quantification.Rmd  (section 8.7 Effect Sharing)
# Chunks     : mashr-sharing / mashr-sharing-plot  (plot p_share)
#
# NOTE ON SELF-CONTAINMENT:
#   The plot needs the fitted mashr model object `m` (get_pairwise_sharing(m)).
#   The model RDS is NOT in the data bundle, so this script REBUILDS the model
#   from the filtered + TMM-normalized DGEList (dge_isoform_longread_filtered),
#   replicating Steps 6-8 of the source Rmd: voom -> duplicateCorrelation ->
#   lmFit -> per-cell-type contrasts -> eBayes -> mashr fit -> pairwise sharing.
#   This is compute-heavy (full mashr fit) but reproduces the author's logic.
# =============================================================================
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(reshape2)
})

DATA <- "nmd_fig_data"
P    <- function(f) file.path(DATA, f)
HERE <- "supplement_figures"
dir.create(HERE, recursive = TRUE, showWarnings = FALSE)

DGE_FILT_PATH <- P("dge_isoform_longread_filtered_2026.3.3.rds")

# ---- Load filtered + normalized DGEList (already filterByExpr + TMM) ------
dge <- readRDS(DGE_FILT_PATH)
# Cell-type factor was saved with LAE as reference and treatment DMSO as reference.

# ---- Step 5/6 — voom + model fit (design as authored) --------------------
design_all <- model.matrix(~ ct + treatment + ct:treatment, data = dge$samples)
v <- voom(dge, design = design_all)

corfit <- duplicateCorrelation(v, design_all, block = dge$samples$id)
fit <- lmFit(v, design_all,
             block       = dge$samples$id,
             correlation = corfit$consensus.correlation)
fit <- eBayes(fit)

# ---- Step 7 — per-cell-type Smg1i contrasts ------------------------------
colnames(fit$design)       <- make.names(colnames(fit$design))
colnames(fit$coefficients) <- colnames(fit$design)

contrast_matrix <- makeContrasts(
  Smg1i_in_LAE = treatmentSmg1i,
  Smg1i_in_AT2 = treatmentSmg1i + ctAT2.treatmentSmg1i,
  Smg1i_in_FB  = treatmentSmg1i + ctFB.treatmentSmg1i,
  Smg1i_in_MV  = treatmentSmg1i + ctMV.treatmentSmg1i,
  levels = fit$design
)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# ---- Step 8 — mashr ------------------------------------------------------
if (!require("ashr",     quietly = TRUE)) install.packages("ashr")
if (!require("mashr",    quietly = TRUE)) BiocManager::install("mashr")
if (!require("reshape2", quietly = TRUE)) install.packages("reshape2")
library(mashr); library(ashr)

bhat <- fit2$coefficients
shat <- fit2$coefficients / fit2$t
finite_rows <- apply(shat, 1, function(x) all(is.finite(x) & x > 0))
bhat <- bhat[finite_rows, ]
shat <- shat[finite_rows, ]

mash_data_init <- mash_set_data(Bhat = bhat, Shat = shat)
m_1by1         <- mash_1by1(mash_data_init)
strong         <- get_significant_results(m_1by1, thresh = 0.05)

set.seed(42)
n_random   <- min(10000, nrow(bhat))
random_idx <- sample(seq_len(nrow(bhat)), n_random)
mash_data_random <- mash_set_data(Bhat = bhat[random_idx, ], Shat = shat[random_idx, ])
Vhat <- estimate_null_correlation_simple(mash_data_random)

mash_data <- mash_set_data(Bhat = bhat, Shat = shat, V = Vhat)

if (length(strong) >= 5) {
  U_pca <- cov_pca(mash_data, npc = min(5, ncol(bhat) - 1), subset = strong)
} else {
  U_pca <- list()
}
U_c   <- cov_canonical(mash_data)
U_all <- c(U_pca, U_c)

m <- mash(mash_data, Ulist = U_all, outputlevel = 2)

# ---- 8.7 Effect sharing --------------------------------------------------
n_sig_total <- length(get_significant_results(m, thresh = 0.05))
stopifnot("Too few significant isoforms for sharing analysis" = n_sig_total >= 2)
sharing <- get_pairwise_sharing(m, factor = 0.5, lfsr_thresh = 0.05)

# ---- 8.8 Sharing heatmap (the figure) ------------------------------------
CT_ORDER <- c("LAE", "AT2", "FB", "MV")
sharing_clean <- sharing
rownames(sharing_clean) <- gsub("Smg1i_in_", "", rownames(sharing_clean))
colnames(sharing_clean) <- gsub("Smg1i_in_", "", colnames(sharing_clean))

sharing_df <- melt(sharing_clean)
colnames(sharing_df) <- c("CellType1", "CellType2", "Sharing")
sharing_df$CellType1 <- factor(sharing_df$CellType1, levels = CT_ORDER)
sharing_df$CellType2 <- factor(sharing_df$CellType2, levels = CT_ORDER)

p_share <- ggplot(sharing_df, aes(x = CellType1, y = CellType2, fill = Sharing)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(Sharing, 2)),
            size = 4, fontface = "bold",
            color = ifelse(sharing_df$Sharing > 0.6, "white", "grey20")) +
  scale_fill_gradient(low = "white", high = "#185FA5",
                      limits = c(0, 1), name = "Sharing") +
  labs(title    = "Pairwise Effect Sharing",
       subtitle = "Same sign & within 2-fold magnitude (lfsr < 0.05)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(color = "grey40", hjust = 0.5))

ggsave(file.path(HERE, "sharing_heatmap_lr_iso.png"),
       p_share, width = 8, height = 7, dpi = 300, bg = "white")
