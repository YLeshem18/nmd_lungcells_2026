# =============================================================================
# Supplemental figure: SR<->LR gene-expression correlation, by cell type
# Source Rmd : correlation_analysis.Rmd
# Chunk      : corr-by-celltype-panel  (plot p_cor_panel)
# Upstream deps inlined: load-libs, config, short-read, long-read,
#                        combine-short-long, corr-by-celltype
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
count  <- dplyr::count;  slice  <- dplyr::slice;  first  <- dplyr::first

# ---- Relative paths -------------------------------------------------------
DATA <- "nmd_fig_data"
P    <- function(f) file.path(DATA, f)
HERE <- "supplement_figures"
dir.create(HERE, recursive = TRUE, showWarnings = FALSE)

DGELIST_GENE  <- P("dge_gene_unfiltered_2026.1.2.rds")
COUNTS_PATH   <- P("salmon.merged.gene_counts_length_scaled.FULL.rds")   # (FULL variant also in bundle)
ISO_COUNT_MAT <- P("nmd_lungcells_filtered.count_matrix.txt")
DGELIST_ISO   <- P("dge_isoform_longread_2026.3.3.rds")

CT_KEEP    <- c("AT", "DD", "FB", "MV")
CT_DISPLAY <- c("AT" = "AT2", "DD" = "LAE", "FB" = "FB", "MV" = "MV")

extract_ct <- function(sample_names) {
  sapply(strsplit(sample_names, "_"), function(x) {
    paste(x[3:length(x)], collapse = "_")
  })
}

# ---- Short-Read -----------------------------------------------------------
dge_gene  <- readRDS(DGELIST_GENE)
se_counts <- readRDS(COUNTS_PATH)
counts_matrix <- assay(se_counts, "counts")

parse_colname <- function(colname) {
  parts      <- strsplit(colname, "_")[[1]]
  first_part <- parts[1]
  if (startsWith(first_part, "AT2")) {
    ct        <- "AT"
    sample_id <- substring(first_part, 4)
  } else {
    ct        <- substring(first_part, 1, 2)
    sample_id <- substring(first_part, 3)
  }
  treatment <- parts[2]
  paste(sample_id, treatment, ct, sep = "_")
}
colnames(counts_matrix) <- sapply(colnames(counts_matrix), parse_colname)

modified_ct   <- gsub("DO_ALI", "DO", dge_gene$samples$ct)
modified_ct   <- gsub("AT2",    "AT", modified_ct)
canonical_ids <- paste(dge_gene$samples$sample_id,
                       dge_gene$samples$treatment,
                       modified_ct, sep = "_")
dd_ali_in_sr <- sapply(colnames(counts_matrix), function(nm) {
  if (!grepl("_DD$", nm)) return(FALSE)
  potential_ali <- gsub("_DD$", "_DD_ALI", nm)
  potential_ali %in% canonical_ids
})
if (sum(dd_ali_in_sr) > 0) counts_matrix <- counts_matrix[, !dd_ali_in_sr]

ct_per_sample <- extract_ct(colnames(counts_matrix))
keep_idx      <- ct_per_sample %in% CT_KEEP
if (sum(!keep_idx) > 0) counts_matrix <- counts_matrix[, keep_idx]

rownames(counts_matrix) <- sub("\\..*", "", rownames(counts_matrix))
counts_combined <- counts_matrix

# ---- Long-Read ------------------------------------------------------------
counts_iso <- fread(ISO_COUNT_MAT, data.table = FALSE)
dge_iso    <- readRDS(DGELIST_ISO)
tx2gene    <- unique(dge_iso$genes[!is.na(dge_iso$genes$txid), c("txid", "gene_id")])

rownames(counts_iso) <- counts_iso$id
counts_iso$id        <- NULL
count_mat            <- round(as.matrix(counts_iso))

nms     <- gsub("^Sample\\d+_", "", colnames(count_mat))
new_nms <- sapply(nms, function(x) {
  parts     <- strsplit(x, "_")[[1]]
  treatment <- parts[length(parts)]
  sample_id <- parts[length(parts) - 1]
  ct        <- paste(parts[1:(length(parts) - 2)], collapse = "_")
  paste(sample_id, treatment, ct, sep = "_")
})
colnames(count_mat) <- new_nms

ct_per_sample_lr <- extract_ct(colnames(count_mat))
keep_idx_lr      <- ct_per_sample_lr %in% CT_KEEP
if (sum(!keep_idx_lr) > 0) count_mat <- count_mat[, keep_idx_lr]

idx              <- match(rownames(count_mat), tx2gene$txid)
gene_ids_for_iso <- tx2gene$gene_id[idx]
has_gene         <- !is.na(gene_ids_for_iso)
counts_with_gene <- count_mat[has_gene, ]
gene_ids_for_iso <- sub("\\..*", "", gene_ids_for_iso[has_gene])
counts_genelevel <- rowsum(counts_with_gene, group = gene_ids_for_iso)

# ---- Combine short + long -------------------------------------------------
common_genes   <- intersect(rownames(counts_combined), rownames(counts_genelevel))
common_samples <- intersect(colnames(counts_combined), colnames(counts_genelevel))
stopifnot(length(common_genes) > 0, length(common_samples) > 0)

counts_short_matched <- counts_combined[common_genes, common_samples]
counts_long_matched  <- counts_genelevel[common_genes, common_samples]

expressed_short <- rowSums(counts_short_matched) > 0
expressed_long  <- rowSums(counts_long_matched)  > 0
keep_genes      <- expressed_short | expressed_long
counts_short_matched <- counts_short_matched[keep_genes, ]
counts_long_matched  <- counts_long_matched[keep_genes, ]

dge_short_matched <- calcNormFactors(DGEList(counts = counts_short_matched), method = "TMM")
dge_long_matched  <- calcNormFactors(DGEList(counts = counts_long_matched),  method = "TMM")
cpm_short_matched <- cpm(dge_short_matched, normalized.lib.sizes = TRUE, log = TRUE)
cpm_long_matched  <- cpm(dge_long_matched,  normalized.lib.sizes = TRUE, log = TRUE)

# ---- Per-cell-type correlation -------------------------------------------
sample_names <- colnames(cpm_short_matched)
cell_types   <- extract_ct(sample_names)
unique_cts   <- intersect(CT_KEEP, unique(cell_types))

ct_results <- list()
for (ct in unique_cts) {
  ct_samples   <- sample_names[cell_types == ct]
  cpm_short_ct <- cpm_short_matched[, ct_samples]
  cpm_long_ct  <- cpm_long_matched[, ct_samples]
  spear <- sapply(seq_len(ncol(cpm_short_ct)), function(i)
    cor(cpm_short_ct[, i], cpm_long_ct[, i], method = "spearman"))
  pears <- sapply(seq_len(ncol(cpm_short_ct)), function(i)
    cor(cpm_short_ct[, i], cpm_long_ct[, i], method = "pearson"))
  ct_results[[ct]] <- list(samples = ct_samples,
                           spearman = mean(spear), pearson = mean(pears))
}

# ---- Multi-panel hex scatter (the figure) --------------------------------
ct_plots <- list()
for (ct in unique_cts) {
  ct_samples   <- ct_results[[ct]]$samples
  cpm_short_ct <- cpm_short_matched[, ct_samples]
  cpm_long_ct  <- cpm_long_matched[, ct_samples]
  plot_data <- data.frame(short_read = as.vector(cpm_short_ct),
                          long_read  = as.vector(cpm_long_ct))
  ct_plots[[ct]] <- ggplot(plot_data, aes(x = short_read, y = long_read)) +
    geom_hex(bins = 50) +
    scale_fill_gradientn(colors = c("white", "lightblue", "blue", "darkblue", "red"),
                         trans = "log10", name = "Count") +
    annotate("text", x = min(plot_data$short_read), y = max(plot_data$long_read),
             hjust = 0, vjust = 1,
             label = paste0("r = ",       sprintf("%.3f", ct_results[[ct]]$pearson),
                            "\nρ = ", sprintf("%.3f", ct_results[[ct]]$spearman)),
             size = 4, fontface = "bold") +
    labs(title = CT_DISPLAY[ct], x = "Short-Read log2(CPM)", y = "Long-Read log2(CPM)") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(size = 14, face = "bold"),
          legend.position  = "none")
}

p_cor_panel <- wrap_plots(ct_plots, ncol = 2) +
  plot_annotation(
    title    = "Short-Read Vs Long-Read Gene Expression By Cell Type",
    subtitle = "log2(CPM) with TMM normalization; Pearson r and Spearman ρ per cell type",
    theme    = theme(plot.title    = element_text(size = 16, face = "bold"),
                     plot.subtitle = element_text(size = 12, color = "grey40")))

ggsave(file.path(HERE, "sr_vs_lr_expression_by_celltype_panel.png"),
       p_cor_panel, width = 12, height = 10, dpi = 300, bg = "white")
