# =============================================================================
# Supplemental figure: Isoform proportion vs expression (all + NMD-only), DMSO
# Source Rmd : interpret_isoform_patterns_mashr_2026.3.10.Rmd
# Chunks     : q3-all-iso (plot p_all)  and  q3-nmd-dmso (plot p_nmd)
# Upstream deps inlined: config, load-libs, load-dgelist, calculate-proportions,
#                        apply-filter, load-die, q3-identify, q3-baseline
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(scales)
  library(data.table)
  library(edgeR)
})
rename <- dplyr::rename; select <- dplyr::select; filter <- dplyr::filter
mutate <- dplyr::mutate; first  <- dplyr::first;  count  <- dplyr::count

DATA <- "nmd_fig_data"
P    <- function(f) file.path(DATA, f)
HERE <- "supplement_figures"
dir.create(HERE, recursive = TRUE, showWarnings = FALSE)

DATE_STAMP             <- "2026.3.10"
INPUT_DGELIST          <- P("dge_isoform_longread_2026.3.3.rds")
INPUT_DGELIST_FILTERED <- P("dge_isoform_longread_filtered_2026.3.3.rds")
INPUT_DGE_DIR          <- DATA
DGE_FILE_PATTERN       <- paste0("^nmd_mashr_die_.*_", gsub("\\.", "\\\\.", DATE_STAMP), "\\.csv$")

CPM_THRESHOLD          <- 0.25
SIGNIFICANCE_THRESHOLD <- 0.05

# Author map is at2/lae/fb/mv; the bundle ships die_{at2,dd,fb,mv}, so a "dd"
# alias is included defensively so LAE is not silently dropped.
CELLTYPE_MAP     <- c("at2" = "AT2", "lae" = "LAE", "dd" = "LAE", "fb" = "FB", "mv" = "MV")
CT_DISPLAY_ORDER <- c("AT2", "LAE", "FB", "MV")

# ---- Load DGELists + normalized CPM --------------------------------------
dge_isoform  <- readRDS(INPUT_DGELIST)
dge_filtered <- readRDS(INPUT_DGELIST_FILTERED)

fix_ct <- function(samples_df) {
  samples_df$ct <- sub("^.*_(DD_ALI|DD|AT|DO|FB|MV)$", "\\1", samples_df$sample_name)
  samples_df$ct <- ifelse(samples_df$ct == "DO", "DO_ALI", samples_df$ct)
  samples_df$ct[samples_df$ct == "DD"] <- "LAE"
  samples_df$ct[samples_df$ct == "AT"] <- "AT2"
  samples_df
}
dge_isoform$samples  <- fix_ct(dge_isoform$samples)
dge_filtered$samples <- fix_ct(dge_filtered$samples)

cpm_normalized <- cpm(dge_filtered, normalized = TRUE, log = FALSE)

# ---- Isoform proportions -------------------------------------------------
sample_meta            <- dge_filtered$samples
sample_meta$cpm_col_id <- rownames(sample_meta)

cpm_long <- as.data.frame(cpm_normalized) %>%
  mutate(transcript_id = dge_filtered$genes$txid) %>%
  pivot_longer(cols = -transcript_id, names_to = "cpm_col_id", values_to = "cpm") %>%
  left_join(sample_meta %>% select(cpm_col_id, sample_name, treatment, ct),
            by = "cpm_col_id") %>%
  rename(cell_type = ct)

cpm_long <- cpm_long %>%
  left_join(dge_filtered$genes %>% select(txid, gene_id),
            by = c("transcript_id" = "txid")) %>%
  rename(gene_group_id = gene_id) %>%
  filter(cell_type %in% CT_DISPLAY_ORDER)

all_proportions <- cpm_long %>%
  group_by(cell_type, treatment, cpm_col_id, gene_group_id) %>%
  mutate(
    gene_total_cpm     = sum(cpm),
    isoform_proportion = cpm / gene_total_cpm,
    isoform_proportion = ifelse(is.nan(isoform_proportion), 0, isoform_proportion),
    n_isoforms         = n_distinct(transcript_id)
  ) %>%
  ungroup() %>%
  rename(normalized_cpm = cpm)

# ---- Abundance filter ----------------------------------------------------
filtered_dict <- all_proportions %>%
  group_by(cell_type, transcript_id) %>%
  summarize(max_cpm = max(normalized_cpm), .groups = "drop") %>%
  filter(max_cpm >= CPM_THRESHOLD)

all_proportions <- all_proportions %>%
  inner_join(filtered_dict %>% select(cell_type, transcript_id),
             by = c("cell_type", "transcript_id"))

# ---- Load DIE (mashr) results --------------------------------------------
dge_files <- list.files(INPUT_DGE_DIR, pattern = DGE_FILE_PATTERN, full.names = TRUE)
all_dge <- as.data.frame(data.table::rbindlist(lapply(dge_files, function(.x){
  df    <- as.data.frame(data.table::fread(.x))
  parts <- strsplit(basename(.x), "_")[[1]]
  df$cell_type_file <- parts[4]
  df
}), fill = TRUE))
all_dge$cell_type <- CELLTYPE_MAP[all_dge$cell_type_file]
all_dge  <- all_dge %>% filter(!is.na(cell_type))
stopifnot(nrow(all_dge) > 0)
if ("txid" %in% colnames(all_dge)) {
  colnames(all_dge)[colnames(all_dge) == "txid"] <- "transcript_id"
}

# ---- NMD susceptible transcripts -----------------------------------------
nmd_responsive <- all_dge %>%
  filter(logFC > 0, adj.P.Val < SIGNIFICANCE_THRESHOLD) %>%
  select(cell_type, transcript_id, logFC, adj.P.Val)

baseline_proportions <- all_proportions %>%
  filter(treatment == "DMSO") %>%
  inner_join(nmd_responsive, by = c("cell_type", "transcript_id")) %>%
  select(cell_type, cpm_col_id, gene_group_id, transcript_id,
         normalized_cpm, isoform_proportion, n_isoforms, logFC, adj.P.Val)

# ---- Figure A: All isoforms — proportion vs expression in DMSO -----------
all_dmso_proportions <- all_proportions %>%
  filter(treatment == "DMSO") %>%
  group_by(transcript_id, gene_group_id) %>%
  summarize(mean_cpm        = mean(normalized_cpm),
            mean_proportion = mean(isoform_proportion),
            .groups = "drop") %>%
  filter(mean_cpm > 0)

p_all <- ggplot(all_dmso_proportions, aes(x = mean_cpm, y = mean_proportion)) +
  geom_point(alpha = 0.3, size = 0.5, color = "gray40") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "orange") +
  scale_x_log10() +
  labs(title    = "All Isoforms — Proportion Vs Expression In DMSO",
       subtitle = "Lines: 50% (red), 10% (orange)",
       x = "Mean Normalized CPM In DMSO (log10)",
       y = "Mean Isoform Proportion In DMSO") +
  theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40"))

ggsave(file.path(HERE, paste0("q3_all_iso_dmso_", DATE_STAMP, ".png")),
       p_all, width = 8, height = 5, dpi = 300, bg = "white")

# ---- Figure B: NMD susceptible isoforms — proportion vs expression in DMSO
nmd_dmso_proportions <- baseline_proportions %>%
  group_by(transcript_id, gene_group_id) %>%
  summarize(mean_cpm        = mean(normalized_cpm),
            mean_proportion = mean(isoform_proportion),
            .groups = "drop") %>%
  filter(mean_cpm > 0)

p_nmd <- ggplot(nmd_dmso_proportions, aes(x = mean_cpm, y = mean_proportion)) +
  geom_point(alpha = 0.4, size = 0.8, color = "gray70") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "orange") +
  scale_x_log10() +
  labs(title    = "NMD Susceptible Isoforms — Proportion Vs Expression In DMSO",
       subtitle = "Lines: 50% (red), 10% (orange)",
       x = "Mean Normalized CPM In DMSO (log10)",
       y = "Mean Isoform Proportion In DMSO") +
  theme_minimal(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40"))

ggsave(file.path(HERE, paste0("q3_nmd_dmso_", DATE_STAMP, ".png")),
       p_nmd, width = 8, height = 5, dpi = 300, bg = "white")
