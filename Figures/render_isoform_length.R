# =============================================================================
# Supplemental figure: Isoform length distribution (known vs novel)
# Source Rmd : Isoform_Landscape.Rmd
# Chunk      : q3-overall-density  (plot p_len_overall)
# Upstream deps inlined: libraries, config, load-dgelist, load-sqanti
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(edgeR)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
count  <- dplyr::count;  slice  <- dplyr::slice;  first  <- dplyr::first

DATA <- "nmd_fig_data"
P    <- function(f) file.path(DATA, f)
HERE <- "supplement_figures"
dir.create(HERE, recursive = TRUE, showWarnings = FALSE)

DGE_FILT_PATH <- P("dge_isoform_longread_filtered_2026.3.3.rds")
SQ_CLASS_PATH <- P("nmd_lungcells_classification.txt")

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

# ---- SQANTI3 classification ----------------------------------------------
sq <- fread(SQ_CLASS_PATH)
sq_filt <- sq %>%
  filter(isoform %in% rownames(dge_filt$counts)) %>%
  mutate(
    known_novel = case_when(
      structural_category %in% c("full-splice_match", "incomplete-splice_match") ~ "Known",
      structural_category %in% c("novel_in_catalog",  "novel_not_in_catalog")    ~ "Novel",
      TRUE                                                                       ~ "Other"
    )
  )

# ---- Figure: overall length distribution ---------------------------------
sq_filt_len <- sq_filt %>% filter(!is.na(length))

p_len_overall <- ggplot(sq_filt_len, aes(x = log10(length), color = known_novel)) +
  geom_freqpoly(bins = 80, linewidth = 1) +
  scale_color_manual(values = c("Known" = "steelblue",
                                "Novel" = "coral",
                                "Other" = "gray60")) +
  labs(title    = "Isoform Length Distribution",
       subtitle = paste0("All ", comma(nrow(sq_filt_len)),
                         " isoforms — SQANTI3 transcript length"),
       x = "log10(Isoform Length, bp)", y = "Number Of Isoforms", color = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(color = "grey40"),
        legend.position = "top")

ggsave(file.path(HERE, "length_distribution.png"),
       p_len_overall, width = 12, height = 6, dpi = 300, bg = "white")
