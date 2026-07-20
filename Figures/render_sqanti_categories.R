# =============================================================================
# Supplemental figure: SQANTI3 structural categories
# Source Rmd : Isoform_Landscape.Rmd
# Chunk      : q1-barplot  (plot p_cat)
# Upstream deps inlined: libraries, config, load-dgelist, load-sqanti, q1-total
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

# ---- SQANTI3 classification (matched to filtered isoforms) ---------------
sq <- fread(SQ_CLASS_PATH)
sq_filt <- sq %>%
  filter(isoform %in% rownames(dge_filt$counts)) %>%
  mutate(
    novelty = case_when(
      structural_category == "full-splice_match"       ~ "Known (FSM)",
      structural_category == "incomplete-splice_match" ~ "Known Variant (ISM)",
      structural_category == "novel_in_catalog"        ~ "Novel (NIC)",
      structural_category == "novel_not_in_catalog"    ~ "Novel (NNC)",
      TRUE                                             ~ structural_category
    ),
    known_novel = case_when(
      structural_category %in% c("full-splice_match", "incomplete-splice_match") ~ "Known",
      structural_category %in% c("novel_in_catalog",  "novel_not_in_catalog")    ~ "Novel",
      TRUE                                                                       ~ "Other"
    )
  )

# ---- q1-total: category table --------------------------------------------
cat_table <- sq_filt %>%
  count(structural_category, novelty, known_novel) %>%
  arrange(desc(n)) %>%
  mutate(pct = round(100 * n / sum(n), 2))

# ---- Figure --------------------------------------------------------------
p_cat <- ggplot(cat_table, aes(x = reorder(novelty, n), y = n, fill = known_novel)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0(comma(n), " (", pct, "%)")), hjust = -0.05, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("Known" = "steelblue",
                               "Novel" = "coral",
                               "Other" = "gray60")) +
  scale_y_continuous(labels = comma) +
  labs(title    = "SQANTI3 Structural Categories",
       subtitle = paste0("Total: ", comma(nrow(sq_filt)),
                         " isoforms from ",
                         comma(length(unique(sq_filt$associated_gene))),
                         " genes"),
       x = NULL, y = "Number Of Isoforms", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(color = "grey40"),
        legend.position = "top") +
  expand_limits(y = max(cat_table$n) * 1.25)

ggsave(file.path(HERE, "structural_categories.png"),
       p_cat, width = 10, height = 6, dpi = 300, bg = "white")
