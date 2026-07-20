#!/usr/bin/env Rscript
## =============================================================================
## Tan et al. (2025) reanalysis -- TRANSCRIPT (isoform) level
## =============================================================================
## Reproduces the manuscript claim (Section 2):
##   "the transcript-level overlap of UPF2-dependent NMD targets between hESCs
##    and NPCs increased to 41-50% (3,069 shared transcripts; 49.8% of hESC and
##    40.6% of NPC targets)."
##
## Method
##   1. Read the 8 transcript-level DET tables from Tan et al. (2025) Supplementary
##      Tables S1, S2, S4, S6 (external, published; see README).
##   2. Part A -- Tan-style thresholded overlap (their own cutoffs: PPDE > 0.95,
##      PostFC > 2 [stringent] or > 1.5 [looser]). Direct mirror of Tan Fig 6C/8D,
##      minus the pUPF1 + downstream-EJ filters that cannot be reconstructed from
##      the supplements alone.
##   3. Part B -- 8-condition transcript-level mashr:
##        bhat = log2(PostFC)
##        shat = |bhat| / qnorm(PPEE/2, lower.tail = FALSE)   (PPEE as p-value proxy)
##      Binary NMD-target call per condition: posterior mean > 0 AND lfsr < 0.05.
##      UPF2-dependent target set = UNION of the 3 UPF2 perturbations
##      (KD, iKO, iKO+KD). The headline overlap is ESC-union vs NPC-union.
##
## Determinism: set.seed(42) before the random subset used for the mashr null
## correlation estimate. With the supplied inputs this reproduces 3,069 / 49.8% /
## 40.6% exactly.
##
## The optional UPF3B 2-condition gene-level sanity check (Part C in the original
## working script) is NOT included here -- it depends on a separately produced
## gene-level mashr object and does not feed any manuscript number.
## =============================================================================

## ---- Paths (edit DATA_DIR to point at the Tan supplementary tables) ---------
## Tan supplementary .xlsx are NOT redistributed with this repo (see README).
## Download them from the source in README and set DATA_DIR accordingly, or set
## the TAN_DATA_DIR environment variable.
DATA_DIR <- Sys.getenv("TAN_DATA_DIR", unset = "data")
OUT_DIR  <- Sys.getenv("TAN_OUT_DIR",  unset = "output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---- Optional: user library on off-PATH Windows R --------------------------
local({
  user_lib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
  if (nzchar(user_lib) && dir.exists(user_lib) && !user_lib %in% .libPaths()) {
    .libPaths(c(user_lib, .libPaths()))
  }
})

suppressPackageStartupMessages({
  library(readxl); library(data.table); library(mashr); library(ashr)
})

cond_order <- c("ESC_UPF2_KD", "ESC_UPF2_iKO", "ESC_UPF2_iKO_KD", "ESC_UPF3B_KO",
                "NPC_UPF2_KD", "NPC_UPF2_iKO", "NPC_UPF2_iKO_KD", "NPC_UPF3B_KO")

## Which sheet of which Tan supplementary file holds each transcript-level DET.
tan_manifest <- data.frame(
  condition = cond_order,
  file = c("Supplementary_Table_S1 (1).xlsx", "Supplementary_Table_S1 (1).xlsx",
           "Supplementary_Table_S1 (1).xlsx", "Supplementary Table S2.xlsx",
           "Supplementary_Table_S4 (1).xlsx", "Supplementary_Table_S4 (1).xlsx",
           "Supplementary_Table_S4 (1).xlsx", "Supplementary_Table_S6 (1).xlsx"),
  sheet = c("Sheet4", "Sheet5", "Sheet6", "Sheet2",
            "Sheet4", "Sheet5", "Sheet6", "Sheet2"),
  stringsAsFactors = FALSE
)

read_tx_sheet <- function(file, sheet, cond) {
  ## Row 1 = title; row 2 = header (TXID, PPEE, PPDE, PostFC, RealFC, C1Mean, C2Mean)
  path <- file.path(DATA_DIR, file)
  if (!file.exists(path))
    stop(sprintf("Missing Tan supplementary file: %s\n  Set DATA_DIR / TAN_DATA_DIR (see README).", path))
  raw <- suppressWarnings(suppressMessages(
    read_excel(path, sheet = sheet, skip = 1, .name_repair = "minimal")))
  raw <- as.data.frame(raw)
  stopifnot(all(c("TXID", "PPEE", "PPDE", "PostFC") %in% colnames(raw)))
  raw <- raw[, c("TXID", "PPEE", "PPDE", "PostFC", "RealFC", "C1Mean", "C2Mean")]
  for (col in c("PPEE", "PPDE", "PostFC", "RealFC", "C1Mean", "C2Mean"))
    if (!is.numeric(raw[[col]])) raw[[col]] <- suppressWarnings(as.numeric(raw[[col]]))
  raw <- raw[!is.na(raw$TXID) & nchar(raw$TXID) > 0, , drop = FALSE]
  raw <- raw[is.finite(raw$PostFC) & is.finite(raw$PPEE) & is.finite(raw$PPDE), , drop = FALSE]
  raw$condition <- cond
  raw
}

cat("Reading 8 transcript-level DET tables from:", normalizePath(DATA_DIR, mustWork = FALSE), "\n")
tx_data <- lapply(seq_len(nrow(tan_manifest)), function(i)
  read_tx_sheet(tan_manifest$file[i], tan_manifest$sheet[i], tan_manifest$condition[i]))
names(tx_data) <- tan_manifest$condition
cat("Rows per condition (transcripts):\n"); print(sapply(tx_data, nrow))

## =============================================================================
## A. Thresholded transcript-level overlap (Tan's own cutoffs)
## =============================================================================
call_targets <- function(df, postfc_cut, ppde_cut = 0.95)
  df$TXID[df$PPDE > ppde_cut & df$PostFC > postfc_cut]

overlap_summary <- function(thr) {
  cat(sprintf("\n----- Thresholded: PPDE>0.95 & PostFC>%.1f -----\n", thr))
  sets <- lapply(tx_data, call_targets, postfc_cut = thr); names(sets) <- cond_order
  factors <- c("UPF2_KD", "UPF2_iKO", "UPF2_iKO_KD", "UPF3B_KO")
  df_per <- data.frame(factor = factors, esc = NA_integer_, npc = NA_integer_,
                       overlap = NA_integer_, pct_esc = NA_real_, pct_npc = NA_real_,
                       jaccard = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_along(factors)) {
    f <- factors[i]; A <- sets[[paste0("ESC_", f)]]; B <- sets[[paste0("NPC_", f)]]
    nI <- length(intersect(A, B)); nU <- length(union(A, B))
    df_per$esc[i] <- length(A); df_per$npc[i] <- length(B); df_per$overlap[i] <- nI
    df_per$pct_esc[i] <- round(100*nI/length(A), 1); df_per$pct_npc[i] <- round(100*nI/length(B), 1)
    df_per$jaccard[i] <- round(nI/nU, 3)
  }
  print(df_per); df_per
}
per_factor_stringent <- overlap_summary(2)
per_factor_looser    <- overlap_summary(1.5)
fwrite(per_factor_stringent, file.path(OUT_DIR, "tan_tx_overlap_per_factor_stringent.csv"))
fwrite(per_factor_looser,    file.path(OUT_DIR, "tan_tx_overlap_per_factor_looser.csv"))

## =============================================================================
## B. 8-condition transcript-level mashr
## =============================================================================
derive_tx <- function(df) {
  ## Drop PostFC == 1 (log2 = 0 -> SE undefined) and the PPEE >= 0.99 pure-noise
  ## tail (runaway SEs, no inferential signal). Tan's own DE call (PPDE > 0.95,
  ## i.e. PPEE < 0.05) is far stricter than this filter.
  df <- df[df$PostFC > 0 & df$PostFC != 1 & df$PPEE < 0.99, , drop = FALSE]
  df$bhat <- log2(df$PostFC)
  z_abs   <- qnorm(pmax(df$PPEE, .Machine$double.xmin) / 2, lower.tail = FALSE)
  df$se   <- abs(df$bhat) / z_abs
  df <- df[is.finite(df$bhat) & is.finite(df$se) & df$se > 0, , drop = FALSE]
  if (anyDuplicated(df$TXID)) {          # keep the most significant PPEE per TXID
    df <- df[order(df$PPEE), ]; df <- df[!duplicated(df$TXID), ]
  }
  df[, c("TXID", "bhat", "se")]
}

tx_dge    <- lapply(tx_data, derive_tx)
common_tx <- Reduce(intersect, lapply(tx_dge, function(d) d$TXID))
cat("\nTranscripts present in all 8 conditions:", length(common_tx), "\n")

bhat_tx <- matrix(NA_real_, length(common_tx), 8, dimnames = list(common_tx, cond_order))
shat_tx <- bhat_tx
for (cond in cond_order) {
  d <- tx_dge[[cond]]; idx <- match(common_tx, d$TXID)
  bhat_tx[, cond] <- d$bhat[idx]; shat_tx[, cond] <- d$se[idx]
}
ok <- apply(shat_tx, 1, function(x) all(is.finite(x) & x > 0)) &
      apply(bhat_tx, 1, function(x) all(is.finite(x)))
bhat_tx <- bhat_tx[ok, , drop = FALSE]; shat_tx <- shat_tx[ok, , drop = FALSE]
cat("Final tx matrices:", nrow(bhat_tx), "transcripts x", ncol(bhat_tx), "conditions\n")

cat("\n--- mashr fit (transcript level) ---\n")
mash_init <- mash_set_data(Bhat = bhat_tx, Shat = shat_tx)
m_1by1    <- mash_1by1(mash_init)
strong    <- get_significant_results(m_1by1, thresh = 0.05)

set.seed(42)                              # determinism for the null-correlation subset
ri      <- sample(seq_len(nrow(bhat_tx)), min(20000, nrow(bhat_tx)))
Vhat_tx <- estimate_null_correlation_simple(mash_set_data(bhat_tx[ri, ], shat_tx[ri, ]))
md_full <- mash_set_data(Bhat = bhat_tx, Shat = shat_tx, V = Vhat_tx)
U_pca   <- if (length(strong) >= 5) cov_pca(md_full, npc = 5, subset = strong) else list()
m_tx    <- mash(md_full, Ulist = c(U_pca, cov_canonical(md_full)), outputlevel = 2)

pm_tx   <- get_pm(m_tx);  colnames(pm_tx)   <- cond_order
lfsr_tx <- get_lfsr(m_tx); colnames(lfsr_tx) <- cond_order
is_target <- (pm_tx > 0) & (lfsr_tx < 0.05)     # binary NMD-target call
cat("\nNMD-target counts per condition (mashr, pm>0 & lfsr<0.05):\n"); print(colSums(is_target))

saveRDS(m_tx, file.path(OUT_DIR, "tan_tx_mashr_model.rds"))

## =============================================================================
## Headline: UPF2-union ESC vs NPC overlap (manuscript 3,069 / 49.8% / 40.6%)
## =============================================================================
upf2 <- c("UPF2_KD", "UPF2_iKO", "UPF2_iKO_KD")
esc_union <- rowSums(is_target[, paste0("ESC_", upf2), drop = FALSE]) >= 1
npc_union <- rowSums(is_target[, paste0("NPC_", upf2), drop = FALSE]) >= 1
nE <- sum(esc_union); nN <- sum(npc_union)
nI <- sum(esc_union & npc_union); nU <- sum(esc_union | npc_union)

cat("\n=== UPF2-dependent (union of 3 perturbations) ESC vs NPC ===\n")
cat(sprintf("  hESC UPF2 targets: %d\n  NPC  UPF2 targets: %d\n  Shared: %d\n", nE, nN, nI))
cat(sprintf("  %.1f%% of hESC also in NPC | %.1f%% of NPC also in hESC | Jaccard %.3f\n",
            100*nI/nE, 100*nI/nN, nI/nU))
cat("  (Manuscript: 3,069 shared; 49.8% of hESC, 40.6% of NPC)\n")

fwrite(data.frame(hESC_targets = nE, NPC_targets = nN, shared = nI,
                  pct_of_hESC = round(100*nI/nE, 1), pct_of_NPC = round(100*nI/nN, 1),
                  jaccard = round(nI/nU, 3)),
       file.path(OUT_DIR, "tan_tx_upf2_union_overlap.csv"))
cat("\nDone. Outputs in:", normalizePath(OUT_DIR, mustWork = FALSE), "\n")
