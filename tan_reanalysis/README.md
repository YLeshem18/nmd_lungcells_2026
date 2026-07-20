# Tan et al. (2025) reanalysis — transcript (isoform) level

Reanalysis of the Tan et al. (2025) UPF2/UPF3B perturbation data used in Section 2
of the manuscript. Joint adaptive-shrinkage (`mashr`) fit across the eight
perturbation conditions, then a Tan-style binary NMD-target overlap between hESC
and NPC.

Reproduces the manuscript claim:

> the transcript-level overlap of UPF2-dependent NMD targets between hESCs and NPCs
> increased to **41–50% (3,069 shared transcripts; 49.8% of hESC and 40.6% of NPC targets)**.

## Inputs (not redistributed here)

The four Tan et al. (2025) Supplementary Tables — **S1, S2, S4, S6** (`.xlsx`).
These are the authors' published data; download them from the paper's supplementary
material and place them in a directory of your choosing. The eight transcript-level
DET tables are read per the manifest in `tan_transcript_reanalysis.R`:

| Condition          | File                              | Sheet  |
|--------------------|-----------------------------------|--------|
| ESC_UPF2_KD        | Supplementary_Table_S1 (1).xlsx   | Sheet4 |
| ESC_UPF2_iKO       | Supplementary_Table_S1 (1).xlsx   | Sheet5 |
| ESC_UPF2_iKO_KD    | Supplementary_Table_S1 (1).xlsx   | Sheet6 |
| ESC_UPF3B_KO       | Supplementary Table S2.xlsx       | Sheet2 |
| NPC_UPF2_KD        | Supplementary_Table_S4 (1).xlsx   | Sheet4 |
| NPC_UPF2_iKO       | Supplementary_Table_S4 (1).xlsx   | Sheet5 |
| NPC_UPF2_iKO_KD    | Supplementary_Table_S4 (1).xlsx   | Sheet6 |
| NPC_UPF3B_KO       | Supplementary_Table_S6 (1).xlsx   | Sheet2 |

Each sheet has a one-row title above the header; columns used are
`TXID, PPEE, PPDE, PostFC`.

## Running

```bash
# point at the folder holding the Tan .xlsx files; output dir is created if absent
TAN_DATA_DIR=/path/to/tan_supplementary \
TAN_OUT_DIR=./output \
Rscript tan_transcript_reanalysis.R
```

Defaults: `DATA_DIR=data`, `OUT_DIR=output` (relative to the working directory).

## Method

1. **Thresholded overlap (Part A)** — Tan's own cutoffs (`PPDE > 0.95`,
   `PostFC > 2` stringent / `> 1.5` looser). Mirrors Tan Fig 6C/8D minus the
   pUPF1 and downstream-EJ filters, which cannot be reconstructed from the
   supplements alone.
2. **8-condition mashr (Part B)** — `bhat = log2(PostFC)`,
   `shat = |bhat| / qnorm(PPEE/2, lower.tail = FALSE)` (PPEE used as a p-value
   proxy). Filter `PostFC != 1` and `PPEE < 0.99`. `mash_1by1` → strong signals
   at `lfsr < 0.05` → `cov_pca(npc = 5)` + `cov_canonical` → `mash`.
   `set.seed(42)` before the 20,000-transcript null-correlation subset makes the
   fit deterministic.
3. **Binary NMD-target call** — per condition, `posterior mean > 0 AND lfsr < 0.05`.
   UPF2-dependent set = **union** of the three UPF2 perturbations (KD, iKO, iKO+KD).
   Headline overlap = hESC-union ∩ NPC-union.

## Outputs (`OUT_DIR/`)

- `tan_tx_upf2_union_overlap.csv` — the headline numbers (3,069 / 49.8% / 40.6%)
- `tan_tx_overlap_per_factor_stringent.csv`, `..._looser.csv` — per-factor thresholded overlaps
- `tan_tx_mashr_model.rds` — the fitted 8-condition mashr model

## Notes

- The optional UPF3B 2-condition gene-level sanity check from the original working
  script is omitted here; it depends on a separately produced gene-level mashr
  object and feeds no manuscript number.
- PPEE is a Bayesian posterior probability, not a frequentist p-value; using it as
  a p-value proxy to derive standard errors is a deliberate approximation (see the
  script header).
