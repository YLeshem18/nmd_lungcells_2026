"""
SR/LR overview composite — 2x3 uniform multipanel layout.
Mirrors figures/multipanel/figure3_isopair_and_ptc/figure3_composite.py:
each panel PNG is rendered at CELL_W x CELL_H (6.0 x 4.0 in) by make_panels.R,
then laid out in a GridSpec(3, 2) of identical cells (zero aspect distortion),
panel letters at 22 pt bold, no frames.

  Row 1: A SR vs LR correlation   | B isoform landscape sharing
  Row 2: C DIE volcano            | D productive logFC SR vs LR
  Row 3: E NMD-susceptible sharing| F KEGG ER pathway (LAE)
"""

import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
# shared figure lib from the repo clone
LIB = Path("../repo_clone_test/nmd_isopair_analysis/figures/lib").resolve()
if LIB.exists():
    sys.path.insert(0, str(LIB))

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

LBL_FS = 22
LBL_COLOR = "#111111"
CELL_W, CELL_H = 6.0, 4.0
N_COLS, N_ROWS = 2, 3

PANELS = {
    "A": "panelA_correlation.png",
    "B": "panelB_landscape.png",
    "C": "panelC_volcano.png",
    "D": "panelD_logfc_box.png",
    "E": "panelE_sharing.png",
    "F": "panelF_pathway.png",
}
ROWS = [["A", "B"], ["C", "D"], ["E", "F"]]


def build_figure():
    paths = {k: HERE / v for k, v in PANELS.items()}
    fig = plt.figure(figsize=(CELL_W * N_COLS, CELL_H * N_ROWS))
    gs = GridSpec(N_ROWS, N_COLS, figure=fig,
                  hspace=0.10, wspace=0.02,
                  left=0.005, right=0.995, top=0.96, bottom=0.005)
    for ri, row in enumerate(ROWS):
        for ci, letter in enumerate(row):
            ax = fig.add_subplot(gs[ri, ci])
            ax.imshow(mpimg.imread(str(paths[letter])), aspect="auto")
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
            ax.text(-0.005, 1.02, letter, transform=ax.transAxes,
                    fontsize=LBL_FS, fontweight="bold", color=LBL_COLOR,
                    ha="left", va="bottom")
    return fig


def main():
    fig = build_figure()
    pdf = HERE / "figure_composite.pdf"
    png = HERE / "figure_composite.png"
    fig.savefig(pdf, facecolor="white")
    fig.savefig(png, dpi=300, facecolor="white")
    print(f"Saved: {pdf.name} and {png.name}")


if __name__ == "__main__":
    main()
