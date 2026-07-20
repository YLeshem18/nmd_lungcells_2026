"""
Productive figure composite.
Row 1 (two 6x4 cells):  A % output lost (isoform)   | B NMD vs non-NMD productive
Rows 2-4 (full-width 12x3.4): C SHMT2 | D SRSF2 | E PCNA  (structure + donor-paired CI)
imshow of pre-rendered PNGs (each at its cell aspect) -> no distortion; frameless; 22pt tags.
"""
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path

HERE = Path(__file__).resolve().parent
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial","Helvetica Neue","Helvetica","DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
LBL_FS = 22

TOP = {"A": "panelA_outputlost.png", "B": "panelB_nmd_vs_nonmd.png"}
ROWS = {"C": "gene_row_SHMT2.png", "D": "gene_row_SRSF2.png", "E": "gene_row_PCNA.png"}

def place(ax, png, letter):
    ax.imshow(mpimg.imread(str(HERE / png)), aspect="auto")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    ax.text(-0.005, 1.02, letter, transform=ax.transAxes,
            fontsize=LBL_FS, fontweight="bold", color="#111111", ha="left", va="bottom")

def main():
    fig = plt.figure(figsize=(12, 14.7))
    gs = GridSpec(4, 2, figure=fig, height_ratios=[4.0, 3.4, 3.4, 3.4],
                  hspace=0.14, wspace=0.03, left=0.01, right=0.99, top=0.945, bottom=0.01)
    place(fig.add_subplot(gs[0, 0]), TOP["A"], "A")
    place(fig.add_subplot(gs[0, 1]), TOP["B"], "B")
    place(fig.add_subplot(gs[1, :]), ROWS["C"], "C")
    place(fig.add_subplot(gs[2, :]), ROWS["D"], "D")
    place(fig.add_subplot(gs[3, :]), ROWS["E"], "E")
    fig.savefig(HERE / "figure_productive_composite.pdf", facecolor="white")
    fig.savefig(HERE / "figure_productive_composite.png", dpi=300, facecolor="white")
    print("Saved figure_productive_composite.{pdf,png}")

if __name__ == "__main__":
    main()
