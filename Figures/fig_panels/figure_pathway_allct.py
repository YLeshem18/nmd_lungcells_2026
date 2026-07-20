"""
All-cell-type KEGG ER protein-processing pathway (hsa04141), colored by SR gene logFC.
2x2 grid (AT2, LAE, FB, MV), each the full pathway diagram with its baked-in logFC color key.
"""
from pathlib import Path
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
CT = ["AT2", "LAE", "FB", "MV"]

def main():
    fig = plt.figure(figsize=(16, 11))
    gs = GridSpec(2, 2, figure=fig, hspace=0.08, wspace=0.02,
                  left=0.005, right=0.995, top=0.95, bottom=0.01)
    for i, ct in enumerate(CT):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        ax.imshow(mpimg.imread(str(HERE / f"hsa04141.{ct}.png")), aspect="auto")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
        ax.set_title(ct, fontsize=22, fontweight="bold", color="#111111", pad=4)
    fig.savefig(HERE / "figure_pathway_allct.pdf", facecolor="white")
    fig.savefig(HERE / "figure_pathway_allct.png", dpi=200, facecolor="white")
    print("Saved figure_pathway_allct.{pdf,png}")

if __name__ == "__main__":
    main()
