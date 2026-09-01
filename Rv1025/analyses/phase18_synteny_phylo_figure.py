#!/usr/bin/env python3
"""phase18_synteny_phylo_figure.py — P3.2.b: the conserved operon along the Actinobacteria.

Companion / phylogenetic-context version of Figure 2: a schematic Actinobacteria
taxonomy (drawn as a cladogram bracket, NOT an inferred tree) with the four-gene
block eno-divIC-Rv1025-ppx2 rendered at each representative tip, the Rv1025 amino-acid
identity to H37Rv shown per tip. The point made visually: the ordered block is retained
from M. tuberculosis down to Bifidobacterium (identity 100% -> 51%), across orders whose
genomes are otherwise extensively rearranged. Data: résultats/synteny_operon.tsv.

Writes: article/supplementary_materials/figureS_synteny_phylo.{png,pdf}
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "article/supplementary_materials/figureS_synteny_phylo"

# Representative tips ordered by accepted Actinobacteria taxonomy (top = M. tuberculosis).
# (label, Rv1025 %id to H37Rv, order/clade). Identities from résultats/synteny_operon.tsv.
TIPS = [
    ("M. tuberculosis (Rv1025)", 100, "Mycobacteriaceae"),
    ("M. avium", 86, "Mycobacteriaceae"),
    ("M. abscessus", 83, "Mycobacteriaceae"),
    ("M. leprae (reductive)", 80, "Mycobacteriaceae"),
    ("Rhodococcus jostii", 76, "Nocardiaceae"),
    ("Nocardia farcinica", 73, "Nocardiaceae"),
    ("Tsukamurella paurometabola", 69, "Tsukamurellaceae"),
    ("Gordonia bronchialis", 68, "Gordoniaceae"),
    ("Corynebacterium glutamicum", 59, "Corynebacteriaceae"),
    ("Corynebacterium diphtheriae", 57, "Corynebacteriaceae"),
    ("Streptomyces coelicolor", 68, "Streptomycetales"),
    ("Cutibacterium acnes", 58, "Propionibacteriales"),
    ("Micrococcus luteus", 59, "Micrococcales"),
    ("Bifidobacterium longum", 51, "Bifidobacteriales"),
    ("Escherichia coli (out-group)", None, "out-of-phylum"),
    ("Bacillus subtilis (out-group)", None, "out-of-phylum"),
]
# clade brackets: (label, first_tip_idx, last_tip_idx, depth-level)
CLADES = [
    ("Corynebacteriales", 0, 9, 1),
    ("Actinobacteria (Actinomycetota)", 0, 13, 2),
]
GENES = ["eno", "divIC", "Rv1025", "ppx2"]
GCOL = {"eno": "#8da0cb", "divIC": "#fc8d62", "Rv1025": "#e41a1c", "ppx2": "#a6d854"}


def ident_color(pid):
    if pid is None:
        return "#dddddd"
    return matplotlib.colormaps["YlGnBu"](0.15 + 0.8 * (pid - 45) / 55.0)


fig, ax = plt.subplots(figsize=(11, 8.5))
n = len(TIPS)
y_of = {i: n - i for i in range(n)}          # top-down
label_x = 0.75                                 # species labels, left-aligned
block_x0 = 4.6                                 # operon block region
gene_w, gene_gap = 1.05, 0.18

for i, (label, pid, clade) in enumerate(TIPS):
    y = y_of[i]
    present = pid is not None
    # species label (left-aligned) + short connector to the block
    ax.text(label_x, y, label, ha="left", va="center", fontsize=9,
            style=("italic" if "out-group" not in label else "normal"))
    ax.plot([block_x0 - 0.55, block_x0 - 0.3], [y, y], color="#bbbbbb", lw=0.6, zorder=1)
    if not present:
        ax.text(block_x0 + 1.5, y, "operon absent", ha="left", va="center",
                fontsize=9, color="#999999", style="italic")
        continue
    # operon block: 4 arrows
    x = block_x0
    for g in GENES:
        ax.add_patch(FancyArrow(x, y, gene_w, 0, width=0.34, head_width=0.34,
                                head_length=0.26, length_includes_head=True,
                                color=GCOL[g], ec="black", lw=0.5, zorder=3))
        x += gene_w + gene_gap
    # Rv1025 identity chip at the right
    ax.add_patch(Rectangle((x + 0.2, y - 0.2), 0.5, 0.4,
                 color=ident_color(pid), ec="black", lw=0.4, zorder=3))
    ax.text(x + 0.9, y, f"{pid}%", ha="left", va="center", fontsize=8.5)

# cladogram brackets (schematic taxonomy), far left of the species labels
for label, a, b, lvl in CLADES:
    x = 0.15 + 0.28 * (2 - lvl)               # Actinobacteria (lvl2) outermost, Corynebacteriales (lvl1) inner
    ya, yb = y_of[b], y_of[a]
    ax.plot([x, x], [ya - 0.35, yb + 0.35], color="#333333", lw=1.4, zorder=1)
    ax.text(x - 0.10, (ya + yb) / 2, label, ha="right", va="center", rotation=90,
            fontsize=8.5, color="#333333")

# gene legend + identity note
handles = [Line2D([0], [0], marker=">", color="w", markerfacecolor=GCOL[g],
                  markeredgecolor="k", markersize=11, label=g) for g in GENES]
ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.30, -0.02),
          ncol=4, frameon=False, fontsize=9, handletextpad=0.3, columnspacing=1.0)
ax.text(block_x0, n + 1.2, "eno–divIC–Rv1025–ppx2 operon", fontsize=11, weight="bold")
ax.text(block_x0, n + 0.6, "Rv1025 a.a. identity to H37Rv shown per tip "
        "(gene order normalised left→right)", fontsize=8.5, color="#555555")

ax.set_xlim(-1.4, block_x0 + 4 * (gene_w + gene_gap) + 2.2)
ax.set_ylim(0.2, n + 1.8)
ax.axis("off")
fig.suptitle("The eno–divIC–Rv1025–ppx2 block is retained across the Actinobacteria phylum",
             y=0.98, fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(f"{OUT}.png", dpi=600, bbox_inches="tight")
fig.savefig(f"{OUT}.pdf", bbox_inches="tight")
print(f"Written: {OUT}.png / .pdf")
