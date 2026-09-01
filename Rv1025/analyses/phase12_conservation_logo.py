#!/usr/bin/env python3
"""
P2.7 / P6.3 - Conservation du site Cys-His-Glu sur l'alignement Pfam CURÉ (PF04417 seed, 39 séq),
indépendant du MSA AF3 redondant. Piste de conservation par colonne (% identité sur non-gap) le long de
l'alignement, avec les 3 résidus du site (Glu59, Cys113, His115) mis en évidence : ils sont parmi les
positions invariantes de la famille curée. Sortie : article/figures/figureS_conservation.{png,pdf}.
"""
import os
from collections import Counter
from Bio import AlignIO
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/christophe/docs/codes/mtbc/Rv1025"
STO = f"{ROOT}/résultats/pfam/PF04417_seed.sto"
OUT = f"{ROOT}/article/figures/figureS_conservation"
# colonnes seed du site (repérées : Glu conservé col68 ; motif C-x-H col184/186)
# (label, aa, dx label, y label, ha)  -- dx/ha pour éviter le chevauchement Cys/His adjacents
SITE = {68: ("Glu59", "E", 0, 120, "center"),
        184: ("Cys113", "C", -9, 132, "right"),
        186: ("His115", "H", 9, 120, "left")}

def main():
    aln = AlignIO.read(STO, "stockholm")
    n = len(aln); L = aln.get_alignment_length()
    cons, cover = [], []
    for j in range(L):
        col = [aln[i, j] for i in range(n)]
        ng = [c for c in col if c not in "-."]
        cover.append(len(ng) / n)
        if not ng:
            cons.append(0); continue
        res, cnt = Counter(ng).most_common(1)[0]
        cons.append(100 * cnt / len(ng))
    # ne tracer que les colonnes bien couvertes (>=50% des séq) pour la lisibilité
    xs = [j for j in range(L) if cover[j] >= 0.5]
    ys = [cons[j] for j in xs]
    idx = {j: k for k, j in enumerate(xs)}
    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.bar(range(len(xs)), ys, color="#bbb", width=1.0)
    for j, (lab, aa, dx, ytxt, ha) in SITE.items():
        if j in idx:
            k = idx[j]
            ax.bar(k, cons[j], color="#C44E52", width=1.2, zorder=3)
            ax.annotate(f"{lab} ({aa}, {cons[j]:.0f}%)", (k, cons[j]), (k + dx, ytxt),
                        ha=ha, fontsize=9, color="#C44E52", fontweight="bold",
                        arrowprops=dict(arrowstyle="-", color="#C44E52", lw=1))
    ax.set_ylim(0, 145); ax.set_xlim(-1, len(xs))
    ax.set_ylabel("column conservation (% of non-gap)")
    ax.set_xlabel("alignment position (well-covered columns of the PF04417 seed, 39 sequences)")
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("The Cys113–His115–Glu59 metal-site residues are invariant in the curated DUF501 family", fontsize=10.5)
    plt.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT + ".pdf", bbox_inches="tight")
    print("écrit:", OUT + ".png/.pdf ; site:", {SITE[j][0]: round(cons[j]) for j in SITE if j < L})

if __name__ == "__main__":
    main()
