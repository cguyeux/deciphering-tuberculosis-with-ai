#!/usr/bin/env python3
"""
P6.3 - Figure de synténie : opéron eno-divIC-Rv1025-ppx2 conservé à travers les Actinobactéries.
Lit résultats/synteny_operon.tsv ; dessine un ladder (1 génome/ligne) d'arrows de gènes, orientation
normalisée (opéron lu eno→divIC→Rv1025→ppx2 ; le brin réel est noté en légende), Rv1025 mis en évidence,
% identité Rv1025 à droite ; témoins hors-phylum (E. coli, B. subtilis) = DUF501 absent.
Sortie : article/figures/figure2_synteny.{png,pdf} (600 dpi).
"""
import csv, re, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Patch

ROOT = "/home/christophe/docs/codes/mtbc/Rv1025"
TSV = f"{ROOT}/résultats/synteny_operon.tsv"
OUT = f"{ROOT}/article/figures/figure2_synteny"

# (clé TSV, label affiché, clade)
ROWS = [
    ("M_bovis(NC_002945)", "M. bovis", "MTBC"),
    ("M_canettii(CP007299)", "M. canettii", "MTBC (outgroup)"),
    ("M_leprae", "M. leprae", "Mycobacterium"),
    ("M_avium", "M. avium", "Mycobacterium"),
    ("M_kansasii", "M. kansasii", "Mycobacterium"),
    ("M_xenopi", "M. xenopi", "Mycobacterium"),
    ("Mabscessus_ATCC19977", "M. abscessus", "Mycobacteriaceae"),
    ("Cdiphtheriae_NCTC13129", "C. diphtheriae", "Corynebacteriales"),
    ("Nfarcinica_IFM10152", "N. farcinica", "Corynebacteriales"),
    ("Rjostii_RHA1", "R. jostii", "Corynebacteriales"),
    ("Scoelicolor_A3-2", "S. coelicolor", "Streptomycetales"),
    ("Blongum_NCC2705", "B. longum", "Bifidobacteriales"),
    ("Mluteus_NCTC2665", "M. luteus", "Micrococcales"),
    ("Ecoli_K12_MG1655", "E. coli", "Proteobacteria (control)"),
    ("Bsubtilis_168", "B. subtilis", "Firmicutes (control)"),
]
GENES = ["eno", "divIC", "Rv1025", "ppx2"]
COL = {"eno": "#4C72B0", "divIC": "#DD8452", "Rv1025": "#C44E52", "ppx2": "#55A868"}

def parse_order(s):
    """'eno@123 divIC@456 ...' -> {gene: midpoint} (gene names as in GENES)."""
    d = {}
    for m in re.finditer(r"(\w+)@(\d+)", s):
        g = m.group(1)
        d[g] = int(m.group(2))
    return d

def main():
    data = {r["genome"]: r for r in csv.DictReader(open(TSV), delimiter="\t")}
    fig, ax = plt.subplots(figsize=(9.5, 7.2))
    n = len(ROWS)
    aw = 0.34  # demi-hauteur d'arrow
    for i, (key, label, clade) in enumerate(ROWS):
        y = n - i
        r = data.get(key, {})
        present = r.get("Rv1025_present", "NON").startswith("OUI")
        ax.text(-0.02, y, label, ha="right", va="center", fontstyle="italic", fontsize=10.5)
        ax.text(1.14, y, clade, ha="left", va="center", fontsize=8, color="grey")
        if not present:
            ax.text(0.5, y, "DUF501 absent", ha="center", va="center", fontsize=9.5,
                    color="#888", fontstyle="italic",
                    bbox=dict(boxstyle="round,pad=0.3", fc="#f0f0f0", ec="#ccc"))
            continue
        pos = parse_order(r["order_on_contig"])
        if not all(g in pos for g in GENES):
            continue
        # normaliser l'orientation : opéron lu eno(min)->ppx2(max)
        reverse = pos["eno"] > pos["ppx2"]
        xs = {g: pos[g] for g in GENES}
        lo, hi = min(xs.values()), max(xs.values())
        span = hi - lo or 1
        def X(g):
            x = (xs[g] - lo) / span
            return (1 - x) if reverse else x
        # dessiner les arrows (largeur fixe, sens = orientation de lecture)
        w = 0.07
        for g in GENES:
            xc = 0.06 + X(g) * 0.78
            hl = w * 0.55
            ax.add_patch(FancyArrow(xc - w/2, y, w, 0, width=aw*0.9, head_width=aw*1.35,
                                    head_length=hl, length_includes_head=True,
                                    fc=COL[g], ec="black" if g == "Rv1025" else "none",
                                    lw=1.6 if g == "Rv1025" else 0, zorder=3))
        # % identité Rv1025
        m = re.search(r"OUI \((\d+)%", r["Rv1025_present"])
        if m:
            ax.text(1.00, y, f"{m.group(1)}%", ha="center", va="center", fontsize=8.5,
                    color="#C44E52", fontweight="bold")
    ax.set_xlim(-0.30, 1.55)
    ax.set_ylim(0.3, n + 1.4)
    ax.axis("off")
    # séparateur avant les témoins
    ax.axhline(2.5, xmin=0.16, xmax=0.83, color="#bbb", ls="--", lw=0.8)
    # légende
    handles = [Patch(fc=COL[g], ec="black" if g == "Rv1025" else "none",
                     lw=1.4 if g == "Rv1025" else 0,
                     label={"eno": "eno (enolase)", "divIC": "divIC (FtsB/DivIC)",
                            "Rv1025": "Rv1025 (DUF501)", "ppx2": "ppx2 (exopolyphosphatase)"}[g])
               for g in GENES]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.02),
              ncol=4, frameon=False, fontsize=9, handlelength=1.4, columnspacing=1.3)
    ax.text(1.00, n + 0.55, "Rv1025 id.", ha="center", va="center", fontsize=8, color="#C44E52")
    plt.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT + ".pdf", bbox_inches="tight")
    print("écrit :", OUT + ".png / .pdf")

if __name__ == "__main__":
    main()
