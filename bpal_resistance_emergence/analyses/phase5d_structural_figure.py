#!/usr/bin/env python3
"""Phase 5d (réponse review R09) — figure composite : filtre structural + contrôle positif de convergence.

Deux panneaux, données curées des expériences structurales (experiments/2026-06-14_structural_f420/RESULTS.md)
et du scan de convergence (phase5_convergence_positions.tsv) :

  (A) DISSOCIATION STRUCTURALE. Distance atome-atome au cofacteur F420 (x) vs RSA (y) pour chaque résidu
      candidat/marqueur, sur les structures holo Fgd1 (3B4Y) et Ddn (3R5W). Zones ombrées : poche
      (<8 Å) et enfoui (RSA<0.20). Couleur = classe (convergent / phénotype-cas-témoins / marqueur de
      lignée / ancre cataloguée). Forme = protéine. Montre que les candidats convergents bordent la poche
      enfouie tandis que les marqueurs de lignée sont distaux/exposés.

  (B) CONTRÔLE POSITIF DE CONVERGENCE. Nombre de lignées majeures indépendantes frappées (barres) et
      nombre de résidus alternatifs distincts (annotation) pour les sites FQ connus (contrôle gyrA/gyrB)
      vs les positions convergentes F420 du projet. Le détecteur récupère la convergence FQ établie à un
      niveau comparable, validant la méthode avant d'interpréter les hits F420.

Sortie : article/figures/phase5d_structural_convergence.{png,pdf} (600 DPI).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths

# --- données panneau A : (label, protéine, dist_F420_Å, RSA, classe) ---------
# classes : conv = convergent candidat ; pheno = cas-témoins DLM ; marker = marqueur de lignée ; anchor = catalogué-R
RESIDUES = [
    ("fgd1 A10",  "Fgd1", 7.6, 0.00, "conv"),
    ("fgd1 A210", "Fgd1", 7.8, 0.04, "conv"),
    ("ddn A111",  "Ddn",  13.1, 0.06, "conv"),
    ("ddn Y65",   "Ddn",  3.1, 0.27, "pheno"),
    ("ddn G81",   "Ddn",  7.0, 0.19, "pheno"),
    ("ddn L49",   "Ddn",  7.7, 0.10, "anchor"),
    ("fgd1 R64",  "Fgd1", 21.4, 0.41, "marker"),
    ("fgd1 K270", "Fgd1", 18.6, 0.59, "marker"),
    ("fgd1 K296", "Fgd1", 17.2, 0.29, "marker"),
    ("ddn D113",  "Ddn",  12.1, 0.27, "marker"),
]
CLS = {
    "conv":   ("#d62728", "convergent candidate"),
    "pheno":  ("#ff7f0e", "phenotype case-control"),
    "marker": ("#7f7f7f", "lineage marker (benign)"),
    "anchor": ("#1f77b4", "catalogued DLM-R anchor"),
}
SHAPE = {"Fgd1": "o", "Ddn": "^"}

# --- données panneau B : (label, n_lignées_indép, n_alts_distincts, classe) ---
# contrôles FQ (phase5 convergence) + candidats F420
CONV = [
    ("gyrA A90 (control)", 7, 4, "control"),
    ("gyrB QRDR (control)", 4, 5, "control"),
    ("fgd1 210", 4, 5, "f420"),
    ("fgd1 10",  4, 3, "f420"),
    ("ddn 111",  4, 3, "f420"),
]


def main():
    figdir = paths.ARTICLE / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.35, 1]})

    # ===== Panneau A : distance F420 vs RSA =====
    axA.axvspan(0, 8, color="#fde0dc", alpha=0.5, zorder=0)          # poche (1re couronne)
    axA.axhspan(0, 0.20, color="#e8e8e8", alpha=0.5, zorder=0)        # enfoui
    axA.axvline(4, ls=":", color="#555", lw=0.8, zorder=1)
    axA.axvline(8, ls="--", color="#999", lw=0.8, zorder=1)
    axA.axvline(12, ls="--", color="#ccc", lw=0.8, zorder=1)
    for label, prot, d, rsa, cls in RESIDUES:
        col = CLS[cls][0]
        axA.scatter(d, rsa, s=95, marker=SHAPE[prot], color=col, edgecolor="black",
                    linewidth=0.6, zorder=3)
        dy = 0.025 if label not in ("ddn L49", "fgd1 A210") else -0.045
        ha = "left" if d < 18 else "right"
        axA.annotate(label, (d, rsa), xytext=(d + (0.5 if ha == "left" else -0.5), rsa + dy),
                     fontsize=7, ha=ha, va="bottom", zorder=4)
    axA.set_xlim(0, 24)
    axA.set_ylim(-0.03, 0.66)
    axA.set_xlabel("distance to F420 cofactor (Å)", fontsize=10)
    axA.set_ylabel("relative solvent accessibility (RSA)", fontsize=10)
    axA.set_title("A  Structural dissociation (holo Fgd1 3B4Y, Ddn 3R5W)", fontsize=10, loc="left")
    axA.text(4.0, 0.63, "pocket\n(1st shell)", fontsize=7, color="#b03020", ha="center")
    axA.text(20, 0.02, "buried", fontsize=7, color="#555")
    cls_handles = [Line2D([0], [0], marker="s", ls="", color=c, markeredgecolor="black",
                          label=lab, markersize=8) for c, lab in CLS.values()]
    shp_handles = [Line2D([0], [0], marker=m, ls="", color="white", markeredgecolor="black",
                          label=p, markersize=8) for p, m in SHAPE.items()]
    leg1 = axA.legend(handles=cls_handles, fontsize=6.8, loc="upper left", frameon=True, framealpha=0.9)
    axA.add_artist(leg1)
    axA.legend(handles=shp_handles, fontsize=6.8, loc="lower right", frameon=True, framealpha=0.9,
               title="protein", title_fontsize=6.8)

    # ===== Panneau B : convergence positive control =====
    labels = [c[0] for c in CONV]
    nlin = [c[1] for c in CONV]
    nalt = [c[2] for c in CONV]
    cols = ["#2ca02c" if c[3] == "control" else "#d62728" for c in CONV]
    y = range(len(CONV))
    axB.barh(list(y), nlin, color=cols, edgecolor="black", linewidth=0.6, height=0.62, zorder=3)
    for i, (n, a) in enumerate(zip(nlin, nalt)):
        axB.text(n + 0.08, i, f"{a} alt. residues", va="center", fontsize=7)
    axB.set_yticks(list(y))
    axB.set_yticklabels(labels, fontsize=8)
    axB.invert_yaxis()
    axB.set_xlim(0, 9)
    axB.set_xlabel("independent major lineages struck", fontsize=10)
    axB.set_title("B  Positional convergence (positive control)", fontsize=10, loc="left")
    axB.legend(handles=[Line2D([0], [0], marker="s", ls="", color="#2ca02c", markeredgecolor="black",
                               label="known FQ sites (control)", markersize=8),
                        Line2D([0], [0], marker="s", ls="", color="#d62728", markeredgecolor="black",
                               label="F420 candidates", markersize=8)],
               fontsize=6.8, loc="lower right", frameon=True)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(figdir / f"phase5d_structural_convergence.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"figure -> {figdir}/phase5d_structural_convergence.png")


if __name__ == "__main__":
    main()
