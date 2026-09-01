#!/usr/bin/env python3
"""
Objet       : S1 -- demontrer par une figure le mecanisme derriere le retrait
              de R4 (densite d'echantillonnage fabrique un faux signal de
              selection) : a Ne vrai identique (meme lignee, meme sous-clade),
              un pool densement echantillonne deplace le pN/pS terminal, plus
              que tout l'IQR inter-lignees observe par ailleurs dans le
              projet (A26). fig_plan.md, Supplementary S1.
Entrees     : résultats/phase6_p101_dense_disperse.tsv (colonnes lignee,
              sous_clade_dense, pn_ps_term_dense, pn_ps_term_disperse, ...
              8 lignees appariees dense/disperse)
Sorties     : article/figures/fig_S1_dense_disperse.pdf
              article/figures/fig_S1_dense_disperse.png
Reutilisable: non -- specifique a ce jeu de resultats et a cette figure
Projet      : GC_par_lignee
Date        : 2026-09-01
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT.parent.parent / "claude_plugins" / "bio_pathogens" / "skills" / "sci-figure" / "scripts"))
import figstyle as fs  # noqa: E402

PRESET = "nature_single"


def main() -> int:
    df = pd.read_csv(ROOT / "résultats" / "phase6_p101_dense_disperse.tsv", sep="\t")
    df["delta"] = df.pn_ps_term_dense - df.pn_ps_term_disperse
    df = df.sort_values("delta").reset_index(drop=True)
    n = len(df)
    n_same_sign = int((df.delta > 0).sum())

    fig, ax = fs.new_figure(PRESET, height_mm=100, width_mm=89)

    y = np.arange(n)
    colors = fs.PALETTE_CATEGORICAL

    for i, row in df.iterrows():
        color = "0.35" if row.delta > 0 else colors[3]
        ax.plot(
            [row.pn_ps_term_disperse, row.pn_ps_term_dense], [i, i],
            color=color, linewidth=1.0, zorder=2, alpha=0.8,
        )
    ax.scatter(
        df.pn_ps_term_disperse, y, marker="o", s=22, color=colors[0],
        edgecolor="black", linewidth=0.4, zorder=3, label="dispersed pool",
    )
    ax.scatter(
        df.pn_ps_term_dense, y, marker="D", s=22, color=colors[1],
        edgecolor="black", linewidth=0.4, zorder=3, label="dense pool (same sub-clade)",
    )

    ax.set_yticks(y)
    ax.set_yticklabels(df.lignee)
    ax.set_xlabel(r"Terminal-branch $\pi_N/\pi_S$")
    ax.axvline(1.0, color="0.7", linestyle=":", linewidth=0.8, zorder=0)
    ax.set_ylim(-0.7, n - 0.3)

    ax.text(
        0.02, -0.16,
        f"{n_same_sign}/{n} pairs: dense pool inflates $\\pi_N/\\pi_S$ "
        f"(median $\\Delta$ = {df.delta.median():+.3f})",
        transform=ax.transAxes, fontsize=6.5, color="0.25", va="top", ha="left",
    )

    ax.legend(loc="lower right", fontsize=6, frameon=False)

    fs.save(fig, str(ROOT / "article" / "figures" / "fig_S1_dense_disperse"), PRESET, raster=True)
    print("Ecrit : article/figures/fig_S1_dense_disperse.pdf (+ .png)")
    print(f"n paires = {n} ; {n_same_sign}/{n} dans le meme sens ; "
          f"delta median = {df.delta.median():+.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
