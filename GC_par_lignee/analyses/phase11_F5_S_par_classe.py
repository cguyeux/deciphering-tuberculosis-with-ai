#!/usr/bin/env python3
"""
Objet       : F5 -- spaghetti plot de la force de maintien du GC, S =
              ln[(GC/(1-GC)) * v/u], a travers quatre classes de sites de
              contrainte selective croissante (4-fold, intergenique, non
              degenere, 2e position de codon). Resultat le plus fort et le
              moins attendu du projet (fig_plan.md, F5) : S est MAXIMAL aux
              sites 4-fois degeneres (aucune contrainte proteique) et MINIMAL
              aux 2emes positions (contrainte maximale), 12 lignees sur 12
              dans ce sens -- l'inverse de ce que predirait une selection sur
              l'acide amine (A37-A39). Falsifie directement le modele nul
              "S est porte par la selection proteique".
Entrees     : résultats/phase9_p42_force_par_classe.tsv (colonnes classe,
              lignee, S -- restreint aux 12 lignees ayant les 4 classes,
              La4 exclue faute de donnees 4fold/intergenique, cf. n_lignees=12
              dans phase9_p42_qualification.tsv)
Sorties     : article/figures/fig_F5_S_par_classe.pdf
              article/figures/fig_F5_S_par_classe.png
Reutilisable: non -- specifique a ce jeu de resultats et a cette figure
Projet      : GC_par_lignee
Date        : 2026-08-31
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT.parent.parent / "claude_plugins" / "bio_pathogens" / "skills" / "sci-figure" / "scripts"))
import figstyle as fs  # noqa: E402

PRESET = "nature_double"

# Ordre d'abscisse = ordre decroissant de S median (contrainte selective
# croissante), tel que decrit dans fig_plan.md F5.
CLASSES = ["4fold", "intergenique", "0fold", "0fold_pos2"]
LABELS = {
    "4fold": "four-fold\ndegenerate",
    "intergenique": "intergenic",
    "0fold": "non-\ndegenerate",
    "0fold_pos2": "2nd\nposition",
}

MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<", ">", "p"]


def main() -> int:
    df = pd.read_csv(ROOT / "résultats" / "phase9_p42_force_par_classe.tsv", sep="\t")
    piv = df.pivot(index="lignee", columns="classe", values="S")
    piv = piv.dropna(subset=CLASSES)  # La4 exclue (pas de 4fold/intergenique)
    piv = piv[CLASSES]
    n = len(piv)

    fig, ax = fs.new_figure(PRESET, height_mm=118, width_mm=140)

    x = np.arange(len(CLASSES))

    # Etendue (min-max) en bande grise, ancre visuelle avant les lignes.
    lo = piv.min(axis=0).values
    hi = piv.max(axis=0).values
    ax.fill_between(x, lo, hi, color="0.85", zorder=1, label="range (min-max)")

    # Une ligne par lignee, couleur x forme cyclees pour rester lisible en N&B.
    colors = fs.PALETTE_CATEGORICAL
    for i, (lignee, row) in enumerate(piv.iterrows()):
        color = colors[i % len(colors)]
        linestyle = "-" if (i // len(colors)) % 2 == 0 else "--"
        marker = MARKERS[i % len(MARKERS)]
        ax.plot(
            x, row.values, color=color, linestyle=linestyle, linewidth=1.1,
            marker=marker, markersize=4.0, alpha=0.85, zorder=3,
            label=lignee,
        )

    # Mediane en gras, au-dessus des lignes individuelles.
    med = piv.median(axis=0).values
    ax.plot(
        x, med, color="black", linewidth=2.2, marker="D", markersize=5.5,
        zorder=4, label="median (n=%d)" % n,
    )

    ax.axhline(0, color="0.6", linestyle=":", linewidth=0.8, zorder=0)

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[c] for c in CLASSES])
    ax.set_ylabel(r"GC maintenance force, $S = \ln\left[\frac{GC}{1-GC}\cdot\frac{v}{u}\right]$")
    ax.set_xlim(-0.3, len(CLASSES) - 0.7)

    ax.text(
        0.02, 0.03,
        f"{n}/{n} lineages: S(4-fold) > S(2nd position)",
        transform=ax.transAxes, fontsize=6.5, color="0.25", va="bottom", ha="left",
    )

    handles, labels_ = ax.get_legend_handles_labels()
    order = [labels_.index(l) for l in labels_ if l not in ("median (n=%d)" % n, "range (min-max)")]
    lignee_handles = [handles[i] for i in order]
    lignee_labels = [labels_[i] for i in order]
    special_handles = [handles[labels_.index("range (min-max)")], handles[labels_.index("median (n=%d)" % n)]]
    special_labels = ["range (min-max)", "median (n=%d)" % n]

    leg1 = ax.legend(
        lignee_handles, lignee_labels, loc="upper left",
        bbox_to_anchor=(1.01, 1.0), fontsize=5.5, ncol=1, frameon=False,
        title="lineage", title_fontsize=6,
    )
    ax.add_artist(leg1)
    ax.legend(
        special_handles, special_labels, loc="lower left",
        bbox_to_anchor=(1.01, 0.0), fontsize=5.5, frameon=False,
    )

    fs.save(fig, str(ROOT / "article" / "figures" / "fig_F5_S_par_classe"), PRESET, raster=True)
    print("Ecrit : article/figures/fig_F5_S_par_classe.pdf (+ .png)")
    print(f"n lignees (4 classes completes) = {n} ; classes exclues faute de "
          f"donnees : {sorted(set(df.lignee) - set(piv.index))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
