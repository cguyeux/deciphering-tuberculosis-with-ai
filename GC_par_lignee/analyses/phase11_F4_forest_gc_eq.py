#!/usr/bin/env python3
"""
Objet       : F4 -- forest plot du GC d'equilibre mutationnel (GC_eq) par
              lignee, figure de resultat principal de l'article (fig_plan.md).
              17 lignees triees par GC_eq croissant, IC bootstrap hierarchique
              (A23), ligne de reference au GC observe de H37Rv (65,6147 %).
              Marqueur neutre uniforme (pas de code humain/animal, tranche par
              I4 : le statut d'hote n'organise pas le classement).
Entrees     : résultats/phase5_p52_ic_bootstrap.tsv (colonnes gc_eq,
              gc_eq_boot_lo, gc_eq_boot_hi, n_souches -- panel final 17
              lignees / 519 souches, P5.3)
Sorties     : article/figures/fig_F4_gc_eq_forest.pdf
              article/figures/fig_F4_gc_eq_forest.png
Reutilisable: non -- specifique a ce jeu de resultats et a cette figure
Projet      : GC_par_lignee
Date        : 2026-08-31
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT.parent.parent / "claude_plugins" / "bio_pathogens" / "skills" / "sci-figure" / "scripts"))
import figstyle as fs  # noqa: E402

GC_OBSERVE_H37RV = 65.6147  # %, cite dans etat_des_decouvertes.md et CLAUDE.md

PRESET = "nature_double"


def main() -> int:
    df = pd.read_csv(ROOT / "résultats" / "phase5_p52_ic_bootstrap.tsv", sep="\t")
    df = df.sort_values("gc_eq").reset_index(drop=True)

    gc_eq = df["gc_eq"] * 100
    lo = df["gc_eq_boot_lo"] * 100
    hi = df["gc_eq_boot_hi"] * 100
    y = range(len(df))

    fig, ax = fs.new_figure(PRESET, height_mm=140)

    ax.axvline(GC_OBSERVE_H37RV, color="0.55", linestyle="--", linewidth=0.9, zorder=1)
    ax.text(
        GC_OBSERVE_H37RV, len(df) - 0.3, "Observed GC\n(H37Rv, 65.61 %)",
        ha="center", va="bottom", fontsize=6.5, color="0.35",
    )

    ax.errorbar(
        gc_eq, y,
        xerr=[gc_eq - lo, hi - gc_eq],
        fmt="o", color=fs.PALETTE_CATEGORICAL[0], ecolor="0.25",
        elinewidth=1.0, capsize=2.2, markersize=4.2, zorder=3,
    )

    labels = [f"{l} (n={n})" for l, n in zip(df["lignee"], df["n_souches"])]
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.7, len(df) - 0.3)
    ax.set_xlabel("Mutational equilibrium GC, GC$_{eq}$ = u/(u+v) (%)")
    ax.set_xlim(15, 70)

    fs.save(fig, str(ROOT / "article" / "figures" / "fig_F4_gc_eq_forest"), PRESET, raster=True)
    print("Ecrit : article/figures/fig_F4_gc_eq_forest.pdf (+ .png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
