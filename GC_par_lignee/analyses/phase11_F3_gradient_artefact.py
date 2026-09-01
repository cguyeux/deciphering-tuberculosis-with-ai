#!/usr/bin/env python3
"""
Objet       : F3 -- figure avant/apres de l'artefact du gradient naif (fig_plan.md,
              archetype C1). Panneau gauche : rapport pertes/gains NAIF, compte
              par souche (pseudo-replique, pre-polarisation), contre la distance
              a H37Rv (nsub) -- reproduit le gradient apparent de A7 (0,89 a
              1,38). Panneau droit : rapport pertes/gains POLARISE par evenement
              contre la profondeur d'arbre (proxy k = nombre de porteurs), pour
              trois lignees reperes (L1, L6.1, L9) -- la courbe en marche
              d'escalier de A9, plate une fois l'orientation corrigee et la
              pseudo-replication levee.
Entrees     : résultats/phase11_F3_naif_par_souche.tsv (ce phase, panneau gauche)
              résultats/phase2_polarisation_mtbc0_n40.tsv (P2.4, panneau droit)
Sorties     : article/figures/fig_F3_gradient_artefact.pdf
              article/figures/fig_F3_gradient_artefact.png
Reutilisable: non -- specifique a cette figure
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
LIGNEES_REPERES = ["L1", "L6.1", "L9"]
DEPTH_ORDER = ["terminal", "k2-8", "k9-15", "k16-22", "k23-29", "k30-36", "k37-39", "racine"]
DEPTH_ORDER_L9 = ["terminal", "k2-7", "k8-13", "k14-19", "k20-25", "k26-31", "k32-35", "racine"]


def depth_rank(clade: str, k_class: str) -> int:
    order = DEPTH_ORDER_L9 if clade == "L9" else DEPTH_ORDER
    return order.index(k_class) if k_class in order else len(order)


def main() -> int:
    naif = pd.read_csv(ROOT / "résultats" / "phase11_F3_naif_par_souche.tsv", sep="\t")
    pol = pd.read_csv(ROOT / "résultats" / "phase2_polarisation_mtbc0_n40.tsv", sep="\t")
    pol = pol[pol["clade"].isin(LIGNEES_REPERES)].copy()
    pol["depth_rank"] = pol.apply(lambda r: depth_rank(r["clade"], r["k_class"]), axis=1)
    pol = pol.sort_values(["clade", "depth_rank"])

    fig, axes = fs.panel_grid(PRESET, nrows=1, ncols=2, width_mm=183, wspace=0.32)
    ax_left, ax_right = axes

    # --- Panneau gauche : gradient naif, comptage par souche (pseudo-replique)
    color_naif = "0.35"
    ax_left.scatter(
        naif["nsub"], naif["ratio"],
        s=6, alpha=0.35, color=color_naif, edgecolors="none", zorder=2,
    )
    coeffs = np.polyfit(naif["nsub"], naif["ratio"], 1)
    xs = np.linspace(naif["nsub"].min(), naif["nsub"].max(), 100)
    ax_left.plot(xs, np.polyval(coeffs, xs), color=fs.PALETTE_CATEGORICAL[3],
                 linewidth=1.4, zorder=3)
    ax_left.axhline(1.0, color="0.75", linestyle=":", linewidth=0.8, zorder=1)
    ax_left.set_xlabel("Distance from H37Rv (substitutions per strain)")
    ax_left.set_ylabel("G:C loss/gain ratio (per strain, naive)")
    ax_left.set_title("Naive count: apparent gradient", fontsize=fs.get_preset(PRESET)["subtitle_size"])

    # --- Panneau droit : rapport polarise par evenement, vs profondeur d'arbre
    for i, clade in enumerate(LIGNEES_REPERES):
        sub = pol[pol["clade"] == clade]
        color = fs.PALETTE_CATEGORICAL[i]
        ax_right.errorbar(
            sub["depth_rank"], sub["ratio_pol"],
            yerr=[sub["ratio_pol"] - sub["ci95_lo"], sub["ci95_hi"] - sub["ratio_pol"]],
            fmt="o-", color=color, ecolor=color, elinewidth=0.9, capsize=2.0,
            markersize=3.4, linewidth=1.1, label=clade, zorder=3,
        )
    ax_right.axhline(1.0, color="0.75", linestyle=":", linewidth=0.8, zorder=1)
    n_ticks = len(DEPTH_ORDER)
    ax_right.set_xticks(range(n_ticks))
    ax_right.set_xticklabels(
        ["terminal", "", "", "", "", "", "", "root"], fontsize=fs.get_preset(PRESET)["note_size"],
    )
    ax_right.set_xlabel("Tree depth (terminal → pool root)")
    ax_right.set_ylabel("G:C loss/gain ratio (per event, polarized)")
    ax_right.set_title("Polarized count: flat ratio", fontsize=fs.get_preset(PRESET)["subtitle_size"])
    ax_right.legend(frameon=False, fontsize=fs.get_preset(PRESET)["note_size"], loc="lower left")

    fs.panel_labels(axes)
    fs.save(fig, str(ROOT / "article" / "figures" / "fig_F3_gradient_artefact"), PRESET, raster=True)
    print("Ecrit : article/figures/fig_F3_gradient_artefact.pdf (+ .png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
