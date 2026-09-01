#!/usr/bin/env python3
"""
Objet       : S2 -- heatmap complete des 96 canaux (6 classes de substitution
              x 16 contextes trinucleotidiques 5'.3') de leur contribution a
              l'heterogeneite contextuelle inter-lignees (G2 d'independance,
              P8.1/P8.2). Mention breve dans le corps (§3.12), figure complete
              releguee en supplementary (fig_plan.md, Supplementary S2) : la
              dissociation gains/pertes reste INDETERMINEE (P8.3, A32) et la
              piste est refermee sans reouverture prevue.
Entrees     : résultats/phase7_p81_canaux_classes.tsv (index = canal
              "CLASSE@5'.3'", colonnes g2_contribution, part, reportable, n)
Sorties     : article/figures/fig_S2_heatmap_contextuel.pdf
              article/figures/fig_S2_heatmap_contextuel.png
Reutilisable: non -- specifique a ce jeu de resultats et a cette figure
Projet      : GC_par_lignee
Date        : 2026-09-01
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT.parent.parent / "claude_plugins" / "bio_pathogens" / "skills" / "sci-figure" / "scripts"))
import figstyle as fs  # noqa: E402

PRESET = "nature_double"
BASES = "ACGT"
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
LOSS = {"C>A", "C>T"}
CTX = [f"{a}.{b}" for a in BASES for b in BASES]  # 5'.3', ordre fixe


def main() -> int:
    df = pd.read_csv(ROOT / "résultats" / "phase7_p81_canaux_classes.tsv", sep="\t")
    df = df.rename(columns={"Unnamed: 0": "canal"})
    df[["classe", "contexte"]] = df["canal"].str.split("@", expand=True)

    mat = df.pivot(index="classe", columns="contexte", values="g2_contribution")
    mat = mat.loc[CLASSES, CTX]
    rep = df.pivot(index="classe", columns="contexte", values="reportable").loc[CLASSES, CTX]

    fig, ax = fs.new_figure(PRESET, height_mm=95, width_mm=183)

    im = ax.imshow(mat.values, cmap="magma", aspect="auto", vmin=0)

    # Croix sur les canaux NON reportables individuellement (C1 : < 100
    # evenements toutes lignees, ou < 5 dans au moins 13/17 lignees).
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if not bool(rep.values[i, j]):
                ax.plot(j, i, marker="x", color="white", markersize=4, markeredgewidth=0.8)

    ax.set_xticks(range(len(CTX)))
    ax.set_xticklabels(CTX, rotation=90, fontsize=5.5, family="monospace")
    ax.set_yticks(range(len(CLASSES)))
    ylabels = [f"{c}{'  (loss)' if c in LOSS else '  (gain)'}" for c in CLASSES]
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel("trinucleotide context (5' flank . 3' flank)")

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
    cbar.set_label(r"$G^2$ contribution", fontsize=6.5)
    cbar.ax.tick_params(labelsize=5.5)

    ax.text(
        0.0, -0.30,
        "x = channel not individually reportable (C1: < 100 events pooled, "
        "or < 5 events in $\\geq$13/17 lineages)",
        transform=ax.transAxes, fontsize=6, color="0.25", va="top", ha="left",
    )

    fs.save(fig, str(ROOT / "article" / "figures" / "fig_S2_heatmap_contextuel"), PRESET, raster=True)
    print("Ecrit : article/figures/fig_S2_heatmap_contextuel.pdf (+ .png)")
    top5 = df.sort_values("g2_contribution", ascending=False).head(5)
    print("5 canaux les plus contributeurs :")
    print(top5[["canal", "g2_contribution", "part", "reportable", "n"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
