#!/usr/bin/env python3
"""
Objet       : F6 -- encart scatter du controle par replichore (fig_plan.md,
              archetype M1 + encart de donnees, P11.2.6). v/u(replichore 1)
              contre v/u(replichore 2), une lignee = un point, pour les
              lignees reportables dans les deux replichores (>= 20 pertes ET
              >= 20 gains, C1 de P3.4). Destine a etre embarque comme encart
              central dans le schema circulaire fig_F6_mecanisme_biais_brin.tex
              (P11.2.6), pas publie seul.
Entrees     : résultats/phase10_p34_vu_replichore.tsv (v/u par lignee x
              replichore, 3.4)
              résultats/phase10_p34_test_replichore.tsv (rho, p, verdict, 3.4)
Sorties     : article/figures/fig_F6_encart_scatter.pdf
              article/figures/fig_F6_encart_scatter.png (revue locale, non
              embarque dans le PDF final)
Reutilisable: non -- specifique a cette figure
Projet      : GC_par_lignee
Date        : 2026-08-31
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT.parent.parent / "claude_plugins" / "bio_pathogens" / "skills" / "sci-figure" / "scripts"))
import figstyle as fs  # noqa: E402

PRESET = "nature_single"
INSET_SIDE_MM = 50.0


def main() -> int:
    vu = pd.read_csv(ROOT / "résultats" / "phase10_p34_vu_replichore.tsv", sep="\t")
    test = pd.read_csv(ROOT / "résultats" / "phase10_p34_test_replichore.tsv", sep="\t").iloc[0]

    vu = vu[vu["reportable"] == True]  # noqa: E712 -- valeur bool ecrite "True" en TSV
    piv = vu.pivot(index="lignee", columns="replichore", values="vu").dropna()
    piv.columns = ["repl1", "repl2"]
    n = len(piv)
    assert n == int(test["n_lignees"]), (
        f"decompte encart ({n}) != decompte du test primaire P3.4 ({test['n_lignees']})"
    )

    fig, ax = fs.panel_grid(PRESET, nrows=1, ncols=1,
                             width_mm=INSET_SIDE_MM, height_mm=INSET_SIDE_MM)

    lo = min(piv["repl1"].min(), piv["repl2"].min()) * 0.92
    hi = max(piv["repl1"].max(), piv["repl2"].max()) * 1.08
    ax.plot([lo, hi], [lo, hi], color="0.7", linestyle=":", linewidth=0.9, zorder=1)
    ax.scatter(piv["repl1"], piv["repl2"], s=13, color=fs.PALETTE_CATEGORICAL[0],
               edgecolors="white", linewidths=0.3, zorder=3)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("v/u (replichore 1)", fontsize=fs.get_preset(PRESET)["label_size"])
    ax.set_ylabel("v/u (replichore 2)", fontsize=fs.get_preset(PRESET)["label_size"])
    ax.tick_params(labelsize=fs.get_preset(PRESET)["note_size"])
    rho = test["rho"]
    p = test["p_rho"]
    mantissa, exponent = f"{p:.1e}".split("e")
    p_str = f"{mantissa}\\times10^{{{int(exponent)}}}"
    ax.text(
        0.05, 0.95,
        f"$\\rho$ = {rho:.3f}\n$n$ = {n}\n$p$ = ${p_str}$",
        transform=ax.transAxes, ha="left", va="top",
        fontsize=fs.get_preset(PRESET)["note_size"],
    )

    stem = ROOT / "article" / "figures" / "fig_F6_encart_scatter"
    fs.save(fig, str(stem), PRESET, raster=True)
    print(f"Ecrit : {stem.relative_to(ROOT)}.pdf/.png ({n} lignees, rho={rho:.4f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
