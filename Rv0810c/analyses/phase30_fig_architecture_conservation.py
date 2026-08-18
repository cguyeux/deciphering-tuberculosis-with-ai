#!/usr/bin/env python3
"""Figures 1 et 2 du manuscrit : architecture bipartite (pLDDT/Rg) et
conservation MSA par position. Piste P11.4."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RESULTATS = ROOT / "résultats"
FIGDIR = ROOT / "article" / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

SKILL_SCRIPTS = Path(
    "/home/christophe/docs/codes/claude_plugins/bio_pathogens/skills/sci-figure/scripts"
)
sys.path.insert(0, str(SKILL_SCRIPTS))
import figstyle as fs  # noqa: E402

PRESET = "generic"
LANG = sys.argv[1] if len(sys.argv) > 1 else "en"

TXT = {
    "en": dict(
        windows=[(1, 19, "rigid module (1–19)"), (20, 33, "hinge (20–33)"),
                  (34, 60, "acidic disordered tail (34–60)")],
        xlabel_a="Residue position", ylabel_a="AlphaFold pLDDT",
        bar_labels=["Rv0810c\n(observed)", "compact globular\n60-residue protein\n(expected)"],
        ylabel_b="Radius of gyration (Å)",
        suffix="",
        xlabel2="Residue position (H37Rv numbering)",
        ylabel2="Majority-residue frequency (%)\n1,930 PF11273 sequences",
        boundary_label="order→disorder\nboundary (E32/L33)",
        legend2=["non-basic residue", "basic residue (K/R)",
                  "perfectly invariant (100%, Shannon = 0)"],
    ),
    "fr": dict(
        windows=[(1, 19, "module rigide (1–19)"), (20, 33, "charnière (20–33)"),
                  (34, 60, "queue désordonnée acide (34–60)")],
        xlabel_a="Position du résidu", ylabel_a="pLDDT AlphaFold",
        bar_labels=["Rv0810c\n(observé)", "protéine globulaire\ncompacte de 60 résidus\n(attendu)"],
        ylabel_b="Rayon de giration (Å)",
        suffix="_fr",
        xlabel2="Position du résidu (numérotation H37Rv)",
        ylabel2="Fréquence du résidu majoritaire (%)\n1930 séquences PF11273",
        boundary_label="frontière\nordre→désordre (E32/L33)",
        legend2=["résidu non basique", "résidu basique (K/R)",
                  "parfaitement invariant (100 %, Shannon = 0)"],
    ),
}[LANG]


def fmt1(x):
    s = f"{x:.1f}"
    return s.replace(".", ",") if LANG == "fr" else s


def save_600dpi(fig, stem):
    written = fs.save(fig, stem, PRESET)
    png = Path(stem).with_suffix(".png")
    fig.savefig(png, format="png", dpi=600)
    written.append(png)
    return written


def load_pdb_plddt(pdb_path):
    plddt, coords = {}, []
    for line in open(pdb_path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            resnum = int(line[22:26])
            b = float(line[60:66])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            plddt[resnum] = b
            coords.append((resnum, x, y, z))
    return plddt, coords


def radius_of_gyration(coords):
    n = len(coords)
    cx = sum(c[1] for c in coords) / n
    cy = sum(c[2] for c in coords) / n
    cz = sum(c[3] for c in coords) / n
    rg2 = sum((c[1] - cx) ** 2 + (c[2] - cy) ** 2 + (c[3] - cz) ** 2 for c in coords) / n
    return rg2 ** 0.5


def fig1_architecture():
    pdb_path = RESULTATS / "p6_3_pocket_detection" / "rv0810c.pdb"
    plddt, coords = load_pdb_plddt(pdb_path)
    residues = sorted(plddt)
    values = [plddt[r] for r in residues]
    rg = radius_of_gyration(coords)

    fig, axes = fs.panel_grid(PRESET, nrows=1, ncols=2, width_ratios=[2.3, 1])
    ax_p, ax_rg = axes

    windows = TXT["windows"]
    colors = [fs.PALETTE_CATEGORICAL[0], fs.PALETTE_CATEGORICAL[2], fs.PALETTE_CATEGORICAL[1]]
    for (start, end, label), col in zip(windows, colors):
        ax_p.axvspan(start - 0.5, end + 0.5, color=col, alpha=0.15, lw=0)
        xs = [r for r in residues if start <= r <= end]
        ys = [plddt[r] for r in xs]
        mean_y = float(np.mean(ys))
        ax_p.plot(xs, ys, color=col, lw=1.4)
        ax_p.hlines(mean_y, start, end, color=col, lw=1.0, linestyle="--")
        # Écrit dans la teinte de la bande, le label de moyenne se lisait mal sur
        # l'axvspan du même ton (le point décimal disparaissait : « 91 9 » au lieu
        # de « 91,9 », /fig-check du 2026-08-17). On garde l'identité de couleur
        # mais on l'assombrit pour retrouver du contraste sur le fond teinté.
        # Ancré en FIN de fenêtre, pas au début : au début, le label de la charnière
        # tombait sous le triangle T24 qui le chevauchait (/fig-check du 2026-08-17).
        ax_p.annotate(fmt1(mean_y), xy=(end, mean_y), xytext=(-2, -16),
                       textcoords="offset points", fontsize=7, ha="right",
                       color=_darken(col))

    phosphosites = {20: ("S20", 14), 21: ("S21", -24), 24: ("T24", 8), 51: ("S51", 8)}
    for r, (lab, dy) in phosphosites.items():
        ax_p.plot(r, plddt[r], marker="v", color="black", ms=4, zorder=5)
        ax_p.annotate(lab, xy=(r, plddt[r]), xytext=(0, dy), textcoords="offset points",
                       fontsize=6.5, ha="center")

    ax_p.set_xlabel(TXT["xlabel_a"])
    ax_p.set_ylabel(TXT["ylabel_a"])
    ax_p.set_ylim(30, 100)
    ax_p.set_xlim(1, 60)
    ax_p.legend(handles=[
        plt_line(col, lab) for (_, _, lab), col in zip(windows, colors)
    ], loc="lower left", fontsize=6.5, frameon=False)

    labels = TXT["bar_labels"]
    means = [rg, 11.5]
    errs = [0, 0.5]
    bars = ax_rg.bar(labels, means, yerr=errs, color=[fs.PALETTE_CATEGORICAL[0], "0.6"],
                      width=0.6, capsize=3)
    ax_rg.set_ylabel(TXT["ylabel_b"])
    # Ancrer au SOMMET DE LA MOUSTACHE (v + e), pas au sommet de la barre : sinon
    # la barre d'erreur traverse le label (« 11 5 » au lieu de « 11,5 » sur le
    # témoin globulaire, dont l'incertitude vaut 0,5 — /fig-check du 2026-08-17).
    for b, v, e in zip(bars, means, errs):
        ax_rg.annotate(fmt1(v), xy=(b.get_x() + b.get_width() / 2, v + e), xytext=(0, 3),
                        textcoords="offset points", ha="center", fontsize=7)
    ax_rg.set_ylim(0, max(means) * 1.35)

    fs.panel_labels(axes)
    out = save_600dpi(fig, FIGDIR / f"fig1_architecture{TXT['suffix']}")
    print("fig1:", [str(p) for p in out], "Rg=", rg)


def _darken(color, factor=0.6):
    """Assombrit une couleur matplotlib en gardant sa teinte (contraste sur fond teinté)."""
    import matplotlib.colors as mcolors
    r, g, b = mcolors.to_rgb(color)
    return (r * factor, g * factor, b * factor)


def plt_line(color, label):
    from matplotlib.lines import Line2D
    return Line2D([0], [0], color=color, lw=3, label=label, alpha=0.5)


def fig2_conservation():
    d = json.load(open(RESULTATS / "p3_3_msa_conservation.json"))
    pdp = sorted(d["per_domain_position"], key=lambda e: e["resnum"])

    resnums = [e["resnum"] for e in pdp]
    freqs = [100.0 * e["freq_majority"] for e in pdp]
    is_basic = [bool(e.get("is_basic_majority")) for e in pdp]
    aa = [e["h37rv_aa"] for e in pdp]

    fig, ax = fs.new_figure(PRESET, height_mm=95)

    perfect = [e["resnum"] for e in pdp if e["freq_majority"] == 1.0]
    for r in perfect:
        ax.axvspan(r - 0.5, r + 0.5, color=fs.PALETTE_CATEGORICAL[0], alpha=0.10, lw=0)

    basic_col = fs.PALETTE_CATEGORICAL[1]
    other_col = "0.35"
    for r, f, b in zip(resnums, freqs, is_basic):
        ax.bar(r, f, width=0.8, color=basic_col if b else other_col,
               edgecolor="none")

    module_end = 33
    ax.axvline(module_end + 0.5, color="black", lw=0.8, linestyle=":")
    ax.annotate(TXT["boundary_label"], xy=(module_end + 0.5, 60),
                xytext=(module_end + 2, 45), fontsize=6.5, ha="left")

    phosphosites = {20: ("S20", 18), 21: ("S21", 4), 24: ("T24", 11), 51: ("S51", 11)}
    for r, (lab, dy) in phosphosites.items():
        ax.plot(r, 106, marker="v", color="black", ms=4, clip_on=False)
        ax.annotate(lab, xy=(r, 106), xytext=(0, dy), textcoords="offset points",
                     fontsize=6, ha="center")

    # Résidus 2 à 5 consécutifs : un décalage vertical PAR RÉSIDU, sinon deux
    # labels adjacents partageant la même hauteur se touchent (R3/G4 rendus
    # « R3G4 », relevé par /fig-check le 2026-08-17).
    named = {2: ("G2", 0), 4: ("G4", 18), 8: ("A8", 0), 14: ("A14", 0),
             29: ("L29", 0), 32: ("E32", 0)}
    for r, (lab, dy) in named.items():
        ax.annotate(lab, xy=(r, 101), xytext=(0, dy), textcoords="offset points",
                     fontsize=6, ha="center", color=fs.PALETTE_CATEGORICAL[0])
    basic_labels = {3: ("R3", 9), 5: ("R5", 27), 18: ("K18", 0)}
    for r, (lab, dy) in basic_labels.items():
        ax.annotate(lab, xy=(r, 101), xytext=(0, dy), textcoords="offset points",
                     fontsize=6, ha="center", color=basic_col)

    ax.set_xlabel(TXT["xlabel2"])
    ax.set_ylabel(TXT["ylabel2"])
    ax.set_xlim(1, 61)
    ax.set_ylim(0, 128)
    ax.set_yticks([0, 25, 50, 75, 100])

    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor=other_col, label=TXT["legend2"][0]),
        Patch(facecolor=basic_col, label=TXT["legend2"][1]),
        Patch(facecolor=fs.PALETTE_CATEGORICAL[0], alpha=0.25, label=TXT["legend2"][2]),
    # `lower left` + `frameon=False` dessinait la légende À MÊME les barres des
    # résidus 2-20 : texte noir sur barres gris foncé, et pastilles indiscernables
    # des barres qu'elles expliquent (/fig-check du 2026-08-17). Aucune zone
    # INTÉRIEURE n'est libre : le coin supérieur droit paraissait vide mais porte
    # le label S51 (masqué au premier essai de correction). La légende sort donc
    # du cadre, au-dessus du panneau, sur une seule ligne.
    ], loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3,
       fontsize=6.5, frameon=False)

    out = save_600dpi(fig, FIGDIR / f"fig2_msa_conservation{TXT['suffix']}")
    print("fig2:", [str(p) for p in out])
    print("perfectly invariant positions:", perfect,
          [aa[resnums.index(r)] for r in perfect])


if __name__ == "__main__":
    fig1_architecture()
    fig2_conservation()
