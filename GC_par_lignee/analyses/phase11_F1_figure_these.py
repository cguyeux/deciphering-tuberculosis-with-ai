#!/usr/bin/env python3
"""
Objet       : F1 -- figure-these (archetype C2) : contraste entre le GC
              observe (identique a 70 ppm pres entre lignees, un point
              unique) et le GC d'equilibre mutationnel GC_eq = u/(u+v) par
              lignee (32-55 %, un nuage de 17 points), relies par une fleche
              dont la position d'arrivee restitue l'ecart reel (A15, A18).
              Genere les coordonnees depuis les donnees, ecrit le .tex TikZ
              standalone complet (fig_plan.md F1, P11.2.1).
Entrees     : résultats/phase5_p53_force_requise.tsv (panel final 17 lignees
              / 519 souches, P5.3 -- meme panel que F4/F5 pour coherence
              inter-figures)
Sorties     : article/figures/fig_F1_figure_these.tex
Reutilisable: non -- specifique a ce jeu de resultats et a cette figure
Projet      : GC_par_lignee
Date        : 2026-08-31
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "résultats" / "phase5_p53_force_requise.tsv"
OUT = ROOT / "article" / "figures" / "fig_F1_figure_these.tex"

GC_OBSERVE_H37RV = 65.6147  # %, identique (± 70 ppm) pour toutes les lignees, A1

AXIS_X0_PCT = 0.0
AXIS_X1_PCT = 100.0
AXIS_LEN_CM = 14.0
AXIS_Y = 0.0
ARROW_Y0_CM = 0.6  # marge entre l'axe et la premiere rangee de fleches
DY_CM = 0.34  # ecart vertical entre deux fleches consecutives du nuage

LANDMARKS = {"L9", "Microti", "La4"}
# Microti (rang 15) et La4 (rang 16) sont deux rangees adjacentes en haut du
# nuage (0,34 cm d'ecart) : leurs etiquettes deux-lignes se chevauchent sans
# ce nudge vertical manuel, constate a l'inspection du PNG.
LABEL_DY_CM = {"Microti": -0.22, "La4": 0.22}

# libelles d'affichage (les codes de la table sont ceux du panel P5.3)
DISPLAY_NAME = {
    "Caprae_La2": "\\emph{M.\\,caprae} (La2)",
    "Orygis_La3": "\\emph{M.\\,orygis} (La3)",
}


def x_of(gc_pct: float) -> float:
    return AXIS_LEN_CM * (gc_pct - AXIS_X0_PCT) / (AXIS_X1_PCT - AXIS_X0_PCT)


def main() -> int:
    with open(SRC, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    rows = [r for r in rows if r["lignee"]]
    rows.sort(key=lambda r: float(r["gc_eq"]))
    n = len(rows)

    x_obs = x_of(GC_OBSERVE_H37RV)

    arrow_lines = []
    dot_lines = []
    label_lines = []
    for rank, r in enumerate(rows):
        lignee = r["lignee"]
        gc_eq_pct = float(r["gc_eq"]) * 100.0
        n_souches = r["n_souches"]
        sous_puissante = r["sous_puissante"] == "True"
        x_eq = x_of(gc_eq_pct)
        y = ARROW_Y0_CM + rank * DY_CM

        style = "pullweak" if sous_puissante else "pull"
        arrow_lines.append(
            f"  \\draw[{style}] ({x_obs:.3f},{y:.3f}) -- ({x_eq:.3f},{y:.3f});"
        )
        dot_lines.append(
            f"  \\node[eqpoint] at ({x_eq:.3f},{y:.3f}) {{}};"
        )

        if lignee in LANDMARKS:
            disp = DISPLAY_NAME.get(lignee, lignee)
            gc_str = f"{gc_eq_pct:.1f}\\,\\%"
            # Chaque fleche va de x_eq (le point) vers x_obs, toujours vers la
            # droite (tous les GC_eq sont < GC observe) : la zone a GAUCHE du
            # point est donc toujours libre sur cette rangee, quel que soit le
            # lignee -- ancrage a l'est systematique, jamais de recouvrement.
            y_label = y + LABEL_DY_CM.get(lignee, 0.0)
            label_lines.append(
                f"  \\node[landmark, anchor=east] at "
                f"({x_eq - 0.28:.3f},{y_label:.3f}) "
                f"{{\\textbf{{{disp}}}\\\\{gc_str} (n={n_souches})}};"
            )

    y_top = ARROW_Y0_CM + (n - 1) * DY_CM

    ticks = []
    for pct in (0, 20, 40, 60, 80, 100):
        xt = x_of(pct)
        ticks.append(
            f"  \\draw[axisline] ({xt:.3f},{AXIS_Y - 0.55:.3f}) -- "
            f"({xt:.3f},{AXIS_Y - 0.45:.3f});"
        )
        ticks.append(
            f"  \\node[axistick] at ({xt:.3f},{AXIS_Y - 0.75:.3f}) {{{pct}}};"
        )

    body = "\n".join(
        [
            "% F1 -- figure-these (archetype C2) : GC observe (point unique, identique",
            "% a 70 ppm pres) contre GC d'equilibre mutationnel GC_eq = u/(u+v) par",
            "% lignee (nuage, 32-55 %), fig_plan.md P11.2.1.",
            "% Genere par analyses/phase11_F1_figure_these.py depuis",
            "% résultats/phase5_p53_force_requise.tsv -- NE PAS EDITER A LA MAIN.",
            "\\documentclass[border=6pt]{standalone}",
            "\\usepackage[T1]{fontenc}",
            "\\usepackage{lmodern}",
            "\\usepackage{xcolor}",
            "\\usepackage{tikz}",
            "\\usetikzlibrary{positioning,arrows.meta,calc,backgrounds}",
            "",
            "\\definecolor{keepblue}{HTML}{0072B2}",
            "\\definecolor{dropred}{HTML}{D55E00}",
            "\\definecolor{neutralgrey}{HTML}{4D4D4D}",
            "\\definecolor{obsdot}{HTML}{000000}",
            "",
            "\\begin{document}",
            "\\begin{tikzpicture}[",
            "    every node/.style={font=\\scriptsize, align=center},",
            "    axisline/.style={line width=0.8pt, color=neutralgrey},",
            "    axistick/.style={font=\\tiny, text=neutralgrey},",
            "    pull/.style={-{Latex[length=3.2pt,width=2.2pt]}, line width=0.55pt,",
            "      color=keepblue, opacity=0.75, shorten <=1.5pt, shorten >=0.5pt},",
            "    pullweak/.style={-{Latex[length=3.2pt,width=2.2pt]}, line width=0.55pt,",
            "      color=keepblue, opacity=0.35, dashed, shorten <=1.5pt, shorten >=0.5pt},",
            "    eqpoint/.style={circle, fill=keepblue, draw=white, line width=0.3pt,",
            "      inner sep=0.9pt},",
            "    landmark/.style={font=\\tiny, text=neutralgrey},",
            "  ]",
            "",
            "  % ---- axe GC, 0 a 100 % -------------------------------------------",
            f"  \\draw[axisline] ({x_of(0):.3f},{AXIS_Y:.3f}) -- "
            f"({x_of(100):.3f},{AXIS_Y:.3f});",
            *ticks,
            f"  \\node[axistick] at ({x_of(50):.3f},{AXIS_Y - 1.15:.3f}) "
            "{\\normalsize GC (\\%), whole genome};",
            "",
            "  % ---- single point: observed GC, identical for the 17 lineages ----",
            f"  \\draw[dropred, thick, dashed] ({x_obs:.3f},{AXIS_Y:.3f}) -- "
            f"({x_obs:.3f},{y_top + 0.55:.3f});",
            f"  \\node[circle, fill=dropred, draw=white, line width=0.4pt, "
            f"inner sep=1.6pt] (obs) at ({x_obs:.3f},{y_top + 0.55:.3f}) {{}};",
            f"  \\node[font=\\scriptsize, text=dropred, anchor=south, align=center] "
            f"at ({x_obs:.3f},{y_top + 0.7:.3f}) "
            "{\\textbf{Observed GC (H37Rv)}\\\\65.6147\\,\\% -- a single point,\\\\"
            "max.\\ spread between lineages~: 70~ppm};",
            "",
            "  % ---- cloud: mutational equilibrium GC by lineage ------------------",
            "  % arrow = pull of the mutational flux (GC:AT loss / AT:GC gain,",
            "  % A15) from the observed GC toward the lineage's own GC_eq;",
            "  % dashed + reduced opacity = underpowered lineage (n<10, A23).",
            *arrow_lines,
            *dot_lines,
            *label_lines,
            "",
            "  % ---- inset: the claim this figure makes ----------------------------",
            "  \\node[draw=neutralgrey, dotted, thick, rounded corners=2pt, "
            "inner sep=5pt,",
            "        text width=12.4cm, font=\\scriptsize, align=left, "
            "text=neutralgrey]",
            f"    at ({x_of(50):.3f},{AXIS_Y - 2.6:.3f})",
            "    {\\textbf{What distinguishes MTBC lineages is not GC content "
            "itself, but its mutational target.} The 17 MTBC lineages share an "
            "observed GC content identical to within 70~ppm (single red point); "
            "their polarized mutational flux (GC:AT losses versus AT:GC gains, "
            "A15) pulls them toward equilibria GC$_{eq}$ = u/(u+v) spanning 23 "
            "percentage points (32.1~\\% to 55.4~\\%, blue cloud). Arrow length = "
            "true spread on the GC axis, not a scale artifact: both quantities "
            "are plotted on the same axis.};",
            "",
            "\\end{tikzpicture}",
            "\\end{document}",
        ]
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(body + "\n")
    print(f"Ecrit : {OUT.relative_to(ROOT)} ({n} lignees, {sum(1 for r in rows if r['sous_puissante']=='True')} sous-puissantes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
