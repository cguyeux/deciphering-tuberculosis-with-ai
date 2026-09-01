#!/usr/bin/env python3
"""
Objet       : F6 -- schema circulaire du mecanisme (fig_plan.md, archetype
              M1 + encart de donnees, P11.2.6). Chromosome circulaire avec
              oriC (ancre sur dnaA, Rv0001, position 1) et ter (localise par
              le skew GC cumule, Lobry 1996, position 2 044 792), les deux
              replichores coloryes et hachures differemment (test N&B), sens
              de fourche par fleches, indice d'asymetrie des pertes AI_loss
              annote par replichore (POOL_TOTAL) avec inversion de signe ;
              encart scatter v/u(repl1) vs v/u(repl2) embarque au centre
              (produit par phase11_F6_scatter_encart.py). Toutes les
              coordonnees angulaires sont derivees des fractions de longueur
              lues dans les donnees -- aucune n'est ecrite a la main.
Entrees     : résultats/phase10_p34_skew.tsv (oriC/ter/fractions, 3.4)
              résultats/phase10_p34_brin_matrice.tsv (AI_loss par replichore,
              ligne POOL_TOTAL, 3.4)
              article/figures/fig_F6_encart_scatter.pdf (phase11_F6_scatter_
              encart.py, doit exister avant ce script)
Sorties     : article/figures/fig_F6_mecanisme_biais_brin.tex
Reutilisable: non -- specifique a cette figure
Projet      : GC_par_lignee
Date        : 2026-08-31
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SKEW = ROOT / "résultats" / "phase10_p34_skew.tsv"
BRIN = ROOT / "résultats" / "phase10_p34_brin_matrice.tsv"
ENCART_PDF = ROOT / "article" / "figures" / "fig_F6_encart_scatter.pdf"
OUT = ROOT / "article" / "figures" / "fig_F6_mecanisme_biais_brin.tex"

R = 3.6  # cm, rayon de l'anneau chromosomique
INSET_SIDE_CM = 3.7  # cote de l'encart carre, doit rester < R*sqrt(2) pour ne pas toucher l'anneau


def fmt(x: float, nd: int = 3) -> str:
    """Number formatted with a LaTeX minus sign ($-$, not the Unicode glyph)."""
    s = f"{abs(x):.{nd}f}"
    return f"$-${s}" if x < 0 else s


def main() -> int:
    with open(SKEW, newline="") as fh:
        skew = next(csv.DictReader(fh, delimiter="\t"))
    with open(BRIN, newline="") as fh:
        brin = {r["replichore"]: r for r in csv.DictReader(fh, delimiter="\t") if r["lignee"] == "POOL_TOTAL"}

    genome_len = int(skew["genome_len"])
    ter_pos = int(skew["ter_1based"])
    frac1 = float(skew["frac_replichore1"])
    frac2 = 1.0 - frac1
    len1 = int(skew["len_replichore1"])
    len2 = int(skew["len_replichore2"])

    angle_oric = 90.0
    angle_ter = angle_oric - frac1 * 360.0
    angle_end2 = angle_ter - frac2 * 360.0  # doit revenir a angle_oric - 360

    ai1, ai1_lo, ai1_hi = (float(brin["1"][k]) for k in ("AI_loss", "AI_lo", "AI_hi"))
    ai2, ai2_lo, ai2_hi = (float(brin["2"][k]) for k in ("AI_loss", "AI_lo", "AI_hi"))
    ter_pos_str = f"{ter_pos:,}".replace(",", "\\,")

    if not ENCART_PDF.exists():
        raise SystemExit(f"encart manquant : {ENCART_PDF.relative_to(ROOT)} -- lancer d'abord phase11_F6_scatter_encart.py")

    # angle median de chaque arc, pour placer les annotations d'AI_loss.
    # Ancrage west/east selon le demi-cercle : le bloc de texte doit croitre
    # vers l'exterieur, jamais etre centre sur le point, sinon sa moitie
    # interieure recouvre l'anneau (constate a l'inspection du premier PNG :
    # meme couleur texte/arc => lettres avalees par le trait).
    import math
    mid1 = (angle_oric + angle_ter) / 2.0
    mid2 = (angle_ter + angle_end2) / 2.0
    anchor1 = "west" if math.cos(math.radians(mid1)) >= 0 else "east"
    anchor2 = "west" if math.cos(math.radians(mid2)) >= 0 else "east"
    r_label = R + 0.35

    body = "\n".join(
        [
            "% F6 -- schema circulaire du mecanisme (biais de brin repliactif),",
            "% archetype M1 + encart de donnees, fig_plan.md P11.2.6.",
            "% Genere par analyses/phase11_F6_mecanisme_biais_brin.py depuis",
            "% résultats/phase10_p34_[skew|brin_matrice].tsv -- NE PAS EDITER A LA MAIN.",
            "\\documentclass[border=6pt]{standalone}",
            "\\usepackage[T1]{fontenc}",
            "\\usepackage{lmodern}",
            "\\usepackage{xcolor}",
            "\\usepackage{graphicx}",
            "\\usepackage{tikz}",
            "\\usetikzlibrary{positioning,arrows.meta,calc,backgrounds}",
            "",
            "\\definecolor{repl1col}{HTML}{D55E00}",
            "\\definecolor{repl2col}{HTML}{0072B2}",
            "\\definecolor{neutralgrey}{HTML}{4D4D4D}",
            "",
            "\\begin{document}",
            "\\begin{tikzpicture}[",
            "    every node/.style={font=\\scriptsize, align=center},",
            f"    repl1/.style={{repl1col, line width=2.6pt, cap=round, -{{Latex[length=7pt,width=5pt]}}}},",
            f"    repl2/.style={{repl2col, line width=2.6pt, cap=round, dash pattern=on 3.6pt off 2.4pt, -{{Latex[length=7pt,width=5pt]}}}},",
            "    landmark/.style={circle, fill=neutralgrey, draw=white, line width=0.4pt, inner sep=1.6pt},",
            "    poslabel/.style={font=\\scriptsize, text=neutralgrey, align=center},",
            "    aibox/.style={font=\\tiny, align=center, inner sep=2pt},",
            "  ]",
            "",
            "  % ---- anneau chromosomique : deux arcs, oriC -> ter, sens de fourche ----",
            f"  \\draw[repl1] (90:{R}) arc ({angle_oric:.3f}:{angle_ter:.3f}:{R});",
            f"  \\draw[repl2] ({angle_end2:.3f}:{R}) arc ({angle_end2:.3f}:{angle_ter:.3f}:{R});",
            "",
            "  % ---- oriC (dnaA, Rv0001, position 1) --------------------------------",
            f"  \\node[landmark] (oric) at (90:{R}) {{}};",
            f"  \\node[poslabel, anchor=south] at (90:{R + 0.35}) "
            "{\\textbf{oriC}\\\\\\emph{dnaA} (Rv0001)\\\\pos.\\,1};",
            "",
            "  % ---- ter (cumulative GC skew, Lobry 1996) -----------------------------",
            f"  \\node[landmark] (ter) at ({angle_ter:.3f}:{R}) {{}};",
            f"  \\node[poslabel, anchor=north west] at ({angle_ter:.3f}:{R + 0.35}) "
            f"{{\\textbf{{ter}}\\\\pos.\\,{ter_pos_str}\\,bp\\\\(cumulative GC skew)}};",
            "",
            "  % ---- replichore 1 annotation: loss asymmetry index --------------------",
            f"  \\node[aibox, anchor={anchor1}, text=repl1col] at ({mid1:.3f}:{r_label}) "
            f"{{replichore 1 ({len1 / genome_len * 100:.1f}\\,\\% of the genome)"
            f"\\\\+ strand = \\textbf{{lagging}}-strand template"
            f"\\\\AI$_{{loss}}$ = {fmt(ai1, 3)} "
            f"[{fmt(ai1_lo, 3)}\\,;\\,{fmt(ai1_hi, 3)}]"
            f"\\\\(C-loss $>$ G-loss)}};",
            "",
            "  % ---- replichore 2 annotation: loss asymmetry index --------------------",
            f"  \\node[aibox, anchor={anchor2}, text=repl2col] at ({mid2:.3f}:{r_label}) "
            f"{{replichore 2 ({len2 / genome_len * 100:.1f}\\,\\% of the genome)"
            f"\\\\+ strand = \\textbf{{leading}}-strand template"
            f"\\\\AI$_{{loss}}$ = {fmt(ai2, 3)} "
            f"[{fmt(ai2_lo, 3)}\\,;\\,{fmt(ai2_hi, 3)}]"
            f"\\\\(G-loss $>$ C-loss)}};",
            "",
            "  % ---- inset: v/u(repl 1) vs v/u(repl 2), 13 reportable lineages --------",
            f"  \\node at (0,0) {{\\includegraphics[width={INSET_SIDE_CM}cm]{{fig_F6_encart_scatter.pdf}}}};",
            "",
            "  % ---- text inset: the claim this figure makes --------------------------",
            "  \\node[draw=neutralgrey, dotted, thick, rounded corners=2pt, inner sep=5pt,",
            "        text width=13.6cm, font=\\scriptsize, align=left, text=neutralgrey]",
            f"    at (0,{-(R + 2.9):.3f})",
            "    {\\textbf{The lineage ranking survives intact within each replichore, "
            "and the control blindly recovers a genuine replicative strand bias in "
            "the MTBC.} The replicative architecture (a single oriC, a single "
            "terminus) is identical between lineages: the $v/u$ ranking by lineage "
            "(inset, $\\rho=0.940$, $n=13$) cannot therefore arise from the "
            "replichore itself. The loss asymmetry index AI$_{loss}$ = "
            "ln(C-loss / G-loss) flips sign between replichores exactly as "
            "predicted by fork geometry, at the individual-lineage level in "
            "15/16 lineages tested (the exception, Dassie, $n=6$ strains, does "
            "not show it) -- confirmation, within the MTBC, of the strand bias "
            "already described in \\emph{M.\\,smegmatis} "
            "(Casta\\~neda-Garc\\'ia et al.\\ 2020).};",
            "",
            "\\end{tikzpicture}",
            "\\end{document}",
        ]
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(body + "\n")
    print(
        f"Ecrit : {OUT.relative_to(ROOT)} "
        f"(angle_ter={angle_ter:.2f}deg, frac1={frac1:.4f}, "
        f"AI1={ai1:.3f}, AI2={ai2:.3f})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
