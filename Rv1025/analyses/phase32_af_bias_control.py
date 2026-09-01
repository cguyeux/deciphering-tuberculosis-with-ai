#!/usr/bin/env python3
"""phase32_af_bias_control.py — P5.2 : le site de Rv1025 est-il VRAIMENT plus ouvert, ou est-ce AlphaFold ?

phase31 a mesuré que le métal de Rv1025 est plus accessible (19,8 % de sphère ouverte,
14,3 % débouchant sur le solvant) que trois zincs catalytiques de référence mesurés sur
CRISTAL (10-14 % et 0,5-4,9 %), et infiniment plus qu'un zinc structural (0,0 %).

BIAIS ÉVIDENT, et fatal s'il n'est pas mesuré : Rv1025 est un MODÈLE AlphaFold, les
contrôles sont des structures CRISTALLOGRAPHIQUES. Un modèle prédit peut être moins
compacté qu'un cristal à 2 A, ce qui gonflerait mécaniquement l'espace libre. Comparer
les deux sans quantifier ce biais serait une faute de méthode.

CONTRÔLE : refaire la mesure sur le modèle AlphaFold DB des MÊMES protéines de contrôle.
La position du métal, absente des modèles AFDB (qui sont apo), est importée du cristal
après superposition Kabsch sur les carbones alpha appariés par numéro de résidu.
  - Si AF-anhydrase montre ~20 % comme Rv1025, alors 19,8 % est la LIGNE DE BASE d'un
    modèle AlphaFold et la comparaison de phase31 ne vaut rien.
  - Si AF-anhydrase reste proche de son cristal, le biais est petit et l'ouverture de
    Rv1025 est une propriété du repli, pas de la méthode.

SECOND CONTRÔLE, indépendant : la dispersion entre les 5 modèles AF3 de Rv1025 et entre
les trois métaux. Un chiffre qui varie fortement d'un modèle à l'autre ne se rapporte pas.

Lit   : data/reference_metalloproteins/*.pdb, résultats/af3_metal_out/
Écrit : résultats/metal_accessibility_af_bias.tsv, .txt
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fast_blas import ensure  # noqa: E402  # type: ignore[import-not-found]

ensure()

from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

from metal_geometry import (kabsch, measure, parse_cif,  # noqa: E402  # type: ignore[import-not-found]
                            parse_pdb, pick_site)

ROOT = Path(__file__).resolve().parent.parent
PDBDIR = ROOT / "data/reference_metalloproteins"
AF3 = ROOT / "résultats/af3_metal_out"
OUT = ROOT / "résultats"

# (PDB cristal, UniProt du modèle AFDB, nom)
PAIRS = [("1CA2", "P00918", "anhydrase carbonique II"),
         ("8TLN", "P00800", "thermolysine")]

out = []


def say(s=""):
    print(s)
    out.append(s)


def ca_by_resnum(poly):
    """{numéro de résidu : (coord, nom du résidu)} pour les carbones alpha."""
    d = {}
    for a in poly:
        if len(a) >= 7 and a[4] == "CA":
            try:
                d[int(a[6])] = (np.array(a[1:4], dtype=float), a[5])
            except ValueError:
                continue
    return d


rows = [("structure", "origine", "ouvert_2.4A", "ouvert_2.6A", "ouvert_2.8A",
         "debouche_solvant", "rmsd_superposition", "n_CA_apparies")]

say("CONTRÔLE DE BIAIS — même protéine, cristal contre modèle AlphaFold DB")
say("(métal importé du cristal après superposition Kabsch sur les CA appariés)")
say()
say(f"{'protéine':<28} {'origine':<10} {'ouvert 2,4 A':>13} {'débouché':>10} {'RMSD':>7}")
for pid, upid, nom in PAIRS:
    xtal = PDBDIR / f"{pid}.pdb"
    afp = PDBDIR / f"AF-{upid}.pdb"
    # Un téléchargement raté laisse un fichier (corps du 404) : re-tenter si trop petit,
    # sinon l'échec se fige définitivement dans le cache.
    if not xtal.exists() or xtal.stat().st_size < 5000:
        subprocess.run(["curl", "-sS", "--max-time", "120", "-o", str(xtal),
                        f"https://files.rcsb.org/download/{pid}.pdb"], check=False)
    if not afp.exists() or afp.stat().st_size < 5000:
        # NE PAS coder la version du modèle en dur (v4 -> v6 en 2026) : la demander à l'API.
        r = subprocess.run(["curl", "-sS", "--max-time", "60",
                            f"https://alphafold.ebi.ac.uk/api/prediction/{upid}"],
                           capture_output=True, text=True)
        url = ""
        try:
            import json
            url = json.loads(r.stdout)[0].get("pdbUrl", "")
        except Exception:
            url = ""
        if url:
            subprocess.run(["curl", "-sS", "--max-time", "180", "-o", str(afp), url],
                           check=False)
    if not afp.exists() or afp.stat().st_size < 5000:
        say(f"  {nom} : modèle AFDB indisponible — contrôle impossible")
        continue

    px, mx = parse_pdb(xtal)
    site = pick_site(px, mx)
    assert site is not None, f"aucun métal dans {pid}"
    _, centre_x, _ = site
    op_x, esc_x = measure(px, centre_x)

    pa, _ = parse_pdb(afp)
    cx, ca = ca_by_resnum(px), ca_by_resnum(pa)
    # La numérotation d'un cristal et celle d'UniProt diffèrent souvent d'un DÉCALAGE
    # constant (propeptide clivé : 8TLN numérote la protéine mature, P00800 inclut le
    # pro-domaine). Chercher le décalage qui maximise les résidus IDENTIQUES appariés,
    # au lieu de supposer une numérotation commune.
    best_off, common = 0, []
    for off in range(-400, 401):
        c = [r for r in cx if (r + off) in ca and cx[r][1] == ca[r + off][1]]
        if len(c) > len(common):
            best_off, common = off, c
    if best_off:
        say(f"    (décalage de numérotation détecté : {best_off:+d} résidus)")
    ca = {r: ca[r + best_off] for r in common}
    if len(common) < 50:
        say(f"  {nom} : appariement insuffisant ({len(common)} CA) — contrôle abandonné")
        continue
    M = np.array([ca[r][0] for r in common])
    T = np.array([cx[r][0] for r in common])
    R, cm, ct, rmsd = kabsch(M, T)
    moved = [(a[0], *(((np.array(a[1:4], dtype=float) - cm) @ R.T) + ct)) for a in pa]
    op_a, esc_a = measure(moved, centre_x)

    say(f"  {nom:<26} {'cristal':<10} {op_x[2.4]:>12.1%} {esc_x:>9.1%} {'-':>7}")
    say(f"  {'':<26} {'AlphaFold':<10} {op_a[2.4]:>12.1%} {esc_a:>9.1%} {rmsd:>6.2f} A"
        f"   ({len(common)} CA appariés)")
    rows.append((nom, "cristal", f"{op_x[2.4]:.3f}", f"{op_x[2.6]:.3f}", f"{op_x[2.8]:.3f}",
                 f"{esc_x:.3f}", "-", "-"))
    rows.append((nom, "AlphaFold", f"{op_a[2.4]:.3f}", f"{op_a[2.6]:.3f}", f"{op_a[2.8]:.3f}",
                 f"{esc_a:.3f}", f"{rmsd:.2f}", str(len(common))))

# ------------------------------------------------ dispersion entre modèles Rv1025
say()
say("DISPERSION DE Rv1025 — 5 modèles AF3 par métal")
say(f"{'jeu':<12} {'ouvert 2,4 A (moy ± et)':>28} {'débouché solvant':>22}")
for metal in ("fe", "zn", "mn"):
    vals, escs = [], []
    for k in range(5):
        f = AF3 / f"fold_rv1025_holo_{metal}/fold_rv1025_holo_{metal}_model_{k}.cif"
        if not f.exists():
            continue
        poly, metals = parse_cif(f)
        site = pick_site(poly, metals)
        if site is None:
            continue
        op, esc = measure(poly, site[1])
        vals.append(op[2.4])
        escs.append(esc)
    if vals:
        say(f"  {metal.upper():<10} {np.mean(vals):>15.1%} ± {np.std(vals):<8.1%} "
            f"{np.mean(escs):>13.1%} ± {np.std(escs):.1%}   (n={len(vals)})")
        rows.append((f"Rv1025 holo {metal.upper()} (5 modèles)", "AlphaFold3",
                     f"{np.mean(vals):.3f}", "-", "-", f"{np.mean(escs):.3f}",
                     f"et={np.std(vals):.3f}", str(len(vals))))

(OUT / "metal_accessibility_af_bias.tsv").write_text(
    "\n".join("\t".join(r) for r in rows) + "\n")
(OUT / "metal_accessibility_af_bias.txt").write_text("\n".join(out) + "\n")
say()
say("Écrit : résultats/metal_accessibility_af_bias.tsv et .txt")
