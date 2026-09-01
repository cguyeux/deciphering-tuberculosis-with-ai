#!/usr/bin/env python3
"""phase31_metal_accessibility.py — P5.2 REFORMULÉE : le métal est-il chimiquement ACCESSIBLE ?

POURQUOI PAS DU DOCKING (challenge appliqué avant de coder). La piste demandait un docking
de fragment chélateur sur le site métal. Trois objections dirimantes :
  1. le site est PRÉDIT et le métal physiologique INCONNU (Zn/Fe/Mn tous accommodés) —
     un score de docking serait une prédiction de troisième ordre ;
  2. les fonctions de score usuelles (Vina et dérivées) ne modélisent PAS la coordination
     métallique ; elles placent volontiers n'importe quoi dans une cavité ;
  3. MODÈLE NUL écrasant : une cavité de 255 A^3 accepte presque tout fragment, et les
     écarts de score entre chimiotypes sont inférieurs à l'erreur de la méthode.
Un chiffre de docking serait donc non falsifiable — exactement ce que ce projet s'interdit.

CE QUI EST MESURABLE, et qui répond à la question SOUS-JACENTE (« ce métal est-il attaquable
par un ligand exogène ? ») : la GÉOMÉTRIE de la sphère de coordination. Un métal catalytique
garde des positions de coordination libres et ouvertes sur le solvant, puisqu'il doit lier
son substrat ; un métal purement STRUCTURAL est saturé par la protéine et enfoui. C'est une
mesure d'espace, pas un score empirique, et surtout elle se CALIBRE.

CALIBRATION (le coeur du dispositif) : la même mesure est appliquée à
  - trois métalloenzymes à zinc CATALYTIQUE, toutes cibles de médicaments coordinants :
    anhydrase carbonique II (1CA2, sulfonamides), thermolysine (8TLN, hydroxamates),
    LpxC (4MDT, hydroxamates) ;
  - un zinc STRUCTURAL saturé : doigt de zinc (1ZNF).
Les eaux et ligands sont EXCLUS partout : on mesure l'espace disponible POUR un ligand,
donc on ne compte pas celui qui l'occupe déjà. Rv1025 doit se comparer aux premières et se
distinguer du second, sinon l'argument de druggabilité du site tombe.

Lit   : résultats/af3_metal_out/fold_rv1025_holo_fe/..._model_0.cif + PDB de contrôle
Écrit : résultats/metal_accessibility.tsv, .txt
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fast_blas import ensure  # noqa: E402  # type: ignore[import-not-found]

ensure()

from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CIF = ROOT / "résultats/af3_metal_out/fold_rv1025_holo_fe/fold_rv1025_holo_fe_model_0.cif"
PDBDIR = ROOT / "data/reference_metalloproteins"
OUT = ROOT / "résultats"

CONTROLS = [("1CA2", "anhydrase carbonique II", "catalytique (sulfonamides)"),
            ("8TLN", "thermolysine", "catalytique (hydroxamates)"),
            ("4MDT", "LpxC", "catalytique (hydroxamates)"),
            ("1ZNF", "doigt de zinc", "STRUCTURAL (contrôle négatif)")]
METALS = {"ZN", "FE", "MN", "CO", "NI", "CU", "MG"}
DONOR = {"N", "O", "S"}
COORD_R = 2.1          # rayon où l'on sonde les positions de coordination (A)
THRESHOLDS = (2.4, 2.6, 2.8)
N_DIR = 4000
ESCAPE_R = 10.0        # distance à atteindre pour considérer qu'on rejoint le solvant
ESCAPE_CLEAR = 1.9     # dégagement minimal le long du trajet de fuite (A)

out = []


def say(s=""):
    print(s)
    out.append(s)


def fib_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)]


def parse_pdb(path):
    """Renvoie (atomes polymère, métaux) ; MODEL 1 seulement, eaux et ligands exclus."""
    poly, metals = [], []
    for line in path.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM") and line[17:20].strip() != "HOH":
            el = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            poly.append((el, float(line[30:38]), float(line[38:46]), float(line[46:54])))
        elif line.startswith("HETATM"):
            resn = line[17:20].strip().upper()
            if resn in METALS:
                metals.append((resn, float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return poly, metals


def parse_cif(path):
    """Parseur minimal d'atom_site mmCIF (sortie AlphaFold3)."""
    lines = path.read_text().splitlines()
    cols, start = [], None
    for i, line in enumerate(lines):
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".")[1])
            start = i + 1
        elif cols and start is not None and i >= start:
            break
    idx = {c: k for k, c in enumerate(cols)}
    poly, metals = [], []
    for line in lines[start:]:
        if line.startswith("#") or not line.strip():
            break
        f = line.split()
        if len(f) < len(cols):
            continue
        comp = f[idx["label_comp_id"]].upper()
        el = f[idx["type_symbol"]].upper()
        xyz = (float(f[idx["Cartn_x"]]), float(f[idx["Cartn_y"]]), float(f[idx["Cartn_z"]]))
        if comp == "HOH":
            continue
        if comp in METALS or el in METALS:
            metals.append((comp, *xyz))
        elif f[idx["group_PDB"]] == "ATOM":
            poly.append((el, *xyz))
    return poly, metals


def analyse(name, poly, metals, label=""):
    if not metals:
        say(f"  {name} : AUCUN métal trouvé — ignoré")
        return None
    P = np.array([[a[1], a[2], a[3]] for a in poly])
    El = [a[0] for a in poly]
    # choisir le métal le mieux coordonné par la protéine (= le site, pas un ion de surface)
    best = None
    for m in metals:
        c = np.array(m[1:])
        d = np.linalg.norm(P - c, axis=1)
        nd = sum(1 for k in np.where(d < 2.8)[0] if El[k] in DONOR)
        if best is None or nd > best[0]:
            best = (nd, c, m[0])
    assert best is not None
    ndon, centre, sym = best

    dirs = fib_sphere(N_DIR)
    probes = centre + dirs * COORD_R
    # distance de chaque sonde à l'atome polymère le plus proche
    near = np.array([np.min(np.linalg.norm(P - p, axis=1)) for p in probes])
    openfrac = {t: float((near >= t).mean()) for t in THRESHOLDS}

    # fuite vers le solvant : le long des directions ouvertes, marcher jusqu'à ESCAPE_R
    steps = np.arange(COORD_R, ESCAPE_R, 0.5)
    esc = 0
    opened = np.where(near >= THRESHOLDS[1])[0]
    for k in opened:
        path_pts = centre + dirs[k] * steps[:, None]
        if all(np.min(np.linalg.norm(P - q, axis=1)) >= ESCAPE_CLEAR for q in path_pts):
            esc += 1
    escfrac = esc / N_DIR

    say(f"  {name:<22} {sym:<3} donneurs {ndon}  "
        f"ouvert {openfrac[2.4]:5.1%}/{openfrac[2.6]:5.1%}/{openfrac[2.8]:5.1%}  "
        f"débouché solvant {escfrac:5.1%}   {label}")
    res: dict = dict(name=name, metal=sym, donors=ndon, esc=escfrac, kind=label or "candidat")
    res.update({f"open{t}": openfrac[t] for t in THRESHOLDS})
    return res


say("Accessibilité de la sphère de coordination — eaux et ligands EXCLUS partout")
say(f"sonde à {COORD_R} A du métal, {N_DIR} directions ; « ouvert » = aucun atome de la")
say(f"protéine à moins de 2,4 / 2,6 / 2,8 A de la sonde (trois seuils = test de robustesse)")
say(f"« débouché solvant » = direction ouverte le long de laquelle on atteint {ESCAPE_R:.0f} A")
say(f"sans jamais passer à moins de {ESCAPE_CLEAR} A d'un atome.")
say()

PDBDIR.mkdir(parents=True, exist_ok=True)
rows = []
poly, metals = parse_cif(CIF)
say(f"Rv1025 : {len(poly)} atomes polymère, {len(metals)} métal(aux)")
r = analyse("Rv1025 (holo Fe)", poly, metals, "<- LE CANDIDAT")
if r:
    rows.append(r)
for pid, nom, kind in CONTROLS:
    f = PDBDIR / f"{pid}.pdb"
    if not f.exists():
        subprocess.run(["curl", "-sS", "--max-time", "120", "-o", str(f),
                        f"https://files.rcsb.org/download/{pid}.pdb"], check=False)
    if not f.exists() or f.stat().st_size < 1000:
        say(f"  {pid} : téléchargement échoué — ignoré")
        continue
    poly, metals = parse_pdb(f)
    r = analyse(f"{pid} {nom}", poly, metals, kind)
    if r:
        rows.append(r)

say()
say("LECTURE : un métal CATALYTIQUE garde des positions de coordination libres ET débouchant")
say("sur le solvant, puisqu'il doit lier son substrat. Un métal STRUCTURAL est saturé et")
say("enfoui. Si Rv1025 se range avec les catalytiques, son site est chimiquement attaquable")
say("par un ligand coordinant — ce qui est l'argument de druggabilité, SANS prétendre prédire")
say("une affinité. Si Rv1025 ressemble au doigt de zinc, l'argument tombe.")

hdr = ["structure", "metal", "donneurs_proteine", "ouvert_2.4A", "ouvert_2.6A",
       "ouvert_2.8A", "debouche_solvant", "type"]
tab = ["\t".join(hdr)]
for r in rows:
    tab.append("\t".join([r["name"], r["metal"], str(r["donors"]),
                          f"{r['open2.4']:.3f}", f"{r['open2.6']:.3f}", f"{r['open2.8']:.3f}",
                          f"{r['esc']:.3f}", r.get("kind", "candidat")]))
(OUT / "metal_accessibility.tsv").write_text("\n".join(tab) + "\n")
(OUT / "metal_accessibility.txt").write_text("\n".join(out) + "\n")
say()
say("Écrit : résultats/metal_accessibility.tsv et .txt")
