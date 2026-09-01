"""metal_geometry.py — mesure de l'accessibilité d'un site métallique.

Extrait de phase31 pour être partagé avec phase32 (contrôle de biais AlphaFold vs
cristal) : deux scripts qui compareraient des sites avec deux implémentations
légèrement différentes ne compareraient rien du tout.

Principe : sonder la sphère de coordination autour d'un ion métallique. Une direction
est OUVERTE si une sonde placée à distance de coordination n'y rencontre aucun atome
de la protéine ; elle DÉBOUCHE si l'on peut ensuite s'éloigner jusqu'au solvant sans
jamais frôler la protéine. Eaux et ligands sont toujours exclus : on mesure la place
disponible POUR un ligand, pas celle qu'occupe déjà celui qui est là.
"""
import numpy as np

METALS = {"ZN", "FE", "MN", "CO", "NI", "CU", "MG"}
DONOR = {"N", "O", "S"}
COORD_R = 2.1          # rayon de sondage des positions de coordination (A)
THRESHOLDS = (2.4, 2.6, 2.8)
N_DIR = 4000
ESCAPE_R = 10.0        # distance à atteindre pour considérer qu'on rejoint le solvant
ESCAPE_CLEAR = 1.9     # dégagement minimal le long du trajet de fuite (A)


def fib_sphere(n=N_DIR):
    """n directions quasi uniformes sur la sphère (spirale de Fibonacci)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)]


def parse_pdb(path, model=1):
    """(atomes polymère, métaux) d'un PDB ; un seul modèle, eaux exclues."""
    poly, metals = [], []
    seen = 0
    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            seen += 1
            if seen > model:
                break
        if line.startswith("ENDMDL") and seen >= model:
            break
        if line.startswith("ATOM") and line[17:20].strip() != "HOH":
            el = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            poly.append((el, float(line[30:38]), float(line[38:46]), float(line[46:54]),
                         line[12:16].strip(), line[17:20].strip(), line[22:26].strip()))
        elif line.startswith("HETATM") and line[17:20].strip().upper() in METALS:
            metals.append((line[17:20].strip().upper(),
                           float(line[30:38]), float(line[38:46]), float(line[46:54])))
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
            poly.append((el, *xyz, f[idx["label_atom_id"]], comp, f[idx["label_seq_id"]]))
    return poly, metals


def pick_site(poly, metals):
    """Le métal le MIEUX coordonné par la protéine (et non un ion de surface)."""
    if not metals:
        return None
    P = np.array([[a[1], a[2], a[3]] for a in poly])
    el = [a[0] for a in poly]
    best = None
    for m in metals:
        c = np.array(m[1:4], dtype=float)
        d = np.linalg.norm(P - c, axis=1)
        nd = sum(1 for k in np.where(d < 2.8)[0] if el[k] in DONOR)
        if best is None or nd > best[0]:
            best = (nd, c, m[0])
    return best


def measure(poly, centre):
    """Fractions ouvertes (par seuil) et fraction débouchant sur le solvant."""
    P = np.array([[a[1], a[2], a[3]] for a in poly])
    dirs = fib_sphere()
    probes = centre + dirs * COORD_R
    near = np.array([np.min(np.linalg.norm(P - p, axis=1)) for p in probes])
    openfrac = {t: float((near >= t).mean()) for t in THRESHOLDS}

    steps = np.arange(COORD_R, ESCAPE_R, 0.5)
    esc = 0
    for k in np.where(near >= THRESHOLDS[1])[0]:
        pts = centre + dirs[k] * steps[:, None]
        if all(np.min(np.linalg.norm(P - q, axis=1)) >= ESCAPE_CLEAR for q in pts):
            esc += 1
    return openfrac, esc / len(dirs)


def kabsch(mobile, target):
    """Rotation+translation optimale amenant `mobile` sur `target` (RMSD minimal)."""
    cm, ct = mobile.mean(0), target.mean(0)
    H = (mobile - cm).T @ (target - ct)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    rmsd = float(np.sqrt((((mobile - cm) @ R.T + ct - target) ** 2).sum(1).mean()))
    return R, cm, ct, rmsd
