#!/usr/bin/env python3
"""
P2.4 - Analyse des jobs AF3 holo (Rv1025 + ion) : l'ion se loge-t-il dans le site Cys113-His115-Glu59 ?

Pour chaque métal (Zn/Fe/Mn) : meilleur modèle (ranking_score), localise l'ion, mesure ses distances aux
donneurs candidats (Cys113-Sγ, His115-Nδ1/Nε2, Glu59-Oε1/Oε2) + 2e sphère (Lys112-Nζ, Arg92), et lit le
pLDDT de l'ion. Coordination si distances ~2.0-2.6 Å à >=2-3 donneurs + pLDDT ion élevé.

Usage : python3 phase7_metal_parse.py <dir_af3_metal_out>
"""
import glob, json, os, sys, warnings
warnings.simplefilter("ignore")
import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser

DONORS = {  # résidu : liste d'atomes donneurs à tester
    ("CYS", 113): ["SG"], ("HIS", 115): ["ND1", "NE2"], ("GLU", 59): ["OE1", "OE2"],
    ("LYS", 112): ["NZ"], ("ARG", 92): ["NH1", "NH2", "NE"], ("ARG", 8): ["NH1", "NH2"],
}
METALS = {"ZN", "FE", "MN", "NI", "CO", "MG", "CA", "CU"}

def best_model(job_dir):
    best = None
    for js in glob.glob(f"{job_dir}/*summary_confidences_*.json"):
        d = json.load(open(js))
        idx = js.split("summary_confidences_")[1].split(".")[0]
        rank = d.get("ranking_score", 0)
        if best is None or rank > best[0]:
            cif = js.replace("summary_confidences", "model").replace(".json", ".cif")
            best = (rank, idx, cif if os.path.exists(cif) else None)
    return best

def analyse_job(job_dir, label):
    rank, idx, cif = best_model(job_dir)
    if not cif:
        print(f"[{label}] pas de CIF"); return
    st = MMCIFParser(QUIET=True).get_structure("m", cif)
    model = list(st)[0]
    # collecte atomes : donneurs protéiques + ion
    donor_atoms, ion = {}, None
    for ch in model:
        for res in ch:
            rn = res.resname.strip().upper()
            resid = res.id[1]
            if rn in METALS:
                for at in res:
                    ion = (rn, at.coord, at.get_bfactor())
            key = (rn, resid)
            if key in DONORS:
                for at in res:
                    if at.name in DONORS[key]:
                        donor_atoms[(rn, resid, at.name)] = at.coord
    if ion is None:
        print(f"[{label}] AUCUN ion trouvé dans le modèle {idx}"); return
    mname, mcoord, mpl = ion
    # distances ion -> donneurs
    dists = sorted(((np.linalg.norm(mcoord - c), k) for k, c in donor_atoms.items()), key=lambda x: x[0])
    coord_shell = [(d, k) for d, k in dists if d <= 2.8]
    print(f"[{label}] modèle {idx} (rank {rank:.2f}) | ion {mname} pLDDT={mpl:.1f}")
    print(f"    coordination (<=2.8 Å) : {len(coord_shell)} donneur(s)")
    for d, (rn, ri, an) in dists[:6]:
        flag = "  <== COORD" if d <= 2.8 else ""
        print(f"      {rn}{ri}:{an:4s} {d:5.2f} Å{flag}")
    site = {(rn, ri) for _, (rn, ri, _) in coord_shell}
    triad = {("CYS", 113), ("HIS", 115), ("GLU", 59)}
    n_triad = len(site & triad)
    verdict = ("SITE CONFIRMÉ" if n_triad >= 2 and mpl >= 50 else
               "partiel/incertain" if n_triad >= 1 else "ion NON logé dans le site prédit")
    print(f"    → {n_triad}/3 de la triade Cys113/His115/Glu59 coordonnent | pLDDT ion {mpl:.0f} → {verdict}\n")
    return dict(metal=mname, plddt=mpl, n_triad=n_triad, coord=len(coord_shell), verdict=verdict)

def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "résultats/af3_metal_out"
    jobs = sorted(d for d in glob.glob(f"{root}/fold_*") if os.path.isdir(d))
    print(f"Site candidat : Cys113 / His115 / Glu59 (2e sphère Lys112, Arg92, Arg8)\n")
    for j in jobs:
        analyse_job(j, os.path.basename(j).replace("fold_rv1025_holo_", "").upper())

if __name__ == "__main__":
    main()
