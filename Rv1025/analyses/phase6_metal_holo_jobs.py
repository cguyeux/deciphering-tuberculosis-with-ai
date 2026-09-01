#!/usr/bin/env python3
"""
P2.4 - Génère les jobs AlphaFold3 Server holo (Rv1025 + ion métallique) pour tester le site
Cys113-His115-Glu59 candidat (P2.3 : donneurs convergents à ~2,3 Å = site métal mononucléaire prédit).

Recette KB (cas Rv3577, tuberculosis.md) : AF3 Server accepte des ions via l'entité
{"ion":{"ion":"ZN","count":1}}. On teste Zn puis Fe (les deux plus probables), + Mn comme comparateur.
Lecture (parseur dédié) : l'ion doit se placer à ~2.0-2.3 Å de Cys113-SG, His115-N et Glu59-Oε,
avec pLDDT d'ion élevé → confirme le site métallique et donc DUF501 = métalloprotéine.

Sortie : résultats/af3_metal_jobs.json (à importer sur alphafoldserver.com).
"""
import json, os
from Bio import SeqIO

ROOT = "/home/christophe/docs/codes/mtbc/Rv1025"
OPERON_FAA = f"{ROOT}/résultats/operon_proteins.faa"
OUT = f"{ROOT}/résultats/af3_metal_jobs.json"

def rv1025_seq():
    for rec in SeqIO.parse(OPERON_FAA, "fasta"):
        if rec.id.split("_")[0] == "Rv1025":
            return str(rec.seq)
    raise SystemExit("Rv1025 introuvable dans operon_proteins.faa")

def job(name, seq, ion):
    return {
        "name": name,
        "modelSeeds": [1],
        "sequences": [
            {"proteinChain": {"sequence": seq, "count": 1}},
            {"ion": {"ion": ion, "count": 1}},
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }

def main():
    s = rv1025_seq()
    jobs = [
        job("Rv1025_holo_Zn", s, "ZN"),
        job("Rv1025_holo_Fe", s, "FE"),
        job("Rv1025_holo_Mn", s, "MN"),
    ]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(jobs, fh, indent=2)
    print(f"{len(jobs)} jobs holo écrits -> {OUT}")
    for j in jobs:
        print(f"  {j['name']}  (Rv1025 155 aa + 1 {j['sequences'][1]['ion']['ion']})")
    print("Site candidat à vérifier : Cys113 / His115 / Glu59 (2e sphère K112, R92).")

if __name__ == "__main__":
    main()
