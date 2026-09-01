#!/usr/bin/env python3
"""
P4.1 - Génère les jobs AlphaFold3 Server pour tester le couplage Rv1025 <-> divIC.

Design (test + contrôles de spécificité) :
  J1 Rv1025 + divIC (Rv1024)  -> TEST (hypothèse : partenaire de division)
  J2 Rv1025 + ppx2  (Rv1026)  -> contrôle voisin d'opéron (voie polyphosphate)
  J3 Rv1025 + eno   (Rv1023)  -> contrôle voisin d'opéron (glycolyse)
  J4 divIC  + ftsQ  (Rv2151c) -> contrôle POSITIF (interface divisome FtsB-FtsQ connue)
  J5 Rv1025 + ftsQ  (Rv2151c) -> Rv1025 touche-t-il le divisome élargi ?

Interprétation : ipTM(J1) >> ipTM(J2,J3) et PAE inter-chaînes basse = interface Rv1025-divIC
spécifique. J4 calibre le "vrai" ipTM d'une interface de division. Garde-fou : AF peut produire
un faux positif d'interface ; la spécificité relative (divIC vs eno/ppx2) et la cohérence des 5
modèles sont les critères, pas la valeur absolue seule.

Sortie : résultats/af3_jobs_p4.json (uploadable sur https://alphafoldserver.com, "Add jobs from JSON").
"""
import json, os
from Bio import SeqIO

ROOT = "/home/christophe/docs/codes/mtbc/Rv1025"
OPERON_FAA = f"{ROOT}/résultats/operon_proteins.faa"
OUT = f"{ROOT}/résultats/af3_jobs_p4.json"

# ftsQ (Rv2151c) extrait du CDS H37Rv (contrôle positif divisome)
FTSQ = ("MTEHNEDPQIERVADDAADEEAVTEPLATESKDEPAEHPEFEGPRRRARRERAERRAAQARATAIEQARRAAKRRARGQIVSEQNPAKP"
        "AARGVVRGLKALLATVVLAVVGIGLGLALYFTPAMSAREIVIIGIGAVSREEVLDAARVRPATPLLQIDTQQVADRVATIRRVASARVQ"
        "RQYPSALRITIVERVPVVVKDFSDGPHLFDRDGVDFATDPPPPALPYFDVDNPGPSDPTTKAALQVLTALHPEVASQVGRIAAPSVASIT"
        "LTLADGRVVIWGTTDRCEEKAEKLAALLTQPGRTYDVSSPDLPTVK")

def load_operon():
    seq = {}
    for rec in SeqIO.parse(OPERON_FAA, "fasta"):
        # header : Rv1023_eno, Rv1024_divIC, Rv1025_Rv1025, Rv1026_ppx2
        lt = rec.id.split("_")[0]
        seq[lt] = str(rec.seq)
    seq["Rv2151c"] = FTSQ
    return seq

def dimer(name, s1, s2):
    return {
        "name": name,
        "modelSeeds": [1],
        "sequences": [
            {"proteinChain": {"sequence": s1, "count": 1}},
            {"proteinChain": {"sequence": s2, "count": 1}},
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }

def main():
    s = load_operon()
    jobs = [
        dimer("Rv1025_divIC_TEST", s["Rv1025"], s["Rv1024"]),
        dimer("Rv1025_ppx2_ctrl", s["Rv1025"], s["Rv1026"]),
        dimer("Rv1025_eno_ctrl", s["Rv1025"], s["Rv1023"]),
        dimer("divIC_ftsQ_posctrl", s["Rv1024"], s["Rv2151c"]),
        dimer("Rv1025_ftsQ", s["Rv1025"], s["Rv2151c"]),
    ]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(jobs, fh, indent=2)
    print(f"{len(jobs)} jobs écrits -> {OUT}")
    for j in jobs:
        lens = [len(c["proteinChain"]["sequence"]) for c in j["sequences"]]
        print(f"  {j['name']:22s} chaînes {lens} (total {sum(lens)} aa)")

if __name__ == "__main__":
    main()
