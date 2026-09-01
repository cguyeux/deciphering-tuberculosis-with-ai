#!/usr/bin/env python3
"""
P2.6 - Contrôle NÉGATIF du site métal : triple mutant C113A/H115A/E59A + ion (Fe/Zn), AF3 holo.

But (lever la circularité, review 2026-07-06) : si l'ion se loge ENCORE dans la région du site avec un
pLDDT aussi élevé qu'en WT, le signal holo était un artefact ; s'il ne se loge plus / pLDDT effondré /
distances aberrantes, le site Cys113-His115-Glu59 est spécifique.

Vérifie les résidus WT aux positions avant de muter (garde-fou anti-décalage).
Sortie : résultats/af3_mutant_control_jobs.json (à importer sur alphafoldserver.com).
"""
import json, os
from Bio import SeqIO

ROOT = "/home/christophe/docs/codes/mtbc/Rv1025"
OPERON_FAA = f"{ROOT}/résultats/operon_proteins.faa"
OUT = f"{ROOT}/résultats/af3_mutant_control_jobs.json"
MUTS = [(59, "E", "A"), (113, "C", "A"), (115, "H", "A")]  # (pos 1-indexée, WT, mut)

def rv1025_seq():
    for rec in SeqIO.parse(OPERON_FAA, "fasta"):
        if rec.id.split("_")[0] == "Rv1025":
            return str(rec.seq)
    raise SystemExit("Rv1025 introuvable")

def mutate(seq):
    s = list(seq)
    for pos, wt, mut in MUTS:
        assert s[pos - 1] == wt, f"attendu {wt}{pos}, trouvé {s[pos-1]}{pos} — DÉCALAGE, on n'écrit rien"
        s[pos - 1] = mut
    return "".join(s)

def job(name, seq, ion):
    return {"name": name, "modelSeeds": [1],
            "sequences": [{"proteinChain": {"sequence": seq, "count": 1}},
                          {"ion": {"ion": ion, "count": 1}}],
            "dialect": "alphafoldserver", "version": 1}

def main():
    wt = rv1025_seq()
    mut = mutate(wt)
    ndiff = sum(a != b for a, b in zip(wt, mut))
    print(f"WT 155 aa ; mutant {ndiff} substitutions : " + ", ".join(f"{w}{p}{m}" for p, w, m in MUTS))
    jobs = [job("Rv1025_mut3A_Fe_ctrl", mut, "FE"), job("Rv1025_mut3A_Zn_ctrl", mut, "ZN")]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(jobs, open(OUT, "w"), indent=2)
    print(f"{len(jobs)} jobs écrits -> {OUT}")
    print("Lecture (parseur phase7 réutilisable) : comparer au WT (Fe pLDDT 98,3, triade 3/3 à 2,3 Å). Attendu si "
          "site spécifique : ion NON logé dans C113/H115/E59 (devenus Ala) et/ou pLDDT ion effondré.")

if __name__ == "__main__":
    main()
