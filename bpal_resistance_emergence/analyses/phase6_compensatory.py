#!/usr/bin/env python3
"""Phase 6, Volet A — Fond compensatoire rpoC/rpoA et émergence de la résistance BPaL.

Hypothèse (épidémiologique). Les mutations rpoB (RIF-R) imposent un coût de fitness ; des mutations
compensatoires rpoC/rpoA le restaurent → clones MDR « fit », transmissibles, épidémiquement réussis.
La résistance aux nouvelles drogues (BDQ/DLM/LZD/FQ) émergeant sur CE fond compensé est plus
préoccupante (elle se propage). Question : la résistance BPaL est-elle enrichie sur fond compensé ?

Méthode (cohorte rifampicine du consensus Resistance_antibio : ~10 920 R / 41 232 S).
1. Identifier empiriquement les variants COMPENSATOIRES = rpoC/rpoA non-syn ENRICHIS chez les RIF-R
   (ils co-occurrent avec rpoB muté ; Fisher carrier × RIF R/S). Cross-check positions connues.
2. « Compensé » = souche portant >=1 variant compensatoire.
3. Émergence : parmi les souches phénotypées pour chaque drogue BPaL, comparer le taux de R entre
   compensé et non-compensé (Fisher), en GÉNÉRAL puis RESTREINT aux RIF-R (dénominateur pertinent).
4. Stratification lignée (le fond compensé est-il concentré en L2/Beijing ? cf. Phase 3b).

Garde-fous. (a) Le fond compensé marque un MDR ANCIEN/avancé (plus d'histoire de traitement → plus
d'exposition BDQ/PMD) : l'association est ÉPIDÉMIOLOGIQUE, pas un mécanisme moléculaire direct.
(b) Confond lignée (L2). (c) rpoC accumule des SNP d'artefact de contamination (KB 2025) → seuil de
porteurs + enrichissement filtrent ; on rapporte aussi le nb de BioProjects si suspicion. p nominaux.

Sortie : résultats/phase6_compensatory.tsv + récap console.
"""
import csv
import pickle
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths

BPAL = ["bedaquiline", "delamanid", "linezolid", "levofloxacin", "moxifloxacin"]
COMP_GENES = {"rpoC", "rpoA"}
# positions compensatoires rpoC/rpoA connues (de Vos 2013, Comas 2012 ; pour cross-check, non exhaustif)
KNOWN_COMP_RPOC = {164, 433, 442, 449, 452, 483, 491, 521, 525, 527, 1040, 1041, 1075}


def aa_pos(s):
    m = re.search(r"(\d+)", s or "")
    return int(m.group(1)) if m else None


def main():
    # --- phénotypes : rif + drogues BPaL ---
    pheno = defaultdict(dict)
    drugs_keep = {"rifampicin"} | set(BPAL)
    with open(paths.PHENOTYPES_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d = r["drug"].lower()
            p = r["phenotype"].strip().upper()
            if d in drugs_keep and p in ("R", "S"):
                pheno[r["strain_id"]][d] = p
    RIFR = {s for s, dd in pheno.items() if dd.get("rifampicin") == "R"}
    RIFS = {s for s, dd in pheno.items() if dd.get("rifampicin") == "S"}
    print(f"Cohorte rifampicine : {len(RIFR)} R / {len(RIFS)} S")

    # --- lignées ---
    lin = {}
    with open(paths.LINEAGE_SNAPSHOT) as fh:
        rd = csv.DictReader(fh)
        sk, lk = rd.fieldnames[0], rd.fieldnames[-1]
        for r in rd:
            lin[r[sk]] = r[lk]

    # --- variants rpoC/rpoA (phase1) ---
    comp_vars = []
    with open(paths.RESULTATS / "phase1_feasibility.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene"] in COMP_GENES and r["effect"] in ("missense", "stop_gained"):
                comp_vars.append({"gene": r["gene"], "aa": r["aa_change"],
                                  "pos": aa_pos(r["aa_change"]), "spdi": r["spdi"]})

    print("chargement pan_spdi…", flush=True)
    pan = pickle.load(open(paths.PAN_SPDI, "rb"))

    def carr(spdi):
        return set(pan.get(spdi, []))

    # --- (1) identification empirique des compensatoires (enrichis chez RIF-R) ---
    nR, nS = len(RIFR), len(RIFS)
    comp_tested = []
    for v in comp_vars:
        c = carr(v["spdi"])
        a = len(c & RIFR); b = len(c & RIFS)
        if a + b < 10:           # seuil de porteurs phénotypés rif
            continue
        cc = nR - a; d = nS - b
        orr, p = fisher_exact([[a, b], [cc, d]], alternative="greater")
        comp_tested.append({**v, "rifR": a, "rifS": b,
                            "OR": (round(orr, 2) if orr not in (float("inf"),) else 9999.0),
                            "p": p, "frac_rifR": round(a / (a + b), 3),
                            "known": v["pos"] in KNOWN_COMP_RPOC if v["gene"] == "rpoC" else False})
    comp_tested.sort(key=lambda x: x["p"])
    # ensemble compensatoire = enrichi RIF-R (OR>=3, p<0.01) -> co-occurrence rpoB forte
    compensatory = [v for v in comp_tested if (v["OR"] >= 3 and v["p"] < 1e-2)]
    comp_spdis = {v["spdi"] for v in compensatory}

    print(f"\n=== (1) Variants rpoC/rpoA enrichis chez les RIF-R (= compensatoires empiriques) ===")
    print(f"{len(compensatory)} variants compensatoires (OR>=3, p<1e-2) / {len(comp_tested)} testés")
    print(f"{'gene':6}{'aa':9}{'rifR':>6}{'rifS':>6}{'fracR':>7}{'OR':>8}{'p':>10}  connu?")
    for v in compensatory[:18]:
        print(f"{v['gene']:6}{v['aa']:9}{v['rifR']:6}{v['rifS']:6}{v['frac_rifR']:7.2f}"
              f"{v['OR']:8.1f}{v['p']:10.1e}  {'OUI' if v['known'] else ''}")

    # --- (2) souches compensées ---
    compensated = set()
    for spdi in comp_spdis:
        compensated |= carr(spdi)
    print(f"\nSouches portant >=1 variant compensatoire : {len(compensated)}")
    # part compensée parmi RIF-R
    compRIFR = compensated & RIFR
    print(f"  dont RIF-R : {len(compRIFR)} ({100*len(compRIFR)/max(len(RIFR),1):.1f}% des RIF-R) "
          f"-- validation : le compensatoire doit être ~exclusif aux RIF-R")
    print(f"  dont RIF-S : {len(compensated & RIFS)} (doit être faible)")

    # --- (3) émergence BPaL : compensé vs non-compensé ---
    rows = []
    print(f"\n=== (3) Résistance BPaL selon le fond compensatoire ===")
    for restrict, label in [(None, "TOUTES souches phénotypées"), (RIFR, "RESTREINT aux RIF-R")]:
        print(f"\n--- {label} ---")
        print(f"{'drogue':13}{'comp:R/N':>12}{'noncomp:R/N':>14}{'OR':>8}{'p':>10}")
        for drug in BPAL:
            tested = {s for s, dd in pheno.items() if drug in dd}
            if restrict is not None:
                tested &= restrict
            comp_t = tested & compensated
            ncomp_t = tested - compensated
            cr = sum(1 for s in comp_t if pheno[s][drug] == "R"); cn = len(comp_t)
            nr = sum(1 for s in ncomp_t if pheno[s][drug] == "R"); nn = len(ncomp_t)
            if cn == 0 or nn == 0:
                continue
            orr, p = fisher_exact([[cr, cn - cr], [nr, nn - nr]], alternative="greater")
            print(f"{drug:13}{f'{cr}/{cn}':>12}{f'{nr}/{nn}':>14}"
                  f"{(orr if orr!=float('inf') else 9999):8.2f}{p:10.1e}")
            rows.append({"restriction": label, "drug": drug, "comp_R": cr, "comp_N": cn,
                         "noncomp_R": nr, "noncomp_N": nn,
                         "OR": round(orr, 2) if orr != float("inf") else 9999.0, "p": p})

    # --- (4) lignées du fond compensé ---
    lc = Counter(lin.get(s, "?") for s in compensated)
    print(f"\n=== (4) Lignées des souches compensées (concentration L2 attendue) ===")
    for l, n in lc.most_common(10):
        print(f"  {l:10} {n:6} ({100*n/max(len(compensated),1):.1f}%)")

    with open(paths.RESULTATS / "phase6_compensatory.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=["restriction", "drug", "comp_R", "comp_N",
                                           "noncomp_R", "noncomp_N", "OR", "p"], delimiter="\t")
        w.writeheader(); w.writerows(rows)
    # table des compensatoires
    with open(paths.RESULTATS / "phase6_compensatory_variants.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "aa", "pos", "spdi", "rifR", "rifS",
                                           "frac_rifR", "OR", "p", "known"], delimiter="\t",
                           extrasaction="ignore")
        w.writeheader(); w.writerows(comp_tested)

    print("\nGarde-fous : fond compensé = MDR avancé (confond exposition/traitement + lignée L2) -> "
          "association ÉPIDÉMIOLOGIQUE, pas mécanisme. p nominaux. La validation (1)/(2) : les "
          "compensatoires doivent être quasi exclusifs aux RIF-R (sinon ce ne sont pas des compensatoires).")
    print(f"\nsortie -> {paths.RESULTATS}/phase6_compensatory.tsv (+ _variants.tsv)")


if __name__ == "__main__":
    main()
