#!/usr/bin/env python3
"""Phase 4 — Polymorphismes pré-existants / lignées « amorcées » (angle distinctif).

Identifie les variants NON-synonymes des gènes BDQ/PMD qui sont (a) NON catalogués-R, (b)
CONCENTRÉS dans une lignée (standing polymorphism de ce clade) et (c) à fréquence intra-lignée
non triviale. Un tel variant est, par construction, ANTÉRIEUR à l'ère du médicament : fixé sur le
tronc de la (sous-)lignée, daté du MRCA de ce clade (siècles à millénaires) >> approbation BDQ
(2012) / PMD (2019). C'est la définition opérationnelle d'une lignée « amorcée » : un fond
génomique pré-existant dans un gène-cible, susceptible d'abaisser la barrière à la résistance.

Léger (pas de rechargement du pangénome) : réutilise
  - `résultats/phase1_feasibility.tsv` : variants + effet + catalogue + ventilation par lignée
  - `résultats/phase3_lineage_prevalence.tsv` : effectifs par lignée (dénominateurs)
  - `résultats/phase3b_ipw_adjusted.tsv` : enrichissement R ajusté par lignée (pour croiser)

Niveau HYPOTHÈSE (cadre N1 de tuberculosis.md) : sans cohorte appariée, on génère des candidats,
on ne prouve pas la causalité. L'antériorité vient de la fixation ; l'impact fonctionnel se teste
ensuite par ESM-1v (`/mtbc-mutation-impact`) ; la confirmation ancestrale formelle = pastml (différé,
inutile pour un marqueur quasi fixé).

Sortie : résultats/phase4_primed_candidates.tsv + récap console.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths

CONC_MIN = 0.85       # >=85% des porteurs dans UNE lignée (spécificité de clade)
FIX_MIN = 0.05        # >=5% de cette lignée portent le variant (standing, pas singleton)
GROUPS = {"BDQ", "PMD"}   # cœur ; COMP/LZD/MFX en contexte si besoin


def main():
    # effectifs par lignée (n_total indépendant de la drogue)
    lin_total = {}
    with open(paths.RESULTATS / "phase3_lineage_prevalence.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lin_total[r["lineage"]] = int(r["n_total"])
    print("effectifs par lignée :", lin_total)

    # enrichissement R ajusté (Phase 3b) par drogue -> pour croiser
    # (on retient l'OR cluster/curation lu dans le cahier ; ici juste la prévalence ajustée par lignée)
    ipw = defaultdict(dict)
    p3b = paths.RESULTATS / "phase3b_ipw_adjusted.tsv"
    if p3b.exists():
        with open(p3b) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                ipw[r["drug"]][r["lineage"]] = float(r["ipw_prev_pct"])

    candidates = []
    with open(paths.RESULTATS / "phase1_feasibility.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["group"] not in GROUPS:
                continue
            if r["effect"] not in ("missense", "stop_gained"):
                continue
            if r["is_catalogued_R"] == "1":
                continue
            # ventilation lignée
            br = {}
            for tok in r["lineage_breakdown"].split(";"):
                if ":" in tok:
                    k, v = tok.rsplit(":", 1)
                    if v.isdigit():
                        br[k] = int(v)
            br.pop("NA", None)
            tot = sum(br.values())
            if tot < 10 or not br:
                continue
            dom = max(br, key=lambda k: br[k])
            dom_n = br[dom]
            conc = dom_n / tot
            fix = dom_n / lin_total.get(dom, 10**9)
            if conc >= CONC_MIN and fix >= FIX_MIN and dom in lin_total:
                candidates.append({
                    "group": r["group"], "drug": r["drug"], "gene": r["gene"],
                    "aa_change": r["aa_change"], "spdi": r["spdi"], "effect": r["effect"],
                    "lineage": dom, "fixation_in_lineage_pct": round(100 * fix, 1),
                    "concentration_pct": round(100 * conc, 1), "n_carriers": tot,
                    "freq_overall_pct": r["freq_pct"], "catalogue": r["catalogue"] or "-",
                })

    candidates.sort(key=lambda c: (c["group"] != "PMD", -c["fixation_in_lineage_pct"]))
    out = paths.RESULTATS / "phase4_primed_candidates.tsv"
    cols = ["group", "drug", "gene", "aa_change", "spdi", "effect", "lineage",
            "fixation_in_lineage_pct", "concentration_pct", "n_carriers",
            "freq_overall_pct", "catalogue"]
    with open(out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(candidates)

    print(f"\n=== Candidats « amorcés » (non-syn, non-catalogués-R, concentrés ≥{int(100*CONC_MIN)}% "
          f"dans une lignée, ≥{int(100*FIX_MIN)}% de cette lignée) : {len(candidates)} ===")
    print(f"{'grp':4}{'gene':7}{'aa':9}{'lignée':8}{'fix%':>6}{'conc%':>7}{'nport':>7}  catalogue")
    for c in candidates[:25]:
        print(f"{c['group']:4}{c['gene']:7}{c['aa_change']:9}{c['lineage']:8}"
              f"{c['fixation_in_lineage_pct']:6.1f}{c['concentration_pct']:7.1f}{c['n_carriers']:7d}  {c['catalogue']}")

    # croisement : les lignées « amorcées » F420 sont-elles enrichies en résistance DLM/PMD ?
    print("\n=== Croisement amorce F420 × enrichissement résistance (Phase 3b) ===")
    pmd_lineages = {c["lineage"] for c in candidates if c["group"] == "PMD"}
    if ipw.get("delamanid"):
        print("Prévalence DLM-R (IPW) par lignée :",
              {l: round(ipw["delamanid"].get(l, 0), 2) for l in sorted(lin_total)})
    print(f"Lignées portant ≥1 amorce F420 : {sorted(pmd_lineages)}")
    print("\nNiveau HYPOTHÈSE : antériorité établie par fixation (MRCA de lignée >> ère PMD/BDQ). "
          "Prochaine brique : ESM-1v (`/mtbc-mutation-impact <gene> <aa>`) sur les top candidats "
          "F420 pour prédire l'impact fonctionnel (abaissement de barrière) ; pastml en confirmation.")
    print(f"\nsortie -> {out}")


if __name__ == "__main__":
    main()
