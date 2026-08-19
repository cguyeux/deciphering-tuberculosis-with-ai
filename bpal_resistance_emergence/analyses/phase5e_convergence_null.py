#!/usr/bin/env python3
"""Phase 5e (réponse review R05) — null de permutation pour la convergence positionnelle.

La force du signal convergent fgd1 210 = « 5 résidus alternatifs distincts à UN codon ». Est-ce
surprenant compte tenu de la densité de variation du gène (fgd1 vs ddn vs fbiC, longueurs différentes) ?
Null : pour chaque gène, on place aléatoirement ses variants non-synonymes sur ses codons (uniforme,
densité = nombre de variants préservé, identité du résidu alt préservée), et on mesure la multiplicité
maximale (n résidus alt distincts à un même codon). Sur n_perm tirages, p = P(un codon atteint >= la
multiplicité observée du candidat). Complément : rang empirique des candidats parmi TOUTES les positions
scannées (phase5_convergence_positions.tsv).

Caveat (énoncé) : null d'OPPORTUNITÉ UNIFORME par codon (ignore l'accessibilité mutationnelle réelle d'un
codon donné) → conservateur ; la vraie convergence est au moins aussi significative.

Sortie : résultats/phase5e_convergence_null.tsv + récap console.
"""
import csv
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths

random.seed(42)
NPERM = 20000
F420 = {"ddn", "fgd1", "fbiA", "fbiB", "fbiC", "fbiD"}
# positions candidates et leur multiplicité observée (n résidus alt distincts)
CANDIDATES = {("fgd1", 210): 5, ("fgd1", 10): 3, ("ddn", 111): 3}
AA = re.compile(r"p?\.?([A-Z*])(\d+)([A-Z*])")


def main():
    # --- longueurs aa des gènes du panel ---
    aalen = {}
    with open(paths.GENE_PANEL) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                aalen[r["gene"]] = max(1, (int(r["end"]) - int(r["start"]) + 1) // 3 - 1)
            except Exception:
                pass

    # --- variants non-syn par gène : (pos, alt) ---
    gene_vars = defaultdict(list)        # gene -> [(pos, alt)]
    obs_alts = defaultdict(lambda: defaultdict(set))  # gene -> pos -> {alt}
    with open(paths.RESULTATS / "phase1_feasibility.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene"] not in F420 or r["effect"] != "missense":
                continue
            m = AA.search(r["aa_change"] or "")
            if not m:
                continue
            pos, alt = int(m.group(2)), m.group(3)
            gene_vars[r["gene"]].append((pos, alt))
            obs_alts[r["gene"]][pos].add(alt)

    rows = []
    print(f"=== Null de permutation ({NPERM} tirages) — multiplicité de résidus alt par codon ===")
    print(f"{'gene':6}{'L_aa':>6}{'n_var':>7}{'obs_max':>8}  candidat(pos:obs->p)")
    for gene in sorted(F420):
        V = gene_vars.get(gene, [])
        L = aalen.get(gene, 0)
        if not V or L < 2:
            continue
        alts = [a for _, a in V]
        obs_max = max((len(s) for s in obs_alts[gene].values()), default=0)
        # null : placer chaque variant (avec son alt) sur un codon aléatoire
        targets = {pos: k for (g, pos), k in CANDIDATES.items() if g == gene}
        hits_max = 0
        hits_target = {pos: 0 for pos in targets}
        # distribution du max
        for _ in range(NPERM):
            codon_alts = defaultdict(set)
            for a in alts:
                codon_alts[random.randint(1, L)].add(a)
            m = max((len(s) for s in codon_alts.values()), default=0)
            if m >= obs_max:
                hits_max += 1
            for pos, k in targets.items():
                if m >= k:
                    hits_target[pos] += 1
        p_max = (hits_max + 1) / (NPERM + 1)
        cand_str = ", ".join(f"{pos}:{k}->p={(hits_target[pos]+1)/(NPERM+1):.4f}"
                             for pos, k in targets.items()) or "—"
        print(f"{gene:6}{L:6}{len(V):7}{obs_max:8}  {cand_str}")
        for pos, k in targets.items():
            rows.append({"gene": gene, "L_aa": L, "n_nonsyn_var": len(V), "obs_max_alts": obs_max,
                         "candidate_pos": pos, "candidate_obs_alts": k,
                         "p_codon_reaches_obs": round((hits_target[pos] + 1) / (NPERM + 1), 5)})
        if not targets:
            rows.append({"gene": gene, "L_aa": L, "n_nonsyn_var": len(V), "obs_max_alts": obs_max,
                         "candidate_pos": "", "candidate_obs_alts": "", "p_codon_reaches_obs": round(p_max, 5)})

    # --- complément : rang empirique parmi les positions scannées ---
    print("\n=== Rang empirique parmi toutes les positions convergentes scannées ===")
    allpos = []
    with open(paths.RESULTATS / "phase5_convergence_positions.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                allpos.append((r["gene"], int(r["pos"]), int(r["n_distinct_alt"]), int(r["n_lineages_sporadic"])))
            except Exception:
                pass
    ntot = len(allpos)
    for (g, pos), k in CANDIDATES.items():
        match = [x for x in allpos if x[0] == g and x[1] == pos]
        if not match:
            continue
        _, _, nalt, nlin = match[0]
        ge_alt = sum(1 for x in allpos if x[2] >= nalt)
        ge_both = sum(1 for x in allpos if x[2] >= nalt and x[3] >= nlin)
        print(f"  {g} {pos}: {nalt} alts / {nlin} lignées -> {ge_alt}/{ntot} positions ont >= {nalt} alts ; "
              f"{ge_both}/{ntot} ont >= {nalt} alts ET >= {nlin} lignées")

    with open(paths.RESULTATS / "phase5e_convergence_null.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "L_aa", "n_nonsyn_var", "obs_max_alts",
                                           "candidate_pos", "candidate_obs_alts", "p_codon_reaches_obs"],
                           delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print("\nLecture : p = probabilité qu'un codon du gène atteigne la multiplicité observée du candidat "
          "sous placement aléatoire des variants (densité préservée). p faible = convergence positionnelle "
          "non due au hasard / à la densité. Null d'opportunité uniforme (conservateur). "
          f"\nsortie -> {paths.RESULTATS}/phase5e_convergence_null.tsv")


if __name__ == "__main__":
    main()
