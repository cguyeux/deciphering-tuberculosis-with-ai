#!/usr/bin/env python3
"""
P2.3 - Conservation par résidu de Rv1025 (DUF501) + cluster spatial = site fonctionnel candidat.

Exploite le MSA profond généré par AF3 (8700 séq, chaîne A = Rv1025) pour calculer la conservation
par position, puis mappe les résidus quasi invariants sur le modèle AlphaFold pour voir s'ils forment
un amas spatial (signature d'un site fonctionnel : liaison / catalyse), même sans hit de repli (Foldseek nul).

Entrées : MSA a3m (chaîne A), modèle AF (AF-P96375-F1.pdb).
Sortie : liste des résidus les plus conservés + amas spatial + résidus « catalytiques » (C/H/D/E/S/K/R/Y).
"""
import math, os
from collections import Counter

ROOT = "/home/christophe/docs/codes/mtbc/Rv1025"
A3M = f"{ROOT}/résultats/af3_out/fold_rv1025_divic_test/msas/fold_rv1025_divic_test_unpaired_msa_chains_a.a3m"
PDB = f"{ROOT}/résultats/structure/AF-P96375-F1.pdb"

def read_a3m_matchcols(path):
    """Retourne (query, list_of_aligned_seqs) sur les colonnes de match (longueur = len(query))."""
    seqs = []
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs.append("".join(buf))
                name, buf = line, []
            else:
                buf.append(line)
    if name is not None:
        seqs.append("".join(buf))
    # a3m : minuscules = insertions vs query -> les retirer donne l'alignement aux colonnes query
    aligned = ["".join(c for c in s if not c.islower()) for s in seqs]
    query = aligned[0]
    L = len(query)
    aligned = [a for a in aligned if len(a) == L]  # garde-fou
    return query, aligned

def conservation(query, aligned):
    L = len(query)
    rows = []
    for i in range(L):
        col = [a[i] for a in aligned]
        nongap = [c for c in col if c != "-" and c != "."]
        n = len(nongap)
        if n == 0:
            rows.append((i + 1, query[i], 0.0, 0.0)); continue
        cnt = Counter(nongap)
        # identité à la query (fraction des séq non-gap portant le résidu query)
        idq = cnt.get(query[i], 0) / n
        # entropie de Shannon (bits), normalisée
        ent = -sum((c / n) * math.log2(c / n) for c in cnt.values())
        rows.append((i + 1, query[i], idq, ent))
    return rows

def load_ca(pdb):
    """resnum -> coord CA (numpy)."""
    import numpy as np
    ca = {}
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                resnum = int(line[22:26])
                ca[resnum] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return ca

def main():
    import numpy as np
    query, aligned = read_a3m_matchcols(A3M)
    print(f"MSA : {len(aligned)} séquences alignées, query {len(query)} aa (Rv1025)")
    rows = conservation(query, aligned)
    # résidus quasi invariants : identité à la query >= 0.90
    conserved = [r for r in rows if r[2] >= 0.90]
    conserved.sort(key=lambda r: -r[2])
    print(f"\nRésidus conservés >=90% identité à la query : {len(conserved)}")
    catal = set("CHDESKRY")
    print(f"{'pos':>4s} {'aa':>2s} {'%idQ':>5s} {'entropie':>8s} {'catal?':>6s}")
    for pos, aa, idq, ent in conserved:
        print(f"{pos:4d} {aa:>2s} {idq*100:5.0f} {ent:8.2f} {'*' if aa in catal else '':>6s}")
    # amas spatial des conservés (CA < 10 Å) — cherche le plus gros cluster
    ca = load_ca(PDB)
    cpos = [p for p, a, idq, e in conserved if p in ca]
    import itertools
    # graphe de proximité
    adj = {p: set() for p in cpos}
    for a, b in itertools.combinations(cpos, 2):
        if np.linalg.norm(ca[a] - ca[b]) < 10.0:
            adj[a].add(b); adj[b].add(a)
    # composantes connexes
    seen, comps = set(), []
    for p in cpos:
        if p in seen: continue
        stack, comp = [p], []
        while stack:
            x = stack.pop()
            if x in seen: continue
            seen.add(x); comp.append(x)
            stack.extend(adj[x] - seen)
        comps.append(sorted(comp))
    comps.sort(key=len, reverse=True)
    print(f"\nAmas spatiaux de résidus conservés (CA<10 Å) :")
    for comp in comps[:4]:
        types = [(p, dict((r[0], r[1]) for r in rows)[p]) for p in comp]
        cat = [f"{a}{p}" for p, a in types if a in catal]
        print(f"  cluster n={len(comp)} : {['%s%d'%(a,p) for p,a in types]}")
        if cat:
            print(f"      dont catalytiques potentiels : {cat}")
    # surface conservée globale : top-12 conservés catalytiques
    topcat = [(p, a, idq) for p, a, idq, e in conserved if a in catal][:15]
    print(f"\nRésidus conservés de type catalytique/fonctionnel (C/H/D/E/S/K/R/Y) : "
          f"{['%s%d(%.0f%%)'%(a,p,idq*100) for p,a,idq in topcat]}")

if __name__ == "__main__":
    main()
