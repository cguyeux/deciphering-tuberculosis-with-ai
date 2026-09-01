#!/usr/bin/env python3
# pyright: reportOptionalSubscript=false, reportAttributeAccessIssue=false
# (stubs Bio.PDB : `residue["CB"]` est typé Optional et `.sasa`, posé par ShrakeRupley
#  au moment du calcul, est inconnu du stub. Idiome documenté, vérifié à l'exécution.)
"""phase23_invariant_clusters.py — P2.3.a: what are the OTHER near-invariant residues?

P2.3 found 23 near-invariant residues (>=90% identity over the 8,700-sequence family
alignment) but the project only ever exploited THREE of them (the Cys113-His115-Glu59
metal triad). The remaining twenty were never looked at. Two very different
possibilities, and the answer changes what they mean:

  * they are BURIED and scattered -> ordinary fold-core constraint, no functional
    information (this is the null, and the most likely outcome for a small protein);
  * they form an EXPOSED spatial cluster -> a conserved surface patch, i.e. a
    candidate protein-interaction site — which would tell P4.3 (divisome panel)
    where to look instead of screening partners blind.

Discriminant = relative solvent accessibility (Shrake-Rupley ASA normalised by the
Tien et al. 2013 theoretical maxima), combined with single-linkage spatial clustering.

Reads : résultats/af3_out/.../..._unpaired_msa_chains_a.a3m   (conservation)
        résultats/structure/AF-P96375-F1.pdb                  (apo AF model)
        résultats/druggability/AF-P96375-F1_out/pockets/pocket{1,3}_atm.pdb
Writes: résultats/invariant_clusters.tsv
"""
from pathlib import Path
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.SASA import ShrakeRupley

ROOT = Path(__file__).resolve().parent.parent
A3M = ROOT / "résultats/af3_out/fold_rv1025_divic_test/msas/fold_rv1025_divic_test_unpaired_msa_chains_a.a3m"
PDB = ROOT / "résultats/structure/AF-P96375-F1.pdb"
POCK = ROOT / "résultats/druggability/AF-P96375-F1_out/pockets"
OUT = ROOT / "résultats/invariant_clusters.tsv"

CONS_CUT = 90.0     # "near-invariant", same threshold as P2.3 / the manuscript
LINK_CUT = 10.0     # single-linkage on CB (CA for Gly), as in the site-clustering step
BURIED, EXPOSED = 0.20, 0.40
TRIAD = {59, 113, 115}

# Tien et al. 2013 (theoretical) maximum accessible surface areas, A^2
MAXASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225,
          "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,
          "MET": 224, "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285,
          "TYR": 263, "VAL": 174}


def conservation():
    seqs, cur = [], []
    for line in A3M.read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur)); cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    aligned = ["".join(c for c in s if not c.islower()) for s in seqs]
    q = aligned[0]
    out = []
    for i, a in enumerate(q):
        col = [s[i] for s in aligned if i < len(s)]
        ng = [c for c in col if c not in "-."]
        out.append(100.0 * sum(1 for c in ng if c == a) / len(ng) if ng else 0.0)
    return q, out


def pocket_res(n):
    f = POCK / f"pocket{n}_atm.pdb"
    return {int(l[22:26]) for l in f.read_text().splitlines() if l.startswith("ATOM")}


query, cons = conservation()
struct = PDBParser(QUIET=True).get_structure("m", str(PDB))
ShrakeRupley().compute(struct[0], level="R")
chain = list(struct[0])[0]
res_by_id = {r.id[1]: r for r in chain if r.id[0] == " "}

inv = [i + 1 for i, c in enumerate(cons) if c >= CONS_CUT]
print(f"Résidus quasi invariants (>={CONS_CUT:.0f}%) : {len(inv)}  -> {inv}")

# RSA per invariant residue
rsa = {}
for p in inv:
    r = res_by_id.get(p)
    if r is None:
        continue
    rsa[p] = r.sasa / MAXASA.get(r.resname.strip(), 200)

# single-linkage clustering on CB (CA for Gly)
def anchor(p):
    r = res_by_id[p]
    return (r["CB"] if "CB" in r else r["CA"]).coord


def dist(a, b):
    pa, pb = anchor(a), anchor(b)
    return sum((pa[k] - pb[k]) ** 2 for k in range(3)) ** 0.5


clusters: list = []
for p in inv:
    if p not in res_by_id:
        continue
    joined = [c for c in clusters if any(dist(p, q) <= LINK_CUT for q in c)]
    if not joined:
        clusters.append([p])
        continue
    merged = {p} | {x for c in joined for x in c}
    for c in joined:
        clusters.remove(c)
    clusters.append(sorted(merged))

clusters.sort(key=len, reverse=True)
p1, p3 = pocket_res(1), pocket_res(3)

rows: list = [("cluster", "n", "residues", "mean_RSA", "burial", "contains_triad",
               "overlap_pocket1", "overlap_pocket3", "interpretation")]
print(f"\n{'cl':>3} {'n':>3} {'RSA moy':>8} {'enfoui/exposé':>14}  résidus")
print("-" * 100)
for i, c in enumerate(clusters, 1):
    vals = [rsa[p] for p in c if p in rsa]
    m = sum(vals) / len(vals) if vals else float("nan")
    burial = "enfoui" if m < BURIED else ("exposé" if m > EXPOSED else "intermédiaire")
    tri = sorted(set(c) & TRIAD)
    o1, o3 = sorted(set(c) & p1), sorted(set(c) & p3)
    if tri:
        interp = "site métal (connu)"
    elif burial == "enfoui":
        interp = "coeur du repli — contrainte structurale, PAS fonctionnelle"
    elif burial == "exposé" and len(c) >= 3:
        interp = "PATCH DE SURFACE conservé — site d'interaction candidat"
    else:
        interp = "isolé / intermédiaire — non concluant"
    rows.append((i, len(c), ";".join(map(str, c)), f"{m:.2f}", burial,
                 ";".join(map(str, tri)) or "-", ";".join(map(str, o1)) or "-",
                 ";".join(map(str, o3)) or "-", interp))
    print(f"{i:>3} {len(c):>3} {m:>8.2f} {burial:>14}  {c}"
          f"{'  [triade ' + str(tri) + ']' if tri else ''}"
          f"{'  [poche1 ' + str(o1) + ']' if o1 else ''}")

OUT.write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")

n_bur = sum(1 for p in rsa if rsa[p] < BURIED)
print(f"\nBilan : {n_bur}/{len(rsa)} invariants sont ENFOUIS (RSA<{BURIED}) ; "
      f"{sum(1 for p in rsa if rsa[p] > EXPOSED)} exposés (RSA>{EXPOSED}).")
print(f"Écrit : {OUT.relative_to(ROOT)}")
print("\nLecture : un amas d'invariants ENFOUIS = cœur du repli, sans portée fonctionnelle (modèle nul).")
print("Seul un amas EXPOSÉ de >=3 invariants désignerait une surface fonctionnelle candidate,")
print("qui orienterait le panel de partenaires de P4.3 au lieu de le chercher à l'aveugle.")
