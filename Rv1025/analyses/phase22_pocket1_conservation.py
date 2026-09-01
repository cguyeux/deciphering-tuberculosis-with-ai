#!/usr/bin/env python3
"""phase22_pocket1_conservation.py — P5.1.a: is the druggable Pocket 1 conserved, or incidental?

fpocket (P5.1) found that the ONLY clearly druggable cavity of Rv1025 (Pocket 1,
druggability 0.645, 361 A^3) is NOT the metal site: it sits 12.6 A away, lined by
N-terminal residues. The manuscript leans on it ("the fold additionally presents a
separate, clearly druggable pocket"), yet nothing told us whether that cavity is a
conserved feature of the DUF501 family (a candidate second functional site, and a
real alternative target) or an incidental surface dent of this particular fold.

Test: per-residue conservation of the Pocket-1 lining residues, compared with
  (a) the metal-site pocket (Pocket 3, positive control: known to be invariant),
  (b) the WHOLE-PROTEIN background — the guard-rail: a pocket must be judged
      against the protein's own conservation level, never in absolute terms,
  (c) an empirical null: 10,000 random residue sets of the same size, giving a
      p-value for "this set is more conserved than chance".

Both alignments are used, exactly as in phase5/phase12:
  - the deep AF3 MSA (8,700 sequences), and
  - the curated Pfam PF04417 seed (39 sequences), independent of MSA redundancy.

Reads : résultats/af3_out/fold_rv1025_divic_test/msas/..._unpaired_msa_chains_a.a3m
        résultats/pfam/PF04417_seed.sto
        résultats/druggability/AF-P96375-F1_out/pockets/pocket{1,3}_atm.pdb
Writes: résultats/druggability/pocket1_conservation.tsv
"""
import random
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
A3M = ROOT / "résultats/af3_out/fold_rv1025_divic_test/msas/fold_rv1025_divic_test_unpaired_msa_chains_a.a3m"
STO = ROOT / "résultats/pfam/PF04417_seed.sto"
POCK = ROOT / "résultats/druggability/AF-P96375-F1_out/pockets"
OUT = ROOT / "résultats/druggability/pocket1_conservation.tsv"
N_NULL = 10000
random.seed(42)


def read_a3m(path):
    """a3m -> (query, aligned seqs on query columns). Lowercase = insertions, dropped."""
    seqs, cur = [], []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur)); cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    aligned = ["".join(c for c in s if not c.islower()) for s in seqs]
    return aligned[0], aligned


def read_stockholm(path):
    """Pfam seed -> aligned sequences (list of gapped strings)."""
    rows = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith(("#", "//")):
            continue
        parts = line.split()
        if len(parts) == 2:
            rows.setdefault(parts[0], []).append(parts[1])
    return ["".join(v) for v in rows.values()]


def conservation_per_query_pos(query, aligned):
    """% of non-gap sequences matching the query residue, per query position."""
    cons = []
    for i, q in enumerate(query):
        col = [s[i] for s in aligned if i < len(s)]
        nongap = [c for c in col if c not in "-."]
        cons.append(100.0 * sum(1 for c in nongap if c == q) / len(nongap) if nongap else float("nan"))
    return cons


def pocket_residues(n):
    f = POCK / f"pocket{n}_atm.pdb"
    return sorted({int(l[22:26]) for l in f.read_text().splitlines() if l.startswith("ATOM")})


# ---- conservation from the deep MSA -------------------------------------------------
query, aligned = read_a3m(A3M)
cons = conservation_per_query_pos(query, aligned)
L = len(query)
print(f"MSA : {len(aligned)} séquences, {L} colonnes de requête")

p1, p3 = pocket_residues(1), pocket_residues(3)
bg = [c for c in cons if c == c]                       # drop NaN
bg_mean = sum(bg) / len(bg)


def stats(res):
    vals = [cons[r - 1] for r in res if 1 <= r <= L]
    return sum(vals) / len(vals), vals


m1, v1 = stats(p1)
m3, v3 = stats(p3)

# ---- empirical null: random residue sets of the same size ---------------------------
def null_p(observed_mean, k):
    idx = list(range(L))
    hits = 0
    for _ in range(N_NULL):
        s = random.sample(idx, k)
        if sum(cons[i] for i in s) / k >= observed_mean:
            hits += 1
    return (hits + 1) / (N_NULL + 1)


p1_p = null_p(m1, len(p1))
p3_p = null_p(m3, len(p3))

# ---- cross-check on the curated Pfam seed (query = Rv1025 row, ungapped mapping) -----
seed_note = ""
try:
    seed = read_stockholm(STO)
    # the query row is the one matching the Rv1025 sequence once ungapped
    qrow = next((s for s in seed if re.sub(r"[-.]", "", s).replace(".", "")[:15] == query[:15]), None)
    if qrow:
        cols = [i for i, c in enumerate(qrow) if c not in "-."]   # query positions -> columns
        seed_cons = []
        for i in cols:
            col = [s[i] for s in seed]
            nongap = [c for c in col if c not in "-."]
            seed_cons.append(100.0 * sum(1 for c in nongap if c == qrow[i]) / len(nongap) if nongap else 0.0)
        s1 = [seed_cons[r - 1] for r in p1 if r <= len(seed_cons)]
        s3 = [seed_cons[r - 1] for r in p3 if r <= len(seed_cons)]
        sbg = sum(seed_cons) / len(seed_cons)
        seed_note = (f"seed PF04417 ({len(seed)} séq) : Poche1 {sum(s1)/len(s1):.1f}% | "
                     f"Poche3 {sum(s3)/len(s3):.1f}% | fond {sbg:.1f}%")
    else:
        seed_note = "seed PF04417 : ligne requête non appariée (mapping ignoré)"
except Exception as e:  # noqa: BLE001
    seed_note = f"seed PF04417 : non exploité ({e})"

# ---- report --------------------------------------------------------------------------
rows: list = [("set", "n_residues", "residues", "mean_conservation_pct", "empirical_p_vs_random")]
rows.append(("pocket1_druggable", len(p1), ";".join(map(str, p1)), f"{m1:.1f}", f"{p1_p:.4f}"))
rows.append(("pocket3_metal_site", len(p3), ";".join(map(str, p3)), f"{m3:.1f}", f"{p3_p:.4f}"))
rows.append(("whole_protein_background", L, "-", f"{bg_mean:.1f}", "-"))
OUT.write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")

print(f"\n{'set':<26} {'n':>3} {'moy. cons.':>11} {'p (vs aléatoire)':>17}")
print("-" * 62)
print(f"{'Poche 1 (druggable)':<26} {len(p1):>3} {m1:>10.1f}% {p1_p:>17.4f}")
print(f"{'Poche 3 (site métal)':<26} {len(p3):>3} {m3:>10.1f}% {p3_p:>17.4f}   <- contrôle positif")
print(f"{'fond protéine entière':<26} {L:>3} {bg_mean:>10.1f}% {'-':>17}")
print(f"\n{seed_note}")
print(f"\nÉcrit : {OUT.relative_to(ROOT)}")
print("\nLecture : la Poche 1 n'est un site fonctionnel candidat que si sa conservation")
print("dépasse NETTEMENT le fond de la protéine (et p faible). Sinon = cavité incidente,")
print("druggable mais sans argument fonctionnel — à ne pas présenter comme une 2e cible.")
