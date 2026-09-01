#!/usr/bin/env python3
"""phase25_conservation_weighted.py — P2.3.a (deepening): redundancy-corrected conservation.

Two limitations were left open by phase22/phase23:
  (a) conservation was computed on the 8,700-sequence AF3 MSA, which is redundant;
  (b) the intended cross-check on the curated Pfam seed was impossible because
      Rv1025 is not among the 39 seed sequences.

Recon then showed something more decisive: the PF04417 domain spans residues 17-136
of Rv1025, so the exposed conserved patch (6-7-8) and part of Pocket 1 lie OUTSIDE the
annotated domain — a seed-based check could never have covered them, whatever the
alignment trick. The right fix is therefore not the seed but a redundancy correction
applied to the full-length MSA itself.

Method: Henikoff & Henikoff (1994) position-based sequence weights, the standard
remedy for over-represented sequences, then weighted per-column conservation. If the
Pocket-1 vs background contrast (74% vs 54%, p=0.005 unweighted) survives, the P5.1.a
conclusion is robust to redundancy; if it collapses, it was an artefact.

Reads : résultats/af3_out/.../..._unpaired_msa_chains_a.a3m
        résultats/druggability/AF-P96375-F1_out/pockets/pocket{1,3}_atm.pdb
Writes: résultats/conservation_weighted.tsv
"""
import random
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
A3M = ROOT / "résultats/af3_out/fold_rv1025_divic_test/msas/fold_rv1025_divic_test_unpaired_msa_chains_a.a3m"
POCK = ROOT / "résultats/druggability/AF-P96375-F1_out/pockets"
OUT = ROOT / "résultats/conservation_weighted.tsv"
DOMAIN = (17, 136)          # PF04417 on Rv1025 (atlas: ali_from/ali_to)
TRIAD = {59, 113, 115}
PATCH = {6, 7, 8}
N_NULL = 10000
random.seed(42)


def read_a3m(path):
    seqs, cur = [], []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur)); cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return ["".join(c for c in s if not c.islower()) for s in seqs]


def henikoff_weights(aligned, L):
    """Position-based sequence weights (Henikoff & Henikoff 1994)."""
    w = [0.0] * len(aligned)
    for i in range(L):
        col = [s[i] if i < len(s) else "-" for s in aligned]
        counts = Counter(c for c in col if c not in "-.")
        r = len(counts)                       # distinct residue types in this column
        if r == 0:
            continue
        for k, c in enumerate(col):
            if c in "-.":
                continue
            w[k] += 1.0 / (r * counts[c])
    tot = sum(w)
    return [x / tot * len(w) for x in w] if tot else w   # mean weight = 1


def weighted_conservation(query, aligned, w):
    cons = []
    for i, q in enumerate(query):
        num = den = 0.0
        for k, s in enumerate(aligned):
            if i >= len(s):
                continue
            c = s[i]
            if c in "-.":
                continue
            den += w[k]
            if c == q:
                num += w[k]
        cons.append(100.0 * num / den if den else 0.0)
    return cons


def pocket_res(n):
    f = POCK / f"pocket{n}_atm.pdb"
    return sorted({int(l[22:26]) for l in f.read_text().splitlines() if l.startswith("ATOM")})


aligned = read_a3m(A3M)
query = aligned[0]
L = len(query)
print(f"MSA : {len(aligned)} séquences, {L} colonnes")
print("Calcul des poids de Henikoff ...")
w = henikoff_weights(aligned, L)
eff = sum(w) ** 2 / sum(x * x for x in w)      # effective number of sequences
print(f"  nombre EFFECTIF de séquences : {eff:.0f} (sur {len(aligned)}) "
      f"-> redondance ~{len(aligned)/eff:.1f}x")

cons_w = weighted_conservation(query, aligned, w)
cons_u = []                                     # unweighted, for comparison
for i, q in enumerate(query):
    col = [s[i] for s in aligned if i < len(s)]
    ng = [c for c in col if c not in "-."]
    cons_u.append(100.0 * sum(1 for c in ng if c == q) / len(ng) if ng else 0.0)

p1, p3 = pocket_res(1), pocket_res(3)
bg_w = sum(cons_w) / L
bg_u = sum(cons_u) / L


def mean(res, cons):
    v = [cons[r - 1] for r in res if 1 <= r <= L]
    return sum(v) / len(v)


def null_p(obs, k, cons):
    idx = list(range(L))
    hits = sum(1 for _ in range(N_NULL)
               if sum(cons[i] for i in random.sample(idx, k)) / k >= obs)
    return (hits + 1) / (N_NULL + 1)


rows: list = [("set", "n", "mean_unweighted", "mean_henikoff", "p_henikoff", "in_PF04417_domain")]
print(f"\n{'set':<24} {'non pondéré':>12} {'Henikoff':>10} {'p':>8}  domaine PF04417 (17-136)")
print("-" * 82)
for name, res in (("Poche 1 (druggable)", p1), ("Poche 3 (site métal)", p3),
                  ("triade métal", sorted(TRIAD)), ("patch exposé 6-7-8", sorted(PATCH))):
    mu, mw = mean(res, cons_u), mean(res, cons_w)
    p = null_p(mw, len(res), cons_w)
    inside = sum(1 for r in res if DOMAIN[0] <= r <= DOMAIN[1])
    print(f"{name:<24} {mu:>11.1f}% {mw:>9.1f}% {p:>8.4f}  {inside}/{len(res)} dans le domaine")
    rows.append((name, len(res), f"{mu:.1f}", f"{mw:.1f}", f"{p:.4f}", f"{inside}/{len(res)}"))
print(f"{'fond protéine':<24} {bg_u:>11.1f}% {bg_w:>9.1f}% {'-':>8}")
rows.append(("whole_protein_background", L, f"{bg_u:.1f}", f"{bg_w:.1f}", "-", "-"))

OUT.write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")
print(f"\nÉcrit : {OUT.relative_to(ROOT)}")
print("\nLecture : si l'écart poche/fond SURVIT à la pondération, la conclusion de P5.1.a")
print("est robuste à la redondance du MSA. La colonne « domaine » rappelle qu'une partie")
print("des résidus est HORS PF04417 (17-136), donc hors de portée d'un contrôle sur le seed.")
