#!/usr/bin/env python3
"""phase28_duf501_family.py — P8.1 : la triade métal est-elle UNIVERSELLE à DUF501 ?

Question : le site Cys113-His115-Glu59 est-il une propriété de TOUTE la famille PF04417
(donc : le site EST la fonction de DUF501), ou d'un seul clade (donc : il existe des
sous-familles métal-indépendantes, ce qui serait une découverte en soi) ?

Données : alignement Pfam FULL de PF04417.18 (2 228 séquences, 471 colonnes), récupéré
depuis l'API InterPro. Il est DÉJÀ aligné — aucun aligneur requis, ce qui compte ici
puisque ni HMMER ni MAFFT ne sont installés. Rv1025 y figure (P96375_MYCTU/16-137), donc
les colonnes de la triade s'ancrent directement, sans réalignement ni hypothèse.

MODÈLES NULS ET PIÈGES ANTICIPÉS :
  1. « 100 % conservé » serait attendu et peu informatif ; ce sont les EXCEPTIONS qui
     portent l'information. On les caractérise donc au lieu de les traiter en bruit.
  2. Un résidu ABSENT peut être un vrai remplacement OU un simple trou d'alignement en
     bord de domaine. On distingue explicitement gap et substitution.
  3. Une substitution CONSERVATIVE (Glu->Asp, autre carboxylate) ne détruit pas un site
     métal. On compte donc séparément « triade stricte » et « triade chimiquement
     compatible » — conclure à une perte sur un E->D serait une faute.
  4. Une perte concentrée dans un clade ne prouve rien si ce clade est aussi celui qui
     est sur-représenté : on teste l'association perte x taxonomie, et l'association
     perte x sous-structure de séquence, avec un null de permutation.

Lit   : data/PF04417_full.sto
Écrit : résultats/duf501_family_triad.tsv, résultats/duf501_family_summary.txt,
        article/supplementary_materials/table_S9_duf501_family.tsv
"""
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fast_blas import ensure  # noqa: E402  # type: ignore[import-not-found]

ensure()

from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.cluster.hierarchy import fcluster, linkage  # noqa: E402
from scipy.spatial.distance import squareform  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
STO = ROOT / "data/PF04417_full.sto"
OUT = ROOT / "résultats"
SUPP = ROOT / "article/supplementary_materials"

ANCHOR = "P96375_MYCTU"          # Rv1025
TRIAD = {59: "E", 113: "C", 115: "H"}
# substitutions qui préservent la chimie de coordination
COMPATIBLE = {"E": set("ED"), "C": set("C"), "H": set("H")}
N_PERM = 10000
rng = np.random.default_rng(42)

out = []


def say(s=""):
    print(s)
    out.append(s)


# ------------------------------------------------------------------ lecture
seqs: dict = defaultdict(str)
order: list = []
for line in STO.read_text().splitlines():
    if not line or line.startswith("#") or line.startswith("//"):
        continue
    parts = line.split()
    if len(parts) != 2:
        continue
    name, frag = parts
    if name not in seqs:
        order.append(name)
    seqs[name] += frag

L = len(seqs[order[0]])
assert all(len(s) == L for s in seqs.values()), "alignement non rectangulaire"
say(f"Alignement PF04417 full : {len(order)} séquences x {L} colonnes")

# --------------------------------------------- ancrage des colonnes de la triade
anchor_name = next(n for n in order if n.startswith(ANCHOR))
start = int(anchor_name.split("/")[1].split("-")[0])
row = seqs[anchor_name]
col_of = {}
pos = start - 1
for c, ch in enumerate(row):
    if ch not in "-.":
        pos += 1
        if pos in TRIAD:
            col_of[pos] = c
say(f"Ancre : {anchor_name} (domaine {start}-{anchor_name.split('-')[-1]})")
for p, aa in TRIAD.items():
    assert p in col_of, f"résidu {p} hors du domaine aligné"
    got = row[col_of[p]]
    assert got.upper() == aa, f"colonne {col_of[p]} = {got}, attendu {aa} pour la position {p}"
    say(f"  {aa}{p} -> colonne d'alignement {col_of[p]}")

# ------------------------------------------------- état de la triade par séquence
def group_of(name):
    """Mnémonique d'organisme UniProt : '9ACTN' = Actinomycetota non classé, etc."""
    return name.split("_")[-1].split("/")[0]


COARSE = {"9ACTN": "Actinomycetota", "MYCTU": "Actinomycetota", "MYCTA": "Actinomycetota",
          "9MICO": "Actinomycetota", "9MICC": "Actinomycetota", "9PSEU": "Actinomycetota",
          "9NOCA": "Actinomycetota", "9CORY": "Actinomycetota", "9GAMM": "Pseudomonadota",
          "9BACT": "autre/inconnu"}


def coarse(name):
    g = group_of(name)
    if g in COARSE:
        return COARSE[g]
    if g.startswith(("MYC", "COR", "STR", "NOC", "MIC", "RHO", "ACT", "PSEU", "BIF", "ARTH")):
        return "Actinomycetota"
    if g.startswith(("PSE", "ESC", "VIB", "XAN", "9GAM", "9PRO", "9ENT", "9BUR", "9ALP")):
        return "Pseudomonadota"
    return "autre/inconnu"


rows = [("sequence", "groupe", "res59", "res113", "res115", "triade_stricte", "triade_compatible")]
status, compat, obs = {}, {}, {p: Counter() for p in TRIAD}
for n in order:
    s = seqs[n]
    aas = {p: s[col_of[p]].upper() for p in TRIAD}
    for p in TRIAD:
        obs[p][aas[p]] += 1
    strict = all(aas[p] == TRIAD[p] for p in TRIAD)
    ok = all(aas[p] in COMPATIBLE[TRIAD[p]] for p in TRIAD)
    status[n], compat[n] = strict, ok
    rows.append((n, coarse(n), aas[59], aas[113], aas[115],
                 "oui" if strict else "non", "oui" if ok else "non"))
(OUT / "duf501_family_triad.tsv").write_text("\n".join("\t".join(r) for r in rows) + "\n")

N = len(order)
say()
say("=" * 78)
say("A. LA TRIADE EST-ELLE UNIVERSELLE ?")
say("=" * 78)
for p in sorted(TRIAD):
    c = obs[p]
    gaps = c.get("-", 0) + c.get(".", 0)
    top = [f"{a}:{n} ({n/N:.1%})" for a, n in c.most_common(4) if a not in "-."]
    say(f"  position {TRIAD[p]}{p} (col {col_of[p]}) : {' | '.join(top)}"
        + (f" | gaps : {gaps} ({gaps/N:.1%})" if gaps else ""))
n_strict = sum(status.values())
n_compat = sum(compat.values())
say(f"TRIADE STRICTE (E/C/H exacts)          : {n_strict}/{N} = {n_strict/N:.1%}")
say(f"TRIADE CHIMIQUEMENT COMPATIBLE (E ou D) : {n_compat}/{N} = {n_compat/N:.1%}")
say(f"MEMBRES SANS triade compatible          : {N - n_compat} ({(N-n_compat)/N:.1%})")

# -------------------------------------------------- B. les exceptions : qui sont-elles ?
say()
say("=" * 78)
say("B. LES EXCEPTIONS : ARTEFACT D'ALIGNEMENT, OU VRAIE PERTE ?")
say("=" * 78)
losers = [n for n in order if not compat[n]]
gap_only = [n for n in losers
            if all(seqs[n][col_of[p]] in "-." or seqs[n][col_of[p]].upper() in COMPATIBLE[TRIAD[p]]
                   for p in TRIAD)]
say(f"{len(losers)} membres sans triade compatible, dont {len(gap_only)} où l'écart tient "
    f"UNIQUEMENT à des gaps (bord de domaine, séquence partielle)")
say(f"-> vraies SUBSTITUTIONS non conservatives : {len(losers) - len(gap_only)}")
if losers:
    say("Composition des positions perdues (hors gaps) :")
    for p in sorted(TRIAD):
        sub = Counter(seqs[n][col_of[p]].upper() for n in losers
                      if seqs[n][col_of[p]] not in "-."
                      and seqs[n][col_of[p]].upper() not in COMPATIBLE[TRIAD[p]])
        if sub:
            say(f"  {TRIAD[p]}{p} remplacé par : {', '.join(f'{a}x{c}' for a, c in sub.most_common(5))}")

# ------------------------------------- C. la perte est-elle liée à la taxonomie ?
say()
say("=" * 78)
say("C. LA PERTE EST-ELLE LIÉE À LA TAXONOMIE ?")
say("=" * 78)
tax = {n: coarse(n) for n in order}
say(f"{'groupe':<18} {'n':>6} {'triade compatible':>20}")
for g, cnt in Counter(tax.values()).most_common():
    ok = sum(1 for n in order if tax[n] == g and compat[n])
    say(f"  {g:<16} {cnt:>6} {ok:>10}/{cnt} = {ok/cnt:6.1%}")
acti = [n for n in order if tax[n] == "Actinomycetota"]
other = [n for n in order if tax[n] != "Actinomycetota"]
if other:
    table = [[sum(1 for n in acti if compat[n]), sum(1 for n in acti if not compat[n])],
             [sum(1 for n in other if compat[n]), sum(1 for n in other if not compat[n])]]
    odds, pv = fisher_exact(table)
    say(f"Fisher (Actinomycetota vs reste) : OR = {odds:.2f}, p = {pv:.3g}")

# -------------------------- D. sous-structure de séquence : perte groupée ou dispersée ?
say()
say("=" * 78)
say("D. SOUS-STRUCTURE : LA PERTE EST-ELLE CONCENTRÉE DANS UN CLADE ?")
say("=" * 78)
AA = "ACDEFGHIKLMNPQRSTVWY"
lut = np.full(256, 20, dtype=np.uint8)
for k, a in enumerate(AA):
    lut[ord(a)] = k
M = lut[np.frombuffer("".join(seqs[n].upper() for n in order).encode(), dtype=np.uint8)].reshape(N, L)
X = np.zeros((N, L * 21), dtype=np.float32)
X[np.arange(N)[:, None], np.arange(L)[None, :] * 21 + M] = 1.0
ident = (X @ X.T) / L
np.fill_diagonal(ident, 1.0)
dist = np.clip(1.0 - ident, 0, None)
Z = linkage(squareform(dist, checks=False), method="average")
lost = np.array([not compat[n] for n in order])
say(f"{'k':>3} {'plus gros amas de pertes':>26} {'p (permutation)':>18}")
for k in (5, 10, 20):
    lab = fcluster(Z, k, criterion="maxclust")
    best = max(Counter(lab[lost]).values()) if lost.any() else 0
    perm = np.array([max(Counter(lab[rng.permutation(lost)]).values()) for _ in range(N_PERM)])
    pv = float(((perm >= best).sum() + 1) / (N_PERM + 1))
    say(f"{k:>3} {best:>15}/{int(lost.sum()):<10} {pv:>18.4f}")
say("Lecture : p faible = les pertes se concentrent dans un même sous-groupe de séquences")
say("(sous-famille candidate) ; p élevé = elles sont dispersées, donc plutôt du bruit")
say("d'annotation ou des séquences partielles, pas un clade métal-indépendant.")

(OUT / "duf501_family_summary.txt").write_text("\n".join(out) + "\n")
SUPP.mkdir(parents=True, exist_ok=True)
tab = [("position", "residu_Rv1025", "colonne", "identique_%", "gaps_%", "principal_substituant")]
for p in sorted(TRIAD):
    c = obs[p]
    sub = [(a, n) for a, n in c.most_common() if a not in "-." and a != TRIAD[p]]
    tab.append((str(p), TRIAD[p], str(col_of[p]), f"{c[TRIAD[p]]/N*100:.1f}",
                f"{(c.get('-',0)+c.get('.',0))/N*100:.1f}",
                f"{sub[0][0]} ({sub[0][1]})" if sub else "-"))
tab.append(("triade", "E/C/H", "-", f"{n_strict/N*100:.1f}", "-",
            f"compatible E/D : {n_compat/N*100:.1f}%"))
(SUPP / "table_S9_duf501_family.tsv").write_text("\n".join("\t".join(r) for r in tab) + "\n")
say()
say("Écrit : résultats/duf501_family_triad.tsv, duf501_family_summary.txt,")
say("        article/supplementary_materials/table_S9_duf501_family.tsv")
