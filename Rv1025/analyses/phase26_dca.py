#!/usr/bin/env python3
"""phase26_dca.py — P8.2 : co-évolution (mean-field DCA) sur le MSA profond DUF501.

But (tel que reformulé dans pistes.md après challenge) : cibler les résidus qui
co-évoluent AVEC les ligands métal (Glu59, Cys113, His115) pour désigner le site
fonctionnel candidat au-delà de la triade — jamais pour assigner un substrat.

MODÈLE NUL EXPLICITE (le contre-argument à battre) : la co-évolution récupère
d'abord les CONTACTS STRUCTURAUX du repli. Une paire co-évoluante proche en 3D
n'apprend donc rien de fonctionnel. Deux garde-fous en conséquence :
  1. la précision top-L/2 contre la carte de contacts du modèle AF sert à VALIDER
     l'implémentation (si elle est proche du hasard, la formule est fausse) ;
  2. le signal réellement informatif est la fraction de couplages forts NON
     expliqués par le monomère (distance Cbeta > 15 A) : interface d'homo-oligomère
     (cf. P4.4) ou couplage fonctionnel/allostérique.

CIRCULARITÉ à ne pas oublier : AF3 a été nourri de CE MSA. L'accord DCA/modèle ne
vaut donc PAS comme preuve indépendante du repli — il ne vaut que comme contrôle
de l'implémentation. Les couplages non expliqués, eux, restent informatifs.

Méthode : Morcos et al. 2011 (mfDCA) — repondération à 80 % d'identité, pseudocompte
0.5, inversion de la matrice de covariance, jauge zéro-somme, norme de Frobenius,
correction APC (Dunn 2008).

Lit   : résultats/af3_out/.../..._unpaired_msa_chains_a.a3m
        résultats/structure/AF-P96375-F1.pdb
Écrit : résultats/dca_pairs.tsv, résultats/dca_triad_partners.tsv, résultats/dca_summary.txt
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "8")       # laisser du CPU au run Boltz en cours

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fast_blas import ensure                        # noqa: E402  # type: ignore[import-not-found]

ensure()            # le numpy système est lié à la BLAS netlib mono-thread (~38x plus lent)

from pathlib import Path                            # noqa: E402

import numpy as np                                  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
A3M = ROOT / "résultats/af3_out/fold_rv1025_divic_test/msas/fold_rv1025_divic_test_unpaired_msa_chains_a.a3m"
PDB = ROOT / "résultats/structure/AF-P96375-F1.pdb"
OUT = ROOT / "résultats"

AA = "ACDEFGHIKLMNPQRSTVWY"
Q = 21                       # 20 acides aminés + gap (état 20)
IDT = 0.8                    # seuil de repondération (standard DCA)
PC = 0.5                     # pseudocompte (standard mfDCA)
SEP = 5                      # |i-j| minimal pour la prédiction de contacts
CONTACT = 8.0                # seuil Cbeta-Cbeta d'un contact vrai (A)
FAR = 15.0                   # au-delà : couplage NON expliqué par le monomère
TRIAD = {59: "E", 113: "C", 115: "H"}


def read_a3m(path):
    seqs, cur = [], []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return ["".join(c for c in s if not c.islower()) for s in seqs]


def encode(seqs, L):
    """MSA -> entiers 0..20 (20 = gap ou acide aminé non standard)."""
    lut = np.full(256, 20, dtype=np.uint8)
    for k, a in enumerate(AA):
        lut[ord(a)] = k
    arr = np.frombuffer("".join(s.ljust(L, "-")[:L] for s in seqs).encode(), dtype=np.uint8)
    return lut[arr].reshape(len(seqs), L)


def onehot(msa):
    M, L = msa.shape
    X = np.zeros((M, L * Q), dtype=np.float32)
    X[np.arange(M)[:, None], np.arange(L)[None, :] * Q + msa] = 1.0
    return X


def weights_80(X, L, block=1000):
    """w_a = 1 / |{b : identité(a,b) >= 0.8}|, calculé par blocs (X @ X.T = # positions identiques)."""
    M = X.shape[0]
    thr = IDT * L
    n = np.zeros(M, dtype=np.int32)
    for s in range(0, M, block):
        sim = X[s:s + block] @ X.T                 # (block, M) : nombre de colonnes identiques
        n[s:s + block] = (sim >= thr).sum(axis=1)
    return 1.0 / n


def cbeta(pdb):
    """Coordonnées Cbeta (CA pour Gly) par numéro de résidu."""
    best = {}
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        name, resn, resi = line[12:16].strip(), line[17:20].strip(), int(line[22:26])
        if name == "CB" or (name == "CA" and resn == "GLY"):
            best[resi] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    idx = sorted(best)
    return idx, np.array([best[i] for i in idx])


# ---------------------------------------------------------------- MSA & poids
seqs = read_a3m(A3M)
L = len(seqs[0])
query = seqs[0]
for pos, aa in TRIAD.items():
    assert query[pos - 1] == aa, f"numérotation cassée : position {pos} = {query[pos-1]}, attendu {aa}"
print(f"MSA : {len(seqs)} séquences x {L} colonnes ; triade vérifiée sur la requête")

msa = encode(seqs, L)
X = onehot(msa)
print("Repondération à 80 % d'identité ...")
w = weights_80(X, L).astype(np.float32)
meff = float(w.sum())
print(f"  Meff = {meff:.0f} (brut {len(seqs)}) -> redondance ~{len(seqs)/meff:.1f}x ; Meff/L = {meff/L:.1f}")

gapfrac = (msa == 20).mean(axis=0)
noisy = np.where(gapfrac > 0.5)[0]
print(f"  colonnes à > 50 % de gaps : {len(noisy)}"
      + (f" (positions {[int(x) + 1 for x in noisy]}) — exclues du classement" if len(noisy) else ""))

# ------------------------------------------------------- fréquences (mfDCA)
fij = ((X * w[:, None]).T @ X).astype(np.float64) / meff        # (L*Q, L*Q)
fi = np.diag(fij).reshape(L, Q).copy()

fi_pc = (1 - PC) * fi + PC / Q
fij_pc = (1 - PC) * fij + PC / (Q * Q)
for i in range(L):                                             # blocs diagonaux : Pii(a,b) = fi(a) d_ab
    s = i * Q
    fij_pc[s:s + Q, s:s + Q] = np.diag(fi_pc[i])

keep = (np.arange(L * Q).reshape(L, Q)[:, :Q - 1]).ravel()     # on retire l'état gap (référence)
C = fij_pc[np.ix_(keep, keep)] - np.outer(fi_pc[:, :Q - 1].ravel(), fi_pc[:, :Q - 1].ravel())
print(f"Inversion de la matrice de covariance {C.shape[0]}x{C.shape[0]} ...")
e = -np.linalg.inv(C)

# ------------------------------------------ Frobenius (jauge zéro-somme) + APC
q1 = Q - 1
E = e.reshape(L, q1, L, q1)
J = E.transpose(0, 2, 1, 3)                                    # (L, L, 20, 20)
J = J - J.mean(axis=2, keepdims=True) - J.mean(axis=3, keepdims=True) + J.mean(axis=(2, 3), keepdims=True)
FN = np.sqrt((J ** 2).sum(axis=(2, 3)))
np.fill_diagonal(FN, 0.0)
mean_row = FN.sum(axis=1) / (L - 1)
apc = FN - np.outer(mean_row, mean_row) / mean_row.mean()
np.fill_diagonal(apc, 0.0)

# ------------------------------------------------------------- structure 3D
resi, cb = cbeta(PDB)
assert len(resi) == L, f"modèle {len(resi)} résidus vs MSA {L}"
D = np.linalg.norm(cb[:, None, :] - cb[None, :, :], axis=2)

# ------------------------------------------------- classement des paires
ii, jj = np.triu_indices(L, k=1)
ok = np.ones(len(ii), dtype=bool)
for c in noisy:
    ok &= (ii != c) & (jj != c)
sep = ok & (np.abs(ii - jj) >= SEP)
order = np.argsort(-apc[ii, jj])
order_sep = order[sep[order]]

# VALIDATION de l'implémentation : précision top-L/2 et top-L contre la carte de contacts
lines = []


def precision(n):
    sel = order_sep[:n]
    return float((D[ii[sel], jj[sel]] < CONTACT).mean())


base = float((D[ii[sep], jj[sep]] < CONTACT).mean())
p_half, p_full = precision(L // 2), precision(L)
lines.append(f"Meff (80 % id) = {meff:.0f} sur {len(seqs)} brutes ; Meff/L = {meff/L:.1f}")
lines.append(f"VALIDATION top-L/2 ({L//2} paires) : précision {p_half:.1%} de vrais contacts "
             f"(Cbeta < {CONTACT} A) ; top-L : {p_full:.1%} ; fond = {base:.1%}")
lines.append(f"  facteur d'enrichissement top-L/2 = {p_half/base:.1f}x")

# ------------------------------- couplages NON expliqués par le monomère
sel = order_sep[:L]
far = D[ii[sel], jj[sel]] > FAR
lines.append(f"Parmi les {L} paires les mieux couplées : {int(far.sum())} ont une distance "
             f"Cbeta > {FAR} A dans le modèle monomère (candidats interface/allostérie)")

# ----------------------------------------------- paires de la triade
lines.append("")
lines.append("PAIRES DE LA TRIADE (rang parmi les paires classées, |i-j| >= 5 sauf mention) :")
allpairs = list(zip(ii[order], jj[order]))
rank_of = {(a, b): r for r, (a, b) in enumerate(allpairs, 1)}
tri = sorted(TRIAD)
for a in range(len(tri)):
    for b in range(a + 1, len(tri)):
        p, q_ = tri[a] - 1, tri[b] - 1
        r = rank_of[(min(p, q_), max(p, q_))]
        note = "" if abs(p - q_) >= SEP else "  [local, hors filtre |i-j|>=5]"
        lines.append(f"  {TRIAD[tri[a]]}{tri[a]}-{TRIAD[tri[b]]}{tri[b]} : APC {apc[p, q_]:+.3f} "
                     f"rang {r}/{len(allpairs)}  d(Cbeta) = {D[p, q_]:.1f} A{note}")

# --------------------------- partenaires co-évoluants de chaque ligand métal
rows = [("ligand", "partenaire", "aa", "apc", "rang_global", "d_cbeta", "eloigne_>15A")]
lines.append("")
lines.append("PARTENAIRES CO-ÉVOLUANTS des ligands métal (top 8 chacun, |i-j| >= 5) :")
for pos in tri:
    p = pos - 1
    part = [(apc[p, k], k) for k in range(L) if abs(k - p) >= SEP and k not in noisy]
    part.sort(reverse=True)
    lines.append(f"  {TRIAD[pos]}{pos} :")
    for score, k in part[:8]:
        r = rank_of[(min(p, k), max(p, k))]
        d = D[p, k]
        lines.append(f"      {query[k]}{k+1:<4} APC {score:+.3f}  rang {r:>5}  d = {d:5.1f} A"
                     + ("   <- non expliqué par le monomère" if d > FAR else ""))
        rows.append((f"{TRIAD[pos]}{pos}", f"{query[k]}{k+1}", query[k], f"{score:.4f}",
                     str(r), f"{d:.1f}", "oui" if d > FAR else "non"))

(OUT / "dca_triad_partners.tsv").write_text("\n".join("\t".join(r) for r in rows) + "\n")

prs: list = [("rang", "pos_i", "aa_i", "pos_j", "aa_j", "apc", "fn", "d_cbeta", "contact_<8A")]
for r, s in enumerate(order_sep[:500], 1):
    a, b = ii[s], jj[s]
    prs.append((str(r), str(a + 1), query[a], str(b + 1), query[b], f"{apc[a, b]:.4f}",
                f"{FN[a, b]:.4f}", f"{D[a, b]:.1f}", "oui" if D[a, b] < CONTACT else "non"))
(OUT / "dca_pairs.tsv").write_text("\n".join("\t".join(r) for r in prs) + "\n")

np.savez_compressed(OUT / "dca_matrices.npz", apc=apc, fn=FN, weights=w, meff=meff,
                    gapfrac=gapfrac, dist=D, noisy=noisy,
                    query=np.array(list(query)))

txt = "\n".join(lines)
(OUT / "dca_summary.txt").write_text(txt + "\n")
print("\n" + txt)
print(f"\nÉcrit : résultats/dca_pairs.tsv (500 paires), dca_triad_partners.tsv, dca_summary.txt")
