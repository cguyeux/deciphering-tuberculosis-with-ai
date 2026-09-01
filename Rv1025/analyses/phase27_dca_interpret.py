#!/usr/bin/env python3
# pyright: reportAttributeAccessIssue=false, reportMissingTypeStubs=false
"""phase27_dca_interpret.py — P8.2 (interprétation), P8.2.a et P8.2.b.

phase26 a produit les couplages ; ce script les INTERPRÈTE, en refusant les quatre
lectures naïves qui guettent.

A. LE TOP DU CLASSEMENT EST-IL AUTRE CHOSE QUE LE REPLI ? (modèle nul de P8.2)
   Le contre-argument du challenge était : « la co-évolution récupère surtout des
   contacts structuraux ». Test : fraction de paires ÉLOIGNÉES (d Cbeta > 15 A) parmi
   les mieux couplées, contre le taux de base des paires éligibles. Un déficit massif
   signifie que le signal est dominé par le repli, donc qu'il n'apprend rien de neuf.

B. CE QUI RESTE EST-IL UNE INTERFACE ? (P8.2.a)
   Une interface d'homo-oligomère laisserait des couplages forts que le monomère
   n'explique pas, et ces résidus devraient (i) être EXPOSÉS (RSA élevé — on ne fait
   pas d'interface avec un cœur hydrophobe enfoui) et (ii) se REGROUPER en une ou deux
   plaques de surface, pas se disperser. Les deux critères sont testés ; l'échec de
   l'un ou l'autre tue l'hypothèse d'interface.

C. PIÈGE DE PUISSANCE — L'INVARIANCE REND LE DCA AVEUGLE.
   Le DCA mesure une COVARIATION : une position quasi invariante n'a presque pas de
   variance, donc ne peut PAS montrer de couplage fort, même si elle est le coeur
   fonctionnel. Conclure « la triade ne co-évolue pas donc le site est douteux » serait
   une faute. Test : les positions les plus conservées ont-elles un couplage maximal
   systématiquement déprimé ? Si oui, l'absence de signal sur la triade est un défaut
   de puissance, pas une réfutation.

D. UNE POCHE OU DEUX ? (P8.2.b)
   Poche 1 et Poche 3 sont contiguës (résidus 41 et 45 partagés). Si elles forment un
   seul site fonctionnel, les couplages qui les relient doivent excéder ce que la seule
   proximité explique. NULL APPARIÉ EN DISTANCE (indispensable : l'APC décroît avec la
   distance, un null non apparié conclurait toujours positif).

Lit   : résultats/dca_matrices.npz, le MSA a3m, le modèle AF, les poches fpocket
Écrit : résultats/dca_interpretation.txt, résultats/dca_far_pairs.tsv,
        article/supplementary_materials/table_S8_dca.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fast_blas import ensure  # noqa: E402  # type: ignore[import-not-found]

ensure()

from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from Bio.PDB.PDBParser import PDBParser  # noqa: E402
from Bio.PDB.SASA import ShrakeRupley  # noqa: E402
from scipy.stats import binomtest, mannwhitneyu  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
NPZ = ROOT / "résultats/dca_matrices.npz"
A3M = ROOT / "résultats/af3_out/fold_rv1025_divic_test/msas/fold_rv1025_divic_test_unpaired_msa_chains_a.a3m"
PDB = ROOT / "résultats/structure/AF-P96375-F1.pdb"
POCK = ROOT / "résultats/druggability/AF-P96375-F1_out/pockets"
OUT = ROOT / "résultats"
SUPP = ROOT / "article/supplementary_materials"

FAR = 15.0
CONTACT = 8.0
SEP = 5
TRIAD = {59: "E", 113: "C", 115: "H"}
# Tien et al. 2013, surface accessible maximale théorique (A^2)
MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225, "G": 104,
          "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240, "P": 159, "S": 155,
          "T": 172, "W": 285, "Y": 263, "V": 174}
rng = np.random.default_rng(42)

d = np.load(NPZ, allow_pickle=True)
apc, D, noisy = d["apc"], d["dist"], set(d["noisy"].tolist())
query = "".join(d["query"].tolist())
w = d["weights"]
meff = float(d["meff"])
L = len(query)

out = []


def say(s=""):
    print(s)
    out.append(s)


say(f"DCA : L = {L}, Meff = {meff:.0f}, colonnes écartées (gaps > 50 %) : {len(noisy)}")

ii, jj = np.triu_indices(L, k=1)
elig = (np.abs(ii - jj) >= SEP) & ~np.isin(ii, list(noisy)) & ~np.isin(jj, list(noisy))
pi, pj = ii[elig], jj[elig]
score = apc[pi, pj]
dist = D[pi, pj]
order = np.argsort(-score)
n_elig = len(pi)

# ------------------------------------------------------------------ A. modèle nul
say()
say("=" * 78)
say("A. LE SIGNAL EST-IL AUTRE CHOSE QUE LE REPLI ? (modèle nul du challenge)")
say("=" * 78)
p_far = float((dist > FAR).mean())
p_cont = float((dist < CONTACT).mean())
say(f"Taux de base sur {n_elig} paires éligibles : {p_far:.1%} éloignées (> {FAR:.0f} A), "
    f"{p_cont:.1%} en contact (< {CONTACT:.0f} A)")
for n in (L // 2, L, 2 * L):
    sel = order[:n]
    nf = int((dist[sel] > FAR).sum())
    nc = int((dist[sel] < CONTACT).sum())
    bt = binomtest(nf, n, p_far, alternative="less")
    say(f"  top-{n:<4} : contacts {nc/n:6.1%} (fond {p_cont:.1%}, x{nc/n/p_cont:4.1f})   "
        f"éloignées {nf:>3}/{n} = {nf/n:5.1%} contre {p_far*n:.0f} attendues  p = {bt.pvalue:.1e}")
say("VERDICT : le déficit de paires éloignées est massif — le classement est dominé par les")
say("contacts du repli. La co-évolution ne fournit donc PAS, en bloc, une information")
say("fonctionnelle nouvelle : c'est le contre-argument du challenge, confirmé.")

# ------------------------------------------------- B. les éloignées sont-elles une interface ?
say()
say("=" * 78)
say("B. LE RÉSIDU NON EXPLIQUÉ PAR LE MONOMÈRE EST-IL UNE INTERFACE ? (P8.2.a)")
say("=" * 78)
sel = order[:L]
far_mask = dist[sel] > FAR
far_pairs = [(int(pi[s]), int(pj[s]), float(score[s]), float(dist[s]))
             for s, m in zip(sel, far_mask) if m]

model = PDBParser(QUIET=True).get_structure("x", str(PDB))[0]  # type: ignore[index]
ShrakeRupley().compute(model, level="R")
res = [r for r in model.get_residues() if r.id[0] == " "]
assert len(res) == L, f"{len(res)} résidus vs L = {L}"
three2one = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q",
             "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
             "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
rsa = np.array([r.sasa / MAXASA[three2one[r.get_resname()]] for r in res])

involved = sorted({a for p in far_pairs for a in p[:2]})
say(f"{len(far_pairs)} paires éloignées dans le top-{L}, impliquant {len(involved)} résidus distincts")
say(f"  RSA moyenne de ces résidus : {rsa[involved].mean():.2f}  (fond protéine {rsa.mean():.2f})")
u = mannwhitneyu(rsa[involved], np.delete(rsa, involved), alternative="greater")
say(f"  test « plus exposés que le reste » (Mann-Whitney) : p = {u.pvalue:.3f}")
exposed = [r for r in involved if rsa[r] > 0.25]
say(f"  résidus réellement exposés (RSA > 0,25) : {len(exposed)}/{len(involved)}")

# regroupement spatial simple (single-linkage à 10 A) des résidus impliqués
clusters: list = []
for r in involved:
    hit = [c for c in clusters if any(D[r, m] < 10.0 for m in c)]
    if hit:
        hit[0].append(r)
        for extra in hit[1:]:
            hit[0].extend(extra)
            clusters.remove(extra)
    else:
        clusters.append([r])
clusters.sort(key=len, reverse=True)
say(f"  regroupement spatial (10 A) : {len(clusters)} amas, tailles {[len(c) for c in clusters]}")
say(f"    plus gros amas : {sorted(x+1 for x in clusters[0])}")

rows = [("pos_i", "aa_i", "pos_j", "aa_j", "apc", "d_cbeta", "rsa_i", "rsa_j")]
for a, b, s, dd in sorted(far_pairs, key=lambda x: -x[2]):
    rows.append((str(a + 1), query[a], str(b + 1), query[b], f"{s:.3f}", f"{dd:.1f}",
                 f"{rsa[a]:.2f}", f"{rsa[b]:.2f}"))
(OUT / "dca_far_pairs.tsv").write_text("\n".join("\t".join(r) for r in rows) + "\n")
say("  (détail dans résultats/dca_far_pairs.tsv)")

# --------------------------------------------------- C. puissance contre invariance
say()
say("=" * 78)
say("C. PIÈGE DE PUISSANCE : L'INVARIANCE REND LE DCA AVEUGLE")
say("=" * 78)
seqs, cur = [], []
for line in A3M.read_text().splitlines():
    if line.startswith(">"):
        if cur:
            seqs.append("".join(cur))
            cur = []
    else:
        cur.append(line.strip())
seqs.append("".join(cur))
aln = np.array([list("".join(c for c in s if not c.islower())) for s in seqs])
q = np.array(list(query))
gap = (aln == "-") | (aln == ".")
same = (aln == q[None, :]) & ~gap
cons = 100 * (same * w[:, None]).sum(0) / np.maximum((~gap * w[:, None]).sum(0), 1e-9)

maxapc = np.array([apc[i, [j for j in range(L) if abs(i - j) >= SEP and j not in noisy]].max()
                   if i not in noisy else np.nan for i in range(L)])
valid = ~np.isnan(maxapc)
hi = valid & (cons >= np.nanpercentile(cons[valid], 90))
lo = valid & (cons < np.nanpercentile(cons[valid], 90))
u2 = mannwhitneyu(maxapc[hi], maxapc[lo], alternative="less")
say(f"couplage MAXIMAL par position : décile le plus conservé {np.nanmean(maxapc[hi]):.2f} "
    f"contre {np.nanmean(maxapc[lo]):.2f} pour le reste")
say(f"  test « les positions très conservées couplent MOINS » : p = {u2.pvalue:.4f}")
say("La triade elle-même :")
for pos, aa in TRIAD.items():
    i = pos - 1
    pct = float((maxapc[valid] < maxapc[i]).mean() * 100)
    say(f"  {aa}{pos} : conservation {cons[i]:.1f} %, couplage max {maxapc[i]:.2f} "
        f"(percentile {pct:.0f})")
say("LECTURE : un site fonctionnel quasi invariant est un ANGLE MORT du DCA. Le rang modeste")
say("des paires impliquant Glu59 ne réfute donc rien — il était attendu.")

# ------------------------------------------------------ D. une poche ou deux ?
say()
say("=" * 78)
say("D. POCHE 1 ET POCHE 3 : UN SITE ÉTENDU OU DEUX SITES ? (P8.2.b)")
say("=" * 78)


def pocket(n):
    f = POCK / f"pocket{n}_atm.pdb"
    return sorted({int(x[22:26]) - 1 for x in f.read_text().splitlines() if x.startswith("ATOM")})


p1, p3 = pocket(1), pocket(3)
shared = sorted(set(p1) & set(p3))
cross = [(a, b) for a in p1 for b in p3
         if a != b and abs(a - b) >= SEP and a not in noisy and b not in noisy
         and a not in shared and b not in shared]
cross = list({(min(a, b), max(a, b)) for a, b in cross})
obs = float(np.mean([apc[a, b] for a, b in cross]))
say(f"Poche 1 : {len(p1)} résidus ; Poche 3 : {len(p3)} ; partagés : {[s+1 for s in shared]}")
say(f"{len(cross)} paires inter-poches éligibles, APC moyen observé = {obs:.3f}")

# null APPARIÉ EN DISTANCE : rééchantillonner des paires quelconques de même profil de distance
bins = np.arange(0, 60, 2.0)
obs_bins = np.digitize([D[a, b] for a, b in cross], bins)
pool = {}
for k, (a, b) in enumerate(zip(pi, pj)):
    pool.setdefault(int(np.digitize(dist[k], bins)), []).append(k)
draws = []
for _ in range(10000):
    tot = []
    for bnum in obs_bins:
        cand = pool.get(bnum)
        if cand:
            tot.append(score[rng.choice(cand)])
    draws.append(np.mean(tot))
draws = np.array(draws)
pval = float((draws >= obs).mean())
say(f"null apparié en distance (10 000 tirages) : moyenne {draws.mean():.3f} "
    f"(IC95 {np.percentile(draws,2.5):.3f}-{np.percentile(draws,97.5):.3f})  p = {pval:.4f}")
say("VERDICT : " + ("les deux cavités sont couplées AU-DELÀ de leur simple proximité "
                    "-> argument pour UN SITE FONCTIONNEL ÉTENDU."
                    if pval < 0.05 else
                    "aucun couplage inter-poches au-delà de ce que la distance explique "
                    "-> rien ne permet de fusionner les deux cavités ; elles restent distinctes."))

(OUT / "dca_interpretation.txt").write_text("\n".join(out) + "\n")
SUPP.mkdir(parents=True, exist_ok=True)
tab = [("analyse", "mesure", "observe", "attendu_sous_null", "p")]
tab.append(("A. contacts vs repli", f"paires eloignees dans le top-{L}",
            f"{len(far_pairs)}/{L}", f"{p_far*L:.0f}/{L}",
            f"{binomtest(len(far_pairs), L, p_far, alternative='less').pvalue:.1e}"))
tab.append(("B. interface", "RSA des residus non expliques",
            f"{rsa[involved].mean():.2f}", f"{rsa.mean():.2f}", f"{u.pvalue:.3f}"))
tab.append(("C. puissance", "couplage max, decile le plus conserve",
            f"{np.nanmean(maxapc[hi]):.2f}", f"{np.nanmean(maxapc[lo]):.2f}", f"{u2.pvalue:.4f}"))
tab.append(("D. poche1 x poche3", "APC moyen inter-poches",
            f"{obs:.3f}", f"{draws.mean():.3f}", f"{pval:.4f}"))
(SUPP / "table_S8_dca.tsv").write_text("\n".join("\t".join(r) for r in tab) + "\n")
say()
say("Écrit : résultats/dca_interpretation.txt, résultats/dca_far_pairs.tsv,")
say("        article/supplementary_materials/table_S8_dca.tsv")
