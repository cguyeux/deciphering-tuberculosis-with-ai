#!/usr/bin/env python3
"""phase34_eukaryote_origin.py — P8.1.a.2 : d'où vient la branche eucaryote de DUF501 ?

Deux origines possibles pour les 102 membres eucaryotes de PF04417 (phase29) :
  (a) acquisition MITOCHONDRIALE, l'endosymbiote fondateur étant une alphaprotéobactérie
      (classe de Pseudomonadota) ;
  (b) transfert horizontal depuis les ACTINOBACTÉRIES (le phylum qui porte 87 % de la
      famille bactérienne, cf. P6.6).

TEST (déjà outillé par la piste) : comparer, pour chaque séquence eucaryote de
l'alignement Pfam full (data/PF04417_full.sto, colonnes déjà alignées), son identité au
plus proche homologue ACTINOBACTÉRIEN contre son identité au plus proche homologue
PSEUDOMONADOTA. Une origine mitochondriale prédit une proximité systématique aux
Pseudomonadota ; un transfert actinobactérien prédit l'inverse.

GARDE-FOU (à énoncer avant de calculer, pas après) : le meilleur-hit (max) sur un groupe
est mécaniquement tiré vers le haut par la PROFONDEUR du groupe échantillonné —
Actinomycetota (1645 séq.) vs Pseudomonadota (104 séq.) dans cet alignement, un facteur
~16x. Comparer les deux max bruts confondrait donc "plus proche parent" avec "groupe plus
échantillonné". PARADE : sous-échantillonner Actinomycetota à N=104 (taille de
Pseudomonadota), répéter 1000 fois, et comparer le max PSEUDOMONADOTA à la distribution
NULLE des max Actinomycetota-sous-échantillonnés de même taille — comparaison à effectif
égal. Le résultat brut (groupes complets) est rapporté à titre descriptif uniquement.

NULL FALSIFIANT : si les deux identités (à taille égale) sont statistiquement
indiscernables, l'alignement ne tranche pas l'origine et il faut une vraie phylogénie
(plus coûteuse) — conclusion consignée telle quelle, pas habillée en résultat positif.

Lit   : data/PF04417_full.sto, résultats/duf501_family_triad.tsv (groupe bactérien déjà
        assigné), article/supplementary_materials/table_S10_duf501_eukaryotes.tsv (102
        membres eucaryotes)
Écrit : résultats/phase34_eukaryote_origin.tsv, résultats/phase34_eukaryote_origin.txt
"""
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fast_blas import ensure  # noqa: E402  # type: ignore[import-not-found]

ensure()

from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.stats import mannwhitneyu, wilcoxon  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
STO = ROOT / "data/PF04417_full.sto"
TRIAD_TSV = ROOT / "résultats/duf501_family_triad.tsv"
S10 = ROOT / "article/supplementary_materials/table_S10_duf501_eukaryotes.tsv"
OUT = ROOT / "résultats"
N_BOOT = 1000
rng = np.random.default_rng(42)

out = []


def say(s=""):
    print(s)
    out.append(s)


# ------------------------------------------------------------------ lecture alignement
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
assert all(len(s) == L for s in seqs.values())
say(f"Alignement PF04417 full : {len(order)} séquences x {L} colonnes")

# ------------------------------------------------------------------ groupes
groupe = {}
with open(TRIAD_TSV, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        groupe[row["sequence"]] = row["groupe"]

euk = []
with open(S10, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        euk.append(row["sequence"])
euk_set = set(euk)
say(f"Membres eucaryotes (table S10) : {len(euk)}")

acti = [n for n in order if n not in euk_set and groupe.get(n) == "Actinomycetota"]
pseudo = [n for n in order if n not in euk_set and groupe.get(n) == "Pseudomonadota"]
say(f"Référence bactérienne : Actinomycetota n={len(acti)}, Pseudomonadota n={len(pseudo)}"
    f" (ratio {len(acti)/len(pseudo):.1f}x — motive le contrôle à taille égale ci-dessous)")

# ------------------------------------------------------------------ matrice encodée
AA_IDX = {c: i for i, c in enumerate("ACDEFGHIKLMNPQRSTVWY")}


def encode(name):
    s = seqs[name]
    arr = np.full(L, -1, dtype=np.int8)
    for i, ch in enumerate(s):
        u = ch.upper()
        if u in AA_IDX:
            arr[i] = AA_IDX[u]
    return arr


names_ref = acti + pseudo
enc_ref = np.stack([encode(n) for n in names_ref])          # (Nref, L)
enc_euk = np.stack([encode(n) for n in euk])                 # (Neuk, L)
is_acti = np.array([1] * len(acti) + [0] * len(pseudo), dtype=bool)

valid_ref = enc_ref >= 0
valid_euk = enc_euk >= 0


def identity_matrix(euk_arr, euk_valid, ref_arr, ref_valid):
    """% identité par paire sur les colonnes non-gap dans LES DEUX séquences."""
    both_valid = euk_valid[:, None, :] & ref_valid[None, :, :]        # (Ne,Nr,L)
    n_cols = both_valid.sum(axis=2)
    match = (euk_arr[:, None, :] == ref_arr[None, :, :]) & both_valid
    n_match = match.sum(axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        pid = np.where(n_cols > 0, n_match / np.maximum(n_cols, 1) * 100.0, 0.0)
    return pid, n_cols


say("Calcul des identités par paire (colonnes co-alignées, non-gap des deux côtés)...")
pid, ncols = identity_matrix(enc_euk, valid_euk, enc_ref, valid_ref)   # (Neuk, Nref)

best_acti_full = pid[:, is_acti].max(axis=1)
best_pseudo_full = pid[:, ~is_acti].max(axis=1)

say()
say("=" * 78)
say("A. COMPARAISON BRUTE (groupes complets, PAS à taille égale — descriptif seulement)")
say("=" * 78)
say(f"  meilleur-hit Actinomycetota (n={len(acti)}) : médiane {np.median(best_acti_full):.1f}%,"
    f" moyenne {best_acti_full.mean():.1f}%")
say(f"  meilleur-hit Pseudomonadota (n={len(pseudo)}) : médiane {np.median(best_pseudo_full):.1f}%,"
    f" moyenne {best_pseudo_full.mean():.1f}%")
n_acti_wins_raw = int((best_acti_full > best_pseudo_full).sum())
say(f"  Actinomycetota gagne (brut) pour {n_acti_wins_raw}/{len(euk)} séquences eucaryotes")

# ------------------------------------------------------------------ B. contrôle à taille égale
say()
say("=" * 78)
say(f"B. TEST DÉCISIF : sous-échantillonnage Actinomycetota à N={len(pseudo)} "
    f"(taille Pseudomonadota), {N_BOOT} tirages")
say("=" * 78)
acti_idx = np.where(is_acti)[0]
pid_acti_only = pid[:, acti_idx]              # (Neuk, Nacti)
boot_best_acti = np.empty((N_BOOT, len(euk)))
for b in range(N_BOOT):
    sub = rng.choice(len(acti_idx), size=len(pseudo), replace=False)
    boot_best_acti[b] = pid_acti_only[:, sub].max(axis=1)

acti_matched_mean_per_seq = boot_best_acti.mean(axis=0)   # (Neuk,) moyenne sur les 1000 tirages
say(f"  meilleur-hit Actinomycetota (sous-échantillon N={len(pseudo)}, moyenne sur {N_BOOT} tirages) :"
    f" médiane {np.median(acti_matched_mean_per_seq):.1f}%, moyenne {acti_matched_mean_per_seq.mean():.1f}%")
say(f"  meilleur-hit Pseudomonadota (groupe complet, N={len(pseudo)}) :"
    f" médiane {np.median(best_pseudo_full):.1f}%, moyenne {best_pseudo_full.mean():.1f}%")

# test apparié par séquence eucaryote : Wilcoxon signé sur (acti_matched - pseudo)
diff = acti_matched_mean_per_seq - best_pseudo_full
w_stat, w_p = wilcoxon(diff)
n_acti_wins_matched = int((diff > 0).sum())
say(f"  Wilcoxon signé apparié (acti sous-échant. vs pseudo, par séquence eucaryote) :"
    f" statistique={w_stat:.1f}, p={w_p:.3g}")
say(f"  Actinomycetota (taille égale) gagne pour {n_acti_wins_matched}/{len(euk)} séquences,"
    f" médiane de l'écart = {np.median(diff):+.1f} points")

# empirical p-value par séquence : fraction des tirages où le sous-échantillon actino
# dépasse le pseudomonadota observé (permutation directe, sans hypothèse de distribution)
emp_p = (boot_best_acti >= best_pseudo_full[None, :]).mean(axis=0)
say(f"  p empirique médian (fraction des tirages actino >= pseudo observé) : {np.median(emp_p):.3g}")

u_stat, u_p = mannwhitneyu(acti_matched_mean_per_seq, best_pseudo_full, alternative="two-sided")
say(f"  Mann-Whitney (distributions, non apparié, contrôle) : U={u_stat:.1f}, p={u_p:.3g}")

# ------------------------------------------------------------------ C. cas du rumen (P8.1.a.3, contexte)
say()
say("=" * 78)
say("C. CONTEXTE — les 3 séquences à plus haute identité bactérienne (candidats HGT récents)")
say("=" * 78)
order_idx = np.argsort(-best_acti_full)
for i in order_idx[:5]:
    say(f"  {euk[i]:<30} acti_max={best_acti_full[i]:5.1f}%  pseudo_max={best_pseudo_full[i]:5.1f}%"
        f"  acti_matched={acti_matched_mean_per_seq[i]:5.1f}%")

# ------------------------------------------------------------------ D. verdict
say()
say("=" * 78)
say("D. VERDICT")
say("=" * 78)
if w_p < 0.05 and np.median(diff) > 0:
    say("  Actinomycetota SIGNIFICATIVEMENT plus proche que Pseudomonadota, À TAILLE ÉGALE."
        " Compatible avec un transfert horizontal depuis les Actinobactéries (ou une origine"
        " ancestrale commune que l'alignement seul ne peut distinguer d'un HGT).")
elif w_p < 0.05 and np.median(diff) < 0:
    say("  Pseudomonadota SIGNIFICATIVEMENT plus proche que Actinomycetota, À TAILLE ÉGALE."
        " Compatible avec une acquisition mitochondriale (endosymbiote alphaprotéobactérien).")
else:
    say("  AUCUNE différence significative à taille égale (NULL NON REJETÉ) : l'identité de"
        " séquence seule NE TRANCHE PAS l'origine. Une phylogénie (positionnement des séquences"
        " eucaryotes dans un arbre ML avec les deux groupes bactériens) serait nécessaire pour"
        " aller plus loin — non entreprise ici (coût plus élevé, gain incertain).")
say()
say("LIMITE MÉTHODOLOGIQUE (à ne pas dissimuler) : 'Pseudomonadota' ici est le groupe COMPLET du"
    " phylum tel que bucketé par le code mnémonique UniProt (phase28), PAS restreint aux"
    " Alphaproteobacteria — le clade spécifiquement pertinent pour une origine mitochondriale."
    " Un résultat 'Pseudomonadota gagne' devrait donc être creusé au niveau de la classe avant"
    " d'être lu comme 'mitochondrial' ; ce script ne fait pas cette résolution plus fine.")

# ------------------------------------------------------------------ écriture
rows = [("sequence", "best_acti_full_%", "best_pseudo_full_%",
         "best_acti_matched_N104_mean_%", "diff_matched_minus_pseudo", "p_empirique")]
for i, n in enumerate(euk):
    rows.append((n, f"{best_acti_full[i]:.2f}", f"{best_pseudo_full[i]:.2f}",
                 f"{acti_matched_mean_per_seq[i]:.2f}", f"{diff[i]:+.2f}", f"{emp_p[i]:.4g}"))
(OUT / "phase34_eukaryote_origin.tsv").write_text("\n".join("\t".join(r) for r in rows) + "\n")
(OUT / "phase34_eukaryote_origin.txt").write_text("\n".join(out) + "\n")
say()
say("Écrit : résultats/phase34_eukaryote_origin.tsv, résultats/phase34_eukaryote_origin.txt")
