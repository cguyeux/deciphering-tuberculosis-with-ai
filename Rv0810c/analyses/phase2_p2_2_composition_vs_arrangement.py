#!/usr/bin/env python3
"""
P2.2 (suite) — Separer COMPOSITION et ARRANGEMENT.

La phase 1 a etabli deux choses :
  - le garde-fou invoque est REEL dans ce proteome (FCR mediane 0,260 pour les
    proteines < 75 aa contre 0,206-0,224 au-dela : les courtes sont bien plus
    chargees). Le controle apparie en longueur etait donc necessaire.
  - kappa depasse 1 sur les sequences pauvres en charges : le parametre de
    Das & Pappu n'est defini que pour les polyampholytes (FCR appreciable).
    Le comparer a une population non appariee en FCR gonfle le nul.

Deux questions distinctes, donc deux nuls distincts :

  (A) COMPOSITION — Rv0810c est-elle anormalement chargee pour une proteine de
      60 aa de H37Rv ? Nul : population appariee en LONGUEUR.

  (B) ARRANGEMENT — a composition donnee, les charges sont-elles anormalement
      SEGREGEES (tete basique / queue acide) ? Nul exact : permutation des
      PROPRES residus de Rv0810c (10^5 melanges). Ce nul controle parfaitement
      la composition : meme longueur, memes acides amines, seul l'ordre change.
      Complete par une comparaison a une population appariee en longueur ET en
      FCR, ou kappa redevient interpretable.

Controle negatif de l'instrument : le meme test de permutation applique a
CHAQUE temoin apparie. Si l'instrument declarait significatif n'importe quelle
petite proteine, il ne vaudrait rien. On mesure son taux de faux positifs.

Statistiques d'arrangement : kappa (Das & Pappu 2013) et SCD (Sawle & Ghosh
2015), cette derniere sans normalisation donc insensible au probleme de kappa.

Sortie : résultats/p2_2_composition_vs_arrangement.json
"""

from __future__ import annotations

import json
import random
from collections import Counter
from pathlib import Path

import numpy as np
from Bio import SeqIO

PROJ = Path(__file__).resolve().parent.parent
PROTEOME = (
    PROJ.parent / "annotation_mtbc" / "résultats" / "phase2d_eggnog" / "proteome.faa"
)
OUTDIR = PROJ / "résultats"
OUTDIR.mkdir(exist_ok=True)

TARGET = "Rv0810c"
N_PERM = 50_000
SEED = 20260810

POS, NEG = set("KR"), set("DE")


def net_charge(s: str) -> int:
    c = Counter(s)
    return sum(c[a] for a in POS) - sum(c[a] for a in NEG)


def fcr(s: str) -> float:
    c = Counter(s)
    return (sum(c[a] for a in POS) + sum(c[a] for a in NEG)) / len(s) if s else 0.0


def ncpr(s: str) -> float:
    return net_charge(s) / len(s) if s else 0.0


def sigma(s: str) -> float:
    f = fcr(s)
    return 0.0 if f == 0 else ncpr(s) ** 2 / f


def _delta(s: str, blob: int) -> float:
    n = len(s)
    if n < blob:
        return 0.0
    g = sigma(s)
    return float(np.mean([(sigma(s[i : i + blob]) - g) ** 2 for i in range(n - blob + 1)]))


def _delta_max(s: str, blob: int) -> float:
    c = Counter(s)
    npos, nneg = sum(c[a] for a in POS), sum(c[a] for a in NEG)
    worst = "D" * nneg + "K" * npos + "G" * (len(s) - npos - nneg)
    return _delta(worst, blob)


def kappa(s: str) -> float:
    c = Counter(s)
    if sum(c[a] for a in POS) == 0 or sum(c[a] for a in NEG) == 0:
        return float("nan")
    vals = []
    for blob in (5, 6):
        dm = _delta_max(s, blob)
        if dm > 0:
            vals.append(_delta(s, blob) / dm)
    return float(np.mean(vals)) if vals else float("nan")


def scd(s: str) -> float:
    """Sequence Charge Decoration (Sawle & Ghosh 2015, J Chem Phys 143:085101).
    SCD = (1/N) * sum_{m>n} q_m q_n * sqrt(m-n).
    Plus SCD est NEGATIF, plus les charges opposees sont segregees en blocs.
    Sans normalisation : pas de pathologie a faible FCR (contrairement a kappa).
    """
    q = np.array([1.0 if a in POS else (-1.0 if a in NEG else 0.0) for a in s])
    idx = np.nonzero(q)[0]
    if len(idx) < 2:
        return 0.0
    tot = 0.0
    for a in range(1, len(idx)):
        m = idx[a]
        n = idx[:a]
        tot += float(np.sum(q[m] * q[n] * np.sqrt(m - n)))
    return tot / len(s)


def best_split_contrast(s: str, min_seg: int = 10) -> tuple[float, int]:
    n = len(s)
    best, bk = -np.inf, -1
    for k in range(min_seg, n - min_seg + 1):
        v = ncpr(s[:k]) - ncpr(s[k:])
        if v > best:
            best, bk = v, k
    return float(best), bk


def perm_test(seq: str, n_perm: int, rng: random.Random) -> dict:
    """Nul exact par permutation des propres residus : composition figee."""
    chars = list(seq)
    obs_k, obs_s, obs_c = kappa(seq), scd(seq), best_split_contrast(seq)[0]
    ge_k = ge_c = le_s = 0
    ks, ss, cs = [], [], []
    for _ in range(n_perm):
        rng.shuffle(chars)
        p = "".join(chars)
        k, sc, c = kappa(p), scd(p), best_split_contrast(p)[0]
        ks.append(k)
        ss.append(sc)
        cs.append(c)
        ge_k += k >= obs_k
        le_s += sc <= obs_s
        ge_c += c >= obs_c
    ks, ss, cs = np.array(ks), np.array(ss), np.array(cs)
    return {
        "n_permutations": n_perm,
        "kappa": {
            "observe": round(obs_k, 4),
            "nul_median": round(float(np.median(ks)), 4),
            "nul_q95": round(float(np.percentile(ks, 95)), 4),
            "p": round((ge_k + 1) / (n_perm + 1), 6),
            "z": round(float((obs_k - ks.mean()) / ks.std()), 2),
        },
        "SCD": {
            "observe": round(obs_s, 4),
            "nul_median": round(float(np.median(ss)), 4),
            "nul_q05": round(float(np.percentile(ss, 5)), 4),
            "p": round((le_s + 1) / (n_perm + 1), 6),
            "z": round(float((obs_s - ss.mean()) / ss.std()), 2),
        },
        "contraste_max_N_moins_C": {
            "observe": round(obs_c, 4),
            "nul_median": round(float(np.median(cs)), 4),
            "nul_q95": round(float(np.percentile(cs, 95)), 4),
            "p": round((ge_c + 1) / (n_perm + 1), 6),
            "z": round(float((obs_c - cs.mean()) / cs.std()), 2),
        },
    }


def rank_in(value: float, pop: list[float], tail: str) -> dict:
    a = np.array([v for v in pop if not np.isnan(v)])
    k = int(np.sum(a >= value)) if tail == "greater" else int(np.sum(a <= value))
    return {
        "valeur": round(float(value), 4),
        "n_temoins": len(a),
        "rang": k,
        "p_empirique": round(k / len(a), 4),
        "mediane_temoins": round(float(np.median(a)), 4),
        "q95_temoins": round(float(np.percentile(a, 95)), 4),
        "q05_temoins": round(float(np.percentile(a, 5)), 4),
    }


def main() -> None:
    rng = random.Random(SEED)
    prots = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(PROTEOME, "fasta")}
    tgt = prots[TARGET]

    rep: dict = {
        "piste": "P2.2 (phase 2)",
        "cible": {
            "rv": TARGET,
            "len": len(tgt),
            "FCR": round(fcr(tgt), 4),
            "NCPR": round(ncpr(tgt), 4),
            "charge_nette": net_charge(tgt),
            "kappa": round(kappa(tgt), 4),
            "SCD": round(scd(tgt), 4),
        },
    }

    # ------------------------------------------------------------------ #
    # (A) COMPOSITION : appariement en LONGUEUR seul
    # ------------------------------------------------------------------ #
    ctrlL = {i: s for i, s in prots.items() if 50 <= len(s) <= 70}
    rep["A_composition"] = {
        "question": "Rv0810c est-elle anormalement chargee pour 60 aa de H37Rv ?",
        "n_temoins_50_70aa": len(ctrlL),
        "FCR": rank_in(fcr(tgt), [fcr(s) for s in ctrlL.values()], "greater"),
        "charge_nette_absolue": rank_in(
            abs(net_charge(tgt)),
            [abs(net_charge(s)) for s in ctrlL.values()],
            "greater",
        ),
    }

    # ------------------------------------------------------------------ #
    # (B1) ARRANGEMENT : nul exact par permutation intra-sequence
    # ------------------------------------------------------------------ #
    rep["B1_arrangement_permutation_intra_sequence"] = {
        "question": (
            "A composition RIGOUREUSEMENT identique (memes 60 acides amines, "
            "seul l'ordre change), les charges de Rv0810c sont-elles anormalement "
            "segregees ?"
        ),
        "resultat": perm_test(tgt, N_PERM, rng),
    }

    # ------------------------------------------------------------------ #
    # (B2) ARRANGEMENT : population appariee en LONGUEUR **ET** en FCR
    #      (seul domaine ou kappa est interpretable)
    # ------------------------------------------------------------------ #
    f_t = fcr(tgt)
    ctrlLF = {
        i: s
        for i, s in prots.items()
        if 45 <= len(s) <= 90 and abs(fcr(s) - f_t) <= 0.08
    }
    rep["B2_arrangement_temoins_apparies_longueur_ET_FCR"] = {
        "critere": f"45-90 aa et |FCR - {f_t:.3f}| <= 0.08",
        "n_temoins": len(ctrlLF),
        "kappa": rank_in(kappa(tgt), [kappa(s) for s in ctrlLF.values()], "greater"),
        "SCD": rank_in(scd(tgt), [scd(s) for s in ctrlLF.values()], "less"),
        "contraste_max": rank_in(
            best_split_contrast(tgt)[0],
            [best_split_contrast(s)[0] for s in ctrlLF.values()],
            "greater",
        ),
        "liste_temoins": sorted(ctrlLF),
    }

    # ------------------------------------------------------------------ #
    # (C) CONTROLE NEGATIF DE L'INSTRUMENT
    #     Le test de permutation est-il trivialement significatif ?
    # ------------------------------------------------------------------ #
    rng2 = random.Random(SEED + 1)
    sub = sorted(ctrlL)  # les 75 temoins apparies en longueur
    pvals = []
    for i in sub:
        r = perm_test(prots[i], 1000, rng2)
        pvals.append((i, r["kappa"]["p"], r["SCD"]["p"], r["contraste_max_N_moins_C"]["p"]))
    pk = np.array([p[1] for p in pvals])
    ps = np.array([p[2] for p in pvals])
    pc = np.array([p[3] for p in pvals])
    rep["C_controle_negatif_instrument"] = {
        "commentaire": (
            "Meme test de permutation (1000 melanges) applique aux 75 temoins "
            "apparies en longueur. Si l'instrument declarait significative "
            "n'importe quelle petite proteine, il ne vaudrait rien."
        ),
        "n_temoins_testes": len(sub),
        "taux_p<0.05_kappa": round(float(np.mean(pk < 0.05)), 3),
        "taux_p<0.05_SCD": round(float(np.mean(ps < 0.05)), 3),
        "taux_p<0.05_contraste": round(float(np.mean(pc < 0.05)), 3),
        "temoins_significatifs_sur_les_trois": [
            p[0] for p in pvals if p[1] < 0.05 and p[2] < 0.05 and p[3] < 0.05
        ],
    }

    out = OUTDIR / "p2_2_composition_vs_arrangement.json"
    out.write_text(json.dumps(rep, indent=2, ensure_ascii=False))
    print(json.dumps(rep, indent=2, ensure_ascii=False))
    print(f"\n[ecrit] {out}")


if __name__ == "__main__":
    main()
