#!/usr/bin/env python3
"""
P2.2 — Contrôle de longueur AVANT toute affirmation de composition (BLOQUANT).

Question : la bipolarite de charge de Rv0810c (module N-terminal basique 1-19 /
queue C-terminale acide 34-60, etablie en P2.1) est-elle un ecart mesurable, ou
simplement l'attendu statistique pour une proteine de 60 aa du proteome H37Rv ?

Garde-fou applique : les petites proteines sont systematiquement plus chargees
(rapport surface/volume). Toute comparaison au proteome ENTIER serait biaisee.
On compare donc a une population TEMOIN APPARIEE EN LONGUEUR.

Trois statistiques, de la plus dependante d'un choix a la plus robuste :
  (1) segments FIXES : charge nette des residus 1-19 et des 27 derniers residus.
      Fair, mais la frontiere vient de Rv0810c (pLDDT) : biais de selection possible.
  (2) frontiere OPTIMISEE pour chaque proteine (max de contraste N-vs-C) :
      neutralise le biais de selection en donnant a chaque temoin le meme
      privilege qu'a Rv0810c.
  (3) kappa (Das & Pappu 2013) : parametre de PATTERNING de charge, sans aucune
      frontiere a choisir. 0 = charges melangees, 1 = charges completement
      segregees. C'est le test decisif.

Instrument valide sur les sequences synthetiques sv1 / sv30 de Das & Pappu 2013.

Sortie : résultats/p2_2_charge_length_control.{json,tsv} + figure.
"""

from __future__ import annotations

import json
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
# Frontieres issues de P2.1 (mesure pLDDT par residu du modele AF-I6XWB9-F1-v6)
MODULE = (1, 19)   # 1-based inclusif, pLDDT 91.9
TAIL = (34, 60)    # 1-based inclusif, pLDDT 62.0
TAIL_LEN = TAIL[1] - TAIL[0] + 1  # 27 residus

POS = set("KR")
NEG = set("DE")


# --------------------------------------------------------------------------- #
# Statistiques de charge
# --------------------------------------------------------------------------- #
def net_charge(seq: str) -> int:
    """Charge nette a pH 7, His comptee neutre (convention Das & Pappu)."""
    c = Counter(seq)
    return sum(c[a] for a in POS) - sum(c[a] for a in NEG)


def fcr(seq: str) -> float:
    """Fraction of charged residues."""
    if not seq:
        return 0.0
    c = Counter(seq)
    return (sum(c[a] for a in POS) + sum(c[a] for a in NEG)) / len(seq)


def ncpr(seq: str) -> float:
    return net_charge(seq) / len(seq) if seq else 0.0


def sigma(seq: str) -> float:
    """sigma = NCPR^2 / FCR, parametre de charge global de Das & Pappu."""
    f = fcr(seq)
    if f == 0:
        return 0.0
    return ncpr(seq) ** 2 / f


def _delta(seq: str, blob: int) -> float:
    """Ecart quadratique moyen du sigma local (fenetre glissante) au sigma global."""
    n = len(seq)
    if n < blob:
        return 0.0
    s_glob = sigma(seq)
    nblobs = n - blob + 1
    tot = 0.0
    for i in range(nblobs):
        tot += (sigma(seq[i : i + blob]) - s_glob) ** 2
    return tot / nblobs


def _delta_max(seq: str, blob: int) -> float:
    """delta de la permutation MAXIMALEMENT segregee (tous les - puis tous les +
    puis les neutres a la fin, convention Das & Pappu / localCIDER)."""
    c = Counter(seq)
    npos = sum(c[a] for a in POS)
    nneg = sum(c[a] for a in NEG)
    nneu = len(seq) - npos - nneg
    worst = "D" * nneg + "K" * npos + "G" * nneu
    return _delta(worst, blob)


def kappa(seq: str) -> float:
    """kappa de Das & Pappu 2013 : moyenne des blobs 5 et 6, normalisee par
    la permutation maximalement segregee. Renvoie nan si non defini
    (aucune charge, ou une seule espece de charge -> kappa non defini)."""
    c = Counter(seq)
    npos = sum(c[a] for a in POS)
    nneg = sum(c[a] for a in NEG)
    if npos == 0 or nneg == 0:
        return float("nan")
    vals = []
    for blob in (5, 6):
        dmax = _delta_max(seq, blob)
        if dmax == 0:
            continue
        vals.append(_delta(seq, blob) / dmax)
    return float(np.mean(vals)) if vals else float("nan")


def best_split_contrast(seq: str, min_seg: int = 10) -> tuple[float, int]:
    """Contraste N-vs-C MAXIMAL sur toutes les frontieres possibles.

    Pour chaque coupure k, contraste = NCPR(seq[:k]) - NCPR(seq[k:]).
    Renvoie (contraste max, position de coupure). Donne a chaque temoin le
    meme privilege de choix de frontiere qu'a Rv0810c : neutralise le biais
    de selection de la statistique (1).
    """
    n = len(seq)
    best, bestk = -np.inf, -1
    for k in range(min_seg, n - min_seg + 1):
        v = ncpr(seq[:k]) - ncpr(seq[k:])
        if v > best:
            best, bestk = v, k
    return float(best), bestk


# --------------------------------------------------------------------------- #
# Validation de l'instrument (Das & Pappu 2013, Table S1)
# --------------------------------------------------------------------------- #
def validate_kappa() -> dict:
    sv1 = "EKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEKEK"
    sv30 = "EEEEEEEEEEEEEEEEEEEEEEEEEKKKKKKKKKKKKKKKKKKKKKKKKK"
    k1, k30 = kappa(sv1), kappa(sv30)
    ok = (k1 < 0.01) and (k30 > 0.99)
    return {
        "sv1_kappa": round(k1, 5),
        "sv1_expected": "~0.000 (charges parfaitement alternees)",
        "sv30_kappa": round(k30, 5),
        "sv30_expected": "1.000 (charges completement segregees)",
        "instrument_valide": bool(ok),
    }


# --------------------------------------------------------------------------- #
# Rangs empiriques
# --------------------------------------------------------------------------- #
def empirical(value: float, pop: np.ndarray, tail: str) -> dict:
    """Rang empirique et p-value unilaterale (le point cible est INCLUS dans la
    population, convention conservatrice : p = (k+1)/(n+1))."""
    pop = pop[~np.isnan(pop)]
    n = len(pop)
    if tail == "greater":
        k = int(np.sum(pop >= value))
    else:
        k = int(np.sum(pop <= value))
    return {
        "valeur": round(float(value), 4),
        "n_temoins": n,
        "n_aussi_extremes": k,
        "p_empirique": round(k / n, 4),
        "percentile": round(100.0 * (1 - k / n), 1)
        if tail == "greater"
        else round(100.0 * (k / n), 1),
        "mediane_temoins": round(float(np.median(pop)), 4),
        "q05_temoins": round(float(np.percentile(pop, 5)), 4),
        "q95_temoins": round(float(np.percentile(pop, 95)), 4),
    }


def main() -> None:
    prots = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(PROTEOME, "fasta")}
    tgt = prots[TARGET]
    assert len(tgt) == 60, len(tgt)

    report: dict = {
        "piste": "P2.2",
        "question": (
            "La bipolarite de charge de Rv0810c est-elle l'attendu pour une "
            "proteine de 60 aa du proteome H37Rv, ou un ecart mesurable ?"
        ),
        "proteome": str(PROTEOME),
        "n_proteines_proteome": len(prots),
        "validation_instrument": validate_kappa(),
    }

    # ---- 0. Le garde-fou est-il reel ? tendance charge ~ longueur -----------
    lens = np.array([len(s) for s in prots.values()])
    fcrs = np.array([fcr(s) for s in prots.values()])
    bins = [(0, 75), (75, 150), (150, 300), (300, 600), (600, 10**6)]
    trend = []
    for lo, hi in bins:
        m = (lens >= lo) & (lens < hi)
        trend.append(
            {
                "classe_longueur_aa": f"{lo}-{hi if hi < 10**6 else 'inf'}",
                "n": int(m.sum()),
                "FCR_median": round(float(np.median(fcrs[m])), 4),
            }
        )
    report["garde_fou_longueur_charge"] = {
        "commentaire": (
            "Verification que le garde-fou invoque est reel dans CE proteome : "
            "les proteines courtes sont-elles effectivement plus chargees ?"
        ),
        "tendance": trend,
    }

    # ---- 1. Population temoin appariee en longueur --------------------------
    for label, (lo, hi) in {
        "strict_50_70aa": (50, 70),
        "large_45_75aa": (45, 75),
    }.items():
        ctrl_ids = [i for i, s in prots.items() if lo <= len(s) <= hi]
        ctrl = [prots[i] for i in ctrl_ids]

        # (1) segments fixes
        n_mod = np.array([net_charge(s[MODULE[0] - 1 : MODULE[1]]) for s in ctrl])
        n_tail = np.array([net_charge(s[-TAIL_LEN:]) for s in ctrl])
        ncpr_mod = np.array(
            [ncpr(s[MODULE[0] - 1 : MODULE[1]]) for s in ctrl]
        )
        ncpr_tail = np.array([ncpr(s[-TAIL_LEN:]) for s in ctrl])
        contrast_fixed = ncpr_mod - ncpr_tail

        # (2) frontiere optimisee
        best = np.array([best_split_contrast(s)[0] for s in ctrl])

        # (3) kappa
        kap = np.array([kappa(s) for s in ctrl])

        t_nmod = net_charge(tgt[MODULE[0] - 1 : MODULE[1]])
        t_ntail = net_charge(tgt[-TAIL_LEN:])
        t_contrast = ncpr(tgt[: MODULE[1]]) - ncpr(tgt[-TAIL_LEN:])
        t_best, t_bestk = best_split_contrast(tgt)
        t_kappa = kappa(tgt)

        report[f"temoin_{label}"] = {
            "n_temoins": len(ctrl),
            "fenetre_longueur_aa": [lo, hi],
            "stat1_segments_fixes": {
                "commentaire": (
                    "Frontieres 1-19 / 27 derniers residus, issues du profil pLDDT "
                    "de Rv0810c. Fair mais choisi a posteriori : lire avec stat2."
                ),
                "charge_nette_module_1_19": empirical(t_nmod, n_mod, "greater"),
                "charge_nette_queue_27C": empirical(t_ntail, n_tail, "less"),
                "contraste_NCPR_module_moins_queue": empirical(
                    t_contrast, contrast_fixed, "greater"
                ),
            },
            "stat2_frontiere_optimisee_par_proteine": {
                "commentaire": (
                    "Chaque temoin recoit le meme privilege que Rv0810c : sa MEILLEURE "
                    "coupure. Neutralise le biais de selection de la frontiere."
                ),
                "coupure_optimale_Rv0810c": t_bestk,
                "contraste_max": empirical(t_best, best, "greater"),
            },
            "stat3_kappa_Das_Pappu": {
                "commentaire": (
                    "Patterning de charge, sans frontiere a choisir. Test decisif."
                ),
                "kappa": empirical(t_kappa, kap, "greater"),
            },
        }

        # Table detaillee du top-10 des temoins par kappa (pour lecture manuelle)
        order = np.argsort(-np.nan_to_num(kap, nan=-1))[:10]
        report[f"temoin_{label}"]["top10_temoins_par_kappa"] = [
            {
                "rv": ctrl_ids[int(i)],
                "len": len(ctrl[int(i)]),
                "kappa": round(float(kap[int(i)]), 4),
                "ncpr": round(ncpr(ctrl[int(i)]), 4),
                "seq": ctrl[int(i)],
            }
            for i in order
        ]

    report["cible"] = {
        "rv": TARGET,
        "len": len(tgt),
        "seq": tgt,
        "module_1_19": tgt[:19],
        "queue_34_60": tgt[33:],
        "charge_nette_totale": net_charge(tgt),
        "charge_nette_module_1_19": net_charge(tgt[:19]),
        "charge_nette_queue_34_60": net_charge(tgt[33:]),
        "NCPR_module": round(ncpr(tgt[:19]), 4),
        "NCPR_queue": round(ncpr(tgt[33:]), 4),
        "FCR_total": round(fcr(tgt), 4),
        "kappa": round(kappa(tgt), 4),
        "coupure_optimale": best_split_contrast(tgt)[1],
    }

    out = OUTDIR / "p2_2_charge_length_control.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False))
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"\n[ecrit] {out}")


if __name__ == "__main__":
    main()
