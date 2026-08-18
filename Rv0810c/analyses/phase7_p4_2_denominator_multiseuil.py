#!/usr/bin/env python3
"""
P4.2 (suite) — Le denominateur, lu a plusieurs seuils et STRATIFIE EN LONGUEUR.

Ce que la phase 3 a revele en cours de route, et qui commande tout le reste :
un seuil unique de E <= 1e-5 sur une recherche PROTEOME -> GENOME ENTIER est
INATTEIGNABLE pour une proteine de 60 aa a distance de *Tropheryma*. Mesure sur
la cible elle-meme :

    Rv0810c vs genome de T. whipplei (tblastn)      E = 8,8e-04   -> "absent"
    Rv0810c vs proteine TWT_722 (blastp)            E = 1,1e-12   -> orthologue
                                                    57,6 % id sur les residus 1-33

L'orthologie est donc certaine, et c'est le SEUIL qui est sous-puissant, pas le
gene qui est perdu. En tirer un taux de retention a seuil unique reviendrait a
mesurer la sensibilite de l'instrument et a l'appeler evolution reductive.

D'ou ce script : le meme calcul a trois seuils, stratifie en longueur, avec le
taux CONDITIONNEL a la detectabilite dans un temoin non reduit du meme clade.

Entree : les sorties tblastn brutes conservees par phase3 (mappees par mtime).
Sortie : résultats/p4_2_denominateur_multiseuil.json
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from Bio import SeqIO
from scipy.stats import beta as beta_dist

PROJ = Path(__file__).resolve().parent.parent
MTBC = PROJ.parent
PROTEOME = MTBC / "annotation_mtbc" / "résultats" / "phase2d_eggnog" / "proteome.faa"
ATLAS = MTBC / "annotation_mtbc" / "site" / "content" / "genes"
WORK = PROJ / "experiments" / "2026-08-10_P4.2_work"
OUTDIR = PROJ / "résultats"

TARGET = "Rv0810c"
# Les quatre premiers genomes ont ete cribles sur le proteome ENTIER (fichiers
# temporaires de phase3, ordonnes par mtime). Le cinquieme, M. luteus, a du etre
# restreint : a 73 % de GC et avec -seg no, le tblastn proteome-entier produisait
# 5 ko/min (>2 h de plus). Il a ete relance sur le sous-ensemble des proteines
# <= 200 aa, qui est exactement la strate qui porte l'argument (Rv0810c fait
# 60 aa). RESTRICTION EXPLICITEMENT DOCUMENTEE : aucun chiffre M. luteus n'est
# rapporte au-dela de 200 aa, plutot qu'un plafond silencieux.
ORDER = ["M_leprae", "M_lepromatosis", "M_abscessus", "T_whipplei"]
LUTEUS_TSV = WORK / "short200_vs_Mluteus.tsv"
LUTEUS_DONE = WORK / "short200_vs_Mluteus.done"
LUTEUS_MAXLEN = 200
REDUIT = {"M_leprae": "M_abscessus", "M_lepromatosis": "M_abscessus",
          "T_whipplei": "M_luteus"}
SEUILS = [("strict", 1e-5), ("intermediaire", 1e-3), ("permissif", 1e-2)]
QCOV_MIN = 40.0


def _parse(f: Path) -> dict[str, dict]:
    best: dict[str, dict] = {}
    for line in f.read_text().splitlines():
        c = line.split("\t")
        if len(c) < 7:
            continue
        q, ev, bs, qc = c[0], float(c[4]), float(c[5]), float(c[6])
        cur = best.get(q)
        if cur is None or bs > cur["bitscore"]:
            best[q] = {"evalue": ev, "bitscore": bs, "qcovs": qc}
        else:
            cur["qcovs"] = max(cur["qcovs"], qc)
            cur["evalue"] = min(cur["evalue"], ev)
    return best


def load_hits() -> dict[str, dict[str, dict]]:
    tsvs = sorted(
        (p for p in WORK.glob("tmp*.tsv")), key=lambda p: p.stat().st_mtime
    )[: len(ORDER)]
    if len(tsvs) != len(ORDER):
        raise SystemExit(f"{len(tsvs)} fichiers tblastn pour {len(ORDER)} genomes")
    out: dict[str, dict[str, dict]] = {}
    for name, f in zip(ORDER, tsvs):
        out[name] = _parse(f)
        print(f"  {name:16} <- {f.name}  ({len(out[name])} requetes avec hit)")
    # M_luteus n'est integre que si son criblage est ACHEVE (fichier sentinelle).
    # Un tblastn encore en cours donnerait un taux de retention artificiellement
    # bas, indiscernable d'une perte reelle : mieux vaut l'absence que le faux.
    if LUTEUS_DONE.exists():
        out["M_luteus"] = _parse(LUTEUS_TSV)
        print(
            f"  {'M_luteus':16} <- {LUTEUS_TSV.name}  "
            f"({len(out['M_luteus'])} requetes ; RESTREINT a <= {LUTEUS_MAXLEN} aa)"
        )
    else:
        print("  M_luteus         ABSENT : criblage non acheve, temoin non rapporte")
    return out


def ci(k: int, n: int) -> list[float]:
    lo = beta_dist.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = beta_dist.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return [round(100 * float(lo), 1), round(100 * float(hi), 1)]


def rate(mask: np.ndarray, ret: np.ndarray) -> dict:
    n = int(mask.sum())
    if n == 0:
        return {"n": 0, "taux_pct": None}
    k = int(ret[mask].sum())
    return {"n": n, "retenus": k, "taux_pct": round(100 * k / n, 1), "IC95": ci(k, n)}


def main() -> None:
    prots = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(PROTEOME, "fasta")}
    ids = sorted(prots)
    lens = np.array([len(prots[i]) for i in ids])
    ess, hyp = [], []
    for i in ids:
        f = ATLAS / f"{i}.json"
        d = json.loads(f.read_text()) if f.exists() else {}
        ess.append((d.get("essentiality") or {}).get("dejesus2017") == "ES")
        hyp.append((d.get("funccat") or {}).get("category") == "conserved hypotheticals")
    ess, hyp = np.array(ess), np.array(hyp)
    ti = ids.index(TARGET)
    hits = load_hits()

    rep: dict = {
        "piste": "P4.2 (multi-seuils)",
        "probleme_resolu": {
            "constat": (
                "Un seuil tblastn unique E<=1e-5 proteome->genome est inatteignable "
                "pour une proteine de 60 aa a distance de Tropheryma."
            ),
            "mesure_sur_la_cible": {
                "Rv0810c_vs_genome_T_whipplei_tblastn": {"evalue": 8.8e-4, "qcov_pct": 75},
                "Rv0810c_vs_proteine_TWT_722_blastp": {
                    "evalue": 1.1e-12, "identite_pct": 57.6,
                    "residus_alignes": "1-33 (le module ordonne)",
                    "longueur_TWT_722_aa": 41,
                },
            },
            "consequence": (
                "L'orthologie est certaine ; c'est le seuil qui est sous-puissant. "
                "Tout taux de retention lu a seuil unique melange perte reelle et "
                "limite de detection, et cette limite depend de la LONGUEUR."
            ),
        },
        "seuils": {},
    }

    bandes = [(1, 100), (101, 200), (201, 400), (401, 10_000)]
    for lbl, e_max in SEUILS:
        pres = {
            g: np.array(
                [
                    (i in h) and h[i]["evalue"] <= e_max and h[i]["qcovs"] >= QCOV_MIN
                    for i in ids
                ]
            )
            for g, h in hits.items()
        }
        # M. luteus n'a ete crible que sur les proteines <= 200 aa : hors de cette
        # strate, "non detecte" signifierait "non teste". On masque explicitement.
        testable_luteus = lens <= LUTEUS_MAXLEN
        bloc: dict = {"E_max": e_max, "qcov_min_pct": QCOV_MIN, "genomes": {}}
        for g in [x for x in ORDER + ["M_luteus"] if x in hits]:
            restreint = g == "M_luteus"
            dom = testable_luteus if restreint else np.ones(len(ids), bool)
            gb: dict = {
                "criblage": (
                    f"proteines <= {LUTEUS_MAXLEN} aa UNIQUEMENT"
                    if restreint
                    else "proteome entier"
                ),
                "global": rate(dom, pres[g]),
                "par_longueur": {
                    f"{lo}-{hi if hi < 10_000 else 'inf'} aa": rate(
                        dom & (lens >= lo) & (lens <= hi), pres[g]
                    )
                    for lo, hi in bandes
                    if not (restreint and lo > LUTEUS_MAXLEN)
                },
                "essentiels_ES": rate(dom & ess, pres[g]),
                "conserved_hypotheticals": rate(dom & hyp, pres[g]),
                "Rv0810c_detecte": bool(pres[g][ti]),
                "Rv0810c_evalue": hits[g].get(TARGET, {}).get("evalue"),
                "Rv0810c_qcov": hits[g].get(TARGET, {}).get("qcovs"),
            }
            ctrl = REDUIT.get(g)
            if ctrl and ctrl in hits:
                # Le temoin M. luteus n'existe que sous 200 aa : le conditionnement
                # est donc lui aussi restreint a cette strate, et c'est dit.
                base = pres[ctrl] & (testable_luteus if ctrl == "M_luteus" else True)
                gb["retention_conditionnelle_au_temoin"] = {
                    "temoin_non_reduit": ctrl,
                    "domaine": (
                        f"proteines <= {LUTEUS_MAXLEN} aa"
                        if ctrl == "M_luteus"
                        else "proteome entier"
                    ),
                    "tous": rate(base, pres[g]),
                    "longueur_1_100aa": rate(base & (lens <= 100), pres[g]),
                    "longueur_50_70aa": rate(base & (lens >= 50) & (lens <= 70), pres[g]),
                    "essentiels_ES": rate(base & ess, pres[g]),
                    "essentiels_ES_et_<=100aa": rate(
                        base & ess & (lens <= 100), pres[g]
                    ),
                    "conserved_hypotheticals": rate(base & hyp, pres[g]),
                }
            bloc["genomes"][g] = gb
        rep["seuils"][lbl] = bloc

    out = OUTDIR / "p4_2_denominateur_multiseuil.json"
    out.write_text(json.dumps(rep, indent=2, ensure_ascii=False))
    print(json.dumps(rep["seuils"], indent=2, ensure_ascii=False))
    print(f"\n[ecrit] {out}")


if __name__ == "__main__":
    main()
