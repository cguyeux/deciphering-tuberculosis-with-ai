#!/usr/bin/env python3
"""
Objet       : P5.1 -- refaire la comparaison inter-lignees en DESCENDANT
              recursivement dans les sous-repertoires de clade. Le sondage
              initial s'etait arrete au conteneur nu <L>/ et lisait donc n = 1
              souche pour L6.1, L9, Bovis et Orygis ; le panel de onze lignees
              de A15 descend deja par prefixe, mais il a ete choisi sur les
              repertoires les plus peuples et non a la maille lignee, et il
              manque un alias (la lignee 1 vit dans DEUX conteneurs disjoints,
              `L1/` nu et `L_1.*`, que le prefixe `L1` ne reunit pas). Ce script
              (A) audite la couverture de chaque lignee, premier niveau contre
              descente recursive ; (B) mesure v/u et GC_eq sur un panel elargi a
              la maille lignee, y compris L5, L8, L10 et Bovis absents de A15 ;
              (C) teste l'EMBOITEMENT mere/fille, qui dit si la valeur d'une
              lignee depend du niveau taxonomique auquel on la lit ; (D) teste
              la sensibilite a l'effectif n du pool, qui dit si comparer des
              lignees a effectifs inegaux est legitime.
Entrees     : bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt
              data/MTBC0/ancestral_on_H37Rv.bin, data/MTBC0/Mask.files/*.bed
Sorties     : résultats/phase5_p51_couverture.tsv       (audit A)
              résultats/phase5_p51_counts_par_souche.tsv (B, une ligne/souche)
              résultats/phase5_p51_vu_panel.tsv          (B, une ligne/lignee)
              résultats/phase5_p51_emboitement.tsv       (C)
              résultats/phase5_p51_sensibilite_n.tsv     (D)
Reutilisable: oui -- l'audit de couverture et la resolution d'alias valent pour
              toute analyse par lignee lisant bdd/actuelle
Projet      : GC_par_lignee
Date        : 2026-08-29
"""
import argparse
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import BDD, read_subs, flux, canonical  # noqa: E402
from phase2_polarisation_mtbc0 import build_ancestral, load_mask  # noqa: E402
from phase3_sueoka_gc_eq import opportunity, vu_ci  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
MASKS = sorted((ROOT / "data" / "MTBC0" / "Mask.files").glob("*.bed"))
MASK_CACHE = ROOT / "data" / "mask_h37rv_positions.npy"
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

# Panel a la maille LIGNEE. La valeur est la liste des PREFIXES de repertoire a
# reunir : la lignee 1 en demande deux, verifies porteurs des memes marqueurs
# positifs de L1 du barcode central (5/5 sur L1 nu comme sur L_1.A/B/C).
PANEL = {
    "L1": ["L1", "L_1"], "L2": ["L2"], "L3": ["L3"], "L4": ["L4"],
    "L5": ["L5"], "L6": ["L6"], "L7": ["L7"], "L8": ["L8"], "L9": ["L9"],
    "L10": ["L10"], "Bovis": ["Bovis"], "Caprae_La2": ["Caprae_La2"],
    "Orygis_La3": ["Orygis_La3"], "Microti": ["Microti"],
    "Pinnipedii": ["Pinnipedii"], "La4": ["La4"], "Dassie": ["Dassie"],
    "Suricattae": ["Suricattae"], "Mungi": ["Mungi"], "Chimpanze": ["Chimpanze"],
    "Borstel": ["Borstel"],
}
# Pools du panel A15, gardes pour l'emboitement mere / fille.
FILLES = {"L2.2.1": ["L2.2.1"], "L4.1.2": ["L4.1.2"], "L4.3": ["L4.3"],
          "L6.1": ["L6.1"], "L1_nu_seul": ["L1"], "L1_alias_seul": ["L_1"]}


def pool_dirs(prefixes):
    """Repertoires de clade reunis par un ou plusieurs prefixes (convention
    dir-mixte : conteneur nu <p>/ + sous-clades <p>.*)."""
    out = []
    for d in sorted(BDD.iterdir()):
        if not d.is_dir():
            continue
        if any(d.name == p or d.name.startswith(p + ".") for p in prefixes):
            out.append(d)
    return out


def strains_of(prefixes, first_level_only=False):
    out = []
    for d in pool_dirs(prefixes):
        if first_level_only and d.name not in prefixes:
            continue
        for s in sorted(d.iterdir()):
            if s.is_dir() and (s / "NC_000962.3" / "spdi.txt").exists():
                out.append(s)
    return out


def masked_positions():
    if MASK_CACHE.exists():
        return set(np.load(MASK_CACHE).tolist())
    m = load_mask(MASKS, in_mtbc0_coords=True)
    np.save(MASK_CACHE, np.fromiter(sorted(m), dtype=np.int64))
    return m


def counts_par_souche(strains, anc, masked, n, seed=0):
    """Comptes de branche terminale (variants portes par UNE souche du pool)
    pour un pool tire au sort. Meme mecanique que phase3_counts_par_souche."""
    rng = random.Random(seed)
    pool = list(strains)
    rng.shuffle(pool)
    subsets, names = [], []
    for s in pool[:n]:
        v = read_subs(s / "NC_000962.3" / "spdi.txt")
        if v:
            subsets.append(v)
            names.append(s.name)
    k = len(subsets)
    if k < 4:
        return None, k
    support = defaultdict(int)
    for i, subs in enumerate(subsets):
        for v in subs:
            support[v] |= 1 << i
    per = defaultdict(Counter)
    for (pos, ref, alt), mask in support.items():
        if mask & (mask - 1) or pos in masked:
            continue
        i = mask.bit_length() - 1
        c = per[i]
        c["n_term_raw"] += 1
        a = chr(anc[pos]) if pos < len(anc) else "N"
        if a == "N":
            continue
        if a == alt:
            c["inverse"] += 1
            continue
        if a != ref:
            c["tierce"] += 1
            ref = a
        c[flux(ref, alt)] += 1
        c[canonical(ref, alt)] += 1
    rows = []
    for i, name in enumerate(names):
        c = per[i]
        rows.append(dict(sra=name, n_pool=k, n_term_raw=c["n_term_raw"],
                         loss=c["loss"], gain=c["gain"], neutral=c["neutral"],
                         inverse=c["inverse"], tierce=c["tierce"],
                         **{cl: c[cl] for cl in CLASSES}))
    return pd.DataFrame(rows), k


def vu_row(name, df, n_gc, n_at, extra=None):
    loss, gain = int(df["loss"].sum()), int(df["gain"].sum())
    if gain == 0 or loss == 0:
        return dict(pool=name, n_souches=len(df), loss=loss, gain=gain,
                    ratio_brut=float("nan"), vu=float("nan"),
                    vu_lo=float("nan"), vu_hi=float("nan"),
                    gc_eq=float("nan"), **(extra or {}))
    vu, lo, hi = vu_ci(loss, gain, n_gc, n_at)
    return dict(pool=name, n_souches=len(df), loss=loss, gain=gain,
                ratio_brut=loss / gain, vu=vu, vu_lo=lo, vu_hi=hi,
                gc_eq=1 / (1 + vu), gc_eq_lo=1 / (1 + hi), gc_eq_hi=1 / (1 + lo),
                **(extra or {}))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--n-per-clade", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-strains", type=int, default=4)
    ap.add_argument("--out-prefix",
                    default=str(ROOT / "résultats" / "phase5_p51_"))
    args = ap.parse_args()

    anc = build_ancestral()
    masked = masked_positions()
    (n_gc, n_at), _ = opportunity()
    print(f"# masque : {len(masked)} positions ; opportunite : {n_gc} G:C, "
          f"{n_at} A:T", file=sys.stderr)

    # ---- A. audit de couverture : premier niveau contre descente recursive
    audit = []
    for lin, prefixes in PANEL.items():
        rec = strains_of(prefixes)
        nu = strains_of(prefixes, first_level_only=True)
        audit.append(dict(lignee=lin, prefixes="+".join(prefixes),
                          n_dirs=len(pool_dirs(prefixes)),
                          n_souches_premier_niveau=len(nu),
                          n_souches_recursif=len(rec),
                          gain_absolu=len(rec) - len(nu),
                          facteur=len(rec) / max(len(nu), 1)))
    a = pd.DataFrame(audit).sort_values("n_souches_recursif", ascending=False)
    a.to_csv(args.out_prefix + "couverture.tsv", sep="\t", index=False)
    print("\n=== A. couverture : ce que la descente recursive ajoute ===")
    print(a.to_string(index=False))

    # ---- B. panel elargi a la maille lignee
    rows_strain, rows_vu = [], []
    for lin, prefixes in PANEL.items():
        st = strains_of(prefixes)
        if len(st) < args.min_strains:
            print(f"# {lin} : {len(st)} souche(s), sous le seuil", file=sys.stderr)
            continue
        df, k = counts_par_souche(st, anc, masked, args.n_per_clade, args.seed)
        if df is None:
            print(f"# {lin} : pool trop petit ({k})", file=sys.stderr)
            continue
        df.insert(0, "lignee", lin)
        rows_strain.append(df)
        rows_vu.append(vu_row(lin, df, n_gc, n_at,
                              dict(n_dispo=len(st), n_dirs=len(pool_dirs(prefixes)))))
    cs = pd.concat(rows_strain, ignore_index=True)
    cs.to_csv(args.out_prefix + "counts_par_souche.tsv", sep="\t", index=False)
    v = pd.DataFrame(rows_vu).sort_values("vu", ascending=False)
    v.to_csv(args.out_prefix + "vu_panel.tsv", sep="\t", index=False)
    print("\n=== B. panel elargi a la maille lignee (n = "
          f"{args.n_per_clade}) ===")
    print(v[["pool", "n_dispo", "n_souches", "loss", "gain", "vu", "vu_lo",
             "vu_hi", "gc_eq"]].round(4).to_string(index=False))

    # ---- C. emboitement mere / fille
    emb = []
    for name, prefixes in FILLES.items():
        st = strains_of(prefixes)
        if len(st) < args.min_strains:
            continue
        df, k = counts_par_souche(st, anc, masked, args.n_per_clade, args.seed)
        if df is None:
            continue
        emb.append(vu_row(name, df, n_gc, n_at, dict(n_dispo=len(st))))
    e = pd.DataFrame(emb)
    e.to_csv(args.out_prefix + "emboitement.tsv", sep="\t", index=False)
    print("\n=== C. emboitement : la valeur depend-elle du niveau lu ? ===")
    print(e[["pool", "n_dispo", "n_souches", "loss", "gain", "vu", "vu_lo",
             "vu_hi", "gc_eq"]].round(4).to_string(index=False))
    for mere, filles in [("L1", ["L1_nu_seul", "L1_alias_seul"]),
                         ("L2", ["L2.2.1"]), ("L4", ["L4.1.2", "L4.3"]),
                         ("L6", ["L6.1"])]:
        m = v[v["pool"] == mere]
        if m.empty:
            continue
        mv = m["vu"].iloc[0]
        for f in filles:
            row = e[e["pool"] == f]
            if row.empty:
                continue
            fv = row["vu"].iloc[0]
            chevauche = (row["vu_lo"].iloc[0] <= m["vu_hi"].iloc[0] and
                         m["vu_lo"].iloc[0] <= row["vu_hi"].iloc[0])
            print(f"  {mere} (v/u {mv:.3f}) vs {f} (v/u {fv:.3f}) : "
                  f"ecart {100*(fv-mv)/mv:+.1f} %, IC95 "
                  f"{'chevauchants' if chevauche else 'DISJOINTS'}")

    # ---- D. sensibilite a l'effectif du pool
    sens = []
    for lin in ["L1", "L2", "L4", "L6", "Bovis", "L9"]:
        if lin not in PANEL:
            continue
        st = strains_of(PANEL[lin])
        for n in [4, 8, 16, 40, 80]:
            if len(st) < n:
                continue
            df, k = counts_par_souche(st, anc, masked, n, args.seed)
            if df is None:
                continue
            sens.append(vu_row(lin, df, n_gc, n_at, dict(n_cible=n)))
    s = pd.DataFrame(sens)
    s.to_csv(args.out_prefix + "sensibilite_n.tsv", sep="\t", index=False)
    print("\n=== D. sensibilite a l'effectif du pool ===")
    print(s.pivot_table(index="pool", columns="n_cible",
                        values="vu").round(3).to_string())
    print("\n# ecrits : " + args.out_prefix + "{couverture,counts_par_souche,"
          "vu_panel,emboitement,sensibilite_n}.tsv")


if __name__ == "__main__":
    main()
