#!/usr/bin/env python3
"""
Objet       : P3.5 -- le dernier confondant BIOLOGIQUE du projet. Les souches
              nucS-deletees recensees par ../nucs_deletion_mutators/ (54 souches,
              9 evenements independants) ont un spectre mutationnel deplace vers
              A:T->G:C : leur presence dans un pool tire au sort abaisserait
              mecaniquement le v/u de la lignee hote. Le script fait trois choses
              distinctes : (A) EXPOSITION, combien de ces souches se trouvent dans
              la population de chaque pool et combien ont ete effectivement
              tirees ; (B) SENSIBILITE, de combien le v/u de la lignee hote bouge
              si on les injecte TOUTES dans le pool (borne superieure du
              confondant, volontairement pessimiste) ; (C) CONTROLE POSITIF,
              mesurer le v/u et le GC_eq des souches nucS-deletees elles-memes
              contre des temoins apparies du meme repertoire de clade, ce qui dit
              si l'effet redoute existe et de quelle taille il est.
Entrees     : ../nucs_deletion_mutators/résultats/phase26_reconciliation_54.tsv
              bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt
              data/MTBC0/ancestral_on_H37Rv.bin (via phase2)
Sorties     : résultats/phase4_p35_exposition.tsv
              résultats/phase4_p35_sensibilite.tsv
              résultats/phase4_p35_cas_temoins.tsv
Reutilisable: oui -- le motif « exposition / sensibilite / controle positif »
              vaut pour tout confondant biologique porte par un sous-ensemble
              identifie de souches
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
from phase1_events_vs_strains import BDD, read_subs, strain_dirs, flux  # noqa: E402
from phase2_polarisation_mtbc0 import build_ancestral, load_mask  # noqa: E402
from phase3_counts_par_souche import MASKS  # noqa: E402
from phase3_sueoka_gc_eq import opportunity, vu_ci  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
NUCS_TSV = (ROOT.parent / "nucs_deletion_mutators" / "résultats" /
            "phase26_reconciliation_54.tsv")
CLADES = ["L1", "L2.2.1", "L3", "L4.1.2", "L4.3", "L6.1", "L7", "L9",
          "Orygis_La3", "Caprae_La2", "Microti"]


def load_nucs():
    """{SRA: (clade_declare_par_le_voisin, evenement)} pour les 54 confirmees."""
    out = {}
    for line in NUCS_TSV.read_text().splitlines()[1:]:
        f = line.split("\t")
        if len(f) >= 5 and f[4].startswith("confirme"):
            out[f[0]] = (f[1], f[3])
    return out


def index_bdd():
    """{SRA: clade_du_repertoire} sur toute la base, pour retrouver les souches
    du voisin la ou elles vivent reellement (sa colonne clade vient d'un autre
    systeme de classification)."""
    loc = {}
    for d in BDD.iterdir():
        if not d.is_dir():
            continue
        for s in d.iterdir():
            if s.is_dir() and (s / "NC_000962.3" / "spdi.txt").exists():
                loc[s.name] = d.name
    return loc


def terminal_counts(paths, anc, masked):
    """Pour un pool de souches, compte de substitutions TERMINALES polarisees
    (singletons dans le pool), par souche. Meme geste que phase3_counts."""
    subsets, names = [], []
    for p in paths:
        v = read_subs(p / "NC_000962.3" / "spdi.txt")
        if v:
            subsets.append(v)
            names.append(p.name)
    support = defaultdict(int)
    for i, subs in enumerate(subsets):
        for v in subs:
            support[v] |= 1 << i
    per = defaultdict(Counter)
    for (pos, ref, alt), mask in support.items():
        if mask & (mask - 1):          # k >= 2 : pas un singleton
            continue
        if pos in masked:
            continue
        a = chr(anc[pos]) if pos < len(anc) else "N"
        if a == "N" or a == alt:       # non liftee, ou H37Rv derive
            continue
        if a != ref:
            ref = a                    # allele ancestral tierce
        per[names[mask.bit_length() - 1]][flux(ref, alt)] += 1
    return {n: per[n] for n in names}


def sample_seed0(clade, n):
    st = strain_dirs(clade)
    rng = random.Random(0)
    rng.shuffle(st)
    return st[:n], st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--n-per-clade", type=int, default=40)
    ap.add_argument("--n-temoins", type=int, default=150,
                    help="temoins tires par repertoire de clade en section C")
    ap.add_argument("--out-prefix", default=str(ROOT / "résultats" / "phase4_p35_"))
    args = ap.parse_args()

    anc = build_ancestral()
    masked = load_mask(MASKS, in_mtbc0_coords=True)
    (n_gc, n_at), _ = opportunity()
    print(f"# masque {len(masked)} positions ; opportunite {n_gc} G:C / {n_at} A:T",
          file=sys.stderr)

    nucs = load_nucs()
    loc = index_bdd()
    present = {k: loc[k] for k in nucs if k in loc}
    print(f"# {len(nucs)} souches nucS-deletees confirmees, {len(present)} "
          f"presentes dans bdd/actuelle avec spdi.txt", file=sys.stderr)

    # ---------- A. EXPOSITION -------------------------------------------------
    rows, exposees = [], {}
    for clade in CLADES:
        sample, pop = sample_seed0(clade, args.n_per_clade)
        popn = {p.name for p in pop}
        dans_pop = sorted(popn & set(present))
        dans_ech = sorted({p.name for p in sample} & set(present))
        exposees[clade] = dans_pop
        rows.append(dict(clade=clade, n_population=len(popn),
                         n_nucs_population=len(dans_pop),
                         prevalence_ppm=1e6 * len(dans_pop) / max(1, len(popn)),
                         n_echantillon=len(sample), n_nucs_echantillon=len(dans_ech),
                         souches_nucs=",".join(dans_pop)))
    ta = pd.DataFrame(rows)
    ta.to_csv(args.out_prefix + "exposition.tsv", sep="\t", index=False)
    print("\n=== A. EXPOSITION : ou sont les souches nucS-deletees ===")
    print(ta.drop(columns=["souches_nucs"]).to_string(index=False))

    # ---------- B. SENSIBILITE (borne superieure pessimiste) ------------------
    rows = []
    for clade in CLADES:
        if not exposees[clade]:
            continue
        sample, pop = sample_seed0(clade, args.n_per_clade)
        byname = {p.name: p for p in pop}
        base = terminal_counts(sample, anc, masked)
        lb = sum(c["loss"] for c in base.values())
        gb = sum(c["gain"] for c in base.values())
        # injection : les souches nucS remplacent les dernieres du tirage,
        # l'effectif du pool reste n (le v/u depend du pool par la definition
        # meme de « terminal »)
        add = [byname[s] for s in exposees[clade] if s not in {p.name for p in sample}]
        spiked = sample[:args.n_per_clade - len(add)] + add
        sp = terminal_counts(spiked, anc, masked)
        ls = sum(c["loss"] for c in sp.values())
        gs = sum(c["gain"] for c in sp.values())
        vb, _, _ = vu_ci(lb, gb, n_gc, n_at)
        vs, _, _ = vu_ci(ls, gs, n_gc, n_at)
        rows.append(dict(clade=clade, n_injectees=len(add),
                         frac_pool=len(add) / args.n_per_clade,
                         loss_base=lb, gain_base=gb, vu_base=vb,
                         gc_eq_base=1 / (1 + vb),
                         loss_spike=ls, gain_spike=gs, vu_spike=vs,
                         gc_eq_spike=1 / (1 + vs),
                         delta_vu_pct=100 * (vs - vb) / vb,
                         delta_gc_eq_pts=100 * (1 / (1 + vs) - 1 / (1 + vb))))
    tb = pd.DataFrame(rows)
    tb.to_csv(args.out_prefix + "sensibilite.tsv", sep="\t", index=False)
    print("\n=== B. SENSIBILITE : injection FORCEE de toutes les nucS du pool ===")
    print(tb.round(4).to_string(index=False))

    # ---------- C. CONTROLE POSITIF : cas contre temoins apparies -------------
    par_dir = defaultdict(list)
    for s, d in present.items():
        par_dir[d].append(s)
    rows = []
    for d, cas in sorted(par_dir.items()):
        pool_dir = sorted(p for p in (BDD / d).iterdir()
                          if p.is_dir() and (p / "NC_000962.3" / "spdi.txt").exists())
        byname = {p.name: p for p in pool_dir}
        cas = [c for c in cas if c in byname]
        if not cas:
            continue
        temoins = [p for p in pool_dir if p.name not in set(cas)]
        rng = random.Random(0)
        rng.shuffle(temoins)
        temoins = temoins[:args.n_temoins]
        pool = [byname[c] for c in cas] + temoins
        cnt = terminal_counts(pool, anc, masked)
        for name, c in cnt.items():
            rows.append(dict(clade_dir=d, sra=name,
                             groupe="cas" if name in set(cas) else "temoin",
                             n_pool=len(pool), loss=c["loss"], gain=c["gain"]))
    tc = pd.DataFrame(rows)
    tc.to_csv(args.out_prefix + "cas_temoins.tsv", sep="\t", index=False)

    print("\n=== C. CONTROLE POSITIF : nucS-deletees contre temoins du meme "
          "repertoire de clade ===")
    print(f"{'repertoire':28s} {'groupe':7s} {'n':>4s} {'loss':>7s} {'gain':>7s} "
          f"{'v/u':>7s} {'IC95':>15s} {'GC_eq':>7s}")
    for d, g in tc.groupby("clade_dir"):
        for grp in ("cas", "temoin"):
            s = g[g.groupe == grp]
            L, G = int(s["loss"].sum()), int(s["gain"].sum())
            if G == 0 or L == 0:
                continue
            vu, lo, hi = vu_ci(L, G, n_gc, n_at)
            print(f"{d:28s} {grp:7s} {len(s):4d} {L:7d} {G:7d} {vu:7.3f} "
                  f"{lo:6.3f}-{hi:6.3f} {100/(1+vu):6.1f}%")
    print("-" * 90)
    for grp in ("cas", "temoin"):
        s = tc[tc.groupe == grp]
        L, G = int(s["loss"].sum()), int(s["gain"].sum())
        vu, lo, hi = vu_ci(L, G, n_gc, n_at)
        print(f"{'TOUS REPERTOIRES':28s} {grp:7s} {len(s):4d} {L:7d} {G:7d} "
              f"{vu:7.3f} {lo:6.3f}-{hi:6.3f} {100/(1+vu):6.1f}%")

    # Test STRATIFIE par repertoire de clade (Mantel-Haenszel) : le test agrege
    # melangerait des repertoires dont le v/u temoin va de 0,98 a 1,55, ce qui
    # est exactement la situation ou un 2x2 poole se laisse retourner (Simpson).
    from scipy import stats
    strates, num, den, mh_num, mh_den = [], 0.0, 0.0, 0.0, 0.0
    for d, g in tc.groupby("clade_dir"):
        a = int(g[g.groupe == "cas"]["loss"].sum())
        b = int(g[g.groupe == "cas"]["gain"].sum())
        c = int(g[g.groupe == "temoin"]["loss"].sum())
        e = int(g[g.groupe == "temoin"]["gain"].sum())
        n = a + b + c + e
        if min(a + b, c + e) == 0 or n == 0:
            continue
        strates.append((d, a, b, c, e))
        num += a - (a + b) * (a + c) / n
        den += (a + b) * (c + e) * (a + c) * (b + e) / (n ** 2 * (n - 1))
        mh_num += a * e / n
        mh_den += b * c / n
    cmh = (abs(num) - 0.5) ** 2 / den
    p_cmh = 1 - stats.chi2.cdf(cmh, 1)
    or_mh = mh_num / mh_den
    print(f"\nMantel-Haenszel stratifie par repertoire ({len(strates)} strates "
          f"informatives) : chi2 = {cmh:.2f}, p = {p_cmh:.3g}, "
          f"OR commun cas/temoin = {or_mh:.3f}")
    print("  (OR < 1 = les cas perdent MOINS de paires G:C par gain que leurs "
          "temoins, sens attendu d'une perte de NucS)")
    for d, a, b, c, e in strates:
        orr = (a * e) / (b * c) if b * c else float("nan")
        print(f"    {d:26s} cas {a:4d}/{b:4d}  temoins {c:6d}/{e:6d}  OR = {orr:.3f}")
    a = tc[tc.groupe == "cas"][["loss", "gain"]].sum()
    b = tc[tc.groupe == "temoin"][["loss", "gain"]].sum()
    tab = np.array([[a["loss"], a["gain"]], [b["loss"], b["gain"]]])
    chi2, p, _, _ = stats.chi2_contingency(tab)
    print(f"  pour memoire, agrege non stratifie : chi2 = {chi2:.2f}, p = {p:.3g}, "
          f"OR = {(tab[0,0]*tab[1,1])/(tab[0,1]*tab[1,0]):.3f}")

    print(f"\n# ecrit dans {args.out_prefix}[exposition|sensibilite|cas_temoins].tsv",
          file=sys.stderr)


if __name__ == "__main__":
    main()
