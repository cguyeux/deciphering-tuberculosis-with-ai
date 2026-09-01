#!/usr/bin/env python3
"""
Objet       : P5.2 -- donner a la comparaison inter-lignees l'incertitude qui
              lui manque et annoncer sa puissance. Les IC de A15 sont des IC
              BINOMIAUX : ils traitent chaque substitution comme une observation
              independante, alors que les substitutions sont groupees par souche
              et que les souches d'une lignee partagent une histoire. Ce script
              (A) refait les IC par BOOTSTRAP SUR SOUCHES, seule unite
              reellement echangeable, et chiffre la sur-dispersion que l'IC
              binomial ignorait ; (B) teste l'homogeneite du spectre a 6 classes
              entre lignees par permutation des etiquettes AU NIVEAU DES SOUCHES,
              le seul nul qui respecte le groupement ; (C) fait les comparaisons
              par paires avec correction FDR ; (D) annonce la puissance, c'est-a-
              dire l'ecart minimal detectable, lignee par lignee -- ce qui dit
              d'avance ce que L8, L9, L10 et les ecotypes rares ne pourront
              jamais trancher, quel que soit le soin de la mesure.
Entrees     : résultats/phase5_p51_counts_par_souche.tsv (P5.1)
Sorties     : résultats/phase5_p52_ic_bootstrap.tsv
              résultats/phase5_p52_homogeneite_spectre.tsv
              résultats/phase5_p52_paires_fdr.tsv
              résultats/phase5_p52_puissance.tsv
Reutilisable: oui -- le bootstrap sur unite groupee et le nul par permutation
              d'etiquettes valent pour toute comparaison de spectres par clade
Projet      : GC_par_lignee
Date        : 2026-08-29
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase3_sueoka_gc_eq import opportunity, vu_ci  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


def boot_vu(loss, gain, n_gc, n_at, B, rng, hierarchique=True):
    """Bootstrap de v/u. Deux niveaux quand `hierarchique` : on retire d'abord
    les SOUCHES avec remise (variance entre souches, la seule unite echangeable),
    puis on retire les COMPTES de chaque souche tiree sous une loi de Poisson
    (variance de comptage a l'interieur d'une souche). Le second niveau n'est pas
    un raffinement : sans lui, un pool de quatre souches ne peut produire qu'une
    poignee de sommes distinctes et le bootstrap rend un IC plus etroit que l'IC
    binomial -- c'est-a-dire faux dans le sens le plus dangereux."""
    k = len(loss)
    idx = rng.integers(0, k, size=(B, k))
    L, G = loss[idx], gain[idx]
    if hierarchique:
        L, G = rng.poisson(L), rng.poisson(G)
    L = L.sum(axis=1).astype(float)
    G = G.sum(axis=1).astype(float)
    ok = (L > 0) & (G > 0)
    return (L[ok] / n_gc) / (G[ok] / n_at)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default=str(ROOT / "résultats" /
                                            "phase5_p51_counts_par_souche.tsv"))
    ap.add_argument("-B", "--boot", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-prefix",
                    default=str(ROOT / "résultats" / "phase5_p52_"))
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    (n_gc, n_at), _ = opportunity()
    df = pd.read_csv(args.counts, sep="\t")
    lignees = list(df.groupby("lignee")["loss"].sum().index)

    # ---- A. IC bootstrap sur souches contre IC binomial
    rows = []
    dist = {}
    for lin in lignees:
        d = df[df["lignee"] == lin]
        loss, gain = d["loss"].to_numpy(), d["gain"].to_numpy()
        L, G = int(loss.sum()), int(gain.sum())
        if L == 0 or G == 0:
            continue
        vu, blo, bhi = vu_ci(L, G, n_gc, n_at)
        bs = boot_vu(loss, gain, n_gc, n_at, args.boot, rng)
        bs1 = boot_vu(loss, gain, n_gc, n_at, args.boot, rng, hierarchique=False)
        dist[lin] = bs
        lo, hi = np.percentile(bs, [2.5, 97.5])
        se_log = np.std(np.log(bs), ddof=1)
        se_log1 = np.std(np.log(bs1), ddof=1)
        rows.append(dict(
            lignee=lin, n_souches=len(d), loss=L, gain=G, vu=vu,
            ic_binom_lo=blo, ic_binom_hi=bhi,
            ic_boot_lo=lo, ic_boot_hi=hi,
            gc_eq=1 / (1 + vu), gc_eq_boot_lo=1 / (1 + hi),
            gc_eq_boot_hi=1 / (1 + lo),
            se_log_boot=se_log,
            se_log_boot_souches_seules=se_log1,
            se_log_binom=np.sqrt(1 / L + 1 / G),
            surdispersion=se_log / np.sqrt(1 / L + 1 / G),
            largeur_binom=np.log(bhi) - np.log(blo),
            largeur_boot=np.log(hi) - np.log(lo)))
    ic = pd.DataFrame(rows).sort_values("vu", ascending=False)
    ic.to_csv(args.out_prefix + "ic_bootstrap.tsv", sep="\t", index=False)
    print("=== A. IC bootstrap sur souches contre IC binomial ===")
    print(ic[["lignee", "n_souches", "vu", "ic_binom_lo", "ic_binom_hi",
              "ic_boot_lo", "ic_boot_hi", "gc_eq", "gc_eq_boot_lo",
              "gc_eq_boot_hi", "surdispersion"]].round(4).to_string(index=False))
    print("\n  se(log v/u) : binomial, bootstrap sur souches seules, "
          "bootstrap hierarchique")
    print(ic[["lignee", "n_souches", "se_log_binom",
              "se_log_boot_souches_seules", "se_log_boot"]]
          .round(4).to_string(index=False))
    print(f"\n  sur-dispersion mediane : x{ic['surdispersion'].median():.2f} "
          f"(min {ic['surdispersion'].min():.2f}, "
          f"max {ic['surdispersion'].max():.2f})")
    print("  l'IC binomial de A15 sous-estime donc l'incertitude d'un facteur "
          f"{ic['surdispersion'].median():.1f} en moyenne.")

    # ---- B. homogeneite du spectre a 6 classes, nul par permutation de souches
    sub = df[df[CLASSES].sum(axis=1) > 0].copy()
    tab = sub.groupby("lignee")[CLASSES].sum()
    chi2_obs = stats.chi2_contingency(tab.to_numpy())[0]
    labels = sub["lignee"].to_numpy()
    mat = sub[CLASSES].to_numpy()
    perm = np.empty(args.boot)
    for b in range(args.boot):
        lab = rng.permutation(labels)
        t = pd.DataFrame(mat).groupby(lab).sum().to_numpy()
        perm[b] = stats.chi2_contingency(t)[0]
    p_perm = (1 + (perm >= chi2_obs).sum()) / (1 + args.boot)
    exp = stats.chi2_contingency(tab.to_numpy())[3]
    resid = (tab.to_numpy() - exp) / np.sqrt(exp)
    r = pd.DataFrame(resid, index=tab.index, columns=CLASSES)
    r.to_csv(args.out_prefix + "homogeneite_spectre.tsv", sep="\t")
    print("\n=== B. homogeneite du spectre a 6 classes ===")
    print(f"  chi2 observe = {chi2_obs:.1f} sur {tab.shape[0]} lignees x 6 "
          f"classes ({(tab.shape[0]-1)*5} ddl)")
    print(f"  p asymptotique = {stats.chi2_contingency(tab.to_numpy())[1]:.3g} "
          f"(NON valide : substitutions groupees par souche)")
    print(f"  p par permutation des etiquettes de lignee entre souches "
          f"({args.boot} tirages) = {p_perm:.4g}")
    print(f"  nul par permutation : chi2 median {np.median(perm):.1f}, "
          f"q99 {np.percentile(perm, 99):.1f}")
    print("\n  residus de Pearson standardises (|r| > 3 : contribution forte)")
    print(r.round(1).to_string())

    # ---- B2. au-dela du v/u : le spectre differe-t-il AUTREMENT que par
    # l'intensite du flux ? On conditionne a la categorie (perte ou gain) et on
    # teste la composition INTERNE : chez les pertes de paires G:C, la part de
    # la transition C>T contre la transversion C>A ; chez les gains, la part de
    # T>C contre T>G. Un ecart ici ne peut pas etre un simple decalage de v/u.
    print("\n=== B2. structure interne du spectre, a categorie fixee ===")
    cond = []
    for nom, cols in [("pertes G:C (C>T vs C>A)", ["C>T", "C>A"]),
                      ("gains G:C (T>C vs T>G)", ["T>C", "T>G"])]:
        sub2 = df[df[cols].sum(axis=1) > 0]
        t = sub2.groupby("lignee")[cols].sum()
        chi_obs = stats.chi2_contingency(t.to_numpy())[0]
        lab0 = sub2["lignee"].to_numpy()
        m2 = sub2[cols].to_numpy()
        pm = np.empty(args.boot)
        for b in range(args.boot):
            lb = rng.permutation(lab0)
            pm[b] = stats.chi2_contingency(
                pd.DataFrame(m2).groupby(lb).sum().to_numpy())[0]
        pp = (1 + (pm >= chi_obs).sum()) / (1 + args.boot)
        frac = (t[cols[0]] / t.sum(axis=1)).sort_values()
        print(f"  {nom} : chi2 = {chi_obs:.1f}, p permutation = {pp:.4g} ; "
              f"part de {cols[0]} de {frac.min():.3f} ({frac.index[0]}) a "
              f"{frac.max():.3f} ({frac.index[-1]})")
        for lin, val in frac.items():
            cond.append(dict(categorie=nom, lignee=lin, part=val,
                             n=int(t.loc[lin].sum()), chi2=chi_obs, p_perm=pp))
    pd.DataFrame(cond).to_csv(args.out_prefix + "structure_conditionnelle.tsv",
                              sep="\t", index=False)

    # ---- C. paires, bootstrap de la difference de log(v/u), FDR
    keys = list(dist.keys())
    pr = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = keys[i], keys[j]
            n = min(len(dist[a]), len(dist[b]))
            d = np.log(dist[a][:n]) - np.log(dist[b][:n])
            p = 2 * min((d <= 0).mean(), (d >= 0).mean())
            p = max(p, 1 / n)
            pr.append(dict(a=a, b=b, vu_a=float(np.median(dist[a])),
                           vu_b=float(np.median(dist[b])),
                           rapport=float(np.median(dist[a]) /
                                         np.median(dist[b])), p=p))
    pdf = pd.DataFrame(pr)
    order = np.argsort(pdf["p"].to_numpy())
    m = len(pdf)
    q = np.empty(m)
    prev = 1.0
    for rank, k in enumerate(order[::-1]):
        val = pdf["p"].to_numpy()[k] * m / (m - rank)
        prev = min(prev, val)
        q[k] = prev
    pdf["q_bh"] = q
    pdf = pdf.sort_values("q_bh")
    pdf.to_csv(args.out_prefix + "paires_fdr.tsv", sep="\t", index=False)
    nsig = int((pdf["q_bh"] < 0.05).sum())
    print(f"\n=== C. comparaisons par paires ({m} paires) ===")
    print(f"  {nsig} paires sur {m} ({100*nsig/m:.0f} %) differentes a "
          f"q_BH < 0,05")
    print("  dix paires les plus separees :")
    print(pdf[["a", "b", "vu_a", "vu_b", "rapport", "q_bh"]].head(10)
          .round(4).to_string(index=False))
    print("\n  paires NON significatives :")
    ns = pdf[pdf["q_bh"] >= 0.05]
    print(ns[["a", "b", "rapport", "p", "q_bh"]].head(15).round(4)
          .to_string(index=False))

    # ---- D. puissance : ecart minimal detectable, lignee par lignee
    z = stats.norm.ppf(0.975) + stats.norm.ppf(0.80)
    ref = ic.loc[ic["se_log_boot"].idxmin()]
    pw = []
    for _, row in ic.iterrows():
        se = np.sqrt(row["se_log_boot"] ** 2 + ref["se_log_boot"] ** 2)
        pw.append(dict(lignee=row["lignee"], n_souches=int(row["n_souches"]),
                       n_evenements=int(row["loss"] + row["gain"]),
                       vu=row["vu"], se_log_boot=row["se_log_boot"],
                       mde_facteur=float(np.exp(z * se)),
                       mde_gc_eq_points=100 * abs(
                           1 / (1 + row["vu"] * np.exp(z * se)) -
                           1 / (1 + row["vu"]))))
    p = pd.DataFrame(pw).sort_values("mde_facteur")
    p.to_csv(args.out_prefix + "puissance.tsv", sep="\t", index=False)
    print(f"\n=== D. puissance (alpha 0,05 bilateral, 1-beta 0,80 ; "
          f"reference = {ref['lignee']}) ===")
    print(p.round(3).to_string(index=False))
    etendue = ic["vu"].max() / ic["vu"].min()
    print(f"\n  etendue inter-lignees observee : facteur {etendue:.2f}")
    for _, row in p.iterrows():
        verdict = ("suffisante" if row["mde_facteur"] < etendue
                   else "INSUFFISANTE")
        if row["mde_facteur"] > 1.5:
            print(f"  {row['lignee']:<12} MDE x{row['mde_facteur']:.2f} : "
                  f"{verdict} pour l'etendue observee")


if __name__ == "__main__":
    main()
