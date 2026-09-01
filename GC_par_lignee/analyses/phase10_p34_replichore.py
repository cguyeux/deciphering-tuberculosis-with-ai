#!/usr/bin/env python3
"""
Objet       : P3.4 -- CONTROLE INTERNE par replichore et par brin de matrice.
              L'architecture replicative (un seul oriC, un seul terminus, pas de
              rearrangement majeur) est identique entre lignees du MTBC : un
              ecart INTER-lignees du flux v/u ne peut donc pas naitre du
              replichore lui-meme. S'il se reproduit a l'identique des deux
              cotes, c'est une replique interne qui renforce le resultat ; s'il
              ne se voit que d'un cote, c'est un artefact positionnel
              (couverture, densite de genes essentiels).

              GEOMETRIE, DERIVEE SANS RIEN EMPRUNTER A LA LITTERATURE POUR LE
              SENS DU BRIN. oriC est ancre sur dnaA (Rv0001, position 1, brin +
              du GFF3) -- convention d'assemblage standard, deja verifiee sur ce
              depot. Le terminus est localise par le skew GC cumule (Lobry
              1996) : extremum de la somme cumulee de (+1 si G, -1 si C) en
              partant d'oriC, loin des bords pour eviter un artefact de bord.
              Le replichore 1 va d'oriC a ter dans le sens des coordonnees
              croissantes ; le replichore 2 le sens inverse. Le brin + (celui
              de la reference H37Rv, sur lequel ref/alt sont toujours exprimes)
              est LU 5'->3' dans le sens des coordonnees croissantes par pure
              convention GenBank -- ceci ne depend d'aucune donnee, seulement de
              la geometrie de la fourche : dans le replichore 1, le brin avance
              (leading) elonge 5'->3' dans le sens croissant, donc son identite
              de sequence est celle du brin + ; son GABARIT (template) est donc
              le brin -, et le brin + porte le gabarit RETARDE (lagging),
              expose plus longtemps en simple brin. Dans le replichore 2 les
              roles s'inversent : le brin + porte le gabarit AVANCE.
              PREDICTION TESTABLE, SANS RIEN SUPPOSER DU SIGNE DU BIAIS. Si la
              desamination de cytosine simple-brin (Castañeda-García et al.
              2020, biais de brin marque chez la mycobacterie ; Niccum et al.
              2024, ces biais positionnels sont normalement corriges par le MMR
              canonique, absent du MTBC) alimente une part du flux pertes
              (G:C -> A:T), alors le taux de pertes dont la base ancestrale
              CORRIGEE est un C (mutant sur le brin +, donc gabarit retarde en
              replichore 1) doit exceder celui dont la base est un G (mutant
              sur le brin -, gabarit avance en replichore 1) -- et ce SENS DOIT
              S'INVERSER dans le replichore 2. Ce test ne suppose rien de plus
              que la geometrie ; il valide ou refute lui-meme le decoupage.
              CRITERES PRE-ENREGISTRES (avant tout regard des resultats).
                C1 REPORTABILITE : une cellule lignee x replichore n'est
                   interpretable en v/u que si elle porte >= 20 pertes ET
                   >= 20 gains (seuil de P4.2 divise par ~1.5, la scission en
                   deux replichores reduisant environ de moitie chaque compte).
                C2 TEST PRIMAIRE (le verdict de P3.4 sur le replichore) : sur
                   les lignees reportables dans les DEUX replichores, rho de
                   Spearman entre v/u(replichore 1) et v/u(replichore 2), et
                   rapport des etendues (ecart type de log v/u, le plus petit
                   sur le plus grand). VERDICT "robuste" si rho >= 0,70 ET
                   rapport >= 0,70 -- exactement les seuils de P4.2/C2, pour
                   rester comparable.
                C3 BRIN DE MATRICE : indice d'asymetrie des pertes
                   AI = ln[(loss_C/n_C) / (loss_G/n_G)] par replichore, sur le
                   pool total (toutes lignees reportables reunies, puissance
                   maximale) et lignee par lignee ou le compte le permet
                   (>= 20 pertes). VERDICT "mecanisme replicatif confirme" si
                   AI change de signe entre replichore 1 et replichore 2.
Entrees     : Canettii/NC_000962.3.fasta (H37Rv, brin +, coordonnees de ref.)
              investigate_phylo/resources/NC_000962.3.gff3 (dnaA, ancrage oriC)
              data/MTBC0/ancestral_on_H37Rv.bin (ancestral, P2.4)
              data/mask_h37rv_positions.npy (masque, P3.1)
              bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt
Sorties     : résultats/phase10_p34_skew.tsv           (oriC, ter, replichores)
              résultats/phase10_p34_vu_replichore.tsv   (par lignee x replichore)
              résultats/phase10_p34_test_replichore.tsv (verdict C2)
              résultats/phase10_p34_brin_matrice.tsv    (indice d'asymetrie C3)
Reutilisable: oui -- localisation oriC/ter par skew GC cumule et decoupage en
              replichores valent pour tout genome bacterien circulaire annote
Projet      : GC_par_lignee
Date        : 2026-08-31
"""
import argparse
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import read_subs, flux  # noqa: E402
from phase2_polarisation_mtbc0 import build_ancestral, read_fasta, H37RV  # noqa: E402
from phase3_sueoka_gc_eq import vu_ci  # noqa: E402
from phase5_p51_panel_recursif import PANEL, strains_of, masked_positions  # noqa: E402
from phase5_p52_bootstrap_puissance import boot_vu  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
MIN_STRAINS = 4
MIN_EVENTS_VU = 20   # C1
MIN_EVENTS_AI = 20   # C3, par classe (pertes)
RHO_SEUIL = 0.70     # C2
RATIO_SEUIL = 0.70   # C2


def locate_ter(seq, oric=0):
    """Terminus par skew GC cumule (Lobry 1996), oriC deja connu (dnaA)."""
    arr = np.frombuffer(seq.encode(), dtype=np.uint8)
    step = np.where(arr == ord("G"), 1, np.where(arr == ord("C"), -1, 0))
    step = step.astype(np.int64)
    rolled = np.roll(step, -oric)
    cum = np.cumsum(rolled)
    n = len(cum)
    lo, hi = int(0.1 * n), int(0.9 * n)
    i_min = lo + int(np.argmin(cum[lo:hi]))
    i_max = lo + int(np.argmax(cum[lo:hi]))
    if abs(cum[i_min]) >= abs(cum[i_max]):
        rel, extreme = i_min, float(cum[i_min])
    else:
        rel, extreme = i_max, float(cum[i_max])
    ter = (oric + rel) % n
    return ter, extreme, cum


def replichore(pos, oric, ter, length):
    d = (pos - oric) % length
    d_ter = (ter - oric) % length
    return 1 if d < d_ter else 2


def opportunity_by_replichore(seq, masked, oric, ter):
    length = len(seq)
    arr = np.frombuffer(seq.encode(), dtype=np.uint8)
    keep = np.ones(length, bool)
    if masked:
        idx = np.fromiter(masked, dtype=np.int64, count=len(masked))
        idx = idx[idx < length]
        keep[idx] = False
    pos = np.arange(length)
    d = (pos - oric) % length
    d_ter = (ter - oric) % length
    repl = np.where(d < d_ter, 1, 2)
    out = {}
    for r in (1, 2):
        sel = keep & (repl == r)
        sub = arr[sel]
        n_c = int((sub == ord("C")).sum())
        n_g = int((sub == ord("G")).sum())
        n_a = int((sub == ord("A")).sum())
        n_t = int((sub == ord("T")).sum())
        out[r] = dict(length=int(sel.sum()), n_c=n_c, n_g=n_g, n_a=n_a, n_t=n_t,
                      n_gc=n_c + n_g, n_at=n_a + n_t)
    return out


def events_par_souche_replichore(strains, anc, masked, oric, ter, length, n,
                                 seed=0):
    """Comptes de branche terminale par souche ET par replichore, avec le
    detail C/G (pertes) et A/T (gains) de la base ANCESTRALE CORRIGEE -- c'est
    elle qui porte l'identite du brin + physiquement mutant, condition du
    controle brin de matrice. Meme mecanique de detection de branche terminale
    (bit unique du masque de support) que phase5_p51.counts_par_souche."""
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
    per = defaultdict(lambda: defaultdict(int))
    for (pos, ref, alt), mask in support.items():
        if mask & (mask - 1) or pos in masked:
            continue
        i = mask.bit_length() - 1
        a = chr(anc[pos]) if pos < len(anc) else "N"
        if a == "N" or a == alt:
            continue
        r = replichore(pos, oric, ter, length)
        c = per[(i, r)]
        f = flux(a, alt)
        if f == "loss":
            c["loss"] += 1
            c["loss_C" if a == "C" else "loss_G"] += 1
        elif f == "gain":
            c["gain"] += 1
            c["gain_A" if a == "A" else "gain_T"] += 1
    rows = []
    for i, name in enumerate(names):
        for r in (1, 2):
            c = per.get((i, r), {})
            rows.append(dict(sra=name, replichore=r,
                             loss=c.get("loss", 0), gain=c.get("gain", 0),
                             loss_C=c.get("loss_C", 0), loss_G=c.get("loss_G", 0),
                             gain_A=c.get("gain_A", 0), gain_T=c.get("gain_T", 0)))
    return pd.DataFrame(rows), k


def asym_index(n_hit, n_opp_hit, n_miss, n_opp_miss):
    """ln[(n_hit/n_opp_hit) / (n_miss/n_opp_miss)], IC95 delta-methode."""
    if n_hit == 0 or n_miss == 0:
        return float("nan"), float("nan"), float("nan")
    ai = np.log((n_hit / n_opp_hit) / (n_miss / n_opp_miss))
    se = np.sqrt(1 / n_hit + 1 / n_miss)
    z = stats.norm.ppf(0.975)
    return ai, ai - z * se, ai + z * se


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--n-per-clade", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-B", "--boot", type=int, default=5000)
    ap.add_argument("--out-prefix",
                    default=str(ROOT / "résultats" / "phase10_p34_"))
    args = ap.parse_args()

    seq = read_fasta(H37RV)
    length = len(seq)
    oric = 0  # dnaA (Rv0001), position 1 (1-based) = index 0
    ter, extreme, cum = locate_ter(seq, oric)
    len1 = (ter - oric) % length
    len2 = length - len1
    print(f"# genome H37Rv NC_000962.3 : {length} pb", file=sys.stderr)
    print(f"# oriC (dnaA, Rv0001) : position 1-based {oric + 1}", file=sys.stderr)
    print(f"# terminus (extremum du skew GC cumule depuis oriC) : "
          f"position 1-based {ter + 1}, skew cumule extreme {extreme:+.0f}",
          file=sys.stderr)
    print(f"# replichore 1 (oriC -> ter, sens croissant) : {len1} pb "
          f"({100 * len1 / length:.1f} %)", file=sys.stderr)
    print(f"# replichore 2 (ter -> oriC, sens croissant) : {len2} pb "
          f"({100 * len2 / length:.1f} %)", file=sys.stderr)

    pd.DataFrame([dict(genome_len=length, oric_1based=oric + 1,
                       ter_1based=ter + 1, skew_extreme=extreme,
                       len_replichore1=len1, len_replichore2=len2,
                       frac_replichore1=len1 / length)]
                ).to_csv(args.out_prefix + "skew.tsv", sep="\t", index=False)

    anc = build_ancestral()
    masked = masked_positions()
    opp = opportunity_by_replichore(seq, masked, oric, ter)
    for r in (1, 2):
        o = opp[r]
        print(f"# replichore {r} hors masque : {o['length']} pb, "
              f"GC = {100 * o['n_gc'] / o['length']:.2f} % "
              f"(C {o['n_c']}, G {o['n_g']}, A {o['n_a']}, T {o['n_t']})",
              file=sys.stderr)

    rng = np.random.default_rng(args.seed)
    vu_rows, cell_dfs = [], {}
    for lin, prefixes in PANEL.items():
        st = strains_of(prefixes)
        if len(st) < MIN_STRAINS:
            continue
        df, k = events_par_souche_replichore(st, anc, masked, oric, ter,
                                             length, args.n_per_clade,
                                             args.seed)
        if df is None:
            continue
        cell_dfs[lin] = df
        for r in (1, 2):
            sub = df[df["replichore"] == r]
            o = opp[r]
            L, G = int(sub["loss"].sum()), int(sub["gain"].sum())
            reportable = L >= MIN_EVENTS_VU and G >= MIN_EVENTS_VU
            row = dict(lignee=lin, replichore=r, n_souches=k, loss=L, gain=G,
                      loss_C=int(sub["loss_C"].sum()),
                      loss_G=int(sub["loss_G"].sum()),
                      gain_A=int(sub["gain_A"].sum()),
                      gain_T=int(sub["gain_T"].sum()),
                      reportable=reportable)
            if reportable:
                vu, lo, hi = vu_ci(L, G, o["n_gc"], o["n_at"])
                bs = boot_vu(sub["loss"].to_numpy(), sub["gain"].to_numpy(),
                            o["n_gc"], o["n_at"], args.boot, rng)
                blo, bhi = (np.percentile(bs, [2.5, 97.5])
                           if len(bs) else (np.nan, np.nan))
                row.update(vu=vu, vu_lo=lo, vu_hi=hi, vu_boot_lo=blo,
                          vu_boot_hi=bhi, gc_eq=1 / (1 + vu))
            else:
                row.update(vu=np.nan, vu_lo=np.nan, vu_hi=np.nan,
                          vu_boot_lo=np.nan, vu_boot_hi=np.nan, gc_eq=np.nan)
            vu_rows.append(row)

    v = pd.DataFrame(vu_rows)
    v.to_csv(args.out_prefix + "vu_replichore.tsv", sep="\t", index=False)
    print("\n=== v/u et GC_eq par lignee x replichore (C1 : "
          f">= {MIN_EVENTS_VU} pertes ET gains) ===")
    print(v[["lignee", "replichore", "n_souches", "loss", "gain", "vu",
             "vu_boot_lo", "vu_boot_hi", "gc_eq",
             "reportable"]].round(4).to_string(index=False))

    # ---- C2 : test primaire, robustesse inter-replichores
    piv = v[v["reportable"]].pivot(index="lignee", columns="replichore",
                                   values="vu")
    piv = piv.dropna()
    if len(piv) >= 4:
        rho, p = stats.spearmanr(piv[1], piv[2])
        s1, s2 = np.std(np.log(piv[1])), np.std(np.log(piv[2]))
        ratio = min(s1, s2) / max(s1, s2)
        verdict = ("ROBUSTE" if rho >= RHO_SEUIL and ratio >= RATIO_SEUIL
                   else "ARTEFACT POSITIONNEL SUSPECTE")
    else:
        rho, p, s1, s2, ratio, verdict = (np.nan,) * 5 + ("SOUS-PUISSANT",)
    test = pd.DataFrame([dict(n_lignees=len(piv), rho=rho, p_rho=p,
                             etendue_log_vu_repl1=s1, etendue_log_vu_repl2=s2,
                             rapport_etendues=ratio, seuil_rho=RHO_SEUIL,
                             seuil_rapport=RATIO_SEUIL, verdict=verdict)])
    test.to_csv(args.out_prefix + "test_replichore.tsv", sep="\t", index=False)
    print(f"\n=== C2 : v/u(replichore 1) vs v/u(replichore 2), "
          f"{len(piv)} lignees reportables dans les deux ===")
    if len(piv):
        print(piv.round(3).to_string())
    print(f"rho de Spearman = {rho:.3f} (p = {p:.3g})"
          if not np.isnan(rho) else "rho : sous-puissant")
    print(f"rapport des etendues (log v/u) = {ratio:.3f}"
          if not np.isnan(ratio) else "rapport des etendues : sous-puissant")
    print(f"VERDICT C2 : {verdict}")

    # ---- C3 : brin de matrice, indice d'asymetrie des pertes C vs G
    ai_rows = []
    for r in (1, 2):
        o = opp[r]
        tot_C = v.loc[v["replichore"] == r, "loss_C"].sum()
        tot_G = v.loc[v["replichore"] == r, "loss_G"].sum()
        ai, lo, hi = asym_index(tot_C, o["n_c"], tot_G, o["n_g"])
        ai_rows.append(dict(lignee="POOL_TOTAL", replichore=r,
                            loss_C=int(tot_C), loss_G=int(tot_G),
                            n_c_site=o["n_c"], n_g_site=o["n_g"],
                            AI_loss=ai, AI_lo=lo, AI_hi=hi))
        for lin in cell_dfs:
            row = v[(v["lignee"] == lin) & (v["replichore"] == r)]
            if row.empty:
                continue
            lc, lg = int(row["loss_C"].iloc[0]), int(row["loss_G"].iloc[0])
            if lc + lg < MIN_EVENTS_AI:
                continue
            ai, lo, hi = asym_index(lc, o["n_c"], lg, o["n_g"])
            ai_rows.append(dict(lignee=lin, replichore=r, loss_C=lc,
                                loss_G=lg, n_c_site=o["n_c"],
                                n_g_site=o["n_g"], AI_loss=ai, AI_lo=lo,
                                AI_hi=hi))
    ai_df = pd.DataFrame(ai_rows)
    ai_df.to_csv(args.out_prefix + "brin_matrice.tsv", sep="\t", index=False)
    print("\n=== C3 : indice d'asymetrie des pertes AI = ln[(C/n_C)/(G/n_G)] "
          f"(>= {MIN_EVENTS_AI} pertes C+G) ===")
    print(ai_df.round(4).to_string(index=False))
    pool = ai_df[ai_df["lignee"] == "POOL_TOTAL"].set_index("replichore")
    if 1 in pool.index and 2 in pool.index:
        ai1, ai2 = pool.loc[1, "AI_loss"], pool.loc[2, "AI_loss"]
        inv = (ai1 * ai2) < 0
        print(f"\nAI(replichore 1) = {ai1:+.4f} [{pool.loc[1,'AI_lo']:+.4f}, "
              f"{pool.loc[1,'AI_hi']:+.4f}]")
        print(f"AI(replichore 2) = {ai2:+.4f} [{pool.loc[2,'AI_lo']:+.4f}, "
              f"{pool.loc[2,'AI_hi']:+.4f}]")
        print(f"VERDICT C3 : {'MECANISME REPLICATIF CONFIRME (signe inverse)' if inv else 'PAS D INVERSION DE SIGNE -- a interpreter avec prudence'}")
        per_lin = ai_df[ai_df["lignee"] != "POOL_TOTAL"]
        if not per_lin.empty:
            pl = per_lin.pivot(index="lignee", columns="replichore",
                               values="AI_loss").dropna()
            if len(pl):
                n_inv = int((pl[1] * pl[2] < 0).sum())
                print(f"  au niveau lignee : {n_inv}/{len(pl)} lignees "
                      "montrent l'inversion de signe attendue")

    print(f"\n# ecrit : {args.out_prefix}"
          "{skew,vu_replichore,test_replichore,brin_matrice}.tsv")


if __name__ == "__main__":
    main()
