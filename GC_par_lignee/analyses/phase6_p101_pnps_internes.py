#!/usr/bin/env python3
"""
Objet       : P10.1 -- chercher un estimateur de Ne INSENSIBLE a la densite
              d'echantillonnage, apres l'effondrement de A20 en P5.3. Le pN/pS
              de A20 etait mesure sur les branches TERMINALES du pool (variants
              portes par une seule souche) ; or la longueur de ces branches est
              une propriete de l'echantillonnage autant que de la biologie, et
              A25 a montre qu'elle covarie avec S (rho = +0,654). Ce script
              mesure le pN/pS sur les branches INTERNES du pool (variants
              partages par 2 a n-1 souches, donc mutations plus anciennes,
              filtrees par la selection sur une duree plus longue), le decompose
              par bande de frequence (profondeur de branche appariee entre
              lignees), et surtout VALIDE ou DISQUALIFIE chaque estimateur par
              deux epreuves INTERNES a chaque lignee, ou le Ne vrai est constant
              par construction et ou seul l'echantillonnage varie.

              CRITERES DE QUALIFICATION, FIXES AVANT DE REGARDER S
              (pre-enregistres ici pour interdire le mecanisme qui a fabrique
              A20 : essayer des proxys jusqu'a ce que l'un rende quelque chose)
                Q1 rarefaction : l'estimateur recalcule sur n' = 10, 20, 30
                   souches tirees du meme pool de 40 ne doit pas deriver de plus
                   que l'ecart interquartile inter-lignees de l'estimateur.
                Q2 dense contre disperse : sur une meme lignee, un pool de 40
                   souches tirees d'UN sous-clade profond (branches terminales
                   courtes) et un pool de 40 souches reparties sur tous les
                   sous-clades (branches longues) doivent donner la meme valeur,
                   a moins d'un ecart interquartile inter-lignees, et sans biais
                   de signe systematique (test des signes).
                Q3 (DIAGNOSTIC, NON QUALIFIANT) correlation inter-lignees avec
                   la longueur de branche terminale. Non qualifiant parce que le
                   Ne vrai et la densite de sequencage sont biologiquement
                   confondus entre lignees : une lignee epidemiologiquement
                   dominante est a la fois plus sequencee et plus grande. Q1 et
                   Q2 tranchent, eux, a Ne constant.
              Un estimateur qui passe Q1 et Q2 est ensuite soumis a UN SEUL test
              de l'hypothese de A20 : Spearman(S, estimateur), rendu quel que
              soit son resultat, sur les dix-sept lignees et sur les treize bien
              pourvues. Aucun autre proxy ne sera essaye dans ce script.
Entrees     : résultats/phase5_p51_counts_par_souche.tsv (pools de P5.1)
              résultats/phase5_p53_force_requise.tsv     (S par lignee, A18)
              bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt
              data/MTBC0/ancestral_on_H37Rv.bin, data/cds_opportunity.json
              investigate_phylo/resources/NC_000962.3.gff3 (CDS H37Rv)
Sorties     : résultats/phase6_p101_classes_support.tsv  (A, une ligne/lignee)
              résultats/phase6_p101_bandes_frequence.tsv (A, lignee x bande)
              résultats/phase6_p101_rarefaction.tsv      (Q1)
              résultats/phase6_p101_dense_disperse.tsv   (Q2)
              résultats/phase6_p101_qualification.tsv    (verdict + test unique)
Reutilisable: oui -- la classe Pool (matrice souche x variant, polarisation,
              effet code, pN/pS corrige du spectre par classe de support) et les
              deux epreuves d'invariance a l'echantillonnage valent pour toute
              bacterie clonale echantillonnee par commodite
Projet      : GC_par_lignee
Date        : 2026-08-30
"""
import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import read_subs  # noqa: E402
from phase2_polarisation_mtbc0 import build_ancestral, read_fasta, H37RV  # noqa: E402
from phase3_sueoka_gc_eq import opportunity  # noqa: E402
from phase4_p93_force_maintien import effect, load_cds  # noqa: E402
from phase5_p51_panel_recursif import (PANEL, masked_positions,  # noqa: E402
                                       pool_dirs, strains_of)

ROOT = Path(__file__).resolve().parent.parent
OPP_CACHE = ROOT / "data" / "cds_opportunity.json"
FAIBLES = {"L8", "Suricattae", "Dassie", "Pinnipedii"}
BASES = "ACGT"
PAIRES = [(a, b) for a in BASES for b in BASES if a != b]
IPAIRE = {p: i for i, p in enumerate(PAIRES)}
# Bandes de frequence intra-pool : « profondeur de branche appariee », lue au
# meme f = k/n chez toutes les lignees. La premiere bande est la branche
# terminale quand n = 40, la derniere la branche menant a la lignee.
# Effet publie en A20 sur les onze pools d'origine, que le test unique doit
# pouvoir exclure ou non : c'est la seule facon de distinguer « rien trouve »
# de « pas assez de puissance pour trouver quoi que ce soit ».
RHO_A20 = -0.736
BANDES = [("f<=0.05", 0.0, 0.05), ("0.05-0.15", 0.05, 0.15),
          ("0.15-0.35", 0.15, 0.35), ("0.35-0.65", 0.35, 0.65),
          ("0.65-0.95", 0.65, 0.95), ("f>0.95", 0.95, 1.0001)]


class Pool:
    """Un pool de souches d'une lignee, sous forme matricielle.

    `M` : matrice booleenne variant x souche des variants non masques (pour pi).
    `A` : sous-matrice des variants POLARISES (allele de H37Rv = etat ancestral
    MTBC0) et CODANTS, avec leur classe de substitution `cls` et leur effet
    `nonsyn`. Toutes les mesures se derivent d'un sous-ensemble de colonnes, ce
    qui rend la rarefaction vectorisable."""

    def __init__(self, subsets, masked, anc, seq, geno, cache):
        owner, off, strand, cds = geno
        n = len(subsets)
        uniq = {}
        for j, s in enumerate(subsets):
            for v in s:
                if v[0] in masked:
                    continue
                uniq.setdefault(v, []).append(j)
        var = list(uniq)
        self.n = n
        self.M = np.zeros((len(var), n), bool)
        for i, v in enumerate(var):
            self.M[i, uniq[v]] = True
        keep, cls, nonsyn = [], [], []
        for i, (pos, ref, alt) in enumerate(var):
            if (chr(anc[pos]) if pos < len(anc) else "N") != ref:
                continue                      # non polarise vers l'ancestral
            key = (pos, alt)
            if key not in cache:
                cache[key] = effect(seq, pos, alt, owner, off, strand, cds)
            eff = cache[key]
            if eff is None:
                continue                      # hors CDS ou codon indecidable
            keep.append(i)
            cls.append(IPAIRE[(ref, alt)])
            nonsyn.append(eff == "nonsyn")
        self.A = self.M[keep]
        self.cls = np.array(cls, np.int8)
        self.nonsyn = np.array(nonsyn, bool)

    def counts(self, idx):
        """(k sur tous les variants, k sur les variants annotes) pour un
        sous-ensemble de souches."""
        return self.M[:, idx].sum(1), self.A[:, idx].sum(1)


def fisher_ci(rho, n, alpha=0.05):
    """IC 95 % d'un rho de Spearman par transformation z de Fisher."""
    if n < 5 or abs(rho) >= 1:
        return np.nan, np.nan
    z, se = np.arctanh(rho), 1 / np.sqrt(n - 3)
    q = stats.norm.ppf(1 - alpha / 2)
    return np.tanh(z - q * se), np.tanh(z + q * se)


def pnps(pool, sel, opp_syn, opp_non):
    """pN/pS corrige du spectre sur le sous-jeu d'evenements `sel` : l'attendu
    N/S est recalcule sous les frequences de classes de substitution observees
    dans CE jeu, ce qui empeche un spectre mutationnel particulier de se faire
    passer pour un regime selectif particulier."""
    cls, ns = pool.cls[sel], pool.nonsyn[sel]
    n_non = int(ns.sum())
    n_syn = int(len(ns) - n_non)
    if not n_syn or not n_non:
        return dict(n_syn=n_syn, n_nonsyn=n_non, pn_ps=np.nan,
                    pn_ps_lo=np.nan, pn_ps_hi=np.nan, attendu_N_sur_S=np.nan)
    obs = np.bincount(cls, minlength=len(PAIRES)).astype(float)
    obs /= obs.sum()
    r = (n_non / n_syn) / ((obs @ opp_non) / (obs @ opp_syn))
    se = np.sqrt(1 / n_non + 1 / n_syn)
    return dict(n_syn=n_syn, n_nonsyn=n_non, pn_ps=r,
                pn_ps_lo=r * np.exp(-1.96 * se), pn_ps_hi=r * np.exp(1.96 * se),
                attendu_N_sur_S=(obs @ opp_non) / (obs @ opp_syn))


def mesure(pool, opp_syn, opp_non, n_sites, idx=None, bandes=False):
    """Toutes les quantites d'un pool (ou d'un sous-echantillon de souches)."""
    if idx is None:
        idx = np.arange(pool.n)
    n = len(idx)
    k_all, k = pool.counts(idx)
    pi_num = float((2 * k_all * (n - k_all)).sum())
    out = dict(n=n, pi=pi_num / (n * (n - 1)) / n_sites,
               branche_term=int((k == 1).sum()) / n,
               n_interne=int(((k >= 2) & (k <= n - 1)).sum()),
               n_fixe=int((k == n).sum()))
    for lab, sel in (("term", k == 1), ("int", (k >= 2) & (k <= n - 1)),
                     ("fixe", k == n)):
        for key, v in pnps(pool, sel, opp_syn, opp_non).items():
            out[f"{key}_{lab}"] = v
    if not bandes:
        return out
    f = k / n
    return out, {lab: pnps(pool, (f > lo) & (f <= hi), opp_syn, opp_non)
                 for lab, lo, hi in BANDES}


def charge(paths, masked, anc, seq, geno, cache):
    subsets = []
    for p in paths:
        v = read_subs(p / "NC_000962.3" / "spdi.txt")
        if v:
            subsets.append(v)
    if len(subsets) < 4:
        return None
    return Pool(subsets, masked, anc, seq, geno, cache)


def pools_dense_disperse(prefixes, n, rng):
    """Deux pools de n souches dans la MEME lignee, a Ne vrai constant :
    « dense » = n souches d'un seul sous-clade profond (branches terminales
    courtes), « disperse » = n souches reparties sur le plus de sous-clades
    possible (branches longues). Rend (dense, disperse, nom_du_sous_clade)."""
    par_dir = {}
    for d in pool_dirs(prefixes):
        s = [x for x in sorted(d.iterdir())
             if x.is_dir() and (x / "NC_000962.3" / "spdi.txt").exists()]
        if s:
            par_dir[d.name] = s
    if len(par_dir) < 8:
        return None
    cands = [k for k, v in par_dir.items() if len(v) >= n]
    if not cands:
        return None
    prof = max(cands, key=lambda k: (k.count("."), len(par_dir[k])))
    dense = rng.sample(par_dir[prof], n)
    noms = sorted(par_dir, key=lambda k: -len(par_dir[k]))
    disperse, tour = [], 0
    while len(disperse) < n:
        pris = False
        for k in noms:
            if tour < len(par_dir[k]):
                disperse.append(par_dir[k][tour])
                pris = True
                if len(disperse) == n:
                    break
        if not pris:
            break
        tour += 1
    if len(disperse) < n:
        return None
    return dense, disperse, prof


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default=str(ROOT / "résultats" /
                                            "phase5_p51_counts_par_souche.tsv"))
    ap.add_argument("--force", default=str(ROOT / "résultats" /
                                           "phase5_p53_force_requise.tsv"))
    ap.add_argument("-B", "--boot", type=int, default=200,
                    help="replicats de rarefaction (Q1)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-prefix", default=str(ROOT / "résultats" /
                                                "phase6_p101_"))
    args = ap.parse_args()

    rng = random.Random(args.seed)
    nprng = np.random.default_rng(args.seed)
    seq = read_fasta(H37RV)
    anc = build_ancestral()
    masked = masked_positions()
    (n_gc, n_at), _ = opportunity()
    n_sites = n_gc + n_at
    geno = load_cds(len(seq))
    raw = json.loads(OPP_CACHE.read_text())
    opp = {tuple(k.split(">")): tuple(v) for k, v in raw["opp"].items()}
    opp_syn = np.array([opp.get(p, (0, 0))[0] for p in PAIRES], float)
    opp_non = np.array([opp.get(p, (0, 0))[1] for p in PAIRES], float)
    print(f"# {len(geno[3])} CDS ; opportunite codante lue du cache "
          f"({opp_syn.sum():.0f} syn, {opp_non.sum():.0f} nonsyn)",
          file=sys.stderr)

    cs = pd.read_csv(args.counts, sep="\t")
    pools_sra = {lin: list(d["sra"]) for lin, d in cs.groupby("lignee")}
    lignees = [l for l in PANEL if l in pools_sra]
    force = pd.read_csv(args.force, sep="\t")[["lignee", "S", "vu", "gc_eq"]]
    cache = {}

    # ---------- A. pN/pS PAR CLASSE DE SUPPORT ---------------------------------
    lignes, bandes_rows, pools = [], [], {}
    for lin in lignees:
        index = {s.name: s for s in strains_of(PANEL[lin])}
        p = charge([index[s] for s in pools_sra[lin] if s in index],
                   masked, anc, seq, geno, cache)
        if p is None:
            continue
        pools[lin] = p
        row, bandes = mesure(p, opp_syn, opp_non, n_sites, bandes=True)
        lignes.append(dict(lignee=lin, **row))
        for lab, v in bandes.items():
            bandes_rows.append(dict(lignee=lin, bande=lab, **v))
        print(f"# {lin} : n = {p.n}, {p.M.shape[0]} variants non masques, "
              f"{p.A.shape[0]} polarises codants", file=sys.stderr)
    ta = pd.DataFrame(lignes).merge(force, on="lignee")
    ta["sous_puissante"] = ta.lignee.isin(FAIBLES)
    ta = ta.sort_values("S", ascending=False)
    ta.to_csv(args.out_prefix + "classes_support.tsv", sep="\t", index=False)
    tb = pd.DataFrame(bandes_rows)
    tb.to_csv(args.out_prefix + "bandes_frequence.tsv", sep="\t", index=False)

    print("\n=== A. pN/pS CORRIGE DU SPECTRE, PAR CLASSE DE SUPPORT ===")
    print("terminal = variants d'une seule souche (mutations jeunes, ce que "
          "mesurait A20) ;\ninterne = 2 a n-1 souches (mutations plus "
          "anciennes) ; fixe = les n souches\n(branche menant a la lignee, "
          "insensible par construction a la densite intra-lignee).")
    print(ta[["lignee", "n", "branche_term", "n_syn_term", "n_nonsyn_term",
              "pn_ps_term", "n_syn_int", "n_nonsyn_int", "pn_ps_int",
              "pn_ps_lo_int", "pn_ps_hi_int", "n_syn_fixe", "n_nonsyn_fixe",
              "pn_ps_fixe", "S"]].round(4).to_string(index=False))
    iqr = {}
    for c in ("pn_ps_term", "pn_ps_int", "pn_ps_fixe", "pi"):
        v = ta[c].dropna()
        iqr[c] = v.quantile(.75) - v.quantile(.25)
        print(f"  {c:12s} : mediane {v.median():.4g}, etendue "
              f"[{v.min():.4g} ; {v.max():.4g}], IQR inter-lignees {iqr[c]:.4g}")

    print("\n  -- pN/pS par bande de frequence f = k/n (profondeur appariee) --")
    piv = tb.pivot(index="lignee", columns="bande",
                   values="pn_ps")[[b[0] for b in BANDES]].reindex(ta.lignee)
    print(piv.round(3).to_string())
    print("  mediane inter-lignees par bande : " +
          "  ".join(f"{b}={piv[b].median():.3f}" for b in piv.columns))

    # ---------- Q1. RAREFACTION INTRA-POOL -------------------------------------
    print("\n=== Q1. RAREFACTION : l'estimateur suit-il l'effectif du pool ? ===")
    print("  A Ne vrai constant (meme lignee, meme pool), on retire des souches.")
    q1 = []
    for lin, p in pools.items():
        if p.n < 30:
            continue
        for m in (10, 20, 30, p.n):
            reps = 1 if m == p.n else args.boot
            vt, vi = [], []
            for _ in range(reps):
                idx = (np.arange(p.n) if reps == 1
                       else nprng.choice(p.n, m, replace=False))
                r = mesure(p, opp_syn, opp_non, n_sites, idx=idx)
                vt.append(r["pn_ps_term"])
                vi.append(r["pn_ps_int"])
            q1.append(dict(lignee=lin, n_sub=m, pn_ps_term=np.nanmedian(vt),
                           pn_ps_int=np.nanmedian(vi)))
    tq1 = pd.DataFrame(q1)
    tq1.to_csv(args.out_prefix + "rarefaction.tsv", sep="\t", index=False)
    print(tq1.pivot(index="lignee", columns="n_sub",
                    values="pn_ps_int").round(3).to_string())
    derive = {}
    for c in ("pn_ps_term", "pn_ps_int"):
        pv = tq1.pivot(index="lignee", columns="n_sub", values=c)
        d = (pv[10] - pv[pv.columns.max()]).abs()
        derive[c] = d.median()
        print(f"  derive mediane |n'=10 - pool complet| de {c} = {d.median():.3f}"
              f" (IQR inter-lignees {iqr[c]:.3f}) -> "
              f"{'PASSE' if d.median() < iqr[c] else 'ECHOUE'} Q1")

    # ---------- Q2. DENSE CONTRE DISPERSE --------------------------------------
    print("\n=== Q2. DENSE CONTRE DISPERSE : l'epreuve decisive ===")
    print("  Deux pools de 40 souches de la MEME lignee, l'un tire d'un seul "
          "sous-clade profond,\n  l'autre reparti sur tous les sous-clades. Le "
          "Ne vrai est le meme ; seule la densite\n  d'echantillonnage change.")
    q2 = []
    for lin in lignees:
        r = pools_dense_disperse(PANEL[lin], 40, rng)
        if r is None:
            continue
        dense, disperse, prof = r
        pd_, ps_ = (charge(dense, masked, anc, seq, geno, cache),
                    charge(disperse, masked, anc, seq, geno, cache))
        if pd_ is None or ps_ is None:
            continue
        md = mesure(pd_, opp_syn, opp_non, n_sites)
        ms = mesure(ps_, opp_syn, opp_non, n_sites)
        q2.append(dict(lignee=lin, sous_clade_dense=prof,
                       **{f"{k}_dense": md[k] for k in
                          ("branche_term", "pn_ps_term", "pn_ps_int", "pi")},
                       **{f"{k}_disperse": ms[k] for k in
                          ("branche_term", "pn_ps_term", "pn_ps_int", "pi")}))
        print(f"# Q2 {lin} : dense {prof} ({md['branche_term']:.1f} subs/souche)"
              f" contre disperse ({ms['branche_term']:.1f})", file=sys.stderr)
    tq2 = pd.DataFrame(q2)
    tq2.to_csv(args.out_prefix + "dense_disperse.tsv", sep="\t", index=False)
    if len(tq2):
        aff = tq2.copy()
        aff["ratio_long"] = aff.branche_term_disperse / aff.branche_term_dense
        print(aff[["lignee", "sous_clade_dense", "branche_term_dense",
                   "branche_term_disperse", "ratio_long", "pn_ps_term_dense",
                   "pn_ps_term_disperse", "pn_ps_int_dense",
                   "pn_ps_int_disperse"]].round(3).to_string(index=False))
        for c, lab in (("pn_ps_term", "pN/pS terminal"),
                       ("pn_ps_int", "pN/pS interne"), ("pi", "diversite pi")):
            d = (tq2[f"{c}_dense"] - tq2[f"{c}_disperse"]).dropna()
            npos = int((d > 0).sum())
            p_signe = stats.binomtest(npos, len(d), 0.5).pvalue
            ok = d.abs().median() < iqr[c] and p_signe > 0.05
            print(f"  {lab} : ecart median |dense - disperse| = "
                  f"{d.abs().median():.4g} (IQR inter-lignees {iqr[c]:.4g}), "
                  f"signes {npos}/{len(d)} p = {p_signe:.3g} -> "
                  f"{'PASSE' if ok else 'ECHOUE'} Q2")
            derive[c + "_q2"] = d.abs().median()
            derive[c + "_q2_p"] = p_signe

    # ---------- Q3. DIAGNOSTIC INTER-LIGNEES (NON QUALIFIANT) -------------------
    print("\n=== Q3 (diagnostic, non qualifiant) : correlation inter-lignees "
          "avec la densite ===")
    bien = ta[~ta.sous_puissante]
    for c in ("pn_ps_term", "pn_ps_int", "pn_ps_fixe", "pi"):
        for t, lab in ((ta, "17"), (bien, "13")):
            u = t.dropna(subset=[c])
            rho, p = stats.spearmanr(u[c], u.branche_term)
            print(f"  Spearman({c:12s}, longueur de branche terminale) n={lab} :"
                  f" {rho:+.3f} (p = {p:.3g})")
    print("  Non qualifiant : entre lignees, Ne vrai et densite de sequencage "
          "sont confondus\n  (une lignee dominante est a la fois plus grande et "
          "plus sequencee). Q1 et Q2\n  tranchent a Ne constant, ce que cette "
          "correlation ne peut pas faire.")

    # ---------- TEST UNIQUE PRE-ENREGISTRE --------------------------------------
    print("\n=== TEST UNIQUE de l'hypothese de A20 avec l'estimateur qualifie ===")
    print("  Attendu si S = 2.Ne.s a s constant : correlation NEGATIVE entre S "
          "et pN/pS.")
    verdict = []
    for c, lab in (("pn_ps_int", "pN/pS interne"),
                   ("pn_ps_term", "pN/pS terminal (rappel A20)")):
        for t, etiq in ((ta, "dix-sept lignees"),
                        (bien, "treize bien pourvues")):
            u = t.dropna(subset=[c])
            rho, p = stats.spearmanr(u.S, u[c])
            lo, hi = fisher_ci(rho, len(u))
            print(f"  Spearman(S, {lab}) sur {etiq} (n = {len(u)}) : "
                  f"{rho:+.3f} (p = {p:.3g}), IC 95 % [{lo:+.3f} ; {hi:+.3f}]"
                  f" -> {'exclut' if not lo <= RHO_A20 <= hi else 'compatible '
                        'avec'} le rho = {RHO_A20:+.3f} de A20")
            verdict.append(dict(estimateur=c, jeu=etiq, n=len(u), rho=rho, p=p,
                                rho_ic_lo=lo, rho_ic_hi=hi,
                                exclut_effet_a20=not lo <= RHO_A20 <= hi,
                                derive_q1=derive.get(c),
                                ecart_q2=derive.get(c + "_q2"),
                                p_signes_q2=derive.get(c + "_q2_p"),
                                iqr_inter_lignees=iqr.get(c)))
    pd.DataFrame(verdict).to_csv(args.out_prefix + "qualification.tsv",
                                 sep="\t", index=False)
    print("\n  Aucun autre proxy n'est essaye ici : P10.3 interdit d'en chercher "
          "un troisieme\n  pour sauver l'hypothese.")


if __name__ == "__main__":
    main()
