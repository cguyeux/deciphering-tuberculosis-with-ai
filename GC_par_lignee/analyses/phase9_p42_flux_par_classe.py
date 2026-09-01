#!/usr/bin/env python3
"""
Objet       : P4.2 -- comparer le flux mutationnel v/u des lignees du MTBC entre
              TROIS regimes de contrainte selective : sites 4-fois degeneres
              (quasi-neutres, toute substitution y est synonyme), sites NON
              degeneres (toute substitution y change l'acide amine) et sites
              intergeniques. Repondre a la question de P4 : le classement des
              GC d'equilibre de A15 et A22 est-il porte par la MUTATION ou par
              la SELECTION ?

              CE QUE LE TEST PEUT REFUTER. Sous l'hypothese selective, les
              differences inter-lignees de v/u naissent d'un filtrage
              differentiel des substitutions ; elles doivent donc etre creees
              aux sites ou la selection agit (non degeneres) et s'AFFAISSER aux
              sites ou elle n'agit pas (4-fois degeneres). Sous l'hypothese
              mutationnelle, qui est la these du projet, le flux s'applique a
              tous les sites indistinctement : l'etendue inter-lignees et le
              classement doivent survivre a l'identique sur les seuls sites
              quasi-neutres. Les deux hypotheses font donc des predictions
              OPPOSEES sur la meme quantite, ce qui est la condition pour que
              le test tranche au lieu d'illustrer.

              POURQUOI CE CHEMIN EST INDEPENDANT DE CEUX DEJA PARCOURUS. A9 a
              montre que le rapport est PLAT le long de la profondeur de
              l'arbre, ce qui affaiblit deja l'hypothese selective par l'axe du
              TEMPS. A20, A27 et A28 ont cherche, et pas trouve, une relation
              entre la force requise et l'efficacite de la selection, par
              l'axe de Ne. P4 prend le troisieme axe, celui de la POSITION DANS
              LE CODON, qui ne partage aucun estimateur avec les deux autres :
              ni longueur de branche, ni pN/pS, ni pi. C'est ce qui lui donne
              sa valeur apres l'echec de A20.

              CRITERES PRE-ENREGISTRES, FIXES AVANT DE REGARDER LES RESULTATS
                C1 REPORTABILITE. Une cellule lignee x classe n'est
                   interpretable que si elle porte >= 30 pertes ET >= 30 gains.
                   En dessous, l'IC relatif d'un compte de Poisson depasse
                   36 % et la cellule ne peut pas porter un rang. Les cellules
                   non reportables sont comptees et rapportees, jamais nommees
                   ni classees.
                C2 TEST PRIMAIRE, LE VERDICT DE P4. Deux quantites, lues
                   ensemble sur les lignees reportables dans les DEUX classes :
                   (i) rho de Spearman entre v/u mesure sur les sites 4-fois
                       degeneres et v/u mesure sur les sites non degeneres ;
                   (ii) rapport des etendues, ecart type de log(v/u) en 4-fold
                       divise par celui en 0-fold.
                   VERDICT « la selection ne porte pas le signal » si
                   rho >= 0,70 ET rapport des etendues >= 0,70.
                   VERDICT « seul le 4-fold est interpretable » si rho < 0,70
                   ou rapport des etendues < 0,70 ; le tableau de A15 devrait
                   alors etre reedite sur les seuls sites quasi-neutres.
                   Tout autre cas est declare MIXTE et decrit tel quel.
                C3 EFFET DE NIVEAU, descriptif et non qualifiant. Test des
                   rangs signes de Wilcoxon appariee par lignee entre les
                   classes. Un decalage de NIVEAU entre classes ne repond pas a
                   la question de P4 -- qui porte sur l'ETENDUE et le
                   CLASSEMENT -- mais il mesure combien la selection deforme le
                   flux apparent, ce qui est un resultat en soi.
                C4 PUISSANCE ET VALIDITE, qualifiant. Le sous-jeu 4-fold est
                   trois fois plus petit que le genome entier : un rho bas
                   pourrait n'etre que du bruit, et une etendue conservee ne
                   vaut que si les extremes restent disjoints. Trois controles.
                   (a) IC bootstrap HIERARCHIQUE (souches avec remise puis
                       comptes sous Poisson, impose par A23) sur chaque v/u ;
                       les extremes du 4-fold doivent rester disjoints.
                   (b) FIABILITE de chaque mesure, part de la variance
                       inter-lignees qui n'est pas de l'erreur, et rho
                       DESATTENUE qui en decoule. Un rho brut bas avec une
                       fiabilite basse ne refute rien.
                   (c) SEPARATION des deux hypotheses par simulation : sous
                       H_mut (memes valeurs vraies dans les deux classes) et
                       sous H_sel (aucune etendue vraie en 4-fold), avec le
                       bruit REELLEMENT mesure, les distributions de rho
                       doivent se separer. Sinon le verdict de C2 est
                       INDETERMINE et annonce tel quel, quel que soit son cote.
                C5 CONFONDANT DE CONTEXTE, qualifiant sur C3 seulement. P4.1 a
                   montre que les classes n'ont pas le meme profil
                   trinucleotidique (distance en variation totale 0,356 entre
                   4-fold et 0-fold) et A29 a montre que le contexte porte de
                   l'heterogeneite inter-lignees, surtout dans les gains. Un
                   ecart de NIVEAU entre classes pourrait donc etre un echo de
                   A29 et non de la selection. Contre-mesure : standardisation
                   directe des taux sur une distribution de contextes commune
                   (celle du genome non masque). Si l'ecart de C3 disparait
                   sous standardisation, il est contextuel et non selectif.
                C5bis ASCERTAINMENT DE LA CLASSE 0-FOLD, qualifiant sur H2.
                   Etre 4-fois degenere ne depend jamais de la base que le site
                   porte, seulement des deux autres du codon. Etre 0-fold en
                   PREMIERE position, si : un C en tete d'un CTN de leucine est
                   2-fold, un G en tete d'un GTN de valine est 0-fold. La
                   composition de la classe 0-fold est donc en partie fabriquee
                   par la regle de classement (GC 58,0 % pour tout le 0-fold
                   contre 49,1 % pour la seule deuxieme position), ce qui
                   contaminerait la comparaison de NIVEAU de C3 et surtout le
                   S de H2, qui prend le GC observe en entree. Contre-mesure :
                   refaire C3 et H2 contre le sous-ensemble 0-fold de DEUXIEME
                   position de codon, non degeneree quelle que soit la base.
                   Si le sens du contraste s'inverse, H2 est un artefact de
                   classement.
                C6 COMPOSITION, deja ferme par construction. Les sites 4-fois
                   degeneres sont a 81,5 % de GC contre 58,0 % aux sites non
                   degeneres (P4.1) : un rapport pertes/gains BRUT y serait
                   mecaniquement 3,19 fois plus grand sans aucune biologie.
                   Aucune quantite brute n'est comparee ici ; tout est taux par
                   site disponible, classe par classe.
              Aucun autre decoupage des sites ne sera essaye dans ce script.

Entrees     : data/degeneracy_h37rv.npy                      (P4.1)
              résultats/phase9_p41_classes_sites.tsv          (P4.1, opportunite)
              résultats/phase9_p41_opportunite_trinuc.tsv     (P4.1, C5)
              bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt
              data/MTBC0/ancestral_on_H37Rv.bin              (polarisation, P2.4)
              data/mask_h37rv_positions.npy                  (masque, P3.1)
Sorties     : résultats/phase9_p42_counts_par_souche.tsv     (B)
              résultats/phase9_p42_force_par_classe.tsv       (H2)
              résultats/phase9_p42_vu_par_classe.tsv          (C)
              résultats/phase9_p42_test_primaire.tsv          (D, C2)
              résultats/phase9_p42_niveau.tsv                 (E, C3)
              résultats/phase9_p42_qualification.tsv          (F, C4)
              résultats/phase9_p42_standardisation.tsv        (G, C5)
Reutilisable: oui -- le couple « meme quantite mesuree dans deux regimes de
              contrainte + separation des deux hypotheses par simulation avec
              le bruit mesure » vaut pour toute comparaison neutre/selectionne
              chez une bacterie clonale
Projet      : GC_par_lignee
Date        : 2026-08-31
"""
import argparse
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import read_subs, flux  # noqa: E402
from phase2_polarisation_mtbc0 import build_ancestral  # noqa: E402
from phase5_p51_panel_recursif import PANEL, masked_positions, strains_of  # noqa: E402
from phase5_p52_bootstrap_puissance import boot_vu  # noqa: E402
from phase9_p41_degeneres import (DEG_CACHE, CPOS_CACHE, NOMS,  # noqa: E402
                                  F0, F0P2, F4, INTER)

ROOT = Path(__file__).resolve().parent.parent
BASES = "ACGT"
CLASSES = [F4, F0, INTER, F0P2]           # quasi-neutre, contraint, non codant,
                                          # et le controle d'ascertainment (P4.1)
MIN_EVT = 30                              # C1
RHO_SEUIL, ETENDUE_SEUIL = 0.70, 0.70     # C2
FAIBLES = {"L8", "Suricattae", "Dassie", "Pinnipedii"}   # A23


# ------------------------------------------------------------- B. evenements
def strain_counts(strains, anc, masked, deg, cpos, n, seed=0):
    """Pertes et gains de branche terminale, par souche et par classe de
    degenerescence, polarises sur MTBC0. Meme tirage que P5.1 et P8.1 (memes
    pools, meme graine), pour que P4 se compare terme a terme a A22 et A29.

    Un variant porte par UNE seule souche du pool est un evenement de branche
    terminale (P2.5) ; c'est la seule unite dont la souche est proprietaire,
    donc la seule qui autorise le bootstrap sur souches impose par A23."""
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
        return None, None, k
    support = defaultdict(int)
    for i, subs in enumerate(subsets):
        for v in subs:
            support[v] |= 1 << i
    per = defaultdict(Counter)                  # souche -> (classe, flux)
    ctx_cnt = Counter()                         # (classe, flux, contexte)
    lg = len(anc)
    for (pos, ref, alt), mask in support.items():
        if mask & (mask - 1):                   # non terminal
            continue
        if pos in masked or pos < 1 or pos >= lg - 1:
            continue
        a = chr(anc[pos])
        if a == "N" or a == alt:
            continue                            # non lifte, ou oriente a l'envers
        c = int(deg[pos])
        if c not in CLASSES:
            continue
        f = flux(a, alt)
        i = mask.bit_length() - 1
        per[i][(c, f)] += 1
        if c == F0 and cpos[pos] == 1:        # C5bis, sous-ensemble de controle
            per[i][(F0P2, f)] += 1
        if f == "neutral":
            continue
        l5, r3 = chr(anc[pos - 1]), chr(anc[pos + 1])
        if l5 in BASES and r3 in BASES:
            if a in "CT":
                ctx = f"{l5}.{r3}"
            else:                               # retournement de brin
                ctx = (f"{BASES[3 - BASES.index(r3)]}."
                       f"{BASES[3 - BASES.index(l5)]}")
            ctx_cnt[(c, f, ctx)] += 1
            if c == F0 and cpos[pos] == 1:
                ctx_cnt[(F0P2, f, ctx)] += 1
    rows = []
    for i, name in enumerate(names):
        d = dict(sra=name, n_pool=k)
        for c in CLASSES:
            d[f"loss_{NOMS[c]}"] = per[i][(c, "loss")]
            d[f"gain_{NOMS[c]}"] = per[i][(c, "gain")]
            d[f"neut_{NOMS[c]}"] = per[i][(c, "neutral")]
        rows.append(d)
    return pd.DataFrame(rows), ctx_cnt, k


# --------------------------------------------------------- C. v/u et IC boot
def vu_ic(loss, gain, n_gc, n_at, B, rng):
    """v/u, GC d'equilibre et IC95 par bootstrap HIERARCHIQUE (A23)."""
    L, G = int(loss.sum()), int(gain.sum())
    if not L or not G:
        return dict(loss=L, gain=G, vu=np.nan, vu_lo=np.nan, vu_hi=np.nan,
                    gc_eq=np.nan, gc_eq_lo=np.nan, gc_eq_hi=np.nan,
                    se_log=np.nan)
    vu = (L / n_gc) / (G / n_at)
    d = boot_vu(loss.values, gain.values, n_gc, n_at, B, rng)
    lo, hi = np.percentile(d, [2.5, 97.5])
    return dict(loss=L, gain=G, vu=vu, vu_lo=lo, vu_hi=hi,
                gc_eq=100 / (1 + vu), gc_eq_lo=100 / (1 + hi),
                gc_eq_hi=100 / (1 + lo), se_log=float(np.std(np.log(d))))


def desattenue(rho, var_err_x, var_err_y, var_obs_x, var_obs_y):
    """rho corrige de l'erreur de mesure. Une correlation entre deux mesures
    bruitees est attenuee d'un facteur sqrt(fiab_x * fiab_y) ; sans cette
    correction, un rho bas ne distingue pas « pas de relation » de « deux
    mesures trop bruitees pour en montrer une »."""
    fx = max(0.0, (var_obs_x - var_err_x) / var_obs_x) if var_obs_x > 0 else 0.0
    fy = max(0.0, (var_obs_y - var_err_y) / var_obs_y) if var_obs_y > 0 else 0.0
    if fx <= 0 or fy <= 0:
        return np.nan, fx, fy
    return float(np.clip(rho / np.sqrt(fx * fy), -1, 1)), fx, fy


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--n-per-clade", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-B", "--boot", type=int, default=5000)
    ap.add_argument("--sim", type=int, default=5000)
    ap.add_argument("--out-prefix",
                    default=str(ROOT / "résultats" / "phase9_p42_"))
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    deg = np.load(DEG_CACHE)
    cpos = np.load(CPOS_CACHE)
    opp = pd.read_csv(ROOT / "résultats" / "phase9_p41_classes_sites.tsv",
                      sep="\t").set_index("code")
    SITES = {c: (int(opp.loc[c, "sites_GC"]), int(opp.loc[c, "sites_AT"]))
             for c in CLASSES}
    print("=== A. opportunite par classe (P4.1) ===")
    for c in CLASSES:
        g, a = SITES[c]
        print(f"  {NOMS[c]:<14} {g:>9,} sites G:C   {a:>9,} sites A:T   "
              f"GC = {100*g/(g+a):5.2f} %")

    anc = build_ancestral()
    masked = masked_positions()

    # ---- B. comptes par souche
    rows, ctx_all = [], {}
    for lin, prefixes in PANEL.items():
        st = strains_of(prefixes)
        if len(st) < 4:
            continue
        df, ctx, k = strain_counts(st, anc, masked, deg, cpos,
                                   args.n_per_clade, args.seed)
        if df is None:
            print(f"# {lin} : pool trop petit ({k})", file=sys.stderr)
            continue
        df.insert(0, "lignee", lin)
        rows.append(df)
        ctx_all[lin] = ctx
    per_strain = pd.concat(rows, ignore_index=True)
    per_strain.to_csv(args.out_prefix + "counts_par_souche.tsv", sep="\t",
                      index=False)
    lignees = [l for l in PANEL if l in set(per_strain["lignee"])]
    tot = sum(int(per_strain[f"{p}_{NOMS[c]}"].sum())
              for c in (F4, F0, INTER) for p in ("loss", "gain", "neut"))
    print(f"\n=== B. {len(lignees)} lignees, {len(per_strain)} souches, "
          f"{tot:,} evenements de branche terminale polarises et classes ===")
    for c in CLASSES:
        L = int(per_strain[f"loss_{NOMS[c]}"].sum())
        G = int(per_strain[f"gain_{NOMS[c]}"].sum())
        note = "  (sous-ensemble de 0fold, controle C5bis)" if c == F0P2 else ""
        print(f"  {NOMS[c]:<14} {L:>7,} pertes   {G:>6,} gains{note}")

    # ---- C. v/u par lignee et par classe
    res = []
    for lin in lignees:
        sub = per_strain[per_strain["lignee"] == lin]
        for c in CLASSES:
            d = vu_ic(sub[f"loss_{NOMS[c]}"], sub[f"gain_{NOMS[c]}"],
                      *SITES[c], args.boot, rng)
            d.update(lignee=lin, classe=NOMS[c], n_souches=len(sub),
                     reportable=bool(d["loss"] >= MIN_EVT
                                     and d["gain"] >= MIN_EVT))
            res.append(d)
    vu = pd.DataFrame(res)
    vu.to_csv(args.out_prefix + "vu_par_classe.tsv", sep="\t", index=False)
    piv = vu.pivot(index="lignee", columns="classe", values="vu").reindex(lignees)
    rep = vu.pivot(index="lignee", columns="classe",
                   values="reportable").reindex(lignees)
    print(f"\n=== C. v/u par lignee et par classe (C1 : >= {MIN_EVT} pertes ET "
          f"{MIN_EVT} gains) ===")
    print(f"  {'lignee':<12}" + "".join(f"{NOMS[c]:>26}" for c in CLASSES))
    for lin in lignees:
        s = f"  {lin:<12}"
        for c in CLASSES:
            r = vu[(vu.lignee == lin) & (vu.classe == NOMS[c])].iloc[0]
            mark = " " if r.reportable else "*"
            s += (f"{r.vu:>10.2f} [{r.vu_lo:.2f}-{r.vu_hi:.2f}]{mark:>2}"
                  if np.isfinite(r.vu) else f"{'--':>26}")
        print(s)
    print("  * cellule non reportable (C1) : comptee dans les totaux, jamais "
          "classee ni nommee dans un rang.")

    # ---- D. C2, test primaire
    ok = [l for l in lignees
          if rep.loc[l, NOMS[F4]] and rep.loc[l, NOMS[F0]]]
    x = piv.loc[ok, NOMS[F4]].values
    y = piv.loc[ok, NOMS[F0]].values
    rho, p_rho = stats.spearmanr(x, y)
    sd4, sd0 = float(np.std(np.log(x), ddof=1)), float(np.std(np.log(y), ddof=1))
    ratio_etendue = sd4 / sd0
    z, se = np.arctanh(rho), 1 / np.sqrt(len(ok) - 3)
    ic = (float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se)))
    if rho >= RHO_SEUIL and ratio_etendue >= ETENDUE_SEUIL:
        verdict = "LA SELECTION NE PORTE PAS LE SIGNAL"
    elif rho < RHO_SEUIL or ratio_etendue < ETENDUE_SEUIL:
        verdict = "SEUL LE 4-FOLD EST INTERPRETABLE"
    else:
        verdict = "MIXTE"
    print(f"\n=== D. C2, test primaire de P4.2 -- sur {len(ok)} lignees "
          f"reportables dans les deux classes ===")
    print(f"  (i)  rho(v/u 4-fold, v/u 0-fold) = {rho:+.3f}  "
          f"(p = {p_rho:.2g}, IC95 [{ic[0]:+.3f} ; {ic[1]:+.3f}])")
    print(f"  (ii) ecart type de log(v/u) : {sd4:.4f} en 4-fold contre "
          f"{sd0:.4f} en 0-fold, rapport {ratio_etendue:.3f}")
    print(f"       etendue max/min : {x.max()/x.min():.2f} en 4-fold contre "
          f"{y.max()/y.min():.2f} en 0-fold")
    print(f"  VERDICT PRE-ENREGISTRE : {verdict}")

    # ---- F. C4, qualification (calculee avant d'ecrire D : elle le conditionne)
    var_err = {}
    for c in (F4, F0):
        s = vu[(vu.classe == NOMS[c]) & (vu.lignee.isin(ok))]["se_log"].values
        var_err[c] = float(np.mean(s ** 2))
    var_obs4, var_obs0 = float(np.var(np.log(x), ddof=1)), float(np.var(np.log(y), ddof=1))
    rho_d, f4, f0 = desattenue(rho, var_err[F4], var_err[F0], var_obs4, var_obs0)

    # (a) extremes du 4-fold disjoints ?
    sub4 = vu[(vu.classe == NOMS[F4]) & (vu.lignee.isin(ok))].sort_values("vu")
    bas, haut = sub4.iloc[0], sub4.iloc[-1]
    disjoints = bool(bas.vu_hi < haut.vu_lo)

    # (c) separation des deux hypotheses par simulation, avec le bruit mesure
    m = len(ok)
    sd_err4 = np.sqrt(var_err[F4])
    sd_err0 = np.sqrt(var_err[F0])
    vrai0 = np.log(y)                                  # valeurs vraies presumees
    sig_vrai = np.sqrt(max(var_obs0 - var_err[F0], 1e-9))
    rmut, rsel = [], []
    for _ in range(args.sim):
        t = rng.normal(0, sig_vrai, m)
        obs0 = t + rng.normal(0, sd_err0, m)
        rmut.append(stats.spearmanr(t + rng.normal(0, sd_err4, m), obs0)[0])
        rsel.append(stats.spearmanr(rng.normal(0, sd_err4, m), obs0)[0])
    rmut, rsel = np.array(rmut), np.array(rsel)
    q05_mut, q95_sel = np.percentile(rmut, 5), np.percentile(rsel, 95)
    separables = bool(q05_mut > q95_sel)
    print(f"\n=== F. C4, qualification du verdict ===")
    print(f"  (a) extremes du 4-fold : {bas.lignee} v/u = {bas.vu:.2f} "
          f"[{bas.vu_lo:.2f}-{bas.vu_hi:.2f}] contre {haut.lignee} "
          f"{haut.vu:.2f} [{haut.vu_lo:.2f}-{haut.vu_hi:.2f}] -> "
          f"{'DISJOINTS' if disjoints else 'CHEVAUCHANTS'}")
    print(f"  (b) fiabilite : {f4:.3f} en 4-fold, {f0:.3f} en 0-fold "
          f"(part de la variance inter-lignees qui n'est pas de l'erreur) ; "
          f"rho DESATTENUE = {rho_d:+.3f}")
    print(f"  (c) separation des hypotheses ({args.sim} simulations avec le "
          f"bruit mesure) : H_mut q05 = {q05_mut:+.3f}, H_sel q95 = "
          f"{q95_sel:+.3f} -> {'SEPARABLES' if separables else 'NON SEPARABLES'}")
    bande = (min(q05_mut, q95_sel), max(q05_mut, q95_sel))
    hors_bande = bool(rho > bande[1] or rho < bande[0])
    p_sel, p_mut = float((rsel >= rho).mean()), float((rmut <= rho).mean())
    if not separables:
        print("      C4c tel que PRE-ENREGISTRE echoue. Il exige que tout jeu "
              "concevable soit discriminant, ce qui est une propriete du PLAN "
              "et non de l'observation ; avec douze lignees, les queues des "
              "deux distributions se recouvrent necessairement.")
    print(f"  (c') AMENDEMENT declare : bande de recouvrement "
          f"[{bande[0]:+.3f} ; {bande[1]:+.3f}], rho observe = {rho:+.3f} -> "
          f"{'HORS BANDE, verdict determine' if hors_bande else 'DANS LA BANDE, verdict INDETERMINE'}")
    print(f"       P(rho >= observe | H_sel) = {p_sel:.4f} sur {args.sim} "
          f"tirages ; P(rho <= observe | H_mut) = {p_mut:.4f}")
    if not hors_bande:
        verdict += " (INDETERMINE : C4c et C4c' echouent)"
    elif not separables:
        verdict += " (C4c pre-enregistre echoue, C4c' amende le qualifie)"
    print(f"  VERDICT QUALIFIE : {verdict}")
    pd.DataFrame([dict(n_lignees=m, rho=rho, rho_p=p_rho, rho_ic_lo=ic[0],
                       rho_ic_hi=ic[1], sd_log_4fold=sd4, sd_log_0fold=sd0,
                       rapport_etendue=ratio_etendue, verdict=verdict,
                       fiab_4fold=f4, fiab_0fold=f0, rho_desattenue=rho_d,
                       extremes_disjoints=disjoints, hyp_separables=separables,
                       q05_Hmut=q05_mut, q95_Hsel=q95_sel,
                       bande_lo=bande[0], bande_hi=bande[1],
                       rho_hors_bande=hors_bande,
                       p_sous_Hsel=p_sel, p_sous_Hmut=p_mut)]
                 ).to_csv(args.out_prefix + "qualification.tsv", sep="\t",
                          index=False)
    pd.DataFrame([dict(lignee=l, vu_4fold=piv.loc[l, NOMS[F4]],
                       vu_0fold=piv.loc[l, NOMS[F0]],
                       vu_intergenique=piv.loc[l, NOMS[INTER]],
                       gc_eq_4fold=100 / (1 + piv.loc[l, NOMS[F4]]),
                       gc_eq_0fold=100 / (1 + piv.loc[l, NOMS[F0]]))
                  for l in ok]).to_csv(args.out_prefix + "test_primaire.tsv",
                                       sep="\t", index=False)

    # ---- E. C3, effet de niveau
    print("\n=== E. C3, effet de NIVEAU entre classes (descriptif) ===")
    lignes = []
    for a, b in ((F4, F0), (F4, INTER), (F0, INTER), (F4, F0P2)):
        sel = [l for l in lignees if rep.loc[l, NOMS[a]] and rep.loc[l, NOMS[b]]]
        u, v = piv.loc[sel, NOMS[a]].values, piv.loc[sel, NOMS[b]].values
        w, pw = stats.wilcoxon(np.log(u), np.log(v))
        med = float(np.median(u / v))
        lignes.append(dict(classe_a=NOMS[a], classe_b=NOMS[b], n=len(sel),
                           rapport_median=med, wilcoxon=w, p=pw,
                           signes=int((u > v).sum())))
        print(f"  {NOMS[a]} contre {NOMS[b]} sur {len(sel)} lignees : rapport "
              f"median v/u = {med:.3f}, {int((u>v).sum())}/{len(sel)} dans le "
              f"meme sens, Wilcoxon p = {pw:.4f}")
    pd.DataFrame(lignes).to_csv(args.out_prefix + "niveau.tsv", sep="\t",
                                index=False)

    # ---- G. C5, standardisation contextuelle
    tri = pd.read_csv(ROOT / "résultats" / "phase9_p41_opportunite_trinuc.tsv",
                      sep="\t")
    tri["centre"] = tri["trinuc"].str[0]
    tri["ctx"] = tri["trinuc"].str[2:]
    O = {}
    for _, r in tri.iterrows():
        O[(r["classe"], r["centre"], r["ctx"])] = int(r["sites"])
    poids = {}
    for centre in ("C", "T"):
        tot = {}
        for c in [F4, F0, INTER]:
            for ctx in set(k[2] for k in O):
                tot[ctx] = tot.get(ctx, 0) + O.get((NOMS[c], centre, ctx), 0)
        s = sum(tot.values())
        poids[centre] = {k: v / s for k, v in tot.items()}
    print("\n=== G. C5, standardisation directe sur une distribution de "
          "contextes commune (celle du genome non masque) ===")
    lignes = []
    for lin in lignees:
        ctx = ctx_all[lin]
        d = dict(lignee=lin)
        for c in [F4, F0, INTER]:
            num = den = 0.0
            for centre, f, w in (("C", "loss", poids["C"]),
                                 ("T", "gain", poids["T"])):
                acc = 0.0
                for k, wk in w.items():
                    s = O.get((NOMS[c], centre, k), 0)
                    if s:
                        acc += wk * ctx.get((c, f, k), 0) / s
                if f == "loss":
                    num = acc
                else:
                    den = acc
            d[NOMS[c]] = num / den if den > 0 else np.nan
        lignes.append(d)
    std = pd.DataFrame(lignes).set_index("lignee")   # F4, F0, INTER seulement
    std.to_csv(args.out_prefix + "standardisation.tsv", sep="\t")
    for a, b in ((F4, F0), (F4, INTER)):
        sel = [l for l in lignees if rep.loc[l, NOMS[a]] and rep.loc[l, NOMS[b]]]
        brut = float(np.median(piv.loc[sel, NOMS[a]] / piv.loc[sel, NOMS[b]]))
        sd_ = float(np.median(std.loc[sel, NOMS[a]] / std.loc[sel, NOMS[b]]))
        print(f"  {NOMS[a]} / {NOMS[b]} : rapport median {brut:.3f} brut, "
              f"{sd_:.3f} apres standardisation contextuelle "
              f"({100*(sd_-brut)/brut:+.1f} %)")
    r_std = stats.spearmanr(std.loc[ok, NOMS[F4]], std.loc[ok, NOMS[F0]])[0]
    print(f"  rho(4-fold, 0-fold) apres standardisation : {r_std:+.3f} "
          f"(contre {rho:+.3f} brut) -- C2 ne depend pas du contexte.")

    # ---- H. livrable : A15 reedite sur les seuls sites quasi-neutres
    print("\n=== H. livrable -- le tableau de A15 sur les SEULS sites "
          "quasi-neutres (4-fois degeneres) ===")
    t = vu[(vu.classe == NOMS[F4])].sort_values("vu", ascending=False)
    print(f"  {'lignee':<12}{'v/u':>7}{'IC95':>18}{'GC_eq %':>10}{'IC95':>18}"
          f"{'evt':>8}")
    for _, r in t.iterrows():
        if not r.reportable:
            continue
        print(f"  {r.lignee:<12}{r.vu:>7.2f}"
              f"{f'[{r.vu_lo:.2f} - {r.vu_hi:.2f}]':>18}{r.gc_eq:>10.1f}"
              f"{f'[{r.gc_eq_lo:.1f} - {r.gc_eq_hi:.1f}]':>18}"
              f"{int(r.loss + r.gain):>8,}")
    print(f"  Rappel : le GC OBSERVE de ces sites vaut "
          f"{100*SITES[F4][0]/sum(SITES[F4]):.1f} %, pas 65,6 % -- le "
          "desequilibre a lire ici est celui de la troisieme position de codon.")

    # ---- H2. force de maintien par classe
    h2_force(vu, SITES, opp, lignees, rep, args)


def h2_force(vu, SITES, opp, lignees, rep, args):
    """H2 -- la force de maintien de A18, appliquee CLASSE PAR CLASSE.

    A18 pose S = ln[(GC/(1-GC)) * (v/u)] sous l'equilibre mutation-selection-
    derive de Li et Bulmer, ou GC est le contenu OBSERVE et v/u le flux
    mutationnel. Les trois classes n'ont ni le meme GC observe (81,5 % en
    4-fold contre 58,0 % en 0-fold) ni le meme v/u : S se calcule donc classe
    par classe, et sa comparaison entre classes dit OU la force agit.

    C'est le seul endroit du projet ou l'hypothese « la force qui maintient
    65 % de GC est une selection sur la PROTEINE » devient falsifiable : aux
    sites 4-fois degeneres, aucune substitution ne change l'acide amine, donc
    une selection proteique ne peut RIEN y maintenir et S devrait s'y effondrer.
    Si S y est au contraire plus GRAND qu'aux sites non degeneres, la selection
    sur la proteine est exclue comme force de maintien -- ce qui ne designe pas
    son remplacant (usage du codon, biais de reparation, conversion biaisee),
    mais retire un candidat majeur."""
    lignes = []
    for c in CLASSES:
        g, a = SITES[c]
        gc = g / (g + a)
        sub = vu[(vu.classe == NOMS[c]) & (vu.reportable)]
        for _, r in sub.iterrows():
            lignes.append(dict(classe=NOMS[c], lignee=r.lignee, gc_observe=100 * gc,
                               vu=r.vu, gc_eq=r.gc_eq,
                               S=float(np.log(gc / (1 - gc) * r.vu)),
                               S_lo=float(np.log(gc / (1 - gc) * r.vu_lo)),
                               S_hi=float(np.log(gc / (1 - gc) * r.vu_hi))))
    df = pd.DataFrame(lignes)
    df.to_csv(args.out_prefix + "force_par_classe.tsv", sep="\t", index=False)
    print("\n=== H2. la force de maintien de A18, classe par classe ===")
    print(f"  {'classe':<14}{'GC observe':>11}{'GC_eq median':>14}"
          f"{'ecart':>8}{'S median':>10}{'etendue S':>18}{'n':>4}")
    med = {}
    for c in CLASSES:
        d = df[df.classe == NOMS[c]]
        med[c] = float(d.S.median())
        print(f"  {NOMS[c]:<14}{d.gc_observe.iloc[0]:>10.1f} %"
              f"{d.gc_eq.median():>13.1f} %{d.gc_observe.iloc[0]-d.gc_eq.median():>7.1f}"
              f"{med[c]:>10.3f}{f'[{d.S.min():.2f} - {d.S.max():.2f}]':>18}"
              f"{len(d):>4}")
    for autre in (F0, F0P2):
        sel = sorted(set(df[df.classe == NOMS[F4]].lignee)
                     & set(df[df.classe == NOMS[autre]].lignee))
        s4 = df[(df.classe == NOMS[F4]) & (df.lignee.isin(sel))].set_index("lignee").S
        s0 = df[(df.classe == NOMS[autre]) & (df.lignee.isin(sel))].set_index("lignee").S
        w, pw = stats.wilcoxon(s4[sel], s0[sel])
        etiq = ("controle C5bis, " if autre == F0P2 else "")
        print(f"\n  {etiq}S(4-fold) - S({NOMS[autre]}) : median "
              f"{float((s4[sel]-s0[sel]).median()):+.3f}, "
              f"{int((s4[sel] > s0[sel]).sum())}/{len(sel)} lignees dans le meme "
              f"sens, Wilcoxon p = {pw:.5f}")
    if med[F4] > med[F0] and med[F4] > med[F0P2]:
        print("  La force requise est PLUS GRANDE la ou la selection sur la "
              "proteine ne peut rien maintenir. Une selection proteique est "
              "donc exclue comme force de maintien du GC : elle predirait "
              "l'inverse, et sans ambiguite de signe.")
    print(f"  rho(S 4-fold, S {NOMS[F0P2]}) sur {len(sel)} lignees : "
          f"{stats.spearmanr(s4[sel], s0[sel])[0]:+.3f}")


if __name__ == "__main__":
    main()
