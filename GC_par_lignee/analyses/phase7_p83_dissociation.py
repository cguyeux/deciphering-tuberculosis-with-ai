#!/usr/bin/env python3
"""
Objet       : P8.3 -- trancher si la DISSOCIATION mise au jour par P8.1 est un
              fait ou un artefact de puissance. P8.1 a etabli deux choses qui ne
              se recouvrent pas. (i) L'heterogeneite contextuelle inter-lignees
              est portee par les canaux de GAIN de paires G:C : T>C et T>G
              portent 274,9 points d'exces de G2 sur 7 878 evenements, contre
              136,5 sur 20 628 pour les pertes, soit cinq fois plus par
              evenement (A29). (ii) Or les deux seuls canaux dont la correlation
              a la force de maintien S survit au controle phylogenetique sont
              des canaux de PERTE, C>T@C.T et C>T@A.G, tandis que le plus fort
              canal de gain, T>C@G.C, s'est revele n'etre que le contraste
              humain/animal deguise en treize observations (lacune 12).
              La composante qui DISTINGUE les lignees entre elles ne serait donc
              pas celle qui SUIT leur force de maintien. Avant de chercher un
              sens a cela, il faut ecarter l'explication ennuyeuse : les gains
              sont trois fois moins nombreux que les pertes, donc leurs
              correlations sont mecaniquement plus bruitees, et leur silence
              pourrait n'etre qu'un manque de resolution.

              LA FAUTE A NE PAS COMMETTRE, et qui est le defaut par defaut de
              ce genre de comparaison : conclure a une difference parce qu'un
              cote est significatif et l'autre non. Une difference de
              significativite n'est pas la significativite d'une difference. Le
              script ne compare donc jamais deux p-valeurs ; il EGALISE la
              puissance et regarde ce qui reste.

              CRITERES PRE-ENREGISTRES, FIXES AVANT TOUTE EXECUTION
                E1 ATTENUATION (mesure directe du bruit, non qualifiant seul).
                   Par canal reportable, erreur type de la part dans la classe
                   estimee par bootstrap SUR LES SOUCHES (l'unite groupee, cf.
                   A23), puis fiabilite R = (var inter-lignees - erreur de
                   mesure moyenne) / var inter-lignees. Le rho maximal
                   detectable vaut approximativement sqrt(R). Si les canaux de
                   gain ont un sqrt(R) tres inferieur a celui des canaux de
                   perte, leur silence n'est pas informatif et la dissociation
                   est suspecte d'etre un artefact de puissance.
                E2 EGALISATION DE PUISSANCE (l'epreuve decisive). Les classes
                   s'apparient naturellement par type de substitution :
                   C>T (transition, perte) avec T>C (transition, gain), et
                   C>A (transversion, perte) avec T>G (transversion, gain). Pour
                   chaque lignee, la classe de perte est AMINCIE par tirage
                   hypergeometrique multivarie a l'effectif exact de la classe
                   de gain appariee, ce qui preserve son profil contextuel en
                   esperance et lui donne le meme bruit d'echantillonnage que la
                   classe de gain. 500 repliques.
                TEST UNIQUE, pre-enregistre. A puissance egalisee, et avec le
                   controle phylogenetique applique d'emblee (les cinq lignees
                   animales forment UN clade, cf. lacune 12) :
                     la dissociation est REELLE si au moins un canal de PERTE
                     aminci garde q_BH <= 0,05 sur les treize lignees dans au
                     moins 50 % des repliques ET un rho de meme signe chez les
                     huit lignees humaines seules dans au moins 50 % des
                     repliques ;
                     elle est un ARTEFACT DE PUISSANCE sinon.
                E3 ACCORD ENTRE EFFECTIFS. Le verdict doit etre le meme a n = 40
                   et a n = 80. S'ils different, le verdict est INDETERMINE et
                   ecrit tel quel. Les deux effectifs sont analyses et rendus,
                   quel que soit leur resultat ; il n'est pas permis d'en
                   choisir un apres coup.
              Aucun autre proxy, aucun autre decoupage, aucune autre classe
              d'appariement ne seront essayes. Si le verdict tombe du cote de
              l'artefact, P8.3 se referme sur cet enonce et ne cherche pas de
              repli, exactement comme P10.3 l'a impose a P10.

Entrees     : résultats/phase7_p81_counts_par_souche.tsv (n = 40, de P8.1)
              bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt (pour n = 80)
              résultats/phase5_p53_force_requise.tsv (S par lignee, A18)
Sorties     : résultats/phase7_p83_attenuation.tsv    (E1)
              résultats/phase7_p83_amincissement.tsv  (E2, canal x effectif)
              résultats/phase7_p83_verdict.tsv        (test unique + E3)
Reutilisable: oui -- l'amincissement hypergeometrique pour egaliser la puissance
              entre deux categories d'effectifs inegaux, et la fiabilite par
              bootstrap sur l'unite groupee, valent pour toute comparaison ou
              l'une des deux categories est structurellement plus rare
Projet      : GC_par_lignee
Date        : 2026-08-30
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase2_polarisation_mtbc0 import build_ancestral  # noqa: E402
from phase5_p51_panel_recursif import PANEL, masked_positions, strains_of  # noqa: E402
from phase7_p81_spectre_contextuel import CANAUX, CTX, strain_channels  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
# appariement par TYPE de substitution : la perte et le gain de meme nature
PAIRES = [("C>T", "T>C"), ("C>A", "T>G")]
PERTES, GAINS = ["C>T", "C>A"], ["T>C", "T>G"]
HUMAINES = ["L1", "L2", "L3", "L4", "L5", "L6", "L7", "L9"]
ANIMALES = ["La4", "Bovis", "Caprae_La2", "Orygis_La3", "Microti"]
TREIZE = HUMAINES + ANIMALES


def charger(n, seed=0):
    """Comptes par souche a l'effectif n. n = 40 est relu de P8.1 ; les autres
    effectifs sont recalcules avec la meme fonction et le meme tirage."""
    cache = ROOT / "résultats" / (f"phase7_p81_counts_par_souche.tsv" if n == 40
                                  else f"phase7_p83_counts_par_souche_n{n}.tsv")
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    anc, masked = build_ancestral(verbose=False), masked_positions()
    rows = []
    for lin, prefixes in PANEL.items():
        st = strains_of(prefixes)
        if len(st) < 4:
            continue
        df, k = strain_channels(st, anc, masked, n, seed)
        if df is None:
            continue
        df.insert(0, "lignee", lin)
        rows.append(df)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(cache, sep="\t", index=False)
    return out


def reportables(agg):
    """C1 de P8.1, reconduit tel quel : >= 100 evenements au total et >= 5 dans
    au moins treize lignees."""
    return (agg.sum(0) >= 100) & ((agg >= 5).sum(0) >= 13)


def parts_classe(agg, classe):
    """Part de chaque canal dans sa classe, par lignee."""
    cols = [f"{classe}@{c}" for c in CTX]
    sub = agg[cols]
    return sub.div(sub.sum(1).replace(0, np.nan), axis=0)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty(len(p))
    q[o] = np.minimum.accumulate(
        (p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def rhos(x, S):
    """Spearman sur les treize, sur les huit humaines, sur les cinq animales,
    et ecart entre groupes. Le controle phylogenetique est applique d'emblee,
    jamais apres coup."""
    r13, p13 = stats.spearmanr(x[TREIZE], S[TREIZE])
    rh, ph = stats.spearmanr(x[HUMAINES], S[HUMAINES])
    ra, pa = stats.spearmanr(x[ANIMALES], S[ANIMALES])
    return r13, p13, rh, ph, ra, pa, float(x[ANIMALES].mean() - x[HUMAINES].mean())


# ------------------------------------------------------------------ E1 bruit
def attenuation(ps, agg, rep, S, boot=400, seed=1):
    """Fiabilite d'un canal : part de sa variance inter-lignees qui n'est pas
    de l'erreur de mesure. Le bootstrap porte sur les SOUCHES, unite groupee
    (A23), et non sur les substitutions."""
    rng = np.random.default_rng(seed)
    se = {}
    for lin in TREIZE:
        sub = ps[ps["lignee"] == lin][CANAUX].to_numpy(float)
        idx = rng.integers(0, len(sub), size=(boot, len(sub)))
        tir = sub[idx].sum(1)                                # boot x 96
        for ci, ca in enumerate(CANAUX):
            k = ca.split("@")[0]
            j = [CANAUX.index(f"{k}@{c}") for c in CTX]
            den = tir[:, j].sum(1)
            with np.errstate(invalid="ignore", divide="ignore"):
                se.setdefault(ca, {})[lin] = float(np.nanstd(
                    np.where(den > 0, tir[:, ci] / den, np.nan)))
    rows = []
    for ca in CANAUX:
        if not rep[ca]:
            continue
        k = ca.split("@")[0]
        x = parts_classe(agg, k)[ca]
        v = float(np.var(x[TREIZE], ddof=1))
        err = float(np.mean([se[ca][l] ** 2 for l in TREIZE]))
        R = max(0.0, (v - err) / v) if v > 0 else 0.0
        r13 = stats.spearmanr(x[TREIZE], S[TREIZE])[0]
        rows.append(dict(canal=ca, classe=k,
                         categorie="perte" if k in PERTES else
                                   ("gain" if k in GAINS else "neutre"),
                         n=int(agg[ca].sum()), var_inter=v, err_mesure=err,
                         fiabilite=R, rho_max_detectable=np.sqrt(R), rho_13=r13))
    return pd.DataFrame(rows)


# ---------------------------------------------------- E2 egalisation par amincissement
def amincir(agg, perte, gain, rng):
    """Tirage hypergeometrique multivarie : dans chaque lignee, la classe de
    perte est ramenee a l'effectif exact de la classe de gain appariee. Le
    profil contextuel est preserve en esperance, le bruit devient celui de la
    classe de gain."""
    cp = [f"{perte}@{c}" for c in CTX]
    cg = [f"{gain}@{c}" for c in CTX]
    out = {}
    for lin in TREIZE:
        obs = agg.loc[lin, cp].to_numpy(int)
        m = int(agg.loc[lin, cg].sum())
        tot = int(obs.sum())
        out[lin] = (rng.multivariate_hypergeometric(obs, min(m, tot))
                    if tot > 0 and m > 0 else np.zeros(16, int))
    return pd.DataFrame(out, index=cp).T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--effectifs", type=int, nargs="+", default=[40, 80])
    ap.add_argument("--rep", type=int, default=500)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    S = pd.read_csv(ROOT / "résultats" / "phase5_p53_force_requise.tsv",
                    sep="\t").set_index("lignee")["S"]
    verdicts, tous_att, tous_am = {}, [], []

    for n in args.effectifs:
        ps = charger(n, args.seed)
        agg = ps.groupby("lignee")[CANAUX].sum()
        rep = reportables(agg)
        print(f"\n{'='*72}\n=== EFFECTIF n = {n} : "
              f"{int(ps['n_evt'].sum()):,} evenements, "
              f"{int(rep.sum())} canaux reportables ===")

        # ---- E1
        att = attenuation(ps, agg, rep, S, seed=args.seed + n)
        att.insert(0, "n_pool", n)
        tous_att.append(att)
        print("\n--- E1. bruit de mesure : le silence des gains est-il informatif ? ---")
        r = att.groupby("categorie").agg(
            n_canaux=("canal", "size"), evt_median=("n", "median"),
            fiabilite_med=("fiabilite", "median"),
            rho_max_med=("rho_max_detectable", "median"),
            rho13_abs_max=("rho_13", lambda v: np.abs(v).max()))
        print(r.round(3).to_string())
        g = att[att.categorie == "gain"]["rho_max_detectable"].median()
        pmed = att[att.categorie == "perte"]["rho_max_detectable"].median()
        print(f"  rho maximal detectable : {g:.3f} (gains) contre {pmed:.3f} "
              f"(pertes) -- rapport {g/max(pmed,1e-9):.2f}")

        # ---- E2
        rng = np.random.default_rng(args.seed + n)
        canaux_p = [c for c in CANAUX if rep[c] and c.split("@")[0] in PERTES]
        canaux_g = [c for c in CANAUX if rep[c] and c.split("@")[0] in GAINS]
        acc = {c: dict(q05=0, sh=0, r13=[], rh=[]) for c in canaux_p}
        for _ in range(args.rep):
            thin = pd.concat([amincir(agg, p, g, rng) for p, g in PAIRES], axis=1)
            stats_rep = {}
            for ca in canaux_p:
                k = ca.split("@")[0]
                cols = [f"{k}@{c}" for c in CTX]
                x = thin[ca] / thin[cols].sum(1).replace(0, np.nan)
                if x.isna().any():
                    continue
                stats_rep[ca] = rhos(x, S)
            if not stats_rep:
                continue
            qs = dict(zip(stats_rep, bh([v[1] for v in stats_rep.values()])))
            for ca, v in stats_rep.items():
                acc[ca]["q05"] += qs[ca] <= 0.05
                acc[ca]["sh"] += np.sign(v[0]) == np.sign(v[2]) and v[2] != 0
                acc[ca]["r13"].append(v[0])
                acc[ca]["rh"].append(v[2])
        rows = []
        for ca, a in acc.items():
            if not a["r13"]:
                continue
            rows.append(dict(n_pool=n, canal=ca, categorie="perte_amincie",
                             n_avant=int(agg[ca].sum()),
                             rho13_median=float(np.median(a["r13"])),
                             rho_humaines_median=float(np.median(a["rh"])),
                             frac_qBH05=a["q05"] / args.rep,
                             frac_signe_concordant=a["sh"] / args.rep))
        am = pd.DataFrame(rows).sort_values("frac_qBH05", ascending=False)
        # les gains, tels quels, pour comparaison
        rows_g = []
        for ca in canaux_g:
            k = ca.split("@")[0]
            x = parts_classe(agg, k)[ca]
            v = rhos(x, S)
            rows_g.append(dict(n_pool=n, canal=ca, categorie="gain_observe",
                               n_avant=int(agg[ca].sum()), rho13_median=v[0],
                               rho_humaines_median=v[2],
                               frac_qBH05=np.nan, frac_signe_concordant=np.nan,
                               p13=v[1], p_humaines=v[3], ecart_groupe=v[6]))
        am = pd.concat([am, pd.DataFrame(rows_g)], ignore_index=True)
        tous_am.append(am)
        print(f"\n--- E2. pertes AMINCIES a l'effectif des gains "
              f"({args.rep} repliques) ---")
        print(am[am.categorie == "perte_amincie"].head(6)[
            ["canal", "n_avant", "rho13_median", "rho_humaines_median",
             "frac_qBH05", "frac_signe_concordant"]].round(3).to_string(index=False))

        # ---- TEST UNIQUE
        ok = am[(am.categorie == "perte_amincie") &
                (am.frac_qBH05 >= 0.5) & (am.frac_signe_concordant >= 0.5)]
        verdicts[n] = "REELLE" if len(ok) else "ARTEFACT DE PUISSANCE"
        print(f"\n  TEST UNIQUE : {len(ok)} canal(aux) de perte gardent "
              f"q_BH <= 0,05 dans >= 50 % des repliques ET un signe concordant "
              f"chez les humaines dans >= 50 %")
        if len(ok):
            print("   " + ", ".join(f"{r.canal} (q05 {r.frac_qBH05:.2f}, "
                                    f"signe {r.frac_signe_concordant:.2f})"
                                    for r in ok.itertuples()))
        print(f"  -> dissociation {verdicts[n]} a n = {n}")

    pd.concat(tous_att).to_csv(ROOT / "résultats" / "phase7_p83_attenuation.tsv",
                               sep="\t", index=False)
    pd.concat(tous_am).to_csv(ROOT / "résultats" / "phase7_p83_amincissement.tsv",
                              sep="\t", index=False)
    fin = ("INDETERMINE (E3 : les effectifs ne s'accordent pas)"
           if len(set(verdicts.values())) > 1 else list(verdicts.values())[0])
    pd.DataFrame([dict(verdict_final=fin,
                       **{f"verdict_n{n}": v for n, v in verdicts.items()},
                       repliques=args.rep)]).to_csv(
        ROOT / "résultats" / "phase7_p83_verdict.tsv", sep="\t", index=False)
    print(f"\n{'='*72}\nE3, accord entre effectifs : "
          + " ; ".join(f"n = {n} -> {v}" for n, v in verdicts.items()))
    print(f"VERDICT P8.3 : la dissociation gains / pertes est {fin}")


if __name__ == "__main__":
    main()
