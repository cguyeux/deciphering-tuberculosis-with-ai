#!/usr/bin/env python3
"""
Objet       : P10.2 -- re-tester l'hypothese de l'efficacite de la selection
              (A20) en changeant non plus l'ESTIMATEUR (ce qu'a fait P10.1) mais
              le PLAN D'ECHANTILLONNAGE. P5 comparait les lignees a EFFECTIF
              egal ; or P10.1 a montre par experience (Q1 contre Q2) que
              l'effectif n'est PAS le confondant -- la derive de rarefaction vaut
              0,029 -- alors que la PROFONDEUR du pool, mesuree par la longueur
              moyenne de branche terminale, deplace le pN/pS terminal de 0,101,
              soit plus que tout l'ecart interquartile inter-lignees. Ce script
              inverse donc le plan : il apparie les lignees sur la profondeur
              L* (variable confondante) et laisse flotter l'effectif (variable
              inoffensive), ce qui revient a lire toutes les lignees a un meme
              AGE MUTATIONNEL des variants prives. A age egal, la selection a eu
              le meme temps pour filtrer, et un ecart residuel de pN/pS est
              alors attribuable a l'efficacite de la selection, c'est-a-dire a
              Ne.s -- exactement la quantite que A20 pretendait mesurer.

              CRITERES, FIXES AVANT DE REGARDER S (meme discipline qu'en P10.1 :
              interdire le mecanisme qui a fabrique A20, essayer des mesures
              jusqu'a ce que l'une rende quelque chose)
                D1 CHOIX DE L*. La profondeur cible est choisie par un critere
                   MECANIQUE portant sur les seules profondeurs, jamais sur le
                   pN/pS ni sur S : L* maximise le nombre de lignees disposant
                   d'au moins R_MIN pools dans la tolerance, les ex aequo etant
                   departages par le nombre minimal de pools qualifiants. Aucun
                   arbitrage humain.
                D2 EXCLUSIONS. Une lignee de moins de N_MIN souches ne peut pas
                   produire de pool apparie ; son exclusion est une propriete du
                   plan, connue avant tout resultat, et non un tri sur donnees.
                Q4 RESOLUTION DU PLAN (remplace le Q2 de P10.1). A profondeur
                   appariee, la LARGEUR du pool (distance moyenne entre souches,
                   qui mesure la part de la lignee balayee) ne doit plus deplacer
                   l'estimateur : ecart median entre pools larges et pools
                   etroits < 0,69 x IQR inter-lignees -- l'artefact residuel de
                   P10.1, que P10.2 doit strictement battre -- et sans biais de
                   signe (test des signes, p > 0,05).
                Q5 BRUIT D'ESTIMATION. L'erreur type de l'estimation d'une
                   lignee -- SE de la mediane sur ses replicats, obtenue par
                   bootstrap, l'ecart-type d'un tirage unique etant rapporte a
                   part -- doit rester sous 0,5 x IQR inter-lignees, sans quoi
                   l'attenuation par erreur de mesure empecherait de distinguer
                   « pas d'effet » de « pas de resolution ».
                Q6 (DIAGNOSTIC, NON QUALIFIANT) dependance residuelle a
                   l'effectif moyen des pools apparies.
              TEST UNIQUE : Spearman(S, pN/pS terminal appariee a L*) sur
              TOUTES les lignees couvertes par le plan (le sous-jeu des lignees
              bien pourvues n'est rapporte que si les deux jeux different),
              rendu quel que soit son resultat, avec IC de Fisher et IC bootstrap
              incluant le bruit de tirage, et verdict d'exclusion du
              rho = -0,736 de A20. SECONDAIRE DECLARE : la diversite pi
              appariee, second proxy de A20, rapportee elle aussi quel que soit
              son resultat mais ne pouvant PAS a elle seule rouvrir A20 (ce
              serait le mecanisme meme que P10.3 interdit). PUISSANCE : sous
              l'effet publie en A20 et le bruit reellement mesure, probabilite
              que ce plan l'aurait vu.
Entrees     : résultats/phase5_p51_counts_par_souche.tsv (liste des lignees)
              résultats/phase5_p53_force_requise.tsv     (S par lignee, A18)
              bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt
              data/MTBC0/ancestral_on_H37Rv.bin, data/cds_opportunity.json
              investigate_phylo/resources/NC_000962.3.gff3 (CDS H37Rv)
Sorties     : résultats/phase6_p102_faisabilite.tsv   (D1, plage par lignee)
              résultats/phase6_p102_pools_apparies.tsv (un pool apparie/ligne)
              résultats/phase6_p102_estimations.tsv   (une lignee/ligne)
              résultats/phase6_p102_qualification.tsv (Q4, Q5, test, puissance)
Reutilisable: oui -- l'appariement d'echantillons de commodite sur une variable
              de plan (ici la profondeur genealogique) et la mesure de sa
              resolution residuelle valent pour toute bacterie clonale
Projet      : GC_par_lignee
Date        : 2026-08-30
"""
import argparse
import json
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
from phase6_p101_pnps_internes import (FAIBLES, PAIRES, RHO_A20,  # noqa: E402
                                       Pool, fisher_ci, mesure)

ROOT = Path(__file__).resolve().parent.parent
OPP_CACHE = ROOT / "data" / "cds_opportunity.json"
N_MIN = 15          # effectif minimal d'un pool apparie (D2)
N_MAX = 50          # effectif maximal ; la plage remplace l'effectif fixe de P5
R_MIN = 10          # pools apparies exiges pour qu'une lignee entre dans le test
R_MAX = 40          # pools apparies retenus par lignee
TOL = (0.05, 0.10, 0.15)   # tolerances relatives essayees dans cet ordre (D1)
CAND = 400          # souches candidates chargees par lignee
ARTEFACT_P101 = 0.69       # artefact residuel de P10.1, en unites d'IQR


def candidats(prefixes, rng, cap=CAND):
    """Jeu de souches candidates d'une lignee, borne a `cap`. Tour de role entre
    repertoires de sous-clade (couverture large de la lignee) puis un bloc dense
    pris dans le repertoire le plus peuple et le plus profond (extremite basse de
    la plage de profondeur). Les lignees a repertoire unique n'ont que le tirage
    aleatoire ; leurs sous-clades sont retrouves plus loin par la distance
    genetique, ce qui rend le plan independant de la taxonomie des repertoires."""
    par_dir = {}
    for d in pool_dirs(prefixes):
        s = [x for x in sorted(d.iterdir())
             if x.is_dir() and (x / "NC_000962.3" / "spdi.txt").exists()]
        if s:
            par_dir[d.name] = s
    total = sum(len(v) for v in par_dir.values())
    if total <= cap:
        return [x for v in par_dir.values() for x in v]
    for v in par_dir.values():
        rng.shuffle(v)
    out, tour = [], 0
    noms = sorted(par_dir, key=lambda k: -len(par_dir[k]))
    while len(out) < int(cap * 0.75):
        pris = False
        for k in noms:
            if tour < len(par_dir[k]):
                out.append(par_dir[k][tour])
                pris = True
                if len(out) >= int(cap * 0.75):
                    break
        if not pris:
            break
        tour += 1
    dense = max(par_dir, key=lambda k: (len(par_dir[k]) >= 40, k.count(".")))
    vus = set(map(str, out))
    for x in par_dir[dense]:
        if len(out) >= cap:
            break
        if str(x) not in vus:
            out.append(x)
    return out


def charge(paths, masked, anc, seq, geno, cache):
    subsets, noms = [], []
    for p in paths:
        v = read_subs(p / "NC_000962.3" / "spdi.txt")
        if v:
            subsets.append(v)
            noms.append(p.name)
    if len(subsets) < N_MIN:
        return None, noms
    return Pool(subsets, masked, anc, seq, geno, cache), noms


def distances(pool):
    """Distances de Hamming entre souches candidates, sur tous les variants non
    masques. Sert a definir des VOISINAGES genetiques, donc des sous-clades, sans
    passer par les noms de repertoires."""
    X = pool.M.astype(np.float32)
    k = X.sum(0)
    D = k[:, None] + k[None, :] - 2.0 * (X.T @ X)
    np.fill_diagonal(D, 0.0)
    return np.maximum(D, 0.0)


def bibliotheque(pool, D, order, ndraw, rng):
    """Bibliotheque de pools candidats couvrant le plan (profondeur x largeur).

    Un tirage : une souche germe, une PORTEE m tiree log-uniformement entre n et
    l'effectif candidat (m = n donne les n plus proches voisins, donc un pool
    etroit et peu profond ; m maximal donne un pool balayant toute la lignee),
    puis n souches tirees dans cette portee. La profondeur est la longueur
    moyenne de branche terminale (variants polarises codants portes par une
    seule souche du pool, par souche), exactement la quantite que Q2 de P10.1 a
    montree responsable de l'artefact. La largeur est la distance moyenne entre
    souches du pool."""
    nc = pool.n
    nmax = min(N_MAX, nc)
    if nmax < N_MIN:
        return []
    lib = []
    for _ in range(ndraw):
        n = int(rng.integers(N_MIN, nmax + 1))
        lo, hi = np.log(n), np.log(nc)
        m = int(round(np.exp(rng.uniform(lo, hi)))) if hi > lo else n
        m = min(max(m, n), nc)
        seed = int(rng.integers(nc))
        scope = order[seed, :m]
        idx = rng.choice(scope, n, replace=False)
        k = pool.A[:, idx].sum(1)
        sub = D[np.ix_(idx, idx)]
        lib.append(dict(idx=idx, n=n, m=m,
                        profondeur=float((k == 1).sum()) / n,
                        largeur=float(sub.sum()) / (n * (n - 1))))
    return lib


def choisir_Lstar(libs, grille, tol):
    """D1 : L* choisi sur les seules PROFONDEURS. Maximise le nombre de lignees
    ayant au moins R_MIN pools dans la tolerance ; ex aequo departages par le
    nombre minimal de pools qualifiants parmi les lignees couvertes."""
    best = None
    for L in grille:
        cnt = {lin: sum(1 for r in lib
                        if abs(np.log(r["profondeur"] / L)) <= np.log(1 + tol))
               for lin, lib in libs.items()}
        couvertes = [lin for lin, c in cnt.items() if c >= R_MIN]
        if not couvertes:
            continue
        cle = (len(couvertes), min(cnt[l] for l in couvertes))
        if best is None or cle > best[0]:
            best = (cle, L, couvertes, cnt)
    return best


def calibre_pente(z, cible, rng, nsim=400):
    """Pente a telle que Spearman(S, -a.z + bruit) vaille en moyenne `cible`."""
    for a in np.linspace(0.05, 8.0, 160):
        r = [stats.spearmanr(z, -a * z + rng.normal(size=len(z)))[0]
             for _ in range(nsim)]
        if np.mean(r) <= cible:
            return a
    return 8.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default=str(ROOT / "résultats" /
                                            "phase5_p51_counts_par_souche.tsv"))
    ap.add_argument("--force", default=str(ROOT / "résultats" /
                                           "phase5_p53_force_requise.tsv"))
    ap.add_argument("--draws", type=int, default=4000,
                    help="pools candidats tires par lignee")
    ap.add_argument("--cand", type=int, default=CAND)
    ap.add_argument("-B", "--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-prefix", default=str(ROOT / "résultats" /
                                                "phase6_p102_"))
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
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

    cs = pd.read_csv(args.counts, sep="\t")
    lignees = [l for l in PANEL if l in set(cs.lignee)]
    force = pd.read_csv(args.force, sep="\t")[["lignee", "S", "vu", "gc_eq"]]
    cache = {}

    # ---------- 1. BIBLIOTHEQUE DE POOLS PAR LIGNEE ---------------------------
    print("=== 1. PLAN : bibliotheque de pools (profondeur x largeur) ===")
    print("  Le plan de P5 fixait l'EFFECTIF ; Q1 de P10.1 a montre que "
          "l'effectif n'est pas\n  le confondant. Ici on fixe la PROFONDEUR, "
          "que Q2 a montree responsable de\n  l'artefact, et on laisse flotter "
          f"l'effectif entre {N_MIN} et {N_MAX}.")
    pools, libs, exclues = {}, {}, {}
    for lin in lignees:
        paths = candidats(PANEL[lin], np.random.default_rng(args.seed), args.cand)
        p, noms = charge(paths, masked, anc, seq, geno, cache)
        if p is None:
            exclues[lin] = f"{len(noms)} souches < N_MIN = {N_MIN} (D2)"
            print(f"  {lin:12s} EXCLUE : {exclues[lin]}")
            continue
        D = distances(p)
        order = np.argsort(D, axis=1, kind="stable")
        lib = bibliotheque(p, D, order, args.draws, rng)
        pools[lin], libs[lin] = p, lib
        prof = np.array([r["profondeur"] for r in lib])
        print(f"  {lin:12s} {p.n:4d} candidates, {p.M.shape[0]:6d} variants, "
              f"{p.A.shape[0]:6d} polarises codants ; profondeur atteignable "
              f"[{np.percentile(prof, 2):7.2f} ; {np.percentile(prof, 98):8.2f}]")
    faisab = pd.DataFrame([
        dict(lignee=lin, n_candidates=pools[lin].n, n_pools=len(libs[lin]),
             prof_min=float(np.min([r["profondeur"] for r in libs[lin]])),
             prof_p02=float(np.percentile([r["profondeur"] for r in libs[lin]], 2)),
             prof_med=float(np.median([r["profondeur"] for r in libs[lin]])),
             prof_p98=float(np.percentile([r["profondeur"] for r in libs[lin]], 98)),
             prof_max=float(np.max([r["profondeur"] for r in libs[lin]])))
        for lin in libs] + [dict(lignee=l, n_candidates=np.nan, n_pools=0)
                            for l in exclues])
    faisab.to_csv(args.out_prefix + "faisabilite.tsv", sep="\t", index=False)

    # ---------- 2. D1 : CHOIX MECANIQUE DE L* ---------------------------------
    print("\n=== 2. D1 : profondeur cible L*, choisie sur les seules "
          "profondeurs ===")
    tous = np.concatenate([[r["profondeur"] for r in lib] for lib in libs.values()])
    grille = np.exp(np.linspace(np.log(max(tous.min(), 0.5)),
                                np.log(tous.max()), 400))
    for tol in TOL:
        best = choisir_Lstar(libs, grille, tol)
        if best and best[0][0] >= 8:
            break
    (nlin, nmin_pools), Lstar, couvertes, cnt = best
    print(f"  L* = {Lstar:.2f} variants prives par souche, tolerance "
          f"+/-{tol:.0%} -> {nlin} lignees couvertes "
          f"(au moins {nmin_pools} pools chacune)")
    for lin in libs:
        etat = "couverte" if lin in couvertes else "hors plage"
        print(f"    {lin:12s} {cnt[lin]:5d} pools dans la tolerance  {etat}")
    non_couv = [l for l in libs if l not in couvertes]

    # ---------- 3. POOLS APPARIES ET ESTIMATIONS ------------------------------
    print("\n=== 3. Mesure sur les pools apparies ===")
    rows = []
    for lin in couvertes:
        q = [r for r in libs[lin]
             if abs(np.log(r["profondeur"] / Lstar)) <= np.log(1 + tol)]
        if len(q) > R_MAX:
            q = [q[i] for i in rng.choice(len(q), R_MAX, replace=False)]
        for j, r in enumerate(q):
            m = mesure(pools[lin], opp_syn, opp_non, n_sites, idx=r["idx"])
            rows.append(dict(lignee=lin, replicat=j, n=r["n"], portee=r["m"],
                             profondeur=r["profondeur"], largeur=r["largeur"],
                             pn_ps_term=m["pn_ps_term"], pn_ps_int=m["pn_ps_int"],
                             pi=m["pi"], n_syn_term=m["n_syn_term"],
                             n_nonsyn_term=m["n_nonsyn_term"]))
    tp = pd.DataFrame(rows)
    tp.to_csv(args.out_prefix + "pools_apparies.tsv", sep="\t", index=False)

    est = tp.groupby("lignee").agg(
        n_replicats=("replicat", "size"), n_moyen=("n", "mean"),
        profondeur=("profondeur", "mean"), largeur=("largeur", "mean"),
        pn_ps_term=("pn_ps_term", "median"), sd_term=("pn_ps_term", "std"),
        pn_ps_int=("pn_ps_int", "median"), pi=("pi", "median"),
        sd_pi=("pi", "std"), n_term_moyen=("n_syn_term", "mean")).reset_index()
    est = est.merge(force, on="lignee")
    est["sous_puissante"] = est.lignee.isin(FAIBLES)
    est = est.sort_values("S", ascending=False)
    est.to_csv(args.out_prefix + "estimations.tsv", sep="\t", index=False)
    print(est[["lignee", "n_replicats", "n_moyen", "profondeur", "largeur",
               "pn_ps_term", "sd_term", "pn_ps_int", "pi", "S"]]
          .round(4).to_string(index=False))
    iqr = {c: est[c].quantile(.75) - est[c].quantile(.25)
           for c in ("pn_ps_term", "pn_ps_int", "pi")}
    for c, v in iqr.items():
        print(f"  IQR inter-lignees de {c:11s} = {v:.4g} "
              f"(mediane {est[c].median():.4g})")

    # ---------- 4. Q4 : RESOLUTION DU PLAN ------------------------------------
    print("\n=== 4. Q4 : a profondeur appariee, la LARGEUR deplace-t-elle "
          "encore l'estimateur ? ===")
    print("  Meme lignee, meme profondeur L*, mais des pools qui balaient une "
          "part plus ou\n  moins large de la lignee. C'est le Q2 de P10.1 refait "
          "sur le plan et non sur\n  l'estimateur : si l'appariement suffit, "
          "l'ecart doit tomber sous "
          f"{ARTEFACT_P101} x IQR.")
    q4 = []
    for lin, g in tp.groupby("lignee"):
        med = g.largeur.median()
        large, etroit = g[g.largeur > med], g[g.largeur <= med]
        if len(large) < 3 or len(etroit) < 3:
            continue
        q4.append(dict(lignee=lin, largeur_large=large.largeur.median(),
                       largeur_etroit=etroit.largeur.median(),
                       **{f"{c}_large": large[c].median() for c in
                          ("pn_ps_term", "pn_ps_int", "pi")},
                       **{f"{c}_etroit": etroit[c].median() for c in
                          ("pn_ps_term", "pn_ps_int", "pi")}))
    tq4 = pd.DataFrame(q4)
    print(tq4.round(4).to_string(index=False))
    verdict_q4 = {}
    for c, lab in (("pn_ps_term", "pN/pS terminal"),
                   ("pn_ps_int", "pN/pS interne"), ("pi", "diversite pi")):
        d = (tq4[f"{c}_large"] - tq4[f"{c}_etroit"]).dropna()
        npos = int((d > 0).sum())
        ps = stats.binomtest(npos, len(d), 0.5).pvalue
        ratio = d.abs().median() / iqr[c]
        ok = ratio < ARTEFACT_P101 and ps > 0.05
        verdict_q4[c] = dict(ecart_median=d.abs().median(), ratio_iqr=ratio,
                             signes=f"{npos}/{len(d)}", p_signes=ps, passe=ok)
        print(f"  {lab:15s} : ecart median |large - etroit| = "
              f"{d.abs().median():.4g} = {ratio:.2f} x IQR, signes "
              f"{npos}/{len(d)} p = {ps:.3g} -> "
              f"{'PASSE' if ok else 'ECHOUE'} Q4 "
              f"(P10.1 laissait {ARTEFACT_P101:.2f} x IQR)")

    # ---------- 5. BOOTSTRAP, Q5, Q6 ------------------------------------------
    reps = {lin: g.pn_ps_term.dropna().values for lin, g in tp.groupby("lignee")}
    reps_pi = {lin: g.pi.dropna().values for lin, g in tp.groupby("lignee")}
    lins = list(est.lignee)
    Sv = est.S.values
    boot_med = {lin: np.empty(args.boot) for lin in lins}
    rho_boot, rho_boot_pi = np.empty(args.boot), np.empty(args.boot)
    for b in range(args.boot):
        v, vp = [], []
        for lin in lins:
            r = reps[lin]
            m = np.median(r[rng.integers(0, len(r), len(r))])
            boot_med[lin][b] = m
            v.append(m)
            rp = reps_pi[lin]
            vp.append(np.median(rp[rng.integers(0, len(rp), len(rp))]))
        rho_boot[b] = stats.spearmanr(Sv, v)[0]
        rho_boot_pi[b] = stats.spearmanr(Sv, vp)[0]
    se_med = np.array([boot_med[lin].std() for lin in lins])
    est["se_mediane"] = se_med
    print("\n=== 5. Q5 : bruit d'estimation contre dispersion inter-lignees ===")
    print(f"  SE mediane de l'estimation par lignee = {np.median(se_med):.4g} "
          f"= {np.median(se_med) / iqr['pn_ps_term']:.2f} x IQR "
          f"(SD d'un tirage unique : {est.sd_term.median():.4g} = "
          f"{est.sd_term.median() / iqr['pn_ps_term']:.2f} x IQR)")
    q5 = np.median(se_med) < 0.5 * iqr["pn_ps_term"]
    print(f"  -> {'PASSE' if q5 else 'ECHOUE'} Q5 (seuil 0,50 x IQR)")

    print("\n=== 6. Q6 (diagnostic, non qualifiant) : effectif residuel ===")
    rho_n, p_n = stats.spearmanr(est.n_moyen, est.pn_ps_term)
    print(f"  Spearman(effectif moyen des pools, pN/pS terminal) = {rho_n:+.3f} "
          f"(p = {p_n:.3g})")
    rk = lambda x: stats.rankdata(x)
    rS, rY, rN = rk(Sv), rk(est.pn_ps_term.values), rk(est.n_moyen.values)
    resS = rS - np.polyval(np.polyfit(rN, rS, 1), rN)
    resY = rY - np.polyval(np.polyfit(rN, rY, 1), rN)
    rho_part, p_part = stats.pearsonr(resS, resY)
    print(f"  Spearman partiel(S, pN/pS terminal | effectif) = {rho_part:+.3f} "
          f"(p = {p_part:.3g})")

    # ---------- 7. TEST UNIQUE PRE-ENREGISTRE ----------------------------------
    print("\n=== 7. TEST UNIQUE de l'hypothese de A20 sur le plan apparie ===")
    print("  Attendu si S = 2.Ne.s a s constant : correlation NEGATIVE entre S "
          "et pN/pS.")
    lignes_v = []
    for c, lab, rb in (("pn_ps_term", "pN/pS terminal apparie (PRIMAIRE)",
                        rho_boot),
                       ("pi", "diversite pi appariee (secondaire declare)",
                        rho_boot_pi)):
        for sous, etiq in ((est, f"{len(est)} lignees couvertes"),
                           (est[~est.sous_puissante],
                            f"{int((~est.sous_puissante).sum())} bien pourvues")):
            if len(sous) == len(est) and etiq.endswith("bien pourvues"):
                continue
            u = sous.dropna(subset=[c])
            rho, p = stats.spearmanr(u.S, u[c])
            lo, hi = fisher_ci(rho, len(u))
            excl = not lo <= RHO_A20 <= hi
            print(f"  Spearman(S, {lab}) sur {etiq} (n = {len(u)}) : "
                  f"{rho:+.3f} (p = {p:.3g}), IC 95 % [{lo:+.3f} ; {hi:+.3f}] "
                  f"-> {'exclut' if excl else 'compatible avec'} "
                  f"le rho = {RHO_A20:+.3f} de A20")
            lignes_v.append(dict(quantite=c, jeu=etiq, n=len(u), rho=rho, p=p,
                                 ic_lo=lo, ic_hi=hi, exclut_a20=excl))
        blo, bhi = np.percentile(rb, [2.5, 97.5])
        print(f"    Bootstrap sur les replicats de pools (STABILITE du rho "
              f"vis-a-vis du tirage,\n    non incertitude inter-lignees) : "
              f"[{blo:+.3f} ; {bhi:+.3f}], mediane {np.median(rb):+.3f}")
        lignes_v.append(dict(quantite=c, jeu="bootstrap replicats (stabilite)",
                             n=len(est), rho=float(np.median(rb)), p=np.nan,
                             ic_lo=blo, ic_hi=bhi,
                             exclut_a20=not blo <= RHO_A20 <= bhi))

    # ---------- 8. PUISSANCE : ce plan aurait-il vu l'effet de A20 ? -----------
    print("\n=== 8. PUISSANCE : sous l'effet publie en A20 et le bruit "
          "reellement mesure ===")
    k = len(est)
    z = stats.zscore(stats.rankdata(Sv))
    a = calibre_pente(z, RHO_A20, rng)
    sigma_b = float(est.pn_ps_term.std())
    med0 = float(est.pn_ps_term.median())
    nsim, det, rhos = 5000, 0, np.empty(5000)
    for i in range(nsim):
        lat = -a * z + rng.normal(size=k)
        y = med0 + sigma_b * (lat - lat.mean()) / lat.std()
        obs = y + rng.normal(0, se_med)
        r, pv = stats.spearmanr(Sv, obs)
        rhos[i] = r
        det += (pv < 0.05) and (r < 0)
    puiss = det / nsim
    print(f"  Effet simule : rho vrai = {RHO_A20:+.3f} sur {k} lignees, "
          f"dispersion inter-lignees observee ({sigma_b:.4g}),\n  bruit "
          f"d'estimation observe (SE mediane {np.median(se_med):.4g}). "
          f"Rho observe simule : {np.median(rhos):+.3f}\n  (attenuation "
          f"{abs(np.median(rhos) / RHO_A20):.0%} du signal conserve). "
          f"PUISSANCE = {puiss:.1%}")

    # ---------- 9. VERDICT ------------------------------------------------------
    prim = [r for r in lignes_v if r["quantite"] == "pn_ps_term"
            and r["jeu"].endswith("couvertes")][0]
    q4_ok = verdict_q4["pn_ps_term"]["passe"]
    print("\n=== 9. VERDICT pour P10.3 ===")
    if not (q4_ok and q5):
        print("  Le plan N'A PAS la resolution exigee : son silence ne "
              "trancherait rien.")
    elif prim["rho"] < 0 and prim["p"] < 0.05:
        print("  L'effet de A20 REAPPARAIT une fois les lignees appariees sur "
              "la profondeur :\n  A20 renait sous condition d'appariement.")
    else:
        print("  Deuxieme mesure independante, resolue, qui ne rend rien : "
              "le critere P10.3 est\n  atteint. L'explication par l'efficacite "
              "de la selection est ABANDONNEE, et le\n  manuscrit ecrit que le "
              "classement des GC d'equilibre est etabli mais inexplique.")
    print("  Aucun troisieme proxy n'est essaye : P10.3 l'interdit.")

    pd.DataFrame(lignes_v + [
        dict(quantite="Q4_" + c, jeu=f"ecart large/etroit ({v['signes']})",
             n=len(tq4), rho=v["ratio_iqr"], p=v["p_signes"],
             ic_lo=np.nan, ic_hi=np.nan, exclut_a20=v["passe"])
        for c, v in verdict_q4.items()] + [
        dict(quantite="Q5_se_mediane", jeu="SE mediane / IQR", n=k,
             rho=float(np.median(se_med) / iqr["pn_ps_term"]), p=np.nan,
             ic_lo=np.nan, ic_hi=np.nan, exclut_a20=bool(q5)),
        dict(quantite="Q6_partiel", jeu="S vs pN/pS | effectif", n=k,
             rho=rho_part, p=p_part, ic_lo=np.nan, ic_hi=np.nan,
             exclut_a20=False),
        dict(quantite="puissance", jeu=f"rho vrai = {RHO_A20}", n=k,
             rho=float(np.median(rhos)), p=puiss, ic_lo=np.nan, ic_hi=np.nan,
             exclut_a20=False),
        dict(quantite="Lstar", jeu=f"tolerance {tol}", n=len(couvertes),
             rho=Lstar, p=np.nan, ic_lo=np.nan, ic_hi=np.nan,
             exclut_a20=False)]).to_csv(
        args.out_prefix + "qualification.tsv", sep="\t", index=False)
    if non_couv:
        print("\n  Lignees hors plage de L* (exclusion de plan, pas de "
              "resultat) : " + ", ".join(non_couv))
    if exclues:
        print("  Lignees sans pool possible : " +
              ", ".join(f"{k_}" for k_ in exclues))


if __name__ == "__main__":
    main()
