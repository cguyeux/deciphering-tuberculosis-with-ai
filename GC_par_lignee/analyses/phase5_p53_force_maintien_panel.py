#!/usr/bin/env python3
"""
Objet       : P5.3 -- refaire A18 (force requise S), A19 (cinetique et ΔGC deja
              accumule depuis le MRCA) et A20 (pN/pS corrige du spectre, proxys
              de Ne) sur le PANEL ELARGI de P5.1 : dix-sept lignees mesurables
              au lieu de onze, a la maille lignee, et lignee 1 COMPLETE (les
              deux conteneurs `L1/` et `L_1.*`, A15 la mesurait sur 84,5 %).
              C'est le point 9 des angles morts de l'etat : A18/A19/A20 tenaient
              encore sur onze points et sur une L1 amputee.
              Deux ecarts assumes avec phase4_p93, tous deux imposes par A23 :
              (a) les IC de S et de v/u sont BOOTSTRAP HIERARCHIQUES (souches
              avec remise puis comptes sous Poisson), l'IC binomial sous-estimant
              l'incertitude d'un facteur 1,61 en mediane ; (b) le ΔGC par lignee
              recoit lui aussi un IC bootstrap sur souches, et les correlations
              de la section D sont rendues sur les dix-sept lignees ET sur les
              treize bien pourvues (hors L8, Suricattae, Dassie, Pinnipedii, que
              A23 declare sous-puissantes), parce qu'une correlation de rang est
              sensible a quatre points bruites.
              Le pool de chaque lignee n'est PAS retire : il est lu depuis
              `phase5_p51_counts_par_souche.tsv` (colonne sra), pour que les
              quatre sections portent exactement sur les memes souches que A22
              et A23.
Entrees     : résultats/phase5_p51_counts_par_souche.tsv (pools et comptes P5.1)
              résultats/phase4_p93_[force_requise|selection].tsv (comparaison n=11)
              bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt
              data/MTBC0/ancestral_on_H37Rv.bin, data/MTBC0/Mask.files/*.bed
              investigate_phylo/resources/NC_000962.3.gff3 (CDS H37Rv)
Sorties     : résultats/phase5_p53_force_requise.tsv    (A)
              résultats/phase5_p53_cinetique.tsv        (B)
              résultats/phase5_p53_dgc_depuis_mrca.tsv  (C, une ligne/souche)
              résultats/phase5_p53_dgc_par_lignee.tsv   (C, une ligne/lignee)
              résultats/phase5_p53_selection.tsv        (D)
              résultats/phase5_p53_comparaison_n11.tsv  (E, A18/A20 avant/apres)
              résultats/phase5_p53_[diagnostic_a20|loo_a20].tsv (F)
              data/cds_opportunity.json                 (cache d'opportunite codante)
Reutilisable: oui -- (A) et (B) valent pour toute bacterie dont on a mesure v/u ;
              (D) est un pN/pS corrige du spectre ; le cache d'opportunite codante
              est reutilisable par tout script travaillant sur H37Rv
Projet      : GC_par_lignee
Date        : 2026-08-30
"""
import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import read_subs, flux  # noqa: E402
from phase2_polarisation_mtbc0 import build_ancestral, read_fasta, H37RV  # noqa: E402
from phase3_sueoka_gc_eq import opportunity, vu_ci  # noqa: E402
from phase4_p93_force_maintien import (cds_opportunity, effect,  # noqa: E402
                                       load_cds)
from phase5_p51_panel_recursif import PANEL, masked_positions, strains_of  # noqa: E402
from phase5_p52_bootstrap_puissance import boot_vu  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OPP_CACHE = ROOT / "data" / "cds_opportunity.json"
# Sous-puissantes selon A23 (MDE >= x2,8) : gardees dans les tableaux, ecartees
# des correlations de rang, ou quatre points bruites pesent lourd sur n = 17.
FAIBLES = {"L8", "Suricattae", "Dassie", "Pinnipedii"}
# Correspondance de maille avec le panel n = 11 de A18/A20 : meme lignee lue au
# meme niveau taxonomique d'un cote et de l'autre.
MEME_MAILLE = {"L1": "L1", "L3": "L3", "L7": "L7", "L9": "L9",
               "Orygis_La3": "Orygis_La3", "Caprae_La2": "Caprae_La2",
               "Microti": "Microti"}


def pools_de_p51(counts_path):
    """{lignee: [sra, ...]} exactement tel que P5.1 les a tires."""
    df = pd.read_csv(counts_path, sep="\t")
    return df, {lin: list(d["sra"]) for lin, d in df.groupby("lignee")}


def chemins(lignee, sras):
    """Resout chaque SRA du pool en son repertoire de souche."""
    index = {s.name: s for s in strains_of(PANEL[lignee])}
    out, manquants = [], []
    for sra in sras:
        if sra in index:
            out.append(index[sra])
        else:
            manquants.append(sra)
    if manquants:
        print(f"# {lignee} : {len(manquants)} SRA introuvables ({manquants[:3]})",
              file=sys.stderr)
    return out


def opportunite_codante(seq, masked, owner, off, strand, cds):
    """Opportunite codante syn/nonsyn par paire (ref, alt), avec cache disque :
    le calcul balaie 4,4 Mb x 3 alleles et ne depend que de H37Rv et du masque."""
    if OPP_CACHE.exists():
        raw = json.loads(OPP_CACHE.read_text())
        if raw.get("n_masked") == len(masked):
            print(f"# opportunite codante lue du cache {OPP_CACHE.name}",
                  file=sys.stderr)
            return {tuple(k.split(">")): tuple(v) for k, v in raw["opp"].items()}
    print("# opportunite codante : calcul complet (quelques minutes)...",
          file=sys.stderr)
    opp = cds_opportunity(seq, masked, owner, off, strand, cds)
    OPP_CACHE.write_text(json.dumps(
        {"n_masked": len(masked),
         "opp": {f"{r}>{a}": list(v) for (r, a), v in opp.items()}}))
    return opp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default=str(ROOT / "résultats" /
                                            "phase5_p51_counts_par_souche.tsv"))
    ap.add_argument("--ancien-force", default=str(ROOT / "résultats" /
                                                  "phase4_p93_force_requise.tsv"))
    ap.add_argument("--ancien-selection", default=str(ROOT / "résultats" /
                                                      "phase4_p93_selection.tsv"))
    ap.add_argument("--clock", type=float, default=1e-7,
                    help="horloge MTBC en substitutions/site/an (Menardo 2019)")
    ap.add_argument("--clock-lo", type=float, default=5e-8)
    ap.add_argument("--clock-hi", type=float, default=2e-7)
    ap.add_argument("--ages", type=float, nargs="*", default=[6000, 20000, 1e5, 1e6])
    ap.add_argument("-B", "--boot", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-prefix",
                    default=str(ROOT / "résultats" / "phase5_p53_"))
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    seq = read_fasta(H37RV)
    anc = build_ancestral()
    masked = masked_positions()
    (n_gc, n_at), _ = opportunity()
    n_sites = n_gc + n_at
    gc_obs = n_gc / n_sites
    print(f"# genome non masque : {n_gc} G:C, {n_at} A:T, GC = {100*gc_obs:.4f} %",
          file=sys.stderr)

    cs, pools = pools_de_p51(args.counts)
    lignees = [l for l in PANEL if l in pools]
    print(f"# panel : {len(lignees)} lignees, "
          f"{sum(len(v) for v in pools.values())} souches", file=sys.stderr)

    # ---------- A. FORCE REQUISE (Li 1987 / Bulmer 1991), IC hierarchiques -----
    f_odds = gc_obs / (1 - gc_obs)
    rows = []
    for lin in lignees:
        d = cs[cs["lignee"] == lin]
        loss, gain = d["loss"].to_numpy(), d["gain"].to_numpy()
        L, G = int(loss.sum()), int(gain.sum())
        vu, blo, bhi = vu_ci(L, G, n_gc, n_at)
        bs = boot_vu(loss, gain, n_gc, n_at, args.boot, rng)
        lo, hi = np.percentile(bs, [2.5, 97.5])
        rows.append(dict(lignee=lin, n_souches=len(d), loss=L, gain=G, vu=vu,
                         vu_boot_lo=lo, vu_boot_hi=hi, gc_eq=1 / (1 + vu),
                         gc_eq_boot_lo=1 / (1 + hi), gc_eq_boot_hi=1 / (1 + lo),
                         S=np.log(f_odds * vu),
                         S_binom_lo=np.log(f_odds * blo),
                         S_binom_hi=np.log(f_odds * bhi),
                         S_boot_lo=np.log(f_odds * lo),
                         S_boot_hi=np.log(f_odds * hi),
                         ecart_gc_pts=100 * (gc_obs - 1 / (1 + vu)),
                         sous_puissante=lin in FAIBLES))
    ta = pd.DataFrame(rows).sort_values("S", ascending=False)
    ta.to_csv(args.out_prefix + "force_requise.tsv", sep="\t", index=False)
    print(f"\n=== A. FORCE REQUISE pour tenir GC = {100*gc_obs:.3f} %, "
          f"panel a {len(ta)} lignees ===")
    print("S = ln[(GC/(1-GC)) . v/u] : coefficient d'echelle (2.Ne.s en haploide) "
          "que devrait avoir\n    toute force favorisant G:C pour que 65,4 % soit "
          "un equilibre mutation-selection-derive.\n    IC bootstrap HIERARCHIQUE "
          "(A23), l'IC binomial est rappele pour memoire.")
    print(ta[["lignee", "n_souches", "loss", "gain", "vu", "S", "S_boot_lo",
              "S_boot_hi", "S_binom_lo", "S_binom_hi", "ecart_gc_pts",
              "sous_puissante"]].round(4).to_string(index=False))
    print(f"  etendue : S = {ta.S.min():.3f} ({ta.lignee.iloc[-1]}) a "
          f"{ta.S.max():.3f} ({ta.lignee.iloc[0]}), rapport "
          f"x{ta.S.max()/ta.S.min():.2f}")
    print(f"  en probabilite de fixation relative, e^S va de "
          f"{np.exp(ta.S.min()):.2f} a {np.exp(ta.S.max()):.2f}")
    fort = ta[~ta.sous_puissante]
    print(f"  hors les quatre sous-puissantes de A23 : S = {fort.S.min():.3f} "
          f"({fort.lignee.iloc[-1]}) a {fort.S.max():.3f} ({fort.lignee.iloc[0]}), "
          f"rapport x{fort.S.max()/fort.S.min():.2f}")
    larg = (ta.S_boot_hi - ta.S_boot_lo) / (ta.S_binom_hi - ta.S_binom_lo)
    print(f"  elargissement median de l'IC de S par le bootstrap hierarchique : "
          f"x{larg.median():.2f} (min x{larg.min():.2f}, max x{larg.max():.2f})")

    # ---------- B. CINETIQUE ---------------------------------------------------
    rows = []
    for lin in lignees:
        d = cs[cs["lignee"] == lin]
        L, G, N = int(d["loss"].sum()), int(d["gain"].sum()), int(d["neutral"].sum())
        frac_gc_changing = (L + G) / (L + G + N)
        for tag, clk in (("central", args.clock), ("lo", args.clock_lo),
                         ("hi", args.clock_hi)):
            R = clk * frac_gc_changing                 # subs GC-changeantes/site/an
            v = R * (L / (L + G)) * n_sites / n_gc     # par site G:C et par an
            u = R * (G / (L + G)) * n_sites / n_at     # par site A:T et par an
            lam = u + v
            gceq = u / (u + v)
            d2 = dict(lignee=lin, horloge=tag, clock=clk,
                      frac_gc_changing=frac_gc_changing, u=u, v=v, lambda_=lam,
                      demi_vie_ans=np.log(2) / lam, gc_eq=gceq,
                      subs_par_genome_par_an=clk * n_sites)
            for t in args.ages:
                # signe negatif = perte de GC (la composition descend vers gceq)
                d2[f"dGC_ppm_{int(t)}ans"] = 1e6 * (gceq - gc_obs) * -np.expm1(-lam * t)
            rows.append(d2)
    tb = pd.DataFrame(rows)
    tb.to_csv(args.out_prefix + "cinetique.tsv", sep="\t", index=False)
    show = tb[tb.horloge == "central"].copy()
    show["u_e8"], show["v_e8"] = show.u * 1e8, show.v * 1e8
    show["demi_vie_Man"] = show.demi_vie_ans / 1e6
    print("\n=== B. CINETIQUE : a quelle vitesse la composition repondrait-elle ? ===")
    print(show[["lignee", "u_e8", "v_e8", "demi_vie_Man"] +
               [f"dGC_ppm_{int(t)}ans" for t in args.ages]]
          .round(2).to_string(index=False))
    print("  u_e8, v_e8 : taux par site A:T et par site G:C, en 1e-8/an ; "
          "demi_vie_Man en millions d'annees")
    print(f"  dGC en ppm, signe NEGATIF = perte de GC ; horloge "
          f"{args.clock:.0e} subs/site/an, soit "
          f"{args.clock*n_sites:.2f} substitution/genome/an")
    col0 = f"dGC_ppm_{int(args.ages[0])}ans"
    print(f"  demi-vie : {show.demi_vie_Man.min():.2f} a "
          f"{show.demi_vie_Man.max():.2f} millions d'annees")
    print(f"  etendue INTER-LIGNEES du deplacement attendu sur "
          f"{int(args.ages[0])} ans : {abs(show[col0].max()-show[col0].min()):.1f} ppm "
          f"(A1 bornait l'ecart observe a 70 ppm)")

    # ---------- C. ΔGC DEJA ACCUMULE DEPUIS LE MRCA -----------------------------
    div = {p: chr(anc[p]) for p in range(len(seq))
           if anc[p] not in (ord("N"), ord(seq[p])) and p not in masked}
    print(f"\n# {len(div)} positions non masquees ou MTBC0 diverge de H37Rv",
          file=sys.stderr)
    rows = []
    for lin in lignees:
        for s in chemins(lin, pools[lin]):
            subs = read_subs(s / "NC_000962.3" / "spdi.txt")
            if not subs:
                continue
            L = G = 0
            vus = set()
            for pos, ref, alt in subs:
                if pos in masked:
                    continue
                a = chr(anc[pos]) if pos < len(anc) else "N"
                if a == "N":
                    continue
                if pos in div:
                    vus.add(pos)
                    if a == alt:
                        continue                  # la souche est ancestrale
                    ref = a
                f = flux(ref, alt)
                L += f == "loss"
                G += f == "gain"
            for pos, a in div.items():            # H37Rv derive, souche non variante
                if pos in vus:
                    continue
                f = flux(a, seq[pos])
                L += f == "loss"
                G += f == "gain"
            rows.append(dict(lignee=lin, sra=s.name, loss=L, gain=G,
                             dgc_ppm=1e6 * (G - L) / n_sites))
    tc = pd.DataFrame(rows)
    tc.to_csv(args.out_prefix + "dgc_depuis_mrca.tsv", sep="\t", index=False)
    rows = []
    fr = tb[tb.horloge == "central"].set_index("lignee").frac_gc_changing
    for lin, d in tc.groupby("lignee"):
        x = d["dgc_ppm"].to_numpy()
        bs = np.array([x[rng.integers(0, len(x), len(x))].mean()
                       for _ in range(args.boot)])
        lo, hi = np.percentile(bs, [2.5, 97.5])
        age = (d["loss"].mean() + d["gain"].mean()) / fr[lin] / (args.clock * n_sites)
        rows.append(dict(lignee=lin, n=len(d), loss=d["loss"].mean(),
                         gain=d["gain"].mean(), dgc_ppm=x.mean(),
                         sd=x.std(ddof=1), boot_lo=lo, boot_hi=hi,
                         age_branche_ans=age, sous_puissante=lin in FAIBLES))
    g = pd.DataFrame(rows).sort_values("dgc_ppm")
    g.to_csv(args.out_prefix + "dgc_par_lignee.tsv", sep="\t", index=False)
    print("\n=== C. ΔGC DEJA ACCUMULE depuis le MRCA du MTBC (branche entiere, "
          "polarisee) ===")
    print(g.round(2).to_string(index=False))
    print(f"  etendue inter-lignees du ΔGC moyen : "
          f"{g.dgc_ppm.max() - g.dgc_ppm.min():.1f} ppm "
          f"(de {g.dgc_ppm.min():.1f} a {g.dgc_ppm.max():.1f})")
    gf = g[~g.sous_puissante]
    print(f"  hors sous-puissantes : {gf.dgc_ppm.max() - gf.dgc_ppm.min():.1f} ppm "
          f"(de {gf.dgc_ppm.min():.1f} a {gf.dgc_ppm.max():.1f})")
    print(f"  duree impliquee de la branche MRCA -> souche, sous l'horloge "
          f"{args.clock:.0e} : {g.age_branche_ans.min():.0f} a "
          f"{g.age_branche_ans.max():.0f} ans "
          f"(mediane {g.age_branche_ans.median():.0f})")
    m = g.merge(ta[["lignee", "S", "vu"]], on="lignee")
    rho, p = stats.spearmanr(m.vu, m.dgc_ppm)
    print(f"  Spearman(v/u terminal, ΔGC branche entiere) = {rho:.3f} "
          f"(p = {p:.2g}) sur {len(m)} lignees -- controle de coherence interne, "
          f"les deux mesures portent sur des jeux d'evenements largement disjoints")
    mf = m[~m.sous_puissante]
    rho2, p2 = stats.spearmanr(mf.vu, mf.dgc_ppm)
    print(f"    hors sous-puissantes ({len(mf)} lignees) : rho = {rho2:.3f} "
          f"(p = {p2:.2g})")

    # ---------- D. EFFICACITE DE LA SELECTION -----------------------------------
    owner, off, strand, cds = load_cds(len(seq))
    print(f"# {len(cds)} CDS charges", file=sys.stderr)
    opp = opportunite_codante(seq, masked, owner, off, strand, cds)
    print(f"# opportunite codante non masquee : {sum(v[0] for v in opp.values())} "
          f"syn, {sum(v[1] for v in opp.values())} nonsyn", file=sys.stderr)

    rows = []
    for lin in lignees:
        subsets = []
        for s in chemins(lin, pools[lin]):
            v = read_subs(s / "NC_000962.3" / "spdi.txt")
            if v:
                subsets.append(v)
        support = defaultdict(int)
        for i, subs in enumerate(subsets):
            for v in subs:
                support[v] |= 1 << i
        obs = Counter()          # (ref, alt) -> compte, sur branches terminales
        n_syn = n_non = 0
        pi_num = 0
        n = len(subsets)
        for (pos, ref, alt), mask in support.items():
            if pos in masked:
                continue
            k = bin(mask).count("1")
            pi_num += 2 * k * (n - k)
            if k != 1:
                continue
            a = chr(anc[pos]) if pos < len(anc) else "N"
            if a != ref:          # non liftee, inverse, ou tierce : ecarte ici
                continue
            eff = effect(seq, pos, alt, owner, off, strand, cds)
            if eff is None:
                continue
            obs[(ref, alt)] += 1
            n_syn += eff == "syn"
            n_non += eff == "nonsyn"
        tot = sum(obs.values())
        if not tot or not n_syn or not n_non:
            print(f"# {lin} : pas assez de singletons codants "
                  f"({n_syn} syn, {n_non} nonsyn), ecartee de D", file=sys.stderr)
            continue
        exp_s = sum(obs[k] / tot * opp.get(k, (0, 0))[0] for k in obs)
        exp_n = sum(obs[k] / tot * opp.get(k, (0, 0))[1] for k in obs)
        pnps = (n_non / n_syn) / (exp_n / exp_s)
        se = np.sqrt(1 / n_non + 1 / n_syn)
        rows.append(dict(lignee=lin, n=n, n_syn=n_syn, n_nonsyn=n_non,
                         attendu_N_sur_S=exp_n / exp_s,
                         observe_N_sur_S=n_non / n_syn, pn_ps=pnps,
                         pn_ps_lo=pnps * np.exp(-1.96 * se),
                         pn_ps_hi=pnps * np.exp(1.96 * se),
                         branche_term=(n_syn + n_non) / n,
                         pi=pi_num / (n * (n - 1)) / n_sites))
    td = pd.DataFrame(rows).merge(ta[["lignee", "S", "vu", "gc_eq",
                                      "sous_puissante"]], on="lignee")
    td = td.sort_values("S", ascending=False)
    td.to_csv(args.out_prefix + "selection.tsv", sep="\t", index=False)
    print("\n=== D. EFFICACITE DE LA SELECTION : la force requise suit-elle un "
          "proxy de Ne ? ===")
    print("pN/pS sur branches terminales, CORRIGE du spectre mutationnel propre "
          "a chaque lignee\n(l'attendu N/S est recalcule sous les frequences de "
          "classes de substitution observees chez elle).")
    aff = td[["lignee", "n", "n_syn", "n_nonsyn", "attendu_N_sur_S",
              "observe_N_sur_S", "pn_ps", "pn_ps_lo", "pn_ps_hi", "branche_term",
              "S", "sous_puissante"]].copy()
    aff["pi_e4"] = td.pi * 1e4
    print(aff.round(4).to_string(index=False))

    def correlations(t, etiquette):
        print(f"\n  -- {etiquette} (n = {len(t)} lignees)")
        for x, lab in (("pn_ps", "pN/pS corrige"), ("pi", "diversite pi")):
            rho, p = stats.spearmanr(t.S, t[x])
            print(f"    Spearman(S requis, {lab}) = {rho:+.3f} (p = {p:.3g})")
        rho_b, p_b = stats.spearmanr(t.pn_ps, t.branche_term)
        rho_sb, p_sb = stats.spearmanr(t.S, t.branche_term)
        print(f"    CONFONDANT longueur de branche terminale "
              f"({t.branche_term.min():.0f} a {t.branche_term.max():.0f} subs "
              f"codantes/souche) :")
        rho_pb, p_pb = stats.spearmanr(t.pi, t.branche_term)
        print(f"      Spearman(pN/pS, longueur) = {rho_b:+.3f} (p = {p_b:.3g}) ; "
              f"Spearman(S, longueur) = {rho_sb:+.3f} (p = {p_sb:.3g}) ; "
              f"Spearman(pi, longueur) = {rho_pb:+.3f} (p = {p_pb:.3g})")
        r = {c: stats.rankdata(t[c]) for c in ("S", "pn_ps", "pi", "branche_term")}

        def resid(y, x):
            return y - np.polyval(np.polyfit(x, y, 1), x)
        # La longueur de branche terminale est une mesure de la DENSITE
        # d'echantillonnage du pool : les deux proxys de Ne en heritent, il faut
        # donc la tenir constante pour chacun des deux, pas seulement pour pN/pS.
        out = dict(jeu=etiquette, n=len(t))
        for x, lab in (("pn_ps", "pN/pS"), ("pi", "pi")):
            pr, pp2 = stats.pearsonr(resid(r["S"], r["branche_term"]),
                                     resid(r[x], r["branche_term"]))
            print(f"      correlation PARTIELLE de rang S x {lab} a longueur "
                  f"tenue constante : {pr:+.3f} (p = {pp2:.3g}, "
                  f"ddl = {len(t)-3})")
            out[f"partielle_S_{x}"] = pr
            out[f"p_partielle_S_{x}"] = pp2
        return out

    print("  Attendu si S = 2.Ne.s a s constant : correlation NEGATIVE avec "
          "pN/pS (fort Ne = selection efficace = pN/pS bas) et POSITIVE avec pi.")
    partielles = pd.DataFrame([
        correlations(td, "panel complet"),
        correlations(td[~td.sous_puissante],
                     "hors les quatre sous-puissantes de A23")])
    partielles.to_csv(args.out_prefix + "partielles.tsv", sep="\t", index=False)

    # ---------- E. COMPARAISON AVEC LE PANEL n = 11 ------------------------------
    try:
        old_f = pd.read_csv(args.ancien_force, sep="\t").rename(
            columns={"clade": "lignee"})
        old_d = pd.read_csv(args.ancien_selection, sep="\t").rename(
            columns={"clade": "lignee"})
    except FileNotFoundError as e:
        print(f"\n# comparaison n=11 impossible : {e}", file=sys.stderr)
        return
    cmp_rows = []
    for new, old in MEME_MAILLE.items():
        a = ta[ta.lignee == new]
        b = old_f[old_f.lignee == old]
        if a.empty or b.empty:
            continue
        r = dict(lignee=new, S_n11=b.S.iloc[0], S_n17=a.S.iloc[0],
                 ecart_S_pct=100 * (a.S.iloc[0] - b.S.iloc[0]) / b.S.iloc[0],
                 vu_n11=b.vu.iloc[0], vu_n17=a.vu.iloc[0])
        c, d = td[td.lignee == new], old_d[old_d.lignee == old]
        if not c.empty and not d.empty:
            r.update(pn_ps_n11=d.pn_ps.iloc[0], pn_ps_n17=c.pn_ps.iloc[0])
        cmp_rows.append(r)
    tcmp = pd.DataFrame(cmp_rows)
    tcmp.to_csv(args.out_prefix + "comparaison_n11.tsv", sep="\t", index=False)
    print("\n=== E. CE QUE LE PASSAGE DE 11 A 17 LIGNEES CHANGE ===")
    print("  lignees lues a la MEME maille des deux cotes (L1 y est complete a "
          "droite, amputee de 15,5 % a gauche) :")
    print(tcmp.round(4).to_string(index=False))
    o = old_f.sort_values("S", ascending=False)
    print(f"  etendue de S : n=11 {o.S.min():.3f}-{o.S.max():.3f} "
          f"(x{o.S.max()/o.S.min():.2f}) contre n=17 {ta.S.min():.3f}-"
          f"{ta.S.max():.3f} (x{ta.S.max()/ta.S.min():.2f})")
    rho_o, p_o = stats.spearmanr(old_d.S, old_d.pn_ps)
    rho_n, p_n = stats.spearmanr(td.S, td.pn_ps)
    print(f"  Spearman(S, pN/pS) : n=11 {rho_o:+.3f} (p = {p_o:.3g}) contre "
          f"n={len(td)} {rho_n:+.3f} (p = {p_n:.3g})")

    # ---------- F. DIAGNOSTIC : d'ou vient l'ecart entre les deux verdicts ? -----
    # La question n'est pas rhetorique. Si A20 s'effondre en passant de onze a
    # dix-sept points, il faut savoir si c'est (i) la L1 amputee, (ii) le
    # changement de MAILLE (le panel n=11 lisait L2.2.1, L4.1.2, L4.3, L6.1 la ou
    # le panel n=17 lit L2, L4, L6), (iii) l'entree des lignees nouvelles, ou
    # (iv) la fragilite intrinseque d'un rho de rang a onze points.
    print("\n=== F. DIAGNOSTIC de l'ecart entre A20 (n=11) et le panel elargi ===")
    diag = []

    def rho_of(t, x="pn_ps"):
        r, p = stats.spearmanr(t.S, t[x])
        return r, p

    # (i) la L1 amputee, seule difference de valeur entre les deux panels
    fix = old_d.copy()
    l1 = td[td.lignee == "L1"]
    if not l1.empty:
        fix.loc[fix.lignee == "L1", ["S", "pn_ps", "pi"]] = [
            l1.S.iloc[0], l1.pn_ps.iloc[0], l1.pi.iloc[0]]
    for lab, t in (("n=11 tel que publie en A20", old_d),
                   ("n=11 avec L1 complete", fix),
                   ("n=17 panel elargi", td),
                   ("n=13 bien pourvues", td[~td.sous_puissante])):
        r, p = rho_of(t)
        rp, pp = rho_of(t, "pi")
        diag.append(dict(jeu=lab, n=len(t), rho_S_pnps=r, p_S_pnps=p,
                         rho_S_pi=rp, p_S_pi=pp))
    d = pd.DataFrame(diag)
    print(d.round(4).to_string(index=False))

    # (ii) et (iii) : leave-one-out sur le panel elargi, pour voir si l'absence de
    # correlation tient a quelques points ou est generale
    loo = []
    for lin in td.lignee:
        t = td[td.lignee != lin]
        r, p = rho_of(t)
        loo.append(dict(lignee_retiree=lin, rho=r, p=p))
    lo = pd.DataFrame(loo).sort_values("rho")
    print("\n  leave-one-out sur les dix-sept (rho de S x pN/pS sans la lignee) :")
    print(lo.round(4).to_string(index=False))
    print(f"    rho reste dans [{lo.rho.min():+.3f}, {lo.rho.max():+.3f}] : "
          f"aucune lignee unique ne restaure la correlation de A20")

    # (iv) fragilite d'un rho a onze points : quelle etendue de rho obtient-on en
    # tirant au hasard onze des dix-sept lignees ?
    idx = np.arange(len(td))
    sub = []
    for _ in range(2000):
        k = rng.choice(idx, size=11, replace=False)
        r, _ = rho_of(td.iloc[k])
        sub.append(r)
    sub = np.array(sub)
    print(f"\n  en tirant 11 lignees au hasard parmi les 17 (2 000 tirages) : "
          f"rho median {np.median(sub):+.3f}, "
          f"IC {np.percentile(sub, 2.5):+.3f} a {np.percentile(sub, 97.5):+.3f} ; "
          f"P(rho <= -0.736) = {(sub <= rho_o).mean():.4f}")

    # confondant de longueur de branche, dans les deux panels
    for lab, t in (("n=11", old_d), ("n=17", td),
                   ("n=13", td[~td.sous_puissante])):
        rb, pb = stats.spearmanr(t.S, t.branche_term)
        print(f"  Spearman(S, longueur de branche terminale) sur {lab} : "
              f"{rb:+.3f} (p = {pb:.3g})")
    d.to_csv(args.out_prefix + "diagnostic_a20.tsv", sep="\t", index=False)
    lo.to_csv(args.out_prefix + "loo_a20.tsv", sep="\t", index=False)

    print(f"\n# ecrits : {args.out_prefix}[force_requise|cinetique|"
          f"dgc_depuis_mrca|dgc_par_lignee|selection|comparaison_n11].tsv",
          file=sys.stderr)


if __name__ == "__main__":
    main()
