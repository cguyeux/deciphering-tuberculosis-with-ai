#!/usr/bin/env python3
"""
Objet       : P9.3 -- confronter le GC d'EQUILIBRE mutationnel (A15 : 32 a 54 %
              selon la lignee) au GC OBSERVE (65,4 %, identique a 70 ppm pres
              entre lignees, A1). Quatre sections.
              (A) FORCE REQUISE : sous l'equilibre mutation-selection-derive de
              Li/Bulmer, la frequence attendue de GC vaut k.e^S/(1+k.e^S) avec
              k = u/v ; le coefficient d'echelle requis pour tenir 65,4 % vaut
              donc S = ln[(GC/(1-GC)).(v/u)]. On le calcule par lignee, avec IC.
              (B) CINETIQUE : la composition ne repond pas instantanement. En
              convertissant v/u en taux absolus par l'horloge du MTBC, on obtient
              la constante de relaxation u+v, la demi-vie de l'approche a
              l'equilibre, et le deplacement de GC attendu sur l'age des lignees.
              (C) MESURE DIRECTE, qui reconcilie A1 et A15 : le ΔGC deja
              accumule depuis le MRCA du MTBC, branche entiere, souche par
              souche, polarise sur MTBC0. C'est la quantite que A1 bornait a
              70 ppm et que A15 predit desormais.
              (D) EFFICACITE DE LA SELECTION : la force requise varie-t-elle
              comme un proxy de Ne ? pN/pS sur branches terminales CORRIGE du
              spectre mutationnel propre a chaque lignee (sans quoi une lignee
              plus riche en transitions parait sous selection differente pour
              une raison purement mutationnelle), et diversite nucleotidique pi.
Entrees     : résultats/phase3_counts_par_souche_n40.tsv (L, G, neutre par souche)
              bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt
              data/MTBC0/ancestral_on_H37Rv.bin, Mask.files/*.bed
              investigate_phylo/resources/NC_000962.3.gff3 (CDS H37Rv)
Sorties     : résultats/phase4_p93_force_requise.tsv
              résultats/phase4_p93_cinetique.tsv
              résultats/phase4_p93_dgc_depuis_mrca.tsv
              résultats/phase4_p93_selection.tsv
Reutilisable: oui -- (A) et (B) valent pour toute bacterie dont on a mesure v/u ;
              (D) est un pN/pS corrige du spectre, reutilisable tel quel
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
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import COMPL, read_subs, strain_dirs, flux  # noqa: E402
from phase2_polarisation_mtbc0 import (build_ancestral, load_mask, read_fasta,  # noqa: E402
                                       H37RV)
from phase3_counts_par_souche import MASKS  # noqa: E402
from phase3_sueoka_gc_eq import opportunity, vu_ci  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
GFF3 = Path("/home/christophe/docs/codes/mtbc/investigate_phylo/resources/"
            "NC_000962.3.gff3")
CLADES = ["L1", "L2.2.1", "L3", "L4.1.2", "L4.3", "L6.1", "L7", "L9",
          "Orygis_La3", "Caprae_La2", "Microti"]
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[i * 16 + j * 4 + k]
        for i, a in enumerate(BASES)
        for j, b in enumerate(BASES)
        for k, c in enumerate(BASES)}


def revcomp(s):
    return "".join(COMPL[c] for c in reversed(s))


def load_cds(seq_len):
    """Pour chaque position 0-based : (index_cds, offset_dans_le_codon, brin),
    ou None. Un CDS unique par position (le premier rencontre en cas de
    chevauchement, cas marginal chez H37Rv)."""
    owner = np.full(seq_len, -1, np.int32)
    off = np.zeros(seq_len, np.int8)
    strand = np.zeros(seq_len, np.int8)
    cds = []
    for line in GFF3.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 8 or f[2] != "CDS":
            continue
        s, e, st = int(f[3]) - 1, int(f[4]), f[6]
        if (e - s) % 3:
            continue
        idx = len(cds)
        cds.append((s, e, st))
        for p in range(s, e):
            if owner[p] != -1:
                continue
            owner[p] = idx
            off[p] = (p - s) % 3 if st == "+" else (e - 1 - p) % 3
            strand[p] = 1 if st == "+" else -1
    return owner, off, strand, cds


def effect(seq, pos, alt, owner, off, strand, cds):
    """« syn » / « nonsyn » / None (hors CDS). SPDI : ref et alt sont sur le brin
    +, on substitue directement dans le codon forward puis on revcomp si besoin
    (piege documente dans la KB : ne PAS complementer alt)."""
    i = owner[pos]
    if i < 0:
        return None
    o = int(off[pos])
    start = pos - o if strand[pos] == 1 else pos - (2 - o)
    s, e, _ = cds[i]
    if start < s or start + 3 > e:
        return None
    cod = seq[start:start + 3]
    mut = cod[:pos - start] + alt + cod[pos - start + 1:]
    if strand[pos] == -1:
        cod, mut = revcomp(cod), revcomp(mut)
    a1, a2 = CODE.get(cod), CODE.get(mut)
    if a1 is None or a2 is None:
        return None
    return "syn" if a1 == a2 else "nonsyn"


def cds_opportunity(seq, masked, owner, off, strand, cds):
    """opp[(ref, alt)] = (n_syn, n_nonsyn) sur tout le codant non masque."""
    opp = defaultdict(lambda: [0, 0])
    for pos in range(len(seq)):
        if owner[pos] < 0 or pos in masked:
            continue
        ref = seq[pos]
        if ref not in COMPL:
            continue
        for alt in "ACGT":
            if alt == ref:
                continue
            eff = effect(seq, pos, alt, owner, off, strand, cds)
            if eff:
                opp[(ref, alt)][0 if eff == "syn" else 1] += 1
    return {k: tuple(v) for k, v in opp.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--n-per-clade", type=int, default=40)
    ap.add_argument("--counts", default=str(ROOT / "résultats" /
                                            "phase3_counts_par_souche_n40.tsv"))
    ap.add_argument("--clock", type=float, default=1e-7,
                    help="horloge MTBC en substitutions/site/an (Menardo 2019)")
    ap.add_argument("--clock-lo", type=float, default=5e-8)
    ap.add_argument("--clock-hi", type=float, default=2e-7)
    ap.add_argument("--ages", type=float, nargs="*", default=[6000, 20000, 1e5, 1e6])
    ap.add_argument("--out-prefix", default=str(ROOT / "résultats" / "phase4_p93_"))
    args = ap.parse_args()

    seq = read_fasta(H37RV)
    anc = build_ancestral()
    masked = load_mask(MASKS, in_mtbc0_coords=True)
    (n_gc, n_at), _ = opportunity()
    n_sites = n_gc + n_at
    gc_obs = n_gc / n_sites

    c = pd.read_csv(args.counts, sep="\t")
    agg = c.groupby("clade")[["loss", "gain", "neutral"]].sum().reset_index()

    # ---------- A. FORCE REQUISE (Li 1987 / Bulmer 1991) ----------------------
    rows = []
    for _, r in agg.iterrows():
        vu, lo, hi = vu_ci(r["loss"], r["gain"], n_gc, n_at)
        f = gc_obs / (1 - gc_obs)
        rows.append(dict(clade=r["clade"], vu=vu, gc_eq=1 / (1 + vu),
                         S=np.log(f * vu), S_lo=np.log(f * lo), S_hi=np.log(f * hi),
                         ecart_gc_pts=100 * (gc_obs - 1 / (1 + vu))))
    ta = pd.DataFrame(rows).sort_values("S", ascending=False)
    ta.to_csv(args.out_prefix + "force_requise.tsv", sep="\t", index=False)
    print(f"\n=== A. FORCE REQUISE pour tenir GC = {100*gc_obs:.3f} % ===")
    print("S = ln[(GC/(1-GC)) . v/u] : coefficient d'echelle (2.Ne.s en haploide) "
          "que devrait avoir\n    toute force favorisant G:C pour que 65,4 % soit "
          "un equilibre mutation-selection-derive.")
    print(ta.round(4).to_string(index=False))
    r_ = ta.S.max() / ta.S.min()
    print(f"  etendue : S = {ta.S.min():.3f} ({ta.clade.iloc[-1]}) a "
          f"{ta.S.max():.3f} ({ta.clade.iloc[0]}), rapport x{r_:.2f}")
    print(f"  en probabilite de fixation relative, e^S va de {np.exp(ta.S.min()):.2f} "
          f"a {np.exp(ta.S.max()):.2f}")

    # ---------- B. CINETIQUE ---------------------------------------------------
    rows = []
    for _, r in agg.iterrows():
        L, G, N = r["loss"], r["gain"], r["neutral"]
        frac_gc_changing = (L + G) / (L + G + N)
        for tag, clk in (("central", args.clock), ("lo", args.clock_lo),
                         ("hi", args.clock_hi)):
            R = clk * frac_gc_changing                 # subs GC-changeantes/site/an
            v = R * (L / (L + G)) * n_sites / n_gc     # par site G:C et par an
            u = R * (G / (L + G)) * n_sites / n_at     # par site A:T et par an
            lam = u + v
            gceq = u / (u + v)
            d = dict(clade=r["clade"], horloge=tag, clock=clk,
                     frac_gc_changing=frac_gc_changing, u=u, v=v, lambda_=lam,
                     demi_vie_ans=np.log(2) / lam, gc_eq=gceq)
            d["subs_par_genome_par_an"] = clk * n_sites
            for t in args.ages:
                # signe negatif = perte de GC (la composition descend vers gceq)
                d[f"dGC_ppm_{int(t)}ans"] = 1e6 * (gceq - gc_obs) * -np.expm1(-lam * t)
            rows.append(d)
    tb = pd.DataFrame(rows)
    tb.to_csv(args.out_prefix + "cinetique.tsv", sep="\t", index=False)
    print("\n=== B. CINETIQUE : a quelle vitesse la composition repondrait-elle ? ===")
    show = tb[tb.horloge == "central"].copy()
    show["u_e8"] = show.u * 1e8
    show["v_e8"] = show.v * 1e8
    show["demi_vie_Man"] = show.demi_vie_ans / 1e6
    print(show[["clade", "u_e8", "v_e8", "demi_vie_Man"] +
               [f"dGC_ppm_{int(t)}ans" for t in args.ages]]
          .round(2).to_string(index=False))
    print(f"  u_e8, v_e8 : taux par site A:T et par site G:C, en 1e-8/an ; "
          f"demi_vie_Man en millions d'annees")
    print(f"  dGC en ppm, signe NEGATIF = perte de GC ; horloge "
          f"{args.clock:.0e} subs/site/an, soit "
          f"{args.clock*n_sites:.2f} substitution/genome/an")
    dmin, dmax = show[f"dGC_ppm_{int(args.ages[0])}ans"].min(), \
        show[f"dGC_ppm_{int(args.ages[0])}ans"].max()
    print(f"  etendue INTER-LIGNEES du deplacement attendu sur "
          f"{int(args.ages[0])} ans : {abs(dmax-dmin):.1f} ppm "
          f"(A1 bornait l'ecart observe a 70 ppm)")

    # ---------- C. ΔGC DEJA ACCUMULE DEPUIS LE MRCA -----------------------------
    div = {p: chr(anc[p]) for p in range(len(seq))
           if anc[p] not in (ord("N"), ord(seq[p])) and p not in masked}
    print(f"\n# {len(div)} positions non masquees ou MTBC0 diverge de H37Rv",
          file=sys.stderr)
    rows = []
    for clade in CLADES:
        st = strain_dirs(clade)
        rng = random.Random(0)
        rng.shuffle(st)
        for s in st[:args.n_per_clade]:
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
            rows.append(dict(clade=clade, sra=s.name, loss=L, gain=G,
                             dgc_ppm=1e6 * (G - L) / n_sites))
    tc = pd.DataFrame(rows)
    tc.to_csv(args.out_prefix + "dgc_depuis_mrca.tsv", sep="\t", index=False)
    g = tc.groupby("clade").agg(n=("sra", "size"), loss=("loss", "mean"),
                                gain=("gain", "mean"),
                                dgc_ppm=("dgc_ppm", "mean"),
                                sd=("dgc_ppm", "std")).reset_index()
    g = g.sort_values("dgc_ppm")
    print("\n=== C. ΔGC DEJA ACCUMULE depuis le MRCA du MTBC (branche entiere, "
          "polarisee) ===")
    print(g.round(2).to_string(index=False))
    print(f"  etendue inter-lignees du ΔGC moyen : "
          f"{g.dgc_ppm.max() - g.dgc_ppm.min():.1f} ppm "
          f"(de {g.dgc_ppm.min():.1f} a {g.dgc_ppm.max():.1f})")
    fr = tb[tb.horloge == "central"].set_index("clade").frac_gc_changing
    age = ((g.set_index("clade").loss + g.set_index("clade").gain) /
           fr / (args.clock * n_sites))
    print(f"  duree impliquee de la branche MRCA -> souche, sous l'horloge "
          f"{args.clock:.0e} : {age.min():.0f} a {age.max():.0f} ans "
          f"(mediane {age.median():.0f}) -- c'est l'echelle a laquelle lire la "
          f"colonne dGC_ppm_{int(args.ages[0])}ans de la section B")
    m = g.merge(ta[["clade", "S", "vu"]], on="clade")
    rho, p = stats.spearmanr(m.vu, m.dgc_ppm)
    print(f"  Spearman(v/u terminal, ΔGC branche entiere) = {rho:.3f} "
          f"(p = {p:.2g}) -- controle de coherence interne, les deux mesures "
          f"portent sur des jeux d'evenements largement disjoints")

    # ---------- D. EFFICACITE DE LA SELECTION -----------------------------------
    owner, off, strand, cds = load_cds(len(seq))
    print(f"# {len(cds)} CDS charges ; opportunite codante en cours...",
          file=sys.stderr)
    opp = cds_opportunity(seq, masked, owner, off, strand, cds)
    tot_syn = sum(v[0] for v in opp.values())
    tot_non = sum(v[1] for v in opp.values())
    print(f"# opportunite codante non masquee : {tot_syn} syn, {tot_non} nonsyn",
          file=sys.stderr)

    rows = []
    for clade in CLADES:
        st = strain_dirs(clade)
        rng = random.Random(0)
        rng.shuffle(st)
        sample = st[:args.n_per_clade]
        subsets, names = [], []
        for s in sample:
            v = read_subs(s / "NC_000962.3" / "spdi.txt")
            if v:
                subsets.append(v)
                names.append(s.name)
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
        exp_s = sum(obs[k] / tot * opp.get(k, (0, 0))[0] for k in obs)
        exp_n = sum(obs[k] / tot * opp.get(k, (0, 0))[1] for k in obs)
        pnps = (n_non / n_syn) / (exp_n / exp_s)
        se = np.sqrt(1 / n_non + 1 / n_syn)
        rows.append(dict(clade=clade, n=n, n_syn=n_syn, n_nonsyn=n_non,
                         attendu_N_sur_S=exp_n / exp_s,
                         observe_N_sur_S=n_non / n_syn, pn_ps=pnps,
                         pn_ps_lo=pnps * np.exp(-1.96 * se),
                         pn_ps_hi=pnps * np.exp(1.96 * se),
                         branche_term=(n_syn + n_non) / n,
                         pi=pi_num / (n * (n - 1)) / n_sites))
    td = pd.DataFrame(rows).merge(ta[["clade", "S", "vu", "gc_eq"]], on="clade")
    td = td.sort_values("S", ascending=False)
    td.to_csv(args.out_prefix + "selection.tsv", sep="\t", index=False)
    print("\n=== D. EFFICACITE DE LA SELECTION : la force requise suit-elle un "
          "proxy de Ne ? ===")
    print("pN/pS sur branches terminales, CORRIGE du spectre mutationnel propre "
          "a chaque lignee\n(l'attendu N/S est recalcule sous les frequences de "
          "classes de substitution observees chez elle).")
    aff = td[["clade", "n_syn", "n_nonsyn", "attendu_N_sur_S", "observe_N_sur_S",
              "pn_ps", "pn_ps_lo", "pn_ps_hi", "branche_term", "S"]].copy()
    aff["pi_e4"] = td.pi * 1e4
    print(aff.round(4).to_string(index=False))
    for x, lab in (("pn_ps", "pN/pS corrige"), ("pi", "diversite pi")):
        rho, p = stats.spearmanr(td.S, td[x])
        print(f"  Spearman(S requis, {lab}) = {rho:+.3f} (p = {p:.3g})")
    print("  Attendu si S = 2.Ne.s a s constant : correlation NEGATIVE avec "
          "pN/pS (fort Ne = selection efficace = pN/pS bas)\n  et POSITIVE avec pi.")

    # CONFONDANT : une lignee densement echantillonnee a des branches terminales
    # COURTES, donc des mutations JEUNES que la selection purificatrice n'a pas
    # encore eliminees, donc un pN/pS mecaniquement plus haut. Il faut verifier
    # que la correlation ci-dessus n'est pas ce seul effet.
    rho_b, p_b = stats.spearmanr(td.pn_ps, td.branche_term)
    rho_sb, p_sb = stats.spearmanr(td.S, td.branche_term)
    print(f"\n  CONFONDANT teste -- longueur des branches terminales "
          f"(subs codantes par souche, {td.branche_term.min():.0f} a "
          f"{td.branche_term.max():.0f}) :")
    print(f"    Spearman(pN/pS, longueur de branche) = {rho_b:+.3f} (p = {p_b:.3g})")
    print(f"    Spearman(S requis, longueur de branche) = {rho_sb:+.3f} "
          f"(p = {p_sb:.3g})")
    r = {c: stats.rankdata(td[c]) for c in ("S", "pn_ps", "branche_term")}
    def resid(y, x):
        b = np.polyfit(x, y, 1)
        return y - np.polyval(b, x)
    pr, pp2 = stats.pearsonr(resid(r["S"], r["branche_term"]),
                             resid(r["pn_ps"], r["branche_term"]))
    print(f"    correlation PARTIELLE de rang S x pN/pS a longueur de branche "
          f"tenue constante : {pr:+.3f} (p = {pp2:.3g}, ddl = {len(td)-3})")

    print(f"\n# ecrit dans {args.out_prefix}"
          f"[force_requise|cinetique|dgc_depuis_mrca|selection].tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
