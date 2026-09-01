#!/usr/bin/env python3
"""
Objet       : P2.1 -- verification croisee INDEPENDANTE de la polarisation par
              MTBC0 (P2.4), en reconstruisant l'etat ancestral par l'outgroup
              *M. canettii*. Logique de parcimonie a un outgroup : a un site
              polymorphe DANS le MTBC (deux alleles, ref H37Rv et alt), l'allele
              porte par canettii est l'ancestral. Le geste est independant de
              MTBC0 par sa source (149 genomes de canettii de la base locale
              contre un genome ancestral impute par ML sur le MTBC) et par sa
              logique (outgroup contre reconstruction interne). Trois sorties :
              une matrice de confusion des deux verdicts d'orientation, le test
              a preuve POSITIVE sur les positions ou MTBC0 diverge de H37Rv
              (celles qui portaient tout l'artefact), et le v/u par lignee
              recalcule sous polarisation canettii.
              LIMITE assumee : chez canettii, absence de variant est lue comme
              « porte l'allele H37Rv », ce qui confond identite et non-couverture
              (le piege de P3.3). Le test a preuve positive y echappe, la matrice
              de confusion non ; les deux sont rapportes separement pour cette
              raison.
Entrees     : bdd/actuelle/Canettii/<SRA>/NC_000962.3/spdi.txt (149 genomes)
              bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt (11 lignees)
              data/MTBC0/ancestral_on_H37Rv.bin (via phase2)
Sorties     : résultats/phase4_p21_confusion.tsv
              résultats/phase4_p21_vu_canettii.tsv
              cache data/canettii_consensus.tsv.gz (allele canettii par position
              H37Rv ou au moins une souche canettii porte un variant)
Reutilisable: oui -- la table consensus canettii sert a toute polarisation MTBC
              par outgroup, independamment du GC
Projet      : GC_par_lignee
Date        : 2026-08-29
"""
import argparse
import gzip
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import BDD, read_subs, strain_dirs, flux  # noqa: E402
from phase2_polarisation_mtbc0 import (build_ancestral, load_mask, read_fasta,  # noqa: E402
                                       H37RV)
from phase3_counts_par_souche import MASKS  # noqa: E402
from phase3_sueoka_gc_eq import opportunity, vu_ci  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "canettii_consensus.tsv.gz"
CLADES = ["L1", "L2.2.1", "L3", "L4.1.2", "L4.3", "L6.1", "L7", "L9",
          "Orygis_La3", "Caprae_La2", "Microti"]


def canettii_consensus(force=False):
    """{pos_0based: (allele_majoritaire, n_porteurs, n_canettii_total)} pour
    toute position ou au moins une souche canettii porte une substitution
    simple. Les positions absentes de la table sont celles ou aucune canettii
    n'a de variant : lues comme « allele H37Rv », avec la reserve ci-dessus."""
    strains = sorted(p for p in (BDD / "Canettii").iterdir()
                     if p.is_dir() and (p / "NC_000962.3" / "spdi.txt").exists())
    if CACHE.exists() and not force:
        tab = {}
        with gzip.open(CACHE, "rt") as fh:
            n_tot = int(next(fh).split("=")[1])
            next(fh)
            for line in fh:
                p, a, k = line.split("\t")
                tab[int(p)] = (a, int(k))
        return tab, n_tot
    votes = defaultdict(Counter)
    for s in strains:
        for pos, _ref, alt in read_subs(s / "NC_000962.3" / "spdi.txt"):
            votes[pos][alt] += 1
    tab = {p: c.most_common(1)[0] for p, c in votes.items()}
    with gzip.open(CACHE, "wt") as fh:
        fh.write(f"# n_canettii={len(strains)}\n")
        fh.write("pos0\tallele\tn_porteurs\n")
        for p in sorted(tab):
            fh.write(f"{p}\t{tab[p][0]}\t{tab[p][1]}\n")
    print(f"# consensus canettii : {len(strains)} genomes, {len(tab)} positions "
          f"variables", file=sys.stderr)
    return tab, len(strains)


def canettii_allele(pos, h37, cons, n_can, seuil):
    """Allele canettii a une position H37Rv 0-based, ou None si le site est
    polymorphe DANS canettii (verdict refuse plutot que force)."""
    hit = cons.get(pos)
    if hit is None:
        return h37[pos]                       # aucune canettii ne varie
    allele, k = hit
    if k >= seuil * n_can:
        return allele                         # quasi-fixe chez canettii
    if k <= (1 - seuil) * n_can:
        return h37[pos]                       # variant marginal chez canettii
    return None                               # canettii polymorphe : ambigu


def verdict(anc_allele, ref, alt):
    if anc_allele is None:
        return "ambigu"
    if anc_allele == ref:
        return "correcte"
    if anc_allele == alt:
        return "inverse"
    return "tierce"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--n-per-clade", type=int, default=40)
    ap.add_argument("--seuil", type=float, default=0.9,
                    help="fraction de canettii requise pour un allele quasi-fixe")
    ap.add_argument("--force", action="store_true", help="recalculer le consensus")
    ap.add_argument("--out-prefix", default=str(ROOT / "résultats" / "phase4_p21_"))
    args = ap.parse_args()

    h37 = read_fasta(H37RV)
    anc = build_ancestral()
    masked = load_mask(MASKS, in_mtbc0_coords=True)
    (n_gc, n_at), _ = opportunity()
    cons, n_can = canettii_consensus(args.force)
    print(f"# consensus canettii : {n_can} genomes, {len(cons)} positions "
          f"variables ; seuil de quasi-fixation {args.seuil}", file=sys.stderr)

    # ---- test 1 : preuve POSITIVE sur les positions ou MTBC0 diverge de H37Rv
    div = [p for p in range(len(h37)) if anc[p] not in (ord("N"), ord(h37[p]))]
    non_masque = [p for p in div if p not in masked]
    d = Counter()
    for p in non_masque:
        hit = cons.get(p)
        if hit is None:
            d["canettii_comme_H37Rv"] += 1
        elif hit[0] == chr(anc[p]):
            d["canettii_confirme_MTBC0"] += 1
            d["canettii_confirme_MTBC0_quasifixe"] += hit[1] >= args.seuil * n_can
        else:
            d["canettii_troisieme_allele"] += 1
    tot = sum(v for k, v in d.items() if not k.endswith("quasifixe"))
    print("\n=== TEST 1 (preuve positive) : les positions ou MTBC0 dit que "
          "H37Rv porte l'etat DERIVE ===")
    print(f"  positions MTBC0 != H37Rv, hors masque : {tot} "
          f"(sur {len(div)} au total, {len(div)-tot} dans le masque)")
    print(f"  canettii porte le MEME allele que MTBC0        : "
          f"{d['canettii_confirme_MTBC0']:6d} ({100*d['canettii_confirme_MTBC0']/tot:5.2f} %)"
          f"  dont quasi-fixe : {d['canettii_confirme_MTBC0_quasifixe']}")
    print(f"  canettii porte un TROISIEME allele              : "
          f"{d['canettii_troisieme_allele']:6d} ({100*d['canettii_troisieme_allele']/tot:5.2f} %)")
    print(f"  aucune canettii ne varie (= allele H37Rv, ou non couvert) : "
          f"{d['canettii_comme_H37Rv']:6d} ({100*d['canettii_comme_H37Rv']/tot:5.2f} %)")

    # ---- test 2 : matrice de confusion + v/u par lignee sous les deux polarisations
    conf = Counter()
    rows = []
    for clade in CLADES:
        st = strain_dirs(clade)
        rng = random.Random(0)
        rng.shuffle(st)
        subsets, names = [], []
        for s in st[:args.n_per_clade]:
            v = read_subs(s / "NC_000962.3" / "spdi.txt")
            if v:
                subsets.append(v)
                names.append(s.name)
        support = defaultdict(int)
        for i, subs in enumerate(subsets):
            for v in subs:
                support[v] |= 1 << i
        c = Counter()
        for (pos, ref, alt), mask in support.items():
            if mask & (mask - 1) or pos in masked:
                continue                       # branches terminales uniquement
            m = chr(anc[pos]) if pos < len(anc) else "N"
            vm = "non_liftee" if m == "N" else verdict(m, ref, alt)
            ca = canettii_allele(pos, h37, cons, n_can, args.seuil)
            vc = verdict(ca, ref, alt)
            conf[(vm, vc)] += 1
            for tag, v, a in (("mtbc0", vm, m), ("canettii", vc, ca)):
                if v in ("correcte", "tierce"):
                    r = ref if v == "correcte" else a
                    c[f"{tag}_{flux(r, alt)}"] += 1
        rows.append(dict(clade=clade, n=len(subsets),
                         **{k: c[k] for k in ("mtbc0_loss", "mtbc0_gain",
                                              "canettii_loss", "canettii_gain")}))
    t = pd.DataFrame(rows)
    for tag in ("mtbc0", "canettii"):
        vu = [vu_ci(r[f"{tag}_loss"], r[f"{tag}_gain"], n_gc, n_at) for _, r in t.iterrows()]
        t[f"vu_{tag}"] = [v[0] for v in vu]
        t[f"vu_{tag}_lo"] = [v[1] for v in vu]
        t[f"vu_{tag}_hi"] = [v[2] for v in vu]
        t[f"gc_eq_{tag}"] = 1 / (1 + t[f"vu_{tag}"])
    t = t.sort_values("vu_mtbc0", ascending=False)
    t.to_csv(args.out_prefix + "vu_canettii.tsv", sep="\t", index=False)

    print("\n=== TEST 2 : matrice de confusion des deux verdicts d'orientation "
          "(variants terminaux des 11 pools) ===")
    labels = ["correcte", "inverse", "tierce", "ambigu", "non_liftee"]
    cm = pd.DataFrame(0, index=[l for l in labels], columns=[l for l in labels])
    for (a, b), v in conf.items():
        cm.loc[a, b] += v
    cm = cm.loc[(cm.sum(1) > 0), (cm.sum(0) > 0)]
    cm.index.name = "MTBC0 \\ canettii"
    print(cm.to_string())
    n = cm.values.sum()
    diag = sum(cm.loc[l, l] for l in cm.index if l in cm.columns)
    print(f"  accord sur {diag}/{n} = {100*diag/n:.3f} % des variants terminaux")

    print("\n=== TEST 3 : v/u par lignee sous les deux polarisations ===")
    show = t[["clade", "n", "vu_mtbc0", "vu_canettii", "gc_eq_mtbc0",
              "gc_eq_canettii"]].copy()
    show["ecart_pct"] = 100 * (show.vu_canettii - show.vu_mtbc0) / show.vu_mtbc0
    print(show.round(4).to_string(index=False))
    rho, p = stats.spearmanr(t.vu_mtbc0, t.vu_canettii)
    r, pp = stats.pearsonr(t.vu_mtbc0, t.vu_canettii)
    print(f"\n  Spearman des classements : rho = {rho:.4f} (p = {p:.2g}) ; "
          f"Pearson r = {r:.4f} (p = {pp:.2g})")
    print(f"  ecart relatif median |v/u_canettii - v/u_MTBC0| / v/u_MTBC0 : "
          f"{show.ecart_pct.abs().median():.2f} %")

    cm.to_csv(args.out_prefix + "confusion.tsv", sep="\t")
    print(f"\n# ecrit dans {args.out_prefix}[confusion|vu_canettii].tsv",
          file=sys.stderr)


if __name__ == "__main__":
    main()
