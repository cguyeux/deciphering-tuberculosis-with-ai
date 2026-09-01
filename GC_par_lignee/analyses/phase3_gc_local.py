#!/usr/bin/env python3
"""
Objet       : P3.2 (volet mecaniste) -- attaquer le biais de couverture
              GC-dependant par son mecanisme et non par ses proxys techniques.
              La couverture d'un site depend du GC du FRAGMENT qui le porte,
              pas du nucleotide lui-meme : en stratifiant les substitutions par
              le GC local (fenetre de la taille d'un fragment autour du site),
              pertes et gains de paires G:C sont compares a couverture attendue
              constante. Si le classement des lignees survit dans chaque strate
              de GC local, un biais de couverture GC-dependant ne peut pas
              l'avoir fabrique.
Entrees     : bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt
              H37Rv NC_000962.3.fasta, data/MTBC0/ancestral_on_H37Rv.bin
Sorties     : TSV clade x strate de GC local (rapport pertes/gains, IC95)
Reutilisable: oui -- la stratification par GC local vaut pour tout signal de
              composition mesure sur donnees de lecture courte
Projet      : GC_par_lignee
Date        : 2026-08-29
"""
import argparse
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import read_subs, strain_dirs, flux  # noqa: E402
from phase1_ratio_par_frequence import ratio_ci  # noqa: E402
from phase2_polarisation_mtbc0 import build_ancestral, load_mask, read_fasta, H37RV  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
MASKS = sorted((ROOT / "data" / "MTBC0" / "Mask.files").glob("*.bed"))


def gc_profile(seq, win):
    """GC en fenetre glissante centree, par position (cumsum, O(L))."""
    isgc = np.frombuffer(seq.encode(), dtype=np.uint8)
    isgc = ((isgc == ord("G")) | (isgc == ord("C"))).astype(np.int32)
    cs = np.concatenate([[0], np.cumsum(isgc)])
    h = win // 2
    n = len(seq)
    lo = np.clip(np.arange(n) - h, 0, n)
    hi = np.clip(np.arange(n) + h + 1, 0, n)
    return (cs[hi] - cs[lo]) / (hi - lo)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("clades", nargs="*")
    ap.add_argument("-n", "--n-per-clade", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--window", type=int, default=300,
                    help="fenetre de GC local, ~taille d'un fragment Illumina")
    ap.add_argument("--strates", type=int, default=4)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    anc = build_ancestral()
    masked = load_mask(MASKS, in_mtbc0_coords=True)
    gc = gc_profile(read_fasta(H37RV), args.window)
    qs = np.quantile(gc, np.linspace(0, 1, args.strates + 1))
    qs[0], qs[-1] = -1, 2
    print(f"# GC local en fenetre {args.window} pb ; bornes des strates : "
          + ", ".join(f"{q:.4f}" for q in qs[1:-1]), file=sys.stderr)

    clades = args.clades or ["L1", "L2.2.1", "L3", "L4.1.2", "L4.3", "L6.1",
                             "L7", "L9", "Orygis_La3", "Caprae_La2", "Microti"]

    out = open(args.out, "w") if args.out else sys.stdout
    print("clade\tn\tstrate_gc\tgc_lo\tgc_hi\tloss\tgain\tratio\tci95_lo\tci95_hi",
          file=out)

    for clade in clades:
        strains = strain_dirs(clade)
        if not strains:
            continue
        rng = random.Random(args.seed)
        rng.shuffle(strains)
        subsets = [read_subs(s / "NC_000962.3" / "spdi.txt")
                   for s in strains[:args.n_per_clade]]
        subsets = [x for x in subsets if x]
        n = len(subsets)
        if n < 4:
            continue

        support = defaultdict(int)
        for i, subs in enumerate(subsets):
            for v in subs:
                support[v] |= 1 << i

        per_strate = defaultdict(Counter)
        for (pos, ref, alt), mask in support.items():
            if mask & (mask - 1):      # branches terminales seulement
                continue
            if pos in masked or pos >= len(anc):
                continue
            a = chr(anc[pos])
            if a == "N" or a == alt:
                continue
            if a != ref:
                ref = a
            s = int(np.searchsorted(qs, gc[pos], side="right")) - 1
            s = min(max(s, 0), args.strates - 1)
            per_strate[s][flux(ref, alt)] += 1

        for s in range(args.strates):
            c = per_strate[s]
            l, g = c["loss"], c["gain"]
            if g == 0:
                continue
            lo, hi = ratio_ci(l, g)
            print(f"{clade}\t{n}\tQ{s+1}\t{qs[s] if s else 0:.4f}\t"
                  f"{qs[s+1] if s+1 < args.strates else 1:.4f}\t{l}\t{g}\t"
                  f"{l/g:.3f}\t{lo:.3f}\t{hi:.3f}", file=out)

    if args.out:
        out.close()
        print(f"# ecrit dans {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
