#!/usr/bin/env python3
"""
Objet       : P3.2 -- produire, POUR CHAQUE SOUCHE, le compte de substitutions
              terminales polarisees (branche propre a la souche : variants
              qu'elle seule porte dans le pool de sa lignee), reparties en
              pertes de paires G:C, gains, et neutres, plus le spectre a 6
              classes canoniques. C'est l'unite d'observation qui permet la
              stratification par BioProject / plateforme / profondeur : le
              rapport pertes/gains cesse d'etre un scalaire par lignee et
              devient une variable mesuree sur chaque souche, joignable aux
              metadonnees de sequencage.
Entrees     : bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt
              data/MTBC0/ancestral_on_H37Rv.bin (via phase2, construit si absent)
              data/MTBC0/Mask.files/*.bed (coordonnees MTBC0)
Sorties     : TSV une ligne par souche : clade, sra, n_term, loss, gain, neutral,
              inverse, tierce, plus les 6 classes canoniques
Reutilisable: oui -- meme geste pour toute analyse de spectre par souche du depot
Projet      : GC_par_lignee
Date        : 2026-08-29
"""
import argparse
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import read_subs, strain_dirs, flux, canonical  # noqa: E402
from phase2_polarisation_mtbc0 import build_ancestral, load_mask  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
MASKS = sorted((ROOT / "data" / "MTBC0" / "Mask.files").glob("*.bed"))
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("clades", nargs="*")
    ap.add_argument("-n", "--n-per-clade", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--no-mask", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    anc = build_ancestral()
    masked = set() if args.no_mask else load_mask(MASKS, in_mtbc0_coords=True)
    print(f"# masque : {len(masked)} positions H37Rv exclues", file=sys.stderr)

    clades = args.clades or ["L1", "L2.2.1", "L3", "L4.1.2", "L4.3", "L6.1",
                             "L7", "L9", "Orygis_La3", "Caprae_La2", "Microti"]

    out = open(args.out, "w") if args.out else sys.stdout
    print("clade\tsra\tn_pool\tn_term_raw\tloss\tgain\tneutral\tinverse\ttierce"
          "\t" + "\t".join(CLASSES), file=out)

    for clade in clades:
        strains = strain_dirs(clade)
        if not strains:
            print(f"{clade}\tABSENT", file=sys.stderr)
            continue
        rng = random.Random(args.seed)
        rng.shuffle(strains)
        sample = strains[:args.n_per_clade]
        subsets = []
        names = []
        for s in sample:
            v = read_subs(s / "NC_000962.3" / "spdi.txt")
            if v:
                subsets.append(v)
                names.append(s.name)
        n = len(subsets)
        if n < 4:
            print(f"{clade}\tTROP PEU ({n})", file=sys.stderr)
            continue

        support = defaultdict(int)
        for i, subs in enumerate(subsets):
            for v in subs:
                support[v] |= 1 << i

        per_strain = defaultdict(Counter)
        for (pos, ref, alt), mask in support.items():
            if mask & (mask - 1):      # k >= 2, pas un singleton
                continue
            if pos in masked:
                continue
            i = mask.bit_length() - 1
            c = per_strain[i]
            c["n_term_raw"] += 1
            a = chr(anc[pos]) if pos < len(anc) else "N"
            if a == "N":
                continue
            if a == alt:
                c["inverse"] += 1
                continue
            if a != ref:
                c["tierce"] += 1
                ref = a                # allele ancestral tierce : on polarise dessus
            c[flux(ref, alt)] += 1
            c[canonical(ref, alt)] += 1

        for i, name in enumerate(names):
            c = per_strain[i]
            print(f"{clade}\t{name}\t{n}\t{c['n_term_raw']}\t{c['loss']}\t"
                  f"{c['gain']}\t{c['neutral']}\t{c['inverse']}\t{c['tierce']}\t"
                  + "\t".join(str(c[k]) for k in CLASSES), file=out)

    if args.out:
        out.close()
        print(f"# ecrit dans {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
