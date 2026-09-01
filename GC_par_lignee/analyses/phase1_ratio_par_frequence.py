#!/usr/bin/env python3
"""
Objet       : P2.5 (diagnostic) -- le rapport pertes/gains de paires G:C n'est pas
              un scalaire par lignee, il depend de la PROFONDEUR de la branche qui
              porte l'evenement. Ce script sort le rapport en fonction de k, la
              frequence du variant dans le pool de n souches echantillonnees
              (k = 1 : branche terminale, k = n : branche racine du pool), avec
              IC95 binomial. Sert a trois choses : (a) montrer que le "gradient
              inter-lignees" du sondage n'etait que le poids relatif de la branche
              racine ; (b) separer ce qui releve du filtre de la selection de ce
              qui releve du biais d'orientation vers H37Rv ; (c) fournir le
              contraste terminal / profond que P2.4 devra confirmer apres
              polarisation par MTBC0.
Entrees     : bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt (0-based, vs H37Rv)
Sorties     : TSV par clade x classe de frequence (stdout ou --out)
Reutilisable: oui -- generique sur tout clade de bdd/actuelle
Projet      : GC_par_lignee
Date        : 2026-08-29
"""
import argparse
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import read_subs, strain_dirs, flux  # noqa: E402


def wilson(k, n, z=1.96):
    """IC95 de Wilson sur une proportion (robuste aux petits effectifs)."""
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def ratio_ci(loss, gain):
    """IC95 du rapport pertes/gains, derive de l'IC de la proportion
    p = loss / (loss + gain) par r = p / (1 - p)."""
    tot = loss + gain
    if tot == 0 or gain == 0:
        return float("nan"), float("nan")
    lo, hi = wilson(loss, tot)
    return lo / (1 - lo), hi / (1 - hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("clades", nargs="*")
    ap.add_argument("-n", "--n-per-clade", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--bins", type=int, default=8,
                    help="nombre de classes de frequence entre k=2 et k=n-1")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    clades = args.clades or ["L1", "L2.2.1", "L3", "L4.1.2", "L5.1", "L6.1",
                             "L7", "L8", "L9", "L10", "Bovis.1.1",
                             "Orygis_La3", "Caprae_La2", "Microti", "Canettii"]

    out = open(args.out, "w") if args.out else sys.stdout
    print("clade\tn\tk_class\tk_lo\tk_hi\tn_events\tloss\tgain\tratio\t"
          "ci95_lo\tci95_hi", file=out)

    for clade in clades:
        strains = strain_dirs(clade)
        if not strains:
            print(f"{clade}\tABSENT", file=sys.stderr)
            continue
        rng = random.Random(args.seed)
        rng.shuffle(strains)
        subsets = [read_subs(s / "NC_000962.3" / "spdi.txt")
                   for s in strains[:args.n_per_clade]]
        subsets = [x for x in subsets if x]
        n = len(subsets)
        if n < 4:
            print(f"{clade}\tTROP PEU ({n})", file=sys.stderr)
            continue

        support = defaultdict(int)
        for i, subs in enumerate(subsets):
            for v in subs:
                support[v] |= 1 << i

        by_k = defaultdict(Counter)
        for (_p, r, a), mask in support.items():  # noqa: F841
            by_k[bin(mask).count("1")][(r, a)] += 1

        # classes : k=1 seul, k=n seul, et le milieu decoupe en --bins classes
        edges = [(1, 1)]
        mid = list(range(2, n))
        if mid:
            step = max(1, len(mid) // args.bins)
            for i in range(0, len(mid), step):
                chunk = mid[i:i + step]
                edges.append((chunk[0], chunk[-1]))
        edges.append((n, n))

        for lo, hi in edges:
            cnt = Counter()
            for k in range(lo, hi + 1):
                cnt.update(by_k.get(k, Counter()))
            loss = sum(v for (r, a), v in cnt.items() if flux(r, a) == "loss")
            gain = sum(v for (r, a), v in cnt.items() if flux(r, a) == "gain")
            tot = sum(cnt.values())
            if tot == 0:
                continue
            r = loss / gain if gain else float("nan")
            clo, chi = ratio_ci(loss, gain)
            label = ("terminal" if lo == hi == 1 else
                     "racine" if lo == hi == n else f"k{lo}-{hi}")
            print(f"{clade}\t{n}\t{label}\t{lo}\t{hi}\t{tot}\t{loss}\t{gain}\t"
                  f"{r:.3f}\t{clo:.3f}\t{chi:.3f}", file=out)

    if args.out:
        out.close()
        print(f"# ecrit dans {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
