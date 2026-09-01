#!/usr/bin/env python3
"""
Objet       : P2.4 -- polariser les substitutions par MTBC0, le genome ancestral
              impute du MRCA du MTBC (Harrison et al. 2024, Microb Genom), au lieu
              de les lire comme un ecart a H37Rv. Le depot MTBC0 fournit un
              liftover position par position vers H37Rv : la polarisation devient
              une jointure sur table. Trois classes en sortent pour chaque site
              substitue : orientation CORRECTE (allele MTBC0 = allele H37Rv, le
              variant est bien derive), orientation INVERSEE (allele MTBC0 =
              allele "variant", c'est H37Rv qui porte l'etat derive et la souche
              qui est ancestrale), et TIERCE (trois alleles en jeu). Le script
              recalcule le rapport pertes/gains de paires G:C par classe de
              frequence en ne comptant que les evenements correctement orientes,
              et chiffre la contamination par classe.
Entrees     : bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt (0-based, vs H37Rv)
              data/MTBC0/MTBC0_v1.1.fasta
              data/MTBC0/liftovers/H37Rv_NC_000962.3.on.MTBC0.all.positions.tsv.gz
              data/MTBC0/Mask.files/*.bed (optionnel, --mask)
Sorties     : TSV par clade x classe de frequence, avant / apres polarisation
              cache data/MTBC0/ancestral_on_H37Rv.bin (allele ancestral par
              position H37Rv 0-based, 'N' si non liftee)
Reutilisable: oui -- la table ancestrale sert a toute analyse de polarisation MTBC
Projet      : GC_par_lignee
Date        : 2026-08-29
"""
import argparse
import gzip
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import read_subs, strain_dirs, flux  # noqa: E402
from phase1_ratio_par_frequence import ratio_ci  # noqa: E402  # type: ignore

ROOT = Path(__file__).resolve().parent.parent
MTBC0_DIR = ROOT / "data" / "MTBC0"
H37RV = Path("/home/christophe/docs/codes/mtbc/Canettii/NC_000962.3.fasta")
CACHE = MTBC0_DIR / "ancestral_on_H37Rv.bin"


def read_fasta(path):
    return "".join(l.strip() for l in path.read_text().splitlines()
                   if not l.startswith(">")).upper()


def build_ancestral(verbose=True):
    """bytearray indexe par position H37Rv 0-based -> allele MTBC0 (b'N' si la
    position n'est pas liftee). Le liftover est 1-based des deux cotes."""
    if CACHE.exists():
        return bytearray(CACHE.read_bytes())
    h37 = read_fasta(H37RV)
    anc_seq = read_fasta(MTBC0_DIR / "MTBC0_v1.1.fasta")
    out = bytearray(b"N" * len(h37))
    lift = MTBC0_DIR / "liftovers" / "H37Rv_NC_000962.3.on.MTBC0.all.positions.tsv.gz"
    same = diff = 0
    with gzip.open(lift, "rt") as fh:
        next(fh)  # entete
        for line in fh:
            a, b = line.split("\t")
            pm, ph = int(a), int(b)
            if not (1 <= pm <= len(anc_seq)) or not (1 <= ph <= len(h37)):
                continue
            base = anc_seq[pm - 1]
            out[ph - 1] = ord(base)
            if base == h37[ph - 1]:
                same += 1
            else:
                diff += 1
    if verbose:
        tot = same + diff
        print(f"# liftover : {tot} positions H37Rv sur {len(h37)} "
              f"({100*tot/len(h37):.2f} %) ; MTBC0 identique a H37Rv sur "
              f"{100*same/tot:.3f} %, divergent sur {diff} positions "
              f"({100*diff/tot:.3f} %)", file=sys.stderr)
    CACHE.write_bytes(bytes(out))
    return out


def mtbc0_to_h37rv():
    """dict position MTBC0 1-based -> position H37Rv 0-based, pour lifter les
    BED de masquage du depot, qui sont en coordonnees MTBC0."""
    lift = MTBC0_DIR / "liftovers" / "H37Rv_NC_000962.3.on.MTBC0.all.positions.tsv.gz"
    m = {}
    with gzip.open(lift, "rt") as fh:
        next(fh)
        for line in fh:
            a, b = line.split("\t")
            m[int(a)] = int(b) - 1
    return m


def load_mask(paths, in_mtbc0_coords=False):
    """Positions H37Rv 0-based a exclure, depuis des BED (0-based half-open).
    Les BED du depot MTBC0 sont en coordonnees MTBC0 : on les lifte."""
    conv = mtbc0_to_h37rv() if in_mtbc0_coords else None
    masked = set()
    for p in paths:
        for line in Path(p).read_text().splitlines():
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            if len(f) < 3:
                continue
            try:
                lo, hi = int(f[1]), int(f[2])
            except ValueError:
                continue
            if conv is None:
                masked.update(range(lo, hi))
            else:
                # BED 0-based half-open sur MTBC0 -> positions MTBC0 1-based
                for pm in range(lo + 1, hi + 1):
                    ph = conv.get(pm)
                    if ph is not None:
                        masked.add(ph)
    return masked


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("clades", nargs="*")
    ap.add_argument("-n", "--n-per-clade", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--bins", type=int, default=6)
    ap.add_argument("--mask", nargs="*", default=None,
                    help="fichiers BED de regions a exclure (coord. H37Rv)")
    ap.add_argument("--mask-mtbc0", nargs="*", default=None,
                    help="idem, mais BED en coordonnees MTBC0 (depot Harrison)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    anc = build_ancestral()
    masked = load_mask(args.mask) if args.mask else set()
    if args.mask_mtbc0:
        masked |= load_mask(args.mask_mtbc0, in_mtbc0_coords=True)
    if masked:
        print(f"# masque : {len(masked)} positions H37Rv exclues", file=sys.stderr)

    clades = args.clades or ["L1", "L2.2.1", "L3", "L4.1.2", "L4.3", "L6.1",
                             "L7", "L9", "Orygis_La3", "Caprae_La2", "Microti"]

    out = open(args.out, "w") if args.out else sys.stdout
    print("clade\tn\tk_class\tk_lo\tk_hi\tn_sites\tfrac_inverse\tfrac_tierce\t"
          "frac_non_liftee\tloss_pol\tgain_pol\tratio_pol\tci95_lo\tci95_hi\t"
          "ratio_brut", file=out)

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

        # classement de chaque variant distinct selon l'allele ancestral MTBC0
        by_k_pol = defaultdict(Counter)    # evenements correctement orientes
        by_k_brut = defaultdict(Counter)   # tous, non polarises (reference)
        by_k_stat = defaultdict(Counter)   # comptes de classes d'orientation
        for (pos, ref, alt), mask in support.items():
            if pos in masked:
                continue
            k = bin(mask).count("1")
            by_k_brut[k][(ref, alt)] += 1
            a = chr(anc[pos]) if pos < len(anc) else "N"
            by_k_stat[k]["total"] += 1
            if a == "N":
                by_k_stat[k]["non_liftee"] += 1
            elif a == ref:
                by_k_stat[k]["correcte"] += 1
                by_k_pol[k][(ref, alt)] += 1
            elif a == alt:
                # H37Rv porte l'etat derive : la souche est ancestrale, il n'y a
                # pas d'evenement a compter dans ce clade
                by_k_stat[k]["inverse"] += 1
            else:
                by_k_stat[k]["tierce"] += 1
                by_k_pol[k][(a, alt)] += 1

        edges = [(1, 1)]
        mid = list(range(2, n))
        if mid:
            step = max(1, len(mid) // args.bins)
            for i in range(0, len(mid), step):
                chunk = mid[i:i + step]
                edges.append((chunk[0], chunk[-1]))
        edges.append((n, n))

        for lo, hi in edges:
            pol, brut, st = Counter(), Counter(), Counter()
            for k in range(lo, hi + 1):
                pol.update(by_k_pol.get(k, Counter()))
                brut.update(by_k_brut.get(k, Counter()))
                st.update(by_k_stat.get(k, Counter()))
            tot = st["total"]
            if tot == 0:
                continue
            lp = sum(v for (r, a), v in pol.items() if flux(r, a) == "loss")
            gp = sum(v for (r, a), v in pol.items() if flux(r, a) == "gain")
            lb = sum(v for (r, a), v in brut.items() if flux(r, a) == "loss")
            gb = sum(v for (r, a), v in brut.items() if flux(r, a) == "gain")
            rp = lp / gp if gp else float("nan")
            rb = lb / gb if gb else float("nan")
            clo, chi = ratio_ci(lp, gp)
            label = ("terminal" if lo == hi == 1 else
                     "racine" if lo == hi == n else f"k{lo}-{hi}")
            print(f"{clade}\t{n}\t{label}\t{lo}\t{hi}\t{tot}\t"
                  f"{st['inverse']/tot:.4f}\t{st['tierce']/tot:.4f}\t"
                  f"{st['non_liftee']/tot:.4f}\t{lp}\t{gp}\t{rp:.3f}\t"
                  f"{clo:.3f}\t{chi:.3f}\t{rb:.3f}", file=out)

    if args.out:
        out.close()
        print(f"# ecrit dans {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
