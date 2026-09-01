#!/usr/bin/env python3
"""
Objet       : P2.5 -- remplacer le comptage des substitutions PAR SOUCHE (qui
              compte n fois chaque substitution ancestrale d'un clade) par un
              comptage d'EVENEMENTS mutationnels places sur l'arbre, et
              rapporter separement le spectre des branches racine / internes /
              terminales. Sous phylogenie parfaite (clonalite MTBC), l'ensemble
              des souches portant un variant identifie sa branche : les variants
              partageant le meme "pattern" de presence sont sur la meme branche,
              et chaque variant distinct vaut UN evenement. Aucun arbre explicite
              n'est donc necessaire.
Entrees     : bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt (0-based, vs H37Rv)
Sorties     : TSV sur stdout (une ligne par clade x mode de comptage) + un TSV
              detaille par classe de substitution si --spectrum
Reutilisable: oui -- s'applique a tout clade de bdd/actuelle, tout preixe
Projet      : GC_par_lignee
Date        : 2026-08-29
"""
import argparse
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

BDD = Path("/home/christophe/docs/codes/mtbc/bdd/actuelle")
GC, AT = set("GC"), set("AT")
COMPL = {"A": "T", "T": "A", "G": "C", "C": "G"}


def strain_dirs(prefix):
    """Toutes les souches d'une lignee : conteneur nu <prefix>/ + sous-clades
    <prefix>.* (convention dir-mixte du depot, les basales vivent dans le nu)."""
    out = []
    for d in sorted(BDD.iterdir()):
        if not d.is_dir():
            continue
        if d.name == prefix or d.name.startswith(prefix + "."):
            for s in sorted(d.iterdir()):
                if s.is_dir() and (s / "NC_000962.3" / "spdi.txt").exists():
                    out.append(s)
    return out


def read_subs(spdi_path):
    """Substitutions simples d'une souche : {(pos, ref, alt)}. Les indels et les
    variants multi-nucleotidiques sont ecartes (ils ne portent pas de flux GC
    interpretable en 6 classes)."""
    subs = set()
    for line in spdi_path.read_text().splitlines():
        p = line.strip().split(":")
        if len(p) != 4:
            continue
        ref, alt = p[2].upper(), p[3].upper()
        if len(ref) != 1 or len(alt) != 1 or ref == alt:
            continue
        if ref not in COMPL or alt not in COMPL:
            continue
        subs.add((int(p[1]), ref, alt))
    return subs


def canonical(ref, alt):
    """Classe de substitution referencee pyrimidine (6 classes canoniques)."""
    if ref in "CT":
        return f"{ref}>{alt}"
    return f"{COMPL[ref]}>{COMPL[alt]}"


def flux(ref, alt):
    if ref in GC and alt in AT:
        return "loss"
    if ref in AT and alt in GC:
        return "gain"
    return "neutral"  # GC<->GC (C<->G) ou AT<->AT (A<->T) : ne change pas le GC


def summarize(counter):
    loss = sum(v for (r, a), v in counter.items() if flux(r, a) == "loss")
    gain = sum(v for (r, a), v in counter.items() if flux(r, a) == "gain")
    neut = sum(v for (r, a), v in counter.items() if flux(r, a) == "neutral")
    ratio = loss / gain if gain else float("nan")
    return loss, gain, neut, ratio


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("clades", nargs="*", default=None)
    ap.add_argument("-n", "--n-per-clade", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--spectrum", metavar="TSV", default=None)
    args = ap.parse_args()

    clades = args.clades or ["L1", "L2.2.1", "L3", "L4.1.2", "L4.3.4", "L5.1",
                             "L6.1", "L7", "L8", "L9", "L10", "Bovis.1.1",
                             "Orygis", "Canettii"]

    print("clade\tn\tmode\tn_units\tloss\tgain\tneutral\tratio_loss_gain")
    spectrum_rows = []

    for clade in clades:
        strains = strain_dirs(clade)
        if not strains:
            print(f"{clade}\tABSENT", file=sys.stderr)
            continue
        rng = random.Random(args.seed)
        rng.shuffle(strains)
        sample = strains[:args.n_per_clade]
        subsets = [read_subs(s / "NC_000962.3" / "spdi.txt") for s in sample]
        subsets = [x for x in subsets if x]
        n = len(subsets)
        if n < 2:
            print(f"{clade}\tTROP PEU ({n})", file=sys.stderr)
            continue

        # --- mode A : comptage PAR SOUCHE (design du sondage, pseudo-replique)
        per_strain = Counter()
        for subs in subsets:
            for _, r, a in subs:
                per_strain[(r, a)] += 1
        # normalise par le nombre de souches pour rester lisible
        loss, gain, neut, ratio = summarize(per_strain)
        print(f"{clade}\t{n}\tA_par_souche\t{loss+gain+neut}\t{loss}\t{gain}\t"
              f"{neut}\t{ratio:.3f}")

        # --- variant -> pattern de presence (bitmask des souches porteuses)
        support = defaultdict(int)
        for i, subs in enumerate(subsets):
            bit = 1 << i
            for v in subs:
                support[v] |= bit

        # --- mode B : comptage PAR EVENEMENT, toutes branches confondues
        by_branch = defaultdict(Counter)   # pattern -> Counter[(ref,alt)]
        for (_pos, r, a), mask in support.items():
            by_branch[mask][(r, a)] += 1

        cat_counter = {"racine": Counter(), "interne": Counter(),
                       "terminal": Counter()}
        cat_branches = Counter()
        for mask, cnt in by_branch.items():
            k = bin(mask).count("1")
            cat = "racine" if k == n else ("terminal" if k == 1 else "interne")
            cat_counter[cat].update(cnt)
            cat_branches[cat] += 1

        allev = Counter()
        for c in cat_counter.values():
            allev.update(c)
        loss, gain, neut, ratio = summarize(allev)
        print(f"{clade}\t{n}\tB_evenements_tous\t{len(support)}\t{loss}\t{gain}\t"
              f"{neut}\t{ratio:.3f}")

        # --- mode C : evenements HORS branche racine (histoire interne au clade)
        interne_plus_term = Counter()
        interne_plus_term.update(cat_counter["interne"])
        interne_plus_term.update(cat_counter["terminal"])
        loss, gain, neut, ratio = summarize(interne_plus_term)
        print(f"{clade}\t{n}\tC_evenements_hors_racine\t"
              f"{sum(interne_plus_term.values())}\t{loss}\t{gain}\t{neut}\t"
              f"{ratio:.3f}")

        # --- modes D/E/F : par categorie de branche
        for cat in ("racine", "interne", "terminal"):
            loss, gain, neut, ratio = summarize(cat_counter[cat])
            print(f"{clade}\t{n}\tD_{cat}\t{sum(cat_counter[cat].values())}\t"
                  f"{loss}\t{gain}\t{neut}\t{ratio:.3f}"
                  f"\t# {cat_branches[cat]} branches")
            if args.spectrum:
                tot = sum(cat_counter[cat].values()) or 1
                canon = Counter()
                for (r, a), v in cat_counter[cat].items():
                    canon[canonical(r, a)] += v
                for cls, v in sorted(canon.items()):
                    spectrum_rows.append(
                        f"{clade}\t{n}\t{cat}\t{cls}\t{v}\t{v/tot:.5f}")

    if args.spectrum:
        with open(args.spectrum, "w") as fh:
            fh.write("clade\tn\tbranch_class\tsubstitution\tcount\tfraction\n")
            fh.write("\n".join(spectrum_rows) + "\n")
        print(f"# spectre canonique ecrit dans {args.spectrum}", file=sys.stderr)


if __name__ == "__main__":
    main()
