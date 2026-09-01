#!/usr/bin/env python3
"""
Objet       : P2.2 -- verification independante de la polarisation MTBC0
              (P2.4) et canettii (P2.1), fondee sur la PARCIMONIE plutot que
              sur un outgroup, comme le prescrit la piste P2.2 (« parcimonie
              ou ML sur l'arbre MTBC »). Principe : un variant fixe a la
              RACINE du pool d'une lignee (support = les 40 souches
              echantillonnees, classe « racine » de phase2) est, sous
              l'hypothese d'origine unique que toute la chaine de mesure
              suppose implicitement (P2.5), un evenement mutationnel PROPRE
              a cette lignee. La parcimonie la plus simple teste cette
              hypothese sans reconstruire d'arbre explicite : si le MEME
              variant (meme position, meme ref, meme alt) est egalement
              present chez une AUTRE lignee, ce n'est plus un evenement isole
              sur SA branche -- soit il est ancestral au MTBC (mal polarise
              si MTBC0 le dit derive chez H37Rv), soit c'est de l'homoplasie.

              DECOUVERTE DE PARCOURS (section A, variants BRUTS non
              polarises, calcul entierement LOCAL) : le partage inter-
              lignees est enorme, 50 a 84 %, et SYSTEMATIQUE -- une lignee
              partage ses variants de racine avec la quasi-totalite des dix
              autres A LA FOIS, jamais avec une seule (ex. L1 partage 577 de
              ses 890 variants racine avec L9, 571 avec Caprae, 534 avec L7
              ET Orygis...). C'est exactement la signature du biais de
              reference deja quantifie par P2.4 sous un angle different :
              H37Rv n'est pas l'ancetre du MTBC, donc l'ecrasante majorite
              des « variants de racine » bruts de toute lignee non-L4 sont
              des sites ou H37Rv porte l'allele DERIVE et toutes les autres
              lignees partagent l'allele ancestral -- ce n'est pas de
              l'homoplasie, c'est le meme artefact que la bascule du rapport
              racine (0,85 -> 2,87 apres polarisation, P2.4) donne a voir
              autrement.

              SECTION B, LE VRAI TEST DE LA PISTE : ne garder que les
              variants de racine pour lesquels MTBC0 CONFIRME que H37Rv porte
              l'allele ancestral (l'alt est alors un veritable etat derive
              propre a la lignee, exactement la classe comptee dans le reste
              du pipeline). Si la polarisation MTBC0 est correcte, le partage
              inter-lignees doit s'effondrer vers zero sur cette classe -- un
              signal residuel serait, lui, une vraie homoplasie ou un residu
              de biais que P2.1 (canettii) n'aurait pas vu. C'est la
              comparaison A vs B, par une methode totalement independante de
              P2.1, qui repond a la piste.
Entrees     : bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt (memes panels
              n=40 seed=0 que tout le pipeline)
              data/MTBC0/ (ancestral MTBC0 liftee sur H37Rv, P2.4)
Sorties     : resultats/phase13_p22_partage_racine.tsv (section A, brut) et
              resultats/phase13_p22_partage_racine_polarise.tsv (section B,
              apres filtre MTBC0) -- memes colonnes, comparables ligne a ligne
Reutilisable: oui -- le test « meme variant, deux lignees, support total dans
              chacune, avant/apres polarisation » est une preuve d'homoplasie
              ou de residu de biais de reference transposable a tout clade
              MTBC etudie par flux mutationnel
Projet      : GC_par_lignee
Date        : 2026-09-01
"""
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from phase3_counts_par_souche import MASKS  # noqa: E402
from phase2_polarisation_mtbc0 import load_mask, build_ancestral  # noqa: E402
from phase12_p33_couverture_singletons import build_panel  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CLADES = ["L1", "L2.2.1", "L3", "L4.1.2", "L4.3", "L6.1", "L7", "L9",
          "Orygis_La3", "Caprae_La2", "Microti"]


def root_variants(clade, masked):
    names, subsets = build_panel(clade)
    n = len(subsets)
    if n < 4:
        return None, 0
    full = (1 << n) - 1
    support = defaultdict(int)
    for i, subs in enumerate(subsets):
        for v in subs:
            support[v] |= 1 << i
    root = {v for v, mask in support.items() if mask == full and v[0] not in masked}
    return root, n


def write_sharing(per_clade, out_path):
    with open(out_path, "w") as f:
        print("clade\tn_racine\tn_partages\tpct_partages\tpartenaires", file=f)
        for clade, root in per_clade.items():
            partners = defaultdict(int)
            shared = 0
            for v in root:
                others = [c2 for c2, r2 in per_clade.items() if c2 != clade and v in r2]
                if others:
                    shared += 1
                    for c2 in others:
                        partners[c2] += 1
            pct = 100 * shared / len(root) if root else 0
            partner_str = ",".join(f"{c}:{k}" for c, k in sorted(partners.items(),
                                    key=lambda x: -x[1]))
            print(f"{clade}\t{len(root)}\t{shared}\t{pct:.3f}\t{partner_str}", file=f)
    print(f"# ecrit {out_path}", file=sys.stderr)


if __name__ == "__main__":
    masked = load_mask(MASKS, in_mtbc0_coords=True)
    anc = build_ancestral()

    raw_per_clade = {}
    polarized_per_clade = {}
    for clade in CLADES:
        root, n = root_variants(clade, masked)
        if root is None:
            print(f"{clade}\tTROP PEU", file=sys.stderr)
            continue
        raw_per_clade[clade] = root
        polarized = set()
        for (pos, ref, alt) in root:
            a = chr(anc[pos]) if pos < len(anc) else "N"
            if a == ref:
                polarized.add((pos, ref, alt))
        polarized_per_clade[clade] = polarized
        print(f"# {clade} : n={n}, racine brute={len(root)}, "
              f"racine polarisee (H37Rv=ancestral, alt=derive propre)={len(polarized)}",
              file=sys.stderr)

    write_sharing(raw_per_clade, ROOT / "résultats" / "phase13_p22_partage_racine.tsv")
    write_sharing(polarized_per_clade,
                   ROOT / "résultats" / "phase13_p22_partage_racine_polarise.tsv")
