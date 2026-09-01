#!/usr/bin/env python3
"""
Objet       : P4.1 -- partitionner le genome H37Rv en CLASSES DE DEGENERESCENCE
              a partir du GFF3 des CDS et de la sequence ANCESTRALE (MTBC0 lifte
              sur H37Rv), puis mesurer l'opportunite de chaque classe : nombre de
              sites G:C et A:T disponibles, composition, opportunite
              trinucleotidique. C'est le socle de P4.2, qui comparera le flux
              mutationnel v/u entre sites quasi-neutres et sites sous selection.

              POURQUOI CETTE PARTITION, ET POURQUOI ELLE EST LE SOCLE. Tout le
              projet mesure v/u, le rapport des taux par site disponible
              (A15), sur le genome non masque pris en bloc. Ce bloc melange des
              sites ou une substitution G:C -> A:T est invisible a la selection
              (troisieme position d'un codon 4-fois degenere) et des sites ou
              elle change systematiquement l'acide amine (site non degenere).
              Si le classement inter-lignees de A15 etait porte par la
              selection plutot que par la mutation, il devrait s'affaisser sur
              les seuls sites 4-fois degeneres. A9 a deja montre que le rapport
              est PLAT le long de la profondeur de l'arbre, ce qui affaiblit
              a priori l'hypothese selective ; P4 la teste par le chemin
              orthogonal, la position dans le codon plutot que le temps.

              CE QUE MESURE LA DEGENERESCENCE, ET CE QU'ELLE NE MESURE PAS. Un
              site n-fois degenere est un site ou n des 4 bases possibles
              codent le meme acide amine. 4-fois degenere = toute substitution
              y est synonyme, donc le site est quasi-neutre au premier ordre.
              Non degenere (dit « 0-fold ») = les trois substitutions changent
              l'acide amine. La classe ne dit rien de la selection sur l'usage
              du codon ni de la structure des ARNm, qui reste la reserve
              standard de la methode et sera ecrite telle quelle.

              TROIS PIEGES FERMES ICI, PAS EN P4.2.
              (a) COMPOSITION. Les sites 4-fois degeneres sont massivement G:C
                  chez une bacterie a 65 % de GC (troisieme position de codon).
                  Un rapport pertes/gains BRUT y sera donc mecaniquement plus
                  eleve qu'ailleurs, sans qu'aucune biologie ne change : c'est
                  exactement le piege que A14 a documente sur le GC local. La
                  seule quantite comparable entre classes est le taux par site
                  DISPONIBLE, d'ou le comptage separe des sites G:C et A:T de
                  chaque classe fait ici.
              (b) CONTEXTE. A29 a etabli que le contexte trinucleotidique porte
                  de l'heterogeneite inter-lignees, surtout dans les GAINS. Les
                  classes de degenerescence n'ont pas le meme profil
                  trinucleotidique (une troisieme position de codon est
                  flanquee des deux premieres). L'opportunite trinucleotidique
                  par classe est donc calculee ici, pour que P4.2 puisse tester
                  si un ecart entre classes n'est pas un echo de A29.
              (c) ARN STRUCTURAUX. tRNA, rRNA et ncRNA ne sont pas codants mais
                  ne sont pas neutres non plus : les verser dans l'intergenique
                  y injecterait de la selection forte. Ils sont ici une classe
                  a part, exclue de l'intergenique et rapportee separement.

              REFERENCE DE SEQUENCE. La degenerescence est lue sur la sequence
              ANCESTRALE (MTBC0 lifte), pas sur H37Rv, pour la meme raison que
              la polarisation : H37Rv est une souche L4 et non l'ancetre. La
              divergence entre les deux etant de l'ordre du millier de
              positions, la section D chiffre le desaccord de classement plutot
              que de le supposer negligeable.

Entrees     : investigate_phylo/resources/NC_000962.3.gff3  (CDS + ARN H37Rv)
              Canettii/NC_000962.3.fasta                    (H37Rv, sensibilite)
              data/MTBC0/ancestral_on_H37Rv.bin             (ancestral, P2.4)
              data/mask_h37rv_positions.npy                 (masque, P3.1)
Sorties     : data/degeneracy_h37rv.npy            (classe par position, int8)
              résultats/phase9_p41_classes_sites.tsv        (A)
              résultats/phase9_p41_opportunite_trinuc.tsv   (B)
              résultats/phase9_p41_sensibilite_reference.tsv (D)
Reutilisable: oui -- la classification de degenerescence vectorisee par CDS et
              le comptage d'opportunite par classe valent pour tout genome
              bacterien annote, independamment du MTBC
Projet      : GC_par_lignee
Date        : 2026-08-31
"""
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from phase2_polarisation_mtbc0 import build_ancestral, read_fasta, H37RV  # noqa: E402
from phase4_p93_force_maintien import CODE, GFF3  # noqa: E402
from phase5_p51_panel_recursif import masked_positions  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
DEG_CACHE = ROOT / "data" / "degeneracy_h37rv.npy"
BASES = "ACGT"

# Classes de sites. L'ordre est celui du gradient de contrainte attendu :
# non degenere (toute substitution est non synonyme) -> 4-fois degenere (aucune
# ne l'est). NC = non classe (base inconnue, codon indecidable, chevauchement).
NC, INTER, ARN, F0, F2, F3, F4 = -1, 0, 1, 2, 3, 4, 5
# F0P2 n'est pas une classe du partitionnement mais un SOUS-ENSEMBLE de F0 : les
# sites non degeneres situes en DEUXIEME position de codon. Il existe pour fermer
# le seul biais d'ascertainment de la partition. Qu'un site soit 4-fois degenere
# ne depend que des deux autres bases de son codon, jamais de la sienne ; mais
# qu'un site de PREMIERE position soit 0-fold ou 2-fold depend, lui, de la base
# qu'il porte (un C en tete d'un codon CTN de leucine est 2-fold, un G en tete
# d'un GTN de valine est 0-fold). La composition de la classe 0-fold est donc en
# partie fabriquee par la regle de classement, ce qui contaminerait toute
# comparaison de NIVEAU avec le 4-fold. La deuxieme position de codon est
# non degeneree quelle que soit la base : le sous-ensemble F0P2 est le meme
# contraste, purge de ce biais.
F0P2 = 6
NOMS = {NC: "non_classe", INTER: "intergenique", ARN: "arn_structural",
        F0: "0fold", F2: "2fold", F3: "3fold", F4: "4fold",
        F0P2: "0fold_pos2"}
# Classes retenues pour la comparaison de P4.2 : les deux extremes de la
# contrainte codante, plus l'intergenique comme troisieme regime.
CLASSES_P42 = [F4, F0, INTER]
CPOS_CACHE = ROOT / "data" / "codonpos_h37rv.npy"


# revcomp local : celui de phase4_p93 leve une KeyError sur 'N', et l'ancestral
# MTBC0 en porte partout ou le liftover ne dit rien. Ici un N doit produire un
# codon indecidable, pas une exception.
_COMPL = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(_COMPL)[::-1]


def table_degenerescence():
    """fold[codon][i] = nombre de bases donnant le MEME acide amine que le
    codon a la position i (1 = non degenere, 4 = 4-fois degenere). Le codon est
    lu sur le brin CODANT ; les codons stop sont traites comme un acide amine
    a part entiere, donc un changement vers stop est non synonyme."""
    fold = {}
    for cod, aa in CODE.items():
        f = []
        for i in range(3):
            n = sum(1 for b in BASES
                    if CODE.get(cod[:i] + b + cod[i + 1:]) == aa)
            f.append(n)
        fold[cod] = tuple(f)
    return fold


def features(kinds):
    """Intervalles 0-based semi-ouverts (start, end, brin) des features du GFF3
    dont le type est dans `kinds`."""
    out = []
    for line in GFF3.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 8 or f[2] not in kinds:
            continue
        out.append((int(f[3]) - 1, int(f[4]), f[6]))
    return out


def classify(seq, verbose=True):
    """Classe de degenerescence par position 0-based, sur la sequence `seq`
    (bytes-like ou str indexable). Les CDS sont traites dans l'ordre du GFF3 et
    le premier arrive garde la position, meme convention que `load_cds` de
    phase4_p93, pour que P4.2 et le pN/pS de P6 parlent des memes sites."""
    n = len(seq)
    cls = np.full(n, INTER, np.int8)
    cpos = np.full(n, -1, np.int8)          # position dans le codon (0, 1, 2)
    fold = table_degenerescence()
    pris = np.zeros(n, bool)

    for s, e, st in features({"tRNA", "rRNA", "ncRNA"}):
        cls[s:e] = ARN
        pris[s:e] = True

    n_cds = n_chev = n_indec = 0
    for s, e, st in features({"CDS"}):
        if (e - s) % 3:
            continue
        n_cds += 1
        brut = seq[s:e] if isinstance(seq, str) else \
            bytes(seq[s:e]).decode("ascii")
        cod_seq = brut if st == "+" else revcomp(brut)
        for j in range(0, len(cod_seq) - 2, 3):
            cod = cod_seq[j:j + 3]
            f = fold.get(cod)
            for i in range(3):
                c = j + i                        # index dans la CDS codante
                p = s + c if st == "+" else e - 1 - c
                if pris[p]:
                    n_chev += 1
                    continue
                pris[p] = True
                cpos[p] = i
                if f is None:
                    cls[p] = NC
                    n_indec += 1
                else:
                    cls[p] = {1: F0, 2: F2, 3: F3, 4: F4}[f[i]]
    # Une base inconnue (N ancestral non lifte) n'est jamais un site disponible.
    arr = np.frombuffer(seq.encode() if isinstance(seq, str) else bytes(seq),
                        dtype=np.uint8)
    inconnu = ~np.isin(arr, np.frombuffer(BASES.encode(), dtype=np.uint8))
    cls[inconnu] = NC
    cpos[inconnu] = -1
    if verbose:
        print(f"# {n_cds} CDS en phase ; {n_chev} positions en chevauchement "
              f"laissees au premier CDS ; {n_indec} positions de codon "
              f"indecidable ; {int(inconnu.sum()):,} bases inconnues",
              file=sys.stderr)
    return cls, cpos


def opportunite(cls, arr, masked_mask, cpos=None):
    """Par classe : nombre de sites G:C et A:T disponibles (non masques, base
    connue). C'est le denominateur de v et de u -- la seule facon de comparer
    des classes dont la composition differe d'un facteur trois."""
    gc = np.isin(arr, np.frombuffer(b"GC", dtype=np.uint8))
    at = np.isin(arr, np.frombuffer(b"AT", dtype=np.uint8))
    ok = (~masked_mask) & (gc | at)
    rows = []
    for c, nom in NOMS.items():
        sel = ok & ((cls == F0) & (cpos == 1) if c == F0P2 else cls == c)
        if c == F0P2 and cpos is None:
            continue
        n_gc, n_at = int((sel & gc).sum()), int((sel & at).sum())
        tot = n_gc + n_at
        rows.append(dict(classe=nom, code=c, sites=tot, sites_GC=n_gc,
                         sites_AT=n_at,
                         gc_pct=100 * n_gc / tot if tot else np.nan,
                         part_genome_pct=100 * tot / int(ok.sum())))
    return pd.DataFrame(rows).sort_values("code")


def opportunite_trinuc(cls, arr, masked_mask):
    """Opportunite trinucleotidique par classe, centre replie sur la pyrimidine
    (meme convention qu'en P8.1, pour que les deux tables se croisent). Le
    CENTRE doit etre non masque et les trois bases connues ; les flancs peuvent
    appartenir a une autre classe, ce qui est la realite d'une troisieme
    position de codon flanquee des deux premieres."""
    n = len(arr)
    code = np.full(n, -1, np.int8)
    for i, b in enumerate(BASES):
        code[arr == ord(b)] = i
    keep = np.zeros(n, bool)
    keep[1:n - 1] = (code[:-2] >= 0) & (code[1:-1] >= 0) & (code[2:] >= 0)
    keep &= ~masked_mask
    rows = []
    for c in CLASSES_P42:
        idx = np.nonzero(keep & (cls == c))[0]
        l5, ce, r3 = (code[idx - 1].astype(int), code[idx].astype(int),
                      code[idx + 1].astype(int))
        pur = (ce == 0) | (ce == 2)                    # A ou G : on retourne
        L = np.where(pur, 3 - r3, l5)
        C = np.where(pur, 3 - ce, ce)
        R = np.where(pur, 3 - l5, r3)
        cnt = np.bincount(C * 16 + L * 4 + R, minlength=64)
        for ci in (1, 3):
            for a in range(4):
                for b in range(4):
                    rows.append(dict(
                        classe=NOMS[c],
                        trinuc=f"{BASES[ci]}@{BASES[a]}.{BASES[b]}",
                        sites=int(cnt[ci * 16 + a * 4 + b]),
                        flancs_gc=(BASES[a] in "GC") + (BASES[b] in "GC")))
    return pd.DataFrame(rows)


def main():
    anc = build_ancestral()
    masked = masked_positions()
    n = len(anc)
    masked_mask = np.zeros(n, bool)
    m = np.fromiter(masked, np.int64, len(masked))
    masked_mask[m[(m >= 0) & (m < n)]] = True

    cls, cpos = classify(anc)
    np.save(DEG_CACHE, cls)
    np.save(CPOS_CACHE, cpos)
    arr = np.frombuffer(bytes(anc), dtype=np.uint8)

    # ---- A. opportunite par classe
    opp = opportunite(cls, arr, masked_mask, cpos)
    opp.to_csv(ROOT / "résultats" / "phase9_p41_classes_sites.tsv",
               sep="\t", index=False)
    print("=== A. classes de degenerescence : sites disponibles (non masques, "
          "base ancestrale connue) ===")
    print(opp.to_string(index=False,
                        float_format=lambda x: f"{x:8.3f}"))
    ref = opp.set_index("classe")
    g4, g0 = ref.loc["4fold", "gc_pct"], ref.loc["0fold", "gc_pct"]
    print(f"\n  GC des sites 4-fois degeneres : {g4:.2f} % contre "
          f"{g0:.2f} % aux sites non degeneres et "
          f"{ref.loc['intergenique','gc_pct']:.2f} % en intergenique.")
    print(f"  Le rapport pertes/gains BRUT est donc attendu "
          f"{(g4/(100-g4))/(g0/(100-g0)):.2f} fois plus grand en 4-fold qu'en "
          "0-fold par pure composition : c'est le piege (a), et c'est pourquoi "
          "P4.2 ne comparera que des taux par site disponible.")
    n2 = ref.loc["0fold_pos2"]
    print(f"  Sous-ensemble de controle : les {int(n2.sites):,} sites 0-fold de "
          f"DEUXIEME position de codon (GC = {n2.gc_pct:.2f} %), non degeneres "
          "quelle que soit la base qu'ils portent, donc exempts du biais "
          "d'ascertainment que la regle de classement introduit en premiere "
          "position.")

    # ---- B. opportunite trinucleotidique par classe
    tri = opportunite_trinuc(cls, arr, masked_mask)
    tri.to_csv(ROOT / "résultats" / "phase9_p41_opportunite_trinuc.tsv",
               sep="\t", index=False)
    print("\n=== B. opportunite trinucleotidique par classe (centre pyrimidine) ===")
    piv = tri.pivot_table(index="trinuc", columns="classe", values="sites")
    part = piv / piv.sum(0)
    d = (part["4fold"] - part["0fold"]).abs().sum() / 2
    print(f"  distance en variation totale entre le profil trinucleotidique du "
          f"4-fold et celui du 0-fold : {d:.3f}")
    print("  contextes les plus discordants (part 4fold - part 0fold) :")
    ecart = (part["4fold"] - part["0fold"]).sort_values()
    for k in list(ecart.index[:3]) + list(ecart.index[-3:]):
        print(f"    {k:<8} {100*part.loc[k,'4fold']:6.2f} % vs "
              f"{100*part.loc[k,'0fold']:6.2f} %  ({100*ecart[k]:+6.2f} pt)")
    print("  Ce profil n'est PAS le meme d'une classe a l'autre : P4.2 devra "
          "montrer qu'un ecart de v/u entre classes n'est pas un echo de A29.")

    # ---- D. sensibilite a la sequence de reference
    h37 = read_fasta(H37RV)
    cls_h, _ = classify(h37, verbose=False)
    comparable = (cls != NC) & (cls_h != NC)
    desaccord = comparable & (cls != cls_h)
    lignes = [dict(positions_comparables=int(comparable.sum()),
                   desaccords=int(desaccord.sum()),
                   pct=100 * int(desaccord.sum()) / int(comparable.sum()))]
    for c in CLASSES_P42:
        sel = comparable & (cls == c)
        lignes.append(dict(classe=NOMS[c], n=int(sel.sum()),
                           change_de_classe=int((sel & desaccord).sum())))
    pd.DataFrame(lignes).to_csv(
        ROOT / "résultats" / "phase9_p41_sensibilite_reference.tsv",
        sep="\t", index=False)
    print(f"\n=== D. sensibilite : ancestral MTBC0 contre H37Rv ===")
    print(f"  {int(desaccord.sum()):,} positions changent de classe sur "
          f"{int(comparable.sum()):,} comparables "
          f"({100*int(desaccord.sum())/int(comparable.sum()):.4f} %). "
          "La classification ne depend pas du choix de reference.")


if __name__ == "__main__":
    main()
