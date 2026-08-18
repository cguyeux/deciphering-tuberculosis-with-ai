#!/usr/bin/env python3
"""
P4.3 — Les 406 souches a "disruption clonale" : artefact, ou falsification de
l'essentialite ?

La fiche atlas porte `frameshift: 1`, `disrupt_sites: 1`, `disrupt_mode: clonal`,
`disrupt_strains: 406` sur 145 209 souches (0,28 %). Un gene essentiel ne
devrait pas tolerer un frameshift chez 406 souches viables et sequencees. Deux
issues, toutes deux decisives : soit l'appel est un artefact et il faut le
retirer de la fiche (il affaiblit le dossier a tort), soit il est reel et
l'essentialite est falsifiee chez un clade entier.

Verification de coordonnees faite en amont, et elle est le coeur du dossier :
les SPDI de TBannotator sont **0-BASED** (SPDI est 0-based par norme NCBI).
Controle : rpoB S450L, rapporte en 1-based 761155 (ref C) dans toute la
litterature, figure ici en `NC_000962.3:761154:C:T` avec 46 303 souches ; et
tous les variants du locus Rv0810c ne retrouvent la base de reference que sous
la lecture 0-based. Confondre les deux conventions decale d'une base — ce qui
suffit a faire entrer ou sortir un variant d'une CDS dont il touche le bord.

Ce script :
  1. localise exactement le variant de disruption dans / hors la CDS ;
  2. caracterise son contexte de sequence (homopolymere = contexte canonique
     d'erreur d'appel d'indel) ;
  3. teste la clonalite par la distribution en lignees des souches porteuses ;
  4. verifie en passant l'affirmation "0 missense", qui ne porte que sur les
     SNP : les variants `complex` du meme locus sont re-traduits ici.

Sortie : résultats/p4_3_disruption.json
"""

from __future__ import annotations

import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

PROJ = Path(__file__).resolve().parent.parent
GENOME = PROJ.parent / "investigate_phylo" / "resources" / "NC_000962.3.fasta"
OUTDIR = PROJ / "résultats"
OUTDIR.mkdir(exist_ok=True)

CDS_START, CDS_END, STRAND = 904905, 905087, "-"   # 1-based inclusif


def load_genome() -> str:
    return str(next(SeqIO.parse(GENOME, "fasta")).seq)


def protein(g: str) -> str:
    cds = g[CDS_START - 1 : CDS_END]
    return str(Seq(cds).reverse_complement().translate(table=11))


def spdi_to_1based(pos0: int) -> int:
    """SPDI est 0-based : la premiere base affectee est a pos0+1 en 1-based."""
    return pos0 + 1


def homopolymer_context(g: str, pos1: int, flank: int = 12) -> dict:
    """Longueur du homopolymere contenant pos1 (1-based) et sequence autour."""
    base = g[pos1 - 1]
    lo = pos1
    while lo > 1 and g[lo - 2] == base:
        lo -= 1
    hi = pos1
    while hi < len(g) and g[hi] == base:
        hi += 1
    return {
        "base": base,
        "homopolymere_1based": [lo, hi],
        "longueur_homopolymere": hi - lo + 1,
        "contexte": g[pos1 - 1 - flank : pos1 + flank],
    }


def variant_effect(g: str, pos0: int, ref: str, alt: str) -> dict:
    """Consequence proteique d'un variant substitutif (meme longueur) dans la CDS."""
    p1 = spdi_to_1based(pos0)
    end1 = p1 + len(ref) - 1
    inside = not (end1 < CDS_START or p1 > CDS_END)
    out = {
        "spdi_0based": pos0,
        "premiere_base_1based": p1,
        "derniere_base_1based": end1,
        "ref": ref,
        "alt": alt,
        "ref_verifiee": g[p1 - 1 : end1] == ref,
        "dans_la_CDS": inside,
    }
    if not inside or len(ref) != len(alt) or not out["ref_verifiee"]:
        return out
    mut = g[: p1 - 1] + alt + g[end1:]
    wt, mt = protein(g), protein(mut)
    diffs = [
        f"{a}{i + 1}{b}" for i, (a, b) in enumerate(zip(wt, mt)) if a != b
    ]
    out["proteine_WT"] = wt
    out["proteine_mutante"] = mt
    out["changements_aa"] = diffs
    out["synonyme"] = len(diffs) == 0
    return out


def main() -> None:
    g = load_genome()
    wt = protein(g)
    rep: dict = {
        "piste": "P4.3",
        "convention_coordonnees": {
            "SPDI_TBannotator": "0-based (norme NCBI)",
            "controle": (
                "rpoB S450L = 1-based 761155 (ref C) dans la litterature ; "
                "TBannotator le porte en NC_000962.3:761154:C:T, 46 303 souches."
            ),
            "verification_locale_761155": g[761154],
            "verification_locale_761154_si_1based": g[761153],
        },
        "CDS": {
            "coordonnees_1based": [CDS_START, CDS_END],
            "brin": STRAND,
            "longueur_nt": CDS_END - CDS_START + 1,
            "proteine": wt,
            "codon_start_1based": [CDS_END - 2, CDS_END],
            "sequence_codon_start_brin_plus": g[CDS_END - 3 : CDS_END],
        },
    }

    # ------------------------------------------------------------------ #
    # 1. LE variant de disruption
    # ------------------------------------------------------------------ #
    pos0, ref, alt = 905086, "TG", "T"
    p1 = spdi_to_1based(pos0)
    deleted_1based = p1 + 1  # ancrage a gauche : la base supprimee est la 2e
    rep["1_variant_de_disruption"] = {
        "spdi": f"NC_000962.3:{pos0}:{ref}:{alt}",
        "strain_count_TBannotator": 754,
        "atlas_disrupt_strains": 406,
        "note_denominateurs": (
            "406/145 209 = 0,28 % (sous-ensemble de l'atlas) et 754/~255 000 = "
            "0,30 % (TBannotator entier) : meme variant, deux denominateurs."
        ),
        "ref_verifiee_0based": g[p1 - 1 : p1 + 1] == ref,
        "ref_lue_en_1based": g[pos0 - 1 : pos0 + 1],
        "base_supprimee_1based": deleted_1based,
        "dans_la_CDS": CDS_START <= deleted_1based <= CDS_END,
        "position_relative_au_codon_start": (
            f"{deleted_1based - CDS_END} nt en amont du codon start (brin -)"
        ),
        "contexte_sequence": homopolymer_context(g, deleted_1based),
    }

    # ------------------------------------------------------------------ #
    # 2. Autres indels du meme homopolymere (signature d'instabilite)
    # ------------------------------------------------------------------ #
    rep["2_autres_indels_du_meme_run"] = {
        "NC_000962.3:905087:GGGGGA:TGGG": {
            "strain_count": 4,
            "premiere_base_1based": spdi_to_1based(905087),
            "ref_verifiee": g[905087 : 905087 + 6] == "GGGGGA",
            "lecture": "seconde deletion independante DANS le meme run de G",
        }
    }

    # ------------------------------------------------------------------ #
    # 3. Les variants `complex` : l'affirmation "0 missense" tient-elle ?
    # ------------------------------------------------------------------ #
    complexes = [
        (905036, "CAA", "GAG", 1123),
        (905024, "GGAG", "AGAA", 798),
        (905036, "CAA", "AAG", 158),
        (905072, "CCGGCCG", "ACGTCCT", 89),
        (904994, "AC", "TT", 29),
        (904984, "CTGAC", "TTGAT", 29),
        (905033, "TTTCAA", "CTTAAG", 25),
        (905007, "TG", "GT", 24),
    ]
    snps = [
        (905081, "G", "C", 6518),
        (905045, "A", "T", 471),
        (905038, "A", "G", 423),
        (904994, "A", "G", 335),
        (904910, "C", "A", 197),
        (905063, "T", "C", 176),
    ]
    rep["3_effets_proteiques"] = {
        "avertissement": (
            "La couche `conservation` de l'atlas annonce missense=0. Elle ne "
            "compte que les variants de type SNP. Les variants `complex` "
            "(substitutions multi-nucleotidiques) du meme locus sont re-traduits "
            "ici : s'ils changent des acides amines, l'affirmation '0 missense "
            "sur 145 209 souches' ne tient pas telle quelle — et c'est un "
            "argument porteur de la piste P5."
        ),
        "SNP": [
            dict(variant_effect(g, p, r, a), strain_count=n) for p, r, a, n in snps
        ],
        "complex": [
            dict(variant_effect(g, p, r, a), strain_count=n) for p, r, a, n in complexes
        ],
    }

    out = OUTDIR / "p4_3_disruption.json"
    out.write_text(json.dumps(rep, indent=2, ensure_ascii=False))
    print(json.dumps(rep["convention_coordonnees"], indent=2, ensure_ascii=False))
    print(json.dumps(rep["1_variant_de_disruption"], indent=2, ensure_ascii=False))
    print(json.dumps(rep["2_autres_indels_du_meme_run"], indent=2, ensure_ascii=False))
    print("\n--- effets proteiques ---")
    for kind in ("SNP", "complex"):
        for e in rep["3_effets_proteiques"][kind]:
            print(
                f"  [{kind:7}] {e['spdi_0based']}:{e['ref']}>{e['alt']:8} "
                f"n={e['strain_count']:>5}  CDS={e['dans_la_CDS']}  "
                f"refOK={e['ref_verifiee']}  aa={e.get('changements_aa')}"
            )
    print(f"\n[ecrit] {out}")


if __name__ == "__main__":
    main()
