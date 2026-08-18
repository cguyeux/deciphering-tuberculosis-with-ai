#!/usr/bin/env python3
"""
P9 -- confirmation par une SECONDE METHODE, independante d'InterPro, des 4 genres
candidats a une perte authentique de DUF3073 soutenus par >=2 assemblages RefSeq
Complete Genome/Chromosome (phase25) : Solwaraspora, Spirillospora, Carbonicoccus,
Oryzobacter. Les 13 autres candidats (un seul assemblage chacun) restent hors de
portee de cette passe -- reportes tels quels, pas verifies ici (cf. resultat).

Pourquoi une seconde methode est necessaire (modele nul de la piste, point c) :
phase25 repose ENTIEREMENT sur l'annotation Pfam precalculee d'InterPro. Le risque
symetrique inverse a deja ete rencontre dans ce projet le jour meme (P9.1) : le
"manque" de *S. coelicolor* dans le pipeline tblastn de `annotation_mtbc` etait un
faux negatif de seuil, corrige seulement par une requete UniProt/InterPro directe.
Ici le risque est l'inverse : si InterProScan n'a jamais tourne sur le proteome
d'un genre, ou a tourne avec un seuil qui rate un orthologue divergent, ce genre
apparait absent dans phase25 sans que le gene soit reellement perdu.

Methode : BLASTp distant (NCBI qblast, meme outil que P0.1 pour la non-homologie
humaine) de la sequence proteique de Rv0810c contre `nr`, restreint par genre via
ENTREZ_QUERY, seuil large (evalue permissif) pour detecter un orthologue meme tres
divergent -- pas une recherche a seuil strict.
"""
import json
import time
from pathlib import Path

from Bio.Blast import NCBIWWW, NCBIXML

PROJ = Path(__file__).resolve().parent.parent
OUT = PROJ / "résultats" / "p9b_blast_confirmation.json"

RV0810C_SEQ = "MGRGRAKAKQTKVARELKYSSPQTDFQRLQRELSGTGTDRLDGDGPSDDDSWNDEDDWRR"

# Les 4 candidats de phase25 soutenus par >=2 assemblages independants
GENERA_TO_CONFIRM = ["Solwaraspora", "Spirillospora", "Carbonicoccus", "Oryzobacter"]


def blast_genus(genus):
    entrez_query = f"{genus}[Organism]"
    print(f"BLASTp distant vs nr, ENTREZ_QUERY={entrez_query!r} ...")
    handle = NCBIWWW.qblast(
        "blastp", "nr", RV0810C_SEQ,
        entrez_query=entrez_query,
        expect=100,  # tres permissif -- on cherche un signal faible, pas un hit fort
        hitlist_size=10,
        word_size=2,
    )
    record = NCBIXML.read(handle)
    hits = []
    for aln in record.alignments:
        for hsp in aln.hsps:
            hits.append({
                "hit_def": aln.hit_def,
                "hit_accession": aln.accession,
                "evalue": hsp.expect,
                "identity": hsp.identities,
                "align_length": hsp.align_length,
                "score": hsp.score,
            })
    return hits


def main():
    result = {}
    for genus in GENERA_TO_CONFIRM:
        try:
            hits = blast_genus(genus)
        except Exception as exc:  # pragma: no cover -- reseau/NCBI cote serveur
            hits = None
            print(f"  ECHEC pour {genus} : {exc}")
        result[genus] = {
            "n_hits": len(hits) if hits is not None else None,
            "hits": hits,
            "verdict": (
                "AUCUN HIT meme a evalue<=100 -- absence corroboree par une seconde methode"
                if hits is not None and len(hits) == 0
                else "HIT TROUVE -- le manque dans InterPro etait un faux negatif de l'outil, PAS une perte"
                if hits
                else "requete echouee, non tranche"
            ),
        }
        print(f"{genus} : {result[genus]['verdict']}")
        time.sleep(3)  # courtoisie NCBI entre deux soumissions qblast

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    print(f"\nEcrit : {OUT}")


if __name__ == "__main__":
    main()
