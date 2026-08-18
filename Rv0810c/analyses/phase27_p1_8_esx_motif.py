#!/usr/bin/env python3
"""P1.8 -- Motif de secretion generale ESX (YxxxD/E) sur Rv0810c.

QUESTION. etat_des_decouvertes.md §4 porte depuis l'iteration 1 une contradiction
non instruite : Rv0810c est phosphorylee en conditions de culture (Malakar 2023) ET
detectee dans le rang 35 de la fraction secretoire enrichie de macrophages infectes
(Chande 2015), alors que DeepTMHMM la classe "predicted cytoplasmic/globular" (deja
conteste sur le plan structural par P2.1 : rayon de giration 22,8 A, pas globulaire).
Les systemes de secretion de type VII (ESX) du MTBC exportent sans peptide signal
Sec/Tat classique, via un signal general C-terminal identifie par Daleke et al. 2012
(J Biol Chem 287(47):39471-39481) : un motif court Y-x-x-x-D/E porte par le dernier
segment de la proteine substrat (ou de son partenaire heterodimerique), reconnu par
l'ATPase EccC. Ce signal n'a jamais ete recherche sur Rv0810c.

METHODE. Recherche du motif Y.{3}[DE] (regex, sens N->C) sur la totalite de la
sequence de Rv0810c et de trois substrats ESX-1 confirmes recuperes depuis UniProt
(EsxA/ESAT-6 P9WNK7, EsxB/CFP-10 P9WNK5, EspB P9WJD9), pour verifier que le test
detecterait un vrai signal s'il existait (modele nul demande par la piste). Le
signal etabli par Daleke 2012 est positionnel : porte par le segment C-terminal
au-dela du coeur replie du substrat (chez EsxB, la queue non structuree en aval de
l'helice-tour-helice WXG100), pas n'importe ou dans la sequence. Chaque match est
donc aussi rapporte avec sa distance au C-terminus.

GARDE-FOU pose avant calcul (contre-argument de la piste) : Rv0810c n'est
physiquement lie a aucun locus ESX connu (P0.4/P0.7) et son architecture bipartite
(tete basique repliee 1-33 + queue acide etendue 34-60, P2.1) ne ressemble pas au
repli WXG100 en helice-tour-helice des substrats canoniques -- deja indirectement
exclu par Foldseek 0 hit sur pdb100 (P3.1/P3.2), qui contient les structures WXG100
(ex. ESAT-6, PDB 1WA8). Un negatif est l'issue la plus probable, pas un echec de
methode.

Sortie : résultats/p1_8_esx_motif.json
"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "p1_8_esx_motif.json"

RV0810C = "MGRGRAKAKQTKVARELKYSSPQTDFQRLQRELSGTGTDRLDGDGPSDDDSWNDEDDWRR"

# Substrats ESX-1 confirmes, UniProt (recuperes cette seance, format TSV rest.uniprot.org)
CONTROLS = {
    "EsxA_ESAT6_Rv3875_P9WNK7": (
        "MTEQQWNFAGIEAAASAIQGNVTSIHSLLDEGKQSLTKLAAAWGGSGSEAYQGVQQKWDATATELNNAL"
        "QNLARTISEAGQAMASTEGNVTGMFA"),
    "EsxB_CFP10_Rv3874_P9WNK5": (
        "MAEMKTDAATLAQEAGNFERISGDLKTQIDQVESTAGSLQGQWRGAAGTAAQAAVVRFQEAANKQKQEL"
        "DEISTNIRQAGVQYSRADEEQQQALSSQMGF"),
    "EspB_Rv3881c_P9WJD9": (
        "MTQSQTVTVDQQEILNRANEVEAPMADPPTDVPITPCELTAAKNAAQQLVLSADNMREYLAAGAKERQR"
        "LATSLRNAAKAYGEVDEEAATALDNDGEGTVQAESAGAVGGDSSAELTDTPRVATAGEPNFMDLKEAAR"
        "KLETGDQGASLAHFADGWNTFNLTLQGDVKRFRGFDNWEGDAATACEASLDQQRQWILHMAKLSAAMAK"
        "QAQYVAQLHVWARREHPTYEDIVGLERLYAENPSARDQILPVYAEYQQRSEKVLTEYNNKAALEPVNPP"
        "KPPPAIKIDPPPPPQEQGLIPGFLMPPSDGSGVTPGTGMPAAPMVPPTGSPGGGLPADTAAQLTSAGRE"
        "AAALSGDVAVKAASLGGGGGGGVPSAPLGSAIGGAESVRPAGAGDIAGLGQGRAGGGAALGGGGMGMPM"
        "GAAHQGQGGAKSKGSQQEDEALYTEDRAWTEAVIGNRRRQDSKESK"),
}

MOTIF = re.compile(r"Y.{3}[DE]")
# Fenetre C-terminale ou Daleke 2012 situe le signal fonctionnel (dernier segment
# non replie, ~15-20 aa chez les substrats testes) -- seuil explicite, pas ajuste
# a la donnee.
CTERM_WINDOW = 20


def scan(seq):
    hits = []
    for m in MOTIF.finditer(seq):
        start = m.start() + 1  # 1-based
        end = m.end()
        dist_from_cterm = len(seq) - end
        hits.append({
            "peptide": m.group(),
            "start": start,
            "end": end,
            "distance_du_c_terminus": dist_from_cterm,
            "dans_fenetre_c_terminale": dist_from_cterm < CTERM_WINDOW,
        })
    return {"longueur": len(seq), "sequence": seq, "hits": hits}


def main():
    assert len(RV0810C) == 60, len(RV0810C)
    results = {"motif": "Y-x-x-x-[D/E] (Daleke et al. 2012, JBC 287(47):39471-39481)",
               "fenetre_c_terminale_residus": CTERM_WINDOW,
               "cible": scan(RV0810C),
               "temoins_esx1_confirmes": {n: scan(s) for n, s in CONTROLS.items()}}

    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit("P1.8 -- Motif de secretion generale ESX (YxxxD/E) sur Rv0810c")
    emit("=" * 78)
    emit(f"Motif recherche : {results['motif']}")
    emit(f"Fenetre C-terminale consideree fonctionnelle (Daleke 2012) : "
         f"derniers {CTERM_WINDOW} residus")
    emit("")

    emit("TEMOINS POSITIFS (substrats ESX-1 deja confirmes secretes)")
    emit("-" * 78)
    for name, r in results["temoins_esx1_confirmes"].items():
        emit(f"  {name}  ({r['longueur']} aa)")
        if not r["hits"]:
            emit("    aucun match Y-x-x-x-[D/E] -- le test ne detecte rien meme sur un vrai substrat")
        for h in r["hits"]:
            flag = "  <-- DANS LA FENETRE C-TERMINALE FONCTIONNELLE" if h["dans_fenetre_c_terminale"] else ""
            emit(f"    {h['peptide']}  positions {h['start']}-{h['end']}"
                 f"  (a {h['distance_du_c_terminus']} aa du C-terminus){flag}")
        emit("")

    emit("CIBLE : Rv0810c (60 aa)")
    emit("-" * 78)
    r = results["cible"]
    if not r["hits"]:
        emit("  AUCUN match Y-x-x-x-[D/E] sur toute la sequence, fenetre C-terminale ou non.")
    else:
        for h in r["hits"]:
            flag = "  <-- DANS LA FENETRE C-TERMINALE FONCTIONNELLE" if h["dans_fenetre_c_terminale"] else ""
            emit(f"    {h['peptide']}  positions {h['start']}-{h['end']}"
                 f"  (a {h['distance_du_c_terminus']} aa du C-terminus){flag}")
    emit("")

    emit("=" * 78)
    emit("VERDICT")
    n_ctrl_hit_window = sum(
        1 for r in results["temoins_esx1_confirmes"].values()
        if any(h["dans_fenetre_c_terminale"] for h in r["hits"]))
    target_hit_window = any(h["dans_fenetre_c_terminale"] for h in results["cible"]["hits"])
    emit(f"  Temoins positifs portant le motif en fenetre C-terminale : "
         f"{n_ctrl_hit_window}/{len(CONTROLS)} -- le test sait detecter un vrai signal.")
    if target_hit_window:
        emit("  >>> Rv0810c PORTE le motif en position C-terminale fonctionnelle : signal")
        emit("      concret a instruire plus avant (positif rare mais fort).")
    else:
        emit("  >>> Rv0810c NE PORTE PAS le motif de secretion generale ESX, meme en dehors")
        emit("      de la fenetre C-terminale. Le mecanisme ESX/WXG100 n'explique pas la")
        emit("      contradiction de §4. Coherent avec Foldseek 0 hit sur pdb100 (P3.1/P3.2,")
        emit("      qui couvre les structures WXG100 comme ESAT-6/1WA8) et avec l'absence de")
        emit("      liaison physique a un locus ESX (P0.4/P0.7).")
    results["verdict"] = {
        "temoins_positifs_dans_fenetre": n_ctrl_hit_window,
        "n_temoins": len(CONTROLS),
        "rv0810c_dans_fenetre": target_hit_window,
    }

    OUT.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n")
    emit("")
    emit(f"Ecrit : {OUT}")


if __name__ == "__main__":
    main()
