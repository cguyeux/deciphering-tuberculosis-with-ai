#!/usr/bin/env python3
"""P1.9 -- Epitopes T/B curates (IEDB) chevauchant Rv0810c.

QUESTION. Comas et al. 2010 montrent que les epitopes T du MTBC sont
paradoxalement HYPERCONSERVES plutot qu'en echappee immunitaire -- signature
d'un pathogene qui "choisit" d'etre reconnu. Rv0810c est l'un des genes les
plus contraints du proteome a sa classe de taille (P8, percentile 2,7 en
NS/S) : la contrainte pourrait recouper une pression immunitaire plutot que
(ou en plus d')une contrainte structurale pure. Jamais verifie si un epitope
T ou B curate chevauche la sequence de Rv0810c.

METHODE. Interrogation de l'API IEDB (query-api.iedb.org, PostgREST),
accession UniProt I6XWB9, volets tcell_search et bcell_search, filtre
parent_source_antigen_iri=eq.UNIPROT:I6XWB9. Les positions rendues par IEDB
se referent a l'antigene source CURE de chaque enregistrement (parfois un
orthologue ou un depot GenPept redondant) : elles ne sont PAS reutilisees
telles quelles, chaque peptide est RELOCALISE ici par recherche exacte dans
la sequence H37Rv de 60 aa.

GARDE-FOU, ecrit avant lecture des resultats (enonce de la piste) : un
"Negative" signifie teste-et-non-reactif dans UN panel d'alleles et de
donneurs, pas l'absence d'epitope pour tout HLA -- un negatif ici est peu
informatif, presque garanti par construction (protéine jamais etudiee pour
elle-meme, cf. garde-fou generique P1.1) plutot que par une propriete reelle
de la proteine. Aucun modele nul supplementaire n'est propose par la piste :
c'est justement la limite du test, a assumer explicitement plutot qu'a
sur-interpreter un eventuel zero.

Sortie : résultats/p1_9_iedb_epitopes/{epitopes_mappes.tsv,rapport.txt}
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "p1_9_iedb_epitopes"

SEQ = "MGRGRAKAKQTKVARELKYSSPQTDFQRLQRELSGTGTDRLDGDGPSDDDSWNDEDDWRR"

MODULE = (1, 33)     # tete rigide, pLDDT 91,9/85,6 (P2.1)
TAIL = (34, 60)       # queue acide desordonnee (P2.1)
PHOSPHOSITES = {20: "S20", 21: "S21", 24: "T24"}
MISSENSE_HOTSPOTS = {24: "p.Thr24Pro (32/145209 souches, clonal L4.3.3, P8.1)",
                      56: "p.Asp56Ala (73/145209 souches, clonal L4.1.1.1, P8.1)"}


def overlap(r1, r2):
    lo, hi = max(r1[0], r2[0]), min(r1[1], r2[1])
    return max(0, hi - lo + 1)


def unit_of(start, end):
    hits = []
    if overlap((start, end), MODULE):
        hits.append("module_1-33")
    if overlap((start, end), TAIL):
        hits.append("queue_34-60")
    return "+".join(hits)


def landmarks_in(start, end):
    hit = [name for pos, name in {**PHOSPHOSITES, **MISSENSE_HOTSPOTS}.items()
           if start <= pos <= end]
    return sorted(set(hit))


def main():
    assert len(SEQ) == 60, len(SEQ)
    records = []
    for kind in ("tcell", "bcell"):
        p = OUT / f"{kind}_raw.json"
        rows = json.loads(p.read_text()) if p.exists() else []
        for r in rows:
            records.append((kind, r))

    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit("P1.9 -- Epitopes T/B references (IEDB) sur Rv0810c / I6XWB9")
    emit("=" * 78)
    n_t = sum(1 for k, _ in records if k == "tcell")
    n_b = sum(1 for k, _ in records if k == "bcell")
    emit(f"Enregistrements IEDB bruts : {n_t} lymphocyte T, {n_b} lymphocyte B.")
    emit("")

    peptides = {}
    for kind, r in records:
        pep = r.get("linear_sequence")
        if not pep:
            continue
        cs = r.get("curated_source_antigen") or {}
        e = peptides.setdefault(pep, {"kind": kind, "assays": [], "iedb_pos": set()})
        e["assays"].append({
            "mesure": r.get("qualitative_measure"),
            "mhc_class": r.get("mhc_class"),
            "mhc": r.get("mhc_restriction"),
            "hote": r.get("host_organism_name"),
            "pmid": r.get("pubmed_id"),
            "antigene_source": cs.get("name"),
            "organisme_source": cs.get("source_organism_name"),
        })
        if cs.get("starting_position"):
            e["iedb_pos"].add((cs["starting_position"], cs["ending_position"]))

    rows_out = []
    emit("PEPTIDES DISTINCTS, RELOCALISES DANS LA SEQUENCE H37Rv (60 aa)")
    emit("-" * 78)
    for pep, e in sorted(peptides.items(), key=lambda kv: SEQ.find(kv[0])):
        idx = SEQ.find(pep)
        if idx < 0:
            emit(f"  {pep}  ABSENT de la sequence H37Rv -- exclu")
            continue
        start, end = idx + 1, idx + len(pep)
        pos_msg = ""
        if e["iedb_pos"]:
            declared = sorted(e["iedb_pos"])
            if all(d[0] != start for d in declared):
                pos_msg = f"   [IEDB annoncait {', '.join(f'{a}-{b}' for a, b in declared)}]"
        pos_all = [a["mesure"] for a in e["assays"] if (a["mesure"] or "").lower() == "positive"]
        neg_all = [a["mesure"] for a in e["assays"] if (a["mesure"] or "").lower() == "negative"]
        lm = landmarks_in(start, end)
        emit(f"  {pep}")
        emit(f"    H37Rv {start}-{end}  unite {unit_of(start, end)}{pos_msg}")
        if lm:
            emit(f"    chevauche : {', '.join(lm)}")
        emit(f"    dosages : {len(pos_all)} POSITIF, {len(neg_all)} negatif")
        for a in e["assays"]:
            emit(f"      - {a['mesure']:<8s} hote {a['hote']}, MHC {a['mhc_class']}/{a['mhc']},"
                 f" PMID {a['pmid']}, antigene source {a['antigene_source']!r} ({a['organisme_source']})")
        emit("")
        rows_out.append({
            "peptide": pep, "h37rv_start": start, "h37rv_end": end,
            "unite": unit_of(start, end), "chevauche": ";".join(lm) or "aucun",
            "n_positif": len(pos_all), "n_negatif": len(neg_all),
            "hotes": ";".join(sorted({a["hote"] or "?" for a in e["assays"]})),
            "pmids": ";".join(sorted({str(a["pmid"]) for a in e["assays"] if a["pmid"]})) or "aucun",
        })

    emit("=" * 78)
    emit("COUVERTURE DE LA SEQUENCE")
    covered = set()
    for r in rows_out:
        covered |= set(range(r["h37rv_start"], r["h37rv_end"] + 1))
    emit(f"  {len(covered)}/60 residus couverts par au moins un dosage IEDB "
         f"({min(covered) if covered else '-'}-{max(covered) if covered else '-'}).")
    uncovered = sorted(set(range(1, 61)) - covered)
    emit(f"  Residus JAMAIS testes : {len(uncovered)}/60.")
    emit("")

    emit("VERDICT SUR LE PARADOXE COMAS 2010 (contrainte = pression immunitaire ?)")
    pos_any = [r for r in rows_out if r["n_positif"] > 0]
    if not rows_out:
        emit("  Aucun enregistrement IEDB chevauchant Rv0810c. Absence de preuve, pas")
        emit("  preuve d'absence : le paradoxe reste totalement non instruit.")
    elif pos_any:
        emit("  >>> AU MOINS UN DOSAGE POSITIF : un confondant immunologique est REEL et")
        emit("      doit etre porte dans l'interpretation de la contrainte P8.")
        for r in pos_any:
            emit(f"      {r['peptide']}  H37Rv {r['h37rv_start']}-{r['h37rv_end']}")
    else:
        emit(f"  {len(rows_out)} peptide(s) distinct(s) testes, {len(covered)}/60 residus couverts,")
        emit("  TOUS LES DOSAGES SONT NEGATIFS en HLA classe II humaine (2 references,")
        emit("  dont Panda/Lindestam Arlehamn et al. 2024, Nat Commun, PMID 38278794).")
        emit("  >>> Le paradoxe Comas 2010 n'est PAS soutenu par la donnee curatee")
        emit("      disponible : la contrainte forte de Rv0810c (P8, percentile 2,7) n'a,")
        emit("      pour l'instant, aucun signal immunologique documente a lui attribuer.")
        emit("      Ceci COUVRE notamment le hotspot missense p.Asp56Ala (73 souches,")
        emit("      expansion clonale L4.1.1.1, P8.1), qui tombe dans le second peptide")
        emit("      teste (45-59) et n'y montre aucune reactivite.")
    emit("")
    emit("  RESERVES, a porter dans toute citation (garde-fou pose avant lecture) :")
    emit("   - 'Negative' = teste et non reactif dans UN panel d'alleles/donneurs, pas")
    emit("     l'absence d'epitope pour tout HLA -- effectif tres faible (2 peptides).")
    emit(f"   - {len(uncovered)}/60 residus n'ont jamais ete testes, dont T24 lui-meme (le")
    emit("     phosphosite nomme) et la region 1-24 en totalite : le paradoxe n'est donc")
    emit("     instruit que sur une fraction de la sequence, pas sur la proteine entiere.")
    emit("   - Deux references seulement, l'une sans PMID (depot IEDB direct) : gisement")
    emit("     etroit, coherent avec 'jamais etudiee pour elle-meme' (garde-fou P1.1).")

    if rows_out:
        hdr = list(rows_out[0].keys())
        (OUT / "epitopes_mappes.tsv").write_text(
            "\t".join(hdr) + "\n"
            + "\n".join("\t".join(str(r[h]) for h in hdr) for r in rows_out) + "\n")
    (OUT / "rapport.txt").write_text("\n".join(lines) + "\n")
    emit("")
    emit(f"Ecrit : {OUT/'epitopes_mappes.tsv'}")
    emit(f"Ecrit : {OUT/'rapport.txt'}")


if __name__ == "__main__":
    main()
