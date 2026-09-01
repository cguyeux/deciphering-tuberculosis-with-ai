#!/usr/bin/env python3
"""P8.1.a.4 (suite) — généralise le tandem-dupliqué-complémentaire trouvé chez
Neospora caninum (F0VAI8, 2026-08-20) aux 5 autres apicomplexes à triade rompue.

Deux sources croisées, indépendantes :
  (a) data/PF04417_full.sto  — alignement Pfam FULL déjà en main (phase28/29/30),
      compte les occurrences PF04417 par accession et leur motif de triade
      (déjà calculé dans résultats/duf501_family_triad.tsv, phase28).
  (b) API InterPro protein/uniprot/<acc>/entry/pfam/PF04417 — matchs Pfam
      indépendants annotés par UniProt/InterPro (entry_protein_locations),
      qui donne les coordonnées RÉELLES du modèle HMM (pas les colonnes
      d'alignement, plus larges car elles incluent le padding autour du hit).

Sortie : résultats/phase38_apicomplexan_tandem.tsv + .txt (verdict).
"""
import csv
import json
import time
import urllib.request

APICOMPLEXES = {
    "A0A9W5T8J3": "Babesia ovis",
    "A0A1D3CV76": "Cyclospora cayetanensis",
    "F0VAI8": "Neospora caninum",
    "A0A2A9M7T6": "Besnoitia besnoiti",
    "U6GA88": "Eimeria acervulina",
    "A0A2C6L2K9": "Cystoisospora suis",
}

STO = "data/PF04417_full.sto"
TRIAD_TSV = "résultats/duf501_family_triad.tsv"


def sto_occurrences():
    """accession -> liste de (label_sto, résidu_range) présentes dans l'alignement full."""
    occ = {}
    with open(STO) as fh:
        for line in fh:
            if not line.startswith("#=GS"):
                continue
            label = line.split()[1]  # ex A0A9W5T8J3_BABOV/437-521
            acc = label.split("_")[0]
            occ.setdefault(acc, []).append(label)
    return occ


def triad_by_label():
    d = {}
    with open(TRIAD_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            d[row["sequence"]] = row
    return d


def interpro_matches(acc):
    url = f"https://www.ebi.ac.uk/interpro/api/protein/uniprot/{acc}/entry/pfam/PF04417"
    req = urllib.request.Request(url, headers={"User-Agent": "Rv1025-research/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        d = json.load(resp)
    locs = []
    for entry in d.get("entries", []):
        for loc in entry.get("entry_protein_locations", []):
            for frag in loc.get("fragments", []):
                locs.append((frag["start"], frag["end"], loc.get("score")))
    locs.sort()
    return locs


def main():
    occ = sto_occurrences()
    triad = triad_by_label()
    rows = []
    for acc, organism in APICOMPLEXES.items():
        sto_labels = sorted(occ.get(acc, []))
        n_sto = len(sto_labels)
        try:
            ip_locs = interpro_matches(acc)
        except Exception as e:
            ip_locs = []
            print(f"[WARN] InterPro échec pour {acc}: {e}")
        time.sleep(0.4)

        motifs = []
        for lab in sto_labels:
            row = triad.get(lab)
            if row:
                motifs.append(f"{row['res59']}{row['res113']}{row['res115']}")
            else:
                motifs.append("?")

        tandem_sto = n_sto >= 2
        tandem_interpro = len(ip_locs) >= 2
        gap = None
        if tandem_interpro:
            gap = ip_locs[1][0] - ip_locs[0][1]

        complementary = False
        if tandem_sto and len(motifs) == 2:
            has_E = any(m[0] == "E" and m[1] == "-" and m[2] == "-" for m in motifs)
            has_CH = any(m[0] == "-" and m[1] == "C" and m[2] == "H" for m in motifs)
            complementary = has_E and has_CH

        rows.append({
            "accession": acc,
            "organisme": organism,
            "n_occurrences_sto": n_sto,
            "n_matches_interpro": len(ip_locs),
            "interpro_coords": ";".join(f"{s}-{e}" for s, e, _ in ip_locs),
            "gap_interpro_aa": gap if gap is not None else "",
            "motifs_sto": ";".join(motifs),
            "tandem_confirme_2_sources": tandem_sto and tandem_interpro,
            "complementaire_E_puis_CH": complementary,
        })

    with open("résultats/phase38_apicomplexan_tandem.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_tandem = sum(r["tandem_confirme_2_sources"] for r in rows)
    n_compl = sum(r["complementaire_E_puis_CH"] for r in rows)
    lines = []
    lines.append(f"Apicomplexes testés : {len(rows)} (dont Neospora, cas index déjà connu)")
    lines.append(f"Tandem PF04417 confirmé par 2 sources indépendantes (sto + InterPro) : {n_tandem}/{len(rows)}")
    lines.append(f"  dont dégénérescence COMPLÉMENTAIRE (une copie E--, l'autre -CH) : {n_compl}/{len(rows)}")
    for r in rows:
        lines.append(
            f"  {r['accession']} ({r['organisme']}): sto={r['n_occurrences_sto']} "
            f"interpro={r['n_matches_interpro']} coords={r['interpro_coords']} "
            f"gap={r['gap_interpro_aa']}aa motifs={r['motifs_sto']} "
            f"tandem={r['tandem_confirme_2_sources']} complementaire={r['complementaire_E_puis_CH']}"
        )
    txt = "\n".join(lines)
    with open("résultats/phase38_apicomplexan_tandem.txt", "w") as fh:
        fh.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
