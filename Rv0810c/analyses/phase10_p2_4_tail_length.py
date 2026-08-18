#!/usr/bin/env python3
"""
P2.4 — longueur de la queue acide C-terminale de DUF3073 : mesure sur les 1930
sequences de l'alignement Pfam PF11273_full, et correlation avec la taxonomie
et un proxy de taille de genome recupere via UniProt/NCBI.
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict

STO = "/home/christophe/docs/codes/mtbc/Rv0810c/data/PF11273_full.sto"
OUT_JSON = "/home/christophe/docs/codes/mtbc/Rv0810c/résultats/p2_4_tail_length.json"
OUT_TSV = "/home/christophe/docs/codes/mtbc/Rv0810c/résultats/p2_4_tail_length_par_sequence.tsv"
CACHE_JSON = "/home/christophe/docs/codes/mtbc/Rv0810c/résultats/p2_4_uniprot_cache.json"

REF_ID = "I6XWB9_MYCTU/2-60"
MODULE_END_RESIDUE = 33  # frontiere ordre/desordre etablie en P2.1 (pLDDT 85.6 -> 62.0)


def parse_stockholm(path):
    seqs = {}       # id -> aligned row (196 col)
    ac_map = {}     # id -> accession (sans version)
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#=GS") and " AC " in line:
                parts = line.split()
                # #=GS <id> AC <accession.version>
                sid = parts[1]
                acc = parts[3]
                ac_map[sid] = acc.split(".")[0]
            elif line.startswith("#") or line.startswith("//") or not line.strip():
                continue
            else:
                parts = line.split(None, 1)
                if len(parts) == 2:
                    seqs[parts[0]] = parts[1]
    return seqs, ac_map


def find_module_boundary_column(ref_row, module_end_residue, ref_start_residue):
    """Colonne d'alignement (0-based) correspondant au dernier residu du module.

    ref_start_residue = residu (numerotation pleine longueur, Met=1) ou commence
    la ligne de reference dans l'alignement (I6XWB9_MYCTU/2-60 -> demarre au residu 2,
    le Met initial n'est jamais inclus dans cette famille Pfam). Sans ce decalage,
    compter 33 caracteres non-gap depuis le debut de la ligne atterrit sur le residu
    34 en numerotation reelle, pas 33 (erreur d'un residu, verifiee sur H37Rv :
    total_len mesure = 59 au lieu de 60 avant correction).
    """
    target_count = module_end_residue - ref_start_residue + 1
    residue_count = 0
    for col, ch in enumerate(ref_row):
        if ch not in (".", "-"):
            residue_count += 1
            if residue_count == target_count:
                return col
    raise ValueError("reference sequence too short")


def measure(seqs, boundary_col):
    out = {}
    for sid, row in seqs.items():
        module_part = row[: boundary_col + 1]
        tail_part = row[boundary_col + 1 :]
        module_len = sum(1 for c in module_part if c not in (".", "-"))
        tail_len = sum(1 for c in tail_part if c not in (".", "-"))
        out[sid] = {
            "module_len": module_len,
            "tail_len": tail_len,
            "total_len": module_len + tail_len,
        }
    return out


def organism_code(sid):
    m = re.match(r"^[A-Za-z0-9]+_([A-Za-z0-9]+)/", sid + "/")
    return m.group(1) if m else None


def batched(iterable, n):
    it = list(iterable)
    for i in range(0, len(it), n):
        yield it[i : i + n]


def fetch_uniprot_metadata(accessions, batch_size=90, sleep=0.34):
    """Retourne acc -> {organism_name, lineage(list), length} via l'API REST UniProt."""
    meta = {}
    fields = "accession,organism_name,lineage,length"
    for batch in batched(sorted(set(accessions)), batch_size):
        query = " OR ".join(f"accession:{a}" for a in batch)
        url = (
            "https://rest.uniprot.org/uniprotkb/search?"
            + urllib.parse.urlencode({"query": query, "fields": fields, "format": "tsv", "size": 500})
        )
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                text = resp.read().decode("utf-8")
        except Exception as e:
            print(f"  [WARN] batch failed ({len(batch)} acc): {e}", file=sys.stderr)
            time.sleep(1.0)
            continue
        lines = text.strip("\n").split("\n")
        if len(lines) < 1:
            continue
        header = lines[0].split("\t")
        for l in lines[1:]:
            cols = l.split("\t")
            if len(cols) != len(header):
                continue
            row = dict(zip(header, cols))
            meta[row["Entry"]] = {
                "organism_name": row.get("Organism", ""),
                "lineage": row.get("Taxonomic lineage", ""),
                "length": row.get("Length", ""),
            }
        time.sleep(sleep)
        print(f"  ... {len(meta)} accessions annotees", file=sys.stderr)
    return meta


LINEAGE_RANKS = [
    "Actinomycetota", "Actinobacteria",  # phylum, deux noms possibles
]


def extract_order(lineage_str):
    """Extrait le taxon marque (order) dans le lignage UniProt 'Nom (rank), Nom (rank), ...'."""
    if not lineage_str:
        return None
    m = re.search(r"([A-Za-z0-9_.-]+) \(order\)", lineage_str)
    return m.group(1) if m else None


def main():
    seqs, ac_map = parse_stockholm(STO)
    print(f"{len(seqs)} sequences alignees, {len(ac_map)} accessions mappees", file=sys.stderr)

    if REF_ID not in seqs:
        raise SystemExit(f"reference {REF_ID} introuvable")
    ref_start_residue = int(REF_ID.split("/")[1].split("-")[0])
    boundary_col = find_module_boundary_column(seqs[REF_ID], MODULE_END_RESIDUE, ref_start_residue)
    print(f"colonne frontiere module/queue (residu {MODULE_END_RESIDUE}) = colonne {boundary_col}", file=sys.stderr)

    measurements = measure(seqs, boundary_col)

    accessions = [ac_map[sid] for sid in seqs if sid in ac_map]
    import os

    if os.path.exists(CACHE_JSON):
        print(f"cache UniProt trouve, lecture {CACHE_JSON}", file=sys.stderr)
        with open(CACHE_JSON) as fh:
            uniprot_meta = json.load(fh)
    else:
        print(f"interrogation UniProt REST pour {len(set(accessions))} accessions uniques...", file=sys.stderr)
        uniprot_meta = fetch_uniprot_metadata(accessions)
        with open(CACHE_JSON, "w") as fh:
            json.dump(uniprot_meta, fh)

    per_seq = []
    for sid, m in measurements.items():
        acc = ac_map.get(sid)
        um = uniprot_meta.get(acc, {})
        order = extract_order(um.get("lineage", ""))
        per_seq.append(
            {
                "id": sid,
                "accession": acc,
                "organism_code": organism_code(sid),
                "organism_name": um.get("organism_name"),
                "order": order,
                "uniprot_length": um.get("length"),
                **m,
            }
        )

    with open(OUT_TSV, "w") as fh:
        fh.write("id\taccession\torganism_code\torganism_name\torder\tuniprot_length\tmodule_len\ttail_len\ttotal_len\n")
        for r in per_seq:
            fh.write(
                f"{r['id']}\t{r['accession']}\t{r['organism_code']}\t{r['organism_name']}\t{r['order']}\t"
                f"{r['uniprot_length']}\t{r['module_len']}\t{r['tail_len']}\t{r['total_len']}\n"
            )

    by_order = defaultdict(list)
    for r in per_seq:
        if r["order"]:
            by_order[r["order"]].append(r["tail_len"])

    order_stats = {}
    for order, vals in by_order.items():
        if len(vals) >= 5:
            vals_sorted = sorted(vals)
            n = len(vals_sorted)
            median = vals_sorted[n // 2] if n % 2 else (vals_sorted[n // 2 - 1] + vals_sorted[n // 2]) / 2
            order_stats[order] = {
                "n": n,
                "median_tail_len": median,
                "min": min(vals_sorted),
                "max": max(vals_sorted),
            }

    n_annotated = sum(1 for r in per_seq if r["order"])
    summary = {
        "n_sequences": len(per_seq),
        "n_uniprot_annotated": len(uniprot_meta),
        "n_with_order": n_annotated,
        "boundary_column": boundary_col,
        "module_end_residue": MODULE_END_RESIDUE,
        "tail_len_distribution": {
            "min": min(r["tail_len"] for r in per_seq),
            "max": max(r["tail_len"] for r in per_seq),
        },
        "order_stats_sorted_by_median": dict(
            sorted(order_stats.items(), key=lambda kv: kv[1]["median_tail_len"])
        ),
    }

    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
