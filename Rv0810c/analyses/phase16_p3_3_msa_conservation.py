#!/usr/bin/env python3
"""
P3.3 -- MSA profond DUF3073 (data/PF11273_full.sto, 1930 sequences curees,
meme fichier que P2.1/P2.4) : conservation par position AU-DELA de la simple
identite basique, avec une mesure de CONTRAINTE RELATIVE (rang parmi les 59
positions du domaine), car le module 1-19/1-33 est deja quasi invariant en
bloc -- un simple comptage d'invariance n'y discrimine plus rien (cf. piste).

Question posee : au-dela des residus bases deja caracterises (P2.2, P2.3),
quels residus du module sont invariants a quasi 100% et NE sont PAS K/R --
et si un tel ensemble existe, forme-t-il une surface/poche exposee (meme
test de contiguite spatiale que P2.3, applique a ce nouvel ensemble) ?

Note de comptage (deja signalee en P3.3 / etat des decouvertes SS7) : le
chiffre de 3947 sequences cite dans la piste vient de l'analyse d'architecture
de domaines InterPro (tous hits, non alignes) ; l'alignement Pfam curated
PF11273_full utilise ici compte 1930 sequences -- c'est la seule des deux
sources qui soit un MSA exploitable pour une conservation par colonne.
"""
import itertools
import json
import math
import random
import statistics
import urllib.request
from collections import Counter
from pathlib import Path

import freesasa
from Bio.PDB import PDBParser

BASE = Path(__file__).resolve().parent.parent
STO = BASE / "data" / "PF11273_full.sto"
PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-I6XWB9-F1-model_v6.pdb"
PDB_LOCAL = BASE / "data" / "AF-I6XWB9-F1-model_v6.pdb"
OUT_JSON = BASE / "résultats" / "p3_3_msa_conservation.json"

REF_ID = "I6XWB9_MYCTU/2-60"
MODULE_END_RESIDUE = 33  # frontiere ordre/desordre P2.1
SEQ = "MGRGRAKAKQTKVARELKYSSPQTDFQRLQRELSGTGTDRLDGDGPSDDDSWNDEDDWRR"
BASIC = set("KR")

MAX_ASA = {
    "ALA": 106, "ARG": 248, "ASN": 157, "ASP": 163, "CYS": 135, "GLN": 198,
    "GLU": 194, "GLY": 84, "HIS": 184, "ILE": 169, "LEU": 164, "LYS": 205,
    "MET": 188, "PHE": 197, "PRO": 136, "SER": 130, "THR": 142, "TRP": 227,
    "TYR": 222, "VAL": 142,
}
AA3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}


def parse_stockholm(path):
    seqs = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or line.startswith("//") or not line.strip():
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                seqs[parts[0]] = parts[1]
    return seqs


def domain_columns(ref_row, ref_start_residue):
    """Liste (colonne, residu_H37Rv) pour chaque position non-gap de la reference."""
    out = []
    residue = ref_start_residue
    for col, ch in enumerate(ref_row):
        if ch not in (".", "-"):
            out.append((col, residue))
            residue += 1
    return out


def column_stats(seqs, col):
    counts = Counter()
    n_total = 0
    for row in seqs.values():
        n_total += 1
        if col < len(row):
            ch = row[col]
            if ch not in (".", "-"):
                counts[ch.upper()] += 1
    occ = sum(counts.values())
    return counts, occ, n_total


def shannon_bits(counts):
    total = sum(counts.values())
    if total == 0:
        return None
    h = 0.0
    for c in counts.values():
        p = c / total
        h -= p * math.log2(p)
    return h


def fetch_pdb():
    if not PDB_LOCAL.exists():
        PDB_LOCAL.parent.mkdir(exist_ok=True)
        urllib.request.urlretrieve(PDB_URL, PDB_LOCAL)
    return PDB_LOCAL


def get_rsa():
    pdb_path = fetch_pdb()
    fs_struct = freesasa.Structure(str(pdb_path))
    fs_result = freesasa.calc(fs_struct)
    residue_areas = fs_result.residueAreas()
    rsa = {}
    for _chain_id, res_map in residue_areas.items():
        for resnum_str, area in res_map.items():
            resnum = int(resnum_str)
            aa = SEQ[resnum - 1] if resnum - 1 < len(SEQ) else None
            max_asa = MAX_ASA.get(AA3.get(aa, ""))
            if max_asa:
                rsa[resnum] = 100.0 * area.total / max_asa
    return rsa, pdb_path


def mean_pairwise_ca_dist(residues_by_num, resnums):
    ca = {i: residues_by_num[i]["CA"] for i in resnums if i in residues_by_num and residues_by_num[i].has_id("CA")}
    pairs = list(itertools.combinations(ca, 2))
    if not pairs:
        return None
    return float(statistics.mean(float(ca[a] - ca[b]) for a, b in pairs))


def main():
    seqs = parse_stockholm(STO)
    if REF_ID not in seqs:
        raise SystemExit(f"reference {REF_ID} introuvable")
    ref_start_residue = int(REF_ID.split("/")[1].split("-")[0])
    cols = domain_columns(seqs[REF_ID], ref_start_residue)
    n_seq = len(seqs)
    print(f"{n_seq} sequences alignees ; {len(cols)} positions de domaine (residus {cols[0][1]}-{cols[-1][1]})")

    per_position = []
    for col, resnum in cols:
        counts, occ, n_total = column_stats(seqs, col)
        if not counts:
            continue
        maj_aa, maj_count = counts.most_common(1)[0]
        freq_maj = maj_count / occ
        entropy = shannon_bits(counts)
        per_position.append({
            "resnum": resnum,
            "h37rv_aa": SEQ[resnum - 1] if resnum - 1 < len(SEQ) else None,
            "majority_aa": maj_aa,
            "freq_majority": freq_maj,
            "occupancy": occ / n_total,
            "n_occupied": occ,
            "shannon_bits": entropy,
            "n_distinct_aa": len(counts),
            "is_module": resnum <= MODULE_END_RESIDUE,
            "is_basic_majority": maj_aa in BASIC,
        })

    # rang de contrainte relative : 1 = position la plus conservee (entropie la plus basse)
    # parmi les 59 positions du domaine -- meme convention que le rang T24 21/59 de P1.5
    ranked = sorted(per_position, key=lambda r: r["shannon_bits"])
    for rank, r in enumerate(ranked, start=1):
        r["conservation_rank"] = rank
    n_domain_positions = len(per_position)
    for r in per_position:
        r["conservation_percentile"] = 100.0 * (n_domain_positions - r["conservation_rank"] + 1) / n_domain_positions

    per_position.sort(key=lambda r: r["resnum"])

    # candidats : dans le module, majorite NON K/R, quasi invariants, bien echantillonnes
    STRICT_FREQ = 0.99
    BROAD_FREQ = 0.90
    MIN_OCC = 0.90
    candidates_strict = [
        r for r in per_position
        if r["is_module"] and not r["is_basic_majority"] and r["freq_majority"] >= STRICT_FREQ and r["occupancy"] >= MIN_OCC
    ]
    candidates_broad = [
        r for r in per_position
        if r["is_module"] and not r["is_basic_majority"] and r["freq_majority"] >= BROAD_FREQ and r["occupancy"] >= MIN_OCC
        and r not in candidates_strict
    ]

    rsa, pdb_path = get_rsa()
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("Rv0810c", str(pdb_path))
    chain = next(structure[0].get_chains())
    residues_by_num = {res.id[1]: res for res in chain if res.id[0] == " "}

    for r in candidates_strict + candidates_broad:
        r["RSA_percent"] = rsa.get(r["resnum"])
        r["exposed"] = rsa.get(r["resnum"], 0) >= 25.0

    exposed_candidates = [r["resnum"] for r in candidates_strict + candidates_broad if r.get("exposed")]

    # test de contiguite spatiale, meme protocole que P2.3 (phase14) : distance CA-CA
    # moyenne observee vs 2000 tirages aleatoires de meme taille parmi les residus
    # exposes du module (RSA >= 25%), toutes identites confondues
    module_exposed_all = [r["resnum"] for r in per_position if r["is_module"] and rsa.get(r["resnum"], 0) >= 25.0]
    observed_dist = mean_pairwise_ca_dist(residues_by_num, exposed_candidates)
    random.seed(2026)
    null_dists = []
    k = len(exposed_candidates)
    if k >= 2 and len(module_exposed_all) > k:
        for _ in range(2000):
            sample = random.sample(module_exposed_all, k)
            d = mean_pairwise_ca_dist(residues_by_num, sample)
            if d is not None:
                null_dists.append(d)
    null_mean = statistics.mean(null_dists) if null_dists else None
    null_sd = statistics.stdev(null_dists) if len(null_dists) > 1 else None
    z = (observed_dist - null_mean) / null_sd if (observed_dist is not None and null_mean and null_sd) else None
    p_le = (
        sum(1 for d in null_dists if d <= observed_dist) / len(null_dists)
        if null_dists and observed_dist is not None
        else None
    )

    # distance des candidats exposes au "patch" basique deja caracterise en P2.3
    basic_module_exposed = [r["resnum"] for r in per_position if r["is_module"] and r["h37rv_aa"] in BASIC and rsa.get(r["resnum"], 0) >= 25.0]
    cross_dists = []
    if exposed_candidates and basic_module_exposed:
        ca = {i: residues_by_num[i]["CA"] for i in set(exposed_candidates) | set(basic_module_exposed) if i in residues_by_num and residues_by_num[i].has_id("CA")}
        for a in exposed_candidates:
            for b in basic_module_exposed:
                if a in ca and b in ca:
                    cross_dists.append(float(ca[a] - ca[b]))
    mean_cross_dist = statistics.mean(cross_dists) if cross_dists else None

    # connu deja caracterise (P1.5) pour contexte, pas un resultat nouveau
    known_phosphosites = {24: "T24 (P1.5, 5 jeux)", 21: "S21 (P1.5, 2 jeux)", 20: "S20 (P1.5, 2 jeux)"}
    for r in per_position:
        r["known_phosphosite"] = known_phosphosites.get(r["resnum"])

    summary = {
        "n_sequences_alignment": n_seq,
        "n_domain_positions": n_domain_positions,
        "module_end_residue": MODULE_END_RESIDUE,
        "note_comptage": (
            "1930 sequences (alignement Pfam curated PF11273_full), a distinguer des 3947 "
            "sequences InterPro/architecture de domaines (tous hits, non alignes) citees dans "
            "la piste et dans etat_des_decouvertes.md SS7 -- seul le premier chiffre est un MSA."
        ),
        "thresholds": {
            "strict_freq_majority": STRICT_FREQ,
            "broad_freq_majority": BROAD_FREQ,
            "min_occupancy": MIN_OCC,
        },
        "candidates_strict_ge99pct_non_basic": [
            {k: v for k, v in r.items()} for r in candidates_strict
        ],
        "candidates_broad_90_99pct_non_basic": [
            {k: v for k, v in r.items()} for r in candidates_broad
        ],
        "n_exposed_candidates_RSA_ge_25pct": len(exposed_candidates),
        "exposed_candidate_resnums": exposed_candidates,
        "spatial_contiguity_test": {
            "method": "identique a P2.3/phase14 : distance CA-CA moyenne observee vs 2000 tirages aleatoires de meme taille parmi les residus exposes (RSA>=25%) du module",
            "observed_mean_pairwise_CA_dist_A": observed_dist,
            "null_n_draws": len(null_dists),
            "null_mean_pairwise_CA_dist_A": null_mean,
            "null_sd_pairwise_CA_dist_A": null_sd,
            "z_score": z,
            "empirical_p_observed_le_null": p_le,
        },
        "distance_to_known_basic_patch": {
            "basic_module_exposed_resnums": basic_module_exposed,
            "mean_CA_dist_candidates_to_basic_A": mean_cross_dist,
            "n_pairs": len(cross_dists),
        },
        "per_domain_position": per_position,
    }

    OUT_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    printable = {k: v for k, v in summary.items() if k != "per_domain_position"}
    print(json.dumps(printable, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
