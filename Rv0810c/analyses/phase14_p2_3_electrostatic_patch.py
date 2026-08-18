#!/usr/bin/env python3
"""
P2.3 — le module invariant 1-19 (et sa charniere 20-33) forme-t-il une SURFACE
electropositive CONTIGUE (signature structurale canonique d'un site de liaison
a un polyanion), ou les residus bases sont-ils exposes mais disperses sur des
faces differentes (pas de patch coherent) ?

Methode : modele AlphaFold local (recupere depuis AlphaFold DB, meme modele que
la fiche atlas et P2.1). Enfouissement = RSA% (freesasa, Shrake-Rupley, normalise
par les ASA max de Sander & Rost 1994), seuil expose = RSA >= 25%. Contiguite =
distance CA-CA moyenne entre residus bases EXPOSES, comparee a des tirages
aleatoires de meme taille parmi les residus exposes du module (null empirique).
"""
import itertools
import json
import random
import statistics
import urllib.request
from pathlib import Path

import freesasa
from Bio.PDB import PDBParser

# ASA max Sander & Rost 1994 (A^2), reference standard pour le RSA
MAX_ASA = {
    "ALA": 106, "ARG": 248, "ASN": 157, "ASP": 163, "CYS": 135, "GLN": 198,
    "GLU": 194, "GLY": 84, "HIS": 184, "ILE": 169, "LEU": 164, "LYS": 205,
    "MET": 188, "PHE": 197, "PRO": 136, "SER": 130, "THR": 142, "TRP": 227,
    "TYR": 222, "VAL": 142,
}

PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-I6XWB9-F1-model_v6.pdb"
PDB_LOCAL = Path(__file__).resolve().parent.parent / "data" / "AF-I6XWB9-F1-model_v6.pdb"
OUT_JSON = Path(__file__).resolve().parent.parent / "résultats" / "p2_3_electrostatic_patch.json"

SEQ = "MGRGRAKAKQTKVARELKYSSPQTDFQRLQRELSGTGTDRLDGDGPSDDDSWNDEDDWRR"
MODULE_RANGE = range(1, 34)  # residus 1-33 (P2.1)
BASIC = set("KR")
ACIDIC = set("DE")
CN_RADIUS = 10.0


def fetch_pdb():
    if not PDB_LOCAL.exists():
        PDB_LOCAL.parent.mkdir(exist_ok=True)
        urllib.request.urlretrieve(PDB_URL, PDB_LOCAL)
    return PDB_LOCAL


def get_side_chain_atom(residue):
    if residue.get_resname() == "GLY":
        return residue["CA"] if residue.has_id("CA") else None
    return residue["CB"] if residue.has_id("CB") else (residue["CA"] if residue.has_id("CA") else None)


def main():
    pdb_path = fetch_pdb()
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("Rv0810c", str(pdb_path))
    chain = next(structure[0].get_chains())
    residues = {r.id[1]: r for r in chain if r.id[0] == " "}

    atoms = {resnum: get_side_chain_atom(res) for resnum, res in residues.items()}
    atoms = {k: v for k, v in atoms.items() if v is not None}

    # SASA reelle (freesasa, algorithme Shrake-Rupley) -> RSA (%) par residu, Sander & Rost 1994
    fs_struct = freesasa.Structure(str(pdb_path))
    fs_result = freesasa.calc(fs_struct)
    residue_areas = fs_result.residueAreas()
    rsa = {}
    for chain_id, res_map in residue_areas.items():
        for resnum_str, area in res_map.items():
            resnum = int(resnum_str)
            resname = residues[resnum].get_resname() if resnum in residues else None
            max_asa = MAX_ASA.get(resname)
            if max_asa:
                rsa[resnum] = 100.0 * area.total / max_asa

    cn = rsa  # conserve le nom de variable en aval ; "cn" = RSA(%) desormais

    per_residue = []
    for i in sorted(atoms):
        aa = SEQ[i - 1] if i - 1 < len(SEQ) else "?"
        per_residue.append({"resnum": i, "aa": aa, "RSA_percent": rsa.get(i)})

    exposed_threshold = 25.0  # seuil standard RSA > 25% = expose (Rost & Sander 1994)

    basic_module = [i for i in MODULE_RANGE if i in cn and SEQ[i - 1] in BASIC]
    basic_module_exposed = [i for i in basic_module if cn[i] >= exposed_threshold]
    acidic_module = [i for i in MODULE_RANGE if i in cn and SEQ[i - 1] in ACIDIC]

    def mean_pairwise_ca_dist(resnums):
        ca = {i: residues[i]["CA"] for i in resnums if residues[i].has_id("CA")}
        pairs = list(itertools.combinations(ca, 2))
        if not pairs:
            return None
        return float(statistics.mean(float(ca[a] - ca[b]) for a, b in pairs))

    observed_dist = mean_pairwise_ca_dist(basic_module_exposed)

    # null : tirages aleatoires de meme taille parmi TOUS les residus exposes du module
    module_exposed_all = [i for i in MODULE_RANGE if i in cn and cn[i] >= exposed_threshold]
    random.seed(2026)
    null_dists = []
    k = len(basic_module_exposed)
    if k >= 2 and len(module_exposed_all) > k:
        for _ in range(2000):
            sample = random.sample(module_exposed_all, k)
            d = mean_pairwise_ca_dist(sample)
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

    summary = {
        "n_residues_modeled": len(atoms),
        "method": "freesasa (Shrake-Rupley) -> RSA%, seuil expose = RSA >= 25% (Rost & Sander 1994)",
        "module_range": [MODULE_RANGE.start, MODULE_RANGE.stop - 1],
        "exposed_threshold_RSA_percent": exposed_threshold,
        "basic_residues_module": [(i, SEQ[i - 1], round(cn.get(i, -1), 1)) for i in basic_module],
        "basic_residues_module_exposed_RSA_ge_25pct": [(i, SEQ[i - 1], round(cn[i], 1)) for i in basic_module_exposed],
        "acidic_residues_module": [(i, SEQ[i - 1], round(cn.get(i, -1), 1)) for i in acidic_module],
        "observed_mean_pairwise_CA_dist_basic_exposed_A": observed_dist,
        "null_n_draws": len(null_dists),
        "null_mean_pairwise_CA_dist_A": null_mean,
        "null_sd_pairwise_CA_dist_A": null_sd,
        "z_score_contiguity": z,
        "empirical_p_observed_le_null": p_le,
        "per_residue_RSA": per_residue,
        "interpretation_note": (
            "z tres negatif / p faible => les residus bases exposes du module sont PLUS PROCHES "
            "les uns des autres qu'un tirage aleatoire parmi les residus exposes du module : patch "
            "electropositif contigu, coherent avec un site de liaison a un polyanion. z proche de 0 "
            "=> pas plus regroupes que le hasard, les residus bases exposes sont disperses sur "
            "plusieurs faces, ce qui affaiblit l'hypothese d'un patch de liaison dedie."
        ),
    }

    OUT_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps({k: v for k, v in summary.items() if k != "per_residue_RSA"}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
