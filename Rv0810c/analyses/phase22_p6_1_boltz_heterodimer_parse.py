#!/usr/bin/env python3
"""
P6.1 -- lecture des sorties Boltz-2 (heterodimeres Rv0810c + partenaires STRING
reels Rv0811c et Rv0812/PabC, + controle de specificite RpmG2/Rv0634B), tous les
modeles (reproductibilite). Meme logique que phase18_p3_4_boltz_parse.py (schema
Boltz != AF3, PAE inter-chaines recalcule depuis le .npz).
"""
import json
import warnings
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent.parent / "résultats" / "p6_1_boltz_heterodimers"
OUT_JSON = BASE / "p6_1_boltz_parsed.json"

JOBS = {
    "rv0810c_rpmg2_control_heterodimer": "controle de specificite -- RpmG2/Rv0634B, 50S ribosomal L33, 56 aa, sans raison biologique de s'associer",
    "rv0810c_rv0812_heterodimer": "test -- Rv0812/PabC, partenaire STRING (P0.5), voie PABA/folate",
    "rv0810c_rv0811c_heterodimer": "test -- Rv0811c, partenaire STRING (P0.5), voisin direct hypothetique",
}


def chain_tokens(cif_path):
    from Bio.PDB.MMCIFParser import MMCIFParser
    warnings.simplefilter("ignore")
    st = MMCIFParser(QUIET=True).get_structure("m", cif_path)
    toks = []
    for ch in list(st)[0]:
        for res in ch:
            if res.id[0] == " ":
                toks.append(ch.id)
            else:
                toks.extend([ch.id] * len(list(res)))
    return toks


def pae_cross_min(npz_path, cif_path):
    z = np.load(npz_path)
    pae = z["pae"] if "pae" in z.files else z[z.files[0]]
    pae = np.asarray(pae)
    if pae.ndim == 3:
        pae = pae[0]
    toks = np.array(chain_tokens(cif_path))
    if len(toks) != pae.shape[0]:
        return f"NA (tokens {len(toks)} != PAE {pae.shape[0]})"
    uniq = list(dict.fromkeys(toks.tolist()))
    if len(uniq) < 2:
        return None
    ia = np.where(toks == uniq[0])[0]
    ib = np.where(toks == uniq[1])[0]
    return round(float(pae[np.ix_(ia, ib)].min()), 3)


def interface_contacts(cif_path, cutoff=5.0):
    from Bio.PDB.MMCIFParser import MMCIFParser
    warnings.simplefilter("ignore")
    st = MMCIFParser(QUIET=True).get_structure("m", cif_path)
    chains = list(list(st)[0])
    if len(chains) < 2:
        return 0
    coords, labels = [], []
    for ch in chains:
        for res in ch:
            for at in res:
                coords.append(at.coord)
                labels.append((ch.id, res.id[1]))
    coords = np.array(coords)
    pairs = set()
    for i in range(len(coords)):
        ci, ri = labels[i]
        d = np.sqrt(((coords - coords[i]) ** 2).sum(1))
        for j in np.where(d < cutoff)[0]:
            cj, rj = labels[j]
            if ci != cj:
                pairs.add(tuple(sorted([(ci, ri), (cj, rj)])))
    return len(pairs)


def parse_job(job_name):
    pred_dir = BASE / f"out_{job_name}" / f"boltz_results_{job_name}" / "predictions" / job_name
    if not pred_dir.exists():
        return None
    models = []
    for conf_json in sorted(pred_dir.glob(f"confidence_{job_name}_model_*.json")):
        rank = conf_json.stem.split("_model_")[-1]
        d = json.load(open(conf_json))
        cif = pred_dir / f"{job_name}_model_{rank}.cif"
        npz = pred_dir / f"pae_{job_name}_model_{rank}.npz"
        cross = pae_cross_min(str(npz), str(cif)) if npz.exists() and cif.exists() else "NA"
        iface = interface_contacts(str(cif)) if cif.exists() else "NA"
        models.append({
            "model": int(rank),
            "iptm": d.get("iptm"),
            "ptm": d.get("ptm"),
            "complex_plddt": d.get("complex_plddt"),
            "confidence_score": d.get("confidence_score"),
            "pair_chains_iptm": d.get("pair_chains_iptm"),
            "pae_cross_chain_min_A": cross,
            "n_interface_contacts_5A": iface,
        })
    return models


def main():
    result = {}
    for job, desc in JOBS.items():
        models = parse_job(job)
        result[job] = {"description": desc, "models": models}
        if models:
            iptms = [m["iptm"] for m in models if isinstance(m["iptm"], (int, float))]
            print(f"\n=== {job} ({desc}) ===")
            for m in models:
                print(f"  model {m['model']}: ipTM={m['iptm']} pTM={m['ptm']} "
                      f"complex_plddt={m['complex_plddt']} PAE_cross_min={m['pae_cross_chain_min_A']} "
                      f"iface_contacts(<5A)={m['n_interface_contacts_5A']}")
            if iptms:
                print(f"  ipTM range: {min(iptms):.3f} - {max(iptms):.3f} (n={len(iptms)})")
        else:
            print(f"\n=== {job} : PAS ENCORE DE SORTIE (job non termine ou absent) ===")

    OUT_JSON.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    print(f"\nJSON: {OUT_JSON}")


if __name__ == "__main__":
    main()
