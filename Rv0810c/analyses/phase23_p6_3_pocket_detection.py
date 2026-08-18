#!/usr/bin/env python3
"""
P6.3 -- scan geometrique de cavites (P2Rank, fpocket) sur le modele AlphaFold de
Rv0810c (AF-I6XWB9-F1-model_v6.pdb, meme modele que P2.1/P2.3/P3.2), avec temoin
interne du projet (RpmG2/Rv0634B, deja utilise comme controle de specificite en
P3.4) et placement par rapport a la calibration ABSOLUE deja etablie ailleurs
(annotation_mtbc, P16.3b/P16.3h) sur un panel apparie de 25 enzymes prouvees vs
25 non-catalytiques, pLDDT_AF>=70, <200 aa -- panel qui couvre deja des longueurs
aussi courtes que 64 aa (Rv1642), donc directement applicable a Rv0810c (60 aa),
pas une extrapolation.

Garde-fou central du skill pocket-detection, verifie a nouveau ici : un score de
poche brut n'est interpretable qu'avec temoins positif ET negatif apparies en
taille. Sous 200 aa, P2Rank est mesure a 3,8% de detection confiante meme chez
des enzymes PROUVEES (contre 60,5% au-dela) ; fpocket et AE-PocketMiner ne
discriminent PAS non plus enzymes vs non-catalytiques sur ce segment (Mann-
Whitney p=0.82 et p=0.20/0.55, annotation_mtbc/résultats/phase99_p16_3h/rapport.md).//
Les deux outils sont donc a leur plancher documente avant meme d'etre lances sur
Rv0810c -- ce script consigne le resultat pour la completude de l'inventaire
methodologique (comme prevu par la piste), PAS comme un test au pouvoir statistique
attendu.

Entrees deja produites manuellement (outils installes dans annotation_mtbc/tools/,
reutilises tels quels) :
  résultats/p6_3_pocket_detection/p2rank/{rv0810c_out,rpmg2_out}/*_predictions.csv
  résultats/p6_3_pocket_detection/{rv0810c_out,rpmg2_rv0634b_out}/*_info.txt
Calibration externe lue directement dans annotation_mtbc (jamais copiee/dupliquee) :
  annotation_mtbc/résultats/phase99_p16_3h/{fpocket.tsv,rapport.md}
  annotation_mtbc/analyses/phase79_pocket_blind_zone.py (chiffre P2Rank <200aa)
"""
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "résultats" / "p6_3_pocket_detection"
CALIB_TSV = Path("/home/christophe/docs/codes/mtbc/annotation_mtbc/résultats/phase99_p16_3h/fpocket.tsv")
CALIB_RAPPORT = Path("/home/christophe/docs/codes/mtbc/annotation_mtbc/résultats/phase99_p16_3h/rapport.md")
OUT_JSON = OUT_DIR / "p6_3_pocket_detection.json"


def parse_p2rank(csv_path):
    if not csv_path.exists():
        return None
    rows = list(csv.DictReader(open(csv_path)))
    return {"n_pockets": len(rows), "pockets": rows}


def parse_fpocket_info(info_path):
    if not info_path.exists():
        return None
    text = info_path.read_text()
    pockets = []
    for block in text.split("Pocket ")[1:]:
        d = {}
        for line in block.splitlines():
            m = re.match(r"\s*(.+?)\s*:\s*\t?\s*([\-0-9.]+)\s*$", line)
            if m:
                d[m.group(1).strip()] = float(m.group(2))
        pockets.append(d)
    return {"n_pockets": len(pockets), "pockets": pockets}


def load_calibration():
    rows = list(csv.DictReader(open(CALIB_TSV))) if CALIB_TSV.exists() else []
    enz = [float(r["best_drug_score"]) for r in rows if r["group"] == "enzyme"]
    non = [float(r["best_drug_score"]) for r in rows if r["group"] == "noncat"]
    short_matches = [r for r in rows if int(r["len_aa"]) <= 100]
    return {
        "n_enzyme": len(enz),
        "n_noncat": len(non),
        "enzyme_drug_score_median": round(sorted(enz)[len(enz) // 2], 3) if enz else None,
        "noncat_drug_score_median": round(sorted(non)[len(non) // 2], 3) if non else None,
        "enzyme_drug_score_range": [round(min(enz), 3), round(max(enz), 3)] if enz else None,
        "noncat_drug_score_range": [round(min(non), 3), round(max(non), 3)] if non else None,
        "closest_length_matches_le100aa": short_matches,
        "source": str(CALIB_TSV),
        "p2rank_short_detection_rate": "3.8% (1/26) enzymes <200aa vs 60.5% (75/124) >=200aa -- annotation_mtbc/analyses/phase79_pocket_blind_zone.py",
        "fpocket_verdict": "NE DISCRIMINE PAS enzymes courtes vs non-catalytiques courtes (Mann-Whitney p=0.82 sur best_drug_score, Fisher p=1.0 sur presence/absence)",
    }


def main():
    result = {
        "rv0810c": {
            "p2rank": parse_p2rank(OUT_DIR / "p2rank" / "rv0810c_out" / "rv0810c.pdb_predictions.csv"),
            "fpocket": parse_fpocket_info(OUT_DIR / "rv0810c_out" / "rv0810c_info.txt"),
        },
        "rpmg2_control": {
            "p2rank": parse_p2rank(OUT_DIR / "p2rank" / "rpmg2_out" / "rpmg2_rv0634b.pdb_predictions.csv"),
            "fpocket": parse_fpocket_info(OUT_DIR / "rpmg2_rv0634b_out" / "rpmg2_rv0634b_info.txt"),
        },
        "calibration_externe_annotation_mtbc": load_calibration(),
    }

    r = result["rv0810c"]
    c = result["rpmg2_control"]
    print(f"Rv0810c    : P2Rank {r['p2rank']['n_pockets']} poche(s) | fpocket {r['fpocket']['n_pockets']} poche(s), "
          f"meilleur DS={max(p['Druggability Score'] for p in r['fpocket']['pockets']):.3f}")
    print(f"RpmG2/Rv0634B (controle) : P2Rank {c['p2rank']['n_pockets']} poche(s) | fpocket {c['fpocket']['n_pockets']} poche(s), "
          f"meilleur DS={max(p['Druggability Score'] for p in c['fpocket']['pockets']):.3f}")
    calib = result["calibration_externe_annotation_mtbc"]
    print(f"\nCalibration externe (annotation_mtbc, n={calib['n_enzyme']} enzymes / {calib['n_noncat']} non-cat, <200aa) :")
    print(f"  best_drug_score median enzymes={calib['enzyme_drug_score_median']} (range {calib['enzyme_drug_score_range']})")
    print(f"  best_drug_score median non-cat={calib['noncat_drug_score_median']} (range {calib['noncat_drug_score_range']})")
    print(f"  {calib['fpocket_verdict']}")
    print(f"  P2Rank : {calib['p2rank_short_detection_rate']}")

    OUT_JSON.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    print(f"\nJSON: {OUT_JSON}")


if __name__ == "__main__":
    main()
