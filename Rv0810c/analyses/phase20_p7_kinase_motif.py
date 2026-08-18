"""P7 -- Attribution kinase de T24 par comparaison de motif consensus publie (Prisic et al. 2010).

Prisic S, Dankwa S, Schwartz D, Chou MF, Locasale JW, Kang CM, Bemis G, Church GM, Steen H,
Husson RN. "Extensive phosphorylation with overlapping specificity by Mycobacterium tuberculosis
serine/threonine protein kinases." PNAS 2010;107:7521-6. PMID 20368441, PMC2867705.

Motif consensus partage par les 6 kinases testees (PknA, PknB, PknD, PknE, PknF, PknH), cite
verbatim depuis le texte integral (section Results, autour de la Figure 2 / pLogo) :

    "the shared motif X-alpha-alpha-alpha-alpha-T-X-(X/V)-phi-(P/R)-I
    (where alpha is an acidic residue and phi a large hydrophobic residue)"

Positions relatives au T phosphoaccepteur (0) : -5 -4 -3 -2 -1 0 +1 +2 +3 +4 +5
Residus attendus                              :  X  a  a  a  a T  X X/V phi P/R  I

Nuances tirees du texte (citees dans pistes.md) :
- "acidic residues at positions N-terminal to the phosphoacceptor (-1 to -4)" = la meme region alpha.
- +4 : Pro prefere par PknA ET PknB ; Arg prefere par PknD, PknE, PknH.
- +5 : Ile predominant (mais la position tolere d'autres hydrophobes selon la figure).
- Le texte etablit explicitement que le motif partage LIMITE l'attribution directe d'un substrat a
  une kinase cognate ("the shared motif ... limits direct mapping of individual substrates to a
  cognate kinase").
"""
import json
from pathlib import Path

RV0810C = "MGRGRAKAKQTKVARELKYSSPQTDFQRLQRELSGTGTDRLDGDGPSDDDSWNDEDDWRR"
ACIDIC = set("DE")
LARGE_HYDROPHOBIC = set("LIVFMWY")  # phi : large hydrophobic side chain

OUT = Path(__file__).resolve().parent.parent / "résultats" / "p7_kinase_motif.json"


def window(seq, center_1based, half=5):
    """Fenetre -half..+half autour d'une position 1-based (numerotation Rv0810c)."""
    i = center_1based - 1  # 0-based
    out = {}
    for rel in range(-half, half + 1):
        idx = i + rel
        out[rel] = seq[idx] if 0 <= idx < len(seq) else None
    return out


def score_against_consensus(win):
    """Score position par position contre le motif publie. Renvoie le detail, pas juste un total."""
    checks = []

    def add(pos, expected_label, test_fn, residue):
        ok = test_fn(residue) if residue is not None else None
        checks.append({
            "position": pos,
            "residue": residue,
            "expected": expected_label,
            "match": ok,
        })

    add(-5, "X (libre)", lambda r: True, win[-5])
    add(-4, "acide (D/E)", lambda r: r in ACIDIC, win[-4])
    add(-3, "acide (D/E)", lambda r: r in ACIDIC, win[-3])
    add(-2, "acide (D/E)", lambda r: r in ACIDIC, win[-2])
    add(-1, "acide (D/E)", lambda r: r in ACIDIC, win[-1])
    add(0, "T (phosphoaccepteur)", lambda r: r == "T", win[0])
    add(1, "X (libre)", lambda r: True, win[1])
    add(2, "X ou V", lambda r: True, win[2])  # X/V : wildcard avec bonus V, teste comme wildcard
    add(3, "hydrophobe volumineux", lambda r: r in LARGE_HYDROPHOBIC, win[3])
    add(4, "Pro (PknA/PknB) ou Arg (PknD/E/H)", lambda r: r in ("P", "R"), win[4])
    add(5, "Ile (predominant)", lambda r: r == "I", win[5])

    # Sous-scores : la region la plus discriminante du motif publie est le run acide -1 a -4
    acidic_run = [c for c in checks if c["position"] in (-4, -3, -2, -1)]
    n_acidic_ok = sum(1 for c in acidic_run if c["match"])

    plus4 = next(c for c in checks if c["position"] == 4)
    plus4_subgroup = None
    if plus4["residue"] == "P":
        plus4_subgroup = "PknA/PknB (Pro)"
    elif plus4["residue"] == "R":
        plus4_subgroup = "PknD/PknE/PknH (Arg)"

    return {
        "checks": checks,
        "acidic_run_-4_to_-1_matches": f"{n_acidic_ok}/4",
        "plus4_subgroup_lean": plus4_subgroup,
        "plus3_hydrophobic_match": next(c for c in checks if c["position"] == 3)["match"],
        "plus5_ile_exact_match": next(c for c in checks if c["position"] == 5)["match"],
    }


def main():
    sites = {
        "T24": {"pos": 24, "n_datasets": 5, "note": "Verma2017, Carette2018, Malakar S2/S3b/S8 (P1.5)"},
        "S21": {"pos": 21, "n_datasets": 2, "note": "charniere du module, meme cluster que T24"},
        "S20": {"pos": 20, "n_datasets": 2, "note": "charniere du module, meme cluster que T24"},
        "S51": {"pos": 51, "n_datasets": 1, "note": "Verma seul, dans la queue acide desordonnee"},
    }

    result = {"consensus_source": "Prisic et al. 2010 PNAS, PMID 20368441, PMC2867705", "sites": {}}

    for name, meta in sites.items():
        win = window(RV0810C, meta["pos"])
        seq_str = "".join(r if r else "-" for r in [win[k] for k in range(-5, 6)])
        scored = score_against_consensus(win)
        result["sites"][name] = {
            "position_1based": meta["pos"],
            "n_datasets": meta["n_datasets"],
            "note": meta["note"],
            "window_-5_to_+5": seq_str,
            **scored,
        }

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    print(f"Ecrit : {OUT}\n")

    for name, d in result["sites"].items():
        print(f"=== {name} (position {d['position_1based']}, {d['n_datasets']} jeux) ===")
        print("fenetre -5..+5 :", d["window_-5_to_+5"])
        print("run acide -4..-1 :", d["acidic_run_-4_to_-1_matches"])
        print("+3 hydrophobe volumineux :", d["plus3_hydrophobic_match"])
        print("+4 penche vers :", d["plus4_subgroup_lean"])
        print("+5 Ile exact :", d["plus5_ile_exact_match"])
        print()


if __name__ == "__main__":
    main()
