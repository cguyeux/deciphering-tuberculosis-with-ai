#!/usr/bin/env python3
"""Phase 5b — Validation DST (délamanide = proxy F420 du PMD) des candidats convergents F420.

Trois volets, sur la cohorte délamanide du consensus Resistance_antibio (187 R / 11 739 S) :
  (A) DIRECT : les porteurs des variants convergents fgd1 10/210 ont-ils un DST délamanide, R ou S ?
      (attendu sous-puissant : ces variants standing sont quasi absents des cohortes DST cliniques).
  (B) CAS-TÉMOINS voie F420 : parmi les souches DLM-phénotypées, la PORTANCE d'un variant des gènes
      F420 (ddn/fgd1/fbiA-D) est-elle enrichie chez les R vs S ? Fisher exact unilatéral.
      - CONTRÔLE POSITIF : les variants catalogués DLM-R doivent être enrichis -> valide cohorte+méthode.
      - puis : un variant fgd1 de la poche (pos 10/210), hors-catalogue, montre-t-il un signal ?
  (C) INVERSE : quelles positions des gènes F420 les 187 souches DLM-R portent-elles (histogramme) ?
      Test : la région de la poche (incl. 10/210) est-elle frappée indépendamment chez les vrais R ?

Garde-fou (cadre N2, tuberculosis.md) : DLM-R est enrichi en L2 (Phase 3b) -> confond lignée. Effectifs
faibles -> on RAPPORTE les lignées des porteurs-R, on ne prétend pas à une association à lignée contrôlée
puissante. Cette phase CORROBORE/réfute, elle ne prouve pas. Jamais de fréquence à partir de cette cohorte.

Sortie : résultats/phase5b_dst_validation.tsv + récap console.
"""
import csv
import pickle
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths

F420 = {"ddn", "fgd1", "fbiA", "fbiB", "fbiC", "fbiD"}
CAND = {  # candidats convergents Phase 5 (fgd1 poche)
    "NC_000962.3:490810:C:T": "fgd1 A10V", "NC_000962.3:490809:G:C": "fgd1 A10P",
    "NC_000962.3:490809:G:A": "fgd1 A10T", "NC_000962.3:491409:G:A": "fgd1 A210T",
    "NC_000962.3:491409:G:C": "fgd1 A210P", "NC_000962.3:491410:C:A": "fgd1 A210E",
    "NC_000962.3:491410:C:G": "fgd1 A210G", "NC_000962.3:491410:C:T": "fgd1 A210V",
}


def aa_pos(s):
    m = re.search(r"(\d+)", s or "")
    return int(m.group(1)) if m else None


def main():
    # --- phénotypes délamanide ---
    dlm = {}
    with open(paths.PHENOTYPES_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if "delamanid" in r["drug"].lower():
                p = r["phenotype"].strip().upper()
                if p in ("R", "S"):
                    dlm[r["strain_id"]] = p
    R = {s for s, p in dlm.items() if p == "R"}
    S = {s for s, p in dlm.items() if p == "S"}
    print(f"Cohorte délamanide : {len(R)} R / {len(S)} S")

    # --- lignées (pour annoter les porteurs-R) ---
    lin = {}
    with open(paths.LINEAGE_SNAPSHOT) as fh:
        rd = csv.DictReader(fh)
        f = rd.fieldnames
        sk = "strain_name" if "strain_name" in f else f[0]
        lk = "lineage_level_1" if "lineage_level_1" in f else ("lineage_code" if "lineage_code" in f else f[-1])
        for r in rd:
            lin[r[sk]] = r[lk]

    # --- variants F420 du panel (phase1) ---
    f420_vars = []
    with open(paths.RESULTATS / "phase1_feasibility.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene"] in F420 and r["effect"] in ("missense", "stop_gained"):
                cat = r["catalogue"] or ""
                f420_vars.append({
                    "gene": r["gene"], "aa": r["aa_change"], "pos": aa_pos(r["aa_change"]),
                    "spdi": r["spdi"], "effect": r["effect"],
                    "dlm_R_assoc": "delamanid:R-associated" in cat, "catalogue": cat or "-",
                })

    print("chargement pan_spdi…", flush=True)
    pan = pickle.load(open(paths.PAN_SPDI, "rb"))

    def carriers(spdi):
        return set(pan.get(spdi, []))

    # ===== (A) DIRECT — candidats fgd1 10/210 =====
    print("\n=== (A) DIRECT : porteurs des candidats convergents fgd1 10/210 dans la cohorte DLM ===")
    print(f"{'variant':12}{'porteurs':>9}{'DLM-testés':>11}{'R':>4}{'S':>4}  lignées des R")
    for spdi, name in CAND.items():
        car = carriers(spdi)
        tested = car & (R | S)
        rr = car & R
        ss = car & S
        rl = ",".join(sorted({lin.get(s, "?") for s in rr})) if rr else ""
        print(f"{name:12}{len(car):9}{len(tested):11}{len(rr):4}{len(ss):4}  {rl}")

    # ===== (B) CAS-TÉMOINS voie F420 =====
    nR, nS = len(R), len(S)
    rows = []
    for v in f420_vars:
        car = carriers(v["spdi"])
        a = len(car & R); b = len(car & S)
        if a + b < 2:   # au moins 2 porteurs DLM-testés
            continue
        c = nR - a; d = nS - b
        try:
            orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        except Exception:
            orr, p = float("nan"), 1.0
        rl = ",".join(sorted(Counter(lin.get(s, "?") for s in (car & R)).keys())) if a else ""
        rows.append({**v, "carr_R": a, "carr_S": b, "OR": round(orr, 2) if orr == orr else "inf",
                     "p": p, "R_lineages": rl})
    rows.sort(key=lambda x: x["p"])

    with open(paths.RESULTATS / "phase5b_dst_validation.tsv", "w") as fh:
        cols = ["gene", "aa", "pos", "effect", "spdi", "dlm_R_assoc", "catalogue",
                "carr_R", "carr_S", "OR", "p", "R_lineages"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    pos_ctrl = [x for x in rows if x["dlm_R_assoc"]]
    print(f"\n=== (B) CONTRÔLE POSITIF : variants F420 catalogués DLM-R, enrichissement chez les R ===")
    print(f"{'gene':6}{'aa':9}{'cR':>4}{'cS':>4}{'OR':>8}{'p':>10}  lignées-R")
    for x in pos_ctrl[:12]:
        print(f"{x['gene']:6}{x['aa']:9}{x['carr_R']:4}{x['carr_S']:4}{str(x['OR']):>8}{x['p']:10.1e}  {x['R_lineages']}")

    print(f"\n=== (B) CANDIDATS hors-catalogue les plus enrichis chez les DLM-R (top par p) ===")
    print(f"{'gene':6}{'aa':9}{'pos':>5}{'cR':>4}{'cS':>4}{'OR':>8}{'p':>10}  lignées-R")
    for x in [r for r in rows if not r["dlm_R_assoc"]][:15]:
        star = " *POCHE10/210*" if (x["gene"] == "fgd1" and x["pos"] in (10, 210)) else ""
        print(f"{x['gene']:6}{x['aa']:9}{str(x['pos']):>5}{x['carr_R']:4}{x['carr_S']:4}"
              f"{str(x['OR']):>8}{x['p']:10.1e}  {x['R_lineages']}{star}")

    # ===== (C) INVERSE — positions F420 portées par les 187 DLM-R =====
    spdi2var = {v["spdi"]: v for v in f420_vars}
    pos_hits = Counter()
    gene_pos_hits = defaultdict(Counter)
    for v in f420_vars:
        car = carriers(v["spdi"])
        nr = len(car & R)
        if nr and v["pos"] is not None:
            pos_hits[(v["gene"], v["pos"])] += nr
            gene_pos_hits[v["gene"]][v["pos"]] += nr
    print(f"\n=== (C) INVERSE : positions F420 les plus portées par les {nR} souches DLM-R ===")
    print(f"{'gene':6}{'pos':>5}{'n_R_porteurs':>14}  (poche fgd1 10/210 ?)")
    for (g, p), n in pos_hits.most_common(15):
        mark = " <-- POCHE" if (g == "fgd1" and p in (10, 210)) else ""
        print(f"{g:6}{p:5}{n:14}{mark}")

    print("\nLecture : (A) sous-puissant attendu (candidats standing ~absents des cohortes DST). "
          "(B) le contrôle positif catalogué doit s'allumer (sinon cohorte/jointure cassée) ; un candidat "
          "hors-catalogue enrichi chez les R sur PLUSIEURS lignées serait corroborant. (C) si la poche "
          "fgd1 10/210 est frappée chez des R indépendants, c'est convergent ET phénotype-associé. "
          "Garde-fou : confond L2 ; effectifs faibles ; corroboration, pas preuve.")
    print(f"\nsortie -> {paths.RESULTATS}/phase5b_dst_validation.tsv")


if __name__ == "__main__":
    main()
