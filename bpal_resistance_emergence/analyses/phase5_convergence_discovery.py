#!/usr/bin/env python3
"""Phase 5 — Découverte de déterminants par CONVERGENCE (homoplasie inter-lignées), phénotype-free.

Principe. Un variant non-synonyme qui apparaît de façon RÉCURRENTE dans plusieurs lignées majeures
phylogénétiquement distinctes, à FRÉQUENCE INTRA-LIGNÉE FAIBLE (sporadique, pas fixé), est le signe
d'une émergence indépendante répétée — signature classique de la sélection médicamenteuse (ex.
gyrA D94G/N/A sous fluoroquinolone). C'est la seule voie d'attaque pour le prétomanide (0 phénotype) :
on cherche les positions des gènes d'activation F420 (ddn/fgd1/fbiA-D) frappées indépendamment.

Discrimination (garde-fou trichotomie, cf. Phase 4) :
  - marqueur de lignée  : concentré >=85% dans UNE lignée, fixé (>=5% de cette lignée) -> PAS convergence
  - convergent récurrent: présent dans >=3 lignées majeures à freq intra-lignée < 10% (sporadique)
  - LoF convergente     : stop_gained/frameshift récurrent -> perte de fonction (mécanisme DLM/PMD établi,
                          Beckert 2020 / Zhang 2026), INTERPRÉTABLE SANS localisation structurale

Unité d'indépendance : lignée majeure (L1,L2,L3,L4,L5,L6,L7,BOV,BOV_AFRI). On rapporte AUSSI
`n_independent_clades` (fonction validée du skill convergent-evolution : regroupe le clade humain
strict), mais on CLASSE/range par `n_lin_sporadic` car la convergence de résistance est surtout
INTRA-humaine (que l'encodage 5-clades écrase).

Contrôle positif OBLIGATOIRE : sans masquer le catalogue, gyrA pos 94 / 90 (+ gyrB) doivent
dominer la convergence FQ. Si oui -> méthode validée, on interprète alors les hits F420/BDQ
hors-catalogue.

Niveau N1/N2 (cadre tuberculosis.md) : proxy d'homoplasie par étalement-lignée + fréquence ; la
confirmation per-émergence (pastml/SNPPar sur l'arbre) reste un N2 différé. Source de la prudence :
l'étalement-lignée à granularité majeure ne distingue pas une émergence sous-clonale unique d'une
récurrence ; la freq intra-lignée FAIBLE est le garde-fou (un marqueur sous-clonal unique serait à
freq plus élevée dans sa lignée).

Sorties : résultats/phase5_convergence_variants.tsv, phase5_convergence_positions.tsv + récap console.
"""
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths

# fonction d'indépendance validée du skill (réutilisation, pas réimplémentation)
sys.path.insert(0, "/home/christophe/docs/codes/claude_plugins/bio_pathogens/"
                    "skills/convergent-evolution/scripts")
from convergent_evolution import count_independent_lineages  # noqa: E402

FIX_CEIL = 0.10          # < 10% intra-lignée = présence sporadique/récurrente
MARKER_CONC = 0.85       # >=85% des porteurs dans 1 lignée = marqueur concentré
MARKER_FIX = 0.05        # ET >=5% de cette lignée = fixé -> marqueur de lignée
MIN_SPORADIC = 3         # >=3 lignées majeures sporadiques = convergent récurrent
F420 = {"ddn", "fgd1", "fbiA", "fbiB", "fbiC", "fbiD"}
BDQ = {"Rv0678", "mmpR5", "atpE", "pepQ"}


def parse_breakdown(s):
    br = {}
    for tok in (s or "").split(";"):
        if ":" in tok:
            k, v = tok.rsplit(":", 1)
            if v.strip().isdigit():
                br[k] = int(v)
    br.pop("NA", None)
    return br


def aa_pos(aa_change):
    m = re.search(r"(\d+)", aa_change or "")
    return int(m.group(1)) if m else None


def main():
    lin_total = {}
    with open(paths.RESULTATS / "phase3_lineage_prevalence.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lin_total[r["lineage"]] = int(r["n_total"])

    rows = []
    # agrégation position-level (gène, position) : union des lignées sporadiques, résidus alt distincts
    pos_agg = defaultdict(lambda: {"lin_spor": set(), "lin_pres": set(), "alts": set(),
                                   "group": None, "drug": None, "n_carriers": 0, "catR_any": False})

    with open(paths.RESULTATS / "phase1_feasibility.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            eff = r["effect"]
            if eff not in ("missense", "stop_gained"):
                continue
            br = parse_breakdown(r["lineage_breakdown"])
            if not br:
                continue
            tot = sum(br.values())
            freqs = {l: br[l] / lin_total[l] for l in br if lin_total.get(l, 0) > 0}
            if not freqs:
                continue
            present = set(freqs)
            sporadic = {l for l, f in freqs.items() if 0 < f < FIX_CEIL}
            maxf = max(freqs.values())
            dom = max(freqs, key=lambda k: freqs[k])
            conc = br[dom] / tot if tot else 0
            n_indep = count_independent_lineages(present)
            catR = (r["is_catalogued_R"] == "1")

            if conc >= MARKER_CONC and (br[dom] / lin_total.get(dom, 1)) >= MARKER_FIX:
                cls = "lineage_marker"
            elif eff in ("stop_gained",) and len(sporadic) >= 2:
                cls = "LoF_convergent"
            elif len(sporadic) >= MIN_SPORADIC:
                cls = "convergent_recurrent"
            elif len(present) >= 2 and maxf < FIX_CEIL:
                cls = "sporadic_lowspread"
            else:
                cls = "single_or_fixed"

            pos = aa_pos(r["aa_change"])
            rows.append({
                "group": r["group"], "drug": r["drug"], "gene": r["gene"],
                "aa_change": r["aa_change"], "pos": pos, "effect": eff, "spdi": r["spdi"],
                "n_carriers": tot, "n_lin_present": len(present), "n_lin_sporadic": len(sporadic),
                "n_independent_clades": n_indep, "max_intra_freq_pct": round(100 * maxf, 3),
                "dom_lineage": dom, "lineages_sporadic": ",".join(sorted(sporadic)),
                "class": cls, "is_catalogued_R": int(catR), "catalogue": r["catalogue"] or "-",
            })

            if pos is not None:
                key = (r["gene"], pos)
                a = pos_agg[key]
                a["group"] = r["group"]; a["drug"] = r["drug"]
                a["lin_spor"] |= sporadic; a["lin_pres"] |= present
                alt = r["aa_change"][-1] if r["aa_change"] else "?"
                a["alts"].add(alt); a["n_carriers"] += tot
                a["catR_any"] = a["catR_any"] or catR

    # --- sortie variant-level ---
    rows.sort(key=lambda x: (-x["n_lin_sporadic"], -x["n_carriers"]))
    vcols = ["group", "drug", "gene", "aa_change", "pos", "effect", "spdi", "n_carriers",
             "n_lin_present", "n_lin_sporadic", "n_independent_clades", "max_intra_freq_pct",
             "dom_lineage", "lineages_sporadic", "class", "is_catalogued_R", "catalogue"]
    with open(paths.RESULTATS / "phase5_convergence_variants.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=vcols, delimiter="\t")
        w.writeheader(); w.writerows(rows)

    # --- sortie position-level ---
    pos_rows = []
    for (gene, pos), a in pos_agg.items():
        n_spor = len(a["lin_spor"])
        if n_spor < 2 and len(a["alts"]) < 2:
            continue
        pos_rows.append({
            "group": a["group"], "drug": a["drug"], "gene": gene, "pos": pos,
            "n_distinct_alt": len(a["alts"]), "alts": "".join(sorted(a["alts"])),
            "n_lin_sporadic": n_spor, "lineages_sporadic": ",".join(sorted(a["lin_spor"])),
            "n_independent_clades": count_independent_lineages(a["lin_spor"]) if a["lin_spor"] else 0,
            "n_carriers": a["n_carriers"], "is_catalogued_R": int(a["catR_any"]),
        })
    pos_rows.sort(key=lambda x: (-x["n_lin_sporadic"], -x["n_distinct_alt"], -x["n_carriers"]))
    pcols = ["group", "drug", "gene", "pos", "n_distinct_alt", "alts", "n_lin_sporadic",
             "lineages_sporadic", "n_independent_clades", "n_carriers", "is_catalogued_R"]
    with open(paths.RESULTATS / "phase5_convergence_positions.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=pcols, delimiter="\t")
        w.writeheader(); w.writerows(pos_rows)

    # ============ RÉCAP CONSOLE ============
    def show(title, sel, key="var", n=12):
        print(f"\n=== {title} ===")
        if key == "var":
            print(f"{'gene':7}{'aa':9}{'eff':12}{'nspor':>6}{'nclade':>7}{'maxf%':>8}{'nport':>7}  class / lignées")
            for x in sel[:n]:
                print(f"{x['gene']:7}{x['aa_change']:9}{x['effect']:12}{x['n_lin_sporadic']:6}"
                      f"{x['n_independent_clades']:7}{x['max_intra_freq_pct']:8.2f}{x['n_carriers']:7}"
                      f"  {x['class']} [{x['lineages_sporadic']}]")
        else:
            print(f"{'gene':7}{'pos':>5}{'nalt':>5} {'alts':6}{'nspor':>6}{'nclade':>7}{'nport':>7}  lignées")
            for x in sel[:n]:
                print(f"{x['gene']:7}{x['pos']:5}{x['n_distinct_alt']:5} {x['alts']:6}"
                      f"{x['n_lin_sporadic']:6}{x['n_independent_clades']:7}{x['n_carriers']:7}"
                      f"  [{x['lineages_sporadic']}]")

    # CONTRÔLE POSITIF : sites FQ catalogués (gyrA/gyrB) — doivent dominer
    ctrl = [x for x in pos_rows if x["gene"] in ("gyrA", "gyrB") and x["is_catalogued_R"]]
    show("CONTRÔLE POSITIF — positions gyrA/gyrB catalogués-R (convergence FQ attendue en tête)",
         ctrl, key="pos")

    # DÉCOUVERTE : F420 (PMD/DLM) hors-catalogue, convergent
    f420v = [x for x in rows if x["gene"] in F420 and not x["is_catalogued_R"]
             and x["class"] in ("convergent_recurrent", "LoF_convergent")]
    show("DÉCOUVERTE — F420 (ddn/fgd1/fbiA-D) hors-catalogue, convergent [missense=à localiser, "
         "stop=LoF directe]", f420v, key="var")

    f420p = [x for x in pos_rows if x["gene"] in F420 and not x["is_catalogued_R"]]
    show("DÉCOUVERTE — positions F420 hors-catalogue frappées dans plusieurs lignées", f420p, key="pos")

    # BDQ hors-catalogue convergent
    bdqv = [x for x in rows if x["gene"] in BDQ and not x["is_catalogued_R"]
            and x["class"] in ("convergent_recurrent", "LoF_convergent")]
    show("DÉCOUVERTE — BDQ (Rv0678/atpE/pepQ) hors-catalogue, convergent", bdqv, key="var")

    # bilan compteurs
    from collections import Counter
    cc = Counter(x["class"] for x in rows)
    print(f"\nClasses (tous variants non-syn du panel) : {dict(cc)}")
    print(f"Variants -> {paths.RESULTATS}/phase5_convergence_variants.tsv")
    print(f"Positions -> {paths.RESULTATS}/phase5_convergence_positions.tsv")
    print("\nProchaine brique : localisation structurale (distance F420) des hits F420 MISSENSE "
          "convergents (un site convergent ET < 8 Å du cofacteur = candidat de résistance fort) ; "
          "les stop_gained convergents sont déjà interprétables comme LoF.")


if __name__ == "__main__":
    main()
