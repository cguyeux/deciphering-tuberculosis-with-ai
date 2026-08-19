#!/usr/bin/env python3
"""Phase 3b — Correction du biais d'échantillonnage de l'enrichissement L2 (Phase 3).

La prévalence BRUTE de Phase 3 (L2/Beijing OR 1.8-5.3 pour les déterminants R) est confondue par
(i) la SÉLECTION : les études MDR (BioProjects cliniques) sur-échantillonnent Beijing et déposent
de gros lots résistants ; (ii) la PSEUDO-RÉPLICATION : des souches d'un même BioProject ne sont pas
indépendantes ; (iii) la CURATION hétérogène de bdd/actuelle par lignée (KB : L1/Bovis peu curés,
L2/L3 moyens, L4/L5/L6/animales fortement curés -> dénominateurs non comparables).

Trois corrections complémentaires (aucune n'est causale en soi ; ensemble = sensibilité honnête) :
  A. IPW par BioProject : poids = 1/(taille du BioProject) -> chaque étude pèse ~1 unité, ce qui
     neutralise les gros lots MDR. Prévalence + OR L2-vs-reste pondérés.
  B. Curation comparable : L2 vs L3 (toutes deux « moyennes ») = le contraste le moins confondu.
  C. Logistique cluster-robuste sur BioProject (statsmodels) : OR L2 vs L4 avec IC corrigés de la
     corrélation intra-étude.
  + Fragilité : nombre de BioProjects DISTINCTS derrière les porteurs R de chaque lignée.

Sorties : résultats/phase3b_ipw_adjusted.tsv + article/figures/phase3b_crude_vs_ipw.{png,pdf}.
Usage : python phase3b_ipw.py [--limit N]
"""
import argparse, csv, glob, pickle, sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths

RDRUGS = ["bedaquiline", "clofazimine", "delamanid", "moxifloxacin", "linezolid", "levofloxacin"]
MAJORS = ["L1", "L2", "L3", "L4", "L5", "L6", "L7", "BOV", "BOV_AFRI"]
# Intensité de curation par lignée (KB reference_bdd_curation_intensity_by_lineage)
CURATION = {"L1": "low", "BOV": "low", "BOV_AFRI": "low", "L2": "med", "L3": "med",
            "L4": "high", "L5": "high", "L6": "high", "L7": "high"}


def lin_label(x):
    if x in (None, "", "NA"):
        return "NA"
    return x if x in ("BOV", "BOV_AFRI") else "L" + str(x)


def weighted_or(a, b, c, d):
    """OR d'un 2x2 (poss. pondéré, fractionnaire) avec correction de Haldane si une case nulle."""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + .5, b + .5, c + .5, d + .5
    return (a * d) / (b * c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int)
    args = ap.parse_args()
    paths.ensure_dirs()
    figdir = paths.ARTICLE / "figures"; figdir.mkdir(parents=True, exist_ok=True)

    # ensembles catalogués-R par drogue (depuis phase1)
    catR_drug = defaultdict(set)
    with open(paths.RESULTATS / "phase1_feasibility.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for tok in r["catalogue"].split(";"):
                if tok.endswith(":R-associated"):
                    catR_drug[tok.split(":")[0]].add(r["spdi"])
    catR_drug = {d: s for d, s in catR_drug.items() if d in RDRUGS}

    strain2lin = {}
    with open(paths.LINEAGE_SNAPSHOT) as fh:
        for row in csv.DictReader(fh):
            strain2lin[row["strain_name"]] = lin_label(row["lineage_level_1"])
    strain2bp = {}
    for f in glob.glob(str(paths.BIOPROJECT_GEO / "consolidated_geo_*.tsv")):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                sra = row.get("strain")
                bp = (row.get("bioproject") or "").strip()
                if sra and bp:
                    strain2bp[sra] = bp
    print(f"snapshot lignée {len(strain2lin)} | bioproject {len(strain2bp)} souches")

    print("chargement pan_strains.pkl…")
    with open(paths.PAN_STRAINS, "rb") as f:
        pan = pickle.load(f)

    # passage : retenir souches lignée-majeure ET bioproject connus
    rows = []                       # (lineage, bp, frozenset R-drugs)
    for i, (strain, spdis) in enumerate(pan.items()):
        if args.limit and i >= args.limit:
            break
        lab = strain2lin.get(strain, "NA")
        bp = strain2bp.get(strain)
        if lab not in MAJORS or not bp:
            continue
        carried = spdis if isinstance(spdis, set) else set(spdis)
        rd = frozenset(d for d, sset in catR_drug.items() if carried & sset)
        rows.append((lab, bp, rd))
    print(f"souches analysables (lignée majeure + bioproject) : {len(rows)}")
    bp_size = Counter(bp for _, bp, _ in rows)

    # --- par drogue : crude / IPW / curation / fragilité ---------------------
    out = paths.RESULTATS / "phase3b_ipw_adjusted.tsv"
    fh = open(out, "w")
    fh.write("drug\tlineage\tn\tn_R\tcrude_prev_pct\tipw_prev_pct\tn_bioprojects_R\n")
    summary = {}
    prev_crude = defaultdict(dict); prev_ipw = defaultdict(dict)
    for drug in sorted(catR_drug):
        # comptes par lignée
        n = Counter(); nR = Counter()
        wn = defaultdict(float); wnR = defaultdict(float)
        bpR = defaultdict(set)
        for lab, bp, rd in rows:
            w = 1.0 / bp_size[bp]
            isR = drug in rd
            n[lab] += 1; wn[lab] += w
            if isR:
                nR[lab] += 1; wnR[lab] += w; bpR[lab].add(bp)
        for lab in MAJORS:
            if n[lab] == 0:
                continue
            cp = 100 * nR[lab] / n[lab]
            ip = 100 * wnR[lab] / wn[lab] if wn[lab] else 0
            prev_crude[drug][lab] = cp; prev_ipw[drug][lab] = ip
            fh.write(f"{drug}\t{lab}\t{n[lab]}\t{nR[lab]}\t{cp:.4f}\t{ip:.4f}\t{len(bpR[lab])}\n")
        # OR L2 vs reste : crude + IPW
        l2, oth = "L2", [l for l in MAJORS if l != "L2" and n[l] > 0]
        a, b = nR[l2], n[l2] - nR[l2]
        c = sum(nR[l] for l in oth); d = sum(n[l] - nR[l] for l in oth)
        or_crude = weighted_or(a, b, c, d)
        wa, wb = wnR[l2], wn[l2] - wnR[l2]
        wc = sum(wnR[l] for l in oth); wd = sum(wn[l] - wnR[l] for l in oth)
        or_ipw = weighted_or(wa, wb, wc, wd)
        # OR L2 vs L3 (curation comparable « med »)
        a3, b3 = nR["L2"], n["L2"] - nR["L2"]
        c3, d3 = nR.get("L3", 0), n.get("L3", 0) - nR.get("L3", 0)
        or_l2l3 = weighted_or(a3, b3, c3, d3) if n.get("L3", 0) else float("nan")
        summary[drug] = (or_crude, or_ipw, or_l2l3, len(bpR["L2"]))
    fh.close()

    # --- logistique cluster-robuste (statsmodels) : OR L2 vs L4 --------------
    cluster = {}
    try:
        import pandas as pd, statsmodels.formula.api as smf
        df = pd.DataFrame({"lin": [r[0] for r in rows], "bp": [r[1] for r in rows]})
        for drug in sorted(catR_drug):
            df["y"] = [1 if drug in r[2] else 0 for r in rows]
            sub = df[df["lin"].isin(["L1", "L2", "L3", "L4"])].copy()  # humaines majeures
            try:
                m = smf.logit("y ~ C(lin, Treatment(reference='L4'))", sub).fit(
                    disp=0, cov_type="cluster", cov_kwds={"groups": sub["bp"]})
                import numpy as np
                coef = m.params.get("C(lin, Treatment(reference='L4'))[T.L2]", np.nan)
                ci = m.conf_int().loc["C(lin, Treatment(reference='L4'))[T.L2]"] \
                    if "C(lin, Treatment(reference='L4'))[T.L2]" in m.params else [np.nan, np.nan]
                cluster[drug] = (np.exp(coef), np.exp(ci[0]), np.exp(ci[1]))
            except Exception as e:
                cluster[drug] = None
    except Exception as e:
        print(f"[logit] statsmodels indisponible : {e}")

    # --- rapport -------------------------------------------------------------
    print("\n=== L2 enrichissement : CRUDE -> AJUSTÉ (OR L2 vs reste) ===")
    print(f"{'drug':13} {'OR_crude':>9} {'OR_IPW':>8} {'OR_L2vsL3':>10} {'nBioProj_L2R':>12} "
          f"{'OR_clust(vsL4)[IC95]':>26}")
    for drug in sorted(catR_drug):
        oc, oi, ol, nbp = summary[drug]
        cl = cluster.get(drug)
        clstr = f"{cl[0]:.2f} [{cl[1]:.2f}-{cl[2]:.2f}]" if cl else "n/a"
        print(f"{drug:13} {oc:9.2f} {oi:8.2f} {ol:10.2f} {nbp:12d} {clstr:>26}")

    print("\nLecture : OR_crude = Phase 3 (confondu). OR_IPW = études pondérées à poids égal "
          "(neutralise gros lots MDR). OR_L2vsL3 = curation comparable (le plus honnête). "
          "OR_clust = logistique vs L4, IC corrigés de la corrélation intra-BioProject.")

    try:
        make_fig(prev_crude, prev_ipw, catR_drug, figdir)
        print(f"figure -> {figdir}/phase3b_crude_vs_ipw.png")
    except Exception as e:
        print(f"[fig] échec : {e}")
    print(f"sortie -> {out.name}")


def make_fig(prev_crude, prev_ipw, catR_drug, figdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    drugs = [d for d in RDRUGS if d in catR_drug]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
    for ax, drug in zip(axes.ravel(), drugs):
        labs = [l for l in MAJORS if l in prev_crude[drug]]
        x = np.arange(len(labs))
        ax.bar(x - 0.2, [prev_crude[drug][l] for l in labs], 0.4, label="brut", color="#c0392b")
        ax.bar(x + 0.2, [prev_ipw[drug][l] for l in labs], 0.4, label="IPW BioProject", color="#2980b9")
        ax.set_title(drug, fontsize=9)
        ax.set_xticks(x); ax.set_xticklabels(labs, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("prévalence R (%)", fontsize=7)
        if drug == drugs[0]:
            ax.legend(fontsize=7)
    fig.suptitle("Prévalence des déterminants catalogués-R : brut vs pondéré IPW (par BioProject)", fontsize=11)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(figdir / f"phase3b_crude_vs_ipw.{ext}", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    main()
