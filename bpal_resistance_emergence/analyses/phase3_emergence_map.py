#!/usr/bin/env python3
"""Phase 3 — Cartographie descriptive de l'émergence de la résistance BPaLM.

UN seul passage sur le pangénome (`pan_strains.pkl`, 132 609 souches) qui produit :
  1. Stratification par lignée : prévalence des déterminants catalogués-R par drogue et par
     lignée Coll, avec IC95 (Clopper-Pearson) et test de Fisher d'hétérogénéité inter-lignées.
  2. Émergence temporelle : porteurs de déterminants R par année + bilan pré/post-2022 (BPaLM).
  3. Phylogéographie : porteurs R par pays.
  4. Trichotomie par variant (catalogué-R / polymorphisme de lignée / candidat-nouveau), réutilise
     `résultats/phase1_feasibility.tsv`.

Garde-fous (cf. cahier) : (a) prévalence CRUE — le biais d'échantillonnage clinique (cliniques MDR
sur-représentées) N'est PAS corrigé par IPW ici (à faire en phase 3b) ; rapporter en « brut »,
ne comparer que lignées de curation comparable. (b) « variant dans gène DR » ≠ « variant de
résistance » : la prévalence n'utilise QUE les SPDI catalogués-R (parse du champ catalogue de
phase1). (c) Dénominateur = souches de la lignée présentes dans le pangénome.

Sorties : résultats/phase3_{lineage_prevalence,temporal,geo}.tsv + figures dans
article/figures/phase3_{lineage_heatmap,emergence_curve}.{png,pdf}.

Usage : python phase3_emergence_map.py [--limit N]
"""
import argparse, csv, glob, os, pickle, sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths

PREFIX = paths.REF_ACC + ":"
PLEN = len(PREFIX)
RDRUGS = ["bedaquiline", "clofazimine", "delamanid", "moxifloxacin", "linezolid", "levofloxacin"]


def lin_label(x):
    if x in (None, "", "NA"):
        return "NA"
    return x if x in ("BOV", "BOV_AFRI") else "L" + str(x)


def clopper_pearson(k, n, alpha=0.05):
    from scipy.stats import beta
    if n == 0:
        return (0.0, 0.0)
    lo = 0.0 if k == 0 else beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - alpha / 2, k + 1, n - k)
    return (lo, hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int)
    args = ap.parse_args()
    paths.ensure_dirs()
    figdir = paths.ARTICLE / "figures"; figdir.mkdir(parents=True, exist_ok=True)

    # --- 1. ensembles de SPDI depuis phase1_feasibility.tsv -------------------
    feas = paths.RESULTATS / "phase1_feasibility.tsv"
    catR_drug = defaultdict(set)        # drogue R -> {spdi catalogués-R}
    candidate = defaultdict(set)        # panel_group -> {spdi non-syn hors-catalogue}
    lineage_marker = []                 # variants à freq élevée = polymorphismes de lignée
    n_var = 0
    with open(feas) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n_var += 1
            spdi = r["spdi"]
            cat = r["catalogue"]
            for tok in cat.split(";"):
                if tok.endswith(":R-associated"):
                    catR_drug[tok.split(":")[0]].add(spdi)
            nonsyn = r["effect"] in ("missense", "stop_gained")
            if nonsyn and not cat:
                candidate[r["group"]].add(spdi)
            # polymorphisme de lignée = non-syn, fréquent, mais NON catalogué-R (ex. gyrA E21Q
            # est catalogué « uncertain », pas R -> il faut filtrer sur is_catalogued_R, pas sur 'cat')
            if nonsyn and r["is_catalogued_R"] == "0" and float(r["freq_pct"]) > 2.0:
                lineage_marker.append((r["gene"], r["aa_change"], float(r["freq_pct"]), r["n_lineages"]))
    allR = set().union(*catR_drug.values()) if catR_drug else set()
    print(f"phase1 : {n_var} variants | catalogués-R par drogue : "
          f"{ {d: len(s) for d, s in catR_drug.items()} }")
    print(f"candidats non-syn hors-catalogue par groupe : { {g: len(s) for g, s in candidate.items()} }")

    # --- 2. métadonnées per-souche -------------------------------------------
    strain2lin = {}
    with open(paths.LINEAGE_SNAPSHOT) as fh:
        for row in csv.DictReader(fh):
            strain2lin[row["strain_name"]] = lin_label(row["lineage_level_1"])
    print(f"snapshot lignée : {len(strain2lin)} souches")

    strain2geo = {}                     # SRA -> (country, year)
    for f in glob.glob(str(paths.BIOPROJECT_GEO / "consolidated_geo_*.tsv")):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                sra = row.get("strain") or row.get("strain_id")
                if not sra:
                    continue
                yr = (row.get("date_min") or "")[:4]
                yr = int(yr) if yr.isdigit() else None
                strain2geo[sra] = (row.get("country") or "", yr)
    print(f"géo+date consolidées : {len(strain2geo)} souches")

    # --- 3. passage pangénome -------------------------------------------------
    print("chargement pan_strains.pkl…")
    with open(paths.PAN_STRAINS, "rb") as f:
        pan = pickle.load(f)

    lin_tot = Counter()                                  # lignée -> n souches
    lin_R = defaultdict(Counter)                         # lignée -> {drogue: n porteurs R}
    lin_cand = defaultdict(Counter)                      # lignée -> {groupe: n porteurs candidat}
    year_R = defaultdict(Counter)                        # drogue -> {année: n porteurs R}
    geo_R = defaultdict(Counter)                         # drogue -> {pays: n porteurs R}
    prepost = defaultdict(lambda: Counter())             # drogue -> {'<=2021','>=2022': n}
    N = 0
    for i, (strain, spdis) in enumerate(pan.items()):
        if args.limit and i >= args.limit:
            break
        N += 1
        lab = strain2lin.get(strain, "NA")
        lin_tot[lab] += 1
        country, year = strain2geo.get(strain, ("", None))
        carried = spdis if isinstance(spdis, set) else set(spdis)
        for drug, sset in catR_drug.items():
            if carried & sset:
                lin_R[lab][drug] += 1
                if drug in RDRUGS and year:
                    year_R[drug][year] += 1
                    prepost[drug]["<=2021" if year <= 2021 else ">=2022"] += 1
                if drug in RDRUGS and country:
                    geo_R[drug][country] += 1
        for grp, sset in candidate.items():
            if carried & sset:
                lin_cand[lab][grp] += 1
        if (i + 1) % 30000 == 0:
            print(f"  …{i+1} souches")
    print(f"pangénome scanné : {N} souches")

    # --- 4. sorties stratification par lignée + Fisher -----------------------
    from scipy.stats import fisher_exact
    majors = ["L1", "L2", "L3", "L4", "L5", "L6", "L7", "BOV", "BOV_AFRI"]
    out = paths.RESULTATS / "phase3_lineage_prevalence.tsv"
    with open(out, "w") as fh:
        fh.write("drug\tlineage\tn_total\tn_R\tprevalence_pct\tci95_low\tci95_high\n")
        for drug in sorted(catR_drug):
            for lab in majors:
                nt = lin_tot.get(lab, 0); nr = lin_R[lab].get(drug, 0)
                if nt == 0:
                    continue
                lo, hi = clopper_pearson(nr, nt)
                fh.write(f"{drug}\t{lab}\t{nt}\t{nr}\t{100*nr/nt:.4f}\t{100*lo:.4f}\t{100*hi:.4f}\n")
    # Fisher 2x2 L2 vs reste (Beijing = drapeau MDR connu), par drogue
    print("\n=== Prévalence catalogués-R par lignée (brut) + Fisher L2-vs-reste ===")
    for drug in sorted(catR_drug):
        tot = sum(lin_tot[l] for l in majors if l in lin_tot)
        rtot = sum(lin_R[l].get(drug, 0) for l in majors)
        l2n, l2r = lin_tot.get("L2", 0), lin_R["L2"].get(drug, 0)
        oth_n, oth_r = tot - l2n, rtot - l2r
        try:
            odr, p = fisher_exact([[l2r, l2n - l2r], [oth_r, oth_n - oth_r]])
        except Exception:
            odr, p = float("nan"), float("nan")
        print(f"  {drug:13}: total R {rtot}/{tot} ({100*rtot/max(tot,1):.3f}%) | "
              f"L2 {l2r}/{l2n} ({100*l2r/max(l2n,1):.3f}%) vs reste {oth_r}/{oth_n} "
              f"({100*oth_r/max(oth_n,1):.3f}%) | OR={odr:.2f} p={p:.2e}")

    # temporel
    tout = paths.RESULTATS / "phase3_temporal.tsv"
    with open(tout, "w") as fh:
        fh.write("drug\tyear\tn_R_carriers\n")
        for drug in sorted(year_R):
            for yr in sorted(year_R[drug]):
                if 1990 <= yr <= 2026:
                    fh.write(f"{drug}\t{yr}\t{year_R[drug][yr]}\n")
    print("\n=== Émergence temporelle (porteurs R, pré/post-2022 BPaLM) ===")
    for drug in sorted(prepost):
        pp = prepost[drug]
        print(f"  {drug:13}: <=2021 {pp.get('<=2021',0)} | >=2022 {pp.get('>=2022',0)}")

    # géo
    gout = paths.RESULTATS / "phase3_geo.tsv"
    with open(gout, "w") as fh:
        fh.write("drug\tcountry\tn_R_carriers\n")
        for drug in sorted(geo_R):
            for c, n in geo_R[drug].most_common():
                fh.write(f"{drug}\t{c}\t{n}\n")

    # trichotomie
    print("\n=== Trichotomie (réponse au garde-fou « gène DR ≠ résistance ») ===")
    print(f"  catalogués-R (SPDI distincts) : {len(allR)}")
    print(f"  candidats non-syn hors-catalogue : {sum(len(s) for s in candidate.values())}")
    print(f"  polymorphismes de lignée probables (non-syn freq>2%) : {len(lineage_marker)}")
    for g, aa, fp, nl in sorted(lineage_marker, key=lambda x: -x[2])[:8]:
        print(f"      {g} {aa}  freq {fp:.1f}%  ({nl} lignées)")

    # --- 5. figures -----------------------------------------------------------
    try:
        make_figures(lin_tot, lin_R, year_R, majors, catR_drug, figdir)
        print(f"\nfigures -> {figdir}")
    except Exception as e:
        print(f"[fig] échec (données écrites quand même) : {e}")

    print(f"\nsorties -> {out.name}, {tout.name}, {gout.name}")


def make_figures(lin_tot, lin_R, year_R, majors, catR_drug, figdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    drugs = [d for d in RDRUGS if d in catR_drug]
    labs = [l for l in majors if lin_tot.get(l, 0) > 0]
    M = np.zeros((len(labs), len(drugs)))
    for i, l in enumerate(labs):
        for j, d in enumerate(drugs):
            nt = lin_tot.get(l, 0)
            M[i, j] = 100 * lin_R[l].get(d, 0) / nt if nt else 0
    fig, ax = plt.subplots(figsize=(1.3 * len(drugs) + 2, 0.5 * len(labs) + 2))
    im = ax.imshow(M, cmap="OrRd", aspect="auto")
    ax.set_xticks(range(len(drugs))); ax.set_xticklabels(drugs, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs, fontsize=8)
    for i in range(len(labs)):
        for j in range(len(drugs)):
            ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=6,
                    color="black" if M[i, j] < M.max() * 0.6 else "white")
    ax.set_title("Prévalence (%) des déterminants catalogués-R par lignée (brut)", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.7, label="% porteurs")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(figdir / f"phase3_lineage_heatmap.{ext}", dpi=600)
    plt.close(fig)

    # courbe d'émergence (bédaquiline + clofazimine)
    fig, ax = plt.subplots(figsize=(7, 4))
    for d in [x for x in ("bedaquiline", "clofazimine", "moxifloxacin") if x in year_R]:
        ys = sorted(y for y in year_R[d] if 1995 <= y <= 2026)
        ax.plot(ys, [year_R[d][y] for y in ys], marker="o", ms=3, label=d)
    ax.axvline(2022, ls="--", color="grey", lw=1)
    ax.text(2022.1, ax.get_ylim()[1] * 0.9, "BPaLM (OMS 2022)", fontsize=7, color="grey")
    ax.set_xlabel("Année de collecte"); ax.set_ylabel("Porteurs de déterminants R (n)")
    ax.set_title("Émergence temporelle des déterminants catalogués-R", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(figdir / f"phase3_emergence_curve.{ext}", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    main()
