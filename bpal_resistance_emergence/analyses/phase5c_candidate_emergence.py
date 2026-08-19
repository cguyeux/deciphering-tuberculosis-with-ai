#!/usr/bin/env python3
"""Phase 5c — Profil d'émergence des candidats convergents F420 (temporel / géo / lignée).

Phase 3 a cartographié l'émergence (temps, pays, lignée) des SEULS déterminants catalogués-R.
Les candidats DÉCOUVERTS en Phase 5 (convergence phénotype-free + structure + ESM-1v) n'ont
jamais reçu ce profil. Ce script comble le trou, sur le pangénome complet (`pan_strains.pkl`),
pour la grappe de candidats F420 retenue :

  convergents (poche cofacteur, ESM-1v damaging) :
    - fgd1 210  (A210E/G/P/T/V)  — NOUVEAU, hors-catalogue, 1re couronne F420 (7.8 Å), enfoui
    - fgd1 10   (A10V/P/T)       — delamanid:uncertain, 1re couronne F420 (7.6 Å), enfoui
    - ddn  111  (A111T/V/G)      — convergent_recurrent, 2de couronne (13 Å), enfoui
  phénotype-associés (cas-témoins DLM, Phase 5b volet B) :
    - ddn Y65S   — contact DIRECT F420 (3.1 Å), delamanid:uncertain, OR 42
    - ddn G81S   — 1re couronne (7 Å), delamanid:uncertain
    - ddn G34R   — N-term, secondaire/faible
  ANCRES cataloguées DLM-R (référence interne « signature d'émergence ère-médicament ») :
    - ddn L49P (R-assoc, OR 456), ddn W88* , ddn W139*  (stop, R-assoc)

Pour chaque POSITION (on agrège les SPDI alternatifs d'une même position : la convergence est
POSITIONNELLE), on calcule sur l'ensemble des porteurs :
  - distribution temporelle (année de collecte) : min/médiane/max + bins <2014 / 2014-2021 / >=2022 ;
  - TEST DE RÉCENCE vs FOND : Mann-Whitney U (unilatéral « carriers plus récents que le fond »),
    le fond = années de TOUTES les souches datées du pangénome. C'EST LE GARDE-FOU central :
    les cohortes cliniques récentes dominent le dénominateur → « tout paraît récent » ; le test
    dit si les porteurs sont PLUS récents que cette tendance de fond, pas juste récents.
  - étalement lignée : n lignées majeures distinctes, n sous-lignées (lineage_code) distinctes —
    discrimine convergence multi-clade vs expansion clonale unique. Par lignée majeure : n porteurs,
    n sous-lignées, fenêtre d'années (clonal si 1 sous-lignée resserrée, convergent sinon) ;
  - géographie : pays des porteurs (top) ;
  - fond MDR : fraction des porteurs co-portant un déterminant RIF-R catalogué (lien Phase 6).

Garde-fous (cf. cahier N1/N2). (a) Date de COLLECTE ≠ date d'ACQUISITION de la mutation : une
collecte récente d'un variant standing ancien est possible ; l'antécédence formelle (le variant
précède l'ère du médicament) exige pastml sur arbre daté (différé). Ici = profil descriptif +
test de récence relatif. (b) N faibles (fgd1 210/10 ≈ 9 porteurs) → on RAPPORTE, on ne teste
puissamment que là où n le permet. (c) confond lignée/échantillonnage non corrigé par IPW ici.

Sorties : résultats/phase5c_candidate_emergence.tsv (1 ligne/position),
          résultats/phase5c_candidate_carriers.tsv (1 ligne/porteur),
          article/figures/phase5c_candidate_timeline.{png,pdf}.
"""
import csv
import glob
import pickle
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths

DLM_CLIN = 2014   # délamanide : autorisation/usage clinique (proxy ère F420 sous pression)
BPALM = 2022      # recommandation OMS BPaLM

# --- candidats : position -> {class, label, spdis:{spdi: alt_label}} ---------
CANDIDATES = {
    "fgd1:210": {
        "cls": "convergent_novel", "gene": "fgd1", "pos": 210,
        "note": "1re couronne F420 7.8A enfoui; ESM A210E -8.4; hors-catalogue",
        "spdis": {
            "NC_000962.3:491409:G:A": "A210T", "NC_000962.3:491409:G:C": "A210P",
            "NC_000962.3:491410:C:A": "A210E", "NC_000962.3:491410:C:G": "A210G",
            "NC_000962.3:491410:C:T": "A210V",
        }},
    "fgd1:10": {
        "cls": "convergent_dlm_uncertain", "gene": "fgd1", "pos": 10,
        "note": "1re couronne F420 7.6A enfoui; ESM A10P -8.2; delamanid:uncertain",
        "spdis": {
            "NC_000962.3:490810:C:T": "A10V", "NC_000962.3:490809:G:C": "A10P",
            "NC_000962.3:490809:G:A": "A10T",
        }},
    "ddn:111": {
        "cls": "convergent_recurrent", "gene": "ddn", "pos": 111,
        "note": "2de couronne 13A enfoui; ESM A111T -5.0; PMD convergent_recurrent",
        "spdis": {
            "NC_000962.3:3987173:G:A": "A111T", "NC_000962.3:3987174:C:T": "A111V",
            "NC_000962.3:3987174:C:G": "A111G",
        }},
    "ddn:65": {
        "cls": "phenotype_casecontrol", "gene": "ddn", "pos": 65,
        "note": "contact DIRECT F420 3.1A; cas-temoins DLM OR42; delamanid:uncertain",
        "spdis": {"NC_000962.3:3987036:A:C": "Y65S"}},
    "ddn:81": {
        "cls": "phenotype_casecontrol", "gene": "ddn", "pos": 81,
        "note": "1re couronne 7A; cas-temoins DLM (2R/0S); delamanid:uncertain",
        "spdis": {"NC_000962.3:3987083:G:A": "G81S"}},
    "ddn:34": {
        "cls": "phenotype_casecontrol_weak", "gene": "ddn", "pos": 34,
        "note": "N-term; cas-temoins OR127 mais peripherique/desordonne",
        "spdis": {"NC_000962.3:3986942:G:A": "G34R", "NC_000962.3:3986942:G:C": "G34R(b)"}},
    "ddn:49": {
        "cls": "catalogued_dlm_R_anchor", "gene": "ddn", "pos": 49,
        "note": "ANCRE delamanid:R-associated (OR456) — reference ere-medicament",
        "spdis": {"NC_000962.3:3986988:T:C": "L49P"}},
    "ddn:88": {
        "cls": "catalogued_dlm_R_anchor", "gene": "ddn", "pos": 88,
        "note": "ANCRE delamanid:R-associated — stop_gained W88*",
        "spdis": {"NC_000962.3:3987105:G:A": "W88*"}},
    "ddn:139": {
        "cls": "catalogued_dlm_R_anchor", "gene": "ddn", "pos": 139,
        "note": "ANCRE delamanid:R-associated — stop_gained W139*",
        "spdis": {"NC_000962.3:3987258:G:A": "W139*"}},
}


def lin_label(x):
    if x in (None, "", "NA"):
        return "NA"
    return x if x in ("BOV", "BOV_AFRI") else "L" + str(x)


def main():
    paths.ensure_dirs()
    figdir = paths.ARTICLE / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    # index spdi -> position-key
    spdi2pos = {}
    all_cand_spdis = set()
    for key, d in CANDIDATES.items():
        for spdi in d["spdis"]:
            spdi2pos[spdi] = key
            all_cand_spdis.add(spdi)

    # --- déterminants RIF-R catalogués (fond MDR), depuis le catalogue consolidé ---
    # rpoB n'est PAS dans le panel BPaL → phase1_feasibility n'a pas de token rifampicin ;
    # on lit donc directement le catalogue (drug_code RIF, call R-associated).
    rifR = set()
    with open(paths.CATALOGUE_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("drug_code") == "RIF" and r.get("call") == "R-associated":
                rifR.add(r["spdi"])
    print(f"déterminants RIF-R catalogués (fond MDR) : {len(rifR)} SPDI")

    # --- lignées : strain -> (major, sublineage code) -------------------------
    strain2lin = {}
    with open(paths.LINEAGE_SNAPSHOT) as fh:
        rd = csv.DictReader(fh)
        f = rd.fieldnames
        sk = "strain_name" if "strain_name" in f else f[0]
        codek = "lineage_code" if "lineage_code" in f else None
        majk = "lineage_level_1" if "lineage_level_1" in f else f[-1]
        for r in rd:
            strain2lin[r[sk]] = (lin_label(r[majk]), (r[codek] if codek else r[majk]) or "NA")
    print(f"snapshot lignée : {len(strain2lin)} souches")

    # --- géo + date par SRA ---------------------------------------------------
    strain2geo = {}
    for fpath in glob.glob(str(paths.BIOPROJECT_GEO / "consolidated_geo_*.tsv")):
        with open(fpath) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                sra = row.get("strain") or row.get("strain_id")
                if not sra:
                    continue
                yr = (row.get("date_min") or "")[:4]
                yr = int(yr) if yr.isdigit() else None
                strain2geo[sra] = (row.get("country") or "", yr)
    print(f"géo+date consolidées : {len(strain2geo)} souches")

    # --- UN passage sur pan_strains -------------------------------------------
    print("chargement pan_strains.pkl…", flush=True)
    with open(paths.PAN_STRAINS, "rb") as f:
        pan = pickle.load(f)
    print(f"pangénome : {len(pan)} souches")

    background_years = []                       # années de TOUTES les souches datées (fond)
    carriers = defaultdict(list)                # poskey -> [carrier dict]
    for strain, spdis in pan.items():
        s = spdis if isinstance(spdis, set) else set(spdis)
        country, year = strain2geo.get(strain, ("", None))
        if year is not None:
            background_years.append(year)
        hit = s & all_cand_spdis
        if not hit:
            continue
        maj, code = strain2lin.get(strain, ("NA", "NA"))
        is_rifr = bool(s & rifR)
        for spdi in hit:
            poskey = spdi2pos[spdi]
            carriers[poskey].append({
                "strain": strain, "spdi": spdi, "alt": CANDIDATES[poskey]["spdis"][spdi],
                "year": year, "country": country, "major": maj, "subcode": code,
                "rifR": int(is_rifr),
            })
    bg_n = len(background_years)
    bg_med = statistics.median(background_years) if bg_n else None
    print(f"fond temporel : {bg_n} souches datées (médiane {bg_med})")

    # --- agrégation par position ----------------------------------------------
    rows = []
    for key, d in CANDIDATES.items():
        cc = carriers.get(key, [])
        # déduplique au niveau souche (une souche porte 1 alt à cette position en pratique)
        by_strain = {}
        for c in cc:
            by_strain.setdefault(c["strain"], c)
        cc = list(by_strain.values())
        n = len(cc)
        years = sorted(c["year"] for c in cc if c["year"] is not None)
        majs = Counter(c["major"] for c in cc)
        subs = {c["subcode"] for c in cc}
        rifr_n = sum(c["rifR"] for c in cc)
        # bins temporels
        b_pre = sum(1 for y in years if y < DLM_CLIN)
        b_mid = sum(1 for y in years if DLM_CLIN <= y < BPALM)
        b_post = sum(1 for y in years if y >= BPALM)
        # test de récence vs fond
        u_p = None
        if len(years) >= 3 and bg_n:
            try:
                _, u_p = mannwhitneyu(years, background_years, alternative="greater")
            except Exception:
                u_p = None
        # étalement intra-lignée (clonal vs convergent) : par lignée majeure, n sous-lignées
        per_major = {}
        for m in majs:
            mc = [c for c in cc if c["major"] == m]
            msub = {c["subcode"] for c in mc}
            myr = sorted(c["year"] for c in mc if c["year"] is not None)
            per_major[m] = (len(mc), len(msub), (myr[0] if myr else None), (myr[-1] if myr else None))
        countries = Counter(c["country"] for c in cc if c["country"])

        rows.append({
            "position": key, "class": d["cls"], "gene": d["gene"], "pos": d["pos"],
            "n_carriers": n, "n_dated": len(years),
            "n_major_lineages": len(majs), "n_sublineages": len(subs),
            "major_lineages": ",".join(f"{m}:{c}" for m, c in majs.most_common()),
            "year_min": years[0] if years else "", "year_med": (statistics.median(years) if years else ""),
            "year_max": years[-1] if years else "",
            "n_pre2014": b_pre, "n_2014_2021": b_mid, "n_post2022": b_post,
            "recency_vs_bg_p": (f"{u_p:.2e}" if u_p is not None else "n/a"),
            "rifR_co_n": rifr_n, "rifR_co_pct": (f"{100*rifr_n/n:.0f}" if n else "0"),
            "top_countries": ",".join(f"{c}:{k}" for c, k in countries.most_common(4)),
            "per_major_subspread": " | ".join(
                f"{m}(n{a},sub{b},{lo}-{hi})" for m, (a, b, lo, hi) in per_major.items()),
            "note": d["note"],
        })

    # --- écriture résumé -------------------------------------------------------
    out = paths.RESULTATS / "phase5c_candidate_emergence.tsv"
    cols = ["position", "class", "gene", "pos", "n_carriers", "n_dated",
            "n_major_lineages", "n_sublineages", "major_lineages",
            "year_min", "year_med", "year_max", "n_pre2014", "n_2014_2021", "n_post2022",
            "recency_vs_bg_p", "rifR_co_n", "rifR_co_pct", "top_countries",
            "per_major_subspread", "note"]
    with open(out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    # --- écriture porteurs -----------------------------------------------------
    cout = paths.RESULTATS / "phase5c_candidate_carriers.tsv"
    with open(cout, "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["position", "class", "alt", "strain", "year", "country", "major", "subcode", "rifR"])
        for key, d in CANDIDATES.items():
            seen = set()
            for c in carriers.get(key, []):
                if c["strain"] in seen:
                    continue
                seen.add(c["strain"])
                w.writerow([key, d["cls"], c["alt"], c["strain"], c["year"] or "",
                            c["country"], c["major"], c["subcode"], c["rifR"]])

    # --- récap console ---------------------------------------------------------
    print("\n=== Profil d'émergence des candidats (fond médiane = "
          f"{bg_med}, {bg_n} souches datées) ===")
    print(f"{'position':11}{'class':26}{'n':>4}{'dat':>4}{'maj':>4}{'sub':>4}"
          f"{'med':>6}{'<14':>5}{'14-21':>6}{'>=22':>5}{'recency_p':>11}{'rifR%':>6}")
    for r in rows:
        print(f"{r['position']:11}{r['class']:26}{r['n_carriers']:4}{r['n_dated']:4}"
              f"{r['n_major_lineages']:4}{r['n_sublineages']:4}{str(r['year_med']):>6}"
              f"{r['n_pre2014']:5}{r['n_2014_2021']:6}{r['n_post2022']:5}"
              f"{r['recency_vs_bg_p']:>11}{r['rifR_co_pct']:>6}")
    print("\nLecture : 'recency_p' = Mann-Whitney unilatéral (porteurs PLUS récents que le fond) ; "
          "significatif (<0.05) = émergence ère-médicament au-delà du biais d'échantillonnage. "
          "Comparer les candidats aux ANCRES ddn L49P/W88*/W139* (signature DLM-R connue). "
          "'sub' >> 1 sur plusieurs 'maj' = convergence multi-clade (≠ expansion clonale). "
          "'rifR%' = co-occurrence fond MDR (lien Phase 6). Garde-fou : collecte ≠ acquisition ; "
          "antécédence formelle = pastml différé.")

    # --- figure ---------------------------------------------------------------
    try:
        make_figure(rows, carriers, background_years, figdir)
        print(f"\nfigure -> {figdir}/phase5c_candidate_timeline.png")
    except Exception as e:
        print(f"[fig] échec (TSV écrits) : {e}")

    print(f"\nsorties -> {out.name}, {cout.name}")


def make_figure(rows, carriers, background_years, figdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # ordre : convergents/phénotype d'abord, ancres en bas
    order = [r for r in rows if not r["class"].startswith("catalogued")] + \
            [r for r in rows if r["class"].startswith("catalogued")]
    order = [r for r in order if r["n_dated"] > 0]
    palette = {"L1": "#1b9e77", "L2": "#d95f02", "L3": "#7570b3", "L4": "#e7298a",
               "L5": "#66a61e", "L6": "#e6ab02", "L7": "#a6761d", "BOV": "#666666",
               "BOV_AFRI": "#999999", "NA": "#cccccc"}

    fig, (axb, ax) = plt.subplots(
        2, 1, figsize=(8.2, 0.46 * len(order) + 2.4), height_ratios=[1, max(len(order), 3)],
        gridspec_kw={"hspace": 0.06})

    # bandeau fond (densité des années d'échantillonnage)
    by = [y for y in background_years if 1995 <= y <= 2026]
    axb.hist(by, bins=range(1995, 2027), color="#dddddd", edgecolor="none")
    axb.axvline(DLM_CLIN, ls="--", color="#444", lw=0.8)
    axb.axvline(BPALM, ls="--", color="crimson", lw=0.9)
    axb.set_xlim(1995, 2026)
    axb.set_yticks([])
    axb.set_xticklabels([])
    axb.set_ylabel("fond\n(échant.)", fontsize=7, rotation=0, ha="right", va="center")
    axb.text(DLM_CLIN + 0.2, axb.get_ylim()[1] * 0.55, "DLM 2014", fontsize=6.5, color="#444")
    axb.text(BPALM + 0.2, axb.get_ylim()[1] * 0.55, "BPaLM 2022", fontsize=6.5, color="crimson")

    seen_lin = set()
    for i, r in enumerate(order):
        cc = {c["strain"]: c for c in carriers.get(r["position"], [])}.values()
        pts = [(c["year"], c["major"]) for c in cc if c["year"] is not None]
        y0 = len(order) - 1 - i
        for (yr, maj) in pts:
            jit = (hash((r["position"], yr, maj)) % 7 - 3) * 0.045
            axb_col = palette.get(maj, "#cccccc")
            ax.scatter(yr, y0 + jit, s=34, color=axb_col, edgecolor="white",
                       linewidth=0.4, zorder=3)
            seen_lin.add(maj)
        med = r["year_med"]
        if med != "":
            ax.scatter(med, y0, marker="|", s=260, color="black", zorder=4, linewidth=1.4)
    ax.axvline(DLM_CLIN, ls="--", color="#444", lw=0.8, zorder=1)
    ax.axvline(BPALM, ls="--", color="crimson", lw=0.9, zorder=1)
    ax.set_xlim(1995, 2026)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_yticks(range(len(order)))
    labels = []
    for i, r in enumerate(reversed(order)):
        tag = "★" if r["class"].startswith("catalogued") else ""
        sig = " *" if (r["recency_vs_bg_p"] not in ("n/a",) and
                       _pf(r["recency_vs_bg_p"]) is not None and _pf(r["recency_vs_bg_p"]) < 0.05) else ""
        labels.append(f"{tag}{r['position']} (n={r['n_carriers']}, {r['n_major_lineages']}maj){sig}")
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("Année de collecte du porteur", fontsize=9)
    fig.suptitle("Phase 5c — émergence des candidats F420 convergents vs ancres DLM-R (★)",
                 fontsize=9.5, y=1.005)
    handles = [plt.Line2D([0], [0], marker="o", ls="", color=palette[l],
                          markeredgecolor="white", label=l)
               for l in ["L1", "L2", "L3", "L4", "L5", "L6", "BOV"] if l in seen_lin]
    handles.append(plt.Line2D([0], [0], marker="|", ls="", color="black", label="médiane"))
    ax.legend(handles=handles, fontsize=6.8, ncol=8, loc="lower center",
              bbox_to_anchor=(0.5, -0.16), frameon=False)
    for ext in ("png", "pdf"):
        fig.savefig(figdir / f"phase5c_candidate_timeline.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


def _pf(s):
    try:
        return float(s)
    except Exception:
        return None


if __name__ == "__main__":
    main()
