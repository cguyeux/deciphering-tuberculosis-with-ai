#!/usr/bin/env python3
"""Phase 5f (réponse review R07) — test de récence STRATIFIÉ (anti-artefact géographique/lignée).

Le test de récence (Phase 5c) compare les années des porteurs au fond GLOBAL des souches datées. Risque
(reviewer R07) : un variant restreint à un pays/lignée récemment sur-séquencé paraîtrait « émergent » par
artefact d'échantillonnage non capté par le fond global. On refait donc le test avec un fond RESTREINT :
  (a) global (rappel Phase 5c) ;
  (b) lignée-apparié : fond = souches des lignées majeures présentes chez les porteurs ;
  (c) pays-apparié : fond = souches des pays présents chez les porteurs.
Conclusion robuste si les ancres cataloguées DLM-R (L49P/W139*) restent significativement récentes ET les
candidats convergents restent non-récents, sous les 3 fonds.

Sortie : résultats/phase5f_recency_stratified.tsv + récap console. (Pas de pan_strains : snapshot + géo +
porteurs phase5c.)
"""
import csv
import glob
import pickle
import sys
from collections import defaultdict
from pathlib import Path

from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths


def lin_label(x):
    if x in (None, "", "NA"):
        return "NA"
    return x if x in ("BOV", "BOV_AFRI") else "L" + str(x)


def mw(carrier_years, bg_years):
    if len(carrier_years) < 3 or len(bg_years) < 10:
        return None
    try:
        _, p = mannwhitneyu(carrier_years, bg_years, alternative="greater")
        return p
    except Exception:
        return None


def main():
    # --- porteurs par position (phase5c) ---
    carriers = defaultdict(list)     # position -> [(year, country, major)]
    with open(paths.RESULTATS / "phase5c_candidate_carriers.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            y = r.get("year", "")
            yr = int(y) if y.isdigit() else None
            carriers[r["position"]].append((yr, r.get("country", ""), r.get("major", "NA")))

    # --- lignées (strain -> major) ---
    strain2lin = {}
    with open(paths.LINEAGE_SNAPSHOT) as fh:
        rd = csv.DictReader(fh)
        sk, lk = rd.fieldnames[0], rd.fieldnames[-1]
        for r in rd:
            strain2lin[r[sk]] = lin_label(r[lk])

    # --- restreindre le fond aux souches DU PANGÉNOME (= même univers que phase5c, cohérence) ---
    print("chargement pan_strains.pkl (restriction du fond = univers phase5c)…", flush=True)
    pan_keys = set(pickle.load(open(paths.PAN_STRAINS, "rb")).keys())

    # --- fond : souches datées du pangénome (strain -> (year, major, country)) ---
    bg_year, bg_major, bg_country = [], [], []
    for f in glob.glob(str(paths.BIOPROJECT_GEO / "consolidated_geo_*.tsv")):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                s = r.get("strain") or r.get("strain_id")
                if not s or s not in pan_keys:
                    continue
                y = (r.get("date_min") or "")[:4]
                if not y.isdigit():
                    continue
                bg_year.append(int(y))
                bg_major.append(strain2lin.get(s, "NA"))
                bg_country.append(r.get("country") or "")
    n_bg = len(bg_year)
    print(f"fond daté : {n_bg} souches (médiane globale {sorted(bg_year)[n_bg//2]})")

    rows = []
    print(f"\n{'position':12}{'n_dat':>6}{'p_global':>11}{'p_lignée':>11}{'p_pays':>11}  fonds (n)")
    for pos, cc in carriers.items():
        years = [y for y, _, _ in cc if y is not None]
        majs = {m for _, _, m in cc if m not in ("", "NA")}
        ctrys = {c for _, c, _ in cc if c}
        bg_lin = [bg_year[i] for i in range(n_bg) if bg_major[i] in majs] if majs else []
        bg_ctry = [bg_year[i] for i in range(n_bg) if bg_country[i] in ctrys] if ctrys else []
        p_g = mw(years, bg_year)
        p_l = mw(years, bg_lin)
        p_c = mw(years, bg_ctry)
        fmt = lambda p: f"{p:.3g}" if p is not None else "n/a"
        print(f"{pos:12}{len(years):6}{fmt(p_g):>11}{fmt(p_l):>11}{fmt(p_c):>11}  "
              f"lin={len(bg_lin)} pays={len(bg_ctry)}")
        rows.append({"position": pos, "n_dated": len(years),
                     "p_recency_global": fmt(p_g), "p_recency_lineage_matched": fmt(p_l),
                     "p_recency_country_matched": fmt(p_c),
                     "n_bg_lineage": len(bg_lin), "n_bg_country": len(bg_ctry),
                     "carrier_lineages": ",".join(sorted(majs)), "carrier_countries": ",".join(sorted(ctrys))})

    with open(paths.RESULTATS / "phase5f_recency_stratified.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=["position", "n_dated", "p_recency_global",
                                           "p_recency_lineage_matched", "p_recency_country_matched",
                                           "n_bg_lineage", "n_bg_country", "carrier_lineages",
                                           "carrier_countries"], delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print("\nLecture : robuste si les ancres ddn:49/139 restent significatives (<0.05) ET les candidats "
          "convergents (fgd1:210/10, ddn:111) restent non-significatifs sous les fonds lignée/pays-appariés. "
          "Si une « récence » disparaît au fond apparié = c'était un artefact géo/lignée.")
    print(f"sortie -> {paths.RESULTATS}/phase5f_recency_stratified.tsv")


if __name__ == "__main__":
    main()
