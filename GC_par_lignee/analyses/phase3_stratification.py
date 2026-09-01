#!/usr/bin/env python3
"""
Objet       : P3.2 -- lever (ou confirmer) le dernier confondant de l'ecart
              inter-lignees du rapport pertes/gains de paires G:C : le biais de
              couverture GC-dependant, corrole au kit, a la profondeur et au
              BioProject, eux-memes corroles a la lignee. Cinq sections :
              (A) structure de confusion lignee x BioProject x plateforme ;
              (B) rapport par BioProject A L'INTERIEUR de chaque lignee, pour
              comparer l'amplitude intra-lignee a l'amplitude inter-lignees ;
              (C) comparaison INTRA-STRATE dans les BioProjects qui contiennent
              plusieurs lignees (la seule stratification vraiment identifiante) ;
              (D) GLM binomial loss ~ lignee + covariables techniques, erreurs
              groupees par BioProject ; (E) correlations du rapport avec la
              profondeur, le GC des reads, le taux Q30 et la duplication.
Entrees     : resultats/phase3_counts_par_souche_n40.tsv (phase3_counts_par_souche.py)
              data/metadata_tbannotator_n40.csv  (mv_strain_metadata + tb_report_quality)
              data/depth_tbannotator_n40.csv     (tb_report_strain_spdi.depth agrege)
Sorties     : rapport texte sur stdout + TSV par section (--out-prefix)
Reutilisable: oui -- le geste "stratifier un signal genomique par BioProject et
              profondeur" vaut pour tout projet du depot sur base de convenance
Projet      : GC_par_lignee
Date        : 2026-08-29
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent


def ratio_ci(loss, gain, conf=0.95):
    """IC de log(loss/gain) par delta-methode sur des comptes de Poisson."""
    if loss == 0 or gain == 0:
        return float("nan"), float("nan")
    z = stats.norm.ppf(0.5 + conf / 2)
    lr = np.log(loss / gain)
    se = np.sqrt(1 / loss + 1 / gain)
    return float(np.exp(lr - z * se)), float(np.exp(lr + z * se))


def load(args):
    c = pd.read_csv(args.counts, sep="\t")
    m = pd.read_csv(args.meta)
    d = pd.read_csv(args.depth)
    df = c.merge(m, left_on="sra", right_on="strain_name", how="left")
    df = df.merge(d, on="strain_name", how="left")
    df["ratio"] = df["loss"] / df["gain"].replace(0, np.nan)
    df["cov_depth"] = df["med_depth"]
    df["gc_reads"] = df["after_gc_content"]
    df["bp"] = df["ncbi_bioproject"].fillna("NA")
    df["model"] = df["ncbi_model"].fillna("NA")
    return df


def section_a(df, out):
    rows = []
    for cl, g in df.groupby("clade"):
        n = len(g)
        bp = g["bp"].value_counts()
        loss, gain = g["loss"].sum(), g["gain"].sum()
        lo, hi = ratio_ci(loss, gain)
        rows.append(dict(
            clade=cl, n=n, ratio=loss / gain, ci_lo=lo, ci_hi=hi,
            n_bioprojects=g["bp"].nunique(),
            top_bp=bp.index[0], top_bp_frac=bp.iloc[0] / n,
            simpson=float((bp / n).pow(2).sum()),
            n_models=g["model"].nunique(),
            top_model=g["model"].value_counts().index[0],
            depth_med=g["cov_depth"].median(),
            gc_reads_med=g["gc_reads"].median(),
            q30_med=g["after_q30_rate"].median(),
            dup_med=g["duplication_rate"].median(),
            n_missing_meta=int(g["ncbi_bioproject"].isna().sum())))
    t = pd.DataFrame(rows).sort_values("ratio", ascending=False)
    t.to_csv(out + "A_structure_confusion.tsv", sep="\t", index=False)
    return t


def section_b(df, out, min_n=5):
    rows = []
    for (cl, bp), g in df.groupby(["clade", "bp"]):
        if len(g) < min_n:
            continue
        loss, gain = g["loss"].sum(), g["gain"].sum()
        if gain == 0:
            continue
        lo, hi = ratio_ci(loss, gain)
        rows.append(dict(clade=cl, bioproject=bp, n=len(g), loss=loss, gain=gain,
                         ratio=loss / gain, ci_lo=lo, ci_hi=hi,
                         depth_med=g["cov_depth"].median(),
                         gc_reads_med=g["gc_reads"].median(),
                         model=g["model"].value_counts().index[0]))
    t = pd.DataFrame(rows).sort_values(["clade", "ratio"])
    t.to_csv(out + "B_ratio_par_bioproject.tsv", sep="\t", index=False)
    return t


def section_c(df, out, min_n=5):
    """BioProjects contenant au moins deux lignees : la seule strate ou l'effet
    lignee est identifiable independamment du protocole de sequencage."""
    rows = []
    for bp, g in df.groupby("bp"):
        if bp == "NA":
            continue
        sub = g.groupby("clade").filter(lambda x: len(x) >= min_n)
        if sub["clade"].nunique() < 2:
            continue
        for cl, gg in sub.groupby("clade"):
            loss, gain = gg["loss"].sum(), gg["gain"].sum()
            lo, hi = ratio_ci(loss, gain)
            rows.append(dict(bioproject=bp, clade=cl, n=len(gg), loss=loss,
                             gain=gain, ratio=loss / gain, ci_lo=lo, ci_hi=hi,
                             depth_med=gg["cov_depth"].median(),
                             gc_reads_med=gg["gc_reads"].median()))
    t = pd.DataFrame(rows)
    if not t.empty:
        t = t.sort_values(["bioproject", "ratio"], ascending=[True, False])
    t.to_csv(out + "C_intra_bioproject.tsv", sep="\t", index=False)
    return t


def glm(df, formula_terms, ref_clade="L2.2.1", cluster=True):
    """GLM binomial : succes = perte de paire G:C parmi (perte + gain), erreurs
    groupees par BioProject (une strate technique = un cluster)."""
    d = df.dropna(subset=["loss", "gain"]).copy()
    d = d[(d["loss"] + d["gain"]) > 0]
    X = pd.get_dummies(d["clade"], prefix="clade", drop_first=False, dtype=float)
    X = X.drop(columns=[f"clade_{ref_clade}"])
    for t in formula_terms:
        v = d[t].astype(float)
        X[t] = (v - v.mean()) / v.std()
    X = sm.add_constant(X, has_constant="add")
    keep = X.notna().all(axis=1)
    X, d = X[keep], d[keep]
    endog = np.column_stack([d["loss"].values, d["gain"].values])
    m = sm.GLM(endog, X.values, family=sm.families.Binomial())
    if cluster:
        r = m.fit(cov_type="cluster",
                  cov_kwds={"groups": pd.factorize(d["bp"])[0]})
    else:
        r = m.fit()
    return r, X.columns.tolist(), d


def section_d(df, out, ref="L2.2.1"):
    lines = []
    tabs = {}
    for name, terms in [("brut", []),
                        ("ajuste", ["cov_depth", "gc_reads", "after_q30_rate",
                                    "duplication_rate"])]:
        r, cols, d = glm(df, terms, ref_clade=ref)
        t = pd.DataFrame(dict(term=cols, coef=r.params, se=r.bse,
                              z=r.tvalues, p=r.pvalues))
        t["odds_ratio"] = np.exp(t["coef"])
        t["model"] = name
        t["n_strains"] = len(d)
        tabs[name] = t
        lines.append(f"[{name}] n={len(d)} souches, "
                     f"df_resid={r.df_resid}, "
                     f"pearson_chi2/df={r.pearson_chi2 / r.df_resid:.2f}")
    t = pd.concat(tabs.values())
    t.to_csv(out + "D_glm.tsv", sep="\t", index=False)
    return t, lines


def section_e(df, out):
    rows = []
    d = df.dropna(subset=["ratio"])
    for v in ["cov_depth", "gc_reads", "after_q30_rate", "duplication_rate",
              "n_var", "insert_size_peak"]:
        s = d.dropna(subset=[v])
        if len(s) < 20:
            continue
        rho, p = stats.spearmanr(s[v], s["ratio"])
        # intra-lignee : correlation partielle par rangs centres par clade
        z = s.copy()
        z["rv"] = z.groupby("clade")[v].rank()
        z["rr"] = z.groupby("clade")["ratio"].rank()
        z["rv"] -= z.groupby("clade")["rv"].transform("mean")
        z["rr"] -= z.groupby("clade")["rr"].transform("mean")
        rho_w, p_w = stats.pearsonr(z["rv"], z["rr"])
        rows.append(dict(variable=v, n=len(s), spearman_global=rho, p_global=p,
                         spearman_intra_lignee=rho_w, p_intra_lignee=p_w))
    t = pd.DataFrame(rows)
    t.to_csv(out + "E_correlations.tsv", sep="\t", index=False)
    return t



def section_f(df, out, min_n=5, min_clades=3):
    """Comparaison INTRA-PLATEFORME. Contrairement au BioProject (presque
    toujours mono-lignee), le modele de sequenceur est partage entre lignees :
    l'effet lignee y est identifiable a protocole constant."""
    rows = []
    for mod, g in df.groupby("model"):
        if mod == "NA":
            continue
        sub = g.groupby("clade").filter(lambda x: len(x) >= min_n)
        if sub["clade"].nunique() < min_clades:
            continue
        for cl, gg in sub.groupby("clade"):
            loss, gain = gg["loss"].sum(), gg["gain"].sum()
            if gain == 0:
                continue
            lo, hi = ratio_ci(loss, gain)
            rows.append(dict(strate=f"model:{mod}", clade=cl, n=len(gg),
                             loss=loss, gain=gain, ratio=loss / gain,
                             ci_lo=lo, ci_hi=hi,
                             depth_med=gg["cov_depth"].median()))
    # strates de profondeur (terciles globaux), a modele libre
    d = df.dropna(subset=["cov_depth"]).copy()
    d["tercile"] = pd.qcut(d["cov_depth"], 3, labels=["Q1", "Q2", "Q3"])
    for terc, g in d.groupby("tercile", observed=True):
        sub = g.groupby("clade").filter(lambda x: len(x) >= min_n)
        for cl, gg in sub.groupby("clade"):
            loss, gain = gg["loss"].sum(), gg["gain"].sum()
            if gain == 0:
                continue
            lo, hi = ratio_ci(loss, gain)
            rows.append(dict(strate=f"depth:{terc}", clade=cl, n=len(gg),
                             loss=loss, gain=gain, ratio=loss / gain,
                             ci_lo=lo, ci_hi=hi,
                             depth_med=gg["cov_depth"].median()))
    t = pd.DataFrame(rows)
    t.to_csv(out + "F_intra_strate_technique.tsv", sep="\t", index=False)
    return t


def section_g(df, out, ref="L2.2.1"):
    """GLM avec effets fixes de plateforme et de tercile de profondeur : le
    classement des lignees est-il le meme a protocole tenu constant ?"""
    d = df.dropna(subset=["cov_depth", "model"]).copy()
    d = d[(d["loss"] + d["gain"]) > 0]
    d["tercile"] = pd.qcut(d["cov_depth"], 3, labels=["Q1", "Q2", "Q3"])
    X = pd.get_dummies(d["clade"], prefix="clade", dtype=float)
    X = X.drop(columns=[f"clade_{ref}"])
    Xm = pd.get_dummies(d["model"], prefix="mod", drop_first=True, dtype=float)
    Xt = pd.get_dummies(d["tercile"], prefix="dep", drop_first=True, dtype=float)
    X = pd.concat([X, Xm, Xt], axis=1)
    X = sm.add_constant(X, has_constant="add")
    endog = np.column_stack([d["loss"].values, d["gain"].values])
    r = sm.GLM(endog, X.values, family=sm.families.Binomial()).fit(
        cov_type="cluster", cov_kwds={"groups": pd.factorize(d["bp"])[0]})
    cols = X.columns.tolist()
    t = pd.DataFrame(dict(term=cols, coef=r.params, se=r.bse, z=r.tvalues,
                          p=r.pvalues, odds_ratio=np.exp(r.params)))
    t = t[t["term"].str.startswith(("clade_", "const"))]
    idx = [i for i, c in enumerate(cols) if c.startswith("clade_")]
    R = np.zeros((len(idx), len(cols)))
    for j, i in enumerate(idx):
        R[j, i] = 1
    w = r.wald_test(R, scalar=True)
    t.to_csv(out + "G_glm_plateforme.tsv", sep="\t", index=False)
    return t, w, len(d), Xm.shape[1] + 1


def section_h(df, out, min_n=5):
    """Unite d'observation = le BioProject, pas la souche. Test le plus
    conservateur possible sur la non-independance technique."""
    rows = []
    for (cl, bp), g in df.groupby(["clade", "bp"]):
        if len(g) < min_n or bp == "NA":
            continue
        loss, gain = g["loss"].sum(), g["gain"].sum()
        if gain < 5:
            continue
        rows.append(dict(clade=cl, bioproject=bp, n=len(g),
                         ratio=loss / gain))
    t = pd.DataFrame(rows)
    t.to_csv(out + "H_unite_bioproject.tsv", sep="\t", index=False)
    groups = [g["ratio"].values for _, g in t.groupby("clade")
              if len(g) >= 1]
    kw = stats.kruskal(*groups) if len(groups) >= 2 else None
    return t, kw



def section_i(df, out):
    """Decomposition de variance : la variabilite du rapport entre souches
    tient-elle davantage a la lignee ou au BioProject ? Modele mixte sur le
    logit du taux de perte par souche, lignee en effet fixe, BioProject en
    intercept aleatoire, poids = nombre de sites informatifs."""
    import statsmodels.formula.api as smf
    d = df.dropna(subset=["loss", "gain"]).copy()
    d = d[(d["loss"] + d["gain"]) >= 20]
    d = d[d["bp"] != "NA"]
    d["y"] = np.log((d["loss"] + 0.5) / (d["gain"] + 0.5))
    d["w"] = d["loss"] + d["gain"]
    m = smf.mixedlm("y ~ C(clade)", d, groups=d["bp"]).fit(reml=True)
    var_bp = float(m.cov_re.iloc[0, 0])
    var_res = float(m.scale)
    # variance inter-lignees : dispersion des moyennes ajustees, ponderee
    mu = d.groupby("clade").apply(lambda g: np.log(g["loss"].sum() /
                                                   g["gain"].sum()),
                                  include_groups=False)
    w = d.groupby("clade")["w"].sum()
    var_lin = float(np.average((mu - np.average(mu, weights=w)) ** 2, weights=w))
    t = pd.DataFrame([dict(composante="inter-lignees (log-ratio)", variance=var_lin,
                           ecart_type=np.sqrt(var_lin)),
                      dict(composante="inter-BioProjects intra-lignee",
                           variance=var_bp, ecart_type=np.sqrt(var_bp)),
                      dict(composante="residuelle inter-souches",
                           variance=var_res, ecart_type=np.sqrt(var_res))])
    t["part"] = t["variance"] / t["variance"].sum()
    t.to_csv(out + "I_composantes_variance.tsv", sep="\t", index=False)
    return t, len(d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default=str(ROOT / "résultats" /
                                            "phase3_counts_par_souche_n40.tsv"))
    ap.add_argument("--meta", default=str(ROOT / "data" /
                                          "metadata_tbannotator_n40.csv"))
    ap.add_argument("--depth", default=str(ROOT / "data" /
                                           "depth_tbannotator_n40.csv"))
    ap.add_argument("--out-prefix", default=str(ROOT / "résultats" / "phase3_"))
    args = ap.parse_args()

    pd.set_option("display.width", 200, "display.max_columns", 40,
                  "display.float_format", lambda x: f"{x:.4g}")
    df = load(args)
    out = args.out_prefix

    print("=" * 78)
    print(f"P3.2 -- stratification. {len(df)} souches, "
          f"{df['ncbi_bioproject'].notna().sum()} avec metadonnees NCBI, "
          f"{df['cov_depth'].notna().sum()} avec profondeur.")
    print("=" * 78)

    print("\n--- A. Structure de confusion lignee x BioProject x plateforme ---")
    print(section_a(df, out).to_string(index=False))

    print("\n--- B. Rapport par BioProject A L'INTERIEUR de chaque lignee "
          "(BioProjects a n>=5) ---")
    b = section_b(df, out)
    print(b.to_string(index=False))
    if not b.empty:
        amp = b.groupby("clade")["ratio"].agg(["min", "max", "count"])
        amp["amplitude"] = amp["max"] / amp["min"]
        print("\n  amplitude intra-lignee entre BioProjects :")
        print(amp[amp["count"] >= 2].to_string())

    print("\n--- C. Comparaison INTRA-STRATE : BioProjects multi-lignees ---")
    c = section_c(df, out)
    print(c.to_string(index=False) if not c.empty else
          "  aucun BioProject ne contient deux lignees a n>=5")

    print("\n--- D. GLM binomial (succes = perte de paire G:C), "
          "erreurs groupees par BioProject ---")
    dtab, lines = section_d(df, out)
    for l in lines:
        print("  " + l)
    print(dtab.to_string(index=False))

    print("\n--- E. Correlations du rapport avec les covariables techniques ---")
    print(section_e(df, out).to_string(index=False))

    print("\n--- F. Comparaison INTRA-STRATE TECHNIQUE (plateforme, profondeur) ---")
    f = section_f(df, out)
    print(f.to_string(index=False))
    glob = df.groupby("clade").apply(
        lambda g: g["loss"].sum() / g["gain"].sum(), include_groups=False)
    print("\n  conservation du classement global dans chaque strate "
          "(Spearman sur les lignees communes) :")
    for st, g in f.groupby("strate"):
        common = [c for c in g["clade"] if c in glob.index]
        if len(common) < 4:
            continue
        a = g.set_index("clade").loc[common, "ratio"]
        b = glob.loc[common]
        rho, p = stats.spearmanr(a, b)
        print(f"    {st:28s} k={len(common)}  rho={rho:+.3f}  p={p:.2g}")

    print("\n--- G. GLM a effets fixes de plateforme et de profondeur ---")
    gt, w, ng, npar = section_g(df, out)
    print(f"  n={ng} souches, {npar} indicatrices techniques ; "
          f"test de Wald joint sur les 10 coefficients de lignee : "
          f"chi2={float(w.statistic):.1f}, ddl={int(w.df_denom) if w.df_denom else 10}, "
          f"p={float(w.pvalue):.3g}")
    print(gt.to_string(index=False))

    print("\n--- H. Unite d'observation = le BioProject ---")
    ht, kw = section_h(df, out)
    print(ht.to_string(index=False))
    if kw:
        print(f"  Kruskal-Wallis sur les rapports par BioProject : "
              f"H={kw.statistic:.2f}, p={kw.pvalue:.3g}")
    ht3, kw3 = section_h(df, out.replace("phase3_", "phase3_minN3_"), min_n=3)
    if kw3:
        print(f"  sensibilite (BioProjects a n>=3, {len(ht3)} strates sur "
              f"{ht3['clade'].nunique()} lignees) : H={kw3.statistic:.2f}, "
              f"p={kw3.pvalue:.3g}")

    print("\n--- I. Composantes de variance du log-rapport ---")
    it, ni = section_i(df, out)
    print(f"  n={ni} souches (>=20 sites informatifs, BioProject connu)")
    print(it.to_string(index=False))
    print(f"\n# TSV ecrits sous {out}[A-E]_*.tsv")


if __name__ == "__main__":
    main()
