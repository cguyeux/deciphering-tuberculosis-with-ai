#!/usr/bin/env python3
"""P1.10 -- Expression de Rv0810c dans un compendium transcriptomique in vivo/ex vivo.

QUESTION. MtbTnDB (P0.6) couvre 146 conditions mais uniquement en FITNESS
(deplétion de mutants), jamais en EXPRESSION -- un gene peut etre induit en
conditions d'infection sans phenotype de fitness (les deux mesures sont
independantes). Aucun axe d'EXPRESSION n'a ete mobilise sur ce projet, alors
qu'il pourrait eclairer directement la question de localisation ouverte en
§4 (Rv0810c detectee dans des macrophages infectes, Chande 2015 ; phosphorylee
en culture, Malakar 2023 ; DeepTMHMM la classe cytoplasmique).

GATE POSE PAR LA PISTE ELLE-MEME avant tout calcul : "l'accès à un compendium
transcriptomique MTBC interrogeable n'est pas clairement outillé dans cet
ecosysteme... risque reel de cout disproportionne pour construire l'acces
avant meme de poser la question. A instruire seulement si un compendium
interrogeable est identifie a faible cout, pas a construire specifiquement
pour cette question." Ce script commence donc APRES cette verification :
Rychel et al. 2022 (mSphere, PMID cf. references.bib -- "Machine Learning of
All Mycobacterium tuberculosis H37Rv RNA-seq Data") a deja construit et publie
un compendium de 647 profils / 231 conditions (dont infection macrophage
THP-1/BMDM/alveolaire murin, neutrophile murin, hypoxie, dormance,
reactivation, stress nitrate/redox, sputum), avec la matrice normalisee
(log-TPM) et les metadonnees telechargeables librement depuis le depot GitHub
associe (Reosu/modulome_mtb, data/processed_data/). Aucune construction n'a
ete necessaire : telechargement direct + lecture pandas, cout quasi nul,
exactement le type de ressource que la piste demandait sans le nommer.

METHODE. Extraction de la ligne Rv0810c + les 74 temoins appariés en longueur
de P2.2/P8 (meme fichier `data/p8_ns_s_sites_75_controls.csv` que partout
ailleurs dans ce dossier -- meme discipline, pour eviter de lire comme
"remarquable" une variation qui serait en realite banale pour un petit gene
de cette classe). Les 231 conditions sont groupees en trois blocs a partir de
la colonne `project` des metadonnees : BASE (croissance standard, temoin de
reference), INFECTION (THP-1, macrophages alveolaires/BMDM murins,
neutrophiles murins), STRESS_HOTE (hypoxie, dormance, reactivation, resus-
citation, nitrate, redox, sputum, cellules chargees en lipides, biofilm --
stress rencontres in vivo sans etre une infection cellulaire a proprement
parler). Pour chaque gene (Rv0810c + 74 temoins), delta = moyenne(INFECTION
UNION STRESS_HOTE) - moyenne(BASE) en log-TPM ; rang percentile de Rv0810c
dans la distribution des 75 deltas.

GARDE-FOU : un delta positif documenterait une induction plausible en
conditions d'infection, coherente avec Chande 2015, mais ne prouverait ni
secretion ni fonction -- seulement que le gene est transcrit davantage dans
ces conditions. Un delta nul ou negatif serait un negatif informatif de plus
(le contraire de ce qu'on attendrait d'un gene specifiquement induit en
infection), a mettre en balance avec le fait que Rv0810c est deja un gene
CONSTITUTIVEMENT contraint (P8) : un gene sous forte selection purificatrice
peut tres bien etre exprime de facon stable plutot qu'inductible.

Sortie : résultats/p1_10_transcriptomic_compendium.json
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "résultats" / "p1_10_transcriptomic_compendium.json"

TARGET = "Rv0810c"

INFECTION_PROJECTS = {"THP-1", "miceAV_macro", "miceBMDM", "miceNF"}
STRESS_HOTE_PROJECTS = {"hypoxia", "dormancy", "reactivation", "resus",
                         "nitrate", "redox", "sputum", "fat_cells", "biofilm"}
BASE_PROJECTS = {"base", "growth"}


def main():
    controls_df = pd.read_csv(DATA / "p8_ns_s_sites_75_controls.csv")
    panel = sorted(controls_df["locus_tag"].unique())
    assert TARGET in panel and len(panel) == 75, (TARGET in panel, len(panel))

    expr = pd.read_csv(DATA / "log_tpm_norm.csv", index_col=0)
    meta = pd.read_csv(DATA / "metadata_final.csv").rename(
        columns={"Unnamed: 0": "sample_id"}).set_index("sample_id")
    meta = meta.loc[meta.index.intersection(expr.columns)]

    present = [g for g in panel if g in expr.index]
    missing = sorted(set(panel) - set(present))

    infection_samples = meta.index[meta["project"].isin(INFECTION_PROJECTS)]
    stress_samples = meta.index[meta["project"].isin(STRESS_HOTE_PROJECTS)]
    host_samples = infection_samples.union(stress_samples)
    base_samples = meta.index[meta["project"].isin(BASE_PROJECTS)]

    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit("P1.10 -- Expression de Rv0810c dans le compendium RNA-seq H37Rv "
         "(Rychel et al. 2022)")
    emit("=" * 78)
    emit(f"Compendium : {expr.shape[1]} profils, {meta['project'].nunique()} projets "
         f"regroupes en {expr.shape[0]} genes (source : GitHub Reosu/modulome_mtb).")
    emit(f"Panel de temoins (P2.2/P8, 50-70 aa) : {len(present)}/{len(panel)} genes "
         f"retrouves dans le compendium ({len(missing)} absents : {missing}).")
    emit(f"Echantillons BASE (reference) : {len(base_samples)}")
    emit(f"Echantillons INFECTION (THP-1/macrophage-neutrophile murin) : {len(infection_samples)}")
    emit(f"Echantillons STRESS_HOTE (hypoxie/dormance/reactivation/nitrate/redox/"
         f"sputum/lipides/biofilm) : {len(stress_samples)}")
    emit(f"Union INFECTION + STRESS_HOTE : {len(host_samples)}")
    emit("")

    deltas = {}
    detail = {}
    for gene in present:
        row = expr.loc[gene]
        base_mean = row[base_samples].mean()
        host_mean = row[host_samples].mean()
        infection_mean = row[infection_samples].mean() if len(infection_samples) else float("nan")
        stress_mean = row[stress_samples].mean() if len(stress_samples) else float("nan")
        deltas[gene] = host_mean - base_mean
        detail[gene] = {
            "base_mean_log_tpm": round(float(base_mean), 4),
            "infection_mean_log_tpm": round(float(infection_mean), 4),
            "stress_hote_mean_log_tpm": round(float(stress_mean), 4),
            "host_union_mean_log_tpm": round(float(host_mean), 4),
            "delta_host_minus_base": round(float(host_mean - base_mean), 4),
            "variance_across_all_conditions": round(float(row.var()), 4),
        }

    deltas_series = pd.Series(deltas).sort_values()
    rank = deltas_series.rank(pct=True)[TARGET]
    n = len(deltas_series)

    emit("DISTRIBUTION DES DELTAS (hote - base), Rv0810c + 74 temoins")
    emit("-" * 78)
    emit(f"  Rv0810c : delta = {deltas[TARGET]:+.4f} log-TPM "
         f"(base={detail[TARGET]['base_mean_log_tpm']}, "
         f"infection={detail[TARGET]['infection_mean_log_tpm']}, "
         f"stress_hote={detail[TARGET]['stress_hote_mean_log_tpm']})")
    emit(f"  Rang percentile parmi les {n} genes (temoins + cible) : {rank*100:.1f}e percentile")
    emit(f"  Distribution des temoins : min={deltas_series.min():+.4f}, "
         f"mediane={deltas_series.median():+.4f}, max={deltas_series.max():+.4f}")
    emit("")

    top5 = deltas_series.tail(5)
    bot5 = deltas_series.head(5)
    emit("  5 genes du panel les PLUS induits (hote - base) :")
    for g, d in top5.items():
        emit(f"    {g:10s} {d:+.4f}" + ("  <-- Rv0810c" if g == TARGET else ""))
    emit("  5 genes du panel les PLUS reprimes (hote - base) :")
    for g, d in bot5.items():
        emit(f"    {g:10s} {d:+.4f}" + ("  <-- Rv0810c" if g == TARGET else ""))
    emit("")

    emit("=" * 78)
    emit("VERDICT")
    if rank >= 0.90:
        emit(f"  >>> Rv0810c figure parmi les genes les PLUS INDUITS de son panel de taille "
             f"({rank*100:.1f}e percentile) en conditions hote (infection + stress). Signal")
        emit("      positif net, coherent avec la detection en macrophage de Chande 2015 --")
        emit("      a instruire plus avant (candidat rare et fort).")
    elif rank <= 0.10:
        emit(f"  >>> Rv0810c figure parmi les genes les PLUS REPRIMES de son panel "
             f"({rank*100:.1f}e percentile) en conditions hote. Ceci va A L'ENCONTRE d'une")
        emit("      induction specifique en infection -- negatif informatif, a mettre en")
        emit("      balance avec la detection en macrophage de Chande 2015 (protéomique, pas")
        emit("      transcriptomique : les deux mesures peuvent diverger).")
    else:
        emit(f"  >>> Rv0810c ({rank*100:.1f}e percentile) n'est ni particulierement induit ni")
        emit("      particulierement reprime en conditions d'infection/stress hote,")
        emit("      relativement a son panel de taille -- EXPRESSION STABLE, cohérente avec")
        emit("      un gene sous forte contrainte purificatrice constitutive (P8, percentile")
        emit("      2,7 en NS/S) plutot qu'un gene specifiquement inductible. La detection en")
        emit("      macrophage de Chande 2015 (proteomique) n'est donc pas expliquee par une")
        emit("      induction transcriptionnelle en infection -- compatible avec une proteine")
        emit("      constitutivement presente, detectee quelle que soit la condition, plutot")
        emit("      qu'un gene de reponse a l'hote.")
    emit("")
    emit("  RESERVE : ce test mesure une INDUCTION relative moyenne sur des blocs de")
    emit("  conditions heterogenes (agregation de plusieurs etudes/souches/protocoles),")
    emit("  pas une cinetique fine par condition individuelle. Un signal localise a une")
    emit("  seule sous-condition pourrait etre dilue dans la moyenne du bloc.")

    results = {
        "compendium": "Rychel et al. 2022 (mSphere), GitHub Reosu/modulome_mtb",
        "n_profils": int(expr.shape[1]),
        "n_projets": int(meta["project"].nunique()),
        "panel_taille": len(panel),
        "panel_retrouve": len(present),
        "panel_absents": missing,
        "n_base": int(len(base_samples)),
        "n_infection": int(len(infection_samples)),
        "n_stress_hote": int(len(stress_samples)),
        "cible": TARGET,
        "rang_percentile_delta": round(float(rank), 4),
        "detail_par_gene": detail,
    }
    OUT.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n")
    emit("")
    emit(f"Ecrit : {OUT}")


if __name__ == "__main__":
    main()
