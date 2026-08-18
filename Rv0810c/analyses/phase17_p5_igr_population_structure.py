#!/usr/bin/env python3
"""
P5 -- IGR purM-Rv0810c (904820-904904, 85 nt) : modele nul obligatoire avant
tout investissement (cf. piste P5). Deux sources reinstruites :

(1) Zeng et al. 2018 (BMC Genomics, PMCID PMC5956929) -- fichiers supplementaires
    recuperes via l'API Europe PMC (supplementaryFiles), PAS relus a la main sur
    une page web (garde-fou CLAUDE.md sur les tables de supplement mal alignees,
    deja rencontre trois fois sur ce projet : Zeng Table S11, Malakar Data Set S8).
    MOESM2 (.docx, 22 tables natives) contient les tables region-marker.

(2) TBannotator (PostgreSQL, mv_spdi_mutations x mv_strain_metadata x
    mv_strain_classification) -- tous les variants SPDI observes dans l'IGR sur
    ~255000 souches, croises avec antibiogramme RIF et lignee.

Verdict recherche : le "hit GWAS RIF/RFB" (CRyPTIC 2022, Table 1, rang 13/20 RIF
et 4/20 RFB) survit-il a un controle direct sur les variants reellement presents
dans l'IGR, et le signal de resistance de Zeng 2018 concerne-t-il vraiment RIF ?
"""
import json
from pathlib import Path

OUT_JSON = Path(__file__).resolve().parent.parent / "résultats" / "p5_igr_population_structure.json"

# --- (1) Zeng et al. 2018, MOESM2, tables natives (extraites par python-docx) ---
# Table 6 : purM-Rv0810c est bien l'un des 20 "region markers" (IGR), confirme.
# Table 3 (Makers x Chi2|Fisher p-value par medicament) :
ZENG_TABLE3_PVALUES = {
    "INH": "0.47|1", "RIF": "0.38|1", "PZA": "--", "STR": "--", "EMB": "0.64|0.62",
    "OFX": "0.24|0.38", "MOX": "--", "ETH": "0.00|0.02", "KAN": "--", "AMI": "--",
    "CAP": "--", "PRO": "--",
}
# Table 14 (p-value ajustee finale, un seul medicament survit par marqueur) :
ZENG_TABLE14_ADJUSTED_P = {"ETH": 3.25e-03}  # tous les autres medicaments = NaN
# Table 20 (frequence du marqueur dans 3 cohortes de validation independantes) :
ZENG_TABLE20_VALIDATION = {
    "Casali_resistant": None, "Casali_sensitive": 0.0,
    "Zhang_resistant": 0.056, "Zhang_sensitive": 0.0,
    "Farhat_resistant": 0.0, "Farhat_sensitive": 0.0,
}

# --- (2) TBannotator, requetes verifiees directement (mcp__tbannotator) ---
# IGR purM-Rv0810c = 904820-904904 (purM/Rv0809 904725-904819 brin + ; Rv0810c 904905-905087 brin -)
IGR_TOP10_VARIANTS = [
    # (spdi, mutation, strain_count, rif_r, rif_s, rif_na)
    ("NC_000962.3:904884:CG:C", "del_CG", 526, 0, 33, 493),
    ("NC_000962.3:904821:AC:A", "del_AC", 126, 0, 15, 111),
    ("NC_000962.3:904884:C:CG", "ins_CG", 59, 0, 2, 57),
    ("NC_000962.3:904885:G:A", "G>A", 23, 0, 3, 20),
    ("NC_000962.3:904855:A:G", "A>G", 17, 0, 0, 17),
    ("NC_000962.3:904891:C:T", "C>T", 16, 0, 0, 16),
    ("NC_000962.3:904889:G:T", "G>T", 14, 0, 0, 14),
    ("NC_000962.3:904893:G:C", "G>C", 12, 0, 0, 12),
    ("NC_000962.3:904896:C:T", "C>T", 11, 0, 3, 8),
    ("NC_000962.3:904822:C:A", "C>A", 10, 0, 3, 7),
]
LINEAGE_PURITY = {
    "NC_000962.3:904821:AC:A (del_AC, 126 souches)": {"lineage_level_1": {"4": 126}, "purity": 1.0,
        "note": "marqueur EXCLUSIF de L4 -- signature de confusion par lignee, pas de causalite"},
    "NC_000962.3:904885:G:A (23 souches)": {"lineage_level_1": {"1": 23}, "purity": 1.0,
        "note": "marqueur EXCLUSIF de L1"},
    "NC_000962.3:904884:CG:C (del_CG, 526 souches)": {
        "lineage_level_1": {"1": 340, "2": 152, "4": 22, "Caprae_La2+La4": 6, "Caprae_La2": 6,
                             "3": 4, "Bovis_sans_La1": 2, "Bovis_La1": 2},
        "purity": 340 / 534, "note": "polyphyletique (8 lignees), pas un marqueur de lignee unique"},
}
RIF_BASELINE = {"rif_r_total": 388, "rif_s_total": 23035}  # mv_strain_metadata, species_group='M. tuberculosis'


def fisher_pooled():
    from scipy.stats import fisher_exact
    igr_rif_r = sum(v[3] for v in IGR_TOP10_VARIANTS)
    igr_rif_s = sum(v[4] for v in IGR_TOP10_VARIANTS)
    table = [[igr_rif_r, igr_rif_s], [RIF_BASELINE["rif_r_total"], RIF_BASELINE["rif_s_total"]]]
    odds, p_two = fisher_exact(table, alternative="two-sided")
    _, p_greater = fisher_exact(table, alternative="greater")
    return {
        "igr_rif_r_pooled_top10": igr_rif_r,
        "igr_rif_s_pooled_top10": igr_rif_s,
        "contingency_table_[[igr_R,igr_S],[baseline_R,baseline_S]]": table,
        "odds_ratio": odds,
        "fisher_p_two_sided": p_two,
        "fisher_p_one_sided_enrichment": p_greater,
    }


def main():
    fisher = fisher_pooled()

    verdict = (
        "MODELE NUL : L'ASSOCIATION RIF/RFB NE SURVIT PAS AU CONTROLE. Trois constats independants "
        "convergent. (1) Zeng 2018 lui-meme, relu sur les tables natives du supplement (pas un resume "
        "web) : la SEULE association medicamenteuse significative de purM-Rv0810c est l'ETHIONAMIDE "
        "(p_ajuste=3.25e-3, table 14 ; chi2/Fisher 0.00/0.02, table 3) -- RIF affiche p=0.38/1, non "
        "significatif. La 'convergence a trois signaux' de l'enonce de piste etait donc en partie une "
        "confusion de medicament : Zeng documente une resistance a l'ETH, pas au RIF/RFB. Validation "
        "croisee faible : frequence du marqueur proche de 0 dans 2 des 3 cohortes independantes "
        "(Casali, Farhat), 5.6% seulement chez Zhang_resistant. (2) TBannotator (recroise "
        "independamment, ~255000 souches) : ZERO souche resistante a la rifampicine parmi TOUTES les "
        "porteuses des 10 variants les plus frequents de l'IGR (0/59 phenotypees), alors que 388/23423 "
        "souches (1,66%) le sont dans la base entiere -- Fisher exact one-sided p=1.0, OR=0, aucune "
        "trace d'enrichissement. (3) Les deux variants les plus discriminants sont des MARQUEURS DE "
        "LIGNEE PURS : 904821:AC:A (126 souches) est EXCLUSIF a L4 (126/126), 904885:G:A (23 souches) "
        "est EXCLUSIF a L1 (23/23) -- la signature exacte de confusion par structure de population que "
        "le contre-argument de la piste anticipait, sur un evenement indel qui cree naturellement un "
        "k-mer/oligonucleotide propre a un clade, exactement le type de variant que la methode "
        "GWAS de CRyPTIC (test oligonucleotide/oligopeptide) est susceptible de capturer sans que cela "
        "implique une causalite. NUANCE A GARDER : CRyPTIC 2022 utilise un modele lineaire mixte "
        "(LMM/GEMMA) cense corriger la parente ; ce controle standard n'elimine pas totalement la "
        "confusion residuelle quand la resistance elle-meme est distribuee de façon clonale (limite "
        "documentee de la methode, pas une erreur des auteurs). Le variant CRyPTIC exact n'a pas ete "
        "retrouve (Table 1 ne donne qu'un rang gene-niveau, pas de position) : cette verification porte "
        "sur les variants IGR reellement observes et les plus frequents, donc les plus capables de "
        "porter un signal a l'echelle de 10228 genomes -- pas sur le variant exact de CRyPTIC, "
        "non identifie precisement malgre la recherche."
    )

    summary = {
        "igr_coordinates": "NC_000962.3:904820-904904 (purM/Rv0809 903725-904819 brin +, Rv0810c 904905-905087 brin -)",
        "zeng2018_source": "Europe PMC supplementaryFiles API, PMCID PMC5956929, MOESM2 (.docx, 22 tables natives, python-docx)",
        "zeng2018_table3_pvalues_by_drug": ZENG_TABLE3_PVALUES,
        "zeng2018_table14_adjusted_p_significant_only": ZENG_TABLE14_ADJUSTED_P,
        "zeng2018_table20_cross_cohort_validation_frequency": ZENG_TABLE20_VALIDATION,
        "tbannotator_igr_top10_variants": [
            {"spdi": v[0], "mutation": v[1], "strain_count": v[2], "rif_r": v[3], "rif_s": v[4], "rif_na": v[5]}
            for v in IGR_TOP10_VARIANTS
        ],
        "lineage_purity_top_variants": LINEAGE_PURITY,
        "rif_baseline_genome_wide": RIF_BASELINE,
        "fisher_test_pooled_igr_vs_baseline": fisher,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
