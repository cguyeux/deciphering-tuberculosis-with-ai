"""P8 -- Formaliser la contrainte purificatrice en dN/dS (test type McDonald-Kreitman),
via le moteur statistique du skill mk-ascertainment.

Modele nul de la piste (pistes.md P8) : comparer le ratio NS:S de Rv0810c a celui de la
MEME classe de temoins appariee en longueur deja utilisee en P2.2 (proteines H37Rv 50-70
aa), plutot qu'au proteome entier.

Donnees : requete TBannotator (tb_report_spdi_annotations, TOUS les sites polymorphes
distincts, sans seuil de frequence) pour Rv0810c + les 74 temoins de P2.2 (meme fichier
proteome.faa, meme fenetre 50-70 aa). Sauvegarde brute : data/p8_ns_s_sites_75_controls.csv

IMPORTANT -- ecart decouvert avec P4.3 : la requete exhaustive ci-dessous recense 84 sites
missense pour Rv0810c (dont p.Thr24Pro a 32 souches, p.Asp56Ala a 73 souches), tres au-dela
des "3 sites missense a l'etat de trace" que P4.3 avait consignes. P4.3 n'avait PAS interroge
tb_report_spdi_annotations de facon exhaustive (sa liste "SNP"/"complex" etait un
sous-ensemble construit a la main pour une autre question, la disruption clonale). Ce
script utilise donc la requete exhaustive comme source de verite pour le test MK, et
signale l'ecart pour correction ulterieure de P4.3 / etat_des_decouvertes.md.

Moteur statistique reutilise du skill mk-ascertainment (mk_per_gene.mk_stats), pour
beneficier des memes formules DoS/NI/alpha/Fisher que le reste de l'ecosysteme MTBC,
adaptees ici a une comparaison POLYMORPHISME vs POLYMORPHISME (Rv0810c vs classe temoin),
pas DIVERGENCE vs POLYMORPHISME (design normalement vise par ce skill pour une
sous-lignee) -- adaptation explicitement justifiee et documentee dans pistes.md.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(
    0,
    "/home/christophe/docs/codes/claude_plugins/bio_pathogens/skills/mk-ascertainment/scripts",
)
from mk_per_gene import mk_stats  # noqa: E402

PROJ = Path(__file__).resolve().parent.parent
DATA = PROJ / "data" / "p8_ns_s_sites_75_controls.csv"
OUT = PROJ / "résultats" / "p8_mk_test.json"

NS_TYPES = {"missense_variant", "stop_gained", "stop_lost"}
S_TYPES = {"synonymous_variant"}


def load_counts():
    counts = defaultdict(lambda: {"NS": 0, "S": 0})
    with open(DATA) as f:
        for row in csv.DictReader(f):
            gene, ann, n = row["locus_tag"], row["annotation_type"], int(row["n_sites"])
            if ann in NS_TYPES:
                counts[gene]["NS"] += n
            elif ann in S_TYPES:
                counts[gene]["S"] += n
    return counts


def main():
    import json

    counts = load_counts()
    target = counts.pop("Rv0810c")
    controls = counts  # 74 genes

    ns_t, s_t = target["NS"], target["S"]
    ratio_t = ns_t / s_t if s_t else float("inf")
    prop_t = ns_t / (ns_t + s_t)

    # ---- 1. Rang empirique parmi les 75 (les 74 temoins + Rv0810c), a la P2.2 ----------
    all_props = [prop_t] + [c["NS"] / (c["NS"] + c["S"]) for c in controls.values()]
    # rang par comptage explicite (gere les ex-aequo proprement)
    n_leq = sum(1 for p in all_props if p <= prop_t)
    percentile = 100 * n_leq / len(all_props)

    # ---- 2. Fisher exact : Rv0810c (NS,S) vs POOL des 74 temoins (NS,S) ------------------
    ns_pool = sum(c["NS"] for c in controls.values())
    s_pool = sum(c["S"] for c in controls.values())
    dos, ni, alpha, p_fisher = mk_stats(ns_t, s_t, ns_pool, s_pool)
    # ici Dn/Ds = Rv0810c (classe "test"), Pn/Ps = pool des temoins (classe "reference") ;
    # relecture : DoS>0 et alpha>0 signifient Rv0810c PLUS charge en NS que le pool
    # (signature d'evolution moins contrainte, PAS de purifying selection) ; DoS<0 et
    # alpha<0 signifient l'inverse (Rv0810c MOINS charge en NS => plus contraint).

    # ---- 3. Meme test, technique 2 : chaque temoin prend le role de comparateur seul ----
    # (pour ne pas laisser le pool ecraser par les 4 genes extremes Rv1576c/Rv2082/
    # Rv3402c/Rv3222c) -- distribution des p-values et DoS individuels
    per_control_tests = []
    for gene, c in controls.items():
        dos_i, _ni_i, _alpha_i, p_i = mk_stats(ns_t, s_t, c["NS"], c["S"])
        per_control_tests.append({
            "gene": gene, "NS": c["NS"], "S": c["S"],
            "NS_over_S": round(c["NS"] / c["S"], 4) if c["S"] else None,
            "DoS_vs_Rv0810c": round(dos_i, 4), "p_fisher": p_i,
        })
    n_sig_more_constrained = sum(
        1 for t in per_control_tests
        if t["p_fisher"] is not None and t["p_fisher"] < 0.05 and t["DoS_vs_Rv0810c"] < 0
    )
    n_sig_less_constrained = sum(
        1 for t in per_control_tests
        if t["p_fisher"] is not None and t["p_fisher"] < 0.05 and t["DoS_vs_Rv0810c"] > 0
    )

    result = {
        "piste": "P8",
        "source_donnees": str(DATA),
        "avertissement_ecart_P4_3": (
            "P4.3 recensait 3 sites missense (R31Q, G35S, Q27T) et 11 sites synonymes pour "
            "Rv0810c, un sous-ensemble construit a la main pour une autre question (la "
            "disruption clonale), PAS une requete exhaustive de tb_report_spdi_annotations. "
            "La requete exhaustive de ce script trouve 84 sites missense (dont p.Thr24Pro a "
            "32 souches et p.Asp56Ala a 73 souches, jamais mentionnes dans P4.3) et 100 sites "
            "synonymes, TOUS seuils de frequence confondus, references verifiees sur le "
            "genome. A signaler pour correction de etat_des_decouvertes.md §2 / P4.3."
        ),
        "Rv0810c": {
            "NS_missense_stop": ns_t, "S_synonymous": s_t,
            "NS_over_S": round(ratio_t, 4), "NS_proportion": round(prop_t, 4),
        },
        "rang_empirique_parmi_75": {
            "n_temoins": len(controls), "rang_1_est_le_moins_charge_NS": n_leq,
            "percentile": round(percentile, 1),
            "mediane_temoins_NS_proportion": round(
                sorted(c["NS"] / (c["NS"] + c["S"]) for c in controls.values())[len(controls) // 2], 4
            ),
        },
        "fisher_vs_pool_74_temoins": {
            "NS_pool": ns_pool, "S_pool": s_pool,
            "NS_over_S_pool": round(ns_pool / s_pool, 4),
            "DoS": round(dos, 4), "NI": round(ni, 4), "alpha_1_moins_NI": round(alpha, 4),
            "p_fisher": p_fisher,
            "lecture": (
                "DoS<0 et alpha<0 : Rv0810c est MOINS charge en NS/S que le pool des temoins "
                "(plus contraint) ; DoS>0 : l'inverse."
                if dos < 0 else
                "DoS>0 et alpha>0 : Rv0810c est PLUS charge en NS/S que le pool des temoins "
                "(moins contraint que la classe de reference), CONTRAIRE a l'hypothese de "
                "contrainte purificatrice accrue."
            ),
        },
        "tests_individuels_vs_chaque_temoin": {
            "n_temoins_testes": len(per_control_tests),
            "n_significatifs_Rv0810c_plus_contraint_p<0.05": n_sig_more_constrained,
            "n_significatifs_Rv0810c_moins_contraint_p<0.05": n_sig_less_constrained,
            "detail": sorted(per_control_tests, key=lambda t: (t["p_fisher"] if t["p_fisher"] is not None else 1)),
        },
    }

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    print(f"Ecrit : {OUT}\n")
    print("Rv0810c : NS =", ns_t, " S =", s_t, " NS/S =", round(ratio_t, 4))
    print("Rang parmi 75 (percentile, NS proportion) :", round(percentile, 1))
    print("Pool 74 temoins : NS =", ns_pool, " S =", s_pool, " NS/S =", round(ns_pool / s_pool, 4))
    print("DoS =", round(dos, 4), " NI =", round(ni, 4), " alpha =", round(alpha, 4), " p_fisher =", p_fisher)
    print(result["fisher_vs_pool_74_temoins"]["lecture"])
    print("Tests individuels significatifs (Rv0810c plus contraint) :", n_sig_more_constrained, "/", len(per_control_tests))
    print("Tests individuels significatifs (Rv0810c moins contraint) :", n_sig_less_constrained, "/", len(per_control_tests))


if __name__ == "__main__":
    main()
