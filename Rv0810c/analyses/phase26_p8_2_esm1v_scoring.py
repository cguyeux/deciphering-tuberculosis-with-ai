#!/usr/bin/env python3
"""
P8.2 -- Noter par ESM-1v (Meier 2021, LLR masque, seule mesure citable, cf.
Evidence status du skill mtbc-gene) les 84 sites missense du catalogue P8
(data/p8_2_missense_catalogue_84.csv, extrait le 2026-08-11 via TBannotator MCP :
tb_report_spdi_annotations x tb_report_spdi x mv_snp_frequency pour les SNP
simples, + tb_report_strain_spdi pour les 18 variants multi-nucleotidiques non
couverts par la vue materialisee mv_snp_frequency ; 3 hgvs_p multi-residus exclus
du scoring, 81/84 sites scorables).

Deux tests prescrits par la piste (pistes.md P8.2), chacun avec son modele nul :

(a) La frequence observee decroit-elle avec le score de dommage predit (sous
    selection purificatrice dose-dependante, un variant plus frequent devrait
    etre note MOINS deletere) ?
    Modele nul : refaire le test hors queue acide (34-60), ou P2.4 a deja montre
    une contrainte de COMPOSITION sans conservation positionnelle stricte qui
    pourrait produire une correlation triviale indépendante de toute selection
    sur le site precis.

(b) Les quatre phosphosites (P1.5 : T24, S21, S20 dans le module 1-33 ; S51 dans
    la queue) sont-ils enrichis en substitutions delateres par rapport au reste
    du module structure, plutot qu'a la proteine entiere (P3.3 a deja montre le
    module entier tres conserve : comparer a la proteine entiere confondrait
    "le module est conserve" avec "les phosphosites le sont plus que le reste
    du module") ?
    Modele nul, implemente EXACTEMENT ici : saturation complete (19 substitutions
    possibles) des 33 positions du module, puis test de permutation EXACT
    (les C(33,3)=5456 triplets de positions, pas un bootstrap) de la moyenne des
    scores de saturation aux positions {20,21,24} contre tous les triplets
    possibles de positions du module. S51 (queue) est hors du perimetre de ce
    test par construction -- rapporte a part, jamais mele au module.

Contre-argument de la piste, garde en tete dans la lecture des resultats : ESM-1v
est un score de conservation de sequence, donc correle par construction avec ce
que P4.3/P8 mesurent deja (NS/S) -- un resultat positif au test (a) confirme la
coherence de la methode plus qu'il n'apporte un fait totalement independant.
"""
import csv
import itertools
import json
import os
import sys
from pathlib import Path

# Machine partagee, souvent en forte contention (2026-08-11 : load average 35-42 sur
# 16 coeurs, swap quasi plein, plusieurs sessions MTBC soeurs actives) -- brider les
# threads BLAS/torch reduit l'empreinte de CE job au lieu d'aggraver le thrashing en
# tentant de saturer des coeurs deja sursouscrits.
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

sys.path.insert(0, "/home/christophe/docs/codes/claude_plugins/bio_pathogens/skills/mtbc-gene/src")
sys.path.insert(
    0, "/home/christophe/docs/codes/claude_plugins/bio_population_genetics/skills/esm-atlas-cli/src"
)
from mtbc_mutation_impact.llr import load_model, llr_of  # noqa: E402

import torch  # noqa: E402

torch.set_num_threads(2)

PROJ = Path(__file__).resolve().parent.parent
CATALOGUE = PROJ / "data" / "p8_2_missense_catalogue_84.csv"
OUT = PROJ / "résultats" / "p8_2_esm1v_scoring.json"

GENE = "Rv0810c"
SEQ = "MGRGRAKAKQTKVARELKYSSPQTDFQRLQRELSGTGTDRLDGDGPSDDDSWNDEDDWRR"
MODULE_POSITIONS = list(range(1, 34))  # 1-33, module rigide (P2.1)
PHOSPHOSITES_MODULE = [20, 21, 24]  # S20, S21, T24 -- P1.5, dans le module
PHOSPHOSITE_TAIL = 51  # S51 -- P1.5, dans la queue desordonnee, hors module
AA20 = "ACDEFGHIKLMNPQRSTVWY"


def load_catalogue():
    rows = list(csv.DictReader(open(CATALOGUE)))
    for r in rows:
        r["protein_position"] = int(r["protein_position"])
        r["frequency"] = int(r["frequency"])
        r["scoreable"] = r["scoreable"] == "1"
    return rows


def build_module_saturation_pairs():
    pairs = []
    for pos in MODULE_POSITIONS:
        wt = SEQ[pos - 1]
        for alt in AA20:
            if alt == wt:
                continue
            pairs.append(f"{wt}{pos}{alt}")
    return pairs


CHECKPOINT = PROJ / "résultats" / "p8_2_esm1v_checkpoint.jsonl"


def load_checkpoint():
    """Reprend un run precedent tue avant la fin (ex. 2026-08-11 : job perdu apres
    ~1h de calcul CPU faute de sauvegarde incrementale -- corrige ici une fois pour
    toutes plutot que de relancer aveuglement)."""
    cache = {}
    if CHECKPOINT.exists():
        with open(CHECKPOINT) as f:
            for line in f:
                rec = json.loads(line)
                cache[rec["mutation"]] = rec["llr"]
        print(f"Checkpoint trouve : {len(cache)} mutations deja notees, reprise.")
    return cache


def score_all(mutation_codes, model, alphabet):
    """Score chaque mutation une seule fois (dedoublonnage), reutilise le modele
    deja charge -- score_many() du module recharge le modele a chaque appel, on
    reimplemente donc la boucle ici avec le modele partage. Sauvegarde CHAQUE score
    immediatement sur disque (append + flush) : un kill externe (timeout, OOM) ne
    doit plus jamais faire perdre le calcul deja fait."""
    cache = load_checkpoint()
    todo = [m for m in sorted(set(mutation_codes)) if m not in cache]
    print(f"{len(cache)} deja en cache, {len(todo)} restant a noter.")
    CHECKPOINT.parent.mkdir(exist_ok=True)
    with open(CHECKPOINT, "a") as ckpt:
        for i, mut in enumerate(todo):
            try:
                llr = llr_of(GENE, mut, model, alphabet)["llr"]
            except Exception as exc:
                llr = None
                print(f"  echec {mut} : {exc}")
            cache[mut] = llr
            ckpt.write(json.dumps({"mutation": mut, "llr": llr}) + "\n")
            ckpt.flush()
            if (i + 1) % 20 == 0:
                print(f"  ... {i + 1}/{len(todo)} mutations notees (cette session)")
    return cache


def spearman(xs, ys):
    from scipy.stats import spearmanr

    rho, p = spearmanr(xs, ys)
    return {"rho": round(float(rho), 4), "p": float(p), "n": len(xs)}


def main():
    rows = load_catalogue()
    observed_codes = [r["mutation_code"] for r in rows if r["scoreable"]]
    saturation_codes = build_module_saturation_pairs()
    all_codes = set(observed_codes) | set(saturation_codes)
    print(f"Sites catalogue scorables : {len(observed_codes)} ; saturation module (33x19) : "
          f"{len(saturation_codes)} ; total unique a noter : {len(all_codes)}")

    print("Chargement du modele ESM-1v (une seule fois)...")
    model, alphabet = load_model()
    print("Modele charge. Notation en cours...")
    scores = score_all(all_codes, model, alphabet)

    # ---- catalogue annote -------------------------------------------------
    for r in rows:
        r["llr"] = scores.get(r["mutation_code"]) if r["scoreable"] else None

    scored_rows = [r for r in rows if r["scoreable"] and r["llr"] is not None]

    # ---- test (a) : frequence vs LLR, protéine entiere puis hors queue ----
    def freq_llr_test(subset, label):
        xs = [r["frequency"] for r in subset]
        ys = [r["llr"] for r in subset]
        res = spearman(xs, ys)
        res["label"] = label
        res["lecture"] = (
            "rho>0 attendu sous selection purificatrice dose-dependante "
            "(frequence elevee <-> LLR moins negatif, moins delatere)"
        )
        return res

    test_a_full = freq_llr_test(scored_rows, "proteine entiere (81 sites)")
    test_a_module = freq_llr_test(
        [r for r in scored_rows if r["protein_position"] <= 33], "module seul (1-33)"
    )
    test_a_tail = freq_llr_test(
        [r for r in scored_rows if r["protein_position"] > 33], "queue seule (34-60, modele nul)"
    )

    # ---- test (b) : phosphosites module vs permutation exacte -------------
    sat = {pos: [scores[f"{SEQ[pos - 1]}{pos}{alt}"] for alt in AA20 if alt != SEQ[pos - 1]]
           for pos in MODULE_POSITIONS}
    sat_mean = {pos: sum(v) / len(v) for pos, v in sat.items() if all(x is not None for x in v)}

    observed_stat = sum(sat_mean[p] for p in PHOSPHOSITES_MODULE) / len(PHOSPHOSITES_MODULE)
    all_triplets = list(itertools.combinations(MODULE_POSITIONS, 3))
    null_means = [sum(sat_mean[p] for p in trip) / 3 for trip in all_triplets]
    n_leq = sum(1 for m in null_means if m <= observed_stat)
    p_exact_one_sided = n_leq / len(null_means)  # H1 : phosphosites PLUS delateres (LLR plus bas)

    phosphosite_observed_llr = {
        r["mutation_code"]: r["llr"]
        for r in scored_rows
        if r["protein_position"] in PHOSPHOSITES_MODULE
    }
    tail_phosphosite_llr = {
        r["mutation_code"]: r["llr"]
        for r in scored_rows
        if r["protein_position"] == PHOSPHOSITE_TAIL
    }

    test_b = {
        "positions_module_testees": PHOSPHOSITES_MODULE,
        "position_queue_hors_test": PHOSPHOSITE_TAIL,
        "methode": (
            "Moyenne des scores de SATURATION (19 substitutions possibles/position) aux "
            "positions {20,21,24}, comparee par permutation EXACTE (les 5456 triplets de "
            "positions du module 1-33, pas un bootstrap) a la meme statistique pour tout "
            "triplet de positions du module. H1 unilaterale : les phosphosites sont PLUS "
            "delateres (LLR moyen plus negatif) que des positions aleatoires du module."
        ),
        "llr_moyen_saturation_phosphosites_module": round(observed_stat, 4),
        "llr_moyen_saturation_module_ensemble": round(
            sum(sat_mean.values()) / len(sat_mean), 4
        ),
        "n_triplets_possibles": len(all_triplets),
        "n_triplets_aussi_ou_plus_delateres": n_leq,
        "p_exact_unilateral": round(p_exact_one_sided, 4),
        "llr_saturation_par_position": {p: round(m, 4) for p, m in sat_mean.items()},
        "llr_substitutions_observees_aux_phosphosites_module": phosphosite_observed_llr,
        "llr_substitutions_observees_S51_queue_hors_test": tail_phosphosite_llr,
        "avertissement_S51": (
            "S51 est dans la queue desordonnee (34-60), pas dans le module (1-33) : "
            "exclu du test de permutation par construction (le modele nul de la piste ne "
            "s'applique qu'au module). Valeur(s) rapportee(s) a titre descriptif seulement."
        ),
    }

    result = {
        "piste": "P8.2",
        "gene": GENE, "sequence": SEQ,
        "n_sites_catalogue": len(rows), "n_scorables": len(observed_codes),
        "n_non_scorables_multi_residus": len(rows) - len(observed_codes),
        "modele_esm1v": "esm1v_t33_650M_UR90S_1 (Meier et al. 2021)",
        "catalogue_note": [
            {k: r[k] for k in ("hgvs_p", "protein_position", "frequency", "mutation_code", "llr", "source")}
            for r in rows
        ],
        "test_a_frequence_vs_llr": {
            "proteine_entiere": test_a_full,
            "module_1_33": test_a_module,
            "queue_34_60": test_a_tail,
        },
        "test_b_phosphosites_vs_module": test_b,
    }

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    print(f"\nEcrit : {OUT}")
    print("Test (a) proteine entiere :", test_a_full)
    print("Test (a) module seul :", test_a_module)
    print("Test (a) queue seule :", test_a_tail)
    print("Test (b) p exact unilateral :", p_exact_one_sided)


if __name__ == "__main__":
    main()
