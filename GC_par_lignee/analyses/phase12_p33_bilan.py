#!/usr/bin/env python3
"""
Objet       : P3.3 -- construire, pour chaque lignee, un echantillon de 300
              singletons terminaux verifiables (situes dans un CDS) et
              generer la requete SQL qui, pour chacun, compte combien des 39
              (ou moins) AUTRES souches du meme pool ont ce gene dans leur
              liste TBannotator de genes manquants (tb_report_strain_missing_
              gene). Un singleton dont AUCUN autre membre du pool n'a le gene
              manquant est CLEAN (l'absence de variant chez les autres se lit
              bien comme « porte H37Rv ») ; un singleton dont au moins un
              autre membre a le gene manquant est SUSPECT (leur silence peut
              signifier « non couvert », pas « refere »).
              Design en echantillon plutot qu'exhaustif : P2.6 a deja etabli
              que la puissance de comptage n'est pas le facteur limitant du
              projet (IC95 bien plus etroit que l'ecart inter-lignees) ; un
              echantillon de 300 par lignee donne un IC95 de quelques points
              sur le taux de suspicion, largement suffisant pour statuer.
              Requete poussee cote SQL (CTE + agregation), pas transferee
              ligne a ligne : le jeu de missing genes complet d'un pool
              (des dizaines de milliers de lignes) n'a pas besoin de quitter
              PostgreSQL, seul le resultat agrege (une ligne par lignee)
              revient au client.
Entrees     : data/p33_sample_query_input.json (produit par
              phase12_p33_couverture_singletons.py + un tirage aleatoire
              seed=0, 300 singletons CDS par lignee) : pour chaque lignee, la
              liste des strain_id du pool et les (locus_tag, carrier_strain_id)
              du sous-echantillon.
              TBannotator PostgreSQL (tb_report_gene, tb_report_strain_missing_
              gene), interroge ICI VIA LE MCP tbannotator DANS LA SESSION
              CLAUDE -- ce script documente et peut regenerer le texte SQL
              exact envoye pour chaque lignee (fonction build_sql), mais ne
              l'execute pas lui-meme (pas de client Postgres direct dans ce
              depot ; toute requete transite par le MCP, jamais par un mot de
              passe stocke localement).
Sorties     : resultats/phase12_p33_bilan_clade.tsv (n_sample, n_suspect,
              % suspect, pire cas par lignee) et
              resultats/phase12_p33_suspect_loci.tsv (detail des singletons
              suspects, un par ligne, avec le nombre d'autres souches du pool
              ou le gene est marque manquant) -- ces deux TSV sont ecrits a la
              main a partir des resultats de la session du 2026-09-01, ils ne
              sont pas regeneres par ce script.
Reutilisable: oui -- le patron (VALUES inline + CTE d'agregation, un aller-
              retour par lignee au lieu d'un transfert ligne a ligne) vaut
              pour toute verification de couverture genique via TBannotator
              sur un sous-ensemble cible de positions, dans n'importe quel
              projet MTBC du depot.
Projet      : GC_par_lignee
Date        : 2026-09-01
"""
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from phase12_p33_couverture_singletons import build_panel, terminal_singletons  # noqa: E402
from phase2_polarisation_mtbc0 import build_ancestral, load_mask  # noqa: E402
from phase3_counts_par_souche import MASKS  # noqa: E402
from phase12_p33_couverture_singletons import locus_tag_by_position, CLADES  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SAMPLE_SIZE = 300
SEED = 0


def build_sample_input():
    """Reproduit data/p33_sample_query_input.json depuis zero (traçabilite).
    Necessite le meme panel n=40/seed=0 que P3.2/P3.3 et la meme table
    strain_name -> strain_id que data/p33_strain_ids.csv (produite par une
    requete MCP separee, non regeneree ici)."""
    import csv
    names_to_id = {}
    id_csv = ROOT / "data" / "p33_strain_ids.csv"
    for row in csv.DictReader(id_csv.open()):
        names_to_id[row["strain_name"]] = row["strain_id"]

    anc = build_ancestral()
    masked = load_mask(MASKS, in_mtbc0_coords=True)
    owner = locus_tag_by_position()
    rng = random.Random(SEED)

    out = {}
    for clade in CLADES:
        names, subsets = build_panel(clade)
        if len(subsets) < 4:
            continue
        term = terminal_singletons(subsets, masked, anc)
        pool_ids = [names_to_id[n] for n in names if n in names_to_id]
        checkable = [(owner[pos], names[i]) for (pos, _ref, _alt), (i, _status) in term.items()
                     if pos in owner]
        sample = checkable if len(checkable) <= SAMPLE_SIZE else rng.sample(checkable, SAMPLE_SIZE)
        sample_ids = [(lt, names_to_id[carrier]) for lt, carrier in sample
                      if carrier in names_to_id]
        out[clade] = {"pool_ids": pool_ids, "sample": sample_ids}
    return out


def build_sql(pool_ids, sample_ids):
    """Le texte SQL exact envoye au MCP tbannotator pour une lignee : les
    singletons cibles et le pool en VALUES inline, jointure et agregation
    cote serveur, un seul resultat (n_sample, n_suspect, detail) qui revient."""
    targets = ",".join(f"('{lt}',{cid})" for lt, cid in sample_ids)
    pool = ",".join(f"({sid})" for sid in pool_ids)
    return (
        f"WITH targets(locus_tag, carrier_id) AS (VALUES {targets}), "
        f"pool(strain_id) AS (VALUES {pool}), "
        "per_target AS (SELECT t.locus_tag, t.carrier_id, "
        "count(DISTINCT smg.strain_id) AS n_others_missing "
        "FROM targets t JOIN pool p ON p.strain_id <> t.carrier_id "
        "LEFT JOIN tb_report_gene g ON g.locus_tag = t.locus_tag "
        "LEFT JOIN tb_report_strain_missing_gene smg "
        "ON smg.gene_id = g.gene_id AND smg.strain_id = p.strain_id "
        "GROUP BY t.locus_tag, t.carrier_id) "
        "SELECT count(*) AS n_sample, "
        "count(*) FILTER (WHERE n_others_missing > 0) AS n_suspect, "
        "string_agg(DISTINCT (locus_tag || ':' || n_others_missing), ',') "
        "FILTER (WHERE n_others_missing > 0) AS suspect_detail FROM per_target"
    )


if __name__ == "__main__":
    data = build_sample_input()
    (ROOT / "data" / "p33_sample_query_input.json").write_text(json.dumps(data))
    for clade, p in data.items():
        sql = build_sql(p["pool_ids"], p["sample"])
        print(f"# {clade} : {len(sql)} caracteres, {len(p['sample'])} singletons, "
              f"{len(p['pool_ids'])} souches de pool", file=sys.stderr)
    print("# Requetes reconstruites. Resultats deja consignes dans "
          "resultats/phase12_p33_bilan_clade.tsv (session 2026-09-01) : "
          "exec via le MCP tbannotator, pas ce script.", file=sys.stderr)
