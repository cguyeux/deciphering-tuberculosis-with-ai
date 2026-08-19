"""Chemins centralisés du projet bpal_resistance_emergence.

Ce projet est CONSOMMATEUR de l'infrastructure du projet voisin
`Resistance_antibio/` (catalogue WHO consolidé, pangénome SPDI, phénotypes
consensus, snapshot de lignée). Il ne réimplémente rien de ce socle ; son
angle propre est l'émergence évolutive / épidémiologique de la résistance au
régime BPaLM.

    from analyses import paths   # ou : import paths (si lancé depuis analyses/)
"""
from pathlib import Path

# --- Racines ---------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent          # bpal_resistance_emergence/
MTBC = ROOT.parent                                     # mtbc/
RA = MTBC / "Resistance_antibio"                       # projet voisin (à consommer)

# --- Local -----------------------------------------------------------------
ANALYSES = ROOT / "analyses"
DATA = ROOT / "data"
RESULTATS = ROOT / "résultats"
ARTICLE = ROOT / "article"
LITREVIEW = ROOT / "litterature_review"
GENE_PANEL = DATA / "gene_panel.tsv"

# --- Actifs réutilisés de Resistance_antibio (read-only, ne pas modifier) --
CATALOGUE_TSV = RA / "catalogue" / "catalogue_consolide.tsv"
PHENOTYPES_TSV = RA / "catalogue" / "phenotypes_souches_consensus.tsv"
WHO_COORDS = RA / "data/sources/who_catalogue/WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"
PAN_SPDI = RA / "data/pangenome/pan_spdi.pkl"          # {SPDI: [souches]}
PAN_STRAINS = RA / "data/pangenome/pan_strains.pkl"    # {souche: {SPDI}}
LINEAGE_SNAPSHOT = RA / "data/sources/tbannotator/lineage_coll_snapshot.csv"  # strain_name,lineage_code,lineage_level_1
RA_MODEL = RA / "résultats/model"                      # X.npz, features.txt, rows.tsv, lineage_mask.tsv

# --- Références génomiques (pipeline central) ------------------------------
RESOURCES = MTBC / "investigate_phylo/resources"
H37RV_FASTA = RESOURCES / "NC_000962.3.fasta"
H37RV_GFF3 = RESOURCES / "NC_000962.3.gff3"
H37RV_CDS = RESOURCES / "NC_000962.3_CDS.fasta"

# --- Données partagées MTBC ------------------------------------------------
BDD = MTBC / "bdd/actuelle"                            # {clade}/{SRA}/NC_000962.3/{spdi.txt,report.json}
GLOBAL_SUPP = MTBC / "global_supplementary"
BARCODE_COMPLETE = GLOBAL_SUPP / "barcoding_v2/barcode_complete.tsv"
BIOPROJECT_GEO = GLOBAL_SUPP / "bioproject_geo"        # consolidated_geo_*.tsv (pays + dates)

# --- Scripts packagés des skills (canonical, ne pas éditer via symlink) ----
SKILLS = MTBC.parent / "claude_plugins/bio_pathogens/skills"
SK_RES_CATALOGUE = SKILLS / "resistance-catalogue/scripts"     # query_catalogue.py, hgvs_to_spdi.py
SK_RES_DISCOVERY = SKILLS / "resistance-discovery/scripts"     # variant_match.py (réconciliation SNP+MNV)
SK_RES_PREDICT = SKILLS / "resistance-predict/scripts"         # resistance_profile.py
SK_CONVERGENT = SKILLS / "convergent-evolution/scripts"        # convergent_evolution.py

REF_ACC = "NC_000962.3"


def ensure_dirs():
    for d in (DATA, RESULTATS, ARTICLE, LITREVIEW):
        d.mkdir(parents=True, exist_ok=True)
