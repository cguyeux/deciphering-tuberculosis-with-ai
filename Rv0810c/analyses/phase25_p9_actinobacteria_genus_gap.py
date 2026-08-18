#!/usr/bin/env python3
"""
P9 -- DUF3073 est-il vraiment UNIVERSEL chez les Actinobacteries, ou existe-t-il
des pertes authentiques ? Voie alternative identifiee en P9.1 point (5) : au lieu
d'un calcul d'orthologie neuf (tblastn/eggNOG a l'echelle de dizaines de milliers
de genomes, cout refuse en P4.4 pour un calcul comparable), comparer directement
la liste des taxons InterPro/PF11273 a une liste de reference des genres
Actinobacteries reconnus -- sans lancer un seul BLAST.

Modele nul de la piste (pistes.md P9), applique ici a l'echelle du GENRE (pas du
genome individuel, hors de portee sans calcul lourd) :
  (a) ne retenir comme candidat un genre "manquant" que s'il possede au moins un
      assemblage RefSeq de qualite "Complete Genome" ou "Chromosome" -- un genre
      jamais bien sequence n'est pas un candidat a une perte authentique, juste
      une absence de donnee (cf. avertissement de la piste elle-meme).
  (b) et (c) (synteny, confirmation par seconde methode) restent hors de portee de
      CETTE passe : elles s'appliquent aux candidats individuels qui survivent au
      filtre (a), pas a l'ensemble -- a instruire separement si des candidats
      emergent (explicitement le cas ici, cf. resultat).

Trois sources, aucune necessitant de calcul BLAST :
  1. NCBI Taxonomy E-utils : tous les taxons de rang "genus" sous Actinomycetia
     (classe, PAS le phylum Actinobacteria/Actinomycetota, cf. controle documente ci-dessous)
     (esearch, puis esummary par lots de 200).
  2. InterPro API (deja utilisee en P0.2/P9.1) : liste complete et paginee des
     taxons UniProt portant PF11273 (`taxonomy/uniprot/entry/pfam/PF11273`,
     page_size=200) -- genre extrait du premier mot du nom d'espece.
  3. assembly_summary_refseq.txt (NCBI FTP, ~240 Mo, cache local) : genre extrait
     du premier mot d'organism_name, restreint a assembly_level in
     {Complete Genome, Chromosome} -- filtre qualite du modele nul (a).

Sortie : genres Actinobacteries reconnus, bien sequences (>=1 assemblage complet/
chromosome), et pourtant ABSENTS de la liste des taxons PF11273 -- candidats a une
perte authentique, a verifier individuellement (synteny + seconde methode) avant
toute revendication.
"""
import gzip
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree

PROJ = Path(__file__).resolve().parent.parent
DATA = PROJ / "data"
OUT = PROJ / "résultats" / "p9_actinobacteria_genus_gap.json"

ASSEMBLY_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
ASSEMBLY_SUMMARY_LOCAL = DATA / "assembly_summary_refseq.txt.gz"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
INTERPRO_TAXA_URL = "https://www.ebi.ac.uk/interpro/api/taxonomy/uniprot/entry/pfam/PF11273/?page_size=200"

UA = {"User-Agent": "Rv0810c-P9-research-script/1.0 (guyeux@gmail.com)"}


def fetch(url, retries=5, sleep=2.0):
    last_err = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except (urllib.error.URLError, TimeoutError) as e:
            last_err = e
            time.sleep(sleep * (i + 1))
    raise RuntimeError(f"echec fetch apres {retries} tentatives : {url}") from last_err


def genus_of(name):
    """Premier mot taxonomiquement informatif d'un nom d'organisme."""
    words = name.strip().split()
    if not words:
        return None
    if words[0] in ("Candidatus", "[Candidatus]", "uncultured", "unclassified"):
        return words[1] if len(words) > 1 else None
    w = words[0].strip("[]")
    if not w or not w[0].isupper() or not w.isalpha():
        return None
    return w


# ---------------------------------------------------------------------------
# 1. Genres NCBI reconnus sous Actinobacteria (E-utils, esearch + esummary)
# ---------------------------------------------------------------------------
def ncbi_actinobacteria_genera():
    # "Actinobacteria" (l'ancien phylum, aujourd'hui Actinomycetota) englobe aussi les
    # classes tres divergentes Coriobacteriia, Rubrobacteria, Acidimicrobiia,
    # Thermoleophilia, Nitriliruptoria -- controle direct (2026-08-11) : AUCUN genre de
    # ces classes (Coriobacterium, Atopobium, Collinsella, Rubrobacter, Conexibacter,
    # Acidimicrobium...) n'est jamais couvert par PF11273, alors que le domaine y est
    # cherche depuis 2006 (Gao). Restreindre a la classe Actinomycetia (Mycobacteriales,
    # Streptomycetales, Corynebacteriales...), la seule ou le domaine a jamais ete trouve
    # -- sinon la totalite de l'ecart genre-par-genre n'est qu'un artefact de perimetre
    # taxonomique trop large (phylum vs classe), pas un signal de perte authentique.
    term = "Actinomycetia[Subtree] AND genus[Rank]"
    esearch_url = f"{EUTILS}/esearch.fcgi?db=taxonomy&term={urllib.parse.quote(term)}&retmax=2000"
    xml = fetch(esearch_url)
    root = ElementTree.fromstring(xml)
    ids = [e.text for e in root.findall(".//Id") if e.text]

    genera = {}  # taxid -> nom de genre
    batch_size = 200
    for i in range(0, len(ids), batch_size):
        batch = ids[i : i + batch_size]
        esum_url = f"{EUTILS}/esummary.fcgi?db=taxonomy&id={','.join(batch)}"
        xml = fetch(esum_url)
        root = ElementTree.fromstring(xml)
        for doc in root.findall(".//DocSum"):
            id_el = doc.find("Id")
            if id_el is None or not id_el.text:
                continue
            taxid = id_el.text
            name_el = [it for it in doc.findall("Item") if it.get("Name") == "ScientificName"]
            if name_el:
                genera[taxid] = name_el[0].text
        time.sleep(0.34)  # E-utils sans cle API : max 3 req/s
    return genera


# ---------------------------------------------------------------------------
# 2. Genres portant PF11273 selon InterPro (pagination complete)
# ---------------------------------------------------------------------------
def interpro_pf11273_genera():
    url = INTERPRO_TAXA_URL
    genera = set()
    species_names = []
    while url:
        payload = json.loads(fetch(url))
        for res in payload["results"]:
            name = res["metadata"]["name"]
            species_names.append(name)
            g = genus_of(name)
            if g:
                genera.add(g)
        url = payload.get("next")
        time.sleep(0.2)
    return genera, species_names


# ---------------------------------------------------------------------------
# 3. Genres Actinobacteries avec >=1 assemblage RefSeq Complete Genome/Chromosome
# ---------------------------------------------------------------------------
def download_assembly_summary():
    if ASSEMBLY_SUMMARY_LOCAL.exists():
        return
    DATA.mkdir(exist_ok=True)
    print("Telechargement assembly_summary_refseq.txt (~240 Mo, une seule fois)...")
    raw = fetch(ASSEMBLY_SUMMARY_URL, retries=3)
    with gzip.open(ASSEMBLY_SUMMARY_LOCAL, "wb") as f:
        f.write(raw)
    print("OK :", ASSEMBLY_SUMMARY_LOCAL, f"({ASSEMBLY_SUMMARY_LOCAL.stat().st_size / 1e6:.1f} Mo compresse)")


def well_sequenced_genera(candidate_genus_names):
    """Genres (parmi candidate_genus_names) avec >=1 assemblage RefSeq de
    niveau Complete Genome ou Chromosome, comptage inclus."""
    download_assembly_summary()
    counts = defaultdict(int)
    with gzip.open(ASSEMBLY_SUMMARY_LOCAL, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            organism_name = cols[7]
            assembly_level = cols[11]
            if assembly_level not in ("Complete Genome", "Chromosome"):
                continue
            g = genus_of(organism_name)
            if g in candidate_genus_names:
                counts[g] += 1
    return counts


def main():
    print("1/3 -- Genres NCBI sous Actinomycetia (E-utils)...")
    ncbi_genera = ncbi_actinobacteria_genera()  # taxid -> nom
    ncbi_genus_names = set(ncbi_genera.values())
    print(f"   {len(ncbi_genus_names)} genres NCBI trouves")

    print("2/3 -- Genres couverts par PF11273 (InterPro, pagination)...")
    ip_genera, ip_species_names = interpro_pf11273_genera()
    print(f"   {len(ip_species_names)} taxons InterPro lus -> {len(ip_genera)} genres distincts")

    print("3/3 -- Filtrage qualite : genres NCBI avec assemblage Complete/Chromosome...")
    quality_counts = well_sequenced_genera(ncbi_genus_names)
    well_sequenced = set(quality_counts)
    print(f"   {len(well_sequenced)} / {len(ncbi_genus_names)} genres NCBI ont >=1 assemblage de qualite")

    missing = sorted(well_sequenced - ip_genera)
    covered = sorted(well_sequenced & ip_genera)
    not_sequenced = sorted(ncbi_genus_names - well_sequenced)

    result = {
        "piste": "P9",
        "methode": (
            "Comparaison genre-par-genre, sans BLAST : (1) genres NCBI sous Actinomycetia "
            "(E-utils esearch+esummary), (2) genres couverts par InterPro/PF11273 (pagination "
            "complete de taxonomy/uniprot/entry/pfam/PF11273), (3) filtre qualite -- ne retenir "
            "comme candidat que les genres avec >=1 assemblage RefSeq Complete Genome/Chromosome "
            "(assembly_summary_refseq.txt), pour respecter le modele nul de la piste (une absence "
            "chez un genre jamais bien sequence n'est pas un candidat)."
        ),
        "n_genres_ncbi_actinobacteria": len(ncbi_genus_names),
        "n_genres_interpro_pf11273": len(ip_genera),
        "n_taxons_interpro_lus": len(ip_species_names),
        "n_genres_bien_sequences_complete_ou_chromosome": len(well_sequenced),
        "n_genres_bien_sequences_ET_couverts_pf11273": len(covered),
        "n_genres_bien_sequences_NON_couverts_pf11273": len(missing),
        "n_genres_jamais_bien_sequences_exclus_du_test": len(not_sequenced),
        "genres_candidats_perte_authentique": [
            {"genre": g, "n_assemblages_complete_ou_chromosome": quality_counts[g]} for g in missing
        ],
        "genres_couverts_pour_controle": covered,
    }

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    print(f"\nEcrit : {OUT}")
    print(f"Genres NCBI Actinomycetia : {len(ncbi_genus_names)}")
    print(f"Genres bien sequences (Complete/Chromosome) : {len(well_sequenced)}")
    print(f"  -> couverts par PF11273 : {len(covered)}")
    print(f"  -> CANDIDATS perte authentique (bien sequences, PF11273 absent) : {len(missing)}")
    if missing:
        print("  Liste :", ", ".join(missing))


if __name__ == "__main__":
    main()
