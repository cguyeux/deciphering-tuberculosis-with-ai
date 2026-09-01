#!/usr/bin/env python3
"""
Objet       : P3.3 -- distinguer, pour chaque substitution TERMINALE (portee
              par une seule souche du pool de sa lignee), « la souche porte
              vraiment l'allele derive et les 39 autres portent vraiment H37Rv »
              de « les 39 autres ne sont simplement pas couvertes a cette
              position », piege deja paye ailleurs dans le depot (KB
              tuberculosis.md). Le disque local ne porte que spdi.txt (la liste
              des variants appeles), jamais un depth par position ; et le
              PostgreSQL de TBannotator (verifie ici par tool_get_schema) ne
              stocke lui non plus la profondeur qu'aux positions ou un variant
              a ete appele (table tb_report_strain_spdi), jamais un depth
              genome entier. Le seul proxy de couverture reellement peuple est
              GENIQUE : tb_report_strain_missing_gene, la liste des genes que
              le pipeline TBannotator a lui-meme juge non fiables pour cette
              souche. Ce script l'utilise comme le fait la piste : pour chaque
              singleton terminal situe dans un CDS, verifier si le gene qui le
              porte est marque « manquant » chez au moins une des 39 autres
              souches du pool -- auquel cas leur « absence de variant » ne peut
              pas etre lue comme « porte H37Rv » et le singleton est SUSPECT.
Entrees     : bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt (panel n=40,
              seed=0, identique a P3.2/P3.3 pour comparabilite)
              data/MTBC0/ (ancestral + masques, P2.4/P3.1)
              investigate_phylo/resources/NC_000962.3.gff3 (locus_tag par CDS)
              TBannotator PostgreSQL (tb_report_strain, tb_report_gene,
              tb_report_strain_missing_gene) via MCP
Sorties     : data/p33_panel_singletons.json (par clade : panel de souches et
              liste des singletons terminaux avec locus_tag), consomme par
              phase12_p33_bilan.py qui croise avec les genes manquants
              TBannotator et ecrit resultats/phase12_p33_*.tsv
Reutilisable: oui -- le geste (position -> locus_tag -> missing_gene comme
              proxy de couverture) vaut pour toute etude MTBC fondee sur
              spdi.txt qui doit statuer sur un singleton
Projet      : GC_par_lignee
Date        : 2026-09-01
"""
import json
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import read_subs, strain_dirs, flux  # noqa: E402
from phase2_polarisation_mtbc0 import build_ancestral, load_mask  # noqa: E402
from phase3_counts_par_souche import MASKS  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
GFF3 = Path("/home/christophe/docs/codes/mtbc/investigate_phylo/resources/"
            "NC_000962.3.gff3")
CLADES = ["L1", "L2.2.1", "L3", "L4.1.2", "L4.3", "L6.1", "L7", "L9",
          "Orygis_La3", "Caprae_La2", "Microti"]
SEED = 0
N_PER_CLADE = 40


def locus_tag_by_position():
    """pos (0-based) -> locus_tag, un CDS par position (le premier rencontre
    en cas de chevauchement, comme load_cds ailleurs dans le projet)."""
    owner = {}
    pat = re.compile(r"locus_tag=([^;]+)")
    for line in GFF3.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        m = pat.search(f[8])
        if not m:
            continue
        tag = m.group(1)
        s, e = int(f[3]) - 1, int(f[4])
        for p in range(s, e):
            if p not in owner:
                owner[p] = tag
    return owner


def build_panel(clade):
    strains = strain_dirs(clade)
    rng = random.Random(SEED)
    rng.shuffle(strains)
    sample = strains[:N_PER_CLADE]
    subsets, names = [], []
    for s in sample:
        v = read_subs(s / "NC_000962.3" / "spdi.txt")
        if v:
            subsets.append(v)
            names.append(s.name)
    return names, subsets


def terminal_singletons(subsets, masked, anc):
    """(pos, ref, alt) -> index de la souche porteuse, pour les singletons non
    masques, avec la classification loss/gain/neutral/inverse/tierce (identique
    a phase3_counts_par_souche)."""
    support = defaultdict(int)
    for i, subs in enumerate(subsets):
        for v in subs:
            support[v] |= 1 << i
    out = {}
    for (pos, ref, alt), mask in support.items():
        if mask & (mask - 1):
            continue
        if pos in masked:
            continue
        i = mask.bit_length() - 1
        a = chr(anc[pos]) if pos < len(anc) else "N"
        if a == "N":
            status = "na"
        elif a == alt:
            status = "inverse"
        elif a != ref:
            status = flux(a, alt)
        else:
            status = flux(ref, alt)
        out[(pos, ref, alt)] = (i, status)
    return out


if __name__ == "__main__":
    anc = build_ancestral()
    masked = load_mask(MASKS, in_mtbc0_coords=True)
    owner = locus_tag_by_position()
    print(f"# masque : {len(masked)} positions ; CDS positions : {len(owner)}",
          file=sys.stderr)

    panels = {}
    for clade in CLADES:
        names, subsets = build_panel(clade)
        if len(subsets) < 4:
            print(f"{clade}\tTROP PEU", file=sys.stderr)
            continue
        term = terminal_singletons(subsets, masked, anc)
        panels[clade] = (names, term)
        print(f"{clade}\tn={len(names)}\tsingletons_non_masques={len(term)}",
              file=sys.stderr)

    payload = {}
    for clade, (names, term) in panels.items():
        rows = []
        for (pos, ref, alt), (i, status) in term.items():
            tag = owner.get(pos)
            rows.append({"pos": pos, "ref": ref, "alt": alt,
                         "carrier": names[i], "status": status, "locus_tag": tag})
        payload[clade] = {"names": names, "singletons": rows}

    out_path = ROOT / "data" / "p33_panel_singletons.json"
    out_path.write_text(json.dumps(payload))
    print(f"# ecrit {out_path} ({sum(len(p['singletons']) for p in payload.values())} "
          f"singletons au total)", file=sys.stderr)
