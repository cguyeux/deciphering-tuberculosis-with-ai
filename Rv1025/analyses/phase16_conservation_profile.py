#!/usr/bin/env python3
"""phase16_conservation_profile.py — P3.1 multi-scale conservation profile of Rv1025.

Consolidates the conservation of Rv1025 (DUF501) at nested scales, from data already
in hand (no new heavy computation; the /challenge control "déjà fait ?" showed the
pieces exist). HONEST FRAMING (KB tuberculosis.md, entries 2026-06-02/06-09): on a
near-invariant gene with a handful of segregating sites, the pN/pS ratio is
underpowered noise; report the NUMBER of segregating sites and their nature, not the
ratio. The atlas fiche labels the intra-MTBC pN/pS 0.581 "relaxed/neutral" — that is
exactly the underpowered verdict; we do NOT propagate it.

Reads : data/Rv1025_atlas_fiche.json      (intra-MTBC: n_strains, snp_sites, syn/missense, disruptive)
        data/Rv1025_canettii_conservation.json  (immediate outgroup)
        data/Rv1025_ntm_conservation.json        (genus Mycobacterium)
Writes: résultats/conservation_profile.tsv
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "data"
OUT = ROOT / "résultats/conservation_profile.tsv"

fiche = json.load(open(D / "Rv1025_atlas_fiche.json"))
can = json.load(open(D / "Rv1025_canettii_conservation.json"))
ntm = json.load(open(D / "Rv1025_ntm_conservation.json"))
c = fiche["conservation"]

rows = [("scale", "sample", "signal", "honest_readout")]

# 1. intra-MTBC — report site count, NOT pN/pS (underpowered)
rows.append((
    "intra-MTBC", f"{c['n_strains']:,} genomes",
    f"{c['snp_sites']} segregating sites ({c['syn']} syn, {c['missense']} missense, "
    f"{c['nonsense']} nonsense, {c['frameshift']} frameshift); {c['disrupt_sites']} disruptive",
    "near-invariant; NOT reported as pN/pS "
    f"(atlas label 'relaxed/neutral' pN/pS={c['pN_pS']} is underpowered on {c['snp_sites']} sites)"))

# 2. M. canettii — immediate outgroup
rows.append((
    "M. canettii (immediate outgroup)", "153 genomes",
    f"{can['n_substitutions']} substitution ({can['n_syn']} syn, {can['n_nonsyn']} non-syn, "
    f"{can['n_disruptive']} disruptive)",
    "invariant at the protein level vs the immediate outgroup; the single change is synonymous "
    "(dN/dS not reliable, low power — reported as a substitution count)"))

# 3. genus Mycobacterium
rows.append((
    "genus Mycobacterium", f"{ntm['n_total']} species",
    f"present in {ntm['n_present']}/{ntm['n_total']} ({ntm['frac']*100:.0f}%), "
    f"mean identity {ntm['mean_pident']:.0f}%; classification '{ntm['classification']}'",
    "genus-core: retained across the whole genus, including the reductive genomes "
    "(M. leprae, M. lepromatosis)"))

# 4. family / phylum (established elsewhere: P3.2 synteny + P2 Pfam)
rows.append((
    "DUF501 family (PF04417)", "4,370 proteins / 5,484 taxa",
    "Actinobacteria-restricted; operon block intact in 61/63 genomes (P3.2); "
    "no GO term, no solved structure",
    "deep family conservation of an ordered gene block (signature of selection on the arrangement)"))

# 5. functional site (triad) — established P2.3/P2.7
rows.append((
    "metal-site triad (Cys113/His115/Glu59)", "8,700-seq MSA + Pfam seed (39)",
    "Cys113 100%, His115 100%, Glu59 99% (MSA) / 100% on curated Pfam seed",
    "the functional site is the most conserved feature at every scale; the 3 intra-MTBC "
    "missense are non-disruptive and fall outside this invariant triad"))

with OUT.open("w") as fh:
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")

print("Multi-scale conservation profile of Rv1025 (DUF501):\n")
for scale, sample, signal, readout in rows[1:]:
    print(f"  [{scale}]  ({sample})")
    print(f"     signal : {signal}")
    print(f"     readout: {readout}\n")
print(f"Written: {OUT.relative_to(ROOT)}")
print("\nNested invariance: intra-MTBC (5 sites, 0 disruptive) -> M. canettii (1 synonymous) ->")
print("genus-core (53/53) -> pan-Actinobacterial family -> invariant metal-site triad.")
