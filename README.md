# Deciphering Tuberculosis with AI

An open, AI-assisted research program on the *Mycobacterium tuberculosis* complex (MTBC), developed in the group of Christophe Guyeux at FEMTO-ST (CNRS UMR 6174, Université de Franche-Comté, Besançon, France). The program combines phylogenomics, population genetics and host-pathogen evolution, with an analysis environment heavily augmented by AI.

This page is a curated index. New items are added one at a time, as each reaches sufficient maturity.

Last updated: 2026-08-26 (grouped manuscripts into thematic categories: lineage phylogenomics and evolution, gene characterisation, antimicrobial resistance).

## About this repository

The software and online tools listed here are released openly. The reference datasets, when published, are provisional, dated and versioned, and may be revised. The manuscripts are preprints or working drafts pending peer review. None of the scientific claims should be cited as established results before peer-reviewed publication. We welcome collaborators willing to help validate, complete and submit this work.

## Applications and web tools

| Tool | Description | Access | Status |
|------|-------------|--------|--------|
| MTBC gene annotation atlas | A successor to the discontinued Mycobrowser reference database, providing curated, AI-assisted functional annotations of MTBC genes, with a focus on requalifying the large stock of genes still labelled as hypothetical proteins. | [Live site](https://mtbc127bab2a-mtbc.functions.fnc.fr-par.scw.cloud) | Online, in development |
| MTBC lineage atlas | A wiki-like dictionary of the taxonomic diversity of the complex, with one normalised fiche per lineage and sub-lineage (phylogeny, geography, dating, drug resistance, selection, codivergence, host and epidemiology), a clickable phylogenetic tree and a world map. A parallel publication channel for the dozens of sub-lineages characterised in this program. | [Live site](https://mtbc127bab2a-atlas-lineage.functions.fnc.fr-par.scw.cloud) | Online, in development |
| MTBC drug-resistance tester | A companion tool that predicts anti-tuberculosis drug resistance from a strain's variants, using the WHO catalogue of mutations (2nd ed., 2023) consolidated with the tb-profiler database and CRyPTIC empirical signal. Variant matching is reconciled across representations (single-nucleotide, multi-nucleotide and indel) to avoid systematic under-calling. A research tool, not a clinical diagnostic. | [Live site](https://mtbc127bab2a-mtbc.functions.fnc.fr-par.scw.cloud/resistance) | Online, in development |

## Manuscripts

### Lineage phylogenomics and evolution

| Title | Description | Access | Status |
|-------|-------------|--------|--------|
| The re-humanisation of L6/L9/L10: opportunistic re-colonisation of *H. sapiens* by an animal-adapted lineage, without genomic reversion | The MTBC lineages L6, L9 and L10 branch within the animal-associated part of the complex yet are found almost exclusively in humans today. The manuscript argues that they re-colonised *Homo sapiens* opportunistically, without reversing the genomic signatures inherited from their animal-adapted ancestry. | [PDF](manuscripts/rehumanisation-L6L9L10.pdf) | Draft, pending peer review |
| Characterisation of L4.13, a cryptic European sub-lineage of *Mycobacterium tuberculosis* with candidate positive selection on accessory ESX-2/4 type VII secretion modules | A phylogenomic characterisation of L4.13, a European L4 sub-lineage defined by Freschi et al. (2021) but absent from the Coll, Stucki and Napier barcoding schemes and undiscussed since. From 258 genomes it proposes validated SNP markers for its bipartite structure (L4.13.1 and L4.13.2), places its origin in Western Europe, and reports candidate positive selection on the accessory ESX-2/4 type VII secretion modules. | [PDF](manuscripts/L4.13-characterisation.pdf) | Draft, pending peer review |

### Gene characterisation (dark genome)

| Title | Description | Access | Status |
|-------|-------------|--------|--------|
| Rv3222c, the unnamed gene between *sigH* and *rshA*, encodes an essential intrinsically disordered protein and is the last intact member of a degenerate gene family in *Mycobacterium tuberculosis* | Rv3222c surfaced from an atlas-wide vulnerability screen as an uncharacterised hypothetical protein physically inserted (4-bp overlaps on both flanks) between the *sigH*-*rshA* oxidative-stress operon. Systematic re-derivation of every inherited claim shows it is essential (Tn-seq, robust to mappability and polarity objections, with a phylogenetic pattern of unequal requirement across lineages), intrinsically disordered (44.8%), under strong purifying selection, has no detectable sequence or structural homolog, and is the last coding survivor of a small gene family degraded across the genus by an ancestral IS110-family insertion. An existing PROSPECT chemogenomic depletion strain makes it experimentally tractable despite the absence of a druggability claim. | [PDF](manuscripts/Rv3222c-characterisation.pdf) | Submitted to Journal of Bacteriology, pending peer review |
| Rv0810c: a genus-conserved, structurally ordered small protein of unknown function carrying DUF3073 in *Mycobacterium tuberculosis* | Rv0810c is a 60-residue protein carrying the unknown-function domain DUF3073 (Pfam PF11273), flagged as an Actinobacteria-signature protein in 2006 and never studied since. This manuscript establishes that it is genuinely translated (detected in 11 of 16 proteomic datasets), has no significant human homologue, no overlap with a neighbouring gene, a bipartite architecture (a rigid conserved module followed by an intrinsically disordered acidic tail), strong purifying selection, and conservation without a single confirmed loss across 260 Actinomycetia genera. Eight independent computational strategies converge on the same negative: no assignable fold, binding site, or functional neighbour. Presented as a rigorous negative result rather than a manufactured function. | [PDF](manuscripts/Rv0810c-characterisation.pdf) | Submitted to Journal of Bacteriology, pending peer review |

### Antimicrobial resistance

| Title | Description | Access | Status |
|-------|-------------|--------|--------|
| Genomic emergence of resistance to the BPaL/BPaLM regimen across the *Mycobacterium tuberculosis* complex: standing pocket variation versus drug-era determinants | The WHO-recommended BPaLM regimen depends on two drugs whose resistance is poorly catalogued, pretomanid especially. Using a consolidated mutation catalogue and a 132,609-strain variant pangenome, this manuscript maps where, when and in which lineages BPaLM resistance arises, discovers F420-pathway pretomanid-resistance candidates phenotype-free by phylogenetic convergence (triangulated with structural localisation in the F420 cofactor pocket and ESM-1v scoring), and shows they are standing polymorphisms that predate the drugs rather than drug-era-selected determinants. Literature mining independently corroborates the F420 pathway as a drug-inducible resistance route and surfaces a fourth, DST-anchored candidate; a compensatory-background analysis serves as a negative control. | [PDF](manuscripts/bpal-resistance-emergence.pdf) | Submitted to Antimicrobial Agents and Chemotherapy, pending peer review |

## Code and reproducibility bundles

| Project | Contents | Access |
|---------|----------|--------|
| Rv0810c | Analysis scripts and derived data reproducing every figure and every reported statistic of the Rv0810c manuscript (DUF3073, *M. tuberculosis*). | [Rv0810c/](Rv0810c/) |
| bpal_resistance_emergence | Phased analysis pipeline and result tables reproducing the numbers, tables and figures of the BPaL/BPaLM resistance-emergence manuscript (feasibility scan, IPW-adjusted lineage prevalence, primed-polymorphism/ESM-1v scoring, DST validation, convergence discovery with permutation null, candidate emergence timing, structural figure, compensatory-background negative control), plus the curated literature-mining evidence tables. Lightweight, project-specific deposit: excludes the shared upstream AMR reference infrastructure (mutation catalogue, pangenome, phenotype table), which belongs to a separate internal project. | [bpal_resistance_emergence/](bpal_resistance_emergence/) |

## Contributing

If you would like to help validate, finish or co-author any of these works, please open an issue on this repository or contact Christophe Guyeux (cguyeux@femto-st.fr). Microbiologists, clinicians and population geneticists are especially welcome.

## Acknowledgements

This project would never have come to life without the many discussions and reflections shared with Guislaine Refrégier, Christophe Sola, Clément Lecarpentier and Gaétan Senelle.
