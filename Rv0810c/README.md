# Rv0810c reproducibility bundle

Analysis scripts and derived data underlying:

> Guyeux C. *Rv0810c: a genus-conserved, structurally ordered small protein of unknown function
> carrying DUF3073 in Mycobacterium tuberculosis.* (submitted)

Rv0810c is a 60-residue *M. tuberculosis* protein carrying the unknown-function domain DUF3073
(Pfam PF11273). It is genuinely translated, has no significant human homologue, no overlap with a
neighbouring gene, and is conserved without a single confirmed loss across 260 Actinomycetia
genera, yet eight independent computational strategies (sequence, structure, electrostatic-patch,
embedding-similarity and homo-oligomerisation searches) converge on the same negative: no
assignable fold, binding site, or functional neighbour. This bundle reproduces every figure and
every reported statistic in the manuscript.

## Contents

- `analyses/` — 31 standard-library-plus-scientific-stack Python scripts, one per investigative
  step (`phaseN_pX_Y_description.py`), covering: proteomic detection controls, structural
  architecture (AlphaFold pLDDT profile, radius of gyration, electrostatic patch), DUF3073
  conservation across Actinomycetia (HHsearch, MSA), essentiality and CRISPR-interference
  vulnerability, McDonald-Kreitman selection tests, population variant structure and ESM-1v
  variant scoring, kinase-motif and phosphosite analysis, homo-/hetero-oligomerisation with
  Boltz-2, pocket detection, and the macrophage-detection contradiction (ESX secretion motif,
  IEDB epitope overlap, transcriptomic compendium). `phase30_fig_architecture_conservation.py`
  generates the two manuscript figures.
- `data/` — input data: the H37Rv reference annotation and gene table, the AlphaFold model of
  Rv0810c, the curated DUF3073 alignment (Pfam PF11273, `PF11273_full.sto`), RefSeq assembly
  summaries, the normalised transcriptomic compendium (`log_tpm_norm.csv`), the missense-variant
  and synonymous/non-synonymous control catalogues, and the atlas fiche for Rv0810c.
- `résultats/` — derived/intermediate results (JSON, TSV) produced by the `analyses/` scripts,
  including cached outputs of computationally expensive steps (Boltz-2 homo- and
  hetero-oligomer predictions, HHsearch, pocket detection, ESM Atlas queries, IEDB epitope
  scan) so that downstream statistics and figures can be regenerated without repeating GPU- or
  API-bound computation.

## Reproducing the figures

Figure 1 (bipartite architecture and pLDDT profile) and Figure 2 (DUF3073 conservation across
Actinomycetia) are both produced by:

```
python3 analyses/phase30_fig_architecture_conservation.py
```

which reads `data/AF-I6XWB9-F1-model_v6.pdb` and `résultats/p3_3_msa_conservation.json` (itself
produced by `analyses/phase16_p3_3_msa_conservation.py` from `data/PF11273_full.sto`).

## Reproducing the reported statistics

Each `analyses/phaseN_*.py` script is self-contained and writes its result to the matching
`résultats/pX_*.json` (or `.tsv`) file consumed elsewhere in the pipeline or cited directly in the
manuscript text. Scripts assume Python 3 with numpy, scipy, pandas, biopython and matplotlib.

## Data provenance

All external data are public and cited with accession numbers, DOIs or supplementary-data sources
in the manuscript's Materials and Methods. Population variant data were queried from a curated
MTBC variant database not distributed here; the derived catalogues needed to reproduce the
selection and variant-structure statistics (`data/p8_2_missense_catalogue_84.csv`,
`data/p8_ns_s_sites_75_controls.csv`) are included.

## Citation

If you use this bundle, please cite the accompanying manuscript and, where relevant, this
repository's Zenodo archive (DOI in the manuscript's Data, Metadata, and Code Availability
section).
