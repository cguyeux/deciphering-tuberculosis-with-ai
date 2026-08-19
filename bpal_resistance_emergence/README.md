# Analysis scripts and result tables for: Genomic emergence of resistance to the BPaL/BPaLM regimen across the *Mycobacterium tuberculosis* complex

Christophe Guyeux, FEMTO-ST Institute (UMR 6174 CNRS), Université Marie et Louis Pasteur, Besançon, France.

## Contents

- `analyses/` — the phased Python scripts (`phase1_...py` to `phase6_...py`) that produced
  every number, table and figure reported in the manuscript: feasibility scan and
  trichotomy, IPW-adjusted lineage prevalence, primed-polymorphism / ESM-1v scoring,
  DST validation, convergence discovery and its permutation null, candidate emergence
  timing (raw and lineage/country-stratified), the structural figure, and the
  compensatory-background negative control. `paths.py` documents every external input
  path the scripts consume (see "Upstream dependency" below).
- `résultats/` — the corresponding output tables (TSV), including
  `phase5c_candidate_carriers.tsv`, the per-candidate carrier table (accession, year,
  country, major lineage/sub-lineage, RIF-R status) for every convergent F420-pocket
  candidate discussed in the manuscript.
- `data/` — the 15-gene resistance panel (`gene_panel.tsv`) and the curated
  literature-mining evidence tables (`litcohort_*`) used for the DST-cohort corroboration
  in the Discussion.

## Upstream dependency (not bundled here)

This project is a **consumer**, not a producer, of a shared antimicrobial-resistance
reference infrastructure maintained by a sibling internal project
(`Resistance_antibio/`): the mutation-to-resistance catalogue (WHO 2nd ed. 2023 +
tb-profiler + CRyPTIC, 62,020 graded assertions), the 132,609-strain MTBC variant
pangenome, the phenotype-consensus table, and the Coll-scheme lineage classification of
248,771 strains. These resources are read-only inputs (see `analyses/paths.py`) and are
**not redistributed in this deposit**: they are substantially heavier, belong to a
separate research infrastructure with its own publication timeline, and are not yet
independently deposited. Anyone wishing to re-run these scripts from scratch needs
access to that infrastructure; this deposit instead makes the analysis logic and the
per-candidate result tables auditable and citable independently of that dependency.

## Also not included

`résultats/litcohort_cache/` (a local cache of mined full-text/JATS content per PMID,
used transiently to build `data/litcohort_manual_evidence.tsv`) is excluded: it is a
working cache of third-party copyrighted article text, not an analysis result.

## Reproducing a figure or table

Each `analyses/phaseN_*.py` script reads its inputs (local `résultats/`/`data/` files
and/or the upstream `Resistance_antibio/` paths above) and writes one or more TSV files
into `résultats/`, matching the files provided here. Scripts are numbered in the order
they were run; later phases consume tables produced by earlier ones.

## License

CC-BY 4.0 (code and data tables).
