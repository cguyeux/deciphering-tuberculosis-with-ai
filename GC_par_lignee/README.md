# GC_par_lignee reproducibility bundle

Analysis scripts and derived results underlying:

> Guyeux C. *GC content and mutational GC-to-AT bias across Mycobacterium tuberculosis complex
> lineages.* (in preparation, not yet submitted)

The genomic GC content of the *Mycobacterium tuberculosis* complex (MTBC) is routinely cited as a
fixed, 65% property of the species. This work shows that while genomic GC content itself is
essentially invariant across the seventeen MTBC lineages (differing by at most 70 parts per
million), the underlying Sueoka mutational equilibrium $GC_{eq}=u/(u+v)$ ranges from 32.1% to
55.4% between lineages — a heterogeneity that survives sequencing-platform and coverage-bias
stratification, an independent polarization check, exclusion of repair-deficient strains, a
replichore/strand-of-synthesis control, and restriction to four-fold degenerate sites that rule
out selection at the amino-acid level. Polarized substitution events are counted without tree
reconstruction, oriented against the MTBC0 ancestral genome, from 142,518 surveillance genomes and
a 519-strain, seventeen-lineage panel drawn from a curated local database.

## Contents

- `analyses/` — 36 standard-library-plus-scientific-stack Python scripts, one per investigative
  step (`phaseN_pX_description.py` for numbered analyses, `phaseN_FX_description.py` for the
  scripts generating the manuscript's main and supplementary figures), covering: event counting
  by strain vs. by polarized tree-depth event, polarization against MTBC0 (with an independent
  *M. canettii* outgroup check), coverage-bias and sequencing-platform stratification, the
  Sueoka $GC_{eq}=u/(u+v)$ conversion, hierarchical-bootstrap confidence intervals, sampling-density
  confound tests and density-matched resampling, *nucS*-deletion confound screening, replichore/
  strand-of-synthesis control, four-fold-degenerate-site restriction and codon-degeneracy-class
  selection coefficients, contextual (trinucleotide) spectrum decomposition, and the mutT3/dinP
  DNA-repair inactivation screen.
- `résultats/` — derived/intermediate results (TSV) produced by the `analyses/` scripts: per-strain
  and per-lineage substitution counts, polarized event tables, bootstrap confidence intervals,
  degeneracy-class site classifications and flux ratios, contextual-spectrum tables, and the
  data underlying every figure and table of the manuscript.

## Reproducing the figures and tables

Each `analyses/phase11_FX_*.py` script reproduces one figure (F1–F6 main text, S1–S2
supplementary) from the corresponding `résultats/*.tsv` file. Table 1 ($GC_{eq}$ by lineage) and
Table 2 (selection coefficient by degeneracy class) are produced by the phase5 and phase9 scripts
respectively. Scripts assume Python 3 with numpy, scipy, pandas and matplotlib.

## Data provenance

Source variant data (SPDI calls per strain) were queried from the group's curated local MTBC
surveillance database (TBannotator pipeline output), not distributed here. The `résultats/` files
included are sufficient to regenerate every figure, table and reported statistic without
re-querying the source database.

## Status

This manuscript is not yet submitted to a peer-reviewed venue or deposited as a preprint. This
bundle is released for transparency and reproducibility ahead of that decision.
