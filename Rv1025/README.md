# Rv1025 reproducibility bundle

Analysis scripts and derived data underlying:

> Guyeux C. *A conserved cysteine-histidine-glutamate metal site identifies DUF501 (Rv1025), an
> essential uncharacterised protein family of Mycobacterium tuberculosis, as a candidate
> metalloenzyme and drug target.* Submitted to *Antonie van Leeuwenhoek*.

Rv1025 is a 155-residue essential, vulnerable, uncharacterised *M. tuberculosis* protein carrying
DUF501 (Pfam PF04417), a family with no Gene Ontology term and no solved structure anywhere. A
structural search (Foldseek against PDB100, AlphaFold-DB/SwissProt and CATH50) finds no homolog,
indicating a novel fold. Conservation across 8,700 homologous sequences reveals a near-invariant
Cys113-His115-Glu59 triad; holo AlphaFold3 modelling places a divalent metal (Zn/Fe/Mn) on this
triad at ~2.3 A, made specific by a triad-mutant control that expels the ion and corroborated
blind by an independent geometry predictor (BioMetAll). The site is universal across 1,472
near-complete bacterial genomes of the family. The conserved *eno*-*divIC*-*Rv1025*-*ppx2* operon
(61/63 actinobacterial genomes) does not correspond to a stable Rv1025-DivIC complex, nor to a
complex with five other divisome components (AlphaFold3 and Boltz-2, positive control
DivIC-FtsQ recovered). This bundle reproduces every figure and every reported statistic in the
manuscript.

## Contents

- `analyses/` — 39 standard-library-plus-scientific-stack Python scripts (plus two PyMOL `.pml`
  scenes and one polling shell script), one per investigative step
  (`phaseN_description.py`), covering: operon synteny and its null control, AlphaFold3/Boltz-2
  complex modelling and parsing for the divisome panel (DivIC, FtsQ, FtsZ, FtsK, FtsW, PBP-b/FtsI,
  SepF), metal-site geometry and holo modelling, conservation logo and weighted conservation
  across the DUF501 family, BioMetAll orthogonal corroboration, druggability and pocket detection,
  direct-coupling analysis (DCA) and invariant-cluster interpretation, homodimer modelling, an
  AlphaFold-bias control, *M. leprae* retention, and the eukaryotic branch of DUF501 (apicomplexan
  tandem-domain architecture, rumen-fungi horizontal transfer). `phase10_synteny_figure.py` and
  `phase11_afmultimer_figure.py` generate Figures 2 and 3; the active-site close-up (Figure 1) is
  rendered headless by the PyMOL scene `phase8_active_site_figure.pml`.
- `data/` — input data: the curated DUF501 alignment (Pfam PF04417, `PF04417_full.sto`), the
  Rv1025 sequence (`Rv1025.faa`), reference metalloprotein structures used as geometric controls
  (carbonic anhydrase 1CA2, zinc-finger 1ZNF, thermolysin 8TLN, matrix metalloproteinase 4MDT,
  AlphaFold models of human ADH1B/P00325-family and catalase P00918), genomic context of the
  rumen-fungi horizontal-transfer candidates, and the atlas fiche and dossier for Rv1025 (with
  dated seeds tracing provenance).
- `résultats/` — derived/intermediate results (JSON, TSV, PDB/CIF, npz) produced by the
  `analyses/` scripts, including cached AlphaFold3 and Boltz-2 structural predictions (apo, holo
  with each candidate metal, triad-mutant controls, the full divisome interaction panel, homodimer
  jobs, apicomplexan triad positions), DCA coupling matrices, synteny and its permutation null, and
  druggability/pocket-detection output. Raw MSA search caches (`msas/`, `msa_cache/`, `msa/`
  subdirectories under the AlphaFold3/Boltz-2 job folders) are excluded as regenerable,
  non-informative bulk (multi-megabyte `.a3m` files); everything needed to re-derive the reported
  statistics and figures from the cached model outputs is included.

## Reproducing the figures

- Figure 1 (active-site close-up, metal triad geometry): PyMOL scene
  `analyses/phase8_active_site_figure.pml`, run headless with `pymol -cq phase8_active_site_figure.pml`.
- Figure 2 (*eno*-*divIC*-*Rv1025*-*ppx2* operon synteny): `python3 analyses/phase10_synteny_figure.py`.
- Figure 3 (AlphaFold3/Boltz-2 divisome interaction panel): `python3 analyses/phase11_afmultimer_figure.py`.

## Reproducing the reported statistics

Each `analyses/phaseN_*.py` script is self-contained and reads from `data/` and/or the cached
`résultats/` produced by an earlier phase, writing its result back to `résultats/` (JSON, TSV or
npz). Scripts assume Python 3 with numpy, scipy, pandas, biopython and matplotlib; the AlphaFold3
and Boltz-2 job scripts (`phase2_afmultimer_jobs.py`, `phase6_metal_holo_jobs.py`,
`phase19_homodimer_jobs.py`, `phase33_divisome_panel_jobs.py`, `phase21_boltz_homodimer.py`)
require the corresponding local or remote inference setup and are provided for full pipeline
transparency; their cached outputs under `résultats/` let every downstream statistic and figure be
regenerated without repeating GPU- or API-bound inference.

## Data provenance

All external data are public and cited with accession numbers, DOIs or supplementary-data sources
in the manuscript's Materials and Methods. The DUF501 family alignment (1,472 near-complete
bacterial sequences plus the eukaryotic branch) was assembled from Pfam PF04417 and public genome
assemblies; the full per-sequence provenance is in `résultats/duf501_taxonomy.tsv` and
`résultats/duf501_family_triad.tsv`.

## Citation

If you use this bundle, please cite the accompanying manuscript and, where relevant, this
repository's Zenodo archive (DOI in the manuscript's Data, Metadata, and Code Availability
section).
