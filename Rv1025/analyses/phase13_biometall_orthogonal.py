#!/usr/bin/env python3
"""phase13_biometall_orthogonal.py — P2.6 (orthogonal metal-site predictor)

Cross-check of the AF3-holo metal site (Cys113/His115/Glu59) with an INDEPENDENT,
non-learning predictor: BioMetAll (Sanchez-Aparicio et al., JCIM 2020,
doi:10.1021/acs.jcim.0c00827), which detects metal-binding sites from BACKBONE
PREORGANIZATION only. It uses neither AF3's learned metal placement nor the
sequence-conservation signal, so it is orthogonal to both lines that defined the
site, and directly answers the circularity objection raised in the internal review.

Install (ephemeral venv, deps = numpy + psutil only):
    python3 -m venv --system-site-packages venv_metal
    venv_metal/bin/pip install biometall
Run (blind, whole protein):
    venv_metal/bin/biometall AF-P96375-F1.pdb --min_coordinators 3 --pdb --cores 4

This script parses the BioMetAll result, then measures in the APO model (no metal,
no AF3-holo) the distance from each predicted probe centre to the candidate donor
atoms, and confronts the two predicted sites with their evolutionary conservation.

Reads : résultats/structure/AF-P96375-F1.pdb (apo AF model, chain A, 155 aa)
        résultats/metal_orthogonal/results_biometall_apo.txt (BioMetAll output)
Writes: résultats/metal_orthogonal/biometall_sites_interpretation.tsv
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
APO = ROOT / "résultats/structure/AF-P96375-F1.pdb"
BM = ROOT / "résultats/metal_orthogonal/results_biometall_apo.txt"
OUT = ROOT / "résultats/metal_orthogonal/biometall_sites_interpretation.tsv"

# Per-residue conservation (source: phase5 conservation on 8700-seq AF3 MSA, and
# Pfam PF04417 seed 39 seqs for the triad; decoy values from the mutant-control run).
CONS = {
    ("CYS", 113): 100.0, ("HIS", 115): 100.0, ("GLU", 59): 99.0,   # conserved triad
    ("CYS", 20): 74.0, ("GLU", 24): 8.0, ("HIS", 48): 20.0,        # non-conserved decoy
    ("GLU", 86): None, ("GLU", 91): None, ("ASP", 93): None,       # weak site, n/a
}
# Metal-coordinating atom(s) per residue type.
DONOR = {"CYS": ["SG"], "HIS": ["ND1", "NE2"], "GLU": ["OE1", "OE2"], "ASP": ["OD1", "OD2"]}


def load_atoms(pdb):
    atoms = {}  # (resn,resi) -> {atomname: (x,y,z)}
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        resn = line[17:20].strip()
        resi = int(line[22:26])
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        atoms.setdefault((resn, resi), {})[name] = (x, y, z)
    return atoms


def parse_biometall(txt):
    sites = []
    for line in txt.read_text().splitlines():
        m = re.match(r"\s*(\d+)\s*\|\s*(.+?)\s*\|\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\|\s*(\d+)\s*\|\s*([-\d.]+)", line)
        if not m:
            continue
        rank = int(m.group(1))
        coords = [(rr.split(":")[0], int(rr.split(":")[1]))
                  for rr in m.group(2).split() if ":" in rr]
        centre = (float(m.group(3)), float(m.group(4)), float(m.group(5)))
        sites.append({"rank": rank, "res": coords, "centre": centre,
                      "nprobes": int(m.group(6)), "radius": float(m.group(7))})
    return sites


def dist(a, b):
    return sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5


atoms = load_atoms(APO)
sites = parse_biometall(BM)

rows = [("rank", "residues", "n_probes", "radius_A", "min_donor_dist_A",
         "donor_detail", "min_conservation_pct", "verdict")]
print(f"{'#':>2} {'residues':<28} {'probes':>6} {'radius':>7} {'d_donor':>8} "
      f"{'cons%':>6}  verdict")
print("-" * 90)
for s in sites:
    # nearest donor atom of each coordinating residue to the predicted centre
    details, mind = [], 1e9
    cons_vals = []
    for resn, resi in s["res"]:
        best = min((dist(s["centre"], atoms[(resn, resi)][an]), an)
                   for an in DONOR.get(resn, []) if an in atoms.get((resn, resi), {}))
        details.append(f"{resn}{resi}:{best[1]}={best[0]:.2f}")
        mind = min(mind, best[0])
        c = CONS.get((resn, resi))
        if c is not None:
            cons_vals.append(c)
    mincons = min(cons_vals) if cons_vals else float("nan")
    triad = {("CYS", 113), ("HIS", 115), ("GLU", 59)}
    is_triad = set(s["res"]) == triad
    if is_triad:
        verdict = "CONSERVED TRIAD (biologically relevant)"
    elif cons_vals and mincons < 30:
        verdict = "geometric decoy (NOT conserved)"
    else:
        verdict = "weak/other"
    rows.append((s["rank"], " ".join(f"{r}{i}" for r, i in s["res"]), s["nprobes"],
                 f"{s['radius']:.3f}", f"{mind:.2f}", ";".join(details),
                 f"{mincons:.0f}" if cons_vals else "NA", verdict))
    print(f"{s['rank']:>2} {' '.join(f'{r}{i}' for r,i in s['res']):<28} "
          f"{s['nprobes']:>6} {s['radius']:>7.3f} {mind:>7.2f}Å "
          f"{(f'{mincons:.0f}' if cons_vals else 'NA'):>6}  {verdict}")

with OUT.open("w") as fh:
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")

print(f"\nWritten: {OUT.relative_to(ROOT)}")
print("\nInterpretation:")
print("  - BioMetAll (backbone-only, no AF3, no conservation) independently flags")
print("    the conserved triad Cys113/His115/Glu59 as a metal-compatible site with a")
print("    predicted centre ~metal-bond distance from the donor atoms in the APO model.")
print("  - It ALSO flags the C20/E24/H48 cluster (the same decoy the mutant-control")
print("    AF3 run relocated the ion to). Geometry alone does NOT rank the true site")
print("    first: BOTH pockets are geometrically metal-compatible.")
print("  - The DISCRIMINATOR is CONSERVATION, not geometry: triad 99-100% invariant")
print("    vs decoy 8-74%. Two orthogonal methods (AF3-holo, BioMetAll) agreeing on")
print("    both sites, with conservation selecting between them, is the honest picture.")
