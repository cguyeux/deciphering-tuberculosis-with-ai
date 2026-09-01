#!/usr/bin/env python3
# pyright: reportAttributeAccessIssue=false, reportOptionalMemberAccess=false, reportOptionalSubscript=false
# (Bio.PDB stubs type Structure/Entity iteration loosely: `atom.element`, `residue.resname`
#  and `structure[0]` are the documented idiom but are flagged. Code verified at runtime;
#  suppressing here rather than contorting correct Biopython usage — cf. phase17.)
"""phase20_homodimer_parse.py — P4.4 read-out: is the metal site completed IN TRANS?

Companion to phase19 (job generation). Answers the ONE question phase3 cannot:
for each ion, do the coordinating donor atoms come from ONE chain (site internal to
a protomer, as in the holo monomer) or from TWO chains (interfacial, shared site)?

Interface quality (ipTM, pTM, inter-chain PAE, contact count) is NOT re-implemented
here: run the existing tool on the same directory, which applies the same bar as P4.1:
    python3 analyses/phase3_afmultimer_parse.py résultats/phase19_homodimer

Decision rule (fixed in phase19, before seeing the result):
  SHARED   -> confident reproducible interface AND >=1 donor from the other chain
              within COORD_CUT of an ion.
  NOT SHARED -> every ion keeps its 3 intra-protomer ligands (Glu59/Cys113/His115),
              no cross-chain donor: the open positions face solvent/substrate.
  ARTEFACT -> interface present only in the holo job (absent in the apo control) or
              ipTM collapsing across models.

Reads : résultats/phase19_homodimer/**/fold_rv1025_homodimer_*_model_*.cif
Writes: résultats/phase19_homodimer/homodimer_metal_readout.tsv
"""
from pathlib import Path
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch

ROOT = Path(__file__).resolve().parent.parent
DIR = ROOT / "résultats/phase19_homodimer"
OUT = DIR / "homodimer_metal_readout.tsv"
COORD_CUT = 2.8      # metal-ligand bond (same cut as phase7)
SHELL_CUT = 3.4
TRIAD = {59, 113, 115}


def load(cif):
    s = MMCIFParser(QUIET=True).get_structure("x", str(cif))[0]
    metals, prot = [], []
    for ch in s:
        for res in ch:
            rn = res.resname.strip()
            if rn in ("ZN", "FE", "FE2", "FE3", "MN"):
                metals += [(ch.id, a) for a in res]
            elif res.id[0] == " ":
                prot += [(ch.id, res, a) for a in res if a.element != "H"]
    return metals, prot


rows: list = [("model", "ion", "ion_chain", "n_donor_chains", "cross_chain_donor",
               "donors", "verdict")]
print(f"{'model':<34} {'ion':>4} {'chains':>7}  donors")
print("-" * 96)
for cif in sorted(DIR.rglob("fold_rv1025_homodimer_*_model_*.cif")):
    metals, prot = load(cif)
    if not metals:
        print(f"{cif.stem:<34}  (apo: no ion — interface control, see phase3)")
        continue
    atoms = [a for _, _, a in prot]
    idx = {id(a): (cid, res) for cid, res, a in prot}
    ns = NeighborSearch(atoms)
    for i, (mchain, m) in enumerate(metals, 1):
        donors, chains = [], set()
        for a in ns.search(m.coord, SHELL_CUT):
            if a.element not in ("N", "O", "S"):
                continue
            cid, res = idx[id(a)]
            d = m - a
            if d <= COORD_CUT:
                donors.append(f"{cid}/{res.resname.strip()}{res.id[1]}:{a.name}={d:.2f}")
                chains.add(cid)
        cross = len(chains) > 1
        triad_only = all(int(x.split("/")[1].split(":")[0][3:]) in TRIAD
                         for x in donors) if donors else False
        verdict = ("SHARED (cross-chain donor)" if cross else
                   ("intra-protomer triad" if triad_only else "intra-protomer (other)"))
        donors.sort(key=lambda s: float(s.split("=")[1]))
        rows.append((cif.stem, f"{m.get_parent().resname.strip()}{i}", mchain,
                     len(chains), "YES" if cross else "no",
                     ";".join(donors) or "none", verdict))
        print(f"{cif.stem:<34} {i:>4} {len(chains):>7}  {';'.join(donors[:4]) or 'none'}"
              f"   <- {verdict}")

with OUT.open("w") as fh:
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")
print(f"\nWritten: {OUT.relative_to(ROOT)}")
print("Reminder: pair this with `phase3_afmultimer_parse.py résultats/phase19_homodimer`")
print("for ipTM / inter-chain PAE / reproducibility across the five models.")
