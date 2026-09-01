#!/usr/bin/env python3
"""phase19_homodimer_jobs.py — P4.4: is the Rv1025 metal site completed IN TRANS?

The holo AF3 runs (P2.4) coordinate the ion with only THREE protein ligands
(Cys113-Sgamma, His115-Ndelta1, Glu59-Oepsilon), leaving open coordination
positions. We read that as a catalytic, solvent-exposed metal. The competing and
untested explanation is mundane in metalloenzymology: the free positions are filled
by ligands from a SECOND SUBUNIT (an interfacial metal site of a homodimer).

All five AF-Multimer jobs run so far (P4.1) were HETEROdimers (operon neighbours +
controls); the most obvious partner, the protein itself, was never tested. This
script writes the two jobs that settle it:

  A) Rv1025 x2 + 2 Fe   -> if a confident interface forms AND chain B donors reach
                           the ion, the site is shared/interfacial (structuring
                           result). If each ion stays in its own protomer on the
                           same three ligands, the site is NOT completed in trans
                           (negative -> strengthens the current catalytic reading,
                           which is so far asserted without a control).
  B) Rv1025 x2, apo     -> control: does the same interface form WITHOUT the metal?
                           Guards against a metal-driven docking artefact.

Same sequence and seed as the holo jobs, so results are directly comparable.
Read-out (KB guard, same bar as P4.1): ipTM + minimum inter-chain PAE +
reproducibility across the five models, AND check that the interface actually puts
donor atoms in contact with the ion, not merely that an interface exists.

Reads : résultats/af3_metal_out/fold_rv1025_holo_fe/..._job_request.json (sequence + seed)
Writes: résultats/phase19_homodimer/fold_rv1025_homodimer_{2fe,apo}_job_request.json
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "résultats/af3_metal_out/fold_rv1025_holo_fe/fold_rv1025_holo_fe_job_request.json"
OUT = ROOT / "résultats/phase19_homodimer"
LIGANDS = {59: "E", 113: "C", 115: "H"}

ref = json.load(REF.open())[0]
seq = ref["sequences"][0]["proteinChain"]["sequence"]
seed = ref["modelSeeds"]
assert len(seq) == 155, f"unexpected length {len(seq)}"
for pos, aa in LIGANDS.items():
    assert seq[pos - 1] == aa, f"residue {pos} is {seq[pos-1]}, expected {aa}"

OUT.mkdir(parents=True, exist_ok=True)


def chain(count):
    return {"proteinChain": {"sequence": seq, "count": count, "useStructureTemplate": True}}


jobs = {
    "fold_rv1025_homodimer_2fe": [{
        "name": "Rv1025_homodimer_2Fe",
        "modelSeeds": seed,
        "sequences": [chain(2), {"ion": {"ion": "FE", "count": 2}}],
        "dialect": "alphafoldserver", "version": 3,
    }],
    "fold_rv1025_homodimer_apo": [{
        "name": "Rv1025_homodimer_apo",
        "modelSeeds": seed,
        "sequences": [chain(2)],
        "dialect": "alphafoldserver", "version": 3,
    }],
}
for name, job in jobs.items():
    (OUT / f"{name}_job_request.json").write_text(json.dumps(job, indent=1))
    print(f"  written: {(OUT / (name + '_job_request.json')).relative_to(ROOT)}")

print(f"\nSequence: {len(seq)} aa, ligands verified (Glu59/Cys113/His115), seed {seed[0]}"
      " (identical to the holo monomer jobs -> directly comparable).")
print("\nHAND-OFF: submit both JSONs on the AlphaFold Server, deposit the result zips in")
print(f"  {OUT.relative_to(ROOT)}/ , then parse with phase20 (interface + ion coordination).")
print("\nDecision rule set IN ADVANCE (so the result cannot be rationalised after the fact):")
print("  * SHARED SITE  : ipTM confident and reproducible across the 5 models, low inter-chain")
print("    PAE at the interface, AND >=1 donor atom from chain B within ~2.5 A of an ion.")
print("  * NOT SHARED   : each ion keeps the same 3 intra-protomer ligands, no chain-B donor")
print("    in contact -> the open positions face solvent/substrate, as currently argued.")
print("  * ARTEFACT     : interface only in the holo job and absent in the apo control, or")
print("    ipTM collapsing across models (cf. Rv1025-DivIC: 0.25 -> 0.08) -> discard.")
