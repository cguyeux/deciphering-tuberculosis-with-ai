#!/usr/bin/env python3
"""phase21_boltz_homodimer.py — P4.4 run LOCALLY with Boltz-2 (no AlphaFold Server hand-off).

Why a local run at all: AlphaFold Server has no public API, its SPA cannot be driven
by the browser extension (KB, 2026-07-05), and it requires a Google login plus
accepting terms of use — actions the assistant must not perform. Boltz-2 is
open-weights, needs no account, and handles multimers AND ions, so the P4.4 question
(is the metal site completed in trans?) can be screened locally.

Scope discipline (decided with CG, 2026-07-31): Boltz-2 is used for SCREENING /
hypothesis testing. AlphaFold3 Server remains the reference method for numbers that
go into the manuscript, so published values stay same-model comparable. Any Boltz
result reported must be labelled as such.

Setup actually used (reproducible):
    uv python install 3.11                       # deps pin numpy<2, scipy==1.13.1 -> no 3.14 wheels
    uv venv --python 3.11 ~/venvs/boltz
    uv pip install --python ~/venvs/boltz/bin/python \
        --index-url https://download.pytorch.org/whl/cpu torch   # CPU-only index ONLY (KB: else ~3 GB of nvidia-*)
    uv pip install --python ~/venvs/boltz/bin/python boltz

MSA: reuses the 8,700-sequence alignment AlphaFold3 already produced for Rv1025
(résultats/af3_out/.../msas/..._unpaired_msa_chains_a.a3m). This avoids sending the
sequence to an external MSA server and gives Boltz the same evolutionary input as AF3.

Writes: résultats/phase21_boltz_homodimer/{homodimer_2fe,homodimer_apo}.yaml
        + run_boltz.sh (the exact commands)
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "résultats/af3_metal_out/fold_rv1025_holo_fe/fold_rv1025_holo_fe_job_request.json"
MSA = ROOT / "résultats/af3_out/fold_rv1025_divic_test/msas/fold_rv1025_divic_test_unpaired_msa_chains_a.a3m"
OUT = ROOT / "résultats/phase21_boltz_homodimer"
VENV = Path.home() / "venvs/boltz/bin/boltz"
LIGANDS = {59: "E", 113: "C", 115: "H"}

seq = json.load(REF.open())[0]["sequences"][0]["proteinChain"]["sequence"]
assert len(seq) == 155
for pos, aa in LIGANDS.items():
    assert seq[pos - 1] == aa, f"residue {pos} is {seq[pos-1]}, expected {aa}"
assert MSA.exists(), f"MSA absent: {MSA}"
with MSA.open() as fh:
    fh.readline()
    assert fh.readline().strip().startswith(seq[:20]), "MSA query != Rv1025 sequence"

OUT.mkdir(parents=True, exist_ok=True)

HOLO = f"""version: 1
sequences:
  - protein:
      id: [A, B]
      sequence: {seq}
      msa: {MSA}
  - ligand:
      id: C
      ccd: FE
  - ligand:
      id: D
      ccd: FE
"""
APO = f"""version: 1
sequences:
  - protein:
      id: [A, B]
      sequence: {seq}
      msa: {MSA}
"""
(OUT / "homodimer_2fe.yaml").write_text(HOLO)
(OUT / "homodimer_apo.yaml").write_text(APO)

run = f"""#!/usr/bin/env bash
# P4.4 local run (Boltz-2, CPU). Long: expect hours on 16 cores for a 310-residue dimer.
# Weights (~1-2 GB) download into ~/.boltz on first call.
set -euo pipefail
BOLTZ="{VENV}"
OUT="{OUT}"
for job in homodimer_2fe homodimer_apo; do
  echo "=== $job ==="
  "$BOLTZ" predict "$OUT/$job.yaml" \\
      --out_dir "$OUT" \\
      --accelerator cpu --devices 1 \\
      --output_format mmcif \\
      --diffusion_samples 1 \\
      --num_workers 4
done
"""
(OUT / "run_boltz.sh").write_text(run)
(OUT / "run_boltz.sh").chmod(0o755)

print(f"Inputs written in {OUT.relative_to(ROOT)}/:")
print("  homodimer_2fe.yaml  (test: 2 chains + 2 Fe)")
print("  homodimer_apo.yaml  (control: 2 chains, no metal)")
print("  run_boltz.sh        (exact commands, CPU)")
print(f"\nMSA reused: {MSA.name} ({8700} sequences, from the AF3 run) -> no external MSA call.")
print("\nRead-out: same phase20 parser (cross-chain donors) + phase3 (interface metrics);")
print("Boltz writes mmCIF, so both parsers apply unchanged.")
