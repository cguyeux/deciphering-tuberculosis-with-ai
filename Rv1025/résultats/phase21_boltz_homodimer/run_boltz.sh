#!/usr/bin/env bash
# P4.4 local run (Boltz-2, CPU). Long: expect hours on 16 cores for a 310-residue dimer.
# Weights (~1-2 GB) download into ~/.boltz on first call.
set -euo pipefail
BOLTZ="/home/christophe/venvs/boltz/bin/boltz"
OUT="/home/christophe/docs/codes/mtbc/Rv1025/résultats/phase21_boltz_homodimer"
for job in homodimer_2fe homodimer_apo; do
  echo "=== $job ==="
  "$BOLTZ" predict "$OUT/$job.yaml" \
      --out_dir "$OUT" \
      --accelerator cpu --devices 1 \
      --output_format mmcif \
      --diffusion_samples 1 \
      --num_workers 4
done
