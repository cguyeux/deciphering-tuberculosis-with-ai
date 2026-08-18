#!/usr/bin/env bash
# P3.4 -- Rv0810c homodimere (hypothese tete-beche) + controle de specificite
# RpmG2/Rv0634B (50S ribosomal L33, 55 aa, taille comparable, ne s'auto-associe pas
# hors ribosome -- garde-fou #4 du skill boltz : Boltz replie volontiers N'IMPORTE
# QUEL petit domaine globulaire en homodimere symetrique confiant, artefact deja
# observe sur Rv2516c. Le score de Rv0810c n'est interpretable QUE s'il depasse
# nettement ce controle.
set -u
trap '' HUP
cd "$(dirname "$0")"
B=$HOME/venvs/boltz/bin/boltz
STATUS=boltz_status.log

if [ ! -x "$B" ]; then
  echo "ERREUR : binaire Boltz introuvable ou non executable : $B" | tee -a "$STATUS" >&2
  exit 1
fi

MEM_C=5367
MEM_A=0.02639

check_resources() {
  local tokens="$1" other avail_mib need_mib
  other=$(pgrep -f "^${B} predict" 2>/dev/null || true)
  if [ -n "$other" ]; then
    echo "  PRE-VOL ECHEC : un autre processus Boltz tourne deja (PID $other)." >&2
    return 1
  fi
  avail_mib=$(free -m | awk '/^Mem:/{print $7}')
  need_mib=$(awk -v n="$tokens" -v c="$MEM_C" -v a="$MEM_A" 'BEGIN{printf "%d", c + a*n*n}')
  if [ "${avail_mib:-0}" -lt $((need_mib * 6 / 5)) ]; then
    echo "  PRE-VOL ECHEC : RAM disponible ${avail_mib:-?} Mio < 1,2x besoin estime ${need_mib} Mio" >&2
    return 1
  fi
  echo "  pre-vol OK : ${avail_mib} Mio dispo, besoin estime ${need_mib} Mio ($tokens tokens)"
  return 0
}

# Attente courtoise d'un job Boltz d'une AUTRE session deja en cours au demarrage
# (constate : e2_unitB_dimer_dna_brouille, PID actif au lancement de cette piste) --
# ne jamais le tuer, seulement patienter avant de commencer la boucle.
waited=0
while pgrep -f "^${B} predict" >/dev/null 2>&1; do
  if [ "$waited" -eq 0 ]; then
    echo "$(date +%H:%M:%S) : un autre job Boltz tourne deja, attente..." | tee -a "$STATUS"
  fi
  sleep 120
  waited=$((waited + 1))
  if [ "$waited" -gt 180 ]; then  # 6h de garde, au-dela on abandonne l'attente
    echo "$(date +%H:%M:%S) : attente > 6h, abandon de l'attente initiale" | tee -a "$STATUS"
    break
  fi
done

for spec in "rv0810c_homodimer 120" "rpmg2_control_homodimer 110"; do
  set -- $spec; y=$1; tokens=$2
  if ! check_resources "$tokens"; then
    echo "  $y SAUTE (pre-vol ressources echoue, $(date +%H:%M:%S)) -- a relancer manuellement plus tard" | tee -a "$STATUS"
    continue
  fi
  echo "=== $y $(date +%H:%M:%S) ===" | tee -a "$STATUS"
  "$B" predict "$y.yaml" --out_dir "out_$y" --accelerator cpu --devices 1 \
      --output_format mmcif --diffusion_samples 5 --num_workers 0 --use_msa_server \
      >> boltz.log 2>&1 && echo "  $y OK $(date +%H:%M:%S)" | tee -a "$STATUS" || echo "  $y ECHEC $(date +%H:%M:%S)" | tee -a "$STATUS"
done
echo "=== TERMINE $(date +%H:%M:%S) ===" | tee -a "$STATUS"
