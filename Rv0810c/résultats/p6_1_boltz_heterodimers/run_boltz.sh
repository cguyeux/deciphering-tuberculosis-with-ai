#!/usr/bin/env bash
# P6.1 -- Boltz-2 heterodimeres Rv0810c + VRAIS partenaires STRING (P0.5) : Rv0811c
# (voisin direct, hypothetique) et Rv0812/PabC (4-amino-4-deoxychorismate lyase, voie
# PABA/folate), plus un controle de specificite Rv0810c + RpmG2/Rv0634B (50S ribosomal
# L33, 56 aa, meme controle deja calibre en P3.4, sans raison biologique de s'associer
# a Rv0810c ni a ses partenaires STRING). Les canaux STRING de P0.5 sont
# context_driven/coexpression uniquement (zero experimental) : un signal positif sur
# Rv0811c ou Rv0812 n'est interpretable que s'il depasse nettement ce controle
# (garde-fou #4 du skill boltz).
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
  # NB : pgrep -f matche la ligne de commande COMPLETE. Un entry-point Python (boltz)
  # est execute comme ".../python .../boltz predict ...", donc un ancrage "^${B} predict"
  # ne matche JAMAIS (l'interpreteur precede le chemin de boltz dans argv). Corrige le
  # 2026-08-11 (P6.1) apres avoir observe deux jobs Boltz tourner simultanement malgre ce
  # garde-fou cense l'empecher -- bug herite du gabarit du skill, a corriger aussi dedans.
  other=$(pgrep -f "${B} predict" 2>/dev/null || true)
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
# (constate au lancement de cette piste : job e4_bldc_dimer_dna d'un autre projet actif) --
# ne jamais le tuer, seulement patienter avant de commencer la boucle.
waited=0
while pgrep -f "${B} predict" >/dev/null 2>&1; do
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

# Ordre : controle de specificite d'abord (116 tokens, le plus rapide, valide le
# pipeline), puis Rv0812/PabC (375 tokens), puis Rv0811c (428 tokens, le plus lent).
for spec in "rv0810c_rpmg2_control_heterodimer 116" "rv0810c_rv0812_heterodimer 375" "rv0810c_rv0811c_heterodimer 428"; do
  set -- $spec; y=$1; tokens=$2
  # Reevalue avant CHAQUE job (pas seulement au demarrage), cf. skill boltz.
  waited_job=0
  while pgrep -f "${B} predict" >/dev/null 2>&1; do
    sleep 120
    waited_job=$((waited_job + 1))
    if [ "$waited_job" -gt 180 ]; then
      break
    fi
  done
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
