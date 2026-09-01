#!/usr/bin/env bash
# Runner Boltz auto-défensif, gabarit canonique du skill `boltz` (révision 2026-08-17),
# adapté pour ce lot combiné : P4.3 (panel divisome, 5 paires + contrôle positif) et
# P8.1.a.4 (monomère apicomplexe, test aveugle du site métal).
#
# Correction P4.3 apportée ici : les 6 YAML du panel avaient chacun leur chaîne B en
# "msa: empty" (repli individuel dégradé, asymétrie consciente consignée au cahier
# 2026-08-01). --use_msa_server est disponible dans ce sandbox (curl direct vérifié
# 2026-08-19, cf. cahier). PREMIER ESSAI (échoué, découvert par un "OK" trompeur du
# gabarit -- le job sepF s'est terminé en quelques secondes en "OK" sans AUCUN fichier
# de sortie) : `--use_msa_server` global + chaîne A à MSA explicite (Rv1025 réutilisé) +
# chaîne B sans champ msa -> Boltz refuse net ("Cannot mix custom and auto-generated MSAs
# in the same input!"), traite l'entrée comme invalide et la SAUTE, code de sortie 0 --
# indiscernable d'un succès pour `wait`. CORRIGÉ : MSA du partenaire calculé À L'AVANCE via
# `boltz.data.msa.mmseqs2.run_mmseqs2()` (appel direct, pas de `boltz predict`), écrit en
# .a3m dans msa_cache/, et pointé en `msa:` explicite -- les DEUX chaînes sont désormais
# "custom", plus de --use_msa_server pour ces jobs. Le job apicomplexe (P8.1.a.4) reste
# seul avec --use_msa_server : une seule chaîne, pas de mélange possible.
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
  other=$(pgrep -f "${B} predict" 2>/dev/null || true)
  if [ -n "$other" ]; then
    echo "  PRE-VOL ECHEC : un autre processus Boltz tourne deja (PID $other)." >&2
    return 1
  fi
  avail_mib=$(free -m | awk '/^Mem:/{print $7}')
  need_mib=$(awk -v n="$tokens" -v c="$MEM_C" -v a="$MEM_A" 'BEGIN{printf "%d", c + a*n*n}')
  if [ "${avail_mib:-0}" -lt $((need_mib * 6 / 5)) ]; then
    echo "  PRE-VOL ECHEC : RAM disponible ${avail_mib:-?} Mio < 1,2x besoin estime ${need_mib} Mio ($tokens tokens)" >&2
    return 1
  fi
  echo "  pre-vol OK : ${avail_mib} Mio dispo, besoin estime ${need_mib} Mio ($tokens tokens)"
  return 0
}

MAX_CYCLES=150
CYCLE_SLEEP=300
STALL_TIMEOUT=1800

run_job() {
  local y="$1" yaml="$2" outdir="$3" use_server="$4" extra=()
  echo "=== $y $(date +%H:%M:%S) ===" | tee -a "$STATUS"
  [ "$use_server" = "server" ] && extra=(--use_msa_server)
  "$B" predict "$yaml" --out_dir "$outdir" --accelerator cpu --devices 1 \
      --output_format mmcif --diffusion_samples 1 --num_workers 0 "${extra[@]}" \
      >> boltz.log 2>&1 &
  local bpid=$! last_size=-1 last_change
  last_change=$(date +%s)
  while kill -0 "$bpid" 2>/dev/null; do
    sleep 30
    local cur_size now
    cur_size=$(stat -c%s boltz.log 2>/dev/null || echo -1)
    now=$(date +%s)
    if [ "$cur_size" != "$last_size" ]; then
      last_size=$cur_size
      last_change=$now
    elif [ $((now - last_change)) -ge "$STALL_TIMEOUT" ]; then
      echo "  $y BLOQUE : boltz.log immobile depuis ${STALL_TIMEOUT}s, process $bpid toujours vivant -> SIGKILL" | tee -a "$STATUS"
      kill -KILL "$bpid" 2>/dev/null
      wait "$bpid" 2>/dev/null
      echo "  $y ECHEC (stall)" | tee -a "$STATUS"
      return 1
    fi
  done
  wait "$bpid" && echo "  $y OK" | tee -a "$STATUS" || echo "  $y ECHEC" | tee -a "$STATUS"
}

is_done() {
  local y="$1" outdir="$2"
  find "$outdir" -name "confidence_${y}_model_0.json" 2>/dev/null | grep -q .
}

R=/home/christophe/docs/codes/mtbc/Rv1025/résultats
# "<nom_job>:<yaml>:<out_dir>:<n_tokens>:<msa_mode>" -- msa_mode "custom" = MSA explicite
# deja fournie pour TOUTES les chaines du YAML (pas de --use_msa_server, cf. incident
# "Cannot mix" ci-dessus) ; "server" = --use_msa_server (chaine(s) sans champ msa).
jobs="\
posctrl_divIC_ftsQ:${R}/phase33_divisome_panel/posctrl_divIC_ftsQ.yaml:${R}/phase33_divisome_panel/out_posctrl_divIC_ftsQ:542:custom \
rv1025_ftsZ_Rv2150c:${R}/phase33_divisome_panel/rv1025_ftsZ_Rv2150c.yaml:${R}/phase33_divisome_panel/out_rv1025_ftsZ_Rv2150c:534:custom \
rv1025_sepF_Rv2147c:${R}/phase33_divisome_panel/rv1025_sepF_Rv2147c.yaml:${R}/phase33_divisome_panel/out_rv1025_sepF_Rv2147c:373:custom \
apicomplex_neospora_F0VAI8:${R}/phase36_apicomplexan/apicomplex_neospora_F0VAI8.yaml:${R}/phase36_apicomplexan/out_apicomplex_neospora_F0VAI8:472:server \
rv1025_ftsW_Rv2154c:${R}/phase33_divisome_panel/rv1025_ftsW_Rv2154c.yaml:${R}/phase33_divisome_panel/out_rv1025_ftsW_Rv2154c:679:custom \
rv1025_pbpB_ftsI_Rv2163c:${R}/phase33_divisome_panel/rv1025_pbpB_ftsI_Rv2163c.yaml:${R}/phase33_divisome_panel/out_rv1025_pbpB_ftsI_Rv2163c:834:custom \
rv1025_ftsK_Rv2748c:${R}/phase33_divisome_panel/rv1025_ftsK_Rv2748c.yaml:${R}/phase33_divisome_panel/out_rv1025_ftsK_Rv2748c:986:custom"
# ordre : du plus leger au plus lourd (posctrl d'abord pour validation, puis petit -> gros)
# afin qu'un maximum de jobs passent des la premiere fenetre de memoire disponible.

remaining="$jobs"
for cycle in $(seq 1 $MAX_CYCLES); do
  [ -z "$remaining" ] && break
  next_remaining=""
  for spec in $remaining; do
    y="${spec%%:*}"; rest="${spec#*:}"
    yaml="${rest%%:*}"; rest="${rest#*:}"
    outdir="${rest%%:*}"; rest="${rest#*:}"
    tokens="${rest%%:*}"; msa_mode="${rest##*:}"
    if is_done "$y" "$outdir"; then
      continue
    fi
    if check_resources "$tokens"; then
      run_job "$y" "$yaml" "$outdir" "$msa_mode"
    else
      next_remaining="$next_remaining $spec"
    fi
  done
  remaining="$next_remaining"
  if [ -n "$remaining" ]; then
    echo "  cycle $cycle/$MAX_CYCLES : $(echo $remaining | wc -w) job(s) restant(s), attente ${CYCLE_SLEEP}s ($(date +%H:%M:%S))" | tee -a "$STATUS"
    sleep "$CYCLE_SLEEP"
  fi
done
if [ -n "$remaining" ]; then
  echo "  TIMEOUT : $(echo $remaining | wc -w) job(s) jamais lances apres $MAX_CYCLES cycles : $remaining" | tee -a "$STATUS"
fi
echo "=== TERMINE $(date +%H:%M:%S) ===" | tee -a "$STATUS"
