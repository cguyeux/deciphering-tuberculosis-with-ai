#!/usr/bin/env bash
# Poll le job Foldseek web, télécharge et décompresse les résultats m8 quand COMPLETE.
set -u
TID="0JewIei9FAdNwmaLOlzhfGojnuoXNeQjqrex5g"
OUT="/home/christophe/docs/codes/mtbc/Rv1025/résultats/structure/foldseek_out"
mkdir -p "$OUT"
for i in $(seq 1 90); do   # ~15 min max (90 x 10 s)
  st=$(curl -s "https://search.foldseek.com/api/ticket/$TID" | python3 -c "import sys,json;print(json.load(sys.stdin).get('status','?'))" 2>/dev/null)
  echo "[$i] status=$st"
  if [ "$st" = "COMPLETE" ]; then
    curl -s -o "$OUT/results.tar.gz" "https://search.foldseek.com/api/result/download/$TID"
    tar -xzf "$OUT/results.tar.gz" -C "$OUT" 2>/dev/null
    echo "TÉLÉCHARGÉ ; fichiers :"; ls -la "$OUT"
    exit 0
  fi
  if [ "$st" = "ERROR" ] || [ -z "$st" ] || [ "$st" = "?" ]; then
    echo "ERREUR ou statut inconnu ($st)"; exit 1
  fi
  sleep 10
done
echo "TIMEOUT après 15 min"; exit 2
