#!/usr/bin/env python3
"""
Parseur RÉUTILISABLE des sorties de la recherche Foldseek web (search.foldseek.com).

Reforgé (2026-07-05) après le piège rencontré : le téléchargement contient, PAR BASE, deux fichiers :
  - `alis_<db>.m8`        = LES VRAIS HITS (format m8 étendu ci-dessous)   <-- à parser
  - `alis_<db>_report.m8` = un RAPPORT DE TAXONOMIE (colonnes pident/rang/taxid/lignée), PAS des hits <-- à IGNORER
Ce script lit uniquement les `alis_<db>.m8` (hors *_report.m8), et sort un tableau propre par base.

Colonnes de `alis_<db>.m8` (observées, foldseek web mode 3diaa) :
  1 query  2 target(+description)  3 pident  4 alnlen  5 mismatch  6 gapopen
  7 qstart 8 qend  9 tstart  10 tend  11 evalue  12 bits  13 ?  14 qlen  15 tlen  16 qaln  17 taln  ...
Couverture query = alnlen / qlen. NB : le format par défaut ne fournit PAS le TM-score (mode tmalign requis) ;
juger la significativité sur evalue + couverture + convergence de plusieurs hits vers un même repli.

Usage : python3 foldseek_parse.py <dir_foldseek_out> [evalue_max=0.01]
"""
import glob, os, sys

def parse_alis(path):
    hits = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        tgt_full = f[1]
        tgt_id = tgt_full.split()[0]
        tgt_desc = tgt_full[len(tgt_id):].strip()
        try:
            pident, alnlen = float(f[2]), int(f[3])
            qstart, qend = int(f[6]), int(f[7])
            evalue, bits = float(f[10]), float(f[11])
            qlen = int(f[13]) if len(f) > 13 and f[13].isdigit() else None
            tlen = int(f[14]) if len(f) > 14 and f[14].isdigit() else None
        except ValueError:
            continue
        cov = (alnlen / qlen) if qlen else None
        hits.append(dict(target=tgt_id, desc=tgt_desc, pident=pident, alnlen=alnlen,
                         qstart=qstart, qend=qend, evalue=evalue, bits=bits,
                         qlen=qlen, tlen=tlen, qcov=cov))
    hits.sort(key=lambda h: h["evalue"])
    return hits

def main():
    if len(sys.argv) < 2:
        print("usage: foldseek_parse.py <dir_foldseek_out> [evalue_max]"); sys.exit(1)
    root = sys.argv[1]
    emax = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
    alis = [p for p in glob.glob(f"{root}/alis_*.m8") if not p.endswith("_report.m8")]
    if not alis:
        print(f"Aucun alis_<db>.m8 (hors _report) sous {root}"); sys.exit(1)
    any_sig = False
    for path in sorted(alis):
        db = os.path.basename(path)[len("alis_"):-len(".m8")]
        hits = parse_alis(path)
        sig = [h for h in hits if h["evalue"] <= emax]
        any_sig = any_sig or bool(sig)
        print(f"\n=== {db} : {len(hits)} hits ({len(sig)} avec e<= {emax}) ===")
        if not hits:
            print("  (aucun)"); continue
        print(f"  {'target':22s} {'pid%':>5s} {'aln':>4s} {'qcov%':>6s} {'eval':>9s} {'bits':>7s}  description")
        for h in hits[:8]:
            cov = f"{h['qcov']*100:5.0f}" if h["qcov"] is not None else "   NA"
            print(f"  {h['target'][:22]:22s} {h['pident']:5.0f} {h['alnlen']:4d} {cov:>6s} "
                  f"{h['evalue']:9.2g} {h['bits']:7.1f}  {h['desc'][:48]}")
    print(f"\nVerdict : {'au moins un hit significatif' if any_sig else 'AUCUN hit significatif = repli sans homologue connu'} "
          f"(seuil e<= {emax}).")
    print("Rappel : convergence de plusieurs hits vers un MÊME repli > un hit isolé ; couverture query faible = match partiel.")

if __name__ == "__main__":
    main()
