#!/usr/bin/env python3
"""phase29_duf501_eukaryotes.py — P8.1 (suite) : la triade est-elle perdue chez les EUCARYOTES ?

phase28 a établi que la triade Cys-His-Glu est présente chez 96,1 % des 2 228 membres de
PF04417, avec ZÉRO substitution non conservative : les 86 exceptions sont toutes des gaps.
Mais 75 de ces 86 sont des fragments (longueur médiane 76 résidus contre 124), tandis que
**11 sont quasi complètes** — et leurs mnémoniques sont eucaryotes, avec un motif partagé
« E-- » (glutamate gardé, paire Cys-His perdue).

Ce script teste si c'est une VRAIE sous-famille eucaryote sans site métal, ou un mirage.
Trois façons dont ce résultat pourrait être faux, testées explicitement :

  1. CONFONDANT DE COMPLÉTUDE. Les eucaryotes pourraient simplement être moins bien
     séquencés/annotés. Parade : restreindre l'analyse aux séquences QUASI COMPLÈTES et
     refaire le test dessus.
  2. TAXONOMIE DEVINÉE. Les mnémoniques UniProt ne sont pas une taxonomie. Parade :
     interroger UniProt pour la lignée réelle de chaque accession.
  3. HITS HMM MARGINAUX. Un domaine Pfam trouvé dans une grande protéine eucaryote à
     une position aberrante peut être un faux positif du modèle. Parade : mesurer
     l'identité à Rv1025 des membres eucaryotes contre celle des membres bactériens.
     Si les eucaryotes sont nettement en dessous, l'homologie elle-même est douteuse et
     la « perte » ne veut rien dire.

Lit   : data/PF04417_full.sto  (+ cache résultats/duf501_taxonomy.tsv)
Écrit : résultats/duf501_taxonomy.tsv, résultats/duf501_eukaryotes.txt,
        article/supplementary_materials/table_S10_duf501_eukaryotes.tsv
"""
import os
import subprocess
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fast_blas import ensure  # noqa: E402  # type: ignore[import-not-found]

ensure()

from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.stats import fisher_exact, mannwhitneyu  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
STO = ROOT / "data/PF04417_full.sto"
TAX = ROOT / "résultats/duf501_taxonomy.tsv"
OUT = ROOT / "résultats"
SUPP = ROOT / "article/supplementary_materials"
COLS = {59: 172, 113: 395, 115: 397}       # colonnes ancrées sur Rv1025 (phase28)
EXPECT = {59: "E", 113: "C", 115: "H"}
ANCHOR = "P96375_MYCTU"
BATCH = 50

out = []


def say(s=""):
    print(s)
    out.append(s)


# ------------------------------------------------------------------ alignement
seqs: dict = defaultdict(str)
order: list = []
acc: dict = {}
for line in STO.read_text().splitlines():
    if line.startswith("#=GS") and " AC " in line:
        parts = line.split()
        acc[parts[1]] = parts[3].split(".")[0]
        continue
    if not line or line.startswith("#") or line.startswith("//"):
        continue
    parts = line.split()
    if len(parts) != 2:
        continue
    if parts[0] not in seqs:
        order.append(parts[0])
    seqs[parts[0]] += parts[1]

length = {n: sum(1 for c in seqs[n] if c not in "-.") for n in order}
triad = {n: "".join(seqs[n][c].upper() if seqs[n][c] not in "-." else "-"
                    for c in COLS.values()) for n in order}
complete = {n: all(seqs[n][c].upper() == EXPECT[p] for p, c in COLS.items()) for n in order}
med = float(np.median([length[n] for n in order if complete[n]]))
near_full = {n for n in order if length[n] >= 0.9 * med}
say(f"{len(order)} séquences ; longueur médiane des membres à triade complète : {med:.0f}")
say(f"{len(near_full)} séquences QUASI COMPLÈTES (>= 90 % de cette médiane)")

# ------------------------------------------------------- taxonomie (UniProt)
if TAX.exists():
    lineage = dict(l.split("\t")[:2] for l in TAX.read_text().splitlines()[1:] if "\t" in l)
    say(f"Taxonomie lue depuis le cache ({len(lineage)} accessions)")
else:
    lineage = {}
    accs = sorted({acc[n] for n in order if n in acc})
    say(f"Interrogation d'UniProt pour {len(accs)} accessions ({-(-len(accs)//BATCH)} lots) ...")
    for i in range(0, len(accs), BATCH):
        chunk = accs[i:i + BATCH]
        q = "+OR+".join(f"accession:{a}" for a in chunk)
        url = (f"https://rest.uniprot.org/uniprotkb/search?query={q}"
               f"&fields=accession,lineage,protein_name&format=tsv&size=500")
        r = subprocess.run(["curl", "-sS", "--max-time", "120", url],
                           capture_output=True, text=True)
        for line in r.stdout.splitlines()[1:]:
            cells = line.split("\t")
            if len(cells) < 2:
                continue
            a, lin = cells[0], cells[1]
            name = cells[2] if len(cells) > 2 else ""
            # Une part des entrées de l'alignement Pfam a été SUPPRIMÉE d'UniProt depuis :
            # lignée vide + protein_name « deleted ». Le dire explicitement vaut mieux que
            # « inconnu », qui laisserait croire à un échec de requête.
            dom = "supprimé d'UniProt" if (not lin.strip() or name.strip() == "deleted") else "inconnu"
            for token in ("Eukaryota", "Archaea", "Bacteria"):
                if f"{token} (domain)" in lin:
                    dom = token
                    break
            lineage[a] = dom
    for a in {acc[n] for n in order if n in acc} - set(lineage):
        lineage[a] = "non retourné par UniProt"
    TAX.write_text("accession\tdomaine\n"
                   + "\n".join(f"{a}\t{d}" for a, d in sorted(lineage.items())) + "\n")
    say(f"Écrit : {TAX.relative_to(ROOT)} ({len(lineage)} accessions)")

dom = {n: lineage.get(acc.get(n, ""), "inconnu") for n in order}

# ------------------------------- A. domaine taxonomique x triade, sur les QUASI COMPLÈTES
say()
say("=" * 78)
say("A. TRIADE PAR DOMAINE DU VIVANT (séquences quasi complètes seulement)")
say("=" * 78)
say(f"{'domaine':<12} {'n':>6} {'triade complète':>18}")
sub = [n for n in order if n in near_full]
for d, cnt in Counter(dom[n] for n in sub).most_common():
    ok = sum(1 for n in sub if dom[n] == d and complete[n])
    say(f"  {d:<10} {cnt:>6} {ok:>10}/{cnt} = {ok/cnt:6.1%}")
euk = [n for n in sub if dom[n] == "Eukaryota"]
bac = [n for n in sub if dom[n] == "Bacteria"]
if euk and bac:
    table = [[sum(complete[n] for n in bac), sum(not complete[n] for n in bac)],
             [sum(complete[n] for n in euk), sum(not complete[n] for n in euk)]]
    odds, pv = fisher_exact(table)
    say(f"Fisher bactéries vs eucaryotes : OR = {odds:.1f}, p = {pv:.3g}")
    say(f"  -> le confondant de complétude est ÉCARTÉ : le test ne porte que sur des")
    say(f"     séquences de longueur comparable.")

# ------------------------------------------ B. quel motif chez les eucaryotes ?
say()
say("=" * 78)
say("B. QUE RESTE-T-IL DU SITE CHEZ LES EUCARYOTES ?")
say("=" * 78)
say(f"motifs observés (E59/C113/H115) sur {len(euk)} eucaryotes quasi complets :")
for pat, n in Counter(triad[x] for x in euk).most_common():
    say(f"  {pat} : {n}")
say(f"motifs chez les {len(bac)} bactéries quasi complètes :")
for pat, n in Counter(triad[x] for x in bac).most_common(4):
    say(f"  {pat} : {n}")

# --------------------------- C. ces hits eucaryotes sont-ils de vrais homologues ?
say()
say("=" * 78)
say("C. HOMOLOGIE RÉELLE OU HIT HMM MARGINAL ?")
say("=" * 78)
ref = seqs[next(n for n in order if n.startswith(ANCHOR))]


def ident_to_ref(n):
    s = seqs[n]
    both = [(a.upper(), b.upper()) for a, b in zip(ref, s) if a not in "-." and b not in "-."]
    return 100.0 * sum(a == b for a, b in both) / len(both) if both else 0.0


ie = np.array([ident_to_ref(n) for n in euk])
ib = np.array([ident_to_ref(n) for n in bac])
say(f"identité à Rv1025 : eucaryotes {ie.mean():.1f} % (médiane {np.median(ie):.1f}, "
    f"min {ie.min():.1f}, max {ie.max():.1f})")
say(f"                    bactéries  {ib.mean():.1f} % (médiane {np.median(ib):.1f})")
u = mannwhitneyu(ie, ib, alternative="less")
say(f"Mann-Whitney (eucaryotes MOINS identiques) : p = {u.pvalue:.3g}")  # type: ignore[attr-defined]
say("Seuil de prudence : sous ~20 % d'identité sur un domaine court, un hit HMM n'est pas")
say("une homologie assurée, et une « perte de site » n'y est pas interprétable.")
frac_low = float((ie < 20).mean())
say(f"  eucaryotes sous 20 % d'identité : {frac_low:.0%}")

rows = [("sequence", "accession", "domaine", "longueur", "motif_E59_C113_H115",
         "identite_Rv1025_%")]
for n in sorted(euk, key=lambda x: -ident_to_ref(x)):
    rows.append((n, acc.get(n, "?"), dom[n], str(length[n]), triad[n], f"{ident_to_ref(n):.1f}"))
(SUPP / "table_S10_duf501_eukaryotes.tsv").write_text(
    "\n".join("\t".join(r) for r in rows) + "\n")
say()
say(f"Détail des {len(euk)} membres eucaryotes : table_S10_duf501_eukaryotes.tsv")
(OUT / "duf501_eukaryotes.txt").write_text("\n".join(out) + "\n")
