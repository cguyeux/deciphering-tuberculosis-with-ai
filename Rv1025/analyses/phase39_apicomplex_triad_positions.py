#!/usr/bin/env python3
"""phase39_apicomplex_triad_positions.py — P9.4 : positions exactes des résidus de
triade (Glu59/Cys113/His115, numérotation ancrée sur Rv1025) dans les deux copies
PF04417 en tandem de Babesia ovis (A0A9W5T8J3) et Besnoitia besnoiti (A0A2A9M7T6),
par le même ancrage de colonnes que phase28 (alignement Pfam full déjà en main).
Nécessaire pour préparer les jobs Boltz (protéine entière) et mesurer, comme pour
Neospora caninum (F0VAI8, P8.1.a.4), le regroupement spatial du triplet reconstitué
par duplication complémentaire.

Lit   : data/PF04417_full.sto
Écrit : résultats/phase39_apicomplex_triad_positions.tsv
"""
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STO = ROOT / "data/PF04417_full.sto"
OUT = ROOT / "résultats/phase39_apicomplex_triad_positions.tsv"

ANCHOR = "P96375_MYCTU"
TRIAD = {59: "E", 113: "C", 115: "H"}
TARGETS = {"A0A9W5T8J3": "Babesia ovis", "A0A2A9M7T6": "Besnoitia besnoiti"}

seqs: dict = defaultdict(str)
order: list = []
for line in STO.read_text().splitlines():
    if not line or line.startswith("#") or line.startswith("//"):
        continue
    parts = line.split()
    if len(parts) != 2:
        continue
    name, frag = parts
    if name not in seqs:
        order.append(name)
    seqs[name] += frag

anchor_name = next(n for n in order if n.startswith(ANCHOR))
start = int(anchor_name.split("/")[1].split("-")[0])
row = seqs[anchor_name]
col_of = {}
pos = start - 1
for c, ch in enumerate(row):
    if ch not in "-.":
        pos += 1
        if pos in TRIAD:
            col_of[pos] = c
for p, aa in TRIAD.items():
    assert row[col_of[p]].upper() == aa

rows = [("accession_organisme", "label_sto", "domaine_start_end",
         "res59_aa", "res59_pos", "res113_aa", "res113_pos", "res115_aa", "res115_pos")]
for acc, organism in TARGETS.items():
    labels = sorted(n for n in order if n.split("_")[0] == acc)
    for label in labels:
        s = seqs[label]
        dom_start, dom_end = label.split("/")[1].split("-")
        dom_start = int(dom_start)
        out_res = {}
        p = dom_start - 1
        for c, ch in enumerate(s):
            if ch not in "-.":
                p += 1
                for pcol in TRIAD:
                    if c == col_of[pcol]:
                        out_res[pcol] = (ch.upper(), p)
        r = [f"{acc} ({organism})", label, f"{dom_start}-{dom_end}"]
        for pcol in (59, 113, 115):
            aa, seqpos = out_res.get(pcol, ("-", None))
            r += [aa, seqpos if seqpos is not None else ""]
        rows.append(tuple(r))

with open(OUT, "w") as fh:
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")

for r in rows:
    print("\t".join(str(x) for x in r))
