#!/usr/bin/env python3
"""phase24_leprae_retention.py — P3.1.a: does retention in M. leprae add information?

The manuscript notes in passing that the eno-divIC-Rv1025-ppx2 block is retained by
the reductive genomes of M. leprae. Surviving a genome that pseudogenised roughly half
its coding capacity is a far stronger functional filter than mere presence — but the
obvious counter-argument, written into the piste before running anything, is
CO-RETENTION: the block contains eno (glycolysis) and two essential genes, so their
retention is expected, and Rv1025 could simply be dragged along.

The test therefore measures what retention is worth GIVEN essentiality:
  * background retention rate of random H37Rv proteins in M. leprae, split by
    essentiality class (from the atlas TnSeq layer, DeJesus 2017);
  * then situates the four operon genes, one of which (ppx2/Rv1026) is NON-essential.
    If non-essential genes are usually lost, the retention of ppx2 alongside the block
    is informative; if they are usually kept, the whole argument is worth little.

Retention = tblastn hit with query coverage >= 80% and E <= 1e-5 (an intact ortholog,
not a fragment) — a pseudogenised or eroded gene fails the coverage criterion.

Reads : Canettii/NC_000962.3_CDS.fasta, bdd/hors_mtbc/M_leprae/ref/genome.fna,
        annotation_mtbc/site/content/genes/<Rv>.json (essentiality)
Writes: résultats/leprae_retention.tsv
"""
import json
import random
import re
import subprocess
import tempfile
from pathlib import Path
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent.parent.parent
CDS = ROOT / "Canettii/NC_000962.3_CDS.fasta"
GENOME = ROOT / "bdd/hors_mtbc/M_leprae/ref/genome.fna"
ATLAS = ROOT / "annotation_mtbc/site/content/genes"
OUT = ROOT / "Rv1025/résultats/leprae_retention.tsv"
OPERON = ["Rv1023", "Rv1024", "Rv1025", "Rv1026"]
N_SAMPLE = 400
COV_MIN, EVAL = 80.0, 1e-5
random.seed(42)


def load_proteins():
    prot, loc, buf = {}, None, []
    for line in CDS.read_text().splitlines():
        if line.startswith(">"):
            if loc and buf:
                prot[loc] = "".join(buf)
            buf = []
            m = re.search(r"\[locus_tag=([^\]]+)\]", line)
            loc = m.group(1) if m else None
        else:
            buf.append(line.strip())
    if loc and buf:
        prot[loc] = "".join(buf)
    out = {}
    for lt, nt in prot.items():
        try:
            p = str(Seq(nt).translate(table=11, to_stop=True))  # type: ignore[arg-type]
        except Exception:  # noqa: BLE001
            continue
        if len(p) >= 50:
            out[lt] = p
    return out


def essential(lt):
    f = ATLAS / f"{lt}.json"
    if not f.exists():
        return None
    e = (json.load(f.open()).get("essentiality") or {})
    return e.get("essential")


prot = load_proteins()
pool = [lt for lt in prot if lt not in OPERON]
sample = random.sample(pool, min(N_SAMPLE, len(pool)))
queries = sample + OPERON
print(f"H37Rv : {len(prot)} protéines ; échantillon {len(sample)} + opéron {len(OPERON)}")

with tempfile.TemporaryDirectory() as td:
    faa = Path(td) / "q.faa"
    faa.write_text("".join(f">{lt}\n{prot[lt]}\n" for lt in queries))
    db = Path(td) / "leprae"
    subprocess.run(["makeblastdb", "-in", str(GENOME), "-dbtype", "nucl", "-out", str(db)],
                   capture_output=True, check=True)
    r = subprocess.run(["tblastn", "-query", str(faa), "-db", str(db), "-evalue", str(EVAL),
                        "-max_target_seqs", "1", "-outfmt", "6 qseqid pident qcovs evalue"],
                       capture_output=True, text=True, check=True)

best = {}
for ln in r.stdout.splitlines():
    q, pid, cov, ev = ln.split("\t")
    cov, ev = float(cov), float(ev)
    if q not in best or cov > best[q][0]:
        best[q] = (cov, float(pid), ev)

retained = {q: (q in best and best[q][0] >= COV_MIN) for q in queries}

# background rates by essentiality class
cls = {lt: essential(lt) for lt in sample}
groups = {"essential": [lt for lt in sample if cls[lt] is True],
          "non-essential": [lt for lt in sample if cls[lt] is False],
          "unknown": [lt for lt in sample if cls[lt] is None]}

rows: list = [("set", "n", "retained", "rate_pct", "note")]
print(f"\n{'classe':<16} {'n':>4} {'retenus':>8} {'taux':>7}")
print("-" * 42)
rates = {}
for name, lts in groups.items():
    if not lts:
        continue
    k = sum(1 for lt in lts if retained[lt])
    rate = 100.0 * k / len(lts)
    rates[name] = k / len(lts)
    rows.append((f"background_{name}", len(lts), k, f"{rate:.1f}", "H37Rv random sample"))
    print(f"{name:<16} {len(lts):>4} {k:>8} {rate:>6.1f}%")

print(f"\n{'gène opéron':<12} {'essentiel':>10} {'couv.':>7} {'%id':>6}  retenu")
print("-" * 52)
for lt in OPERON:
    cov, pid, ev = best.get(lt, (0.0, 0.0, 1.0))
    e = essential(lt)
    print(f"{lt:<12} {str(e):>10} {cov:>6.0f}% {pid:>5.0f}%  {'OUI' if retained[lt] else 'non'}")
    rows.append((lt, 1, int(retained[lt]), f"{cov:.0f}", f"essential={e}; pident={pid:.0f}"))

# the informative question: the NON-essential member of the block
ne_rate = rates.get("non-essential")
if ne_rate is not None:
    print(f"\nTaux de rétention des gènes NON essentiels (fond) : {100*ne_rate:.1f}%")
    print(f"-> ppx2/Rv1026 est NON essentiel et pourtant retenu : probabilité sous le fond = {ne_rate:.2f}")
    print("   (c'est CE membre, pas Rv1025 lui-même, qui rend la rétention du bloc informative :")
    print("    Rv1025 étant essentiel, sa propre rétention est déjà attendue.)")

OUT.write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")
print(f"\nÉcrit : {OUT.relative_to(ROOT)}")
