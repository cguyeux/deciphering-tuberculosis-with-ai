#!/usr/bin/env python3
"""phase17_synteny_null.py — P3.2.a null model: is the operon's synteny atypical?

The manuscript claims the eno-divIC-Rv1025-ppx2 block conservation across the
Actinobacteria is "a signature of selection on the arrangement rather than passive
co-inheritance". This script tests that against a null, on the DISTANT Actinobacteria
(non-Mycobacterium; genome order otherwise extensively rearranged vs H37Rv).

KEY DESIGN (avoids a naive trap): a raw "block conserved?" null conflates ARRANGEMENT
conservation with GENE conservation (Mtb-specific PE/PPE genes would score 0 for the
wrong reason). We therefore CONDITION ON GENE PRESENCE: for each 4-gene block and each
distant genome where ALL 4 orthologs are found, we ask whether they are SYNTENIC
(same contig, small window, consistent H37Rv order). The per-block metric is
    syntenic_fraction = (# genomes syntenic) / (# genomes with all 4 orthologs present).
The operon (syntenic_fraction ~ 1) is then placed in the null distribution of random
4-consecutive-gene blocks. This isolates the arrangement signal from gene conservation.

Reads : Canettii/NC_000962.3_CDS.fasta  (H37Rv CDS, headers carry [locus_tag] [location])
        bdd/hors_mycobacterium/<sp>/ref/genome.fna  (11 distant Actinobacteria)
Writes: résultats/synteny_null/synteny_null.tsv, null_summary.json
Needs  : tblastn, makeblastdb (BLAST+), Biopython
"""
import json
import random
import re
import subprocess
import tempfile
from pathlib import Path
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent.parent.parent          # .../mtbc
CDS = ROOT / "Canettii/NC_000962.3_CDS.fasta"
GDIR = ROOT / "bdd/hors_mycobacterium"
OUT = ROOT / "Rv1025/résultats/synteny_null"
OUT.mkdir(parents=True, exist_ok=True)

# 11 distant Actinobacteria (exclude the 2 out-of-phylum negatives E. coli, B. subtilis)
DISTANT = ["Blongum_NCC2705", "Cacnes_KPA171202", "Cdiphtheriae_NCTC13129",
           "Cglutamicum_ATCC13032", "Gbronchialis_DSM43247", "Mabscessus_ATCC19977",
           "Mluteus_NCTC2665", "Nfarcinica_IFM10152", "Rjostii_RHA1",
           "Scoelicolor_A3-2", "Tpaurometabola_DSM20162"]
OPERON = ["Rv1023", "Rv1024", "Rv1025", "Rv1026"]
N_RANDOM = 90
WINDOW_BP = 8000          # syntenic if all 4 hits fall within this span on one contig
EVAL = 1e-5
random.seed(42)


def parse_cds():
    """Return ordered list of (locus, strand, protein_seq) by genome start coordinate."""
    genes = {}
    loc, seqbuf = None, []
    for line in CDS.read_text().splitlines():
        if line.startswith(">"):
            if loc:
                genes[loc[0]] = (loc[1], loc[2], "".join(seqbuf))
            seqbuf = []
            lt = re.search(r"\[locus_tag=([^\]]+)\]", line)
            m = re.search(r"\[location=(complement\()?(?:join\()?<?(\d+)\.\.>?(\d+)", line)
            if lt and m:
                strand = "-" if m.group(1) else "+"
                loc = (lt.group(1), int(m.group(2)), strand)
            else:
                loc = None
        else:
            seqbuf.append(line.strip())
    if loc:
        genes[loc[0]] = (loc[1], loc[2], "".join(seqbuf))
    # translate, order by start
    prot = {}
    for lt, (start, strand, nt) in genes.items():
        try:
            p = str(Seq(nt).translate(table=11, to_stop=True))  # type: ignore[arg-type]  # NCBI table 11 (int is valid Biopython)
        except Exception:  # noqa: BLE001
            continue
        if len(p) >= 30:
            prot[lt] = (start, strand, p)
    ordered = sorted(prot, key=lambda k: prot[k][0])
    return prot, ordered


def pick_blocks(prot, ordered):
    """Operon + N random 4-consecutive co-directional gene blocks (disjoint, exclude operon)."""
    blocks = {"OPERON": OPERON}
    used = set(OPERON)
    tries = 0
    while len([b for b in blocks if b != "OPERON"]) < N_RANDOM and tries < 5000:
        tries += 1
        i = random.randint(0, len(ordered) - 4)
        win = ordered[i:i + 4]
        if any(g in used for g in win):
            continue
        strands = {prot[g][1] for g in win}
        if len(strands) != 1:            # co-directional, like the operon
            continue
        span = prot[win[-1]][0] - prot[win[0]][0]
        if span > 12000:                 # avoid windows straddling a huge gene/gap
            continue
        blocks[win[0]] = win
        used.update(win)
    return blocks


def run():
    prot, ordered = parse_cds()
    blocks = pick_blocks(prot, ordered)
    print(f"H37Rv genes parsed: {len(prot)} ; blocks: operon + {len(blocks)-1} random")

    # one multi-FASTA with every block's 4 proteins, id = <blockname>__<locus>
    qfaa = OUT / "block_queries.faa"
    with qfaa.open("w") as fh:
        for bname, gs in blocks.items():
            for g in gs:
                fh.write(f">{bname}__{g}\n{prot[g][2]}\n")

    # tblastn all queries vs each distant genome (11 calls total)
    hits = {b: {} for b in blocks}       # block -> genome -> {locus: (contig, mid, ok)}
    with tempfile.TemporaryDirectory() as td:
        for sp in DISTANT:
            fna = GDIR / sp / "ref/genome.fna"
            if not fna.exists():
                print(f"  (absent) {sp}"); continue
            db = Path(td) / sp
            subprocess.run(["makeblastdb", "-in", str(fna), "-dbtype", "nucl",
                            "-out", str(db)], capture_output=True, check=True)
            r = subprocess.run(
                ["tblastn", "-query", str(qfaa), "-db", str(db), "-evalue", str(EVAL),
                 "-max_target_seqs", "1", "-outfmt", "6 qseqid sseqid sstart send evalue pident"],
                capture_output=True, text=True, check=True)
            best = {}   # qseqid -> (contig, mid, evalue)
            for ln in r.stdout.splitlines():
                q, sseq, ss, se, ev, _pid = ln.split("\t")
                ev = float(ev)
                if q not in best or ev < best[q][2]:
                    best[q] = (sseq, (int(ss) + int(se)) // 2, ev)
            for q, (contig, mid, ev) in best.items():
                bname, locus = q.split("__")
                hits[bname].setdefault(sp, {})[locus] = (contig, mid)
            print(f"  tblastn done: {sp} ({len(best)} hits)")

    # score each block: present (all 4 orthologs) and syntenic (1 contig, window, order)
    rows: list = [("block", "genes", "n_present", "n_syntenic", "syntenic_fraction")]
    results = {}
    for bname, gs in blocks.items():
        n_present = n_syn = 0
        for sp, hh in hits[bname].items():
            if len(hh) < 4:
                continue
            n_present += 1
            contigs = {hh[g][0] for g in gs}
            if len(contigs) != 1:
                continue
            mids = [hh[g][1] for g in gs]
            if max(mids) - min(mids) > WINDOW_BP:
                continue
            order = [hh[g][1] for g in gs]           # H37Rv order of the 4 genes
            if order == sorted(order) or order == sorted(order, reverse=True):
                n_syn += 1
        frac = (n_syn / n_present) if n_present else float("nan")
        results[bname] = {"genes": gs, "n_present": n_present, "n_syntenic": n_syn,
                          "syntenic_fraction": frac}
        rows.append((bname, "-".join(gs), n_present, n_syn,
                     f"{frac:.3f}" if n_present else "NA"))

    (OUT / "synteny_null.tsv").write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")

    # operon vs null distribution (random blocks with >=3 genomes present, for a meaningful fraction)
    op = results["OPERON"]
    null = [v["syntenic_fraction"] for b, v in results.items()
            if b != "OPERON" and v["n_present"] >= 3]
    n_ge = sum(1 for x in null if x >= op["syntenic_fraction"])
    pctile = 100 * (1 - n_ge / len(null)) if null else float("nan")
    summary = {
        "operon": op,
        "n_random_blocks_evaluable": len(null),
        "null_mean_syntenic_fraction": round(sum(null) / len(null), 3) if null else None,
        "null_median": round(sorted(null)[len(null) // 2], 3) if null else None,
        "operon_syntenic_fraction": round(op["syntenic_fraction"], 3) if op["n_present"] else None,
        "n_random_ge_operon": n_ge,
        "empirical_p": round((n_ge + 1) / (len(null) + 1), 4) if null else None,
        "operon_percentile": round(pctile, 1) if null else None,
        "distant_genomes": DISTANT,
    }
    (OUT / "null_summary.json").write_text(json.dumps(summary, indent=2))

    print("\n=== P3.2.a null-model result (synteny GIVEN gene presence, distant Actinobacteria) ===")
    print(f"Operon: present in {op['n_present']}/{len(DISTANT)} distant genomes, "
          f"syntenic in {op['n_syntenic']} -> fraction {op['syntenic_fraction']:.2f}")
    if null:
        print(f"Random blocks (n={len(null)} evaluable): mean syntenic_fraction "
              f"{summary['null_mean_syntenic_fraction']}, median {summary['null_median']}")
        print(f"Random blocks with fraction >= operon: {n_ge}/{len(null)}  "
              f"-> empirical p = {summary['empirical_p']}, operon at ~{summary['operon_percentile']}th percentile")
    print(f"\nWritten: {(OUT/'synteny_null.tsv').relative_to(ROOT)} , null_summary.json")


if __name__ == "__main__":
    run()
