#!/usr/bin/env python3
"""
P1.7 — ncRv10810c comme element du contexte transcriptionnel de Rv0810c.
Deux tests independants des trois lectures (leader/5'UTR, sRNA trans, artefact) :
(1) conservation ADN de la fenetre Rv0810c-Rv0811c chez M. bovis / M. marinum / M. smegmatis
    (blastn local, genomes en cache dans le depot) ;
(2) repliement du transcrit de 51 nt (ncRv10810c/candidate_1074, 905113-905163) et
    z-score contre 1000 permutations mononucleotidiques (ViennaRNA doit etre installe :
    venv utilise cette seance, cf. experiments/2026-08-10_P1_7_ncRv10810c/).

Necessite RNA (ViennaRNA python bindings) — installe dans un venv dedie si absent
du systeme (PEP 668 sur Arch empeche pip install --system).
"""
import json
import random
import statistics
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

REPO = Path(__file__).resolve().parent.parent
H37RV_FASTA = REPO.parent / "investigate_phylo" / "resources" / "NC_000962.3.fasta"
BOVIS_FASTA = REPO.parent / "investigate_phylo" / "resources" / "LT708304.1.fasta"
MARINUM_FASTA = REPO.parent / "bdd" / "hors_mtbc" / "M_marinum" / "ref" / "genome.fna"
SMEGMATIS_FASTA = REPO / "experiments" / "2026-08-10_P1_7_ncRv10810c" / "M_smegmatis_NC_008596.1.fasta"

WORKDIR = REPO / "experiments" / "2026-08-10_P1_7_ncRv10810c"
OUT_JSON = REPO / "résultats" / "p1_7_ncrv10810c_context.json"

# coordonnees H37Rv 1-based
WINDOW_START, WINDOW_END = 904950, 905350          # Rv0810c CDS (fin) + IGR 146nt + debut Rv0811c
CORE_START, CORE_END = 905113, 905163              # ncRv10810c / candidate_1074 (Miotto 2012, Ami 2020)


def run(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd), file=sys.stderr)
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def blast_window(genome_fasta, db_prefix, evalue="1e-5"):
    query = WORKDIR / "query_window.fasta"
    run(["makeblastdb", "-in", str(genome_fasta), "-dbtype", "nucl", "-out", str(db_prefix)])
    r = run(
        [
            "blastn",
            "-query",
            str(query),
            "-db",
            str(db_prefix),
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-evalue",
            evalue,
        ]
    )
    hits = []
    for line in r.stdout.strip().splitlines():
        cols = line.split("\t")
        hits.append(
            dict(
                zip(
                    ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                     "qstart", "qend", "sstart", "send", "evalue", "bitscore"],
                    cols,
                )
            )
        )
    return hits


def mfe_shuffle_zscore(seq_rna, n=1000, seed=12345):
    import RNA

    def mfe(s):
        fc = RNA.fold_compound(s)
        structure, e = fc.mfe()
        return structure, e

    real_structure, real_mfe = mfe(seq_rna)
    random.seed(seed)
    chars = list(seq_rna)
    shuffles = []
    for _ in range(n):
        random.shuffle(chars)
        _, e = mfe("".join(chars))
        shuffles.append(e)
    mean_s = statistics.mean(shuffles)
    sd_s = statistics.stdev(shuffles)
    z = (real_mfe - mean_s) / sd_s
    n_as_or_more_stable = sum(1 for x in shuffles if x <= real_mfe)
    return {
        "structure": real_structure,
        "mfe_kcal_mol": real_mfe,
        "shuffle_n": n,
        "shuffle_mfe_mean": mean_s,
        "shuffle_mfe_sd": sd_s,
        "z_score": z,
        "empirical_p_le": n_as_or_more_stable / n,
    }


def main():
    WORKDIR.mkdir(parents=True, exist_ok=True)

    rec = next(SeqIO.parse(str(H37RV_FASTA), "fasta"))
    window = rec.seq[WINDOW_START - 1 : WINDOW_END]
    core_plus = rec.seq[CORE_START - 1 : CORE_END]
    core_transcription_sense = str(Seq(str(core_plus)).reverse_complement())  # gene sur brin -

    query_fa = WORKDIR / "query_window.fasta"
    query_fa.write_text(f">Rv0810c_Rv0811c_IGR_window_{WINDOW_START}_{WINDOW_END}\n{window}\n")

    genomes = {
        "M_bovis_AF2122_97": (BOVIS_FASTA, WORKDIR / "bovis_db"),
        "M_marinum_E11": (MARINUM_FASTA, WORKDIR / "marinum_db"),
        "M_smegmatis_MC2_155": (SMEGMATIS_FASTA, WORKDIR / "smegmatis_db"),
    }

    conservation = {}
    for name, (fasta, db) in genomes.items():
        hits = blast_window(fasta, db)
        best = max(hits, key=lambda h: float(h["bitscore"])) if hits else None
        conservation[name] = {"n_hits_evalue_1e-5": len(hits), "best_hit": best}

    seq_rna = core_transcription_sense.replace("T", "U")
    fold_result = mfe_shuffle_zscore(seq_rna)

    out = {
        "window_coords_H37Rv_1based": [WINDOW_START, WINDOW_END],
        "core_ncRv10810c_coords_1based": [CORE_START, CORE_END],
        "core_sequence_transcription_sense_RNA": seq_rna,
        "conservation_blastn_window_vs_genome": conservation,
        "rna_fold_vs_mononucleotide_shuffle": fold_result,
    }

    OUT_JSON.parent.mkdir(exist_ok=True)
    OUT_JSON.write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
