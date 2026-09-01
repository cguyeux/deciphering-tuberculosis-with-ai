#!/usr/bin/env python3
"""
P3.2 - Test de conservation de synténie de l'opéron eno(Rv1023)-divIC(Rv1024)-Rv1025-ppx2(Rv1026).

Question : l'opéron est-il conservé, et l'adjacence divIC-Rv1025 tient-elle à travers le genre
Mycobacterium et l'outgroup M. canettii ? (P3.2)

Méthode : tblastn des 4 protéines H37Rv contre chaque génome de référence, meilleur hit par requête,
puis analyse de l'ordre et des distances autour de l'ancre Rv1025 (DUF501, spécifique Actinobactéries,
donc hit propre). eno/ppx2 étant des familles répandues, on ne compte comme "syntenique" que la copie
PROCHE de l'ancre.

Entrées : bdd/hors_mtbc/<sp>/ref/genome.fna, Canettii/References/{CP007299.1 (canettii), NC_002945.3 (bovis)}.
Sortie : résultats/synteny_operon.tsv + résumé stdout.
"""
import json, os, re, subprocess
from Bio import SeqIO
from Bio.Seq import Seq

ROOT = "/home/christophe/docs/codes/mtbc"
CDS_H37RV = f"{ROOT}/Canettii/NC_000962.3_CDS.fasta"
HORS = f"{ROOT}/bdd/hors_mtbc"
HORS_MYCO = f"{ROOT}/bdd/hors_mycobacterium"
CANETTII_REF = f"{ROOT}/Canettii/References"
OUTDIR = f"{ROOT}/Rv1025/résultats"
os.makedirs(OUTDIR, exist_ok=True)

# Opéron H37Rv : locus -> (gene, ordre attendu)
OPERON = ["Rv1023", "Rv1024", "Rv1025", "Rv1026"]
GENE = {"Rv1023": "eno", "Rv1024": "divIC", "Rv1025": "Rv1025(DUF501)", "Rv1026": "ppx2"}
# coords H37Rv (start,end,strand) extraites du CDS fasta pour référence
H37RV_COORDS = {
    "Rv1023": (1144564, 1145853), "Rv1024": (1145858, 1146544),
    "Rv1025": (1146561, 1147028), "Rv1026": (1147019, 1147978),
}

def extract_query_proteins():
    """Extrait et traduit les 4 CDS H37Rv -> fichier .faa des protéines requêtes."""
    want = set(OPERON)
    prots = {}
    for rec in SeqIO.parse(CDS_H37RV, "fasta"):
        m = re.search(r"\[locus_tag=(Rv\d+[A-Za-z]?)\]", rec.description)
        if m and m.group(1) in want:
            aa = str(Seq(str(rec.seq)).translate(table=11, cds=False)).rstrip("*")  # type: ignore[arg-type]
            prots[m.group(1)] = aa
    faa = f"{OUTDIR}/operon_proteins.faa"
    with open(faa, "w") as fh:
        for lt in OPERON:
            fh.write(f">{lt}_{GENE[lt].split('(')[0]}\n{prots[lt]}\n")
    return faa, {lt: len(prots[lt]) for lt in OPERON}

def collect_targets():
    """Liste (label, groupe, chemin fna) des génomes cibles."""
    targets = []
    # Outgroup + MTBC contrôle
    targets.append(("M_canettii(CP007299)", "outgroup", f"{CANETTII_REF}/CP007299.1.fasta"))
    targets.append(("M_bovis(NC_002945)", "MTBC", f"{CANETTII_REF}/NC_002945.3.fasta"))
    # Toutes les espèces NTM avec un génome de référence
    for sp in sorted(os.listdir(HORS)):
        fna = f"{HORS}/{sp}/ref/genome.fna"
        if os.path.isfile(fna):
            targets.append((sp, "NTM", fna))
    # Actinobactéries hors genre Mycobacterium + témoins (clade dans meta.json)
    if os.path.isdir(HORS_MYCO):
        for sp in sorted(os.listdir(HORS_MYCO)):
            fna = f"{HORS_MYCO}/{sp}/ref/genome.fna"
            if os.path.isfile(fna):
                clade = "hors-Myco"
                try:
                    clade = json.load(open(f"{HORS_MYCO}/{sp}/ref/meta.json")).get("clade", clade)
                except Exception:
                    pass
                targets.append((sp, clade, fna))
    return targets

def run_tblastn(query_faa, subject_fna):
    """tblastn query vs subject ; retourne dict locus -> best hit."""
    cols = "qseqid sseqid pident length evalue bitscore qcovs sstart send"
    cmd = ["tblastn", "-query", query_faa, "-subject", subject_fna,
           "-evalue", "1e-5", "-max_target_seqs", "20", "-seg", "yes",
           "-outfmt", f"6 {cols}"]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=600).stdout
    except subprocess.TimeoutExpired:
        return {}
    best = {}
    for line in out.strip().splitlines():
        f = line.split("\t")
        if len(f) < 9:
            continue
        q = f[0].split("_")[0]  # Rv1023_eno -> Rv1023
        pident, evalue, bit, qcov = float(f[2]), float(f[4]), float(f[5]), float(f[6])
        sstart, send = int(f[7]), int(f[8])
        strand = "+" if send >= sstart else "-"
        lo, hi = min(sstart, send), max(sstart, send)
        hit = dict(contig=f[1], pident=pident, qcov=qcov, evalue=evalue, bit=bit,
                   lo=lo, hi=hi, mid=(lo + hi) / 2, strand=strand)
        if q not in best or bit > best[q]["bit"]:
            best[q] = hit
    return best

def analyse(label, group, best):
    """Analyse la synténie autour de l'ancre Rv1025."""
    row = {"genome": label, "group": group}
    for lt in OPERON:
        h = best.get(lt)
        row[f"{lt}_id"] = f"{h['pident']:.0f}%/{h['qcov']:.0f}%cov" if h else "ABSENT"
    anchor = best.get("Rv1025")
    if not anchor:
        row["Rv1025_present"] = "NON"
        row["divIC_adjacent"] = "NA"
        row["operon_intact"] = "NA"
        row["order_on_contig"] = "Rv1025 absent"
        return row
    row["Rv1025_present"] = f"OUI ({anchor['pident']:.0f}% id, {anchor['qcov']:.0f}% cov)"
    ac = anchor["contig"]
    # hits sur le même contig que l'ancre
    same = {lt: best[lt] for lt in OPERON if lt in best and best[lt]["contig"] == ac}
    ordered = sorted(same, key=lambda lt: same[lt]["mid"])
    row["order_on_contig"] = " ".join(f"{GENE[lt].split('(')[0]}@{int(same[lt]['mid'])}" for lt in ordered)
    # adjacence divIC-Rv1025 : même contig ET aucun autre gène de l'opéron entre eux
    div = best.get("Rv1024")
    if div and div["contig"] == ac:
        gap = same["Rv1025"]["lo"] - div["hi"] if div["mid"] < anchor["mid"] else div["lo"] - same["Rv1025"]["hi"]
        # y a-t-il un autre hit de l'opéron strictement entre divIC et Rv1025 ?
        between = [lt for lt in ("Rv1023", "Rv1026") if lt in same and
                   min(div["mid"], anchor["mid"]) < same[lt]["mid"] < max(div["mid"], anchor["mid"])]
        adj = (abs(gap) < 2000) and not between
        row["divIC_adjacent"] = f"OUI (gap {gap:+d} nt)" if adj else f"non (gap {gap:+d} nt, entre: {between})"
    else:
        row["divIC_adjacent"] = "divIC ABSENT" if not div else "autre contig"
    # opéron intact : 4 gènes, même contig, ordre eno<divIC<Rv1025<ppx2 (ou inverse)
    if len(same) == 4:
        pos = [ordered.index(lt) for lt in OPERON]
        intact = (pos == sorted(pos)) or (pos == sorted(pos, reverse=True))
        row["operon_intact"] = "OUI (4/4, ordre conservé)" if intact else f"4/4 réarrangé ({row['order_on_contig']})"
    else:
        row["operon_intact"] = f"{len(same)}/4 sur contig ancre"
    return row

def main():
    faa, qlen = extract_query_proteins()
    print(f"Protéines requêtes (aa) : " + ", ".join(f"{GENE[lt].split('(')[0]}={qlen[lt]}" for lt in OPERON))
    print(f"Opéron H37Rv (réf) : eno {H37RV_COORDS['Rv1023']} < divIC {H37RV_COORDS['Rv1024']} < "
          f"Rv1025 {H37RV_COORDS['Rv1025']} < ppx2 {H37RV_COORDS['Rv1026']}, tous brin +\n")
    targets = collect_targets()
    rows = []
    for label, group, fna in targets:
        best = run_tblastn(faa, fna)
        rows.append(analyse(label, group, best))
        r = rows[-1]
        print(f"[{group:8s}] {label:28s} Rv1025:{r['Rv1025_present']:24s} "
              f"divIC-adj:{r['divIC_adjacent']:28s} opéron:{r['operon_intact']}")
    # TSV
    cols = ["genome", "group", "Rv1025_present", "divIC_adjacent", "operon_intact",
            "order_on_contig", "Rv1023_id", "Rv1024_id", "Rv1025_id", "Rv1026_id"]
    tsv = f"{OUTDIR}/synteny_operon.tsv"
    with open(tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    # Résumé
    n = len(rows)
    present = sum(1 for r in rows if r["Rv1025_present"].startswith("OUI"))
    adj = sum(1 for r in rows if r["divIC_adjacent"].startswith("OUI"))
    intact = sum(1 for r in rows if r["operon_intact"].startswith("OUI"))
    print(f"\n=== RÉSUMÉ ({n} génomes) ===")
    print(f"Rv1025/DUF501 détecté        : {present}/{n}")
    print(f"divIC immédiatement adjacent : {adj}/{n}")
    print(f"opéron 4/4 intact & ordonné  : {intact}/{n}")
    print(f"TSV : {tsv}")

if __name__ == "__main__":
    main()
