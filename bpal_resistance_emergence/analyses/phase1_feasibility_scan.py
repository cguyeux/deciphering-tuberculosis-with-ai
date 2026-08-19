#!/usr/bin/env python3
"""Phase 1 — Sonde de faisabilité BPaLM.

Balaye le pangénome SPDI (souche -> {SPDI}, projet voisin Resistance_antibio) en
le restreignant aux variants tombant dans les gènes du panel BPaLM
(data/gene_panel.tsv). Pour chaque variant : effet codon (annoté de façon
déterministe depuis H37Rv), nombre de porteurs, ventilation par lignée Coll
(level_1), et statut dans le catalogue WHO consolidé.

But : go/no-go sur le substrat avant d'engager le pipeline complet. On vérifie
notamment (1) que les déterminants BDQ catalogués (Rv0678) ressortent avec des
comptes plausibles, (2) que gyrA D94 est présent (contrôle positif FQ),
(3) que la voie F420 (PMD) contient des variants NON catalogués = la place pour
la découverte.

Sortie : résultats/phase1_feasibility.tsv + récap par gène à l'écran.

Usage :
    python phase1_feasibility_scan.py [--limit N] [--min-carriers K]
    --limit N      : ne scanner que N souches (validation rapide de la logique)
    --min-carriers : seuil d'affichage dans le TSV (défaut 1)
"""
import argparse
import pickle
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths

REF = paths.REF_ACC
PREFIX = REF + ":"
PLEN = len(PREFIX)

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def load_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


def load_panel():
    """Renvoie : pos2gene {pos_0based: locus}, genes {locus: dict}."""
    pos2gene = {}
    genes = {}
    with open(paths.GENE_PANEL) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            locus = r["locus"]
            s, e = int(r["start"]), int(r["end"])
            genes[locus] = {
                "gene": r["gene"], "start": s, "end": e, "strand": r["strand"],
                "biotype": r["biotype"], "drug": r["drug"], "group": r["panel_group"],
            }
            # SPDI = 0-based ; gène 1-based [s,e] -> positions SPDI [s-1, e-1]
            for p in range(s - 1, e):
                pos2gene[p] = locus
    return pos2gene, genes


def annotate_effect(genes, locus, spdi_pos, ref, alt, fasta):
    """Effet d'un SNV (SPDI 0-based). Renvoie (vtype, effect, aa_change, gene_nt_pos)."""
    g = genes[locus]
    G = spdi_pos + 1                       # position génomique 1-based
    nt_pos = (G - g["start"] + 1) if g["strand"] == "+" else (g["end"] - G + 1)
    if len(ref) != 1 or len(alt) != 1:
        return ("indel/MNV", "complex", "", nt_pos)
    if g["biotype"] == "rRNA":
        return ("SNV", "rRNA", f"r.{nt_pos}{ref}>{alt}", nt_pos)
    # garde-fou : la base ref du SPDI doit coïncider avec H37Rv
    if fasta[spdi_pos].upper() != ref.upper():
        return ("SNV", "ref_mismatch", "", nt_pos)
    s, e, strand = g["start"], g["end"], g["strand"]
    if strand == "+":
        offset = G - s
        cstart = s + (offset // 3) * 3          # 1-based début codon
        plus = fasta[cstart - 1:cstart + 2]
        i = offset % 3
        ref_codon = plus
        alt_codon = plus[:i] + alt + plus[i + 1:]
    else:
        offset = e - G
        cend = e - (offset // 3) * 3             # 1-based base 5' (coord la + haute)
        plus = fasta[cend - 3:cend]
        i = G - (cend - 2)                       # index 0-based dans 'plus' (+ strand)
        alt_plus = plus[:i] + alt + plus[i + 1:]
        ref_codon = revcomp(plus)
        alt_codon = revcomp(alt_plus)
    aa_n = offset // 3 + 1
    aa_ref = CODON_TABLE.get(ref_codon.upper())
    aa_alt = CODON_TABLE.get(alt_codon.upper())
    if aa_ref is None or aa_alt is None:
        return ("SNV", "complex", "", nt_pos)
    if aa_ref == aa_alt:
        eff = "synonymous"
    elif aa_alt == "*":
        eff = "stop_gained"
    elif aa_ref == "*":
        eff = "stop_lost"
    elif aa_n == 1:
        eff = "start_related"
    else:
        eff = "missense"
    return ("SNV", eff, f"p.{aa_ref}{aa_n}{aa_alt}", nt_pos)


def load_catalogue():
    """spdi -> 'drug:call' agrégé (drogues pertinentes BPaLM + proxies)."""
    import csv
    keep = {"bedaquiline", "pretomanid", "delamanid", "clofazimine",
            "linezolid", "moxifloxacin", "levofloxacin", "rifampicin", "rifabutin"}
    cat = defaultdict(list)
    with open(paths.CATALOGUE_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if row["drug"] in keep and row["spdi"]:
                cat[row["spdi"]].append((row["drug"], row["call"]))
    return cat


def load_lineages():
    import csv
    m = {}
    with open(paths.LINEAGE_SNAPSHOT) as fh:
        r = csv.DictReader(fh)
        for row in r:
            m[row["strain_name"]] = row["lineage_level_1"]
    return m


def lin_label(x):
    if x in (None, "", "NA"):
        return "NA"
    if x in ("BOV", "BOV_AFRI"):
        return x
    return "L" + x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--min-carriers", type=int, default=1)
    args = ap.parse_args()

    paths.ensure_dirs()
    print("panel :", end=" ")
    pos2gene, genes = load_panel()
    print(f"{len(genes)} gènes, {len(pos2gene)} positions couvertes")

    print("chargement H37Rv…")
    fasta = load_fasta(paths.H37RV_FASTA)
    print(f"  génome {len(fasta)} pb")

    print("chargement catalogue consolidé…")
    cat = load_catalogue()
    print(f"  {len(cat)} SPDI catalogués (drogues pertinentes)")

    print("chargement snapshot lignées Coll…")
    strain2lin = load_lineages()
    print(f"  {len(strain2lin)} souches typées")

    print("chargement pan_strains.pkl (souche -> {SPDI})… [volumineux]")
    with open(paths.PAN_STRAINS, "rb") as f:
        pan = pickle.load(f)
    print(f"  {len(pan)} souches dans le pangénome")

    # scan
    var_total = Counter()                       # spdi -> n porteurs
    var_lin = defaultdict(Counter)              # spdi -> {lignée: n}
    N = 0
    for i, (strain, spdis) in enumerate(pan.items()):
        if args.limit and i >= args.limit:
            break
        N += 1
        lab = lin_label(strain2lin.get(strain))
        for spdi in spdis:
            # parse rapide de la position (0-based)
            j = spdi.find(":", PLEN)
            if j < 0:
                continue
            try:
                pos = int(spdi[PLEN:j])
            except ValueError:
                continue
            if pos in pos2gene:
                var_total[spdi] += 1
                var_lin[spdi][lab] += 1
        if (i + 1) % 20000 == 0:
            print(f"  …{i + 1} souches, {len(var_total)} variants panel")
    print(f"scan terminé : {N} souches, {len(var_total)} variants distincts dans le panel")

    # annotation + écriture
    rows = []
    for spdi, n in var_total.items():
        parts = spdi.split(":")
        pos = int(parts[1]); ref = parts[2]; alt = parts[3]
        locus = pos2gene[pos]
        g = genes[locus]
        vtype, eff, aa, ntpos = annotate_effect(genes, locus, pos, ref, alt, fasta)
        catst = cat.get(spdi, [])
        cat_str = ";".join(f"{d}:{c}" for d, c in catst)
        is_R = any(c == "R-associated" for _, c in catst)
        linbk = ";".join(f"{k}:{v}" for k, v in var_lin[spdi].most_common())
        rows.append({
            "group": g["group"], "drug": g["drug"], "locus": locus, "gene": g["gene"],
            "spdi": spdi, "ref": ref, "alt": alt, "vtype": vtype, "effect": eff,
            "aa_change": aa, "gene_nt_pos": ntpos, "n_carriers": n,
            "freq_pct": round(100 * n / N, 4), "n_lineages": len(var_lin[spdi]),
            "catalogue": cat_str, "is_catalogued_R": int(is_R),
            "lineage_breakdown": linbk,
        })
    cols = ["group", "drug", "locus", "gene", "spdi", "ref", "alt", "vtype", "effect",
            "aa_change", "gene_nt_pos", "n_carriers", "freq_pct", "n_lineages",
            "catalogue", "is_catalogued_R", "lineage_breakdown"]
    out = paths.RESULTATS / "phase1_feasibility.tsv"
    rows.sort(key=lambda r: (r["group"], r["locus"], -r["n_carriers"]))
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            if r["n_carriers"] >= args.min_carriers:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\nécrit -> {out}  ({sum(1 for r in rows if r['n_carriers']>=args.min_carriers)} variants)")

    # récap par gène
    print("\n=== RÉCAP PAR GÈNE ===")
    print(f"{'gène':8} {'grp':4} {'var':>5} {'SNV':>5} {'miss':>5} {'syn':>5} "
          f"{'cat':>5} {'catR':>5} {'novNS':>6} {'porteurs':>9}")
    by_g = defaultdict(list)
    for r in rows:
        by_g[r["locus"]].append(r)
    for locus, g in sorted(genes.items(), key=lambda kv: (kv[1]["group"], kv[0])):
        rs = by_g.get(locus, [])
        nvar = len(rs)
        nsnv = sum(1 for r in rs if r["vtype"] == "SNV")
        nmiss = sum(1 for r in rs if r["effect"] in ("missense", "stop_gained"))
        nsyn = sum(1 for r in rs if r["effect"] == "synonymous")
        ncat = sum(1 for r in rs if r["catalogue"])
        ncatR = sum(1 for r in rs if r["is_catalogued_R"])
        nnov = sum(1 for r in rs if r["effect"] in ("missense", "stop_gained") and not r["catalogue"])
        ncarr = sum(r["n_carriers"] for r in rs)
        print(f"{g['gene']:8} {g['group']:4} {nvar:5d} {nsnv:5d} {nmiss:5d} {nsyn:5d} "
              f"{ncat:5d} {ncatR:5d} {nnov:6d} {ncarr:9d}")

    # sanity checks
    print("\n=== SANITY CHECKS ===")
    rv0678 = [r for r in rows if r["locus"] == "Rv0678"]
    print(f"Rv0678 : {len(rv0678)} variants, "
          f"{sum(1 for r in rv0678 if r['is_catalogued_R'])} catalogués-R, "
          f"{sum(r['n_carriers'] for r in rv0678 if r['is_catalogued_R'])} porteurs catalogués-R")
    gyrA94 = [r for r in rows if r["locus"] == "Rv0006" and r["aa_change"].startswith("p.D94")]
    print(f"gyrA D94* (contrôle FQ) : {len(gyrA94)} variants -> "
          + ", ".join(f"{r['aa_change']}({r['n_carriers']})" for r in gyrA94[:6]))
    f420 = [r for r in rows if r["group"] == "PMD" and r["effect"] in ("missense", "stop_gained")]
    f420_nov = [r for r in f420 if not r["catalogue"]]
    print(f"voie F420 (PMD) : {len(f420)} variants non-syn, dont {len(f420_nov)} NON catalogués "
          f"(place pour la découverte)")


if __name__ == "__main__":
    main()
