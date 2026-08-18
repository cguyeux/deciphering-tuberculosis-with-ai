#!/usr/bin/env python3
"""
P1.4 — L'IGR purM-Rv0810c exprime-t-il un sRNA ?

Zeng et al. 2018 (BMC Genomics, PMID 29769016) affirment en une phrase du texte
principal : « Three (purM-Rv0810c, Rv3848-espR, PPE36-prcA) of the 14 drug IGRs
express small RNAs (sRNAs). » Aucune source primaire n'est appelee a cet endroit ;
la seule reference sRNA du papier renvoie a un sRNA de *Vibrio cholerae*.

Le script instruit l'affirmation sur pieces :

  (1) geometrie du locus depuis l'annotation RefSeq (GFF3 GCF_000195955.2) ;
  (2) ncRNA annotes RefSeq : presence dans l'IGR, distance au plus proche ;
  (3) balayage de SIX catalogues sRNA primaires par COORDONNEES et par NOM ;
  (4) CONTROLE POSITIF sur le meme pipeline : les 20 ncRNA RefSeq et les 10 IGR
      que Zeng lui-meme designe en Table S11 comme « expressing sRNA » ;
  (5) parsing STRUCTURE des trois catalogues qui couvrent la region :
      Miotto 2012 (candidats sRNA avec leurs colonnes), Ami 2020 (ncRv nommes,
      nombre de conditions), Arnvig 2011 (couverture RNA-seq par IGR et par CDS) ;
  (6) TEST DU DENOMINATEUR — un IGR entre deux genes CONVERGENTS recoit
      mecaniquement la lecture-through des deux transcrits. La couverture d'un tel
      IGR n'est donc pas une preuve de sRNA tant qu'elle n'est pas comparee aux
      IGR de MEME ORIENTATION. Le rang est calcule dans cette strate.

Sortie : résultats/p1_4_srna_igr.json
"""

import json
import re
import subprocess
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from statistics import median

PROJ = Path(__file__).resolve().parent.parent
EXP = PROJ / "experiments" / "2026-08-10_P1_4_P1_5"
DATA = EXP / "data"
CAT = DATA / "srna_catalogues"
OUT = PROJ / "résultats" / "p1_4_srna_igr.json"
GFF = DATA / "GCF_000195955.2_genomic.gff"

WINDOW = (904_300, 905_500)   # fenetre de balayage, genereuse a dessein
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


# --------------------------------------------------------------------------
# lecture de fichiers
# --------------------------------------------------------------------------
def xlsx_rows(path):
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        r = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in r.findall(NS + "si"):
            shared.append("".join(t.text or "" for t in si.iter(NS + "t")))
    rows = []
    for name in sorted(z.namelist()):
        if not re.match(r"xl/worksheets/sheet\d+\.xml$", name):
            continue
        root = ET.fromstring(z.read(name))
        for row in root.iter(NS + "row"):
            vals = []
            for c in row.findall(NS + "c"):
                t, v, isel = c.get("t"), c.find(NS + "v"), c.find(NS + "is")
                if t == "s" and v is not None:
                    val = shared[int(v.text or "0")]
                elif isel is not None:
                    val = "".join(x.text or "" for x in isel.iter(NS + "t"))
                elif v is not None:
                    val = v.text or ""
                else:
                    val = ""
                vals.append(val)
            rows.append(vals)
    return rows


def read_pdf_lines(path):
    try:
        r = subprocess.run(["pdftotext", "-layout", str(path), "-"],
                           capture_output=True, text=True, timeout=600)
        return r.stdout.splitlines()
    except Exception as exc:
        print(f"  [warn] pdftotext {path.name}: {exc}", file=sys.stderr)
        return []


def read_lines(path):
    suf = path.suffix.lower()
    if suf == ".xlsx":
        return ["\t".join(r) for r in xlsx_rows(path) if any(r)]
    if suf == ".pdf":
        return read_pdf_lines(path)
    if suf in (".txt", ".csv", ".tsv"):
        return path.read_text(errors="replace").splitlines()
    return []


def to_float(x, default=None):
    try:
        return float(str(x).strip())
    except (TypeError, ValueError):
        return default


# --------------------------------------------------------------------------
# geometrie du locus
# --------------------------------------------------------------------------
def load_gff():
    genes, ncrna, ordered = {}, [], []
    for line in GFF.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        typ, start, end, strand, attr = f[2], int(f[3]), int(f[4]), f[6], f[8]
        d = dict(kv.split("=", 1) for kv in attr.split(";") if "=" in kv)
        if typ in ("gene", "pseudogene"):
            tag = d.get("locus_tag")
            if not tag:
                continue
            name = d.get("gene", tag)
            g = dict(start=start, end=end, strand=strand, name=name, locus_tag=tag)
            genes[tag] = g
            if name != tag:
                genes.setdefault(name, g)
            ordered.append(g)
        elif typ == "ncRNA":
            ncrna.append(dict(start=start, end=end, strand=strand,
                              name=d.get("gene", d.get("ID", "?")),
                              product=d.get("product", "")))
    ordered.sort(key=lambda g: g["start"])
    return genes, ncrna, ordered


def orientation(left, right):
    if left["strand"] == "+" and right["strand"] == "-":
        return "convergent"
    if left["strand"] == "-" and right["strand"] == "+":
        return "divergent"
    return "tandem"


def igr_between(genes, a, b):
    if a not in genes or b not in genes:
        return None
    ga, gb = genes[a], genes[b]
    left, right = (ga, gb) if ga["start"] <= gb["start"] else (gb, ga)
    return dict(pair=f"{a}-{b}",
                start=left["end"] + 1, end=right["start"] - 1,
                length=right["start"] - left["end"] - 1,
                left=f'{left["locus_tag"]}({left["strand"]})',
                right=f'{right["locus_tag"]}({right["strand"]})',
                orientation=orientation(left, right))


def all_igrs(ordered):
    """Tous les IGR non chevauchants du genome, avec leur orientation."""
    out = []
    for a, b in zip(ordered, ordered[1:]):
        if b["start"] - a["end"] - 1 <= 0:
            continue
        out.append(dict(start=a["end"] + 1, end=b["start"] - 1,
                        length=b["start"] - a["end"] - 1,
                        orientation=orientation(a, b),
                        pair=f'{a["locus_tag"]}-{b["locus_tag"]}'))
    return out


# --------------------------------------------------------------------------
# balayage generique
# --------------------------------------------------------------------------
COORD_RE = re.compile(r"(?<![\d.,])(\d{6,7})(?![\d])")


def scan(lines, lo, hi, patterns):
    by_coord, by_name = [], []
    for ln in lines:
        if any(lo <= int(m.group(1)) <= hi for m in COORD_RE.finditer(ln)):
            by_coord.append(ln.strip()[:300])
        if any(p.search(ln) for p in patterns):
            by_name.append(ln.strip()[:300])
    return by_coord, by_name


CATALOGUES = {
    "arnvig2009_molmicrobiol_PMID19555452": ["19555452/mmi0073-0397-SD1.pdf"],
    "dichiara2010_nar_PMID20181675": ["converted/gkq101_nar-02696-a-2009-File007.txt"],
    "arnvig2011_plospathog_PMID22072964": [
        "converted/ppat.1002342.s008.xlsx", "converted/ppat.1002342.s009.xlsx",
        "converted/ppat.1002342.s011.xlsx", "converted/ppat.1002342.s012.xlsx",
        "converted/ppat.1002342.s007.txt", "converted/ppat.1002342.s010.txt"],
    "miotto2012_plosone_PMID23284830": [
        "23284830/pone.0051950.s007.xlsx", "23284830/pone.0051950.s008.xlsx",
        "converted/pone.0051950.s009.xlsx", "23284830/pone.0051950.s010.xlsx",
        "23284830/pone.0051950.s011.xlsx", "23284830/pone.0051950.s012.xlsx"],
    "gerrick2018_pnas_PMID29871950": [
        "29871950/pnas.1718003115.sd01.xlsx", "29871950/pnas.1718003115.sd02.xlsx",
        "29871950/pnas.1718003115.sd03.xlsx", "29871950/pnas.1718003115.sd04.xlsx",
        "29871950/pnas.1718003115.sapp.pdf"],
    "ami2020_bmcgenomics_PMID32070281": [
        "32070281/12864_2020_6573_MOESM8_ESM.pdf",
        "32070281/12864_2020_6573_MOESM9_ESM.pdf",
        "32070281/12864_2020_6573_MOESM10_ESM.pdf",
        "32070281/12864_2020_6573_MOESM11_ESM.pdf"],
}

ZENG_S11 = [("ctpE", "Rv0909"), ("lipX", "mshB"), ("nrdH", "Rv3054c"),
            ("PE12", "fbiC"), ("rplM", "esxT"), ("Rv0922", "Rv0923c"),
            ("Rv1045", "Rv1047"), ("Rv1179c", "pks3"), ("Rv3402c", "Rv3403c"),
            ("serB1", "mmpS2")]


# --------------------------------------------------------------------------
# parsing structure des trois catalogues qui couvrent la region
# --------------------------------------------------------------------------
def parse_miotto(lo, hi):
    """Table S1 de Miotto 2012 : candidats sRNA avec colonnes nommees."""
    rows = xlsx_rows(CAT / "23284830/pone.0051950.s007.xlsx")
    hdr_i = next(i for i, r in enumerate(rows)
                 if any("start" == (c or "").strip().lower() for c in r))
    hdr = [c.strip() for c in rows[hdr_i]]
    hits, total = [], 0
    for r in rows[hdr_i + 1:]:
        if len(r) < 4 or not r[0].startswith("candidate_"):
            continue
        total += 1
        s, e = to_float(r[1]), to_float(r[2])
        if s is None or e is None:
            continue
        if e >= lo and s <= hi:
            hits.append(dict(zip(hdr, r)))
    return {"n_candidats_total": total, "hits_fenetre": hits}


def parse_ami(lo, hi):
    """Tables S4 d'Ami 2020 : sRNA nommes ncRv..., coordonnees + conditions."""
    lines = read_pdf_lines(CAT / "32070281/12864_2020_6573_MOESM10_ESM.pdf")
    coords, conds = {}, {}
    INT = re.compile(r"\d{4,7}$")
    for ln in lines:
        f = ln.split()
        if not f or not f[0].startswith("ncRv"):
            continue
        name = f[0]
        # Table S4a/S4b : nom, start, end, localisation TEXTUELLE [, n conditions].
        # Le 4e champ doit etre non numerique : sans cette garde, la table S4d
        # (15 valeurs d'expression) passe pour des coordonnees et les ecrase.
        if (len(f) >= 4 and INT.fullmatch(f[1]) and INT.fullmatch(f[2])
                and to_float(f[3]) is None):
            s, e = int(f[1]), int(f[2])
            loc = " ".join(f[3:]).strip()
            m = re.search(r"\s(\d+)$", loc)
            n_cond = int(m.group(1)) if m else None
            if n_cond is not None:
                loc = loc[:m.start()].strip()
            prev = coords.get(name, {})
            coords[name] = dict(start=s, end=e, location=loc or prev.get("location", ""),
                                n_conditions=n_cond if n_cond is not None
                                else prev.get("n_conditions"))
        # Table S4c : matrice binaire (nom puis 15 x 0/1 puis total)
        elif len(f) == 17 and all(x in ("0", "1") for x in f[1:16]):
            conds[name] = dict(binaire=f[1:16], total=int(f[16]))
    for k, v in conds.items():
        coords.setdefault(k, {}).update(v)
    hits = {k: v for k, v in coords.items()
            if v.get("start") and v["end"] >= lo and v["start"] <= hi}
    return {"n_srna_nommes": len(coords), "hits_fenetre": hits}


def parse_arnvig(igr_purm, igr_811, genome_igrs):
    """Tables S2 (CDS) et S3 (IGR) d'Arnvig 2011 + test du denominateur."""
    ig_rows = xlsx_rows(CAT / "converted/ppat.1002342.s009.xlsx")
    hdr = [c.strip() for c in ig_rows[0]]
    igs = []
    for r in ig_rows[1:]:
        if len(r) < 6 or not r[0].startswith("IG"):
            continue
        s, e, ln = to_float(r[1]), to_float(r[2]), to_float(r[3])
        exp, sta = to_float(r[4]), to_float(r[5])
        if None in (s, e, ln, exp) or ln <= 0:
            continue
        igs.append(dict(region=r[0], start=int(s), end=int(e), length=int(ln),
                        exp_reads=exp, sta_reads=sta,
                        exp_per_nt=exp / ln,
                        sta_per_nt=(sta / ln) if sta is not None else None))
    # orientation de chaque IG d'Arnvig, par appariement de coordonnees
    by_span = {(g["start"], g["end"]): g["orientation"] for g in genome_igrs}
    for ig in igs:
        o = by_span.get((ig["start"], ig["end"]))
        if o is None:  # tolerance de +-2 nt (conventions de bornes)
            for (s, e), oo in by_span.items():
                if abs(s - ig["start"]) <= 2 and abs(e - ig["end"]) <= 2:
                    o = oo
                    break
        ig["orientation"] = o or "inconnue"

    def rank_of(target, strata=None):
        pool = [g for g in igs if strata is None or g["orientation"] == strata]
        pool = [g for g in pool if g["exp_per_nt"] is not None]
        vals = sorted((g["exp_per_nt"] for g in pool), reverse=True)
        r = vals.index(target["exp_per_nt"]) + 1
        return dict(rang=r, n=len(vals),
                    percentile=round(100 * (1 - (r - 1) / len(vals)), 1),
                    mediane_strate=round(median(vals), 3))

    def find(start, end):
        for g in igs:
            if abs(g["start"] - start) <= 3 and abs(g["end"] - end) <= 3:
                return g
        return None

    ig624 = find(igr_purm["start"], igr_purm["end"])
    ig625 = find(igr_811["start"], igr_811["end"])
    res = {"n_IGR_table": len(igs),
           "orientations": {o: sum(1 for g in igs if g["orientation"] == o)
                            for o in ("convergent", "divergent", "tandem", "inconnue")}}
    for label, ig in (("IGR_purM_Rv0810c", ig624), ("IGR_Rv0810c_Rv0811c", ig625)):
        if ig is None:
            res[label] = None
            continue
        res[label] = {**{k: v for k, v in ig.items()},
                      "rang_global": rank_of(ig),
                      "rang_dans_sa_strate": rank_of(ig, ig["orientation"])}
    # couverture par strate d'orientation : l'attendu geometrique
    res["mediane_exp_par_nt_par_orientation"] = {
        o: round(median([g["exp_per_nt"] for g in igs if g["orientation"] == o]), 3)
        for o in ("convergent", "divergent", "tandem")
        if any(g["orientation"] == o for g in igs)}

    # table S2 : CDS, transcription sens et antisens
    cds_rows = xlsx_rows(CAT / "converted/ppat.1002342.s008.xlsx")
    chdr = [c.strip() for c in cds_rows[0]]
    cds = {}
    for r in cds_rows[1:]:
        if len(r) < len(chdr) or not r[0].strip().startswith("Rv"):
            continue
        cds[r[0].strip()] = dict(zip(chdr, [c.strip() for c in r]))
    out_cds = {}
    for tag in ("Rv0810c", "Rv0809", "Rv0811c"):
        if tag not in cds:
            continue
        row = cds[tag]
        as_rpkm = to_float(row.get("exp_AS_RPKM"))
        sense = to_float(row.get("exp_sense_RPKM"))
        pool = sorted((to_float(v.get("exp_AS_RPKM"), 0) or 0 for v in cds.values()),
                      reverse=True)
        rank = pool.index(as_rpkm) + 1 if as_rpkm in pool else None
        out_cds[tag] = {
            **row,
            "ratio_AS_sur_sens_exp": round(as_rpkm / sense, 3)
            if (as_rpkm and sense) else None,
            "rang_AS_RPKM": rank, "n_CDS": len(pool),
            "percentile_AS_RPKM": round(100 * (1 - (rank - 1) / len(pool)), 1)
            if rank else None,
        }
    res["CDS"] = out_cds
    return res


# --------------------------------------------------------------------------
def main():
    genes, refseq_ncrna, ordered = load_gff()
    genome_igrs = all_igrs(ordered)
    res = {"question": "L'IGR purM-Rv0810c exprime-t-il un sRNA ? (Zeng et al. 2018)"}

    igr_purm = igr_between(genes, "purM", "Rv0810c")
    igr_811 = igr_between(genes, "Rv0810c", "Rv0811c")
    assert igr_purm and igr_811
    res["geometrie"] = {
        "Rv0810c": {k: genes["Rv0810c"][k] for k in ("start", "end", "strand")},
        "purM_Rv0809": {k: genes["purM"][k] for k in ("start", "end", "strand")},
        "Rv0811c": {k: genes["Rv0811c"][k] for k in ("start", "end", "strand")},
        "IGR_purM_Rv0810c": igr_purm,
        "IGR_Rv0810c_Rv0811c": igr_811,
    }
    print(f"IGR purM-Rv0810c   : {igr_purm['start']}-{igr_purm['end']} "
          f"({igr_purm['length']} nt, {igr_purm['orientation']})")
    print(f"IGR Rv0810c-Rv0811c: {igr_811['start']}-{igr_811['end']} "
          f"({igr_811['length']} nt, {igr_811['orientation']})")

    # (2) ncRNA RefSeq
    inside = [n for n in refseq_ncrna
              if n["end"] >= igr_purm["start"] and n["start"] <= igr_purm["end"]]
    nearest = min(refseq_ncrna, key=lambda n: min(abs(n["start"] - igr_purm["end"]),
                                                  abs(igr_purm["start"] - n["end"])))
    dist = min(abs(nearest["start"] - igr_purm["end"]),
               abs(igr_purm["start"] - nearest["end"]))
    res["refseq_ncrna"] = {"n_total_genome": len(refseq_ncrna),
                           "dans_IGR_purM_Rv0810c": inside,
                           "plus_proche": {**nearest, "distance_nt": dist},
                           "liste": sorted(refseq_ncrna, key=lambda n: n["start"])}
    print(f"ncRNA RefSeq : {len(refseq_ncrna)} au genome, {len(inside)} dans l'IGR ; "
          f"le plus proche ({nearest['name']}) est a {dist} nt")

    # (3) balayage des catalogues
    pats = [re.compile(p, re.I) for p in
            (r"\bRv0810", r"\bRv0809\b", r"\bpurM\b", r"ncRv0809", r"ncRv0810")]
    lo, hi = WINDOW
    cat_res, all_lines = {}, {}
    for cat, files in CATALOGUES.items():
        bc, bn, missing, nl = [], [], [], 0
        for rel in files:
            p = CAT / rel
            if not p.exists():
                missing.append(rel)
                continue
            lines = read_lines(p)
            all_lines[rel] = lines
            nl += len(lines)
            c, n = scan(lines, lo, hi, pats)
            bc += [f"{p.name}: {x}" for x in c]
            bn += [f"{p.name}: {x}" for x in n]
        cat_res[cat] = {"n_fichiers": len(files) - len(missing),
                        "fichiers_absents": missing, "n_lignes_scannees": nl,
                        "hits_par_coordonnee": bc, "hits_par_nom": bn}
        print(f"{cat}: {nl} lignes, {len(bc)} coord, {len(bn)} nom")
    res["catalogues"] = cat_res

    # (4) controle positif
    ctrl_igr = {}
    for a, b in ZENG_S11:
        igr = igr_between(genes, a, b)
        if igr is None:
            ctrl_igr[f"{a}-{b}"] = {"erreur": "gene absent du GFF"}
            continue
        clo, chi = igr["start"] - 250, igr["end"] + 250
        ipats = [re.compile(rf"\b{re.escape(a)}\b", re.I),
                 re.compile(rf"\b{re.escape(b)}\b", re.I)]
        tc = tn = 0
        cats_hit = []
        for cat, files in CATALOGUES.items():
            hc = hn = 0
            for rel in files:
                if rel in all_lines:
                    c, n = scan(all_lines[rel], clo, chi, ipats)
                    hc += len(c)
                    hn += len(n)
            tc += hc
            tn += hn
            if hc or hn:
                cats_hit.append(cat.split("_")[0])
        ctrl_igr[igr["pair"]] = {
            "intervalle": [igr["start"], igr["end"]], "longueur_nt": igr["length"],
            "orientation": igr["orientation"],
            "hits_coord": tc, "hits_nom": tn, "catalogues_touches": cats_hit,
            "ncRNA_RefSeq_dans_IGR": [n["name"] for n in refseq_ncrna
                                      if n["end"] >= igr["start"]
                                      and n["start"] <= igr["end"]]}
    ctrl_nc = {}
    for n in refseq_ncrna:
        clo, chi = n["start"] - 200, n["end"] + 200
        npats = [re.compile(rf"\b{re.escape(n['name'])}\b")]
        hits = []
        for cat, files in CATALOGUES.items():
            hc = hn = 0
            for rel in files:
                if rel in all_lines:
                    c, nn = scan(all_lines[rel], clo, chi, npats)
                    hc += len(c)
                    hn += len(nn)
            if hc or hn:
                hits.append(f"{cat.split('_')[0]}({hc}c/{hn}n)")
        ctrl_nc[n["name"]] = {"intervalle": [n["start"], n["end"]],
                              "catalogues_touches": hits, "n_catalogues": len(hits)}
    n_ok = sum(1 for v in ctrl_nc.values() if v["n_catalogues"])
    res["controle_positif"] = {
        "igr_zeng_table_s11": ctrl_igr, "ncrna_refseq": ctrl_nc,
        "resume": {
            "ncRNA_RefSeq_retrouves": f"{n_ok}/{len(ctrl_nc)}",
            "IGR_S11_avec_hit":
                f"{sum(1 for v in ctrl_igr.values() if v.get('hits_coord') or v.get('hits_nom'))}"
                f"/{len(ctrl_igr)}"}}
    print(f"CONTROLE : {n_ok}/{len(ctrl_nc)} ncRNA RefSeq retrouves")

    # (5) parsing structure
    res["miotto2012_structure"] = parse_miotto(lo, hi)
    res["ami2020_structure"] = parse_ami(lo, hi)
    res["arnvig2011_structure"] = parse_arnvig(igr_purm, igr_811, genome_igrs)

    print(f"Miotto 2012 : {res['miotto2012_structure']['n_candidats_total']} candidats, "
          f"{len(res['miotto2012_structure']['hits_fenetre'])} dans la fenetre")
    print(f"Ami 2020    : {res['ami2020_structure']['n_srna_nommes']} sRNA nommes, "
          f"{len(res['ami2020_structure']['hits_fenetre'])} dans la fenetre "
          f"-> {list(res['ami2020_structure']['hits_fenetre'])}")
    a = res["arnvig2011_structure"]
    for lab in ("IGR_purM_Rv0810c", "IGR_Rv0810c_Rv0811c"):
        if a.get(lab):
            print(f"Arnvig {lab}: {a[lab]['exp_reads']} reads / {a[lab]['length']} nt, "
                  f"rang global {a[lab]['rang_global']['rang']}/{a[lab]['rang_global']['n']}, "
                  f"rang dans la strate {a[lab]['orientation']} "
                  f"{a[lab]['rang_dans_sa_strate']['rang']}/"
                  f"{a[lab]['rang_dans_sa_strate']['n']}")
    print(f"Mediane couverture/nt par orientation : "
          f"{a['mediane_exp_par_nt_par_orientation']}")

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(res, indent=2, ensure_ascii=False))
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
