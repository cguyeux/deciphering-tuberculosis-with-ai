#!/usr/bin/env python3
"""Phase 2c — minage des tables du CORPS d'article (JATS <table-wrap>) du corpus litcohort.

Pourquoi ce script existe (cf. cahier 2026-06-12 « double panne » et la reprise P3.5).
`phase2b_litcohort_dst.py` ne lit que les fichiers SUPPLEMENTARY. Or pour Frontiers,
MDPI, BMC et AAC, l'endpoint EPMC `/supplementaryFiles` ne renvoie que des IMAGES de
figures : les tables de DST sont dans le corps de l'article (`fullTextXML`, balises
`<table-wrap>`), déjà en cache et jamais parsées. Le zéro-rendement du batch du 12 juin
n'était donc pas seulement une panne EBI, c'était surtout un angle mort du parseur.

Constat dirimant mesuré avant d'écrire ce script : sur 40 tables du corps, 15 portent
une colonne médicament et 2 portent des accessions, mais AUCUNE ne porte les deux. Le
DST par souche du PMD/DLM n'est pas publié sous une forme reliable à une accession SRA
(identifiants internes + dépôt en BioProject global). L'objectif initial (densifier
`res.csv` en paires accession↔phénotype) est donc structurellement hors d'atteinte pour
ces molécules ; ce script vise l'objectif atteignable : une table de preuves au niveau
VARIANT (gène, mutation, médicament, MIC/R-S, dispositif), citable dans le manuscrit.

Sorties :
  résultats/phase2c_f420_phenotype_evidence.tsv   preuves (gène panel × médicament)
  résultats/phase2c_table_inventory.tsv           audit de toutes les tables + liabilité

Lancer avec le venv qui porte xlrd/python-docx (formats supp non-XLSX) :
  /home/christophe/venvs/litcohort/bin/python analyses/phase2c_litcohort_jats.py
"""
import csv
import io
import re
import subprocess
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import pandas as pd
from lxml import etree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths

CACHE = paths.RESULTATS / "litcohort_cache"
OUT_EVID = paths.RESULTATS / "phase2c_f420_phenotype_evidence.tsv"
OUT_INV = paths.RESULTATS / "phase2c_table_inventory.tsv"
OUT_XCHK = paths.RESULTATS / "phase2c_candidate_crosscheck.tsv"
MANUAL_EVID = paths.DATA / "litcohort_manual_evidence.tsv"

ACC_RE = re.compile(r"\b([SED]RR\d{5,})\b")
PROJ_RE = re.compile(r"\b(PRJ[EDN][A-Z]\d+)\b")

# --- médicaments -----------------------------------------------------------
# Deux niveaux : codes courts (frontière de mot obligatoire, sinon faux positifs)
# et noms longs (sous-chaîne suffisante). « Pa » = prétomanide dans les tables
# éthiopiennes ; « DEL » = délamanide chez Frontiers.
DRUG_LONG = {
    "bedaquiline": "bedaquiline", "pretomanid": "pretomanid", "delamanid": "delamanid",
    "clofazim": "clofazimine", "linezolid": "linezolid", "moxiflox": "moxifloxacin",
    "levofloxac": "levofloxacin", "pa-824": "pretomanid", "pa824": "pretomanid",
    "opc-67683": "delamanid",
}
DRUG_SHORT = {
    "bdq": "bedaquiline", "pmd": "pretomanid", "pa": "pretomanid", "dlm": "delamanid",
    "del": "delamanid", "del amanid": "delamanid", "cfz": "clofazimine", "clf": "clofazimine",
    "cfx": "clofazimine", "lzd": "linezolid", "mfx": "moxifloxacin", "mox": "moxifloxacin",
    "lfx": "levofloxacin",
}

# --- gènes du panel --------------------------------------------------------
PANEL = {}
with open(paths.GENE_PANEL) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        PANEL[r["gene"].lower()] = (r["gene"], r["panel_group"])
        PANEL[r["locus"].lower()] = (r["gene"], r["panel_group"])
PANEL["mmpr5"] = PANEL.get("mmpr5", ("mmpR5", "BDQ"))
GENE_RE = re.compile(r"\b(" + "|".join(sorted(map(re.escape, PANEL), key=len, reverse=True)) + r")\b",
                     re.I)

AA3 = ("Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val")
MUT_RE = re.compile(
    r"(?:p\.)?(?:" + AA3 + r")\d+(?:" + AA3 + r"|\*|fs|Ter|del|ins|dup)"      # Arg280Cys, Ala855fs
    r"|(?:p\.)?[ACDEFGHIKLMNPQRSTVWY]\d+(?:[ACDEFGHIKLMNPQRSTVWY]|\*|fs)"     # A210E, W88*, L49P
    r"|\d+_\d+(?:del|ins|dup)[A-Za-z]*"                                        # 559_561delTGC
    r"|\d+(?:del|ins|dup)[A-Za-z]+"                                            # 344delA, 214dupC
    r"|\d+_del_\d+_[A-Za-z]+_[A-Za-z]+"                                        # 398_del_1_TC_T (ACS Omega)
    r"|c\.\d+[_\d]*[A-Za-z>\s]*[ACGT]",                                        # c.1134C > T
    re.I)
WT_TOKENS = {"", "-", "nan", "na", "nd", "wt", "gwt", "h37rv", "ref", "reference",
             "none", "wild-type", "wildtype", "0", "no", "s"}

# Les sauts de ligne HTML perdus collent le gène suivant à la mutation précédente :
# « fgd1: Arg280CysfbiC: Val318IlerplC: Cys154Arg ». Sans frontière de mot, fbiC et
# rplC ne sont plus détectés et TOUTES les mutations retombent sur fgd1 (faux
# « fgd1 Cys154Arg », alors que Cys154Arg est le déterminant LZD canonique de rplC).
UNGLUE_RE = re.compile(r"(?<=[A-Za-z0-9])(?=(?:" +
                       "|".join(sorted(map(re.escape, PANEL), key=len, reverse=True)) +
                       r")\s*:)", re.I)


def unglue(s):
    return UNGLUE_RE.sub(" ", s)


def norm_hdr(h):
    return re.sub(r"[^a-z0-9]+", " ", str(h).lower()).strip()


def drug_of_header(h):
    """Médicament d'un en-tête, ou None.

    GARDE-FOU : un en-tête qui contient un token de mutation est une colonne de
    GÉNOTYPE, jamais une colonne de médicament. Sans cela, l'en-tête fusionné
    « mutation 281_del_21_CAACCCC… » fait matcher l'alias court « del » et la
    colonne passe pour du délamanide (faux positif constaté au 2e passage)."""
    if MUT_RE.search(str(h)):
        return None
    n = norm_hdr(h)
    for k, d in DRUG_LONG.items():
        if k in n:
            return d
    for k, d in DRUG_SHORT.items():
        if re.search(r"\b" + re.escape(k) + r"\b", n):
            return d
    return None


MIC_CELL = re.compile(
    r"^\s*([<>≤≥]?=?\s*\d*\.?\d+)\s*(?:µg/?m?l|ug/?m?l|mg/?l)?\s*"     # la MIC
    r"\(?\s*(r|s|siu|riu)?\s*\)?\s*$", re.I)                            # label optionnel
RS_CELL = re.compile(r"^\s*\(?\s*(r|s|siu|riu)\s*\)?\s*$", re.I)


def parse_pheno(raw):
    """(rs, mic, raw) depuis une cellule.

    GARDE-FOU (bug vécu, cf. cahier 2026-06-12 et 2e passage de ce script) : une MIC
    doit valider la cellule ENTIÈRE, jamais être un chiffre glané en tête. Sinon
    « Rv0678 M146T » donne un faux 146, et la mutation nucléotidique
    « 281_del_21_CAACCCC… » un faux 281 — les deux lus comme des MIC résistantes."""
    s = str(raw).strip()
    if s.lower() in WT_TOKENS - {"s"}:
        return "", "", s
    if MUT_RE.search(s):                       # cellule de génotype : pas un phénotype
        return "", "", s
    # marqueurs de réplicat « (x3) » / « (2x) » et valeurs multiples « a;b » -> 1re valeur
    s2 = re.sub(r"\((?:x\s*\d+|\d+\s*x)\)", "", s).split(";")[0].strip()
    m = MIC_CELL.match(s2)
    if m:
        mic = re.sub(r"\s+", "", m.group(1)).replace("≤", "<=").replace("≥", ">=")
        lab = (m.group(2) or "").lower()
        return ({"siu": "S", "riu": "R"}.get(lab, lab.upper())), mic, s
    m = RS_CELL.match(s2)
    if m:
        lab = m.group(1).lower()
        return {"siu": "S", "riu": "R"}.get(lab, lab.upper()), "", s
    return "", "", s


def merge_subheaders(df, max_rounds=3):
    """Fusionne les en-têtes sur plusieurs lignes (fréquent en JATS : la ligne 0 du
    corps porte les vrais noms de médicaments sous un en-tête « pDST (MIC) »)."""
    for _ in range(max_rounds):
        if df.empty:
            break
        r0 = df.iloc[0].astype(str)
        has_drug = any(drug_of_header(v) for v in r0)
        has_meta = sum(bool(re.search(r"baseline|passage|breakpoint|critical|µg|ug/ml|mic|"
                                      r"sensitive|resistant", str(v), re.I)) for v in r0) >= 2
        if not (has_drug or has_meta):
            break
        df = df.copy()
        df.columns = [re.sub(r"\s+", " ", f"{c} {v}").strip()
                      for c, v in zip(df.columns.astype(str), r0)]
        df = df.iloc[1:].reset_index(drop=True)
    return df


def genes_muts_in_row(cells, caption_gene=None):
    """[(gène, mutation)] trouvés dans une ligne.

    Règle d'appariement : chaque mutation revient au gène qui la PRÉCÈDE
    IMMÉDIATEMENT (appariement 1:1), jamais à tous les gènes d'une fenêtre. Les
    cellules listent en effet plusieurs couples à la suite, parfois sans séparateur
    quand les sauts de ligne HTML sont perdus :
        « Rv0678 214dupC (1.0), fbiC Glu142Asp (0.1), fbiD 559_561delTGC (0.5) »
        « fgd1: Arg280CysfbiC: Val318IlerplC: Cys154Arg »
    Une règle de fenêtre attribuait 559_561delTGC à la fois à fbiC et à fbiD (faux
    positif constaté au premier passage). À défaut de gène précédent, on retombe sur
    le gène annoncé par la légende de la table."""
    out = []
    for cell in cells:
        s = unglue(str(cell))
        if s.strip().lower() in WT_TOKENS:
            continue
        hits = list(GENE_RE.finditer(s))
        muts = list(MUT_RE.finditer(s))
        for m in muts:
            prev = [g for g in hits if g.end() <= m.start()]
            if prev:
                g = max(prev, key=lambda x: x.end())
                if m.start() - g.end() > 30:      # trop loin : appariement non fiable
                    continue
                gene, grp = PANEL[g.group(1).lower()]
            elif caption_gene and not hits:
                gene, grp = caption_gene
            else:
                continue
            out.append((gene, grp, m.group(0).strip()))
    # dédup en conservant l'ordre
    seen, uniq = set(), []
    for t in out:
        if t not in seen:
            seen.add(t); uniq.append(t)
    return uniq


def caption_gene_of(text):
    m = GENE_RE.search(text or "")
    if not m:
        return None
    g = PANEL[m.group(1).lower()]
    return g


# --- collecte des tables ---------------------------------------------------
def jats_tables(pmcid):
    """[(label, caption, DataFrame)] depuis le fullTextXML en cache."""
    ft = CACHE / pmcid / "fulltext.xml"
    if not ft.exists():
        return []
    out = []
    try:
        root = etree.fromstring(ft.read_bytes())
    except Exception as e:
        print(f"    [jats] XML illisible {pmcid}: {e}")
        return []
    for tw in root.findall(".//table-wrap"):
        tab = tw.find(".//table")
        if tab is None:
            continue
        lab = (tw.findtext("label") or "?").strip()
        cap = " ".join(tw.find("caption").itertext()).strip() if tw.find("caption") is not None else ""
        try:
            html = etree.tostring(tab, encoding="unicode", method="html")
            df = pd.read_html(io.StringIO(html), header=0)[0]
        except Exception:
            continue
        out.append((lab, re.sub(r"\s+", " ", cap), merge_subheaders(df)))
    return out


def supp_tables(pmcid):
    """[(label, caption, DataFrame)] depuis les fichiers supplementary lisibles."""
    d = CACHE / pmcid / "supp"
    if not d.exists():
        return []
    out = []
    for f in sorted(d.rglob("*")):
        if not f.is_file() or f.suffix.lower() in (".gif", ".jpg", ".jpeg", ".png", ".tif"):
            continue
        ext = f.suffix.lower()
        try:
            if ext in (".xlsx", ".xls"):
                xl = pd.ExcelFile(f)
                for sh in xl.sheet_names:
                    out.append((f"{f.name}::{sh}", "", merge_subheaders(
                        pd.read_excel(f, sheet_name=sh, dtype=str))))
            elif ext in (".csv", ".tsv"):
                sep = "\t" if ext == ".tsv" else ","
                out.append((f.name, "", merge_subheaders(
                    pd.read_csv(f, sep=sep, dtype=str, on_bad_lines="skip"))))
            elif ext == ".docx":
                import docx
                for i, t in enumerate(docx.Document(str(f)).tables):
                    rows = [[c.text for c in r.cells] for r in t.rows]
                    if len(rows) > 1:
                        out.append((f"{f.name}::t{i}", "", merge_subheaders(
                            pd.DataFrame(rows[1:], columns=rows[0]))))
            elif ext == ".pdf":
                txt = subprocess.run(["pdftotext", "-layout", str(f), "-"],
                                     capture_output=True, text=True, timeout=120).stdout
                lines = [l for l in txt.splitlines() if GENE_RE.search(l) and MUT_RE.search(l)]
                if lines:
                    out.append((f"{f.name}::pdftext", "", pd.DataFrame({"_line": lines})))
        except Exception as e:
            print(f"    [supp] parse échoué {f.name}: {type(e).__name__} {e}")
    return out


def main():
    papers = []
    for cand in csv.DictReader(open(paths.DATA / "litcohort_candidates.tsv"), delimiter="\t"):
        pmid = cand["pmid"]
        core = CACHE / pmid / "core.json"
        if not core.exists():
            continue
        import json
        try:
            r = json.load(open(core))["resultList"]["result"][0]
        except Exception:
            continue
        if r.get("pmcid"):
            papers.append((pmid, r["pmcid"], r.get("title", "")[:70],
                           r.get("journalInfo", {}).get("journal", {}).get("title", ""),
                           cand.get("drug_focus", "")))

    evid, inv = [], []
    for pmid, pmcid, title, journal, focus in papers:
        tables = [("body", l, c, d) for l, c, d in jats_tables(pmcid)] + \
                 [("supp", l, c, d) for l, c, d in supp_tables(pmcid)]
        ft = CACHE / pmcid / "fulltext.xml"
        blob = ft.read_text(errors="replace") if ft.exists() else ""
        da_acc = sorted(set(ACC_RE.findall(blob)))
        da_prj = sorted(set(PROJ_RE.findall(blob)))
        # médicament unique du papier (titre), pour le repli « colonne MIC nue »
        tdrugs = {d for k, d in DRUG_LONG.items() if k in title.lower()}
        paper_drug = tdrugs.pop() if len(tdrugs) == 1 else None
        n_ev_paper = 0
        for src, lab, cap, df in tables:
            if df is None or df.empty:
                continue
            df = df.dropna(how="all")
            cg = caption_gene_of(cap)
            dcols = {c: d for c in df.columns if (d := drug_of_header(c))}
            # Repli : une colonne « MIC » NUE (sans nom de médicament) est attribuée au
            # médicament unique du papier, quand le titre n'en désigne qu'un. Cas ACS
            # Omega (ddn, prétomanide) : colonne « MIC », gène dans la légende.
            if paper_drug and not dcols:
                for c in df.columns:
                    if re.fullmatch(r"mic(?:\s*(?:µg|ug|mg)?\s*/?\s*(?:ml|l)?)?", norm_hdr(c)):
                        dcols[c] = paper_drug
            acc_in_tab = len(set(ACC_RE.findall(df.astype(str).to_string())))
            n_ev = 0
            for _, row in df.iterrows():
                cells = [row[c] for c in df.columns]
                gm = genes_muts_in_row(cells, cg)
                if not gm:
                    continue
                phen = []
                for c, drug in dcols.items():
                    rs, mic, raw = parse_pheno(row.get(c))
                    if rs or mic:
                        phen.append((drug, rs, mic, raw, str(c)[:60]))
                if not phen:
                    continue
                for gene, grp, mut in gm:
                    for drug, rs, mic, raw, col in phen:
                        evid.append(dict(
                            pmid=pmid, pmcid=pmcid, journal=journal, source=src,
                            table=lab, gene=gene, panel_group=grp, mutation=mut,
                            drug=drug, rs=rs, mic=mic, raw_value=raw, column=col,
                            caption=cap[:150]))
                        n_ev += 1
            n_ev_paper += n_ev
            inv.append(dict(pmid=pmid, pmcid=pmcid, source=src, table=lab,
                            rows=len(df), cols=len(df.columns),
                            drug_cols=";".join(sorted(set(dcols.values()))) or "-",
                            acc_in_table=acc_in_tab, evidence_rows=n_ev,
                            caption=cap[:110]))
        print(f"  {pmid} {pmcid:12s} {len(tables):2d} tables | preuves={n_ev_paper:3d} | "
              f"data-avail {len(da_acc)} SRA / {len(da_prj)} PRJ | {title}")

    ev = pd.DataFrame(evid).drop_duplicates()
    # Canal manuel : preuves curées à la main depuis les PDF paywall que CG fournit
    # (piste P3.5.b — ces papiers n'ont ni fullTextXML ni supp récupérables via EPMC).
    # Même schéma que l'extraction automatique, source = « manual_pdf ».
    if MANUAL_EVID.exists():
        man = pd.read_csv(MANUAL_EVID, sep="\t", dtype=str).fillna("")
        n_before = len(ev)
        ev = pd.concat([ev, man], ignore_index=True).drop_duplicates()
        print(f"\n  + {len(ev) - n_before} preuves manuelles (PDF fournis par CG) "
              f"depuis {MANUAL_EVID.name} : {', '.join(sorted(set(man.pmid)))}")
    iv = pd.DataFrame(inv)
    ev.to_csv(OUT_EVID, sep="\t", index=False)
    iv.to_csv(OUT_INV, sep="\t", index=False)

    print(f"\n=== {len(ev)} lignes de preuve, {len(iv)} tables inventoriées")
    if not ev.empty:
        print("\n--- preuves par gène du panel ---")
        print(ev.groupby(["panel_group", "gene"]).size().to_string())
        f420 = ev[ev.gene.isin(["ddn", "fgd1", "fbiA", "fbiB", "fbiC", "fbiD"])]
        print(f"\n--- voie F420 : {len(f420)} preuves, "
              f"{f420.mutation.nunique()} mutations distinctes ---")
        cols = ["pmid", "gene", "mutation", "drug", "rs", "mic", "table"]
        print(f420[cols].drop_duplicates().to_string(index=False))
    crosscheck(ev)
    print(f"\nSorties : {OUT_EVID.name}, {OUT_INV.name}, {OUT_XCHK.name}")


AA3_POS = re.compile(r"^(?:p\.)?(?:" + AA3 + r")(\d+)", re.I)
AA1_POS = re.compile(r"^(?:p\.)?[ACDEFGHIKLMNPQRSTVWY](\d+)")


def prot_pos(mut):
    """Position PROTÉIQUE d'une mutation, ou None si la mutation est en coordonnées
    NUCLÉOTIDIQUES (« 559_561delTGC », « 398_del_1_TC_T ») : les deux systèmes ne sont
    pas comparables et confondre les deux fabriquerait de faux appariements."""
    for rgx in (AA3_POS, AA1_POS):
        m = rgx.match(str(mut).strip())
        if m:
            return int(m.group(1))
    return None


def crosscheck(ev):
    """Les candidats du projet sont-ils phénotypés dans la littérature minée ?

    C'est la question qui décide de la valeur de cette phase pour le manuscrit :
    l'absence de DST sur les candidats convergents cesse d'être une hypothèse pour
    devenir un constat MESURÉ sur un corpus BPaL récent."""
    src = paths.RESULTATS / "phase5c_candidate_emergence.tsv"
    if not src.exists():
        return
    cand = pd.read_csv(src, sep="\t")
    lit = ev[ev.gene.isin(["ddn", "fgd1", "fbiA", "fbiB", "fbiC", "fbiD"])].copy()
    lit["prot_pos"] = lit.mutation.map(prot_pos)
    rows = []
    for _, c in cand.iterrows():
        hit = lit[(lit.gene == c.gene) & (lit.prot_pos == c.pos)]
        drugs = sorted(set(hit[hit.rs.isin(["R"]) | hit.mic.notna()].drug)) if not hit.empty else []
        rows.append(dict(
            candidate=f"{c.gene}:{c.pos}", classe=c["class"], n_carriers=c.n_carriers,
            lit_phenotyped="OUI" if not hit.empty else "non",
            lit_mutations=";".join(sorted(set(hit.mutation))) if not hit.empty else "",
            lit_drugs=";".join(d for d in drugs if d in ("pretomanid", "delamanid")),
            lit_pmids=";".join(sorted(set(hit.pmid.astype(str)))) if not hit.empty else ""))
    xc = pd.DataFrame(rows)
    xc.to_csv(OUT_XCHK, sep="\t", index=False)
    print("\n--- CROISEMENT candidats du projet × phénotypes de la littérature ---")
    print(xc.to_string(index=False))
    # variants F420 phénotypés que le projet n'a PAS dans ses candidats
    known = {(r.gene, r.pos) for _, r in cand.iterrows()}
    extra = {(g, int(p)) for g, p in zip(lit.gene, lit.prot_pos)
             if pd.notna(p) and (g, int(p)) not in known}
    if extra:
        print(f"\n--- positions F420 phénotypées ABSENTES des candidats du projet : "
              f"{', '.join(f'{g}:{int(p)}' for g, p in sorted(extra))}")


if __name__ == "__main__":
    main()
