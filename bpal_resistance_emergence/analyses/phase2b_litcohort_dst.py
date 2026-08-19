#!/usr/bin/env python3
"""Phase 2.5 — Cohorte DST par minage de littérature (étages B + C).

Étage A (découverte des papiers) est fait côté agent via tbmonitor et déposé dans
`data/litcohort_candidates.tsv` (pmid, drug_focus, access_hint, note). Ce script fait :

  Étage B — pour chaque PMID : Europe PMC core (pmcid, OA), fullTextXML (accessions de la
    déclaration « Availability of data » ET tables du CORPS d'article, cf. `jats_tables`),
    supplementaryFiles (zip) → parse des tables (.xlsx/.csv ; .docx/.pdf via
    `uvx markitdown` si dispo) → extraction heuristique de (accession, drogue+valeur,
    gène+mutation). Non-OA → fallback Unpaywall, puis routage vers `need_manual`
    (PDF à fournir par l'utilisateur).

PÉRIMÈTRE ET SCRIPT FRÈRE (établi le 2026-07-20). Ce script produit des paires
**accession × phénotype** (schéma `res.csv`). C'est atteignable pour la bédaquiline et
la clofazimine, mais PAS pour le prétomanide ni le délamanide : mesuré sur le corpus,
1 seule table sur 51 porte à la fois une colonne médicament et une accession, car le
DST par souche de ces molécules récentes est publié avec des identifiants internes
contre un dépôt BioProject global. Un rendement nul sur un papier PMD/DLM n'est donc
PAS un bug — voir `phase2c_litcohort_jats.py`, qui vise la granularité atteignable :
la preuve au niveau **variant** (gène, mutation, médicament, MIC/R-S, dispositif).

  Étage C — binarisation MIC→R/S (breakpoints documentés ; PMD sans breakpoint = R/S sur
    label explicite seulement, MIC conservée), normalisation d'ID, jointure inter-tables,
    déduplication contre `phenotypes_souches_consensus.tsv` (ne garder que l'incrément),
    map SRA→lignée (snapshot Coll) + validation des mutations panel contre le SPDI de
    `bdd/actuelle` (réutilise l'annotateur codon de phase1).

Sorties :
  data/litcohort_dst.csv                 schéma res.csv (SRA, link, <drogue>=R/S) -> ingest_literature
  résultats/litcohort_raw_extraction.tsv toutes les paires extraites (traçabilité)
  résultats/litcohort_validation.tsv     mutations panel rapportées vs SPDI bdd/actuelle
  résultats/litcohort_need_manual.tsv    papiers non-OA (PDF à fournir) + URL Unpaywall

Usage :
  python phase2b_litcohort_dst.py --candidates data/litcohort_candidates.tsv [--limit N]
  python phase2b_litcohort_dst.py --paper 33239092        # un seul PMID
  python phase2b_litcohort_dst.py --test-beckert          # régression (cache OK)
"""
import argparse, csv, json, os, re, shutil, subprocess, sys, tarfile, time, urllib.parse, urllib.request, zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths
from phase1_feasibility_scan import load_panel, load_fasta, annotate_effect

EMAIL = "christophe.guyeux@univ-fcomte.fr"
CACHE = paths.RESULTATS / "litcohort_cache"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
ACC_RE = re.compile(r"\b([DES]RR\d{4,})\b")
PROJ_RE = re.compile(r"\b(PRJ[END][ABJ]\d+|[DES]RP\d{5,})\b")
MUT_RE = re.compile(r"^[A-Z]\d{1,4}[A-Z*]$|del|ins|fs|\d+_\d+|>")   # token "mutation-like"

# Alias de drogues pour reconnaître les en-têtes de colonnes (codes courts = tokens).
DRUG_ALIASES = {
    "bedaquiline":  ["bedaquilin", "bdq"],
    "pretomanid":   ["pretomanid", "pa-824", "pa824", "ptm", "pmd"],
    "delamanid":    ["delamanid", "dlm"],
    "clofazimine":  ["clofazimin", "cfz", "cfm"],
    "linezolid":    ["linezolid", "lzd"],
    "moxifloxacin": ["moxifloxacin", "mfx", "mxf", "moxi"],
}
# Gènes du panel reconnus comme colonnes "mutation" (le contenu de la cellule = la mutation).
PANEL_GENES = {"rv0678": "Rv0678", "mmpr5": "Rv0678", "atpe": "atpE", "pepq": "pepQ",
               "ddn": "ddn", "fgd1": "fgd1", "fbia": "fbiA", "fbib": "fbiB", "fbic": "fbiC",
               "fbid": "fbiD", "rrl": "rrl", "rplc": "rplC", "gyra": "gyrA", "gyrb": "gyrB",
               "rpoc": "rpoC", "rpoa": "rpoA"}
GENERIC_MUT = {"mutation", "variant", "snp", "aachange", "aminoacid", "rv0678mutation"}
# Breakpoints MIC (mg/L) — APPROXIMATIFS, dépendants du milieu ; flag mic_binarized.
# PMD : pas de breakpoint WHO consolidé -> None (jamais binariser, garder la MIC).
BREAKPOINTS = {"bedaquiline": 1.0, "clofazimine": 1.0, "delamanid": 0.06,
               "linezolid": 1.0, "moxifloxacin": 1.0, "pretomanid": None}
# NB : "1"/"0" volontairement EXCLUS — ils collisionnent avec une MIC de 1 ou 0 mg/L.
# Une cellule "1" dans une table DST est presque toujours une MIC, pas un flag binaire R/S
# (les tables binaires 0/1 seraient mal binarisées, cas rare à signaler en revue manuelle).
RS_LABELS = {"r": "R", "resistant": "R", "res": "R", "resistance": "R",
             "s": "S", "susceptible": "S", "sensitive": "S"}

norm_id = lambda s: re.sub(r"[^0-9a-z]", "", str(s).lower())
norm_hdr = lambda s: re.sub(r"[^0-9a-z]", "", str(s).lower())


# ------------------------------------------------------------------ HTTP + cache
def _get(url, dest=None, binary=False, tries=3):
    if dest and dest.exists() and dest.stat().st_size > 0:
        return dest.read_bytes() if binary else dest.read_text(errors="replace")
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": f"bpal-litcohort ({EMAIL})"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
            if dest:
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(data)
            return data if binary else data.decode("utf-8", "replace")
        except Exception as e:
            if k == tries - 1:
                print(f"    [http] échec {url[:70]} : {e}")
                return None
            time.sleep(1.0 + k)


def epmc_core(pmid):
    q = urllib.parse.quote(f"EXT_ID:{pmid} AND SRC:MED")
    txt = _get(f"{EPMC}/search?query={q}&format=json&resultType=core",
               CACHE / pmid / "core.json")
    if not txt:
        return {}
    try:
        L = json.loads(txt).get("resultList", {}).get("result", [])
        return L[0] if L else {}
    except Exception:
        return {}


def unpaywall(doi):
    if not doi:
        return None
    txt = _get(f"https://api.unpaywall.org/v2/{urllib.parse.quote(doi)}?email={EMAIL}",
               CACHE / "unpaywall" / (norm_id(doi) + ".json"))
    try:
        loc = (json.loads(txt) or {}).get("best_oa_location") or {}
        return loc.get("url_for_pdf") or loc.get("url")
    except Exception:
        return None


# ------------------------------------------------------------------ étage B : parse tables
def fetch_supp_dir(pmcid):
    """Dossier des fichiers supplementary. 1) EBI supplementaryFiles (zip) ; si EBI renvoie
    du HTML/500, 2) bascule sur le PACKAGE OA NCBI (oa.fcgi -> tar.gz), souvent UP quand EBI
    est en panne. Valide is_zipfile/is_tarfile, purge le cache empoisonné, retry. None si rien."""
    d = CACHE / pmcid / "supp"
    if d.exists() and any(d.rglob("*")):
        return d
    d.parent.mkdir(parents=True, exist_ok=True)
    # 1) EBI supplementaryFiles (zip)
    zf = CACHE / pmcid / "supp.zip"
    if zf.exists() and not zipfile.is_zipfile(zf):
        zf.unlink()
    for attempt in range(2):
        if zf.exists() and zipfile.is_zipfile(zf):
            break
        try:
            subprocess.run(["curl", "-sL", "--max-time", "240", "-o", str(zf),
                            f"{EPMC}/{pmcid}/supplementaryFiles"], timeout=260, check=False)
        except Exception:
            pass
        if zf.exists() and not zipfile.is_zipfile(zf):
            zf.unlink(); time.sleep(2 * (attempt + 1))
    if zf.exists() and zipfile.is_zipfile(zf):
        try:
            with zipfile.ZipFile(zf) as z:
                z.extractall(d)
            return d
        except Exception:
            pass
    # 2) Fallback NCBI PMC OA package (tar.gz)
    oa = _get(f"https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmcid}", CACHE / pmcid / "oa.xml")
    m = re.search(r'href="(ftp://\S+?\.tar\.gz)"', oa or "")
    if m:
        url = m.group(1).replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
        tgz = CACHE / pmcid / "oa.tar.gz"
        if not (tgz.exists() and tarfile.is_tarfile(tgz)):
            try:
                subprocess.run(["curl", "-sL", "--max-time", "240", "-o", str(tgz), url],
                               timeout=260, check=False)
            except Exception:
                pass
        if tgz.exists() and tarfile.is_tarfile(tgz):
            try:
                with tarfile.open(tgz) as t:
                    t.extractall(d)
                print(f"    [supp] récupéré via NCBI OA package (fallback EBI 500) pour {pmcid}")
                return d
            except Exception:
                pass
    print(f"    [supp] aucun supp récupérable pour {pmcid} (EBI 500 + NCBI OA indispo) -> à re-tenter")
    return None


def supp_tables(pmcid):
    """{nom_fichier::feuille: DataFrame} depuis le dossier supp (EBI zip ou NCBI OA tgz)."""
    out = {}
    d = fetch_supp_dir(pmcid)
    if d is None:
        return out
    for f in sorted(d.rglob("*")):
        if not f.is_file():
            continue
        ext = f.suffix.lower()
        try:
            if ext in (".xlsx", ".xls"):
                xl = pd.ExcelFile(f)
                for sh in xl.sheet_names:
                    out[f"{f.name}::{sh}"] = pd.read_excel(f, sheet_name=sh, dtype=str)
            elif ext in (".csv", ".tsv"):
                sep = "\t" if ext == ".tsv" else ","
                out[f.name] = pd.read_csv(f, sep=sep, dtype=str, on_bad_lines="skip")
            elif ext in (".docx", ".pdf") and shutil.which("uvx") and \
                    re.search(r"supp|mmc|additional|electronic|_s\d|table[_\s]?s?\d", f.name.lower()):
                # restreint aux fichiers SUPPLEMENTARY (le package OA NCBI inclut le PDF principal)
                md = subprocess.run(["uvx", "markitdown", str(f)], capture_output=True,
                                    text=True, timeout=120).stdout
                # extraire les tables markdown grossièrement -> DataFrame par bloc de lignes |...|
                rows = [ln for ln in md.splitlines() if ln.count("|") >= 2]
                if rows:
                    out[f.name + "::md"] = pd.DataFrame({"_raw": rows})
        except Exception as e:
            print(f"    [supp] parse échoué {f.name}: {e}")
    return out


def _merge_subheaders(df, max_rounds=3):
    """Fusionne les en-têtes JATS étalés sur 2-3 lignes : la ligne 0 du CORPS porte
    souvent les vrais noms de médicaments sous un en-tête générique « pDST (MIC) »."""
    for _ in range(max_rounds):
        if df.empty:
            break
        r0 = df.iloc[0].astype(str)
        hdrish = any(any(a in re.sub(r"[^a-z0-9]+", " ", str(v).lower())
                         for al in DRUG_ALIASES.values() for a in al) for v in r0) or \
            sum(bool(re.search(r"baseline|passage|breakpoint|critical|µg|ug/ml|mic", str(v), re.I))
                for v in r0) >= 2
        if not hdrish:
            break
        df = df.copy()
        df.columns = [re.sub(r"\s+", " ", f"{c} {v}").strip()
                      for c, v in zip(df.columns.astype(str), r0)]
        df = df.iloc[1:].reset_index(drop=True)
    return df


def jats_tables(pmcid):
    """{'JATS::<label>': DataFrame} — tables du CORPS de l'article (fullTextXML).

    POURQUOI (défaut corrigé le 2026-07-20, cf. cahier). Ce script ne lisait que les
    fichiers SUPPLEMENTARY. Or pour Frontiers, MDPI, BMC et AAC, l'endpoint EPMC
    `/supplementaryFiles` ne renvoie QUE des images de figures : les tables de DST
    sont dans le corps de l'article, balises `<table-wrap>`. Le batch du 12/06 avait
    donc un rendement quasi nul (« Beckert seul »), imputé à tort à une panne EBI —
    l'endpoint était bien en 500, mais le vrai verrou était cet angle mort. Mesure
    faite sur le corpus : 40 tables exploitables dans le corps contre 11 en supp.

    NB : `pd.read_html` appliqué DIRECTEMENT au fichier JATS échoue (« No tables
    found ») ; il faut extraire le `<table>` puis le sérialiser en HTML."""
    ft = CACHE / pmcid / "fulltext.xml"
    if not ft.exists():
        return {}
    try:
        import io as _io
        from lxml import etree
    except ImportError:
        print("    [jats] lxml absent -> tables du corps non lues")
        return {}
    out = {}
    try:
        root = etree.fromstring(ft.read_bytes())
    except Exception as e:
        print(f"    [jats] XML illisible {pmcid}: {e}")
        return {}
    for tw in root.findall(".//table-wrap"):
        tab = tw.find(".//table")
        if tab is None:
            continue
        lab = (tw.findtext("label") or "?").strip()
        try:
            html = etree.tostring(tab, encoding="unicode", method="html")
            df = pd.read_html(_io.StringIO(html), header=0)[0]
        except Exception:
            continue
        out[f"JATS::{lab}"] = _merge_subheaders(df)
    return out


def classify_cols(df):
    cols = list(df.columns)
    id_col = acc_col = None
    drug_cols, mut_cols = {}, {}
    for c in cols:
        h = norm_hdr(c)
        # accession : valeurs SRR/ERR/DRR
        if acc_col is None and df[c].astype(str).map(lambda v: bool(ACC_RE.search(v))).any():
            acc_col = c
        for drug, al in DRUG_ALIASES.items():
            if any(a.replace("-", "") in h for a in al) and drug not in drug_cols:
                drug_cols[c] = drug
        if h in PANEL_GENES:
            mut_cols[c] = PANEL_GENES[h]
        elif h in GENERIC_MUT:
            mut_cols[c] = "?"
    for c in cols:                                   # id : en-tête explicite, sinon 1re colonne
        h = norm_hdr(c)
        if h in ("id", "isolate", "strain", "sample", "sampleid", "strainid", "originalid",
                 "isolateid", "id1", "run", "ena", "enarun"):
            id_col = c; break
    if id_col is None and cols:
        id_col = cols[0]
    return id_col, acc_col, drug_cols, mut_cols


def extract_paper(pmid, core):
    """Renvoie (rows, dataavail_accessions). rows = liste de dicts d'extraction."""
    pmcid = core.get("pmcid")
    rows = []
    da_acc = set()
    if pmcid:
        ft = _get(f"{EPMC}/{pmcid}/fullTextXML", CACHE / pmcid / "fulltext.xml")
        if ft:
            da_acc = set(PROJ_RE.findall(ft)) | set(ACC_RE.findall(ft))
        tables = supp_tables(pmcid)
        tables.update(jats_tables(pmcid))     # + tables du CORPS d'article (cf. jats_tables)
    else:
        tables = {}
    # collecte
    id2acc = defaultdict(set)        # norm_id -> accessions
    pheno = []                        # (norm_id|None, accession|None, drug, value, sheet)
    muts = []                         # (norm_id|None, accession|None, gene, mutation, sheet)
    for sheet, df in tables.items():
        if df is None or df.empty:
            continue
        df = df.dropna(how="all")
        idc, accc, dcols, mcols = classify_cols(df)
        for _, r in df.iterrows():
            rid = norm_id(r[idc]) if idc and pd.notna(r.get(idc)) else None
            racc = None
            if accc and pd.notna(r.get(accc)):
                m = ACC_RE.search(str(r[accc]))
                racc = m.group(1) if m else None
            if rid and racc:
                id2acc[rid].add(racc)
            for c, drug in dcols.items():
                v = r.get(c)
                if pd.isna(v) or str(v).strip() in ("", "-", "nan", "NA", "ND"):
                    continue
                pheno.append((rid, racc, drug, str(v).strip(), sheet))
            for c, gene in mcols.items():
                v = r.get(c)
                # marqueurs wild-type / référence = PAS une mutation (gWT, H37Rv dans la table S3)
                if pd.isna(v) or str(v).strip().lower() in (
                        "", "-", "nan", "na", "nd", "wt", "gwt", "h37rv", "ref",
                        "reference", "none", "wild-type", "wildtype", "0"):
                    continue
                muts.append((rid, racc, gene, str(v).strip(), sheet))
    # jointure par norm_id quand l'accession n'est pas dans la même ligne
    def resolve(rid, racc):
        if racc:                              # accession dans la même ligne -> sûr
            return {racc}
        if not rid or len(rid) < 4:           # id trop court/vide -> pas de jointure (collisions)
            return set()
        accs = id2acc.get(rid, set())
        return accs if len(accs) == 1 else set()   # jointure 1:1 SEULEMENT (évite le cross-product)
    for rid, racc, drug, val, sheet in pheno:
        for acc in resolve(rid, racc):
            rows.append(dict(pmid=pmid, pmcid=pmcid or "", accession=acc, kind="pheno",
                             key=drug, value=val, sheet=sheet))
    for rid, racc, gene, mut, sheet in muts:
        for acc in resolve(rid, racc):
            # ne garder que des tokens mutation-like (évite le bruit des colonnes génériques)
            if gene != "?" or MUT_RE.search(mut):
                rows.append(dict(pmid=pmid, pmcid=pmcid or "", accession=acc, kind="mutation",
                                 key=gene, value=mut, sheet=sheet))
    return rows, da_acc


# ------------------------------------------------------------------ étage C : binarisation + RS
def to_rs(drug, raw):
    """(phenotype R/S|'', mic|'', source). Label R/S explicite prioritaire, sinon MIC.

    GARDE-FOU : une MIC doit être un nombre PUR (option. comparateur/unité). On REJETTE
    tout token contenant des lettres alphanumériques (ex. 'Rv0678 M146T', 'gyrA D94N')
    — sinon une mutation comme M146T donnerait un faux 146 → faux R (bug vécu, Table S2
    = table de mutations PAR drogue, pas de MIC)."""
    s = str(raw).strip().lower()
    if s in RS_LABELS:
        return RS_LABELS[s], "", "label"
    s2 = re.sub(r"(mg/?l|µg/?ml|ug/?ml|µg/?l|µg|\s)", "", s).replace("≤", "<=").replace("≥", ">=").replace(",", ".")
    if not re.fullmatch(r"[<>]?=?\d*\.?\d+", s2):     # MIC pure uniquement
        return "", "", ""
    mic = float(re.search(r"\d*\.?\d+", s2).group())
    bp = BREAKPOINTS.get(drug)
    if bp is None:                       # PMD : pas de breakpoint -> garder MIC, pas de R/S auto
        return "", str(mic), "mic_only"
    return ("R" if mic > bp else "S"), str(mic), "mic_binarized"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates")
    ap.add_argument("--paper")
    ap.add_argument("--test-beckert", action="store_true")
    ap.add_argument("--limit", type=int)
    args = ap.parse_args()
    paths.ensure_dirs(); CACHE.mkdir(parents=True, exist_ok=True)

    if args.test_beckert:
        pmids = [("33239092", "bedaquiline", "oa")]
    elif args.paper:
        pmids = [(args.paper, "?", "?")]
    else:
        rows = list(csv.DictReader(open(args.candidates), delimiter="\t"))
        pmids = [(r["pmid"], r.get("drug_focus", "?"), r.get("access_hint", "?")) for r in rows]
    if args.limit:
        pmids = pmids[:args.limit]

    all_rows, need_manual = [], []
    for i, (pmid, focus, hint) in enumerate(pmids, 1):
        core = epmc_core(pmid)
        oa = core.get("isOpenAccess") == "Y" and core.get("pmcid")
        title = (core.get("title") or "")[:60]
        print(f"[{i}/{len(pmids)}] {pmid} OA={'Y' if oa else 'N'} pmcid={core.get('pmcid')} | {title}")
        if not oa:
            need_manual.append(dict(pmid=pmid, doi=core.get("doi", ""), focus=focus,
                                    journal=(core.get("journalInfo") or {}).get("journal", {}).get("title", ""),
                                    unpaywall=unpaywall(core.get("doi")) or "", title=title))
            continue
        rows, da = extract_paper(pmid, core)
        for r in rows:
            r["drug_focus"] = focus
        all_rows.extend(rows)
        accs = sorted({r["accession"] for r in rows})
        print(f"      -> {len(rows)} paires, {len(accs)} accessions, {len(da)} accessions data-availability")
        time.sleep(0.2)

    # ---- raw extraction
    raw_path = paths.RESULTATS / "litcohort_raw_extraction.tsv"
    with open(raw_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pmid", "pmcid", "drug_focus", "accession", "kind", "key", "value", "sheet"], delimiter="\t")
        w.writeheader()
        for r in all_rows:
            w.writerow({k: r.get(k, "") for k in w.fieldnames})

    # ---- need_manual
    nm_path = paths.RESULTATS / "litcohort_need_manual.tsv"
    with open(nm_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pmid", "doi", "focus", "journal", "unpaywall", "title"], delimiter="\t")
        w.writeheader(); w.writerows(need_manual)

    # ---- agrégation par accession : phénotype R/S par drogue + MIC + mutations
    per = defaultdict(lambda: {"drugs": {}, "mic": {}, "muts": set()})
    for r in all_rows:
        acc = r["accession"]
        if r["kind"] == "pheno":
            rs, mic, src = to_rs(r["key"], r["value"])
            if rs and r["key"] not in per[acc]["drugs"]:
                per[acc]["drugs"][r["key"]] = rs
            if mic:
                per[acc]["mic"][r["key"]] = mic
        else:
            per[acc]["muts"].add((r["key"], r["value"]))

    # ---- dédup contre le consensus (par paire accession×drogue)
    consensus = {}
    if paths.PHENOTYPES_TSV.exists():
        cdf = pd.read_csv(paths.PHENOTYPES_TSV, sep="\t", low_memory=False)
        cdf["drug"] = cdf["drug"].astype(str).str.lower()
        consensus = cdf.groupby("strain_id")["drug"].agg(set).to_dict()   # vectorisé (évite iterrows 317k)

    # ---- res.csv (incrément seulement)
    drugs_order = list(DRUG_ALIASES.keys())
    res_path = paths.DATA / "litcohort_dst.csv"
    n_inc_pairs = 0
    with open(res_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["SRA", "link"] + drugs_order)
        for acc, d in sorted(per.items()):
            row_drugs = {}
            for drug, rs in d["drugs"].items():
                if drug not in consensus.get(acc, set()):     # incrément uniquement
                    row_drugs[drug] = rs; n_inc_pairs += 1
            if not row_drugs:
                continue
            w.writerow([acc, ""] + [row_drugs.get(dr, "") for dr in drugs_order])

    # ---- validation des mutations panel contre le SPDI bdd/actuelle
    pos2gene, genes = load_panel()
    fasta = load_fasta(paths.H37RV_FASTA)
    bdd = Path(paths.BDD)
    idx = {}
    for clade in os.scandir(bdd):
        if clade.is_dir():
            try:
                for sra in os.scandir(clade.path):
                    if sra.is_dir():
                        idx.setdefault(sra.name, f"{clade.path}/{sra.name}/NC_000962.3/spdi.txt")
            except OSError:
                pass

    def bdd_panel_aas(sra, locus):
        p = idx.get(sra)
        if not p or not os.path.exists(p):
            return None
        g = genes[locus]; lo, hi = g["start"] - 1, g["end"] - 1
        seen = []
        for line in open(p):
            line = line.strip()
            if not line:
                continue
            pr = line.split(":")
            try:
                pos = int(pr[1])
            except (IndexError, ValueError):
                continue
            if lo <= pos <= hi:
                _, _, aa, _ = annotate_effect(genes, locus, pos, pr[2], pr[3], fasta)
                if aa.startswith("p."):
                    seen.append(aa[2:])
        return seen

    val_path = paths.RESULTATS / "litcohort_validation.tsv"
    n_val_ok = n_val_tot = 0
    with open(val_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "gene", "reported", "in_bdd", "bdd_panel_variants", "match"])
        for acc, d in sorted(per.items()):
            for gene, mut in sorted(d["muts"]):
                if gene not in genes:            # seulement les gènes du panel
                    continue
                seen = bdd_panel_aas(acc, gene)
                if seen is None:
                    w.writerow([acc, gene, mut, "no", "", ""]); continue
                n_val_tot += 1
                clean = re.sub(r"[^A-Za-z0-9*]", "", mut)
                ok = clean in seen
                n_val_ok += ok
                w.writerow([acc, gene, mut, "yes", ";".join(seen), "OK" if ok else "MISMATCH"])

    # ---- rapport
    print("\n=== RÉCAP Phase 2.5 ===")
    print(f"papiers traités : {len(pmids)} | OA exploités : {len(pmids)-len(need_manual)} | non-OA (need_manual) : {len(need_manual)}")
    print(f"accessions distinctes avec extraction : {len(per)}")
    print(f"paires accession×drogue d'incrément (vs consensus) écrites -> {res_path.name} : {n_inc_pairs}")
    bydrug = defaultdict(int)
    for acc, d in per.items():
        for dr in d["drugs"]:
            bydrug[dr] += 1
    print("phénotypes R/S par drogue (avant dédup) :", dict(bydrug))
    print(f"validation mutations panel (souches dans bdd) : {n_val_ok}/{n_val_tot} OK -> {val_path.name}")
    print(f"need_manual (PDF à fournir) : {len(need_manual)} -> {nm_path.name}")
    if args.test_beckert:
        ok = len(per) >= 30 and n_val_ok >= 3
        print(f"\n[RÉGRESSION BECKERT] {'PASS' if ok else 'FAIL'} "
              f"(attendu : >=30 accessions, >=3 validations Rv0678 OK)")


if __name__ == "__main__":
    main()
