#!/usr/bin/env python3
"""bib-check Phase 3 — vérification anti-hallucination des entrées BibTeX via CrossRef (par DOI).

Source primaire déterministe (≫ WebSearch) : https://api.crossref.org/works/{doi}. Pour chaque
entrée du .bib avec un DOI, compare title / year / 1er auteur (famille) / journal / volume / issue /
pages à la vérité CrossRef, et classe OK / MISMATCH(champs) / NOT_FOUND. Écrit un rapport TSV et
imprime les écarts. N'édite PAS le .bib (corrections appliquées à la main après revue).

Usage : python bibcheck_crossref.py [--bib ../article/references.bib] [--only-unverified]
"""
import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

UA = "bibcheck/1.0 (mailto:guyeux@gmail.com)"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def strip_braces(s):
    # enlève les accolades LaTeX et escapes simples {\c{c}}->c, {\"o}->o pour comparaison
    s = re.sub(r"\{\\[a-zA-Z]\{(\w)\}\}", r"\1", s)
    s = re.sub(r"\{\\.(\w)\}", r"\1", s)
    s = s.replace("{", "").replace("}", "").replace("\\", "")
    return s


def parse_bib(text):
    """parse minimal : retourne [{key, type, fields:{}}] (gère champs multi-lignes, accolades imbriquées)."""
    entries = []
    i = 0
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,]+),", text):
        typ, key = m.group(1), m.group(2).strip()
        # corps = du début après la 1re virgule jusqu'à l'accolade fermante équilibrée
        start = m.end()
        depth = 1
        j = m.start(0) + text[m.start(0):].index("{")
        depth = 0
        k = j
        while k < len(text):
            if text[k] == "{":
                depth += 1
            elif text[k] == "}":
                depth -= 1
                if depth == 0:
                    break
            k += 1
        body = text[start:k]
        fields = {}
        for fm in re.finditer(r"(\w+)\s*=\s*", body):
            fname = fm.group(1).lower()
            p = fm.end()
            if p >= len(body):
                break
            if body[p] == "{":
                d = 0
                q = p
                while q < len(body):
                    if body[q] == "{":
                        d += 1
                    elif body[q] == "}":
                        d -= 1
                        if d == 0:
                            break
                    q += 1
                val = body[p + 1:q]
            else:
                mm = re.match(r"([^,]+)", body[p:])
                val = mm.group(1).strip() if mm else ""
            fields[fname] = val.strip()
        entries.append({"key": key, "type": typ.lower(), "fields": fields})
    return entries


def crossref(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=20) as r:
        return json.load(r)["message"]


def first_family(author_field):
    # BibTeX "Last, First and Last, First" -> 1re famille ; ou "First Last and ..."
    a = strip_braces(author_field).split(" and ")[0].strip()
    if "," in a:
        return a.split(",")[0].strip()
    return a.split()[-1] if a else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bib", default=str(Path(__file__).resolve().parent.parent / "article/references.bib"))
    ap.add_argument("--only-unverified", action="store_true")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent.parent / "résultats/bibcheck_crossref.tsv"))
    args = ap.parse_args()

    text = Path(args.bib).read_text(encoding="utf-8")
    entries = parse_bib(text)
    todo = [e for e in entries if e["fields"].get("doi")
            and not (args.only_unverified and "verified" in e["fields"])]
    print(f"{len(entries)} entrées, {len(todo)} avec DOI à vérifier (CrossRef)\n")

    rows = []
    n_ok = n_mis = n_nf = 0
    for idx, e in enumerate(todo, 1):
        f = e["fields"]
        doi = f["doi"].strip()
        try:
            cr = crossref(doi)
        except Exception as ex:
            n_nf += 1
            rows.append((e["key"], doi, "NOT_FOUND", f"crossref error: {ex}"))
            print(f"[{idx}/{len(todo)}] {e['key']:34} NOT_FOUND ({ex})")
            time.sleep(0.2)
            continue
        problems = []
        # titre
        ct = (cr.get("title") or [""])[0]
        if norm(strip_braces(f.get("title", "")))[:60] != norm(ct)[:60]:
            # tolère sous-titres : exige inclusion du titre bib dans crossref ou inverse
            a, b = norm(strip_braces(f.get("title", ""))), norm(ct)
            if not (a in b or b in a or a[:40] == b[:40]):
                problems.append(f"title: bib='{strip_braces(f.get('title',''))[:45]}' cr='{ct[:45]}'")
        # année
        cy = None
        for kk in ("published", "published-print", "published-online", "issued"):
            if cr.get(kk, {}).get("date-parts"):
                cy = cr[kk]["date-parts"][0][0]
                break
        if f.get("year") and cy and str(cy) != str(f["year"]):
            problems.append(f"year: bib={f['year']} cr={cy}")
        # 1er auteur
        crfam = ""
        if cr.get("author"):
            crfam = cr["author"][0].get("family", "")
        bibfam = first_family(f.get("author", ""))
        if bibfam and crfam and norm(bibfam) != norm(crfam):
            problems.append(f"author1: bib='{bibfam}' cr='{crfam}'")
        # journal
        cj = (cr.get("container-title") or [""])[0]
        if f.get("journal") and cj and norm(f["journal"])[:12] not in norm(cj) and norm(cj)[:12] not in norm(f["journal"]):
            problems.append(f"journal: bib='{f['journal'][:30]}' cr='{cj[:30]}'")
        # volume
        if f.get("volume") and cr.get("volume") and str(f["volume"]) != str(cr["volume"]):
            problems.append(f"volume: bib={f['volume']} cr={cr['volume']}")
        # issue
        if f.get("number") and cr.get("issue") and str(f["number"]) != str(cr["issue"]):
            problems.append(f"issue: bib={f['number']} cr={cr['issue']}")
        # pages / article-number
        cr_page = cr.get("page") or cr.get("article-number") or ""
        if f.get("pages") and cr_page:
            bp = norm(f["pages"]); cp = norm(cr_page)
            if bp != cp and bp not in cp and cp not in bp:
                problems.append(f"pages: bib='{f['pages']}' cr='{cr_page}'")
        if problems:
            n_mis += 1
            rows.append((e["key"], doi, "MISMATCH", " | ".join(problems)))
            print(f"[{idx}/{len(todo)}] {e['key']:34} MISMATCH  {' | '.join(problems)}")
        else:
            n_ok += 1
            rows.append((e["key"], doi, "OK", ""))
        time.sleep(0.2)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("key\tdoi\tstatus\tproblems\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"\n=== {n_ok} OK | {n_mis} MISMATCH | {n_nf} NOT_FOUND === -> {args.out}")


if __name__ == "__main__":
    main()
