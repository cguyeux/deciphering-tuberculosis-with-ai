#!/usr/bin/env python3
"""
P1.5 — Quel est le SITE de phosphorylation exact de Rv0810c ?

Malakar et al. 2023 (mBio, PMID 37791794) rangent Rv0810c parmi treize proteines
« phosphorylees supplementaires » identifiees par phospho-secretome, sans nommer
de residu dans le corps du texte. Les residus sont dans les Data Sets S2, S3, S7
et S8, telecharges via Europe PMC.

PIEGE MAJEUR, et c'est la raison d'etre du controle ci-dessous : le Data Set S8
est un onglet ou DEUX listes de longueurs differentes ont ete collees cote a cote.
La colonne A (Rv_number) est la liste des 38 proteines uniques ; les colonnes B-G
sont la liste des 57 phosphopeptides. Lues ligne a ligne, elles s'apparient
faussement : la ligne « Rv0810c » y porte un peptide en positions 801-829, ce qui
est impossible pour une proteine de 60 aa. La colonne G (Protein_phospho_site)
est la seule qui soit interne a la table des peptides.

Le script tranche par la SEQUENCE : un phosphopeptide n'est attribue a Rv0810c que
s'il se retrouve exactement dans sa sequence, aux positions annoncees, et que le
residu modifie y est bien un S, T ou Y. Puis il croise les quatre jeux et la
comparaison inter-etudes du Data Set S3b (Fortuin, Prisic, Verma, Carette), et
situe chaque site dans l'architecture bipartite etablie en P2.1.

Sortie : résultats/p1_5_phosphosites.json
"""

import json
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
DATA = PROJ / "experiments" / "2026-08-10_P1_4_P1_5" / "data"
OUT = PROJ / "résultats" / "p1_5_phosphosites.json"
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

SEQ = "MGRGRAKAKQTKVARELKYSSPQTDFQRLQRELSGTGTDRLDGDGPSDDDSWNDEDDWRR"

# Architecture etablie en P2.1 (pLDDT par residu sur AF-I6XWB9-F1-model_v6).
SEGMENTS = [(1, 19, "module invariant (pLDDT 91,9)"),
            (20, 33, "charniere du module (pLDDT 85,6)"),
            (34, 60, "queue acide desordonnee (pLDDT 62,0)")]

FILES = {
    "S1_proteines_phosphorylees": "mbio.01232-23-s0001.xlsx",
    "S2_phosphopeptides_lysat": "mbio.01232-23-s0002.xlsx",
    "S3_comparaison_etudes": "mbio.01232-23-s0003.xlsx",
    "S4_secretome": "mbio.01232-23-s0004.xlsx",
    "S5_secretome_multi_etudes": "mbio.01232-23-s0005.xlsx",
    "S7_phosphoproteines_filtrat": "mbio.01232-23-s0007.xlsx",
    "S8_phosphosites_filtrat": "mbio.01232-23-s0008.xlsx",
}


def xlsx_sheets(path):
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        r = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in r.findall(NS + "si"):
            shared.append("".join(t.text or "" for t in si.iter(NS + "t")))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    names = [s.get("name") for s in wb.iter(NS + "sheet")]
    out = {}
    for i, name in enumerate(names, 1):
        fn = f"xl/worksheets/sheet{i}.xml"
        if fn not in z.namelist():
            continue
        root = ET.fromstring(z.read(fn))
        rows = []
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
                vals.append(val.strip())
            rows.append(vals)
        out[name or f"sheet{i}"] = rows
    return out


def segment_of(pos):
    for a, b, label in SEGMENTS:
        if a <= pos <= b:
            return label
    return "hors sequence"


def verify_peptide(peptide, start, stop, mod, seq=SEQ):
    """Le peptide est-il reellement celui de Rv0810c, aux positions annoncees ?"""
    checks = {"peptide": peptide, "start_annonce": start, "stop_annonce": stop,
              "modification": mod}
    idx = seq.find(peptide) if peptide else -1
    checks["present_dans_la_sequence"] = idx >= 0
    checks["position_reelle"] = (idx + 1) if idx >= 0 else None
    checks["positions_coherentes"] = (idx >= 0 and start is not None
                                      and idx + 1 == start
                                      and stop == start + len(peptide) - 1)
    m = re.match(r"([STY])(\d+)", mod or "")
    if m and idx >= 0:
        aa, off = m.group(1), int(m.group(2))
        site = idx + off               # position 1-based dans la proteine
        checks["site_calcule"] = f"{aa}{site}"
        checks["residu_a_cette_position"] = seq[site - 1] if site <= len(seq) else None
        checks["residu_concorde"] = (seq[site - 1] == aa) if site <= len(seq) else False
        checks["segment"] = segment_of(site)
        checks["site_position"] = site
        checks["site_aa"] = aa
    else:
        checks["site_calcule"] = None
        checks["residu_concorde"] = False
    checks["tryptique"] = bool(peptide) and (
        (idx <= 0 or seq[idx - 1] in "KR") and peptide[-1] in "KR")
    return checks


def conservation(synth, sto=None):
    """Conservation des positions phosphorylees dans l'alignement Pfam complet
    de DUF3073. Le rang est donne parmi les 59 positions du domaine : sans ce
    denominateur, « 81 % conserve » ne dit pas si la position est remarquable."""
    from collections import Counter
    sto = sto or (PROJ / "data" / "PF11273_full.sto")
    if not sto.exists():
        return {"erreur": f"alignement absent : {sto}"}
    seqs, sq_declared = {}, None
    for l in sto.read_text(encoding="utf-8", errors="replace").splitlines():
        if l.startswith("#=GF SQ"):
            sq_declared = int(l.split()[-1])
        if not l or l.startswith("#") or l.strip() == "//":
            continue
        p = l.split(None, 1)
        if len(p) == 2:
            seqs[p[0]] = p[1]
    ref_id = next((k for k in seqs if k.startswith("I6XWB9")), None)
    if ref_id is None:
        return {"erreur": "Rv0810c (I6XWB9) absent de l'alignement"}
    ref = seqs[ref_id]
    first = int(ref_id.split("/")[1].split("-")[0])
    col_of, pos = {}, first
    for i, ch in enumerate(ref):
        if ch not in ".-":
            col_of[pos] = i
            pos += 1

    def profile(p):
        c = col_of[p]
        col = [s[c] for s in seqs.values() if len(s) > c]
        gaps = sum(1 for x in col if x in ".-")
        aas = [x.upper() for x in col if x not in ".-"]
        cnt = Counter(aas)
        n = len(aas)
        top = cnt.most_common(3)
        return {
            "position": p, "residu_Rv0810c": ref[c].upper(), "colonne": c,
            "n_sequences_alignees": n, "n_gaps": gaps,
            "taux_gaps_pct": round(100 * gaps / len(col), 2),
            "identite_max_pct": round(100 * top[0][1] / n, 1) if n else None,
            "residu_majoritaire": top[0][0] if n else None,
            "top3": [[a, c_, round(100 * c_ / n, 1)] for a, c_ in top],
            "phosphorylable_STY_pct": round(
                100 * sum(cnt.get(a, 0) for a in "STY") / n, 1) if n else None,
        }

    profils = {p: profile(p) for p in sorted(col_of)}
    ident = sorted((v["identite_max_pct"] for v in profils.values()), reverse=True)
    gaps = sorted(v["taux_gaps_pct"] for v in profils.values())
    out = {"alignement": sto.name, "n_sequences": len(seqs),
           "n_sequences_declarees_GF_SQ": sq_declared,
           "n_positions_mappees": len(col_of), "sites": {}}
    for site, e in synth.items():
        p = e["position"]
        if p not in profils:
            out["sites"][site] = {"erreur": "position hors alignement"}
            continue
        pr = profils[p]
        out["sites"][site] = {
            **pr,
            "rang_identite": ident.index(pr["identite_max_pct"]) + 1,
            "n_positions": len(ident),
            "percentile_identite": round(
                100 * (1 - ident.index(pr["identite_max_pct"]) / len(ident)), 1),
            "rang_absence_de_gaps": gaps.index(pr["taux_gaps_pct"]) + 1,
        }
    # profil complet, pour situer les sites dans le domaine
    out["profil_par_position"] = profils
    out["positions_sans_aucun_gap"] = [p for p, v in profils.items() if v["n_gaps"] == 0]
    out["Thr_de_la_proteine"] = [i + 1 for i, a in enumerate(SEQ) if a == "T"]
    return out


def tryptic_coverage(synth, seq=SEQ, min_len=6, max_len=40):
    """Contre-argument a instruire : le regroupement des sites dans 20-24
    est-il un fait biologique, ou le simple reflet des peptides que la
    trypsine rend detectables ? Une proteine de 60 aa riche en K/R n'offre
    que quelques peptides de taille exploitable en LC-MS/MS."""
    cuts = [0] + [i + 1 for i, a in enumerate(seq)
                  if a in "KR" and i + 1 < len(seq)] + [len(seq)]
    peps = []
    for a, b in zip(cuts, cuts[1:]):
        p = seq[a:b]
        sty = [(a + i + 1, c) for i, c in enumerate(p) if c in "STY"]
        peps.append({
            "peptide": p, "start": a + 1, "stop": b, "longueur": len(p),
            "detectable_LCMS": min_len <= len(p) <= max_len,
            "sites_STY": [f"{c}{pos}" for pos, c in sty],
            "sites_detectes": sorted(s for s in synth if s in
                                     {f"{c}{pos}" for pos, c in sty}),
        })
    det = [p for p in peps if p["detectable_LCMS"]]
    sty_all = [f"{c}{i+1}" for i, c in enumerate(seq) if c in "STY"]
    sty_det = [s for p in det for s in p["sites_STY"]]
    return {
        "n_peptides_tryptiques": len(peps),
        "n_peptides_detectables": len(det),
        "peptides": peps,
        "n_STY_total": len(sty_all),
        "n_STY_dans_un_peptide_detectable": len(sty_det),
        "STY_hors_de_portee": [s for s in sty_all if s not in sty_det],
        "reserve": "Les sites ne peuvent etre appeles que dans les peptides "
                   "detectables ; l'absence de site ailleurs n'est pas une "
                   "absence de phosphorylation.",
    }


def main():
    res = {"question": "Quel est le site de phosphorylation exact de Rv0810c ?",
           "sequence": SEQ, "longueur_aa": len(SEQ),
           "source": "Malakar et al. 2023, mBio, PMID 37791794, "
                     "Data Sets S1-S8 (Europe PMC)"}
    sheets = {k: xlsx_sheets(DATA / v) for k, v in FILES.items()}

    # ---- 1. S2 : phosphopeptides du lysat, verifies sur la sequence ----------
    s2 = next(iter(sheets["S2_phosphopeptides_lysat"].values()))
    hdr2 = next(r for r in s2 if r and r[0].lower().startswith("rv number"))
    sites = []
    for r in s2:
        if not r or r[0] != "Rv0810c":
            continue
        d = dict(zip(hdr2, r))
        v = verify_peptide(d.get("Peptide sequence", ""),
                           int(d["Peptide start"]) if d.get("Peptide start") else None,
                           int(d["Peptide stop"]) if d.get("Peptide stop") else None,
                           d.get("Modification", ""))
        v["etiquette_publiee"] = d.get("Protein phospho site")
        v["jeu"] = "S2 (lysat, cette etude)"
        sites.append(v)
    res["S2_lysat"] = {"n_lignes_Rv0810c": len(sites), "sites": sites}

    # ---- 2. S8 : phosphosites du filtrat de culture --------------------------
    s8 = next(iter(sheets["S8_phosphosites_filtrat"].values()))
    hdr8 = next(r for r in s8 if r and r[0] == "Rv_number")
    naif, corrige = [], []
    for r in s8:
        d = dict(zip(hdr8, r))
        pep = d.get("Peptide_seq", "")
        start = int(d["peptide_start"]) if (d.get("peptide_start") or "").isdigit() else None
        stop = int(d["peptide_end"]) if (d.get("peptide_end") or "").isdigit() else None
        mod = (d.get("Modification") or "").replace("(Phospho)", "")
        label = d.get("Protein_phospho_site", "")
        # (a) lecture NAIVE : la colonne Rv_number de la meme ligne
        if d.get("Rv_number") == "Rv0810c":
            v = verify_peptide(pep, start, stop, mod)
            v["etiquette_publiee"] = label
            naif.append(v)
        # (b) lecture CORRIGEE : la colonne Protein_phospho_site, interne a la table
        if label.startswith("Rv0810c_"):
            v = verify_peptide(pep, start, stop, mod)
            v["etiquette_publiee"] = label
            v["jeu"] = "S8 (filtrat de culture)"
            corrige.append(v)
    res["S8_filtrat"] = {
        "lecture_naive_colonne_Rv_number": naif,
        "lecture_corrigee_colonne_Protein_phospho_site": corrige,
        "diagnostic_desalignement": {
            "n_lignes_avec_Rv_number": sum(1 for r in s8[2:] if r and r[0].startswith("Rv")),
            "n_lignes_de_peptides": sum(1 for r in s8[2:] if len(r) >= 5),
            "explication": "Colonne A = 38 proteines uniques ; colonnes B-G = 57 "
                           "phosphopeptides. Deux listes collees cote a cote, donc "
                           "desalignees des la 2e ligne de donnees.",
        },
    }

    # ---- 3. S3b : le meme site chez les autres etudes ------------------------
    inter = {}
    for name, rows in sheets["S3_comparaison_etudes"].items():
        # l'en-tete n'est pas forcement la 1re ligne (titre de l'onglet au-dessus)
        hdr = next((r for r in rows[:5] if any("Fortuin" in c for c in r)), None)
        if hdr is None:
            continue
        par_onglet = "site" if any("_" in r[0] for r in rows[3:8] if r) else "proteine"
        for r in rows:
            if r and r[0].startswith("Rv0810c"):
                inter[r[0] if par_onglet == "site" else f"{r[0]} ({name})"] = \
                    dict(zip(hdr, r))
    res["S3_comparaison_inter_etudes"] = inter

    # ---- 4. presence dans les autres jeux ------------------------------------
    presence = {}
    for key in ("S1_proteines_phosphorylees", "S4_secretome",
                "S5_secretome_multi_etudes", "S7_phosphoproteines_filtrat"):
        hits = []
        for name, rows in sheets[key].items():
            for r in rows:
                if any(c == "Rv0810c" or c == "I6XWB9" for c in r):
                    hits.append({"onglet": name, "ligne": r})
        presence[key] = hits
    res["presence_autres_jeux"] = presence

    # ---- 5. synthese : un site par residu ------------------------------------
    synth = {}
    for v in sites + corrige:
        if not v.get("site_calcule") or not v.get("residu_concorde"):
            continue
        s = v["site_calcule"]
        e = synth.setdefault(s, {"position": v["site_position"],
                                 "residu": v["site_aa"],
                                 "segment": v["segment"],
                                 "peptides": set(), "jeux": set()})
        e["peptides"].add(v["peptide"])
        e["jeux"].add(v["jeu"])
    # apport des autres etudes, depuis S3b
    for key, row in inter.items():
        m = re.match(r"Rv0810c_([STY])(\d+)$", key)
        if not m:
            continue
        s = f"{m.group(1)}{m.group(2)}"
        pos = int(m.group(2))
        e = synth.setdefault(s, {"position": pos, "residu": m.group(1),
                                 "segment": segment_of(pos),
                                 "peptides": set(), "jeux": set()})
        for etude in ("Fortuin", "Prisic", "Verma", "Carette", "This study_any two rep"):
            if row.get(etude) == "1":
                e["jeux"].add(etude)
        e["residu_concorde_sequence"] = (SEQ[pos - 1] == m.group(1)
                                         if pos <= len(SEQ) else False)
    for s, e in synth.items():
        e["peptides"] = sorted(e["peptides"])
        e["jeux"] = sorted(e["jeux"])
        e["n_jeux_independants"] = len(e["jeux"])
    res["synthese_sites"] = dict(sorted(synth.items(), key=lambda kv: kv[1]["position"]))

    # ---- 6. conservation des sites dans DUF3073 (PF11273) --------------------
    res["conservation_PF11273"] = conservation(synth)

    # ---- 7. detectabilite : quels S/T/Y sont seulement ATTEIGNABLES ? --------
    res["detectabilite_tryptique"] = tryptic_coverage(synth)

    best = max(synth.items(), key=lambda kv: kv[1]["n_jeux_independants"])
    res["site_principal"] = {
        "site": best[0], **best[1],
        "justification": "site retenu par le plus grand nombre de jeux independants",
    }

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(res, indent=2, ensure_ascii=False))

    # ---- affichage ----------------------------------------------------------
    print(f"Sequence Rv0810c ({len(SEQ)} aa)\n{SEQ}\n")
    print("S2 (lysat) — phosphopeptides attribues a Rv0810c :")
    for v in sites:
        print(f"  {v['peptide']:16s} {v['start_annonce']}-{v['stop_annonce']} "
              f"{v['modification']:4s} -> {v['site_calcule']} "
              f"| dans la sequence: {v['present_dans_la_sequence']} "
              f"| positions coherentes: {v['positions_coherentes']} "
              f"| residu concorde: {v['residu_concorde']} | {v.get('segment')}")
    print("\nS8 (filtrat) — lecture NAIVE de la colonne Rv_number :")
    for v in naif:
        print(f"  {v['peptide'][:32]:34s} {v['start_annonce']}-{v['stop_annonce']} "
              f"-> present dans la sequence: {v['present_dans_la_sequence']} "
              f"(etiquette publiee: {v['etiquette_publiee']})")
    print("S8 (filtrat) — lecture CORRIGEE (colonne Protein_phospho_site) :")
    for v in corrige:
        print(f"  {v['peptide']:16s} {v['start_annonce']}-{v['stop_annonce']} "
              f"{v['modification']:4s} -> {v['site_calcule']} "
              f"| dans la sequence: {v['present_dans_la_sequence']} "
              f"| tryptique: {v['tryptique']} | {v.get('segment')}")
    print("\nSites, tous jeux confondus :")
    for s, e in res["synthese_sites"].items():
        print(f"  {s:5s} pos {e['position']:2d}  {e['n_jeux_independants']} jeu(x) : "
              f"{', '.join(e['jeux'])}   [{e['segment']}]")
    cons = res["conservation_PF11273"]
    if "sites" in cons:
        print(f"\nConservation dans DUF3073 ({cons['n_sequences']} sequences, "
              f"{cons['n_positions_mappees']} positions) :")
        for site, v in cons["sites"].items():
            if "erreur" in v:
                print(f"  {site}: {v['erreur']}")
                continue
            print(f"  {site:5s} {v['residu_majoritaire']} majoritaire a "
                  f"{v['identite_max_pct']:5.1f}% | gaps {v['taux_gaps_pct']:5.2f}% | "
                  f"S/T/Y {v['phosphorylable_STY_pct']:5.1f}% | "
                  f"rang identite {v['rang_identite']}/{v['n_positions']} "
                  f"(percentile {v['percentile_identite']})")
        print(f"  positions sans aucun gap : {cons['positions_sans_aucun_gap']}")
        print(f"  Thr de la proteine : {cons['Thr_de_la_proteine']}")
    tc = res["detectabilite_tryptique"]
    print(f"\nDetectabilite tryptique : {tc['n_peptides_detectables']}/"
          f"{tc['n_peptides_tryptiques']} peptides exploitables en LC-MS/MS ; "
          f"{tc['n_STY_dans_un_peptide_detectable']}/{tc['n_STY_total']} S/T/Y "
          f"atteignables")
    for p in tc["peptides"]:
        if p["detectable_LCMS"]:
            print(f"  {p['peptide']:22s} {p['start']:2d}-{p['stop']:2d} "
                  f"S/T/Y {p['sites_STY']} -> detectes {p['sites_detectes']}")
    print(f"  S/T/Y hors de portee : {tc['STY_hors_de_portee']}")
    print(f"\nSITE PRINCIPAL : {res['site_principal']['site']} "
          f"({res['site_principal']['n_jeux_independants']} jeux)")
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
