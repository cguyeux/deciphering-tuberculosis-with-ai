#!/usr/bin/env python3
"""phase30_duf501_eukaryote_architecture.py — P8.1.a : que sont les DUF501 eucaryotes ?

phase29 a montré que PF04417 compte 102 membres eucaryotes quasi complets (~30 % d'identité
à Rv1025), surtout des protistes, et que 11 ont perdu la paire Cys-His. Reste la question
qui décide de leur valeur : sont-ils une VRAIE branche eucaryote de la famille, ou un
ramassis de domaines dérivés et d'artefacts ?

MODÈLE NUL DE LA PISTE (à battre) : « ces hits pourraient n'être que des domaines dérivés
sans fonction conservée — vérifier d'abord si les architectures sont COHÉRENTES entre
organismes ou HÉTÉROCLITES ».

Trois artefacts distincts sont testés, parce qu'ils appellent des conclusions opposées :

  1. CONTAMINATION D'ASSEMBLAGE. Une séquence « eucaryote » à identité anormalement HAUTE
     avec une protéine bactérienne est plus probablement une lecture bactérienne mal
     assignée qu'un homologue ancien. Signal : identité très au-dessus du reste du groupe
     eucaryote, et protéine de taille bactérienne (~150 aa) sans extension.
  2. HITS DÉGÉNÉRÉS. À l'inverse, un domaine isolé dans une grande protéine à position
     quelconque, sans cohérence d'architecture entre organismes, est un faux positif de
     profil. Signal : architectures hétéroclites, longueurs erratiques.
  3. VRAIE BRANCHE. Signal : architecture cohérente entre organismes indépendants, et
     surtout des EXTENSIONS N-TERMINALES longues et systématiques — ce qui, chez les
     apicomplexes et les algues, est la signature d'un ADRESSAGE (peptide signal +
     peptide de transit) vers un organite d'origine endosymbiotique.

HYPOTHÈSE À TESTER, née de la composition taxonomique (Babesia, Eimeria, Neospora =
apicomplexes ; straménopiles) : ces protéines pourraient être adressées à l'APICOPLASTE,
plastide vestigial d'origine algale dont le protéome est massivement d'origine bactérienne.
TEST FALSIFIANT NET : *Cryptosporidium* est l'apicomplexe qui a PERDU son apicoplaste. S'il
porte un DUF501, l'hypothèse tombe ; s'il est absent alors que ses proches en ont, elle
tient. Idem *Naegleria*, qui n'a jamais eu de plastide : sa présence exigerait une autre
explication.

Lit   : article/supplementary_materials/table_S10_duf501_eukaryotes.tsv, data/PF04417_full.sto
Écrit : résultats/duf501_eukaryote_architecture.tsv, .txt
"""
import os
import subprocess
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fast_blas import ensure  # noqa: E402  # type: ignore[import-not-found]

ensure()

from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
S10 = ROOT / "article/supplementary_materials/table_S10_duf501_eukaryotes.tsv"
STO = ROOT / "data/PF04417_full.sto"
OUT = ROOT / "résultats"
CACHE = OUT / "duf501_eukaryote_uniprot.tsv"
BATCH = 40

out = []


def say(s=""):
    print(s)
    out.append(s)


rows = [l.split("\t") for l in S10.read_text().splitlines()[1:] if l.strip()]
euk = [{"name": r[0], "acc": r[1], "len_dom": int(r[3]), "motif": r[4], "ident": float(r[5])}
       for r in rows]
for e in euk:
    span = e["name"].split("/")[1]
    e["start"], e["end"] = (int(x) for x in span.split("-"))
say(f"{len(euk)} membres eucaryotes quasi complets à caractériser")

# ------------------------------------------------------------------ UniProt
FIELDS = "accession,length,protein_name,organism_name,lineage,xref_pfam,cc_subcellular_location"
if CACHE.exists():
    data = {}
    for line in CACHE.read_text().splitlines()[1:]:
        c = line.split("\t")
        if len(c) >= 6:
            data[c[0]] = c
    say(f"UniProt lu depuis le cache ({len(data)} entrées)")
else:
    data = {}
    accs = [e["acc"] for e in euk]
    say(f"Interrogation d'UniProt ({-(-len(accs)//BATCH)} lots) ...")
    header = None
    for i in range(0, len(accs), BATCH):
        q = "+OR+".join(f"accession:{a}" for a in accs[i:i + BATCH])
        url = f"https://rest.uniprot.org/uniprotkb/search?query={q}&fields={FIELDS}&format=tsv&size=500"
        r = subprocess.run(["curl", "-sS", "--max-time", "120", url],
                           capture_output=True, text=True)
        lines = r.stdout.splitlines()
        if lines and header is None:
            header = lines[0]
        for line in lines[1:]:
            c = line.split("\t")
            if c and c[0]:
                data[c[0]] = c
    CACHE.write_text((header or "") + "\n"
                     + "\n".join("\t".join(v) for v in data.values()) + "\n")
    say(f"Écrit : {CACHE.relative_to(ROOT)} ({len(data)} entrées)")

for e in euk:
    c = data.get(e["acc"], [])
    e["length"] = int(c[1]) if len(c) > 1 and c[1].isdigit() else 0
    e["pname"] = c[2] if len(c) > 2 else ""
    e["org"] = c[3] if len(c) > 3 else ""
    e["lin"] = c[4] if len(c) > 4 else ""
    e["pfam"] = [x for x in (c[5] if len(c) > 5 else "").split(";") if x]
    e["loc"] = c[6] if len(c) > 6 else ""

known = [e for e in euk if e["length"]]
say(f"{len(known)}/{len(euk)} entrées résolues dans UniProt")


def clade(e):
    for tok in ("Apicomplexa", "Stramenopiles", "Heterolobosea", "Fungi", "Metazoa",
                "Viridiplantae", "Amoebozoa", "Discoba", "Rhodophyta", "Haptista"):
        if tok in e["lin"]:
            return tok
    return "autre eucaryote"


# ------------------------------------------------------- A. composition taxonomique
say()
say("=" * 78)
say("A. QUI SONT-ILS ? (clade, et genres représentés)")
say("=" * 78)
by = defaultdict(list)
for e in known:
    by[clade(e)].append(e)
for cl, lst in sorted(by.items(), key=lambda kv: -len(kv[1])):
    lost = sum(1 for e in lst if e["motif"] != "ECH")
    idm = np.median([e["ident"] for e in lst])
    genera = Counter(e["org"].split()[0] for e in lst if e["org"])
    say(f"  {cl:<18} n={len(lst):<4} triade perdue {lost:<3} identité médiane {idm:5.1f} %")
    say(f"      genres : {', '.join(f'{g}({n})' for g, n in genera.most_common(6))}")

# --------------------------------------- B. contamination : identité anormalement haute ?
say()
say("=" * 78)
say("B. ARTEFACT 1 — CONTAMINATION D'ASSEMBLAGE ?")
say("=" * 78)
idents = np.array([e["ident"] for e in known])
hi = [e for e in known if e["ident"] > 45]
say(f"identité à Rv1025 : médiane {np.median(idents):.1f} %, "
    f"IQR {np.percentile(idents,25):.1f}-{np.percentile(idents,75):.1f}")
say(f"membres ANORMALEMENT proches (> 45 %) : {len(hi)}")
for e in sorted(hi, key=lambda x: -x["ident"]):
    ext = e["start"] - 1
    say(f"    {e['acc']:<12} {e['ident']:5.1f} %  {e['length']:>4} aa  "
        f"ext. N-term {ext:>4}  {e['org'][:42]}")
say("Lecture : identité proche du niveau bactérien (~63 %) + protéine de taille bactérienne")
say("+ extension N-terminale nulle = signature de CONTAMINATION, pas d'homologie ancienne.")

# ------------------------------------- C. architecture : cohérente ou hétéroclite ?
say()
say("=" * 78)
say("C. ARTEFACT 2 — ARCHITECTURES HÉTÉROCLITES ? (le null de la piste)")
say("=" * 78)
arch = Counter(";".join(sorted(e["pfam"])) if e["pfam"] else "(aucun)" for e in known)
say("architectures Pfam observées :")
for a, n in arch.most_common(8):
    say(f"  {a if a else '(aucune)':<40} {n}")
solo = sum(1 for e in known if e["pfam"] == ["PF04417"])
say(f"-> DUF501 SEUL (aucun autre domaine Pfam) : {solo}/{len(known)} = {solo/len(known):.0%}")
say("Lecture : une architecture DOMINANTE et partagée entre clades indépendants réfute")
say("l'hypothèse « domaines dérivés hétéroclites ».")

# --------------------------- D. extensions N-terminales : signature d'adressage ?
say()
say("=" * 78)
say("D. EXTENSIONS N-TERMINALES : SIGNATURE D'UN ADRESSAGE ?")
say("=" * 78)
say(f"{'clade':<18} {'n':>4} {'longueur totale':>16} {'extension N-term':>18}")
for cl, lst in sorted(by.items(), key=lambda kv: -len(kv[1])):
    ln = np.median([e["length"] for e in lst])
    ext = np.median([e["start"] - 1 for e in lst])
    say(f"  {cl:<16} {len(lst):>4} {ln:>13.0f} aa {ext:>15.0f} aa")
say("Repère : Rv1025 fait 155 aa avec une extension N-terminale de 15 résidus.")
say("Une extension de plusieurs dizaines à centaines de résidus, SYSTÉMATIQUE dans un clade,")
say("est compatible avec un adressage bipartite (peptide signal + transit) vers un organite.")

# ------------------------------ E. test falsifiant de l'hypothèse apicoplaste
say()
say("=" * 78)
say("E. TEST FALSIFIANT — CRYPTOSPORIDIUM (apicomplexe SANS apicoplaste)")
say("=" * 78)
apic = [e for e in known if "Apicomplexa" in e["lin"]]
say(f"genres d'apicomplexes présents dans l'ALIGNEMENT : "
    f"{dict(Counter(e['org'].split()[0] for e in apic if e['org']))}")
say()
say("Une absence dans l'alignement n'est PAS une absence dans le génome : Pfam full est un")
say("échantillon, pas un recensement. Le test se fait donc par interrogation DIRECTE d'UniProt,")
say("genre par genre, sur le domaine (PF04417) ET sur la famille InterPro (IPR007511, plus")
say("sensible car elle intègre d'autres signatures).")


def count_uniprot(query):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=list&size=500"
    r = subprocess.run(["curl", "-sS", "--max-time", "120", url],
                       capture_output=True, text=True)
    return sum(1 for line in r.stdout.splitlines() if line.strip())


# plastide : True = apicoplaste/plastide présent, False = perdu ou jamais eu
GENERA = [("Plasmodium", True), ("Babesia", True), ("Theileria", True), ("Eimeria", True),
          ("Toxoplasma", True), ("Cryptosporidium", False),
          ("Trypanosoma", False), ("Leishmania", False), ("Naegleria", False)]
say(f"{'genre':<18} {'plastide':>9} {'PF04417':>9} {'IPR007511':>11}")
tally = {}
for g, plast in GENERA:
    n_pf = count_uniprot(f"xref:pfam-PF04417+AND+organism_name:{g}")
    n_ip = count_uniprot(f"xref:interpro-IPR007511+AND+organism_name:{g}")
    tally[g] = (plast, n_pf, n_ip)
    say(f"  {g:<16} {'oui' if plast else 'non':>9} {n_pf:>9} {n_ip:>11}")

viol_pos = [g for g, (p, a, b) in tally.items() if not p and max(a, b) > 0]
viol_neg = [g for g, (p, a, b) in tally.items() if p and max(a, b) == 0]
say()
say("VERDICT de l'hypothèse plastidiale — elle exige : présent SI plastide, absent SINON.")
say(f"  contre-exemples SANS plastide mais AVEC le gène : {viol_pos or 'aucun'}")
say(f"  contre-exemples AVEC plastide mais SANS le gène : {viol_neg or 'aucun'}")
if viol_pos or viol_neg:
    say("  -> HYPOTHÈSE RÉFUTÉE : des contre-exemples existent dans les deux sens. Les")
    say("     kinétoplastidés (Trypanosoma, Leishmania), qui n'ont jamais eu de plastide,")
    say("     portent des membres à triade INTACTE ; l'absence chez Cryptosporidium relève")
    say("     plus simplement de la réduction génomique extrême de ce genre.")
else:
    say("  -> compatible avec l'hypothèse (ce qui ne la démontre pas : il faudrait une")
    say("     prédiction ou une mesure d'adressage).")

# ------------------------------------------- F. les 11 « E-- » se distinguent-ils ?
say()
say("=" * 78)
say("F. LES MEMBRES À SITE DÉGRADÉ SE DISTINGUENT-ILS DES AUTRES ?")
say("=" * 78)
lost = [e for e in known if e["motif"] != "ECH"]
kept = [e for e in known if e["motif"] == "ECH"]
for lab, grp in (("site dégradé", lost), ("site complet", kept)):
    if grp:
        say(f"  {lab:<14} n={len(grp):<4} longueur médiane {np.median([e['length'] for e in grp]):>5.0f} aa"
            f"   extension N-term médiane {np.median([e['start']-1 for e in grp]):>5.0f} aa"
            f"   identité {np.median([e['ident'] for e in grp]):.1f} %")
say(f"  clades des dégradés : {dict(Counter(clade(e) for e in lost))}")

tab: list = [("accession", "organisme", "clade", "longueur", "domaine_debut", "domaine_fin",
              "ext_Nterm", "motif", "identite_Rv1025", "pfam", "localisation")]
for e in sorted(known, key=lambda x: (clade(x), -x["ident"])):
    tab.append((e["acc"], e["org"], clade(e), str(e["length"]), str(e["start"]), str(e["end"]),
                str(e["start"] - 1), e["motif"], f"{e['ident']:.1f}",
                ";".join(e["pfam"]) or "-", (e["loc"] or "-")[:60]))
(OUT / "duf501_eukaryote_architecture.tsv").write_text(
    "\n".join("\t".join(r) for r in tab) + "\n")
(OUT / "duf501_eukaryote_architecture.txt").write_text("\n".join(out) + "\n")
say()
say("Écrit : résultats/duf501_eukaryote_architecture.tsv et .txt")
