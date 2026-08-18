#!/usr/bin/env python3
"""
P3.1 — HHsearch profil-profil sur Rv0810c / DUF3073.

Pourquoi ce chemin plutot que HHpred web. La KB documente deux echecs
reproduits (2026-06-06, 2026-07-04) de l'API REST du MPI Toolkit : la
soumission passe, le MSA se construit, puis hhsearch plante pour les jobs
anonymes. Seul le chemin navigateur aboutit — indisponible ici (extension
Chrome deconnectee). On fait donc tourner hhsearch EN LOCAL, ce qui est
en outre plus reproductible pour un manuscrit qu'une soumission web.

Difference methodologique a assumer, et elle joue en notre faveur : HHpred
construit le profil requete par HHblits contre UniRef30 depuis la sequence
seule. Ici le profil vient directement de l'ALIGNEMENT PFAM CURE de la famille
(PF11273, 1930 sequences, Rv0810c y figure sous I6XWB9_MYCTU/2-60). Le profil
est donc au moins aussi profond, et il n'est pas exposé au risque de derive
d'une recherche iterative.

Conversion Stockholm -> A3M : dans un alignement Pfam, les colonnes match sont
en MAJUSCULES ou '-', les colonnes d'insertion en minuscules ou '.'. C'est du
A2M ; retirer les '.' donne du A3M valide, sans passer par reformat.pl.

Controle positif OBLIGATOIRE de la chaine : la recherche contre Pfam DOIT
retrouver PF11273 (DUF3073) lui-meme a une probabilite proche de 100 %. Si ce
n'est pas le cas, l'instrument est casse et aucun negatif n'est interpretable.

Garde-fou de lecture (KB, 2026-07-04) : ne JAMAIS conclure sur un hit sans SA
probabilite et son E-value — le hit n°1 peut n'etre que le moins mauvais des
non-hits. Et nommer le REPLI, pas la sous-famille.

Sortie : résultats/p3_1_hhsearch/*.hhr + p3_1_hhsearch.json
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
RES = PROJ.parent / "investigate_phylo" / "resources"
HH = RES / "hhsuite_bin"
DB = RES / "hhsuite"
STO = PROJ / "data" / "PF11273_full.sto"
OUT = PROJ / "résultats" / "p3_1_hhsearch"
OUT.mkdir(parents=True, exist_ok=True)

QUERY_ID = "I6XWB9_MYCTU"       # Rv0810c dans l'alignement Pfam
# Region ORDONNEE etablie en P2.1 : residus 1-33 (pLDDT 91,9 sur 1-19, 85,6 sur
# 20-33), contre 62,0 sur la queue 34-60. La queue desordonnee de 27 residus
# dilue mecaniquement le signal profil-profil : on interroge donc aussi le
# module seul. L'alignement Pfam demarre au residu 2 de la proteine.
MODULE_RES = (2, 33)             # en coordonnees proteine, borne basse = 2


def sto_to_a3m(sto: Path, first: str) -> tuple[list[tuple[str, str]], str]:
    """Stockholm Pfam -> A3M. Place `first` en tete (il definit les etats match)."""
    seqs: dict[str, list[str]] = {}
    order: list[str] = []
    for line in sto.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        name, block = parts
        if name not in seqs:
            seqs[name] = []
            order.append(name)
        seqs[name].append(block.strip())
    a3m = {n: "".join(v).replace(".", "") for n, v in seqs.items()}
    hits = [n for n in order if n.startswith(first)]
    if not hits:
        raise SystemExit(f"{first} absent de {sto}")
    q = hits[0]
    ordered = [(q, a3m[q])] + [(n, a3m[n]) for n in order if n != q]
    return ordered, a3m[q]


def write_a3m(pairs: list[tuple[str, str]], path: Path) -> None:
    path.write_text("".join(f">{n}\n{s}\n" for n, s in pairs))


def slice_a3m(seq: str, col_lo: int, col_hi: int) -> str:
    """Extrait les colonnes MATCH [col_lo, col_hi] (1-based) d'une ligne A3M.

    Un A3M n'est PAS colonne-alignable : les insertions (minuscules) decalent
    les positions d'une sequence a l'autre. Un slice par indice brut melangerait
    des colonnes differentes selon les sequences. On compte donc les etats match
    (majuscule ou '-') et on conserve les insertions situees a l'interieur.
    """
    out: list[str] = []
    m = 0
    for ch in seq:
        is_match = ch.isupper() or ch == "-"
        if is_match:
            m += 1
            if col_lo <= m <= col_hi:
                out.append(ch)
        elif col_lo <= m < col_hi:  # insertion interne au segment
            out.append(ch)
    return "".join(out)


def parse_hhr(path: Path, top: int = 25) -> list[dict]:
    """Extrait le tableau de hits d'un .hhr (colonnes fixes de hhsearch)."""
    lines = path.read_text().splitlines()
    try:
        i = next(k for k, l in enumerate(lines) if l.startswith(" No Hit"))
    except StopIteration:
        return []
    out = []
    for l in lines[i + 1 :]:
        if not l.strip():
            break
        m = re.match(
            r"\s*(\d+)\s+(.{30})\s+([\d.]+)\s+([\dEe.+-]+)\s+([\dEe.+-]+)\s+"
            r"([-\d.]+)\s+([-\d.]+)\s+(\d+)\s+(\S+)\s+(\S+)",
            l,
        )
        if not m:
            continue
        out.append(
            {
                "rang": int(m.group(1)),
                "hit": m.group(2).strip(),
                "Prob": float(m.group(3)),
                "E_value": float(m.group(4)),
                "P_value": float(m.group(5)),
                "Score": float(m.group(6)),
                "SS": float(m.group(7)),
                "Cols": int(m.group(8)),
                "Query_HMM": m.group(9),
                "Template_HMM": m.group(10),
            }
        )
        if len(out) >= top:
            break
    return out


def run(cmd: list[str]) -> None:
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd[:3])} -> {r.returncode}\n{r.stderr[-3000:]}")


def main() -> None:
    pairs, qseq = sto_to_a3m(STO, QUERY_ID)
    rep: dict = {
        "piste": "P3.1",
        "profil_requete": {
            "source": "alignement Pfam full PF11273 (release courante InterPro)",
            "n_sequences": len(pairs),
            "requete_en_tete": pairs[0][0],
            "sequence_requete_a3m": qseq,
        },
        "garde_fou_lecture": (
            "Prob > 95 % = hit confiant ; 50-95 % = a signaler ; sous 50 %, la "
            "CONVERGENCE de plusieurs hits vers un meme repli compte autant qu'un "
            "hit isole. Ne jamais lire le rang 1 sans sa Prob et son E-value."
        ),
        "recherches": {},
    }

    # --- profil complet -----------------------------------------------------
    # NB : sur une entree A3M, hhmake NE DOIT PAS recevoir -M — les etats match
    # y sont deja encodes par la casse. Passer -M first lui fait lire le fichier
    # comme un FASTA aligne et il exige alors des longueurs egales, ce qu'un A3M
    # n'a jamais (les insertions decalent chaque ligne).
    a3m = OUT / "duf3073_full.a3m"
    write_a3m(pairs, a3m)
    hhm = OUT / "duf3073_full.hhm"
    run([str(HH / "hhmake"), "-i", str(a3m), "-o", str(hhm), "-v", "0"])

    # --- profil de la REGION ORDONNEE SEULE --------------------------------
    col_lo = MODULE_RES[0] - 1  # l'alignement demarre au residu 2 -> colonne 1
    col_hi = MODULE_RES[1] - 1
    keep = [(n, slice_a3m(s, col_lo, col_hi)) for n, s in pairs]
    keep = [(n, s) for n, s in keep if s.replace("-", "").strip()]
    a3m_m = OUT / "duf3073_module.a3m"
    write_a3m(keep, a3m_m)
    hhm_m = OUT / "duf3073_module.hhm"
    run([str(HH / "hhmake"), "-i", str(a3m_m), "-o", str(hhm_m), "-v", "0"])
    rep["profil_module"] = {
        "residus_proteine": list(MODULE_RES),
        "colonnes_match_alignement": [col_lo, col_hi],
        "n_sequences": len(keep),
        "sequence_module_requete": keep[0][1],
    }

    dbs = {"pfam": DB / "pfam"}
    if (DB / "scop70_1.75_hhm.ffdata").exists():
        dbs["scop70"] = DB / "scop70_1.75"

    for prof_name, prof in (("proteine_entiere", hhm), ("module_2_33", hhm_m)):
        for db_name, dbp in dbs.items():
            hhr = OUT / f"{prof_name}__{db_name}.hhr"
            run(
                [
                    str(HH / "hhsearch"), "-i", str(prof), "-d", str(dbp),
                    "-o", str(hhr), "-cpu", "8", "-v", "0",
                ]
            )
            rep["recherches"][f"{prof_name}__{db_name}"] = parse_hhr(hhr)

    # --- controle positif ---------------------------------------------------
    top_pfam = rep["recherches"].get("proteine_entiere__pfam", [])
    self_hit = next(
        (h for h in top_pfam if "DUF3073" in h["hit"] or "PF11273" in h["hit"]), None
    )
    rep["controle_positif"] = {
        "attendu": "la recherche Pfam doit retrouver DUF3073 (PF11273) a Prob ~100 %",
        "auto_hit_trouve": self_hit,
        "instrument_valide": bool(self_hit and self_hit["Prob"] >= 95.0),
    }

    (PROJ / "résultats" / "p3_1_hhsearch.json").write_text(
        json.dumps(rep, indent=2, ensure_ascii=False)
    )
    print(json.dumps(rep["controle_positif"], indent=2, ensure_ascii=False))
    for k, v in rep["recherches"].items():
        print(f"\n===== {k} =====")
        for h in v[:12]:
            print(
                f"  {h['rang']:>2}. Prob={h['Prob']:>6.1f}  E={h['E_value']:<10.2g} "
                f"Cols={h['Cols']:>3}  {h['hit'][:70]}"
            )


if __name__ == "__main__":
    main()
