#!/usr/bin/env python3
"""
P4.2 — Le DENOMINATEUR de la retention reductive (BLOQUANT pour P4.1).

P4.1 a etabli que *Tropheryma whipplei*, l'actinobacterie au genome le plus
reduit, conserve DUF3073 (TWT_722). Tel quel, ce fait ne prouve rien : si les
genomes reduits conservent l'essentiel de leur coeur, la retention de Rv0810c
est *coherente avec* l'essentialite sans la demontrer.

Ce script produit le denominateur manquant, et corrige au passage un defaut de
raisonnement du plan initial : **un gene absent d'un genome reduit peut l'etre
par PERTE ou par simple DIVERGENCE a cette distance phylogenetique.** Sans
temoin non reduit a distance comparable, les deux sont indiscernables.

Dispositif apparie :
    reduit                       temoin NON reduit, meme clade
    ------------------------     ---------------------------------
    M. leprae      3,31 Mb   <-> M. abscessus   5,14 Mb  (genre Mycobacterium)
    M. lepromatosis 3,31 Mb  <-> M. abscessus   5,14 Mb
    T. whipplei    0,93 Mb   <-> M. luteus      2,54 Mb  (ordre Micrococcales)

La quantite qui mord est donc CONDITIONNELLE : parmi les proteines de H37Rv
encore detectables dans le temoin non reduit (donc assez anciennes et assez
conservees pour etre vues a cette distance), quelle fraction survit a la
reduction ? C'est le vrai denominateur.

Stratification obligatoire : essentialite (DeJesus 2017), LONGUEUR (une
proteine de 60 aa offre bien moins de signal tblastn qu'une de 400 : sans
appariement en longueur, la comparaison est un artefact de sensibilite), et
categorie fonctionnelle (conserved hypotheticals).

Controles de l'instrument, dans les deux sens :
  - positif : proteines ribosomiques (rpl/rps/rpm) — doivent etre retenues
    quasi partout. Si elles ne le sont pas, le seuil est trop severe.
  - negatif : familles PE/PPE, propres aux mycobacteries — doivent etre
    massivement absentes hors du genre. Si elles sont "retenues", le seuil
    est trop laxiste.

Sortie : résultats/p4_2_reductive_denominator.{json,tsv}
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from Bio import SeqIO

PROJ = Path(__file__).resolve().parent.parent
MTBC = PROJ.parent
PROTEOME = MTBC / "annotation_mtbc" / "résultats" / "phase2d_eggnog" / "proteome.faa"
ATLAS = MTBC / "annotation_mtbc" / "site" / "content" / "genes"
# NB : /tmp est un tmpfs de 32 Go monte en RAM sur cette machine. Les fichiers
# intermediaires de blast vont donc sur DISQUE, sinon un tblastn proteome-entier
# peut saturer la RAM et produire une sortie tronquee (constate le 2026-08-10 :
# ligne tblastn coupee en plein champ -> IndexError).
SCRATCH = PROJ / "experiments" / "2026-08-10_P4.2_work"
SCRATCH.mkdir(parents=True, exist_ok=True)
GENOME_DL = SCRATCH / "Twhipplei_Twist.fna"
OUTDIR = PROJ / "résultats"
OUTDIR.mkdir(exist_ok=True)

TARGET = "Rv0810c"
EVALUE = 1e-3          # recherche large, filtrage ensuite
E_KEEP = 1e-5
QCOV_PRESENT = 40.0    # homologue detectable
QCOV_INTACT = 80.0     # orthologue de pleine longueur
THREADS = 12

GENOMES = {
    "M_leprae": {
        "path": MTBC / "bdd/hors_mtbc/M_leprae/ref/genome.fna",
        "reduit": True,
        "clade": "Mycobacterium",
        "temoin_non_reduit": "M_abscessus",
        "taille_mb": 3.31,
    },
    "M_lepromatosis": {
        "path": MTBC / "bdd/hors_mtbc/M_lepromatosis/ref/genome.fna",
        "reduit": True,
        "clade": "Mycobacterium",
        "temoin_non_reduit": "M_abscessus",
        "taille_mb": 3.31,
    },
    "M_abscessus": {
        "path": MTBC / "bdd/hors_mycobacterium/Mabscessus_ATCC19977/ref/genome.fna",
        "reduit": False,
        "clade": "Mycobacterium",
        "temoin_non_reduit": None,
        "taille_mb": 5.14,
    },
    "T_whipplei": {
        "path": GENOME_DL,
        "reduit": True,
        "clade": "Micrococcales",
        "temoin_non_reduit": "M_luteus",
        "taille_mb": 0.93,
    },
    "M_luteus": {
        "path": MTBC / "bdd/hors_mycobacterium/Mluteus_NCTC2665/ref/genome.fna",
        "reduit": False,
        "clade": "Micrococcales",
        "temoin_non_reduit": None,
        "taille_mb": 2.54,
    },
}

FMT = "6 qseqid sseqid pident length evalue bitscore qcovs"


def run_tblastn(query: Path, genome: Path, tag: str) -> dict[str, dict]:
    """tblastn proteome -> genome. -seg no : sinon le filtrage basse complexite
    ampute la couverture des proteines a region de faible complexite (la queue
    acide de Rv0810c en est une), ce qui biaiserait la cible a la baisse."""
    db = SCRATCH / f"db_{tag}"
    db.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["makeblastdb", "-in", str(genome), "-dbtype", "nucl", "-out", str(db)],
        check=True,
        capture_output=True,
    )
    with tempfile.NamedTemporaryFile(
        "w+", suffix=".tsv", delete=False, dir=SCRATCH
    ) as tmp:
        out = Path(tmp.name)
    subprocess.run(
        [
            "tblastn", "-query", str(query), "-db", str(db),
            "-outfmt", FMT, "-evalue", str(EVALUE), "-seg", "no",
            "-num_threads", str(THREADS), "-max_target_seqs", "20",
            "-out", str(out),
        ],
        check=True,
        capture_output=True,
    )
    best: dict[str, dict] = {}
    n_bad = 0
    for line in out.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 7:  # ligne tronquee : ne JAMAIS l'avaler en silence
            n_bad += 1
            continue
        q, pid, ev, bs, qc = f[0], float(f[2]), float(f[4]), float(f[5]), float(f[6])
        cur = best.get(q)
        if cur is None or bs > cur["bitscore"]:
            best[q] = {"pident": pid, "evalue": ev, "bitscore": bs, "qcovs": qc}
        elif qc > cur["qcovs"]:
            cur["qcovs"] = qc
    if n_bad:
        raise RuntimeError(
            f"{tag}: {n_bad} lignes tblastn tronquees — sortie non fiable, "
            "verifier l'espace disque avant d'interpreter quoi que ce soit."
        )
    return best


def load_annotation() -> dict[str, dict]:
    ann: dict[str, dict] = {}
    for f in ATLAS.glob("Rv*.json"):
        try:
            d = json.loads(f.read_text())
        except Exception:
            continue
        ess = (d.get("essentiality") or {}).get("dejesus2017")
        ann[d["rv"]] = {
            "len_aa": d.get("len_aa"),
            "dejesus": ess,
            "essential": bool((d.get("essentiality") or {}).get("essential")),
            "funccat": (d.get("funccat") or {}).get("category"),
            "gene": d.get("gene") or "",
            "product": d.get("product_h37rv") or "",
        }
    return ann


def rate(mask: np.ndarray, ret: np.ndarray) -> dict:
    n = int(mask.sum())
    k = int(ret[mask].sum()) if n else 0
    if n == 0:
        return {"n": 0, "retenus": 0, "taux_pct": None}
    # IC binomial exact (Clopper-Pearson)
    from scipy.stats import beta

    lo = beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return {
        "n": n,
        "retenus": k,
        "taux_pct": round(100.0 * k / n, 1),
        "IC95_pct": [round(100 * lo, 1), round(100 * hi, 1)],
    }


def main() -> None:
    prots = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(PROTEOME, "fasta")}
    ann = load_annotation()
    ids = sorted(prots)
    lens = np.array([len(prots[i]) for i in ids])
    ess = np.array([ann.get(i, {}).get("dejesus") == "ES" for i in ids])
    hyp = np.array(
        [ann.get(i, {}).get("funccat") == "conserved hypotheticals" for i in ids]
    )
    gene = [ann.get(i, {}).get("gene", "") or "" for i in ids]
    prod = [ann.get(i, {}).get("product", "") or "" for i in ids]
    ribo = np.array(
        [g[:3] in ("rpl", "rps", "rpm") for g in gene]
    )
    pepp = np.array(
        [("PE-PGRS" in p) or ("PPE" in p) or p.startswith("PE family") for p in prod]
    )
    ti = ids.index(TARGET)

    hits: dict[str, dict] = {}
    for name, meta in GENOMES.items():
        print(f"[tblastn] {name} ({meta['taille_mb']} Mb) ...", flush=True)
        hits[name] = run_tblastn(PROTEOME, meta["path"], name)
        print(f"    {len(hits[name])} requetes avec >=1 hit", flush=True)

    present: dict[str, np.ndarray] = {}
    intact: dict[str, np.ndarray] = {}
    for name in GENOMES:
        h = hits[name]
        present[name] = np.array(
            [
                (i in h) and h[i]["evalue"] <= E_KEEP and h[i]["qcovs"] >= QCOV_PRESENT
                for i in ids
            ]
        )
        intact[name] = np.array(
            [
                (i in h) and h[i]["evalue"] <= E_KEEP and h[i]["qcovs"] >= QCOV_INTACT
                for i in ids
            ]
        )

    rep: dict = {
        "piste": "P4.2",
        "n_proteines_H37Rv": len(ids),
        "criteres": {
            "present": f"E<={E_KEEP} et qcov>={QCOV_PRESENT}%",
            "intact": f"E<={E_KEEP} et qcov>={QCOV_INTACT}%",
            "tblastn": "-seg no (evite d'amputer la couverture des regions de faible complexite)",
        },
        "genomes": {
            k: {kk: (str(vv) if kk == "path" else vv) for kk, vv in v.items()}
            for k, v in GENOMES.items()
        },
    }

    # ---------------- controles de l'instrument ----------------
    rep["controles_instrument"] = {}
    for name in GENOMES:
        rep["controles_instrument"][name] = {
            "positif_ribosomiques": rate(ribo, present[name]),
            "negatif_PE_PPE": rate(pepp, present[name]),
            "global": rate(np.ones(len(ids), bool), present[name]),
        }

    # ---------------- taux stratifies ----------------
    strata = {
        "tous": np.ones(len(ids), bool),
        "essentiels_ES": ess,
        "non_essentiels": ~ess,
        "conserved_hypotheticals": hyp,
        "longueur_50_70aa": (lens >= 50) & (lens <= 70),
        "longueur_<=100aa": lens <= 100,
        "longueur_>100aa": lens > 100,
        "ES_ET_50_70aa": ess & (lens >= 50) & (lens <= 70),
        "hypothetiques_50_70aa": hyp & (lens >= 50) & (lens <= 70),
        "hypothetiques_<=100aa": hyp & (lens <= 100),
    }
    rep["taux_par_strate"] = {}
    for name in GENOMES:
        rep["taux_par_strate"][name] = {
            s: {"present": rate(m, present[name]), "intact": rate(m, intact[name])}
            for s, m in strata.items()
        }

    # -------- LE denominateur : retention CONDITIONNELLE a la detectabilite --------
    rep["retention_conditionnelle"] = {
        "principe": (
            "Parmi les proteines de H37Rv encore detectables dans le temoin NON "
            "reduit du meme clade (donc assez anciennes/conservees pour etre vues "
            "a cette distance), quelle fraction survit a la reduction ? "
            "C'est le seul denominateur qui separe PERTE de DIVERGENCE."
        )
    }
    for name, meta in GENOMES.items():
        ctrl = meta["temoin_non_reduit"]
        if not ctrl:
            continue
        base = present[ctrl]
        sub = {
            "tous_detectables_dans_temoin": base,
            "essentiels_ES": base & ess,
            "non_essentiels": base & ~ess,
            "conserved_hypotheticals": base & hyp,
            "longueur_50_70aa": base & (lens >= 50) & (lens <= 70),
            "ES_ET_50_70aa": base & ess & (lens >= 50) & (lens <= 70),
            "hypothetiques_ET_ES": base & hyp & ess,
        }
        rep["retention_conditionnelle"][name] = {
            "temoin_non_reduit": ctrl,
            "strates": {s: rate(m, present[name]) for s, m in sub.items()},
        }

    # ---------------- la cible ----------------
    rep["cible_Rv0810c"] = {
        "len_aa": int(lens[ti]),
        "dejesus": ann[TARGET]["dejesus"],
        "funccat": ann[TARGET]["funccat"],
        "par_genome": {
            name: {
                "hit": hits[name].get(TARGET),
                "present": bool(present[name][ti]),
                "intact": bool(intact[name][ti]),
            }
            for name in GENOMES
        },
    }

    # table par gene pour reutilisation
    tsv = ["\t".join(["rv", "len_aa", "dejesus", "funccat"] + [
        f"{n}_{c}" for n in GENOMES for c in ("qcov", "pident", "present")
    ])]
    for j, i in enumerate(ids):
        row = [i, str(lens[j]), str(ann.get(i, {}).get("dejesus")), str(ann.get(i, {}).get("funccat"))]
        for n in GENOMES:
            h = hits[n].get(i)
            row += [
                f"{h['qcovs']:.0f}" if h else "0",
                f"{h['pident']:.1f}" if h else "0",
                "1" if present[n][j] else "0",
            ]
        tsv.append("\t".join(row))
    (OUTDIR / "p4_2_retention_par_gene.tsv").write_text("\n".join(tsv))

    out = OUTDIR / "p4_2_reductive_denominator.json"
    out.write_text(json.dumps(rep, indent=2, ensure_ascii=False))
    print(json.dumps(rep["controles_instrument"], indent=2, ensure_ascii=False))
    print(json.dumps(rep["retention_conditionnelle"], indent=2, ensure_ascii=False))
    print(json.dumps(rep["cible_Rv0810c"], indent=2, ensure_ascii=False))
    print(f"\n[ecrit] {out}")


if __name__ == "__main__":
    main()
