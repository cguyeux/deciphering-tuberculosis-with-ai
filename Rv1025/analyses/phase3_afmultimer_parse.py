#!/usr/bin/env python3
# pyright: reportArgumentType=false
# (stubs Bio.PDB : `get_structure` est typé Structure|None alors que l'échec lève ;
#  code vérifié à l'exécution — on marque le faux positif au lieu de tordre l'idiome.)
"""
P4.1 - Analyse les sorties de prédiction de complexes (AlphaFold3 Server OU Boltz-2).

Usage : python3 phase3_afmultimer_parse.py <dir_resultats>

Deux schémas de sortie sont reconnus automatiquement, par dossier de job :

  * AlphaFold3 Server : *_summary_confidences_*.json (iptm, ptm, ranking_score,
    chain_pair_iptm, chain_pair_pae_min) + *_model_*.cif.
  * Boltz-2 : confidence_<id>_model_<rank>.json (iptm, ptm, confidence_score,
    pair_chains_iptm, complex_plddt) + pae_<id>_model_<rank>.npz + <id>_model_<rank>.cif.
    ATTENTION : Boltz n'expose PAS `chain_pair_pae_min` ; le PAE minimum inter-chaînes
    est donc RECALCULÉ ici depuis la matrice brute du .npz, en affectant chaque token
    à sa chaîne d'après le mmCIF (convention : 1 token par résidu polymère, 1 token par
    atome de ligand). Si la somme des tokens ne correspond pas à la dimension de la
    matrice, on renvoie NA plutôt qu'un chiffre faux.

Sort un tableau : job | ipTM | pTM | ranking | PAE_min inter-chaînes | #contacts interface.
Critère de lecture : comparer au CONTRÔLE POSITIF du même système (le plafond d'ipTM est
spécifique au système, un seuil générique 0,6 rejette de vraies interfaces), lire surtout
le PAE inter-chaînes et la cohérence entre modèles.
"""
import glob, json, os, sys
from collections import defaultdict


def interface_contacts(cif_path, cutoff=5.0):
    """Nombre de paires de résidus inter-chaînes avec un atome à < cutoff Å."""
    try:
        from Bio.PDB.MMCIFParser import MMCIFParser
        import warnings
        warnings.simplefilter("ignore")
        st = MMCIFParser(QUIET=True).get_structure("m", cif_path)
    except Exception as e:  # noqa: BLE001
        return f"NA ({e})"
    chains = list(list(st)[0])
    if len(chains) < 2:
        return 0
    import numpy as np
    coords, labels = [], []
    for ch in chains:
        for res in ch:
            for at in res:
                coords.append(at.coord); labels.append((ch.id, res.id[1]))
    coords = np.array(coords)
    pairs = set()
    for i in range(len(coords)):
        ci, ri = labels[i]
        d = np.sqrt(((coords - coords[i]) ** 2).sum(1))
        for j in np.where(d < cutoff)[0]:
            cj, rj = labels[j]
            if ci != cj:
                pairs.add(tuple(sorted([(ci, ri), (cj, rj)])))
    return len(pairs)


def chain_tokens(cif_path):
    """Liste des chaînes token par token, convention AF3/Boltz.

    1 token par résidu polymère, 1 token par atome de ligand/hétéro-résidu.
    Renvoie None si la structure est illisible.
    """
    try:
        from Bio.PDB.MMCIFParser import MMCIFParser
        import warnings
        warnings.simplefilter("ignore")
        st = MMCIFParser(QUIET=True).get_structure("m", cif_path)
    except Exception:  # noqa: BLE001
        return None
    toks = []
    for ch in list(st)[0]:
        for res in ch:
            if res.id[0] == " ":                 # résidu polymère -> 1 token
                toks.append(ch.id)
            else:                                 # ligand/ion -> 1 token par atome
                toks.extend([ch.id] * len(list(res)))
    return toks


def pae_cross_min_from_npz(npz_path, cif_path):
    """PAE minimum inter-chaînes recalculé (Boltz). NA explicite si incohérent."""
    import numpy as np
    try:
        z = np.load(npz_path)
        pae = z[z.files[0]] if "pae" not in z.files else z["pae"]
    except Exception as e:  # noqa: BLE001
        return f"NA ({e})"
    pae = np.asarray(pae)
    if pae.ndim == 3:            # (1, N, N) éventuel
        pae = pae[0]
    toks = chain_tokens(cif_path) if cif_path else None
    if toks is None:
        return "NA (cif illisible)"
    if len(toks) != pae.shape[0]:
        # Ne PAS deviner : un mauvais découpage donnerait un chiffre faux.
        return f"NA (tokens {len(toks)} != PAE {pae.shape[0]})"
    toks = np.array(toks)
    uniq = list(dict.fromkeys(toks.tolist()))
    if len(uniq) < 2:
        return None
    best = None
    for a in range(len(uniq)):
        for b in range(a + 1, len(uniq)):
            ia = np.where(toks == uniq[a])[0]
            ib = np.where(toks == uniq[b])[0]
            block = pae[np.ix_(ia, ib)]
            m = float(block.min())
            best = m if best is None else min(best, m)
    return round(best, 2) if best is not None else None


def parse_af3(js):
    d = json.load(open(js))
    iptm, ptm, rank = d.get("iptm"), d.get("ptm"), d.get("ranking_score")
    pae_min = d.get("chain_pair_pae_min")
    cross = None
    if isinstance(pae_min, list) and len(pae_min) >= 2:
        cross = min(pae_min[0][1], pae_min[1][0])
    cif = js.replace("summary_confidences", "model").replace(".json", ".cif")
    if not os.path.exists(cif):
        cif = (glob.glob(f"{os.path.dirname(js)}/*model*.cif") or [None])[0]
    return iptm, ptm, rank, cross, cif


def parse_boltz(js):
    """confidence_<id>_model_<r>.json -> ipTM/pTM/score ; PAE recalculé depuis le .npz."""
    d = json.load(open(js))
    iptm, ptm = d.get("iptm"), d.get("ptm")
    rank = d.get("confidence_score")
    base = os.path.basename(js)[len("confidence_"):-len(".json")]   # <id>_model_<r>
    dirn = os.path.dirname(js)
    cif = os.path.join(dirn, f"{base}.cif")
    if not os.path.exists(cif):
        cif = (glob.glob(f"{dirn}/{base}.*cif") or glob.glob(f"{dirn}/*model*.cif") or [None])[0]
    npz = os.path.join(dirn, f"pae_{base}.npz")
    cross = pae_cross_min_from_npz(npz, cif) if os.path.exists(npz) else "NA (pas de .npz)"
    return iptm, ptm, rank, cross, cif


def main():
    if len(sys.argv) < 2:
        print("usage: phase3_afmultimer_parse.py <dir_resultats>"); sys.exit(1)
    root = sys.argv[1]
    by_job = defaultdict(list)          # job -> [(format, json)]
    for js in glob.glob(f"{root}/**/*.json", recursive=True):
        b = os.path.basename(js)
        if "summary_confidences" in b:
            by_job[os.path.basename(os.path.dirname(js))].append(("af3", js))
        elif b.startswith("confidence_"):
            by_job[os.path.basename(os.path.dirname(js))].append(("boltz", js))
    if not by_job:
        print(f"Aucune sortie de confiance (AF3 ou Boltz) sous {root}"); sys.exit(1)

    print(f"{'job':30s} {'src':>5s} {'ipTM':>6s} {'pTM':>6s} {'rank':>6s} {'PAExchain':>10s} {'iface':>7s}")
    rows = []
    for job, entries in sorted(by_job.items()):
        best, src = None, entries[0][0]
        for fmt, js in entries:
            iptm, ptm, rank, cross, cif = parse_af3(js) if fmt == "af3" else parse_boltz(js)
            key = rank if rank is not None else (iptm or 0)
            if best is None or key > best[0]:
                best, src = (key, iptm, ptm, rank, cross, cif), fmt
        if best is None:            # dossier sans sortie de confiance exploitable
            continue
        _, iptm, ptm, rank, cross, cif = best
        iface = interface_contacts(cif) if cif else "NA"
        print(f"{job:30s} {src:>5s} {iptm!s:>6} {ptm!s:>6} {rank!s:>6} {cross!s:>10} {iface!s:>7}")
        rows.append(dict(job=job, source=src, iptm=iptm, ptm=ptm, ranking=rank,
                         pae_cross=cross, iface_contacts=iface))
    out = f"{root}/afmultimer_summary.tsv"
    with open(out, "w") as fh:
        fh.write("job\tsource\tiptm\tptm\tranking\tpae_cross_min\tiface_contacts\n")
        for r in rows:
            fh.write(f"{r['job']}\t{r['source']}\t{r['iptm']}\t{r['ptm']}\t{r['ranking']}"
                     f"\t{r['pae_cross']}\t{r['iface_contacts']}\n")
    print(f"\nTSV : {out}")
    print("Lecture : calibrer sur le CONTRÔLE POSITIF du même système (pas de seuil ipTM générique),")
    print("privilégier le PAE inter-chaînes et la cohérence entre modèles ; comparer aux contrôles de spécificité.")
    print("NB : une valeur 'boltz' ne se compare pas directement à une valeur 'af3' (modèles différents).")


if __name__ == "__main__":
    main()
