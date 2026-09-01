#!/usr/bin/env python3
"""
Objet       : P7.5 -- DinB2/dinP (Rv3056) a des roles ETABLIS dans la
              mutagenese par SUBSTITUTION et par DECALAGE DE CADRE chez la
              mycobacterie (Ordonez et al., eLife 2023, 10.7554/eLife.83094,
              polymerase de translesion). P7.2 n'avait mesure QUE les
              substitutions (statistique C>A/(C>A+C>T)) sur ce meme gene,
              en EXPLORATOIRE (bilateral, p = 0,44, puissance 42 % a cap600 --
              silence de resolution, pas absence d'effet). La prediction sur
              les decalages de cadre n'avait jamais ete mesuree : ce script la
              teste, comme piste separee et pre-enregistree, PAS comme une
              nouvelle tentative sur le meme signal apres coup.

              CRIBLE /challenge du 2026-09-01 (avant tout code) : reutiliser
              les 5 EVENEMENTS deja qualifies par P7.2 (appariement cap600),
              PAS d'evenements nouveaux -- miner de nouveaux evenements parmi
              les 54 decalages catalogues par P7.1 romprait le principe que
              P7.1 avait lui-meme pose (les indels recidivent aux
              homopolymeres et fabriquent de la fausse convergence, c'est
              pour cela qu'ils sont exclus de Q2/Q3 la-bas). Le plan neuf
              porte donc sur une MESURE neuve (decalages de cadre prives),
              pas sur un CHOIX neuf d'evenements.

              REAPPARIEMENT INDEPENDANT, PAS REUTILISATION LITTERALE. Les
              identites des temoins de P7.2 n'ont pas ete serialisees sur
              disque (seuls les agregats le sont) ; ce script rejoue
              l'appariement porteur/temoin AVEC LES MEMES REGLES (memes
              repertoires de clade candidats, meme cap 600, meme tolerance
              20 % sur la profondeur en variants prives, fonction apparier()
              importee telle quelle de phase8_p72), mais avec son propre
              generateur aleatoire. Le jeu de temoins obtenu n'est donc pas
              bit-a-bit identique a celui de P7.2, seule la METHODE l'est.
              A DECLARER, pas a masquer.

              PREDICTION DIRIGEE, ECRITE AVANT LA MESURE. DinB2 CONTRIBUE a la
              mutagenese (substitutions ET decalages de cadre) sous stress
              genotoxique in vivo (eLife 2023). Sa PERTE doit donc REDUIRE la
              part des decalages de cadre parmi les variants prives, pas
              l'augmenter -- a l'oppose des genes de reparation dont la perte
              AUGMENTE le taux de mutation. Test UNILATERAL a gauche.

              STATISTIQUE PRIMAIRE : part des decalages de cadre parmi les
              variants prives totaux, indels_prives / (indels_prives +
              substitutions_prives), par cote (porteurs, temoins).
              Substitutions au denominateur, comme mesure de
              l'opportunite/profondeur de branche, pour ne pas dependre d'un
              modele d'opportunite des indels (tres heterogene, domine par le
              contenu en homopolymeres) -- meme principe que la statistique de
              P7.2.

              Q2 CALCULEE ET IMPRIMEE AVANT LE TEST, comme P7.2. Effet par
              defaut : cote x0.6 sur la part de decalages parmi les porteurs
              (reduction attendue, cf. prediction dirigee). Si la puissance
              est sous 50 %, un silence sera declare SILENCE DE RESOLUTION,
              pas absence d'effet -- ecrit ici, avant le calcul.

Entrees     : résultats/phase8_p72_evenements_qualifies_cap600.tsv (filtre :
              locus == Rv3056, qualifie == True -- les 5 evenements)
              résultats/phase8_p71_evenements.tsv, data/p71_variants_3r.tsv.gz
              bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt
              data/MTBC0/ancestral_on_H37Rv.bin, data/mask_h37rv_positions.npy
Sorties     : résultats/phase8_p75_evenements.tsv
              résultats/phase8_p75_verdict.tsv
Reutilisable: partiellement -- la mecanique d'appariement est celle de
              phase8_p72 (importee), la mesure des decalages prives est neuve
              et vaut pour toute question opposant substitutions et indels
              prives sur un meme pool apparie
Projet      : GC_par_lignee
Date        : 2026-09-01
"""
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase2_polarisation_mtbc0 import build_ancestral  # noqa: E402
from phase5_p51_panel_recursif import masked_positions  # noqa: E402
from phase8_p72_temoins_apparies import (  # noqa: E402
    GENES, BDD, ROOT, apparier, subs_cache, prives,
)

_CACHE_INDELS = {}
CAP = 600


def read_indels(spdi_path):
    """Variants INDELS d'une souche (longueur ref != longueur alt), tous
    confondus (insertion/deletion, sans filtre de modulo-3 : la mesure porte
    sur le taux d'accident de replication genome-entier, pas sur la
    consequence proteique d'un decalage particulier)."""
    out = set()
    for line in spdi_path.read_text().splitlines():
        p = line.strip().split(":")
        if len(p) != 4:
            continue
        ref, alt = p[2].upper(), p[3].upper()
        if len(ref) == len(alt) or not ref.isalpha() or not alt.isalpha():
            continue
        out.add((int(p[1]), ref, alt))
    return out


def indels_cache(f):
    k = str(f)
    if k not in _CACHE_INDELS:
        _CACHE_INDELS[k] = read_indels(f)
    return _CACHE_INDELS[k]


def prives_count(subsets):
    """Compte, pour chaque souche, le nombre de variants PRIVES du pool
    (presents chez une seule souche) -- meme logique de privaute que
    prives()/phase8_p72, generalisee a n'importe quel type de variant."""
    sup = defaultdict(list)
    for j, s in enumerate(subsets):
        for v in s:
            sup[v].append(j)
    n = np.zeros(len(subsets), int)
    for v, js in sup.items():
        if len(js) == 1:
            n[js[0]] += 1
    return n


def assemble_pool(porteurs, clades, cap, seed):
    """Reproduit l'algorithme de pool_local() (phase8_p72) mais conserve le
    couple (clade, sra) de chaque souche retenue, necessaire pour relire ses
    indels -- pool_local() ne renvoie que le nom de souche, pas son chemin."""
    rng = np.random.default_rng(seed)
    cand = []
    for c in sorted(clades):
        d = BDD / c
        if not d.is_dir():
            continue
        s = [x.name for x in sorted(d.iterdir())
             if x.is_dir() and (x / "NC_000962.3" / "spdi.txt").exists()]
        cand.append((c, s))
    sra_p = {s for _, s in porteurs}
    temoins = []
    for c, s in cand:
        rng.shuffle(s)
    i = 0
    while len(temoins) < cap and any(i < len(s) for _, s in cand):
        for c, s in cand:
            if i < len(s) and s[i] not in sra_p:
                temoins.append((c, s[i]))
                if len(temoins) >= cap:
                    break
        i += 1
    membres, subsets, est_porteur = [], [], []
    for c, sra in list(porteurs) + temoins:
        f = BDD / c / sra / "NC_000962.3" / "spdi.txt"
        if not f.exists():
            continue
        v = subs_cache(f, masked_positions())
        if not v:
            continue
        membres.append((c, sra))
        subsets.append(v)
        est_porteur.append(sra in sra_p)
    return membres, subsets, np.array(est_porteur, bool)


def main():
    rng = np.random.default_rng(0)
    anc, masked = build_ancestral(verbose=False), masked_positions()

    qual = pd.read_csv(ROOT / "résultats" / "phase8_p72_evenements_qualifies_cap600.tsv",
                       sep="\t")
    qual = qual[(qual.locus == "Rv3056") & (qual.qualifie == True)]  # noqa: E712
    print(f"Evenements dinP repris de P7.2 (cap{CAP}, qualifiants) : {len(qual)}")
    assert len(qual) == 5, "attendu : les 5 evenements deja qualifies par P7.2"

    var = pd.read_csv(ROOT / "data" / "p71_variants_3r.tsv.gz", sep="\t")
    ev71 = pd.read_csv(ROOT / "résultats" / "phase8_p71_evenements.tsv", sep="\t")
    ev71 = ev71[(ev71.locus == "Rv3056") & (ev71.type == "non-sens")]

    lignes = []
    for r in qual.itertuples():
        e = ev71[ev71.pos == r.pos].iloc[0]
        sel = var[(var.pos == e.pos) & (var.ref == e.ref) & (var.alt == e.alt)]
        porteurs = [(x.clade, x.sra) for x in sel.itertuples()]
        clades = {c for c, _ in porteurs}
        membres, subsets, est_p = assemble_pool(porteurs, clades, CAP, seed=0)
        if est_p.sum() < 2 or (~est_p).sum() < 10:
            print(f"  pos {r.pos} : reappariement echoue (trop peu de porteurs/"
                  f"temoins lus), evenement ecarte")
            continue
        sp, tot_sub = prives(subsets, anc)
        idx_t, ecart = apparier(tot_sub, est_p, rng)
        if idx_t is None:
            print(f"  pos {r.pos} : reappariement echoue (profondeur non "
                  f"appariable, ecart = {ecart:.3f}), evenement ecarte")
            continue

        subsets_idl = [indels_cache(BDD / c / sra / "NC_000962.3" / "spdi.txt")
                       for c, sra in membres]
        n_idl = prives_count(subsets_idl)

        ip = np.flatnonzero(est_p)
        idl_p, idl_t = int(n_idl[ip].sum()), int(n_idl[idx_t].sum())
        sub_p, sub_t = int(tot_sub[ip].sum()), int(tot_sub[idx_t].sum())
        part_p = idl_p / (idl_p + sub_p) if (idl_p + sub_p) else np.nan
        part_t = idl_t / (idl_t + sub_t) if (idl_t + sub_t) else np.nan
        lignes.append(dict(
            pos=r.pos, n_porteurs_lus=int(est_p.sum()),
            n_temoins_lus=len(idx_t), ecart_profondeur=ecart,
            idl_porteurs=idl_p, sub_porteurs=sub_p, part_porteurs=part_p,
            idl_temoins=idl_t, sub_temoins=sub_t, part_temoins=part_t,
            delta=part_p - part_t))

    df = pd.DataFrame(lignes)
    df.to_csv(ROOT / "résultats" / "phase8_p75_evenements.tsv", sep="\t", index=False)
    print()
    print(df.round(4).to_string(index=False))

    verdict = dict(n=len(df))
    if len(df) < 4:
        print(f"\n{len(df)} evenements appariables, trop peu pour un test. "
              f"Aucun verdict.")
        verdict["verdict"] = "non testable"
    else:
        d_prof = df.ecart_profondeur.to_numpy()  # deja un ecart relatif, pas un delta signe
        rho, p_q1 = stats.spearmanr(np.abs(d_prof), np.abs(df.delta.to_numpy()))
        print(f"\nQ1 appariement : rho(|ecart profondeur|, |delta statistique|) = "
              f"{rho:+.3f} (p = {p_q1:.3f}) -> "
              f"{'PASSE' if p_q1 > 0.05 else 'ECHOUE, test degrade'}")

        effet = 0.6  # cote x0.6 : reduction attendue chez les porteurs
        sim_n = 4000
        sim = np.zeros(sim_n)
        for b in range(sim_n):
            dd = []
            for row in df.itertuples():
                base = (row.part_porteurs * (row.idl_porteurs + row.sub_porteurs) +
                        row.part_temoins * (row.idl_temoins + row.sub_temoins)) / \
                       max((row.idl_porteurs + row.sub_porteurs) +
                           (row.idl_temoins + row.sub_temoins), 1)
                o = base / max(1 - base, 1e-9)
                pp = np.clip(o * effet / (1 + o * effet), 0, 1)
                n_p = max(row.idl_porteurs + row.sub_porteurs, 1)
                n_t = max(row.idl_temoins + row.sub_temoins, 1)
                a = rng.binomial(n_p, pp) / n_p
                c = rng.binomial(n_t, base) / n_t
                dd.append(a - c)
            dd = np.array(dd)
            try:
                sim[b] = stats.wilcoxon(dd, alternative="less").pvalue
            except ValueError:
                sim[b] = 1.0
        puissance = float((sim <= 0.05).mean())
        print(f"Q2 PUISSANCE, calculee avant le test : sous une cote de x{effet:g} "
              f"sur la part de decalages parmi les porteurs (reduction attendue), "
              f"ce plan verrait l'effet dans {100*puissance:.1f} % des cas")
        if puissance < 0.50:
            print("   -> sous 50 % : un silence sera declare SILENCE DE "
                  "RESOLUTION, pas absence d'effet")

        try:
            st = stats.wilcoxon(df.delta.to_numpy(), alternative="less")
            w, p = float(st.statistic), float(st.pvalue)
        except ValueError:
            w, p = np.nan, 1.0
        med = float(np.median(df.delta))
        signes = int((df.delta < 0).sum())
        print(f"\nTEST UNIQUE (unilateral a gauche, prediction dirigee : reduction "
              f"des decalages parmi les porteurs) :")
        print(f"  mediane des differences appariees = {med:+.4f}, "
              f"{signes}/{len(df)} evenements de signe negatif (sens predit)")
        print(f"  Wilcoxon W = {w:.1f}, p = {p:.4g}")
        v = ("EFFET DETECTE" if p <= 0.05 else
             "silence de resolution" if puissance < 0.50 else
             "PAS D'EFFET, a resolution suffisante")
        print(f"  verdict : {v}")
        verdict.update(rho_q1=rho, p_q1=p_q1, puissance=puissance,
                       mediane_delta=med, signes_sens_predit=signes,
                       W=w, p=p, verdict=v)

    pd.DataFrame([verdict]).to_csv(ROOT / "résultats" / "phase8_p75_verdict.tsv",
                                   sep="\t", index=False)


if __name__ == "__main__":
    main()
