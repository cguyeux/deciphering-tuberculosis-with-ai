#!/usr/bin/env python3
"""
Objet       : P7.2 -- comparer chaque porteur d'une inactivation convergente de
              gene 3R a des TEMOINS APPARIES phylogenetiquement et en
              PROFONDEUR, et non a la moyenne du complexe. P7.1 a livre le
              dispositif que P7 n'avait pas : mutT3 (Rv0413) et dinP (Rv3056)
              sont frappes par 14 et 15 evenements non-sens independants, a
              ensembles de porteurs deux a deux non emboites, repartis sur cinq
              et huit lignees. Chaque evenement est donc un REPLICAT, et
              l'allele se separe du fond genetique -- ce que les alleles nommes
              par la litterature ne permettaient pas (A33).

              L'UNITE DE REPLICATION EST L'EVENEMENT, PAS LA SOUCHE. Les
              porteurs d'un meme evenement descendent de la branche ou il est
              apparu : ils sont clonaux et ne comptent que pour une observation.
              Toutes les mesures sont donc agregees a l'interieur d'un
              evenement, et le test porte sur les differences appariees ENTRE
              evenements.

              CE QUI EST MESURE, ET POURQUOI C'EST POST-INACTIVATION. Les
              variants PRIVES d'une souche dans son pool local sont ceux de sa
              branche terminale, donc en aval de la branche ou l'inactivation
              est apparue. Le spectre de ces variants prives est donc accumule
              APRES la perte de fonction chez les porteurs, et sur une duree
              comparable chez les temoins puisque l'appariement fixe la
              profondeur.

              LA PREDICTION, DIRIGEE ET ECRITE AVANT LA MESURE. mutT3 est une
              8-oxo-dGTP diphosphatase : sa perte laisse le 8-oxo-dGTP entrer
              dans l'ADN et produit des TRANSVERSIONS G:C -> T:A, c'est-a-dire
              le canal C>A. Prediction : chez les porteurs, la part de C>A parmi
              les pertes de paires G:C AUGMENTE. dinP est une polymerase de
              translesion et ne porte AUCUNE prediction dirigee ; elle est
              testee en bilateral et declaree exploratoire.
              ET IL FAUT LE DIRE MAINTENANT : cette prediction porte sur un
              canal de PERTE, alors que A29 et A31 designent les GAINS comme la
              composante qui distingue les lignees. Un succes de P7.2 ne
              confirmerait donc PAS le pointeur de P8, et un echec ne
              l'infirmerait pas non plus. Ecrit ici pour qu'aucun des deux
              resultats ne soit relu apres coup comme la confirmation de ce
              qu'on preferait.

              STATISTIQUE PRIMAIRE, CHOISIE POUR ETRE EXEMPTE D'OPPORTUNITE :
              C>A / (C>A + C>T), la part de la transversion parmi les pertes de
              paires G:C. Numerateur et denominateur portent sur les MEMES sites
              G:C, donc l'opportunite se simplifie exactement et aucun modele de
              composition n'entre dans la mesure.

              CRITERES PRE-ENREGISTRES, FIXES AVANT TOUTE EXECUTION
                D1 QUALIFICATION D'UN EVENEMENT. >= 2 porteurs lisibles, >= 10
                   temoins candidats non porteurs dans les memes repertoires de
                   clade, et un jeu de temoins apparie tel que l'ecart relatif
                   des profondeurs moyennes (variants prives par souche) reste
                   <= 20 %. Un evenement non appariable est exclu, et
                   l'exclusion est une propriete du plan.
                D2 TRONCATURE UTILE. Un codon stop dans les 10 % terminaux du
                   CDS ne detruit pas forcement la proteine ; l'analyse primaire
                   ne retient que les stops situes dans les 90 % premiers, regle
                   standard d'appel de perte de fonction. Les autres sont
                   rapportes a part.
                Q1 RESOLUTION DE L'APPARIEMENT. Apres appariement, la difference
                   de statistique primaire entre porteurs et temoins ne doit
                   plus correler a la difference de profondeur residuelle
                   (Spearman, p > 0,05). Sinon l'appariement n'a pas fait son
                   office et le test n'est pas rendu.
                Q2 PUISSANCE, CALCULEE ET IMPRIMEE AVANT LE TEST. Sous un effet
                   de +50 % en cote sur la part de C>A parmi les pertes, avec
                   les effectifs de substitutions reellement obtenus, quelle est
                   la probabilite que ce plan le voie ? Si elle est inferieure a
                   50 %, le silence sera declare SILENCE DE RESOLUTION et non
                   absence d'effet, quel que soit le resultat.
                TEST UNIQUE : test des rangs signes de Wilcoxon sur les
                   differences appariees porteurs - temoins de la part de C>A
                   parmi les pertes, sur les evenements qualifiants de mutT3,
                   UNILATERAL a droite puisque la prediction est dirigee. Rendu
                   quel que soit son resultat. dinP est bilateral et
                   exploratoire, et ne peut PAS a lui seul soutenir une
                   conclusion mecaniste.
              Aucun autre gene, aucune autre statistique, aucune autre
              definition de temoin ne seront essayes.

Entrees     : résultats/phase8_p71_evenements.tsv, data/p71_variants_3r.tsv.gz
              bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt
              data/MTBC0/ancestral_on_H37Rv.bin, data/mask_h37rv_positions.npy
Sorties     : résultats/phase8_p72_evenements_qualifies.tsv (D1, D2)
              résultats/phase8_p72_appariement.tsv          (Q1, profondeurs)
              résultats/phase8_p72_spectres.tsv             (spectre des deux cotes)
              résultats/phase8_p72_verdict.tsv              (Q2 + test unique)
Reutilisable: oui -- l'appariement porteur / temoin sur voisinage genetique ET
              profondeur de branche terminale, avec l'evenement comme unite de
              replication, vaut pour toute mutation convergente d'une bacterie
              clonale
Projet      : GC_par_lignee
Date        : 2026-08-30
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import read_subs, canonical  # noqa: E402
from phase2_polarisation_mtbc0 import build_ancestral  # noqa: E402
from phase5_p51_panel_recursif import masked_positions  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
BDD = Path("/home/christophe/docs/codes/mtbc/bdd/actuelle")
def _coords():
    """Coordonnees LUES du GFF3, jamais codees en dur : un premier jet les avait
    ecrites de memoire et se trompait de 141 pb sur Rv0413, ce qui aurait fausse
    le critere D2 sans rien signaler."""
    import re as _re
    g = {"Rv0413": dict(nom="mutT3", dirige=True),
         "Rv3056": dict(nom="dinP", dirige=False)}
    for line in Path("/home/christophe/docs/codes/mtbc/investigate_phylo/"
                     "resources/NC_000962.3.gff3").read_text().splitlines():
        if "\tCDS\t" not in line:
            continue
        f = line.split("\t")
        m = _re.search(r"locus_tag=(Rv0413|Rv3056);", f[8])
        if m:
            g[m.group(1)].update(debut=int(f[3]) - 1, fin=int(f[4]), brin=f[6])
    return g


GENES = _coords()
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
CAP_TEMOINS = 180
_CACHE_SPDI = {}


def subs_cache(f, masked):
    """Les 29 evenements puisent leurs temoins dans des repertoires de clade qui
    se recoupent largement ; sans cache, le meme genome est relu jusqu'a une
    dizaine de fois."""
    k = str(f)
    if k not in _CACHE_SPDI:
        _CACHE_SPDI[k] = {x for x in read_subs(f) if x[0] not in masked}
    return _CACHE_SPDI[k]


def coord_cds(locus, pos):
    """Fraction du CDS parcourue au moment du stop, orientation du gene prise en
    compte (D2)."""
    g = GENES[locus]
    L = g["fin"] - g["debut"]
    x = (pos - g["debut"]) / L if g["brin"] == "+" else (g["fin"] - 1 - pos) / L
    return float(np.clip(x, 0, 1))


def pool_local(porteurs, clades, anc, masked, cap=CAP_TEMOINS, seed=0):
    """`porteurs` : liste de couples (clade, SRA)."""
    """Pool local d'un evenement : tous les porteurs, plus des candidats temoins
    tires dans les MEMES repertoires de clade -- donc le voisinage phylogenetique
    le plus proche que la taxonomie du depot permette de nommer."""
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
    for c, s in cand:                       # tour de role entre repertoires
        rng.shuffle(s)
    i = 0
    while len(temoins) < cap and any(i < len(s) for _, s in cand):
        for c, s in cand:
            if i < len(s) and s[i] not in sra_p:
                temoins.append((c, s[i]))
                if len(temoins) >= cap:
                    break
        i += 1
    noms, subsets, est_porteur = [], [], []
    for c, sra in list(porteurs) + temoins:
        f = BDD / c / sra / "NC_000962.3" / "spdi.txt"
        if not f.exists():
            continue
        v = subs_cache(f, masked)
        if not v:
            continue
        noms.append(sra)
        subsets.append(v)
        est_porteur.append(sra in sra_p)
    return noms, subsets, np.array(est_porteur, bool)


def prives(subsets, anc):
    """Pour chaque souche, les variants PRIVES du pool (branche terminale),
    polarises sur MTBC0, classes en six classes canoniques."""
    sup = defaultdict(list)
    for j, s in enumerate(subsets):
        for v in s:
            sup[v].append(j)
    out = [np.zeros(6, int) for _ in subsets]
    tot = np.zeros(len(subsets), int)
    for (pos, ref, alt), js in sup.items():
        if len(js) != 1:
            continue
        a = chr(anc[pos]) if pos < len(anc) else "N"
        if a == "N" or a == alt:
            continue
        if a != ref:
            ref = a
        tot[js[0]] += 1
        out[js[0]][CLASSES.index(canonical(ref, alt))] += 1
    return np.array(out), tot


def apparier(prof, est_porteur, rng, tol=0.20, essais=3000):
    """Jeu de temoins de meme effectif que les porteurs, choisi pour egaliser la
    PROFONDEUR moyenne (variants prives par souche). C'est le confondant que
    A25 et A26 ont montre actif : un groupe plus dense a des branches terminales
    plus courtes, donc des mutations plus jeunes."""
    ip = np.flatnonzero(est_porteur)
    it = np.flatnonzero(~est_porteur)
    if len(ip) < 2 or len(it) < 10:
        return None, np.nan
    cible = prof[ip].mean()
    best, ecart = None, np.inf
    for _ in range(essais):
        pick = rng.choice(it, min(len(ip), len(it)), replace=False)
        e = abs(prof[pick].mean() - cible) / max(cible, 1e-9)
        if e < ecart:
            best, ecart = pick, e
            if e < 0.01:
                break
    return (best, ecart) if ecart <= tol else (None, ecart)


def part_ca(sp):
    """Statistique primaire : C>A / (C>A + C>T), part de la transversion parmi
    les pertes de paires G:C. Exempte d'opportunite par construction."""
    ca, ct = sp[CLASSES.index("C>A")], sp[CLASSES.index("C>T")]
    return (ca / (ca + ct), ca + ct) if (ca + ct) else (np.nan, 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--sim", type=int, default=4000)
    ap.add_argument("--effet", type=float, default=1.5, help="cote, pour Q2")
    ap.add_argument("--cap", type=int, default=CAP_TEMOINS,
                    help="plafond de temoins candidats par evenement. PARAMETRE "
                         "D'IMPLEMENTATION, absent des criteres pre-enregistres : "
                         "D1 n'exige que >= 10 temoins candidats et un ecart de "
                         "profondeur <= 20 pour cent. Le premier run, a 180, a echoue "
                         "l'appariement sur 10 des 29 evenements, dont les deux "
                         "plus peuples (77 et 72 porteurs) : avec 77 porteurs il "
                         "faut tirer 77 temoins parmi 180, ce qui laisse trop peu "
                         "de latitude pour deplacer la profondeur moyenne. "
                         "L'elargir est AVEUGLE AU RESULTAT (l'exclusion se decide "
                         "avant toute statistique) et les DEUX runs sont rapportes.")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    anc, masked = build_ancestral(verbose=False), masked_positions()
    ev = pd.read_csv(ROOT / "résultats" / "phase8_p71_evenements.tsv", sep="\t")
    var = pd.read_csv(ROOT / "data" / "p71_variants_3r.tsv.gz", sep="\t")
    ev = ev[(ev.locus.isin(GENES)) & (ev.type == "non-sens") & (ev.n_porteurs >= 2)]

    lignes, spectres = [], []
    for e in ev.itertuples():
        sel = var[(var.pos == e.pos) & (var.ref == e.ref) & (var.alt == e.alt)]
        porteurs = [(r.clade, r.sra) for r in sel.itertuples()]   # (clade, SRA)
        clades = {c for c, _ in porteurs}
        frac = coord_cds(e.locus, e.pos)
        noms, subsets, est_p = pool_local(porteurs, clades, anc, masked,
                                          cap=args.cap, seed=args.seed)
        if est_p.sum() < 2 or (~est_p).sum() < 10:
            lignes.append(dict(locus=e.locus, nom=GENES[e.locus]["nom"],
                               pos=e.pos, frac_cds=frac, n_porteurs=int(e.n_porteurs),
                               n_porteurs_lus=int(est_p.sum()),
                               n_temoins_cand=int((~est_p).sum()),
                               qualifie=False, motif="trop peu de porteurs/temoins"))
            continue
        sp, tot = prives(subsets, anc)
        idx_t, ecart = apparier(tot, est_p, rng)
        if idx_t is None:
            lignes.append(dict(locus=e.locus, nom=GENES[e.locus]["nom"],
                               pos=e.pos, frac_cds=frac, n_porteurs=int(e.n_porteurs),
                               n_porteurs_lus=int(est_p.sum()),
                               n_temoins_cand=int((~est_p).sum()),
                               ecart_profondeur=ecart, qualifie=False,
                               motif="profondeur non appariable"))
            continue
        ip = np.flatnonzero(est_p)
        sp_p, sp_t = sp[ip].sum(0), sp[idx_t].sum(0)
        f_p, n_p = part_ca(sp_p)
        f_t, n_t = part_ca(sp_t)
        lignes.append(dict(
            locus=e.locus, nom=GENES[e.locus]["nom"], pos=e.pos, frac_cds=frac,
            n_porteurs=int(e.n_porteurs), n_porteurs_lus=int(est_p.sum()),
            n_temoins_cand=int((~est_p).sum()),
            prof_porteurs=float(tot[ip].mean()), prof_temoins=float(tot[idx_t].mean()),
            ecart_profondeur=ecart, n_pertes_porteurs=int(n_p),
            n_pertes_temoins=int(n_t), part_ca_porteurs=f_p, part_ca_temoins=f_t,
            delta=f_p - f_t, qualifie=bool(n_p >= 5 and n_t >= 5 and frac <= 0.90),
            motif="" if (n_p >= 5 and n_t >= 5) else "moins de 5 pertes d'un cote"))
        for cote, s in (("porteurs", sp_p), ("temoins", sp_t)):
            spectres.append(dict(locus=e.locus, nom=GENES[e.locus]["nom"],
                                 pos=e.pos, cote=cote,
                                 **{c: int(v) for c, v in zip(CLASSES, s)}))

    df = pd.DataFrame(lignes)
    df.to_csv(ROOT / "résultats" /
              f"phase8_p72_evenements_qualifies_cap{args.cap}.tsv",
              sep="\t", index=False)
    pd.DataFrame(spectres).to_csv(ROOT / "résultats" /
                                  f"phase8_p72_spectres_cap{args.cap}.tsv",
                                  sep="\t", index=False)
    print("=== D1 et D2. qualification des evenements ===")
    for locus, g in GENES.items():
        s = df[df.locus == locus]
        q = s[s.qualifie]
        print(f"  {g['nom']:6} : {len(s)} evenements non-sens a >= 2 porteurs, "
              f"{len(q)} qualifiants "
              f"({(s.frac_cds > 0.90).sum()} exclus par D2, stop dans les 10 % "
              f"terminaux du CDS)")
    print()
    print(df[df.qualifie][["nom", "pos", "frac_cds", "n_porteurs_lus",
                           "prof_porteurs", "prof_temoins", "ecart_profondeur",
                           "n_pertes_porteurs", "n_pertes_temoins",
                           "part_ca_porteurs", "part_ca_temoins", "delta"]]
          .round(4).to_string(index=False))

    res = {}
    for locus, g in GENES.items():
        q = df[(df.locus == locus) & df.qualifie].copy()
        if len(q) < 4:
            print(f"\n  {g['nom']} : {len(q)} evenements qualifiants, "
                  f"trop peu pour un test. Aucun verdict.")
            res[g["nom"]] = dict(n=len(q), verdict="non testable")
            continue
        # ---- Q1 resolution de l'appariement
        d_prof = (q.prof_porteurs - q.prof_temoins).to_numpy()
        rho, p_q1 = stats.spearmanr(d_prof, q.delta.to_numpy())
        # ---- Q2 puissance, AVANT le test
        odds = args.effet
        sim = np.zeros(args.sim)
        for b in range(args.sim):
            dd = []
            for r in q.itertuples():
                base = (r.part_ca_porteurs * r.n_pertes_porteurs +
                        r.part_ca_temoins * r.n_pertes_temoins) / \
                       max(r.n_pertes_porteurs + r.n_pertes_temoins, 1)
                o = base / max(1 - base, 1e-9)
                pp = np.clip(o * odds / (1 + o * odds), 0, 1)
                a = rng.binomial(r.n_pertes_porteurs, pp) / max(r.n_pertes_porteurs, 1)
                c = rng.binomial(r.n_pertes_temoins, base) / max(r.n_pertes_temoins, 1)
                dd.append(a - c)
            dd = np.array(dd)
            try:
                sim[b] = stats.wilcoxon(dd, alternative="greater").pvalue
            except ValueError:
                sim[b] = 1.0
        puissance = float((sim <= 0.05).mean())
        # ---- TEST UNIQUE
        alt = "greater" if g["dirige"] else "two-sided"
        try:
            st = stats.wilcoxon(q.delta.to_numpy(), alternative=alt)
            w, p = float(st.statistic), float(st.pvalue)
        except ValueError:
            w, p = np.nan, 1.0
        med = float(np.median(q.delta))
        signes = int((q.delta > 0).sum())
        print(f"\n=== {g['nom']} ({locus}) ===")
        print(f"  Q1 appariement : rho(delta profondeur, delta statistique) = "
              f"{rho:+.3f} (p = {p_q1:.3f}) -> "
              f"{'PASSE' if p_q1 > 0.05 else 'ECHOUE, test non rendu'}")
        print(f"  Q2 PUISSANCE, calculee avant le test : sous une cote de "
              f"x{odds:g} sur la part de C>A parmi les pertes, ce plan verrait "
              f"l'effet dans {100*puissance:.1f} % des cas")
        if puissance < 0.50:
            print(f"     -> sous 50 % : un silence sera declare SILENCE DE "
                  f"RESOLUTION, pas absence d'effet")
        if p_q1 <= 0.05:
            res[g["nom"]] = dict(n=len(q), verdict="Q1 echoue, test non rendu",
                                 puissance=puissance)
            continue
        print(f"  TEST UNIQUE ({'unilateral a droite, prediction dirigee' if g['dirige'] else 'bilateral, EXPLORATOIRE'}) :")
        print(f"     mediane des differences appariees = {med:+.4f}, "
              f"{signes}/{len(q)} evenements de signe positif")
        print(f"     Wilcoxon W = {w:.1f}, p = {p:.4g}")
        verdict = ("EFFET DETECTE" if p <= 0.05 else
                   "silence de resolution" if puissance < 0.50 else
                   "PAS D'EFFET, a resolution suffisante")
        print(f"     verdict : {verdict}")
        res[g["nom"]] = dict(n=len(q), rho_q1=rho, p_q1=p_q1,
                             puissance=puissance, mediane_delta=med,
                             signes_positifs=signes, W=w, p=p, verdict=verdict)
    pd.DataFrame([dict(gene=k, **v) for k, v in res.items()]).to_csv(
        ROOT / "résultats" / f"phase8_p72_verdict_cap{args.cap}.tsv",
        sep="\t", index=False)
    print(f"\n=== VERDICT P7.2 ===")
    for k, v in res.items():
        print(f"  {k:6} : {v['verdict']}"
              + (f" (n = {v['n']} evenements, puissance "
                 f"{100*v.get('puissance', float('nan')):.0f} %)"
                 if "puissance" in v else f" (n = {v['n']})"))


if __name__ == "__main__":
    main()
