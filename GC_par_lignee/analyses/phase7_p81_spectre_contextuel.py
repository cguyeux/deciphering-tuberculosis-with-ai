#!/usr/bin/env python3
"""
Objet       : P8.1 -- decomposer le spectre mutationnel de chaque lignee en
              6 classes canoniques x 16 contextes trinucleotidiques (96 canaux),
              sur EVENEMENTS polarises sur MTBC0 (donc apres P2.4 et P2.5), et
              P8.2 -- dire si l'ecart inter-lignees, une fois retiree la part
              deja etablie a la maille des 6 classes (A24), se CONCENTRE dans
              quelques contextes ou se DISTRIBUE uniformement. Un ecart
              concentre oriente vers un processus enzymatique identifie ; un
              ecart diffus vers un biais global de replication ou un artefact.

              POURQUOI LE CONTEXTE, ET POURQUOI MAINTENANT. A24 a etabli deux
              choses : le spectre a 6 classes n'est pas homogene entre lignees,
              et il porte DEUX axes orthogonaux, le flux GC (objet du projet) et
              un axe Ts/Tv d'amplitude comparable (1,35 chez Microti a 3,08 chez
              Bovis) que le projet constate sans l'expliquer -- lacune n°11 de
              l'etat. Le contexte trinucleotidique est la premiere resolution
              superieure disponible, et c'est celle ou les mecanismes de
              reparation candidats (Ogt sur O6-methylguanine, MutT sur 8-oxo-dG,
              Ung sur l'uracile de desamination) laissent des empreintes
              distinctes. Le volet METHYLOME de la proposition P8 initiale reste
              ABANDONNE (pas de 5-methylcytosine chez les mycobacteries,
              Hemavathy 1995 ; methylation N6-adenine, Zhu 2016) : la section F1
              en fait un CONTROLE POSITIF plutot qu'une hypothese.

              CRITERES PRE-ENREGISTRES, FIXES AVANT DE REGARDER LES RESULTATS
              (meme discipline qu'en P10.1 et P10.2, et pour la meme raison :
              96 canaux offrent 96 occasions de trouver quelque chose)
                C1 REPORTABILITE. Un canal n'est interpretable INDIVIDUELLEMENT
                   que s'il porte >= 100 evenements toutes lignees confondues ET
                   >= 5 evenements dans au moins 13 des 17 lignees. Les autres
                   entrent dans les tests globaux mais ne sont jamais nommes.
                C2 EXISTENCE D'UN ECART CONTEXTUEL. Le nul valide est la
                   permutation des etiquettes de lignee ENTRE SOUCHES (jamais
                   entre substitutions : A24 a etabli que le p asymptotique,
                   qui traite chaque substitution comme independante, ne doit
                   pas etre cite). Statistique : somme sur les 6 classes du G2
                   d'independance de la table lignee x contexte. Le modele nul
                   M_class donne donc a chaque lignee son propre taux par
                   CLASSE -- toute la part deja acquise en A24 est absorbee --
                   et ne teste QUE le profil contextuel interne a chaque classe.
                   Seuil : p_permutation <= 0,05. Si C2 echoue, P8.2 est sans
                   objet et le verdict est ecrit tel quel.
                C3 CONCENTRATION (le test de P8.2). Part du G2 total portee par
                   les 5 canaux les plus contributeurs sur 96. Nul DIFFUS :
                   meme quantite totale d'heterogeneite, repartie egalement sur
                   les 96 canaux, puis bruit de Poisson. Verdict CONCENTRE si la
                   part observee depasse le q95 du nul diffus, DIFFUS sinon.
                C4 PUISSANCE DU TEST C3, epreuve de validite. Une simulation
                   CONCENTREE (toute l'heterogeneite dans 5 canaux) doit se
                   separer du nul diffus (q05 concentre > q95 diffus). Sinon le
                   verdict de C3 est INDETERMINE et doit etre annonce comme tel,
                   quel que soit son cote.
                C5 CONFONDANT DE COMPOSITION. L'heterogeneite par canal ne doit
                   pas etre monotone dans le nombre de G/C des flancs (0, 1, 2) :
                   une monotonie signerait un echo du gradient de GC local, donc
                   du biais de couverture (residu de +7 % laisse ouvert par P9.4),
                   et non un processus de contexte. Diagnostic, non qualifiant.
              Aucun autre decoupage du spectre ne sera essaye dans ce script.

              CE QUE LE SCRIPT NE FAIT PAS. Il ne teste pas les canaux contre S
              en confirmatoire : apres A20, toute correlation entre une mesure
              du spectre et la force requise est DESCRIPTIVE et declaree telle
              (section F3). Rouvrir une inference causale demanderait une piste
              neuve, pas un 97e canal.

Entrees     : bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt
              data/MTBC0/ancestral_on_H37Rv.bin (polarisation, P2.4)
              data/mask_h37rv_positions.npy    (masque MTBC0, P3.1)
              résultats/phase5_p53_force_requise.tsv (S par lignee, A18)
Sorties     : résultats/phase7_p81_opportunite_trinuc.tsv  (A)
              résultats/phase7_p81_counts_par_souche.tsv    (B, souche x canal)
              résultats/phase7_p81_spectre_96.tsv           (C, lignee x canal)
              résultats/phase7_p81_ecart_contextuel.tsv     (D, C2)
              résultats/phase7_p81_concentration.tsv        (E, C3 + C4)
              résultats/phase7_p81_controles.tsv            (F, C5 + F1/F2/F3)
Reutilisable: oui -- l'opportunite trinucleotidique vectorisee, la collapse
              pyrimidine du contexte et le couple « existence par permutation de
              souches / concentration contre nul diffus » valent pour tout
              spectre mutationnel de bacterie clonale
Projet      : GC_par_lignee
Date        : 2026-08-30
"""
import argparse
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import read_subs, flux, canonical  # noqa: E402
from phase2_polarisation_mtbc0 import build_ancestral  # noqa: E402
from phase5_p51_panel_recursif import PANEL, masked_positions, strains_of  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
BASES = "ACGT"                     # index i, complement = 3 - i (A<->T, C<->G)
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
LOSS, GAIN = {"C>A", "C>T"}, {"T>C", "T>G"}
CTX = [f"{a}.{b}" for a in BASES for b in BASES]      # 5' . 3'
CANAUX = [f"{k}@{c}" for k in CLASSES for c in CTX]   # 96 canaux
FAIBLES = {"L8", "Suricattae", "Dassie", "Pinnipedii"}


# --------------------------------------------------------------- A. opportunite
def trinuc_opportunity(anc, masked):
    """Nombre de sites disponibles par trinucleotide reference pyrimidine, sur
    le genome ancestral non masque. Le centre doit etre non masque et les trois
    bases connues ; meme critere que pour les evenements, sans quoi le taux par
    site n'est pas un taux."""
    arr = np.frombuffer(bytes(anc), dtype=np.uint8)
    n = len(arr)
    code = np.full(n, -1, np.int8)
    for i, b in enumerate(BASES):
        code[arr == ord(b)] = i
    keep = np.zeros(n, bool)
    keep[1:n - 1] = (code[:-2] >= 0) & (code[1:-1] >= 0) & (code[2:] >= 0)
    m = np.fromiter(masked, np.int64, len(masked))
    keep[m[(m >= 0) & (m < n)]] = False
    idx = np.nonzero(keep)[0]
    l5, ce, r3 = code[idx - 1].astype(int), code[idx].astype(int), code[idx + 1].astype(int)
    pur = (ce == 0) | (ce == 2)                        # A ou G : on retourne
    L = np.where(pur, 3 - r3, l5)
    C = np.where(pur, 3 - ce, ce)
    R = np.where(pur, 3 - l5, r3)
    cnt = np.bincount(C * 16 + L * 4 + R, minlength=64)
    opp = {}
    for c in (1, 3):                                   # C, T
        for a in range(4):
            for b in range(4):
                opp[f"{BASES[c]}@{BASES[a]}.{BASES[b]}"] = int(cnt[c * 16 + a * 4 + b])
    return opp


# ------------------------------------------------------- B. evenements polarises
def strain_channels(strains, anc, masked, n, seed=0):
    """Comptes par canal pour chaque souche d'un pool de n souches. Un variant
    porte par une seule souche = un evenement de branche terminale (P2.5).

    POURQUOI SEULEMENT LES BRANCHES TERMINALES, alors que les branches internes
    tripleraient l'effectif. Le seul nul valide pour ce test est la permutation
    des etiquettes de lignee ENTRE SOUCHES (A24), qui suppose la souche
    echangeable. Un evenement de branche interne n'appartient pas a une souche
    mais a un GROUPE de souches : il n'a pas de place dans la table souche x
    canal et le permuter reviendrait a casser precisement le groupement dont le
    nul doit tenir compte. Y adjoindre les branches internes demanderait un nul
    au niveau de la BRANCHE, donc un plan different, pas une option. La
    puissance se gagne ici en faisant varier n (sensibilite ci-dessous), pas en
    melangeant deux unites d'observation."""
    rng = random.Random(seed)      # meme tirage que P5.1 : pools IDENTIQUES
    pool = list(strains)
    rng.shuffle(pool)
    subsets, names = [], []
    for s in pool[:n]:
        v = read_subs(s / "NC_000962.3" / "spdi.txt")
        if v:
            subsets.append(v)
            names.append(s.name)
    k = len(subsets)
    if k < 4:
        return None, k
    support = defaultdict(int)
    for i, subs in enumerate(subsets):
        for v in subs:
            support[v] |= 1 << i
    per = defaultdict(Counter)
    lg = len(anc)
    for (pos, ref, alt), mask in support.items():
        if mask & (mask - 1):                          # non terminal
            continue
        if pos in masked or pos < 1 or pos >= lg - 1:
            continue
        a = chr(anc[pos])
        if a == "N" or a == alt:
            continue                                   # non lifte, ou inverse
        tierce = a != ref
        ref = a
        l5, r3 = chr(anc[pos - 1]), chr(anc[pos + 1])
        if l5 not in BASES or r3 not in BASES:
            continue
        if ref in "CT":
            ctx = f"{l5}.{r3}"
        else:                                          # retournement de brin
            ctx = f"{BASES[3 - BASES.index(r3)]}.{BASES[3 - BASES.index(l5)]}"
        c = per[mask.bit_length() - 1]
        c[f"{canonical(ref, alt)}@{ctx}"] += 1
        c["_n"] += 1
        c["_tierce"] += tierce
    rows = []
    for i, name in enumerate(names):
        c = per[i]
        rows.append(dict(sra=name, n_pool=k, n_evt=c["_n"], n_tierce=c["_tierce"],
                         **{ca: c[ca] for ca in CANAUX}))
    return pd.DataFrame(rows), k


# ------------------------------------------------- D/E. modele M_class et nul
def g2_independance(tab):
    """G2 d'independance d'une table lignee x contexte, et contribution de
    chaque contexte. C'est la deviance du modele M_class : chaque lignee a son
    propre taux pour la CLASSE, le profil contextuel interne est universel.
    L'opportunite trinucleotidique est absorbee par l'effet de colonne, donc le
    test d'ecart contextuel n'en depend pas -- elle ne sert qu'aux taux."""
    tab = np.asarray(tab, float)
    tot = tab.sum()
    if tot <= 0:
        return 0.0, np.zeros(tab.shape[1]), np.zeros_like(tab)
    exp = np.outer(tab.sum(1), tab.sum(0)) / tot
    ok = exp > 0
    term = np.zeros_like(tab)
    nz = ok & (tab > 0)
    term[nz] = 2 * tab[nz] * np.log(tab[nz] / exp[nz])
    return float(term.sum()), term.sum(0), exp


def g2_from_agg(agg):
    """Somme sur les 6 classes du G2 d'independance lignee x contexte, et
    contribution des 96 canaux. `agg` : matrice (lignee x 96 canaux)."""
    tot, contrib = 0.0, np.zeros(96)
    for ki in range(6):
        sub = agg[:, ki * 16:(ki + 1) * 16]
        g, per_col, _ = g2_independance(sub)
        tot += g
        contrib[ki * 16:(ki + 1) * 16] = per_col
    return tot, contrib


def agreger(mat, lab_idx, n_lin):
    out = np.zeros((n_lin, mat.shape[1]))
    np.add.at(out, lab_idx, mat)
    return out


def simule(mat, lab_idx, n_lin, cible_g2, mode, rng, n_top=5):
    """Nul calibre. On part d'une table obtenue par PERMUTATION des etiquettes
    de lignee entre souches -- elle porte donc exactement la sur-dispersion de
    la donnee reelle, que le groupement des substitutions par souche fabrique --
    puis on y INJECTE une deviation systematique de G2 total cible_g2,
    repartie egalement sur les 96 canaux (mode 'diffus') ou concentree sur
    n_top canaux (mode 'concentre'). Le G2 total des deux regimes est donc le
    meme, et seule leur REPARTITION entre canaux differe : c'est exactement la
    question de P8.2, et rien d'autre."""
    agg = agreger(mat, rng.permutation(lab_idx), n_lin)
    exp = np.zeros_like(agg)
    for ki in range(6):
        sub = agg[:, ki * 16:(ki + 1) * 16]
        t = sub.sum()
        if t > 0:
            exp[:, ki * 16:(ki + 1) * 16] = np.outer(sub.sum(1), sub.sum(0)) / t
    cible = np.zeros(96)
    if mode == "diffus":
        cible[:] = cible_g2 / 96
    else:
        cible[rng.choice(96, n_top, replace=False)] = cible_g2 / n_top
    pert = rng.normal(size=agg.shape)
    pert -= pert.mean(0)                        # deviation d'interaction, centree
    for c in range(96):
        norm = (pert[:, c] ** 2).sum()
        if exp[:, c].sum() <= 0 or cible[c] <= 0 or norm <= 0:
            pert[:, c] = 0
        else:
            pert[:, c] *= np.sqrt(cible[c] / norm)
    sim = np.maximum(agg + pert * np.sqrt(np.maximum(exp, 1e-9)), 0.0)
    g, contrib = g2_from_agg(sim)
    if g <= 0:
        return np.nan, np.nan
    return float(np.sort(contrib)[::-1][:n_top].sum() / g), g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--n-per-clade", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--perm", type=int, default=5000)
    ap.add_argument("--sim", type=int, default=2000)
    ap.add_argument("--out-prefix",
                    default=str(ROOT / "résultats" / "phase7_p81_"))
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    anc = build_ancestral()
    masked = masked_positions()

    # ---- A. opportunite trinucleotidique
    opp = trinuc_opportunity(anc, masked)
    oppd = pd.DataFrame([dict(trinuc=k, sites=v) for k, v in opp.items()])
    oppd["centre"] = oppd["trinuc"].str[0]
    oppd["flancs_gc"] = oppd["trinuc"].str[2:].str.count("[GC]")
    oppd.to_csv(args.out_prefix + "opportunite_trinuc.tsv", sep="\t", index=False)
    print("=== A. opportunite : sites disponibles par trinucleotide "
          "(centre pyrimidine, genome ancestral non masque) ===")
    print(f"  total centre C : {oppd[oppd.centre=='C'].sites.sum():,} ; "
          f"centre T : {oppd[oppd.centre=='T'].sites.sum():,}")
    print(f"  etendue par contexte : {oppd.sites.min():,} (" 
          f"{oppd.loc[oppd.sites.idxmin(),'trinuc']}) a {oppd.sites.max():,} ("
          f"{oppd.loc[oppd.sites.idxmax(),'trinuc']}), facteur "
          f"{oppd.sites.max()/oppd.sites.min():.1f} -- c'est ce facteur que la "
          "normalisation retire et qu'un spectre brut confondrait avec un taux")

    # ---- B. comptes par souche
    rows = []
    for lin, prefixes in PANEL.items():
        st = strains_of(prefixes)
        if len(st) < 4:
            continue
        df, k = strain_channels(st, anc, masked, args.n_per_clade, args.seed)
        if df is None:
            print(f"# {lin} : pool trop petit ({k})", file=sys.stderr)
            continue
        df.insert(0, "lignee", lin)
        rows.append(df)
    per_strain = pd.concat(rows, ignore_index=True)
    per_strain.to_csv(args.out_prefix + "counts_par_souche.tsv", sep="\t",
                      index=False)
    lignees = [l for l in PANEL if l in set(per_strain["lignee"])]
    agg = per_strain.groupby("lignee")[CANAUX].sum().reindex(lignees)
    print(f"\n=== B. {len(lignees)} lignees, "
          f"{int(per_strain['n_evt'].sum()):,} evenements de branche terminale "
          f"polarises, {int(per_strain['n_tierce'].sum()):,} etats tierces "
          f"({100*per_strain['n_tierce'].sum()/max(per_strain['n_evt'].sum(),1):.2f} %) ===")

    # ---- C. spectre 96 canaux, brut et normalise par l'opportunite
    tot_evt = agg.sum(1)
    spec = []
    for lin in lignees:
        for ca in CANAUX:
            k, c = ca.split("@")
            o = opp[f"{k[0]}@{c}"]
            n = int(agg.loc[lin, ca])
            spec.append(dict(lignee=lin, canal=ca, classe=k, contexte=c,
                             n=n, part=n / max(tot_evt[lin], 1),
                             taux_par_site=n / o,
                             taux_relatif=(n / o) / max(tot_evt[lin], 1),
                             flancs_gc=c.count("G") + c.count("C")))
    spec = pd.DataFrame(spec)
    spec.to_csv(args.out_prefix + "spectre_96.tsv", sep="\t", index=False)
    # C1 : reportabilite
    tot_par_canal = agg.sum(0)
    n_lin_ok = (agg >= 5).sum(0)
    reportable = (tot_par_canal >= 100) & (n_lin_ok >= 13)
    print(f"\n=== C. spectre a 96 canaux ===")
    print(f"  C1 : {int(reportable.sum())} canaux sur 96 reportables "
          f"(>= 100 evenements au total et >= 5 dans >= 13 lignees)")
    piv = spec.pivot_table(index="lignee", columns="canal", values="taux_relatif")
    piv = piv.reindex(lignees)[CANAUX]
    v = piv.to_numpy()
    v = v / np.linalg.norm(v, axis=1, keepdims=True)
    cos = pd.DataFrame(v @ v.T, index=lignees, columns=lignees)
    off = cos.to_numpy()[~np.eye(len(lignees), dtype=bool)]
    print(f"  similarite cosinus entre spectres normalises : mediane "
          f"{np.median(off):.4f}, min {off.min():.4f} "
          f"({cos.stack().idxmin()})")

    # ---- D. C2 : y a-t-il un ecart CONTEXTUEL, une fois A24 absorbe ?
    mat = per_strain[CANAUX].to_numpy(float)
    ilin = {l: i for i, l in enumerate(lignees)}
    lab_idx = per_strain["lignee"].map(ilin).to_numpy()
    aggn = agreger(mat, lab_idx, len(lignees))
    g_obs, contrib_a = g2_from_agg(aggn)
    contrib = pd.Series(contrib_a, index=CANAUX)
    perm = np.empty(args.perm)
    for b in range(args.perm):
        perm[b] = g2_from_agg(agreger(mat, rng.permutation(lab_idx),
                                      len(lignees)))[0]
    p_perm = (1 + (perm >= g_obs).sum()) / (1 + args.perm)
    # par classe, pour voir QUELLE classe porte l'heterogeneite contextuelle
    par_classe = []
    for k in CLASSES:
        cols = [f"{k}@{c}" for c in CTX]
        g, _, _ = g2_independance(agg[cols].to_numpy())
        par_classe.append(dict(classe=k, n=int(agg[cols].to_numpy().sum()),
                               g2=g, ddl=(len(lignees) - 1) * 15,
                               g2_par_ddl=g / ((len(lignees) - 1) * 15)))
    par_classe = pd.DataFrame(par_classe).sort_values("g2_par_ddl", ascending=False)
    d = pd.DataFrame([dict(g2_observe=g_obs, ddl_total=6 * (len(lignees) - 1) * 15,
                           p_permutation=p_perm, perm_median=float(np.median(perm)),
                           perm_q95=float(np.percentile(perm, 95)),
                           perm_q99=float(np.percentile(perm, 99)),
                           n_permutations=args.perm)])
    d.to_csv(args.out_prefix + "ecart_contextuel.tsv", sep="\t", index=False)
    par_classe.to_csv(args.out_prefix + "ecart_contextuel_par_classe.tsv",
                      sep="\t", index=False)
    print("\n=== D. C2 : existe-t-il un ecart CONTEXTUEL inter-lignees ? ===")
    print("  modele nul M_class : chaque lignee a son propre taux par classe "
          "(toute la part A24 absorbee), profil contextuel universel")
    print(f"  G2 observe = {g_obs:.1f} sur {6*(len(lignees)-1)*15} ddl")
    print(f"  nul par permutation des etiquettes ENTRE SOUCHES "
          f"({args.perm} tirages) : median {np.median(perm):.1f}, "
          f"q95 {np.percentile(perm,95):.1f}, q99 {np.percentile(perm,99):.1f}")
    print(f"  p = {p_perm:.4g}  -> C2 {'PASSE' if p_perm <= 0.05 else 'ECHOUE'}")
    print("\n  heterogeneite contextuelle par classe (G2 par ddl) :")
    print(par_classe.round(3).to_string(index=False))

    # ---- E. C3/C4 : concentration contre diffusion (P8.2)
    ordre = contrib.sort_values(ascending=False)
    top5 = float(ordre.iloc[:5].sum() / g_obs) if g_obs > 0 else np.nan
    if p_perm <= 0.05:
        exces = max(g_obs - float(np.median(perm)), 0.0)
        rd = np.array([simule(mat, lab_idx, len(lignees), exces, "diffus", rng)
                       for _ in range(args.sim)])
        rc = np.array([simule(mat, lab_idx, len(lignees), exces, "concentre", rng)
                       for _ in range(args.sim)])
        sd, gd = rd[:, 0], rd[:, 1]
        sc, gc_ = rc[:, 0], rc[:, 1]
        sd, sc = sd[~np.isnan(sd)], sc[~np.isnan(sc)]
        print(f"  calibration : G2 observe {g_obs:.1f} ; G2 des simulations, "
              f"median {np.nanmedian(gd):.1f} (diffus) et {np.nanmedian(gc_):.1f} "
              f"(concentre) -- les trois doivent coincider, seule la "
              f"REPARTITION entre canaux differe")
        q95d, q05c = np.percentile(sd, 95), np.percentile(sc, 5)
        c4 = q05c > q95d
        verdict = ("INDETERMINE (C4 echoue : le test ne separe pas les deux "
                   "regimes)" if not c4 else
                   "CONCENTRE" if top5 > q95d else "DIFFUS")
    else:
        sd = sc = np.array([np.nan])
        q95d = q05c = np.nan
        c4 = False
        verdict = "SANS OBJET (C2 echoue)"
    e = pd.DataFrame([dict(part_top5_observee=top5, exces_g2=g_obs - float(np.median(perm)),
                           nul_diffus_median=float(np.nanmedian(sd)),
                           nul_diffus_q95=float(q95d),
                           alt_concentree_median=float(np.nanmedian(sc)),
                           alt_concentree_q05=float(q05c),
                           C4_puissance=bool(c4), verdict=verdict,
                           n_simulations=args.sim)])
    e.to_csv(args.out_prefix + "concentration.tsv", sep="\t", index=False)
    ordre.rename("g2_contribution").to_frame().assign(
        part=lambda x: x.g2_contribution / g_obs,
        reportable=reportable.reindex(ordre.index).to_numpy(),
        n=tot_par_canal.reindex(ordre.index).to_numpy()
    ).to_csv(args.out_prefix + "canaux_classes.tsv", sep="\t")
    print("\n=== E. C3/C4 : l'ecart contextuel est-il CONCENTRE ou DIFFUS ? ===")
    print(f"  part du G2 portee par les 5 premiers canaux : {top5:.3f}")
    print(f"  nul DIFFUS (meme total, reparti sur 96 canaux) : median "
          f"{np.nanmedian(sd):.3f}, q95 {q95d:.3f}")
    print(f"  alternative CONCENTREE (tout dans 5 canaux) : median "
          f"{np.nanmedian(sc):.3f}, q05 {q05c:.3f}")
    print(f"  C4 (le test separe-t-il les deux regimes ?) : "
          f"{'PASSE' if c4 else 'ECHOUE'}")
    print(f"  VERDICT P8.2 : {verdict}")
    print("\n  dix premiers canaux par contribution au G2 :")
    tete = ordre.iloc[:10].to_frame("g2")
    tete["part"] = tete["g2"] / g_obs
    tete["n_total"] = tot_par_canal.reindex(tete.index)
    tete["reportable"] = reportable.reindex(tete.index)
    print(tete.round(4).to_string())

    # ---- F. controles
    ctrl = []
    # F1 controle positif : pas de 5mC chez les mycobacteries -> pas d'exces
    # de C>T en contexte NCG (la signature de desamination de la 5mC)
    for lin in lignees:
        ncg = sum(agg.loc[lin, f"C>T@{a}.G"] for a in BASES)
        oncg = sum(opp[f"C@{a}.G"] for a in BASES)
        aut = sum(agg.loc[lin, f"C>T@{a}.{b}"] for a in BASES for b in BASES if b != "G")
        oaut = sum(opp[f"C@{a}.{b}"] for a in BASES for b in BASES if b != "G")
        r = (ncg / oncg) / (aut / oaut) if aut and ncg else np.nan
        ctrl.append(dict(controle="F1_CpG", lignee=lin, valeur=r,
                         n=int(ncg + aut)))
    f1 = pd.Series([c["valeur"] for c in ctrl if c["controle"] == "F1_CpG"])
    # F2 : quelle classe porte l'heterogeneite contextuelle (deja en par_classe)
    # F3 DESCRIPTIF : lien des taux par canal avec S (jamais confirmatoire)
    S = pd.read_csv(ROOT / "résultats" / "phase5_p53_force_requise.tsv",
                    sep="\t").set_index("lignee")["S"]
    comm = [l for l in lignees if l in S.index and l not in FAIBLES]
    f3 = []
    for ca in CANAUX:
        if not reportable[ca]:
            continue
        k = ca.split("@")[0]
        den_k = agg.loc[comm, [f"{k}@{c}" for c in CTX]].sum(1)
        x = (agg.loc[comm, ca] / den_k.replace(0, np.nan)).to_numpy()
        if np.isnan(x).any():
            continue
        rho, p = stats.spearmanr(x, S[comm].to_numpy())
        f3.append(dict(controle="F3_descriptif", canal=ca, part_dans_classe=True,
                       rho=rho, p=p, n_lignees=len(comm)))
    f3 = pd.DataFrame(f3)
    if len(f3):
        o = np.argsort(f3["p"].to_numpy())
        q = np.empty(len(f3))
        pv = f3["p"].to_numpy()[o]
        q[o] = np.minimum.accumulate(
            (pv * len(f3) / np.arange(1, len(f3) + 1))[::-1])[::-1]
        f3["q_bh"] = np.clip(q, 0, 1)
        f3 = f3.sort_values("p")
    # F5 (C5) : l'heterogeneite est-elle monotone dans le GC des flancs ?
    gcf = pd.Series({ca: ca.split("@")[1].count("G") + ca.split("@")[1].count("C")
                     for ca in CANAUX})
    xg = gcf.reindex(CANAUX).to_numpy(float)
    yc = contrib.reindex(CANAUX).to_numpy(float)
    zn = np.log1p(tot_par_canal.reindex(CANAUX).to_numpy(float))
    rho5, p5 = stats.spearmanr(xg, yc)
    rx = stats.rankdata(xg); ry = stats.rankdata(yc); rz = stats.rankdata(zn)
    rxz = np.corrcoef(rx, rz)[0, 1]; ryz = np.corrcoef(ry, rz)[0, 1]
    rxy = np.corrcoef(rx, ry)[0, 1]
    den = np.sqrt(max((1 - rxz ** 2) * (1 - ryz ** 2), 1e-12))
    rho5p = (rxy - rxz * ryz) / den
    tstat = rho5p * np.sqrt((96 - 3) / max(1 - rho5p ** 2, 1e-12))
    p5p = 2 * stats.t.sf(abs(tstat), 96 - 3)
    pd.DataFrame(ctrl).to_csv(args.out_prefix + "controles.tsv", sep="\t",
                              index=False)
    f3.to_csv(args.out_prefix + "controles_f3_descriptif.tsv", sep="\t",
              index=False)
    print("\n=== F. controles ===")
    print(f"  F1 CONTROLE POSITIF (pas de 5mC chez les mycobacteries, donc pas "
          f"de signature CpG attendue) :")
    print(f"     taux de C>T en contexte NCG / hors NCG, par lignee : mediane "
          f"{f1.median():.3f} (min {f1.min():.3f}, max {f1.max():.3f})")
    print(f"     attendu ~1 si le mecanisme absent est bien absent ; un exces "
          f"net invaliderait soit la methode soit Hemavathy 1995")
    print(f"  C5 diagnostic, monotonie de l'heterogeneite dans le GC des "
          f"flancs : rho brut = {rho5:+.3f} (p = {p5:.3g}) ; PARTIEL a effectif "
          f"du canal tenu constant = {rho5p:+.3f} (p = {p5p:.3g})")
    print(f"     le brut confond contexte et opportunite (les contextes "
          f"GC-riches sont jusqu'a 19 fois plus abondants) ; c'est le partiel "
          f"qui vaut diagnostic")
    print(f"     un rho fort signerait un echo du gradient de GC local "
          f"(residu P9.4), donc un artefact plutot qu'un contexte")
    if len(f3):
        print(f"  F3 DESCRIPTIF (non confirmatoire, cf. A20) : part du canal "
              f"DANS SA CLASSE contre S, sur {len(comm)} lignees bien pourvues, "
              f"{len(f3)} canaux reportables testes, correction BH :")
        print(f3.head(5)[["canal", "rho", "p", "q_bh"]].round(5)
              .to_string(index=False))
        print(f"     {int((f3['q_bh'] <= 0.05).sum())} canaux a q_BH <= 0,05. "
              f"Descriptif par pre-enregistrement : aucune inference causale "
              f"n'en est tiree ici (cf. A20).")


if __name__ == "__main__":
    main()
