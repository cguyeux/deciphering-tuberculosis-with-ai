#!/usr/bin/env python3
"""
Objet       : P7.1 -- recenser dans la base locale les genes de REPARATION,
              RECOMBINAISON et REPLICATION (3R) frappes plusieurs fois de facon
              PHYLOGENETIQUEMENT INDEPENDANTE, seule construction ou l'allele de
              reparation se separe du fond genetique de sa lignee.

              POURQUOI CETTE CONTRAINTE, ET PAS UNE COMPARAISON PORTEURS /
              NON-PORTEURS. Le crible /challenge du 2026-08-29 avait juge
              redhibitoire la forme naive de P7 : les alleles decrits dans la
              litterature (mutT2 G58R, mutT4 et ogt chez W-Beijing, ogt44 +
              ung501 chez Haarlem, ada/alkA chez M. tuberculosis et M. bovis)
              ont chacun ete acquis UNE SEULE FOIS, a la base de leur lignee.
              « Porteur de ogt44 » est alors synonyme de « Haarlem », et l'effet
              de l'allele est totalement confondu avec le fond genetique et
              l'histoire demographique de la lignee ; un denominateur de
              milliers de souches n'y change rien, le denominateur EFFECTIF vaut
              UN evenement evolutif. Seule la CONVERGENCE fournit des replicats.

              CE QUE P8 VIENT D'AJOUTER, ET QUI ORIENTE LA RECHERCHE. A29 etablit
              que l'heterogeneite contextuelle inter-lignees est portee par les
              GAINS de paires G:C, et A31 que ce cote est aussi le mieux mesure et
              qu'il ne suit pas la force de maintien S. Le mecanisme cherche doit
              donc agir sur les gains, ce qui est precisement la composante que la
              perte de NucS inverse chez la mycobacterie (Castaneda-Garcia 2020).
              P7.1 ne teste pas cette prediction -- elle recense les candidats sur
              lesquels P7.2 pourrait la tester.

              CRITERES PRE-ENREGISTRES, FIXES AVANT TOUTE EXECUTION
                Q1 CONTROLES POSITIFS. Les temoins PRIMAIRES sont les genes dont
                   la PERTE DE FONCTION est le mecanisme de resistance etabli :
                   pncA (Rv2043c) et ethA (Rv3854c), enzymes d'activation de
                   prodrogue dont l'inactivation abolit la conversion du
                   pyrazinamide et de l'ethionamide, et gid (Rv3919c). Chacun
                   doit ressortir a q_BH <= 0,05 et dans le premier decile des
                   rangs. Temoin SECONDAIRE : katG (Rv1908c), attendu seulement
                   significativement enrichi, la resistance a l'isoniazide etant
                   dominee par le faux-sens S315T qui conserve la proteine. Si
                   les temoins primaires echouent, le detecteur est casse et
                   AUCUN verdict n'est rendu.
                Q2 RECURRENCE. Un gene 3R qualifie s'il porte >= 3 evenements
                   d'inactivation distincts (SPDI distincts, non-sens ou
                   decalage de cadre), chacun porte par >= 2 souches -- le seuil
                   a deux souches ecarte l'erreur de sequencage d'une souche
                   isolee -- et sur des ensembles de porteurs disjoints.
                Q3 EXCES SUR L'ATTENDU NEUTRE. La recurrence seule ne suffit pas :
                   un gene long est frappe plus souvent. Le gene doit aussi
                   depasser son attendu de non-sens, calcule depuis son
                   OPPORTUNITE en codons stop ponderee par le spectre mutationnel
                   du MTBC, et calibre sur un jeu TEMOIN de 200 genes tires au
                   hasard hors 3R et hors resistance. q_BH <= 0,05.
                Q4 EXCLUSION DE LA SELECTION MEDICAMENTEUSE. Un gene cible ou
                   determinant de resistance ne peut jamais qualifier comme gene
                   3R : sa convergence est produite par le traitement, pas par la
                   biologie de la reparation. gyrA et gyrB sont donc exclus du
                   jeu 3R bien qu'ils soient replicatifs, et katG et pncA n'y
                   entrent que comme temoins.
                SECONDE PASSE, ET CE QUI A ETE VU AVANT DE L'ECRIRE. La
                   premiere passe a ECHOUE Q1 telle que specifiee et n'a donc
                   rendu aucun verdict. Trois defauts ont ete diagnostiques,
                   dont deux sont des BOGUES DE MISE EN OEUVRE contre la
                   specification ci-dessus, verifiables en confrontant cet
                   en-tete au code, et un est un defaut de modele MESURE dans
                   les donnees.
                   (i) Q3 etait calcule sur les non-sens ET les decalages de
                       cadre, alors que l'en-tete le specifie sur les seuls
                       non-sens et sur l'opportunite en codons stop. La
                       consequence etait absurde et visible : pncA affichait
                       304 evenements pour 45 sites-stop disponibles, soit plus
                       d'evenements que de sites, parce que ses decalages de
                       cadre etaient comptes contre un denominateur qui ne les
                       couvre pas.
                   (ii) L'opportunite n'etait PAS ponderee par le spectre
                       mutationnel, alors que l'en-tete le specifie. Elle l'est
                       desormais, canal par canal, depuis le spectre a 96 canaux
                       de A29 : un site dont le stop demande un C>T n'a pas la
                       meme probabilite qu'un site qui demande un T>A.
                   (iii) A 142 706 genomes le paysage des non-sens est SATURE :
                       22,5 % de tous les sites-stop possibles des genes temoins
                       sont deja frappes chez au moins deux souches. Le compte
                       d'evenements distincts est donc borne par l'opportunite,
                       ce qu'un modele de Poisson ignore ; il est remplace par
                       un modele BINOMIAL sur l'opportunite, qui est le modele
                       correct sous saturation.
                   Les seuils de Q2 ne sont PAS touches. Les resultats 3R de la
                   premiere passe ont ete vus avant l'ecriture de celle-ci, et
                   c'est dit ici parce que cela affaiblit le pre-enregistrement :
                   aucune de ces trois corrections ne deplace un seuil pour faire
                   passer un gene, mais le lecteur doit pouvoir en juger.
                   Q1 est aussi RESPECIFIE, sur une base mecaniste externe aux
                   donnees : un detecteur de perte de fonction doit trouver les
                   genes dont la PERTE est selectionnee, ce qui est le cas des
                   enzymes d'ACTIVATION DE PRODROGUE, pncA pour le pyrazinamide
                   et ethA pour l'ethionamide, et de gid pour la streptomycine.
                   katG n'en fait pas partie : la resistance a l'isoniazide est
                   dominee par la substitution faux-sens S315T, qui conserve la
                   proteine, et sa perte de fonction complete est minoritaire.
                   Le demander dans les cinq premiers rangs etait une erreur de
                   domaine de ma part, pas un defaut du detecteur ; il reste
                   temoin SECONDAIRE, avec l'attendu « significativement
                   enrichi » et non « en tete ».
                CRITERE GO / NO-GO DE P7.4. P7 continue vers P7.2 si et seulement
                   si au moins un gene 3R passe Q2 ET Q3. Sinon P7 se referme, le
                   mecanisme reste une hypothese de discussion et non un resultat,
                   et il est explicitement interdit de le sauver en revenant a la
                   comparaison naive porteurs / non-porteurs.

              LIMITE DECLAREE D'AVANCE, ET ELLE EST STRUCTURELLE. Seuls les
              `spdi.txt` sont sur le disque ; les `report.json` n'y sont pas
              (correction d'enonce de P3.3), et leur champ `percent_missing` est
              de toute facon inutilisable, etant rabattu sur un catalogue fixe et
              non mesure par souche (KB tuberculosis.md, 2026-08-26). Une souche
              qui a DELETE un gene de reparation ne peut donc pas etre vue ici :
              toute chaine referencee lui impute silencieusement l'allele de
              H37Rv, et elle est comptee comme sauvage. Le recensement est de ce
              fait CONSERVATEUR, et il l'est le plus la ou une deletion
              recurrente est justement documentee, le quartier ogt (Rv1316c) /
              alkA (Rv1317c) / nucS (Rv1321), ou 19 evenements independants ont
              ete decrits. FRONTIERE : cette deletion est l'objet du projet voisin
              `../nucs_deletion_mutators/`, qui travaille a la maille du clade
              porteur d'une lesion ; P7.1 travaille a la maille de la LIGNEE et ne
              compte que les inactivations PONCTUELLES. Un gene qui ne qualifierait
              que par ce bloc delete releve du voisin, pas d'ici.

Entrees     : bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt (142 706 souches)
              investigate_phylo/resources/NC_000962.3.gff3
              H37Rv NC_000962.3.fasta
              résultats/phase7_p81_spectre_96.tsv (spectre mutationnel, A29)
Sorties     : résultats/phase8_p71_genes_3r.tsv        (jeu de genes retenu)
              résultats/phase8_p71_evenements.tsv      (un evenement par ligne)
              résultats/phase8_p71_par_gene.tsv        (Q2, Q3, verdict par gene)
              résultats/phase8_p71_verdict.tsv         (Q1 et go/no-go de P7.4)
Reutilisable: oui -- le balayage filtre par ensemble de couples (position, alt)
              creant un stop, qui ramene 142 706 genomes a un fichier de quelques
              centaines de milliers de lignes, vaut pour toute recherche
              d'inactivation convergente dans une bacterie clonale
Projet      : GC_par_lignee
Date        : 2026-08-30
"""
import argparse
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase2_polarisation_mtbc0 import read_fasta, H37RV  # noqa: E402
from phase4_p93_force_maintien import CODE, GFF3, revcomp  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
MTBC = Path("/home/christophe/docs/codes/mtbc")
BDD = MTBC / "bdd" / "actuelle"
TEMOINS_PRIMAIRES = {"Rv2043c": "pncA", "Rv3854c": "ethA", "Rv3919c": "gid"}
TEMOINS_SECONDAIRES = {"Rv1908c": "katG"}
TEMOINS = {**TEMOINS_PRIMAIRES, **TEMOINS_SECONDAIRES}    # Q1
# Q4 : exclus du jeu 3R parce que cibles ou determinants de resistance connus.
RESISTANCE = {"Rv0006": "gyrA", "Rv0005": "gyrB", "Rv1908c": "katG",
              "Rv2043c": "pncA", "Rv0667": "rpoB", "Rv0668": "rpoC",
              "Rv3795": "embB", "Rv3794": "embA", "Rv3793": "embC",
              "Rv1483": "fabG1", "Rv1484": "inhA", "Rv2245": "kasA",
              "Rv2764c": "thyA", "Rv3919c": "gid", "Rv0682": "rpsL",
              "Rv1694": "tlyA", "Rv2416c": "eis", "Rv0678": "mmpR5",
              "Rv3547": "ddn", "Rv1173": "fbiC", "Rv3261": "fbiA",
              "Rv3262": "fbiB", "Rv2983": "fbiD", "Rv0407": "fgd1",
              "Rv2535c": "pepQ", "Rv1305": "atpE", "Rv3854c": "ethA",
              "Rv3855": "ethR", "Rv1267c": "embR", "Rv2447c": "folC",
              "Rv3609c": "ribD", "Rv0341": "iniB", "Rv1592c": "", }
NOMS_3R = """ung udgB mutY mutT1 mutT2 mutT3 mutT4 ogt alkA nth xthA tagA nei fpg
uvrA uvrB uvrC uvrD1 uvrD2 mfd recA recB recC recD recF recG recN recO recR ruvA
ruvB ruvC radA xerC xerD ssb recX mku ligA ligB ligC ligD dnaE1 dnaE2 dinX dinP
polA lexA dnaQ dnaN alkB priA dnaB helY topA""".split()
LOCUS_3R = ["Rv1321", "Rv0937c", "Rv3202c", "Rv3201c", "Rv1537", "Rv3056",
            "Rv3394c", "Rv3395c", "Rv2924c", "Rv0944", "Rv2464c", "Rv3297",
            "Rv1316c", "Rv1317c", "Rv1210", "Rv1259", "Rv2985", "Rv1160",
            "Rv0413", "Rv0427c", "Rv3674c", "Rv3589", "Rv2736c", "Rv2976c",
            "Rv3908", "Rv3370c"]
KW = re.compile(r"DNA (repair|polymerase|ligase|helicase|glycosylase)|excinuclease|"
                r"recombinase|Holliday|mismatch|excision|RecA|endonuclease (III|V|VIII|NucS)|"
                r"exodeoxyribonuclease|8-oxo-dGTP|methylated-DNA|non-homologous end", re.I)


def lire_cds():
    genes = {}
    for line in GFF3.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        a = f[8]
        lt = re.search(r"locus_tag=([^;]+)", a)
        if not lt or (int(f[4]) - int(f[3]) + 1) % 3:
            continue
        nm = re.search(r"gene=([^;]+)", a)
        pr = re.search(r"product=([^;]+)", a)
        genes[lt.group(1)] = dict(
            locus=lt.group(1), nom=nm.group(1) if nm else "",
            debut=int(f[3]) - 1, fin=int(f[4]), brin=f[6],
            produit=pr.group(1) if pr else "")
    return genes


def jeu_de_genes(genes, n_temoins_neutres=200, seed=0):
    """3R (noms curates + locus verifies + mots-cles de produit), moins les genes
    de resistance (Q4) ; plus katG et pncA comme temoins positifs (Q1) ; plus un
    jeu temoin NEUTRE tire au hasard hors 3R et hors resistance (calibration Q3)."""
    trois_r = {lt for lt, d in genes.items()
               if d["nom"] in NOMS_3R or lt in LOCUS_3R or KW.search(d["produit"])}
    trois_r -= set(RESISTANCE)
    rng = np.random.default_rng(seed)
    dispo = sorted(set(genes) - trois_r - set(RESISTANCE) - set(TEMOINS))
    neutres = set(rng.choice(dispo, size=min(n_temoins_neutres, len(dispo)),
                             replace=False))
    out = []
    for lt in sorted(trois_r | set(TEMOINS) | neutres):
        d = dict(genes[lt])
        d["categorie"] = ("temoin_positif" if lt in TEMOINS else
                          "3R" if lt in trois_r else "temoin_neutre")
        d["longueur"] = d["fin"] - d["debut"]
        out.append(d)
    return pd.DataFrame(out)


def taux_par_canal():
    """Taux mutationnel par site et par canal (classe x contexte trinucleotidique),
    agrege sur les dix-sept lignees du spectre de A29. C'est la ponderation que
    l'en-tete specifie et que la premiere passe n'appliquait pas : un site dont le
    stop demande un C>T n'a pas la meme probabilite qu'un site qui demande un T>A."""
    sp = pd.read_csv(ROOT / "résultats" / "phase7_p81_spectre_96.tsv", sep="\t")
    n = sp.groupby("canal")["n"].sum()
    opp = pd.read_csv(ROOT / "résultats" / "phase7_p81_opportunite_trinuc.tsv",
                      sep="\t").set_index("trinuc")["sites"]
    out = {}
    for canal, cnt in n.items():
        k, ctx = canal.split("@")
        o = opp.get(f"{k[0]}@{ctx}", 0)
        out[canal] = cnt / o if o else 0.0
    return out


def canal_de(seq, pos, alt):
    """Canal (classe pyrimidine x contexte) d'une substitution, meme convention
    de repliement de brin que P8.1."""
    ref = seq[pos]
    if pos < 1 or pos + 1 >= len(seq):
        return None
    l5, r3 = seq[pos - 1], seq[pos + 1]
    if any(b not in "ACGT" for b in (ref, alt, l5, r3)):
        return None
    C = {"A": "T", "C": "G", "G": "C", "T": "A"}
    if ref in "CT":
        return f"{ref}>{alt}@{l5}.{r3}"
    return f"{C[ref]}>{C[alt]}@{C[r3]}.{C[l5]}"


def stops_et_fenetres(seq, jeu, taux):
    """Ensemble des couples (position 0-based, base alternative) qui creent un
    codon stop, et ensemble des positions du jeu de genes. C'est ce filtre, calcule
    une fois sur H37Rv, qui ramene 142 706 genomes a quelques centaines de milliers
    de lignes : sans lui il faudrait lire 6,3 Go de texte pour les jeter."""
    stops, fenetres, opp, nsites = set(), set(), {}, {}
    for g in jeu.itertuples():
        n_stop, k_sites = 0.0, 0
        for pos in range(g.debut, g.fin):
            fenetres.add(pos)
            o = (pos - g.debut) % 3 if g.brin == "+" else (g.fin - 1 - pos) % 3
            start = pos - o if g.brin == "+" else pos - (2 - o)
            if start < g.debut or start + 3 > g.fin:
                continue
            cod = seq[start:start + 3]
            for alt in "ACGT":
                if alt == seq[pos]:
                    continue
                mut = cod[:pos - start] + alt + cod[pos - start + 1:]
                c1, c2 = (cod, mut) if g.brin == "+" else (revcomp(cod), revcomp(mut))
                a1, a2 = CODE.get(c1), CODE.get(c2)
                if a1 is None or a2 is None or a1 == "*":
                    continue
                if a2 == "*":
                    stops.add((pos, alt))
                    k_sites += 1
                    c = canal_de(seq, pos, alt)
                    n_stop += taux.get(c, 0.0) if c else 0.0
        opp[g.locus] = n_stop        # opportunite PONDEREE par le spectre
        nsites[g.locus] = k_sites    # nombre brut de sites, borne du binomial
    return stops, fenetres, opp, nsites


def balayer(stops, fenetres, cache, jobs=1):
    """Un seul passage awk sur les 142 706 spdi.txt. Ne sort que les variants
    creant un stop et les indels tombant dans le jeu de genes."""
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    tmp = cache.parent / "_p71_filtres.txt"
    with open(tmp, "w") as fh:
        for p, a in sorted(stops):
            fh.write(f"S{p}:{a}\n")
        for p in sorted(fenetres):
            fh.write(f"W{p}\n")
    # Le programme awk est ecrit dans un FICHIER et charge par `awk -f` :
    # le passer en argument de `bash -c` obligeait a le faire survivre a deux
    # niveaux de quoting (Python puis shell), ou ses retours a la ligne se
    # transforment en `\n` litteraux et cassent la syntaxe awk.
    # Le programme awk est ecrit dans un FICHIER et charge par `awk -f` : le
    # passer en argument de `sh -c` obligerait a le faire survivre a deux
    # niveaux de quoting (Python puis shell), ou ses retours a la ligne
    # deviennent des `\n` litteraux et cassent la syntaxe awk.
    # DEUX OPTIMISATIONS, validees a sortie strictement identique sur 400
    # genomes (2 817 lignes des deux cotes), qui divisent le temps par quatre :
    #   - `split(FILENAME, ...)` ne s'execute plus a chaque ligne mais une fois
    #     par fichier (FNR==1), soit 142 706 fois au lieu de ~285 millions ;
    #   - la fenetre `$2 in w` prefiltre AVANT la concatenation `$2":"$4`, qui
    #     ne se fait donc plus que sur les 7 % de lignes qui tombent dans le
    #     jeu de genes.
    progf = cache.parent / "_p71_prog.awk"
    progf.write_text(
        'BEGIN { FS=":"; while ((getline l < F) > 0) {\n'
        '          if (substr(l,1,1)=="S") s[substr(l,2)]=1; else w[substr(l,2)]=1 } }\n'
        'FNR==1 { n = split(FILENAME, p, "/"); pre = p[n-3] "\\t" p[n-2] }\n'
        '($2 in w) {\n'
        '  if (length($3)==1 && length($4)==1) {\n'
        '      if (($2 ":" $4) in s) print "stop\\t" pre "\\t" $2 "\\t" $3 "\\t" $4 }\n'
        '  else print "indel\\t" pre "\\t" $2 "\\t" $3 "\\t" $4 }\n')
    # UN FICHIER DE SORTIE PAR LOT. Douze processus awk ecrivant dans le meme
    # tube entrelacent leurs lignes des qu'une ecriture depasse PIPE_BUF : le
    # premier essai a produit des lignes a onze champs et perdu 31 minutes de
    # balayage. Et des lots DIX FOIS plus gros : le cout dominant n'etait pas la
    # lecture des genomes mais le rechargement du filtre de 340 014 lignes a
    # chaque invocation, soit 357 fois au lieu de 36.
    sortie = cache.parent / "_p71_lots"
    subprocess.run(["rm", "-rf", str(sortie)], check=False)
    sortie.mkdir(parents=True)
    cmd = (f"find {BDD} -name spdi.txt -path '*NC_000962.3*' -print0 | "
           f"xargs -0 -n 4000 -P {jobs} sh -c "
           f"'awk -v F={tmp} -f {progf} \"$@\" > {sortie}/lot.$$.tsv' sh")
    print(f"# balayage des spdi.txt (awk, {jobs} processus, lots de 4 000)...",
          file=sys.stderr)
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr[-2000:])
    lignes, mauvaises = [], 0
    for f in sorted(sortie.glob("lot.*.tsv")):
        for l in f.read_text().splitlines():
            c = l.split("\t")
            if len(c) == 6:
                lignes.append(c)
            else:
                mauvaises += 1
    if mauvaises:
        raise RuntimeError(f"{mauvaises} lignes malformees : entrelacement non "
                           f"resolu, ne pas poursuivre sur des donnees douteuses")
    df = pd.DataFrame(lignes, columns=["type", "clade", "sra", "pos", "ref", "alt"])
    df["pos"] = df["pos"].astype(int)
    df.to_csv(cache, sep="\t", index=False)
    tmp.unlink()
    progf.unlink()
    subprocess.run(["rm", "-rf", str(sortie)], check=False)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--min-porteurs", type=int, default=2)   # Q2
    ap.add_argument("--min-evenements", type=int, default=3)  # Q2
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    seq = read_fasta(H37RV)
    genes = lire_cds()
    jeu = jeu_de_genes(genes, seed=args.seed)
    jeu.to_csv(ROOT / "résultats" / "phase8_p71_genes_3r.tsv", sep="\t", index=False)
    print(f"=== jeu de genes : {(jeu.categorie=='3R').sum()} genes 3R, "
          f"{(jeu.categorie=='temoin_positif').sum()} temoins positifs, "
          f"{(jeu.categorie=='temoin_neutre').sum()} temoins neutres ===")

    taux = taux_par_canal()
    stops, fenetres, opp, nsites = stops_et_fenetres(seq, jeu, taux)
    print(f"  opportunite PONDEREE par le spectre a 96 canaux de A29")
    print(f"  {len(stops):,} couples (position, alt) creant un stop ; "
          f"{len(fenetres):,} positions balayees "
          f"({100*len(fenetres)/len(seq):.1f} % du genome)")

    df = balayer(stops, fenetres, ROOT / "data" / "p71_variants_3r.tsv.gz",
                 args.jobs)
    print(f"  {len(df):,} variants ramenes ({(df.type=='stop').sum():,} non-sens, "
          f"{(df.type=='indel').sum():,} indels) sur "
          f"{df.sra.nunique():,} souches")

    # ---- attribution au gene et classement des indels
    bornes = jeu.sort_values("debut")
    deb = bornes["debut"].to_numpy()
    idx = np.searchsorted(deb, df["pos"].to_numpy(), side="right") - 1
    ok = (idx >= 0) & (df["pos"].to_numpy() < bornes["fin"].to_numpy()[np.clip(idx, 0, None)])
    df = df[ok].copy()
    df["locus"] = bornes["locus"].to_numpy()[idx[ok]]
    df["decalage"] = [(len(a) - len(r)) % 3 != 0 if t == "indel" else False
                      for t, r, a in zip(df.type, df.ref, df.alt)]
    df["inactivant"] = (df.type == "stop") | df.decalage

    # ---- evenements : un SPDI distinct = un evenement mutationnel
    inact = df[df.inactivant]
    ev = (inact.groupby(["locus", "pos", "ref", "alt"])
          .agg(n_porteurs=("sra", "nunique"), n_clades=("clade", "nunique"),
               clades=("clade", lambda s: ",".join(sorted(set(s))[:6])))
          .reset_index())
    ev["type"] = np.where(ev["ref"].str.len().eq(1) & ev["alt"].str.len().eq(1),
                          "non-sens", "decalage")
    ev = ev.merge(jeu[["locus", "nom", "categorie", "longueur"]], on="locus")
    ev.sort_values(["categorie", "n_porteurs"], ascending=[True, False]).to_csv(
        ROOT / "résultats" / "phase8_p71_evenements.tsv", sep="\t", index=False)

    # ---- Q2 recurrence, Q3 exces sur l'attendu neutre
    qual = ev[ev.n_porteurs >= args.min_porteurs]
    par = (qual.groupby("locus").agg(n_evenements=("pos", "size"),
                                     n_porteurs_tot=("n_porteurs", "sum"),
                                     n_clades_max=("n_clades", "max"))
           .reindex(jeu["locus"]).fillna(0).reset_index())
    # Q3 porte sur les NON-SENS SEULS : c'est ce que l'en-tete specifie, et le
    # denominateur (opportunite en codons stop) ne couvre pas les decalages de
    # cadre. Les compter ensemble donnait a pncA 304 evenements pour 45 sites.
    ns = (qual[qual.type == "non-sens"].groupby("locus")["pos"].size()
          .reindex(jeu["locus"]).fillna(0).to_numpy())
    par["n_non_sens"] = ns
    par["n_decalages"] = par["n_evenements"] - par["n_non_sens"]
    par = par.merge(jeu[["locus", "nom", "categorie", "longueur"]], on="locus")
    par["opp_stop"] = par["locus"].map(opp)
    # calibration du taux neutre sur les SEULS temoins neutres
    tn = par[par.categorie == "temoin_neutre"]
    # MODELE BINOMIAL, pas Poisson : 22,5 % des sites-stop des temoins sont deja
    # frappes, donc le compte est BORNE par le nombre de sites disponibles et un
    # modele non borne comprime tout enrichissement. n = sites-stop du gene,
    # p = probabilite par site calibree sur les seuls temoins neutres.
    par["n_sites_stop"] = par["locus"].map(nsites)
    p_site = tn["n_non_sens"].sum() / max(par.loc[tn.index, "n_sites_stop"].sum(), 1)
    # ponderation par le spectre : chaque gene a sa propre exposition relative
    expo = par["opp_stop"] / par["n_sites_stop"].replace(0, np.nan)
    expo = expo / (tn["opp_stop"].sum() / max(par.loc[tn.index, "n_sites_stop"].sum(), 1))
    par["p_site_gene"] = np.clip(p_site * expo.fillna(1.0), 0, 1)
    par["attendu"] = par["n_sites_stop"] * par["p_site_gene"]
    par["p_binom"] = [stats.binom.sf(int(k) - 1, int(n), pg) if k > 0 and n > 0 else 1.0
                      for k, n, pg in zip(par.n_non_sens, par.n_sites_stop,
                                          par.p_site_gene)]
    o_ = np.argsort(par["p_binom"].to_numpy())
    q = np.empty(len(par))
    pv = par["p_binom"].to_numpy()[o_]
    q[o_] = np.minimum.accumulate((pv * len(par) / np.arange(1, len(par) + 1))[::-1])[::-1]
    par["q_bh"] = np.clip(q, 0, 1)
    par["Q2"] = par["n_evenements"] >= args.min_evenements
    par["Q3"] = par["q_bh"] <= 0.05
    par = par.sort_values(["q_bh", "n_non_sens"], ascending=[True, False])
    par.to_csv(ROOT / "résultats" / "phase8_p71_par_gene.tsv", sep="\t", index=False)

    print(f"\n  modele binomial calibre sur {len(tn)} temoins neutres : "
          f"p = {p_site:.4f} non-sens par site-stop disponible "
          f"({100*p_site:.1f} % des sites deja frappes -> saturation, "
          f"d'ou le binomial et non Poisson)")
    print("\n=== Q1. CONTROLES POSITIFS ===")
    tp = par[par.categorie == "temoin_positif"]
    rangs = {r.locus: i + 1 for i, r in enumerate(par.itertuples())}
    print(tp[["locus", "nom", "n_non_sens", "n_decalages", "n_porteurs_tot",
              "n_sites_stop", "attendu", "q_bh"]].round(4).to_string(index=False))
    decile = max(1, len(par) // 10)
    prim = tp[tp.locus.isin(TEMOINS_PRIMAIRES)]
    sec = tp[tp.locus.isin(TEMOINS_SECONDAIRES)]
    q1 = bool(len(prim) and (prim.q_bh <= 0.05).all() and
              all(rangs[r.locus] <= decile for r in prim.itertuples()))
    print(f"  temoins PRIMAIRES (perte de fonction = mecanisme de resistance) : "
          f"q_BH <= 0,05 et rang <= {decile} (premier decile)")
    for r in prim.itertuples():
        print(f"    {r.nom:5} rang {rangs[r.locus]:>3}/{len(par)}, "
              f"q_BH = {r.q_bh:.3g} -> "
              f"{'OK' if r.q_bh <= 0.05 and rangs[r.locus] <= decile else 'ECHEC'}")
    for r in sec.itertuples():
        print(f"    {r.nom:5} (SECONDAIRE, attendu seulement enrichi) rang "
              f"{rangs[r.locus]:>3}/{len(par)}, q_BH = {r.q_bh:.3g}")
    print(f"  Q1 {'PASSE' if q1 else 'ECHOUE'}")

    print("\n=== Q2 + Q3. genes 3R candidats ===")
    tr = par[par.categorie == "3R"]
    tete = tr.nsmallest(12, "q_bh")
    print(tete[["locus", "nom", "longueur", "n_non_sens", "n_decalages",
                "n_porteurs_tot", "n_clades_max", "n_sites_stop", "attendu",
                "q_bh", "Q2", "Q3"]].round(4).to_string(index=False))
    passe = tr[tr.Q2 & tr.Q3]
    verdict = ("DETECTEUR CASSE, aucun verdict rendu" if not q1 else
               "GO vers P7.2" if len(passe) else
               "NO-GO : P7 se referme (critere P7.4)")
    print(f"\n=== CRITERE GO / NO-GO DE P7.4 ===")
    print(f"  genes 3R passant Q2 (>= {args.min_evenements} evenements a "
          f">= {args.min_porteurs} porteurs) ET Q3 (q_BH <= 0,05) : {len(passe)}")
    if len(passe):
        print(passe[["locus", "nom", "n_non_sens", "n_decalages",
                     "n_porteurs_tot", "n_clades_max", "attendu", "q_bh"]]
              .round(4).to_string(index=False))
    print(f"  VERDICT P7.1 : {verdict}")
    pd.DataFrame([dict(Q1_controles_positifs=q1, n_3r_Q2=int(tr.Q2.sum()),
                       n_3r_Q3=int(tr.Q3.sum()), n_3r_Q2_et_Q3=len(passe),
                       taux_neutre=taux, verdict=verdict)]).to_csv(
        ROOT / "résultats" / "phase8_p71_verdict.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
