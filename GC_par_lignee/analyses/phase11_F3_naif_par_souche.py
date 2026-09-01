#!/usr/bin/env python3
"""
Objet       : F3 (panneau gauche) -- reproduire la quantite naive du sondage
              initial (A7) : pour CHAQUE SOUCHE, le nombre total de substitutions
              simples vs H37Rv (= distance a H37Rv, incluant les substitutions
              partagees par toute la lignee, donc pseudo-repliquees) et le
              rapport pertes/gains de paires G:C compte sur ces memes
              substitutions. C'est le design pre-correction du projet, avant
              polarisation (phase2) et avant passage a l'evenement (phase1
              mode B/C) : il ne doit PAS etre confondu avec phase3_counts_par_
              souche_n40.tsv (deja polarise, branche terminale propre).
Entrees     : bdd/actuelle/<clade>/<SRA>/NC_000962.3/spdi.txt (via phase1_events_
              vs_strains.strain_dirs / read_subs / flux, meme seed et meme n que
              le panel standard du projet)
Sorties     : résultats/phase11_F3_naif_par_souche.tsv (clade, sra, nsub, loss,
              gain, neutral, ratio)
Reutilisable: non -- specifique a la figure F3 (le geste generique vit deja dans
              phase1_events_vs_strains.py)
Projet      : GC_par_lignee
Date        : 2026-08-31
"""
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from phase1_events_vs_strains import read_subs, strain_dirs, flux  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CLADES = ["L1", "L2.2.1", "L3", "L4.1.2", "L4.3", "L5.1", "L6.1", "L7", "L8",
          "L9", "L10", "Bovis.1.1", "Orygis_La3", "Caprae_La2", "Microti"]
N_PER_CLADE = 40
SEED = 0


def main() -> int:
    out = ROOT / "résultats" / "phase11_F3_naif_par_souche.tsv"
    rows = []
    for clade in CLADES:
        strains = strain_dirs(clade)
        if not strains:
            print(f"# {clade} : ABSENT", file=sys.stderr)
            continue
        rng = random.Random(SEED)
        rng.shuffle(strains)
        sample = strains[:N_PER_CLADE]
        for s in sample:
            f = s / "NC_000962.3" / "spdi.txt"
            if not f.exists():
                continue
            subs = read_subs(f)
            if not subs:
                continue
            loss = sum(1 for _, r, a in subs if flux(r, a) == "loss")
            gain = sum(1 for _, r, a in subs if flux(r, a) == "gain")
            neutral = len(subs) - loss - gain
            ratio = loss / gain if gain else float("nan")
            rows.append((clade, s.name, len(subs), loss, gain, neutral, ratio))
        print(f"# {clade} : {len(sample)} souches echantillonnees", file=sys.stderr)

    with open(out, "w") as fh:
        fh.write("clade\tsra\tnsub\tloss\tgain\tneutral\tratio\n")
        for clade, sra, nsub, loss, gain, neutral, ratio in rows:
            fh.write(f"{clade}\t{sra}\t{nsub}\t{loss}\t{gain}\t{neutral}\t{ratio:.4f}\n")
    print(f"Ecrit : {out} ({len(rows)} souches)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
