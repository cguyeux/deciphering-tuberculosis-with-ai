#!/usr/bin/env python3
"""
Objet       : convertir le rapport pertes/gains de paires G:C en la quantite
              que la theorie designe : le rapport des TAUX mutationnels v/u
              (par site disponible) et le GC d'equilibre de Sueoka
              GC_eq = u/(u+v). Le rapport brut confond taux et opportunite --
              un genome a 65 % de GC offre deux fois plus de sites G:C a perdre
              que de sites A:T a gagner. Normaliser par le nombre de sites
              disponibles dans le genome non masque rend la mesure comparable a
              la litterature (Sueoka 1988) et interpretable : chaque lignee
              recoit le GC vers lequel son propre flux mutationnel la pousse.
              Fournit aussi la stratification par GC local, ou la normalisation
              sert de controle interne : le gradient du rapport brut le long du
              GC local doit s'effondrer une fois l'opportunite retiree.
Entrees     : resultats/phase3_counts_par_souche_n40.tsv (par lignee)
              resultats/phase3_gc_local_n40.tsv (par strate de GC local)
              H37Rv NC_000962.3.fasta + data/MTBC0/Mask.files/*.bed
Sorties     : TSV par lignee (v/u, GC_eq, IC95) et TSV par lignee x strate
Reutilisable: oui -- tout projet mesurant un flux de composition sur genome
Projet      : GC_par_lignee
Date        : 2026-08-29
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from phase2_polarisation_mtbc0 import read_fasta, load_mask, H37RV  # noqa: E402
from phase3_gc_local import gc_profile, MASKS  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


def opportunity(window=300, strates=4):
    """Nombre de sites G:C et A:T disponibles, en tout et par strate de GC local,
    sur le genome non masque."""
    seq = read_fasta(H37RV)
    arr = np.frombuffer(seq.encode(), dtype=np.uint8)
    is_gc = (arr == ord("G")) | (arr == ord("C"))
    is_at = (arr == ord("A")) | (arr == ord("T"))
    masked = load_mask(MASKS, in_mtbc0_coords=True)
    keep = np.ones(len(seq), bool)
    keep[np.fromiter(masked, int, len(masked))] = False
    gc = gc_profile(seq, window)
    qs = np.quantile(gc, np.linspace(0, 1, strates + 1))
    qs[0], qs[-1] = -1, 2
    st = np.clip(np.searchsorted(qs, gc, side="right") - 1, 0, strates - 1)
    per = {f"Q{s+1}": (int((is_gc & keep & (st == s)).sum()),
                       int((is_at & keep & (st == s)).sum()))
           for s in range(strates)}
    return (int((is_gc & keep).sum()), int((is_at & keep).sum())), per


def vu_ci(loss, gain, n_gc, n_at, conf=0.95):
    """v/u = (loss/n_gc)/(gain/n_at) et son IC95 par delta-methode sur le log."""
    vu = (loss / n_gc) / (gain / n_at)
    z = stats.norm.ppf(0.5 + conf / 2)
    se = np.sqrt(1 / loss + 1 / gain)
    return vu, vu * np.exp(-z * se), vu * np.exp(z * se)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default=str(ROOT / "résultats" /
                                            "phase3_counts_par_souche_n40.tsv"))
    ap.add_argument("--gc-local", default=str(ROOT / "résultats" /
                                              "phase3_gc_local_n40.tsv"))
    ap.add_argument("--window", type=int, default=300)
    ap.add_argument("--out-prefix", default=str(ROOT / "résultats" / "phase3_"))
    args = ap.parse_args()

    (n_gc, n_at), per = opportunity(args.window)
    print(f"# genome non masque : {n_gc} sites G:C, {n_at} sites A:T "
          f"(GC = {n_gc/(n_gc+n_at):.4f})", file=sys.stderr)

    c = pd.read_csv(args.counts, sep="\t").groupby("clade")[
        ["loss", "gain"]].sum().reset_index()
    rows = []
    for _, r in c.iterrows():
        vu, lo, hi = vu_ci(r["loss"], r["gain"], n_gc, n_at)
        rows.append(dict(clade=r["clade"], loss=int(r["loss"]),
                         gain=int(r["gain"]), ratio_brut=r["loss"] / r["gain"],
                         vu=vu, vu_lo=lo, vu_hi=hi,
                         gc_eq=1 / (1 + vu), gc_eq_lo=1 / (1 + hi),
                         gc_eq_hi=1 / (1 + lo)))
    t = pd.DataFrame(rows).sort_values("vu", ascending=False)
    t.to_csv(args.out_prefix + "J_sueoka_gc_eq.tsv", sep="\t", index=False)
    print(t.to_string(index=False))

    g = pd.read_csv(args.gc_local, sep="\t")
    g["n_gc"] = g["strate_gc"].map({k: v[0] for k, v in per.items()})
    g["n_at"] = g["strate_gc"].map({k: v[1] for k, v in per.items()})
    g["vu"] = (g["loss"] / g["n_gc"]) / (g["gain"] / g["n_at"])
    g.to_csv(args.out_prefix + "K_sueoka_par_gc_local.tsv", sep="\t", index=False)
    tot = g.groupby("strate_gc").apply(
        lambda x: pd.Series(dict(
            ratio_brut=x["loss"].sum() / x["gain"].sum(),
            vu=(x["loss"].sum() / x["n_gc"].iloc[0]) /
               (x["gain"].sum() / x["n_at"].iloc[0]),
            opp_gc_sur_at=x["n_gc"].iloc[0] / x["n_at"].iloc[0])),
        include_groups=False)
    print("\n# controle interne : le gradient le long du GC local est-il un "
          "effet d'opportunite ?")
    print(tot.round(3).to_string())
    print(f"  amplitude Q1->Q4 : rapport brut x"
          f"{tot['ratio_brut'].iloc[-1]/tot['ratio_brut'].iloc[0]:.2f}, "
          f"taux v/u x{tot['vu'].iloc[-1]/tot['vu'].iloc[0]:.2f}")


if __name__ == "__main__":
    main()
