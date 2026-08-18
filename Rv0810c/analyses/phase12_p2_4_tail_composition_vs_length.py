#!/usr/bin/env python3
"""
P2.4 (suite) — la composition (charge nette, fraction D/E) de la queue reste-t-elle
stable quand sa longueur varie, ou la charge s'accumule-t-elle proportionnellement ?
Reutilise le parsing de phase10 (meme frontiere de colonne).
"""
import json
import statistics

from phase10_p2_4_tail_length import (
    STO,
    REF_ID,
    MODULE_END_RESIDUE,
    parse_stockholm,
    find_module_boundary_column,
)

OUT_JSON = "/home/christophe/docs/codes/mtbc/Rv0810c/résultats/p2_4_tail_composition_vs_length.json"

ACIDIC = set("DE")
BASIC = set("KR")


def spearman(x, y):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        ranks = [0] * len(v)
        for pos, i in enumerate(order):
            ranks[i] = pos + 1
        return ranks

    rx, ry = rank(x), rank(y)
    n = len(x)
    d2 = sum((a - b) ** 2 for a, b in zip(rx, ry))
    return 1 - 6 * d2 / (n * (n**2 - 1))


def main():
    seqs, ac_map = parse_stockholm(STO)
    ref_start_residue = int(REF_ID.split("/")[1].split("-")[0])
    boundary_col = find_module_boundary_column(seqs[REF_ID], MODULE_END_RESIDUE, ref_start_residue)

    rows = []
    for sid, row in seqs.items():
        tail = row[boundary_col + 1 :]
        residues = [c.upper() for c in tail if c not in (".", "-")]
        n = len(residues)
        if n < 5:
            continue
        n_acidic = sum(1 for c in residues if c in ACIDIC)
        n_basic = sum(1 for c in residues if c in BASIC)
        net_charge = n_acidic * -1 + n_basic * 1
        rows.append(
            {
                "id": sid,
                "tail_len": n,
                "n_acidic": n_acidic,
                "n_basic": n_basic,
                "net_charge": net_charge,
                "frac_DE": n_acidic / n,
            }
        )

    lengths = [r["tail_len"] for r in rows]
    net_charges = [r["net_charge"] for r in rows]
    frac_DE = [r["frac_DE"] for r in rows]

    rho_charge_vs_len = spearman(lengths, net_charges)
    rho_fracDE_vs_len = spearman(lengths, frac_DE)

    # H37Rv (Rv0810c) pour reference
    h37rv = next((r for r in rows if r["id"] == REF_ID), None)

    summary = {
        "n_sequences": len(rows),
        "spearman_net_charge_vs_tail_len": rho_charge_vs_len,
        "spearman_fracDE_vs_tail_len": rho_fracDE_vs_len,
        "median_frac_DE": statistics.median(frac_DE),
        "median_net_charge": statistics.median(net_charges),
        "H37Rv_Rv0810c": h37rv,
        "note": (
            "rho proche de 0 sur frac_DE => composition maintenue independamment de la longueur "
            "(la queue s'allonge en ajoutant plus de D/E dans la meme proportion, pas en diluant "
            "la charge) ; rho positif fort sur net_charge => la charge negative TOTALE croit "
            "avec la longueur (mecanique, pas informatif en soi)."
        ),
    }

    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
