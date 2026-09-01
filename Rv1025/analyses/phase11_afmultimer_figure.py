#!/usr/bin/env python3
"""
P6.3 - Figure AF-Multimer : ipTM (max sur 5 modèles) et PAE inter-chaînes min par job.
Distingue TEST (Rv1025-divIC), contrôle POSITIF (divIC-ftsQ), contrôles de spécificité (eno/ppx2/ftsQ).
Sortie : article/figures/figure3_afmultimer.{png,pdf}.
"""
import glob, json, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/christophe/docs/codes/mtbc/Rv1025"
AF = f"{ROOT}/résultats/af3_out"
OUT = f"{ROOT}/article/figures/figure3_afmultimer"

LABEL = {"divic_ftsq_posctrl": "DivIC–FtsQ\n(positive control)",
         "rv1025_divic_test": "Rv1025–DivIC\n(test)",
         "rv1025_ftsq": "Rv1025–FtsQ",
         "rv1025_eno_ctrl": "Rv1025–Eno",
         "rv1025_ppx2_ctrl": "Rv1025–Ppx2"}
KIND = {"divic_ftsq_posctrl": "pos", "rv1025_divic_test": "test",
        "rv1025_ftsq": "ctrl", "rv1025_eno_ctrl": "ctrl", "rv1025_ppx2_ctrl": "ctrl"}
COL = {"pos": "#55A868", "test": "#C44E52", "ctrl": "#8C8C8C"}

def collect():
    d = {}
    for job in os.listdir(AF):
        p = f"{AF}/{job}"
        if not os.path.isdir(p): continue
        key = job.replace("fold_", "")
        if key not in LABEL: continue
        ipt, pae = [], []
        for js in glob.glob(f"{p}/*summary_confidences*.json"):
            x = json.load(open(js)); ipt.append(x["iptm"])
            m = x["chain_pair_pae_min"]; pae.append(min(m[0][1], m[1][0]))
        d[key] = (max(ipt), min(pae))
    return d

def main():
    d = collect()
    order = sorted(d, key=lambda k: -d[k][0])
    labels = [LABEL[k] for k in order]
    cols = [COL[KIND[k]] for k in order]
    iptm = [d[k][0] for k in order]
    pae = [d[k][1] for k in order]
    x = range(len(order))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))
    a1.bar(x, iptm, color=cols, edgecolor="black", lw=0.6)
    a1.axhline(0.6, ls="--", color="grey", lw=1); a1.text(len(order)-0.4, 0.61, "ipTM 0.6", fontsize=8, color="grey", ha="right")
    a1.set_ylabel("best ipTM (5 models)"); a1.set_ylim(0, 0.7); a1.set_title("(A) Interface confidence", fontsize=11)
    a2.bar(x, pae, color=cols, edgecolor="black", lw=0.6)
    a2.set_ylabel("min inter-chain PAE (Å)"); a2.set_title("(B) Inter-chain error", fontsize=11)
    for a in (a1, a2):
        a.set_xticks(list(x)); a.set_xticklabels(labels, fontsize=8.5)
        a.spines[["top", "right"]].set_visible(False)
    for xi, v in zip(x, iptm): a1.text(xi, v+0.012, f"{v:.2f}", ha="center", fontsize=8.5)
    for xi, v in zip(x, pae): a2.text(xi, v+0.3, f"{v:.1f}", ha="center", fontsize=8.5)
    plt.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT + ".pdf", bbox_inches="tight")
    print("écrit:", OUT + ".png/.pdf ;", {k: (round(d[k][0],2), round(d[k][1],1)) for k in order})

if __name__ == "__main__":
    main()
