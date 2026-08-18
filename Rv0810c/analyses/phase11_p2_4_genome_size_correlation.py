#!/usr/bin/env python3
"""
P2.4 (suite) — correle la longueur mediane de queue acide par ORDRE (phase10)
a une taille de genome representative par ordre (NCBI Datasets, genomes
"reference_only", median de assembly_stats.total_sequence_length).
"""
import json
import statistics
import time
import urllib.parse
import urllib.request

TAIL_JSON = "/home/christophe/docs/codes/mtbc/Rv0810c/résultats/p2_4_tail_length.json"
OUT_JSON = "/home/christophe/docs/codes/mtbc/Rv0810c/résultats/p2_4_genome_size_correlation.json"


def fetch_genome_sizes(order, page_size=20):
    url = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{urllib.parse.quote(order)}/dataset_report"
        f"?filters.reference_only=true&page_size={page_size}"
    )
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            d = json.load(resp)
    except Exception as e:
        print(f"  [WARN] {order}: {e}")
        return []
    sizes = []
    for r in d.get("reports", []):
        L = r.get("assembly_stats", {}).get("total_sequence_length")
        if L:
            sizes.append(int(L))
    return sizes


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
    tail_data = json.load(open(TAIL_JSON))
    order_stats = tail_data["order_stats_sorted_by_median"]

    rows = []
    for order, stats in order_stats.items():
        if stats["n"] < 5:
            continue
        sizes = fetch_genome_sizes(order)
        time.sleep(0.4)
        if not sizes:
            print(f"{order}: AUCUN genome recupere, exclu")
            continue
        median_size = statistics.median(sizes)
        rows.append(
            {
                "order": order,
                "n_seq": stats["n"],
                "median_tail_len": stats["median_tail_len"],
                "n_genomes_ref": len(sizes),
                "median_genome_size_bp": median_size,
            }
        )
        print(f"{order:22s} n_seq={stats['n']:4d}  median_tail={stats['median_tail_len']:5.1f}  "
              f"n_genomes={len(sizes):3d}  median_genome_size={median_size/1e6:.2f} Mb")

    tail_vals = [r["median_tail_len"] for r in rows]
    size_vals = [r["median_genome_size_bp"] for r in rows]
    rho = spearman(tail_vals, size_vals) if len(rows) >= 4 else None

    out = {
        "n_orders": len(rows),
        "spearman_rho_tail_vs_genome_size": rho,
        "rows": sorted(rows, key=lambda r: r["median_genome_size_bp"]),
    }
    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    print()
    print(f"n_orders = {len(rows)}")
    print(f"Spearman rho (median_tail_len vs median_genome_size) = {rho}")


if __name__ == "__main__":
    main()
