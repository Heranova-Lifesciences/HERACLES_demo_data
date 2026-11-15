#!/usr/bin/env python3
"""
FigS21: cross-parent clusters in the 0.94 clustering (R3 Major 1-3).

A. parents-per-cluster distribution over all 16,931 clusters (log y).
B. The 679 cross-parent clusters split by mechanism:
   isodecoder / same anticodon different aa / cross-anticodon,
   cluster counts (left axis) and representative-read weight (right axis).
C. Effect of blocking cross-parent edges: 16,931 -> 17,559 clusters (+628),
   read from a real re-clustering with the --block-cross-parent flag
   (cluster_summary_blocked.csv). The earlier post-hoc parent split gave
   17,728 and the first external approximation 17,784; both differ from the
   flag because they approximated the edge rules outside the pipeline.

Data: TableS14_cluster_members.tsv / TableS18_cross_parent_clusters.csv.

Output: figures/FigS21_cross_parent.pdf
"""
import csv, os
from collections import defaultdict, Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MEMBERS = os.path.join(BASE, "output", "r1text20261020", "tables", "TableS14_cluster_members.tsv")
OUT = os.path.join(BASE, "figures", "FigS21_cross_parent.pdf")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "pdf.fonttype": 42,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
C_MAIN = "#4a7fb5"
C_ORANGE = "#c08552"
C_TEAL = "#5a9e8f"
C_RED = "#b5544a"
C_EDGE = "#7a7f84"
C_DARK = "#3f3f3f"

CROSS = "#e8a33d"      # cross-parent
SINGLE = "#b9c4cf"     # single-parent


def load():
    c2parents = defaultdict(set)
    c2n = {}
    c2repreads = {}
    with open(MEMBERS, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cid = r["cluster_id"]
            c2parents[cid].add(r["parent_tRNA"])
            c2n[cid] = int(r["cluster_size"])
            if r["is_representative"] == "yes":
                c2repreads[cid] = int(r["read_count_total"])
    return c2parents, c2n, c2repreads


def main():
    c2parents, c2n, c2repreads = load()
    total = len(c2parents)
    dist = Counter(len(p) for p in c2parents.values())
    cross = {c: p for c, p in c2parents.items() if len(p) > 1}

    # mechanism categories
    def fam(p):
        m = p.split("-")
        return m[1], m[2]          # aa, anticodon
    cats = {}
    reads_cats = defaultdict(int)
    for c, ps in cross.items():
        fams = {fam(p) for p in ps}
        if len(fams) == 1:
            k = "isodecoder\nfamily"
        elif len({f[1] for f in fams}) == 1:
            k = "same anticodon\ndiff. aa"
        else:
            k = "cross-anticodon"
        cats[c] = k
        reads_cats[k] += c2repreads[c]
    cat_counts = Counter(cats.values())
    total_reads = sum(c2repreads.values())
    cross_reads = sum(c2repreads[c] for c in cross)

    # block view
    mem_groups = defaultdict(set)
    with open(MEMBERS, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["cluster_id"] in cross:
                mem_groups[r["cluster_id"]].add(r["parent_tRNA"])
    posthoc = total - len(cross) + sum(len(g) for g in mem_groups.values())
    n_block = posthoc
    blocked_csv = os.path.join(BASE, "cluster_summary_blocked.csv")
    if os.path.exists(blocked_csv):
        import pandas as _pd
        _b = _pd.read_csv(blocked_csv, usecols=[0], low_memory=False)
        n_block = int(_b.iloc[:, 0].nunique())
        print(f"blocked clusters from re-clustering: {n_block:,} (post-hoc split would give {posthoc:,})")

    fig, axes = plt.subplots(1, 3, figsize=(8.0, 2.7), constrained_layout=True)

    # ---- A: parents per cluster -------------------------------------
    ax = axes[0]
    ks = [1, 2, 3, 4, 5]
    vals = [dist[k] for k in ks]
    cols = [SINGLE] + [CROSS] * 4
    ax.bar(range(len(ks)), vals, color=cols, edgecolor=C_EDGE, lw=0.5)
    ax.set_yscale("log")
    ax.set_xticks(range(len(ks)))
    ax.set_xticklabels([str(k) for k in ks])
    ax.set_xlabel("Parental tRNAs per cluster")
    ax.set_ylabel("Clusters (log scale)")
    for i, v in enumerate(vals):
        ax.text(i, v * 1.7, f"{v:,}", ha="center", fontsize=6)
    ax.set_ylim(1, 3e5)
    ax.text(0.04, 0.98,
            f"{len(cross):,} / {total:,} clusters\n({len(cross)/total*100:.1f}%) contain >1 parent",
            transform=ax.transAxes, va="top", fontsize=6.2)
    ax.text(0.04, 0.74,
            f"rep. reads in\ncross-parent: {cross_reads/total_reads*100:.1f}%",
            transform=ax.transAxes, va="top", fontsize=6.2, color=C_DARK)
    ax.set_title("A", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    # ---- B: mechanism categories ------------------------------------
    ax = axes[1]
    order = ["isodecoder\nfamily", "same anticodon\ndiff. aa", "cross-anticodon"]
    n_vals = [cat_counts[k] for k in order]
    r_vals = [reads_cats[k] for k in order]
    colors = [C_TEAL, C_MAIN, C_RED]
    ax.bar(range(3), n_vals, color=colors, edgecolor=C_EDGE, lw=0.5)
    ax.set_xticks(range(3))
    ax.set_xticklabels(order, fontsize=6.0)
    ax.set_ylabel("Cross-parent clusters")
    for i, (n, rr) in enumerate(zip(n_vals, r_vals)):
        ax.text(i, n + 14, f"{n}", ha="center", fontsize=6.4)
        if n > 100:
            ax.text(i, n / 2, f"{rr/total_reads*100:.1f}%\nof all reads",
                    ha="center", fontsize=6.2, color="white")
        else:
            ax.text(i, n + 62, f"{rr/total_reads*100:.2f}% of reads",
                    ha="center", fontsize=5.8, color=C_DARK)
    ax.set_ylim(0, max(n_vals) * 1.30)
    ax.set_xlim(-0.55, 2.55)
    ax.set_title("B", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    # ---- C: blocking cross-parent edges -----------------------------
    ax = axes[2]
    ax.bar([0], [total], color=SINGLE, edgecolor=C_EDGE, lw=0.5, width=0.55)
    ax.bar([1], [n_block], color=CROSS, edgecolor=C_EDGE, lw=0.5, width=0.55)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["default\n(cross-parent allowed)", "cross-parent\nblocked (re-clustered)"], fontsize=6.5)
    ax.set_ylabel("Clusters")
    ax.text(0, total + 250, f"{total:,}", ha="center", fontsize=6.5)
    ax.text(1, n_block + 250, f"{n_block:,}", ha="center", fontsize=6.5)
    ax.annotate(f"+{n_block-total:,}", xy=(0.5, (total + n_block) / 2), ha="center",
                fontsize=7, color=C_DARK)
    ax.set_ylim(0, n_block * 1.12)
    ax.set_title("C", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    fig.savefig(OUT)
    # PNG preview alongside the PDF
    fig.savefig(OUT.replace(".pdf", ".png"), dpi=300)
    print(f"saved -> {OUT}")
    print(f"check: cross {len(cross)}  reads% {cross_reads/total_reads*100:.2f}  "
          f"cats {dict(cat_counts)}  block {total} -> {n_block}")


if __name__ == "__main__":
    main()
