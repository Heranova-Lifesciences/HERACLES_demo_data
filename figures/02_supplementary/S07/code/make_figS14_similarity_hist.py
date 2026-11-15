# -*- coding: utf-8 -*-
"""Redraw Figure S14 (Letter Figure R2-1) with "similarity" instead of "identity".

Reviewer comment: the metric is interval overlap divided by interval union, which is
a sequence similarity, not a sequence identity (it needs no alignment and treats
partial overlaps as partial similarity). The numbers are unchanged; only the wording
on the axes changes.

A  pairwise similarity between representatives of different clusters on the same
   parent tRNA (interval overlap / union, no alignment step)
B  own similarity minus similarity to the nearest representative of another cluster

Writes figures/FigS14_representative_identity.{png,pdf} so existing references
keep working, plus copies named *_similarity.* for clarity.
"""
import os
import re
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "figures")
PNG = os.path.join(OUT, "_build_png")

C_A = "#c98a4b"
C_B = "#5b7c99"
GRID = "#e3e3e3"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.6,
    "axes.labelsize": 9.4,
    "xtick.labelsize": 8.4,
    "ytick.labelsize": 8.4,
    "pdf.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04,
})

PAT = re.compile(r"^(tRNA-[^:]+):(\d+)-(\d+):")


def main():
    cs = pd.read_csv(os.path.join(BASE, "r1_diag", "cluster_summary_detailed_46.csv"),
                     low_memory=False)
    frag = cs[cs.is_cluster_total == False].copy()
    frag["is_rep"] = frag.tsRNA_id.astype(str) == frag.rep_tsRNA_id.astype(str)
    size = frag.groupby("cluster_id").size()
    multi = set(size[size >= 2].index)

    rep = frag[frag.is_rep & frag.cluster_id.isin(multi)].copy()

    def parse(s):
        m = PAT.match(str(s))
        return (m.group(1), int(m.group(2)), int(m.group(3))) if m else (None, None, None)

    rep[["parent", "start", "end"]] = rep.tsRNA_id.apply(lambda x: pd.Series(parse(x)))
    rep = rep.dropna(subset=["parent"])

    vals = []
    for _, g in rep.groupby("parent"):
        iv = list(zip(g.start.astype(int), g.end.astype(int)))
        for (a1, a2), (b1, b2) in combinations(iv, 2):
            inter = max(0, min(a2, b2) - max(a1, b1) + 1)
            union = (a2 - a1 + 1) + (b2 - b1 + 1) - inter
            vals.append(inter / union if union > 0 else 0.0)
    v = np.array(vals)

    sep = pd.read_csv(os.path.join(BASE, "r1_diag", "cluster_separation.csv"))
    margin = sep.margin.dropna().to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.55),
                             gridspec_kw={"wspace": 0.22})

    ax = axes[0]
    w = np.ones(len(v)) / len(v) * 100.0            # percent of pairs, not raw counts
    counts, _, patches = ax.hist(v, bins=np.arange(0, 1.02, 0.02), weights=w,
                                 color=C_A, edgecolor="white", linewidth=0.35,
                                 zorder=3)
    # the zero bin pools every non-overlapping pair; colour it apart and label it
    patches[0].set_facecolor("#9a9a9a")
    ax.set_xlabel("Sequence similarity between representatives")
    ax.set_ylabel("Pairs (%)")
    ax.set_title("A  Representatives of different clusters", loc="left",
                 fontsize=9.6, pad=6)
    ax.set_xlim(-0.025, 1.02)
    ax.set_ylim(0, counts.max() * 1.28)
    ax.grid(axis="y", color=GRID, zorder=0)
    ax.set_axisbelow(True)

    med = np.median(v)
    ax.axvline(med, color="#444444", lw=0.9, ls=(0, (4, 3)), zorder=4)
    ax.text(med + 0.02, counts.max() * 1.16, f"median {med:.2f}",
            fontsize=8.2, color="#333333")
    ax.annotate(f"{100 * (v == 0).mean():.1f}%\nno overlap",
                xy=(0.01, counts[0]), xytext=(0.15, counts[0] * 0.92),
                fontsize=8.0, color="#5a5a5a",
                arrowprops=dict(arrowstyle="-", color="#9a9a9a", lw=0.8))
    ax.text(0.985, counts.max() * 1.06, f"{100 * (v > 0.90).mean():.1f}% above 0.90",
            ha="right", fontsize=8.2, color="#333333")

    ax = axes[1]
    ax.hist(margin, bins=np.arange(-0.6, 0.62, 0.02), color=C_B, edgecolor="white",
            linewidth=0.3, zorder=3)
    ax.axvline(0, color="#555555", lw=0.9, ls=(0, (4, 3)), zorder=4)
    ax.set_xlabel("Own minus nearest-other similarity")
    ax.set_ylabel("Members")
    ax.set_title("B  Which representative is closer?", loc="left",
                 fontsize=9.6, pad=6)
    ax.set_xlim(-0.68, 0.68)
    ax.grid(axis="y", color=GRID, zorder=0)
    ax.set_axisbelow(True)
    ax.text(0.03, 0.95, f"{100 * (margin > 0).mean():.0f}% own-closer",
            transform=ax.transAxes, fontsize=8.4)
    ax.text(0.03, 0.87, f"margin {np.median(margin):+.2f}",
            transform=ax.transAxes, fontsize=8.4)

    for a in axes:
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)

    targets = [
        os.path.join(OUT, "FigS14_representative_identity"),
        os.path.join(OUT, "FigS14_representative_similarity"),
    ]
    for stem in targets:
        fig.savefig(stem + ".pdf")
    # only the canonical PNG name is used by the letter
    fig.savefig(os.path.join(PNG, "FigS14_representative_identity.png"), dpi=300)
    fig.savefig(os.path.join(OUT, "FigS14_representative_similarity.png"), dpi=300)
    print(f"pairs {len(v):,} | median {np.median(v):.3f} | >0.90 {100*(v>0.90).mean():.1f}% "
          f"| non-overlap {100*(v==0).mean():.1f}%")
    print(f"margin median {np.median(margin):+.3f} | own-closer {100*(margin>0).mean():.1f}%")
    print("saved", targets[0] + ".{png,pdf} and *_similarity copies")


if __name__ == "__main__":
    main()
