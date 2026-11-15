# -*- coding: utf-8 -*-
"""Letter figure for R3 Major 2.2: FDR-positive fraction against the number of tests.

One point per analysis (three datasets x two aggregation levels).  x is the
number of features that entered the test, y is the percentage of those features
reaching Benjamini-Hochberg < 0.05, and the arrow runs from the fragment level to
the cluster level of the same dataset.

All values are read from Table S25 (observed rows, perm = -1); nothing is typed
in by hand.

Writes figures/FigR3_fdr_vs_tests.{png,pdf}.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(BASE, "output", "补充材料_提交版", "TableS25_permuted_label_controls.csv")
OUT = os.path.join(BASE, "figures")

COL = {"cohort": "#8c5b9e", "paper1": "#5b7c99", "paper2": "#c0743c"}
NAME = {"cohort": "Endometriosis cohort", "paper1": "Xiong 2024", "paper2": "Lu 2021"}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.8,
    "axes.labelsize": 9.4,
    "xtick.labelsize": 8.4,
    "ytick.labelsize": 8.4,
    "pdf.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})


def load():
    d = pd.read_csv(TAB)
    d = d[d.perm == -1].copy()
    frag = d[d.workflow.isin(["frag", "nonclust"])].copy()
    clus = d[d.workflow == "clust"].copy()

    def key(s):
        for k in COL:
            if s.startswith(k):
                return k
        return s

    frag["ds"] = frag.dataset.map(key)
    clus["ds"] = clus.dataset.map(key)
    out = {}
    for k in COL:
        f = frag[frag.ds == k].iloc[0]
        c = clus[clus.ds == k].iloc[0]
        out[k] = dict(
            f_tests=int(f.n_tested), f_sig=int(f.n_sig_fdr),
            c_tests=int(c.n_tested), c_sig=int(c.n_sig_fdr),
        )
    return out


def pct(sig, tests):
    return 100.0 * sig / tests


def main():
    v = load()
    fig, ax = plt.subplots(figsize=(6.6, 4.1))

    NAME_AT = {"cohort": (15500, 1.00), "paper1": (1850, 2.55), "paper2": (2760, 8.55)}
    FRAG_AT = {"cohort": (21288, -0.80), "paper1": (1652, -0.80), "paper2": (3108, 6.15)}
    CLUS_AT = {"cohort": (10345, 0.55), "paper1": (1084, 2.98), "paper2": (2002, 10.55)}

    for k, d in v.items():
        xf, xc = d["f_tests"], d["c_tests"]
        yf = pct(d["f_sig"], xf)
        yc = pct(d["c_sig"], xc)
        col = COL[k]

        ax.annotate("", xy=(xc, yc), xytext=(xf, yf),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.1,
                                    shrinkA=6, shrinkB=8, alpha=0.75,
                                    connectionstyle="arc3,rad=0.06"))
        ax.plot([xf], [yf], "o", ms=6.6, color="white", mec=col, mew=1.5, zorder=5)
        ax.plot([xc], [yc], "s", ms=6.6, color=col, mec=col, zorder=5)

        ax.text(*NAME_AT[k], NAME[k], color=col, fontsize=8.4, fontweight="bold",
                ha="left" if k != "cohort" else "center")
        ax.text(*FRAG_AT[k], f"{d['f_sig']} / {xf:,}", ha="center", fontsize=7.4,
                color="#555555")
        ax.text(*CLUS_AT[k], f"{d['c_sig']} / {xc:,}", ha="center", fontsize=7.4,
                color="#555555")

    ax.set_xscale("log")
    ax.set_xlim(800, 40000)
    ax.set_ylim(-1.6, 11.5)
    ax.set_xlabel("Features tested (log scale)")
    ax.set_ylabel("Features reaching BH < 0.05 (%)")
    ax.grid(axis="both", color="#e8e8e8", lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # marker key, no legend box
    ax.plot([], [], "o", ms=6.2, color="white", mec="#666666", mew=1.4, label="fragment level")
    ax.plot([], [], "s", ms=6.2, color="#666666", label="cluster level")
    ax.legend(loc="upper right", frameon=False, fontsize=8.0, handletextpad=0.5,
              borderaxespad=0.2)

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"FigR3_fdr_vs_tests.{ext}"),
                    dpi=300 if ext == "png" else None)

    for k, d in v.items():
        print(f"{k:7s} fragment {d['f_sig']}/{d['f_tests']:,} = {pct(d['f_sig'], d['f_tests']):.2f}%"
              f" | cluster {d['c_sig']}/{d['c_tests']:,} = {pct(d['c_sig'], d['c_tests']):.2f}%")
    print("saved", os.path.join(OUT, "FigR3_fdr_vs_tests.png"))


if __name__ == "__main__":
    main()
