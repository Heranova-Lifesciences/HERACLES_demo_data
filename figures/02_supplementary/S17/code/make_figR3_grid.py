# -*- coding: utf-8 -*-
"""Letter figure for R3 Major 3.2: the one-parameter-at-a-time grid.

Panel A  cluster count per setting, coloured by the parameter family, with the
         percentage change against the default annotated.
Panel B  nominally significant features under the true labels (filled circle)
         against the median of ten permuted-label replicates (open circle);
         the grey bar spans the permuted minimum to maximum.

Input  output/R3_回应包/03_tables/TableS31_param_grid_topology.csv
Output figures/FigR3_grid_sensitivity.pdf / .png
"""
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "output", "R3_回应包", "03_tables")
FIG = os.path.join(ROOT, "figures")
BLUE, ORANGE, GREY, RED = "#31688e", "#e07b39", "#9a9a9a", "#b2182b"
TEAL, PURPLE, DARK = "#1f9e89", "#7b52ab", "#333333"

plt.rcParams.update({
    "font.size": 8.0, "font.family": "DejaVu Sans",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
})

ORDER = ["default", "abund1.5", "abund3", "kmer8", "kmer12",
         "smk3", "smk8", "ratio090", "ratio098", "standard"]
LABEL = {
    "default": ("default", DARK),
    "abund1.5": ("abundance \u00d71.5", ORANGE),
    "abund3": ("abundance \u00d73.0", ORANGE),
    "kmer8": ("k-mer 8", BLUE),
    "kmer12": ("k-mer 12", BLUE),
    "smk3": ("shared k-mers 3", TEAL),
    "smk8": ("shared k-mers 8", TEAL),
    "ratio090": ("length ratio 0.90", PURPLE),
    "ratio098": ("length ratio 0.98", PURPLE),
    "standard": ("standard", RED),
}


def main():
    d = pd.read_csv(os.path.join(TAB, "TableS31_param_grid_topology.csv"))
    d["o"] = d.setting.map({s: i for i, s in enumerate(ORDER)})
    d = d.sort_values("o").reset_index(drop=True)
    base = int(d.loc[d.setting == "default", "n_clusters"].iloc[0])
    y = np.arange(len(d))[::-1]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.05), sharey=True,
                             gridspec_kw=dict(width_ratios=[1.0, 1.25], wspace=0.12,
                                              left=0.235, right=0.975, top=0.855, bottom=0.20))

    # ---- A: cluster count -------------------------------------------------
    ax = axes[0]
    for yy, (_, r) in zip(y, d.iterrows()):
        c = LABEL[r.setting][1]
        ax.barh(yy, r.n_clusters, height=0.62, color=c, alpha=0.9, lw=0)
        pct = 100.0 * (r.n_clusters - base) / base
        txt = "\u2013" if r.setting == "default" else ("%+.0f%%" % pct)
        ax.text(r.n_clusters + 700, yy, txt, va="center", fontsize=6.6,
                color="#444444" if r.setting != "default" else "#999999")
    ax.axvline(base, color=GREY, lw=0.8, ls=":")
    ax.set_xlim(0, 43500)
    ax.set_xlabel("clusters (43,554 fragments)")
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[s][0] for s in d.setting], fontsize=7.2)
    for tick, s in zip(ax.get_yticklabels(), d.setting):
        tick.set_color(LABEL[s][1])
    ax.set_title("A  segmentation changes a lot", loc="left", fontsize=8.6)
    ax.tick_params(axis="y", length=0)

    # ---- B: true vs permuted significant counts ---------------------------
    ax = axes[1]
    for yy, (_, r) in zip(y, d.iterrows()):
        ax.plot([r.perm_sig_p_min, r.perm_sig_p_max], [yy, yy], color="#cccccc",
                lw=2.6, solid_capstyle="round", zorder=1)
        ax.plot(r.perm_sig_p_median, yy, marker="o", ms=4.6, mfc="white",
                mec="#666666", mew=1.0, zorder=2)
        c = LABEL[r.setting][1]
        ax.plot(r.n_sig_p_full, yy, marker="o", ms=4.8, color=c, mec="white",
                mew=0.5, zorder=3)
        ax.text(470, yy, "%.2f\u00d7" % r.obs_over_perm, va="center", fontsize=6.6,
                color="#444444")
    ax.set_xlim(0, 490)
    ax.set_xlabel("features with nominal p < 0.05")
    ax.set_title("B  true vs permuted significant counts", loc="left", fontsize=8.6)
    ax.text(0.03, 0.97, "filled = true labels\n\u25cb = median of 10 permutations\n"
                        "bar = permuted range", transform=ax.transAxes,
            ha="left", va="top", fontsize=6.2, color="#666666", linespacing=1.5)

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, "FigR3_grid_sensitivity.%s" % ext), dpi=300)
    print("wrote FigR3_grid_sensitivity  (%d settings)" % len(d))


if __name__ == "__main__":
    main()
