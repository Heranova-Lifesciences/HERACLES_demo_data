# -*- coding: utf-8 -*-
"""R3 Major 2.1 / 2.3 figures that were quoted in the letter but never drawn.

FigR3_2C_stratified
    A  paired p-values, fragment-level vs cluster-level, split by whether the
       representative inherits any members (singleton vs multi-member).
    B  distribution of the gain in -log10(p), with medians and the share of
       features that improve.

FigR3_LOO_stability
    Nominal p of every Xiong/Lu marker under each leave-one-library-out
    refit, with BH-significant cells marked.

Inputs
    output/R1_回应包/03_tables/TableS20_fig2C_pairs_stratified.csv
    output/R1_回应包/03_tables/TableS22_LOO_marker_stability.csv
"""
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "output", "R1_回应包", "03_tables")
FIG = os.path.join(ROOT, "figures")

# Features that the family-pattern extraction picked up but that are not
# published markers themselves. Plotted, they inflate the Xiong panel to six
# rows, so they are dropped here; they stay in Table S13 for the record.
NOT_MARKERS = {"tRNA-Lys-TTT-3-1:1-35"}
# display-only: row labels follow the revised Figure 4 marker names
# (positional prefixes dropped; the cluster-level rows keep the name of the
# published marker whose cluster-level statistic they carry)
DISPLAY = {
    "tRNA-Lys-TTT-3-1:1-33": "tiRNA-Lys-TTT-3-M2",
    "tRNA-Lys-CTT-1-1:1-34": "tiRNA-Lys-CTT-1-M2",
    "tRNA-Lys-CTT-2-1:1-33": "tiRNA-Lys-CTT-3",
    "tRNA-Val-CAC-2-1:1-34": "tiRNA-Val-CAC-2",
    "tRNA-Val-CAC-1-1:1-34": "tiRNA-Val-CAC-1-M3",
    "tRNA-Val-CAC-2-1:1-29": "tRF-Val-CAC-2",
    "tRNA-Arg-ACG-1-1:54-76": "tRF-Arg-ACG-1-M2",
    "tRNA-Ser-CGA-1-1:1-24": "tRF-Ser-CGA-1-M3",
}
BLUE, ORANGE, GREY, RED = "#31688e", "#e07b39", "#9a9a9a", "#b2182b"

plt.rcParams.update({
    "font.size": 8.2, "font.family": "DejaVu Sans",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
})


def neglog(p):
    p = np.asarray(p, dtype=float)
    floor = np.nanmin(p[p > 0]) / 10 if np.any(p > 0) else 1e-12
    return -np.log10(np.clip(p, floor, None))


def fig2c():
    d = pd.read_csv(os.path.join(TAB, "TableS20_fig2C_pairs_stratified.csv"))
    d["x"] = neglog(d["pvalue_non_clustered_matrix"])
    d["y"] = neglog(d["pvalue_clustered_matrix"])
    d["layer"] = d["layer"].map({"singleton": "singleton (n=%d)" % (d.layer == "singleton").sum(),
                                 "multi-member": "multi-member (n=%d)" % (d.layer == "multi-member").sum()})
    cols = {k: c for k, c in zip(sorted(d["layer"].unique()), [GREY, BLUE])}

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9),
                             gridspec_kw=dict(width_ratios=[1.15, 1], wspace=0.33,
                                              left=0.085, right=0.985, top=0.86, bottom=0.19))

    ax = axes[0]
    lim = [0, max(d.x.max(), d.y.max()) * 1.06]
    ax.plot(lim, lim, color="#cccccc", lw=0.8, ls="--", zorder=0)
    ax.text(lim[1] * 0.62, lim[1] * 0.66, "no change", rotation=45, fontsize=6.2,
            color="#999999", ha="center", va="center")
    for lab in sorted(d["layer"].unique()):
        s = d[d["layer"] == lab]
        ax.scatter(s.x, s.y, s=13, color=cols[lab], alpha=0.85, lw=0, label=lab)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("−log$_{10}$ p, fragment level")
    ax.set_ylabel("−log$_{10}$ p, cluster level")
    ax.set_title("A  paired p-values per tested feature", loc="left", fontsize=8.8)
    ax.legend(fontsize=6.6, frameon=False, loc="lower right",
              handletextpad=0.4, borderpad=0.2)

    ax = axes[1]
    order = ["singleton", "multi-member"]
    groups = [d.loc[d["layer"].str.startswith(o), "delta_neglog10p"].to_numpy() for o in order]
    for i, (g, o) in enumerate(zip(groups, order)):
        lab = [k for k in cols if k.startswith(o)][0]
        x = np.random.default_rng(7 + i).normal(i, 0.055, len(g))
        ax.scatter(x, g, s=12, color=cols[lab], alpha=0.8, lw=0)
        med = np.median(g)
        ax.plot([i - 0.22, i + 0.22], [med, med], color=RED, lw=1.4, zorder=3)
        ax.text(i, 0.02, f"{(g > 0).mean()*100:.0f}% improve · median {med:+.2f}",
                ha="center", fontsize=6.2, color="#555555",
                transform=ax.get_xaxis_transform())
    ax.axhline(0, color="#bbbbbb", lw=0.8, ls=":")
    ax.set_xticks(range(2))
    ax.set_xticklabels(["singleton\nn=%d" % len(groups[0]), "multi-member\nn=%d" % len(groups[1])],
                       fontsize=7)
    ax.set_ylabel("Δ(−log$_{10}$ p), cluster − fragment")
    lo = min(np.concatenate(groups))
    hi = max(np.concatenate(groups))
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(lo * 1.35 - 0.05, hi * 1.18)
    ax.set_title("B  gain per feature", loc="left", fontsize=8.8)

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"FigR3_2C_stratified.{ext}"), dpi=300)
    print("wrote FigR3_2C_stratified (n=%d, %d degenerate)" % (
        len(d), int((d["log2FC_clustered"].abs() > 30).sum())))


def loo():
    d = pd.read_csv(os.path.join(TAB, "TableS22_LOO_marker_stability.csv"))
    d = d[d["workflow"].str.startswith("clustering")]
    # Block labels give the group sizes, not the row count: after dropping the
    # non-marker singleton the Xiong block holds five markers, so "(n = 6)"
    # would read as a sixth marker.
    d["ds"] = np.where(d["dataset"].str.startswith("Xiong"),
                       "Xiong 2024 (3 vs 3)", "Lu 2021 (4 vs 4)")
    # Only rows that carry a reported marker. tRNA-Lys-TTT-3-1:1-35 is a
    # separate 35-nt singleton cluster with no published counterpart; it used
    # to be pulled in through the "Lys-TTT-3" family pattern and shown as a
    # sixth Xiong row, which wrongly read as a sixth marker.
    d = d[~d["feature"].isin(NOT_MARKERS)]
    feats = list(dict.fromkeys(d["feature"]))
    blocks = [("Xiong 2024 (3 vs 3)", "drop_GSM"), ("Lu 2021 (4 vs 4)", "drop_SRR")]

    colmap = {}
    xlabels = {}
    col = 0
    for ds, pref in blocks:
        runs = ["full"] + sorted([r for r in d.run.unique() if r.startswith(pref)])
        for k, r in enumerate(runs):
            colmap[(ds, r)] = col
            # sequential library number within each block: full 1 2 3 ...
            xlabels[(ds, r)] = "full" if r == "full" else str(k)
            col += 1
        col += 0.8

    xiong_rows = [f for f in feats
                  if d[(d.feature == f) & (d.ds.str.startswith("Xiong"))].shape[0]]
    lu_rows = [f for f in feats if f not in xiong_rows]
    allrows = xiong_rows + lu_rows

    fig, ax = plt.subplots(figsize=(7.0, 3.1),
                           gridspec_kw=dict(left=0.28, right=0.885, top=0.87, bottom=0.19))
    for i, f in enumerate(allrows):
        sub = d[d.feature == f]
        for _, r in sub.iterrows():
            c = colmap.get((r.ds, r.run))
            if c is None:
                continue
            p = float(r.nominal_p)
            if pd.isna(p):
                ax.add_patch(plt.Rectangle((c - 0.45, i - 0.45), 0.9, 0.9,
                                           color="#eeeeee", lw=0))
                ax.text(c, i, "NA", ha="center", va="center", fontsize=5.0, color="#888888")
                continue
            v = -np.log10(max(p, 1e-4))
            ax.add_patch(plt.Rectangle((c - 0.45, i - 0.45), 0.9, 0.9,
                                       color=plt.cm.Blues(min(v / 3.2, 1.0)), lw=0))
            if pd.notna(r.padj) and r.padj < 0.05:
                ax.plot(c, i, marker="o", ms=3.4, mfc="none", mec=RED, mew=0.9)
            ax.text(c, i, "%.2f" % p, ha="center", va="center",
                    fontsize=4.6, color="white" if v > 2.0 else "#333333")
    ax.set_xlim(-0.6, col - 0.3)
    ax.set_ylim(len(allrows) - 0.4, -0.55)
    ax.set_yticks(range(len(allrows)))
    ax.set_yticklabels([DISPLAY.get(f, f.replace("tRNA-", "").replace("-1:", " ")) for f in allrows],
                       fontsize=5.8)
    ticks = [c for c in colmap.values()]
    labels = [xlabels[(ds, r)] for (ds, r) in colmap.keys()]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=5.2, rotation=90)
    for ds, _ in blocks:
        cs = [c for (d2, _r), c in colmap.items() if d2 == ds]
        ax.text(np.mean(cs), 1.02, ds, ha="center", fontsize=7,
                transform=ax.get_xaxis_transform())
    for s in ("left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    sm = plt.cm.ScalarMappable(cmap="Blues", norm=plt.Normalize(0, 3.2))
    cax = fig.add_axes([0.905, 0.24, 0.012, 0.56])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("−log$_{10}$ nominal p", fontsize=6.4)
    cb.ax.tick_params(labelsize=5.6, width=0.6)
    cb.outline.set_linewidth(0.5)
    ax.set_title("Leave-one-library-out refit of published markers (cluster level); "
                 "open circle = BH < 0.05", loc="left", fontsize=8.2, pad=16)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"FigR3_LOO_stability.{ext}"), dpi=300)
    print("wrote FigR3_LOO_stability (%d features)" % len(allrows))


if __name__ == "__main__":
    fig2c()
    loo()
