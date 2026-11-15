# -*- coding: utf-8 -*-
"""Journal-size versions of the two main-text panels that replace the withdrawn
"clustering enhances statistical power" claims.

Fig2C_paired_pvalues   <- what was letter Figure R3-3 (stratified paired p values)
Fig4C_marker_reanalysis<- what was letter Figure R3-5 (Xiong markers, log2FC + P/BH)

Differences from the letter versions: 17 cm width, no internal A/B letters, larger
fonts, 600 dpi PNG and vector PDF so the panels can be dropped into the composite
Figure 2 / Figure 4 in Illustrator.

Output
    figures/正文主图/Fig2C_paired_pvalues.png/.pdf
    figures/正文主图/Fig4C_marker_reanalysis.png/.pdf
"""
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "output", "R1_回应包", "03_tables")
OUT = os.path.join(ROOT, "figures", "正文主图")
BLUE, ORANGE, GREY, RED = "#31688e", "#e07b39", "#9a9a9a", "#b2182b"

plt.rcParams.update({
    "font.size": 8.6, "font.family": "DejaVu Sans",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def neglog(p):
    p = np.asarray(p, dtype=float)
    floor = np.nanmin(p[p > 0]) / 10 if np.any(p > 0) else 1e-12
    return -np.log10(np.clip(p, floor, None))


def fig2c():
    d = pd.read_csv(os.path.join(TAB, "TableS20_fig2C_pairs_stratified.csv"))
    d["x"] = neglog(d["pvalue_non_clustered_matrix"])
    d["y"] = neglog(d["pvalue_clustered_matrix"])
    n_s = int((d.layer == "singleton").sum())
    n_m = int((d.layer == "multi-member").sum())
    d["layer"] = d["layer"].map({"singleton": f"singleton (n={n_s})",
                                 "multi-member": f"multi-member (n={n_m})"})
    cols = {f"multi-member (n={n_m})": GREY, f"singleton (n={n_s})": BLUE}

    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.85),
                             gridspec_kw=dict(width_ratios=[1.12, 1], wspace=0.34,
                                              left=0.085, right=0.985, top=0.90, bottom=0.185))

    ax = axes[0]
    lim = [0, max(d.x.max(), d.y.max()) * 1.06]
    ax.plot(lim, lim, color="#cccccc", lw=0.8, ls="--", zorder=0)
    ax.text(lim[1] * 0.60, lim[1] * 0.66, "no change", rotation=45, fontsize=6.6,
            color="#999999", ha="center", va="center")
    for lab in sorted(d["layer"].unique()):
        s = d[d["layer"] == lab]
        ax.scatter(s.x, s.y, s=14, color=cols[lab], alpha=0.85, lw=0, label=lab)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("−log$_{10}$ $P$, fragment level")
    ax.set_ylabel("−log$_{10}$ $P$, cluster level")
    ax.legend(fontsize=7.2, frameon=False, loc="lower right",
              handletextpad=0.4, borderpad=0.2)
    ax.tick_params(labelsize=7.6, length=2.6)

    ax = axes[1]
    order = ["singleton", "multi-member"]
    groups = [d.loc[d["layer"].str.startswith(o), "delta_neglog10p"].to_numpy() for o in order]
    for i, (g, o) in enumerate(zip(groups, order)):
        lab = [k for k in cols if k.startswith(o)][0]
        x = np.random.default_rng(7 + i).normal(i, 0.055, len(g))
        ax.scatter(x, g, s=13, color=cols[lab], alpha=0.8, lw=0)
        med = np.median(g)
        ax.plot([i - 0.24, i + 0.24], [med, med], color=RED, lw=1.5, zorder=3)
        ax.text(i, 0.025, f"{(g > 0).mean()*100:.0f}% improve\nmedian {med:+.2f}",
                ha="center", va="bottom", fontsize=7.0, color="#444444",
                transform=ax.get_xaxis_transform())
    ax.axhline(0, color="#bbbbbb", lw=0.8, ls=":")
    ax.set_xticks(range(2))
    ax.set_xticklabels([f"singleton\nn={len(groups[0])}", f"multi-member\nn={len(groups[1])}"],
                       fontsize=7.6)
    ax.set_ylabel("Δ(−log$_{10}$ $P$), cluster − fragment")
    lo, hi = min(np.concatenate(groups)), max(np.concatenate(groups))
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(lo * 1.45 - 0.1, hi * 1.20)
    ax.tick_params(labelsize=7.6, length=2.6)

    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"Fig2C_paired_pvalues.{ext}"), dpi=600)
    print("wrote Fig2C (n=%d: %d singleton, %d multi-member)" % (len(d), n_s, n_m))


# ---------------------------------------------------------------- Figure 4C
MARKERS = ["tiRNA-Val-CAC-2", "tiRNA-Lys-CTT-1-M2", "tiRNA-Lys-TTT-3-M2",
           "tiRNA-Lys-CTT-3", "tiRNA-Val-CAC-1-M3"]
SHORT = {"tiRNA-Val-CAC-2": "Val-CAC-2", "tiRNA-Lys-CTT-1-M2": "Lys-CTT-1-M2",
         "tiRNA-Lys-TTT-3-M2": "Lys-TTT-3-M2", "tiRNA-Lys-CTT-3": "Lys-CTT-3",
         "tiRNA-Val-CAC-1-M3": "Val-CAC-1-M3"}
C_NON, C_CLU, CGREY = "#3b6ea5", "#b5523f", "#8a8f98"


def fig4c():
    t19 = pd.read_csv(os.path.join(TAB, "TableS19_marker_benchmark.csv"))
    rec = {}
    for m in MARKERS:
        sub = t19[t19.reported_marker == m]
        d = {}
        for _, r in sub.iterrows():
            w = str(r["workflow"]).lower()
            key = "clu" if "clust" in w and "non" not in w else "non"
            d[key] = dict(l2fc=r["log2FoldChange"], p=r["pvalue"], fdr=r["FDR_BH"])
        rec[m] = d

    fig = plt.figure(figsize=(6.9, 6.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.08], hspace=0.62,
                          left=0.135, right=0.995, top=0.965, bottom=0.095)

    ax = fig.add_subplot(gs[0])
    ymax = max(max(v["non"]["l2fc"], v["clu"]["l2fc"]) for v in rec.values())
    ymin = min(min(v["non"]["l2fc"], v["clu"]["l2fc"]) for v in rec.values())
    gap = 0.30

    def place(vals):
        """label y positions, pushed apart so nothing overlaps"""
        order = np.argsort(vals)
        out = np.array(vals, dtype=float)
        for a, b in zip(order[:-1], order[1:]):
            if out[b] - out[a] < gap:
                out[b] = out[a] + gap
        return out

    left_vals = [rec[m]["non"]["l2fc"] for m in MARKERS]
    right_vals = [rec[m]["clu"]["l2fc"] for m in MARKERS]
    left_lab, right_lab = place(left_vals), place(right_vals)

    for m, y0, y1, ly0, ly1 in zip(MARKERS, left_vals, right_vals, left_lab, right_lab):
        hot = m == "tiRNA-Lys-CTT-3"
        col = RED if hot else CGREY
        ax.plot([0, 1], [y0, y1], color=col, lw=2.0 if hot else 1.1,
                zorder=4 if hot else 2, alpha=1 if hot else 0.85)
        ax.scatter([0, 1], [y0, y1], s=26, color=col, zorder=5 if hot else 3,
                   edgecolors="white", linewidths=0.5)
        for xpt, ypt, ylab, ha in ((0, y0, ly0, "right"), (1, y1, ly1, "left")):
            ax.annotate(SHORT[m], xy=(xpt, ypt), xytext=(-0.06 if ha == "right" else 1.06, ylab),
                        textcoords="data", ha=ha, va="center", fontsize=7.3,
                        color=RED if hot else "#444444",
                        fontweight="bold" if hot else "normal",
                        arrowprops=dict(arrowstyle="-", color=col, lw=0.5,
                                        shrinkA=0, shrinkB=3) if abs(ylab - ypt) > 0.02 else None)

    ax.axhline(0, color="#999999", lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.set_xlim(-0.62, 1.62)
    ax.set_ylim(ymin - 1.05, ymax + 0.55)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["non-clustering\n(fragment level)", "clustering\n(representative level)"],
                       fontsize=7.8)
    ax.set_ylabel("log$_2$ fold change\n(disease / control)", fontsize=8.2)
    ax.tick_params(axis="y", labelsize=7.6, length=2.6)
    ax.tick_params(axis="x", length=0, pad=6)
    ax.annotate("sign reverses",
                xy=(0.985, rec["tiRNA-Lys-CTT-3"]["clu"]["l2fc"]),
                xytext=(0.62, ymin - 0.72), fontsize=7.2, color=RED, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.8,
                                connectionstyle="arc3,rad=0.25"))

    ax = fig.add_subplot(gs[1])
    x = np.arange(len(MARKERS))
    w = 0.19
    series = (("non", "p", C_NON, None, "non-clustering, $P$"),
              ("non", "fdr", C_NON, "////", "non-clustering, BH"),
              ("clu", "p", C_CLU, None, "clustering, $P$"),
              ("clu", "fdr", C_CLU, "////", "clustering, BH"))
    for i, m in enumerate(MARKERS):
        for j, (key, stat, col, hatch, _lab) in enumerate(series):
            v = -np.log10(max(rec[m][key][stat], 1e-300))
            xp = i + (j - 1.5) * w
            if hatch:
                ax.bar(xp, v, w * 0.94, facecolor="white", edgecolor=col, hatch=hatch,
                       lw=0.7, zorder=4)
            else:
                ax.bar(xp, v, w * 0.94, color=col, zorder=3)
    ax.axhline(-np.log10(0.05), color="#666666", lw=0.8, ls=":")
    ax.text(-0.98, -np.log10(0.05) + 0.10, "$P$ = 0.05", fontsize=6.8, color="#555555",
            bbox=dict(fc="white", ec="none", pad=0.5))
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT[m] for m in MARKERS], fontsize=7.5)
    ax.set_ylabel("−log$_{10}$ $P$  /  −log$_{10}$ BH", fontsize=8.0)
    ax.set_ylim(0, 5.7)
    ax.set_xlim(-1.02, len(MARKERS) - 0.38)
    ax.tick_params(labelsize=7.6, length=2.6)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=C_NON, label="non-clustering, $P$"),
                       Patch(facecolor="white", edgecolor=C_NON, hatch="////",
                             label="non-clustering, BH"),
                       Patch(facecolor=C_CLU, label="clustering, $P$"),
                       Patch(facecolor="white", edgecolor=C_CLU, hatch="////",
                             label="clustering, BH")],
              fontsize=7.0, frameon=False, ncol=4, loc="upper left",
              bbox_to_anchor=(-0.062, 1.10), handlelength=1.1, columnspacing=1.0,
              handletextpad=0.35)

    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"Fig4C_marker_reanalysis.{ext}"), dpi=600)
    print("wrote Fig4C (5 markers, both workflows)")


if __name__ == "__main__":
    fig2c()
    fig4c()
