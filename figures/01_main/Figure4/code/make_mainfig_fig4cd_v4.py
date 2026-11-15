# -*- coding: utf-8 -*-
"""Figure 4 C/D, both datasets.

C  slope graph of log2 fold change (non-clustering -> clustering) for the five
   Xiong 2024 markers and the three Lu 2021 markers; colour identifies the
   dataset, with a legend.
D  -log10 nominal P and -log10 BH per marker, four bars per marker; colour
   identifies the workflow (blue = non-clustering, brick = clustering) and
   hatching marks BH. Datasets are labelled above their marker groups.

Output: figures/正文主图/Fig4_panelC.*, Fig4_panelD.*,
        figures/正文主图/Fig4C_marker_reanalysis_v4.*
"""
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "output", "R1_回应包", "03_tables")
OUT = os.path.join(ROOT, "figures", "正文主图")

plt.rcParams.update({
    "font.size": 8.6, "font.family": "DejaVu Sans",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

# dataset colours: Xiong blue, Lu warm orange
X_NON, X_CLU = "#2c5985", "#7aa6cc"
L_NON, L_CLU = "#a8452f", "#e39a7f"
# panel D: colour encodes the workflow, so the four legend entries can be
# written out in full (blue = non-clustering, brick = clustering).
D_NON, D_CLU = "#3d6fa8", "#b0502c"

SHORT = {
    "tiRNA-Val-CAC-2": "Val-CAC-2",
    "tiRNA-Lys-CTT-1-M2": "Lys-CTT-1-M2",
    "tiRNA-Lys-TTT-3-M2": "Lys-TTT-3-M2",
    "tiRNA-Lys-CTT-3": "Lys-CTT-3",
    "tiRNA-Val-CAC-1-M3": "Val-CAC-1-M3",
    "tRF-1:28-Val-CAC-2": "tRF-Val-CAC-2",
    "tRF-55:76-Arg-ACG-1-M2": "tRF-Arg-ACG-1-M2",
    "tRF-1:24-Ser-CGA-1-M3": "tRF-Ser-CGA-1-M3",
}
XIONG = ["tiRNA-Val-CAC-2", "tiRNA-Lys-CTT-1-M2", "tiRNA-Lys-TTT-3-M2",
         "tiRNA-Lys-CTT-3", "tiRNA-Val-CAC-1-M3"]
LU = ["tRF-1:28-Val-CAC-2", "tRF-55:76-Arg-ACG-1-M2", "tRF-1:24-Ser-CGA-1-M3"]
ORDER = XIONG + LU


def load():
    t = pd.read_csv(os.path.join(TAB, "TableS19_marker_benchmark.csv"))
    rec = {}
    for m in ORDER:
        sub = t[t.reported_marker == m]
        d = {}
        for _, r in sub.iterrows():
            w = str(r["workflow"]).lower()
            key = "clu" if "clust" in w and "non" not in w else "non"
            d[key] = dict(l2fc=float(r["log2FoldChange"]), p=float(r["pvalue"]),
                          fdr=float(r["FDR_BH"]))
        rec[m] = d
    return rec


rec = load()


def place(vals, gap):
    """push label positions apart, keeping order"""
    order = np.argsort(vals)
    out = np.array(vals, dtype=float)
    for a, b in zip(order[:-1], order[1:]):
        if out[b] - out[a] < gap:
            out[b] = out[a] + gap
    return out


def fig4c():
    fig = plt.figure(figsize=(7.7, 7.1))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.22, 1.0], hspace=0.52,
                          left=0.135, right=0.995, top=0.955, bottom=0.135)

    # ---------------------------------------------------------------- C
    ax = fig.add_subplot(gs[0])
    left_vals = [rec[m]["non"]["l2fc"] for m in ORDER]
    right_vals = [rec[m]["clu"]["l2fc"] for m in ORDER]
    left_lab = place(left_vals, 0.78)
    right_lab = place(right_vals, 0.78)

    for m, y0, y1, ly0, ly1 in zip(ORDER, left_vals, right_vals, left_lab, right_lab):
        is_lu = m in LU
        col = L_NON if is_lu else X_NON
        ax.plot([0, 1], [y0, y1], color=col, lw=1.15, zorder=2, alpha=0.9)
        ax.scatter([0, 1], [y0, y1], s=26, color=col, zorder=3,
                   edgecolors="white", linewidths=0.5)
        for xpt, ypt, ylab, ha in ((0, y0, ly0, "right"), (1, y1, ly1, "left")):
            ax.annotate(SHORT[m], xy=(xpt, ypt), xytext=(-0.06 if ha == "right" else 1.06, ylab),
                        textcoords="data", ha=ha, va="center", fontsize=7.0,
                        color=col,
                        arrowprops=dict(arrowstyle="-", color=col, lw=0.5, alpha=0.55,
                                        shrinkA=0, shrinkB=3)
                        if abs(ylab - ypt) > 0.02 else None)

    ax.axhline(0, color="#999999", lw=0.8, ls=(0, (4, 3)), zorder=1)
    lo = min(min(left_lab), min(right_lab)) - 1.2
    hi = max(max(left_lab), max(right_lab)) + 1.2
    ax.set_xlim(-0.72, 1.72)
    ax.set_ylim(lo, hi)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["non-clustering\n(fragment level)", "clustering\n(representative level)"],
                       fontsize=8.4)
    ax.set_ylabel("log$_2$ fold change\n(disease / control)", fontsize=8.8)
    ax.tick_params(axis="y", labelsize=8.0, length=2.8)
    ax.tick_params(axis="x", length=0, pad=6)
    # "n =" would read as the number of samples here; the parenthetical
    # counts markers, so say so explicitly.
    ax.legend(handles=[Patch(color=X_NON, label="Xiong 2024 (5 markers)"),
                       Patch(color=L_NON, label="Lu 2021 (3 markers)")],
              fontsize=7.8, frameon=False, loc="lower left",
              bbox_to_anchor=(0.005, -0.015), handlelength=1.3, handletextpad=0.5)
    ax_top = ax

    # ---------------------------------------------------------------- D
    ax = fig.add_subplot(gs[1])
    n = len(ORDER)
    w = 0.185
    series = (("non", "p", 0), ("non", "fdr", 1), ("clu", "p", 2), ("clu", "fdr", 3))
    for i, m in enumerate(ORDER):
        for key, stat, j in series:
            v = -np.log10(max(rec[m][key][stat], 1e-300))
            xp = i + (j - 1.5) * w
            col = D_CLU if key == "clu" else D_NON
            if stat == "fdr":
                ax.bar(xp, v, w * 0.94, facecolor="white", edgecolor=col,
                       hatch="////", lw=0.7, zorder=4)
            else:
                ax.bar(xp, v, w * 0.94, color=col, zorder=3)

    ax.axhline(-np.log10(0.05), color="#666666", lw=0.8, ls=":")
    ax.text(-0.95, -np.log10(0.05) + 0.12, "$P$ = 0.05", fontsize=7.2, color="#555555",
            bbox=dict(fc="white", ec="none", pad=0.5))
    # dataset separator
    ax.axvline(len(XIONG) - 0.5, color="#bbbbbb", lw=0.8, ls=(0, (2, 2)), zorder=1)
    ax.set_xticks(range(n))
    ax.set_xticklabels([SHORT[m] for m in ORDER], fontsize=6.7, rotation=32,
                       ha="right", rotation_mode="anchor")
    ax.set_ylabel("$-\\log_{10}$ $P$   /   $-\\log_{10}$ BH", fontsize=8.6)
    ax.set_ylim(0, 7.2)
    ax.set_xlim(-1.05, n - 0.35)
    ax.tick_params(axis="x", length=0, pad=3)
    ax.tick_params(axis="y", labelsize=8.0, length=2.8)

    leg = [Patch(facecolor=D_NON, label="non-clustering, $P$"),
           Patch(facecolor="white", edgecolor=D_NON, hatch="////",
                 label="non-clustering, BH"),
           Patch(facecolor=D_CLU, label="clustering, $P$"),
           Patch(facecolor="white", edgecolor=D_CLU, hatch="////",
                 label="clustering, BH")]
    ax.legend(handles=leg, fontsize=7.6, frameon=False, ncol=4,
              loc="upper left", bbox_to_anchor=(-0.028, 1.185),
              handlelength=1.2, columnspacing=1.0, handletextpad=0.4,
              borderaxespad=0.0)
    # the dataset is no longer encoded by colour in panel D
    ax.text(2.0, 6.85, "Xiong 2024 (5 markers)", ha="center", va="top",
            fontsize=7.2, color="#555555")
    ax.text(6.0, 6.85, "Lu 2021 (3 markers)", ha="center", va="top",
            fontsize=7.2, color="#555555")
    ax_bot = ax

    os.makedirs(OUT, exist_ok=True)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for a_, tag in ((ax_top, "C"), (ax_bot, "D")):
        bb = a_.get_tightbbox(renderer).transformed(fig.dpi_scale_trans.inverted())
        bb = bb.expanded(1.03, 1.05)
        for ext in ("pdf", "png"):
            fig.savefig(os.path.join(OUT, f"Fig4_panel{tag}.{ext}"),
                        dpi=600, bbox_inches=bb, pad_inches=0.03)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"Fig4C_marker_reanalysis_v4.{ext}"), dpi=600)
    print("wrote Fig4 C/D v4 (both datasets) + separate panels")


if __name__ == "__main__":
    fig4c()
