# -*- coding: utf-8 -*-
"""R3 Major 3.3 letter figure, rebuilt on the ablation_v2 rerun.

Same layout as make_figR3_3p3_mapping.py (panel A stacked bars of shared /
unique fragments, panel B Venn of the two nominally significant sets), but all
numbers come from output/ablation_v2/report.json, where the -v 0 arm reproduces
the main-workflow baseline of 43,554 fragments (the old arm used max-len 45 and
gave 37,669).

Writes figures/FigR3_mapping_ablation.{png,pdf}.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(BASE, "output", "ablation_v2", "report.json")
OUT = os.path.join(BASE, "figures")

C_SHARED = "#b9c6d4"
C_V0 = "#5b7c99"
C_V1 = "#c0743c"
DARK = "#333333"
GRID = "#e6e6e6"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.6,
    "axes.labelsize": 9.2,
    "xtick.labelsize": 8.4,
    "ytick.labelsize": 8.4,
    "pdf.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})


def main():
    r = json.load(open(SRC, encoding="utf8"))
    v0, v1 = r["v0"], r["v1"]
    n0, n1 = v0["n_features"], v1["n_features"]
    shared, lost, gained = r["common"], r["only_v0"], r["only_v1"]
    reads_new = 100 * r["reads_share_only_v1"]
    a, b, both = v0["n_sig_p"], v1["n_sig_p"], r["sig_overlap"]

    fig = plt.figure(figsize=(7.4, 3.25))
    axA = fig.add_axes([0.085, 0.155, 0.40, 0.70])
    axB = fig.add_axes([0.545, 0.11, 0.44, 0.76])

    # ---- A: fragment turnover ----
    xs = [0, 1]
    axA.bar(xs[0], shared, width=0.52, color=C_SHARED, edgecolor="white", lw=0.6, zorder=3)
    axA.bar(xs[0], lost, bottom=shared, width=0.52, color=C_V0, edgecolor="white", lw=0.6, zorder=3)
    axA.bar(xs[1], shared, width=0.52, color=C_SHARED, edgecolor="white", lw=0.6, zorder=3)
    axA.bar(xs[1], gained, bottom=shared, width=0.52, color=C_V1, edgecolor="white", lw=0.6, zorder=3)

    axA.set_xticks(xs)
    axA.set_xticklabels(["\u2212v 0\n(default)", "\u2212v 1\n(one mismatch)"])
    axA.set_ylabel("Fragments")
    axA.set_title("A  What the mapping flag changes", loc="left", fontsize=9.6, pad=6)
    axA.set_xlim(-0.55, 1.55)
    axA.set_ylim(0, 64000)
    axA.set_yticks([0, 20000, 40000, 60000])
    axA.set_yticklabels(["0", "20,000", "40,000", "60,000"])
    axA.grid(axis="y", color=GRID, zorder=0)
    axA.set_axisbelow(True)
    for s in ("top", "right"):
        axA.spines[s].set_visible(False)

    axA.text(0, n0 + 1600, f"{n0:,}", ha="center", fontsize=8.6, color=DARK)
    axA.text(1, n1 + 3400, f"{n1:,}", ha="center", fontsize=8.6, color=DARK, fontweight="bold")
    axA.text(1, n1 + 900, f"(+{100 * (n1 / n0 - 1):.0f}%)", ha="center", fontsize=7.6, color="#8a4a20")
    axA.text(0.5, shared * 0.46, f"{shared:,} shared", ha="center", fontsize=7.8, color="#4a5b6b")

    axA.annotate(f"{gained:,} new at \u2212v 1,\ncarrying {reads_new:.1f}% of reads",
                 xy=(1.27, shared + gained * 0.55), xytext=(1.32, 50500),
                 fontsize=7.8, color="#8a4a20", ha="left",
                 arrowprops=dict(arrowstyle="-", color="#c0743c", lw=0.8, shrinkA=0, shrinkB=3))
    axA.annotate(f"{lost:,} lost", xy=(-0.27, shared + lost * 0.5), xytext=(-0.53, 47200),
                 fontsize=7.8, color="#33506b", ha="left",
                 arrowprops=dict(arrowstyle="-", color=C_V0, lw=0.8, shrinkA=0, shrinkB=3))

    # ---- B: overlap of the two call sets ----
    venn2(subsets=(a - both, b - both, both), set_labels=None, ax=axB,
          set_colors=(C_V0, C_V1), alpha=0.35)
    venn2_circles(subsets=(a - both, b - both, both), ax=axB, color="#888888", lw=0.7)
    for t in axB.texts:
        t.set_fontsize(8.0)
        t.set_color(DARK)
        if t.get_text() == str(both):
            t.set_fontweight("bold")
            t.set_fontsize(9.2)

    axB.set_title("B  Overlap of the two significant sets", loc="left", fontsize=9.6, pad=6)
    axB.set_axis_off()
    axB.text(0.02, 0.97, f"\u2212v 0  n = {a}", transform=axB.transAxes, fontsize=7.8,
             color=C_V0, va="top")
    axB.text(0.98, 0.97, f"\u2212v 1  n = {b}", transform=axB.transAxes, fontsize=7.8,
             color=C_V1, va="top", ha="right")
    axB.text(0.5, -0.06, f"Jaccard {r['jaccard_all']:.3f}", transform=axB.transAxes,
             ha="center", fontsize=8.2, color="#333333")

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"FigR3_mapping_ablation.{ext}"),
                    dpi=300 if ext == "png" else None)
    print(f"v0 {n0:,} v1 {n1:,} | shared {shared:,} new {gained:,} lost {lost:,} "
          f"| sig {a}/{b} both {both} | Jaccard {r['jaccard_all']:.3f}")
    print("saved FigR3_mapping_ablation.png/.pdf")


if __name__ == "__main__":
    main()
