#!/usr/bin/env python3
"""
阈值敏感性（R2-3）：长度比阈值对聚类拓扑与差异表达的影响。

数字来自服务器上的重跑（同一计数矩阵、同一 DE 流程）：
  0.90 与 0.94 拓扑完全相同 —— 阈值在 0.94 以下没有约束力
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "figures")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "axes.labelsize": 7.6,
    "xtick.labelsize": 6.8, "ytick.labelsize": 7, "pdf.fonttype": 42,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
C_BLUE = "#5b7c99"
C_GREY = "#c9ced3"
C_ORANGE = "#c08552"

TH = ["0.90", "0.94", "0.95", "0.98"]
# Re-run on the main-workflow count matrix (original row order), 2026-11-13.
# Baseline is 16,931 clusters, consistent with the main workflow / Table S5.
CLUSTERS = [16931, 16931, 18945, 36575]
SINGLETON_PCT = [59.6, 59.6, 63.3, 88.3]
DE_CALLS = [239, 239, 240, 406]
SHARED_WITH_094 = [239, 239, 186, 111]


def main():
    with open(os.path.join(BASE, "r1_diag", "threshold_summary.json"), "w") as fh:
        json.dump({"threshold": TH, "clusters": CLUSTERS,
                   "singleton_pct": SINGLETON_PCT, "de_calls": DE_CALLS,
                   "shared_with_0.94": SHARED_WITH_094}, fh, indent=2)

    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.5),
                             gridspec_kw={"width_ratios": [1.0, 1.0], "wspace": 0.32})
    x = np.arange(len(TH))

    ax = axes[0]
    ax.bar(x, CLUSTERS, color=C_GREY, width=0.6, zorder=3)
    for i, (c, s) in enumerate(zip(CLUSTERS, SINGLETON_PCT)):
        ax.text(i, c + 700, f"{c:,}\n{s:.0f}% single", ha="center", fontsize=5.8)
    ax.set_xticks(x)
    ax.set_xticklabels(TH)
    ax.set_ylim(0, 44000)
    ax.set_xlabel("Minimum length ratio")
    ax.set_ylabel("Clusters")
    ax.set_title("A  Cluster topology", fontsize=7.2, loc="left", pad=4)

    ax = axes[1]
    unique = [d - s for d, s in zip(DE_CALLS, SHARED_WITH_094)]
    ax.bar(x, SHARED_WITH_094, color=C_BLUE, width=0.6, label="also called at 0.94", zorder=3)
    ax.bar(x, unique, bottom=SHARED_WITH_094, color=C_ORANGE, width=0.6,
           label="only at this threshold", zorder=3)
    for i, d in enumerate(DE_CALLS):
        ax.text(i, d + 10, f"{d}", ha="center", fontsize=6.2)
    ax.set_xticks(x)
    ax.set_xticklabels(TH)
    ax.set_xlabel("Minimum length ratio")
    ax.set_ylabel("Clusters with P < 0.05")
    ax.set_ylim(0, 470)
    ax.legend(frameon=False, fontsize=5.8, loc="upper left")
    ax.set_title("B  Differential expression", fontsize=7.2, loc="left", pad=4)

    for ax in axes:
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(axis="y", color="#e8eaec", lw=0.5)
        ax.set_axisbelow(True)

    fig.savefig(os.path.join(OUT, "FigS15_threshold_sensitivity.pdf"))
    fig.savefig(os.path.join(OUT, "_build_png", "_preview_FigS15.png"), dpi=300)
    # Deliverable names used by Figure S8 / the response letter (350 dpi keeps
    # the same pixel dimensions as the previously submitted figure).
    fig.savefig(os.path.join(OUT, "FigS15_min_length_ratio.pdf"))
    fig.savefig(os.path.join(OUT, "_build_png", "FigS15_min_length_ratio.png"), dpi=350)
    plt.close(fig)
    print("saved")


if __name__ == "__main__":
    main()
