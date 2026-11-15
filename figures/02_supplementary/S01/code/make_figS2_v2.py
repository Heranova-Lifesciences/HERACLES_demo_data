#!/usr/bin/env python3
"""
成员相对簇代表的差异构成、偏移量分布与丰度比（R1-2 新口径）。
输出 figures/FigS2_cluster_relationships.pdf
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
    "font.family": "DejaVu Sans",
    "font.size": 7,
    "axes.labelsize": 7.6,
    "xtick.labelsize": 6.6,
    "ytick.labelsize": 7,
    "pdf.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.08,
})

ORDER = ["5p_trimmed", "3p_trimmed", "5p_extended", "3p_extended",
         "both_trimmed", "both_extended", "shift_5p_trim_3p_ext", "shift_5p_ext_3p_trim"]
LAB = {"5p_trimmed": "5\u2032 trimmed", "3p_trimmed": "3\u2032 trimmed",
       "5p_extended": "5\u2032 extended", "3p_extended": "3\u2032 extended",
       "both_trimmed": "both trimmed", "both_extended": "both extended",
       "shift_5p_trim_3p_ext": "3\u2032 extended, 5\u2032 trimmed",
       "shift_5p_ext_3p_trim": "5\u2032 extended, 3\u2032 trimmed"}
C_ONE = "#5b7c99"
C_TWO = "#9db4c8"
C_MID = "#d9c3a5"
C_FAR = "#c08552"


def main():
    d = json.load(open(os.path.join(BASE, "r1_diag", "member_vs_rep.json")))
    n = {k: len(d[k]) for k in ORDER}
    tot = sum(n.values())
    y = np.arange(len(ORDER))[::-1]

    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.9),
                             gridspec_kw={"width_ratios": [1.15, 1.35, 1.0], "wspace": 0.14})

    # A 构成
    ax = axes[0]
    ax.barh(y, [100 * n[k] / tot for k in ORDER], color=C_ONE, height=0.62, zorder=3)
    for yy, k in zip(y, ORDER):
        ax.text(100 * n[k] / tot + 0.4, yy, f"{100*n[k]/tot:.1f}%", va="center", fontsize=6.0)
    ax.set_yticks(y)
    ax.set_yticklabels([LAB[k] for k in ORDER])
    ax.set_xlim(0, 22)
    ax.set_xlabel("Members (%)")
    ax.set_title("A  Relation to representative", fontsize=7.2, loc="left", pad=4)

    # B 偏移量分布
    ax = axes[1]
    left = np.zeros(len(ORDER))
    cats = [("offset 1", lambda o: o == 1, C_ONE),
            ("offset 2", lambda o: o == 2, C_TWO),
            ("offset 3", lambda o: o == 3, C_MID),
            ("offset \u2265 4", lambda o: o >= 4, C_FAR)]
    for name, fn, col in cats:
        w = [100 * float(np.mean(fn(np.array([p[0] for p in d[k]])))) for k in ORDER]
        ax.barh(y, w, left=left, color=col, height=0.62, zorder=3, label=name)
        left += np.array(w)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlim(0, 100)
    ax.set_xlabel("Members (%)")
    ax.legend(frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.55, -0.34),
              fontsize=5.8, handlelength=0.9, columnspacing=1.0)
    ax.set_title("B  Terminal offset", fontsize=7.2, loc="left", pad=4)

    # C 丰度比
    ax = axes[2]
    med = [float(np.median([p[1] for p in d[k]])) for k in ORDER]
    q1 = [float(np.percentile([p[1] for p in d[k]], 25)) for k in ORDER]
    q3 = [float(np.percentile([p[1] for p in d[k]], 75)) for k in ORDER]
    ax.hlines(y, q1, q3, color="#c7ccd1", lw=2.4, zorder=2)
    ax.scatter(med, y, s=16, color="#2c3e50", zorder=3)
    ax.axvline(0, lw=0.6, color="#888888", ls=(0, (2, 2)))
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlim(-12, 1)
    ax.set_xlabel("log$_2$(member / rep.)")
    ax.set_title("C  Abundance ratio", fontsize=7.2, loc="left", pad=4)

    for ax in axes:
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(axis="x", color="#e8eaec", lw=0.5)
        ax.set_axisbelow(True)

    fig.savefig(os.path.join(OUT, "FigS2_cluster_relationships.pdf"))
    fig.savefig(os.path.join(OUT, "_build_png", "_preview_FigS2.png"), dpi=300)
    plt.close(fig)
    print(f"saved; n={tot:,}")


if __name__ == "__main__":
    main()
