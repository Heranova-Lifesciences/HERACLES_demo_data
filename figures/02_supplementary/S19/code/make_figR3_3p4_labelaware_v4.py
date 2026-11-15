# -*- coding: utf-8 -*-
"""Letter Figure R3-10 (final): label-aware re-clustering, conclusion-first.

No legend, no overlapping annotations. Titles state the conclusions; numbers
sit inside the bars; the single quantitative contrast is called out in place.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = r"C:\Users\huyao\Desktop\HERACLES修改"
FIG = os.path.join(ROOT, "figures")

DSETS = ["Xiong 2024", "Lu 2021", "In-house 46"]
SUB = ["3 vs 3 libraries", "4 vs 4 libraries", "20 vs 26 libraries"]
spec_frag = [59.15, 63.25, 49.48]
spec_read = [0.26, 0.28, 3.17]
fold = [228, 226, 16]
true_w = [95.66, 91.51, 77.35]
perm_w = [[96.64, 97.94, 97.94], [89.92, 92.89, 84.10], [77.20, 76.39, 77.55]]
perm_mean = [float(np.mean(p)) for p in perm_w]
delta = [t - m for t, m in zip(true_w, perm_mean)]

BLUE = "#5b87b5"
RED = "#b03a2e"
GREY = "#9a9a9a"

plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.9))
fig.subplots_adjust(left=0.09, right=0.985, top=0.78, bottom=0.25, wspace=0.32)
fig.suptitle("Splitting libraries by condition reveals no condition-specific cluster structure",
             x=0.09, ha="left", fontsize=9.6, fontweight="bold", y=0.965)

x = np.arange(3)
w = 0.36

# ---------------- (A) ----------------
ax = axes[0]
b1 = ax.bar(x - w / 2, spec_frag, w, color=BLUE)
b2 = ax.bar(x + w / 2, spec_read, w, color=RED)
for b, v in zip(b1, spec_frag):
    ax.text(b.get_x() + b.get_width() / 2, v + 1.5, f"{v:.1f}%",
            ha="center", fontsize=8.2, fontweight="bold", color=BLUE)
for b, v, f in zip(b2, spec_read, fold):
    cx = b.get_x() + b.get_width() / 2
    tx = cx + 0.05
    ax.text(tx, 12.5, f"{v:.2f}%", ha="left", fontsize=8.2,
            fontweight="bold", color=RED)
    ax.text(tx, 8.6, f"×{f} less", ha="left", fontsize=7.2, color="#7a3b34")
    ax.annotate("", xy=(cx, max(v, 0.8)), xytext=(cx, 6.8),
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.9))
ax.text(0.01, 0.965, "blue = share of fragments     red = share of reads",
        transform=ax.transAxes, fontsize=6.8, color="#555555", va="top")
ax.set_xlim(-0.6, 2.95)
ax.set_xticks(x)
ax.set_xticklabels([f"{d}\n{s}" for d, s in zip(DSETS, SUB)], fontsize=7)
ax.set_ylabel("Percentage of all fragments / all reads")
ax.set_ylim(0, 78)
ax.set_title("A  Many condition-specific fragments,\n    almost no reads",
             loc="left", fontsize=8.8, fontweight="bold", pad=8)

# ---------------- (B) ----------------
ax = axes[1]
bt = ax.bar(x - w / 2, true_w, w, color=RED)
bp = ax.bar(x + w / 2, perm_mean, w, color=GREY,
            yerr=[[m - min(p) for m, p in zip(perm_mean, perm_w)],
                  [max(p) - m for m, p in zip(perm_mean, perm_w)]],
            error_kw=dict(ecolor="#707070", lw=0.9, capsize=3))
for bars, vals in ((bt, true_w), (bp, perm_mean)):
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v - 9, f"{v:.1f}%",
                ha="center", fontsize=8.0, fontweight="bold", color="white")
ax.text(0.5, 0.965,
        "difference from random (left → right): " +
        " / ".join(f"{d:+.1f}" for d in delta) + " pp",
        transform=ax.transAxes, ha="center", va="top", fontsize=7.0, color="#444444")
ax.set_xticks(x)
ax.set_xticklabels([f"{d}\n{s}" for d, s in zip(DSETS, SUB)], fontsize=7)
ax.set_ylabel("Read-weighted representative agreement (%)")
ax.set_ylim(0, 110)
ax.set_title("B  True labels ≈ random splits\n    (red = true; grey = 3 random splits, ± range)",
             loc="left", fontsize=8.2, fontweight="bold", pad=8)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, f"FigR3_3p4_labelaware.{ext}"), dpi=300,
                facecolor="white")
print("saved FigR3_3p4_labelaware.png / .pdf (final)")
