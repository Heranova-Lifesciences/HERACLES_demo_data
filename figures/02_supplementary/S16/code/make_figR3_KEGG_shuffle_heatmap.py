# -*- coding: utf-8 -*-
"""Major 2.6: KEGG pathway sets of blacklist-filtered shuffled controls.

Same layout as the reference letter figure:
  Panel A  21 x 21 heatmap of pairwise Jaccard between the real (observed) run
           and 20 dinucleotide-shuffled replicates; real row/column boxed red.
  Panel B  Number of shuffles (of 20) recovering each real KEGG term.

The 63 Enrichr runs in server_results/results/shuffle_rep242/enrichr/ were
produced from blacklist-filtered target-gene lists (CDS_black_list.txt removed
before Enrichr), so no new KEGG run is needed.  Because that filtering leaves
~94-118 genes, no term reaches BH-adjusted p < 0.05 (minimum q = 0.178), so the
reference figure's BH criterion yields zero terms.  Two usable criteria are
offered:
  --metric top20  : top 20 terms per run, fixed set size (default, closest to
                    the reference figure)
  --metric nom    : nominal p < 0.05 terms (matches Tables S28-S30)

Usage: python make_figR3_KEGG_shuffle_heatmap.py [--metric top20|nom] [--energy 29]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
E = os.path.join(ROOT, "server_results", "results", "shuffle_rep242", "enrichr")
FIG = os.path.join(ROOT, "figures")
KS = [f"k{k:03d}" for k in range(1, 21)]
TAGS = ["obs"] + KS
BLUE, GREY, RED = "#31688e", "#9a9a9a", "#b2182b"

plt.rcParams.update({
    "font.size": 8.0, "font.family": "DejaVu Sans",
    "axes.linewidth": 0.7, "pdf.fonttype": 42,
})


def terms(tag, e, metric):
    d = pd.read_csv(os.path.join(E, f"{tag}_e{e}.csv"))
    if metric == "nom":
        return set(d.loc[d["P-value"] < 0.05, "Term"])
    return set(d.nsmallest(20, "P-value")["Term"])


def jac(a, b):
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b) if (a | b) else 1.0


def main(metric, energy):
    S = {t: terms(t, energy, metric) for t in TAGS}
    M = np.array([[jac(S[a], S[b]) for b in TAGS] for a in TAGS], dtype=float)
    cnt = {x: sum(1 for t in KS if x in S[t]) for x in S["obs"]}
    # Only terms recovered in at least one shuffle are displayed (0/20 terms are
    # ranking-cutoff artefacts or singleton-gene terms; see analysis notes).
    ranked = sorted((kv for kv in cnt.items() if kv[1] > 0),
                    key=lambda kv: (kv[1], kv[0]))

    fig = plt.figure(figsize=(10.2, 4.4))
    axA = fig.add_axes([0.075, 0.16, 0.33, 0.72])
    axC = fig.add_axes([0.415, 0.16, 0.012, 0.72])
    axB = fig.add_axes([0.615, 0.16, 0.36, 0.72])

    # ---- A: heatmap ----
    Mp = M.copy()
    np.fill_diagonal(Mp, np.nan)
    cmap = plt.get_cmap("Blues").copy()
    cmap.set_bad("white")
    lo = 0.0 if metric == "nom" else 0.10
    hi = 0.60 if metric == "nom" else 0.75
    im = axA.imshow(Mp, cmap=cmap, vmin=lo, vmax=hi, aspect="auto", interpolation="nearest")
    axA.add_patch(Rectangle((-0.5, -0.5), 21, 1, fill=False, edgecolor=RED, lw=1.6, zorder=5))
    axA.add_patch(Rectangle((-0.5, -0.5), 1, 21, fill=False, edgecolor=RED, lw=1.6, zorder=5))
    ticks = [0, 4, 8, 12, 16, 20]
    labs = ["real"] + [f"s{i}" for i in ticks[1:]]
    axA.set_xticks(ticks); axA.set_xticklabels(labs, fontsize=7.4)
    axA.set_yticks(ticks); axA.set_yticklabels(labs, fontsize=7.4)
    axA.tick_params(length=0)
    for sp in axA.spines.values():
        sp.set_linewidth(0.6)
    axA.set_title("A  pathway sets: real vs 20 shuffles", loc="left", fontsize=9.2, pad=8)
    cb = fig.colorbar(im, cax=axC)
    cb.set_label("Jaccard", fontsize=7.8)
    cb.ax.tick_params(labelsize=6.8)
    cb.outline.set_linewidth(0.5)

    # ---- B: recovery of real terms ----
    y = np.arange(len(ranked))
    vals = np.array([v for _, v in ranked], dtype=float)
    cols = [BLUE if v >= 15 else GREY if v >= 10 else RED for v in vals]
    axB.barh(y, vals, color=cols, height=0.72, zorder=3)
    axB.set_yticks(y)
    axB.set_yticklabels([n if len(n) <= 36 else n[:35] + "\u2026" for n, _ in ranked], fontsize=6.2)
    axB.tick_params(axis="y", length=0, pad=3)
    axB.set_xlim(0, 22)
    axB.set_xticks([0, 5, 10, 15, 20])
    axB.set_xlabel("number of shuffles recovering the term (of 20)", fontsize=8.0)
    axB.set_title("B  real KEGG terms recur after shuffling", loc="left", fontsize=9.2, pad=8)
    for yy, v in zip(y, vals):
        axB.text(v + 0.45, yy, f"{int(v)}/20", va="center", fontsize=6.4, color="#333333")
    axB.grid(axis="x", color="#eceef0", lw=0.6)
    axB.set_axisbelow(True)
    for sp in ("top", "right"):
        axB.spines[sp].set_visible(False)

    # No footer note: the figure carries no library/version/criterion annotation.
    stem = f"FigR3_KEGG_shuffle_heatmap_{metric}"
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{stem}.{ext}"), dpi=300)
    plt.close(fig)
    n_zero = sum(1 for v in cnt.values() if v == 0)
    print(f"wrote {stem}.pdf/.png | shown terms={len(ranked)} (0/20 hidden={n_zero}) "
          f"| 20/20={int((vals==20).sum())} >=15={int((vals>=15).sum())}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--metric", choices=["top20", "nom"], default="top20")
    ap.add_argument("--energy", type=int, default=29)
    a = ap.parse_args()
    main(a.metric, a.energy)
