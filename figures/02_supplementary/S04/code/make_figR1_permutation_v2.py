# -*- coding: utf-8 -*-
"""R1 Comment 6 letter figure, version 2.

Only the two workflows the reviewer's comment is about -- graph clustering and no
clustering (fragment level) -- in three datasets: this cohort, Paper 1 (Xiong) and
Paper 2 (Lu).  The parent-tRNA / +5' anchor / representative-only sweeps remain in
Table S15 and are no longer drawn.

Bars are the observed number of nominal calls (P < 0.05); the permuted-label null is
drawn as median + IQR box + full-range whisker.  The split that reproduces the true
assignment (and, in the enumerated six-sample design, its complement, which is the
same test with the labels swapped) is not a null draw and is dropped.

Outputs
    figures/FigR1_permutation.pdf/.png
    output/R1_回应包/03_tables/TableS25b_permutation_two_workflows.csv
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "output", "补充材料_提交版")
FIG = os.path.join(ROOT, "figures")
OUT = os.path.join(ROOT, "output", "R1_回应包", "03_tables")


def null_of(d, workflow):
    d = d[d["workflow"] == workflow]
    true = d[d["perm"] == -1].iloc[0]
    null = d[d["perm"] != -1].copy()
    # 1) the identity partition and its complement reproduce the observed test
    null = null[~np.isclose(null["median_p"], true["median_p"], rtol=0, atol=1e-4)]
    # 2) the six-sample design enumerates each split with its complement, which gives
    #    the same p-values; keep one copy of each distinct assignment
    key = (null["n_tested"].astype(int).astype(str) + "|"
           + null["n_sig_p"].astype(int).astype(str) + "|"
           + null["n_sig_fdr"].astype(int).astype(str) + "|"
           + null["median_p"].round(3).astype(str))
    null = null.loc[~key.duplicated()]
    return true, null


def main():
    s25 = pd.read_csv(os.path.join(TAB, "TableS16_permuted_label_controls.csv"))
    plan = [("cohort (n=46, 26v20)", "A", "In-house cohort (n = 46)", "clust", "frag"),
            ("paper1 (n=6, 3v3)", "B", "Xiong (n = 6)", "clust", "nonclust"),
            ("paper2 (n=8, 4v4)", "C", "Lu (n = 8)", "clust", "nonclust")]
    names = {"clust": "graph\nclustering", "frag": "no clustering\n(fragments)",
             "nonclust": "no clustering\n(fragments)"}

    rows = []
    for ds, tag, title, wc, wn in plan:
        d = s25[s25["dataset"] == ds]
        for wf in (wc, wn):
            true, null = null_of(d, wf)
            v = null["n_sig_p"].to_numpy(float)
            rows.append(dict(panel=tag, title=title, workflow=wf, tests=int(true["n_tested"]),
                             obs=int(true["n_sig_p"]), fdr=int(true["n_sig_fdr"]),
                             n_null=len(v), med=float(np.median(v)),
                             lo=float(np.percentile(v, 25)), hi=float(np.percentile(v, 75)),
                             mn=float(v.min()), mx=float(v.max()),
                             n_at_or_above=int((v >= true["n_sig_p"]).sum())))
    res = pd.DataFrame(rows)
    print(res.to_string(index=False))
    os.makedirs(OUT, exist_ok=True)
    res.to_csv(os.path.join(OUT, "TableS25b_permutation_two_workflows.csv"), index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "font.family": "DejaVu Sans",
                         "axes.linewidth": 0.8})

    blue, grey = "#31688e", "#c2c2c2"
    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.72))
    fig.subplots_adjust(left=0.080, right=0.988, top=0.875, bottom=0.285, wspace=0.30)

    for ax, group in zip(axes, [res[res["panel"] == p] for p in "ABC"]):
        ymax = max(group["obs"].max(), group["mx"].max()) * 1.30
        for xi, (_, r) in enumerate(group.iterrows()):
            ax.bar(xi - 0.17, r["obs"], 0.32, color=blue, zorder=3)
            ax.bar(xi + 0.17, r["med"], 0.32, facecolor=grey, edgecolor="#8a8a8a",
                   lw=0.6, zorder=3)
            ax.errorbar(xi + 0.17, r["med"], yerr=[[r["med"] - r["mn"]], [r["mx"] - r["med"]]],
                        fmt="none", ecolor="#5a5a5a", elinewidth=0.8, capsize=2.2, zorder=4)
            ax.text(xi - 0.17, r["obs"] + ymax * 0.018, f"{r['obs']}", ha="center",
                    va="bottom", fontsize=7.9, color=blue)
            ax.text(xi + 0.17, r["mx"] + ymax * 0.018, f"{r['med']:.0f}", ha="center",
                    va="bottom", fontsize=7.9, color="#6d6d6d")
        ax.set_xticks([0, 1])
        ax.set_xticklabels([names[w] for w in group["workflow"]], fontsize=7.6)
        ax.set_xlim(-0.55, 1.55)
        ax.set_ylim(0, ymax)
        ax.set_title(f"{group['panel'].iloc[0]}   {group['title'].iloc[0]}",
                     loc="left", fontsize=8.3, pad=3)
        ax.tick_params(axis="both", labelsize=7.5, length=2.8)

    axes[0].set_ylabel("nominal calls (P < 0.05)", fontsize=8.0, labelpad=2)
    handles = [plt.Rectangle((0, 0), 1, 1, color=blue),
               plt.Rectangle((0, 0), 1, 1, facecolor=grey, edgecolor="#8a8a8a", lw=0.6)]
    fig.legend(handles, ["observed calls",
                         "permuted labels: median, error bar = full range"],
               fontsize=7.0, frameon=False, loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, -0.015), handlelength=1.1, handleheight=0.75,
               columnspacing=1.6, borderpad=0.1)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"FigR1_permutation.{ext}"), dpi=300)
    print("saved", os.path.join(FIG, "FigR1_permutation.png"))


if __name__ == "__main__":
    main()
