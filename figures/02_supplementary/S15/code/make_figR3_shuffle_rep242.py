# -*- coding: utf-8 -*-
"""Major 2.6 (R3) letter figure, rebuilt on the 242-cluster-representative shuffle.

Identical layout AND identical metric to the original FigR3_shuffle_thresholds:
panel A plots the Jaccard of top-5% target-gene lists (real vs shuffled, and
shuffled vs shuffled) with IQR bands and a red median-gap arrow; panel B plots
per-gene reproducibility of the observed targets across the 20 shuffles.

Only the underlying analysis changed: 242 cluster representatives instead of
444 individual fragments, blacklist genes removed, score ordering aligned with
Figure 3B.  Inputs are the 63 gene lists in
server_results/results/shuffle_rep242/genes/ (one per arm x energy), which also
feed Table S28.

Output figures/FigR3_shuffle_thresholds.pdf / .png
"""
import itertools
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G = os.path.join(ROOT, "server_results", "results", "shuffle_rep242", "genes")
FIG = os.path.join(ROOT, "figures")

KS = [f"k{k:03d}" for k in range(1, 21)]
XS = [-40, -36, -29]
BLUE, GREY, RED = "#31688e", "#9a9a9a", "#b2182b"

plt.rcParams.update({
    "font.size": 8.0, "font.family": "DejaVu Sans",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "pdf.fonttype": 42,
})


def genes(tag, e):
    with open(os.path.join(G, f"{tag}_e{abs(e)}_genes.txt"), encoding="utf8") as f:
        return set(f.read().split())


def jac(a, b):
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b) if (a | b) else 1.0


def collect():
    out = {}
    for e in XS:
        obs = genes("obs", e)
        S = {t: genes(t, e) for t in KS}
        ov = [jac(obs, S[t]) for t in KS]
        pv = [jac(S[a], S[b]) for a, b in itertools.combinations(KS, 2)]
        cnt = {g: sum(1 for t in KS if g in S[t]) for g in obs}
        out[e] = dict(
            n_obs=len(obs),
            obs_med=float(np.median(ov)), shuf_med=float(np.median(pv)),
            oi=np.percentile(ov, [25, 75]), si=np.percentile(pv, [25, 75]),
            gap=float(np.median(pv) - np.median(ov)),   # same sign convention as original
            never=sum(v == 0 for v in cnt.values()),
            once=sum(0 < v < 20 for v in cnt.values()),
            always=sum(v == 20 for v in cnt.values()),
        )
    return out


def main():
    d = collect()
    om = [d[e]["obs_med"] for e in XS]
    sm = [d[e]["shuf_med"] for e in XS]
    oi = np.array([d[e]["oi"] for e in XS])
    si = np.array([d[e]["si"] for e in XS])

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9),
                             gridspec_kw=dict(width_ratios=[1.15, 1], wspace=0.30,
                                              left=0.105, right=0.985, top=0.86, bottom=0.295))

    # ---- A: median Jaccard with IQR bands ----
    ax = axes[0]
    ax.fill_between(XS, oi[:, 0], oi[:, 1], color=BLUE, alpha=0.18, lw=0)
    ax.fill_between(XS, si[:, 0], si[:, 1], color=GREY, alpha=0.22, lw=0)
    ax.plot(XS, om, "-o", color=BLUE, ms=5, lw=1.3, label="real vs shuffled")
    ax.plot(XS, sm, "--o", color=GREY, mec="#666666", mfc="white", ms=5, lw=1.3,
            label="shuffled vs shuffled")
    for x, o, s in zip(XS, om, sm):
        ax.annotate("", xy=(x, s - 0.008), xytext=(x, o + 0.008),
                    arrowprops=dict(arrowstyle="<->", color=RED, lw=0.9))
        ax.text(x + 0.35, (o + s) / 2, "\u0394=%.3f" % (s - o), color=RED,
                fontsize=6.8, va="center")
    ax.set_xticks(XS)
    ax.set_xticklabels(["e = \u221240", "e = \u221236", "e = \u221229\n(default)"])
    ax.set_xlim(-41, -27.4)
    ax.set_ylim(0.10, 0.62)
    ax.set_ylabel("Jaccard of top-5% target lists")
    ax.legend(frameon=False, fontsize=7, loc="lower right", handlelength=1.6)
    ax.set_title("A  real lists leave the shuffle cloud at every threshold",
                 loc="left", fontsize=8.4)
    ax.grid(axis="y", color="#eceef0", lw=0.6)
    ax.set_axisbelow(True)

    # ---- B: per-gene reproducibility ----
    ax = axes[1]
    order = [-29, -36, -40]
    y = np.arange(3)[::-1]
    left = np.zeros(3)
    for key, col, lab in [("always", BLUE, "in all 20 shuffles"),
                          ("once", "#d8d8d8", "in 1\u201319"),
                          ("never", RED, "in none")]:
        vals = np.array([d[e][key] for e in order], dtype=float)
        pct = 100 * vals / np.array([d[e]["n_obs"] for e in order])
        ax.barh(y, pct, left=left, color=col, height=0.58, lw=0, label=lab)
        for yy, l, v, p in zip(y, left, vals, pct):
            if p > 6:
                ax.text(l + p / 2, yy, "%d" % int(v), ha="center", va="center",
                        fontsize=6.6, color="white" if col in (BLUE, RED) else "#333333")
        left += pct
    ax.set_yticks(y)
    ax.set_yticklabels(["e = \u221229", "e = \u221236", "e = \u221240"])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of observed target genes (n = 94 / 118 / 106)")
    ax.legend(frameon=False, fontsize=6.6, loc="upper center", bbox_to_anchor=(0.5, -0.32),
              ncol=3, handlelength=1.1, columnspacing=0.9)
    ax.set_title("B  per-gene reproducibility across shuffles", loc="left", fontsize=8.4)

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"FigR3_shuffle_thresholds.{ext}"), dpi=300)
    plt.close(fig)
    print("wrote FigR3_shuffle_thresholds.pdf/.png (242-rep, gene-list Jaccard)")
    for e in XS:
        k = d[e]
        print(f"  e={e}: n_obs={k['n_obs']} obs_med={k['obs_med']:.3f} "
              f"shuf_med={k['shuf_med']:.3f} gap={k['gap']:.3f} "
              f"recovery never/once/always={k['never']}/{k['once']}/{k['always']}")


if __name__ == "__main__":
    main()
