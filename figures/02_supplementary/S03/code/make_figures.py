#!/usr/bin/env python3
"""
生成 Reviewer 1 回应所需的 6 张图（PDF 矢量）。

输出到 figures/：
  FigS1_read_length_boundary.pdf     R1-3  读长边界与 unsupported 3' 端
  FigS2_cluster_relationships.pdf    R1-2  簇内关系构成、末端偏移、丰度比、跨度
  FigS3_singleton_structure.pdf      R1-5  单例结构与非单例分层
  FigS4_boundary_exclusion.pdf       R1-3  剔除边界片段后的 DE 敏感性
  FigS5_threshold_sensitivity.pdf    R1-9  能量阈值敏感性
  FigS6_rep_vs_members.pdf           R1-7  代表序列与成员靶基因的重叠
"""
import collections
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "pdf.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

C_DARK = "#2c3e50"
C_STEEL = "#5b7c99"
C_TERRA = "#c08552"
C_RED = "#c0392b"
C_GRAY = "#8a9199"
C_LGRAY = "#d5dadd"


def style(ax, grid_axis="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#444444")
    ax.spines["bottom"].set_color("#444444")
    if grid_axis:
        ax.grid(axis=grid_axis, color="#e8eaec", lw=0.5, zorder=0)
        ax.set_axisbelow(True)


def panel(ax, letter, dx=-0.13, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="left", color="black")


# ---------------------------------------------------------------- 数据载入
def load_clusters():
    path = os.path.join(BASE, "cluster_summary_detailed.csv")
    cl = collections.defaultdict(list)
    with open(path) as f:
        for row in csv.DictReader(f):
            if row["is_cluster_total"] == "True":
                continue
            cl[row["cluster_id"]].append(row)
    return cl


def parse_id(tid):
    a = tid.split(":")
    s, e = map(int, a[1].split("-"))
    return a[0], s, e, a[2]


def load_trna_lengths():
    path = os.path.join(BASE, "HERACLES", "tRNA_index", "hg38-tRNA.fa")
    ref, cur, n = {}, None, 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur:
                    ref[cur] = n
                cur = line[1:].split()[0]
                n = 0
            elif cur:
                n += len(line)
        if cur:
            ref[cur] = n
    return ref


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------- Fig S1
def fig_read_length(cl, rlen):
    by_len = collections.Counter()
    at_end = collections.Counter()
    total = collections.Counter()
    dist_l50 = collections.Counter()
    dist_ref = collections.Counter()

    for ms in cl.values():
        for m in ms:
            p, s, e, _ = parse_id(m["tsRNA_id"])
            L = e - s + 1
            by_len[L] += 1
            if rlen.get(p) == e:
                at_end[L] += 1
            d = rlen.get(p, 0) - e
            if d < 0:
                continue
            total[L] += 1
            if L == 50:
                dist_l50[min(d, 7)] += 1
            elif 43 <= L <= 49:
                dist_ref[min(d, 7)] += 1

    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.05),
                             gridspec_kw={"width_ratios": [1.15, 1.15, 1.0],
                                          "wspace": 0.46})

    # A: length distribution
    ax = axes[0]
    xs = sorted(by_len)
    ax.bar(xs, [by_len[x] for x in xs], color=C_GRAY, width=0.85, zorder=3)
    ymax = max(by_len.values())
    for v in (18, 50):
        ax.axvline(v, ls="--", lw=0.8, color=C_RED, zorder=4)
    ax.set_xlabel("Fragment length (nt)")
    ax.set_ylabel("Unique fragments")
    ax.set_ylim(0, ymax * 1.26)
    ax.set_xlim(15, 53)
    ax.text(19.2, ymax * 1.19, "18 nt filter", fontsize=5.4, color=C_RED, ha="left")
    ax.text(49, ymax * 1.09, "50 nt read length", fontsize=5.4, color=C_RED, ha="right")
    style(ax)
    panel(ax, "A")

    # B: fraction ending at tRNA 3' terminus
    ax = axes[1]
    lens = [L for L in range(20, 51) if by_len[L] >= 30]
    frac = [100.0 * at_end[L] / by_len[L] for L in lens]
    ax.plot(lens, frac, "o-", ms=2.6, lw=0.9, color=C_DARK, zorder=3)
    ax.axvline(50, ls="--", lw=0.8, color=C_RED, zorder=4)
    ax.scatter([50], [100.0 * at_end[50] / by_len[50]], s=26, facecolor="white",
               edgecolor=C_RED, lw=1.0, zorder=5)
    ax.set_xlabel("Fragment length (nt)")
    ax.set_ylabel("Ending at tRNA 3' terminus (%)")
    ax.set_xlim(19, 52)
    ax.set_ylim(-0.6, max(frac) * 1.30)
    ax.annotate("0 of 1,235", (50, 0), xytext=(34.5, max(frac) * 0.16),
                fontsize=5.6, color=C_RED,
                arrowprops=dict(arrowstyle="->", color=C_RED, lw=0.7,
                                connectionstyle="arc3,rad=-0.18"))
    style(ax)
    panel(ax, "B")

    # C: distance-to-terminus profile
    ax = axes[2]
    ds = list(range(7))
    n50 = sum(dist_l50.values())
    nr = sum(dist_ref.values())
    a = [100.0 * dist_l50[d] / n50 for d in ds]
    b = [100.0 * dist_ref[d] / nr for d in ds]
    a7 = 100.0 * dist_l50[7] / n50
    b7 = 100.0 * dist_ref[7] / nr
    w = 0.4
    ax.bar([d - w / 2 for d in ds], b, width=w, color=C_STEEL,
           label="43-49 nt", zorder=3)
    ax.bar([d + w / 2 for d in ds], a, width=w, color=C_RED,
           label="50 nt", zorder=3)
    ax.set_xticks(ds)
    ax.set_xlabel(f"Distance from tRNA 3\u2032 terminus (nt)\n(d \u2265 7 omitted: {a7:.0f}% vs {b7:.0f}%)",
                  fontsize=6.2)
    ax.set_ylabel("Fragments (%)")
    ax.set_ylim(0, max(max(a), max(b)) * 1.45)
    ax.legend(frameon=False, loc="upper right", handlelength=1.0,
              borderaxespad=0.15)
    style(ax)
    panel(ax, "C")

    fig.savefig(os.path.join(OUT, "FigS1_read_length_boundary.pdf"))
    fig.savefig(os.path.join(OUT, "_preview_FigS1.png"), dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------- Fig S2
def fig_relationships():
    near = [r for r in read_csv(os.path.join(BASE, "r1_diag", "edges_nearest.csv"))
            if r["same_parent"] == "1"]
    star = read_csv(os.path.join(BASE, "r1_diag", "edges_star.csv"))
    n = len(near)

    labels = {
        "one_end_3p_shorten": "3' end shorter",
        "one_end_3p_lengthen": "3' end longer",
        "one_end_5p_lengthen": "5' end longer",
        "one_end_5p_shorten": "5' end shorter",
        "shift_5p_extend_3p_trim": "5' longer + 3' shorter",
        "shift_5p_trim_3p_extend": "5' shorter + 3' longer",
        "both_ends_in": "both ends shorter",
        "both_ends_out": "both ends longer",
    }
    order = ["one_end_5p_shorten", "one_end_5p_lengthen",
             "one_end_3p_shorten", "one_end_3p_lengthen",
             "shift_5p_extend_3p_trim", "shift_5p_trim_3p_extend",
             "both_ends_in", "both_ends_out"]
    cnt = collections.Counter(r["typ"] for r in near)

    fig = plt.figure(figsize=(7.1, 2.75))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.38, 0.92, 1.0, 0.92], wspace=0.80)

    # A: composition
    ax = fig.add_subplot(gs[0, 0])
    vals = [100.0 * cnt[k] / n for k in order]
    colors = [C_STEEL, "#8aa8c0", C_TERRA, "#d9a877",
              C_GRAY, C_LGRAY, "#b0b7bd", "#cdd3d7"]
    ypos = np.arange(len(order))[::-1]
    ax.barh(ypos, vals, color=colors, height=0.68, zorder=3)
    ax.set_yticks(ypos)
    ax.set_yticklabels([labels[k] for k in order], fontsize=6.2)
    for y, v in zip(ypos, vals):
        ax.text(v + 0.7, y, f"{v:.1f}", va="center", fontsize=5.4, color="#333333")
    ax.set_xlabel("Edges (%)")
    ax.set_xlim(0, max(vals) * 1.30)
    style(ax, grid_axis="x")
    panel(ax, "A", dx=-0.62)

    # B: offsets
    ax = fig.add_subplot(gs[0, 1])
    d5 = [abs(int(r["d5"])) for r in near]
    d3 = [abs(int(r["d3"])) for r in near]
    mx = 6
    a = collections.Counter(min(x, mx) for x in d5)
    b = collections.Counter(min(x, mx) for x in d3)
    xs = list(range(mx + 1))
    ax.bar([x - 0.2 for x in xs], [max(a[x], 0.6) for x in xs], width=0.38,
           color=C_STEEL, label="|\u03945\u2032|", zorder=3)
    ax.bar([x + 0.2 for x in xs], [max(b[x], 0.6) for x in xs], width=0.38,
           color=C_TERRA, label="|\u03943\u2032|", zorder=3)
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels([str(x) for x in xs[:-1]] + [f"\u2265{mx}"])
    ax.set_xlabel("Terminal offset (nt)")
    ax.set_ylabel("Edges (log scale)")
    ax.legend(frameon=False, loc="upper right", handlelength=1.0)
    ax.set_title("93% differ by 1 nt", fontsize=5.4, color="#333333", pad=3)
    style(ax)
    panel(ax, "B", dx=-0.28)

    # C: abundance ratio
    ax = fig.add_subplot(gs[0, 2])
    lr = [np.log2(float(r["ratio"])) for r in star
          if float(r["ratio"]) > 0 and r["same_parent"] == "1"]
    lo = min(-6, np.floor(min(lr)))
    ax.hist(lr, bins=np.arange(lo, 0.4, 0.5), color=C_GRAY, zorder=3)
    ax.axvline(-1, ls="--", lw=0.8, color=C_RED, zorder=4)
    frac_below = 100.0 * np.mean(np.array(lr) <= -1)
    ax.set_xlim(lo - 0.3, 0.6)
    ax.text(0.04, 0.97, f"{frac_below:.0f}% >2-fold\nbelow rep.", transform=ax.transAxes,
            ha="left", va="top", fontsize=5.0, color=C_RED)
    ax.text(-1.15, 0.98, "2-fold", transform=ax.get_xaxis_transform(),
            ha="right", va="top", fontsize=5.0, color=C_RED, rotation=90)
    ax.set_xlabel("log$_2$(member / rep.)")
    ax.set_ylabel("Pairs")
    style(ax)
    panel(ax, "C", dx=-0.30)

    # D: within-cluster span
    ax = fig.add_subplot(gs[0, 3])
    spans = []
    for ms in load_clusters().values():
        if len(ms) < 2:
            continue
        P = [parse_id(m["tsRNA_id"]) for m in ms]
        spans.append(max(p[1] for p in P) - min(p[1] for p in P))
    ax.hist([min(x, 20) for x in spans], bins=range(0, 22), color=C_GRAY, zorder=3)
    ax.set_xticks([0, 5, 10, 15, 20])
    ax.set_xticklabels(["0", "5", "10", "15", "\u226520"])
    ax.set_xlabel("Cluster span (nt)")
    ax.set_ylabel("Clusters")
    n_big = sum(1 for x in spans if x >= 10)
    ax.text(0.97, 0.95, f"{n_big} clusters \u226510 nt\n({100.0 * n_big / len(spans):.1f}%)",
            transform=ax.transAxes, fontsize=5.2, color=C_RED, ha="right", va="top")
    style(ax)
    panel(ax, "D", dx=-0.30)

    fig.savefig(os.path.join(OUT, "FigS2_cluster_relationships.pdf"))
    fig.savefig(os.path.join(OUT, "_preview_FigS2.png"), dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------- Fig S3
def fig_singleton():
    rows = read_csv(os.path.join(BASE, "r1_diag", "cluster_overview.csv"))
    sizes = np.array([int(r["n_members"]) for r in rows])
    totals = np.array([int(float(r["cluster_total"])) for r in rows])
    n_cl = len(rows)
    n_single = int((sizes == 1).sum())

    fig = plt.figure(figsize=(7.1, 2.25))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.30)

    # A: cluster size distribution
    ax = fig.add_subplot(gs[0, 0])
    cnt = collections.Counter(np.minimum(sizes, 21))
    xs = sorted(cnt)
    cols = [C_RED if x == 1 else C_GRAY for x in xs]
    ax.bar(xs, [cnt[x] for x in xs], color=cols, width=0.82, zorder=3)
    ax.set_yscale("log")
    ax.set_xlabel("Cluster size (number of member fragments)")
    ax.set_ylabel("Clusters (log scale)")
    ax.set_xticks([1, 2, 3, 5, 10, 15, 21])
    ax.set_xticklabels(["1", "2", "3", "5", "10", "15", "\u226521"])
    ax.set_ylim(0.7, max(cnt.values()) * 3.2)
    ax.annotate(f"singletons: {n_single:,} / {n_cl:,}\n({100.0 * n_single / n_cl:.1f}%)",
                (1, cnt[1]), xytext=(2.6, cnt[1] * 0.42), fontsize=6.0, color=C_RED,
                arrowprops=dict(arrowstyle="->", color=C_RED, lw=0.7))
    style(ax)
    panel(ax, "A")

    # B: abundance by stratum
    ax = fig.add_subplot(gs[0, 1])
    strata = [("Singleton", sizes == 1), ("2\u20133", (sizes >= 2) & (sizes <= 3)),
              ("\u22654", sizes >= 4)]
    data = [totals[m] for _, m in strata]
    bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False,
                    medianprops=dict(color="white", lw=0.9),
                    whiskerprops=dict(color="#555555", lw=0.6),
                    capprops=dict(color="#555555", lw=0.6))
    for patch, c in zip(bp["boxes"], [C_RED, C_STEEL, C_DARK]):
        patch.set_facecolor(c)
        patch.set_edgecolor("#444444")
        patch.set_linewidth(0.5)
    ax.set_yscale("log")
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels([f"{s[0]}\n(n={(s[1]).sum():,})" for s in strata], fontsize=6.4)
    ax.set_ylabel("Cluster total count (log scale)")
    ax.set_xlabel("Cluster size stratum")
    med = [np.median(d) for d in data]
    for i, m in enumerate(med, 1):
        ax.text(i, m * 1.55, f"{m:,.0f}", ha="center", fontsize=5.8, color="#333333")
    style(ax)
    panel(ax, "B", dx=-0.26)

    fig.savefig(os.path.join(OUT, "FigS3_singleton_structure.pdf"))
    fig.savefig(os.path.join(OUT, "_preview_FigS3.png"), dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------- Fig S4
def fig_boundary_exclusion():
    full = {r["cluster_id"]: r for r in read_csv(
        os.path.join(BASE, "server_results", "results", "three_way", "de_cluster.csv"))}
    nb = {r["cluster_id"]: r for r in read_csv(
        os.path.join(BASE, "server_results", "results", "boundary_excl",
                     "de_cluster_no_boundary.csv"))}

    def pv(r):
        if r is None:
            return np.nan
        try:
            v = float(r["pvalue"])
        except (TypeError, ValueError):
            return np.nan
        return v if np.isfinite(v) and v > 0 else np.nan

    common = [k for k in full if k in nb and np.isfinite(pv(full[k])) and np.isfinite(pv(nb[k]))]
    x = -np.log10([pv(full[k]) for k in common])
    y = -np.log10([pv(nb[k]) for k in common])
    lost = [(a, b, k) for a, b, k in zip(x, y, common) if a > -np.log10(0.05) and b <= -np.log10(0.05)]
    kept = [(a, b) for a, b, k in zip(x, y, common) if not (a > -np.log10(0.05) and b <= -np.log10(0.05))]

    with open(os.path.join(BASE, "server_results", "results", "boundary_excl",
                           "summary.json")) as f:
        bs = json.load(f)

    fig = plt.figure(figsize=(7.1, 2.35))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.34)

    ax = fig.add_subplot(gs[0, 0])
    ax.scatter([k[0] for k in kept], [k[1] for k in kept], s=1.6, c=C_LGRAY,
               edgecolors="none", zorder=3, label="retained")
    ax.scatter([l[0] for l in lost], [l[1] for l in lost], s=4.2, c=C_RED,
               edgecolors="none", zorder=4, label=f"lost after exclusion (n={len(lost)})")
    lim = max(x) * 1.05
    ax.plot([0, lim], [0, lim], ls=":", lw=0.7, color="#777777", zorder=2)
    ax.axhline(-np.log10(0.05), ls="--", lw=0.7, color="#999999", zorder=2)
    ax.axvline(-np.log10(0.05), ls="--", lw=0.7, color="#999999", zorder=2)
    top = sorted(range(len(x)), key=lambda i: -x[i])[:2]
    for i in top:
        ax.scatter([x[i]], [y[i]], s=22, facecolor="white", edgecolor=C_RED, lw=0.9, zorder=6)
    ax.text(0.98, 0.46, "n = 21 of 26 lost calls\n(two strongest: 6.4e-7\u21920.40; 6.1e-6\u21920.36)",
            transform=ax.transAxes, ha="right", va="center", fontsize=5.4, color=C_RED)
    ax.set_xlabel("$-$log$_{10}$ p, all fragments")
    ax.set_ylabel("$-$log$_{10}$ p, boundary removed", labelpad=2.5)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.legend(frameon=False, loc="upper left", handlelength=0.9,
              markerscale=2.2, borderaxespad=0.2)
    style(ax, grid_axis=None)
    panel(ax, "A", dx=-0.16)

    ax = fig.add_subplot(gs[0, 1])
    groups = ["Fragment", "Cluster"]
    nom_full = [bs["de_fragment"]["full"]["nominal"], bs["de_cluster"]["full"]["nominal"]]
    nom_nb = [bs["de_fragment"]["no_boundary"]["nominal"],
              bs["de_cluster"]["no_boundary"]["nominal"]]
    fdr_full = [bs["de_fragment"]["full"]["fdr"], bs["de_cluster"]["full"]["fdr"]]
    fdr_nb = [bs["de_fragment"]["no_boundary"]["fdr"],
              bs["de_cluster"]["no_boundary"]["fdr"]]

    pos = np.arange(2)
    w = 0.34
    ax.bar(pos - w / 2, nom_full, width=w, color=C_STEEL, label="all fragments", zorder=3)
    ax.bar(pos + w / 2, nom_nb, width=w, color=C_TERRA, label="boundary removed", zorder=3)
    for i in range(2):
        ax.text(pos[i] - w / 2, nom_full[i] + 6, f"{nom_full[i]}", ha="center", fontsize=5.8)
        ax.text(pos[i] + w / 2, nom_nb[i] + 6, f"{nom_nb[i]}", ha="center", fontsize=5.8)
    ax.set_xticks(pos)
    ax.set_xticklabels(groups)
    ax.set_ylabel("Nominal DE calls (p < 0.05)")
    ax.set_ylim(0, max(nom_full) * 1.30)
    ax.text(0.03, 0.97,
            f"FDR < 0.05: fragment {fdr_full[0]}\u2192{fdr_nb[0]}, cluster {fdr_full[1]}\u2192{fdr_nb[1]}",
            transform=ax.transAxes, fontsize=5.2, color="#555555", va="top")
    leg = ax.legend(frameon=False, loc="center right", handlelength=1.0,
                    bbox_to_anchor=(0.99, 0.60), borderaxespad=0.0)
    ax.set_xlabel("Analysis level")
    style(ax)
    panel(ax, "B", dx=-0.26)

    fig.savefig(os.path.join(OUT, "FigS4_boundary_exclusion.pdf"))
    fig.savefig(os.path.join(OUT, "_preview_FigS4.png"), dpi=300)
    plt.close(fig)
    return len(lost)


# ---------------------------------------------------------------- Fig S5
def fig_threshold():
    with open(os.path.join(BASE, "server_results", "results", "energy_rep242",
                           "summary.json")) as f:
        sm = json.load(f)
    hp = read_csv(os.path.join(BASE, "server_results", "results", "energy_rep242",
                               "headline_pathways.csv"))

    ths = ["e25", "e29", "e33"]
    tlab = ["\u221225", "\u221229", "\u221233"]
    terms = ["Fanconi anemia", "Endometrial cancer", "Sulfur metabolism"]
    disp = ["Fanconi anaemia", "Endometrial cancer", "Sulfur metabolism"]
    pv = collections.defaultdict(dict)
    for r in hp:
        if r["p"] in ("", "None"):
            continue
        pv[r["term"]][r["threshold"]] = float(r["p"])

    fig = plt.figure(figsize=(7.1, 2.35))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.18, 1.0], wspace=0.34)

    # A: p across thresholds
    ax = fig.add_subplot(gs[0, 0])
    cols = [C_DARK, C_RED, C_TERRA]
    marks = ["o", "s", "^"]
    for t, d, c, m in zip(terms, disp, cols, marks):
        ys = [pv[t].get(x, np.nan) for x in ths]
        ax.plot(range(3), ys, m + "-", ms=3.6, lw=0.9, color=c, label=d, zorder=4)
    ax.axhline(0.05, ls="--", lw=0.8, color="#777777", zorder=2)
    ax.text(0.12, 0.062, "p = 0.05", fontsize=5.4, color="#555555",
            ha="left", va="bottom")
    ax.set_yscale("log")
    ax.set_xticks(range(3))
    ax.set_xticklabels([f"{l} kcal/mol" for l in tlab])
    ax.set_xlim(-0.35, 2.45)
    ax.set_ylim(0.0014, 0.9)
    ax.set_ylabel("Enrichr nominal p (log scale)")
    ax.legend(frameon=False, loc="lower left", handlelength=1.2,
              bbox_to_anchor=(0.01, 0.10), borderaxespad=0.0)
    ax.set_xlabel("RIsearch2 free-energy threshold")
    style(ax)
    panel(ax, "A", dx=-0.20)

    # B: stability metrics
    ax = fig.add_subplot(gs[0, 1])
    labels, gene_j, kegg_j, top10 = [], [], [], []
    for key, lab in (("e29_vs_e25", "\u221225 vs \u221229"), ("e29_vs_e33", "\u221233 vs \u221229")):
        c = sm["comparisons"][key]
        labels.append(lab)
        gene_j.append(c["top5"]["jaccard"])
        k = c["enrichment"]["KEGG_2021_Human"]
        kegg_j.append(k["jaccard_sig"])
        top10.append((k["top10_overlap"] / 10.0,
                      c["enrichment"]["GO_Biological_Process_2025"]["top10_overlap"] / 10.0))

    pos = np.arange(len(labels))
    w = 0.26
    ax.bar(pos - w, gene_j, width=w, color=C_STEEL, label="top-5% gene set", zorder=3)
    ax.bar(pos, kegg_j, width=w, color=C_TERRA, label="KEGG terms (p<0.05)", zorder=3)
    gbp = [t[1] for t in top10]
    ax.bar(pos + w, gbp, width=w, color=C_GRAY, label="GO BP top-10 overlap", zorder=3)
    for i in pos:
        for off, val in ((-w, gene_j[i]), (0, kegg_j[i]), (w, gbp[i])):
            ax.text(pos[i] + off, val + 0.02, f"{val:.2f}", ha="center", fontsize=5.4)
    ax.set_xticks(pos)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Jaccard index / overlap fraction")
    ax.legend(frameon=False, loc="upper right", ncol=1, handlelength=1.0,
              borderaxespad=0.2, fontsize=5.8)
    style(ax)
    panel(ax, "B", dx=-0.24)

    fig.savefig(os.path.join(OUT, "FigS5_threshold_sensitivity.pdf"))
    fig.savefig(os.path.join(OUT, "_preview_FigS5.png"), dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------- Fig S6
def fig_rep_members():
    rows = read_csv(os.path.join(BASE, "server_results", "results", "target_overlap",
                                 "per_cluster_overlap.csv"))
    ev = [r for r in rows
          if int(r["n_members"]) > 1 and int(r["rep_genes"]) > 0 and int(r["member_union_genes"]) > 0]
    jac = np.array([float(r["jaccard"]) for r in ev])
    with open(os.path.join(BASE, "server_results", "results", "target_overlap",
                           "summary.json")) as f:
        ts = json.load(f)
    g = ts["global_top5"]

    fig = plt.figure(figsize=(7.1, 2.15))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.32)

    ax = fig.add_subplot(gs[0, 0])
    ax.hist(jac, bins=np.linspace(0, 1, 26), color=C_GRAY, zorder=3)
    ax.axvline(np.median(jac), ls="--", lw=0.9, color=C_RED, zorder=4)
    ax.axvline(0.5, ls=":", lw=0.8, color="#777777", zorder=4)
    ax.set_xlabel("Jaccard index (representative vs union of members)")
    ax.set_ylabel(f"Clusters (n={len(ev)})")
    ax.text(np.median(jac) - 0.03, ax.get_ylim()[1] * 0.99,
            f"median {np.median(jac):.2f}", fontsize=6.0, color=C_RED,
            ha="right", va="top")
    ax.text(0.515, ax.get_ylim()[1] * 0.72, f"{100.0 * (jac < 0.5).mean():.1f}% < 0.5",
            fontsize=5.4, color="#555555", ha="left", rotation=90, va="center")
    style(ax)
    panel(ax, "A", dx=-0.16)

    ax = fig.add_subplot(gs[0, 1])
    vals = [g["rep_only"], g["overlap"], g["member_only"]]
    labs = ["Representative\nonly", "Shared", "Members\nonly"]
    cols = [C_TERRA, C_STEEL, C_DARK]
    bars = ax.bar(range(3), vals, color=cols, width=0.62, zorder=3)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 14, f"{v:,}", ha="center", fontsize=6.2)
    ax.set_xticks(range(3))
    ax.set_xticklabels(labs, fontsize=6.4)
    ax.set_ylabel("Genes in top 5% (CDS targets)")
    ax.set_ylim(0, max(vals) * 1.20)
    ax.text(0.5, 0.90, f"global Jaccard {g['jaccard']:.2f}", transform=ax.transAxes,
            ha="center", fontsize=6.2, color="#333333")
    style(ax)
    panel(ax, "B", dx=-0.26)

    fig.savefig(os.path.join(OUT, "FigS6_rep_vs_members.pdf"))
    fig.savefig(os.path.join(OUT, "_preview_FigS6.png"), dpi=300)
    plt.close(fig)


def main():
    cl = load_clusters()
    rlen = load_trna_lengths()
    print(f"clusters {len(cl):,}, tRNA refs {len(rlen):,}")
    fig_read_length(cl, rlen)
    print("FigS1 done")
    fig_relationships()
    print("FigS2 done")
    fig_singleton()
    print("FigS3 done")
    n_lost = fig_boundary_exclusion()
    print(f"FigS4 done (lost calls {n_lost})")
    fig_threshold()
    print("FigS5 done")
    fig_rep_members()
    print("FigS6 done")
    for f in sorted(os.listdir(OUT)):
        if f.endswith(".pdf"):
            print("  ", f, f"{os.path.getsize(os.path.join(OUT, f)):,} bytes")


if __name__ == "__main__":
    main()
