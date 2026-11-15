# -*- coding: utf-8 -*-
"""Figure 3 rebuilt in the style of the reference figure (Desktop/fig3.pdf).

Panel A  Flowchart copied from the reference: rounded 242 DE tsRNAs box, two
         dashed arrows to the CDS / 3'UTR references, filtering-criteria
         cartoons, pink 95 / 691 Targets boxes and dashed arrows down to B/C.
Panel B/C  Identical bubble panels: x = GSEA Combined Score, colour =
         -log10(pvalue), marker size = Combined Score, colour bar and size
         legend at the right of each panel.

Data are our own Enrichr results (CDS: energy_rep242 e29; 3'UTR:
energy_rep242/enrichment_3UTR), i.e. the seven nominal-P<0.05 terms per arm
ordered by Combined Score.

Writes figures/正文主图/Figure3_revised.{png,pdf} and Fig3_panelA_workflow.*
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "figures", "正文主图")
CDS = os.path.join(BASE, "server_results", "results", "energy_rep242",
                   "enrichment", "e29_KEGG_2021_Human.csv")
UTR = os.path.join(BASE, "server_results", "results", "energy_rep242",
                   "enrichment_3UTR", "KEGG_enrichment.csv")

INK = "#2f2f2f"
PINK_FC, PINK_EC = "#f7dbe8", "#b97aa0"
TOP_FC, TOP_EC = "#e9edfa", "#7c8fc0"
CMAP = LinearSegmentedColormap.from_list(
    "pval", ["#463a99", "#7747ad", "#bc3b78", "#e52b46"])

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.2,
    "axes.labelsize": 9.0, "xtick.labelsize": 8.0, "ytick.labelsize": 8.0,
    "pdf.fonttype": 42, "savefig.pad_inches": 0.04,
})


def rbox(ax, x, y, w, h, text="", fc="white", ec=INK, lw=0.9, fs=8.2,
         bold=False, round_=0.02, color=INK):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle=f"round,pad=0,rounding_size={round_}",
                                fc=fc, ec=ec, lw=lw, zorder=3))
    if text:
        ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=color,
                fontweight="bold" if bold else "normal", zorder=4)


def darrow(ax, p0, p1, lw=1.0, color="#555555"):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=9,
                                 lw=lw, color=color, linestyle=(0, (3, 3)),
                                 zorder=2, shrinkA=2, shrinkB=2))


def cartoons(ax):
    # free-energy cartoon: a row of short bars with a red/orange energy trace
    cols = ["#5b9bd5", "#ed7d31", "#4472c4", "#ed7d31", "#5b9bd5"]
    for i in range(11):
        ax.plot([0.365 + i * 0.007, 0.365 + i * 0.007], [0.505, 0.525],
                color=cols[i % len(cols)], lw=1.5, solid_capstyle="round")
    ax.plot([0.362, 0.375, 0.390, 0.405, 0.420, 0.435],
            [0.548, 0.532, 0.558, 0.536, 0.550, 0.535],
            color="#e0662d", lw=1.2)
    ax.text(0.418, 0.575, "\u0394G", fontsize=6.2, color="#777777", style="italic")
    ax.add_patch(Circle((0.443, 0.578), 0.005, fill=False, ec="#999999", lw=0.7))
    ax.text(0.400, 0.425, "Free energy", ha="center", fontsize=7.0, style="italic", color=INK)
    ax.text(0.400, 0.395, "\u2264-29 kcal/mol", ha="center", fontsize=6.8,
            style="italic", color=INK)

    # node-connectivity cartoon: coloured nodes around a central dark box
    nodes = [(0.575, 0.585, "#c8a2a2"), (0.605, 0.595, "#8fa8d0"),
             (0.635, 0.580, "#d8cf8a"), (0.565, 0.545, "#d98f8f"),
             (0.645, 0.545, "#9fc8a0"), (0.605, 0.555, "#a9b7d8")]
    for x, y, c in nodes:
        ax.add_patch(Circle((x, y), 0.007, fc=c, ec="none", zorder=3))
        ax.plot([x, 0.605], [y, 0.520], color="#c9c9c9", lw=0.6, zorder=2)
    rbox(ax, 0.605, 0.520, 0.048, 0.022, "tsRNA", fc="#3b4f86", ec="none",
         fs=4.8, color="white", round_=0.006)
    ax.text(0.605, 0.425, "Node connectivity", ha="center", fontsize=7.0,
            style="italic", color=INK)
    ax.text(0.605, 0.395, "top 5%", ha="center", fontsize=6.8,
            style="italic", color=INK)


def panel_a(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.012, 0.95, "A", fontsize=15, fontweight="bold", color=INK,
            ha="left", va="top")

    rbox(ax, 0.5, 0.90, 0.19, 0.105, "242 DE tsRNAs", fc=TOP_FC, ec=TOP_EC,
         fs=8.6, bold=True, round_=0.022)
    darrow(ax, (0.425, 0.862), (0.255, 0.665))
    darrow(ax, (0.575, 0.862), (0.745, 0.665))
    ax.text(0.270, 0.775, "Human CDS reference", ha="center", fontsize=7.0, color=INK)
    ax.text(0.730, 0.775, "Human 3UTR reference", ha="center", fontsize=7.0, color=INK)
    ax.text(0.5, 0.775, "Filtering Criteria", ha="center", fontsize=7.0,
            style="italic", color=INK)

    cartoons(ax)

    rbox(ax, 0.225, 0.215, 0.155, 0.105, "95 Targets", fc=PINK_FC, ec=PINK_EC,
         fs=8.4, bold=True, round_=0.028)
    rbox(ax, 0.775, 0.215, 0.155, 0.105, "691 Targets", fc=PINK_FC, ec=PINK_EC,
         fs=8.4, bold=True, round_=0.028)
    darrow(ax, (0.225, 0.155), (0.225, 0.035))
    darrow(ax, (0.775, 0.155), (0.775, 0.035))


def bubble_panel(fig, ax, csv_path, letter, title, cb_ax, leg_ax):
    d = pd.read_csv(csv_path)
    d = d[d["P-value"] < 0.05].nlargest(7, "Combined Score").iloc[::-1]
    x = d["Combined Score"].astype(float).values
    c = -np.log10(d["P-value"].astype(float).values)
    y = np.arange(len(d))
    smax = max(x.max(), 1)
    sizes = 18 + 190 * (x / smax) ** 1.2

    sc = ax.scatter(x, y, s=sizes, c=c, cmap=CMAP, vmin=c.min(), vmax=c.max(),
                    edgecolors="none", zorder=3)
    ax.set_yticks(y)
    labels = [t if len(t) <= 26 else t[:25] + "\u2026" for t in d["Term"].values]
    ax.set_yticklabels(labels, fontsize=7.2)
    ax.tick_params(axis="y", length=3, width=0.7, pad=4)
    ax.set_ylim(-0.6, len(d) - 0.4)
    ax.set_xlim(0, x.max() * 1.14)
    ax.set_xlabel("Enrichr combined score", fontsize=9.0, labelpad=6)
    ax.set_title(title, fontsize=9.6, fontweight="bold", pad=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.grid(False)
    ax.text(-0.12, 1.10, letter, transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="top", ha="left")

    cb = fig.colorbar(sc, cax=cb_ax)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=6.6, length=2)
    cb_ax.set_title("-log10(nominal P)", fontsize=7.0, pad=13, loc="left")

    leg_ax.axis("off")
    leg_ax.set_title("Combined Score", fontsize=7.0, pad=6, loc="left")
    ticks = [round(v) for v in np.linspace(max(1, x.min() * 0.8), x.max(), 4)]
    for i, v in enumerate(ticks):
        yy = 0.80 - i * 0.235
        sz = 18 + 190 * (v / smax) ** 1.2
        leg_ax.scatter([0.22], [yy], s=sz, c="black", transform=leg_ax.transAxes,
                       clip_on=False)
        leg_ax.text(0.55, yy, str(v), transform=leg_ax.transAxes, va="center",
                    fontsize=6.6)


def main():
    os.makedirs(OUT, exist_ok=True)
    # A alone
    f1 = plt.figure(figsize=(7.0, 2.6))
    panel_a(f1.add_axes([0, 0, 1, 1]))
    for ext in ("pdf", "png"):
        f1.savefig(os.path.join(OUT, f"Fig3_panelA_workflow.{ext}"), dpi=400)
    plt.close(f1)

    # full figure
    fig = plt.figure(figsize=(10.5, 6.2))
    axA = fig.add_axes([0.02, 0.635, 0.96, 0.345])
    axB = fig.add_axes([0.150, 0.115, 0.255, 0.455])
    axC = fig.add_axes([0.645, 0.115, 0.245, 0.455])
    cbB = fig.add_axes([0.415, 0.330, 0.013, 0.175])
    cbC = fig.add_axes([0.900, 0.330, 0.013, 0.175])
    lgB = fig.add_axes([0.405, 0.125, 0.042, 0.185])
    lgC = fig.add_axes([0.890, 0.125, 0.042, 0.185])
    panel_a(axA)
    bubble_panel(fig, axB, CDS, "B", "KEGG over-representation (CDS)", cbB, lgB)
    bubble_panel(fig, axC, UTR, "C", "KEGG over-representation (3\u2032 UTR)", cbC, lgC)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"Figure3_revised.{ext}"), dpi=400)
    plt.close(fig)
    print("saved Figure3_revised.png/.pdf (reference style)")


if __name__ == "__main__":
    main()
