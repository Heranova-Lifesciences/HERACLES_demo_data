# -*- coding: utf-8 -*-
"""Figure 4A/B v2: corrected comparator marker matrix (R3-2 'Comparators').

Panel A (Xiong 2024, 5 markers) and B (Lu 2021, 3 markers).
Rows = markers, columns = HERACLES (non-clustering) / MINTmap / tDRmapper /
tsRNAsearch.  Cell = -log10(nominal p) with log2FC annotated; degenerate
all-zero-group fits hatched; absent rows / parental-tRNA-level NA marked.
Data: TableS17 (comparator rerun) + TableS10 (HERACLES non-clustering).
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T26 = pd.read_csv(os.path.join(
    ROOT, "output", "R1_回应包", "03_tables",
    "TableS26_comparator_marker_rerun.csv"))
T19 = pd.read_csv(os.path.join(
    ROOT, "output", "R1_回应包", "03_tables",
    "TableS19_marker_benchmark.csv"))

TOOLS = ["HERACLES (non-clust.)", "MINTmap", "tDRmapper", "tsRNAsearch"]

# HERACLES non-clustering nominal p / lfc per marker (TableS10 workflow rows)
HER = {}
for _, r in T19[T19.workflow == "non-clustering"].iterrows():
    HER[r.reported_marker] = (float(r.log2FoldChange), float(r.pvalue), "ok")

# per-dataset marker lists (display order)
PANELS = [
    ("A", "Xiong 2024",
     ["tiRNA-Lys-TTT-3-M2", "tiRNA-Val-CAC-2", "tiRNA-Lys-CTT-3",
      "tiRNA-Lys-CTT-1-M2", "tiRNA-Val-CAC-1-M3"]),
    ("B", "Lu 2021",
     ["tRF-1:28-Val-CAC-2", "tRF-55:76-Arg-ACG-1-M2", "tRF-1:24-Ser-CGA-1-M3"]),
]

TITLES = {"A": "Xiong et al. (Oncogene 2024)",
          "B": "Lu et al. (Front Mol Bio 2021)"}

# display-only: drop the positional prefix from Lu marker IDs in the figure
SHORT = {
    "tRF-1:28-Val-CAC-2": "tRF-Val-CAC-2",
    "tRF-55:76-Arg-ACG-1-M2": "tRF-Arg-ACG-1-M2",
    "tRF-1:24-Ser-CGA-1-M3": "tRF-Ser-CGA-1-M3",
}

TOOL_MAP = {"HERACLES (non-clust.)": "HERACLES (non-clust.)",
            "MINTmap": "MINTmap v1", "tDRmapper": "tDRmapper",
            "tsRNAsearch": "tsRNAsearch 1.02"}


def cell(tool, marker):
    """return (neglog10p, lfc, state) with state in ok/degenerate/absent/na"""
    if tool == "HERACLES (non-clust.)":
        lfc, p, _ = HER[marker]
        return -np.log10(p), lfc, "ok"
    sub = T26[(T26.tool == TOOL_MAP[tool]) & (T26.marker == marker)]
    if not len(sub):
        return None
    r = sub.iloc[0]
    if tool == "tsRNAsearch":
        return None, None, "na"
    if r["detection_corrected"].startswith("not detected"):
        return None, None, "absent"
    degen = str(r["degenerate_allzero_group"]) == "True"
    p, lfc = r["rerun_p"], r["rerun_log2FC"]
    if pd.isna(p):
        return None, None, "absent"
    if degen:
        return None, lfc, "degenerate"
    return -np.log10(p), lfc, "ok"


fig, axes = plt.subplots(2, 1, figsize=(8.6, 7.4),
                         gridspec_kw=dict(height_ratios=[5, 3], hspace=0.62))
for ax, (lab, ds, markers) in zip(axes, PANELS):
    nc, nr = len(TOOLS), len(markers)
    grid = np.full((nr, nc), np.nan)
    states = np.empty((nr, nc), dtype=object)
    for i, m in enumerate(markers):
        for j, t in enumerate(TOOLS):
            c = cell(t, m)
            if c is None:
                states[i, j] = "missing"
                continue
            nlp, lfc, st = c
            states[i, j] = st
            if nlp is not None:
                grid[i, j] = nlp
    im = ax.imshow(grid, cmap="Blues", vmin=0, vmax=6, aspect="auto")
    for i, m in enumerate(markers):
        for j, t in enumerate(TOOLS):
            st = states[i, j]
            x, y = j, i
            if st == "ok":
                _, lfc, _ = cell(t, m)
                ax.text(x, y, f"{lfc:+.2f}\np={10**-grid[i, j]:.3g}",
                        ha="center", va="center", fontsize=8.2, fontweight="bold",
                        color="white" if grid[i, j] > 3.0 else "#101010")
            elif st == "degenerate":
                ax.add_patch(Rectangle((x - .5, y - .5), 1, 1, fill=False,
                                       hatch="///", edgecolor="#37474f", lw=1.4))
                ax.text(x, y, "degenerate\nall-zero group", ha="center",
                        va="center", fontsize=7.4, color="#37474f",
                        bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2))
            elif st == "absent":
                ax.text(x, y, "absent\n(no exact\nsequence row)", ha="center",
                        va="center", fontsize=7.2, color="#555555",
                        style="italic", bbox=dict(fc="white", ec="none", pad=1.0))
            elif st == "na":
                ax.text(x, y, "NA\n(parental tRNA\nresolution)", ha="center",
                        va="center", fontsize=7.2, color="#555555",
                        style="italic", bbox=dict(fc="white", ec="none", pad=1.0))
    ax.set_xticks(range(nc))
    ax.set_xticklabels(TOOLS, fontsize=9.2)
    ax.set_yticks(range(nr))
    ax.set_yticklabels([SHORT.get(m, m) for m in markers], fontsize=8.8)
    ax.set_title(TITLES[lab], fontsize=11.5, loc="center", fontweight="bold", pad=9)
    for e in ["top", "right"]:
        ax.spines[e].set_visible(False)
fig.text(0.98, 0.018,
         "colour = $-\\log_{10}$(nominal $p$);  hatched = degenerate fit "
         "(entire control group zero counts)", fontsize=8.4,
         ha="right", va="bottom", color="#333333")
fig.subplots_adjust(left=0.28, right=0.98, top=0.925, bottom=0.10)
out = os.path.join(ROOT, "figures", "Fig4AB_comparator_rerun_v3")
fig.savefig(out + ".pdf")
fig.savefig(out + ".png", dpi=600)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
for a_, tag in ((axes[0], "A"), (axes[1], "B")):
    bb = a_.get_tightbbox(renderer).transformed(fig.dpi_scale_trans.inverted())
    bb = bb.expanded(1.02, 1.06)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(ROOT, "figures", f"Fig4_panel{tag}.{ext}"),
                    dpi=600, bbox_inches=bb, pad_inches=0.03)
print("saved", out + ".pdf/.png + separate A/B panels")
