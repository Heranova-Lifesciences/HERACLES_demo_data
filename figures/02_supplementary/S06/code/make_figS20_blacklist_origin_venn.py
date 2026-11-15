#!/usr/bin/env python3
"""
FigS20：程序自带 CDS 黑名单（1,377 基因）的来源（三队列预测靶基因三方交集 Venn）。

三个队列各自取 top-10% 预测 tsRNA 连通度基因（2026-06-12，
overlap/enrichment/{our_data,paper1,paper2}.txt，2,014 / 2,010 / 2,018 个），
三方交集恰好 1,377 个 = overlap_genes.txt，原样拷为
RIsearch2_index/CDS_black_list.txt（md5 一致）。

风格与 FigS19 一致：经典固定布局三圆、同色橙 alpha 叠加。

区域数（本地校验过，ABC 与黑名单逐基因相等）：
  A only 424 | B only 166 | C only 157
  AB 98 | AC 115 | BC 369 | ABC 1,377

输出：figures/FigS20_blacklist_origin_venn.pdf
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "figures")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "pdf.fonttype": 42,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
C_ORANGE = "#c08552"
ALPHA = 0.38
C_EDGE = "#7a7f84"
C_DARK = "#3f3f3f"

R = 1.0
A = (0.0, 0.60)      # in-house（上）
B = (-0.60, -0.35)   # Xiong（左下）
C = (0.60, -0.35)    # Lu（右下）

# 区域 (in_A, in_B, in_C) -> (数量, 标签位置)
REGIONS = {
    (1, 0, 0): (424, (0.0, 1.05)),
    (0, 1, 0): (166, (-1.18, -0.62)),
    (0, 0, 1): (157, (1.18, -0.62)),
    (1, 1, 0): (98, (-0.52, 0.42)),
    (1, 0, 1): (115, (0.52, 0.42)),
    (0, 1, 1): (369, (0.0, -0.92)),
    (1, 1, 1): (1377, (0.0, -0.12)),
}
N_ALL3 = REGIONS[(1, 1, 1)][0]


def main():
    fig, ax = plt.subplots(figsize=(4.6, 4.0))

    # 三层同色 alpha 叠加：颜色深浅 = 支持队列数
    for ctr in (A, B, C):
        ax.add_patch(mpatches.Circle(ctr, R, facecolor=C_ORANGE, alpha=ALPHA,
                                     edgecolor="none", zorder=2))
    for ctr in (A, B, C):
        ax.add_patch(mpatches.Circle(ctr, R, fill=False, ec=C_EDGE, lw=0.7,
                                     zorder=4))

    for (a, b, c), (n, (x, y)) in REGIONS.items():
        k = a + b + c
        ax.text(x, y, f"{n:,}", ha="center", va="center", fontsize=7.0,
                zorder=5, color="white" if k == 3 else C_DARK)

    # 集合标签
    ax.text(A[0], A[1] + R + 0.22, "In-house cohort\n(n = 2,014)", ha="center",
            fontsize=6.8, color=C_DARK)
    ax.text(B[0] - 1.05, B[1] - 0.92, "Xiong et al.\n(n = 2,010)", ha="center",
            fontsize=6.8, color=C_DARK)
    ax.text(C[0] + 1.05, C[1] - 0.92, "Lu et al.\n(n = 2,018)", ha="center",
            fontsize=6.8, color=C_DARK)

    # 图例：按合成色画 swatch（1/2/3 层 alpha 在白底上的合成结果）
    comp = {k: 1 - (1 - ALPHA) ** k for k in (1, 2, 3)}
    handles = [
        mpatches.Patch(fc=C_ORANGE, alpha=comp[1], ec=C_EDGE, lw=0.4,
                       label="target in 1 cohort"),
        mpatches.Patch(fc=C_ORANGE, alpha=comp[2], ec=C_EDGE, lw=0.4,
                       label="target in 2 cohorts"),
        mpatches.Patch(fc=C_ORANGE, alpha=comp[3], ec=C_EDGE, lw=0.4,
                       label="target in all 3"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=6.0, loc="upper left",
              bbox_to_anchor=(-0.08, 1.0), handlelength=1.2, borderaxespad=0)

    ax.text(0, -1.86,
            f"{N_ALL3:,} genes predicted as tsRNA targets in all 3 cohorts "
            "\u2192 CDS blacklist",
            ha="center", fontsize=6.8, color=C_DARK)
    ax.text(0, -2.14,
            "top-10% predicted tsRNA connectivity per cohort "
            "(RIsearch2, \u221229 kcal/mol)",
            ha="center", fontsize=5.6, color="#8a8f94")

    ax.set_xlim(-2.15, 2.15)
    ax.set_ylim(-2.25, 1.95)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.savefig(os.path.join(OUT, "FigS20_blacklist_origin_venn.pdf"))
    fig.savefig(os.path.join(OUT, "_build_png", "_preview_FigS20.png"), dpi=300)
    plt.close(fig)
    print("saved")


if __name__ == "__main__":
    main()
