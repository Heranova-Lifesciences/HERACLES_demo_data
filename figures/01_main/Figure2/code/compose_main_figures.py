# -*- coding: utf-8 -*-
"""Compose revised main figures (v2, fixed crop boxes measured from preview).
Fig2 = original A + original B + NEW Fig2C (old C/D removed)
Fig4 = rerun Fig4AB heatmaps + NEW Fig4C marker reanalysis
"""
import os
import numpy as np
import pymupdf
from PIL import Image, ImageDraw, ImageFont

SRC = r"C:/Users/huyao/Desktop/figure"
FIG = r"C:/Users/huyao/Desktop/HERACLES修改/figures"
OUT = os.path.join(FIG, "主图重拼")
os.makedirs(OUT, exist_ok=True)
DPI = 600
PT2PX = DPI / 72.0
FONT_BOLD = r"C:/Windows/Fonts/arialbd.ttf"


def render_pt(pdf, dpi=DPI):
    d = pymupdf.open(pdf)
    pix = d[0].get_pixmap(dpi=dpi)
    img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    return img, img.size[0] / d[0].rect.width  # px per pt


def add_letter(canvas, letter, xy, size_pt=20):
    f = ImageFont.truetype(FONT_BOLD, int(size_pt * PT2PX))
    ImageDraw.Draw(canvas).text(xy, letter, font=f, fill=(0, 0, 0))


def save(img, name):
    img.save(os.path.join(OUT, name + ".png"), dpi=(DPI, DPI))
    # vector-wrapper PDF: page sized to image at 72 dpi scale
    w_pt, h_pt = img.size[0] / PT2PX, img.size[1] / PT2PX
    doc = pymupdf.open()
    pg = doc.new_page(width=w_pt, height=h_pt)
    tmp = os.path.join(OUT, "_tmp.png")
    img.save(tmp)
    pg.insert_image(pymupdf.Rect(0, 0, w_pt, h_pt), filename=tmp)
    doc.save(os.path.join(OUT, name + ".pdf"))
    doc.close()
    os.remove(tmp)
    print("saved", name, img.size)


# ---------------------------------------------------------------- Fig 2
base, k = render_pt(os.path.join(SRC, "fig2.pdf"))  # k px per pt
W, H = base.size

# measured from 150dpi preview (display 836x1030, factor 1.525 to 150dpi):
# A: y 6..228 pt full width | B: x 14..302 pt, y 230..412 pt
a_bot = int(229 * k)
A = base.crop((0, 0, W, a_bot))

b_x0, b_x1 = 0, int(304 * k)
b_y0, b_y1 = int(230 * k), int(436 * k)
B = base.crop((b_x0, b_y0, b_x1, b_y1))

# new C box: x 316..604 pt, vertically centered in bottom row
cx0, cx1 = int(316 * k), int(604 * k)
row_h = max(B.size[1], (b_y1 - b_y0))
newC = Image.open(os.path.join(FIG, "正文主图", "Fig2C_paired_pvalues.png")).convert("RGB")
box_w, box_h = cx1 - cx0, b_y1 - b_y0
scale = min(box_w / newC.size[0], box_h / newC.size[1])
newC2 = newC.resize((int(newC.size[0] * scale), int(newC.size[1] * scale)),
                    Image.LANCZOS)

canvas = Image.new("RGB", (W, a_bot + int(10 * k) + row_h + int(8 * k)), "white")
canvas.paste(A, (0, 0))
y_bot = a_bot + int(10 * k)
canvas.paste(B, (b_x0, y_bot))
px = cx0 + (box_w - newC2.size[0]) // 2
py = y_bot + (row_h - newC2.size[1]) // 2
canvas.paste(newC2, (px, py))
add_letter(canvas, "C", (cx0, y_bot - int(2 * k)), size_pt=18)
save(canvas, "Figure2_revised")

# ---------------------------------------------------------------- Fig 4
def trim(img, tol=250):
    a = np.asarray(img.convert("L"))
    rows_ = ~(a > tol).all(axis=1)
    cols_ = ~(a > tol).all(axis=0)
    r0 = int(np.argmax(rows_)); r1 = len(rows_) - int(np.argmax(rows_[::-1]))
    c0 = int(np.argmax(cols_)); c1 = len(cols_) - int(np.argmax(cols_[::-1]))
    return img.crop((c0, r0, c1, r1))


ab = trim(Image.open(os.path.join(FIG, "Fig4AB_comparator_rerun.png")).convert("RGB"))
c4 = trim(Image.open(os.path.join(FIG, "正文主图", "Fig4C_marker_reanalysis.png")).convert("RGB"))

target_h = 2400
ab2 = ab.resize((int(ab.size[0] * target_h / ab.size[1]), target_h), Image.LANCZOS)
c42 = c4.resize((int(c4.size[0] * target_h / c4.size[1]), target_h), Image.LANCZOS)

gap = int(0.06 * target_h)
marg = int(0.03 * target_h)
Wc = marg + ab2.size[0] + gap + c42.size[0] + marg
canvas4 = Image.new("RGB", (Wc, target_h + 2 * marg), "white")
canvas4.paste(ab2, (marg, marg))
canvas4.paste(c42, (marg + ab2.size[0] + gap, marg))
add_letter(canvas4, "C", (marg + ab2.size[0] + gap, marg - int(0.010 * target_h)),
           size_pt=18)
save(canvas4, "Figure4_revised")
