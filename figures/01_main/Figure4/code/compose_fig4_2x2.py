# -*- coding: utf-8 -*-
"""Figure 4 final layout: 2x2  (A | B on top, C | D below).
A/B = restyled comparator heatmaps; C = slope graph (no red highlight);
D = P / BH bar panel, now labelled D.
"""
import os
import numpy as np
import pymupdf
from PIL import Image, ImageDraw, ImageFont

FIG = r"C:/Users/huyao/Desktop/HERACLES修改/figures"
OUT = os.path.join(FIG, "主图重拼")
os.makedirs(OUT, exist_ok=True)
FONT_BOLD = r"C:/Windows/Fonts/arialbd.ttf"
FONT = r"C:/Windows/Fonts/arial.ttf"

P = {
    "A": os.path.join(FIG, "Fig4_panelA.png"),
    "B": os.path.join(FIG, "Fig4_panelB.png"),
    "C": os.path.join(FIG, "正文主图", "Fig4_panelC.png"),
    "D": os.path.join(FIG, "正文主图", "Fig4_panelD.png"),
}


def trim(img, tol=250):
    a = np.asarray(img.convert("L"))
    rows_ = ~(a > tol).all(axis=1)
    cols_ = ~(a > tol).all(axis=0)
    r0 = int(np.argmax(rows_)); r1 = len(rows_) - int(np.argmax(rows_[::-1]))
    c0 = int(np.argmax(cols_)); c1 = len(cols_) - int(np.argmax(cols_[::-1]))
    return img.crop((c0, r0, c1, r1))


panels = {k: trim(Image.open(v).convert("RGB")) for k, v in P.items()}

CELL = 3600                      # common cell width
gapx = int(0.045 * CELL)
marg = int(0.025 * CELL)
head = int(0.065 * CELL)         # room above each row for the panel letters
gapy = int(0.075 * CELL)
foot = int(0.030 * CELL)

scaled = {}
for k, im in panels.items():
    h = int(im.size[1] * CELL / im.size[0])
    scaled[k] = im.resize((CELL, h), Image.LANCZOS)

row1_h = max(scaled["A"].size[1], scaled["B"].size[1])
row2_h = max(scaled["C"].size[1], scaled["D"].size[1])
W = marg + CELL + gapx + CELL + marg
y1 = head
y2 = y1 + row1_h + gapy + head
H = y2 + row2_h + foot

canvas = Image.new("RGB", (W, H), "white")
pos = {"A": (marg, y1), "B": (marg + CELL + gapx, y1),
       "C": (marg, y2), "D": (marg + CELL + gapx, y2)}
for k, (x, y) in pos.items():
    canvas.paste(scaled[k], (x, y))

d = ImageDraw.Draw(canvas)
f_letter = ImageFont.truetype(FONT_BOLD, int(0.031 * CELL))
for k, (x, y) in pos.items():
    d.text((x, y - int(0.058 * CELL)), k, font=f_letter, fill=(12, 12, 12))

f_note = ImageFont.truetype(FONT, int(0.012 * CELL))
note = ("colour = \u2212log10(nominal P);  hatched = degenerate fit "
        "(entire control group zero counts)")
# 作为 A/B 热图的共用脚注，放在第一行下面
d.text((marg, y1 + row1_h + int(0.020 * CELL)), note, font=f_note,
       fill=(70, 70, 70))

name = "Figure4_revised_v3"
canvas.save(os.path.join(OUT, name + ".png"), dpi=(600, 600))

w_pt, h_pt = W / (600 / 72), H / (600 / 72)
tmp = os.path.join(OUT, "_tmp.jpg")
canvas.save(tmp, quality=92)
doc = pymupdf.open(); pg = doc.new_page(width=w_pt, height=h_pt)
pg.insert_image(pymupdf.Rect(0, 0, w_pt, h_pt), filename=tmp)
doc.save(os.path.join(OUT, name + ".pdf"), deflate=True, garbage=4)
doc.close(); os.remove(tmp)
print("saved", name, canvas.size)
for k in "ABCD":
    print(" ", k, panels[k].size, "->", scaled[k].size)
