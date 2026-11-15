# -*- coding: utf-8 -*-
"""Rebuild Figure 1 (HERACLES workflow overview) as a true vector PDF.
Layout and colours replicate the original bitmap version (Desktop/figure/fig1.pdf).
"""
import math
import os
import pymupdf

OUT = r"C:/Users/huyao/Desktop/HERACLES修改/figures/主图重拼"
os.makedirs(OUT, exist_ok=True)

W, H = 500, 760
NAVY = (0.114, 0.239, 0.42)
TXT = (0.10, 0.10, 0.12)
DESC = (0.36, 0.38, 0.43)

doc = pymupdf.open()
page = doc.new_page(width=W, height=H)
shape = page.new_shape()


def rr(rect, fill, stroke, width=1.3, radius_pt=9):
    frac = min(0.5, radius_pt / min(rect.width, rect.height))
    shape.draw_rect(rect, radius=frac)
    shape.finish(fill=fill, color=stroke, width=width, closePath=True)


# 文字统一缓存，等所有图形 commit 之后再画，保证文字压在色块/图标之上
TEXT_BUFFER = []


def ctext(pt, s, size, font="helv", color=TXT):
    w = pymupdf.get_text_length(s, fontname=font, fontsize=size)
    TEXT_BUFFER.append(((pt[0] - w / 2, pt[1]), s, size, font, color))


def ltext(pt, s, size, font="helv", color=TXT):
    TEXT_BUFFER.append((pt, s, size, font, color))


def arrow(points):
    shape.draw_polyline(points)
    shape.finish(color=NAVY, width=2.2, closePath=False)
    (x1, y1), (x2, y2) = points[-2], points[-1]
    ang = math.atan2(y2 - y1, x2 - x1)
    L = 8.0
    p1 = (x2 - L * math.cos(ang - 0.42), y2 - L * math.sin(ang - 0.42))
    p2 = (x2 - L * math.cos(ang + 0.42), y2 - L * math.sin(ang + 0.42))
    shape.draw_polyline([p1, (x2, y2), p2])
    shape.finish(fill=NAVY, color=NAVY, width=0.5, closePath=True)


def elbow(points):
    shape.draw_polyline(points)
    shape.finish(color=NAVY, width=2.2, closePath=False)


def head_at(pt, ang):
    x2, y2 = pt
    L = 8.0
    p1 = (x2 - L * math.cos(ang - 0.42), y2 - L * math.sin(ang - 0.42))
    p2 = (x2 - L * math.cos(ang + 0.42), y2 - L * math.sin(ang + 0.42))
    shape.draw_polyline([p1, (x2, y2), p2])
    shape.finish(fill=NAVY, color=NAVY, width=0.5, closePath=True)


def circle(center, r, fill, text=None, tcol=(1, 1, 1), size=9.0):
    shape.draw_circle(center, r)
    shape.finish(fill=fill, color=fill, width=0)
    if text:
        ctext((center[0], center[1] + size * 0.36), text, size, "hebo", tcol)


# ------------------------------------------------------------- icons ----
def icon_doc(c):
    x, y = c[0] - 12, c[1] - 16
    shape.draw_rect(pymupdf.Rect(x, y, x + 22, y + 28))
    shape.finish(fill=(1, 1, 1), color=(0.35, 0.5, 0.75), width=1.1)
    shape.draw_polyline([(x + 14, y), (x + 14, y + 8), (x + 22, y + 8)])
    shape.finish(fill=(1, 1, 1), color=(0.35, 0.5, 0.75), width=1.1, closePath=False)
    for dy in (5, 10, 15):
        shape.draw_line((x + 4, y + dy), (x + 18, y + dy))
        shape.finish(color=(0.55, 0.65, 0.8), width=1.2)
    shape.draw_rect(pymupdf.Rect(x - 2, y + 24, x + 18, y + 33))
    shape.finish(fill=(0.23, 0.42, 0.72), color=(0.23, 0.42, 0.72), width=0)
    ctext((x + 8, y + 31), "FASTQ", 5.4, "hebo", (1, 1, 1))


def icon_scissors(c):
    x, y = c
    for dx in (-6, 6):
        shape.draw_circle((x + dx, y + 8), 3.4)
        shape.finish(fill=(0.16, 0.55, 0.45), color=(0.16, 0.55, 0.45), width=0)
    shape.draw_polyline([(x - 6, y + 4), (x + 9, y - 10)])
    shape.finish(color=(0.55, 0.6, 0.65), width=2.6)
    shape.draw_polyline([(x + 6, y + 4), (x - 9, y - 10)])
    shape.finish(color=(0.55, 0.6, 0.65), width=2.6)
    for i, dx in enumerate((-3, 1, 5)):
        shape.draw_line((x + dx, y - 14), (x + dx, y - 9 - (i % 2) * 2))
        shape.finish(color=(0.45, 0.72, 0.62), width=1.4)


def icon_table(c):
    x, y = c[0] - 13, c[1] - 12
    shape.draw_rect(pymupdf.Rect(x, y, x + 26, y + 20))
    shape.finish(fill=(1, 1, 1), color=(0.62, 0.7, 0.42), width=1.1)
    for dy in (5, 10, 15):
        shape.draw_line((x, y + dy), (x + 26, y + dy))
        shape.finish(color=(0.62, 0.7, 0.42), width=0.8)
    for dx in (7, 14, 21):
        shape.draw_line((x + dx, y), (x + dx, y + 20))
        shape.finish(color=(0.62, 0.7, 0.42), width=0.8)
    shape.draw_rect(pymupdf.Rect(x + 14, y + 4, x + 20, y + 10))
    shape.finish(fill=(0.85, 0.9, 0.68), color=(0.85, 0.9, 0.68), width=0)
    shape.draw_circle((x + 24, y + 20), 5.4)
    shape.finish(fill=(0.72, 0.78, 0.45), color=(0.72, 0.78, 0.45), width=0)
    ctext((x + 24, y + 22.6), "\u03a3", 7.5, "hebo", (1, 1, 1))


def icon_branch(c):
    x, y = c
    pts = [(x - 10, y - 9), (x + 10, y - 9), (x, y), (x - 8, y + 10), (x + 8, y + 10)]
    cols = [(0.85, 0.45, 0.25), (0.3, 0.55, 0.8), (0.4, 0.42, 0.75),
            (0.9, 0.72, 0.3), (0.5, 0.75, 0.5)]
    for i, a in enumerate(pts):
        for b in pts[i + 1:]:
            shape.draw_line(a, b)
            shape.finish(color=(0.6, 0.65, 0.7), width=1.0)
    for p, col in zip(pts, cols):
        shape.draw_circle(p, 3.6)
        shape.finish(fill=col, color=col, width=0)


def icon_align(c, teal=False):
    if not teal:
        bars = [(-14, 10, 24, (0.85, 0.35, 0.3)), (-8, 2, 22, (0.3, 0.5, 0.8)),
                (-14, -6, 16, (0.55, 0.4, 0.75)), (-2, -14, 14, (0.85, 0.55, 0.3))]
    else:
        bars = [(-14, 6, 26, (0.55, 0.78, 0.76)), (-14, -2, 20, (0.3, 0.62, 0.6)),
                (-8, -10, 22, (0.75, 0.88, 0.87))]
    for x0, dy, ln, col in bars:
        shape.draw_rect(pymupdf.Rect(c[0] + x0, c[1] + dy - 2.2,
                                     c[0] + x0 + ln, c[1] + dy + 2.2))
        shape.finish(fill=col, color=col, width=0)


def icon_netgraph(c):
    x, y = c
    nodes = [(0, -11), (10, -4), (7, 9), (-7, 9), (-10, -4), (0, 1)]
    for i, a in enumerate(nodes):
        for b in nodes[i + 1:]:
            shape.draw_line((x + a[0], y + a[1]), (x + b[0], y + b[1]))
            shape.finish(color=(0.62, 0.5, 0.82), width=0.9)
    for a in nodes:
        shape.draw_circle((x + a[0], y + a[1]), 3.2)
        shape.finish(fill=(0.45, 0.32, 0.68), color=(0.45, 0.32, 0.68), width=0)


def icon_volcano(c):
    x, y = c
    shape.draw_polyline([(x - 13, y + 10), (x + 13, y + 10)])
    shape.finish(color=(0.45, 0.45, 0.5), width=1.1)
    shape.draw_polyline([(x - 13, y + 10), (x - 13, y - 10)])
    shape.finish(color=(0.45, 0.45, 0.5), width=1.1)
    blue = [(-9, 6), (-6, 3), (-10, 1), (-4, 7), (-7, -2), (-3, 4), (-11, -4)]
    red = [(4, 7), (6, 4), (8, 8), (5, -1), (9, 3), (3, 9), (10, 7)]
    for dx, dy in blue:
        shape.draw_circle((x + dx, y + 6 - dy), 1.7)
        shape.finish(fill=(0.3, 0.5, 0.8), color=(0.3, 0.5, 0.8), width=0)
    for dx, dy in red:
        shape.draw_circle((x + dx, y + 6 - dy), 1.7)
        shape.finish(fill=(0.85, 0.3, 0.3), color=(0.85, 0.3, 0.3), width=0)


def icon_target(c):
    x, y = c
    for r, col in ((11, (0.88, 0.42, 0.5)), (7.4, (1, 1, 1)), (4.2, (0.88, 0.3, 0.4))):
        shape.draw_circle((x, y), r)
        shape.finish(fill=col, color=(0.75, 0.35, 0.42),
                     width=0.8 if col == (1, 1, 1) else 0)
    shape.draw_polyline([(x + 1, y - 1), (x + 13, y - 13)])
    shape.finish(color=(0.3, 0.3, 0.35), width=1.6)
    shape.draw_polyline([(x + 13, y - 13), (x + 9, y - 13), (x + 13, y - 9)])
    shape.finish(color=(0.3, 0.3, 0.35), width=1.6, closePath=False)


def icon_charts(c):
    x, y = c
    for i, (h, col) in enumerate(((7, (0.85, 0.45, 0.3)), (11, (0.9, 0.62, 0.3)),
                                  (5, (0.93, 0.75, 0.4)))):
        shape.draw_rect(pymupdf.Rect(x - 13 + i * 5, y + 9 - h,
                                     x - 9.6 + i * 5, y + 9))
        shape.finish(fill=col, color=col, width=0)
    px, py, r = x + 9, y + 3, 6.4
    shape.draw_circle((px, py), r)
    shape.finish(fill=(0.5, 0.62, 0.82), color=(0.5, 0.62, 0.82), width=0)
    shape.draw_sector((px, py), (px + r, py), 120)
    shape.finish(fill=(0.92, 0.6, 0.32), color=(0.92, 0.6, 0.32), width=0)
    shape.draw_circle((px, py), r)
    shape.finish(fill=None, color=(0.6, 0.5, 0.4), width=0.9)


def icon_monitor(c):
    x, y = c
    shape.draw_rect(pymupdf.Rect(x - 14, y - 11, x + 14, y + 7))
    shape.finish(fill=(1, 1, 1), color=(0.4, 0.5, 0.68), width=1.2)
    px, py, r = x + 5, y - 2, 5
    shape.draw_circle((px, py), r)
    shape.finish(fill=(0.4, 0.55, 0.8), color=(0.4, 0.55, 0.8), width=0)
    shape.draw_sector((px, py), (px + r, py), 130)
    shape.finish(fill=(0.85, 0.45, 0.35), color=(0.85, 0.45, 0.35), width=0)
    shape.draw_circle((px, py), r)
    shape.finish(fill=None, color=(0.4, 0.5, 0.68), width=0.8)
    for i in range(3):
        shape.draw_line((x - 11, y - 6 + i * 3.2), (x - 2 + i * 2, y - 6 + i * 3.2))
        shape.finish(color=(0.6, 0.7, 0.85), width=1.1)
    shape.draw_rect(pymupdf.Rect(x - 5, y + 7, x + 5, y + 10))
    shape.finish(fill=(0.55, 0.62, 0.75), color=(0.55, 0.62, 0.75), width=0)
    shape.draw_rect(pymupdf.Rect(x - 9, y + 10, x + 9, y + 12))
    shape.finish(fill=(0.55, 0.62, 0.75), color=(0.55, 0.62, 0.75), width=0)


# ---------------------------------------------------------- module ----
def module(x0, y0, w, h, bg, border, num, numcol, title, desc, iconf):
    rr(pymupdf.Rect(x0, y0, x0 + w, y0 + h), bg, border)
    iconf((x0 + 26, y0 + h / 2))
    circle((x0 + 55, y0 + 16), 9, numcol, num)
    tx = x0 + 70
    y = y0 + 19
    for s in title:
        ltext((tx, y), s, 10.2, "hebo", TXT)
        y += 12
    y += 3
    for s in desc:
        ltext((tx, y), s, 8.2, "helv", DESC)
        y += 9.2


# ------------------------------------------------------------ build ----
module(100, 18, 290, 50, (0.918, 0.949, 0.988), (0.49, 0.65, 0.85),
       "1", (0.24, 0.42, 0.70), ["Raw FASTQ"],
       ["Input raw sequencing data"], icon_doc)
arrow([(245, 68), (245, 84)])

module(100, 84, 290, 50, (0.914, 0.961, 0.937), (0.49, 0.79, 0.71),
       "2", (0.16, 0.55, 0.45), ["Quality control"],
       ["Adapter and quality trimming"], icon_scissors)
arrow([(245, 134), (245, 150)])

module(100, 150, 290, 56, (0.933, 0.953, 0.878), (0.72, 0.80, 0.50),
       "3", (0.63, 0.68, 0.30), ["Preprocessing"],
       ["Collapse identical sequences", "and sum counts"], icon_table)
arrow([(245, 206), (245, 222)])

module(100, 222, 290, 56, (0.992, 0.965, 0.890), (0.89, 0.77, 0.42),
       "4", (0.91, 0.72, 0.24), ["Mapping to tRNA reference"],
       ["Map reads to tRNA reference", "database"], icon_branch)

elbow([(245, 278), (245, 292), (125, 292), (125, 308)])
head_at((125, 310), math.pi / 2)
elbow([(245, 278), (245, 292), (375, 292), (375, 308)])
head_at((375, 310), math.pi / 2)

module(10, 310, 230, 80, (0.937, 0.918, 0.984), (0.65, 0.60, 0.88),
       "5A", (0.50, 0.37, 0.77), ["Local alignment of", "tsRNA fragments"],
       ["Perform pairwise alignment of", "all tsRNA fragments to identify",
        "read families"], icon_align)
module(260, 310, 230, 80, (0.918, 0.945, 0.988), (0.56, 0.69, 0.86),
       "5B", (0.29, 0.47, 0.76), ["Non-clustering of", "tsRNAs"],
       ["Quantify tsRNAs", "without clustering"],
       lambda c: icon_align(c, teal=True))

arrow([(125, 390), (125, 404)])

module(28, 406, 222, 54, (0.925, 0.906, 0.969), (0.70, 0.64, 0.89),
       "6", (0.55, 0.44, 0.78), ["Clustering of tsRNAs"],
       ["Using network-based", "directional algorithm"], icon_netgraph)
arrow([(139, 460), (139, 474)])

module(28, 476, 222, 66, (0.925, 0.937, 0.980), (0.64, 0.68, 0.85),
       "7", (0.36, 0.44, 0.75), ["Differential expression", "analysis"],
       ["Identify differentially", "expressed tsRNAs"], icon_volcano)
elbow([(375, 390), (375, 509), (254, 509)])
head_at((252, 509), math.pi)
arrow([(139, 542), (139, 556)])

module(28, 558, 222, 50, (0.984, 0.925, 0.937), (0.92, 0.66, 0.73),
       "8", (0.84, 0.43, 0.55), ["tsRNA target prediction"],
       ["Predict targets in CDS", "and/or 3\u2032UTR"], icon_target)
arrow([(139, 608), (139, 622)])

module(28, 624, 222, 52, (0.992, 0.941, 0.867), (0.93, 0.78, 0.58),
       "9", (0.90, 0.60, 0.24), ["Enrichment analysis"],
       ["Using predicted", "gene targets"], icon_charts)
arrow([(139, 676), (139, 690)])

module(28, 692, 222, 50, (0.918, 0.949, 0.988), (0.56, 0.69, 0.86),
       "10", (0.29, 0.47, 0.76), ["Visualization"],
       ["Results visualization", "and interpretation"], icon_monitor)

shape.commit()
# 文字最后画，避免被不透明色块盖住
for pt, s, size, font, color in TEXT_BUFFER:
    page.insert_text(pt, s, fontname=font, fontsize=size, color=color)
pdf_path = os.path.join(OUT, "Figure1_vector.pdf")
doc.save(pdf_path, deflate=True, garbage=4)
doc.close()
print("saved", pdf_path)

# 600 dpi PNG preview
d = pymupdf.open(pdf_path)
d[0].get_pixmap(dpi=150).save(os.path.join(OUT, "Figure1_vector_preview.png"))
print("preview saved")
