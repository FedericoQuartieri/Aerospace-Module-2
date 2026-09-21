#!/usr/bin/env python3
# extracts the hy2Foam curves from the figures of the paper (Casseau et al., Aerospace 2016)
# and saves them in csv files "t,value", one per curve
#
# the figures of the pdf are vector graphics: each page is converted to svg with
# mutool, then the strokes (paths) of the curves are read and mapped to
# physical coordinates using the axis ticks and the numeric labels
#
# usage: python3 extract-paper-curves.py [fig3a fig3b ...]   (without arguments: all)

import os
import re
import sys
import math
import subprocess
import tempfile

PDF = os.path.join(os.path.dirname(__file__), "..", "aerospace-03-00034-1.pdf")
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ---- description of the figures
# page   : page of the pdf
# group  : panel (svg group with clip-path), found by looking at the svg file
# style  : (colour, width, dash pattern) of the hy2Foam paths in that panel
# yscale : the y-axis labels are in K x 10^3 -> multiply by 1000
# curves : name of the curve and expected initial value, used to recognise
#          which stroke is which; if two curves start from the same value
#          the first in the list is the one with the lower mean value
# match  : "final" -> the curves are recognised by their final value, in
#          decreasing order (used for the number densities)
# out    : name of the output files, if different from the name of the figure
FIGURES = {
    "fig3a": dict(page=10, group="clip_5", style=("#000000", "5", None),
                  yscale=1000, curves=[("Ttr", 10000), ("Tv", 1000)]),
    "fig3b": dict(page=10, group="clip_9", style=("#000000", "12", None),
                  yscale=1000, curves=[("Ttr", 3000), ("Tv", 10000)]),
    "fig4-noEl": dict(page=12, group="clip_4", style=("#000000", "12", None),
                      yscale=1000, curves=[("Ttr", 30000), ("Tv", 1000)]),
    "fig4-el": dict(page=12, group="clip_4", style=("#000000", "9", "54,54,18,54"),
                    yscale=1000, curves=[("Ttr", 30000), ("Tv", 1000)]),
    "fig5": dict(page=12, group="clip_10", style=("#000000", "12", None),
                 yscale=1000, curves=[("Ttr", 30000), ("Tv", 1000)]),
    # red curve: hy2Foam with the LeMANS convention (mixture n, sigma 1e-20)
    "fig5-lemans": dict(page=12, group="clip_10", style=("#ff0000", "5", None),
                        yscale=1000, curves=[("Ttr", 30000), ("Tv", 1000)]),
    "fig6-noVV": dict(page=13, group="clip_4", style=("#000000", "20", None),
                      yscale=1000, curves=[("Ttr", 5000), ("Tv_O2", 30000), ("Tv_N2", 30000)]),
    "fig6-VV": dict(page=13, group="clip_4", style=("#000000", "7", None),
                    yscale=1000, curves=[("Ttr", 5000), ("Tv_O2", 30000), ("Tv_N2", 30000)]),
    "fig7a": dict(page=14, group="clip_7", style=("#000000", "12", None),
                  yscale=1000, curves=[("Ttr", 30000), ("Tv", 1000)]),
    "fig7b": dict(page=14, group="clip_11", style=("#000000", "12", None),
                  yscale=1, match="final", curves=[("N", None), ("N2", None)]),
    # in fig 8a Ttr and Tv of hy2Foam-Park have different colours: two entries, same output file
    "fig8a-Ttr": dict(page=15, group="clip_5", style=("#000000", "12", None),
                      yscale=1000, curves=[("Ttr", 30000)], out="fig8a"),
    "fig8a-Tv": dict(page=15, group="clip_5", style=("#0000ff", "12", None),
                     yscale=1000, curves=[("Tv", 30000)], out="fig8a"),
    "fig8b": dict(page=15, group="clip_9", style=("#000000", "12", None),
                  yscale=1, match="final", curves=[("N", None), ("N2", None)]),
    "fig9a": dict(page=16, group="clip_4", style=("#000000", "12", None),
                  yscale=1000, curves=[("T", 10000)]),
    "fig9b": dict(page=16, group="clip_8", style=("#000000", "12", None),
                  yscale=1, ylog=True, match="final",
                  curves=[("N2", None), ("O", None), ("N", None), ("NO", None), ("O2", None)]),
}
# ---- end of description of the figures


def parse_path(d):
    # reads the d attribute of an svg path (only M, L, H, V, Z and implicit lineto)
    # and returns a list of subpaths, each a list of points (x, y)
    toks = re.findall(r"[MLHVZ]|-?[\d.]+(?:e-?\d+)?", d)
    subs = []
    cur = []
    cmd = None
    x = y = 0.0
    i = 0
    while i < len(toks):
        t = toks[i]
        if t in "MLHVZ":
            cmd = t
            i += 1
            if t == "Z" and cur:
                cur.append(cur[0])
                subs.append(cur)
                cur = []
            continue
        if cmd == "M":
            if cur:
                subs.append(cur)
            x, y = float(t), float(toks[i + 1])
            i += 2
            cur = [(x, y)]
            cmd = "L"
        elif cmd == "L":
            x, y = float(t), float(toks[i + 1])
            i += 2
            cur.append((x, y))
        elif cmd == "H":
            x = float(t)
            i += 1
            cur.append((x, y))
        elif cmd == "V":
            y = float(t)
            i += 1
            cur.append((x, y))
        else:
            raise ValueError("unsupported svg command in: " + d[:40])
    if cur:
        subs.append(cur)
    return subs


def page_to_svg(page):
    # converts one page of the pdf to svg (in a temporary folder)
    tmp = tempfile.mkdtemp()
    out = os.path.join(tmp, "page%d.svg")
    subprocess.run(["mutool", "draw", "-o", out, "-F", "svg", PDF, str(page)],
                   check=True, capture_output=True)
    return open(out % page).read()


def read_panel(svg, group):
    # returns the transformation of the panel and its paths with style and points
    body = svg.split("</defs>", 1)[1]
    groups = dict(re.findall(r'<g clip-path="url\(#(clip_\d+)\)">(.*?)</g>', body, re.S))
    g = groups[group]
    tr = re.search(r'transform="matrix\(([^)]*)\)"', g).group(1)
    tr = [float(v) for v in tr.split(",")]
    paths = []
    for p in re.findall(r"<path([^>]*)/>", g):
        d = re.search(r' d="([^"]*)"', p).group(1)
        col = re.search(r'stroke="([^"]*)"', p)
        wid = re.search(r'stroke-width="([^"]*)"', p)
        dash = re.search(r'stroke-dasharray="([^"]*)"', p)
        style = (col.group(1) if col else None,
                 wid.group(1) if wid else None,
                 dash.group(1) if dash else None)
        paths.append((style, parse_path(d)))
    return tr, paths, body


def read_words(body):
    # rebuilds the words from the individual glyphs (data-text) with their position
    uses = re.findall(
        r'<use data-text="([^"]*)" xlink:href="#font_\d+_\d+" '
        r'transform="matrix\(([^)]*)\)"/>', body)
    words = []
    cur = None
    for ch, m in uses:
        a = [float(v) for v in m.split(",")]
        x, y = a[4], a[5]
        size = abs(a[0]) or abs(a[1])
        # rotated text (y-axis title): a[0] = 0
        rotated = a[0] == 0
        # same line and right after the previous letter: same word
        if cur and abs(y - cur["y"]) < 0.5 and abs(x - cur["xend"]) < size * 1.2:
            cur["text"] += ch
            cur["xend"] = x + size * 0.55
        else:
            cur = dict(text=ch, x=x, y=y, xend=x + size * 0.55, size=size, rotated=rotated)
            words.append(cur)
    return words


def find_frame(paths):
    # the frame of the plot is the largest closed rectangle
    # returns (x0, x1, y0, y1) and the style it is drawn with
    best = None
    for style, subs in paths:
        for sp in subs:
            if len(sp) == 5 and sp[0] == sp[-1]:
                xs = [p[0] for p in sp]
                ys = [p[1] for p in sp]
                area = (max(xs) - min(xs)) * (max(ys) - min(ys))
                if best is None or area > best[0]:
                    best = (area, min(xs), max(xs), min(ys), max(ys), style)
    return best[1:5], best[5]


def major_ticks(paths, frame, frame_style):
    # major ticks: the longest segments that touch the frame,
    # drawn with the same stroke as the frame
    x0, x1, y0, y1 = frame
    xt = {}
    yt = {}
    for style, subs in paths:
        if style != frame_style:
            continue
        for sp in subs:
            if len(sp) != 2:
                continue
            (ax, ay), (bx, by) = sp
            if ax == bx and (ay in (y0, y1) or by in (y0, y1)):
                xt[ax] = max(xt.get(ax, 0), abs(by - ay))
            if ay == by and (ax in (x0, x1) or bx in (x0, x1)):
                yt[ay] = max(yt.get(ay, 0), abs(bx - ax))

    def keep_long(t):
        if not t:
            return []
        lmax = max(t.values())
        return sorted(v for v, l in t.items() if l > 0.75 * lmax)
    # the edges of the frame can be major ticks too
    return keep_long(xt) + [x0, x1], keep_long(yt) + [y0, y1]


def axis_labels(words, tr, frame):
    # numeric labels of the axes, already in the local coordinates of the panel
    # returns lists of (value, local coordinate, log?) for x and for y
    sc, e, f = tr[0], tr[4], tr[5]
    x0, x1, y0, y1 = frame
    px0, px1 = e + sc * x0, e + sc * x1
    py_top, py_bot = f - sc * y1, f - sc * y0
    xlab = []
    ylab = []
    skip = set()
    for i, w in enumerate(words):
        t = w["text"].strip()
        if w["rotated"] or i in skip:
            continue
        lab = None
        if t == "10":
            # logarithmic label: "10" followed by the smaller exponent
            nxt = words[i + 1] if i + 1 < len(words) else None
            if nxt and nxt["size"] < w["size"] and 0 < nxt["x"] - w["x"] < 2 * w["size"]:
                try:
                    expo = int(nxt["text"].replace("−", "-"))
                    width = 2 * 0.556 * w["size"] + len(nxt["text"]) * 0.56 * nxt["size"]
                    lab = (expo, True)
                    skip.add(i + 1)
                except ValueError:
                    pass
        if lab is None:
            try:
                val = float(t)
            except ValueError:
                continue
            width = len(t) * 0.556 * w["size"]
            lab = (val, False)
        xc = w["x"] + width / 2
        yc = w["y"] - 0.35 * w["size"]
        # left of the frame and within its height: y axis; below it: x axis
        tol = 0.6 * w["size"]
        if w["xend"] < px0 + 2 and py_top - tol < yc < py_bot + tol:
            ylab.append((lab, (f - yc) / sc))
        elif py_bot < w["y"] < py_bot + 3 * w["size"] and px0 - 20 < xc < px1 + 20:
            xlab.append((lab, (xc - e) / sc))
    return xlab, ylab


def calibrate(labels, ticks, name):
    # matches each label to the nearest major tick and makes a linear fit
    # value = a + b * local coordinate (value = log10 if the axis is logarithmic)
    pts = []
    for (val, is_log), pos in labels:
        near = min(ticks, key=lambda t: abs(t - pos))
        if abs(near - pos) < 90:
            pts.append((near, val))
    pts = sorted(set(pts))
    if len(pts) < 2:
        raise RuntimeError("axis %s: only %d usable labels found" % (name, len(pts)))
    n = len(pts)
    sx = sum(p[0] for p in pts)
    sy = sum(p[1] for p in pts)
    sxx = sum(p[0] * p[0] for p in pts)
    sxy = sum(p[0] * p[1] for p in pts)
    b = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    a = (sy - b * sx) / n
    is_log = labels[0][0][1]
    return a, b, is_log, pts


def chain_subpaths(subs, tol):
    # joins consecutive strokes into continuous curves: the next stroke is the
    # one that starts close to the end of the current stroke and in the same
    # direction (needed for dashed curves, broken into many pieces)
    subs = [sp if sp[0][0] <= sp[-1][0] else sp[::-1] for sp in subs]
    subs.sort(key=lambda sp: sp[0][0])
    used = [False] * len(subs)
    chains = []
    for i in range(len(subs)):
        if used[i]:
            continue
        used[i] = True
        chain = list(subs[i])
        while True:
            ex, ey = chain[-1]
            # local direction of the curve at its end
            px, py = chain[-2] if len(chain) > 1 else (ex - 1, ey)
            dx, dy = ex - px, ey - py
            norm = math.hypot(dx, dy) or 1.0
            dx, dy = dx / norm, dy / norm
            best = None
            for j in range(len(subs)):
                if used[j]:
                    continue
                sx, sy = subs[j][0]
                dist = math.hypot(sx - ex, sy - ey)
                if dist > tol or sx < ex - 1:
                    continue
                # distance from the straight line that extends the curve
                dev = abs((sx - ex) * dy - (sy - ey) * dx)
                score = dev + 0.2 * dist
                if best is None or score < best[0]:
                    best = (score, j)
            if best is None:
                break
            used[best[1]] = True
            chain.extend(subs[best[1]])
        chains.append(chain)
    # discards the short pieces: line samples in the legend
    return [c for c in chains if c[-1][0] - c[0][0] > 600]


def extract(name, cfg):
    svg = page_to_svg(cfg["page"])
    tr, paths, body = read_panel(svg, cfg["group"])
    frame, frame_style = find_frame(paths)
    xticks, yticks = major_ticks(paths, frame, frame_style)
    xlab, ylab = axis_labels(read_words(body), tr, frame)
    ax, bx, xlog, xpts = calibrate(xlab, xticks, "x")
    ay, by, ylog, ypts = calibrate(ylab, yticks, "y")
    print("%s: x axis %s, y axis %s" % (name, xpts, ypts))

    # strokes of the hy2Foam curves (same style) joined into continuous curves
    subs = [sp for style, ss in paths if style == cfg["style"] for sp in ss if len(sp) >= 2]
    tol = 2.0 if cfg["style"][2] is None else 250.0
    chains = chain_subpaths(subs, tol)

    def to_data(pt):
        x = ax + bx * pt[0]
        y = ay + by * pt[1]
        x = 10 ** x if xlog else x
        y = 10 ** y if ylog else y
        return x, y * cfg["yscale"]

    curves = []
    for ch in chains:
        pts = [to_data(p) for p in ch if frame[0] - 1 <= p[0] <= frame[1] + 1]
        curves.append(pts)

    # ---- recognises the curves
    names = [c[0] for c in cfg["curves"]]
    if cfg.get("match") == "final":
        # sorted by decreasing final value
        curves.sort(key=lambda c: -c[-1][1])
        assigned = list(zip(names, curves))
    else:
        assigned = []
        for cname, y0 in cfg["curves"]:
            if not curves:
                raise RuntimeError("%s: no curve left for %s" % (name, cname))
            # candidates: curves that start from the expected initial value,
            # otherwise the one that starts closest to it
            cand = [c for c in curves if abs(c[0][1] - y0) < 0.08 * y0]
            if not cand:
                cand = [min(curves, key=lambda c: abs(c[0][1] - y0))]
            # for the same starting value, first the one with the lower mean value
            cand.sort(key=lambda c: sum(p[1] for p in c) / len(c))
            curves.remove(cand[0])
            assigned.append((cname, cand[0]))
    # ---- end of recognition

    for cname, pts in assigned:
        fn = os.path.join(OUT_DIR, "%s-%s.csv" % (cfg.get("out", name), cname))
        with open(fn, "w") as f:
            f.write("t,%s\n" % cname)
            for x, y in pts:
                f.write("%.6g,%.6g\n" % (x, y))
        print("   %-6s %4d points, from (%.3g, %.5g) to (%.3g, %.5g) -> %s"
              % (cname, len(pts), pts[0][0], pts[0][1], pts[-1][0], pts[-1][1],
                 os.path.basename(fn)))
    if curves and cfg.get("match") != "final":
        print("   warning: %d curves with the same style not assigned, starting from %s"
              % (len(curves), [round(c[0][1]) for c in curves]))


if __name__ == "__main__":
    wanted = sys.argv[1:] or list(FIGURES)
    for name in wanted:
        extract(name, FIGURES[name])
