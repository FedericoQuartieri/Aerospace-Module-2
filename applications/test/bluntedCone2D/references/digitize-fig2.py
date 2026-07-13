"""Vector-extraction digitisation of Figure 2 (Casseau et al., Aerospace
2016, 3, 45 - Part Two, CC-BY 4.0), page 9 of the published PDF.

Unlike Engauge-style manual picking, this reads the curve coordinates
directly from the PDF vector drawing commands (PyMuPDF get_drawings), so
the accuracy is limited only by the axis calibration, not by hand-eye
work.

Panels (3 rows x 2 cols):
  (a) T/Tinf vs stagnation-line position [mm]   x [-3, 0],  y [0, 35]
  (b) rho/rhoinf vs position [mm]               x [-3, 0],  y [0, 15]
  (c) U/Uinf vs position [mm]                   x [-3, 0],  y [0, 1]
  (d) pressure coefficient vs axial dist [cm]   x [0, 4],   y [0, 1]
  (e) friction coefficient                      x [0, 4],   y [0, 0.15]
  (f) Stanton number                            x [0, 4],   y [0, 0.2]

Series classification inside each panel:
  - blue stroke               -> hy2Foam, 10 um first spacing (e, f only)
  - black long polylines      -> hy2Foam (solid)
  - black short dash segments -> CFD Michigan (dash-dot)
  - filled tiny paths         -> DSMC MONACO triangles (use centroid)
  - remaining short verticals with caps -> experiment error bars (d, f):
    the experiment marker centre is taken as the midpoint of each bar

Legend boxes (inner rectangles) are masked out.

Outputs: fig2<panel>-<series>.csv (x,y in data units) + a QA overlay
image fig2-overlay.png rendered from the page with extracted points on
top - inspect it before trusting the CSVs.

Usage: python3 digitize-fig2.py
"""

import os
import fitz
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "aerospace-03-00045.pdf")
PAGE = 8  # 0-based: page 9 of the PDF

# panel name -> (data x range, data y range)
PANELS = {
    "a": ((-3.0, 0.0), (0.0, 35.0)),
    "b": ((-3.0, 0.0), (0.0, 15.0)),
    "c": ((-3.0, 0.0), (0.0, 1.0)),
    "d": ((0.0, 4.0), (0.0, 1.0)),
    "e": ((0.0, 4.0), (0.0, 0.15)),
    "f": ((0.0, 4.0), (0.0, 0.2)),
}

doc = fitz.open(PDF)
page = doc[PAGE]
drawings = page.get_drawings()

# ---------------------------------------------------------------- panels ---
# the six axis boxes are the six largest stroked rectangles
rects = []
for d in drawings:
    for item in d["items"]:
        if item[0] == "re":
            r = item[1]
            if r.width > 100 and r.height > 80:
                rects.append(fitz.Rect(r))
# dedupe (stroke+fill can duplicate)
uniq = []
for r in rects:
    if not any(abs(r.x0 - u.x0) < 2 and abs(r.y0 - u.y0) < 2 for u in uniq):
        uniq.append(r)
uniq.sort(key=lambda r: (round(r.y0), round(r.x0)))
assert len(uniq) == 6, f"expected 6 panel boxes, got {len(uniq)}"
panelRect = dict(zip("abcdef", uniq))

# legend boxes: smaller stroked rectangles fully inside a panel
legendRects = []
for d in drawings:
    for item in d["items"]:
        if item[0] == "re":
            r = fitz.Rect(item[1])
            for p in uniq:
                if p.contains(r) and r.width < 0.9*p.width and r.height < 0.7*p.height:
                    legendRects.append(r)

def in_legend(pt):
    return any(r.contains(fitz.Point(pt)) for r in legendRects)

# ------------------------------------------------------------ transforms ---
def make_maps(panel):
    pr = panelRect[panel]
    (x0, x1), (y0, y1) = PANELS[panel]

    def to_data(px, py):
        xd = x0 + (px - pr.x0)/(pr.x1 - pr.x0)*(x1 - x0)
        yd = y0 + (pr.y1 - py)/(pr.y1 - pr.y0)*(y1 - y0)  # pdf y grows down
        return xd, yd

    return to_data

# ------------------------------------------------------------- classify ---
series = {p: {"hy2foam": [], "hy2foam10um": [], "michigan": [],
              "dsmc": [], "experiments": []} for p in PANELS}

def which_panel(rect):
    c = fitz.Point((rect.x0 + rect.x1)/2, (rect.y0 + rect.y1)/2)
    for name, pr in panelRect.items():
        if pr.contains(c):
            return name
    return None

def is_tick_segment(p1, p2, pr, B=10.0, maxlen=13.0):
    """A tick: short segment with both endpoints hugging the same border.
    gnuplot emits mirror ticks (bottom+top / left+right) in a single path,
    so filtering must happen per segment, not per path."""
    if p1.distance_to(p2) > maxlen:
        return False
    for lo, hi, coord in (
        (pr.x0 - B, pr.x0 + B, 0), (pr.x1 - B, pr.x1 + B, 0),
        (pr.y0 - B, pr.y0 + B, 1), (pr.y1 - B, pr.y1 + B, 1),
    ):
        a = (p1.x, p1.y)[coord]
        b = (p2.x, p2.y)[coord]
        if lo <= a <= hi and lo <= b <= hi:
            return True
    return False


for d in drawings:
    pn = which_panel(d["rect"])
    if pn is None:
        continue
    color = d.get("color")
    fill = d.get("fill")
    pr = panelRect[pn]

    pts = []
    segs = []          # (p1, p2) after tick filtering, in path order
    seglen = 0.0
    for item in d["items"]:
        if item[0] == "l":
            p1, p2 = fitz.Point(item[1]), fitz.Point(item[2])
        elif item[0] == "c":
            p1, p2 = fitz.Point(item[1]), fitz.Point(item[4])
        else:
            continue
        if is_tick_segment(p1, p2, pr):
            continue
        pts.append((p1.x, p1.y))
        pts.append((p2.x, p2.y))
        segs.append((p1, p2))
        seglen += p1.distance_to(p2)
    if not pts:
        continue
    if any(in_legend(p) for p in pts):
        continue
    arr = np.array(pts)
    w = arr[:, 0].max() - arr[:, 0].min()
    h = arr[:, 1].max() - arr[:, 1].min()

    # contiguity: fraction of consecutive segments that join end-to-start.
    # A solid curve is one connected polyline (high contiguity); a
    # dash-dot path carries many detached little strokes (low contiguity).
    if len(segs) > 1:
        joined = sum(
            1 for i in range(len(segs) - 1)
            if segs[i][1].distance_to(segs[i + 1][0]) < 0.5
        )
        contiguity = joined/(len(segs) - 1)
    else:
        contiguity = 1.0

    if color == (0.0, 0.0, 1.0):
        series[pn]["hy2foam10um"].append(arr)
    elif fill is not None and color is None:
        # filled symbol (DSMC triangle): keep centroid
        if w < 8 and h < 8:
            series[pn]["dsmc"].append(arr.mean(axis=0, keepdims=True))
    elif color == (0.0, 0.0, 0.0):
        if contiguity > 0.8 and seglen > 25:
            # connected long polyline: solid hy2Foam curve
            series[pn]["hy2foam"].append(arr)
        elif len(segs) <= 3 and w < 10 and h < 18 and seglen < 25:
            # error bar: vertical stem +/- horizontal caps, tiny bbox.
            # The experiment point is the centre of the bar.
            series[pn]["experiments"].append(
                arr.mean(axis=0, keepdims=True))
        else:
            # detached strokes: Michigan dash-dot (keep stroke midpoints
            # to avoid double-weighting dash endpoints)
            mids = np.array([[(a.x + b.x)/2, (a.y + b.y)/2]
                             for a, b in segs])
            series[pn]["michigan"].append(mids)

# ------------------------------------------------------------- outputs ----
def save_csv(panel, name, pts_pdf):
    to_data = make_maps(panel)
    data = np.array([to_data(px, py) for px, py in pts_pdf])
    order = np.argsort(data[:, 0])
    data = data[order]
    path = os.path.join(HERE, f"fig2{panel}-{name}.csv")
    np.savetxt(path, data, delimiter=",", header="x,y", comments="")
    return data

summary = []
extracted = {}
for pn in PANELS:
    for name, chunks in series[pn].items():
        if not chunks:
            continue
        allpts = np.vstack(chunks)
        # thin out duplicated vertices in polylines
        if name in ("hy2foam", "hy2foam10um", "michigan"):
            allpts = np.unique(np.round(allpts, 2), axis=0)
        data = save_csv(pn, name, allpts)
        extracted[(pn, name)] = allpts
        summary.append(f"fig2{pn}-{name}.csv: {len(data)} punti, "
                       f"x [{data[:,0].min():.3g}, {data[:,0].max():.3g}], "
                       f"y [{data[:,1].min():.3g}, {data[:,1].max():.3g}]")

print("\n".join(summary))

# ------------------------------------------------------------ QA overlay --
zoom = 3
pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom))
img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)

fig, ax = plt.subplots(figsize=(14, 20))
ax.imshow(img)
cols = {"hy2foam": "red", "hy2foam10um": "cyan", "michigan": "lime",
        "dsmc": "orange", "experiments": "magenta"}
for (pn, name), pts in extracted.items():
    ax.scatter(pts[:, 0]*zoom, pts[:, 1]*zoom, s=1.5,
               c=cols[name], label=f"{pn}-{name}")
ax.set_title("QA overlay: punti estratti sopra il render della pagina")
ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig2-overlay.png"), dpi=110)
print("QA overlay salvato in fig2-overlay.png - CONTROLLARLO")
