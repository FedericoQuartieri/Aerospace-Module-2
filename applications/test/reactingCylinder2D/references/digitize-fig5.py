"""Vector-extraction digitisation of Figure 5 (Casseau et al., Aerospace
2016, 3, 45 - Part Two, CC-BY 4.0), page 12 of the published PDF: the
Mach 20 reacting cylinder references.

Same technique as bluntedCone2D/references/digitize-fig2.py (see the
provenance and caveat notes there): PyMuPDF get_drawings, panel boxes
from the six largest stroked rectangles, per-segment tick filtering,
legend masking. The PDF is NOT duplicated: this script reads the copy
vendored with the cone case.

Panels (3 rows x 2 cols):
  (a) Mach vs stagnation-line position [m]      x [-1.5, -1], y [0, 20]
  (b) T vs position [K x 10^3]                  x [-1.5, -1], y [0, 15]
  (c) number density [1/m^3], LOG y             x [-1.5, -1], y [1e18, 1e22]
  (d) pressure coefficient vs theta [deg]       x [0, 180],   y [0, 2]
  (e) skin-friction coefficient vs theta        x [0, 180],   y [0, 0.06]
  (f) surface heat flux [W/cm^2] vs theta       x [0, 180],   y [0, 15]

Series (Fig. 5 caption: run 1 = black, run 2 = red, run 3 = blue; the
symbols are dsmcFoam (+ -> short 2-segment crosses) and QK circles):
  - blue strokes  -> run3   (Park TTv + Park rates = OUR setup)
  - red strokes   -> run2   (SSH, CVDV-QK)
  - black strokes -> run1 / hy2Foam non-reacting in (d-f)
  - tiny black crosses  -> dsmc (marker centre)
  - small circles/curves -> qk  (marker centre)
Gray Kn_GLL band fills are ignored (fill colour not black/white).

Outputs: fig5<panel>-<series>.csv + QA overlay fig5-overlay.png -
inspect the overlay after every regeneration.

Usage: python3 digitize-fig5.py
"""

import os
import math
import fitz
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(
    HERE, "..", "..", "bluntedCone2D", "references", "aerospace-03-00045.pdf")
PAGE = 12  # 0-based

# panel name -> (x range, y range, logy)
PANELS = {
    "a": ((-1.5, -1.0), (0.0, 20.0), False),
    "b": ((-1.5, -1.0), (0.0, 15.0), False),
    "c": ((-1.5, -1.0), (18.0, 22.0), True),   # log10(n)
    "d": ((0.0, 180.0), (0.0, 2.0), False),
    "e": ((0.0, 180.0), (0.0, 0.06), False),
    "f": ((0.0, 180.0), (0.0, 15.0), False),
}

doc = fitz.open(PDF)
page = doc[PAGE]
drawings = page.get_drawings()

# ---------------------------------------------------------------- panels ---
rects = []
for d in drawings:
    for item in d["items"]:
        if item[0] == "re":
            r = item[1]
            if r.width > 100 and r.height > 80:
                rects.append(fitz.Rect(r))
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
                if p.contains(r) and r.width < 0.9*p.width \
                        and r.height < 0.7*p.height:
                    legendRects.append(r)

def in_legend(pt):
    return any(r.contains(fitz.Point(pt)) for r in legendRects)

# ------------------------------------------------------------ transforms ---
def make_maps(panel):
    pr = panelRect[panel]
    (x0, x1), (y0, y1), logy = PANELS[panel]

    def to_data(px, py):
        xd = x0 + (px - pr.x0)/(pr.x1 - pr.x0)*(x1 - x0)
        yd = y0 + (pr.y1 - py)/(pr.y1 - pr.y0)*(y1 - y0)
        if logy:
            yd = 10.0**yd
        return xd, yd

    return to_data

# ------------------------------------------------------------- classify ---
series = {p: {"run1": [], "run2": [], "run3": [],
              "dsmc": [], "qk": []} for p in PANELS}

def which_panel(rect):
    c = fitz.Point((rect.x0 + rect.x1)/2, (rect.y0 + rect.y1)/2)
    for name, pr in panelRect.items():
        if pr.contains(c):
            return name
    return None

def is_tick_segment(p1, p2, pr, B=10.0, maxlen=13.0):
    """Tick = short segment in a border band AND roughly PERPENDICULAR
    to that border. The perpendicularity test is what saves the many
    legitimate curve runs that hug the x axis in this figure (post-shock
    Mach ~ 0, Cp/Cf/q -> 0 at high theta): those are parallel to the
    border and must survive (the Fig 2 version of this filter ate
    them)."""
    if p1.distance_to(p2) > maxlen:
        return False
    dx = abs(p2.x - p1.x)
    dy = abs(p2.y - p1.y)
    for lo, hi, coord, perp in (
        (pr.x0 - B, pr.x0 + B, 0, "h"), (pr.x1 - B, pr.x1 + B, 0, "h"),
        (pr.y0 - B, pr.y0 + B, 1, "v"), (pr.y1 - B, pr.y1 + B, 1, "v"),
    ):
        a = (p1.x, p1.y)[coord]
        b = (p2.x, p2.y)[coord]
        if lo <= a <= hi and lo <= b <= hi:
            # vertical border (x band) -> tick is horizontal; horizontal
            # border (y band) -> tick is vertical
            if perp == "h" and dx > 2*dy:
                return True
            if perp == "v" and dy > 2*dx:
                return True
    return False

def colour_of(c):
    if c is None:
        return None
    r, g, b = c
    if b > 0.5 and r < 0.5 and g < 0.5:
        return "blue"
    if r > 0.5 and g < 0.5 and b < 0.5:
        return "red"
    if r < 0.35 and g < 0.35 and b < 0.35:
        return "black"
    return "other"

for d in drawings:
    pn = which_panel(d["rect"])
    if pn is None:
        continue
    col = colour_of(d.get("color"))
    fill = d.get("fill")
    pr = panelRect[pn]

    pts = []
    segs = []
    seglen = 0.0
    hasCurve = False
    for item in d["items"]:
        if item[0] == "l":
            p1, p2 = fitz.Point(item[1]), fitz.Point(item[2])
        elif item[0] == "c":
            p1, p2 = fitz.Point(item[1]), fitz.Point(item[4])
            hasCurve = True
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

    if col == "blue":
        series[pn]["run3"].append(arr)
    elif col == "red":
        series[pn]["run2"].append(arr)
    elif col == "black":
        # markers: crosses (2 straight segments, tiny bbox) or circles
        # (bezier curves, tiny bbox); everything else is the run-1 /
        # non-reacting reference line work
        if hasCurve and w < 8 and h < 8:
            series[pn]["qk"].append(arr.mean(axis=0, keepdims=True))
        elif len(segs) <= 3 and w < 8 and h < 8:
            series[pn]["dsmc"].append(arr.mean(axis=0, keepdims=True))
        else:
            series[pn]["run1"].append(arr)

# ------------------------------------------------------------- outputs ----
def save_csv(panel, name, pts_pdf):
    to_data = make_maps(panel)
    data = np.array([to_data(px, py) for px, py in pts_pdf])
    order = np.argsort(data[:, 0])
    data = data[order]
    path = os.path.join(HERE, f"fig5{panel}-{name}.csv")
    np.savetxt(path, data, delimiter=",", header="x,y", comments="")
    return data

summary = []
extracted = {}
for pn in PANELS:
    for name, chunks in series[pn].items():
        if not chunks:
            continue
        allpts = np.vstack(chunks)
        if name in ("run1", "run2", "run3"):
            allpts = np.unique(np.round(allpts, 2), axis=0)
        data = save_csv(pn, name, allpts)
        extracted[(pn, name)] = allpts
        summary.append(f"fig5{pn}-{name}.csv: {len(data)} punti, "
                       f"x [{data[:,0].min():.3g}, {data[:,0].max():.3g}], "
                       f"y [{data[:,1].min():.3g}, {data[:,1].max():.3g}]")

print("\n".join(summary))

# ------------------------------------------------------------ QA overlay --
zoom = 3
pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom))
img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(
    pix.height, pix.width, pix.n)

fig, ax = plt.subplots(figsize=(14, 20))
ax.imshow(img)
cols = {"run1": "magenta", "run2": "orange", "run3": "cyan",
        "dsmc": "lime", "qk": "red"}
for (pn, name), pts in extracted.items():
    ax.scatter(pts[:, 0]*zoom, pts[:, 1]*zoom, s=1.5,
               c=cols[name], label=f"{pn}-{name}")
ax.set_title("QA overlay: punti estratti sopra il render della pagina")
ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig5-overlay.png"), dpi=110)
print("QA overlay salvato in fig5-overlay.png - CONTROLLARLO")
