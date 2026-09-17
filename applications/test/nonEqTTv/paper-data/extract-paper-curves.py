#!/usr/bin/env python3
# estrae le curve di hy2Foam dalle figure del paper (Casseau et al., Aerospace 2016)
# e le salva in file csv "t,valore", uno per curva
#
# le figure del pdf sono vettoriali: ogni pagina viene convertita in svg con
# mutool, poi si leggono i tratti (path) delle curve e li si riporta alle
# coordinate fisiche usando i tick degli assi e le etichette numeriche
#
# uso: python3 extract-paper-curves.py [fig3a fig3b ...]   (senza argomenti: tutte)

import os
import re
import sys
import math
import subprocess
import tempfile

PDF = os.path.join(os.path.dirname(__file__), "..", "aerospace-03-00034-1.pdf")
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ---- descrizione delle figure
# page   : pagina del pdf
# group  : pannello (gruppo svg con clip-path), trovato guardando il file svg
# style  : (colore, spessore, tratteggio) dei path di hy2Foam in quel pannello
# yscale : le etichette dell'asse y sono in K x 10^3 -> moltiplico per 1000
# curves : nome della curva e valore iniziale atteso, serve per riconoscere
#          quale tratto e' quale; se due curve partono dallo stesso valore
#          la prima della lista e' quella con valore medio piu' basso
# match  : "final" -> le curve si riconoscono dal valore finale, in ordine
#          decrescente (usato per le densita' numeriche)
# out    : nome dei file di uscita, se diverso dal nome della figura
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
    # curva rossa: hy2Foam con la convenzione di LeMANS (n miscela, sigma 1e-20)
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
    # in fig 8a Ttr e Tv di hy2Foam-Park hanno colori diversi: due voci, stesso file di uscita
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
# ---- fine descrizione delle figure


def parse_path(d):
    # legge l'attributo d di un path svg (solo M, L, H, V, Z e lineto impliciti)
    # e ritorna una lista di sottopercorsi, ognuno lista di punti (x, y)
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
            raise ValueError("comando svg non gestito in: " + d[:40])
    if cur:
        subs.append(cur)
    return subs


def page_to_svg(page):
    # converte una pagina del pdf in svg (in una cartella temporanea)
    tmp = tempfile.mkdtemp()
    out = os.path.join(tmp, "page%d.svg")
    subprocess.run(["mutool", "draw", "-o", out, "-F", "svg", PDF, str(page)],
                   check=True, capture_output=True)
    return open(out % page).read()


def read_panel(svg, group):
    # ritorna la trasformazione del pannello e i suoi path con stile e punti
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
    # ricompone le parole dai singoli glifi (data-text) con la loro posizione
    uses = re.findall(
        r'<use data-text="([^"]*)" xlink:href="#font_\d+_\d+" '
        r'transform="matrix\(([^)]*)\)"/>', body)
    words = []
    cur = None
    for ch, m in uses:
        a = [float(v) for v in m.split(",")]
        x, y = a[4], a[5]
        size = abs(a[0]) or abs(a[1])
        # testo ruotato (titolo dell'asse y): a[0] = 0
        rotated = a[0] == 0
        # stessa riga e subito dopo la lettera precedente: stessa parola
        if cur and abs(y - cur["y"]) < 0.5 and abs(x - cur["xend"]) < size * 1.2:
            cur["text"] += ch
            cur["xend"] = x + size * 0.55
        else:
            cur = dict(text=ch, x=x, y=y, xend=x + size * 0.55, size=size, rotated=rotated)
            words.append(cur)
    return words


def find_frame(paths):
    # la cornice del grafico e' il rettangolo chiuso piu' grande
    # ritorna (x0, x1, y0, y1) e lo stile con cui e' disegnata
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
    # tick maggiori: i segmenti piu' lunghi che toccano la cornice,
    # disegnati con lo stesso tratto della cornice
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
    # anche i bordi della cornice possono essere tick maggiori
    return keep_long(xt) + [x0, x1], keep_long(yt) + [y0, y1]


def axis_labels(words, tr, frame):
    # etichette numeriche degli assi, gia' in coordinate locali del pannello
    # ritorna liste di (valore, coordinata locale, log?) per x e per y
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
            # etichetta logaritmica: "10" seguito dall'esponente piu' piccolo
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
        # a sinistra della cornice e alla sua altezza: asse y; sotto: asse x
        tol = 0.6 * w["size"]
        if w["xend"] < px0 + 2 and py_top - tol < yc < py_bot + tol:
            ylab.append((lab, (f - yc) / sc))
        elif py_bot < w["y"] < py_bot + 3 * w["size"] and px0 - 20 < xc < px1 + 20:
            xlab.append((lab, (xc - e) / sc))
    return xlab, ylab


def calibrate(labels, ticks, name):
    # associa ogni etichetta al tick maggiore piu' vicino e fa un fit lineare
    # valore = a + b * coordinata locale (valore = log10 se asse logaritmico)
    pts = []
    for (val, is_log), pos in labels:
        near = min(ticks, key=lambda t: abs(t - pos))
        if abs(near - pos) < 90:
            pts.append((near, val))
    pts = sorted(set(pts))
    if len(pts) < 2:
        raise RuntimeError("asse %s: trovate solo %d etichette utili" % (name, len(pts)))
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
    # unisce i tratti consecutivi in curve continue: il tratto seguente e'
    # quello che inizia vicino alla fine del tratto corrente e nella stessa
    # direzione (serve per le curve tratteggiate, spezzate in tanti pezzi)
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
            # direzione locale della curva alla sua fine
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
                # distanza dalla retta che prolunga la curva
                dev = abs((sx - ex) * dy - (sy - ey) * dx)
                score = dev + 0.2 * dist
                if best is None or score < best[0]:
                    best = (score, j)
            if best is None:
                break
            used[best[1]] = True
            chain.extend(subs[best[1]])
        chains.append(chain)
    # scarta i pezzi corti: campioni di linea nella legenda
    return [c for c in chains if c[-1][0] - c[0][0] > 600]


def extract(name, cfg):
    svg = page_to_svg(cfg["page"])
    tr, paths, body = read_panel(svg, cfg["group"])
    frame, frame_style = find_frame(paths)
    xticks, yticks = major_ticks(paths, frame, frame_style)
    xlab, ylab = axis_labels(read_words(body), tr, frame)
    ax, bx, xlog, xpts = calibrate(xlab, xticks, "x")
    ay, by, ylog, ypts = calibrate(ylab, yticks, "y")
    print("%s: asse x %s, asse y %s" % (name, xpts, ypts))

    # tratti delle curve di hy2Foam (stesso stile) uniti in curve continue
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

    # ---- riconosce le curve
    names = [c[0] for c in cfg["curves"]]
    if cfg.get("match") == "final":
        # ordinate per valore finale decrescente
        curves.sort(key=lambda c: -c[-1][1])
        assigned = list(zip(names, curves))
    else:
        assigned = []
        for cname, y0 in cfg["curves"]:
            if not curves:
                raise RuntimeError("%s: nessuna curva rimasta per %s" % (name, cname))
            # candidate: curve che partono dal valore iniziale atteso,
            # altrimenti quella che parte piu' vicino
            cand = [c for c in curves if abs(c[0][1] - y0) < 0.08 * y0]
            if not cand:
                cand = [min(curves, key=lambda c: abs(c[0][1] - y0))]
            # a parita' di partenza, prima quella con valore medio piu' basso
            cand.sort(key=lambda c: sum(p[1] for p in c) / len(c))
            curves.remove(cand[0])
            assigned.append((cname, cand[0]))
    # ---- fine riconoscimento

    for cname, pts in assigned:
        fn = os.path.join(OUT_DIR, "%s-%s.csv" % (cfg.get("out", name), cname))
        with open(fn, "w") as f:
            f.write("t,%s\n" % cname)
            for x, y in pts:
                f.write("%.6g,%.6g\n" % (x, y))
        print("   %-6s %4d punti, da (%.3g, %.5g) a (%.3g, %.5g) -> %s"
              % (cname, len(pts), pts[0][0], pts[0][1], pts[-1][0], pts[-1][1],
                 os.path.basename(fn)))
    if curves and cfg.get("match") != "final":
        print("   attenzione: %d curve con lo stesso stile non assegnate, partono da %s"
              % (len(curves), [round(c[0][1]) for c in curves]))


if __name__ == "__main__":
    wanted = sys.argv[1:] or list(FIGURES)
    for name in wanted:
        extract(name, FIGURES[name])
