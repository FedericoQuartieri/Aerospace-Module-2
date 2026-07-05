"""Isothermal-wall unit test: the decoded wall temperature must equal the
imposed 300 K (datum-consistent energy BC) and the near-wall profile must
cool monotonically towards the wall."""
import os, re, sys
import numpy as np

NC = 50
def read(tdir, name):
    txt = open(os.path.join(tdir, name)).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\((.*?)\)\s*;", txt, re.S)
    vals = np.array([float(q) for q in m.group(1).split()]) if m else \
        np.full(NC, float(re.search(r"internalField\s+uniform\s+([^;]+);", txt).group(1)))
    mb = re.search(r"wall\s*\{[^}]*?value\s+uniform\s+([0-9.eE+-]+)", txt, re.S)
    wall = float(mb.group(1)) if mb else None
    if wall is None:
        mb = re.search(r"wall\s*\{[^}]*?value\s+nonuniform\s+List<scalar>\s*\d+\s*\(([^)]*)\)", txt, re.S)
        wall = float(mb.group(1).split()[0]) if mb else float('nan')
    return vals, wall

times = sorted((d for d in os.listdir(".") if re.fullmatch(r"[0-9.e+-]+", d) and float(d) > 0), key=float)
t = times[-1]
T, Tw = read(t, "T")

print(f"t = {t} s")
print(f"wall T (decoded)      = {Tw:.4f} K   (imposed: 300)")
print(f"first 5 cells T       = {np.array2string(T[:5], precision=1)}")
mono = np.all(np.diff(T[:10]) >= -1e-6)
print(f"monotone cooling near wall: {mono}")
ok = abs(Tw - 300) < 0.5 and mono and T[0] < 999
print("WALL TEST:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
