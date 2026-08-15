#!/usr/bin/env python3
"""
analyze.py -- Auswertung der seriellen Kaskaden-Ausgaben.
Erzeugt aus den drei CSVs eines Laufs (Prefix) drei Grafiken:
  1) Energiebilanz ueber Zeit  (Korrektheits-Check, sichtbar)
  2) Schadenskarte             (permanente Defekte = gerissene Bindungen)
  3) Defektcluster-Groessenverteilung (log-log; Vorstufe zum Power-Law-Fit)

Aufruf:  python3 analyze.py [prefix]      (default: run)
Liest:   <prefix>_energy.csv, <prefix>_state_final.csv, <prefix>_broken_bonds.csv
"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

prefix = sys.argv[1] if len(sys.argv) > 1 else "run"
L0 = 1.0
DISP_THRESH = 0.5 * L0     # Atom gilt als "verlagert" ab dieser Verschiebung

# ---------------------------------------------------------------- Energie
e = np.genfromtxt(f"{prefix}_energy.csv", delimiter=",", names=True)

# ---------------------------------------------------------------- Zustand
d = np.genfromtxt(f"{prefix}_state_final.csv", delimiter=",", names=True)
x, y, disp = d["x"], d["y"], d["disp"]
moved = disp > DISP_THRESH                  # permanent verlagerte Atome

# ------------------------------------------- Defektcluster (Union-Find)
# Verlagerte Atome, die raeumlich benachbart sind (< 1.5*L0), bilden ein Cluster.
def cluster_sizes(px, py, link=1.5*L0):
    n = len(px)
    if n == 0: return np.array([])
    parent = list(range(n))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    # Nachbarsuche ueber ein grobes Gitter (Bucketing), O(n) statt O(n^2)
    cell = link
    buckets = {}
    for i in range(n):
        key = (int(px[i]//cell), int(py[i]//cell))
        buckets.setdefault(key, []).append(i)
    link2 = link*link
    for i in range(n):
        cx, cy = int(px[i]//cell), int(py[i]//cell)
        for ax in (cx-1,cx,cx+1):
            for ay in (cy-1,cy,cy+1):
                for j in buckets.get((ax,ay), ()):
                    if j <= i: continue
                    if (px[i]-px[j])**2 + (py[i]-py[j])**2 < link2:
                        union(i, j)
    roots = {}
    for i in range(n):
        r = find(i); roots[r] = roots.get(r, 0) + 1
    return np.array(sorted(roots.values(), reverse=True))

sizes = cluster_sizes(x[moved], y[moved])

# ================================================================ Plot
fig = plt.figure(figsize=(15, 4.6))

# (1) Energie
ax1 = fig.add_subplot(1, 3, 1)
ax1.plot(e["t"], e["E_kin"],    label="E_kin")
ax1.plot(e["t"], e["E_spring"], label="E_Feder")
ax1.plot(e["t"], e["E_rep"],    label="E_Abstoss")
ax1.plot(e["t"], e["E_broken"], label="E_Bruch")
ax1.plot(e["t"], e["E_total"],  "k--", lw=2, label="E_total")
ax1.set_xlabel("Zeit"); ax1.set_ylabel("Energie")
ax1.set_title("1) Energiebilanz (E_total flach = korrekt)")
ax1.legend(fontsize=8); ax1.grid(alpha=.3)

# (2) Schadenskarte aus gerissenen Bindungen, auf den Kern gezoomt
ax2 = fig.add_subplot(1, 3, 2)
try:
    bb = np.genfromtxt(f"{prefix}_broken_bonds.csv", delimiter=",", names=True)
    xa, ya = np.atleast_1d(bb["xa"]), np.atleast_1d(bb["ya"])
    xb, yb = np.atleast_1d(bb["xb"]), np.atleast_1d(bb["yb"])
    nseg = len(xa)
    if nseg > 0:
        segs = [[[xa[i],ya[i]],[xb[i],yb[i]]] for i in range(nseg)]
        ax2.add_collection(LineCollection(segs, colors="crimson", linewidths=1.0))
        cx, cy = np.r_[xa,xb].mean(), np.r_[ya,yb].mean()
        R = max(15.0, 1.3*max(np.ptp(np.r_[xa,xb]), np.ptp(np.r_[ya,yb]))/2)
        ax2.set_xlim(cx-R, cx+R); ax2.set_ylim(cy-R, cy+R)
    ax2.set_aspect("equal")
    ax2.set_title(f"2) Schaden: {nseg} gerissene Bindungen")
    ax2.set_xlabel("x"); ax2.set_ylabel("y")
except (OSError, ValueError):
    ax2.set_title("2) keine broken_bonds.csv")

# (3) Clustergroessen log-log
ax3 = fig.add_subplot(1, 3, 3)
if len(sizes) > 0:
    vals, counts = np.unique(sizes, return_counts=True)
    ax3.loglog(vals, counts, "o", ms=6)
    ax3.set_xlabel("Clustergroesse n"); ax3.set_ylabel("Haeufigkeit")
    ax3.set_title(f"3) {len(sizes)} Cluster (groesster: {sizes.max()})")
    ax3.grid(alpha=.3, which="both")
else:
    ax3.set_title("3) keine Cluster ueber Schwelle")

fig.suptitle(f"Lauf: {prefix}   |   verlagerte Atome: {int(moved.sum())}", fontsize=11)
fig.tight_layout()
out = f"{prefix}_analysis.png"
fig.savefig(out, dpi=130)
print(f"geschrieben: {out}")
print(f"verlagerte Atome (disp>{DISP_THRESH}): {int(moved.sum())}")
print(f"Cluster: {len(sizes)}  Groessen (Top10): {sizes[:10].tolist()}")
