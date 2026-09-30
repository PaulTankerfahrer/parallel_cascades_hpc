#!/usr/bin/env python3
"""
analyze_ensemble.py -- poolt die Defektcluster ALLER Laeufe eines Ensembles
und prueft auf Power-Law F(n) ~ n^-S.

  1) jeder Lauf -> verlagerte Atome -> Cluster (Union-Find)
  2) alle Clustergroessen zusammenwerfen
  3) log-log-Histogramm + CCDF + MLE-Exponent (Clauset-Schaetzer)

Aufruf:  python3 analyze_ensemble.py [outdir]
"""
import sys, os, glob, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

args = [a for a in sys.argv[1:]]
FIX_XMIN = None
if "--xmin" in args:
    i = args.index("--xmin"); FIX_XMIN = float(args[i+1]); del args[i:i+2]
OUTDIR = args[0] if args else "ensemble"
L0 = 1.0
DISP_THRESH = 0.5 * L0
LINK = 1.5 * L0

def cluster_sizes(px, py, link=LINK):
    n = len(px)
    if n == 0: return []
    parent = list(range(n))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    cell = link; buckets = {}
    for i in range(n):
        buckets.setdefault((int(px[i]//cell), int(py[i]//cell)), []).append(i)
    link2 = link*link
    for i in range(n):
        cx, cy = int(px[i]//cell), int(py[i]//cell)
        for ax in (cx-1, cx, cx+1):
            for ay in (cy-1, cy, cy+1):
                for j in buckets.get((ax, ay), ()):
                    if j <= i: continue
                    if (px[i]-px[j])**2 + (py[i]-py[j])**2 < link2:
                        ra, rb = find(i), find(j)
                        if ra != rb: parent[ra] = rb
    roots = {}
    for i in range(n):
        r = find(i); roots[r] = roots.get(r, 0) + 1
    return list(roots.values())

# ---- alle Laeufe einsammeln ----
all_sizes = []
files = sorted(glob.glob(os.path.join(OUTDIR, "run_*_state_final.csv")))
for fn in files:
    d = np.genfromtxt(fn, delimiter=",", names=True)
    if d.size == 0: continue
    x, y, disp = np.atleast_1d(d["x"]), np.atleast_1d(d["y"]), np.atleast_1d(d["disp"])
    moved = disp > DISP_THRESH
    all_sizes += cluster_sizes(x[moved], y[moved])

sizes = np.array(sorted(all_sizes, reverse=True), dtype=float)
print(f"Laeufe: {len(files)}   Cluster gesamt: {len(sizes)}")
if len(sizes) == 0:
    sys.exit("keine Cluster gefunden")
print(f"Groessen: min={sizes.min():.0f} max={sizes.max():.0f} median={np.median(sizes):.0f}")

# ---- MLE-Exponent (Clauset, kontinuierlich): S = 1 + n / sum ln(x/xmin) ----
def mle_alpha(data, xmin):
    d = data[data >= xmin]
    if len(d) < 5: return None, 0
    return 1.0 + len(d) / np.sum(np.log(d / xmin)), len(d)

# xmin per einfachem KS-Scan ueber Kandidaten
def ccdf(d):
    s = np.sort(d); c = 1.0 - np.arange(len(s)) / len(s)
    return s, c
def ks_for(xmin):
    a, ndat = mle_alpha(sizes, xmin)
    if a is None: return None
    d = sizes[sizes >= xmin]; xs, emp = ccdf(d)
    theo = (xs / xmin) ** (1.0 - a)
    return a, np.max(np.abs(emp - theo)), ndat

# Schwanz-Schutz: mind. 5% der Cluster ODER 50 Punkte muessen ueber xmin bleiben
min_keep = max(10, min(max(50, int(0.05*len(sizes))), len(sizes)//3))
if FIX_XMIN is not None:
    res = ks_for(FIX_XMIN)
    alpha, ks, ndat = (res[0], res[1], res[2]) if res else (None,1e9,0)
    xmin = FIX_XMIN
    print(f"[festes xmin={xmin:.0f}]")
else:
    best = (None, None, 1e9, 0)
    for xmin in np.unique(sizes):
        if (sizes >= xmin).sum() < min_keep: break   # nicht in den Schwanz fliehen
        res = ks_for(xmin)
        if res is None: continue
        a, ks, ndat = res
        if ks < best[2]: best = (a, xmin, ks, ndat)
    alpha, xmin, ks, ndat = best

print(f"MLE: S={alpha:.2f}  xmin={xmin:.0f}  (n>=xmin: {ndat})  KS={ks:.3f}")
if ks > 0.10:
    print("  WARNUNG: KS>0.10 -> schlechter Power-Law-Fit (Verteilung evtl. kein Power Law)")
print(f"Referenz (3D-Wolfram, Sand et al.): S = 1.63 +/- 0.07  (3D; 2D darf abweichen)")

# pooled Clustergroessen sichern (fuer spaeteres Re-Fitten ohne Re-Clustern)
np.savetxt(os.path.join(OUTDIR,"cluster_sizes.csv"), sizes, fmt="%d", header="size", comments="")

# ---- Plot: Histogramm (log-log) + CCDF mit Fit ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))

vals, counts = np.unique(sizes, return_counts=True)
ax1.loglog(vals, counts, "o", ms=5, color="#185FA5")
ax1.set_xlabel("Clustergroesse n"); ax1.set_ylabel("Haeufigkeit")
ax1.set_title(f"Gepoolte Clustergroessen ({len(sizes)} Cluster, {len(files)} Laeufe)")
ax1.grid(alpha=.3, which="both")

xs, emp = ccdf(sizes)
ax2.loglog(xs, emp, ".", ms=4, color="#185FA5", label="Daten (CCDF)")
if alpha is not None:
    xx = np.array([xmin, sizes.max()])
    yy = (sizes >= xmin).mean() * (xx / xmin) ** (1.0 - alpha)
    ax2.loglog(xx, yy, "-", lw=2, color="#D85A30",
               label=f"Power-Law-Fit  S={alpha:.2f}")
    ax2.axvline(xmin, ls=":", color="gray", lw=1, label=f"xmin={xmin:.0f}")
ax2.set_xlabel("Clustergroesse n"); ax2.set_ylabel("P(N \u2265 n)")
ax2.set_title("CCDF mit MLE-Fit (Clauset-Methode)")
ax2.legend(fontsize=8); ax2.grid(alpha=.3, which="both")

fig.suptitle(f"Ensemble-Auswertung: {OUTDIR}", fontsize=11)
fig.tight_layout()
out = os.path.join(OUTDIR, "ensemble_powerlaw.png")
fig.savefig(out, dpi=130)
print("geschrieben:", out)
