#!/usr/bin/env python3
"""
overlay_ccdf.py -- legt zwei Cluster-Verteilungen (z. B. mit/ohne Healing)
in EINER CCDF-Grafik uebereinander und fittet beide mit DEMSELBEN xmin.
Zeigt visuell, dass der Power-Law-Exponent robust ist.

Braucht die cluster_sizes.csv, die analyze_ensemble.py jetzt schreibt.
Aufruf:
  python3 overlay_ccdf.py ohne/cluster_sizes.csv "ohne Healing" \
                          mit/cluster_sizes.csv  "mit Healing" [xmin]
"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

a = sys.argv[1:]
XMIN = 4.0
if len(a) in (5, 7) and a[-1].replace(".", "", 1).isdigit():
    XMIN = float(a[-1]); a = a[:-1]
pairs = [(a[i], a[i+1]) for i in range(0, len(a), 2)]
colors = ["#185FA5", "#D85A30", "#2E8B57"]

def load(fn):
    s = np.atleast_1d(np.genfromtxt(fn, delimiter=",", skip_header=1))
    return np.sort(s[s > 0].astype(float))

def ccdf(s):
    xs = np.sort(s); c = 1.0 - np.arange(len(xs))/len(xs); return xs, c

def mle(s, xmin):
    d = s[s >= xmin]
    if len(d) < 5: return None
    return 1.0 + len(d)/np.sum(np.log(d/xmin))

fig, ax = plt.subplots(figsize=(7.2, 5.4))
print(f"festes xmin = {XMIN:.0f}")
for k, (fn, lbl) in enumerate(pairs):
    s = load(fn); xs, c = ccdf(s); alpha = mle(s, XMIN)
    col = colors[k % len(colors)]
    ax.loglog(xs, c, ".", ms=4, color=col, alpha=.55)
    xx = np.array([XMIN, s.max()])
    yy = (s >= XMIN).mean() * (xx/XMIN)**(1.0-alpha)
    ax.loglog(xx, yy, "-", lw=2, color=col,
              label=f"{lbl}: S={alpha:.2f}  ({len(s)} Cluster)")
    print(f"  {lbl}: S={alpha:.2f}  (n={len(s)})")

ax.axvline(XMIN, ls=":", color="gray", lw=1)
ax.text(XMIN*1.1, 1.3e-3, f"xmin={XMIN:.0f}", color="gray", fontsize=9)
ax.set_xlabel("Clustergroesse n"); ax.set_ylabel("P(N \u2265 n)")
ax.set_title("Cluster-Verteilung: mit vs. ohne Healing")
ax.legend(fontsize=10); ax.grid(alpha=.3, which="both")
fig.tight_layout(); fig.savefig("ccdf_overlay.png", dpi=140)
print("geschrieben: ccdf_overlay.png")
