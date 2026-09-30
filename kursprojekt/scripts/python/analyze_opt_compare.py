#!/usr/bin/env python3
"""
analyze_opt_compare.py -- Strong-Scaling-Varianten vergleichen.
Legt die Speedup-Kurven uebereinander (Optimierung aendert die Kurve oder
nicht?) und zeigt die Zeitersparnis pro Rankzahl gegen die Baseline.

Aufruf:
  python3 analyze_opt_compare.py off.csv "AUS" on.csv "non-blocking + resort"
  (weitere Varianten anhaengbar: ... v2.csv "nur resort" )
"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

args = sys.argv[1:]
if len(args) < 4 or len(args) % 2 != 0:
    print("Aufruf: analyze_opt_compare.py base.csv 'Label' var1.csv 'Label1' [var2.csv 'Label2' ...]")
    sys.exit(1)

base_csv, base_lbl = args[0], args[1]
variants = [(args[i], args[i+1]) for i in range(2, len(args), 2)]
colors  = ["#E08A1E", "#185FA5", "#2E8B57", "#A0408A"]
markers = ["s", "^", "o", "D"]

def load(fn):
    d = np.genfromtxt(fn, delimiter=",", names=True)
    r = np.atleast_1d(d["ranks"]).astype(float)
    w = np.atleast_1d(d["walltime"])
    o = np.argsort(r); return r[o], w[o]

rb, wb = load(base_csv); spb = wb[0]/wb
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.7))
x  = np.arange(len(rb)); nv = len(variants); wbar = 0.8/max(nv, 1)

ax1.plot(rb, rb, "--", color="gray", lw=1, label="ideal")
ax1.plot(rb, spb, "o-", color="#999999", label=base_lbl)

print(f"Baseline = {base_lbl}")
for k, (fn, lbl) in enumerate(variants):
    rv, wv = load(fn); spv = wv[0]/wv
    ax1.plot(rv, spv, marker=markers[k % 4], ls="-", color=colors[k % 4], label=lbl)
    gains = []
    for i, rr in enumerate(rb):
        j = np.where(rv == rr)[0]
        gains.append(100*(wb[i]-wv[j[0]])/wb[i] if len(j) else np.nan)
    ax2.bar(x + (k-(nv-1)/2)*wbar, gains, wbar, color=colors[k % 4], label=lbl)
    g = np.array([g for g in gains if not np.isnan(g)])
    print(f"  {lbl}: mittlerer Gewinn {g.mean():+.1f} %  (Spanne {g.min():+.1f} .. {g.max():+.1f})")

ax1.set_xscale("log", base=2); ax1.set_yscale("log", base=2)
ax1.set_xlabel("Kerne"); ax1.set_ylabel("Speedup T(1)/T(p)")
ax1.set_title("Strong Scaling: Varianten im Vergleich")
ax1.legend(fontsize=8); ax1.grid(alpha=.3, which="both")

ax2.axhspan(-2, 2, color="gray", alpha=.15, label="Rauschband ±2%")
ax2.axhline(0, color="gray", lw=.8)
ax2.set_xticks(x); ax2.set_xticklabels([str(int(v)) for v in rb])
ax2.set_xlabel("Kerne"); ax2.set_ylabel(f"Zeitersparnis vs. {base_lbl} [%]")
ax2.set_title("Gewinn pro Rankzahl (>0 = schneller)")
ax2.legend(fontsize=8); ax2.grid(alpha=.3, axis="y")

fig.suptitle("Optimierungs-Vergleich (Strong Scaling)", fontsize=11)
fig.tight_layout(); fig.savefig("opt_compare.png", dpi=130)
print("geschrieben: opt_compare.png")
