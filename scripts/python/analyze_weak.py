#!/usr/bin/env python3
"""
analyze_weak.py -- Weak-Scaling-Auswertung.
Bei perfektem Weak Scaling bleibt die Wandzeit konstant (jeder Rank hat gleich
viel zu tun). Der Anstieg zeigt den Kommunikations-/Overhead-Anteil.
Aufruf:  python3 analyze_weak.py weak.csv
"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fn = sys.argv[1] if len(sys.argv) > 1 else "weak.csv"
d = np.genfromtxt(fn, delimiter=",", names=True)
r = np.atleast_1d(d["ranks"]).astype(float)
wall = np.atleast_1d(d["walltime"])
comm = np.atleast_1d(d["comm"])
o = np.argsort(r); r, wall, comm = r[o], wall[o], comm[o]
eff = wall[0] / wall * 100        # Weak-Effizienz: t(1)/t(P)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

ax1.plot(r, wall, "o-", color="#185FA5", label="Wandzeit")
ax1.axhline(wall[0], ls="--", color="gray", label="ideal (flach)")
ax1.set_xscale("log", base=2)
ax1.set_xlabel("Kerne"); ax1.set_ylabel("Wandzeit [s]")
ax1.set_title("Weak Scaling: Wandzeit (flach = ideal)")
ax1.legend(fontsize=9); ax1.grid(alpha=.3, which="both")
ax1.set_ylim(0, max(wall)*1.15)

ax2.plot(r, eff, "s-", color="#2E8B57")
ax2.axhline(100, ls="--", color="gray")
ax2.set_xscale("log", base=2)
ax2.set_xlabel("Kerne"); ax2.set_ylabel("Weak-Effizienz [%]")
ax2.set_title("Effizienz t(1)/t(P)")
ax2.grid(alpha=.3, which="both"); ax2.set_ylim(0, 110)

fig.suptitle(f"Weak Scaling: {fn}", fontsize=11)
fig.tight_layout()
out = fn.replace(".csv", "") + "_weak.png"
fig.savefig(out, dpi=130)
print("geschrieben:", out)
print(f"Weak-Effizienz bei {int(r[-1])} Kernen: {eff[-1]:.0f} %")
