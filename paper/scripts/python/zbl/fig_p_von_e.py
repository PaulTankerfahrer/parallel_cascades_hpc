#!/usr/bin/env python3
"""fig_p_von_e.py -- Defektwahrscheinlichkeit P(E) aus kalibrierung_ed.py.

Aufruf:
  .venv-powerlaw/bin/python scripts/python/zbl/fig_p_von_e.py \
      results/zbl/kalibrierung/p_von_e.npz results/images/D2_p_von_e.png
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CAT = ["#2a78d6", "#eb6834"]
INK, INK2, SURFACE, GRID = "#1f1f1e", "#5c5c58", "#fcfcfb", "#ececE6"

d = np.load(sys.argv[1])
E = d["E"]
fig, ax = plt.subplots(figsize=(7.2, 4.0), facecolor=SURFACE)
ax.set_facecolor(SURFACE)
ax.grid(True, color=GRID, lw=0.8)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(INK2)
ax.tick_params(colors=INK2)

ax.step(E, d["r12"], where="mid", color=CAT[0], lw=2, label="alt (r⁻¹²)")
ax.step(E, d["zbl_eps0.4"], where="mid", color=CAT[1], lw=2, label="neu (ZBL, ε = 0,4 eV)")
ax.axhspan(0.4, 0.6, xmin=0, xmax=1, color=GRID, zorder=0)
ax.text(44, 0.80, "grau: Plateau bei P ≈ 0,5 –\ndie Hälfte der Richtungen braucht\n3–4× mehr Energie",
        color=INK2, fontsize=9, va="center")
ax.set_xscale("log")
ax.set_xlim(40, 1500)
ticks = [50, 100, 200, 500, 1000]
ax.set_xticks(ticks)
ax.set_xticklabels([str(t) for t in ticks])
ax.minorticks_off()
ax.set_ylim(-0.02, 1.05)
ax.set_xlabel("PKA-Energie E [Modelleinheiten]", color=INK)
ax.set_ylabel("Anteil der Richtungen mit\nbleibendem Leerplatz, P(E)", color=INK)
ax.set_title("Verlagerungswahrscheinlichkeit, 31 Richtungen, 40×40-Gitter",
             loc="left", fontsize=10.5, color=INK)
ax.legend(frameon=False, loc="lower right", fontsize=9)
fig.tight_layout()
fig.savefig(sys.argv[2], dpi=160, facecolor=SURFACE)
print("->", sys.argv[2])
