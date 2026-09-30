#!/usr/bin/env python3
"""
plot_cuda_opt.py -- publication-ready Plot der GPU-Optimierung.

Kernaussage in zwei Panels:
  links  : die Occupancy variiert stark (45 % <-> 66 %, per launch_bounds)
  rechts : der Durchsatz bleibt flach (alle Varianten ~3800 Matom/s)
-> weder Occupancy noch Launch-Overhead (CUDA Graphs) heben die Leistung;
   der Kernel ist speicher-latenz-/bandbreitengebunden.

Aufruf:  python3 plot_cuda_opt.py [cuda_opt_results.csv]
Ausgabe: cuda_opt_plot.png (300 dpi)
"""
import sys, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CSV = sys.argv[1] if len(sys.argv) > 1 else "cuda_opt_results.csv"
rows = list(csv.DictReader(open(CSV)))
order = ["base", "lb", "graph", "lbgraph"]
rows = [next(r for r in rows if r["variant"] == v) for v in order]
labels = ["base", "lb\n(launch_bounds)", "graph\n(CUDA Graphs)", "lb+graph"]
occ    = [float(r["occupancy_pct"]) for r in rows]
regs   = [int(r["registers"]) for r in rows]
thr    = [float(r["throughput_matom_s"]) for r in rows]
kf     = [float(r["kforce_us"]) for r in rows]

BG, FG, GRID = "white", "#1a1a1a", "#c8c8c8"
C_OCC, C_THR, C_REF = "#185FA5", "#D85A30", "#8a95b0"
plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
    "text.color": FG, "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG,
    "axes.edgecolor": GRID, "font.size": 11,
})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.8))
x = np.arange(len(labels))

# --- links: Occupancy ---
b1 = ax1.bar(x, occ, 0.62, color=C_OCC)
ax1.set_ylim(0, 78); ax1.set_ylabel("Occupancy (%)  — ncu")
ax1.set_title("Occupancy variiert stark", color=FG, fontsize=12, weight="bold")
ax1.set_xticks(x); ax1.set_xticklabels(labels, fontsize=9)
for r, v, reg in zip(b1, occ, regs):
    ax1.text(r.get_x()+r.get_width()/2, v+1.5, f"{v:.0f}%", ha="center", va="bottom",
             color=C_OCC, fontsize=10, weight="bold")
    ax1.text(r.get_x()+r.get_width()/2, 3, f"{reg} Reg", ha="center", va="bottom",
             color="white", fontsize=8.5)
ax1.text(0.5, 0.93, "k_force-Laufzeit dabei konstant: 307–313 µs",
         transform=ax1.transAxes, ha="center", color=FG, fontsize=9.5, style="italic")

# --- rechts: Durchsatz (Achse bei 0 -> ehrliche Flachheit) ---
b2 = ax2.bar(x, thr, 0.62, color=C_THR)
ax2.set_ylim(0, 4300); ax2.set_ylabel("Durchsatz (Matom/s)")
ax2.set_title("Durchsatz bleibt flach", color=FG, fontsize=12, weight="bold")
ax2.set_xticks(x); ax2.set_xticklabels(labels, fontsize=9)
ax2.axhline(thr[0], color=C_REF, lw=0.9, ls="--", alpha=0.6)
for r, v in zip(b2, thr):
    d = 100*(v-thr[0])/thr[0]
    ax2.text(r.get_x()+r.get_width()/2, v+60, f"{v:.0f}\n({d:+.1f}%)", ha="center",
             va="bottom", color=C_THR, fontsize=9)

for ax in (ax1, ax2):
    ax.grid(axis="y", color=GRID, lw=0.6, alpha=0.5); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)

fig.suptitle("GPU-Optimierung (fp64, 2 Mio Atome, A100): Occupancy hebbar — Leistung nicht",
             color=FG, fontsize=13, weight="bold", y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig("cuda_opt_plot.png", dpi=300)
print("geschrieben: cuda_opt_plot.png")
