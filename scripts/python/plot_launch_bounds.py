#!/usr/bin/env python3
"""
plot_launch_bounds.py -- visualisiert das __launch_bounds__-Experiment.

Kernaussage: die Occupancy laesst sich per Register-Begrenzung von 45 % auf
66 % heben, der Durchsatz folgt aber NICHT -- der Kraft-Kernel ist
latenz-/bandbreitengebunden. Bei 32 Registern (theor. 100 % Occupancy) kippt
es durch Register-Spilling.

Aufruf:  python3 plot_launch_bounds.py [launch_bounds_results.csv]
Ausgabe: launch_bounds_plot.png
"""
import sys, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CSV = sys.argv[1] if len(sys.argv) > 1 else "launch_bounds_results.csv"

rows = list(csv.DictReader(open(CSV)))
# lb4 ist identisch zu base (gleiche Register) -> fuer die Aussage base/lb6/lb8 zeigen
order = ["base", "lb6", "lb8"]
rows = [next(r for r in rows if r["variant"] == v) for v in order]

labels   = [r["variant"] for r in rows]
regs     = [int(r["registers"]) for r in rows]
occ      = [float(r["occupancy_pct"]) for r in rows]
occ_src  = [r["occupancy_source"] for r in rows]
thr_abs  = [float(r["throughput_matom_s"]) for r in rows]
base_thr = thr_abs[0]
thr_rel  = [100.0 * t / base_thr for t in thr_abs]   # Durchsatz relativ zu base

# --- helles Theme (einheitlich mit den übrigen Abbildungen) -----------------
BG, FG, GRID = "white", "#1a1a1a", "#c8c8c8"
C_OCC, C_THR = "#185FA5", "#D85A30"     # blau = Occupancy, orange = Durchsatz
plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
    "text.color": FG, "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG,
    "axes.edgecolor": GRID, "font.size": 11,
})

fig, ax = plt.subplots(figsize=(8.4, 5.0))
x = np.arange(len(labels)); w = 0.38

b1 = ax.bar(x - w/2, occ,     w, color=C_OCC, label="Occupancy (%)")
b2 = ax.bar(x + w/2, thr_rel, w, color=C_THR, label="Durchsatz (% von base)")

ax.axhline(100, color=FG, lw=0.8, ls="--", alpha=0.4)   # base-Durchsatz-Referenz
ax.set_ylim(0, 115)
ax.set_ylabel("Prozent")
ax.set_title("Occupancy-Hebel (__launch_bounds__): mehr Occupancy \u2260 mehr Durchsatz",
             color=FG, pad=14, fontsize=12.5, weight="bold")

# x-Beschriftung mit Registern + Occupancy-Quelle
xt = []
for lab, reg, src in zip(labels, regs, occ_src):
    tag = "ncu" if src == "ncu" else "theor."
    xt.append(f"{lab}\n{reg} Reg/Thread")
ax.set_xticks(x); ax.set_xticklabels(xt)

# Werte ueber die Balken
for rect, v, src in zip(b1, occ, occ_src):
    note = "" if src == "ncu" else " (theor.)"
    ax.text(rect.get_x()+rect.get_width()/2, v+1.5, f"{v:.0f}%{note}",
            ha="center", va="bottom", color=C_OCC, fontsize=9.5)
for rect, v, a in zip(b2, thr_rel, thr_abs):
    ax.text(rect.get_x()+rect.get_width()/2, v+1.5, f"{v:.0f}%\n{a:.0f} Matom/s",
            ha="center", va="bottom", color=C_THR, fontsize=9.0)

ax.grid(axis="y", color=GRID, lw=0.6, alpha=0.5)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(loc="lower left", framealpha=0.0, fontsize=10)

# Annotation der Pointe
ax.annotate("Occupancy +21 Pp,\nDurchsatz flach",
            xy=(1+w/2, thr_rel[1]), xytext=(1.05, 60),
            color=FG, fontsize=9.5, ha="center",
            arrowprops=dict(arrowstyle="->", color=FG, alpha=0.6))
ax.annotate("32 Reg \u2192 Spilling\n\u221217 % Durchsatz",
            xy=(2+w/2, thr_rel[2]), xytext=(2.0, 55),
            color=FG, fontsize=9.5, ha="center",
            arrowprops=dict(arrowstyle="->", color=FG, alpha=0.6))

fig.tight_layout()
fig.savefig("launch_bounds_plot.png", dpi=200)
print("geschrieben: launch_bounds_plot.png")
