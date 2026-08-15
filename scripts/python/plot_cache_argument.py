#!/usr/bin/env python3
"""
plot_cache_argument.py -- Cache-Argument: Datensatz/Rank [MB] + Effizienz [%] vs. Kerne.

Dual-Achsen-Plot:
  - Links: Datensatz pro Rank in MB (log-Skala)
  - Rechts: parallele Effizienz in %
  - Horizontale Linie: L3-Cache-Schwelle (4 MB = 32 MB L3 / 8 Kerne pro EPYC-Chip)

Daten aus mpi_strong_scaling_1node.csv + mpi_strong_scaling_2node.csv.
Erzeugt: cache_argument.png
"""
import sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

f1 = sys.argv[1] if len(sys.argv) > 1 else "mpi_strong_scaling_1node.csv"
f2 = sys.argv[2] if len(sys.argv) > 2 else "mpi_strong_scaling_2node.csv"

def load(fn):
    d = np.genfromtxt(fn, delimiter=",", names=True)
    r = np.atleast_1d(d["ranks"]).astype(float)
    w = np.atleast_1d(d["walltime"])
    o = np.argsort(r)
    return r[o], w[o]

r1, w1 = load(f1)
T1 = w1[0]  # Referenz: 1 Kern
sp1 = T1 / w1
eff1 = sp1 / (r1 / r1[0]) * 100

try:
    r2, w2 = load(f2)
    sp2 = T1 / w2
    eff2 = sp2 / (r2 / r1[0]) * 100
    has2 = True
except Exception:
    has2 = False

# Datensatz pro Rank: 76 Byte/Atom × N_total / ranks → in MB
# N_total = 1999396 für die 1-Knoten-Skalierung
N_TOTAL  = 1_999_396
BYTES_PER_ATOM = 76   # px,py,vx,vy,fx,fy,nbr[12],cell → 8+8+8+8+8+8+12*4+4 = 96... aber ~76 Nutzlast
MB_PER_RANK_1 = N_TOTAL * BYTES_PER_ATOM / r1 / (1024**2)

L3_THRESHOLD = 4.0   # MB pro Rank (32 MB L3 / 8 Kerne/EPYC-Chip)

fig, ax1 = plt.subplots(figsize=(10, 6))

# ---- Datensatz/Rank (linke Achse) -------------------------------------------
color_ds = "#185FA5"
color_ef = "#D85A30"

ax1.semilogy(r1, MB_PER_RANK_1, "s--", color=color_ds, ms=8, lw=1.5,
             label="Datensatz/Rank (1 Knoten)")
ax1.axhline(L3_THRESHOLD, ls=":", color="black", lw=1.8,
            label=f"L3-Schwelle = {L3_THRESHOLD:.0f} MB/Rank\n(32 MB ÷ 8 Kerne/EPYC-Chip)")
ax1.fill_between(r1, MB_PER_RANK_1, L3_THRESHOLD,
                 where=(MB_PER_RANK_1 > L3_THRESHOLD),
                 alpha=0.12, color=color_ds, label="Cache-miss-Regime")
ax1.fill_between(r1, MB_PER_RANK_1, L3_THRESHOLD,
                 where=(MB_PER_RANK_1 <= L3_THRESHOLD),
                 alpha=0.15, color="green", label="Cache-hit-Regime")

ax1.set_xlabel("MPI-Ranks (Kerne)", fontsize=11)
ax1.set_ylabel("Datensatz pro Rank  [MB, log]", color=color_ds, fontsize=11)
ax1.tick_params(axis="y", labelcolor=color_ds)
ax1.set_xscale("log", base=2)
ax1.set_xticks(r1)
ax1.set_xticklabels([f"{int(r)}" for r in r1])
ax1.grid(alpha=0.25, which="both")

for r, mb in zip(r1, MB_PER_RANK_1):
    ax1.annotate(f"{mb:.1f} MB", (r, mb), textcoords="offset points",
                 xytext=(0, 7), ha="center", fontsize=8, color=color_ds)

# ---- Effizienz (rechte Achse) ------------------------------------------------
ax2 = ax1.twinx()
ax2.plot(r1, eff1, "o-", color=color_ef, ms=9, lw=2, label="Effizienz 1 Knoten")
if has2:
    ax2.plot(r2, eff2, "D--", color="#7A3B8A", ms=8, lw=2, label="Effizienz 2 Knoten")
ax2.axhline(100, ls="--", color="gray", lw=1)
for r, e in zip(r1, eff1):
    ax2.annotate(f"{e:.0f}%", (r, e), textcoords="offset points",
                 xytext=(0, -14), ha="center", fontsize=8.5, color=color_ef,
                 fontweight="bold")
ax2.set_ylabel("Parallele Effizienz  [%]", color=color_ef, fontsize=11)
ax2.tick_params(axis="y", labelcolor=color_ef)
ax2.set_ylim(0, 160)

# Superlinear-Hinweis
super_idx = np.argmax(eff1)
if eff1[super_idx] > 110:
    ax2.annotate("Superlinear!\n(Cache-Effekt)",
                 xy=(r1[super_idx], eff1[super_idx]),
                 xytext=(r1[super_idx] / 2, eff1[super_idx] + 15),
                 arrowprops=dict(arrowstyle="->", color=color_ef, lw=1.5),
                 fontsize=9, color=color_ef, ha="center")

# ---- Legenden zusammenführen -------------------------------------------------
# Außerhalb der Achse platziert, damit sie nicht die Superlinear-Annotation
# oder die Datenpunkte im rechten oberen Bereich überdeckt.
lines1, labs1 = ax1.get_legend_handles_labels()
lines2, labs2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labs1 + labs2, fontsize=8.5,
           loc="upper left", bbox_to_anchor=(1.12, 1.0), borderaxespad=0)
ax1.set_title("Cache-Argument: L3-Schwelle erklärt superlinearen Speedup\n"
              f"(N={N_TOTAL:,} Atome, {BYTES_PER_ATOM} B/Atom → {N_TOTAL*BYTES_PER_ATOM/1e6:.0f} MB gesamt)",
              fontsize=11)

fig.tight_layout()
fig.savefig("cache_argument.png", dpi=150, bbox_inches="tight")
print("geschrieben: cache_argument.png")
# Crossing point
crossing = r1[MB_PER_RANK_1 <= L3_THRESHOLD]
if len(crossing):
    print(f"Cache-Hit-Regime ab {int(crossing[0])} Ranks "
          f"(Datensatz/Rank = {MB_PER_RANK_1[np.where(r1==crossing[0])[0][0]]:.1f} MB ≤ {L3_THRESHOLD} MB)")
