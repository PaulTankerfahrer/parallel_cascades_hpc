#!/usr/bin/env python3
"""
plot_ensemble_scaling.py -- Machbarkeitsstudie "skalierter Ensemble-Lauf".

Zwei Befunde der Voruntersuchung in einer Abbildung:
  Links  -- MPI-Strong-Scaling eines einzelnen 1000x1000-Laufs (1 Mio Atome).
            Reproduziert die superlineare Skalierung (Cache-Effekt) bei 64/128
            Kernen -- derselbe rote Faden wie im Hauptteil.
  Rechts -- Energie-Kalibrierung: verschobene Atome vs. PKA-Energie. Ab ~5000
            perkoliert der Schaden zu EINEM Cluster (groesster Cluster ~ alle
            verschobenen Atome); ab ~15000 ueberrennt er das Gitter. Das ist
            der wissenschaftliche Grund gegen den Grosslauf.

Daten hardcoded aus den Voruntersuchungs-Laeufen (s. Ensemble-Skalierung-Bericht).
Aufruf:  python3 plot_ensemble_scaling.py
Ausgabe: ensemble_scaling.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Messdaten --------------------------------------------------------------
# MPI-Strong-Scaling, 1 Mio Atome (1000x1000), 3000 Schritte
cores    = np.array([1, 16, 32, 64, 128])
speedup  = np.array([1.0, 9.0, 27.1, 87.3, 171.0])

# Energie-Kalibrierung (1000x1000, PKA zentral, 40000 Schritte)
energy    = np.array([3000, 5000, 8000, 15000])
displaced = np.array([4983, 44626, 59347, 67995])
overrun   = np.array([False, False, False, True])   # 15000 ueberrennt das Gitter
percolate = np.array([False, True, True, True])      # groesster Cluster ~ alle Defekte

# --- helles Theme (einheitlich mit den übrigen Abbildungen) -----------------
BG, FG, GRID = "white", "#1a1a1a", "#c8c8c8"
C_COMP, C_COMM, C_IDEAL = "#185FA5", "#D85A30", "#9AA7B2"
plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
    "text.color": FG, "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG,
    "axes.edgecolor": GRID, "font.size": 11,
})

fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.4, 4.8))

# === Links: MPI strong scaling, 1 Mio Atome =================================
axL.plot(cores, cores, ls="--", lw=1.2, color=C_IDEAL, label="ideal ($S=P$)")
axL.plot(cores, speedup, "o-", lw=1.8, ms=7, color=C_COMP, label="gemessen")
# superlinearen Bereich markieren
for c, s in zip(cores, speedup):
    if s > c * 1.05:
        axL.annotate(f"{s:.0f}$\\times$", (c, s), textcoords="offset points",
                     xytext=(-4, 8), color=C_COMP, fontsize=9.5, ha="right")
axL.set_xscale("log", base=2); axL.set_yscale("log", base=2)
axL.set_xticks(cores); axL.set_xticklabels(cores)
axL.set_yticks([1, 4, 16, 64, 171]); axL.set_yticklabels([1, 4, 16, 64, 171])
axL.set_xlabel("MPI-Ranks (Kerne)"); axL.set_ylabel("Speedup $S = T(1)/T(P)$")
axL.set_title("MPI-Strong-Scaling, 1 Mio Atome\n(superlinear ab 64 Kernen, Cache-Effekt)",
              color=FG, fontsize=11.5, pad=10)
axL.grid(True, which="both", color=GRID, lw=0.6, alpha=0.5)
axL.set_axisbelow(True)
for s in ("top", "right"):
    axL.spines[s].set_visible(False)
axL.legend(loc="upper left", framealpha=0.0, fontsize=10)

# === Rechts: Energie-Kalibrierung ==========================================
colors = [C_COMM if (p or o) else C_COMP for p, o in zip(percolate, overrun)]
x = np.arange(len(energy))
bars = axR.bar(x, displaced, color=colors, width=0.6)
for rect, d, p, o in zip(bars, displaced, percolate, overrun):
    tag = ""
    if o:
        tag = "\nGitter\nueberrannt"
    elif p:
        tag = "\n1 Blob"
    axR.text(rect.get_x()+rect.get_width()/2, d+1500, f"{d:,}".replace(",", ".")+tag,
             ha="center", va="bottom", color=FG, fontsize=8.5)
axR.set_xticks(x); axR.set_xticklabels([f"{e}" for e in energy])
axR.set_ylim(0, 82000)
axR.set_xlabel("PKA-Energie (dimensionslos)")
axR.set_ylabel("verschobene Atome")
axR.set_title("Energie-Kalibrierung: ein Grosslauf liefert\nkeine breite Verteilung, sondern einen Blob",
              color=FG, fontsize=11.5, pad=10)
# Legende von Hand
from matplotlib.patches import Patch
axR.legend(handles=[
    Patch(color=C_COMP, label="breite Verteilung (gewuenscht)"),
    Patch(color=C_COMM, label="perkoliert / ueberrannt"),
], loc="upper left", framealpha=0.0, fontsize=9.5)
axR.grid(axis="y", color=GRID, lw=0.6, alpha=0.5)
axR.set_axisbelow(True)
for s in ("top", "right"):
    axR.spines[s].set_visible(False)

fig.tight_layout()
fig.savefig("ensemble_scaling.png", dpi=200)
print("geschrieben: ensemble_scaling.png")
