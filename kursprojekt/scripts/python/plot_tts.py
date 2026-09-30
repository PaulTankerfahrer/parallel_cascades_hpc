#!/usr/bin/env python3
"""
plot_tts.py -- Time-to-Solution Vergleich: CPU (1 Kern), CPU (64 Kerne),
               Tesla T4, A100 fp64, A100 fp32.

Erzeugt: tts.png
Daten sind Messwerte aus den Benchmark-Ausgaben (hardcoded).
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Messwerte: N=~2M Atome, 3000 Schritte
systems = [
    ("CPU\n1 Kern",    1152.9, "#9AA7B2"),
    ("CPU\n64 Kerne",  22.5,   "#3E8AC0"),
    ("Tesla T4\n(fp64)", 19.09, "#E8915F"),
    ("A100\n(fp64)",   1.5765, "#185FA5"),
    ("A100\n(fp32)",   0.999,  "#D85A30"),
]

labels  = [s[0] for s in systems]
times   = [s[1] for s in systems]
colors  = [s[2] for s in systems]
T1      = 1152.9
speedups = [T1 / t for t in times]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))

# ---- (1) Laufzeit (log-Skala) -----------------------------------------------
bars = ax1.bar(range(len(labels)), times, color=colors, zorder=3)
ax1.set_yscale("log")
ax1.set_xticks(range(len(labels)))
ax1.set_xticklabels(labels, fontsize=10)
ax1.set_ylabel("Laufzeit [s]  (log-Skala, kleiner = besser)")
ax1.set_title("Time-to-Solution\n(N ≈ 2 Mio Atome, 3000 Schritte)")
ax1.grid(axis="y", alpha=0.35, which="both", zorder=0)
for i, (t, sp) in enumerate(zip(times, speedups)):
    label = f"{t:.3f} s" if t < 10 else f"{t:.1f} s"
    ax1.text(i, t * 1.6, label, ha="center", va="bottom", fontsize=8.5, fontweight="bold")

# Speedup-Anmerkungen relativ zu CPU 1 Kern.
# Position ist relativ zum Balken-TOP (nicht zur Achsen-Unterkante) verankert,
# damit das Label auch bei sehr kurzen Balken (A100) innerhalb des Balkens
# bleibt und nicht mit den x-Achsen-Beschriftungen kollidiert.
for i, sp in enumerate(speedups):
    if sp > 1.5:
        ax1.text(i, times[i] / 1.15, f"×{sp:.0f}", ha="center", va="top",
                 fontsize=8, color="white", fontweight="bold")

# ---- (2) Speedup-Balken -------------------------------------------------------
bars2 = ax2.bar(range(len(labels)), speedups, color=colors, zorder=3)
ax2.set_xticks(range(len(labels)))
ax2.set_xticklabels(labels, fontsize=10)
ax2.set_ylabel("Speedup gegenüber 1 CPU-Kern")
ax2.set_title("Speedup  T(1-Kern) / T(System)")
ax2.grid(axis="y", alpha=0.35, zorder=0)
for i, sp in enumerate(speedups):
    ax2.text(i, sp + 5, f"{sp:.0f}×", ha="center", va="bottom", fontsize=9, fontweight="bold")

fig.suptitle("GPU-Beschleunigung der Kollisionskaskaden-Simulation", fontsize=13)
fig.tight_layout()
fig.savefig("tts.png", dpi=150)
print("geschrieben: tts.png")
print(f"A100 fp64 Speedup vs 64-Kern-CPU: {22.5/1.5765:.1f}x  (= {1.5765:.3f}s vs {22.5:.1f}s)")
