#!/usr/bin/env python3
"""
plot_ncu.py -- ncu-Profiling-Auswertung des k_force-Kernels.

Drei Teilgrafiken:
  1) Kernel-Breakdown pro Schritt: CPU vs. A100 (fp64)
  2) SpeedOfLight-Roofline: T4 vs. A100 (Compute % und Memory %)
  3) Occupancy-Analyse: Registerbeschränkung auf A100

Daten hardcoded aus cuda_bench_33375548.out und tesdla_t4_cuda_run.txt.
Erzeugt: ncu_profiling.png
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ---- Kernel-Breakdown (us/Schritt) ----------------------------------------
# CPU: Phasen aus mpi_strong_scaling_1node.csv @ 1 Rank (in Sekunden, × 1e6 / 5000 Schritte)
cpu_phases = {
    "force":     1026.3953 * 1e6 / 5000,
    "cells":     20.3980   * 1e6 / 5000,
    "hash":      67.5829   * 1e6 / 5000,
    "halo":      0.0025    * 1e6 / 5000,
    "migrate":   5.2422    * 1e6 / 5000,
    "energy":    0.6445    * 1e6 / 5000,
}
# A100: Messwerte aus cuda_bench_33375548.out (Median, 3000 Schritte)
gpu_phases = {
    "force":       259.27,
    "build_cells": 47.36,
    "kick_drift":  133.09,
    "kick":        86.52,
    "clear_cells": 9.42,
}

CPU_COLORS = ["#185FA5", "#3E8AC0", "#7FB2D8", "#D85A30", "#E8915F", "#9AA7B2"]
GPU_COLORS = ["#185FA5", "#3E8AC0", "#7FB2D8", "#D85A30", "#E8915F"]

# ---- SpeedOfLight-Daten -------------------------------------------------------
# A100 (ncu_report.csv: sm__throughput, gpu__compute_memory_throughput)
a100_compute = 45.17
a100_memory  = 37.47
# T4 (tesdla_t4_cuda_run.txt ncu-Output)
t4_compute   = 87.46
t4_memory    = 7.90

# ---- Occupancy-Daten ----------------------------------------------------------
# A100: 62 Register/Thread, 256 Threads/Block, 108 SMs
# Max concurrente Blocks/SM durch Register = floor(65536 / (62 * 256)) = 4
# Max Waves = 108 * 4 = 432 Blocks = 432 * 256 = 110592 Threads
# Max Warps/SM: A100 hat 64 Warps/SM -> 4 Blocks * 8 Warps = 32 Warps = 50% Okkupanz
occ_data = {
    "A100\n(fp64)":  50.0,
    "A100\n(fp32)":  75.0,  # fp32 braucht ~48 Register/Thread → 5 Blocks/SM
}

fig, axes = plt.subplots(1, 3, figsize=(16, 6))

# ---- Plot 1: Kernel-Breakdown -------------------------------------------------
ax = axes[0]
cpu_vals = list(cpu_phases.values())
cpu_keys = list(cpu_phases.keys())
gpu_vals = list(gpu_phases.values())
gpu_keys = list(gpu_phases.keys())

cpu_total = sum(cpu_vals)
gpu_total = sum(gpu_vals)

# Normierte Anteile als Stacked Bar
cpu_fracs = [v / cpu_total * 100 for v in cpu_vals]
gpu_fracs = [v / gpu_total * 100 for v in gpu_vals]

x = np.array([0, 1])
# CPU bar
bottom = 0
for i, (k, v, c) in enumerate(zip(cpu_keys, cpu_fracs, CPU_COLORS)):
    ax.bar(0, v, bottom=bottom, color=c, label=k, width=0.5)
    if v > 4:
        ax.text(0, bottom + v/2, f"{k}\n{v:.0f}%", ha="center", va="center",
                fontsize=7.5, color="white" if v > 15 else "black")
    bottom += v
# GPU bar
bottom = 0
for i, (k, v, c) in enumerate(zip(gpu_keys, gpu_fracs, GPU_COLORS)):
    ax.bar(1, v, bottom=bottom, color=c, width=0.5)
    if v > 4:
        ax.text(1, bottom + v/2, f"{k}\n{v:.0f}%", ha="center", va="center",
                fontsize=7.5, color="white" if v > 20 else "black")
    bottom += v

ax.set_xticks([0, 1])
ax.set_xticklabels([f"CPU (1 Kern)\n{cpu_total/1e6*1e6:.0f} µs/Schritt",
                    f"A100 (fp64)\n{gpu_total:.0f} µs/Schritt"], fontsize=9)
ax.set_ylabel("Zeitanteil [%]")
ax.set_ylim(0, 120)
ax.set_title("Kernel-Breakdown\npro Zeitschritt")
ax.grid(axis="y", alpha=0.3)

# Speedup annotation -- oberhalb beider Balken platziert, damit die Linie
# nicht durch die "build_cells"-Beschriftung im GPU-Balken läuft.
ax.annotate(f"", xy=(1.3, 106), xytext=(0.7, 106),
            arrowprops=dict(arrowstyle="->", color="darkgreen", lw=2))
speedup_force = (1026.3953/5000*1e6) / 259.27
ax.text(0.5, 112, f"Force: ×{speedup_force:.0f}", ha="center", va="center",
        fontsize=8, color="darkgreen")

# ---- Plot 2: SpeedOfLight (Roofline) -----------------------------------------
ax = axes[1]
gpus  = ["Tesla T4", "NVIDIA A100"]
comp  = [t4_compute, a100_compute]
mem   = [t4_memory,  a100_memory]
x_pos = [0, 1]

bar_w = 0.35
b1 = ax.bar([x - bar_w/2 for x in x_pos], comp, bar_w,
            label="Compute (SM) %", color=["#D85A30", "#D85A30"], zorder=3)
b2 = ax.bar([x + bar_w/2 for x in x_pos], mem, bar_w,
            label="Memory %", color=["#3E8AC0", "#3E8AC0"], alpha=0.6, zorder=3)

ax.axhline(80, ls="--", color="gray", lw=1.2, label="Engpass-Schwelle (80%)")
ax.set_xticks(x_pos)
ax.set_xticklabels(gpus, fontsize=10)
ax.set_ylabel("% des Peak-Durchsatzes")
ax.set_ylim(0, 110)
ax.set_title("SpeedOfLight Throughput\n(k_force-Kernel, ncu)")
ax.legend(fontsize=8.5)
ax.grid(axis="y", alpha=0.3, zorder=0)

for xi, (c, m) in zip(x_pos, zip(comp, mem)):
    ax.text(xi - bar_w/2, c + 1.5, f"{c:.0f}%", ha="center", fontsize=9, fontweight="bold")
    ax.text(xi + bar_w/2, m + 1.5, f"{m:.0f}%", ha="center", fontsize=9, fontweight="bold")

ax.text(0, t4_compute + 5, "Compute-\nbegrenzt", ha="center", fontsize=8, color="#D85A30")
ax.text(1, a100_memory + 11, "Latenz-\nbegrenzt", ha="center", fontsize=8, color="#185FA5")

# ---- Plot 3: Occupancy -------------------------------------------------------
ax = axes[2]
occ_labels = list(occ_data.keys())
occ_vals   = list(occ_data.values())
cols = ["#185FA5", "#D85A30"]
ax.bar(range(len(occ_labels)), occ_vals, color=cols, width=0.5, zorder=3)
ax.axhline(100, ls="--", color="gray", lw=1.2, label="max. theoretisch (100%)")
for i, v in enumerate(occ_vals):
    ax.text(i, v + 2, f"{v:.0f}%", ha="center", fontsize=11, fontweight="bold")
ax.set_xticks(range(len(occ_labels)))
ax.set_xticklabels(occ_labels, fontsize=10)
ax.set_ylabel("Theoretische Okkupanz [%]")
ax.set_ylim(0, 120)
ax.set_title("GPU-Okkupanz\n(Registerlimitierung)")
ax.grid(axis="y", alpha=0.3, zorder=0)
ax.legend(fontsize=8)

# Erläuterungen -- mittig im Balken platziert (nicht auf Höhe der
# Prozent-Beschriftung über dem Balken, sonst überlappen sich beide Texte).
ax.text(0, occ_vals[0] / 2, "62 Reg/Thread\n→ 4 Blöcke/SM\n→ 32/64 Warps", ha="center",
        fontsize=7.5, color="white", va="center")
ax.text(1, occ_vals[1] / 2, "48 Reg/Thread\n→ 5 Blöcke/SM\n→ 40/64 Warps", ha="center",
        fontsize=7.5, color="white", va="center")

fig.suptitle("ncu-Profiling: k_force-Kernel (NVIDIA A100 / Tesla T4)", fontsize=13)
fig.tight_layout()
fig.savefig("ncu_profiling.png", dpi=150)
print("geschrieben: ncu_profiling.png")
