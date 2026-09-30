#!/usr/bin/env python3
"""
plot_gpu_saturation.py -- GPU-Sättigungskurve: Matom/s vs. N.

Zeigt, dass der A100 bereits ab N≈250k voll gesättigt ist (Plateau bei ~3800 Matom/s).
Daten aus cuda_bench_33375548.out (Sättigungs-Sweep, fp64, T=256).
Erzeugt: gpu_saturation.png
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Messwerte aus Sättigungs-Sweep (fp64, T=256, 400 Schritte)
data_sweep = [
    (250_000,    3834.5),
    (499_849,    4074.4),
    (1_000_000,  3690.5),
    (1_999_396,  3804.8),   # Hauptergebnis (3000 Schritte)
    (4_000_000,  3791.9),
    (7_997_584,  3870.6),
    (16_000_000, 3842.4),
]
N_arr    = np.array([d[0] for d in data_sweep], dtype=float)
mats_arr = np.array([d[1] for d in data_sweep])

fig, ax = plt.subplots(figsize=(9, 5.5))

ax.semilogx(N_arr, mats_arr, "o-", color="#185FA5", ms=8, lw=2, zorder=4, label="A100 fp64")

# Plateaulinie
plateau = np.mean(mats_arr)
ax.axhline(plateau, ls="--", color="#D85A30", lw=1.5,
           label=f"Plateau ≈ {plateau:.0f} Matom/s")

# Markierung: Arbeitspunkt 2 Mio
main_idx = 3
ax.plot(N_arr[main_idx], mats_arr[main_idx], "*", color="#D85A30", ms=18, zorder=5,
        label=f"Arbeitspunkt (N≈2M): {mats_arr[main_idx]:.0f} Matom/s")

# Beschriftung aller Punkte
for n, m in zip(N_arr, mats_arr):
    label = f"{int(n/1e6):.0f}M" if n >= 1e6 else f"{int(n/1e3):.0f}k"
    ax.annotate(f"{m:.0f}", (n, m), textcoords="offset points", xytext=(0, 8),
                ha="center", fontsize=8.5)

ax.set_xlabel("Problemgröße N  [Atome]")
ax.set_ylabel("Durchsatz  [Matom/s]")
ax.set_title("GPU-Sättigung: NVIDIA A100 (fp64, T=256)\nDurchsatz als Funktion der Systemgröße")
ax.set_ylim(0, 5000)
ax.legend(fontsize=10)
ax.grid(alpha=0.3, which="both")

# x-Ticks anpassen
ax.set_xticks(N_arr)
ax.set_xticklabels(["250k", "500k", "1M", "2M", "4M", "8M", "16M"], fontsize=9)

fig.tight_layout()
fig.savefig("gpu_saturation.png", dpi=150)
print(f"geschrieben: gpu_saturation.png")
print(f"Plateau-Mittelwert: {plateau:.1f} Matom/s")
print(f"Min: {mats_arr.min():.1f}, Max: {mats_arr.max():.1f}, Spanne: {(mats_arr.max()-mats_arr.min())/plateau*100:.1f}%")
