#!/usr/bin/env python3
"""
plot_energy.py -- Energiebilanz EINES Laufs plotten (Korrektheits-Check).
Funktioniert mit der Energie-CSV von cascade_serial UND mpi_cascade.
Nimmt den DATEINAMEN direkt (nicht den Prefix).

Aufruf:  python3 plot_energy.py np1_energy.csv
"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fn = sys.argv[1] if len(sys.argv) > 1 else "run_energy.csv"
e = np.genfromtxt(fn, delimiter=",", names=True)
t = e["t"]; E0 = e["E_total"][0]
drift = 100.0*(e["E_total"] - E0)/abs(E0)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

# links: Energiekomponenten + Summe
ax1.plot(t, e["E_kin"],    label="E_kin")
ax1.plot(t, e["E_spring"], label="E_Feder")
ax1.plot(t, e["E_rep"],    label="E_Abstoss")
ax1.plot(t, e["E_broken"], label="E_Bruch")
ax1.plot(t, e["E_total"],  "k--", lw=2, label="E_total")
ax1.set_xlabel("Zeit"); ax1.set_ylabel("Energie")
ax1.set_title("Energiebilanz (E_total flach = korrekt)")
ax1.legend(fontsize=8); ax1.grid(alpha=.3)

# rechts: Drift in Prozent
ax2.plot(t, drift, "o-", color="#D85A30")
ax2.axhline(0, color="gray", lw=.8)
ax2.set_xlabel("Zeit"); ax2.set_ylabel("Drift [%]")
ax2.set_title(f"Energie-Drift (final {drift[-1]:.3f} %)")
ax2.grid(alpha=.3)

fig.suptitle(f"Energie: {fn}", fontsize=11)
fig.tight_layout()
out = fn.replace(".csv", "") + "_plot.png"
fig.savefig(out, dpi=130)
print("geschrieben:", out)
print(f"final E_total={e['E_total'][-1]:.4f}  Drift={drift[-1]:.4f} %")
