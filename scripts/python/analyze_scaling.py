#!/usr/bin/env python3
"""
analyze_scaling.py -- Strong-Scaling-Auswertung mit Effizienz + 2-Knoten-Vergleich.

Erzeugt vier Grafiken (2x2):
  1) Speedup            (ideal + 1 Knoten + optional 2 Knoten)
  2) Parallele Effizienz(100%-Referenz + 1 Knoten + optional 2 Knoten)
  3) Phasen-Breakdown   (wo geht die Zeit hin; 1 Knoten + 2-Knoten-Vergleichsbalken)
  4) Kommunikationsanteil vs. Kerne (erklaert, WO das Scaling einbricht)

CSV-Spalten:  ranks,walltime,force,cells,hash,halo,migrate,energy,comm
Aufruf:
  python3 analyze_scaling.py mpi_strong_scaling_1node.csv [mpi_strong_scaling_2node.csv]

Die 2-Knoten-Speedups/Effizienzen werden gegen DIESELBE Baseline T(1) (1 Kern,
1 Knoten) berechnet -> 1- und 2-Knoten-Punkte sind direkt vergleichbar.
"""
import sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

f1 = sys.argv[1] if len(sys.argv) > 1 else "mpi_strong_scaling_1node.csv"
f2 = sys.argv[2] if len(sys.argv) > 2 else None

PHASES  = ["force", "cells", "hash", "halo", "migrate", "energy"]
PCOLORS = ["#185FA5", "#3E8AC0", "#7FB2D8", "#D85A30", "#E8915F", "#9AA7B2"]

def load(fn):
    d = np.genfromtxt(fn, delimiter=",", names=True)
    r = np.atleast_1d(d["ranks"]).astype(float)
    o = np.argsort(r)
    out = {"ranks": r[o], "walltime": np.atleast_1d(d["walltime"])[o],
           "comm": np.atleast_1d(d["comm"])[o]}
    for p in PHASES:
        out[p] = np.atleast_1d(d[p])[o]
    return out

A = load(f1)
B = load(f2) if f2 else None

T1   = A["walltime"][0]                 # Referenz: kleinste Rankzahl 1 Knoten (=1 Kern)
r1   = A["ranks"]; w1 = A["walltime"]
sp1  = T1 / w1
eff1 = sp1 / (r1 / r1[0])
if B:
    r2, w2 = B["ranks"], B["walltime"]
    sp2  = T1 / w2                       # SELBE Baseline
    eff2 = sp2 / (r2 / r1[0])

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

# ---- (1) Speedup --------------------------------------------------------
ax1.plot(r1, r1 / r1[0], "--", color="gray", label="ideal (linear)")
ax1.plot(r1, sp1, "o-", color="#185FA5", label="1 Knoten")
if B:
    ax1.plot(r2, sp2, "s--", color="#D85A30", ms=9, label="2 Knoten")
ax1.set_xscale("log", base=2); ax1.set_yscale("log", base=2)
ax1.set_xlabel("Kerne"); ax1.set_ylabel("Speedup  T(1)/T(p)")
ax1.set_title("1) Strong-Scaling-Speedup")
ax1.legend(fontsize=9); ax1.grid(alpha=.3, which="both")

# ---- (2) Effizienz (eigene Kurve, wie in der Anforderung) ---------------
ax2.axhline(100, ls="--", color="gray", label="ideal (100%)")
ax2.plot(r1, eff1 * 100, "o-", color="#185FA5", label="1 Knoten")
if B:
    ax2.plot(r2, eff2 * 100, "s--", color="#D85A30", ms=9, label="2 Knoten")
for x, e in zip(r1, eff1):
    ax2.annotate(f"{e*100:.0f}%", (x, e*100), fontsize=7, color="#185FA5",
                 xytext=(0, 6), textcoords="offset points", ha="center")
ax2.set_xscale("log", base=2)
ax2.set_xlabel("Kerne"); ax2.set_ylabel("parallele Effizienz [%]")
ax2.set_title("2) Parallele Effizienz")
ax2.legend(fontsize=9); ax2.grid(alpha=.3, which="both")

# ---- (3) Phasen-Zusammensetzung (normiert -> bei hohen Ranks sichtbar) --
#  Absolute Zeiten werden vom 1-Kern-Lauf erdrueckt; normiert sieht man bei
#  JEDER Rankzahl, welcher Teil dominiert (Flaschenhals-Verschiebung).
def comp(D):
    tot = np.sum([D[p] for p in PHASES], axis=0)
    return {p: 100.0 * D[p] / tot for p in PHASES}
xp1 = list(range(len(r1)))
labels = [f"{int(x)}" for x in r1]
C1 = comp(A)
bottom = np.zeros(len(r1))
for p, c in zip(PHASES, PCOLORS):
    ax3.bar(xp1, C1[p], bottom=bottom, color=c, label=p, width=0.82)
    bottom += C1[p]
xp2 = []
if B:
    gap = 1
    xp2 = [len(r1) + gap + i for i in range(len(r2))]
    C2 = comp(B)
    bottom2 = np.zeros(len(r2))
    for p, c in zip(PHASES, PCOLORS):
        ax3.bar(xp2, C2[p], bottom=bottom2, color=c, width=0.82)
        bottom2 += C2[p]
    labels += [f"{int(x)}\n(2N)" for x in r2]
    ax3.axvline(len(r1) + gap/2 - 0.5, ls=":", color="gray", lw=1)
ax3.set_xticks(xp1 + xp2); ax3.set_xticklabels(labels, fontsize=8)
ax3.set_xlabel("Kerne   (rechts der Linie: 2 Knoten)")
ax3.set_ylabel("Phasen-Anteil [%]")
ax3.set_title("3) Phasen-Zusammensetzung (wo dominiert was?)")
ax3.set_ylim(0, 100)
ax3.legend(fontsize=8, ncol=2)

# ---- (4) Kommunikationsanteil ------------------------------------------
ax4.axhline(50, ls=":", color="gray")
ax4.plot(r1, 100 * A["comm"] / w1, "s-", color="#D85A30", label="1 Knoten")
if B:
    ax4.plot(r2, 100 * B["comm"] / w2, "D--", color="#7A3B8A", ms=9, label="2 Knoten")
ax4.set_xscale("log", base=2)
ax4.set_xlabel("Kerne"); ax4.set_ylabel("Kommunikation / Wandzeit [%]")
ax4.set_title("4) Kommunikationsanteil (wo bricht es ein?)")
ax4.legend(fontsize=9); ax4.grid(alpha=.3, which="both")

fig.suptitle("Strong Scaling der MPI-Version", fontsize=14, y=0.995)
fig.tight_layout()
fig.savefig("strong_scaling.png", dpi=140)

# ---- Konsolen-Zusammenfassung ------------------------------------------
def at(ranks, vals, k):
    idx = np.where(ranks == k)[0]
    return float(vals[idx[0]]) if len(idx) else None
print(f"Baseline T(1) = {T1:.3f} s (1 Kern, 1 Knoten)")
print(f"1 Knoten  @128: Speedup {at(r1,sp1,128):.1f}x  Effizienz {at(r1,eff1,128)*100:.0f}%"
      f"  comm {at(r1,100*A['comm']/w1,128):.0f}%")
if B:
    s128_2 = at(r2, sp2, 128)
    if s128_2:
        print(f"2 Knoten  @128: Speedup {s128_2:.1f}x  Effizienz {at(r2,eff2,128)*100:.0f}%"
              f"  -> 1N/2N-Verhaeltnis {at(r1,w1,128)/at(r2,w2,128):.3f} (≈1 = Inter-Node vernachlaessigbar)")
    s256 = at(r2, sp2, 256)
    if s256:
        print(f"2 Knoten  @256: Speedup {s256:.1f}x  Effizienz {at(r2,eff2,256)*100:.0f}%")
print("geschrieben: strong_scaling.png")