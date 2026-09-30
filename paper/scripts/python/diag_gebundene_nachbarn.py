#!/usr/bin/env python3
"""diag_gebundene_nachbarn.py -- Fliegen schnelle Atome durch gebundene Nachbarn?

Reproduziert den Befund in docs/befund_gebundene_nachbarn.md.

Aufbau: 21x21-Dreiecksgitter, PKA im Zentrum (Atom 220), keine Daempfung,
kein Healing. Verfolgt wird der Abstand des PKA zu seinen 6 gebundenen
Nachbarn aus den XYZ-Frames von cascade_serial.

  1. Schwellen frontal (90 Grad, entlang der Bindung zum rechten Nachbarn),
     per Bisektion: (a) PKA erreicht den Ort des Nachbarn (Abstand < 0,05),
     (b) PKA fliegt ganz hindurch und landet jenseits des Nachbarn.
  2. Winkelscan: Mindestabstand zu einem gebundenen Nachbarn ueber den
     Abschusswinkel. Vergleich mit dem Abstand einer geraden, unabgelenkten
     Linie, sin(Winkelabstand zur naechsten Bindung).
  3. Beide Modelle: rep_model = r12 (alt) und zbl (neu).

Aufruf:
  .venv-powerlaw/bin/python scripts/python/diag_gebundene_nachbarn.py \
      --out results/images/D1_gebundene_nachbarn.png
"""
import argparse
import os
import shutil
import subprocess
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(REPO, "src", "cascade_serial.c")
NX = 21
PKA = 10 * NX + 10
RIGHT = PKA + 1

CAT = ["#2a78d6", "#eb6834"]          # r12, zbl (validiert, siehe dataviz)
INK, INK2, SURFACE, GRID = "#1f1f1e", "#5c5c58", "#fcfcfb", "#ececE6"


def build(tmp):
    exe = os.path.join(tmp, "cascade_serial")
    subprocess.run(["gcc", "-O2", "-o", exe, SRC, "-lm"], check=True)
    return exe


def run(exe, tmp, model, energy, angle):
    """Ein Lauf; liefert Liste von (t, pos[Nx2]) und die Ruhelage."""
    tag = os.path.join(tmp, f"{model}_{energy:.1f}_{angle:.1f}")
    if model == "r12":
        dyn = "dt=2e-5\nn_steps=3000\n"
        dump = 5
    else:
        dyn = "dt=2e-5\ndt_adapt=true\nt_max=0.06\nn_steps=1000000\n"
        dump = 20
    with open(tag + ".ini", "w") as f:
        f.write(f"""[grid]
NX = {NX}
NY = {NX}
[potential]
rep_model={model}
[pka]
pka_x=0.5
pka_y=0.5
pka_energy={energy}
pka_angle={angle}
[dynamics]
{dyn}[output]
log_every=1000000
dump_every={dump}
out_prefix={tag}
[ensemble]
seed=1
healing=false
""")
    subprocess.run([exe, tag + ".ini"], stdout=subprocess.DEVNULL, check=True)
    frames = []
    with open(tag + ".xyz") as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            n = int(line)
            t = float(fh.readline().split("Time=")[1])
            pos = np.array([[float(v) for v in fh.readline().split()[1:3]]
                            for _ in range(n)])
            frames.append((t, pos))
    os.remove(tag + ".xyz")
    return frames


def bonded_neighbours(pos0):
    d = np.hypot(*(pos0 - pos0[PKA]).T)
    return [q for q in range(len(pos0)) if q != PKA and abs(d[q] - 1.0) < 1e-2]


def min_dist(frames, nbs):
    return min(np.hypot(*(pos[nbs] - pos[PKA]).T).min() for _, pos in frames)


def passes_through(frames):
    """Frontal nach rechts: PKA landet rechts vom Nachbarn."""
    return any(pos[PKA, 0] > pos[RIGHT, 0] + 0.05 for _, pos in frames)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(REPO, "results", "images",
                                                   "D1_gebundene_nachbarn.png"))
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix="bonddiag_")
    try:
        exe = build(tmp)

        # 1. Schwellen frontal (r12)
        def bisect(crit):
            lo, hi = 50.0, 400.0
            for _ in range(10):
                m = 0.5 * (lo + hi)
                if crit(run(exe, tmp, "r12", m, 90)):
                    hi = m
                else:
                    lo = m
            return lo, hi
        near = lambda fr: min(np.hypot(*(p[RIGHT] - p[PKA])) for _, p in fr) < 0.05
        lo, hi = bisect(near)
        print(f"r12, frontal: erreicht den Ort des Nachbarn ab E = {lo:.0f}-{hi:.0f}")
        lo, hi = bisect(passes_through)
        print(f"r12, frontal: fliegt ganz hindurch ab E = {lo:.0f}-{hi:.0f}")

        # Zeitverlauf frontal, E = 1000
        E_head = 1000
        traces = {}
        for model in ("r12", "zbl"):
            fr = run(exe, tmp, model, E_head, 90)
            t = np.array([f[0] for f in fr])
            dx = np.array([f[1][RIGHT, 0] - f[1][PKA, 0] for f in fr])
            dy = np.array([f[1][RIGHT, 1] - f[1][PKA, 1] for f in fr])
            traces[model] = (t, np.sign(dx) * np.hypot(dx, dy))
            print(f"{model}, frontal E={E_head}: min Abstand {np.abs(traces[model][1]).min():.3f}, "
                  f"Durchflug: {'ja' if passes_through(fr) else 'nein'}")

        # 2. Winkelscan, E = 1300 (kleinste Energie des Fixed-Energy-Ensembles)
        E_scan = 1300
        angles = np.arange(0, 61, 5)
        nbs = None
        scan = {"r12": [], "zbl": []}
        for model in scan:
            for ang in angles:
                fr = run(exe, tmp, model, E_scan, float(ang))
                if nbs is None:
                    nbs = bonded_neighbours(fr[0][1])
                scan[model].append(min_dist(fr, nbs))
        # Bindungen liegen bei 30 und 90 Grad (Konvention: 0 = -y, 90 = +x)
        geo = np.abs(np.sin(np.radians(angles - 30)))
        print(f"\nWinkelscan E={E_scan}: Winkel | gerade Linie | r12 | zbl")
        for ang, g, r, z in zip(angles, geo, scan["r12"], scan["zbl"]):
            print(f"  {ang:3d}°  {g:.3f}  {r:.3f}  {z:.3f}")
    finally:
        shutil.rmtree(tmp)

    # ---- Abbildung: zwei Panels, je eine eigene y-Achse -----------------
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2,
                         "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "text.color": INK})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.9), facecolor=SURFACE)
    for ax in (ax1, ax2):
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    lab = {"r12": "alt (r⁻¹², nur ungebunden)", "zbl": "neu (ZBL, alle Paare)"}
    for k, (model, (t, d)) in enumerate(traces.items()):
        ax1.plot(t, d, color=CAT[k], lw=2, label=lab[model])
    ax1.axhline(0, color=INK2, lw=1, ls="--")
    ax1.text(0.0005, 0.04, "gleicher Ort", color=INK2, fontsize=9, va="bottom")
    ax1.set_xlim(0, 0.03)
    ax1.set_xlabel("Zeit t")
    ax1.set_ylabel("Abstand PKA → Nachbar  [L0]\n(< 0: PKA ist durchgeflogen)")
    ax1.set_title(f"a) Frontal auf gebundenen Nachbarn, E = {E_head}",
                  loc="left", fontsize=10.5)
    ax1.legend(frameon=False, loc="upper right", fontsize=9)

    ax2.plot(angles, geo, color=INK2, lw=1.5, ls="--",
             label="gerade Linie ohne Ablenkung")
    for k, model in enumerate(("r12", "zbl")):
        ax2.plot(angles, scan[model], color=CAT[k], lw=2, marker="o", ms=5,
                 mec=SURFACE, mew=1.5, label=lab[model])
    ax2.set_xlabel("Abschusswinkel [°]   (Bindung bei 30°)")
    ax2.set_ylabel("kleinster Abstand zu einem\ngebundenen Nachbarn [L0]")
    ax2.set_ylim(0, 0.6)
    ax2.set_xticks(range(0, 61, 10))
    ax2.set_title(f"b) Alle Richtungen, E = {E_scan}", loc="left", fontsize=10.5)
    ax2.legend(frameon=False, loc="upper center", fontsize=9, ncol=1)

    fig.tight_layout()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=160, facecolor=SURFACE)
    print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()
