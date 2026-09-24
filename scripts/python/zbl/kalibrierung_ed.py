#!/usr/bin/env python3
"""kalibrierung_ed.py -- Energieskala eps [eV/Einheit] ueber die Verlagerungsschwelle E_d.

Phase 2 von docs/plan_zbl.md.

Defekt:      nach t_max bleibt mindestens ein Leerplatz (Wigner-Seitz, defekte.py).
P(E):        Anteil der Richtungen, in denen bei Energie E ein Defekt bleibt.
             Richtungen 0..30 Grad in 1-Grad-Schritten; das Dreiecksgitter ist
             60-Grad-periodisch und bei 30 Grad gespiegelt, das deckt alle ab.
E_d:         mittlere Schwelle  E_d = INT_0^inf (1 - P(E)) dE  (Flaeche ueber P).
             Fuer eine Richtung mit scharfer Schwelle ist das genau die Schwelle;
             ueber Richtungen gemittelt ist es deren Mittelwert.  Anders als
             "erste Energie mit Defekt" zaehlen auch Energien oberhalb der ersten
             Schwelle, in denen der Defekt wieder rekombiniert (E_d ist in MD nicht
             streng monoton).  Das ist deutlich robuster: Eine Iteration mit
             "erster Defekt" pendelte um +-15 % (siehe docs/plan_zbl.md).
Kalibrierung: eps * E_d(eps) = 90 eV (ASTM, Wolfram).  Statt zu iterieren wird
             E_d bei mehreren eps gemessen und f(eps) = eps*E_d(eps) linear
             interpoliert.

Aufruf:
  .venv-powerlaw/bin/python scripts/python/zbl/kalibrierung_ed.py \
      --outdir results/zbl/kalibrierung --eps 0.35,0.40,0.45 --jobs 8
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from defekte import ws_analyse, lies_endzustand  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
SRC = os.path.join(REPO, "src", "cascade_serial.c")
ED_EV = 90.0
NX = 40
T_MAX = 3.0
ANGLES = np.arange(0.0, 30.5, 1.0)
E_GRID = 30.0 * 1.12 ** np.arange(0, 41)      # 30 ... ~2800

_CFG = {}


def _init(cfg):
    _CFG.update(cfg)


def one(arg):
    """Ein Lauf: (variante, energie, winkel) -> Anzahl Leerplaetze."""
    variant, block, energy, angle = arg
    scratch = tempfile.mkdtemp(prefix="ed_")
    tag = os.path.join(scratch, "run")
    try:
        with open(tag + ".ini", "w") as f:
            f.write(f"""[grid]
NX={NX}
NY={NX}
[potential]
{block}
[pka]
pka_energy={energy}
pka_angle={angle}
[dynamics]
dt=3e-4
dt_adapt=true
t_max={T_MAX}
n_steps=100000000
[output]
log_every=100000000
out_prefix={tag}
[ensemble]
healing=true
seed=1
""")
        subprocess.run([_CFG["bin"], tag + ".ini"], stdout=subprocess.DEVNULL, check=True)
        leer, _ = ws_analyse(lies_endzustand(tag + "_state_final.csv"), NX, NX)
    finally:
        shutil.rmtree(scratch)
    return variant, energy, angle, leer


def mean_threshold(E, P):
    """INT (1-P) dE; unterhalb des Gitters P = 0, oberhalb P = 1 angenommen."""
    trapz = getattr(np, "trapezoid", None) or np.trapz
    return float(E[0] + trapz(1.0 - P, E))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(REPO, "results", "zbl", "kalibrierung"))
    ap.add_argument("--eps", default="0.35,0.40,0.45")
    ap.add_argument("--keep-r12-at", type=float, default=None,
                    help="zusaetzlich zbl_keep_r12=true bei diesem eps")
    ap.add_argument("--r12", action="store_true", help="altes Modell als Referenz")
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    variants = {f"zbl eps={e}": f"rep_model=zbl\neps_eV={e}" for e in map(float, a.eps.split(","))}
    if a.keep_r12_at:
        variants[f"zbl+r12 eps={a.keep_r12_at}"] = \
            f"rep_model=zbl\neps_eV={a.keep_r12_at}\nzbl_keep_r12=true"
    if a.r12:
        variants["r12"] = "rep_model=r12"

    tmp = tempfile.mkdtemp(prefix="edbin_")
    exe = os.path.join(tmp, "cascade_serial")
    subprocess.run(["gcc", "-O3", "-march=native", "-o", exe, SRC, "-lm"], check=True)
    tasks = [(v, b, float(E), float(ang)) for v, b in variants.items()
             for E in E_GRID for ang in ANGLES]
    print(f"{len(tasks)} Laeufe, {len(variants)} Varianten", flush=True)
    res = []
    with Pool(a.jobs, initializer=_init, initargs=({"bin": exe},)) as pool:
        for k, r in enumerate(pool.imap_unordered(one, tasks, chunksize=4)):
            res.append(r)
            if (k + 1) % 250 == 0:
                print(f"  {k+1}/{len(tasks)}", flush=True)
    shutil.rmtree(tmp)

    with open(os.path.join(a.outdir, "ed_laeufe.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variante", "energie", "winkel", "leerplaetze"])
        w.writerows(sorted(res))

    # Auswertung
    summary = {}
    for v in variants:
        P = np.array([np.mean([r[3] > 0 for r in res if r[0] == v and r[1] == E])
                      for E in E_GRID])
        first = []
        for ang in ANGLES:
            hit = sorted(r[1] for r in res if r[0] == v and r[2] == ang and r[3] > 0)
            first.append(hit[0] if hit else np.nan)
        summary[v] = dict(P=P, ed=mean_threshold(E_GRID, P),
                          e50=float(np.interp(0.5, np.maximum.accumulate(P), E_GRID)),
                          first=np.array(first))
        print(f"{v:22s} E_d = INT(1-P) = {summary[v]['ed']:7.1f}  "
              f"P=50% bei {summary[v]['e50']:7.1f}  "
              f"erste Schwelle min {np.nanmin(first):.0f} max {np.nanmax(first):.0f}", flush=True)

    for v, s in summary.items():
        if s["P"][-1] < 1.0:
            print(f"WARNUNG: {v}: P(E_max) = {s['P'][-1]:.2f} < 1, E_d unterschaetzt", flush=True)
    zbl = sorted((float(v.split("=")[1]), s["ed"]) for v, s in summary.items()
                 if v.startswith("zbl eps"))
    eps_arr = np.array([e for e, _ in zbl])
    f_arr = np.array([e * ed for e, ed in zbl])
    if np.all(np.diff(f_arr) > 0) and f_arr[0] <= ED_EV <= f_arr[-1]:
        eps_cal = float(np.interp(ED_EV, f_arr, eps_arr))
    else:
        eps_cal = float("nan")
        print(f"WARNUNG: {ED_EV} eV nicht von eps*E_d eingeschlossen oder nicht monoton, "
              f"eps-Bereich erweitern", flush=True)
    print(f"f(eps) = eps*E_d: " + ", ".join(f"{e}: {fv:.1f} eV" for e, fv in zip(eps_arr, f_arr)))
    print(f"-> eps = {eps_cal:.4f} eV/Einheit  (E_d = {ED_EV/eps_cal:.1f} Einheiten)")

    with open(os.path.join(a.outdir, "ergebnis.md"), "w") as f:
        f.write("# Kalibrierung der Energieskala über E_d\n\n")
        f.write(f"Gitter {NX}×{NX}, t_max = {T_MAX}, Healing an, keine Bremsung. "
                f"{len(ANGLES)} Richtungen (0–30°, 1°), {len(E_GRID)} Energien "
                f"({E_GRID[0]:.0f}–{E_GRID[-1]:.0f}, Faktor 1,12). "
                f"Defekt = mindestens ein Leerplatz nach t_max (Wigner-Seitz).\n\n")
        f.write("E_d = ∫(1 − P(E)) dE, P = Anteil der Richtungen mit Defekt.\n\n")
        f.write("| Variante | E_d [Einh.] | P = 50 % bei [Einh.] | erste Schwelle min–max [Einh.] | ε·E_d [eV] |\n")
        f.write("|---|--:|--:|--:|--:|\n")
        for v, s in summary.items():
            ev = f"{float(v.split('=')[1])*s['ed']:.1f}" if "eps=" in v else "–"
            f.write(f"| {v} | {s['ed']:.1f} | {s['e50']:.1f} | "
                    f"{np.nanmin(s['first']):.0f}–{np.nanmax(s['first']):.0f} | {ev} |\n")
        f.write(f"\n**Ergebnis: ε = {eps_cal:.3f} eV pro Energieeinheit** "
                f"(lineare Interpolation von ε·E_d(ε) auf {ED_EV:.0f} eV).\n\n")
        f.write("## P(E)\n\n| E [Einh.] | " + " | ".join(summary) + " |\n|--:|" +
                "--:|" * len(summary) + "\n")
        for k, E in enumerate(E_GRID):
            f.write(f"| {E:.0f} | " + " | ".join(f"{summary[v]['P'][k]:.2f}" for v in summary) + " |\n")
    np.savez(os.path.join(a.outdir, "p_von_e.npz"), E=E_GRID,
             **{v.replace(" ", "_").replace("=", ""): s["P"] for v, s in summary.items()})
    print(f"-> {a.outdir}/ergebnis.md")


if __name__ == "__main__":
    main()
