#!/usr/bin/env python3
"""
analyze_compiler_opt.py -- wertet den Compiler-Flag-Sweep aus.

Zwei Fragen gleichzeitig:
  1) WIE SCHNELL ist jede Flag-Kombination? (Laufzeit, Speedup gegen -O0)
  2) IST SIE NOCH KORREKT? (Energie-Drift -- aggressive Flags wie -ffast-math/
     -Ofast koennen die Energieerhaltung verletzen)

Gruppiert mehrere Wiederholungen je Label automatisch (Median der Zeit,
schlechtester |Drift| als sichere Obergrenze).

Aufruf:  python3 analyze_compiler_opt.py compiler_opt.csv [drift_tol_pct]
         drift_tol_pct = erlaubter |Drift| in % (Default 0.1)
"""
import sys, csv
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

FN  = sys.argv[1] if len(sys.argv) > 1 else "compiler_opt.csv"
TOL = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5

def fnum(x):
    try: return float(x)
    except: return float("nan")

# ---- einlesen + nach Label gruppieren -----------------------------------
groups = defaultdict(list)
with open(FN) as f:
    for r in csv.DictReader(f):
        groups[r["label"]].append(r)

rows = []
for label, rs in groups.items():
    flags = rs[0].get("flags", "")
    ok = all(r.get("compile_ok") == "1" for r in rs)
    if not ok:
        rows.append({"label": label, "flags": flags, "ok": False,
                     "wt": float("nan"), "drift": float("nan"), "ad": float("nan")})
        continue
    wts = [fnum(r.get("walltime_s")) for r in rs]
    ads = [fnum(r.get("abs_drift_pct")) for r in rs]
    drs = [fnum(r.get("drift_pct")) for r in rs]
    k = int(np.nanargmax(ads)) if np.any(~np.isnan(ads)) else 0   # schlechtester Drift
    rows.append({"label": label, "flags": flags, "ok": True,
                 "wt": float(np.nanmedian(wts)),
                 "ad": float(np.nanmax(ads)), "drift": drs[k]})

# ---- Baseline = -O0 (sonst langsamste lauffaehige) ----------------------
runnable = [d for d in rows if d["ok"] and d["wt"] == d["wt"]]
base = next((d for d in rows if d["label"] == "O0" and d["ok"]), None)
base_wt = base["wt"] if base else (max((d["wt"] for d in runnable), default=float("nan")))

# ---- Tabelle ------------------------------------------------------------
print(f"\nBaseline (Speedup-Bezug) = {'O0' if base else 'langsamste'}  "
      f"({base_wt:.3f}s) | Drift-Toleranz = {TOL}%\n")
print(f"{'Label':<18}{'Flags':<34}{'Zeit[s]':>9}{'Speedup':>9}{'Drift%':>10}  Energie")
print("-"*92)
order = sorted(rows, key=lambda d: (not d["ok"],
               d["wt"] if d["wt"] == d["wt"] else 1e18))
for d in order:
    if not d["ok"]:
        print(f"{d['label']:<18}{d['flags']:<34}{'--':>9}{'--':>9}{'--':>10}  KOMPILIERUNG FEHLT")
        continue
    sp = base_wt / d["wt"] if d["wt"] else float("nan")
    estat = "OK" if d["ad"] <= TOL else f"VERLETZT (|{d['drift']:.3f}|>{TOL})"
    print(f"{d['label']:<18}{d['flags']:<34}{d['wt']:>9.3f}{sp:>8.2f}x{d['drift']:>10.4f}  {estat}")

# ---- Empfehlung ---------------------------------------------------------
correct = [d for d in runnable if d["ad"] <= TOL]
fastest = min(runnable, key=lambda d: d["wt"]) if runnable else None
print()
if correct:
    best = min(correct, key=lambda d: d["wt"])
    print(f"-> Schnellstes KORREKTES Setup: {best['label']}  ({best['flags']})")
    print(f"   Zeit={best['wt']:.3f}s  Drift={best['drift']:.4f}%  "
          f"Speedup vs O0={base_wt/best['wt']:.2f}x")
if fastest and fastest["ad"] > TOL:
    print(f"-> ACHTUNG: das absolut schnellste Setup '{fastest['label']}' "
          f"VERLETZT die Energieerhaltung (Drift {fastest['drift']:.3f}%).")
    print(f"   Schneller, aber physikalisch falsch -> nicht fuer Produktion nutzen.")

# ---- Plot: Zeit-Balken, gruen=korrekt / rot=Energie verletzt ------------
rd = sorted(runnable, key=lambda d: d["wt"])
labels = [d["label"] for d in rd]; wts = [d["wt"] for d in rd]
cols = ["#2E8B57" if d["ad"] <= TOL else "#C0392B" for d in rd]
fig, ax = plt.subplots(figsize=(10, 5.2))
ax.bar(range(len(labels)), wts, color=cols)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=9)
ax.set_ylabel("Wandzeit [s]  (kleiner = schneller)")
ax.set_title(f"Compiler-Optimierungen: Zeit vs. Korrektheit\n"
             f"gruen = Energie ok, rot = Drift > {TOL}%")
for i, d in enumerate(rd):
    ax.text(i, wts[i], f"{wts[i]:.1f}", ha="center", va="bottom", fontsize=8)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#2E8B57", label="Energie ok"),
                   Patch(color="#C0392B", label=f"Drift > {TOL}%")], fontsize=9)
fig.tight_layout(); fig.savefig("compiler_opt.png", dpi=140)
print("\ngeschrieben: compiler_opt.png")
