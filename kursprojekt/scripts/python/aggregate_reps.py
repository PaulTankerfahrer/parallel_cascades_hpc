#!/usr/bin/env python3
"""
aggregate_reps.py -- fasst mehrere Wiederholungen pro Rankzahl zusammen.
Liest eine CSV mit mehreren Zeilen je Rankzahl (z. B. 3 Wiederholungen) und
schreibt eine CSV mit EINER Zeile je Rankzahl (Median, robust gegen Ausreisser).

Aufruf:  python3 aggregate_reps.py scaling_raw.csv scaling.csv [median|mean]
Danach:  python3 analyze_scaling.py scaling.csv
"""
import sys, numpy as np

fn_in  = sys.argv[1]
fn_out = sys.argv[2] if len(sys.argv) > 2 else "scaling_agg.csv"
stat   = sys.argv[3] if len(sys.argv) > 3 else "median"

d = np.genfromtxt(fn_in, delimiter=",", names=True)
names = d.dtype.names
ranks_all = np.atleast_1d(d["ranks"])

f = np.median if stat == "median" else np.mean
out_rows = []
print(f"Aggregation = {stat}")
print(f"{'ranks':>6} {'n_rep':>6} {'walltime':>10} {'min':>9} {'max':>9} {'Spanne%':>8}")
for r in sorted(np.unique(ranks_all)):
    mask = ranks_all == r
    sub = d[mask]
    row = [r]
    for nm in names[1:]:
        vals = np.atleast_1d(sub[nm])
        row.append(float(f(vals)))
    out_rows.append(row)
    w = np.atleast_1d(sub["walltime"])
    spread = 100*(w.max()-w.min())/np.median(w) if np.median(w) else 0
    print(f"{int(r):6d} {len(w):6d} {f(w):10.4f} {w.min():9.4f} {w.max():9.4f} {spread:8.1f}")

with open(fn_out, "w") as fout:
    fout.write(",".join(names) + "\n")
    for row in out_rows:
        fout.write(",".join((f"{int(v)}" if i == 0 else f"{v:.4f}")
                            for i, v in enumerate(row)) + "\n")
print("geschrieben:", fn_out)
