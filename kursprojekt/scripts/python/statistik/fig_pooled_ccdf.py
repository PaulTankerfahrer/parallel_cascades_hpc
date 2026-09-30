#!/usr/bin/env python3
"""fig_pooled_ccdf.py -- CCDF des urspruenglichen 1000-Kaskaden-Ensembles.

Fuer die Folie "Projektergebnis": die gepoolten Clustergroessen mit und ohne
Defektheilung, dazu die angepasste Gerade. Grosse Schrift und deutsche
Beschriftung, weil die Abbildung auf der Folie stark verkleinert wird.

Aufruf:
  .venv-powerlaw/bin/python scripts/python/statistik/fig_pooled_ccdf.py \
      results/statistik/ensemble_gepoolt presentation/figures/P5_ccdf_gepoolt.png
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clauset_refit import load, fit_of, pl_stats  # noqa: E402

CAT = ["#2a78d6", "#eb6834"]
INK, INK2, SURFACE = "#1f1f1e", "#5c5c58", "#fcfcfb"
SETS = [("heal", "mit Defektheilung"), ("no_heal", "ohne Defektheilung")]


def ccdf(d):
    x = np.sort(np.unique(d))
    return x, np.array([(d >= v).sum() / len(d) for v in x])


def main():
    ensdir, out = sys.argv[1], sys.argv[2]
    fig, ax = plt.subplots(figsize=(5.2, 4.4), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for i, (sub, label) in enumerate(SETS):
        d = load(os.path.join(ensdir, sub, "cluster_sizes.csv"))
        alpha, xmin, sigma, D, ntail = pl_stats(fit_of(d))
        x, y = ccdf(d)
        ax.plot(x, y, "o", ms=4.4, color=CAT[i], mew=0,
                label=f"{label}:  $S={alpha:.2f}$".replace(".", ","))
        xf = np.logspace(np.log10(max(xmin, 1)), np.log10(d.max()), 60)
        y0 = (d >= max(xmin, 1)).sum() / len(d)
        ax.plot(xf, y0 * (xf / max(xmin, 1)) ** (1 - alpha), "-",
                color=CAT[i], lw=1.8, alpha=0.55)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Clustergröße $n$", color=INK2, fontsize=15)
    ax.set_ylabel("Anteil Cluster mit\nmindestens $n$ Atomen",
                  color=INK2, fontsize=15)
    ax.tick_params(colors=INK2, labelsize=13)
    ax.grid(True, which="major", color="#ececE6", lw=0.7)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_color("#dcdcd6")
    leg = ax.legend(frameon=True, fontsize=12.5, loc="lower left")
    leg.get_frame().set_edgecolor("#dcdcd6")
    leg.get_frame().set_facecolor(SURFACE)
    for t in leg.get_texts():
        t.set_color(INK)
    fig.tight_layout()
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print("geschrieben:", out)


if __name__ == "__main__":
    main()
