#!/usr/bin/env python3
"""
fig_fixed_energy_ccdf.py -- CCDF-Vergleich des Fixed-Energy-Ensembles.

Links  die bisherige Clusterdefinition (disp > 0.5, Linkradius 1.5),
rechts die feinere (disp > 1.0, Linkradius 1.0) -- jeweils dieselben
900 Laeufe, dieselben Atome, nur anders zu Clustern zusammengefasst.

Aufruf:
  .venv-powerlaw/bin/python scripts/python/statistik/fig_fixed_energy_ccdf.py \
      results/statistik/ensemble_fest results/statistik/fit_fest/fixed_energy_ccdf.png

Optionen:
  --panel=alt|fein   nur eine der beiden Definitionen (Standard: beide)
  --folie            groessere Schrift, kein Obertitel -- fuer Beamer-Folien,
                     auf denen die Abbildung stark verkleinert wird
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import powerlaw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clauset_refit import load, fit_of, pl_stats  # noqa: E402

CAT = ["#2a78d6", "#eb6834", "#1baf7a"]        # kategoriale Slots 1-3
INK, INK2, SURFACE = "#1f1f1e", "#5c5c58", "#fcfcfb"
ENERGIES = [1300, 1800, 2400]
PANELS = [("clusters_thr050_link150",
           "alte Definition: verschoben ab 0,5, zusammen unter 1,5"),
          ("clusters_thr100_link100",
           "neue Definition: verschoben ab 1,0, zusammen unter 1,0")]


def ccdf(d):
    x = np.sort(np.unique(d))
    y = np.array([(d >= v).sum() / len(d) for v in x])
    return x, y


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = [a for a in sys.argv[1:] if a.startswith("--")]
    ensdir, out = args[0], args[1]
    folie = "--folie" in opts
    which = next((o.split("=", 1)[1] for o in opts if o.startswith("--panel=")), "both")
    panels = {"alt": PANELS[:1], "fein": PANELS[1:], "both": PANELS}[which]
    if "--folie" in opts:   # auf der Folie steht die Definition im Text daneben
        panels = [(tag, t.split(":")[0]) for tag, t in panels]

    # Folienvariante: groessere Schrift, weil die Abbildung auf der Folie
    # auf etwa die Haelfte verkleinert wird.
    F = dict(title=16, label=15, tick=13, leg=13.5) if folie \
        else dict(title=12, label=11, tick=9, leg=9)
    w = 5.6 * len(panels) if folie else 5.5 * len(panels)
    h = 5.0 if folie else 4.6

    fig, axes = plt.subplots(1, len(panels), figsize=(w, h), facecolor=SURFACE,
                             squeeze=False)
    axes = axes[0]
    for ax, (tag, title) in zip(axes, panels):
        ax.set_facecolor(SURFACE)
        for i, E in enumerate(ENERGIES):
            d = load(os.path.join(ensdir, f"E{E:04d}", tag + ".csv"))
            fit = fit_of(d)
            alpha, xmin, sigma, D, ntail = pl_stats(fit)
            x, y = ccdf(d)
            ax.plot(x, y, "o", ms=4.6 if folie else 3.4, color=CAT[i], mew=0,
                    label=(f"$E={E}$:  $S={alpha:.2f}$,  Abw. ${D:.3f}$"
                           .replace(".", ",") if folie
                           else f"$E={E}$:  $S={alpha:.2f}$,  $D={D:.3f}$"))
            # angepasstes Potenzgesetz ab x_min, auf die Daten dort normiert
            xf = np.logspace(np.log10(max(xmin, 1)), np.log10(d.max()), 60)
            y0 = (d >= max(xmin, 1)).sum() / len(d)
            ax.plot(xf, y0 * (xf / max(xmin, 1)) ** (1 - alpha),
                    "-", color=CAT[i], lw=1.6, alpha=0.55)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Clustergröße $n$", color=INK2, fontsize=F["label"])
        ax.set_title(title, color=INK, fontsize=F["title"], pad=9)
        ax.tick_params(colors=INK2, labelsize=F["tick"])
        ax.grid(True, which="major", color="#ececE6", lw=0.7)
        ax.set_axisbelow(True)
        for s in ax.spines.values():
            s.set_color("#dcdcd6")
        leg = ax.legend(frameon=True, fontsize=F["leg"], loc="lower left")
        leg.get_frame().set_edgecolor("#dcdcd6")
        leg.get_frame().set_facecolor(SURFACE)
        for t in leg.get_texts():
            t.set_color(INK)
    axes[0].set_ylabel("Anteil Cluster mit\nmindestens $n$ Atomen" if folie
                       else "$P(N \\geq n)$",
                       color=INK2, fontsize=F["label"])
    if folie:
        fig.tight_layout()
    else:
        fig.suptitle("Feste PKA-Energie: dieselben 900 Läufe, zwei Clusterdefinitionen",
                     color=INK, fontsize=13)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print("geschrieben:", out)


if __name__ == "__main__":
    main()
