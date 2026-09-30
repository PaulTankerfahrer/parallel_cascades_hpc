#!/usr/bin/env python3
"""fig_force_law.py -- Die Kraft zwischen zwei Nachbaratomen, eine Kurve.

Die alte Abbildung M2_force_law.png zeigte Kraefte und Potentiale in zwei
Panels mit sich ueberlappender Legende. Fuer den Vortrag genuegt eine einzige
Kurve: die Gesamtkraft zwischen zwei benachbarten Atomen als Funktion ihres
Abstands. Daran laesst sich alles ablesen, was das Modell ausmacht.

Damit die Kurve ohne Erklaerung lesbar ist:
  * die y-Achse ist mit "anziehend"/"abstoszend" beschriftet statt mit Zahlen,
  * die drei Bereiche (zu nah / Feder / gerissen) sind hinterlegt und benannt,
  * der Sprung beim Reissen ist als gestrichelte Linie eingezeichnet.

Parameter wie in src/params.ini:
  L0 = 1.0, K_SPRING = 100, MAX_STRETCH = 1.15, RCUT = 0.9, K_REP = 400, n = 12

Aufruf:
  python3 scripts/python/fig_force_law.py presentation/figures/P2_kraft.png
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L0, K, MAX_STRETCH = 1.0, 100.0, 1.15
RCUT, KREP, NREP = 0.9, 400.0, 12.0

CAT = ["#2a78d6", "#eb6834"]
INK, INK2, SURFACE = "#1f1f1e", "#5c5c58", "#fcfcfb"
GRID = "#ececE6"


def kraft(r):
    """Positiv = die Atome stossen sich ab, negativ = sie ziehen sich an."""
    f = -K * (r - L0)                      # Feder
    f = np.where(r > MAX_STRETCH * L0, 0.0, f)   # gerissen: keine Feder mehr
    f = f + np.where(r < RCUT, KREP * (RCUT / r) ** NREP, 0.0)
    return f


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "presentation/figures/P2_kraft.png"
    fig, ax = plt.subplots(figsize=(5.6, 3.9), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    xlo, xhi, ylo, yhi = 0.855, 1.33, -36, 54

    # Die drei Bereiche hinterlegen
    ax.axvspan(xlo, RCUT, color=CAT[1], alpha=0.09, lw=0)
    ax.axvspan(MAX_STRETCH, xhi, color=INK2, alpha=0.07, lw=0)

    # Kurve in drei Stuecken, damit der Sprung beim Reissen sichtbar bleibt
    r1 = np.linspace(0.80, MAX_STRETCH, 600)
    ax.plot(r1, kraft(r1), color=CAT[0], lw=3.2, solid_capstyle="round")
    r2 = np.linspace(MAX_STRETCH, xhi, 60)
    ax.plot(r2, np.zeros_like(r2), color=CAT[0], lw=3.2, alpha=0.30)
    f_riss = kraft(np.array([MAX_STRETCH]))[0]
    ax.plot([MAX_STRETCH, MAX_STRETCH], [f_riss, 0], ls=(0, (2, 2)),
            color=CAT[1], lw=1.8, zorder=4)
    ax.plot([MAX_STRETCH], [f_riss], "o", ms=9, mfc=SURFACE, mec=CAT[1],
            mew=2.6, zorder=5)

    ax.axhline(0, color=INK2, lw=1.0)
    ax.axvline(L0, color=INK2, lw=1.0, ls=":")

    ax.annotate("zu nah:\nstoßen sich stark ab", xy=(0.905, 40),
                xytext=(0.965, 44), fontsize=13, color=CAT[0],
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=CAT[0], lw=1.1))
    ax.annotate("Feder", xy=(1.045, -2.5), fontsize=14, color=CAT[0],
                ha="left", va="bottom", rotation=-30, rotation_mode="anchor")
    ax.annotate("Bindung gerissen:\nkeine Kraft mehr", xy=(1.238, 24),
                fontsize=13, color=CAT[1], ha="center", va="center")

    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.set_xlabel("Abstand zweier Nachbaratome, in Gitterweiten",
                  fontsize=14, color=INK2)
    ax.set_ylabel("Kraft", fontsize=14, color=INK2)
    ax.set_yticks([-22, 0, 38])
    ax.set_yticklabels(["anziehend", "0", "abstoßend"])
    ax.set_xticks([0.9, 1.0, MAX_STRETCH, 1.3])
    ax.set_xticklabels(["0,9", "1,0\nRuhelage", "1,15\nreißt", "1,3"])
    ax.tick_params(colors=INK2, labelsize=12.5)
    ax.grid(True, axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#dcdcd6")

    fig.tight_layout()
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    print("geschrieben:", out)


if __name__ == "__main__":
    main()
