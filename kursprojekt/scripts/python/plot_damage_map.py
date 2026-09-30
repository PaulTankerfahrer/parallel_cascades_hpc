#!/usr/bin/env python3
"""
plot_damage_map.py -- Schadenskarte einer Einzelkaskade, mit vs. ohne Healing.

Liest die per-Atom-Endzustaende (x, y, disp, broken) aus
results/ensemble_run_0999_state_final.csv (mit Healing) und
results/ensemble_no_heal_run_0999_state_final.csv (ohne Healing) -- siehe
Abschnitt "Machbarkeitsstudie: Skalierung des Ensembles" im Bericht.

Die CSVs enthalten keine explizite Bindungsliste (keine *_broken_bonds.csv),
nur einen kumulativen Bruchzaehler pro Atom. Um trotzdem Bindungen als
Liniensegmente zeichnen zu koennen, wird das Dreiecksgitter geometrisch
rekonstruiert:
  1) Aus (x, y) werden die ganzzahligen Gitterindizes (i, j) zurueckgerechnet
     (Zeilenabstand sqrt(3)/2 * L0, Zeilen-Versatz 0.5 * L0 bei ungeradem j).
  2) Für jedes Atom werden die 6 nominellen Gitternachbarn ueber die
     bekannte Dreiecksgitter-Konnektivitaet bestimmt.
  3) Eine Bindung gilt als gerissen, wenn der AKTUELLE Abstand der beiden
     Atome die im Bericht dokumentierte Bruchschwelle r > 1.15 L0
     ueberschreitet -- dieselbe Regel wie im Simulationscode
     (vgl. Abschnitt "Physikalisches Modell").

Damit ist jede gezeichnete Bindung direkt aus den Rohdaten ableitbar, ohne
Bindungsidentitaeten zu erfinden.

Aufruf:
  python3 plot_damage_map.py [no_heal_csv] [heal_csv] [--zoom R]
Default-Pfade passen zum Projektlayout (results/ zwei Verzeichnisebenen
oberhalb von scripts/python/ oder als Pfad-Argument uebergeben).
Ausgabe: vid_S2_damage.png
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L0 = 1.0
ROW_DY = np.sqrt(3) / 2 * L0
BREAK_R = 1.15 * L0

# Nominelle Nachbar-Offsets (di, dj) im Dreiecksgitter, getrennt nach
# Zeilenparitaet (siehe Docstring / Herleitung ueber die Zeilen-Versatz-Formel).
NBR_EVEN_J = [(-1, 0), (1, 0), (0, 1), (-1, 1), (0, -1), (-1, -1)]
NBR_ODD_J  = [(-1, 0), (1, 0), (0, 1), (1, 1), (0, -1), (1, -1)]


def load(fn):
    d = np.genfromtxt(fn, delimiter=",", names=True)
    return d


def reconstruct_bonds(d, cx, cy, zoom):
    """Baut (i,j)->Position-Lookup und gibt intakte/gerissene Segmente
    innerhalb des Zoom-Fensters um (cx, cy) zurueck."""
    x, y, disp = d["x"], d["y"], d["disp"]
    j = np.round(y / ROW_DY).astype(int)
    i = np.round(x - 0.5 * (j % 2)).astype(int)

    # Nur Atome im (grosszuegigen) Zoom-Fenster betrachten -- spart Zeit bei
    # 62500 Atomen und reicht, um auch Randbindungen des Fensters zu zeichnen.
    pad = 2.0
    sel = (np.abs(x - cx) <= zoom + pad) & (np.abs(y - cy) <= zoom + pad)
    idx_sel = np.nonzero(sel)[0]

    pos = {}       # (i,j) -> (x, y)
    for k in idx_sel:
        pos[(i[k], j[k])] = (x[k], y[k])

    intact, broken = [], []
    seen = set()
    for k in idx_sel:
        ii, jj = i[k], j[k]
        offsets = NBR_EVEN_J if jj % 2 == 0 else NBR_ODD_J
        for di, dj in offsets:
            key = (ii + di, jj + dj)
            bond_key = frozenset({(ii, jj), key})
            if bond_key in seen:
                continue
            seen.add(bond_key)
            if key not in pos:
                continue
            x0, y0 = pos[(ii, jj)]
            x1, y1 = pos[key]
            r = np.hypot(x1 - x0, y1 - y0)
            seg = ((x0, y0), (x1, y1))
            if r > BREAK_R:
                broken.append(seg)
            else:
                intact.append(seg)
    return intact, broken


def panel(ax, d, title, cx, cy, zoom, vmin, vmax, cmap):
    intact, broken = reconstruct_bonds(d, cx, cy, zoom)

    from matplotlib.collections import LineCollection
    if intact:
        ax.add_collection(LineCollection(intact, colors="#b0b0b0", linewidths=0.4,
                                          alpha=0.6, zorder=1))
    if broken:
        ax.add_collection(LineCollection(broken, colors="#d62728", linewidths=0.9,
                                          alpha=0.85, zorder=2))

    x, y, disp = d["x"], d["y"], d["disp"]
    sel = (np.abs(x - cx) <= zoom) & (np.abs(y - cy) <= zoom)
    logdisp = np.log10(np.clip(disp[sel], 1e-3, None))
    sc = ax.scatter(x[sel], y[sel], c=logdisp, cmap=cmap, vmin=vmin, vmax=vmax,
                     s=10, zorder=3, edgecolors="none")
    ax.scatter([cx], [cy], marker="x", s=140, c="black", linewidths=2.5, zorder=4)

    ax.set_xlim(cx - zoom, cx + zoom)
    ax.set_ylim(cy - zoom, cy + zoom)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_title(title)
    ax.grid(alpha=0.25)
    return sc, len(broken)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("no_heal_csv", nargs="?",
                     default="results/ensemble_no_heal_run_0999_state_final.csv")
    ap.add_argument("heal_csv", nargs="?",
                     default="results/ensemble_run_0999_state_final.csv")
    ap.add_argument("--zoom", type=float, default=15.0,
                     help="Zoom-Radius um den Einschlag in L0 (Default: 15)")
    args = ap.parse_args()

    d_no  = load(args.no_heal_csv)
    d_yes = load(args.heal_csv)

    # Einschlagzentrum: Position des am staerksten verschobenen Atoms. Ein
    # disp^2-gewichteter Schwerpunkt ueber die GESAMTE Domaene wurde verworfen,
    # da schon schwaches Hintergrundrauschen weit ausserhalb der Kaskade
    # (viele Atome mit kleinem, aber nicht-null disp) das Zentrum spuerbar
    # verzieht -- das Maximum liegt dagegen robust im Kaskadenkern.
    def centroid(d):
        k = np.argmax(d["disp"])
        return d["x"][k], d["y"][k]

    cx_no, cy_no = centroid(d_no)
    cx_yes, cy_yes = centroid(d_yes)
    # gemeinsames Zentrum (beide Laeufe starten am selben PKA-Ort)
    cx, cy = 0.5 * (cx_no + cx_yes), 0.5 * (cy_no + cy_yes)

    vmax = max(np.log10(max(d_no["disp"].max(), 1e-3)),
               np.log10(max(d_yes["disp"].max(), 1e-3)))
    vmin = -3.0
    cmap = "inferno"

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4), sharey=True)
    sc, n_broken_no = panel(axes[0], d_no, "(a) ohne Healing", cx, cy, args.zoom,
                             vmin, vmax, cmap)
    _, n_broken_yes = panel(axes[1], d_yes, "(b) mit Healing", cx, cy, args.zoom,
                             vmin, vmax, cmap)
    axes[0].set_ylabel("y")

    fig.suptitle(f"S2 -- Schadenskarte  ({n_broken_no} vs. {n_broken_yes} "
                 f"gerissene Bindungen im Zoom-Fenster)", fontsize=13)

    # Farbbalken ausserhalb beider Achsen -- verhindert die Ueberlappung mit
    # den Streudaten, die im urspruenglichen Plot auftrat.
    fig.tight_layout(rect=[0, 0, 0.9, 0.95])
    cax = fig.add_axes([0.915, 0.12, 0.02, 0.72])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label(r"Verschiebung $|u|$ (log)")
    ticks = np.arange(np.ceil(vmin), np.floor(vmax) + 1)
    cb.set_ticks(ticks)
    cb.set_ticklabels([f"$10^{{{int(t)}}}$" for t in ticks])

    fig.savefig("vid_S2_damage.png", dpi=150)
    print("geschrieben: vid_S2_damage.png")
    print(f"Einschlagzentrum (gemittelt): ({cx:.2f}, {cy:.2f})")
    print(f"gerissene Bindungen im Zoom-Fenster: ohne Healing={n_broken_no}, "
          f"mit Healing={n_broken_yes}")


if __name__ == "__main__":
    main()
