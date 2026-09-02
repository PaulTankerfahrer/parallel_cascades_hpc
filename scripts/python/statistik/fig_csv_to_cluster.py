#!/usr/bin/env python3
"""
fig_csv_to_cluster.py -- Abbildung fuer die Physik-Praesentation:
zeigt den Weg von der Endzustands-CSV eines einzelnen Laufs zu den
Clustergroessen, die spaeter ueber 1000 Laeufe gepoolt werden.

Links : Spalte `disp` aller Atome (Verschiebung vom Gitterplatz).
Rechts: nach Schwelle disp > 0.5*L0 und Nachbarschaftsverknuepfung
        (< 1.5*L0) -- jeder Cluster eine Farbe, Groesse direkt angeschrieben.

Aufruf:
  python3 fig_csv_to_cluster.py ARCHIV.zip ensemble/run_0116_state_final.csv AUSGABE.png
"""
import io, sys, zipfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

L0 = 1.0
DISP_THRESH = 0.5 * L0
LINK = 1.5 * L0

# Kategoriale Slots 1-3 (validiert fuer all-pairs), Rest neutral.
CAT = ["#2a78d6", "#eb6834", "#1baf7a"]
REST = "#8c8c88"
INK = "#1f1f1e"
INK2 = "#5c5c58"
SURFACE = "#fcfcfb"
# Sequenzieller Ein-Ton-Verlauf (Blau-Rampe, hell -> dunkel)
SEQ = LinearSegmentedColormap.from_list(
    "blues1", ["#e8f0fc", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#123a70"])


def labels(px, py, link=LINK):
    """Union-Find ueber Nachbarpaare; gibt fuer jedes Atom eine Cluster-ID."""
    n = len(px)
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    buckets = {}
    for i in range(n):
        buckets.setdefault((int(px[i] // link), int(py[i] // link)), []).append(i)
    link2 = link * link
    for i in range(n):
        cx, cy = int(px[i] // link), int(py[i] // link)
        for ax in (cx - 1, cx, cx + 1):
            for ay in (cy - 1, cy, cy + 1):
                for j in buckets.get((ax, ay), ()):
                    if j <= i:
                        continue
                    if (px[i] - px[j]) ** 2 + (py[i] - py[j]) ** 2 < link2:
                        ra, rb = find(i), find(j)
                        if ra != rb:
                            parent[ra] = rb
    return np.array([find(i) for i in range(n)])


def read(zpath, member):
    x, y, d = [], [], []
    with zipfile.ZipFile(zpath) as z, z.open(member) as fh:
        txt = io.TextIOWrapper(fh, encoding="utf-8", errors="replace")
        txt.readline()
        for line in txt:
            p = line.split(",")
            x.append(float(p[0])); y.append(float(p[1])); d.append(float(p[2]))
    return np.array(x), np.array(y), np.array(d)


def main():
    zpath, member, out = sys.argv[1], sys.argv[2], sys.argv[3]
    x, y, d = read(zpath, member)
    m = d > DISP_THRESH
    lab = labels(x[m], y[m])
    ids, counts = np.unique(lab, return_counts=True)
    order = np.argsort(-counts)
    ids, counts = ids[order], counts[order]

    # Ausschnitt um die Schadenszone
    pad = 10.0
    x0, x1 = x[m].min() - pad, x[m].max() + pad
    y0, y1 = y[m].min() - pad, y[m].max() + pad
    win = (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)

    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.6), facecolor=SURFACE)
    for a in ax:
        a.set_facecolor(SURFACE)
        a.set_xlim(x0, x1); a.set_ylim(y0, y1)
        a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
        for s in a.spines.values():
            s.set_color("#dcdcd6")

    # --- links: rohe CSV-Spalte disp ---------------------------------
    dv = d[win]
    sc = ax[0].scatter(x[win], y[win], c=dv, cmap=SEQ, s=3.2,
                       vmin=0.0, vmax=np.percentile(d[m], 97), linewidths=0)
    cb = fig.colorbar(sc, ax=ax[0], fraction=0.046, pad=0.02)
    cb.set_label("disp  (Verschiebung in $L_0$)", color=INK2, fontsize=11)
    cb.ax.tick_params(colors=INK2, labelsize=10)
    cb.outline.set_edgecolor("#dcdcd6")
    cb.ax.axhline(DISP_THRESH, color=INK, lw=1.4)
    cb.ax.annotate("Schwelle", xy=(0.5, DISP_THRESH), xycoords=("axes fraction", "data"),
                   xytext=(0, 7), textcoords="offset points", ha="center",
                   color=INK, fontsize=10)
    ax[0].set_title("Spalte $\\mathtt{disp}$ aus der Endzustands-CSV",
                    color=INK, fontsize=14, pad=9)
    ax[0].text(0.02, 0.97, f"{len(x):d} Atome – je eine Zeile".replace("62500", "62 500"),
               transform=ax[0].transAxes, va="top", ha="left",
               color=INK2, fontsize=11)

    # --- rechts: Cluster ---------------------------------------------
    ax[1].scatter(x[win], y[win], c="#e6e6e0", s=2.0, linewidths=0)
    xs, ys = x[m], y[m]
    mx0, my0 = xs.mean(), ys.mean()
    rest_sizes = []
    for k, (cid, cnt) in enumerate(zip(ids, counts)):
        sel = lab == cid
        col = CAT[k] if k < len(CAT) else REST
        ax[1].scatter(xs[sel], ys[sel], c=col, s=4.5, linewidths=0,
                      zorder=3 if k < len(CAT) else 2)
        cx, cy = xs[sel].mean(), ys[sel].mean()
        if k < len(CAT):
            # Label radial nach aussen vom Schadenszentrum, mit Fuehrungslinie
            vx, vy = cx - mx0, cy - my0
            norm = max((vx * vx + vy * vy) ** 0.5, 1e-9)
            reach = 0.19 * min(x1 - x0, y1 - y0)
            lx, ly = mx0 + vx / norm * (norm + reach), my0 + vy / norm * (norm + reach)
            ax[1].annotate(f"$n={cnt}$", xy=(cx, cy), xytext=(lx, ly),
                           ha="center", va="center", color=INK, fontsize=13,
                           arrowprops=dict(arrowstyle="-", color="#9a9a95", lw=0.9,
                                           shrinkA=6, shrinkB=4),
                           bbox=dict(boxstyle="round,pad=0.25", fc=SURFACE,
                                     ec="#dcdcd6", lw=0.8), zorder=5)
        else:
            rest_sizes.append(int(cnt))
    ax[1].set_title("Nach Schwelle und Nachbarschaft: Cluster",
                    color=INK, fontsize=14, pad=9)
    txt = (f"{int(m.sum())} verlagerte Atome  →  {len(ids)} Cluster\n"
           f"Größen: " + ", ".join(str(int(c)) for c in counts))
    ax[1].text(0.02, 0.97, txt, transform=ax[1].transAxes, va="top", ha="left",
               color=INK2, fontsize=11)

    fig.tight_layout()
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print("geschrieben:", out, "| Cluster:", list(map(int, counts)))


if __name__ == "__main__":
    main()
