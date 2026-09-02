#!/usr/bin/env python3
"""make_cascade_videos.py -- rendert die vier Physik-Videos der Kollisionskaskade
mit OVITO (headless, Tachyon) und kodiert sie via ffmpeg zu MP4.

Voraussetzungen:
  * OVITO-Python-Modul (pip install ovito) -- die System-GUI reicht NICHT.
  * ffmpeg im PATH.
  * Zwei XYZ-Dumps, erzeugt vom seriellen Simulator (dump_every>0):
      single.xyz  -- ein PKA  (Videos 1-3, nur andere Einfaerbung)
      multi.xyz   -- fuenf PKAs (Video 4)

Aufruf:
  ovito-venv/bin/python make_cascade_videos.py <xyz-verzeichnis> <ausgabe-verzeichnis>
  [--only=3_cluster,1_geschwindigkeit]   nur einzelne Videos rendern

Die vier Videos:
  1. Geschwindigkeitsfeld  (vmag)         -- der heisse Stoss-Spike expandiert und kuehlt ab
  2. Schadensbildung       (broken bonds) -- Defektcluster friert ins Gitter ein
  3. Defektcluster         (disp + Cluster) -- die Objekte, die ausgewertet werden
  4. Mehrere Einschlaege   (vmag, 5 PKAs) -- ueberlappende Kaskaden + Channeling

Video 3 setzt die Clusterdefinition der Auswertung um: Atome mit
disp > 1,0 L0 gelten als verlagert, und alles, was naeher als 1,0 L0
beieinander liegt, ist ein Cluster (results/statistik/fit_fest/BEFUND_fest.md).
Jeder Cluster bekommt eine eigene Farbe, der Rest des Gitters bleibt grau.
"""
import os, sys, subprocess, tempfile, shutil
import numpy as np

from ovito.io import import_file
from ovito.vis import (Viewport, TachyonRenderer,
                       ColorLegendOverlay, TextLabelOverlay)
from ovito.modifiers import (ColorCodingModifier, ExpressionSelectionModifier,
                             ClusterAnalysisModifier)
from ovito.qt_compat import QtCore

# Der Simulator schreibt "disp" als Skalar; ohne explizite Zuordnung legt OVITO
# den Wert in die X-Komponente der Standard-Eigenschaft "Displacement".
COLUMNS = ["Particle Type", "Position.X", "Position.Y", "Position.Z",
           "disp", "vmag", "broken"]

# Clusterdefinition -- identisch zur Auswertung
DISP_THRESH = 1.0
LINK = 1.0

# Kategoriale Farben fuer die Cluster, hell genug fuer dunklen Hintergrund.
CLUSTER_COLORS = np.array([
    (0x4a, 0x9e, 0xff), (0xff, 0x8a, 0x4a), (0x3f, 0xd9, 0x9b),
    (0xff, 0xd2, 0x4a), (0xc7, 0x7d, 0xff), (0xff, 0x6b, 0x8a),
    (0x4a, 0xe0, 0xe0), (0xa3, 0xe6, 0x35), (0xff, 0x5f, 0x5f),
    (0x7a, 0xa5, 0xff), (0xff, 0xb3, 0xd9), (0xb8, 0xf2, 0x7a),
], dtype=float) / 255.0
LATTICE_COLOR = (0.17, 0.17, 0.21)     # nicht verlagerte Atome
R_LATTICE, R_CLUSTER = 0.30, 0.62      # Radien: Cluster treten hervor


def damage_view(pl, thresh, margin):
    """Kameraausschnitt aus der Schadenszone des Endzustands.

    Bei NX=250 fuellt die Kaskade nur einen kleinen Teil der Box -- zoom_all()
    wuerde sie auf wenige Pixel schrumpfen. Ausschnitt und Mitte werden daher
    aus den verlagerten Atomen des letzten Frames bestimmt.
    """
    d = pl.compute(pl.source.num_frames - 1)
    disp = np.asarray(d.particles["disp"])
    pos = np.asarray(d.particles["Position"])
    q = pos[disp > thresh]
    if len(q) == 0:
        return None, None
    cx = (q[:, 0].min() + q[:, 0].max()) / 2
    cy = (q[:, 1].min() + q[:, 1].max()) / 2
    half = max(q[:, 0].max() - q[:, 0].min(),
               q[:, 1].max() - q[:, 1].min()) / 2 * margin
    return (cx, cy), half


def color_by_cluster(frame, data):
    """Faerbt jeden Cluster eigenstaendig ein, alles darunter bleibt grau."""
    cl = np.asarray(data.particles["Cluster"])
    col = np.empty((len(cl), 3), dtype=float)
    col[:] = LATTICE_COLOR
    rad = np.full(len(cl), R_LATTICE)
    m = cl > 0
    col[m] = CLUSTER_COLORS[(cl[m] - 1) % len(CLUSTER_COLORS)]
    rad[m] = R_CLUSTER
    data.particles_.create_property("Color", data=col)
    data.particles_.create_property("Radius", data=rad)

A = QtCore.Qt.AlignmentFlag
H = QtCore.Qt.Orientation.Horizontal

# --- Job-Definitionen -------------------------------------------------------
# (xyz, property, lo, hi, gradient, dateiname, titel, legenden-titel, hintergrund)
def jobs():
    dark = (0.02, 0.02, 0.06)
    black = (0.0, 0.0, 0.0)
    return [
        dict(xyz="single.xyz", prop="vmag", lo=0.0, hi=0.6,
             grad=ColorCodingModifier.Jet(), name="1_geschwindigkeit",
             title="Stosskaskade - Geschwindigkeitsfeld",
             legend="Geschwindigkeit |v|", bg=dark, view=(0.2, 1.35)),
        dict(xyz="single.xyz", prop="broken", lo=0.0, hi=6.0,
             grad=ColorCodingModifier.Hot(), name="2_schaden",
             title="Schadensbildung - gerissene Bindungen",
             legend="fehlende Bindungen", bg=black, view=(0.2, 1.35)),
        dict(xyz="single.xyz", prop=None, name="3_cluster",
             title="Defektcluster - disp > 1,0, Linkradius 1,0",
             legend=None, bg=dark, cluster=True, view=(1.0, 2.0)),
        dict(xyz="multi.xyz", prop="vmag", lo=0.0, hi=0.6,
             grad=ColorCodingModifier.Jet(), name="4_mehrere_pka",
             title="Mehrere Einschlaege (5 PKAs)",
             legend="Geschwindigkeit |v|", bg=dark),
    ]

SIZE = (800, 800)
RADIUS = 0.55
FPS = 20
TOTAL_SECONDS = 15.0     # Gesamtlaenge je Video; das Endbild wird so lange
                         # eingefroren, dass Animation + Freeze == TOTAL_SECONDS


def set_view(vp, pl, j):
    """Kamera auf die Schadenszone setzen, sonst die ganze Box zeigen."""
    if "view" not in j:
        vp.zoom_all()
        return
    thresh, margin = j["view"]
    center, half = damage_view(pl, thresh, margin)
    if center is None:
        vp.zoom_all()
        return
    vp.camera_pos = (center[0], center[1], 0.0)
    vp.camera_dir = (0.0, 0.0, -1.0)
    vp.fov = half


def render_job(j, xyz_dir, out_dir, renderer):
    pl = import_file(os.path.join(xyz_dir, j["xyz"]), columns=COLUMNS)
    pl.add_to_scene()
    pl.source.data.particles_.vis.radius = RADIUS

    vp = Viewport(type=Viewport.Type.Top)

    if j.get("cluster"):
        # Clusterdefinition der Auswertung, Schritt fuer Schritt
        pl.modifiers.append(ExpressionSelectionModifier(
            expression=f"disp > {DISP_THRESH}"))
        pl.modifiers.append(ClusterAnalysisModifier(
            cutoff=LINK, only_selected=True, sort_by_size=True))
        pl.modifiers.append(color_by_cluster)
        set_view(vp, pl, j)
        vp.overlays.append(TextLabelOverlay(
            text="[ClusterAnalysis.cluster_count] getrennte Cluster  ·  groesster: [ClusterAnalysis.largest_size]",
            alignment=A.AlignBottom | A.AlignHCenter,
            text_color=(1, 1, 1), font_size=0.038,
            offset_y=0.03, source_pipeline=pl))
    else:
        cc = ColorCodingModifier(property=j["prop"],
                                 start_value=j["lo"], end_value=j["hi"])
        cc.gradient = j["grad"]
        pl.modifiers.append(cc)
        set_view(vp, pl, j)
        vp.overlays.append(ColorLegendOverlay(
            modifier=cc, title=j["legend"],
            alignment=A.AlignBottom | A.AlignHCenter, orientation=H,
            text_color=(1, 1, 1), font_size=0.04, legend_size=0.28,
            offset_y=0.04))

    vp.overlays.append(TextLabelOverlay(
        text=j["title"], alignment=A.AlignTop | A.AlignLeft,
        text_color=(1, 1, 1), font_size=0.045, offset_x=0.02, offset_y=0.02))

    nframes = pl.source.num_frames
    tmp = tempfile.mkdtemp(prefix="frames_")
    for fr in range(nframes):
        vp.render_image(filename=os.path.join(tmp, f"f_{fr:04d}.png"),
                        size=SIZE, frame=fr, renderer=renderer, background=j["bg"])
    pl.remove_from_scene()

    out = os.path.join(out_dir, f"cascade_{j['name']}.mp4")
    # H.264, yuv420p (breit kompatibel). Endbild so lange einfrieren, dass
    # Animation + Freeze == TOTAL_SECONDS (Animation selbst bleibt unveraendert).
    freeze = max(0.0, TOTAL_SECONDS - nframes / FPS)
    subprocess.run([
        "ffmpeg", "-y", "-loglevel", "error",
        "-framerate", str(FPS), "-i", os.path.join(tmp, "f_%04d.png"),
        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
        "-vf", f"tpad=stop_mode=clone:stop_duration={freeze:.3f}",
        out], check=True)
    shutil.rmtree(tmp)
    print(f"  -> {out}  ({nframes} Frames)")
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    only = next((a.split("=", 1)[1].split(",") for a in sys.argv[1:]
                 if a.startswith("--only=")), None)
    xyz_dir = args[0] if len(args) > 0 else "."
    out_dir = args[1] if len(args) > 1 else "videos"
    os.makedirs(out_dir, exist_ok=True)
    renderer = TachyonRenderer(ambient_occlusion=False, shadows=False)
    for j in jobs():
        if only and j["name"] not in only:
            continue
        print(f"[{j['name']}] rendere {j['xyz']} ({j['prop']}) ...")
        render_job(j, xyz_dir, out_dir, renderer)
    print("fertig.")


if __name__ == "__main__":
    main()
