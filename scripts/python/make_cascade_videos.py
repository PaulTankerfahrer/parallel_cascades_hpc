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

Die vier Videos:
  1. Geschwindigkeitsfeld  (vmag)         -- der heisse Stoss-Spike expandiert und kuehlt ab
  2. Schadensbildung       (broken bonds) -- Defektcluster friert ins Gitter ein
  3. Verschiebungsfeld     (disp)         -- bleibende Verformung / Heat-Spike-Halo
  4. Mehrere Einschlaege   (vmag, 5 PKAs) -- ueberlappende Kaskaden + Channeling
"""
import os, sys, subprocess, tempfile, shutil

from ovito.io import import_file
from ovito.vis import (Viewport, TachyonRenderer,
                       ColorLegendOverlay, TextLabelOverlay)
from ovito.modifiers import ColorCodingModifier
from ovito.qt_compat import QtCore

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
             legend="Geschwindigkeit |v|", bg=dark),
        dict(xyz="single.xyz", prop="broken", lo=0.0, hi=6.0,
             grad=ColorCodingModifier.Hot(), name="2_schaden",
             title="Schadensbildung - gerissene Bindungen",
             legend="fehlende Bindungen", bg=black),
        dict(xyz="single.xyz", prop="Displacement", lo=0.0, hi=3.0,
             grad=ColorCodingModifier.Viridis(), name="3_verschiebung",
             title="Verschiebungsfeld - Heat-Spike",
             legend="Verschiebung |u|", bg=dark),
        dict(xyz="multi.xyz", prop="vmag", lo=0.0, hi=0.6,
             grad=ColorCodingModifier.Jet(), name="4_mehrere_pka",
             title="Mehrere Einschlaege (5 PKAs)",
             legend="Geschwindigkeit |v|", bg=dark),
    ]

SIZE = (800, 800)
RADIUS = 0.55
FPS = 20


def render_job(j, xyz_dir, out_dir, renderer):
    pl = import_file(os.path.join(xyz_dir, j["xyz"]))
    pl.add_to_scene()
    pl.source.data.particles_.vis.radius = RADIUS
    cc = ColorCodingModifier(property=j["prop"],
                             start_value=j["lo"], end_value=j["hi"])
    cc.gradient = j["grad"]
    pl.modifiers.append(cc)

    vp = Viewport(type=Viewport.Type.Top)
    vp.zoom_all()
    vp.overlays.append(ColorLegendOverlay(
        modifier=cc, title=j["legend"],
        alignment=A.AlignBottom | A.AlignHCenter, orientation=H,
        text_color=(1, 1, 1), font_size=0.04, legend_size=0.28, offset_y=0.04))
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
    # H.264, yuv420p (breit kompatibel), letztes Bild 1,5 s einfrieren
    subprocess.run([
        "ffmpeg", "-y", "-loglevel", "error",
        "-framerate", str(FPS), "-i", os.path.join(tmp, "f_%04d.png"),
        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
        "-vf", f"tpad=stop_mode=clone:stop_duration=1.5",
        out], check=True)
    shutil.rmtree(tmp)
    print(f"  -> {out}  ({nframes} Frames)")
    return out


def main():
    xyz_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "videos"
    os.makedirs(out_dir, exist_ok=True)
    renderer = TachyonRenderer(ambient_occlusion=False, shadows=False)
    for j in jobs():
        print(f"[{j['name']}] rendere {j['xyz']} ({j['prop']}) ...")
        render_job(j, xyz_dir, out_dir, renderer)
    print("fertig.")


if __name__ == "__main__":
    main()
