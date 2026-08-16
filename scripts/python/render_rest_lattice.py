#!/usr/bin/env python3
"""render_rest_lattice.py -- rendert das Modell in Ruhelage: ein Ausschnitt des
Dreiecksgitters mit Federbindungen (Frame 0 eines XYZ-Dumps).

Erzeugt das Standbild fuer die "Das Modell"-Folie der Physik-Praesentation.
Braucht das OVITO-Python-Modul (pip install ovito).

Aufruf:
  ovito-venv/bin/python render_rest_lattice.py <single.xyz> <ausgabe.png>
"""
import os, sys
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer
from ovito.modifiers import (ExpressionSelectionModifier, DeleteSelectedModifier,
                             CreateBondsModifier, ComputePropertyModifier)

xyz = sys.argv[1] if len(sys.argv) > 1 else "single.xyz"
out = sys.argv[2] if len(sys.argv) > 2 else "modell_ruhelage.png"

pl = import_file(xyz)
pl.add_to_scene()
pl.source.data.particles_.vis.radius = 0.34

# zentralen Gitterausschnitt (~16 x 14 L0) behalten
pl.modifiers.append(ExpressionSelectionModifier(
    expression="Position.X<62 || Position.X>78 || Position.Y<52 || Position.Y>66"))
pl.modifiers.append(DeleteSelectedModifier())

# einheitliche Atomfarbe (Stahlblau)
pl.modifiers.append(ComputePropertyModifier(
    operate_on="particles", output_property="Color", expressions=["0.27", "0.45", "0.78"]))

# Federbindungen (naechste Nachbarn, r=1.0); Breite ueber cb.vis, Farbe als
# Per-Bond-Property erzwingen (cb.vis.color wird beim Rendern ignoriert).
cb = CreateBondsModifier(cutoff=1.15)
pl.modifiers.append(cb)
cb.vis.width = 0.09
pl.modifiers.append(ComputePropertyModifier(
    operate_on="bonds", output_property="Color", expressions=["0.34", "0.34", "0.40"]))

vp = Viewport(type=Viewport.Type.Top)
vp.zoom_all()
r = TachyonRenderer(ambient_occlusion=True, ambient_occlusion_brightness=0.85, shadows=False)
vp.render_image(filename=out, size=(900, 760), frame=0, renderer=r, background=(1, 1, 1))
print("geschrieben:", out)
