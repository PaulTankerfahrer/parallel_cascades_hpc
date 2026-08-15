# Physik-Präsentation & Kaskaden-Videos

Kurze, physikorientierte Vorstellung der Kollisionskaskaden-Simulation
(ohne den HPC-/Programmierteil) samt vier OVITO-Videos.

## Inhalt

| Datei | Zweck |
|---|---|
| `physik_praesentation.tex` / `.pdf` | Beamer-Folien (Physik, ~9 Folien) |
| `videos/cascade_1_geschwindigkeit.mp4` | Geschwindigkeitsfeld — thermischer Spike |
| `videos/cascade_2_schaden.mp4` | Schadensbildung — gerissene Bindungen |
| `videos/cascade_3_verschiebung.mp4` | Verschiebungsfeld — Heat-Spike |
| `videos/cascade_4_mehrere_pka.mp4` | Fünf gleichzeitige Einschläge |
| `stills/` | Standbilder für die Folien (aus den Videos extrahiert) |
| `single.ini`, `multi.ini` | Simulations-Konfigurationen der beiden Läufe |
| `build_videos.sh` | Erzeugt alles neu (Sim → Render → Standbilder) |

## Videos abspielen

Die Videos sind **eigenständige MP4-Dateien**. In der PDF ist jede Video-Folie
mit einem Standbild hinterlegt; ein Klick darauf (bzw. auf „▶ Video abspielen``)
startet in den meisten PDF-Viewern den Systemplayer. Alternativ die MP4s direkt
aus `videos/` öffnen.

## Neu erzeugen

```bash
# OVITO-Python bereitstellen (einmalig):
python3 -m venv ovito-venv && ovito-venv/bin/pip install ovito

# Videos + Standbilder erzeugen:
OVITO_PY=ovito-venv/bin/python ./build_videos.sh

# Folien bauen:
pdflatex physik_praesentation.tex   # 2x für Referenzen
```

> Hinweis: Es genügt das PyPI-Paket `ovito` (headless-Rendering).
> Die OVITO-Desktop-GUI allein bringt **kein** skriptbares Python mit.

## Physikalische Kernaussagen

- Ein PKA (Primäres Rückstoßatom) startet eine verzweigende Kollisionskaskade.
- Sichtbar werden: heißer Spike (Geschwindigkeit), Defektbildung (gerissene
  Bindungen) und bleibende Verschiebung.
- Über 1000 Kaskaden folgt die Verteilung der Defekt-Clustergrößen einem
  **Potenzgesetz** $F(n)\propto n^{-S}$ mit $S\approx 1{,}4$ — skalenfrei und
  robust gegen Defektheilung. Vergleichswert 3D-Wolfram (Sand et al. 2013):
  $S=1{,}63$.
