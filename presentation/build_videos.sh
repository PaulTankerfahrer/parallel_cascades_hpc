#!/usr/bin/env bash
# build_videos.sh -- erzeugt die vier Physik-Videos komplett neu.
#
# Ablauf:  Simulator bauen -> zwei Kaskaden-Laeufe (single/multi) mit XYZ-Dump
#          -> OVITO-Rendering + ffmpeg-Kodierung -> Standbilder fuer die Folien.
#
# Voraussetzungen: gcc, ffmpeg und ein Python mit OVITO-Modul (pip install ovito).
# Den Pfad zum OVITO-Python ggf. unten anpassen (OVITO_PY).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
WORK="$HERE/_work"
OVITO_PY="${OVITO_PY:-python3}"     # z.B. /pfad/zu/ovito-venv/bin/python

mkdir -p "$WORK" "$HERE/videos" "$HERE/stills"

echo "[1/4] Simulator bauen"
gcc -O3 -march=native -o "$WORK/cascade_serial" "$ROOT/src/cascade_serial.c" -lm

echo "[2/4] Kaskaden simulieren (XYZ-Dumps)"
( cd "$WORK" && ./cascade_serial "$HERE/single.ini" )
( cd "$WORK" && ./cascade_serial "$HERE/multi.ini" )

echo "[3/4] Videos rendern (OVITO + ffmpeg)"
"$OVITO_PY" "$ROOT/scripts/python/make_cascade_videos.py" "$WORK" "$HERE/videos"

echo "[4/4] Standbilder fuer die Folien extrahieren (15-s-Videos)"
ffmpeg -y -loglevel error -ss 6.6  -i "$HERE/videos/cascade_1_geschwindigkeit.mp4" -frames:v 1 "$HERE/stills/still_1.png"
ffmpeg -y -loglevel error -ss 14.3 -i "$HERE/videos/cascade_2_schaden.mp4"        -frames:v 1 "$HERE/stills/still_2.png"
ffmpeg -y -loglevel error -ss 14.3 -i "$HERE/videos/cascade_3_verschiebung.mp4"   -frames:v 1 "$HERE/stills/still_3.png"
ffmpeg -y -loglevel error -ss 14.3 -i "$HERE/videos/cascade_4_mehrere_pka.mp4"    -frames:v 1 "$HERE/stills/still_4.png"

echo "fertig. Videos in videos/, Standbilder in stills/."
echo "Praesentation bauen:  pdflatex physik_praesentation.tex"
