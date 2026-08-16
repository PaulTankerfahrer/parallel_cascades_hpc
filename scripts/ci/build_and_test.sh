#!/usr/bin/env bash
# build_and_test.sh -- CI: Simulator bauen + Physik-Smoke-Test.
#
# Prueft billig und schnell (kein Cluster/GPU noetig):
#   1. cascade_serial.c und mpi_cascade.c kompilieren
#   2. Energieerhaltung: |drift%| < DRIFT_MAX auf kleinem Gitter
#   3. Determinismus: gleicher Seed -> gleiche broken_bonds (seriell, 2 Laeufe)
#   4. Sanity: broken_bonds > 0 (es ist ueberhaupt eine Kaskade passiert)
#   5. MPI-Smoke: 2 Ranks laufen ohne Absturz, geben RESULT aus
#
# Lokal ausfuehrbar:  bash scripts/ci/build_and_test.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CFG="$ROOT/src/params_valgrind.ini"   # 80x80, 300 Schritte -- klein & schnell
DRIFT_MAX=0.5                          # Prozent; Lauf liegt real bei ~0.04 %
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

echo "== [1/5] Kompilieren =="
gcc -O2 -o cascade_serial "$ROOT/src/cascade_serial.c" -lm
mpicc -O2 -o mpi_cascade  "$ROOT/src/mpi_cascade.c"    -lm
echo "   ok: seriell + MPI gebaut"

echo "== [2/5] Energieerhaltung =="
out1="$(./cascade_serial "$CFG")"
drift="$(printf '%s\n' "$out1" | awk '/^ *[0-9]/{d=$NF} END{print d+0}')"
echo "   drift% = $drift  (Grenze +/-$DRIFT_MAX)"
awk -v d="$drift" -v m="$DRIFT_MAX" 'BEGIN{ if (d<0) d=-d; exit !(d<m) }' \
  || { echo "   FEHLER: Energiedrift zu gross"; exit 1; }
echo "   ok"

echo "== [3/5] Determinismus =="
bb1="$(printf '%s\n' "$out1" | grep -oE 'broken_bonds=[0-9]+' | head -1)"
out2="$(./cascade_serial "$CFG")"
bb2="$(printf '%s\n' "$out2" | grep -oE 'broken_bonds=[0-9]+' | head -1)"
echo "   Lauf1: $bb1 | Lauf2: $bb2"
[ "$bb1" = "$bb2" ] || { echo "   FEHLER: nicht deterministisch"; exit 1; }
echo "   ok"

echo "== [4/5] Sanity (broken_bonds > 0) =="
n="${bb1#broken_bonds=}"
[ "${n:-0}" -gt 0 ] || { echo "   FEHLER: keine gerissenen Bindungen"; exit 1; }
echo "   ok: $n gerissene Bindungen"

echo "== [5/5] MPI-Smoke (2 Ranks) =="
# CI-Container laufen als root -> --allow-run-as-root; --oversubscribe fuer <2 Slots
mout="$(mpirun --allow-run-as-root --oversubscribe -np 2 ./mpi_cascade "$CFG" 2>&1)"
printf '%s\n' "$mout" | grep -q 'RESULT' \
  || { echo "   FEHLER: MPI-Lauf ohne RESULT"; printf '%s\n' "$mout" | tail -5; exit 1; }
printf '%s\n' "$mout" | grep -oE 'RESULT.*'
echo "   ok"

echo "== Alle Physik-Smoke-Tests bestanden =="
