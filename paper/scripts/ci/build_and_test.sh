#!/usr/bin/env bash
# build_and_test.sh -- CI: Simulator bauen + Physik-Smoke-Test.
#
# Prueft billig und schnell (kein Cluster/GPU noetig):
#   1. cascade_serial.c kompilieren (MPI testet kursprojekt/scripts/ci)
#   2. Energieerhaltung: |drift%| < DRIFT_MAX auf kleinem Gitter
#   3. Determinismus: gleicher Seed -> gleiche broken_bonds (seriell, 2 Laeufe)
#   4. Sanity: broken_bonds > 0 (es ist ueberhaupt eine Kaskade passiert)
#   5. ZBL-Streutest: Zwei-Koerper-Streuung gegen die klassische Theorie
#   6. ZBL im Gitter: Energieerhaltung mit adaptivem Zeitschritt bei E = 10^4
#   7. Elektronische Bremsung (Lindhard): E_total + E_elec bleibt erhalten,
#      und es wird ueberhaupt Energie abgegeben
#
# Lokal ausfuehrbar:  bash paper/scripts/ci/build_and_test.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CFG="$ROOT/src/params_valgrind.ini"   # 80x80, 300 Schritte -- klein & schnell
DRIFT_MAX=0.5                          # Prozent; Lauf liegt real bei ~0.04 %
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

echo "== [1/7] Kompilieren =="
gcc -O2 -o cascade_serial "$ROOT/src/cascade_serial.c" -lm
echo "   ok: seriell gebaut"

echo "== [2/7] Energieerhaltung =="
out1="$(./cascade_serial "$CFG")"
drift="$(printf '%s\n' "$out1" | awk '/^ *[0-9]/{d=$NF} END{print d+0}')"
echo "   drift% = $drift  (Grenze +/-$DRIFT_MAX)"
awk -v d="$drift" -v m="$DRIFT_MAX" 'BEGIN{ if (d<0) d=-d; exit !(d<m) }' \
  || { echo "   FEHLER: Energiedrift zu gross"; exit 1; }
echo "   ok"

echo "== [3/7] Determinismus =="
bb1="$(printf '%s\n' "$out1" | grep -oE 'broken_bonds=[0-9]+' | head -1)"
out2="$(./cascade_serial "$CFG")"
bb2="$(printf '%s\n' "$out2" | grep -oE 'broken_bonds=[0-9]+' | head -1)"
echo "   Lauf1: $bb1 | Lauf2: $bb2"
[ "$bb1" = "$bb2" ] || { echo "   FEHLER: nicht deterministisch"; exit 1; }
echo "   ok"

echo "== [4/7] Sanity (broken_bonds > 0) =="
n="${bb1#broken_bonds=}"
[ "${n:-0}" -gt 0 ] || { echo "   FEHLER: keine gerissenen Bindungen"; exit 1; }
echo "   ok: $n gerissene Bindungen"

echo "== [5/7] ZBL-Streutest =="
gcc -O2 -I"$ROOT/src" -o zbl_scatter "$ROOT/scripts/ci/zbl_scatter_test.c" -lm
./zbl_scatter > zbl_scatter.out || { cat zbl_scatter.out; exit 1; }
tail -1 zbl_scatter.out
echo "   ok"

echo "== [6/7] ZBL im Gitter: Energieerhaltung (adaptiver Zeitschritt) =="
cat > zbl.ini <<'EOF'
[grid]
NX = 80
NY = 80
[potential]
rep_model = zbl
[pka]
pka_energy = 10000
pka_angle  = 17.3
[dynamics]
dt       = 3e-4
dt_adapt = true
t_max    = 0.3
n_steps  = 1000000
[output]
log_every  = 100000
out_prefix = zbl
EOF
out3="$(./cascade_serial zbl.ini)"
drift3="$(printf '%s\n' "$out3" | awk '/^ *[0-9]/{d=$NF} END{print d+0}')"
echo "   drift% = $drift3  (Grenze +/-$DRIFT_MAX)"
awk -v d="$drift3" -v m="$DRIFT_MAX" 'BEGIN{ if (d<0) d=-d; exit !(d<m) }' \
  || { echo "   FEHLER: Energiedrift zu gross"; exit 1; }
echo "   ok"

echo "== [7/7] Elektronische Bremsung: Energiebilanz =="
sed 's/^out_prefix = zbl/out_prefix = zbl_el\nel_stopping = lindhard/' zbl.ini > zbl_el.ini
out4="$(./cascade_serial zbl_el.ini)"
drift4="$(printf '%s\n' "$out4" | awk '/^ *[0-9]/{d=$NF} END{print d+0}')"
elec="$(printf '%s\n' "$out4" | grep -oE 'E_elec=[0-9.]+' | cut -d= -f2)"
echo "   drift% = $drift4  E_elec = ${elec:-fehlt}"
awk -v d="$drift4" -v m="$DRIFT_MAX" -v e="${elec:-0}" 'BEGIN{ if (d<0) d=-d; exit !(d<m && e>0) }' \
  || { echo "   FEHLER: Energiebilanz mit Bremsung verletzt oder keine Bremsung"; exit 1; }
echo "   ok"

echo "== Alle Physik-Smoke-Tests bestanden =="
