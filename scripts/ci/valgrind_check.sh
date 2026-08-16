#!/usr/bin/env bash
# valgrind_check.sh -- CI: Speicherkorrektheit des seriellen Simulators.
#
# Baut mit Debug-Infos und laeuft unter Valgrind auf einem kleinen Gitter
# (80x80, 300 Schritte). --error-exitcode=1 laesst den Job bei jedem
# Speicherfehler oder definitiven Leak rot werden.
#
# Hinweis: Auf sehr aktuellem glibc (z.B. Rolling-Release-Distros) kann
# Valgrind beim Start scheitern ("function redirection ... mandatory").
# Auf dem Ubuntu-CI-Runner (stabiles glibc) laeuft es normal.
#
# Lokal ausfuehrbar:  bash scripts/ci/valgrind_check.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CFG="$ROOT/src/params_valgrind.ini"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

echo "== Kompilieren (mit -g) =="
gcc -O1 -g -o cascade_serial "$ROOT/src/cascade_serial.c" -lm

echo "== Valgrind (seriell, strict) =="
valgrind --error-exitcode=1 \
         --leak-check=full \
         --errors-for-leak-kinds=definite \
         --track-origins=yes \
         ./cascade_serial "$CFG"

echo "== Keine Speicherfehler / definitiven Leaks =="
