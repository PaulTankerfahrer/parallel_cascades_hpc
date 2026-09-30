#!/bin/bash
#SBATCH -J cascade_ensemble
#SBATCH -N 1
#SBATCH -n 128                 # 1 Knoten, 128 Kerne -> 128 Kaskaden gleichzeitig
#SBATCH --mem-per-cpu 1G
#SBATCH -t 00:30:00
#SBATCH -p normal
##SBATCH -A hpc-prf-simvine
# ============================================================================
#  job_ensemble.sh -- 1000 unabhaengige Kaskaden auf EINEM Knoten.
#  xargs -P 128 haelt immer 128 Laeufe gleichzeitig am Laufen (embarrassingly
#  parallel). Jeder Lauf ist EINE serielle Kaskade auf 1 Kern.
#
#  VORHER auf dem Login-Knoten EINMAL bauen:
#      gcc -O3 -march=znver3 -o cascade_serial cascade_serial.c -lm
#  Submit:
#      sbatch job_ensemble.sh           # 1000 Laeufe (Default)
#      sbatch job_ensemble.sh 500       # oder weniger zum Testen
#  Wenn fertig (auf dem Login-Knoten):
#      python3 analyze_ensemble.py ensemble
# ============================================================================
module reset
module load lang/Python/3.13.5-GCCcore-14.3.0       # nur fuer run_one.py (Standardbibliothek reicht)

NRUNS=${1:-1000}
mkdir -p ensemble

echo "Starte $NRUNS Kaskaden, 128 gleichzeitig ..."
seq 0 $((NRUNS-1)) | xargs -P 128 -I{} python3 run_one.py {} ensemble

echo "fertig: $(ls ensemble/run_*_state_final.csv 2>/dev/null | wc -l) Laeufe abgeschlossen"
