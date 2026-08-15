#!/bin/bash
#SBATCH -J "cascade_strong_rep"
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 02:00:00
#SBATCH --mem-per-cpu 1G
#SBATCH -p normal
##SBATCH -A hpc-prf-simvine
# ============================================================================
#  job_scaling_rep.sh -- Strong Scaling MIT WIEDERHOLUNGEN (gegen Rauschen).
#  Faehrt den kompletten Sweep REPS-mal; jede Rankzahl bekommt REPS Zeilen.
#  Danach Median bilden und auswerten:
#      python3 aggregate_reps.py scaling_raw.csv scaling.csv
#      python3 analyze_scaling.py scaling.csv
#  Fuer den Off-vs-On-Vergleich: zweimal laufen (params mit Opt. aus / an),
#  z. B. -> scaling_off_raw.csv und scaling_on_raw.csv, beide aggregieren,
#  dann:  python3 analyze_opt_compare.py scaling_off.csv "AUS" \
#                                        scaling_on.csv "non-blocking + resort"
#  Aufruf:  sbatch job_scaling_rep.sh [params] [out_raw.csv] [reps]
# ============================================================================
module reset
module load mpi/OpenMPI/5.0.8-GCC-14.3.0

CFG=${1:-params_scaling.ini}
OUT=${2:-scaling_raw.csv}
REPS=${3:-3}

echo "== Build (znver3) =="
mpicc -O3 -march=znver3 -o mpi_cascade mpi_cascade.c -lm

echo "ranks,walltime,force,cells,hash,halo,migrate,energy,comm" > "$OUT"
for rep in $(seq 1 "$REPS"); do
    echo "=== Wiederholung $rep / $REPS ==="
    for NP in 1 2 4 8 16 32 64 128; do
        echo "--- ranks = $NP (rep $rep) ---"
        srun -N 1 -n $NP --cpu-bind=cores ./mpi_cascade "$CFG" \
            | grep "^# DATA" | sed 's/^# DATA //' >> "$OUT"
    done
done
echo "fertig -> $OUT  ($REPS Wiederholungen je Rankzahl)"
