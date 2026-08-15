#!/bin/bash
#SBATCH -J "cascade_strong"
#SBATCH -N 2
#SBATCH -n 256
#SBATCH --mem-per-cpu 1G
#SBATCH -t 01:00:00
#SBATCH -p normal
##SBATCH -A hpc-prf-simvine
# ============================================================================
#  job_scaling.sh -- STRONG SCALING: feste Problemgroesse (immer gleiche Paramter),
#  wachsende Rankzahl.
#  Sammelt die "# DATA"-Zeilen in scaling.csv.
#  Aufruf (im sbatch):  sbatch job_scaling.sh [params] [out.csv]
#  Danach:              python3 analyze_scaling.py scaling.csv
# ============================================================================
module reset
module load mpi/OpenMPI/5.0.8-GCC-14.3.0    # ggf. anpassen ('ml av' zeigt Module)
 
CFG=${1:-params_scaling.ini}
OUT=${2:-scaling.csv}
OUT2="${OUT%.csv}_2node.csv"
HDR="ranks,walltime,force,cells,hash,halo,migrate,energy,comm"
 
echo "== Build (znver3) =="
mpicc -O3 -march=znver3 -o mpi_cascade mpi_cascade.c -lm
 
# --- 1 Knoten: 1..128 ---
echo "$HDR" > "$OUT"
for NP in 1 2 4 8 16 32 64 128; do
    echo "--- ranks = $NP (1 Knoten) ---"
    srun -N 1 -n $NP --cpu-bind=cores ./mpi_cascade "$CFG" \
        | grep "^# DATA" | sed 's/^# DATA //' >> "$OUT"
done
 
# --- 2 Knoten: 64/128/256 (Inter-Node) ---
echo "$HDR" > "$OUT2"
for NP in 64 128 256; do
    echo "--- ranks = $NP (2 Knoten) ---"
    srun -N 2 -n $NP --cpu-bind=cores --distribution=block ./mpi_cascade "$CFG" \
        | grep "^# DATA" | sed 's/^# DATA //' >> "$OUT2"
done
 
echo "fertig -> $OUT  und  $OUT2"

