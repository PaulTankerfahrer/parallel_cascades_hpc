#!/bin/bash
#SBATCH -J "cascade_weak"
#SBATCH -N 2
#SBATCH -n 256
#SBATCH --mem-per-cpu 1G
#SBATCH -t 01:00:00
#SBATCH -p normal
##SBATCH -A hpc-prf-simvine
# ============================================================================
#  job_weak.sh -- WEAK SCALING: konstante Last PRO Rank, Gesamtgroesse waechst
#  mit der Rankzahl. Ideal: walltime bleibt flach. Der Anstieg = reiner
#  Kommunikations-/Overhead-Anteil. 1D-Zerlegung: NX fest, NY = rows*P.
#  Danach:  python3 analyze_weak.py weak.csv
# ============================================================================
module reset
module load mpi/OpenMPI/5.0.8-GCC-14.3.0

OUT=weak.csv
NX_FIX=1024          # feste Breite
ROWS=128             # Zeilen PRO Rank  -> 1024*128 = 131072 Atome/Rank (konstant)

echo "== Build (znver3) =="
mpicc -O3 -march=znver3 -mprefer-vector-width=256 -o mpi_cascade mpi_cascade.c -lm

echo "ranks,walltime,force,cells,hash,halo,migrate,energy,comm" > "$OUT"
for NP in 8 16 32 64 128 256; do
    NY=$(( ROWS * NP ))
    cat > _weak_$NP.ini <<EOF
[grid]
NX = $NX_FIX
NY = $NY
L0 = 1.0
[potential]
K_SPRING=100.0
MAX_STRETCH=1.15
RCUT=0.9
K_REP=400.0
REP_N=12.0
[pka]
pka_x=0.5
pka_y=0.5
pka_energy=5000.0
pka_angle=15.0
pka_mass=1.0
n_pka=1
timing=true
[dynamics]
dt=2e-4
n_steps=3000
[output]
log_every=3000
dump_every=0
out_prefix=weak_$NP
[ensemble]
seed=3
EOF
    echo "--- ranks=$NP  NY=$NY  (N=$(( NX_FIX*NY )), $(( NX_FIX*ROWS )) Atome/Rank) ---"
    srun -n $NP -N 2 \
         ./mpi_cascade _weak_$NP.ini \
         | grep "^# DATA" | sed 's/^# DATA //' >> "$OUT"
    rm -f _weak_$NP.ini
done
echo "fertig -> $OUT"
