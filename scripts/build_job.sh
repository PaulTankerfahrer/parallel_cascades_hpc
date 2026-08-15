#!/bin/bash
#SBATCH -J build_executable
#SBATCH -N 1
#SBATCH -n 1                 # 1 Knoten, 128 Kerne -> 128 Kaskaden gleichzeitig
#SBATCH --mem-per-cpu 2G
#SBATCH -t 00:05:00
#SBATCH -p normal
##SBATCH -A hpc-prf-simvine

##### Dieser Job dient nur dazu, die serielle Version zu builden.

module reset
module load mpi/OpenMPI/5.0.8-GCC-14.3.0 

echo "Building serial cascade build_executable"

mpicc -O3 -lm -march=native -o cascade_serial cascade_serial.c

echo "build done"