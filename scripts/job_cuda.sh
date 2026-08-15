#!/bin/bash
#SBATCH --job-name=cuda_cascade
#SBATCH -A hpc-prf-simvine
#SBATCH -p gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 00:30:00
#SBATCH --mem 40G
#SBATCH --output=cuda_cascade_%j.out
#SBATCH --error=cuda_cascade_%j.err

set -euo pipefail

# --- Module ---------------------------------------------------------------
# PRUEFEN mit: module spider CUDA   (Modulname/Version auf Noctua 2 anpassen)
module reset
module load system/CUDA/13.2.0
nvcc --version
nvidia-smi --query-gpu=name,memory.total --format=csv

# --- Build ----------------------------------------------------------------
nvcc -O3 -arch=sm_80 -o cascade_cuda cascade_cuda.cu

# --- 1) Timing-Lauf (liefert den ERGEBNIS-BLOCK) --------------------------
echo "=============== TIMING-LAUF ==============="
./cascade_cuda_step5

# --- 2) ncu-Profiling des Force-Kernels -----------------------------------
# Nur 2 repraesentative k_force-Launches (Warmup ueberspringen) -> schnell.
# SpeedOfLight = % Peak Compute/Memory; MemoryWorkloadAnalysis = DRAM-Throughput.
echo "=============== NCU-PROFILING ==============="
ncu --target-processes all \
    --kernel-name k_force \
    --launch-skip 60 --launch-count 2 \
    --section SpeedOfLight \
    --section MemoryWorkloadAnalysis \
    --section Occupancy \
    -o ncu_force_%j \
    ./cascade_cuda | tee ncu_force_summary.txt

echo "fertig. Output: cuda_cascade_<jobid>.out  +  ncu_force_<jobid>.ncu-rep"
