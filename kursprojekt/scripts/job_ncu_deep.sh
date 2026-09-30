#!/bin/bash
#SBATCH --job-name=cuda_cascade
#SBATCH -A hpc-prf-simvine
#SBATCH -p gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 00:30:00
#SBATCH --mem 40G
#SBATCH --output=cuda_ncu_%j.out
#SBATCH --error=cuda_ncu_%j.err

set -euo pipefail
module reset
module load system/CUDA/13.2.0
nvcc --version

# fp64-Binary bauen (kleiner/kuerzer fuer ncu: 600 Schritte reichen, Profiling
# zielt ohnehin auf 2 repraesentative k_force-Launches nach Warmup)
nvcc -O3 -arch=sm_80 -o cascade_cuda_bench_fp64 cascade_cuda_bench.cu

# Tiefes Profiling: WARUM ist der Kernel so ausgelastet?
#   SpeedOfLight            -> Compute% vs Memory% (Bound-Typ)
#   ComputeWorkloadAnalysis -> welche Pipes (FP64? MUFU/transzendent?)
#   InstructionStats        -> Instruktionsmix (pow/sqrt-Anteil)
#   Occupancy               -> Auslastung der SMs
ncu --target-processes all \
    --kernel-name k_force \
    --launch-skip 60 --launch-count 2 \
    --section SpeedOfLight \
    --section ComputeWorkloadAnalysis \
    --section InstructionStats \
    --section Occupancy \
    --section MemoryWorkloadAnalysis \
    -o ncu_force_deep_${SLURM_JOB_ID} \
    ./cascade_cuda_bench_fp64 2000000 256 600 1 | tee ncu_force_deep_summary.txt

echo "fertig. Report: ncu_force_deep_${SLURM_JOB_ID}.ncu-rep  +  ncu_force_deep_summary.txt"
