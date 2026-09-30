#!/bin/bash
#SBATCH --job-name=cuda_cascade
#SBATCH -A hpc-prf-simvine
#SBATCH -p gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 01:00:00
#SBATCH --mem 40G
#SBATCH --output=cuda_bench_%j.out
#SBATCH --error=cuda_bench_%j.err

set -euo pipefail
module reset
module load system/CUDA/13.2.0
nvcc --version
nvidia-smi --query-gpu=name,compute_cap,memory.total --format=csv

# --- Build: double und float ----------------------------------------------
nvcc -O3 -arch=sm_80            -o cascade_cuda_bench_fp64 cascade_cuda_bench.cu
nvcc -O3 -arch=sm_80 -DUSE_FP32 -o cascade_cuda_bench_fp32 cascade_cuda_bench.cu

N2M=2000000

echo "############## (1) HEADLINE fp64 @ 2 Mio, 3000 Schritte, 5 Wdh ##############"
./cascade_cuda_bench_fp64 $N2M 256 3000 5

echo "############## (2) HEADLINE fp32 @ 2 Mio, 3000 Schritte, 5 Wdh ##############"
./cascade_cuda_bench_fp32 $N2M 256 3000 5

echo "############## (3) SAETTIGUNGS-SWEEP fp64 (400 Schritte, 3 Wdh) #############"
for N in 250000 500000 1000000 2000000 4000000 8000000 16000000; do
  ./cascade_cuda_bench_fp64 $N 256 400 3
done

echo "############## (4) BLOCKGROESSEN-SWEEP fp64 @ 2 Mio (400 Schritte) ##########"
for T in 128 256 512 1024; do
  ./cascade_cuda_bench_fp64 $N2M $T 400 3
done

echo "############## CSV-Zeilen gesammelt (fuer Plots) ##########################"
grep -h "# CSV " cuda_bench_${SLURM_JOB_ID}.out | sed 's/# CSV  //' || true
echo "fertig."
