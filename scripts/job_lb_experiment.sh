#!/bin/bash
#SBATCH --job-name=lb_experiment
#SBATCH --account=hpc-prf-simvine
#SBATCH --partition=gpu                 # PRUEFEN: GPU-Partition auf Noctua 2
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=lb_experiment_%j.out
#SBATCH --error=lb_experiment_%j.err
# ============================================================================
#  job_lb_experiment.sh -- EIN Experiment: hilft __launch_bounds__?
#
#  Vergleicht den Force-Kernel ohne/mit Register-Begrenzung. __launch_bounds__
#  (LB_T, MINBLK) zwingt ptxas, die Register je Thread so zu druecken, dass
#  MINBLK Bloecke gleichzeitig je SM laufen -> hoehere Occupancy. Gemessen wird,
#  ob das den DURCHSATZ hebt (oder ob die Latenz-/Bandbreitengrenze haelt).
#
#  Erwartung laut Diagnose: Register sinken, Occupancy steigt -- der Durchsatz
#  aendert sich aber evtl. kaum, weil der Kernel latenz-/bandbreiten-limitiert
#  ist. BEIDE Ausgaenge sind berichtfaehig.
#
#  Blockgroesse fest T=256 (== LB_T), sonst Launch-Fehler.
#  Submit:  sbatch job_lb_experiment.sh
# ============================================================================
set -u
module reset 2>/dev/null || true
module load CUDA || ml CUDA
nvcc --version | tail -1
nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader

SRC=cascade_cuda_bench_lb.cu
N=${1:-2000000}        # ~2 Mio Atome (Produktionsgroesse)
T=256                  # MUSS 256 sein (== LB_T)
NREP=${2:-5}
ARCH=sm_80
RES=lb_results.txt; : > "$RES"

build_run(){   # name  "extraflags"
    local name=$1 flags=$2
    echo ">> build $name ($flags)"
    nvcc -O3 -arch=$ARCH --ptxas-options=-v $flags -o "bench_$name" "$SRC" 2> "ptxas_$name.log"
    if [ ! -x "./bench_$name" ]; then echo "   BUILD FEHLER ($name)"; return; fi
    # Register je Thread fuer k_force aus der ptxas-Ausgabe
    local reg
    reg=$(awk '/entry function.*k_force/{f=1} f&&/Used [0-9]+ registers/{
              for(i=1;i<=NF;i++) if($i=="registers,"||$i=="registers") print $(i-1); f=0}' "ptxas_$name.log" | head -1)
    [ -z "$reg" ] && reg="?"
    ./bench_$name "$N" "$T" "$NREP" > "out_$name.txt" 2>&1
    local csv force matom
    csv=$(grep "^# CSV  fp64" "out_$name.txt" | tail -1)
    matom=$(echo "$csv" | awk -F, '{print $6}')
    force=$(echo "$csv" | awk -F, '{print $7}')
    echo "$name $reg $force $matom" >> "$RES"
    echo "   k_force: ${reg} Register/Thread | force=${force} us | ${matom} Matom/s"
}

build_run base ""
build_run lb4  "-DUSE_LB -DLB_MINBLK=4"
build_run lb6  "-DUSE_LB -DLB_MINBLK=6"
build_run lb8  "-DUSE_LB -DLB_MINBLK=8"

echo ""
python3 - "$RES" <<'PY'
import sys
rows=[l.split() for l in open(sys.argv[1]) if l.strip()]
base=next((r for r in rows if r[0]=="base"), None)
bforce=float(base[2]) if base and base[2] not in("","?") else None
bmat=float(base[3]) if base else None
print("="*70)
print("  __launch_bounds__-EXPERIMENT (k_force, fp64, A100)")
print("="*70)
print(f"{'Variante':<8}{'Reg/Thread':>11}{'force[us]':>11}{'Matom/s':>11}{'Speedup':>9}")
print("-"*70)
for r in rows:
    name,reg,force,mat=r[0],r[1],r[2],r[3]
    sp=""
    try:
        if bforce: sp=f"{bforce/float(force):.3f}x"
    except: pass
    print(f"{name:<8}{reg:>11}{force:>11}{mat:>11}{sp:>9}")
print("-"*70)
print("Lies es so: sinken die Register und steigt die Occupancy, aber der")
print("Durchsatz bleibt ~gleich -> die Latenz-/Bandbreitengrenze haelt (= Ergebnis).")
print("Steigt der Durchsatz -> __launch_bounds__ hat geholfen.")
print("="*70)
PY

# --- optional: echte erreichte Occupancy via ncu (base vs lb6) --------------
if command -v ncu >/dev/null 2>&1; then
    echo ""; echo "### ncu: erreichte Occupancy + Register (base vs lb6) ###"
    for v in base lb6; do
        [ -x "./bench_$v" ] || continue
        echo "--- $v ---"
        ncu --launch-count 1 --launch-skip 50 --kernel-name regex:k_force \
            --metrics sm__warps_active.avg.pct_of_peak_sustained_active,launch__registers_per_thread \
            ./bench_$v "$N" "$T" 1 2>/dev/null \
          | grep -E "warps_active|registers_per_thread" || echo "  (ncu lieferte nichts)"
    done
else
    echo "ncu nicht verfuegbar -- Register/Durchsatz oben reichen fuer die Aussage."
fi
