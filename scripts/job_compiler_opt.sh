#!/bin/bash
#SBATCH -J cc_opt
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
# ============================================================================
#  job_compiler_opt.sh -- Compiler-Flag-Sweep der MPI-Version.
#  Fuer JEDE Flag-Kombination: kompilieren -> laufen -> Zeit + Energie-Drift
#  in EINE CSV schreiben. Der Energie-Drift zeigt, ob eine Optimierung die
#  Energieerhaltung verletzt (typisch bei -ffast-math / -Ofast!).
#
#  Aufruf:  sbatch job_compiler_opt.sh [params.ini] [ranks] [out.csv] [reps]
#  Default: params_scaling.ini  16  compiler_opt.csv  1
#
#  WICHTIG: log_every sollte n_steps teilen/gleich sein, damit Anfang UND Ende
#  in der Energie-CSV stehen (params_scaling.ini erfuellt das: 3000 == 3000).
# ============================================================================
set -u
SRC=mpi_cascade.c
PARAMS=${1:-params_scaling.ini}
RANKS=${2:-16}
OUT=${3:-compiler_opt.csv}
REPS=${4:-1}
MPICC=${MPICC:-mpicc}
MPIRUN=${MPIRUN:-mpirun}

# --- Flag-Kombinationen:  label | Flags --------------------------------------
#  znver3 = AMD EPYC 7763 (Noctua 2). Bei anderer CPU -march anpassen.
FLAGSETS=(
  "O0|-O0"
  "O1|-O1"
  "O2|-O2"
  "O3|-O3"
  "O3_march|-O3 -march=znver3"
  "O3_march_unroll|-O3 -march=znver3 -funroll-loops"
  "O3_march_lto|-O3 -march=znver3 -flto"
  "O3_native|-O3 -march=native"
  "O3_ffast|-O3 -march=znver3 -ffast-math"
  "Ofast|-Ofast -march=znver3"
)

# Praefix der Energie-CSV aus der params.ini lesen (Fallback: scale)
PREFIX=$(grep -E '^[[:space:]]*out_prefix' "$PARAMS" 2>/dev/null \
         | sed 's/;.*//' | awk -F= '{gsub(/[ \t]/,"",$2);print $2}')
PREFIX=${PREFIX:-scale}

echo "label,rep,flags,compile_ok,walltime_s,E_start,E_end,drift_pct,abs_drift_pct" > "$OUT"

for entry in "${FLAGSETS[@]}"; do
    label="${entry%%|*}"; flags="${entry#*|}"
    bin="./mpi_cc_${label}"
    echo ">>> $label : $flags"

    # ---- kompilieren ----
    if ! $MPICC $flags -o "$bin" "$SRC" -lm 2> "cc_${label}.log"; then
        echo "    KOMPILIERUNG FEHLGESCHLAGEN (siehe cc_${label}.log)"
        echo "${label},1,\"${flags}\",0,NA,NA,NA,NA,NA" >> "$OUT"
        continue
    fi

    # ---- REPS-mal laufen ----
    for ((rep=1; rep<=REPS; rep++)); do
        runlog=$($MPIRUN -np "$RANKS" "$bin" "$PARAMS" 2>&1)
        wt=$(printf '%s\n' "$runlog" | grep -oE 'max_walltime=[0-9.]+' | head -1 | cut -d= -f2)

        ecsv="${PREFIX}_energy.csv"
        if [ -f "$ecsv" ]; then
            es=$(awk -F, 'NR==2{print $7}' "$ecsv")          # E_total Anfang
            ee=$(awk -F, 'END{print $7}'   "$ecsv")          # E_total Ende
            dr=$(awk -v a="$es" -v b="$ee" 'BEGIN{if(a+0!=0)printf"%.6f",100*(b-a)/a;else print"NA"}')
            ad=$(awk -v d="$dr" 'BEGIN{if(d=="NA")print"NA";else printf"%.6f",(d<0?-d:d)}')
        else
            es=NA; ee=NA; dr=NA; ad=NA
        fi
        echo "    rep $rep: walltime=${wt:-NA}s  drift=${dr}%"
        echo "${label},${rep},\"${flags}\",1,${wt:-NA},${es},${ee},${dr},${ad}" >> "$OUT"
    done

    rm -f "$bin"     # Platz sparen; Logs cc_*.log bleiben fuer die Fehlersuche
done

echo "fertig -> $OUT"
echo "Auswertung:  python3 analyze_compiler_opt.py $OUT"
