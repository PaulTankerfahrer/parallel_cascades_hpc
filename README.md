# Kollisionskaskaden-Simulator — HPC-Endprojekt

Ein konfigurierbarer, paralleler 2D-Masse-Feder-Simulator für Kollisions-
kaskaden (Strahlenschaden). Drei Implementierungen teilen sich dieselbe
Konfigurationsdatei und dasselbe Modell: **seriell**, **MPI** und **CUDA**.

[![Simulation CI](https://forgejo.paultankerfahrer.org/PaulTankerfahrer/parallel_cascades_hpc/actions/workflows/simulation.yml/badge.svg)](https://forgejo.paultankerfahrer.org/PaulTankerfahrer/parallel_cascades_hpc/actions)

---

## Projektstruktur

```
project/
├── src/                             # Quellcode
│   ├── cascade_serial.c             # Serielle Referenzimplementierung
│   ├── mpi_cascade.c                # MPI-Version (1D-Domänenzerlegung)
│   ├── cascade_cuda.cu              # CUDA-Version (Gather, ein Thread/Atom)
│   ├── cascade_cuda_bench.cu        # CUDA-Benchmark: fp32/fp64, N-Sweep, Kernel-Breakdown
│   ├── cascade_cuda_bench_opt.cu    # CUDA-Bench + __launch_bounds__ & CUDA Graphs (4 Varianten)
│   ├── params.ini                   # Konfigurationsdatei (alle Versionen)
│   ├── params_multipka.ini          # Konfiguration mit mehreren PKAs
│   └── params_valgrind.ini          # Kleines Gitter (80×80) für die Valgrind-Analyse
│
├── scripts/                         # Slurm-Job-Skripte
│   ├── build_job.sh
│   ├── job_ensemble.sh              # 1000-Kaskaden-Ensemble
│   ├── job_strong_scaling_mpi.sh
│   ├── job_scaling_rep.sh           # Scaling mit Wiederholungen
│   ├── job_weak.sh
│   ├── job_compiler_opt.sh
│   ├── job_bench.sh                 # CUDA-Benchmark
│   ├── job_cuda.sh                  # CUDA-Lauf
│   ├── job_ncu_deep.sh              # ncu-Profiling
│   └── job_lb_experiment.sh         # __launch_bounds__ (base/lb4/lb6/lb8) — veraltet, siehe unten
│
├── scripts/python/                  # Python-Auswertungs- und Plot-Skripte
│   ├── plot_energy.py               # Energiebilanz + Drift-Check
│   ├── analyze.py                   # Schadenskarte + Cluster (seriell)
│   ├── preview_xyz.py               # GIF-Vorschau aus .xyz-Dump
│   ├── run_one.py                   # randomisierte Parameter für Ensemble-Job
│   ├── aggregate_reps.py            # Median über Wiederholungen
│   ├── analyze_scaling.py           # Strong-Scaling: Speedup, Effizienz, Phasen
│   ├── analyze_weak.py              # Weak-Scaling: Wandzeit + Effizienz
│   ├── analyze_compiler_opt.py      # Compiler-Flag-Sweep
│   ├── analyze_opt_compare.py       # Blocking vs. non-blocking + resort
│   ├── analyze_ensemble.py          # Union-Find + MLE-Fit über alle Läufe
│   ├── overlay_ccdf.py              # Overlay CCDF: mit vs. ohne Healing
│   ├── plot_damage_map.py           # Schadenskarte einer Einzelkaskade (heal vs. no-heal)
│   ├── plot_tts.py                  # Time-to-Solution: CPU vs. T4 vs. A100
│   ├── plot_ncu.py                  # ncu-Profiling: Kernel-Breakdown, SpeedOfLight, Okkupanz
│   ├── plot_launch_bounds.py        # __launch_bounds__: Occupancy vs. Durchsatz
│   ├── plot_cuda_opt.py             # __launch_bounds__ + CUDA Graphs: 4-Varianten-Vergleich
│   ├── plot_gpu_saturation.py       # Durchsatz [Matom/s] vs. Problemgröße N
│   ├── plot_cache_argument.py       # Datensatz/Rank + Effizienz (Cache-Schwelle)
│   └── plot_ensemble_scaling.py     # Ensemble-Machbarkeit: MPI-Scaling + Energie-Kalibrierung
│
├── results/                         # Messdaten + erzeugte Plots
│   ├── KATALOG.md                   # Vollständige Bestandsaufnahme aller Messdateien
│   ├── mpi_strong_scaling_1node.csv
│   ├── mpi_strong_scaling_2node.csv
│   ├── mpi_weak_scaling_2nodes.csv
│   ├── mpi_compiler_opt.csv
│   ├── opt_compare_opt_off_raw.csv
│   ├── opt_compare_opt_on_raw.csv
│   ├── opt_compare_n128.csv
│   ├── ncu_report.csv               # Nsight Compute A100 (detailliert)
│   ├── ncu_report2.csv
│   ├── ncu_force_deep_summary.txt
│   ├── launch_bounds_results.csv    # base/lb4/lb6/lb8: Register, Occupancy, Durchsatz
│   ├── cuda_opt_results.csv         # base/lb/graph/lbgraph: Occupancy, Durchsatz
│   ├── cuda_bench_33375548.out      # A100 fp64+fp32 Benchmark + Sweep (Slurm-Log)
│   ├── cuda_cascade_33377028.out    # A100 voller Kaskaden-Lauf (Slurm-Log)
│   ├── cuda_ncu_33375549.out        # ncu-Profiling-Job (Slurm-Log)
│   ├── tesdla_t4_cuda_run.txt       # Tesla T4 Lauf + ncu-Output
│   ├── ensemble_run_0999_state_final.csv          # Einzellauf mit Healing (Schadenskarte)
│   ├── ensemble_no_heal_run_0999_state_final.csv  # Einzellauf ohne Healing (Schadenskarte)
│   ├── valgrind/                    # Valgrind-Memcheck-Logs (seriell + MPI, 2 Ränge)
│   │   ├── valgrind_serial.txt
│   │   ├── valgrind_mpi_rank_0.txt
│   │   └── valgrind_mpi_rank_1.txt
│   └── images/                      # alle Plots (PNG) für report.tex
│       ├── — Modell —
│       ├── M1_model_schema.png
│       ├── M2_force_law.png
│       ├── M3_cell_list.png
│       ├── M4_pka_velocity.png
│       ├── — Korrektheit —
│       ├── vid_S1_energy.png
│       ├── vid_S2_damage.png
│       ├── vid_S6_timeseries.png
│       ├── — MPI-Skalierung —
│       ├── H1_strong_scaling.png
│       ├── H2_phase_breakdown.png
│       ├── H3_compute_vs_comm.png
│       ├── H5_weak_scaling.png
│       ├── H6_opt_compare.png
│       ├── H7_cache_argument.png
│       ├── H8_node_comparison.png
│       ├── H9_compiler_opt.png
│       ├── — CUDA / GPU —
│       ├── tts.png
│       ├── ncu_profiling.png
│       ├── launch_bounds_plot.png
│       ├── cuda_opt_plot.png
│       ├── gpu_saturation.png
│       ├── cache_argument.png
│       ├── — Power Law —
│       ├── E1_pooling_progression.png
│       ├── E2_powerlaw_overlay.png
│       ├── E3_powerlaw_pdf_ccdf.png
│       ├── ccdf_overlay.png
│       ├── — Ensemble-Machbarkeit —
│       └── ensemble_scaling.png
│
├── report.tex                       # LaTeX-Bericht (2-teilig: HPC + Power Law)
├── report.pdf                       # kompiliertes PDF
└── Präsentation.html                # eigenständige HTML/CSS/JS-Präsentation (im Browser öffnen)
```

`results/images/` enthält teils zwei Versionen desselben Plots: eine
unpräfigierte Rohversion direkt aus dem Python-Skript (z. B.
`strong_scaling.png`) und eine für den Bericht umbenannte/neu erzeugte
Version mit `H`/`E`/`M`/`vid`-Präfix (z. B. `H1_strong_scaling.png`). Details
zur Herkunft jeder einzelnen Datei in `results/`: siehe `results/KATALOG.md`.

---

## Allgemein

Dies ist das Endprojekt für den Wahlpflichtkurs HPC. Es verfolgt zwei Ziele:

- **Primär (HPC):** ein konfigurierbarer, paralleler Kaskaden-Simulator.
  Nachweisbare **MPI-Strong-Scaling-Kurve** und **Single-GPU-CUDA-Port** mit
  vollständigem CPU-vs-GPU-Vergleich (A100: 731× gegenüber 1 EPYC-Kern),
  ergänzt um eine ncu-gestützte Bottleneck-Diagnose und drei gezielte
  GPU-Optimierungsexperimente (Occupancy, Launch-Overhead, Präzision), die
  den Speicherbandbreiten-Engpass des `k_force`-Kernels experimentell
  bestätigen — sowie eine Valgrind-Speicherkorrektheitsprüfung beider
  CPU-Implementierungen.
- **Sekundär (Physik):** Aus einem **Ensemble** von 1000 Kaskaden die
  Defektcluster-Größenverteilung bestimmen und prüfen, ob ein **Power Law**
  F(n) ∝ n⁻ˢ emergiert. Ergebnis: S = 1,38–1,41, robust gegen Healing.

---

## Konfiguration: `src/params.ini`

Alle drei Implementierungen lesen dieselbe Datei. Keine Neukompilierung bei
Parameteränderungen nötig. Die Datei selbst ist inzwischen ausführlich
kommentiert (Einheitensystem, Bedeutung jedes Feldes); hier nur die
aktuellen Standardwerte (Strong-Scaling-Benchmarkkonfiguration, ~2 Mio. Atome):

```ini
[model]
model         = A          ; A = harmonische Feder + Potenzgesetz-Abstoßung
healing       = false      ; true = gerissene Bindungen können rekombinieren
healing_dist  = 1.10
healing_vrel  = 0.5

[grid]
NX            = 1414       ; 1414 × 1414 ≈ 2,0 Mio. Atome
NY            = 1414
L0            = 1.0
lattice       = triangular

[potential]
K_SPRING      = 100.0
MAX_STRETCH   = 1.15       ; Bruch bei r > L0 × 1,15
RCUT          = 0.9
K_REP         = 400.0
REP_N         = 12.0

[pka]
pka_x         = 0.5
pka_y         = 0.5
pka_energy    = 5000.0
pka_angle     = 15.0
pka_mass      = 1.0
n_pka         = 1
timing        = true
halo_nonblocking = false   ; NUR MPI: Isend/Irecv statt Sendrecv
resort_every  = 0          ; >0: alle N Schritte cell-linear umsortieren

[dynamics]
dt            = 2e-4
n_steps       = 3000
damping       = 0.0
v_thresh      = 1e30
absorb_border = 0

[output]
log_every     = 3000
dump_every    = 0          ; >0: XYZ-Frames für OVITO
enable_vtk    = false
out_prefix    = scale

[ensemble]
seed          = 3
```

Für die Valgrind-Analyse gibt es ein deutlich kleineres Gitter
(`src/params_valgrind.ini`, 80×80 = 6400 Atome, 300 Schritte) — Memcheck
instrumentiert jeden Speicherzugriff und läuft 20–50× langsamer, das
2-Mio.-Atome-Gitter wäre damit stundenlang unterwegs.

**Gemessene Phasen (MPI-Version, `timing=true`):**

| Phase | Typ | Anteil bei 1 Kern |
|---|---|---|
| `force` | Rechnung | 89 % |
| `cells` | Rechnung | 6 % |
| `hash` | Rechnung | 5 % |
| `halo` | Kommunikation | < 1 % |
| `migrate` | Kommunikation | < 1 % |
| `energy` | Diagnose | < 1 % |

---

## Build-Anleitungen

### Serielle Version

```bash
gcc -O3 -march=native -o cascade_serial src/cascade_serial.c -lm
./cascade_serial src/params.ini
```

### MPI-Version

```bash
mpicc -O3 -march=znver3 -o mpi_cascade src/mpi_cascade.c -lm
mpirun -np 4 ./mpi_cascade src/params.ini
```

### CUDA-Version

```bash
# A100 (sm_80), fp64
nvcc -O3 -arch=sm_80 -o cascade_cuda src/cascade_cuda.cu

# Benchmark-Binary (Sweep + Block-Größen-Test, fp64/fp32 per Flag)
nvcc -O3 -arch=sm_80 -o cascade_cuda_bench src/cascade_cuda_bench.cu
nvcc -O3 -arch=sm_80 -DUSE_FP32 -o cascade_cuda_bench_fp32 src/cascade_cuda_bench.cu

# GPU-Optimierungsexperiment: __launch_bounds__ + CUDA Graphs, 4 unabhängige
# und kombinierbare Varianten (Details im Kopfkommentar der Datei):
nvcc -O3 -arch=sm_80                      -o bench_base    src/cascade_cuda_bench_opt.cu
nvcc -O3 -arch=sm_80 -DUSE_LB             -o bench_lb      src/cascade_cuda_bench_opt.cu
nvcc -O3 -arch=sm_80 -DUSE_GRAPH          -o bench_graph   src/cascade_cuda_bench_opt.cu
nvcc -O3 -arch=sm_80 -DUSE_LB -DUSE_GRAPH -o bench_lbgraph src/cascade_cuda_bench_opt.cu
# Aufruf wie cascade_cuda_bench.cu: ./bench_<variante> [N] [blocksize=256] [n_steps] [n_rep]
# -> Ergebnisse nach results/cuda_opt_results.csv, Plot: plot_cuda_opt.py

# Tesla T4 (sm_75)
nvcc -O3 -arch=sm_75 -o cascade_cuda src/cascade_cuda.cu
```

### Speicherkorrektheit: Valgrind

```bash
# Serielle Version (mit gcc gebaut, damit keine MPI-Bibliothek stört):
gcc -O1 -g -march=native -o cascade_serial src/cascade_serial.c -lm
valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes \
         ./cascade_serial src/params_valgrind.ini

# MPI-Version (2 Ränge, ein Logfile pro Rang):
mpicc -O1 -g -march=native -o mpi_cascade src/mpi_cascade.c -lm
mpirun --oversubscribe -np 2 valgrind --leak-check=full \
       --show-leak-kinds=definite,indirect --track-origins=yes \
       --log-file=vg_rank_%q{OMPI_COMM_WORLD_RANK}.txt \
       ./mpi_cascade src/params_valgrind.ini
```

Ergebnis (siehe `results/valgrind/` und Abschnitt 2.3 im Bericht): 0 ungültige
Zugriffe und 0 uninitialisierte Werte in beiden Versionen; die einzigen
Funde sind ≈1 MB harmlose `still reachable`/`definitely lost`-Allokationen
(App) sowie eine bekannte OpenMPI-Fehlalarmklasse aus `MPI_Init`.

---

## SLURM-Jobs (Noctua 2)

Alle Messungen laufen über Slurm-Batch-Jobs in `scripts/`. Die `#SBATCH`-Header
sind im jeweiligen Skript gesetzt; hier die Submit-Befehle. **Vor dem ersten
Lauf** in den Skripten `-A <projekt>` / `--account` (Compute-Kontingent) und die
Partition (`-p normal` für CPU, `--partition=gpu` für A100) an den eigenen
Cluster-Account anpassen.

```bash
cd scripts/

# ---- Setup -------------------------------------------------------------
sbatch build_job.sh                       # serielle Version bauen (1 Kern, 5 min)

# ---- MPI-Skalierung (CPU, normal-Partition) ----------------------------
sbatch job_strong_scaling_mpi.sh          # Strong Scaling, 2 Knoten / 256 Ranks
sbatch job_scaling_rep.sh                 # Strong Scaling mit Wiederholungen (Median)
sbatch job_weak.sh                        # Weak Scaling, 2 Knoten / 256 Ranks
sbatch job_compiler_opt.sh                # Compiler-Flag-Sweep (-O0 … -O3 …)

# ---- Ensemble (CPU, 1 Knoten, xargs -P 128) ----------------------------
sbatch job_ensemble.sh                     # 1000 Kaskaden (Default)
sbatch job_ensemble.sh 500                 # oder weniger zum Testen

# ---- GPU (A100, gpu-Partition, --gres=gpu:a100:1) ----------------------
sbatch job_bench.sh                        # CUDA-Benchmark fp64/fp32 + Sättigungs-Sweep
sbatch job_cuda.sh                         # voller CUDA-Kaskaden-Lauf
sbatch job_ncu_deep.sh                     # Nsight-Compute-Profiling des force-Kernels
sbatch job_lb_experiment.sh                # __launch_bounds__-Experiment (base/lb4/lb6/lb8)

# ---- Monitoring --------------------------------------------------------
squeue --me                                # eigene Jobs anzeigen
sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,State   # Ressourcen nach Lauf
```

| Job-Skript | Ressourcen (`#SBATCH`) | Erzeugt |
|---|---|---|
| `build_job.sh` | `-N1 -n1 -p normal -t 5min` | `cascade_serial` |
| `job_strong_scaling_mpi.sh` | `-N2 -n256 -p normal -t 1h` | `scaling.csv` |
| `job_scaling_rep.sh` | `-N2 -n256 -p normal` | `*_raw.csv` (mehrere Reps) |
| `job_weak.sh` | `-N2 -n256 -p normal -t 1h` | `weak.csv` |
| `job_compiler_opt.sh` | `-N1 -n… -p normal` | `mpi_compiler_opt.csv` |
| `job_ensemble.sh` | `-N1 -n128 -p normal -t 30min` | `ensemble/run_*_state_final.csv` |
| `job_bench.sh` | `--gres=gpu:a100:1` | `cuda_bench_*.out` |
| `job_cuda.sh` | `--gres=gpu:a100:1` | `cuda_cascade_*.out` |
| `job_ncu_deep.sh` | `--gres=gpu:a100:1` | `ncu_report*.csv` |
| `job_lb_experiment.sh` | `--partition=gpu --gres=gpu:a100:1 -t 30min` | `lb_results.txt`, `ptxas_*.log` |

> **Hinweis zu `job_lb_experiment.sh`:** Das Skript baut
> `cascade_cuda_bench_lb.cu`, das im Zuge der GPU-Optimierungsversuche durch
> `cascade_cuda_bench_opt.cu` ersetzt wurde (dasselbe `__launch_bounds__`-
> Experiment, zusätzlich CUDA Graphs) — die alte Quelldatei existiert nicht
> mehr, das Skript ist damit aktuell nicht lauffähig. Die vier Varianten
> (`base`/`lb`/`graph`/`lbgraph`) für `results/cuda_opt_results.csv` wurden
> manuell mit den Befehlen aus dem Build-Abschnitt oben gebaut und
> ausgeführt; ein passendes Slurm-Skript dafür fehlt noch.

Die meisten Jobs nehmen optionale Positionsargumente (Parameterdatei, Ausgabe-CSV,
Rankzahl, Wiederholungen) — siehe Kopfkommentar im jeweiligen Skript, z. B.
`sbatch job_compiler_opt.sh params.ini 128 out.csv 3`.

---

## Python-Skripte ausführen

Alle Skripte werden aus `results/images/` aufgerufen, damit die PNGs dort landen:

```bash
cd results/images/
P=../../scripts/python   # Kürzel
R=..                     # Kürzel für results/
```

### MPI-Skalierung

```bash
python3 $P/analyze_scaling.py \
    $R/mpi_strong_scaling_1node.csv $R/mpi_strong_scaling_2node.csv
# → strong_scaling.png

python3 $P/analyze_weak.py $R/mpi_weak_scaling_2nodes.csv
# → mpi_weak_scaling_2nodes_weak.png

python3 $P/analyze_compiler_opt.py $R/mpi_compiler_opt.csv
# → compiler_opt.png

python3 $P/analyze_opt_compare.py \
    $R/opt_compare_opt_off_raw.csv "blocking (Baseline)" \
    $R/opt_compare_opt_on_raw.csv  "non-blocking + resort"
# → opt_compare.png

python3 $P/plot_cache_argument.py \
    $R/mpi_strong_scaling_1node.csv $R/mpi_strong_scaling_2node.csv
# → cache_argument.png
```

### CUDA / GPU

```bash
python3 $P/plot_tts.py             # → tts.png
python3 $P/plot_ncu.py             # → ncu_profiling.png
python3 $P/plot_gpu_saturation.py  # → gpu_saturation.png

# __launch_bounds__-Experiment (CSV aus job_lb_experiment.sh + ncu-Occupancy):
python3 $P/plot_launch_bounds.py $R/launch_bounds_results.csv
# → launch_bounds_plot.png

# GPU-Optimierung: __launch_bounds__ + CUDA Graphs, 4 Varianten (siehe Build-Anleitung):
python3 $P/plot_cuda_opt.py $R/cuda_opt_results.csv
# → cuda_opt_plot.png
```

### Power Law / Ensemble

```bash
python3 $P/analyze_ensemble.py $R/ensemble
python3 $P/analyze_ensemble.py $R/ensemble_no_heal
# → ensemble/cluster_sizes.csv + ensemble_powerlaw.png (je Ordner)
# Hinweis: die vollen 1000-Läufe-Rohdaten (ensemble/, ensemble_no_heal/) werden
# von job_ensemble.sh erzeugt, sind aber wegen der Datenmenge nicht Teil des
# Repos — nur je ein repräsentativer Einzellauf (run_0999) ist eingecheckt.

python3 $P/overlay_ccdf.py \
    $R/ensemble/cluster_sizes.csv         "mit Healing" \
    $R/ensemble_no_heal/cluster_sizes.csv "ohne Healing"
# → ccdf_overlay.png

# Schadenskarte einer Einzelkaskade (heal vs. no-heal, rekonstruiert Bindungen
# geometrisch aus den beiden eingecheckten End-CSVs):
python3 $P/plot_damage_map.py \
    $R/ensemble_no_heal_run_0999_state_final.csv \
    $R/ensemble_run_0999_state_final.csv
# → vid_S2_damage.png

# Machbarkeitsstudie skalierter Ensemble-Lauf (Daten hardcoded):
python3 $P/plot_ensemble_scaling.py   # → ensemble_scaling.png
```

### Einzel-Lauf

```bash
python3 $P/plot_energy.py $R/run_energy.csv   # Energiebilanz
python3 $P/analyze.py     $R/run              # Schadenskarte + Cluster
```

---

## Skript-Referenz

| Skript | Eingabe | Ausgabe |
|---|---|---|
| `analyze_scaling.py` | `mpi_strong_scaling_*.csv` | `strong_scaling.png` |
| `analyze_weak.py` | `mpi_weak_scaling_*.csv` | `*_weak.png` |
| `analyze_compiler_opt.py` | `mpi_compiler_opt.csv` | `compiler_opt.png` |
| `analyze_opt_compare.py` | `opt_compare_opt_{off,on}_raw.csv` | `opt_compare.png` |
| `plot_cache_argument.py` | `mpi_strong_scaling_*.csv` | `cache_argument.png` |
| `plot_tts.py` | *(hardcoded)* | `tts.png` |
| `plot_ncu.py` | *(hardcoded)* | `ncu_profiling.png` |
| `plot_launch_bounds.py` | `launch_bounds_results.csv` | `launch_bounds_plot.png` |
| `plot_cuda_opt.py` | `cuda_opt_results.csv` | `cuda_opt_plot.png` |
| `plot_gpu_saturation.py` | *(hardcoded)* | `gpu_saturation.png` |
| `plot_ensemble_scaling.py` | *(hardcoded)* | `ensemble_scaling.png` |
| `analyze_ensemble.py` | `ensemble/run_*_state_final.csv` | `cluster_sizes.csv` + PNG |
| `plot_damage_map.py` | `ensemble*_run_0999_state_final.csv` | `vid_S2_damage.png` |
| `overlay_ccdf.py` | `*/cluster_sizes.csv` | `ccdf_overlay.png` |
| `aggregate_reps.py` | `*_raw.csv` | aggregierte CSV |
| `plot_energy.py` | `*_energy.csv` | Energiebilanz-Plot |
| `analyze.py` | Lauf-Präfix | Schadenskarte + Cluster |

---

## Wichtige Messergebnisse

### GPU-Beschleunigung (N = 2 Mio Atome, 3000 Schritte)

| System | Laufzeit | Speedup vs. 1 Kern |
|---|---|---|
| CPU 1 Kern (EPYC 7763) | 1152,9 s | 1× |
| CPU 64 Kerne | 22,5 s | 51× |
| Tesla T4 (fp64) | 19,1 s | 60× |
| **A100 (fp64)** | **1,58 s** | **731×** |
| A100 (fp32) | 1,00 s | 1154× |

*A100 läuft mit 45 % Compute- und 37 % Speicher-Auslastung (latenz-limitiert
durch 50 % Okkupanz bei 62 Registern/Thread, fp64).*

### GPU-Optimierungsversuche (base/lb/graph/lbgraph, N ≈ 2 Mio Atome, fp64)

| Variante | Reg/Thread | Occupancy | Durchsatz |
|---|---|---|---|
| `base` | 62 | 45 % | 3808 Matom/s |
| `lb` (`__launch_bounds__`) | 40 | 66 % | 3760 Matom/s |
| `graph` (CUDA Graphs) | 62 | 45 % | 3847 Matom/s |
| `lb+graph` | 40 | 66 % | 3830 Matom/s |

*Occupancy fast verdoppelt (45 %→66 %), Durchsatz bleibt im Messrauschen
(±1 %) — bestätigt experimentell, dass `k_force` speicherlatenz-/
bandbreitengebunden ist, nicht occupancy- oder launch-overhead-limitiert.
Ein vollwertiges Negativergebnis, siehe Abschnitt 11.4 im Bericht.*

### Power Law (1000 Kaskaden, MLE-Fit mit x\_min = 4)

| Variante | Cluster | Exponent S |
|---|---|---|
| Mit Healing | 2391 | 1,41 |
| Ohne Healing | 2335 | 1,38 |
| Sand et al. 2013 (3D-Wolfram) | — | 1,63 ± 0,07 |

---

## Bericht

- **`report.tex`** / **`report.pdf`** — zweiteiliger LaTeX-Bericht:
  - Teil I: HPC (Modell → MPI-Parallelisierung + Valgrind-Speicherkorrektheit
    → Compiler-Optimierungen → Strong/Weak Scaling → Cache-Argument →
    1-/2-Knoten-Vergleich → CUDA-Portierung → GPU-Sättigung → ncu-Profiling
    → `__launch_bounds__`-Experiment → CUDA Graphs & Präzision → Fazit)
  - Teil II: Power Law (Theorie → Ensemble-Setup → MLE-Fit-Methodik →
    Ergebnisse → Diskussion → Ensemble-Machbarkeitsstudie → Fazit)
  - Bau: `latexmk -pdf report.tex` (oder `pdflatex report.tex && pdflatex report.tex`)
- **`Präsentation.html`** — eigenständige HTML/CSS/JS-Präsentation, keine
  externen Abhängigkeiten, Tastatur-Navigation (Pfeiltasten/Space/PageDown);
  einfach im Browser öffnen.

---

## Continuous Integration

Bei jedem Push laufen über **Forgejo Actions** (`.forgejo/workflows/`) billige,
cluster-freie Checks — die schwere Rechnung (Scaling, GPU, 1000er-Ensembles)
bleibt dem Cluster vorbehalten, nicht der CI.

### Simulation CI (`simulation.yml`) — läuft bei jedem Push

- **Physik-Smoke-Test** (`scripts/ci/build_and_test.sh`): baut serielle und
  MPI-Version, dann auf einem kleinen Gitter (`src/params_valgrind.ini`):
  1. **Energieerhaltung** — der relative Energiedrift muss `|Δ| < 0,5 %`
     bleiben (real ~0,04 %). Fängt physikalisch kaputte Änderungen an der
     Integration/Kraftberechnung, bevor sie Ergebnisse verfälschen.
  2. **Determinismus** — zweimal derselbe Seed muss dieselbe Zahl gerissener
     Bindungen (`broken_bonds`) liefern.
  3. **Sanity** — es ist überhaupt eine Kaskade passiert (`broken_bonds > 0`).
  4. **MPI-Smoke** — 2 Ränge laufen ohne Absturz und geben `RESULT` aus.
- **Valgrind-Speichercheck** (`scripts/ci/valgrind_check.sh`): baut mit `-g`
  und lässt die serielle Version streng unter Valgrind laufen
  (`--error-exitcode=1`, `--leak-check=full`, definitive Leaks = rot). Der
  vollständige Report wird als Artefakt **`valgrind-report`** hochgeladen (auch
  bei gefundenen Fehlern, `if: always()`) und lässt sich vom jeweiligen
  CI-Lauf herunterladen.

Beide Skripte sind auch **lokal** ausführbar:
`bash scripts/ci/build_and_test.sh` bzw. `bash scripts/ci/valgrind_check.sh`.
(Hinweis: Auf sehr aktuellem glibc — Rolling-Release — kann Valgrind lokal
beim Start scheitern; auf dem CI-Runner mit stabilem glibc läuft es.)

### Dokumente (`documents.yml`) — nur auf Anforderung

Baut `report.pdf` und die Präsentation mit TeXLive und lädt beide PDFs als
Artefakt **`pdfs`** hoch. Läuft **nicht** bei jedem Push, sondern nur, wenn die
Commit-Message das Schlüsselwort **`document`** enthält — oder manuell über
„Run workflow" (`workflow_dispatch`).

> Beide Jobs laufen im Container `node:20-bookworm`: Der podman-Runner führt
> JS-Actions (`checkout`, `upload-artifact`) mit dem node *im Container* aus,
> darum ein Image mit node an Bord (plus git/gcc); OpenMPI, Valgrind bzw. eine
> TeXLive-Teilmenge werden per `apt` nachgezogen.

---

## Über das Projekt

Dieses Projekt entstand als Kursprojekt im Modul *High Performance Computing*
an der **TH Lübeck**, betreut von **Dr. Maurice Maurer**. Die Messungen wurden
auf dem Cluster *Noctua 2* (PC2 Paderborn, EPYC 7763, A100-SXM4-40GB)
durchgeführt.

## Lizenz

- **Code** (`src/`, `scripts/`) steht unter der **MIT-Lizenz** — siehe
  [`LICENSE`](LICENSE).
- **Bericht und Präsentation** (`report.tex`, `report.pdf`,
  `Präsentation.html`) stehen unter
  [**CC BY 4.0**](https://creativecommons.org/licenses/by/4.0/) — freie
  Weiterverwendung bei Namensnennung.

© 2026 Paul Kamm
