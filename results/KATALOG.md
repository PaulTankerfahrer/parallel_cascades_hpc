# Katalog der Ergebnisse (`results/`)

Vollständige Auflistung aller in `results/` liegenden Mess- und Ausgabedateien.
Reine Bestandsaufnahme (was wurde gemessen, womit, womit), keine Analyse.

---

## 1. Messumgebung

Aus den Slurm-`.out`-Dateien und `scripts/*.sh` rekonstruiert:

| Komponente     | Wert                                                      |
|----------------|-----------------------------------------------------------|
| System         | Noctua 2 (PC2 Paderborn)                                  |
| CPU-Knoten     | AMD EPYC 7763, 128 Kerne/Knoten (`-march=znver3`)         |
| GPU-Knoten     | NVIDIA A100-SXM4-40GB, CC 8.0, 108 SMs, 40 GB             |
| CUDA           | 13.2 (nvcc 13.2.51)                                       |
| MPI            | OpenMPI 5.0.8 / GCC 14.3.0                                |
| Python         | 3.13.5 (für `run_one.py`, Auswerteskripte)               |
| MPI-Code       | `mpi_cascade.c`, 1D-Zerlegung, `-O3 -march=znver3`        |
| CUDA-Code      | `cascade_cuda_bench.cu`, `-O3 -arch=sm_80` (fp32 via `-DUSE_FP32`) |

## 2. Bedeutung der MPI-CSV-Spalten

`ranks, walltime, force, cells, hash, halo, migrate, energy, comm`

- `walltime` = maximale Wandzeit über alle Ranks (bestimmt die Laufzeit).
- `force, cells, hash, halo, migrate, energy` = phasenaufgelöste Zeiten
  (Maximum über alle Ranks), ausgegeben in `src/mpi_cascade.c` (Abschnitt
  `# DATA`, Zeilen ~606–623).
- `comm` = `halo + migrate` (Kommunikationsanteil).

---

## 3. HPC-Performance-Messungen (MPI)

### 3.1 MPI Strong Scaling – 1 Knoten
- **Datei:** `mpi_strong_scaling_1node.csv`
- **Skript:** `job_strong_scaling_mpi.sh` (1-Knoten-Teil), `-N 1`.
- **Methode:** Feste Problemgröße, wachsende Rankzahl.
- **Parameter:** `params.ini` – NX=NY=1414 (~2,0 Mio. Atome), Modell A,
  n_steps=3000, dt=2e-4, 1 PKA, seed=3, `halo_nonblocking=false`,
  `resort_every=0`.
- **Build:** `mpicc -O3 -march=znver3 mpi_cascade.c -lm`.
- **Ranks:** 1, 2, 4, 8, 16, 32, 64, 128 (je 1 Lauf).
- **Wertebereich walltime:** 1152,9 s (1 Rank) → 6,7 s (128 Ranks).

### 3.2 MPI Strong Scaling – 2 Knoten (Inter-Node)
- **Datei:** `mpi_strong_scaling_2node.csv`
- **Skript:** `job_strong_scaling_mpi.sh` (2-Knoten-Teil), `-N 2`,
  `--distribution=block`.
- **Methode:** Gleiche Problemgröße wie 3.1, über 2 Knoten verteilt zum
  Messen der Inter-Node-Kommunikation.
- **Ranks:** 64, 128, 256 (je 1 Lauf).
- **Wertebereich walltime:** 22,3 s (64) → 3,6 s (256).

### 3.3 MPI Weak Scaling – 2 Knoten
- **Datei:** `mpi_weak_scaling_2nodes.csv`
- **Skript:** `job_weak.sh`, `-N 2`.
- **Methode:** Konstante Last pro Rank: NX_fix=1024, 128 Zeilen/Rank
  (= 131 072 Atome/Rank), NY wächst mit der Rankzahl. Pro Rankzahl wird eine
  eigene `.ini` generiert (gleiche Physik-Parameter wie `params.ini`, seed=3).
- **Build:** `mpicc -O3 -march=znver3 -mprefer-vector-width=256`.
- **Ranks:** 8, 16, 32, 64, 128, 256 (je 1 Lauf).
- **Wertebereich walltime:** 37,9 s (8 Ranks) → 161,8 s (256 Ranks).

### 3.4 MPI Compiler-Flag-Sweep
- **Datei:** `mpi_compiler_opt.csv`
- **Skript:** `job_compiler_opt.sh`, `-N 1`, `-n 16`.
- **Methode:** 10 Flag-Kombinationen, jeweils neu kompiliert, 1 Lauf mit
  `params.ini` (NX=NY=1414). Pro Flag: Wandzeit + Energie-Drift
  (E_start/E_end aus `scale_energy.csv`) zur Erkennung von
  Energieneverletzungen (z. B. durch `-ffast-math` / `-Ofast`).
- **Spalten:** `label, rep, flags, compile_ok, walltime_s, E_start, E_end,
  drift_pct, abs_drift_pct`.
- **Flag-Sets:** O0, O1, O2, O3, O3_march, O3_march_unroll, O3_march_lto,
  O3_native, O3_ffast, Ofast.
- **Hinweis:** Alle 10 Sets kompilierten erfolgreich (`compile_ok=1`);
  `drift_pct` ist über alle Flags identisch (0,215 %).

### 3.5 MPI Opt-Vergleich (Off vs. On) – Strong Scaling mit Wiederholungen
- **Dateien:** `opt_compare_opt_off_raw.csv`, `opt_compare_opt_on_raw.csv`
- **Skript:** `job_scaling_rep.sh` (zweimal gelaufen), `-N 1`, REPS=3.
- **Methode:** Strong Scaling wie 3.1, aber **3 Wiederholungen** pro Rankzahl
  gegen Messrauschen. "opt off" = `halo_nonblocking=false`,
  `resort_every=0`; "opt on" = non-blocking + resort. Vorgesehene Auswertung
  über `aggregate_reps.py` (Median) + `analyze_opt_compare.py`.
- **Ranks:** 1, 2, 4, 8, 16, 32, 64, 128 (je 3 Wdh. = 24 Datenzeilen/Datei).
- **Spalten:** wie in Abschnitt 2.

### 3.6 MPI Opt-Vergleich – Konfigurationsmatrix @ 128 Ranks
- **Datei:** `opt_compare_n128.csv`
- **Skript:** manuell zusammengesetzte Datei (Header + 7 Blöcke), fest 128 Ranks.
- **Methode:** Systematische Variation von Kommunikationsmodus (blocking /
  non-blocking), Resort-Intervall (0 / 300) und PKA-Anzahl (1 / 8) bei
  festem N=128, um deren Einfluss zu isolieren.
- **Konfigurationen (Label):**
  - `blocking_no_resort_1pka`
  - `blocking_no_resort_8pka`
  - `blocking_resort300_1pka`
  - `blocking_resort300_8pka`
  - `non_blocking_no_resort_1pka`
  - `non_blocking_no_resort_8pka`
  - `non_blocking_resort300_1pka`
- **Spalten:** wie Abschnitt 2, je 1 Zeile (`ranks=128`) pro Konfiguration.

---

## 4. HPC-Performance-Messungen (CUDA / A100)

### 4.1 CUDA Haupt-Timing-Lauf (Headline-Punkt)
- **Datei:** `cuda_cascade_33377028.out` (Slurm-Job 33377028)
- **Skript:** `job_cuda.sh`, `-p gpu`, `--gres=gpu:a100:1`.
- **Methode:** Einzelner Timing-Lauf des Produktions-Binarys (`cascade_cuda`)
  auf dem Referenzpunkt. Der nachfolgende ncu-Teil des Skripts ist fehlge-
  schlagen (`%j` wurde nicht expandiert → `==ERROR== Unknown macro '%j'`),
  daher nur Timing hier; tiefes ncu-Profiling siehe 4.3.
- **Parameter:** N=1 999 396 (1414×1414), PKA_E=5000, dt=2e-4, 3000 Schritte,
  Zellen=2 405 340.
- **Erfasste Metriken:** Energie-Drift-Verlauf (alle 500 Schritte:
  0,0700 % → 0,2164 %), Zeitschleife gesamt (1,594 s), µs/Schritt (531,2),
  Durchsatz (3763,7 Matom/s), Force-Kernel isoliert (259,1 µs, 586,4 GB/s
  = 37,7 % Peak), H2D (184 MB, 14 GB/s) / D2H (16 MB, 13 GB/s), Speedup
  vs. 1 EPYC-Kern (723,4×) / 64-Kern-Knoten (14,13×).
- **Hinweis:** gleiche Konfiguration wie Block (1) des Sättigungs-Benchmarks
  (4.2), aber 1 Wiederholung statt 5; daher leicht abweichende Wandzeit
  (1,594 s vs. 1,5765 s Median). ERGEBNIS-BLOCK ist beschriftet mit
  "fuer Folie 9b".

### 4.2 CUDA Sättigungs- & Blockgrößen-Benchmark
- **Datei:** `cuda_bench_33375548.out` (Slurm-Job 33375548)
- **Skript:** `job_bench.sh`, `-p gpu`, `--gres=gpu:a100:1`.
- **Methode:** `cascade_cuda_bench` in 4 Blöcken (fp64- und fp32-Build via
  `-DUSE_FP32`), `nvcc -O3 -arch=sm_80`:
  1. **Headline fp64** – N=2 Mio (1414×1414), T=256, 3000 Schritte, 5 Wdh.
  2. **Headline fp32** – gleiche Konfiguration.
  3. **Sättigungs-Sweep fp64** – N ∈ {250k, 500k, 1 Mio, 2 Mio, 4 Mio,
     8 Mio, 16 Mio}, T=256, 400 Schritte, 3 Wdh.
  4. **Blockgrößen-Sweep fp64** – T ∈ {128, 256, 512, 1024}, N=2 Mio,
     400 Schritte, 3 Wdh.
- **Erfasste Metriken pro Lauf:** Reduktions-Drift, Median-Wandzeit (min/max),
  µs/Schritt, Durchsatz [Matom/s], Kernel-Breakdown (clear_cells, build_cells,
  force, energy, kick_drift, kick in µs/Aufruf), force-Anteil, H2D/D2H-Transfer
  (MB, GB/s), force-Bandbreite [% Peak], Speedup vs. 1 EPYC-Kern / 64-Kern-Knoten.
- **Maschinenlesbar:** `# CSV prec,N,T,steps,per_step_us,matom_s,force_us,
  step_us`-Zeilen, am Ende gesammelt.

### 4.3 CUDA ncu-Profiling (tief)
- **Dateien:** `cuda_ncu_33375549.out` (Job 33375549),
  `ncu_force_deep_summary.txt`, `ncu_report.csv`, `ncu_report2.csv`.
- **Skript:** `job_ncu_deep.sh`, `-p gpu`, `--gres=gpu:a100:1`.
- **Methode:** Nsight Compute (`ncu`) auf den `k_force`-Kernel des
  fp64-Benchmarks, N=2 Mio, T=256, 600 Schritte, 1 Wdh. Profilierung mit
  `--launch-skip 60 --launch-count 2` (2 repräsentative Launches nach Warmup),
  11 Passes.
- **Sections:** SpeedOfLight, ComputeWorkloadAnalysis, InstructionStats,
  Occupancy, MemoryWorkloadAnalysis.
- **`ncu_force_deep_summary.txt`:** Enthält den ERGEBNIS-BLOCK unter
  Profiling-Last (z. B. 244,5 Matom/s, force-Anteil 38,1 %, 581,4 GB/s
  = 37,4 % Peak) – bewusst langsamer als der unkonditionierte Benchmark
  (4.1 / 4.2) wegen Profiling-Overhead.
- **`ncu_report.csv` / `ncu_report2.csv`:** Vollständiger Metrik-Dump
  (transponiert: 2 Launches als Spalten 0/1, Metriken als Zeilen). Beide
  Dateien inhaltlich identisch (derselbe Lauf). Enthalten u. a. Estimated
  Speedup ~9 %, SM-Throughput ~45 %, DRAM-Throughput ~37 %, FP64-Pipe ~44 %,
  62 Register/Thread, Occupancy-Limit durch Register (4 Blöcke/SM),
  L1-/L2-Hitraten, Instruction-Mix, keine Register-/Shared-Memory-Spills.

---

## 5. Science-Schicht: Ensemble (Defektmuster / Power Law)

### 5.1 Ensemble – einzelne Kaskaden-Snapshots (Stichprobe)
- **Dateien:** `ensemble_run_0999_state_final.csv`,
  `ensemble_no_heal_run_0999_state_final.csv`.
- **Skript:** `job_ensemble.sh` (1000 unabhängige Kaskaden, `xargs -P 128`)
  + `scripts/python/run_one.py` + serieller Code `cascade_serial.c`.
- **Erzeugung:** `run_one.py` generiert pro Task-ID eine `.ini` mit
  deterministisch variiertem Einschlag (Zentrum 0,40–0,60), Winkel (0–360°)
  und Energie (log-uniform 200–2500), seed=1000+TASK, und ruft `cascade_serial`
  auf. Dies ist der Finalzustand von Task 999 (der 1000. Lauf).
- **Parameter:** NX=NY=250 (62 500 Atome, Dreiecksgitter), n_steps=12 000,
  dt=3e-4, K_SPRING=100, MAX_STRETCH=1.15, K_REP=400, REP_N=12, 1 PKA.
- **Spalten:** `x, y, disp, broken` (disp = Verschiebungsbetrag,
  broken = Anzahl gebrochener Bindungen an diesem Atom).
- **Variante:** `ensemble_run_0999_*` mit Healing, `ensemble_no_heal_run_0999_*`
  ohne Healing. Defekt-Atoms (broken>0): 1031 (heal) vs. 1073 (no_heal);
  Unterschiede im Verschiebungsfeld nur im Kaskadenkern (4. Nachkommastelle).
- **Hinweis:** Dies ist **kein** HPC-Performance-Datensatz, sondern eine
  Science-Ausgabe (Defektcluster-Snapshot). Die vollständigen Ensemble-Ordner
  (1000 Läufe + Power-Law-Auswertung via `analyze_ensemble.py` /
  `overlay_ccdf.py`) liegen nicht in `results/`; hier nur eine einzelne
  Stichprobe als Beispiel.

---

## 6. Bilder (`results/images/`)

5 PNG-Plots, erzeugt von den `scripts/python/analyze_*.py`-Skripten:

| Datei                                | gehört zu | Inhalt                                        |
|--------------------------------------|-----------|-----------------------------------------------|
| `compiler_opt.png`                   | 3.4       | Compiler-Flag-Sweep                           |
| `mpi_weak_scaling_2nodes_weak.png`   | 3.3       | Weak Scaling                                   |
| `opt_compare.png`                    | 3.5       | Opt on vs. off (Strong Scaling, 3 Wdh.)       |
| `opt_compare_n128.png`               | 3.6       | Konfigurationsmatrix @ 128 Ranks              |
| `strong_scaling_2nodes.png`          | 3.1 + 3.2 | Strong Scaling (1 und 2 Knoten kombiniert)    |

---

## 7. Querverweise & Hinweise

- **Problemgröße MPI (3.1–3.6):** durchgehend NX=NY=1414 (~2 Mio. Atome),
  Modell A, 1 PKA, E=5000, dt=2e-4, 3000 Schritte, seed=3 (Ausnahme: Weak
  Scaling 3.3 mit variierender NY).
- **CUDA-Referenzpunkt (4.x):** N=2 Mio (1414×1414), T=256 als Standard;
  fp64 und fp32 verglichen.
- **Wiederholungen:** Strong Scaling (3.1, 3.2) und Weak Scaling (3.3) je
  1 Lauf pro Punkt; Opt-Vergleich (3.5) und CUDA-Benchmark (4.1) mit 3–5 Wdh.
- **Leere/abgebrochene Job-Teile:** Der Timing-Teil von `job_cuda.sh`
  (4.1, Job 33377028) lief erfolgreich; der angehängte ncu-Aufruf desselben
  Jobs schlug fehl (`%j` nicht expandiert). Die tiefen ncu-Daten wurden
  stattdessen über `job_ncu_deep.sh` (4.3, Job 33375549) erhoben.
  Fruehere, leere `cuda_cascade_*.out`-Laeufe sind aus `results/` entfernt.
- **Trennung HPC / Science:** Abschnitt 3–4 = HPC-Kern (Performance);
  Abschnitt 5 = Science-Schicht (Defektmuster, Power Law). Vollständige
  Ensemble-Auswertungen und Power-Law-Fits liegen nicht in `results/`.
