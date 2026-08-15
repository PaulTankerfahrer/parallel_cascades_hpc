#!/usr/bin/env python3
"""
run_one.py -- EINE Kaskade fuer eine Slurm-Job-Array-Task.
Variiert deterministisch (aus der Task-ID) Einschlagort, Winkel und Energie
(log-uniform = Rueckstoss-Spektrum-Proxy) und ruft die serielle Version auf.
Ausgabe: <outdir>/run_<id>_state_final.csv  (+ _energy, _broken_bonds)

Aufruf:  python3 run_one.py <task_id> [outdir]
Groesse/Schritte via Umgebungsvariablen ENS_NX / ENS_STEPS / ENS_DT
(Default = Produktionswerte; fuer schnelle Tests kleiner setzen).
"""
import sys, os, subprocess, random, math

TASK   = int(sys.argv[1])
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else "ensemble"
os.makedirs(OUTDIR, exist_ok=True)

NX     = int(os.environ.get("ENS_NX", 250))     # gross genug -> Kaskade bleibt lokal
NSTEPS = int(os.environ.get("ENS_STEPS", 12000))
DT     = os.environ.get("ENS_DT", "3e-4")       # klein genug gegen Blow-up
E_MIN, E_MAX = 200.0, 2500.0                     # Energiebereich (log-uniform)

rng = random.Random(1000 + TASK)                 # reproduzierbar pro Task
energy = math.exp(rng.uniform(math.log(E_MIN), math.log(E_MAX)))
px  = rng.uniform(0.40, 0.60)                    # zentrales Fenster -> Rand frei
py  = rng.uniform(0.40, 0.60)
ang = rng.uniform(0.0, 360.0)

prefix = os.path.join(OUTDIR, f"run_{TASK:04d}")
cfg    = prefix + ".ini"
with open(cfg, "w") as f:
    f.write(f"""[grid]
NX = {NX}
NY = {NX}
L0 = 1.0
[potential]
K_SPRING=100.0
MAX_STRETCH=1.15
RCUT=0.9
K_REP=400.0
REP_N=12.0
[pka]
pka_x={px:.4f}
pka_y={py:.4f}
pka_energy={energy:.2f}
pka_angle={ang:.2f}
pka_mass=1.0
[dynamics]
dt={DT}
n_steps={NSTEPS}
[output]
log_every={NSTEPS}
dump_every=0
out_prefix={prefix}
[ensemble]
seed={TASK}
""")
subprocess.run(["./cascade_serial", cfg],
               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
os.remove(cfg)
