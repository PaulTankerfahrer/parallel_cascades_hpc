#!/usr/bin/env python3
"""
run_fixed_energy.py -- Ensemble bei FESTER PKA-Energie.

Gegenprobe zum urspruenglichen Ensemble: dort wurde die Energie log-uniform
gewuerfelt, wodurch die Stichprobenwahl selbst ein Potenzgesetz erzeugt.
Hier ist die Energie konstant; variiert werden nur Einschlagort, Winkel und
Seed. Damit misst man die Fragmentierungsstatistik der Kaskade und nicht die
Energieverteilung.

Die Einschlagparameter haengen nur von der Task-ID ab, nicht von der Energie:
Lauf 42 hat bei allen Energien dieselbe Geometrie (gepaarter Vergleich).

Pro Lauf wird ausgewertet und danach die grosse Rohdatei geloescht; behalten
wird eine getrimmte Fassung mit allen Atomen ueber disp > 0.4 (Faktor ~30
kleiner, verlustfrei fuer jede Schwelle >= 0.4).

Aufruf:
  python3 run_fixed_energy.py --outdir results/statistik/ensemble_fest \
      --energies 1300,1800,2400 --runs 300 --jobs 8
"""
import argparse, math, os, random, shutil, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(REPO, "src", "cascade_serial.c")

# (Schwelle auf disp, Linkradius) -- identisch benannt wie in analyze_ensemble.py
COMBOS = [(0.5, 1.5), (0.5, 1.0), (1.0, 1.5), (1.0, 1.0), (1.5, 1.0), (2.0, 1.0)]
TRIM = 0.4          # Atome darunter werden nicht aufbewahrt
L0 = 1.0

_CFG = {}


def _init(cfg):
    _CFG.update(cfg)


def build(binpath):
    os.makedirs(os.path.dirname(binpath), exist_ok=True)
    cc = shutil.which("gcc") or shutil.which("cc")
    subprocess.run([cc, "-O3", "-march=native", "-o", binpath, SRC, "-lm"], check=True)
    return binpath


def cluster_sizes(px, py, link):
    n = len(px)
    if n == 0:
        return []
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    buckets = {}
    for i in range(n):
        buckets.setdefault((int(px[i] // link), int(py[i] // link)), []).append(i)
    link2 = link * link
    for i in range(n):
        cx, cy = int(px[i] // link), int(py[i] // link)
        for ax in (cx - 1, cx, cx + 1):
            for ay in (cy - 1, cy, cy + 1):
                for j in buckets.get((ax, ay), ()):
                    if j <= i:
                        continue
                    if (px[i] - px[j]) ** 2 + (py[i] - py[j]) ** 2 < link2:
                        ra, rb = find(i), find(j)
                        if ra != rb:
                            parent[ra] = rb
    roots = {}
    for i in range(n):
        r = find(i)
        roots[r] = roots.get(r, 0) + 1
    return sorted(roots.values(), reverse=True)


def impact(task):
    """Einschlagparameter -- bewusst unabhaengig von der Energie."""
    rng = random.Random(50000 + int(task))
    return (rng.uniform(0.40, 0.60), rng.uniform(0.40, 0.60), rng.uniform(0.0, 360.0))


def one_run(arg):
    energy, task = arg
    c = _CFG
    edir = os.path.join(c["outdir"], f"E{int(round(energy)):04d}")
    dispdir = os.path.join(edir, "disp")
    px, py, ang = impact(task)
    scratch = tempfile.mkdtemp(prefix="fixE_")
    prefix = os.path.join(scratch, f"run_{task:04d}")
    cfg = prefix + ".ini"
    with open(cfg, "w") as f:
        f.write(f"""[grid]
NX = {c['nx']}
NY = {c['nx']}
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
dt={c['dt']}
n_steps={c['steps']}
[output]
log_every={c['steps']}
dump_every=0
out_prefix={prefix}
[ensemble]
seed={task}
healing=true
""")
    try:
        subprocess.run([c["bin"], cfg], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=True)
        xs, ys, ds, bs = [], [], [], []
        with open(prefix + "_state_final.csv") as fh:
            fh.readline()
            for line in fh:
                p = line.split(",")
                d = float(p[2])
                if d > TRIM:
                    xs.append(float(p[0])); ys.append(float(p[1]))
                    ds.append(d); bs.append(int(p[3]))
        with open(os.path.join(dispdir, f"run_{task:04d}.csv"), "w") as fh:
            fh.write("x,y,disp,broken\n")
            for i in range(len(xs)):
                fh.write(f"{xs[i]:.4f},{ys[i]:.4f},{ds[i]:.4f},{bs[i]}\n")
        out = {}
        for thr, link in COMBOS:
            sel = [i for i in range(len(xs)) if ds[i] > thr]
            out[(thr, link)] = cluster_sizes([xs[i] for i in sel],
                                             [ys[i] for i in sel], link)
        return task, energy, px, py, ang, out
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def tag(thr, link):
    return f"thr{int(thr*100):03d}_link{int(link*100):03d}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--energies", default="1300,1800,2400")
    ap.add_argument("--runs", type=int, default=300)
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--nx", type=int, default=250)
    ap.add_argument("--steps", type=int, default=12000)
    ap.add_argument("--dt", default="3e-4")
    ap.add_argument("--bin", default=None)
    a = ap.parse_args()

    binpath = a.bin or build(os.path.join(tempfile.gettempdir(),
                                          "cascade_fixed_bin", "cascade_serial"))
    energies = [float(e) for e in a.energies.split(",")]
    os.makedirs(a.outdir, exist_ok=True)
    for e in energies:
        d = os.path.join(a.outdir, f"E{int(round(e)):04d}")
        os.makedirs(os.path.join(d, "disp"), exist_ok=True)

    cfg = dict(outdir=a.outdir, nx=a.nx, steps=a.steps, dt=a.dt, bin=binpath)
    jobs = [(e, t) for e in energies for t in range(a.runs)]
    print(f"{len(jobs)} Laeufe ({len(energies)} Energien x {a.runs}), "
          f"{a.jobs} parallel, healing=true, NX={a.nx}", flush=True)

    files, meta = {}, {}
    for e in energies:
        d = os.path.join(a.outdir, f"E{int(round(e)):04d}")
        meta[e] = open(os.path.join(d, "run_meta.csv"), "w")
        meta[e].write("run,energy,px,py,angle,n_moved,n_clusters,max_cluster\n")
        for thr, link in COMBOS:
            fh = open(os.path.join(d, f"clusters_{tag(thr, link)}.csv"), "w")
            fh.write("task,size\n")
            files[(e, thr, link)] = fh

    done = 0
    with ProcessPoolExecutor(a.jobs, initializer=_init, initargs=(cfg,)) as ex:
        for task, e, px, py, ang, out in ex.map(one_run, jobs, chunksize=1):
            base = out[(0.5, 1.5)]
            meta[e].write(f"{task},{e:.1f},{px:.4f},{py:.4f},{ang:.2f},"
                          f"{sum(base)},{len(base)},{base[0] if base else 0}\n")
            for (thr, link), sizes in out.items():
                fh = files[(e, thr, link)]
                for s in sizes:
                    fh.write(f"{task},{s}\n")
            done += 1
            if done % 25 == 0:
                for fh in list(files.values()) + list(meta.values()):
                    fh.flush()
                print(f"  {done}/{len(jobs)}", flush=True)
    for fh in list(files.values()) + list(meta.values()):
        fh.close()
    print("fertig ->", a.outdir, flush=True)


if __name__ == "__main__":
    main()
