#!/usr/bin/env python3
"""
cluster_from_zip.py -- liest die run_*_state_final.csv eines Ensembles direkt
aus einem ZIP-Archiv und erzeugt die gepoolte cluster_sizes.csv.

Identische Clusterdefinition wie scripts/python/analyze_ensemble.py:
  verlagert  := disp > 0.5 * L0
  Cluster    := Union-Find ueber Paare mit Abstand < 1.5 * L0

Aufruf:
  python3 cluster_from_zip.py ARCHIV.zip PRAEFIX_IM_ZIP AUSGABEVERZEICHNIS [--jobs 8]
Beispiel:
  python3 cluster_from_zip.py Messungen_USW.zip ensemble_no_heal/ results/statistik/ensemble_gepoolt/no_heal
"""
import argparse, io, os, sys, zipfile
from concurrent.futures import ProcessPoolExecutor

L0 = 1.0
DISP_THRESH = 0.5 * L0
LINK = 1.5 * L0

_Z = {}


def _init(zpath):
    _Z["z"] = zipfile.ZipFile(zpath)


def cluster_sizes(px, py, link=LINK):
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
    return list(roots.values())


def one_file(name):
    xs, ys = [], []
    with _Z["z"].open(name) as fh:
        txt = io.TextIOWrapper(fh, encoding="utf-8", errors="replace")
        txt.readline()                       # Header: x,y,disp,broken
        for line in txt:
            p = line.split(",")
            if len(p) < 3:
                continue
            if float(p[2]) > DISP_THRESH:
                xs.append(float(p[0]))
                ys.append(float(p[1]))
    return name, len(xs), cluster_sizes(xs, ys)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("zip")
    ap.add_argument("prefix")
    ap.add_argument("outdir")
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    a = ap.parse_args()

    with zipfile.ZipFile(a.zip) as z:
        names = sorted(n for n in z.namelist()
                       if n.startswith(a.prefix) and n.endswith("_state_final.csv"))
    if not names:
        sys.exit(f"keine *_state_final.csv unter '{a.prefix}' in {a.zip}")
    print(f"{len(names)} Laeufe unter '{a.prefix}', {a.jobs} Prozesse")
    os.makedirs(a.outdir, exist_ok=True)

    pooled, meta, per_run, done = [], [], [], 0
    with ProcessPoolExecutor(max_workers=a.jobs, initializer=_init,
                             initargs=(a.zip,)) as ex:
        for name, nmoved, sizes in ex.map(one_file, names, chunksize=4):
            pooled.extend(sizes)
            per_run.append((name, sizes))
            meta.append((os.path.basename(name), nmoved, len(sizes),
                         max(sizes) if sizes else 0))
            done += 1
            if done % 100 == 0 or done == len(names):
                print(f"  {done}/{len(names)}  ->  {len(pooled)} Cluster", flush=True)

    pooled.sort(reverse=True)
    out = os.path.join(a.outdir, "cluster_sizes.csv")
    with open(out, "w") as f:
        f.write("size\n")
        for s in pooled:
            f.write(f"{s}\n")
    # Cluster einzeln mit Lauf-Zuordnung sichern (fuer energieaufgeloeste Analysen)
    with open(os.path.join(a.outdir, "cluster_by_run.csv"), "w") as f:
        f.write("task,size\n")
        for name, sizes in per_run:
            task = int(os.path.basename(name).split("_")[1])
            for sz in sizes:
                f.write(f"{task},{sz}\n")

    meta.sort()
    with open(os.path.join(a.outdir, "run_meta.csv"), "w") as f:
        f.write("run,n_moved,n_clusters,max_cluster\n")
        for r in meta:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]}\n")
    print(f"\n{len(pooled)} Cluster (min={min(pooled)}, max={max(pooled)})")
    print("geschrieben:", out)


if __name__ == "__main__":
    main()
