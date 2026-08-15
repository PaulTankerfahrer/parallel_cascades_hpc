#!/usr/bin/env python3
"""preview_xyz.py -- schnelle Vorschau-Animation aus einer XYZ-Frame-Datei.
Faerbt nach Tempo (vmag) -> man sieht die Schockwelle laufen.
Aufruf: python3 preview_xyz.py vid.xyz"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

fn = sys.argv[1] if len(sys.argv) > 1 else "vid.xyz"

# --- XYZ-Frames parsen ---
frames = []
with open(fn) as f:
    lines = f.read().splitlines()
i = 0
while i < len(lines):
    n = int(lines[i]); i += 2                      # Atomzahl + Properties-Zeile
    block = lines[i:i+n]; i += n
    arr = np.array([ln.split()[1:] for ln in block], dtype=float)  # ohne Typ-Spalte
    frames.append(arr)   # Spalten: x y z disp vmag broken
print(f"{len(frames)} Frames, {frames[0].shape[0]} Atome")

# jeden 2. Frame nehmen (kleinere GIF)
frames = frames[::2]
x0, y0 = frames[0][:,0], frames[0][:,1]

fig, ax = plt.subplots(figsize=(5.5,5.5))
ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
sc = ax.scatter(x0, y0, c=frames[0][:,4], s=2, cmap="inferno", vmin=0, vmax=6)
ax.set_xlim(x0.min()-2, x0.max()+2); ax.set_ylim(y0.min()-2, y0.max()+2)
ttl = ax.set_title("t = 0.00")

def update(k):
    fr = frames[k]
    sc.set_offsets(np.c_[fr[:,0], fr[:,1]])
    sc.set_array(fr[:,4])           # vmag
    ttl.set_text(f"Frame {k*2}")
    return sc, ttl

anim = FuncAnimation(fig, update, frames=len(frames), interval=80, blit=False)
anim.save("vid_preview.gif", writer=PillowWriter(fps=12), dpi=90)
print("geschrieben: vid_preview.gif")

# 3 Standbilder (frueh / mitte / spaet) zur Kontrolle, gefaerbt nach disp
idx = [len(frames)//6, len(frames)//2, len(frames)-1]
fig2, axes = plt.subplots(1,3,figsize=(13,4.4))
for ax,k in zip(axes, idx):
    fr = frames[k]
    ax.scatter(fr[:,0], fr[:,1], c=fr[:,3], s=2, cmap="viridis", vmin=0, vmax=3)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"Frame {k*2} (Farbe = Verschiebung)")
fig2.tight_layout(); fig2.savefig("vid_snapshots.png", dpi=110)
print("geschrieben: vid_snapshots.png")
