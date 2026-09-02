# Physik-Präsentation & Kaskaden-Videos

Kurze, physikorientierte Vorstellung der Kollisionskaskaden-Simulation
(ohne den HPC-/Programmierteil) samt vier OVITO-Videos.

## Inhalt

| Datei | Zweck |
|---|---|
| `physik_praesentation.tex` / `.pdf` | Beamer-Folien: 8 Folien + 5 Reservefolien im Anhang |
| `videos/cascade_1_geschwindigkeit.mp4` | Geschwindigkeitsfeld — thermischer Spike |
| `videos/cascade_2_schaden.mp4` | Schadensbildung — gerissene Bindungen |
| `videos/cascade_3_cluster.mp4` | Defektcluster nach der Definition der Auswertung (`disp > 1,0`, Linkradius `1,0`) |
| `videos/cascade_4_mehrere_pka.mp4` | Fünf gleichzeitige Einschläge |
| `stills/` | Standbilder für die Folien (aus den Videos extrahiert) |
| `single.ini`, `multi.ini` | Simulations-Konfigurationen der beiden Läufe. `single.ini` entspricht dem Fixed-Energy-Ensemble (NX=250, `healing=true`, E=2400, Einschlagparameter von Lauf 276). |
| `build_videos.sh` | Erzeugt alles neu (Sim → Render → Standbilder) |

## Videos abspielen

Die Videos sind **eigenständige MP4-Dateien** in `videos/`. Die PDF bindet sie
auf zwei Wegen ein: als Movie-Annotation über dem Standbild (Inline-Wiedergabe)
und als Launch-Action in der Zeile darunter (öffnet den Systemplayer).

**Nicht jeder Betrachter kann das.** Firefox bzw. pdf.js ignoriert beides
stillschweigend — der Link ist sichtbar, der Klick bewirkt nichts.

| Weg | Wiedergabe |
|---|---|
| `pdfpc physik_praesentation.pdf` | **Empfohlen zum Vortragen.** Spielt die Videos inline ab, dazu Referentenansicht mit Timer. |
| Okular, Präsentationsmodus | Spielt die Videos inline ab. |
| Acrobat Reader | „extern öffnen`` startet den Systemplayer. |
| Firefox / pdf.js | Keine Wiedergabe. Stattdessen `videos/index.html` öffnen. |
| beliebig | MP4s direkt aus `videos/` mit mpv oder VLC. |

Ohne Installation: **`videos/index.html`** im Browser öffnen — eine Seite mit
allen vier Videos und ihren Bildunterschriften.

```bash
firefox presentation/videos/index.html
```

## Neu erzeugen

```bash
# OVITO-Python bereitstellen (einmalig):
python3 -m venv ovito-venv && ovito-venv/bin/pip install ovito

# Videos + Standbilder erzeugen:
OVITO_PY=ovito-venv/bin/python ./build_videos.sh

# Folien bauen:
pdflatex physik_praesentation.tex   # 2x für Referenzen
```

> Hinweis: Es genügt das PyPI-Paket `ovito` (headless-Rendering).
> Die OVITO-Desktop-GUI allein bringt **kein** skriptbares Python mit.

## Aufbau der Folien

1. Warum ein Minimalmodell? — Wolfram, Neutronen, Kaskaden; Verteilungen
   brauchen Ensembles, volle MD ist dafür zu teuer.
2. Das Modell — Dreiecksgitter, brechbare Feder + steile Abstoßung.
3. Video: Der thermische Spike.
4. Video: Die Defektcluster entstehen — Atome mit `disp > 1,0`, nach Cluster eingefärbt.
5. Ein Ergebnis, das zu schön war — der alte Exponent und seine Prüfung.
6. Das saubere Ensemble — feste Energie, 3 × 300 Läufe.
7. Was die Zahlen sagen — und was nicht.
8. Von der Korrektur zur Bachelorarbeit.

Danach folgen nach `\appendix` fünf Reservefolien für Rückfragen (Weg von der
CSV zum Cluster, Inhalt der Endzustandstabelle, beide Clusterdefinitionen im
Vergleich, zwei weitere Videos). Sie werden nicht mitgezählt.

## Physikalische Kernaussagen

- Ein PKA (Primäres Rückstoßatom) startet eine verzweigende Kollisionskaskade.
- Sichtbar werden: heißer Spike (Geschwindigkeit), Defektbildung (gerissene
  Bindungen) und bleibende Verschiebung.
- Der frühere Befund „Clustergrößen folgen einem Potenzgesetz mit $S\approx1{,}4$"
  ist **nicht haltbar**: Der Goodness-of-Fit-Test verwirft ihn, und der Exponent
  entsteht aus der log-uniform gezogenen PKA-Energie zusammen mit $n\sim E^{2,9}$,
  nicht aus der Kaskadenphysik.
- Bei **fester** PKA-Energie und einer Clusterdefinition, die getrennte
  Defektcluster auflöst (`disp > 1,0`, Linkradius `1,0`), ergibt sich
  $S = 2{,}00 \pm 0{,}01$ bei $E=2400$ — mit fallendem Trend über die
  Kaskadengröße (2,50 → 2,19 → 2,00) und daher nicht konvergiert.
- Vergleichswert 3D-Wolfram (Sand et al. 2013): $S=1{,}63$ — der 2D-Wert liegt
  **darüber**, konsistent mit der geringeren Konnektivität in zwei Dimensionen.

Ausführlich: `../results/statistik/fit_gepoolt/BEFUND_gepoolt.md` und
`../results/statistik/fit_fest/BEFUND_fest.md`.
