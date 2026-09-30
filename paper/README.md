# Paper: Weiterentwicklung des Modells und echte Läufe

Aktiver Teil des Projekts. Hier wird das Minimalmodell aus dem Kursprojekt
physikalisch weiterentwickelt. Später laufen hier die großen Rechnungen auf dem
PC2, mit denen sich die Ergebnisse des Agentenlabors (Bachelorarbeit) prüfen
lassen. Hintergrund: [`../docs/UEBERBLICK.md`](../docs/UEBERBLICK.md).

**Stand (30.09.):** Der ZBL-Umbau ist gebaut und getestet. Die Kalibrierung der
Energieskala wartet auf eine Entscheidung zur Bindungsstärke, siehe
[`docs/haltepunkt_2026-09-24.md`](docs/haltepunkt_2026-09-24.md).

## Inhalt

| Pfad | Inhalt |
|---|---|
| `src/cascade_serial.c` | Serieller Simulator. `rep_model = r12` ist das alte Modell (bitgleich zum Kursprojekt), `rep_model = zbl` das neue: ZBL für alle Paare, adaptiver Zeitschritt, Lindhard-Bremsung, absorbierender Rand, Frühabbruch |
| `src/potential.h` | r⁻¹²- und ZBL-Potenzial mit glatter Abschaltung |
| `src/params.ini` | alle Schlüssel mit Erklärung (Werte: Benchmark aus dem Kursprojekt) |
| `scripts/ci/` | Tests: Energieerhaltung, Determinismus, Zweikörperstreuung gegen Theorie, ZBL im Gitter, Energiebilanz mit Bremsung, Valgrind |
| `scripts/python/zbl/` | Kalibrierung: P(E)-Scan der Verlagerungsschwelle, Wigner-Seitz-Defekte |
| `scripts/python/statistik/` | Ensemble-Treiber und Clauset-Fit, aus dem Kursprojekt übernommen. Wird für die Produktion erweitert |
| `scripts/python/diag_gebundene_nachbarn.py` | Diagnose des Durchflug-Fehlers im alten Modell |
| `results/` | Kalibrierung, Abbildungen D1, D2 |
| `docs/plan_zbl.md` | Plan und Fortschritt (Phasen 0–6) |
| `docs/befund_gebundene_nachbarn.md` | Der Modellfehler mit den durchlässigen Bindungen |
| `docs/haltepunkt_2026-09-24.md` | aktueller Stand, offene Entscheidung |

## Bauen und testen

Alle Pfade in Code-Kommentaren und Doku unter `paper/` gelten relativ zu
`paper/`. Befehle deshalb aus diesem Ordner starten:

```bash
cd paper
bash scripts/ci/build_and_test.sh
gcc -O3 -march=native -o cascade_serial src/cascade_serial.c -lm
```

Kalibrierung (Python mit numpy/scipy, das venv liegt im Wurzelverzeichnis):

```bash
cd paper
../.venv-powerlaw/bin/python scripts/python/zbl/kalibrierung_ed.py \
    --outdir results/zbl/kalibrierung --eps 0.3,0.4,0.5 --jobs 8
```
