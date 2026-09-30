# Clusterstatistik: Potenzgesetz-Analyse

Skripte für die Auswertung der Defekt-Clustergrößen nach Clauset, Shalizi &
Newman (2009). Sie brauchen das Paket `powerlaw`; im Repo liegt dafür ein
eigenes venv:

```bash
python3 -m venv .venv-powerlaw && .venv-powerlaw/bin/pip install powerlaw
```

Aufgerufen wird immer aus dem Projektwurzelverzeichnis.

## Zweig 1 — gepooltes Ensemble (Originaldaten des Berichts)

| Skript | Zweck |
|---|---|
| `cluster_from_zip.py` | Liest die 1000 Kaskaden direkt aus dem Nextcloud-Zip und clustert sie, ohne es zu entpacken. → `results/statistik/ensemble_gepoolt/` |
| `clauset_refit.py` | Der vollständige Fit: diskrete MLE, x_min-Scan, Bootstrap-CI, Goodness-of-Fit, Likelihood-Ratio-Tests gegen vier Alternativen. |
| `energy_control.py` | Kontrolltest: derselbe Fit gepoolt gegen einzelne Energiebänder. Zeigt, dass der Exponent aus der Energieziehung stammt. |

## Zweig 2 — Ensemble bei fester Energie

| Skript | Zweck |
|---|---|
| `run_fixed_energy.py` | Treiber: baut das serielle Binary und rechnet N Läufe je Energie, variiert nur Einschlagort, Winkel und Seed. Clustert jeden Lauf gleich in sechs Definitionen. |
| `analyze_fixed_energy.py` | Schnellübersicht über alle Energien × Clusterdefinitionen (Sekunden statt Minuten, ohne Bootstrap). |
| `fig_fixed_energy_ccdf.py` | Die CCDF-Vergleichsabbildung. `--panel=alt\|fein` für ein einzelnes Panel, `--folie` für die Beamer-Variante mit größerer Schrift. |

## Abbildungen für die Präsentation

| Skript | Zweck |
|---|---|
| `fig_csv_to_cluster.py` | Zweiteilige Abbildung: rohes `disp`-Feld eines Laufs neben den daraus gebildeten Clustern. |

## Zwei Fallstricke des `powerlaw`-Pakets

- Es sucht `alpha` nur in **[0, 3]** und meldet bei steileren Verteilungen
  kommentarlos exakt `3.000`. `clauset_refit.py` und `analyze_fixed_energy.py`
  rechnen jeden Wert nahe 3,000 mit `mle_unbounded()` frei nach und markieren
  ihn mit ⚠.
- Der `truncated_power_law`-Fit konvergiert nicht immer und meldet dann ein
  positives R gegen das Potenzgesetz — für ein genestetes Modell unmöglich.
  Solche Zellen erscheinen als „n. k.".

Erklärung aller Fachbegriffe: [`../../../docs/clauset_glossar.md`](../../../../docs/clauset_glossar.md).
