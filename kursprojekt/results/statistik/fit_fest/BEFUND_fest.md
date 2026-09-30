# Befund: Ensemble bei fester PKA-Energie

**Datum:** 2026-09-01 · **Daten:** `results/statistik/ensemble_fest` · 900 Läufe
(3 × 300 bei E = 1300 / 1800 / 2400), `NX = 250`, `healing = true`,
`dt = 3e-4`, 12 000 Schritte.

Einschlagort, Winkel und Seed hängen **nur von der Lauf-Nummer** ab, nicht von
der Energie. Lauf 42 hat bei allen drei Energien exakt dieselbe Geometrie — der
Vergleich über die Energien ist damit gepaart und der Energieeffekt isoliert.

---

## 1. Warum dieser Lauf nötig war

Im ursprünglichen Ensemble wurde die PKA-Energie log-uniform aus [200, 2500]
gezogen. Zusammen mit `n_max ~ E^2,9` (r = 0,97) erzeugt allein diese Wahl der
Stichprobe ein `p(n) ~ n^-1`, unabhängig von der Physik. Der Exponent S = 1,41
maß daher die Energieverteilung, nicht die Kaskade. Details:
`results/statistik/fit_gepoolt/BEFUND_gepoolt.md`.

Bei fester Energie fällt dieser Mechanismus weg. Was übrig bleibt, ist die
Fragmentierungsstatistik der Kaskade selbst.

---

## 2. Die bisherige Clusterdefinition trägt bei fester Energie nicht

`disp > 0,5 · L0`, Linkradius `1,5 · L0`:

| E | Cluster | pro Lauf | S (95 %-Bootstrap-CI) | KS-D | GoF p |
|---|--:|--:|--:|--:|--:|
| 1300 | 1131 | 3,8 | 1,400 [1,379, 1,426] | 0,173 | 0,000 |
| 1800 | 1100 | 3,7 | 1,354 [1,333, 1,378] | 0,201 | 0,000 |
| 2400 | 675 | 2,2 | 1,234 [1,220, 1,250] | 0,300 | 0,000 |

Die 1,4 taucht wieder auf, aber der KS-Abstand ist 4- bis 7-mal so groß wie im
gepoolten Ensemble (dort D = 0,041). Das ist kein Fit.

**Der Grund ist geometrisch.** Mit diesen Schwellen verschmilzt die gesamte
Schadenszone zu einem einzigen zusammenhängenden Klumpen: über die 1000 Läufe
des alten Ensembles liegen **92 % aller verlagerten Atome im größten Cluster**,
und 365 Läufe haben überhaupt nur einen. Die „Verteilung" besteht aus einem
Haufen Einer und Zweier plus einem Riesenblob — zweigipflig, und ein
Potenzgesetz-Schätzer darauf liefert eine Zahl ohne Bedeutung.

---

## 3. Mit einer feineren Definition entsteht eine echte Verteilung

`disp > 1,0 · L0`, Linkradius `1,0 · L0` — dieselben Läufe, dieselben Atome:

| E | Cluster | pro Lauf | größter | x_min | S (95 %-CI) | KS-D | GoF p |
|---|--:|--:|--:|--:|--:|--:|--:|
| 1300 | 9 281 | 30,9 | 43 | 2 | 2,498 [2,230, 2,950] | 0,032 | 0,000 |
| 1800 | 14 588 | 48,6 | 113 | 2 | 2,191 [2,076, 2,279] | 0,017 | 0,000 |
| 2400 | 21 580 | 71,9 | 330 | 1 | 1,996 [1,982, 2,010] | 0,013 | 0,000 |

30–70 getrennte Cluster pro Lauf, KS-Abstand um Faktor 10–20 besser. **Das ist
die Größe, die Sand et al. zählen** — getrennte Defektcluster, nicht die
Ausdehnung der Schadenszone.

Likelihood-Ratio-Tests bei E = 2400 (oberhalb x_min = 1, also auf allen Daten):

| Alternative | R | p | wer gewinnt |
|---|--:|--:|---|
| Exponential | +41,63 | < 1e-4 | Potenzgesetz |
| Stretched Exponential | +12,41 | < 1e-4 | Potenzgesetz |
| Lognormal | +1,80 | 0,072 | Potenzgesetz, aber nicht signifikant |
| Potenzgesetz + Cutoff | −4,30 | 0,003 | Cutoff — **aber bei n ≈ 2946** |

Der Cutoff liegt damit fast eine Größenordnung jenseits des größten je
beobachteten Clusters (330). Im Datenbereich ist er wirkungslos. Bei den
niedrigeren Energien ist das anders: dort liegt der Cutoff bei n ≈ 13 (E = 1300)
bzw. n ≈ 71 (E = 1800), also **mitten in den Daten**, und Lognormal und
Stretched Exponential schlagen das Potenzgesetz.

---

## 4. Der Exponent ist keine Konstante des Modells

Dieselben 900 Läufe, dieselben Atome, sechs Clusterdefinitionen:

| Definition | E=1300 | E=1800 | E=2400 |
|---|--:|--:|--:|
| `disp>0,5` / link 1,5 | 1,40 | 1,35 | 1,23 |
| `disp>1,0` / link 1,5 | 1,57 | 1,50 | 1,45 |
| `disp>0,5` / link 1,0 | 1,98 | 2,10 | 2,25 |
| `disp>1,0` / link 1,0 | 2,50 | 2,19 | 2,00 |
| `disp>1,5` / link 1,0 | 2,78 | 2,51 | 2,28 |
| `disp>2,0` / link 1,0 | 3,30 | 3,02 | 2,77 |

S läuft von **1,23 bis 3,30**. Die Konvention entscheidet stärker als die
Physik. Eine Zahl „S = 1,4" ohne Angabe der Clusterdefinition ist bedeutungslos —
und genau so steht sie derzeit auf der Folie.

---

## 5. Was ehrlich gesagt werden muss: auch die feinere Definition ist kein Potenzgesetz

Der GoF-Test verwirft alle sechs Datensätze (p = 0,000). Zwei Einordnungen dazu:

- **Entlastend:** bei n = 21 580 verwirft der Test schon winzige Abweichungen
  (D_beobachtet = 0,0127 gegen D_synthetisch = 0,0024). Clauset warnt selbst
  davor, p bei sehr großem n als Urteil zu lesen.
- **Belastend:** die CCDF in `fixed_energy_ccdf.png` ist sichtbar gekrümmt, nicht
  gerade. Das ist kein Rauschen, das ist Form.

Der brauchbare Befund steckt im **Trend**: mit wachsender Kaskade fallen S und D
gemeinsam und monoton.

| E | verlagerte Atome / Lauf | S | KS-D |
|---|--:|--:|--:|
| 1300 | 474 | 2,50 | 0,032 |
| 1800 | 1399 | 2,19 | 0,017 |
| 2400 | 2646 | 2,00 | 0,013 |

Die Verteilung **nähert sich** einem Potenzgesetz, je größer die Kaskade wird,
hat es bei den erreichbaren Kaskadengrößen aber nicht erreicht. Konvergiert ist
S ebenfalls nicht — 2,50 → 2,19 → 2,00 fällt noch.

---

## 6. Konsequenz für Präsentation und Arbeit

**Der Vergleich mit Sand et al. wird wieder zulässig** — feste Energie und eine
Clusterdefinition, die getrennte Defektcluster zählt. Er fällt aber anders aus
als behauptet: der 2D-Wert liegt bei S ≈ 2,0, also **über** dem 3D-Wert 1,63.
Die Folienaussage „der 2D-Wert liegt konsistent darunter" ist damit in zwei
Punkten falsch — Zahl und Richtung.

Formulierung, die trägt:

> Bei fester PKA-Energie und einer Clusterdefinition, die getrennte
> Defektcluster auflöst (`disp > 1,0 L0`, Linkradius `1,0 L0`), folgt die
> Clustergrößenverteilung näherungsweise einem Potenzgesetz mit S = 2,00 ± 0,01
> (E = 2400, 21 580 Cluster aus 300 Kaskaden). Der Exponent fällt systematisch
> mit der Kaskadengröße (2,50 → 2,19 → 2,00) und ist nicht konvergiert; der
> strenge Goodness-of-Fit-Test verwirft das reine Potenzgesetz. Der 3D-Wert von
> Sand et al. (1,63) liegt unterhalb — konsistent mit der geringeren
> Konnektivität in zwei Dimensionen.

**Nächster Schritt** für die Arbeit: größere Box und höhere Energie, um zu
prüfen, ob S konvergiert. Der Trend legt nahe, dass die 2D-Kaskade bei diesen
Größen noch nicht im asymptotischen Regime ist.

---

## 7. Zwei Werkzeugfallen, die unterwegs auffielen

- `powerlaw` sucht `alpha` nur in **[0, 3]** und meldet bei steileren
  Verteilungen kommentarlos exakt `3.000`. Bei `disp > 2,0` war die freie MLE
  tatsächlich 3,30 / 3,02 / 2,77. `clauset_refit.py` und
  `analyze_fixed_energy.py` rechnen jetzt jeden Wert nahe 3,000 mit
  `mle_unbounded()` frei nach.
- Bei zwei sehr großen Datensätzen konvergierte der `truncated_power_law`-Fit
  nicht und meldete R = +757 gegen das Potenzgesetz. Das ist unmöglich — das
  genestete Modell kann nie schlechter sein. Solche Zellen stehen jetzt auf
  „n. k.".

---

## 8. Reproduktion

```bash
# Ensemble (900 Laeufe, ~35 min auf 8 Kernen)
python3 scripts/python/statistik/run_fixed_energy.py \
    --outdir results/statistik/ensemble_fest --energies 1300,1800,2400 \
    --runs 300 --jobs 8

# Uebersicht ueber alle Energien x Clusterdefinitionen (Sekunden)
.venv-powerlaw/bin/python scripts/python/statistik/analyze_fixed_energy.py \
    results/statistik/ensemble_fest --outdir results/statistik/fit_fest

# Voller Clauset-Fit mit Bootstrap und GoF (~7 min)
.venv-powerlaw/bin/python scripts/python/statistik/clauset_refit.py \
    results/statistik/ensemble_fest/E1300/clusters_thr050_link150.csv "E=1300, alte Definition" \
    ... \
    --outdir results/statistik/fit_fest --bootstrap 300 --gof 200 --jobs 8

# CCDF-Vergleichsabbildung
.venv-powerlaw/bin/python scripts/python/statistik/fig_fixed_energy_ccdf.py \
    results/statistik/ensemble_fest results/statistik/fit_fest/fixed_energy_ccdf.png
```

Dateien: `fixed_energy_uebersicht.md` (alle Kombinationen),
`clauset_refit.md` / `.json` (voller Fit), `refit.log` (Rohausgabe),
`fixed_energy_ccdf.png`, `clauset_ccdf_alternatives.png`,
`clauset_xmin_scan.png`, `clauset_bootstrap_alpha.png`.
