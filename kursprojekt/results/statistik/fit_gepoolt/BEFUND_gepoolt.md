# Schritt 0 -- Clauset-Refit: Ergebnis und Konsequenz für die Präsentation

Datenbasis: die **Originaldaten** des Projektberichts, je 1000 Kaskaden mit und
ohne Healing, aus `Messungen_USW.zip` (Nextcloud). Die Clusterzahlen 2391 (heal)
und 2335 (no_heal) wurden aus den Rohläufen exakt reproduziert.
Werkzeug: Python-Paket `powerlaw` 2.0.0 (Alstott et al. 2014), Methodik nach
Clauset, Shalizi & Newman (2009).

---

## Kurzfassung

Der Exponent des Berichts ist **reproduzierbar und präzise**. Aber der strenge
Test zeigt: Er ist **kein Ergebnis der Kaskadenphysik**, sondern überwiegend ein
Abbild der log-uniformen PKA-Energieverteilung, aus der das Ensemble gezogen
wurde. Die Folie in ihrer bisherigen Form ist nicht haltbar.

---

## 1. Der Exponent ist reproduzierbar

| | Bericht | Refit (diskretes MLE) |
|---|---|---|
| Mit Healing, x_min = 4 | S = 1,41 | S = 1,391 ± 0,012 |
| Ohne Healing, x_min = 4 | S = 1,38 | S = 1,361 ± 0,011 |

Die Berichtswerte 1,41 / 1,38 lassen sich **exakt** rekonstruieren, wenn man den
**kontinuierlichen** MLE-Schätzer `S = 1 + n / sum(ln(x/x_min))` verwendet
(nachgerechnet: 1,4114 und 1,3790). Clustergrößen sind aber ganze Zahlen; der
korrekte **diskrete** Schätzer liegt rund 0,02 tiefer. Das ist kein Fehler des
Berichts, sondern ein bekannter kleiner Bias der kontinuierlichen Näherung.

**Ohne** die Beschränkung auf x_min = 4 wählt die KS-Minimierung x_min = 1 und
liefert S = 1,409 (heal) / 1,390 (no_heal), Bootstrap-95-%-CI [1,393; 1,426]
bzw. [1,374; 1,406]. Der Fehler ±0,03 des Berichts war eher zu konservativ --
die statistische Präzision ist besser als angegeben.

---

## 2. Aber: Das reine Potenzgesetz wird verworfen

**Goodness-of-Fit-Test** (1000 semiparametrische Bootstrap-Datensätze):
**p = 0,000** für beide Varianten. Clauset-Kriterium ist p > 0,1. Das reine
Potenzgesetz wird als Erzeugungsmodell abgelehnt.

**Likelihood-Ratio-Tests** (alle oberhalb desselben x_min):

| Alternative | R (heal) | p | Ergebnis |
|---|---:|---:|---|
| Exponential | +50,2 | < 10⁻⁴ | Potenzgesetz klar besser |
| Lognormal | −9,5 | < 10⁻⁴ | **Lognormal besser** |
| Stretched Exponential | −2,4 | 0,015 | Alternative besser |
| Potenzgesetz + Cutoff | −45,0 | < 10⁻⁴ | **Cutoff-Variante besser** |

Die einzige Alternative, die das Potenzgesetz schlägt *nicht*, ist die
Exponentialverteilung -- also die Mindesthürde. Gegen die ernsthaften
schwerschwänzigen Konkurrenten verliert es.

Die `x_min`-Diagnose (`clauset_xmin_scan.png`) zeigt warum: Der KS-Abstand
**steigt monoton** mit x_min. Es gibt keinen Bereich, in dem das Potenzgesetz
besonders gut passt -- der Fit wird nach oben hin nur schlechter.

---

## 3. Der eigentliche Befund: Das Potenzgesetz kommt aus dem Sampling

`run_one.py` zieht die PKA-Energie **log-uniform** aus [200, 2500]. Und:

- Der größte Cluster eines Laufs skaliert wie **n_max ~ E^2,9** mit **r = 0,97**.
  Die Clustergröße ist damit fast vollständig durch die Anregungsenergie bestimmt.
- Ist E log-uniform verteilt und gilt n ~ E^k, dann ist ln(n) gleichverteilt.
  Eine gleichverteilte Log-Größe bedeutet aber **p(n) ~ n⁻¹** -- ein Potenzgesetz,
  das rein aus der Wahl der Energieverteilung folgt, unabhängig von k und von
  jeder Kaskadenphysik.

**Kontrolltest** (`energy_control.py`): derselbe Fit in drei Energiebändern statt
gepoolt.

| Datensatz | gepoolt | Band 1 | Band 2 | Band 3 | KS D gepoolt | KS D Bänder |
|---|---:|---:|---:|---:|---:|---|
| Mit Healing | **1,41** | 2,96 | 1,80 | 2,45 | 0,041 | 0,093 -- 0,186 |
| Ohne Healing | **1,39** | 2,92 | 1,76 | 2,16 | 0,051 | 0,107 -- 0,176 |

In **keinem** einzelnen Energieband taucht der Wert 1,4 als stabiler Exponent
auf, und in **jedem** Band passt das Potenzgesetz deutlich schlechter als
gepoolt. Der gepoolte Exponent ist damit eine Eigenschaft der
Energie-Ziehungsverteilung, nicht der Defektbildung.

### Warum das den Vergleich mit Sand et al. trifft

Sand et al. (2013) simulieren ~100 Kaskaden bei **fester** Energie (150 keV) und
finden dort S = 1,63 ± 0,07. Das ist genau die Größe, die hier in den
Energiebändern gemessen wird -- und dort liegt der Wert bei 1,8 bis 3,0 mit
schlechtem Fit. Der bisherige Vergleich stellt eine gepoolte Spektrumsgröße einer
Festenergie-Größe gegenüber. Die Differenz "2D vs. 3D" erklärt das nicht.

---

## 4. Was auf die Folie kann -- und was nicht

**Nicht haltbar:**
> "Die Clustergrößenverteilung folgt einem Potenzgesetz mit S = 1,4 ± σ,
> gegen Alternativen getestet."

Der Test wurde gemacht -- und ist negativ ausgefallen.

**Haltbar und methodisch stärker:**
> "Der Exponent des gepoolten Ensembles ist präzise bestimmbar
> (S = 1,41, 95-%-CI [1,39; 1,43], 1000 Bootstrap-Resamples). Der strenge
> Clauset-Test verwirft das reine Potenzgesetz jedoch (p < 0,001); Lognormal
> und Potenzgesetz-mit-Cutoff beschreiben die Daten besser. Ein Kontrolltest
> in Energiebändern zeigt die Ursache: Bei n ~ E^2,9 und log-uniform gezogener
> PKA-Energie erzeugt bereits das Sampling ein n⁻¹-Verhalten. Der Exponent
> misst die Energieverteilung, nicht die Kaskadenfragmentierung."

Das ist die ehrlichere Folie -- und die, die einer Nachfrage standhält. Sie
zeigt, dass das Verfahren beherrscht wird, statt ein Ergebnis zu behaupten,
das der erste kritische Blick zerlegt.

---

## 5. Konsequenz für die Bachelorarbeit

Der Test bei **fester** PKA-Energie ist der eigentlich interessante -- und er
ist mit dem vorhandenen Code sofort machbar: `run_one.py` so ändern, dass
`pka_energy` konstant bleibt und nur Ort/Winkel/Seed variieren. Dann misst man
die Fragmentierungsstatistik statt der Energieziehung, und der Vergleich mit
Sand et al. wird zulässig.

Zweiter offener Punkt: Die Clusterdefinition (verlagert ab `disp > 0,5·L0`,
Verknüpfungsradius `1,5·L0`) ist eine Konvention. Eine Sensitivitätsanalyse
über beide Schwellwerte gehört in die Arbeit.

---

## Dateien

| Datei | Inhalt |
|---|---|
| `clauset_refit.md` / `.json` | Vollständige Fit-Ergebnisse, LR-Tests, x_min-Scan |
| `clauset_ccdf_alternatives.png` | CCDF mit Potenzgesetz und allen vier Alternativen |
| `clauset_bootstrap_alpha.png` | Bootstrap-Verteilung des Exponenten |
| `clauset_xmin_scan.png` | KS-Abstand und S als Funktion von x_min |
| `energy_control.png` | Exponent gepoolt vs. pro Energieband |
| `../../../../docs/clauset_glossar.md` | Erklärung aller Fachbegriffe |

Reproduktion:
```bash
python3 scripts/python/statistik/cluster_from_zip.py Messungen_USW.zip ensemble/ results/statistik/ensemble_gepoolt/heal
.venv-powerlaw/bin/python scripts/python/statistik/clauset_refit.py \
    results/statistik/ensemble_gepoolt/heal/cluster_sizes.csv "Mit Healing" \
    results/statistik/ensemble_gepoolt/no_heal/cluster_sizes.csv "Ohne Healing" \
    --outdir results/statistik/fit_gepoolt --bootstrap 1000 --gof 1000
.venv-powerlaw/bin/python scripts/python/statistik/energy_control.py \
    results/statistik/ensemble_gepoolt/heal/cluster_by_run.csv "Mit Healing" \
    results/statistik/ensemble_gepoolt/no_heal/cluster_by_run.csv "Ohne Healing"
```
