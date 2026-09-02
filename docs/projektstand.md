# Projektstand: Kollisionskaskaden-HPC — Kontext für einen neuen Chat

*Stand: 2026-09-02. Zum Einfügen als Kontext am Anfang einer neuen Sitzung.*

---

## Kontext

Repo: `/home/paul/Dokumente/parallel_cascades_hpc` (git, Branch `main`, Remote Forgejo
+ GitHub-Mirror). Ursprung ist ein abgeschlossenes HPC-Semesterprojekt
(`report.pdf`, `report.tex`): 2D-Molekulardynamik von Kollisionskaskaden,
seriell / MPI / CUDA. Darauf baut jetzt eine **Bachelorarbeit** auf.

Sprache: Deutsch. Keine Attribution an KI-Werkzeuge in Commits, Code oder Dokumenten.

Nebenprodukte aus dem Semesterprojekt, die weiterlaufen:
- Beamer-Präsentation `presentation/physik_praesentation.tex` (9 Seiten) mit vier
  OVITO-Videos. Headless-Rendering braucht das PyPI-Paket `ovito`, nicht das GUI-Paket.
- CI in drei Workflows (Build/Physik-Smoke-Test, Valgrind, LaTeX-PDFs), gespiegelt in
  `.forgejo/workflows/` und `.github/workflows/`. Valgrind schlägt **lokal** auf neuer
  glibc fehl, in CI läuft es sauber.

---

## Die zentrale physikalische Erkenntnis

Der Bericht behauptet: Clustergrößen folgen einem Potenzgesetz `p(n) ~ n^-S` mit
S ≈ 1,4. Zwei aufeinander aufbauende Analysen haben das zerlegt und neu aufgebaut.

### Teil 1 — Das alte Ergebnis misst das Sampling, nicht die Physik

Datenbasis: die Originaldaten, je 1000 Kaskaden mit und ohne Healing (liegen als Zip
in der Nextcloud, **nicht** entpackt im Repo; `cluster_from_zip.py` liest direkt daraus).
Methodik: Clauset, Shalizi & Newman (2009), Paket `powerlaw` 2.0.0 in `.venv-powerlaw`.

- **Reproduzierbar:** S = 1,391 ± 0,012 (heal) / 1,361 ± 0,011 (no_heal). Die
  Berichtswerte 1,41/1,38 entstehen exakt mit dem *kontinuierlichen* MLE; Clustergrößen
  sind ganzzahlig, der korrekte *diskrete* Schätzer liegt ~0,02 tiefer. Bekannter Bias,
  kein Fehler. Der angegebene Fehler ±0,03 war zu konservativ.
- **Aber verworfen:** Goodness-of-Fit p = 0,000. Im Likelihood-Ratio-Test verliert das
  Potenzgesetz gegen Lognormal (R = −9,5), Stretched Exponential (−2,4) und
  Potenzgesetz-mit-Cutoff (−45,0). Es schlägt nur die Exponentialverteilung.
- **Ursache:** `run_one.py` zieht die PKA-Energie **log-uniform** aus [200, 2500]. Der
  größte Cluster skaliert wie `n_max ~ E^2,9` (r = 0,97). Ist E log-uniform und gilt
  `n ~ E^k`, dann ist ln(n) gleichverteilt — und das heißt `p(n) ~ n^-1`, unabhängig
  von k und von jeder Kaskadenphysik.
- **Kontrolltest:** In drei Energiebändern statt gepoolt taucht die 1,4 in **keinem**
  Band auf (dort 1,8–3,0) und der Fit wird in **jedem** Band schlechter
  (KS-D 0,09–0,19 gegen 0,041 gepoolt).

Folge: Der Vergleich mit Sand et al. (2013, 3D-Wolfram, ~100 Kaskaden bei **fester**
Energie 150 keV, S = 1,63 ± 0,07) war unzulässig — gepoolte Spektrumsgröße gegen
Festenergie-Größe.

### Teil 2 — Ensemble bei fester Energie

900 Läufe (3 × 300 bei E = 1300 / 1800 / 2400), NX = 250, `healing = true`,
dt = 3e-4, 12 000 Schritte, ~35 min auf 8 Kernen. Einschlagort, Winkel und Seed hängen
**nur von der Lauf-Nummer** ab, nicht von der Energie — der Vergleich über die Energien
ist gepaart, der Energieeffekt isoliert.

**a) Die alte Clusterdefinition trägt nicht.** Mit `disp > 0,5·L0` / Linkradius `1,5·L0`
kommt S = 1,40 / 1,35 / 1,23 heraus, aber KS-D = 0,173 / 0,201 / 0,300 — vier- bis
siebenmal schlechter als gepoolt. Der Grund ist geometrisch: die ganze Schadenszone
verschmilzt zu **einem** Klumpen. 92 % aller verlagerten Atome liegen im größten
Cluster; 365 der 1000 Altläufe haben überhaupt nur einen. Die „Verteilung" ist
zweigipflig — Einer und Zweier plus ein Riesenblob.

**b) Eine feinere Definition liefert eine echte Verteilung.** `disp > 1,0·L0` /
Linkradius `1,0·L0`, dieselben Läufe, dieselben Atome:

| E | Cluster | pro Lauf | größter | x_min | S (95 %-CI) | KS-D |
|---|--:|--:|--:|--:|--:|--:|
| 1300 | 9 281 | 30,9 | 43 | 2 | 2,498 [2,230, 2,950] | 0,032 |
| 1800 | 14 588 | 48,6 | 113 | 2 | 2,191 [2,076, 2,279] | 0,017 |
| 2400 | 21 580 | 71,9 | 330 | 1 | 1,996 [1,982, 2,010] | 0,013 |

30–70 getrennte Cluster pro Lauf, KS-Abstand 10–20-fach besser. **Das sind die Objekte,
die Sand et al. zählen** — getrennte Defektcluster, nicht die Ausdehnung der Schadenszone.

LR-Tests bei E = 2400: Exponential +41,6 und Stretched Exponential +12,4 (Potenzgesetz
gewinnt klar), Lognormal +1,8 (p = 0,072, gewinnt nicht mehr), Cutoff −4,3 — der Cutoff
liegt aber bei n ≈ 2946, fast eine Größenordnung jenseits des größten beobachteten
Clusters (330), im Datenbereich also wirkungslos. Bei E = 1300/1800 liegt der Cutoff
dagegen bei n ≈ 13 bzw. 71, mitten in den Daten.

**c) S ist keine Konstante des Modells.** Dieselben 900 Läufe, sechs Clusterdefinitionen:

| Definition | E=1300 | E=1800 | E=2400 |
|---|--:|--:|--:|
| `disp>0,5` / link 1,5 | 1,40 | 1,35 | 1,23 |
| `disp>1,0` / link 1,5 | 1,57 | 1,50 | 1,45 |
| `disp>0,5` / link 1,0 | 1,98 | 2,10 | 2,25 |
| `disp>1,0` / link 1,0 | 2,50 | 2,19 | 2,00 |
| `disp>1,5` / link 1,0 | 2,78 | 2,51 | 2,28 |
| `disp>2,0` / link 1,0 | 3,30 | 3,02 | 2,77 |

S läuft von **1,23 bis 3,30**. Die Konvention entscheidet stärker als die Physik. Eine
Zahl „S = 1,4" ohne Angabe der Clusterdefinition ist bedeutungslos.

**d) Der ehrliche Vorbehalt.** Der GoF-Test verwirft auch die feine Definition
(p = 0,000). Entlastend: bei n = 21 580 verwirft der Test schon winzige Abweichungen
(D_beob = 0,0127 gegen D_synth = 0,0024). Belastend: die CCDF ist sichtbar gekrümmt —
das ist Form, nicht Rauschen. Der brauchbare Befund ist der **Trend**: S und KS-D fallen
gemeinsam und monoton mit der Kaskadengröße (474 → 1399 → 2646 verlagerte Atome pro Lauf
ergeben S = 2,50 → 2,19 → 2,00 und D = 0,032 → 0,017 → 0,013). Die Verteilung nähert
sich einem Potenzgesetz, hat es bei diesen Kaskadengrößen aber nicht erreicht; S ist
nicht konvergiert.

**e) Sand-Vergleich.** Wieder zulässig (feste Energie + fragmentierende
Clusterdefinition), fällt aber umgekehrt aus: 2D ≈ 2,0 liegt **über** dem 3D-Wert 1,63 —
konsistent mit der geringeren Konnektivität in zwei Dimensionen.

---

## Zwei Werkzeugfallen

- **`powerlaw` klemmt `alpha` bei 3.** `DEFAULT_PARAMETER_RANGES = {'alpha': [0, 3]}`
  in `distributions.py`; steilere Verteilungen werden kommentarlos als exakt `3.000`
  gemeldet. Der Kwarg heißt `parameter_ranges` (nicht `parameter_range`), hebt die
  Klemme aber **nicht** auf. Lösung: `mle_unbounded()` in `clauset_refit.py` und
  `analyze_fixed_energy.py` rechnet jeden Wert nahe 3,000 frei nach, markiert mit ⚠.
- **`truncated_power_law` konvergiert nicht immer** und meldete R = +757 gegen das
  Potenzgesetz. Für ein genestetes Modell unmöglich. Solche Zellen stehen auf „n. k.".

---

## Dateien

| Pfad | Inhalt |
|---|---|
| `results/statistik/fit_gepoolt/BEFUND_gepoolt.md` | Teil 1: Refit der Originaldaten, Energie-Kontrolltest |
| `results/statistik/fit_fest/BEFUND_fest.md` | Teil 2: Fixed-Energy-Ensemble, 8 Abschnitte |
| `results/statistik/fit_fest/fixed_energy_uebersicht.md` | Alle Energien × Clusterdefinitionen |
| `results/statistik/fit_fest/fixed_energy_ccdf.png` | CCDF alte vs. feine Definition, drei Energien |
| `docs/clauset_glossar.md` | Alle Fachbegriffe erklärt (9 Kapitel) |
| `results/statistik/ensemble_fest/E*/` | Rohdaten der 900 Läufe (`disp/` ist gitignored) |
| `scripts/python/statistik/run_fixed_energy.py` | Ensemble-Treiber |
| `scripts/python/statistik/clauset_refit.py` | Voller Clauset-Fit: Bootstrap, GoF, LR, x_min-Scan |
| `scripts/python/statistik/analyze_fixed_energy.py` | Schnellübersicht über alle Definitionen |
| `scripts/python/statistik/fig_fixed_energy_ccdf.py` | Die CCDF-Vergleichsabbildung |
| `scripts/python/statistik/energy_control.py` | Fit gepoolt vs. pro Energieband |
| `scripts/python/statistik/cluster_from_zip.py` | Clustering direkt aus dem Nextcloud-Zip |
| `scripts/python/statistik/fig_csv_to_cluster.py` | Präsentationsabbildung CSV → Cluster |

Reproduktion:
```bash
python3 scripts/python/statistik/run_fixed_energy.py \
    --outdir results/statistik/ensemble_fest --energies 1300,1800,2400 --runs 300 --jobs 8

.venv-powerlaw/bin/python scripts/python/statistik/analyze_fixed_energy.py \
    results/statistik/ensemble_fest --outdir results/statistik/fit_fest

.venv-powerlaw/bin/python scripts/python/statistik/clauset_refit.py \
    results/statistik/ensemble_fest/E2400/clusters_thr100_link100.csv "E=2400, feine Definition" \
    --outdir results/statistik/fit_fest --bootstrap 300 --gof 200 --jobs 8
```

---

## Offene Punkte

1. **Folie 9** (`Das physikalische Ergebnis: ein Potenzgesetz`) ist an drei Stellen
   widerlegt und noch nicht angepasst: „S ≈ 1,4, robust, 1,36 vs. 1,39"; „skalenfrei /
   selbstorganisierte Kritikalität"; „Sand: 1,63 — der 2D-Wert liegt konsistent
   darunter" (falsch in Zahl **und** Richtung). Eine fertige Ersatzformulierung steht in
   §6 von `BEFUND_fest.md`.
2. **Nichts vom Physik-Teil ist committet.** Letzter Commit `ad5be04` (17.08.), alles
   danach liegt unversioniert im Arbeitsverzeichnis.
3. **Für die Arbeit:** größere Box und höhere Energie, um zu prüfen, ob S konvergiert.
   Der Trend legt nahe, dass die 2D-Kaskade bei 250×250 noch nicht im asymptotischen
   Regime ist.
4. **Für die Arbeit:** Sensitivitätsanalyse über die Clusterdefinition — durch die
   Sechs-Definitionen-Tabelle bereits halb geschrieben.
