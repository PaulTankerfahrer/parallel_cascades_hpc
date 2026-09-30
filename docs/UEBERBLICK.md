# Überblick: Was ist was, und warum

*Stand: 2026-09-24. Einstieg ins Projekt. Alles Weitere ist von hier verlinkt.*

Zwei Begleitdateien:
- [`ENTSCHEIDUNGEN.md`](ENTSCHEIDUNGEN.md): Welche Entscheidung wann und warum
  gefallen ist, welche Alternativen es gab, was sie an alten Ergebnissen ändert.
- [`PARAMETER.md`](PARAMETER.md): Jeder Modellparameter mit Wert, Bedeutung,
  Herkunft und Status.

---

## Die Frage in einem Satz

Folgt die Größenverteilung der Defektcluster, die eine Kollisionskaskade in
einem 2D-Gitter hinterlässt, einem Potenzgesetz p(n) ~ n^−S, und wenn ja, mit
welchem S?

## Drei Projektteile, die im selben Repo liegen

| Teil | Zeitraum | Worum es ging | Wo es liegt | Status |
|---|---|---|---|---|
| **1. HPC-Semesterprojekt** | bis 17.08. | Dieselbe Simulation seriell, mit MPI und mit CUDA. Hauptthema war die Laufzeit (Skalierung, GPU-Optimierung). Nebenbei ein erstes Potenzgesetz-Ergebnis, S ≈ 1,4. | `report.pdf`/`report.tex`, `src/mpi_cascade.c`, `src/cascade_cuda*.cu`, `scripts/job_*.sh`, CSVs direkt in `results/`, [`results/KATALOG.md`](../results/KATALOG.md) | abgeschlossen |
| **2. Clusterstatistik** | bis 02.09. | Prüfen, ob S ≈ 1,4 stimmt. Ergebnis: nein. Dann ein neues Ensemble bei fester Energie. Dazu die Physik-Präsentation. | `results/statistik/`, `scripts/python/statistik/`, `presentation/` | abgeschlossen, Präsentation gehalten und fertig |
| **3. ZBL-Umbau** | seit 24.09. | Das Modell physikalisch realistischer machen: energieabhängiger Wirkungsquerschnitt, Wolfram-Einheiten, größere Gitter. Dann die Frage aus Teil 2 neu stellen. | `src/cascade_serial.c`, `src/potential.h`, `scripts/python/zbl/`, `results/zbl/`, [`plan_zbl.md`](plan_zbl.md) | **in Arbeit, wartet auf eine Entscheidung** |

Wichtig für die Orientierung:
- `src/params.ini` ist die **Benchmark-Konfiguration aus Teil 1** (1414², kein
  Healing). Die Physik-Läufe aus Teil 2 und 3 schreiben ihre eigene `.ini` pro
  Lauf aus den Python-Skripten. Die Datei dokumentiert aber alle Schlüssel.
- Nur `cascade_serial.c` kennt den Umbau aus Teil 3. MPI und CUDA rechnen
  weiter das alte Modell. Das alte Modell ist in `cascade_serial.c` als
  `rep_model = r12` bitgleich erhalten, der Standard.

---

## Die Geschichte als Kette

Jeder Schritt folgt aus dem vorigen.

1. **Behauptung aus dem Bericht:** S = 1,41 (mit Healing) bzw. 1,38 (ohne),
   aus 1000 Kaskaden, Vergleich mit Sand et al. (3D-Wolfram, S = 1,63).

2. **Nachgeprüft (Teil 2): Die Zahl ist reproduzierbar, aber kein
   Potenzgesetz.** Der Goodness-of-Fit-Test verwirft sie (p = 0,000), eine
   Lognormalverteilung passt besser.
   → [`BEFUND_gepoolt.md`](../results/statistik/fit_gepoolt/BEFUND_gepoolt.md)

3. **Ursache gefunden: Das Energie-Sampling erzeugt das Potenzgesetz.** Die
   PKA-Energie wurde log-uniform gewürfelt, und der größte Cluster wächst wie
   E^2,9. Daraus folgt rein rechnerisch p(n) ~ n^−1, ohne jede Kaskadenphysik.
   Bei fester Energie taucht die 1,4 nicht mehr auf.

4. **Deshalb neues Ensemble bei fester Energie** (900 Läufe, E = 1300 / 1800 /
   2400, 250²). Zweites Problem: Mit der alten Clusterdefinition verschmilzt
   der ganze Schaden zu einem Klumpen, da gibt es keine Verteilung.

5. **Feinere Clusterdefinition:** 30–70 getrennte Cluster pro Lauf, S = 2,50 →
   2,19 → 2,00 mit steigender Energie. Aber: S hängt stark von der Definition ab
   (1,23 bis 3,30 auf denselben Daten), und auch die feine Definition ist
   statistisch **kein** Potenzgesetz. Belastbar ist nur der Trend: Je größer die
   Kaskade, desto näher am Potenzgesetz.
   → [`BEFUND_fest.md`](../results/statistik/fit_fest/BEFUND_fest.md)

6. **Neue Idee (24.09.):** Realistisch wäre ein Wirkungsradius, der mit der
   Geschwindigkeit schrumpft, und größere Gitter. Gemessen:
   - Im alten Modell schrumpft der Radius kaum (0,64 L0 bei E = 2400, noch
     0,49 L0 bei 50 000; real 0,07 Nachbarabstände).
   - Keine Kaskade kommt bei 250² auch nur in die Nähe des Rands. Größere Gitter
     allein bringen also nichts, erst zusammen mit einem Modell, in dem Atome
     weit fliegen.

   Der physikalische Gedanke: Schnelle Atome fliegen weit, stoßen selten hart
   und zerfallen in Subkaskaden verschiedener Größe. Das gilt als wichtigster
   Kandidat für ein Potenzgesetz, und das alte Modell kann es nicht.

7. **Plan: ZBL-Potenzial**, das Standardpotenzial für kurze Abstände in der
   Kaskaden-MD, dazu Wolfram-Einheiten und elektronische Bremsung.
   → [`plan_zbl.md`](plan_zbl.md)

8. **Beim Vorbereiten gefunden: ein Modellfehler.** Im alten Modell spüren
   gebundene Nachbarn keine Abstoßung, nur ihre Feder, und die gibt höchstens
   50 Energieeinheiten her. Ab E ≈ 200 fliegt ein Atom einfach durch seinen
   Nachbarn hindurch, schräg sogar ohne jede Ablenkung.
   → [`befund_gebundene_nachbarn.md`](befund_gebundene_nachbarn.md)

9. **Umbau gebaut und getestet:** ZBL für alle Paare, adaptiver Zeitschritt,
   elektronische Bremsung, absorbierender Rand, Frühabbruch. Alle 8 CI-Stufen
   laufen durch.

10. **Kalibrierung hängt (Stand jetzt).** Um „E = 2400“ in eV zu übersetzen,
    braucht es den Umrechnungsfaktor ε. Dabei zeigt sich: Die Bindungen des
    Federmodells sind 5–13× schwächer als in Wolfram. Das Gitter schmilzt in
    der Kaskade großflächig. Entscheidung offen, siehe
    [`haltepunkt_2026-09-24.md`](haltepunkt_2026-09-24.md) und
    [`ENTSCHEIDUNGEN.md`](ENTSCHEIDUNGEN.md) E17.

---

## Was von den Ergebnissen noch gilt

| Ergebnis | Gilt noch? | Warum |
|---|---|---|
| Bericht: S ≈ 1,4, Potenzgesetz | **nein** | Artefakt des Energie-Samplings (Schritt 3) |
| Sampling-Mechanismus: log-uniforme Energie erzeugt n^−1 | **ja** | reine Mathematik, hängt nicht am Modell |
| S hängt stark von der Clusterdefinition ab | **ja**, als methodischer Befund | Die Zahlenwerte stammen aus dem alten Modell |
| S = 2,50 → 2,00, Trend zum Potenzgesetz | **nur für das alte Modell** | Das Modell hat durchlässige Bindungen (Schritt 8) |
| Sand-Vergleich (2D ≈ 2,0 über 3D 1,63) | **nicht belastbar** | Altes Modell, keine physikalischen Einheiten |
| „Größeres Gitter allein hilft nicht“ | **ja, für das alte Modell** | Mit ZBL werden die Kaskaden größer, dann schon |
| HPC-Ergebnisse (Skalierung, GPU) | **ja** | Laufzeitmessungen, von der Stoßphysik unabhängig |

---

## Wo steht was

| Ich will wissen … | Datei |
|---|---|
| … wo wir gerade stehen und was als Nächstes ansteht | [`haltepunkt_2026-09-24.md`](haltepunkt_2026-09-24.md) |
| … warum etwas so ist, wie es ist | [`ENTSCHEIDUNGEN.md`](ENTSCHEIDUNGEN.md) |
| … was ein Parameter bedeutet und woher sein Wert kommt | [`PARAMETER.md`](PARAMETER.md) |
| … wie der Code aufgebaut ist, wie man baut und startet | [`../README.md`](../README.md) |
| … was die Fachbegriffe der Potenzgesetz-Statistik bedeuten | [`clauset_glossar.md`](clauset_glossar.md) |
| … alle Zahlen aus Teil 2 kompakt | [`projektstand.md`](projektstand.md) |
| … Details zum alten Sampling-Befund | [`BEFUND_gepoolt.md`](../results/statistik/fit_gepoolt/BEFUND_gepoolt.md) |
| … Details zum Festenergie-Ensemble | [`BEFUND_fest.md`](../results/statistik/fit_fest/BEFUND_fest.md) |
| … den Modellfehler mit den Bindungen | [`befund_gebundene_nachbarn.md`](befund_gebundene_nachbarn.md) |
| … den Plan und die Tests des Umbaus | [`plan_zbl.md`](plan_zbl.md) |
| … die Kalibrierungszahlen | [`../results/zbl/kalibrierung/ergebnis.md`](../results/zbl/kalibrierung/ergebnis.md) |
| … was die HPC-Messdateien in `results/` sind | [`../results/KATALOG.md`](../results/KATALOG.md) |

---

## Offene Punkte, nach Dringlichkeit

1. **Entscheidung Bindungsstärke** (E17 in `ENTSCHEIDUNGEN.md`). Davon hängt
   alles Weitere im Umbau ab.
2. ~~Präsentation~~: gehalten und fertig (30.09.). Bekannte Ungenauigkeiten
   (Kraftkurve `P2_kraft.png` mit Abstoßung auch für gebundene Nachbarn,
   Werte aus dem alten Modell) bleiben so, wie sie gezeigt wurden.
3. Umbau fertigstellen: Pilotläufe, Produktion (700²), Auswertung
   (Phasen 4–6 in `plan_zbl.md`).
4. Für die Arbeit: Sensitivität über die Clusterdefinition als eigenes Kapitel.
