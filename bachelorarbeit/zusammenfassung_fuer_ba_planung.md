# Kontext für die Planung der Bachelorarbeit

*Stand: 2026-09-30. Geschrieben als Einstieg für einen anderen KI-Chat, in dem
die Bachelorarbeit geplant wird. Sie ist in sich abgeschlossen, also ohne
Zugriff auf das Repo lesbar.*

> **Nicht an den Agenten der Bachelorarbeit weitergeben.** Abschnitt 4 listet
> Fehler und Erkenntnisse, die der Agent selbst finden soll. Kennt er sie
> vorher, ist das Experiment wertlos.

---

## 1. Worum es geht

Paul Kamm, Bachelorarbeit. Ausgangspunkt ist ein abgeschlossenes
HPC-Semesterprojekt: eine 2D-Molekulardynamik-Simulation von
Kollisionskaskaden (Strahlenschaden in Kristallen). Darauf bauen jetzt zwei
getrennte Vorhaben auf:

| | Bachelorarbeit | Paper |
|---|---|---|
| **Frage** | Kann ein agentenbetriebenes HPC-Labor ein bestehendes Minimalmodell eigenständig so weiterentwickeln, dass es MD-Referenzen näher kommt, ohne seine Einfachheit zu verlieren? (Arbeitstitel, noch nicht final) | Physik der Kaskaden: Folgt die Größenverteilung der Defektcluster einem Potenzgesetz, und unter welchen Modellannahmen? |
| **Wer arbeitet** | ein KI-Agent, möglichst autonom | Paul, mit KI-Hilfe, interaktiv |
| **Rechner** | minimaler HPC-Cluster auf Pauls Kubernetes-Cluster zu Hause | PC2 Paderborn (Noctua 2: AMD EPYC 7763 mit 128 Kernen pro Knoten, NVIDIA A100) |
| **Rolle des anderen** | Das Paper prüft mit echten Läufen, ob die Ideen und Ergebnisse des Agenten tragen | Die Arbeit liefert Ideen und Modellvarianten, die das Paper bestätigt oder widerlegt |

Betreuung/Kontakt: Dr. Maurice Maurer (TH Lübeck).

---

## 2. Das Minimalmodell (Stand Kursprojekt)

Der Code wurde **KI-generiert und nie systematisch validiert**. Paul hat das
selbst so eingeordnet. Die Modellannahmen sind frei gesetzt oder ohne bewusste
Entscheidung entstanden.

**Physik:**
- 2D-Dreiecksgitter, jedes Atom hat 6 Nachbarn. Gleiche Massen m = 1, Abstand L0 = 1.
- Nachbarn sind durch harmonische Federn verbunden, K = 100.
- Eine Feder reißt, wenn sie über 1,15·L0 gedehnt wird (`MAX_STRETCH`).
  Ursprünglich war das als Stellschraube für die Duktilität verschiedener
  Materialien gedacht.
- Abstoßung auf kurze Distanz: V ~ r⁻¹² mit Reichweite 0,9·L0 und Stärke 400,
  glatt auf 0 geführt.
- „Healing“: Gerissene Bindungen schließen sich wieder, wenn zwei Atome näher
  als 1,10·L0 und langsamer als 0,5 relativ zueinander sind. Die Schwellen
  sind frei gesetzt.
- Ein „Primary Knock-on Atom“ (PKA) bekommt eine Startenergie und löst die
  Kaskade aus. Am Ende bleiben verlagerte Atome und gerissene Bindungen.
- Velocity-Verlet mit festem Zeitschritt (2–3·10⁻⁴). Die Nachbarsuche läuft
  über eine Zellliste.
- Alles ist dimensionslos, es gibt keine Kalibrierung auf ein reales Material.
- Konfiguration über eine `.ini`-Datei, die alle drei Implementierungen lesen.

**Implementierungen:**
- seriell in C (~1200 Zeilen),
- MPI mit 1D-Streifenzerlegung (~1100 Zeilen),
- CUDA mit Benchmark- und Optimierungsvarianten (~2000 Zeilen).

Dazu kommen SLURM-Jobskripte für Noctua 2, Python-Auswertung, eine CI mit
Build, Physik-Smoke-Test und Valgrind, und Visualisierung mit OVITO.

**Ergebnisse des Kursprojekts (Bericht):**
- **HPC:** Der Code ist durchgehend durch die Speicherbandbreite begrenzt.
  - Eine A100 ist 731× schneller als ein EPYC-Kern.
  - Die MPI-Skalierung verläuft U-förmig, mit superlinearer Effizienz (133 %
    bei 128 Kernen) durch Cache-Effekte.
  - Optimierungen der Kommunikation bringen nichts.
  - Diese Ergebnisse gelten weiter.
- **Physik:** Die Clustergrößen aus 1000 Kaskaden folgen angeblich einem
  Potenzgesetz p(n) ~ n^−S mit S ≈ 1,4. Der Bericht verglich das mit Sand et
  al. 2013 (3D-Wolfram, S = 1,63). **Diese Aussage ist inzwischen widerlegt**,
  siehe Abschnitt 4.

Eine Physik-Präsentation mit Videos für Maurices Professor ist gehalten und
abgeschlossen.

---

## 3. Was nach dem Kursprojekt passiert ist (Aug.–Sept. 2026)

Mit KI-Hilfe interaktiv gearbeitet, nicht autonom:

1. **Nachanalyse der Clusterstatistik** mit der Methode von Clauset, Shalizi
   und Newman (2009): diskreter Maximum-Likelihood-Schätzer, Bootstrap,
   Goodness-of-Fit, Likelihood-Ratio-Tests.
2. **Neues Ensemble bei fester Energie:** 900 Läufe, 250×250 Atome,
   3 Energien, gepaartes Design.
3. **ZBL-Umbau der seriellen Version,** begonnen, nicht abgeschlossen:
   - Die r⁻¹²-Abstoßung wird durch das ZBL-Potenzial ersetzt, also abgeschirmte
     Coulomb-Abstoßung, Standard in der Kaskaden-MD. Damit schrumpft der
     Wirkungsquerschnitt mit der Energie, wie in der Realität.
   - Dazu kommen ein adaptiver Zeitschritt, elektronische Bremsung nach
     Lindhard-Scharff, ein absorbierender Rand und ein Frühabbruch.
   - Die Einheiten sollen auf Wolfram kalibriert werden
     (L0 = 2,74 Å, Verlagerungsschwelle E_d = 90 eV als Anker).
   - Das alte Modell bleibt per Schalter bitgleich erhalten.
   - Stand: Die Kalibrierung hängt an einer Modellentscheidung (Abschnitt 4,
     Punkt 5). Dieser Umbau gehört ab jetzt zum **Paper**.

---

## 4. Bekannte Fehler und Erkenntnisse = Prüfliste für den Agenten

Mensch und KI haben das interaktiv gefunden. Für die Bachelorarbeit ist es eine
fertige Referenz: Findet der Agent dasselbe selbst?

| # | Befund | Art | Wie gefunden |
|---|---|---|---|
| 1 | **Energie-Sampling erzeugt ein Scheinpotenzgesetz.** Die PKA-Energie wurde log-uniform aus [200, 2500] gezogen, und der größte Cluster wächst wie E^2,9. Daraus folgt rein mathematisch p(n) ~ n^−1. Das S ≈ 1,4 misst die Stichprobe, nicht die Physik. In einzelnen Energiebändern taucht 1,4 nicht auf. | methodischer Fehler | Fit-Tests (GoF p = 0, Lognormal passt besser), dann Zusammenhang zwischen Energie und Clustergröße |
| 2 | **Kontinuierlicher statt diskreter Schätzer:** Der Bericht liegt ~0,02 zu hoch, der Fehlerbalken war zu groß. | kleiner methodischer Fehler | Nachrechnen |
| 3 | **Die Clusterdefinition bestimmt das Ergebnis.** Mit der alten Definition (verlagert ab 0,5·L0, Linkradius 1,5·L0) verschmilzt der Schaden zu einem Klumpen, 92 % der Atome liegen im größten Cluster. Mit feinerer Definition (1,0/1,0) gibt es 30–70 Cluster pro Lauf. Über sechs Definitionen läuft S von 1,23 bis 3,30 auf denselben Daten. Auch die feine Definition ist statistisch kein Potenzgesetz. Es gibt nur einen Trend dorthin mit wachsender Kaskadengröße (S = 2,50 → 2,19 → 2,00). | methodische Erkenntnis | Analyse bei fester Energie |
| 4 | **Gebundene Nachbarn sind durchlässig.** Die Abstoßung wirkt nur zwischen *ungebundenen* Paaren, gebundene spüren nur ihre Feder, und die liefert höchstens ½·K·L0² = 50. Ab E ≈ 200 fliegt ein Atom durch seinen Nachbarn hindurch, schräg ganz ohne Ablenkung. Die erste Nachbarschale ist für schnelle Atome unsichtbar. Die Energie bleibt dabei erhalten, deshalb hat kein Test es bemerkt. | **Modellfehler im Code** | gezielte Diagnose (Mindestabstand PKA–Nachbar gegen E) |
| 5 | **Die Bindungen sind 5–13× zu schwach für ein reales Metall.** Kohäsionsenergie 3,4 Modelleinheiten, je nach Kalibrierung 0,7–2,0 eV, Wolfram hat 8,9 eV. In der Kaskade schmilzt das Gitter großflächig. | Modellgrenze | Kalibrierung auf Wolfram |
| 6 | **Der Wirkungsquerschnitt hängt kaum von der Energie ab.** Mindestabstand beim zentralen Stoß 0,64·L0 bei E = 2400, noch 0,49·L0 bei 50 000. Real (W, 150 keV) sind es 0,07 Nachbarabstände. Schnelle Atome können deshalb nicht weit fliegen und in Subkaskaden zerfallen, dem wichtigsten Kandidaten für ein Potenzgesetz. | Modellgrenze | Zweikörper-Rechnung |
| 7 | **Größere Gitter allein bringen nichts.** Bei 250² kommt keine Kaskade in die Nähe des Rands (Ausdehnung 17–35 L0). | Erkenntnis | Messung |
| 8 | **Keine elektronische Bremsung, keine Einheiten.** Die Bremsung ist im Code vorhanden, war aber in allen Läufen aus. „E = 2400“ lässt sich nicht in eV übersetzen. | Modellgrenze | Code-Lektüre |
| 9 | **Werkzeugfalle:** Das Python-Paket `powerlaw` klemmt den Exponenten stillschweigend bei 3,000. | Werkzeugfehler | auffällige Werte |
| 10 | Kommentare zu den Einheiten in der `.ini` sind irreführend („Zeiteinheit ~0,1 s“). | Doku-Fehler | Lektüre |

Referenz für die ZBL-Richtung (ein möglicher „Zielpunkt“ des Agenten): Der
Mensch-plus-KI-Umbau aus Abschnitt 3 mit Tests: Zweikörperstreuung gegen das
klassische Ablenkintegral, Energieerhaltung auf 10⁻⁴, Byte-Gleichheit des
alten Modells.

---

## 5. Was der Agent bekommt (beschlossen)

- **Ein eigenes Repo**, getrennt vom Repo des Kursprojekts und des Papers. Er
  sieht weder unsere Doku noch die Befunde noch den ZBL-Umbau.
- **Startpunkt: das Minimalmodell genau wie im Kursprojekt, mit allen
  Fehlern.** Die Fehler aus Abschnitt 4 sind absichtlich nicht behoben. Sie
  dienen als Messinstrument.
- **Rechner:** ein minimaler HPC-Cluster auf Pauls Kubernetes-Cluster zu Hause.

Noch **nicht** festgelegt, und Thema für die Planung:
- Welche MD-Referenzen er bekommt und wie „näher kommen“ gemessen wird.
  Kandidaten: Verlagerungsschwelle, Defektzahl über Energie (NRT/arc-dpa),
  Clustergrößenverteilung (Sand et al. 2013), Schallgeschwindigkeit,
  Kohäsionsenergie.
- Wie „ohne seine Einfachheit zu verlieren“ gemessen wird, etwa über Zahl der
  Parameter, Codeumfang oder Laufzeit pro Atom und Schritt.
- Welche Code-Teile er bekommt: nur seriell, oder auch MPI und CUDA?
- Aufgabenformulierung, Autonomiegrad, Budget, Protokollierung seiner
  Entscheidungen, Wiederholbarkeit (mehrere Agentenläufe?).
- Ob und wie der Bericht des Kursprojekts (mit der widerlegten S ≈ 1,4-Aussage)
  Teil seiner Eingabe ist.

## 6. Was Paul im Paper macht

- Echte Läufe auf dem PC2 mit großen Gittern und hohen Energien.
- (a) Ideen und Ergebnisse des Agenten mit echten Läufen bestätigen oder
  widerlegen.
- (b) Eigene Ideen umsetzen. Wichtigste: der energieabhängige
  Wirkungsquerschnitt, also der ZBL-Umbau aus Abschnitt 3.
  - Nächster Schritt dort: entscheiden, ob die Bruchdehnung angehoben wird,
    damit die Kohäsionsenergie zu Wolfram passt (Empfehlung: ≈ 1,35, Energieskala
    aus der Schallgeschwindigkeit, ≈ 0,44 eV pro Einheit).
  - Danach Pilotläufe, Produktionsensemble, Clauset-Auswertung.
- Physikalische Kernfrage: Wird die Clustergrößenverteilung mit realistischem
  Querschnitt, Subkaskaden und großen Gittern zu einem Potenzgesetz? Wie
  verhält sich S in 2D zu 3D (Sand et al.: 1,63)?

## 7. Organisatorisches

- Git-Repo auf Pauls Forgejo, GitHub-Mirror. Sprache Deutsch.
- Keine Attribution an KI-Werkzeuge in Commits, Code oder Dokumenten.
- Geplante Neuordnung dieses Repos: `kursprojekt/` (eingefroren, mit
  Clusterstatistik und Präsentation), `paper/` (ZBL-Weiterentwicklung,
  PC2-Läufe), `bachelorarbeit/` (Konzept und Auswertung des Agentenlabors,
  nicht das Arbeitsrepo des Agenten), `docs/`.
- Literatur: Sand, Dudarev, Nordlund 2013 (EPL 103, 46003; Clustergrößen in
  W-Kaskaden); Clauset, Shalizi, Newman 2009 (SIAM Review, Potenzgesetz-Fits);
  De Backer et al. 2016 (Subkaskaden); Ziegler, Biersack, Littmark (ZBL);
  Lindhard-Scharff (elektronische Bremsung).
