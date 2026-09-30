# Entscheidungslog

*Stand: 2026-09-30. Zurück zum [Überblick](UEBERBLICK.md).*

Hier steht jede Entscheidung, die das Modell oder die Auswertung prägt, mit
Begründung, verworfenen Alternativen und Folgen.

**Wichtig für Teil 1:** Der Code des Semesterprojekts ist KI-generiert und wurde
nie systematisch validiert. Seine Modellannahmen (S1–S5) sind entweder frei
gesetzt oder ohne bewusste Entscheidung entstanden. Sie sind ein Ausgangspunkt,
keine begründeten Modellwahlen.

Status: **gilt** · **ersetzt** (durch spätere Entscheidung) · **offen**

---

## Teil 1: HPC-Semesterprojekt (bis 17.08.)

Diese Entscheidungen stammen aus dem Semesterprojekt. Sie waren auf Laufzeit
ausgerichtet, nicht auf physikalischen Realismus.

### S1 Dreiecksgitter mit harmonischen Federn · gilt
- **Was:** 2D-Dreiecksgitter, 6 Nachbarn, Federn mit K = 100 und L0 = 1.
- **Warum:** Ein Quadratgitter mit Federn nur entlang der Bindungen ist nicht
  schubstabil, das Dreiecksgitter schon (Kommentar in `params.ini`).
- **Folge:** Das Modell ist ein Minimalmodell und kein Wolfram-Potenzial.
  Kohäsion und Thermodynamik sind deshalb nicht realistisch (siehe E17).

### S2 Bindungsbruch bei 15 % Dehnung (MAX_STRETCH = 1,15) · gilt, steht zur Diskussion
- **Warum:** frei gesetzt (Paul, 30.09.). Gedacht als Stellschraube, um
  verschiedene Materialien und ihre Duktilität zu simulieren. Kein Bezug zu
  einem realen Material.
- **Folge:** Legt fest, wie stark das Gitter zusammenhält. Pro Bindung
  ½·K·(0,15)² = 1,125 Einheiten, pro Atom 3,4 Einheiten. Das ist der Kern von E17.

### S3 Abstoßung r⁻¹² nur zwischen ungebundenen Paaren · ersetzt durch E10
- **Was:** Die Abstoßung wird für gebundene Paare übersprungen.
- **Warum:** keine bewusste Entscheidung. Der Code wurde KI-generiert und nie
  validiert (Paul, 30.09.). Sehr wahrscheinlich ein Fehler.
- **Folge:** Gebundene Nachbarn sind für schnelle Atome durchlässig (E09).

### S4 PKA-Energie log-uniform aus [200, 2500] · ersetzt durch E02
- **Warum:** keine bewusste Entscheidung, vom KI-generierten Code übernommen
  (Paul, 30.09.). Ein tieferer Grund ist nicht bekannt.
- **Folge:** Genau diese Wahl erzeugt das scheinbare Potenzgesetz mit S ≈ 1,4 (E02).

### S5 Healing ein und aus als zwei Varianten · gilt
- **Was:** Gerissene Bindungen können sich wieder schließen, wenn zwei Atome
  nah und langsam genug sind (`healing_dist` = 1,10, `healing_vrel` = 0,5).
- **Warum:** Beide Varianten wurden als Robustheitstest gerechnet und ergaben
  fast gleich viele Cluster (2391 gegen 2335). Die Schwellenwerte selbst sind
  frei gesetzt, wie S2 (Paul, 30.09.).
- **Folge:** Ab Teil 2 wird immer mit Healing gerechnet.

---

## Teil 2: Clusterstatistik (bis 02.09.)

### E01 Diskreter Maximum-Likelihood-Schätzer nach Clauset · gilt
- **Was:** S wird mit dem diskreten MLE geschätzt, dazu Bootstrap,
  Goodness-of-Fit und Likelihood-Ratio-Tests (Clauset, Shalizi & Newman 2009).
- **Warum:** Clustergrößen sind ganze Zahlen. Der Bericht benutzte den
  kontinuierlichen Schätzer, der liegt ~0,02 zu hoch. Und ein Fit allein zeigt
  nicht, ob überhaupt ein Potenzgesetz vorliegt, das leisten erst die Tests.
- **Alternative:** Gerade durch ein Log-log-Histogramm legen. Das ist als
  unzuverlässig bekannt.

### E02 Energie-Sampling ist die Ursache für S ≈ 1,4 · gilt
- **Befund:** Der größte Cluster wächst wie E^2,9. Ist E log-uniform verteilt,
  folgt daraus p(n) ~ n^−1, ganz ohne Kaskadenphysik. In einzelnen
  Energiebändern kommt 1,4 nicht vor.
- **Entscheidung:** Ab jetzt nur noch Ensembles bei **fester** Energie.
- **Folge:** Das Ergebnis im Bericht und der Vergleich mit Sand et al. sind in
  dieser Form nicht haltbar.

### E03 Design des Festenergie-Ensembles · gilt für Teil 2
- **Was:** E = 1300 / 1800 / 2400, je 300 Läufe, 250², dt = 3·10⁻⁴,
  12 000 Schritte, Healing an.
- **Gepaart:** Einschlagort (x, y aus 0,4–0,6), Winkel (0–360°) und Seed hängen
  nur von der Lauf-Nummer ab. Lauf 42 hat bei jeder Energie dieselbe
  Geometrie, so bleibt beim Energievergleich nur die Energie als Unterschied.
- **Warum diese drei Energien:** keine besondere Begründung, relativ zufällig
  gewählt (Paul, 30.09.).

### E04 Feine Clusterdefinition, sechs Definitionen als Sensitivitätstest · gilt
- **Was:** Ein Atom zählt als verlagert bei `disp > 1,0 L0`, zwei Atome gehören
  zum selben Cluster bei Abstand < 1,0 L0. Die alte Definition war
  0,5 / 1,5.
- **Warum:** Mit der alten Definition verschmilzt der Schaden zu einem Klumpen
  (92 % der Atome im größten Cluster). Die feine Definition ergibt getrennte
  Cluster. So zählen auch Sand et al.
- **Folge:** S hängt stark von der Definition ab (1,23–3,30). Jede S-Angabe
  braucht deshalb die Definition dazu.

### E05 Eigene Nachrechnung von α nahe 3 · gilt
- **Warum:** Das Paket `powerlaw` klemmt α stillschweigend bei 3,000.
  `mle_unbounded()` rechnet solche Werte frei nach.

---

## Teil 3: ZBL-Umbau (seit 24.09.)

### E06 Wirkungsquerschnitt über ein ZBL-Potenzial („Variante 2“) · gilt
- **Was:** Die Abstoßung auf kurze Distanz ist das abgeschirmte
  Coulomb-Potenzial nach Ziegler, Biersack und Littmark. Schnelle Atome kommen
  damit näher heran, der Querschnitt schrumpft von selbst mit der Energie.
- **Warum:** So arbeitet die Kaskaden-MD in der Forschung, auch Sand et al. Die
  Energie bleibt erhalten, und alle Physik-Tests bleiben möglich.
- **Verworfen:**
  - Binary-Collision-Approximation (wie SRIM): kennt kein Healing, keine
    Relaxation und keine Clusterbildung, also gerade das nicht, was wir messen.
  - Geschwindigkeitsabhängiges RCUT: künstlich, verletzt die Energieerhaltung.
  - Nur ein kleineres REP_N: richtige Richtung, aber nicht kalibrierbar.

### E07 Material Wolfram · gilt
- **Warum:** Vergleich mit Sand et al. 2013 (3D-W). Deine Wahl.
- **Folge:** L0 = 2,74 Å, m = 183,84 u, Z = 74.

### E08 Energieanker E_d = 90 eV · gilt, Umsetzung offen (E17)
- **Was:** Die Verlagerungsschwelle E_d des Modells soll dem Wolfram-Wert
  entsprechen (ASTM-Mittel 90 eV). Daraus folgt ε, also eV pro Modellenergie.
- **Warum:** Übliche Kalibriergröße für Kaskadensimulationen. Deine Wahl.

### E09 Befund: Gebundene Nachbarn sind durchlässig · Befund
- Ab E ≈ 200 fliegt ein Atom durch seinen gebundenen Nachbarn. Folge von S3.
  → [`befund_gebundene_nachbarn.md`](../paper/docs/befund_gebundene_nachbarn.md)

### E10 ZBL für alle Paare, r⁻¹²-Wand im ZBL-Modus aus · gilt
- **Warum:** Behebt E09. Die r⁻¹²-Wand ist außerdem viel zu hart:
  V(0,8 L0) ≈ 50 Einheiten, das wären bei ε = 10 eV rund 500 eV bei 2,2 Å.
  ZBL gibt dort ~8 eV. Bliebe sie stehen, würde sie den Querschnitt weiter
  bestimmen.
- **Getestet:** Zusätzlich eingeschaltet (`zbl_keep_r12`) hebt sie E_d nur um
  ~20 %. Sie bleibt deshalb aus.
- **Übergang:** ZBL wird zwischen 0,5 und 0,85 L0 glatt auf 0 geführt. Die
  Ruhelage bei L0 bleibt dadurch unverändert.

### E11 Altes Modell bleibt bitgleich erhalten · gilt
- **Was:** `rep_model = r12` ist Standard und liefert byte-identische Ergebnisse
  zum Stand vor dem Umbau. Das wird in der CI geprüft.
- **Warum:** Alle alten Ergebnisse bleiben reproduzierbar. Ein r12-Kontrollensemble
  trennt später den Effekt des Potenzials von allem anderen.

### E12 Nur Laptop, nur serielle Version · gilt
- **Was:** 8 Kerne, 30 GB, abends und nachts. Kein MPI, kein CUDA.
- **Warum:** CUDA vorerst weglassen war deine Wahl. MPI habe ich zusätzlich
  weggelassen: Auf 8 Kernen sind 8 parallele serielle Läufe effizienter als ein
  MPI-Lauf, weil die Kommunikation zwischen den Prozessen entfällt.
- **Folge:** MPI und CUDA rechnen weiter das alte Modell.

### E13 Größere Gitter nur zusammen mit dem neuen Potenzial · gilt
- **Befund:** Bei 250² erreicht keine Kaskade den Rand (Ausdehnung
  17–35 L0). Im alten Modell begrenzt die Physik, nicht die Box.
- **Entscheidung:** Das Gitter wird erst größer, wenn die Atome weit fliegen
  können. Ziel für die Produktion: 700², etwa 13 h über Nacht. Die endgültige
  Größe kommt aus den Pilotläufen.

### E14 Adaptiver Zeitschritt, Δx_max = 0,0025 L0 · gilt
- **Warum:** ZBL ist auf kurze Distanz sehr steif. Ein fester Zeitschritt wäre
  entweder zu grob für den harten Stoß oder zu fein für den Rest.
- **Wert:** Beim zentralen Stoß mit E = 10⁵ ergibt 0,005 einen Energiefehler von
  7,6·10⁻⁴, 0,0025 nur 1,1·10⁻⁴.

### E15 Elektronische Bremsung nach Lindhard-Scharff mit 3D-Dichte · gilt
- **Warum:** Schnelle Atome verlieren real einen großen Teil ihrer Energie an
  die Elektronen. Mit der 3D-Dichte verliert das Modellatom pro Strecke so viel
  wie in echtem Wolfram. Unter 10 eV wird nicht gebremst, das ist die übliche
  Konvention.

### E16 Absorbierender Rand ist Pflicht · gilt
- **Warum:** Ohne Wärmesenke heizt die Kaskade die ganze Box auf. Bei 150 keV
  und 700² wären es 0,77 Einheiten pro Atom, fast eine Bindungsenergie. Die
  abgeführte Energie wird mitgezählt.

### E17 Bindungsstärke und Energieskala ε · **offen**
- **Problem:** Die Kalibrierung über E_d ist mehrdeutig. P(E) ist
  zweigipflig, und je nach Definition von E_d kommt ε = 0,21 / 0,40 / 0,6 eV
  pro Einheit heraus. Die Schallgeschwindigkeit ergibt unabhängig davon
  ε ≈ 0,44. In jedem Fall sind die Bindungen 5–13× zu schwach
  (Kohäsionsenergie 0,7–2,0 eV, Wolfram 8,9 eV), und das Gitter schmilzt
  großflächig.
- **Optionen:**
  - **A (empfohlen):** MAX_STRETCH auf ≈ 1,35 anheben (ändert S2), ε ≈ 0,44
    aus der Schallgeschwindigkeit, E_d als Kontrolle nachmessen.
  - **B:** So lassen und ε aus E_d wählen. Die Schmelzzone dominiert dann.
  - **C:** Ohne eV-Skala, alles in E/E_d.
- **Details:** [`haltepunkt_2026-09-24.md`](../paper/docs/haltepunkt_2026-09-24.md)

### E18 Kalibrierung per P(E)-Scan statt Iteration · gilt
- **Warum:** Der erste Versuch passte ε iterativ an die jeweils erste
  Defektenergie an. Das pendelte (0,37 → 0,41 → 0,40 → 0,33), weil die erste
  Schwelle chaotisch auf kleine Änderungen reagiert. Der Scan über 31 Richtungen
  und 41 Energien misst stattdessen die ganze Kurve.

---

## Teil 4: Neuordnung (30.09.)

### E19 Bachelorarbeit und Paper getrennt · gilt
- **Was:** Die Bachelorarbeit untersucht, ob ein agentenbetriebenes HPC-Labor
  das Minimalmodell eigenständig in Richtung MD-Referenzen weiterentwickeln kann
  (Arbeitstitel). Das Paper rechnet echte Läufe auf dem PC2: Es prüft die
  Ergebnisse des Agenten und setzt Pauls eigene Ideen um.
- **Folge:** Der ZBL-Umbau gehört zum Paper, er ist interaktiv mit KI-Hilfe
  entstanden und nicht autonom.

### E20 Der Agent startet vom unveränderten Kursprojekt-Modell, in einem eigenen Repo · gilt
- **Was:** Startpunkt ist der Code im Stand von Tag `kursprojekt-code`, mit
  allen bekannten Fehlern (S3, S4 usw.). Der Agent sieht dieses Repo nicht.
- **Warum:** Kennt der Agent unsere Befunde, misst das Experiment nur, ob er
  abschreiben kann. Die bekannten Fehler dienen als Prüfliste für die
  Auswertung ([`zusammenfassung_fuer_ba_planung.md`](../bachelorarbeit/zusammenfassung_fuer_ba_planung.md), Abschnitt 4).
- **Offen:** Inhalt des Agenten-Repos, Referenzen und Messgrößen. Wird angelegt,
  wenn der Plan der Arbeit steht.

### E21 Kursprojekt eingefroren · gilt
- **Was:** Alles aus Teil 1 und 2 inklusive Präsentation liegt in
  `kursprojekt/` und wird nicht mehr geändert. Tag `kursprojekt-final`.
  Ausnahme: Pfadkorrekturen, damit es weiter baut.
