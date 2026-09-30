# Vortragsnotizen

Stichpunkte zum Sprechen, eine Sektion je Folie. Nichts davon vorlesen --
die Folie zeigt das Bild, du sagst den Satz dazu.
Richtwert: **8 bis 9 Minuten** ohne Rückfragen. Die Reservefolien im Anhang
sind nur für Nachfragen.

---

## 1 · Titel (~20 s)

- Projektarbeit bei Dr. Maurer: Kollisionskaskaden in einem 2D-Modellkristall.
- Warum das Thema: die erste Wand im Fusionsreaktor steht im Neutronenfluss.
  Ein Neutron stößt ein Gitteratom an, das schlägt weitere an -- eine Kaskade.
  Übrig bleiben Defektcluster, das Material versprödet.
- Interessant ist nicht die einzelne Kaskade, sondern die **Verteilung der
  Clustergrößen**. Dafür braucht man Hunderte Läufe.
- Volle Molekulardynamik ist dafür zu teuer -> Minimalmodell.
- Ansage fürs Publikum: "Ich zeige das Modell, zwei Videos, ein Ergebnis,
  warum es falsch war, und was jetzt rauskommt."

## 2 · Das Modell (~1 min)

- Dreiecksgitter, jedes Atom hat sechs Nachbarn. 250 × 250 = 62.500 Atome.
- **Linkes Bild, die Kraft zwischen zwei Nachbarn** -- das ist die ganze Physik:
  - bei 1,0 Ruhelage, keine Kraft;
  - gedehnt zieht die Feder zurück (linear, Hooke);
  - unter 0,9 sehr steile Abstoßung -- das ist der eigentliche Stoß;
  - über 1,15 (15 % Dehnung) **reißt** die Bindung, Kraft springt auf null und
    bleibt dort. Das ist der bleibende Schaden.
- Rechtes Bild: derselbe Kristall in Ruhe, ein Ausschnitt.
- Das PKA ist einfach ein Atom, das eine hohe Startgeschwindigkeit bekommt;
  danach wird nur noch Newton integriert.
- Kernsatz: **nur Massen und brechbare Federn** -- eine Kaskade in Sekunden
  statt Stunden. Die Frage ist, wie viel Statistik davon übrig bleibt.

## 3 · Video "Geschwindigkeitspeak" (~45 s)

Was zu sehen ist (Farbe = Betrag der Geschwindigkeit, blau = ruhig,
rot = schnell; Blick von oben, das Video läuft 15 s):

- Anfang: ein kaltes, dunkelblaues Gitter, in der Mitte eine Handvoll bunter
  Punkte -- das PKA und die ersten angestoßenen Nachbarn.
- Dann verzweigt es sich: aus wenigen schnellen Atomen wird ein wachsender,
  tiefroter Fleck. Die Energie verteilt sich auf immer mehr Atome.
- Die Zone wird rund und weitet sich; am Rand läuft eine helle grün-türkise
  Front nach außen ins kalte Gitter -- die Druckwelle.
- Zum Schluss kühlt die Mitte ab, die Farben werden dunkler und unruhig
  gesprenkelt: geordnetes Gitter wird zu ungeordneter, "flüssiger" Zone.
- Sprechsatz: das ist der **thermische Spike** -- lokal schmilzt der Kristall
  auf und friert dann wieder ein. Was beim Einfrieren nicht auf seinen Platz
  zurückfindet, ist der Schaden.

## 4 · Video "Defektcluster" (~45 s)

Was zu sehen ist (derselbe Lauf, andere Färbung):

- Grau/dunkel: Atome, die auf ihrem Gitterplatz sitzen. Farbig: Atome, die
  mehr als **eine Gitterweite** von ihrem eigenen Platz weg sind.
- **Jede Farbe ist ein eigener Cluster** -- alles, was näher als eine
  Gitterweite beieinander liegt, gehört zusammen.
- Unten im Bild läuft der Zähler mit: erst 33 Cluster, größter 3 Atome
  (nur vereinzelte Atome), dann 51 Cluster / 17 Atome, am Ende
  **72 Cluster, größter 239 Atome**.
- Man sieht also, wie aus einzelnen verschobenen Atomen ein zusammenhängendes
  Schadensgebiet mit einem großen und vielen kleinen Clustern wird.
- Sprechsatz: **genau diese farbigen Objekte werden gezählt** -- ihre
  Größenverteilung ist das Ergebnis der Arbeit.

## 5 · Projektergebnis (~45 s)

- 1000 Kaskaden gerechnet, alle Cluster in einen Topf, doppelt logarithmisch
  aufgetragen: Anteil der Cluster mit mindestens n Atomen.
- Die Punkte liegen auf einer Geraden -> **Potenzgesetz**, S ≈ 1,4.
- Bedeutet: es gibt keine typische Clustergröße; kleine häufig, große selten,
  im selben Verhältnis über alle Größen.
- Mit und ohne Defektheilung praktisch derselbe Wert -- sah robust aus.
- Und: 1,4 liegt nah an dem, was in der Literatur für Wolfram steht (1,63).
  Das sah nach einem schönen Ergebnis aus.

## 6 · Projekt überprüft (~1 min)

- Statt es zu glauben, durchgetestet (Verfahren nach Clauset):
  - **Anpassungsgüte**: Passt das Potenzgesetz überhaupt? Nein.
  - **Vergleich mit anderen Verteilungen**: Lognormal und Potenzgesetz mit
    Cutoff beschreiben die Daten besser.
- Der Grund (das ist der Punkt der Folie): Die Einschlagsenergie war in jedem
  Lauf **zufällig gewürfelt**, jeder Zehnerbereich gleich oft. Die Clustergröße
  hängt fast nur von der Energie ab -- die Gerade kommt also aus meiner
  Energieauswahl, nicht aus der Physik der Kaskade.
- Deshalb ist auch der Vergleich mit Sand et al. hinfällig: dort ~100 Kaskaden
  bei **einer** festen Energie, hier 1000 bei gestreuten Energien.
- Ruhig offen sagen: das erste Ergebnis war ein Artefakt des Aufbaus.

## 7 · Clusterdefinition angepasst (~1 min)

- Konsequenz: neues Ensemble, **900 Läufe bei fester Energie**, drei Stufen
  (1300 / 1800 / 2400) mit je 300 Läufen an denselben Stellen und Winkeln.
- Zweiter Fehler, der dabei auffiel: die Clusterdefinition. Alt (verschoben ab
  0,5, zusammen unter 1,5) klebt die ganze Schadenszone zu **einem** Klumpen
  zusammen -- 92 % aller verschobenen Atome in einem einzigen Cluster.
  Damit gibt es gar keine Verteilung zu messen.
- Neu (1,0 / 1,0): dieselben Atome, aber 30 bis 70 getrennte Cluster pro Lauf.
  Das ist auch die Definition, die Sand et al. verwenden.
- Im Bild: die Verteilung geht jetzt über zwei Größenordnungen und die Gerade
  passt sichtbar besser ("Abw." = Abstand der Punkte von der Geraden).

## 8 · Neues Ergebnis (~1 min)

- Wichtigste Erkenntnis: **S ist keine feste Zahl des Modells, sondern hängt
  an der Definition.** Sechs Definitionen auf denselben 900 Läufen ergeben
  Werte zwischen 1,23 und 3,30 (Tabelle).
- Wer einen Exponenten angibt, muss die Clusterdefinition dazu angeben.
- Bei der Sand-Definition wird S mit steigender Energie kleiner:
  2,50 -> 2,19 -> 2,00, und die Gerade passt dabei immer besser.
  Ob sich das auf einen Wert einpendelt, ist offen.
- Größenvergleich: Wolfram in 3D hat 1,63, hier in 2D rund 2,0 -- also steiler,
  große Cluster sind seltener. Passt qualitativ dazu, dass jedes Atom in 2D
  weniger Nachbarn hat.

## 9 · Plan für die Bachelorarbeit (~1 min)

- Kurz durchgehen, nicht ausformulieren:
  1. **Pendelt sich S ein?** Verdacht: 250 × 250 ist zu klein, die größten
     Cluster stoßen an den Rand. Größeres Gitter, höhere Energien, Rechencluster.
  2. **Definition systematisch scannen** -- gibt es einen Bereich, in dem S
     stabil bleibt? Das wäre die sinnvolle Definition.
  3. **Parameter an Wolfram anpassen** (Bindungsenergie, Verlagerungsschwelle).
  4. **3D** -- dann sind die Literaturwerte direkt vergleichbar.
  5. **Hauptfrage:** taugt das billige Modell als Vorfilter, um viele
     Parameter grob zu durchsuchen und teure MD nur gezielt einzusetzen?
- Schlusssatz: Aus dem Projekt ist vor allem eine methodische Lehre geworden --
  ein Potenzgesetz muss man testen, nicht nur ansehen.

---

### Wenn die Zeit knapp wird

- Folie 5 und 6 zusammenziehen: "Erst kam 1,4 raus, der Test hat gezeigt,
  dass das an der Energieauswahl lag."
- Das Geschwindigkeits-Video nach ~8 s abbrechen, das Cluster-Video ganz
  laufen lassen -- daran hängt die Clusterdefinition.

### Fragen, die kommen könnten

- *Warum 2D?* Kosten. 3D ist Punkt 4 des Plans.
- *Warum reißen die Federn bei 15 %?* Freier Modellparameter, an keinem
  Material geeicht -- das ist Punkt 3 des Plans.
- *Wie viele Cluster stecken in der Statistik?* 900 Läufe, 30-70 Cluster je
  Lauf. Reservefolie mit beiden Definitionen liegt hinten.
