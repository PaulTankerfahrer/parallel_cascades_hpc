# Kleines Glossar zur Präsentation

Einfach erklärt, zum Nachschlagen vor dem Vortrag. Die technische
Langfassung steht in [`../docs/clauset_glossar.md`](../docs/clauset_glossar.md).

---

## 1 · Die linke Grafik auf „Das Modell": die Kraftkurve

**Was aufgetragen ist:** waagerecht der Abstand zweier benachbarter Atome
(1,0 = normaler Gitterabstand), senkrecht die Kraft zwischen ihnen.
Oben = die Atome stoßen sich ab, unten = sie ziehen sich an, auf der Nulllinie
= keine Kraft.

**Was man ablesen soll** — die Kurve *ist* das ganze Modell, von links nach
rechts gelesen:

| Bereich | Was passiert |
|---|---|
| unter 0,9 | Die Kurve schießt steil nach oben: kommen sich zwei Atome zu nahe, drücken sie sich sehr stark auseinander. Das ist der eigentliche **Stoß**. |
| 0,9 bis 1,0 | Zusammengedrückte Feder, drückt zurück (Kraft leicht positiv). |
| bei 1,0 | Ruhelage, Kraft null. Hier sitzt jedes Atom im ungestörten Kristall. |
| 1,0 bis 1,15 | Gedehnte Feder, zieht die Atome zusammen. Der Zug wächst gleichmäßig mit der Dehnung (Hookesches Gesetz, deshalb eine Gerade). |
| bei 1,15 | Die Feder **reißt**: 15 % Dehnung ist die Grenze. Die Kraft springt auf null — das ist der orange Kreis mit der gestrichelten Linie. |
| über 1,15 | Nichts mehr. Die Bindung bleibt gerissen, auch wenn die Atome später wieder zusammenkommen. Das ist der **bleibende Schaden**. |

**Der Satz dazu:** Mehr Physik steckt nicht drin — Massen, eine Feder, die
reißt, und eine harte Abstoßung. Alles, was in den Videos passiert, folgt aus
dieser einen Kurve.

---

## 2 · Die Zeile `Häufigkeit ∝ n^(−S)`, S ≈ 1,4

Zeichen für Zeichen:

- **n** — die Größe eines Clusters, also die Zahl der Atome darin.
- **∝** — „ist proportional zu", gesprochen „geht wie". Kein Gleichheitszeichen,
  weil der Vorfaktor egal ist: es geht nur darum, *wie* die Häufigkeit mit der
  Größe abfällt, nicht um absolute Zahlen.
- **n^(−S)** — n hoch minus S, also `1 / n^S`. Das Minus heißt: je größer n,
  desto seltener. Große Cluster sind selten, kleine häufig.
- **S** — der Exponent, die einzige interessante Zahl. Er sagt, *wie schnell*
  die Häufigkeit abfällt. Großes S = große Cluster praktisch ausgeschlossen,
  kleines S = auch sehr große Cluster kommen noch vor.

**In Worten:** „Wie oft ein Cluster mit n Atomen vorkommt, fällt mit einer festen
Potenz von n ab."

**Warum das interessant ist:** Der Faktor hängt nur vom *Verhältnis* der Größen
ab, nicht von der Größe selbst. Von 10 auf 100 Atome wird ein Cluster um
denselben Faktor seltener wie von 100 auf 1000. Es gibt also **keine typische
Clustergröße** — man nennt das *skalenfrei*. Bei einer Normalverteilung wäre es
umgekehrt: ein typischer Wert, und alles weit davon weg ist praktisch
unmöglich.

**Zahlenbeispiel für S = 1,4:** zehnmal so großer Cluster → rund 25-mal
seltener (10^1,4 ≈ 25).

---

## 3 · Wie die Grafik auf „Projektergebnis" entsteht

Vier Schritte:

1. **Zählen.** Jeder Lauf liefert seine Cluster, jeder Cluster eine Zahl n.
   1000 Läufe zusammengeworfen ergeben eine lange Liste von Clustergrößen.
2. **CCDF bilden** (siehe unten) statt eines Histogramms.
3. **Beide Achsen logarithmisch** auftragen (siehe unten).
4. **Gerade anpassen** — nicht mit dem Lineal, sondern mit dem
   Maximum-Likelihood-Verfahren nach Clauset (siehe Punkt 5). Deren Steigung
   liefert S.

### Logarithmische Skala

Auf einer normalen Achse bedeutet gleicher Abstand gleiche *Differenz*
(1, 2, 3, 4 …). Auf einer logarithmischen Achse bedeutet gleicher Abstand
gleichen *Faktor* (1, 10, 100, 1000 …). Nötig, weil hier alles über mehrere
Größenordnungen läuft: Cluster von 1 bis fast 3000 Atomen, Häufigkeiten von
100 % bis 0,04 %. Auf einer normalen Achse würde alles Interessante am linken
Rand zusammenkleben.

**Der Trick:** Ein Potenzgesetz wird auf einer log-log-Achse zu einer
**Geraden**. Eine Gerade erkennt das Auge sofort — deshalb plottet man das so.

### CCDF

Ausgeschrieben *complementary cumulative distribution function*, deutsch etwa
„Überlebensfunktion". Sie beantwortet die Frage:

> **Welcher Anteil aller Cluster hat mindestens n Atome?**

Beispiel: Bei n = 10 steht 0,3 → 30 % aller Cluster haben 10 oder mehr Atome.
Ganz links ist sie immer 1 (100 % der Cluster sind mindestens so groß wie der
kleinste), nach rechts fällt sie ab.

**Warum nicht einfach ein Histogramm?** Ein Histogramm zählt „genau n Atome".
Bei großen Clustern liegt dann in jedem Balken eine 1 oder eine 0 — reines
Rauschen, und das Ergebnis hängt auch noch von der Balkenbreite ab. Die CCDF
summiert auf, ist dadurch glatt und braucht keine Balkenwahl.

**Ein Detail, falls jemand nachrechnet:** Die CCDF fällt um genau 1 flacher ab
als die Häufigkeit selbst. Zu S = 1,4 in der Formel gehört also eine Gerade mit
Steigung −0,4 im Bild. Auf der Folie steht die Zahl aus der Formel, nicht die
gemessene Steigung.

### Was in der Grafik zu sehen ist

Punkte = Messwerte, blasse Linie = angepasstes Potenzgesetz. Blau und orange
sind derselbe Datensatz mit und ohne Defektheilung — praktisch deckungsgleich.
Der senkrechte Absturz ganz rechts ist kein Effekt der Physik, sondern die
endliche Gittergröße: größer als die Schadenszone kann kein Cluster werden.

---

## 4 · Die Grafik auf „Clusterdefinition angepasst"

**Dieselbe Art Bild wie auf „Projektergebnis"** — CCDF, beide Achsen
logarithmisch, Punkte = Messwerte, blasse Linie = angepasstes Potenzgesetz.
Zwei Unterschiede: es ist die **neue Clusterdefinition** (verschoben ab 1,0,
zusammengehörig unter 1,0), und es sind **drei getrennte Kurven** statt einer.

**Die drei Farben sind die drei Energiestufen**, je 300 Läufe:
blau E = 1300, orange E = 1800, grün E = 2400. Alle drei starten links oben im
selben Punkt (1 | 1) — 100 % der Cluster haben mindestens ein Atom; dort liegen
die Punkte übereinander.

**Die Legende:**
- **S** — der Exponent, also die Steilheit der angepassten Geraden.
- **Abw.** — der größte senkrechte Abstand zwischen den Messpunkten und der
  angepassten Geraden (in der Statistik: Kolmogorov-Smirnov-Abstand *D*).
  Ein Maß dafür, wie gut die Gerade zu den Punkten passt: **kleiner ist besser**.

**Was man sehen soll** — dreimal derselbe Trend, von blau über orange nach grün:

1. **Die Kurven reichen weiter nach rechts.** Größere Kaskade, größere Cluster:
   der größte Cluster wächst von 43 über 113 auf 330 Atome; die Statistik
   umfasst 9 281 / 14 588 / 21 580 Cluster.
2. **Die Kurven werden flacher:** S = 2,50 → 2,19 → 2,00. Große Cluster werden
   relativ häufiger.
3. **Die Punkte legen sich besser auf die Gerade:** Abw. 0,032 → 0,017 → 0,013.
   Die Verteilung nähert sich mit wachsender Kaskade einem Potenzgesetz — hat es
   aber bei keiner der drei Energien erreicht.

**Zwei Details, falls jemand nachfragt:**
- Der **Absturz am rechten Ende** jeder Kurve ist keine Physik, sondern die
  endliche Kaskadengröße: größere Cluster gibt es in diesen Läufen schlicht nicht.
- Die Gerade beginnt nicht am linken Rand, sondern erst bei der Größe, ab der
  gefittet wird (`x_min`, hier 1 bzw. 2 Atome). Unterhalb davon wird das
  Potenzgesetz gar nicht behauptet.

**Der Kontrast zur alten Definition** steckt auf der Reservefolie: dort bricht
die CCDF fast senkrecht ab, weil praktisch alle verschobenen Atome in einem
einzigen Riesencluster stecken. Hier läuft die Verteilung durchgehend über zwei
Größenordnungen — erst damit gibt es überhaupt etwas zu fitten.

---

## 5 · Die drei Vergleichsverteilungen

Auf „Projekt überprüft" wird das Potenzgesetz gegen drei andere Kurvenformen
gehalten. Alle drei beschreiben ebenfalls „viele kleine, wenige große", nur mit
anderem Abfall. Im log-log-Bild:

| Verteilung | Alltagsbeispiel | Form im log-log-Bild |
|---|---|---|
| **Potenzgesetz** | Erdbebenstärken, Stadtgrößen | **Gerade** über den ganzen Bereich |
| **Exponentialverteilung** | Zerfall radioaktiver Kerne, Wartezeiten an der Kasse | Beginnt flach und **knickt dann immer steiler nach unten ab** — stürzt viel schneller ab als eine Gerade. Große Werte sind praktisch ausgeschlossen. |
| **Lognormalverteilung** | Einkommen, Teilchengrößen beim Mahlen | Eine **sanft gebogene Kurve** statt einer Geraden: nach rechts wird sie allmählich immer steiler. Über einen kurzen Bereich sieht sie fast gerade aus — genau deshalb wird sie so oft mit einem Potenzgesetz verwechselt. |
| **Potenzgesetz mit Cutoff** | alles, was eine natürliche Obergrenze hat | Erst eine **Gerade**, die am oberen Ende **abknickt** und weggeht. Das Potenzgesetz gilt, aber nur bis zu einer Grenze. |

Kurzfassung der Unterschiede:

- **Exponential:** Es gibt eine typische Größe, alles Größere fällt schlagartig weg.
- **Lognormal:** Kein sauberes Potenzgesetz, aber ein langer Schwanz — entsteht
  typischerweise, wenn viele zufällige Faktoren sich *multiplizieren*.
- **Potenzgesetz mit Cutoff:** Potenzgesetz plus eine Obergrenze. Hier plausibel,
  weil das Gitter endlich ist — die Cluster *können* gar nicht beliebig groß werden.

Ergebnis im Projekt: Gegen die Exponentialverteilung gewinnt das Potenzgesetz,
gegen Lognormal und Cutoff verliert es.

---

## 6 · Die beiden Referenzen in zwei Sätzen

**Sand et al. (2013)** — *High-energy collision cascades in tungsten:
dislocation loops structure and clustering scaling laws.*
Molekulardynamik-Kaskaden in echtem **Wolfram in 3D**: rund 100 Kaskaden bei
**einer festen** Energie (150 keV). Ergebnis: die Größenverteilung der
Defektcluster folgt einem Potenzgesetz mit **S = 1,63 ± 0,07**. Das ist der
Wert, an dem sich dieses Modell messen lässt — und der Grund, warum die
Läufe hier ebenfalls bei fester Energie gerechnet werden müssen.

**Clauset, Shalizi & Newman (2009)** — *Power-law distributions in empirical
data.* Das Standardrezept dafür, ob Daten wirklich einem Potenzgesetz folgen,
statt nur im log-log-Plot geradeaus auszusehen. Drei Zutaten:
**(1)** Exponent per Maximum-Likelihood schätzen statt per Geradenfit,
**(2)** einen **Anpassungsgüte-Test**, der sagt, ob das Potenzgesetz überhaupt
passt, **(3)** einen **Vergleich mit Alternativen** (Lognormal, Exponential,
Cutoff) — denn „passt" heißt nichts, wenn eine andere Kurve besser passt.
Die Kernaussage des Papiers: die meisten publizierten Potenzgesetze halten
diesen Test nicht aus. Genau das ist hier auch passiert.
