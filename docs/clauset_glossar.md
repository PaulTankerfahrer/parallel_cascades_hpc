# Glossar: Die Begriffe hinter dem Clauset-Refit

Alles, was in `results/statistik/fit_gepoolt/clauset_refit.md` an Fachbegriffen
vorkommt -- in der Reihenfolge, in der man es beim Lesen des Ergebnisses braucht.

---

## 1. Was überhaupt gefittet wird

### Potenzgesetz (power law)
Eine Verteilung, bei der die Wahrscheinlichkeit, einen Cluster der Größe `n`
zu finden, wie

    p(n) ~ n^(-S)

abfällt. `S` heißt **Exponent** (in der Strahlenschadensliteratur `S`, in der
Statistikliteratur meist `alpha` -- gleiche Größe, das `powerlaw`-Paket nennt
sie `alpha`).

Der Witz daran: Die Verteilung hat **keine ausgezeichnete Größenskala**. Wenn du
von Clustergröße 10 auf 100 gehst, ändert sich die Häufigkeit um denselben
Faktor wie von 100 auf 1000. Man nennt das **skalenfrei**. Das ist der
physikalisch interessante Teil: Eine Kaskade produziert nicht "typisch große"
Defektcluster, sondern Cluster aller Größen nach einer festen Regel. Bei einer
Normalverteilung wäre das anders -- da gibt es einen typischen Wert und
Abweichungen davon sind exponentiell unwahrscheinlich.

### PDF und CCDF
- **PDF** (probability density function, Wahrscheinlichkeitsdichte): "Wie viele
  Cluster haben *genau* Größe n?" Im Histogramm gezählt. Problem: bei großen `n`
  sind die Bins fast leer, das Rauschen dominiert.
- **CCDF** (complementary cumulative distribution function, komplementäre
  Verteilungsfunktion): "Wie viele Cluster sind *mindestens* so groß wie n?"
  Also `P(N >= n)`. Weil aufsummiert wird, ist die CCDF glatt und rauscharm --
  deshalb plottet man Potenzgesetze fast immer als CCDF.
  Wenn die PDF mit Exponent `S` abfällt, fällt die CCDF mit `S - 1`.

### Log-log-Plot
Beide Achsen logarithmisch. Ein Potenzgesetz wird dort zur **Geraden**, deren
Steigung der Exponent ist. Das ist der Grund, warum man Potenzgesetze so
plottet -- und gleichzeitig die Quelle des klassischen Fehlers (siehe 2).

---

## 2. Warum nicht einfach eine Gerade durch den log-log-Plot legen?

Genau das ist der Punkt von **Clauset, Shalizi & Newman (2009)**, dem
meistzitierten Methodenpapier zu diesem Thema. Ein linearer Fit
(kleinste Quadrate) im log-log-Histogramm ist statistisch **falsch**:

1. **Binning-Bias.** Das Ergebnis hängt davon ab, wie breit man die Histogramm-
   Balken wählt. Andere Balkenbreite, anderer Exponent.
2. **Falsche Fehlerstruktur.** Die kleinste-Quadrate-Methode setzt voraus, dass
   die Messfehler in allen Punkten gleich groß und normalverteilt sind. Im
   log-log-Histogramm ist beides falsch -- die Punkte im Schwanz beruhen auf
   1--2 Zählungen, die Punkte vorne auf Hunderten.
3. **Kein Test.** Eine Gerade lässt sich durch fast jede krumme Punktwolke
   legen. Dass ein Fit existiert, heißt nicht, dass die Verteilung ein
   Potenzgesetz *ist*.

Clauset et al. ersetzen das durch drei Schritte: MLE für den Exponenten,
KS-Minimierung für `x_min`, und einen echten Hypothesentest gegen Alternativen.

---

## 3. Die Schätzung des Exponenten

### Maximum-Likelihood-Schätzung (MLE)
"Likelihood" = die Wahrscheinlichkeit, **genau die beobachteten Daten** zu
sehen, wenn ein bestimmtes Modell mit bestimmten Parametern zutrifft. Die
MLE-Antwort auf "welches `S` ist das richtige?" lautet: dasjenige `S`, das
diese Wahrscheinlichkeit maximal macht. Kein Binning, keine Annahme über
Fehlerbalken, und statistisch beweisbar der beste Schätzer für große
Stichproben (er ist *konsistent* und *effizient*).

Praktisch maximiert man die **Log-Likelihood** (Logarithmus davon), weil sich
das Produkt vieler kleiner Wahrscheinlichkeiten in eine Summe verwandelt --
numerisch stabiler und leichter abzuleiten.

### Diskret vs. kontinuierlich
Clustergrößen sind **ganze Zahlen** (3 Atome, 4 Atome -- nichts dazwischen).
Es gibt deshalb zwei Varianten des Schätzers:

- **Kontinuierlich** (Näherung): geschlossene Formel
  `S = 1 + n / sum(ln(x_i / x_min))`. Schnell, aber für ganzzahlige Daten
  leicht verzerrt. Eine gebräuchliche Korrektur ersetzt `x_min` durch
  `x_min - 1/2` -- so steht es auch im Projektbericht (Gleichung 17).
- **Diskret** (exakt): erfordert die **Hurwitz-Zeta-Funktion**
  `zeta(S, x_min) = sum_{k>=0} (k + x_min)^(-S)` als Normierung und muss
  numerisch optimiert werden. Etwas teurer, aber unverzerrt.

Der Refit verwendet den **exakten diskreten** Schätzer. Das ist die
Hauptursache dafür, dass die neuen Zahlen etwas unter den Berichtswerten
liegen (siehe Ergebnisdokument).

### Standardfehler sigma
Die Unsicherheit des MLE-Exponenten, aus der Krümmung der Likelihood-Funktion
hergeleitet (Clauset Gl. 3.2):

    sigma = (S - 1) / sqrt(n_tail)

`n_tail` ist die Zahl der Datenpunkte oberhalb `x_min` -- **nicht** die Zahl
aller Cluster. Wer `x_min` hoch setzt, wirft Daten weg und bekommt einen
größeren Fehlerbalken. Das ist der eigentliche Preis der `x_min`-Wahl.

---

## 4. Der untere Cutoff x_min

### Warum es ihn gibt
Reale Verteilungen folgen dem Potenzgesetz meist nur **im Schwanz**, oberhalb
einer Schwelle `x_min`. Darunter gilt andere Physik. Bei den Kaskaden: ein
Cluster aus 1 oder 2 Atomen ist kein "Subkaskadenfragment", sondern ein
einzelnes Frenkel-Paar -- eine andere Sorte Objekt.

### KS-Minimierung: wie x_min bestimmt wird
- **Kolmogorov-Smirnov-Abstand (KS-Abstand, `D`):** der größte senkrechte
  Abstand zwischen der empirischen CCDF (den Daten) und der theoretischen CCDF
  (dem Fit). Ein einziger Zahlenwert für "wie weit ist der Fit von den Daten
  weg". `D = 0` wäre perfekt.
- **Das Verfahren:** Für *jeden* Kandidaten `x_min` wird ein eigener MLE-Fit
  gemacht und `D` berechnet. Gewählt wird das `x_min` mit dem **kleinsten `D`**.
  Idee: dort beginnt der Bereich, in dem das Potenzgesetz am besten passt.

### Die "x_min-Falle"
Weil mit steigendem `x_min` immer weniger Punkte übrig bleiben, kann `D`
künstlich klein werden -- drei Punkte lassen sich durch jede Gerade legen. Die
KS-Minimierung "flieht dann in den Schwanz" und liefert absurde Exponenten
(im ersten Auswertungsversuch des Projekts: `S ~ 29.6` bei `x_min ~ 2571`).
Der Bericht hat das mit einer Untergrenze (mindestens 5 % der Punkte müssen
über `x_min` bleiben) abgefangen und `x_min = 4` gesetzt.

Der Refit macht den Scan **uneingeschränkt**, prüft aber zusätzlich, was `D`
über den ganzen `x_min`-Bereich tut (Plot `clauset_xmin_scan.png`). Das ist
die ehrlichere Variante: man sieht, ob es überhaupt einen Bereich gibt, in dem
das Potenzgesetz besser passt.

---

## 5. Die Unsicherheit: Bootstrap

### Bootstrap
Verfahren, um die Streuung einer Schätzung zu bestimmen, ohne die wahre
Verteilung zu kennen. Prinzip: Man tut so, als *wäre* die Stichprobe die
Grundgesamtheit. Aus den 2391 gemessenen Clustergrößen zieht man 1000-mal je
2391 Werte **mit Zurücklegen** (manche Werte kommen mehrfach vor, andere gar
nicht), fittet jedes Mal neu und schaut, wie stark `S` dabei schwankt.

Wichtig: In jedem Resample wird **auch `x_min` neu geschätzt**. Sonst
unterschätzt man die Unsicherheit, weil die `x_min`-Wahl selbst eine
Fehlerquelle ist.

### Perzentil-Konfidenzintervall
Aus den 1000 Bootstrap-Werten von `S` nimmt man das 2,5-%- und das
97,5-%-Perzentil. Dazwischen liegen 95 % der Werte -- das ist das
**95-%-Konfidenzintervall**. Es ist verlässlicher als `S ± 2 sigma`, weil es
keine Symmetrie oder Normalverteilung voraussetzt.

---

## 6. Der Härtetest: Passt das Potenzgesetz überhaupt?

### Goodness-of-Fit-p-Wert (semiparametrischer Bootstrap)
Die entscheidende Frage: Ist der gemessene KS-Abstand `D` klein genug, dass die
Daten *plausibel* aus einem Potenzgesetz stammen könnten? Clauset et al.
beantworten das so (Abschnitt 4.1):

1. Fitte das Potenzgesetz an die echten Daten -> `S`, `x_min`, `D_beobachtet`.
2. Erzeuge 1000 **synthetische** Datensätze derselben Größe, die *garantiert*
   aus diesem Potenzgesetz stammen: oberhalb `x_min` echte Zufallszahlen aus
   dem gefitteten Potenzgesetz, unterhalb `x_min` gezogen aus den echten Daten
   (deshalb "semiparametrisch" -- halb Modell, halb Empirie).
3. Fitte jeden synthetischen Datensatz mit demselben Verfahren -> 1000 Werte
   `D_synthetisch`.
4. **p = Anteil der synthetischen Datensätze mit `D_synthetisch >= D_beobachtet`.**

Interpretation:
- **p > 0,1**: Die echten Daten weichen nicht stärker vom Potenzgesetz ab als
  echte Potenzgesetz-Daten es zufällig tun. Das Potenzgesetz ist **plausibel**.
- **p <= 0,1**: Die Abweichung ist zu groß. Das Potenzgesetz wird **verworfen**.

Zwei Fallstricke:
- Ein großes `p` **beweist** kein Potenzgesetz. Es heißt nur "nicht widerlegt".
  Andere Verteilungen können genauso gut passen -- dafür gibt es die
  LR-Tests (7).
- Bei sehr großen Stichproben wird *jedes* idealisierte Modell irgendwann
  verworfen, weil selbst winzige systematische Abweichungen statistisch
  signifikant werden.

---

## 7. Der Vergleich: Potenzgesetz gegen Alternativen

Ein Modell allein zu testen reicht nicht. Die relevante Frage ist:
**Beschreibt eine andere schwerschwänzige Verteilung die Daten besser?**

### Likelihood-Ratio-Test (LR-Test)
Man berechnet für jeden Datenpunkt die Log-Likelihood unter Modell A
(Potenzgesetz) und unter Modell B (Alternative) und bildet die Differenz-Summe

    R = sum( ln L_A(x_i) - ln L_B(x_i) )

- `R > 0`: Das Potenzgesetz erklärt die Daten besser.
- `R < 0`: Die Alternative erklärt sie besser.
- `R ~ 0`: Unentscheidbar.

Beide Modelle werden dabei **oberhalb desselben `x_min`** gefittet -- sonst
vergleicht man Äpfel mit Birnen.

### Vuong-Test und das normierte R
`R` allein sagt nichts, solange man nicht weiß, ob es größer ist als die
zufällige Schwankung. Vuong (1989) zeigt: Bei nicht-verwandten Modellen ist
`R` näherungsweise normalverteilt, wenn man es auf seine Streuung normiert:

    R_normiert = R / (sqrt(n) * sd(Einzeldifferenzen))

Daraus folgt ein **p-Wert**. `p > 0,1` bedeutet: Das Vorzeichen von `R` ist
statistisch nicht abgesichert -- die beiden Modelle sind auf diesen Daten
**nicht unterscheidbar**. Das ist ein sehr häufiges und völlig legitimes
Ergebnis, gerade bei Lognormal vs. Potenzgesetz.

### Genestete Modelle
Zwei Modelle sind **genestet** (nested), wenn eines ein Spezialfall des anderen
ist. Das Potenzgesetz ist der Grenzfall `lambda -> 0` des Potenzgesetzes mit
Cutoff. Hier ist der Vuong-Test nicht anwendbar; stattdessen ist `2R`
chi-quadrat-verteilt mit einem Freiheitsgrad (Wilks-Theorem). Das Skript
behandelt diesen Fall separat -- daran ist `[genestet, chi^2]` in der Ausgabe
zu erkennen. Wichtig: Ein größeres Modell passt **immer** mindestens so gut;
die Frage ist nur, ob die Verbesserung den zusätzlichen Parameter rechtfertigt.

### Die vier Alternativen
| Modell | Form | Physikalische Bedeutung |
|---|---|---|
| **Lognormal** | Logarithmus der Größe ist normalverteilt | Entsteht aus **multiplikativen** Zufallsprozessen (viele Faktoren multiplizieren sich). Sieht im log-log-Plot über begrenzte Bereiche fast wie eine Gerade aus -- der klassische Konkurrent, extrem schwer zu unterscheiden. |
| **Exponential** | `p(n) ~ exp(-lambda n)` | Die *leichtschwänzige* Referenz: eine charakteristische Größenskala `1/lambda`. Wenn das Potenzgesetz hiergegen nicht gewinnt, ist die Skalenfreiheit gar nicht belegt. Das ist die Mindesthürde. |
| **Potenzgesetz mit exponentiellem Cutoff** | `p(n) ~ n^(-S) exp(-lambda n)` | Potenzgesetz, das oberhalb `n ~ 1/lambda` abgeschnitten wird. Physikalisch genau das, was ein **endliches Gitter** erzwingt -- ein 250x250-System kann keine beliebig großen Cluster tragen. Der wichtigste Kandidat für dieses Projekt. |
| **Stretched Exponential** (Weibull) | `p(n) ~ n^(b-1) exp(-lambda n^b)` | Zwischen Potenzgesetz und Exponential; typisch für Systeme mit einer Verteilung von Relaxationszeiten. |

---

## 8. Begriffe aus der Datenaufbereitung

### Cluster / Union-Find
Ein **Defektcluster** ist eine zusammenhängende Gruppe verlagerter Atome.
Die Definition im Projekt:
- Ein Atom gilt als **verlagert**, wenn seine Verschiebung `disp > 0,5 * L0`
  ist (`L0` = Gitterkonstante).
- Zwei verlagerte Atome gehören zum selben Cluster, wenn ihr Abstand
  `< 1,5 * L0` ist (**Verknüpfungsradius**, link distance).
- **Union-Find** ist der Algorithmus, der daraus die Zusammenhangskomponenten
  berechnet: jedes Atom zeigt auf einen "Vertreter", verbundene Atome werden
  zusammengelegt. Läuft in nahezu linearer Zeit.

Beide Schwellwerte sind **Konventionen**, keine Naturkonstanten -- eine
Sensitivitätsanalyse darüber wäre ein sinnvoller Schritt für die Bachelorarbeit.

### Pooling
Eine einzelne Kaskade liefert nur wenige Cluster (hier: 2--3 im Schnitt). Für
eine Verteilung braucht man Statistik, also wirft man die Cluster **aller
1000 Läufe zusammen** ("poolen"). Die Läufe unterscheiden sich in PKA-Energie
(log-uniform 200--2500), Einschlagort und -winkel; gepoolt wird also über ein
**Rückstoßspektrum**, nicht über identische Bedingungen.

### Finite-Size-Cutoff
Das simulierte Gitter ist endlich (250x250 = 62 500 Atome). Cluster können
deshalb nicht beliebig groß werden -- die Verteilung bricht am oberen Ende ab
("Klippe" im log-log-Plot). Das ist ein **Artefakt der Simulationsgröße**,
keine Physik, und der Hauptgrund, warum das Modell "Potenzgesetz mit Cutoff"
ernst zu nehmen ist.

### Healing (Frenkel-Paar-Rekombination)
Nach der ballistischen Phase können sich ausgeschlagene Atome und Leerstellen
wieder zusammenfinden, wenn sie nah genug (`r < 1,10 * L0`) und langsam genug
(`|v_rel| < 0,5`) sind. Der Code kann das an- und abschalten -- daher die
beiden Datensätze `heal` und `no_heal`.

---

## 9. Der Referenzwert aus der Literatur

**Sand et al. (2013)**, "High-energy collision cascades in tungsten:
dislocation loops structure and clustering scaling laws", finden für
**3D-Wolfram** einen Exponenten `S = 1,63 ± 0,07` über fast drei
Größenordnungen. Das ist der Vergleichspunkt des Projekts.

Ein Unterschied zum 2D-Wert ist **erwartet**, nicht überraschend: In 2D hat
eine Kaskade weniger Raumrichtungen, in die sie fragmentieren kann, und das
harmonische Federmodell bildet weder ZBL-Streuung noch Einbettungsenergie ab.
Ein direkter Zahlenvergleich 2D-Modell gegen 3D-Experiment trägt daher nur
begrenzt.

---

## Kurzreferenz: Was steht in der Ergebnistabelle?

| Spalte | Bedeutung |
|---|---|
| `n` | Anzahl gepoolter Cluster insgesamt |
| `x_min` | Untere Grenze des gefitteten Bereichs (per KS-Minimierung) |
| `n_tail` | Datenpunkte oberhalb `x_min` -- nur die gehen in den Fit ein |
| `S` | Exponent des Potenzgesetzes (diskretes MLE) |
| `sigma_MLE` | Standardfehler `(S-1)/sqrt(n_tail)` |
| `95 %-CI` | Bootstrap-Konfidenzintervall (1000 Resamples, `x_min` je neu) |
| `KS D` | Kolmogorov-Smirnov-Abstand des Fits |
| `GoF p` | Goodness-of-Fit-p-Wert; `> 0,1` = Potenzgesetz plausibel |
| `R` | Log-Likelihood-Differenz: `> 0` Potenzgesetz besser, `< 0` Alternative besser |
| `p` (LR) | Signifikanz von `R`; `> 0,1` = nicht unterscheidbar |
