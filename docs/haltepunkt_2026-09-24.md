# Haltepunkt 2026-09-24: Stand des ZBL-Umbaus

Festgehalten, bevor wir das Projekt einmal sortieren. Hier steht, wo die Arbeit
steht und was als Nächstes vorgeschlagen war. Details stehen in
[`plan_zbl.md`](plan_zbl.md).

## Wo wir stehen

Die Ausgangsfrage war: Wird die Clustergrößenverteilung zu einem Potenzgesetz,
wenn (1) der Wirkungsradius mit der Geschwindigkeit schrumpft und (2) die Gitter
viel größer werden? Dafür wird das Modell auf ein Wolfram-Analogon umgebaut.

| Phase | Inhalt | Stand |
|---|---|---|
| 0 | Referenz und Diagnose des alten Modells | erledigt. Dabei gefunden: gebundene Nachbarn sind durchlässig ([Befund](befund_gebundene_nachbarn.md)) |
| 1 | ZBL-Abstoßung (energieabhängiger Wirkungsquerschnitt), adaptiver Zeitschritt | erledigt, CI-Stufen 6–7 |
| 2 | Energieskala ε [eV pro Modelleinheit] kalibrieren | **hängt an einer Entscheidung**, siehe unten |
| 3 | Elektronische Bremsung (Lindhard-Scharff) | erledigt, CI-Stufe 8. Offen: Vergleich des Elektronenanteils mit Lindhard in den Pilotläufen |
| 4 | Absorbierender Rand, Frühabbruch, Pilotläufe | Rand und Frühabbruch eingebaut. Offen: Frühabbruch gegen vollen Lauf prüfen, dt-Obergrenze anheben, Pilotläufe |
| 5 | Produktionsensemble (Ziel 700²) plus r12-Kontrolle | nicht begonnen |
| 6 | Clauset-Auswertung, BEFUND_zbl.md | nicht begonnen |

Commits auf `clusterstatistik`, nicht gepusht: a8dd985, 55ea542, d49feca,
e41e8c6, eddcd20. Deine Präsentationsänderungen (presentation/, zwei Skripte)
sind uncommittet und unangetastet.

## Ergebnis der Kalibrierung (Phase 2)

Scan: 40×40, 31 Richtungen, 41 Energien, Defekt = Leerplatz nach t = 3.
Zahlen in [`../results/zbl/kalibrierung/ergebnis.md`](../results/zbl/kalibrierung/ergebnis.md),
Grafik `results/images/D2_p_von_e.png`.

| Variante | E_d = ∫(1−P) dE | P = 50 % bei | erste Schwelle min–max |
|---|--:|--:|--:|
| zbl ε = 0,3 | 408 | 212 | 164–717 |
| zbl ε = 0,4 | 382 | 182 | 147–571 |
| zbl ε = 0,5 | 365 | 162 | 131–510 |
| zbl + r12, ε = 0,4 | 454 | 201 | 147–717 |
| r12 (altes Modell) | 159 | 124 | 66–258 |

Drei Befunde:

1. **E_d ist nicht eindeutig.** P(E) ist zweigipflig: Etwa die Hälfte der
   Richtungen macht ab ~170 Einheiten einen Defekt, die andere erst ab ~700. Je
   nach Definition folgt aus E_d = 90 eV ε ≈ 0,21 (∫(1−P), extrapoliert), ≈ 0,40
   (erste Schwelle) oder ≈ 0,6 (Median).
2. **Die Schallgeschwindigkeit ist ein unabhängiger Anker.** c_L = 7,7·√ε km/s,
   c_T = 4,4·√ε km/s. Beide passen bei ε ≈ 0,43–0,46 zu Wolfram.
3. **Die Bindungen sind 5–13× zu schwach.** Bruch bei 15 % Dehnung
   (MAX_STRETCH = 1,15). Das ergibt eine Kohäsionsenergie von 3,4 Einheiten,
   also 0,7–2,0 eV. Wolfram hat 8,9 eV. Test (120², E = 10⁴, t = 20): 10 % aller
   Bindungen reißen, 45 % der Energie steckt in gebrochenen Bindungen, der Rand
   schluckt 31 %. Das Gitter schmilzt großflächig.

## Vorschlag, der zur Entscheidung stand

| Option | Inhalt |
|---|---|
| **A (empfohlen)** | MAX_STRETCH auf ≈ 1,33–1,39 anheben, damit die Kohäsionsenergie bei ε ≈ 0,4 zu Wolfram passt. ε aus der Schallgeschwindigkeit (≈ 0,44). E_d danach als Kontrolle nachmessen (~1 h Rechenzeit). |
| B | Schwache Bindungen behalten und ε über eine der E_d-Definitionen wählen. Die Schmelzzone dominiert, nur kleine Energien sind machbar. |
| C | Ohne eV-Skala arbeiten, alles in E/E_d. Das Schmelzproblem bleibt. |

Offene Rückfrage: A ändert eine Grundannahme aus dem Semesterprojekt. Eventuell
vorher mit Maurice klären.

## Danach (unabhängig von der Wahl)

1. ε festlegen, Kalibrierung mit neuem MAX_STRETCH wiederholen (bei A).
2. Phase 4: Frühabbruch gegen vollen Lauf prüfen, dt-Obergrenze auf ~2·10⁻³
   anheben und Energieerhaltung testen, Pilotläufe (~10 Läufe × 3–4 Energien)
   für Boxgröße und Rechenbudget.
3. Phase 5: `run_fixed_energy.py` erweitern (Potenzial, Bremsung, Rand,
   fortsetzbar), Produktion über Nacht.
4. Phase 6: Clauset-Auswertung.

Randnotiz: `presentation/figures/P2_kraft.png` (`fig_force_law.py`) zeigt die
Abstoßung auch für gebundene Nachbarn. Im alten Modell stimmt das nicht. Nicht
geändert, weil es deine Präsentation ist.
