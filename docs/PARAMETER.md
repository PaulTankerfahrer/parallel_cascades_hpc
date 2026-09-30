# Modellparameter: Wert, Bedeutung, Herkunft

*Stand: 2026-09-30. Zurück zum [Überblick](UEBERBLICK.md). Die Begründungen
stehen in [`ENTSCHEIDUNGEN.md`](ENTSCHEIDUNGEN.md), die Kürzel (S2, E14, …)
verweisen dorthin.*

**Herkunft**
- **Physik:** Messwert von Wolfram oder Naturkonstante
- **Modell:** Modellannahme aus dem Semesterprojekt
- **Wahl:** bewusst gewählt, mit Begründung
- **Kalibriert:** aus einer Messung im Modell bestimmt
- **Numerik:** nur für Genauigkeit oder Laufzeit, ohne physikalische Bedeutung

**Wo die Werte herkommen:** Die Standardwerte stehen im Code
(`paper/src/cascade_serial.c` bzw. `kursprojekt/src/cascade_serial.c`, oben).
`src/params.ini` ist die Benchmark-Konfiguration aus dem Semesterprojekt
(1414², kein Healing). Die Physik-Läufe schreiben ihre eigene `.ini` pro Lauf:
`kursprojekt/scripts/python/run_one.py` (altes Ensemble),
`…/statistik/run_fixed_energy.py` (Teil 2) und
`paper/scripts/python/zbl/kalibrierung_ed.py`.

## Einheiten

Das Modell rechnet dimensionslos: L0 = 1, m = 1, K = 100. Eine Energieeinheit
ist m·L0²/t², eine Zeiteinheit die Zeit, in der ein Atom mit v = 1 einen
Gitterabstand zurücklegt. Die Schallgeschwindigkeit ist damit ≈ 10 L0 pro
Zeiteinheit.

Umrechnung auf Wolfram: 1 L0 = 2,74 Å, m = 183,84 u, 1 Energieeinheit =
**ε eV (noch offen, E17)**. Die Zeiteinheit folgt daraus: t = L0·√(m/ε).

> Der Kopfkommentar in `src/params.ini` zu den Einheiten ist verwirrend. Er
> misst Energien in K·L0² (1 davon = 100 Modelleinheiten) und nennt „eine
> Zeiteinheit ~0,1 s“, was physikalisch nichts bedeutet. Das stammt aus dem
> Semesterprojekt und ist nicht korrigiert. Alle Energien im Code, in den
> Skripten und in dieser Doku sind Modelleinheiten wie oben.

## Gitter und Bindungen

| Parameter | Wert | Bedeutung | Herkunft | Status |
|---|---|---|---|---|
| `L0` | 1 | Gleichgewichtsabstand der Federn = Gitterkonstante | Modell | fest |
| `MASS` | 1 | Masse jedes Atoms | Modell | fest |
| `lattice` | Dreieck | 6 Nachbarn, schubstabil (S1) | Modell | fest |
| `NX`, `NY` | 250 (Teil 2), Ziel 700 (Teil 3) | Gittergröße in Atomen | Wahl (E13) | 700 vorläufig, kommt aus Pilotläufen |
| `K_SPRING` | 100 | Federkonstante. Bestimmt die Schallgeschwindigkeit | Modell | fest |
| `MAX_STRETCH` | 1,15 | Feder reißt bei r > 1,15·L0 | Modell (S2) | **zur Diskussion (E17)** |
| `healing` | true (ab Teil 2) | Gerissene Bindungen dürfen sich wieder schließen | Wahl (S5) | fest |
| `healing_dist` | 1,10 | … wenn der Abstand kleiner als 1,10·L0 ist | Modell, frei gesetzt | fest |
| `healing_vrel` | 0,5 | … und die Relativgeschwindigkeit kleiner als 0,5 | Modell, frei gesetzt | fest |

**Abgeleitet:** Bindungsbruch kostet ½·K·(0,15)² = 1,125 Einheiten. Die
Kohäsionsenergie (3 Bindungen pro Atom) ist 3,4 Einheiten. Wolfram hat 8,9 eV.
Schall: c_L = 7,7·√ε km/s, c_T = 4,4·√ε km/s. Für Wolfram passt das bei ε ≈ 0,44.

## Abstoßung, altes Modell (`rep_model = r12`)

Wirkt nur zwischen ungebundenen Paaren (S3, E09).

| Parameter | Wert | Bedeutung | Herkunft | Status |
|---|---|---|---|---|
| `rep_model` | r12 | Standard, bitgleich zum alten Stand (E11) | Wahl | fest |
| `RCUT` | 0,9 | Reichweite der Abstoßung, zugleich Zellgröße der Nachbarsuche | Modell | fest |
| `K_REP` | 400 | Stärke der Abstoßung | Modell | fest |
| `REP_N` | 12 | Steilheit | Modell | fest |

## Abstoßung, neues Modell (`rep_model = zbl`)

Wirkt zwischen allen Paaren (E10).

| Parameter | Wert | Bedeutung | Herkunft | Status |
|---|---|---|---|---|
| `zbl_Z` | 74 | Kernladungszahl | Physik (W) | fest |
| `zbl_L0_A` | 2,74 Å | L0 in Ångström | Physik (W) | fest |
| `eps_eV` | 10 (Platzhalter) | eV pro Energieeinheit. Steckt in der ZBL-Stärke und in der Bremsung | Kalibriert (E08) | **offen (E17)**, Kandidat 0,44 |
| `zbl_rs1`, `zbl_rs2` | 0,50 / 0,85 | ZBL wird zwischen diesen Abständen glatt auf 0 geführt. Die Ruhelage bei L0 bleibt unberührt, `rs2 ≤ RCUT` hält die Nachbarsuche gültig | Numerik | fest |
| `zbl_keep_r12` | false | r⁻¹²-Wand zusätzlich, nur zum Vergleich | Wahl (E10) | fest |

## Zeitintegration

| Parameter | Wert | Bedeutung | Herkunft | Status |
|---|---|---|---|---|
| `dt` | 3·10⁻⁴ (Teil 2) | fester Zeitschritt, bei `dt_adapt` die Obergrenze | Numerik | **anheben auf ~2·10⁻³ geplant** (Phase 4) |
| `dt_adapt` | true (Teil 3) | Zeitschritt so, dass kein Atom pro Schritt weiter als `dx_max` fliegt | Numerik (E14) | fest |
| `dx_max` | 0,0025 | größte Strecke pro Schritt | Numerik (E14) | fest |
| `n_steps` | 12 000 (Teil 2) | Schrittzahl. Mit `dt_adapt` nur noch Notbremse | Numerik | – |
| `t_max` | 3 (Kalibrierung) | Simulationszeit | Numerik | Produktion offen (Pilotläufe) |

## Elektronische Bremsung

| Parameter | Wert | Bedeutung | Herkunft | Status |
|---|---|---|---|---|
| `el_stopping` | lindhard | Reibung aus der Lindhard-Scharff-Theorie (E15) | Wahl | fest |
| `mass_u` | 183,84 | Atommasse in u | Physik (W) | fest |
| `dens_A3` | 0,0632 Å⁻³ | 3D-Atomdichte | Physik (W) | fest |
| `el_cutoff_eV` | 10 | darunter keine Bremsung | Wahl, Konvention | fest |
| `damping`, `v_thresh` | berechnet | Bei `lindhard` aus den Werten oben bestimmt. Im alten Modell 0 bzw. aus | – | – |

## Rand und Abbruch

| Parameter | Wert | Bedeutung | Herkunft | Status |
|---|---|---|---|---|
| `absorb_border` | – | Breite der Randschicht mit Reibung (E16) | Numerik | Wert aus Pilotläufen |
| `absorb_gamma` | 10 | Reibung am äußeren Rand, nach innen linear auf 0 | Numerik, nicht optimiert | vorläufig |
| `stop_ekin` | – | Abbruch, wenn kein Atom mehr so viel kinetische Energie hat | Numerik | **noch nicht gegen Volllauf geprüft** |
| `stop_tmin` | – | frühester Abbruchzeitpunkt | Numerik | wie oben |

## Einschlag (PKA)

| Parameter | Wert | Bedeutung | Herkunft | Status |
|---|---|---|---|---|
| `pka_energy` | alt: log-uniform 200–2500; Teil 2: fest 1300 / 1800 / 2400 | Startenergie | Wahl (S4 → E02, E03) | Teil 3: Energien aus Pilotläufen |
| `pka_x`, `pka_y` | gleichverteilt 0,4–0,6 | Einschlagort relativ zur Box | Wahl (E03) | fest |
| `pka_angle` | gleichverteilt 0–360° | Flugrichtung | Wahl (E03) | fest |
| `pka_mass` | 1 | PKA gleich schwer wie Gitteratome | Modell | fest |
| `seed` | = Lauf-Nummer | gepaartes Design (E03) | Wahl | fest |

## Auswertung

| Größe | Wert | Bedeutung | Herkunft | Status |
|---|---|---|---|---|
| Verlagerungsschwelle `disp` | 1,0 L0 (Haupt), 0,5–2,0 als Test | ab wann ein Atom als verlagert zählt | Wahl (E04) | fest |
| Linkradius | 1,0 L0 (Haupt), 1,5 als Test | Abstand, bis zu dem zwei Atome zum selben Cluster gehören | Wahl (E04) | fest |
| E_d-Anker | 90 eV | Verlagerungsschwelle Wolfram, ASTM | Physik (E08) | fest |
| Defekt (Kalibrierung) | ≥ 1 Leerplatz nach Wigner-Seitz, gezählt nur ≥ 3 Zeilen/Spalten vom Rand | Kriterium im P(E)-Scan | Wahl (E18) | fest |
| Potenzgesetz-Fit | diskreter MLE, Bootstrap, GoF, LR | Clauset 2009 | Wahl (E01, E05) | fest |
