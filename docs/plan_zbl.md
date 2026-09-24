# Plan: Energieabhängiger Wirkungsquerschnitt (ZBL) und große Gitter

*Stand: 2026-09-24. Status: Plan beschlossen, Umsetzung noch nicht begonnen.*

**Entscheidungen (2026-09-24):** Material Wolfram · Energieanker E_d = 90 eV ·
Rechnen nur auf dem Laptop (8 Kerne, 30 GB RAM) · CUDA vorerst nicht.

## Ziel

Klären, ob die Clustergrößenverteilung ein Potenzgesetz wird, wenn das Modell
das zeigt, was reale Kaskaden haben: Schnelle Atome haben einen kleinen
Wirkungsquerschnitt, fliegen lange Strecken und zerfallen in Subkaskaden.
Erst am Ende der Kaskade wird der Stoß dicht und nachbarweise.

Ausgangslage (siehe `projektstand.md`): Bei fester Energie und feiner
Clusterdefinition fallen S und KS-D mit der Kaskadengröße (S = 2,50 → 2,00),
aber der GoF-Test verwirft das Potenzgesetz. Bei 250×250 erreicht keine Kaskade
den Rand (Ausdehnung 17–35 L0). Die Box begrenzt also nicht, sondern die Physik.

## Warum das bisherige Modell das nicht leisten kann

1. **Querschnitt zu groß und zu wenig energieabhängig.** Mindestabstand beim
   zentralen Stoß mit der r⁻¹²-Abstoßung: 0,64 L0 bei E = 2400, 0,49 L0 bei
   E = 50 000. Real (W→W, 150 keV, ZBL) sind es ≈ 0,2 Å = 0,07 Nachbarabstände.
2. **Gebundene Paare sind „durchlässig".** `compute_forces()` überspringt die
   Abstoßung für gebundene Paare, dort wirkt nur die Feder. Deren Energie ist auf
   0,5·K·L0² = 50 beschränkt. Ein schnelles Atom kann daher durch einen
   gebundenen Nachbarn hindurchfliegen, statt an ihm zu streuen.
   → Wird in Phase 0 gemessen, bevor wir es ändern.
3. **Keine elektronische Bremsung.** `damping` ist vorhanden, war in allen
   Ensemble-Läufen aber 0.
4. **Keine physikalischen Einheiten.** „E = 2400" ist nicht in eV übersetzbar.
   Damit ist auch der Vergleich mit Sand et al. (150 keV) nur qualitativ.

## Material und Einheiten

Wolfram, wegen des Vergleichs mit Sand et al. (2013).

| Größe | Modell | Real (W) |
|---|---|---|
| Länge | L0 = 1 | 2,74 Å (Nachbarabstand bcc) |
| Masse | m = 1 | 183,84 u |
| Energie | 1 Einheit = ? eV | **wird in Phase 2 kalibriert** |

Anker für die Energieskala ist die **Verlagerungsschwellenenergie** E_d: Wir
messen sie im Modell und setzen sie gleich dem Wolfram-Wert (ASTM-Mittel 90 eV,
richtungsabhängig ca. 40–90 eV). Grobe Vorabschätzung: Der Bindungsbruch kostet im
Modell 0,5·K·(0,15 L0)² ≈ 1,1 Einheiten. Die Skala könnte daher so liegen, dass
E = 2400 schon im Bereich von 10–100 keV liegt. Das zu prüfen ist Teil von Phase 2.

## Phase 0: Referenz und Diagnose (erledigt 2026-09-24)

- [x] Schalter `rep_model = r12 | zbl` in `cascade_serial.c`, Standard `r12`.
      `zbl` bricht vorerst mit Fehlermeldung ab.
- [x] Abstoßung nach `src/potential.h` ausgelagert (`r12_force`, `r12_pot`).
      Vorerst nur von `cascade_serial.c` benutzt.
- [x] **Bitgleichheit geprüft:** 120×120, E = 2400, 4000 Schritte, Healing an.
      `state_final`, `energy` und `broken_bonds` sind byte-identisch zum Stand
      vor dem Umbau, mit `-O2` und mit `-O3 -march=native` (Ensemble-Flags).
- [x] **Diagnose „durchlässige Bindung": bestätigt.** Ausführlich mit
      Abbildung: [`befund_gebundene_nachbarn.md`](befund_gebundene_nachbarn.md).
      21×21-Gitter, PKA aus der Mitte, Mindestabstand zu einem der 6 gebundenen
      Nachbarn:
  - Frontal (Flugrichtung = Bindungsrichtung): Ab **E ≈ 163** erreicht das PKA
    den Ort des Nachbarn, ab **E ≈ 198** fliegt es ganz hindurch. Der Nachbar
    bewegt sich dabei kaum. Das passt zur Obergrenze der Federenergie im
    Schwerpunktsystem (E/2 gegen ½·K·L0² = 50).
  - Schräg: Der Mindestabstand ist für E = 300 bis 2400 **gleich dem
    geometrischen Stoßparameter der geraden Linie**, sin(Δθ): 0,088 bei 5°,
    0,26 bei 15°, 0,50 bei 30° Abstand zur Bindung. Die Bahn wird also gar nicht
    abgelenkt, die erste Nachbarschale ist für schnelle Atome **unsichtbar**.
  - Das gilt nicht nur für das PKA: Jedes angestoßene Atom ist anfangs an 6
    Nachbarn gebunden und fliegt oberhalb von E ≈ 160 durch diese hindurch.
    Stöße finden nur mit Atomen statt, zu denen keine Bindung besteht.
  - Die Energie bleibt erhalten, deshalb hat die CI das nie bemerkt.
  - **Bedeutung für die bisherigen Ergebnisse:** Der Befund zum Energie-Sampling
    (Teil 1) und die Abhängigkeit von der Clusterdefinition (Teil 2) bleiben
    gültig, sie hängen nicht an der Stoßphysik. Die absoluten Werte von S bei
    fester Energie gelten aber für ein Modell mit diesem Artefakt.

**Folgerung für Phase 1 (neu):** Die r⁻¹²-Wand ist für die schnelle Phase
nicht nur zu dick, sondern mit K_REP = 400 auch viel zu hart. V(0,8 L0) ≈ 50
Modelleinheiten; bei ε ~ 10 eV/Einheit wären das ~500 eV bei 2,2 Å, ZBL liefert
dort ~8 eV. Bleibt die r⁻¹²-Wand neben ZBL stehen, dominiert sie weiter den
Querschnitt ungebundener Paare, und es ändert sich wenig. Im `zbl`-Modus ersetzt
ZBL daher die r⁻¹²-Wand für **alle** Paare. Die r⁻¹²-Wand ist per Schalter
`zbl_keep_r12` zuschaltbar, als Vergleich in Phase 2.

## Phase 1: ZBL-Potenzial (erledigt 2026-09-24)

**Form.** Für alle Paare mit r < r_s2 (gebunden **und** ungebunden):

    V_ZBL(r) = A/r · φ(r/a),
    φ(x) = 0,1818 e^(-3,2x) + 0,5099 e^(-0,9423x) + 0,2802 e^(-0,4029x) + 0,02817 e^(-0,2016x)

mit a = 0,8854·a_B / (2·Z^0,23) = 0,087 Å = 0,032 L0 und
A = Z²e² in Modelleinheiten (Z = 74; hängt von der Energiekalibrierung ab).

**Übergang.** Zwischen r_s1 und r_s2 wird V glatt auf 0 geführt, mit einem
Polynom, sodass V, F und F' stetig sind (wie `pair zbl` in LAMMPS).
Vorschlag: r_s1 ≈ 0,5 L0, r_s2 ≈ 0,85 L0 < RCUT. Damit gilt:
- Ruhelage und Phononen bei r ≈ L0 bleiben unverändert (Federn wie bisher).
- `cell_size = RCUT` bleibt gültig.
- Das Healing-Fenster `RCUT ≤ r < HEALING_DIST` bleibt unberührt.
- Die r⁻¹²-Wand ist im `zbl`-Modus aus (siehe Folgerung aus Phase 0),
  optional zuschaltbar mit `zbl_keep_r12 = true`.

**Gebundene Paare.** ZBL wirkt zusätzlich zur Feder. Das behebt Punkt 2 oben.
In der Paarschleife heißt das: `is_bonded` überspringt nur noch den r⁻¹²-Anteil,
nicht den ZBL-Anteil.

**Zeitschritt.** ZBL ist bei kleinem r sehr steif, und schnelle Atome kommen nah
heran. Wir brauchen einen **adaptiven Zeitschritt** wie in Kaskaden-MD üblich:
dt = min(dt_max, Δx_max / v_max, ΔE_max / (F_max·v_max)), mit Δx_max = 0,0025 L0 (siehe Tests).
Ausgabe ist dann nach Simulationszeit, nicht nach Schrittzahl.

**Parameter (`params.ini`, Sektion `[potential]`).**

    rep_model    = zbl     ; r12 (alt, Standard) | zbl
    zbl_Z        = 74
    zbl_L0_A     = 2.74    ; Å
    eps_eV       = 10.0    ; vorläufig, wird in Phase 2 kalibriert
    zbl_rs1      = 0.50
    zbl_rs2      = 0.85
    zbl_keep_r12 = false
    ; [dynamics]
    dt_adapt     = true
    dx_max       = 0.0025
    t_max        = ...     ; Simulationszeit, n_steps nur Notbremse

**Tests (alle grün, in `scripts/ci/build_and_test.sh` als Stufen 6 und 7).**
- [x] Zwei-Körper-Streuung (`scripts/ci/zbl_scatter_test.c`): E = 10 … 10⁵,
      Stoßparameter 0 … 0,4 L0. r_min auf < 0,1 %, Laborwinkel auf < 0,1 % genau.
      Energiefehler pro Stoß: 7,6·10⁻⁴ mit Δx_max = 0,005 (schlimmster Fall: zentral,
      E = 10⁵, r_min = 0,027 L0), 1,1·10⁻⁴ mit 0,0025, 8·10⁻⁶ mit 0,001.
      → Standard 0,0025.
- [x] Durchflug-Diagnose im ZBL-Modus: frontal auf gebundenen Nachbarn, kein
      Durchflug mehr. r_min = 0,39 / 0,20 / 0,087 bei E = 100 / 1000 / 10⁴,
      Theorie ohne Feder 0,35 / 0,19 / 0,086.
- [x] Energieerhaltung 120×120, t = 0,5, ohne Dämpfung: Drift 0,10 % / 0,04 % /
      0,005 % bei E = 10³ / 10⁴ / 10⁵. Das r12-Modell hatte bei E = 2400 0,33 %.
      Schritte: 4 300 / 13 100 / 27 900 (fest bei dt = 3·10⁻⁴ wären es 1 670).
- [x] Ruhelage: 60×60 ohne PKA, t = 1: Verschiebung exakt 0.
- [x] `rep_model = r12`: byte-identisch zum Stand vor dem Umbau (-O2 und -O3 -march=native).

**Gefundener und behobener Fehler:** Die Lehrbuchform
dt = (−v + √(v² + 2a·Δx))/a löscht bei a ≈ 10⁻¹³ (Rundungsrauschen der Federn)
zu exakt 0 aus, die Simulation stand dann still. Stabil ist 2Δx/(v + √(v² + 2a·Δx)).

**Offen für Phase 2:** `zbl_keep_r12 = true` als Vergleich. Bei ε = 10 eV ist die
reine ZBL-Wand bei 0,6 L0 nur ~2 Einheiten hoch. Ob das für die dichte Endphase
reicht, zeigt sich an E_d und am Endzustand.

## Phase 2: Kalibrierung der Energieskala

- [ ] E_d im Modell messen: PKA-Energie in kleinen Schritten erhöhen, für ca. 20
      Richtungen. E_d ist die kleinste Energie, bei der ein stabiles
      Frenkel-Paar übrig bleibt (mit Healing an).
- [ ] Energieeinheit festlegen: ε = 90 eV / E_d,Modell. Daraus A für ZBL.
      **Achtung, Zirkularität:** E_d hängt selbst leicht vom ZBL-Anteil ab.
      → 1–2 Iterationen, bis ε stabil ist.
- [ ] Plausibilitätscheck: Schallgeschwindigkeit bzw. Debye-Frequenz im Modell,
      mit ε umgerechnet, gegen W vergleichen (c_W ≈ 5,2 km/s). Nicht als Anker,
      nur als Größenordnungskontrolle. Das Federmodell ist kein W-Potenzial.

## Phase 3: Elektronische Bremsung

Die vorhandene Implementierung F = −γ·v für v > v_thresh hat schon die
Lindhard-Scharff-Form (S_e ∝ v).
- [ ] γ aus dem Lindhard-Scharff-Koeffizienten für W→W bestimmen und in
      Modelleinheiten umrechnen.
- [ ] v_thresh entsprechend E_kin ≈ 10 eV setzen (übliche Konvention).
      Darunter wirkt die Bremsung nicht auf die thermische Phase.
- [ ] Test: Anteil der elektronisch verlorenen Energie über E. Mit der
      Lindhard-Partition (Damage Energy) vergleichen.

## Phase 4: Pilotläufe, dann Boxgröße festlegen

**Rechenbudget Laptop.** Referenz: Das bisherige Ensemble kostete ≈ 19 s pro Lauf
(62 500 Atome, 12 000 Schritte, 8 Läufe parallel), also ≈ 2,5·10⁻⁸ s pro
Atom und Schritt und Kern. Hochgerechnet mit dem Faktor ~3 für den adaptiven
Zeitschritt:

| Box | Atome | RAM/Lauf | Zeit/Lauf | 800 Läufe (8 parallel) |
|---|--:|--:|--:|--:|
| 250² | 62 500 | ~5 MB | ~1 min | ~2 h |
| 700² | 490 000 | ~40 MB | ~8 min | ~13 h |
| 1500² | 2,25 Mio. | ~170 MB | ~35 min | ~2,5 Tage |

**Zielgröße Produktion: 700²** (Laptop läuft abends und nachts durch, ~13 h).

Die Box ist damit **der** Kostentreiber. Hebel, in dieser Reihenfolge:
1. **Absorbierender Rand** statt großer Pufferzone. Die Box muss dann nur die
   Kaskade selbst fassen, nicht die Druckwelle.
2. Box aus den Pilotläufen so klein wie zulässig wählen, nicht pauschal groß.
3. Frühzeitig abbrechen, sobald die Kaskade abgeklungen ist (E_kin der Atome mit
   v > v_thresh ≈ 0 über mehrere 100 Schritte), statt fester Schrittzahl.
4. Erst wenn das nicht reicht: Energien reduzieren oder weniger Läufe (≥ 150 pro
   Energie).

Nicht die Box raten, sondern messen:
- [ ] Je 10 Läufe bei 3–4 Energien (z. B. 10, 30, 100, 150 keV nach Kalibrierung)
      auf einer großzügigen Box (z. B. 1000×1000, einmalig, über Nacht).
- [ ] Messen: Ausdehnung der Schadenszone, Anzahl Subkaskaden, Laufzeit,
      Simulationszeit bis zum Abklingen, Energieerhaltung.
- [ ] Box so festlegen, dass die größte Kaskade ≥ 20 % Abstand zum Rand hat.
- [ ] **Absorbierender Rand:** `absorb_border` ist ein Platzhalter. Eine
      gedämpfte Randschicht (Langevin oder linear ansteigende Reibung) einbauen,
      damit die Druckwelle nicht vom freien Rand zurückläuft. Wird bei großen
      Energien nötig.
- [ ] Frühabbruch-Kriterium einbauen und prüfen (gleiches Endergebnis wie Volllauf).
- [ ] Laufzeitbudget aus den Pilotläufen neu hochrechnen. Produktion so zuschneiden,
      dass sie in ≲ 1–2 Nächten durchläuft.

## Phase 5: Produktionsensemble

- [ ] `run_fixed_energy.py` erweitern: `--rep-model`, Parameter für die Bremsung
      und den Rand. Fortsetzbar machen: fertige Läufe überspringen, damit der
      Laptop zwischendurch zugeklappt werden kann.
- [ ] Gepaartes Design wie bisher: Einschlagort, Winkel und Seed hängen nur von der
      Lauf-Nummer ab.
- [ ] Energien: 3–4 feste Werte, je ≥ 200 Läufe (Clauset braucht ≳ 10⁴ Cluster
      pro Energie für enge Konfidenzintervalle).
- [ ] Zusätzlich ein kleines Ensemble mit `r12` auf gleicher Box und Energie, um
      den Effekt des Potenzials allein zu isolieren.

## Phase 6: Auswertung

- [ ] Clusteranalyse mit allen sechs Definitionen (`analyze_fixed_energy.py`).
- [ ] Voller Clauset-Fit (`clauset_refit.py`): Bootstrap, GoF, LR-Tests. Die
      `mle_unbounded()`-Prüfung bei α ≈ 3 bleibt Pflicht.
- [ ] Kernfragen:
  1. Wird die CCDF gerade? Geht der GoF-p-Wert von 0 weg?
  2. Konvergiert S mit der Energie, und ist er weniger definitionsabhängig als bisher?
  3. Zahl und Größenverteilung der Subkaskaden (Clustering der Schadenszonen
     auf gröberer Skala).
  4. Vergleich mit Sand et al.: S = 1,63 in 3D. Wo liegt 2D?
- [ ] Befund in `results/statistik/fit_zbl/BEFUND_zbl.md`, analog zu den bisherigen.

## Reihenfolge und Aufwand (grob)

| Phase | Inhalt | Aufwand |
|---|---|---|
| 0 | Referenz, Diagnose, `potential.h` | klein |
| 1 | ZBL + adaptiver Zeitschritt + Tests | **groß** (Kern) |
| 2 | Kalibrierung | mittel |
| 3 | Elektronische Bremsung | klein |
| 4 | Pilotläufe, absorbierender Rand | mittel |
| 5 | Produktion | Rechenzeit |
| 6 | Auswertung | mittel (Werkzeuge existieren) |

Alles wird in `cascade_serial.c` umgesetzt. Auf einem Laptop mit 8 Kernen ist
ein Ensemble aus 8 parallelen seriellen Läufen effizienter als MPI innerhalb
eines Laufs (keine Halo-Kommunikation). **MPI und CUDA werden deshalb vorerst
nicht angefasst.** Sie behalten `rep_model = r12` und bleiben so, wie sie im Bericht
stehen.

## Bekannte Grenzen

- 2D: Der Querschnitt ist eine Länge, keine Fläche. Kanalisierung und
  Konnektivität sind anders als in 3D. Ein quantitativer Vergleich mit 3D-S ist
  nur mit Vorbehalt möglich, qualitativ (Potenzgesetz ja/nein) aber sauber.
- Das Federmodell ist kein Wolfram-Potenzial. Die Kalibrierung über E_d macht die
  Energieskala plausibel, nicht die Thermodynamik der Endphase.

## Offene Entscheidungen

Keine für Phase 0–3. Nach Phase 4: Boxgröße und Energien der Produktion, abhängig
vom gemessenen Laufzeitbudget.
