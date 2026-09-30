# Kollisionskaskaden in 2D: Kursprojekt, Paper, Bachelorarbeit

2D-Molekulardynamik von Kollisionskaskaden (Strahlenschaden) in einem
Masse-Feder-Gitter. Das Projekt begann als HPC-Semesterprojekt. Daraus sind
zwei getrennte Vorhaben entstanden.

**Einstieg: [`docs/UEBERBLICK.md`](docs/UEBERBLICK.md).** Dort stehen die ganze
Geschichte, was von den Ergebnissen noch gilt und wo was liegt.

| Ordner | Inhalt | Status |
|---|---|---|
| [`kursprojekt/`](kursprojekt/README.md) | HPC-Semesterprojekt: Bericht, Code seriell/MPI/CUDA, Messungen auf Noctua 2, Clusterstatistik, Physik-Präsentation | eingefroren |
| [`paper/`](paper/README.md) | Weiterentwicklung des Modells (ZBL-Potenzial, Wolfram-Einheiten, elektronische Bremsung), später echte Läufe auf dem PC2 | aktiv |
| [`bachelorarbeit/`](bachelorarbeit/README.md) | Agentenbetriebenes HPC-Labor: Planung, später Auswertung | Planung |
| [`docs/`](docs/README.md) | Überblick, Entscheidungslog, Parametertabelle, Glossar | gilt für alles |

**CI (Forgejo):**
[![Simulation CI](https://forgejo.paultankerfahrer.org/PaulTankerfahrer/parallel_cascades_hpc/actions/workflows/simulation.yml/badge.svg)](https://forgejo.paultankerfahrer.org/PaulTankerfahrer/parallel_cascades_hpc/actions?workflow=simulation.yml)
[![Valgrind](https://forgejo.paultankerfahrer.org/PaulTankerfahrer/parallel_cascades_hpc/actions/workflows/valgrind.yml/badge.svg)](https://forgejo.paultankerfahrer.org/PaulTankerfahrer/parallel_cascades_hpc/actions?workflow=valgrind.yml)

Die CI testet beide Codestände: `kursprojekt/scripts/ci/` (seriell + MPI, altes
Modell) und `paper/scripts/ci/` (aktueller Code mit ZBL und Bremsung).

## Tags

| Tag | Bedeutung |
|---|---|
| `kursprojekt-code` | Modellcode genau wie im Kursprojekt, vor jeder Änderung. Startpunkt für das Agentenlabor |
| `kursprojekt-final` | Ende des Kursprojekts: Bericht, Clusterstatistik, gehaltene Präsentation |

## Lizenz

Code: MIT ([`LICENSE`](LICENSE)). Bericht und Folien: CC BY 4.0.
