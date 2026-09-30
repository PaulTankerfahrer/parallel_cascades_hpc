# Clauset-Refit der Kaskaden-Clustergroessen

Paket: `powerlaw` 2.0.0 (Alstott et al. 2014), Methodik: Clauset, Shalizi & Newman (2009).
Bootstrap: 1000 Resamples, Goodness-of-Fit: 1000 synthetische Datensaetze, Seed 12345.

## Hauptergebnis

| Datensatz | n | x_min | n_tail | S | sigma_MLE | 95 %-CI (Bootstrap) | KS D | GoF p |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| Mit Healing | 2391 | 1 | 2391 | 1.409 | 0.008 | [1.393, 1.426] | 0.0414 | 0.000 |
| Ohne Healing | 2335 | 1 | 2335 | 1.390 | 0.008 | [1.374, 1.406] | 0.0514 | 0.000 |

## Likelihood-Ratio-Tests (Potenzgesetz vs. Alternative)

R > 0 = Potenzgesetz bevorzugt, R < 0 = Alternative bevorzugt; p > 0.1 = Unterschied statistisch nicht signifikant.
Nicht genestete Alternativen: Vuong-Test mit normiertem R. Genestet (Cutoff): R = rohe Log-Likelihood-Differenz, p aus chi^2.

| Datensatz | Alternative | Test | R | p | Interpretation |
|---|---|---|---:|---:|---|
| Mit Healing | Lognormal | Vuong | -9.498 | 0.0000 | Alternative bevorzugt |
| Mit Healing | Exponential | Vuong | 50.209 | 0.0000 | Potenzgesetz bevorzugt |
| Mit Healing | Stretched Exponential | Vuong | -2.445 | 0.0145 | Alternative bevorzugt |
| Mit Healing | Potenzgesetz mit Cutoff | genestet | -44.997 | 0.0000 | Alternative bevorzugt |
| Ohne Healing | Lognormal | Vuong | -11.719 | 0.0000 | Alternative bevorzugt |
| Ohne Healing | Exponential | Vuong | 49.740 | 0.0000 | Potenzgesetz bevorzugt |
| Ohne Healing | Stretched Exponential | Vuong | -4.872 | 0.0000 | Alternative bevorzugt |
| Ohne Healing | Potenzgesetz mit Cutoff | genestet | -53.965 | 0.0000 | Alternative bevorzugt |

## Kontrollfit bei festem x_min = 4.0 (Wert aus dem Projektbericht)

| Datensatz | S | sigma | n_tail | KS D |
|---|---:|---:|---:|---:|
| Mit Healing | 1.391 | 0.012 | 1131 | 0.0759 |
| Ohne Healing | 1.361 | 0.011 | 1117 | 0.0858 |

## x_min-Diagnose (KS-Abstand als Funktion von x_min)


**Mit Healing**

| x_min | n_tail | S | KS D |
|---:|---:|---:|---:|
| 1 | 2391 | 1.409 | 0.0414 |
| 2 | 1513 | 1.381 | 0.0556 |
| 3 | 1247 | 1.379 | 0.0681 |
| 4 | 1131 | 1.391 | 0.0759 |
| 5 | 1031 | 1.393 | 0.0839 |
| 6 | 945 | 1.389 | 0.0904 |
| 7 | 884 | 1.389 | 0.0964 |
| 9 | 818 | 1.400 | 0.1074 |
| 10 | 783 | 1.400 | 0.1123 |
| 12 | 719 | 1.396 | 0.1214 |
| 13 | 702 | 1.400 | 0.1252 |
| 16 | 647 | 1.402 | 0.1363 |
| 18 | 626 | 1.408 | 0.1419 |
| 21 | 594 | 1.413 | 0.1500 |

**Ohne Healing**

| x_min | n_tail | S | KS D |
|---:|---:|---:|---:|
| 1 | 2335 | 1.390 | 0.0514 |
| 2 | 1478 | 1.356 | 0.0662 |
| 3 | 1230 | 1.353 | 0.0805 |
| 4 | 1117 | 1.361 | 0.0858 |
| 5 | 1025 | 1.363 | 0.0929 |
| 6 | 976 | 1.372 | 0.0944 |
| 7 | 934 | 1.379 | 0.0986 |
| 9 | 871 | 1.392 | 0.1095 |
| 10 | 832 | 1.391 | 0.1144 |
| 12 | 771 | 1.390 | 0.1232 |
| 14 | 727 | 1.392 | 0.1310 |
| 16 | 693 | 1.394 | 0.1379 |
| 19 | 643 | 1.392 | 0.1483 |
| 22 | 613 | 1.397 | 0.1559 |

## Foliensatz-Formulierung

- **Mit Healing:** S = 1.41 +/- 0.01 (95 %-Bootstrap-CI [1.39, 1.43]), x_min = 1 per KS-Minimierung, Goodness-of-Fit p = 0.00
- **Ohne Healing:** S = 1.39 +/- 0.01 (95 %-Bootstrap-CI [1.37, 1.41]), x_min = 1 per KS-Minimierung, Goodness-of-Fit p = 0.00
