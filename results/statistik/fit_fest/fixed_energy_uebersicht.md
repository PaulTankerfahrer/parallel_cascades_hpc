# Fixed-Energy-Ensemble: Uebersicht

Feste PKA-Energie, variiert werden nur Einschlagort, Winkel und Seed. Die Einschlagparameter haengen nur von der Task-ID ab, sind also bei allen drei Energien identisch (gepaarter Vergleich).

`R` < 0 heisst: die Alternative beschreibt die Daten besser als das Potenzgesetz.

⚠ = das powerlaw-Paket sucht `alpha` nur in [0, 3] und haette 3.000 gemeldet; angegeben ist die freie MLE.
`n. k.` = der Cutoff-Fit ist nicht konvergiert (das genestete Modell kann nie schlechter sein als das Potenzgesetz, ein positives R ist daher ein numerischer Fehlschlag).


## E = 1300  (im Mittel 474 verlagerte Atome pro Lauf)

| Clusterdefinition | Cluster | pro Lauf | groesster | S | sigma | x_min | n_tail | KS-D | R Lognormal | R Cutoff |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| disp>0.5 / link 1.0 | 15826 | 52.8 | 632 | 1.980 | 0.008 | 1 | 15826 | 0.0321 | -0.5 | -0.0 |
| disp>0.5 / link 1.5 | 1131 | 3.8 | 751 | 1.400 | 0.012 | 1 | 1131 | 0.1734 | -11.7 | -31.6 |
| disp>1.0 / link 1.0 | 9281 | 30.9 | 43 | 2.498 | 0.026 | 2 | 3273 | 0.0320 | -5.6 | -43.4 |
| disp>1.0 / link 1.5 | 2005 | 6.7 | 96 | 1.565 | 0.013 | 1 | 2005 | 0.0691 | -9.9 | -149.2 |
| disp>1.5 / link 1.0 | 6481 | 21.6 | 16 | 2.778 | 0.022 | 1 | 6481 | 0.0272 | -7.8 | -80.3 |
| disp>2.0 / link 1.0 | 3770 | 12.6 | 9 | 3.296 ⚠ | 0.033 | 1 | 3770 | 0.0259 | -7.5 | -49.2 |

## E = 1800  (im Mittel 1399 verlagerte Atome pro Lauf)

| Clusterdefinition | Cluster | pro Lauf | groesster | S | sigma | x_min | n_tail | KS-D | R Lognormal | R Cutoff |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| disp>0.5 / link 1.0 | 39522 | 131.7 | 1421 | 2.098 | 0.006 | 1 | 39522 | 0.0229 | +7.2 | n. k. |
| disp>0.5 / link 1.5 | 1100 | 3.7 | 1778 | 1.354 | 0.011 | 1 | 1100 | 0.2013 | -4.0 | -29.2 |
| disp>1.0 / link 1.0 | 14588 | 48.6 | 113 | 2.191 | 0.016 | 2 | 5728 | 0.0173 | -4.9 | -36.7 |
| disp>1.0 / link 1.5 | 1905 | 6.3 | 205 | 1.499 | 0.011 | 1 | 1905 | 0.0669 | -8.0 | -82.8 |
| disp>1.5 / link 1.0 | 11430 | 38.1 | 28 | 2.506 | 0.014 | 1 | 11430 | 0.0306 | -9.9 | -132.6 |
| disp>2.0 / link 1.0 | 7088 | 23.6 | 13 | 3.020 ⚠ | 0.024 | 1 | 7088 | 0.0149 | -6.6 | -55.2 |

## E = 2400  (im Mittel 2646 verlagerte Atome pro Lauf)

| Clusterdefinition | Cluster | pro Lauf | groesster | S | sigma | x_min | n_tail | KS-D | R Lognormal | R Cutoff |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| disp>0.5 / link 1.0 | 55595 | 185.3 | 2596 | 2.251 | 0.005 | 1 | 55595 | 0.0179 | +6.6 | n. k. |
| disp>0.5 / link 1.5 | 675 | 2.2 | 2906 | 1.234 | 0.009 | 1 | 675 | 0.3003 | -7.7 | -92.4 |
| disp>1.0 / link 1.0 | 21580 | 71.9 | 330 | 1.996 | 0.007 | 1 | 21580 | 0.0127 | +1.8 | -4.3 |
| disp>1.0 / link 1.5 | 1415 | 4.7 | 528 | 1.454 | 0.012 | 1 | 1415 | 0.1366 | -3.9 | -25.5 |
| disp>1.5 / link 1.0 | 18150 | 60.5 | 58 | 2.279 | 0.009 | 1 | 18150 | 0.0215 | -8.1 | -90.4 |
| disp>2.0 / link 1.0 | 12318 | 41.1 | 15 | 2.769 | 0.016 | 1 | 12318 | 0.0212 | -8.4 | -96.5 |


## Dieselbe Clusterdefinition ueber die Energien

- **disp>0.5 / link 1.0**: S = 1.98 (E=1300), 2.10 (E=1800), 2.25 (E=2400)  -- Spannweite 0.27
- **disp>0.5 / link 1.5**: S = 1.40 (E=1300), 1.35 (E=1800), 1.23 (E=2400)  -- Spannweite 0.17
- **disp>1.0 / link 1.0**: S = 2.50 (E=1300), 2.19 (E=1800), 2.00 (E=2400)  -- Spannweite 0.50
- **disp>1.0 / link 1.5**: S = 1.57 (E=1300), 1.50 (E=1800), 1.45 (E=2400)  -- Spannweite 0.11
- **disp>1.5 / link 1.0**: S = 2.78 (E=1300), 2.51 (E=1800), 2.28 (E=2400)  -- Spannweite 0.50
- **disp>2.0 / link 1.0**: S = 3.30 (E=1300), 3.02 (E=1800), 2.77 (E=2400)  -- Spannweite 0.53
