# Clauset-Refit der Kaskaden-Clustergroessen

Paket: `powerlaw` 2.0.0 (Alstott et al. 2014), Methodik: Clauset, Shalizi & Newman (2009).
Bootstrap: 300 Resamples, Goodness-of-Fit: 200 synthetische Datensaetze, Seed 12345.

## Hauptergebnis

| Datensatz | n | x_min | n_tail | S | sigma_MLE | 95 %-CI (Bootstrap) | KS D | GoF p |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| E=1300, disp>0,5 / link 1,5 (alte Definition) | 1131 | 1 | 1131 | 1.400 | 0.012 | [1.379, 1.426] | 0.1734 | 0.000 |
| E=1800, disp>0,5 / link 1,5 (alte Definition) | 1100 | 1 | 1100 | 1.354 | 0.011 | [1.333, 1.378] | 0.2013 | 0.000 |
| E=2400, disp>0,5 / link 1,5 (alte Definition) | 675 | 1 | 675 | 1.234 | 0.009 | [1.220, 1.250] | 0.3003 | 0.000 |
| E=1300, disp>1,0 / link 1,0 | 9281 | 2 | 3273 | 2.498 | 0.026 | [2.230, 2.950] | 0.0320 | 0.000 |
| E=1800, disp>1,0 / link 1,0 | 14588 | 2 | 5728 | 2.191 | 0.016 | [2.076, 2.279] | 0.0173 | 0.000 |
| E=2400, disp>1,0 / link 1,0 | 21580 | 1 | 21580 | 1.996 | 0.007 | [1.982, 2.010] | 0.0127 | 0.000 |

## Likelihood-Ratio-Tests (Potenzgesetz vs. Alternative)

R > 0 = Potenzgesetz bevorzugt, R < 0 = Alternative bevorzugt; p > 0.1 = Unterschied statistisch nicht signifikant.
Nicht genestete Alternativen: Vuong-Test mit normiertem R. Genestet (Cutoff): R = rohe Log-Likelihood-Differenz, p aus chi^2.

| Datensatz | Alternative | Test | R | p | Interpretation |
|---|---|---|---:|---:|---|
| E=1300, disp>0,5 / link 1,5 (alte Definition) | Lognormal | Vuong | -11.673 | 0.0000 | Alternative bevorzugt |
| E=1300, disp>0,5 / link 1,5 (alte Definition) | Exponential | Vuong | 31.237 | 0.0000 | Potenzgesetz bevorzugt |
| E=1300, disp>0,5 / link 1,5 (alte Definition) | Stretched Exponential | Vuong | -0.845 | 0.3983 | nicht unterscheidbar |
| E=1300, disp>0,5 / link 1,5 (alte Definition) | Potenzgesetz mit Cutoff | genestet | -31.625 | 0.0000 | Alternative bevorzugt |
| E=1800, disp>0,5 / link 1,5 (alte Definition) | Lognormal | Vuong | -3.988 | 0.0001 | Alternative bevorzugt |
| E=1800, disp>0,5 / link 1,5 (alte Definition) | Exponential | Vuong | 32.929 | 0.0000 | Potenzgesetz bevorzugt |
| E=1800, disp>0,5 / link 1,5 (alte Definition) | Stretched Exponential | Vuong | 0.536 | 0.5918 | nicht unterscheidbar |
| E=1800, disp>0,5 / link 1,5 (alte Definition) | Potenzgesetz mit Cutoff | genestet | -29.228 | 0.0000 | Alternative bevorzugt |
| E=2400, disp>0,5 / link 1,5 (alte Definition) | Lognormal | Vuong | -7.670 | 0.0000 | Alternative bevorzugt |
| E=2400, disp>0,5 / link 1,5 (alte Definition) | Exponential | Vuong | 13.858 | 0.0000 | Potenzgesetz bevorzugt |
| E=2400, disp>0,5 / link 1,5 (alte Definition) | Stretched Exponential | Vuong | -4.825 | 0.0000 | Alternative bevorzugt |
| E=2400, disp>0,5 / link 1,5 (alte Definition) | Potenzgesetz mit Cutoff | genestet | -92.359 | 0.0000 | Alternative bevorzugt |
| E=1300, disp>1,0 / link 1,0 | Lognormal | Vuong | -5.585 | 0.0000 | Alternative bevorzugt |
| E=1300, disp>1,0 / link 1,0 | Exponential | Vuong | 6.213 | 0.0000 | Potenzgesetz bevorzugt |
| E=1300, disp>1,0 / link 1,0 | Stretched Exponential | Vuong | -5.772 | 0.0000 | Alternative bevorzugt |
| E=1300, disp>1,0 / link 1,0 | Potenzgesetz mit Cutoff | genestet | -43.396 | 0.0000 | Alternative bevorzugt |
| E=1800, disp>1,0 / link 1,0 | Lognormal | Vuong | -4.915 | 0.0000 | Alternative bevorzugt |
| E=1800, disp>1,0 / link 1,0 | Exponential | Vuong | 17.347 | 0.0000 | Potenzgesetz bevorzugt |
| E=1800, disp>1,0 / link 1,0 | Stretched Exponential | Vuong | -4.910 | 0.0000 | Alternative bevorzugt |
| E=1800, disp>1,0 / link 1,0 | Potenzgesetz mit Cutoff | genestet | -36.693 | 0.0000 | Alternative bevorzugt |
| E=2400, disp>1,0 / link 1,0 | Lognormal | Vuong | 1.802 | 0.0716 | Potenzgesetz bevorzugt |
| E=2400, disp>1,0 / link 1,0 | Exponential | Vuong | 41.629 | 0.0000 | Potenzgesetz bevorzugt |
| E=2400, disp>1,0 / link 1,0 | Stretched Exponential | Vuong | 12.413 | 0.0000 | Potenzgesetz bevorzugt |
| E=2400, disp>1,0 / link 1,0 | Potenzgesetz mit Cutoff | genestet | -4.298 | 0.0034 | Alternative bevorzugt |

## Kontrollfit bei festem x_min = 4.0 (Wert aus dem Projektbericht)

| Datensatz | S | sigma | n_tail | KS D |
|---|---:|---:|---:|---:|
| E=1300, disp>0,5 / link 1,5 (alte Definition) | 1.278 | 0.013 | 431 | 0.3806 |
| E=1800, disp>0,5 / link 1,5 (alte Definition) | 1.225 | 0.011 | 422 | 0.4305 |
| E=2400, disp>0,5 / link 1,5 (alte Definition) | 1.174 | 0.009 | 351 | 0.5294 |
| E=1300, disp>1,0 / link 1,0 | 2.815 | 0.056 | 1060 | 0.0440 |
| E=1800, disp>1,0 / link 1,0 | 2.278 | 0.027 | 2245 | 0.0320 |
| E=2400, disp>1,0 / link 1,0 | 1.877 | 0.015 | 3484 | 0.0339 |

## x_min-Diagnose (KS-Abstand als Funktion von x_min)


**E=1300, disp>0,5 / link 1,5 (alte Definition)**

| x_min | n_tail | S | KS D |
|---:|---:|---:|---:|
| 1 | 1131 | 1.400 | 0.1734 |
| 2 | 626 | 1.317 | 0.2733 |
| 3 | 482 | 1.283 | 0.3415 |
| 4 | 431 | 1.278 | 0.3806 |
| 5 | 399 | 1.276 | 0.4106 |
| 6 | 376 | 1.274 | 0.4356 |
| 7 | 357 | 1.273 | 0.4587 |
| 8 | 345 | 1.274 | 0.4746 |
| 9 | 337 | 1.277 | 0.4858 |
| 10 | 333 | 1.282 | 0.4912 |
| 12 | 323 | 1.289 | 0.5056 |
| 13 | 320 | 1.294 | 0.5097 |
| 16 | 316 | 1.309 | 0.5136 |
| 19 | 312 | 1.323 | 0.5178 |

**E=1800, disp>0,5 / link 1,5 (alte Definition)**

| x_min | n_tail | S | KS D |
|---:|---:|---:|---:|
| 1 | 1100 | 1.354 | 0.2013 |
| 2 | 586 | 1.256 | 0.3221 |
| 3 | 471 | 1.233 | 0.3885 |
| 4 | 422 | 1.225 | 0.4305 |
| 5 | 393 | 1.221 | 0.4610 |
| 6 | 371 | 1.218 | 0.4876 |
| 7 | 355 | 1.217 | 0.5093 |
| 8 | 347 | 1.218 | 0.5212 |
| 9 | 339 | 1.219 | 0.5337 |
| 11 | 328 | 1.222 | 0.5518 |
| 12 | 321 | 1.222 | 0.5638 |
| 14 | 313 | 1.224 | 0.5781 |
| 16 | 310 | 1.229 | 0.5836 |
| 21 | 306 | 1.242 | 0.5905 |

**E=2400, disp>0,5 / link 1,5 (alte Definition)**

| x_min | n_tail | S | KS D |
|---:|---:|---:|---:|
| 1 | 675 | 1.234 | 0.3003 |
| 2 | 428 | 1.183 | 0.4375 |
| 3 | 372 | 1.174 | 0.4995 |
| 4 | 351 | 1.174 | 0.5294 |
| 5 | 336 | 1.174 | 0.5530 |
| 6 | 329 | 1.176 | 0.5650 |
| 7 | 326 | 1.180 | 0.5707 |
| 8 | 318 | 1.180 | 0.5850 |
| 9 | 313 | 1.181 | 0.5945 |
| 11 | 309 | 1.186 | 0.6021 |
| 16 | 306 | 1.198 | 0.6081 |
| 24 | 303 | 1.214 | 0.6139 |
| 2364 | 299 | 9.983 | 0.6584 |
| 2384 | 297 | 10.648 | 0.6696 |

**E=1300, disp>1,0 / link 1,0**

| x_min | n_tail | S | KS D |
|---:|---:|---:|---:|
| 1 | 9281 | 2.245 | 0.0360 |
| 2 | 3273 | 2.498 | 0.0320 |
| 3 | 1685 | 2.642 | 0.0461 |
| 4 | 1060 | 2.815 | 0.0440 |
| 5 | 720 | 2.964 | 0.0426 |
| 6 | 518 | 3.116 | 0.0301 |
| 7 | 378 | 3.209 | 0.0425 |
| 8 | 285 | 3.294 | 0.0566 |
| 9 | 226 | 3.442 | 0.0737 |
| 10 | 169 | 3.387 | 0.0869 |
| 11 | 149 | 3.711 | 0.1128 |
| 12 | 116 | 3.698 | 0.1277 |
| 13 | 97 | 3.840 | 0.1492 |
| 14 | 82 | 4.004 | 0.1718 |

**E=1800, disp>1,0 / link 1,0**

| x_min | n_tail | S | KS D |
|---:|---:|---:|---:|
| 1 | 14588 | 2.074 | 0.0252 |
| 2 | 5728 | 2.191 | 0.0173 |
| 3 | 3299 | 2.241 | 0.0194 |
| 4 | 2245 | 2.278 | 0.0320 |
| 5 | 1615 | 2.260 | 0.0441 |
| 6 | 1302 | 2.311 | 0.0510 |
| 7 | 1060 | 2.330 | 0.0601 |
| 8 | 903 | 2.373 | 0.0614 |
| 9 | 791 | 2.433 | 0.0592 |
| 10 | 704 | 2.500 | 0.0686 |
| 11 | 625 | 2.551 | 0.0775 |
| 12 | 551 | 2.575 | 0.0871 |
| 14 | 443 | 2.641 | 0.1029 |
| 15 | 407 | 2.699 | 0.1028 |

**E=2400, disp>1,0 / link 1,0**

| x_min | n_tail | S | KS D |
|---:|---:|---:|---:|
| 1 | 21580 | 1.996 | 0.0127 |
| 2 | 8423 | 1.983 | 0.0324 |
| 3 | 4938 | 1.922 | 0.0384 |
| 4 | 3484 | 1.877 | 0.0339 |
| 5 | 2704 | 1.843 | 0.0379 |
| 6 | 2260 | 1.833 | 0.0440 |
| 7 | 1943 | 1.822 | 0.0501 |
| 8 | 1711 | 1.813 | 0.0565 |
| 9 | 1556 | 1.819 | 0.0621 |
| 10 | 1424 | 1.821 | 0.0678 |
| 11 | 1336 | 1.837 | 0.0712 |
| 13 | 1188 | 1.862 | 0.0804 |
| 15 | 1047 | 1.863 | 0.0910 |
| 17 | 957 | 1.884 | 0.0986 |

## Foliensatz-Formulierung

- **E=1300, disp>0,5 / link 1,5 (alte Definition):** S = 1.40 +/- 0.01 (95 %-Bootstrap-CI [1.38, 1.43]), x_min = 1 per KS-Minimierung, Goodness-of-Fit p = 0.00
- **E=1800, disp>0,5 / link 1,5 (alte Definition):** S = 1.35 +/- 0.01 (95 %-Bootstrap-CI [1.33, 1.38]), x_min = 1 per KS-Minimierung, Goodness-of-Fit p = 0.00
- **E=2400, disp>0,5 / link 1,5 (alte Definition):** S = 1.23 +/- 0.01 (95 %-Bootstrap-CI [1.22, 1.25]), x_min = 1 per KS-Minimierung, Goodness-of-Fit p = 0.00
- **E=1300, disp>1,0 / link 1,0:** S = 2.50 +/- 0.03 (95 %-Bootstrap-CI [2.23, 2.95]), x_min = 2 per KS-Minimierung, Goodness-of-Fit p = 0.00
- **E=1800, disp>1,0 / link 1,0:** S = 2.19 +/- 0.02 (95 %-Bootstrap-CI [2.08, 2.28]), x_min = 2 per KS-Minimierung, Goodness-of-Fit p = 0.00
- **E=2400, disp>1,0 / link 1,0:** S = 2.00 +/- 0.01 (95 %-Bootstrap-CI [1.98, 2.01]), x_min = 1 per KS-Minimierung, Goodness-of-Fit p = 0.00
