# Swift-Hohenberg — Multigraph vs Modified Multigraph

Equation: $f_t = (r-1)f - 2\Delta f - \Delta^2 f - f^3$

| Library | True terms |
|---|---|
| Multigraph | f, f³, Δf, Δ²f |
| Modified MG | f, f³, Δf, Δ²f |

## Best per configuration (STRidge)

| Noise | Backend | m | IBP | MG fp/fn | MG E | Mod fp/fn | Mod E | Winner |
|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | 0/0 | 0.0255 | 0/0 | 0.0255 | ⚪ |
| clean | polynomial | 12 | yes | 0/0 | 0.0230 | 0/7 | 1.0000 | 🔴 |
| clean | polynomial | 4 | no | 0/0 | 0.1387 | 0/0 | 0.1387 | ⚪ |
| clean | polynomial | 4 | yes | 0/7 | 1.0000 | 0/7 | 1.0000 | ⚪ |
| clean | polynomial | 8 | no | 0/0 | 0.0249 | 0/0 | 0.0249 | ⚪ |
| clean | polynomial | 8 | yes | 0/0 | 0.1719 | 0/0 | 0.1719 | ⚪ |
| clean | spectral | 12 | no | 0/7 | 1.0000 | 0/0 | 0.0258 | 🟢 |
| clean | spectral | 12 | yes | 0/7 | 1.0000 | 0/7 | 1.0000 | ⚪ |
| clean | spectral | 4 | no | 0/0 | 0.0252 | 0/0 | 0.0252 | ⚪ |
| clean | spectral | 4 | yes | 0/7 | 1.0000 | 0/7 | 1.0000 | ⚪ |
| clean | spectral | 8 | no | 0/0 | 0.0281 | 0/0 | 0.0281 | ⚪ |
| clean | spectral | 8 | yes | 0/0 | 0.1719 | 0/0 | 0.1719 | ⚪ |
| 1% | polynomial | 12 | no | 0/0 | 0.0123 | 0/0 | 0.0123 | ⚪ |
| 1% | polynomial | 12 | yes | 0/0 | 0.0098 | 0/0 | 0.0098 | ⚪ |
| 1% | polynomial | 4 | no | 0/7 | 1.0000 | 0/7 | 1.0000 | ⚪ |
| 1% | polynomial | 4 | yes | 0/7 | 1.0000 | 0/7 | 1.0000 | ⚪ |
| 1% | polynomial | 8 | no | 0/7 | 1.0000 | 0/7 | 1.0000 | ⚪ |
| 1% | polynomial | 8 | yes | 0/7 | 1.0000 | 8/0 | 0.1242 | 🔴 |
| 1% | spectral | 12 | no | 0/0 | 0.0115 | 0/0 | 0.0115 | ⚪ |
| 1% | spectral | 12 | yes | 0/0 | 0.0098 | 0/0 | 0.0098 | ⚪ |
| 1% | spectral | 4 | no | 3/5 | 1.0008 | 0/7 | 1.0000 | 🟢 |
| 1% | spectral | 4 | yes | 3/5 | 1.0008 | 0/7 | 1.0000 | 🟢 |
| 1% | spectral | 8 | no | 0/0 | 0.1039 | 0/0 | 0.1039 | ⚪ |
| 1% | spectral | 8 | yes | 0/0 | 0.1148 | 0/0 | 0.1148 | ⚪ |
| 2% | polynomial | 12 | no | 0/0 | 0.0050 | 0/0 | 0.0050 | ⚪ |
| 2% | polynomial | 12 | yes | 0/0 | 0.0023 | 0/0 | 0.0023 | ⚪ |
| 2% | polynomial | 4 | no | 8/5 | 1.0010 | 0/7 | 1.0000 | 🟢 |
| 2% | polynomial | 4 | yes | 0/7 | 1.0000 | 0/7 | 1.0000 | ⚪ |
| 2% | polynomial | 8 | no | 2/1 | 0.3072 | 3/1 | 0.3790 | 🔴 |
| 2% | polynomial | 8 | yes | 0/7 | 1.0000 | 0/7 | 1.0000 | ⚪ |
| 2% | spectral | 12 | no | 0/0 | 0.0040 | 0/0 | 0.0040 | ⚪ |
| 2% | spectral | 12 | yes | 0/0 | 0.0023 | 0/0 | 0.0023 | ⚪ |
| 2% | spectral | 4 | no | 2/5 | 1.0000 | 0/7 | 1.0000 | 🟥 |
| 2% | spectral | 4 | yes | 2/5 | 1.0000 | 0/7 | 1.0000 | 🟥 |
| 2% | spectral | 8 | no | 0/0 | 0.3083 | 0/0 | 0.3083 | ⚪ |
| 2% | spectral | 8 | yes | 0/0 | 0.3128 | 2/0 | 0.3653 | 🔴 |
| 5% | polynomial | 12 | no | 0/0 | 0.0497 | 0/0 | 0.0497 | ⚪ |
| 5% | polynomial | 12 | yes | 0/0 | 0.0536 | 0/0 | 0.0536 | ⚪ |
| 5% | polynomial | 4 | no | 0/7 | 1.0000 | 0/7 | 1.0000 | ⚪ |
| 5% | polynomial | 4 | yes | 0/7 | 1.0000 | 0/7 | 1.0000 | ⚪ |
| 5% | polynomial | 8 | no | 0/7 | 1.0000 | 0/7 | 1.0000 | ⚪ |
| 5% | polynomial | 8 | yes | 0/7 | 1.0000 | 0/7 | 1.0000 | ⚪ |
| 5% | spectral | 12 | no | 0/0 | 0.0521 | 0/0 | 0.0521 | ⚪ |
| 5% | spectral | 12 | yes | 0/0 | 0.0536 | 0/0 | 0.0536 | ⚪ |
| 5% | spectral | 4 | no | 2/5 | 1.0000 | 0/5 | 1.0000 | 🟢 |
| 5% | spectral | 4 | yes | 2/5 | 1.0000 | 0/5 | 1.0000 | 🟢 |
| 5% | spectral | 8 | no | 2/2 | 0.8715 | 2/5 | 1.0027 | 🔴 |
| 5% | spectral | 8 | yes | 2/2 | 0.8714 | 2/0 | 0.7608 | 🟢 |

## Detailed results

| Noise | Backend | m | IBP | Library | fp/fn | n_act | E | Exact? | d_tol | cond |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | Modified MG | 0/0 | 7 | 0.0255 | ✓ | 5.0 | 42835.5 |
| clean | polynomial | 12 | no | Multigraph | 0/0 | 7 | 0.0255 | ✓ | 2.0 | 28643.3 |
| clean | polynomial | 12 | yes | Modified MG | 0/7 | 0 | 1.0000 |  | 10.0 | 43141.6 |
| clean | polynomial | 12 | yes | Multigraph | 0/0 | 7 | 0.0230 | ✓ | 5.0 | 28691.7 |
| clean | polynomial | 4 | no | Modified MG | 0/0 | 7 | 0.1387 | ✓ | 0.5 | 23433.0 |
| clean | polynomial | 4 | no | Multigraph | 0/0 | 7 | 0.1387 | ✓ | 2.0 | 10430.3 |
| clean | polynomial | 4 | yes | Modified MG | 0/7 | 0 | 1.0000 |  | 5.0 | 27406.6 |
| clean | polynomial | 4 | yes | Multigraph | 0/7 | 0 | 1.0000 |  | 5.0 | 10448.2 |
| clean | polynomial | 8 | no | Modified MG | 0/0 | 7 | 0.0249 | ✓ | 0.5 | 27385.2 |
| clean | polynomial | 8 | no | Multigraph | 0/0 | 7 | 0.0249 | ✓ | 2.0 | 15952.5 |
| clean | polynomial | 8 | yes | Modified MG | 0/0 | 7 | 0.1719 | ✓ | 2.0 | 27810.9 |
| clean | polynomial | 8 | yes | Multigraph | 0/0 | 7 | 0.1719 | ✓ | 5.0 | 15962.4 |
| clean | spectral | 12 | no | Modified MG | 0/0 | 7 | 0.0258 | ✓ | 5.0 | 42761.4 |
| clean | spectral | 12 | no | Multigraph | 0/7 | 0 | 1.0000 |  | 10.0 | 28641.5 |
| clean | spectral | 12 | yes | Modified MG | 0/7 | 0 | 1.0000 |  | 10.0 | 42767.6 |
| clean | spectral | 12 | yes | Multigraph | 0/7 | 0 | 1.0000 |  | 10.0 | 28695.1 |
| clean | spectral | 4 | no | Modified MG | 0/0 | 7 | 0.0252 | ✓ | 0.5 | 23276.8 |
| clean | spectral | 4 | no | Multigraph | 0/0 | 7 | 0.0252 | ✓ | 1.0 | 10372.4 |
| clean | spectral | 4 | yes | Modified MG | 0/7 | 0 | 1.0000 |  | 2.0 | 27461.9 |
| clean | spectral | 4 | yes | Multigraph | 0/7 | 0 | 1.0000 |  | 5.0 | 10380.5 |
| clean | spectral | 8 | no | Modified MG | 0/0 | 7 | 0.0281 | ✓ | 0.2 | 27294.5 |
| clean | spectral | 8 | no | Multigraph | 0/0 | 7 | 0.0281 | ✓ | 2.0 | 15893.5 |
| clean | spectral | 8 | yes | Modified MG | 0/0 | 7 | 0.1719 | ✓ | 2.0 | 27481.0 |
| clean | spectral | 8 | yes | Multigraph | 0/0 | 7 | 0.1719 | ✓ | 5.0 | 15897.9 |
| 1% | polynomial | 12 | no | Modified MG | 0/0 | 7 | 0.0123 | ✓ | 2.0 | 37593.6 |
| 1% | polynomial | 12 | no | Multigraph | 0/0 | 7 | 0.0123 | ✓ | 1.0 | 27057.2 |
| 1% | polynomial | 12 | yes | Modified MG | 0/0 | 7 | 0.0098 | ✓ | 2.0 | 38241.2 |
| 1% | polynomial | 12 | yes | Multigraph | 0/0 | 7 | 0.0098 | ✓ | 5.0 | 27091.5 |
| 1% | polynomial | 4 | no | Modified MG | 0/7 | 0 | 1.0000 |  | 10.0 | 20650.7 |
| 1% | polynomial | 4 | no | Multigraph | 0/7 | 0 | 1.0000 |  | 2.0 | 9335.6 |
| 1% | polynomial | 4 | yes | Modified MG | 0/7 | 0 | 1.0000 |  | 2.0 | 7676.0 |
| 1% | polynomial | 4 | yes | Multigraph | 0/7 | 0 | 1.0000 |  | 10.0 | 9290.9 |
| 1% | polynomial | 8 | no | Modified MG | 0/7 | 0 | 1.0000 |  | 10.0 | 24037.6 |
| 1% | polynomial | 8 | no | Multigraph | 0/7 | 0 | 1.0000 |  | 10.0 | 14472.2 |
| 1% | polynomial | 8 | yes | Modified MG | 8/0 | 15 | 0.1242 |  | 2.0 | 23417.9 |
| 1% | polynomial | 8 | yes | Multigraph | 0/7 | 0 | 1.0000 |  | 10.0 | 14477.3 |
| 1% | spectral | 12 | no | Modified MG | 0/0 | 7 | 0.0115 | ✓ | 0.5 | 10410.6 |
| 1% | spectral | 12 | no | Multigraph | 0/0 | 7 | 0.0115 | ✓ | 1.0 | 4693.5 |
| 1% | spectral | 12 | yes | Modified MG | 0/0 | 7 | 0.0098 | ✓ | 2.0 | 10067.7 |
| 1% | spectral | 12 | yes | Multigraph | 0/0 | 7 | 0.0098 | ✓ | 0.5 | 4693.5 |
| 1% | spectral | 4 | no | Modified MG | 0/7 | 0 | 1.0000 |  | 10.0 | 1712.8 |
| 1% | spectral | 4 | no | Multigraph | 3/5 | 5 | 1.0008 |  | 5.0 | 892.8 |
| 1% | spectral | 4 | yes | Modified MG | 0/7 | 0 | 1.0000 |  | 1.0 | 2223.0 |
| 1% | spectral | 4 | yes | Multigraph | 3/5 | 5 | 1.0008 |  | 5.0 | 893.7 |
| 1% | spectral | 8 | no | Modified MG | 0/0 | 7 | 0.1039 | ✓ | 1.0 | 3643.0 |
| 1% | spectral | 8 | no | Multigraph | 0/0 | 7 | 0.1039 | ✓ | 1.0 | 1600.9 |
| 1% | spectral | 8 | yes | Modified MG | 0/0 | 7 | 0.1148 | ✓ | 2.0 | 4180.4 |
| 1% | spectral | 8 | yes | Multigraph | 0/0 | 7 | 0.1148 | ✓ | 1.0 | 1601.9 |
| 2% | polynomial | 12 | no | Modified MG | 0/0 | 7 | 0.0050 | ✓ | 2.0 | 34554.7 |
| 2% | polynomial | 12 | no | Multigraph | 0/0 | 7 | 0.0050 | ✓ | 2.0 | 25224.6 |
| 2% | polynomial | 12 | yes | Modified MG | 0/0 | 7 | 0.0023 | ✓ | 2.0 | 34914.8 |
| 2% | polynomial | 12 | yes | Multigraph | 0/0 | 7 | 0.0023 | ✓ | 5.0 | 25249.1 |
| 2% | polynomial | 4 | no | Modified MG | 0/7 | 0 | 1.0000 |  | 10.0 | 17745.0 |
| 2% | polynomial | 4 | no | Multigraph | 8/5 | 10 | 1.0010 |  | 2.0 | 8016.7 |
| 2% | polynomial | 4 | yes | Modified MG | 0/7 | 0 | 1.0000 |  | 10.0 | 5680.9 |
| 2% | polynomial | 4 | yes | Multigraph | 0/7 | 0 | 1.0000 |  | 5.0 | 8011.6 |
| 2% | polynomial | 8 | no | Modified MG | 3/1 | 9 | 0.3790 |  | 5.0 | 22169.0 |
| 2% | polynomial | 8 | no | Multigraph | 2/1 | 8 | 0.3072 |  | 5.0 | 13410.6 |
| 2% | polynomial | 8 | yes | Modified MG | 0/7 | 0 | 1.0000 |  | 5.0 | 20533.1 |
| 2% | polynomial | 8 | yes | Multigraph | 0/7 | 0 | 1.0000 |  | 10.0 | 13414.0 |
| 2% | spectral | 12 | no | Modified MG | 0/0 | 7 | 0.0040 | ✓ | 0.5 | 7715.4 |
| 2% | spectral | 12 | no | Multigraph | 0/0 | 7 | 0.0040 | ✓ | 1.0 | 3350.2 |
| 2% | spectral | 12 | yes | Modified MG | 0/0 | 7 | 0.0023 | ✓ | 0.5 | 9090.4 |
| 2% | spectral | 12 | yes | Multigraph | 0/0 | 7 | 0.0023 | ✓ | 1.0 | 3350.1 |
| 2% | spectral | 4 | no | Modified MG | 0/7 | 0 | 1.0000 |  | 2.0 | 1384.1 |
| 2% | spectral | 4 | no | Multigraph | 2/5 | 4 | 1.0000 |  | 1.0 | 692.9 |
| 2% | spectral | 4 | yes | Modified MG | 0/7 | 0 | 1.0000 |  | 1.0 | 1720.3 |
| 2% | spectral | 4 | yes | Multigraph | 2/5 | 4 | 1.0000 |  | 1.0 | 691.0 |
| 2% | spectral | 8 | no | Modified MG | 0/0 | 7 | 0.3083 | ✓ | 2.0 | 2585.0 |
| 2% | spectral | 8 | no | Multigraph | 0/0 | 7 | 0.3083 | ✓ | 2.0 | 1432.1 |
| 2% | spectral | 8 | yes | Modified MG | 2/0 | 9 | 0.3653 |  | 1.0 | 3561.6 |
| 2% | spectral | 8 | yes | Multigraph | 0/0 | 7 | 0.3128 | ✓ | 2.0 | 1432.2 |
| 5% | polynomial | 12 | no | Modified MG | 0/0 | 7 | 0.0497 | ✓ | 5.0 | 13863.7 |
| 5% | polynomial | 12 | no | Multigraph | 0/0 | 7 | 0.0497 | ✓ | 2.0 | 13259.2 |
| 5% | polynomial | 12 | yes | Modified MG | 0/0 | 7 | 0.0536 | ✓ | 1.0 | 14719.0 |
| 5% | polynomial | 12 | yes | Multigraph | 0/0 | 7 | 0.0536 | ✓ | 2.0 | 13263.5 |
| 5% | polynomial | 4 | no | Modified MG | 0/7 | 0 | 1.0000 |  | 2.0 | 6688.2 |
| 5% | polynomial | 4 | no | Multigraph | 0/7 | 0 | 1.0000 |  | 2.0 | 3051.7 |
| 5% | polynomial | 4 | yes | Modified MG | 0/7 | 0 | 1.0000 |  | 1.0 | 2721.5 |
| 5% | polynomial | 4 | yes | Multigraph | 0/7 | 0 | 1.0000 |  | 2.0 | 3037.9 |
| 5% | polynomial | 8 | no | Modified MG | 0/7 | 0 | 1.0000 |  | 5.0 | 9148.8 |
| 5% | polynomial | 8 | no | Multigraph | 0/7 | 0 | 1.0000 |  | 10.0 | 6492.0 |
| 5% | polynomial | 8 | yes | Modified MG | 0/7 | 0 | 1.0000 |  | 10.0 | 8809.0 |
| 5% | polynomial | 8 | yes | Multigraph | 0/7 | 0 | 1.0000 |  | 10.0 | 6488.1 |
| 5% | spectral | 12 | no | Modified MG | 0/0 | 7 | 0.0521 | ✓ | 2.0 | 4162.5 |
| 5% | spectral | 12 | no | Multigraph | 0/0 | 7 | 0.0521 | ✓ | 5.0 | 4258.3 |
| 5% | spectral | 12 | yes | Modified MG | 0/0 | 7 | 0.0536 | ✓ | 1.0 | 6436.1 |
| 5% | spectral | 12 | yes | Multigraph | 0/0 | 7 | 0.0536 | ✓ | 5.0 | 4258.1 |
| 5% | spectral | 4 | no | Modified MG | 0/5 | 2 | 1.0000 |  | 1.0 | 1535.1 |
| 5% | spectral | 4 | no | Multigraph | 2/5 | 4 | 1.0000 |  | 0.5 | 744.2 |
| 5% | spectral | 4 | yes | Modified MG | 0/5 | 2 | 1.0000 |  | 1.0 | 1049.4 |
| 5% | spectral | 4 | yes | Multigraph | 2/5 | 4 | 1.0000 |  | 0.5 | 741.0 |
| 5% | spectral | 8 | no | Modified MG | 2/5 | 4 | 1.0027 |  | 5.0 | 2122.1 |
| 5% | spectral | 8 | no | Multigraph | 2/2 | 7 | 0.8715 |  | 5.0 | 1578.6 |
| 5% | spectral | 8 | yes | Modified MG | 2/0 | 9 | 0.7608 |  | 2.0 | 2746.5 |
| 5% | spectral | 8 | yes | Multigraph | 2/2 | 7 | 0.8714 |  | 5.0 | 1579.0 |

---
*Source: `mg_comparison_20260629_165010_full.csv`*
