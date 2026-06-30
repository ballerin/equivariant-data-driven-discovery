# KS (weak-chaotic) — Multigraph vs Modified Multigraph

Equation: $f_t = -\Delta f - \Delta^2 f - \|\nabla f\|^2/2$

| Library | True terms |
|---|---|
| Multigraph | Δf, ‖∇f‖², Δ²f |
| Modified MG | Δf, ‖∇f‖², Δ²f |

## Best per configuration (STRidge)

| Noise | Backend | m | IBP | MG fp/fn | MG E | Mod fp/fn | Mod E | Winner |
|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | 9/0 | 0.0311 | 2/0 | 0.0067 | 🟢 |
| clean | polynomial | 12 | yes | 10/0 | 0.0337 | 1/0 | 0.0296 | 🟢 |
| clean | polynomial | 4 | no | 0/0 | 0.0080 | 0/0 | 0.0080 | ⚪ |
| clean | polynomial | 4 | yes | 37/0 | 121.5561 | 36/0 | 91.8140 | 🟢 |
| clean | polynomial | 8 | no | 0/0 | 0.0040 | 0/0 | 0.0040 | ⚪ |
| clean | polynomial | 8 | yes | 1/0 | 0.7472 | 1/0 | 0.7472 | ⚪ |
| clean | spectral | 12 | no | 11/0 | 0.0178 | 2/0 | 0.0067 | 🟢 |
| clean | spectral | 12 | yes | 10/0 | 0.0336 | 1/0 | 0.0297 | 🟢 |
| clean | spectral | 4 | no | 0/0 | 0.0040 | 0/0 | 0.0040 | ⚪ |
| clean | spectral | 4 | yes | 37/0 | 121.3384 | 36/0 | 91.7938 | 🟢 |
| clean | spectral | 8 | no | 0/0 | 0.0031 | 0/0 | 0.0031 | ⚪ |
| clean | spectral | 8 | yes | 1/0 | 0.7470 | 1/0 | 0.7470 | ⚪ |
| 2% | polynomial | 12 | no | 32/0 | 0.2283 | 28/0 | 0.2281 | 🟢 |
| 2% | polynomial | 12 | yes | 32/0 | 0.2304 | 37/0 | 0.3724 | 🔴 |
| 2% | polynomial | 4 | no | 11/5 | 1.0012 | 11/5 | 1.0012 | 🟥 |
| 2% | polynomial | 4 | yes | 12/5 | 1.0012 | 31/5 | 1.0009 | 🔴 |
| 2% | polynomial | 8 | no | 32/0 | 0.7626 | 37/0 | 0.7626 | 🔴 |
| 2% | polynomial | 8 | yes | 32/3 | 0.7814 | 33/0 | 0.8621 | 🟢 |
| 2% | spectral | 12 | no | 30/0 | 0.7504 | 35/0 | 0.7484 | 🔴 |
| 2% | spectral | 12 | yes | 30/0 | 0.7502 | 35/3 | 0.8673 | 🔴 |
| 2% | spectral | 4 | no | 37/0 | 1.0001 | 37/0 | 1.0001 | 🟩 |
| 2% | spectral | 4 | yes | 8/5 | 1.0000 | 31/5 | 1.0001 | 🔴 |
| 2% | spectral | 8 | no | 29/5 | 0.9360 | 34/5 | 0.9375 | 🔴 |
| 2% | spectral | 8 | yes | 29/5 | 0.9360 | 32/5 | 0.9756 | 🔴 |
| 5% | polynomial | 12 | no | 30/5 | 0.8089 | 35/5 | 0.8095 | 🔴 |
| 5% | polynomial | 12 | yes | 30/5 | 0.8089 | 32/0 | 0.8866 | 🟢 |
| 5% | polynomial | 4 | no | 11/5 | 1.0004 | 11/5 | 1.0004 | 🟥 |
| 5% | polynomial | 4 | yes | 11/5 | 1.0004 | 12/5 | 1.0005 | 🔴 |
| 5% | polynomial | 8 | no | 26/5 | 0.9734 | 22/5 | 0.9734 | 🟢 |
| 5% | polynomial | 8 | yes | 31/0 | 0.9748 | 36/5 | 0.9921 | 🔴 |
| 5% | spectral | 12 | no | 37/0 | 0.9991 | 35/3 | 0.9889 | 🔴 |
| 5% | spectral | 12 | yes | 37/0 | 0.9991 | 34/0 | 1.0026 | 🟢 |
| 5% | spectral | 4 | no | 8/5 | 1.0000 | 4/7 | 1.0000 | 🟢 |
| 5% | spectral | 4 | yes | 8/5 | 1.0000 | 4/7 | 1.0000 | 🟢 |
| 5% | spectral | 8 | no | 28/3 | 0.9961 | 33/5 | 0.9976 | 🔴 |
| 5% | spectral | 8 | yes | 28/3 | 0.9961 | 37/0 | 0.9993 | 🔴 |

## Detailed results

| Noise | Backend | m | IBP | Library | fp/fn | n_act | E | Exact? | d_tol | cond |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | Modified MG | 2/0 | 9 | 0.0067 |  | 10.0 | 33163493.6 |
| clean | polynomial | 12 | no | Multigraph | 9/0 | 16 | 0.0311 |  | 10.0 | 34063590.5 |
| clean | polynomial | 12 | yes | Modified MG | 1/0 | 8 | 0.0296 |  | 10.0 | 33238548.3 |
| clean | polynomial | 12 | yes | Multigraph | 10/0 | 17 | 0.0337 |  | 10.0 | 34199127.9 |
| clean | polynomial | 4 | no | Modified MG | 0/0 | 7 | 0.0080 | ✓ | 5.0 | 15581625.2 |
| clean | polynomial | 4 | no | Multigraph | 0/0 | 7 | 0.0080 | ✓ | 5.0 | 15564349.4 |
| clean | polynomial | 4 | yes | Modified MG | 36/0 | 43 | 91.8140 |  | 5.0 | 378560254.5 |
| clean | polynomial | 4 | yes | Multigraph | 37/0 | 44 | 121.5561 |  | 0.0001 | 59112570.3 |
| clean | polynomial | 8 | no | Modified MG | 0/0 | 7 | 0.0040 | ✓ | 10.0 | 21242633.9 |
| clean | polynomial | 8 | no | Multigraph | 0/0 | 7 | 0.0040 | ✓ | 10.0 | 21542667.8 |
| clean | polynomial | 8 | yes | Modified MG | 1/0 | 8 | 0.7472 |  | 10.0 | 21996367.3 |
| clean | polynomial | 8 | yes | Multigraph | 1/0 | 8 | 0.7472 |  | 10.0 | 22378946.1 |
| clean | spectral | 12 | no | Modified MG | 2/0 | 9 | 0.0067 |  | 10.0 | 33105025.5 |
| clean | spectral | 12 | no | Multigraph | 11/0 | 18 | 0.0178 |  | 10.0 | 34013918.5 |
| clean | spectral | 12 | yes | Modified MG | 1/0 | 8 | 0.0297 |  | 10.0 | 33225928.9 |
| clean | spectral | 12 | yes | Multigraph | 10/0 | 17 | 0.0336 |  | 10.0 | 34140137.1 |
| clean | spectral | 4 | no | Modified MG | 0/0 | 7 | 0.0040 | ✓ | 5.0 | 15569745.3 |
| clean | spectral | 4 | no | Multigraph | 0/0 | 7 | 0.0040 | ✓ | 5.0 | 15559622.6 |
| clean | spectral | 4 | yes | Modified MG | 36/0 | 43 | 91.7938 |  | 5.0 | 378328508.5 |
| clean | spectral | 4 | yes | Multigraph | 37/0 | 44 | 121.3384 |  | 0.0001 | 59713062.9 |
| clean | spectral | 8 | no | Modified MG | 0/0 | 7 | 0.0031 | ✓ | 10.0 | 21210481.5 |
| clean | spectral | 8 | no | Multigraph | 0/0 | 7 | 0.0031 | ✓ | 10.0 | 21531316.4 |
| clean | spectral | 8 | yes | Modified MG | 1/0 | 8 | 0.7470 |  | 10.0 | 21985357.5 |
| clean | spectral | 8 | yes | Multigraph | 1/0 | 8 | 0.7470 |  | 10.0 | 22310889.7 |
| 2% | polynomial | 12 | no | Modified MG | 28/0 | 35 | 0.2281 |  | 5.0 | 31917642.3 |
| 2% | polynomial | 12 | no | Multigraph | 32/0 | 39 | 0.2283 |  | 5.0 | 32039332.3 |
| 2% | polynomial | 12 | yes | Modified MG | 37/0 | 44 | 0.3724 |  | 10.0 | 30679279.8 |
| 2% | polynomial | 12 | yes | Multigraph | 32/0 | 39 | 0.2304 |  | 10.0 | 32035475.0 |
| 2% | polynomial | 4 | no | Modified MG | 11/5 | 13 | 1.0012 |  | 10.0 | 18300302.2 |
| 2% | polynomial | 4 | no | Multigraph | 11/5 | 13 | 1.0012 |  | 10.0 | 18300991.0 |
| 2% | polynomial | 4 | yes | Modified MG | 31/5 | 33 | 1.0009 |  | 10.0 | 87630380.8 |
| 2% | polynomial | 4 | yes | Multigraph | 12/5 | 14 | 1.0012 |  | 2.0 | 18291220.6 |
| 2% | polynomial | 8 | no | Modified MG | 37/0 | 44 | 0.7626 |  | 10.0 | 24684616.5 |
| 2% | polynomial | 8 | no | Multigraph | 32/0 | 39 | 0.7626 |  | 10.0 | 24727456.6 |
| 2% | polynomial | 8 | yes | Modified MG | 33/0 | 40 | 0.8621 |  | 10.0 | 24030483.1 |
| 2% | polynomial | 8 | yes | Multigraph | 32/3 | 36 | 0.7814 |  | 10.0 | 24702480.8 |
| 2% | spectral | 12 | no | Modified MG | 35/0 | 42 | 0.7484 |  | 0.5 | 43722243.4 |
| 2% | spectral | 12 | no | Multigraph | 30/0 | 37 | 0.7504 |  | 10.0 | 45377815.8 |
| 2% | spectral | 12 | yes | Modified MG | 35/3 | 39 | 0.8673 |  | 5.0 | 54338651.1 |
| 2% | spectral | 12 | yes | Multigraph | 30/0 | 37 | 0.7502 |  | 10.0 | 45383276.7 |
| 2% | spectral | 4 | no | Modified MG | 37/0 | 44 | 1.0001 |  | 0.0001 | 34077902.0 |
| 2% | spectral | 4 | no | Multigraph | 37/0 | 44 | 1.0001 |  | 0.0001 | 21474809.9 |
| 2% | spectral | 4 | yes | Modified MG | 31/5 | 33 | 1.0001 |  | 10.0 | 109778594.8 |
| 2% | spectral | 4 | yes | Multigraph | 8/5 | 10 | 1.0000 |  | 10.0 | 21421560.4 |
| 2% | spectral | 8 | no | Modified MG | 34/5 | 36 | 0.9375 |  | 10.0 | 32234645.8 |
| 2% | spectral | 8 | no | Multigraph | 29/5 | 31 | 0.9360 |  | 10.0 | 33213146.7 |
| 2% | spectral | 8 | yes | Modified MG | 32/5 | 34 | 0.9756 |  | 10.0 | 38364902.9 |
| 2% | spectral | 8 | yes | Multigraph | 29/5 | 31 | 0.9360 |  | 10.0 | 33218011.5 |
| 5% | polynomial | 12 | no | Modified MG | 35/5 | 37 | 0.8095 |  | 5.0 | 36156961.9 |
| 5% | polynomial | 12 | no | Multigraph | 30/5 | 32 | 0.8089 |  | 5.0 | 36171050.9 |
| 5% | polynomial | 12 | yes | Modified MG | 32/0 | 39 | 0.8866 |  | 10.0 | 36004061.2 |
| 5% | polynomial | 12 | yes | Multigraph | 30/5 | 32 | 0.8089 |  | 5.0 | 36169050.1 |
| 5% | polynomial | 4 | no | Modified MG | 11/5 | 13 | 1.0004 |  | 10.0 | 19026937.1 |
| 5% | polynomial | 4 | no | Multigraph | 11/5 | 13 | 1.0004 |  | 10.0 | 18973257.7 |
| 5% | polynomial | 4 | yes | Modified MG | 12/5 | 14 | 1.0005 |  | 10.0 | 97575906.7 |
| 5% | polynomial | 4 | yes | Multigraph | 11/5 | 13 | 1.0004 |  | 5.0 | 18941693.6 |
| 5% | polynomial | 8 | no | Modified MG | 22/5 | 24 | 0.9734 |  | 10.0 | 27932865.4 |
| 5% | polynomial | 8 | no | Multigraph | 26/5 | 28 | 0.9734 |  | 10.0 | 27941601.0 |
| 5% | polynomial | 8 | yes | Modified MG | 36/5 | 38 | 0.9921 |  | 10.0 | 27127344.3 |
| 5% | polynomial | 8 | yes | Multigraph | 31/0 | 38 | 0.9748 |  | 10.0 | 27917652.5 |
| 5% | spectral | 12 | no | Modified MG | 35/3 | 39 | 0.9889 |  | 10.0 | 111486078.0 |
| 5% | spectral | 12 | no | Multigraph | 37/0 | 44 | 0.9991 |  | 0.0001 | 125136634.2 |
| 5% | spectral | 12 | yes | Modified MG | 34/0 | 41 | 1.0026 |  | 10.0 | 93454652.8 |
| 5% | spectral | 12 | yes | Multigraph | 37/0 | 44 | 0.9991 |  | 0.0001 | 125136813.2 |
| 5% | spectral | 4 | no | Modified MG | 4/7 | 4 | 1.0000 |  | 10.0 | 192109333.9 |
| 5% | spectral | 4 | no | Multigraph | 8/5 | 10 | 1.0000 |  | 10.0 | 64231689.3 |
| 5% | spectral | 4 | yes | Modified MG | 4/7 | 4 | 1.0000 |  | 10.0 | 118628576.9 |
| 5% | spectral | 4 | yes | Multigraph | 8/5 | 10 | 1.0000 |  | 10.0 | 63892845.2 |
| 5% | spectral | 8 | no | Modified MG | 33/5 | 35 | 0.9976 |  | 10.0 | 121753269.5 |
| 5% | spectral | 8 | no | Multigraph | 28/3 | 32 | 0.9961 |  | 5.0 | 93763895.2 |
| 5% | spectral | 8 | yes | Modified MG | 37/0 | 44 | 0.9993 |  | 0.0001 | 65761915.7 |
| 5% | spectral | 8 | yes | Multigraph | 28/3 | 32 | 0.9961 |  | 5.0 | 93755748.7 |

---
*Source: `mg_comparison_20260629_165010_full.csv`*
