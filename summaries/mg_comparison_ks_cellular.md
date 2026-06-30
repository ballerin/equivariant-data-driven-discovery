# KS (cellular) — Multigraph vs Modified Multigraph

Equation: $f_t = -\Delta f - \Delta^2 f - \|\nabla f\|^2/2$

| Library | True terms |
|---|---|
| Multigraph | Δf, ‖∇f‖², Δ²f |
| Modified MG | Δf, ‖∇f‖², Δ²f |

## Best per configuration (STRidge)

| Noise | Backend | m | IBP | MG fp/fn | MG E | Mod fp/fn | Mod E | Winner |
|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | 0/0 | 0.0029 | 0/0 | 0.0029 | ⚪ |
| clean | polynomial | 12 | yes | 1/0 | 0.1388 | 1/0 | 0.1388 | ⚪ |
| clean | polynomial | 4 | no | 0/0 | 0.0038 | 0/0 | 0.0038 | ⚪ |
| clean | polynomial | 4 | yes | 15/2 | 4.8068 | 36/0 | 108.9949 | 🔴 |
| clean | polynomial | 8 | no | 0/0 | 0.0035 | 0/0 | 0.0035 | ⚪ |
| clean | polynomial | 8 | yes | 1/0 | 3.5865 | 1/0 | 3.5865 | ⚪ |
| clean | spectral | 12 | no | 0/0 | 0.0029 | 0/0 | 0.0029 | ⚪ |
| clean | spectral | 12 | yes | 1/0 | 0.1388 | 1/0 | 0.1388 | ⚪ |
| clean | spectral | 4 | no | 0/0 | 0.0038 | 0/0 | 0.0038 | ⚪ |
| clean | spectral | 4 | yes | 15/2 | 4.8068 | 36/0 | 108.9949 | 🔴 |
| clean | spectral | 8 | no | 0/0 | 0.0035 | 0/0 | 0.0035 | ⚪ |
| clean | spectral | 8 | yes | 1/0 | 3.5865 | 1/0 | 3.5865 | ⚪ |
| 2% | polynomial | 12 | no | 30/0 | 0.6454 | 35/0 | 0.6436 | 🔴 |
| 2% | polynomial | 12 | yes | 30/0 | 0.6509 | 37/0 | 0.7191 | 🔴 |
| 2% | polynomial | 4 | no | 3/5 | 1.0000 | 3/5 | 1.0000 | ⚪ |
| 2% | polynomial | 4 | yes | 3/5 | 1.0000 | 3/5 | 1.0000 | ⚪ |
| 2% | polynomial | 8 | no | 30/3 | 0.9843 | 32/3 | 0.9843 | 🔴 |
| 2% | polynomial | 8 | yes | 30/3 | 0.9840 | 26/3 | 0.9999 | 🟢 |
| 2% | spectral | 12 | no | 37/0 | 0.6510 | 37/0 | 0.6510 | 🟥 |
| 2% | spectral | 12 | yes | 37/0 | 0.6517 | 33/0 | 0.6878 | 🟢 |
| 2% | spectral | 4 | no | 7/5 | 1.0000 | 12/5 | 1.0000 | 🔴 |
| 2% | spectral | 4 | yes | 7/5 | 1.0000 | 5/7 | 1.0000 | 🟥 |
| 2% | spectral | 8 | no | 29/3 | 0.9868 | 35/3 | 0.9874 | 🔴 |
| 2% | spectral | 8 | yes | 29/3 | 0.9868 | 35/0 | 0.9842 | 🔴 |
| 5% | polynomial | 12 | no | 31/0 | 0.8713 | 22/0 | 0.8684 | 🟢 |
| 5% | polynomial | 12 | yes | 31/0 | 0.8733 | 35/0 | 0.8969 | 🔴 |
| 5% | polynomial | 4 | no | 7/5 | 1.0000 | 3/5 | 1.0000 | 🟢 |
| 5% | polynomial | 4 | yes | 7/5 | 1.0000 | 6/5 | 1.0000 | 🟢 |
| 5% | polynomial | 8 | no | 23/3 | 1.0001 | 19/3 | 1.0000 | 🟢 |
| 5% | polynomial | 8 | yes | 21/5 | 1.0000 | 27/5 | 1.0002 | 🔴 |
| 5% | spectral | 12 | no | 32/0 | 0.9201 | 37/0 | 0.9193 | 🔴 |
| 5% | spectral | 12 | yes | 32/0 | 0.9201 | 37/0 | 0.9352 | 🔴 |
| 5% | spectral | 4 | no | 13/5 | 1.0000 | 31/5 | 1.0001 | 🔴 |
| 5% | spectral | 4 | yes | 13/5 | 1.0000 | 9/7 | 1.0000 | 🟢 |
| 5% | spectral | 8 | no | 28/3 | 1.0003 | 35/0 | 0.9975 | 🔴 |
| 5% | spectral | 8 | yes | 28/3 | 1.0003 | 22/0 | 0.9968 | 🟢 |

## Detailed results

| Noise | Backend | m | IBP | Library | fp/fn | n_act | E | Exact? | d_tol | cond |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | Modified MG | 0/0 | 7 | 0.0029 | ✓ | 5.0 | 130678.0 |
| clean | polynomial | 12 | no | Multigraph | 0/0 | 7 | 0.0029 | ✓ | 5.0 | 125064.7 |
| clean | polynomial | 12 | yes | Modified MG | 1/0 | 8 | 0.1388 |  | 5.0 | 132114.0 |
| clean | polynomial | 12 | yes | Multigraph | 1/0 | 8 | 0.1388 |  | 5.0 | 126842.9 |
| clean | polynomial | 4 | no | Modified MG | 0/0 | 7 | 0.0038 | ✓ | 2.0 | 107674.9 |
| clean | polynomial | 4 | no | Multigraph | 0/0 | 7 | 0.0038 | ✓ | 2.0 | 105504.4 |
| clean | polynomial | 4 | yes | Modified MG | 36/0 | 43 | 108.9949 |  | 10.0 | 513374280.7 |
| clean | polynomial | 4 | yes | Multigraph | 15/2 | 20 | 4.8068 |  | 5.0 | 6541669.1 |
| clean | polynomial | 8 | no | Modified MG | 0/0 | 7 | 0.0035 | ✓ | 5.0 | 118284.0 |
| clean | polynomial | 8 | no | Multigraph | 0/0 | 7 | 0.0035 | ✓ | 5.0 | 114785.6 |
| clean | polynomial | 8 | yes | Modified MG | 1/0 | 8 | 3.5865 |  | 10.0 | 142215.5 |
| clean | polynomial | 8 | yes | Multigraph | 1/0 | 8 | 3.5865 |  | 5.0 | 133850.3 |
| clean | spectral | 12 | no | Modified MG | 0/0 | 7 | 0.0029 | ✓ | 5.0 | 130676.6 |
| clean | spectral | 12 | no | Multigraph | 0/0 | 7 | 0.0029 | ✓ | 5.0 | 125064.5 |
| clean | spectral | 12 | yes | Modified MG | 1/0 | 8 | 0.1388 |  | 5.0 | 132113.2 |
| clean | spectral | 12 | yes | Multigraph | 1/0 | 8 | 0.1388 |  | 5.0 | 126840.5 |
| clean | spectral | 4 | no | Modified MG | 0/0 | 7 | 0.0038 | ✓ | 2.0 | 107672.7 |
| clean | spectral | 4 | no | Multigraph | 0/0 | 7 | 0.0038 | ✓ | 2.0 | 105502.8 |
| clean | spectral | 4 | yes | Modified MG | 36/0 | 43 | 108.9949 |  | 10.0 | 513373848.7 |
| clean | spectral | 4 | yes | Multigraph | 15/2 | 20 | 4.8068 |  | 5.0 | 6541786.5 |
| clean | spectral | 8 | no | Modified MG | 0/0 | 7 | 0.0035 | ✓ | 5.0 | 118282.2 |
| clean | spectral | 8 | no | Multigraph | 0/0 | 7 | 0.0035 | ✓ | 5.0 | 114784.0 |
| clean | spectral | 8 | yes | Modified MG | 1/0 | 8 | 3.5865 |  | 10.0 | 142215.0 |
| clean | spectral | 8 | yes | Multigraph | 1/0 | 8 | 3.5865 |  | 5.0 | 133848.0 |
| 2% | polynomial | 12 | no | Modified MG | 35/0 | 42 | 0.6436 |  | 10.0 | 73531.1 |
| 2% | polynomial | 12 | no | Multigraph | 30/0 | 37 | 0.6454 |  | 10.0 | 75641.8 |
| 2% | polynomial | 12 | yes | Modified MG | 37/0 | 44 | 0.7191 |  | 5.0 | 71304.4 |
| 2% | polynomial | 12 | yes | Multigraph | 30/0 | 37 | 0.6509 |  | 10.0 | 75666.7 |
| 2% | polynomial | 4 | no | Modified MG | 3/5 | 5 | 1.0000 |  | 5.0 | 166584.3 |
| 2% | polynomial | 4 | no | Multigraph | 3/5 | 5 | 1.0000 |  | 10.0 | 55740.0 |
| 2% | polynomial | 4 | yes | Modified MG | 3/5 | 5 | 1.0000 |  | 5.0 | 2672571.6 |
| 2% | polynomial | 4 | yes | Multigraph | 3/5 | 5 | 1.0000 |  | 10.0 | 55767.2 |
| 2% | polynomial | 8 | no | Modified MG | 32/3 | 36 | 0.9843 |  | 10.0 | 53833.5 |
| 2% | polynomial | 8 | no | Multigraph | 30/3 | 34 | 0.9843 |  | 10.0 | 54983.2 |
| 2% | polynomial | 8 | yes | Modified MG | 26/3 | 30 | 0.9999 |  | 10.0 | 52184.4 |
| 2% | polynomial | 8 | yes | Multigraph | 30/3 | 34 | 0.9840 |  | 10.0 | 54914.0 |
| 2% | spectral | 12 | no | Modified MG | 37/0 | 44 | 0.6510 |  | 0.0001 | 645445.9 |
| 2% | spectral | 12 | no | Multigraph | 37/0 | 44 | 0.6510 |  | 0.0001 | 693750.9 |
| 2% | spectral | 12 | yes | Modified MG | 33/0 | 40 | 0.6878 |  | 10.0 | 444353.3 |
| 2% | spectral | 12 | yes | Multigraph | 37/0 | 44 | 0.6517 |  | 0.0001 | 693845.0 |
| 2% | spectral | 4 | no | Modified MG | 12/5 | 14 | 1.0000 |  | 5.0 | 1050694.1 |
| 2% | spectral | 4 | no | Multigraph | 7/5 | 9 | 1.0000 |  | 10.0 | 336798.9 |
| 2% | spectral | 4 | yes | Modified MG | 5/7 | 5 | 1.0000 |  | 10.0 | 3828793.7 |
| 2% | spectral | 4 | yes | Multigraph | 7/5 | 9 | 1.0000 |  | 10.0 | 338696.5 |
| 2% | spectral | 8 | no | Modified MG | 35/3 | 39 | 0.9874 |  | 10.0 | 633740.1 |
| 2% | spectral | 8 | no | Multigraph | 29/3 | 33 | 0.9868 |  | 10.0 | 471736.7 |
| 2% | spectral | 8 | yes | Modified MG | 35/0 | 42 | 0.9842 |  | 5.0 | 308500.2 |
| 2% | spectral | 8 | yes | Multigraph | 29/3 | 33 | 0.9868 |  | 10.0 | 472036.4 |
| 5% | polynomial | 12 | no | Modified MG | 22/0 | 29 | 0.8684 |  | 10.0 | 81617.7 |
| 5% | polynomial | 12 | no | Multigraph | 31/0 | 38 | 0.8713 |  | 5.0 | 85906.5 |
| 5% | polynomial | 12 | yes | Modified MG | 35/0 | 42 | 0.8969 |  | 5.0 | 72857.3 |
| 5% | polynomial | 12 | yes | Multigraph | 31/0 | 38 | 0.8733 |  | 5.0 | 85896.2 |
| 5% | polynomial | 4 | no | Modified MG | 3/5 | 5 | 1.0000 |  | 10.0 | 444605.1 |
| 5% | polynomial | 4 | no | Multigraph | 7/5 | 9 | 1.0000 |  | 10.0 | 149255.4 |
| 5% | polynomial | 4 | yes | Modified MG | 6/5 | 8 | 1.0000 |  | 10.0 | 3092163.6 |
| 5% | polynomial | 4 | yes | Multigraph | 7/5 | 9 | 1.0000 |  | 10.0 | 149255.1 |
| 5% | polynomial | 8 | no | Modified MG | 19/3 | 23 | 1.0000 |  | 10.0 | 60733.4 |
| 5% | polynomial | 8 | no | Multigraph | 23/3 | 27 | 1.0001 |  | 10.0 | 63243.2 |
| 5% | polynomial | 8 | yes | Modified MG | 27/5 | 29 | 1.0002 |  | 10.0 | 65941.4 |
| 5% | polynomial | 8 | yes | Multigraph | 21/5 | 23 | 1.0000 |  | 10.0 | 63144.3 |
| 5% | spectral | 12 | no | Modified MG | 37/0 | 44 | 0.9193 |  | 5.0 | 4537347.8 |
| 5% | spectral | 12 | no | Multigraph | 32/0 | 39 | 0.9201 |  | 10.0 | 4850904.7 |
| 5% | spectral | 12 | yes | Modified MG | 37/0 | 44 | 0.9352 |  | 5.0 | 3269545.0 |
| 5% | spectral | 12 | yes | Multigraph | 32/0 | 39 | 0.9201 |  | 10.0 | 4851117.8 |
| 5% | spectral | 4 | no | Modified MG | 31/5 | 33 | 1.0001 |  | 10.0 | 6175523.2 |
| 5% | spectral | 4 | no | Multigraph | 13/5 | 15 | 1.0000 |  | 5.0 | 2029412.3 |
| 5% | spectral | 4 | yes | Modified MG | 9/7 | 9 | 1.0000 |  | 10.0 | 4617059.9 |
| 5% | spectral | 4 | yes | Multigraph | 13/5 | 15 | 1.0000 |  | 5.0 | 2031277.0 |
| 5% | spectral | 8 | no | Modified MG | 35/0 | 42 | 0.9975 |  | 2.0 | 4638490.4 |
| 5% | spectral | 8 | no | Multigraph | 28/3 | 32 | 1.0003 |  | 5.0 | 3418856.0 |
| 5% | spectral | 8 | yes | Modified MG | 22/0 | 29 | 0.9968 |  | 10.0 | 2209880.3 |
| 5% | spectral | 8 | yes | Multigraph | 28/3 | 32 | 1.0003 |  | 10.0 | 3419390.6 |

---
*Source: `mg_comparison_20260629_165010_full.csv`*
