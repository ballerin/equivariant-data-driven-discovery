# Allen-Cahn — Multigraph vs Modified Multigraph

Equation: $f_t = \varepsilon^2 \Delta f + f - f^3$

| Library | True terms |
|---|---|
| Multigraph | f, f³, Δf |
| Modified MG | f, f³, Δf |

## Best per configuration (STRidge)

| Noise | Backend | m | IBP | MG fp/fn | MG E | Mod fp/fn | Mod E | Winner |
|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | 0/0 | 0.0048 | 0/0 | 0.0048 | ⚪ |
| clean | polynomial | 12 | yes | 0/0 | 0.0048 | 0/0 | 0.0048 | ⚪ |
| clean | polynomial | 4 | no | 0/0 | 0.0063 | 0/0 | 0.0063 | ⚪ |
| clean | polynomial | 4 | yes | 0/0 | 0.0157 | 0/0 | 0.0157 | ⚪ |
| clean | polynomial | 8 | no | 0/0 | 0.0052 | 0/0 | 0.0052 | ⚪ |
| clean | polynomial | 8 | yes | 0/0 | 0.0051 | 0/0 | 0.0051 | ⚪ |
| clean | spectral | 12 | no | 0/0 | 0.0048 | 0/0 | 0.0048 | ⚪ |
| clean | spectral | 12 | yes | 0/0 | 0.0048 | 0/0 | 0.0048 | ⚪ |
| clean | spectral | 4 | no | 0/0 | 0.0055 | 0/0 | 0.0055 | ⚪ |
| clean | spectral | 4 | yes | 0/0 | 0.0157 | 0/0 | 0.0157 | ⚪ |
| clean | spectral | 8 | no | 0/0 | 0.0051 | 0/0 | 0.0051 | ⚪ |
| clean | spectral | 8 | yes | 0/0 | 0.0051 | 0/0 | 0.0051 | ⚪ |
| 2% | polynomial | 12 | no | 0/0 | 0.0045 | 0/0 | 0.0045 | ⚪ |
| 2% | polynomial | 12 | yes | 0/0 | 0.0045 | 0/0 | 0.0045 | ⚪ |
| 2% | polynomial | 4 | no | 0/0 | 0.0098 | 0/0 | 0.0098 | ⚪ |
| 2% | polynomial | 4 | yes | 0/0 | 0.0283 | 0/0 | 0.0283 | ⚪ |
| 2% | polynomial | 8 | no | 0/0 | 0.0048 | 0/0 | 0.0048 | ⚪ |
| 2% | polynomial | 8 | yes | 0/0 | 0.0052 | 0/0 | 0.0052 | ⚪ |
| 2% | spectral | 12 | no | 0/0 | 0.0045 | 0/0 | 0.0045 | ⚪ |
| 2% | spectral | 12 | yes | 0/0 | 0.0045 | 0/0 | 0.0045 | ⚪ |
| 2% | spectral | 4 | no | 0/0 | 0.0200 | 0/0 | 0.0200 | ⚪ |
| 2% | spectral | 4 | yes | 0/0 | 0.0283 | 0/0 | 0.0283 | ⚪ |
| 2% | spectral | 8 | no | 0/0 | 0.0051 | 0/0 | 0.0051 | ⚪ |
| 2% | spectral | 8 | yes | 0/0 | 0.0052 | 0/0 | 0.0052 | ⚪ |
| 5% | polynomial | 12 | no | 0/0 | 0.0048 | 0/0 | 0.0048 | ⚪ |
| 5% | polynomial | 12 | yes | 0/0 | 0.0049 | 0/0 | 0.0049 | ⚪ |
| 5% | polynomial | 4 | no | 0/0 | 0.0506 | 0/0 | 0.0506 | ⚪ |
| 5% | polynomial | 4 | yes | 0/0 | 0.1170 | 0/0 | 0.1170 | ⚪ |
| 5% | polynomial | 8 | no | 0/0 | 0.0075 | 0/0 | 0.0075 | ⚪ |
| 5% | polynomial | 8 | yes | 0/0 | 0.0082 | 0/0 | 0.0082 | ⚪ |
| 5% | spectral | 12 | no | 0/0 | 0.0049 | 0/0 | 0.0049 | ⚪ |
| 5% | spectral | 12 | yes | 0/0 | 0.0049 | 0/0 | 0.0049 | ⚪ |
| 5% | spectral | 4 | no | 0/0 | 0.1105 | 0/0 | 0.1105 | ⚪ |
| 5% | spectral | 4 | yes | 0/0 | 0.1170 | 0/0 | 0.1170 | ⚪ |
| 5% | spectral | 8 | no | 0/0 | 0.0081 | 0/0 | 0.0081 | ⚪ |
| 5% | spectral | 8 | yes | 0/0 | 0.0082 | 0/0 | 0.0082 | ⚪ |

## Detailed results

| Noise | Backend | m | IBP | Library | fp/fn | n_act | E | Exact? | d_tol | cond |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | Modified MG | 0/0 | 4 | 0.0048 | ✓ | 0.05 | 1171.3 |
| clean | polynomial | 12 | no | Multigraph | 0/0 | 4 | 0.0048 | ✓ | 0.2 | 2568.6 |
| clean | polynomial | 12 | yes | Modified MG | 0/0 | 4 | 0.0048 | ✓ | 0.1 | 1163.8 |
| clean | polynomial | 12 | yes | Multigraph | 0/0 | 4 | 0.0048 | ✓ | 0.05 | 2563.8 |
| clean | polynomial | 4 | no | Modified MG | 0/0 | 4 | 0.0063 | ✓ | 0.05 | 3985.4 |
| clean | polynomial | 4 | no | Multigraph | 0/0 | 4 | 0.0063 | ✓ | 0.05 | 1704.3 |
| clean | polynomial | 4 | yes | Modified MG | 0/0 | 4 | 0.0157 | ✓ | 0.05 | 37566.9 |
| clean | polynomial | 4 | yes | Multigraph | 0/0 | 4 | 0.0157 | ✓ | 0.05 | 25885.5 |
| clean | polynomial | 8 | no | Modified MG | 0/0 | 4 | 0.0052 | ✓ | 0.1 | 1501.1 |
| clean | polynomial | 8 | no | Multigraph | 0/0 | 4 | 0.0052 | ✓ | 0.1 | 1418.8 |
| clean | polynomial | 8 | yes | Modified MG | 0/0 | 4 | 0.0051 | ✓ | 0.05 | 1505.1 |
| clean | polynomial | 8 | yes | Multigraph | 0/0 | 4 | 0.0051 | ✓ | 0.2 | 1406.6 |
| clean | spectral | 12 | no | Modified MG | 0/0 | 4 | 0.0048 | ✓ | 0.05 | 1256.9 |
| clean | spectral | 12 | no | Multigraph | 0/0 | 4 | 0.0048 | ✓ | 0.2 | 2683.0 |
| clean | spectral | 12 | yes | Modified MG | 0/0 | 4 | 0.0048 | ✓ | 0.05 | 1258.3 |
| clean | spectral | 12 | yes | Multigraph | 0/0 | 4 | 0.0048 | ✓ | 0.2 | 2682.2 |
| clean | spectral | 4 | no | Modified MG | 0/0 | 4 | 0.0055 | ✓ | 0.05 | 3847.1 |
| clean | spectral | 4 | no | Multigraph | 0/0 | 4 | 0.0055 | ✓ | 0.01 | 1742.0 |
| clean | spectral | 4 | yes | Modified MG | 0/0 | 4 | 0.0157 | ✓ | 0.05 | 36624.2 |
| clean | spectral | 4 | yes | Multigraph | 0/0 | 4 | 0.0157 | ✓ | 0.05 | 25684.2 |
| clean | spectral | 8 | no | Modified MG | 0/0 | 4 | 0.0051 | ✓ | 0.05 | 1485.1 |
| clean | spectral | 8 | no | Multigraph | 0/0 | 4 | 0.0051 | ✓ | 0.5 | 1475.2 |
| clean | spectral | 8 | yes | Modified MG | 0/0 | 4 | 0.0051 | ✓ | 0.05 | 1508.9 |
| clean | spectral | 8 | yes | Multigraph | 0/0 | 4 | 0.0051 | ✓ | 0.5 | 1477.8 |
| 2% | polynomial | 12 | no | Modified MG | 0/0 | 4 | 0.0045 | ✓ | 0.05 | 1160.8 |
| 2% | polynomial | 12 | no | Multigraph | 0/0 | 4 | 0.0045 | ✓ | 0.2 | 2582.2 |
| 2% | polynomial | 12 | yes | Modified MG | 0/0 | 4 | 0.0045 | ✓ | 0.5 | 1154.8 |
| 2% | polynomial | 12 | yes | Multigraph | 0/0 | 4 | 0.0045 | ✓ | 0.2 | 2577.5 |
| 2% | polynomial | 4 | no | Modified MG | 0/0 | 4 | 0.0098 | ✓ | 0.1 | 3890.9 |
| 2% | polynomial | 4 | no | Multigraph | 0/0 | 4 | 0.0098 | ✓ | 0.1 | 1731.3 |
| 2% | polynomial | 4 | yes | Modified MG | 0/0 | 4 | 0.0283 | ✓ | 0.5 | 19196.1 |
| 2% | polynomial | 4 | yes | Multigraph | 0/0 | 4 | 0.0283 | ✓ | 0.5 | 11217.1 |
| 2% | polynomial | 8 | no | Modified MG | 0/0 | 4 | 0.0048 | ✓ | 0.05 | 1491.0 |
| 2% | polynomial | 8 | no | Multigraph | 0/0 | 4 | 0.0048 | ✓ | 0.05 | 1427.5 |
| 2% | polynomial | 8 | yes | Modified MG | 0/0 | 4 | 0.0052 | ✓ | 0.5 | 1490.3 |
| 2% | polynomial | 8 | yes | Multigraph | 0/0 | 4 | 0.0052 | ✓ | 0.05 | 1414.7 |
| 2% | spectral | 12 | no | Modified MG | 0/0 | 4 | 0.0045 | ✓ | 0.1 | 2304.7 |
| 2% | spectral | 12 | no | Multigraph | 0/0 | 4 | 0.0045 | ✓ | 0.1 | 4037.5 |
| 2% | spectral | 12 | yes | Modified MG | 0/0 | 4 | 0.0045 | ✓ | 0.1 | 2733.1 |
| 2% | spectral | 12 | yes | Multigraph | 0/0 | 4 | 0.0045 | ✓ | 0.1 | 4037.9 |
| 2% | spectral | 4 | no | Modified MG | 0/0 | 4 | 0.0200 | ✓ | 0.2 | 4861.0 |
| 2% | spectral | 4 | no | Multigraph | 0/0 | 4 | 0.0200 | ✓ | 0.2 | 2330.4 |
| 2% | spectral | 4 | yes | Modified MG | 0/0 | 4 | 0.0283 | ✓ | 0.5 | 18138.9 |
| 2% | spectral | 4 | yes | Multigraph | 0/0 | 4 | 0.0283 | ✓ | 0.5 | 11816.9 |
| 2% | spectral | 8 | no | Modified MG | 0/0 | 4 | 0.0051 | ✓ | 0.1 | 2056.2 |
| 2% | spectral | 8 | no | Multigraph | 0/0 | 4 | 0.0051 | ✓ | 0.1 | 2821.4 |
| 2% | spectral | 8 | yes | Modified MG | 0/0 | 4 | 0.0052 | ✓ | 0.1 | 1934.6 |
| 2% | spectral | 8 | yes | Multigraph | 0/0 | 4 | 0.0052 | ✓ | 0.1 | 2841.3 |
| 5% | polynomial | 12 | no | Modified MG | 0/0 | 4 | 0.0048 | ✓ | 1.0 | 1134.5 |
| 5% | polynomial | 12 | no | Multigraph | 0/0 | 4 | 0.0048 | ✓ | 0.5 | 2531.0 |
| 5% | polynomial | 12 | yes | Modified MG | 0/0 | 4 | 0.0049 | ✓ | 0.5 | 1161.2 |
| 5% | polynomial | 12 | yes | Multigraph | 0/0 | 4 | 0.0049 | ✓ | 0.5 | 2527.3 |
| 5% | polynomial | 4 | no | Modified MG | 0/0 | 4 | 0.0506 | ✓ | 0.2 | 3941.4 |
| 5% | polynomial | 4 | no | Multigraph | 0/0 | 4 | 0.0506 | ✓ | 0.2 | 2046.5 |
| 5% | polynomial | 4 | yes | Modified MG | 0/0 | 4 | 0.1170 | ✓ | 1.0 | 14240.4 |
| 5% | polynomial | 4 | yes | Multigraph | 0/0 | 4 | 0.1170 | ✓ | 1.0 | 9239.9 |
| 5% | polynomial | 8 | no | Modified MG | 0/0 | 4 | 0.0075 | ✓ | 0.2 | 1455.3 |
| 5% | polynomial | 8 | no | Multigraph | 0/0 | 4 | 0.0075 | ✓ | 0.2 | 1462.6 |
| 5% | polynomial | 8 | yes | Modified MG | 0/0 | 4 | 0.0082 | ✓ | 1.0 | 1444.7 |
| 5% | polynomial | 8 | yes | Multigraph | 0/0 | 4 | 0.0082 | ✓ | 0.2 | 1449.7 |
| 5% | spectral | 12 | no | Modified MG | 0/0 | 4 | 0.0049 | ✓ | 0.2 | 11483.3 |
| 5% | spectral | 12 | no | Multigraph | 0/0 | 4 | 0.0049 | ✓ | 2.0 | 19563.0 |
| 5% | spectral | 12 | yes | Modified MG | 0/0 | 4 | 0.0049 | ✓ | 0.5 | 13027.6 |
| 5% | spectral | 12 | yes | Multigraph | 0/0 | 4 | 0.0049 | ✓ | 2.0 | 19566.7 |
| 5% | spectral | 4 | no | Modified MG | 0/0 | 4 | 0.1105 | ✓ | 0.5 | 19374.6 |
| 5% | spectral | 4 | no | Multigraph | 0/0 | 4 | 0.1105 | ✓ | 0.5 | 10382.8 |
| 5% | spectral | 4 | yes | Modified MG | 0/0 | 4 | 0.1170 | ✓ | 1.0 | 14721.9 |
| 5% | spectral | 4 | yes | Multigraph | 0/0 | 4 | 0.1170 | ✓ | 1.0 | 12502.2 |
| 5% | spectral | 8 | no | Modified MG | 0/0 | 4 | 0.0081 | ✓ | 0.2 | 10706.3 |
| 5% | spectral | 8 | no | Multigraph | 0/0 | 4 | 0.0081 | ✓ | 0.5 | 13252.1 |
| 5% | spectral | 8 | yes | Modified MG | 0/0 | 4 | 0.0082 | ✓ | 0.2 | 8908.1 |
| 5% | spectral | 8 | yes | Multigraph | 0/0 | 4 | 0.0082 | ✓ | 0.5 | 13321.0 |

---
*Source: `mg_comparison_20260629_165010_full.csv`*
