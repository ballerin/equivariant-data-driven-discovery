# Porous Medium — Multigraph vs Modified Multigraph

Equation: $f_t = \Delta(f^3)/3$ (mod) = $f^2\Delta f + 2f\|\nabla f\|^2$ (mg)

| Library | True terms |
|---|---|
| Multigraph | f²·Δf, f·‖∇f‖² |
| Modified MG | Δf³ |

## Best per configuration (STRidge)

| Noise | Backend | m | IBP | MG fp/fn | MG E | Mod fp/fn | Mod E | Winner |
|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | 0/0 | 0.0007 | 0/0 | 0.0007 | 🟩 |
| clean | polynomial | 12 | yes | 0/0 | 0.0007 | 0/0 | 0.0007 | 🟩 |
| clean | polynomial | 4 | no | 0/0 | 0.0010 | 0/0 | 0.0010 | 🟩 |
| clean | polynomial | 4 | yes | 0/0 | 0.0010 | 0/0 | 0.0279 | 🟥 |
| clean | polynomial | 8 | no | 0/0 | 0.0009 | 0/0 | 0.0009 | 🟩 |
| clean | polynomial | 8 | yes | 0/0 | 0.0009 | 0/0 | 0.0009 | 🟩 |
| clean | spectral | 12 | no | 0/0 | 0.0008 | 0/0 | 0.0007 | 🟩 |
| clean | spectral | 12 | yes | 0/0 | 0.0008 | 0/0 | 0.0007 | 🟩 |
| clean | spectral | 4 | no | 0/0 | 0.0012 | 0/0 | 0.0011 | 🟩 |
| clean | spectral | 4 | yes | 0/0 | 0.0012 | 0/0 | 0.0279 | 🟥 |
| clean | spectral | 8 | no | 0/0 | 0.0010 | 0/0 | 0.0009 | 🟩 |
| clean | spectral | 8 | yes | 0/0 | 0.0010 | 0/0 | 0.0009 | 🟩 |
| 2% | polynomial | 12 | no | 6/2 | 1.5258 | 2/2 | 0.9266 | 🟢 |
| 2% | polynomial | 12 | yes | 6/2 | 1.5258 | 0/0 | 0.0014 | 🟢 |
| 2% | polynomial | 4 | no | 0/4 | 1.0000 | 0/4 | 1.0000 | ⚪ |
| 2% | polynomial | 4 | yes | 0/4 | 1.0000 | 0/4 | 1.0000 | ⚪ |
| 2% | polynomial | 8 | no | 0/0 | 0.0035 | 0/0 | 0.0031 | 🟩 |
| 2% | polynomial | 8 | yes | 0/0 | 0.0035 | 0/0 | 0.0028 | 🟩 |
| 2% | spectral | 12 | no | 0/0 | 0.0022 | 0/0 | 0.0004 | 🟩 |
| 2% | spectral | 12 | yes | 0/0 | 0.0022 | 0/0 | 0.0014 | 🟩 |
| 2% | spectral | 4 | no | 0/0 | 0.0152 | 0/4 | 1.0000 | 🔴 |
| 2% | spectral | 4 | yes | 0/4 | 1.0000 | 0/0 | 0.0109 | 🟢 |
| 2% | spectral | 8 | no | 0/4 | 1.0000 | 0/4 | 1.0000 | ⚪ |
| 2% | spectral | 8 | yes | 0/4 | 1.0000 | 0/0 | 0.0028 | 🟢 |
| 5% | polynomial | 12 | no | 0/4 | 1.0000 | 4/0 | 0.0430 | 🟩 |
| 5% | polynomial | 12 | yes | 0/4 | 1.0000 | 0/0 | 0.0022 | 🟢 |
| 5% | polynomial | 4 | no | 0/4 | 1.0000 | 0/0 | 0.0367 | 🟢 |
| 5% | polynomial | 4 | yes | 0/4 | 1.0000 | 0/4 | 1.0000 | ⚪ |
| 5% | polynomial | 8 | no | 0/4 | 1.0000 | 0/0 | 0.0070 | 🟢 |
| 5% | polynomial | 8 | yes | 0/4 | 1.0000 | 0/0 | 0.0047 | 🟢 |
| 5% | spectral | 12 | no | 0/0 | 0.0025 | 0/0 | 0.0046 | 🟥 |
| 5% | spectral | 12 | yes | 0/0 | 0.0025 | 0/0 | 0.0022 | 🟩 |
| 5% | spectral | 4 | no | 0/4 | 1.0000 | 0/4 | 1.0000 | ⚪ |
| 5% | spectral | 4 | yes | 0/4 | 1.0000 | 0/4 | 1.0000 | ⚪ |
| 5% | spectral | 8 | no | 0/4 | 1.0000 | 0/0 | 0.0006 | 🟢 |
| 5% | spectral | 8 | yes | 0/4 | 1.0000 | 0/0 | 0.0047 | 🟢 |

## Detailed results

| Noise | Backend | m | IBP | Library | fp/fn | n_act | E | Exact? | d_tol | cond |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | Modified MG | 0/0 | 4 | 0.0007 | ✓ | 0.1 | 1870424.8 |
| clean | polynomial | 12 | no | Multigraph | 0/0 | 4 | 0.0007 | ✓ | 1.0 | 1896428.8 |
| clean | polynomial | 12 | yes | Modified MG | 0/0 | 4 | 0.0007 | ✓ | 0.05 | 1981014.5 |
| clean | polynomial | 12 | yes | Multigraph | 0/0 | 4 | 0.0007 | ✓ | 1.0 | 1884770.8 |
| clean | polynomial | 4 | no | Modified MG | 0/0 | 4 | 0.0010 | ✓ | 0.05 | 393092.6 |
| clean | polynomial | 4 | no | Multigraph | 0/0 | 4 | 0.0010 | ✓ | 0.1 | 399891.4 |
| clean | polynomial | 4 | yes | Modified MG | 0/0 | 4 | 0.0279 | ✓ | 0.2 | 2631720.0 |
| clean | polynomial | 4 | yes | Multigraph | 0/0 | 4 | 0.0010 | ✓ | 0.05 | 1438261.1 |
| clean | polynomial | 8 | no | Modified MG | 0/0 | 4 | 0.0009 | ✓ | 0.05 | 1089047.0 |
| clean | polynomial | 8 | no | Multigraph | 0/0 | 4 | 0.0009 | ✓ | 0.05 | 1103259.0 |
| clean | polynomial | 8 | yes | Modified MG | 0/0 | 4 | 0.0009 | ✓ | 0.05 | 1448759.9 |
| clean | polynomial | 8 | yes | Multigraph | 0/0 | 4 | 0.0009 | ✓ | 0.1 | 1060615.2 |
| clean | spectral | 12 | no | Modified MG | 0/0 | 4 | 0.0007 | ✓ | 0.5 | 3029383.0 |
| clean | spectral | 12 | no | Multigraph | 0/0 | 4 | 0.0008 | ✓ | 0.1 | 3063426.0 |
| clean | spectral | 12 | yes | Modified MG | 0/0 | 4 | 0.0007 | ✓ | 0.5 | 3031653.1 |
| clean | spectral | 12 | yes | Multigraph | 0/0 | 4 | 0.0008 | ✓ | 0.1 | 3065203.8 |
| clean | spectral | 4 | no | Modified MG | 0/0 | 4 | 0.0011 | ✓ | 0.1 | 1119543.9 |
| clean | spectral | 4 | no | Multigraph | 0/0 | 4 | 0.0012 | ✓ | 0.005 | 1144526.1 |
| clean | spectral | 4 | yes | Modified MG | 0/0 | 4 | 0.0279 | ✓ | 0.1 | 3208058.9 |
| clean | spectral | 4 | yes | Multigraph | 0/0 | 4 | 0.0012 | ✓ | 0.1 | 3534221.1 |
| clean | spectral | 8 | no | Modified MG | 0/0 | 4 | 0.0009 | ✓ | 0.05 | 1908107.5 |
| clean | spectral | 8 | no | Multigraph | 0/0 | 4 | 0.0010 | ✓ | 0.2 | 1933914.7 |
| clean | spectral | 8 | yes | Modified MG | 0/0 | 4 | 0.0009 | ✓ | 0.05 | 1700172.7 |
| clean | spectral | 8 | yes | Multigraph | 0/0 | 4 | 0.0010 | ✓ | 0.2 | 1940044.0 |
| 2% | polynomial | 12 | no | Modified MG | 2/2 | 4 | 0.9266 |  | 5.0 | 179587.7 |
| 2% | polynomial | 12 | no | Multigraph | 6/2 | 8 | 1.5258 |  | 1.0 | 174841.5 |
| 2% | polynomial | 12 | yes | Modified MG | 0/0 | 4 | 0.0014 | ✓ | 0.5 | 201035.1 |
| 2% | polynomial | 12 | yes | Multigraph | 6/2 | 8 | 1.5258 |  | 1.0 | 174824.0 |
| 2% | polynomial | 4 | no | Modified MG | 0/4 | 0 | 1.0000 |  | 0.5 | 27361.5 |
| 2% | polynomial | 4 | no | Multigraph | 0/4 | 0 | 1.0000 |  | 2.0 | 26996.4 |
| 2% | polynomial | 4 | yes | Modified MG | 0/4 | 0 | 1.0000 |  | 1.0 | 94248.1 |
| 2% | polynomial | 4 | yes | Multigraph | 0/4 | 0 | 1.0000 |  | 1.0 | 88760.8 |
| 2% | polynomial | 8 | no | Modified MG | 0/0 | 4 | 0.0031 | ✓ | 0.5 | 74881.4 |
| 2% | polynomial | 8 | no | Multigraph | 0/0 | 4 | 0.0035 | ✓ | 0.5 | 74434.3 |
| 2% | polynomial | 8 | yes | Modified MG | 0/0 | 4 | 0.0028 | ✓ | 0.5 | 71326.3 |
| 2% | polynomial | 8 | yes | Multigraph | 0/0 | 4 | 0.0035 | ✓ | 0.5 | 72730.0 |
| 2% | spectral | 12 | no | Modified MG | 0/0 | 4 | 0.0004 | ✓ | 2.0 | 41406.5 |
| 2% | spectral | 12 | no | Multigraph | 0/0 | 4 | 0.0022 | ✓ | 0.5 | 25197.9 |
| 2% | spectral | 12 | yes | Modified MG | 0/0 | 4 | 0.0014 | ✓ | 1.0 | 27346.0 |
| 2% | spectral | 12 | yes | Multigraph | 0/0 | 4 | 0.0022 | ✓ | 0.5 | 25197.6 |
| 2% | spectral | 4 | no | Modified MG | 0/4 | 0 | 1.0000 |  | 0.5 | 5574.8 |
| 2% | spectral | 4 | no | Multigraph | 0/0 | 4 | 0.0152 | ✓ | 0.5 | 5461.2 |
| 2% | spectral | 4 | yes | Modified MG | 0/0 | 4 | 0.0109 | ✓ | 0.5 | 22060.2 |
| 2% | spectral | 4 | yes | Multigraph | 0/4 | 0 | 1.0000 |  | 0.5 | 18061.6 |
| 2% | spectral | 8 | no | Modified MG | 0/4 | 0 | 1.0000 |  | 2.0 | 16318.0 |
| 2% | spectral | 8 | no | Multigraph | 0/4 | 0 | 1.0000 |  | 2.0 | 14793.5 |
| 2% | spectral | 8 | yes | Modified MG | 0/0 | 4 | 0.0028 | ✓ | 0.2 | 15291.4 |
| 2% | spectral | 8 | yes | Multigraph | 0/4 | 0 | 1.0000 |  | 2.0 | 14794.1 |
| 5% | polynomial | 12 | no | Modified MG | 4/0 | 8 | 0.0430 |  | 1.0 | 62155.6 |
| 5% | polynomial | 12 | no | Multigraph | 0/4 | 0 | 1.0000 |  | 5.0 | 62737.0 |
| 5% | polynomial | 12 | yes | Modified MG | 0/0 | 4 | 0.0022 | ✓ | 1.0 | 69783.4 |
| 5% | polynomial | 12 | yes | Multigraph | 0/4 | 0 | 1.0000 |  | 5.0 | 62734.9 |
| 5% | polynomial | 4 | no | Modified MG | 0/0 | 4 | 0.0367 | ✓ | 0.2 | 9255.3 |
| 5% | polynomial | 4 | no | Multigraph | 0/4 | 0 | 1.0000 |  | 5.0 | 9053.8 |
| 5% | polynomial | 4 | yes | Modified MG | 0/4 | 0 | 1.0000 |  | 5.0 | 30506.3 |
| 5% | polynomial | 4 | yes | Multigraph | 0/4 | 0 | 1.0000 |  | 2.0 | 31368.5 |
| 5% | polynomial | 8 | no | Modified MG | 0/0 | 4 | 0.0070 | ✓ | 1.0 | 31271.7 |
| 5% | polynomial | 8 | no | Multigraph | 0/4 | 0 | 1.0000 |  | 5.0 | 30716.2 |
| 5% | polynomial | 8 | yes | Modified MG | 0/0 | 4 | 0.0047 | ✓ | 1.0 | 30453.9 |
| 5% | polynomial | 8 | yes | Multigraph | 0/4 | 0 | 1.0000 |  | 10.0 | 30454.4 |
| 5% | spectral | 12 | no | Modified MG | 0/0 | 4 | 0.0046 | ✓ | 1.0 | 14021.9 |
| 5% | spectral | 12 | no | Multigraph | 0/0 | 4 | 0.0025 | ✓ | 10.0 | 11056.3 |
| 5% | spectral | 12 | yes | Modified MG | 0/0 | 4 | 0.0022 | ✓ | 1.0 | 11869.6 |
| 5% | spectral | 12 | yes | Multigraph | 0/0 | 4 | 0.0025 | ✓ | 10.0 | 11050.3 |
| 5% | spectral | 4 | no | Modified MG | 0/4 | 0 | 1.0000 |  | 2.0 | 4894.6 |
| 5% | spectral | 4 | no | Multigraph | 0/4 | 0 | 1.0000 |  | 2.0 | 2218.6 |
| 5% | spectral | 4 | yes | Modified MG | 0/4 | 0 | 1.0000 |  | 1.0 | 7837.9 |
| 5% | spectral | 4 | yes | Multigraph | 0/4 | 0 | 1.0000 |  | 2.0 | 4890.3 |
| 5% | spectral | 8 | no | Modified MG | 0/0 | 4 | 0.0006 | ✓ | 1.0 | 6260.1 |
| 5% | spectral | 8 | no | Multigraph | 0/4 | 0 | 1.0000 |  | 5.0 | 5152.9 |
| 5% | spectral | 8 | yes | Modified MG | 0/0 | 4 | 0.0047 | ✓ | 0.5 | 6475.9 |
| 5% | spectral | 8 | yes | Multigraph | 0/4 | 0 | 1.0000 |  | 5.0 | 5133.2 |

---
*Source: `mg_comparison_20260629_165010_full.csv`*
