# Cahn-Hilliard — Multigraph vs Modified Multigraph

Equation: $f_t = -\Delta^2 f + \Delta(f^3 - f)$

| Library | True terms |
|---|---|
| Multigraph | Δf, f²·Δf, f·‖∇f‖², Δ²f |
| Modified MG | Δf, Δf³, Δ²f |

## Best per configuration (STRidge)

| Noise | Backend | m | IBP | MG fp/fn | MG E | Mod fp/fn | Mod E | Winner |
|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | 5/0 | 1.2376 | 18/0 | 2.2088 | 🔴 |
| clean | polynomial | 12 | yes | 5/0 | 1.1910 | 6/0 | 0.4364 | 🔴 |
| clean | polynomial | 4 | no | 7/0 | 1.3089 | 7/0 | 1.3093 | 🟥 |
| clean | polynomial | 4 | yes | 7/0 | 14.8645 | 14/0 | 21.4485 | 🔴 |
| clean | polynomial | 8 | no | 15/0 | 3.8627 | 16/0 | 3.9438 | 🔴 |
| clean | polynomial | 8 | yes | 8/0 | 1.4129 | 14/0 | 0.6548 | 🔴 |
| clean | spectral | 12 | no | 3/0 | 0.5489 | 6/0 | 0.3653 | 🔴 |
| clean | spectral | 12 | yes | 3/0 | 0.6183 | 6/0 | 0.4291 | 🔴 |
| clean | spectral | 4 | no | 3/0 | 0.1363 | 0/0 | 0.0834 | 🟢 |
| clean | spectral | 4 | yes | 13/0 | 8.8112 | 18/0 | 24.3715 | 🔴 |
| clean | spectral | 8 | no | 0/0 | 0.0577 | 0/0 | 0.0559 | 🟩 |
| clean | spectral | 8 | yes | 1/0 | 0.5784 | 14/0 | 0.6600 | 🔴 |
| 0.1% | polynomial | 12 | no | 5/0 | 0.1063 | 18/0 | 0.1269 | 🔴 |
| 0.1% | polynomial | 12 | yes | 5/0 | 0.1037 | 8/0 | 0.0639 | 🔴 |
| 0.1% | polynomial | 4 | no | 10/0 | 0.9588 | 10/0 | 0.9587 | 🟩 |
| 0.1% | polynomial | 4 | yes | 10/0 | 0.9351 | 17/0 | 0.9519 | 🔴 |
| 0.1% | polynomial | 8 | no | 15/0 | 0.2803 | 16/0 | 0.2727 | 🔴 |
| 0.1% | polynomial | 8 | yes | 8/0 | 0.1380 | 14/0 | 0.0709 | 🔴 |
| 0.1% | spectral | 12 | no | 35/0 | 0.0740 | 8/0 | 0.0662 | 🟢 |
| 0.1% | spectral | 12 | yes | 35/0 | 0.0681 | 0/0 | 0.0401 | 🟢 |
| 0.1% | spectral | 4 | no | 9/0 | 0.8855 | 18/0 | 0.8947 | 🔴 |
| 0.1% | spectral | 4 | yes | 8/0 | 0.9690 | 14/0 | 0.9632 | 🔴 |
| 0.1% | spectral | 8 | no | 0/0 | 0.0646 | 0/0 | 0.0616 | 🟩 |
| 0.1% | spectral | 8 | yes | 1/0 | 0.0856 | 14/0 | 0.0734 | 🔴 |
| 0.5% | polynomial | 12 | no | 12/0 | 0.1621 | 18/0 | 0.1627 | 🔴 |
| 0.5% | polynomial | 12 | yes | 9/0 | 0.1454 | 6/0 | 0.0601 | 🟢 |
| 0.5% | polynomial | 4 | no | 10/3 | 0.9977 | 9/7 | 0.9986 | 🔴 |
| 0.5% | polynomial | 4 | yes | 10/0 | 0.9886 | 17/0 | 0.9939 | 🔴 |
| 0.5% | polynomial | 8 | no | 8/0 | 0.4578 | 16/0 | 0.4471 | 🔴 |
| 0.5% | polynomial | 8 | yes | 11/0 | 0.3807 | 18/0 | 0.3790 | 🔴 |
| 0.5% | spectral | 12 | no | 2/0 | 0.0748 | 18/0 | 0.1085 | 🔴 |
| 0.5% | spectral | 12 | yes | 0/0 | 0.0545 | 0/0 | 0.0493 | 🟩 |
| 0.5% | spectral | 4 | no | 9/5 | 0.9982 | 9/7 | 0.9989 | 🔴 |
| 0.5% | spectral | 4 | yes | 0/9 | 1.0000 | 5/9 | 0.9995 | 🔴 |
| 0.5% | spectral | 8 | no | 11/0 | 0.3648 | 18/0 | 0.3888 | 🔴 |
| 0.5% | spectral | 8 | yes | 12/0 | 0.3659 | 18/0 | 0.3713 | 🔴 |
| 1% | polynomial | 12 | no | 12/0 | 0.2328 | 18/0 | 0.2340 | 🔴 |
| 1% | polynomial | 12 | yes | 13/0 | 0.2311 | 18/0 | 0.1963 | 🔴 |
| 1% | polynomial | 4 | no | 9/7 | 0.9986 | 9/7 | 0.9986 | ⚪ |
| 1% | polynomial | 4 | yes | 9/3 | 0.9984 | 9/7 | 0.9986 | 🔴 |
| 1% | polynomial | 8 | no | 13/0 | 0.6684 | 18/0 | 0.6672 | 🔴 |
| 1% | polynomial | 8 | yes | 11/0 | 0.6512 | 17/0 | 0.6532 | 🔴 |
| 1% | spectral | 12 | no | 17/0 | 0.2047 | 18/0 | 0.1946 | 🔴 |
| 1% | spectral | 12 | yes | 17/0 | 0.2019 | 5/0 | 0.1524 | 🟢 |
| 1% | spectral | 4 | no | 6/2 | 1.0042 | 0/9 | 1.0000 | 🔴 |
| 1% | spectral | 4 | yes | 0/9 | 1.0000 | 0/9 | 1.0000 | ⚪ |
| 1% | spectral | 8 | no | 14/0 | 0.6549 | 16/0 | 0.6487 | 🔴 |
| 1% | spectral | 8 | yes | 14/0 | 0.6552 | 18/0 | 0.6475 | 🔴 |
| 2% | polynomial | 12 | no | 12/0 | 0.3843 | 18/0 | 0.3860 | 🔴 |
| 2% | polynomial | 12 | yes | 14/0 | 0.3795 | 18/0 | 0.3659 | 🔴 |
| 2% | polynomial | 4 | no | 9/3 | 0.9985 | 9/7 | 0.9986 | 🔴 |
| 2% | polynomial | 4 | yes | 9/3 | 0.9985 | 0/9 | 1.0000 | 🟢 |
| 2% | polynomial | 8 | no | 14/0 | 0.8656 | 18/0 | 0.8652 | 🔴 |
| 2% | polynomial | 8 | yes | 13/0 | 0.8755 | 18/0 | 0.8711 | 🔴 |
| 2% | spectral | 12 | no | 18/0 | 0.3526 | 18/0 | 0.3417 | 🟩 |
| 2% | spectral | 12 | yes | 18/0 | 0.3508 | 7/0 | 0.2787 | 🟢 |
| 2% | spectral | 4 | no | 6/2 | 1.0051 | 6/4 | 1.0053 | 🔴 |
| 2% | spectral | 4 | yes | 0/9 | 1.0000 | 2/4 | 1.0035 | 🟢 |
| 2% | spectral | 8 | no | 17/0 | 0.8556 | 18/0 | 0.8563 | 🔴 |
| 2% | spectral | 8 | yes | 17/0 | 0.8561 | 7/0 | 0.8613 | 🟢 |
| 5% | polynomial | 12 | no | 14/0 | 0.7374 | 14/0 | 0.7305 | 🟩 |
| 5% | polynomial | 12 | yes | 14/0 | 0.7385 | 18/0 | 0.7132 | 🔴 |
| 5% | polynomial | 4 | no | 6/4 | 1.0090 | 0/9 | 1.0000 | 🟢 |
| 5% | polynomial | 4 | yes | 0/9 | 1.0000 | 0/9 | 1.0000 | ⚪ |
| 5% | polynomial | 8 | no | 11/3 | 0.9970 | 11/3 | 0.9970 | 🟩 |
| 5% | polynomial | 8 | yes | 11/3 | 0.9970 | 11/3 | 0.9970 | 🟩 |
| 5% | spectral | 12 | no | 35/0 | 0.7659 | 35/0 | 0.7659 | 🟥 |
| 5% | spectral | 12 | yes | 35/0 | 0.7656 | 18/0 | 0.7025 | 🟢 |
| 5% | spectral | 4 | no | 0/9 | 1.0000 | 0/9 | 1.0000 | ⚪ |
| 5% | spectral | 4 | yes | 0/9 | 1.0000 | 0/9 | 1.0000 | ⚪ |
| 5% | spectral | 8 | no | 35/0 | 0.9942 | 35/0 | 0.9942 | 🟩 |
| 5% | spectral | 8 | yes | 35/0 | 0.9944 | 2/7 | 0.9983 | 🟢 |

## Detailed results

| Noise | Backend | m | IBP | Library | fp/fn | n_act | E | Exact? | d_tol | cond |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | Modified MG | 18/0 | 27 | 2.2088 |  | 1.0 | 4588.1 |
| clean | polynomial | 12 | no | Multigraph | 5/0 | 14 | 1.2376 |  | 10.0 | 6897.0 |
| clean | polynomial | 12 | yes | Modified MG | 6/0 | 15 | 0.4364 |  | 10.0 | 4594.3 |
| clean | polynomial | 12 | yes | Multigraph | 5/0 | 14 | 1.1910 |  | 10.0 | 6897.3 |
| clean | polynomial | 4 | no | Modified MG | 7/0 | 16 | 1.3093 |  | 10.0 | 25129.7 |
| clean | polynomial | 4 | no | Multigraph | 7/0 | 16 | 1.3089 |  | 10.0 | 8893.6 |
| clean | polynomial | 4 | yes | Modified MG | 14/0 | 23 | 21.4485 |  | 10.0 | 176565.3 |
| clean | polynomial | 4 | yes | Multigraph | 7/0 | 16 | 14.8645 |  | 10.0 | 151189.3 |
| clean | polynomial | 8 | no | Modified MG | 16/0 | 25 | 3.9438 |  | 5.0 | 8993.5 |
| clean | polynomial | 8 | no | Multigraph | 15/0 | 24 | 3.8627 |  | 10.0 | 6603.5 |
| clean | polynomial | 8 | yes | Modified MG | 14/0 | 23 | 0.6548 |  | 10.0 | 9241.2 |
| clean | polynomial | 8 | yes | Multigraph | 8/0 | 17 | 1.4129 |  | 10.0 | 6798.2 |
| clean | spectral | 12 | no | Modified MG | 6/0 | 15 | 0.3653 |  | 10.0 | 5058.5 |
| clean | spectral | 12 | no | Multigraph | 3/0 | 12 | 0.5489 |  | 10.0 | 7501.5 |
| clean | spectral | 12 | yes | Modified MG | 6/0 | 15 | 0.4291 |  | 5.0 | 5065.0 |
| clean | spectral | 12 | yes | Multigraph | 3/0 | 12 | 0.6183 |  | 10.0 | 7505.1 |
| clean | spectral | 4 | no | Modified MG | 0/0 | 9 | 0.0834 | ✓ | 10.0 | 30569.3 |
| clean | spectral | 4 | no | Multigraph | 3/0 | 12 | 0.1363 |  | 2.0 | 11574.8 |
| clean | spectral | 4 | yes | Modified MG | 18/0 | 27 | 24.3715 |  | 1.0 | 188308.1 |
| clean | spectral | 4 | yes | Multigraph | 13/0 | 22 | 8.8112 |  | 5.0 | 115242.0 |
| clean | spectral | 8 | no | Modified MG | 0/0 | 9 | 0.0559 | ✓ | 10.0 | 10226.5 |
| clean | spectral | 8 | no | Multigraph | 0/0 | 9 | 0.0577 | ✓ | 5.0 | 7441.2 |
| clean | spectral | 8 | yes | Modified MG | 14/0 | 23 | 0.6600 |  | 5.0 | 10328.0 |
| clean | spectral | 8 | yes | Multigraph | 1/0 | 10 | 0.5784 |  | 10.0 | 7494.7 |
| 0.1% | polynomial | 12 | no | Modified MG | 18/0 | 27 | 0.1269 |  | 10.0 | 4586.3 |
| 0.1% | polynomial | 12 | no | Multigraph | 5/0 | 14 | 0.1063 |  | 10.0 | 6895.9 |
| 0.1% | polynomial | 12 | yes | Modified MG | 8/0 | 17 | 0.0639 |  | 5.0 | 4592.0 |
| 0.1% | polynomial | 12 | yes | Multigraph | 5/0 | 14 | 0.1037 |  | 10.0 | 6896.1 |
| 0.1% | polynomial | 4 | no | Modified MG | 10/0 | 19 | 0.9587 |  | 5.0 | 24989.7 |
| 0.1% | polynomial | 4 | no | Multigraph | 10/0 | 19 | 0.9588 |  | 5.0 | 8849.8 |
| 0.1% | polynomial | 4 | yes | Modified MG | 17/0 | 26 | 0.9519 |  | 10.0 | 163134.9 |
| 0.1% | polynomial | 4 | yes | Multigraph | 10/0 | 19 | 0.9351 |  | 5.0 | 124946.6 |
| 0.1% | polynomial | 8 | no | Modified MG | 16/0 | 25 | 0.2727 |  | 5.0 | 8984.9 |
| 0.1% | polynomial | 8 | no | Multigraph | 15/0 | 24 | 0.2803 |  | 10.0 | 6596.5 |
| 0.1% | polynomial | 8 | yes | Modified MG | 14/0 | 23 | 0.0709 |  | 10.0 | 9226.3 |
| 0.1% | polynomial | 8 | yes | Multigraph | 8/0 | 17 | 0.1380 |  | 10.0 | 6787.4 |
| 0.1% | spectral | 12 | no | Modified MG | 8/0 | 17 | 0.0662 |  | 5.0 | 5068.3 |
| 0.1% | spectral | 12 | no | Multigraph | 35/0 | 44 | 0.0740 |  | 0.0001 | 7540.6 |
| 0.1% | spectral | 12 | yes | Modified MG | 0/0 | 9 | 0.0401 | ✓ | 10.0 | 5067.4 |
| 0.1% | spectral | 12 | yes | Multigraph | 35/0 | 44 | 0.0681 |  | 0.0001 | 7544.3 |
| 0.1% | spectral | 4 | no | Modified MG | 18/0 | 27 | 0.8947 |  | 1.0 | 25888.5 |
| 0.1% | spectral | 4 | no | Multigraph | 9/0 | 18 | 0.8855 |  | 10.0 | 9870.8 |
| 0.1% | spectral | 4 | yes | Modified MG | 14/0 | 23 | 0.9632 |  | 10.0 | 126574.2 |
| 0.1% | spectral | 4 | yes | Multigraph | 8/0 | 17 | 0.9690 |  | 10.0 | 107861.7 |
| 0.1% | spectral | 8 | no | Modified MG | 0/0 | 9 | 0.0616 | ✓ | 10.0 | 10133.9 |
| 0.1% | spectral | 8 | no | Multigraph | 0/0 | 9 | 0.0646 | ✓ | 10.0 | 7426.0 |
| 0.1% | spectral | 8 | yes | Modified MG | 14/0 | 23 | 0.0734 |  | 5.0 | 10224.7 |
| 0.1% | spectral | 8 | yes | Multigraph | 1/0 | 10 | 0.0856 |  | 10.0 | 7482.1 |
| 0.5% | polynomial | 12 | no | Modified MG | 18/0 | 27 | 0.1627 |  | 2.0 | 4540.9 |
| 0.5% | polynomial | 12 | no | Multigraph | 12/0 | 21 | 0.1621 |  | 5.0 | 6885.3 |
| 0.5% | polynomial | 12 | yes | Modified MG | 6/0 | 15 | 0.0601 |  | 10.0 | 4532.8 |
| 0.5% | polynomial | 12 | yes | Multigraph | 9/0 | 18 | 0.1454 |  | 10.0 | 6884.4 |
| 0.5% | polynomial | 4 | no | Modified MG | 9/7 | 11 | 0.9986 |  | 10.0 | 23876.6 |
| 0.5% | polynomial | 4 | no | Multigraph | 10/3 | 16 | 0.9977 |  | 5.0 | 8503.5 |
| 0.5% | polynomial | 4 | yes | Modified MG | 17/0 | 26 | 0.9939 |  | 5.0 | 114143.3 |
| 0.5% | polynomial | 4 | yes | Multigraph | 10/0 | 19 | 0.9886 |  | 5.0 | 62209.9 |
| 0.5% | polynomial | 8 | no | Modified MG | 16/0 | 25 | 0.4471 |  | 2.0 | 8869.1 |
| 0.5% | polynomial | 8 | no | Multigraph | 8/0 | 17 | 0.4578 |  | 10.0 | 6497.2 |
| 0.5% | polynomial | 8 | yes | Modified MG | 18/0 | 27 | 0.3790 |  | 1.0 | 9015.3 |
| 0.5% | polynomial | 8 | yes | Multigraph | 11/0 | 20 | 0.3807 |  | 10.0 | 6633.4 |
| 0.5% | spectral | 12 | no | Modified MG | 18/0 | 27 | 0.1085 |  | 10.0 | 5751.2 |
| 0.5% | spectral | 12 | no | Multigraph | 2/0 | 11 | 0.0748 |  | 10.0 | 8797.1 |
| 0.5% | spectral | 12 | yes | Modified MG | 0/0 | 9 | 0.0493 | ✓ | 10.0 | 4843.6 |
| 0.5% | spectral | 12 | yes | Multigraph | 0/0 | 9 | 0.0545 | ✓ | 10.0 | 8802.7 |
| 0.5% | spectral | 4 | no | Modified MG | 9/7 | 11 | 0.9989 |  | 5.0 | 21717.2 |
| 0.5% | spectral | 4 | no | Multigraph | 9/5 | 13 | 0.9982 |  | 10.0 | 8277.0 |
| 0.5% | spectral | 4 | yes | Modified MG | 5/9 | 5 | 0.9995 |  | 10.0 | 82339.1 |
| 0.5% | spectral | 4 | yes | Multigraph | 0/9 | 0 | 1.0000 |  | 10.0 | 65284.0 |
| 0.5% | spectral | 8 | no | Modified MG | 18/0 | 27 | 0.3888 |  | 0.5 | 8997.4 |
| 0.5% | spectral | 8 | no | Multigraph | 11/0 | 20 | 0.3648 |  | 10.0 | 7785.3 |
| 0.5% | spectral | 8 | yes | Modified MG | 18/0 | 27 | 0.3713 |  | 2.0 | 7829.8 |
| 0.5% | spectral | 8 | yes | Multigraph | 12/0 | 21 | 0.3659 |  | 10.0 | 7825.5 |
| 1% | polynomial | 12 | no | Modified MG | 18/0 | 27 | 0.2340 |  | 5.0 | 4460.5 |
| 1% | polynomial | 12 | no | Multigraph | 12/0 | 21 | 0.2328 |  | 10.0 | 6870.8 |
| 1% | polynomial | 12 | yes | Modified MG | 18/0 | 27 | 0.1963 |  | 5.0 | 4431.8 |
| 1% | polynomial | 12 | yes | Multigraph | 13/0 | 22 | 0.2311 |  | 10.0 | 6868.4 |
| 1% | polynomial | 4 | no | Modified MG | 9/7 | 11 | 0.9986 |  | 10.0 | 23067.4 |
| 1% | polynomial | 4 | no | Multigraph | 9/7 | 11 | 0.9986 |  | 10.0 | 8285.2 |
| 1% | polynomial | 4 | yes | Modified MG | 9/7 | 11 | 0.9986 |  | 10.0 | 86872.8 |
| 1% | polynomial | 4 | yes | Multigraph | 9/3 | 15 | 0.9984 |  | 10.0 | 53167.9 |
| 1% | polynomial | 8 | no | Modified MG | 18/0 | 27 | 0.6672 |  | 1.0 | 8597.5 |
| 1% | polynomial | 8 | no | Multigraph | 13/0 | 22 | 0.6684 |  | 1.0 | 6293.6 |
| 1% | polynomial | 8 | yes | Modified MG | 17/0 | 26 | 0.6532 |  | 2.0 | 8617.2 |
| 1% | polynomial | 8 | yes | Multigraph | 11/0 | 20 | 0.6512 |  | 10.0 | 6343.7 |
| 1% | spectral | 12 | no | Modified MG | 18/0 | 27 | 0.1946 |  | 10.0 | 9983.3 |
| 1% | spectral | 12 | no | Multigraph | 17/0 | 26 | 0.2047 |  | 10.0 | 14597.8 |
| 1% | spectral | 12 | yes | Modified MG | 5/0 | 14 | 0.1524 |  | 10.0 | 8165.0 |
| 1% | spectral | 12 | yes | Multigraph | 17/0 | 26 | 0.2019 |  | 10.0 | 14610.5 |
| 1% | spectral | 4 | no | Modified MG | 0/9 | 0 | 1.0000 |  | 10.0 | 25820.3 |
| 1% | spectral | 4 | no | Multigraph | 6/2 | 13 | 1.0042 |  | 2.0 | 8664.7 |
| 1% | spectral | 4 | yes | Modified MG | 0/9 | 0 | 1.0000 |  | 10.0 | 65331.2 |
| 1% | spectral | 4 | yes | Multigraph | 0/9 | 0 | 1.0000 |  | 10.0 | 48729.1 |
| 1% | spectral | 8 | no | Modified MG | 16/0 | 25 | 0.6487 |  | 10.0 | 11663.8 |
| 1% | spectral | 8 | no | Multigraph | 14/0 | 23 | 0.6549 |  | 5.0 | 11349.1 |
| 1% | spectral | 8 | yes | Modified MG | 18/0 | 27 | 0.6475 |  | 10.0 | 7231.5 |
| 1% | spectral | 8 | yes | Multigraph | 14/0 | 23 | 0.6552 |  | 5.0 | 11392.3 |
| 2% | polynomial | 12 | no | Modified MG | 18/0 | 27 | 0.3860 |  | 5.0 | 4308.1 |
| 2% | polynomial | 12 | no | Multigraph | 12/0 | 21 | 0.3843 |  | 10.0 | 6863.3 |
| 2% | polynomial | 12 | yes | Modified MG | 18/0 | 27 | 0.3659 |  | 5.0 | 4280.4 |
| 2% | polynomial | 12 | yes | Multigraph | 14/0 | 23 | 0.3795 |  | 5.0 | 6858.7 |
| 2% | polynomial | 4 | no | Modified MG | 9/7 | 11 | 0.9986 |  | 10.0 | 22399.2 |
| 2% | polynomial | 4 | no | Multigraph | 9/3 | 15 | 0.9985 |  | 5.0 | 8383.3 |
| 2% | polynomial | 4 | yes | Modified MG | 0/9 | 0 | 1.0000 |  | 10.0 | 71225.6 |
| 2% | polynomial | 4 | yes | Multigraph | 9/3 | 15 | 0.9985 |  | 5.0 | 47663.2 |
| 2% | polynomial | 8 | no | Modified MG | 18/0 | 27 | 0.8652 |  | 2.0 | 8089.1 |
| 2% | polynomial | 8 | no | Multigraph | 14/0 | 23 | 0.8656 |  | 2.0 | 5957.7 |
| 2% | polynomial | 8 | yes | Modified MG | 18/0 | 27 | 0.8711 |  | 10.0 | 8016.8 |
| 2% | polynomial | 8 | yes | Multigraph | 13/0 | 22 | 0.8755 |  | 10.0 | 5947.1 |
| 2% | spectral | 12 | no | Modified MG | 18/0 | 27 | 0.3417 |  | 10.0 | 30683.0 |
| 2% | spectral | 12 | no | Multigraph | 18/0 | 27 | 0.3526 |  | 10.0 | 38098.6 |
| 2% | spectral | 12 | yes | Modified MG | 7/0 | 16 | 0.2787 |  | 10.0 | 22738.9 |
| 2% | spectral | 12 | yes | Multigraph | 18/0 | 27 | 0.3508 |  | 10.0 | 38140.9 |
| 2% | spectral | 4 | no | Modified MG | 6/4 | 11 | 1.0053 |  | 5.0 | 57493.2 |
| 2% | spectral | 4 | no | Multigraph | 6/2 | 13 | 1.0051 |  | 5.0 | 17640.7 |
| 2% | spectral | 4 | yes | Modified MG | 2/4 | 7 | 1.0035 |  | 10.0 | 55919.9 |
| 2% | spectral | 4 | yes | Multigraph | 0/9 | 0 | 1.0000 |  | 10.0 | 40979.0 |
| 2% | spectral | 8 | no | Modified MG | 18/0 | 27 | 0.8563 |  | 2.0 | 35744.6 |
| 2% | spectral | 8 | no | Multigraph | 17/0 | 26 | 0.8556 |  | 10.0 | 26461.4 |
| 2% | spectral | 8 | yes | Modified MG | 7/0 | 16 | 0.8613 |  | 10.0 | 16757.6 |
| 2% | spectral | 8 | yes | Multigraph | 17/0 | 26 | 0.8561 |  | 10.0 | 26587.8 |
| 5% | polynomial | 12 | no | Modified MG | 14/0 | 23 | 0.7305 |  | 10.0 | 4557.8 |
| 5% | polynomial | 12 | no | Multigraph | 14/0 | 23 | 0.7374 |  | 5.0 | 7321.3 |
| 5% | polynomial | 12 | yes | Modified MG | 18/0 | 27 | 0.7132 |  | 10.0 | 4561.2 |
| 5% | polynomial | 12 | yes | Multigraph | 14/0 | 23 | 0.7385 |  | 5.0 | 7321.0 |
| 5% | polynomial | 4 | no | Modified MG | 0/9 | 0 | 1.0000 |  | 10.0 | 22432.0 |
| 5% | polynomial | 4 | no | Multigraph | 6/4 | 11 | 1.0090 |  | 5.0 | 9678.6 |
| 5% | polynomial | 4 | yes | Modified MG | 0/9 | 0 | 1.0000 |  | 10.0 | 64539.5 |
| 5% | polynomial | 4 | yes | Multigraph | 0/9 | 0 | 1.0000 |  | 10.0 | 38330.7 |
| 5% | polynomial | 8 | no | Modified MG | 11/3 | 17 | 0.9970 |  | 5.0 | 7256.2 |
| 5% | polynomial | 8 | no | Multigraph | 11/3 | 17 | 0.9970 |  | 5.0 | 5943.7 |
| 5% | polynomial | 8 | yes | Modified MG | 11/3 | 17 | 0.9970 |  | 5.0 | 7183.1 |
| 5% | polynomial | 8 | yes | Multigraph | 11/3 | 17 | 0.9970 |  | 5.0 | 5945.2 |
| 5% | spectral | 12 | no | Modified MG | 35/0 | 44 | 0.7659 |  | 0.0001 | 120972.7 |
| 5% | spectral | 12 | no | Multigraph | 35/0 | 44 | 0.7659 |  | 0.0001 | 141216.4 |
| 5% | spectral | 12 | yes | Modified MG | 18/0 | 27 | 0.7025 |  | 5.0 | 94530.9 |
| 5% | spectral | 12 | yes | Multigraph | 35/0 | 44 | 0.7656 |  | 0.0001 | 141328.9 |
| 5% | spectral | 4 | no | Modified MG | 0/9 | 0 | 1.0000 |  | 10.0 | 240792.5 |
| 5% | spectral | 4 | no | Multigraph | 0/9 | 0 | 1.0000 |  | 10.0 | 73926.3 |
| 5% | spectral | 4 | yes | Modified MG | 0/9 | 0 | 1.0000 |  | 10.0 | 64422.5 |
| 5% | spectral | 4 | yes | Multigraph | 0/9 | 0 | 1.0000 |  | 10.0 | 79797.7 |
| 5% | spectral | 8 | no | Modified MG | 35/0 | 44 | 0.9942 |  | 0.0001 | 144478.5 |
| 5% | spectral | 8 | no | Multigraph | 35/0 | 44 | 0.9942 |  | 0.0001 | 98139.8 |
| 5% | spectral | 8 | yes | Modified MG | 2/7 | 4 | 0.9983 |  | 5.0 | 68587.1 |
| 5% | spectral | 8 | yes | Multigraph | 35/0 | 44 | 0.9944 |  | 0.0001 | 98410.2 |

---
*Source: `mg_comparison_20260629_165010_full.csv`*
