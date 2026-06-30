# KS (chaotic) — Multigraph vs Modified Multigraph

Equation: $f_t = -\Delta f - \Delta^2 f - \|\nabla f\|^2/2$

| Library | True terms |
|---|---|
| Multigraph | Δf, ‖∇f‖², Δ²f |
| Modified MG | Δf, ‖∇f‖², Δ²f |

## Best per configuration (STRidge)

| Noise | Backend | m | IBP | MG fp/fn | MG E | Mod fp/fn | Mod E | Winner |
|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | 15/0 | 0.0589 | 4/0 | 0.0579 | 🟢 |
| clean | polynomial | 12 | yes | 11/0 | 0.0552 | 4/0 | 0.0549 | 🟢 |
| clean | polynomial | 4 | no | 7/0 | 0.2787 | 7/0 | 0.2667 | 🟩 |
| clean | polynomial | 4 | yes | 11/0 | 11.7479 | 21/0 | 14.8292 | 🔴 |
| clean | polynomial | 8 | no | 7/0 | 0.1336 | 0/0 | 0.0485 | 🟢 |
| clean | polynomial | 8 | yes | 8/0 | 0.1019 | 1/0 | 0.0578 | 🟢 |
| clean | spectral | 12 | no | 7/0 | 0.1753 | 1/0 | 0.0068 | 🟢 |
| clean | spectral | 12 | yes | 7/0 | 0.1555 | 0/0 | 0.0059 | 🟢 |
| clean | spectral | 4 | no | 0/0 | 0.0037 | 0/0 | 0.0037 | ⚪ |
| clean | spectral | 4 | yes | 25/0 | 7.5028 | 30/0 | 14.6208 | 🔴 |
| clean | spectral | 8 | no | 37/0 | 0.0534 | 0/0 | 0.0017 | 🟢 |
| clean | spectral | 8 | yes | 15/0 | 0.0706 | 4/0 | 0.0509 | 🟢 |
| 2% | polynomial | 12 | no | 27/0 | 0.4204 | 22/0 | 0.4184 | 🟢 |
| 2% | polynomial | 12 | yes | 27/0 | 0.4221 | 16/0 | 0.4192 | 🟢 |
| 2% | polynomial | 4 | no | 8/0 | 0.9135 | 9/0 | 0.9055 | 🔴 |
| 2% | polynomial | 4 | yes | 18/5 | 0.9220 | 9/2 | 1.0719 | 🟢 |
| 2% | polynomial | 8 | no | 19/0 | 0.4929 | 19/0 | 0.4931 | 🟥 |
| 2% | polynomial | 8 | yes | 37/0 | 0.4900 | 20/0 | 0.6454 | 🟢 |
| 2% | spectral | 12 | no | 32/0 | 0.4032 | 18/0 | 0.3806 | 🟢 |
| 2% | spectral | 12 | yes | 32/0 | 0.4031 | 37/0 | 0.6291 | 🔴 |
| 2% | spectral | 4 | no | 8/5 | 0.9996 | 11/5 | 0.9996 | 🔴 |
| 2% | spectral | 4 | yes | 37/0 | 0.9689 | 12/5 | 1.0044 | 🟢 |
| 2% | spectral | 8 | no | 30/0 | 0.4555 | 32/0 | 0.4546 | 🔴 |
| 2% | spectral | 8 | yes | 16/0 | 0.4412 | 37/0 | 0.6910 | 🔴 |
| 5% | polynomial | 12 | no | 37/0 | 1.0160 | 37/0 | 1.0160 | 🟥 |
| 5% | polynomial | 12 | yes | 37/0 | 1.0139 | 37/0 | 1.6102 | 🟥 |
| 5% | polynomial | 4 | no | 37/0 | 1.0352 | 37/0 | 1.0352 | 🟩 |
| 5% | polynomial | 4 | yes | 37/0 | 1.0383 | 11/5 | 1.0724 | 🟢 |
| 5% | polynomial | 8 | no | 37/0 | 0.8154 | 37/0 | 0.8154 | 🟩 |
| 5% | polynomial | 8 | yes | 37/0 | 0.8089 | 23/0 | 1.0909 | 🟢 |
| 5% | spectral | 12 | no | 35/2 | 1.0753 | 34/0 | 1.0591 | 🟢 |
| 5% | spectral | 12 | yes | 35/2 | 1.0758 | 37/0 | 1.1915 | 🟥 |
| 5% | spectral | 4 | no | 8/5 | 0.9997 | 11/5 | 0.9997 | 🔴 |
| 5% | spectral | 4 | yes | 14/5 | 0.9997 | 4/5 | 0.9997 | 🟢 |
| 5% | spectral | 8 | no | 30/5 | 0.9468 | 37/0 | 0.9480 | 🔴 |
| 5% | spectral | 8 | yes | 30/5 | 0.9300 | 32/2 | 0.9588 | 🟢 |

## Detailed results

| Noise | Backend | m | IBP | Library | fp/fn | n_act | E | Exact? | d_tol | cond |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | polynomial | 12 | no | Modified MG | 4/0 | 11 | 0.0579 |  | 10.0 | 317272517.4 |
| clean | polynomial | 12 | no | Multigraph | 15/0 | 22 | 0.0589 |  | 5.0 | 318159697.9 |
| clean | polynomial | 12 | yes | Modified MG | 4/0 | 11 | 0.0549 |  | 10.0 | 315815967.5 |
| clean | polynomial | 12 | yes | Multigraph | 11/0 | 18 | 0.0552 |  | 5.0 | 315952930.1 |
| clean | polynomial | 4 | no | Modified MG | 7/0 | 14 | 0.2667 |  | 10.0 | 39747181.8 |
| clean | polynomial | 4 | no | Multigraph | 7/0 | 14 | 0.2787 |  | 5.0 | 39807836.4 |
| clean | polynomial | 4 | yes | Modified MG | 21/0 | 28 | 14.8292 |  | 5.0 | 63961636.1 |
| clean | polynomial | 4 | yes | Multigraph | 11/0 | 18 | 11.7479 |  | 10.0 | 40706924.2 |
| clean | polynomial | 8 | no | Modified MG | 0/0 | 7 | 0.0485 | ✓ | 5.0 | 107952639.0 |
| clean | polynomial | 8 | no | Multigraph | 7/0 | 14 | 0.1336 |  | 5.0 | 108519130.9 |
| clean | polynomial | 8 | yes | Modified MG | 1/0 | 8 | 0.0578 |  | 10.0 | 105364877.7 |
| clean | polynomial | 8 | yes | Multigraph | 8/0 | 15 | 0.1019 |  | 5.0 | 104746960.5 |
| clean | spectral | 12 | no | Modified MG | 1/0 | 8 | 0.0068 |  | 10.0 | 315469761.6 |
| clean | spectral | 12 | no | Multigraph | 7/0 | 14 | 0.1753 |  | 10.0 | 317175677.7 |
| clean | spectral | 12 | yes | Modified MG | 0/0 | 7 | 0.0059 | ✓ | 10.0 | 315772580.2 |
| clean | spectral | 12 | yes | Multigraph | 7/0 | 14 | 0.1555 |  | 10.0 | 317477220.6 |
| clean | spectral | 4 | no | Modified MG | 0/0 | 7 | 0.0037 | ✓ | 5.0 | 36725677.3 |
| clean | spectral | 4 | no | Multigraph | 0/0 | 7 | 0.0037 | ✓ | 5.0 | 37240999.0 |
| clean | spectral | 4 | yes | Modified MG | 30/0 | 37 | 14.6208 |  | 5.0 | 62445876.5 |
| clean | spectral | 4 | yes | Multigraph | 25/0 | 32 | 7.5028 |  | 5.0 | 46753006.8 |
| clean | spectral | 8 | no | Modified MG | 0/0 | 7 | 0.0017 | ✓ | 10.0 | 103575094.2 |
| clean | spectral | 8 | no | Multigraph | 37/0 | 44 | 0.0534 |  | 0.0001 | 103883766.1 |
| clean | spectral | 8 | yes | Modified MG | 4/0 | 11 | 0.0509 |  | 10.0 | 104713347.1 |
| clean | spectral | 8 | yes | Multigraph | 15/0 | 22 | 0.0706 |  | 10.0 | 105059979.5 |
| 2% | polynomial | 12 | no | Modified MG | 22/0 | 29 | 0.4184 |  | 10.0 | 304180159.7 |
| 2% | polynomial | 12 | no | Multigraph | 27/0 | 34 | 0.4204 |  | 10.0 | 304727897.1 |
| 2% | polynomial | 12 | yes | Modified MG | 16/0 | 23 | 0.4192 |  | 10.0 | 312794249.7 |
| 2% | polynomial | 12 | yes | Multigraph | 27/0 | 34 | 0.4221 |  | 10.0 | 302805744.9 |
| 2% | polynomial | 4 | no | Modified MG | 9/0 | 16 | 0.9055 |  | 10.0 | 36801244.3 |
| 2% | polynomial | 4 | no | Multigraph | 8/0 | 15 | 0.9135 |  | 10.0 | 36907183.3 |
| 2% | polynomial | 4 | yes | Modified MG | 9/2 | 14 | 1.0719 |  | 10.0 | 38413184.7 |
| 2% | polynomial | 4 | yes | Multigraph | 18/5 | 20 | 0.9220 |  | 10.0 | 36906281.6 |
| 2% | polynomial | 8 | no | Modified MG | 19/0 | 26 | 0.4931 |  | 10.0 | 100530058.1 |
| 2% | polynomial | 8 | no | Multigraph | 19/0 | 26 | 0.4929 |  | 10.0 | 100564243.2 |
| 2% | polynomial | 8 | yes | Modified MG | 20/0 | 27 | 0.6454 |  | 10.0 | 98330150.5 |
| 2% | polynomial | 8 | yes | Multigraph | 37/0 | 44 | 0.4900 |  | 0.0001 | 91189739.9 |
| 2% | spectral | 12 | no | Modified MG | 18/0 | 25 | 0.3806 |  | 10.0 | 273987602.4 |
| 2% | spectral | 12 | no | Multigraph | 32/0 | 39 | 0.4032 |  | 5.0 | 273651905.4 |
| 2% | spectral | 12 | yes | Modified MG | 37/0 | 44 | 0.6291 |  | 0.0001 | 311903338.0 |
| 2% | spectral | 12 | yes | Multigraph | 32/0 | 39 | 0.4031 |  | 5.0 | 273959604.1 |
| 2% | spectral | 4 | no | Modified MG | 11/5 | 13 | 0.9996 |  | 10.0 | 39790864.7 |
| 2% | spectral | 4 | no | Multigraph | 8/5 | 10 | 0.9996 |  | 10.0 | 39795395.5 |
| 2% | spectral | 4 | yes | Modified MG | 12/5 | 14 | 1.0044 |  | 10.0 | 46486696.7 |
| 2% | spectral | 4 | yes | Multigraph | 37/0 | 44 | 0.9689 |  | 0.0001 | 39797489.0 |
| 2% | spectral | 8 | no | Modified MG | 32/0 | 39 | 0.4546 |  | 5.0 | 90927192.0 |
| 2% | spectral | 8 | no | Multigraph | 30/0 | 37 | 0.4555 |  | 2.0 | 90842207.1 |
| 2% | spectral | 8 | yes | Modified MG | 37/0 | 44 | 0.6910 |  | 0.0001 | 96236633.1 |
| 2% | spectral | 8 | yes | Multigraph | 16/0 | 23 | 0.4412 |  | 10.0 | 91745285.8 |
| 5% | polynomial | 12 | no | Modified MG | 37/0 | 44 | 1.0160 |  | 0.0001 | 179003676.5 |
| 5% | polynomial | 12 | no | Multigraph | 37/0 | 44 | 1.0160 |  | 0.0001 | 179198493.1 |
| 5% | polynomial | 12 | yes | Modified MG | 37/0 | 44 | 1.6102 |  | 0.0001 | 302560144.7 |
| 5% | polynomial | 12 | yes | Multigraph | 37/0 | 44 | 1.0139 |  | 0.0001 | 178146133.7 |
| 5% | polynomial | 4 | no | Modified MG | 37/0 | 44 | 1.0352 |  | 0.0001 | 37585479.7 |
| 5% | polynomial | 4 | no | Multigraph | 37/0 | 44 | 1.0352 |  | 0.0001 | 37606173.9 |
| 5% | polynomial | 4 | yes | Modified MG | 11/5 | 13 | 1.0724 |  | 5.0 | 39698520.8 |
| 5% | polynomial | 4 | yes | Multigraph | 37/0 | 44 | 1.0383 |  | 0.0001 | 37557011.2 |
| 5% | polynomial | 8 | no | Modified MG | 37/0 | 44 | 0.8154 |  | 0.0001 | 75876089.7 |
| 5% | polynomial | 8 | no | Multigraph | 37/0 | 44 | 0.8154 |  | 0.0001 | 75976318.7 |
| 5% | polynomial | 8 | yes | Modified MG | 23/0 | 30 | 1.0909 |  | 10.0 | 82931219.8 |
| 5% | polynomial | 8 | yes | Multigraph | 37/0 | 44 | 0.8089 |  | 0.0001 | 64171012.7 |
| 5% | spectral | 12 | no | Modified MG | 34/0 | 41 | 1.0591 |  | 0.5 | 123485427.5 |
| 5% | spectral | 12 | no | Multigraph | 35/2 | 40 | 1.0753 |  | 2.0 | 124183762.7 |
| 5% | spectral | 12 | yes | Modified MG | 37/0 | 44 | 1.1915 |  | 5.0 | 303223049.1 |
| 5% | spectral | 12 | yes | Multigraph | 35/2 | 40 | 1.0758 |  | 2.0 | 124299300.8 |
| 5% | spectral | 4 | no | Modified MG | 11/5 | 13 | 0.9997 |  | 5.0 | 39643195.3 |
| 5% | spectral | 4 | no | Multigraph | 8/5 | 10 | 0.9997 |  | 10.0 | 39785550.0 |
| 5% | spectral | 4 | yes | Modified MG | 4/5 | 6 | 0.9997 |  | 10.0 | 44846864.9 |
| 5% | spectral | 4 | yes | Multigraph | 14/5 | 16 | 0.9997 |  | 10.0 | 39566356.0 |
| 5% | spectral | 8 | no | Modified MG | 37/0 | 44 | 0.9480 |  | 0.0001 | 58621260.5 |
| 5% | spectral | 8 | no | Multigraph | 30/5 | 32 | 0.9468 |  | 10.0 | 58952277.8 |
| 5% | spectral | 8 | yes | Modified MG | 32/2 | 37 | 0.9588 |  | 10.0 | 82479846.7 |
| 5% | spectral | 8 | yes | Multigraph | 30/5 | 32 | 0.9300 |  | 10.0 | 59009394.7 |

---
*Source: `mg_comparison_20260629_165010_full.csv`*
