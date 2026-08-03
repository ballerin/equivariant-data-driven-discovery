# Kuramoto-Sivashinsky equation (weakly chaotic, $L=8\pi$)

The equation under consideration is

$$f_t = -\Delta f - \Delta^2 f - \tfrac12\|\nabla f\|^2$$

which can be written in terms of $O(2)$-equivariant terms

| Operator         | True coefficient |
|:-----------------|:-----------------|
| $\Delta f$       | $-1.0$           |
| $\Delta^2 f$     | $-1.0$           |
| $\|\nabla f\|^2$ | $-0.5$           |

or standard terms (PolyJet basis)

| Operator            | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|:--------------------|:-----------|:---------|:-------|:---------|:-----------|
| $\partial_{xx} f$   |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yy} f$   |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxxx} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxyy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yyyy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f_x^2$             |            | 🔴        | 🔴      | 🟢        | 🟢          |
| $f_y^2$             |            | 🔴        | 🔴      | 🟢        | 🟢          |

KS contains the gradient norm $\|\nabla f\|^2 = f_x^2 + f_y^2$ which neither PDE-FIND nor WSINDy can represent as a sparse column (PDE-FIND lacks products of derivatives; in WSINDy, $f_x^2$ is locked with $f\cdot f_{xx}$ inside $\partial_{xx}(f^2)$).  Weakly chaotic regime ($L=8\pi$).

### Legend for summary table

| Symbol | Meaning |
|---|---|
| 🟢 | Equivariant **clear** win (fewer total fp+fn) |
| 🟩 | Equivariant **marginal** win (same fp+fn, lower E) |
| ⚪ | Tie (same fp+fn, same E) |
| 🟥 | Non-equivariant **marginal** win (same fp+fn, lower E) |
| 🔴 | Non-equivariant **clear** win (fewer total fp+fn) |
| — | Both methods failed, no comparison |
| ✓ | Exact recovery (fp = fn = 0) |
| ✗ | Failed |


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best         | Mg E_data | Std best        | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-------------------|:----------|:----------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | ⚪        | 0/0 ✓             | 0.0047     | multigraph/RFE     | 0.0176    | PolyJet/RFE     | 0/0 ✓     | 0.0047      | 0.0175     | —            |
| Clean | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.0117     | multigraph/RFE     | 0.0325    | PolyJet/RFE     | 0/0 ✓     | 0.0128      | 0.0320     | —            |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0021     | multigraph/STRidge | 0.0112    | PolyJet/STRidge | 0/0 ✓     | 0.0023      | 0.0112     | —            |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0022     | multigraph/MSTLS   | 0.0112    | PolyJet/MSTLS   | 0/0 ✓     | 0.0023      | 0.0112     | —            |
| Clean | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0029     | multigraph/STRidge | 0.0133    | PolyJet/STRidge | 0/0 ✓     | 0.0030      | 0.0133     | —            |
| Clean | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0029     | multigraph/STRidge | 0.0133    | PolyJet/STRidge | 0/0 ✓     | 0.0029      | 0.0133     | —            |
| 2%    | pointwise     | spectral   | —        | 18/0 ✗            | 0.4266     | —                  | 0.6378    | —               | 83/0 ✗    | 0.4249      | 0.6367     | —            |
| 2%    | pointwise     | polynomial | —        | 28/0 ✗            | 0.3663     | —                  | 0.6210    | —               | 79/0 ✗    | 0.4291      | 0.6371     | —            |
| 2%    | weak (no IBP) | spectral   | —        | 6/0 ✗             | 0.5485     | —                  | 0.3476    | —               | 62/0 ✗    | 0.4969      | 0.3239     | —            |
| 2%    | weak (no IBP) | polynomial | —        | 7/0 ✗             | 0.1307     | —                  | 0.1168    | —               | 9/0 ✗     | 0.0954      | 0.1099     | —            |
| 2%    | weak (IBP)    | spectral   | —        | 30/0 ✗            | 0.4404     | —                  | 0.3931    | —               | 11/0 ✗    | 0.5365      | 0.3255     | —            |
| 2%    | weak (IBP)    | polynomial | —        | 9/0 ✗             | 0.3100     | —                  | 0.2168    | —               | 33/0 ✗    | 0.2693      | 0.2026     | —            |
| 5%    | pointwise     | spectral   | —        | 28/0 ✗            | 0.4318     | —                  | 0.8725    | —               | 89/0 ✗    | 0.4310      | 0.8722     | —            |
| 5%    | pointwise     | polynomial | —        | 33/0 ✗            | 0.4222     | —                  | 0.8717    | —               | 93/0 ✗    | 0.4194      | 0.8702     | —            |
| 5%    | weak (no IBP) | spectral   | —        | 30/0 ✗            | 0.4313     | —                  | 0.5369    | —               | 59/0 ✗    | 0.4459      | 0.4115     | —            |
| 5%    | weak (no IBP) | polynomial | —        | 8/0 ✗             | 0.3897     | —                  | 0.2940    | —               | 13/0 ✗    | 0.2512      | 0.2611     | —            |
| 5%    | weak (IBP)    | spectral   | —        | 37/0 ✗            | 0.4345     | —                  | 0.4875    | —               | 78/0 ✗    | 0.4340      | 0.4833     | —            |
| 5%    | weak (IBP)    | polynomial | —        | 37/0 ✗            | 0.3803     | —                  | 0.4640    | —               | 17/0 ✗    | 0.3406      | 0.3414     | —            |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                | Best conf.                     | fp/fn  | E_coeff     | Time (s) | vs. best |
|:----------------------|:-------------------------------|:-------|------------:|---------:|---------:|
| Multigraphs + STRidge | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0021      | 0.3      | —        |
| Multigraphs + MSTLS   | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0022      | 0.3      | +3.6%    |
| Multigraphs + RFE     | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0022      | 0.1      | +3.6%    |
| PolyJet + STRidge     | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0023      | 1.3      | +6.9%    |
| PolyJet + MSTLS       | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0023      | 3.7      | +8.3%    |
| PolyJet + RFE         | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0023      | 0.1      | +8.3%    |
| WSINDy + STRidge      | spectral, weak (no IBP), m=8   | 76/0 ✗ | 0.4334      | 1.7      | ✗        |
| PDE-FIND + STRidge    | polynomial, pointwise          | 24/2 ✗ | 0.4288      | 2.9      | ✗        |
| PDE-FIND + MSTLS      | spectral, weak (IBP), m=12     | 7/5 ✗  | 0.4412      | 1.5      | ✗        |
| PDE-FIND + RFE        | spectral, pointwise            | 1/7 ✗  | 0.4409      | 0.1      | ✗        |
| WSINDy + RFE          | spectral, pointwise            | 1/7 ✗  | 0.4409      | 0.1      | ✗        |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=4 | 4/7 ✗  | 444514.6420 | 2.1      | ✗        |

### 2% noise

| Method                | Best conf.                      | fp/fn  | E_coeff     | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|------------:|---------:|---------:|
| Multigraphs + MSTLS   | spectral, weak (no IBP), m=8    | 6/0 ✗  | 0.5485      | 0.4      | ✗        |
| Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 7/0 ✗  | 0.1341      | 0.4      | ✗        |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=12 | 9/0 ✗  | 0.0954      | 4.0      | ✗        |
| PolyJet + STRidge     | polynomial, weak (no IBP), m=8  | 10/0 ✗ | 0.2682      | 2.4      | ✗        |
| WSINDy + STRidge      | polynomial, pointwise           | 79/0 ✗ | 0.4291      | 2.9      | ✗        |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=4    | 23/2 ✗ | 0.4283      | 2.7      | ✗        |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.5578      | 0.1      | ✗        |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12   | 7/5 ✗  | 0.4447      | 1.6      | ✗        |
| PolyJet + RFE         | polynomial, weak (no IBP), m=12 | 0/6 ✗  | 0.7043      | 0.1      | ✗        |
| PDE-FIND + RFE        | spectral, pointwise             | 1/7 ✗  | 0.4409      | 0.1      | ✗        |
| WSINDy + RFE          | spectral, pointwise             | 1/7 ✗  | 0.4409      | 0.1      | ✗        |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=4  | 4/7 ✗  | 438169.9763 | 2.2      | ✗        |

### 5% noise

| Method                | Best conf.                      | fp/fn  | E_coeff     | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|------------:|---------:|---------:|
| Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 8/0 ✗  | 0.3897      | 0.3      | ✗        |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=12 | 13/0 ✗ | 0.2512      | 3.9      | ✗        |
| Multigraphs + STRidge | polynomial, weak (no IBP), m=4  | 22/0 ✗ | 0.3894      | 0.9      | ✗        |
| PolyJet + STRidge     | polynomial, weak (no IBP), m=12 | 37/0 ✗ | 0.2682      | 2.6      | ✗        |
| WSINDy + STRidge      | polynomial, weak (no IBP), m=4  | 85/0 ✗ | 0.4292      | 2.8      | ✗        |
| PDE-FIND + STRidge    | polynomial, weak (IBP), m=12    | 11/2 ✗ | 0.4347      | 1.0      | ✗        |
| Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/5 ✗  | 0.8602      | 0.1      | ✗        |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=4    | 1/7 ✗  | 0.4409      | 3.3      | ✗        |
| PDE-FIND + RFE        | spectral, pointwise             | 1/7 ✗  | 0.5796      | 0.1      | ✗        |
| WSINDy + RFE          | spectral, pointwise             | 1/7 ✗  | 0.5796      | 0.1      | ✗        |
| PolyJet + RFE         | spectral, pointwise             | 1/7 ✗  | 0.5796      | 0.2      | ✗        |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=4  | 4/7 ✗  | 408758.5378 | 2.1      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.001   | 55/2  | 60       | 0.4272  | 0.4146 | 5.0      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 1/7   | 1        | 0.4409  | 0.4294 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.001   | 54/2  | 59       | 0.4292  | 0.6377 | 4.7      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.01    | 47/2  | 52       | 0.4308  | 0.8720 | 4.5      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.001   | 83/0  | 40       | 0.4280  | 0.4157 | 2.8      | 🔴      |
| WSINDy      | Clean | RFE       |         | 1/7   | 1        | 0.4409  | 0.4294 | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.01    | 94/0  | 45       | 0.4291  | 0.6378 | 2.9      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.01    | 89/0  | 44       | 0.4310  | 0.8722 | 2.9      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0047  | 0.0175 | 3.1      | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.0047  | 0.0175 | 0.2      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.01    | 83/0  | 90       | 0.4249  | 0.6367 | 11.5     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.01    | 91/0  | 98       | 0.4309  | 0.8714 | 12.2     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.2      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0047  | 0.0176 | 0.5      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0047  | 0.0176 | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.1     | 18/0  | 12       | 0.4266  | 0.6378 | 0.7      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.01    | 28/0  | 16       | 0.4318  | 0.8725 | 0.8      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.1      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 24/2  | 29       | 0.4288  | 0.4158 | 2.9      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 1/7   | 1        | 0.4416  | 0.4294 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.001   | 53/2  | 58       | 0.4294  | 0.6374 | 4.8      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/7   | 1        | 0.4417  | 0.6433 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 55/2  | 60       | 0.4304  | 0.8716 | 4.6      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/7   | 1        | 0.5800  | 0.8979 | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.001   | 83/0  | 40       | 0.4288  | 0.4157 | 2.9      | 🔴      |
| WSINDy      | Clean | RFE       |         | 1/7   | 1        | 0.4416  | 0.4294 | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.001   | 79/0  | 41       | 0.4291  | 0.6371 | 2.9      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 1/7   | 1        | 0.4417  | 0.6433 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.01    | 93/0  | 44       | 0.4305  | 0.8718 | 2.9      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/7   | 1        | 0.5800  | 0.8979 | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0130  | 0.0320 | 3.3      | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.0128  | 0.0320 | 0.2      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.001   | 95/0  | 102      | 0.3606  | 0.6182 | 12.6     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/7   | 1        | 0.4417  | 0.6433 | 0.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.01    | 93/0  | 100      | 0.4194  | 0.8702 | 13.5     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/7   | 1        | 0.5800  | 0.8979 | 0.2      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0119  | 0.0325 | 0.5      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0117  | 0.0325 | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.01    | 28/0  | 16       | 0.3663  | 0.6210 | 0.8      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.4417  | 0.6433 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.001   | 33/0  | 18       | 0.4222  | 0.8717 | 0.9      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.5800  | 0.8979 | 0.1      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff     | E_data      | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:------------|:------------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 4  | no  | 43/2  | 48       | 0.4284      | 0.4002      | 3.9      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.4173      | 3.2      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5615      | 0.1      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 0.001   | 8  | yes | 55/2  | 60       | 0.4322      | 0.3888      | 3.2      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 12 | yes | 7/5   | 9        | 0.4412      | 0.3678      | 1.5      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5526      | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 4  | no  | 23/2  | 28       | 0.4283      | 0.4295      | 2.7      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | no  | 7/5   | 9        | 0.4447      | 0.3361      | 1.6      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5778      | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 8  | yes | 28/2  | 33       | 0.4311      | 0.4095      | 2.0      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4409      | 0.4266      | 2.3      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5630      | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | no  | 34/2  | 39       | 0.4432      | 0.3506      | 1.5      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.5449      | 3.3      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 8  | yes | 22/2  | 27       | 0.4310      | 0.4851      | 1.5      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4409      | 0.4974      | 2.2      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.6098      | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 1.0     | 8  | no  | 76/0  | 38       | 0.4334      | 0.3644      | 1.7      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 444689.2112 | 448874.7916 | 2.1      | 🔴      |
| WSINDy      | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5615      | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 8  | yes | 90/0  | 44       | 0.4318      | 0.3903      | 1.9      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 663202.9871 | 668683.6885 | 1.5      | 🔴      |
| WSINDy      | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5526      | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.01    | 12 | no  | 88/0  | 43       | 0.4506      | 0.3247      | 1.3      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 438342.0538 | 436503.6332 | 2.1      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5778      | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | yes | 85/0  | 43       | 0.4359      | 0.3652      | 1.3      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 653069.4089 | 653081.9030 | 1.5      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5630      | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 8  | no  | 85/0  | 43       | 0.4302      | 0.4131      | 1.9      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 408919.0648 | 381256.0793 | 2.1      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.01    | 8  | yes | 94/0  | 45       | 0.4298      | 0.4844      | 2.0      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 618999.3579 | 594279.6067 | 1.5      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.6098      | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 1.3      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 3.7      | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 0.1      | 🟢      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0030      | 0.0133      | 1.3      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 0.0030      | 0.0133      | 3.6      | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0030      | 0.0133      | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.1     | 8  | no  | 62/0  | 69       | 0.4969      | 0.3239      | 6.1      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 4  | no  | 12/1  | 18       | 0.4233      | 0.4182      | 8.8      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 0/6   | 1        | 0.9275      | 0.3548      | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 0.1     | 8  | yes | 37/0  | 44       | 0.4143      | 0.3812      | 4.3      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.5     | 12 | yes | 11/0  | 18       | 0.5365      | 0.3255      | 3.9      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 0.9288      | 0.4103      | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 8  | no  | 59/0  | 66       | 0.4459      | 0.4115      | 5.6      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | no  | 9/3   | 13       | 0.4430      | 0.3575      | 3.6      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 8  | yes | 78/0  | 85       | 0.4340      | 0.4833      | 5.8      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4409      | 0.4974      | 5.4      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.6098      | 0.2      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0021      | 0.0112      | 0.3      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.3      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0029      | 0.0133      | 0.3      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0030      | 0.0133      | 0.3      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0030      | 0.0133      | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.1     | 4  | no  | 30/0  | 17       | 0.4212      | 0.4176      | 0.8      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 6/0   | 6        | 0.5485      | 0.3476      | 0.4      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9250      | 0.3445      | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 0.1     | 8  | yes | 30/0  | 17       | 0.4404      | 0.3931      | 0.6      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 8  | yes | 6/3   | 5        | 0.4485      | 0.3981      | 0.4      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9259      | 0.3931      | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.1     | 4  | no  | 30/0  | 17       | 0.4313      | 0.5369      | 0.8      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.5449      | 0.6      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.001   | 8  | yes | 37/0  | 20       | 0.4345      | 0.4875      | 0.6      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4409      | 0.4974      | 0.4      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.6098      | 0.1      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff     | E_data      | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:------------|:------------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 4  | no  | 43/2  | 48       | 0.4291      | 0.4002      | 4.2      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4416      | 0.4173      | 3.3      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.5615      | 0.1      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 0.001   | 8  | yes | 54/2  | 59       | 0.4327      | 0.3889      | 3.2      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 12 | yes | 7/5   | 9        | 0.4419      | 0.3678      | 1.6      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.5526      | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | no  | 27/2  | 32       | 0.4398      | 0.3266      | 1.4      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4416      | 0.4432      | 3.1      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.5778      | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 12 | yes | 41/2  | 46       | 0.4368      | 0.3649      | 1.7      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4416      | 0.4266      | 2.2      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.5630      | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | no  | 27/2  | 32       | 0.4374      | 0.3491      | 1.4      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4416      | 0.5449      | 3.4      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.6460      | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | yes | 11/2  | 16       | 0.4347      | 0.4060      | 1.0      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4416      | 0.4974      | 2.2      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.6098      | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 1.0     | 8  | no  | 76/0  | 38       | 0.4336      | 0.3644      | 1.8      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 444514.6420 | 448874.7916 | 2.1      | 🔴      |
| WSINDy      | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.5615      | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.01    | 8  | yes | 83/0  | 41       | 0.4321      | 0.3906      | 2.0      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 662942.6372 | 668683.6885 | 1.6      | 🔴      |
| WSINDy      | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.5526      | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.01    | 12 | no  | 88/0  | 43       | 0.4503      | 0.3247      | 1.3      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 438169.9763 | 436503.6332 | 2.2      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.5778      | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.01    | 8  | yes | 95/0  | 45       | 0.4313      | 0.4096      | 2.1      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 652813.0372 | 653081.9030 | 1.5      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.5630      | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.01    | 4  | no  | 85/0  | 44       | 0.4292      | 0.5353      | 2.8      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 408758.5378 | 381256.0793 | 2.1      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.6460      | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.01    | 12 | yes | 90/0  | 44       | 0.4343      | 0.4018      | 1.3      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 618756.3608 | 594279.6067 | 1.5      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.6098      | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0024      | 0.0112      | 1.4      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 3.7      | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 0.1      | 🟢      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0029      | 0.0133      | 1.3      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 0.0030      | 0.0133      | 3.8      | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0030      | 0.0133      | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 10.0    | 8  | no  | 10/0  | 17       | 0.2682      | 0.2197      | 2.4      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.5     | 12 | no  | 9/0   | 16       | 0.0954      | 0.1099      | 4.0      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 0/6   | 1        | 0.7043      | 0.3248      | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | yes | 33/0  | 40       | 0.2693      | 0.2026      | 2.9      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 10/1  | 16       | 0.2877      | 0.2122      | 3.9      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 0.7062      | 0.4060      | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 10.0    | 12 | no  | 37/0  | 44       | 0.2682      | 0.2533      | 2.6      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.5     | 12 | no  | 13/0  | 20       | 0.2512      | 0.2611      | 3.9      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.6460      | 0.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | yes | 58/0  | 65       | 0.3374      | 0.3361      | 3.2      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | yes | 17/0  | 24       | 0.3406      | 0.3414      | 3.7      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.6098      | 0.2      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0039      | 0.0141      | 0.4      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.3      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0029      | 0.0133      | 0.3      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0029      | 0.0133      | 0.3      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0029      | 0.0133      | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 10.0    | 12 | no  | 7/0   | 5        | 0.1341      | 0.1168      | 0.4      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 12 | no  | 7/0   | 5        | 0.1307      | 0.1168      | 0.3      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.5578      | 0.2228      | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | yes | 33/0  | 18       | 0.3000      | 0.2119      | 0.5      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 9/0   | 7        | 0.3100      | 0.2168      | 0.3      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.5617      | 0.2773      | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 4  | no  | 22/0  | 14       | 0.3894      | 0.5225      | 0.9      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | no  | 8/0   | 6        | 0.3897      | 0.2940      | 0.3      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.6460      | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.001   | 8  | yes | 37/0  | 20       | 0.3803      | 0.4640      | 0.6      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | yes | 9/2   | 6        | 0.4026      | 0.3677      | 0.3      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.8602      | 0.4354      | 0.1      | 🔴      |

---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ks_wchaos`.

| Source                                     | Wall time | CPU                                      | Threads |
|:-------------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260802_081855_ks_wchaos_full.csv` | 163.3 min | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | 1       |
