# Cahn-Hilliard equation

The equation under consideration is

$$f_t = -\Delta^2 f + \Delta(f^3 - f) = -\Delta f + 3f^2\Delta f + 6f\|\nabla f\|^2 - \varepsilon^2\Delta^2 f, \qquad \varepsilon = 0.2$$

which can be written in terms of $O(2)$-equivariant terms

| Operator               | True coefficient         |
|:-----------------------|:-------------------------|
| $\Delta f$             | $-1.0$                   |
| $f^2\cdot\Delta f$     | $+3.0$                   |
| $f\cdot\|\nabla f\|^2$ | $+6.0$                   |
| $\Delta^2 f$           | $-\varepsilon^2 = -0.04$ |

or standard terms (PolyJet basis)

| Operator                  | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|:--------------------------|:-----------|:---------|:-------|:---------|:-----------|
| $\partial_{xx} f$         |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yy} f$         |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f^2\cdot\partial_{xx} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f^2\cdot\partial_{yy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f\cdot f_x^2$            |            | 🔴        | 🟢      | 🟢        | 🟢          |
| $f\cdot f_y^2$            |            | 🔴        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxxx} f$       |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxyy} f$       |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yyyy} f$       |            | 🟢        | 🟢      | 🟢        | 🟢          |

Cahn-Hilliard combines second and fourth-order terms with nonlinear products.  PDE-FIND cannot represent $f\cdot f_x^2$ (products of derivatives).  WSINDy cannot represent either $f^2\cdot f_{xx}$ or $f\cdot f_x^2$ individually, but the combination $3f^2\Delta f + 6f\|\nabla f\|^2$ equals $\Delta(f^3)$, which IS representable as $(f^3)_{xx} + (f^3)_{yy}$ in the WSINDy basis.

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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best             | Mg E_data | Std best       | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-----------------------|:----------|:---------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | 🔴        | 17/0 ✗            | 0.9209     | —                      | 0.6452    | WSINDy/STRidge | 0/0 ✓     | 0.1205      | 0.2837     | 2.9863       |
| Clean | pointwise     | polynomial | —        | 18/0 ✗            | 0.9994     | —                      | 0.7537    | —              | 2/0 ✗     | 0.8250      | 0.8429     | 39.4467      |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0136     | mod_multigraph/STRidge | 0.0535    | WSINDy/RFE     | 0/0 ✓     | 0.0162      | 0.0247     | 0.0608       |
| Clean | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.1935     | mod_multigraph/RFE     | 0.4138    | WSINDy/STRidge | 0/0 ✓     | 0.0111      | 0.0248     | 0.8057       |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0246     | mod_multigraph/RFE     | 0.0248    | WSINDy/RFE     | 0/0 ✓     | 0.0281      | 0.0247     | 0.0602       |
| Clean | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0111     | mod_multigraph/STRidge | 0.0248    | WSINDy/RFE     | 0/0 ✓     | 0.0121      | 0.0248     | 0.8047       |
| 2%    | pointwise     | spectral   | —        | 17/0 ✗            | 0.9775     | —                      | 0.8748    | —              | 18/0 ✗    | 0.9782      | 0.8863     | 216.5872     |
| 2%    | pointwise     | polynomial | —        | 18/0 ✗            | 0.9995     | —                      | 0.8168    | —              | 18/0 ✗    | 1.0458      | 0.9283     | 36.9446      |
| 2%    | weak (no IBP) | spectral   | 🔴        | 2/0 ✗             | 0.4197     | —                      | 0.4608    | WSINDy/STRidge | 0/0 ✓     | 0.0344      | 0.1886     | 0.4861       |
| 2%    | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.2288     | mod_multigraph/RFE     | 0.4642    | WSINDy/STRidge | 0/0 ✓     | 0.0301      | 0.1884     | 0.8072       |
| 2%    | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0361     | mod_multigraph/STRidge | 0.1886    | WSINDy/STRidge | 0/0 ✓     | 0.0412      | 0.1886     | 0.4863       |
| 2%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0301     | mod_multigraph/STRidge | 0.1886    | WSINDy/STRidge | 0/0 ✓     | 0.0308      | 0.1884     | 0.8061       |
| 5%    | pointwise     | spectral   | —        | 2/0 ✗             | 0.9979     | —                      | 0.9797    | —              | 37/0 ✗    | 0.9890      | 0.9651     | 404.9209     |
| 5%    | pointwise     | polynomial | —        | 22/0 ✗            | 0.9997     | —                      | 0.9213    | —              | 62/0 ✗    | 0.9999      | 0.9189     | 29.9851      |
| 5%    | weak (no IBP) | spectral   | —        | 10/0 ✗            | 1.6238     | —                      | 0.5544    | —              | 2/0 ✗     | 0.2005      | 0.4331     | 2.5795       |
| 5%    | weak (no IBP) | polynomial | 🟢        | 0/0 ✓             | 0.3896     | mod_multigraph/RFE     | 0.6366    | —              | 2/0 ✗     | 0.1894      | 0.4327     | 0.9214       |
| 5%    | weak (IBP)    | spectral   | 🟢        | 0/0 ✓             | 0.1417     | mod_multigraph/STRidge | 0.4384    | —              | 2/0 ✗     | 0.2022      | 0.4331     | 2.5797       |
| 5%    | weak (IBP)    | polynomial | 🟢        | 0/0 ✓             | 0.1271     | mod_multigraph/RFE     | 0.4384    | —              | 2/0 ✗     | 0.1897      | 0.4328     | 0.9199       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                      | fp/fn | E_coeff  | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:------|---------:|---------:|---------:|
| WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0111   | 2.1      | —        |
| Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0111   | 1.2      | +0.2%    |
| WSINDy + RFE                   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0113   | 0.9      | +2.1%    |
| Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0113   | 0.7      | +2.2%    |
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0113   | 1.4      | +2.2%    |
| Multigraphs + RFE              | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0349   | 0.7      | +215.9%  |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0349   | 1.4      | +215.9%  |
| Multigraphs + STRidge          | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0351   | 1.2      | +217.9%  |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0519   | 13.4     | +369.9%  |
| PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0519   | 0.9      | +369.9%  |
| PolyJet + STRidge              | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0525   | 4.6      | +374.8%  |
| WSINDy + MSTLS                 | spectral, weak (IBP), m=12      | 2/0 ✗ | 320.2403 | 2.9      | ✗        |
| PDE-FIND + STRidge             | spectral, weak (no IBP), m=4    | 5/2 ✗ | 1.1966   | 3.7      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=12   | 9/2 ✗ | 4.9996   | 8.0      | ✗        |
| PDE-FIND + RFE                 | polynomial, pointwise           | 2/6 ✗ | 0.9945   | 0.8      | ✗        |

### 2% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0301  | 1.8      | —        |
| Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓  | 0.0301  | 1.2      | +0.1%    |
| WSINDy + RFE                   | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0304  | 0.7      | +1.0%    |
| Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓  | 0.0305  | 0.7      | +1.2%    |
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=12    | 0/0 ✓  | 0.0491  | 1.4      | +63.3%   |
| Multigraphs + RFE              | spectral, weak (IBP), m=8       | 0/0 ✓  | 0.2910  | 0.7      | +866.8%  |
| PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/0 ✓  | 0.2987  | 0.7      | +892.2%  |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=8    | 2/0 ✗  | 0.4197  | 1.6      | ✗        |
| PolyJet + STRidge              | spectral, weak (no IBP), m=12   | 4/0 ✗  | 1.5168  | 4.3      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=8  | 7/0 ✗  | 0.4980  | 1.1      | ✗        |
| WSINDy + MSTLS                 | polynomial, weak (IBP), m=8     | 10/0 ✗ | 22.1400 | 3.6      | ✗        |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=8    | 11/0 ✗ | 1.9260  | 12.8     | ✗        |
| PDE-FIND + STRidge             | spectral, pointwise             | 5/2 ✗  | 0.9877  | 3.6      | ✗        |
| PDE-FIND + MSTLS               | polynomial, weak (no IBP), m=4  | 29/2 ✗ | 0.9901  | 6.3      | ✗        |
| PDE-FIND + RFE                 | spectral, pointwise             | 2/5 ✗  | 0.9774  | 0.7      | ✗        |

### 5% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓  | 0.1271  | 0.7      | —        |
| Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓  | 0.1276  | 1.2      | +0.4%    |
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=12    | 0/0 ✓  | 0.1755  | 1.5      | +38.1%   |
| Multigraphs + RFE              | polynomial, weak (IBP), m=8     | 0/0 ✓  | 0.5795  | 0.7      | +356.0%  |
| WSINDy + RFE                   | polynomial, weak (no IBP), m=12 | 2/0 ✗  | 0.1894  | 0.9      | ✗        |
| WSINDy + STRidge               | polynomial, weak (IBP), m=12    | 2/0 ✗  | 0.1916  | 2.0      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=8  | 7/0 ✗  | 0.7128  | 1.2      | ✗        |
| PolyJet + STRidge              | polynomial, weak (IBP), m=8     | 8/0 ✗  | 0.6992  | 4.6      | ✗        |
| PolyJet + RFE                  | spectral, weak (no IBP), m=12   | 8/0 ✗  | 1.5522  | 1.4      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=4    | 8/0 ✗  | 40.4156 | 3.2      | ✗        |
| Multigraphs + MSTLS            | polynomial, weak (IBP), m=8     | 9/0 ✗  | 0.6860  | 1.5      | ✗        |
| PolyJet + MSTLS                | polynomial, weak (IBP), m=12    | 10/0 ✗ | 0.5369  | 21.8     | ✗        |
| PDE-FIND + STRidge             | polynomial, weak (IBP), m=12    | 10/2 ✗ | 1.0172  | 2.9      | ✗        |
| PDE-FIND + MSTLS               | polynomial, weak (no IBP), m=4  | 31/2 ✗ | 0.9951  | 4.9      | ✗        |
| PDE-FIND + RFE                 | spectral, pointwise             | 3/4 ✗  | 0.9929  | 0.7      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 36/2  | 43       | 0.9731  | 0.7703 | 7.4      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 2/7   | 4        | 0.9769  | 0.7734 | 0.8      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.01    | 5/2   | 12       | 0.9877  | 0.8922 | 3.6      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 2/5   | 6        | 0.9774  | 0.9309 | 0.7      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.01    | 15/2  | 22       | 0.9948  | 0.9722 | 4.5      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 3/4   | 8        | 0.9929  | 0.9775 | 0.7      | 🔴      |
| WSINDy           | Clean | STRidge   | 1.0     | 0/0   | 7        | 0.1205  | 0.2837 | 2.0      | 🟢      |
| WSINDy           | Clean | RFE       |         | 0/0   | 7        | 0.1215  | 0.2837 | 0.9      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.01    | 18/0  | 12       | 0.9782  | 0.8863 | 2.9      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 2/5   | 6        | 0.9774  | 0.9309 | 1.2      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.001   | 96/0  | 46       | 0.9795  | 0.9654 | 3.9      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 6/3   | 8        | 0.9829  | 0.9762 | 0.7      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.1     | 15/0  | 24       | 0.9220  | 0.6098 | 5.5      | 🔴      |
| PolyJet          | Clean | RFE       |         | 2/2   | 9        | 0.9017  | 0.7424 | 1.0      | 🔴      |
| PolyJet          | 2%    | STRidge   | 0.01    | 78/0  | 87       | 0.9754  | 0.8639 | 13.1     | 🔴      |
| PolyJet          | 2%    | RFE       |         | 2/5   | 6        | 0.9774  | 0.9309 | 0.7      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.01    | 37/0  | 46       | 0.9890  | 0.9651 | 17.8     | 🔴      |
| PolyJet          | 5%    | RFE       |         | 3/4   | 8        | 0.9929  | 0.9775 | 1.4      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 17/0  | 11       | 0.9209  | 0.6452 | 1.4      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 2/2   | 5        | 0.9012  | 0.7430 | 0.7      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 0.01    | 17/0  | 11       | 0.9775  | 0.8748 | 1.7      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 7/7   | 5        | 0.9873  | 0.9138 | 0.7      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.01    | 17/0  | 11       | 0.9942  | 0.9697 | 1.9      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 5/2   | 6        | 0.9977  | 0.9766 | 0.9      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.1     | 18/0  | 11       | 0.9205  | 0.6496 | 1.5      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 2/2   | 4        | 0.9032  | 0.7431 | 0.7      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 0.01    | 18/0  | 11       | 0.9782  | 0.8751 | 1.7      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 18/6  | 5        | 0.9884  | 0.8956 | 0.9      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 18/0  | 9        | 0.9941  | 0.9702 | 1.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 2/0   | 5        | 0.9979  | 0.9797 | 0.7      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 19/2  | 26       | 0.9975  | 0.8344 | 4.3      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 2/6   | 5        | 0.9945  | 0.9004 | 0.8      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.001   | 43/2  | 50       | 0.9982  | 0.8857 | 6.6      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 2/6   | 5        | 0.9954  | 0.9266 | 0.9      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.01    | 17/2  | 24       | 0.9994  | 0.9638 | 3.8      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 2/6   | 5        | 0.9975  | 0.9744 | 0.7      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 2/0   | 9        | 0.8288  | 0.8430 | 2.2      | 🔴      |
| WSINDy           | Clean | RFE       |         | 2/0   | 9        | 0.8250  | 0.8429 | 1.0      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 18/0  | 12       | 1.0458  | 0.9283 | 2.9      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 2/8   | 3        | 0.9973  | 0.9756 | 1.2      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.001   | 83/0  | 42       | 1.0288  | 0.9763 | 4.5      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 2/7   | 4        | 0.9960  | 0.9813 | 0.9      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.01    | 35/0  | 44       | 0.9994  | 0.7258 | 8.8      | 🔴      |
| PolyJet          | Clean | RFE       |         | 4/7   | 6        | 0.9996  | 0.8642 | 0.8      | 🔴      |
| PolyJet          | 2%    | STRidge   | 0.1     | 31/0  | 40       | 0.9997  | 0.8027 | 6.2      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 3/9   | 3        | 0.9999  | 0.9211 | 0.8      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.01    | 62/0  | 71       | 0.9999  | 0.9189 | 14.3     | 🔴      |
| PolyJet          | 5%    | RFE       |         | 3/9   | 3        | 1.0000  | 0.9627 | 0.9      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.01    | 18/0  | 12       | 0.9994  | 0.7537 | 1.7      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 6/9   | 3        | 1.0036  | 0.8544 | 1.2      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 0.01    | 18/0  | 12       | 0.9995  | 0.8168 | 1.5      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 6/9   | 3        | 1.0034  | 0.8879 | 0.6      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 35/0  | 20       | 0.9997  | 0.9213 | 1.7      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 6/9   | 3        | 1.0027  | 0.9565 | 0.7      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.01    | 18/0  | 12       | 0.9994  | 0.7537 | 1.8      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 6/9   | 3        | 1.0036  | 0.8544 | 0.7      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 0.01    | 18/0  | 12       | 0.9995  | 0.8168 | 1.5      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 6/9   | 3        | 1.0034  | 0.8879 | 0.7      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.01    | 22/0  | 14       | 0.9997  | 0.9213 | 1.6      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 6/9   | 3        | 1.0027  | 0.9565 | 0.7      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff  | E_data   | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:---------|:---------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 4  | no  | 5/2   | 12       | 1.1966   | 0.7449   | 3.7      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.5     | 12 | no  | 9/2   | 16       | 4.9996   | 0.6064   | 8.0      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 8  | no  | 2/7   | 4        | 0.9654   | 0.7371   | 0.8      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 10.0    | 8  | yes | 5/2   | 12       | 3.8366   | 0.6752   | 2.9      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.5     | 12 | yes | 9/2   | 16       | 5.0016   | 0.6065   | 7.5      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 8  | yes | 2/7   | 4        | 0.9654   | 0.7371   | 4.5      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 8  | no  | 14/2  | 21       | 1.2058   | 0.7283   | 3.0      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.75    | 8  | no  | 31/2  | 38       | 1.2307   | 0.7230   | 4.9      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 8  | no  | 2/7   | 4        | 0.9658   | 0.7441   | 0.9      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 8  | yes | 14/2  | 21       | 1.2068   | 0.7284   | 3.0      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.75    | 8  | yes | 31/2  | 38       | 1.2317   | 0.7231   | 4.9      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 8  | yes | 2/7   | 4        | 0.9658   | 0.7441   | 0.7      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.1     | 8  | no  | 10/2  | 17       | 1.1240   | 0.7631   | 2.7      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 4  | no  | 2/5   | 6        | 1.0010   | 0.8351   | 6.2      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 2/7   | 4        | 0.9783   | 0.8424   | 0.9      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.1     | 8  | yes | 10/2  | 17       | 1.1240   | 0.7631   | 2.8      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.75    | 12 | yes | 8/5   | 12       | 1.0078   | 0.7191   | 5.7      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 2/7   | 4        | 0.9784   | 0.7680   | 0.7      | 🔴      |
| WSINDy           | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0163   | 0.0247   | 1.8      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 8  | no  | 2/0   | 9        | 704.2663 | 209.3462 | 3.2      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0162   | 0.0247   | 0.7      | 🟢      |
| WSINDy           | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0284   | 0.0247   | 1.8      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | yes | 2/0   | 9        | 320.2403 | 41.2621  | 2.9      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0281   | 0.0247   | 0.7      | 🟢      |
| WSINDy           | 2%    | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0344   | 0.1886   | 2.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 18/0  | 12       | 426.8456 | 40.0176  | 4.2      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0356   | 0.1886   | 1.2      | 🟢      |
| WSINDy           | 2%    | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0412   | 0.1886   | 2.1      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | yes | 18/0  | 12       | 409.8011 | 39.9805  | 4.2      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0432   | 0.1886   | 1.2      | 🟢      |
| WSINDy           | 5%    | STRidge   | 0.1     | 4  | no  | 2/0   | 9        | 1.0081   | 0.8300   | 2.1      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 4  | no  | 8/0   | 15       | 40.4156  | 37.5187  | 3.2      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | no  | 2/0   | 9        | 0.2005   | 0.4331   | 0.7      | 🔴      |
| WSINDy           | 5%    | STRidge   | 10.0    | 12 | yes | 2/0   | 9        | 0.2022   | 0.4331   | 1.8      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 8  | yes | 29/0  | 20       | 613.3703 | 94.6662  | 3.4      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | yes | 2/0   | 9        | 0.2043   | 0.4331   | 0.7      | 🔴      |
| PolyJet          | Clean | STRidge   | 1.0     | 8  | no  | 0/0   | 9        | 0.0525   | 0.0726   | 4.6      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 8  | no  | 0/0   | 9        | 0.0519   | 0.0726   | 13.4     | 🟢      |
| PolyJet          | Clean | RFE       |         | 8  | no  | 0/0   | 9        | 0.0519   | 0.0726   | 0.9      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 9        | 0.0553   | 0.0515   | 4.4      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 9        | 0.0559   | 0.0515   | 16.1     | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 9        | 0.0559   | 0.0515   | 1.5      | 🟢      |
| PolyJet          | 2%    | STRidge   | 10.0    | 12 | no  | 4/0   | 13       | 1.5168   | 0.2518   | 4.3      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 8  | no  | 11/0  | 20       | 1.9260   | 0.2965   | 12.8     | 🔴      |
| PolyJet          | 2%    | RFE       |         | 8  | no  | 0/0   | 9        | 0.2987   | 0.4864   | 0.7      | 🟢      |
| PolyJet          | 2%    | STRidge   | 10.0    | 12 | yes | 4/0   | 13       | 1.5282   | 0.2517   | 5.8      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 8  | yes | 11/0  | 20       | 1.9968   | 0.2977   | 13.4     | 🔴      |
| PolyJet          | 2%    | RFE       |         | 8  | yes | 0/0   | 9        | 0.3003   | 0.4869   | 0.8      | 🟢      |
| PolyJet          | 5%    | STRidge   | 10.0    | 8  | no  | 8/0   | 17       | 1.6224   | 0.5896   | 8.4      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | no  | 14/0  | 23       | 1.8560   | 0.4103   | 31.1     | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | no  | 8/0   | 17       | 1.5522   | 0.4480   | 1.4      | 🔴      |
| PolyJet          | 5%    | STRidge   | 10.0    | 8  | yes | 8/0   | 17       | 1.6554   | 0.5904   | 8.7      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | yes | 14/0  | 23       | 1.8653   | 0.4103   | 25.4     | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | yes | 8/0   | 17       | 1.5611   | 0.4480   | 1.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 4        | 0.0351   | 0.0734   | 1.2      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0349   | 0.0734   | 1.4      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 8  | no  | 0/0   | 4        | 0.0349   | 0.0734   | 0.7      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 4        | 0.0453   | 0.0522   | 1.2      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0458   | 0.0522   | 1.4      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0458   | 0.0522   | 0.7      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | no  | 7/0   | 6        | 1.6244   | 0.2144   | 1.3      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 8  | no  | 2/0   | 6        | 0.4197   | 0.4608   | 1.6      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 8  | no  | 2/0   | 6        | 0.4197   | 0.4608   | 0.7      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | yes | 7/0   | 6        | 1.6342   | 0.2145   | 1.4      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 8  | yes | 2/0   | 6        | 0.4515   | 0.4615   | 1.8      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 8  | yes | 0/0   | 4        | 0.2910   | 0.4873   | 0.7      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | no  | 10/0  | 7        | 1.6238   | 0.5544   | 1.6      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 12 | no  | 14/0  | 9        | 2.1135   | 0.5425   | 1.8      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | no  | 10/0  | 7        | 1.6248   | 0.6324   | 0.9      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | yes | 10/0  | 7        | 1.6283   | 0.5545   | 1.6      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 12 | yes | 18/0  | 12       | 1.8655   | 0.5236   | 1.8      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 10/0  | 7        | 1.6490   | 0.6328   | 0.9      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 3        | 0.0136   | 0.0535   | 1.4      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0140   | 0.0535   | 2.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0140   | 0.0535   | 1.5      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 3        | 0.0246   | 0.0248   | 1.3      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0246   | 0.0248   | 1.5      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0246   | 0.0248   | 0.8      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | no  | 6/0   | 7        | 1.8199   | 0.3166   | 1.5      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.5     | 8  | no  | 14/0  | 7        | 1.5510   | 0.3723   | 2.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 8  | no  | 2/0   | 5        | 0.4291   | 0.4633   | 0.9      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0361   | 0.1886   | 1.4      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.1809   | 0.1853   | 1.6      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0371   | 0.1886   | 0.8      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 18/0  | 6        | 1.5954   | 0.5747   | 1.3      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 4  | no  | 19/0  | 13       | 1.0154   | 0.8058   | 1.5      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | no  | 16/0  | 8        | 1.6448   | 0.6160   | 0.7      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.1417   | 0.4384   | 1.2      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 2/0   | 6        | 0.2215   | 0.4238   | 1.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.1428   | 0.4384   | 0.7      | 🟢      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data   | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:---------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | no  | 21/2  | 28       | 1.0248  | 0.6667   | 4.2      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.75    | 4  | no  | 29/2  | 36       | 0.9885  | 0.7547   | 4.9      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 8  | no  | 2/7   | 4        | 0.9725  | 0.7366   | 0.8      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | yes | 21/2  | 28       | 1.0243  | 0.6670   | 3.5      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.75    | 8  | yes | 38/2  | 45       | 0.9886  | 0.7123   | 5.7      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 8  | yes | 2/7   | 4        | 0.9725  | 0.7366   | 1.2      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 12 | no  | 23/2  | 30       | 1.0229  | 0.6740   | 4.7      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.75    | 4  | no  | 29/2  | 36       | 0.9901  | 0.7715   | 6.3      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 8  | no  | 2/7   | 4        | 0.9728  | 0.7439   | 1.2      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 12 | yes | 23/2  | 30       | 1.0225  | 0.6742   | 5.1      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.75    | 12 | yes | 22/3  | 28       | 1.0172  | 0.6703   | 7.2      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 8  | yes | 2/7   | 4        | 0.9728  | 0.7439   | 1.2      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 10/2  | 17       | 1.0175  | 0.7111   | 2.9      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.75    | 4  | no  | 31/2  | 38       | 0.9951  | 0.8291   | 4.9      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 8  | no  | 2/7   | 4        | 0.9798  | 0.7744   | 0.7      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 10/2  | 17       | 1.0172  | 0.7112   | 2.9      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.75    | 8  | yes | 35/2  | 42       | 0.9956  | 0.7527   | 4.9      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 2/7   | 4        | 0.9798  | 0.7742   | 0.8      | 🔴      |
| WSINDy           | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0111  | 0.0248   | 2.1      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | no  | 18/0  | 12       | 12.8940 | 41.4123  | 3.1      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0113  | 0.0248   | 0.9      | 🟢      |
| WSINDy           | Clean | STRidge   | 1.0     | 12 | yes | 10/0  | 10       | 0.0965  | 0.0230   | 2.1      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | yes | 10/0  | 10       | 13.2636 | 41.4248  | 3.0      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0121  | 0.0248   | 0.9      | 🟢      |
| WSINDy           | 2%    | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0301  | 0.1884   | 1.8      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 18/0  | 12       | 13.5319 | 40.0414  | 3.6      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0304  | 0.1884   | 0.7      | 🟢      |
| WSINDy           | 2%    | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0308  | 0.1884   | 2.0      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 8  | yes | 10/0  | 10       | 22.1400 | 176.3373 | 3.6      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0311  | 0.1884   | 0.7      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 4  | no  | 6/0   | 10       | 1.0178  | 0.8383   | 2.0      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 8  | no  | 25/0  | 18       | 29.3716 | 99.1766  | 3.7      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | no  | 2/0   | 9        | 0.1894  | 0.4327   | 0.9      | 🔴      |
| WSINDy           | 5%    | STRidge   | 10.0    | 12 | yes | 2/0   | 9        | 0.1916  | 0.4328   | 2.0      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 8  | yes | 29/0  | 20       | 21.4300 | 95.4550  | 3.9      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | yes | 2/0   | 9        | 0.1897  | 0.4328   | 0.9      | 🔴      |
| PolyJet          | Clean | STRidge   | 10.0    | 8  | no  | 7/0   | 16       | 0.4577  | 0.3432   | 3.9      | 🔴      |
| PolyJet          | Clean | MSTLS     | 0.25    | 8  | no  | 9/0   | 18       | 0.4899  | 0.3296   | 11.9     | 🔴      |
| PolyJet          | Clean | RFE       |         | 8  | no  | 0/0   | 9        | 0.3098  | 0.4598   | 0.8      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 9        | 0.2798  | 0.4152   | 4.2      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 8  | yes | 5/0   | 14       | 0.3045  | 0.3320   | 11.5     | 🔴      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 2/0   | 11       | 0.3115  | 0.3496   | 0.8      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 8  | no  | 9/0   | 18       | 0.4746  | 0.3910   | 5.2      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | no  | 14/0  | 23       | 0.4445  | 0.2441   | 13.8     | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/0   | 9        | 0.3007  | 0.4180   | 1.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 10.0    | 8  | yes | 9/0   | 18       | 0.4694  | 0.3528   | 5.3      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | yes | 14/0  | 23       | 0.4434  | 0.2437   | 13.9     | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/0   | 9        | 0.3010  | 0.4180   | 1.2      | 🟢      |
| PolyJet          | 5%    | STRidge   | 10.0    | 8  | no  | 8/0   | 17       | 0.7098  | 0.5704   | 4.6      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | no  | 10/0  | 19       | 0.5376  | 0.4459   | 23.6     | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | no  | 0/5   | 4        | 0.7467  | 0.7313   | 0.8      | 🔴      |
| PolyJet          | 5%    | STRidge   | 10.0    | 8  | yes | 8/0   | 17       | 0.6992  | 0.5682   | 4.6      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | yes | 10/0  | 19       | 0.5369  | 0.4458   | 21.8     | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | yes | 0/5   | 4        | 0.7467  | 0.7313   | 0.8      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | no  | 15/0  | 9        | 0.4528  | 0.1873   | 1.8      | 🔴      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.2643  | 0.3654   | 2.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.2643  | 0.3654   | 1.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 1.0     | 8  | yes | 15/0  | 11       | 0.4396  | 0.3049   | 1.8      | 🔴      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.2647  | 0.3654   | 2.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.2647  | 0.3654   | 1.2      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 8  | no  | 7/0   | 6        | 0.4980  | 0.4629   | 1.1      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 8  | no  | 2/0   | 6        | 0.4692  | 0.4872   | 1.5      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.2981  | 0.4204   | 0.6      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 8  | yes | 15/0  | 11       | 0.5105  | 0.4048   | 1.2      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | yes | 16/0  | 10       | 0.4334  | 0.2732   | 2.0      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.2985  | 0.4205   | 1.2      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 8  | no  | 7/0   | 6        | 0.7128  | 0.6195   | 1.2      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 8  | no  | 9/0   | 8        | 0.6970  | 0.6056   | 1.5      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | no  | 0/0   | 4        | 0.5894  | 0.7306   | 0.7      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 8  | yes | 17/0  | 11       | 0.7361  | 0.5861   | 1.2      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 8  | yes | 9/0   | 8        | 0.6860  | 0.6028   | 1.5      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 0/0   | 4        | 0.5795  | 0.7218   | 0.7      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 16/0  | 8        | 0.6206  | 0.3081   | 1.8      | 🔴      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 12 | no  | 2/0   | 6        | 0.3128  | 0.3514   | 2.2      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.1935  | 0.4138   | 1.2      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 3        | 0.0111  | 0.0248   | 1.2      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0113  | 0.0248   | 1.4      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0113  | 0.0248   | 0.7      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 12 | no  | 15/0  | 11       | 0.4199  | 0.3167   | 1.3      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.2981  | 0.4204   | 1.6      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 3        | 0.2288  | 0.4642   | 0.7      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0301  | 0.1886   | 1.2      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0491  | 0.1857   | 1.4      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0305  | 0.1886   | 0.7      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 8  | no  | 7/0   | 5        | 0.7215  | 0.6209   | 1.2      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.5     | 8  | no  | 9/0   | 8        | 0.6970  | 0.6056   | 1.5      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | no  | 0/0   | 3        | 0.3896  | 0.6366   | 0.7      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.1276  | 0.4384   | 1.2      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.1755  | 0.4218   | 1.5      | 🟢      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.1271  | 0.4384   | 0.7      | 🟢      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ) |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:--------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.2604       | 6829.7  |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0839       | 8665.4  |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0608       | 13607.8 |
| Clean | weak (no IBP) | polynomial | 4  | no  | 7.8874       | 3660.2  |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.9336       | 4661.8  |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.8057       | 7174.7  |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0897       | 8662.7  |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0602       | 13605.4 |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.9050       | 4719.8  |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.8047       | 7167.9  |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 3.4975       | 5305.4  |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 0.5875       | 7937.5  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 0.4861       | 13214.6 |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 7.9794       | 3555.6  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.9941       | 4466.0  |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.8072       | 6899.8  |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 0.5886       | 7958.9  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 0.4863       | 13218.4 |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 0.9691       | 4506.8  |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.8061       | 6897.0  |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 8.9563       | 8724.2  |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 2.7055       | 15618.3 |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 2.5795       | 25685.6 |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 8.4374       | 3198.1  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 1.2927       | 4026.9  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 0.9214       | 6294.1  |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 2.7062       | 15643.2 |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 2.5797       | 25692.9 |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 1.2794       | 4047.7  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 0.9199       | 6298.3  |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ch`.

| Source                                 | Wall time | CPU | Threads |
|:---------------------------------------|----------:|:----|:-------:|
| `sweep_ch_full_ch_20260803_181412.csv` | ? min     | ?   | ?       |
