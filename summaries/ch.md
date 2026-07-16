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
| $f^2\cdot\partial_{xx} f$ |            | 🟢        | 🔴      | 🟢        | 🟢          |
| $f^2\cdot\partial_{yy} f$ |            | 🟢        | 🔴      | 🟢        | 🟢          |
| $f\cdot f_x^2$            |            | 🔴        | 🔴      | 🟢        | 🟢          |
| $f\cdot f_y^2$            |            | 🔴        | 🔴      | 🟢        | 🟢          |
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
| Clean | pointwise     | polynomial | —        | 18/0 ✗            | 0.9986     | —                      | 0.7526    | —              | 2/0 ✗     | 0.8250      | 0.8429     | 39.4467      |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0136     | mod_multigraph/STRidge | 0.0535    | WSINDy/RFE     | 0/0 ✓     | 0.0162      | 0.0247     | 0.0608       |
| Clean | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.1935     | mod_multigraph/RFE     | 0.4138    | WSINDy/STRidge | 0/0 ✓     | 0.0111      | 0.0248     | 0.8057       |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0238     | mod_multigraph/STRidge | 0.0227    | WSINDy/RFE     | 0/0 ✓     | 0.0281      | 0.0247     | 0.0602       |
| Clean | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0111     | mod_multigraph/STRidge | 0.0248    | WSINDy/RFE     | 0/0 ✓     | 0.0121      | 0.0248     | 0.8047       |
| 2%    | pointwise     | spectral   | —        | 35/0 ✗            | 1.0000     | —                      | 0.8748    | —              | 18/0 ✗    | 0.9999      | 0.8863     | 216.5872     |
| 2%    | pointwise     | polynomial | —        | 18/0 ✗            | 0.9985     | —                      | 0.8164    | —              | 18/0 ✗    | 1.0441      | 0.9283     | 36.9446      |
| 2%    | weak (no IBP) | spectral   | 🔴        | 2/0 ✗             | 0.3742     | —                      | 0.4633    | WSINDy/STRidge | 0/0 ✓     | 0.0296      | 0.1886     | 0.4861       |
| 2%    | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.2307     | mod_multigraph/RFE     | 0.4642    | WSINDy/STRidge | 0/0 ✓     | 0.0301      | 0.1884     | 0.8072       |
| 2%    | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0292     | mod_multigraph/RFE     | 0.1886    | WSINDy/STRidge | 0/0 ✓     | 0.0292      | 0.1886     | 0.4863       |
| 2%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0305     | mod_multigraph/RFE     | 0.1886    | WSINDy/STRidge | 0/0 ✓     | 0.0307      | 0.1884     | 0.8061       |
| 5%    | pointwise     | spectral   | —        | 2/0 ✗             | 1.0000     | —                      | 0.9797    | —              | 10/0 ✗    | 1.0000      | 0.9669     | 404.9209     |
| 5%    | pointwise     | polynomial | —        | 18/0 ✗            | 0.9985     | —                      | 0.9218    | —              | 77/0 ✗    | 1.0241      | 0.9763     | 29.9851      |
| 5%    | weak (no IBP) | spectral   | —        | 10/0 ✗            | 0.6939     | —                      | 0.6324    | —              | 2/0 ✗     | 0.1892      | 0.4331     | 2.5795       |
| 5%    | weak (no IBP) | polynomial | 🟢        | 0/0 ✓             | 0.5804     | mod_multigraph/RFE     | 0.7523    | —              | 2/0 ✗     | 0.1888      | 0.4327     | 0.9214       |
| 5%    | weak (IBP)    | spectral   | 🟢        | 0/0 ✓             | 0.1234     | mod_multigraph/RFE     | 0.4384    | —              | 2/0 ✗     | 0.1885      | 0.4331     | 2.5797       |
| 5%    | weak (IBP)    | polynomial | 🟢        | 0/0 ✓             | 0.1268     | mod_multigraph/RFE     | 0.4384    | —              | 2/0 ✗     | 0.1888      | 0.4328     | 0.9199       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                      | fp/fn | E_coeff  | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:------|---------:|---------:|---------:|
| WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0111   | 1.2      | —        |
| Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0111   | 1.2      | +0.2%    |
| WSINDy + RFE                   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0113   | 0.9      | +2.1%    |
| Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0113   | 1.1      | +2.2%    |
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0113   | 1.6      | +2.2%    |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0349   | 1.3      | +215.9%  |
| Multigraphs + RFE              | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0349   | 0.9      | +215.9%  |
| PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0519   | 1.3      | +369.9%  |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0519   | 6.4      | +369.9%  |
| PolyJet + STRidge              | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0525   | 2.0      | +374.8%  |
| Multigraphs + STRidge          | spectral, weak (no IBP), m=4    | 0/0 ✓ | 0.1207   | 1.0      | +991.5%  |
| WSINDy + MSTLS                 | spectral, weak (IBP), m=12      | 2/0 ✗ | 320.2403 | 2.4      | ✗        |
| PDE-FIND + STRidge             | polynomial, weak (no IBP), m=4  | 5/2 ✗ | 0.9882   | 1.4      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=12   | 7/2 ✗ | 5.0162   | 4.2      | ✗        |
| PDE-FIND + RFE                 | polynomial, pointwise           | 2/6 ✗ | 0.9945   | 1.0      | ✗        |

### 2% noise

| Method                         | Best conf.                     | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:-------------------------------|:-------|--------:|---------:|---------:|
| Modified Multigraphs + RFE     | spectral, weak (IBP), m=12     | 0/0 ✓  | 0.0292  | 1.0      | —        |
| WSINDy + STRidge               | spectral, weak (IBP), m=12     | 0/0 ✓  | 0.0292  | 1.2      | +0.1%    |
| WSINDy + RFE                   | spectral, weak (IBP), m=12     | 0/0 ✓  | 0.0295  | 0.9      | +1.3%    |
| Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12   | 0/0 ✓  | 0.0493  | 1.3      | +68.8%   |
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=12   | 0/0 ✓  | 0.0493  | 1.6      | +69.1%   |
| Multigraphs + RFE              | spectral, weak (IBP), m=8      | 0/0 ✓  | 0.2964  | 0.9      | +916.0%  |
| PolyJet + RFE                  | spectral, weak (IBP), m=8      | 0/0 ✓  | 0.2986  | 1.3      | +923.5%  |
| Multigraphs + MSTLS            | spectral, weak (IBP), m=8      | 2/0 ✗  | 0.3701  | 1.5      | ✗        |
| PolyJet + STRidge              | spectral, weak (no IBP), m=8   | 5/0 ✗  | 0.4680  | 2.1      | ✗        |
| WSINDy + MSTLS                 | polynomial, weak (IBP), m=8    | 10/0 ✗ | 22.7095 | 2.7      | ✗        |
| PolyJet + MSTLS                | spectral, weak (IBP), m=12     | 11/0 ✗ | 0.1754  | 6.5      | ✗        |
| Multigraphs + STRidge          | spectral, weak (no IBP), m=4   | 11/0 ✗ | 0.9712  | 1.0      | ✗        |
| PDE-FIND + STRidge             | spectral, pointwise            | 5/2 ✗  | 0.9999  | 1.4      | ✗        |
| PDE-FIND + MSTLS               | polynomial, weak (no IBP), m=4 | 21/2 ✗ | 0.9894  | 3.7      | ✗        |
| PDE-FIND + RFE                 | spectral, pointwise            | 2/5 ✗  | 0.9997  | 1.0      | ✗        |

### 5% noise

| Method                         | Best conf.                     | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:-------------------------------|:-------|--------:|---------:|---------:|
| Modified Multigraphs + RFE     | spectral, weak (IBP), m=12     | 0/0 ✓  | 0.1234  | 0.9      | —        |
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=12   | 0/0 ✓  | 0.1795  | 1.4      | +45.5%   |
| Modified Multigraphs + STRidge | spectral, weak (IBP), m=12     | 0/0 ✓  | 0.1965  | 1.0      | +59.3%   |
| Multigraphs + RFE              | polynomial, weak (IBP), m=8    | 0/0 ✓  | 0.5848  | 0.9      | +374.1%  |
| WSINDy + RFE                   | spectral, weak (IBP), m=12     | 2/0 ✗  | 0.1885  | 1.0      | ✗        |
| WSINDy + STRidge               | spectral, weak (IBP), m=12     | 2/0 ✗  | 0.1908  | 1.3      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (IBP), m=4    | 7/0 ✗  | 1.0024  | 1.0      | ✗        |
| PolyJet + RFE                  | spectral, weak (IBP), m=12     | 8/0 ✗  | 0.3244  | 1.3      | ✗        |
| PolyJet + STRidge              | spectral, weak (IBP), m=8      | 8/0 ✗  | 0.6135  | 2.6      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=4   | 8/0 ✗  | 3.2675  | 3.3      | ✗        |
| PolyJet + MSTLS                | polynomial, weak (IBP), m=12   | 10/0 ✗ | 0.5318  | 8.0      | ✗        |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=12  | 14/0 ✗ | 0.5736  | 1.5      | ✗        |
| PDE-FIND + STRidge             | spectral, pointwise            | 3/2 ✗  | 1.0000  | 1.6      | ✗        |
| PDE-FIND + MSTLS               | polynomial, weak (no IBP), m=4 | 23/2 ✗ | 0.9926  | 4.1      | ✗        |
| PDE-FIND + RFE                 | spectral, pointwise            | 3/4 ✗  | 1.0000  | 1.1      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm   | d_tol/α   | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:------------|:----------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge     | 0.001     | 36/2  | 57       | 0.9737  | 0.7702 | 1.6      | 🔴      |
| PDE-FIND         | Clean | RFE         |           | 2/7   | 4        | 0.9769  | 0.7734 | 1.0      | 🔴      |
| PDE-FIND         | 2%    | STRidge     | 0.1       | 5/2   | 12       | 0.9999  | 0.8922 | 1.4      | 🔴      |
| PDE-FIND         | 2%    | RFE         |           | 2/5   | 6        | 0.9997  | 0.9309 | 1.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge     | 0.1       | 3/2   | 10       | 1.0000  | 0.9763 | 1.6      | 🔴      |
| PDE-FIND         | 5%    | RFE         |           | 3/4   | 8        | 1.0000  | 0.9775 | 1.1      | 🔴      |
| WSINDy           | Clean | STRidge     | 10.0      | 0/0   | 7        | 0.1205  | 0.2837 | 1.2      | 🟢      |
| WSINDy           | Clean | RFE         |           | 0/0   | 7        | 0.1215  | 0.2837 | 0.9      | 🟢      |
| WSINDy           | 2%    | STRidge     | 0.1       | 18/0  | 12       | 0.9999  | 0.8863 | 1.3      | 🔴      |
| WSINDy           | 2%    | RFE         |           | 2/5   | 6        | 0.9997  | 0.9309 | 0.9      | 🔴      |
| WSINDy           | 5%    | STRidge     | 0.1       | 18/0  | 12       | 1.0000  | 0.9658 | 1.4      | 🔴      |
| WSINDy           | 5%    | RFE         |           | 6/3   | 8        | 1.0000  | 0.9762 | 1.1      | 🔴      |
| PolyJet          | Clean | STRidge     | 1.0       | 4/0   | 13       | 0.9071  | 0.8115 | 2.1      | 🔴      |
| PolyJet          | Clean | RFE         |           | 2/2   | 9        | 0.9017  | 0.7424 | 1.3      | 🔴      |
| PolyJet          | 2%    | STRidge     | 0.001     | 77/0  | 86       | 1.0000  | 0.8642 | 3.0      | 🔴      |
| PolyJet          | 2%    | RFE         |           | 2/5   | 6        | 0.9997  | 0.9309 | 1.3      | 🔴      |
| PolyJet          | 5%    | STRidge     | 0.1       | 10/0  | 19       | 1.0000  | 0.9669 | 2.4      | 🔴      |
| PolyJet          | 5%    | RFE         |           | 3/4   | 8        | 1.0000  | 0.9775 | 1.6      | 🔴      |
| Multigraphs      | Clean | STRidge     | 1.0       | 17/0  | 11       | 0.9209  | 0.6452 | 1.0      | 🔴      |
| Multigraphs      | Clean | RFE         |           | 2/2   | 5        | 0.9012  | 0.7430 | 0.9      | 🔴      |
| Multigraphs      | Clean | Brute-force | L₀=0.01   | 17/0  | 11       | 0.9207  | 0.6450 | 11.9     | 🔴      |
| Multigraphs      | 2%    | STRidge     | 0.001     | 35/0  | 20       | 1.0000  | 0.8748 | 1.1      | 🔴      |
| Multigraphs      | 2%    | RFE         |           | 7/7   | 5        | 0.9996  | 0.9138 | 0.9      | 🔴      |
| Multigraphs      | 2%    | Brute-force | L₀=0.01   | 11/2  | 8        | 1.0000  | 0.8789 | 11.3     | 🔴      |
| Multigraphs      | 5%    | STRidge     | 0.1       | 14/0  | 10       | 1.0000  | 0.9700 | 1.3      | 🔴      |
| Multigraphs      | 5%    | RFE         |           | 5/2   | 6        | 1.0000  | 0.9766 | 1.0      | 🔴      |
| Multigraphs      | 5%    | Brute-force | L₀=0.0001 | 21/0  | 12       | 1.0000  | 0.9696 | 14.0     | 🔴      |
| Mod. Multigraphs | Clean | STRidge     | 1.0       | 18/0  | 10       | 0.9204  | 0.6533 | 1.0      | 🔴      |
| Mod. Multigraphs | Clean | RFE         |           | 2/2   | 4        | 0.9032  | 0.7431 | 0.9      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge     | 0.001     | 35/0  | 20       | 1.0000  | 0.8748 | 1.1      | 🔴      |
| Mod. Multigraphs | 2%    | RFE         |           | 18/6  | 5        | 1.0000  | 0.8956 | 0.9      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge     | 0.1       | 18/0  | 9        | 1.0000  | 0.9702 | 1.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE         |           | 2/0   | 5        | 1.0000  | 0.9797 | 0.9      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm   | d_tol/α   | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:------------|:----------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge     | 0.1       | 5/2   | 12       | 0.9975  | 0.8347 | 1.4      | 🔴      |
| PDE-FIND         | Clean | RFE         |           | 2/6   | 5        | 0.9945  | 0.9004 | 1.0      | 🔴      |
| PDE-FIND         | 2%    | STRidge     | 0.001     | 39/2  | 60       | 0.9978  | 0.8857 | 1.7      | 🔴      |
| PDE-FIND         | 2%    | RFE         |           | 2/6   | 5        | 0.9951  | 0.9266 | 1.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge     | 0.1       | 5/2   | 14       | 0.9988  | 0.9640 | 1.7      | 🔴      |
| PDE-FIND         | 5%    | RFE         |           | 2/6   | 5        | 0.9968  | 0.9744 | 1.0      | 🔴      |
| WSINDy           | Clean | STRidge     | 1.0       | 2/0   | 9        | 0.8288  | 0.8430 | 1.2      | 🔴      |
| WSINDy           | Clean | RFE         |           | 2/0   | 9        | 0.8250  | 0.8429 | 0.9      | 🔴      |
| WSINDy           | 2%    | STRidge     | 0.1       | 18/0  | 12       | 1.0441  | 0.9283 | 1.3      | 🔴      |
| WSINDy           | 2%    | RFE         |           | 2/8   | 3        | 0.9973  | 0.9756 | 0.9      | 🔴      |
| WSINDy           | 5%    | STRidge     | 0.01      | 77/0  | 39       | 1.0241  | 0.9763 | 1.6      | 🔴      |
| WSINDy           | 5%    | RFE         |           | 2/7   | 4        | 0.9955  | 0.9813 | 1.0      | 🔴      |
| PolyJet          | Clean | STRidge     | 0.1       | 15/0  | 24       | 0.9994  | 0.7276 | 2.4      | 🔴      |
| PolyJet          | Clean | RFE         |           | 4/7   | 6        | 0.9996  | 0.8642 | 1.3      | 🔴      |
| PolyJet          | 2%    | STRidge     | 0.1       | 31/0  | 40       | 0.9994  | 0.8027 | 2.4      | 🔴      |
| PolyJet          | 2%    | RFE         |           | 3/9   | 3        | 0.9998  | 0.9211 | 1.3      | 🔴      |
| PolyJet          | 5%    | STRidge     | 0.01      | 80/0  | 89       | 0.9994  | 0.9187 | 3.4      | 🔴      |
| PolyJet          | 5%    | RFE         |           | 3/9   | 3        | 0.9995  | 0.9627 | 1.4      | 🔴      |
| Multigraphs      | Clean | STRidge     | 0.1       | 18/0  | 12       | 0.9986  | 0.7526 | 1.0      | 🔴      |
| Multigraphs      | Clean | RFE         |           | 6/9   | 3        | 1.0036  | 0.8544 | 0.9      | 🔴      |
| Multigraphs      | Clean | Brute-force | L₀=0.0001 | 18/2  | 11       | 0.9986  | 0.7525 | 11.9     | 🔴      |
| Multigraphs      | 2%    | STRidge     | 0.1       | 18/0  | 12       | 0.9985  | 0.8164 | 1.0      | 🔴      |
| Multigraphs      | 2%    | RFE         |           | 6/9   | 3        | 1.0030  | 0.8879 | 0.9      | 🔴      |
| Multigraphs      | 2%    | Brute-force | L₀=0.0001 | 18/0  | 12       | 0.9986  | 0.8163 | 11.7     | 🔴      |
| Multigraphs      | 5%    | STRidge     | 0.1       | 18/0  | 12       | 0.9985  | 0.9218 | 1.2      | 🔴      |
| Multigraphs      | 5%    | RFE         |           | 6/9   | 3        | 1.0013  | 0.9565 | 1.1      | 🔴      |
| Multigraphs      | 5%    | Brute-force | L₀=0.0001 | 18/0  | 12       | 0.9985  | 0.9218 | 10.9     | 🔴      |
| Mod. Multigraphs | Clean | STRidge     | 0.1       | 18/0  | 12       | 0.9986  | 0.7526 | 1.1      | 🔴      |
| Mod. Multigraphs | Clean | RFE         |           | 6/9   | 3        | 1.0036  | 0.8544 | 1.0      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge     | 0.1       | 18/0  | 12       | 0.9985  | 0.8164 | 1.1      | 🔴      |
| Mod. Multigraphs | 2%    | RFE         |           | 6/9   | 3        | 1.0030  | 0.8879 | 1.0      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge     | 0.1       | 18/0  | 12       | 0.9985  | 0.9218 | 1.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE         |           | 6/9   | 3        | 1.0013  | 0.9565 | 1.0      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm   | d_tol/α | m    | IBP | fp/fn | n_active | E_coeff  | E_data   | Time (s) | Exact? |
|:-----------------|:------|:------------|:--------|:-----|:----|:------|:---------|:---------|:---------|:---------|:-------|
| PDE-FIND         | Clean | STRidge     | 10.0    | 8.0  | no  | 5/2   | 12       | 3.8398   | 0.6742   | 1.3      | 🔴      |
| PDE-FIND         | Clean | MSTLS       | 0.5     | 12.0 | no  | 7/2   | 16       | 5.0162   | 0.6064   | 4.2      | 🔴      |
| PDE-FIND         | Clean | RFE         |         | 8.0  | no  | 2/7   | 4        | 0.9654   | 0.7371   | 1.0      | 🔴      |
| PDE-FIND         | Clean | STRidge     | 1.0     | 4.0  | yes | 5/2   | 12       | 1.0947   | 0.7482   | 1.3      | 🔴      |
| PDE-FIND         | Clean | MSTLS       | 0.5     | 12.0 | yes | 7/2   | 16       | 5.0183   | 0.6065   | 4.2      | 🔴      |
| PDE-FIND         | Clean | RFE         |         | 8.0  | yes | 2/7   | 4        | 0.9654   | 0.7371   | 1.0      | 🔴      |
| PDE-FIND         | 2%    | STRidge     | 1.0     | 8.0  | no  | 6/2   | 14       | 0.9622   | 0.7296   | 1.4      | 🔴      |
| PDE-FIND         | 2%    | MSTLS       | 0.75    | 8.0  | no  | 22/2  | 38       | 0.9609   | 0.7230   | 3.5      | 🔴      |
| PDE-FIND         | 2%    | RFE         |         | 4.0  | no  | 2/7   | 4        | 0.9935   | 0.7728   | 1.0      | 🔴      |
| PDE-FIND         | 2%    | STRidge     | 1.0     | 8.0  | yes | 6/2   | 14       | 0.9621   | 0.7297   | 1.4      | 🔴      |
| PDE-FIND         | 2%    | MSTLS       | 0.75    | 8.0  | yes | 22/2  | 38       | 0.9607   | 0.7231   | 3.5      | 🔴      |
| PDE-FIND         | 2%    | RFE         |         | 4.0  | yes | 2/7   | 4        | 0.9935   | 0.7731   | 1.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge     | 1.0     | 8.0  | no  | 12/2  | 22       | 0.9711   | 0.7605   | 1.7      | 🔴      |
| PDE-FIND         | 5%    | MSTLS       | 0.25    | 4.0  | no  | 2/5   | 6        | 0.9956   | 0.8351   | 4.0      | 🔴      |
| PDE-FIND         | 5%    | RFE         |         | 4.0  | no  | 2/7   | 4        | 0.9949   | 0.8424   | 1.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge     | 1.0     | 8.0  | yes | 12/2  | 22       | 0.9711   | 0.7605   | 1.6      | 🔴      |
| PDE-FIND         | 5%    | MSTLS       | 0.25    | 4.0  | yes | 2/5   | 6        | 0.9956   | 0.8350   | 3.7      | 🔴      |
| PDE-FIND         | 5%    | RFE         |         | 4.0  | yes | 2/7   | 4        | 0.9950   | 0.8424   | 1.0      | 🔴      |
| WSINDy           | Clean | STRidge     | 10.0    | 12.0 | no  | 0/0   | 7        | 0.0163   | 0.0247   | 1.2      | 🟢      |
| WSINDy           | Clean | MSTLS       | 0.25    | 8.0  | no  | 2/0   | 9        | 704.2663 | 209.3462 | 2.5      | 🔴      |
| WSINDy           | Clean | RFE         |         | 12.0 | no  | 0/0   | 7        | 0.0162   | 0.0247   | 0.9      | 🟢      |
| WSINDy           | Clean | STRidge     | 10.0    | 12.0 | yes | 0/0   | 7        | 0.0284   | 0.0247   | 1.2      | 🟢      |
| WSINDy           | Clean | MSTLS       | 0.25    | 12.0 | yes | 2/0   | 9        | 320.2403 | 41.2621  | 2.4      | 🔴      |
| WSINDy           | Clean | RFE         |         | 12.0 | yes | 0/0   | 7        | 0.0281   | 0.0247   | 0.9      | 🟢      |
| WSINDy           | 2%    | STRidge     | 10.0    | 12.0 | no  | 0/0   | 7        | 0.0296   | 0.1886   | 1.3      | 🟢      |
| WSINDy           | 2%    | MSTLS       | 0.25    | 12.0 | no  | 18/0  | 12       | 22.5244  | 40.0176  | 2.8      | 🔴      |
| WSINDy           | 2%    | RFE         |         | 12.0 | no  | 0/0   | 7        | 0.0299   | 0.1886   | 0.9      | 🟢      |
| WSINDy           | 2%    | STRidge     | 10.0    | 12.0 | yes | 0/0   | 7        | 0.0292   | 0.1886   | 1.2      | 🟢      |
| WSINDy           | 2%    | MSTLS       | 0.25    | 12.0 | yes | 18/0  | 12       | 22.6549  | 39.9805  | 2.8      | 🔴      |
| WSINDy           | 2%    | RFE         |         | 12.0 | yes | 0/0   | 7        | 0.0295   | 0.1886   | 0.9      | 🟢      |
| WSINDy           | 5%    | STRidge     | 10.0    | 8.0  | no  | 6/0   | 10       | 0.6603   | 0.6360   | 1.4      | 🔴      |
| WSINDy           | 5%    | MSTLS       | 0.25    | 4.0  | no  | 8/0   | 15       | 3.2675   | 37.5187  | 3.3      | 🔴      |
| WSINDy           | 5%    | RFE         |         | 12.0 | no  | 2/0   | 9        | 0.1892   | 0.4331   | 0.9      | 🔴      |
| WSINDy           | 5%    | STRidge     | 10.0    | 12.0 | yes | 2/0   | 9        | 0.1908   | 0.4331   | 1.3      | 🔴      |
| WSINDy           | 5%    | MSTLS       | 0.25    | 4.0  | yes | 16/0  | 17       | 4.6420   | 71.2490  | 2.9      | 🔴      |
| WSINDy           | 5%    | RFE         |         | 12.0 | yes | 2/0   | 9        | 0.1885   | 0.4331   | 1.0      | 🔴      |
| PolyJet          | Clean | STRidge     | 10.0    | 8.0  | no  | 0/0   | 9        | 0.0525   | 0.0726   | 2.0      | 🟢      |
| PolyJet          | Clean | MSTLS       | 0.25    | 8.0  | no  | 0/0   | 9        | 0.0519   | 0.0726   | 6.4      | 🟢      |
| PolyJet          | Clean | RFE         |         | 8.0  | no  | 0/0   | 9        | 0.0519   | 0.0726   | 1.3      | 🟢      |
| PolyJet          | Clean | STRidge     | 10.0    | 12.0 | yes | 0/0   | 9        | 0.0553   | 0.0515   | 2.0      | 🟢      |
| PolyJet          | Clean | MSTLS       | 0.25    | 12.0 | yes | 0/0   | 9        | 0.0559   | 0.0515   | 7.4      | 🟢      |
| PolyJet          | Clean | RFE         |         | 12.0 | yes | 0/0   | 9        | 0.0559   | 0.0515   | 1.3      | 🟢      |
| PolyJet          | 2%    | STRidge     | 10.0    | 8.0  | no  | 5/0   | 14       | 0.4680   | 0.3914   | 2.1      | 🔴      |
| PolyJet          | 2%    | MSTLS       | 0.25    | 12.0 | no  | 11/0  | 20       | 0.1761   | 0.1773   | 6.7      | 🔴      |
| PolyJet          | 2%    | RFE         |         | 8.0  | no  | 0/0   | 9        | 0.2988   | 0.4864   | 1.3      | 🟢      |
| PolyJet          | 2%    | STRidge     | 10.0    | 8.0  | yes | 5/0   | 14       | 0.4682   | 0.3943   | 2.1      | 🔴      |
| PolyJet          | 2%    | MSTLS       | 0.25    | 12.0 | yes | 11/0  | 20       | 0.1754   | 0.1773   | 6.5      | 🔴      |
| PolyJet          | 2%    | RFE         |         | 8.0  | yes | 0/0   | 9        | 0.2986   | 0.4869   | 1.3      | 🟢      |
| PolyJet          | 5%    | STRidge     | 10.0    | 8.0  | no  | 8/0   | 17       | 0.6155   | 0.5896   | 2.3      | 🔴      |
| PolyJet          | 5%    | MSTLS       | 0.25    | 12.0 | no  | 14/0  | 23       | 0.3016   | 0.4103   | 7.8      | 🔴      |
| PolyJet          | 5%    | RFE         |         | 12.0 | no  | 8/0   | 17       | 0.3250   | 0.4480   | 1.6      | 🔴      |
| PolyJet          | 5%    | STRidge     | 10.0    | 8.0  | yes | 8/0   | 17       | 0.6135   | 0.5904   | 2.6      | 🔴      |
| PolyJet          | 5%    | MSTLS       | 0.25    | 12.0 | yes | 14/0  | 23       | 0.3010   | 0.4103   | 7.2      | 🔴      |
| PolyJet          | 5%    | RFE         |         | 12.0 | yes | 8/0   | 17       | 0.3244   | 0.4480   | 1.3      | 🔴      |
| Multigraphs      | Clean | STRidge     | 10.0    | 4.0  | no  | 0/0   | 4        | 0.1207   | 0.2305   | 1.0      | 🟢      |
| Multigraphs      | Clean | MSTLS       | 0.25    | 8.0  | no  | 0/0   | 4        | 0.0349   | 0.0734   | 1.3      | 🟢      |
| Multigraphs      | Clean | RFE         |         | 8.0  | no  | 0/0   | 4        | 0.0349   | 0.0734   | 0.9      | 🟢      |
| Multigraphs      | Clean | Brute-force | L₀=0.01 | 8.0  | no  | 0/0   | 4        | 0.0349   | 0.0734   | 11.5     | 🟢      |
| Multigraphs      | Clean | STRidge     | 10.0    | 8.0  | yes | 13/0  | 8        | 0.9062   | 0.0662   | 1.0      | 🔴      |
| Multigraphs      | Clean | MSTLS       | 0.25    | 12.0 | yes | 0/0   | 4        | 0.0458   | 0.0522   | 1.3      | 🟢      |
| Multigraphs      | Clean | RFE         |         | 12.0 | yes | 0/0   | 4        | 0.0458   | 0.0522   | 0.9      | 🟢      |
| Multigraphs      | Clean | Brute-force | L₀=0.01 | 12.0 | yes | 0/0   | 4        | 0.0458   | 0.0522   | 11.2     | 🟢      |
| Multigraphs      | 2%    | STRidge     | 1.0     | 4.0  | no  | 11/0  | 10       | 0.9712   | 0.7050   | 1.0      | 🔴      |
| Multigraphs      | 2%    | MSTLS       | 0.25    | 8.0  | no  | 2/0   | 6        | 0.3744   | 0.4608   | 1.4      | 🔴      |
| Multigraphs      | 2%    | RFE         |         | 8.0  | no  | 2/0   | 6        | 0.3744   | 0.4608   | 0.9      | 🔴      |
| Multigraphs      | 2%    | Brute-force | L₀=0.01 | 12.0 | no  | 7/0   | 6        | 0.1431   | 0.2144   | 10.5     | 🔴      |
| Multigraphs      | 2%    | STRidge     | 10.0    | 8.0  | yes | 17/0  | 11       | 0.3393   | 0.3020   | 1.0      | 🔴      |
| Multigraphs      | 2%    | MSTLS       | 0.25    | 8.0  | yes | 2/0   | 6        | 0.3701   | 0.4615   | 1.5      | 🔴      |
| Multigraphs      | 2%    | RFE         |         | 8.0  | yes | 0/0   | 4        | 0.2964   | 0.4873   | 0.9      | 🟢      |
| Multigraphs      | 2%    | Brute-force | L₀=0.01 | 12.0 | yes | 7/0   | 6        | 0.1425   | 0.2145   | 13.2     | 🔴      |
| Multigraphs      | 5%    | STRidge     | 10.0    | 8.0  | no  | 11/0  | 8        | 0.6867   | 0.6246   | 1.1      | 🔴      |
| Multigraphs      | 5%    | MSTLS       | 0.5     | 12.0 | no  | 14/0  | 9        | 0.5736   | 0.5425   | 1.5      | 🔴      |
| Multigraphs      | 5%    | RFE         |         | 8.0  | no  | 10/0  | 7        | 0.6939   | 0.6324   | 0.9      | 🔴      |
| Multigraphs      | 5%    | Brute-force | L₀=0.01 | 12.0 | no  | 14/0  | 8        | 0.5125   | 0.5371   | 11.0     | 🔴      |
| Multigraphs      | 5%    | STRidge     | 10.0    | 12.0 | yes | 14/0  | 10       | 0.5397   | 0.5571   | 1.1      | 🔴      |
| Multigraphs      | 5%    | MSTLS       | 0.5     | 12.0 | yes | 18/0  | 12       | 0.5210   | 0.5236   | 1.4      | 🔴      |
| Multigraphs      | 5%    | RFE         |         | 8.0  | yes | 10/0  | 7        | 0.6915   | 0.6328   | 0.9      | 🔴      |
| Multigraphs      | 5%    | Brute-force | L₀=0.01 | 12.0 | yes | 14/0  | 8        | 0.5121   | 0.5372   | 11.3     | 🔴      |
| Mod. Multigraphs | Clean | STRidge     | 10.0    | 12.0 | no  | 0/0   | 3        | 0.0136   | 0.0535   | 1.0      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS       | 0.25    | 12.0 | no  | 0/0   | 3        | 0.0140   | 0.0535   | 1.4      | 🟢      |
| Mod. Multigraphs | Clean | RFE         |         | 12.0 | no  | 0/0   | 3        | 0.0140   | 0.0535   | 0.9      | 🟢      |
| Mod. Multigraphs | Clean | STRidge     | 1.0     | 12.0 | yes | 0/0   | 4        | 0.0238   | 0.0227   | 1.3      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS       | 0.25    | 12.0 | yes | 0/0   | 3        | 0.0246   | 0.0248   | 1.4      | 🟢      |
| Mod. Multigraphs | Clean | RFE         |         | 12.0 | yes | 0/0   | 3        | 0.0246   | 0.0248   | 1.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge     | 10.0    | 12.0 | no  | 6/0   | 7        | 0.3881   | 0.3166   | 1.2      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS       | 0.5     | 8.0  | no  | 14/0  | 7        | 0.2769   | 0.3723   | 1.6      | 🔴      |
| Mod. Multigraphs | 2%    | RFE         |         | 8.0  | no  | 2/0   | 5        | 0.3742   | 0.4633   | 1.0      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge     | 10.0    | 12.0 | yes | 0/0   | 4        | 0.0645   | 0.1853   | 1.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS       | 0.25    | 12.0 | yes | 0/0   | 4        | 0.0647   | 0.1853   | 1.6      | 🟢      |
| Mod. Multigraphs | 2%    | RFE         |         | 12.0 | yes | 0/0   | 3        | 0.0292   | 0.1886   | 1.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge     | 10.0    | 12.0 | no  | 18/0  | 11       | 0.5243   | 0.5236   | 1.2      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS       | 0.5     | 12.0 | no  | 19/0  | 12       | 0.5230   | 0.5230   | 1.8      | 🔴      |
| Mod. Multigraphs | 5%    | RFE         |         | 12.0 | no  | 16/0  | 8        | 0.5158   | 0.5305   | 1.0      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge     | 10.0    | 12.0 | yes | 0/0   | 4        | 0.1965   | 0.4300   | 1.0      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS       | 0.25    | 12.0 | yes | 2/0   | 6        | 0.2008   | 0.4238   | 1.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE         |         | 12.0 | yes | 0/0   | 3        | 0.1234   | 0.4384   | 0.9      | 🟢      |

### Polynomial stencils

| Library          | Noise | Algorithm   | d_tol/α | m    | IBP | fp/fn | n_active | E_coeff | E_data   | Time (s) | Exact? |
|:-----------------|:------|:------------|:--------|:-----|:----|:------|:---------|:--------|:---------|:---------|:-------|
| PDE-FIND         | Clean | STRidge     | 1.0     | 4.0  | no  | 5/2   | 12       | 0.9882  | 0.7575   | 1.4      | 🔴      |
| PDE-FIND         | Clean | MSTLS       | 0.75    | 4.0  | no  | 21/2  | 36       | 0.9885  | 0.7547   | 3.5      | 🔴      |
| PDE-FIND         | Clean | RFE         |         | 8.0  | no  | 2/7   | 4        | 0.9725  | 0.7366   | 1.0      | 🔴      |
| PDE-FIND         | Clean | STRidge     | 0.1     | 12.0 | yes | 21/2  | 37       | 1.0245  | 0.6661   | 1.6      | 🔴      |
| PDE-FIND         | Clean | MSTLS       | 0.75    | 8.0  | yes | 29/2  | 45       | 0.9888  | 0.7123   | 3.5      | 🔴      |
| PDE-FIND         | Clean | RFE         |         | 8.0  | yes | 2/7   | 4        | 0.9725  | 0.7366   | 1.0      | 🔴      |
| PDE-FIND         | 2%    | STRidge     | 1.0     | 12.0 | no  | 21/2  | 34       | 1.0196  | 0.6684   | 1.6      | 🔴      |
| PDE-FIND         | 2%    | MSTLS       | 0.75    | 4.0  | no  | 21/2  | 36       | 0.9894  | 0.7715   | 3.7      | 🔴      |
| PDE-FIND         | 2%    | RFE         |         | 8.0  | no  | 2/7   | 4        | 0.9728  | 0.7439   | 1.0      | 🔴      |
| PDE-FIND         | 2%    | STRidge     | 1.0     | 12.0 | yes | 21/2  | 34       | 1.0192  | 0.6686   | 1.5      | 🔴      |
| PDE-FIND         | 2%    | MSTLS       | 0.75    | 12.0 | yes | 16/3  | 28       | 1.0147  | 0.6703   | 3.7      | 🔴      |
| PDE-FIND         | 2%    | RFE         |         | 8.0  | yes | 2/7   | 4        | 0.9728  | 0.7439   | 1.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge     | 1.0     | 12.0 | no  | 7/2   | 17       | 1.0039  | 0.7111   | 1.6      | 🔴      |
| PDE-FIND         | 5%    | MSTLS       | 0.75    | 4.0  | no  | 23/2  | 38       | 0.9926  | 0.8291   | 4.1      | 🔴      |
| PDE-FIND         | 5%    | RFE         |         | 8.0  | no  | 2/7   | 4        | 0.9793  | 0.7744   | 1.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge     | 1.0     | 8.0  | yes | 7/2   | 16       | 0.9915  | 0.7584   | 1.5      | 🔴      |
| PDE-FIND         | 5%    | MSTLS       | 0.75    | 8.0  | yes | 28/2  | 42       | 0.9905  | 0.7527   | 3.6      | 🔴      |
| PDE-FIND         | 5%    | RFE         |         | 8.0  | yes | 2/7   | 4        | 0.9793  | 0.7742   | 1.0      | 🔴      |
| WSINDy           | Clean | STRidge     | 10.0    | 12.0 | no  | 0/0   | 7        | 0.0111  | 0.0248   | 1.2      | 🟢      |
| WSINDy           | Clean | MSTLS       | 0.25    | 12.0 | no  | 18/0  | 12       | 12.8940 | 41.4123  | 2.3      | 🔴      |
| WSINDy           | Clean | RFE         |         | 12.0 | no  | 0/0   | 7        | 0.0113  | 0.0248   | 0.9      | 🟢      |
| WSINDy           | Clean | STRidge     | 1.0     | 4.0  | yes | 2/0   | 9        | 0.9278  | 0.7343   | 1.2      | 🔴      |
| WSINDy           | Clean | MSTLS       | 0.25    | 12.0 | yes | 10/0  | 10       | 13.2636 | 41.4248  | 2.4      | 🔴      |
| WSINDy           | Clean | RFE         |         | 12.0 | yes | 0/0   | 7        | 0.0121  | 0.0248   | 0.9      | 🟢      |
| WSINDy           | 2%    | STRidge     | 10.0    | 12.0 | no  | 0/0   | 7        | 0.0301  | 0.1884   | 1.3      | 🟢      |
| WSINDy           | 2%    | MSTLS       | 0.25    | 12.0 | no  | 18/0  | 12       | 13.9053 | 40.0414  | 2.7      | 🔴      |
| WSINDy           | 2%    | RFE         |         | 12.0 | no  | 0/0   | 7        | 0.0304  | 0.1884   | 0.9      | 🟢      |
| WSINDy           | 2%    | STRidge     | 10.0    | 12.0 | yes | 0/0   | 7        | 0.0307  | 0.1884   | 1.2      | 🟢      |
| WSINDy           | 2%    | MSTLS       | 0.25    | 8.0  | yes | 10/0  | 10       | 22.7095 | 176.3373 | 2.7      | 🔴      |
| WSINDy           | 2%    | RFE         |         | 12.0 | yes | 0/0   | 7        | 0.0310  | 0.1884   | 0.9      | 🟢      |
| WSINDy           | 5%    | STRidge     | 10.0    | 4.0  | no  | 10/0  | 11       | 1.0298  | 0.8364   | 1.3      | 🔴      |
| WSINDy           | 5%    | MSTLS       | 0.25    | 8.0  | no  | 25/0  | 18       | 27.9617 | 99.1766  | 2.9      | 🔴      |
| WSINDy           | 5%    | RFE         |         | 12.0 | no  | 2/0   | 9        | 0.1888  | 0.4327   | 1.0      | 🔴      |
| WSINDy           | 5%    | STRidge     | 10.0    | 12.0 | yes | 2/0   | 9        | 0.1908  | 0.4328   | 1.2      | 🔴      |
| WSINDy           | 5%    | MSTLS       | 0.25    | 8.0  | yes | 29/0  | 20       | 20.7747 | 95.4550  | 2.8      | 🔴      |
| WSINDy           | 5%    | RFE         |         | 12.0 | yes | 2/0   | 9        | 0.1888  | 0.4328   | 0.9      | 🔴      |
| PolyJet          | Clean | STRidge     | 10.0    | 8.0  | no  | 7/0   | 16       | 0.4577  | 0.3432   | 2.1      | 🔴      |
| PolyJet          | Clean | MSTLS       | 0.25    | 8.0  | no  | 9/0   | 18       | 0.4899  | 0.3296   | 7.3      | 🔴      |
| PolyJet          | Clean | RFE         |         | 8.0  | no  | 0/0   | 9        | 0.3098  | 0.4598   | 1.3      | 🟢      |
| PolyJet          | Clean | STRidge     | 10.0    | 8.0  | yes | 0/0   | 9        | 0.2798  | 0.4152   | 2.1      | 🟢      |
| PolyJet          | Clean | MSTLS       | 0.25    | 8.0  | yes | 5/0   | 14       | 0.3045  | 0.3320   | 7.1      | 🔴      |
| PolyJet          | Clean | RFE         |         | 4.0  | yes | 1/0   | 10       | 0.7633  | 0.6982   | 1.3      | 🔴      |
| PolyJet          | 2%    | STRidge     | 10.0    | 8.0  | no  | 9/0   | 18       | 0.4756  | 0.3910   | 2.1      | 🔴      |
| PolyJet          | 2%    | MSTLS       | 0.25    | 12.0 | no  | 14/0  | 23       | 0.4438  | 0.2441   | 7.5      | 🔴      |
| PolyJet          | 2%    | RFE         |         | 12.0 | no  | 0/0   | 9        | 0.3041  | 0.4180   | 1.3      | 🟢      |
| PolyJet          | 2%    | STRidge     | 10.0    | 8.0  | yes | 9/0   | 18       | 0.4707  | 0.3528   | 2.1      | 🔴      |
| PolyJet          | 2%    | MSTLS       | 0.25    | 12.0 | yes | 14/0  | 23       | 0.4427  | 0.2437   | 7.7      | 🔴      |
| PolyJet          | 2%    | RFE         |         | 12.0 | yes | 0/0   | 9        | 0.3044  | 0.4180   | 1.3      | 🟢      |
| PolyJet          | 5%    | STRidge     | 10.0    | 8.0  | no  | 8/0   | 17       | 0.7159  | 0.5704   | 2.2      | 🔴      |
| PolyJet          | 5%    | MSTLS       | 0.25    | 12.0 | no  | 10/0  | 19       | 0.5325  | 0.4459   | 8.5      | 🔴      |
| PolyJet          | 5%    | RFE         |         | 12.0 | no  | 0/5   | 4        | 0.7759  | 0.7313   | 1.3      | 🔴      |
| PolyJet          | 5%    | STRidge     | 10.0    | 8.0  | yes | 8/0   | 17       | 0.7055  | 0.5682   | 2.2      | 🔴      |
| PolyJet          | 5%    | MSTLS       | 0.25    | 12.0 | yes | 10/0  | 19       | 0.5318  | 0.4458   | 8.0      | 🔴      |
| PolyJet          | 5%    | RFE         |         | 12.0 | yes | 0/5   | 4        | 0.7759  | 0.7313   | 1.3      | 🔴      |
| Multigraphs      | Clean | STRidge     | 1.0     | 4.0  | no  | 12/0  | 9        | 0.9954  | 0.7267   | 1.0      | 🔴      |
| Multigraphs      | Clean | MSTLS       | 0.5     | 8.0  | no  | 16/0  | 10       | 1.0649  | 0.2959   | 1.4      | 🔴      |
| Multigraphs      | Clean | RFE         |         | 12.0 | no  | 0/0   | 4        | 0.2643  | 0.3654   | 0.9      | 🟢      |
| Multigraphs      | Clean | Brute-force | L₀=0.01 | 12.0 | no  | 15/0  | 9        | 0.8537  | 0.2328   | 11.4     | 🔴      |
| Multigraphs      | Clean | STRidge     | 10.0    | 4.0  | yes | 2/0   | 6        | 0.8303  | 0.6739   | 1.0      | 🔴      |
| Multigraphs      | Clean | MSTLS       | 0.25    | 8.0  | yes | 2/0   | 6        | 0.3398  | 0.3751   | 1.4      | 🔴      |
| Multigraphs      | Clean | RFE         |         | 12.0 | yes | 0/0   | 4        | 0.2647  | 0.3654   | 0.9      | 🟢      |
| Multigraphs      | Clean | Brute-force | L₀=0.01 | 4.0  | yes | 14/0  | 9        | 0.7141  | 0.6402   | 11.7     | 🔴      |
| Multigraphs      | 2%    | STRidge     | 1.0     | 4.0  | no  | 12/0  | 9        | 0.9970  | 0.7394   | 1.0      | 🔴      |
| Multigraphs      | 2%    | MSTLS       | 0.5     | 8.0  | no  | 16/0  | 10       | 1.0769  | 0.3808   | 1.4      | 🔴      |
| Multigraphs      | 2%    | RFE         |         | 12.0 | no  | 0/0   | 4        | 0.3009  | 0.4204   | 0.9      | 🟢      |
| Multigraphs      | 2%    | Brute-force | L₀=0.01 | 12.0 | no  | 15/0  | 9        | 0.9111  | 0.2929   | 11.2     | 🔴      |
| Multigraphs      | 2%    | STRidge     | 10.0    | 8.0  | yes | 16/0  | 10       | 0.8997  | 0.3700   | 1.0      | 🔴      |
| Multigraphs      | 2%    | MSTLS       | 0.5     | 8.0  | yes | 16/0  | 10       | 0.8996  | 0.3700   | 1.5      | 🔴      |
| Multigraphs      | 2%    | RFE         |         | 12.0 | yes | 0/0   | 4        | 0.3012  | 0.4205   | 0.9      | 🟢      |
| Multigraphs      | 2%    | Brute-force | L₀=0.01 | 8.0  | yes | 15/0  | 9        | 0.7957  | 0.3742   | 10.5     | 🔴      |
| Multigraphs      | 5%    | STRidge     | 1.0     | 4.0  | no  | 9/0   | 8        | 1.0016  | 0.7846   | 1.1      | 🔴      |
| Multigraphs      | 5%    | MSTLS       | 0.5     | 12.0 | no  | 17/0  | 11       | 0.9365  | 0.4689   | 1.7      | 🔴      |
| Multigraphs      | 5%    | RFE         |         | 8.0  | no  | 0/0   | 4        | 0.5945  | 0.7306   | 1.0      | 🟢      |
| Multigraphs      | 5%    | Brute-force | L₀=0.01 | 12.0 | no  | 15/0  | 9        | 0.9815  | 0.4792   | 11.7     | 🔴      |
| Multigraphs      | 5%    | STRidge     | 1.0     | 4.0  | yes | 7/0   | 8        | 1.0024  | 0.7883   | 1.0      | 🔴      |
| Multigraphs      | 5%    | MSTLS       | 0.5     | 8.0  | yes | 16/0  | 10       | 0.9401  | 0.5699   | 1.5      | 🔴      |
| Multigraphs      | 5%    | RFE         |         | 8.0  | yes | 0/0   | 4        | 0.5848  | 0.7218   | 0.9      | 🟢      |
| Multigraphs      | 5%    | Brute-force | L₀=0.01 | 8.0  | yes | 15/0  | 9        | 0.9740  | 0.5732   | 11.0     | 🔴      |
| Mod. Multigraphs | Clean | STRidge     | 10.0    | 12.0 | no  | 16/0  | 9        | 1.0598  | 0.2220   | 1.1      | 🔴      |
| Mod. Multigraphs | Clean | MSTLS       | 0.25    | 8.0  | no  | 0/0   | 4        | 0.3150  | 0.4693   | 1.7      | 🟢      |
| Mod. Multigraphs | Clean | RFE         |         | 12.0 | no  | 0/0   | 3        | 0.1935  | 0.4138   | 1.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge     | 10.0    | 12.0 | yes | 0/0   | 3        | 0.0111  | 0.0248   | 1.2      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS       | 0.25    | 12.0 | yes | 0/0   | 3        | 0.0113  | 0.0248   | 1.6      | 🟢      |
| Mod. Multigraphs | Clean | RFE         |         | 12.0 | yes | 0/0   | 3        | 0.0113  | 0.0248   | 1.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge     | 10.0    | 12.0 | no  | 10/0  | 9        | 0.3479  | 0.3394   | 1.1      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS       | 0.25    | 8.0  | no  | 0/0   | 4        | 0.3808  | 0.5422   | 1.6      | 🟢      |
| Mod. Multigraphs | 2%    | RFE         |         | 12.0 | no  | 0/0   | 3        | 0.2307  | 0.4642   | 0.9      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge     | 10.0    | 12.0 | yes | 0/0   | 4        | 0.0493  | 0.1857   | 1.3      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS       | 0.25    | 12.0 | yes | 0/0   | 4        | 0.0493  | 0.1857   | 1.6      | 🟢      |
| Mod. Multigraphs | 2%    | RFE         |         | 12.0 | yes | 0/0   | 3        | 0.0305  | 0.1886   | 1.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge     | 10.0    | 12.0 | no  | 16/0  | 8        | 0.9295  | 0.4787   | 1.3      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS       | 0.5     | 8.0  | no  | 9/0   | 8        | 0.6959  | 0.6056   | 1.7      | 🔴      |
| Mod. Multigraphs | 5%    | RFE         |         | 8.0  | no  | 0/0   | 3        | 0.5804  | 0.7523   | 1.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge     | 10.0    | 8.0  | yes | 8/0   | 6        | 0.6449  | 0.5856   | 1.1      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS       | 0.25    | 12.0 | yes | 0/0   | 4        | 0.1795  | 0.4218   | 1.4      | 🟢      |
| Mod. Multigraphs | 5%    | RFE         |         | 12.0 | yes | 0/0   | 3        | 0.1268  | 0.4384   | 0.9      | 🟢      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation | Backend | m | IBP | E_data floor | cond(Θ) |
|:------|:------------|:--------|:-:|:---:|-------------:|:--------|
---
| Source              | Wall time | CPU | Threads |
|:--------------------|----------:|:----|:-------:|
| `ch_merged_all.csv` | ? min     | ?   | ?       |
