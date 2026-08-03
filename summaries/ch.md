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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best         | Mg E_data | Std best       | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-------------------|:----------|:---------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | 🔴        | 17/0 ✗            | 0.9209     | —                  | 0.6452    | WSINDy/STRidge | 0/0 ✓     | 0.1205      | 0.2837     | —            |
| Clean | pointwise     | polynomial | —        | 35/0 ✗            | 0.9994     | —                  | 0.7536    | —              | 2/0 ✗     | 0.8250      | 0.8429     | —            |
| Clean | weak (no IBP) | spectral   | 🟥        | 0/0 ✓             | 0.0349     | multigraph/MSTLS   | 0.0734    | WSINDy/RFE     | 0/0 ✓     | 0.0162      | 0.0247     | —            |
| Clean | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.2643     | multigraph/MSTLS   | 0.3654    | WSINDy/STRidge | 0/0 ✓     | 0.0111      | 0.0248     | —            |
| Clean | weak (IBP)    | spectral   | 🟥        | 0/0 ✓             | 0.0332     | multigraph/STRidge | 0.0696    | WSINDy/STRidge | 0/0 ✓     | 0.0258      | 0.0440     | —            |
| Clean | weak (IBP)    | polynomial | 🟥        | 0/0 ✓             | 0.2572     | multigraph/RFE     | 0.3984    | WSINDy/STRidge | 0/0 ✓     | 0.0260      | 0.0700     | —            |
| 2%    | pointwise     | spectral   | —        | 17/0 ✗            | 0.9775     | —                  | 0.8748    | —              | 18/0 ✗    | 0.9782      | 0.8863     | —            |
| 2%    | pointwise     | polynomial | —        | 18/0 ✗            | 0.9995     | —                  | 0.8168    | —              | 18/0 ✗    | 1.0458      | 0.9283     | —            |
| 2%    | weak (no IBP) | spectral   | 🔴        | 2/0 ✗             | 0.4197     | —                  | 0.4608    | WSINDy/STRidge | 0/0 ✓     | 0.0344      | 0.1886     | —            |
| 2%    | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.2981     | multigraph/RFE     | 0.4204    | WSINDy/STRidge | 0/0 ✓     | 0.0301      | 0.1884     | —            |
| 2%    | weak (IBP)    | spectral   | 🔴        | 2/0 ✗             | 0.3935     | —                  | 0.4387    | PolyJet/RFE    | 0/0 ✓     | 0.2871      | 0.4617     | —            |
| 2%    | weak (IBP)    | polynomial | 🟢        | 0/0 ✓             | 0.3210     | multigraph/RFE     | 0.4815    | —              | 2/0 ✗     | 0.1613      | 0.3059     | —            |
| 5%    | pointwise     | spectral   | —        | 17/0 ✗            | 0.9942     | —                  | 0.9697    | —              | 37/0 ✗    | 0.9890      | 0.9651     | —            |
| 5%    | pointwise     | polynomial | —        | 35/0 ✗            | 0.9997     | —                  | 0.9213    | —              | 62/0 ✗    | 0.9999      | 0.9189     | —            |
| 5%    | weak (no IBP) | spectral   | —        | 10/0 ✗            | 1.6238     | —                  | 0.5544    | —              | 2/0 ✗     | 0.2005      | 0.4331     | —            |
| 5%    | weak (no IBP) | polynomial | 🟢        | 0/0 ✓             | 0.5894     | multigraph/RFE     | 0.7306    | —              | 2/0 ✗     | 0.1894      | 0.4327     | —            |
| 5%    | weak (IBP)    | spectral   | —        | 10/0 ✗            | 1.6645     | —                  | 0.6086    | —              | 2/0 ✗     | 0.5394      | 0.6054     | —            |
| 5%    | weak (IBP)    | polynomial | 🟢        | 0/0 ✓             | 0.5387     | multigraph/RFE     | 0.6926    | —              | 2/0 ✗     | 0.5143      | 0.6004     | —            |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                | Best conf.                      | fp/fn | E_coeff  | Time (s) | vs. best |
|:----------------------|:--------------------------------|:------|---------:|---------:|---------:|
| WSINDy + STRidge      | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0111   | 2.1      | —        |
| WSINDy + RFE          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0113   | 0.6      | +2.1%    |
| Multigraphs + STRidge | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0332   | 1.2      | +200.7%  |
| Multigraphs + RFE     | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0334   | 0.6      | +201.7%  |
| Multigraphs + MSTLS   | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0334   | 1.6      | +201.7%  |
| PolyJet + MSTLS       | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0519   | 16.5     | +369.9%  |
| PolyJet + RFE         | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0519   | 0.8      | +369.9%  |
| PolyJet + STRidge     | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0525   | 4.8      | +374.8%  |
| WSINDy + MSTLS        | spectral, weak (IBP), m=12      | 2/0 ✗ | 636.2848 | 3.5      | ✗        |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=4    | 5/2 ✗ | 1.1966   | 4.1      | ✗        |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12   | 9/2 ✗ | 4.9996   | 8.0      | ✗        |
| PDE-FIND + RFE        | polynomial, pointwise           | 2/6 ✗ | 0.9945   | 0.6      | ✗        |

### 2% noise

| Method                | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| WSINDy + STRidge      | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0301  | 1.9      | —        |
| WSINDy + RFE          | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0304  | 0.6      | +1.0%    |
| PolyJet + RFE         | spectral, weak (IBP), m=12      | 0/0 ✓  | 0.2871  | 0.8      | +853.8%  |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.2981  | 0.6      | +890.2%  |
| PolyJet + STRidge     | polynomial, weak (IBP), m=12    | 2/0 ✗  | 0.3924  | 4.3      | ✗        |
| Multigraphs + MSTLS   | spectral, weak (IBP), m=12      | 2/0 ✗  | 0.3935  | 1.6      | ✗        |
| Multigraphs + STRidge | polynomial, weak (IBP), m=8     | 2/0 ✗  | 0.8645  | 1.1      | ✗        |
| PolyJet + MSTLS       | spectral, weak (IBP), m=12      | 6/0 ✗  | 1.8419  | 17.5     | ✗        |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=12 | 18/0 ✗ | 13.5319 | 3.7      | ✗        |
| PDE-FIND + STRidge    | spectral, pointwise             | 5/2 ✗  | 0.9877  | 4.3      | ✗        |
| PDE-FIND + MSTLS      | polynomial, weak (IBP), m=12    | 28/2 ✗ | 0.9909  | 6.5      | ✗        |
| PDE-FIND + RFE        | spectral, pointwise             | 2/5 ✗  | 0.9774  | 0.6      | ✗        |

### 5% noise

| Method                | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓  | 0.5387  | 0.6      | —        |
| WSINDy + RFE          | polynomial, weak (no IBP), m=12 | 2/0 ✗  | 0.1894  | 0.6      | ✗        |
| WSINDy + STRidge      | spectral, weak (no IBP), m=4    | 2/0 ✗  | 1.0081  | 2.4      | ✗        |
| Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 7/0 ✗  | 0.6634  | 1.2      | ✗        |
| PolyJet + STRidge     | polynomial, weak (IBP), m=12    | 8/0 ✗  | 0.6686  | 4.5      | ✗        |
| PolyJet + RFE         | spectral, weak (no IBP), m=12   | 8/0 ✗  | 1.5522  | 0.8      | ✗        |
| WSINDy + MSTLS        | spectral, weak (no IBP), m=4    | 8/0 ✗  | 40.4156 | 3.8      | ✗        |
| Multigraphs + MSTLS   | polynomial, weak (IBP), m=12    | 9/0 ✗  | 0.6451  | 1.4      | ✗        |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=12 | 10/0 ✗ | 0.5376  | 19.3     | ✗        |
| PDE-FIND + STRidge    | polynomial, weak (no IBP), m=4  | 5/2 ✗  | 0.9951  | 4.7      | ✗        |
| PDE-FIND + MSTLS      | polynomial, weak (IBP), m=12    | 28/2 ✗ | 0.9974  | 6.6      | ✗        |
| PDE-FIND + RFE        | spectral, pointwise             | 3/4 ✗  | 0.9929  | 0.7      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 36/2  | 43       | 0.9731  | 0.7703 | 5.7      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 2/7   | 4        | 0.9769  | 0.7734 | 0.6      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.01    | 5/2   | 12       | 0.9877  | 0.8922 | 4.3      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 2/5   | 6        | 0.9774  | 0.9309 | 0.6      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.01    | 15/2  | 22       | 0.9948  | 0.9722 | 4.6      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 3/4   | 8        | 0.9929  | 0.9775 | 0.7      | 🔴      |
| WSINDy      | Clean | STRidge   | 1.0     | 0/0   | 7        | 0.1205  | 0.2837 | 2.1      | 🟢      |
| WSINDy      | Clean | RFE       |         | 0/0   | 7        | 0.1215  | 0.2837 | 0.6      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.01    | 18/0  | 12       | 0.9782  | 0.8863 | 2.7      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 2/5   | 6        | 0.9774  | 0.9309 | 0.6      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 96/0  | 46       | 0.9795  | 0.9654 | 4.6      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 6/3   | 8        | 0.9829  | 0.9762 | 0.6      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.1     | 15/0  | 24       | 0.9220  | 0.6098 | 6.2      | 🔴      |
| PolyJet     | Clean | RFE       |         | 2/2   | 9        | 0.9017  | 0.7424 | 0.8      | 🔴      |
| PolyJet     | 2%    | STRidge   | 0.01    | 78/0  | 87       | 0.9754  | 0.8639 | 16.2     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 2/5   | 6        | 0.9774  | 0.9309 | 0.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.01    | 37/0  | 46       | 0.9890  | 0.9651 | 9.4      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 3/4   | 8        | 0.9929  | 0.9775 | 0.7      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 17/0  | 11       | 0.9209  | 0.6452 | 1.5      | 🔴      |
| Multigraphs | Clean | RFE       |         | 2/2   | 5        | 0.9012  | 0.7430 | 0.6      | 🔴      |
| Multigraphs | 2%    | STRidge   | 0.01    | 17/0  | 11       | 0.9775  | 0.8748 | 1.5      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 7/7   | 5        | 0.9873  | 0.9138 | 0.6      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.01    | 17/0  | 11       | 0.9942  | 0.9697 | 1.5      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 5/2   | 6        | 0.9977  | 0.9766 | 0.6      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.001   | 53/2  | 60       | 0.9975  | 0.8338 | 8.1      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 2/6   | 5        | 0.9945  | 0.9004 | 0.6      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.001   | 53/2  | 60       | 0.9982  | 0.8857 | 8.0      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 2/6   | 5        | 0.9954  | 0.9266 | 0.6      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.01    | 17/2  | 24       | 0.9994  | 0.9638 | 4.5      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 2/6   | 5        | 0.9975  | 0.9744 | 0.6      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 2/0   | 9        | 0.8288  | 0.8430 | 2.2      | 🔴      |
| WSINDy      | Clean | RFE       |         | 2/0   | 9        | 0.8250  | 0.8429 | 0.6      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.01    | 18/0  | 12       | 1.0458  | 0.9283 | 2.8      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 2/8   | 3        | 0.9973  | 0.9756 | 0.6      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 83/0  | 42       | 1.0288  | 0.9763 | 4.7      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 2/7   | 4        | 0.9960  | 0.9813 | 0.6      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.01    | 35/0  | 44       | 0.9994  | 0.7258 | 11.6     | 🔴      |
| PolyJet     | Clean | RFE       |         | 4/7   | 6        | 0.9996  | 0.8642 | 0.8      | 🔴      |
| PolyJet     | 2%    | STRidge   | 0.1     | 31/0  | 40       | 0.9997  | 0.8027 | 7.0      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 3/9   | 3        | 0.9999  | 0.9211 | 0.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.01    | 62/0  | 71       | 0.9999  | 0.9189 | 12.7     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 3/9   | 3        | 1.0000  | 0.9627 | 0.7      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.001   | 35/0  | 20       | 0.9994  | 0.7536 | 2.0      | 🔴      |
| Multigraphs | Clean | RFE       |         | 6/9   | 3        | 1.0036  | 0.8544 | 0.6      | 🔴      |
| Multigraphs | 2%    | STRidge   | 0.01    | 18/0  | 12       | 0.9995  | 0.8168 | 1.5      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 6/9   | 3        | 1.0034  | 0.8879 | 0.6      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.001   | 35/0  | 20       | 0.9997  | 0.9213 | 1.9      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 6/9   | 3        | 1.0027  | 0.9565 | 0.6      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff  | E_data   | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:---------|:---------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 4  | no  | 5/2   | 12       | 1.1966   | 0.7449   | 4.1      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.5     | 12 | no  | 9/2   | 16       | 4.9996   | 0.6064   | 8.0      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | no  | 2/7   | 4        | 0.9654   | 0.7371   | 0.6      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 10.0    | 12 | yes | 5/2   | 12       | 4.1175   | 0.6603   | 2.5      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.5     | 12 | yes | 12/2  | 19       | 4.0909   | 0.6566   | 8.0      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 4/7   | 6        | 1.4373   | 0.7980   | 0.6      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 8  | no  | 14/2  | 21       | 1.2058   | 0.7283   | 3.5      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 8  | no  | 31/2  | 38       | 1.2307   | 0.7230   | 6.6      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | no  | 2/7   | 4        | 0.9658   | 0.7441   | 0.6      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 12 | yes | 29/2  | 36       | 1.3175   | 0.7134   | 5.2      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | yes | 31/2  | 38       | 1.3207   | 0.7125   | 6.3      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | yes | 2/7   | 4        | 0.9764   | 0.7606   | 0.6      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 8  | no  | 10/2  | 17       | 1.1240   | 0.7631   | 3.1      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 4  | no  | 2/5   | 6        | 1.0010   | 0.8351   | 6.2      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 4  | no  | 2/7   | 4        | 0.9783   | 0.8424   | 0.6      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.01    | 12 | yes | 26/2  | 33       | 1.1791   | 0.7522   | 5.4      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 8  | yes | 4/5   | 8        | 0.9933   | 0.8046   | 6.3      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 2/7   | 4        | 0.9776   | 0.8073   | 0.6      | 🔴      |
| WSINDy      | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0163   | 0.0247   | 2.0      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 8  | no  | 2/0   | 9        | 704.2663 | 209.3462 | 3.4      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0162   | 0.0247   | 0.6      | 🟢      |
| WSINDy      | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0258   | 0.0440   | 1.9      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 2/0   | 9        | 636.2848 | 170.4945 | 3.5      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0260   | 0.0440   | 0.7      | 🟢      |
| WSINDy      | 2%    | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0344   | 0.1886   | 1.9      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | no  | 18/0  | 12       | 426.8456 | 40.0176  | 3.9      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0356   | 0.1886   | 0.6      | 🟢      |
| WSINDy      | 2%    | STRidge   | 1.0     | 12 | yes | 2/0   | 9        | 0.1742   | 0.3077   | 2.1      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | yes | 18/0  | 12       | 641.9611 | 148.2076 | 3.5      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 2/0   | 9        | 0.1709   | 0.3077   | 0.6      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 4  | no  | 2/0   | 9        | 1.0081   | 0.8300   | 2.4      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 4  | no  | 8/0   | 15       | 40.4156  | 37.5187  | 3.8      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 2/0   | 9        | 0.2005   | 0.4331   | 0.6      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | yes | 10/0  | 11       | 10.4241  | 0.5973   | 2.2      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | yes | 26/0  | 19       | 643.0872 | 87.2160  | 3.8      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | yes | 2/0   | 9        | 0.5394   | 0.6054   | 0.6      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 8  | no  | 0/0   | 9        | 0.0525   | 0.0726   | 4.8      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | no  | 0/0   | 9        | 0.0519   | 0.0726   | 16.5     | 🟢      |
| PolyJet     | Clean | RFE       |         | 8  | no  | 0/0   | 9        | 0.0519   | 0.0726   | 0.8      | 🟢      |
| PolyJet     | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 9        | 0.0536   | 0.0686   | 4.9      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 9        | 0.0534   | 0.0686   | 16.7     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 9        | 0.0534   | 0.0686   | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 10.0    | 12 | no  | 4/0   | 13       | 1.5168   | 0.2518   | 4.6      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 8  | no  | 11/0  | 20       | 1.9260   | 0.2965   | 16.7     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 8  | no  | 0/0   | 9        | 0.2987   | 0.4864   | 0.7      | 🟢      |
| PolyJet     | 2%    | STRidge   | 10.0    | 12 | yes | 2/0   | 11       | 1.0652   | 0.3902   | 4.4      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 12 | yes | 6/0   | 15       | 1.8419   | 0.3049   | 17.5     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 0/0   | 9        | 0.2871   | 0.4617   | 0.8      | 🟢      |
| PolyJet     | 5%    | STRidge   | 10.0    | 8  | no  | 8/0   | 17       | 1.6224   | 0.5896   | 4.6      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | no  | 14/0  | 23       | 1.8560   | 0.4103   | 17.6     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | no  | 8/0   | 17       | 1.5522   | 0.4480   | 0.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 10.0    | 12 | yes | 10/0  | 19       | 1.6524   | 0.5435   | 5.0      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | yes | 13/0  | 22       | 1.6643   | 0.5436   | 17.5     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | yes | 8/0   | 17       | 1.6702   | 0.5597   | 0.7      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 4        | 0.0351   | 0.0734   | 1.1      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0349   | 0.0734   | 1.4      | 🟢      |
| Multigraphs | Clean | RFE       |         | 8  | no  | 0/0   | 4        | 0.0349   | 0.0734   | 0.6      | 🟢      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 4        | 0.0332   | 0.0696   | 1.2      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0334   | 0.0696   | 1.6      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0334   | 0.0696   | 0.6      | 🟢      |
| Multigraphs | 2%    | STRidge   | 10.0    | 12 | no  | 7/0   | 6        | 1.6244   | 0.2144   | 1.2      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 8  | no  | 2/0   | 6        | 0.4197   | 0.4608   | 1.6      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 8  | no  | 2/0   | 6        | 0.4197   | 0.4608   | 0.6      | 🔴      |
| Multigraphs | 2%    | STRidge   | 10.0    | 12 | yes | 5/0   | 6        | 0.9862   | 0.3700   | 1.1      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 2/0   | 6        | 0.3935   | 0.4387   | 1.6      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 2/0   | 6        | 0.3935   | 0.4387   | 0.6      | 🔴      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 10/0  | 7        | 1.6238   | 0.5544   | 1.3      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.5     | 12 | no  | 14/0  | 9        | 2.1135   | 0.5425   | 1.4      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 8  | no  | 10/0  | 7        | 1.6248   | 0.6324   | 0.6      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 8  | yes | 12/0  | 9        | 1.0733   | 0.7703   | 1.3      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 8  | yes | 17/0  | 13       | 1.0691   | 0.7683   | 1.5      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | yes | 10/0  | 7        | 1.6645   | 0.6086   | 0.6      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data   | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:---------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 12 | no  | 21/2  | 28       | 1.0248  | 0.6667   | 4.1      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 4  | no  | 29/2  | 36       | 0.9885  | 0.7547   | 6.5      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | no  | 2/7   | 4        | 0.9725  | 0.7366   | 0.6      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 0.1     | 8  | yes | 29/2  | 36       | 0.9764  | 0.7373   | 5.5      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 12 | yes | 28/3  | 34       | 0.9896  | 0.7034   | 7.1      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 2/7   | 4        | 0.9806  | 0.7503   | 0.6      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 12 | no  | 23/2  | 30       | 1.0229  | 0.6740   | 4.2      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 4  | no  | 29/2  | 36       | 0.9901  | 0.7715   | 6.9      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | no  | 2/7   | 4        | 0.9728  | 0.7439   | 0.6      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 8  | yes | 28/2  | 35       | 0.9781  | 0.7509   | 5.1      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | yes | 28/2  | 35       | 0.9909  | 0.7102   | 6.5      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | yes | 2/7   | 4        | 0.9811  | 0.7606   | 0.6      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 4  | no  | 5/2   | 12       | 0.9951  | 0.8307   | 4.7      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.75    | 4  | no  | 31/2  | 38       | 0.9951  | 0.8291   | 6.5      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | no  | 2/7   | 4        | 0.9798  | 0.7744   | 0.6      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 8  | yes | 46/2  | 53       | 0.9845  | 0.8034   | 7.0      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.75    | 12 | yes | 28/2  | 35       | 0.9974  | 0.7419   | 6.6      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 2/7   | 4        | 0.9831  | 0.8073   | 0.6      | 🔴      |
| WSINDy      | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0111  | 0.0248   | 2.1      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | no  | 18/0  | 12       | 12.8940 | 41.4123  | 3.3      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0113  | 0.0248   | 0.6      | 🟢      |
| WSINDy      | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0260  | 0.0700   | 2.0      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 18/0  | 12       | 23.8549 | 169.1508 | 4.0      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0264  | 0.0700   | 0.6      | 🟢      |
| WSINDy      | 2%    | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0301  | 0.1884   | 1.9      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | no  | 18/0  | 12       | 13.5319 | 40.0414  | 3.7      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0304  | 0.1884   | 0.6      | 🟢      |
| WSINDy      | 2%    | STRidge   | 1.0     | 12 | yes | 2/0   | 9        | 0.1638  | 0.3060   | 2.1      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 8  | yes | 20/0  | 17       | 10.3382 | 184.3150 | 4.0      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 2/0   | 9        | 0.1613  | 0.3059   | 0.6      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 4  | no  | 6/0   | 10       | 1.0178  | 0.8383   | 1.9      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 8  | no  | 25/0  | 18       | 29.3716 | 99.1766  | 4.0      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 2/0   | 9        | 0.1894  | 0.4327   | 0.6      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | yes | 10/0  | 11       | 0.9801  | 0.5947   | 2.2      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | yes | 26/0  | 19       | 22.5246 | 87.7398  | 3.9      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | yes | 2/0   | 9        | 0.5143  | 0.6004   | 0.6      | 🔴      |
| PolyJet     | Clean | STRidge   | 10.0    | 8  | no  | 7/0   | 16       | 0.4577  | 0.3432   | 4.5      | 🔴      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | no  | 9/0   | 18       | 0.4899  | 0.3296   | 19.7     | 🔴      |
| PolyJet     | Clean | RFE       |         | 8  | no  | 0/0   | 9        | 0.3098  | 0.4598   | 0.8      | 🟢      |
| PolyJet     | Clean | STRidge   | 1.0     | 8  | yes | 7/0   | 16       | 0.6284  | 0.4652   | 5.7      | 🔴      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 12/0  | 21       | 0.4424  | 0.1950   | 19.3     | 🔴      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 9        | 0.2548  | 0.3919   | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 10.0    | 8  | no  | 9/0   | 18       | 0.4746  | 0.3910   | 4.9      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 12 | no  | 14/0  | 23       | 0.4445  | 0.2441   | 19.2     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 0/0   | 9        | 0.3007  | 0.4180   | 0.7      | 🟢      |
| PolyJet     | 2%    | STRidge   | 10.0    | 12 | yes | 2/0   | 11       | 0.3924  | 0.4339   | 4.3      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 12 | yes | 7/0   | 16       | 0.4094  | 0.3995   | 17.2     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 2/0   | 11       | 0.3912  | 0.4339   | 0.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 10.0    | 8  | no  | 8/0   | 17       | 0.7098  | 0.5704   | 4.7      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | no  | 10/0  | 19       | 0.5376  | 0.4459   | 19.3     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | no  | 0/5   | 4        | 0.7467  | 0.7313   | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 10.0    | 12 | yes | 8/0   | 17       | 0.6686  | 0.5456   | 4.5      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.5     | 12 | yes | 16/0  | 25       | 0.6853  | 0.5175   | 17.4     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | yes | 2/6   | 5        | 0.9612  | 0.7856   | 0.7      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 15/0  | 9        | 0.4528  | 0.1873   | 1.4      | 🔴      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.2643  | 0.3654   | 1.6      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.2643  | 0.3654   | 0.6      | 🟢      |
| Multigraphs | Clean | STRidge   | 1.0     | 12 | yes | 15/0  | 11       | 0.4300  | 0.2855   | 1.3      | 🔴      |
| Multigraphs | Clean | MSTLS     | 0.5     | 8  | yes | 16/0  | 10       | 0.7119  | 0.5211   | 1.5      | 🔴      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.2572  | 0.3984   | 0.6      | 🟢      |
| Multigraphs | 2%    | STRidge   | 10.0    | 8  | no  | 7/0   | 6        | 0.4980  | 0.4629   | 1.2      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 8  | no  | 2/0   | 6        | 0.4692  | 0.4872   | 1.6      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.2981  | 0.4204   | 0.6      | 🟢      |
| Multigraphs | 2%    | STRidge   | 10.0    | 8  | yes | 2/0   | 6        | 0.8645  | 0.7105   | 1.1      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.5     | 12 | yes | 17/0  | 11       | 0.5221  | 0.3394   | 1.5      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.3210  | 0.4815   | 0.6      | 🟢      |
| Multigraphs | 5%    | STRidge   | 10.0    | 8  | no  | 7/0   | 6        | 0.7128  | 0.6195   | 1.2      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.5     | 8  | no  | 9/0   | 8        | 0.6970  | 0.6056   | 1.6      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 8  | no  | 0/0   | 4        | 0.5894  | 0.7306   | 0.6      | 🟢      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | yes | 7/0   | 6        | 0.6634  | 0.5983   | 1.2      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.5     | 12 | yes | 9/0   | 8        | 0.6451  | 0.5815   | 1.4      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.5387  | 0.6926   | 0.6      | 🟢      |

---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ch`.

| Source                              | Wall time | CPU                                      | Threads |
|:------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260802_005618_ch_full.csv` | 364.6 min | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | 1       |
