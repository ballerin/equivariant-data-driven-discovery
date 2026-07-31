# Kuramoto-Sivashinsky equation

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

KS contains the gradient norm $\|\nabla f\|^2 = f_x^2 + f_y^2$ which neither PDE-FIND nor WSINDy can represent as a sparse column (PDE-FIND lacks products of derivatives; in WSINDy, $f_x^2$ is locked with $f\cdot f_{xx}$ inside $\partial_{xx}(f^2)$).  Three KS regimes are simulated: cellular ($L=4\pi$), weak-chaotic ($L=8\pi$), chaotic ($L=16\pi$).

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
| Clean | pointwise     | spectral   | ⚪        | 0/0 ✓             | 0.0047     | multigraph/RFE     | 0.0176    | PolyJet/RFE     | 0/0 ✓     | 0.0047      | 0.0175     | 0.0182       |
| Clean | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.0117     | multigraph/RFE     | 0.0325    | PolyJet/RFE     | 0/0 ✓     | 0.0128      | 0.0320     | 0.0345       |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0022     | multigraph/RFE     | 0.0112    | PolyJet/MSTLS   | 0/0 ✓     | 0.0023      | 0.0112     | 0.0114       |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0022     | multigraph/MSTLS   | 0.0112    | PolyJet/MSTLS   | 0/0 ✓     | 0.0023      | 0.0112     | 0.0114       |
| Clean | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0029     | multigraph/STRidge | 0.0133    | PolyJet/STRidge | 0/0 ✓     | 0.0030      | 0.0133     | 0.0136       |
| Clean | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0029     | multigraph/STRidge | 0.0133    | PolyJet/STRidge | 0/0 ✓     | 0.0029      | 0.0133     | 0.0136       |
| 2%    | pointwise     | spectral   | —        | 18/0 ✗            | 1.0000     | —                  | 0.6378    | —               | 73/0 ✗    | 1.0000      | 0.6364     | 1381.7721    |
| 2%    | pointwise     | polynomial | —        | 25/0 ✗            | 0.9990     | —                  | 0.6210    | —               | 86/0 ✗    | 0.9992      | 0.6375     | 29.8498      |
| 2%    | weak (no IBP) | spectral   | —        | 6/0 ✗             | 0.9307     | —                  | 0.3476    | —               | 66/0 ✗    | 0.7213      | 0.2564     | 11.3180      |
| 2%    | weak (no IBP) | polynomial | —        | 7/0 ✗             | 0.1789     | —                  | 0.1168    | —               | 9/0 ✗     | 0.1883      | 0.1099     | 0.8842       |
| 2%    | weak (IBP)    | spectral   | —        | 37/0 ✗            | 0.9279     | —                  | 0.3341    | —               | 11/0 ✗    | 0.9039      | 0.3255     | 11.3156      |
| 2%    | weak (IBP)    | polynomial | —        | 9/0 ✗             | 0.7271     | —                  | 0.2168    | —               | 67/0 ✗    | 0.9911      | 0.3119     | 0.9437       |
| 5%    | pointwise     | spectral   | —        | 34/0 ✗            | 1.0000     | —                  | 0.8725    | —               | 86/0 ✗    | 1.0000      | 0.8714     | 2190.0271    |
| 5%    | pointwise     | polynomial | —        | 37/0 ✗            | 0.9994     | —                  | 0.8717    | —               | 86/0 ✗    | 0.9994      | 0.8702     | 47.4252      |
| 5%    | weak (no IBP) | spectral   | —        | 28/0 ✗            | 1.0001     | —                  | 0.5369    | —               | 42/0 ✗    | 0.9995      | 0.4133     | 69.4094      |
| 5%    | weak (no IBP) | polynomial | —        | 8/0 ✗             | 0.7278     | —                  | 0.2940    | —               | 13/0 ✗    | 0.6399      | 0.2611     | 5.4493       |
| 5%    | weak (IBP)    | spectral   | —        | 27/0 ✗            | 1.0005     | —                  | 0.4893    | —               | 78/0 ✗    | 1.0001      | 0.4834     | 68.8768      |
| 5%    | weak (IBP)    | polynomial | —        | 36/0 ✗            | 0.9995     | —                  | 0.4641    | —               | 17/0 ✗    | 0.9582      | 0.3414     | 5.4607       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                | Best conf.                     | fp/fn  | E_coeff     | Time (s) | vs. best |
|:----------------------|:-------------------------------|:-------|------------:|---------:|---------:|
| Multigraphs + RFE     | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0022      | 0.1      | —        |
| Multigraphs + MSTLS   | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0022      | 0.9      | +0.0%    |
| PolyJet + MSTLS       | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0023      | 41.6     | +4.6%    |
| PolyJet + RFE         | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0023      | 0.7      | +4.6%    |
| Multigraphs + STRidge | polynomial, weak (IBP), m=12   | 0/0 ✓  | 0.0029      | 1.6      | +29.9%   |
| PolyJet + STRidge     | polynomial, weak (IBP), m=12   | 0/0 ✓  | 0.0029      | 11.8     | +32.2%   |
| WSINDy + STRidge      | spectral, weak (IBP), m=8      | 75/0 ✗ | 0.4316      | 4.6      | ✗        |
| PDE-FIND + STRidge    | spectral, pointwise            | 25/2 ✗ | 0.5503      | 7.0      | ✗        |
| PDE-FIND + MSTLS      | spectral, weak (IBP), m=12     | 7/5 ✗  | 0.4412      | 3.2      | ✗        |
| PDE-FIND + RFE        | spectral, pointwise            | 1/7 ✗  | 0.4409      | 0.6      | ✗        |
| WSINDy + RFE          | spectral, pointwise            | 1/7 ✗  | 0.4409      | 0.3      | ✗        |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=4 | 4/7 ✗  | 444514.6421 | 3.5      | ✗        |

### 2% noise

| Method                | Best conf.                      | fp/fn  | E_coeff  | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|---------:|---------:|---------:|
| Multigraphs + MSTLS   | spectral, weak (no IBP), m=8    | 6/0 ✗  | 0.9307   | 1.0      | ✗        |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=12 | 9/0 ✗  | 0.1883   | 39.1     | ✗        |
| PolyJet + STRidge     | polynomial, weak (no IBP), m=12 | 17/0 ✗ | 0.2182   | 12.8     | ✗        |
| Multigraphs + STRidge | spectral, pointwise             | 18/0 ✗ | 1.0000   | 1.8      | ✗        |
| WSINDy + STRidge      | spectral, weak (IBP), m=12      | 83/0 ✗ | 1.0027   | 4.3      | ✗        |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=4    | 11/2 ✗ | 0.9998   | 5.9      | ✗        |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.9984   | 0.1      | ✗        |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12   | 7/5 ✗  | 1.0113   | 3.8      | ✗        |
| PolyJet + RFE         | polynomial, weak (no IBP), m=8  | 0/6 ✗  | 0.9987   | 0.8      | ✗        |
| PDE-FIND + RFE        | polynomial, pointwise           | 1/7 ✗  | 0.9991   | 0.4      | ✗        |
| WSINDy + RFE          | polynomial, pointwise           | 1/7 ✗  | 0.9991   | 0.3      | ✗        |
| WSINDy + MSTLS        | spectral, weak (no IBP), m=4    | 4/7 ✗  | 268.7026 | 3.0      | ✗        |

### 5% noise

| Method                | Best conf.                      | fp/fn  | E_coeff  | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|---------:|---------:|---------:|
| Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 8/0 ✗  | 0.7278   | 0.8      | ✗        |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=12 | 13/0 ✗ | 0.6399   | 37.3     | ✗        |
| Multigraphs + STRidge | spectral, weak (IBP), m=8       | 27/0 ✗ | 1.0005   | 1.8      | ✗        |
| PolyJet + STRidge     | spectral, weak (no IBP), m=8    | 42/0 ✗ | 0.9995   | 15.8     | ✗        |
| WSINDy + STRidge      | spectral, weak (no IBP), m=8    | 79/0 ✗ | 1.0015   | 4.8      | ✗        |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=4    | 12/2 ✗ | 0.9999   | 5.8      | ✗        |
| Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/5 ✗  | 0.9988   | 0.2      | ✗        |
| PDE-FIND + MSTLS      | polynomial, weak (IBP), m=12    | 1/7 ✗  | 0.9994   | 3.5      | ✗        |
| PDE-FIND + RFE        | polynomial, weak (IBP), m=12    | 1/7 ✗  | 0.9996   | 0.4      | ✗        |
| WSINDy + RFE          | polynomial, weak (IBP), m=12    | 1/7 ✗  | 0.9996   | 0.3      | ✗        |
| PolyJet + RFE         | polynomial, weak (IBP), m=12    | 1/7 ✗  | 0.9996   | 0.7      | ✗        |
| WSINDy + MSTLS        | spectral, weak (no IBP), m=4    | 4/7 ✗  | 100.5003 | 3.3      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 25/2  | 39       | 0.5503  | 0.4154 | 7.0      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 1/7   | 1        | 0.4409  | 0.4294 | 0.6      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 36/2  | 54       | 1.0000  | 0.6377 | 7.1      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6433 | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.01    | 35/2  | 48       | 1.0000  | 0.8721 | 6.7      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.8979 | 0.4      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.001   | 98/0  | 46       | 0.4276  | 0.4152 | 6.3      | 🔴      |
| WSINDy      | Clean | RFE       |         | 1/7   | 1        | 0.4409  | 0.4294 | 0.3      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 91/0  | 43       | 1.0000  | 0.6378 | 4.8      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6433 | 0.4      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 88/0  | 43       | 1.0000  | 0.8722 | 5.2      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.8979 | 0.5      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 0/0   | 7        | 0.0047  | 0.0175 | 13.6     | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.0047  | 0.0175 | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.1     | 73/0  | 80       | 1.0000  | 0.6364 | 19.5     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6433 | 0.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 86/0  | 93       | 1.0000  | 0.8714 | 22.1     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.8979 | 0.8      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 0/0   | 3        | 0.0047  | 0.0176 | 2.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0047  | 0.0176 | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 18/0  | 12       | 1.0000  | 0.6378 | 1.8      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6433 | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.1     | 34/0  | 18       | 1.0000  | 0.8725 | 1.7      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.8979 | 0.2      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 25/2  | 39       | 0.5560  | 0.4154 | 7.4      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 1/7   | 1        | 0.4416  | 0.4294 | 0.5      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 38/2  | 56       | 0.9980  | 0.6370 | 7.1      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/7   | 1        | 0.9991  | 0.6433 | 0.4      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 35/2  | 53       | 0.9971  | 0.8716 | 6.9      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/7   | 1        | 0.9996  | 0.8979 | 0.4      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.01    | 94/0  | 44       | 0.4284  | 0.4152 | 5.9      | 🔴      |
| WSINDy      | Clean | RFE       |         | 1/7   | 1        | 0.4416  | 0.4294 | 0.4      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.01    | 86/0  | 41       | 0.9992  | 0.6375 | 5.1      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 1/7   | 1        | 0.9991  | 0.6433 | 0.3      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 98/0  | 46       | 0.9995  | 0.8717 | 4.6      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/7   | 1        | 0.9996  | 0.8979 | 0.3      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 0/0   | 7        | 0.0130  | 0.0320 | 13.4     | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.0128  | 0.0320 | 0.7      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.01    | 91/0  | 98       | 0.9991  | 0.6182 | 22.5     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/7   | 1        | 0.9991  | 0.6433 | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 86/0  | 93       | 0.9994  | 0.8702 | 22.0     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/7   | 1        | 0.9996  | 0.8979 | 0.9      | 🔴      |
| Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0119  | 0.0325 | 1.7      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0117  | 0.0325 | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.1     | 25/0  | 15       | 0.9990  | 0.6210 | 1.8      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.9991  | 0.6433 | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.001   | 37/0  | 20       | 0.9994  | 0.8717 | 1.9      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.9996  | 0.8979 | 0.2      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff     | E_data      | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:------------|:------------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 8  | no  | 36/2  | 54       | 1.1506      | 0.3629      | 7.8      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.4173      | 7.1      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5615      | 0.4      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 0.01    | 8  | yes | 40/2  | 58       | 0.8958      | 0.3889      | 6.5      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 12 | yes | 7/5   | 9        | 0.4412      | 0.3678      | 3.2      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5526      | 0.4      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 4  | no  | 11/2  | 19       | 0.9998      | 0.4312      | 5.9      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | no  | 7/5   | 9        | 1.0113      | 0.3361      | 3.8      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 1.0000      | 0.5778      | 0.5      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 8  | yes | 18/2  | 27       | 1.0003      | 0.4097      | 5.7      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 1.0000      | 0.3879      | 3.1      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 1.0000      | 0.5630      | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 4  | no  | 12/2  | 20       | 0.9999      | 0.5369      | 5.8      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 1.0000      | 0.3701      | 3.8      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 1.0000      | 0.6460      | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 8  | yes | 19/2  | 29       | 0.9998      | 0.4851      | 6.0      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 1.0000      | 0.4206      | 3.5      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 1.0000      | 0.6098      | 0.4      | 🔴      |
| WSINDy      | Clean | STRidge   | 1.0     | 4  | no  | 78/0  | 41       | 0.4286      | 0.4008      | 4.7      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 444689.2112 | 448874.7917 | 3.4      | 🔴      |
| WSINDy      | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5615      | 0.4      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 8  | yes | 75/0  | 37       | 0.4316      | 0.3912      | 4.6      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 663202.9872 | 668683.6885 | 2.5      | 🔴      |
| WSINDy      | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5526      | 0.4      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | no  | 87/0  | 42       | 1.0404      | 0.3256      | 4.9      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 268.7026    | 436503.6333 | 3.0      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 1.0000      | 0.5778      | 0.3      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | yes | 83/0  | 42       | 1.0027      | 0.3656      | 4.3      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 400.3246    | 653081.9031 | 2.9      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 1.0000      | 0.5630      | 0.3      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 8  | no  | 79/0  | 42       | 1.0015      | 0.4131      | 4.8      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 100.5003    | 381256.0794 | 3.3      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 1.0000      | 0.6460      | 0.5      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 12 | yes | 85/0  | 42       | 0.9999      | 0.4019      | 4.1      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 152.1176    | 594279.6068 | 3.1      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 1.0000      | 0.6098      | 0.3      | 🔴      |
| PolyJet     | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 7        | 0.0032      | 0.0138      | 12.8     | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 41.6     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 0.7      | 🟢      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0030      | 0.0133      | 11.9     | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 0.0030      | 0.0133      | 34.5     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0030      | 0.0133      | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 10.0    | 12 | no  | 66/0  | 73       | 0.7213      | 0.2564      | 17.7     | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 8  | no  | 12/1  | 18       | 0.9369      | 0.3362      | 35.4     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 0/6   | 1        | 1.0000      | 0.3548      | 0.7      | 🔴      |
| PolyJet     | 2%    | STRidge   | 1.0     | 8  | yes | 70/0  | 77       | 0.9989      | 0.3805      | 18.8     | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.5     | 12 | yes | 11/0  | 18       | 0.9039      | 0.3255      | 36.4     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 1.0000      | 0.4103      | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 8  | no  | 42/0  | 49       | 0.9995      | 0.4133      | 15.8     | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | no  | 9/3   | 13       | 1.0009      | 0.3575      | 36.1     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 1.0000      | 0.6460      | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 8  | yes | 78/0  | 85       | 1.0001      | 0.4834      | 20.0     | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 1.0000      | 0.4206      | 32.5     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 1.0000      | 0.6098      | 0.8      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0031      | 0.0139      | 1.9      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.9      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0029      | 0.0133      | 1.6      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0030      | 0.0133      | 0.7      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0030      | 0.0133      | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 8  | no  | 32/0  | 18       | 0.9388      | 0.3425      | 1.8      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 6/0   | 6        | 0.9307      | 0.3476      | 1.0      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 1.0000      | 0.3445      | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 0.001   | 12 | yes | 37/0  | 20       | 0.9279      | 0.3341      | 1.7      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 8  | yes | 6/3   | 5        | 1.0000      | 0.3981      | 0.8      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 1.0000      | 0.3931      | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 4  | no  | 28/0  | 17       | 1.0001      | 0.5369      | 1.7      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 1.0000      | 0.3701      | 0.8      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 1.0000      | 0.6460      | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 8  | yes | 27/0  | 15       | 1.0005      | 0.4893      | 1.8      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 1.0000      | 0.4206      | 0.7      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 1.0000      | 0.6098      | 0.1      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff     | E_data      | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:------------|:------------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 4  | no  | 27/2  | 38       | 0.5246      | 0.4007      | 8.0      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4416      | 0.4173      | 6.7      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.5615      | 0.4      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 0.01    | 8  | yes | 40/2  | 58       | 0.8988      | 0.3889      | 7.2      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 12 | yes | 7/5   | 9        | 0.4419      | 0.3678      | 3.6      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.5526      | 0.5      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 10.0    | 12 | no  | 12/2  | 20       | 0.9253      | 0.3312      | 5.3      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9991      | 0.3501      | 4.2      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9993      | 0.5127      | 0.4      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | yes | 24/2  | 36       | 1.0416      | 0.3663      | 5.9      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9991      | 0.3879      | 3.6      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9993      | 0.5346      | 0.4      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 10.0    | 12 | no  | 12/2  | 18       | 1.0510      | 0.3529      | 5.3      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9994      | 0.3701      | 3.7      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9996      | 0.5246      | 0.4      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | yes | 22/2  | 33       | 1.0107      | 0.4027      | 6.0      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9994      | 0.4206      | 3.5      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9996      | 0.5552      | 0.4      | 🔴      |
| WSINDy      | Clean | STRidge   | 1.0     | 4  | no  | 78/0  | 41       | 0.4292      | 0.4008      | 5.4      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 444514.6421 | 448874.7917 | 3.5      | 🔴      |
| WSINDy      | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.5615      | 0.3      | 🔴      |
| WSINDy      | Clean | STRidge   | 1.0     | 8  | yes | 88/0  | 43       | 0.4322      | 0.3903      | 4.8      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 662942.6373 | 668683.6885 | 2.9      | 🔴      |
| WSINDy      | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.5526      | 0.3      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | no  | 87/0  | 42       | 1.0758      | 0.3256      | 4.4      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 12438.5549  | 436503.6333 | 3.4      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9993      | 0.5127      | 0.4      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.001   | 8  | yes | 98/0  | 46       | 1.0037      | 0.4095      | 5.2      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 18531.7329  | 653081.9031 | 2.6      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9993      | 0.5346      | 0.3      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 12 | no  | 88/0  | 42       | 1.0577      | 0.3479      | 4.8      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 4641.1233   | 381256.0794 | 3.4      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9996      | 0.5246      | 0.4      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 8  | yes | 98/0  | 46       | 1.0016      | 0.4845      | 5.1      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 7025.4698   | 594279.6068 | 2.8      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9996      | 0.5552      | 0.3      | 🔴      |
| PolyJet     | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 7        | 0.0042      | 0.0141      | 23.3     | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 34.5     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 0.7      | 🟢      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0029      | 0.0133      | 11.8     | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 0.0030      | 0.0133      | 37.0     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0030      | 0.0133      | 0.9      | 🟢      |
| PolyJet     | 2%    | STRidge   | 10.0    | 12 | no  | 17/0  | 24       | 0.2182      | 0.1089      | 12.8     | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.5     | 12 | no  | 9/0   | 16       | 0.1883      | 0.1099      | 39.1     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 8  | no  | 0/6   | 1        | 0.9987      | 0.4262      | 0.8      | 🔴      |
| PolyJet     | 2%    | STRidge   | 1.0     | 8  | yes | 67/0  | 74       | 0.9911      | 0.3119      | 18.7     | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 10/1  | 16       | 0.7003      | 0.2122      | 38.2     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 0.9988      | 0.4060      | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 10.0    | 12 | no  | 55/0  | 62       | 0.6891      | 0.2526      | 17.2     | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.5     | 12 | no  | 13/0  | 20       | 0.6399      | 0.2611      | 37.3     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9996      | 0.5246      | 0.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | yes | 58/0  | 65       | 0.9599      | 0.3361      | 15.9     | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | yes | 17/0  | 24       | 0.9582      | 0.3414      | 36.1     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9996      | 0.5552      | 0.7      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0039      | 0.0141      | 1.9      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.8      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0029      | 0.0133      | 1.6      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0029      | 0.0133      | 0.8      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0029      | 0.0133      | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 4  | no  | 27/0  | 16       | 0.9974      | 0.3618      | 1.8      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 12 | no  | 7/0   | 5        | 0.1789      | 0.1168      | 0.7      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9984      | 0.2228      | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 10.0    | 12 | yes | 28/0  | 16       | 0.7288      | 0.2114      | 1.6      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 9/0   | 7        | 0.7271      | 0.2168      | 0.7      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9984      | 0.2773      | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.01    | 4  | no  | 34/0  | 19       | 0.9995      | 0.5225      | 2.0      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | no  | 8/0   | 6        | 0.7278      | 0.2940      | 0.8      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.9996      | 0.6460      | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 8  | yes | 36/0  | 19       | 0.9995      | 0.4641      | 1.9      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | yes | 9/2   | 6        | 0.9918      | 0.3677      | 0.6      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9988      | 0.4354      | 0.2      | 🔴      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ)     |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:------------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0171       | 15559622.6  |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0142       | 21531316.4  |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0114       | 34013918.5  |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.0227       | 15558268.0  |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0145       | 21531863.1  |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0114       | 34040234.6  |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0161       | 16584038.6  |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0136       | 22977365.0  |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.0161       | 16618686.5  |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0136       | 23050006.8  |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 18.4247      | 21474809.9  |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 11.3180      | 33213146.7  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 11.4029      | 45377815.8  |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 5.3848       | 17403912.2  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.9636       | 21997505.8  |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.8842       | 28567629.8  |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 12.1722      | 24334095.7  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 11.3156      | 34520993.6  |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 4.5355       | 18536933.6  |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.9437       | 22570644.0  |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 75.2191      | 64231689.3  |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 69.4094      | 93763895.2  |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 70.6711      | 125136634.2 |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 13.4584      | 17372921.2  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 5.4893       | 23903276.1  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 5.4493       | 31307599.2  |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 68.8768      | 71034902.5  |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 69.5782      | 97513806.4  |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 11.8995      | 18751098.2  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 5.4607       | 24719677.4  |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ks_wchaos`.

| Source                                     | Wall time | CPU                                      | Threads |
|:-------------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260730_182404_ks_wchaos_full.csv` | 486.2 min | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | not set |
