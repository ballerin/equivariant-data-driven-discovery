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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best         | Mg E_data | Std best      | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-------------------|:----------|:--------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | 🟩        | 0/0 ✓             | 0.0040     | multigraph/RFE     | 0.0123    | PolyJet/RFE   | 0/0 ✓     | 0.0041      | 0.0122     | 0.0129       |
| Clean | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.0040     | multigraph/RFE     | 0.0123    | PolyJet/RFE   | 0/0 ✓     | 0.0041      | 0.0122     | 0.0129       |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0027     | multigraph/STRidge | 0.0088    | PolyJet/MSTLS | 0/0 ✓     | 0.0030      | 0.0088     | 0.0092       |
| Clean | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0027     | multigraph/STRidge | 0.0088    | PolyJet/MSTLS | 0/0 ✓     | 0.0030      | 0.0088     | 0.0092       |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0023     | multigraph/MSTLS   | 0.0116    | PolyJet/MSTLS | 0/0 ✓     | 0.0026      | 0.0116     | 0.0097       |
| Clean | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0023     | multigraph/RFE     | 0.0116    | PolyJet/MSTLS | 0/0 ✓     | 0.0026      | 0.0116     | 0.0097       |
| 2%    | pointwise     | spectral   | —        | 19/0 ✗            | 1.0000     | —                  | 0.5297    | —             | 88/0 ✗    | 1.0000      | 0.4941     | 9249.6733    |
| 2%    | pointwise     | polynomial | —        | 37/0 ✗            | 0.9999     | —                  | 0.4488    | —             | 93/0 ✗    | 1.0000      | 0.4002     | 204.6417     |
| 2%    | weak (no IBP) | spectral   | —        | 9/0 ✗             | 0.9806     | —                  | 0.3719    | —             | 35/0 ✗    | 0.8094      | 0.2031     | 2.9497       |
| 2%    | weak (no IBP) | polynomial | 🟢        | 0/0 ✓             | 0.5015     | multigraph/MSTLS   | 0.2656    | —             | 9/0 ✗     | 0.4473      | 0.2086     | 0.4322       |
| 2%    | weak (IBP)    | spectral   | —        | 11/0 ✗            | 0.9666     | —                  | 0.3636    | —             | 43/0 ✗    | 0.9845      | 0.2872     | 3.7960       |
| 2%    | weak (IBP)    | polynomial | —        | 10/0 ✗            | 0.9670     | —                  | 0.3582    | —             | 65/0 ✗    | 0.9598      | 0.2606     | 2.4255       |
| 5%    | pointwise     | spectral   | —        | 24/0 ✗            | 1.0000     | —                  | 0.6742    | —             | 49/0 ✗    | 1.0000      | 0.6532     | 20670.3966   |
| 5%    | pointwise     | polynomial | —        | 37/0 ✗            | 1.0000     | —                  | 0.6387    | —             | 93/0 ✗    | 1.0000      | 0.6247     | 457.3188     |
| 5%    | weak (no IBP) | spectral   | —        | 18/0 ✗            | 0.9281     | —                  | 0.4005    | —             | 38/0 ✗    | 0.9975      | 0.3841     | 18.2763      |
| 5%    | weak (no IBP) | polynomial | —        | 7/0 ✗             | 0.8607     | —                  | 0.3259    | —             | 6/0 ✗     | 0.7824      | 0.3192     | 1.6768       |
| 5%    | weak (IBP)    | spectral   | —        | 27/0 ✗            | 0.9910     | —                  | 0.4391    | —             | 36/0 ✗    | 1.0000      | 0.4179     | 19.2240      |
| 5%    | weak (IBP)    | polynomial | —        | 7/0 ✗             | 0.9886     | —                  | 0.3949    | —             | 53/0 ✗    | 1.0004      | 0.4299     | 6.1829       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                | Best conf.                   | fp/fn  | E_coeff   | Time (s) | vs. best |
|:----------------------|:-----------------------------|:-------|----------:|---------:|---------:|
| Multigraphs + RFE     | polynomial, weak (IBP), m=8  | 0/0 ✓  | 0.0023    | 0.1      | —        |
| Multigraphs + MSTLS   | polynomial, weak (IBP), m=8  | 0/0 ✓  | 0.0023    | 0.3      | +0.0%    |
| Multigraphs + STRidge | polynomial, weak (IBP), m=8  | 0/0 ✓  | 0.0023    | 0.2      | +0.5%    |
| PolyJet + MSTLS       | polynomial, weak (IBP), m=8  | 0/0 ✓  | 0.0026    | 52.1     | +12.1%   |
| PolyJet + RFE         | polynomial, weak (IBP), m=8  | 0/0 ✓  | 0.0026    | 0.5      | +12.1%   |
| PolyJet + STRidge     | polynomial, weak (IBP), m=8  | 0/0 ✓  | 0.0027    | 5.0      | +17.3%   |
| WSINDy + MSTLS        | polynomial, weak (IBP), m=8  | 44/0 ✗ | 1985.5275 | 1.2      | ✗        |
| WSINDy + STRidge      | spectral, weak (IBP), m=8    | 48/0 ✗ | 0.6248    | 0.4      | ✗        |
| PDE-FIND + MSTLS      | polynomial, weak (IBP), m=12 | 7/2 ✗  | 0.6539    | 3.5      | ✗        |
| PDE-FIND + STRidge    | polynomial, pointwise        | 27/2 ✗ | 0.8327    | 8.7      | ✗        |
| PDE-FIND + RFE        | spectral, weak (IBP), m=12   | 14/4 ✗ | 4.8515    | 0.1      | ✗        |
| WSINDy + RFE          | polynomial, weak (IBP), m=12 | 4/5 ✗  | 0.7367    | 0.1      | ✗        |

### 2% noise

| Method                | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.5015  | 0.3      | —        |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=12 | 9/0 ✗  | 0.4473  | 49.3     | ✗        |
| Multigraphs + STRidge | spectral, pointwise             | 19/0 ✗ | 1.0000  | 2.3      | ✗        |
| PolyJet + STRidge     | spectral, weak (IBP), m=8       | 48/0 ✗ | 0.9995  | 10.1     | ✗        |
| WSINDy + STRidge      | spectral, weak (IBP), m=12      | 64/0 ✗ | 0.9808  | 0.3      | ✗        |
| WSINDy + MSTLS        | spectral, weak (IBP), m=12      | 21/2 ✗ | 41.4239 | 0.8      | ✗        |
| PDE-FIND + STRidge    | polynomial, weak (IBP), m=8     | 22/2 ✗ | 1.0001  | 0.9      | ✗        |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12   | 23/2 ✗ | 0.9481  | 4.5      | ✗        |
| Multigraphs + RFE     | spectral, weak (no IBP), m=8    | 0/5 ✗  | 1.0000  | 0.1      | ✗        |
| WSINDy + RFE          | spectral, weak (IBP), m=12      | 4/5 ✗  | 1.0011  | 0.1      | ✗        |
| PolyJet + RFE         | spectral, weak (no IBP), m=12   | 1/7 ✗  | 1.0000  | 0.4      | ✗        |
| PDE-FIND + RFE        | spectral, weak (IBP), m=12      | 1/7 ✗  | 1.0000  | 0.1      | ✗        |

### 5% noise

| Method                | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=12 | 6/0 ✗  | 0.7824  | 37.9     | ✗        |
| Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 7/0 ✗  | 0.8607  | 0.2      | ✗        |
| Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 18/0 ✗ | 0.8779  | 0.1      | ✗        |
| PolyJet + STRidge     | polynomial, weak (no IBP), m=12 | 41/0 ✗ | 0.8247  | 5.6      | ✗        |
| WSINDy + STRidge      | spectral, weak (IBP), m=8       | 43/0 ✗ | 1.0000  | 0.3      | ✗        |
| PDE-FIND + STRidge    | polynomial, weak (no IBP), m=8  | 20/2 ✗ | 1.0027  | 0.9      | ✗        |
| PDE-FIND + MSTLS      | polynomial, weak (no IBP), m=12 | 22/2 ✗ | 0.9549  | 2.9      | ✗        |
| WSINDy + MSTLS        | spectral, weak (IBP), m=12      | 29/2 ✗ | 47.4136 | 0.8      | ✗        |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=8  | 0/5 ✗  | 1.0000  | 0.0      | ✗        |
| WSINDy + RFE          | spectral, weak (no IBP), m=12   | 4/5 ✗  | 1.0012  | 0.1      | ✗        |
| PolyJet + RFE         | spectral, weak (no IBP), m=8    | 0/6 ✗  | 1.0000  | 0.4      | ✗        |
| PDE-FIND + RFE        | spectral, weak (IBP), m=12      | 1/7 ✗  | 1.0000  | 0.1      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 27/2  | 39       | 0.8327  | 0.4141 | 7.5      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 1/7   | 1        | 0.6070  | 0.6047 | 0.4      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.01    | 37/2  | 56       | 1.0000  | 0.4924 | 8.6      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 0.9      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 31/2  | 50       | 1.0000  | 0.6544 | 7.3      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 3/7   | 5        | 1.0000  | 0.7065 | 0.5      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 62/0  | 33       | 0.4222  | 0.4158 | 12.8     | 🔴      |
| WSINDy      | Clean | RFE       |         | 1/7   | 1        | 0.6070  | 0.6047 | 1.2      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.01    | 88/0  | 43       | 1.0000  | 0.4941 | 5.6      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 0.3      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.01    | 97/0  | 45       | 1.0000  | 0.6555 | 5.3      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.7219 | 0.3      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 0/0   | 7        | 0.0041  | 0.0122 | 33.1     | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.0041  | 0.0122 | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.01    | 88/0  | 95       | 1.0000  | 0.4713 | 27.0     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 49/0  | 56       | 1.0000  | 0.6532 | 19.2     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.7219 | 0.7      | 🔴      |
| Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0040  | 0.0123 | 1.6      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.1     | 19/0  | 12       | 1.0000  | 0.5297 | 2.3      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.1     | 24/0  | 14       | 1.0000  | 0.6742 | 1.7      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.7219 | 0.1      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 27/2  | 39       | 0.8327  | 0.4141 | 8.7      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 1/7   | 1        | 0.6070  | 0.6047 | 0.9      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.001   | 41/2  | 60       | 0.9995  | 0.4709 | 9.2      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 40/2  | 59       | 0.9998  | 0.6237 | 8.2      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.7219 | 0.4      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 62/0  | 33       | 0.4222  | 0.4158 | 13.8     | 🔴      |
| WSINDy      | Clean | RFE       |         | 1/7   | 1        | 0.6070  | 0.6047 | 0.7      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.001   | 98/0  | 46       | 0.9999  | 0.4722 | 5.6      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 0.3      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 93/0  | 45       | 1.0000  | 0.6247 | 5.5      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.7219 | 0.3      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 0/0   | 7        | 0.0041  | 0.0122 | 17.2     | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.0041  | 0.0122 | 0.7      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.01    | 93/0  | 100      | 1.0000  | 0.4002 | 27.8     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 1.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.001   | 97/0  | 104      | 1.0000  | 0.6009 | 24.2     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.7219 | 0.9      | 🔴      |
| Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0040  | 0.0123 | 2.1      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.001   | 37/0  | 20       | 0.9999  | 0.4488 | 2.1      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.001   | 37/0  | 20       | 1.0000  | 0.6387 | 2.2      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.7219 | 0.1      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff   | E_data    | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:----------|:----------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 4  | no  | 27/2  | 41       | 1.1391    | 0.4127    | 8.5      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 4  | no  | 27/2  | 33       | 0.7140    | 0.4195    | 9.9      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | no  | 16/4  | 26       | 5.4960    | 0.3547    | 0.1      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 1.0     | 12 | yes | 29/2  | 44       | 2.8852    | 0.3315    | 0.7      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | yes | 7/2   | 12       | 0.6539    | 0.4650    | 3.6      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | yes | 14/4  | 23       | 4.8515    | 0.3720    | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | no  | 39/2  | 57       | 0.9380    | 0.3454    | 1.3      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | no  | 23/2  | 28       | 0.9481    | 0.3543    | 4.5      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 1.0000    | 0.5668    | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 12 | yes | 39/2  | 57       | 0.9885    | 0.3630    | 0.7      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | yes | 24/2  | 29       | 0.9886    | 0.3700    | 5.0      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000    | 0.5779    | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 12 | no  | 41/2  | 60       | 0.9789    | 0.3599    | 1.2      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.75    | 12 | no  | 27/2  | 33       | 0.9867    | 0.3650    | 3.2      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | no  | 1/7   | 1        | 1.0000    | 0.5955    | 0.2      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 12 | yes | 41/2  | 60       | 0.9954    | 0.3806    | 0.9      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.75    | 8  | yes | 17/3  | 21       | 1.0002    | 0.4384    | 5.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000    | 0.5812    | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 10.0    | 12 | no  | 51/0  | 30       | 0.7601    | 0.3146    | 0.4      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.75    | 8  | no  | 29/1  | 20       | 566.3077  | 622.3344  | 1.9      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 4/5   | 4        | 0.7577    | 0.4809    | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 10.0    | 8  | yes | 48/0  | 31       | 0.6248    | 0.3608    | 0.4      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.75    | 8  | yes | 44/0  | 30       | 1985.5287 | 2170.3700 | 1.3      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 4/5   | 4        | 0.7367    | 0.4959    | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 4  | no  | 89/0  | 42       | 1.0000    | 0.4289    | 5.8      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.75    | 12 | no  | 44/2  | 28       | 98.7468   | 1879.8458 | 0.9      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 4/5   | 4        | 1.0014    | 0.4838    | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 1.0     | 12 | yes | 64/0  | 34       | 0.9808    | 0.3507    | 0.3      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.75    | 12 | yes | 21/2  | 17       | 41.4239   | 837.4022  | 0.8      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 4/5   | 4        | 1.0011    | 0.5081    | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 8  | no  | 60/0  | 37       | 0.9953    | 0.4058    | 0.4      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.75    | 12 | no  | 17/4  | 14       | 8.3290    | 763.0814  | 0.7      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 4/5   | 4        | 1.0012    | 0.4954    | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 8  | yes | 43/0  | 29       | 1.0000    | 0.4323    | 0.3      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.75    | 12 | yes | 29/2  | 22       | 47.4136   | 2485.6412 | 0.8      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000    | 0.5812    | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 10.0     | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 56.9     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 0.4      | 🟢      |
| PolyJet     | Clean | STRidge   | 1.0     | 8  | yes | 0/0   | 7        | 0.0027    | 0.0116    | 5.0      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 7        | 0.0026    | 0.0116    | 49.2     | 🟢      |
| PolyJet     | Clean | RFE       |         | 8  | yes | 0/0   | 7        | 0.0026    | 0.0116    | 0.6      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 8  | no  | 58/0  | 65       | 0.9938    | 0.3054    | 11.7     | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | no  | 35/0  | 42       | 0.8094    | 0.2031    | 39.1     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 1.0000    | 0.6415    | 0.4      | 🔴      |
| PolyJet     | 2%    | STRidge   | 1.0     | 8  | yes | 48/0  | 55       | 0.9995    | 0.3354    | 10.1     | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 43/0  | 50       | 0.9845    | 0.2872    | 42.4     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000    | 0.5779    | 0.5      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 4  | no  | 81/0  | 88       | 1.0000    | 0.4365    | 24.4     | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 8  | no  | 38/0  | 45       | 0.9975    | 0.3841    | 40.6     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 8  | no  | 0/6   | 1        | 1.0000    | 0.6513    | 0.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | yes | 68/0  | 75       | 0.9915    | 0.3473    | 8.9      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 8  | yes | 36/0  | 43       | 1.0000    | 0.4179    | 38.7     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 1.0000    | 0.6067    | 0.4      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.3      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.0      | 🟢      |
| Multigraphs | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 3        | 0.0023    | 0.0116    | 0.2      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 3        | 0.0023    | 0.0116    | 0.3      | 🟢      |
| Multigraphs | Clean | RFE       |         | 8  | yes | 0/0   | 3        | 0.0023    | 0.0116    | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 10.0    | 12 | no  | 25/0  | 14       | 0.5878    | 0.2571    | 0.2      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 9/0   | 8        | 0.9806    | 0.3719    | 0.4      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 8  | no  | 0/5   | 1        | 1.0000    | 0.5441    | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 1.0     | 8  | yes | 31/0  | 17       | 0.9995    | 0.3959    | 0.2      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 11/0  | 9        | 0.9666    | 0.3636    | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 8  | yes | 0/5   | 1        | 1.0000    | 0.5679    | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 31/0  | 17       | 0.9245    | 0.3920    | 0.1      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | no  | 18/0  | 12       | 0.9281    | 0.4005    | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 1.0000    | 0.6154    | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 8  | yes | 31/0  | 17       | 1.0001    | 0.5015    | 0.1      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | yes | 27/0  | 15       | 0.9910    | 0.4391    | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 1.0000    | 0.6067    | 0.1      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff   | E_data    | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:----------|:----------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 4  | no  | 27/2  | 41       | 1.1391    | 0.4127    | 8.3      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 4  | no  | 27/2  | 33       | 0.7140    | 0.4195    | 11.3     | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | no  | 16/4  | 26       | 5.4961    | 0.3547    | 0.1      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 1.0     | 12 | yes | 29/2  | 44       | 2.8852    | 0.3315    | 0.6      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | yes | 7/2   | 12       | 0.6539    | 0.4650    | 3.5      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | yes | 14/4  | 23       | 4.8515    | 0.3720    | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.001   | 12 | no  | 41/2  | 60       | 1.0005    | 0.3323    | 1.7      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | no  | 25/2  | 30       | 0.9324    | 0.3425    | 3.5      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | no  | 1/7   | 1        | 1.0000    | 0.5924    | 0.3      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 8  | yes | 22/2  | 34       | 1.0001    | 0.4159    | 0.9      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | yes | 28/2  | 33       | 0.9821    | 0.3636    | 2.9      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 1.0000    | 0.5996    | 0.3      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 8  | no  | 20/2  | 30       | 1.0027    | 0.4207    | 0.9      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.75    | 12 | no  | 22/2  | 27       | 0.9549    | 0.3632    | 2.9      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | no  | 1/7   | 1        | 1.0000    | 0.5955    | 0.2      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | yes | 25/2  | 37       | 0.9938    | 0.3811    | 0.5      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.75    | 12 | yes | 15/3  | 19       | 0.9920    | 0.3921    | 2.3      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000    | 0.5812    | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 10.0    | 8  | no  | 51/0  | 32       | 0.6387    | 0.3554    | 0.4      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.75    | 8  | no  | 29/1  | 20       | 566.3111  | 622.3367  | 1.3      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 4/5   | 4        | 0.7577    | 0.4809    | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 10.0    | 8  | yes | 48/0  | 31       | 0.6248    | 0.3608    | 0.3      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.75    | 8  | yes | 44/0  | 30       | 1985.5275 | 2170.3680 | 1.2      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 4/5   | 4        | 0.7367    | 0.4959    | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 1.0     | 12 | no  | 90/0  | 45       | 0.9239    | 0.3220    | 0.3      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.75    | 12 | no  | 44/2  | 28       | 231.2112  | 1877.6884 | 0.9      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 4/5   | 4        | 1.0055    | 0.4838    | 0.2      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.001   | 12 | yes | 98/0  | 46       | 0.9836    | 0.3485    | 0.3      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.75    | 12 | yes | 11/3  | 12       | 6.8674    | 770.7203  | 0.7      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 4/5   | 4        | 1.0044    | 0.5081    | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 4  | no  | 54/0  | 30       | 1.0002    | 0.4502    | 4.6      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.75    | 12 | no  | 17/4  | 14       | 44.0980   | 763.2070  | 0.8      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 4/5   | 4        | 1.0049    | 0.4953    | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 8  | yes | 53/0  | 30       | 1.0004    | 0.4299    | 0.4      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.75    | 12 | yes | 29/2  | 22       | 98.4861   | 2526.6223 | 0.8      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000    | 0.5812    | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 4.7      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 51.1     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 0.5      | 🟢      |
| PolyJet     | Clean | STRidge   | 1.0     | 8  | yes | 0/0   | 7        | 0.0027    | 0.0116    | 5.0      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 7        | 0.0026    | 0.0116    | 52.1     | 🟢      |
| PolyJet     | Clean | RFE       |         | 8  | yes | 0/0   | 7        | 0.0026    | 0.0116    | 0.5      | 🟢      |
| PolyJet     | 2%    | STRidge   | 10.0    | 12 | no  | 63/0  | 70       | 0.5965    | 0.1642    | 11.2     | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | no  | 9/0   | 16       | 0.4473    | 0.2086    | 49.3     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 1.0000    | 0.5668    | 0.4      | 🔴      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | yes | 77/0  | 84       | 0.9637    | 0.2556    | 9.7      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 65/0  | 72       | 0.9598    | 0.2606    | 38.8     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000    | 0.5779    | 0.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 10.0    | 12 | no  | 41/0  | 48       | 0.8247    | 0.2458    | 5.6      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.5     | 12 | no  | 6/0   | 13       | 0.7824    | 0.3192    | 37.9     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 1.0000    | 0.5687    | 0.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | yes | 66/0  | 73       | 0.9882    | 0.2988    | 10.2     | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | yes | 8/2   | 13       | 0.9887    | 0.3693    | 34.4     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000    | 0.5812    | 0.4      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.3      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 3        | 0.0023    | 0.0116    | 0.2      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 3        | 0.0023    | 0.0116    | 0.3      | 🟢      |
| Multigraphs | Clean | RFE       |         | 8  | yes | 0/0   | 3        | 0.0023    | 0.0116    | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 10.0    | 12 | no  | 20/0  | 12       | 0.6173    | 0.2255    | 0.1      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.5015    | 0.2656    | 0.3      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 1.0000    | 0.5668    | 0.0      | 🔴      |
| Multigraphs | 2%    | STRidge   | 1.0     | 8  | yes | 31/0  | 16       | 0.9991    | 0.3810    | 0.2      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 10/0  | 8        | 0.9670    | 0.3582    | 0.3      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000    | 0.5779    | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 18/0  | 12       | 0.8779    | 0.3078    | 0.1      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | no  | 7/0   | 7        | 0.8607    | 0.3259    | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 8  | no  | 0/5   | 1        | 1.0000    | 0.4886    | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 8  | yes | 29/0  | 17       | 1.0004    | 0.4129    | 0.2      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | yes | 7/0   | 7        | 0.9886    | 0.3949    | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000    | 0.5812    | 0.0      | 🔴      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ)   |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:----------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0126       | 105502.8  |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0114       | 114784.0  |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0092       | 125064.5  |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.0126       | 105504.4  |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0114       | 114785.5  |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0092       | 125064.6  |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0118       | 108826.2  |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0097       | 112577.6  |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.0118       | 108828.0  |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0097       | 112580.0  |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 88.9904      | 336798.9  |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 3.9765       | 471736.7  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 2.9497       | 693750.9  |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 33.2300      | 55772.9   |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 2.5916       | 49341.4   |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.4322       | 68936.2   |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 28.9149      | 352903.8  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 3.7960       | 507486.5  |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 28.7695      | 39179.9   |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 2.4255       | 55894.4   |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 220.7844     | 2029412.3 |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 18.7254      | 3418856.0 |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 18.2763      | 4850904.7 |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 82.2495      | 142599.5  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 6.5727       | 58296.0   |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 1.6768       | 79176.8   |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 73.6059      | 2400411.1 |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 19.2240      | 3655517.7 |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 71.4408      | 78089.2   |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 6.1829       | 59364.4   |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ks_cell`.

| Source                                   | Wall time | CPU                                      | Threads |
|:-----------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260730_115725_ks_cell_full.csv` | 386.7 min | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | not set |
