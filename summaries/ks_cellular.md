# Kuramoto-Sivashinsky equation (cellular, $L=4\pi$)

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

KS contains the gradient norm $\|\nabla f\|^2 = f_x^2 + f_y^2$ which neither PDE-FIND nor WSINDy can represent as a sparse column (PDE-FIND lacks products of derivatives; in WSINDy, $f_x^2$ is locked with $f\cdot f_{xx}$ inside $\partial_{xx}(f^2)$).  Cellular regime ($L=4\pi$).

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
| Clean | pointwise     | spectral   | 🟩        | 0/0 ✓             | 0.0040     | multigraph/RFE     | 0.0123    | PolyJet/RFE   | 0/0 ✓     | 0.0041      | 0.0122     | —            |
| Clean | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.0040     | multigraph/RFE     | 0.0123    | PolyJet/RFE   | 0/0 ✓     | 0.0041      | 0.0122     | —            |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0027     | multigraph/STRidge | 0.0088    | PolyJet/MSTLS | 0/0 ✓     | 0.0030      | 0.0088     | —            |
| Clean | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0027     | multigraph/STRidge | 0.0088    | PolyJet/RFE   | 0/0 ✓     | 0.0030      | 0.0088     | —            |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0023     | multigraph/MSTLS   | 0.0116    | PolyJet/MSTLS | 0/0 ✓     | 0.0026      | 0.0116     | —            |
| Clean | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0023     | multigraph/MSTLS   | 0.0116    | PolyJet/MSTLS | 0/0 ✓     | 0.0026      | 0.0116     | —            |
| 2%    | pointwise     | spectral   | —        | 37/0 ✗            | 0.4675     | —                  | 0.5297    | —             | 60/0 ✗    | 0.3941      | 0.4714     | —            |
| 2%    | pointwise     | polynomial | —        | 34/0 ✗            | 0.3850     | —                  | 0.4488    | —             | 98/0 ✗    | 0.3351      | 0.4002     | —            |
| 2%    | weak (no IBP) | spectral   | —        | 9/0 ✗             | 0.4585     | —                  | 0.3719    | —             | 35/0 ✗    | 0.4372      | 0.2031     | —            |
| 2%    | weak (no IBP) | polynomial | 🟢        | 0/0 ✓             | 0.2711     | multigraph/STRidge | 0.2656    | —             | 9/0 ✗     | 0.3141      | 0.2086     | —            |
| 2%    | weak (IBP)    | spectral   | —        | 11/0 ✗            | 0.5216     | —                  | 0.3636    | —             | 43/0 ✗    | 0.4585      | 0.2872     | —            |
| 2%    | weak (IBP)    | polynomial | —        | 10/0 ✗            | 0.5112     | —                  | 0.3582    | —             | 61/0 ✗    | 0.5088      | 0.2563     | —            |
| 5%    | pointwise     | spectral   | —        | 32/0 ✗            | 0.5283     | —                  | 0.6742    | —             | 38/0 ✗    | 0.4714      | 0.6538     | —            |
| 5%    | pointwise     | polynomial | —        | 34/0 ✗            | 0.4307     | —                  | 0.6387    | —             | 91/0 ✗    | 0.3660      | 0.6009     | —            |
| 5%    | weak (no IBP) | spectral   | —        | 18/0 ✗            | 0.9587     | —                  | 0.4005    | —             | 38/0 ✗    | 0.4873      | 0.3841     | —            |
| 5%    | weak (no IBP) | polynomial | —        | 7/0 ✗             | 0.4445     | —                  | 0.3947    | —             | 6/0 ✗     | 0.4636      | 0.3192     | —            |
| 5%    | weak (IBP)    | spectral   | —        | 27/0 ✗            | 0.5433     | —                  | 0.5017    | —             | 36/0 ✗    | 0.4477      | 0.4179     | —            |
| 5%    | weak (IBP)    | polynomial | —        | 7/0 ✗             | 0.4906     | —                  | 0.3949    | —             | 53/0 ✗    | 0.4551      | 0.4299     | —            |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                | Best conf.                   | fp/fn  | E_coeff   | Time (s) | vs. best |
|:----------------------|:-----------------------------|:-------|----------:|---------:|---------:|
| Multigraphs + MSTLS   | polynomial, weak (IBP), m=8  | 0/0 ✓  | 0.0023    | 0.2      | —        |
| Multigraphs + RFE     | polynomial, weak (IBP), m=8  | 0/0 ✓  | 0.0023    | 0.0      | +0.0%    |
| Multigraphs + STRidge | polynomial, weak (IBP), m=8  | 0/0 ✓  | 0.0023    | 0.2      | +0.5%    |
| PolyJet + MSTLS       | polynomial, weak (IBP), m=8  | 0/0 ✓  | 0.0026    | 2.8      | +12.1%   |
| PolyJet + RFE         | polynomial, weak (IBP), m=8  | 0/0 ✓  | 0.0026    | 0.1      | +12.1%   |
| PolyJet + STRidge     | polynomial, weak (IBP), m=8  | 0/0 ✓  | 0.0027    | 1.1      | +17.3%   |
| WSINDy + STRidge      | spectral, weak (IBP), m=12   | 42/0 ✗ | 0.7245    | 0.5      | ✗        |
| WSINDy + MSTLS        | polynomial, weak (IBP), m=8  | 44/0 ✗ | 1985.5275 | 0.7      | ✗        |
| PDE-FIND + MSTLS      | polynomial, weak (IBP), m=12 | 7/2 ✗  | 0.6539    | 0.8      | ✗        |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=8 | 31/2 ✗ | 0.6789    | 1.2      | ✗        |
| PDE-FIND + RFE        | spectral, weak (IBP), m=12   | 20/4 ✗ | 0.9472    | 0.1      | ✗        |
| WSINDy + RFE          | polynomial, weak (IBP), m=12 | 4/5 ✗  | 0.7367    | 0.1      | ✗        |

### 2% noise

| Method                | Best conf.                      | fp/fn  | E_coeff  | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|---------:|---------:|---------:|
| Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.2711   | 0.2      | —        |
| Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.2712   | 0.2      | +0.1%    |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=12 | 9/0 ✗  | 0.3141   | 2.1      | ✗        |
| PolyJet + STRidge     | spectral, weak (no IBP), m=8    | 38/0 ✗ | 0.3767   | 1.9      | ✗        |
| WSINDy + STRidge      | polynomial, weak (no IBP), m=4  | 39/0 ✗ | 0.4346   | 1.0      | ✗        |
| WSINDy + MSTLS        | spectral, weak (IBP), m=12      | 21/2 ✗ | 662.2014 | 0.5      | ✗        |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12   | 23/2 ✗ | 0.6171   | 0.8      | ✗        |
| PDE-FIND + STRidge    | polynomial, pointwise           | 38/2 ✗ | 0.4278   | 2.0      | ✗        |
| Multigraphs + RFE     | spectral, weak (no IBP), m=8    | 0/5 ✗  | 0.7796   | 0.1      | ✗        |
| WSINDy + RFE          | polynomial, weak (IBP), m=12    | 4/5 ✗  | 0.7045   | 0.1      | ✗        |
| PDE-FIND + RFE        | polynomial, pointwise           | 1/7 ✗  | 0.6070   | 0.1      | ✗        |
| PolyJet + RFE         | polynomial, pointwise           | 1/7 ✗  | 0.6070   | 0.1      | ✗        |

### 5% noise

| Method                | Best conf.                      | fp/fn  | E_coeff   | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|----------:|---------:|---------:|
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=12 | 6/0 ✗  | 0.4636    | 1.9      | ✗        |
| Multigraphs + MSTLS   | polynomial, weak (no IBP), m=8  | 7/0 ✗  | 0.4445    | 0.2      | ✗        |
| Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 13/0 ✗ | 0.4499    | 0.2      | ✗        |
| PolyJet + STRidge     | spectral, pointwise             | 38/0 ✗ | 0.4714    | 2.7      | ✗        |
| WSINDy + STRidge      | polynomial, weak (IBP), m=8     | 53/0 ✗ | 0.4551    | 0.7      | ✗        |
| PDE-FIND + MSTLS      | polynomial, weak (no IBP), m=12 | 22/2 ✗ | 0.6175    | 0.7      | ✗        |
| PDE-FIND + STRidge    | polynomial, weak (IBP), m=12    | 27/2 ✗ | 0.5556    | 0.7      | ✗        |
| WSINDy + MSTLS        | spectral, weak (IBP), m=12      | 29/2 ✗ | 2741.0990 | 0.5      | ✗        |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=8  | 0/5 ✗  | 0.6489    | 0.1      | ✗        |
| WSINDy + RFE          | spectral, weak (no IBP), m=12   | 4/5 ✗  | 0.7226    | 0.1      | ✗        |
| PolyJet + RFE         | spectral, weak (no IBP), m=8    | 0/6 ✗  | 0.9571    | 0.1      | ✗        |
| PDE-FIND + RFE        | polynomial, pointwise           | 1/7 ✗  | 0.6070    | 0.1      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 45/2  | 50       | 0.4166  | 0.4099 | 1.9      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 1/7   | 1        | 0.6070  | 0.6047 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.001   | 55/2  | 60       | 0.4349  | 0.4922 | 2.2      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 42/2  | 47       | 0.4732  | 0.6544 | 1.7      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 5/7   | 5        | 0.5678  | 0.7065 | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.01    | 62/0  | 33       | 0.4222  | 0.4158 | 1.2      | 🔴      |
| WSINDy      | Clean | RFE       |         | 1/7   | 1        | 0.6070  | 0.6047 | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.01    | 79/0  | 39       | 0.4365  | 0.4941 | 1.3      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.01    | 85/0  | 42       | 0.4760  | 0.6557 | 1.4      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0041  | 0.0122 | 1.4      | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.0041  | 0.0122 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.01    | 60/0  | 67       | 0.3941  | 0.4714 | 4.4      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 38/0  | 45       | 0.4714  | 0.6538 | 2.7      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0040  | 0.0123 | 0.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.001   | 37/0  | 20       | 0.4675  | 0.5297 | 0.4      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.01    | 32/0  | 18       | 0.5283  | 0.6742 | 0.4      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 45/2  | 50       | 0.4166  | 0.4099 | 1.8      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 1/7   | 1        | 0.6070  | 0.6047 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.01    | 38/2  | 43       | 0.4278  | 0.4712 | 2.0      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 53/2  | 58       | 0.4306  | 0.6237 | 2.3      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.01    | 62/0  | 33       | 0.4222  | 0.4158 | 1.2      | 🔴      |
| WSINDy      | Clean | RFE       |         | 1/7   | 1        | 0.6070  | 0.6047 | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.001   | 98/1  | 45       | 0.4290  | 0.4722 | 1.4      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 72/1  | 38       | 0.4318  | 0.6248 | 1.2      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0041  | 0.0122 | 1.4      | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.0041  | 0.0122 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.001   | 98/0  | 105      | 0.3351  | 0.4002 | 5.5      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.001   | 91/0  | 98       | 0.3660  | 0.6009 | 5.9      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0040  | 0.0123 | 0.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.01    | 34/0  | 17       | 0.3850  | 0.4488 | 0.4      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.001   | 34/0  | 19       | 0.4307  | 0.6387 | 0.4      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff   | E_data    | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:----------|:----------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 10.0    | 8  | no  | 31/2  | 36       | 0.6789    | 0.3706    | 1.2      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 4  | no  | 28/2  | 33       | 0.4584    | 0.4195    | 1.6      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | no  | 23/4  | 26       | 0.9619    | 0.3547    | 0.1      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 1.0     | 8  | yes | 37/2  | 42       | 0.6819    | 0.3662    | 1.5      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | yes | 7/2   | 12       | 0.6539    | 0.4650    | 0.9      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | yes | 20/4  | 23       | 0.9472    | 0.3720    | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 12 | no  | 53/2  | 58       | 0.6193    | 0.3454    | 1.1      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | no  | 23/2  | 28       | 0.6171    | 0.3543    | 0.8      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6051    | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.01    | 8  | yes | 47/2  | 52       | 0.4719    | 0.4164    | 1.5      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | yes | 24/2  | 29       | 0.5836    | 0.3700    | 0.8      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6085    | 0.5996    | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 4  | no  | 42/2  | 47       | 0.4355    | 0.4546    | 1.7      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.75    | 12 | no  | 28/2  | 33       | 0.5562    | 0.3650    | 0.8      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6154    | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.01    | 8  | yes | 42/2  | 47       | 0.4556    | 0.4363    | 1.4      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.75    | 8  | yes | 17/3  | 21       | 0.4499    | 0.4384    | 1.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6085    | 0.6067    | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 10.0    | 12 | no  | 59/0  | 33       | 0.7731    | 0.3104    | 0.6      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.75    | 8  | no  | 29/1  | 20       | 566.3077  | 622.3344  | 0.7      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 4/5   | 4        | 0.7577    | 0.4809    | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 10.0    | 12 | yes | 42/0  | 27       | 0.7245    | 0.3266    | 0.5      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.75    | 8  | yes | 44/0  | 30       | 1985.5287 | 2170.3700 | 0.7      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 4/5   | 4        | 0.7367    | 0.4959    | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 8  | no  | 77/0  | 41       | 0.5043    | 0.3894    | 0.9      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.75    | 12 | no  | 44/2  | 28       | 1625.2772 | 1879.8458 | 0.5      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 4/5   | 4        | 0.7504    | 0.4838    | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 8  | yes | 75/0  | 40       | 0.4711    | 0.4088    | 0.9      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.75    | 12 | yes | 21/2  | 17       | 662.2014  | 837.4022  | 0.5      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 4/5   | 4        | 0.7047    | 0.5081    | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 4  | no  | 98/0  | 46       | 0.4337    | 0.4533    | 1.3      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.75    | 12 | no  | 17/4  | 14       | 654.7380  | 763.0814  | 0.5      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 4/5   | 4        | 0.7226    | 0.4954    | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 8  | yes | 65/0  | 38       | 0.4543    | 0.4314    | 0.9      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.75    | 12 | yes | 29/2  | 22       | 2741.0990 | 2485.6412 | 0.5      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6085    | 0.6067    | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 8  | no  | 0/0   | 7        | 0.0037    | 0.0108    | 1.1      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 1.9      | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 0.1      | 🟢      |
| PolyJet     | Clean | STRidge   | 1.0     | 8  | yes | 0/0   | 7        | 0.0027    | 0.0116    | 1.0      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 7        | 0.0026    | 0.0116    | 2.9      | 🟢      |
| PolyJet     | Clean | RFE       |         | 8  | yes | 0/0   | 7        | 0.0026    | 0.0116    | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 8  | no  | 38/0  | 45       | 0.3767    | 0.3079    | 1.9      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | no  | 35/0  | 42       | 0.4372    | 0.2031    | 2.0      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6051    | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 1.0     | 8  | yes | 48/0  | 55       | 0.3600    | 0.3354    | 2.1      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 43/0  | 50       | 0.4585    | 0.2872    | 1.9      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6085    | 0.5996    | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.01    | 4  | no  | 80/0  | 87       | 0.4181    | 0.4365    | 5.2      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 8  | no  | 38/0  | 45       | 0.4873    | 0.3841    | 2.8      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 8  | no  | 0/6   | 1        | 0.9571    | 0.6513    | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 12 | yes | 78/0  | 85       | 0.5573    | 0.3471    | 2.4      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 8  | yes | 36/0  | 43       | 0.4477    | 0.4179    | 2.6      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6085    | 0.6067    | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.2      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 1.0     | 8  | yes | 0/0   | 3        | 0.0023    | 0.0116    | 0.2      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 3        | 0.0023    | 0.0116    | 0.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 8  | yes | 0/0   | 3        | 0.0023    | 0.0116    | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 10.0    | 12 | no  | 11/0  | 7        | 0.3770    | 0.2649    | 0.2      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 9/0   | 8        | 0.4585    | 0.3719    | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.7796    | 0.5441    | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 1.0     | 8  | yes | 27/0  | 15       | 0.4348    | 0.3959    | 0.3      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 11/0  | 9        | 0.5216    | 0.3636    | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 8  | yes | 0/5   | 1        | 0.7857    | 0.5679    | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 28/0  | 16       | 0.9366    | 0.3920    | 0.2      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | no  | 18/0  | 12       | 0.9587    | 0.4005    | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6154    | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 12 | yes | 31/0  | 17       | 0.6273    | 0.4393    | 0.3      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 8  | yes | 27/0  | 15       | 0.5433    | 0.5017    | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6085    | 0.6067    | 0.1      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff   | E_data    | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:----------|:----------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 10.0    | 8  | no  | 31/2  | 36       | 0.6789    | 0.3706    | 1.2      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 4  | no  | 28/2  | 33       | 0.4584    | 0.4195    | 1.6      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | no  | 23/4  | 26       | 0.9619    | 0.3547    | 0.1      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 1.0     | 8  | yes | 37/2  | 42       | 0.6819    | 0.3662    | 1.5      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | yes | 7/2   | 12       | 0.6539    | 0.4650    | 0.8      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | yes | 20/4  | 23       | 0.9472    | 0.3720    | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | no  | 43/2  | 48       | 0.6880    | 0.3395    | 1.0      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | no  | 25/2  | 30       | 0.6531    | 0.3425    | 0.8      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6051    | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 8  | yes | 40/2  | 45       | 0.4700    | 0.4149    | 1.4      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | yes | 28/2  | 33       | 0.6195    | 0.3636    | 0.8      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6085    | 0.5996    | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 8  | no  | 34/2  | 39       | 0.4654    | 0.4194    | 1.2      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.75    | 12 | no  | 22/2  | 27       | 0.6175    | 0.3632    | 0.7      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6154    | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | yes | 27/2  | 32       | 0.5556    | 0.3815    | 0.7      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.75    | 12 | yes | 15/3  | 19       | 0.5740    | 0.3921    | 0.7      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6085    | 0.6067    | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 10.0    | 12 | no  | 59/0  | 33       | 0.7731    | 0.3104    | 0.6      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.75    | 8  | no  | 29/1  | 20       | 566.3111  | 622.3367  | 0.7      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 4/5   | 4        | 0.7577    | 0.4809    | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 10.0    | 8  | yes | 48/0  | 31       | 0.6248    | 0.3608    | 0.7      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.75    | 8  | yes | 44/0  | 30       | 1985.5275 | 2170.3680 | 0.7      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 4/5   | 4        | 0.7367    | 0.4959    | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 4  | no  | 39/0  | 25       | 0.4346    | 0.4281    | 1.0      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.75    | 12 | no  | 44/2  | 28       | 1623.1074 | 1877.6884 | 0.5      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 4/5   | 4        | 0.7504    | 0.4838    | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 1.0     | 12 | yes | 64/0  | 34       | 0.6099    | 0.3505    | 0.5      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.75    | 12 | yes | 11/3  | 12       | 602.9416  | 770.7203  | 0.5      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 4/5   | 4        | 0.7045    | 0.5081    | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 4  | no  | 54/0  | 30       | 0.4326    | 0.4502    | 1.0      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.75    | 12 | no  | 17/4  | 14       | 654.8884  | 763.2070  | 0.5      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 4/5   | 4        | 0.7226    | 0.4953    | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 8  | yes | 53/0  | 30       | 0.4551    | 0.4299    | 0.7      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.75    | 12 | yes | 29/2  | 22       | 2788.3221 | 2526.6223 | 0.5      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6085    | 0.6067    | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 4  | no  | 0/0   | 7        | 0.0041    | 0.0120    | 1.4      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 1.8      | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 0.1      | 🟢      |
| PolyJet     | Clean | STRidge   | 1.0     | 8  | yes | 0/0   | 7        | 0.0027    | 0.0116    | 1.1      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 7        | 0.0026    | 0.0116    | 2.8      | 🟢      |
| PolyJet     | Clean | RFE       |         | 8  | yes | 0/0   | 7        | 0.0026    | 0.0116    | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 8  | no  | 54/0  | 61       | 0.3781    | 0.2821    | 2.4      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | no  | 9/0   | 16       | 0.3141    | 0.2086    | 2.1      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6051    | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | yes | 61/0  | 68       | 0.5088    | 0.2563    | 1.8      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 65/0  | 72       | 0.5043    | 0.2606    | 1.9      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6085    | 0.5996    | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 8  | no  | 41/0  | 48       | 0.3551    | 0.3177    | 2.2      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.5     | 12 | no  | 6/0   | 13       | 0.4636    | 0.3192    | 1.9      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6154    | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | yes | 64/0  | 71       | 0.4050    | 0.2989    | 1.8      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | yes | 8/2   | 13       | 0.5071    | 0.3693    | 1.8      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6085    | 0.6067    | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.2      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 1.0     | 8  | yes | 0/0   | 3        | 0.0023    | 0.0116    | 0.2      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 3        | 0.0023    | 0.0116    | 0.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 8  | yes | 0/0   | 3        | 0.0023    | 0.0116    | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.2711    | 0.2656    | 0.2      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.2712    | 0.2656    | 0.2      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6051    | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | yes | 30/0  | 15       | 0.4513    | 0.3454    | 0.2      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 10/0  | 8        | 0.5112    | 0.3582    | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6085    | 0.5996    | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 13/0  | 9        | 0.4499    | 0.3141    | 0.2      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 8  | no  | 7/0   | 7        | 0.4445    | 0.3947    | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.6489    | 0.4886    | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.1     | 8  | yes | 32/0  | 18       | 0.4034    | 0.4129    | 0.3      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | yes | 7/0   | 7        | 0.4906    | 0.3949    | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6085    | 0.6067    | 0.0      | 🔴      |

---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ks_cell`.

| Source                                   | Wall time | CPU                                      | Threads |
|:-----------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260802_070102_ks_cell_full.csv` | 77.8 min  | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | 1       |
