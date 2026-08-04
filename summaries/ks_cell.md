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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best             | Mg E_data | Std best      | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-----------------------|:----------|:--------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | 🟩        | 0/0 ✓             | 0.0040     | mod_multigraph/RFE     | 0.0123    | PolyJet/RFE   | 0/0 ✓     | 0.0041      | 0.0122     | 0.0129       |
| Clean | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.0040     | mod_multigraph/RFE     | 0.0123    | PolyJet/RFE   | 0/0 ✓     | 0.0041      | 0.0122     | 0.0129       |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0027     | mod_multigraph/STRidge | 0.0088    | PolyJet/MSTLS | 0/0 ✓     | 0.0030      | 0.0088     | 0.0092       |
| Clean | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0027     | mod_multigraph/STRidge | 0.0088    | PolyJet/MSTLS | 0/0 ✓     | 0.0030      | 0.0088     | 0.0092       |
| Clean | weak (IBP)    | spectral   | 🟢        | 0/0 ✓             | 0.1234     | multigraph/STRidge     | 0.0760    | —             | 1/0 ✗     | 0.1389      | 0.0081     | 0.1572       |
| Clean | weak (IBP)    | polynomial | 🟢        | 0/0 ✓             | 0.1234     | multigraph/STRidge     | 0.0760    | —             | 1/0 ✗     | 0.1389      | 0.0081     | 0.1572       |
| 2%    | pointwise     | spectral   | —        | 37/0 ✗            | 0.4675     | —                      | 0.5297    | —             | 60/0 ✗    | 0.3941      | 0.4714     | 9249.6733    |
| 2%    | pointwise     | polynomial | —        | 37/0 ✗            | 0.3850     | —                      | 0.4488    | —             | 98/0 ✗    | 0.3351      | 0.4002     | 204.6417     |
| 2%    | weak (no IBP) | spectral   | —        | 9/0 ✗             | 0.4585     | —                      | 0.3719    | —             | 35/0 ✗    | 0.4372      | 0.2031     | 2.9497       |
| 2%    | weak (no IBP) | polynomial | 🟢        | 0/0 ✓             | 0.2711     | multigraph/STRidge     | 0.2656    | —             | 9/0 ✗     | 0.3141      | 0.2086     | 0.4322       |
| 2%    | weak (IBP)    | spectral   | —        | 8/0 ✗             | 1.8988     | —                      | 0.3222    | —             | 35/0 ✗    | 0.4718      | 0.2033     | 3.0651       |
| 2%    | weak (IBP)    | polynomial | 🟢        | 0/0 ✓             | 0.3032     | multigraph/MSTLS       | 0.2617    | —             | 4/0 ✗     | 0.3327      | 0.2149     | 0.5150       |
| 5%    | pointwise     | spectral   | —        | 37/0 ✗            | 0.5283     | —                      | 0.6742    | —             | 38/0 ✗    | 0.4714      | 0.6538     | 20670.3966   |
| 5%    | pointwise     | polynomial | —        | 34/0 ✗            | 0.4307     | —                      | 0.6387    | —             | 91/0 ✗    | 0.3660      | 0.6009     | 457.3188     |
| 5%    | weak (no IBP) | spectral   | —        | 15/0 ✗            | 0.6630     | —                      | 0.4714    | —             | 38/0 ✗    | 0.4873      | 0.3841     | 18.2763      |
| 5%    | weak (no IBP) | polynomial | —        | 7/0 ✗             | 0.4445     | —                      | 0.3947    | —             | 6/0 ✗     | 0.4636      | 0.3192     | 1.6768       |
| 5%    | weak (IBP)    | spectral   | —        | 4/0 ✗             | 0.6365     | —                      | 0.4616    | —             | 38/0 ✗    | 0.4963      | 0.3841     | 18.3906      |
| 5%    | weak (IBP)    | polynomial | —        | 7/0 ✗             | 0.4525     | —                      | 0.3949    | —             | 2/0 ✗     | 0.4684      | 0.3447     | 1.7805       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                      | fp/fn  | E_coeff   | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|----------:|---------:|---------:|
| Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0027    | 0.3      | —        |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0027    | 0.2      | +0.0%    |
| Modified Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0027    | 0.3      | +0.2%    |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0027    | 0.1      | +0.2%    |
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0027    | 0.3      | +0.2%    |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0027    | 0.1      | +0.2%    |
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0030    | 1.9      | +13.9%   |
| PolyJet + RFE                  | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0030    | 0.1      | +13.9%   |
| PolyJet + STRidge              | spectral, weak (no IBP), m=8    | 0/0 ✓  | 0.0037    | 0.9      | +38.3%   |
| WSINDy + STRidge               | spectral, weak (IBP), m=8       | 29/0 ✗ | 1.1328    | 0.6      | ✗        |
| WSINDy + MSTLS                 | polynomial, weak (IBP), m=8     | 50/0 ✗ | 2055.6319 | 0.8      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=4    | 28/2 ✗ | 0.4584    | 1.7      | ✗        |
| PDE-FIND + STRidge             | spectral, weak (IBP), m=8       | 28/2 ✗ | 1.4997    | 1.0      | ✗        |
| PDE-FIND + RFE                 | polynomial, weak (no IBP), m=12 | 23/4 ✗ | 0.9619    | 0.1      | ✗        |
| WSINDy + RFE                   | polynomial, weak (IBP), m=12    | 4/5 ✗  | 0.7576    | 0.1      | ✗        |

### 2% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff   | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|----------:|---------:|---------:|
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.2711    | 0.2      | —        |
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.2712    | 0.2      | +0.1%    |
| Modified Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.2712    | 0.2      | +0.1%    |
| PolyJet + MSTLS                | polynomial, weak (IBP), m=12    | 4/0 ✗  | 0.3327    | 2.3      | ✗        |
| PolyJet + STRidge              | polynomial, weak (IBP), m=12    | 25/0 ✗ | 0.3054    | 1.1      | ✗        |
| Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=8  | 28/0 ✗ | 0.4081    | 0.3      | ✗        |
| WSINDy + STRidge               | polynomial, weak (no IBP), m=4  | 39/0 ✗ | 0.4346    | 1.2      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=12   | 23/2 ✗ | 0.6171    | 0.8      | ✗        |
| PDE-FIND + STRidge             | polynomial, weak (IBP), m=8     | 33/2 ✗ | 0.4869    | 1.8      | ✗        |
| WSINDy + MSTLS                 | polynomial, weak (IBP), m=12    | 44/2 ✗ | 1618.5991 | 0.5      | ✗        |
| Multigraphs + RFE              | spectral, weak (no IBP), m=8    | 0/5 ✗  | 0.7796    | 0.1      | ✗        |
| Modified Multigraphs + RFE     | spectral, weak (no IBP), m=8    | 0/5 ✗  | 0.7796    | 0.1      | ✗        |
| WSINDy + RFE                   | polynomial, weak (IBP), m=12    | 4/5 ✗  | 0.7504    | 0.1      | ✗        |
| PDE-FIND + RFE                 | polynomial, pointwise           | 1/7 ✗  | 0.6070    | 0.1      | ✗        |
| PolyJet + RFE                  | polynomial, pointwise           | 1/7 ✗  | 0.6070    | 0.1      | ✗        |

### 5% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff  | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|---------:|---------:|---------:|
| PolyJet + MSTLS                | polynomial, weak (IBP), m=12    | 2/0 ✗  | 0.4684   | 2.1      | ✗        |
| Modified Multigraphs + MSTLS   | spectral, weak (IBP), m=12      | 4/0 ✗  | 0.6365   | 0.2      | ✗        |
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=8  | 7/0 ✗  | 0.4445   | 0.3      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 13/0 ✗ | 0.4499   | 0.3      | ✗        |
| Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 13/0 ✗ | 0.4953   | 0.2      | ✗        |
| PolyJet + STRidge              | spectral, pointwise             | 38/0 ✗ | 0.4714   | 2.3      | ✗        |
| WSINDy + STRidge               | polynomial, weak (no IBP), m=4  | 54/0 ✗ | 0.4326   | 1.2      | ✗        |
| PDE-FIND + MSTLS               | polynomial, weak (no IBP), m=12 | 22/2 ✗ | 0.6175   | 0.9      | ✗        |
| PDE-FIND + STRidge             | spectral, weak (IBP), m=12      | 31/2 ✗ | 0.5766   | 0.8      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=12   | 17/4 ✗ | 654.7380 | 0.6      | ✗        |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=8  | 0/5 ✗  | 0.6489   | 0.1      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=8  | 0/5 ✗  | 0.6489   | 0.1      | ✗        |
| WSINDy + RFE                   | polynomial, weak (IBP), m=12    | 4/5 ✗  | 0.7225   | 0.1      | ✗        |
| PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/6 ✗  | 0.9571   | 0.1      | ✗        |
| PDE-FIND + RFE                 | polynomial, pointwise           | 1/7 ✗  | 0.6070   | 0.1      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 45/2  | 50       | 0.4166  | 0.4099 | 3.0      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 1/7   | 1        | 0.6070  | 0.6047 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.01    | 44/2  | 49       | 0.4349  | 0.4923 | 2.0      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.1     | 42/2  | 47       | 0.4732  | 0.6544 | 1.8      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 5/7   | 5        | 0.5678  | 0.7065 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.01    | 62/0  | 33       | 0.4222  | 0.4158 | 2.2      | 🔴      |
| WSINDy           | Clean | RFE       |         | 1/7   | 1        | 0.6070  | 0.6047 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 79/0  | 39       | 0.4365  | 0.4941 | 1.8      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 85/0  | 42       | 0.4760  | 0.6557 | 2.4      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0041  | 0.0122 | 1.1      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0041  | 0.0122 | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.01    | 60/0  | 67       | 0.3941  | 0.4714 | 4.0      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 38/0  | 45       | 0.4714  | 0.6538 | 2.3      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0040  | 0.0123 | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.001   | 37/0  | 20       | 0.4675  | 0.5297 | 0.5      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 37/0  | 20       | 0.5283  | 0.6742 | 0.7      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0040  | 0.0123 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.001   | 37/0  | 20       | 0.4675  | 0.5297 | 0.5      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 37/0  | 20       | 0.5283  | 0.6742 | 0.5      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 45/2  | 50       | 0.4166  | 0.4099 | 1.8      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 1/7   | 1        | 0.6070  | 0.6047 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.01    | 38/2  | 43       | 0.4278  | 0.4712 | 3.4      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.001   | 53/2  | 58       | 0.4306  | 0.6237 | 3.3      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.01    | 62/0  | 33       | 0.4222  | 0.4158 | 1.4      | 🔴      |
| WSINDy           | Clean | RFE       |         | 1/7   | 1        | 0.6070  | 0.6047 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.001   | 98/1  | 45       | 0.4290  | 0.4722 | 1.8      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.001   | 72/1  | 38       | 0.4318  | 0.6248 | 2.0      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0041  | 0.0122 | 1.3      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0041  | 0.0122 | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.001   | 98/0  | 105      | 0.3351  | 0.4002 | 6.8      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.001   | 91/0  | 98       | 0.3660  | 0.6009 | 6.2      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0040  | 0.0123 | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.001   | 37/0  | 20       | 0.3850  | 0.4488 | 0.5      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 37/0  | 20       | 0.4307  | 0.6387 | 0.6      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0040  | 0.0123 | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.01    | 37/0  | 17       | 0.3850  | 0.4488 | 0.5      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.6309 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 34/0  | 19       | 0.4307  | 0.6387 | 0.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.6070  | 0.7219 | 0.1      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff   | E_data    | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:----------|:----------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 10.0    | 8  | no  | 31/2  | 36       | 0.6789    | 0.3706    | 1.2      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.75    | 4  | no  | 28/2  | 33       | 0.4584    | 0.4195    | 1.7      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | no  | 23/4  | 26       | 0.9619    | 0.3547    | 0.1      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 10.0    | 8  | yes | 28/2  | 33       | 1.4997    | 0.3710    | 1.0      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.75    | 8  | yes | 46/2  | 51       | 0.9119    | 0.3521    | 1.3      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | yes | 10/6  | 11       | 0.9062    | 0.4740    | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.01    | 4  | no  | 54/2  | 59       | 0.4347    | 0.4304    | 2.6      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.75    | 12 | no  | 23/2  | 28       | 0.6171    | 0.3543    | 0.8      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6051    | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 8  | yes | 39/2  | 44       | 0.4690    | 0.4098    | 1.2      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.75    | 12 | yes | 23/2  | 28       | 0.6196    | 0.3543    | 0.7      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6084    | 0.5924    | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.1     | 4  | no  | 42/2  | 47       | 0.4355    | 0.4546    | 1.9      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.75    | 12 | no  | 28/2  | 33       | 0.5562    | 0.3650    | 0.8      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6154    | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 31/2  | 36       | 0.5766    | 0.3657    | 0.8      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.75    | 12 | yes | 28/2  | 33       | 0.5574    | 0.3650    | 0.9      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6084    | 0.5955    | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 10.0    | 12 | no  | 59/0  | 33       | 0.7731    | 0.3104    | 0.6      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.75    | 8  | no  | 29/1  | 20       | 566.3077  | 622.3344  | 0.8      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 4/5   | 4        | 0.7577    | 0.4809    | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 10.0    | 8  | yes | 29/0  | 21       | 1.1328    | 0.3691    | 0.6      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.75    | 8  | yes | 50/0  | 31       | 2055.6337 | 2251.0317 | 0.8      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 4/5   | 4        | 0.7576    | 0.4809    | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.001   | 4  | no  | 98/0  | 46       | 0.4344    | 0.4287    | 1.9      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.75    | 12 | no  | 44/2  | 28       | 1625.2772 | 1879.8458 | 0.5      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 4/5   | 4        | 0.7504    | 0.4838    | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.1     | 8  | yes | 50/0  | 32       | 0.5014    | 0.3936    | 0.8      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.75    | 12 | yes | 44/2  | 28       | 1619.7628 | 1874.2926 | 0.5      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 4/5   | 4        | 0.7504    | 0.4838    | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.001   | 4  | no  | 98/0  | 46       | 0.4337    | 0.4533    | 2.3      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.75    | 12 | no  | 17/4  | 14       | 654.7380  | 763.0814  | 0.6      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | no  | 4/5   | 4        | 0.7226    | 0.4954    | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.001   | 8  | yes | 98/0  | 46       | 0.4897    | 0.4053    | 1.0      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.75    | 12 | yes | 17/4  | 14       | 654.7436  | 763.0839  | 0.6      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | yes | 4/5   | 4        | 0.7226    | 0.4954    | 0.0      | 🔴      |
| PolyJet          | Clean | STRidge   | 1.0     | 8  | no  | 0/0   | 7        | 0.0037    | 0.0108    | 0.9      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 1.6      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 1.0     | 8  | yes | 1/0   | 8        | 3.5954    | 0.0101    | 0.9      | 🔴      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 1/0   | 8        | 0.1389    | 0.0081    | 1.7      | 🔴      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 1/0   | 8        | 0.1389    | 0.0081    | 0.1      | 🔴      |
| PolyJet          | 2%    | STRidge   | 1.0     | 8  | no  | 38/0  | 45       | 0.3767    | 0.3079    | 1.7      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 12 | no  | 35/0  | 42       | 0.4372    | 0.2031    | 1.8      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6051    | 0.1      | 🔴      |
| PolyJet          | 2%    | STRidge   | 1.0     | 8  | yes | 45/0  | 52       | 0.3923    | 0.3075    | 1.7      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 12 | yes | 35/0  | 42       | 0.4718    | 0.2033    | 1.8      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6084    | 0.5924    | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 8  | no  | 82/0  | 89       | 0.4893    | 0.3783    | 3.0      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 8  | no  | 38/0  | 45       | 0.4873    | 0.3841    | 2.3      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 8  | no  | 0/6   | 1        | 0.9571    | 0.6513    | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 12 | yes | 58/0  | 65       | 0.7440    | 0.3088    | 1.7      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 8  | yes | 38/0  | 45       | 0.4963    | 0.3841    | 2.3      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 0/6   | 1        | 0.9571    | 0.6513    | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.2      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.2      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.1234    | 0.0760    | 0.2      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.1234    | 0.0760    | 0.2      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 1/0   | 4        | 0.1388    | 0.0082    | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | no  | 11/0  | 7        | 0.3770    | 0.2649    | 0.2      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 8  | no  | 9/0   | 8        | 0.4585    | 0.3719    | 0.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.7796    | 0.5441    | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | yes | 10/0  | 7        | 0.4550    | 0.2708    | 0.2      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | yes | 13/0  | 8        | 0.5681    | 0.2535    | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 8  | yes | 0/5   | 1        | 0.7796    | 0.5441    | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | no  | 28/0  | 16       | 0.9366    | 0.3920    | 0.3      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 12 | no  | 18/0  | 12       | 0.9587    | 0.4005    | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6154    | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | yes | 28/0  | 16       | 0.9427    | 0.3920    | 0.3      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 12 | yes | 18/0  | 12       | 0.9650    | 0.4005    | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 2/7   | 1        | 0.9483    | 0.5977    | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.1234    | 0.0760    | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.1234    | 0.0760    | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 1/0   | 4        | 0.1388    | 0.0082    | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 0.001   | 4  | no  | 37/0  | 20       | 0.4157    | 0.4085    | 0.5      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 12/0  | 9        | 0.4602    | 0.3719    | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.7796    | 0.5441    | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 12 | yes | 35/0  | 17       | 0.8495    | 0.2394    | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.5     | 12 | yes | 8/0   | 5        | 1.8988    | 0.3222    | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 8  | yes | 0/5   | 1        | 0.7796    | 0.5441    | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 8  | no  | 35/0  | 16       | 0.6373    | 0.4633    | 0.4      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 8  | no  | 15/0  | 11       | 0.6630    | 0.4714    | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6154    | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | yes | 36/0  | 17       | 1.1666    | 0.4086    | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.5     | 12 | yes | 4/0   | 5        | 0.6365    | 0.4616    | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6084    | 0.5955    | 0.1      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff   | E_data    | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:----------|:----------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 10.0    | 8  | no  | 31/2  | 36       | 0.6789    | 0.3706    | 1.1      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.75    | 4  | no  | 28/2  | 33       | 0.4584    | 0.4195    | 1.5      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | no  | 23/4  | 26       | 0.9619    | 0.3547    | 0.1      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 10.0    | 8  | yes | 28/2  | 33       | 1.4999    | 0.3710    | 0.9      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.75    | 8  | yes | 46/2  | 51       | 0.9121    | 0.3521    | 1.2      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | yes | 10/6  | 11       | 0.9061    | 0.4740    | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | no  | 43/2  | 48       | 0.6880    | 0.3395    | 1.4      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.75    | 12 | no  | 25/2  | 30       | 0.6531    | 0.3425    | 1.1      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6051    | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 8  | yes | 33/2  | 38       | 0.4869    | 0.4062    | 1.8      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.75    | 12 | yes | 25/2  | 30       | 0.6553    | 0.3425    | 1.1      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6084    | 0.5924    | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.1     | 8  | no  | 34/2  | 39       | 0.4654    | 0.4194    | 1.4      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.75    | 12 | no  | 22/2  | 27       | 0.6175    | 0.3632    | 0.9      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6154    | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.1     | 8  | yes | 34/2  | 39       | 0.4680    | 0.4194    | 1.3      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.75    | 12 | yes | 25/2  | 30       | 0.6212    | 0.3595    | 0.9      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6084    | 0.5955    | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 1.0     | 8  | no  | 51/0  | 32       | 0.6387    | 0.3554    | 0.9      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.75    | 8  | no  | 29/1  | 20       | 566.3111  | 622.3367  | 0.8      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 4/5   | 4        | 0.7577    | 0.4809    | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 12 | yes | 98/0  | 46       | 0.7742    | 0.3076    | 0.8      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.75    | 8  | yes | 50/0  | 31       | 2055.6319 | 2251.0293 | 0.8      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 4/5   | 4        | 0.7576    | 0.4809    | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.1     | 4  | no  | 39/0  | 25       | 0.4346    | 0.4281    | 1.2      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.75    | 12 | no  | 44/2  | 28       | 1623.1074 | 1877.6884 | 0.5      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 4/5   | 4        | 0.7504    | 0.4838    | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.001   | 8  | yes | 98/0  | 46       | 0.5232    | 0.3887    | 0.9      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.75    | 12 | yes | 44/2  | 28       | 1618.5991 | 1873.1182 | 0.5      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 4/5   | 4        | 0.7504    | 0.4838    | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 1.0     | 4  | no  | 54/0  | 30       | 0.4326    | 0.4502    | 1.2      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.75    | 12 | no  | 17/4  | 14       | 654.8884  | 763.2070  | 0.6      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | no  | 4/5   | 4        | 0.7226    | 0.4953    | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.001   | 8  | yes | 98/0  | 46       | 0.4895    | 0.4053    | 1.2      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.75    | 12 | yes | 17/4  | 14       | 654.7410  | 763.0583  | 0.6      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | yes | 4/5   | 4        | 0.7225    | 0.4954    | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 1.0     | 4  | no  | 0/0   | 7        | 0.0041    | 0.0120    | 1.3      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 1.9      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0030    | 0.0088    | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 1.0     | 8  | yes | 1/0   | 8        | 3.5954    | 0.0101    | 1.3      | 🔴      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 1/0   | 8        | 0.1389    | 0.0081    | 2.0      | 🔴      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 1/0   | 8        | 0.1389    | 0.0081    | 0.1      | 🔴      |
| PolyJet          | 2%    | STRidge   | 1.0     | 8  | no  | 54/0  | 61       | 0.3781    | 0.2821    | 2.7      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 12 | no  | 9/0   | 16       | 0.3141    | 0.2086    | 2.3      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6051    | 0.1      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 12 | yes | 25/0  | 32       | 0.3054    | 0.1792    | 1.1      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 12 | yes | 4/0   | 11       | 0.3327    | 0.2149    | 2.3      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6084    | 0.5924    | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | no  | 41/0  | 48       | 0.3551    | 0.3177    | 2.2      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | no  | 6/0   | 13       | 0.4636    | 0.3192    | 2.1      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6154    | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 10.0    | 12 | yes | 41/0  | 48       | 0.3900    | 0.2456    | 1.3      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | yes | 2/0   | 9        | 0.4684    | 0.3447    | 2.1      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6084    | 0.5955    | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.2      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.1234    | 0.0760    | 0.3      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.1234    | 0.0760    | 0.2      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 1/0   | 4        | 0.1388    | 0.0082    | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.2711    | 0.2656    | 0.2      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.2712    | 0.2656    | 0.2      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6051    | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.3034    | 0.2617    | 0.2      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.3032    | 0.2617    | 0.2      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6084    | 0.5924    | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | no  | 13/0  | 9        | 0.4499    | 0.3141    | 0.3      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 8  | no  | 7/0   | 7        | 0.4445    | 0.3947    | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.6489    | 0.4886    | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | yes | 20/0  | 13       | 0.4493    | 0.3068    | 0.3      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 8  | yes | 7/0   | 7        | 0.4525    | 0.3949    | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 0/5   | 1        | 0.6489    | 0.4886    | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0027    | 0.0088    | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 4        | 0.1387    | 0.0082    | 0.3      | 🔴      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.1234    | 0.0760    | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 1/0   | 4        | 0.1388    | 0.0082    | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 8  | no  | 28/0  | 12       | 0.4081    | 0.3469    | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.2712    | 0.2656    | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.6073    | 0.6051    | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 0.001   | 12 | yes | 37/0  | 20       | 0.2640    | 0.2022    | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.5     | 12 | yes | 8/0   | 5        | 0.3001    | 0.2338    | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.6084    | 0.5924    | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 13/0  | 8        | 0.4953    | 0.3252    | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 8  | no  | 7/0   | 7        | 0.4506    | 0.3949    | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.6489    | 0.4886    | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | yes | 35/0  | 17       | 0.5477    | 0.2830    | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 12 | yes | 9/0   | 7        | 0.5593    | 0.3147    | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | yes | 0/5   | 1        | 0.6489    | 0.4886    | 0.1      | 🔴      |


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
| Clean | weak (IBP)    | spectral   | 8  | yes | 3.7624       | 133848.0  |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.1572       | 126840.5  |
| Clean | weak (IBP)    | polynomial | 8  | yes | 3.7624       | 133850.2  |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.1572       | 126842.8  |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 88.9904      | 336798.9  |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 3.9765       | 471736.7  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 2.9497       | 693750.9  |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 33.2300      | 55772.9   |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 2.5916       | 49341.4   |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.4322       | 68936.2   |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 6.7089       | 472036.4  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 3.0651       | 693845.0  |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 4.8434       | 49347.7   |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.5150       | 68934.5   |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 220.7844     | 2029412.3 |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 18.7254      | 3418856.0 |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 18.2763      | 4850904.7 |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 82.2495      | 142599.5  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 6.5727       | 58296.0   |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 1.6768       | 79176.8   |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 21.4031      | 3419390.6 |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 18.3906      | 4851117.8 |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 8.5625       | 58182.9   |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 1.7805       | 79178.9   |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ks_cell`.

| Source                                       | Wall time | CPU | Threads |
|:---------------------------------------------|----------:|:----|:-------:|
| `sweep_no_ac_ch_ks_cell_20260803_170832.csv` | ? min     | ?   | ?       |
