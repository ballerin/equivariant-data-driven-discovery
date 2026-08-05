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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best             | Mg E_data | Std best        | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-----------------------|:----------|:----------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | ⚪        | 0/0 ✓             | 0.0047     | multigraph/RFE         | 0.0176    | PolyJet/RFE     | 0/0 ✓     | 0.0047      | 0.0175     | 0.0182       |
| Clean | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.0119     | multigraph/RFE         | 0.0325    | PolyJet/RFE     | 0/0 ✓     | 0.0130      | 0.0320     | 0.0345       |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0021     | mod_multigraph/STRidge | 0.0112    | PolyJet/STRidge | 0/0 ✓     | 0.0023      | 0.0112     | 0.0114       |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0022     | mod_multigraph/RFE     | 0.0112    | PolyJet/MSTLS   | 0/0 ✓     | 0.0023      | 0.0112     | 0.0114       |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0278     | multigraph/STRidge     | 0.0158    | PolyJet/STRidge | 0/0 ✓     | 0.0280      | 0.0157     | 0.0326       |
| Clean | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0278     | mod_multigraph/STRidge | 0.0158    | PolyJet/STRidge | 0/0 ✓     | 0.0279      | 0.0157     | 0.0326       |
| 2%    | pointwise     | spectral   | —        | 17/0 ✗            | 0.4267     | —                      | 0.6377    | —               | 83/0 ✗    | 0.4249      | 0.6367     | 1381.7721    |
| 2%    | pointwise     | polynomial | —        | 28/0 ✗            | 0.3683     | —                      | 0.6210    | —               | 79/0 ✗    | 0.4282      | 0.6371     | 29.8498      |
| 2%    | weak (no IBP) | spectral   | —        | 6/0 ✗             | 0.5485     | —                      | 0.3476    | —               | 62/0 ✗    | 0.4969      | 0.3239     | 11.3180      |
| 2%    | weak (no IBP) | polynomial | —        | 3/0 ✗             | 0.1979     | —                      | 0.1320    | —               | 9/0 ✗     | 0.0962      | 0.1099     | 0.8842       |
| 2%    | weak (IBP)    | spectral   | —        | 9/0 ✗             | 0.8931     | —                      | 0.3082    | —               | 11/0 ✗    | 0.5803      | 0.3359     | 11.4293      |
| 2%    | weak (IBP)    | polynomial | —        | 1/0 ✗             | 0.6262     | —                      | 0.1858    | —               | 9/0 ✗     | 0.1000      | 0.1098     | 0.9103       |
| 5%    | pointwise     | spectral   | —        | 26/0 ✗            | 0.4316     | —                      | 0.8725    | —               | 89/0 ✗    | 0.4310      | 0.8722     | 2190.0271    |
| 5%    | pointwise     | polynomial | —        | 33/0 ✗            | 0.4217     | —                      | 0.8717    | —               | 93/0 ✗    | 0.4190      | 0.8702     | 47.4252      |
| 5%    | weak (no IBP) | spectral   | —        | 30/0 ✗            | 0.4313     | —                      | 0.5369    | —               | 59/0 ✗    | 0.4459      | 0.4115     | 69.4094      |
| 5%    | weak (no IBP) | polynomial | —        | 8/0 ✗             | 0.3919     | —                      | 0.2940    | —               | 13/0 ✗    | 0.2558      | 0.2611     | 5.4493       |
| 5%    | weak (IBP)    | spectral   | —        | 31/0 ✗            | 0.4804     | —                      | 0.3565    | —               | 70/0 ✗    | 0.4452      | 0.4134     | 70.0482      |
| 5%    | weak (IBP)    | polynomial | —        | 8/0 ✗             | 0.3937     | —                      | 0.2934    | —               | 13/0 ✗    | 0.2518      | 0.2604     | 5.4753       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                      | fp/fn  | E_coeff     | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|------------:|---------:|---------:|
| Modified Multigraphs + STRidge | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0021      | 0.3      | —        |
| Multigraphs + STRidge          | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0021      | 0.3      | +0.0%    |
| Modified Multigraphs + MSTLS   | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0022      | 0.4      | +3.6%    |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0022      | 0.3      | +3.6%    |
| Multigraphs + RFE              | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0022      | 0.1      | +3.6%    |
| Modified Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0022      | 0.1      | +3.6%    |
| PolyJet + STRidge              | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0023      | 1.3      | +6.9%    |
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0023      | 3.9      | +8.3%    |
| PolyJet + RFE                  | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0023      | 0.1      | +8.3%    |
| WSINDy + STRidge               | polynomial, weak (IBP), m=8     | 63/0 ✗ | 0.4328      | 1.6      | ✗        |
| PDE-FIND + STRidge             | spectral, weak (IBP), m=8       | 20/2 ✗ | 0.4489      | 1.5      | ✗        |
| PDE-FIND + MSTLS               | polynomial, weak (IBP), m=8     | 1/5 ✗  | 0.4535      | 1.9      | ✗        |
| PDE-FIND + RFE                 | spectral, pointwise             | 1/7 ✗  | 0.4409      | 0.1      | ✗        |
| WSINDy + RFE                   | spectral, pointwise             | 1/7 ✗  | 0.4409      | 0.1      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=4    | 4/7 ✗  | 444689.2112 | 1.7      | ✗        |

### 2% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff     | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|------------:|---------:|---------:|
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=12    | 1/0 ✗  | 0.6262      | 0.6      | ✗        |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=8    | 6/0 ✗  | 0.5485      | 0.5      | ✗        |
| Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 8/0 ✗  | 0.1951      | 0.7      | ✗        |
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=12 | 9/0 ✗  | 0.0962      | 4.2      | ✗        |
| PolyJet + STRidge              | polynomial, weak (no IBP), m=8  | 10/0 ✗ | 0.2711      | 2.3      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (IBP), m=12    | 12/0 ✗ | 0.1460      | 0.4      | ✗        |
| WSINDy + STRidge               | polynomial, pointwise           | 79/0 ✗ | 0.4282      | 3.2      | ✗        |
| PDE-FIND + STRidge             | spectral, weak (no IBP), m=4    | 23/2 ✗ | 0.4283      | 4.9      | ✗        |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.5636      | 0.1      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.5636      | 0.2      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=12   | 7/5 ✗  | 0.4447      | 2.7      | ✗        |
| PolyJet + RFE                  | polynomial, weak (no IBP), m=12 | 0/6 ✗  | 0.7087      | 0.1      | ✗        |
| PDE-FIND + RFE                 | polynomial, pointwise           | 1/7 ✗  | 0.4409      | 0.1      | ✗        |
| WSINDy + RFE                   | spectral, pointwise             | 1/7 ✗  | 0.4409      | 0.2      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=4    | 4/7 ✗  | 438342.0538 | 1.8      | ✗        |

### 5% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff     | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|------------:|---------:|---------:|
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=12 | 8/0 ✗  | 0.3919      | 0.3      | ✗        |
| Modified Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 9/0 ✗  | 0.4122      | 0.4      | ✗        |
| PolyJet + MSTLS                | polynomial, weak (IBP), m=12    | 13/0 ✗ | 0.2518      | 4.2      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 20/0 ✗ | 0.4188      | 0.4      | ✗        |
| Modified Multigraphs + STRidge | spectral, pointwise             | 26/0 ✗ | 0.4316      | 1.6      | ✗        |
| PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 37/0 ✗ | 0.2726      | 2.9      | ✗        |
| WSINDy + STRidge               | polynomial, weak (no IBP), m=4  | 85/0 ✗ | 0.4284      | 3.9      | ✗        |
| PDE-FIND + STRidge             | polynomial, weak (IBP), m=12    | 27/2 ✗ | 0.4370      | 1.2      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/5 ✗  | 0.8584      | 0.1      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=4    | 1/7 ✗  | 0.4409      | 2.7      | ✗        |
| Multigraphs + RFE              | spectral, pointwise             | 1/7 ✗  | 0.5796      | 0.2      | ✗        |
| WSINDy + RFE                   | polynomial, pointwise           | 1/7 ✗  | 0.5796      | 0.1      | ✗        |
| PolyJet + RFE                  | spectral, pointwise             | 1/7 ✗  | 0.5796      | 0.2      | ✗        |
| PDE-FIND + RFE                 | spectral, pointwise             | 1/7 ✗  | 0.5796      | 0.1      | ✗        |
| WSINDy + MSTLS                 | polynomial, weak (no IBP), m=4  | 4/7 ✗  | 408919.0648 | 2.0      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.001   | 47/2  | 52       | 0.4274  | 0.4150 | 4.9      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 1/7   | 1        | 0.4409  | 0.4294 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 44/2  | 49       | 0.4292  | 0.6377 | 8.6      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.3      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.01    | 47/2  | 52       | 0.4308  | 0.8720 | 4.8      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 98/0  | 46       | 0.4276  | 0.4152 | 3.5      | 🔴      |
| WSINDy           | Clean | RFE       |         | 1/7   | 1        | 0.4409  | 0.4294 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 94/0  | 45       | 0.4291  | 0.6378 | 3.4      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.2      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 89/0  | 44       | 0.4310  | 0.8722 | 3.3      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.2      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0047  | 0.0175 | 2.5      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0047  | 0.0175 | 0.3      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.01    | 83/0  | 90       | 0.4249  | 0.6367 | 9.5      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.01    | 91/0  | 98       | 0.4309  | 0.8714 | 11.6     | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0047  | 0.0176 | 0.5      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0047  | 0.0176 | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 18/0  | 12       | 0.4266  | 0.6378 | 1.0      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.01    | 28/0  | 16       | 0.4318  | 0.8725 | 1.8      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0047  | 0.0176 | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0047  | 0.0176 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 17/0  | 11       | 0.4267  | 0.6377 | 1.0      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.01    | 26/0  | 13       | 0.4316  | 0.8725 | 1.6      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.2      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 24/2  | 29       | 0.4279  | 0.4158 | 2.3      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 1/7   | 1        | 0.4409  | 0.4294 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.001   | 53/2  | 58       | 0.4285  | 0.6374 | 6.2      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.01    | 50/2  | 55       | 0.4296  | 0.8716 | 4.7      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 83/0  | 40       | 0.4280  | 0.4157 | 3.1      | 🔴      |
| WSINDy           | Clean | RFE       |         | 1/7   | 1        | 0.4409  | 0.4294 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.001   | 79/0  | 41       | 0.4282  | 0.6371 | 3.2      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 93/0  | 44       | 0.4297  | 0.8718 | 3.8      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0132  | 0.0320 | 2.5      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0130  | 0.0320 | 0.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.001   | 98/0  | 105      | 0.3624  | 0.6182 | 11.9     | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.01    | 93/0  | 100      | 0.4190  | 0.8702 | 11.7     | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.3      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0120  | 0.0325 | 1.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0119  | 0.0325 | 0.2      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.01    | 28/0  | 16       | 0.3683  | 0.6210 | 1.4      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 33/0  | 18       | 0.4217  | 0.8717 | 1.6      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0120  | 0.0325 | 0.8      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0119  | 0.0325 | 0.3      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 33/0  | 14       | 0.3683  | 0.6210 | 2.0      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 35/0  | 19       | 0.4218  | 0.8717 | 1.6      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.1      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff     | E_data      | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:------------|:------------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 4  | no  | 43/2  | 48       | 0.4284      | 0.4002      | 6.5      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.4173      | 2.9      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5615      | 0.1      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 1.0     | 8  | yes | 20/2  | 25       | 0.4489      | 0.3656      | 1.5      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.75    | 8  | yes | 1/5   | 3        | 0.4535      | 0.3811      | 2.1      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5374      | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 4  | no  | 23/2  | 28       | 0.4283      | 0.4295      | 4.9      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.75    | 12 | no  | 7/5   | 9        | 0.4447      | 0.3361      | 2.7      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5778      | 0.3      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | yes | 38/2  | 43       | 0.4435      | 0.3253      | 3.2      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.75    | 12 | yes | 7/5   | 9        | 0.4448      | 0.3361      | 3.0      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5422      | 0.3      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 34/2  | 39       | 0.4432      | 0.3506      | 1.5      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.5449      | 2.7      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.2      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 34/2  | 39       | 0.4430      | 0.3506      | 1.5      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4409      | 0.4298      | 2.4      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5652      | 0.2      | 🔴      |
| WSINDy           | Clean | STRidge   | 1.0     | 8  | no  | 76/0  | 38       | 0.4334      | 0.3644      | 1.6      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 444689.2112 | 448874.7916 | 1.7      | 🔴      |
| WSINDy           | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5615      | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 8  | yes | 98/0  | 46       | 0.4345      | 0.3638      | 1.8      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 662524.5252 | 675775.7402 | 1.3      | 🔴      |
| WSINDy           | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5374      | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 12 | no  | 88/0  | 43       | 0.4506      | 0.3247      | 1.3      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 438342.0538 | 436503.6332 | 1.8      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5778      | 0.2      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 12 | yes | 88/0  | 43       | 0.4507      | 0.3247      | 1.3      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 658615.2958 | 669444.4221 | 1.4      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5422      | 0.2      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 8  | no  | 85/0  | 43       | 0.4302      | 0.4131      | 1.9      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 408919.0648 | 381256.0793 | 1.8      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.2      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 12 | yes | 90/0  | 45       | 0.4459      | 0.3469      | 1.3      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 644131.2469 | 643358.1219 | 1.4      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5652      | 0.2      | 🔴      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 1.3      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 3.9      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0280      | 0.0157      | 1.3      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 0.0280      | 0.0157      | 4.0      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0280      | 0.0157      | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.1     | 8  | no  | 62/0  | 69       | 0.4969      | 0.3239      | 5.4      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 4  | no  | 12/1  | 18       | 0.4233      | 0.4182      | 5.7      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/6   | 1        | 0.9275      | 0.3548      | 0.1      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 12 | yes | 66/0  | 73       | 0.6562      | 0.2564      | 3.1      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 8  | yes | 11/0  | 18       | 0.5803      | 0.3359      | 4.5      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 0.9275      | 0.3548      | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 8  | no  | 59/0  | 66       | 0.4459      | 0.4115      | 5.8      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 12 | no  | 9/3   | 13       | 0.4430      | 0.3575      | 3.8      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | yes | 70/0  | 77       | 0.4452      | 0.4134      | 5.8      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 12 | yes | 9/3   | 13       | 0.4430      | 0.3575      | 3.8      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5652      | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0021      | 0.0112      | 0.3      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.3      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0279      | 0.0158      | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0279      | 0.0158      | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 4  | no  | 30/0  | 17       | 0.4212      | 0.4176      | 2.0      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 8  | no  | 6/0   | 6        | 0.5485      | 0.3476      | 0.5      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9250      | 0.3445      | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 0.001   | 8  | yes | 37/0  | 20       | 0.5555      | 0.3424      | 0.7      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 12 | yes | 14/0  | 9        | 0.7705      | 0.2872      | 0.4      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9250      | 0.3445      | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.1     | 4  | no  | 30/0  | 17       | 0.4313      | 0.5369      | 2.0      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.5449      | 0.6      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 12 | yes | 31/0  | 17       | 0.4804      | 0.3565      | 0.5      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4409      | 0.4298      | 0.5      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5652      | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0021      | 0.0112      | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0279      | 0.0158      | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0279      | 0.0158      | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 8  | no  | 33/0  | 15       | 0.5413      | 0.3428      | 0.6      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 9/0   | 7        | 0.5490      | 0.3476      | 0.4      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9250      | 0.3445      | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 12 | yes | 33/0  | 18       | 0.8892      | 0.3010      | 0.5      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 9/0   | 7        | 0.8931      | 0.3082      | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9250      | 0.3445      | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | no  | 33/0  | 18       | 0.4683      | 0.3542      | 0.5      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.5449      | 0.7      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 8  | yes | 33/0  | 18       | 0.4440      | 0.4177      | 0.7      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4409      | 0.4298      | 0.5      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5652      | 0.2      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff     | E_data      | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:------------|:------------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | no  | 48/2  | 53       | 0.4648      | 0.3177      | 1.8      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.4173      | 2.6      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5615      | 0.1      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 1.0     | 8  | yes | 20/2  | 25       | 0.4490      | 0.3656      | 1.4      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.75    | 8  | yes | 1/5   | 3        | 0.4535      | 0.3811      | 1.9      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5374      | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | no  | 27/2  | 32       | 0.4395      | 0.3266      | 1.4      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.4432      | 3.0      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5778      | 0.2      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | yes | 27/2  | 32       | 0.4395      | 0.3266      | 1.4      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4409      | 0.3929      | 2.3      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5422      | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 27/2  | 32       | 0.4370      | 0.3491      | 1.2      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.5449      | 2.6      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 27/2  | 32       | 0.4370      | 0.3491      | 1.2      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4409      | 0.4298      | 1.8      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5652      | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 4  | no  | 89/0  | 44       | 0.4285      | 0.4008      | 3.2      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 444689.2112 | 448874.7916 | 1.7      | 🔴      |
| WSINDy           | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5615      | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 8  | yes | 63/0  | 34       | 0.4328      | 0.3653      | 1.6      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 662524.5252 | 675775.7402 | 1.3      | 🔴      |
| WSINDy           | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5374      | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 4  | no  | 97/0  | 45       | 0.4285      | 0.4291      | 3.1      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 438342.0538 | 436503.6332 | 1.7      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5778      | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 12 | yes | 88/0  | 43       | 0.4507      | 0.3247      | 1.2      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 658615.2958 | 669444.4221 | 1.4      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5422      | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 4  | no  | 85/0  | 44       | 0.4284      | 0.5353      | 3.9      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 408919.0648 | 381256.0793 | 2.0      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 12 | yes | 88/0  | 44       | 0.4455      | 0.3471      | 1.4      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 644131.2469 | 643358.1219 | 1.5      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5652      | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0024      | 0.0112      | 1.3      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 3.9      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0279      | 0.0157      | 1.3      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 0.0279      | 0.0157      | 4.0      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0279      | 0.0157      | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 10.0    | 8  | no  | 10/0  | 17       | 0.2711      | 0.2197      | 2.3      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 12 | no  | 9/0   | 16       | 0.0962      | 0.1099      | 4.2      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/6   | 1        | 0.7087      | 0.3248      | 0.1      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 12 | yes | 13/0  | 20       | 0.1047      | 0.1095      | 1.7      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 12 | yes | 9/0   | 16       | 0.1000      | 0.1098      | 4.3      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 0.7087      | 0.3248      | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 10.0    | 12 | no  | 37/0  | 44       | 0.2726      | 0.2533      | 2.9      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | no  | 13/0  | 20       | 0.2558      | 0.2611      | 4.2      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.3      | 🔴      |
| PolyJet          | 5%    | STRidge   | 10.0    | 12 | yes | 47/0  | 54       | 0.2802      | 0.2529      | 3.3      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | yes | 13/0  | 20       | 0.2518      | 0.2604      | 4.2      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5652      | 0.3      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0039      | 0.0141      | 0.7      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.6      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.2      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.5      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.6      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.3      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 4  | no  | 31/0  | 18       | 0.3397      | 0.3618      | 1.5      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 12 | no  | 7/0   | 5        | 0.1322      | 0.1168      | 0.4      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.5636      | 0.2228      | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | yes | 12/0  | 8        | 0.1460      | 0.1165      | 0.4      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 12 | yes | 7/0   | 5        | 0.1416      | 0.1166      | 0.4      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.5636      | 0.2228      | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | no  | 20/0  | 13       | 0.4188      | 0.2859      | 0.4      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 12 | no  | 8/0   | 6        | 0.3919      | 0.2940      | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.1     | 8  | yes | 34/0  | 19       | 0.3961      | 0.3708      | 0.6      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 12 | yes | 8/0   | 6        | 0.3937      | 0.2934      | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 3/7   | 1        | 0.9852      | 0.5084      | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0024      | 0.0112      | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.3      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 8  | no  | 11/0  | 9        | 0.3121      | 0.2279      | 0.9      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 12 | no  | 3/0   | 4        | 0.1979      | 0.1320      | 0.6      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.5636      | 0.2228      | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | yes | 8/0   | 7        | 0.1951      | 0.1309      | 0.7      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 1/0   | 4        | 0.6262      | 0.1858      | 0.6      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.5636      | 0.2228      | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 8  | no  | 30/0  | 17       | 0.3937      | 0.3711      | 0.7      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 12 | no  | 9/0   | 6        | 0.4122      | 0.2941      | 0.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 8  | yes | 37/0  | 20       | 0.5152      | 0.3847      | 0.7      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/5   | 1        | 0.8584      | 0.3690      | 0.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.8584      | 0.3690      | 0.1      | 🔴      |


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
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.7600       | 22310889.7  |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0326       | 34140137.1  |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.7599       | 22365392.2  |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0326       | 34175096.6  |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 18.4247      | 21474809.9  |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 11.3180      | 33213146.7  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 11.4029      | 45377815.8  |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 5.3848       | 17403912.2  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.9636       | 21997505.8  |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.8842       | 28567629.8  |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 11.9736      | 33218011.5  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 11.4293      | 45383276.7  |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 1.6306       | 21967269.0  |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.9103       | 28568312.5  |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 75.2191      | 64231689.3  |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 69.4094      | 93763895.2  |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 70.6711      | 125136634.2 |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 13.4584      | 17372921.2  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 5.4893       | 23903276.1  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 5.4493       | 31307599.2  |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 70.0482      | 93755748.7  |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 70.6973      | 125136813.2 |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 6.1362       | 23791682.7  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 5.4753       | 31307505.2  |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ks_wchaos`.

| Source                                         | Wall time | CPU | Threads |
|:-----------------------------------------------|----------:|:----|:-------:|
| `sweep_paper_v2_ks_wchaos_20260804_223951.csv` | ? min     | ?   | ?       |
