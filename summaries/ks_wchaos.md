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
| Clean | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.0117     | multigraph/RFE         | 0.0325    | PolyJet/RFE     | 0/0 ✓     | 0.0128      | 0.0320     | 0.0345       |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0021     | mod_multigraph/STRidge | 0.0112    | PolyJet/STRidge | 0/0 ✓     | 0.0023      | 0.0112     | 0.0114       |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0022     | mod_multigraph/RFE     | 0.0112    | PolyJet/MSTLS   | 0/0 ✓     | 0.0023      | 0.0112     | 0.0114       |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0278     | multigraph/STRidge     | 0.0158    | PolyJet/STRidge | 0/0 ✓     | 0.0280      | 0.0157     | 0.0326       |
| Clean | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0277     | mod_multigraph/STRidge | 0.0158    | PolyJet/STRidge | 0/0 ✓     | 0.0279      | 0.0157     | 0.0326       |
| 2%    | pointwise     | spectral   | —        | 17/0 ✗            | 0.4267     | —                      | 0.6377    | —               | 83/0 ✗    | 0.4249      | 0.6367     | 1381.7721    |
| 2%    | pointwise     | polynomial | —        | 28/0 ✗            | 0.3663     | —                      | 0.6210    | —               | 79/0 ✗    | 0.4291      | 0.6371     | 29.8498      |
| 2%    | weak (no IBP) | spectral   | —        | 6/0 ✗             | 0.5485     | —                      | 0.3476    | —               | 62/0 ✗    | 0.4969      | 0.3239     | 11.3180      |
| 2%    | weak (no IBP) | polynomial | —        | 3/0 ✗             | 0.1957     | —                      | 0.1320    | —               | 9/0 ✗     | 0.0954      | 0.1099     | 0.8842       |
| 2%    | weak (IBP)    | spectral   | —        | 9/0 ✗             | 0.8931     | —                      | 0.3082    | —               | 11/0 ✗    | 0.5803      | 0.3359     | 11.4293      |
| 2%    | weak (IBP)    | polynomial | —        | 1/0 ✗             | 0.6258     | —                      | 0.1858    | —               | 9/0 ✗     | 0.0992      | 0.1098     | 0.9103       |
| 5%    | pointwise     | spectral   | —        | 26/0 ✗            | 0.4316     | —                      | 0.8725    | —               | 89/0 ✗    | 0.4310      | 0.8722     | 2190.0271    |
| 5%    | pointwise     | polynomial | —        | 35/0 ✗            | 0.4222     | —                      | 0.8717    | —               | 93/0 ✗    | 0.4194      | 0.8702     | 47.4252      |
| 5%    | weak (no IBP) | spectral   | —        | 30/0 ✗            | 0.4313     | —                      | 0.5369    | —               | 59/0 ✗    | 0.4459      | 0.4115     | 69.4094      |
| 5%    | weak (no IBP) | polynomial | —        | 8/0 ✗             | 0.3897     | —                      | 0.2940    | —               | 13/0 ✗    | 0.2512      | 0.2611     | 5.4493       |
| 5%    | weak (IBP)    | spectral   | —        | 31/0 ✗            | 0.4804     | —                      | 0.3565    | —               | 70/0 ✗    | 0.4452      | 0.4134     | 70.0482      |
| 5%    | weak (IBP)    | polynomial | —        | 8/0 ✗             | 0.3918     | —                      | 0.2934    | —               | 13/0 ✗    | 0.2474      | 0.2604     | 5.4753       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                     | fp/fn  | E_coeff     | Time (s) | vs. best |
|:-------------------------------|:-------------------------------|:-------|------------:|---------:|---------:|
| Modified Multigraphs + STRidge | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0021      | 0.4      | —        |
| Multigraphs + STRidge          | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0021      | 0.3      | +0.0%    |
| Modified Multigraphs + MSTLS   | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0022      | 0.5      | +3.6%    |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0022      | 0.3      | +3.6%    |
| Multigraphs + RFE              | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0022      | 0.1      | +3.6%    |
| Modified Multigraphs + RFE     | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0022      | 0.2      | +3.6%    |
| PolyJet + STRidge              | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0023      | 1.4      | +6.9%    |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0023      | 4.1      | +8.3%    |
| PolyJet + RFE                  | spectral, weak (no IBP), m=12  | 0/0 ✓  | 0.0023      | 0.2      | +8.3%    |
| WSINDy + STRidge               | polynomial, weak (IBP), m=8    | 63/0 ✗ | 0.4331      | 1.8      | ✗        |
| PDE-FIND + STRidge             | spectral, weak (IBP), m=8      | 20/2 ✗ | 0.4489      | 1.5      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (IBP), m=8      | 1/5 ✗  | 0.4535      | 2.3      | ✗        |
| PDE-FIND + RFE                 | spectral, pointwise            | 1/7 ✗  | 0.4409      | 0.1      | ✗        |
| WSINDy + RFE                   | spectral, pointwise            | 1/7 ✗  | 0.4409      | 0.2      | ✗        |
| WSINDy + MSTLS                 | polynomial, weak (no IBP), m=4 | 4/7 ✗  | 444514.6420 | 2.0      | ✗        |

### 2% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff     | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|------------:|---------:|---------:|
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=12    | 1/0 ✗  | 0.6258      | 0.3      | ✗        |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=8    | 6/0 ✗  | 0.5485      | 0.5      | ✗        |
| Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 8/0 ✗  | 0.1929      | 0.4      | ✗        |
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=12 | 9/0 ✗  | 0.0954      | 4.2      | ✗        |
| PolyJet + STRidge              | polynomial, weak (no IBP), m=8  | 10/0 ✗ | 0.2682      | 2.3      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (IBP), m=12    | 12/0 ✗ | 0.1444      | 0.4      | ✗        |
| WSINDy + STRidge               | polynomial, pointwise           | 79/0 ✗ | 0.4291      | 5.7      | ✗        |
| PDE-FIND + STRidge             | spectral, weak (no IBP), m=4    | 23/2 ✗ | 0.4283      | 5.4      | ✗        |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.5578      | 0.1      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.5578      | 0.1      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=12   | 7/5 ✗  | 0.4447      | 2.5      | ✗        |
| PolyJet + RFE                  | polynomial, weak (no IBP), m=12 | 0/6 ✗  | 0.7043      | 0.2      | ✗        |
| WSINDy + RFE                   | spectral, pointwise             | 1/7 ✗  | 0.4409      | 0.1      | ✗        |
| PDE-FIND + RFE                 | spectral, pointwise             | 1/7 ✗  | 0.4409      | 0.2      | ✗        |
| WSINDy + MSTLS                 | polynomial, weak (no IBP), m=4  | 4/7 ✗  | 438169.9763 | 2.8      | ✗        |

### 5% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff     | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|------------:|---------:|---------:|
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=12 | 8/0 ✗  | 0.3897      | 0.4      | ✗        |
| Modified Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 9/0 ✗  | 0.4100      | 0.3      | ✗        |
| PolyJet + MSTLS                | polynomial, weak (IBP), m=12    | 13/0 ✗ | 0.2474      | 3.8      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 20/0 ✗ | 0.4165      | 0.5      | ✗        |
| Modified Multigraphs + STRidge | spectral, pointwise             | 26/0 ✗ | 0.4316      | 2.0      | ✗        |
| PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 37/0 ✗ | 0.2682      | 2.8      | ✗        |
| WSINDy + STRidge               | polynomial, weak (no IBP), m=4  | 85/0 ✗ | 0.4292      | 4.0      | ✗        |
| PDE-FIND + STRidge             | polynomial, weak (IBP), m=12    | 27/2 ✗ | 0.4373      | 1.3      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/5 ✗  | 0.8574      | 0.1      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=4    | 1/7 ✗  | 0.4409      | 2.6      | ✗        |
| Multigraphs + RFE              | spectral, pointwise             | 1/7 ✗  | 0.5796      | 0.1      | ✗        |
| WSINDy + RFE                   | spectral, pointwise             | 1/7 ✗  | 0.5796      | 0.1      | ✗        |
| PolyJet + RFE                  | spectral, pointwise             | 1/7 ✗  | 0.5796      | 0.2      | ✗        |
| PDE-FIND + RFE                 | spectral, pointwise             | 1/7 ✗  | 0.5796      | 0.1      | ✗        |
| WSINDy + MSTLS                 | polynomial, weak (no IBP), m=4  | 4/7 ✗  | 408758.5378 | 2.0      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.001   | 47/2  | 52       | 0.4274  | 0.4150 | 4.9      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 1/7   | 1        | 0.4409  | 0.4294 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 44/2  | 49       | 0.4292  | 0.6377 | 8.2      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.2      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.01    | 47/2  | 52       | 0.4308  | 0.8720 | 4.6      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 98/0  | 46       | 0.4276  | 0.4152 | 3.4      | 🔴      |
| WSINDy           | Clean | RFE       |         | 1/7   | 1        | 0.4409  | 0.4294 | 0.2      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 94/0  | 45       | 0.4291  | 0.6378 | 3.6      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 89/0  | 44       | 0.4310  | 0.8722 | 3.4      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0047  | 0.0175 | 2.4      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0047  | 0.0175 | 0.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.01    | 83/0  | 90       | 0.4249  | 0.6367 | 12.1     | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.01    | 91/0  | 98       | 0.4309  | 0.8714 | 11.6     | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0047  | 0.0176 | 0.5      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0047  | 0.0176 | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 18/0  | 12       | 0.4266  | 0.6378 | 1.0      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.01    | 28/0  | 16       | 0.4318  | 0.8725 | 1.9      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0047  | 0.0176 | 0.7      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0047  | 0.0176 | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 17/0  | 11       | 0.4267  | 0.6377 | 1.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.4409  | 0.6433 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.01    | 26/0  | 13       | 0.4316  | 0.8725 | 2.0      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.5796  | 0.8979 | 0.2      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 24/2  | 29       | 0.4288  | 0.4158 | 3.4      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 1/7   | 1        | 0.4416  | 0.4294 | 0.2      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.001   | 53/2  | 58       | 0.4294  | 0.6374 | 11.4     | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 1/7   | 1        | 0.4417  | 0.6433 | 0.3      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.01    | 50/2  | 55       | 0.4304  | 0.8716 | 4.8      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 1/7   | 1        | 0.5800  | 0.8979 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 83/0  | 40       | 0.4288  | 0.4157 | 3.6      | 🔴      |
| WSINDy           | Clean | RFE       |         | 1/7   | 1        | 0.4416  | 0.4294 | 0.2      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.001   | 79/0  | 41       | 0.4291  | 0.6371 | 5.7      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 1/7   | 1        | 0.4417  | 0.6433 | 0.2      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 93/0  | 44       | 0.4305  | 0.8718 | 3.7      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 1/7   | 1        | 0.5800  | 0.8979 | 0.2      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0130  | 0.0320 | 2.6      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0128  | 0.0320 | 0.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.001   | 98/0  | 105      | 0.3606  | 0.6182 | 14.3     | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 0.4417  | 0.6433 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.01    | 93/0  | 100      | 0.4194  | 0.8702 | 12.5     | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 0.5800  | 0.8979 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0119  | 0.0325 | 1.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0117  | 0.0325 | 0.2      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.01    | 28/0  | 16       | 0.3663  | 0.6210 | 1.5      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 0.4417  | 0.6433 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 37/0  | 20       | 0.4222  | 0.8717 | 1.9      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 0.5800  | 0.8979 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0119  | 0.0325 | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0117  | 0.0325 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 33/0  | 14       | 0.3663  | 0.6210 | 0.7      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.4417  | 0.6433 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 35/0  | 19       | 0.4222  | 0.8717 | 1.6      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.5800  | 0.8979 | 0.1      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff     | E_data      | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:------------|:------------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 4  | no  | 43/2  | 48       | 0.4284      | 0.4002      | 4.7      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.4173      | 3.1      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5615      | 0.1      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 1.0     | 8  | yes | 20/2  | 25       | 0.4489      | 0.3656      | 1.5      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.75    | 8  | yes | 1/5   | 3        | 0.4535      | 0.3811      | 2.3      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5374      | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 4  | no  | 23/2  | 28       | 0.4283      | 0.4295      | 5.4      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.75    | 12 | no  | 7/5   | 9        | 0.4447      | 0.3361      | 2.5      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5778      | 0.3      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | yes | 38/2  | 43       | 0.4435      | 0.3253      | 2.9      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.75    | 12 | yes | 7/5   | 9        | 0.4448      | 0.3361      | 2.9      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5422      | 0.3      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 34/2  | 39       | 0.4432      | 0.3506      | 1.4      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.5449      | 2.6      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.2      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 34/2  | 39       | 0.4430      | 0.3506      | 1.4      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4409      | 0.4298      | 1.9      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5652      | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 1.0     | 8  | no  | 76/0  | 38       | 0.4334      | 0.3644      | 1.7      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 444689.2112 | 448874.7916 | 1.8      | 🔴      |
| WSINDy           | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5615      | 0.2      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 8  | yes | 98/0  | 46       | 0.4345      | 0.3638      | 1.8      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 662524.5252 | 675775.7402 | 1.4      | 🔴      |
| WSINDy           | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5374      | 0.2      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 12 | no  | 88/0  | 43       | 0.4506      | 0.3247      | 1.3      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 438342.0538 | 436503.6332 | 1.9      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.5778      | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 12 | yes | 88/0  | 43       | 0.4507      | 0.3247      | 1.4      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 658615.2958 | 669444.4221 | 1.5      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5422      | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 8  | no  | 85/0  | 43       | 0.4302      | 0.4131      | 2.0      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 408919.0648 | 381256.0793 | 1.9      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 12 | yes | 90/0  | 45       | 0.4459      | 0.3469      | 1.5      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 644131.2469 | 643358.1219 | 1.6      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5652      | 0.2      | 🔴      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 1.4      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 4.1      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 0.2      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0280      | 0.0157      | 1.4      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 0.0280      | 0.0157      | 4.1      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0280      | 0.0157      | 0.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.1     | 8  | no  | 62/0  | 69       | 0.4969      | 0.3239      | 6.7      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 4  | no  | 12/1  | 18       | 0.4233      | 0.4182      | 7.2      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/6   | 1        | 0.9275      | 0.3548      | 0.1      | 🔴      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | yes | 64/0  | 71       | 0.6538      | 0.2569      | 3.9      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 8  | yes | 11/0  | 18       | 0.5803      | 0.3359      | 5.8      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 0.9275      | 0.3548      | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 8  | no  | 59/0  | 66       | 0.4459      | 0.4115      | 5.7      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 12 | no  | 9/3   | 13       | 0.4430      | 0.3575      | 3.6      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | yes | 70/0  | 77       | 0.4452      | 0.4134      | 5.8      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 12 | yes | 9/3   | 13       | 0.4430      | 0.3575      | 3.9      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5652      | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0021      | 0.0112      | 0.3      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.3      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0279      | 0.0158      | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0279      | 0.0158      | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 4  | no  | 30/0  | 17       | 0.4212      | 0.4176      | 1.6      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 8  | no  | 6/0   | 6        | 0.5485      | 0.3476      | 0.5      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9250      | 0.3445      | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 12 | yes | 31/0  | 17       | 0.7730      | 0.2842      | 0.5      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 12 | yes | 14/0  | 9        | 0.7705      | 0.2872      | 0.4      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9250      | 0.3445      | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.1     | 4  | no  | 30/0  | 17       | 0.4313      | 0.5369      | 1.9      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.5449      | 0.6      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 12 | yes | 31/0  | 17       | 0.4804      | 0.3565      | 0.5      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4409      | 0.4298      | 0.5      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5652      | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0021      | 0.0112      | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0279      | 0.0158      | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0279      | 0.0158      | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 8  | no  | 33/0  | 15       | 0.5413      | 0.3428      | 0.7      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 9/0   | 7        | 0.5490      | 0.3476      | 0.5      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9250      | 0.3445      | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 12 | yes | 33/0  | 18       | 0.8892      | 0.3010      | 0.5      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 9/0   | 7        | 0.8931      | 0.3082      | 0.4      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9250      | 0.3445      | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.01    | 4  | no  | 37/0  | 19       | 0.4313      | 0.5369      | 2.1      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4409      | 0.5449      | 0.8      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5796      | 0.6460      | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 12 | yes | 32/0  | 18       | 0.4729      | 0.3543      | 0.6      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4409      | 0.4298      | 0.6      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5796      | 0.5652      | 0.2      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff     | E_data      | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:------------|:------------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | no  | 48/2  | 53       | 0.4647      | 0.3177      | 2.9      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4416      | 0.4173      | 4.1      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.5615      | 0.2      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 1.0     | 8  | yes | 20/2  | 25       | 0.4495      | 0.3656      | 1.9      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.75    | 8  | yes | 1/5   | 3        | 0.4551      | 0.3811      | 3.0      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.5374      | 0.2      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | no  | 27/2  | 32       | 0.4398      | 0.3266      | 2.7      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4416      | 0.4432      | 5.4      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.5778      | 0.2      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | yes | 27/2  | 32       | 0.4399      | 0.3266      | 2.9      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4416      | 0.3929      | 4.3      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.5422      | 0.2      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 27/2  | 32       | 0.4374      | 0.3491      | 1.2      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.4416      | 0.5449      | 2.8      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.6460      | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 27/2  | 32       | 0.4373      | 0.3491      | 1.3      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.4416      | 0.4298      | 2.1      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.5652      | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 1.0     | 8  | no  | 76/0  | 38       | 0.4336      | 0.3644      | 1.8      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 444514.6420 | 448874.7916 | 2.0      | 🔴      |
| WSINDy           | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.5615      | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 8  | yes | 63/0  | 34       | 0.4331      | 0.3653      | 1.8      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 662264.4417 | 675775.7402 | 1.5      | 🔴      |
| WSINDy           | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.5374      | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 4  | no  | 97/0  | 45       | 0.4293      | 0.4291      | 6.1      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 438169.9763 | 436503.6332 | 2.8      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.5778      | 0.2      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 12 | yes | 88/0  | 43       | 0.4503      | 0.3247      | 1.8      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 658356.7470 | 669444.4221 | 2.2      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.5422      | 0.2      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 4  | no  | 85/0  | 44       | 0.4292      | 0.5353      | 4.0      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 4  | no  | 4/7   | 4        | 408758.5378 | 381256.0793 | 2.0      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.6460      | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 12 | yes | 88/0  | 44       | 0.4455      | 0.3471      | 1.4      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 8  | yes | 4/7   | 4        | 643878.3840 | 643358.1219 | 1.6      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.5652      | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0024      | 0.0112      | 1.2      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 3.7      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0023      | 0.0112      | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0279      | 0.0157      | 1.3      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 0.0279      | 0.0157      | 3.9      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0279      | 0.0157      | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 10.0    | 8  | no  | 10/0  | 17       | 0.2682      | 0.2197      | 2.3      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 12 | no  | 9/0   | 16       | 0.0954      | 0.1099      | 4.2      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/6   | 1        | 0.7043      | 0.3248      | 0.2      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 12 | yes | 13/0  | 20       | 0.1030      | 0.1095      | 1.7      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 12 | yes | 9/0   | 16       | 0.0992      | 0.1098      | 4.2      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 0.7043      | 0.3248      | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 10.0    | 12 | no  | 37/0  | 44       | 0.2682      | 0.2533      | 2.8      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | no  | 13/0  | 20       | 0.2512      | 0.2611      | 3.9      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.6460      | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 10.0    | 12 | yes | 47/0  | 54       | 0.2757      | 0.2529      | 2.9      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | yes | 13/0  | 20       | 0.2474      | 0.2604      | 3.8      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5801      | 0.5652      | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0039      | 0.0141      | 0.7      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.6      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.2      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0277      | 0.0158      | 0.5      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.6      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.2      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 4  | no  | 31/0  | 18       | 0.3347      | 0.3618      | 1.4      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 12 | no  | 7/0   | 5        | 0.1307      | 0.1168      | 0.4      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.5578      | 0.2228      | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | yes | 12/0  | 8        | 0.1444      | 0.1165      | 0.4      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 12 | yes | 7/0   | 5        | 0.1401      | 0.1166      | 0.4      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.5578      | 0.2228      | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | no  | 20/0  | 13       | 0.4165      | 0.2859      | 0.5      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 12 | no  | 8/0   | 6        | 0.3897      | 0.2940      | 0.4      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.6460      | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 8  | yes | 37/0  | 20       | 0.3910      | 0.3708      | 0.7      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 12 | yes | 8/0   | 6        | 0.3918      | 0.2934      | 0.4      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 3/7   | 1        | 0.9854      | 0.5084      | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0024      | 0.0112      | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022      | 0.0112      | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0277      | 0.0158      | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0278      | 0.0158      | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 8  | no  | 11/0  | 9        | 0.3054      | 0.2279      | 0.5      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 12 | no  | 3/0   | 4        | 0.1957      | 0.1320      | 0.4      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.5578      | 0.2228      | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | yes | 8/0   | 7        | 0.1929      | 0.1309      | 0.4      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 1/0   | 4        | 0.6258      | 0.1858      | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.5578      | 0.2228      | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 8  | no  | 30/0  | 17       | 0.3899      | 0.3711      | 0.6      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 12 | no  | 9/0   | 6        | 0.4100      | 0.2941      | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5800      | 0.6460      | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 8  | yes | 37/0  | 20       | 0.5162      | 0.3847      | 0.6      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/5   | 1        | 0.8574      | 0.3690      | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.8574      | 0.3690      | 0.1      | 🔴      |


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
| `sweep_no_ac_ch_ks_wchaos_20260803_170832.csv` | ? min     | ?   | ?       |
