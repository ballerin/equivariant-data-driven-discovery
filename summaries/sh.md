# Swift-Hohenberg equation

The equation under consideration is

$$f_t = (r-1)f - 2\Delta f - \Delta^2 f - f^3, \qquad r = 0.3$$

which can be written in terms of $O(2)$-equivariant terms

| Operator     | True coefficient |
|:-------------|:-----------------|
| $f$          | $r-1 = -0.7$     |
| $\Delta f$   | $-2.0$           |
| $\Delta^2 f$ | $-1.0$           |
| $f^3$        | $-1.0$           |

or standard terms (PolyJet basis)

| Operator            | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|:--------------------|:-----------|:---------|:-------|:---------|:-----------|
| $f$                 |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f^3$               |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xx} f$   |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yy} f$   |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxxx} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxyy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yyyy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |

Swift-Hohenberg contains a 4th-order biharmonic term $\Delta^2 f = f_{xxxx} + 2f_{xxyy} + f_{yyyy}$ which is the main source of difficulty.  All four libraries can represent every term individually, but the equivariant basis encodes the isotropic combination as a single operator, giving it a sparsity advantage

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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best         | Mg E_data | Std best         | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-------------------|:----------|:-----------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | 🟩        | 0/0 ✓             | 0.0273     | multigraph/RFE     | 0.0209    | WSINDy/RFE       | 0/0 ✓     | 0.0274      | 0.0208     | —            |
| Clean | pointwise     | polynomial | 🔴        | 5/0 ✗             | 0.1902     | —                  | 0.1428    | PolyJet/RFE      | 0/0 ✓     | 0.1848      | 0.1484     | —            |
| Clean | weak (no IBP) | spectral   | 🟥        | 0/0 ✓             | 0.0260     | multigraph/RFE     | 0.0138    | WSINDy/STRidge   | 0/0 ✓     | 0.0257      | 0.0138     | —            |
| Clean | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.0260     | multigraph/RFE     | 0.0139    | WSINDy/STRidge   | 0/0 ✓     | 0.0257      | 0.0138     | —            |
| Clean | weak (IBP)    | spectral   | 🟥        | 0/0 ✓             | 0.0265     | multigraph/STRidge | 0.0130    | PDE-FIND/STRidge | 0/0 ✓     | 0.0263      | 0.0129     | —            |
| Clean | weak (IBP)    | polynomial | 🟥        | 0/0 ✓             | 0.0267     | multigraph/STRidge | 0.0130    | WSINDy/STRidge   | 0/0 ✓     | 0.0264      | 0.0129     | —            |
| 2%    | pointwise     | spectral   | —        | 20/0 ✗            | 0.5341     | —                  | 0.6620    | —                | 55/0 ✗    | 0.4476      | 0.6093     | —            |
| 2%    | pointwise     | polynomial | —        | 16/0 ✗            | 0.4884     | —                  | 0.6076    | —                | 18/0 ✗    | 0.7243      | 0.8462     | —            |
| 2%    | weak (no IBP) | spectral   | 🟥        | 0/0 ✓             | 0.0299     | multigraph/STRidge | 0.1183    | PolyJet/STRidge  | 0/0 ✓     | 0.0298      | 0.1183     | —            |
| 2%    | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0294     | multigraph/STRidge | 0.1180    | WSINDy/STRidge   | 0/0 ✓     | 0.0293      | 0.1180     | —            |
| 2%    | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.3024     | multigraph/STRidge | 0.2871    | WSINDy/STRidge   | 0/0 ✓     | 0.3129      | 0.2815     | —            |
| 2%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.2557     | multigraph/STRidge | 0.2871    | WSINDy/STRidge   | 0/0 ✓     | 0.2665      | 0.2815     | —            |
| 5%    | pointwise     | spectral   | —        | 14/0 ✗            | 0.5627     | —                  | 0.8619    | —                | 28/0 ✗    | 0.5768      | 0.8683     | —            |
| 5%    | pointwise     | polynomial | —        | 20/0 ✗            | 0.5093     | —                  | 0.8457    | —                | 41/0 ✗    | 0.7317      | 0.9274     | —            |
| 5%    | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0923     | multigraph/STRidge | 0.2848    | PDE-FIND/STRidge | 0/0 ✓     | 0.0956      | 0.2845     | —            |
| 5%    | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0779     | multigraph/STRidge | 0.2840    | WSINDy/STRidge   | 0/0 ✓     | 0.0811      | 0.2837     | —            |
| 5%    | weak (IBP)    | spectral   | 🔴        | 4/0 ✗             | 0.8446     | —                  | 0.5024    | WSINDy/STRidge   | 0/0 ✓     | 0.8497      | 0.4825     | —            |
| 5%    | weak (IBP)    | polynomial | 🔴        | 14/0 ✗            | 0.7896     | —                  | 0.5789    | WSINDy/STRidge   | 0/0 ✓     | 0.7173      | 0.4825     | —            |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best   |
|:----------------------|:--------------------------------|:------|--------:|---------:|-----------:|
| WSINDy + STRidge      | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0257  | 0.4      | —          |
| PDE-FIND + RFE        | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0259  | 0.1      | +0.6%      |
| WSINDy + RFE          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0259  | 0.1      | +0.6%      |
| PolyJet + RFE         | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0259  | 0.1      | +0.6%      |
| Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0260  | 0.1      | +1.1%      |
| PDE-FIND + STRidge    | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0263  | 0.5      | +2.3%      |
| Multigraphs + STRidge | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0265  | 0.2      | +2.8%      |
| PolyJet + STRidge     | polynomial, weak (IBP), m=8     | 0/0 ✓ | 0.0268  | 1.3      | +4.2%      |
| WSINDy + MSTLS        | polynomial, weak (IBP), m=12    | 0/0 ✓ | 71.2362 | 0.5      | +276867.2% |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=4    | 2/0 ✗ | 0.0301  | 1.9      | ✗          |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=8  | 4/0 ✗ | 0.1248  | 4.0      | ✗          |
| Multigraphs + MSTLS   | polynomial, weak (IBP), m=8     | 8/0 ✗ | 0.0272  | 0.3      | ✗          |

### 2% noise

| Method                | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| WSINDy + STRidge      | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0293  | 0.4      | —        |
| Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0294  | 0.2      | +0.3%    |
| PolyJet + STRidge     | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0298  | 0.9      | +1.8%    |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0298  | 0.5      | +1.8%    |
| WSINDy + RFE          | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0299  | 0.1      | +2.2%    |
| PDE-FIND + RFE        | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0299  | 0.1      | +2.2%    |
| PolyJet + RFE         | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0299  | 0.1      | +2.2%    |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0301  | 0.1      | +2.7%    |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0304  | 1.0      | +3.9%    |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=12 | 6/0 ✗  | 71.1891 | 0.7      | ✗        |
| PolyJet + MSTLS       | spectral, weak (no IBP), m=12   | 8/0 ✗  | 0.3180  | 2.5      | ✗        |
| Multigraphs + MSTLS   | spectral, weak (no IBP), m=12   | 10/0 ✗ | 0.2035  | 0.2      | ✗        |

### 5% noise

| Method                | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0779  | 0.3      | —        |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0799  | 0.1      | +2.5%    |
| WSINDy + STRidge      | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0811  | 0.4      | +4.0%    |
| WSINDy + RFE          | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0832  | 0.1      | +6.8%    |
| PDE-FIND + RFE        | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0832  | 0.1      | +6.8%    |
| PolyJet + RFE         | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0832  | 0.1      | +6.8%    |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0956  | 0.5      | +22.7%   |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12   | 3/0 ✗  | 0.1730  | 1.1      | ✗        |
| PolyJet + STRidge     | spectral, weak (no IBP), m=12   | 5/0 ✗  | 0.1905  | 1.0      | ✗        |
| Multigraphs + MSTLS   | spectral, weak (no IBP), m=12   | 8/0 ✗  | 0.4818  | 0.2      | ✗        |
| PolyJet + MSTLS       | spectral, weak (no IBP), m=12   | 21/0 ✗ | 0.6253  | 2.5      | ✗        |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=12 | 31/0 ✗ | 86.9394 | 0.7      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 0/0   | 7        | 0.0276  | 0.0208 | 0.9      | 🟢      |
| PDE-FIND    | Clean | RFE       |         | 0/0   | 7        | 0.0274  | 0.0208 | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 0.001   | 5/1   | 11       | 0.8356  | 0.8680 | 1.1      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 0/5   | 2        | 0.8437  | 0.8727 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 10/1  | 16       | 0.8260  | 0.9354 | 1.5      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 0/5   | 2        | 0.8436  | 0.9404 | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.001   | 0/0   | 7        | 0.0276  | 0.0208 | 0.7      | 🟢      |
| WSINDy      | Clean | RFE       |         | 0/0   | 7        | 0.0274  | 0.0208 | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.001   | 79/0  | 38       | 0.5160  | 0.6568 | 1.6      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 2/5   | 3        | 0.6324  | 0.7532 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 28/0  | 14       | 0.5768  | 0.8683 | 1.0      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 6/5   | 4        | 0.6164  | 0.8806 | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.01    | 0/0   | 7        | 0.0276  | 0.0208 | 1.8      | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.0274  | 0.0208 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.001   | 55/0  | 62       | 0.4476  | 0.6093 | 4.4      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/5   | 3        | 0.5788  | 0.6951 | 0.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.001   | 29/1  | 35       | 0.4965  | 0.8491 | 3.4      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/5   | 3        | 0.6099  | 0.8872 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.01    | 5/0   | 6        | 0.0278  | 0.0206 | 0.3      | 🔴      |
| Multigraphs | Clean | RFE       |         | 0/0   | 4        | 0.0273  | 0.0209 | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.01    | 20/0  | 12       | 0.5341  | 0.6620 | 0.5      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 2/5   | 3        | 0.5616  | 0.6850 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.001   | 14/0  | 9        | 0.5627  | 0.8619 | 0.5      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 2/5   | 3        | 0.6034  | 0.8860 | 0.1      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 0/0   | 7        | 0.1858  | 0.1484 | 1.0      | 🟢      |
| PDE-FIND    | Clean | RFE       |         | 0/0   | 7        | 0.1848  | 0.1484 | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 0.001   | 18/0  | 25       | 0.7243  | 0.8462 | 1.7      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/6   | 2        | 0.7549  | 0.8649 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 41/0  | 48       | 0.7317  | 0.9274 | 2.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/6   | 2        | 0.7638  | 0.9391 | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.01    | 0/0   | 7        | 0.1858  | 0.1484 | 0.7      | 🟢      |
| WSINDy      | Clean | RFE       |         | 0/0   | 7        | 0.1848  | 0.1484 | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.001   | 41/0  | 21       | 0.5682  | 0.6089 | 1.2      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 2/5   | 3        | 0.6626  | 0.7105 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 77/0  | 36       | 0.6141  | 0.8500 | 1.5      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 2/5   | 3        | 0.6582  | 0.8799 | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.01    | 19/0  | 26       | 0.2384  | 0.0206 | 2.3      | 🔴      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.1848  | 0.1484 | 0.2      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.001   | 44/1  | 50       | 0.4029  | 0.5483 | 4.2      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 2/5   | 4        | 0.5451  | 0.6451 | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.001   | 43/0  | 50       | 0.4331  | 0.8225 | 4.1      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 2/5   | 4        | 0.5526  | 0.8591 | 0.2      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.01    | 5/0   | 6        | 0.1902  | 0.1428 | 0.4      | 🔴      |
| Multigraphs | Clean | RFE       |         | 2/1   | 4        | 0.2336  | 0.2349 | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 0.001   | 16/0  | 11       | 0.4884  | 0.6076 | 0.5      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 13/6  | 7        | 0.5832  | 0.6651 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.001   | 20/0  | 12       | 0.5093  | 0.8457 | 0.5      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 2/5   | 3        | 0.5988  | 0.8640 | 0.1      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data  | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:--------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 8  | no  | 0/0   | 7        | 0.0270  | 0.0127  | 0.6      | 🟢      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 4  | no  | 2/0   | 9        | 0.0301  | 0.0158  | 1.9      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| PDE-FIND    | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 7        | 0.0263  | 0.0129  | 0.5      | 🟢      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 8  | yes | 4/0   | 11       | 0.0340  | 0.0133  | 1.4      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0265  | 0.0129  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0298  | 0.1183  | 0.5      | 🟢      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0304  | 0.1183  | 1.0      | 🟢      |
| PDE-FIND    | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0304  | 0.1183  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | yes | 3/0   | 10       | 0.4297  | 0.2706  | 0.5      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.5     | 12 | yes | 9/0   | 16       | 0.4771  | 0.2671  | 1.1      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.3179  | 0.2815  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0956  | 0.2845  | 0.5      | 🟢      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | no  | 3/0   | 10       | 0.1730  | 0.2825  | 1.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0974  | 0.2845  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 8  | yes | 23/0  | 30       | 1.0886  | 0.6358  | 1.3      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | yes | 3/1   | 9        | 1.0608  | 0.4514  | 1.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | yes | 0/4   | 3        | 0.8261  | 0.6601  | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 7        | 0.0257  | 0.0138  | 0.4      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 8  | no  | 0/0   | 7        | 76.5564 | 70.8534 | 0.7      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| WSINDy      | Clean | STRidge   | 0.01    | 12 | yes | 0/0   | 7        | 0.0263  | 0.0129  | 0.4      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 73.9016 | 67.4830 | 0.5      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0265  | 0.0129  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 7        | 0.0298  | 0.1183  | 0.4      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | no  | 6/0   | 9        | 73.9646 | 63.1279 | 0.6      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0304  | 0.1183  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.3129  | 0.2815  | 0.3      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.75    | 12 | yes | 70/0  | 32       | 78.1935 | 63.5057 | 0.7      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.3179  | 0.2815  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0956  | 0.2845  | 0.4      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.5     | 12 | no  | 31/0  | 18       | 90.3314 | 58.9798 | 0.7      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0974  | 0.2845  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.8497  | 0.4825  | 0.4      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.75    | 8  | yes | 37/0  | 18       | 95.9453 | 63.8236 | 1.0      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 8  | yes | 0/5   | 2        | 1.0654  | 0.7087  | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.1     | 4  | no  | 0/0   | 7        | 0.0271  | 0.0169  | 1.8      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 4  | no  | 5/0   | 12       | 0.0302  | 0.0156  | 5.7      | 🔴      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| PolyJet     | Clean | STRidge   | 0.1     | 8  | yes | 0/0   | 7        | 0.0269  | 0.0149  | 1.4      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | yes | 7/0   | 14       | 0.0344  | 0.0129  | 4.2      | 🔴      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0265  | 0.0129  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0298  | 0.1183  | 0.9      | 🟢      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 12 | no  | 8/0   | 15       | 0.3180  | 0.1173  | 2.5      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0304  | 0.1183  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | yes | 18/0  | 25       | 0.9499  | 0.2406  | 1.2      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 43/0  | 50       | 1.1217  | 0.2223  | 2.5      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.3179  | 0.2815  | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | no  | 5/0   | 12       | 0.1905  | 0.2829  | 1.0      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.5     | 12 | no  | 21/0  | 28       | 0.6253  | 0.2795  | 2.5      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0974  | 0.2845  | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge   | 0.1     | 12 | yes | 27/0  | 34       | 0.7169  | 0.4022  | 1.4      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | yes | 31/0  | 38       | 0.9883  | 0.3759  | 2.7      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 0/5   | 2        | 1.0654  | 0.7087  | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge   | 1.0     | 8  | no  | 0/0   | 4        | 0.0271  | 0.0129  | 0.2      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 8  | no  | 8/0   | 7        | 0.0277  | 0.0124  | 0.3      | 🔴      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0260  | 0.0138  | 0.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0265  | 0.0130  | 0.2      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 12 | yes | 8/0   | 7        | 0.0273  | 0.0124  | 0.2      | 🔴      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0266  | 0.0130  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0299  | 0.1183  | 0.2      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.5     | 12 | no  | 10/0  | 8        | 0.2035  | 0.1176  | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0306  | 0.1183  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.3024  | 0.2871  | 0.2      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.5     | 12 | yes | 15/0  | 10       | 0.7636  | 0.2774  | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.3066  | 0.2871  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 4        | 0.0923  | 0.2848  | 0.2      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.5     | 12 | no  | 8/0   | 7        | 0.4818  | 0.2830  | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0941  | 0.2848  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 0.1     | 12 | yes | 4/0   | 5        | 0.8446  | 0.5024  | 0.2      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.5     | 12 | yes | 10/1  | 7        | 1.2042  | 0.4778  | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 8  | yes | 0/5   | 2        | 1.0654  | 0.7087  | 0.1      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data  | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:--------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 12 | no  | 2/0   | 9        | 0.0261  | 0.0135  | 0.5      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | no  | 4/0   | 11       | 0.0371  | 0.0110  | 1.1      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| PDE-FIND    | Clean | STRidge   | 0.1     | 8  | yes | 0/0   | 7        | 0.0268  | 0.0149  | 0.7      | 🟢      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 8  | yes | 4/0   | 11       | 0.0325  | 0.0134  | 1.4      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0266  | 0.0129  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | no  | 2/0   | 9        | 0.0370  | 0.1178  | 0.6      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | no  | 4/0   | 11       | 0.0511  | 0.1169  | 1.1      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0299  | 0.1180  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | yes | 3/0   | 10       | 0.3718  | 0.2694  | 0.5      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | yes | 23/0  | 30       | 0.4502  | 0.2574  | 1.1      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.2713  | 0.2815  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | no  | 3/0   | 10       | 0.1252  | 0.2826  | 0.5      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | no  | 5/0   | 12       | 0.1655  | 0.2806  | 1.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0832  | 0.2837  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 8  | yes | 17/0  | 24       | 0.8841  | 0.6349  | 1.3      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.5     | 12 | yes | 7/1   | 13       | 0.9078  | 0.4387  | 1.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | yes | 0/4   | 3        | 0.7218  | 0.6601  | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 7        | 0.0257  | 0.0138  | 0.4      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 8  | no  | 0/0   | 7        | 73.8555 | 70.8725 | 0.7      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 7        | 0.0264  | 0.0129  | 0.4      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 71.2362 | 67.4830 | 0.5      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0266  | 0.0129  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 7        | 0.0293  | 0.1180  | 0.4      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | no  | 6/0   | 9        | 71.1891 | 63.1179 | 0.7      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0299  | 0.1180  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 7        | 0.2665  | 0.2815  | 0.4      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.5     | 12 | yes | 66/0  | 31       | 75.8712 | 63.6146 | 0.7      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.2713  | 0.2815  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0811  | 0.2837  | 0.4      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.5     | 12 | no  | 31/0  | 18       | 86.9394 | 59.0284 | 0.7      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0832  | 0.2837  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.7173  | 0.4825  | 0.4      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.75    | 8  | yes | 37/0  | 18       | 98.6183 | 66.4441 | 1.0      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 8  | yes | 0/5   | 2        | 0.8802  | 0.7087  | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 4  | no  | 4/0   | 11       | 0.2187  | 0.0618  | 1.9      | 🔴      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | no  | 4/0   | 11       | 0.1248  | 0.0147  | 4.0      | 🔴      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| PolyJet     | Clean | STRidge   | 0.1     | 8  | yes | 0/0   | 7        | 0.0268  | 0.0149  | 1.3      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | yes | 6/0   | 13       | 0.0297  | 0.0136  | 4.1      | 🔴      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0266  | 0.0129  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | no  | 10/0  | 17       | 2.0901  | 0.1169  | 1.1      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 8  | no  | 30/0  | 37       | 1.2410  | 0.2199  | 4.0      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0299  | 0.1180  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | yes | 23/0  | 30       | 1.5924  | 0.2232  | 1.2      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 43/1  | 49       | 1.6784  | 0.2095  | 2.5      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 1/5   | 3        | 0.8601  | 0.5330  | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | no  | 15/0  | 22       | 2.6815  | 0.2791  | 1.1      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | no  | 42/0  | 49       | 4.2615  | 0.2713  | 2.7      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0832  | 0.2837  | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge   | 0.01    | 12 | yes | 51/0  | 58       | 2.3655  | 0.3428  | 2.0      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | yes | 28/2  | 33       | 1.9773  | 0.3409  | 2.6      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 0/5   | 2        | 0.8802  | 0.7087  | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge   | 1.0     | 8  | no  | 0/0   | 4        | 0.0266  | 0.0178  | 0.2      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 4  | no  | 18/0  | 10       | 0.2299  | 0.0550  | 0.4      | 🔴      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0260  | 0.0139  | 0.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0267  | 0.0130  | 0.2      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 8  | yes | 8/0   | 7        | 0.0272  | 0.0145  | 0.3      | 🔴      |
| Multigraphs | Clean | RFE       |         | 8  | yes | 0/0   | 4        | 0.0268  | 0.0151  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0294  | 0.1180  | 0.2      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.5     | 12 | no  | 20/0  | 12       | 2.5016  | 0.1165  | 0.3      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0301  | 0.1180  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.2557  | 0.2871  | 0.2      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 27/0  | 14       | 0.9855  | 0.2724  | 0.3      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 2/1   | 4        | 0.3542  | 0.2927  | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0779  | 0.2840  | 0.3      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.5     | 12 | no  | 17/0  | 11       | 4.9908  | 0.2796  | 0.3      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0799  | 0.2839  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 0.001   | 8  | yes | 14/0  | 9        | 0.7896  | 0.5789  | 0.3      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | yes | 24/0  | 14       | 1.8536  | 0.4590  | 0.3      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | yes | 2/1   | 4        | 0.7573  | 0.4949  | 0.1      | 🔴      |

---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation sh`.

| Source                              | Wall time | CPU                                      | Threads |
|:------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260801_230839_sh_full.csv` | 76.1 min  | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | 1       |
