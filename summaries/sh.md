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
| Clean | pointwise     | spectral   | 🟩        | 0/0 ✓             | 0.0273     | multigraph/RFE     | 0.0209    | PDE-FIND/RFE     | 0/0 ✓     | 0.0274      | 0.0208     | 0.0338       |
| Clean | pointwise     | polynomial | 🔴        | 5/0 ✗             | 0.1902     | —                  | 0.1428    | PDE-FIND/RFE     | 0/0 ✓     | 0.1848      | 0.1484     | 0.2347       |
| Clean | weak (no IBP) | spectral   | 🟥        | 0/0 ✓             | 0.0260     | multigraph/RFE     | 0.0138    | WSINDy/STRidge   | 0/0 ✓     | 0.0257      | 0.0138     | 0.0288       |
| Clean | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.0260     | multigraph/RFE     | 0.0139    | WSINDy/STRidge   | 0/0 ✓     | 0.0257      | 0.0138     | 0.0290       |
| Clean | weak (IBP)    | spectral   | 🟥        | 0/0 ✓             | 0.0266     | multigraph/RFE     | 0.0130    | PolyJet/STRidge  | 0/0 ✓     | 0.0263      | 0.0129     | 0.0290       |
| Clean | weak (IBP)    | polynomial | 🟥        | 0/0 ✓             | 0.0268     | multigraph/RFE     | 0.0151    | PolyJet/STRidge  | 0/0 ✓     | 0.0264      | 0.0129     | 0.0290       |
| 2%    | pointwise     | spectral   | —        | 12/0 ✗            | 1.0000     | —                  | 0.8696    | —                | 55/0 ✗    | 1.0000      | 0.6093     | 1150.0123    |
| 2%    | pointwise     | polynomial | —        | 14/0 ✗            | 0.9995     | —                  | 0.6653    | —                | 27/0 ✗    | 0.9991      | 0.8454     | 23.6646      |
| 2%    | weak (no IBP) | spectral   | 🟥        | 0/0 ✓             | 0.0096     | multigraph/STRidge | 0.1183    | PDE-FIND/RFE     | 0/0 ✓     | 0.0085      | 0.1183     | 0.1208       |
| 2%    | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.0110     | multigraph/RFE     | 0.1180    | WSINDy/RFE       | 0/0 ✓     | 0.0098      | 0.1180     | 0.1204       |
| 2%    | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.2130     | multigraph/RFE     | 0.2871    | PolyJet/RFE      | 0/0 ✓     | 0.2447      | 0.2815     | 0.3186       |
| 2%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.2148     | multigraph/STRidge | 0.2871    | PDE-FIND/RFE     | 0/0 ✓     | 0.2563      | 0.2815     | 0.3186       |
| 5%    | pointwise     | spectral   | —        | 18/0 ✗            | 1.0000     | —                  | 0.8608    | —                | 17/0 ✗    | 1.0000      | 0.8600     | 2032.3997    |
| 5%    | pointwise     | polynomial | —        | 9/0 ✗             | 1.0001     | —                  | 0.8738    | —                | 24/0 ✗    | 1.0000      | 0.9276     | 41.8218      |
| 5%    | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0472     | multigraph/RFE     | 0.2848    | PDE-FIND/STRidge | 0/0 ✓     | 0.0559      | 0.2845     | 0.2874       |
| 5%    | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0450     | multigraph/RFE     | 0.2839    | WSINDy/STRidge   | 0/0 ✓     | 0.0553      | 0.2837     | 0.2864       |
| 5%    | weak (IBP)    | spectral   | 🔴        | 4/0 ✗             | 0.6375     | —                  | 0.5024    | WSINDy/STRidge   | 0/0 ✓     | 0.6666      | 0.4825     | 0.7756       |
| 5%    | weak (IBP)    | polynomial | 🔴        | 14/0 ✗            | 0.9954     | —                  | 0.5789    | WSINDy/STRidge   | 0/0 ✓     | 0.6920      | 0.4825     | 0.7756       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best   |
|:----------------------|:--------------------------------|:------|--------:|---------:|-----------:|
| WSINDy + STRidge      | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0257  | 0.2      | —          |
| PDE-FIND + RFE        | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0259  | 0.2      | +0.6%      |
| WSINDy + RFE          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0259  | 0.1      | +0.6%      |
| PolyJet + RFE         | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0259  | 0.4      | +0.6%      |
| PDE-FIND + STRidge    | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0259  | 0.3      | +0.7%      |
| Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0260  | 0.1      | +1.1%      |
| PolyJet + STRidge     | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0263  | 3.2      | +2.3%      |
| Multigraphs + STRidge | polynomial, weak (IBP), m=8     | 0/0 ✓ | 0.0269  | 1.6      | +4.6%      |
| WSINDy + MSTLS        | polynomial, weak (IBP), m=12    | 0/0 ✓ | 71.2362 | 0.8      | +276867.2% |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=4    | 2/0 ✗ | 0.0301  | 4.2      | ✗          |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=8  | 4/0 ✗ | 0.1248  | 36.4     | ✗          |
| Multigraphs + MSTLS   | polynomial, weak (IBP), m=8     | 8/0 ✗ | 0.0272  | 0.9      | ✗          |

### 2% noise

| Method                | Best conf.                    | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:------------------------------|:-------|--------:|---------:|---------:|
| PDE-FIND + RFE        | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0085  | 0.3      | —        |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0085  | 2.4      | +0.0%    |
| WSINDy + RFE          | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0085  | 0.1      | +0.0%    |
| PolyJet + RFE         | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0085  | 0.4      | +0.0%    |
| PolyJet + STRidge     | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0095  | 3.1      | +11.8%   |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0095  | 0.2      | +11.8%   |
| WSINDy + STRidge      | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0095  | 0.2      | +11.8%   |
| Multigraphs + STRidge | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0096  | 0.1      | +13.2%   |
| Multigraphs + RFE     | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0100  | 0.1      | +17.1%   |
| WSINDy + MSTLS        | spectral, weak (no IBP), m=12 | 6/0 ✗  | 6.0366  | 0.9      | ✗        |
| PolyJet + MSTLS       | spectral, weak (no IBP), m=12 | 8/0 ✗  | 0.0162  | 33.6     | ✗        |
| Multigraphs + MSTLS   | spectral, weak (no IBP), m=12 | 10/0 ✗ | 0.0108  | 0.2      | ✗        |

### 5% noise

| Method                | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0450  | 0.1      | —        |
| Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0456  | 0.1      | +1.5%    |
| WSINDy + STRidge      | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0553  | 0.2      | +23.0%   |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0559  | 0.2      | +24.4%   |
| PDE-FIND + RFE        | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0584  | 0.2      | +29.9%   |
| WSINDy + RFE          | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0584  | 0.1      | +29.9%   |
| PolyJet + RFE         | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0584  | 0.4      | +29.9%   |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12   | 2/0 ✗  | 0.1000  | 2.4      | ✗        |
| PolyJet + STRidge     | spectral, weak (no IBP), m=12   | 3/0 ✗  | 0.0895  | 3.6      | ✗        |
| Multigraphs + MSTLS   | spectral, weak (no IBP), m=12   | 8/0 ✗  | 0.0532  | 0.3      | ✗        |
| PolyJet + MSTLS       | spectral, weak (no IBP), m=12   | 21/0 ✗ | 0.1237  | 35.8     | ✗        |
| WSINDy + MSTLS        | spectral, weak (no IBP), m=12   | 31/0 ✗ | 4.8985  | 1.0      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 0/0   | 7        | 0.0276  | 0.0208 | 5.1      | 🟢      |
| PDE-FIND    | Clean | RFE       |         | 0/0   | 7        | 0.0274  | 0.0208 | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 0.01    | 4/1   | 10       | 1.0000  | 0.8582 | 5.4      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 0/5   | 2        | 1.0000  | 0.8727 | 0.4      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 28/0  | 46       | 1.0000  | 0.9345 | 6.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 0/5   | 2        | 1.0000  | 0.9404 | 0.4      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.01    | 0/0   | 7        | 0.0276  | 0.0208 | 3.7      | 🟢      |
| WSINDy      | Clean | RFE       |         | 0/0   | 7        | 0.0274  | 0.0208 | 0.3      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.001   | 60/0  | 31       | 1.0000  | 0.6569 | 4.0      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 2/5   | 3        | 1.0001  | 0.7532 | 0.3      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 64/0  | 31       | 1.0000  | 0.8664 | 4.8      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 6/5   | 4        | 1.0000  | 0.8806 | 0.3      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0276  | 0.0208 | 12.0     | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.0274  | 0.0208 | 0.7      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.001   | 55/0  | 62       | 1.0000  | 0.6093 | 16.6     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/5   | 3        | 1.0000  | 0.6951 | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.001   | 17/0  | 24       | 1.0000  | 0.8600 | 14.9     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/5   | 3        | 1.0000  | 0.8872 | 0.7      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 4        | 0.0275  | 0.0209 | 1.6      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 4        | 0.0273  | 0.0209 | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.001   | 12/0  | 8        | 1.0000  | 0.8696 | 2.2      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 2/5   | 3        | 1.0000  | 0.6850 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.01    | 18/0  | 11       | 1.0000  | 0.8608 | 1.7      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 2/5   | 3        | 1.0000  | 0.8860 | 0.1      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.1858  | 0.1484 | 4.9      | 🟢      |
| PDE-FIND    | Clean | RFE       |         | 0/0   | 7        | 0.1848  | 0.1484 | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 0.001   | 27/0  | 44       | 0.9991  | 0.8454 | 5.6      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/6   | 2        | 1.0009  | 0.8649 | 0.4      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 24/0  | 41       | 1.0000  | 0.9276 | 5.5      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/6   | 2        | 1.0011  | 0.9391 | 0.4      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.1858  | 0.1484 | 3.7      | 🟢      |
| WSINDy      | Clean | RFE       |         | 0/0   | 7        | 0.1848  | 0.1484 | 0.3      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.001   | 64/0  | 31       | 0.9995  | 0.6076 | 4.8      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 2/5   | 3        | 1.0002  | 0.7105 | 0.3      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 75/0  | 35       | 1.0001  | 0.8500 | 4.5      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 2/5   | 3        | 1.0004  | 0.8799 | 0.3      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.1858  | 0.1484 | 11.2     | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.1848  | 0.1484 | 0.7      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.001   | 76/0  | 83       | 0.9994  | 0.5419 | 18.9     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 2/5   | 4        | 1.0003  | 0.6451 | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.001   | 37/0  | 44       | 1.0001  | 0.8292 | 14.9     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 2/5   | 4        | 1.0005  | 0.8591 | 0.8      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 5/0   | 6        | 0.1902  | 0.1428 | 1.6      | 🔴      |
| Multigraphs | Clean | RFE       |         | 2/1   | 4        | 0.2336  | 0.2349 | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 0.001   | 14/0  | 9        | 0.9995  | 0.6653 | 1.6      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 13/6  | 7        | 0.9983  | 0.6651 | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.001   | 9/0   | 7        | 1.0001  | 0.8738 | 2.0      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 2/5   | 3        | 0.9999  | 0.8640 | 0.1      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data  | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:--------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 8  | no  | 0/0   | 7        | 0.0270  | 0.0127  | 5.2      | 🟢      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 4  | no  | 2/0   | 9        | 0.0301  | 0.0158  | 4.2      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.3      | 🟢      |
| PDE-FIND    | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0263  | 0.0129  | 0.2      | 🟢      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 8  | yes | 4/0   | 11       | 0.0340  | 0.0133  | 3.7      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0265  | 0.0129  | 0.2      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0095  | 0.1183  | 0.2      | 🟢      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0085  | 0.1183  | 2.4      | 🟢      |
| PDE-FIND    | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0085  | 0.1183  | 0.3      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.2457  | 0.2815  | 0.2      | 🟢      |
| PDE-FIND    | 2%    | MSTLS     | 0.5     | 12 | yes | 7/0   | 16       | 0.3246  | 0.2671  | 3.5      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.2447  | 0.2815  | 0.2      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0559  | 0.2845  | 0.2      | 🟢      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | no  | 2/0   | 10       | 0.1000  | 0.2825  | 2.4      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0584  | 0.2845  | 0.2      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 12 | yes | 9/0   | 18       | 0.7522  | 0.4342  | 0.3      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | yes | 2/1   | 9        | 0.7539  | 0.4514  | 2.6      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | yes | 0/4   | 3        | 0.9114  | 0.6601  | 0.3      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 7        | 0.0257  | 0.0138  | 0.2      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 8  | no  | 0/0   | 7        | 76.5564 | 70.8534 | 1.8      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 7        | 0.0263  | 0.0129  | 0.2      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 73.9016 | 67.4830 | 0.7      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0265  | 0.0129  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0095  | 0.1183  | 0.2      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | no  | 6/0   | 9        | 6.0366  | 63.1279 | 0.9      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0085  | 0.1183  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.2457  | 0.2815  | 0.2      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.75    | 12 | yes | 70/0  | 32       | 3.9330  | 63.5057 | 1.0      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.2447  | 0.2815  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0559  | 0.2845  | 0.2      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.5     | 12 | no  | 31/0  | 18       | 4.8985  | 58.9798 | 1.0      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0584  | 0.2845  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.6666  | 0.4825  | 0.2      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.75    | 8  | yes | 37/0  | 18       | 1.0110  | 63.8236 | 2.3      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | yes | 0/5   | 2        | 1.0000  | 0.6044  | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 8  | no  | 0/0   | 7        | 0.0270  | 0.0127  | 12.6     | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 4  | no  | 5/0   | 12       | 0.0302  | 0.0156  | 35.2     | 🔴      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.4      | 🟢      |
| PolyJet     | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0263  | 0.0129  | 3.2      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | yes | 7/0   | 14       | 0.0344  | 0.0129  | 37.8     | 🔴      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0265  | 0.0129  | 0.4      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0095  | 0.1183  | 3.1      | 🟢      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 12 | no  | 8/0   | 15       | 0.0162  | 0.1173  | 33.6     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0085  | 0.1183  | 0.4      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | yes | 18/0  | 25       | 0.4687  | 0.2406  | 4.1      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 43/0  | 50       | 0.5594  | 0.2223  | 36.2     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.2447  | 0.2815  | 0.4      | 🟢      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | no  | 3/0   | 10       | 0.0895  | 0.2829  | 3.6      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.5     | 12 | no  | 21/0  | 28       | 0.1237  | 0.2795  | 35.8     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0584  | 0.2845  | 0.4      | 🟢      |
| PolyJet     | 5%    | STRidge   | 0.1     | 12 | yes | 27/0  | 34       | 0.7909  | 0.4022  | 5.5      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | yes | 31/0  | 38       | 0.8280  | 0.3759  | 36.5     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 0/5   | 2        | 1.0000  | 0.7087  | 0.7      | 🔴      |
| Multigraphs | Clean | STRidge   | 1.0     | 4  | no  | 0/0   | 4        | 0.0270  | 0.0170  | 1.6      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 8  | no  | 8/0   | 7        | 0.0277  | 0.0124  | 0.8      | 🔴      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0260  | 0.0138  | 0.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 1.0     | 8  | yes | 0/0   | 4        | 0.0269  | 0.0151  | 1.6      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 12 | yes | 8/0   | 7        | 0.0273  | 0.0124  | 0.3      | 🔴      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0266  | 0.0130  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0096  | 0.1183  | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.5     | 12 | no  | 10/0  | 8        | 0.0108  | 0.1176  | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0100  | 0.1183  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.2143  | 0.2871  | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.5     | 12 | yes | 15/0  | 10       | 0.2713  | 0.2774  | 0.3      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.2130  | 0.2871  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 4        | 0.0480  | 0.2848  | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.5     | 12 | no  | 8/0   | 7        | 0.0532  | 0.2830  | 0.3      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0472  | 0.2848  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 0.1     | 12 | yes | 4/0   | 5        | 0.6375  | 0.5024  | 0.1      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.5     | 12 | yes | 10/1  | 7        | 0.6747  | 0.4778  | 0.3      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 8  | yes | 0/5   | 2        | 1.0000  | 0.7087  | 0.1      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data  | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:--------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 9        | 0.0259  | 0.0135  | 0.3      | 🟢      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 8  | no  | 4/0   | 13       | 0.0361  | 0.0157  | 3.8      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.2      | 🟢      |
| PDE-FIND    | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0264  | 0.0129  | 0.2      | 🟢      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | yes | 2/0   | 11       | 0.0358  | 0.0107  | 2.6      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0266  | 0.0129  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 8  | no  | 0/0   | 7        | 0.2780  | 0.2910  | 4.9      | 🟢      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | no  | 3/0   | 11       | 0.0085  | 0.1169  | 2.8      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0098  | 0.1180  | 0.3      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | yes | 2/0   | 10       | 0.3336  | 0.2694  | 0.3      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.75    | 12 | yes | 18/0  | 30       | 0.4052  | 0.2574  | 3.1      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.2563  | 0.2815  | 0.3      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | no  | 2/0   | 10       | 0.0814  | 0.2826  | 0.3      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | no  | 4/0   | 12       | 0.1113  | 0.2806  | 3.3      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0589  | 0.2837  | 0.2      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 12 | yes | 16/0  | 26       | 0.7946  | 0.4250  | 0.4      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.5     | 12 | yes | 3/1   | 13       | 0.7870  | 0.4387  | 2.8      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | yes | 0/4   | 3        | 0.8577  | 0.6601  | 0.3      | 🔴      |
| WSINDy      | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0257  | 0.0138  | 0.2      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 8  | no  | 0/0   | 7        | 73.8555 | 70.8725 | 2.0      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 7        | 0.0264  | 0.0129  | 0.2      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 71.2362 | 67.4830 | 0.8      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0266  | 0.0129  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0107  | 0.1180  | 0.2      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | no  | 6/0   | 9        | 11.0707 | 63.1179 | 0.9      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0098  | 0.1180  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.2570  | 0.2815  | 0.2      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.5     | 12 | yes | 66/0  | 31       | 7.8182  | 63.6146 | 1.0      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.2563  | 0.2815  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0553  | 0.2837  | 0.2      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.5     | 12 | no  | 31/0  | 18       | 9.0052  | 59.0284 | 0.9      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0589  | 0.2837  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.6920  | 0.4825  | 0.3      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.75    | 8  | yes | 37/0  | 18       | 2.1074  | 66.4441 | 2.3      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | yes | 0/5   | 2        | 0.9999  | 0.6044  | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 4  | no  | 4/0   | 11       | 0.2187  | 0.0618  | 11.5     | 🔴      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | no  | 4/0   | 11       | 0.1248  | 0.0147  | 36.4     | 🔴      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.4      | 🟢      |
| PolyJet     | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0264  | 0.0129  | 3.6      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | yes | 6/0   | 13       | 0.0297  | 0.0136  | 35.4     | 🔴      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0266  | 0.0129  | 0.4      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | no  | 10/0  | 17       | 0.0793  | 0.1169  | 3.3      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.5     | 12 | no  | 30/0  | 37       | 0.1504  | 0.1153  | 39.3     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0098  | 0.1180  | 0.4      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | yes | 22/0  | 29       | 0.5711  | 0.2258  | 4.0      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 43/1  | 49       | 0.6494  | 0.2095  | 37.6     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 1/5   | 3        | 0.9033  | 0.5330  | 0.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | no  | 15/0  | 22       | 0.1308  | 0.2791  | 3.6      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | no  | 42/0  | 49       | 0.3284  | 0.2713  | 39.5     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0589  | 0.2837  | 0.4      | 🟢      |
| PolyJet     | 5%    | STRidge   | 0.01    | 12 | yes | 51/0  | 58       | 0.9042  | 0.3428  | 7.5      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | yes | 28/2  | 33       | 0.9413  | 0.3409  | 37.4     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 0/5   | 2        | 0.9999  | 0.7087  | 0.6      | 🔴      |
| Multigraphs | Clean | STRidge   | 1.0     | 12 | no  | 3/0   | 5        | 0.0258  | 0.0139  | 0.1      | 🔴      |
| Multigraphs | Clean | MSTLS     | 0.5     | 4  | no  | 18/0  | 10       | 0.2299  | 0.0550  | 0.9      | 🔴      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0260  | 0.0139  | 0.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 1.0     | 8  | yes | 0/0   | 4        | 0.0269  | 0.0151  | 1.6      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 8  | yes | 8/0   | 7        | 0.0272  | 0.0145  | 0.9      | 🔴      |
| Multigraphs | Clean | RFE       |         | 8  | yes | 0/0   | 4        | 0.0268  | 0.0151  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | no  | 5/0   | 6        | 0.0103  | 0.1180  | 0.1      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.5     | 12 | no  | 20/0  | 12       | 0.0983  | 0.1165  | 0.3      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0110  | 0.1180  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 10.0    | 12 | yes | 0/0   | 4        | 0.2148  | 0.2871  | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 27/0  | 14       | 0.3093  | 0.2724  | 0.3      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 2/1   | 4        | 0.2703  | 0.2927  | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 4        | 0.0456  | 0.2840  | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.5     | 12 | no  | 17/0  | 11       | 0.1196  | 0.2796  | 0.3      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0450  | 0.2839  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 0.01    | 8  | yes | 14/0  | 9        | 0.9954  | 0.5789  | 1.5      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | yes | 24/0  | 14       | 0.7069  | 0.4590  | 0.3      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | yes | 2/1   | 4        | 0.6571  | 0.4949  | 0.1      | 🔴      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ) |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:--------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0314       | 10372.4 |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0291       | 15893.5 |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0288       | 28641.5 |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.1275       | 10429.5 |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0316       | 15964.3 |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0290       | 28227.6 |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0301       | 11248.4 |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0290       | 17099.7 |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.0301       | 11323.6 |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0290       | 17129.6 |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 11.8754      | 692.9   |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 0.3767       | 1432.1  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 0.1208       | 3350.2  |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 4.1394       | 8018.2  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.3327       | 13590.1 |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.1204       | 23171.2 |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 3.7798       | 842.6   |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 0.3186       | 1499.3  |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 3.7798       | 9167.5  |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.3186       | 14775.9 |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 28.1988      | 744.2   |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 0.9184       | 1578.6  |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 0.2874       | 4258.3  |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 9.8233       | 3075.9  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 0.8096       | 6209.9  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 0.2864       | 11568.1 |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 9.1095       | 1015.5  |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 0.7756       | 1695.2  |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 9.1095       | 3797.9  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 0.7756       | 6686.9  |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation sh`.

| Source                              | Wall time | CPU                                      | Threads |
|:------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260729_200349_sh_full.csv` | 423.7 min | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | not set |
