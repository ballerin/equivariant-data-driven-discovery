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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best             | Mg E_data | Std best         | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-----------------------|:----------|:-----------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | 🟩        | 0/0 ✓             | 0.0273     | multigraph/RFE         | 0.0209    | PolyJet/RFE      | 0/0 ✓     | 0.0274      | 0.0208     | 0.0338       |
| Clean | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.1826     | mod_multigraph/RFE     | 0.1511    | PolyJet/RFE      | 0/0 ✓     | 0.1848      | 0.1484     | 0.2347       |
| Clean | weak (no IBP) | spectral   | 🟥        | 0/0 ✓             | 0.0259     | mod_multigraph/STRidge | 0.0138    | WSINDy/STRidge   | 0/0 ✓     | 0.0257      | 0.0138     | 0.0288       |
| Clean | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.0259     | mod_multigraph/STRidge | 0.0139    | WSINDy/STRidge   | 0/0 ✓     | 0.0257      | 0.0138     | 0.0290       |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0224     | mod_multigraph/STRidge | 0.0153    | PDE-FIND/STRidge | 0/0 ✓     | 0.0230      | 0.0138     | 0.0268       |
| Clean | weak (IBP)    | polynomial | 🟥        | 0/0 ✓             | 0.0231     | mod_multigraph/STRidge | 0.0153    | PDE-FIND/STRidge | 0/0 ✓     | 0.0225      | 0.0138     | 0.0268       |
| 2%    | pointwise     | spectral   | —        | 12/0 ✗            | 1.0000     | —                      | 0.8696    | —                | 55/0 ✗    | 1.0000      | 0.6093     | 1150.0123    |
| 2%    | pointwise     | polynomial | —        | 9/0 ✗             | 0.9994     | —                      | 0.6679    | —                | 27/0 ✗    | 0.9991      | 0.8454     | 23.6646      |
| 2%    | weak (no IBP) | spectral   | 🟥        | 0/0 ✓             | 0.0096     | multigraph/STRidge     | 0.1183    | PDE-FIND/RFE     | 0/0 ✓     | 0.0085      | 0.1183     | 0.1208       |
| 2%    | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.0108     | mod_multigraph/STRidge | 0.1180    | PolyJet/RFE      | 0/0 ✓     | 0.0098      | 0.1180     | 0.1204       |
| 2%    | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0079     | multigraph/STRidge     | 0.1185    | PolyJet/RFE      | 0/0 ✓     | 0.0099      | 0.1184     | 0.1202       |
| 2%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0080     | mod_multigraph/STRidge | 0.1185    | PolyJet/RFE      | 0/0 ✓     | 0.0115      | 0.1184     | 0.1202       |
| 5%    | pointwise     | spectral   | —        | 13/0 ✗            | 1.0000     | —                      | 0.9396    | —                | 17/0 ✗    | 1.0000      | 0.8600     | 2032.3997    |
| 5%    | pointwise     | polynomial | —        | 7/0 ✗             | 1.0001     | —                      | 0.8740    | —                | 24/0 ✗    | 1.0000      | 0.9276     | 41.8218      |
| 5%    | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0472     | multigraph/RFE         | 0.2848    | PDE-FIND/STRidge | 0/0 ✓     | 0.0559      | 0.2845     | 0.2874       |
| 5%    | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0450     | multigraph/RFE         | 0.2839    | WSINDy/STRidge   | 0/0 ✓     | 0.0553      | 0.2837     | 0.2864       |
| 5%    | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0487     | multigraph/RFE         | 0.2847    | PDE-FIND/STRidge | 0/0 ✓     | 0.0548      | 0.2846     | 0.2869       |
| 5%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0489     | multigraph/RFE         | 0.2847    | WSINDy/STRidge   | 0/0 ✓     | 0.0559      | 0.2846     | 0.2869       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                   | fp/fn  | E_coeff | Time (s) | vs. best   |
|:-------------------------------|:-----------------------------|:-------|--------:|---------:|-----------:|
| Modified Multigraphs + STRidge | spectral, weak (IBP), m=12   | 0/0 ✓  | 0.0224  | 0.1      | —          |
| Multigraphs + RFE              | spectral, weak (IBP), m=12   | 0/0 ✓  | 0.0225  | 0.1      | +0.3%      |
| Modified Multigraphs + RFE     | spectral, weak (IBP), m=12   | 0/0 ✓  | 0.0225  | 0.1      | +0.3%      |
| PDE-FIND + STRidge             | polynomial, weak (IBP), m=12 | 0/0 ✓  | 0.0225  | 0.2      | +0.5%      |
| WSINDy + STRidge               | polynomial, weak (IBP), m=12 | 0/0 ✓  | 0.0225  | 0.2      | +0.5%      |
| PolyJet + RFE                  | polynomial, weak (IBP), m=12 | 0/0 ✓  | 0.0227  | 0.2      | +1.2%      |
| PDE-FIND + RFE                 | polynomial, weak (IBP), m=12 | 0/0 ✓  | 0.0227  | 0.1      | +1.2%      |
| WSINDy + RFE                   | polynomial, weak (IBP), m=12 | 0/0 ✓  | 0.0227  | 0.1      | +1.2%      |
| PolyJet + STRidge              | spectral, weak (no IBP), m=8 | 0/0 ✓  | 0.0270  | 0.5      | +20.5%     |
| Multigraphs + STRidge          | spectral, weak (no IBP), m=4 | 0/0 ✓  | 0.0270  | 0.1      | +20.8%     |
| WSINDy + MSTLS                 | polynomial, weak (IBP), m=8  | 0/0 ✓  | 71.4884 | 0.5      | +319211.2% |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=4 | 2/0 ✗  | 0.0301  | 1.2      | ✗          |
| PolyJet + MSTLS                | polynomial, weak (IBP), m=8  | 2/0 ✗  | 0.1729  | 2.0      | ✗          |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=8 | 8/0 ✗  | 0.0277  | 0.3      | ✗          |
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=8  | 17/0 ✗ | 0.1794  | 0.3      | ✗          |

### 2% noise

| Method                         | Best conf.                    | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + STRidge          | spectral, weak (IBP), m=12    | 0/0 ✓  | 0.0079  | 0.1      | —        |
| Modified Multigraphs + STRidge | spectral, weak (IBP), m=12    | 0/0 ✓  | 0.0079  | 0.1      | +0.0%    |
| Multigraphs + RFE              | spectral, weak (IBP), m=12    | 0/0 ✓  | 0.0083  | 0.1      | +4.5%    |
| Modified Multigraphs + RFE     | spectral, weak (IBP), m=12    | 0/0 ✓  | 0.0083  | 0.1      | +4.5%    |
| PDE-FIND + RFE                 | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0085  | 0.1      | +7.8%    |
| PolyJet + RFE                  | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0085  | 0.1      | +7.8%    |
| WSINDy + RFE                   | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0085  | 0.1      | +7.8%    |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0085  | 0.5      | +7.8%    |
| PDE-FIND + STRidge             | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0095  | 0.2      | +20.4%   |
| WSINDy + STRidge               | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0095  | 0.2      | +20.4%   |
| PolyJet + STRidge              | spectral, weak (no IBP), m=12 | 0/0 ✓  | 0.0095  | 0.4      | +20.4%   |
| Modified Multigraphs + MSTLS   | spectral, weak (IBP), m=12    | 6/0 ✗  | 0.0090  | 0.2      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=12 | 6/0 ✗  | 6.0366  | 0.4      | ✗        |
| PolyJet + MSTLS                | spectral, weak (IBP), m=12    | 8/0 ✗  | 0.0150  | 1.3      | ✗        |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=12 | 10/0 ✗ | 0.0108  | 0.2      | ✗        |

### 5% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0450  | 0.1      | —        |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0450  | 0.1      | +0.0%    |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0456  | 0.1      | +1.5%    |
| Modified Multigraphs + STRidge | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0480  | 0.1      | +6.8%    |
| PDE-FIND + STRidge             | spectral, weak (IBP), m=12      | 0/0 ✓  | 0.0548  | 0.2      | +21.9%   |
| WSINDy + STRidge               | spectral, weak (IBP), m=12      | 0/0 ✓  | 0.0548  | 0.2      | +21.9%   |
| PDE-FIND + RFE                 | spectral, weak (IBP), m=12      | 0/0 ✓  | 0.0573  | 0.1      | +27.4%   |
| PolyJet + RFE                  | spectral, weak (IBP), m=12      | 0/0 ✓  | 0.0573  | 0.1      | +27.4%   |
| WSINDy + RFE                   | spectral, weak (IBP), m=12      | 0/0 ✓  | 0.0573  | 0.1      | +27.4%   |
| Modified Multigraphs + MSTLS   | spectral, weak (no IBP), m=12   | 2/0 ✗  | 0.0468  | 0.2      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (IBP), m=12      | 2/0 ✗  | 0.0992  | 0.5      | ✗        |
| PolyJet + STRidge              | spectral, weak (IBP), m=12      | 3/0 ✗  | 0.0887  | 0.4      | ✗        |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=12   | 8/0 ✗  | 0.0532  | 0.2      | ✗        |
| PolyJet + MSTLS                | spectral, weak (IBP), m=12      | 21/0 ✗ | 0.1229  | 1.2      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=12   | 31/0 ✗ | 4.8985  | 0.5      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 0/0   | 7        | 0.0276  | 0.0208 | 0.3      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 0/0   | 7        | 0.0274  | 0.0208 | 0.1      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 0.01    | 4/1   | 10       | 1.0000  | 0.8582 | 0.3      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 0/5   | 2        | 1.0000  | 0.8727 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.001   | 28/0  | 46       | 1.0000  | 0.9345 | 0.4      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 0/5   | 2        | 1.0000  | 0.9404 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.01    | 0/0   | 7        | 0.0276  | 0.0208 | 0.3      | 🟢      |
| WSINDy           | Clean | RFE       |         | 0/0   | 7        | 0.0274  | 0.0208 | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.001   | 60/0  | 31       | 1.0000  | 0.6569 | 0.4      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 2/5   | 3        | 1.0001  | 0.7532 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.001   | 64/0  | 31       | 1.0000  | 0.8664 | 0.3      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 6/5   | 4        | 1.0000  | 0.8806 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0276  | 0.0208 | 0.6      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0274  | 0.0208 | 0.3      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.001   | 55/0  | 62       | 1.0000  | 0.6093 | 1.1      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/5   | 3        | 1.0000  | 0.6951 | 0.3      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.001   | 17/0  | 24       | 1.0000  | 0.8600 | 0.9      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/5   | 3        | 1.0000  | 0.8872 | 0.3      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 0/0   | 4        | 0.0275  | 0.0209 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 4        | 0.0273  | 0.0209 | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.001   | 12/0  | 8        | 1.0000  | 0.8696 | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 2/5   | 3        | 1.0000  | 0.6850 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 13/0  | 9        | 1.0000  | 0.9396 | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 2/5   | 3        | 1.0000  | 0.8860 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 4        | 0.0275  | 0.0209 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 4        | 0.0273  | 0.0209 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.001   | 16/0  | 6        | 1.0000  | 0.8676 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 2/5   | 3        | 1.0000  | 0.6850 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 30/0  | 18       | 1.0000  | 0.8603 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 2/5   | 3        | 1.0000  | 0.8860 | 0.1      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.1858  | 0.1484 | 0.3      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 0/0   | 7        | 0.1848  | 0.1484 | 0.1      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 0.001   | 27/0  | 44       | 0.9991  | 0.8454 | 0.5      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 1/6   | 2        | 1.0009  | 0.8649 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.001   | 24/0  | 41       | 1.0000  | 0.9276 | 0.4      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 1/6   | 2        | 1.0011  | 0.9391 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.1858  | 0.1484 | 0.3      | 🟢      |
| WSINDy           | Clean | RFE       |         | 0/0   | 7        | 0.1848  | 0.1484 | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.001   | 64/0  | 31       | 0.9995  | 0.6076 | 0.4      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 2/5   | 3        | 1.0002  | 0.7105 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.001   | 75/0  | 35       | 1.0001  | 0.8500 | 0.4      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 2/5   | 3        | 1.0004  | 0.8799 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.1858  | 0.1484 | 0.7      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.1848  | 0.1484 | 0.3      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.001   | 76/0  | 83       | 0.9994  | 0.5419 | 1.1      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 2/5   | 4        | 1.0003  | 0.6451 | 0.3      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.001   | 37/0  | 44       | 1.0001  | 0.8292 | 0.9      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 2/5   | 4        | 1.0005  | 0.8591 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 5/0   | 6        | 0.1902  | 0.1428 | 0.1      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 2/1   | 4        | 0.2336  | 0.2349 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 0.001   | 14/0  | 9        | 0.9995  | 0.6653 | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 13/6  | 7        | 0.9983  | 0.6651 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 9/0   | 7        | 1.0001  | 0.8738 | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 2/5   | 3        | 0.9999  | 0.8640 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 4        | 0.1837  | 0.1511 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 4        | 0.1826  | 0.1511 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.001   | 9/0   | 7        | 0.9994  | 0.6679 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12/6  | 7        | 0.9969  | 0.6698 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 7/0   | 6        | 1.0001  | 0.8740 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 2/5   | 3        | 0.9999  | 0.8640 | 0.1      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data  | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:--------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 8  | no  | 0/0   | 7        | 0.0270  | 0.0127  | 0.2      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 4  | no  | 2/0   | 9        | 0.0301  | 0.0158  | 1.2      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| PDE-FIND         | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0230  | 0.0138  | 0.2      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | yes | 4/0   | 11       | 0.0366  | 0.0108  | 0.6      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.1      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0095  | 0.1183  | 0.2      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0085  | 0.1183  | 0.5      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0085  | 0.1183  | 0.1      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0109  | 0.1184  | 0.2      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 0.0099  | 0.1184  | 0.5      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0099  | 0.1184  | 0.1      | 🟢      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0559  | 0.2845  | 0.2      | 🟢      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | no  | 2/0   | 10       | 0.1000  | 0.2825  | 0.6      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0584  | 0.2845  | 0.1      | 🟢      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0548  | 0.2846  | 0.2      | 🟢      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | yes | 2/0   | 10       | 0.0992  | 0.2826  | 0.5      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0573  | 0.2846  | 0.1      | 🟢      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 7        | 0.0257  | 0.0138  | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 8  | no  | 0/0   | 7        | 76.5564 | 70.8534 | 0.5      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 7        | 0.0230  | 0.0138  | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 7        | 74.2567 | 64.4902 | 0.5      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0095  | 0.1183  | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 6/0   | 9        | 6.0366  | 63.1279 | 0.4      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0085  | 0.1183  | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0109  | 0.1184  | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | yes | 6/0   | 9        | 6.0511  | 63.2057 | 0.4      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0099  | 0.1184  | 0.1      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0559  | 0.2845  | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.5     | 12 | no  | 31/0  | 18       | 4.8985  | 58.9798 | 0.5      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0584  | 0.2845  | 0.1      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0548  | 0.2846  | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.5     | 12 | yes | 31/0  | 18       | 4.9068  | 59.0455 | 0.4      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0573  | 0.2846  | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 1.0     | 8  | no  | 0/0   | 7        | 0.0270  | 0.0127  | 0.5      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 4  | no  | 5/0   | 12       | 0.0302  | 0.0156  | 2.7      | 🔴      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 1.0     | 8  | yes | 0/0   | 7        | 0.1780  | 0.0134  | 0.5      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 8  | yes | 4/0   | 11       | 0.1818  | 0.0115  | 2.0      | 🔴      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0095  | 0.1183  | 0.4      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | no  | 8/0   | 15       | 0.0162  | 0.1173  | 1.3      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0085  | 0.1183  | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0109  | 0.1184  | 0.4      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | yes | 8/0   | 15       | 0.0150  | 0.1173  | 1.3      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0099  | 0.1184  | 0.1      | 🟢      |
| PolyJet          | 5%    | STRidge   | 1.0     | 12 | no  | 3/0   | 10       | 0.0895  | 0.2829  | 0.4      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | no  | 21/0  | 28       | 0.1237  | 0.2795  | 1.3      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0584  | 0.2845  | 0.1      | 🟢      |
| PolyJet          | 5%    | STRidge   | 1.0     | 12 | yes | 3/0   | 10       | 0.0887  | 0.2830  | 0.4      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | yes | 21/0  | 28       | 0.1229  | 0.2796  | 1.2      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0573  | 0.2846  | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 1.0     | 4  | no  | 0/0   | 4        | 0.0270  | 0.0170  | 0.1      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 8  | no  | 8/0   | 7        | 0.0277  | 0.0124  | 0.3      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0260  | 0.0138  | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 1.0     | 12 | yes | 3/0   | 5        | 0.0224  | 0.0153  | 0.1      | 🔴      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 8  | yes | 14/0  | 9        | 0.2600  | 0.0365  | 0.3      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0225  | 0.0153  | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0096  | 0.1183  | 0.1      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | no  | 10/0  | 8        | 0.0108  | 0.1176  | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0100  | 0.1183  | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0079  | 0.1185  | 0.1      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | yes | 15/0  | 10       | 0.0087  | 0.1174  | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0083  | 0.1185  | 0.1      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 4        | 0.0480  | 0.2848  | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 12 | no  | 8/0   | 7        | 0.0532  | 0.2830  | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0472  | 0.2848  | 0.1      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | yes | 0/0   | 4        | 0.0496  | 0.2848  | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 12 | yes | 11/0  | 8        | 0.0618  | 0.2826  | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0487  | 0.2847  | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 4  | no  | 20/0  | 11       | 0.0319  | 0.0151  | 0.3      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0260  | 0.0138  | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0224  | 0.0153  | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 8  | yes | 17/0  | 10       | 0.1887  | 0.0364  | 0.3      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0225  | 0.0153  | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0096  | 0.1183  | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 6/0   | 6        | 0.0109  | 0.1181  | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0100  | 0.1183  | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0079  | 0.1185  | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 6/0   | 6        | 0.0090  | 0.1182  | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0083  | 0.1185  | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 4        | 0.0480  | 0.2848  | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 2/0   | 5        | 0.0468  | 0.2843  | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0472  | 0.2848  | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0496  | 0.2848  | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 8  | yes | 5/0   | 6        | 0.7655  | 0.4798  | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0487  | 0.2847  | 0.1      | 🟢      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data  | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:--------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 9        | 0.0259  | 0.0135  | 0.2      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 8  | no  | 4/0   | 13       | 0.0361  | 0.0157  | 0.8      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| PDE-FIND         | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0225  | 0.0138  | 0.2      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 8  | yes | 2/0   | 11       | 0.1734  | 0.0114  | 0.8      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0227  | 0.0138  | 0.1      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 8  | no  | 0/0   | 7        | 0.2780  | 0.2910  | 0.2      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | no  | 3/0   | 11       | 0.0085  | 0.1169  | 0.6      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0098  | 0.1180  | 0.1      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0124  | 0.1184  | 0.2      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | yes | 3/0   | 11       | 0.0082  | 0.1173  | 0.6      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0115  | 0.1184  | 0.1      | 🟢      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 2/0   | 10       | 0.0814  | 0.2826  | 0.2      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | no  | 4/0   | 12       | 0.1113  | 0.2806  | 0.6      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0589  | 0.2837  | 0.1      | 🟢      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 2/0   | 10       | 0.0828  | 0.2835  | 0.2      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.5     | 12 | yes | 7/0   | 17       | 0.1194  | 0.2811  | 0.7      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0594  | 0.2846  | 0.1      | 🟢      |
| WSINDy           | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0257  | 0.0138  | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 8  | no  | 0/0   | 7        | 73.8555 | 70.8725 | 0.6      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 7        | 0.0225  | 0.0138  | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 7        | 71.4884 | 64.4902 | 0.5      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0227  | 0.0138  | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0107  | 0.1180  | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 6/0   | 9        | 11.0707 | 63.1179 | 0.4      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0098  | 0.1180  | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0124  | 0.1184  | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | yes | 6/0   | 9        | 11.0826 | 63.2065 | 0.4      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0115  | 0.1184  | 0.1      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0553  | 0.2837  | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.5     | 12 | no  | 31/0  | 18       | 9.0052  | 59.0284 | 0.4      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0589  | 0.2837  | 0.1      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0559  | 0.2846  | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.5     | 12 | yes | 31/0  | 18       | 8.9877  | 59.0564 | 0.4      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0594  | 0.2846  | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 1.0     | 4  | no  | 4/0   | 11       | 0.2187  | 0.0618  | 0.6      | 🔴      |
| PolyJet          | Clean | MSTLS     | 0.25    | 8  | no  | 4/0   | 11       | 0.1248  | 0.0147  | 2.1      | 🔴      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.2      | 🟢      |
| PolyJet          | Clean | STRidge   | 1.0     | 8  | yes | 0/0   | 7        | 0.1717  | 0.0134  | 0.5      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 8  | yes | 2/0   | 9        | 0.1729  | 0.0122  | 2.0      | 🔴      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0227  | 0.0138  | 0.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | no  | 10/0  | 17       | 0.0793  | 0.1169  | 0.4      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 12 | no  | 30/0  | 37       | 0.1504  | 0.1153  | 1.5      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0098  | 0.1180  | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | yes | 5/0   | 12       | 0.0838  | 0.1173  | 0.4      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 12 | yes | 30/0  | 37       | 0.1449  | 0.1156  | 1.6      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0115  | 0.1184  | 0.1      | 🟢      |
| PolyJet          | 5%    | STRidge   | 1.0     | 12 | no  | 15/0  | 22       | 0.1308  | 0.2791  | 0.4      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 12 | no  | 42/0  | 49       | 0.3284  | 0.2713  | 1.4      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0589  | 0.2837  | 0.1      | 🟢      |
| PolyJet          | 5%    | STRidge   | 1.0     | 12 | yes | 15/0  | 22       | 0.1333  | 0.2798  | 0.4      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 12 | yes | 42/0  | 49       | 0.3337  | 0.2718  | 1.5      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0594  | 0.2846  | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 1.0     | 12 | no  | 3/0   | 5        | 0.0258  | 0.0139  | 0.1      | 🔴      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 4  | no  | 18/0  | 10       | 0.2300  | 0.0550  | 0.3      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0260  | 0.0139  | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 1.0     | 12 | yes | 3/0   | 5        | 0.0231  | 0.0153  | 0.1      | 🔴      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 8  | yes | 16/0  | 10       | 0.2519  | 0.0364  | 0.2      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0232  | 0.0153  | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 12 | no  | 5/0   | 6        | 0.0146  | 0.1173  | 0.1      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | no  | 18/0  | 11       | 0.0964  | 0.1166  | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0110  | 0.1180  | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 12 | yes | 5/0   | 6        | 0.0121  | 0.1179  | 0.1      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | yes | 20/0  | 12       | 0.0969  | 0.1169  | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0084  | 0.1185  | 0.1      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 4        | 0.0456  | 0.2840  | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 12 | no  | 18/0  | 11       | 0.1749  | 0.2790  | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0450  | 0.2839  | 0.1      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | yes | 0/0   | 4        | 0.0497  | 0.2848  | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 12 | yes | 18/0  | 11       | 0.1786  | 0.2796  | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0489  | 0.2847  | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0259  | 0.0139  | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 4  | no  | 20/0  | 11       | 0.2312  | 0.0548  | 0.4      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0260  | 0.0139  | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0231  | 0.0153  | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 8  | yes | 17/0  | 10       | 0.1794  | 0.0363  | 0.3      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0232  | 0.0153  | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | no  | 0/0   | 4        | 0.0108  | 0.1180  | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.5     | 12 | no  | 20/0  | 11       | 0.0967  | 0.1166  | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0110  | 0.1180  | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | yes | 0/0   | 4        | 0.0080  | 0.1185  | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.5     | 12 | yes | 20/0  | 11       | 0.0891  | 0.1169  | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0084  | 0.1185  | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 8  | no  | 2/0   | 5        | 0.7172  | 0.4683  | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.5     | 12 | no  | 20/0  | 12       | 0.1745  | 0.2790  | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0450  | 0.2839  | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0497  | 0.2848  | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.5     | 12 | yes | 20/0  | 11       | 0.2054  | 0.2798  | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0489  | 0.2847  | 0.1      | 🟢      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ) |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:--------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0314       | 10372.4 |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0291       | 15893.5 |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0288       | 28641.5 |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.1275       | 10430.3 |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0316       | 15952.5 |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0290       | 28643.3 |
| Clean | weak (IBP)    | spectral   | 4  | yes | 33.2439      | 10380.5 |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.1679       | 15897.9 |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0268       | 28695.1 |
| Clean | weak (IBP)    | polynomial | 4  | yes | 33.2439      | 10448.2 |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.1679       | 15962.4 |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0268       | 28691.7 |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 11.8754      | 692.9   |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 0.3767       | 1432.1  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 0.1208       | 3350.2  |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 4.1394       | 8016.7  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.3327       | 13410.6 |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.1204       | 25224.6 |
| 2%    | weak (IBP)    | spectral   | 4  | yes | 35.0283      | 691.0   |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 0.4121       | 1432.2  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 0.1202       | 3350.1  |
| 2%    | weak (IBP)    | polynomial | 4  | yes | 35.0283      | 8011.6  |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 0.4121       | 13414.0 |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.1202       | 25249.1 |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 28.1988      | 744.2   |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 0.9184       | 1578.6  |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 0.2874       | 4258.3  |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 9.8233       | 3051.7  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 0.8096       | 6492.0  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 0.2864       | 13259.2 |
| 5%    | weak (IBP)    | spectral   | 4  | yes | 42.2120      | 741.0   |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 0.9352       | 1579.0  |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 0.2869       | 4258.1  |
| 5%    | weak (IBP)    | polynomial | 4  | yes | 42.2120      | 3037.9  |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 0.9352       | 6488.1  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 0.2869       | 13263.5 |
---
| Source                              | Wall time | CPU                                           | Threads |
|:------------------------------------|----------:|:----------------------------------------------|:-------:|
| `sweep_20260715_150518_sh_full.csv` | 40.5 min  | AMD Ryzen 7 PRO 7840U w/ Radeon 780M Graphics | 1       |
