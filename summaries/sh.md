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
| Clean | pointwise     | spectral   | 🟩        | 0/0 ✓             | 0.0273     | mod_multigraph/RFE     | 0.0209    | WSINDy/RFE       | 0/0 ✓     | 0.0274      | 0.0208     | 0.0338       |
| Clean | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.2077     | mod_multigraph/RFE     | 0.1511    | WSINDy/RFE       | 0/0 ✓     | 0.2110      | 0.1484     | 0.2347       |
| Clean | weak (no IBP) | spectral   | 🟥        | 0/0 ✓             | 0.0259     | mod_multigraph/STRidge | 0.0138    | WSINDy/STRidge   | 0/0 ✓     | 0.0257      | 0.0138     | 0.0288       |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0259     | mod_multigraph/STRidge | 0.0139    | WSINDy/STRidge   | 0/0 ✓     | 0.0258      | 0.0138     | 0.0290       |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0224     | mod_multigraph/STRidge | 0.0153    | PDE-FIND/STRidge | 0/0 ✓     | 0.0230      | 0.0138     | 0.0268       |
| Clean | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0224     | multigraph/STRidge     | 0.0153    | PDE-FIND/STRidge | 0/0 ✓     | 0.0230      | 0.0138     | 0.0268       |
| 2%    | pointwise     | spectral   | —        | 20/0 ✗            | 0.5341     | —                      | 0.6620    | —                | 55/0 ✗    | 0.4476      | 0.6093     | 1150.0123    |
| 2%    | pointwise     | polynomial | —        | 16/0 ✗            | 0.4413     | —                      | 0.6076    | —                | 18/0 ✗    | 0.8080      | 0.8462     | 23.6646      |
| 2%    | weak (no IBP) | spectral   | 🟥        | 0/0 ✓             | 0.0299     | multigraph/STRidge     | 0.1183    | PolyJet/STRidge  | 0/0 ✓     | 0.0298      | 0.1183     | 0.1208       |
| 2%    | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0289     | multigraph/STRidge     | 0.1180    | WSINDy/STRidge   | 0/0 ✓     | 0.0288      | 0.1180     | 0.1204       |
| 2%    | weak (IBP)    | spectral   | 🟥        | 0/0 ✓             | 0.0278     | multigraph/STRidge     | 0.1185    | PolyJet/STRidge  | 0/0 ✓     | 0.0237      | 0.1184     | 0.1202       |
| 2%    | weak (IBP)    | polynomial | 🟥        | 0/0 ✓             | 0.0278     | multigraph/STRidge     | 0.1185    | PDE-FIND/STRidge | 0/0 ✓     | 0.0237      | 0.1184     | 0.1202       |
| 5%    | pointwise     | spectral   | —        | 14/0 ✗            | 0.5627     | —                      | 0.8619    | —                | 28/0 ✗    | 0.5768      | 0.8683     | 2032.3997    |
| 5%    | pointwise     | polynomial | —        | 20/0 ✗            | 0.4641     | —                      | 0.8457    | —                | 41/0 ✗    | 0.8135      | 0.9274     | 41.8218      |
| 5%    | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0923     | mod_multigraph/STRidge | 0.2848    | PDE-FIND/STRidge | 0/0 ✓     | 0.0956      | 0.2845     | 0.2874       |
| 5%    | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0885     | multigraph/STRidge     | 0.2840    | WSINDy/STRidge   | 0/0 ✓     | 0.0918      | 0.2837     | 0.2864       |
| 5%    | weak (IBP)    | spectral   | 🟥        | 0/0 ✓             | 0.0922     | multigraph/STRidge     | 0.2848    | PDE-FIND/STRidge | 0/0 ✓     | 0.0891      | 0.2846     | 0.2869       |
| 5%    | weak (IBP)    | polynomial | 🟥        | 0/0 ✓             | 0.0922     | multigraph/STRidge     | 0.2848    | WSINDy/STRidge   | 0/0 ✓     | 0.0891      | 0.2846     | 0.2869       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                   | fp/fn  | E_coeff | Time (s) | vs. best   |
|:-------------------------------|:-----------------------------|:-------|--------:|---------:|-----------:|
| Multigraphs + STRidge          | polynomial, weak (IBP), m=12 | 0/0 ✓  | 0.0224  | 0.2      | —          |
| Modified Multigraphs + STRidge | spectral, weak (IBP), m=12   | 0/0 ✓  | 0.0224  | 0.3      | +0.0%      |
| Multigraphs + RFE              | spectral, weak (IBP), m=12   | 0/0 ✓  | 0.0225  | 0.1      | +0.3%      |
| Modified Multigraphs + RFE     | spectral, weak (IBP), m=12   | 0/0 ✓  | 0.0225  | 0.1      | +0.3%      |
| PDE-FIND + STRidge             | polynomial, weak (IBP), m=12 | 0/0 ✓  | 0.0230  | 0.5      | +2.9%      |
| WSINDy + STRidge               | polynomial, weak (IBP), m=12 | 0/0 ✓  | 0.0230  | 0.4      | +2.9%      |
| PDE-FIND + RFE                 | spectral, weak (IBP), m=12   | 0/0 ✓  | 0.0233  | 0.2      | +4.1%      |
| PolyJet + RFE                  | spectral, weak (IBP), m=12   | 0/0 ✓  | 0.0233  | 0.2      | +4.1%      |
| WSINDy + RFE                   | spectral, weak (IBP), m=12   | 0/0 ✓  | 0.0233  | 0.1      | +4.1%      |
| PolyJet + STRidge              | spectral, weak (no IBP), m=4 | 0/0 ✓  | 0.0271  | 3.4      | +21.0%     |
| WSINDy + MSTLS                 | polynomial, weak (IBP), m=8  | 0/0 ✓  | 74.2567 | 0.8      | +331576.1% |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=4 | 2/0 ✗  | 0.0301  | 3.4      | ✗          |
| PolyJet + MSTLS                | polynomial, weak (IBP), m=8  | 2/0 ✗  | 0.1779  | 3.8      | ✗          |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=8 | 8/0 ✗  | 0.0277  | 0.4      | ✗          |
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=8  | 17/0 ✗ | 0.1848  | 0.4      | ✗          |

### 2% noise

| Method                         | Best conf.                    | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:------------------------------|:-------|--------:|---------:|---------:|
| PolyJet + STRidge              | spectral, weak (IBP), m=12    | 0/0 ✓  | 0.0237  | 1.3      | —        |
| PDE-FIND + STRidge             | spectral, weak (IBP), m=12    | 0/0 ✓  | 0.0237  | 0.5      | +0.0%    |
| WSINDy + STRidge               | spectral, weak (IBP), m=12    | 0/0 ✓  | 0.0237  | 0.3      | +0.0%    |
| WSINDy + RFE                   | spectral, weak (IBP), m=12    | 0/0 ✓  | 0.0244  | 0.1      | +3.1%    |
| PolyJet + RFE                  | spectral, weak (IBP), m=12    | 0/0 ✓  | 0.0244  | 0.2      | +3.1%    |
| PDE-FIND + MSTLS               | spectral, weak (IBP), m=12    | 0/0 ✓  | 0.0244  | 1.1      | +3.1%    |
| PDE-FIND + RFE                 | polynomial, weak (IBP), m=12  | 0/0 ✓  | 0.0244  | 0.1      | +3.1%    |
| Multigraphs + STRidge          | polynomial, weak (IBP), m=12  | 0/0 ✓  | 0.0278  | 0.2      | +17.5%   |
| Modified Multigraphs + STRidge | spectral, weak (IBP), m=12    | 0/0 ✓  | 0.0278  | 0.2      | +17.5%   |
| Multigraphs + RFE              | spectral, weak (IBP), m=12    | 0/0 ✓  | 0.0287  | 0.1      | +21.3%   |
| Modified Multigraphs + RFE     | spectral, weak (IBP), m=12    | 0/0 ✓  | 0.0287  | 0.1      | +21.3%   |
| Modified Multigraphs + MSTLS   | spectral, weak (no IBP), m=12 | 6/0 ✗  | 0.1479  | 0.3      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (IBP), m=12    | 6/0 ✗  | 73.6692 | 0.6      | ✗        |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=12 | 8/0 ✗  | 0.3180  | 4.2      | ✗        |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=12 | 10/0 ✗ | 0.2035  | 0.2      | ✗        |

### 5% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0885  | 0.2      | —        |
| Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0885  | 0.3      | +0.0%    |
| WSINDy + STRidge               | polynomial, weak (IBP), m=12    | 0/0 ✓  | 0.0891  | 0.4      | +0.7%    |
| PDE-FIND + STRidge             | spectral, weak (IBP), m=12      | 0/0 ✓  | 0.0891  | 0.8      | +0.7%    |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0904  | 0.1      | +2.2%    |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0904  | 0.1      | +2.2%    |
| WSINDy + RFE                   | polynomial, weak (IBP), m=12    | 0/0 ✓  | 0.0910  | 0.1      | +2.8%    |
| PDE-FIND + RFE                 | polynomial, weak (IBP), m=12    | 0/0 ✓  | 0.0910  | 0.1      | +2.8%    |
| PolyJet + RFE                  | spectral, weak (IBP), m=12      | 0/0 ✓  | 0.0910  | 0.1      | +2.8%    |
| Modified Multigraphs + MSTLS   | spectral, weak (no IBP), m=12   | 2/0 ✗  | 0.4611  | 0.3      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (IBP), m=12      | 3/0 ✗  | 0.1693  | 1.9      | ✗        |
| PolyJet + STRidge              | spectral, weak (IBP), m=12      | 5/0 ✗  | 0.1868  | 1.0      | ✗        |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=12   | 8/0 ✗  | 0.4818  | 0.3      | ✗        |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=12   | 21/0 ✗ | 0.6253  | 2.7      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (IBP), m=12      | 31/0 ✗ | 89.9817 | 1.2      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 0/0   | 7        | 0.0276  | 0.0208 | 1.5      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 0/0   | 7        | 0.0274  | 0.0208 | 0.2      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 0.001   | 5/1   | 11       | 0.8356  | 0.8680 | 1.0      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 0/5   | 2        | 0.8437  | 0.8727 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.001   | 10/1  | 16       | 0.8260  | 0.9354 | 2.9      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 0/5   | 2        | 0.8436  | 0.9404 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 0/0   | 7        | 0.0276  | 0.0208 | 1.9      | 🟢      |
| WSINDy           | Clean | RFE       |         | 0/0   | 7        | 0.0274  | 0.0208 | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.001   | 79/0  | 38       | 0.5160  | 0.6568 | 2.2      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 2/5   | 3        | 0.6324  | 0.7532 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.001   | 28/0  | 14       | 0.5768  | 0.8683 | 1.8      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 6/5   | 4        | 0.6164  | 0.8806 | 0.2      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.01    | 0/0   | 7        | 0.0276  | 0.0208 | 3.4      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0274  | 0.0208 | 0.3      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.001   | 55/0  | 62       | 0.4476  | 0.6093 | 6.8      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/5   | 3        | 0.5788  | 0.6951 | 0.3      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.001   | 29/1  | 35       | 0.4965  | 0.8491 | 3.6      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/5   | 3        | 0.6099  | 0.8872 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.01    | 5/0   | 6        | 0.0278  | 0.0206 | 0.4      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 4        | 0.0273  | 0.0209 | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.01    | 20/0  | 12       | 0.5341  | 0.6620 | 0.4      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 2/5   | 3        | 0.5616  | 0.6850 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 14/0  | 9        | 0.5627  | 0.8619 | 0.6      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 2/5   | 3        | 0.6034  | 0.8860 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.01    | 7/0   | 7        | 0.0280  | 0.0203 | 0.7      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 4        | 0.0273  | 0.0209 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.001   | 21/0  | 13       | 0.5302  | 0.6592 | 0.5      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 2/5   | 3        | 0.5616  | 0.6850 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 20/0  | 10       | 0.5598  | 0.8620 | 0.5      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 2/5   | 3        | 0.6034  | 0.8860 | 0.1      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 0/0   | 7        | 0.2122  | 0.1484 | 0.9      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 0/0   | 7        | 0.2110  | 0.1484 | 0.1      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 0.001   | 18/0  | 25       | 0.8080  | 0.8462 | 1.6      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 1/6   | 2        | 0.8318  | 0.8649 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.001   | 41/0  | 48       | 0.8135  | 0.9274 | 2.2      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 1/6   | 2        | 0.8322  | 0.9391 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.01    | 0/0   | 7        | 0.2122  | 0.1484 | 0.6      | 🟢      |
| WSINDy           | Clean | RFE       |         | 0/0   | 7        | 0.2110  | 0.1484 | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.001   | 41/0  | 21       | 0.4667  | 0.6089 | 1.0      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 2/5   | 3        | 0.6227  | 0.7105 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.001   | 77/0  | 36       | 0.5128  | 0.8500 | 1.6      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 2/5   | 3        | 0.6230  | 0.8799 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.01    | 19/0  | 26       | 0.2479  | 0.0206 | 2.0      | 🔴      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.2110  | 0.1484 | 0.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.001   | 44/1  | 50       | 0.3412  | 0.5483 | 3.8      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 2/5   | 4        | 0.5127  | 0.6451 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.001   | 43/0  | 50       | 0.3772  | 0.8225 | 4.4      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 2/5   | 4        | 0.5204  | 0.8591 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.01    | 5/0   | 6        | 0.2168  | 0.1428 | 0.4      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 2/1   | 4        | 0.2472  | 0.2349 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 0.001   | 16/0  | 11       | 0.4413  | 0.6076 | 0.4      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 13/6  | 7        | 0.5445  | 0.6651 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 20/0  | 12       | 0.4641  | 0.8457 | 0.4      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 2/5   | 3        | 0.5626  | 0.8640 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.01    | 14/0  | 8        | 0.2177  | 0.0941 | 0.5      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 4        | 0.2077  | 0.1511 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.001   | 30/0  | 18       | 0.4313  | 0.6006 | 0.7      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12/6  | 7        | 0.5380  | 0.6698 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 20/0  | 12       | 0.4641  | 0.8457 | 0.8      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 2/5   | 3        | 0.5626  | 0.8640 | 0.1      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data  | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:--------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 8  | no  | 0/0   | 7        | 0.0270  | 0.0127  | 1.1      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 4  | no  | 2/0   | 9        | 0.0301  | 0.0158  | 3.4      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.2      | 🟢      |
| PDE-FIND         | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0230  | 0.0138  | 0.7      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | yes | 4/0   | 11       | 0.0366  | 0.0108  | 1.7      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.2      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0298  | 0.1183  | 0.5      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0304  | 0.1183  | 1.1      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0304  | 0.1183  | 0.1      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0237  | 0.1184  | 0.5      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 0.0244  | 0.1184  | 1.1      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0244  | 0.1184  | 0.1      | 🟢      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0956  | 0.2845  | 0.8      | 🟢      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | no  | 3/0   | 10       | 0.1730  | 0.2825  | 1.9      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0974  | 0.2845  | 0.1      | 🟢      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0891  | 0.2846  | 0.8      | 🟢      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | yes | 3/0   | 10       | 0.1693  | 0.2826  | 1.9      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0910  | 0.2846  | 0.1      | 🟢      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 7        | 0.0257  | 0.0138  | 0.6      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 8  | no  | 0/0   | 7        | 76.5564 | 70.8534 | 1.4      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 7        | 0.0230  | 0.0138  | 0.6      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 7        | 74.2567 | 64.4902 | 1.4      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 7        | 0.0298  | 0.1183  | 0.3      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 6/0   | 9        | 73.9646 | 63.1279 | 0.6      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0304  | 0.1183  | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 7        | 0.0237  | 0.1184  | 0.3      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | yes | 6/0   | 9        | 73.6692 | 63.2057 | 0.6      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0244  | 0.1184  | 0.1      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0956  | 0.2845  | 0.6      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.5     | 12 | no  | 31/0  | 18       | 90.3314 | 58.9798 | 1.3      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0974  | 0.2845  | 0.2      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0891  | 0.2846  | 0.6      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.5     | 12 | yes | 31/0  | 18       | 89.9817 | 59.0455 | 1.2      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0910  | 0.2846  | 0.2      | 🟢      |
| PolyJet          | Clean | STRidge   | 0.1     | 4  | no  | 0/0   | 7        | 0.0271  | 0.0169  | 3.4      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 4  | no  | 5/0   | 12       | 0.0302  | 0.0156  | 10.4     | 🔴      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.2      | 🟢      |
| PolyJet          | Clean | STRidge   | 1.0     | 12 | yes | 5/0   | 12       | 0.1607  | 0.0118  | 1.7      | 🔴      |
| PolyJet          | Clean | MSTLS     | 0.25    | 8  | yes | 4/0   | 11       | 0.1818  | 0.0115  | 7.1      | 🔴      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0298  | 0.1183  | 1.1      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | no  | 8/0   | 15       | 0.3180  | 0.1173  | 4.2      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0304  | 0.1183  | 0.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0237  | 0.1184  | 1.3      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | yes | 8/0   | 15       | 0.3257  | 0.1173  | 4.5      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0244  | 0.1184  | 0.2      | 🟢      |
| PolyJet          | 5%    | STRidge   | 1.0     | 12 | no  | 5/0   | 12       | 0.1905  | 0.2829  | 1.0      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | no  | 21/0  | 28       | 0.6253  | 0.2795  | 2.7      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0974  | 0.2845  | 0.1      | 🟢      |
| PolyJet          | 5%    | STRidge   | 1.0     | 12 | yes | 5/0   | 12       | 0.1868  | 0.2830  | 1.0      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | yes | 21/0  | 28       | 0.6326  | 0.2796  | 2.6      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0910  | 0.2846  | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 1.0     | 8  | no  | 0/0   | 4        | 0.0271  | 0.0129  | 0.3      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 8  | no  | 8/0   | 7        | 0.0277  | 0.0124  | 0.4      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0260  | 0.0138  | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0224  | 0.0153  | 0.2      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 8  | yes | 14/0  | 9        | 0.2600  | 0.0365  | 0.4      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0225  | 0.0153  | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0299  | 0.1183  | 0.2      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | no  | 10/0  | 8        | 0.2035  | 0.1176  | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0306  | 0.1183  | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0278  | 0.1185  | 0.2      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | yes | 15/0  | 10       | 0.0389  | 0.1174  | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0287  | 0.1185  | 0.1      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 4        | 0.0923  | 0.2848  | 0.2      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 12 | no  | 8/0   | 7        | 0.4818  | 0.2830  | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0941  | 0.2848  | 0.1      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | yes | 0/0   | 4        | 0.0922  | 0.2848  | 0.2      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 12 | yes | 11/0  | 8        | 0.4376  | 0.2826  | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0940  | 0.2847  | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0259  | 0.0138  | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 4  | no  | 20/0  | 11       | 0.0319  | 0.0151  | 0.8      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0260  | 0.0138  | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0224  | 0.0153  | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 8  | yes | 17/0  | 10       | 0.1887  | 0.0364  | 0.6      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0225  | 0.0153  | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0299  | 0.1183  | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 6/0   | 6        | 0.1479  | 0.1181  | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0306  | 0.1183  | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0278  | 0.1185  | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 6/0   | 6        | 0.2112  | 0.1182  | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0287  | 0.1185  | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0923  | 0.2848  | 0.2      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 2/0   | 5        | 0.4611  | 0.2843  | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0941  | 0.2848  | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0922  | 0.2848  | 0.2      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 8  | yes | 5/0   | 6        | 1.8727  | 0.4798  | 0.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0940  | 0.2847  | 0.1      | 🟢      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data  | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:--------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | no  | 2/0   | 9        | 0.0260  | 0.0135  | 0.5      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | no  | 4/0   | 11       | 0.0386  | 0.0110  | 1.2      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0260  | 0.0138  | 0.1      | 🟢      |
| PDE-FIND         | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0230  | 0.0138  | 0.5      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | yes | 4/0   | 11       | 0.0365  | 0.0109  | 1.1      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.1      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | no  | 2/0   | 9        | 0.0373  | 0.1178  | 0.5      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | no  | 4/0   | 11       | 0.0539  | 0.1169  | 0.9      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0295  | 0.1180  | 0.1      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0237  | 0.1184  | 0.4      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | yes | 4/0   | 11       | 0.0513  | 0.1173  | 1.0      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0244  | 0.1184  | 0.1      | 🟢      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 3/0   | 10       | 0.1416  | 0.2826  | 0.6      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | no  | 5/0   | 12       | 0.1892  | 0.2806  | 1.3      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0937  | 0.2837  | 0.1      | 🟢      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 3/0   | 10       | 0.1413  | 0.2835  | 0.6      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.5     | 12 | yes | 10/0  | 17       | 0.2077  | 0.2811  | 1.4      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0910  | 0.2846  | 0.1      | 🟢      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 7        | 0.0258  | 0.0138  | 0.4      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 8  | no  | 0/0   | 7        | 76.6159 | 70.8725 | 0.8      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0260  | 0.0138  | 0.1      | 🟢      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 7        | 0.0230  | 0.0138  | 0.4      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 7        | 74.2567 | 64.4902 | 0.8      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 7        | 0.0288  | 0.1180  | 0.3      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 6/0   | 9        | 74.1547 | 63.1179 | 0.6      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0295  | 0.1180  | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 7        | 0.0237  | 0.1184  | 0.3      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | yes | 6/0   | 9        | 73.7051 | 63.2065 | 0.6      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0244  | 0.1184  | 0.1      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 7        | 0.0918  | 0.2837  | 0.4      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.5     | 12 | no  | 31/0  | 18       | 90.9074 | 59.0284 | 0.8      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0937  | 0.2837  | 0.1      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 7        | 0.0891  | 0.2846  | 0.4      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.5     | 12 | yes | 31/0  | 18       | 90.7376 | 59.0564 | 0.8      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0910  | 0.2846  | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 1.0     | 4  | no  | 4/0   | 11       | 0.2188  | 0.0618  | 1.5      | 🔴      |
| PolyJet          | Clean | MSTLS     | 0.25    | 8  | no  | 4/0   | 11       | 0.1316  | 0.0147  | 4.0      | 🔴      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0260  | 0.0138  | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 1.0     | 12 | yes | 5/0   | 12       | 0.1617  | 0.0118  | 1.0      | 🔴      |
| PolyJet          | Clean | MSTLS     | 0.25    | 8  | yes | 2/0   | 9        | 0.1779  | 0.0122  | 3.8      | 🔴      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | no  | 10/0  | 17       | 2.1529  | 0.1169  | 1.1      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 8  | no  | 30/0  | 37       | 1.1538  | 0.2199  | 4.1      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0295  | 0.1180  | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | yes | 5/0   | 12       | 2.2502  | 0.1173  | 1.0      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 12 | yes | 30/0  | 37       | 3.3853  | 0.1156  | 3.3      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0244  | 0.1184  | 0.1      | 🟢      |
| PolyJet          | 5%    | STRidge   | 1.0     | 12 | no  | 15/0  | 22       | 2.7400  | 0.2791  | 1.2      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 12 | no  | 42/0  | 49       | 4.3316  | 0.2713  | 3.1      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | no  | 0/0   | 7        | 0.0937  | 0.2837  | 0.1      | 🟢      |
| PolyJet          | 5%    | STRidge   | 1.0     | 12 | yes | 15/0  | 22       | 2.7099  | 0.2798  | 1.2      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 12 | yes | 42/0  | 49       | 4.2928  | 0.2718  | 3.4      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | yes | 0/0   | 7        | 0.0910  | 0.2846  | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 1.0     | 8  | no  | 0/0   | 4        | 0.0279  | 0.0178  | 0.3      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 4  | no  | 18/0  | 10       | 0.2357  | 0.0550  | 0.5      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0261  | 0.0139  | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0224  | 0.0153  | 0.2      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 8  | yes | 14/0  | 9        | 0.2598  | 0.0365  | 0.4      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0225  | 0.0153  | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0289  | 0.1180  | 0.2      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | no  | 20/0  | 12       | 2.5713  | 0.1165  | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0296  | 0.1180  | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0278  | 0.1185  | 0.2      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | yes | 20/0  | 12       | 2.4854  | 0.1168  | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0287  | 0.1185  | 0.1      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0885  | 0.2840  | 0.2      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 12 | no  | 17/0  | 11       | 5.1180  | 0.2796  | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0904  | 0.2839  | 0.1      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0922  | 0.2848  | 0.2      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 12 | yes | 17/0  | 11       | 5.0452  | 0.2802  | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0940  | 0.2847  | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0259  | 0.0139  | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 4  | no  | 20/0  | 11       | 0.2366  | 0.0548  | 0.5      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0261  | 0.0139  | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0224  | 0.0153  | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 8  | yes | 17/0  | 10       | 0.1848  | 0.0363  | 0.4      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0225  | 0.0153  | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0289  | 0.1180  | 0.3      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.5     | 12 | no  | 20/0  | 10       | 2.6473  | 0.1165  | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0296  | 0.1180  | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0278  | 0.1185  | 0.3      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.5     | 12 | yes | 20/0  | 11       | 2.3250  | 0.1169  | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0287  | 0.1185  | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0885  | 0.2840  | 0.3      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.5     | 8  | no  | 12/0  | 9        | 2.1626  | 0.4565  | 0.6      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0904  | 0.2839  | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0922  | 0.2848  | 0.4      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.5     | 12 | yes | 20/0  | 11       | 9.7093  | 0.2798  | 0.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0940  | 0.2847  | 0.3      | 🟢      |


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
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.1679       | 15897.9 |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0268       | 28695.1 |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.1679       | 15976.1 |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0268       | 28279.7 |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 11.8754      | 692.9   |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 0.3767       | 1432.1  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 0.1208       | 3350.2  |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 4.1394       | 8018.2  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.3327       | 13590.1 |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.1204       | 23171.2 |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 0.4121       | 1432.2  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 0.1202       | 3350.1  |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 0.4121       | 13595.3 |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.1202       | 23187.1 |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 28.1988      | 744.2   |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 0.9184       | 1578.6  |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 0.2874       | 4258.3  |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 9.8233       | 3075.9  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 0.8096       | 6209.9  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 0.2864       | 11568.1 |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 0.9352       | 1579.0  |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 0.2869       | 4258.1  |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 0.9352       | 6209.1  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 0.2869       | 11573.6 |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation sh`.

| Source                                  | Wall time | CPU | Threads |
|:----------------------------------------|----------:|:----|:-------:|
| `sweep_paper_v2_sh_20260804_223951.csv` | ? min     | ?   | ?       |
