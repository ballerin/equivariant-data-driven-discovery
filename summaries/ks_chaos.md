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
| Clean | pointwise     | spectral   | 🟩        | 0/0 ✓             | 0.0055     | multigraph/RFE     | 0.0197    | PolyJet/RFE     | 0/0 ✓     | 0.0056      | 0.0197     | 0.0205       |
| Clean | pointwise     | polynomial | —        | 10/0 ✗            | 0.2774     | —                  | 0.1129    | —               | 2/0 ✗     | 0.2174      | 0.2119     | 0.3112       |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0017     | multigraph/STRidge | 0.0095    | PolyJet/RFE     | 0/0 ✓     | 0.0019      | 0.0095     | 0.0066       |
| Clean | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.0350     | multigraph/MSTLS   | 0.0071    | PolyJet/STRidge | 0/0 ✓     | 0.0335      | 0.0071     | 0.0115       |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0020     | multigraph/MSTLS   | 0.0087    | PolyJet/MSTLS   | 0/0 ✓     | 0.0024      | 0.0087     | 0.0089       |
| Clean | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0277     | multigraph/STRidge | 0.0156    | PolyJet/STRidge | 0/0 ✓     | 0.0287      | 0.0146     | 0.0150       |
| 2%    | pointwise     | spectral   | —        | 17/0 ✗            | 0.9997     | —                  | 0.6535    | —               | 65/0 ✗    | 0.9997      | 0.6516     | 85.8534      |
| 2%    | pointwise     | polynomial | —        | 33/0 ✗            | 0.8591     | —                  | 0.6178    | —               | 62/0 ✗    | 0.8601      | 0.6126     | 1.8548       |
| 2%    | weak (no IBP) | spectral   | —        | 5/0 ✗             | 0.3943     | —                  | 0.1551    | —               | 3/0 ✗     | 0.4115      | 0.2057     | 3.3691       |
| 2%    | weak (no IBP) | polynomial | —        | 5/0 ✗             | 0.2663     | —                  | 0.1483    | —               | 8/0 ✗     | 0.7204      | 0.2721     | 0.2790       |
| 2%    | weak (IBP)    | spectral   | —        | 4/0 ✗             | 0.4568     | —                  | 0.2039    | —               | 4/0 ✗     | 0.4053      | 0.1904     | 3.3618       |
| 2%    | weak (IBP)    | polynomial | —        | 5/0 ✗             | 0.2721     | —                  | 0.1395    | —               | 9/0 ✗     | 0.3918      | 0.1359     | 0.2983       |
| 5%    | pointwise     | spectral   | —        | 18/0 ✗            | 0.9998     | —                  | 0.8839    | —               | 87/0 ✗    | 0.9998      | 0.8824     | 133.9474     |
| 5%    | pointwise     | polynomial | —        | 37/0 ✗            | 0.9432     | —                  | 0.8801    | —               | 84/0 ✗    | 0.9437      | 0.8784     | 3.0277       |
| 5%    | weak (no IBP) | spectral   | —        | 30/0 ✗            | 0.9985     | —                  | 0.5367    | —               | 36/0 ✗    | 0.9689      | 0.3587     | 19.5646      |
| 5%    | weak (no IBP) | polynomial | —        | 32/0 ✗            | 0.8969     | —                  | 0.5044    | —               | 69/0 ✗    | 0.7383      | 0.3011     | 1.6411       |
| 5%    | weak (IBP)    | spectral   | —        | 32/0 ✗            | 0.9953     | —                  | 0.4681    | —               | 38/0 ✗    | 0.9862      | 0.3397     | 19.9290      |
| 5%    | weak (IBP)    | polynomial | —        | 37/0 ✗            | 0.6972     | —                  | 0.2916    | —               | 64/0 ✗    | 0.7462      | 0.2840     | 1.6435       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                | Best conf.                     | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:-------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + STRidge | spectral, weak (no IBP), m=8   | 0/0 ✓  | 0.0017  | 1.7      | —        |
| Multigraphs + MSTLS   | spectral, weak (no IBP), m=8   | 0/0 ✓  | 0.0017  | 0.8      | +1.0%    |
| PolyJet + RFE         | spectral, weak (no IBP), m=8   | 0/0 ✓  | 0.0019  | 0.8      | +14.5%   |
| PolyJet + MSTLS       | spectral, weak (no IBP), m=8   | 0/0 ✓  | 0.0019  | 39.4     | +14.5%   |
| PolyJet + STRidge     | spectral, weak (no IBP), m=8   | 0/0 ✓  | 0.0021  | 11.5     | +23.6%   |
| Multigraphs + RFE     | spectral, weak (IBP), m=8      | 0/0 ✓  | 0.0028  | 0.1      | +66.2%   |
| WSINDy + STRidge      | spectral, weak (IBP), m=12     | 71/0 ✗ | 0.4838  | 4.6      | ✗        |
| PDE-FIND + STRidge    | polynomial, weak (no IBP), m=8 | 12/2 ✗ | 4.0849  | 5.7      | ✗        |
| PDE-FIND + RFE        | spectral, weak (no IBP), m=12  | 1/7 ✗  | 0.4515  | 0.4      | ✗        |
| WSINDy + RFE          | spectral, weak (no IBP), m=12  | 1/7 ✗  | 0.4515  | 0.3      | ✗        |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12  | 1/7 ✗  | 0.4515  | 4.1      | ✗        |
| WSINDy + MSTLS        | spectral, weak (no IBP), m=12  | 1/7 ✗  | 0.7177  | 2.4      | ✗        |

### 2% noise

| Method                | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| PolyJet + MSTLS       | spectral, weak (no IBP), m=8    | 3/0 ✗  | 0.4115  | 39.3     | ✗        |
| Multigraphs + MSTLS   | spectral, weak (IBP), m=12      | 4/0 ✗  | 0.4568  | 0.7      | ✗        |
| Multigraphs + STRidge | polynomial, weak (no IBP), m=4  | 9/0 ✗  | 0.6901  | 1.6      | ✗        |
| PolyJet + STRidge     | polynomial, weak (no IBP), m=8  | 13/0 ✗ | 0.3970  | 13.2     | ✗        |
| WSINDy + STRidge      | polynomial, weak (no IBP), m=4  | 64/0 ✗ | 0.8885  | 5.2      | ✗        |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=4    | 10/2 ✗ | 0.9752  | 5.4      | ✗        |
| Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/5 ✗  | 0.8261  | 0.1      | ✗        |
| PolyJet + RFE         | polynomial, weak (IBP), m=8     | 0/6 ✗  | 0.8870  | 0.7      | ✗        |
| PDE-FIND + MSTLS      | polynomial, weak (IBP), m=12    | 1/7 ✗  | 0.8908  | 4.3      | ✗        |
| PDE-FIND + RFE        | polynomial, weak (IBP), m=12    | 1/7 ✗  | 0.8908  | 0.5      | ✗        |
| WSINDy + RFE          | polynomial, weak (IBP), m=12    | 1/7 ✗  | 0.8908  | 0.3      | ✗        |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=12 | 1/7 ✗  | 0.9439  | 2.7      | ✗        |

### 5% noise

| Method                | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + STRidge | spectral, pointwise             | 18/0 ✗ | 0.9998  | 1.7      | ✗        |
| PolyJet + STRidge     | spectral, weak (no IBP), m=8    | 36/0 ✗ | 0.9689  | 15.2     | ✗        |
| WSINDy + STRidge      | polynomial, weak (no IBP), m=12 | 83/0 ✗ | 1.2103  | 4.3      | ✗        |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=4  | 8/2 ✗  | 0.9190  | 38.7     | ✗        |
| PDE-FIND + STRidge    | polynomial, weak (no IBP), m=12 | 10/2 ✗ | 1.4514  | 5.4      | ✗        |
| Multigraphs + MSTLS   | polynomial, weak (no IBP), m=8  | 9/3 ✗  | 0.9463  | 1.1      | ✗        |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.9016  | 0.1      | ✗        |
| PDE-FIND + RFE        | polynomial, pointwise           | 1/7 ✗  | 0.9490  | 0.4      | ✗        |
| WSINDy + RFE          | polynomial, pointwise           | 1/7 ✗  | 0.9490  | 0.4      | ✗        |
| PolyJet + RFE         | polynomial, pointwise           | 1/7 ✗  | 0.9490  | 0.7      | ✗        |
| PDE-FIND + MSTLS      | polynomial, weak (IBP), m=12    | 1/7 ✗  | 0.9490  | 3.8      | ✗        |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=12 | 1/7 ✗  | 0.9787  | 1.8      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 30/2  | 43       | 0.6827  | 0.4383 | 6.1      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 1/7   | 1        | 0.4515  | 0.4404 | 0.4      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.01    | 38/2  | 57       | 0.9990  | 0.6593 | 7.3      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.6602 | 0.4      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 34/2  | 51       | 0.9991  | 0.8833 | 7.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.8839 | 0.5      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.01    | 83/0  | 43       | 0.4503  | 0.4387 | 4.9      | 🔴      |
| WSINDy      | Clean | RFE       |         | 1/7   | 1        | 0.4515  | 0.4404 | 0.4      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 86/0  | 42       | 0.9997  | 0.6594 | 5.4      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.6602 | 0.3      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 89/0  | 44       | 0.9998  | 0.8833 | 4.8      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.8839 | 0.3      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 0/0   | 7        | 0.0058  | 0.0197 | 12.7     | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.0056  | 0.0197 | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.1     | 65/0  | 72       | 0.9997  | 0.6516 | 18.9     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.6602 | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 87/0  | 94       | 0.9998  | 0.8824 | 22.3     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.8839 | 0.9      | 🔴      |
| Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0057  | 0.0197 | 1.6      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0055  | 0.0197 | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 17/0  | 11       | 0.9997  | 0.6535 | 1.6      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.6602 | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 18/0  | 12       | 0.9998  | 0.8839 | 1.7      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.8839 | 0.2      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 28/2  | 42       | 0.7101  | 0.4384 | 6.7      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 1/7   | 1        | 0.5227  | 0.4404 | 0.7      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.01    | 35/2  | 52       | 0.9056  | 0.6593 | 6.9      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/7   | 1        | 0.8908  | 0.6602 | 0.4      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 32/2  | 49       | 0.9096  | 0.8831 | 6.6      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/7   | 1        | 0.9490  | 0.8839 | 0.4      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.01    | 94/0  | 43       | 0.5211  | 0.4387 | 4.5      | 🔴      |
| WSINDy      | Clean | RFE       |         | 1/7   | 1        | 0.5227  | 0.4404 | 0.3      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 96/0  | 44       | 0.8915  | 0.6594 | 5.2      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 1/7   | 1        | 0.8908  | 0.6602 | 0.3      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.01    | 89/0  | 44       | 0.9508  | 0.8832 | 4.9      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/7   | 1        | 0.9490  | 0.8839 | 0.4      | 🔴      |
| PolyJet     | Clean | STRidge   | 10.0    | 2/0   | 9        | 0.2174  | 0.2119 | 12.4     | 🔴      |
| PolyJet     | Clean | RFE       |         | 0/6   | 1        | 0.6498  | 0.7076 | 0.7      | 🔴      |
| PolyJet     | 2%    | STRidge   | 0.1     | 62/0  | 69       | 0.8601  | 0.6126 | 18.5     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/7   | 1        | 0.8908  | 0.6602 | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 84/0  | 91       | 0.9437  | 0.8784 | 21.0     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/7   | 1        | 0.9490  | 0.8839 | 0.7      | 🔴      |
| Multigraphs | Clean | STRidge   | 1.0     | 10/0  | 6        | 0.2774  | 0.1129 | 1.6      | 🔴      |
| Multigraphs | Clean | RFE       |         | 0/5   | 1        | 0.3549  | 0.5044 | 0.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 0.01    | 33/0  | 19       | 0.8591  | 0.6178 | 1.9      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.8908  | 0.6602 | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.001   | 37/0  | 20       | 0.9432  | 0.8801 | 1.9      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.9490  | 0.8839 | 0.2      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 1.0     | 8  | no  | 12/2  | 20       | 5.1566  | 0.2947 | 5.4      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 4.1      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 0.4      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 1.0     | 12 | yes | 16/2  | 24       | 6.1202  | 0.2818 | 5.3      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2889 | 3.6      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2889 | 0.5      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 4  | no  | 10/2  | 18       | 0.9752  | 0.4183 | 5.4      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9997  | 0.2431 | 3.1      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9997  | 0.2431 | 0.5      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 8  | yes | 11/2  | 20       | 0.9846  | 0.3797 | 5.2      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9997  | 0.3004 | 3.3      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9997  | 0.3004 | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 4  | no  | 12/2  | 19       | 0.9988  | 0.5393 | 5.7      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9997  | 0.2788 | 3.5      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9997  | 0.2788 | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 8  | yes | 14/2  | 21       | 0.9930  | 0.4719 | 5.4      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9997  | 0.3519 | 3.6      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9997  | 0.3519 | 0.6      | 🔴      |
| WSINDy      | Clean | STRidge   | 1.0     | 4  | no  | 84/0  | 39       | 0.4528  | 0.3840 | 4.7      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7177  | 0.6516 | 2.4      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 0.3      | 🔴      |
| WSINDy      | Clean | STRidge   | 1.0     | 12 | yes | 71/0  | 34       | 0.4838  | 0.2838 | 4.6      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7216  | 0.6696 | 2.7      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2889 | 0.3      | 🔴      |
| WSINDy      | 2%    | STRidge   | 1.0     | 4  | no  | 87/0  | 41       | 1.0004  | 0.4179 | 4.9      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9999  | 0.6540 | 2.0      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9997  | 0.2431 | 0.3      | 🔴      |
| WSINDy      | 2%    | STRidge   | 1.0     | 12 | yes | 69/0  | 34       | 1.0959  | 0.2957 | 4.6      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9999  | 0.6738 | 2.5      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9997  | 0.3004 | 0.3      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 12 | no  | 83/0  | 42       | 1.2747  | 0.2731 | 4.3      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9999  | 0.6660 | 2.5      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9997  | 0.2788 | 0.3      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 8  | yes | 98/0  | 46       | 0.9991  | 0.4714 | 4.7      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9999  | 0.6942 | 2.7      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9997  | 0.3519 | 0.3      | 🔴      |
| PolyJet     | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 7        | 0.0021  | 0.0095 | 11.5     | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 7        | 0.0019  | 0.0095 | 39.4     | 🟢      |
| PolyJet     | Clean | RFE       |         | 8  | no  | 0/0   | 7        | 0.0019  | 0.0095 | 0.8      | 🟢      |
| PolyJet     | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 7        | 0.0029  | 0.0127 | 11.9     | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 7        | 0.0024  | 0.0087 | 38.4     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0024  | 0.0087 | 0.7      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | no  | 57/0  | 64       | 0.3901  | 0.1400 | 16.7     | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 8  | no  | 3/0   | 10       | 0.4115  | 0.2057 | 39.3     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 0/6   | 1        | 0.9995  | 0.2326 | 0.8      | 🔴      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | yes | 41/0  | 48       | 0.4254  | 0.1878 | 15.0     | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 4/0   | 11       | 0.4053  | 0.1904 | 38.2     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 0.9995  | 0.3011 | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 8  | no  | 36/0  | 43       | 0.9689  | 0.3587 | 15.2     | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9997  | 0.2788 | 37.7     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9997  | 0.2788 | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | yes | 38/0  | 45       | 0.9862  | 0.3397 | 14.5     | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9997  | 0.3519 | 35.9     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9997  | 0.3519 | 0.8      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 1.7      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.8      | 🟢      |
| Multigraphs | Clean | RFE       |         | 4  | no  | 0/0   | 3        | 0.0037  | 0.0152 | 0.2      | 🟢      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0022  | 0.0087 | 1.6      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0020  | 0.0087 | 0.8      | 🟢      |
| Multigraphs | Clean | RFE       |         | 8  | yes | 0/0   | 3        | 0.0028  | 0.0127 | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 10.0    | 8  | no  | 25/0  | 14       | 0.4829  | 0.2148 | 1.7      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 12 | no  | 5/0   | 5        | 0.3943  | 0.1551 | 0.8      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9994  | 0.2197 | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 10.0    | 8  | yes | 14/0  | 10       | 0.8389  | 0.2983 | 1.9      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.5     | 12 | yes | 4/0   | 4        | 0.4568  | 0.2039 | 0.7      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9994  | 0.2785 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 4  | no  | 30/0  | 17       | 0.9985  | 0.5367 | 1.7      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9997  | 0.2788 | 0.6      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.9995  | 0.3968 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 8  | yes | 32/0  | 18       | 0.9953  | 0.4681 | 1.8      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9997  | 0.3519 | 0.6      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9997  | 0.3519 | 0.2      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 1.0     | 8  | no  | 12/2  | 20       | 4.0849  | 0.2950 | 5.7      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.5226  | 0.3865 | 4.8      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.3865 | 0.4      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 1.0     | 12 | yes | 16/2  | 24       | 4.8686  | 0.2822 | 5.7      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.5226  | 0.3592 | 4.3      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5226  | 0.3592 | 0.4      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 8  | no  | 12/2  | 21       | 5.0461  | 0.3083 | 5.7      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.8908  | 0.2431 | 3.6      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.8908  | 0.2431 | 0.4      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | yes | 15/2  | 24       | 5.6037  | 0.2946 | 5.7      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.8908  | 0.3004 | 4.3      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.8908  | 0.3004 | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 10.0    | 12 | no  | 10/2  | 19       | 1.4514  | 0.2734 | 5.4      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9491  | 0.2788 | 3.2      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9491  | 0.2788 | 0.4      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | yes | 16/2  | 26       | 1.0279  | 0.3480 | 5.4      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9490  | 0.3519 | 3.8      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9490  | 0.3519 | 0.5      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.01    | 4  | no  | 94/0  | 45       | 0.5223  | 0.3840 | 4.9      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7444  | 0.6516 | 2.4      | 🔴      |
| WSINDy      | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.3865 | 0.3      | 🔴      |
| WSINDy      | Clean | STRidge   | 1.0     | 12 | yes | 84/0  | 41       | 0.5188  | 0.2838 | 4.2      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7479  | 0.6696 | 2.2      | 🔴      |
| WSINDy      | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5226  | 0.3592 | 0.4      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 4  | no  | 64/0  | 32       | 0.8885  | 0.4179 | 5.2      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9439  | 0.6540 | 2.7      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.8908  | 0.2431 | 0.4      | 🔴      |
| WSINDy      | 2%    | STRidge   | 1.0     | 12 | yes | 65/0  | 33       | 0.9801  | 0.2956 | 4.3      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9447  | 0.6738 | 2.4      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.8908  | 0.3004 | 0.3      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 12 | no  | 83/0  | 42       | 1.2103  | 0.2731 | 4.3      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9787  | 0.6660 | 1.8      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9491  | 0.2788 | 0.4      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 8  | yes | 98/0  | 46       | 0.9534  | 0.4714 | 4.6      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9792  | 0.6942 | 1.9      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9490  | 0.3519 | 0.3      | 🔴      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0335  | 0.0071 | 12.0     | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 7        | 0.0392  | 0.0110 | 36.0     | 🟢      |
| PolyJet     | Clean | RFE       |         | 8  | no  | 0/0   | 7        | 0.0392  | 0.0110 | 0.8      | 🟢      |
| PolyJet     | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 7        | 0.0287  | 0.0146 | 11.9     | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 7        | 0.0288  | 0.0146 | 37.2     | 🟢      |
| PolyJet     | Clean | RFE       |         | 8  | yes | 0/0   | 7        | 0.0288  | 0.0146 | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 10.0    | 8  | no  | 13/0  | 20       | 0.3970  | 0.1447 | 13.2     | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.5     | 4  | no  | 8/0   | 15       | 0.7204  | 0.2721 | 40.1     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 8  | no  | 0/6   | 1        | 0.8912  | 0.3243 | 0.7      | 🔴      |
| PolyJet     | 2%    | STRidge   | 10.0    | 12 | yes | 19/0  | 26       | 0.3949  | 0.1351 | 13.2     | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 9/0   | 16       | 0.3918  | 0.1359 | 39.4     | 🔴      |
| PolyJet     | 2%    | RFE       |         | 8  | yes | 0/6   | 1        | 0.8870  | 0.4623 | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 8  | no  | 69/0  | 76       | 0.7383  | 0.3011 | 19.4     | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 4  | no  | 8/2   | 13       | 0.9190  | 0.4999 | 38.7     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9491  | 0.2788 | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 10.0    | 12 | yes | 64/0  | 71       | 0.7462  | 0.2840 | 15.6     | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | yes | 9/2   | 14       | 0.8041  | 0.2869 | 40.7     | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9490  | 0.3519 | 0.7      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0386  | 0.0115 | 1.7      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.75    | 12 | no  | 0/0   | 3        | 0.0350  | 0.0071 | 0.6      | 🟢      |
| Multigraphs | Clean | RFE       |         | 4  | no  | 0/5   | 1        | 0.3492  | 0.2638 | 0.2      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 3        | 0.0277  | 0.0156 | 1.7      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 8  | yes | 0/0   | 3        | 0.0277  | 0.0156 | 0.8      | 🟢      |
| Multigraphs | Clean | RFE       |         | 8  | yes | 0/0   | 3        | 0.0277  | 0.0156 | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 10.0    | 4  | no  | 9/0   | 7        | 0.6901  | 0.2792 | 1.6      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 5/0   | 5        | 0.2663  | 0.1483 | 1.0      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.8264  | 0.1912 | 0.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 10.0    | 12 | yes | 22/0  | 12       | 0.3015  | 0.1368 | 1.8      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 5/0   | 5        | 0.2721  | 0.1395 | 0.8      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.8261  | 0.1803 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 4  | no  | 32/0  | 17       | 0.8969  | 0.5044 | 1.7      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 8  | no  | 9/3   | 6        | 0.9463  | 0.3114 | 1.1      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9016  | 0.2652 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.001   | 12 | yes | 37/0  | 20       | 0.6972  | 0.2916 | 1.7      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | yes | 9/3   | 6        | 0.9468  | 0.2945 | 0.7      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9026  | 0.3432 | 0.1      | 🔴      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ)     |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:------------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0156       | 37240999.0  |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0096       | 103883766.1 |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0066       | 317175677.7 |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.0891       | 39715936.0  |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0172       | 107952857.1 |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0115       | 320266523.5 |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0130       | 50038749.7  |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0089       | 118142901.4 |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.0239       | 46089832.6  |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0150       | 117691652.7 |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 3.4737       | 39795395.5  |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 3.3691       | 90842207.1  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 3.4010       | 273651905.4 |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 0.4838       | 35806846.8  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.3048       | 100771904.5 |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.2790       | 293455352.7 |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 3.3618       | 44766235.9  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 3.3725       | 104180806.6 |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 0.4215       | 38412858.0  |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.2983       | 106432775.4 |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 19.5646      | 39785550.0  |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 20.5720      | 58952277.8  |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 21.0208      | 124183762.7 |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 1.7861       | 35711721.6  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 1.6500       | 72025688.0  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 1.6411       | 153225710.5 |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 19.9290      | 44841209.3  |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 20.6511      | 62264823.2  |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 1.7507       | 38262359.9  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 1.6435       | 72801128.4  |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ks_chaos`.

| Source                                    | Wall time | CPU                                      | Threads |
|:------------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260731_023014_ks_chaos_full.csv` | 474.9 min | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | not set |
