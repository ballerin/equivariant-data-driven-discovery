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
| Clean | pointwise     | polynomial | —        | 7/0 ✗             | 0.2846     | —                  | 0.1198    | —               | 2/0 ✗     | 0.2174      | 0.2119     | 0.3112       |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0017     | multigraph/STRidge | 0.0095    | PolyJet/MSTLS   | 0/0 ✓     | 0.0019      | 0.0095     | 0.0066       |
| Clean | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.0350     | multigraph/MSTLS   | 0.0071    | PolyJet/STRidge | 0/0 ✓     | 0.0335      | 0.0071     | 0.0115       |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0057     | multigraph/MSTLS   | 0.0063    | PolyJet/STRidge | 0/0 ✓     | 0.0061      | 0.0063     | 0.0069       |
| Clean | weak (IBP)    | polynomial | 🟥        | 0/0 ✓             | 0.0340     | multigraph/MSTLS   | 0.0070    | PolyJet/MSTLS   | 0/0 ✓     | 0.0324      | 0.0070     | 0.0101       |
| 2%    | pointwise     | spectral   | —        | 17/0 ✗            | 0.9997     | —                  | 0.6535    | —               | 65/0 ✗    | 0.9997      | 0.6516     | 85.8534      |
| 2%    | pointwise     | polynomial | —        | 37/0 ✗            | 0.8954     | —                  | 0.6216    | —               | 62/0 ✗    | 0.8601      | 0.6126     | 1.8548       |
| 2%    | weak (no IBP) | spectral   | —        | 5/0 ✗             | 0.3943     | —                  | 0.1551    | —               | 3/0 ✗     | 0.4115      | 0.2057     | 3.3691       |
| 2%    | weak (no IBP) | polynomial | —        | 7/0 ✗             | 0.6507     | —                  | 0.1478    | —               | 8/0 ✗     | 0.7204      | 0.2721     | 0.2790       |
| 2%    | weak (IBP)    | spectral   | —        | 4/0 ✗             | 0.4601     | —                  | 0.2151    | —               | 4/0 ✗     | 0.4155      | 0.2022     | 3.4030       |
| 2%    | weak (IBP)    | polynomial | —        | 3/0 ✗             | 0.2289     | —                  | 0.1552    | —               | 9/0 ✗     | 0.4393      | 0.1446     | 0.2807       |
| 5%    | pointwise     | spectral   | —        | 18/0 ✗            | 0.9998     | —                  | 0.8839    | —               | 87/0 ✗    | 0.9998      | 0.8824     | 133.9474     |
| 5%    | pointwise     | polynomial | —        | 24/0 ✗            | 0.9629     | —                  | 0.8812    | —               | 84/0 ✗    | 0.9437      | 0.8784     | 3.0277       |
| 5%    | weak (no IBP) | spectral   | —        | 30/0 ✗            | 0.9985     | —                  | 0.5367    | —               | 36/0 ✗    | 0.9689      | 0.3587     | 19.5646      |
| 5%    | weak (no IBP) | polynomial | —        | 30/0 ✗            | 1.0068     | —                  | 0.5053    | —               | 69/0 ✗    | 0.7383      | 0.3011     | 1.6411       |
| 5%    | weak (IBP)    | spectral   | —        | 25/0 ✗            | 0.9999     | —                  | 0.5382    | —               | 36/0 ✗    | 0.9700      | 0.3587     | 20.6187      |
| 5%    | weak (IBP)    | polynomial | —        | 32/0 ✗            | 1.2217     | —                  | 0.3194    | —               | 38/0 ✗    | 0.7384      | 0.3015     | 1.6430       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-------------------------------|:-----------------------------|:------|--------:|---------:|---------:|
| Multigraphs + STRidge          | spectral, weak (no IBP), m=8 | 0/0 ✓ | 0.0017  | 0.2      | —        |
| Modified Multigraphs + STRidge | spectral, weak (no IBP), m=8 | 0/0 ✓ | 0.0017  | 0.2      | +0.0%    |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=8 | 0/0 ✓ | 0.0017  | 0.4      | +1.0%    |
| Modified Multigraphs + MSTLS   | spectral, weak (no IBP), m=8 | 0/0 ✓ | 0.0017  | 0.4      | +1.0%    |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=8 | 0/0 ✓ | 0.0019  | 3.4      | +14.5%   |
| PolyJet + RFE                  | spectral, weak (no IBP), m=8 | 0/0 ✓ | 0.0019  | 0.4      | +14.5%   |
| PolyJet + STRidge              | spectral, weak (no IBP), m=8 | 0/0 ✓ | 0.0021  | 0.8      | +23.6%   |
| Multigraphs + RFE              | spectral, weak (no IBP), m=4 | 0/0 ✓ | 0.0037  | 0.2      | +121.5%  |
| Modified Multigraphs + RFE     | spectral, weak (no IBP), m=4 | 0/0 ✓ | 0.0037  | 0.2      | +121.5%  |

### 2% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=8     | 3/0 ✗  | 0.2289  | 0.4      | ✗        |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=8    | 3/0 ✗  | 0.4115  | 2.9      | ✗        |
| Multigraphs + MSTLS            | spectral, weak (IBP), m=8       | 4/0 ✗  | 0.4601  | 0.4      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=4  | 8/0 ✗  | 0.9584  | 0.2      | ✗        |
| Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=4  | 8/0 ✗  | 0.9584  | 0.2      | ✗        |
| PolyJet + STRidge              | polynomial, weak (IBP), m=8     | 13/0 ✗ | 0.3144  | 0.8      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.8258  | 0.2      | ✗        |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=8  | 0/5 ✗  | 0.8264  | 0.2      | ✗        |
| PolyJet + RFE                  | polynomial, weak (no IBP), m=8  | 0/6 ✗  | 0.8912  | 0.3      | ✗        |

### 5% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + STRidge          | spectral, pointwise             | 18/0 ✗ | 0.9998  | 0.3      | ✗        |
| Modified Multigraphs + STRidge | spectral, weak (IBP), m=4       | 25/0 ✗ | 0.9999  | 0.3      | ✗        |
| PolyJet + STRidge              | spectral, weak (no IBP), m=8    | 36/0 ✗ | 0.9689  | 1.0      | ✗        |
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=4  | 8/2 ✗  | 0.9190  | 3.9      | ✗        |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.9016  | 0.2      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.9016  | 0.2      | ✗        |
| Modified Multigraphs + MSTLS   | polynomial, weak (no IBP), m=8  | 1/5 ✗  | 0.9204  | 0.4      | ✗        |
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=4  | 5/5 ✗  | 0.9951  | 0.4      | ✗        |
| PolyJet + RFE                  | polynomial, pointwise           | 1/7 ✗  | 0.9490  | 0.4      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PolyJet          | Clean | STRidge   | 1.0     | 0/0   | 7        | 0.0058  | 0.0197 | 1.2      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0056  | 0.0197 | 0.5      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.1     | 65/0  | 72       | 0.9997  | 0.6516 | 1.6      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.6602 | 0.4      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 87/0  | 94       | 0.9998  | 0.8824 | 1.7      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.8839 | 0.4      | 🔴      |
| Multigraphs      | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0057  | 0.0197 | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0055  | 0.0197 | 0.2      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 17/0  | 11       | 0.9997  | 0.6535 | 0.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.6602 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 18/0  | 12       | 0.9998  | 0.8839 | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.8839 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0057  | 0.0197 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0055  | 0.0197 | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 31/0  | 17       | 0.9997  | 0.6534 | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.6602 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.01    | 34/0  | 19       | 0.9998  | 0.8838 | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.9997  | 0.8839 | 0.2      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PolyJet          | Clean | STRidge   | 10.0    | 2/0   | 9        | 0.2174  | 0.2119 | 1.1      | 🔴      |
| PolyJet          | Clean | RFE       |         | 0/6   | 1        | 0.6498  | 0.7076 | 0.5      | 🔴      |
| PolyJet          | 2%    | STRidge   | 0.1     | 62/0  | 69       | 0.8601  | 0.6126 | 1.5      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 0.8908  | 0.6602 | 0.5      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 84/0  | 91       | 0.9437  | 0.8784 | 1.6      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 0.9490  | 0.8839 | 0.4      | 🔴      |
| Multigraphs      | Clean | STRidge   | 1.0     | 7/0   | 5        | 0.2856  | 0.1196 | 0.3      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 0/5   | 1        | 0.3549  | 0.5044 | 0.2      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 0.001   | 37/0  | 20       | 0.8954  | 0.6216 | 0.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 0.8908  | 0.6602 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 24/0  | 14       | 0.9629  | 0.8812 | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 0.9490  | 0.8839 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 7/0   | 4        | 0.2846  | 0.1198 | 0.3      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 7/0   | 4        | 0.2851  | 0.1198 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 0.001   | 37/0  | 20       | 0.8954  | 0.6216 | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.8908  | 0.6602 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 37/0  | 20       | 0.9628  | 0.8810 | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.9490  | 0.8839 | 0.2      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PolyJet          | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 7        | 0.0021  | 0.0095 | 0.8      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 7        | 0.0019  | 0.0095 | 3.4      | 🟢      |
| PolyJet          | Clean | RFE       |         | 8  | no  | 0/0   | 7        | 0.0019  | 0.0095 | 0.4      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0061  | 0.0063 | 0.6      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.5     | 8  | yes | 0/0   | 7        | 0.0836  | 0.0222 | 3.0      | 🟢      |
| PolyJet          | Clean | RFE       |         | 4  | yes | 1/0   | 8        | 17.2773 | 0.0204 | 0.4      | 🔴      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | no  | 57/0  | 64       | 0.3901  | 0.1400 | 0.7      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 8  | no  | 3/0   | 10       | 0.4115  | 0.2057 | 2.9      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/6   | 1        | 0.9995  | 0.2326 | 0.2      | 🔴      |
| PolyJet          | 2%    | STRidge   | 1.0     | 4  | yes | 54/0  | 61       | 0.9761  | 0.3460 | 1.4      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 8  | yes | 4/0   | 11       | 0.4155  | 0.2022 | 2.8      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 0.9995  | 0.2326 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | no  | 36/0  | 43       | 0.9689  | 0.3587 | 1.0      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9997  | 0.2788 | 1.9      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9997  | 0.2788 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | yes | 36/0  | 43       | 0.9700  | 0.3587 | 1.0      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9997  | 0.2788 | 1.9      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9997  | 0.2788 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.2      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.4      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 4  | no  | 0/0   | 3        | 0.0037  | 0.0152 | 0.2      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | yes | 1/0   | 4        | 0.0509  | 0.0090 | 0.2      | 🔴      |
| Multigraphs      | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 3        | 0.0057  | 0.0063 | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 4  | yes | 1/0   | 4        | 15.0047 | 0.0435 | 0.2      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 8  | no  | 25/0  | 14       | 0.4829  | 0.2148 | 0.2      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 12 | no  | 5/0   | 5        | 0.3943  | 0.1551 | 0.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9994  | 0.2197 | 0.2      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 8  | yes | 25/0  | 14       | 0.4787  | 0.2148 | 0.2      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 8  | yes | 4/0   | 4        | 0.4601  | 0.2151 | 0.4      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9994  | 0.2197 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 4  | no  | 30/0  | 17       | 0.9985  | 0.5367 | 0.3      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9997  | 0.2788 | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.9995  | 0.3968 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 4  | yes | 27/0  | 16       | 0.9982  | 0.5367 | 0.3      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9997  | 0.2788 | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 0/5   | 1        | 0.9995  | 0.3968 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 4  | no  | 0/0   | 3        | 0.0037  | 0.0152 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 8  | yes | 1/0   | 4        | 0.0509  | 0.0090 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 3        | 0.0057  | 0.0063 | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 4  | yes | 1/0   | 4        | 15.0047 | 0.0435 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | no  | 29/0  | 14       | 0.3404  | 0.1535 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 8/0   | 6        | 0.4590  | 0.2153 | 0.4      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9994  | 0.2197 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | yes | 17/0  | 12       | 0.5748  | 0.1733 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 9/0   | 7        | 0.6217  | 0.1749 | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9994  | 0.2197 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | no  | 33/0  | 17       | 1.0235  | 0.2696 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.9997  | 0.2788 | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9997  | 0.2788 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 4  | yes | 25/0  | 16       | 0.9999  | 0.5382 | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.9997  | 0.2788 | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9997  | 0.2788 | 0.2      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PolyJet          | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0335  | 0.0071 | 0.6      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 7        | 0.0392  | 0.0110 | 2.9      | 🟢      |
| PolyJet          | Clean | RFE       |         | 8  | no  | 0/0   | 7        | 0.0392  | 0.0110 | 0.3      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0326  | 0.0070 | 0.6      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 7        | 0.0324  | 0.0070 | 2.1      | 🟢      |
| PolyJet          | Clean | RFE       |         | 4  | yes | 1/0   | 8        | 16.6548 | 0.0229 | 0.4      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 8  | no  | 13/0  | 20       | 0.3970  | 0.1447 | 0.9      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 4  | no  | 8/0   | 15       | 0.7204  | 0.2721 | 4.1      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 8  | no  | 0/6   | 1        | 0.8912  | 0.3243 | 0.3      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 8  | yes | 13/0  | 20       | 0.3144  | 0.1474 | 0.8      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 8  | yes | 9/0   | 16       | 0.4393  | 0.1446 | 2.9      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 8  | yes | 0/6   | 1        | 0.8912  | 0.3243 | 0.3      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | no  | 69/0  | 76       | 0.7383  | 0.3011 | 1.1      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 4  | no  | 8/2   | 13       | 0.9190  | 0.4999 | 3.9      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9491  | 0.2788 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | yes | 38/0  | 45       | 0.7384  | 0.3015 | 1.0      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 8  | yes | 9/2   | 14       | 0.9563  | 0.3036 | 3.1      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9491  | 0.2788 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | no  | 4/0   | 5        | 0.0498  | 0.0066 | 0.2      | 🔴      |
| Multigraphs      | Clean | MSTLS     | 0.75    | 12 | no  | 0/0   | 3        | 0.0350  | 0.0071 | 0.4      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 4  | no  | 0/5   | 1        | 0.3492  | 0.2638 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | yes | 1/0   | 4        | 0.0524  | 0.0104 | 0.2      | 🔴      |
| Multigraphs      | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 3        | 0.0340  | 0.0070 | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 4  | yes | 1/0   | 4        | 14.8344 | 0.0395 | 0.2      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 4  | no  | 8/0   | 6        | 0.9584  | 0.2665 | 0.2      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 8  | no  | 7/0   | 5        | 0.6507  | 0.1478 | 0.4      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.8264  | 0.1912 | 0.2      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 8  | yes | 19/0  | 12       | 0.6449  | 0.1462 | 0.2      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | yes | 0/5   | 1        | 0.8258  | 0.1301 | 0.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 8  | yes | 0/5   | 1        | 0.8264  | 0.1912 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 4  | no  | 31/0  | 18       | 1.0068  | 0.5053 | 0.3      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 4  | no  | 5/5   | 3        | 0.9951  | 0.5080 | 0.4      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9016  | 0.2652 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 8  | yes | 34/0  | 19       | 0.8454  | 0.3160 | 0.2      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.5     | 4  | yes | 5/5   | 3        | 0.9951  | 0.5080 | 0.5      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9016  | 0.2652 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0386  | 0.0115 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.75    | 12 | no  | 0/0   | 3        | 0.0350  | 0.0071 | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 4  | no  | 0/5   | 1        | 0.3492  | 0.2638 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 8  | yes | 1/0   | 4        | 0.0524  | 0.0104 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 3        | 0.0340  | 0.0070 | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 4  | yes | 1/0   | 4        | 14.8344 | 0.0395 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 4  | no  | 8/0   | 6        | 0.9584  | 0.2665 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 8/0   | 6        | 0.6386  | 0.1475 | 0.4      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.8258  | 0.1301 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 8  | yes | 20/0  | 13       | 0.9092  | 0.1455 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 8  | yes | 3/0   | 4        | 0.2289  | 0.1552 | 0.4      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.8258  | 0.1301 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 4  | no  | 30/0  | 17       | 1.0068  | 0.5053 | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.5     | 8  | no  | 1/5   | 2        | 0.9204  | 0.3449 | 0.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9016  | 0.2652 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 8  | yes | 32/0  | 18       | 1.2217  | 0.3194 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.5     | 8  | yes | 4/5   | 3        | 0.8982  | 0.3256 | 0.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9016  | 0.2652 | 0.2      | 🔴      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ)     |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:------------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0156       | 37240999.0  |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0096       | 103883766.1 |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0066       | 317175677.7 |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.0891       | 39807836.4  |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0172       | 108519130.8 |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0115       | 318159697.5 |
| Clean | weak (IBP)    | spectral   | 4  | yes | 12.4492      | 46753006.8  |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0548       | 105059979.5 |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0069       | 317477220.6 |
| Clean | weak (IBP)    | polynomial | 4  | yes | 12.4415      | 40706924.2  |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.0481       | 104746960.5 |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0101       | 315952930.8 |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 3.4737       | 39795395.5  |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 3.3691       | 90842207.1  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 3.4010       | 273651905.4 |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 0.4838       | 36907183.3  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.3048       | 100564243.2 |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.2790       | 304727896.7 |
| 2%    | weak (IBP)    | spectral   | 4  | yes | 15.2311      | 39797489.0  |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 3.4168       | 91745285.8  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 3.4030       | 273959604.1 |
| 2%    | weak (IBP)    | polynomial | 4  | yes | 12.4899      | 36906281.6  |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 0.3460       | 91189739.9  |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.2807       | 302805745.6 |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 19.5646      | 39785550.0  |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 20.5720      | 58952277.8  |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 21.0208      | 124183762.7 |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 1.7861       | 37606173.9  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 1.6500       | 75976318.7  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 1.6411       | 179198493.1 |
| 5%    | weak (IBP)    | spectral   | 4  | yes | 29.7857      | 39566356.0  |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 20.6187      | 59009394.7  |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 21.0227      | 124299300.8 |
| 5%    | weak (IBP)    | polynomial | 4  | yes | 12.7976      | 37557011.2  |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 1.6947       | 64171012.7  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 1.6430       | 178146133.8 |
---
| Source                                       | Wall time | CPU                                           | Threads |
|:---------------------------------------------|----------:|:----------------------------------------------|:-------:|
| `sweep_20260715_212141_ks_chaos_trimmed.csv` | 43.2 min  | AMD Ryzen 7 PRO 7840U w/ Radeon 780M Graphics | 1       |
