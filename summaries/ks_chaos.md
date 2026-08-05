# Kuramoto-Sivashinsky equation (strongly chaotic, $L=16\pi$)

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

KS contains the gradient norm $\|\nabla f\|^2 = f_x^2 + f_y^2$ which neither PDE-FIND nor WSINDy can represent as a sparse column (PDE-FIND lacks products of derivatives; in WSINDy, $f_x^2$ is locked with $f\cdot f_{xx}$ inside $\partial_{xx}(f^2)$).  Chaotic regime ($L=16\pi$).

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
| Clean | pointwise     | spectral   | 🟩        | 0/0 ✓             | 0.0055     | mod_multigraph/RFE     | 0.0197    | PolyJet/RFE     | 0/0 ✓     | 0.0056      | 0.0197     | 0.0205       |
| Clean | pointwise     | polynomial | 🟢        | 0/0 ✓             | 0.1871     | mod_multigraph/STRidge | 0.2729    | —               | 2/0 ✗     | 0.2946      | 0.2119     | 0.3112       |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0017     | multigraph/STRidge     | 0.0095    | PolyJet/MSTLS   | 0/0 ✓     | 0.0019      | 0.0095     | 0.0066       |
| Clean | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0423     | multigraph/MSTLS       | 0.0071    | PolyJet/MSTLS   | 0/0 ✓     | 0.0523      | 0.0110     | 0.0115       |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0057     | multigraph/MSTLS       | 0.0063    | PolyJet/STRidge | 0/0 ✓     | 0.0061      | 0.0063     | 0.0069       |
| Clean | weak (IBP)    | polynomial | 🟥        | 0/0 ✓             | 0.0409     | multigraph/MSTLS       | 0.0070    | PolyJet/MSTLS   | 0/0 ✓     | 0.0408      | 0.0070     | 0.0101       |
| 2%    | pointwise     | spectral   | —        | 17/0 ✗            | 0.4234     | —                      | 0.6535    | —               | 47/0 ✗    | 0.4212      | 0.6518     | 85.8534      |
| 2%    | pointwise     | polynomial | —        | 31/0 ✗            | 0.3084     | —                      | 0.6179    | —               | 62/0 ✗    | 0.2942      | 0.6127     | 1.8548       |
| 2%    | weak (no IBP) | spectral   | —        | 5/0 ✗             | 0.3918     | —                      | 0.3575    | —               | 3/0 ✗     | 0.3483      | 0.2057     | 3.3691       |
| 2%    | weak (no IBP) | polynomial | —        | 3/0 ✗             | 0.3515     | —                      | 0.1548    | —               | 8/0 ✗     | 0.2652      | 0.2721     | 0.2790       |
| 2%    | weak (IBP)    | spectral   | —        | 4/0 ✗             | 0.4356     | —                      | 0.2151    | —               | 4/0 ✗     | 0.3596      | 0.2022     | 3.4030       |
| 2%    | weak (IBP)    | polynomial | —        | 3/0 ✗             | 0.3170     | —                      | 0.1552    | —               | 9/0 ✗     | 0.2786      | 0.1446     | 0.2807       |
| 5%    | pointwise     | spectral   | —        | 24/0 ✗            | 0.4506     | —                      | 0.8838    | —               | 77/0 ✗    | 0.4500      | 0.8824     | 133.9474     |
| 5%    | pointwise     | polynomial | —        | 37/0 ✗            | 0.4056     | —                      | 0.8801    | —               | 54/0 ✗    | 0.4004      | 0.8785     | 3.0277       |
| 5%    | weak (no IBP) | spectral   | —        | 30/0 ✗            | 0.4511     | —                      | 0.5367    | —               | 27/0 ✗    | 0.4451      | 0.5352     | 19.5646      |
| 5%    | weak (no IBP) | polynomial | —        | 9/0 ✗             | 0.3619     | —                      | 0.5053    | —               | 59/0 ✗    | 0.6333      | 0.2733     | 1.6411       |
| 5%    | weak (IBP)    | spectral   | —        | 32/0 ✗            | 0.4804     | —                      | 0.3616    | —               | 53/0 ✗    | 0.4781      | 0.3584     | 20.6187      |
| 5%    | weak (IBP)    | polynomial | —        | 25/0 ✗            | 0.9582     | —                      | 0.2327    | —               | 62/0 ✗    | 0.3252      | 0.3012     | 1.6430       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + STRidge          | spectral, weak (no IBP), m=8    | 0/0 ✓  | 0.0017  | 0.4      | —        |
| Modified Multigraphs + STRidge | spectral, weak (no IBP), m=8    | 0/0 ✓  | 0.0017  | 0.6      | +0.0%    |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=8    | 0/0 ✓  | 0.0017  | 0.6      | +1.0%    |
| Modified Multigraphs + MSTLS   | spectral, weak (no IBP), m=8    | 0/0 ✓  | 0.0017  | 0.8      | +1.0%    |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=8    | 0/0 ✓  | 0.0019  | 5.2      | +14.5%   |
| PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/0 ✓  | 0.0019  | 0.2      | +14.5%   |
| Modified Multigraphs + RFE     | spectral, weak (no IBP), m=4    | 0/0 ✓  | 0.0037  | 0.2      | +121.5%  |
| Multigraphs + RFE              | spectral, weak (no IBP), m=4    | 0/0 ✓  | 0.0037  | 0.1      | +121.5%  |
| PolyJet + STRidge              | spectral, weak (no IBP), m=4    | 0/0 ✓  | 0.0038  | 2.2      | +126.3%  |
| WSINDy + STRidge               | spectral, weak (IBP), m=8       | 67/0 ✗ | 0.4819  | 1.6      | ✗        |
| PDE-FIND + STRidge             | polynomial, weak (no IBP), m=12 | 19/2 ✗ | 2.1454  | 1.0      | ✗        |
| PDE-FIND + RFE                 | polynomial, weak (no IBP), m=12 | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| PDE-FIND + MSTLS               | polynomial, weak (no IBP), m=12 | 1/7 ✗  | 0.4515  | 1.8      | ✗        |
| WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.7177  | 0.9      | ✗        |

### 2% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=8     | 3/0 ✗  | 0.3170  | 0.6      | ✗        |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=8    | 3/0 ✗  | 0.3483  | 11.4     | ✗        |
| Multigraphs + MSTLS            | spectral, weak (IBP), m=8       | 4/0 ✗  | 0.4356  | 0.7      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=4  | 11/0 ✗ | 0.2188  | 0.7      | ✗        |
| Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=4  | 11/0 ✗ | 0.2190  | 0.8      | ✗        |
| PolyJet + STRidge              | polynomial, weak (no IBP), m=8  | 14/0 ✗ | 0.2742  | 2.2      | ✗        |
| WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 59/0 ✗ | 0.6308  | 1.1      | ✗        |
| PDE-FIND + STRidge             | spectral, weak (no IBP), m=4    | 17/2 ✗ | 0.4508  | 2.3      | ✗        |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=4  | 0/5 ✗  | 0.4990  | 0.1      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=4  | 0/5 ✗  | 0.4990  | 0.1      | ✗        |
| PolyJet + RFE                  | polynomial, weak (no IBP), m=8  | 0/6 ✗  | 0.7758  | 0.2      | ✗        |
| PDE-FIND + RFE                 | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 1.7      | ✗        |
| WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.2      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.7182  | 1.0      | ✗        |

### 5% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Modified Multigraphs + MSTLS   | polynomial, weak (no IBP), m=4  | 9/0 ✗  | 0.3619  | 1.2      | ✗        |
| Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=4  | 11/0 ✗ | 0.3608  | 1.4      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=4  | 19/0 ✗ | 0.3609  | 1.1      | ✗        |
| PolyJet + STRidge              | spectral, weak (no IBP), m=4    | 27/0 ✗ | 0.4451  | 5.5      | ✗        |
| WSINDy + STRidge               | spectral, weak (no IBP), m=12   | 59/0 ✗ | 0.6317  | 1.1      | ✗        |
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=4  | 8/2 ✗  | 0.3481  | 7.8      | ✗        |
| PDE-FIND + STRidge             | spectral, weak (no IBP), m=12   | 16/2 ✗ | 0.5161  | 1.0      | ✗        |
| Multigraphs + MSTLS            | polynomial, weak (IBP), m=8     | 9/3 ✗  | 0.4782  | 0.5      | ✗        |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.6626  | 0.1      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.6626  | 0.2      | ✗        |
| PDE-FIND + MSTLS               | polynomial, weak (no IBP), m=12 | 1/7 ✗  | 0.4515  | 1.8      | ✗        |
| PDE-FIND + RFE                 | polynomial, weak (no IBP), m=12 | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| PolyJet + RFE                  | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.2      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.7208  | 0.9      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.001   | 49/2  | 54       | 0.4502  | 0.4382 | 9.5      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 1/7   | 1        | 0.4515  | 0.4404 | 0.3      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.01    | 43/2  | 48       | 0.4505  | 0.6593 | 5.2      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.001   | 49/2  | 54       | 0.4516  | 0.8832 | 4.6      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 84/0  | 43       | 0.4503  | 0.4387 | 3.3      | 🔴      |
| WSINDy           | Clean | RFE       |         | 1/7   | 1        | 0.4515  | 0.4404 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 86/0  | 42       | 0.4506  | 0.6594 | 3.2      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 88/0  | 43       | 0.4517  | 0.8833 | 4.0      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0058  | 0.0197 | 2.3      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0056  | 0.0197 | 0.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.1     | 47/0  | 54       | 0.4212  | 0.6518 | 12.9     | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.4      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 77/0  | 84       | 0.4500  | 0.8824 | 9.6      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0057  | 0.0197 | 0.6      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0055  | 0.0197 | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 17/0  | 11       | 0.4234  | 0.6535 | 1.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.1     | 24/0  | 14       | 0.4506  | 0.8838 | 1.6      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0057  | 0.0197 | 0.8      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0055  | 0.0197 | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 36/0  | 12       | 0.4236  | 0.6535 | 1.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.01    | 30/0  | 16       | 0.4513  | 0.8838 | 1.7      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.1      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 1.0     | 45/2  | 50       | 0.4503  | 0.4380 | 4.6      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 1/7   | 1        | 0.4515  | 0.4404 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.001   | 55/2  | 60       | 0.4510  | 0.6592 | 5.8      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.001   | 55/2  | 60       | 0.4527  | 0.8830 | 5.3      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 84/0  | 41       | 0.4505  | 0.4387 | 4.1      | 🔴      |
| WSINDy           | Clean | RFE       |         | 1/7   | 1        | 0.4515  | 0.4404 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.001   | 97/0  | 45       | 0.4508  | 0.6594 | 3.1      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 94/0  | 45       | 0.4521  | 0.8832 | 3.6      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 10.0    | 2/0   | 9        | 0.2946  | 0.2119 | 2.6      | 🔴      |
| PolyJet          | Clean | RFE       |         | 0/6   | 1        | 0.7068  | 0.7076 | 0.2      | 🔴      |
| PolyJet          | 2%    | STRidge   | 0.01    | 62/0  | 69       | 0.2942  | 0.6127 | 8.6      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 54/0  | 61       | 0.4004  | 0.8785 | 7.7      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 1.0     | 7/0   | 5        | 0.2740  | 0.1194 | 0.5      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 0/5   | 1        | 0.4989  | 0.5044 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 31/0  | 16       | 0.3084  | 0.6179 | 1.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 37/0  | 20       | 0.4056  | 0.8801 | 1.5      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.1871  | 0.2729 | 0.7      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 7/0   | 4        | 0.2740  | 0.1194 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 0.01    | 32/0  | 18       | 0.3083  | 0.6178 | 1.8      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.01    | 37/0  | 19       | 0.4056  | 0.8801 | 2.9      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.2      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 10.0    | 12 | no  | 19/2  | 24       | 2.3130  | 0.2164 | 1.6      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 2.8      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 0.3      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 10.0    | 12 | yes | 25/2  | 30       | 2.5018  | 0.2151 | 2.2      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2352 | 3.1      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2352 | 0.2      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 4  | no  | 17/2  | 22       | 0.4508  | 0.4182 | 2.3      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2431 | 1.7      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2431 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | yes | 17/2  | 22       | 0.5180  | 0.2394 | 1.3      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2431 | 1.8      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2431 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 16/2  | 21       | 0.5161  | 0.2758 | 1.0      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 1.4      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 16/2  | 21       | 0.5162  | 0.2758 | 1.0      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 1.4      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | no  | 70/0  | 36       | 0.7050  | 0.2284 | 1.1      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7177  | 0.6516 | 0.9      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 8  | yes | 67/0  | 35       | 0.4819  | 0.2963 | 1.6      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7177  | 0.6516 | 0.9      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2352 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 8  | no  | 88/0  | 43       | 0.4725  | 0.3092 | 1.8      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7182  | 0.6540 | 1.0      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2431 | 0.2      | 🔴      |
| WSINDy           | 2%    | STRidge   | 1.0     | 12 | yes | 93/0  | 41       | 0.7240  | 0.2364 | 1.3      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7182  | 0.6540 | 1.0      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2431 | 0.2      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 12 | no  | 59/0  | 34       | 0.6317  | 0.2733 | 1.1      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7208  | 0.6660 | 0.9      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 8  | yes | 88/0  | 43       | 0.4747  | 0.3659 | 1.8      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7208  | 0.6660 | 0.9      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 1.0     | 4  | no  | 0/0   | 7        | 0.0038  | 0.0152 | 2.2      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 7        | 0.0019  | 0.0095 | 5.2      | 🟢      |
| PolyJet          | Clean | RFE       |         | 8  | no  | 0/0   | 7        | 0.0019  | 0.0095 | 0.2      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0061  | 0.0063 | 1.4      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.5     | 8  | yes | 0/0   | 7        | 0.0836  | 0.0222 | 5.1      | 🟢      |
| PolyJet          | Clean | RFE       |         | 8  | yes | 0/6   | 1        | 0.8898  | 0.3298 | 0.2      | 🔴      |
| PolyJet          | 2%    | STRidge   | 0.1     | 4  | no  | 68/0  | 75       | 0.3650  | 0.3469 | 18.6     | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 8  | no  | 3/0   | 10       | 0.3483  | 0.2057 | 11.4     | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/6   | 1        | 0.8091  | 0.2326 | 0.3      | 🔴      |
| PolyJet          | 2%    | STRidge   | 1.0     | 8  | yes | 61/0  | 68       | 0.3246  | 0.1995 | 9.6      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 8  | yes | 4/0   | 11       | 0.3596  | 0.2022 | 11.1     | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 0.8091  | 0.2326 | 0.3      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 4  | no  | 27/0  | 34       | 0.4451  | 0.5352 | 5.5      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 3.9      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 8  | yes | 53/0  | 60       | 0.4781  | 0.3584 | 5.7      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 3.9      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.4      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.6      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 4  | no  | 0/0   | 3        | 0.0037  | 0.0152 | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 3        | 0.0817  | 0.0227 | 0.4      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 3        | 0.0057  | 0.0063 | 0.5      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 8  | yes | 0/5   | 1        | 0.5661  | 0.1266 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 4  | no  | 14/0  | 10       | 0.3910  | 0.3568 | 1.4      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 4  | no  | 5/0   | 5        | 0.3918  | 0.3575 | 0.8      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.7874  | 0.2197 | 0.2      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 8  | yes | 18/0  | 12       | 0.4259  | 0.2146 | 0.7      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 8  | yes | 4/0   | 4        | 0.4356  | 0.2151 | 0.7      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.7874  | 0.2197 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.1     | 4  | no  | 30/0  | 17       | 0.4511  | 0.5367 | 1.8      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.4      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.9576  | 0.3968 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.1     | 8  | yes | 32/0  | 18       | 0.4804  | 0.3616 | 0.7      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.4      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 0/5   | 1        | 0.9576  | 0.3968 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.6      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.8      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 4  | no  | 0/0   | 3        | 0.0037  | 0.0152 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 3        | 0.0817  | 0.0227 | 0.6      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 3        | 0.0057  | 0.0063 | 0.6      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 8  | yes | 0/5   | 1        | 0.5661  | 0.1266 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | no  | 17/0  | 12       | 0.4487  | 0.1535 | 0.5      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 4  | no  | 8/0   | 6        | 0.3896  | 0.3575 | 0.7      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.7874  | 0.2197 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 8  | yes | 37/0  | 19       | 0.7057  | 0.2449 | 0.7      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 8  | yes | 9/0   | 7        | 0.7548  | 0.2450 | 0.6      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.7874  | 0.2197 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 8  | no  | 35/0  | 17       | 0.4839  | 0.3616 | 0.7      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 8  | yes | 32/0  | 17       | 0.4889  | 0.3640 | 0.7      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 10.0    | 12 | no  | 19/2  | 24       | 2.1454  | 0.2182 | 1.0      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 1.8      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 0.1      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 10.0    | 12 | yes | 19/2  | 24       | 2.1460  | 0.2182 | 1.0      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2352 | 1.8      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2352 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | no  | 34/2  | 39       | 1.6831  | 0.2292 | 1.8      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2431 | 1.8      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2431 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 8  | yes | 25/2  | 30       | 0.5055  | 0.3079 | 2.5      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2431 | 1.8      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2431 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.1     | 8  | no  | 27/2  | 32       | 0.4599  | 0.3656 | 2.7      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 1.8      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 28/2  | 33       | 0.7093  | 0.2716 | 1.3      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 1.7      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | no  | 70/0  | 36       | 0.7075  | 0.2284 | 1.3      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7177  | 0.6516 | 1.0      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 8  | yes | 92/0  | 44       | 0.4962  | 0.2962 | 2.0      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7177  | 0.6516 | 1.1      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2352 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | no  | 59/0  | 34       | 0.6308  | 0.2366 | 1.1      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7182  | 0.6540 | 0.9      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2431 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 8  | yes | 95/0  | 45       | 0.4883  | 0.3091 | 1.8      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7182  | 0.6540 | 0.9      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2431 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 12 | no  | 59/0  | 34       | 0.6333  | 0.2733 | 1.1      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7208  | 0.6660 | 1.0      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 8  | yes | 83/0  | 43       | 0.4714  | 0.3660 | 1.9      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7208  | 0.6660 | 1.0      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 10.0    | 4  | no  | 0/0   | 7        | 0.1376  | 0.0719 | 2.7      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 7        | 0.0523  | 0.0110 | 5.9      | 🟢      |
| PolyJet          | Clean | RFE       |         | 8  | no  | 0/0   | 7        | 0.0523  | 0.0110 | 0.2      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 8  | yes | 3/0   | 10       | 0.0641  | 0.0085 | 2.0      | 🔴      |
| PolyJet          | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 7        | 0.0408  | 0.0070 | 4.2      | 🟢      |
| PolyJet          | Clean | RFE       |         | 8  | yes | 0/6   | 1        | 0.8976  | 0.3287 | 0.2      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 8  | no  | 14/0  | 21       | 0.2742  | 0.1447 | 2.2      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 4  | no  | 8/0   | 15       | 0.2652  | 0.2721 | 7.9      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 8  | no  | 0/6   | 1        | 0.7758  | 0.3243 | 0.2      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 8  | yes | 18/0  | 25       | 0.2718  | 0.1439 | 2.3      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 8  | yes | 9/0   | 16       | 0.2786  | 0.1446 | 5.7      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 8  | yes | 0/6   | 1        | 0.7758  | 0.3243 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 4  | no  | 65/0  | 72       | 0.3428  | 0.4983 | 10.3     | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 4  | no  | 8/2   | 13       | 0.3481  | 0.4999 | 7.8      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.3      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 8  | yes | 62/0  | 69       | 0.3252  | 0.3012 | 6.5      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 8  | yes | 9/2   | 14       | 0.4180  | 0.3036 | 6.3      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0424  | 0.0071 | 0.3      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.75    | 12 | no  | 0/0   | 3        | 0.0423  | 0.0071 | 0.4      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 4  | no  | 0/5   | 1        | 0.5437  | 0.2638 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 3        | 0.0626  | 0.0241 | 0.4      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 3        | 0.0409  | 0.0070 | 0.4      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 8  | yes | 0/5   | 1        | 0.5721  | 0.1208 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 4  | no  | 11/0  | 9        | 0.2188  | 0.2789 | 0.7      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 8  | no  | 5/0   | 5        | 0.1016  | 0.1483 | 0.5      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 4  | no  | 0/5   | 1        | 0.4990  | 0.3563 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | yes | 26/0  | 15       | 0.2288  | 0.1020 | 0.4      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 8  | yes | 5/0   | 5        | 0.1215  | 0.1487 | 0.5      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 8  | yes | 0/5   | 1        | 0.5006  | 0.1912 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 4  | no  | 19/0  | 12       | 0.3609  | 0.5045 | 1.1      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 8  | no  | 9/3   | 6        | 0.4792  | 0.3114 | 0.5      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.6626  | 0.2652 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 12 | yes | 35/0  | 19       | 0.3603  | 0.2233 | 0.5      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 8  | yes | 9/3   | 6        | 0.4782  | 0.3115 | 0.5      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.6626  | 0.2652 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0481  | 0.0115 | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.75    | 12 | no  | 0/0   | 3        | 0.0423  | 0.0071 | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 4  | no  | 0/5   | 1        | 0.5437  | 0.2638 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 3        | 0.0626  | 0.0241 | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 3        | 0.0409  | 0.0070 | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 8  | yes | 0/5   | 1        | 0.5721  | 0.1208 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 4  | no  | 11/0  | 8        | 0.2190  | 0.2790 | 0.8      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 3/0   | 4        | 0.3515  | 0.1548 | 0.6      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 4  | no  | 0/5   | 1        | 0.4990  | 0.3563 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 8  | yes | 11/0  | 9        | 0.6071  | 0.1473 | 0.6      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 8  | yes | 3/0   | 4        | 0.3170  | 0.1552 | 0.6      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 8  | yes | 0/5   | 1        | 0.5006  | 0.1912 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 4  | no  | 11/0  | 9        | 0.3608  | 0.5046 | 1.4      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 4  | no  | 9/0   | 6        | 0.3619  | 0.5053 | 1.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.6626  | 0.2652 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 12 | yes | 25/0  | 16       | 0.9582  | 0.2327 | 0.9      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.5     | 8  | yes | 4/5   | 3        | 0.4798  | 0.3256 | 0.9      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.6626  | 0.2652 | 0.2      | 🔴      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ)     |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:------------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0156       | 37240999.1  |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0096       | 103883766.2 |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0066       | 317175677.6 |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.0891       | 39715936.0  |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0172       | 107952857.1 |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0115       | 320266523.4 |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0548       | 105059979.6 |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0069       | 317477220.5 |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.0481       | 103982777.7 |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0101       | 317959968.0 |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 3.4737       | 39795395.5  |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 3.3691       | 90842207.1  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 3.4010       | 273651905.3 |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 0.4838       | 35806846.8  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.3048       | 100771904.4 |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.2790       | 293455352.0 |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 3.4168       | 91745285.8  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 3.4030       | 273959604.0 |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 0.3460       | 91260093.3  |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.2807       | 291754076.4 |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 19.5646      | 39785550.0  |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 20.5720      | 58952277.8  |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 21.0208      | 124183762.7 |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 1.7861       | 35711721.6  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 1.6500       | 72025688.0  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 1.6411       | 153225710.3 |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 20.6187      | 59009394.7  |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 21.0227      | 124299300.8 |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 1.6947       | 61871523.2  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 1.6430       | 152400614.3 |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ks_chaos`.

| Source                                        | Wall time | CPU | Threads |
|:----------------------------------------------|----------:|:----|:-------:|
| `sweep_paper_v2_ks_chaos_20260804_223951.csv` | ? min     | ?   | ?       |
