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
| Clean | pointwise     | spectral   | 🟩        | 0/0 ✓             | 0.0055     | multigraph/RFE         | 0.0197    | PolyJet/RFE     | 0/0 ✓     | 0.0056      | 0.0197     | 0.0205       |
| Clean | pointwise     | polynomial | 🟢        | 0/0 ✓             | 0.1421     | mod_multigraph/STRidge | 0.2729    | —               | 2/0 ✗     | 0.2174      | 0.2119     | 0.3112       |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0017     | multigraph/STRidge     | 0.0095    | PolyJet/MSTLS   | 0/0 ✓     | 0.0019      | 0.0095     | 0.0066       |
| Clean | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0350     | multigraph/MSTLS       | 0.0071    | PolyJet/RFE     | 0/0 ✓     | 0.0392      | 0.0110     | 0.0115       |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0057     | multigraph/MSTLS       | 0.0063    | PolyJet/STRidge | 0/0 ✓     | 0.0061      | 0.0063     | 0.0069       |
| Clean | weak (IBP)    | polynomial | 🟥        | 0/0 ✓             | 0.0340     | mod_multigraph/MSTLS   | 0.0070    | PolyJet/MSTLS   | 0/0 ✓     | 0.0324      | 0.0070     | 0.0101       |
| 2%    | pointwise     | spectral   | —        | 17/0 ✗            | 0.4234     | —                      | 0.6535    | —               | 47/0 ✗    | 0.4212      | 0.6518     | 85.8534      |
| 2%    | pointwise     | polynomial | —        | 31/0 ✗            | 0.3559     | —                      | 0.6179    | —               | 62/0 ✗    | 0.3628      | 0.6127     | 1.8548       |
| 2%    | weak (no IBP) | spectral   | —        | 5/0 ✗             | 0.3918     | —                      | 0.3575    | —               | 3/0 ✗     | 0.3483      | 0.2057     | 3.3691       |
| 2%    | weak (no IBP) | polynomial | —        | 3/0 ✗             | 0.2723     | —                      | 0.1548    | —               | 8/0 ✗     | 0.3319      | 0.2721     | 0.2790       |
| 2%    | weak (IBP)    | spectral   | —        | 4/0 ✗             | 0.4356     | —                      | 0.2151    | —               | 4/0 ✗     | 0.3596      | 0.2022     | 3.4030       |
| 2%    | weak (IBP)    | polynomial | —        | 3/0 ✗             | 0.2359     | —                      | 0.1552    | —               | 9/0 ✗     | 0.2710      | 0.1446     | 0.2807       |
| 5%    | pointwise     | spectral   | —        | 24/0 ✗            | 0.4506     | —                      | 0.8838    | —               | 77/0 ✗    | 0.4500      | 0.8824     | 133.9474     |
| 5%    | pointwise     | polynomial | —        | 26/0 ✗            | 0.4593     | —                      | 0.8801    | —               | 54/0 ✗    | 0.4604      | 0.8785     | 3.0277       |
| 5%    | weak (no IBP) | spectral   | —        | 30/0 ✗            | 0.4511     | —                      | 0.5367    | —               | 27/0 ✗    | 0.4451      | 0.5352     | 19.5646      |
| 5%    | weak (no IBP) | polynomial | —        | 9/0 ✗             | 0.3606     | —                      | 0.5053    | —               | 59/0 ✗    | 0.5998      | 0.2733     | 1.6411       |
| 5%    | weak (IBP)    | spectral   | —        | 32/0 ✗            | 0.4804     | —                      | 0.3616    | —               | 53/0 ✗    | 0.4781      | 0.3584     | 20.6187      |
| 5%    | weak (IBP)    | polynomial | —        | 25/0 ✗            | 0.9377     | —                      | 0.2327    | —               | 62/0 ✗    | 0.3036      | 0.3012     | 1.6430       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                    | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + STRidge          | spectral, weak (no IBP), m=8  | 0/0 ✓  | 0.0017  | 0.4      | —        |
| Modified Multigraphs + STRidge | spectral, weak (no IBP), m=8  | 0/0 ✓  | 0.0017  | 0.5      | +0.0%    |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=8  | 0/0 ✓  | 0.0017  | 0.6      | +1.0%    |
| Modified Multigraphs + MSTLS   | spectral, weak (no IBP), m=8  | 0/0 ✓  | 0.0017  | 0.6      | +1.0%    |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=8  | 0/0 ✓  | 0.0019  | 5.2      | +14.5%   |
| PolyJet + RFE                  | spectral, weak (no IBP), m=8  | 0/0 ✓  | 0.0019  | 0.2      | +14.5%   |
| Multigraphs + RFE              | spectral, weak (no IBP), m=4  | 0/0 ✓  | 0.0037  | 0.1      | +121.5%  |
| Modified Multigraphs + RFE     | spectral, weak (no IBP), m=4  | 0/0 ✓  | 0.0037  | 0.3      | +121.5%  |
| PolyJet + STRidge              | spectral, weak (no IBP), m=4  | 0/0 ✓  | 0.0038  | 2.7      | +126.3%  |
| WSINDy + STRidge               | spectral, weak (IBP), m=8     | 67/0 ✗ | 0.4819  | 1.5      | ✗        |
| PDE-FIND + STRidge             | polynomial, weak (IBP), m=12  | 19/2 ✗ | 1.7161  | 1.1      | ✗        |
| PDE-FIND + RFE                 | spectral, weak (no IBP), m=12 | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=12 | 1/7 ✗  | 0.4515  | 1.8      | ✗        |
| WSINDy + RFE                   | spectral, weak (no IBP), m=12 | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=12 | 1/7 ✗  | 0.7177  | 0.9      | ✗        |

### 2% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=8     | 3/0 ✗  | 0.2359  | 0.5      | ✗        |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=8    | 3/0 ✗  | 0.3483  | 6.0      | ✗        |
| Multigraphs + MSTLS            | spectral, weak (IBP), m=8       | 4/0 ✗  | 0.4356  | 0.6      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=4  | 11/0 ✗ | 0.2738  | 1.3      | ✗        |
| Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=4  | 11/0 ✗ | 0.2761  | 0.7      | ✗        |
| PolyJet + STRidge              | polynomial, weak (no IBP), m=8  | 14/0 ✗ | 0.2636  | 2.3      | ✗        |
| WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 59/0 ✗ | 0.6002  | 1.3      | ✗        |
| PDE-FIND + STRidge             | spectral, weak (no IBP), m=4    | 17/2 ✗ | 0.4508  | 2.2      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.3415  | 0.1      | ✗        |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=8  | 0/5 ✗  | 0.3436  | 0.3      | ✗        |
| PolyJet + RFE                  | polynomial, weak (no IBP), m=8  | 0/6 ✗  | 0.6793  | 0.2      | ✗        |
| PDE-FIND + RFE                 | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 1.5      | ✗        |
| WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.7182  | 1.1      | ✗        |

### 5% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Modified Multigraphs + MSTLS   | polynomial, weak (no IBP), m=4  | 9/0 ✗  | 0.3606  | 0.9      | ✗        |
| Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=4  | 11/0 ✗ | 0.3568  | 1.1      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=4  | 19/0 ✗ | 0.3569  | 1.2      | ✗        |
| PolyJet + STRidge              | spectral, weak (no IBP), m=4    | 27/0 ✗ | 0.4451  | 5.4      | ✗        |
| WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 59/0 ✗ | 0.5998  | 1.2      | ✗        |
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=4  | 8/2 ✗  | 0.3742  | 7.1      | ✗        |
| PDE-FIND + STRidge             | spectral, weak (no IBP), m=12   | 16/2 ✗ | 0.5161  | 1.1      | ✗        |
| Multigraphs + MSTLS            | polynomial, weak (IBP), m=8     | 9/3 ✗  | 0.3348  | 0.6      | ✗        |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.6270  | 0.2      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.6270  | 0.2      | ✗        |
| PolyJet + RFE                  | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| PDE-FIND + RFE                 | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 1.4      | ✗        |
| WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.7208  | 0.9      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.001   | 49/2  | 54       | 0.4502  | 0.4382 | 5.6      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 1/7   | 1        | 0.4515  | 0.4404 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.01    | 43/2  | 48       | 0.4505  | 0.6593 | 4.8      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.001   | 49/2  | 54       | 0.4516  | 0.8832 | 4.7      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 84/0  | 43       | 0.4503  | 0.4387 | 3.2      | 🔴      |
| WSINDy           | Clean | RFE       |         | 1/7   | 1        | 0.4515  | 0.4404 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 86/0  | 42       | 0.4506  | 0.6594 | 3.8      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 88/0  | 43       | 0.4517  | 0.8833 | 3.1      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0058  | 0.0197 | 2.9      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0056  | 0.0197 | 0.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.1     | 47/0  | 54       | 0.4212  | 0.6518 | 6.9      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 77/0  | 84       | 0.4500  | 0.8824 | 9.0      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0057  | 0.0197 | 0.6      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0055  | 0.0197 | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 17/0  | 11       | 0.4234  | 0.6535 | 1.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.1     | 24/0  | 14       | 0.4506  | 0.8838 | 1.8      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0057  | 0.0197 | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0055  | 0.0197 | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 36/0  | 12       | 0.4236  | 0.6535 | 1.1      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.01    | 30/0  | 16       | 0.4513  | 0.8838 | 1.5      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.2      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 1.0     | 45/2  | 50       | 0.5217  | 0.4380 | 4.6      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 1/7   | 1        | 0.5227  | 0.4404 | 0.3      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.001   | 55/2  | 60       | 0.5226  | 0.6592 | 5.7      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.6602 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.001   | 51/2  | 56       | 0.5235  | 0.8831 | 4.7      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.8839 | 0.2      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 84/0  | 41       | 0.5211  | 0.4387 | 3.3      | 🔴      |
| WSINDy           | Clean | RFE       |         | 1/7   | 1        | 0.5227  | 0.4404 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.001   | 97/0  | 45       | 0.5225  | 0.6594 | 3.6      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.6602 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 94/0  | 45       | 0.5233  | 0.8832 | 3.8      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.8839 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 10.0    | 2/0   | 9        | 0.2174  | 0.2119 | 2.7      | 🔴      |
| PolyJet          | Clean | RFE       |         | 0/6   | 1        | 0.6498  | 0.7076 | 0.2      | 🔴      |
| PolyJet          | 2%    | STRidge   | 0.01    | 62/0  | 69       | 0.3628  | 0.6127 | 8.4      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.6602 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 54/0  | 61       | 0.4604  | 0.8785 | 7.8      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.8839 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 1.0     | 7/0   | 5        | 0.2750  | 0.1194 | 0.5      | 🔴      |
| Multigraphs      | Clean | RFE       |         | 0/5   | 1        | 0.3549  | 0.5044 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 31/0  | 16       | 0.3559  | 0.6179 | 2.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.6602 | 0.3      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.1     | 26/0  | 15       | 0.4593  | 0.8801 | 1.6      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.8839 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.1421  | 0.2729 | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 7/0   | 4        | 0.2755  | 0.1194 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 0.01    | 32/0  | 18       | 0.3560  | 0.6178 | 1.7      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.6602 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.01    | 37/0  | 19       | 0.4605  | 0.8801 | 2.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.8839 | 0.2      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 10.0    | 12 | no  | 19/2  | 24       | 2.3130  | 0.2164 | 1.0      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 1.8      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 0.1      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 10.0    | 12 | yes | 25/2  | 30       | 2.5018  | 0.2151 | 1.2      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2352 | 1.9      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2352 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 4  | no  | 17/2  | 22       | 0.4508  | 0.4182 | 2.2      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2431 | 1.5      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2431 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | yes | 17/2  | 22       | 0.5180  | 0.2394 | 1.1      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2431 | 1.5      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2431 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 16/2  | 21       | 0.5161  | 0.2758 | 1.1      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 1.4      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 16/2  | 21       | 0.5162  | 0.2758 | 1.1      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 1.4      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | no  | 70/0  | 36       | 0.7050  | 0.2284 | 1.1      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7177  | 0.6516 | 0.9      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 8  | yes | 67/0  | 35       | 0.4819  | 0.2963 | 1.5      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7177  | 0.6516 | 0.9      | 🔴      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2352 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.01    | 8  | no  | 88/0  | 43       | 0.4725  | 0.3092 | 1.9      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7182  | 0.6540 | 1.1      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2431 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 1.0     | 12 | yes | 93/0  | 41       | 0.7240  | 0.2364 | 1.4      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7182  | 0.6540 | 1.1      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2431 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 12 | no  | 59/0  | 34       | 0.6317  | 0.2733 | 1.1      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7208  | 0.6660 | 0.9      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 8  | yes | 88/0  | 43       | 0.4747  | 0.3659 | 1.7      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7208  | 0.6660 | 0.9      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 1.0     | 4  | no  | 0/0   | 7        | 0.0038  | 0.0152 | 2.7      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 7        | 0.0019  | 0.0095 | 5.2      | 🟢      |
| PolyJet          | Clean | RFE       |         | 8  | no  | 0/0   | 7        | 0.0019  | 0.0095 | 0.2      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 7        | 0.0061  | 0.0063 | 1.3      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.5     | 8  | yes | 0/0   | 7        | 0.0836  | 0.0222 | 5.2      | 🟢      |
| PolyJet          | Clean | RFE       |         | 8  | yes | 0/6   | 1        | 0.8898  | 0.3298 | 0.2      | 🔴      |
| PolyJet          | 2%    | STRidge   | 0.1     | 8  | no  | 54/0  | 61       | 0.3142  | 0.1996 | 6.1      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 8  | no  | 3/0   | 10       | 0.3483  | 0.2057 | 6.0      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/6   | 1        | 0.8091  | 0.2326 | 0.3      | 🔴      |
| PolyJet          | 2%    | STRidge   | 0.1     | 12 | yes | 72/0  | 79       | 0.4007  | 0.1399 | 4.7      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 8  | yes | 4/0   | 11       | 0.3596  | 0.2022 | 5.9      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 0.8091  | 0.2326 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 4  | no  | 27/0  | 34       | 0.4451  | 0.5352 | 5.4      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 3.7      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 8  | yes | 53/0  | 60       | 0.4781  | 0.3584 | 5.2      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 3.7      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.4      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.6      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 4  | no  | 0/0   | 3        | 0.0037  | 0.0152 | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 3        | 0.0817  | 0.0227 | 0.4      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 3        | 0.0057  | 0.0063 | 0.5      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 8  | yes | 0/5   | 1        | 0.5661  | 0.1266 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 4  | no  | 14/0  | 10       | 0.3910  | 0.3568 | 1.3      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 4  | no  | 5/0   | 5        | 0.3918  | 0.3575 | 0.7      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.7874  | 0.2197 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 8  | yes | 18/0  | 12       | 0.4259  | 0.2146 | 0.6      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 8  | yes | 4/0   | 4        | 0.4356  | 0.2151 | 0.6      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.7874  | 0.2197 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.1     | 4  | no  | 30/0  | 17       | 0.4511  | 0.5367 | 1.9      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.4      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.9576  | 0.3968 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.1     | 8  | yes | 32/0  | 18       | 0.4804  | 0.3616 | 0.7      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.6      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 0/5   | 1        | 0.9576  | 0.3968 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.6      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 4  | no  | 0/0   | 3        | 0.0037  | 0.0152 | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 3        | 0.0817  | 0.0227 | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 3        | 0.0057  | 0.0063 | 0.5      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 8  | yes | 0/5   | 1        | 0.5661  | 0.1266 | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | no  | 17/0  | 12       | 0.4487  | 0.1535 | 0.5      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 4  | no  | 8/0   | 6        | 0.3896  | 0.3575 | 0.7      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.7874  | 0.2197 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 8  | yes | 37/0  | 19       | 0.7057  | 0.2449 | 0.7      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 8  | yes | 9/0   | 7        | 0.7548  | 0.2450 | 0.6      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.7874  | 0.2197 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 4  | no  | 35/0  | 17       | 0.4510  | 0.5367 | 1.5      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.5      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 8  | yes | 32/0  | 17       | 0.4889  | 0.3640 | 0.8      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2788 | 0.2      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 10.0    | 12 | no  | 19/2  | 24       | 1.7170  | 0.2182 | 1.1      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.5226  | 0.3865 | 3.3      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.3865 | 0.2      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 10.0    | 12 | yes | 19/2  | 24       | 1.7161  | 0.2182 | 1.1      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.5226  | 0.3001 | 2.8      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5226  | 0.3001 | 0.2      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | no  | 34/2  | 39       | 1.3403  | 0.2292 | 1.8      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.5226  | 0.4199 | 3.0      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.4199 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 8  | yes | 25/2  | 30       | 0.5480  | 0.3079 | 2.5      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.5226  | 0.3127 | 2.4      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5226  | 0.3127 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 28/2  | 33       | 0.6529  | 0.2716 | 1.3      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.5226  | 0.5403 | 2.7      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.5403 | 0.2      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 28/2  | 33       | 0.6521  | 0.2716 | 1.4      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.5227  | 0.3688 | 1.9      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5227  | 0.3688 | 0.2      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | no  | 70/0  | 36       | 0.6364  | 0.2284 | 1.4      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7444  | 0.6516 | 1.1      | 🔴      |
| WSINDy           | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.3865 | 0.2      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 8  | yes | 92/0  | 44       | 0.5345  | 0.2962 | 2.0      | 🔴      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7444  | 0.6516 | 1.1      | 🔴      |
| WSINDy           | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5226  | 0.3001 | 0.2      | 🔴      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | no  | 59/0  | 34       | 0.6002  | 0.2366 | 1.3      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7449  | 0.6540 | 1.0      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.4199 | 0.1      | 🔴      |
| WSINDy           | 2%    | STRidge   | 1.0     | 12 | yes | 95/0  | 43       | 0.6506  | 0.2364 | 1.4      | 🔴      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7449  | 0.6540 | 1.0      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5226  | 0.3127 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 12 | no  | 59/0  | 34       | 0.5998  | 0.2733 | 1.2      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7472  | 0.6660 | 1.0      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.5403 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 8  | yes | 83/0  | 43       | 0.5246  | 0.3660 | 1.9      | 🔴      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7472  | 0.6660 | 1.1      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5227  | 0.3688 | 0.1      | 🔴      |
| PolyJet          | Clean | STRidge   | 10.0    | 4  | no  | 0/0   | 7        | 0.0884  | 0.0719 | 2.7      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 7        | 0.0392  | 0.0110 | 4.9      | 🟢      |
| PolyJet          | Clean | RFE       |         | 8  | no  | 0/0   | 7        | 0.0392  | 0.0110 | 0.2      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 8  | yes | 3/0   | 10       | 0.0551  | 0.0085 | 1.9      | 🔴      |
| PolyJet          | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 7        | 0.0324  | 0.0070 | 3.6      | 🟢      |
| PolyJet          | Clean | RFE       |         | 8  | yes | 0/6   | 1        | 0.7813  | 0.3287 | 0.2      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 8  | no  | 14/0  | 21       | 0.2636  | 0.1447 | 2.3      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 4  | no  | 8/0   | 15       | 0.3319  | 0.2721 | 8.0      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 8  | no  | 0/6   | 1        | 0.6793  | 0.3243 | 0.2      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 8  | yes | 18/0  | 25       | 0.2619  | 0.1439 | 2.5      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 8  | yes | 9/0   | 16       | 0.2710  | 0.1446 | 5.5      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 8  | yes | 0/6   | 1        | 0.6793  | 0.3243 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 4  | no  | 65/0  | 72       | 0.3700  | 0.4983 | 9.2      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 4  | no  | 8/2   | 13       | 0.3742  | 0.4999 | 7.1      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.5403 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 8  | yes | 62/0  | 69       | 0.3036  | 0.3012 | 5.6      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 8  | yes | 9/2   | 14       | 0.3521  | 0.3036 | 5.7      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5227  | 0.3688 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0351  | 0.0071 | 0.4      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.75    | 12 | no  | 0/0   | 3        | 0.0350  | 0.0071 | 0.4      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 4  | no  | 0/5   | 1        | 0.3492  | 0.2638 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 3        | 0.0453  | 0.0241 | 0.4      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 3        | 0.0340  | 0.0070 | 0.4      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 8  | yes | 0/5   | 1        | 0.3706  | 0.1208 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 4  | no  | 11/0  | 9        | 0.2738  | 0.2789 | 1.3      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 8  | no  | 5/0   | 5        | 0.1162  | 0.1483 | 1.0      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.3436  | 0.1912 | 0.3      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | yes | 26/0  | 15       | 0.1926  | 0.1020 | 0.7      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 8  | yes | 5/0   | 5        | 0.1228  | 0.1487 | 0.9      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 8  | yes | 0/5   | 1        | 0.3436  | 0.1912 | 0.3      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 4  | no  | 19/0  | 12       | 0.3569  | 0.5045 | 1.2      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 8  | no  | 9/3   | 6        | 0.3354  | 0.3114 | 0.5      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.6270  | 0.2652 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 8  | yes | 37/0  | 20       | 0.2853  | 0.3092 | 0.8      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 8  | yes | 9/3   | 6        | 0.3348  | 0.3115 | 0.6      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.6270  | 0.2652 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0386  | 0.0115 | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.75    | 12 | no  | 0/0   | 3        | 0.0350  | 0.0071 | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 4  | no  | 0/5   | 1        | 0.3492  | 0.2638 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 3        | 0.0453  | 0.0241 | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.75    | 12 | yes | 0/0   | 3        | 0.0340  | 0.0070 | 0.4      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 8  | yes | 0/5   | 1        | 0.3706  | 0.1208 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 4  | no  | 11/0  | 8        | 0.2761  | 0.2790 | 0.7      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 3/0   | 4        | 0.2723  | 0.1548 | 0.5      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.3415  | 0.1301 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 8  | yes | 11/0  | 9        | 0.5953  | 0.1473 | 0.5      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 8  | yes | 3/0   | 4        | 0.2359  | 0.1552 | 0.5      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.3415  | 0.1301 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 4  | no  | 11/0  | 9        | 0.3568  | 0.5046 | 1.1      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 4  | no  | 9/0   | 6        | 0.3606  | 0.5053 | 0.9      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.6270  | 0.2652 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 12 | yes | 25/0  | 16       | 0.9377  | 0.2327 | 0.7      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.5     | 8  | yes | 4/5   | 3        | 0.3434  | 0.3256 | 0.8      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.6270  | 0.2652 | 0.2      | 🔴      |


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
| 5%    | weak (no IBP) | polynomial | 12 | no  | 1.6411       | 153225710.5 |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 20.6187      | 59009394.7  |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 21.0227      | 124299300.8 |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 1.6947       | 61871523.2  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 1.6430       | 152400614.1 |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ks_chaos`.

| Source                                        | Wall time | CPU | Threads |
|:----------------------------------------------|----------:|:----|:-------:|
| `sweep_no_ac_ch_ks_chaos_20260803_170832.csv` | ? min     | ?   | ?       |
