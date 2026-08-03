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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best         | Mg E_data | Std best      | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-------------------|:----------|:--------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | 🟩        | 0/0 ✓             | 0.0055     | multigraph/RFE     | 0.0197    | PolyJet/RFE   | 0/0 ✓     | 0.0056      | 0.0197     | —            |
| Clean | pointwise     | polynomial | —        | 7/0 ✗             | 0.2750     | —                  | 0.1194    | —             | 2/0 ✗     | 0.2174      | 0.2119     | —            |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0017     | multigraph/STRidge | 0.0095    | PolyJet/MSTLS | 0/0 ✓     | 0.0019      | 0.0095     | —            |
| Clean | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0350     | multigraph/MSTLS   | 0.0071    | PolyJet/MSTLS | 0/0 ✓     | 0.0392      | 0.0110     | —            |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0020     | multigraph/MSTLS   | 0.0087    | PolyJet/MSTLS | 0/0 ✓     | 0.0024      | 0.0087     | —            |
| Clean | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0277     | multigraph/STRidge | 0.0156    | PolyJet/RFE   | 0/0 ✓     | 0.0288      | 0.0146     | —            |
| 2%    | pointwise     | spectral   | —        | 17/0 ✗            | 0.4234     | —                  | 0.6535    | —             | 47/0 ✗    | 0.4212      | 0.6518     | —            |
| 2%    | pointwise     | polynomial | —        | 31/0 ✗            | 0.3559     | —                  | 0.6179    | —             | 62/0 ✗    | 0.3628      | 0.6127     | —            |
| 2%    | weak (no IBP) | spectral   | —        | 5/0 ✗             | 0.3918     | —                  | 0.3575    | —             | 3/0 ✗     | 0.3483      | 0.2057     | —            |
| 2%    | weak (no IBP) | polynomial | —        | 5/0 ✗             | 0.1162     | —                  | 0.1483    | —             | 8/0 ✗     | 0.3319      | 0.2721     | —            |
| 2%    | weak (IBP)    | spectral   | —        | 4/0 ✗             | 0.4217     | —                  | 0.2039    | —             | 4/0 ✗     | 0.3380      | 0.1904     | —            |
| 2%    | weak (IBP)    | polynomial | —        | 5/0 ✗             | 0.1101     | —                  | 0.1395    | —             | 9/0 ✗     | 0.2682      | 0.1359     | —            |
| 5%    | pointwise     | spectral   | —        | 24/0 ✗            | 0.4506     | —                  | 0.8838    | —             | 77/0 ✗    | 0.4500      | 0.8824     | —            |
| 5%    | pointwise     | polynomial | —        | 26/0 ✗            | 0.4593     | —                  | 0.8801    | —             | 54/0 ✗    | 0.4604      | 0.8785     | —            |
| 5%    | weak (no IBP) | spectral   | —        | 30/0 ✗            | 0.4511     | —                  | 0.5367    | —             | 27/0 ✗    | 0.4451      | 0.5352     | —            |
| 5%    | weak (no IBP) | polynomial | —        | 19/0 ✗            | 0.3569     | —                  | 0.5045    | —             | 59/0 ✗    | 0.5998      | 0.2733     | —            |
| 5%    | weak (IBP)    | spectral   | —        | 32/0 ✗            | 0.4558     | —                  | 0.4681    | —             | 54/0 ✗    | 0.4502      | 0.4655     | —            |
| 5%    | weak (IBP)    | polynomial | —        | 32/0 ✗            | 0.3162     | —                  | 0.4245    | —             | 42/0 ✗    | 0.2997      | 0.2840     | —            |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + STRidge | spectral, weak (no IBP), m=8    | 0/0 ✓  | 0.0017  | 0.3      | —        |
| Multigraphs + MSTLS   | spectral, weak (no IBP), m=8    | 0/0 ✓  | 0.0017  | 0.5      | +1.0%    |
| PolyJet + MSTLS       | spectral, weak (no IBP), m=8    | 0/0 ✓  | 0.0019  | 6.7      | +14.5%   |
| PolyJet + RFE         | spectral, weak (no IBP), m=8    | 0/0 ✓  | 0.0019  | 0.2      | +14.5%   |
| Multigraphs + RFE     | spectral, weak (IBP), m=8       | 0/0 ✓  | 0.0028  | 0.1      | +66.2%   |
| PolyJet + STRidge     | spectral, weak (IBP), m=8       | 0/0 ✓  | 0.0029  | 2.0      | +77.5%   |
| WSINDy + STRidge      | polynomial, weak (no IBP), m=12 | 70/0 ✗ | 0.6364  | 1.2      | ✗        |
| PDE-FIND + STRidge    | polynomial, weak (no IBP), m=12 | 19/2 ✗ | 1.7170  | 1.0      | ✗        |
| PDE-FIND + RFE        | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 1.6      | ✗        |
| WSINDy + RFE          | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| WSINDy + MSTLS        | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.7177  | 0.9      | ✗        |

### 2% noise

| Method                | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| PolyJet + MSTLS       | spectral, weak (no IBP), m=8    | 3/0 ✗  | 0.3483  | 6.5      | ✗        |
| Multigraphs + MSTLS   | spectral, weak (IBP), m=12      | 4/0 ✗  | 0.4217  | 0.3      | ✗        |
| Multigraphs + STRidge | polynomial, weak (no IBP), m=4  | 11/0 ✗ | 0.2738  | 0.7      | ✗        |
| PolyJet + STRidge     | polynomial, weak (no IBP), m=8  | 14/0 ✗ | 0.2636  | 2.3      | ✗        |
| WSINDy + STRidge      | polynomial, weak (no IBP), m=12 | 59/0 ✗ | 0.6002  | 1.3      | ✗        |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=4    | 17/2 ✗ | 0.4508  | 2.3      | ✗        |
| Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/5 ✗  | 0.3425  | 0.1      | ✗        |
| PolyJet + RFE         | polynomial, weak (IBP), m=8     | 0/6 ✗  | 0.6527  | 0.2      | ✗        |
| PDE-FIND + RFE        | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 1.5      | ✗        |
| WSINDy + RFE          | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| WSINDy + MSTLS        | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.7182  | 1.0      | ✗        |

### 5% noise

| Method                | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + STRidge | polynomial, weak (no IBP), m=4  | 19/0 ✗ | 0.3569  | 0.7      | ✗        |
| PolyJet + STRidge     | spectral, weak (no IBP), m=4    | 27/0 ✗ | 0.4451  | 5.8      | ✗        |
| WSINDy + STRidge      | polynomial, weak (no IBP), m=12 | 59/0 ✗ | 0.5998  | 1.1      | ✗        |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=4  | 8/2 ✗  | 0.3742  | 9.8      | ✗        |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=4    | 19/2 ✗ | 0.4512  | 2.5      | ✗        |
| Multigraphs + MSTLS   | polynomial, weak (no IBP), m=8  | 9/3 ✗  | 0.3354  | 0.5      | ✗        |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.6270  | 0.1      | ✗        |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 1.6      | ✗        |
| PDE-FIND + RFE        | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| WSINDy + RFE          | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.1      | ✗        |
| PolyJet + RFE         | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.4515  | 0.2      | ✗        |
| WSINDy + MSTLS        | spectral, weak (no IBP), m=12   | 1/7 ✗  | 0.7208  | 1.0      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.001   | 49/2  | 54       | 0.4502  | 0.4382 | 4.9      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 1/7   | 1        | 0.4515  | 0.4404 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.01    | 43/2  | 48       | 0.4505  | 0.6593 | 4.3      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 49/2  | 54       | 0.4516  | 0.8832 | 4.7      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.001   | 84/0  | 43       | 0.4503  | 0.4387 | 2.9      | 🔴      |
| WSINDy      | Clean | RFE       |         | 1/7   | 1        | 0.4515  | 0.4404 | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.01    | 86/0  | 42       | 0.4506  | 0.6594 | 2.9      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 88/0  | 43       | 0.4517  | 0.8833 | 2.7      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.1     | 0/0   | 7        | 0.0058  | 0.0197 | 3.1      | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 7        | 0.0056  | 0.0197 | 0.2      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.1     | 47/0  | 54       | 0.4212  | 0.6518 | 7.4      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 77/0  | 84       | 0.4500  | 0.8824 | 9.4      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.2      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0057  | 0.0197 | 0.5      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0055  | 0.0197 | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.1     | 17/0  | 11       | 0.4234  | 0.6535 | 0.8      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.6602 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.1     | 24/0  | 14       | 0.4506  | 0.8838 | 0.8      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.4515  | 0.8839 | 0.1      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 1.0     | 45/2  | 50       | 0.5217  | 0.4380 | 4.0      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 1/7   | 1        | 0.5227  | 0.4404 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.001   | 52/2  | 57       | 0.5226  | 0.6592 | 4.8      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.6602 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 51/2  | 56       | 0.5235  | 0.8831 | 4.4      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.8839 | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.001   | 84/0  | 41       | 0.5211  | 0.4387 | 2.7      | 🔴      |
| WSINDy      | Clean | RFE       |         | 1/7   | 1        | 0.5227  | 0.4404 | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.001   | 97/0  | 45       | 0.5225  | 0.6594 | 3.0      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.6602 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.01    | 94/0  | 45       | 0.5233  | 0.8832 | 2.8      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.8839 | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 10.0    | 2/0   | 9        | 0.2174  | 0.2119 | 2.8      | 🔴      |
| PolyJet     | Clean | RFE       |         | 0/6   | 1        | 0.6498  | 0.7076 | 0.3      | 🔴      |
| PolyJet     | 2%    | STRidge   | 0.01    | 62/0  | 69       | 0.3628  | 0.6127 | 9.1      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.6602 | 0.3      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 54/0  | 61       | 0.4604  | 0.8785 | 7.8      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.8839 | 0.2      | 🔴      |
| Multigraphs | Clean | STRidge   | 1.0     | 7/0   | 5        | 0.2750  | 0.1194 | 0.5      | 🔴      |
| Multigraphs | Clean | RFE       |         | 0/5   | 1        | 0.3549  | 0.5044 | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 0.1     | 31/0  | 16       | 0.3559  | 0.6179 | 0.9      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.6602 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.1     | 26/0  | 15       | 0.4593  | 0.8801 | 0.8      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.5227  | 0.8839 | 0.1      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 10.0    | 12 | no  | 19/2  | 24       | 2.3130  | 0.2164 | 1.0      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 1.6      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 0.1      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 1.0     | 12 | yes | 50/2  | 55       | 0.5938  | 0.2807 | 2.1      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.2889 | 1.6      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2889 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 4  | no  | 17/2  | 22       | 0.4508  | 0.4182 | 2.3      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2431 | 1.5      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2431 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 8  | yes | 18/2  | 23       | 0.4511  | 0.3796 | 2.0      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.3004 | 1.5      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.3004 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 4  | no  | 19/2  | 24       | 0.4512  | 0.5391 | 2.5      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 1.6      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 8  | yes | 35/2  | 40       | 0.4515  | 0.4715 | 2.3      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.3519 | 1.6      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.3519 | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | no  | 70/0  | 36       | 0.7050  | 0.2284 | 1.2      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7177  | 0.6516 | 0.9      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2352 | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | yes | 87/0  | 43       | 0.4833  | 0.2838 | 1.3      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7216  | 0.6696 | 1.0      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.2889 | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.01    | 8  | no  | 88/0  | 43       | 0.4725  | 0.3092 | 2.0      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7182  | 0.6540 | 1.0      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2431 | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.01    | 8  | yes | 89/0  | 43       | 0.4526  | 0.3791 | 1.9      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7225  | 0.6738 | 1.0      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.3004 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 12 | no  | 59/0  | 34       | 0.6317  | 0.2733 | 1.2      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7208  | 0.6660 | 1.0      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 12 | yes | 83/0  | 42       | 0.4784  | 0.3481 | 1.3      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7273  | 0.6942 | 0.9      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.3519 | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 4  | no  | 0/0   | 7        | 0.0038  | 0.0152 | 2.6      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 7        | 0.0019  | 0.0095 | 6.7      | 🟢      |
| PolyJet     | Clean | RFE       |         | 8  | no  | 0/0   | 7        | 0.0019  | 0.0095 | 0.2      | 🟢      |
| PolyJet     | Clean | STRidge   | 1.0     | 8  | yes | 0/0   | 7        | 0.0029  | 0.0127 | 2.0      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 7        | 0.0024  | 0.0087 | 3.8      | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0024  | 0.0087 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.1     | 4  | no  | 68/0  | 75       | 0.3650  | 0.3469 | 10.4     | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 8  | no  | 3/0   | 10       | 0.3483  | 0.2057 | 6.5      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 0/6   | 1        | 0.8091  | 0.2326 | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | yes | 41/0  | 48       | 0.3188  | 0.1871 | 3.0      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 4/0   | 11       | 0.3380  | 0.1904 | 4.1      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 0.8110  | 0.3011 | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 4  | no  | 27/0  | 34       | 0.4451  | 0.5352 | 5.8      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 3.8      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 8  | yes | 54/0  | 61       | 0.4502  | 0.4655 | 4.8      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.3519 | 3.9      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.3519 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.3      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 3        | 0.0017  | 0.0095 | 0.5      | 🟢      |
| Multigraphs | Clean | RFE       |         | 4  | no  | 0/0   | 3        | 0.0037  | 0.0152 | 0.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0022  | 0.0087 | 0.3      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0020  | 0.0087 | 0.3      | 🟢      |
| Multigraphs | Clean | RFE       |         | 8  | yes | 0/0   | 3        | 0.0028  | 0.0127 | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 4  | no  | 14/0  | 10       | 0.3910  | 0.3568 | 0.7      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.5     | 4  | no  | 5/0   | 5        | 0.3918  | 0.3575 | 0.6      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.7874  | 0.2197 | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 10.0    | 8  | yes | 14/0  | 10       | 0.3935  | 0.2983 | 0.5      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.5     | 12 | yes | 4/0   | 4        | 0.4217  | 0.2039 | 0.3      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.7883  | 0.2785 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.1     | 4  | no  | 30/0  | 17       | 0.4511  | 0.5367 | 0.9      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.4515  | 0.2788 | 0.3      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.9576  | 0.3968 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.1     | 8  | yes | 32/0  | 18       | 0.4558  | 0.4681 | 0.6      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.4515  | 0.3519 | 0.3      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.4515  | 0.3519 | 0.1      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 10.0    | 12 | no  | 19/2  | 24       | 1.7170  | 0.2182 | 1.0      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.5226  | 0.3865 | 3.3      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.3865 | 0.1      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 0.1     | 12 | yes | 30/2  | 35       | 0.5699  | 0.2817 | 1.7      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.5226  | 0.3592 | 2.2      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5226  | 0.3592 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | no  | 34/2  | 39       | 1.3403  | 0.2292 | 1.8      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.5226  | 0.4199 | 3.2      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.4199 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 8  | yes | 51/2  | 56       | 0.5211  | 0.3784 | 3.2      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.5226  | 0.3821 | 2.3      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5226  | 0.3821 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 8  | no  | 27/2  | 32       | 0.5260  | 0.3656 | 2.7      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 4  | no  | 1/7   | 1        | 0.5226  | 0.5403 | 3.3      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.5403 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 12 | yes | 37/2  | 42       | 0.5230  | 0.3476 | 1.7      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 8  | yes | 1/7   | 1        | 0.5226  | 0.4735 | 2.3      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5226  | 0.4735 | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | no  | 70/0  | 36       | 0.6364  | 0.2284 | 1.2      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7444  | 0.6516 | 0.9      | 🔴      |
| WSINDy      | Clean | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.3865 | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | yes | 84/0  | 41       | 0.5188  | 0.2838 | 1.3      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7479  | 0.6696 | 1.0      | 🔴      |
| WSINDy      | Clean | RFE       |         | 8  | yes | 1/7   | 1        | 0.5226  | 0.3592 | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | no  | 59/0  | 34       | 0.6002  | 0.2366 | 1.3      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7449  | 0.6540 | 1.0      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.4199 | 0.1      | 🔴      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | yes | 89/0  | 42       | 0.5189  | 0.2956 | 1.3      | 🔴      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7487  | 0.6738 | 0.9      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5226  | 0.3821 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 12 | no  | 59/0  | 34       | 0.5998  | 0.2733 | 1.1      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 0.7472  | 0.6660 | 0.9      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.5403 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 8  | yes | 92/0  | 44       | 0.5212  | 0.4714 | 1.9      | 🔴      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 0.7529  | 0.6942 | 0.9      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5226  | 0.4735 | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 10.0    | 4  | no  | 0/0   | 7        | 0.0884  | 0.0719 | 2.6      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.5     | 8  | no  | 0/0   | 7        | 0.0392  | 0.0110 | 6.1      | 🟢      |
| PolyJet     | Clean | RFE       |         | 8  | no  | 0/0   | 7        | 0.0392  | 0.0110 | 0.2      | 🟢      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 8        | 0.0304  | 0.0088 | 1.3      | 🔴      |
| PolyJet     | Clean | MSTLS     | 0.25    | 8  | yes | 0/0   | 7        | 0.0288  | 0.0146 | 6.1      | 🟢      |
| PolyJet     | Clean | RFE       |         | 8  | yes | 0/0   | 7        | 0.0288  | 0.0146 | 0.2      | 🟢      |
| PolyJet     | 2%    | STRidge   | 10.0    | 8  | no  | 14/0  | 21       | 0.2636  | 0.1447 | 2.3      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.5     | 4  | no  | 8/0   | 15       | 0.3319  | 0.2721 | 9.9      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 8  | no  | 0/6   | 1        | 0.6793  | 0.3243 | 0.2      | 🔴      |
| PolyJet     | 2%    | STRidge   | 10.0    | 12 | yes | 19/0  | 26       | 0.2630  | 0.1351 | 1.7      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.75    | 12 | yes | 9/0   | 16       | 0.2682  | 0.1359 | 4.2      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 8  | yes | 0/6   | 1        | 0.6527  | 0.4623 | 0.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 4  | no  | 65/0  | 72       | 0.3700  | 0.4983 | 10.3     | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 4  | no  | 8/2   | 13       | 0.3742  | 0.4999 | 9.8      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 0.5226  | 0.5403 | 0.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | yes | 42/0  | 49       | 0.2997  | 0.2840 | 3.0      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.75    | 12 | yes | 9/2   | 14       | 0.3160  | 0.2869 | 4.4      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 1/7   | 1        | 0.5226  | 0.4735 | 0.2      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0351  | 0.0071 | 0.3      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.75    | 12 | no  | 0/0   | 3        | 0.0350  | 0.0071 | 0.3      | 🟢      |
| Multigraphs | Clean | RFE       |         | 4  | no  | 0/5   | 1        | 0.3492  | 0.2638 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge   | 10.0    | 8  | yes | 0/0   | 3        | 0.0277  | 0.0156 | 0.3      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.5     | 8  | yes | 0/0   | 3        | 0.0277  | 0.0156 | 0.5      | 🟢      |
| Multigraphs | Clean | RFE       |         | 8  | yes | 0/0   | 3        | 0.0277  | 0.0156 | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 4  | no  | 11/0  | 9        | 0.2738  | 0.2789 | 0.7      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 5/0   | 5        | 0.1162  | 0.1483 | 0.5      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.3436  | 0.1912 | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 1.0     | 8  | yes | 27/0  | 13       | 0.1971  | 0.2226 | 0.6      | 🔴      |
| Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 5/0   | 5        | 0.1101  | 0.1395 | 0.4      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.3425  | 0.1803 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 1.0     | 4  | no  | 19/0  | 12       | 0.3569  | 0.5045 | 0.7      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 8  | no  | 9/3   | 6        | 0.3354  | 0.3114 | 0.5      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.6270  | 0.2652 | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.1     | 8  | yes | 32/0  | 18       | 0.3162  | 0.4245 | 0.6      | 🔴      |
| Multigraphs | 5%    | MSTLS     | 0.75    | 12 | yes | 9/3   | 6        | 0.3358  | 0.2945 | 0.3      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.6313  | 0.3432 | 0.1      | 🔴      |

---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ks_chaos`.

| Source                                    | Wall time | CPU                                      | Threads |
|:------------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260802_110219_ks_chaos_full.csv` | 162.2 min | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | 1       |
