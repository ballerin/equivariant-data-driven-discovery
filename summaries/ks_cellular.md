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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best         | Mg E_data | Std best      | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-------------------|:----------|:--------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | 🟩        | 0/0 ✓             | 0.0040     | multigraph/RFE     | 0.0123    | PolyJet/RFE   | 0/0 ✓     | 0.0041      | 0.0122     | 0.0129       |
| Clean | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.0040     | multigraph/RFE     | 0.0123    | PolyJet/RFE   | 0/0 ✓     | 0.0041      | 0.0122     | 0.0129       |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0027     | multigraph/STRidge | 0.0088    | PolyJet/MSTLS | 0/0 ✓     | 0.0030      | 0.0088     | 0.0092       |
| Clean | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0027     | multigraph/STRidge | 0.0088    | PolyJet/RFE   | 0/0 ✓     | 0.0030      | 0.0088     | 0.0092       |
| Clean | weak (IBP)    | spectral   | 🟢        | 0/0 ✓             | 0.1234     | multigraph/MSTLS   | 0.0760    | —             | 1/0 ✗     | 0.1388      | 0.0081     | 0.1572       |
| Clean | weak (IBP)    | polynomial | 🟢        | 0/0 ✓             | 0.1234     | multigraph/MSTLS   | 0.0760    | —             | 1/0 ✗     | 0.1388      | 0.0081     | 0.1572       |
| 2%    | pointwise     | spectral   | —        | 19/0 ✗            | 1.0000     | —                  | 0.5297    | —             | 88/0 ✗    | 1.0000      | 0.4713     | 9249.6733    |
| 2%    | pointwise     | polynomial | —        | 37/0 ✗            | 0.9999     | —                  | 0.4552    | —             | 93/0 ✗    | 1.0000      | 0.4002     | 204.6417     |
| 2%    | weak (no IBP) | spectral   | —        | 9/0 ✗             | 0.9806     | —                  | 0.3719    | —             | 35/0 ✗    | 0.8094      | 0.2031     | 2.9497       |
| 2%    | weak (no IBP) | polynomial | 🟢        | 0/0 ✓             | 0.5015     | multigraph/MSTLS   | 0.2656    | —             | 9/0 ✗     | 0.4473      | 0.2086     | 0.4322       |
| 2%    | weak (IBP)    | spectral   | —        | 8/0 ✗             | 0.6349     | —                  | 0.3222    | —             | 35/0 ✗    | 0.8090      | 0.2033     | 3.0651       |
| 2%    | weak (IBP)    | polynomial | 🟢        | 0/0 ✓             | 0.5254     | multigraph/MSTLS   | 0.2617    | —             | 4/0 ✗     | 0.4364      | 0.2149     | 0.5150       |
| 5%    | pointwise     | spectral   | —        | 24/0 ✗            | 1.0000     | —                  | 0.6742    | —             | 49/0 ✗    | 1.0000      | 0.6532     | 20670.3966   |
| 5%    | pointwise     | polynomial | —        | 37/0 ✗            | 1.0000     | —                  | 0.6474    | —             | 89/0 ✗    | 1.0000      | 0.6009     | 457.3188     |
| 5%    | weak (no IBP) | spectral   | —        | 15/0 ✗            | 0.9980     | —                  | 0.4714    | —             | 38/0 ✗    | 0.9975      | 0.3841     | 18.2763      |
| 5%    | weak (no IBP) | polynomial | —        | 7/0 ✗             | 0.8607     | —                  | 0.3259    | —             | 6/0 ✗     | 0.7824      | 0.3192     | 1.6768       |
| 5%    | weak (IBP)    | spectral   | —        | 4/0 ✗             | 0.9427     | —                  | 0.4545    | —             | 38/0 ✗    | 0.9975      | 0.3841     | 18.3906      |
| 5%    | weak (IBP)    | polynomial | —        | 7/0 ✗             | 0.8630     | —                  | 0.3265    | —             | 2/0 ✗     | 0.7414      | 0.3447     | 1.7805       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:------|--------:|---------:|---------:|
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0027  | 0.1      | —        |
| Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0027  | 0.1      | +0.0%    |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0027  | 0.1      | +0.2%    |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0027  | 0.1      | +0.2%    |
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0027  | 0.2      | +0.2%    |
| Modified Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0027  | 0.2      | +0.2%    |
| PolyJet + RFE                  | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0030  | 0.1      | +13.9%   |
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0030  | 1.0      | +13.9%   |
| PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0030  | 0.3      | +14.2%   |

### 2% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.5015  | 0.2      | —        |
| Modified Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.5015  | 0.2      | +0.0%    |
| PolyJet + MSTLS                | polynomial, weak (IBP), m=12    | 4/0 ✗  | 0.4364  | 0.8      | ✗        |
| Modified Multigraphs + STRidge | polynomial, weak (IBP), m=8     | 16/0 ✗ | 0.9842  | 0.2      | ✗        |
| Multigraphs + STRidge          | spectral, weak (IBP), m=8       | 17/0 ✗ | 0.9817  | 0.1      | ✗        |
| PolyJet + STRidge              | polynomial, weak (IBP), m=12    | 25/0 ✗ | 0.5573  | 0.3      | ✗        |
| Multigraphs + RFE              | spectral, weak (no IBP), m=8    | 0/5 ✗  | 1.0000  | 0.1      | ✗        |
| Modified Multigraphs + RFE     | spectral, weak (no IBP), m=8    | 0/5 ✗  | 1.0000  | 0.1      | ✗        |
| PolyJet + RFE                  | spectral, weak (no IBP), m=12   | 1/7 ✗  | 1.0000  | 0.1      | ✗        |

### 5% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| PolyJet + MSTLS                | polynomial, weak (IBP), m=12    | 2/0 ✗  | 0.7414  | 0.7      | ✗        |
| Modified Multigraphs + MSTLS   | spectral, weak (IBP), m=12      | 4/0 ✗  | 0.9427  | 0.2      | ✗        |
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=12 | 7/0 ✗  | 0.8607  | 0.2      | ✗        |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 18/0 ✗ | 0.8812  | 0.1      | ✗        |
| Modified Multigraphs + STRidge | polynomial, weak (IBP), m=4     | 29/0 ✗ | 1.0000  | 0.1      | ✗        |
| PolyJet + STRidge              | polynomial, weak (IBP), m=12    | 41/0 ✗ | 0.8246  | 0.4      | ✗        |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=8  | 0/5 ✗  | 1.0000  | 0.1      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=8  | 0/5 ✗  | 1.0000  | 0.1      | ✗        |
| PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/6 ✗  | 1.0000  | 0.1      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PolyJet          | Clean | STRidge   | 1.0     | 0/0   | 7        | 0.0041  | 0.0122 | 0.4      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0041  | 0.0122 | 0.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.01    | 88/0  | 95       | 1.0000  | 0.4713 | 0.7      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 0.3      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 49/0  | 56       | 1.0000  | 0.6532 | 0.6      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.7219 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 19/0  | 12       | 1.0000  | 0.5297 | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.1     | 24/0  | 14       | 1.0000  | 0.6742 | 0.1      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.7219 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.001   | 37/0  | 20       | 1.0000  | 0.5297 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 32/0  | 16       | 1.0000  | 0.6743 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.7219 | 0.1      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PolyJet          | Clean | STRidge   | 1.0     | 0/0   | 7        | 0.0041  | 0.0122 | 0.4      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0041  | 0.0122 | 0.2      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.01    | 93/0  | 100      | 1.0000  | 0.4002 | 0.7      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.01    | 89/0  | 96       | 1.0000  | 0.6009 | 0.6      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.7219 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0040  | 0.0123 | 0.2      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.001   | 37/0  | 20       | 0.9999  | 0.4551 | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 37/0  | 20       | 1.0000  | 0.6474 | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.7219 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0040  | 0.0123 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 37/0  | 15       | 0.9999  | 0.4552 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6309 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 37/0  | 20       | 1.0000  | 0.6474 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.7219 | 0.1      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PolyJet          | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0030  | 0.0088 | 0.3      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0030  | 0.0088 | 0.7      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0030  | 0.0088 | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 8        | 0.1388  | 0.0081 | 0.3      | 🔴      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 1/0   | 8        | 0.1389  | 0.0081 | 1.0      | 🔴      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 1/0   | 8        | 0.1389  | 0.0081 | 0.2      | 🔴      |
| PolyJet          | 2%    | STRidge   | 1.0     | 8  | no  | 58/0  | 65       | 0.9938  | 0.3054 | 0.5      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 12 | no  | 35/0  | 42       | 0.8094  | 0.2031 | 0.8      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 1.0000  | 0.6415 | 0.1      | 🔴      |
| PolyJet          | 2%    | STRidge   | 1.0     | 8  | yes | 45/0  | 52       | 0.9943  | 0.3075 | 0.5      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 12 | yes | 35/0  | 42       | 0.8090  | 0.2033 | 0.8      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000  | 0.6415 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 4  | no  | 81/0  | 88       | 1.0000  | 0.4365 | 0.7      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 8  | no  | 38/0  | 45       | 0.9975  | 0.3841 | 1.0      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 8  | no  | 0/6   | 1        | 1.0000  | 0.6513 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | yes | 79/0  | 86       | 0.9974  | 0.3784 | 0.6      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 8  | yes | 38/0  | 45       | 0.9975  | 0.3841 | 1.2      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 0/6   | 1        | 1.0000  | 0.6513 | 0.2      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0027  | 0.0088 | 0.1      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0027  | 0.0088 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0027  | 0.0088 | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 4        | 0.1387  | 0.0082 | 0.1      | 🔴      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.1234  | 0.0760 | 0.2      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 1/0   | 4        | 0.1388  | 0.0082 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | no  | 25/0  | 14       | 0.5878  | 0.2571 | 0.1      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 8  | no  | 9/0   | 8        | 0.9806  | 0.3719 | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 8  | no  | 0/5   | 1        | 1.0000  | 0.5441 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 8  | yes | 17/0  | 10       | 0.9817  | 0.3723 | 0.1      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | yes | 13/0  | 8        | 0.6111  | 0.2535 | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 8  | yes | 0/5   | 1        | 1.0000  | 0.5441 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | no  | 31/0  | 17       | 0.9245  | 0.3920 | 0.1      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 12 | no  | 18/0  | 12       | 0.9281  | 0.4005 | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 1.0000  | 0.6154 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 4  | yes | 23/0  | 14       | 1.0001  | 0.5219 | 0.1      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 4  | yes | 12/0  | 10       | 1.0001  | 0.5273 | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 4  | yes | 1/7   | 1        | 1.0000  | 0.6154 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0027  | 0.0088 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0027  | 0.0088 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0027  | 0.0088 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 4        | 0.1387  | 0.0082 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.1234  | 0.0760 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 1/0   | 4        | 0.1388  | 0.0082 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | no  | 25/0  | 15       | 0.5827  | 0.2547 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 8  | no  | 12/0  | 9        | 0.9806  | 0.3719 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 8  | no  | 0/5   | 1        | 1.0000  | 0.5441 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | yes | 33/0  | 15       | 0.6876  | 0.2405 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.5     | 12 | yes | 8/0   | 5        | 0.6349  | 0.3222 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 8  | yes | 0/5   | 1        | 1.0000  | 0.5441 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 8  | no  | 35/0  | 16       | 0.9979  | 0.4633 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 8  | no  | 15/0  | 11       | 0.9980  | 0.4714 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | no  | 1/7   | 1        | 1.0000  | 0.5955 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 12 | yes | 37/0  | 20       | 0.9398  | 0.4049 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 12 | yes | 4/0   | 6        | 0.9427  | 0.4545 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000  | 0.5687 | 0.1      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PolyJet          | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 7        | 0.0030  | 0.0088 | 0.3      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0030  | 0.0088 | 1.0      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0030  | 0.0088 | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 8        | 0.1388  | 0.0081 | 0.4      | 🔴      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 1/0   | 8        | 0.1389  | 0.0081 | 1.0      | 🔴      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 1/0   | 8        | 0.1389  | 0.0081 | 0.2      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 12 | no  | 63/0  | 70       | 0.5965  | 0.1642 | 0.6      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 12 | no  | 9/0   | 16       | 0.4473  | 0.2086 | 0.8      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 1.0000  | 0.5668 | 0.1      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 12 | yes | 25/0  | 32       | 0.5573  | 0.1792 | 0.3      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 12 | yes | 4/0   | 11       | 0.4364  | 0.2149 | 0.8      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000  | 0.5668 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 10.0    | 12 | no  | 41/0  | 48       | 0.8247  | 0.2458 | 0.3      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | no  | 6/0   | 13       | 0.7824  | 0.3192 | 0.7      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 1.0000  | 0.5687 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 10.0    | 12 | yes | 41/0  | 48       | 0.8246  | 0.2456 | 0.4      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | yes | 2/0   | 9        | 0.7414  | 0.3447 | 0.7      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000  | 0.5687 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0027  | 0.0088 | 0.1      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0027  | 0.0088 | 0.2      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0027  | 0.0088 | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 4        | 0.1387  | 0.0082 | 0.1      | 🔴      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.1234  | 0.0760 | 0.2      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 1/0   | 4        | 0.1388  | 0.0082 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | no  | 21/0  | 13       | 0.5594  | 0.2409 | 0.1      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.5015  | 0.2656 | 0.2      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 1.0000  | 0.5668 | 0.1      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | yes | 21/0  | 13       | 0.5643  | 0.2425 | 0.1      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.5254  | 0.2617 | 0.2      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000  | 0.5668 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | no  | 18/0  | 12       | 0.8812  | 0.3113 | 0.1      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 12 | no  | 7/0   | 7        | 0.8607  | 0.3259 | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | no  | 0/5   | 1        | 1.0000  | 0.4886 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | yes | 18/0  | 12       | 0.8833  | 0.3117 | 0.1      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 12 | yes | 7/0   | 7        | 0.8630  | 0.3265 | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 0/5   | 1        | 1.0000  | 0.4886 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0027  | 0.0088 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0027  | 0.0088 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0027  | 0.0088 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 4        | 0.1387  | 0.0082 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.1234  | 0.0760 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 1/0   | 4        | 0.1388  | 0.0082 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 0.001   | 12 | no  | 37/0  | 20       | 0.6441  | 0.2182 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.5015  | 0.2656 | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 1/7   | 1        | 1.0000  | 0.5668 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 8  | yes | 16/0  | 9        | 0.9842  | 0.3109 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.5     | 12 | yes | 8/0   | 5        | 0.6079  | 0.2338 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 1/7   | 1        | 1.0000  | 0.5668 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 8  | no  | 33/0  | 16       | 0.9951  | 0.3718 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 12 | no  | 7/0   | 7        | 0.8607  | 0.3259 | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | no  | 0/5   | 1        | 1.0000  | 0.4886 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 4  | yes | 29/0  | 12       | 1.0000  | 0.4343 | 0.1      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 12 | yes | 9/0   | 7        | 0.8732  | 0.3147 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | yes | 0/5   | 1        | 1.0000  | 0.4886 | 0.1      | 🔴      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ)   |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:----------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0126       | 105502.8  |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0114       | 114784.0  |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0092       | 125064.5  |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.0126       | 105504.4  |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0114       | 114785.6  |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0092       | 125064.7  |
| Clean | weak (IBP)    | spectral   | 4  | yes | 864.6296     | 6541786.5 |
| Clean | weak (IBP)    | spectral   | 8  | yes | 3.7624       | 133848.0  |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.1572       | 126840.5  |
| Clean | weak (IBP)    | polynomial | 4  | yes | 864.6296     | 6541669.1 |
| Clean | weak (IBP)    | polynomial | 8  | yes | 3.7624       | 133850.3  |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.1572       | 126842.9  |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 88.9904      | 336798.9  |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 3.9765       | 471736.7  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 2.9497       | 693750.9  |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 33.2300      | 55740.0   |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 2.5916       | 54983.2   |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.4322       | 75641.8   |
| 2%    | weak (IBP)    | spectral   | 4  | yes | 868.8468     | 338696.5  |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 6.7089       | 472036.4  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 3.0651       | 693845.0  |
| 2%    | weak (IBP)    | polynomial | 4  | yes | 867.1098     | 55767.2   |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 4.8434       | 54914.0   |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.5150       | 75666.7   |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 220.7844     | 2029412.3 |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 18.7254      | 3418856.0 |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 18.2763      | 4850904.7 |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 82.2495      | 149255.4  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 6.5727       | 63243.2   |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 1.6768       | 85906.5   |
| 5%    | weak (IBP)    | spectral   | 4  | yes | 892.4608     | 2031277.0 |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 21.4031      | 3419390.6 |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 18.3906      | 4851117.8 |
| 5%    | weak (IBP)    | polynomial | 4  | yes | 881.9141     | 149255.1  |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 8.5625       | 63144.3   |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 1.7805       | 85896.2   |
---
| Source                                      | Wall time | CPU | Threads |
|:--------------------------------------------|----------:|:----|:-------:|
| `sweep_20260715_180039_ks_cell_trimmed.csv` | ? min     | ?   | ?       |
