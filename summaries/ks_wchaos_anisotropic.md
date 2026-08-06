# Kuramoto-Sivashinsky equation (weakly chaotic, $L=8\pi$) — Anisotropic noise

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

All libraries are evaluated in the common polynomial-jet basis (max order 4,
max degree 3).  fp/fn/E are computed against the same ground truth for every
library.  KS contains the gradient norm $\|\nabla f\|^2 = f_x^2 + f_y^2$ which neither PDE-FIND nor WSINDy can represent as a sparse column (PDE-FIND lacks products of derivatives; in WSINDy, $f_x^2$ is locked with $f\cdot f_{xx}$ inside $\partial_{xx}(f^2)$).  Weakly chaotic regime ($L=8\pi$).

**Noise model**: Space-dependent observation noise
$\sigma \cdot P(x,y) \cdot \mathcal{N}(0,1)$ where
$\sigma = \mathrm{obs\_noise} \cdot \mathrm{std}(f)$ and $P(x,y)$ is a
normalised spatial profile ($\mathrm{mean}(P^2)=1$, preserving total noise
power relative to isotropic).

| Profile | $P(x,y)$ | Symmetry |
|---|---|---|
| Radial decay | $\propto \exp(-r^2/2s^2)$, $r=\sqrt{x^2+y^2}$, $s=L/4$ | O(2) |
| Radial growth | $\propto 1 + (1 - \exp(-r^2/2s^2))$ | O(2) |
| Axial decay | $\propto \exp(-x^2/2s^2)$ | Reflection only |
| Axial growth | $\propto 1 + (1 - \exp(-x^2/2s^2))$ | Reflection only |

Sweep configuration: 2 libraries (Multigraphs, PolyJet),
2 algorithm(s) (RFE, STRidge),
2 method(s) (weak (no IBP), weak (IBP)),
2 backend(s) (Polynomial stencils, Spectral derivatives),
2 noise level(s),
4 noise profile(s) (Axial decay, Radial decay, Axial growth, Radial growth).

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


| Noise | Method        | Backend    | Profile       | Eq best | Multigraphs fp/fn | Multigraphs E | Mult. best | Non-eq best | Non-eq fp/fn | Non-eq E |
|:------|:-------------:|:----------:|:-------------:|:-------:|:------------------|:--------------|:-----------|:------------|:-------------|:---------|
| 2%    | weak (no IBP) | spectral   | Radial decay  | —       | 30/0 ✗            | 0.4151        | —          | —           | 81/0 ✗       | 0.4035   |
| 2%    | weak (no IBP) | spectral   | Radial growth | —       | 20/0 ✗            | 0.3217        | —          | —           | 61/0 ✗       | 0.2750   |
| 2%    | weak (no IBP) | spectral   | Axial decay   | —       | 34/0 ✗            | 0.4014        | —          | —           | 80/0 ✗       | 0.3738   |
| 2%    | weak (no IBP) | spectral   | Axial growth  | —       | 20/0 ✗            | 0.3152        | —          | —           | 59/0 ✗       | 0.2722   |
| 2%    | weak (no IBP) | polynomial | Radial decay  | —       | 16/0 ✗            | 0.1586        | —          | —           | 23/0 ✗       | 0.3138   |
| 2%    | weak (no IBP) | polynomial | Radial growth | —       | 9/0 ✗             | 0.0852        | —          | —           | 7/0 ✗        | 0.0704   |
| 2%    | weak (no IBP) | polynomial | Axial decay   | —       | 11/0 ✗            | 0.1303        | —          | —           | 19/0 ✗       | 0.1321   |
| 2%    | weak (no IBP) | polynomial | Axial growth  | —       | 8/0 ✗             | 0.0962        | —          | —           | 9/0 ✗        | 0.2549   |
| 2%    | weak (IBP)    | spectral   | Radial decay  | —       | 27/0 ✗            | 0.4237        | —          | —           | 46/0 ✗       | 0.4047   |
| 2%    | weak (IBP)    | spectral   | Radial growth | —       | 21/0 ✗            | 0.4192        | —          | —           | 61/0 ✗       | 0.2786   |
| 2%    | weak (IBP)    | spectral   | Axial decay   | —       | 27/0 ✗            | 0.4202        | —          | —           | 76/0 ✗       | 0.3839   |
| 2%    | weak (IBP)    | spectral   | Axial growth  | —       | 24/0 ✗            | 0.3165        | —          | —           | 35/0 ✗       | 0.2688   |
| 2%    | weak (IBP)    | polynomial | Radial decay  | —       | 10/0 ✗            | 0.3540        | —          | —           | 42/0 ✗       | 0.3217   |
| 2%    | weak (IBP)    | polynomial | Radial growth | —       | 9/0 ✗             | 0.0965        | —          | —           | 14/0 ✗       | 0.0820   |
| 2%    | weak (IBP)    | polynomial | Axial decay   | —       | 7/0 ✗             | 0.1778        | —          | —           | 14/0 ✗       | 0.1322   |
| 2%    | weak (IBP)    | polynomial | Axial growth  | —       | 7/0 ✗             | 0.1677        | —          | —           | 14/0 ✗       | 0.0818   |
| 5%    | weak (no IBP) | spectral   | Radial decay  | —       | 33/0 ✗            | 0.4325        | —          | —           | 93/0 ✗       | 0.4294   |
| 5%    | weak (no IBP) | spectral   | Radial growth | —       | 35/0 ✗            | 0.4306        | —          | —           | 65/0 ✗       | 0.4284   |
| 5%    | weak (no IBP) | spectral   | Axial decay   | —       | 35/0 ✗            | 0.4328        | —          | —           | 79/0 ✗       | 0.4291   |
| 5%    | weak (no IBP) | spectral   | Axial growth  | —       | 28/0 ✗            | 0.4290        | —          | —           | 58/0 ✗       | 0.4336   |
| 5%    | weak (no IBP) | polynomial | Radial decay  | —       | 27/0 ✗            | 0.3923        | —          | —           | 90/0 ✗       | 0.3649   |
| 5%    | weak (no IBP) | polynomial | Radial growth | —       | 18/0 ✗            | 0.2552        | —          | —           | 29/0 ✗       | 0.2279   |
| 5%    | weak (no IBP) | polynomial | Axial decay   | —       | 30/0 ✗            | 0.3730        | —          | —           | 80/0 ✗       | 0.3055   |
| 5%    | weak (no IBP) | polynomial | Axial growth  | —       | 18/0 ✗            | 0.2593        | —          | —           | 17/0 ✗       | 0.2298   |
| 5%    | weak (IBP)    | spectral   | Radial decay  | —       | 28/0 ✗            | 0.4328        | —          | —           | 75/0 ✗       | 0.4296   |
| 5%    | weak (IBP)    | spectral   | Radial growth | —       | 23/0 ✗            | 0.4327        | —          | —           | 65/0 ✗       | 0.4288   |
| 5%    | weak (IBP)    | spectral   | Axial decay   | —       | 34/0 ✗            | 0.4319        | —          | —           | 44/0 ✗       | 0.4292   |
| 5%    | weak (IBP)    | spectral   | Axial growth  | —       | 28/0 ✗            | 0.4289        | —          | —           | 95/0 ✗       | 0.4359   |
| 5%    | weak (IBP)    | polynomial | Radial decay  | —       | 27/0 ✗            | 0.3923        | —          | —           | 71/0 ✗       | 0.3648   |
| 5%    | weak (IBP)    | polynomial | Radial growth | —       | 18/0 ✗            | 0.2586        | —          | —           | 29/0 ✗       | 0.2296   |
| 5%    | weak (IBP)    | polynomial | Axial decay   | —       | 27/0 ✗            | 0.4005        | —          | —           | 69/0 ✗       | 0.3463   |
| 5%    | weak (IBP)    | polynomial | Axial growth  | —       | 17/0 ✗            | 0.2621        | —          | —           | 37/0 ✗       | 0.2344   |

## Ablation

### 2% — Axial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 7/0   | 0.1778  | 0.5      | ✗        |
| 2    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 11/0  | 0.1303  | 0.4      | ✗        |
| 3    | PolyJet + STRidge              | polynomial, weak (IBP), m=12    | 14/0  | 0.1322  | 1.8      | ✗        |
| 4    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/5   | 0.6336  | 0.1      | ✗        |
| 5    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5   | 0.6336  | 0.2      | ✗        |
| 6    | PolyJet + RFE                  | polynomial, weak (no IBP), m=12 | 0/6   | 0.7243  | 0.2      | ✗        |

### 2% — Radial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=8     | 10/0  | 0.3540  | 0.4      | ✗        |
| 2    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 16/0  | 0.1586  | 0.4      | ✗        |
| 3    | PolyJet + STRidge              | polynomial, weak (no IBP), m=8  | 23/0  | 0.3138  | 2.6      | ✗        |
| 4    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/5   | 0.7211  | 0.1      | ✗        |
| 5    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5   | 0.7211  | 0.1      | ✗        |
| 6    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 1/7   | 0.5796  | 0.2      | ✗        |

### 2% — Axial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 7/0   | 0.1677  | 0.3      | ✗        |
| 2    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 9/0   | 0.0871  | 0.4      | ✗        |
| 3    | PolyJet + STRidge              | polynomial, weak (no IBP), m=8  | 9/0   | 0.2549  | 2.2      | ✗        |
| 4    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/5   | 0.5343  | 0.1      | ✗        |
| 5    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5   | 0.5343  | 0.1      | ✗        |
| 6    | PolyJet + RFE                  | polynomial, weak (no IBP), m=12 | 0/6   | 0.7120  | 0.2      | ✗        |

### 2% — Radial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 7/0   | 0.0704  | 1.5      | ✗        |
| 2    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 9/0   | 0.0852  | 0.4      | ✗        |
| 3    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 10/0  | 0.0847  | 0.3      | ✗        |
| 4    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/5   | 0.5303  | 0.1      | ✗        |
| 5    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5   | 0.5303  | 0.1      | ✗        |
| 6    | PolyJet + RFE                  | polynomial, weak (no IBP), m=8  | 0/6   | 0.7086  | 0.2      | ✗        |

### 5% — Axial decay

| Rank | Library + Algorithm            | Configuration                  | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=8    | 27/0  | 0.4005  | 1.0      | ✗        |
| 2    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=8 | 32/0  | 0.3731  | 0.6      | ✗        |
| 3    | PolyJet + STRidge              | spectral, weak (IBP), m=8      | 44/0  | 0.4292  | 8.6      | ✗        |
| 4    | Modified Multigraphs + RFE     | spectral, weak (no IBP), m=8   | 1/7   | 0.5796  | 0.1      | ✗        |
| 5    | PolyJet + RFE                  | polynomial, weak (no IBP), m=8 | 1/7   | 0.5796  | 0.2      | ✗        |
| 6    | Multigraphs + RFE              | spectral, weak (no IBP), m=8   | 1/7   | 0.5796  | 0.2      | ✗        |

### 5% — Radial decay

| Rank | Library + Algorithm            | Configuration                  | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + STRidge          | polynomial, weak (IBP), m=12   | 27/0  | 0.3923  | 0.4      | ✗        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12     | 28/0  | 0.4328  | 0.5      | ✗        |
| 3    | PolyJet + STRidge              | polynomial, weak (IBP), m=8    | 71/0  | 0.3648  | 6.7      | ✗        |
| 4    | Modified Multigraphs + RFE     | spectral, weak (no IBP), m=8   | 1/7   | 0.5796  | 0.1      | ✗        |
| 5    | PolyJet + RFE                  | polynomial, weak (no IBP), m=8 | 1/7   | 0.5796  | 0.2      | ✗        |
| 6    | Multigraphs + RFE              | spectral, weak (no IBP), m=8   | 1/7   | 0.5796  | 0.2      | ✗        |

### 5% — Axial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 17/0  | 0.2298  | 2.2      | ✗        |
| 2    | Multigraphs + STRidge          | polynomial, weak (IBP), m=12    | 17/0  | 0.2621  | 0.4      | ✗        |
| 3    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 18/0  | 0.2593  | 0.7      | ✗        |
| 4    | Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/5   | 0.8295  | 0.2      | ✗        |
| 5    | PolyJet + RFE                  | polynomial, weak (no IBP), m=8  | 1/7   | 0.5796  | 0.2      | ✗        |
| 6    | Multigraphs + RFE              | spectral, weak (no IBP), m=8    | 1/7   | 0.5796  | 0.2      | ✗        |

### 5% — Radial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 18/0  | 0.2552  | 0.4      | ✗        |
| 2    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 18/0  | 0.2552  | 0.8      | ✗        |
| 3    | PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 29/0  | 0.2279  | 2.5      | ✗        |
| 4    | Multigraphs + RFE              | spectral, weak (no IBP), m=8    | 1/7   | 0.5796  | 0.2      | ✗        |
| 5    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=8  | 1/7   | 0.5796  | 0.2      | ✗        |
| 6    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 1/7   | 0.5796  | 0.4      | ✗        |


## Detailed results

### weak (no IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/6   | 1        | 0.7243  | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/6   | 1        | 0.7120  | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial growth | 0/6   | 1        | 0.7086  | 0.2      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 19/0  | 26       | 0.1321  | 1.9      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Radial decay  | 23/0  | 30       | 0.3138  | 2.6      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Axial growth  | 9/0   | 16       | 0.2549  | 2.2      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial growth | 7/0   | 14       | 0.0704  | 1.5      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 80/0  | 87       | 0.3055  | 5.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Radial decay  | 90/0  | 97       | 0.3649  | 6.9      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 17/0  | 24       | 0.2298  | 2.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Radial growth | 29/0  | 36       | 0.2279  | 2.5      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/5   | 1        | 0.6336  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/5   | 1        | 0.7211  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/5   | 1        | 0.5343  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/5   | 1        | 0.5303  | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 11/0  | 9        | 0.1303  | 0.4      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 16/0  | 12       | 0.1586  | 0.4      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 8/0   | 6        | 0.0962  | 0.4      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial growth | 9/0   | 7        | 0.0852  | 0.4      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.5796  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/7   | 1        | 0.5796  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.5796  | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 1.0       | Axial decay   | 30/0  | 15       | 0.3730  | 1.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 27/0  | 16       | 0.3923  | 0.4      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 18/0  | 13       | 0.2593  | 0.7      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial growth | 18/0  | 13       | 0.2552  | 0.4      | 🔴      |

### weak (no IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/6   | 1        | 0.9118  | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/6   | 1        | 0.9081  | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 80/0  | 87       | 0.3738  | 4.6      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 1.0       | Radial decay  | 81/0  | 88       | 0.4035  | 6.2      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 59/0  | 66       | 0.2722  | 3.3      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial growth | 61/0  | 68       | 0.2750  | 3.5      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.5796  | 0.4      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.3      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/7   | 1        | 0.5796  | 0.4      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.5796  | 0.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 1.0       | Axial decay   | 79/0  | 86       | 0.4291  | 15.4     | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.01      | Radial decay  | 93/0  | 100      | 0.4294  | 15.7     | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 58/0  | 65       | 0.4336  | 6.3      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Radial growth | 65/0  | 72       | 0.4284  | 11.2     | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/5   | 1        | 0.9080  | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/5   | 1        | 0.9042  | 0.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 0.1       | Axial decay   | 34/0  | 18       | 0.4014  | 1.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 1.0       | Radial decay  | 30/0  | 16       | 0.4151  | 1.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 20/0  | 12       | 0.3152  | 0.8      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial growth | 20/0  | 12       | 0.3217  | 0.8      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.5796  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/7   | 1        | 0.5796  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 0.1       | Axial decay   | 35/0  | 19       | 0.4328  | 0.6      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 33/0  | 18       | 0.4325  | 0.8      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Axial growth  | 28/0  | 16       | 0.4290  | 1.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 0.1       | Radial growth | 35/0  | 18       | 0.4306  | 0.5      | 🔴      |

### weak (IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/6   | 1        | 0.7243  | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/6   | 1        | 0.7120  | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial growth | 0/6   | 1        | 0.7086  | 0.2      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 14/0  | 21       | 0.1322  | 1.8      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 1.0       | Radial decay  | 42/0  | 49       | 0.3217  | 3.6      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 14/0  | 21       | 0.0818  | 1.7      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial growth | 14/0  | 21       | 0.0820  | 1.6      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 69/0  | 76       | 0.3463  | 7.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Radial decay  | 71/0  | 78       | 0.3648  | 6.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 37/0  | 44       | 0.2344  | 2.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Radial growth | 29/0  | 36       | 0.2296  | 2.6      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/5   | 1        | 0.6336  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/5   | 1        | 0.7211  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/5   | 1        | 0.5343  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/5   | 1        | 0.5303  | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 7/0   | 8        | 0.1778  | 0.5      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 10.0      | Radial decay  | 10/0  | 10       | 0.3540  | 0.4      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 7/0   | 5        | 0.1677  | 0.3      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial growth | 9/0   | 7        | 0.0965  | 0.4      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.5796  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/5   | 1        | 0.8295  | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.5796  | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Axial decay   | 27/0  | 14       | 0.4005  | 1.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 27/0  | 16       | 0.3923  | 0.4      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 17/0  | 12       | 0.2621  | 0.4      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial growth | 18/0  | 13       | 0.2586  | 0.4      | 🔴      |

### weak (IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/6   | 1        | 0.9118  | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/6   | 1        | 0.9081  | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 0.1       | Axial decay   | 76/0  | 83       | 0.3839  | 7.5      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 1.0       | Radial decay  | 46/0  | 53       | 0.4047  | 4.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 35/0  | 42       | 0.2688  | 2.7      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial growth | 61/0  | 68       | 0.2786  | 3.4      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.5796  | 0.4      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.4      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/7   | 1        | 0.5796  | 0.4      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.5796  | 0.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 1.0       | Axial decay   | 44/0  | 51       | 0.4292  | 8.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Radial decay  | 75/0  | 82       | 0.4296  | 14.1     | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 0.1       | Axial growth  | 95/0  | 102      | 0.4359  | 9.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Radial growth | 65/0  | 72       | 0.4288  | 13.4     | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/5   | 1        | 0.9080  | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/5   | 1        | 0.9042  | 0.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 1.0       | Axial decay   | 27/0  | 14       | 0.4202  | 1.3      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 1.0       | Radial decay  | 27/0  | 14       | 0.4237  | 1.3      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 24/0  | 14       | 0.3165  | 0.7      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 1.0       | Radial growth | 21/0  | 13       | 0.4192  | 1.4      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.5796  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.5796  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/7   | 1        | 0.5796  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.5796  | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 34/0  | 19       | 0.4319  | 0.7      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial decay  | 28/0  | 14       | 0.4328  | 0.5      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Axial growth  | 28/0  | 16       | 0.4289  | 1.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial growth | 23/0  | 11       | 0.4327  | 0.4      | 🔴      |


---
*Source: `sweep_aniso_default_ks_wchaos_20260805_215410.csv`*

