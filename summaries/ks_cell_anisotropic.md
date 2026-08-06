# Kuramoto-Sivashinsky equation (cellular, $L=4\pi$) — Anisotropic noise

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
library.  KS contains the gradient norm $\|\nabla f\|^2 = f_x^2 + f_y^2$ which neither PDE-FIND nor WSINDy can represent as a sparse column (PDE-FIND lacks products of derivatives; in WSINDy, $f_x^2$ is locked with $f\cdot f_{xx}$ inside $\partial_{xx}(f^2)$).  Cellular regime ($L=4\pi$).

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


| Noise | Method        | Backend    | Profile       | Eq best | Multigraphs fp/fn | Multigraphs E | Mult. best             | Non-eq best | Non-eq fp/fn | Non-eq E |
|:------|:-------------:|:----------:|:-------------:|:-------:|:------------------|:--------------|:-----------------------|:------------|:-------------|:---------|
| 2%    | weak (no IBP) | spectral   | Radial decay  | —       | 23/0 ✗            | 0.4001        | —                      | —           | 56/0 ✗       | 0.4058   |
| 2%    | weak (no IBP) | spectral   | Radial growth | —       | 9/0 ✗             | 0.2554        | —                      | —           | 59/0 ✗       | 0.3877   |
| 2%    | weak (no IBP) | spectral   | Axial decay   | —       | 14/0 ✗            | 0.3518        | —                      | —           | 46/0 ✗       | 0.3665   |
| 2%    | weak (no IBP) | spectral   | Axial growth  | —       | 8/0 ✗             | 0.3579        | —                      | —           | 45/0 ✗       | 0.3680   |
| 2%    | weak (no IBP) | polynomial | Radial decay  | —       | 4/0 ✗             | 0.3258        | —                      | —           | 43/0 ✗       | 0.3620   |
| 2%    | weak (no IBP) | polynomial | Radial growth | —       | 3/0 ✗             | 0.2193        | —                      | —           | 30/0 ✗       | 0.2700   |
| 2%    | weak (no IBP) | polynomial | Axial decay   | —       | 3/0 ✗             | 0.3080        | —                      | —           | 38/0 ✗       | 0.3780   |
| 2%    | weak (no IBP) | polynomial | Axial growth  | 🟢       | 0/0 ✓             | 0.2282        | multigraph/STRidge     | —           | 33/0 ✗       | 0.3237   |
| 2%    | weak (IBP)    | spectral   | Radial decay  | —       | 8/0 ✗             | 0.5533        | —                      | —           | 45/0 ✗       | 0.3933   |
| 2%    | weak (IBP)    | spectral   | Radial growth | —       | 7/0 ✗             | 0.2583        | —                      | —           | 30/0 ✗       | 0.3394   |
| 2%    | weak (IBP)    | spectral   | Axial decay   | —       | 8/0 ✗             | 0.5588        | —                      | —           | 35/0 ✗       | 0.3667   |
| 2%    | weak (IBP)    | spectral   | Axial growth  | —       | 7/0 ✗             | 0.2617        | —                      | —           | 27/0 ✗       | 0.3900   |
| 2%    | weak (IBP)    | polynomial | Radial decay  | —       | 3/0 ✗             | 0.3579        | —                      | —           | 54/0 ✗       | 0.3689   |
| 2%    | weak (IBP)    | polynomial | Radial growth | 🟢       | 0/0 ✓             | 0.2644        | mod_multigraph/STRidge | —           | 31/0 ✗       | 0.2776   |
| 2%    | weak (IBP)    | polynomial | Axial decay   | —       | 3/0 ✗             | 0.3163        | —                      | —           | 32/0 ✗       | 0.3739   |
| 2%    | weak (IBP)    | polynomial | Axial growth  | —       | 12/0 ✗            | 0.3806        | —                      | —           | 27/0 ✗       | 0.2996   |
| 5%    | weak (no IBP) | spectral   | Radial decay  | —       | 37/0 ✗            | 0.5184        | —                      | —           | 65/0 ✗       | 0.5411   |
| 5%    | weak (no IBP) | spectral   | Radial growth | —       | 36/0 ✗            | 0.6519        | —                      | —           | 72/0 ✗       | 0.6484   |
| 5%    | weak (no IBP) | spectral   | Axial decay   | —       | 31/0 ✗            | 0.5054        | —                      | —           | 98/0 ✗       | 0.4334   |
| 5%    | weak (no IBP) | spectral   | Axial growth  | —       | 27/0 ✗            | 0.6810        | —                      | —           | 58/0 ✗       | 0.5830   |
| 5%    | weak (no IBP) | polynomial | Radial decay  | —       | 16/0 ✗            | 0.4287        | —                      | —           | 75/0 ✗       | 0.4034   |
| 5%    | weak (no IBP) | polynomial | Radial growth | —       | 17/0 ✗            | 0.4564        | —                      | —           | 38/0 ✗       | 0.3767   |
| 5%    | weak (no IBP) | polynomial | Axial decay   | —       | 26/0 ✗            | 0.4054        | —                      | —           | 40/0 ✗       | 0.3984   |
| 5%    | weak (no IBP) | polynomial | Axial growth  | —       | 37/0 ✗            | 0.4425        | —                      | —           | 31/0 ✗       | 0.3652   |
| 5%    | weak (IBP)    | spectral   | Radial decay  | —       | 32/0 ✗            | 0.5545        | —                      | —           | 69/0 ✗       | 0.5406   |
| 5%    | weak (IBP)    | spectral   | Radial growth | —       | 28/0 ✗            | 0.6569        | —                      | —           | 72/0 ✗       | 0.6608   |
| 5%    | weak (IBP)    | spectral   | Axial decay   | —       | 31/0 ✗            | 0.5034        | —                      | —           | 51/0 ✗       | 0.5276   |
| 5%    | weak (IBP)    | spectral   | Axial growth  | —       | 27/0 ✗            | 0.6901        | —                      | —           | 55/0 ✗       | 0.4719   |
| 5%    | weak (IBP)    | polynomial | Radial decay  | —       | 15/0 ✗            | 0.4054        | —                      | —           | 71/0 ✗       | 0.4043   |
| 5%    | weak (IBP)    | polynomial | Radial growth | —       | 15/0 ✗            | 0.4562        | —                      | —           | 38/0 ✗       | 0.3850   |
| 5%    | weak (IBP)    | polynomial | Axial decay   | —       | 26/0 ✗            | 0.4037        | —                      | —           | 38/0 ✗       | 0.3935   |
| 5%    | weak (IBP)    | polynomial | Axial growth  | —       | 15/0 ✗            | 0.4671        | —                      | —           | 31/0 ✗       | 0.3714   |

## Ablation

### 2% — Axial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 3/0   | 0.3080  | 0.2      | ✗        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=8       | 8/0   | 0.5588  | 0.5      | ✗        |
| 3    | PolyJet + STRidge              | polynomial, weak (IBP), m=8     | 32/0  | 0.3739  | 1.5      | ✗        |
| 4    | Multigraphs + RFE              | spectral, weak (no IBP), m=8    | 0/5   | 0.8460  | 0.1      | ✗        |
| 5    | Modified Multigraphs + RFE     | spectral, weak (no IBP), m=8    | 0/5   | 0.8460  | 0.1      | ✗        |
| 6    | PolyJet + RFE                  | polynomial, weak (no IBP), m=8  | 1/7   | 0.6084  | 0.1      | ✗        |

### 2% — Radial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + STRidge          | polynomial, weak (IBP), m=12    | 3/0   | 0.3579  | 0.2      | ✗        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=8       | 8/0   | 0.5533  | 0.4      | ✗        |
| 3    | PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 43/0  | 0.3620  | 1.2      | ✗        |
| 4    | Modified Multigraphs + RFE     | spectral, weak (no IBP), m=8    | 0/5   | 0.8982  | 0.1      | ✗        |
| 5    | Multigraphs + RFE              | spectral, weak (no IBP), m=8    | 4/5   | 0.5314  | 0.1      | ✗        |
| 6    | PolyJet + RFE                  | polynomial, weak (no IBP), m=8  | 1/7   | 0.6084  | 0.1      | ✗        |

### 2% — Axial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.2282  | 0.2      | —        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (no IBP), m=12   | 8/0   | 0.3579  | 0.4      | ✗        |
| 3    | PolyJet + STRidge              | polynomial, weak (IBP), m=12    | 27/0  | 0.2996  | 1.0      | ✗        |
| 4    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5   | 0.5048  | 0.1      | ✗        |
| 5    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/5   | 0.5048  | 0.1      | ✗        |
| 6    | PolyJet + RFE                  | polynomial, weak (no IBP), m=12 | 0/6   | 0.6848  | 0.1      | ✗        |

### 2% — Radial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.2644  | 0.3      | —        |
| 2    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 3/0   | 0.2193  | 0.2      | ✗        |
| 3    | PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 30/0  | 0.2700  | 1.1      | ✗        |
| 4    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/5   | 0.5048  | 0.1      | ✗        |
| 5    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5   | 0.5048  | 0.1      | ✗        |
| 6    | PolyJet + RFE                  | polynomial, weak (no IBP), m=12 | 0/6   | 0.6851  | 0.1      | ✗        |

### 5% — Axial decay

| Rank | Library + Algorithm            | Configuration                  | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + STRidge          | polynomial, weak (IBP), m=8    | 26/0  | 0.4037  | 0.3      | ✗        |
| 2    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=8 | 34/0  | 0.4057  | 0.3      | ✗        |
| 3    | PolyJet + STRidge              | polynomial, weak (IBP), m=12   | 38/0  | 0.3935  | 1.3      | ✗        |
| 4    | PolyJet + RFE                  | spectral, weak (no IBP), m=8   | 0/6   | 0.9727  | 0.1      | ✗        |
| 5    | Multigraphs + RFE              | polynomial, weak (no IBP), m=8 | 1/7   | 0.6084  | 0.1      | ✗        |
| 6    | Modified Multigraphs + RFE     | spectral, weak (no IBP), m=8   | 1/7   | 0.6084  | 0.1      | ✗        |

### 5% — Radial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + STRidge          | polynomial, weak (IBP), m=8     | 15/0  | 0.4054  | 0.2      | ✗        |
| 2    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 16/0  | 0.4287  | 0.2      | ✗        |
| 3    | PolyJet + STRidge              | spectral, weak (no IBP), m=12   | 65/0  | 0.5411  | 2.1      | ✗        |
| 4    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=8  | 0/5   | 0.7979  | 0.1      | ✗        |
| 5    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 4/5   | 0.5306  | 0.1      | ✗        |
| 6    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 1/7   | 0.6084  | 0.1      | ✗        |

### 5% — Axial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + STRidge          | polynomial, weak (IBP), m=12    | 15/0  | 0.4671  | 0.2      | ✗        |
| 2    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=8     | 28/0  | 0.4713  | 0.4      | ✗        |
| 3    | PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 31/0  | 0.3652  | 1.2      | ✗        |
| 4    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=8  | 0/5   | 0.6280  | 0.1      | ✗        |
| 5    | Multigraphs + RFE              | polynomial, weak (no IBP), m=8  | 0/5   | 0.6280  | 0.1      | ✗        |
| 6    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/6   | 0.9528  | 0.1      | ✗        |

### 5% — Radial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + STRidge          | polynomial, weak (IBP), m=12    | 15/0  | 0.4562  | 0.2      | ✗        |
| 2    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=8     | 27/0  | 0.4806  | 0.4      | ✗        |
| 3    | PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 38/0  | 0.3767  | 1.3      | ✗        |
| 4    | Modified Multigraphs + RFE     | spectral, weak (no IBP), m=8    | 0/5   | 0.9526  | 0.0      | ✗        |
| 5    | Multigraphs + RFE              | polynomial, weak (no IBP), m=8  | 1/7   | 0.6084  | 0.1      | ✗        |
| 6    | PolyJet + RFE                  | polynomial, weak (no IBP), m=8  | 1/7   | 0.6084  | 0.1      | ✗        |


## Detailed results

### weak (no IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/6   | 1        | 0.6848  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/6   | 1        | 0.6851  | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 38/0  | 45       | 0.3780  | 1.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 43/0  | 50       | 0.3620  | 1.2      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 33/0  | 40       | 0.3237  | 1.0      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial growth | 30/0  | 37       | 0.2700  | 1.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 40/0  | 47       | 0.3984  | 1.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 75/0  | 82       | 0.4034  | 2.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 31/0  | 38       | 0.3652  | 1.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Radial growth | 38/0  | 45       | 0.3767  | 1.3      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/5   | 1        | 0.5048  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/5   | 1        | 0.5048  | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 3/0   | 4        | 0.3080  | 0.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 4/0   | 5        | 0.3258  | 0.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 0/0   | 3        | 0.2282  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial growth | 3/0   | 6        | 0.2193  | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/5   | 1        | 0.7979  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/5   | 1        | 0.6280  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 1.0       | Axial decay   | 26/0  | 15       | 0.4054  | 0.3      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial decay  | 16/0  | 9        | 0.4287  | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 0.001     | Axial growth  | 37/0  | 20       | 0.4425  | 0.3      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial growth | 17/0  | 11       | 0.4564  | 0.2      | 🔴      |

### weak (no IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.9994  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.9996  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial growth  | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 1.0       | Axial decay   | 46/0  | 53       | 0.3665  | 2.0      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 56/0  | 63       | 0.4058  | 1.8      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 1.0       | Axial growth  | 45/0  | 52       | 0.3680  | 1.8      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 1.0       | Radial growth | 59/0  | 66       | 0.3877  | 2.2      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/6   | 1        | 0.9727  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/6   | 1        | 0.9528  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.9998  | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.001     | Axial decay   | 98/0  | 105      | 0.4334  | 4.5      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 65/0  | 72       | 0.5411  | 2.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 58/0  | 65       | 0.5830  | 1.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial growth | 72/0  | 79       | 0.6484  | 2.0      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/5   | 1        | 0.8460  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial decay  | 0/5   | 1        | 0.8982  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial growth  | 0/5   | 1        | 0.7593  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial growth | 0/5   | 1        | 0.7550  | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 14/0  | 11       | 0.3518  | 0.4      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 23/0  | 13       | 0.4001  | 0.4      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 8/0   | 7        | 0.3579  | 0.4      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial growth | 9/0   | 6        | 0.2554  | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.6084  | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/7   | 1        | 0.6084  | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/5   | 1        | 0.9526  | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 1.0       | Axial decay   | 31/0  | 17       | 0.5054  | 0.5      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 0.001     | Radial decay  | 37/0  | 20       | 0.5184  | 0.6      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 27/0  | 15       | 0.6810  | 0.4      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial growth | 36/0  | 16       | 0.6519  | 0.3      | 🔴      |

### weak (IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/6   | 1        | 0.6848  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/6   | 1        | 0.6851  | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 1.0       | Axial decay   | 32/0  | 39       | 0.3739  | 1.5      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 54/0  | 61       | 0.3689  | 1.4      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 27/0  | 34       | 0.2996  | 1.0      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial growth | 31/0  | 38       | 0.2776  | 1.0      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 38/0  | 45       | 0.3935  | 1.3      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 71/0  | 78       | 0.4043  | 2.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 31/0  | 38       | 0.3714  | 1.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Radial growth | 38/0  | 45       | 0.3850  | 1.3      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/5   | 1        | 0.5048  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/5   | 1        | 0.5048  | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 3/0   | 4        | 0.3163  | 0.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 3/0   | 4        | 0.3579  | 0.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 10.0      | Axial growth  | 12/0  | 6        | 0.3806  | 0.4      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial growth | 0/0   | 3        | 0.2644  | 0.3      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 4/5   | 3        | 0.5306  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/5   | 1        | 0.6280  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 1.0       | Axial decay   | 26/0  | 15       | 0.4037  | 0.3      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 15/0  | 11       | 0.4054  | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 15/0  | 11       | 0.4671  | 0.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial growth | 15/0  | 11       | 0.4562  | 0.2      | 🔴      |

### weak (IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.9994  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.9996  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial growth  | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Axial decay   | 35/0  | 42       | 0.3667  | 1.4      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 45/0  | 52       | 0.3933  | 1.3      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 1.0       | Axial growth  | 27/0  | 34       | 0.3900  | 1.5      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial growth | 30/0  | 37       | 0.3394  | 1.0      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/6   | 1        | 0.9727  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/6   | 1        | 0.9528  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.9998  | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 51/0  | 58       | 0.5276  | 1.9      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 69/0  | 76       | 0.5406  | 2.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Axial growth  | 55/0  | 62       | 0.4719  | 2.9      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial growth | 72/0  | 79       | 0.6608  | 2.0      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/5   | 1        | 0.8460  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial decay  | 0/5   | 1        | 0.8982  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial growth  | 0/5   | 1        | 0.7593  | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial growth | 0/5   | 1        | 0.7550  | 0.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 1.0       | Axial decay   | 8/0   | 9        | 0.5588  | 0.5      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 1.0       | Radial decay  | 8/0   | 9        | 0.5533  | 0.4      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 7/0   | 5        | 0.2617  | 0.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial growth | 7/0   | 5        | 0.2583  | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 1/7   | 1        | 0.6084  | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial decay  | 1/7   | 1        | 0.6084  | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/5   | 1        | 0.9485  | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/7   | 1        | 0.6084  | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 1.0       | Axial decay   | 31/0  | 17       | 0.5034  | 0.5      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 32/0  | 18       | 0.5545  | 0.4      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 27/0  | 15       | 0.6901  | 0.4      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial growth | 28/0  | 16       | 0.6569  | 0.4      | 🔴      |


---
*Source: `sweep_aniso_default_ks_cell_20260805_215410.csv`*

