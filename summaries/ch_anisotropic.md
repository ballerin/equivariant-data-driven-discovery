# Cahn-Hilliard equation — Anisotropic noise

The equation under consideration is

$$f_t = -\Delta^2 f + \Delta(f^3 - f) = -\Delta f + 3f^2\Delta f + 6f\|\nabla f\|^2 - \varepsilon^2\Delta^2 f, \qquad \varepsilon = 0.2$$

which can be written in terms of $O(2)$-equivariant terms

| Operator               | True coefficient         |
|:-----------------------|:-------------------------|
| $\Delta f$             | $-1.0$                   |
| $f^2\cdot\Delta f$     | $+3.0$                   |
| $f\cdot\|\nabla f\|^2$ | $+6.0$                   |
| $\Delta^2 f$           | $-\varepsilon^2 = -0.04$ |

or standard terms (PolyJet basis)

| Operator                  | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|:--------------------------|:-----------|:---------|:-------|:---------|:-----------|
| $\partial_{xx} f$         |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yy} f$         |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f^2\cdot\partial_{xx} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f^2\cdot\partial_{yy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f\cdot f_x^2$            |            | 🔴        | 🟢      | 🟢        | 🟢          |
| $f\cdot f_y^2$            |            | 🔴        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxxx} f$       |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxyy} f$       |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yyyy} f$       |            | 🟢        | 🟢      | 🟢        | 🟢          |

All libraries are evaluated in the common polynomial-jet basis (max order 4,
max degree 3).  fp/fn/E are computed against the same ground truth for every
library.  Cahn-Hilliard combines second and fourth-order terms with nonlinear products.  PDE-FIND cannot represent $f\cdot f_x^2$ (products of derivatives).  WSINDy cannot represent either $f^2\cdot f_{xx}$ or $f\cdot f_x^2$ individually, but the combination $3f^2\Delta f + 6f\|\nabla f\|^2$ equals $\Delta(f^3)$, which IS representable as $(f^3)_{xx} + (f^3)_{yy}$ in the WSINDy basis.

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

Sweep configuration: 4 libraries (Multigraphs, PolyJet, PDE-FIND, WSINDy),
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


| Noise | Method        | Backend    | Profile       | Eq best | Multigraphs fp/fn | Multigraphs E | Mult. best         | Non-eq best | Non-eq fp/fn | Non-eq E |
|:------|:-------------:|:----------:|:-------------:|:-------:|:------------------|:--------------|:-------------------|:------------|:-------------|:---------|
| 2%    | weak (no IBP) | spectral   | Radial decay  | 🔴       | 7/0 ✗             | 1.2540        | —                  | WSINDy/RFE  | 0/0 ✓        | 0.0634   |
| 2%    | weak (no IBP) | spectral   | Radial growth | 🟥       | 0/0 ✓             | 0.2760        | multigraph/RFE     | WSINDy/RFE  | 0/0 ✓        | 0.0424   |
| 2%    | weak (no IBP) | spectral   | Axial decay   | 🟥       | 0/0 ✓             | 0.3769        | multigraph/RFE     | WSINDy/RFE  | 0/0 ✓        | 0.0539   |
| 2%    | weak (no IBP) | spectral   | Axial growth  | 🟥       | 0/0 ✓             | 0.2798        | multigraph/RFE     | WSINDy/RFE  | 0/0 ✓        | 0.0430   |
| 2%    | weak (no IBP) | polynomial | Radial decay  | 🟥       | 0/0 ✓             | 0.2500        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0424   |
| 2%    | weak (no IBP) | polynomial | Radial growth | 🟥       | 0/0 ✓             | 0.2253        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0298   |
| 2%    | weak (no IBP) | polynomial | Axial decay   | 🟥       | 0/0 ✓             | 0.2383        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0369   |
| 2%    | weak (no IBP) | polynomial | Axial growth  | 🟥       | 0/0 ✓             | 0.2258        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0299   |
| 2%    | weak (IBP)    | spectral   | Radial decay  | 🟩       | 0/0 ✓             | 0.0608        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0735   |
| 2%    | weak (IBP)    | spectral   | Radial growth | 🟩       | 0/0 ✓             | 0.0444        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0523   |
| 2%    | weak (IBP)    | spectral   | Axial decay   | 🟩       | 0/0 ✓             | 0.0564        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0642   |
| 2%    | weak (IBP)    | spectral   | Axial growth  | 🟩       | 0/0 ✓             | 0.0434        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0528   |
| 2%    | weak (IBP)    | polynomial | Radial decay  | 🟩       | 0/0 ✓             | 0.0423        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0430   |
| 2%    | weak (IBP)    | polynomial | Radial growth | 🟩       | 0/0 ✓             | 0.0298        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0305   |
| 2%    | weak (IBP)    | polynomial | Axial decay   | 🟩       | 0/0 ✓             | 0.0370        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0376   |
| 2%    | weak (IBP)    | polynomial | Axial growth  | 🟩       | 0/0 ✓             | 0.0299        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0305   |
| 5%    | weak (no IBP) | spectral   | Radial decay  | —       | 11/0 ✗            | 1.0909        | —                  | —           | 2/0 ✗        | 0.3518   |
| 5%    | weak (no IBP) | spectral   | Radial growth | —       | 10/0 ✗            | 1.4304        | —                  | —           | 2/0 ✗        | 0.1802   |
| 5%    | weak (no IBP) | spectral   | Axial decay   | —       | 17/0 ✗            | 1.5291        | —                  | —           | 2/0 ✗        | 0.2814   |
| 5%    | weak (no IBP) | spectral   | Axial growth  | —       | 10/0 ✗            | 1.3851        | —                  | —           | 2/0 ✗        | 0.1818   |
| 5%    | weak (no IBP) | polynomial | Radial decay  | 🟢       | 0/0 ✓             | 0.5360        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.2683   |
| 5%    | weak (no IBP) | polynomial | Radial growth | 🟢       | 0/0 ✓             | 0.3753        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.1698   |
| 5%    | weak (no IBP) | polynomial | Axial decay   | 🟢       | 0/0 ✓             | 0.4531        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.2241   |
| 5%    | weak (no IBP) | polynomial | Axial growth  | 🟢       | 0/0 ✓             | 0.3793        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.1708   |
| 5%    | weak (IBP)    | spectral   | Radial decay  | 🟢       | 0/0 ✓             | 0.2354        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.3589   |
| 5%    | weak (IBP)    | spectral   | Radial growth | 🟢       | 0/0 ✓             | 0.1461        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.1835   |
| 5%    | weak (IBP)    | spectral   | Axial decay   | 🟢       | 0/0 ✓             | 0.2032        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.2885   |
| 5%    | weak (IBP)    | spectral   | Axial growth  | 🟢       | 0/0 ✓             | 0.1422        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.1844   |
| 5%    | weak (IBP)    | polynomial | Radial decay  | 🟢       | 0/0 ✓             | 0.1916        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.2686   |
| 5%    | weak (IBP)    | polynomial | Radial growth | 🟢       | 0/0 ✓             | 0.1176        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.1701   |
| 5%    | weak (IBP)    | polynomial | Axial decay   | 🟢       | 0/0 ✓             | 0.1569        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.2244   |
| 5%    | weak (IBP)    | polynomial | Axial growth  | 🟢       | 0/0 ✓             | 0.1182        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.1711   |

## Ablation

### 2% — Axial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | WSINDy + RFE                   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0369  | 0.6      | —        |
| 2    | Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0370  | 0.6      | +0.1%    |
| 3    | WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0374  | 1.8      | +1.2%    |
| 4    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0374  | 1.1      | +1.3%    |
| 5    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.3078  | 0.6      | +733.1%  |
| 6    | PolyJet + RFE                  | polynomial, weak (IBP), m=8     | 0/0 ✓ | 0.3569  | 0.8      | +865.9%  |
| 7    | PolyJet + STRidge              | spectral, weak (no IBP), m=8    | 4/0   | 1.4887  | 4.2      | ✗        |
| 8    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=8  | 7/0   | 0.5077  | 1.2      | ✗        |

### 2% — Radial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0423  | 0.6      | —        |
| 2    | WSINDy + RFE                   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0424  | 0.6      | +0.1%    |
| 3    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0427  | 1.1      | +0.9%    |
| 4    | WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0428  | 1.9      | +1.0%    |
| 5    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.3187  | 0.6      | +652.9%  |
| 6    | PolyJet + RFE                  | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.3925  | 0.7      | +827.1%  |
| 7    | PolyJet + STRidge              | spectral, weak (no IBP), m=12   | 4/0   | 1.5508  | 4.2      | ✗        |
| 8    | Multigraphs + STRidge          | spectral, weak (no IBP), m=12   | 7/0   | 1.4968  | 1.3      | ✗        |

### 2% — Axial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | WSINDy + RFE                   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0299  | 0.6      | —        |
| 2    | Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0299  | 0.6      | +0.0%    |
| 3    | WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0302  | 2.0      | +1.0%    |
| 4    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0302  | 1.1      | +1.2%    |
| 5    | Multigraphs + RFE              | spectral, weak (IBP), m=8       | 0/0 ✓ | 0.2763  | 0.6      | +824.7%  |
| 6    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.2849  | 0.8      | +853.3%  |
| 7    | PolyJet + STRidge              | polynomial, weak (IBP), m=8     | 2/0   | 0.4113  | 4.1      | ✗        |
| 8    | Multigraphs + STRidge          | spectral, weak (no IBP), m=8    | 5/0   | 0.9498  | 1.1      | ✗        |

### 2% — Radial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | WSINDy + RFE                   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0298  | 0.6      | —        |
| 2    | Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0298  | 0.6      | +0.1%    |
| 3    | WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0301  | 2.0      | +1.1%    |
| 4    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0302  | 1.1      | +1.4%    |
| 5    | Multigraphs + RFE              | spectral, weak (IBP), m=8       | 0/0 ✓ | 0.2724  | 0.6      | +814.3%  |
| 6    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.2816  | 0.8      | +845.1%  |
| 7    | PolyJet + STRidge              | polynomial, weak (IBP), m=8     | 2/0   | 0.4109  | 4.2      | ✗        |
| 8    | Multigraphs + STRidge          | spectral, weak (no IBP), m=8    | 5/0   | 0.9444  | 1.1      | ✗        |

### 5% — Axial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.1569  | 0.6      | —        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.2059  | 1.1      | +31.2%   |
| 3    | Multigraphs + RFE              | polynomial, weak (IBP), m=8     | 0/0 ✓ | 0.6097  | 0.6      | +288.6%  |
| 4    | WSINDy + RFE                   | polynomial, weak (no IBP), m=12 | 2/0   | 0.2241  | 0.6      | ✗        |
| 5    | WSINDy + STRidge               | polynomial, weak (IBP), m=12    | 6/0   | 0.9690  | 1.9      | ✗        |
| 6    | PolyJet + STRidge              | spectral, weak (no IBP), m=12   | 6/0   | 1.2944  | 4.4      | ✗        |
| 7    | PolyJet + RFE                  | spectral, weak (no IBP), m=12   | 6/0   | 1.2954  | 0.8      | ✗        |
| 8    | Multigraphs + STRidge          | polynomial, weak (IBP), m=8     | 7/0   | 0.7315  | 1.2      | ✗        |

### 5% — Radial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.1916  | 0.6      | —        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.2404  | 1.1      | +25.4%   |
| 3    | WSINDy + RFE                   | polynomial, weak (no IBP), m=12 | 2/0   | 0.2683  | 0.6      | ✗        |
| 4    | WSINDy + STRidge               | spectral, weak (no IBP), m=8    | 6/0   | 6.1287  | 2.0      | ✗        |
| 5    | PolyJet + STRidge              | polynomial, weak (IBP), m=8     | 8/0   | 0.7491  | 4.5      | ✗        |
| 6    | PolyJet + RFE                  | spectral, weak (no IBP), m=12   | 8/0   | 1.1595  | 0.8      | ✗        |
| 7    | Multigraphs + STRidge          | spectral, weak (no IBP), m=8    | 11/0  | 1.0909  | 1.3      | ✗        |
| 8    | Multigraphs + RFE              | polynomial, weak (no IBP), m=8  | 4/3   | 0.9660  | 0.6      | ✗        |

### 5% — Axial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.1182  | 0.6      | —        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.1468  | 1.0      | +24.1%   |
| 3    | Multigraphs + RFE              | polynomial, weak (IBP), m=8     | 0/0 ✓ | 0.5649  | 0.6      | +377.8%  |
| 4    | WSINDy + RFE                   | polynomial, weak (no IBP), m=12 | 2/0   | 0.1708  | 0.6      | ✗        |
| 5    | WSINDy + STRidge               | polynomial, weak (IBP), m=12    | 2/0   | 0.1723  | 1.8      | ✗        |
| 6    | Multigraphs + STRidge          | polynomial, weak (IBP), m=8     | 7/0   | 0.6938  | 1.1      | ✗        |
| 7    | PolyJet + STRidge              | polynomial, weak (IBP), m=8     | 8/0   | 0.6882  | 4.4      | ✗        |
| 8    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 8/0   | 1.5697  | 0.7      | ✗        |

### 5% — Radial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.1176  | 0.6      | —        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.1501  | 1.0      | +27.6%   |
| 3    | Multigraphs + RFE              | polynomial, weak (IBP), m=8     | 0/0 ✓ | 0.5633  | 0.6      | +378.8%  |
| 4    | WSINDy + RFE                   | polynomial, weak (no IBP), m=12 | 2/0   | 0.1698  | 0.6      | ✗        |
| 5    | WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 2/0   | 0.1712  | 2.1      | ✗        |
| 6    | Multigraphs + STRidge          | polynomial, weak (IBP), m=8     | 7/0   | 0.6922  | 1.2      | ✗        |
| 7    | PolyJet + STRidge              | polynomial, weak (IBP), m=8     | 8/0   | 0.6866  | 4.3      | ✗        |
| 8    | PolyJet + RFE                  | spectral, weak (no IBP), m=12   | 8/0   | 1.5883  | 0.8      | ✗        |


## Detailed results

### weak (no IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0369  | 0.6      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0424  | 0.6      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0299  | 0.6      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0298  | 0.6      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 7        | 0.0374  | 1.8      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 7        | 0.0428  | 1.9      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0302  | 2.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0301  | 2.0      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 2/0   | 9        | 0.2241  | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 2/0   | 9        | 0.2683  | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 2/0   | 9        | 0.1708  | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 2/0   | 9        | 0.1698  | 0.6      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 6/0   | 10       | 1.0359  | 2.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 8  | 10.0      | Radial decay  | 10/0  | 11       | 1.1922  | 1.9      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 6/0   | 10       | 0.8706  | 1.9      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 2/0   | 9        | 0.1712  | 2.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 9        | 0.3827  | 0.8      | 🟢      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 0/0   | 9        | 0.3925  | 0.7      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 2/0   | 11       | 0.3502  | 0.8      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 2/0   | 11       | 0.3492  | 0.7      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Axial decay   | 9/0   | 18       | 0.4846  | 4.6      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 6/0   | 15       | 0.3929  | 4.5      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 9/0   | 18       | 0.3855  | 4.4      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial growth | 9/0   | 18       | 0.3852  | 4.3      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/5   | 4        | 0.7630  | 0.7      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/5   | 4        | 0.7909  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/3   | 7        | 0.8343  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/3   | 7        | 0.8343  | 0.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Axial decay   | 8/0   | 17       | 0.7309  | 4.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Radial decay  | 8/0   | 17       | 0.7553  | 4.5      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Axial growth  | 8/0   | 17       | 0.6999  | 4.3      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Radial growth | 8/0   | 17       | 0.6985  | 4.4      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.2383  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.2500  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.2258  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.2253  | 0.6      | 🟢      |
| Multigraphs | 2%    | STRidge   | 8  | 10.0      | Axial decay   | 7/0   | 6        | 0.5077  | 1.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 10/0  | 9        | 0.3762  | 1.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 15/0  | 11       | 0.4209  | 1.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial growth | 15/0  | 11       | 0.4196  | 1.3      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.4531  | 0.6      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.5360  | 0.6      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.3793  | 0.6      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.3753  | 0.6      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Axial decay   | 7/0   | 6        | 0.7361  | 1.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Radial decay  | 7/0   | 5        | 0.7689  | 1.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Axial growth  | 7/0   | 6        | 0.7014  | 1.1      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Radial growth | 7/0   | 6        | 0.7001  | 1.1      | 🔴      |

### weak (no IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0539  | 0.6      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0634  | 0.6      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0430  | 0.6      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0424  | 0.6      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 7        | 0.0546  | 1.8      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 7        | 0.0642  | 1.9      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0433  | 2.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0428  | 2.0      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 2/0   | 9        | 0.2814  | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 2/0   | 9        | 0.3518  | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 2/0   | 9        | 0.1842  | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 2/0   | 9        | 0.1820  | 0.6      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 6/0   | 10       | 14.9189 | 2.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 8  | 10.0      | Radial decay  | 6/0   | 10       | 6.1287  | 2.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 2/0   | 9        | 0.1818  | 2.2      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Radial growth | 2/0   | 9        | 0.1802  | 1.9      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 9        | 0.3813  | 0.8      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 6/0   | 15       | 1.3572  | 0.8      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial growth  | 0/0   | 9        | 0.2849  | 0.8      | 🟢      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial growth | 0/0   | 9        | 0.2816  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Axial decay   | 4/0   | 13       | 1.4887  | 4.2      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 4/0   | 13       | 1.5508  | 4.2      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 6/0   | 15       | 1.4234  | 4.3      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial growth | 5/0   | 14       | 1.3149  | 4.3      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 6/0   | 15       | 1.2954  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 8/0   | 17       | 1.1595  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 8/0   | 17       | 1.5697  | 0.7      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial growth | 8/0   | 17       | 1.5883  | 0.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 6/0   | 15       | 1.2944  | 4.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 10/0  | 19       | 1.3061  | 5.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Axial growth  | 8/0   | 17       | 1.5634  | 4.3      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Radial growth | 8/0   | 17       | 1.5818  | 4.3      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 4        | 0.3769  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 7/0   | 6        | 1.4951  | 0.6      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial growth  | 0/0   | 4        | 0.2798  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial growth | 0/0   | 4        | 0.2760  | 0.6      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 7/0   | 6        | 1.5840  | 1.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 10.0      | Radial decay  | 7/0   | 8        | 1.2540  | 1.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 10.0      | Axial growth  | 5/0   | 6        | 0.9498  | 1.1      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 10.0      | Radial growth | 5/0   | 6        | 0.9444  | 1.1      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 18/0  | 10       | 1.4716  | 0.6      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial decay  | 2/7   | 3        | 0.9787  | 0.6      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial growth  | 10/0  | 7        | 1.5281  | 0.6      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial growth | 10/0  | 7        | 1.5596  | 0.6      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 17/0  | 11       | 1.5291  | 1.5      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 11/0  | 10       | 1.0909  | 1.3      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 10/0  | 7        | 1.3851  | 1.2      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial growth | 10/0  | 7        | 1.4304  | 1.2      | 🔴      |

### weak (IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0376  | 0.6      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0430  | 0.6      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0305  | 0.6      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0305  | 0.6      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 7        | 0.0380  | 1.8      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 7        | 0.0434  | 1.8      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0309  | 2.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0308  | 2.0      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 2/0   | 9        | 0.2244  | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 2/0   | 9        | 0.2686  | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 2/0   | 9        | 0.1711  | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 2/0   | 9        | 0.1701  | 0.6      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 6/0   | 10       | 0.9690  | 1.9      | 🔴      |
| WSINDy      | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 10/0  | 11       | 1.0118  | 2.2      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 2/0   | 9        | 0.1723  | 1.8      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Radial growth | 2/0   | 9        | 0.1714  | 1.8      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 9        | 0.3569  | 0.8      | 🟢      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 2/0   | 11       | 0.4427  | 0.8      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 2/0   | 11       | 0.3502  | 0.8      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 2/0   | 11       | 0.3492  | 0.8      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Axial decay   | 9/0   | 18       | 0.4800  | 4.3      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Radial decay  | 9/0   | 18       | 0.4932  | 4.3      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Axial growth  | 2/0   | 11       | 0.4113  | 4.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Radial growth | 2/0   | 11       | 0.4109  | 4.2      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/5   | 4        | 0.7630  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/5   | 4        | 0.7909  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/3   | 7        | 0.8355  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/3   | 7        | 0.8356  | 0.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Axial decay   | 8/0   | 17       | 0.7225  | 4.5      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Radial decay  | 8/0   | 17       | 0.7491  | 4.5      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Axial growth  | 8/0   | 17       | 0.6882  | 4.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Radial growth | 8/0   | 17       | 0.6866  | 4.3      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.0370  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.0423  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.0299  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.0298  | 0.6      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 3        | 0.0374  | 1.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 3        | 0.0427  | 1.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 0/0   | 3        | 0.0302  | 1.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial growth | 0/0   | 3        | 0.0302  | 1.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.1569  | 0.6      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.1916  | 0.6      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.1182  | 0.6      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.1176  | 0.6      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 4        | 0.2113  | 1.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 4        | 0.2425  | 1.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 0/0   | 4        | 0.1690  | 1.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Radial growth | 7/0   | 6        | 0.6922  | 1.2      | 🔴      |

### weak (IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0642  | 0.6      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0735  | 0.6      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0528  | 0.6      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0523  | 0.6      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 7        | 0.0649  | 1.9      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 7        | 0.0743  | 1.9      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0529  | 2.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0526  | 2.1      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 2/0   | 9        | 0.2885  | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 2/0   | 9        | 0.3589  | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 2/0   | 9        | 0.1877  | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 2/0   | 9        | 0.1860  | 0.6      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 6/0   | 10       | 14.4017 | 2.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 10/0  | 11       | 8.9799  | 2.2      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 2/0   | 9        | 0.1844  | 1.9      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Radial growth | 2/0   | 9        | 0.1835  | 1.9      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 9        | 0.3819  | 0.8      | 🟢      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 2/0   | 11       | 1.2264  | 0.8      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial growth  | 2/0   | 11       | 0.4381  | 0.7      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial growth | 2/0   | 11       | 0.4330  | 0.8      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 4/0   | 13       | 1.5120  | 4.0      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Radial decay  | 4/0   | 13       | 1.5578  | 4.3      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Axial growth  | 4/0   | 13       | 1.4676  | 4.2      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial growth | 5/0   | 14       | 1.3205  | 4.3      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 6/0   | 15       | 1.3011  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 8/0   | 17       | 1.1637  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial growth  | 8/0   | 17       | 1.5815  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial growth | 8/0   | 17       | 1.5976  | 0.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 6/0   | 15       | 1.3001  | 4.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 12/0  | 21       | 1.4241  | 6.9      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Axial growth  | 8/0   | 17       | 1.5850  | 4.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Radial growth | 8/0   | 17       | 1.6066  | 4.4      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.0564  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.0608  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.0434  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.0444  | 0.6      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 3        | 0.0572  | 1.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 3        | 0.0625  | 1.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 0/0   | 3        | 0.0453  | 1.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial growth | 0/0   | 3        | 0.0459  | 1.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.2032  | 0.6      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.2354  | 0.6      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.1422  | 0.6      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.1461  | 0.6      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 3        | 0.2059  | 1.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 3        | 0.2404  | 1.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 0/0   | 3        | 0.1468  | 1.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial growth | 0/0   | 3        | 0.1501  | 1.0      | 🟢      |


---
*Source: `sweep_aniso_20260803_194535_full_local.csv`*

*Run: 775.9 min (Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz), OMP_NUM_THREADS=1*
