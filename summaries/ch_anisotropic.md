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

Sweep configuration: 3 libraries (Multigraphs, PolyJet, WSINDy),
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
| 2%    | weak (no IBP) | polynomial | Radial decay  | 🟥       | 0/0 ✓             | 0.4482        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0796   |
| 2%    | weak (no IBP) | polynomial | Radial growth | 🟥       | 0/0 ✓             | 0.4260        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0602   |
| 2%    | weak (no IBP) | polynomial | Axial decay   | 🟥       | 0/0 ✓             | 0.4385        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0719   |
| 2%    | weak (no IBP) | polynomial | Axial growth  | 🟥       | 0/0 ✓             | 0.4264        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0602   |
| 2%    | weak (IBP)    | spectral   | Radial decay  | 🟩       | 0/0 ✓             | 0.0608        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0735   |
| 2%    | weak (IBP)    | spectral   | Radial growth | 🟩       | 0/0 ✓             | 0.0444        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0523   |
| 2%    | weak (IBP)    | spectral   | Axial decay   | 🟩       | 0/0 ✓             | 0.0564        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0642   |
| 2%    | weak (IBP)    | spectral   | Axial growth  | 🟩       | 0/0 ✓             | 0.0434        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0528   |
| 2%    | weak (IBP)    | polynomial | Radial decay  | 🟩       | 0/0 ✓             | 0.0608        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.1131   |
| 2%    | weak (IBP)    | polynomial | Radial growth | 🟩       | 0/0 ✓             | 0.0444        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0958   |
| 2%    | weak (IBP)    | polynomial | Axial decay   | 🟩       | 0/0 ✓             | 0.0564        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.1070   |
| 2%    | weak (IBP)    | polynomial | Axial growth  | 🟩       | 0/0 ✓             | 0.0434        | mod_multigraph/RFE | WSINDy/RFE  | 0/0 ✓        | 0.0954   |
| 5%    | weak (no IBP) | spectral   | Radial decay  | —       | 11/0 ✗            | 1.0909        | —                  | —           | 2/0 ✗        | 0.3518   |
| 5%    | weak (no IBP) | spectral   | Radial growth | —       | 10/0 ✗            | 1.5596        | —                  | —           | 2/0 ✗        | 0.1802   |
| 5%    | weak (no IBP) | spectral   | Axial decay   | —       | 17/0 ✗            | 1.5291        | —                  | —           | 2/0 ✗        | 0.2814   |
| 5%    | weak (no IBP) | spectral   | Axial growth  | —       | 10/0 ✗            | 1.3851        | —                  | —           | 2/0 ✗        | 0.1818   |
| 5%    | weak (no IBP) | polynomial | Radial decay  | 🟢       | 0/0 ✓             | 1.2588        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.3710   |
| 5%    | weak (no IBP) | polynomial | Radial growth | 🟢       | 0/0 ✓             | 0.6208        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.1900   |
| 5%    | weak (no IBP) | polynomial | Axial decay   | 🟢       | 0/0 ✓             | 1.4953        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.3024   |
| 5%    | weak (no IBP) | polynomial | Axial growth  | 🟢       | 0/0 ✓             | 0.6214        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.1932   |
| 5%    | weak (IBP)    | spectral   | Radial decay  | 🟢       | 0/0 ✓             | 0.2354        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.3589   |
| 5%    | weak (IBP)    | spectral   | Radial growth | 🟢       | 0/0 ✓             | 0.1461        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.1835   |
| 5%    | weak (IBP)    | spectral   | Axial decay   | 🟢       | 0/0 ✓             | 0.2032        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.2885   |
| 5%    | weak (IBP)    | spectral   | Axial growth  | 🟢       | 0/0 ✓             | 0.1422        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.1844   |
| 5%    | weak (IBP)    | polynomial | Radial decay  | 🟢       | 0/0 ✓             | 0.2354        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.3934   |
| 5%    | weak (IBP)    | polynomial | Radial growth | 🟢       | 0/0 ✓             | 0.1461        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.2082   |
| 5%    | weak (IBP)    | polynomial | Axial decay   | 🟢       | 0/0 ✓             | 0.2032        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.3270   |
| 5%    | weak (IBP)    | polynomial | Axial growth  | 🟢       | 0/0 ✓             | 0.1422        | mod_multigraph/RFE | —           | 2/0 ✗        | 0.2043   |

## Ablation

### 2% — Axial decay

| Rank | Library + Algorithm            | Configuration                  | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-------------------------------|:------|:--------|:---------|:---------|
| 1    | WSINDy + RFE                   | spectral, weak (no IBP), m=12  | 0/0 ✓ | 0.0539  | 0.7      | —        |
| 2    | WSINDy + STRidge               | spectral, weak (no IBP), m=12  | 0/0 ✓ | 0.0546  | 2.0      | +1.4%    |
| 3    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0564  | 0.7      | +4.8%    |
| 4    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0572  | 1.2      | +6.2%    |
| 5    | Multigraphs + RFE              | spectral, weak (IBP), m=8      | 0/0 ✓ | 0.3736  | 1.0      | +593.6%  |
| 6    | PolyJet + RFE                  | spectral, weak (no IBP), m=8   | 0/0 ✓ | 0.3813  | 1.2      | +607.8%  |
| 7    | PolyJet + STRidge              | spectral, weak (no IBP), m=8   | 4/0   | 1.4887  | 4.7      | ✗        |
| 8    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=8 | 7/0   | 0.5576  | 1.3      | ✗        |

### 2% — Radial decay

| Rank | Library + Algorithm            | Configuration                  | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0608  | 0.7      | —        |
| 2    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12   | 0/0 ✓ | 0.0625  | 2.3      | +2.7%    |
| 3    | WSINDy + RFE                   | spectral, weak (no IBP), m=12  | 0/0 ✓ | 0.0634  | 0.7      | +4.1%    |
| 4    | WSINDy + STRidge               | spectral, weak (no IBP), m=12  | 0/0 ✓ | 0.0642  | 1.8      | +5.5%    |
| 5    | Multigraphs + RFE              | polynomial, weak (no IBP), m=8 | 0/0 ✓ | 0.4843  | 0.8      | +696.0%  |
| 6    | PolyJet + RFE                  | polynomial, weak (no IBP), m=8 | 0/0 ✓ | 0.5540  | 0.8      | +810.5%  |
| 7    | PolyJet + STRidge              | spectral, weak (no IBP), m=12  | 4/0   | 1.5508  | 5.0      | ✗        |
| 8    | Multigraphs + STRidge          | spectral, weak (no IBP), m=12  | 7/0   | 1.4968  | 1.5      | ✗        |

### 2% — Axial growth

| Rank | Library + Algorithm            | Configuration                 | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:------------------------------|:------|:--------|:---------|:---------|
| 1    | WSINDy + RFE                   | spectral, weak (no IBP), m=12 | 0/0 ✓ | 0.0430  | 0.8      | —        |
| 2    | WSINDy + STRidge               | spectral, weak (no IBP), m=12 | 0/0 ✓ | 0.0433  | 2.4      | +0.6%    |
| 3    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12    | 0/0 ✓ | 0.0434  | 0.7      | +0.8%    |
| 4    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12    | 0/0 ✓ | 0.0453  | 1.1      | +5.2%    |
| 5    | Multigraphs + RFE              | spectral, weak (IBP), m=8     | 0/0 ✓ | 0.2763  | 1.0      | +542.1%  |
| 6    | PolyJet + RFE                  | spectral, weak (no IBP), m=8  | 0/0 ✓ | 0.2849  | 0.8      | +562.0%  |
| 7    | PolyJet + STRidge              | polynomial, weak (IBP), m=8   | 2/0   | 0.5810  | 4.0      | ✗        |
| 8    | Multigraphs + STRidge          | spectral, weak (no IBP), m=8  | 5/0   | 0.9498  | 1.4      | ✗        |

### 2% — Radial growth

| Rank | Library + Algorithm            | Configuration                 | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:------------------------------|:------|:--------|:---------|:---------|
| 1    | WSINDy + RFE                   | spectral, weak (no IBP), m=12 | 0/0 ✓ | 0.0424  | 0.7      | —        |
| 2    | WSINDy + STRidge               | spectral, weak (no IBP), m=12 | 0/0 ✓ | 0.0428  | 2.0      | +0.9%    |
| 3    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12    | 0/0 ✓ | 0.0444  | 0.7      | +4.6%    |
| 4    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12  | 0/0 ✓ | 0.0459  | 2.0      | +8.2%    |
| 5    | Multigraphs + RFE              | spectral, weak (IBP), m=8     | 0/0 ✓ | 0.2724  | 0.9      | +542.5%  |
| 6    | PolyJet + RFE                  | spectral, weak (no IBP), m=8  | 0/0 ✓ | 0.2816  | 1.2      | +564.2%  |
| 7    | PolyJet + STRidge              | polynomial, weak (IBP), m=8   | 2/0   | 0.5813  | 4.3      | ✗        |
| 8    | Multigraphs + STRidge          | spectral, weak (no IBP), m=12 | 7/0   | 1.6434  | 1.5      | ✗        |

### 5% — Axial decay

| Rank | Library + Algorithm            | Configuration                  | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.2032  | 0.7      | —        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.2059  | 1.2      | +1.3%    |
| 3    | Multigraphs + RFE              | polynomial, weak (IBP), m=8    | 0/0 ✓ | 0.7256  | 1.3      | +257.2%  |
| 4    | WSINDy + RFE                   | spectral, weak (no IBP), m=12  | 2/0   | 0.2814  | 1.2      | ✗        |
| 5    | PolyJet + STRidge              | spectral, weak (no IBP), m=12  | 6/0   | 1.2944  | 4.6      | ✗        |
| 6    | PolyJet + RFE                  | spectral, weak (no IBP), m=12  | 6/0   | 1.2954  | 0.9      | ✗        |
| 7    | WSINDy + STRidge               | polynomial, weak (IBP), m=12   | 6/0   | 14.2209 | 1.8      | ✗        |
| 8    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=8 | 7/0   | 0.8452  | 1.9      | ✗        |

### 5% — Radial decay

| Rank | Library + Algorithm            | Configuration                  | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.2354  | 0.7      | —        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.2404  | 1.1      | +2.1%    |
| 3    | WSINDy + RFE                   | spectral, weak (no IBP), m=12  | 2/0   | 0.3518  | 1.3      | ✗        |
| 4    | WSINDy + STRidge               | spectral, weak (no IBP), m=8   | 6/0   | 6.1287  | 4.2      | ✗        |
| 5    | Multigraphs + STRidge          | polynomial, weak (IBP), m=8    | 7/0   | 0.9068  | 1.4      | ✗        |
| 6    | PolyJet + RFE                  | spectral, weak (no IBP), m=12  | 8/0   | 1.1595  | 0.9      | ✗        |
| 7    | PolyJet + STRidge              | polynomial, weak (IBP), m=8    | 8/0   | 2.1081  | 4.8      | ✗        |
| 8    | Multigraphs + RFE              | polynomial, weak (no IBP), m=8 | 4/3   | 1.0553  | 0.7      | ✗        |

### 5% — Axial growth

| Rank | Library + Algorithm            | Configuration                  | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.1422  | 0.7      | —        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.1468  | 1.2      | +3.2%    |
| 3    | Multigraphs + RFE              | polynomial, weak (no IBP), m=8 | 0/0 ✓ | 0.6759  | 0.8      | +375.4%  |
| 4    | WSINDy + STRidge               | spectral, weak (no IBP), m=12  | 2/0   | 0.1818  | 3.4      | ✗        |
| 5    | WSINDy + RFE                   | spectral, weak (no IBP), m=12  | 2/0   | 0.1842  | 1.2      | ✗        |
| 6    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=8 | 7/0   | 0.7993  | 1.7      | ✗        |
| 7    | PolyJet + STRidge              | spectral, weak (no IBP), m=8   | 8/0   | 1.5634  | 4.6      | ✗        |
| 8    | PolyJet + RFE                  | spectral, weak (no IBP), m=8   | 8/0   | 1.5697  | 0.9      | ✗        |

### 5% — Radial growth

| Rank | Library + Algorithm            | Configuration                  | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.1461  | 0.7      | —        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.1501  | 1.1      | +2.7%    |
| 3    | Multigraphs + RFE              | polynomial, weak (no IBP), m=8 | 0/0 ✓ | 0.6763  | 0.8      | +362.9%  |
| 4    | WSINDy + STRidge               | spectral, weak (no IBP), m=12  | 2/0   | 0.1802  | 4.2      | ✗        |
| 5    | WSINDy + RFE                   | spectral, weak (no IBP), m=12  | 2/0   | 0.1820  | 1.1      | ✗        |
| 6    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=8 | 7/0   | 0.7972  | 1.5      | ✗        |
| 7    | PolyJet + STRidge              | spectral, weak (no IBP), m=8   | 8/0   | 1.5818  | 4.5      | ✗        |
| 8    | PolyJet + RFE                  | spectral, weak (no IBP), m=12  | 8/0   | 1.5883  | 0.9      | ✗        |


## Detailed results

### weak (no IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0719  | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0796  | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0602  | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0602  | 0.7      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 7        | 0.0727  | 1.8      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 7        | 0.0808  | 1.9      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0612  | 1.9      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0611  | 1.9      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 2/0   | 9        | 0.3024  | 0.7      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 2/0   | 9        | 0.3710  | 0.7      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 2/0   | 9        | 0.1932  | 0.8      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 2/0   | 9        | 0.1936  | 0.7      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 6/0   | 10       | 15.3829 | 1.9      | 🔴      |
| WSINDy      | 5%    | STRidge   | 8  | 10.0      | Radial decay  | 10/0  | 11       | 11.3127 | 1.8      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 6/0   | 10       | 13.2071 | 1.9      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 2/0   | 9        | 0.1900  | 2.0      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 9        | 0.5446  | 0.8      | 🟢      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 0/0   | 9        | 0.5540  | 0.8      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 2/0   | 11       | 0.6341  | 0.8      | 🔴      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 2/0   | 11       | 0.6375  | 1.2      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Axial decay   | 9/0   | 18       | 2.3601  | 4.4      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 6/0   | 15       | 1.2140  | 4.5      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 9/0   | 18       | 1.7944  | 4.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial growth | 9/0   | 18       | 1.7930  | 4.5      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/5   | 4        | 2.8680  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/5   | 4        | 2.5329  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/3   | 7        | 1.0184  | 1.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/3   | 7        | 1.0184  | 0.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Axial decay   | 8/0   | 17       | 2.1901  | 4.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Radial decay  | 8/0   | 17       | 2.1378  | 4.9      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Axial growth  | 8/0   | 17       | 2.2860  | 4.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Radial growth | 8/0   | 17       | 2.2833  | 4.8      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 3        | 0.4385  | 1.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial decay  | 0/0   | 3        | 0.4482  | 1.3      | 🟢      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial growth  | 0/0   | 3        | 0.4264  | 1.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial growth | 0/0   | 3        | 0.4260  | 1.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 8  | 10.0      | Axial decay   | 7/0   | 6        | 0.5576  | 1.3      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 10/0  | 9        | 0.4120  | 2.4      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 15/0  | 11       | 1.9105  | 2.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial growth | 15/0  | 11       | 1.8929  | 2.4      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 1.4953  | 0.7      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 1.2588  | 0.9      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/0   | 3        | 0.6214  | 0.7      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/0   | 3        | 0.6208  | 0.7      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Axial decay   | 7/0   | 6        | 0.8452  | 1.9      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Radial decay  | 7/0   | 5        | 0.9105  | 1.4      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Axial growth  | 7/0   | 6        | 0.7993  | 1.7      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Radial growth | 7/0   | 6        | 0.7972  | 1.5      | 🔴      |

### weak (no IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0539  | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0634  | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0430  | 0.8      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0424  | 0.7      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 7        | 0.0546  | 2.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 7        | 0.0642  | 1.8      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0433  | 2.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0428  | 2.0      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 2/0   | 9        | 0.2814  | 1.2      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 2/0   | 9        | 0.3518  | 1.3      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 2/0   | 9        | 0.1842  | 1.2      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 2/0   | 9        | 0.1820  | 1.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 6/0   | 10       | 14.9189 | 4.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 8  | 10.0      | Radial decay  | 6/0   | 10       | 6.1287  | 4.2      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 2/0   | 9        | 0.1818  | 3.4      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Radial growth | 2/0   | 9        | 0.1802  | 4.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 9        | 0.3813  | 1.2      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 6/0   | 15       | 1.3572  | 1.2      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial growth  | 0/0   | 9        | 0.2849  | 0.8      | 🟢      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial growth | 0/0   | 9        | 0.2816  | 1.2      | 🟢      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Axial decay   | 4/0   | 13       | 1.4887  | 4.7      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 4/0   | 13       | 1.5508  | 5.0      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 6/0   | 15       | 1.4234  | 4.0      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial growth | 5/0   | 14       | 1.3149  | 4.5      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 6/0   | 15       | 1.2954  | 0.9      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 8/0   | 17       | 1.1595  | 0.9      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 8/0   | 17       | 1.5697  | 0.9      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial growth | 8/0   | 17       | 1.5883  | 0.9      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 6/0   | 15       | 1.2944  | 4.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 10/0  | 19       | 1.3061  | 6.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Axial growth  | 8/0   | 17       | 1.5634  | 4.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Radial growth | 8/0   | 17       | 1.5818  | 4.5      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 4        | 0.3769  | 1.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 7/0   | 6        | 1.4951  | 0.8      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial growth  | 0/0   | 4        | 0.2798  | 1.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial growth | 0/0   | 4        | 0.2760  | 0.9      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 7/0   | 6        | 1.5840  | 1.6      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 10.0      | Radial decay  | 7/0   | 8        | 1.2540  | 1.2      | 🔴      |
| Multigraphs | 2%    | STRidge   | 8  | 10.0      | Axial growth  | 5/0   | 6        | 0.9498  | 1.4      | 🔴      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial growth | 7/0   | 6        | 1.6434  | 1.5      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 18/0  | 10       | 1.4716  | 0.7      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial decay  | 2/7   | 3        | 0.9787  | 0.7      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial growth  | 10/0  | 7        | 1.5281  | 0.6      | 🔴      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial growth | 10/0  | 7        | 1.5596  | 0.6      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 17/0  | 11       | 1.5291  | 1.6      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 11/0  | 10       | 1.0909  | 1.4      | 🔴      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 10/0  | 7        | 1.3851  | 1.8      | 🔴      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Radial growth | 11/0  | 8        | 1.5555  | 1.4      | 🔴      |

### weak (IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.1070  | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.1131  | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0954  | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0958  | 0.7      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 7        | 0.1076  | 1.8      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 7        | 0.1143  | 1.7      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0965  | 1.9      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0967  | 1.9      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 2/0   | 9        | 0.3270  | 0.7      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 2/0   | 9        | 0.3934  | 0.7      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 2/0   | 9        | 0.2109  | 0.7      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 2/0   | 9        | 0.2135  | 0.7      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 6/0   | 10       | 14.2209 | 1.8      | 🔴      |
| WSINDy      | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 10/0  | 11       | 7.7532  | 2.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 2/0   | 9        | 0.2043  | 1.8      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Radial growth | 2/0   | 9        | 0.2082  | 1.8      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 9        | 0.5581  | 0.8      | 🟢      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 2/0   | 11       | 0.5797  | 0.8      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial growth  | 2/0   | 11       | 0.5726  | 0.8      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial growth | 2/0   | 11       | 0.5730  | 1.2      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Axial decay   | 9/0   | 18       | 3.3146  | 4.3      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Radial decay  | 9/0   | 18       | 3.2965  | 4.7      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Axial growth  | 2/0   | 11       | 0.5810  | 4.0      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Radial growth | 2/0   | 11       | 0.5813  | 4.3      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/5   | 4        | 2.8680  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/5   | 4        | 2.5329  | 0.8      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 1/3   | 7        | 1.0088  | 1.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 1/3   | 7        | 1.0088  | 0.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Axial decay   | 8/0   | 17       | 2.1509  | 4.3      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Radial decay  | 8/0   | 17       | 2.1081  | 4.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Axial growth  | 8/0   | 17       | 2.2300  | 4.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Radial growth | 8/0   | 17       | 2.2275  | 5.0      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.0564  | 1.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.0608  | 1.4      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.0434  | 1.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.0444  | 1.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 3        | 0.0572  | 2.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 3        | 0.0625  | 2.3      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 0/0   | 3        | 0.0453  | 2.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial growth | 0/0   | 3        | 0.0459  | 2.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.2032  | 1.2      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.2354  | 0.9      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.1422  | 0.7      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.1461  | 1.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 4        | 0.4325  | 1.8      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 4        | 0.5277  | 1.4      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 0/0   | 3        | 0.1468  | 1.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Radial growth | 7/0   | 6        | 0.8083  | 1.4      | 🔴      |

### weak (IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0642  | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0735  | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0528  | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0523  | 0.7      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 7        | 0.0649  | 1.9      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 7        | 0.0743  | 1.9      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0529  | 2.3      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0526  | 1.9      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 2/0   | 9        | 0.2885  | 1.2      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 2/0   | 9        | 0.3589  | 1.3      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 2/0   | 9        | 0.1877  | 1.5      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 2/0   | 9        | 0.1860  | 1.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 6/0   | 10       | 14.4017 | 4.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 10/0  | 11       | 8.9799  | 4.5      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 2/0   | 9        | 0.1844  | 3.7      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Radial growth | 2/0   | 9        | 0.1835  | 3.6      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 9        | 0.3819  | 1.2      | 🟢      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 2/0   | 11       | 1.2264  | 1.3      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial growth  | 2/0   | 11       | 0.4381  | 0.8      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial growth | 2/0   | 11       | 0.4330  | 1.2      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 4/0   | 13       | 1.5120  | 4.4      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Radial decay  | 4/0   | 13       | 1.5578  | 5.3      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 10.0      | Axial growth  | 4/0   | 13       | 1.4676  | 4.0      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial growth | 5/0   | 14       | 1.3205  | 4.4      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 6/0   | 15       | 1.3011  | 0.9      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 8/0   | 17       | 1.1637  | 0.9      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial growth  | 8/0   | 17       | 1.5815  | 0.9      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial growth | 8/0   | 17       | 1.5976  | 0.9      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 6/0   | 15       | 1.3001  | 4.8      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 12/0  | 21       | 1.4241  | 6.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Axial growth  | 8/0   | 17       | 1.5850  | 4.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 10.0      | Radial growth | 8/0   | 17       | 1.6066  | 4.4      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.0564  | 0.7      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.0608  | 0.7      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.0434  | 0.7      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.0444  | 0.7      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 3        | 0.0572  | 1.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 3        | 0.0625  | 1.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 0/0   | 3        | 0.0453  | 1.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial growth | 0/0   | 3        | 0.0459  | 1.2      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.2032  | 0.7      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.2354  | 0.7      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.1422  | 0.7      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.1461  | 0.7      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 3        | 0.2059  | 1.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 3        | 0.2404  | 1.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial growth  | 0/0   | 3        | 0.1468  | 1.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial growth | 0/0   | 3        | 0.1501  | 1.1      | 🟢      |


---
*Source: `sweep_aniso_default_ch_20260805_215410.csv`*

