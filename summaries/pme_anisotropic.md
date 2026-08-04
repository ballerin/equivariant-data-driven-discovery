# Porous Medium equation ($m = 3$) — Anisotropic noise

The equation under consideration is

$$f_t = \nabla\cdot(f^2\nabla f) = f^2\Delta f + 2f\|\nabla f\|^2 = \tfrac13\Delta(f^3)$$

which can be written in terms of $O(2)$-equivariant terms

| Operator               | True coefficient |
|:-----------------------|:-----------------|
| $f^2\cdot\Delta f$     | $+1.0$           |
| $f\cdot\|\nabla f\|^2$ | $+2.0$           |

or standard terms (PolyJet basis)

| Operator                  | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|:--------------------------|:-----------|:---------|:-------|:---------|:-----------|
| $f^2\cdot\partial_{xx} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f^2\cdot\partial_{yy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f\cdot f_x^2$            |            | 🔴        | 🟢      | 🟢        | 🟢          |
| $f\cdot f_y^2$            |            | 🔴        | 🟢      | 🟢        | 🟢          |

All libraries are evaluated in the common polynomial-jet basis (max order 4,
max degree 3).  fp/fn/E are computed against the same ground truth for every
library.  In the standard convention $u_t = \Delta(u^m)$ this is $m = 3$; the $1/3$ prefactor is inessential and absorbed by rescaling time.  The coefficients below belong to the $\nabla\cdot(f^2\nabla f)$ normalisation actually simulated.  PDE-FIND cannot represent $f\cdot f_x^2$ or $f\cdot f_y^2$ (products of derivatives), so its minimum fn is 2.  WSINDy represents PME exactly via the chain rule: $f_t = \frac13\partial_x^2(f^3) + \frac13\partial_y^2(f^3)$, which expands to all 4 jet terms.

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


| Noise | Method        | Backend    | Profile       | Eq best | Multigraphs fp/fn | Multigraphs E | Mult. best             | Non-eq best    | Non-eq fp/fn | Non-eq E |
|:------|:-------------:|:----------:|:-------------:|:-------:|:------------------|:--------------|:-----------------------|:---------------|:-------------|:---------|
| 2%    | weak (no IBP) | spectral   | Radial decay  | 🟥       | 0/0 ✓             | 0.0031        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0024   |
| 2%    | weak (no IBP) | spectral   | Radial growth | 🟩       | 0/0 ✓             | 0.0002        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0014   |
| 2%    | weak (no IBP) | spectral   | Axial decay   | 🟩       | 0/0 ✓             | 0.0001        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0018   |
| 2%    | weak (no IBP) | spectral   | Axial growth  | 🟩       | 0/0 ✓             | 0.0001        | mod_multigraph/STRidge | WSINDy/STRidge | 0/0 ✓        | 0.0015   |
| 2%    | weak (no IBP) | polynomial | Radial decay  | 🟩       | 0/0 ✓             | 0.0014        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0024   |
| 2%    | weak (no IBP) | polynomial | Radial growth | 🟩       | 0/0 ✓             | 0.0004        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0014   |
| 2%    | weak (no IBP) | polynomial | Axial decay   | 🟩       | 0/0 ✓             | 0.0007        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0018   |
| 2%    | weak (no IBP) | polynomial | Axial growth  | 🟩       | 0/0 ✓             | 0.0004        | mod_multigraph/STRidge | WSINDy/STRidge | 0/0 ✓        | 0.0015   |
| 2%    | weak (IBP)    | spectral   | Radial decay  | 🟩       | 0/0 ✓             | 0.0011        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0024   |
| 2%    | weak (IBP)    | spectral   | Radial growth | 🟩       | 0/0 ✓             | 0.0003        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0014   |
| 2%    | weak (IBP)    | spectral   | Axial decay   | 🟩       | 0/0 ✓             | 0.0003        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0018   |
| 2%    | weak (IBP)    | spectral   | Axial growth  | 🟩       | 0/0 ✓             | 0.0002        | mod_multigraph/STRidge | WSINDy/STRidge | 0/0 ✓        | 0.0015   |
| 2%    | weak (IBP)    | polynomial | Radial decay  | 🟩       | 0/0 ✓             | 0.0011        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0024   |
| 2%    | weak (IBP)    | polynomial | Radial growth | 🟩       | 0/0 ✓             | 0.0003        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0014   |
| 2%    | weak (IBP)    | polynomial | Axial decay   | 🟩       | 0/0 ✓             | 0.0003        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0018   |
| 2%    | weak (IBP)    | polynomial | Axial growth  | 🟩       | 0/0 ✓             | 0.0002        | mod_multigraph/STRidge | WSINDy/STRidge | 0/0 ✓        | 0.0015   |
| 5%    | weak (no IBP) | spectral   | Radial decay  | 🟥       | 0/0 ✓             | 0.0474        | mod_multigraph/RFE     | WSINDy/STRidge | 0/0 ✓        | 0.0108   |
| 5%    | weak (no IBP) | spectral   | Radial growth | 🟩       | 0/0 ✓             | 0.0003        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0025   |
| 5%    | weak (no IBP) | spectral   | Axial decay   | 🟥       | 0/0 ✓             | 0.0146        | mod_multigraph/RFE     | WSINDy/STRidge | 0/0 ✓        | 0.0090   |
| 5%    | weak (no IBP) | spectral   | Axial growth  | 🟩       | 0/0 ✓             | 0.0001        | mod_multigraph/RFE     | WSINDy/RFE     | 0/0 ✓        | 0.0028   |
| 5%    | weak (no IBP) | polynomial | Radial decay  | 🟥       | 0/0 ✓             | 0.0208        | multigraph/STRidge     | WSINDy/STRidge | 0/0 ✓        | 0.0107   |
| 5%    | weak (no IBP) | polynomial | Radial growth | 🟩       | 0/0 ✓             | 0.0006        | multigraph/STRidge     | WSINDy/RFE     | 0/0 ✓        | 0.0025   |
| 5%    | weak (no IBP) | polynomial | Axial decay   | 🟩       | 0/0 ✓             | 0.0075        | mod_multigraph/STRidge | WSINDy/STRidge | 0/0 ✓        | 0.0090   |
| 5%    | weak (no IBP) | polynomial | Axial growth  | 🟩       | 0/0 ✓             | 0.0008        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0028   |
| 5%    | weak (IBP)    | spectral   | Radial decay  | 🟩       | 0/0 ✓             | 0.0013        | mod_multigraph/STRidge | WSINDy/STRidge | 0/0 ✓        | 0.0108   |
| 5%    | weak (IBP)    | spectral   | Radial growth | 🟩       | 0/0 ✓             | 0.0004        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0025   |
| 5%    | weak (IBP)    | spectral   | Axial decay   | 🟩       | 0/0 ✓             | 0.0021        | mod_multigraph/STRidge | WSINDy/STRidge | 0/0 ✓        | 0.0090   |
| 5%    | weak (IBP)    | spectral   | Axial growth  | 🟩       | 0/0 ✓             | 0.0006        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0028   |
| 5%    | weak (IBP)    | polynomial | Radial decay  | 🟩       | 0/0 ✓             | 0.0013        | mod_multigraph/STRidge | WSINDy/STRidge | 0/0 ✓        | 0.0107   |
| 5%    | weak (IBP)    | polynomial | Radial growth | 🟩       | 0/0 ✓             | 0.0004        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0025   |
| 5%    | weak (IBP)    | polynomial | Axial decay   | 🟩       | 0/0 ✓             | 0.0026        | mod_multigraph/STRidge | WSINDy/STRidge | 0/0 ✓        | 0.0090   |
| 5%    | weak (IBP)    | polynomial | Axial growth  | 🟩       | 0/0 ✓             | 0.0006        | mod_multigraph/STRidge | WSINDy/RFE     | 0/0 ✓        | 0.0028   |

## Ablation

### 2% — Axial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best  |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:----------|
| 1    | Modified Multigraphs + STRidge | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0001  | 0.1      | —         |
| 2    | Modified Multigraphs + RFE     | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0002  | 0.0      | +115.4%   |
| 3    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0011  | 0.1      | +1174.2%  |
| 4    | WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0018  | 0.0      | +1944.8%  |
| 5    | WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0035  | 0.2      | +3873.0%  |
| 6    | Multigraphs + RFE              | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0040  | 0.0      | +4438.3%  |
| 7    | PolyJet + STRidge              | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0104  | 0.4      | +11772.6% |
| 8    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/1   | 0.2978  | 0.0      | ✗         |

### 2% — Radial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0011  | 0.1      | —        |
| 2    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0014  | 0.1      | +26.9%   |
| 3    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0022  | 0.0      | +107.7%  |
| 4    | WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0024  | 0.0      | +118.5%  |
| 5    | WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0040  | 0.2      | +269.4%  |
| 6    | Multigraphs + RFE              | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0095  | 0.0      | +784.6%  |
| 7    | PolyJet + STRidge              | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0157  | 0.4      | +1355.2% |
| 8    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/1   | 0.2980  | 0.0      | ✗        |

### 2% — Axial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0001  | 0.1      | —        |
| 2    | Multigraphs + STRidge          | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0009  | 0.1      | +753.3%  |
| 3    | Modified Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0011  | 0.0      | +935.8%  |
| 4    | WSINDy + STRidge               | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0015  | 0.2      | +1354.1% |
| 5    | WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0015  | 0.0      | +1397.3% |
| 6    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0024  | 0.0      | +2279.6% |
| 7    | PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0052  | 0.5      | +5046.9% |
| 8    | PolyJet + RFE                  | polynomial, weak (no IBP), m=8  | 0/1   | 0.2976  | 0.0      | ✗        |

### 2% — Radial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0002  | 0.1      | —        |
| 2    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0009  | 0.1      | +370.3%  |
| 3    | Modified Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0011  | 0.0      | +488.8%  |
| 4    | WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0014  | 0.0      | +618.3%  |
| 5    | WSINDy + STRidge               | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0017  | 0.2      | +763.7%  |
| 6    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0022  | 0.0      | +1027.4% |
| 7    | PolyJet + STRidge              | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0056  | 0.7      | +2801.8% |
| 8    | PolyJet + RFE                  | polynomial, weak (no IBP), m=8  | 0/1   | 0.2976  | 0.1      | ✗        |

### 5% — Axial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=8       | 0/0 ✓ | 0.0021  | 0.1      | —        |
| 2    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=8       | 0/0 ✓ | 0.0024  | 0.0      | +10.8%   |
| 3    | WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0090  | 0.2      | +321.5%  |
| 4    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0112  | 0.1      | +426.0%  |
| 5    | PolyJet + STRidge              | spectral, weak (no IBP), m=12   | 1/0   | 0.0298  | 0.6      | ✗        |
| 6    | WSINDy + RFE                   | polynomial, weak (no IBP), m=8  | 0/2   | 0.5773  | 0.0      | ✗        |
| 7    | Multigraphs + RFE              | spectral, weak (no IBP), m=8    | 6/2   | 0.6270  | 0.0      | ✗        |
| 8    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 2/3   | 0.6545  | 0.1      | ✗        |

### 5% — Radial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=8     | 0/0 ✓ | 0.0013  | 0.2      | —        |
| 2    | WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0107  | 0.2      | +700.9%  |
| 3    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0208  | 0.1      | +1450.6% |
| 4    | Modified Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0474  | 0.0      | +3437.1% |
| 5    | PolyJet + STRidge              | spectral, weak (no IBP), m=12   | 1/0   | 0.0455  | 0.4      | ✗        |
| 6    | Multigraphs + RFE              | spectral, weak (no IBP), m=12   | 6/0   | 0.0359  | 0.0      | ✗        |
| 7    | WSINDy + RFE                   | polynomial, weak (no IBP), m=8  | 0/2   | 0.5773  | 0.0      | ✗        |
| 8    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 4/3   | 0.6400  | 0.0      | ✗        |

### 5% — Axial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0001  | 0.0      | —        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0004  | 0.1      | +522.1%  |
| 3    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0008  | 0.1      | +1253.3% |
| 4    | WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0028  | 0.0      | +4535.2% |
| 5    | WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0037  | 0.2      | +5931.1% |
| 6    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0043  | 0.0      | +6939.8% |
| 7    | PolyJet + STRidge              | spectral, weak (no IBP), m=8    | 1/0   | 0.0156  | 0.6      | ✗        |
| 8    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/1   | 0.2981  | 0.1      | ✗        |

### 5% — Radial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0003  | 0.1      | —        |
| 2    | Multigraphs + STRidge          | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0006  | 0.1      | +66.9%   |
| 3    | Modified Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0011  | 0.0      | +220.7%  |
| 4    | WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0025  | 0.0      | +630.8%  |
| 5    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0036  | 0.0      | +959.0%  |
| 6    | WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0040  | 0.2      | +1099.7% |
| 7    | PolyJet + STRidge              | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0128  | 0.6      | +3696.4% |
| 8    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/1   | 0.2983  | 0.1      | ✗        |


## Detailed results

### weak (no IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 2        | 0.0018  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 2        | 0.0024  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 2        | 0.0015  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 2        | 0.0014  | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 2        | 0.0035  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 2        | 0.0040  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 2        | 0.0015  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.01      | Radial growth | 0/0   | 2        | 0.0017  | 0.3      | 🟢      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 1        | 0.5773  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 1        | 0.5773  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 2        | 0.0028  | 0.0      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 2        | 0.0025  | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 2        | 0.0090  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 2        | 0.0107  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 2        | 0.0037  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 2        | 0.0040  | 0.2      | 🟢      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 1/3   | 2        | 0.3701  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/4   | 1        | 0.6272  | 0.0      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial growth  | 0/1   | 3        | 0.2976  | 0.0      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial growth | 0/1   | 3        | 0.2976  | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 14/0  | 18       | 0.0708  | 0.5      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 0.01      | Radial decay  | 21/0  | 25       | 0.0825  | 0.8      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0052  | 0.5      | 🟢      |
| PolyJet     | 2%    | STRidge   | 8  | 0.1       | Radial growth | 0/0   | 4        | 0.0056  | 0.7      | 🟢      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 2/4   | 2        | 0.2599  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 2/4   | 2        | 0.2644  | 0.0      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/1   | 3        | 0.3040  | 0.0      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/3   | 1        | 0.6696  | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 1.0       | Axial decay   | 2/0   | 6        | 0.0424  | 0.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 1/0   | 5        | 0.0555  | 0.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.01      | Axial growth  | 20/0  | 24       | 0.0463  | 0.9      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Radial growth | 14/0  | 18       | 0.0308  | 0.7      | 🔴      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 1        | 0.0033  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial decay  | 3/2   | 2        | 0.2695  | 0.0      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 1        | 0.0016  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 1        | 0.0014  | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 1        | 0.0007  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 1        | 0.0014  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 1        | 0.0004  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 1        | 0.0004  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 4/4   | 2        | 0.5509  | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 2/4   | 2        | 0.4418  | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 1        | 0.0037  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 1        | 0.0030  | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 0/0   | 1        | 0.0075  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 0.1       | Radial decay  | 0/0   | 2        | 0.0208  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 1        | 0.0008  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 2        | 0.0006  | 0.1      | 🟢      |

### weak (no IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 2        | 0.0018  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 2        | 0.0024  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 2        | 0.0015  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 2        | 0.0014  | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 2        | 0.0035  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 2        | 0.0040  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 2        | 0.0015  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 2        | 0.0017  | 0.2      | 🟢      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 1        | 0.5774  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 1        | 0.5774  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 2        | 0.0028  | 0.0      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 2        | 0.0025  | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 2        | 0.0090  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 2        | 0.0108  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 2        | 0.0037  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 2        | 0.0041  | 0.2      | 🟢      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/1   | 3        | 0.2978  | 0.0      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 0/1   | 3        | 0.2980  | 0.0      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial growth  | 0/1   | 3        | 0.2977  | 0.0      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial growth | 0/1   | 3        | 0.2977  | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 4        | 0.0104  | 0.4      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 4        | 0.0157  | 0.4      | 🟢      |
| PolyJet     | 2%    | STRidge   | 8  | 0.1       | Axial growth  | 0/0   | 4        | 0.0057  | 0.6      | 🟢      |
| PolyJet     | 2%    | STRidge   | 8  | 0.1       | Radial growth | 0/0   | 4        | 0.0056  | 0.7      | 🟢      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 2/3   | 3        | 0.6545  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 4/3   | 5        | 0.6400  | 0.0      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/1   | 3        | 0.2981  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/1   | 3        | 0.2983  | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 0.1       | Axial decay   | 1/0   | 5        | 0.0298  | 0.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 1/0   | 5        | 0.0455  | 0.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Axial growth  | 1/0   | 5        | 0.0156  | 0.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Radial growth | 0/0   | 4        | 0.0128  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 1        | 0.0002  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial decay  | 0/0   | 1        | 0.0035  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 1        | 0.0011  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 1        | 0.0011  | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 8  | 0.1       | Axial decay   | 0/0   | 1        | 0.0001  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 8  | 10.0      | Radial decay  | 0/0   | 1        | 0.0031  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 1        | 0.0001  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 1        | 0.0002  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 1        | 0.0146  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 1        | 0.0474  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 1        | 0.0001  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 1        | 0.0011  | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Axial decay   | 0/0   | 1        | 0.0150  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 1        | 0.0555  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 0.1       | Axial growth  | 0/0   | 1        | 0.0004  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 10.0      | Radial growth | 0/0   | 1        | 0.0003  | 0.1      | 🟢      |

### weak (IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 2        | 0.0018  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 2        | 0.0024  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 2        | 0.0015  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 2        | 0.0014  | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 2        | 0.0035  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 2        | 0.0040  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 2        | 0.0015  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 2        | 0.0017  | 0.2      | 🟢      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 1        | 0.5773  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 1        | 0.5773  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 2        | 0.0028  | 0.0      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 2        | 0.0025  | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 2        | 0.0090  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 2        | 0.0107  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 2        | 0.0037  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 2        | 0.0040  | 0.2      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/2   | 2        | 0.3545  | 0.0      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 1/4   | 1        | 0.6272  | 0.0      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial growth  | 0/1   | 3        | 0.2976  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial growth | 0/1   | 3        | 0.2976  | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 12/0  | 16       | 0.0325  | 0.5      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 0.01      | Radial decay  | 21/0  | 25       | 0.0825  | 0.9      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0052  | 0.5      | 🟢      |
| PolyJet     | 2%    | STRidge   | 8  | 1.0       | Radial growth | 0/0   | 4        | 0.0056  | 0.5      | 🟢      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 2/4   | 2        | 0.2599  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 2/4   | 2        | 0.2644  | 0.0      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/1   | 3        | 0.3040  | 0.0      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/3   | 1        | 0.6696  | 0.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 10/0  | 14       | 0.1363  | 0.7      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 1/0   | 5        | 0.0555  | 0.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Axial growth  | 4/0   | 8        | 0.0215  | 0.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Radial growth | 14/0  | 18       | 0.0308  | 0.7      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 1        | 0.0017  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 8  | 0.0       | Radial decay  | 3/2   | 2        | 0.2695  | 0.0      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 1        | 0.0014  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 1        | 0.0012  | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 1        | 0.0003  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 1        | 0.0011  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 1        | 0.0002  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 1        | 0.0003  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial decay   | 2/4   | 2        | 0.2640  | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 2/4   | 2        | 0.4418  | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 1        | 0.0022  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 1        | 0.0020  | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 1        | 0.0026  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 0.1       | Radial decay  | 0/0   | 1        | 0.0013  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 1        | 0.0006  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 1        | 0.0004  | 0.1      | 🟢      |

### weak (IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 2        | 0.0018  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 2        | 0.0024  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 2        | 0.0015  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 2        | 0.0014  | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 2        | 0.0035  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 2        | 0.0040  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 2        | 0.0015  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 2        | 0.0017  | 0.2      | 🟢      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 1        | 0.5774  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 1        | 0.5774  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 2        | 0.0028  | 0.0      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 2        | 0.0025  | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 2        | 0.0090  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 2        | 0.0108  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 2        | 0.0037  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 2        | 0.0041  | 0.2      | 🟢      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial decay   | 0/1   | 3        | 0.2978  | 0.0      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial decay  | 0/1   | 3        | 0.2980  | 0.1      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Axial growth  | 0/1   | 3        | 0.2977  | 0.0      | 🔴      |
| PolyJet     | 2%    | RFE       | 8  | 0.0       | Radial growth | 0/1   | 3        | 0.2977  | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 4        | 0.0104  | 0.4      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 4        | 0.0157  | 0.4      | 🟢      |
| PolyJet     | 2%    | STRidge   | 8  | 0.1       | Axial growth  | 0/0   | 4        | 0.0057  | 0.6      | 🟢      |
| PolyJet     | 2%    | STRidge   | 8  | 0.1       | Radial growth | 0/0   | 4        | 0.0056  | 0.7      | 🟢      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 2/3   | 3        | 0.6545  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 4/3   | 5        | 0.6400  | 0.0      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/1   | 3        | 0.2981  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/1   | 3        | 0.2983  | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 0.1       | Axial decay   | 1/0   | 5        | 0.0298  | 0.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 1/0   | 5        | 0.0455  | 0.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Axial growth  | 1/0   | 5        | 0.0156  | 0.6      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Radial growth | 0/0   | 4        | 0.0128  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 1        | 0.0017  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 1        | 0.0022  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 1        | 0.0014  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 1        | 0.0012  | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 1        | 0.0003  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 1        | 0.0011  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 1        | 0.0002  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 1        | 0.0003  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 1        | 0.0024  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 6/0   | 4        | 0.0359  | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 1        | 0.0022  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 1        | 0.0020  | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 1.0       | Axial decay   | 0/0   | 1        | 0.0021  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 0/0   | 1        | 0.0013  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 1        | 0.0006  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 1        | 0.0004  | 0.1      | 🟢      |


---
*Source: `sweep_aniso_20260803_194535_full_local.csv`*

*Run: 775.9 min (Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz), OMP_NUM_THREADS=1*
