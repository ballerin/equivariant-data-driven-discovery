# Swift-Hohenberg equation — Anisotropic noise

The equation under consideration is

$$f_t = (r-1)f - 2\Delta f - \Delta^2 f - f^3, \qquad r = 0.3$$

which can be written in terms of $O(2)$-equivariant terms

| Operator     | True coefficient |
|:-------------|:-----------------|
| $f$          | $r-1 = -0.7$     |
| $\Delta f$   | $-2.0$           |
| $\Delta^2 f$ | $-1.0$           |
| $f^3$        | $-1.0$           |

or standard terms (PolyJet basis)

| Operator            | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|:--------------------|:-----------|:---------|:-------|:---------|:-----------|
| $f$                 |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f^3$               |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xx} f$   |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yy} f$   |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxxx} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxyy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yyyy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |

All libraries are evaluated in the common polynomial-jet basis (max order 4,
max degree 3).  fp/fn/E are computed against the same ground truth for every
library.  Swift-Hohenberg contains a 4th-order biharmonic term $\Delta^2 f = f_{xxxx} + 2f_{xxyy} + f_{yyyy}$ which is the main source of difficulty.  All four libraries can represent every term individually, but the equivariant basis encodes the isotropic combination as a single operator, giving it a sparsity advantage

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


| Noise | Method        | Backend    | Profile       | Eq best | Multigraphs fp/fn | Multigraphs E | Mult. best             | Non-eq best      | Non-eq fp/fn | Non-eq E |
|:------|:-------------:|:----------:|:-------------:|:-------:|:------------------|:--------------|:-----------------------|:-----------------|:-------------|:---------|
| 2%    | weak (no IBP) | spectral   | Radial decay  | 🟥       | 0/0 ✓             | 0.0384        | mod_multigraph/STRidge | PolyJet/STRidge  | 0/0 ✓        | 0.0381   |
| 2%    | weak (no IBP) | spectral   | Radial growth | 🟥       | 0/0 ✓             | 0.0285        | mod_multigraph/STRidge | WSINDy/STRidge   | 0/0 ✓        | 0.0280   |
| 2%    | weak (no IBP) | spectral   | Axial decay   | 🟥       | 0/0 ✓             | 0.0329        | multigraph/STRidge     | PolyJet/STRidge  | 0/0 ✓        | 0.0325   |
| 2%    | weak (no IBP) | spectral   | Axial growth  | 🟥       | 0/0 ✓             | 0.0287        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0282   |
| 2%    | weak (no IBP) | polynomial | Radial decay  | 🟥       | 0/0 ✓             | 0.0366        | mod_multigraph/STRidge | PDE-FIND/STRidge | 0/0 ✓        | 0.0363   |
| 2%    | weak (no IBP) | polynomial | Radial growth | 🟥       | 0/0 ✓             | 0.0276        | mod_multigraph/STRidge | PDE-FIND/STRidge | 0/0 ✓        | 0.0272   |
| 2%    | weak (no IBP) | polynomial | Axial decay   | 🟥       | 0/0 ✓             | 0.0316        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0311   |
| 2%    | weak (no IBP) | polynomial | Axial growth  | 🟥       | 0/0 ✓             | 0.0278        | mod_multigraph/STRidge | PDE-FIND/STRidge | 0/0 ✓        | 0.0274   |
| 2%    | weak (IBP)    | spectral   | Radial decay  | 🟥       | 0/0 ✓             | 0.0377        | multigraph/STRidge     | PolyJet/STRidge  | 0/0 ✓        | 0.0316   |
| 2%    | weak (IBP)    | spectral   | Radial growth | 🟥       | 0/0 ✓             | 0.0268        | multigraph/STRidge     | WSINDy/STRidge   | 0/0 ✓        | 0.0230   |
| 2%    | weak (IBP)    | spectral   | Axial decay   | 🟥       | 0/0 ✓             | 0.0320        | multigraph/STRidge     | PolyJet/STRidge  | 0/0 ✓        | 0.0266   |
| 2%    | weak (IBP)    | spectral   | Axial growth  | 🟥       | 0/0 ✓             | 0.0269        | mod_multigraph/RFE     | PolyJet/STRidge  | 0/0 ✓        | 0.0231   |
| 2%    | weak (IBP)    | polynomial | Radial decay  | 🟥       | 0/0 ✓             | 0.0377        | mod_multigraph/STRidge | PDE-FIND/STRidge | 0/0 ✓        | 0.0316   |
| 2%    | weak (IBP)    | polynomial | Radial growth | 🟥       | 0/0 ✓             | 0.0268        | mod_multigraph/STRidge | PDE-FIND/STRidge | 0/0 ✓        | 0.0230   |
| 2%    | weak (IBP)    | polynomial | Axial decay   | 🟥       | 0/0 ✓             | 0.0320        | mod_multigraph/STRidge | WSINDy/STRidge   | 0/0 ✓        | 0.0266   |
| 2%    | weak (IBP)    | polynomial | Axial growth  | 🟥       | 0/0 ✓             | 0.0269        | multigraph/RFE         | WSINDy/STRidge   | 0/0 ✓        | 0.0231   |
| 5%    | weak (no IBP) | spectral   | Radial decay  | 🟢       | 0/0 ✓             | 0.1618        | mod_multigraph/STRidge | —                | 2/0 ✗        | 0.2238   |
| 5%    | weak (no IBP) | spectral   | Radial growth | 🟩       | 0/0 ✓             | 0.0742        | mod_multigraph/STRidge | PDE-FIND/STRidge | 0/0 ✓        | 0.0756   |
| 5%    | weak (no IBP) | spectral   | Axial decay   | 🟩       | 0/0 ✓             | 0.1218        | mod_multigraph/STRidge | WSINDy/STRidge   | 0/0 ✓        | 0.1240   |
| 5%    | weak (no IBP) | spectral   | Axial growth  | 🟩       | 0/0 ✓             | 0.0760        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0776   |
| 5%    | weak (no IBP) | polynomial | Radial decay  | 🟢       | 0/0 ✓             | 0.1561        | mod_multigraph/STRidge | —                | 5/0 ✗        | 0.3145   |
| 5%    | weak (no IBP) | polynomial | Radial growth | 🟩       | 0/0 ✓             | 0.0713        | mod_multigraph/STRidge | WSINDy/STRidge   | 0/0 ✓        | 0.0726   |
| 5%    | weak (no IBP) | polynomial | Axial decay   | 🟩       | 0/0 ✓             | 0.1174        | mod_multigraph/STRidge | WSINDy/STRidge   | 0/0 ✓        | 0.1196   |
| 5%    | weak (no IBP) | polynomial | Axial growth  | 🟩       | 0/0 ✓             | 0.0730        | mod_multigraph/STRidge | WSINDy/STRidge   | 0/0 ✓        | 0.0745   |
| 5%    | weak (IBP)    | spectral   | Radial decay  | 🟢       | 0/0 ✓             | 0.1619        | mod_multigraph/STRidge | —                | 2/0 ✗        | 0.2190   |
| 5%    | weak (IBP)    | spectral   | Radial growth | 🟥       | 0/0 ✓             | 0.0745        | mod_multigraph/STRidge | PDE-FIND/STRidge | 0/0 ✓        | 0.0695   |
| 5%    | weak (IBP)    | spectral   | Axial decay   | 🟥       | 0/0 ✓             | 0.1221        | mod_multigraph/STRidge | WSINDy/STRidge   | 0/0 ✓        | 0.1182   |
| 5%    | weak (IBP)    | spectral   | Axial growth  | 🟥       | 0/0 ✓             | 0.0760        | mod_multigraph/RFE     | PDE-FIND/STRidge | 0/0 ✓        | 0.0715   |
| 5%    | weak (IBP)    | polynomial | Radial decay  | 🟢       | 0/0 ✓             | 0.1619        | multigraph/STRidge     | —                | 5/0 ✗        | 0.3197   |
| 5%    | weak (IBP)    | polynomial | Radial growth | 🟥       | 0/0 ✓             | 0.0745        | multigraph/STRidge     | WSINDy/STRidge   | 0/0 ✓        | 0.0695   |
| 5%    | weak (IBP)    | polynomial | Axial decay   | 🟥       | 0/0 ✓             | 0.1221        | multigraph/STRidge     | WSINDy/STRidge   | 0/0 ✓        | 0.1182   |
| 5%    | weak (IBP)    | polynomial | Axial growth  | 🟥       | 0/0 ✓             | 0.0760        | multigraph/RFE         | PDE-FIND/STRidge | 0/0 ✓        | 0.0715   |

## Ablation

### 2% — Axial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | PolyJet + STRidge              | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0266  | 0.8      | —        |
| 2    | WSINDy + STRidge               | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0266  | 0.4      | +0.0%    |
| 3    | PDE-FIND + STRidge             | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0266  | 0.4      | +0.0%    |
| 4    | PolyJet + RFE                  | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0282  | 0.1      | +6.0%    |
| 5    | WSINDy + RFE                   | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0282  | 0.1      | +6.0%    |
| 6    | PDE-FIND + RFE                 | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0282  | 0.1      | +6.0%    |
| 7    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0316  | 0.2      | +18.6%   |
| 8    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0320  | 0.2      | +20.3%   |
| 9    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0333  | 0.1      | +25.2%   |
| 10   | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0333  | 0.1      | +25.2%   |

### 2% — Radial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | PolyJet + STRidge              | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0316  | 0.8      | —        |
| 2    | PDE-FIND + STRidge             | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0316  | 0.4      | +0.0%    |
| 3    | WSINDy + STRidge               | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0316  | 0.4      | +0.0%    |
| 4    | PolyJet + RFE                  | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0341  | 0.1      | +8.0%    |
| 5    | PDE-FIND + RFE                 | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0341  | 0.1      | +8.0%    |
| 6    | WSINDy + RFE                   | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0341  | 0.1      | +8.0%    |
| 7    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0366  | 0.2      | +15.9%   |
| 8    | Multigraphs + STRidge          | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0377  | 0.2      | +19.3%   |
| 9    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0389  | 0.1      | +23.0%   |
| 10   | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0389  | 0.1      | +23.0%   |

### 2% — Axial growth

| Rank | Library + Algorithm            | Configuration                | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-----------------------------|:------|:--------|:---------|:---------|
| 1    | PolyJet + STRidge              | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0231  | 0.8      | —        |
| 2    | PDE-FIND + STRidge             | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0231  | 0.4      | +0.0%    |
| 3    | WSINDy + STRidge               | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0231  | 0.4      | +0.0%    |
| 4    | PolyJet + RFE                  | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0233  | 0.1      | +0.7%    |
| 5    | PDE-FIND + RFE                 | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0233  | 0.1      | +0.7%    |
| 6    | WSINDy + RFE                   | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0233  | 0.1      | +0.7%    |
| 7    | Multigraphs + RFE              | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0269  | 0.1      | +16.4%   |
| 8    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0269  | 0.1      | +16.4%   |
| 9    | Multigraphs + STRidge          | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0270  | 0.2      | +16.7%   |
| 10   | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0270  | 0.2      | +16.7%   |

### 2% — Radial growth

| Rank | Library + Algorithm            | Configuration                  | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-------------------------------|:------|:--------|:---------|:---------|
| 1    | PDE-FIND + STRidge             | polynomial, weak (IBP), m=12   | 0/0 ✓ | 0.0230  | 0.4      | —        |
| 2    | WSINDy + STRidge               | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0230  | 0.4      | +0.0%    |
| 3    | PolyJet + RFE                  | polynomial, weak (IBP), m=12   | 0/0 ✓ | 0.0232  | 0.1      | +0.9%    |
| 4    | PDE-FIND + RFE                 | polynomial, weak (IBP), m=12   | 0/0 ✓ | 0.0232  | 0.1      | +0.9%    |
| 5    | WSINDy + RFE                   | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0232  | 0.1      | +0.9%    |
| 6    | Multigraphs + STRidge          | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0268  | 0.2      | +16.3%   |
| 7    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12   | 0/0 ✓ | 0.0268  | 0.2      | +16.3%   |
| 8    | Multigraphs + RFE              | polynomial, weak (IBP), m=12   | 0/0 ✓ | 0.0268  | 0.1      | +16.5%   |
| 9    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0268  | 0.1      | +16.5%   |
| 10   | PolyJet + STRidge              | polynomial, weak (no IBP), m=8 | 0/0 ✓ | 0.3053  | 1.2      | +1225.9% |

### 5% — Axial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1174  | 0.2      | —        |
| 2    | WSINDy + STRidge               | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.1182  | 0.5      | +0.7%    |
| 3    | Multigraphs + STRidge          | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.1218  | 0.3      | +3.8%    |
| 4    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1238  | 0.1      | +5.4%    |
| 5    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1238  | 0.1      | +5.4%    |
| 6    | PDE-FIND + STRidge             | spectral, weak (IBP), m=12      | 2/0   | 0.1729  | 0.8      | ✗        |
| 7    | PolyJet + STRidge              | spectral, weak (IBP), m=12      | 9/0   | 0.2663  | 1.1      | ✗        |
| 8    | WSINDy + RFE                   | spectral, weak (IBP), m=12      | 0/4   | 0.8228  | 0.1      | ✗        |
| 9    | PDE-FIND + RFE                 | spectral, weak (IBP), m=12      | 0/4   | 0.8228  | 0.1      | ✗        |
| 10   | PolyJet + RFE                  | spectral, weak (IBP), m=12      | 0/4   | 0.8228  | 0.1      | ✗        |

### 5% — Radial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1561  | 0.2      | —        |
| 2    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1561  | 0.2      | +0.0%    |
| 3    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1638  | 0.1      | +4.9%    |
| 4    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1638  | 0.1      | +4.9%    |
| 5    | PDE-FIND + STRidge             | spectral, weak (IBP), m=12      | 2/0   | 0.2190  | 0.9      | ✗        |
| 6    | PolyJet + STRidge              | spectral, weak (IBP), m=8       | 10/0  | 1.5215  | 1.5      | ✗        |
| 7    | WSINDy + STRidge               | spectral, weak (IBP), m=12      | 18/0  | 0.9991  | 0.5      | ✗        |
| 8    | WSINDy + RFE                   | polynomial, weak (IBP), m=12    | 0/4   | 0.8210  | 0.1      | ✗        |
| 9    | PolyJet + RFE                  | spectral, weak (IBP), m=12      | 0/4   | 0.8210  | 0.1      | ✗        |
| 10   | PDE-FIND + RFE                 | spectral, weak (IBP), m=12      | 0/4   | 0.8210  | 0.1      | ✗        |

### 5% — Axial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | PDE-FIND + STRidge             | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0715  | 0.9      | —        |
| 2    | WSINDy + STRidge               | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0715  | 0.5      | +0.0%    |
| 3    | WSINDy + RFE                   | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0723  | 0.1      | +1.2%    |
| 4    | PolyJet + RFE                  | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0723  | 0.1      | +1.2%    |
| 5    | PDE-FIND + RFE                 | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0723  | 0.1      | +1.2%    |
| 6    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0730  | 0.2      | +2.1%    |
| 7    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0732  | 0.1      | +2.4%    |
| 8    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0732  | 0.1      | +2.4%    |
| 9    | Multigraphs + STRidge          | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0760  | 0.2      | +6.3%    |
| 10   | PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 5/0   | 1.5002  | 1.0      | ✗        |

### 5% — Radial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | WSINDy + STRidge               | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0695  | 0.4      | —        |
| 2    | PDE-FIND + STRidge             | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0695  | 0.7      | +0.0%    |
| 3    | PDE-FIND + RFE                 | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0708  | 0.2      | +1.8%    |
| 4    | PolyJet + RFE                  | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0708  | 0.1      | +1.8%    |
| 5    | WSINDy + RFE                   | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0708  | 0.1      | +1.8%    |
| 6    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0713  | 0.2      | +2.5%    |
| 7    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0713  | 0.2      | +2.5%    |
| 8    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0720  | 0.1      | +3.5%    |
| 9    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0720  | 0.1      | +3.5%    |
| 10   | PolyJet + STRidge              | spectral, weak (IBP), m=12      | 2/0   | 0.0864  | 0.9      | ✗        |


## Detailed results

### weak (no IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0334  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0393  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0282  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0280  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 7        | 0.0311  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 7        | 0.0363  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0274  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0272  | 0.4      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8237  | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8219  | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0759  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0744  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 5/0   | 12       | 0.2426  | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 5/0   | 12       | 0.3145  | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 8  | 0.1       | Axial growth  | 5/0   | 12       | 1.0499  | 0.8      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Radial growth | 15/0  | 22       | 0.1665  | 0.6      | 🔴      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0334  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0393  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0282  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0280  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 7        | 0.0311  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 7        | 0.0363  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 7        | 0.0274  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 7        | 0.0272  | 0.4      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8237  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8219  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0759  | 0.1      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0744  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 7        | 0.1196  | 0.4      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 0.1       | Radial decay  | 18/0  | 13       | 1.0351  | 0.4      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0745  | 0.4      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0726  | 0.4      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0334  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0393  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0282  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0280  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 11/0  | 18       | 2.3894  | 0.9      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 6/0   | 13       | 2.4411  | 1.0      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 3/0   | 10       | 1.2328  | 1.0      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 1.0       | Radial growth | 0/0   | 7        | 0.3053  | 1.2      | 🟢      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8237  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8219  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0759  | 0.1      | 🟢      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0744  | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 21/0  | 28       | 3.8957  | 1.3      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 16/0  | 23       | 2.8793  | 1.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 5/0   | 12       | 1.5002  | 1.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial growth | 7/0   | 14       | 1.6286  | 1.0      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0333  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0389  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0284  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0283  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.0316  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.0366  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Axial growth  | 0/0   | 4        | 0.0278  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0276  | 0.2      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.1238  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.1638  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0732  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0720  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.1174  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.1561  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0730  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0713  | 0.2      | 🟢      |

### weak (no IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0348  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0411  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0290  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0288  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 2/0   | 9        | 0.0448  | 0.5      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 7        | 0.0381  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0282  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Radial growth | 2/0   | 9        | 0.0329  | 0.5      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8238  | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8220  | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0790  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0774  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 2/0   | 9        | 0.1779  | 0.9      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 2/0   | 9        | 0.2238  | 0.9      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0776  | 0.9      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 7        | 0.0756  | 0.9      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0348  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0411  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0290  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0288  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 7        | 0.0325  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 7        | 0.0381  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 7        | 0.0282  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0280  | 0.3      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8238  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8220  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0790  | 0.1      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0774  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 7        | 0.1240  | 0.4      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 0.1       | Radial decay  | 18/0  | 13       | 1.0033  | 0.5      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0776  | 0.5      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0756  | 0.5      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0348  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0411  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0290  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0288  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 7        | 0.0325  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 7        | 0.0381  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0282  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial growth | 6/0   | 13       | 0.1319  | 0.9      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8238  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8220  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0790  | 0.1      | 🟢      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0774  | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 9/0   | 16       | 0.2707  | 1.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 0.1       | Radial decay  | 17/0  | 24       | 0.6322  | 1.3      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 11/0  | 18       | 0.2144  | 1.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial growth | 2/0   | 9        | 0.0923  | 0.9      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0347  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0407  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0292  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0291  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.0329  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.0384  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0287  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 4        | 0.0285  | 0.2      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.1283  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.1694  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0763  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0750  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.1218  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.1618  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0760  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0742  | 0.2      | 🟢      |

### weak (IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0282  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0341  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0233  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0232  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 7        | 0.0266  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 7        | 0.0316  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0231  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0230  | 0.4      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8228  | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8210  | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0723  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0708  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 0.1       | Axial decay   | 4/0   | 11       | 0.2176  | 0.6      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 5/0   | 12       | 0.3197  | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0715  | 0.4      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Radial growth | 15/0  | 22       | 0.1689  | 0.6      | 🔴      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0282  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0341  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0233  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0232  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 7        | 0.0266  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 7        | 0.0316  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 7        | 0.0231  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 7        | 0.0230  | 0.4      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8228  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8210  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0723  | 0.1      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0708  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 7        | 0.1182  | 0.3      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 0.1       | Radial decay  | 18/0  | 13       | 1.0095  | 0.4      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0715  | 0.4      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0695  | 0.4      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0282  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0341  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0233  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0232  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 4/0   | 11       | 2.1134  | 0.9      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 7        | 0.0316  | 0.9      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 3/0   | 10       | 1.1633  | 0.9      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial growth | 6/0   | 13       | 1.1412  | 1.0      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8228  | 0.2      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8210  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0723  | 0.1      | 🟢      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0708  | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 14/0  | 21       | 2.9865  | 1.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 16/0  | 23       | 2.8866  | 1.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 14/0  | 21       | 2.6024  | 1.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial growth | 33/0  | 40       | 3.3528  | 1.4      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0334  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0396  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0269  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0268  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.0320  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.0377  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0270  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0268  | 0.2      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.1282  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.1693  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0760  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0748  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.1221  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.1619  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0762  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0745  | 0.2      | 🟢      |

### weak (IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0282  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0341  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0233  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0232  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 2/0   | 9        | 0.0391  | 0.5      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 7        | 0.0316  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0231  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Radial growth | 2/0   | 9        | 0.0330  | 0.5      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8228  | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8210  | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0723  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0708  | 0.2      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 2/0   | 9        | 0.1729  | 0.8      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 2/0   | 9        | 0.2190  | 0.9      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0715  | 0.9      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0695  | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0282  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0341  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0233  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0232  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 7        | 0.0266  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 7        | 0.0316  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 7        | 0.0231  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 7        | 0.0230  | 0.4      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8228  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8210  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0723  | 0.1      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0708  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 7        | 0.1182  | 0.5      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 0.1       | Radial decay  | 18/0  | 13       | 0.9991  | 0.5      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0715  | 0.5      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0695  | 0.5      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0282  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0341  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0233  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0232  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 7        | 0.0266  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 7        | 0.0316  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0231  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial growth | 6/0   | 13       | 0.1363  | 0.9      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8228  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8210  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0723  | 0.1      | 🟢      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0708  | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 9/0   | 16       | 0.2663  | 1.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Radial decay  | 10/0  | 17       | 1.5215  | 1.5      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 11/0  | 18       | 0.2118  | 1.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial growth | 2/0   | 9        | 0.0864  | 0.9      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0334  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0396  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0269  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0268  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.0320  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 4        | 0.0377  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0270  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0268  | 0.2      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.1282  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.1693  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0760  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0748  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.1221  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.1619  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0762  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0745  | 0.2      | 🟢      |


---
*Source: `sweep_aniso_default_sh_20260805_215410.csv`*

