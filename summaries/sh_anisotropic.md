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
| 2%    | weak (no IBP) | spectral   | Radial decay  | 🟥       | 0/0 ✓             | 0.0384        | multigraph/STRidge     | PolyJet/STRidge  | 0/0 ✓        | 0.0381   |
| 2%    | weak (no IBP) | spectral   | Radial growth | 🟥       | 0/0 ✓             | 0.0285        | multigraph/STRidge     | WSINDy/STRidge   | 0/0 ✓        | 0.0280   |
| 2%    | weak (no IBP) | spectral   | Axial decay   | 🟥       | 0/0 ✓             | 0.0329        | multigraph/STRidge     | WSINDy/STRidge   | 0/0 ✓        | 0.0325   |
| 2%    | weak (no IBP) | spectral   | Axial growth  | 🟥       | 0/0 ✓             | 0.0287        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0282   |
| 2%    | weak (no IBP) | polynomial | Radial decay  | 🟥       | 0/0 ✓             | 0.0356        | mod_multigraph/STRidge | PDE-FIND/STRidge | 0/0 ✓        | 0.0354   |
| 2%    | weak (no IBP) | polynomial | Radial growth | 🟥       | 0/0 ✓             | 0.0283        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0279   |
| 2%    | weak (no IBP) | polynomial | Axial decay   | 🟥       | 0/0 ✓             | 0.0316        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0312   |
| 2%    | weak (no IBP) | polynomial | Axial growth  | 🟥       | 0/0 ✓             | 0.0285        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0280   |
| 2%    | weak (IBP)    | spectral   | Radial decay  | 🟥       | 0/0 ✓             | 0.0377        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0316   |
| 2%    | weak (IBP)    | spectral   | Radial growth | 🟥       | 0/0 ✓             | 0.0268        | multigraph/STRidge     | WSINDy/STRidge   | 0/0 ✓        | 0.0230   |
| 2%    | weak (IBP)    | spectral   | Axial decay   | 🟥       | 0/0 ✓             | 0.0320        | multigraph/STRidge     | PolyJet/STRidge  | 0/0 ✓        | 0.0266   |
| 2%    | weak (IBP)    | spectral   | Axial growth  | 🟥       | 0/0 ✓             | 0.0269        | multigraph/RFE         | PDE-FIND/STRidge | 0/0 ✓        | 0.0231   |
| 2%    | weak (IBP)    | polynomial | Radial decay  | 🟥       | 0/0 ✓             | 0.0362        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0310   |
| 2%    | weak (IBP)    | polynomial | Radial growth | 🟥       | 0/0 ✓             | 0.0275        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0239   |
| 2%    | weak (IBP)    | polynomial | Axial decay   | 🟥       | 0/0 ✓             | 0.0317        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0270   |
| 2%    | weak (IBP)    | polynomial | Axial growth  | 🟥       | 0/0 ✓             | 0.0277        | multigraph/RFE         | PDE-FIND/STRidge | 0/0 ✓        | 0.0240   |
| 5%    | weak (no IBP) | spectral   | Radial decay  | 🟢       | 0/0 ✓             | 0.1618        | multigraph/STRidge     | —                | 2/0 ✗        | 0.2238   |
| 5%    | weak (no IBP) | spectral   | Radial growth | 🟩       | 0/0 ✓             | 0.0742        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0756   |
| 5%    | weak (no IBP) | spectral   | Axial decay   | 🟩       | 0/0 ✓             | 0.1218        | multigraph/STRidge     | WSINDy/STRidge   | 0/0 ✓        | 0.1240   |
| 5%    | weak (no IBP) | spectral   | Axial growth  | 🟩       | 0/0 ✓             | 0.0760        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0776   |
| 5%    | weak (no IBP) | polynomial | Radial decay  | 🟢       | 0/0 ✓             | 0.1339        | multigraph/STRidge     | —                | 5/0 ✗        | 0.2723   |
| 5%    | weak (no IBP) | polynomial | Radial growth | 🟩       | 0/0 ✓             | 0.0640        | multigraph/STRidge     | WSINDy/STRidge   | 0/0 ✓        | 0.0654   |
| 5%    | weak (no IBP) | polynomial | Axial decay   | 🟩       | 0/0 ✓             | 0.1018        | mod_multigraph/STRidge | WSINDy/STRidge   | 0/0 ✓        | 0.1046   |
| 5%    | weak (no IBP) | polynomial | Axial growth  | 🟩       | 0/0 ✓             | 0.0654        | mod_multigraph/STRidge | WSINDy/STRidge   | 0/0 ✓        | 0.0670   |
| 5%    | weak (IBP)    | spectral   | Radial decay  | 🟢       | 0/0 ✓             | 0.1619        | multigraph/STRidge     | —                | 2/0 ✗        | 0.2190   |
| 5%    | weak (IBP)    | spectral   | Radial growth | 🟥       | 0/0 ✓             | 0.0745        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0695   |
| 5%    | weak (IBP)    | spectral   | Axial decay   | 🟥       | 0/0 ✓             | 0.1221        | multigraph/STRidge     | WSINDy/STRidge   | 0/0 ✓        | 0.1182   |
| 5%    | weak (IBP)    | spectral   | Axial growth  | 🟥       | 0/0 ✓             | 0.0760        | multigraph/RFE         | PDE-FIND/STRidge | 0/0 ✓        | 0.0715   |
| 5%    | weak (IBP)    | polynomial | Radial decay  | 🟢       | 0/0 ✓             | 0.1383        | multigraph/STRidge     | —                | 5/0 ✗        | 0.2765   |
| 5%    | weak (IBP)    | polynomial | Radial growth | 🟥       | 0/0 ✓             | 0.0663        | multigraph/STRidge     | WSINDy/STRidge   | 0/0 ✓        | 0.0623   |
| 5%    | weak (IBP)    | polynomial | Axial decay   | 🟥       | 0/0 ✓             | 0.1054        | multigraph/STRidge     | WSINDy/STRidge   | 0/0 ✓        | 0.1030   |
| 5%    | weak (IBP)    | polynomial | Axial growth  | 🟥       | 0/0 ✓             | 0.0677        | multigraph/RFE         | PDE-FIND/STRidge | 0/0 ✓        | 0.0639   |

## Ablation

### 2% — Axial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | PolyJet + STRidge              | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0266  | 0.9      | —        |
| 2    | WSINDy + STRidge               | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0266  | 0.4      | +0.0%    |
| 3    | PDE-FIND + STRidge             | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0270  | 0.4      | +1.5%    |
| 4    | PolyJet + RFE                  | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0282  | 0.1      | +6.0%    |
| 5    | WSINDy + RFE                   | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0282  | 0.1      | +6.0%    |
| 6    | PDE-FIND + RFE                 | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0282  | 0.1      | +6.0%    |
| 7    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0316  | 0.2      | +18.7%   |
| 8    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0317  | 0.2      | +19.2%   |
| 9    | Multigraphs + RFE              | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0329  | 0.1      | +23.6%   |
| 10   | Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0329  | 0.1      | +23.6%   |

### 2% — Radial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | PDE-FIND + STRidge             | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0310  | 0.4      | —        |
| 2    | WSINDy + STRidge               | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0310  | 0.4      | +0.0%    |
| 3    | PolyJet + STRidge              | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0310  | 0.8      | +0.0%    |
| 4    | PDE-FIND + RFE                 | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0331  | 0.1      | +6.8%    |
| 5    | WSINDy + RFE                   | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0331  | 0.1      | +6.8%    |
| 6    | PolyJet + RFE                  | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0331  | 0.1      | +6.8%    |
| 7    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0356  | 0.2      | +14.8%   |
| 8    | Multigraphs + STRidge          | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0362  | 0.2      | +16.6%   |
| 9    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0376  | 0.1      | +21.2%   |
| 10   | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0376  | 0.1      | +21.2%   |

### 2% — Axial growth

| Rank | Library + Algorithm            | Configuration              | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:---------------------------|:------|:--------|:---------|:---------|
| 1    | PDE-FIND + STRidge             | spectral, weak (IBP), m=12 | 0/0 ✓ | 0.0231  | 0.4      | —        |
| 2    | WSINDy + STRidge               | spectral, weak (IBP), m=12 | 0/0 ✓ | 0.0231  | 0.4      | +0.0%    |
| 3    | PolyJet + STRidge              | spectral, weak (IBP), m=12 | 0/0 ✓ | 0.0231  | 0.8      | +0.0%    |
| 4    | PolyJet + RFE                  | spectral, weak (IBP), m=12 | 0/0 ✓ | 0.0233  | 0.1      | +0.7%    |
| 5    | PDE-FIND + RFE                 | spectral, weak (IBP), m=12 | 0/0 ✓ | 0.0233  | 0.1      | +0.7%    |
| 6    | WSINDy + RFE                   | spectral, weak (IBP), m=12 | 0/0 ✓ | 0.0233  | 0.1      | +0.7%    |
| 7    | Multigraphs + RFE              | spectral, weak (IBP), m=12 | 0/0 ✓ | 0.0269  | 0.1      | +16.4%   |
| 8    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12 | 0/0 ✓ | 0.0269  | 0.1      | +16.4%   |
| 9    | Multigraphs + STRidge          | spectral, weak (IBP), m=12 | 0/0 ✓ | 0.0270  | 0.2      | +16.7%   |
| 10   | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12 | 0/0 ✓ | 0.0270  | 0.2      | +16.7%   |

### 2% — Radial growth

| Rank | Library + Algorithm            | Configuration                  | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-------------------------------|:------|:--------|:---------|:---------|
| 1    | WSINDy + STRidge               | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0230  | 0.3      | —        |
| 2    | PolyJet + RFE                  | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0232  | 0.1      | +0.9%    |
| 3    | PDE-FIND + RFE                 | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0232  | 0.1      | +0.9%    |
| 4    | WSINDy + RFE                   | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0232  | 0.1      | +0.9%    |
| 5    | PDE-FIND + STRidge             | polynomial, weak (IBP), m=12   | 0/0 ✓ | 0.0239  | 0.4      | +3.8%    |
| 6    | Multigraphs + STRidge          | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0268  | 0.2      | +16.3%   |
| 7    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0268  | 0.2      | +16.3%   |
| 8    | Multigraphs + RFE              | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0268  | 0.1      | +16.5%   |
| 9    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0268  | 0.1      | +16.5%   |
| 10   | PolyJet + STRidge              | polynomial, weak (no IBP), m=8 | 0/0 ✓ | 0.2597  | 1.2      | +1027.7% |

### 5% — Axial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1018  | 0.2      | —        |
| 2    | WSINDy + STRidge               | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.1030  | 0.3      | +1.2%    |
| 3    | Multigraphs + STRidge          | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.1054  | 0.2      | +3.5%    |
| 4    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1074  | 0.1      | +5.5%    |
| 5    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1074  | 0.1      | +5.5%    |
| 6    | PDE-FIND + STRidge             | spectral, weak (IBP), m=12      | 2/0   | 0.1729  | 0.5      | ✗        |
| 7    | PolyJet + STRidge              | spectral, weak (IBP), m=12      | 9/0   | 0.2663  | 1.0      | ✗        |
| 8    | PDE-FIND + RFE                 | polynomial, weak (IBP), m=12    | 0/4   | 0.7214  | 0.1      | ✗        |
| 9    | WSINDy + RFE                   | polynomial, weak (IBP), m=12    | 0/4   | 0.7214  | 0.1      | ✗        |
| 10   | PolyJet + RFE                  | polynomial, weak (IBP), m=12    | 0/4   | 0.7214  | 0.1      | ✗        |

### 5% — Radial decay

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1339  | 0.2      | —        |
| 2    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1339  | 0.2      | +0.0%    |
| 3    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1406  | 0.1      | +5.0%    |
| 4    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.1406  | 0.1      | +5.0%    |
| 5    | PDE-FIND + STRidge             | spectral, weak (IBP), m=12      | 2/0   | 0.2190  | 0.5      | ✗        |
| 6    | PolyJet + STRidge              | spectral, weak (IBP), m=8       | 10/0  | 1.5215  | 1.5      | ✗        |
| 7    | WSINDy + STRidge               | polynomial, weak (IBP), m=12    | 18/0  | 0.9988  | 0.4      | ✗        |
| 8    | WSINDy + RFE                   | polynomial, weak (IBP), m=12    | 0/4   | 0.7199  | 0.1      | ✗        |
| 9    | PDE-FIND + RFE                 | polynomial, weak (IBP), m=12    | 0/4   | 0.7199  | 0.1      | ✗        |
| 10   | PolyJet + RFE                  | polynomial, weak (IBP), m=12    | 0/4   | 0.7199  | 0.1      | ✗        |

### 5% — Axial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | PDE-FIND + STRidge             | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0639  | 0.4      | —        |
| 2    | WSINDy + STRidge               | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0639  | 0.3      | +0.0%    |
| 3    | WSINDy + RFE                   | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0647  | 0.1      | +1.3%    |
| 4    | PDE-FIND + RFE                 | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0647  | 0.1      | +1.3%    |
| 5    | PolyJet + RFE                  | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0647  | 0.1      | +1.3%    |
| 6    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0654  | 0.2      | +2.3%    |
| 7    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0658  | 0.1      | +2.9%    |
| 8    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0658  | 0.1      | +2.9%    |
| 9    | Multigraphs + STRidge          | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0760  | 0.2      | +18.9%   |
| 10   | PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 5/0   | 1.4673  | 0.9      | ✗        |

### 5% — Radial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | WSINDy + STRidge               | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0623  | 0.3      | —        |
| 2    | PolyJet + RFE                  | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0635  | 0.1      | +1.8%    |
| 3    | WSINDy + RFE                   | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0635  | 0.1      | +1.8%    |
| 4    | PDE-FIND + RFE                 | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0635  | 0.1      | +1.8%    |
| 5    | Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0640  | 0.2      | +2.6%    |
| 6    | Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0640  | 0.2      | +2.6%    |
| 7    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0648  | 0.1      | +3.9%    |
| 8    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0648  | 0.1      | +3.9%    |
| 9    | PDE-FIND + STRidge             | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0695  | 0.4      | +11.5%   |
| 10   | PolyJet + STRidge              | spectral, weak (IBP), m=12      | 2/0   | 0.0864  | 0.9      | ✗        |


## Detailed results

### weak (no IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0332  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0381  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0288  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0287  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 7        | 0.0312  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 7        | 0.0354  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0280  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0279  | 0.4      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.7220  | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.7205  | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0684  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0671  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 5/0   | 12       | 0.2109  | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 5/0   | 12       | 0.2723  | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 8  | 0.1       | Axial growth  | 5/0   | 12       | 0.8883  | 0.8      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 8  | 0.1       | Radial growth | 5/0   | 12       | 0.8849  | 0.8      | 🔴      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0332  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0381  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0288  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0287  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 7        | 0.0312  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 7        | 0.0354  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 7        | 0.0280  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 7        | 0.0279  | 0.3      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.7220  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.7205  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0684  | 0.1      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0671  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 7        | 0.1046  | 0.3      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 0.1       | Radial decay  | 18/0  | 13       | 1.0240  | 0.4      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0670  | 0.3      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0654  | 0.3      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0332  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0381  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0288  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0287  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 11/0  | 18       | 2.3217  | 1.0      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 6/0   | 13       | 2.3733  | 0.9      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 3/0   | 10       | 1.1953  | 0.9      | 🔴      |
| PolyJet     | 2%    | STRidge   | 8  | 1.0       | Radial growth | 0/0   | 7        | 0.2597  | 1.2      | 🟢      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.7220  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.7205  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0684  | 0.1      | 🟢      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0671  | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 21/0  | 28       | 3.8400  | 1.2      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 16/0  | 23       | 2.8504  | 1.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 5/0   | 12       | 1.4673  | 0.9      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial growth | 7/0   | 14       | 1.5905  | 0.9      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0332  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0376  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0291  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0290  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.0316  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.0356  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0285  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0283  | 0.2      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.1074  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.1406  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0658  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0648  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.1018  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.1339  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0654  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0640  | 0.2      | 🟢      |

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
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 2/0   | 9        | 0.1779  | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 2/0   | 9        | 0.2238  | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0776  | 0.5      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 7        | 0.0756  | 0.5      | 🟢      |
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
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0774  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 7        | 0.1240  | 0.3      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 0.1       | Radial decay  | 18/0  | 13       | 1.0033  | 0.4      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0776  | 0.3      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0756  | 0.3      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0348  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0411  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0290  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0288  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 7        | 0.0325  | 0.9      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 7        | 0.0381  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0282  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial growth | 6/0   | 13       | 0.1319  | 0.9      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8238  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8220  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0790  | 0.1      | 🟢      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0774  | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 9/0   | 16       | 0.2707  | 1.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 0.1       | Radial decay  | 17/0  | 24       | 0.6322  | 1.3      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 11/0  | 18       | 0.2144  | 0.9      | 🔴      |
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
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0284  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0331  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0242  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0241  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 7        | 0.0270  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 7        | 0.0310  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0240  | 0.4      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0239  | 0.4      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.7214  | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.7199  | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0647  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0635  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 0.1       | Axial decay   | 4/0   | 11       | 0.1896  | 0.6      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 5/0   | 12       | 0.2765  | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0639  | 0.4      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Radial growth | 15/0  | 22       | 0.1489  | 0.6      | 🔴      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0284  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0331  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0242  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0241  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 7        | 0.0270  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 7        | 0.0310  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 7        | 0.0240  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 7        | 0.0239  | 0.3      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.7214  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.7199  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0647  | 0.1      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0635  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 7        | 0.1030  | 0.3      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 0.1       | Radial decay  | 18/0  | 13       | 0.9988  | 0.4      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0639  | 0.3      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0623  | 0.3      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0284  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0331  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0242  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0241  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 4/0   | 11       | 2.0522  | 0.9      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 7        | 0.0310  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 3/0   | 10       | 1.1276  | 0.9      | 🔴      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial growth | 6/0   | 13       | 1.1063  | 0.9      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.7214  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.7199  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0647  | 0.1      | 🟢      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0635  | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 14/0  | 21       | 2.9418  | 1.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 16/0  | 23       | 2.8584  | 1.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 14/0  | 21       | 2.5531  | 1.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial growth | 33/0  | 40       | 3.2988  | 1.4      | 🔴      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0329  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0378  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0277  | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0276  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.0317  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.0362  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0277  | 0.2      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0275  | 0.2      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.1107  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.1447  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0677  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0667  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.1054  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.1383  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0678  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0663  | 0.2      | 🟢      |

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
| PDE-FIND    | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0708  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 2/0   | 9        | 0.1729  | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Radial decay  | 2/0   | 9        | 0.2190  | 0.5      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0715  | 0.4      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0695  | 0.4      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0282  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0341  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0233  | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0232  | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 7        | 0.0266  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 7        | 0.0316  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 7        | 0.0231  | 0.4      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 7        | 0.0230  | 0.3      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8228  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8210  | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0723  | 0.1      | 🟢      |
| WSINDy      | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0708  | 0.1      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 7        | 0.1182  | 0.4      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 0.1       | Radial decay  | 18/0  | 13       | 0.9991  | 0.4      | 🔴      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0715  | 0.3      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 7        | 0.0695  | 0.3      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 7        | 0.0282  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 7        | 0.0341  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0233  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0232  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 7        | 0.0266  | 0.9      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 7        | 0.0316  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 7        | 0.0231  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial growth | 6/0   | 13       | 0.1363  | 0.9      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial decay   | 0/4   | 3        | 0.8228  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial decay  | 0/4   | 3        | 0.8210  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 7        | 0.0723  | 0.1      | 🟢      |
| PolyJet     | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 7        | 0.0708  | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial decay   | 9/0   | 16       | 0.2663  | 1.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Radial decay  | 10/0  | 17       | 1.5215  | 1.5      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 11/0  | 18       | 0.2118  | 0.9      | 🔴      |
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
*Source: `sweep_aniso_20260803_194535_full_local.csv`*

*Run: 775.9 min (Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz), OMP_NUM_THREADS=1*
