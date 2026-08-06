# Allen-Cahn equation — Anisotropic noise

The equation under consideration is

$$f_t = f - f^3 + \varepsilon^2\Delta f, \qquad \varepsilon = 0.15$$

which can be written in terms of $O(2)$-equivariant terms

| Operator   | True coefficient           |
|:-----------|:---------------------------|
| $f$        | $+1.0$                     |
| $f^3$      | $-1.0$                     |
| $\Delta f$ | $+\varepsilon^2 = +0.0225$ |

or standard terms (PolyJet basis)

| Operator          | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|:------------------|:-----------|:---------|:-------|:---------|:-----------|
| $f$               |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f^3$             |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xx} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |

All libraries are evaluated in the common polynomial-jet basis (max order 4,
max degree 3).  fp/fn/E are computed against the same ground truth for every
library.  Allen-Cahn contains only pure monomials ($f$, $f^3$) and pure second derivatives ($f_{xx}$, $f_{yy}$).  All four libraries can represent every term, and this is considered the simplest equation in the benchmark set.

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
| 2%    | weak (no IBP) | spectral   | Radial decay  | 🟩       | 0/0 ✓             | 0.0049        | mod_multigraph/RFE     | PolyJet/RFE      | 0/0 ✓        | 0.0052   |
| 2%    | weak (no IBP) | spectral   | Radial growth | ⚪       | 0/0 ✓             | 0.0046        | mod_multigraph/RFE     | WSINDy/RFE       | 0/0 ✓        | 0.0046   |
| 2%    | weak (no IBP) | spectral   | Axial decay   | 🟩       | 0/0 ✓             | 0.0049        | multigraph/RFE         | PDE-FIND/RFE     | 0/0 ✓        | 0.0050   |
| 2%    | weak (no IBP) | spectral   | Axial growth  | ⚪       | 0/0 ✓             | 0.0045        | multigraph/RFE         | WSINDy/RFE       | 0/0 ✓        | 0.0046   |
| 2%    | weak (no IBP) | polynomial | Radial decay  | 🟩       | 0/0 ✓             | 0.0049        | multigraph/RFE         | WSINDy/RFE       | 0/0 ✓        | 0.0052   |
| 2%    | weak (no IBP) | polynomial | Radial growth | ⚪       | 0/0 ✓             | 0.0046        | multigraph/RFE         | WSINDy/RFE       | 0/0 ✓        | 0.0047   |
| 2%    | weak (no IBP) | polynomial | Axial decay   | 🟩       | 0/0 ✓             | 0.0049        | multigraph/RFE         | PolyJet/RFE      | 0/0 ✓        | 0.0050   |
| 2%    | weak (no IBP) | polynomial | Axial growth  | ⚪       | 0/0 ✓             | 0.0045        | multigraph/RFE         | WSINDy/RFE       | 0/0 ✓        | 0.0046   |
| 2%    | weak (IBP)    | spectral   | Radial decay  | 🟩       | 0/0 ✓             | 0.0049        | mod_multigraph/RFE     | PolyJet/RFE      | 0/0 ✓        | 0.0052   |
| 2%    | weak (IBP)    | spectral   | Radial growth | ⚪       | 0/0 ✓             | 0.0046        | mod_multigraph/RFE     | PolyJet/RFE      | 0/0 ✓        | 0.0046   |
| 2%    | weak (IBP)    | spectral   | Axial decay   | 🟩       | 0/0 ✓             | 0.0049        | multigraph/RFE         | PDE-FIND/RFE     | 0/0 ✓        | 0.0050   |
| 2%    | weak (IBP)    | spectral   | Axial growth  | ⚪       | 0/0 ✓             | 0.0045        | multigraph/RFE         | WSINDy/RFE       | 0/0 ✓        | 0.0046   |
| 2%    | weak (IBP)    | polynomial | Radial decay  | 🟩       | 0/0 ✓             | 0.0049        | multigraph/RFE         | WSINDy/RFE       | 0/0 ✓        | 0.0052   |
| 2%    | weak (IBP)    | polynomial | Radial growth | ⚪       | 0/0 ✓             | 0.0046        | multigraph/RFE         | WSINDy/RFE       | 0/0 ✓        | 0.0046   |
| 2%    | weak (IBP)    | polynomial | Axial decay   | 🟩       | 0/0 ✓             | 0.0049        | multigraph/RFE         | PolyJet/RFE      | 0/0 ✓        | 0.0050   |
| 2%    | weak (IBP)    | polynomial | Axial growth  | ⚪       | 0/0 ✓             | 0.0045        | multigraph/RFE         | WSINDy/RFE       | 0/0 ✓        | 0.0046   |
| 5%    | weak (no IBP) | spectral   | Radial decay  | 🟥       | 0/0 ✓             | 0.0106        | multigraph/STRidge     | WSINDy/STRidge   | 0/0 ✓        | 0.0102   |
| 5%    | weak (no IBP) | spectral   | Radial growth | 🟩       | 0/0 ✓             | 0.0048        | multigraph/RFE         | WSINDy/STRidge   | 0/0 ✓        | 0.0052   |
| 5%    | weak (no IBP) | spectral   | Axial decay   | 🟩       | 0/0 ✓             | 0.0076        | multigraph/RFE         | PDE-FIND/STRidge | 0/0 ✓        | 0.0080   |
| 5%    | weak (no IBP) | spectral   | Axial growth  | 🟩       | 0/0 ✓             | 0.0047        | multigraph/RFE         | PDE-FIND/STRidge | 0/0 ✓        | 0.0051   |
| 5%    | weak (no IBP) | polynomial | Radial decay  | 🟩       | 0/0 ✓             | 0.0100        | multigraph/STRidge     | PDE-FIND/STRidge | 0/0 ✓        | 0.0102   |
| 5%    | weak (no IBP) | polynomial | Radial growth | 🟩       | 0/0 ✓             | 0.0048        | multigraph/RFE         | PDE-FIND/STRidge | 0/0 ✓        | 0.0052   |
| 5%    | weak (no IBP) | polynomial | Axial decay   | 🟩       | 0/0 ✓             | 0.0076        | multigraph/RFE         | PDE-FIND/STRidge | 0/0 ✓        | 0.0080   |
| 5%    | weak (no IBP) | polynomial | Axial growth  | 🟩       | 0/0 ✓             | 0.0047        | multigraph/RFE         | PDE-FIND/STRidge | 0/0 ✓        | 0.0051   |
| 5%    | weak (IBP)    | spectral   | Radial decay  | 🟩       | 0/0 ✓             | 0.0099        | mod_multigraph/STRidge | WSINDy/STRidge   | 0/0 ✓        | 0.0101   |
| 5%    | weak (IBP)    | spectral   | Radial growth | 🟩       | 0/0 ✓             | 0.0048        | multigraph/RFE         | WSINDy/STRidge   | 0/0 ✓        | 0.0052   |
| 5%    | weak (IBP)    | spectral   | Axial decay   | 🟩       | 0/0 ✓             | 0.0076        | multigraph/RFE         | PDE-FIND/STRidge | 0/0 ✓        | 0.0080   |
| 5%    | weak (IBP)    | spectral   | Axial growth  | 🟩       | 0/0 ✓             | 0.0047        | multigraph/RFE         | PDE-FIND/STRidge | 0/0 ✓        | 0.0051   |
| 5%    | weak (IBP)    | polynomial | Radial decay  | 🟩       | 0/0 ✓             | 0.0099        | multigraph/STRidge     | PolyJet/STRidge  | 0/0 ✓        | 0.0101   |
| 5%    | weak (IBP)    | polynomial | Radial growth | 🟩       | 0/0 ✓             | 0.0048        | multigraph/RFE         | PolyJet/STRidge  | 0/0 ✓        | 0.0052   |
| 5%    | weak (IBP)    | polynomial | Axial decay   | 🟩       | 0/0 ✓             | 0.0076        | multigraph/RFE         | PolyJet/STRidge  | 0/0 ✓        | 0.0080   |
| 5%    | weak (IBP)    | polynomial | Axial growth  | 🟩       | 0/0 ✓             | 0.0047        | multigraph/RFE         | PDE-FIND/STRidge | 0/0 ✓        | 0.0051   |

## Ablation

### 2% — Axial decay

| Rank | Library + Algorithm            | Configuration                | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-----------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + RFE              | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0049  | 0.0      | —        |
| 2    | Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0049  | 0.0      | +0.0%    |
| 3    | PolyJet + RFE                  | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0050  | 0.0      | +2.7%    |
| 4    | PDE-FIND + RFE                 | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0050  | 0.1      | +2.7%    |
| 5    | WSINDy + RFE                   | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0050  | 0.0      | +2.7%    |
| 6    | Multigraphs + STRidge          | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0052  | 0.1      | +6.2%    |
| 7    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0052  | 0.1      | +6.2%    |
| 8    | PolyJet + STRidge              | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0053  | 0.7      | +8.4%    |
| 9    | PDE-FIND + STRidge             | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0053  | 0.2      | +8.4%    |
| 10   | WSINDy + STRidge               | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0053  | 0.2      | +8.4%    |

### 2% — Radial decay

| Rank | Library + Algorithm            | Configuration                | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-----------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0049  | 0.0      | —        |
| 2    | Multigraphs + RFE              | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0049  | 0.0      | +0.0%    |
| 3    | PolyJet + RFE                  | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0052  | 0.1      | +5.8%    |
| 4    | WSINDy + RFE                   | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0052  | 0.0      | +5.8%    |
| 5    | PDE-FIND + RFE                 | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0052  | 0.0      | +5.8%    |
| 6    | Multigraphs + STRidge          | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0052  | 0.1      | +6.8%    |
| 7    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0052  | 0.1      | +6.8%    |
| 8    | PDE-FIND + STRidge             | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0054  | 0.2      | +11.0%   |
| 9    | WSINDy + STRidge               | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0054  | 0.2      | +11.0%   |
| 10   | PolyJet + STRidge              | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0054  | 0.7      | +11.0%   |

### 2% — Axial growth

| Rank | Library + Algorithm            | Configuration                | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-----------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + RFE              | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0045  | 0.0      | —        |
| 2    | Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0045  | 0.0      | +0.0%    |
| 3    | WSINDy + RFE                   | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0046  | 0.0      | +2.0%    |
| 4    | PDE-FIND + RFE                 | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0046  | 0.0      | +2.0%    |
| 5    | PolyJet + RFE                  | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0046  | 0.0      | +2.0%    |
| 6    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0046  | 0.1      | +3.0%    |
| 7    | Multigraphs + STRidge          | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0046  | 0.1      | +3.0%    |
| 8    | PolyJet + STRidge              | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0047  | 0.8      | +4.5%    |
| 9    | PDE-FIND + STRidge             | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0047  | 0.2      | +4.5%    |
| 10   | WSINDy + STRidge               | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0047  | 0.2      | +4.5%    |

### 2% — Radial growth

| Rank | Library + Algorithm            | Configuration                | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-----------------------------|:------|:--------|:---------|:---------|
| 1    | Modified Multigraphs + RFE     | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0046  | 0.0      | —        |
| 2    | Multigraphs + RFE              | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0046  | 0.0      | +0.0%    |
| 3    | PolyJet + RFE                  | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0046  | 0.1      | +1.6%    |
| 4    | WSINDy + RFE                   | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0046  | 0.0      | +1.6%    |
| 5    | PDE-FIND + RFE                 | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0046  | 0.0      | +1.6%    |
| 6    | Multigraphs + STRidge          | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0047  | 0.1      | +3.2%    |
| 7    | Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0047  | 0.1      | +3.2%    |
| 8    | PDE-FIND + STRidge             | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.2      | +4.5%    |
| 9    | WSINDy + STRidge               | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.2      | +4.5%    |
| 10   | PolyJet + STRidge              | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.7      | +4.5%    |

### 5% — Axial decay

| Rank | Library + Algorithm            | Configuration                  | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + RFE              | polynomial, weak (no IBP), m=8 | 0/0 ✓ | 0.0076  | 0.1      | —        |
| 2    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=8 | 0/0 ✓ | 0.0076  | 0.0      | +0.0%    |
| 3    | Multigraphs + STRidge          | spectral, weak (IBP), m=8      | 0/0 ✓ | 0.0079  | 0.2      | +4.0%    |
| 4    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=8      | 0/0 ✓ | 0.0079  | 0.1      | +4.0%    |
| 5    | PolyJet + STRidge              | polynomial, weak (IBP), m=8    | 0/0 ✓ | 0.0080  | 0.6      | +5.6%    |
| 6    | PDE-FIND + STRidge             | spectral, weak (IBP), m=8      | 0/0 ✓ | 0.0080  | 0.3      | +5.6%    |
| 7    | WSINDy + STRidge               | spectral, weak (IBP), m=8      | 0/0 ✓ | 0.0080  | 0.2      | +5.6%    |
| 8    | PDE-FIND + RFE                 | polynomial, weak (no IBP), m=8 | 0/2   | 0.3555  | 0.0      | ✗        |
| 9    | PolyJet + RFE                  | spectral, weak (no IBP), m=8   | 0/2   | 0.3555  | 0.1      | ✗        |
| 10   | WSINDy + RFE                   | spectral, weak (no IBP), m=8   | 0/2   | 0.3555  | 0.0      | ✗        |

### 5% — Radial decay

| Rank | Library + Algorithm            | Configuration                  | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:-------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + STRidge          | polynomial, weak (IBP), m=8    | 0/0 ✓ | 0.0099  | 0.1      | —        |
| 2    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=8      | 0/0 ✓ | 0.0099  | 0.1      | +0.0%    |
| 3    | PolyJet + STRidge              | polynomial, weak (IBP), m=8    | 0/0 ✓ | 0.0101  | 0.6      | +2.6%    |
| 4    | PDE-FIND + STRidge             | polynomial, weak (IBP), m=8    | 0/0 ✓ | 0.0101  | 0.3      | +2.6%    |
| 5    | WSINDy + STRidge               | spectral, weak (IBP), m=8      | 0/0 ✓ | 0.0101  | 0.2      | +2.6%    |
| 6    | WSINDy + RFE                   | spectral, weak (no IBP), m=8   | 0/2   | 0.3556  | 0.0      | ✗        |
| 7    | PDE-FIND + RFE                 | spectral, weak (no IBP), m=8   | 0/2   | 0.3556  | 0.0      | ✗        |
| 8    | PolyJet + RFE                  | polynomial, weak (no IBP), m=8 | 0/2   | 0.3556  | 0.1      | ✗        |
| 9    | Multigraphs + RFE              | spectral, weak (no IBP), m=8   | 0/2   | 0.3556  | 0.0      | ✗        |
| 10   | Modified Multigraphs + RFE     | spectral, weak (no IBP), m=8   | 0/2   | 0.3556  | 0.0      | ✗        |

### 5% — Axial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0047  | 0.0      | —        |
| 2    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0047  | 0.0      | +0.0%    |
| 3    | Multigraphs + STRidge          | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0050  | 0.1      | +6.3%    |
| 4    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0050  | 0.1      | +6.3%    |
| 5    | PDE-FIND + STRidge             | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0051  | 0.3      | +7.2%    |
| 6    | WSINDy + STRidge               | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0051  | 0.2      | +7.2%    |
| 7    | PolyJet + STRidge              | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0051  | 0.5      | +7.2%    |
| 8    | PDE-FIND + RFE                 | polynomial, weak (no IBP), m=8  | 0/2   | 0.3555  | 0.0      | ✗        |
| 9    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/2   | 0.3555  | 0.1      | ✗        |
| 10   | WSINDy + RFE                   | spectral, weak (no IBP), m=8    | 0/2   | 0.3555  | 0.0      | ✗        |

### 5% — Radial growth

| Rank | Library + Algorithm            | Configuration                   | fp/fn | E_coeff | Time (s) | vs. best |
|:-----|:-------------------------------|:--------------------------------|:------|:--------|:---------|:---------|
| 1    | Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0048  | 0.1      | —        |
| 2    | Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0048  | 0.0      | +0.0%    |
| 3    | Multigraphs + STRidge          | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0052  | 0.1      | +6.9%    |
| 4    | Modified Multigraphs + STRidge | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0052  | 0.1      | +6.9%    |
| 5    | PolyJet + STRidge              | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0052  | 0.4      | +7.7%    |
| 6    | WSINDy + STRidge               | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0052  | 0.2      | +7.7%    |
| 7    | PDE-FIND + STRidge             | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0052  | 0.2      | +7.7%    |
| 8    | PDE-FIND + RFE                 | polynomial, weak (no IBP), m=8  | 0/2   | 0.3555  | 0.0      | ✗        |
| 9    | PolyJet + RFE                  | spectral, weak (no IBP), m=8    | 0/2   | 0.3555  | 0.1      | ✗        |
| 10   | WSINDy + RFE                   | spectral, weak (no IBP), m=8    | 0/2   | 0.3555  | 0.0      | ✗        |


## Detailed results

### weak (no IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0050  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0052  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0047  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 4        | 0.0053  | 0.3      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.0054  | 0.3      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 4        | 0.0047  | 0.3      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 4        | 0.0048  | 0.3      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 0/0   | 4        | 0.0080  | 0.3      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 0/0   | 4        | 0.0102  | 0.3      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 4        | 0.0051  | 0.2      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 4        | 0.0052  | 0.2      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0050  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0052  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0047  | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.01      | Axial decay   | 0/0   | 4        | 0.0053  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 4        | 0.0054  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 4        | 0.0047  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 4        | 0.0048  | 0.2      | 🟢      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 0/0   | 4        | 0.0080  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 8  | 0.1       | Radial decay  | 0/0   | 4        | 0.0102  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 4        | 0.0051  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 4        | 0.0052  | 0.2      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0050  | 0.0      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0052  | 0.0      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0047  | 0.0      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.0053  | 0.5      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.0054  | 0.5      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0047  | 0.5      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0048  | 0.5      | 🟢      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 2        | 0.3555  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/2   | 2        | 0.3555  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/2   | 2        | 0.3555  | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 0/0   | 4        | 0.0080  | 0.7      | 🟢      |
| PolyJet     | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 0/0   | 4        | 0.0102  | 0.6      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0051  | 0.5      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0052  | 0.4      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.0049  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.0049  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.0045  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.0046  | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 3        | 0.0052  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 3        | 0.0052  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 3        | 0.0047  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 3        | 0.0047  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 3        | 0.0076  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.0047  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.0048  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 0/0   | 3        | 0.0079  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 0.1       | Radial decay  | 0/0   | 3        | 0.0100  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 3        | 0.0050  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 3        | 0.0052  | 0.1      | 🟢      |

### weak (no IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0050  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0052  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.0053  | 0.2      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 4        | 0.0054  | 0.2      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 4        | 0.0047  | 0.2      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 4        | 0.0048  | 0.2      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 0/0   | 4        | 0.0080  | 0.3      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 4        | 0.0114  | 0.3      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 4        | 0.0051  | 0.2      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 8  | 0.1       | Radial growth | 0/0   | 4        | 0.0055  | 0.3      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0050  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0052  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 4        | 0.0053  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 4        | 0.0054  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.01      | Axial growth  | 0/0   | 4        | 0.0047  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 4        | 0.0048  | 0.2      | 🟢      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 0/0   | 4        | 0.0080  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 8  | 0.1       | Radial decay  | 0/0   | 4        | 0.0102  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 4        | 0.0051  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0052  | 0.2      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0050  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0052  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0046  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0046  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.0053  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.0054  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0047  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0048  | 0.7      | 🟢      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 2        | 0.3555  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/2   | 2        | 0.3555  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/2   | 2        | 0.3555  | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 4        | 0.0086  | 0.5      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 4        | 0.0114  | 0.5      | 🟢      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Axial growth  | 0/0   | 4        | 0.0054  | 0.7      | 🟢      |
| PolyJet     | 5%    | STRidge   | 8  | 1.0       | Radial growth | 0/0   | 4        | 0.0055  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.0049  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.0049  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.0045  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.0046  | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 3        | 0.0052  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 3        | 0.0052  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.01      | Axial growth  | 0/0   | 3        | 0.0046  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 3        | 0.0047  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 3        | 0.0076  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.0047  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.0048  | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 0/0   | 3        | 0.0079  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 3        | 0.0106  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 3        | 0.0050  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 3        | 0.0052  | 0.1      | 🟢      |

### weak (IBP) — Polynomial stencils

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0050  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0052  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 4        | 0.0053  | 0.3      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.0054  | 0.3      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 4        | 0.0047  | 0.3      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0048  | 0.3      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 8  | 1.0       | Axial decay   | 0/0   | 4        | 0.0080  | 0.3      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 8  | 0.1       | Radial decay  | 0/0   | 4        | 0.0101  | 0.3      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 4        | 0.0051  | 0.3      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 4        | 0.0052  | 0.2      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0050  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0052  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 4        | 0.0053  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 4        | 0.0054  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.01      | Axial growth  | 0/0   | 4        | 0.0047  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 4        | 0.0048  | 0.2      | 🟢      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 0/0   | 4        | 0.0080  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 8  | 0.1       | Radial decay  | 0/0   | 4        | 0.0101  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 4        | 0.0051  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0052  | 0.2      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0050  | 0.0      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0052  | 0.0      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 4        | 0.0053  | 0.6      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.0054  | 0.5      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0047  | 0.5      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 4        | 0.0048  | 0.5      | 🟢      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 2        | 0.3555  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/2   | 2        | 0.3555  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/2   | 2        | 0.3555  | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 8  | 1.0       | Axial decay   | 0/0   | 4        | 0.0080  | 0.6      | 🟢      |
| PolyJet     | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 0/0   | 4        | 0.0101  | 0.6      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0051  | 0.5      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0052  | 0.4      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.0049  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.0049  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.0045  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.0046  | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 3        | 0.0052  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 3        | 0.0052  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 3        | 0.0046  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 3        | 0.0047  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 3        | 0.0076  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.0047  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.0048  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 1.0       | Axial decay   | 0/0   | 3        | 0.0079  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 0/0   | 3        | 0.0099  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 3        | 0.0050  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 3        | 0.0052  | 0.1      | 🟢      |

### weak (IBP) — Spectral derivatives

| Library     | Noise | Algorithm | m  | d_tol/lam | Profile       | fp/fn | n_active | E_coeff | Time (s) | Exact? |
|:------------|:------|:----------|:---|:----------|:--------------|:------|:---------|:--------|:---------|:-------|
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0050  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0052  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.0053  | 0.2      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 4        | 0.0054  | 0.2      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 4        | 0.0047  | 0.2      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 4        | 0.0048  | 0.2      | 🟢      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 0/0   | 4        | 0.0080  | 0.3      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 4        | 0.0114  | 0.3      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 4        | 0.0051  | 0.2      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 8  | 0.1       | Radial growth | 0/0   | 4        | 0.0055  | 0.3      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0050  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0052  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| WSINDy      | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0046  | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 4        | 0.0053  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial decay  | 0/0   | 4        | 0.0054  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.01      | Axial growth  | 0/0   | 4        | 0.0047  | 0.2      | 🟢      |
| WSINDy      | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 4        | 0.0048  | 0.2      | 🟢      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| WSINDy      | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/2   | 2        | 0.3555  | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 0/0   | 4        | 0.0080  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 8  | 0.1       | Radial decay  | 0/0   | 4        | 0.0101  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 4        | 0.0051  | 0.2      | 🟢      |
| WSINDy      | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0052  | 0.2      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 4        | 0.0050  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 4        | 0.0052  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 4        | 0.0046  | 0.1      | 🟢      |
| PolyJet     | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 4        | 0.0046  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial decay   | 0/0   | 4        | 0.0053  | 0.7      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 4        | 0.0054  | 0.7      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Axial growth  | 0/0   | 4        | 0.0047  | 0.8      | 🟢      |
| PolyJet     | 2%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 4        | 0.0048  | 0.7      | 🟢      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/2   | 2        | 0.3555  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Axial growth  | 0/2   | 2        | 0.3555  | 0.1      | 🔴      |
| PolyJet     | 5%    | RFE       | 8  | 0.0       | Radial growth | 0/2   | 2        | 0.3555  | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Axial decay   | 0/0   | 4        | 0.0086  | 0.5      | 🟢      |
| PolyJet     | 5%    | STRidge   | 12 | 10.0      | Radial decay  | 0/0   | 4        | 0.0114  | 0.5      | 🟢      |
| PolyJet     | 5%    | STRidge   | 8  | 0.1       | Axial growth  | 0/0   | 4        | 0.0054  | 0.7      | 🟢      |
| PolyJet     | 5%    | STRidge   | 8  | 1.0       | Radial growth | 0/0   | 4        | 0.0055  | 0.6      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial decay   | 0/0   | 3        | 0.0049  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial decay  | 0/0   | 3        | 0.0049  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.0045  | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.0046  | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Axial decay   | 0/0   | 3        | 0.0052  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 1.0       | Radial decay  | 0/0   | 3        | 0.0052  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 3        | 0.0046  | 0.1      | 🟢      |
| Multigraphs | 2%    | STRidge   | 12 | 0.1       | Radial growth | 0/0   | 3        | 0.0047  | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Axial decay   | 0/0   | 3        | 0.0076  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 8  | 0.0       | Radial decay  | 0/2   | 2        | 0.3556  | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Axial growth  | 0/0   | 3        | 0.0047  | 0.0      | 🟢      |
| Multigraphs | 5%    | RFE       | 12 | 0.0       | Radial growth | 0/0   | 3        | 0.0048  | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 0.1       | Axial decay   | 0/0   | 3        | 0.0079  | 0.2      | 🟢      |
| Multigraphs | 5%    | STRidge   | 8  | 1.0       | Radial decay  | 0/0   | 3        | 0.0099  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 0.1       | Axial growth  | 0/0   | 3        | 0.0050  | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 12 | 1.0       | Radial growth | 0/0   | 3        | 0.0052  | 0.1      | 🟢      |


---
*Source: `sweep_aniso_default_ac_20260805_215410.csv`*

