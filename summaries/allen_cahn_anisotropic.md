# Allen-Cahn equation — Anisotropic noise

The equation under consideration is

$$f_t = f - f^3 + \varepsilon^2\Delta f, \qquad \varepsilon = 0.15$$

which can be written in terms of $O(2)$-equivariant terms

| Operator | True coefficient |
|---|---|
| $f$ | $+1.0$ |
| $f^3$ | $-1.0$ |
| $\Delta f$ | $+\varepsilon^2 = +0.0225$ |

or standard terms (PolyJet basis)

| Operator | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|---|---|---|---|---|---|
| $f$ | | 🟢 | 🟢 | 🟢 | 🟢 |
| $f^3$ | | 🟢 | 🟢 | 🟢 | 🟢 |
| $\partial_{xx} f$ | | 🟢 | 🟢 | 🟢 | 🟢 |
| $\partial_{yy} f$ | | 🟢 | 🟢 | 🟢 | 🟢 |

All libraries are evaluated in the common polynomial-jet basis (105 monomials,
max order 4, max degree 3).  fp/fn/E are computed against the same ground truth
for every library.  Allen-Cahn contains only pure monomials ($f$, $f^3$) and pure second derivatives ($f_{xx}$, $f_{yy}$).  All four libraries can represent every term, so the structural advantage of the equivariant basis is less pronounced here; differences come from numerical conditioning, not representability.

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
3 algorithm(s) (mstls, sr3, strridge),
3 method(s) (pointwise, weak (no IBP), weak (IBP)),
2 backend(s) (polynomial, spectral),
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


## Summary: equivariant vs best non-equivariant (anisotropic noise)

| Noise | Method | Backend | Profile | Eq best | Multigraphs fp/fn | Multigraphs E | Non-eq best | Non-eq fp/fn | Non-eq E | MWSINDy fp/fn | MWSINDy E |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2% | pointwise | spectral | Radial decay | ⚪ | 0/2 | 0.7620 | PDE-FIND | 0/2 | 0.7620 | — | — |
| 2% | pointwise | spectral | Radial growth | ⚪ | 0/2 | 0.7642 | PDE-FIND | 0/2 | 0.7642 | — | — |
| 2% | pointwise | spectral | Axial decay | 🔴 | 0/2 | 0.7631 | PolyJet | 0/1 | 0.7395 | — | — |
| 2% | pointwise | spectral | Axial growth | ⚪ | 0/2 | 0.7640 | PDE-FIND | 0/2 | 0.7640 | — | — |
| 2% | pointwise | polynomial | Radial decay | ⚪ | 0/0 ✓ | 0.0656 | PDE-FIND | 0/0 ✓ | 0.0656 | — | — |
| 2% | pointwise | polynomial | Radial growth | ⚪ | 0/0 ✓ | 0.0601 | PDE-FIND | 0/0 ✓ | 0.0601 | — | — |
| 2% | pointwise | polynomial | Axial decay | ⚪ | 0/0 ✓ | 0.0587 | PDE-FIND | 0/0 ✓ | 0.0587 | — | — |
| 2% | pointwise | polynomial | Axial growth | ⚪ | 0/0 ✓ | 0.0622 | PolyJet | 0/0 ✓ | 0.0623 | — | — |
| 2% | weak (no IBP) | spectral | Radial decay | ⚪ | 0/0 ✓ | 0.0044 | PDE-FIND | 0/0 ✓ | 0.0045 | — | — |
| 2% | weak (no IBP) | spectral | Radial growth | ⚪ | 0/0 ✓ | 0.0046 | PDE-FIND | 0/0 ✓ | 0.0046 | — | — |
| 2% | weak (no IBP) | spectral | Axial decay | ⚪ | 0/0 ✓ | 0.0049 | PDE-FIND | 0/0 ✓ | 0.0050 | — | — |
| 2% | weak (no IBP) | spectral | Axial growth | ⚪ | 0/0 ✓ | 0.0044 | PolyJet | 0/0 ✓ | 0.0045 | — | — |
| 2% | weak (no IBP) | polynomial | Radial decay | ⚪ | 0/0 ✓ | 0.0042 | PDE-FIND | 0/0 ✓ | 0.0043 | — | — |
| 2% | weak (no IBP) | polynomial | Radial growth | ⚪ | 0/0 ✓ | 0.0044 | WSINDy | 0/0 ✓ | 0.0045 | — | — |
| 2% | weak (no IBP) | polynomial | Axial decay | ⚪ | 0/0 ✓ | 0.0047 | WSINDy | 0/0 ✓ | 0.0048 | — | — |
| 2% | weak (no IBP) | polynomial | Axial growth | ⚪ | 0/0 ✓ | 0.0042 | WSINDy | 0/0 ✓ | 0.0042 | — | — |
| 2% | weak (IBP) | spectral | Radial decay | ⚪ | 0/0 ✓ | 0.0044 | PDE-FIND | 0/0 ✓ | 0.0045 | — | — |
| 2% | weak (IBP) | spectral | Radial growth | ⚪ | 0/0 ✓ | 0.0046 | PolyJet | 0/0 ✓ | 0.0046 | — | — |
| 2% | weak (IBP) | spectral | Axial decay | ⚪ | 0/0 ✓ | 0.0049 | PolyJet | 0/0 ✓ | 0.0050 | — | — |
| 2% | weak (IBP) | spectral | Axial growth | ⚪ | 0/0 ✓ | 0.0044 | PolyJet | 0/0 ✓ | 0.0045 | — | — |
| 2% | weak (IBP) | polynomial | Radial decay | ⚪ | 0/0 ✓ | 0.0044 | PDE-FIND | 0/0 ✓ | 0.0045 | — | — |
| 2% | weak (IBP) | polynomial | Radial growth | ⚪ | 0/0 ✓ | 0.0046 | PolyJet | 0/0 ✓ | 0.0046 | — | — |
| 2% | weak (IBP) | polynomial | Axial decay | ⚪ | 0/0 ✓ | 0.0049 | PolyJet | 0/0 ✓ | 0.0050 | — | — |
| 2% | weak (IBP) | polynomial | Axial growth | ⚪ | 0/0 ✓ | 0.0044 | PolyJet | 0/0 ✓ | 0.0045 | — | — |
| 5% | pointwise | spectral | Radial decay | ⚪ | 0/2 | 0.9499 | PDE-FIND | 0/2 | 0.9499 | — | — |
| 5% | pointwise | spectral | Radial growth | 🔴 | 0/2 | 0.9543 | PolyJet | 0/1 | 0.9406 | — | — |
| 5% | pointwise | spectral | Axial decay | 🔴 | 0/2 | 0.9526 | WSINDy | 0/1 | 0.9306 | — | — |
| 5% | pointwise | spectral | Axial growth | 🔴 | 0/2 | 0.9540 | PolyJet | 0/1 | 0.9400 | — | — |
| 5% | pointwise | polynomial | Radial decay | 🟩 | 0/0 ✓ | 0.3492 | PolyJet | 0/0 ✓ | 0.3495 | — | — |
| 5% | pointwise | polynomial | Radial growth | 🟩 | 0/0 ✓ | 0.3293 | PDE-FIND | 0/0 ✓ | 0.3297 | — | — |
| 5% | pointwise | polynomial | Axial decay | ⚪ | 0/0 ✓ | 0.3271 | PDE-FIND | 0/0 ✓ | 0.3272 | — | — |
| 5% | pointwise | polynomial | Axial growth | 🟩 | 0/0 ✓ | 0.3372 | PDE-FIND | 0/0 ✓ | 0.3377 | — | — |
| 5% | weak (no IBP) | spectral | Radial decay | 🟩 | 0/0 ✓ | 0.0041 | PolyJet | 0/0 ✓ | 0.0057 | — | — |
| 5% | weak (no IBP) | spectral | Radial growth | ⚪ | 0/0 ✓ | 0.0057 | PolyJet | 0/0 ✓ | 0.0057 | — | — |
| 5% | weak (no IBP) | spectral | Axial decay | 🟩 | 0/0 ✓ | 0.0059 | PDE-FIND | 0/0 ✓ | 0.0073 | — | — |
| 5% | weak (no IBP) | spectral | Axial growth | ⚪ | 0/0 ✓ | 0.0050 | PolyJet | 0/0 ✓ | 0.0050 | — | — |
| 5% | weak (no IBP) | polynomial | Radial decay | 🟩 | 0/0 ✓ | 0.0038 | PDE-FIND | 0/0 ✓ | 0.0053 | — | — |
| 5% | weak (no IBP) | polynomial | Radial growth | ⚪ | 0/0 ✓ | 0.0055 | PDE-FIND | 0/0 ✓ | 0.0055 | — | — |
| 5% | weak (no IBP) | polynomial | Axial decay | 🟩 | 0/0 ✓ | 0.0059 | WSINDy | 0/0 ✓ | 0.0069 | — | — |
| 5% | weak (no IBP) | polynomial | Axial growth | ⚪ | 0/0 ✓ | 0.0044 | PolyJet | 0/0 ✓ | 0.0044 | — | — |
| 5% | weak (IBP) | spectral | Radial decay | 🟩 | 0/0 ✓ | 0.0041 | PDE-FIND | 0/0 ✓ | 0.0057 | — | — |
| 5% | weak (IBP) | spectral | Radial growth | ⚪ | 0/0 ✓ | 0.0057 | PolyJet | 0/0 ✓ | 0.0057 | — | — |
| 5% | weak (IBP) | spectral | Axial decay | 🟩 | 0/0 ✓ | 0.0059 | PolyJet | 0/0 ✓ | 0.0073 | — | — |
| 5% | weak (IBP) | spectral | Axial growth | ⚪ | 0/0 ✓ | 0.0051 | PDE-FIND | 0/0 ✓ | 0.0051 | — | — |
| 5% | weak (IBP) | polynomial | Radial decay | 🟩 | 0/0 ✓ | 0.0041 | PDE-FIND | 0/0 ✓ | 0.0057 | — | — |
| 5% | weak (IBP) | polynomial | Radial growth | ⚪ | 0/0 ✓ | 0.0057 | PolyJet | 0/0 ✓ | 0.0057 | — | — |
| 5% | weak (IBP) | polynomial | Axial decay | 🟩 | 0/0 ✓ | 0.0059 | PolyJet | 0/0 ✓ | 0.0073 | — | — |
| 5% | weak (IBP) | polynomial | Axial growth | ⚪ | 0/0 ✓ | 0.0051 | PDE-FIND | 0/0 ✓ | 0.0051 | — | — |

## Pointwise formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol/α | Profile | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|---|
| PDE-FIND | 2% | STRidge | 0.2 | Radial decay | 0/2 | 2 | 0.7620 | 0.3 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.2 | Radial growth | 0/2 | 2 | 0.7642 | 0.3 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.2 | Axial decay | 0/2 | 2 | 0.7631 | 0.3 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.2 | Axial growth | 0/2 | 2 | 0.7640 | 0.3 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.05 | Radial decay | 0/2 | 2 | 0.7759 | 0.1 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.05 | Radial growth | 0/2 | 2 | 0.7776 | 0.1 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.05 | Axial decay | 0/2 | 2 | 0.7765 | 0.1 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.05 | Axial growth | 0/2 | 2 | 0.7775 | 0.1 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | Radial decay | 0/2 | 2 | 0.7621 | 0.7 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | Radial growth | 0/2 | 2 | 0.7642 | 0.7 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | Axial decay | 0/2 | 2 | 0.7631 | 0.7 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | Axial growth | 0/2 | 2 | 0.7640 | 0.7 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.2 | Radial decay | 0/2 | 2 | 0.9499 | 0.3 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.2 | Radial growth | 0/2 | 2 | 0.9543 | 0.3 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.2 | Axial decay | 0/2 | 2 | 0.9526 | 0.3 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.2 | Axial growth | 0/2 | 2 | 0.9540 | 0.3 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.01 | Radial decay | 0/2 | 2 | 0.9526 | 0.1 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.01 | Radial growth | 0/2 | 2 | 0.9558 | 0.1 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.01 | Axial decay | 0/2 | 2 | 0.9544 | 0.1 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.01 | Axial growth | 0/2 | 2 | 0.9555 | 0.1 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | Radial decay | 0/2 | 2 | 0.9501 | 0.8 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | Radial growth | 0/2 | 2 | 0.9547 | 0.8 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | Axial decay | 0/2 | 2 | 0.9527 | 0.8 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | Axial growth | 0/2 | 2 | 0.9544 | 0.8 | 🔴 |
| WSINDy | 2% | STRidge | 0.2 | Radial decay | 0/2 | 2 | 0.7620 | 0.2 | 🔴 |
| WSINDy | 2% | STRidge | 0.2 | Radial growth | 0/2 | 2 | 0.7642 | 0.2 | 🔴 |
| WSINDy | 2% | STRidge | 0.2 | Axial decay | 0/2 | 2 | 0.7631 | 0.2 | 🔴 |
| WSINDy | 2% | STRidge | 0.2 | Axial growth | 0/2 | 2 | 0.7640 | 0.2 | 🔴 |
| WSINDy | 2% | SR3 | 0.05 | Radial decay | 0/2 | 2 | 0.7721 | 0.1 | 🔴 |
| WSINDy | 2% | SR3 | 0.05 | Radial growth | 0/2 | 2 | 0.7739 | 0.1 | 🔴 |
| WSINDy | 2% | SR3 | 0.05 | Axial decay | 0/2 | 2 | 0.7728 | 0.1 | 🔴 |
| WSINDy | 2% | SR3 | 0.05 | Axial growth | 0/2 | 2 | 0.7738 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | Radial decay | 0/2 | 2 | 0.7621 | 0.4 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | Radial growth | 0/2 | 2 | 0.7642 | 0.4 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | Axial decay | 0/2 | 2 | 0.7631 | 0.4 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | Axial growth | 0/2 | 2 | 0.7640 | 0.4 | 🔴 |
| WSINDy | 5% | STRidge | 0.2 | Radial decay | 0/2 | 2 | 0.9499 | 0.2 | 🔴 |
| WSINDy | 5% | STRidge | 0.2 | Radial growth | 0/2 | 2 | 0.9543 | 0.2 | 🔴 |
| WSINDy | 5% | STRidge | 0.2 | Axial decay | 0/2 | 2 | 0.9526 | 0.2 | 🔴 |
| WSINDy | 5% | STRidge | 0.2 | Axial growth | 0/2 | 2 | 0.9540 | 0.2 | 🔴 |
| WSINDy | 5% | SR3 | 0.01 | Radial decay | 0/2 | 2 | 0.9519 | 0.1 | 🔴 |
| WSINDy | 5% | SR3 | 0.01 | Radial growth | 0/2 | 2 | 0.9553 | 0.1 | 🔴 |
| WSINDy | 5% | SR3 | 0.005 | Axial decay | 0/1 | 3 | 0.9306 | 0.1 | 🔴 |
| WSINDy | 5% | SR3 | 0.01 | Axial growth | 0/2 | 2 | 0.9550 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | Radial decay | 0/2 | 2 | 0.9501 | 0.5 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | Radial growth | 0/2 | 2 | 0.9547 | 0.5 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | Axial decay | 0/2 | 2 | 0.9527 | 0.5 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | Axial growth | 0/2 | 2 | 0.9544 | 0.4 | 🔴 |
| PolyJet | 2% | STRidge | 0.2 | Radial decay | 0/2 | 2 | 0.7620 | 0.5 | 🔴 |
| PolyJet | 2% | STRidge | 0.2 | Radial growth | 0/2 | 2 | 0.7642 | 0.5 | 🔴 |
| PolyJet | 2% | STRidge | 0.2 | Axial decay | 0/2 | 2 | 0.7631 | 0.5 | 🔴 |
| PolyJet | 2% | STRidge | 0.2 | Axial growth | 0/2 | 2 | 0.7640 | 0.5 | 🔴 |
| PolyJet | 2% | SR3 | 0.005 | Radial decay | 0/2 | 2 | 0.7678 | 0.1 | 🔴 |
| PolyJet | 2% | SR3 | 0.005 | Radial growth | 0/2 | 2 | 0.7688 | 0.1 | 🔴 |
| PolyJet | 2% | SR3 | 0.005 | Axial decay | 0/1 | 3 | 0.7395 | 0.1 | 🔴 |
| PolyJet | 2% | SR3 | 0.005 | Axial growth | 0/2 | 2 | 0.7686 | 0.1 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | Radial decay | 0/2 | 2 | 0.7621 | 1.3 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | Radial growth | 0/2 | 2 | 0.7642 | 1.3 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | Axial decay | 0/2 | 2 | 0.7631 | 1.3 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | Axial growth | 0/2 | 2 | 0.7640 | 1.3 | 🔴 |
| PolyJet | 5% | STRidge | 0.2 | Radial decay | 0/2 | 2 | 0.9499 | 0.5 | 🔴 |
| PolyJet | 5% | STRidge | 0.5 | Radial growth | 0/2 | 2 | 0.9543 | 0.5 | 🔴 |
| PolyJet | 5% | STRidge | 0.2 | Axial decay | 0/2 | 2 | 0.9526 | 0.5 | 🔴 |
| PolyJet | 5% | STRidge | 0.2 | Axial growth | 0/2 | 2 | 0.9540 | 0.5 | 🔴 |
| PolyJet | 5% | SR3 | 0.005 | Radial decay | 0/2 | 2 | 0.9511 | 0.1 | 🔴 |
| PolyJet | 5% | SR3 | 0.005 | Radial growth | 0/1 | 3 | 0.9406 | 0.1 | 🔴 |
| PolyJet | 5% | SR3 | 0.005 | Axial decay | 0/1 | 3 | 0.9326 | 0.1 | 🔴 |
| PolyJet | 5% | SR3 | 0.005 | Axial growth | 0/1 | 3 | 0.9400 | 0.1 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | Radial decay | 0/2 | 2 | 0.9501 | 1.6 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | Radial growth | 0/2 | 2 | 0.9547 | 1.5 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | Axial decay | 0/2 | 2 | 0.9527 | 1.5 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | Axial growth | 0/2 | 2 | 0.9544 | 1.5 | 🔴 |
| Multigraphs | 2% | STRidge | 0.2 | Radial decay | 0/2 | 2 | 0.7620 | 0.1 | 🔴 |
| Multigraphs | 2% | STRidge | 0.2 | Radial growth | 0/2 | 2 | 0.7642 | 0.1 | 🔴 |
| Multigraphs | 2% | STRidge | 0.2 | Axial decay | 0/2 | 2 | 0.7631 | 0.1 | 🔴 |
| Multigraphs | 2% | STRidge | 0.2 | Axial growth | 0/2 | 2 | 0.7640 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | Radial decay | 0/2 | 2 | 0.7696 | 0.0 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | Radial growth | 0/2 | 2 | 0.7714 | 0.0 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | Axial decay | 0/2 | 2 | 0.7704 | 0.0 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | Axial growth | 0/2 | 2 | 0.7712 | 0.0 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | Radial decay | 0/2 | 2 | 0.7621 | 0.1 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | Radial growth | 0/2 | 2 | 0.7642 | 0.1 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | Axial decay | 0/2 | 2 | 0.7631 | 0.2 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | Axial growth | 0/2 | 2 | 0.7640 | 0.1 | 🔴 |
| Multigraphs | 5% | STRidge | 0.2 | Radial decay | 0/2 | 2 | 0.9499 | 0.1 | 🔴 |
| Multigraphs | 5% | STRidge | 0.1 | Radial growth | 0/2 | 2 | 0.9543 | 0.1 | 🔴 |
| Multigraphs | 5% | STRidge | 0.2 | Axial decay | 0/2 | 2 | 0.9526 | 0.1 | 🔴 |
| Multigraphs | 5% | STRidge | 0.1 | Axial growth | 0/2 | 2 | 0.9540 | 0.1 | 🔴 |
| Multigraphs | 5% | SR3 | 0.01 | Radial decay | 0/2 | 2 | 0.9517 | 0.0 | 🔴 |
| Multigraphs | 5% | SR3 | 0.05 | Radial growth | 0/2 | 2 | 0.9556 | 0.0 | 🔴 |
| Multigraphs | 5% | SR3 | 0.01 | Axial decay | 0/2 | 2 | 0.9537 | 0.0 | 🔴 |
| Multigraphs | 5% | SR3 | 0.01 | Axial growth | 0/2 | 2 | 0.9545 | 0.0 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | Radial decay | 0/2 | 2 | 0.9501 | 0.1 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | Radial growth | 0/2 | 2 | 0.9547 | 0.1 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | Axial decay | 0/2 | 2 | 0.9527 | 0.2 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | Axial growth | 0/2 | 2 | 0.9544 | 0.1 | 🔴 |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol/α | Profile | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|---|
| PDE-FIND | 2% | STRidge | 0.05 | Radial decay | 0/0 | 4 | 0.0661 | 0.3 | 🟢 |
| PDE-FIND | 2% | STRidge | 0.05 | Radial growth | 0/0 | 4 | 0.0637 | 0.4 | 🟢 |
| PDE-FIND | 2% | STRidge | 0.05 | Axial decay | 0/0 | 4 | 0.0600 | 0.4 | 🟢 |
| PDE-FIND | 2% | STRidge | 0.05 | Axial growth | 0/0 | 4 | 0.0661 | 0.3 | 🟢 |
| PDE-FIND | 2% | SR3 | 0.005 | Radial decay | 0/0 | 4 | 0.5703 | 0.1 | 🟢 |
| PDE-FIND | 2% | SR3 | 0.005 | Radial growth | 0/0 | 4 | 0.5627 | 0.1 | 🟢 |
| PDE-FIND | 2% | SR3 | 0.005 | Axial decay | 0/0 | 4 | 0.5729 | 0.1 | 🟢 |
| PDE-FIND | 2% | SR3 | 0.005 | Axial growth | 0/0 | 4 | 0.5598 | 0.1 | 🟢 |
| PDE-FIND | 2% | MSTLS | 0.25 | Radial decay | 0/0 | 4 | 0.0656 | 0.6 | 🟢 |
| PDE-FIND | 2% | MSTLS | 0.25 | Radial growth | 0/0 | 4 | 0.0601 | 0.7 | 🟢 |
| PDE-FIND | 2% | MSTLS | 0.25 | Axial decay | 0/0 | 4 | 0.0587 | 0.7 | 🟢 |
| PDE-FIND | 2% | MSTLS | 0.25 | Axial growth | 0/0 | 4 | 0.0623 | 0.6 | 🟢 |
| PDE-FIND | 5% | STRidge | 0.05 | Radial decay | 0/0 | 4 | 0.4838 | 0.3 | 🟢 |
| PDE-FIND | 5% | STRidge | 0.05 | Radial growth | 0/0 | 4 | 0.3375 | 0.3 | 🟢 |
| PDE-FIND | 5% | STRidge | 0.05 | Axial decay | 0/0 | 4 | 0.3276 | 0.3 | 🟢 |
| PDE-FIND | 5% | STRidge | 0.05 | Axial growth | 0/0 | 4 | 0.3462 | 0.3 | 🟢 |
| PDE-FIND | 5% | SR3 | 0.001 | Radial decay | 0/0 | 4 | 0.5800 | 0.1 | 🟢 |
| PDE-FIND | 5% | SR3 | 0.001 | Radial growth | 0/0 | 4 | 0.5923 | 0.1 | 🟢 |
| PDE-FIND | 5% | SR3 | 0.0005 | Axial decay | 0/0 | 4 | 0.6006 | 0.1 | 🟢 |
| PDE-FIND | 5% | SR3 | 0.001 | Axial growth | 0/0 | 4 | 0.6113 | 0.1 | 🟢 |
| PDE-FIND | 5% | MSTLS | 0.25 | Radial decay | 0/0 | 4 | 0.3542 | 0.7 | 🟢 |
| PDE-FIND | 5% | MSTLS | 0.25 | Radial growth | 0/0 | 4 | 0.3297 | 0.7 | 🟢 |
| PDE-FIND | 5% | MSTLS | 0.25 | Axial decay | 0/0 | 4 | 0.3272 | 0.7 | 🟢 |
| PDE-FIND | 5% | MSTLS | 0.25 | Axial growth | 0/0 | 4 | 0.3377 | 0.6 | 🟢 |
| WSINDy | 2% | STRidge | 0.05 | Radial decay | 0/0 | 4 | 0.0661 | 0.2 | 🟢 |
| WSINDy | 2% | STRidge | 0.05 | Radial growth | 0/0 | 4 | 0.0637 | 0.2 | 🟢 |
| WSINDy | 2% | STRidge | 0.05 | Axial decay | 0/0 | 4 | 0.0600 | 0.2 | 🟢 |
| WSINDy | 2% | STRidge | 0.05 | Axial growth | 0/0 | 4 | 0.0661 | 0.2 | 🟢 |
| WSINDy | 2% | SR3 | 0.001 | Radial decay | 0/0 | 4 | 0.5156 | 0.1 | 🟢 |
| WSINDy | 2% | SR3 | 0.005 | Radial growth | 0/0 | 4 | 0.5416 | 0.1 | 🟢 |
| WSINDy | 2% | SR3 | 0.001 | Axial decay | 0/0 | 4 | 0.5173 | 0.1 | 🟢 |
| WSINDy | 2% | SR3 | 0.005 | Axial growth | 0/0 | 4 | 0.5411 | 0.1 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | Radial decay | 0/0 | 4 | 0.0656 | 0.4 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | Radial growth | 0/0 | 4 | 0.0601 | 0.4 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | Axial decay | 0/0 | 4 | 0.0587 | 0.4 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | Axial growth | 0/0 | 4 | 0.0623 | 0.4 | 🟢 |
| WSINDy | 5% | STRidge | 0.05 | Radial decay | 0/0 | 4 | 0.3495 | 0.2 | 🟢 |
| WSINDy | 5% | STRidge | 0.05 | Radial growth | 0/0 | 4 | 0.3375 | 0.2 | 🟢 |
| WSINDy | 5% | STRidge | 0.05 | Axial decay | 0/0 | 4 | 0.3276 | 0.2 | 🟢 |
| WSINDy | 5% | STRidge | 0.05 | Axial growth | 0/0 | 4 | 0.3462 | 0.2 | 🟢 |
| WSINDy | 5% | SR3 | 0.001 | Radial decay | 0/0 | 4 | 0.6438 | 0.1 | 🟢 |
| WSINDy | 5% | SR3 | 0.001 | Radial growth | 0/0 | 4 | 0.5980 | 0.1 | 🟢 |
| WSINDy | 5% | SR3 | 0.001 | Axial decay | 0/0 | 4 | 0.5977 | 0.1 | 🟢 |
| WSINDy | 5% | SR3 | 0.005 | Axial growth | 0/0 | 4 | 0.6725 | 0.1 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | Radial decay | 0/0 | 4 | 0.3542 | 0.4 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | Radial growth | 0/0 | 4 | 0.3297 | 0.4 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | Axial decay | 0/0 | 4 | 0.3272 | 0.4 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | Axial growth | 0/0 | 4 | 0.3377 | 0.4 | 🟢 |
| PolyJet | 2% | STRidge | 0.05 | Radial decay | 0/0 | 4 | 0.0661 | 0.5 | 🟢 |
| PolyJet | 2% | STRidge | 0.05 | Radial growth | 0/0 | 4 | 0.0637 | 0.5 | 🟢 |
| PolyJet | 2% | STRidge | 0.05 | Axial decay | 0/0 | 4 | 0.0600 | 0.6 | 🟢 |
| PolyJet | 2% | STRidge | 0.05 | Axial growth | 0/0 | 4 | 0.0661 | 0.6 | 🟢 |
| PolyJet | 2% | SR3 | 0.005 | Radial decay | 0/2 | 2 | 0.7684 | 0.1 | 🔴 |
| PolyJet | 2% | SR3 | 0.005 | Radial growth | 0/2 | 2 | 0.7704 | 0.1 | 🔴 |
| PolyJet | 2% | SR3 | 0.005 | Axial decay | 0/2 | 2 | 0.7693 | 0.1 | 🔴 |
| PolyJet | 2% | SR3 | 0.005 | Axial growth | 0/2 | 2 | 0.7703 | 0.1 | 🔴 |
| PolyJet | 2% | MSTLS | 0.25 | Radial decay | 0/0 | 4 | 0.0656 | 1.2 | 🟢 |
| PolyJet | 2% | MSTLS | 0.25 | Radial growth | 0/0 | 4 | 0.0601 | 1.2 | 🟢 |
| PolyJet | 2% | MSTLS | 0.25 | Axial decay | 0/0 | 4 | 0.0587 | 1.2 | 🟢 |
| PolyJet | 2% | MSTLS | 0.25 | Axial growth | 0/0 | 4 | 0.0623 | 1.2 | 🟢 |
| PolyJet | 5% | STRidge | 0.05 | Radial decay | 0/0 | 4 | 0.3495 | 0.5 | 🟢 |
| PolyJet | 5% | STRidge | 0.05 | Radial growth | 0/0 | 4 | 0.4765 | 0.5 | 🟢 |
| PolyJet | 5% | STRidge | 0.05 | Axial decay | 0/0 | 4 | 0.3276 | 0.5 | 🟢 |
| PolyJet | 5% | STRidge | 0.05 | Axial growth | 0/0 | 4 | 0.4843 | 0.5 | 🟢 |
| PolyJet | 5% | SR3 | 0.001 | Radial decay | 0/2 | 2 | 0.9520 | 0.1 | 🔴 |
| PolyJet | 5% | SR3 | 0.01 | Radial growth | 0/2 | 2 | 0.9564 | 0.1 | 🔴 |
| PolyJet | 5% | SR3 | 0.005 | Axial decay | 0/2 | 2 | 0.9545 | 0.1 | 🔴 |
| PolyJet | 5% | SR3 | 0.005 | Axial growth | 0/2 | 2 | 0.9561 | 0.1 | 🔴 |
| PolyJet | 5% | MSTLS | 0.25 | Radial decay | 0/0 | 4 | 0.3542 | 1.4 | 🟢 |
| PolyJet | 5% | MSTLS | 0.0 | Radial growth | 0/2 | 2 | 0.9547 | 1.3 | 🔴 |
| PolyJet | 5% | MSTLS | 0.25 | Axial decay | 0/0 | 4 | 0.3272 | 1.5 | 🟢 |
| PolyJet | 5% | MSTLS | 0.0 | Axial growth | 0/2 | 2 | 0.9544 | 1.3 | 🔴 |
| Multigraphs | 2% | STRidge | 0.05 | Radial decay | 0/0 | 3 | 0.0660 | 0.1 | 🟢 |
| Multigraphs | 2% | STRidge | 0.05 | Radial growth | 0/0 | 3 | 0.0637 | 0.1 | 🟢 |
| Multigraphs | 2% | STRidge | 0.05 | Axial decay | 0/0 | 3 | 0.0600 | 0.1 | 🟢 |
| Multigraphs | 2% | STRidge | 0.05 | Axial growth | 0/0 | 3 | 0.0661 | 0.1 | 🟢 |
| Multigraphs | 2% | SR3 | 0.005 | Radial decay | 0/2 | 2 | 0.7693 | 0.0 | 🔴 |
| Multigraphs | 2% | SR3 | 0.005 | Radial growth | 0/2 | 2 | 0.7715 | 0.0 | 🔴 |
| Multigraphs | 2% | SR3 | 0.005 | Axial decay | 0/2 | 2 | 0.7704 | 0.0 | 🔴 |
| Multigraphs | 2% | SR3 | 0.005 | Axial growth | 0/0 | 3 | 0.6994 | 0.0 | 🟢 |
| Multigraphs | 2% | MSTLS | 0.25 | Radial decay | 0/0 | 3 | 0.0656 | 0.1 | 🟢 |
| Multigraphs | 2% | MSTLS | 0.25 | Radial growth | 0/0 | 3 | 0.0601 | 0.1 | 🟢 |
| Multigraphs | 2% | MSTLS | 0.25 | Axial decay | 0/0 | 3 | 0.0587 | 0.1 | 🟢 |
| Multigraphs | 2% | MSTLS | 0.25 | Axial growth | 0/0 | 3 | 0.0622 | 0.1 | 🟢 |
| Multigraphs | 5% | STRidge | 0.1 | Radial decay | 0/0 | 3 | 0.3492 | 0.1 | 🟢 |
| Multigraphs | 5% | STRidge | 0.1 | Radial growth | 0/0 | 3 | 0.3372 | 0.1 | 🟢 |
| Multigraphs | 5% | STRidge | 0.05 | Axial decay | 0/0 | 3 | 0.4677 | 0.1 | 🟢 |
| Multigraphs | 5% | STRidge | 0.1 | Axial growth | 0/0 | 3 | 0.3457 | 0.1 | 🟢 |
| Multigraphs | 5% | SR3 | 0.005 | Radial decay | 0/0 | 3 | 0.7993 | 0.0 | 🟢 |
| Multigraphs | 5% | SR3 | 0.005 | Radial growth | 0/0 | 3 | 0.8070 | 0.0 | 🟢 |
| Multigraphs | 5% | SR3 | 0.005 | Axial decay | 0/0 | 3 | 0.8320 | 0.0 | 🟢 |
| Multigraphs | 5% | SR3 | 0.005 | Axial growth | 0/0 | 3 | 0.7950 | 0.0 | 🟢 |
| Multigraphs | 5% | MSTLS | 0.0 | Radial decay | 0/2 | 2 | 0.9501 | 0.1 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.25 | Radial growth | 0/0 | 3 | 0.3293 | 0.1 | 🟢 |
| Multigraphs | 5% | MSTLS | 0.25 | Axial decay | 0/0 | 3 | 0.3271 | 0.1 | 🟢 |
| Multigraphs | 5% | MSTLS | 0.25 | Axial growth | 0/0 | 3 | 0.3372 | 0.1 | 🟢 |

---
*Source: `sweep_20260630_161023_aniso_full_ac_sh.csv`*

*Run: 716.7 min (Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz), OMP_NUM_THREADS=8*
