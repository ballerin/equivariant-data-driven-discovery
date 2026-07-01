# Swift-Hohenberg equation — Anisotropic noise

The equation under consideration is

$$f_t = (r-1)f - 2\Delta f - \Delta^2 f - f^3, \qquad r = 0.3$$

which can be written in terms of $O(2)$-equivariant terms

| Operator | True coefficient |
|---|---|
| $f$ | $r-1 = -0.7$ |
| $\Delta f$ | $-2.0$ |
| $\Delta^2 f$ | $-1.0$ |
| $f^3$ | $-1.0$ |

or standard terms (PolyJet basis)

| Operator | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|---|---|---|---|---|---|
| $f$ | | 🟢 | 🟢 | 🟢 | 🟢 |
| $f^3$ | | 🟢 | 🟢 | 🟢 | 🟢 |
| $\partial_{xx} f$ | | 🟢 | 🟢 | 🟢 | 🟢 |
| $\partial_{yy} f$ | | 🟢 | 🟢 | 🟢 | 🟢 |
| $\partial_{xxxx} f$ | | 🟢 | 🟢 | 🟢 | 🟢 |
| $\partial_{xxyy} f$ | | 🟢 | 🟢 | 🟢 | 🟢 |
| $\partial_{yyyy} f$ | | 🟢 | 🟢 | 🟢 | 🟢 |

All libraries are evaluated in the common polynomial-jet basis (105 monomials,
max order 4, max degree 3).  fp/fn/E are computed against the same ground truth
for every library.  Swift-Hohenberg contains a 4th-order biharmonic term $\Delta^2 f = f_{xxxx} + 2f_{xxyy} + f_{yyyy}$ which is the main source of difficulty.  All four libraries can represent every term individually, but the equivariant basis encodes the isotropic combination as a single operator, giving it a $3\times$ sparsity advantage over coordinate bases for the biharmonic term.

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
| 2% | pointwise | spectral | Radial decay | — | 0/1 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | — | — |
| 2% | pointwise | spectral | Radial growth | — | 0/4 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | — | — |
| 2% | pointwise | spectral | Axial decay | — | 0/1 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | — | — |
| 2% | pointwise | spectral | Axial growth | — | 0/4 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | — | — |
| 2% | pointwise | polynomial | Radial decay | — | 0/1 ✗ | 1.0009 | — | 0/4 ✗ | 1.0000 | — | — |
| 2% | pointwise | polynomial | Radial growth | — | 0/1 ✗ | 1.0008 | — | 0/4 ✗ | 1.0000 | — | — |
| 2% | pointwise | polynomial | Axial decay | — | 0/1 ✗ | 1.0008 | — | 0/4 ✗ | 1.0000 | — | — |
| 2% | pointwise | polynomial | Axial growth | — | 0/1 ✗ | 1.0008 | — | 0/4 ✗ | 1.0000 | — | — |
| 2% | weak (no IBP) | spectral | Radial decay | 🟥 | 0/0 ✓ | 0.0078 | PolyJet | 0/0 ✓ | 0.0068 | — | — |
| 2% | weak (no IBP) | spectral | Radial growth | ⚪ | 0/0 ✓ | 0.0076 | PDE-FIND | 0/0 ✓ | 0.0076 | — | — |
| 2% | weak (no IBP) | spectral | Axial decay | 🟩 | 0/0 ✓ | 0.0008 | PolyJet | 0/0 ✓ | 0.0025 | — | — |
| 2% | weak (no IBP) | spectral | Axial growth | 🟥 | 0/0 ✓ | 0.0074 | PolyJet | 0/0 ✓ | 0.0072 | — | — |
| 2% | weak (no IBP) | polynomial | Radial decay | 🟥 | 0/0 ✓ | 0.0065 | PDE-FIND | 0/0 ✓ | 0.0052 | — | — |
| 2% | weak (no IBP) | polynomial | Radial growth | 🟩 | 0/0 ✓ | 0.0086 | PolyJet | 0/0 ✓ | 0.0087 | — | — |
| 2% | weak (no IBP) | polynomial | Axial decay | 🟩 | 0/0 ✓ | 0.0003 | PDE-FIND | 0/0 ✓ | 0.0019 | — | — |
| 2% | weak (no IBP) | polynomial | Axial growth | ⚪ | 0/0 ✓ | 0.0084 | PDE-FIND | 0/0 ✓ | 0.0083 | — | — |
| 2% | weak (IBP) | spectral | Radial decay | 🟥 | 0/0 ✓ | 0.0094 | PolyJet | 0/0 ✓ | 0.0055 | — | — |
| 2% | weak (IBP) | spectral | Radial growth | 🟩 | 0/0 ✓ | 0.0059 | WSINDy | 0/0 ✓ | 0.0090 | — | — |
| 2% | weak (IBP) | spectral | Axial decay | 🟥 | 0/0 ✓ | 0.0024 | PolyJet | 0/0 ✓ | 0.0020 | — | — |
| 2% | weak (IBP) | spectral | Axial growth | 🟩 | 0/0 ✓ | 0.0057 | PDE-FIND | 0/0 ✓ | 0.0086 | — | — |
| 2% | weak (IBP) | polynomial | Radial decay | 🟥 | 0/0 ✓ | 0.0094 | PDE-FIND | 0/0 ✓ | 0.0055 | — | — |
| 2% | weak (IBP) | polynomial | Radial growth | 🟩 | 0/0 ✓ | 0.0059 | WSINDy | 0/0 ✓ | 0.0090 | — | — |
| 2% | weak (IBP) | polynomial | Axial decay | 🟥 | 0/0 ✓ | 0.0024 | PolyJet | 0/0 ✓ | 0.0020 | — | — |
| 2% | weak (IBP) | polynomial | Axial growth | 🟩 | 0/0 ✓ | 0.0057 | PDE-FIND | 0/0 ✓ | 0.0086 | — | — |
| 5% | pointwise | spectral | Radial decay | — | 0/6 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | — | — |
| 5% | pointwise | spectral | Radial growth | — | 0/5 ✗ | 1.0000 | — | 0/5 ✗ | 1.0000 | — | — |
| 5% | pointwise | spectral | Axial decay | — | 0/6 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | — | — |
| 5% | pointwise | spectral | Axial growth | — | 0/6 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | — | — |
| 5% | pointwise | polynomial | Radial decay | — | 0/1 ✗ | 1.0001 | — | 0/4 ✗ | 1.0000 | — | — |
| 5% | pointwise | polynomial | Radial growth | — | 0/4 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | — | — |
| 5% | pointwise | polynomial | Axial decay | — | 0/4 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | — | — |
| 5% | pointwise | polynomial | Axial growth | — | 0/4 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | — | — |
| 5% | weak (no IBP) | spectral | Radial decay | 🟩 | 0/0 ✓ | 0.1139 | PDE-FIND | 0/0 ✓ | 0.1302 | — | — |
| 5% | weak (no IBP) | spectral | Radial growth | 🟩 | 0/0 ✓ | 0.0333 | PDE-FIND | 0/0 ✓ | 0.0392 | — | — |
| 5% | weak (no IBP) | spectral | Axial decay | 🟩 | 0/0 ✓ | 0.0782 | PDE-FIND | 0/0 ✓ | 0.0909 | — | — |
| 5% | weak (no IBP) | spectral | Axial growth | 🟩 | 0/0 ✓ | 0.0345 | PDE-FIND | 0/0 ✓ | 0.0406 | — | — |
| 5% | weak (no IBP) | polynomial | Radial decay | 🟩 | 0/0 ✓ | 0.1102 | PolyJet | 0/0 ✓ | 0.1301 | — | — |
| 5% | weak (no IBP) | polynomial | Radial growth | 🟩 | 0/0 ✓ | 0.0313 | PDE-FIND | 0/0 ✓ | 0.0369 | — | — |
| 5% | weak (no IBP) | polynomial | Axial decay | 🟩 | 0/0 ✓ | 0.0753 | PDE-FIND | 0/0 ✓ | 0.0874 | — | — |
| 5% | weak (no IBP) | polynomial | Axial growth | 🟩 | 0/0 ✓ | 0.0325 | PDE-FIND | 0/0 ✓ | 0.0382 | — | — |
| 5% | weak (IBP) | spectral | Radial decay | 🟩 | 0/0 ✓ | 0.1151 | PDE-FIND | 0/0 ✓ | 0.1294 | — | — |
| 5% | weak (IBP) | spectral | Radial growth | 🟩 | 0/0 ✓ | 0.0348 | PDE-FIND | 0/0 ✓ | 0.0380 | — | — |
| 5% | weak (IBP) | spectral | Axial decay | 🟩 | 0/0 ✓ | 0.0795 | PDE-FIND | 0/0 ✓ | 0.0900 | — | — |
| 5% | weak (IBP) | spectral | Axial growth | 🟩 | 0/0 ✓ | 0.0360 | PDE-FIND | 0/0 ✓ | 0.0394 | — | — |
| 5% | weak (IBP) | polynomial | Radial decay | 🟩 | 0/0 ✓ | 0.1151 | PolyJet | 0/0 ✓ | 0.1336 | — | — |
| 5% | weak (IBP) | polynomial | Radial growth | 🟩 | 0/0 ✓ | 0.0348 | PDE-FIND | 0/0 ✓ | 0.0380 | — | — |
| 5% | weak (IBP) | polynomial | Axial decay | 🟩 | 0/0 ✓ | 0.0795 | PDE-FIND | 0/0 ✓ | 0.0900 | — | — |
| 5% | weak (IBP) | polynomial | Axial growth | 🟩 | 0/0 ✓ | 0.0360 | PDE-FIND | 0/0 ✓ | 0.0394 | — | — |

## Pointwise formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol/α | Profile | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|---|
| PDE-FIND | 2% | STRidge | 0.05 | Radial decay | 0/5 | 2 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.01 | Radial growth | 0/6 | 1 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.01 | Axial decay | 0/6 | 1 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.01 | Axial growth | 0/6 | 1 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.1 | Radial decay | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.1 | Radial growth | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.1 | Axial decay | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.1 | Axial growth | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | Radial decay | 1/6 | 2 | 0.9999 | 1.6 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | Radial growth | 2/6 | 3 | 1.0000 | 1.7 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | Axial decay | 1/6 | 2 | 0.9999 | 1.8 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | Axial growth | 2/6 | 3 | 1.0000 | 1.5 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.05 | Radial decay | 0/5 | 2 | 1.0000 | 0.7 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.05 | Radial growth | 0/5 | 2 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.05 | Axial decay | 0/5 | 2 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.05 | Axial growth | 0/5 | 2 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.05 | Radial decay | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.05 | Radial growth | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.05 | Axial decay | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.05 | Axial growth | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | Radial decay | 0/5 | 2 | 1.0000 | 1.7 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | Radial growth | 0/5 | 2 | 1.0000 | 1.6 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | Axial decay | 0/5 | 2 | 1.0000 | 1.7 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | Axial growth | 0/5 | 2 | 1.0000 | 1.7 | 🔴 |
| WSINDy | 2% | STRidge | 0.01 | Radial decay | 0/6 | 1 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 2% | STRidge | 0.01 | Radial growth | 0/6 | 1 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 2% | STRidge | 0.01 | Axial decay | 0/6 | 1 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 2% | STRidge | 0.01 | Axial growth | 0/6 | 1 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 2% | SR3 | 0.01 | Radial decay | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 2% | SR3 | 0.01 | Radial growth | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 2% | SR3 | 0.01 | Axial decay | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 2% | SR3 | 0.01 | Axial growth | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | Radial decay | 12/3 | 8 | 1.0000 | 0.9 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | Radial growth | 12/3 | 8 | 1.0000 | 1.0 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | Axial decay | 12/3 | 8 | 1.0000 | 0.9 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | Axial growth | 12/3 | 8 | 1.0000 | 0.9 | 🔴 |
| WSINDy | 5% | STRidge | 0.05 | Radial decay | 0/5 | 2 | 1.0000 | 0.4 | 🔴 |
| WSINDy | 5% | STRidge | 0.05 | Radial growth | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 5% | STRidge | 0.05 | Axial decay | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 5% | STRidge | 0.05 | Axial growth | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 5% | SR3 | 0.01 | Radial decay | 0/5 | 2 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | SR3 | 0.01 | Radial growth | 0/5 | 2 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | SR3 | 0.1 | Axial decay | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | SR3 | 0.01 | Axial growth | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.5 | Radial decay | 12/1 | 10 | 1.0000 | 1.0 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | Radial growth | 12/1 | 10 | 1.0000 | 1.0 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | Axial decay | 12/1 | 10 | 1.0000 | 1.0 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | Axial growth | 12/1 | 10 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 2% | STRidge | 0.01 | Radial decay | 0/6 | 1 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 2% | STRidge | 0.01 | Radial growth | 0/6 | 1 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 2% | STRidge | 0.01 | Axial decay | 0/6 | 1 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 2% | STRidge | 0.01 | Axial growth | 0/6 | 1 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 2% | SR3 | 0.05 | Radial decay | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 2% | SR3 | 0.05 | Radial growth | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 2% | SR3 | 0.01 | Axial decay | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 2% | SR3 | 0.05 | Axial growth | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | Radial decay | 3/3 | 7 | 1.0000 | 3.3 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | Radial growth | 1/5 | 3 | 1.0000 | 3.2 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | Axial decay | 1/5 | 3 | 1.0000 | 3.3 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | Axial growth | 1/5 | 3 | 1.0000 | 3.0 | 🔴 |
| PolyJet | 5% | STRidge | 0.01 | Radial decay | 0/6 | 1 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 5% | STRidge | 0.01 | Radial growth | 0/6 | 1 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 5% | STRidge | 0.01 | Axial decay | 0/6 | 1 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 5% | STRidge | 0.01 | Axial growth | 0/6 | 1 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 5% | SR3 | 0.01 | Radial decay | 0/4 | 3 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 5% | SR3 | 0.005 | Radial growth | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 5% | SR3 | 0.005 | Axial decay | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 5% | SR3 | 0.005 | Axial growth | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | Radial decay | 3/3 | 7 | 1.0000 | 3.5 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | Radial growth | 3/3 | 7 | 1.0000 | 3.2 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | Axial decay | 3/3 | 7 | 1.0000 | 3.2 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | Axial growth | 3/3 | 7 | 1.0000 | 3.0 | 🔴 |
| Multigraphs | 2% | STRidge | 0.005 | Radial decay | 0/4 | 2 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 2% | STRidge | 0.005 | Radial growth | 0/6 | 1 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 2% | STRidge | 0.005 | Axial decay | 0/4 | 2 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 2% | STRidge | 0.005 | Axial growth | 0/6 | 1 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | Radial decay | 0/1 | 3 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.05 | Radial growth | 0/4 | 2 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | Axial decay | 0/1 | 3 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.05 | Axial growth | 0/4 | 2 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | Radial decay | 2/5 | 3 | 1.0000 | 0.4 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | Radial growth | 2/5 | 3 | 1.0000 | 0.4 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | Axial decay | 2/5 | 3 | 1.0000 | 0.4 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | Axial growth | 2/5 | 3 | 1.0000 | 0.4 | 🔴 |
| Multigraphs | 5% | STRidge | 0.005 | Radial decay | 0/6 | 1 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | STRidge | 0.005 | Radial growth | 0/6 | 1 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | STRidge | 0.005 | Axial decay | 0/6 | 1 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | STRidge | 0.01 | Axial growth | 0/6 | 1 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | SR3 | 0.1 | Radial decay | 0/6 | 1 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | SR3 | 0.01 | Radial growth | 0/5 | 2 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | SR3 | 0.1 | Axial decay | 0/6 | 1 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | SR3 | 0.1 | Axial growth | 0/6 | 1 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | Radial decay | 6/5 | 4 | 1.0000 | 0.4 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | Radial growth | 6/5 | 4 | 1.0000 | 0.4 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | Axial decay | 6/5 | 4 | 1.0000 | 0.4 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | Axial growth | 6/5 | 4 | 1.0000 | 0.4 | 🔴 |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol/α | Profile | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|---|
| PDE-FIND | 2% | STRidge | 0.1 | Radial decay | 0/7 | 0 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.1 | Radial growth | 0/7 | 0 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.1 | Axial decay | 0/7 | 0 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.1 | Axial growth | 0/7 | 0 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.01 | Radial decay | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.01 | Radial growth | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.01 | Axial decay | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.01 | Axial growth | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | Radial decay | 1/5 | 3 | 1.0002 | 1.7 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.25 | Radial growth | 1/5 | 3 | 1.0002 | 1.7 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | Axial decay | 1/5 | 3 | 1.0002 | 1.6 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.25 | Axial growth | 1/5 | 3 | 1.0002 | 1.6 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.1 | Radial decay | 0/7 | 0 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.1 | Radial growth | 0/7 | 0 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.1 | Axial decay | 0/7 | 0 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.1 | Axial growth | 0/7 | 0 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.005 | Radial decay | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.005 | Radial growth | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.005 | Axial decay | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.005 | Axial growth | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | Radial decay | 1/5 | 3 | 1.0001 | 1.8 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | Radial growth | 1/5 | 3 | 1.0001 | 1.6 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | Axial decay | 1/5 | 3 | 1.0001 | 1.7 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | Axial growth | 1/5 | 3 | 1.0001 | 1.5 | 🔴 |
| WSINDy | 2% | STRidge | 0.01 | Radial decay | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 2% | STRidge | 0.01 | Radial growth | 0/5 | 2 | 1.0000 | 0.4 | 🔴 |
| WSINDy | 2% | STRidge | 0.01 | Axial decay | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 2% | STRidge | 0.01 | Axial growth | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 2% | SR3 | 0.01 | Radial decay | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 2% | SR3 | 0.01 | Radial growth | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 2% | SR3 | 0.01 | Axial decay | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 2% | SR3 | 0.01 | Axial growth | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | Radial decay | 6/5 | 4 | 0.9992 | 1.0 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | Radial growth | 6/5 | 4 | 0.9992 | 1.0 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | Axial decay | 6/5 | 4 | 0.9992 | 1.0 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | Axial growth | 6/5 | 4 | 0.9992 | 1.0 | 🔴 |
| WSINDy | 5% | STRidge | 0.1 | Radial decay | 0/7 | 0 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 5% | STRidge | 0.01 | Radial growth | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 5% | STRidge | 0.1 | Axial decay | 0/7 | 0 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 5% | STRidge | 0.01 | Axial growth | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 5% | SR3 | 0.01 | Radial decay | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | SR3 | 0.01 | Radial growth | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | SR3 | 0.01 | Axial decay | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | SR3 | 0.01 | Axial growth | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | Radial decay | 2/5 | 3 | 1.0001 | 1.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | Radial growth | 6/5 | 4 | 0.9997 | 1.0 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | Axial decay | 2/5 | 3 | 1.0002 | 1.0 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | Axial growth | 6/5 | 4 | 0.9997 | 1.0 | 🔴 |
| PolyJet | 2% | STRidge | 0.5 | Radial decay | 0/7 | 0 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 2% | STRidge | 0.1 | Radial growth | 0/7 | 0 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 2% | STRidge | 0.1 | Axial decay | 0/7 | 0 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 2% | STRidge | 0.5 | Axial growth | 0/7 | 0 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 2% | SR3 | 0.01 | Radial decay | 0/5 | 2 | 1.0001 | 0.3 | 🔴 |
| PolyJet | 2% | SR3 | 0.01 | Radial growth | 0/5 | 2 | 1.0001 | 0.3 | 🔴 |
| PolyJet | 2% | SR3 | 0.01 | Axial decay | 0/5 | 2 | 1.0001 | 0.3 | 🔴 |
| PolyJet | 2% | SR3 | 0.01 | Axial growth | 0/5 | 2 | 1.0001 | 0.3 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | Radial decay | 9/5 | 11 | 1.0003 | 4.0 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | Radial growth | 9/5 | 11 | 1.0003 | 4.0 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | Axial decay | 9/5 | 11 | 1.0003 | 3.6 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | Axial growth | 9/5 | 11 | 1.0003 | 3.6 | 🔴 |
| PolyJet | 5% | STRidge | 0.1 | Radial decay | 0/7 | 0 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 5% | STRidge | 0.5 | Radial growth | 0/7 | 0 | 1.0000 | 1.0 | 🔴 |
| PolyJet | 5% | STRidge | 0.2 | Axial decay | 0/7 | 0 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 5% | STRidge | 0.2 | Axial growth | 0/7 | 0 | 1.0000 | 0.8 | 🔴 |
| PolyJet | 5% | SR3 | 0.01 | Radial decay | 0/5 | 2 | 1.0001 | 0.3 | 🔴 |
| PolyJet | 5% | SR3 | 0.01 | Radial growth | 0/5 | 2 | 1.0001 | 0.3 | 🔴 |
| PolyJet | 5% | SR3 | 0.01 | Axial decay | 0/5 | 2 | 1.0001 | 0.3 | 🔴 |
| PolyJet | 5% | SR3 | 0.005 | Axial growth | 0/5 | 2 | 1.0001 | 0.2 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | Radial decay | 4/5 | 6 | 1.0001 | 3.5 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | Radial growth | 5/5 | 7 | 1.0001 | 3.5 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | Axial decay | 4/5 | 6 | 1.0001 | 3.4 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | Axial growth | 5/5 | 7 | 1.0001 | 3.4 | 🔴 |
| Multigraphs | 2% | STRidge | 0.005 | Radial decay | 0/4 | 2 | 1.0001 | 0.1 | 🔴 |
| Multigraphs | 2% | STRidge | 0.005 | Radial growth | 0/4 | 2 | 1.0001 | 0.1 | 🔴 |
| Multigraphs | 2% | STRidge | 0.005 | Axial decay | 0/4 | 2 | 1.0001 | 0.1 | 🔴 |
| Multigraphs | 2% | STRidge | 0.005 | Axial growth | 0/4 | 2 | 1.0001 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | Radial decay | 0/1 | 3 | 1.0009 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | Radial growth | 0/1 | 3 | 1.0008 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | Axial decay | 0/1 | 3 | 1.0008 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | Axial growth | 0/1 | 3 | 1.0008 | 0.1 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | Radial decay | 19/4 | 11 | 1.0001 | 0.4 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | Radial growth | 19/3 | 12 | 1.0000 | 0.4 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | Axial decay | 19/4 | 11 | 1.0001 | 0.4 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | Axial growth | 19/4 | 11 | 1.0000 | 0.5 | 🔴 |
| Multigraphs | 5% | STRidge | 0.005 | Radial decay | 0/4 | 2 | 1.0001 | 0.1 | 🔴 |
| Multigraphs | 5% | STRidge | 0.005 | Radial growth | 0/4 | 2 | 1.0001 | 0.1 | 🔴 |
| Multigraphs | 5% | STRidge | 0.005 | Axial decay | 0/4 | 2 | 1.0001 | 0.1 | 🔴 |
| Multigraphs | 5% | STRidge | 0.005 | Axial growth | 0/6 | 1 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | SR3 | 0.01 | Radial decay | 0/1 | 3 | 1.0001 | 0.1 | 🔴 |
| Multigraphs | 5% | SR3 | 0.1 | Radial growth | 0/4 | 2 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | SR3 | 0.1 | Axial decay | 0/4 | 2 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | SR3 | 0.1 | Axial growth | 0/4 | 2 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | Radial decay | 13/5 | 8 | 1.0001 | 0.5 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | Radial growth | 13/5 | 8 | 1.0001 | 0.4 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | Axial decay | 13/5 | 8 | 1.0001 | 0.4 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | Axial growth | 13/5 | 8 | 1.0000 | 0.4 | 🔴 |

---
*Source: `sweep_20260630_161023_aniso_full_ac_sh.csv`*

*Run: 716.7 min (Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz), OMP_NUM_THREADS=8*
