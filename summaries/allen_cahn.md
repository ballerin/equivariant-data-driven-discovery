# Allen-Cahn equation

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
max order 4, max degree 3).  fp/fn/E are computed against the same 4-term
ground truth for every library.  Allen-Cahn contains only pure monomials ($f$, $f^3$) and pure second derivatives ($f_{xx}$, $f_{yy}$).  All four libraries can represent every term, so the structural advantage of the equivariant basis is less pronounced here; differences come from numerical conditioning, not representability.

Sweep configuration: all 4 libraries, 3 algorithms (STRidge, SR3, MSTLS),
3 methods (pointwise, weak-direct, weak-IBP), 2 backends (spectral, polynomial),
3 noise levels (clean, 2%, 5%).

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


## Summary: equivariant vs best non-equivariant

| Noise | Method | Backend | Eq best | Multigraphs fp/fn | Multigraphs E | Non-eq best | Non-eq fp/fn | Non-eq E | MWSINDy fp/fn | MWSINDy E |
|---|---|---|---|---|---|---|---|---|---|---|
| Clean | pointwise | spectral | 🟩 | 0/0 ✓ | 0.0080 | PDE-FIND | 0/0 ✓ | 0.0080 | 0/0 ✓ | 0.0080 |
| Clean | pointwise | polynomial | 🟩 | 0/0 ✓ | 0.0114 | PolyJet | 0/0 ✓ | 0.0115 | 0/0 ✓ | 0.0114 |
| Clean | weak (no IBP) | spectral | 🟩 | 0/0 ✓ | 0.0048 | WSINDy | 0/0 ✓ | 0.0048 | 0/0 ✓ | 0.0048 |
| Clean | weak (no IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0048 | WSINDy | 0/0 ✓ | 0.0048 | 0/0 ✓ | 0.0048 |
| Clean | weak (IBP) | spectral | 🟩 | 0/0 ✓ | 0.0048 | WSINDy | 0/0 ✓ | 0.0048 | 0/0 ✓ | 0.0048 |
| Clean | weak (IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0048 | WSINDy | 0/0 ✓ | 0.0048 | 0/0 ✓ | 0.0048 |
| 2% | pointwise | spectral | ⚪ | 0/2 | 0.7637 | PDE-FIND | 0/2 | 0.7637 | 0/0 ✓ | 0.6763 |
| 2% | pointwise | polynomial | 🟩 | 0/0 ✓ | 0.0540 | PolyJet | 0/0 ✓ | 0.0543 | 0/0 ✓ | 0.0540 |
| 2% | weak (no IBP) | spectral | 🟩 | 0/0 ✓ | 0.0045 | WSINDy | 0/0 ✓ | 0.0046 | 0/0 ✓ | 0.0052 |
| 2% | weak (no IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0045 | PolyJet | 0/0 ✓ | 0.0045 | 0/0 ✓ | 0.0052 |
| 2% | weak (IBP) | spectral | 🟩 | 0/0 ✓ | 0.0045 | PDE-FIND | 0/0 ✓ | 0.0046 | 0/0 ✓ | 0.0052 |
| 2% | weak (IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0045 | PDE-FIND | 0/0 ✓ | 0.0046 | 0/0 ✓ | 0.0052 |
| 5% | pointwise | spectral | 🔴 | 0/2 | 0.9539 | WSINDy | 0/1 | 0.9323 | 0/0 ✓ | 0.9211 |
| 5% | pointwise | polynomial | 🟥 | 0/0 ✓ | 0.3523 | PDE-FIND | 0/0 ✓ | 0.3421 | 0/0 ✓ | 0.3421 |
| 5% | weak (no IBP) | spectral | 🟩 | 0/0 ✓ | 0.0049 | WSINDy | 0/0 ✓ | 0.0054 | 0/0 ✓ | 0.0078 |
| 5% | weak (no IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0048 | PolyJet | 0/0 ✓ | 0.0053 | 0/0 ✓ | 0.0077 |
| 5% | weak (IBP) | spectral | 🟩 | 0/0 ✓ | 0.0049 | WSINDy | 0/0 ✓ | 0.0054 | 0/0 ✓ | 0.0078 |
| 5% | weak (IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0049 | WSINDy | 0/0 ✓ | 0.0054 | 0/0 ✓ | 0.0078 |

## Pointwise formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol/α | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.05 | 0/0 | 4 | 0.0081 | 0.3 | 🟢 |
| PDE-FIND | Clean | SR3 | 0.005 | 0/0 | 4 | 0.3042 | 0.1 | 🟢 |
| PDE-FIND | Clean | MSTLS | 0.25 | 0/0 | 4 | 0.0080 | 0.6 | 🟢 |
| PDE-FIND | 2% | STRidge | 0.2 | 0/2 | 2 | 0.7637 | 0.3 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.05 | 0/2 | 2 | 0.7775 | 0.1 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 0/2 | 2 | 0.7637 | 0.7 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.2 | 0/2 | 2 | 0.9539 | 0.3 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.01 | 0/2 | 2 | 0.9555 | 0.1 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | 0/2 | 2 | 0.9539 | 0.7 | 🔴 |
| WSINDy | Clean | STRidge | 0.05 | 0/0 | 4 | 0.0081 | 0.2 | 🟢 |
| WSINDy | Clean | SR3 | 0.005 | 0/0 | 4 | 0.2678 | 0.1 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 0/0 | 4 | 0.0080 | 0.3 | 🟢 |
| WSINDy | 2% | STRidge | 0.2 | 0/2 | 2 | 0.7637 | 0.2 | 🔴 |
| WSINDy | 2% | SR3 | 0.05 | 0/2 | 2 | 0.7737 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | 0/2 | 2 | 0.7637 | 0.4 | 🔴 |
| WSINDy | 5% | STRidge | 0.2 | 0/2 | 2 | 0.9539 | 0.2 | 🔴 |
| WSINDy | 5% | SR3 | 0.005 | 0/1 | 3 | 0.9323 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | 0/2 | 2 | 0.9539 | 0.4 | 🔴 |
| PolyJet | Clean | STRidge | 0.05 | 0/0 | 4 | 0.0081 | 0.6 | 🟢 |
| PolyJet | Clean | SR3 | 0.005 | 0/2 | 2 | 0.3666 | 0.2 | 🔴 |
| PolyJet | Clean | MSTLS | 0.25 | 0/0 | 4 | 0.0080 | 1.4 | 🟢 |
| PolyJet | 2% | STRidge | 0.2 | 0/2 | 2 | 0.7637 | 0.7 | 🔴 |
| PolyJet | 2% | SR3 | 0.005 | 0/2 | 2 | 0.7681 | 0.2 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | 0/2 | 2 | 0.7637 | 1.7 | 🔴 |
| PolyJet | 5% | STRidge | 0.2 | 0/2 | 2 | 0.9539 | 0.6 | 🔴 |
| PolyJet | 5% | SR3 | 0.005 | 0/1 | 3 | 0.9377 | 0.2 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | 0/2 | 2 | 0.9539 | 1.8 | 🔴 |
| Multigraphs | Clean | STRidge | 0.05 | 0/0 | 3 | 0.0080 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.005 | 0/2 | 2 | 0.3723 | 0.0 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.25 | 0/0 | 3 | 0.0080 | 0.1 | 🟢 |
| Multigraphs | 2% | STRidge | 0.2 | 0/2 | 2 | 0.7637 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | 0/2 | 2 | 0.7715 | 0.0 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 0/2 | 2 | 0.7637 | 0.1 | 🔴 |
| Multigraphs | 5% | STRidge | 0.1 | 0/2 | 2 | 0.9539 | 0.1 | 🔴 |
| Multigraphs | 5% | SR3 | 0.01 | 0/2 | 2 | 0.9545 | 0.0 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 0/2 | 2 | 0.9539 | 0.1 | 🔴 |
| MWSINDy | Clean | STRidge | 0.01 | 0/0 | 6 | 0.0080 | 0.2 | 🟢 |
| MWSINDy | Clean | SR3 | 0.005 | 0/0 | 4 | 0.0080 | 0.1 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.25 | 0/0 | 4 | 0.0080 | 0.4 | 🟢 |
| MWSINDy | 2% | STRidge | 0.05 | 0/0 | 6 | 0.6763 | 0.2 | 🟢 |
| MWSINDy | 2% | SR3 | 0.005 | 0/0 | 6 | 0.6763 | 0.1 | 🟢 |
| MWSINDy | 2% | MSTLS | 0.5 | 0/0 | 6 | 0.6763 | 0.4 | 🟢 |
| MWSINDy | 5% | STRidge | 0.1 | 0/0 | 6 | 0.9211 | 0.2 | 🟢 |
| MWSINDy | 5% | SR3 | 0.005 | 0/2 | 3 | 0.9539 | 0.1 | 🔴 |
| MWSINDy | 5% | MSTLS | 0.0 | 0/2 | 2 | 0.9539 | 0.4 | 🔴 |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol/α | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.05 | 0/0 | 4 | 0.0116 | 0.3 | 🟢 |
| PDE-FIND | Clean | SR3 | 0.001 | 0/1 | 3 | 0.3733 | 0.1 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.25 | 0/0 | 4 | 0.0115 | 0.6 | 🟢 |
| PDE-FIND | 2% | STRidge | 0.05 | 0/0 | 4 | 0.0591 | 0.3 | 🟢 |
| PDE-FIND | 2% | SR3 | 0.005 | 0/0 | 4 | 0.5559 | 0.1 | 🟢 |
| PDE-FIND | 2% | MSTLS | 0.25 | 0/0 | 4 | 0.0543 | 0.6 | 🟢 |
| PDE-FIND | 5% | STRidge | 0.05 | 0/0 | 4 | 0.3522 | 0.3 | 🟢 |
| PDE-FIND | 5% | SR3 | 0.001 | 0/0 | 4 | 0.6270 | 0.1 | 🟢 |
| PDE-FIND | 5% | MSTLS | 0.25 | 0/0 | 4 | 0.3421 | 0.6 | 🟢 |
| WSINDy | Clean | STRidge | 0.05 | 0/0 | 4 | 0.0116 | 0.2 | 🟢 |
| WSINDy | Clean | SR3 | 0.005 | 0/2 | 2 | 0.3797 | 0.1 | 🔴 |
| WSINDy | Clean | MSTLS | 0.25 | 0/0 | 4 | 0.0115 | 0.3 | 🟢 |
| WSINDy | 2% | STRidge | 0.05 | 0/0 | 4 | 0.0591 | 0.2 | 🟢 |
| WSINDy | 2% | SR3 | 0.005 | 0/0 | 4 | 0.5243 | 0.1 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | 0/0 | 4 | 0.0543 | 0.4 | 🟢 |
| WSINDy | 5% | STRidge | 0.05 | 0/0 | 4 | 0.4933 | 0.2 | 🟢 |
| WSINDy | 5% | SR3 | 0.001 | 0/0 | 4 | 0.6252 | 0.1 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | 0/0 | 4 | 0.3421 | 0.4 | 🟢 |
| PolyJet | Clean | STRidge | 0.05 | 0/0 | 4 | 0.0116 | 0.6 | 🟢 |
| PolyJet | Clean | SR3 | 0.005 | 0/2 | 2 | 0.3677 | 0.2 | 🔴 |
| PolyJet | Clean | MSTLS | 0.25 | 0/0 | 4 | 0.0115 | 1.4 | 🟢 |
| PolyJet | 2% | STRidge | 0.05 | 0/0 | 4 | 0.2432 | 0.6 | 🟢 |
| PolyJet | 2% | SR3 | 0.005 | 0/2 | 2 | 0.7710 | 0.2 | 🔴 |
| PolyJet | 2% | MSTLS | 0.25 | 0/0 | 4 | 0.0543 | 1.6 | 🟢 |
| PolyJet | 5% | STRidge | 0.05 | 0/0 | 4 | 0.3522 | 0.7 | 🟢 |
| PolyJet | 5% | SR3 | 0.005 | 0/2 | 2 | 0.9561 | 0.2 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | 0/2 | 2 | 0.9539 | 1.7 | 🔴 |
| Multigraphs | Clean | STRidge | 0.05 | 0/0 | 3 | 0.0114 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.005 | 0/2 | 2 | 0.3730 | 0.0 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.25 | 0/0 | 3 | 0.0114 | 0.1 | 🟢 |
| Multigraphs | 2% | STRidge | 0.05 | 0/0 | 3 | 0.0591 | 0.1 | 🟢 |
| Multigraphs | 2% | SR3 | 0.005 | 0/2 | 2 | 0.7720 | 0.0 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.25 | 0/0 | 3 | 0.0540 | 0.1 | 🟢 |
| Multigraphs | 5% | STRidge | 0.1 | 0/0 | 3 | 0.3523 | 0.1 | 🟢 |
| Multigraphs | 5% | SR3 | 0.005 | 0/0 | 3 | 0.8424 | 0.0 | 🟢 |
| Multigraphs | 5% | MSTLS | 0.0 | 0/2 | 2 | 0.9539 | 0.2 | 🔴 |
| MWSINDy | Clean | STRidge | 0.05 | 0/0 | 4 | 0.0114 | 0.2 | 🟢 |
| MWSINDy | Clean | SR3 | 0.005 | 0/2 | 2 | 0.3395 | 0.1 | 🔴 |
| MWSINDy | Clean | MSTLS | 0.25 | 0/0 | 4 | 0.0114 | 0.4 | 🟢 |
| MWSINDy | 2% | STRidge | 0.01 | 0/0 | 8 | 0.0540 | 0.2 | 🟢 |
| MWSINDy | 2% | SR3 | 0.005 | 0/0 | 4 | 0.0540 | 0.1 | 🟢 |
| MWSINDy | 2% | MSTLS | 0.25 | 0/0 | 4 | 0.0540 | 0.4 | 🟢 |
| MWSINDy | 5% | STRidge | 0.01 | 0/0 | 9 | 0.3421 | 0.2 | 🟢 |
| MWSINDy | 5% | SR3 | 0.001 | 0/0 | 4 | 0.3421 | 0.1 | 🟢 |
| MWSINDy | 5% | MSTLS | 0.25 | 0/0 | 4 | 0.3421 | 0.4 | 🟢 |

## Weak formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol/α | m | IBP | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.05 | 12 | no | 0/0 | 4 | 0.0048 | 0.1 | 🟢 |
| PDE-FIND | Clean | SR3 | 7.0 | 8 | no | 0/1 | 3 | 0.6983 | 0.1 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0048 | 0.3 | 🟢 |
| PDE-FIND | Clean | STRidge | 0.05 | 12 | yes | 0/0 | 4 | 0.0048 | 0.2 | 🟢 |
| PDE-FIND | Clean | SR3 | 7.0 | 8 | yes | 0/1 | 3 | 0.6994 | 0.1 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0048 | 0.3 | 🟢 |
| PDE-FIND | 2% | STRidge | 0.1 | 12 | no | 0/0 | 4 | 0.0046 | 0.2 | 🟢 |
| PDE-FIND | 2% | SR3 | 0.01 | 4 | no | 0/0 | 4 | 0.6350 | 0.1 | 🟢 |
| PDE-FIND | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0052 | 0.3 | 🟢 |
| PDE-FIND | 2% | STRidge | 0.1 | 12 | yes | 0/0 | 4 | 0.0046 | 0.1 | 🟢 |
| PDE-FIND | 2% | SR3 | 0.5 | 8 | yes | 0/2 | 2 | 0.7984 | 0.1 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0052 | 0.3 | 🟢 |
| PDE-FIND | 5% | STRidge | 0.2 | 12 | no | 0/0 | 4 | 0.0054 | 0.2 | 🟢 |
| PDE-FIND | 5% | SR3 | 0.1 | 8 | no | 0/2 | 2 | 0.9548 | 0.1 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0078 | 0.3 | 🟢 |
| PDE-FIND | 5% | STRidge | 0.2 | 12 | yes | 0/0 | 4 | 0.0054 | 0.2 | 🟢 |
| PDE-FIND | 5% | SR3 | 0.1 | 8 | yes | 0/2 | 2 | 0.9550 | 0.1 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0078 | 0.3 | 🟢 |
| WSINDy | Clean | STRidge | 0.05 | 12 | no | 0/0 | 4 | 0.0048 | 0.1 | 🟢 |
| WSINDy | Clean | SR3 | 3.0 | 4 | no | 0/0 | 4 | 0.9618 | 0.1 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0048 | 0.2 | 🟢 |
| WSINDy | Clean | STRidge | 0.05 | 12 | yes | 0/0 | 4 | 0.0048 | 0.1 | 🟢 |
| WSINDy | Clean | SR3 | 0.5 | 4 | yes | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| WSINDy | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0048 | 0.2 | 🟢 |
| WSINDy | 2% | STRidge | 0.05 | 12 | no | 0/0 | 4 | 0.0046 | 0.1 | 🟢 |
| WSINDy | 2% | SR3 | 4.0 | 4 | no | 0/1 | 3 | 0.9937 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0052 | 0.2 | 🟢 |
| WSINDy | 2% | STRidge | 0.05 | 12 | yes | 0/0 | 4 | 0.0046 | 0.1 | 🟢 |
| WSINDy | 2% | SR3 | 0.5 | 4 | yes | 0/2 | 2 | 1.0084 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0052 | 0.2 | 🟢 |
| WSINDy | 5% | STRidge | 0.1 | 12 | no | 0/0 | 4 | 0.0054 | 0.1 | 🟢 |
| WSINDy | 5% | SR3 | 0.01 | 4 | no | 0/0 | 4 | 0.5566 | 0.1 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0078 | 0.2 | 🟢 |
| WSINDy | 5% | STRidge | 0.1 | 12 | yes | 0/0 | 4 | 0.0054 | 0.1 | 🟢 |
| WSINDy | 5% | SR3 | 0.1 | 4 | yes | 0/2 | 2 | 0.9556 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0078 | 0.2 | 🟢 |
| PolyJet | Clean | STRidge | 0.1 | 12 | no | 0/0 | 4 | 0.0048 | 0.3 | 🟢 |
| PolyJet | Clean | SR3 | 3.0 | 8 | no | 0/0 | 4 | 0.7218 | 0.1 | 🟢 |
| PolyJet | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0048 | 0.7 | 🟢 |
| PolyJet | Clean | STRidge | 0.1 | 12 | yes | 0/0 | 4 | 0.0048 | 0.3 | 🟢 |
| PolyJet | Clean | SR3 | 3.0 | 8 | yes | 0/0 | 4 | 0.7226 | 0.1 | 🟢 |
| PolyJet | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0048 | 0.7 | 🟢 |
| PolyJet | 2% | STRidge | 0.1 | 12 | no | 0/0 | 4 | 0.0046 | 0.4 | 🟢 |
| PolyJet | 2% | SR3 | 0.05 | 4 | no | 0/2 | 2 | 0.7670 | 0.2 | 🔴 |
| PolyJet | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0052 | 0.7 | 🟢 |
| PolyJet | 2% | STRidge | 0.1 | 12 | yes | 0/0 | 4 | 0.0046 | 0.4 | 🟢 |
| PolyJet | 2% | SR3 | 0.2 | 8 | yes | 0/2 | 2 | 0.8108 | 0.2 | 🔴 |
| PolyJet | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0052 | 0.7 | 🟢 |
| PolyJet | 5% | STRidge | 1.0 | 12 | no | 0/0 | 4 | 0.0054 | 0.3 | 🟢 |
| PolyJet | 5% | SR3 | 0.05 | 4 | no | 0/2 | 2 | 0.9532 | 0.2 | 🔴 |
| PolyJet | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0078 | 0.7 | 🟢 |
| PolyJet | 5% | STRidge | 1.0 | 12 | yes | 0/0 | 4 | 0.0054 | 0.3 | 🟢 |
| PolyJet | 5% | SR3 | 0.1 | 8 | yes | 0/2 | 2 | 0.9584 | 0.1 | 🔴 |
| PolyJet | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0078 | 0.7 | 🟢 |
| Multigraphs | Clean | STRidge | 0.2 | 12 | no | 0/0 | 3 | 0.0048 | 0.0 | 🟢 |
| Multigraphs | Clean | SR3 | 0.05 | 4 | no | 0/2 | 2 | 0.3708 | 0.0 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 3 | 0.0048 | 0.1 | 🟢 |
| Multigraphs | Clean | STRidge | 0.2 | 12 | yes | 0/0 | 3 | 0.0048 | 0.0 | 🟢 |
| Multigraphs | Clean | SR3 | 0.5 | 4 | yes | 0/4 | 0 | 1.0000 | 0.0 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 3 | 0.0048 | 0.1 | 🟢 |
| Multigraphs | 2% | STRidge | 0.1 | 12 | no | 0/0 | 3 | 0.0045 | 0.0 | 🟢 |
| Multigraphs | 2% | SR3 | 0.05 | 4 | no | 0/2 | 2 | 0.7661 | 0.0 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 3 | 0.0052 | 0.1 | 🟢 |
| Multigraphs | 2% | STRidge | 0.1 | 12 | yes | 0/0 | 3 | 0.0045 | 0.1 | 🟢 |
| Multigraphs | 2% | SR3 | 0.5 | 4 | yes | 0/4 | 0 | 1.0000 | 0.0 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 3 | 0.0052 | 0.1 | 🟢 |
| Multigraphs | 5% | STRidge | 1.5 | 12 | no | 0/0 | 3 | 0.0049 | 0.0 | 🟢 |
| Multigraphs | 5% | SR3 | 0.05 | 4 | no | 0/2 | 2 | 0.9536 | 0.0 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 3 | 0.0078 | 0.1 | 🟢 |
| Multigraphs | 5% | STRidge | 1.5 | 12 | yes | 0/0 | 3 | 0.0049 | 0.0 | 🟢 |
| Multigraphs | 5% | SR3 | 0.05 | 4 | yes | 0/2 | 2 | 0.9537 | 0.0 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 3 | 0.0078 | 0.1 | 🟢 |
| MWSINDy | Clean | STRidge | 0.05 | 12 | no | 0/0 | 4 | 0.0048 | 0.1 | 🟢 |
| MWSINDy | Clean | SR3 | 3.0 | 4 | no | 0/0 | 4 | 0.0055 | 0.1 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0048 | 0.2 | 🟢 |
| MWSINDy | Clean | STRidge | 0.05 | 12 | yes | 0/0 | 4 | 0.0048 | 0.1 | 🟢 |
| MWSINDy | Clean | SR3 | 3.0 | 4 | yes | 0/0 | 6 | 0.0156 | 0.1 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0048 | 0.2 | 🟢 |
| MWSINDy | 2% | STRidge | 0.01 | 12 | no | 0/0 | 17 | 0.0052 | 0.1 | 🟢 |
| MWSINDy | 2% | SR3 | 0.05 | 4 | no | 0/2 | 2 | 0.7644 | 0.1 | 🔴 |
| MWSINDy | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0052 | 0.2 | 🟢 |
| MWSINDy | 2% | STRidge | 0.01 | 12 | yes | 0/0 | 17 | 0.0052 | 0.1 | 🟢 |
| MWSINDy | 2% | SR3 | 3.0 | 4 | yes | 0/0 | 6 | 0.0273 | 0.1 | 🟢 |
| MWSINDy | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0052 | 0.2 | 🟢 |
| MWSINDy | 5% | STRidge | 0.05 | 12 | no | 0/0 | 13 | 0.0078 | 0.1 | 🟢 |
| MWSINDy | 5% | SR3 | 0.01 | 4 | no | 0/0 | 4 | 0.1059 | 0.1 | 🟢 |
| MWSINDy | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0078 | 0.2 | 🟢 |
| MWSINDy | 5% | STRidge | 0.05 | 12 | yes | 0/0 | 13 | 0.0078 | 0.1 | 🟢 |
| MWSINDy | 5% | SR3 | 3.0 | 4 | yes | 0/0 | 6 | 0.1125 | 0.1 | 🟢 |
| MWSINDy | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0078 | 0.2 | 🟢 |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol/α | m | IBP | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.05 | 12 | no | 0/0 | 4 | 0.0048 | 0.2 | 🟢 |
| PDE-FIND | Clean | SR3 | 0.2 | 4 | no | 0/1 | 3 | 0.6251 | 0.1 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0048 | 0.3 | 🟢 |
| PDE-FIND | Clean | STRidge | 0.05 | 12 | yes | 0/0 | 4 | 0.0048 | 0.2 | 🟢 |
| PDE-FIND | Clean | SR3 | 7.0 | 8 | yes | 0/1 | 3 | 0.6930 | 0.1 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0048 | 0.3 | 🟢 |
| PDE-FIND | 2% | STRidge | 0.05 | 12 | no | 0/0 | 4 | 0.0045 | 0.1 | 🟢 |
| PDE-FIND | 2% | SR3 | 0.2 | 4 | no | 0/1 | 3 | 0.9061 | 0.1 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0051 | 0.3 | 🟢 |
| PDE-FIND | 2% | STRidge | 0.1 | 12 | yes | 0/0 | 4 | 0.0046 | 0.2 | 🟢 |
| PDE-FIND | 2% | SR3 | 4.0 | 8 | yes | 0/1 | 3 | 0.9352 | 0.1 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0052 | 0.3 | 🟢 |
| PDE-FIND | 5% | STRidge | 0.5 | 12 | no | 0/0 | 4 | 0.0053 | 0.2 | 🟢 |
| PDE-FIND | 5% | SR3 | 1.5 | 8 | no | 0/0 | 4 | 1.2172 | 0.1 | 🟢 |
| PDE-FIND | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0077 | 0.3 | 🟢 |
| PDE-FIND | 5% | STRidge | 0.5 | 12 | yes | 0/0 | 4 | 0.0054 | 0.2 | 🟢 |
| PDE-FIND | 5% | SR3 | 1.5 | 8 | yes | 0/0 | 4 | 1.2314 | 0.1 | 🟢 |
| PDE-FIND | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0078 | 0.3 | 🟢 |
| WSINDy | Clean | STRidge | 0.05 | 12 | no | 0/0 | 4 | 0.0048 | 0.1 | 🟢 |
| WSINDy | Clean | SR3 | 3.0 | 4 | no | 0/0 | 4 | 0.9619 | 0.1 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0048 | 0.2 | 🟢 |
| WSINDy | Clean | STRidge | 0.05 | 12 | yes | 0/0 | 4 | 0.0048 | 0.1 | 🟢 |
| WSINDy | Clean | SR3 | 0.5 | 4 | yes | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| WSINDy | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0048 | 0.2 | 🟢 |
| WSINDy | 2% | STRidge | 0.05 | 12 | no | 0/0 | 4 | 0.0045 | 0.1 | 🟢 |
| WSINDy | 2% | SR3 | 5.0 | 4 | no | 0/0 | 4 | 0.9999 | 0.1 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0051 | 0.2 | 🟢 |
| WSINDy | 2% | STRidge | 0.05 | 12 | yes | 0/0 | 4 | 0.0046 | 0.1 | 🟢 |
| WSINDy | 2% | SR3 | 0.5 | 4 | yes | 0/2 | 2 | 1.0081 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0052 | 0.2 | 🟢 |
| WSINDy | 5% | STRidge | 0.1 | 12 | no | 0/0 | 4 | 0.0053 | 0.1 | 🟢 |
| WSINDy | 5% | SR3 | 3.0 | 4 | no | 0/0 | 4 | 1.0036 | 0.1 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0077 | 0.2 | 🟢 |
| WSINDy | 5% | STRidge | 0.1 | 12 | yes | 0/0 | 4 | 0.0054 | 0.1 | 🟢 |
| WSINDy | 5% | SR3 | 0.1 | 4 | yes | 0/2 | 2 | 0.9547 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0078 | 0.2 | 🟢 |
| PolyJet | Clean | STRidge | 0.1 | 12 | no | 0/0 | 4 | 0.0048 | 0.3 | 🟢 |
| PolyJet | Clean | SR3 | 0.2 | 4 | no | 0/0 | 4 | 0.6667 | 0.2 | 🟢 |
| PolyJet | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0048 | 0.7 | 🟢 |
| PolyJet | Clean | STRidge | 0.5 | 12 | yes | 0/0 | 4 | 0.0048 | 0.3 | 🟢 |
| PolyJet | Clean | SR3 | 3.0 | 8 | yes | 0/0 | 4 | 0.7220 | 0.1 | 🟢 |
| PolyJet | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0048 | 0.7 | 🟢 |
| PolyJet | 2% | STRidge | 0.5 | 12 | no | 0/0 | 4 | 0.0045 | 0.3 | 🟢 |
| PolyJet | 2% | SR3 | 2.0 | 12 | no | 0/1 | 3 | 0.8118 | 0.1 | 🔴 |
| PolyJet | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0051 | 0.7 | 🟢 |
| PolyJet | 2% | STRidge | 0.2 | 12 | yes | 0/0 | 4 | 0.0046 | 0.3 | 🟢 |
| PolyJet | 2% | SR3 | 2.0 | 12 | yes | 0/1 | 3 | 0.8119 | 0.1 | 🔴 |
| PolyJet | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0052 | 0.7 | 🟢 |
| PolyJet | 5% | STRidge | 0.5 | 12 | no | 0/0 | 4 | 0.0053 | 0.3 | 🟢 |
| PolyJet | 5% | SR3 | 1.5 | 4 | no | 0/1 | 3 | 1.0222 | 0.2 | 🔴 |
| PolyJet | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0077 | 0.7 | 🟢 |
| PolyJet | 5% | STRidge | 0.5 | 12 | yes | 0/0 | 4 | 0.0054 | 0.3 | 🟢 |
| PolyJet | 5% | SR3 | 0.1 | 8 | yes | 0/1 | 3 | 1.0465 | 0.1 | 🔴 |
| PolyJet | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0078 | 0.7 | 🟢 |
| Multigraphs | Clean | STRidge | 0.2 | 12 | no | 0/0 | 3 | 0.0048 | 0.0 | 🟢 |
| Multigraphs | Clean | SR3 | 0.05 | 4 | no | 0/2 | 2 | 0.3769 | 0.0 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 3 | 0.0048 | 0.1 | 🟢 |
| Multigraphs | Clean | STRidge | 0.05 | 12 | yes | 0/0 | 3 | 0.0048 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.5 | 4 | yes | 0/4 | 0 | 1.0000 | 0.0 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 3 | 0.0048 | 0.1 | 🟢 |
| Multigraphs | 2% | STRidge | 0.2 | 12 | no | 0/0 | 3 | 0.0045 | 0.0 | 🟢 |
| Multigraphs | 2% | SR3 | 0.05 | 4 | no | 0/2 | 2 | 0.7706 | 0.0 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 3 | 0.0052 | 0.1 | 🟢 |
| Multigraphs | 2% | STRidge | 0.2 | 12 | yes | 0/0 | 3 | 0.0045 | 0.0 | 🟢 |
| Multigraphs | 2% | SR3 | 0.5 | 4 | yes | 0/2 | 2 | 0.8750 | 0.0 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 3 | 0.0052 | 0.1 | 🟢 |
| Multigraphs | 5% | STRidge | 0.5 | 12 | no | 0/0 | 3 | 0.0048 | 0.1 | 🟢 |
| Multigraphs | 5% | SR3 | 0.05 | 4 | no | 0/2 | 2 | 0.9558 | 0.0 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 3 | 0.0077 | 0.1 | 🟢 |
| Multigraphs | 5% | STRidge | 0.5 | 12 | yes | 0/0 | 3 | 0.0049 | 0.0 | 🟢 |
| Multigraphs | 5% | SR3 | 0.05 | 4 | yes | 0/2 | 2 | 0.9533 | 0.0 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 3 | 0.0078 | 0.1 | 🟢 |
| MWSINDy | Clean | STRidge | 0.01 | 12 | no | 0/0 | 6 | 0.0048 | 0.1 | 🟢 |
| MWSINDy | Clean | SR3 | 3.0 | 4 | no | 0/0 | 4 | 0.0063 | 0.1 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0048 | 0.2 | 🟢 |
| MWSINDy | Clean | STRidge | 0.05 | 12 | yes | 0/0 | 4 | 0.0048 | 0.1 | 🟢 |
| MWSINDy | Clean | SR3 | 3.0 | 4 | yes | 0/0 | 6 | 0.0156 | 0.1 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0048 | 0.2 | 🟢 |
| MWSINDy | 2% | STRidge | 0.01 | 12 | no | 0/0 | 13 | 0.0052 | 0.1 | 🟢 |
| MWSINDy | 2% | SR3 | 3.0 | 4 | no | 0/0 | 4 | 0.0086 | 0.1 | 🟢 |
| MWSINDy | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0052 | 0.2 | 🟢 |
| MWSINDy | 2% | STRidge | 0.01 | 12 | yes | 0/0 | 17 | 0.0052 | 0.1 | 🟢 |
| MWSINDy | 2% | SR3 | 3.0 | 4 | yes | 0/0 | 6 | 0.0273 | 0.1 | 🟢 |
| MWSINDy | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0052 | 0.2 | 🟢 |
| MWSINDy | 5% | STRidge | 0.05 | 12 | no | 0/0 | 13 | 0.0077 | 0.1 | 🟢 |
| MWSINDy | 5% | SR3 | 3.0 | 4 | no | 0/0 | 4 | 0.0453 | 0.1 | 🟢 |
| MWSINDy | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0077 | 0.2 | 🟢 |
| MWSINDy | 5% | STRidge | 0.05 | 12 | yes | 0/0 | 13 | 0.0078 | 0.1 | 🟢 |
| MWSINDy | 5% | SR3 | 3.0 | 4 | yes | 0/0 | 6 | 0.1125 | 0.1 | 🟢 |
| MWSINDy | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0078 | 0.2 | 🟢 |

---
*Source: `sweep_20260623_101642_ac_full.csv`*

*Run: 102.7 min (Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz), OMP_NUM_THREADS=1*
