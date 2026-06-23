# Porous Medium equation

The equation under consideration is

$$f_t = \nabla\cdot(f^2\nabla f) = f^2\Delta f + 2f\|\nabla f\|^2$$

which can be written in terms of $O(2)$-equivariant terms

| Operator | True coefficient |
|---|---|
| $f^2\cdot\Delta f$ | $+1.0$ |
| $f\cdot\|\nabla f\|^2$ | $+2.0$ |

or standard terms (polynomial jet basis)

| Operator | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|---|---|---|---|---|---|
| $f^2\cdot\partial_{xx} f$ | | 🟢 | 🟢 | 🟢 | 🟢 |
| $f^2\cdot\partial_{yy} f$ | | 🟢 | 🟢 | 🟢 | 🟢 |
| $f\cdot f_x^2$ | | 🔴 | 🟢 | 🟢 | 🟢 |
| $f\cdot f_y^2$ | | 🔴 | 🟢 | 🟢 | 🟢 |

All libraries are evaluated in the common polynomial-jet basis (105 monomials,
max order 4, max degree 3).  fp/fn/E are computed against the same 4-term
ground truth for every library.  PDE-FIND cannot represent $f\cdot f_x^2$ or $f\cdot f_y^2$ (products of derivatives), so its minimum fn is 2.  WSINDy represents PME exactly via the chain rule: $f_t = \frac13\partial_x^2(f^3) + \frac13\partial_y^2(f^3)$, which expands to all 4 jet terms.

Sweep configuration: all 4 libraries, 3 algorithms (STRidge, SR3, MSTLS),
3 methods (pointwise, weak-direct, weak-IBP), 2 backends (spectral, polynomial),
3 noise levels (clean, 2%, 5%).

### Legend for summary table

| Symbol | Meaning |
|---|---|
| 🟢 | Equivariant **clear** win |
| 🟩 | Equivariant **marginal** win |
| ⚪ | Tie  |
| 🟥 | Non-equivariant **marginal** win |
| 🔴 | Non-equivariant **clear** win |
| — | Both methods failed, no comparison |
| ✓ | Exact recovery (fp = fn = 0) |
| ✗ | Failed |


## Summary: equivariant vs best non-equivariant

| Noise | Method | Backend | Eq. wins? | Multigraphs fp/fn | Multigraphs E | Non-eq best | Non-eq fp/fn | Non-eq E |
|---|---|---|---|---|---|---|---|---|
| Clean | pointwise | spectral | 🟩 | 0/0 ✓ | 0.0009 | WSINDy | 0/0 ✓ | 0.0010 |
| Clean | pointwise | polynomial | 🟩 | 0/0 ✓ | 0.0006 | Polynomial jet | 0/0 ✓ | 0.0007 |
| Clean | weak (no IBP) | spectral | 🟩 | 0/0 ✓ | 0.0002 | Polynomial jet | 0/0 ✓ | 0.0004 |
| Clean | weak (no IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0002 | Polynomial jet | 0/0 ✓ | 0.0004 |
| Clean | weak (IBP) | spectral | 🟩 | 0/0 ✓ | 0.0002 | Polynomial jet | 0/0 ✓ | 0.0004 |
| Clean | weak (IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0002 | Polynomial jet | 0/0 ✓ | 0.0004 |
| 2% | pointwise | spectral | — | 0/4 ✗ | 1.0000 | — | — | — |
| 2% | pointwise | polynomial | 🔴 | 0/4 ✗ | 1.0000 | WSINDy | 0/0 ✓ | 0.0453 |
| 2% | weak (no IBP) | spectral | 🟥 | 0/0 ✓ | 0.0019 | WSINDy | 0/0 ✓ | 0.0012 |
| 2% | weak (no IBP) | polynomial | 🟥 | 0/0 ✓ | 0.0017 | WSINDy | 0/0 ✓ | 0.0012 |
| 2% | weak (IBP) | spectral | 🟥 | 0/0 ✓ | 0.0019 | WSINDy | 0/0 ✓ | 0.0012 |
| 2% | weak (IBP) | polynomial | 🟥 | 0/0 ✓ | 0.0017 | WSINDy | 0/0 ✓ | 0.0012 |
| 5% | pointwise | spectral | — | 0/4 ✗ | 1.0000 | — | — | — |
| 5% | pointwise | polynomial | — | 0/4 ✗ | 1.0000 | — | — | — |
| 5% | weak (no IBP) | spectral | 🟥 | 0/0 ✓ | 0.0022 | WSINDy | 0/0 ✓ | 0.0017 |
| 5% | weak (no IBP) | polynomial | 🔴 | 0/0 ✓ | 0.0051 | WSINDy | 0/0 ✓ | 0.0018 |
| 5% | weak (IBP) | spectral | ⚪ | 0/0 ✓ | 0.0022 | WSINDy | 0/0 ✓ | 0.0022 |
| 5% | weak (IBP) | polynomial | 🔴 | 0/0 ✓ | 0.0051 | WSINDy | 0/0 ✓ | 0.0018 |

## Pointwise formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol/α | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | Clean | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 2/4 | 2 | 0.2585 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.2 | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | 4/4 | 4 | 0.9946 | 🔴 |
| WSINDy | Clean | STRidge | 0.001 | 0/0 | 2 | 0.0010 | 🟢 |
| WSINDy | Clean | SR3 | 0.001 | 0/0 | 2 | 0.5622 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 0/0 | 2 | 0.0010 | 🟢 |
| WSINDy | 2% | STRidge | 0.2 | 0/4 | 0 | 1.0000 | 🔴 |
| WSINDy | 2% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 🔴 |
| WSINDy | 5% | STRidge | 0.2 | 0/4 | 0 | 1.0000 | 🔴 |
| WSINDy | 5% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | 4/4 | 4 | 0.9946 | 🔴 |
| Polynomial jet | Clean | STRidge | 0.005 | 0/0 | 4 | 0.0010 | 🟢 |
| Polynomial jet | Clean | SR3 | 0.001 | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | Clean | MSTLS | 0.25 | 0/0 | 4 | 0.0011 | 🟢 |
| Polynomial jet | 2% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 2% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 🔴 |
| Polynomial jet | 5% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 5% | SR3 | 0.01 | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 5% | MSTLS | 0.0 | 3/4 | 3 | 1.0083 | 🔴 |
| Multigraphs | Clean | STRidge | 0.05 | 0/0 | 2 | 0.0009 | 🟢 |
| Multigraphs | Clean | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 0/0 | 2 | 0.0011 | 🟢 |
| Multigraphs | 2% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | 2% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 🔴 |
| Multigraphs | 5% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | 5% | SR3 | 0.01 | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 2/4 | 2 | 0.9805 | 🔴 |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol/α | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | Clean | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 2/4 | 2 | 0.2585 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.01 | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | 2/3 | 4 | 0.5502 | 🔴 |
| WSINDy | Clean | STRidge | 0.005 | 0/0 | 2 | 0.0012 | 🟢 |
| WSINDy | Clean | SR3 | 0.001 | 0/0 | 2 | 0.6394 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 0/0 | 2 | 0.0012 | 🟢 |
| WSINDy | 2% | STRidge | 0.1 | 0/0 | 2 | 0.0453 | 🟢 |
| WSINDy | 2% | SR3 | 0.001 | 0/0 | 2 | 0.5716 | 🟢 |
| WSINDy | 2% | MSTLS | 0.0 | 0/4 | 0 | 1.0000 | 🔴 |
| WSINDy | 5% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 🔴 |
| WSINDy | 5% | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | 2/4 | 2 | 0.9805 | 🔴 |
| Polynomial jet | Clean | STRidge | 0.005 | 0/0 | 4 | 0.0007 | 🟢 |
| Polynomial jet | Clean | SR3 | 0.001 | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | Clean | MSTLS | 0.25 | 0/0 | 4 | 0.0008 | 🟢 |
| Polynomial jet | 2% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 2% | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 🔴 |
| Polynomial jet | 5% | STRidge | 0.2 | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 5% | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 5% | MSTLS | 0.0 | 2/4 | 2 | 0.9819 | 🔴 |
| Multigraphs | Clean | STRidge | 0.05 | 0/0 | 2 | 0.0006 | 🟢 |
| Multigraphs | Clean | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.25 | 0/0 | 2 | 0.0008 | 🟢 |
| Multigraphs | 2% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | 2% | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 🔴 |
| Multigraphs | 5% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | 5% | SR3 | 0.01 | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 2/4 | 2 | 0.9805 | 🔴 |

## Weak formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol/α | m | IBP | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.2 | 8 | no | 0/2 | 4 | 0.8389 | 🔴 |
| PDE-FIND | Clean | SR3 | 0.2 | 4 | no | 0/2 | 4 | 0.9712 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 4 | no | 2/4 | 2 | 0.2586 | 🔴 |
| PDE-FIND | Clean | STRidge | 0.5 | 8 | yes | 0/2 | 3 | 0.7949 | 🔴 |
| PDE-FIND | Clean | SR3 | 0.05 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 12 | yes | 5/2 | 7 | 0.3752 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.2 | 8 | no | 0/2 | 4 | 0.8205 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.1 | 4 | no | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 4 | no | 2/3 | 4 | 0.5617 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.2 | 8 | yes | 0/2 | 4 | 0.8205 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.05 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 12 | yes | 5/2 | 7 | 1.3576 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.1 | 4 | no | 0/4 | 2 | 1.0000 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.1 | 4 | no | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.75 | 4 | no | 2/2 | 6 | 0.2743 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.5 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.1 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.5 | 12 | yes | 2/2 | 5 | 1.5476 | 🔴 |
| WSINDy | Clean | STRidge | 0.05 | 12 | no | 0/0 | 2 | 0.0006 | 🟢 |
| WSINDy | Clean | SR3 | 0.01 | 4 | no | 0/0 | 2 | 0.6544 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0006 | 🟢 |
| WSINDy | Clean | STRidge | 0.05 | 12 | yes | 0/0 | 2 | 0.0006 | 🟢 |
| WSINDy | Clean | SR3 | 0.01 | 4 | yes | 0/0 | 2 | 0.6341 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0006 | 🟢 |
| WSINDy | 2% | STRidge | 0.5 | 12 | no | 0/0 | 2 | 0.0012 | 🟢 |
| WSINDy | 2% | SR3 | 0.005 | 4 | no | 0/0 | 2 | 0.4190 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0013 | 🟢 |
| WSINDy | 2% | STRidge | 0.5 | 12 | yes | 0/0 | 2 | 0.0012 | 🟢 |
| WSINDy | 2% | SR3 | 0.005 | 4 | yes | 0/0 | 2 | 0.4083 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0013 | 🟢 |
| WSINDy | 5% | STRidge | 1.0 | 12 | no | 0/0 | 2 | 0.0017 | 🟢 |
| WSINDy | 5% | SR3 | 0.05 | 4 | no | 0/0 | 2 | 0.6330 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0022 | 🟢 |
| WSINDy | 5% | STRidge | 0.2 | 8 | yes | 0/0 | 2 | 0.0043 | 🟢 |
| WSINDy | 5% | SR3 | 0.01 | 4 | yes | 0/0 | 2 | 0.5038 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0022 | 🟢 |
| Polynomial jet | Clean | STRidge | 0.2 | 12 | no | 0/0 | 4 | 0.0004 | 🟢 |
| Polynomial jet | Clean | SR3 | 0.01 | 4 | no | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0007 | 🟢 |
| Polynomial jet | Clean | STRidge | 0.2 | 12 | yes | 0/0 | 4 | 0.0004 | 🟢 |
| Polynomial jet | Clean | SR3 | 0.005 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0007 | 🟢 |
| Polynomial jet | 2% | STRidge | 0.5 | 8 | no | 0/0 | 4 | 0.0037 | 🟢 |
| Polynomial jet | 2% | SR3 | 0.005 | 4 | no | 0/0 | 4 | 0.5623 | 🟢 |
| Polynomial jet | 2% | MSTLS | 0.25 | 8 | no | 0/0 | 4 | 0.0043 | 🟢 |
| Polynomial jet | 2% | STRidge | 0.5 | 8 | yes | 0/0 | 4 | 0.0037 | 🟢 |
| Polynomial jet | 2% | SR3 | 0.01 | 4 | yes | 0/0 | 4 | 0.6677 | 🟢 |
| Polynomial jet | 2% | MSTLS | 0.25 | 8 | yes | 0/0 | 4 | 0.0043 | 🟢 |
| Polynomial jet | 5% | STRidge | 0.5 | 8 | no | 0/0 | 4 | 0.0202 | 🟢 |
| Polynomial jet | 5% | SR3 | 0.05 | 8 | no | 0/0 | 4 | 0.7120 | 🟢 |
| Polynomial jet | 5% | MSTLS | 0.0 | 8 | no | 2/2 | 4 | 50.9075 | 🔴 |
| Polynomial jet | 5% | STRidge | 0.2 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 5% | SR3 | 0.05 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 5% | MSTLS | 0.0 | 4 | yes | 5/4 | 5 | 7.3139 | 🔴 |
| Multigraphs | Clean | STRidge | 0.2 | 12 | no | 0/0 | 2 | 0.0002 | 🟢 |
| Multigraphs | Clean | SR3 | 0.05 | 4 | no | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 2 | 0.0007 | 🟢 |
| Multigraphs | Clean | STRidge | 0.2 | 12 | yes | 0/0 | 2 | 0.0002 | 🟢 |
| Multigraphs | Clean | SR3 | 0.05 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 12 | yes | 0/0 | 2 | 0.0007 | 🟢 |
| Multigraphs | 2% | STRidge | 0.5 | 8 | no | 0/0 | 2 | 0.0026 | 🟢 |
| Multigraphs | 2% | SR3 | 0.05 | 8 | no | 0/2 | 1 | 0.7542 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 12 | no | 0/0 | 2 | 0.0019 | 🟢 |
| Multigraphs | 2% | STRidge | 0.5 | 8 | yes | 0/0 | 2 | 0.0026 | 🟢 |
| Multigraphs | 2% | SR3 | 0.05 | 8 | yes | 0/2 | 1 | 0.7543 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 12 | yes | 0/0 | 2 | 0.0019 | 🟢 |
| Multigraphs | 5% | STRidge | 1.0 | 8 | no | 0/0 | 2 | 0.0375 | 🟢 |
| Multigraphs | 5% | SR3 | 0.05 | 4 | no | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 12 | no | 0/0 | 2 | 0.0022 | 🟢 |
| Multigraphs | 5% | STRidge | 0.5 | 8 | yes | 0/0 | 2 | 0.0375 | 🟢 |
| Multigraphs | 5% | SR3 | 0.05 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 12 | yes | 0/0 | 2 | 0.0022 | 🟢 |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol/α | m | IBP | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.2 | 8 | no | 0/2 | 4 | 0.8389 | 🔴 |
| PDE-FIND | Clean | SR3 | 0.2 | 4 | no | 0/2 | 4 | 0.9713 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 4 | no | 2/4 | 2 | 0.2586 | 🔴 |
| PDE-FIND | Clean | STRidge | 0.5 | 8 | yes | 0/2 | 3 | 0.7949 | 🔴 |
| PDE-FIND | Clean | SR3 | 0.01 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 12 | yes | 5/2 | 7 | 0.3755 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.2 | 8 | no | 0/2 | 4 | 0.8144 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.1 | 4 | no | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 4 | no | 2/4 | 2 | 0.9287 | 🔴 |
| PDE-FIND | 2% | STRidge | 1.0 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.05 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 12 | yes | 5/2 | 8 | 4.6131 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.1 | 4 | no | 0/4 | 2 | 1.0000 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.5 | 8 | no | 0/2 | 4 | 0.9705 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.5 | 4 | no | 2/3 | 4 | 0.5028 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.1 | 4 | yes | 0/2 | 3 | 0.7244 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.1 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | 12 | yes | 5/2 | 8 | 2.8699 | 🔴 |
| WSINDy | Clean | STRidge | 0.05 | 12 | no | 0/0 | 2 | 0.0006 | 🟢 |
| WSINDy | Clean | SR3 | 0.01 | 4 | no | 0/0 | 2 | 0.6575 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0006 | 🟢 |
| WSINDy | Clean | STRidge | 0.05 | 12 | yes | 0/0 | 2 | 0.0006 | 🟢 |
| WSINDy | Clean | SR3 | 0.01 | 4 | yes | 0/0 | 2 | 0.6369 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0006 | 🟢 |
| WSINDy | 2% | STRidge | 0.2 | 12 | no | 0/0 | 2 | 0.0012 | 🟢 |
| WSINDy | 2% | SR3 | 0.01 | 4 | no | 0/0 | 2 | 0.4963 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0013 | 🟢 |
| WSINDy | 2% | STRidge | 1.0 | 12 | yes | 0/0 | 2 | 0.0012 | 🟢 |
| WSINDy | 2% | SR3 | 0.005 | 4 | yes | 0/0 | 2 | 0.4379 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0013 | 🟢 |
| WSINDy | 5% | STRidge | 0.5 | 12 | no | 0/0 | 2 | 0.0018 | 🟢 |
| WSINDy | 5% | SR3 | 0.005 | 4 | no | 0/0 | 2 | 0.3573 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0022 | 🟢 |
| WSINDy | 5% | STRidge | 1.0 | 12 | yes | 0/0 | 2 | 0.0018 | 🟢 |
| WSINDy | 5% | SR3 | 0.005 | 4 | yes | 0/0 | 2 | 0.3113 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0022 | 🟢 |
| Polynomial jet | Clean | STRidge | 0.5 | 12 | no | 0/0 | 4 | 0.0004 | 🟢 |
| Polynomial jet | Clean | SR3 | 0.01 | 4 | no | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0006 | 🟢 |
| Polynomial jet | Clean | STRidge | 0.2 | 12 | yes | 0/0 | 4 | 0.0004 | 🟢 |
| Polynomial jet | Clean | SR3 | 0.01 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0006 | 🟢 |
| Polynomial jet | 2% | STRidge | 0.5 | 8 | no | 0/0 | 4 | 0.0040 | 🟢 |
| Polynomial jet | 2% | SR3 | 0.005 | 4 | no | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 2% | MSTLS | 0.25 | 8 | no | 0/0 | 4 | 0.0048 | 🟢 |
| Polynomial jet | 2% | STRidge | 0.5 | 8 | yes | 0/0 | 4 | 0.0040 | 🟢 |
| Polynomial jet | 2% | SR3 | 0.005 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0209 | 🟢 |
| Polynomial jet | 5% | STRidge | 0.1 | 4 | no | 0/0 | 4 | 0.0443 | 🟢 |
| Polynomial jet | 5% | SR3 | 0.005 | 4 | no | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 5% | MSTLS | 0.25 | 4 | no | 0/0 | 4 | 0.0352 | 🟢 |
| Polynomial jet | 5% | STRidge | 0.1 | 4 | yes | 0/0 | 4 | 0.0413 | 🟢 |
| Polynomial jet | 5% | SR3 | 0.01 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| Polynomial jet | 5% | MSTLS | 0.25 | 4 | yes | 3/0 | 7 | 14.2391 | 🔴 |
| Multigraphs | Clean | STRidge | 0.2 | 12 | no | 0/0 | 2 | 0.0002 | 🟢 |
| Multigraphs | Clean | SR3 | 0.05 | 4 | no | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 2 | 0.0006 | 🟢 |
| Multigraphs | Clean | STRidge | 0.2 | 12 | yes | 0/0 | 2 | 0.0002 | 🟢 |
| Multigraphs | Clean | SR3 | 0.05 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 12 | yes | 0/0 | 2 | 0.0006 | 🟢 |
| Multigraphs | 2% | STRidge | 1.0 | 4 | no | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | 4 | no | 0/2 | 1 | 0.9992 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 4 | no | 0/0 | 2 | 0.0017 | 🟢 |
| Multigraphs | 2% | STRidge | 0.5 | 4 | yes | 0/4 | 0 | 1.0000 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | 4 | yes | 0/2 | 1 | 0.9992 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 4 | yes | 0/0 | 2 | 0.0017 | 🟢 |
| Multigraphs | 5% | STRidge | 1.0 | 12 | no | 0/0 | 2 | 0.0051 | 🟢 |
| Multigraphs | 5% | SR3 | 0.05 | 8 | no | 0/2 | 1 | 0.9188 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 12 | no | 0/0 | 2 | 0.0053 | 🟢 |
| Multigraphs | 5% | STRidge | 1.0 | 12 | yes | 0/0 | 2 | 0.0051 | 🟢 |
| Multigraphs | 5% | SR3 | 0.05 | 8 | yes | 0/2 | 1 | 0.9191 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 12 | yes | 0/0 | 2 | 0.0053 | 🟢 |
