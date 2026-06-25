# Porous Medium equation

The equation under consideration is

$$f_t = \nabla\cdot(f^2\nabla f) = f^2\Delta f + 2f\|\nabla f\|^2$$

which can be written in terms of $O(2)$-equivariant terms

| Operator | True coefficient |
|---|---|
| $f^2\cdot\Delta f$ | $+1.0$ |
| $f\cdot\|\nabla f\|^2$ | $+2.0$ |

or standard terms (PolyJet basis)

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
| Clean | pointwise | spectral | 🟩 | 0/0 ✓ | 0.0009 | WSINDy | 0/0 ✓ | 0.0010 | 0/0 ✓ | 0.0011 |
| Clean | pointwise | polynomial | 🟩 | 0/0 ✓ | 0.0006 | PolyJet | 0/0 ✓ | 0.0007 | 0/0 ✓ | 0.0008 |
| Clean | weak (no IBP) | spectral | 🟩 | 0/0 ✓ | 0.0004 | PolyJet | 0/0 ✓ | 0.0006 | 0/0 ✓ | 0.0008 |
| Clean | weak (no IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0003 | PolyJet | 0/0 ✓ | 0.0005 | 0/0 ✓ | 0.0007 |
| Clean | weak (IBP) | spectral | 🟩 | 0/0 ✓ | 0.0004 | PolyJet | 0/0 ✓ | 0.0006 | 0/0 ✓ | 0.0008 |
| Clean | weak (IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0003 | PolyJet | 0/0 ✓ | 0.0005 | 0/0 ✓ | 0.0007 |
| 2% | pointwise | spectral | — | 0/4 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | 0/4 ✗ | 1.0000 |
| 2% | pointwise | polynomial | 🔴 | 0/4 ✗ | 1.0000 | WSINDy | 0/0 ✓ | 0.0453 | 0/0 ✓ | 0.0664 |
| 2% | weak (no IBP) | spectral | 🟥 | 0/0 ✓ | 0.0019 | WSINDy | 0/0 ✓ | 0.0014 | 0/0 ✓ | 0.0020 |
| 2% | weak (no IBP) | polynomial | 🟥 | 0/0 ✓ | 0.0024 | WSINDy | 0/0 ✓ | 0.0014 | 0/0 ✓ | 0.0024 |
| 2% | weak (IBP) | spectral | 🟥 | 0/0 ✓ | 0.0019 | WSINDy | 0/0 ✓ | 0.0014 | 0/0 ✓ | 0.0020 |
| 2% | weak (IBP) | polynomial | 🟥 | 0/0 ✓ | 0.0024 | WSINDy | 0/0 ✓ | 0.0014 | 0/0 ✓ | 0.0024 |
| 5% | pointwise | spectral | — | 0/4 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | 0/4 ✗ | 1.0000 |
| 5% | pointwise | polynomial | — | 0/4 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | 0/4 ✗ | 1.0000 |
| 5% | weak (no IBP) | spectral | 🟩 | 0/0 ✓ | 0.0022 | WSINDy | 0/0 ✓ | 0.0022 | 0/0 ✓ | 0.0022 |
| 5% | weak (no IBP) | polynomial | 🟥 | 0/0 ✓ | 0.0051 | WSINDy | 0/0 ✓ | 0.0022 | 0/0 ✓ | 0.0055 |
| 5% | weak (IBP) | spectral | 🟩 | 0/0 ✓ | 0.0022 | WSINDy | 0/0 ✓ | 0.0022 | 0/0 ✓ | 0.0022 |
| 5% | weak (IBP) | polynomial | 🟥 | 0/0 ✓ | 0.0051 | WSINDy | 0/0 ✓ | 0.0022 | 0/0 ✓ | 0.0055 |

## Pointwise formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol/α | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | Clean | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 2/4 | 2 | 0.2585 | 2.7 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 2.2 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.2 | 0/4 | 0 | 1.0000 | 0.4 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 0.3 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | 4/4 | 4 | 0.9946 | 3.0 | 🔴 |
| WSINDy | Clean | STRidge | 0.001 | 0/0 | 2 | 0.0010 | 0.3 | 🟢 |
| WSINDy | Clean | SR3 | 0.001 | 0/0 | 2 | 0.5622 | 0.1 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 0/0 | 2 | 0.0010 | 0.4 | 🟢 |
| WSINDy | 2% | STRidge | 0.2 | 0/4 | 0 | 1.0000 | 0.2 | 🔴 |
| WSINDy | 2% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 0.8 | 🔴 |
| WSINDy | 5% | STRidge | 0.2 | 0/4 | 0 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 5% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | 4/4 | 4 | 0.9946 | 0.8 | 🔴 |
| PolyJet | Clean | STRidge | 0.005 | 0/0 | 4 | 0.0010 | 13.3 | 🟢 |
| PolyJet | Clean | SR3 | 0.001 | 0/4 | 0 | 1.0000 | 1.0 | 🔴 |
| PolyJet | Clean | MSTLS | 0.25 | 0/0 | 4 | 0.0011 | 29.1 | 🟢 |
| PolyJet | 2% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 1.0 | 🔴 |
| PolyJet | 2% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 3.9 | 🔴 |
| PolyJet | 5% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 0.8 | 🔴 |
| PolyJet | 5% | SR3 | 0.01 | 0/4 | 0 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | 3/4 | 3 | 1.0083 | 2.9 | 🔴 |
| Multigraphs | Clean | STRidge | 0.05 | 0/0 | 2 | 0.0009 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 0.0 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 0/0 | 2 | 0.0011 | 0.2 | 🟢 |
| Multigraphs | 2% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 0.2 | 🔴 |
| Multigraphs | 5% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 0.2 | 🔴 |
| Multigraphs | 5% | SR3 | 0.01 | 0/4 | 0 | 1.0000 | 0.0 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 2/4 | 2 | 0.9805 | 0.2 | 🔴 |
| MWSINDy | Clean | STRidge | 0.0005 | 0/0 | 2 | 0.0011 | 0.3 | 🟢 |
| MWSINDy | Clean | SR3 | 0.001 | 0/0 | 2 | 0.0011 | 0.1 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.25 | 0/0 | 2 | 0.0011 | 0.5 | 🟢 |
| MWSINDy | 2% | STRidge | 0.2 | 0/4 | 0 | 1.0000 | 0.2 | 🔴 |
| MWSINDy | 2% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| MWSINDy | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 0.8 | 🔴 |
| MWSINDy | 5% | STRidge | 0.2 | 0/4 | 0 | 1.0000 | 0.2 | 🔴 |
| MWSINDy | 5% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| MWSINDy | 5% | MSTLS | 0.0 | 4/4 | 4 | 0.9946 | 0.8 | 🔴 |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol/α | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | Clean | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 2/4 | 2 | 0.2585 | 2.9 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 0.4 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.01 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 1.9 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 0.4 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.05 | 0/4 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | 2/3 | 4 | 0.5502 | 2.3 | 🔴 |
| WSINDy | Clean | STRidge | 0.005 | 0/0 | 2 | 0.0012 | 0.2 | 🟢 |
| WSINDy | Clean | SR3 | 0.001 | 0/0 | 2 | 0.6394 | 0.1 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 0/0 | 2 | 0.0012 | 0.5 | 🟢 |
| WSINDy | 2% | STRidge | 0.1 | 0/0 | 2 | 0.0453 | 0.2 | 🟢 |
| WSINDy | 2% | SR3 | 0.001 | 0/0 | 2 | 0.5716 | 0.1 | 🟢 |
| WSINDy | 2% | MSTLS | 0.0 | 0/4 | 0 | 1.0000 | 0.7 | 🔴 |
| WSINDy | 5% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 0.2 | 🔴 |
| WSINDy | 5% | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | 2/4 | 2 | 0.9805 | 0.7 | 🔴 |
| PolyJet | Clean | STRidge | 0.005 | 0/0 | 4 | 0.0007 | 12.5 | 🟢 |
| PolyJet | Clean | SR3 | 0.001 | 0/4 | 0 | 1.0000 | 0.3 | 🔴 |
| PolyJet | Clean | MSTLS | 0.25 | 0/0 | 4 | 0.0008 | 2.1 | 🟢 |
| PolyJet | 2% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 0.8 | 🔴 |
| PolyJet | 2% | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 0.2 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 3.2 | 🔴 |
| PolyJet | 5% | STRidge | 0.2 | 0/4 | 0 | 1.0000 | 0.8 | 🔴 |
| PolyJet | 5% | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 0.2 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | 2/4 | 2 | 0.9819 | 3.5 | 🔴 |
| Multigraphs | Clean | STRidge | 0.05 | 0/0 | 2 | 0.0006 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 0.0 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.25 | 0/0 | 2 | 0.0008 | 0.3 | 🟢 |
| Multigraphs | 2% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 0.0 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 2/4 | 2 | 0.9285 | 0.2 | 🔴 |
| Multigraphs | 5% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | SR3 | 0.01 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 2/4 | 2 | 0.9805 | 0.2 | 🔴 |
| MWSINDy | Clean | STRidge | 0.001 | 0/0 | 2 | 0.0008 | 0.3 | 🟢 |
| MWSINDy | Clean | SR3 | 0.001 | 0/0 | 2 | 0.0008 | 0.1 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.25 | 0/0 | 2 | 0.0008 | 0.5 | 🟢 |
| MWSINDy | 2% | STRidge | 0.05 | 0/0 | 4 | 0.0664 | 0.3 | 🟢 |
| MWSINDy | 2% | SR3 | 0.001 | 0/0 | 2 | 0.0664 | 0.1 | 🟢 |
| MWSINDy | 2% | MSTLS | 0.0 | 0/4 | 0 | 1.0000 | 0.7 | 🔴 |
| MWSINDy | 5% | STRidge | 0.5 | 0/4 | 0 | 1.0000 | 0.2 | 🔴 |
| MWSINDy | 5% | SR3 | 0.005 | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| MWSINDy | 5% | MSTLS | 0.0 | 2/4 | 2 | 0.9805 | 0.9 | 🔴 |

## Weak formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol/α | m | IBP | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.05 | 4 | no | 0/2 | 4 | 0.8211 | 0.4 | 🔴 |
| PDE-FIND | Clean | SR3 | 1.0 | 8 | no | 0/2 | 4 | 0.9638 | 0.1 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 4 | no | 2/4 | 2 | 0.2585 | 1.3 | 🔴 |
| PDE-FIND | Clean | STRidge | 0.05 | 4 | yes | 0/2 | 4 | 0.8211 | 0.4 | 🔴 |
| PDE-FIND | Clean | SR3 | 0.01 | 4 | yes | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 4 | yes | 2/4 | 2 | 0.2585 | 1.1 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.05 | 4 | no | 0/2 | 4 | 0.7986 | 0.4 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.1 | 4 | no | 0/4 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 4 | no | 2/3 | 4 | 0.5651 | 1.0 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.05 | 4 | yes | 0/2 | 4 | 0.7986 | 0.5 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.05 | 4 | yes | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 4 | yes | 2/3 | 4 | 0.5560 | 1.2 | 🔴 |
| PDE-FIND | 5% | STRidge | 2.0 | 4 | no | 0/4 | 0 | 1.0000 | 0.4 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.1 | 4 | no | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.5 | 8 | no | 2/3 | 4 | 0.5134 | 1.0 | 🔴 |
| PDE-FIND | 5% | STRidge | 2.0 | 4 | yes | 0/4 | 0 | 1.0000 | 0.4 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.05 | 4 | yes | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.5 | 8 | yes | 2/3 | 4 | 0.5134 | 1.0 | 🔴 |
| WSINDy | Clean | STRidge | 0.05 | 12 | no | 0/0 | 2 | 0.0007 | 0.1 | 🟢 |
| WSINDy | Clean | SR3 | 0.01 | 4 | no | 0/0 | 2 | 0.6257 | 0.1 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0007 | 0.3 | 🟢 |
| WSINDy | Clean | STRidge | 0.05 | 12 | yes | 0/0 | 2 | 0.0007 | 0.2 | 🟢 |
| WSINDy | Clean | SR3 | 0.01 | 4 | yes | 0/0 | 2 | 0.6015 | 0.1 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0007 | 0.4 | 🟢 |
| WSINDy | 2% | STRidge | 0.2 | 12 | no | 0/0 | 2 | 0.0014 | 0.1 | 🟢 |
| WSINDy | 2% | SR3 | 0.005 | 4 | no | 0/0 | 2 | 0.4157 | 0.1 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0016 | 0.4 | 🟢 |
| WSINDy | 2% | STRidge | 0.2 | 12 | yes | 0/0 | 2 | 0.0014 | 0.1 | 🟢 |
| WSINDy | 2% | SR3 | 0.005 | 4 | yes | 0/0 | 2 | 0.4120 | 0.1 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0016 | 0.3 | 🟢 |
| WSINDy | 5% | STRidge | 0.5 | 12 | no | 0/0 | 2 | 0.0022 | 0.2 | 🟢 |
| WSINDy | 5% | SR3 | 0.1 | 8 | no | 0/0 | 2 | 0.7370 | 0.1 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0028 | 0.4 | 🟢 |
| WSINDy | 5% | STRidge | 0.5 | 12 | yes | 0/0 | 2 | 0.0022 | 0.1 | 🟢 |
| WSINDy | 5% | SR3 | 0.1 | 8 | yes | 0/0 | 2 | 0.7371 | 0.1 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0028 | 0.4 | 🟢 |
| PolyJet | Clean | STRidge | 0.1 | 12 | no | 0/0 | 4 | 0.0006 | 0.5 | 🟢 |
| PolyJet | Clean | SR3 | 0.005 | 4 | no | 0/4 | 0 | 1.0000 | 0.4 | 🔴 |
| PolyJet | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0008 | 1.1 | 🟢 |
| PolyJet | Clean | STRidge | 0.1 | 12 | yes | 0/0 | 4 | 0.0006 | 1.0 | 🟢 |
| PolyJet | Clean | SR3 | 0.005 | 4 | yes | 0/4 | 0 | 1.0000 | 1.4 | 🔴 |
| PolyJet | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0008 | 1.5 | 🟢 |
| PolyJet | 2% | STRidge | 0.2 | 8 | no | 0/0 | 4 | 0.0085 | 0.8 | 🟢 |
| PolyJet | 2% | SR3 | 0.05 | 12 | no | 0/2 | 2 | 0.7247 | 0.3 | 🔴 |
| PolyJet | 2% | MSTLS | 0.25 | 8 | no | 0/0 | 4 | 0.0023 | 3.9 | 🟢 |
| PolyJet | 2% | STRidge | 0.2 | 8 | yes | 0/0 | 4 | 0.0085 | 1.0 | 🟢 |
| PolyJet | 2% | SR3 | 0.05 | 12 | yes | 0/2 | 2 | 0.7246 | 0.2 | 🔴 |
| PolyJet | 2% | MSTLS | 0.25 | 8 | yes | 0/0 | 4 | 0.0023 | 3.4 | 🟢 |
| PolyJet | 5% | STRidge | 0.2 | 8 | no | 0/0 | 4 | 0.0184 | 0.7 | 🟢 |
| PolyJet | 5% | SR3 | 0.01 | 4 | no | 0/4 | 0 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 5% | MSTLS | 0.25 | 8 | no | 2/0 | 6 | 5.0653 | 3.6 | 🔴 |
| PolyJet | 5% | STRidge | 0.2 | 8 | yes | 0/0 | 4 | 0.0184 | 1.3 | 🟢 |
| PolyJet | 5% | SR3 | 0.01 | 4 | yes | 0/4 | 0 | 1.0000 | 0.2 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | 8 | yes | 2/0 | 6 | 5.0653 | 2.5 | 🔴 |
| Multigraphs | Clean | STRidge | 0.5 | 12 | no | 0/0 | 2 | 0.0004 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.01 | 4 | no | 0/4 | 0 | 1.0000 | 0.0 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 2 | 0.0008 | 0.1 | 🟢 |
| Multigraphs | Clean | STRidge | 0.5 | 12 | yes | 0/0 | 2 | 0.0004 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.01 | 4 | yes | 0/4 | 0 | 1.0000 | 0.0 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 12 | yes | 0/0 | 2 | 0.0008 | 0.1 | 🟢 |
| Multigraphs | 2% | STRidge | 1.5 | 12 | no | 0/0 | 2 | 0.0019 | 0.1 | 🟢 |
| Multigraphs | 2% | SR3 | 0.05 | 8 | no | 0/2 | 1 | 0.7686 | 0.0 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 12 | no | 0/0 | 2 | 0.0020 | 0.2 | 🟢 |
| Multigraphs | 2% | STRidge | 1.5 | 12 | yes | 0/0 | 2 | 0.0019 | 0.1 | 🟢 |
| Multigraphs | 2% | SR3 | 0.05 | 8 | yes | 0/2 | 1 | 0.7685 | 0.0 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 12 | yes | 0/0 | 2 | 0.0020 | 0.2 | 🟢 |
| Multigraphs | 5% | STRidge | 7.0 | 12 | no | 0/0 | 2 | 0.0025 | 0.1 | 🟢 |
| Multigraphs | 5% | SR3 | 0.05 | 8 | no | 0/0 | 2 | 0.7705 | 0.0 | 🟢 |
| Multigraphs | 5% | MSTLS | 0.0 | 12 | no | 0/0 | 2 | 0.0022 | 0.2 | 🟢 |
| Multigraphs | 5% | STRidge | 7.0 | 12 | yes | 0/0 | 2 | 0.0025 | 0.1 | 🟢 |
| Multigraphs | 5% | SR3 | 0.05 | 8 | yes | 0/0 | 2 | 0.7703 | 0.0 | 🟢 |
| Multigraphs | 5% | MSTLS | 0.0 | 12 | yes | 0/0 | 2 | 0.0022 | 0.2 | 🟢 |
| MWSINDy | Clean | STRidge | 0.005 | 12 | no | 0/0 | 2 | 0.0008 | 0.1 | 🟢 |
| MWSINDy | Clean | SR3 | 0.05 | 12 | no | 0/0 | 2 | 0.0008 | 0.1 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0008 | 0.3 | 🟢 |
| MWSINDy | Clean | STRidge | 0.005 | 12 | yes | 0/0 | 2 | 0.0008 | 0.1 | 🟢 |
| MWSINDy | Clean | SR3 | 0.05 | 12 | yes | 0/0 | 2 | 0.0008 | 0.1 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0008 | 0.4 | 🟢 |
| MWSINDy | 2% | STRidge | 0.5 | 12 | no | 0/0 | 2 | 0.0020 | 0.2 | 🟢 |
| MWSINDy | 2% | SR3 | 0.05 | 12 | no | 0/0 | 2 | 0.0020 | 0.1 | 🟢 |
| MWSINDy | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0020 | 0.3 | 🟢 |
| MWSINDy | 2% | STRidge | 0.5 | 12 | yes | 0/0 | 2 | 0.0020 | 0.2 | 🟢 |
| MWSINDy | 2% | SR3 | 0.05 | 12 | yes | 0/0 | 2 | 0.0020 | 0.1 | 🟢 |
| MWSINDy | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0020 | 0.4 | 🟢 |
| MWSINDy | 5% | STRidge | 1.0 | 12 | no | 0/0 | 2 | 0.0022 | 0.2 | 🟢 |
| MWSINDy | 5% | SR3 | 0.1 | 8 | no | 0/0 | 2 | 0.0024 | 0.1 | 🟢 |
| MWSINDy | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0022 | 0.4 | 🟢 |
| MWSINDy | 5% | STRidge | 1.0 | 12 | yes | 0/0 | 2 | 0.0022 | 0.1 | 🟢 |
| MWSINDy | 5% | SR3 | 0.1 | 8 | yes | 0/0 | 2 | 0.0024 | 0.1 | 🟢 |
| MWSINDy | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0022 | 0.3 | 🟢 |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol/α | m | IBP | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.5 | 8 | no | 0/3 | 3 | 0.8943 | 0.5 | 🔴 |
| PDE-FIND | Clean | SR3 | 1.0 | 8 | no | 0/2 | 4 | 0.9638 | 0.1 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 4 | no | 2/4 | 2 | 0.2585 | 1.1 | 🔴 |
| PDE-FIND | Clean | STRidge | 0.05 | 4 | yes | 0/2 | 4 | 0.8210 | 0.5 | 🔴 |
| PDE-FIND | Clean | SR3 | 0.005 | 4 | yes | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 4 | yes | 2/4 | 2 | 0.2585 | 1.0 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.05 | 4 | no | 0/2 | 4 | 0.7918 | 0.4 | 🔴 |
| PDE-FIND | 2% | SR3 | 1.0 | 8 | no | 0/2 | 4 | 0.9605 | 0.1 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.5 | 4 | no | 2/3 | 4 | 0.5627 | 1.2 | 🔴 |
| PDE-FIND | 2% | STRidge | 0.05 | 4 | yes | 0/2 | 4 | 0.7918 | 0.6 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.05 | 4 | yes | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.5 | 4 | yes | 2/3 | 4 | 0.5627 | 1.3 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.2 | 8 | no | 0/2 | 4 | 0.8073 | 0.3 | 🔴 |
| PDE-FIND | 5% | SR3 | 1.0 | 12 | no | 0/2 | 4 | 0.9777 | 0.1 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.75 | 4 | no | 2/3 | 4 | 0.5069 | 1.2 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.2 | 8 | yes | 0/3 | 3 | 0.8824 | 0.3 | 🔴 |
| PDE-FIND | 5% | SR3 | 1.0 | 12 | yes | 0/2 | 4 | 0.9777 | 0.1 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | 4 | yes | 2/4 | 2 | 0.9804 | 1.2 | 🔴 |
| WSINDy | Clean | STRidge | 0.05 | 12 | no | 0/0 | 2 | 0.0007 | 0.2 | 🟢 |
| WSINDy | Clean | SR3 | 0.005 | 4 | no | 0/0 | 2 | 0.6294 | 0.2 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0007 | 0.5 | 🟢 |
| WSINDy | Clean | STRidge | 0.05 | 12 | yes | 0/0 | 2 | 0.0007 | 0.1 | 🟢 |
| WSINDy | Clean | SR3 | 0.01 | 4 | yes | 0/0 | 2 | 0.6023 | 0.1 | 🟢 |
| WSINDy | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0007 | 0.3 | 🟢 |
| WSINDy | 2% | STRidge | 0.2 | 12 | no | 0/0 | 2 | 0.0014 | 0.1 | 🟢 |
| WSINDy | 2% | SR3 | 0.005 | 4 | no | 0/0 | 2 | 0.4295 | 0.1 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0016 | 0.3 | 🟢 |
| WSINDy | 2% | STRidge | 0.2 | 12 | yes | 0/0 | 2 | 0.0014 | 0.2 | 🟢 |
| WSINDy | 2% | SR3 | 0.005 | 4 | yes | 0/0 | 2 | 0.3454 | 0.1 | 🟢 |
| WSINDy | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0016 | 0.4 | 🟢 |
| WSINDy | 5% | STRidge | 1.0 | 12 | no | 0/0 | 2 | 0.0022 | 0.2 | 🟢 |
| WSINDy | 5% | SR3 | 0.005 | 4 | no | 0/0 | 2 | 0.3677 | 0.1 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0028 | 0.4 | 🟢 |
| WSINDy | 5% | STRidge | 1.0 | 12 | yes | 0/0 | 2 | 0.0022 | 0.1 | 🟢 |
| WSINDy | 5% | SR3 | 0.005 | 4 | yes | 0/0 | 2 | 0.3277 | 0.1 | 🟢 |
| WSINDy | 5% | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0028 | 0.4 | 🟢 |
| PolyJet | Clean | STRidge | 0.2 | 12 | no | 0/0 | 4 | 0.0005 | 0.7 | 🟢 |
| PolyJet | Clean | SR3 | 0.005 | 4 | no | 0/4 | 0 | 1.0000 | 0.3 | 🔴 |
| PolyJet | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 4 | 0.0007 | 1.8 | 🟢 |
| PolyJet | Clean | STRidge | 0.2 | 12 | yes | 0/0 | 4 | 0.0005 | 1.0 | 🟢 |
| PolyJet | Clean | SR3 | 0.005 | 4 | yes | 0/4 | 0 | 1.0000 | 0.4 | 🔴 |
| PolyJet | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 4 | 0.0007 | 1.3 | 🟢 |
| PolyJet | 2% | STRidge | 0.2 | 8 | no | 0/0 | 4 | 0.0091 | 1.2 | 🟢 |
| PolyJet | 2% | SR3 | 0.005 | 4 | no | 0/4 | 0 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 2% | MSTLS | 0.25 | 8 | no | 0/0 | 4 | 0.0027 | 2.2 | 🟢 |
| PolyJet | 2% | STRidge | 0.5 | 8 | yes | 0/0 | 4 | 0.0091 | 0.9 | 🟢 |
| PolyJet | 2% | SR3 | 0.005 | 4 | yes | 0/4 | 0 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 2% | MSTLS | 0.25 | 8 | yes | 0/0 | 4 | 0.0027 | 2.3 | 🟢 |
| PolyJet | 5% | STRidge | 0.5 | 8 | no | 0/0 | 4 | 0.0082 | 0.8 | 🟢 |
| PolyJet | 5% | SR3 | 0.005 | 4 | no | 0/4 | 0 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 5% | MSTLS | 0.25 | 8 | no | 0/0 | 4 | 0.0080 | 3.1 | 🟢 |
| PolyJet | 5% | STRidge | 1.0 | 8 | yes | 0/0 | 4 | 0.0082 | 0.6 | 🟢 |
| PolyJet | 5% | SR3 | 0.05 | 4 | yes | 0/4 | 0 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 5% | MSTLS | 0.25 | 8 | yes | 0/0 | 4 | 0.0080 | 2.3 | 🟢 |
| Multigraphs | Clean | STRidge | 0.5 | 12 | no | 0/0 | 2 | 0.0003 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.01 | 4 | no | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 2 | 0.0007 | 0.1 | 🟢 |
| Multigraphs | Clean | STRidge | 0.5 | 12 | yes | 0/0 | 2 | 0.0003 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.005 | 4 | yes | 0/4 | 0 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 12 | yes | 0/0 | 2 | 0.0007 | 0.1 | 🟢 |
| Multigraphs | 2% | STRidge | 0.5 | 8 | no | 0/0 | 2 | 0.0032 | 0.1 | 🟢 |
| Multigraphs | 2% | SR3 | 0.01 | 4 | no | 0/2 | 1 | 0.9995 | 0.0 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 12 | no | 0/0 | 2 | 0.0024 | 0.1 | 🟢 |
| Multigraphs | 2% | STRidge | 0.5 | 8 | yes | 0/0 | 2 | 0.0032 | 0.1 | 🟢 |
| Multigraphs | 2% | SR3 | 0.01 | 4 | yes | 0/2 | 1 | 0.9995 | 0.0 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 12 | yes | 0/0 | 2 | 0.0024 | 0.1 | 🟢 |
| Multigraphs | 5% | STRidge | 1.0 | 12 | no | 0/0 | 2 | 0.0051 | 0.1 | 🟢 |
| Multigraphs | 5% | SR3 | 0.05 | 8 | no | 0/2 | 1 | 0.9028 | 0.0 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.25 | 8 | no | 0/0 | 2 | 0.0055 | 0.2 | 🟢 |
| Multigraphs | 5% | STRidge | 1.0 | 12 | yes | 0/0 | 2 | 0.0051 | 0.1 | 🟢 |
| Multigraphs | 5% | SR3 | 0.05 | 8 | yes | 0/2 | 1 | 0.9029 | 0.0 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 4 | yes | 3/2 | 2 | 7.7843 | 0.2 | 🔴 |
| MWSINDy | Clean | STRidge | 0.01 | 12 | no | 0/0 | 2 | 0.0007 | 0.2 | 🟢 |
| MWSINDy | Clean | SR3 | 0.05 | 12 | no | 0/0 | 2 | 0.0007 | 0.1 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0007 | 0.3 | 🟢 |
| MWSINDy | Clean | STRidge | 0.01 | 12 | yes | 0/0 | 2 | 0.0007 | 0.2 | 🟢 |
| MWSINDy | Clean | SR3 | 0.05 | 12 | yes | 0/0 | 2 | 0.0007 | 0.1 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0007 | 0.3 | 🟢 |
| MWSINDy | 2% | STRidge | 0.2 | 12 | no | 0/0 | 2 | 0.0024 | 0.2 | 🟢 |
| MWSINDy | 2% | SR3 | 0.05 | 12 | no | 0/0 | 2 | 0.0024 | 0.1 | 🟢 |
| MWSINDy | 2% | MSTLS | 0.25 | 12 | no | 0/0 | 2 | 0.0024 | 0.3 | 🟢 |
| MWSINDy | 2% | STRidge | 0.5 | 12 | yes | 0/0 | 2 | 0.0024 | 0.2 | 🟢 |
| MWSINDy | 2% | SR3 | 0.05 | 12 | yes | 0/0 | 2 | 0.0024 | 0.1 | 🟢 |
| MWSINDy | 2% | MSTLS | 0.25 | 12 | yes | 0/0 | 2 | 0.0024 | 0.3 | 🟢 |
| MWSINDy | 5% | STRidge | 0.1 | 8 | no | 0/0 | 2 | 0.0055 | 0.2 | 🟢 |
| MWSINDy | 5% | SR3 | 0.05 | 8 | no | 0/0 | 2 | 0.0055 | 0.1 | 🟢 |
| MWSINDy | 5% | MSTLS | 0.25 | 8 | no | 0/0 | 2 | 0.0055 | 0.4 | 🟢 |
| MWSINDy | 5% | STRidge | 0.5 | 8 | yes | 0/0 | 2 | 0.0055 | 0.2 | 🟢 |
| MWSINDy | 5% | SR3 | 0.1 | 8 | yes | 0/0 | 2 | 0.0055 | 0.1 | 🟢 |
| MWSINDy | 5% | MSTLS | 0.25 | 8 | yes | 0/0 | 2 | 0.0055 | 0.4 | 🟢 |

---
*Source: `sweep_20260623_133611_pme_full.csv`*

*Run: 182.0 min (Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz), OMP_NUM_THREADS=not set*
