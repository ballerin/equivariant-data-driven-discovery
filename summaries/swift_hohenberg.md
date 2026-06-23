# Swift-Hohenberg equation

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
max order 4, max degree 3).  fp/fn/E are computed against the same 4-term
ground truth for every library.  Swift-Hohenberg contains a 4th-order biharmonic term $\Delta^2 f = f_{xxxx} + 2f_{xxyy} + f_{yyyy}$ which is the main source of difficulty.  All four libraries can represent every term individually, but the equivariant basis encodes the isotropic combination as a single operator, giving it a 3\times sparsity advantage over coordinate bases for the biharmonic term.

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
| Clean | pointwise | spectral | 🟩 | 0/0 ✓ | 0.0248 | PolyJet | 0/0 ✓ | 0.0252 | 0/0 ✓ | 0.0273 |
| Clean | pointwise | polynomial | 🟩 | 0/0 ✓ | 0.2052 | WSINDy | 0/0 ✓ | 0.2075 | 0/0 ✓ | 0.2077 |
| Clean | weak (no IBP) | spectral | 🟩 | 0/0 ✓ | 0.0252 | PolyJet | 0/0 ✓ | 0.0255 | 0/0 ✓ | 0.0260 |
| Clean | weak (no IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0249 | PDE-FIND | 0/0 ✓ | 0.0253 | 0/0 ✓ | 0.0261 |
| Clean | weak (IBP) | spectral | 🟥 | 0/0 ✓ | 0.1719 | WSINDy | 0/0 ✓ | 0.0214 | 0/0 ✓ | 0.0225 |
| Clean | weak (IBP) | polynomial | 🟥 | 0/0 ✓ | 0.0230 | WSINDy | 0/0 ✓ | 0.0214 | 0/0 ✓ | 0.0225 |
| 2% | pointwise | spectral | — | 0/4 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | 0/2 ✗ | 1.0000 |
| 2% | pointwise | polynomial | — | 0/1 ✗ | 1.0009 | — | 0/4 ✗ | 1.0000 | 0/3 ✗ | 0.9998 |
| 2% | weak (no IBP) | spectral | 🟥 | 0/0 ✓ | 0.0040 | PolyJet | 0/0 ✓ | 0.0032 | 0/0 ✓ | 0.0100 |
| 2% | weak (no IBP) | polynomial | 🟥 | 0/0 ✓ | 0.0050 | WSINDy | 0/0 ✓ | 0.0044 | 0/0 ✓ | 0.0110 |
| 2% | weak (IBP) | spectral | 🟩 | 0/0 ✓ | 0.0023 | PolyJet | 0/0 ✓ | 0.0045 | 0/0 ✓ | 0.0083 |
| 2% | weak (IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0023 | PolyJet | 0/0 ✓ | 0.0045 | 0/0 ✓ | 0.0083 |
| 5% | pointwise | spectral | — | 0/6 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | 0/3 ✗ | 1.0000 |
| 5% | pointwise | polynomial | — | 0/4 ✗ | 1.0000 | — | 0/4 ✗ | 1.0000 | 0/5 ✗ | 1.0000 |
| 5% | weak (no IBP) | spectral | 🟩 | 0/0 ✓ | 0.0521 | PDE-FIND | 0/0 ✓ | 0.0584 | 0/0 ✓ | 0.0472 |
| 5% | weak (no IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0497 | PDE-FIND | 0/0 ✓ | 0.0556 | 0/0 ✓ | 0.0448 |
| 5% | weak (IBP) | spectral | 🟩 | 0/0 ✓ | 0.0536 | PDE-FIND | 0/0 ✓ | 0.0573 | 0/0 ✓ | 0.0487 |
| 5% | weak (IBP) | polynomial | 🟩 | 0/0 ✓ | 0.0536 | PDE-FIND | 0/0 ✓ | 0.0573 | 0/0 ✓ | 0.0487 |

## Pointwise formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol/α | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.005 | 0/0 | 7 | 0.0252 | 0.5 | 🟢 |
| PDE-FIND | Clean | SR3 | 0.005 | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 0/0 | 7 | 0.0274 | 1.1 | 🟢 |
| PDE-FIND | 2% | STRidge | 0.05 | 0/5 | 2 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.1 | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 1/6 | 2 | 0.9999 | 1.5 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.05 | 0/5 | 2 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.05 | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | 0/5 | 2 | 1.0000 | 1.6 | 🔴 |
| WSINDy | Clean | STRidge | 0.001 | 0/0 | 7 | 0.0252 | 0.3 | 🟢 |
| WSINDy | Clean | SR3 | 0.01 | 0/5 | 2 | 1.1587 | 0.1 | 🔴 |
| WSINDy | Clean | MSTLS | 0.0 | 0/0 | 7 | 0.0274 | 0.7 | 🟢 |
| WSINDy | 2% | STRidge | 0.005 | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 2% | SR3 | 0.01 | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | 12/3 | 8 | 1.0000 | 0.9 | 🔴 |
| WSINDy | 5% | STRidge | 0.05 | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 5% | SR3 | 0.01 | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.5 | 12/1 | 10 | 1.0000 | 0.9 | 🔴 |
| PolyJet | Clean | STRidge | 0.05 | 0/0 | 7 | 0.0252 | 0.9 | 🟢 |
| PolyJet | Clean | SR3 | 0.01 | 0/4 | 3 | 0.9317 | 0.3 | 🔴 |
| PolyJet | Clean | MSTLS | 0.0 | 0/0 | 7 | 0.0274 | 2.4 | 🟢 |
| PolyJet | 2% | STRidge | 0.005 | 0/5 | 2 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 2% | SR3 | 0.05 | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | 1/5 | 3 | 1.0000 | 3.1 | 🔴 |
| PolyJet | 5% | STRidge | 0.1 | 0/7 | 0 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 5% | SR3 | 0.01 | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | 3/3 | 7 | 1.0000 | 3.3 | 🔴 |
| Multigraphs | Clean | STRidge | 0.1 | 0/0 | 4 | 0.0248 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.05 | 0/1 | 3 | 0.9385 | 0.1 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 5/0 | 6 | 0.0276 | 0.3 | 🔴 |
| Multigraphs | 2% | STRidge | 0.005 | 0/4 | 2 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.05 | 0/4 | 2 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 2/5 | 3 | 1.0000 | 0.4 | 🔴 |
| Multigraphs | 5% | STRidge | 0.01 | 0/6 | 1 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | SR3 | 0.1 | 0/6 | 1 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 6/5 | 4 | 1.0000 | 0.4 | 🔴 |
| MWSINDy | Clean | STRidge | 0.0005 | 0/0 | 8 | 0.0273 | 0.4 | 🟢 |
| MWSINDy | Clean | SR3 | 0.0005 | 0/3 | 9 | 0.8403 | 0.2 | 🔴 |
| MWSINDy | Clean | MSTLS | 0.0 | 0/0 | 7 | 0.0273 | 0.7 | 🟢 |
| MWSINDy | 2% | STRidge | 0.0005 | 0/2 | 19 | 1.0000 | 0.4 | 🔴 |
| MWSINDy | 2% | SR3 | 0.001 | 0/2 | 23 | 1.0000 | 0.2 | 🔴 |
| MWSINDy | 2% | MSTLS | 0.0 | 4/5 | 8 | 1.0000 | 0.9 | 🔴 |
| MWSINDy | 5% | STRidge | 0.05 | 0/5 | 4 | 1.0000 | 0.4 | 🔴 |
| MWSINDy | 5% | SR3 | 0.0005 | 0/3 | 20 | 1.0000 | 0.2 | 🔴 |
| MWSINDy | 5% | MSTLS | 0.5 | 4/3 | 10 | 1.0000 | 0.9 | 🔴 |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol/α | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.05 | 0/0 | 7 | 0.2075 | 0.5 | 🟢 |
| PDE-FIND | Clean | SR3 | 0.01 | 0/5 | 2 | 0.9615 | 0.2 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 0/0 | 7 | 0.2110 | 1.6 | 🟢 |
| PDE-FIND | 2% | STRidge | 0.1 | 0/7 | 0 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 2% | SR3 | 0.01 | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 1/5 | 3 | 1.0002 | 1.5 | 🔴 |
| PDE-FIND | 5% | STRidge | 0.1 | 0/7 | 0 | 1.0000 | 0.5 | 🔴 |
| PDE-FIND | 5% | SR3 | 0.005 | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | 1/5 | 3 | 1.0001 | 1.5 | 🔴 |
| WSINDy | Clean | STRidge | 0.05 | 0/0 | 7 | 0.2075 | 0.3 | 🟢 |
| WSINDy | Clean | SR3 | 0.01 | 0/5 | 2 | 1.1561 | 0.1 | 🔴 |
| WSINDy | Clean | MSTLS | 0.0 | 0/0 | 7 | 0.2110 | 0.8 | 🟢 |
| WSINDy | 2% | STRidge | 0.2 | 0/7 | 0 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 2% | SR3 | 0.01 | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | 6/5 | 4 | 0.9992 | 1.0 | 🔴 |
| WSINDy | 5% | STRidge | 0.2 | 0/7 | 0 | 1.0000 | 0.3 | 🔴 |
| WSINDy | 5% | SR3 | 0.01 | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | 6/5 | 4 | 0.9997 | 1.1 | 🔴 |
| PolyJet | Clean | STRidge | 0.1 | 0/0 | 7 | 0.2075 | 0.9 | 🟢 |
| PolyJet | Clean | SR3 | 0.01 | 0/4 | 3 | 0.9317 | 0.3 | 🔴 |
| PolyJet | Clean | MSTLS | 0.0 | 2/0 | 9 | 0.2139 | 2.4 | 🔴 |
| PolyJet | 2% | STRidge | 0.2 | 0/7 | 0 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 2% | SR3 | 0.01 | 0/5 | 2 | 1.0001 | 0.3 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | 8/5 | 10 | 1.0003 | 3.6 | 🔴 |
| PolyJet | 5% | STRidge | 0.1 | 0/7 | 0 | 1.0000 | 0.9 | 🔴 |
| PolyJet | 5% | SR3 | 0.01 | 0/5 | 2 | 1.0001 | 0.3 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | 8/5 | 10 | 1.0002 | 3.4 | 🔴 |
| Multigraphs | Clean | STRidge | 0.2 | 0/0 | 4 | 0.2052 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.05 | 0/1 | 3 | 0.9384 | 0.1 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 21/0 | 12 | 0.2205 | 0.4 | 🔴 |
| Multigraphs | 2% | STRidge | 0.005 | 0/4 | 2 | 1.0001 | 0.1 | 🔴 |
| Multigraphs | 2% | SR3 | 0.01 | 0/1 | 3 | 1.0009 | 0.1 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 19/4 | 11 | 1.0001 | 0.4 | 🔴 |
| Multigraphs | 5% | STRidge | 0.005 | 0/4 | 2 | 1.0001 | 0.1 | 🔴 |
| Multigraphs | 5% | SR3 | 0.1 | 0/4 | 2 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 13/5 | 8 | 1.0001 | 0.5 | 🔴 |
| MWSINDy | Clean | STRidge | 0.005 | 0/0 | 10 | 0.2077 | 0.4 | 🟢 |
| MWSINDy | Clean | SR3 | 0.0005 | 0/0 | 13 | 0.2077 | 0.2 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.0 | 0/0 | 7 | 0.2077 | 0.9 | 🟢 |
| MWSINDy | 2% | STRidge | 0.001 | 0/3 | 11 | 0.9998 | 0.4 | 🔴 |
| MWSINDy | 2% | SR3 | 0.005 | 0/3 | 12 | 0.9998 | 0.2 | 🔴 |
| MWSINDy | 2% | MSTLS | 0.0 | 0/5 | 4 | 1.0000 | 1.0 | 🔴 |
| MWSINDy | 5% | STRidge | 0.05 | 0/5 | 3 | 1.0000 | 0.4 | 🔴 |
| MWSINDy | 5% | SR3 | 0.0005 | 0/5 | 8 | 1.0000 | 0.2 | 🔴 |
| MWSINDy | 5% | MSTLS | 0.0 | 0/5 | 4 | 1.0000 | 1.1 | 🔴 |

## Weak formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol/α | m | IBP | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.2 | 4 | no | 0/0 | 7 | 0.0255 | 0.5 | 🟢 |
| PDE-FIND | Clean | SR3 | 0.1 | 4 | no | 0/5 | 2 | 0.9454 | 0.2 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0259 | 0.8 | 🟢 |
| PDE-FIND | Clean | STRidge | 0.5 | 12 | yes | 0/0 | 7 | 0.0214 | 0.3 | 🟢 |
| PDE-FIND | Clean | SR3 | 0.05 | 4 | yes | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0233 | 0.8 | 🟢 |
| PDE-FIND | 2% | STRidge | 1.0 | 12 | no | 0/0 | 7 | 0.0032 | 0.3 | 🟢 |
| PDE-FIND | 2% | SR3 | 0.05 | 4 | no | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0085 | 0.8 | 🟢 |
| PDE-FIND | 2% | STRidge | 1.0 | 12 | yes | 0/0 | 7 | 0.0045 | 0.3 | 🟢 |
| PDE-FIND | 2% | SR3 | 0.05 | 4 | yes | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0099 | 0.8 | 🟢 |
| PDE-FIND | 5% | STRidge | 0.5 | 12 | no | 0/0 | 7 | 0.0628 | 0.3 | 🟢 |
| PDE-FIND | 5% | SR3 | 0.05 | 4 | no | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0584 | 0.9 | 🟢 |
| PDE-FIND | 5% | STRidge | 0.5 | 12 | yes | 0/0 | 7 | 0.0617 | 0.3 | 🟢 |
| PDE-FIND | 5% | SR3 | 0.05 | 4 | yes | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0573 | 1.0 | 🟢 |
| WSINDy | Clean | STRidge | 0.01 | 4 | no | 0/0 | 7 | 0.0255 | 0.3 | 🟢 |
| WSINDy | Clean | SR3 | 0.05 | 4 | no | 0/3 | 4 | 1.4590 | 0.1 | 🔴 |
| WSINDy | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0259 | 0.5 | 🟢 |
| WSINDy | Clean | STRidge | 0.1 | 12 | yes | 0/0 | 7 | 0.0214 | 0.2 | 🟢 |
| WSINDy | Clean | SR3 | 0.05 | 4 | yes | 0/3 | 4 | 2.2479 | 0.1 | 🔴 |
| WSINDy | Clean | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0233 | 0.5 | 🟢 |
| WSINDy | 2% | STRidge | 0.5 | 12 | no | 0/0 | 7 | 0.0032 | 0.2 | 🟢 |
| WSINDy | 2% | SR3 | 0.2 | 4 | no | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0085 | 0.5 | 🟢 |
| WSINDy | 2% | STRidge | 0.5 | 12 | yes | 0/0 | 7 | 0.0045 | 0.2 | 🟢 |
| WSINDy | 2% | SR3 | 0.1 | 4 | yes | 0/4 | 3 | 1.0001 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0099 | 0.5 | 🟢 |
| WSINDy | 5% | STRidge | 1.0 | 12 | no | 0/0 | 7 | 0.0628 | 0.2 | 🟢 |
| WSINDy | 5% | SR3 | 0.1 | 4 | no | 0/4 | 3 | 1.0001 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0584 | 0.7 | 🟢 |
| WSINDy | 5% | STRidge | 1.0 | 12 | yes | 0/0 | 7 | 0.0617 | 0.2 | 🟢 |
| WSINDy | 5% | SR3 | 0.1 | 4 | yes | 0/3 | 4 | 0.9999 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0573 | 0.7 | 🟢 |
| PolyJet | Clean | STRidge | 0.1 | 4 | no | 0/0 | 7 | 0.0255 | 0.9 | 🟢 |
| PolyJet | Clean | SR3 | 0.1 | 8 | no | 0/4 | 3 | 0.9311 | 0.2 | 🔴 |
| PolyJet | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0259 | 1.8 | 🟢 |
| PolyJet | Clean | STRidge | 0.5 | 12 | yes | 0/0 | 7 | 0.0214 | 0.5 | 🟢 |
| PolyJet | Clean | SR3 | 0.05 | 4 | yes | 0/3 | 4 | 0.9447 | 0.3 | 🔴 |
| PolyJet | Clean | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0233 | 1.9 | 🟢 |
| PolyJet | 2% | STRidge | 0.5 | 12 | no | 0/0 | 7 | 0.0032 | 0.5 | 🟢 |
| PolyJet | 2% | SR3 | 0.2 | 12 | no | 0/4 | 3 | 1.0010 | 0.1 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0085 | 2.2 | 🟢 |
| PolyJet | 2% | STRidge | 0.5 | 12 | yes | 0/0 | 7 | 0.0045 | 0.6 | 🟢 |
| PolyJet | 2% | SR3 | 0.2 | 12 | yes | 0/4 | 3 | 1.0010 | 0.2 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0099 | 2.2 | 🟢 |
| PolyJet | 5% | STRidge | 1.5 | 12 | no | 0/0 | 7 | 0.0628 | 0.6 | 🟢 |
| PolyJet | 5% | SR3 | 0.1 | 4 | no | 0/5 | 2 | 1.0000 | 0.3 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | 4 | no | 3/4 | 6 | 1.0002 | 3.3 | 🔴 |
| PolyJet | 5% | STRidge | 1.5 | 12 | yes | 0/0 | 7 | 0.0617 | 0.6 | 🟢 |
| PolyJet | 5% | SR3 | 0.05 | 12 | yes | 0/6 | 1 | 1.0000 | 0.2 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | 4 | yes | 3/5 | 5 | 1.0000 | 3.4 | 🔴 |
| Multigraphs | Clean | STRidge | 1.0 | 4 | no | 0/0 | 4 | 0.0252 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.1 | 4 | no | 0/1 | 3 | 0.9343 | 0.1 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 4 | 0.0260 | 0.2 | 🟢 |
| Multigraphs | Clean | STRidge | 3.0 | 8 | yes | 0/0 | 4 | 0.1719 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.5 | 8 | yes | 0/1 | 3 | 0.9352 | 0.1 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 12 | yes | 5/0 | 6 | 0.0226 | 0.2 | 🔴 |
| Multigraphs | 2% | STRidge | 1.0 | 12 | no | 0/0 | 4 | 0.0040 | 0.1 | 🟢 |
| Multigraphs | 2% | SR3 | 0.5 | 4 | no | 0/1 | 3 | 1.0002 | 0.1 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 12 | no | 8/0 | 7 | 0.0128 | 0.2 | 🔴 |
| Multigraphs | 2% | STRidge | 1.0 | 12 | yes | 0/0 | 4 | 0.0023 | 0.1 | 🟢 |
| Multigraphs | 2% | SR3 | 2.0 | 12 | yes | 0/1 | 3 | 1.0003 | 0.0 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 12 | yes | 8/0 | 7 | 0.0112 | 0.2 | 🔴 |
| Multigraphs | 5% | STRidge | 3.0 | 12 | no | 0/0 | 4 | 0.0521 | 0.1 | 🟢 |
| Multigraphs | 5% | SR3 | 0.5 | 8 | no | 0/1 | 3 | 1.0011 | 0.1 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 4 | no | 2/5 | 3 | 1.0000 | 0.4 | 🔴 |
| Multigraphs | 5% | STRidge | 3.0 | 12 | yes | 0/0 | 4 | 0.0536 | 0.1 | 🟢 |
| Multigraphs | 5% | SR3 | 0.5 | 8 | yes | 0/1 | 3 | 1.0011 | 0.1 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 4 | yes | 2/5 | 3 | 1.0000 | 0.4 | 🔴 |
| MWSINDy | Clean | STRidge | 0.05 | 12 | no | 0/0 | 9 | 0.0260 | 0.2 | 🟢 |
| MWSINDy | Clean | SR3 | 0.005 | 4 | no | 0/0 | 10 | 0.0272 | 0.2 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0260 | 0.5 | 🟢 |
| MWSINDy | Clean | STRidge | 0.05 | 12 | yes | 0/0 | 9 | 0.0225 | 0.2 | 🟢 |
| MWSINDy | Clean | SR3 | 0.05 | 4 | yes | 0/0 | 11 | 19.5859 | 0.2 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0225 | 0.5 | 🟢 |
| MWSINDy | 2% | STRidge | 0.05 | 12 | no | 0/0 | 9 | 0.0100 | 0.2 | 🟢 |
| MWSINDy | 2% | SR3 | 0.5 | 4 | no | 0/3 | 10 | 0.9995 | 0.2 | 🔴 |
| MWSINDy | 2% | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0100 | 0.5 | 🟢 |
| MWSINDy | 2% | STRidge | 0.1 | 12 | yes | 0/0 | 9 | 0.0083 | 0.2 | 🟢 |
| MWSINDy | 2% | SR3 | 0.2 | 4 | yes | 0/0 | 11 | 0.9987 | 0.1 | 🟢 |
| MWSINDy | 2% | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0083 | 0.5 | 🟢 |
| MWSINDy | 5% | STRidge | 0.2 | 12 | no | 0/0 | 9 | 0.0472 | 0.2 | 🟢 |
| MWSINDy | 5% | SR3 | 0.2 | 8 | no | 0/3 | 5 | 1.0002 | 0.1 | 🔴 |
| MWSINDy | 5% | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0472 | 0.7 | 🟢 |
| MWSINDy | 5% | STRidge | 0.2 | 12 | yes | 0/0 | 9 | 0.0487 | 0.2 | 🟢 |
| MWSINDy | 5% | SR3 | 0.05 | 4 | yes | 0/3 | 8 | 0.9998 | 0.2 | 🔴 |
| MWSINDy | 5% | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0487 | 0.7 | 🟢 |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol/α | m | IBP | fp/fn | n_active | E | Time (s) | Exact? |
|---|---|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.5 | 8 | no | 0/0 | 7 | 0.0253 | 0.4 | 🟢 |
| PDE-FIND | Clean | SR3 | 0.1 | 4 | no | 0/5 | 2 | 0.9471 | 0.2 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0260 | 0.8 | 🟢 |
| PDE-FIND | Clean | STRidge | 2.0 | 12 | yes | 0/0 | 7 | 0.0214 | 0.3 | 🟢 |
| PDE-FIND | Clean | SR3 | 0.05 | 4 | yes | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | Clean | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0233 | 0.8 | 🟢 |
| PDE-FIND | 2% | STRidge | 1.0 | 12 | no | 0/0 | 7 | 0.0044 | 0.3 | 🟢 |
| PDE-FIND | 2% | SR3 | 0.05 | 4 | no | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0097 | 1.0 | 🟢 |
| PDE-FIND | 2% | STRidge | 1.0 | 12 | yes | 0/0 | 7 | 0.0045 | 0.4 | 🟢 |
| PDE-FIND | 2% | SR3 | 0.05 | 4 | yes | 0/7 | 0 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 2% | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0099 | 1.0 | 🟢 |
| PDE-FIND | 5% | STRidge | 1.5 | 12 | no | 0/0 | 7 | 0.0599 | 0.3 | 🟢 |
| PDE-FIND | 5% | SR3 | 0.1 | 4 | no | 0/5 | 2 | 1.0001 | 0.2 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0556 | 1.1 | 🟢 |
| PDE-FIND | 5% | STRidge | 1.5 | 12 | yes | 0/0 | 7 | 0.0617 | 0.3 | 🟢 |
| PDE-FIND | 5% | SR3 | 0.1 | 4 | yes | 0/5 | 2 | 1.0000 | 0.2 | 🔴 |
| PDE-FIND | 5% | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0573 | 1.1 | 🟢 |
| WSINDy | Clean | STRidge | 0.1 | 8 | no | 0/0 | 7 | 0.0253 | 0.3 | 🟢 |
| WSINDy | Clean | SR3 | 0.05 | 4 | no | 0/3 | 4 | 1.4120 | 0.1 | 🔴 |
| WSINDy | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0260 | 0.5 | 🟢 |
| WSINDy | Clean | STRidge | 0.1 | 12 | yes | 0/0 | 7 | 0.0214 | 0.2 | 🟢 |
| WSINDy | Clean | SR3 | 0.05 | 4 | yes | 0/3 | 4 | 2.1889 | 0.1 | 🔴 |
| WSINDy | Clean | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0233 | 0.5 | 🟢 |
| WSINDy | 2% | STRidge | 0.5 | 12 | no | 0/0 | 7 | 0.0044 | 0.2 | 🟢 |
| WSINDy | 2% | SR3 | 0.1 | 4 | no | 0/4 | 3 | 1.0002 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0097 | 0.5 | 🟢 |
| WSINDy | 2% | STRidge | 0.5 | 12 | yes | 0/0 | 7 | 0.0045 | 0.2 | 🟢 |
| WSINDy | 2% | SR3 | 0.1 | 4 | yes | 0/4 | 3 | 1.0001 | 0.1 | 🔴 |
| WSINDy | 2% | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0099 | 0.5 | 🟢 |
| WSINDy | 5% | STRidge | 1.0 | 12 | no | 0/0 | 7 | 0.0599 | 0.2 | 🟢 |
| WSINDy | 5% | SR3 | 0.2 | 4 | no | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0556 | 0.7 | 🟢 |
| WSINDy | 5% | STRidge | 1.0 | 12 | yes | 0/0 | 7 | 0.0617 | 0.2 | 🟢 |
| WSINDy | 5% | SR3 | 0.1 | 4 | yes | 0/4 | 3 | 1.0000 | 0.1 | 🔴 |
| WSINDy | 5% | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0573 | 0.7 | 🟢 |
| PolyJet | Clean | STRidge | 0.5 | 8 | no | 0/0 | 7 | 0.0253 | 0.7 | 🟢 |
| PolyJet | Clean | SR3 | 0.1 | 8 | no | 0/4 | 3 | 0.9311 | 0.2 | 🔴 |
| PolyJet | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0260 | 1.8 | 🟢 |
| PolyJet | Clean | STRidge | 0.5 | 12 | yes | 0/0 | 7 | 0.0214 | 0.6 | 🟢 |
| PolyJet | Clean | SR3 | 0.05 | 4 | yes | 0/3 | 4 | 0.9447 | 0.3 | 🔴 |
| PolyJet | Clean | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0233 | 1.9 | 🟢 |
| PolyJet | 2% | STRidge | 4.0 | 12 | no | 0/0 | 7 | 0.0044 | 0.6 | 🟢 |
| PolyJet | 2% | SR3 | 0.2 | 12 | no | 0/4 | 3 | 1.0009 | 0.1 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | 12 | no | 3/0 | 10 | 0.0051 | 2.7 | 🔴 |
| PolyJet | 2% | STRidge | 4.0 | 12 | yes | 0/0 | 7 | 0.0045 | 0.6 | 🟢 |
| PolyJet | 2% | SR3 | 0.2 | 12 | yes | 0/4 | 3 | 1.0009 | 0.2 | 🔴 |
| PolyJet | 2% | MSTLS | 0.0 | 12 | yes | 3/0 | 10 | 0.0051 | 2.8 | 🔴 |
| PolyJet | 5% | STRidge | 2.0 | 12 | no | 0/0 | 7 | 0.0599 | 0.5 | 🟢 |
| PolyJet | 5% | SR3 | 0.1 | 12 | no | 0/3 | 4 | 1.0018 | 0.1 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | 12 | no | 2/0 | 9 | 0.0720 | 3.2 | 🔴 |
| PolyJet | 5% | STRidge | 2.0 | 12 | yes | 0/0 | 7 | 0.0617 | 0.6 | 🟢 |
| PolyJet | 5% | SR3 | 0.1 | 12 | yes | 0/3 | 4 | 1.0018 | 0.2 | 🔴 |
| PolyJet | 5% | MSTLS | 0.0 | 12 | yes | 2/0 | 9 | 0.0744 | 3.4 | 🔴 |
| Multigraphs | Clean | STRidge | 1.5 | 8 | no | 0/0 | 4 | 0.0249 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.1 | 4 | no | 0/1 | 3 | 0.9343 | 0.1 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 4 | 0.0261 | 0.3 | 🟢 |
| Multigraphs | Clean | STRidge | 5.0 | 12 | yes | 0/0 | 4 | 0.0230 | 0.1 | 🟢 |
| Multigraphs | Clean | SR3 | 0.5 | 8 | yes | 0/1 | 3 | 0.9352 | 0.1 | 🔴 |
| Multigraphs | Clean | MSTLS | 0.0 | 12 | yes | 5/0 | 6 | 0.0225 | 0.2 | 🔴 |
| Multigraphs | 2% | STRidge | 2.0 | 12 | no | 0/0 | 4 | 0.0050 | 0.1 | 🟢 |
| Multigraphs | 2% | SR3 | 1.5 | 4 | no | 0/1 | 3 | 1.0001 | 0.1 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 12 | no | 5/0 | 6 | 0.0145 | 0.2 | 🔴 |
| Multigraphs | 2% | STRidge | 3.0 | 12 | yes | 0/0 | 4 | 0.0023 | 0.1 | 🟢 |
| Multigraphs | 2% | SR3 | 0.05 | 4 | yes | 0/1 | 3 | 0.9986 | 0.1 | 🔴 |
| Multigraphs | 2% | MSTLS | 0.0 | 12 | yes | 5/0 | 6 | 0.0075 | 0.2 | 🔴 |
| Multigraphs | 5% | STRidge | 2.0 | 12 | no | 0/0 | 4 | 0.0497 | 0.1 | 🟢 |
| Multigraphs | 5% | SR3 | 5.0 | 12 | no | 0/1 | 3 | 1.0001 | 0.0 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 12 | no | 18/0 | 11 | 0.1124 | 0.2 | 🔴 |
| Multigraphs | 5% | STRidge | 2.0 | 12 | yes | 0/0 | 4 | 0.0536 | 0.1 | 🟢 |
| Multigraphs | 5% | SR3 | 0.05 | 4 | yes | 0/1 | 3 | 1.0000 | 0.1 | 🔴 |
| Multigraphs | 5% | MSTLS | 0.0 | 12 | yes | 18/0 | 11 | 0.1175 | 0.2 | 🔴 |
| MWSINDy | Clean | STRidge | 0.05 | 12 | no | 0/0 | 9 | 0.0261 | 0.2 | 🟢 |
| MWSINDy | Clean | SR3 | 0.005 | 4 | no | 0/0 | 11 | 0.1434 | 0.2 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0261 | 0.5 | 🟢 |
| MWSINDy | Clean | STRidge | 0.05 | 12 | yes | 0/0 | 9 | 0.0225 | 0.2 | 🟢 |
| MWSINDy | Clean | SR3 | 0.05 | 4 | yes | 0/0 | 11 | 19.5859 | 0.2 | 🟢 |
| MWSINDy | Clean | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0225 | 0.5 | 🟢 |
| MWSINDy | 2% | STRidge | 0.1 | 12 | no | 0/0 | 9 | 0.0110 | 0.2 | 🟢 |
| MWSINDy | 2% | SR3 | 0.005 | 4 | no | 0/3 | 9 | 0.9994 | 0.2 | 🔴 |
| MWSINDy | 2% | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0110 | 0.6 | 🟢 |
| MWSINDy | 2% | STRidge | 0.1 | 12 | yes | 0/0 | 9 | 0.0083 | 0.2 | 🟢 |
| MWSINDy | 2% | SR3 | 0.2 | 4 | yes | 0/0 | 11 | 0.9987 | 0.1 | 🟢 |
| MWSINDy | 2% | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0083 | 0.5 | 🟢 |
| MWSINDy | 5% | STRidge | 0.2 | 12 | no | 0/0 | 9 | 0.0448 | 0.2 | 🟢 |
| MWSINDy | 5% | SR3 | 0.2 | 4 | no | 0/3 | 10 | 0.9996 | 0.2 | 🔴 |
| MWSINDy | 5% | MSTLS | 0.0 | 12 | no | 0/0 | 7 | 0.0448 | 0.7 | 🟢 |
| MWSINDy | 5% | STRidge | 0.2 | 12 | yes | 0/0 | 9 | 0.0487 | 0.2 | 🟢 |
| MWSINDy | 5% | SR3 | 0.05 | 4 | yes | 0/3 | 9 | 0.9998 | 0.2 | 🔴 |
| MWSINDy | 5% | MSTLS | 0.0 | 12 | yes | 0/0 | 7 | 0.0487 | 0.7 | 🟢 |

---
*Source: `sweep_20260622_202544_sh_full_v2.csv`*

*Run: 194.2 min (Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz), OMP_NUM_THREADS=8*
