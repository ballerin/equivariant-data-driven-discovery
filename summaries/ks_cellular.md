# Kuramoto–Sivashinsky equation — Cellular regime ($L=4\pi$)

The equation under consideration is

$$f_t = -\Delta f - \Delta^2 f - \tfrac{1}{2}\|\nabla f\|^2$$

which can be written in terms of $O(2)$-equivariant terms

| Operator | True coefficient |
|---|---|
| $\Delta f$ | $-1.0$ |
| $\Delta^2 f$ | $-1.0$ |
| $\|\nabla f\|^2$ | $-0.5$ |

or standard terms

| Operator | True coefficient | PDE-FIND |
|---|---|---|
| $\partial_{xx} f$ | $-1.0$ | ✅ |
| $\partial_{yy} f$ | $-1.0$ | ✅ |
| $\partial_{xxxx} f$ | $-1.0$ | ✅ |
| $\partial_{xxyy} f$ | $-2.0$ | ✅ |
| $\partial_{yyyy} f$ | $-1.0$ | ✅ |
| $f_x^2$ | $-0.5$ | ❌ |
| $f_y^2$ | $-0.5$ | ❌ |

KS contains $\|\nabla f\|^2 = f_x^2 + f_y^2$, which neither standard library can express as a sparse column. The cellular regime ($L=4\pi$) is quasi-periodic with only ~8 unstable modes — the field lacks spatiotemporal diversity, collapsing the design matrix.

## Standard formulation
### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| Multigraphs | Clean | STRidge | 0.2 | 0/0 | 3 | 0.00395 | ✅ |
| Multigraphs | 2% | SR3 | 1 | 0/3 | 0 | 1.000 | ❌ |
| Multigraphs | 5% | SR3 | 1 | 0/3 | 0 | 1.000 | ❌ |
| PDE-FIND | Clean | SR3 | 1 | 0/5 | 0 | 1.000 | ❌ |
| PDE-FIND | 2% | SR3 | 1 | 0/5 | 0 | 1.000 | ❌ |
| PDE-FIND | 5% | SR3 | 1 | 0/5 | 0 | 1.000 | ❌ |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| Multigraphs | Clean | STRidge | 0.2 | 0/0 | 3 | 0.00395 | ✅ |
| Multigraphs | 2% | SR3 | 1 | 0/2 | 1 | 1.000 | ❌ |
| Multigraphs | 5% | SR3 | 0.5 | 0/3 | 0 | 1.000 | ❌ |
| PDE-FIND | Clean | SR3 | 1 | 0/5 | 0 | 1.000 | ❌ |
| PDE-FIND | 2% | SR3 | 1 | 0/5 | 0 | 1.000 | ❌ |
| PDE-FIND | 5% | SR3 | 1 | 0/5 | 0 | 1.000 | ❌ |

## Weak formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol | m | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND (no IBP) | Clean | SR3 | 1 | 4 | 1/5 | 1 | 1.437 | ❌ |
| PDE-FIND (no IBP) | 2% | SR3 | 1 | 4 | 1/5 | 1 | 1.000 | ❌ |
| PDE-FIND (no IBP) | 5% | SR3 | 1 | 4 | 3/5 | 3 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | Clean | SR3 | 0.5 | 4 | 0/5 | 0 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | 2% | SR3 | 0.5 | 4 | 1/5 | 1 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | 5% | SR3 | 0.5 | 4 | 1/5 | 1 | 1.000 | ❌ |
| Multigraphs (no IBP) | Clean | STRidge | 1 | 4 | 2/0 | 5 | 0.0574 | ❌ |
| Multigraphs (no IBP) | 2% | SR3 | 1 | 4 | 3/2 | 4 | 1.000 | ❌ |
| Multigraphs (no IBP) | 5% | SR3 | 1 | 4 | 2/2 | 3 | 1.000 | ❌ |
| Multigraphs (yes IBP) | Clean | SR3 | 1 | 4 | 0/3 | 0 | 1.000 | ❌ |
| Multigraphs (yes IBP) | 2% | SR3 | 1 | 4 | 3/2 | 4 | 1.000 | ❌ |
| Multigraphs (yes IBP) | 5% | SR3 | 0.5 | 4 | 5/2 | 6 | 1.000 | ❌ |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol | m | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND (no IBP) | Clean | SR3 | 1 | 4 | 1/5 | 1 | 1.437 | ❌ |
| PDE-FIND (no IBP) | 2% | SR3 | 1 | 4 | 1/5 | 1 | 1.001 | ❌ |
| PDE-FIND (no IBP) | 5% | SR3 | 1 | 4 | 1/5 | 1 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | Clean | SR3 | 0.5 | 4 | 0/5 | 0 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | 2% | SR3 | 0.2 | 4 | 1/5 | 1 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | 5% | SR3 | 0.2 | 4 | 1/5 | 1 | 1.000 | ❌ |
| Multigraphs (no IBP) | Clean | STRidge | 1 | 4 | 2/0 | 5 | 0.0574 | ❌ |
| Multigraphs (no IBP) | 2% | SR3 | 1 | 4 | 2/2 | 3 | 1.000 | ❌ |
| Multigraphs (no IBP) | 5% | SR3 | 1 | 8 | 6/1 | 8 | 1.000 | ❌ |
| Multigraphs (yes IBP) | Clean | SR3 | 1 | 4 | 0/3 | 0 | 1.000 | ❌ |
| Multigraphs (yes IBP) | 2% | SR3 | 1 | 4 | 2/2 | 3 | 1.000 | ❌ |
| Multigraphs (yes IBP) | 5% | SR3 | 1 | 4 | 2/2 | 3 | 1.000 | ❌ |
