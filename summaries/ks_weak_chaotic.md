# Kuramoto–Sivashinsky equation — Weak-chaotic regime ($L=8\pi$)

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

KS contains $\|\nabla f\|^2 = f_x^2 + f_y^2$, which neither standard library can express as a sparse column. The weak-chaotic regime ($L=8\pi$) has ~44 unstable modes — sufficient chaotic mixing for good recovery.

## Standard formulation
### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| Multigraphs | Clean | STRidge | 0.5 | 0/0 | 3 | 0.0413 | ✅ |
| Multigraphs | 2% | SR3 | 1 | 0/3 | 0 | 1.000 | ❌ |
| Multigraphs | 5% | SR3 | 1 | 0/3 | 0 | 1.000 | ❌ |
| PDE-FIND | Clean | SR3 | 1 | 0/5 | 0 | 1.000 | ❌ |
| PDE-FIND | 2% | SR3 | 1 | 0/5 | 0 | 1.000 | ❌ |
| PDE-FIND | 5% | SR3 | 1 | 0/5 | 0 | 1.000 | ❌ |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| Multigraphs | Clean | STRidge | 0.5 | 0/0 | 3 | 0.0119 | ✅ |
| Multigraphs | 2% | SR3 | 1 | 0/3 | 0 | 1.000 | ❌ |
| Multigraphs | 5% | SR3 | 1 | 0/3 | 0 | 1.000 | ❌ |
| PDE-FIND | Clean | SR3 | 1 | 0/5 | 0 | 1.000 | ❌ |
| PDE-FIND | 2% | SR3 | 1 | 0/5 | 0 | 1.000 | ❌ |
| PDE-FIND | 5% | SR3 | 1 | 0/5 | 0 | 1.000 | ❌ |

## Weak formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol | m | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND (no IBP) | Clean | SR3 | 1 | 4 | 3/5 | 3 | 2.153 | ❌ |
| PDE-FIND (no IBP) | 2% | SR3 | 1 | 4 | 3/5 | 3 | 1.009 | ❌ |
| PDE-FIND (no IBP) | 5% | SR3 | 0.2 | 4 | 3/5 | 3 | 1.002 | ❌ |
| PDE-FIND (yes IBP) | Clean | SR3 | 1 | 4 | 1/5 | 1 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | 2% | SR3 | 0.5 | 4 | 1/5 | 1 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | 5% | SR3 | 1 | 8 | 2/5 | 2 | 0.9982 | ❌ |
| Multigraphs (no IBP) | Clean | STRidge | 1 | 4 | 4/0 | 7 | 0.0101 | ❌ |
| Multigraphs (no IBP) | 2% | SR3 | 0.5 | 4 | 6/2 | 7 | 1.000 | ❌ |
| Multigraphs (no IBP) | 5% | SR3 | 1 | 4 | 3/3 | 3 | 1.000 | ❌ |
| Multigraphs (yes IBP) | Clean | SR3 | 1 | 4 | 5/0 | 8 | 0.258 | ❌ |
| Multigraphs (yes IBP) | 2% | SR3 | 0.5 | 4 | 6/2 | 7 | 1.000 | ❌ |
| Multigraphs (yes IBP) | 5% | SR3 | 1 | 4 | 3/3 | 3 | 1.000 | ❌ |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol | m | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND (no IBP) | Clean | SR3 | 1 | 4 | 3/5 | 3 | 2.170 | ❌ |
| PDE-FIND (no IBP) | 2% | SR3 | 0.5 | 4 | 3/5 | 3 | 1.037 | ❌ |
| PDE-FIND (no IBP) | 5% | SR3 | 1 | 4 | 3/5 | 3 | 1.007 | ❌ |
| PDE-FIND (yes IBP) | Clean | SR3 | 1 | 4 | 1/5 | 1 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | 2% | SR3 | 1 | 4 | 1/5 | 1 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | 5% | SR3 | 1 | 4 | 1/5 | 1 | 1.000 | ❌ |
| Multigraphs (no IBP) | Clean | SR3 | 0.5 | 4 | 5/0 | 8 | 0.2307 | ❌ |
| Multigraphs (no IBP) | 2% | SR3 | 0.5 | 4 | 8/2 | 9 | 1.001 | ❌ |
| Multigraphs (no IBP) | 5% | SR3 | 1 | 4 | 6/2 | 7 | 1.000 | ❌ |
| Multigraphs (yes IBP) | Clean | SR3 | 1 | 4 | 5/0 | 8 | 0.2546 | ❌ |
| Multigraphs (yes IBP) | 2% | SR3 | 0.5 | 4 | 8/2 | 9 | 1.001 | ❌ |
| Multigraphs (yes IBP) | 5% | SR3 | 1 | 4 | 6/2 | 7 | 1.000 | ❌ |
