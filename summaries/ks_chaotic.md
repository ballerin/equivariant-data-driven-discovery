# Kuramoto–Sivashinsky equation — Chaotic regime ($L=16\pi$)

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

KS contains $\|\nabla f\|^2 = f_x^2 + f_y^2$, which neither standard library can express as a sparse column. The chaotic regime ($L=16\pi$) has ~192 unstable modes — maximal spatiotemporal diversity.

## Standard formulation
### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| Multigraphs | Clean | STRidge | 0.5 | 0/0 | 3 | 0.00549 | ✅ |
| Multigraphs | 2% | SR3 | 0.5 | 1/3 | 1 | 1.000 | ❌ |
| Multigraphs | 5% | SR3 | 0.5 | 1/3 | 1 | 1.000 | ❌ |
| PDE-FIND | Clean | SR3 | 0.5 | 1/5 | 1 | 1.350 | ❌ |
| PDE-FIND | 2% | SR3 | 0.5 | 1/5 | 1 | 1.000 | ❌ |
| PDE-FIND | 5% | SR3 | 0.5 | 1/5 | 1 | 1.000 | ❌ |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| Multigraphs | Clean | SR3 | 0.5 | 0/2 | 1 | 0.7678 | ❌ |
| Multigraphs | 2% | SR3 | 1 | 0/3 | 0 | 1.000 | ❌ |
| Multigraphs | 5% | SR3 | 1 | 0/3 | 0 | 1.000 | ❌ |
| PDE-FIND | Clean | SR3 | 0.5 | 1/5 | 1 | 1.681 | ❌ |
| PDE-FIND | 2% | SR3 | 0.5 | 1/5 | 1 | 1.038 | ❌ |
| PDE-FIND | 5% | SR3 | 0.5 | 1/5 | 1 | 1.007 | ❌ |

## Weak formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol | m | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND (no IBP) | Clean | SR3 | 0.5 | 4 | 3/5 | 3 | 3.917 | ❌ |
| PDE-FIND (no IBP) | 2% | SR3 | 1 | 4 | 3/5 | 3 | 2.183 | ❌ |
| PDE-FIND (no IBP) | 5% | SR3 | 1 | 4 | 2/5 | 2 | 1.318 | ❌ |
| PDE-FIND (yes IBP) | Clean | SR3 | 1 | 4 | 1/5 | 1 | 0.9971 | ❌ |
| PDE-FIND (yes IBP) | 2% | SR3 | 0.2 | 4 | 1/5 | 1 | 0.9973 | ❌ |
| PDE-FIND (yes IBP) | 5% | SR3 | 0.2 | 4 | 1/5 | 1 | 0.9976 | ❌ |
| Multigraphs (no IBP) | Clean | STRidge | 1 | 4 | 4/0 | 7 | 0.0212 | ❌ |
| Multigraphs (no IBP) | 2% | SR3 | 0.5 | 4 | 5/2 | 6 | 1.000 | ❌ |
| Multigraphs (no IBP) | 5% | SR3 | 0.5 | 4 | 6/2 | 7 | 1.000 | ❌ |
| Multigraphs (yes IBP) | Clean | SR3 | 1 | 4 | 1/0 | 4 | 0.5134 | ❌ |
| Multigraphs (yes IBP) | 2% | SR3 | 0.5 | 4 | 5/2 | 6 | 1.000 | ❌ |
| Multigraphs (yes IBP) | 5% | SR3 | 0.5 | 4 | 6/2 | 7 | 1.000 | ❌ |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol | m | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND (no IBP) | Clean | SR3 | 0.5 | 4 | 3/5 | 3 | 4.262 | ❌ |
| PDE-FIND (no IBP) | 2% | SR3 | 0.5 | 4 | 3/5 | 3 | 3.623 | ❌ |
| PDE-FIND (no IBP) | 5% | SR3 | 1 | 4 | 3/5 | 3 | 2.392 | ❌ |
| PDE-FIND (yes IBP) | Clean | SR3 | 1 | 4 | 1/5 | 1 | 0.9971 | ❌ |
| PDE-FIND (yes IBP) | 2% | SR3 | 1 | 8 | 1/5 | 1 | 0.9293 | ❌ |
| PDE-FIND (yes IBP) | 5% | SR3 | 1 | 8 | 1/5 | 1 | 0.9335 | ❌ |
| Multigraphs (no IBP) | Clean | SR3 | 1 | 4 | 6/0 | 9 | 0.2391 | ❌ |
| Multigraphs (no IBP) | 2% | SR3 | 1 | 8 | 3/2 | 4 | 0.999 | ❌ |
| Multigraphs (no IBP) | 5% | SR3 | 1 | 4 | 4/2 | 5 | 1.012 | ❌ |
| Multigraphs (yes IBP) | Clean | SR3 | 1 | 4 | 1/0 | 4 | 0.4383 | ❌ |
| Multigraphs (yes IBP) | 2% | SR3 | 1 | 8 | 3/2 | 4 | 0.999 | ❌ |
| Multigraphs (yes IBP) | 5% | SR3 | 1 | 4 | 4/2 | 5 | 1.012 | ❌ |
