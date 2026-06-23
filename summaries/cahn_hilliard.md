# Cahn–Hilliard equation

The equation under consideration is

$$f_t = -\Delta(\Delta f - f + f^3) = -\varepsilon^2\Delta^2 f + \Delta f - 3f^2\Delta f - 6f\|\nabla f\|^2, \qquad \varepsilon = 0.2$$

which can be written in terms of $O(2)$-equivariant terms

| Operator | True coefficient |
|---|---|
| $\Delta f$ | $-1.0$ |
| $f^2\cdot\Delta f$ | $+3.0$ |
| $f\cdot\|\nabla f\|^2$ | $+6.0$ |
| $\Delta^2 f$ | $-\varepsilon^2 = -0.04$ |

or standard terms

| Operator | True coefficient | PDE-FIND |
|---|---|---|
| $\partial_{xx} f$ | $-1.0$ | ✅ |
| $\partial_{yy} f$ | $-1.0$ | ✅ |
| $f^2 \cdot \partial_{xx} f$ | $+3.0$ | ✅ |
| $f^2 \cdot \partial_{yy} f$ | $+3.0$ | ✅ |
| $f \cdot f_x^2$ | $+6.0$ | ❌ |
| $f \cdot f_y^2$ | $+6.0$ | ❌ |
| $\partial_{xxxx} f$ | $-\varepsilon^2$ | ✅ |
| $\partial_{xxyy} f$ | $-2\varepsilon^2$ | ✅ |
| $\partial_{yyyy} f$ | $-\varepsilon^2$ | ✅ |

CH contains $f\|\nabla f\|^2 = f(f_x^2 + f_y^2)$, a product of two derivatives multiplied by $f$ — this **does not exist** as an isolated column in the standard library.

## Standard formulation
### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| Multigraphs | Clean | STRidge | 0.5 | 0/0 | 4 | 0.1272 | ✅ |
| Multigraphs | 2% | STRidge | 1 | 0/4 | 0 | 1.000 | ❌ |
| Multigraphs | 5% | STRidge | 1 | 0/4 | 0 | 1.000 | ❌ |
| PDE-FIND | Clean | STRidge | 1 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND | 2% | STRidge | 1 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND | 5% | STRidge | 1 | 0/7 | 0 | 1.000 | ❌ |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| Multigraphs | Clean | STRidge | 0.2 | 0/4 | 0 | 1.000 | ❌ |
| Multigraphs | 2% | STRidge | 1 | 0/4 | 0 | 1.000 | ❌ |
| Multigraphs | 5% | STRidge | 1 | 0/4 | 0 | 1.000 | ❌ |
| PDE-FIND | Clean | STRidge | 1 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND | 2% | STRidge | 1 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND | 5% | STRidge | 1 | 0/7 | 0 | 1.000 | ❌ |

## Weak formulation

### Spectral derivatives backend

| Library | Noise | Algorithm | d_tol | m | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND (no IBP) | Clean | STRidge | 1 | 4 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND (no IBP) | 2% | STRidge | 0.5 | 4 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND (no IBP) | 5% | STRidge | 0.5 | 4 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | Clean | STRidge | 0.5 | 4 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | 2% | STRidge | 0.5 | 4 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | 5% | STRidge | 0.5 | 4 | 0/7 | 0 | 1.000 | ❌ |
| Multigraphs (no IBP) | Clean | STRidge | 1 | 4 | 7/0 | 11 | 0.1308 | ❌ |
| Multigraphs (no IBP) | 2% | STRidge | 1 | 4 | 0/4 | 0 | 1.000 | ❌ |
| Multigraphs (no IBP) | 5% | STRidge | 0.5 | 4 | 0/4 | 0 | 1.000 | ❌ |
| Multigraphs (yes IBP) | Clean | SR3 | 0.5 | 4 | 1/2 | 3 | 0.9945 | ❌ |
| Multigraphs (yes IBP) | 2% | STRidge | 1 | 4 | 0/4 | 0 | 1.000 | ❌ |
| Multigraphs (yes IBP) | 5% | STRidge | 0.5 | 4 | 0/4 | 0 | 1.000 | ❌ |

### Polynomial stencils backend

| Library | Noise | Algorithm | d_tol | m | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND (no IBP) | Clean | STRidge | 1 | 4 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND (no IBP) | 2% | STRidge | 1 | 4 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND (no IBP) | 5% | STRidge | 1 | 4 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | Clean | STRidge | 0.5 | 4 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | 2% | STRidge | 0.5 | 4 | 0/7 | 0 | 1.000 | ❌ |
| PDE-FIND (yes IBP) | 5% | STRidge | 0.5 | 4 | 0/7 | 0 | 1.000 | ❌ |
| Multigraphs (no IBP) | Clean | SR3 | 1 | 4 | 0/4 | 0 | 1.000 | ❌ |
| Multigraphs (no IBP) | 2% | STRidge | 1 | 4 | 0/4 | 0 | 1.000 | ❌ |
| Multigraphs (no IBP) | 5% | SR3 | 1 | 4 | 0/4 | 0 | 1.000 | ❌ |
| Multigraphs (yes IBP) | Clean | SR3 | 1 | 4 | 0/4 | 0 | 1.000 | ❌ |
| Multigraphs (yes IBP) | 2% | STRidge | 1 | 4 | 0/4 | 0 | 1.000 | ❌ |
| Multigraphs (yes IBP) | 5% | STRidge | 1 | 4 | 0/4 | 0 | 1.000 | ❌ |
