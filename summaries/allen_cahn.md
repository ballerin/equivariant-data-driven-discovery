# Allen–Cahn equation

$$f_t = \varepsilon^2 \Delta f + f - f^3, \qquad \varepsilon = 0.15$$

Equivariant library:

| Operator | True coefficient |
|---|---|
| $f$ | $+1.0$ |
| $\Delta f$ | $+\varepsilon^2 = 0.0225$ |
| $f^3$ | $-1.0$ |

Standard libraries:
| Operator | True coefficient | PDE-FIND | WSINDy |
|---|---|---|---|
| $f$ | $+1.0$ | ✅ | ✅ |
| $\partial_{xx} f$ | $+\varepsilon^2 = 0.0225$ | ✅ | ✅ |
| $\partial_{yy} f$ | $+\varepsilon^2 = 0.0225$ | ✅ | ✅ |
| $f^3$ | $-1.0$ | ✅ | ✅ |

Allen–Cahn remains the easy benchmark: it is second-order, mildly nonlinear,
and **all its terms are representable in the coordinate-basis standard library**
(no $\|\nabla f\|^2$ or other product-of-derivatives terms).

## Spectral derivatives

| Library | Noise | Algorithm | d_tol | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.05 | 0/0 | 4 | 0.00778 | ✅ |
| PDE-FIND | 2% | STRidge | 0.20 | 0/2 | 2 | 0.771 | ❌ |
| PDE-FIND | 5% | STRidge | 0.20 | 0/2 | 2 | 0.957 | ❌ |
| Equivariant | Clean | STRidge | 0.05 | 0/0 | 3 | 0.00766 | ✅ |
| Equivariant | 2% | STRidge | 0.20 | 0/1 | 2 | 0.764 | ❌ |
| Equivariant | 5% | STRidge | 0.20 | 0/1 | 2 | 0.954 | ❌ |

At 0% noise, both libraries achieve exact recovery — the standard library
needs 4 columns (f, f_xx, f_yy, f³) to represent the 3 equivariant terms
(f, Δf=f_xx+f_yy, f³).  At 2%+ noise, both libraries lose terms because
the small Laplacian coefficient ε²=0.0225 is drowned by noise-amplified
spectral derivatives.  The |k|⁴ amplification makes spectral derivatives
fragile for any library structure on Allen–Cahn.

## Polynomial stencils

| Library | Noise | Algorithm | d_tol | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.05 | 0/0 | 4 | 0.0109 | ✅ |
| PDE-FIND | 2% | STRidge | 0.05 | 0/0 | 4 | 0.0285 | ✅ |
| PDE-FIND | 5% | SR3 | 0.05 | 0/0 | 4 | 0.272 | ✅ |
| Equivariant | Clean | STRidge | 0.05 | 0/0 | 3 | 0.0109 | ✅ |
| Equivariant | 2% | STRidge | 0.05 | 0/0 | 3 | 0.0535 | ✅ |
| Equivariant | 5% | STRidge | 0.10 | 0/0 | 3 | 0.339 | ✅ |

With polynomial stencils, **both libraries achieve exact recovery at all
three noise levels**.  The natural smoothing of the polynomial fit tames
the noise amplification that kills spectral derivatives.  The only
difference is that the standard library uses 4 columns (f, f_xx, f_yy, f³)
while the equivariant library uses 3 (f, Δf, f³).

## Weak formulation

### Spectral backend

| Library | Noise | Algorithm | d_tol | m | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.20 | 12 | 0/0 | 4 | 0.00463 | ✅ |
| PDE-FIND | 2% | STRidge | 0.10 | 12 | 0/0 | 4 | 0.00454 | ✅ |
| PDE-FIND | 5% | STRidge | 0.20 | 8 | 0/0 | 4 | 0.00656 | ✅ |
| WSINDy | Clean | STRidge | 0.10 | 12 | 0/0 | 3 | 0.00473 | ✅ |
| WSINDy | 2% | STRidge | 0.10 | 12 | 0/0 | 3 | 0.00532 | ✅ |
| WSINDy | 5% | STRidge | 5.0 | 8 | 0/0 | 3 | 0.00868 | ✅ |
| WSINDy IBP | Clean | STRidge | 5.0 | 12 | 0/0 | 3 | 0.00589 | ✅ |
| WSINDy IBP | 2% | STRidge | 0.10 | 12 | 0/0 | 3 | 0.00361 | ✅ |
| WSINDy IBP | 5% | STRidge | 5.0 | 8 | 0/0 | 3 | 0.00245 | ✅ |

### Polynomial backend

| Library | Noise | Algorithm | d_tol | m | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.05 | 12 | 0/0 | 4 | 0.00463 | ✅ |
| PDE-FIND | 2% | STRidge | 0.10 | 12 | 0/0 | 4 | 0.00454 | ✅ |
| PDE-FIND | 5% | STRidge | 0.20 | 12 | 0/0 | 4 | 0.00604 | ✅ |
| WSINDy | Clean | STRidge | 5.0 | 12 | 0/0 | 3 | 0.00473 | ✅ |
| WSINDy | 2% | STRidge | 0.05 | 8 | 0/0 | 3 | 0.00511 | ✅ |
| WSINDy | 5% | STRidge | 5.0 | 8 | 0/0 | 3 | 0.00830 | ✅ |
| WSINDy IBP | Clean | STRidge | 0.10 | 12 | 0/0 | 3 | 0.00589 | ✅ |
| WSINDy IBP | 2% | STRidge | 0.20 | 12 | 0/0 | 3 | 0.00361 | ✅ |
| WSINDy IBP | 5% | STRidge | 5.0 | 8 | 0/0 | 3 | 0.00245 | ✅ |

**Row mapping:** PDE-FIND = standard library + WSINDy direct (the baseline
weak-form method from Messenger & Bortz).  WSINDy = equivariant library +
WSINDy direct (our method, no IBP).  WSINDy IBP = equivariant library +
WSINDy with integration-by-parts on Δf and Δ²f.

**Key result:** all methods and both libraries achieve exact recovery at
all noise levels with both backends via weak formulation.  The weak-form
convolution with smooth test functions effectively denoises even the
standard library.  IBP helps most at 5% noise where it halves E compared
to WSINDy direct (0.00245 vs 0.00868/0.00830).  The only structural
difference is that the standard library uses 4 columns (n_active=4) while
the equivariant library uses 3 (n_active=3).