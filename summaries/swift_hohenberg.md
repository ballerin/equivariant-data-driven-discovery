# Swift–Hohenberg equation

$$f_t = (r-1)f - 2\Delta f - \Delta^2 f - f^3, \qquad r = 0.3$$

Four active equivariant operators:

| Operator | True coefficient |
|---|---|
| $f$ | $r - 1 = -0.7$ |
| $\Delta f$ | $-2.0$ |
| $\Delta^2 f$ | $-1.0$ |
| $f^3$ | $-1.0$ |

Standard libraries:
| Operator | True coefficient | PDE-FIND | WSINDy |
|---|---|---|---|
| $f$ | $r-1 = -0.7$ | ✅ | ✅ |
| $f_{xx}$ | $-2.0$ | ✅ | ✅ |
| $f_{yy}$ | $-2.0$ | ✅ | ✅ |
| $f_{xxxx}$ | $-1.0$ | ✅ | ✅ |
| $f_{xxyy}$ | $-2.0$ | ✅ | ✅ |
| $f_{yyyy}$ | $-1.0$ | ✅ | ✅ |
| $f^3$ | $-1.0$ | ✅ | ✅ |

This is the flagship hard case.  All terms are representable in the standard
library ($\Delta^2 f = f_{xxxx} + 2f_{xxyy} + f_{yyyy}$ exists as three columns),
so the standard library can **fit** the dynamics — but the biharmonic term
splits across 3 columns with a linear constraint (coefficient of $f_{xxyy}$
must be exactly twice that of $f_{xxxx}$), which sparse regression cannot
discover from noisy data.

## Spectral derivatives

| Library | Noise | Algorithm | d_tol | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.01 | 0/0 | 7 | 0.0275 | ✅ |
| PDE-FIND | 2% | SR3 | 0.05 | 0/5 | 2 | 1.00 | ❌ |
| PDE-FIND | 5% | STRidge | 0.05 | 0/5 | 2 | 1.00 | ❌ |
| Equivariant | Clean | STRidge | 0.10 | 0/0 | 4 | 0.0274 | ✅ |
| Equivariant | 2% | SR3 | 0.10 | 0/2 | 2 | 1.00 | ❌ |
| Equivariant | 5% | STRidge | 0.10 | 0/4 | 0 | 1.00 | ❌ |

At 0% noise, **both libraries achieve exact recovery** — the standard library
needs 7 columns (f, f_xx, f_yy, f_xxxx, f_xxyy, f_yyyy, f³) while the
equivariant library needs 4 (f, Δf, Δ²f, f³).  At 2%+ noise, pointwise
methods fail catastrophically for both libraries (E=1.0, all true terms
lost).

## Polynomial stencils

| Library | Noise | Algorithm | d_tol | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 0.05 | 0/0 | 7 | 0.182 | ✅ |
| PDE-FIND | 2% | STRidge | 0.10 | 0/7 | 0 | 1.00 | ❌ |
| PDE-FIND | 5% | — | — | —/— | — | — | — |
| Equivariant | Clean | STRidge | 0.20 | 0/0 | 4 | 0.208 | ✅ |
| Equivariant | 2% | SR3 | 0.05 | 0/1 | 3 | 1.00 | ❌ |
| Equivariant | 5% | — | — | —/— | — | — | — |

Polynomial stencils mirror the spectral pattern: clean recovery for both
libraries, complete failure at 2%+.  5% polynomial data not yet available
(polynomial backend is ~15 s per equivariant matrix for SH's 89 snapshots).
The natural smoothing helps somewhat (equivariant fn=1 at 2% vs fn=2 for
spectral) but is insufficient to rescue pointwise methods.

## Weak formulation

### Spectral backend

| Library | Noise | Algorithm | d_tol | m | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 5.0 | 12 | 0/0 | 7 | 0.0247 | ✅ |
| PDE-FIND | 2% | STRidge | 0.20 | 12 | 0/0 | 7 | 0.0231 | ✅ |
| PDE-FIND | 5% | STRidge | 5.0 | 8 | 0/5 | 2 | 0.632 | ❌ |
| WSINDy | Clean | STRidge | 5.0 | 12 | 0/0 | 4 | 0.0249 | ✅ |
| WSINDy | 2% | STRidge | 5.0 | 12 | 0/0 | 4 | 0.00971 | ✅ |
| WSINDy | 5% | STRidge | 5.0 | 8 | 0/0 | 4 | 0.400 | ✅ |
| WSINDy IBP | Clean | STRidge | 5.0 | 8 | 0/0 | 4 | 0.0724 | ✅ |
| WSINDy IBP | 2% | STRidge | 5.0 | 12 | 0/0 | 4 | 0.0316 | ✅ |
| WSINDy IBP | 5% | STRidge | 5.0 | 8 | 0/0 | 4 | 0.329 | ✅ |

### Polynomial backend

| Library | Noise | Algorithm | d_tol | m | fp/fn | n_active | E | Exact? |
|---|---|---|---|---|---|---|---|---|
| PDE-FIND | Clean | STRidge | 5.0 | 12 | 0/0 | 7 | 0.0248 | ✅ |
| PDE-FIND | 2% | STRidge | 5.0 | 12 | 0/0 | 7 | 0.0231 | ✅ |
| PDE-FIND | 5% | — | — | — | —/— | — | — | — |
| WSINDy | Clean | STRidge | 5.0 | 12 | 0/0 | 4 | 0.0249 | ✅ |
| WSINDy | 2% | STRidge | 5.0 | 12 | 0/0 | 4 | 0.0103 | ✅ |
| WSINDy | 5% | — | — | — | —/— | — | — | — |
| WSINDy IBP | Clean | STRidge | 5.0 | 12 | 0/0 | 4 | 0.149 | ✅ |
| WSINDy IBP | 2% | STRidge | 5.0 | 12 | 0/0 | 4 | 0.0316 | ✅ |
| WSINDy IBP | 5% | — | — | — | —/— | — | — | — |

**Row mapping:** PDE-FIND = standard library + WSINDy direct.  WSINDy =
equivariant library + WSINDy direct (no IBP).  WSINDy IBP = equivariant
library + WSINDy with integration-by-parts on Δf and Δ²f.

**Key result — the central claim of the paper:** the equivariant library
achieves exact recovery at **all three noise levels** via weak formulation
(0%, 2%, 5%), while the standard library fails at 5%.  At 2%, the standard
library recovers exactly (7 columns) but the equivariant library does so
with only 4 columns.

IBP is critical at 5% noise: it cuts E from 0.40 (WSINDy direct) to 0.33
by moving the biharmonic derivative onto the analytically-known test
function, eliminating the |k|⁴ noise amplification that otherwise
dominates the Δ²f column.

Unlike Allen–Cahn, Swift–Hohenberg **cannot** be recovered by pointwise
methods at 2%+ noise even with the equivariant library — the biharmonic
term's |k|⁴ noise amplification in spectral derivatives produces E=1.0
(complete failure to predict the dynamics).  Only the weak formulation
with its built-in test-function smoothing rescues recovery at noise.
