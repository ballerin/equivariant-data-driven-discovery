# Cahn–Hilliard equation — recovery summary

## Equation and ground truth

$$f_t = -\Delta(\Delta f - f + f^3) = -\varepsilon^2\Delta^2 f + \Delta f - 3f^2\Delta f - 6f\|\nabla f\|^2$$

Four active equivariant operators (eps=0.2):

| Operator | True coefficient |
|---|---|
| $\Delta f$ | $-1.0$ |
| $f^2\Delta f$ | $+3.0$ |
| $f\|\nabla f\|^2$ | $+6.0$ |
| $\Delta^2 f$ | $-\varepsilon^2 = -0.04$ |

Standard libraries:
| Operator | True coefficient | PDE-FIND | WSINDy |
|---|---|---|---|
| $f_{xx}$ | $+1.0$ | ✅ | ✅ |
| $f_{yy}$ | $+1.0$ | ✅ | ✅ |
| $f^2 \cdot f_{xx}$ | $-3.0$ | ✅ | ✅ |
| $f^2 \cdot f_{yy}$ | $-3.0$ | ✅ | ✅ |
| $f \cdot f_x^2$ | $-6.0$ | ❌ | ❌ |
| $f \cdot f_y^2$ | $-6.0$ | ❌ | ❌ |
| $f_{xxxx}$ | $-\varepsilon^2$ | ✅ | ✅ |
| $f_{xxyy}$ | $-2\varepsilon^2$ | ✅ | ✅ |
| $f_{yyyy}$ | $-\varepsilon^2$ | ✅ | ✅ |

CH contains $f\|\nabla f\|^2 = f(f_x^2 + f_y^2)$, a product of two derivatives
multiplied by $f$ — this **does not exist** in the standard library.

## Current source of truth

`results/unified_sweep_20260607_225223_efit_full.csv` (unified ablation sweep, 2026-06-07).

Simulation: N=128, L=4π, T_end=30, eps=0.2, dt=0.005, save_every=30.

## Results — equivariant library

| Noise | Best method | d_tol | E | n_active |
|---|---|---|---|---|
| 0% | WSINDy direct (m=8) | 15.0 | 0.049 | 4 |
| 2% | Pointwise STRidge | 5.0 | 1.0 | 0 |
| 5% | Pointwise STRidge | 5.0 | 1.0 | 0 |

**Clean data**: exact recovery (fp=0, fn=0) with E=0.049.  The relatively
high E (compared to AC at 0.005) reflects the small $\varepsilon^2=0.04$
coefficient — a 5% relative error on $\Delta^2 f$ is hard to resolve.

**Noisy data**: total failure (E=1.0).  The biharmonic term at N=128 is stiff
and the small $\varepsilon^2$ coefficient makes it the first term lost to noise.

## Results — standard library (structurally missing $f\|\nabla f\|^2$)

| Noise | Best method | n_active | E_fit |
|---|---|---|---|
| 0% | WSINDy direct SR3 (m=8) | 60 | **0.693** |
| 2% | WSINDy direct (m=8) | 60 | **0.707** |
| 5% | WSINDy direct SR3 (m=8) | 60 | **0.736** |

**Key result**: even with ALL 60 columns, the standard library leaves
**69–74% of $f_t$ variance unexplained**.  This is the largest irreducible
error in any benchmark — $f\|\nabla f\|^2$ is a dominant term in CH dynamics
and the standard library simply lacks it.

## Condition numbers

| Noise | Pointwise equivariant | Pointwise standard |
|---|---|---|
| 0% | 1.6×10⁴ | 1.3×10⁴ |
| 2% | 3.6×10⁵ | 1.2×10⁵ |
| 5% | 6.6×10⁵ | 2.1×10⁵ |

## Standard-library decomposition

| Equivariant term | True coeff | Standard-library columns | Representable? |
|---|---|---|---|
| $\Delta f$ | $-1.0$ | $f_{xx}$ (col 5) $+$ $f_{yy}$ (col 3), coeff $-1.0$ each | ✅ 2 columns |
| $-\varepsilon^2\Delta^2 f$ | $-\varepsilon^2$ | $f_{xxxx}$ (col 14) $+$ $2f_{xxyy}$ (col 12) $+$ $f_{yyyy}$ (col 10), coeffs $-\varepsilon^2$, $-2\varepsilon^2$, $-\varepsilon^2$ | ✅ 3 columns with constraint |
| $-3f^2\Delta f$ | $-3.0$ | $f^2\!\cdot\!f_{xx}$ (col 35) $+$ $f^2\!\cdot\!f_{yy}$ (col 33), coeff $-3.0$ each | ✅ 2 columns |
| $-6f\|\nabla f\|^2$ | $-6.0$ | **NOT REPRESENTABLE** — $f\!(f_x)^2$ is $f$ times the square of a derivative; the standard library has $f\!\cdot\!f_{xx}$ (col 20) but not $f\!\cdot\!(f_x)^2$. | ❌ Missing |
| **Total** | — | 7 representable, 2 missing | ❌ Incomplete |

The term $-6f\|\nabla f\|^2$ is the dominant nonlinearity in CH dynamics
and accounts for most of the E_fit=0.69 irreducible error.

## Remaining work

1. Test eps=0.4 and eps=0.1 with WSINDy — may recover at larger eps.
2. Forward-simulate clean recovery to verify dynamics.
3. Spectral lowpass sweep — the biharmonic stiffness may benefit.
3. Treat weak-form regression as the main path for improving noisy CH, especially for small `eps`.
