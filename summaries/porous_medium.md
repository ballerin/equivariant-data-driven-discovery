# Porous Medium equation — recovery summary

## Equation and ground truth

$$f_t = \Delta(f^2) = \nabla\cdot(2f\nabla f) = 2f\Delta f + 2\|\nabla f\|^2$$

Two active equivariant operators:

| Operator | True coefficient |
|---|---|
| $f^2 \cdot \Delta f$ | $1.0$ |
| $f \cdot \|\nabla f\|^2$ | $2.0$ |

Standard libraries:
| Operator | True coefficient | PDE-FIND | WSINDy |
|---|---|---|---|
| $f \cdot f_{xx}$ | $2.0$ | ✅ | ✅ |
| $f \cdot f_{yy}$ | $2.0$ | ✅ | ✅ |
| $f_x^2$ | $2.0$ | ❌ | ❌ |
| $f_y^2$ | $2.0$ | ❌ | ❌ |

PME contains $f\|\nabla f\|^2$ — the standard library cannot represent it.

## Current source of truth

`results/unified_sweep_20260607_225223_efit_full.csv` (unified ablation sweep, 2026-06-07).

Simulation: N=64, L=6π, T_end=2, dt=0.002, save_every=25.  Spatial mask
excludes f < 0.05 (degenerate support boundary where D = f² → 0).

## Results — equivariant library (exact recovery at all noise levels)

| Noise | Best method | d_tol | E | n_active | cond(Θ) |
|---|---|---|---|---|---|
| 0% | WSINDy direct (m=8) | 0.05 | 0.0009 | 2 | 2.4×10⁶ |
| 2% | WSINDy direct (m=8) | 5.0 | 0.0032 | 2 | 1.9×10⁴ |
| 5% | WSINDy direct SR3 (m=8) | 0.01 | 0.14 | 2 | 6.9×10³ |

PME achieves **exact recovery at all noise levels**.  This is the only
equation besides Allen–Cahn to do so.  The second-order dynamics and
only 2 active terms make it robust.

## Results — standard library (structurally missing $f\|\nabla f\|^2$)

| Noise | Best method | n_active | E_fit |
|---|---|---|---|
| 0% | WSINDy direct (m=8) | 2 | **0.109** |
| 2% | WSINDy direct (m=8) | 60 | **0.334** |
| 5% | WSINDy direct (m=8) | 58 | **0.623** |

**Key**: on clean data, the standard library manages E_fit=0.11 with only
2 active terms — it picks two partial-derivative columns that approximate
$f^2\Delta f$ and $f\|\nabla f\|^2$ imperfectly.  At noise, it needs
58–60 columns and still leaves 33–62% of variance unexplained.

## Standard-library decomposition

The simulated PME is $f_t = \nabla\!\cdot\!(f^2\nabla f)$, which expands as:

$$f_t = f^2\Delta f + 2f\|\nabla f\|^2 = f^2(f_{xx} + f_{yy}) + 2f(f_x^2 + f_y^2)$$

(Note: this differs from the more common $f_t = \Delta(f^2) = 2f\Delta f + 2\|\nabla f\|^2$.)

| Equivariant term | True coeff | Standard-library columns | Representable? |
|---|---|---|---|
| $f^2\Delta f$ | $+1.0$ | $f^2\!\cdot\!f_{xx}$ (col 35) $+$ $f^2\!\cdot\!f_{yy}$ (col 33), coeff $+1.0$ each | ✅ 2 columns |
| $f\|\nabla f\|^2$ | $+2.0$ | **NOT REPRESENTABLE** — same reason as KS and CH: $f\!(f_x)^2$ is a product of derivatives that the $f^p\!\cdot\!\partial^a f$ library lacks. | ❌ Missing |
| **Total** | — | 2 representable, 2 missing | ❌ Incomplete |

## Remaining work

1. Forward-simulate recovered PME at 5% noise to verify dynamics.
2. The 5% recovery (E=0.14) has correct sparsity but noisy coefficients.

### Standard library

At small thresholds the standard library selects all 60 columns. fp/fn are not computed (true-coefficient vector not mapped to standard library columns). The standard library is expected to fail structurally because $f^2 f_{xx}$, $f^2 f_{yy}$, $f f_x^2$, $f f_y^2$ appear as separate columns rather than the two equivariant operators $f^2\Delta f$ and $f\|\nabla f\|^2$.

### Standard library

The standard library has not been specifically analysed for PME. Since the true operators ($f^2\Delta f$ and $f\|\nabla f\|^2$) involve nonlinear products with derivatives, their coordinate-basis representations require multiple cross-product columns (e.g., $f^2 f_{xx}$, $f^2 f_{yy}$, $f f_x^2$, $f f_y^2$, ...). The standard library is expected to spread these across many columns, requiring more active terms even in the correctly identified equation.

## Additional notes

**Sharp transition in active set**: the drop from 15 to 2 active operators between d_tol = 0.005 and d_tol = 0.01 (spectral) is striking. It suggests that most of the 15 false-positive operators have normalised coefficients clustered around the threshold range 5–25 × 0.005 = 0.125–0.625, while the two surviving operators (fp=1, fn=1) have larger normalised coefficients. Understanding which operators survive and why would clarify the regression behaviour.

**Mask effectiveness**: the quality of the support mask directly affects the regression. An overly tight mask may exclude valid data; an overly loose mask includes degenerate points that distort the regression. The current mask threshold has not been optimised.

**Short simulation time**: T_end = 2 with save_every = 25 produces only 41 snapshots. The PME evolves quickly and the support expands; after T_end = 2 the solution may be approaching the boundary of the periodic domain. More snapshots from a longer or repeated simulation might improve recovery.

## Noisy data results

The older summary stopped at clean data. The comprehensive sweep now adds 2% and 5% observation noise:

| Noise | Best result | Best recorded configuration |
|---|---|---|
| 2% | `fp=1, fn=0` | polynomial SR3 + cutoff 3.0 (`d_tol=0.005`, `lam=0.1`) |
| 5% | `fp=1, fn=0` | spectral SR3 + cutoff 3.0 (`d_tol=0.01`, `lam=0.1`) |

PME is therefore no longer "untested under noise". It remains close to exact recovery, but not exact, in the saved noisy runs.

## What needs to be done

1. Identify the single recurring false-positive operator in the noisy PME runs.
2. Add a quantitative standard-library comparison once the coordinate-basis mapping is implemented.
3. Test whether mask tuning can turn the `fp=1, fn=0` noisy cases into exact recovery.

## Status

**✅ Clean exact recovery achieved; noisy recovery is close but not exact.** Both spectral and polynomial equivariant backends recover the correct 2 operators on clean data. The saved noisy sweeps reach `fp=1, fn=0` at both 2% and 5% noise.
