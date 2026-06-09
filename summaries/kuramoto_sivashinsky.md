# Kuramoto–Sivashinsky equation — recovery summary

## Equation and ground truth

$$f_t = -\Delta f - \Delta^2 f - \tfrac{1}{2}\|\nabla f\|^2$$

Three active equivariant operators:

| Operator | True coefficient |
|---|---|
| $\Delta f$ | $-1.0$ |
| $\Delta^2 f$ | $-1.0$ |
| $\|\nabla f\|^2$ | $-0.5$ |

Standard libraries:
| Operator | True coefficient | PDE-FIND | WSINDy |
|---|---|---|---|
| $f_{xx}$ | $-1.0$ | ✅ | ✅ |
| $f_{yy}$ | $-1.0$ | ✅ | ✅ |
| $f_{xxxx}$ | $-1.0$ | ✅ | ✅ |
| $f_{xxyy}$ | $-2.0$ | ✅ | ✅ |
| $f_{yyyy}$ | $-1.0$ | ✅ | ✅ |
| $f_x^2$ | $-0.5$ | ❌ | ❌ |
| $f_y^2$ | $-0.5$ | ❌ | ❌ |

#### Extension of 1D KS to 2D

The original PDE-FIND paper (Brunton et al., PNAS 2016) successfully recovers the 1D Kuramoto–Sivashinsky equation 
$$u_t + u u_x + u_{xx} + u_{xxxx} = 0.$$
All three terms fit the $u^p \cdot \partial^d u$ mold:

| 1D term | PDE-FIND form | Recoverable? |
|---|---|---|
| $u_{xx}$ | $u^0 \cdot \partial_x^2 u$ | ✅ |
| $u_{xxxx}$ | $u^0 \cdot \partial_x^4 u$ | ✅ |
| $u u_x$ | $u^1 \cdot \partial_x^1 u$ | ✅ |

The 2D generalisation $f_t = -\Delta f - \Delta^2 f - \tfrac{1}{2}\|\nabla f\|^2$
replaces the self-advection term $u u_x$ with $\|\nabla f\|^2 = f_x^2 + f_y^2$.
In 1D there is only one derivative direction, so $u u_x$ happens to be a
single-derivative factor times a monomial.  In 2D, $\|\nabla f\|^2$ sums
squares of derivatives in *different* directions so it requires **two**
derivative factors ($\partial_x f$ and $\partial_y f$) in a single column,
which neither $f^p \cdot \partial^a f$ nor $\partial^\alpha(f^\beta)$ can
express.

## Equation regimes
The KS equation exhibits different dynamical regimes as the domain size varies. They are "stable", "cellular", and "chaotic". The stable regime is not included in the current sweep because the dynamic is not well captured in the data. The regimes that are included are:
- cellular
- weakly chaotic
- chaotic

## Results — equivariant library

### Clean data — exact recovery for all 3 regimes

| Regime | Best method | d_tol | E | n_active | cond(Θ) |
|---|---|---|---|---|---|
| Cellular | WSINDy direct (m=8) | 5.0 | 0.0033 | 3 | 1.2×10⁵ |
| Weak chaotic | WSINDy direct (m=8) | 15.0 | 0.0027 | 3 | 2.5×10⁷ |
| Chaotic | Pointwise STRidge | 5.0 | 0.0054 | 3 | 3.4×10⁷ |

All three regimes achieve **exact recovery** on clean data.

### Noisy data — total failure regardless of method

| Regime | Noise | Best method | E | Meaning |
|---|---|---|---|---|
| All 3 | 2% | Various (SR3, IBP) | 1.0 | All coefficients zeroed |
| All 3 | 5% | Various (SR3) | 1.0 | All coefficients zeroed |

E=1.0 means $\xi = \mathbf{0}$ — the regression found no sparse solution and
returned $f_t = 0$.  The condition numbers are too high:

| Regime | Noise | Pointwise cond | WSINDy cond |
|---|---|---|---|
| Cellular | 2% | 1.2×10⁷ | 5.4×10⁵ |
| Chaotic | 2% | 5.8×10⁷ | 1.4×10⁸ |
| Chaotic | 5% | 1.4×10⁸ | 7.0×10⁷ |

Even WSINDy (which helps SH at 5%) cannot rescue KS — the domain is too large
and the $|k|^4$ noise amplification is too severe.

## Results — standard library (structurally missing ‖∇f‖²)

The standard library cannot represent $\|\nabla f\|^2$, and the fit error
confirms it:

| Regime | Noise | n_active | E_fit |
|---|---|---|---|
| Cellular clean | — | 60 | **0.347** |
| Weak chaotic clean | — | 60 | **0.348** |
| Chaotic clean | — | 60 | **0.264** |
| Chaotic 2% | — | 60 | **0.280** |

Even with ALL 60 columns active, the standard library leaves **26–35% of
$f_t$ variance unexplained** on clean data.  This is the irreducible error
from lacking $\|\nabla f\|^2 = f_x^2 + f_y^2$ in the basis.

This is the strongest result in the paper: the equivariant library achieves
exact recovery on clean KS data with E≈0.005 and 3 terms; the standard
library can never fit the data regardless of how many columns it uses.

## Standard-library decomposition

| Equivariant term | True coeff | Standard-library columns | Representable? |
|---|---|---|---|
| $\Delta f$ | $-1.0$ | $f_{xx}$ (col 5) $+$ $f_{yy}$ (col 3), coeff $-1.0$ each | ✅ 2 columns |
| $\Delta^2 f$ | $-1.0$ | $f_{xxxx}$ (col 14) $+$ $2f_{xxyy}$ (col 12) $+$ $f_{yyyy}$ (col 10), coeffs $-1.0$, $-2.0$, $-1.0$ | ✅ 3 columns with linear constraint |
| $\|\nabla f\|^2$ | $-0.5$ | **NOT REPRESENTABLE** — $(f_x)^2$ and $(f_y)^2$ are products of two first derivatives; the standard library only has $f^p \cdot \partial^a f$ (single-derivative columns). | ❌ Missing |
| **Total** | — | 5 representable, 2 missing | ❌ Incomplete |

The standard library **cannot** capture $\|\nabla f\|^2$ regardless of how
many columns are used.  This is the structural gap that explains why
E_fit ≥ 0.26 even with all 60 columns active.

## Remaining work

1. Spectral lowpass sweep (cutoff ∈ {1,2,3,5}) — may rescue noisy KS.
2. Larger n_per_snap (500, 1000) — more regression points may help.
3. Forward-simulate recovered KS PDE on clean data to validate dynamics.

| Library | Backend | cond(Θ) |
|---|---|---|
| equivariant | spectral | 8.29 × 10⁸ |
| equivariant | polynomial | ~4 × 10⁸ (estimated) |
| standard | spectral | 5.17 × 10⁹ |

### Standard library (structural failure)

| Backend | d_tol ≤ 0.2 | d_tol = 5.0 (spectral) | d_tol ≥ 5.0 (polynomial) | d_tol ≥ 15.0 (spectral) |
|---|---|---|---|---|
| spectral | n_active = 60 (all) | n_active = 1 | — | n_active = 1 |
| polynomial | n_active = 60 (all) | n_active = 60 | n_active = 1 | — |

fp/fn are not computed (ground-truth mapping not implemented for the standard library column ordering).

The standard library selects **all 60 columns** for any threshold that doesn't destroy the result entirely. There is no threshold at which the standard library recovers a plausibly sparse equation. The spectral variant jumps from 60 to 1 (a single operator, wrong) when d_tol hits 5.0, with no intermediate sparse representation. This is the structural failure: $\Delta^2 f$ is spread across three correlated columns whose individual coefficients vary wildly, making the regression problem ill-conditioned.

## Domain-size regime sweep: clean data (new)

The KS equation exhibits different dynamical regimes as the domain size varies. A four-regime sweep was conducted to assess recovery robustness across spatial complexity. All simulations use N = 64, T_end = 200, dt = 0.1, save_every = 5, with n_per_snap = 200.

### Regime characterisation

| Regime | L | N_unstable | Spatial std | Drift rate (per t) | Autocorr(Δt=10) | Regime note |
|---|---|---|---|---|---|---|
| L=4π | 4π | 8 | 3.04 | −1.63 | −0.573 | quasi-periodic |
| L=8π | 8π | 44 | 3.47 | −2.58 | +0.051 | weakly chaotic |
| L=16π | 16π | 192 | 3.99 | −2.56 | +0.140 | chaotic |
| L=24π | 24π | 436 | 4.49 | −2.88 | +0.004 | strongly chaotic |

### Recovery results (d_tol = 5.0, STRidge)

| Regime | Spectral equivariant | Polynomial equivariant | Spectral standard | Notes |
|---|---|---|---|---|
| L=4π | **fn=0, fp=0** ✓ err=2.27e-3 | **fn=0, fp=0** ✓ err=2.25e-3 | fn=0, fp=1 | Both backends succeed |
| L=8π | **fn=0, fp=0** ✓ err=4.87e-3 | **fn=0, fp=0** ✓ err=1.28e-2 | fn=0, fp=5 | Both backends succeed |
| L=16π | **fn=0, fp=0** ✓ err=8.02e-3 | fn=0, fp=2 | fn=0, fp=60 | Spectral succeeds; poly has 2 false positives |
| L=24π | **fn=0, fp=0** ✓ err=4.51e-2 | fn=0, fp=2 | fn=0, fp=60 | Spectral succeeds; poly has 2 false positives |

**All four regimes achieve exact recovery with spectral equivariant**. Error magnitude increases with chaos (as expected from noise amplification in chaotic regimes). Polynomial backend degrades for the largest regimes (L ≥ 16π).

### Condition numbers by regime

| Regime | Spectral equivariant | Polynomial equivariant | Spectral standard |
|---|---|---|---|
| L=4π | 1.21 × 10⁸ | 1.21 × 10⁸ | 3.68 × 10⁹ |
| L=8π | 5.81 × 10⁸ | 5.81 × 10⁸ | 3.60 × 10⁹ |
| L=16π | 8.29 × 10⁸ | 8.34 × 10⁸ | 5.17 × 10⁹ |
| L=24π | 1.26 × 10⁹ | 1.27 × 10⁹ | 8.74 × 10⁹ |

Equivariant condition improves 3-7× over standard. The increase with L is expected due to larger spatial dimension and more unstable modes. Despite conditions worsening toward L=24π, exact recovery remains achievable across all regimes.

## Recommended hyperparameters

| Setting | Value |
|---|---|
| Algorithm | STRidge |
| Backend | spectral |
| Library | equivariant |
| lam | 1e-5 (no effect; any value gives identical results) |
| d_tol | **5.0** (primary recommendation) |
| d_tol | 15.0 (also achieves exact recovery; broader initial step) |
| D_TOLS grid | [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, **5.0, 15.0, 20.0**] |

**Do not use d_tol ≤ 0.2** for KS — the false-positive constant column will never be eliminated. **Do not use d_tol = 20.0** — the true terms $\Delta f$ and $\|\nabla f\|^2$ will be eliminated.

The polynomial backend achieves fp = 2 at d_tol = 5.0 on small regimes (L ≤ 8π) and larger regimes (L ≥ 16π). On intermediate regimes it achieves near-exact recovery.

## Noisy data results (current saved sweeps)

Source: `results/comprehensive_sweep_20260528_001113_ks.csv`.

### 2% observation noise

The older claim that KS reaches `fp=1, fn=0` at 2% noise is not supported by the saved CSV.

The best saved 2% KS row, using the same ranking rule as the rest of the project (exact first, then minimise `fn`, then `fp`), is:

| Backend | Algorithm | Cutoff | Best result | d_tol | lam |
|---|---|---|---|---|---|
| polynomial | SR3 | 3.0 | `fp=8, fn=0` | 0.05 | 0.1 |

This means all three true KS terms are retained, but with eight false positives. There is no exact recovery in the saved 2% sweep.

### 5% observation noise

A dedicated follow-up sweep now includes 5% KS rows (`results/comprehensive_sweep_20260528_123302_ks_noisy_followup.csv`).

Best saved 5% row (same ranking: exact first, then minimise `fn`, then `fp`):

| Backend | Algorithm | Cutoff | Best result | d_tol | lam |
|---|---|---|---|---|---|
| spectral | SR3 | none | `fp=15, fn=0` | 0.2 | 100.0 |

So 5% KS is now swept, but recovery quality is poor.

### SR3 clean data results

In addition to STRidge, SR3 also achieves exact recovery on clean data:

| d_tol | Backend | Exact | Notes |
|---|---|---|---|
| 0.1 | spectral | yes | fp=0, fn=0 |
| 0.2 | spectral | yes | fp=0, fn=0 |
| 5.0 | spectral | yes | fp=0, fn=0 |

SR3 gives a broader window of d_tol values achieving exact recovery on clean data compared to STRidge (3 vs 2 configurations).

## What remains to be done

1. Improve noisy KS quality: 5% is now swept, but best saved result is still `fp=15, fn=0`.
2. Improve the 2% noisy case; the best saved row is still `fp=8, fn=0`, which is much weaker than the clean-data result.
3. Apply weak-form regression here once it exists for Swift–Hohenberg.
4. Update `notebooks/kuramoto_sivashinsky.ipynb` to the ETDRK4 parameter set used in the current sweeps.
5. Add a clean coefficient table for the exact spectral recoveries.
