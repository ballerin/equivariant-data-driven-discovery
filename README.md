# Equivariant PDE discovery

## Summary
This directory contains the core Python code used to
- simulate benchmark PDEs, optionally with noise;
- construct an O(2)-equivariant library of nonlinear partial differential operators;
- build the standard coordinate libraries (non-equivariant) for comparison;
- apply sparse regression to recover PDEs from data;
- 2D flat domain with periodic boundaries conditions.

The implementation provides several differentiation backends (spectral, local-polynomial stencils, isotropic discretization stencils) so experiments can trade accuracy for noise-robustness.

Both SINDy and WSINDy approaches available.

## Motivation
The code implements the multigraph-operator approach from Ballerin & Grong [1]: a minimal, complete, rotation-equivariant basis of polynomial nonlinear partial differential operators. It is organised to make it straightforward to (re)generate datasets, assemble library matrices, and run hyperparameter sweeps that compare the equivariant basis from the paper against standard non-equivariant libraries such as those used in SINDy/PDE-FIND  and WSINDy.

## Repository structure

TODO

## Error metric

Recovery is measured in terms of how many false positives (fp) and false negatives (fn) are present in the recovered PDE, as well as the relative L2 error between the recovered PDE and the true PDE, evaluated on the same dataset:

$$E = \frac{\|\Theta (\xi - \xi_{\mathrm{true}})\|_2}{\|\Theta \xi_{\mathrm{true}}\|_2}$$

where $\Theta$ is the library matrix, $\xi$ is the recovered coefficient vector, and $\xi_{\mathrm{true}}$ is the ground-truth coefficient vector. The denominator uses $\Theta \xi_{\mathrm{true}}$ (the true dynamics), *not* the observed $f_t$, so the metric is independent of observation noise in the time derivative.

This is a good metric if the PDE to recover is representable in the library, but it can be misleading if the library is missing key operators. In that case, a low E_fit (fit error) can be achieved by using many spurious terms to approximate the missing dynamics, which is what happens sometimes.

## Libraries

We compare classes of libraries with maximum order `max_p=2` and maximum polynomial degree `max_n=3` for the 2D flat settings.

Notice that in general the terms in the library are not linearly independent. In particular, the equivariant library has 20 operators but only 18 are linearly independent as shown in [1].

| Library | Form | # terms | $O(2)$-equivariant | Linear | Spanning |
|---|---|---|---|---|---|
| PDE-FIND | $f^p \cdot \partial^a f$ | 60| ❌ | ❌ | ❌ |
| WSINDy | $\partial^\alpha f^\beta$ | 46 | ❌ | ✅ | ❌ |
| Polynomial jet | monomials in jet variables | 105 | ❌ | ❌ | ✅ |
| Multigraph | multigraph operators | 20 | ✅ | ❌ | ✅ |

### Overlap between PDE-FIND and WSINDy

Both libraries share the pure monomials ($1$, $f$, $f^2$, $f^3$) and the pure
derivatives ($\partial_x f$, $\partial_{xx} f$, …, $\partial_{xxxx} f$).
The structural difference is in how they handle nonlinear derivative terms:

* **PDE-FIND unique** (42 columns): $f^p \cdot \partial^a f$ for $p \ge 1$.
  Multiplies a monomial by a **single** derivative factor.  Simple and
  interpretable, but cannot represent products of derivatives like
  $\|\nabla f\|^2 = f_x^2 + f_y^2$.

* **WSINDy unique** (28 columns): $\partial^a(f^q)$ for $q \ge 2$.
  Differentiates a monomial power.  By the chain rule these expand to
  **mixtures** of terms — e.g.\ $\partial_x^2(f^2) = 2f_x^2 + 2f f_{xx}$
  bundles the product-of-derivatives term $f_x^2$ with the PDE-FIND-style
  term $f f_{xx}$.  WSINDy can *reach* $f_x^2$ but cannot isolate it as a
  sparse column.

Consequently, **neither** standard library can represent $\|\nabla f\|^2$
as an independent sparse term.

## Equivariant operator library — canonical ordering

`get_operators(max_p, max_n)` returns operators in a **fixed three-level
lexicographic order**:

| Level | Key | Meaning |
|-------|-----|---------|
| 1 (primary) | `p` | total edge count in the multigraph (= half the differential order, i.e. order = 2p) |
| 2 (secondary) | `k` | polynomial degree (total vertex count, non-isolated + isolated) |
| 3 (tertiary) | description string | Unicode lexicographic order — deterministic tie-breaker within each (p, k) group |

This ordering is hardcoded, and is **enforced by `tests/test_multigraph_operators.py`**.  Any change to the operator set or ordering must be accompanied by an update to that test.

For `get_operators(max_p=2, max_n=3)` the 20 operators are, in order:

| Index | p | k | Description |
|-------|---|---|-------------|
| 0  | 0 | 0 | `1` |
| 1  | 0 | 1 | `f` |
| 2  | 0 | 2 | `f²` |
| 3  | 0 | 3 | `f³` |
| 4  | 1 | 1 | `Δf` |
| 5  | 1 | 2 | `f·Δf` |
| 6  | 1 | 2 | `‖∇f‖²` |
| 7  | 1 | 3 | `f²·Δf` |
| 8  | 1 | 3 | `f·‖∇f‖²` |
| 9  | 2 | 1 | `Δ²f` |
| 10 | 2 | 2 | `f·Δ²f` |
| 11 | 2 | 2 | `Δf·Δf` |
| 12 | 2 | 2 | `‖∇²f‖²` |
| 13 | 2 | 2 | `⟨∇Δf,∇f⟩` |
| 14 | 2 | 3 | `f²·Δ²f` |
| 15 | 2 | 3 | `f·Δf·Δf` |
| 16 | 2 | 3 | `f·‖∇²f‖²` |
| 17 | 2 | 3 | `f·⟨∇Δf,∇f⟩` |
| 18 | 2 | 3 | `‖∇f‖²·Δf` |
| 19 | 2 | 3 | `∇²f(∇f,∇f)` |


# References
[1]: Ballerin, F. & Grong, E. (2026). Equivariant nonlinear partial differential operators on constant curvature spaces. *arXiv preprint arXiv:2605.16847*. https://arxiv.org/abs/2605.16847