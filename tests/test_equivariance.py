"""
Tests for O(2)-equivariance of multigraph operators and the equivariant library.

The central claim of the paper is that operators built from contractions of
covariant derivatives via the Euclidean metric are O(2)-equivariant, while
coordinate-basis operators (f_x, f_xx, ...) are not.  These tests verify:

1. Rotation equivariance (90°, 180°, 270°) on a square periodic grid.
2. Reflection equivariance (x → −x, y → −y, x ↔ y).
3. Translation equivariance (np.roll).
4. Non-equivariance of coordinate-basis partial derivatives.
5. Library-level equivariance under rotation.

All tests use the spectral backend (exact for band-limited inputs on a
periodic grid), which gives machine-precision agreement for equivariant
operators.  Tests that should fail (non-equivariant operators) use
assert_raises or explicit inequality checks.

Conventions
-----------
- R_90[f][i, j] = f[j, N−1−i]    (np.rot90, 90° CCW)
- R_180[f]       = np.rot90(f, 2)
- R_270[f]       = np.rot90(f, 3)
- Reflect_x[f]   = f[−i, j]       (np.flip(f, axis=0))
- Reflect_y[f]   = f[i, −j]       (np.flip(f, axis=1))
- Swap_xy[f]     = f[j, i]        (f.T, diagonal reflection)

An operator O is O(2)-equivariant if

    O[f ∘ g] = O[f] ∘ g      for all g ∈ O(2).

For a discrete square grid the rotations and reflections above are exact
(no interpolation needed), so the equality can be checked to machine
precision when using the spectral backend.
"""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives as _derivatives
from equivariant_pde_find.multigraph_operators import _op_bilap, _op_Delta, _op_grad_sq, _op_hess_norm_sq


# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------


def make_grid(N, L=2 * math.pi):
    x = np.linspace(0.0, L, N, endpoint=False)
    y = np.linspace(0.0, L, N, endpoint=False)
    return np.meshgrid(x, y, indexing="ij")


def rotation_90(f):
    """Rotate a 2-D field 90° CCW (exact on square grid)."""
    return np.rot90(f, k=1)


def rotation_180(f):
    return np.rot90(f, k=2)


def rotation_270(f):
    return np.rot90(f, k=3)


def reflect_x(f):
    return np.flip(f, axis=0)


def reflect_y(f):
    return np.flip(f, axis=1)


def swap_xy(f):
    return f.T


# Transformations as (name, transform_fn) pairs.
ROTATIONS = [
    ("R_90", rotation_90),
    ("R_180", rotation_180),
    ("R_270", rotation_270),
]

REFLECTIONS = [
    ("refl_x", reflect_x),
    ("refl_y", reflect_y),
    ("swap_xy", swap_xy),
]

ALL_TRANSFORMS = ROTATIONS + REFLECTIONS


# ---------------------------------------------------------------------------
# Fixture: a non-axis-aligned test field (so rotations are non-trivial)
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _set_spectral_backend():
    """All equivariance tests use the spectral backend for exactness."""
    _derivatives.set_backend("spectral")
    yield
    # Leave spectral as the default; other test files may rely on it.


@pytest.fixture
def field_and_domain():
    """A smooth, non-axis-aligned field on a square periodic domain."""
    N = 64
    L = 2 * math.pi
    X, Y = make_grid(N, L)
    mx, my = 2, 3  # non-axis-aligned → rotation changes the field
    f = np.cos(mx * X + my * Y) + 0.5 * np.sin(2 * mx * X - my * Y)
    return f, L, N


# ---------------------------------------------------------------------------
# 1.  Rotation equivariance of key scalar operators
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("op, name", [
    (_op_Delta, "Δf"),
    (_op_bilap, "Δ²f"),
    (_op_grad_sq, "‖∇f‖²"),
    (_op_hess_norm_sq, "‖∇²f‖²"),
])
@pytest.mark.parametrize("tf_name, tf", ROTATIONS)
def test_rotation_equivariance_scalar_ops(op, name, tf_name, tf, field_and_domain):
    """O[f ∘ R] = O[f] ∘ R  for all scalar O(2)-equivariant operators."""
    f, L, _N = field_and_domain

    # LHS: rotate field, then apply operator
    f_rot = tf(f)
    lhs = op(f_rot, L, L)

    # RHS: apply operator, then rotate result
    op_f = op(f, L, L)
    rhs = tf(op_f)

    # atol=1e-9: bilaplacian (k⁴ multiplier) accumulates ~6e-10 round-off
    # on a 64² grid; 1e-9 is still 9 orders below the signal amplitude.
    np.testing.assert_allclose(lhs, rhs, rtol=1e-12, atol=1e-9,
                               err_msg=f"{name} not {tf_name}-equivariant")


# ---------------------------------------------------------------------------
# 2.  Reflection equivariance
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("op, name", [
    (_op_Delta, "Δf"),
    (_op_bilap, "Δ²f"),
    (_op_grad_sq, "‖∇f‖²"),
    (_op_hess_norm_sq, "‖∇²f‖²"),
])
@pytest.mark.parametrize("tf_name, tf", REFLECTIONS)
def test_reflection_equivariance_scalar_ops(op, name, tf_name, tf, field_and_domain):
    """O[f ∘ S] = O[f] ∘ S  for reflections across axes and diagonal."""
    f, L, _N = field_and_domain

    f_refl = tf(f)
    lhs = op(f_refl, L, L)
    rhs = tf(op(f, L, L))

    np.testing.assert_allclose(lhs, rhs, rtol=1e-12, atol=1e-9,
                               err_msg=f"{name} not {tf_name}-equivariant")


# ---------------------------------------------------------------------------
# 3.  Translation equivariance  (periodic shifts)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("op, name", [
    (_op_Delta, "Δf"),
    (_op_bilap, "Δ²f"),
    (_op_grad_sq, "‖∇f‖²"),
    (_op_hess_norm_sq, "‖∇²f‖²"),
])
@pytest.mark.parametrize("dx, dy", [(3, 0), (0, 5), (7, 11), (-2, 4)])
def test_translation_equivariance(op, name, dx, dy, field_and_domain):
    """O[f(· + τ)] = O[f](· + τ) — operators are translation-invariant."""
    f, L, _N = field_and_domain

    f_shifted = np.roll(np.roll(f, dx, axis=0), dy, axis=1)
    lhs = op(f_shifted, L, L)
    rhs = np.roll(np.roll(op(f, L, L), dx, axis=0), dy, axis=1)

    np.testing.assert_allclose(lhs, rhs, rtol=1e-12, atol=1e-9,
                               err_msg=f"{name} not translation-equivariant by ({dx},{dy})")


# ---------------------------------------------------------------------------
# 4.  Non-equivariance of coordinate-basis partial derivatives
# ---------------------------------------------------------------------------


def _partial_x_spectral(f, L):
    """Compute ∂f/∂x via FFT (not O(2)-equivariant as a scalar operator)."""
    from equivariant_pde_find.derivatives import spectral as sd
    fx, _fy = sd.spectral_grad(f, L, L)
    return fx


def _partial_xx_spectral(f, L):
    """Compute ∂²f/∂x² via FFT (not O(2)-equivariant as a scalar operator)."""
    from equivariant_pde_find.derivatives import spectral as sd
    fxx, _fxy, _fyy = sd.spectral_hessian(f, L, L)
    return fxx


def test_fx_not_rotation_equivariant(field_and_domain):
    """f_x is NOT equivariant under 90° rotation: (f∘R)_x ≠ (f_x)∘R.

    Under a 90° CCW rotation, ∂/∂x maps to ∂/∂y (up to a sign depending on
    convention).  This test verifies that the naive equivariance check fails.
    """
    f, L, _N = field_and_domain

    f_rot = rotation_90(f)
    lhs = _partial_x_spectral(f_rot, L)          # (f∘R)_x
    rhs = rotation_90(_partial_x_spectral(f, L))  # (f_x)∘R

    # These must differ substantially (not be close).
    assert not np.allclose(lhs, rhs, rtol=1e-3), (
        "f_x unexpectedly appears 90°-rotation-equivariant "
        "(max |diff| = {:.2e})".format(np.max(np.abs(lhs - rhs)))
    )


def test_fxx_not_rotation_equivariant(field_and_domain):
    """f_xx is NOT equivariant under 90° rotation.  Under R_90, ∂²/∂x² → ∂²/∂y²."""
    f, L, _N = field_and_domain

    f_rot = rotation_90(f)
    lhs = _partial_xx_spectral(f_rot, L)
    rhs = rotation_90(_partial_xx_spectral(f, L))

    assert not np.allclose(lhs, rhs, rtol=1e-3), (
        "f_xx unexpectedly appears 90°-rotation-equivariant "
        "(max |diff| = {:.2e})".format(np.max(np.abs(lhs - rhs)))
    )


def test_fx_maps_to_fy_under_rotation(field_and_domain):
    """Under 90° CCW: (f∘R)_x = −(f_y)∘R.

    np.rot90(k=1) maps physical coordinates (x, y) → (y, −x) (modulo the
    periodic grid).  By the chain rule:

        ∂_x (f∘R)(x,y) = ∂_x [f(y, −x)] = f_y(y, −x) · (−1) = −(f_y ∘ R)(x,y).

    So the x-derivative of the rotated field equals the negative of the
    rotated y-derivative of the original field.
    """
    f, L, _N = field_and_domain

    from equivariant_pde_find.derivatives import spectral as sd

    f_rot = rotation_90(f)
    lhs = _partial_x_spectral(f_rot, L)        # (f∘R)_x

    _fx, fy = sd.spectral_grad(f, L, L)
    rhs = -rotation_90(fy)                       # −(f_y)∘R

    np.testing.assert_allclose(lhs, rhs, rtol=1e-12, atol=1e-12,
                               err_msg="(f∘R)_x ≠ −(f_y)∘R under np.rot90(k=1)")


# ---------------------------------------------------------------------------
# 5.  Library-level equivariance  (regression design matrix)
# ---------------------------------------------------------------------------


def test_equivariant_library_rotated_field(field_and_domain):
    """Each column of build_multigraph_library(f∘R) equals the corresponding
    column of build_multigraph_library(f) rotated by R.

    This is the key structural property: the design matrix Θ transforms
    column-wise under rotations of the input field.
    """
    f, L, N = field_and_domain
    from equivariant_pde_find.library import build_multigraph_library

    Theta, descs = build_multigraph_library(f, L, L, max_p=2, max_n=3)
    assert Theta.shape[1] == 20  # sanity

    f_rot = rotation_90(f)
    Theta_rot, descs_rot = build_multigraph_library(f_rot, L, L, max_p=2, max_n=3)
    assert descs_rot == descs  # same columns, same order

    # Each column of Theta_rot should equal the rotated column of Theta.
    for j in range(Theta.shape[1]):
        col_rot_expected = rotation_90(Theta[:, j].reshape(N, N)).ravel()
        np.testing.assert_allclose(
            Theta_rot[:, j], col_rot_expected,
            rtol=1e-11, atol=1e-8,
            err_msg=f"Library column '{descs[j]}' not 90°-equivariant at library level.",
        )


def test_standard_library_not_columnwise_equivariant(field_and_domain):
    """The standard (coordinate-basis) library does NOT transform column-wise
    under rotation.  At least one column must differ substantially.

    This is the central negative result: the coordinate splitting of
    isotropic derivatives destroys the structural sparsity.
    """
    f, L, N = field_and_domain
    from equivariant_pde_find.library import build_pdefind_library

    Theta, descs = build_pdefind_library(f, L, L)
    f_rot = rotation_90(f)
    Theta_rot, _ = build_pdefind_library(f_rot, L, L)

    # Check how many columns fail the equivariance test.
    n_fail = 0
    for j in range(Theta.shape[1]):
        col_rot_expected = rotation_90(Theta[:, j].reshape(N, N)).ravel()
        if not np.allclose(Theta_rot[:, j], col_rot_expected, rtol=1e-3):
            n_fail += 1

    assert n_fail > 0, (
        "All standard-library columns appear rotation-equivariant — "
        "unexpected for coordinate-basis operators like f_x, f_xx."
    )


# ---------------------------------------------------------------------------
# 6.  Composite equivariance: all operators in the full basis
# ---------------------------------------------------------------------------


def test_all_equivariant_operators_rotation_equivariant(field_and_domain):
    """Every operator returned by get_operators(max_p=2, max_n=3) must be
    equivariant under 90° rotation.  This is the exhaustive check."""
    f, L, N = field_and_domain
    from equivariant_pde_find.multigraph_operators import get_operators

    ops = get_operators(2, 3)
    assert len(ops) == 20

    f_rot = rotation_90(f)

    for func, desc in ops:
        lhs = func(f_rot, L, L)
        rhs = rotation_90(func(f, L, L))
        np.testing.assert_allclose(
            lhs, rhs, rtol=1e-11, atol=1e-8,
            err_msg=f"Operator '{desc}' not 90°-rotation-equivariant.",
        )


def test_all_equivariant_operators_reflection_equivariant(field_and_domain):
    """Every operator in the basis must be equivariant under x-axis reflection."""
    f, L, N = field_and_domain
    from equivariant_pde_find.multigraph_operators import get_operators

    ops = get_operators(2, 3)
    f_refl = reflect_x(f)

    for func, desc in ops:
        lhs = func(f_refl, L, L)
        rhs = reflect_x(func(f, L, L))
        np.testing.assert_allclose(
            lhs, rhs, rtol=1e-11, atol=1e-8,
            err_msg=f"Operator '{desc}' not x-reflection-equivariant.",
        )


# ---------------------------------------------------------------------------
# 7.  Regression equivariance: same coefficients on rotated data
# ---------------------------------------------------------------------------


def test_regression_coefficients_rotation_invariant(field_and_domain):
    """STRidge on f and on R[f] recovers the same coefficient pattern.

    For a PDE whose ground truth is built from equivariant operators
    (here: Allen-Cahn f_t = Δf + f − f³), the regression should select
    the same columns regardless of whether the data is rotated.

    We use a simple synthetic design matrix rather than a full simulation
    to keep the test fast and deterministic.
    """
    f, L, N = field_and_domain
    from equivariant_pde_find.library import build_multigraph_library
    from equivariant_pde_find.regression import STRidge

    # Build a synthetic f_t from the true Allen-Cahn operators.
    Theta, descs = build_multigraph_library(f, L, L, max_p=2, max_n=3)
    # True coefficients: Δf (+1), f (+1), f³ (−1).  Find their column indices.
    idx_Delta = descs.index("Δf")
    idx_f = descs.index("f")
    idx_f3 = descs.index("f³")

    xi_true = np.zeros(Theta.shape[1])
    xi_true[idx_Delta] = 1.0
    xi_true[idx_f] = 1.0
    xi_true[idx_f3] = -1.0

    f_t = Theta @ xi_true  # synthetic time derivative

    # STRidge on original data.
    # Signature: STRidge(X0, y, lam, maxit, tol, normalize=2)
    xi_rec = STRidge(Theta, f_t, lam=1e-5, maxit=25, tol=1.0, normalize=0)

    # Now rotate everything and re-run
    f_rot = rotation_90(f)
    Theta_rot, descs_rot = build_multigraph_library(f_rot, L, L, max_p=2, max_n=3)
    assert descs_rot == descs

    f_t_rot = rotation_90(f_t.reshape(N, N)).ravel()
    xi_rot = STRidge(Theta_rot, f_t_rot, lam=1e-5, maxit=25, tol=1.0, normalize=0)

    # The recovered sparse patterns should be identical.
    active = set(np.where(np.abs(xi_rec) > 1e-10)[0])
    active_rot = set(np.where(np.abs(xi_rot) > 1e-10)[0])
    assert active == active_rot, (
        f"Active column mismatch: original {active}, rotated {active_rot}"
    )

    # Coefficients should be nearly identical.
    np.testing.assert_allclose(xi_rec, xi_rot, rtol=1e-10, atol=1e-10,
                               err_msg="Regression coefficients differ under rotation.")


# ---------------------------------------------------------------------------
# 8.  Identity and constant are trivially equivariant
# ---------------------------------------------------------------------------


def test_identity_and_constant_equivariant(field_and_domain):
    """f and 1 are trivially equivariant under all transformations."""
    f, L, _N = field_and_domain
    from equivariant_pde_find.multigraph_operators import _op_const, _op_identity

    for tf_name, tf in ALL_TRANSFORMS:
        # Identity: id(f∘T) = (id f)∘T  →  f∘T = f∘T  (tautology)
        lhs_id = _op_identity(tf(f), L, L)
        rhs_id = tf(_op_identity(f, L, L))
        np.testing.assert_array_equal(lhs_id, rhs_id)

        # Constant: 1(f∘T) = 1(f)∘T  →  ones = ones  (tautology)
        lhs_c = _op_const(tf(f), L, L)
        rhs_c = tf(_op_const(f, L, L))
        np.testing.assert_array_equal(lhs_c, rhs_c)
