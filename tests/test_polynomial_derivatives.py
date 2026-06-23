"""
Tests for `polynomial_derivatives.py` (local least-squares polynomial fits).

These tests check that low-order polynomials are recovered exactly by the
local fits (within floating-point rounding), and that the polynomial
approximation behaves as expected on low-frequency trigonometric modes.
"""

import math
import os
import sys

import numpy as np
import pytest

# Ensure the code directory is importable when tests are run from the repo root.
CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import polynomial_derivatives as pd  # noqa: E402


def make_grid(Nx, Ny, Lx, Ly):
    x = np.linspace(0.0, Lx, Nx, endpoint=False)
    y = np.linspace(0.0, Ly, Ny, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="ij")
    return X, Y


def test_polynomial_grad_linear_exactness():
    """Degree-1 polynomials should be recovered exactly by the local fit."""
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    a, b, c = 1.234, -0.987, 0.5
    f = a * X + b * Y + c

    fx_p, fy_p = pd.polynomial_grad(f, Lx, Ly)
    # The polynomial backend uses periodic wrapping when building local
    # stencils. A globally non-periodic polynomial creates a jump at the
    # domain boundary which contaminates fits near the edges. Restrict the
    # exactness check to interior points where the stencil does not wrap.
    _, r = pd.get_default_stencil(1)  # use the actual default stencil half-width
    interior_fx = fx_p[r:-r, r:-r]
    interior_fy = fy_p[r:-r, r:-r]
    np.testing.assert_allclose(interior_fx, a * np.ones_like(interior_fx), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(interior_fy, b * np.ones_like(interior_fy), rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize(
    "Lx,Ly,Nx,Ny,mx,my",
    [
        (2 * math.pi, 2 * math.pi, 128, 128, 1, 1),
        (8 * math.pi, 8 * math.pi, 128, 128, 1, 1),
    ],
)
def test_polynomial_derivatives_approximate(Lx, Ly, Nx, Ny, mx, my):
    """Polynomial local fits approximate low-frequency trigonometric modes."""
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    theta = kx * X + ky * Y
    f = np.cos(theta)

    # Gradient
    fx_true = -kx * np.sin(theta)
    fy_true = -ky * np.sin(theta)
    fx_p, fy_p = pd.polynomial_grad(f, Lx, Ly)
    np.testing.assert_allclose(fx_p, fx_true, rtol=1e-3, atol=1e-6)
    np.testing.assert_allclose(fy_p, fy_true, rtol=1e-3, atol=1e-6)

    # Laplacian
    lap_true = -(kx**2 + ky**2) * np.cos(theta)
    lap_p = pd.polynomial_laplacian(f, Lx, Ly)
    np.testing.assert_allclose(lap_p, lap_true, rtol=5e-3, atol=1e-6)

    # Hessian
    fxx_true = -(kx**2) * np.cos(theta)
    fxy_true = -kx * ky * np.cos(theta)
    fyy_true = -(ky**2) * np.cos(theta)
    fxx_p, fxy_p, fyy_p = pd.polynomial_hessian(f, Lx, Ly)
    np.testing.assert_allclose(fxx_p, fxx_true, rtol=5e-3, atol=1e-6)
    np.testing.assert_allclose(fxy_p, fxy_true, rtol=5e-3, atol=1e-6)
    np.testing.assert_allclose(fyy_p, fyy_true, rtol=5e-3, atol=1e-6)


def test_poly_partial_1d_against_analytic():
    L = 2 * math.pi
    N = 256
    m = 3
    x = np.linspace(0.0, L, N, endpoint=False)
    k = 2.0 * math.pi * m / L
    u = np.cos(k * x)

    # 1st derivative analytic
    du_true = -k * np.sin(k * x)
    du_p = pd.poly_partial_1d(u, L, 1)
    np.testing.assert_allclose(du_p, du_true, rtol=2e-3, atol=1e-6)

    # 2nd derivative analytic
    d2_true = -(k**2) * np.cos(k * x)
    d2_p = pd.poly_partial_1d(u, L, 2)
    np.testing.assert_allclose(d2_p, d2_true, rtol=2e-3, atol=1e-6)


def test_polynomial_constant_field_zero():
    """Derivatives of a constant field should be (numerically) zero."""
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 32, 32
    f = np.ones((Nx, Ny))

    fx, fy = pd.polynomial_grad(f, Lx, Ly)
    np.testing.assert_allclose(fx, 0.0, atol=1e-12)
    np.testing.assert_allclose(fy, 0.0, atol=1e-12)


def test_poly_partial_matches_specialised_call():
    """`poly_partial(f,1,0)` should match `polynomial_grad(f)[0]`."""
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 48, 48
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx, ky = 2.0 * math.pi * 2 / Lx, 2.0 * math.pi * 3 / Ly
    f = np.cos(kx * X + ky * Y)

    fx_p, _ = pd.polynomial_grad(f, Lx, Ly)
    fx_generic = pd.poly_partial(f, Lx, Ly, 1, 0)
    np.testing.assert_allclose(fx_p, fx_generic, rtol=1e-3, atol=1e-6)


def test_shape_error_on_wrong_input():
    """Passing a 1-D array to 2-D routines should raise ValueError."""
    u = np.arange(16)
    with pytest.raises(ValueError):
        pd.polynomial_grad(u, 2 * math.pi, 2 * math.pi)


def test_poly_partial_respects_stencil_override():
    """poly_partial must honour set_default_stencil overrides (regression for the
    _default_stencil vs get_default_stencil bug)."""
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    k = 3.0
    f = np.cos(k * X)

    # Compute with the small override (3, 2) and large override (9, 5) and
    # verify that the results differ — proving poly_partial actually uses the
    # override rather than hardcoded defaults.
    pd.set_default_stencil(3, 2, order=1)
    result_small = pd.poly_partial(f, Lx, Ly, 1, 0)

    pd.set_default_stencil(9, 5, order=1)
    result_large = pd.poly_partial(f, Lx, Ly, 1, 0)

    pd.set_default_stencil(None, None, order=1)  # restore

    true_fx = -k * np.sin(k * X)
    err_small = np.max(np.abs(result_small - true_fx))
    err_large = np.max(np.abs(result_large - true_fx))

    # Large stencil must be substantially more accurate.
    assert err_large < err_small * 0.01, (
        f"Large stencil err ({err_large:.2e}) should be <1% of small stencil err "
        f"({err_small:.2e}); override is probably not being applied."
    )
    # And the large stencil should be highly accurate.
    np.testing.assert_allclose(result_large, true_fx, rtol=1e-4, atol=1e-6)


@pytest.mark.parametrize(
    "k,op,true_fn",
    [
        (3, "laplacian", lambda k, X: -(k**2) * np.cos(k * X)),
        (5, "laplacian", lambda k, X: -(k**2) * np.cos(k * X)),
        (3, "bilaplacian", lambda k, X: k**4 * np.cos(k * X)),
    ],
)
def test_default_stencil_accuracy_moderate_wavenumber(k, op, true_fn):
    """Default stencils (updated to deg=9,r=5 / deg=11,r=6) should give <0.1%
    relative error for wavenumbers up to k=5 on a 64-point grid."""
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    f = np.cos(k * X)
    true = true_fn(k, X)
    scale = np.max(np.abs(true))

    if op == "laplacian":
        result = pd.polynomial_laplacian(f, Lx, Ly)
    else:
        result = pd.polynomial_bilaplacian(f, Lx, Ly)

    rel_err = np.max(np.abs(result - true)) / scale
    assert rel_err < 1e-3, (
        f"Default stencil {op} at k={k}: relative error {rel_err:.2e} exceeds 0.1%."
    )


def test_poly_partial_and_named_function_agree_after_override():
    """After a stencil override, poly_partial and polynomial_laplacian must
    both use the overridden (deg, r) and produce identical results."""
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    k = 3.0
    f = np.cos(k * X)

    pd.set_default_stencil(7, 4, order=2)

    lap_named = pd.polynomial_laplacian(f, Lx, Ly)
    lap_generic = pd.poly_partial(f, Lx, Ly, 2, 0) + pd.poly_partial(f, Lx, Ly, 0, 2)

    pd.set_default_stencil(None, None, order=2)  # restore

    np.testing.assert_allclose(lap_named, lap_generic, rtol=1e-12, atol=1e-12)


# ---------------------------------------------------------------------------
# Additional robustness tests
# ---------------------------------------------------------------------------


def test_polynomial_quadratic_exactness_interior():
    """Degree ≥ 2 polynomial fitting recovers quadratic fields exactly in the interior.

    f = a·x² + b·x·y + c·y² has constant second-order derivatives: f_xx = 2a,
    f_xy = b, f_yy = 2c.  The fit is only exact away from the domain boundary
    because the stencil wraps periodically (non-periodic polynomial creates a jump).
    """
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    a, b, c = 1.5, -0.7, 2.3
    f = a * X ** 2 + b * X * Y + c * Y ** 2

    _, r = pd.get_default_stencil(2)
    fxx, fxy, fyy = pd.polynomial_hessian(f, Lx, Ly)

    np.testing.assert_allclose(fxx[r:-r, r:-r], 2 * a, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(fxy[r:-r, r:-r], b,     rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(fyy[r:-r, r:-r], 2 * c, rtol=1e-10, atol=1e-10)


def test_polynomial_linearity():
    """Polynomial differentiation satisfies d(a·f + b·g) = a·d(f) + b·d(g).

    The kernel is a fixed linear filter, so linearity must hold to machine
    precision regardless of the polynomial approximation error.
    """
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx1 = 2.0 * math.pi * 1 / Lx
    ky1 = 2.0 * math.pi * 1 / Ly
    kx2 = 2.0 * math.pi * 2 / Lx
    ky2 = 2.0 * math.pi * 1 / Ly
    f = np.cos(kx1 * X + ky1 * Y)
    g = np.sin(kx2 * X + ky2 * Y)
    a, b = 3.0, -2.0
    h = a * f + b * g

    fx, fy = pd.polynomial_grad(f, Lx, Ly)
    gx, gy = pd.polynomial_grad(g, Lx, Ly)
    hx, hy = pd.polynomial_grad(h, Lx, Ly)
    # atol=1e-12: floating-point accumulation in the kernel correlation
    # (sum of weighted np.roll outputs) introduces ~1e-14 round-off per
    # element; cross-architecture variance can push a few elements to
    # ~1.7e-13, so we use atol=1e-12 for safety.
    np.testing.assert_allclose(hx, a * fx + b * gx, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(hy, a * fy + b * gy, rtol=1e-12, atol=1e-12)

    lap_f = pd.polynomial_laplacian(f, Lx, Ly)
    lap_g = pd.polynomial_laplacian(g, Lx, Ly)
    lap_h = pd.polynomial_laplacian(h, Lx, Ly)
    np.testing.assert_allclose(lap_h, a * lap_f + b * lap_g, rtol=1e-12, atol=1e-12)


def test_polynomial_hessian_trace_equals_laplacian():
    """f_xx + f_yy from polynomial_hessian must equal polynomial_laplacian exactly.

    Both functions call _deriv(f, 2, 0) + _deriv(f, 0, 2) with the same (deg, r)
    so the identity is exact, not just approximate.
    """
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * 2 / Lx
    ky = 2.0 * math.pi * 3 / Ly
    f = np.cos(kx * X + ky * Y)

    fxx, _fxy, fyy = pd.polynomial_hessian(f, Lx, Ly)
    lap = pd.polynomial_laplacian(f, Lx, Ly)
    np.testing.assert_array_equal(fxx + fyy, lap)


def test_polynomial_third_order_approximate():
    """polynomial_third_order gives correct 3rd-order partials on a 2D cosine mode."""
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 128, 128
    mx, my = 1, 1
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    theta = kx * X + ky * Y
    f = np.cos(theta)

    # ∂^3 cos(kx·x + ky·y) by repeated chain rule
    fxxx_true =  kx ** 3 * np.sin(theta)
    fxxy_true =  kx ** 2 * ky * np.sin(theta)
    fxyy_true =  kx * ky ** 2 * np.sin(theta)
    fyyy_true =  ky ** 3 * np.sin(theta)

    fxxx, fxxy, fxyy, fyyy = pd.polynomial_third_order(f, Lx, Ly)
    np.testing.assert_allclose(fxxx, fxxx_true, rtol=5e-3, atol=1e-6)
    np.testing.assert_allclose(fxxy, fxxy_true, rtol=5e-3, atol=1e-6)
    np.testing.assert_allclose(fxyy, fxyy_true, rtol=5e-3, atol=1e-6)
    np.testing.assert_allclose(fyyy, fyyy_true, rtol=5e-3, atol=1e-6)


def test_polynomial_third_order_matches_poly_partial():
    """polynomial_third_order components match poly_partial for each multi-index."""
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * 2 / Lx
    ky = 2.0 * math.pi * 1 / Ly
    f = np.cos(kx * X + ky * Y)

    fxxx, fxxy, fxyy, fyyy = pd.polynomial_third_order(f, Lx, Ly)

    np.testing.assert_array_equal(fxxx, pd.poly_partial(f, Lx, Ly, 3, 0))
    np.testing.assert_array_equal(fxxy, pd.poly_partial(f, Lx, Ly, 2, 1))
    np.testing.assert_array_equal(fxyy, pd.poly_partial(f, Lx, Ly, 1, 2))
    np.testing.assert_array_equal(fyyy, pd.poly_partial(f, Lx, Ly, 0, 3))


def test_polynomial_bilaplacian_2d_diagonal_mode():
    """polynomial_bilaplacian correctly handles a non-axis-aligned mode.

    A mode with kx ≠ 0 and ky ≠ 0 exercises the cross-derivative term 2·f_xxyy
    in the bilaplacian expansion.  The analytic result is (kx²+ky²)² cos(θ).
    """
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 128, 128
    mx, my = 1, 2  # non-axis-aligned
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    theta = kx * X + ky * Y
    f = np.cos(theta)

    bilap_true = (kx ** 2 + ky ** 2) ** 2 * np.cos(theta)
    bilap = pd.polynomial_bilaplacian(f, Lx, Ly)

    scale = np.max(np.abs(bilap_true))
    rel_err = np.max(np.abs(bilap - bilap_true)) / scale
    assert rel_err < 2e-2, (
        f"Bilaplacian 2D diagonal mode: relative error {rel_err:.2e} exceeds 2 %."
    )


def test_polynomial_trilaplacian_approximate():
    """polynomial_trilaplacian gives reasonable accuracy on a low-k 2D cosine.

    Δ³f = −(kx²+ky²)³ f for f = cos(kx·x+ky·y).  We use the lowest wavenumber
    (mx=my=1) to minimise the local-fit truncation error on 6th-order operators.
    """
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 128, 128
    mx, my = 1, 1
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    theta = kx * X + ky * Y
    f = np.cos(theta)

    trilap_true = -(kx ** 2 + ky ** 2) ** 3 * np.cos(theta)
    trilap = pd.polynomial_trilaplacian(f, Lx, Ly)

    scale = np.max(np.abs(trilap_true))
    rel_err = np.max(np.abs(trilap - trilap_true)) / scale
    assert rel_err < 5e-2, (
        f"Trilaplacian 2D mode: relative error {rel_err:.2e} exceeds 5 %."
    )


def test_polynomial_prefilter_changes_result():
    """prefilter_sigma != 0 must change the polynomial gradient output.

    For a high-k mode the Gaussian filter exp(−½σ²|k|²) is far from 1, so the
    filtered and unfiltered results must differ substantially.
    """
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * 4 / Lx
    ky = 2.0 * math.pi * 4 / Ly
    f = np.cos(kx * X + ky * Y)

    fx_plain, _ = pd.polynomial_grad(f, Lx, Ly)
    fx_filtered, _ = pd.polynomial_grad(f, Lx, Ly, prefilter_sigma=0.5)

    assert not np.allclose(fx_plain, fx_filtered, rtol=1e-2), (
        "prefilter_sigma=0.5 should substantially change the gradient output "
        "for a high-k mode."
    )


def test_polynomial_stencil_restored_after_clear():
    """After clearing a per-order stencil override, results return to default.

    Sets a deliberately small stencil (deg=3, r=1) that gives a measurably
    different result, then clears it and verifies the output is bit-for-bit
    identical to the pre-override baseline.
    """
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    f = np.cos(2.0 * math.pi * 2 / Lx * X)

    lap_default = pd.polynomial_laplacian(f, Lx, Ly)

    pd.set_default_stencil(3, 1, order=2)
    lap_overridden = pd.polynomial_laplacian(f, Lx, Ly)

    pd.set_default_stencil(None, None, order=2)
    lap_restored = pd.polynomial_laplacian(f, Lx, Ly)

    assert not np.allclose(lap_overridden, lap_default), (
        "The small stencil override should have changed the Laplacian."
    )
    np.testing.assert_array_equal(lap_restored, lap_default)


def test_spectral_polynomial_backends_agree_low_frequency():
    """Spectral and polynomial backends agree to < 0.1 % on low-frequency data.

    Import spectral_derivatives inline so the test file does not acquire a
    hard module-level dependency on the spectral backend.
    """
    from equivariant_pde_find import spectral_derivatives as sd

    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 128, 128
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    mx, my = 1, 1
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    f = np.cos(kx * X + ky * Y)

    # Gradient
    fx_s, fy_s = sd.spectral_grad(f, Lx, Ly)
    fx_p, fy_p = pd.polynomial_grad(f, Lx, Ly)
    np.testing.assert_allclose(fx_p, fx_s, rtol=1e-3, atol=1e-6)
    np.testing.assert_allclose(fy_p, fy_s, rtol=1e-3, atol=1e-6)

    # Laplacian
    lap_s = sd.spectral_laplacian(f, Lx, Ly)
    lap_p = pd.polynomial_laplacian(f, Lx, Ly)
    np.testing.assert_allclose(lap_p, lap_s, rtol=1e-3, atol=1e-6)

    # Bilaplacian
    bilap_s = sd.spectral_bilaplacian(f, Lx, Ly)
    bilap_p = pd.polynomial_bilaplacian(f, Lx, Ly)
    np.testing.assert_allclose(bilap_p, bilap_s, rtol=1e-3, atol=1e-6)
