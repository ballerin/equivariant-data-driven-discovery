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
CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, CODE_DIR)

import polynomial_derivatives as pd  # noqa: E402


def make_grid(Nx, Ny, Lx, Ly):
    x = np.linspace(0.0, Lx, Nx, endpoint=False)
    y = np.linspace(0.0, Ly, Ny, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="ij")
    return X, Y


def test_polynomial_grad_linear_exactness():
    """Degree-1 polynomials should be recovered exactly by the local fit."""
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 32, 32
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    a, b, c = 1.234, -0.987, 0.5
    f = a * X + b * Y + c

    fx_p, fy_p = pd.polynomial_grad(f, Lx, Ly)
    # The polynomial backend uses periodic wrapping when building local
    # stencils. A globally non-periodic polynomial creates a jump at the
    # domain boundary which contaminates fits near the edges. Restrict the
    # exactness check to interior points where the stencil does not wrap.
    r = 2  # default stencil half-width for order-1 (deg=3, r=2)
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
