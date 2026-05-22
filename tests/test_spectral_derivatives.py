"""
Tests for `spectral_derivatives.py`.

These exercises verify that the FFT-based derivative operators produce the
analytic derivatives for band-limited trigonometric polynomials on periodic
rectangles [0, Lx) × [0, Ly).  Tests include different domain sizes and
rectangular grids, and verify behaviour for complex-valued inputs.

Notes
-----
- Tolerances are tight (near machine precision) for spectral operators on
  band-limited inputs; very-high-order operators use slightly relaxed
  absolute tolerances due to floating-point accumulation.
"""

import math
import os
import sys

import numpy as np
import pytest

# Ensure the code directory is importable when tests are run from the repo root.
CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, CODE_DIR)

import spectral_derivatives as sd  # noqa: E402


def make_grid(Nx, Ny, Lx, Ly):
    x = np.linspace(0.0, Lx, Nx, endpoint=False)
    y = np.linspace(0.0, Ly, Ny, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="ij")
    return X, Y


@pytest.mark.parametrize(
    "Lx,Ly,Nx,Ny,mx,my",
    [
        (2 * math.pi, 2 * math.pi, 32, 24, 3, 4),
        (8 * math.pi, 8 * math.pi, 32, 24, 3, 4),
    ],
)
def test_spectral_grad_and_laplacian_exact(Lx, Ly, Nx, Ny, mx, my):
    """Gradient and Laplacian are exact for a single cosine mode."""
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    theta = kx * X + ky * Y
    f = np.cos(theta)

    fx_true = -kx * np.sin(theta)
    fy_true = -ky * np.sin(theta)

    fx, fy = sd.spectral_grad(f, Lx, Ly)
    np.testing.assert_allclose(fx, fx_true, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(fy, fy_true, rtol=1e-12, atol=1e-12)

    lap_true = -(kx**2 + ky**2) * np.cos(theta)
    lap = sd.spectral_laplacian(f, Lx, Ly)
    np.testing.assert_allclose(lap, lap_true, rtol=1e-12, atol=1e-12)


def test_spectral_hessian_bilap_trilap_and_third_order():
    """Higher-order operators match analytic derivatives on a single mode."""
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    mx, my = 7, 5
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    theta = kx * X + ky * Y
    f = np.cos(theta)

    fxx_true = -(kx**2) * np.cos(theta)
    fxy_true = -kx * ky * np.cos(theta)
    fyy_true = -(ky**2) * np.cos(theta)
    fxx, fxy, fyy = sd.spectral_hessian(f, Lx, Ly)
    np.testing.assert_allclose(fxx, fxx_true, rtol=1e-11, atol=1e-11)
    np.testing.assert_allclose(fxy, fxy_true, rtol=1e-11, atol=1e-11)
    np.testing.assert_allclose(fyy, fyy_true, rtol=1e-11, atol=1e-11)

    bilap_true = (kx**2 + ky**2) ** 2 * np.cos(theta)
    np.testing.assert_allclose(
        sd.spectral_bilaplacian(f, Lx, Ly), bilap_true, rtol=1e-10, atol=1e-8
    )

    trilap_true = -((kx**2 + ky**2) ** 3) * np.cos(theta)
    np.testing.assert_allclose(
        sd.spectral_trilaplacian(f, Lx, Ly), trilap_true, rtol=1e-10, atol=1e-5
    )

    fxxx_true = kx**3 * np.sin(theta)
    fxxy_true = kx**2 * ky * np.sin(theta)
    fxyy_true = kx * ky**2 * np.sin(theta)
    fyyy_true = ky**3 * np.sin(theta)
    fxxx, fxxy, fxyy, fyyy = sd.spectral_third_order(f, Lx, Ly)
    np.testing.assert_allclose(fxxx, fxxx_true, rtol=1e-11, atol=1e-10)
    np.testing.assert_allclose(fxxy, fxxy_true, rtol=1e-11, atol=1e-10)
    np.testing.assert_allclose(fxyy, fxyy_true, rtol=1e-11, atol=1e-10)
    np.testing.assert_allclose(fyyy, fyyy_true, rtol=1e-11, atol=1e-10)


def test_rectangular_domain_scaling():
    """Check correctness on a rectangular domain (Nx != Ny, Lx != Ly)."""
    Lx, Ly, Nx, Ny = 4 * math.pi, 2 * math.pi, 64, 32
    mx, my = 5, 3
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    theta = kx * X + ky * Y
    f = np.cos(theta)

    fx_true = -kx * np.sin(theta)
    fy_true = -ky * np.sin(theta)
    fx, fy = sd.spectral_grad(f, Lx, Ly)
    np.testing.assert_allclose(fx, fx_true, rtol=1e-11, atol=1e-11)
    np.testing.assert_allclose(fy, fy_true, rtol=1e-11, atol=1e-11)

    lap_true = -(kx**2 + ky**2) * np.cos(theta)
    lap = sd.spectral_laplacian(f, Lx, Ly)
    np.testing.assert_allclose(lap, lap_true, rtol=1e-11, atol=1e-11)


def test_complex_input_imag_discarded():
    """If the input has an imaginary part, the routines return the real part
    of the derivative (module uses ``np.real`` on IFFT results).
    """
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 32, 32
    mx, my = 3, 2
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    theta = kx * X + ky * Y
    f = np.cos(theta)
    f_complex = f + 0.5j * f

    fx_real, fy_real = sd.spectral_grad(f, Lx, Ly)
    fx_c, fy_c = sd.spectral_grad(f_complex, Lx, Ly)
    np.testing.assert_allclose(fx_c, fx_real, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(fy_c, fy_real, rtol=1e-12, atol=1e-12)


def test_input_shape_error():
    """Passing a 1-D array should raise a ValueError due to shape unpacking."""
    arr = np.arange(16)
    with pytest.raises(ValueError):
        sd.spectral_grad(arr, 2 * math.pi, 2 * math.pi)
