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
CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.derivatives import spectral as sd  # noqa: E402


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


def test_1d_input_is_valid_ndim_1():
    """A 1-D array is valid input (ndim=1) since the N-D generalization —
    spectral_grad returns a 1-tuple, exact for a single cosine mode."""
    N, L, m = 16, 2 * math.pi, 3
    x = np.linspace(0, L, N, endpoint=False)
    k = 2.0 * math.pi * m / L
    f = np.cos(k * x)
    (fx,) = sd.spectral_grad(f, L, L)
    np.testing.assert_allclose(fx, -k * np.sin(k * x), rtol=1e-12, atol=1e-12)


# ---------------------------------------------------------------------------
# Additional robustness tests
# ---------------------------------------------------------------------------


def test_sine_mode_gradient_and_laplacian():
    """Gradient and Laplacian are exact on a single sine mode (complements cosine tests)."""
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 32, 32
    mx, my = 4, 3
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    theta = kx * X + ky * Y
    f = np.sin(theta)

    fx_true = kx * np.cos(theta)
    fy_true = ky * np.cos(theta)
    fx, fy = sd.spectral_grad(f, Lx, Ly)
    np.testing.assert_allclose(fx, fx_true, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(fy, fy_true, rtol=1e-12, atol=1e-12)

    lap_true = -(kx ** 2 + ky ** 2) * np.sin(theta)
    np.testing.assert_allclose(
        sd.spectral_laplacian(f, Lx, Ly), lap_true, rtol=1e-12, atol=1e-12
    )


def test_hessian_trace_equals_laplacian():
    """f_xx + f_yy from spectral_hessian must equal spectral_laplacian exactly.

    Both compute IFFT(−|k|² f̂) split by component, so the identity is exact
    up to floating-point round-off from the independent IFFT calls.
    """
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 32, 32
    mx, my = 4, 3
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    f = np.cos(kx * X + ky * Y)

    fxx, _fxy, fyy = sd.spectral_hessian(f, Lx, Ly)
    lap = sd.spectral_laplacian(f, Lx, Ly)
    # Same multiplier split over two IFFT calls → machine precision agreement.
    np.testing.assert_allclose(fxx + fyy, lap, rtol=1e-12, atol=1e-12)


def test_bilaplacian_equals_iterated_laplacian():
    """Δ²f = Δ(Δf): applying the spectral Laplacian twice must equal bilaplacian.

    For band-limited input, each FFT→IFFT round-trip is exact, so the two-step
    composition and the direct bilaplacian should agree to near machine precision.
    """
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    mx, my = 3, 5
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    f = np.cos(kx * X + ky * Y)

    bilap_direct = sd.spectral_bilaplacian(f, Lx, Ly)
    bilap_iterated = sd.spectral_laplacian(sd.spectral_laplacian(f, Lx, Ly), Lx, Ly)
    np.testing.assert_allclose(bilap_iterated, bilap_direct, rtol=1e-10, atol=1e-8)


def test_linearity_gradient_and_laplacian():
    """∂(a·f + b·g) = a·∂f + b·∂g for gradient and Laplacian (superposition)."""
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 32, 32
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx1 = 2.0 * math.pi * 2 / Lx
    ky1 = 2.0 * math.pi * 3 / Ly
    kx2 = 2.0 * math.pi * 5 / Lx
    ky2 = 2.0 * math.pi * 1 / Ly
    f = np.cos(kx1 * X + ky1 * Y)
    g = np.sin(kx2 * X + ky2 * Y)
    a, b = 2.5, -1.3
    h = a * f + b * g

    fx, fy = sd.spectral_grad(f, Lx, Ly)
    gx, gy = sd.spectral_grad(g, Lx, Ly)
    hx, hy = sd.spectral_grad(h, Lx, Ly)
    np.testing.assert_allclose(hx, a * fx + b * gx, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(hy, a * fy + b * gy, rtol=1e-12, atol=1e-12)

    lap_f = sd.spectral_laplacian(f, Lx, Ly)
    lap_g = sd.spectral_laplacian(g, Lx, Ly)
    lap_h = sd.spectral_laplacian(h, Lx, Ly)
    np.testing.assert_allclose(lap_h, a * lap_f + b * lap_g, rtol=1e-12, atol=1e-12)


def test_zero_mean_of_laplacian():
    """For a zero-mean periodic field, ∫Δf dx dy = 0 (Gauss / periodicity).

    The spectral multiplier at k=(0,0) is zero, so the Laplacian of any field
    always has zero spatial mean — regardless of the mean of f itself.
    """
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    # Non-trivial field with a non-zero mean constant offset included.
    f = np.cos(3 * X + 2 * Y) + 0.5 * np.sin(5 * X - Y) + 1.7
    lap = sd.spectral_laplacian(f, Lx, Ly)
    np.testing.assert_allclose(np.mean(lap), 0.0, atol=1e-12)


def test_output_shapes_preserved():
    """Every spectral operator returns an array with the same shape as the input."""
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 32, 48  # rectangular on purpose
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    f = np.cos(2 * X + 3 * Y)

    assert sd.spectral_laplacian(f, Lx, Ly).shape == (Nx, Ny)
    assert sd.spectral_bilaplacian(f, Lx, Ly).shape == (Nx, Ny)
    assert sd.spectral_trilaplacian(f, Lx, Ly).shape == (Nx, Ny)
    for arr in sd.spectral_grad(f, Lx, Ly):
        assert arr.shape == (Nx, Ny)
    for arr in sd.spectral_hessian(f, Lx, Ly):
        assert arr.shape == (Nx, Ny)
    for arr in sd.spectral_third_order(f, Lx, Ly):
        assert arr.shape == (Nx, Ny)


def test_prefilter_sigma_attenuates_amplitude():
    """prefilter_sigma=σ suppresses a single-mode field by ≈exp(−½ σ² |k|²).

    The Gaussian filter in Fourier space is exact (it is implemented as a
    spectral multiplier), so the amplitude ratio should match to high precision.
    """
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    mx, my = 5, 5
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    f = np.cos(kx * X + ky * Y)

    # Choose sigma so attenuation is modest (~10 %) — large enough to measure
    # accurately but well above the floating-point noise floor.
    sigma = 0.3  # exp(-0.5 * 0.09 * k²) ≈ 0.105 for k² = (5·2π/L)² × 2 = 50
    k2 = kx ** 2 + ky ** 2
    expected_attenuation = math.exp(-0.5 * sigma ** 2 * k2)

    lap_plain = sd.spectral_laplacian(f, Lx, Ly)
    lap_filtered = sd.spectral_laplacian(f, Lx, Ly, prefilter_sigma=sigma)

    # Ratio of peak amplitudes should match the Gaussian weight.
    ratio = np.max(np.abs(lap_filtered)) / np.max(np.abs(lap_plain))
    np.testing.assert_allclose(ratio, expected_attenuation, rtol=1e-6)


def test_third_order_is_gradient_of_hessian():
    """f_xxy from spectral_third_order equals ∂_y(f_xx) via iterated spectral ops.

    Both paths reduce to the same Fourier multiplier (−i k_x² k_y) applied to
    the common FFT of f, so they agree to near machine precision for band-limited
    input despite involving different numbers of FFT round-trips.
    """
    Lx, Ly, Nx, Ny = 2 * math.pi, 2 * math.pi, 64, 64
    mx, my = 3, 4
    X, Y = make_grid(Nx, Ny, Lx, Ly)
    kx = 2.0 * math.pi * mx / Lx
    ky = 2.0 * math.pi * my / Ly
    f = np.cos(kx * X + ky * Y)

    _fxxx, fxxy, fxyy, _fyyy = sd.spectral_third_order(f, Lx, Ly)

    fxx, _fxy, fyy = sd.spectral_hessian(f, Lx, Ly)
    # ∂_y(f_xx) should equal f_xxy; ∂_x(f_yy) should equal f_xyy.
    d_fxx_dy = sd.spectral_grad(fxx, Lx, Ly)[1]
    d_fyy_dx = sd.spectral_grad(fyy, Lx, Ly)[0]

    np.testing.assert_allclose(fxxy, d_fxx_dy, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(fxyy, d_fyy_dx, rtol=1e-10, atol=1e-10)


# =============================================================================
# 3-D spectral derivative tests
# =============================================================================


def make_grid_3d(Nx, Ny, Nz, Lx, Ly, Lz):
    """Build 3-D coordinate grids."""
    x = np.linspace(0.0, Lx, Nx, endpoint=False)
    y = np.linspace(0.0, Ly, Ny, endpoint=False)
    z = np.linspace(0.0, Lz, Nz, endpoint=False)
    return np.meshgrid(x, y, z, indexing="ij")


class TestSpectralDerivatives3D:
    """3-D spectral derivative correctness on trigonometric polynomials."""

    def test_grad_3d_exact(self):
        """Gradient on a pure cosine mode in 3-D."""
        Nx, Ny, Nz = 16, 16, 16
        Lx, Ly, Lz = 2 * math.pi, 2 * math.pi, 2 * math.pi
        mx, my, mz = 2, 3, 1
        X, Y, Z = make_grid_3d(Nx, Ny, Nz, Lx, Ly, Lz)
        kx = 2 * math.pi * mx / Lx
        ky = 2 * math.pi * my / Ly
        kz = 2 * math.pi * mz / Lz
        theta = kx * X + ky * Y + kz * Z
        f = np.cos(theta)

        fx, fy, fz = sd.spectral_grad(f, Lx, Ly, Lz)
        np.testing.assert_allclose(fx, -kx * np.sin(theta), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(fy, -ky * np.sin(theta), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(fz, -kz * np.sin(theta), rtol=1e-12, atol=1e-12)

    def test_laplacian_3d_exact(self):
        """Laplacian on a pure cosine mode in 3-D."""
        Nx, Ny, Nz = 16, 16, 16
        Lx, Ly, Lz = 2 * math.pi, 2 * math.pi, 2 * math.pi
        mx, my, mz = 2, 3, 1
        X, Y, Z = make_grid_3d(Nx, Ny, Nz, Lx, Ly, Lz)
        kx = 2 * math.pi * mx / Lx
        ky = 2 * math.pi * my / Ly
        kz = 2 * math.pi * mz / Lz
        theta = kx * X + ky * Y + kz * Z
        f = np.cos(theta)

        lap = sd.spectral_laplacian(f, Lx, Ly, Lz)
        lap_true = -(kx**2 + ky**2 + kz**2) * np.cos(theta)
        np.testing.assert_allclose(lap, lap_true, rtol=1e-12, atol=1e-12)

    def test_hessian_3d_exact(self):
        """Hessian on a pure cosine mode in 3-D — 6 components."""
        Nx, Ny, Nz = 16, 16, 16
        Lx, Ly, Lz = 2 * math.pi, 2 * math.pi, 2 * math.pi
        mx, my, mz = 2, 3, 1
        X, Y, Z = make_grid_3d(Nx, Ny, Nz, Lx, Ly, Lz)
        kx = 2 * math.pi * mx / Lx
        ky = 2 * math.pi * my / Ly
        kz = 2 * math.pi * mz / Lz
        theta = kx * X + ky * Y + kz * Z
        f = np.cos(theta)

        fxx, fxy, fxz, fyy, fyz, fzz = sd.spectral_hessian(f, Lx, Ly, Lz)

        np.testing.assert_allclose(fxx, -(kx**2) * np.cos(theta), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(fxy, -kx * ky * np.cos(theta), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(fxz, -kx * kz * np.cos(theta), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(fyy, -(ky**2) * np.cos(theta), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(fyz, -ky * kz * np.cos(theta), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(fzz, -(kz**2) * np.cos(theta), rtol=1e-12, atol=1e-12)

    def test_bilaplacian_3d_exact(self):
        """Bilaplacian on a pure cosine mode in 3-D."""
        Nx, Ny, Nz = 16, 16, 16
        Lx, Ly, Lz = 2 * math.pi, 2 * math.pi, 2 * math.pi
        mx, my, mz = 2, 3, 1
        X, Y, Z = make_grid_3d(Nx, Ny, Nz, Lx, Ly, Lz)
        kx = 2 * math.pi * mx / Lx
        ky = 2 * math.pi * my / Ly
        kz = 2 * math.pi * mz / Lz
        theta = kx * X + ky * Y + kz * Z
        f = np.cos(theta)

        bilap = sd.spectral_bilaplacian(f, Lx, Ly, Lz)
        bilap_true = (kx**2 + ky**2 + kz**2) ** 2 * np.cos(theta)
        np.testing.assert_allclose(bilap, bilap_true, rtol=1e-10, atol=1e-8)

    def test_trilaplacian_3d_exact(self):
        """Trilaplacian on a pure cosine mode in 3-D."""
        Nx, Ny, Nz = 16, 16, 16
        Lx, Ly, Lz = 2 * math.pi, 2 * math.pi, 2 * math.pi
        mx, my, mz = 2, 3, 1
        X, Y, Z = make_grid_3d(Nx, Ny, Nz, Lx, Ly, Lz)
        kx = 2 * math.pi * mx / Lx
        ky = 2 * math.pi * my / Ly
        kz = 2 * math.pi * mz / Lz
        theta = kx * X + ky * Y + kz * Z
        f = np.cos(theta)

        trilap = sd.spectral_trilaplacian(f, Lx, Ly, Lz)
        trilap_true = -((kx**2 + ky**2 + kz**2) ** 3) * np.cos(theta)
        np.testing.assert_allclose(trilap, trilap_true, rtol=1e-10, atol=1e-5)

    def test_third_order_3d_exact(self):
        """Third-order partials on a pure cosine mode in 3-D — 10 components."""
        Nx, Ny, Nz = 16, 16, 16
        Lx, Ly, Lz = 2 * math.pi, 2 * math.pi, 2 * math.pi
        mx, my, mz = 2, 3, 1
        X, Y, Z = make_grid_3d(Nx, Ny, Nz, Lx, Ly, Lz)
        kx = 2 * math.pi * mx / Lx
        ky = 2 * math.pi * my / Ly
        kz = 2 * math.pi * mz / Lz
        theta = kx * X + ky * Y + kz * Z
        f = np.cos(theta)
        # ∂^3 cos(θ) = sin(θ) · (+k_i k_j k_l)  (i³ = −i, −i·−i = −1, wait: i³=−i,
        #  cos→sin has sign −1, so ∂^3 cos = +sin · product)
        # For cos(k·x): ∂_i cos = −k_i sin, ∂_ij cos = −k_i k_j cos,
        # ∂_ijk cos = +k_i k_j k_k sin
        sin_theta = np.sin(theta)

        derivs = sd.spectral_third_order(f, Lx, Ly, Lz)
        names = ["f_xxx", "f_xxy", "f_xxz", "f_xyy", "f_xyz", "f_xzz",
                 "f_yyy", "f_yyz", "f_yzz", "f_zzz"]
        expected_multipliers = [
            kx**3, kx**2 * ky, kx**2 * kz,
            kx * ky**2, kx * ky * kz, kx * kz**2,
            ky**3, ky**2 * kz, ky * kz**2,
            kz**3,
        ]

        assert len(derivs) == 10, f"Expected 10 third-order components, got {len(derivs)}"
        for i, (deriv, mult) in enumerate(zip(derivs, expected_multipliers)):
            expected = mult * sin_theta
            np.testing.assert_allclose(deriv, expected, rtol=1e-11, atol=1e-10,
                                       err_msg=f"Component {names[i]} mismatch")

    def test_hessian_trace_equals_laplacian_3d(self):
        """f_xx + f_yy + f_zz must equal spectral_laplacian."""
        Nx, Ny, Nz = 16, 16, 16
        Lx, Ly, Lz = 2 * math.pi, 2 * math.pi, 2 * math.pi
        X, Y, Z = make_grid_3d(Nx, Ny, Nz, Lx, Ly, Lz)
        f = np.cos(2 * X + 3 * Y + Z)

        fxx, fxy, fxz, fyy, fyz, fzz = sd.spectral_hessian(f, Lx, Ly, Lz)
        lap = sd.spectral_laplacian(f, Lx, Ly, Lz)
        np.testing.assert_allclose(fxx + fyy + fzz, lap, rtol=1e-12, atol=1e-12)

    def test_output_shapes_preserved_3d(self):
        """Every spectral operator returns arrays with the same shape as input."""
        Nx, Ny, Nz = 16, 24, 20  # non-cubic on purpose
        Lx, Ly, Lz = 2 * math.pi, 3 * math.pi, 4 * math.pi
        X, Y, Z = make_grid_3d(Nx, Ny, Nz, Lx, Ly, Lz)
        f = np.cos(2 * X + 3 * Y + Z)

        assert sd.spectral_laplacian(f, Lx, Ly, Lz).shape == (Nx, Ny, Nz)
        assert sd.spectral_bilaplacian(f, Lx, Ly, Lz).shape == (Nx, Ny, Nz)
        assert sd.spectral_trilaplacian(f, Lx, Ly, Lz).shape == (Nx, Ny, Nz)
        for arr in sd.spectral_grad(f, Lx, Ly, Lz):
            assert arr.shape == (Nx, Ny, Nz)
        for arr in sd.spectral_hessian(f, Lx, Ly, Lz):
            assert arr.shape == (Nx, Ny, Nz)
        for arr in sd.spectral_third_order(f, Lx, Ly, Lz):
            assert arr.shape == (Nx, Ny, Nz)

    def test_rectangular_3d_domain(self):
        """Check correctness on a fully non-cubic domain (Nx≠Ny≠Nz, Lx≠Ly≠Lz)."""
        Nx, Ny, Nz = 20, 16, 12
        Lx, Ly, Lz = 4 * math.pi, 3 * math.pi, 2 * math.pi
        mx, my, mz = 3, 2, 4
        X, Y, Z = make_grid_3d(Nx, Ny, Nz, Lx, Ly, Lz)
        kx = 2 * math.pi * mx / Lx
        ky = 2 * math.pi * my / Ly
        kz = 2 * math.pi * mz / Lz
        theta = kx * X + ky * Y + kz * Z
        f = np.cos(theta)

        fx, fy, fz = sd.spectral_grad(f, Lx, Ly, Lz)
        np.testing.assert_allclose(fx, -kx * np.sin(theta), rtol=1e-11, atol=1e-11)
        np.testing.assert_allclose(fy, -ky * np.sin(theta), rtol=1e-11, atol=1e-11)
        np.testing.assert_allclose(fz, -kz * np.sin(theta), rtol=1e-11, atol=1e-11)

        lap = sd.spectral_laplacian(f, Lx, Ly, Lz)
        lap_true = -(kx**2 + ky**2 + kz**2) * np.cos(theta)
        np.testing.assert_allclose(lap, lap_true, rtol=1e-11, atol=1e-11)

    def test_bilaplacian_equals_iterated_laplacian_3d(self):
        """Δ²f = Δ(Δf) in 3-D."""
        Nx, Ny, Nz = 16, 16, 16
        Lx, Ly, Lz = 2 * math.pi, 2 * math.pi, 2 * math.pi
        X, Y, Z = make_grid_3d(Nx, Ny, Nz, Lx, Ly, Lz)
        f = np.cos(3 * X + 2 * Y + 4 * Z)

        bilap_direct = sd.spectral_bilaplacian(f, Lx, Ly, Lz)
        bilap_iter = sd.spectral_laplacian(
            sd.spectral_laplacian(f, Lx, Ly, Lz), Lx, Ly, Lz
        )
        np.testing.assert_allclose(bilap_iter, bilap_direct, rtol=1e-10, atol=1e-8)

    def test_linearity_3d(self):
        """∂(a·f + b·g) = a·∂f + b·∂g in 3-D."""
        Nx, Ny, Nz = 16, 16, 16
        Lx, Ly, Lz = 2 * math.pi, 2 * math.pi, 2 * math.pi
        X, Y, Z = make_grid_3d(Nx, Ny, Nz, Lx, Ly, Lz)
        f = np.cos(2 * X + 3 * Y + Z)
        g = np.sin(X - 2 * Y + 5 * Z)
        a, b = 2.5, -1.3
        h = a * f + b * g

        fx, fy, fz = sd.spectral_grad(f, Lx, Ly, Lz)
        gx, gy, gz = sd.spectral_grad(g, Lx, Ly, Lz)
        hx, hy, hz = sd.spectral_grad(h, Lx, Ly, Lz)
        np.testing.assert_allclose(hx, a * fx + b * gx, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(hy, a * fy + b * gy, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(hz, a * fz + b * gz, rtol=1e-12, atol=1e-12)

        lap_f = sd.spectral_laplacian(f, Lx, Ly, Lz)
        lap_g = sd.spectral_laplacian(g, Lx, Ly, Lz)
        lap_h = sd.spectral_laplacian(h, Lx, Ly, Lz)
        np.testing.assert_allclose(lap_h, a * lap_f + b * lap_g, rtol=1e-12, atol=1e-12)

    def test_4d_input_is_valid(self):
        """A 4-D array is valid input too (isotropic domain beyond ndim=3) —
        exact for a single cosine mode, matching the 2-D/3-D exactness tests."""
        N, L, ndim = 8, 2 * math.pi, 4  # N=8 -> Nyquist=4; keep ks well below it
        ks = np.array([2, 1, 3, 1])
        x = np.linspace(0, L, N, endpoint=False)
        grids = np.meshgrid(*([x] * ndim), indexing="ij")
        kvec = 2.0 * math.pi / L * ks
        theta = sum(k * g for k, g in zip(kvec, grids))
        f = np.cos(theta)

        grad = sd.spectral_grad(f, L, L)
        assert len(grad) == ndim
        for i in range(ndim):
            np.testing.assert_allclose(
                grad[i], -kvec[i] * np.sin(theta), rtol=1e-10, atol=1e-10
            )

        lap = sd.spectral_laplacian(f, L, L)
        np.testing.assert_allclose(
            lap, -np.sum(kvec**2) * f, rtol=1e-10, atol=1e-10
        )

    def test_polynomial_backend_supports_3d(self):
        """Polynomial backend now supports 3-D input (N-D generalization) —
        gradient of a single low-frequency cosine mode is approximately
        correct (polynomial fitting has O(dx^(deg+1-order)) bias, unlike
        the exact spectral backend)."""
        from equivariant_pde_find import derivatives
        N, L, m = 24, 2 * math.pi, 2
        x = np.linspace(0, L, N, endpoint=False)
        X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
        k = 2.0 * math.pi * m / L
        f = np.cos(k * X)
        with derivatives.use_backend("polynomial"):
            fx, fy, fz = derivatives.grad(f, L, L, L)
        np.testing.assert_allclose(fx, -k * np.sin(k * X), rtol=1e-2, atol=1e-2)
        np.testing.assert_allclose(fy, 0.0, atol=1e-2)
        np.testing.assert_allclose(fz, 0.0, atol=1e-2)


# ---------------------------------------------------------------------------
# Spectral-method conventions (Trefethen, "Spectral Methods in MATLAB", ch. 3)
#
# These pin two textbook properties that the N-D generalization relies on but
# that nothing previously asserted:
#
#   1. Nyquist handling.  For even N the mode k = N/2 is shared between +N/2
#      and -N/2 and must be zeroed for odd-order derivatives, kept for even
#      ones.  spectral.py achieves this implicitly (an odd-order multiplier
#      makes the real Nyquist coefficient purely imaginary, which the closing
#      np.real() discards) rather than with an explicit mask -- so a refactor
#      of that np.real could silently change the answer without these tests.
#
#   2. Exactness on trigonometric polynomials.  Pseudospectral
#      differentiation is machine-precision exact for band-limited periodic
#      data; this is the property that makes every analytic N-D test in this
#      suite meaningful.
# ---------------------------------------------------------------------------


class TestNyquistConvention:
    def _nyquist_field(self, ndim, N, L, axis=0):
        """cos(k_Nyquist · x_axis): the alternating +1/-1 mode."""
        x = np.linspace(0.0, L, N, endpoint=False)
        grids = np.meshgrid(*([x] * ndim), indexing="ij")
        k_nyq = (N // 2) * 2.0 * math.pi / L
        return np.cos(k_nyq * grids[axis]), k_nyq

    @pytest.mark.parametrize("ndim", [1, 2, 4])
    @pytest.mark.parametrize("order", [1, 3])
    def test_odd_order_derivative_of_nyquist_mode_is_zero(self, ndim, order):
        N, L = 16, 2 * math.pi
        f, _ = self._nyquist_field(ndim, N, L)
        multi_index = tuple(order if a == 0 else 0 for a in range(ndim))
        d = sd.spectral_partial_nd(f, (L,) * ndim, multi_index)
        np.testing.assert_allclose(d, 0.0, atol=1e-10)

    @pytest.mark.parametrize("ndim", [1, 2, 4])
    @pytest.mark.parametrize("order", [2, 4])
    def test_even_order_derivative_of_nyquist_mode_is_retained(self, ndim, order):
        N, L = 16, 2 * math.pi
        f, k_nyq = self._nyquist_field(ndim, N, L)
        multi_index = tuple(order if a == 0 else 0 for a in range(ndim))
        d = sd.spectral_partial_nd(f, (L,) * ndim, multi_index)
        expected = ((-1) ** (order // 2)) * k_nyq**order * f
        np.testing.assert_allclose(d, expected, rtol=1e-10, atol=1e-8)

    def test_nyquist_field_is_not_trivially_zero(self):
        """Guard: the tests above would pass vacuously on a zero field."""
        f, k_nyq = self._nyquist_field(2, 16, 2 * math.pi)
        assert np.allclose(np.abs(f), 1.0)
        assert k_nyq > 0


class TestPseudospectralExactnessND:
    """Machine-precision exactness on trigonometric polynomials at ndim=4."""

    def _trig_poly(self, grids, modes):
        """Σ_terms a · cos(m·x + φ) — band-limited, hence spectrally exact."""
        out = np.zeros_like(grids[0])
        for amp, ms, phase in modes:
            out = out + amp * np.cos(sum(m * g for m, g in zip(ms, grids)) + phase)
        return out

    def _analytic_partial(self, grids, modes, multi_index):
        """Exact ∂^α of the same trigonometric polynomial."""
        order = sum(multi_index)
        out = np.zeros_like(grids[0])
        for amp, ms, phase in modes:
            factor = 1.0
            for axis, a in enumerate(multi_index):
                factor *= ms[axis] ** a
            arg = sum(m * g for m, g in zip(ms, grids)) + phase
            # d^order/dt^order cos(t) cycles cos, -sin, -cos, sin
            wave = [np.cos(arg), -np.sin(arg), -np.cos(arg), np.sin(arg)][order % 4]
            out = out + amp * factor * wave
        return out

    @pytest.mark.parametrize("multi_index", [
        (1, 0, 0, 0), (0, 2, 0, 0), (1, 1, 0, 0),
        (2, 0, 1, 0), (2, 0, 0, 2), (1, 1, 1, 1),
    ])
    def test_exact_on_trig_polynomial_ndim4(self, multi_index):
        ndim, N, L = 4, 12, 2 * math.pi
        x = np.linspace(0.0, L, N, endpoint=False)
        grids = np.meshgrid(*([x] * ndim), indexing="ij")
        modes = [
            (1.0, (1, 2, 1, 0), 0.0),
            (0.5, (2, 0, 1, 3), 0.7),
            (-0.3, (0, 1, 2, 1), 1.9),
        ]
        f = self._trig_poly(grids, modes)
        got = sd.spectral_partial_nd(f, (L,) * ndim, multi_index)
        expected = self._analytic_partial(grids, modes, multi_index)
        scale = max(np.max(np.abs(expected)), 1.0)
        np.testing.assert_allclose(got, expected, atol=1e-10 * scale, rtol=0)

    def test_batch_matches_individual_calls_ndim4(self):
        """spectral_partials_nd_batch shares one FFT across many multi-indices
        -- it must agree exactly with computing them one at a time."""
        ndim, N, L = 4, 10, 2 * math.pi
        rng = np.random.default_rng(0)
        f = rng.standard_normal((N,) * ndim)
        indices = [(1, 0, 0, 0), (0, 0, 2, 0), (1, 1, 0, 0), (0, 0, 0, 0), (2, 0, 0, 2)]
        batch = sd.spectral_partials_nd_batch(f, (L,) * ndim, indices)
        for mi in indices:
            np.testing.assert_allclose(
                batch[mi], sd.spectral_partial_nd(f, (L,) * ndim, mi), rtol=1e-12
            )
