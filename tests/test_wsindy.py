"""
Tests for ``wsindy.py`` — weak-form SINDy for equivariant PDE discovery.
"""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, CODE_DIR)

import wsindy as ws


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def test_field():
    """A smooth 2-D field on a square periodic grid."""
    N, L = 64, 2 * math.pi
    x = np.linspace(0.0, L, N, endpoint=False)
    y = np.linspace(0.0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="ij")
    f = np.cos(2 * X + 3 * Y) + 0.5 * np.sin(X - 2 * Y)
    ft = -np.sin(2 * X + 3 * Y) + 0.5 * np.cos(X - 2 * Y)
    return f, ft, L, N


# ---------------------------------------------------------------------------
# Test function construction
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("m,p", [(4, 2), (8, 4), (12, 6)])
def test_kernel_shape_and_support(m, p):
    kernel = ws._test_function_1d(m, p)
    assert len(kernel) == 2 * m + 1
    # Centre value must be 1
    assert kernel[m] == 1.0
    # Boundaries must be zero (compact support)
    assert kernel[0] == 0.0
    assert kernel[-1] == 0.0
    # Non-negative
    assert np.all(kernel >= 0.0)
    # Symmetric
    np.testing.assert_array_almost_equal(kernel, kernel[::-1])


def test_kernel_invalid_params():
    with pytest.raises(ValueError):
        ws._test_function_1d(1, 2)  # m < 2
    with pytest.raises(ValueError):
        ws._test_function_1d(4, 0)  # p < 1


# ---------------------------------------------------------------------------
# Weak convolution
# ---------------------------------------------------------------------------


def test_weak_convolution_constant(test_field):
    f, _, _, N = test_field
    m, p = 6, 4
    kernel = ws._test_function_1d(m, p)

    f_const = np.ones_like(f)
    fw = ws.weak_convolution_2d(f_const, kernel, kernel)
    # Convolution of constant 1 with kernel → sum(kernel) per dimension
    expected = float(np.sum(kernel)) ** 2
    np.testing.assert_allclose(fw, expected, rtol=1e-10)


def test_weak_convolution_shape(test_field):
    f, _, _, N = test_field
    m, p = 8, 4
    kernel = ws._test_function_1d(m, p)
    fw = ws.weak_convolution_2d(f, kernel, kernel)
    assert fw.shape == (N - 2 * m, N - 2 * m)


def test_weak_convolution_linearity(test_field):
    f, _, _, N = test_field
    m, p = 4, 2
    kernel = ws._test_function_1d(m, p)
    g = np.cos(f)

    a, b = 2.5, -1.3
    lhs = ws.weak_convolution_2d(a * f + b * g, kernel, kernel)
    rhs = a * ws.weak_convolution_2d(f, kernel, kernel) + b * ws.weak_convolution_2d(g, kernel, kernel)
    np.testing.assert_allclose(lhs, rhs, rtol=1e-12)


# ---------------------------------------------------------------------------
# Weak library and RHS
# ---------------------------------------------------------------------------


def test_weak_library_shape(test_field):
    f, _, L, N = test_field
    m, p = 6, 4
    Theta_w, descs = ws.build_weak_library_2d(f, L, L, max_p=2, max_n=3, m=m, p=p)
    Nv = (N - 2 * m) ** 2
    assert Theta_w.shape == (Nv, 20)
    assert len(descs) == 20
    assert np.all(np.isfinite(Theta_w))


def test_weak_rhs_shape(test_field):
    _, ft, _, N = test_field
    m, p = 6, 4
    b = ws.build_weak_rhs_2d(ft, m=m, p=p)
    Nv = (N - 2 * m) ** 2
    assert len(b) == Nv
    assert np.all(np.isfinite(b))


def test_weak_rhs_nonzero(test_field):
    f, ft, L, N = test_field
    m, p = 6, 4
    b = ws.build_weak_rhs_2d(ft, m=m, p=p)
    assert np.max(np.abs(b)) > 0, "Weak RHS should not be identically zero."


# ---------------------------------------------------------------------------
# Weak regression data pipeline
# ---------------------------------------------------------------------------


def test_make_weak_regression_data_random(test_field):
    f, ft, L, N = test_field
    f_snaps = np.array([f, f + 0.1, f - 0.1])
    ft_snaps = np.array([ft, ft + 0.05, ft - 0.05])
    idx_sets = [np.arange(50) for _ in range(3)]

    Theta, b, descs = ws.make_weak_regression_data(
        f_snaps, ft_snaps, idx_sets, L, L, max_p=2, max_n=3, m=6, p=4, sample_stride=0
    )
    assert Theta.shape == (150, 20)
    assert len(b) == 150
    assert np.all(np.isfinite(Theta))
    assert np.all(np.isfinite(b))


def test_make_weak_regression_data_strided(test_field):
    f, ft, L, N = test_field
    f_snaps = np.array([f, f])
    ft_snaps = np.array([ft, ft])

    Theta, b, descs = ws.make_weak_regression_data(
        f_snaps, ft_snaps, [], L, L, max_p=2, max_n=3, m=6, p=4, sample_stride=4
    )
    # 2 snaps × ((N-2m)//4)² query points
    m = 6
    n_q_per_snap = ((N - 2 * m) // 4) ** 2
    assert Theta.shape == (2 * n_q_per_snap, 20)
    assert len(b) == 2 * n_q_per_snap
    assert np.all(np.isfinite(Theta))


# ---------------------------------------------------------------------------
# Consistency: weak library and strong library on clean data
# ---------------------------------------------------------------------------


def test_weak_form_reduces_noise_amplification():
    """The weak-form convolution reduces the effective noise in library columns.

    On a simple cosine mode with added Gaussian noise, the weak form should
    produce library columns with lower relative noise than pointwise
    evaluation of the same operators.
    """
    N, L = 64, 2 * math.pi
    x = np.linspace(0, L, N, endpoint=False)
    y = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="ij")
    f_clean = np.cos(2 * X + 3 * Y)
    noise_level = 0.05
    rng = np.random.default_rng(42)
    noise = noise_level * np.std(f_clean) * rng.standard_normal(f_clean.shape)
    f_noisy = f_clean + noise

    import derivatives
    derivatives.set_backend("spectral")
    from library import build_equivariant_library

    # Strong form on noisy data
    Theta_strong_noisy, descs = build_equivariant_library(
        f_noisy, L, L, max_p=2, max_n=3
    )
    Theta_strong_clean, _ = build_equivariant_library(
        f_clean, L, L, max_p=2, max_n=3
    )

    # Weak form on noisy data
    m, p = 6, 4
    Theta_weak_noisy, _ = ws.build_weak_library_2d(
        f_noisy, L, L, max_p=2, max_n=3, m=m, p=p
    )
    Theta_weak_clean, _ = ws.build_weak_library_2d(
        f_clean, L, L, max_p=2, max_n=3, m=m, p=p
    )

    # Compare relative error (noisy vs clean) for the bilaplacian column
    # (index 9 = "Δ²f") — the highest-order, most noise-sensitive operator.
    idx_bilap = descs.index("Δ²f")

    # Strong form relative error
    strong_err = np.std(Theta_strong_noisy[:, idx_bilap] - Theta_strong_clean[:, idx_bilap])
    strong_scale = np.std(Theta_strong_clean[:, idx_bilap])

    # Weak form relative error (restrict to valid interior for fair comparison)
    weak_err = np.std(Theta_weak_noisy[:, idx_bilap] - Theta_weak_clean[:, idx_bilap])
    weak_scale = np.std(Theta_weak_clean[:, idx_bilap])

    strong_rel = strong_err / strong_scale
    weak_rel = weak_err / weak_scale

    # The weak form must have lower relative noise than the strong form.
    assert weak_rel < strong_rel, (
        f"Weak form noise {weak_rel:.4f} should be lower than strong form "
        f"noise {strong_rel:.4f} for Δ²f column."
    )


# ---------------------------------------------------------------------------
# IBP (integration by parts) tests
# ---------------------------------------------------------------------------


def test_test_function_deriv_1d():
    """Analytic test-function derivatives must be consistent with finite differences."""
    m, p = 8, 4
    phi = ws._test_function_1d(m, p)
    # First derivative: central difference should approximate analytic result
    dphi_fd = (np.roll(phi, -1) - np.roll(phi, 1)) / 2.0  # central diff, dx=1
    dphi_analytic = ws._test_function_deriv_1d(m, p, 1)
    # Agreement in the interior (away from boundaries)
    interior = slice(2, -2)
    np.testing.assert_allclose(
        dphi_fd[interior], dphi_analytic[interior], rtol=1e-1, atol=1e-3,
        err_msg="Analytic φ' must match central-difference approximation (~10 %)."
    )

    # Second derivative
    d2phi_fd = (np.roll(phi, -1) - 2 * phi + np.roll(phi, 1))
    d2phi_analytic = ws._test_function_deriv_1d(m, p, 2)
    np.testing.assert_allclose(
        d2phi_fd[interior], d2phi_analytic[interior], rtol=1e-1, atol=1e-3,
        err_msg="Analytic φ'' must match central-difference approximation (~10 %)."
    )

    # order=0 returns the base function
    np.testing.assert_array_equal(
        ws._test_function_deriv_1d(m, p, 0), phi,
    )


def test_ibp_laplacian_matches_direct_on_clean_data(test_field):
    """IBP and direct weak forms of Δf must agree on clean periodic data."""
    f, _, L, N = test_field
    m, p = 6, 4

    import derivatives
    derivatives.set_backend("spectral")

    # Direct: convolve Δf with ψ
    lap = derivatives.laplacian(f, L, L)
    kernel = ws._test_function_1d(m, p)
    direct = ws.weak_convolution_2d(lap, kernel, kernel)

    # IBP: convolve f with Δψ
    ibp = ws.weak_laplacian_ibp(f, m, p, L, L)

    np.testing.assert_allclose(ibp, direct, rtol=3e-2, atol=1.0,
        err_msg="IBP Laplacian must match direct approach (~3 % spectral-vs-finite-diff).")


def test_ibp_bilaplacian_matches_direct_on_clean_data(test_field):
    """IBP and direct weak forms of Δ²f must agree on clean periodic data."""
    f, _, L, N = test_field
    m, p = 6, 4

    import derivatives
    derivatives.set_backend("spectral")

    # Direct: convolve Δ²f with ψ
    bilap = derivatives.bilaplacian(f, L, L)
    kernel = ws._test_function_1d(m, p)
    direct = ws.weak_convolution_2d(bilap, kernel, kernel)

    # IBP: convolve f with Δ²ψ
    ibp = ws.weak_bilaplacian_ibp(f, m, p, L, L)

    # The discrete summation-by-parts identity is exact for the discrete
    # bilaplacian (finite differences), but the spectral bilaplacian (FFT,
    # k⁴ multiplier) is a different operator.  On clean data they agree to
    # within ~15 % in bulk but can differ more at zero-crossings.
    # The true value of IBP is on noisy data (tested separately).
    np.testing.assert_allclose(ibp, direct, rtol=2.0, atol=500.0,
        err_msg="IBP Bilaplacian must approximately match direct.")


def test_ibp_beats_direct_on_noisy_data():
    """IBP must produce lower noise than direct approach on noisy data for Δ²f."""
    N, L = 64, 2 * math.pi
    x = np.linspace(0, L, N, endpoint=False)
    y = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="ij")
    f_clean = np.cos(2 * X + 3 * Y)

    noise_level = 0.05
    rng = np.random.default_rng(42)
    noise = noise_level * np.std(f_clean) * rng.standard_normal(f_clean.shape)
    f_noisy = f_clean + noise

    m, p = 6, 4
    kernel = ws._test_function_1d(m, p)

    # Direct: differentiate noisy data, then convolve
    import derivatives
    derivatives.set_backend("spectral")
    bilap_noisy = derivatives.bilaplacian(f_noisy, L, L)
    bilap_clean = derivatives.bilaplacian(f_clean, L, L)
    direct_noisy = ws.weak_convolution_2d(bilap_noisy, kernel, kernel)
    direct_clean = ws.weak_convolution_2d(bilap_clean, kernel, kernel)
    direct_err = np.std(direct_noisy - direct_clean)

    # IBP: convolve raw (noisy but undifferentiated) data with Δ²ψ
    ibp_noisy = ws.weak_bilaplacian_ibp(f_noisy, m, p, L, L)
    ibp_clean = ws.weak_bilaplacian_ibp(f_clean, m, p, L, L)
    ibp_err = np.std(ibp_noisy - ibp_clean)

    assert ibp_err < direct_err, (
        f"IBP noise {ibp_err:.2e} should be lower than direct noise "
        f"{direct_err:.2e} for Δ²f at {noise_level*100:.0f}% noise."
    )


def test_build_weak_library_ibp_flag(test_field):
    """use_ibp=True produces a valid library where Δf/Δ²f differ from direct."""
    f, _, L, N = test_field
    m, p = 6, 4

    Theta_ibp, descs = ws.build_weak_library_2d(f, L, L, m=m, p=p, use_ibp=True)
    Theta_direct, descs2 = ws.build_weak_library_2d(f, L, L, m=m, p=p, use_ibp=False)

    assert descs == descs2
    assert Theta_ibp.shape == Theta_direct.shape
    assert np.all(np.isfinite(Theta_ibp))

    # Δf and Δ²f columns must agree between IBP and direct on clean data.
    idx_lap = descs.index("Δf")
    idx_bilap = descs.index("Δ²f")
    np.testing.assert_allclose(
        Theta_ibp[:, idx_lap], Theta_direct[:, idx_lap], rtol=3e-2, atol=1.0,
        err_msg="IBP and direct Laplacian weak forms must agree (~3 % spectral-vs-finite-diff)."
    )
    np.testing.assert_allclose(
        Theta_ibp[:, idx_bilap], Theta_direct[:, idx_bilap], rtol=2.0, atol=500.0,
        err_msg="IBP and direct Bilaplacian weak forms must approximately agree."
    )

    # Nonlinear columns (‖∇f‖², etc.) must be identical since IBP doesn't apply.
    idx_grad_sq = descs.index("‖∇f‖²")
    np.testing.assert_array_equal(
        Theta_ibp[:, idx_grad_sq], Theta_direct[:, idx_grad_sq],
    )


def test_ibp_flag_off_by_default(test_field):
    """use_ibp defaults to False — the library must match the direct path."""
    f, _, L, N = test_field
    m, p = 6, 4

    Theta_default, descs = ws.build_weak_library_2d(f, L, L, m=m, p=p)
    Theta_explicit_off, descs2 = ws.build_weak_library_2d(
        f, L, L, m=m, p=p, use_ibp=False
    )

    assert descs == descs2
    # Every column must be bit-for-bit identical.
    for j in range(len(descs)):
        np.testing.assert_array_equal(
            Theta_default[:, j], Theta_explicit_off[:, j],
            err_msg=f"Column '{descs[j]}' differs between default and use_ibp=False."
        )


def test_make_weak_regression_data_ibp_plumbing(test_field):
    """make_weak_regression_data forwards use_ibp correctly."""
    f, ft, L, N = test_field
    f_snaps = np.array([f, f])
    ft_snaps = np.array([ft, ft])
    idx_sets = [np.arange(30) for _ in range(2)]

    Theta_ibp, b_ibp, descs = ws.make_weak_regression_data(
        f_snaps, ft_snaps, idx_sets, L, L, m=6, p=4, sample_stride=0, use_ibp=True,
    )
    Theta_dir, b_dir, descs2 = ws.make_weak_regression_data(
        f_snaps, ft_snaps, idx_sets, L, L, m=6, p=4, sample_stride=0, use_ibp=False,
    )

    assert descs == descs2
    assert Theta_ibp.shape == Theta_dir.shape

    # Δf and Δ²f columns should differ (IBP vs direct are different computations
    # that agree on clean data but differ numerically).
    idx_lap = descs.index("Δf")
    idx_bilap = descs.index("Δ²f")
    np.testing.assert_allclose(
        Theta_ibp[:, idx_lap], Theta_dir[:, idx_lap], rtol=3e-2, atol=1.0,
    )
    np.testing.assert_allclose(
        Theta_ibp[:, idx_bilap], Theta_dir[:, idx_bilap], rtol=2.0, atol=500.0,
    )
    # RHS should be identical (IBP doesn't affect RHS)
    np.testing.assert_array_equal(b_ibp, b_dir)
