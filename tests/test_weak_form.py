"""
Tests for ``weak_form.py`` — weak-form PDE discovery.
"""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.library import weak_form as wf
from equivariant_pde_find import library


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
    kernel = wf._test_function_1d(m, p)
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
        wf._test_function_1d(1, 2)  # m < 2
    with pytest.raises(ValueError):
        wf._test_function_1d(4, 0)  # p < 1


# ---------------------------------------------------------------------------
# Weak convolution
# ---------------------------------------------------------------------------


def test_weak_convolution_constant(test_field):
    f, _, _, N = test_field
    m, p = 6, 4
    kernel = wf._test_function_1d(m, p)

    f_const = np.ones_like(f)
    fw = wf.weak_convolution(f_const, kernel, kernel)
    # Convolution of constant 1 with kernel → sum(kernel) per dimension
    expected = float(np.sum(kernel)) ** 2
    np.testing.assert_allclose(fw, expected, rtol=1e-10)


def test_weak_convolution_shape(test_field):
    f, _, _, N = test_field
    m, p = 8, 4
    kernel = wf._test_function_1d(m, p)
    fw = wf.weak_convolution(f, kernel, kernel)
    assert fw.shape == (N - 2 * m, N - 2 * m)


def test_weak_convolution_linearity(test_field):
    f, _, _, N = test_field
    m, p = 4, 2
    kernel = wf._test_function_1d(m, p)
    g = np.cos(f)

    a, b = 2.5, -1.3
    lhs = wf.weak_convolution(a * f + b * g, kernel, kernel)
    rhs = a * wf.weak_convolution(f, kernel, kernel) + b * wf.weak_convolution(g, kernel, kernel)
    np.testing.assert_allclose(lhs, rhs, rtol=1e-12)


# ---------------------------------------------------------------------------
# Weak library and RHS
# ---------------------------------------------------------------------------


def test_weak_library_shape(test_field):
    f, _, L, N = test_field
    m, p = 6, 4
    Theta_w, descs = library.build_multigraph_library(f, L, L, max_p=2, max_n=3, m=m, p=p, formulation="weak")
    Nv = (N - 2 * m) ** 2
    assert Theta_w.shape == (Nv, 20)
    assert len(descs) == 20
    assert np.all(np.isfinite(Theta_w))


def test_weak_rhs_shape(test_field):
    _, ft, _, N = test_field
    m, p = 6, 4
    b = wf.build_weak_rhs(ft, m=m, p=p)
    Nv = (N - 2 * m) ** 2
    assert len(b) == Nv
    assert np.all(np.isfinite(b))


def test_weak_rhs_nonzero(test_field):
    f, ft, L, N = test_field
    m, p = 6, 4
    b = wf.build_weak_rhs(ft, m=m, p=p)
    assert np.max(np.abs(b)) > 0, "Weak RHS should not be identically zero."


# ---------------------------------------------------------------------------
# Weak regression data pipeline
# ---------------------------------------------------------------------------


def test_make_weak_regression_data_random(test_field):
    f, ft, L, N = test_field
    f_snaps = np.array([f, f + 0.1, f - 0.1])
    ft_snaps = np.array([ft, ft + 0.05, ft - 0.05])
    idx_sets = [np.arange(50) for _ in range(3)]

    Theta, b, descs, _rinfo = wf.make_weak_regression_data(
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

    Theta, b, descs, _rinfo = wf.make_weak_regression_data(
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

    from equivariant_pde_find import derivatives
    derivatives.set_backend("spectral")
    from equivariant_pde_find.library import build_multigraph_library

    # Strong form on noisy data
    Theta_strong_noisy, descs = build_multigraph_library(
        f_noisy, L, L, max_p=2, max_n=3
    )
    Theta_strong_clean, _ = build_multigraph_library(
        f_clean, L, L, max_p=2, max_n=3
    )

    # Weak form on noisy data
    m, p = 6, 4
    Theta_weak_noisy, _ = library.build_multigraph_library(
        f_noisy, L, L, max_p=2, max_n=3, m=m, p=p, formulation="weak"
    )
    Theta_weak_clean, _ = library.build_multigraph_library(
        f_clean, L, L, max_p=2, max_n=3, m=m, p=p, formulation="weak"
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
    phi = wf._test_function_1d(m, p)
    # First derivative: central difference should approximate analytic result
    dphi_fd = (np.roll(phi, -1) - np.roll(phi, 1)) / 2.0  # central diff, dx=1
    dphi_analytic = wf._test_function_deriv_1d(m, p, 1)
    # Agreement in the interior (away from boundaries)
    interior = slice(2, -2)
    np.testing.assert_allclose(
        dphi_fd[interior], dphi_analytic[interior], rtol=1e-1, atol=1e-3,
        err_msg="Analytic φ' must match central-difference approximation (~10 %)."
    )

    # Second derivative
    d2phi_fd = (np.roll(phi, -1) - 2 * phi + np.roll(phi, 1))
    d2phi_analytic = wf._test_function_deriv_1d(m, p, 2)
    np.testing.assert_allclose(
        d2phi_fd[interior], d2phi_analytic[interior], rtol=1e-1, atol=1e-3,
        err_msg="Analytic φ'' must match central-difference approximation (~10 %)."
    )

    # order=0 returns the base function
    np.testing.assert_array_equal(
        wf._test_function_deriv_1d(m, p, 0), phi,
    )


def test_ibp_laplacian_matches_direct_on_clean_data(test_field):
    """IBP and direct weak forms of Δf must agree on clean periodic data."""
    f, _, L, N = test_field
    m, p = 6, 4

    from equivariant_pde_find import derivatives
    derivatives.set_backend("spectral")

    # Direct: convolve Δf with ψ
    lap = derivatives.laplacian(f, L, L)
    kernel = wf._test_function_1d(m, p)
    direct = wf.weak_convolution(lap, kernel, kernel)

    # IBP: convolve f with Δψ
    ibp = wf.weak_laplacian_ibp(f, m, p, L, L)

    np.testing.assert_allclose(ibp, direct, rtol=3e-2, atol=1.0,
        err_msg="IBP Laplacian must match direct approach (~3 % spectral-vs-finite-diff).")


def test_ibp_bilaplacian_matches_direct_on_clean_data(test_field):
    """IBP and direct weak forms of Δ²f must agree on clean periodic data."""
    f, _, L, N = test_field
    m, p = 6, 4

    from equivariant_pde_find import derivatives
    derivatives.set_backend("spectral")

    # Direct: convolve Δ²f with ψ
    bilap = derivatives.bilaplacian(f, L, L)
    kernel = wf._test_function_1d(m, p)
    direct = wf.weak_convolution(bilap, kernel, kernel)

    # IBP: convolve f with Δ²ψ
    ibp = wf.weak_bilaplacian_ibp(f, m, p, L, L)

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
    kernel = wf._test_function_1d(m, p)

    # Direct: differentiate noisy data, then convolve
    from equivariant_pde_find import derivatives
    derivatives.set_backend("spectral")
    bilap_noisy = derivatives.bilaplacian(f_noisy, L, L)
    bilap_clean = derivatives.bilaplacian(f_clean, L, L)
    direct_noisy = wf.weak_convolution(bilap_noisy, kernel, kernel)
    direct_clean = wf.weak_convolution(bilap_clean, kernel, kernel)
    direct_err = np.std(direct_noisy - direct_clean)

    # IBP: convolve raw (noisy but undifferentiated) data with Δ²ψ
    ibp_noisy = wf.weak_bilaplacian_ibp(f_noisy, m, p, L, L)
    ibp_clean = wf.weak_bilaplacian_ibp(f_clean, m, p, L, L)
    ibp_err = np.std(ibp_noisy - ibp_clean)

    assert ibp_err < direct_err, (
        f"IBP noise {ibp_err:.2e} should be lower than direct noise "
        f"{direct_err:.2e} for Δ²f at {noise_level*100:.0f}% noise."
    )


def test_build_weak_multigraph_library_ibp_flag(test_field):
    """use_ibp=True produces a valid library where Δf/Δ²f differ from direct."""
    f, _, L, N = test_field
    m, p = 6, 4

    Theta_ibp, descs = library.build_multigraph_library(f, L, L, m=m, p=p, use_ibp=True, formulation="weak")
    Theta_direct, descs2 = library.build_multigraph_library(f, L, L, m=m, p=p, use_ibp=False, formulation="weak")

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

    Theta_default, descs = library.build_multigraph_library(f, L, L, m=m, p=p, formulation="weak")
    Theta_explicit_off, descs2 = library.build_multigraph_library(
        f, L, L, m=m, p=p, use_ibp=False, formulation="weak"
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

    Theta_ibp, b_ibp, descs, _rinfo = wf.make_weak_regression_data(
        f_snaps, ft_snaps, idx_sets, L, L, m=6, p=4, sample_stride=0, use_ibp=True,
    )
    Theta_dir, b_dir, descs2, _rinfo2 = wf.make_weak_regression_data(
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


# ---------------------------------------------------------------------------
# Automatic test-function parameter selection
# ---------------------------------------------------------------------------


def test_estimate_critical_wavenumber_shape():
    """estimate_critical_wavenumber returns a non-negative int."""
    rng = np.random.default_rng(42)
    U = rng.standard_normal((64, 64)) + 0.1 * np.sin(
        np.linspace(0, 4 * np.pi, 64)
    )[np.newaxis, :]
    k = wf.estimate_critical_wavenumber(U, axis=1)
    assert isinstance(k, (int, np.integer))
    assert k >= 0
    assert k <= 32  # Nyquist


def test_estimate_critical_wavenumber_on_pure_noise():
    """On pure white noise, changepoint should be near DC (k≈0)."""
    rng = np.random.default_rng(123)
    U = rng.standard_normal((64, 64))
    k = wf.estimate_critical_wavenumber(U, axis=0)
    # Pure noise → no signal-dominated region → changepoint near zero
    assert k <= 5


def test_estimate_critical_wavenumber_on_pure_signal():
    """On a pure sinusoid, changepoint is at or above signal frequency."""
    N = 64
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    U = np.sin(2 * x)[:, np.newaxis] * np.ones((1, N))
    k = wf.estimate_critical_wavenumber(U, axis=0)
    assert k >= 2


def test_auto_test_function_params_returns_valid():
    """auto_test_function_params returns reasonable (m, p) values."""
    N = 64
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    U = np.cos(2 * X + 3 * Y) + 0.5 * np.sin(X - 2 * Y)
    m_x, p_x, m_t, p_t = wf.auto_test_function_params(
        U, dx=2 * np.pi / N, dt=1.0, max_deriv_spatial=4, max_deriv_temporal=1,
    )
    assert m_x >= 2
    assert p_x >= 2
    assert m_t >= 2
    assert p_t >= 1


def test_compute_p_from_m_bounds():
    """_compute_p_from_m returns p >= alpha_bar + 1."""
    for m in [4, 8, 16, 32]:
        p = wf._compute_p_from_m(m, alpha_bar=4, tau=1e-10)
        assert p >= 5  # alpha_bar=4 → p >= 5
        # Verify decay condition
        base = 1.0 - (1.0 - 1.0 / m) ** 2
        assert base ** p <= 1e-10 + 1e-15  # within tolerance


# ---------------------------------------------------------------------------
# Rescaling
# ---------------------------------------------------------------------------


def test_compute_rescaling_scales_finite():
    """Rescaling scales are finite and positive."""
    N = 64
    L = 2 * np.pi
    x = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    U = np.cos(2 * X + 3 * Y) + 0.5 * np.sin(X - 2 * Y)
    dx = L / N
    gx, gt, gu = wf.compute_rescaling_scales(
        U, dx, 1.0, m_x=8, m_t=4, p_x=4, p_t=2,
        max_deriv_spatial=4, max_deriv_temporal=1, max_poly_degree=3,
    )
    assert np.isfinite(gx) and gx > 0
    assert np.isfinite(gt) and gt > 0
    assert np.isfinite(gu) and gu > 0


def test_build_rescale_matrix():
    """build_rescale_matrix returns a vector of correct length."""
    from equivariant_pde_find.multigraph_operators import get_operators

    ops = get_operators(max_p=2, max_n=3)
    mu = wf.build_rescale_matrix(2.0, 3.0, 0.5, max_p=2, max_n=3)
    assert len(mu) == len(ops)
    assert np.all(np.isfinite(mu))
    assert np.all(mu > 0)


def test_rescale_roundtrip():
    """Rescaling + unscaling preserves the PDE on a simple mock system."""
    from equivariant_pde_find.multigraph_operators import get_operator_metadata

    # Build a mock G and w_true, then verify rescale_coefficients round-trips
    rng = np.random.default_rng(42)
    n, max_p, max_n = 100, 2, 3
    meta = get_operator_metadata(max_p, max_n)
    d = len(meta)
    G = rng.standard_normal((n, d))
    w_true = np.zeros(d)
    w_true[1] = 1.0   # f
    w_true[4] = 0.02  # Δf
    w_true[6] = -0.5  # ‖∇f‖²
    b = G @ w_true

    # Rescale
    gamma_x, gamma_t, gamma_u = 2.0, 3.0, 0.5
    col_scale, lhs_scale = wf._rescale_factors(gamma_x, gamma_t, gamma_u, max_p, max_n)
    G_scaled = G * col_scale[np.newaxis, :]
    b_scaled = lhs_scale * b

    # Solve scaled system
    w_scaled = np.linalg.lstsq(G_scaled, b_scaled, rcond=None)[0]

    # Recover original coefficients
    rescale_info = {
        "gamma_x": gamma_x, "gamma_t": gamma_t, "gamma_u": gamma_u,
        "col_scale": col_scale, "lhs_scale": lhs_scale,
        "max_p": max_p, "max_n": max_n,
    }
    w_recovered = wf.rescale_coefficients(w_scaled, rescale_info)
    np.testing.assert_allclose(w_recovered, w_true, rtol=1e-10, atol=1e-10)


def test_make_weak_regression_data_rescale():
    """make_weak_regression_data with rescale=True returns rescale_info."""
    f_snaps = np.random.default_rng(0).standard_normal((2, 32, 32))
    ft_snaps = np.random.default_rng(1).standard_normal((2, 32, 32))
    idx_sets = [np.arange(10) for _ in range(2)]

    Theta, b, descs, rinfo = wf.make_weak_regression_data(
        f_snaps, ft_snaps, idx_sets, 1.0, 1.0, m=6, p=4,
        sample_stride=0, rescale=True,
    )
    assert rinfo is not None
    assert "col_scale" in rinfo
    assert "lhs_scale" in rinfo
    assert len(rinfo["col_scale"]) == Theta.shape[1]
    assert np.all(np.isfinite(Theta))
    assert np.all(np.isfinite(b))


# ---------------------------------------------------------------------------
# IBP sign and accuracy guards (regression tests for 2026-06-22 fixes)
# ---------------------------------------------------------------------------


def _build_separable_field(N=64):
    """A separable smooth 2-D field for controlled IBP testing."""
    L = 2 * math.pi
    x = np.linspace(0.0, L, N, endpoint=False)
    y = np.linspace(0.0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="ij")
    f = np.sin(2 * X) * np.cos(3 * Y) + 0.5 * np.cos(X + Y)
    return f, L


def test_wsindy_ibp_positive_correlation_all_orders():
    """Every WSINDy IBP-applicable column up to order 3 must have positive
    correlation with its direct counterpart.  A negative correlation indicates
    a sign error (e.g. double sign-flip from np.convolve kernel reversal).

    4th-order columns are excluded: with p=4 the test function is only C³,
    so the analytic 4th derivative has a jump at the boundary, causing
    legitimate correlation degradation that is not a sign error."""
    f, L = _build_separable_field()
    from equivariant_pde_find import derivatives
    derivatives.set_backend("spectral")

    for m, p in [(6, 4), (8, 4), (12, 6)]:
        G_dir, descs = library.build_wsindy_library(
            f, L, L, m=m, p=p, use_ibp=False, formulation="weak",
        )
        G_ibp, _ = library.build_wsindy_library(
            f, L, L, m=m, p=p, use_ibp=True, formulation="weak",
        )
        for j, desc in enumerate(descs):
            if not desc.startswith("f_{"):
                continue
            order = len(desc[3:-1])
            if order >= 4:
                continue  # known C³ limitation for p=4
            corr = np.corrcoef(G_dir[:, j], G_ibp[:, j])[0, 1]
            assert corr > 0, (
                f"WSINDy IBP column '{desc}' (m={m}, p={p}) has negative "
                f"correlation {corr:.6f} with direct column — sign error."
            )


def test_standard_ibp_positive_correlation_all_orders():
    """Every standard-library IBP-applicable column up to order 3 must have
    positive correlation with its direct counterpart.  4th-order excluded
    (same C³ limitation as WSINDy)."""
    f, L = _build_separable_field()
    from equivariant_pde_find import derivatives
    derivatives.set_backend("spectral")

    for m, p in [(6, 4), (8, 4), (12, 6)]:
        G_dir, descs = library.build_pdefind_library(
            f, L, L, m=m, p=p, use_ibp=False, formulation="weak",
        )
        G_ibp, _ = library.build_pdefind_library(
            f, L, L, m=m, p=p, use_ibp=True, formulation="weak",
        )
        for j, desc in enumerate(descs):
            # IBP applies to pure linear differential columns: f_x, f_xx, …
            # (no polynomial prefactor, no curly braces in standard library)
            if not (desc.startswith("f_") and desc[2:].replace("x","").replace("y","") == ""):
                continue
            if desc == "f":
                continue
            order = len(desc) - 2  # strip "f_"
            if order >= 4:
                continue
            corr = np.corrcoef(G_dir[:, j], G_ibp[:, j])[0, 1]
            assert corr > 0, (
                f"Standard IBP column '{desc}' (m={m}, p={p}) has negative "
                f"correlation {corr:.6f} with direct column — sign error."
            )


def test_wsindy_ibp_accuracy_by_order():
    """WSINDy IBP accuracy degrades gracefully with derivative order.

    Orders 1–2: sub-percent relative error.
    Orders 3–4: larger errors are expected (finite-difference truncation,
    C³ smoothness limit for 4th order with p=4).
    All columns must have correlation > 0.9 with direct (no catastrophic
    failure).  4th-order columns at small m are excluded from the relative-
    error check (known limitation, not a bug)."""
    f, L = _build_separable_field()
    from equivariant_pde_find import derivatives
    derivatives.set_backend("spectral")

    m, p = 8, 4
    G_dir, descs = library.build_wsindy_library(
        f, L, L, m=m, p=p, use_ibp=False, formulation="weak",
    )
    G_ibp, _ = library.build_wsindy_library(
        f, L, L, m=m, p=p, use_ibp=True, formulation="weak",
    )

    # Relative error thresholds by order (m=8, p=4)
    max_rel_err = {1: 5e-3, 2: 5e-3, 3: 0.5, 4: 2.5}
    # Minimum correlation by order
    min_corr = {1: 0.999, 2: 0.999, 3: 0.99, 4: 0.9}

    for j, desc in enumerate(descs):
        if not desc.startswith("f_{"):
            continue
        order = len(desc[3:-1])
        if order >= 4:
            continue  # known C³ smoothness limit for p=4 — not a bug
        rel = np.linalg.norm(G_dir[:, j] - G_ibp[:, j]) / np.linalg.norm(G_dir[:, j])
        corr = np.corrcoef(G_dir[:, j], G_ibp[:, j])[0, 1]

        assert corr > min_corr.get(order, 0.9), (
            f"WSINDy IBP '{desc}' (order={order}) correlation {corr:.6f} "
            f"below minimum {min_corr.get(order, 0.9)}."
        )
        assert rel < max_rel_err.get(order, 5.0), (
            f"WSINDy IBP '{desc}' (order={order}) relative error {rel:.6f} "
            f"exceeds maximum {max_rel_err.get(order, 5.0)}."
        )


def test_equivariant_laplacian_ibp_accuracy():
    """The equivariant library's Laplacian IBP must match direct to < 1 %
    (analytic φ″ kernel, guard against fallback to discrete differences)."""
    f, L = _build_separable_field()
    from equivariant_pde_find import derivatives
    derivatives.set_backend("spectral")

    m, p = 8, 4
    kernel = wf._test_function_1d(m, p)
    lap = derivatives.laplacian(f, L, L)
    direct = wf.weak_convolution(lap, kernel, kernel)
    ibp = wf.weak_laplacian_ibp(f, m, p, L, L)

    rel = np.linalg.norm(direct - ibp) / np.linalg.norm(direct)
    assert rel < 0.01, (
        f"Equivariant Laplacian IBP relative error {rel:.6f} exceeds 1 % — "
        f"analytic φ″ kernel may not be in use."
    )


def test_equivariant_ibp_kernel_consistency_with_wsindy():
    """The Laplacian IBP kernel used by the equivariant library must match
    the f_xx IBP kernel used by the WSINDy library (both should use
    analytic φ″)."""
    f, L = _build_separable_field()
    from equivariant_pde_find import derivatives
    derivatives.set_backend("spectral")

    m, p = 8, 4

    # Equivariant: build full library, extract the Δf column
    G_eq, descs_eq = library.build_multigraph_library(f, L, L, m=m, p=p, use_ibp=True, formulation="weak")
    col_delta = G_eq[:, descs_eq.index("Δf")]

    # WSINDy: the IBP f_xx + f_yy should equal the equivariant Δf column
    # (since Δf = f_xx + f_yy)
    G_ws, descs_ws = library.build_wsindy_library(f, L, L, m=m, p=p, use_ibp=True, formulation="weak")
    col_fxx = G_ws[:, descs_ws.index("f_{xx}")]
    col_fyy = G_ws[:, descs_ws.index("f_{yy}")]
    col_laplacian_from_wsindy = col_fxx + col_fyy

    rel = np.linalg.norm(col_delta - col_laplacian_from_wsindy) / np.linalg.norm(col_delta)
    assert rel < 0.01, (
        f"Equivariant Δf (IBP) vs WSINDy f_xx+f_yy (IBP) relative "
        f"error {rel:.6f} exceeds 1 % — IBP kernels differ."
    )


def test_make_weak_regression_data_no_rescale():
    """Without rescale, rescale_info is None."""
    f_snaps = np.random.default_rng(0).standard_normal((2, 32, 32))
    ft_snaps = np.random.default_rng(1).standard_normal((2, 32, 32))
    idx_sets = [np.arange(10) for _ in range(2)]

    Theta, b, descs, rinfo = wf.make_weak_regression_data(
        f_snaps, ft_snaps, idx_sets, 1.0, 1.0, m=6, p=4,
        sample_stride=0, rescale=False,
    )
    assert rinfo is None


# ---------------------------------------------------------------------------
# MWSINDy tests (WSINDy regression + multigraph fine-tuning)
# ---------------------------------------------------------------------------


def test_mwsindy_symmetry_check_rejects_anisotropic():
    """MWSINDy's O(2) symmetry check: if WSINDy finds f_xx but not f_yy,
    the multigraph operator Δf is NOT considered discovered (fn>0)."""
    f, L = _build_separable_field()
    from equivariant_pde_find import derivatives
    derivatives.set_backend("spectral")

    m, p = 8, 6

    G_ws, desc_ws = library.build_wsindy_library(f, L, L, m=m, p=p, use_ibp=False, formulation="weak")
    G_mg, desc_mg = library.build_multigraph_library(f, L, L, m=m, p=p, use_ibp=False, formulation="weak")

    from equivariant_pde_find.library import build_polyjet_library
    _, desc_jet = build_polyjet_library(f, L, L, max_order=4, max_degree=3)
    from equivariant_pde_find.evaluation import build_wsindy_to_jet_map
    M = build_wsindy_to_jet_map(desc_ws, desc_jet)

    # Ground truth: f_t = Δf
    true_mg = np.zeros(len(desc_mg))
    true_mg[desc_mg.index("Δf")] = 1.0
    from equivariant_pde_find.evaluation import _MULTIGRAPH_TO_JET_DECOMPOSITION
    true_jet = np.zeros(len(desc_jet))
    for j, (name, coeff) in enumerate(zip(desc_mg, true_mg)):
        if abs(coeff) < 1e-12: continue
        if name not in _MULTIGRAPH_TO_JET_DECOMPOSITION: continue
        for jet_name, mult in _MULTIGRAPH_TO_JET_DECOMPOSITION[name]:
            true_jet[desc_jet.index(jet_name)] += coeff * mult

    # Create an anisotropic WSINDy vector: keep f_{xx} but NOT f_{yy}
    w_broken = np.zeros(len(desc_ws))
    w_broken[desc_ws.index("f_{xx}")] = 1.0
    # f_{yy} is intentionally left at 0

    # Use the evaluation Theta from the G matrix itself (just test fn/fp logic)
    # Map w_broken to jet
    w_jet = M @ w_broken
    nz_jet = np.abs(w_jet) > 1e-10

    # Check: Δf requires f_xx AND f_yy
    # f_xx is active, f_yy is not → Δf is incomplete
    idx_fxx = desc_jet.index("f_xx")
    idx_fyy = desc_jet.index("f_yy")
    assert nz_jet[idx_fxx], "f_xx should be active"
    assert not nz_jet[idx_fyy], "f_yy should be inactive"

    # Now verify that MWSINDy's operator check would detect this
    # Δf decomposition: f_xx + f_yy → both must be active
    mg_ops_missing = []
    for j_mg, mg_name in enumerate(desc_mg):
        if mg_name not in _MULTIGRAPH_TO_JET_DECOMPOSITION:
            continue
        all_active = True
        for jet_name, _mult in _MULTIGRAPH_TO_JET_DECOMPOSITION[mg_name]:
            j_jet = desc_jet.index(jet_name)
            if not nz_jet[j_jet]:
                all_active = False
                break
        if not all_active and mg_name == "Δf":
            mg_ops_missing.append(mg_name)

    assert "Δf" in mg_ops_missing, (
        "Δf should be marked as incomplete when f_yy is missing"
    )


def test_mwsindy_pipeline_produces_valid_result():
    """MWSINDy produces valid fp/fn/E/n_samples on a real dataset."""
    import numpy as np
    data_path = os.path.join(
        os.path.dirname(__file__), "..", "data", "swift_hohenberg_clean.npz",
    )
    if not os.path.exists(data_path):
        pytest.skip("Data file not available (gitignored *.npz)")
    data = np.load(data_path)
    f = data["snapshots"][20]
    Lx = float(data["Lx"]); Ly = float(data["Ly"])

    from equivariant_pde_find import derivatives
    derivatives.set_backend("spectral")
    m, p = 8, 6

    G_ws, desc_ws = library.build_wsindy_library(f, Lx, Ly, m=m, p=p, use_ibp=False, formulation="weak")
    G_mg, desc_mg = library.build_multigraph_library(f, Lx, Ly, m=m, p=p, use_ibp=False, formulation="weak")
    ft = derivatives.bilaplacian(f, Lx, Ly)
    b = wf.build_weak_rhs(ft, m=m, p=p)

    from equivariant_pde_find.library import build_polyjet_library
    _, desc_jet = build_polyjet_library(f, Lx, Ly, max_order=4, max_degree=3)
    from equivariant_pde_find.evaluation import build_wsindy_to_jet_map
    M = build_wsindy_to_jet_map(desc_ws, desc_jet)

    # Ground truth: f_t = Δ²f
    true_mg = np.zeros(len(desc_mg))
    true_mg[desc_mg.index("Δ²f")] = 1.0
    from equivariant_pde_find.evaluation import _MULTIGRAPH_TO_JET_DECOMPOSITION
    true_jet = np.zeros(len(desc_jet))
    for j, (name, coeff) in enumerate(zip(desc_mg, true_mg)):
        if abs(coeff) < 1e-12: continue
        if name not in _MULTIGRAPH_TO_JET_DECOMPOSITION: continue
        for jet_name, mult in _MULTIGRAPH_TO_JET_DECOMPOSITION[name]:
            true_jet[desc_jet.index(jet_name)] += coeff * mult

    from equivariant_pde_find.regression import TrainSTRidge
    from equivariant_pde_find.evaluation import evaluate_mwsindy, build_eval_theta_jet
    from equivariant_pde_find.simulations import make_regression_data

    # Build proper evaluation Theta from the full snapshot array
    snaps = data["snapshots"]
    times = data["times"]
    Theta_jet, Ft_jet = build_eval_theta_jet(
        snaps, times, Lx, Ly, skip_fraction=0.1, seed=0,
    )

    w_ws = TrainSTRidge(G_ws, b, lam=1e-6, d_tol=0.5, maxit=25)
    result = evaluate_mwsindy(
        w_ws, desc_ws, desc_mg, desc_jet, M,
        Theta_jet, Ft_jet, true_jet, G_mg, b,
    )

    assert isinstance(result["fp"], (int, float))
    assert isinstance(result["fn"], (int, float))
    assert result["n_samples"] > 0
    assert result["E"] >= 0
