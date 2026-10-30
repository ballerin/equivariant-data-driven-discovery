"""
Tests for ``library/spacetime.py`` (weak form in space, time, both, neither)
and the Messenger & Bortz test-function rule in ``library/weak_form.py``.
"""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives  # noqa: E402
from equivariant_pde_find.equations import find_op  # noqa: E402
from equivariant_pde_find.library import (  # noqa: E402
    build_multigraph_library, build_wsindy_library,
)
from equivariant_pde_find.library import weak_form as wf  # noqa: E402
from equivariant_pde_find.library.spacetime import (  # noqa: E402
    FORMULATIONS, build_weak_system, time_kernel, time_weak_convolve,
)
from equivariant_pde_find.regression import centering  # noqa: E402


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

NU, SOURCE = 0.5, 0.3


@pytest.fixture(scope="module")
def heat_data():
    """Exact solution of u_t = ν Δu + c on the periodic square [0, 2π)²."""
    N, L, Nt, dt = 32, 2 * math.pi, 61, 0.02
    x = np.linspace(0.0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    t = np.arange(Nt) * dt
    modes = [((1, 2), 1.0), ((3, 1), 0.6), ((2, 2), -0.4)]
    S = np.zeros((Nt, N, N))
    for (kx, ky), a in modes:
        S += (a * np.exp(-NU * (kx ** 2 + ky ** 2) * t)[:, None, None]
              * np.cos(kx * X + ky * Y)[None])
    S += SOURCE * t[:, None, None] + 1.0
    return S, t, L


def multigraph_builder(L):
    def build(f, formulation, **kw):
        return build_multigraph_library(f, L, L, formulation=formulation, **kw)
    return build


def heat_truth(desc):
    w = np.zeros(len(desc))
    w[find_op(desc, "Δf")] = NU
    w[find_op(desc, "1")] = SOURCE
    return w


def ols(X, y):
    return np.linalg.lstsq(X, y, rcond=None)[0]


# ---------------------------------------------------------------------------
# Temporal convolution
# ---------------------------------------------------------------------------


def test_time_weak_convolve_matches_np_convolve():
    rng = np.random.default_rng(0)
    A = rng.standard_normal((30, 5, 3))
    k = time_kernel(6, 9, order=1, dt=0.1)
    out = time_weak_convolve(A, k)
    assert out.shape == (30 - 12, 5, 3)
    ref = np.apply_along_axis(lambda v: np.convolve(v, k, mode="valid"), 0, A)
    np.testing.assert_allclose(out, ref, rtol=1e-12, atol=1e-12)


def test_time_weak_convolve_rejects_short_series():
    with pytest.raises(ValueError):
        time_weak_convolve(np.zeros((10, 2)), time_kernel(6, 9))


def test_time_ibp_sign_and_scale():
    """φ_t′/dt ∗ u equals φ_t ∗ u_t with no extra sign."""
    dt = 0.05
    t = np.arange(400) * dt
    u, ut = np.sin(1.3 * t), 1.3 * np.cos(1.3 * t)
    m_t = 10
    p_t = wf.select_test_function_degree(m_t, max_deriv=1)
    lhs_ibp = time_weak_convolve(u, time_kernel(m_t, p_t, 1, dt))
    lhs_fd = time_weak_convolve(ut, time_kernel(m_t, p_t, 0, dt))
    assert np.linalg.norm(lhs_ibp - lhs_fd) / np.linalg.norm(lhs_fd) < 1e-8


def test_spacetime_factorises_into_spatial_builder_then_time(heat_data):
    """Stacking per-frame spatial weak columns and convolving in time equals the
    full separable (t, x, y) convolution — for a direct and an IBP column."""
    S, t, L = heat_data
    m, p, m_t, p_t = 5, 7, 6, 9
    dx = L / S.shape[1]
    phi_x, phi_t = wf._test_function_1d(m, p), wf._test_function_1d(m_t, p_t)
    d2 = wf._test_function_deriv_1d(m, p, 2) / dx ** 2

    def build(f, formulation, **kw):
        return build_wsindy_library(f, L, L, formulation=formulation, **kw)

    with derivatives.use_backend("spectral"):
        G, _, desc, _ = build_weak_system(
            S, t, build, "spacetime_weak", m=m, p=p, space_ibp=True,
            m_t=m_t, p_t=p_t)
    n_t, n_x = S.shape[0] - 2 * m_t, S.shape[1] - 2 * m

    full_sq = wf.weak_convolution(S ** 2, [phi_t, phi_x, phi_x])
    np.testing.assert_allclose(
        G[:, desc.index("f²")].reshape(n_t, n_x, n_x), full_sq, rtol=1e-10)

    full_xx = wf.weak_convolution(S, [phi_t, d2, phi_x])
    np.testing.assert_allclose(
        G[:, desc.index("f_{xx}")].reshape(n_t, n_x, n_x), full_xx, rtol=1e-10)


# ---------------------------------------------------------------------------
# The four formulations on an exact solution
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("formulation", sorted(FORMULATIONS))
@pytest.mark.parametrize("time_ibp", [True, False])
def test_oracle_residual_on_exact_solution(heat_data, formulation, time_ibp):
    S, t, L = heat_data
    weak_t = FORMULATIONS[formulation][1]
    if not weak_t and not time_ibp:
        pytest.skip("time_ibp only applies to formulations weak in time")
    with derivatives.use_backend("spectral"):
        G, b, desc, info = build_weak_system(
            S, t, multigraph_builder(L), formulation, m=8,
            p=wf.select_test_function_degree(8, 4), space_ibp=True,
            m_t=8, p_t=wf.select_test_function_degree(8, 1), time_ibp=time_ibp)
    r = np.linalg.norm(G @ heat_truth(desc) - b) / np.linalg.norm(b)
    # Integrating by parts in time removes the O(dt²) central-FD error; with
    # (m, p) meeting the τ = 1e-10 rule the spatial quadrature is as accurate
    assert r < (1e-6 if (weak_t and time_ibp) else 5e-3), (formulation, r)
    assert info["formulation"] == formulation


def test_time_ibp_beats_fd_on_exact_solution(heat_data):
    S, t, L = heat_data
    res = {}
    with derivatives.use_backend("spectral"):
        for tibp in (True, False):
            G, b, desc, _ = build_weak_system(
                S, t, multigraph_builder(L), "time_weak",
                m_t=8, p_t=wf.select_test_function_degree(8, 1), time_ibp=tibp)
            res[tibp] = np.linalg.norm(G @ heat_truth(desc) - b) / np.linalg.norm(b)
    assert res[True] < 1e-3 * res[False]


@pytest.mark.parametrize("formulation", sorted(FORMULATIONS))
def test_recovery_with_and_without_centering(heat_data, formulation):
    """Plain least squares and the centred fit give the same coefficients,
    including the constant source written into the '1' slot."""
    S, t, L = heat_data
    with derivatives.use_backend("spectral"):
        G, b, desc, _ = build_weak_system(
            S, t, multigraph_builder(L), formulation, m=8,
            p=wf.select_test_function_degree(8, 4), space_ibp=True,
            m_t=8, p_t=wf.select_test_function_degree(8, 1))
    keep = [find_op(desc, "1"), find_op(desc, "f"), find_op(desc, "Δf")]
    G, desc = G[:, keep], [desc[j] for j in keep]
    truth = heat_truth(desc)

    w_plain = ols(G, b)
    w_cent = centering(G, b, ols, descs=desc)
    tol = 1e-5 if FORMULATIONS[formulation][1] else 5e-3
    np.testing.assert_allclose(w_plain, truth, atol=tol)
    np.testing.assert_allclose(w_cent, truth, atol=tol)
    np.testing.assert_allclose(G @ w_cent, G @ w_plain, rtol=1e-6, atol=1e-8)


# ---------------------------------------------------------------------------
# Query points, strides, masks, validation
# ---------------------------------------------------------------------------


def test_query_grid_sizes(heat_data):
    S, t, L = heat_data
    Nt, N = S.shape[0], S.shape[1]
    m, m_t = 5, 8
    kw = dict(m=m, p=9, m_t=m_t, p_t=16, spatial_stride=3, time_stride=2)
    expected_t = {
        "pointwise": len(range(0, Nt - 2, 2)),
        "space_weak": len(range(0, Nt - 2, 2)),
        "time_weak": len(range(0, Nt - 2 * m_t, 2)),
        "spacetime_weak": len(range(0, Nt - 2 * m_t, 2)),
    }
    expected_x = {
        "pointwise": len(range(0, N, 3)) ** 2,
        "time_weak": len(range(0, N, 3)) ** 2,
        "space_weak": len(range(0, N - 2 * m, 3)) ** 2,
        "spacetime_weak": len(range(0, N - 2 * m, 3)) ** 2,
    }
    with derivatives.use_backend("spectral"):
        for form in FORMULATIONS:
            G, b, _, info = build_weak_system(S, t, multigraph_builder(L),
                                              form, **kw)
            assert info["n_time_query"] == expected_t[form]
            assert info["n_space_query"] == expected_x[form]
            assert G.shape[0] == b.shape[0] == expected_t[form] * expected_x[form]


def test_point_mask_is_intersected_over_frames(heat_data):
    S, t, L = heat_data
    N = S.shape[1]
    moving = []

    def mask_fn(f):
        mask = np.zeros(f.shape, dtype=bool)
        mask[: N // 2 + len(moving) % 3] = True
        moving.append(None)
        return mask.ravel()

    with derivatives.use_backend("spectral"):
        _, _, _, info = build_weak_system(
            S[:20], t[:20], multigraph_builder(L), "pointwise",
            point_mask_fn=mask_fn)
    assert info["n_space_query"] == (N // 2) * N


def test_validation_errors(heat_data):
    S, t, L = heat_data
    build = multigraph_builder(L)
    with pytest.raises(ValueError, match="Unknown formulation"):
        build_weak_system(S, t, build, "weak")
    with pytest.raises(ValueError, match="uniformly spaced"):
        t_bad = t.copy()
        t_bad[5] += 1e-3
        build_weak_system(S, t_bad, build, "pointwise")
    with pytest.raises(ValueError, match="needs spatial m"):
        build_weak_system(S, t, build, "space_weak")
    with pytest.raises(ValueError, match="needs temporal"):
        build_weak_system(S, t, build, "time_weak")
    with pytest.raises(ValueError, match="frames"):
        build_weak_system(S, t, build, "time_weak", m_t=40, p_t=10)


# ---------------------------------------------------------------------------
# Messenger & Bortz rule
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("m_x, m_t, p_x, p_t", [
    # arXiv:2007.02848, Table 3 (m_x, m_t) -> Table 4 (p_x, p_t), tau = 1e-10
    (60, 60, 7, 7),    # Burgers
    (45, 80, 8, 7),    # KdV
    (23, 22, 10, 10),  # KS
    (19, 25, 11, 10),  # NLS
    (40, 25, 8, 10),   # SG
    (13, 14, 13, 12),  # RD
    (31, 14, 9, 12),   # NS
])
def test_degree_rule_reproduces_paper_tables(m_x, m_t, p_x, p_t):
    assert wf._compute_p_from_m(m_x, 6, 1e-10) == p_x
    assert wf._compute_p_from_m(m_t, 1, 1e-10) == p_t


@pytest.mark.parametrize("k_star, N, tau_hat", [
    (24, 256, 3.0), (5, 64, 2.0), (1, 64, 1.0), (40, 128, 1.0), (3, 41, 2.0),
])
def test_root_m_solves_F(k_star, N, tau_hat):
    tau = 1e-10
    m = wf.mb_root_m(k_star, N, tau, tau_hat)
    assert 2.0 <= m <= N / 2 - 1
    if m < N / 2 - 1:
        F = (math.log((2 * m - 1) / m ** 2)
             * (4 * math.pi ** 2 * k_star ** 2 * m ** 2 - 3 * N ** 2 * tau_hat ** 2)
             - 2 * N ** 2 * tau_hat ** 2 * math.log(tau))
        assert abs(F) < 1e-6 * N ** 2 * tau_hat ** 2
        # The unique root lies in the paper's bracket
        lo = math.sqrt(3) * (N / 2) * tau_hat / (math.pi * k_star)
        assert m >= lo - 1e-9


def test_root_m_matches_paper_ks_width():
    """Paper's KS example (N=256, τ̂=3): the rule's Gaussian-equivalent width
    m/√(2p+3) equals τ̂N/(2πk*) — (C2) of wsindy-test-function-rule.md."""
    k_star, N, tau_hat = 24, 256, 3.0
    m = wf.mb_root_m(k_star, N, 1e-10, tau_hat)
    p = math.log(1e-10) / math.log((2 * m - 1) / m ** 2)
    assert m / math.sqrt(2 * p + 3) == pytest.approx(
        tau_hat * N / (2 * math.pi * k_star), rel=1e-6)


def test_critical_wavenumber_finds_noise_floor():
    rng = np.random.default_rng(1)
    N, Nt = 128, 50
    x = np.arange(N) * 2 * math.pi / N
    U = np.zeros((Nt, N))
    for k in range(1, 21):
        U += (0.8 ** k) * np.cos(k * x[None] + rng.uniform(0, 2 * math.pi, (Nt, 1)))
    U += 1e-3 * rng.standard_normal(U.shape)
    k_star = wf.mb_critical_wavenumber(U, axis=1)
    # 0.8^k meets the 1e-3·√N-scaled noise floor near k ≈ 20; the reference
    # implementation reports 3 modes below the corner
    assert 10 <= k_star <= 22


def test_mb_params_respect_bounds():
    rng = np.random.default_rng(2)
    U = rng.standard_normal((41, 24, 30))       # pure noise: k* → 1
    out = wf.mb_test_function_params(U, spatial_axes=(1, 2), time_axis=0)
    assert 2 <= out["m"] <= (24 - 1) // 2
    assert 2 <= out["m_t"] <= (41 - 1) // 2
    assert out["p"] >= 5 and out["p_t"] >= 2
    assert out["p"] == wf._compute_p_from_m(out["m"], 4)
    assert out["p_t"] == wf._compute_p_from_m(out["m_t"], 1)
    assert len(out["k_x"]) == 2
    no_t = wf.mb_test_function_params(U, spatial_axes=(1, 2))
    assert "m_t" not in no_t and no_t["m"] == out["m"]
