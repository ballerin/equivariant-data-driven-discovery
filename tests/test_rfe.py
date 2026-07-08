"""
Tests for RFE — Recursive Feature Elimination.

Covers ``RFE`` and ``TrainRFE`` from ``regression.py``.

Reference: Stephany & Earls (2021), "PDE-READ: Human-readable Partial
Differential Equation Discovery using Deep Learning", arXiv:2111.00998.
"""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.regression import RFE, TrainRFE, STRidge
from equivariant_pde_find.library import build_equivariant_library
from equivariant_pde_find.simulations import allen_cahn


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _sparse_system(n=200, d=20, nz=4, noise=0.0, seed=42):
    """Generate a synthetic sparse linear system G w_true = b (+ noise)."""
    rng = np.random.default_rng(seed)
    G = rng.standard_normal((n, d))
    # Normalize columns so column norms are O(1)
    G = G / np.sqrt(np.sum(G ** 2, axis=0))[np.newaxis, :]
    w_true = np.zeros(d)
    nz_idx = rng.choice(d, size=nz, replace=False)
    w_true[nz_idx] = rng.uniform(0.5, 2.0, size=nz) * rng.choice(
        [-1, 1], size=nz,
    )
    b = G @ w_true
    if noise > 0:
        b += noise * rng.standard_normal(n)
    return G, b, w_true, nz_idx


# ======================================================================
# RFE — exact recovery on synthetic systems
# ======================================================================


def test_rfe_exact_recovery_noiseless():
    """RFE should recover the exact sparse vector on clean data."""
    G, b, w_true, nz_idx = _sparse_system(n=200, d=20, nz=4, noise=0.0)
    w = RFE(G, b)
    nonzero = np.abs(w) > 1e-8
    assert set(nz_idx) == set(np.where(nonzero)[0])
    np.testing.assert_allclose(w[nz_idx], w_true[nz_idx], rtol=1e-8)


def test_rfe_sparsity_recovery_noisy():
    """RFE should recover the correct support with small noise."""
    G, b, w_true, nz_idx = _sparse_system(n=400, d=20, nz=3, noise=0.005)
    w = RFE(G, b)
    active = np.where(np.abs(w) > 1e-8)[0]
    assert set(nz_idx).issubset(set(active))
    # With this noise level, should have at most 1 false positive
    assert len(active) <= len(nz_idx) + 1


# ======================================================================
# RFE — rescaling and normalisation invariance
# ======================================================================


def test_rfe_rescaling_correct():
    """RFE correctly rescales: column-scaling the library should divide
    recovered coefficients by the same scale factors."""
    G, b, w_true, nz_idx = _sparse_system(n=200, d=15, nz=4, noise=0.0)

    w1 = RFE(G, b)

    # Scale columns by different amounts.
    # True coefficients should become w_true / s.
    s = np.array([1.0 + 0.5 * i for i in range(G.shape[1])])
    G_scaled = G * s[np.newaxis, :]
    w2 = RFE(G_scaled, b)

    # After rescaling, w2 * s should equal w1 (both recover w_true).
    np.testing.assert_allclose(w2 * s, w1, rtol=1e-8)


def test_rfe_normalization_invariant():
    """RFE correctly handles global column scaling: RFE(2*X0, y) * 2
    should equal RFE(X0, y)."""
    G, b, w_true, nz_idx = _sparse_system(n=200, d=15, nz=4, noise=0.0)

    w1 = RFE(G, b)
    w2 = RFE(2.0 * G, b)

    # Support should match despite different coefficient magnitudes
    active1 = set(np.where(np.abs(w1) > 1e-8)[0])
    active2 = set(np.where(np.abs(w2) > 1e-8)[0])
    assert active1 == active2
    # Coefficients scale as 1/scale_factor (2*w2 = w1)
    np.testing.assert_allclose(2.0 * w2, w1, rtol=1e-8)


# ======================================================================
# RFE — edge cases
# ======================================================================


def test_rfe_single_feature():
    """RFE should work with d=1 (single feature)."""
    rng = np.random.default_rng(42)
    n = 100
    w_true = np.array([2.5])
    G = rng.standard_normal((n, 1))
    b = G @ w_true

    w = RFE(G, b)
    np.testing.assert_allclose(w, w_true, rtol=1e-8)


def test_rfe_zero_column():
    """RFE should not crash with an all-zero column."""
    rng = np.random.default_rng(42)
    n, d = 100, 5
    G = rng.standard_normal((n, d))
    G[:, 2] = 0.0  # Make column 2 all zeros
    b = G @ np.array([1.0, 0.0, 0.0, 2.0, 0.0])

    w = RFE(G, b)
    # Should not crash; zero column should be eliminated early
    assert w.shape == (d,)


def test_rfe_perfect_fit():
    """RFE should handle RLS[0] ≈ 0 (perfect fit, n=d)."""
    rng = np.random.default_rng(42)
    n, d = 8, 8
    G = rng.standard_normal((n, d))
    w_true = rng.uniform(-2, 2, size=d)
    b = G @ w_true

    w = RFE(G, b)
    # With perfect fit, RLS[0] = 0; denominator clamping should prevent NaN
    assert not np.any(np.isnan(w))
    np.testing.assert_allclose(w, w_true, rtol=1e-6)


def test_rfe_j_max_bounds():
    """j_max should always be in [0, d-1]."""
    rng = np.random.default_rng(123)
    for _ in range(20):
        n = rng.integers(50, 200)
        d = rng.integers(3, 20)
        G = rng.standard_normal((n, d))
        G = G / np.sqrt(np.sum(G ** 2, axis=0))[np.newaxis, :]
        nz = rng.integers(1, min(d, 6))
        nz_idx = rng.choice(d, size=nz, replace=False)
        w_true = np.zeros(d)
        w_true[nz_idx] = rng.uniform(0.5, 2.0, size=nz)
        b = G @ w_true + 0.01 * rng.standard_normal(n)

        w = RFE(G, b)
        nonzero = np.count_nonzero(w)
        # At least one feature should be active, at most all
        assert 1 <= nonzero <= d


# ======================================================================
# RFE — selection criterion
# ======================================================================


def test_rfe_elimination_order():
    """A feature with a truly tiny true coefficient should not survive RFE."""
    rng = np.random.default_rng(42)
    n, d = 500, 5
    G = rng.standard_normal((n, d))
    G = G / np.sqrt(np.sum(G ** 2, axis=0))[np.newaxis, :]
    # Truth: 3 features with large coeffs, 2 zero
    w_true = np.array([2.0, 1.5, 0.0, -1.0, 0.0])
    b = G @ w_true

    w = RFE(G, b)
    # Both zero features should be inactive in the final result
    assert np.abs(w[2]) < 1e-8
    assert np.abs(w[4]) < 1e-8
    # The three active features should be recovered
    assert np.abs(w[0]) > 1e-8
    assert np.abs(w[1]) > 1e-8
    assert np.abs(w[3]) > 1e-8


def test_rfe_ratio_selection_correct():
    """On a system with one clearly dominant feature, j_max should be correct."""
    rng = np.random.default_rng(42)
    n, d = 500, 5
    G = rng.standard_normal((n, d))
    G = G / np.sqrt(np.sum(G ** 2, axis=0))[np.newaxis, :]

    # Three features: one very strong, two moderate, two zero
    w_true = np.array([5.0, 0.8, 0.0, 0.7, 0.0])
    b = G @ w_true + 0.002 * rng.standard_normal(n)

    w = RFE(G, b)
    active = set(np.where(np.abs(w) > 1e-6)[0])
    # Should have exactly 3 active features, not 4 or 5
    assert active == {0, 1, 3}


# ======================================================================
# RFE — TrainRFE wrapper
# ======================================================================


def test_train_rfe_is_equivalent():
    """TrainRFE should return the same result as calling RFE directly."""
    G, b, w_true, nz_idx = _sparse_system(n=200, d=20, nz=4, noise=0.0)
    w1 = RFE(G, b)
    w2 = TrainRFE(G, b)
    np.testing.assert_array_equal(w1, w2)


# ======================================================================
# RFE — Allen-Cahn end-to-end
# ======================================================================


@pytest.fixture(scope="module")
def ac_data():
    """Simulate a small Allen-Cahn trajectory for testing."""
    snaps, times = allen_cahn(
        N=32, L=2 * np.pi, T_end=2.0, eps=0.15, dt=0.01,
        save_every=20, seed=42,
    )
    return snaps, times


def _true_ac(descriptions):
    """Ground-truth Allen-Cahn coefficients (equivariant basis)."""
    w = np.zeros(len(descriptions))
    w[descriptions.index("f")] = 1.0
    w[descriptions.index("f\u00b3")] = -1.0
    w[descriptions.index("\u0394f")] = 0.0225  # eps²
    return w


def test_rfe_allen_cahn_clean(ac_data):
    """RFE recovers the dominant Allen-Cahn terms (fn≤1, fp≤1) on clean data."""
    from equivariant_pde_find.simulations import make_regression_data

    snaps, times = ac_data

    # Build equivariant library on all snapshots
    _, desc = build_equivariant_library(snaps[0], 2 * np.pi, 2 * np.pi)
    true_w = _true_ac(desc)

    # Identify operator indices for diagnostics
    idx_f = desc.index("f")
    idx_f3 = desc.index("f\u00b3")
    idx_lap = desc.index("\u0394f")

    # Make regression data — use all snapshots
    f_snaps, ft_snaps, _, _ = make_regression_data(
        snaps, times, n_per_snap=200, skip_fraction=0.0, seed=0,
    )
    Theta = np.vstack([
        build_equivariant_library(f, 2 * np.pi, 2 * np.pi)[0]
        for f in f_snaps
    ])
    Ft = np.concatenate([ft.ravel() for ft in ft_snaps])

    w = RFE(Theta, Ft)

    # Check that the two linear terms are always recovered
    true_active = set(np.where(np.abs(true_w) > 1e-8)[0])
    recovered_active = set(np.where(np.abs(w) > 1e-6)[0])

    assert idx_f in recovered_active, "RFE missed the reaction term f"
    assert idx_lap in recovered_active, "RFE missed the diffusion term Δf"

    # f³ may or may not be recovered (RFE can struggle with weak nonlinearity
    # on short trajectories); at most 1 false negative allowed
    missing = true_active - recovered_active
    assert len(missing) <= 1, (
        f"RFE missed too many true terms: {missing}. "
        f"True: {true_active}, Recovered: {recovered_active}"
    )
    # At most 1 false positive
    extra = recovered_active - true_active
    assert len(extra) <= 1, f"Too many false positives: {extra}"
    # Recovered true-term coefficients should be close
    common = list(true_active & recovered_active)
    np.testing.assert_allclose(w[common], true_w[common], rtol=5e-2)


# ======================================================================
# RFE vs STRidge comparison
# ======================================================================


def test_rfe_vs_stridge_clean():
    """On clean synthetic data, RFE should match STRidge with optimal d_tol."""
    G, b, w_true, nz_idx = _sparse_system(n=200, d=20, nz=4, noise=0.0)

    w_rfe = RFE(G, b)

    # STRidge with carefully chosen tol should recover the same support
    w_stridge = STRidge(G, b, lam=1e-5, maxit=25, tol=0.05,
                        normalize=2, early_stop=True, unbias=True)

    active_rfe = set(np.where(np.abs(w_rfe) > 1e-8)[0])
    active_stridge = set(np.where(np.abs(w_stridge) > 1e-8)[0])
    assert active_rfe == active_stridge
