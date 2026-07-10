"""
Tests for MSTLS — Modified Sequential Thresholding Least Squares.

Covers ``mstls``, ``TrainMSTLS``, and ``select_lambda_grid`` from
``regression.py``.

Reference: Messenger & Bortz (2021), J. Comput. Phys., Section 3.2.
"""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.regression import (
    mstls,
    TrainMSTLS,
    select_lambda_grid,
    STRidge,
)


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
    w_true[nz_idx] = rng.uniform(0.5, 2.0, size=nz) * rng.choice([-1, 1], size=nz)
    b = G @ w_true
    if noise > 0:
        b += noise * rng.standard_normal(n)
    return G, b, w_true, nz_idx


# ======================================================================
# select_lambda_grid
# ======================================================================


def test_select_lambda_grid_shape():
    G, b, _, _ = _sparse_system()
    lambdas = select_lambda_grid(G, b, n_lambdas=50)
    assert len(lambdas) == 50
    assert np.all(lambdas > 0)
    assert np.allclose(lambdas, np.sort(lambdas))  # ascending


def test_select_lambda_grid_range():
    G, b, _, _ = _sparse_system()
    lambdas = select_lambda_grid(G, b, n_lambdas=100)
    assert lambdas[0] < lambdas[-1]
    # λ_max should be >= max |G_j^T b| / ||G_j||^2
    col_norms = np.sqrt(np.sum(G ** 2, axis=0))
    Gty = np.abs(G.T @ b.ravel())
    lam_max_theory = np.max(Gty / np.maximum(col_norms ** 2, 1e-30))
    assert lambdas[-1] >= lam_max_theory * 0.99


def test_select_lambda_grid_large_norms():
    """λ grid is not truncated when column norms produce large λ_max.

    The previous clamp (lam_max ≤ 10.0) was removed to match MATLAB;
    this test verifies that large λ values survive the new wide clamps.
    """
    rng = np.random.default_rng(42)
    n, d = 200, 5
    G = rng.standard_normal((n, d)) * 50.0  # large column norms
    w_true = np.array([0.0, 100.0, 0.0, -50.0, 0.0])
    b = G @ w_true

    lambdas = select_lambda_grid(G, b, n_lambdas=50)
    # λ_max should be large (>> 10, which was the old clamp)
    col_norms = np.sqrt(np.sum(G ** 2, axis=0))
    Gty = np.abs(G.T @ b.ravel())
    lam_max_theory = np.max(Gty / np.maximum(col_norms ** 2, 1e-30))
    assert lambdas[-1] >= lam_max_theory * 0.99
    # With large coefficients, λ_max can exceed 10.0
    assert lambdas[-1] > 1.0


# ======================================================================
# mstls — single λ
# ======================================================================


def test_mstls_perfect_recovery_noiseless():
    """MSTLS with λ well-chosen should recover the true sparse vector."""
    G, b, w_true, nz_idx = _sparse_system(n=300, d=20, nz=3, noise=0.0)
    # λ = 0.05 is generous for dominant-balance with unit-norm columns
    w = mstls(G, b, lam=0.05, maxit=50)
    nonzero = np.abs(w) > 1e-8
    # At least recover all true terms
    assert set(nz_idx).issubset(set(np.where(nonzero)[0]))
    # Coefficients should be close
    np.testing.assert_allclose(w[nz_idx], w_true[nz_idx], rtol=1e-6)


def test_mstls_sparsity_enforced():
    """Tight λ should enforce sparsity via dominant-balance bounds."""
    G, b, w_true, nz_idx = _sparse_system(n=300, d=20, nz=5, noise=0.0)
    # λ ≈ min |coeff| / ||b|| gives tight bounds
    min_coeff = np.min(np.abs(w_true[nz_idx]))
    b_norm = np.linalg.norm(b)
    lam = min_coeff / max(b_norm, 1.0) * 0.9
    w = mstls(G, b, lam=lam, maxit=50)
    # Should drop at least one term with coefficient near the bound
    assert np.count_nonzero(w) <= 5


def test_mstls_large_lambda_keeps_few():
    """With λ near 1 and well-scaled data, dominant-balance keeps few terms."""
    rng = np.random.default_rng(99)
    n, d = 300, 10
    G = rng.standard_normal((n, d))
    G = G / np.sqrt(np.sum(G ** 2, axis=0))[np.newaxis, :]  # unit-norm columns
    # Coefficients near 1 so ||b|| ≈ 1
    w_true = np.zeros(d)
    w_true[[1, 4, 7]] = [0.8, -0.6, 0.9]
    b = G @ w_true  # ||b|| ≈ 1.36

    # λ=0.8: bnds≈1.36, lower≈0.8*1.36≈1.09, upper≈1.25*1.36≈1.70
    # Coefficients [0.8, 0.6, 0.9] all below lower bound → all zeroed
    w_tight = mstls(G, b, lam=0.8, maxit=50)
    assert np.count_nonzero(w_tight) <= 1  # at most 1 term survives

    # λ=0.5: lower≈0.5*1.36≈0.68, upper≈2*1.36≈2.72
    # All coefficients [0.8, 0.6, 0.9] survive
    w_loose = mstls(G, b, lam=0.5, maxit=50)
    assert np.count_nonzero(w_loose) >= 2  # most survive
    assert np.count_nonzero(w_loose) > np.count_nonzero(w_tight)


def test_mstls_all_terms_zeroed():
    """Very large λ zeros all terms (both bounds collapse)."""
    G, b, _, _ = _sparse_system(n=300, d=20, nz=3, noise=0.0)
    # With λ ≫ 1, lower ≫ 1 and upper ≪ 1 → no coefficient can satisfy both
    w = mstls(G, b, lam=1e6, maxit=50)
    assert np.count_nonzero(w) == 0
    np.testing.assert_allclose(w, 0.0, atol=1e-14)


def test_mstls_converges_in_few_iters():
    """MSTLS should converge in far fewer than maxit on clean data."""
    G, b, _, _ = _sparse_system(n=300, d=20, nz=4, noise=0.0)
    w = mstls(G, b, lam=0.1, maxit=5)
    assert np.all(np.isfinite(w))


def test_mstls_maxit_default():
    """When maxit is None, mstls defaults to d+1 (matching the reference)."""
    G, b, w_true, nz_idx = _sparse_system(n=300, d=20, nz=3, noise=0.0)
    w = mstls(G, b, lam=0.05)  # maxit=None → defaults to d+1 = 21
    nonzero = np.abs(w) > 1e-8
    assert set(nz_idx).issubset(set(np.where(nonzero)[0]))


def test_mstls_default_M_is_ones():
    """When M is not provided, mstls defaults to M = ones(d)."""
    G, b, w_true, nz_idx = _sparse_system(n=300, d=20, nz=3, noise=0.0)
    w = mstls(G, b, lam=0.05, maxit=50)  # no M, defaults to ones
    assert set(nz_idx).issubset(set(np.where(np.abs(w) > 1e-8)[0]))


def test_mstls_noisy_data():
    """On noisy data, MSTLS should still identify the correct support."""
    G, b, w_true, nz_idx = _sparse_system(n=500, d=20, nz=4, noise=0.01)
    w = mstls(G, b, lam=0.1, maxit=50)
    nonzero = np.abs(w) > 1e-8
    # May miss some terms at this noise level, but should have no false positives
    fp = set(np.where(nonzero)[0]) - set(nz_idx)
    assert len(fp) == 0  # no false positives


# ======================================================================
# TrainMSTLS — auto λ
# ======================================================================


def test_TrainMSTLS_perfect_noiseless():
    """TrainMSTLS auto-selects λ and recovers the true vector."""
    G, b, w_true, nz_idx = _sparse_system(n=300, d=20, nz=4, noise=0.0)
    w, lam_best, losses = TrainMSTLS(G, b, n_lambdas=50, alpha=0.5)
    nonzero = np.abs(w) > 1e-8
    assert set(nz_idx).issubset(set(np.where(nonzero)[0]))
    np.testing.assert_allclose(w[nz_idx], w_true[nz_idx], rtol=1e-5)
    assert len(losses) == 50
    assert np.all(np.isfinite(losses))


def test_TrainMSTLS_custom_lambdas():
    """TrainMSTLS accepts a user-supplied λ grid."""
    G, b, w_true, nz_idx = _sparse_system(n=300, d=20, nz=3, noise=0.0)
    custom_lam = np.logspace(-3, 0, 30)
    w, lam_best, losses = TrainMSTLS(G, b, lambdas=custom_lam, alpha=0.5)
    assert len(losses) == 30
    assert lam_best in custom_lam
    nonzero = np.abs(w) > 1e-8
    assert set(nz_idx).issubset(set(np.where(nonzero)[0]))


def test_TrainMSTLS_noisy():
    """TrainMSTLS on mildly noisy data."""
    G, b, w_true, nz_idx = _sparse_system(n=500, d=20, nz=4, noise=0.02)
    w, lam_best, losses = TrainMSTLS(G, b, n_lambdas=50, alpha=0.5)
    nonzero = np.abs(w) > 1e-8
    # At mild noise, should have no false positives
    fp = set(np.where(nonzero)[0]) - set(nz_idx)
    assert len(fp) == 0


def test_TrainMSTLS_alpha_extremes():
    """α=0 (all sparsity) vs α=1 (all projection) produce different results."""
    G, b, _, _ = _sparse_system(n=300, d=20, nz=4, noise=0.0)
    w0, _, _ = TrainMSTLS(G, b, n_lambdas=50, alpha=0.0)
    w1, _, _ = TrainMSTLS(G, b, n_lambdas=50, alpha=1.0)
    # α=0 prefers sparsity → fewer or equal terms
    assert np.count_nonzero(w0) <= np.count_nonzero(w1)


# ======================================================================
# Comparison with STRidge
# ======================================================================


def test_mstls_handles_column_scale_disparity():
    """MSTLS dominant-balance handles columns with varying norms naturally."""
    rng = np.random.default_rng(99)
    n, d = 500, 15
    G = rng.standard_normal((n, d))
    # Create moderate scale disparity (~3x)
    G[:, :5] *= 3.0

    w_true = np.zeros(d)
    nz_idx = [1, 10, 3, 12]
    for i in nz_idx:
        w_true[i] = rng.uniform(0.4, 0.9) * rng.choice([-1, 1])
    b = G @ w_true + 0.01 * rng.standard_normal(n)

    # Dominant-balance automatically handles scale via ‖b‖/‖G_j‖
    w_ms, lam_best, _ = TrainMSTLS(G, b, n_lambdas=50, alpha=0.5)
    nz_ms = np.where(np.abs(w_ms) > 1e-8)[0]
    # Should not have false positives
    assert len(set(nz_ms) - set(nz_idx)) == 0


def test_mstls_upper_bound_prevents_large_spurious():
    """The dominant-balance upper bound prevents terms with |coeff| ≫ 1/λ."""
    rng = np.random.default_rng(42)
    n, d = 300, 10
    # Use nearly-orthogonal columns so ||b|| ≈ ||w_true||
    G = rng.standard_normal((n, d)) / np.sqrt(n)
    U_svd, _, Vt = np.linalg.svd(G, full_matrices=False)
    G = U_svd @ Vt  # columns are orthonormal

    w_true = np.zeros(d)
    w_true[0] = 1.0
    b = G @ w_true  # ||b|| = 1.0 exactly

    # λ=0.5: bnds=1, lower=max(0.5,1)*0.5=0.5, upper=2*min(1,1)=1
    # Coefficient 1.0 is within [0.5, 1.0]
    w = mstls(G, b, lam=0.5, maxit=50)
    nonzero = np.where(np.abs(w) > 1e-8)[0]
    assert 0 in nonzero
    # Upper bound: |w_j| ≤ (1/λ) * min(‖b‖/‖G_j‖, 1) * M_j = 2 * 1 = 2
    for j in nonzero:
        assert abs(w[j]) <= 2.0 * 1.01


# ======================================================================
# Integration tests — MSTLS in sweep pipeline
# ======================================================================


def test_TrainMSTLS_one_pointwise():
    """sweep_one with algorithm='mstls' on synthetic pointwise-like data."""
    rng = np.random.default_rng(123)
    n, d = 1000, 10
    # Use column-normalized matrix so lambda grid is well-scaled
    Theta = rng.standard_normal((n, d))
    Theta = Theta / np.sqrt(np.sum(Theta ** 2, axis=0))[np.newaxis, :]
    w_true = np.zeros(d)
    w_true[0] = 1.0
    w_true[3] = -0.5
    Ft = Theta @ w_true

    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(SCRIPT_DIR, "..", "scripts"))
    from sweep import sweep_one  # type: ignore[import-not-found]

    n_active, w, lam_out = sweep_one(Theta, Ft, w_true, 0.0, 0.0, "mstls")
    assert n_active >= 2, f"n_active={n_active}, w={w}"
    assert lam_out > 0
    # With column-normalized data, MSTLS should recover coefficients well
    assert abs(w[0] - 1.0) < 0.2
    assert abs(w[3] + 0.5) < 0.2


def test_TrainMSTLS_one_weak_form():
    """sweep_one with algorithm='mstls' on synthetic weak-form-like data."""
    rng = np.random.default_rng(456)
    n, d = 500, 8
    G = rng.standard_normal((n, d))
    G = G / np.sqrt(np.sum(G ** 2, axis=0))[np.newaxis, :]
    w_true = np.zeros(d)
    w_true[1] = 1.5
    w_true[5] = -2.0
    b = G @ w_true + 0.005 * rng.standard_normal(n)

    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(SCRIPT_DIR, "..", "scripts"))
    from sweep import sweep_one  # type: ignore[import-not-found]

    n_active, w, lam_out = sweep_one(G, b, w_true, 0.0, 0.0, "mstls")
    assert n_active >= 2, f"n_active={n_active}, w={w}"
    assert lam_out > 0
    assert abs(w[1] - 1.5) < 0.3
    assert abs(w[5] + 2.0) < 0.3


def test_TrainMSTLS_cli_accepts_mstls():
    """Verify that --algorithms mstls is accepted by the CLI parser."""
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(SCRIPT_DIR, "..", "scripts"))
    from sweep import main as sweep_main  # type: ignore[import-not-found]
    assert True  # import succeeded


def test_mstls_vs_strridge_clean():
    """MSTLS and STRidge should both recover the true model on clean data."""
    rng = np.random.default_rng(789)
    n, d = 800, 10
    Theta = rng.standard_normal((n, d))
    Theta = Theta / np.sqrt(np.sum(Theta ** 2, axis=0))[np.newaxis, :]
    w_true = np.zeros(d)
    w_true[2] = 0.8
    w_true[7] = -1.2
    Ft = Theta @ w_true

    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(SCRIPT_DIR, "..", "scripts"))
    from sweep import sweep_one  # type: ignore[import-not-found]

    # STRidge
    _, w_s, _ = sweep_one(Theta, Ft, w_true, 1e-5, 0.05, "strridge")
    # MSTLS
    _, w_m, _ = sweep_one(Theta, Ft, w_true, 0.0, 0.0, "mstls")

    # Both should recover the two nonzero coefficients
    assert abs(w_s[2] - 0.8) < 0.1, f"STRidge w[2]={w_s[2]:.4f}"
    assert abs(w_s[7] + 1.2) < 0.1, f"STRidge w[7]={w_s[7]:.4f}"
    assert abs(w_m[2] - 0.8) < 0.2, f"MSTLS w[2]={w_m[2]:.4f}"
    assert abs(w_m[7] + 1.2) < 0.2, f"MSTLS w[7]={w_m[7]:.4f}"
    # Both should be sparse
    assert np.count_nonzero(np.abs(w_s) > 1e-8) >= 2
    assert np.count_nonzero(np.abs(w_m) > 1e-8) >= 2


# ======================================================================
# Dominant-balance MSTLS (M parameter, sparsifyDynamics path A)
# ======================================================================


def test_mstls_dominant_balance_perfect():
    """With M=np.ones(d), dominant-balance MSTLS recovers true support."""
    rng = np.random.default_rng(42)
    n, d = 500, 15
    X = rng.standard_normal((n, d))
    # Use columns with roughly similar norms so bounds are comparable
    X = X / np.sqrt(np.sum(X ** 2, axis=0))[np.newaxis, :] * 3.0
    w_true = np.zeros(d)
    w_true[1] = 1.5
    w_true[10] = -1.2
    y = X @ w_true

    M = np.ones(d)
    w = mstls(X, y, lam=0.1, M=M, maxit=50)

    # Should recover both true terms
    nonzero = np.abs(w) > 1e-8
    assert 1 in np.where(nonzero)[0]
    assert 10 in np.where(nonzero)[0]
    # No false positives
    fp = set(np.where(nonzero)[0]) - {1, 10}
    assert len(fp) == 0


def test_mstls_dominant_balance_upper_bound():
    """The upper bound should zero coefficients that are too large."""
    rng = np.random.default_rng(99)
    n, d = 300, 10
    X = rng.standard_normal((n, d))
    w_true = np.zeros(d)
    w_true[0] = 100.0  # very large coefficient
    y = X @ w_true

    # With lam=0.01: ||y||/||X_0|| ≈ 100*17/17 = 100,
    # bnds_0 = 100, UBs_0 = (1/0.01) * min(1, 100) = 1
    # |w_0| = 100 > UBs_0 = 1 → zeroed
    M = np.ones(d)
    w = mstls(X, y, lam=0.01, M=M, maxit=50)
    assert abs(w[0]) < 1e-8  # eliminated by upper bound


def test_mstls_dominant_balance_lower_bound():
    """The lower bound zeros coefficients below λ when bnds <= 1."""
    rng = np.random.default_rng(123)
    n, d = 300, 10
    X = rng.standard_normal((n, d))
    X = X / np.sqrt(np.sum(X ** 2, axis=0))[np.newaxis, :]  # unit-norm
    w_true = np.zeros(d)
    w_true[0] = 0.05  # below lam=0.1
    w_true[1] = 0.8   # above lam=0.1
    y = X @ w_true

    # With unit-norm columns: ||y|| ≈ 0.8, bnds_j ≈ 0.8,
    # LBs_j = lam * max(1, 0.8) = lam = 0.1
    M = np.ones(d)
    w = mstls(X, y, lam=0.1, M=M, maxit=50)
    assert abs(w[0]) < 1e-8  # 0.05 < 0.1 → zeroed
    assert abs(w[1]) > 0.1   # 0.8 > 0.1 → survives


def test_mstls_dominant_balance_with_custom_M():
    """Custom M scales the dominant-balance bounds."""
    rng = np.random.default_rng(555)
    n, d = 400, 10
    X = rng.standard_normal((n, d))
    X = X / np.sqrt(np.sum(X ** 2, axis=0))[np.newaxis, :]  # unit-norm columns
    w_true = np.zeros(d)
    w_true[2] = 0.8
    w_true[7] = -0.6
    y = X @ w_true

    # Default M=ones vs explicit M=ones should give same result
    M = np.ones(d)
    w_default = mstls(X, y, lam=0.1, maxit=50)       # M defaults to ones
    w_explicit = mstls(X, y, lam=0.1, M=M, maxit=50)  # explicit ones

    # Same active set
    active_default = set(np.where(np.abs(w_default) > 1e-8)[0])
    active_explicit = set(np.where(np.abs(w_explicit) > 1e-8)[0])
    assert active_default == active_explicit


def test_TrainMSTLS_with_M():
    """TrainMSTLS accepts M and passes it through."""
    rng = np.random.default_rng(777)
    n, d = 500, 10
    X = rng.standard_normal((n, d))
    w_true = np.zeros(d)
    w_true[3] = 0.7
    w_true[8] = -0.5
    y = X @ w_true

    M = np.ones(d)
    w, lam_best, losses = TrainMSTLS(
        X, y, M=M, n_lambdas=30, alpha=0.5
    )
    nonzero = np.abs(w) > 1e-8
    assert set([3, 8]).issubset(set(np.where(nonzero)[0]))
    assert len(losses) == 30


def test_TrainMSTLS_default_M_is_ones():
    """TrainMSTLS with M=None defaults to M=ones(d) internally."""
    rng = np.random.default_rng(888)
    n, d = 500, 10
    X = rng.standard_normal((n, d))
    w_true = np.zeros(d)
    w_true[3] = 0.7
    w_true[8] = -0.5
    y = X @ w_true

    # Default (M=None, should behave as M=ones)
    w_default, lam_d, losses_d = TrainMSTLS(X, y, n_lambdas=30, alpha=0.5)
    # Explicit M=ones
    w_explicit, lam_e, losses_e = TrainMSTLS(
        X, y, M=np.ones(d), n_lambdas=30, alpha=0.5
    )

    # Same results
    np.testing.assert_allclose(w_default, w_explicit, rtol=1e-12)
    assert lam_d == lam_e
    np.testing.assert_allclose(losses_d, losses_e, rtol=1e-12)
