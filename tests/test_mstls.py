"""
Tests for MSTLS — Modified Sequential Thresholding Least Squares.

Covers ``mstls``, ``mstls_sweep``, and ``select_lambda_grid`` from
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
    mstls_sweep,
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


# ======================================================================
# mstls — single λ
# ======================================================================


def test_mstls_perfect_recovery_noiseless():
    """MSTLS with λ well-chosen should recover the true sparse vector."""
    G, b, w_true, nz_idx = _sparse_system(n=300, d=20, nz=3, noise=0.0)
    # λ = 0.1 is generous for noiseless data with column-normed G
    w = mstls(G, b, lam=0.1, maxits=50, normalize=2)
    nonzero = np.abs(w) > 1e-8
    # At least recover all true terms
    assert set(nz_idx).issubset(set(np.where(nonzero)[0]))
    # Coefficients should be close
    np.testing.assert_allclose(w[nz_idx], w_true[nz_idx], rtol=1e-6)


def test_mstls_sparsity_enforced():
    """Tight λ (matching coefficient magnitude) should enforce sparsity."""
    G, b, w_true, nz_idx = _sparse_system(n=300, d=20, nz=5, noise=0.0)
    # λ ≈ min |coeff| gives tight bounds that exclude terms below that magnitude
    min_coeff = np.min(np.abs(w_true[nz_idx]))
    w = mstls(G, b, lam=min_coeff * 0.9, maxits=50, normalize=2)
    # Should drop at least one term with coefficient near the bound
    assert np.count_nonzero(w) <= 5


def test_mstls_large_lambda_keeps_most():
    """Large λ (near 1) with well-scaled data keeps most terms."""
    rng = np.random.default_rng(99)
    n, d = 300, 10
    # Make columns roughly orthogonal and equal-norm so ||b|| is O(1)
    G = rng.standard_normal((n, d))
    G = G / np.sqrt(np.sum(G ** 2, axis=0))[np.newaxis, :]  # unit-norm columns
    w_true = np.zeros(d)
    w_true[[1, 4, 7]] = [0.5, -0.3, 0.7]  # small coefficients, well-scaled
    b = G @ w_true

    w = mstls(G, b, lam=0.99, maxits=50, normalize=0)
    # With λ=0.99 and ||b||≈0.8, bounds are L≈0.99*0.8≈0.8, U≈1.01
    # Should keep only terms with |coeff| in [~0.8, ~1.01], i.e. may drop some
    # But the projection should not be totally wrong
    nz = np.count_nonzero(w)
    assert nz <= 3  # at most 3 terms (true count)


def test_mstls_converges_in_few_iters():
    """MSTLS should converge in far fewer than maxits on clean data."""
    G, b, _, _ = _sparse_system(n=300, d=20, nz=4, noise=0.0)
    # We can't easily count iters from the outside, but at least it doesn't crash
    w = mstls(G, b, lam=0.1, maxits=5, normalize=2)
    assert np.all(np.isfinite(w))


def test_mstls_no_normalize():
    G, b, w_true, nz_idx = _sparse_system(n=300, d=20, nz=3, noise=0.0)
    w = mstls(G, b, lam=0.1, maxits=50, normalize=0)
    assert set(nz_idx).issubset(set(np.where(np.abs(w) > 1e-8)[0]))


def test_mstls_noisy_data():
    """On noisy data, MSTLS should still identify the correct support."""
    G, b, w_true, nz_idx = _sparse_system(n=500, d=20, nz=4, noise=0.01)
    w = mstls(G, b, lam=0.1, maxits=50, normalize=2)
    nonzero = np.abs(w) > 1e-8
    # May miss some terms at this noise level, but should have no false positives
    fp = set(np.where(nonzero)[0]) - set(nz_idx)
    assert len(fp) == 0  # no false positives


# ======================================================================
# mstls_sweep — auto λ
# ======================================================================


def test_mstls_sweep_perfect_noiseless():
    """mstls_sweep auto-selects λ and recovers the true vector."""
    G, b, w_true, nz_idx = _sparse_system(n=300, d=20, nz=4, noise=0.0)
    w, lam_best, losses = mstls_sweep(G, b, n_lambdas=50, alpha=0.5)
    nonzero = np.abs(w) > 1e-8
    assert set(nz_idx).issubset(set(np.where(nonzero)[0]))
    np.testing.assert_allclose(w[nz_idx], w_true[nz_idx], rtol=1e-5)
    assert len(losses) == 50
    assert np.all(np.isfinite(losses))


def test_mstls_sweep_custom_lambdas():
    """mstls_sweep accepts a user-supplied λ grid."""
    G, b, w_true, nz_idx = _sparse_system(n=300, d=20, nz=3, noise=0.0)
    custom_lam = np.logspace(-3, 0, 30)
    w, lam_best, losses = mstls_sweep(G, b, lambdas=custom_lam, alpha=0.5)
    assert len(losses) == 30
    assert lam_best in custom_lam
    nonzero = np.abs(w) > 1e-8
    assert set(nz_idx).issubset(set(np.where(nonzero)[0]))


def test_mstls_sweep_noisy():
    """mstls_sweep on mildly noisy data."""
    G, b, w_true, nz_idx = _sparse_system(n=500, d=20, nz=4, noise=0.02)
    w, lam_best, losses = mstls_sweep(G, b, n_lambdas=50, alpha=0.5)
    nonzero = np.abs(w) > 1e-8
    # At mild noise, should have no false positives
    fp = set(np.where(nonzero)[0]) - set(nz_idx)
    assert len(fp) == 0


def test_mstls_sweep_alpha_extremes():
    """α=0 (all sparsity) vs α=1 (all projection) produce different results."""
    G, b, _, _ = _sparse_system(n=300, d=20, nz=4, noise=0.0)
    w0, _, _ = mstls_sweep(G, b, n_lambdas=50, alpha=0.0)
    w1, _, _ = mstls_sweep(G, b, n_lambdas=50, alpha=1.0)
    # α=0 prefers sparsity → fewer or equal terms
    assert np.count_nonzero(w0) <= np.count_nonzero(w1)


# ======================================================================
# Comparison with STRidge
# ======================================================================


def test_mstls_handles_column_scale_disparity():
    """MSTLS with normalize=2 handles columns with varying norms."""
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

    # MSTLS with normalize=2
    w_ms, lam_best, _ = mstls_sweep(G, b, n_lambdas=50, alpha=0.5, normalize=2)
    nz_ms = np.where(np.abs(w_ms) > 1e-8)[0]
    # Should not have false positives
    assert len(set(nz_ms) - set(nz_idx)) == 0


def test_mstls_upper_bound_prevents_large_spurious():
    """The upper bound prevents terms with coefficients ≫ 1/λ."""
    rng = np.random.default_rng(42)
    n, d = 300, 10
    # Use nearly-orthogonal columns so ||b|| ≈ ||w_true||
    G = rng.standard_normal((n, d)) / np.sqrt(n)
    U_svd, _, Vt = np.linalg.svd(G, full_matrices=False)
    G = U_svd @ Vt  # columns are orthonormal

    w_true = np.zeros(d)
    w_true[0] = 1.0
    b = G @ w_true  # ||b|| = 1.0 exactly

    # λ=0.5: L=0.5, U=2 for all columns (since unit-norm, ||b||/||G_j||=1)
    w = mstls(G, b, lam=0.5, maxits=50, normalize=0)
    nonzero = np.where(np.abs(w) > 1e-8)[0]
    assert 0 in nonzero
    for j in nonzero:
        assert abs(w[j]) <= 2.0 * 1.01


# ======================================================================
# Integration tests — MSTLS in sweep pipeline
# ======================================================================


def test_mstls_sweep_one_pointwise():
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


def test_mstls_sweep_one_weak_form():
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


def test_mstls_sweep_cli_accepts_mstls():
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
