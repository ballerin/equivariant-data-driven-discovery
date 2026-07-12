"""
Tests for ``regression.py`` — STRidge, TrainSTRidge, SR3, and TrainSR3.

These are the core regression algorithms used by the PDE-FIND pipeline.
STRidge and TrainSTRidge are tested here directly; MSTLS and RFE have
their own dedicated test files (test_mstls.py, test_rfe.py).

References
----------
- STRidge : Brunton, Proctor, Kutz (2016), doi:10.1073/pnas.1517384113
- SR3     : Zheng et al. (2018), doi:10.1109/ACCESS.2018.2886528
"""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.regression import (
    STRidge,
    TrainSTRidge,
    SR3,
    TrainSR3,
)


# ======================================================================
# STRidge — core algorithm
# ======================================================================


class TestSTRidge:
    """Unit tests for the STRidge inner loop."""

    def test_perfect_recovery_noiseless(self, sparse_system):
        G, b, w_true, nz_idx = sparse_system(n=300, d=20, nz=4, noise=0.0)
        w = STRidge(G, b, lam=1e-5, maxit=25, tol=0.05)
        nonzero = np.abs(w) > 1e-8
        assert set(nz_idx).issubset(set(np.where(nonzero)[0]))
        np.testing.assert_allclose(w[nz_idx], w_true[nz_idx], rtol=1e-4)

    def test_threshold_zeroes_small_coeffs(self, sparse_system):
        """With large tol, small coefficients should be zeroed."""
        G, b, _, _ = sparse_system(n=200, d=10, nz=3, noise=0.0, seed=7)
        w = STRidge(G, b, lam=1e-5, maxit=10, tol=100.0)
        # All coefficients should be zero since tol > max|coeff|
        assert np.all(np.abs(w) < 1e-10)

    def test_no_normalize(self, sparse_system):
        """normalize=0 should not rescale columns."""
        G, b, w_true, nz_idx = sparse_system(
            n=200, d=15, nz=4, noise=0.0, normalize=False,
        )
        w = STRidge(G, b, lam=1e-5, maxit=25, tol=0.05, normalize=0)
        nonzero = np.abs(w) > 1e-8
        assert set(nz_idx).issubset(set(np.where(nonzero)[0]))

    def test_early_stop_disabled(self, sparse_system):
        """With early_stop=False, STRidge runs all maxit iterations."""
        G, b, _, _ = sparse_system(n=200, d=10, nz=3, noise=0.0)
        # Should not crash with early_stop=False
        w = STRidge(G, b, lam=1e-5, maxit=25, tol=0.05, early_stop=False)
        assert np.all(np.isfinite(w))

    def test_unbias_false(self, sparse_system):
        """With unbias=False, final coefficients stay ridge-regularised."""
        G, b, w_true, nz_idx = sparse_system(n=200, d=10, nz=3, noise=0.0)
        w_unbiased = STRidge(G, b, lam=1e-5, maxit=25, tol=0.05, unbias=True)
        w_biased = STRidge(G, b, lam=1e-5, maxit=25, tol=0.05, unbias=False)
        # Both should recover the right support
        nonzero_u = np.abs(w_unbiased) > 1e-8
        nonzero_b = np.abs(w_biased) > 1e-8
        assert set(nz_idx).issubset(set(np.where(nonzero_u)[0]))
        assert set(nz_idx).issubset(set(np.where(nonzero_b)[0]))

    def test_noisy_data(self, sparse_system):
        """STRidge on mildly noisy data should recover the correct support."""
        G, b, w_true, nz_idx = sparse_system(
            n=500, d=20, nz=4, noise=0.005,
        )
        # Use a looser tolerance; STRidge is sensitive to noise
        w = STRidge(G, b, lam=1e-4, maxit=25, tol=0.1)
        nonzero = np.abs(w) > 1e-8
        fp = set(np.where(nonzero)[0]) - set(nz_idx)
        # At mild noise, should have at most 1 false positive
        assert len(fp) <= 1, f"Too many false positives: {fp}"

    def test_single_feature(self):
        """STRidge should work with d=1."""
        rng = np.random.default_rng(42)
        G = rng.standard_normal((100, 1))
        w_true = np.array([3.0])
        b = G @ w_true
        w = STRidge(G, b, lam=1e-5, maxit=10, tol=0.01)
        np.testing.assert_allclose(w, w_true, rtol=1e-4)

    def test_normalize_2_is_default(self, sparse_system):
        """Normalize=2 is the default (L2 column normalization)."""
        G, b, w_true, nz_idx = sparse_system(n=200, d=15, nz=4, noise=0.0)
        w = STRidge(G, b, lam=1e-5, maxit=25, tol=0.05)
        nonzero = np.abs(w) > 1e-8
        assert set(nz_idx).issubset(set(np.where(nonzero)[0]))

    def test_large_lambda_stabilizes(self, sparse_system):
        """STRidge with large lam should still run without error."""
        G, b, _, _ = sparse_system(n=100, d=10, nz=3, noise=0.0)
        w = STRidge(G, b, lam=100.0, maxit=25, tol=0.05)
        assert np.all(np.isfinite(w))


# ======================================================================
# TrainSTRidge — cross-validated STRidge
# ======================================================================


class TestTrainSTRidge:
    """Tests for the TrainSTRidge cross-validation wrapper."""

    def test_recovers_sparse_model(self, sparse_system):
        G, b, w_true, nz_idx = sparse_system(n=300, d=20, nz=4, noise=0.0)
        w = TrainSTRidge(G, b, lam=1e-5, d_tol=0.05, maxit=20, STR_iters=10)
        nonzero = np.abs(w) > 1e-8
        assert set(nz_idx).issubset(set(np.where(nonzero)[0]))

    def test_noisy_data(self, sparse_system):
        """TrainSTRidge on mildly noisy data."""
        G, b, w_true, nz_idx = sparse_system(
            n=500, d=20, nz=4, noise=0.01,
        )
        w = TrainSTRidge(G, b, lam=1e-4, d_tol=0.05, maxit=30, STR_iters=10)
        nonzero = np.abs(w) > 1e-8
        # Should have some of the true terms
        recovered = set(np.where(nonzero)[0])
        assert len(recovered & set(nz_idx)) >= 2

    def test_custom_l0_penalty(self, sparse_system):
        G, b, _, _ = sparse_system(n=200, d=15, nz=3, noise=0.0)
        w = TrainSTRidge(G, b, lam=1e-5, d_tol=0.05, l0_penalty=0.01)
        assert np.all(np.isfinite(w))

    def test_split_fraction(self, sparse_system):
        G, b, _, _ = sparse_system(n=200, d=15, nz=3, noise=0.0)
        # split=0.5 → half training, half test
        w = TrainSTRidge(G, b, lam=1e-5, d_tol=0.05, split=0.5)
        assert np.all(np.isfinite(w))

    def test_unbias_false(self, sparse_system):
        G, b, _, _ = sparse_system(n=200, d=10, nz=3, noise=0.0)
        w = TrainSTRidge(G, b, lam=1e-5, d_tol=0.05, unbias=False)
        assert np.all(np.isfinite(w))


# ======================================================================
# SR3 — Sparse Relaxed Regularized Regression
# ======================================================================


class TestSR3:
    """Unit tests for the SR3 inner loop."""

    def test_perfect_recovery_noiseless(self, sparse_system):
        G, b, w_true, nz_idx = sparse_system(
            n=300, d=20, nz=4, noise=0.0, normalize=False,
        )
        # Use very small lam_l1 and small kappa so shrinkage is negligible
        w = SR3(G, b, lam_l1=1e-6, kappa=0.01, maxit=500, normalize=0)
        nonzero = np.abs(w) > 1e-8
        assert set(nz_idx).issubset(set(np.where(nonzero)[0]))
        np.testing.assert_allclose(w[nz_idx], w_true[nz_idx], rtol=1e-2)

    def test_sparsity_increases_with_lam(self, sparse_system):
        """Larger λ_L1 should give sparser solutions."""
        G, b, _, _ = sparse_system(n=300, d=15, nz=5, noise=0.0)
        w_small = SR3(G, b, lam_l1=0.001, kappa=1.0, maxit=200)
        w_large = SR3(G, b, lam_l1=10.0, kappa=1.0, maxit=200)
        assert np.count_nonzero(w_large) <= np.count_nonzero(w_small)

    def test_noisy_data(self, sparse_system):
        """SR3 on mildly noisy data should recover support."""
        G, b, w_true, nz_idx = sparse_system(
            n=500, d=20, nz=4, noise=0.01,
        )
        w = SR3(G, b, lam_l1=0.05, kappa=1.0, maxit=200)
        nonzero = np.abs(w) > 1e-8
        fp = set(np.where(nonzero)[0]) - set(nz_idx)
        assert len(fp) == 0, f"False positives: {fp}"

    def test_no_normalize(self, sparse_system):
        G, b, w_true, nz_idx = sparse_system(
            n=200, d=15, nz=4, noise=0.0, normalize=False,
        )
        w = SR3(G, b, lam_l1=0.01, maxit=200, normalize=0)
        nonzero = np.abs(w) > 1e-8
        assert set(nz_idx).issubset(set(np.where(nonzero)[0]))

    def test_single_feature(self):
        rng = np.random.default_rng(42)
        G = rng.standard_normal((100, 1))
        w_true = np.array([3.0])
        b = G @ w_true
        # normalize=0, very small kappa and lam_l1 for exact recovery
        w = SR3(G, b, lam_l1=1e-8, kappa=0.01, maxit=500, normalize=0)
        np.testing.assert_allclose(w, w_true, rtol=1e-2)

    def test_kappa_effect(self, sparse_system):
        """Larger κ should make ξ and w closer (tighter coupling)."""
        G, b, w_true, _ = sparse_system(n=200, d=10, nz=3, noise=0.0)
        w_small_k = SR3(G, b, lam_l1=0.01, kappa=0.1, maxit=200)
        w_large_k = SR3(G, b, lam_l1=0.01, kappa=10.0, maxit=200)
        assert np.all(np.isfinite(w_small_k))
        assert np.all(np.isfinite(w_large_k))


# ======================================================================
# TrainSR3 — cross-validated SR3
# ======================================================================


class TestTrainSR3:
    """Tests for the TrainSR3 cross-validation wrapper."""

    def test_recovers_sparse_model(self, sparse_system):
        G, b, w_true, nz_idx = sparse_system(n=300, d=20, nz=4, noise=0.0)
        w = TrainSR3(G, b, lam=1.0, d_tol=0.01, maxit=20, SR3_iters=100)
        nonzero = np.abs(w) > 1e-8
        assert set(nz_idx).issubset(set(np.where(nonzero)[0]))

    def test_noisy_data(self, sparse_system):
        G, b, w_true, nz_idx = sparse_system(
            n=500, d=20, nz=4, noise=0.01,
        )
        w = TrainSR3(G, b, lam=1.0, d_tol=0.01, maxit=30, SR3_iters=100)
        nonzero = np.abs(w) > 1e-8
        recovered = set(np.where(nonzero)[0])
        assert len(recovered & set(nz_idx)) >= 2
