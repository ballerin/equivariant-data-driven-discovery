"""
Tests for ``centering`` in ``regression.py`` — intercept partialling for
sparse PDE regression.

``centering`` drops the constant column, centres the target and the
remaining columns, delegates to a sparse solver, then recovers the
intercept in closed form and thresholds it.  The returned vector is
full-length and aligned with the library's columns, so it is a drop-in
replacement for the output of ``TrainSTRidge`` / ``TrainRFE``.

References
----------
- Frisch & Waugh (1933), doi:10.2307/1907330
- Lovell (1963), doi:10.1080/01621459.1963.10480682
"""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.regression import (  # noqa: E402
    centering,
    TrainSTRidge,
    TrainRFE,
    _resolve_const_index,
)


# ======================================================================
# Local fixtures
#
# The shared ``sparse_system`` fixture draws standard-normal columns, which
# are already mean-zero — centring would be a no-op on it.  These systems
# deliberately carry nonzero column means plus an explicit constant column.
# ======================================================================


def offset_library(n=400, d=12, nz=3, w0=0.0, noise=0.0, seed=0,
                   const_at=0, offset_scale=3.0):
    """Build (Theta, ut, w_true, descs) with a constant column and offset columns.

    ``Theta[:, const_at]`` is all ones; every other column has a nonzero mean
    drawn on the scale of ``offset_scale``, so centring genuinely changes the
    design matrix.
    """
    rng = np.random.default_rng(seed)

    n_free = d - 1
    X00 = rng.standard_normal((n, n_free))
    X00 += offset_scale * rng.standard_normal(n_free)[np.newaxis, :]

    w_free = np.zeros(n_free)
    nz_idx = rng.choice(n_free, size=nz, replace=False)
    w_free[nz_idx] = rng.uniform(0.5, 2.0, size=nz) * rng.choice([-1, 1], size=nz)

    # Assemble with the constant column at position ``const_at``
    cols = list(X00.T)
    cols.insert(const_at, np.ones(n))
    Theta = np.column_stack(cols)

    w_true = np.insert(w_free, const_at, w0)
    descs = [f"f{i}" for i in range(n_free)]
    descs.insert(const_at, "1")

    ut = Theta @ w_true
    if noise > 0:
        ut = ut + noise * np.std(ut) * rng.standard_normal(n)

    return Theta, ut, w_true, descs


def ols(X, y):
    """Plain unregularised least squares, in (design, target) order."""
    return np.linalg.lstsq(X, y, rcond=None)[0]


# ======================================================================
# Frisch-Waugh-Lovell exactness
# ======================================================================


class TestFWLIdentity:
    """For OLS, centring must reproduce the full-library solution exactly."""

    def test_matches_full_ols(self):
        Theta, ut, _, descs = offset_library(n=400, d=12, nz=4, w0=1.7,
                                             noise=0.05, seed=1)
        w_direct = np.linalg.lstsq(Theta, ut, rcond=None)[0]
        w_centered = centering(Theta, ut, ols, tol=0.0, descs=descs)
        np.testing.assert_allclose(w_centered, w_direct, atol=1e-9, rtol=1e-7)

    def test_matches_full_ols_zero_intercept(self):
        Theta, ut, _, descs = offset_library(n=300, d=10, nz=3, w0=0.0,
                                             noise=0.0, seed=2)
        w_direct = np.linalg.lstsq(Theta, ut, rcond=None)[0]
        w_centered = centering(Theta, ut, ols, tol=0.0, descs=descs)
        np.testing.assert_allclose(w_centered, w_direct, atol=1e-9, rtol=1e-7)

    def test_intercept_returned_matches_slot(self):
        Theta, ut, _, descs = offset_library(n=300, d=10, nz=3, w0=2.5, seed=3)
        w, w0 = centering(Theta, ut, ols, tol=0.0, descs=descs,
                          return_intercept=True)
        assert w0 == pytest.approx(w[descs.index("1")])


# ======================================================================
# Recovery with the real sparse solvers
# ======================================================================


class TestPlantedIntercept:
    """A planted nonzero intercept must be recovered with its support."""

    def test_stridge_recovers_intercept_and_support(self):
        Theta, ut, w_true, descs = offset_library(n=500, d=12, nz=3, w0=1.25,
                                                  noise=0.0, seed=4)
        w = centering(Theta, ut, lambda X, y: TrainSTRidge(X, y, 1e-6, 1.0),
                      tol=0.0, descs=descs)
        np.testing.assert_allclose(w, w_true, atol=1e-6)

    def test_rfe_ridge_recovers_intercept_and_support(self):
        Theta, ut, w_true, descs = offset_library(n=500, d=12, nz=3, w0=1.25,
                                                  noise=0.0, seed=5)
        w = centering(Theta, ut, lambda X, y: TrainRFE(X, y, lam=1e-6),
                      tol=0.0, descs=descs)
        support_true = np.abs(w_true) > 1e-12
        support_rec = np.abs(w) > 1e-8
        np.testing.assert_array_equal(support_rec, support_true)
        # RFE applies ridge (lam) without an unbiasing step, so the surviving
        # coefficients carry an O(lam) shrinkage bias.
        np.testing.assert_allclose(w, w_true, atol=1e-4)

    def test_noisy_recovery_keeps_support(self):
        Theta, ut, w_true, descs = offset_library(n=800, d=12, nz=3, w0=1.5,
                                                  noise=0.01, seed=6)
        w = centering(Theta, ut, lambda X, y: TrainSTRidge(X, y, 1e-6, 1.0),
                      tol=0.0, descs=descs)
        support_true = np.abs(w_true) > 1e-12
        support_rec = np.abs(w) > 1e-8
        np.testing.assert_array_equal(support_rec, support_true)


# ======================================================================
# The intercept threshold
# ======================================================================


class TestIntercepThreshold:
    """``|w0| < tol * mean(|ut|)`` must zero the intercept, and only then."""

    def test_threshold_fires_on_spurious_intercept(self):
        # True intercept is zero; noise leaves a small spurious one.
        Theta, ut, w_true, descs = offset_library(n=600, d=12, nz=3, w0=0.0,
                                                  noise=0.02, seed=7)
        j = descs.index("1")
        alg = lambda X, y: TrainSTRidge(X, y, 1e-6, 1.0)  # noqa: E731

        w_raw = centering(Theta, ut, alg, tol=0.0, descs=descs)
        assert w_raw[j] != 0.0  # a spurious intercept is present to begin with

        w = centering(Theta, ut, alg, tol=1e-2, descs=descs)
        assert w[j] == 0.0

        # Without the threshold the spurious intercept is the sole false
        # positive; with it, the recovery is exact.
        active_true = np.abs(w_true) > 1e-12
        assert np.sum((np.abs(w_raw) > 1e-10) & ~active_true) == 1
        assert np.sum((np.abs(w) > 1e-10) & ~active_true) == 0

    def test_threshold_does_not_fire_on_genuine_intercept(self):
        Theta, ut, _, descs = offset_library(n=600, d=12, nz=3, w0=2.0,
                                             noise=0.01, seed=8)
        j = descs.index("1")
        w = centering(Theta, ut, ols, tol=0.1, descs=descs)
        assert w[j] == pytest.approx(2.0, abs=1e-2)

    def test_tol_zero_disables_threshold(self):
        Theta, ut, _, descs = offset_library(n=400, d=10, nz=3, w0=0.0,
                                             noise=0.02, seed=9)
        j = descs.index("1")
        w = centering(Theta, ut, ols, tol=0.0, descs=descs)
        assert w[j] != 0.0

    def test_zero_target_does_not_divide_by_zero(self):
        Theta, ut, _, descs = offset_library(n=200, d=8, nz=2, w0=0.0, seed=10)
        ut = np.zeros_like(ut)
        w = centering(Theta, ut, ols, tol=1e-3, descs=descs)
        assert np.all(np.isfinite(w))


# ======================================================================
# Locating the constant column
# ======================================================================


class TestConstantColumnResolution:
    """The constant must be found by name, not assumed to be at index 0."""

    @pytest.mark.parametrize("const_at", [0, 1, 5, 11])
    def test_invariant_to_constant_position(self, const_at):
        Theta, ut, w_true, descs = offset_library(n=400, d=12, nz=3, w0=1.4,
                                                  noise=0.0, seed=11,
                                                  const_at=const_at)
        w = centering(Theta, ut, lambda X, y: TrainSTRidge(X, y, 1e-6, 1.0),
                      tol=0.0, descs=descs)
        assert w[const_at] == pytest.approx(1.4, abs=1e-6)
        np.testing.assert_allclose(w, w_true, atol=1e-6)

    def test_explicit_const_index(self):
        Theta, ut, w_true, _ = offset_library(n=300, d=10, nz=3, w0=1.0,
                                              noise=0.0, seed=12, const_at=4)
        w = centering(Theta, ut, ols, tol=0.0, const_index=4)
        np.testing.assert_allclose(w, w_true, atol=1e-8)

    def test_numeric_fallback_without_descs(self):
        Theta, ut, w_true, _ = offset_library(n=300, d=10, nz=3, w0=1.0,
                                              noise=0.0, seed=13, const_at=3)
        w = centering(Theta, ut, ols, tol=0.0)
        np.testing.assert_allclose(w, w_true, atol=1e-8)

    def test_resolver_returns_none_when_absent(self):
        rng = np.random.default_rng(14)
        Theta = rng.standard_normal((100, 5)) + 2.0
        assert _resolve_const_index(Theta) is None
        assert _resolve_const_index(Theta, descs=list("abcde")) is None

    def test_descs_length_mismatch_raises(self):
        Theta, ut, _, descs = offset_library(n=100, d=8, nz=2, seed=15)
        with pytest.raises(ValueError, match="length"):
            centering(Theta, ut, ols, descs=descs[:-1])


# ======================================================================
# Libraries with no constant column
# ======================================================================


class TestMissingConstantColumn:
    """Equivariant libraries may carry no constant; fail loudly by default."""

    def _no_const(self, seed=16):
        rng = np.random.default_rng(seed)
        X = rng.standard_normal((300, 6)) + 2.0
        w_true = np.array([0.0, 1.5, 0.0, -0.8, 0.0, 0.0])
        y = X @ w_true + 0.7  # a genuine offset the library cannot represent
        return X, y, w_true

    def test_raises_by_default(self):
        X, y, _ = self._no_const()
        with pytest.raises(ValueError, match="No constant column"):
            centering(X, y, ols)

    def test_allow_missing_returns_intercept_separately(self):
        X, y, w_true = self._no_const()
        w, w0 = centering(X, y, ols, tol=0.0, allow_missing=True,
                          return_intercept=True)
        assert len(w) == X.shape[1]
        np.testing.assert_allclose(w, w_true, atol=1e-8)
        assert w0 == pytest.approx(0.7, abs=1e-8)


# ======================================================================
# Degeneracy created by centring
# ======================================================================


class TestDegenerateColumns:
    """A second constant column centres to zero and must not absorb noise."""

    def test_duplicate_constant_is_zeroed(self):
        Theta, ut, w_true, descs = offset_library(n=400, d=10, nz=3, w0=1.0,
                                                  noise=0.01, seed=17)
        # Append a second, differently-scaled constant column
        Theta = np.column_stack([Theta, 5.0 * np.ones(Theta.shape[0])])
        descs = descs + ["const_dup"]
        w_true = np.append(w_true, 0.0)

        w = centering(Theta, ut, lambda X, y: TrainSTRidge(X, y, 1e-6, 1.0),
                      tol=1e-2, descs=descs)
        assert w[-1] == 0.0
        assert np.all(np.isfinite(w))
        support_true = np.abs(w_true) > 1e-12
        support_rec = np.abs(w) > 1e-6
        np.testing.assert_array_equal(support_rec, support_true)

    def test_all_zero_column_is_zeroed(self):
        Theta, ut, _, descs = offset_library(n=300, d=10, nz=3, w0=1.0,
                                             noise=0.0, seed=18)
        Theta = np.column_stack([Theta, np.zeros(Theta.shape[0])])
        descs = descs + ["dead"]
        w = centering(Theta, ut, ols, tol=0.0, descs=descs)
        assert w[-1] == 0.0
        assert np.all(np.isfinite(w))


# ======================================================================
# Interface contract
# ======================================================================


class TestInterface:
    """Shape, alignment, and solver-return handling."""

    def test_returns_full_length_vector(self):
        Theta, ut, _, descs = offset_library(n=200, d=9, nz=2, seed=19)
        w = centering(Theta, ut, ols, descs=descs)
        assert w.shape == (Theta.shape[1],)

    def test_unwraps_tuple_returning_solver(self):
        Theta, ut, w_true, descs = offset_library(n=300, d=10, nz=3, w0=1.0,
                                                  noise=0.0, seed=20)

        def tuple_solver(X, y):
            return ols(X, y), "extra", 42

        w = centering(Theta, ut, tuple_solver, tol=0.0, descs=descs)
        np.testing.assert_allclose(w, w_true, atol=1e-8)

    def test_row_mismatch_raises(self):
        Theta, ut, _, descs = offset_library(n=200, d=8, nz=2, seed=21)
        with pytest.raises(ValueError, match="rows"):
            centering(Theta, ut[:-1], ols, descs=descs)

    def test_does_not_mutate_inputs(self):
        Theta, ut, _, descs = offset_library(n=200, d=8, nz=2, w0=1.0, seed=22)
        Theta_copy, ut_copy = Theta.copy(), ut.copy()
        centering(Theta, ut, ols, tol=0.1, descs=descs)
        np.testing.assert_array_equal(Theta, Theta_copy)
        np.testing.assert_array_equal(ut, ut_copy)

    def test_e_data_formula_still_works(self):
        # The full-length return must plug straight into ||Theta w - ut||/||ut||
        Theta, ut, _, descs = offset_library(n=400, d=12, nz=3, w0=1.5,
                                             noise=0.0, seed=23)
        w = centering(Theta, ut, lambda X, y: TrainSTRidge(X, y, 1e-6, 1.0),
                      tol=0.0, descs=descs)
        E_data = np.linalg.norm(Theta @ w - ut) / max(np.linalg.norm(ut), 1e-14)
        assert E_data < 1e-8

    @pytest.mark.parametrize("c_val", [7.44, 0.2, -3.0])
    def test_non_unit_constant_column(self, c_val):
        # Weak form: the constant column is the test-function mass <psi, 1>,
        # not 1.  The slot must hold the PDE's constant, so Theta @ w is the fit.
        Theta, ut, w_true, descs = offset_library(n=400, d=10, nz=3, w0=1.25,
                                                  noise=0.0, seed=24)
        j = descs.index("1")
        Theta = Theta.copy()
        Theta[:, j] = c_val
        ut = Theta @ w_true
        w, w0 = centering(Theta, ut, ols, descs=descs, return_intercept=True)
        np.testing.assert_allclose(w, w_true, atol=1e-9)
        assert w0 == pytest.approx(1.25 * c_val)
        np.testing.assert_allclose(Theta @ w, ut, atol=1e-8)
