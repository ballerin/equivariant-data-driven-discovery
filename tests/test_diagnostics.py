"""
Tests for `diagnostics.py` — pairwise_max_correlation, vif_factors,
condition_numbers, collinearity_report.
"""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.diagnostics import (                              # noqa: E402
    collinearity_report,
    condition_numbers,
    pairwise_max_correlation,
    vif_factors,
)


# ---------------------------------------------------------------------------
# pairwise_max_correlation
# ---------------------------------------------------------------------------

def test_pairwise_corr_empty():
    """Single column or zero columns — max corr is 0."""
    assert pairwise_max_correlation(np.ones((100, 1))) == 0.0


def test_pairwise_corr_orthogonal():
    """Orthogonal columns — max |corr| ≈ 0."""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(1000, 4))
    A = A - A.mean(axis=0, keepdims=True)
    # Orthogonalise via QR
    Q, _ = np.linalg.qr(A)
    max_corr = pairwise_max_correlation(Q)
    assert max_corr < 1e-12, f"Expected ~0, got {max_corr}"


def test_pairwise_corr_perfect_collinear():
    """col_2 = 2 * col_1 → max |corr| = 1."""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(100, 2))
    A[:, 1] = 2.0 * A[:, 0]  # perfect collinearity
    max_corr = pairwise_max_correlation(A)
    assert np.isclose(max_corr, 1.0, rtol=1e-12)


def test_pairwise_corr_high_correlation():
    """col_2 = col_1 + small noise → max |corr| ≈ 1."""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(1000, 2))
    A[:, 1] = A[:, 0] + 0.001 * rng.normal(size=1000)
    max_corr = pairwise_max_correlation(A)
    assert max_corr > 0.999


def test_pairwise_corr_identity():
    """Identity matrix → max |corr| = 0 (columns are orthogonal)."""
    I = np.eye(100)
    max_corr = pairwise_max_correlation(I)
    # With 100 rows and n-1 normalisation, correlations are exactly -1/(n-1)
    # for any pair.  With 100 rows this is ≈ -0.0101.
    assert max_corr < 0.02


# ---------------------------------------------------------------------------
# vif_factors
# ---------------------------------------------------------------------------

def test_vif_orthogonal():
    """Orthogonal columns → all VIF ≈ 1."""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(200, 4))
    Q, _ = np.linalg.qr(A)
    vif, skipped = vif_factors(Q)
    assert not skipped.any()
    assert np.allclose(vif, 1.0, atol=0.05), f"VIFs: {vif}"


def test_vif_perfect_collinear_pair():
    """col_2 = 2 * col_1 → VIF_1 and VIF_2 both capped (R² → 1)."""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(200, 3))
    A[:, 0] = rng.normal(size=200)  # independent
    A[:, 1] = rng.normal(size=200)  # independent
    A[:, 2] = 3.0 * A[:, 0]         # collinear with col_0
    vif, skipped = vif_factors(A)
    assert not skipped.any()
    # col_0 and col_2 should have huge VIF
    assert vif[0] > 1e5, f"VIF_0 = {vif[0]}, expected very large"
    assert vif[2] > 1e5, f"VIF_2 = {vif[2]}, expected very large"
    # col_1 is independent → VIF ≈ 1
    assert vif[1] < 2.0, f"VIF_1 = {vif[1]}, expected ≈ 1"


def test_vif_constant_column_skipped():
    """Constant column → skipped (VIF = NaN, std ≈ 0)."""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(100, 3))
    A[:, 0] = 5.0  # constant column
    vif, skipped = vif_factors(A)
    assert skipped[0], "Constant column should be skipped"
    assert np.isnan(vif[0]), "VIF for constant column should be NaN"
    assert not skipped[1], "Non-constant column should not be skipped"
    assert not skipped[2], "Non-constant column should not be skipped"


def test_vif_three_column_dependence():
    """col_3 = col_1 + col_2 + small noise → VIF_3 large, VIF_1,2 moderate."""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(500, 4))
    A[:, 3] = A[:, 0] + A[:, 1] + 0.001 * rng.normal(size=500)
    vif, skipped = vif_factors(A)
    assert not skipped.any()
    # col_3 (dependent) should have very high VIF
    assert vif[3] > 1e5, f"VIF_3 = {vif[3]}, expected very large"
    # col_0 and col_1 should have moderate VIF (they're partly explained by others)
    assert vif[0] > 2.0, f"VIF_0 = {vif[0]}, expected > 2"


def test_vif_identity():
    """Identity matrix → all VIF ≈ n/(n-2) (well-conditioned)."""
    I = np.eye(50)
    vif, skipped = vif_factors(I)
    assert not skipped.any()
    # For identity, each column is nearly orthogonal to others
    # VIF should be modest
    assert np.all(vif < 2.0), f"VIFs for identity: {vif}"


def test_vif_single_column():
    """Single column → all NaN (nothing to regress on)."""
    A = np.random.default_rng(42).normal(size=(100, 1))
    vif, skipped = vif_factors(A)
    assert np.isnan(vif[0])
    # A single non-constant column has variance but nothing to regress
    # on (n_cols < 2).  Not "skipped" (that flag means constant).


# ---------------------------------------------------------------------------
# condition_numbers
# ---------------------------------------------------------------------------

def test_cond_identity():
    """Identity matrix → κ = 1."""
    I = np.eye(50)
    cond_raw, cond_norm = condition_numbers(I)
    assert np.isclose(cond_raw, 1.0, rtol=1e-12)
    assert np.isclose(cond_norm, 1.0, rtol=1e-12)


def test_cond_nearly_singular():
    """Nearly rank-deficient → κ very large."""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(200, 3))
    A[:, 2] = A[:, 0] + A[:, 1] + 1e-12 * rng.normal(size=200)
    cond_raw, cond_norm = condition_numbers(A)
    assert cond_raw > 1e8, f"cond_raw = {cond_raw}, expected > 1e8"
    assert cond_norm > 1e8, f"cond_norm = {cond_norm}, expected > 1e8"


def test_cond_normalisation_effect():
    """Column scaling changes raw κ but not normalised κ."""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(100, 3))
    # Scale columns dramatically
    A_scaled = A * np.array([1.0, 1e6, 1e-3])
    cond_raw_A, cond_norm_A = condition_numbers(A)
    cond_raw_S, cond_norm_S = condition_numbers(A_scaled)
    # Raw κ should change due to scaling
    assert not np.isclose(cond_raw_A, cond_raw_S, rtol=0.1)
    # Normalised κ should be (nearly) unchanged
    assert np.isclose(cond_norm_A, cond_norm_S, rtol=1e-10)


# ---------------------------------------------------------------------------
# collinearity_report
# ---------------------------------------------------------------------------

def test_collinearity_report_keys():
    """Report dict has all expected keys."""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(100, 4))
    desc = ["a", "b", "c", "d"]
    report = collinearity_report(A, desc)
    expected_keys = {
        "n_rows", "n_cols", "max_abs_corr", "max_vif", "mean_vif",
        "cond_raw", "cond_normalised", "skipped_columns", "vif_per_column",
    }
    assert set(report.keys()) == expected_keys
    assert report["n_rows"] == 100
    assert report["n_cols"] == 4


def test_collinearity_report_constant_column():
    """Constant column is reported in skipped_columns."""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(100, 3))
    A[:, 0] = 42.0
    desc = ["const", "x", "y"]
    report = collinearity_report(A, desc)
    assert report["skipped_columns"] == ["const"]
    # vif_per_column should have NaN for the constant
    assert np.isnan(report["vif_per_column"][0][1])
    assert not np.isnan(report["vif_per_column"][1][1])
