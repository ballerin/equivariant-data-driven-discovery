"""
Tests for ``evaluation.py`` — evaluate_recovery, evaluate_from_theta,
build_eval_theta, and best_of_sweep.

Also covers a regression bug: ground-truth coefficient vectors must always
be built from **equivariant** library descriptions, never from the standard
(coordinate-basis) library, which uses different operator names (f_xx vs Δf).
"""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.evaluation import (                                        # noqa: E402
    evaluate_recovery, evaluate_from_theta, build_eval_theta, best_of_sweep,
    evaluate_standard,
)
from equivariant_pde_find.library import (                                           # noqa: E402
    build_equivariant_library, build_standard_library,
)
from equivariant_pde_find.simulations import make_regression_data                    # noqa: E402
from equivariant_pde_find.simulations import allen_cahn                              # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _true_allen_cahn(desc_eq, eps=0.15):
    """Ground-truth Allen-Cahn coefficients (equivariant basis)."""
    w = np.zeros(len(desc_eq))
    w[desc_eq.index("f")] = 1.0
    w[desc_eq.index("f\u00b3")] = -1.0   # f³
    w[desc_eq.index("\u0394f")] = eps ** 2   # Δf
    return w


@pytest.fixture(scope="module")
def ac_data():
    """Simulate a small Allen-Cahn trajectory for testing."""
    snaps, times = allen_cahn(
        N=32, L=2 * np.pi, T_end=2.0, eps=0.15, dt=0.01,
        save_every=20, seed=42,
    )
    # Get descriptions from one snapshot
    _, desc = build_equivariant_library(snaps[0], 2 * np.pi, 2 * np.pi)
    true_eq = _true_allen_cahn(desc)
    return snaps, times, true_eq, desc


# ---------------------------------------------------------------------------
# Ground-truth / library description mismatch (regression guard)
# ---------------------------------------------------------------------------


class TestGroundTruthDescriptionMismatch:
    """Guard: ground truth must use equivariant descriptions, not standard."""

    def test_equivariant_descriptions_contain_delta(self, ac_data):
        """Equivariant library must contain 'Δf', 'Δ²f', etc."""
        snaps, times, _, _ = ac_data
        _, desc = build_equivariant_library(snaps[0], 2 * np.pi, 2 * np.pi)
        assert "Δf" in desc, f"Expected 'Δf' in {desc[:5]}..."
        assert "f³" in desc or "f\u00b3" in desc
        assert "Δ²f" in desc

    def test_standard_descriptions_do_not_contain_delta(self, ac_data):
        """Standard library must NOT contain 'Δf' (it has 'f_xx', etc.)."""
        snaps, times, _, _ = ac_data
        _, desc = build_standard_library(snaps[0], 2 * np.pi, 2 * np.pi,
                                          max_order=4, max_poly=3)
        assert "Δf" not in desc
        assert "f_xx" in desc
        assert "f_xxxx" in desc

    def test_true_fn_works_with_equivariant(self, ac_data):
        """_true_allen_cahn must succeed with equivariant descriptions."""
        snaps, times, _, _ = ac_data
        _, desc_eq = build_equivariant_library(snaps[0], 2 * np.pi, 2 * np.pi)
        w = _true_allen_cahn(desc_eq)
        assert len(w) == len(desc_eq)
        assert np.sum(np.abs(w) > 1e-10) == 3

    def test_true_fn_fails_with_standard(self, ac_data):
        """_true_allen_cahn must raise with standard library descriptions."""
        snaps, times, _, _ = ac_data
        _, desc_std = build_standard_library(snaps[0], 2 * np.pi, 2 * np.pi,
                                              max_order=4, max_poly=3)
        with pytest.raises((ValueError, KeyError)):
            _true_allen_cahn(desc_std)

    def test_evaluate_from_theta_dimension_mismatch_raises(self, ac_data):
        """Passing a 60-dim w (standard) with a 20-dim true_eq must raise."""
        snaps, times, true_eq, desc = ac_data
        Theta, Ft = build_eval_theta(snaps, times, 2 * np.pi, 2 * np.pi)
        # Simulate a standard-library coefficient vector (wrong size)
        w_wrong_size = np.zeros(60)
        with pytest.raises(ValueError):
            evaluate_from_theta(w_wrong_size, true_eq, Theta, Ft)


# ---------------------------------------------------------------------------
# Evaluate_recovery
# ---------------------------------------------------------------------------


class TestEvaluateRecovery:
    """Tests for evaluate_recovery()."""

    def test_exact_recovery_gives_small_E(self, ac_data):
        """w = true_eq should give E ≈ 0, exact_recovery = 1."""
        snaps, times, true_eq, desc = ac_data
        result = evaluate_recovery(true_eq, true_eq, snaps, times,
                                    2 * np.pi, 2 * np.pi)
        assert result["exact_recovery"] == 1
        assert result["fp"] == 0
        assert result["fn"] == 0
        assert result["E"] < 1e-8, f"E={result['E']} should be ~0"

    def test_zero_coefficients_gives_E_one(self, ac_data):
        """w = 0 should give E = 1."""
        snaps, times, true_eq, desc = ac_data
        w_zero = np.zeros_like(true_eq)
        result = evaluate_recovery(w_zero, true_eq, snaps, times,
                                    2 * np.pi, 2 * np.pi)
        # E = ||Theta·(0 - true)|| / ||Theta·true|| = 1 exactly
        assert abs(result["E"] - 1.0) < 1e-10

    def test_fp_counting(self, ac_data):
        """Adding a false positive should increase fp."""
        snaps, times, true_eq, desc = ac_data
        # Pick a true-inactive operator and give it a large coefficient
        w_bad = true_eq.copy()
        # Find an inactive operator (not f, f³, Δf)
        active = np.abs(true_eq) > 1e-10
        inactive_idx = np.where(~active)[0][0]
        w_bad[inactive_idx] = 0.5
        result = evaluate_recovery(w_bad, true_eq, snaps, times,
                                    2 * np.pi, 2 * np.pi)
        assert result["fp"] >= 1
        assert result["exact_recovery"] == 0

    def test_fn_counting(self, ac_data):
        """Zeroing a true term should increase fn."""
        snaps, times, true_eq, desc = ac_data
        w_bad = true_eq.copy()
        active = np.abs(true_eq) > 1e-10
        true_idx = np.where(active)[0][0]
        w_bad[true_idx] = 0.0
        result = evaluate_recovery(w_bad, true_eq, snaps, times,
                                    2 * np.pi, 2 * np.pi)
        assert result["fn"] >= 1
        assert result["exact_recovery"] == 0

    def test_n_active_count(self, ac_data):
        """n_active should count nonzero coefficients."""
        snaps, times, true_eq, desc = ac_data
        result = evaluate_recovery(true_eq, true_eq, snaps, times,
                                    2 * np.pi, 2 * np.pi)
        expected = int(np.sum(np.abs(true_eq) > 1e-10))
        assert result["n_active"] == expected

    def test_max_coeff_err_on_exact(self, ac_data):
        """max_coeff_err should be small for exact recovery."""
        snaps, times, true_eq, desc = ac_data
        result = evaluate_recovery(true_eq, true_eq, snaps, times,
                                    2 * np.pi, 2 * np.pi)
        assert result["max_coeff_err"] < 1e-10

    def test_max_coeff_err_on_miss(self, ac_data):
        """max_coeff_err should reflect missing true terms."""
        snaps, times, true_eq, desc = ac_data
        w_bad = true_eq.copy()
        active = np.abs(true_eq) > 1e-10
        true_idx = np.where(active)[0][0]
        w_bad[true_idx] = 0.0
        result = evaluate_recovery(w_bad, true_eq, snaps, times,
                                    2 * np.pi, 2 * np.pi)
        # Error on the zeroed term should be |true_coeff|
        expected_err = abs(true_eq[true_idx])
        assert abs(result["max_coeff_err"] - expected_err) < 1e-10

    def test_E_is_monotone_in_error(self, ac_data):
        """E should be larger when coefficients are worse."""
        snaps, times, true_eq, desc = ac_data
        # Small perturbation
        w_small = true_eq * 1.01
        # Large perturbation
        w_large = true_eq * 2.0
        E_small = evaluate_recovery(w_small, true_eq, snaps, times,
                                     2 * np.pi, 2 * np.pi)["E"]
        E_large = evaluate_recovery(w_large, true_eq, snaps, times,
                                     2 * np.pi, 2 * np.pi)["E"]
        assert E_large > E_small

    def test_n_samples_is_positive(self, ac_data):
        """n_samples should be positive (at least one row)."""
        snaps, times, true_eq, desc = ac_data
        result = evaluate_recovery(true_eq, true_eq, snaps, times,
                                    2 * np.pi, 2 * np.pi)
        assert result["n_samples"] > 0

    def test_skip_fraction_reduces_samples(self, ac_data):
        """skip_fraction=0.5 should give fewer samples than 0."""
        snaps, times, true_eq, desc = ac_data
        r0 = evaluate_recovery(true_eq, true_eq, snaps, times,
                                2 * np.pi, 2 * np.pi, skip_fraction=0.0)
        r5 = evaluate_recovery(true_eq, true_eq, snaps, times,
                                2 * np.pi, 2 * np.pi, skip_fraction=0.5)
        assert r5["n_samples"] < r0["n_samples"]


# ---------------------------------------------------------------------------
# best_of_sweep
# ---------------------------------------------------------------------------


class TestBestOfSweep:
    """Tests for best_of_sweep()."""

    def test_exact_beats_nonexact(self):
        """An exact recovery (fp=0, fn=0) should beat any non-exact."""
        rows = [
            dict(exact_recovery=0, fp=1, fn=0, E=0.001),
            dict(exact_recovery=1, fp=0, fn=0, E=0.01),
        ]
        best = best_of_sweep(rows)
        assert best["exact_recovery"] == 1

    def test_lower_fp_beats_higher_fp(self):
        """When neither is exact, lower fp wins."""
        rows = [
            dict(exact_recovery=0, fp=5, fn=0, E=0.001),
            dict(exact_recovery=0, fp=2, fn=0, E=0.1),
        ]
        best = best_of_sweep(rows)
        assert best["fp"] == 2

    def test_lower_fn_breaks_fp_tie(self):
        """When fp is tied, lower fn wins."""
        rows = [
            dict(exact_recovery=0, fp=2, fn=3, E=0.001),
            dict(exact_recovery=0, fp=2, fn=1, E=0.1),
        ]
        best = best_of_sweep(rows)
        assert best["fn"] == 1

    def test_lower_E_breaks_fn_tie(self):
        """When fp and fn are tied, lower E wins."""
        rows = [
            dict(exact_recovery=0, fp=2, fn=1, E=0.05),
            dict(exact_recovery=0, fp=2, fn=1, E=0.01),
        ]
        best = best_of_sweep(rows)
        assert best["E"] == 0.01

    def test_multiple_exact_all_tied(self):
        """When multiple rows are exact, fp/fn tie, pick lower E."""
        rows = [
            dict(exact_recovery=1, fp=0, fn=0, E=0.002),
            dict(exact_recovery=1, fp=0, fn=0, E=0.001),
        ]
        best = best_of_sweep(rows)
        assert best["E"] == 0.001

    def test_single_row(self):
        """Single row is trivially the best."""
        row = dict(exact_recovery=0, fp=3, fn=1, E=0.05)
        best = best_of_sweep([row])
        assert best is row

    def test_nan_fp_fn_handled(self):
        """Rows with NaN fp/fn are sorted after non-NaN rows.

        Since NaN comparisons in Python are unpredictable, the lexicographic
        criterion uses (exact_score, fp, fn, E).  NaN values in fp/fn
        should not crash the sort — they should be handled by the min()
        built-in, which treats NaN > any number in most implementations,
        but we guard by ensuring standard-library rows (fp=nan) don't beat
        valid rows.
        """
        # NaN fp should lose to valid fp=1
        rows = [
            dict(exact_recovery=0, fp=1, fn=0, E=0.01),
            dict(exact_recovery=0, fp=float("nan"), fn=float("nan"), E=0.001),
        ]
        best = best_of_sweep(rows)
        # The row with fp=1 should win (NaN goes to the end in min)
        assert best["fp"] == 1


# ---------------------------------------------------------------------------
# evaluate_standard
# ---------------------------------------------------------------------------


class TestEvaluateStandard:
    """Tests for evaluate_standard()."""

    def test_zero_coeffs_gives_E_fit_one(self, ac_data):
        """w=0 should give E_fit = 1 (predicts zero, target is nonzero)."""
        snaps, times, true_eq, desc = ac_data
        _, desc_std = build_standard_library(snaps[0], 2 * np.pi, 2 * np.pi,
                                              max_order=4, max_poly=3)
        # Build a small Theta_std on one snapshot
        T_std, _ = build_standard_library(snaps[1], 2 * np.pi, 2 * np.pi,
                                           max_order=4, max_poly=3)
        ft = (snaps[2] - snaps[0]) / (2 * (times[1] - times[0]))
        w_zero = np.zeros(T_std.shape[1])
        result = evaluate_standard(w_zero, T_std, ft.ravel())
        assert abs(result["E_fit"] - 1.0) < 1e-8
        assert result["n_active"] == 0

    def test_least_squares_fit_has_low_E_fit(self, ac_data):
        """The LSQ fit of the standard library should have E_fit ≪ 1
        because Allen-Cahn (f_t = Δf + f − f³) is representable in the
        standard basis."""
        snaps, times, true_eq, desc = ac_data
        T_std, _ = build_standard_library(snaps[1], 2 * np.pi, 2 * np.pi,
                                           max_order=4, max_poly=3)
        ft = (snaps[2] - snaps[0]) / (2 * (times[1] - times[0]))
        w_lsq = np.linalg.lstsq(T_std, ft.ravel(), rcond=None)[0]
        result = evaluate_standard(w_lsq, T_std, ft.ravel())
        # Allen-Cahn IS representable: f, f³, f_xx, f_yy all exist
        assert result["E_fit"] < 0.1, f"E_fit={result['E_fit']} should be small"
        assert result["n_active"] > 0

    def test_n_active_count(self, ac_data):
        """n_active should count nonzeros."""
        snaps, times, true_eq, desc = ac_data
        T_std, _ = build_standard_library(snaps[1], 2 * np.pi, 2 * np.pi,
                                           max_order=4, max_poly=3)
        ft = (snaps[2] - snaps[0]) / (2 * (times[1] - times[0]))
        w = np.zeros(T_std.shape[1])
        w[0] = 1.0
        w[3] = -0.5
        result = evaluate_standard(w, T_std, ft.ravel())
        assert result["n_active"] == 2
