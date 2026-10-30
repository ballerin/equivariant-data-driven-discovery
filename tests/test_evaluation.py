"""
Tests for ``evaluation.py`` — evaluate_recovery, evaluate_from_theta,
build_eval_theta_multigraphs, and best_of_sweep.

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
    evaluate_recovery, evaluate_from_theta, build_eval_theta_multigraphs,
    best_of_sweep, evaluate_standard, build_eval_theta_jet, evaluate_unified,
)
from equivariant_pde_find.library import (                                           # noqa: E402
    build_multigraph_library, build_pdefind_library,
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
    _, desc = build_multigraph_library(snaps[0], 2 * np.pi, 2 * np.pi)
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
        _, desc = build_multigraph_library(snaps[0], 2 * np.pi, 2 * np.pi)
        assert "Δf" in desc, f"Expected 'Δf' in {desc[:5]}..."
        assert "f³" in desc or "f\u00b3" in desc
        assert "Δ²f" in desc

    def test_standard_descriptions_do_not_contain_delta(self, ac_data):
        """Standard library must NOT contain 'Δf' (it has 'f_xx', etc.)."""
        snaps, times, _, _ = ac_data
        _, desc = build_pdefind_library(snaps[0], 2 * np.pi, 2 * np.pi,
                                          max_order=4, max_poly=3)
        assert "Δf" not in desc
        assert "f_xx" in desc
        assert "f_xxxx" in desc

    def test_true_fn_works_with_equivariant(self, ac_data):
        """_true_allen_cahn must succeed with equivariant descriptions."""
        snaps, times, _, _ = ac_data
        _, desc_eq = build_multigraph_library(snaps[0], 2 * np.pi, 2 * np.pi)
        w = _true_allen_cahn(desc_eq)
        assert len(w) == len(desc_eq)
        assert np.sum(np.abs(w) > 1e-10) == 3

    def test_true_fn_fails_with_standard(self, ac_data):
        """_true_allen_cahn must raise with standard library descriptions."""
        snaps, times, _, _ = ac_data
        _, desc_std = build_pdefind_library(snaps[0], 2 * np.pi, 2 * np.pi,
                                              max_order=4, max_poly=3)
        with pytest.raises((ValueError, KeyError)):
            _true_allen_cahn(desc_std)

    def test_evaluate_from_theta_dimension_mismatch_raises(self, ac_data):
        """Passing a 60-dim w (standard) with a 20-dim true_eq must raise."""
        snaps, times, true_eq, desc = ac_data
        Theta, Ft = build_eval_theta_multigraphs(
            snaps, times, 2 * np.pi, 2 * np.pi,
        )
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
        assert result["E_coeff"] < 1e-8, f"E={result['E_coeff']} should be ~0"

    def test_zero_coefficients_gives_E_one(self, ac_data):
        """w = 0 should give E = 1."""
        snaps, times, true_eq, desc = ac_data
        w_zero = np.zeros_like(true_eq)
        result = evaluate_recovery(w_zero, true_eq, snaps, times,
                                    2 * np.pi, 2 * np.pi)
        # E = ||Theta·(0 - true)|| / ||Theta·true|| = 1 exactly
        assert abs(result["E_coeff"] - 1.0) < 1e-10

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
                                     2 * np.pi, 2 * np.pi)["E_coeff"]
        E_large = evaluate_recovery(w_large, true_eq, snaps, times,
                                     2 * np.pi, 2 * np.pi)["E_coeff"]
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
            dict(exact_recovery=0, fp=1, fn=0, E_coeff=0.001),
            dict(exact_recovery=1, fp=0, fn=0, E_coeff=0.01),
        ]
        best = best_of_sweep(rows)
        assert best["exact_recovery"] == 1

    def test_lower_fp_beats_higher_fp(self):
        """When neither is exact, lower fp wins."""
        rows = [
            dict(exact_recovery=0, fp=5, fn=0, E_coeff=0.001),
            dict(exact_recovery=0, fp=2, fn=0, E_coeff=0.1),
        ]
        best = best_of_sweep(rows)
        assert best["fp"] == 2

    def test_lower_fn_breaks_fp_tie(self):
        """When fp is tied, lower fn wins."""
        rows = [
            dict(exact_recovery=0, fp=2, fn=3, E_coeff=0.001),
            dict(exact_recovery=0, fp=2, fn=1, E_coeff=0.1),
        ]
        best = best_of_sweep(rows)
        assert best["fn"] == 1

    def test_lower_E_breaks_fn_tie(self):
        """When fp and fn are tied, lower E wins."""
        rows = [
            dict(exact_recovery=0, fp=2, fn=1, E_coeff=0.05),
            dict(exact_recovery=0, fp=2, fn=1, E_coeff=0.01),
        ]
        best = best_of_sweep(rows)
        assert best["E_coeff"] == 0.01

    def test_multiple_exact_all_tied(self):
        """When multiple rows are exact, fp/fn tie, pick lower E."""
        rows = [
            dict(exact_recovery=1, fp=0, fn=0, E_coeff=0.002),
            dict(exact_recovery=1, fp=0, fn=0, E_coeff=0.001),
        ]
        best = best_of_sweep(rows)
        assert best["E_coeff"] == 0.001

    def test_single_row(self):
        """Single row is trivially the best."""
        row = dict(exact_recovery=0, fp=3, fn=1, E_coeff=0.05)
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
            dict(exact_recovery=0, fp=1, fn=0, E_coeff=0.01),
            dict(exact_recovery=0, fp=float("nan"), fn=float("nan"), E_coeff=0.001),
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
        _, desc_std = build_pdefind_library(snaps[0], 2 * np.pi, 2 * np.pi,
                                              max_order=4, max_poly=3)
        # Build a small Theta_std on one snapshot
        T_std, _ = build_pdefind_library(snaps[1], 2 * np.pi, 2 * np.pi,
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
        T_std, _ = build_pdefind_library(snaps[1], 2 * np.pi, 2 * np.pi,
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
        T_std, _ = build_pdefind_library(snaps[1], 2 * np.pi, 2 * np.pi,
                                           max_order=4, max_poly=3)
        ft = (snaps[2] - snaps[0]) / (2 * (times[1] - times[0]))
        w = np.zeros(T_std.shape[1])
        w[0] = 1.0
        w[3] = -0.5
        result = evaluate_standard(w, T_std, ft.ravel())
        assert result["n_active"] == 2


# ======================================================================
# Dynamic polyjet evaluation integration (M4)
# ======================================================================


class TestDynamicEvalTheta:
    """Tests for build_eval_theta_jet with jet_monomials parameter."""

    def test_legacy_behavior_unchanged(self, ac_data):
        """jet_monomials=None preserves old 105-column behavior."""
        snaps, times, true_eq, desc = ac_data
        T, Ft = build_eval_theta_jet(snaps, times, 2 * np.pi, 2 * np.pi,
                                      jet_monomials=None)
        assert T.shape[1] == 105, f"Expected 105 columns, got {T.shape[1]}"

    def test_dynamic_with_short_list(self, ac_data):
        """With explicit monomial list, only those columns are built."""
        snaps, times, true_eq, desc = ac_data
        monomials = ["1", "f", "f_x"]
        T, Ft = build_eval_theta_jet(snaps, times, 2 * np.pi, 2 * np.pi,
                                      jet_monomials=monomials)
        assert T.shape[1] == 3
        # Check column 0 is all ones
        assert np.allclose(T[:, 0], 1.0, atol=1e-12)

    def test_dynamic_spans_same_n_samples(self, ac_data):
        """Dynamic and legacy Thetas have the same number of rows."""
        snaps, times, true_eq, desc = ac_data
        monomials = ["1", "f", "f_x", "f_xx", "f³"]
        T_legacy, _ = build_eval_theta_jet(snaps, times, 2 * np.pi, 2 * np.pi)
        T_dynamic, _ = build_eval_theta_jet(snaps, times, 2 * np.pi, 2 * np.pi,
                                             jet_monomials=monomials)
        assert T_legacy.shape[0] == T_dynamic.shape[0]


class TestEvaluateUnifiedPolyjet:
    """Tests for evaluate_unified with library='polyjet'."""

    def test_polyjet_name_mapping_old_to_new(self, ac_data):
        """Regression in old 105-col polyjet, evaluation in dynamic subset:
        name-based mapping produces correct fp/fn/E."""
        snaps, times, true_eq, desc = ac_data

        # Build old regression polyjet (105 cols)
        from equivariant_pde_find.library import build_polyjet_library
        T_reg, desc_reg = build_polyjet_library(
            snaps[1], 2 * np.pi, 2 * np.pi, max_order=4, max_degree=3,
        )

        # Build dynamic evaluation jet (subset)
        from equivariant_pde_find.library.conversions import (
            collect_needed_jet_monomials_from_multigraph,
        )
        desc_jet = collect_needed_jet_monomials_from_multigraph(desc)
        T_jet, _ = build_eval_theta_jet(
            snaps, times, 2 * np.pi, 2 * np.pi, jet_monomials=desc_jet,
        )

        # Ground truth in jet
        from equivariant_pde_find.library.conversions import map_multigraph_to_jet
        true_jet = map_multigraph_to_jet(true_eq, desc, desc_jet)

        # AC in old polyjet: coefficients for f, f³, f_xx, f_yy
        w_reg = np.zeros(len(desc_reg))
        eps = 0.15
        w_reg[desc_reg.index("f")] = 1.0
        w_reg[desc_reg.index("f³")] = -1.0
        w_reg[desc_reg.index("f_xx")] = eps ** 2
        w_reg[desc_reg.index("f_yy")] = eps ** 2

        result = evaluate_unified(
            w_reg, "polyjet", desc_reg, desc_jet, T_jet, T_jet @ true_jet,
            true_jet,
        )
        assert result["exact_recovery"] == 1
        assert result["fp"] == 0
        assert result["fn"] == 0

    def test_polyjet_identity_when_same_sets(self, ac_data):
        """When regression and evaluation use the same dynamic set,
        the name-based mapping degenerates to identity."""
        snaps, times, true_eq, desc = ac_data

        from equivariant_pde_find.library.conversions import (
            collect_needed_jet_monomials_from_multigraph,
        )
        desc_jet = collect_needed_jet_monomials_from_multigraph(desc)
        T_jet, _ = build_eval_theta_jet(
            snaps, times, 2 * np.pi, 2 * np.pi, jet_monomials=desc_jet,
        )

        from equivariant_pde_find.library.conversions import map_multigraph_to_jet
        true_jet = map_multigraph_to_jet(true_eq, desc, desc_jet)

        # Recovered w is already in the dynamic jet
        w = true_jet.copy()

        result = evaluate_unified(
            w, "polyjet", desc_jet, desc_jet, T_jet, T_jet @ true_jet,
            true_jet,
        )
        assert result["exact_recovery"] == 1
        assert result["fp"] == 0
        assert result["fn"] == 0


# ======================================================================


class TestSpectralDefaultEvaluation:
    """Verify that evaluation builders always use spectral derivatives
    regardless of the globally configured backend."""

    def test_jet_defaults_to_spectral_when_polynomial_is_active(self, ac_data):
        """build_eval_theta_jet uses spectral even when global backend
        is polynomial."""
        snaps, times, true_eq, desc = ac_data
        from equivariant_pde_find import derivatives

        # Force polynomial as the global backend
        derivatives.set_backend("polynomial")
        try:
            # Build WITHOUT explicit backend — should still use spectral
            T_jet, _ = build_eval_theta_jet(
                snaps, times, 2 * np.pi, 2 * np.pi, jet_monomials=None,
            )
            # Build WITH explicit spectral for comparison
            T_spec, _ = build_eval_theta_jet(
                snaps, times, 2 * np.pi, 2 * np.pi,
                backend="spectral", jet_monomials=None,
            )
        finally:
            derivatives.set_backend("spectral")

        # They should match (both spectral)
        assert np.allclose(T_jet, T_spec, atol=1e-12), (
            "build_eval_theta_jet did not default to spectral"
        )

    def test_jet_default_spectral_differs_from_polynomial(self, ac_data):
        """Explicitly confirm that spectral vs polynomial eval Thetas
        differ somewhere (higher-order derivatives), to guard against
        a degenerate test where both backends accidentally match."""
        snaps, times, true_eq, desc = ac_data

        T_spec, _ = build_eval_theta_jet(
            snaps, times, 2 * np.pi, 2 * np.pi,
            backend="spectral", jet_monomials=None,
        )
        T_poly, _ = build_eval_theta_jet(
            snaps, times, 2 * np.pi, 2 * np.pi,
            backend="polynomial", jet_monomials=None,
        )

        # Same shape
        assert T_spec.shape == T_poly.shape

        # The "f" column is backend-independent (no derivatives taken)
        assert np.allclose(T_spec[:, 0], T_poly[:, 0], atol=1e-12)

        # Higher-order derivative columns should differ at machine precision
        # because the polynomial stencil is not exact.  Check the full
        # matrices differ (excluding the field-only column).
        max_diff = np.max(np.abs(T_spec[:, 1:] - T_poly[:, 1:]))
        assert max_diff > 1e-15, (
            f"Polynomial and spectral Thetas are identical to machine "
            f"precision (max_diff={max_diff}); test is degenerate"
        )

    def test_multigraphs_defaults_to_spectral(self, ac_data):
        """build_eval_theta_multigraphs defaults to spectral."""
        snaps, times, true_eq, desc = ac_data
        from equivariant_pde_find import derivatives

        derivatives.set_backend("polynomial")
        try:
            T_def, _ = build_eval_theta_multigraphs(
                snaps, times, 2 * np.pi, 2 * np.pi,
            )
            T_spec, _ = build_eval_theta_multigraphs(
                snaps, times, 2 * np.pi, 2 * np.pi, backend="spectral",
            )
        finally:
            derivatives.set_backend("spectral")

        assert np.allclose(T_def, T_spec, atol=1e-12), (
            "build_eval_theta_multigraphs did not default to spectral"
        )

    def test_evaluate_recovery_uses_spectral_eval_theta(self, ac_data):
        """evaluate_recovery computes E_coeff against spectral eval Theta
        by default, so exact ground-truth coefficients yield E_coeff ≈ 0."""
        snaps, times, true_eq, desc = ac_data
        from equivariant_pde_find import derivatives

        # Force polynomial as global backend — should not affect evaluation
        derivatives.set_backend("polynomial")
        try:
            result = evaluate_recovery(
                true_eq, true_eq, snaps, times, 2 * np.pi, 2 * np.pi,
            )
        finally:
            derivatives.set_backend("spectral")

        assert result["exact_recovery"] == 1
        assert result["fp"] == 0
        assert result["fn"] == 0
        assert result["E_coeff"] < 1e-8, (
            f"E_coeff={result['E_coeff']} should be near 0 for exact recovery"
        )

    def test_evaluate_recovery_spectral_vs_polynomial_yield_different_E(
        self, ac_data,
    ):
        """E_coeff for the SAME perturbed coefficients differs between
        spectral and polynomial eval Thetas, confirming the backend
        choice matters for the metric."""
        snaps, times, true_eq, desc = ac_data

        # Perturb the coefficients slightly so E_coeff > 0
        w_perturbed = true_eq.copy()
        # Add 1 % error to the first non-zero coefficient
        nz = np.where(np.abs(true_eq) > 1e-10)[0]
        w_perturbed[nz[0]] *= 1.01

        result_spec = evaluate_recovery(
            w_perturbed, true_eq, snaps, times, 2 * np.pi, 2 * np.pi,
            backend="spectral",
        )
        result_poly = evaluate_recovery(
            w_perturbed, true_eq, snaps, times, 2 * np.pi, 2 * np.pi,
            backend="polynomial",
        )

        # Both should have E_coeff > 0 (coefficients are wrong)
        assert result_spec["E_coeff"] > 0
        assert result_poly["E_coeff"] > 0

        # They should differ because the eval Thetas differ
        assert not np.isclose(
            result_spec["E_coeff"], result_poly["E_coeff"], rtol=0, atol=0,
        ), (
            f"Spectral E_coeff ({result_spec['E_coeff']:.6e}) and polynomial "
            f"E_coeff ({result_poly['E_coeff']:.6e}) should differ"
        )


# ---------------------------------------------------------------------------
# FactoredEval — QR stand-in for a large evaluation (Theta, Ft)
# ---------------------------------------------------------------------------

from equivariant_pde_find.evaluation import FactoredEval                             # noqa: E402


def _ill_conditioned_pair(n, K, seed=0):
    rng = np.random.default_rng(seed)
    Theta = rng.standard_normal((n, K)) * np.logspace(0, 6, K)
    Theta[:, K // 2] = 3.0 * Theta[:, 0] + 1e-6 * rng.standard_normal(n)
    true = np.zeros(K)
    true[[1, 5, K // 2]] = [1.0, -0.5, 2.0]
    Ft = Theta @ true
    Ft = Ft + 0.01 * np.linalg.norm(Ft) / np.sqrt(n) * rng.standard_normal(n)
    return Theta, Ft, true


@pytest.mark.parametrize("n,K,chunk_rows", [(5000, 40, 1 << 18), (5000, 40, 700),
                                            (30, 40, 8)])
def test_factored_eval_matches_dense(n, K, chunk_rows):
    """Same metrics as the dense evaluation, to round-off; structural metrics
    exactly.  Covers several chunk sizes, including fewer rows than columns."""
    Theta, Ft, true = _ill_conditioned_pair(n, K)
    fe = FactoredEval(Theta, Ft, chunk_rows=chunk_rows)
    assert fe.shape == (n, K)
    rng = np.random.default_rng(1)
    candidates = [true.copy(), true + 1e-8 * np.eye(K)[3],
                  true * (1 + 1e-3 * rng.standard_normal(K)),
                  rng.standard_normal(K) * (rng.random(K) < 0.2)]
    for w in candidates:
        dense = evaluate_from_theta(w, true, Theta, Ft)
        fact = evaluate_from_theta(w, true, fe, None)
        assert dense.keys() == fact.keys()
        for key in ("fp", "fn", "exact_recovery", "n_active", "n_samples",
                    "max_coeff_err"):
            assert dense[key] == fact[key], key
        for key in ("E_coeff", "E_data"):
            assert fact[key] == pytest.approx(dense[key], rel=1e-9, abs=1e-15), key
        assert evaluate_standard(w, fe, None)["E_fit"] == pytest.approx(
            evaluate_standard(w, Theta, Ft)["E_fit"], rel=1e-9)


def test_factored_eval_unrepresentable_truth_branch():
    Theta, Ft, _ = _ill_conditioned_pair(400, 12)
    w = np.random.default_rng(2).standard_normal(12)
    zero_truth = np.zeros(12)
    dense = evaluate_from_theta(w, zero_truth, Theta, Ft)
    fact = evaluate_from_theta(w, zero_truth, FactoredEval(Theta, Ft), None)
    assert fact["E_fit"] == pytest.approx(dense["E_fit"], rel=1e-9)
    assert fact["n_samples"] == dense["n_samples"] == 400


def test_factored_eval_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        FactoredEval(np.ones((10, 3)), np.ones(9))
