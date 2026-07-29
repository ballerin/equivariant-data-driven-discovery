"""Tests for the N-D generalization of multigraph_operators.py (max_p<=2).

Covers:
- All max_p<=2 multigraph/mod-multigraph operators run and match analytic
  single-Fourier-mode values at ndim=4,5.
- Operator count is dimension-independent (20 for max_p=2, max_n=3).
- max_p=3 operators raise a clear NotImplementedError outside ndim in (2,3).
- Numeric evaluation matches the symbolic jet_algebra formulas directly
  (cross-check between the two "modes" of the same formula).
"""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.multigraph_operators import (
    get_operators,
    get_mod_multigraph_operators,
    _evaluate_jet_expr,
    _sym_hess_norm_sq,
)


def _single_mode_field(ndim, N, L, ks):
    x = np.linspace(0, L, N, endpoint=False)
    grids = np.meshgrid(*([x] * ndim), indexing="ij")
    kvec = 2.0 * math.pi / L * np.asarray(ks, dtype=float)
    theta = sum(k * g for k, g in zip(kvec, grids))
    return np.cos(theta), kvec, theta, grids


class TestOperatorCountIsDimensionIndependent:
    @pytest.mark.parametrize("ndim", [2, 3, 4, 5])
    def test_multigraph_count(self, ndim):
        assert len(get_operators(max_p=2, max_n=3)) == 20

    @pytest.mark.parametrize("ndim", [2, 3, 4, 5])
    def test_mod_multigraph_count(self, ndim):
        assert len(get_mod_multigraph_operators(max_p=2, max_n=3)) == 20


class TestAllMaxP2OperatorsRunAtHighNdim:
    @pytest.mark.parametrize("ndim", [4, 5])
    def test_multigraph_operators_finite_and_shaped(self, ndim):
        N, L = 8, 2 * math.pi
        rng = np.random.default_rng(0)
        f = 0.3 * rng.standard_normal((N,) * ndim)
        for func, desc in get_operators(max_p=2, max_n=3):
            val = func(f, L, L)
            assert val.shape == f.shape, desc
            assert np.all(np.isfinite(val)), desc

    @pytest.mark.parametrize("ndim", [4, 5])
    def test_mod_multigraph_operators_finite_and_shaped(self, ndim):
        N, L = 8, 2 * math.pi
        rng = np.random.default_rng(1)
        f = 0.3 * rng.standard_normal((N,) * ndim)
        for func, desc in get_mod_multigraph_operators(max_p=2, max_n=3):
            val = func(f, L, L)
            assert val.shape == f.shape, desc
            assert np.all(np.isfinite(val)), desc


class TestAnalyticCorrectnessAtNdim5:
    """Single-mode f = cos(k.x) has closed-form values for every p<=2 operator."""

    def setup_field(self):
        ndim, N, L = 5, 10, 2 * math.pi
        ks = [1, 1, 0, 1, 0]
        f, kvec, theta, grids = _single_mode_field(ndim, N, L, ks)
        return f, kvec, theta, L

    def test_Delta_f(self):
        f, kvec, theta, L = self.setup_field()
        ops = dict((d, fn) for fn, d in get_operators(2, 3))
        val = ops["Δf"](f, L, L)
        expected = -np.sum(kvec**2) * f
        np.testing.assert_allclose(val, expected, atol=1e-9)

    def test_grad_sq(self):
        f, kvec, theta, L = self.setup_field()
        ops = dict((d, fn) for fn, d in get_operators(2, 3))
        val = ops["‖∇f‖²"](f, L, L)
        expected = np.sum(kvec**2) * np.sin(theta) ** 2
        np.testing.assert_allclose(val, expected, atol=1e-9)

    def test_bilaplacian(self):
        f, kvec, theta, L = self.setup_field()
        ops = dict((d, fn) for fn, d in get_operators(2, 3))
        val = ops["Δ²f"](f, L, L)
        expected = (np.sum(kvec**2) ** 2) * f
        np.testing.assert_allclose(val, expected, atol=1e-8)

    def test_hess_norm_sq(self):
        f, kvec, theta, L = self.setup_field()
        ops = dict((d, fn) for fn, d in get_operators(2, 3))
        val = ops["‖∇²f‖²"](f, L, L)
        # f_ij = -k_i k_j cos(theta) for a single mode -> sum_ij f_ij^2 = |k|^4 cos^2
        expected = (np.sum(kvec**2) ** 2) * np.cos(theta) ** 2
        np.testing.assert_allclose(val, expected, atol=1e-8)


class TestAllOperatorsAnalyticAtHighNdim:
    """Closed-form check of **all 20** max_p<=2 operators at ndim=4 and 5.

    ``TestAnalyticCorrectnessAtNdim5`` above covers only 4 operators, and
    beyond ndim=3 the domain is isotropic — every axis has the same length —
    so an axis-ordering/transposition bug in the multi-index bookkeeping is
    invisible unless the *field* distinguishes the axes.  Here the single
    mode uses **distinct, non-symmetric** wavenumbers per axis, so any axis
    permutation inside the derivative machinery changes the answer.

    Ground truth is assembled from the exact single-mode derivative
    ``∂^α cos(k·x) = (∏ k_i^{α_i}) · cos^{(|α|)}(k·x)`` and the operators'
    mathematical definitions — it never calls the derivative backend or
    ``_evaluate_jet_expr``, so agreement is real evidence rather than a
    self-comparison.
    """

    @staticmethod
    def _analytic_ops(kvec, theta):
        """{operator description: exact array} for a field f = cos(k·x)."""
        ndim = len(kvec)
        # cos^{(m)}: m = 0,1,2,3 (mod 4)
        cycle = [np.cos(theta), -np.sin(theta), -np.cos(theta), np.sin(theta)]

        def d(*axes):
            """Exact ∂^α cos(k·x) for α given as a list of (repeated) axes."""
            factor = 1.0
            for a in axes:
                factor *= kvec[a]
            return factor * cycle[len(axes) % 4]

        f = np.cos(theta)
        lap = sum(d(i, i) for i in range(ndim))
        bilap = sum(d(i, i, j, j) for i in range(ndim) for j in range(ndim))
        grad_sq = sum(d(i) ** 2 for i in range(ndim))
        grad_lap_grad = sum(
            sum(d(j, j, i) for j in range(ndim)) * d(i) for i in range(ndim)
        )
        hess_norm_sq = sum(
            d(i, j) ** 2 for i in range(ndim) for j in range(ndim)
        )
        hess_grad_grad = sum(
            d(i, j) * d(i) * d(j) for i in range(ndim) for j in range(ndim)
        )
        return {
            "1": np.ones_like(f),
            "f": f,
            "f²": f ** 2,
            "f³": f ** 3,
            "Δf": lap,
            "‖∇f‖²": grad_sq,
            "f·‖∇f‖²": f * grad_sq,
            "Δ²f": bilap,
            "Δf·Δf": lap * lap,
            "‖∇²f‖²": hess_norm_sq,
            "f·Δf·Δf": f * lap * lap,
            "f·‖∇²f‖²": f * hess_norm_sq,
            "∇²f(∇f,∇f)": hess_grad_grad,
            "‖∇f‖²·Δf": grad_sq * lap,
            "f·Δf": f * lap,
            "f²·Δf": f ** 2 * lap,
            "f·Δ²f": f * bilap,
            "⟨∇Δf,∇f⟩": grad_lap_grad,
            "f²·Δ²f": f ** 2 * bilap,
            "f·⟨∇Δf,∇f⟩": f * grad_lap_grad,
        }

    @pytest.mark.parametrize("ndim,ks", [
        (4, [1, 2, 3, 1]),
        (5, [1, 3, 2, 1, 2]),
    ])
    def test_all_20_multigraph_operators(self, ndim, ks):
        N, L = 8, 2 * math.pi
        f, kvec, theta, _ = _single_mode_field(ndim, N, L, ks)
        expected = self._analytic_ops(kvec, theta)

        ops = get_operators(2, 3)
        assert len(ops) == 20
        for func, desc in ops:
            assert desc in expected, f"no analytic form recorded for {desc!r}"
            value = func(f, L, L)
            scale = max(np.max(np.abs(expected[desc])), 1.0)
            np.testing.assert_allclose(
                value, expected[desc], atol=1e-8 * scale, rtol=0,
                err_msg=f"ndim={ndim} operator {desc!r}",
            )

    @pytest.mark.parametrize("ndim,ks", [
        (4, [1, 2, 3, 1]),
        (5, [1, 3, 2, 1, 2]),
    ])
    def test_all_20_mod_multigraph_operators(self, ndim, ks):
        """The 6 replacement operators are Δ(·)/Δ²(·) forms; their exact
        values follow from the same single-mode derivative rule applied to
        f², f³ and ‖∇f‖²."""
        N, L = 8, 2 * math.pi
        f, kvec, theta, _ = _single_mode_field(ndim, N, L, ks)
        shared = self._analytic_ops(kvec, theta)

        K = float(np.sum(kvec ** 2))
        cos, sin = np.cos(theta), np.sin(theta)
        # Δ(f²) = 2fΔf + 2‖∇f‖²;  Δ(f³) = 3f²Δf + 6f‖∇f‖²
        lap_f2 = 2 * cos * (-K * cos) + 2 * K * sin ** 2
        lap_f3 = 3 * cos ** 2 * (-K * cos) + 6 * cos * K * sin ** 2
        # Δ‖∇f‖² = 2⟨∇Δf,∇f⟩ + 2‖∇²f‖²
        lap_grad_sq = 2 * shared["⟨∇Δf,∇f⟩"] + 2 * shared["‖∇²f‖²"]
        # Δ(fg) = fΔg + 2⟨∇f,∇g⟩ + gΔf, with g = ‖∇f‖² and
        # ∂_i g = 2K k_i sinθ cosθ  =>  ⟨∇f,∇g⟩ = -2K² sin²θ cosθ
        lap_f_grad_sq = (
            cos * lap_grad_sq
            + 2 * (-2 * K ** 2 * sin ** 2 * cos)
            + (K * sin ** 2) * (-K * cos)
        )
        # Δ² acts on harmonics: Δ(cos mθ) = -m²K cos mθ.
        # Δ(f²) = -2K cos2θ  =>  Δ²(f²) = 8K² cos2θ
        bilap_f2 = 8 * K ** 2 * np.cos(2 * theta)
        # Δ(f³) = -(3K/4)(cosθ + 3cos3θ)  =>  Δ²(f³) = (3K²/4)(cosθ + 27cos3θ)
        bilap_f3 = (3 * K ** 2 / 4) * (cos + 27 * np.cos(3 * theta))

        expected = {k: shared[k] for k in [
            "1", "f", "f²", "f³", "Δf", "‖∇f‖²", "f·‖∇f‖²", "Δ²f", "Δf·Δf",
            "‖∇²f‖²", "f·Δf·Δf", "f·‖∇²f‖²", "∇²f(∇f,∇f)", "‖∇f‖²·Δf",
        ]}
        expected.update({
            "Δf²": lap_f2,
            "Δf³": lap_f3,
            "Δ‖∇f‖²": lap_grad_sq,
            "Δ(f·‖∇f‖²)": lap_f_grad_sq,
            "Δ²f²": bilap_f2,
            "Δ²f³": bilap_f3,
        })

        ops = get_mod_multigraph_operators(2, 3)
        assert len(ops) == 20
        for func, desc in ops:
            assert desc in expected, f"no analytic form recorded for {desc!r}"
            value = func(f, L, L)
            scale = max(np.max(np.abs(expected[desc])), 1.0)
            np.testing.assert_allclose(
                value, expected[desc], atol=1e-8 * scale, rtol=0,
                err_msg=f"ndim={ndim} operator {desc!r}",
            )

    def test_wavenumbers_are_distinct_enough_to_expose_axis_swaps(self):
        """Guard: the chosen wavenumbers must make the field sensitive to an
        axis permutation, otherwise this whole class proves nothing about
        axis bookkeeping."""
        ndim, N, L, ks = 4, 8, 2 * math.pi, [1, 2, 3, 1]
        f, _, _, _ = _single_mode_field(ndim, N, L, ks)
        swapped = np.transpose(f, (1, 0, 2, 3))
        assert not np.allclose(f, swapped, atol=1e-12)


class TestMaxP3RaisesOutsideNdim2Or3:
    @pytest.mark.parametrize("ndim", [1, 4, 5])
    def test_p3_operator_raises(self, ndim):
        N, L = 6, 2 * math.pi
        rng = np.random.default_rng(0)
        f = rng.standard_normal((N,) * ndim)
        ops = dict((d, fn) for fn, d in get_operators(max_p=3, max_n=4))
        with pytest.raises(NotImplementedError, match="only implemented for"):
            ops["Δ³f"](f, L, L)

    def test_p3_operator_still_works_at_2d_3d(self):
        for ndim, L in [(2, 2 * math.pi), (3, 2 * math.pi)]:
            N = 8
            rng = np.random.default_rng(0)
            f = rng.standard_normal((N,) * ndim)
            ops = dict((d, fn) for fn, d in get_operators(max_p=3, max_n=4))
            val = ops["Δ³f"](f, L, L, L if ndim == 3 else None)
            assert val.shape == f.shape
            assert np.all(np.isfinite(val))


class TestNumericMatchesSymbolic:
    """Cross-check: evaluating the symbolic JetExpr numerically for a random
    polynomial-in-jet-variables input reproduces the array from evaluating
    factor-by-factor by hand -- guards against a mismatch between JetExpr's
    monomial bookkeeping and _evaluate_jet_expr's substitution."""

    def test_hess_norm_sq_matches_manual_sum(self):
        ndim, N, L = 4, 8, 2 * math.pi
        rng = np.random.default_rng(0)
        f = 0.2 * rng.standard_normal((N,) * ndim)

        expr = _sym_hess_norm_sq(ndim)
        val = _evaluate_jet_expr(expr, f, L, L, L)

        from equivariant_pde_find import derivatives
        Ls = derivatives.axis_lengths(f, L, L, L)

        def mi_ij(i, j):
            m = [0] * ndim
            m[i] += 1
            m[j] += 1
            return tuple(m)

        manual = np.zeros_like(f)
        for i in range(ndim):
            for j in range(ndim):
                d = derivatives.partial_nd(f, Ls, mi_ij(i, j))
                manual += d * d
        np.testing.assert_allclose(val, manual, atol=1e-10)
