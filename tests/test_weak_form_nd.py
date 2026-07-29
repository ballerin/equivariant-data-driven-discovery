"""Tests for the N-D generalization of library/weak_form.py.

Covers:
- weak_convolution accepts a list of ndim kernels (new calling convention)
  in addition to the historical positional style.
- build_weak_rhs works at ndim=4 with default (m, p) fill-in.
- weak_laplacian_ibp / weak_bilaplacian_ibp agree with the direct
  (derivative-then-convolve) approach at ndim=4, within the same tolerance
  bands already accepted for the existing 2-D (<10%) and 3-D (<30%) tests.
"""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives
from equivariant_pde_find.library import weak_form as wf


class TestWeakConvolutionND:
    def test_list_and_positional_conventions_agree(self):
        N, m = 16, 3
        rng = np.random.default_rng(0)
        f = rng.standard_normal((N, N, N))
        k = wf._test_function_1d(m, 4)
        via_positional = wf.weak_convolution(f, k, k, k)
        via_list = wf.weak_convolution(f, [k, k, k])
        np.testing.assert_allclose(via_positional, via_list)

    def test_ndim_4(self):
        N, m = 10, 2
        rng = np.random.default_rng(0)
        f = rng.standard_normal((N,) * 4)
        k = wf._test_function_1d(m, 4)
        out = wf.weak_convolution(f, [k] * 4)
        assert out.shape == (N - 2 * m,) * 4
        assert np.all(np.isfinite(out))

    def test_wrong_kernel_count_raises(self):
        f = np.ones((8, 8, 8, 8))
        k = wf._test_function_1d(3, 4)
        with pytest.raises(ValueError, match="need exactly f.ndim kernels"):
            wf.weak_convolution(f, [k, k])


class TestBuildWeakRhsND:
    def test_ndim_4_default_fill(self):
        N, m, p = 12, 3, 4
        rng = np.random.default_rng(0)
        ft = rng.standard_normal((N,) * 4)
        b = wf.build_weak_rhs(ft, m=m, p=p)
        assert b.shape == ((N - 2 * m) ** 4,)
        assert np.all(np.isfinite(b))

    def test_explicit_kernels_list(self):
        N, m, p = 12, 3, 4
        rng = np.random.default_rng(0)
        ft = rng.standard_normal((N,) * 4)
        k = wf._cached_kernel_1d(m, p)
        b1 = wf.build_weak_rhs(ft, m=m, p=p)
        b2 = wf.build_weak_rhs(ft, kernels=[k, k, k, k])
        np.testing.assert_allclose(b1, b2)


class TestIBPAgreesWithDirectAtNdim4:
    def setup_field(self):
        derivatives.set_backend("spectral")
        ndim, N, L = 4, 16, 2 * math.pi
        rng = np.random.default_rng(0)
        f = 0.3 * rng.standard_normal((N,) * ndim)
        return f, L

    def test_laplacian_ibp(self):
        f, L = self.setup_field()
        m, p = 6, 6
        lap = derivatives.laplacian(f, L, L)
        kernel = wf._cached_kernel_1d(m, p)
        direct = wf.weak_convolution(lap, [kernel] * f.ndim)
        ibp = wf.weak_laplacian_ibp(f, m, p, L, L)
        rel = np.linalg.norm(ibp - direct) / np.linalg.norm(direct)
        assert rel < 0.1, f"relative error {rel:.4f} exceeds 10%"

    def test_bilaplacian_ibp(self):
        f, L = self.setup_field()
        m, p = 6, 6
        bilap = derivatives.bilaplacian(f, L, L)
        kernel = wf._cached_kernel_1d(m, p)
        direct = wf.weak_convolution(bilap, [kernel] * f.ndim)
        ibp = wf.weak_bilaplacian_ibp(f, m, p, L, L)
        rel = np.linalg.norm(ibp - direct) / np.linalg.norm(direct)
        assert rel < 0.3, f"relative error {rel:.4f} exceeds 30%"


class TestIBPKernelStructureIsExact:
    """Structural (machine-precision) checks of the N-D IBP kernel expansions.

    The agreement tests above allow 10-30% relative error, because the IBP
    weak form and the direct derivative-then-convolve form are genuinely
    different discretizations.  That tolerance is wide enough to hide a
    wrong *term count* in the ndim-generic expansion — e.g. dropping the
    mixed ``φ″(i)φ″(j)`` terms of Δ²ψ, or getting their coefficient wrong.

    Here the IBP result is instead compared against a separable kernel sum
    assembled independently from the definitions

        Δψ  = Σ_i φ″(i)/dx_i² · Π_{j≠i} φ(j)
        Δ²ψ = Σ_i φ⁗(i)/dx_i⁴ · Π_{j≠i} φ(j)
              + 2 Σ_{i<j} φ″(i)φ″(j)/(dx_i²dx_j²) · Π_{k∉{i,j}} φ(k)

    so agreement must be exact.  At ndim=4 that is 4 pure + 6 mixed terms;
    a missing or mis-weighted mixed term fails immediately.
    """

    def _setup(self, ndim=4, N=12, L=2 * math.pi, seed=5):
        derivatives.set_backend("spectral")
        rng = np.random.default_rng(seed)
        return rng.standard_normal((N,) * ndim), L

    def test_laplacian_kernel_expansion(self):
        f, L = self._setup()
        ndim, m, p = f.ndim, 3, 6
        Ls = derivatives.axis_lengths(f, L, L)
        dx2 = [(Ls[a] / f.shape[a]) ** 2 for a in range(ndim)]
        phi = wf._test_function_1d(m, p)
        d2 = wf._test_function_deriv_1d(m, p, 2)

        expected = None
        for i in range(ndim):
            kernels = [phi] * ndim
            kernels[i] = d2 / dx2[i]
            term = wf.weak_convolution(f, kernels)
            expected = term if expected is None else expected + term

        got = wf.weak_laplacian_ibp(f, m, p, L, L)
        np.testing.assert_allclose(got, expected, rtol=1e-12, atol=0)

    def test_bilaplacian_kernel_expansion_has_all_mixed_terms(self):
        f, L = self._setup()
        ndim, m, p = f.ndim, 3, 6
        Ls = derivatives.axis_lengths(f, L, L)
        dx2 = [(Ls[a] / f.shape[a]) ** 2 for a in range(ndim)]
        phi = wf._test_function_1d(m, p)
        d2 = wf._test_function_deriv_1d(m, p, 2)
        d4 = wf._test_function_deriv_1d(m, p, 4)

        expected = None
        n_terms = 0
        for i in range(ndim):  # pure φ⁗ terms
            kernels = [phi] * ndim
            kernels[i] = d4 / dx2[i] ** 2
            term = wf.weak_convolution(f, kernels)
            expected = term if expected is None else expected + term
            n_terms += 1
        for i in range(ndim):  # mixed 2·φ″φ″ terms
            for j in range(i + 1, ndim):
                kernels = [phi] * ndim
                kernels[i] = d2 / dx2[i]
                kernels[j] = d2 / dx2[j]
                expected = expected + 2.0 * wf.weak_convolution(f, kernels)
                n_terms += 1
        assert n_terms == ndim + ndim * (ndim - 1) // 2 == 10

        got = wf.weak_bilaplacian_ibp(f, m, p, L, L)
        np.testing.assert_allclose(got, expected, rtol=1e-12, atol=0)

    def test_dropping_a_mixed_term_would_be_detected(self):
        """Guard: the exact comparison above must be sensitive to the mixed
        terms it is meant to protect (they are not numerically negligible)."""
        f, L = self._setup()
        ndim, m, p = f.ndim, 3, 6
        Ls = derivatives.axis_lengths(f, L, L)
        dx2 = [(Ls[a] / f.shape[a]) ** 2 for a in range(ndim)]
        phi = wf._test_function_1d(m, p)
        d2 = wf._test_function_deriv_1d(m, p, 2)

        kernels = [phi] * ndim
        kernels[0] = d2 / dx2[0]
        kernels[1] = d2 / dx2[1]
        one_mixed_term = 2.0 * wf.weak_convolution(f, kernels)
        full = wf.weak_bilaplacian_ibp(f, m, p, L, L)
        assert np.linalg.norm(one_mixed_term) > 1e-3 * np.linalg.norm(full)


class TestWeakFormExactOnSingleMode:
    """On a band-limited field the *direct* weak form is machine-exact, so it
    is a genuine analytic reference for the whole ndim-generic pipeline."""

    def _mode(self, ndim, N, L, ks):
        x = np.linspace(0.0, L, N, endpoint=False)
        grids = np.meshgrid(*([x] * ndim), indexing="ij")
        kvec = 2.0 * math.pi / L * np.asarray(ks, dtype=float)
        theta = sum(k * g for k, g in zip(kvec, grids))
        return np.cos(theta), float(np.sum(kvec ** 2))

    def test_laplacian_weak_form_matches_analytic_delta_f(self):
        """∫ψ·Δf for f=cos(k·x) must equal −|k|²·∫ψ·f, since Δf = −|k|²f
        exactly.  This checks the convolution path end to end at ndim=4."""
        derivatives.set_backend("spectral")
        ndim, N, L, m, p = 4, 16, 2 * math.pi, 4, 6
        f, K = self._mode(ndim, N, L, [1, 2, 1, 2])
        kernel = wf._cached_kernel_1d(m, p)

        lap = derivatives.laplacian(f, L, L)
        np.testing.assert_allclose(lap, -K * f, atol=1e-9)

        direct = wf.weak_convolution(lap, [kernel] * ndim)
        analytic = -K * wf.weak_convolution(f, [kernel] * ndim)
        np.testing.assert_allclose(direct, analytic, rtol=1e-10, atol=1e-12)

    def test_ibp_converges_to_analytic_as_m_grows(self):
        """The IBP form approximates the same integral; its relative error
        against the analytic value must shrink as the test function support
        widens (more grid points resolving φ)."""
        derivatives.set_backend("spectral")
        ndim, N, L, p = 4, 20, 2 * math.pi, 6
        f, K = self._mode(ndim, N, L, [1, 1, 1, 1])

        errors = []
        for m in (3, 6):
            kernel = wf._cached_kernel_1d(m, p)
            analytic = -K * wf.weak_convolution(f, [kernel] * ndim)
            ibp = wf.weak_laplacian_ibp(f, m, p, L, L)
            errors.append(
                np.linalg.norm(ibp - analytic) / np.linalg.norm(analytic)
            )
        assert errors[1] < errors[0], f"IBP error did not improve with m: {errors}"
        assert errors[1] < 0.02, f"IBP error at m=6 too large: {errors[1]:.4f}"


class TestComputeMScale:
    """First tests for compute_m_scale (previously untested at any ndim).

    Per Messenger & Bortz's WSINDy pipeline it is M_j = ‖G_j‖ / ‖b‖, the
    per-column scale feeding MSTLS's dominant-balance bounds.
    """

    def test_matches_definition(self):
        rng = np.random.default_rng(0)
        G = rng.standard_normal((50, 7))
        b = rng.standard_normal(50)
        expected = np.linalg.norm(G, axis=0) / np.linalg.norm(b)
        np.testing.assert_allclose(wf.compute_m_scale(G, b), expected)

    def test_shape_and_positivity(self):
        rng = np.random.default_rng(1)
        G = rng.standard_normal((30, 5))
        b = rng.standard_normal(30)
        M = wf.compute_m_scale(G, b)
        assert M.shape == (5,)
        assert np.all(M > 0)

    def test_scales_linearly_in_columns_and_inversely_in_b(self):
        rng = np.random.default_rng(2)
        G = rng.standard_normal((40, 4))
        b = rng.standard_normal(40)
        base = wf.compute_m_scale(G, b)
        np.testing.assert_allclose(wf.compute_m_scale(3.0 * G, b), 3.0 * base)
        np.testing.assert_allclose(wf.compute_m_scale(G, 2.0 * b), base / 2.0)

    def test_zero_rhs_does_not_divide_by_zero(self):
        G = np.ones((10, 3))
        M = wf.compute_m_scale(G, np.zeros(10))
        assert np.all(np.isfinite(M))

    def test_on_a_real_ndim4_weak_library(self):
        """Smoke check on an actual 4-D weak library rather than random data."""
        from equivariant_pde_find.library import build_multigraph_library

        derivatives.set_backend("spectral")
        ndim, N, L, m, p = 4, 12, 2 * math.pi, 3, 4
        rng = np.random.default_rng(0)
        f = 0.3 * rng.standard_normal((N,) * ndim)
        ft = 0.3 * rng.standard_normal((N,) * ndim)

        G, desc = build_multigraph_library(
            f, L, L, formulation="weak", m=m, p=p, use_ibp=True,
        )
        b = wf.build_weak_rhs(ft, m=m, p=p)
        M = wf.compute_m_scale(G, b)
        assert M.shape == (len(desc),)
        assert np.all(np.isfinite(M))
