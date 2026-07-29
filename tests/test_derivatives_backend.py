"""
Tests for ``derivatives/__init__.py`` — backend switching and dispatch.

Covers:
- set_backend / get_backend
- use_backend context manager (including exception safety)
- Invalid backend error handling
- All six backend-neutral functions dispatch correctly
- After switching backends, the same function calls different implementations
"""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives


def _make_field(N=32, L=2 * math.pi):
    x = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    return np.cos(2 * X + 3 * Y), L


# ======================================================================
# set_backend / get_backend
# ======================================================================


class TestBackendSwitching:
    def test_default_is_spectral(self):
        # conftest.py sets spectral before each test
        assert derivatives.get_backend() == "spectral"

    def test_set_to_polynomial(self):
        derivatives.set_backend("polynomial")
        assert derivatives.get_backend() == "polynomial"
        derivatives.set_backend("spectral")

    def test_set_to_spectral(self):
        derivatives.set_backend("polynomial")
        derivatives.set_backend("spectral")
        assert derivatives.get_backend() == "spectral"

    def test_invalid_backend_raises(self):
        with pytest.raises(ValueError, match="Unknown backend"):
            derivatives.set_backend("finite_difference")


# ======================================================================
# use_backend context manager
# ======================================================================


class TestUseBackend:
    def test_temporary_switch(self):
        assert derivatives.get_backend() == "spectral"
        with derivatives.use_backend("polynomial"):
            assert derivatives.get_backend() == "polynomial"
        assert derivatives.get_backend() == "spectral"

    def test_nested_contexts(self):
        with derivatives.use_backend("polynomial"):
            assert derivatives.get_backend() == "polynomial"
            with derivatives.use_backend("spectral"):
                assert derivatives.get_backend() == "spectral"
            assert derivatives.get_backend() == "polynomial"
        assert derivatives.get_backend() == "spectral"

    def test_exception_restores_backend(self):
        assert derivatives.get_backend() == "spectral"
        try:
            with derivatives.use_backend("polynomial"):
                raise RuntimeError("test")
        except RuntimeError:
            pass
        assert derivatives.get_backend() == "spectral"


# ======================================================================
# All six backend-neutral functions exist and are callable
# ======================================================================


class TestDispatchFunctions:
    def test_grad_exists(self):
        f, L = _make_field()
        fx, fy = derivatives.grad(f, L, L)
        assert fx.shape == f.shape
        assert fy.shape == f.shape

    def test_laplacian_exists(self):
        f, L = _make_field()
        lap = derivatives.laplacian(f, L, L)
        assert lap.shape == f.shape

    def test_hessian_exists(self):
        f, L = _make_field()
        fxx, fxy, fyy = derivatives.hessian(f, L, L)
        assert fxx.shape == f.shape
        assert fxy.shape == f.shape
        assert fyy.shape == f.shape

    def test_bilaplacian_exists(self):
        f, L = _make_field()
        bilap = derivatives.bilaplacian(f, L, L)
        assert bilap.shape == f.shape

    def test_trilaplacian_exists(self):
        f, L = _make_field()
        trilap = derivatives.trilaplacian(f, L, L)
        assert trilap.shape == f.shape

    def test_third_order_exists(self):
        f, L = _make_field()
        fxxx, fxxy, fxyy, fyyy = derivatives.third_order(f, L, L)
        assert fxxx.shape == f.shape
        assert fxxy.shape == f.shape
        assert fxyy.shape == f.shape
        assert fyyy.shape == f.shape


# ======================================================================
# Backend-neutral functions give different results with different backends
# (on a smooth field, both should give similar results)
# ======================================================================


class TestBackendAgreement:
    def test_spectral_and_polynomial_agree_on_gradient(self):
        f, L = _make_field(N=64)
        with derivatives.use_backend("spectral"):
            fx_s, fy_s = derivatives.grad(f, L, L)
        with derivatives.use_backend("polynomial"):
            fx_p, fy_p = derivatives.grad(f, L, L)
        # Interior agreement should be very good for low-frequency mode
        r = 5
        np.testing.assert_allclose(
            fx_s[r:-r, r:-r], fx_p[r:-r, r:-r], rtol=1e-3, atol=1e-4,
        )
        np.testing.assert_allclose(
            fy_s[r:-r, r:-r], fy_p[r:-r, r:-r], rtol=1e-3, atol=1e-4,
        )

    def test_spectral_and_polynomial_agree_on_laplacian(self):
        f, L = _make_field(N=64)
        with derivatives.use_backend("spectral"):
            lap_s = derivatives.laplacian(f, L, L)
        with derivatives.use_backend("polynomial"):
            lap_p = derivatives.laplacian(f, L, L)
        r = 5
        np.testing.assert_allclose(
            lap_s[r:-r, r:-r], lap_p[r:-r, r:-r], rtol=5e-3, atol=1e-3,
        )


# ======================================================================
# prefilter_sigma kwarg passthrough
# ======================================================================


class TestPrefilterPassthrough:
    def test_prefilter_sigma_is_accepted(self):
        f, L = _make_field()
        # Should not raise
        lap = derivatives.laplacian(f, L, L, prefilter_sigma=0.5)
        assert lap.shape == f.shape

    def test_prefilter_sigma_none(self):
        f, L = _make_field()
        lap_no = derivatives.laplacian(f, L, L, prefilter_sigma=None)
        lap_0 = derivatives.laplacian(f, L, L, prefilter_sigma=0.0)
        np.testing.assert_allclose(lap_no, lap_0, rtol=1e-12, atol=1e-12)


# ======================================================================
# set_polynomial_stencil / get_polynomial_stencil
# ======================================================================


class TestPolynomialStencilConfig:
    def test_get_default_stencil(self):
        deg, r = derivatives.get_polynomial_stencil(2)
        assert isinstance(deg, int)
        assert isinstance(r, int)
        assert deg > 0
        assert r > 0

    def test_set_and_get_stencil(self):
        original = derivatives.get_polynomial_stencil(1)
        try:
            derivatives.set_polynomial_stencil(deg=7, r=3, order=1)
            deg, r = derivatives.get_polynomial_stencil(1)
            assert deg == 7
            assert r == 3
        finally:
            # Restore
            derivatives.set_polynomial_stencil(
                deg=original[0], r=original[1], order=1,
            )

    def test_clear_stencil_override(self):
        derivatives.set_polynomial_stencil(deg=None, r=None, order=1)
        deg, r = derivatives.get_polynomial_stencil(1)
        assert isinstance(deg, int)
        assert isinstance(r, int)


class TestAxisLengthsIsSingleSourced:
    """``axis_lengths`` used to exist as three independent copies (spectral.py,
    polynomial.py, derivatives/__init__.py), with nothing asserting they
    agreed -- a silent-divergence risk, since ``partial_nd`` used one copy
    while the named operators used the others.  They now all delegate to
    ``derivatives/_common.py``; these tests pin that."""

    def test_all_three_entry_points_are_the_same_function(self):
        from equivariant_pde_find.derivatives import polynomial as poly
        from equivariant_pde_find.derivatives import spectral as spec

        assert derivatives.axis_lengths is spec._axis_lengths
        assert derivatives.axis_lengths is poly._axis_lengths

    @pytest.mark.parametrize("ndim", [1, 2, 3, 4, 5, 6])
    def test_convention_by_dimension(self, ndim):
        f = np.zeros((3,) * ndim)
        Lx, Ly, Lz = 1.0, 2.0, 3.0
        got = derivatives.axis_lengths(f, Lx, Ly, Lz)
        assert len(got) == ndim
        if ndim == 1:
            assert got == (Lx,)
        elif ndim == 2:
            assert got == (Lx, Ly)  # anisotropic 2-D preserved
        elif ndim == 3:
            assert got == (Lx, Ly, Lz)  # anisotropic 3-D preserved
        else:
            assert got == (Lx,) * ndim  # isotropic beyond 3-D

    def test_missing_lz_in_3d_raises_a_clear_error(self):
        """Previously this surfaced as a TypeError deep inside the FFT
        wavenumber construction ('unsupported operand type(s) for /: float
        and NoneType')."""
        f = np.zeros((4, 4, 4))
        with pytest.raises(ValueError, match="Lz is required for 3-D fields"):
            derivatives.axis_lengths(f, 1.0, 2.0)

    def test_lz_is_ignored_beyond_3d(self):
        """Beyond 3-D the domain is isotropic, so a stray Lz must not leak in."""
        f = np.zeros((3,) * 5)
        assert derivatives.axis_lengths(f, 1.0, 2.0, 99.0) == (1.0,) * 5
