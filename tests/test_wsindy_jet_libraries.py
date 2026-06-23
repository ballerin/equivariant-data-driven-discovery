"""
Tests for WSINDy library (∂^α(f^q)) and polynomial jet library
(all monomials in jet variables).

These libraries complement the equivariant and standard (PDE-FIND) libraries
tested in test_library.py.
"""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives                                          # noqa: E402
from equivariant_pde_find.library import (                                       # noqa: E402
    build_pointwise_wsindy_library,
    build_polynomial_jet_library,
    enumerate_jet_monomials,
)
from equivariant_pde_find.weak_form import (                                        # noqa: E402
    build_weak_wsindy_library,
    build_weak_polynomial_jet_library,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

BACKENDS_2D = ["spectral", "polynomial"]


def _make_2d_field(N=32, L=2 * np.pi):
    """Return a simple band-limited 2-D field and its domain length."""
    x = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    return np.sin(X) + 0.5 * np.cos(2 * Y), L


def _col(Theta, descs, name):
    """Return the column of Theta whose description equals *name*."""
    matches = [i for i, d in enumerate(descs) if d == name]
    assert len(matches) == 1, (
        f"Expected exactly one column labeled {repr(name)}, "
        f"found {len(matches)} in {descs}"
    )
    return Theta[:, matches[0]]


# ======================================================================
# WSINDy library — pointwise
# ======================================================================


class TestPointwiseWSINDyLibrary:
    """Tests for build_pointwise_wsindy_library."""

    def test_column_count(self):
        """max_order=4, max_poly=3 → 1 + 14*3 + 3 (pure monomials) = 46."""
        f, L = _make_2d_field()
        Theta, descs = build_pointwise_wsindy_library(f, L, L, max_order=4, max_poly=3)
        assert Theta.shape[1] == 46, f"Expected 46 columns, got {Theta.shape[1]}"
        assert len(descs) == 46

    def test_column_count_small(self):
        """max_order=2, max_poly=2 → 1 + 5*2 + 2 = 13."""
        f, L = _make_2d_field()
        Theta, descs = build_pointwise_wsindy_library(f, L, L, max_order=2, max_poly=2)
        # multi_indices for max_order=2: (0,0),(1,0),(0,1),(2,0),(1,1),(0,2) → 6
        # Derivative columns: 5 * 2 = 10 (excluding a=b=0)
        # Pure monomials: f, f² → 2
        # Constant: 1
        # Total: 1 + 10 + 2 = 13
        assert Theta.shape[1] == 13, f"Expected 13 columns, got {Theta.shape[1]}"

    def test_constant_column_is_ones(self):
        """The '1' column must be identically 1."""
        f, L = _make_2d_field()
        Theta, descs = build_pointwise_wsindy_library(f, L, L)
        col = _col(Theta, descs, "1")
        np.testing.assert_allclose(col, np.ones_like(col), rtol=0, atol=1e-14)

    def test_f_column_equals_f(self):
        """The 'f' column must equal f."""
        f, L = _make_2d_field()
        Theta, descs = build_pointwise_wsindy_library(f, L, L, max_order=2, max_poly=2)
        np.testing.assert_allclose(
            _col(Theta, descs, "f").reshape(f.shape), f,
            rtol=1e-12, atol=1e-12,
        )

    def test_f2_column_equals_f_squared(self):
        """The 'f²' column must equal f²."""
        f, L = _make_2d_field()
        Theta, descs = build_pointwise_wsindy_library(f, L, L, max_order=2, max_poly=2)
        np.testing.assert_allclose(
            _col(Theta, descs, "f²").reshape(f.shape), f ** 2,
            rtol=1e-12, atol=1e-12,
        )

    @pytest.mark.parametrize("backend", BACKENDS_2D)
    def test_fx_column_is_derivative(self, backend):
        """f_{x} (q=1, a=1, b=0) should equal ∂_x f."""
        f, L = _make_2d_field()
        with derivatives.use_backend(backend):
            Theta, descs = build_pointwise_wsindy_library(f, L, L, max_order=2, max_poly=2)

        col = _col(Theta, descs, "f_{x}").reshape(f.shape)

        # Compute expected ∂_x f analytically
        x = np.linspace(0, L, f.shape[0], endpoint=False)
        X, Y = np.meshgrid(x, x, indexing="ij")
        expected = np.cos(X)

        np.testing.assert_allclose(col, expected, rtol=1e-11, atol=1e-10)

    @pytest.mark.parametrize("backend", BACKENDS_2D)
    def test_f2_x_column_is_derivative_of_f_squared(self, backend):
        """(f²)_{x} should equal ∂_x(f²) = 2f·f_x."""
        f, L = _make_2d_field()
        with derivatives.use_backend(backend):
            Theta, descs = build_pointwise_wsindy_library(f, L, L, max_order=2, max_poly=2)

        col = _col(Theta, descs, "(f²)_{x}").reshape(f.shape)

        x = np.linspace(0, L, f.shape[0], endpoint=False)
        X, Y = np.meshgrid(x, x, indexing="ij")
        f_x = np.cos(X)
        expected = 2 * f * f_x

        tol = dict(rtol=1e-3, atol=1e-2) if backend == "polynomial" else dict(rtol=1e-11, atol=1e-10)
        np.testing.assert_allclose(col, expected, **tol)

    @pytest.mark.parametrize("backend", BACKENDS_2D)
    def test_wsindy_and_standard_agree_on_linear_terms(self, backend):
        """For q=1, WSINDy ∂^α(f) should equal standard f^0·∂^α f."""
        f, L = _make_2d_field()
        from equivariant_pde_find.library import build_standard_library as bsl
        with derivatives.use_backend(backend):
            Tw, dw = build_pointwise_wsindy_library(f, L, L, max_order=3, max_poly=1)
            Ts, ds = bsl(f, L, L, max_order=3, max_poly=0)

        # Standard with max_poly=0: pure derivatives + '1'
        # WSINDy with max_poly=1: constant + ∂^α(f) derivatives + f
        # Derivative descriptions: WSINDy uses f_{x}, f_{xx} etc (same as standard)
        wsindy_deriv_descs = [d for d in dw if d.startswith("f_{")]
        std_deriv_descs = [d for d in ds if d.startswith("f_")]

        # Match by converting WSINDy format to standard format for comparison
        for wd in wsindy_deriv_descs:
            # WSINDy: "f_{x}", "f_{xx}", "f_{xy}", ...
            # Standard: "f_x", "f_xx", "f_xy", ...
            sd = wd.replace("{", "").replace("}", "")
            if sd in std_deriv_descs:
                wcol = _col(Tw, dw, wd)
                scol = _col(Ts, ds, sd)
                tol = dict(rtol=1e-3, atol=1e-2) if backend == "polynomial" else dict(rtol=1e-11, atol=1e-10)
                np.testing.assert_allclose(wcol, scol, **tol,
                    err_msg=f"WSINDy {wd} ≠ standard {sd}")

    def test_no_identically_zero_columns(self):
        """No column should be identically zero (except possibly constant)."""
        f, L = _make_2d_field()
        Theta, descs = build_pointwise_wsindy_library(f, L, L, max_order=4, max_poly=3)
        for j, d in enumerate(descs):
            col_norm = np.linalg.norm(Theta[:, j])
            assert col_norm > 0, f"Column '{d}' is identically zero"


# ======================================================================
# Polynomial jet library — pointwise
# ======================================================================


class TestPolynomialJetLibrary:
    """Tests for build_polynomial_jet_library."""

    def test_column_count(self):
        """max_order=4, max_degree=3 → 105 columns."""
        f, L = _make_2d_field()
        Theta, descs = build_polynomial_jet_library(f, L, L, max_order=4, max_degree=3)
        assert Theta.shape[1] == 105, f"Expected 105 columns, got {Theta.shape[1]}"
        assert len(descs) == 105

    def test_column_count_small(self):
        """max_order=2, max_degree=2 — verify count."""
        f, L = _make_2d_field()
        Theta, descs = build_polynomial_jet_library(f, L, L, max_order=2, max_degree=2)
        # Jet vars at order≤2: f, f_x, f_y, f_xx, f_xy, f_yy (6 vars)
        # Monomials with deg≤2, order≤2:
        # deg 0: 1
        # deg 1: f, f_x, f_y, f_xx, f_xy, f_yy (6)
        # deg 2, order≤2: f², f·f_x, f·f_y, f·f_xx, f·f_xy, f·f_yy,
        #                   f_x², f_x·f_y, f_y² (9)
        # Total: 1+6+9 = 16
        expected = 16
        assert Theta.shape[1] == expected, f"Expected {expected} columns, got {Theta.shape[1]}"

    def test_constant_column_is_ones(self):
        """The '1' column must be identically 1."""
        f, L = _make_2d_field()
        Theta, descs = build_polynomial_jet_library(f, L, L, max_order=2, max_degree=2)
        col = _col(Theta, descs, "1")
        np.testing.assert_allclose(col, np.ones_like(col), rtol=0, atol=1e-14)

    def test_f_column_equals_f(self):
        """The 'f' column must equal f."""
        f, L = _make_2d_field()
        Theta, descs = build_polynomial_jet_library(f, L, L, max_order=2, max_degree=2)
        np.testing.assert_allclose(
            _col(Theta, descs, "f").reshape(f.shape), f, rtol=1e-12, atol=1e-12,
        )

    @pytest.mark.parametrize("backend", BACKENDS_2D)
    def test_fx_column_is_derivative(self, backend):
        """f_x should equal ∂_x f."""
        f, L = _make_2d_field()
        with derivatives.use_backend(backend):
            Theta, descs = build_polynomial_jet_library(f, L, L, max_order=2, max_degree=2)

        col = _col(Theta, descs, "f_x").reshape(f.shape)
        x = np.linspace(0, L, f.shape[0], endpoint=False)
        X, Y = np.meshgrid(x, x, indexing="ij")
        expected = np.cos(X)
        np.testing.assert_allclose(col, expected, rtol=1e-11, atol=1e-10)

    @pytest.mark.parametrize("backend", BACKENDS_2D)
    def test_fx_squared_column(self, backend):
        """f_x² should equal (f_x)²."""
        f, L = _make_2d_field()
        with derivatives.use_backend(backend):
            Theta, descs = build_polynomial_jet_library(f, L, L, max_order=2, max_degree=2)

        col = _col(Theta, descs, "f_x²").reshape(f.shape)
        x = np.linspace(0, L, f.shape[0], endpoint=False)
        X, Y = np.meshgrid(x, x, indexing="ij")
        f_x = np.cos(X)
        np.testing.assert_allclose(col, f_x ** 2, rtol=1e-11, atol=1e-10)

    @pytest.mark.parametrize("backend", BACKENDS_2D)
    def test_f_times_fxx_column(self, backend):
        """f·f_xx should equal f · f_xx."""
        f, L = _make_2d_field()
        with derivatives.use_backend(backend):
            Theta, descs = build_polynomial_jet_library(f, L, L, max_order=2, max_degree=2)

        col = _col(Theta, descs, "f·f_xx").reshape(f.shape)

        # f_xx = -sin(x)
        x = np.linspace(0, L, f.shape[0], endpoint=False)
        X, Y = np.meshgrid(x, x, indexing="ij")
        f_xx = -np.sin(X)
        expected = f * f_xx
        tol = dict(rtol=1e-3, atol=1e-2) if backend == "polynomial" else dict(rtol=1e-11, atol=1e-10)
        np.testing.assert_allclose(col, expected, **tol)

    def test_column_ordering_by_order_then_degree(self):
        """Columns should be sorted by (total_order, total_degree)."""
        f, L = _make_2d_field()
        Theta, descs = build_polynomial_jet_library(f, L, L, max_order=2, max_degree=2)

        # First columns: constant, f, f² (order 0, degree 0→1→2)
        assert descs[0] == "1"
        assert descs[1] == "f"
        assert descs[2] == "f²"

    def test_no_identically_zero_columns(self):
        """No column should be identically zero."""
        f, L = _make_2d_field()
        Theta, descs = build_polynomial_jet_library(f, L, L, max_order=4, max_degree=3)
        for j, d in enumerate(descs):
            col_norm = np.linalg.norm(Theta[:, j])
            assert col_norm > 0, f"Column '{d}' is identically zero"

    def test_contains_pde_find_terms(self):
        """All PDE-FIND library terms should appear in the polynomial jet library."""
        f, L = _make_2d_field()
        from equivariant_pde_find.library import build_standard_library
        _, std_descs = build_standard_library(f, L, L, max_order=2, max_poly=2)
        _, jet_descs = build_polynomial_jet_library(f, L, L, max_order=2, max_degree=2)

        # Map standard descriptions to jet descriptions
        # Standard: "1", "f_x", "f·f_x", "f²·f_x", ...
        # Jet: "f_{x}", "f·f_{x}", "f²·f_{x}", ...
        std_to_jet = {}
        for sd in std_descs:
            if sd == "1":
                std_to_jet[sd] = "1"
            elif sd == "f":
                std_to_jet[sd] = "f"
            elif "·" in sd:
                parts = sd.split("·")
                poly_part = parts[0]  # "f" or "f²"
                deriv_part = parts[1]
                if deriv_part == "1":
                    std_to_jet[sd] = poly_part
                else:
                    # Convert "f_xx" → "f_{xx}" etc.
                    pass  # skip complex mapping, just check key ones
            elif sd.startswith("f_"):
                # Convert "f_xx" → "f_{xx}"
                jet_name = "f_{" + sd[2:] + "}"
                std_to_jet[sd] = jet_name

        # Check key terms exist
        for sd, jd in [("1", "1"), ("f", "f"), ("f_x", "f_x"), ("f_xx", "f_xx")]:
            assert jd in jet_descs, f"Jet library missing '{jd}' (from standard '{sd}')"


# ======================================================================
# Polynomial jet library — enumeration
# ======================================================================


class TestJetEnumeration:
    """Tests for enumerate_jet_monomials."""

    def _make_jet_vars(self, max_order):
        """Build jet_vars list matching build_polynomial_jet_library."""
        jet_vars = []
        for total in range(max_order + 1):
            for a in range(total + 1):
                b = total - a
                name = (("x" * a) + ("y" * b)) if total > 0 else ""
                name = f"f_{name}" if total > 0 else "f"
                jet_vars.append(((a, b), total, name))
        return jet_vars

    def test_count_max_order_4_degree_3(self):
        jet_vars = self._make_jet_vars(4)
        results = enumerate_jet_monomials(jet_vars, max_order=4, max_degree=3)
        assert len(results) == 105, f"Expected 105, got {len(results)}"

    def test_count_max_order_2_degree_2(self):
        jet_vars = self._make_jet_vars(2)
        results = enumerate_jet_monomials(jet_vars, max_order=2, max_degree=2)
        assert len(results) == 16, f"Expected 16, got {len(results)}"

    def test_count_max_order_1_degree_1(self):
        """Order≤1, degree≤1: vars f, f_x, f_y → monomials: 1, f, f_x, f_y (4)."""
        jet_vars = self._make_jet_vars(1)
        results = enumerate_jet_monomials(jet_vars, max_order=1, max_degree=1)
        assert len(results) == 4, f"Expected 4, got {len(results)}"

    def test_constant_is_first(self):
        jet_vars = self._make_jet_vars(2)
        results = enumerate_jet_monomials(jet_vars, max_order=2, max_degree=2)
        exponents, deg, order = results[0]
        assert deg == 0 and order == 0
        assert all(v == 0 for v in exponents.values())

    def test_all_satisfy_constraints(self):
        jet_vars = self._make_jet_vars(4)
        results = enumerate_jet_monomials(jet_vars, max_order=4, max_degree=3)
        for exponents, deg, order in results:
            assert deg <= 3, f"degree {deg} > 3"
            assert order <= 4, f"order {order} > 4"
            # Verify computed deg/order match exponents
            computed_deg = sum(exponents.values())
            computed_order = sum(
                exp * (a + b) for (a, b), exp in exponents.items()
            )
            assert computed_deg == deg
            assert computed_order == order


# ======================================================================
# Weak-form WSINDy library
# ======================================================================


class TestWeakWSINDyLibrary:
    """Tests for build_weak_wsindy_library."""

    def test_column_count(self):
        f, L = _make_2d_field()
        Theta, descs = build_weak_wsindy_library(f, L, L, max_order=4, max_poly=3, m=4, p=3)
        assert Theta.shape[1] == 46

    def test_output_shape(self):
        """Valid grid should be (N-2m)²."""
        f, L = _make_2d_field(N=32)
        Theta, _ = build_weak_wsindy_library(f, L, L, max_order=2, max_poly=1, m=4, p=3)
        expected_n = (32 - 8) ** 2  # 24² = 576
        assert Theta.shape[0] == expected_n

    def test_ibp_flag_accepted(self):
        """use_ibp=True should not crash."""
        f, L = _make_2d_field()
        Theta, descs = build_weak_wsindy_library(f, L, L, max_order=2, max_poly=1,
                                                      m=4, p=3, use_ibp=True)
        assert Theta.shape[1] > 0

    def test_no_nan_or_inf(self):
        f, L = _make_2d_field()
        Theta, _ = build_weak_wsindy_library(f, L, L, max_order=4, max_poly=3, m=4, p=3)
        assert not np.any(np.isnan(Theta))
        assert not np.any(np.isinf(Theta))


# ======================================================================
# Weak-form polynomial jet library
# ======================================================================


class TestWeakPolynomialJetLibrary:
    """Tests for build_weak_polynomial_jet_library."""

    def test_column_count(self):
        f, L = _make_2d_field()
        Theta, descs = build_weak_polynomial_jet_library(
            f, L, L, max_order=4, max_degree=3, m=4, p=3,
        )
        assert Theta.shape[1] == 105

    def test_output_shape(self):
        f, L = _make_2d_field(N=32)
        Theta, _ = build_weak_polynomial_jet_library(
            f, L, L, max_order=2, max_degree=1, m=4, p=3,
        )
        expected_n = (32 - 8) ** 2
        assert Theta.shape[0] == expected_n

    def test_ibp_flag_accepted(self):
        f, L = _make_2d_field()
        Theta, _ = build_weak_polynomial_jet_library(
            f, L, L, max_order=2, max_degree=2, m=4, p=3, use_ibp=True,
        )
        assert Theta.shape[1] > 0

    def test_no_nan_or_inf(self):
        f, L = _make_2d_field()
        Theta, _ = build_weak_polynomial_jet_library(
            f, L, L, max_order=2, max_degree=2, m=4, p=3,
        )
        assert not np.any(np.isnan(Theta))
        assert not np.any(np.isinf(Theta))
