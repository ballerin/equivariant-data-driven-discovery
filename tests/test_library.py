"""
Tests for `library.py` — build_standard_library, build_standard_library_1d,
and build_equivariant_library.

The most critical class of bugs to guard against is a *polynomial-degree shift*
in the standard library: a column labeled "f^p" must contain f^p multiplied by
the derivative operator, not f^{p+1}.  This happened before because the
(a=0, b=0) derivative entry stored f instead of np.ones_like(f).  The 1D
variant is tested here too, even though it was correct, to prevent regressions.
"""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, CODE_DIR)

import derivatives                                          # noqa: E402
from library import (                                       # noqa: E402
    build_equivariant_library,
    build_standard_library,
    build_standard_library_1d,
)
from multigraph_operators import get_operators              # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_2d_field(N=32, L=2 * np.pi):
    """Return a simple band-limited 2-D field and its domain length."""
    x = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    return np.sin(X) + 0.5 * np.cos(2 * Y), L


def _make_1d_field(N=64, L=2 * np.pi):
    """Return a simple band-limited 1-D field and its domain length."""
    x = np.linspace(0, L, N, endpoint=False)
    return np.sin(x) + 0.3 * np.cos(2 * x), L


def _col(Theta, descs, name):
    """Return the column of Theta whose description equals *name*."""
    matches = [i for i, d in enumerate(descs) if d == name]
    assert len(matches) == 1, (
        f"Expected exactly one column labeled {repr(name)}, "
        f"found {len(matches)} in {descs}"
    )
    return Theta[:, matches[0]]


# ---------------------------------------------------------------------------
# 2-D standard library — polynomial degree correctness
# ---------------------------------------------------------------------------

BACKENDS_2D = ["spectral", "polynomial"]


@pytest.mark.parametrize("backend", BACKENDS_2D)
def test_standard_library_constant_column_is_ones(backend):
    """The '1' column must be identically 1 for both backends.

    Before the fix the (a=0,b=0) derivative entry stored f instead of
    np.ones_like(f), so the '1' column silently contained f and every
    polynomial column was shifted up by one degree.
    """
    f, L = _make_2d_field()
    with derivatives.use_backend(backend):
        Theta, descs = build_standard_library(f, L, L, max_order=2, max_poly=3)
    col = _col(Theta, descs, "1")
    np.testing.assert_allclose(
        col, np.ones_like(col),
        rtol=0, atol=1e-14,
        err_msg=(
            f"[backend={backend}] The '1' column should be all-ones; "
            "if it equals f, the degree-shift bug has reappeared."
        ),
    )


@pytest.mark.parametrize("backend", BACKENDS_2D)
def test_standard_library_f_column_equals_f(backend):
    """The 'f' column must equal f (not f²)."""
    f, L = _make_2d_field()
    with derivatives.use_backend(backend):
        Theta, descs = build_standard_library(f, L, L, max_order=2, max_poly=3)
    np.testing.assert_allclose(
        _col(Theta, descs, "f").reshape(f.shape), f,
        rtol=1e-12, atol=1e-12,
        err_msg=f"[backend={backend}] The 'f' column should equal f, not f².",
    )


@pytest.mark.parametrize("backend", BACKENDS_2D)
def test_standard_library_f2_column_equals_f_squared(backend):
    """The 'f²' column must equal f² (not f³)."""
    f, L = _make_2d_field()
    with derivatives.use_backend(backend):
        Theta, descs = build_standard_library(f, L, L, max_order=2, max_poly=3)
    np.testing.assert_allclose(
        _col(Theta, descs, "f²").reshape(f.shape), f ** 2,
        rtol=1e-12, atol=1e-12,
        err_msg=f"[backend={backend}] The 'f²' column should equal f², not f³.",
    )


@pytest.mark.parametrize("backend", BACKENDS_2D)
def test_standard_library_f3_column_equals_f_cubed(backend):
    """The 'f³' column must equal f³ (not f⁴)."""
    f, L = _make_2d_field()
    with derivatives.use_backend(backend):
        Theta, descs = build_standard_library(f, L, L, max_order=2, max_poly=3)
    np.testing.assert_allclose(
        _col(Theta, descs, "f³").reshape(f.shape), f ** 3,
        rtol=1e-12, atol=1e-12,
        err_msg=f"[backend={backend}] The 'f³' column should equal f³, not f⁴.",
    )


# ---------------------------------------------------------------------------
# 2-D standard library — product term correctness
# ---------------------------------------------------------------------------

def test_standard_library_product_term_equals_f_times_derivative():
    """A product column 'f·f_xx' must equal f multiplied by the f_xx column.

    This test catches off-by-one bugs in the polynomial factor of product terms:
    if the library has degree-shift, 'f·f_xx' would contain f²·f_xx instead.
    Uses the spectral backend for exact derivatives.
    """
    f, L = _make_2d_field(N=32)
    with derivatives.use_backend("spectral"):
        Theta, descs = build_standard_library(f, L, L, max_order=2, max_poly=2)

    f_xx_col   = _col(Theta, descs, "f_xx").reshape(f.shape)
    f_f_xx_col = _col(Theta, descs, "f·f_xx").reshape(f.shape)

    np.testing.assert_allclose(
        f_f_xx_col, f * f_xx_col,
        rtol=1e-12, atol=1e-12,
        err_msg=(
            "'f·f_xx' column should equal f * f_xx column element-wise; "
            "a degree-shift bug would make it equal f² * f_xx."
        ),
    )


def test_standard_library_product_terms_are_consistent():
    """For every product column 'f^p·D', verify it equals (f^p) * (D column).

    Covers all polynomial degrees and all derivative columns so that a degree
    shift in any factor is detected.
    """
    f, L = _make_2d_field(N=32)
    max_order, max_poly = 2, 3
    with derivatives.use_backend("spectral"):
        Theta, descs = build_standard_library(f, L, L, max_order=max_order, max_poly=max_poly)

    # Gather the pure-derivative columns (p=0 factor)
    pure_deriv_descs = [d for d in descs if not d.startswith("f^") and not d.startswith("f·")]
    # Remove the plain "f" which is p=1, d="1"
    pure_deriv_descs = [d for d in pure_deriv_descs if d != "f"]

    for d_name in pure_deriv_descs:
        D = _col(Theta, descs, d_name).reshape(f.shape)
        for p in range(1, max_poly + 1):
            if p == 1:
                prod_name = f"f·{d_name}" if d_name != "1" else "f"
            else:
                sup = {2: "²", 3: "³", 4: "⁴"}.get(p, f"^{p}")
                prod_name = f"f{sup}·{d_name}" if d_name != "1" else f"f{sup}"
            if prod_name not in descs:
                continue
            col = _col(Theta, descs, prod_name).reshape(f.shape)
            np.testing.assert_allclose(
                col, f ** p * D,
                rtol=1e-12, atol=1e-12,
                err_msg=(
                    f"Column {repr(prod_name)} should equal f^{p} * {repr(d_name)} column "
                    f"element-wise. Got max deviation "
                    f"{np.max(np.abs(col - f**p * D)):.3e}."
                ),
            )


# ---------------------------------------------------------------------------
# 2-D standard library — derivative column correctness (spectral only)
# ---------------------------------------------------------------------------

def test_standard_library_derivative_columns_spectral():
    """Derivative columns must match the analytic derivatives of a known field.

    Uses f(x,y) = sin(x) + 0.5*cos(2y) whose partials are analytic:
      f_x  = cos(x),  f_y  = -sin(2y),
      f_xx = -sin(x), f_yy = -2*cos(2y).
    Only tests the spectral backend (which is exact for band-limited data).
    """
    N = 64
    L = 2 * np.pi
    x = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    f = np.sin(X) + 0.5 * np.cos(2 * Y)

    analytic = {
        "f_x":  np.cos(X),
        "f_y":  -np.sin(2 * Y),
        "f_xx": -np.sin(X),
        "f_yy": -2 * np.cos(2 * Y),
    }

    with derivatives.use_backend("spectral"):
        Theta, descs = build_standard_library(f, L, L, max_order=2, max_poly=1)

    for name, expected in analytic.items():
        col = _col(Theta, descs, name).reshape(f.shape)
        np.testing.assert_allclose(
            col, expected,
            rtol=1e-10, atol=1e-10,
            err_msg=(
                f"Derivative column {repr(name)} does not match analytic value. "
                "This may indicate a backend wiring or index-ordering bug."
            ),
        )


# ---------------------------------------------------------------------------
# 2-D standard library — column count
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("max_order,max_poly", [(2, 2), (2, 3), (4, 3)])
def test_standard_library_column_count(max_order, max_poly):
    """Column count must equal (number of 2-D multi-indices) * (max_poly + 1).

    The number of 2-D multi-indices (a, b) with a+b <= max_order is
    (max_order+1)*(max_order+2)/2.
    """
    f, L = _make_2d_field()
    with derivatives.use_backend("spectral"):
        Theta, descs = build_standard_library(f, L, L, max_order=max_order, max_poly=max_poly)
    expected_n_derivs = (max_order + 1) * (max_order + 2) // 2
    expected_n_cols   = expected_n_derivs * (max_poly + 1)
    assert Theta.shape[1] == expected_n_cols, (
        f"max_order={max_order}, max_poly={max_poly}: "
        f"expected {expected_n_cols} columns ({expected_n_derivs} derivative terms "
        f"× {max_poly+1} polynomial factors), got {Theta.shape[1]}."
    )
    assert len(descs) == expected_n_cols, (
        f"len(descs)={len(descs)} != Theta.shape[1]={Theta.shape[1]}."
    )


# ---------------------------------------------------------------------------
# 1-D standard library — polynomial degree correctness
# ---------------------------------------------------------------------------

def test_standard_library_1d_constant_column_is_ones():
    """The '1' column of build_standard_library_1d must be identically 1."""
    u, L = _make_1d_field()
    with derivatives.use_backend("spectral"):
        Theta, descs = build_standard_library_1d(u, L, max_deriv=3, max_poly=3)
    col = _col(Theta, descs, "1")
    np.testing.assert_allclose(
        col, np.ones_like(col),
        rtol=0, atol=1e-14,
        err_msg="1D library: the '1' column should be identically 1.",
    )


def test_standard_library_1d_u_column_equals_u():
    """The 'u' column of build_standard_library_1d must equal u (not u²)."""
    u, L = _make_1d_field()
    with derivatives.use_backend("spectral"):
        Theta, descs = build_standard_library_1d(u, L, max_deriv=3, max_poly=3)
    np.testing.assert_allclose(
        _col(Theta, descs, "u"), u,
        rtol=1e-12, atol=1e-12,
        err_msg="1D library: the 'u' column should equal u, not u².",
    )


def test_standard_library_1d_u2_column_equals_u_squared():
    """The 'u²' column of build_standard_library_1d must equal u²."""
    u, L = _make_1d_field()
    with derivatives.use_backend("spectral"):
        Theta, descs = build_standard_library_1d(u, L, max_deriv=3, max_poly=3)
    np.testing.assert_allclose(
        _col(Theta, descs, "u²"), u ** 2,
        rtol=1e-12, atol=1e-12,
        err_msg="1D library: the 'u²' column should equal u², not u³.",
    )


def test_standard_library_1d_column_count():
    """Column count must equal (max_deriv+1) * (max_poly+1)."""
    u, L = _make_1d_field()
    max_deriv, max_poly = 3, 3
    with derivatives.use_backend("spectral"):
        Theta, descs = build_standard_library_1d(u, L, max_deriv=max_deriv, max_poly=max_poly)
    expected = (max_deriv + 1) * (max_poly + 1)
    assert Theta.shape[1] == expected, (
        f"Expected {expected} columns ((max_deriv+1)*(max_poly+1) = "
        f"{max_deriv+1}×{max_poly+1}), got {Theta.shape[1]}."
    )


# ---------------------------------------------------------------------------
# Equivariant library — structural checks
# ---------------------------------------------------------------------------

def test_equivariant_library_column_count():
    """build_equivariant_library column count must match get_operators count."""
    f, L = _make_2d_field()
    for max_p, max_n in [(2, 3), (1, 2), (2, 4)]:
        with derivatives.use_backend("spectral"):
            Theta, descs = build_equivariant_library(f, L, L, max_p=max_p, max_n=max_n)
        expected = len(get_operators(max_p, max_n))
        assert Theta.shape[1] == expected, (
            f"max_p={max_p}, max_n={max_n}: expected {expected} columns from "
            f"get_operators, got {Theta.shape[1]}."
        )
        assert len(descs) == expected, (
            f"len(descs)={len(descs)} != expected {expected}."
        )


def test_equivariant_library_descriptions_match_operators():
    """Descriptions returned by build_equivariant_library must match get_operators exactly."""
    f, L = _make_2d_field()
    with derivatives.use_backend("spectral"):
        _, descs = build_equivariant_library(f, L, L, max_p=2, max_n=3)
    expected_descs = [name for _, name in get_operators(2, 3)]
    assert descs == expected_descs, (
        "Equivariant library descriptions do not match get_operators order. Diff:\n"
        + "\n".join(
            f"  [{i}] expected {repr(e)}, got {repr(a)}"
            for i, (e, a) in enumerate(zip(expected_descs, descs))
            if e != a
        )
    )


def test_equivariant_library_constant_column_is_ones():
    """The '1' operator column must be identically 1 for equivariant library."""
    f, L = _make_2d_field()
    with derivatives.use_backend("spectral"):
        Theta, descs = build_equivariant_library(f, L, L, max_p=2, max_n=3)
    col = _col(Theta, descs, "1")
    np.testing.assert_allclose(
        col, np.ones_like(col),
        rtol=0, atol=1e-14,
        err_msg="Equivariant library: the '1' column should be identically 1.",
    )


def test_equivariant_library_f_column_equals_f():
    """The 'f' operator column must equal f for equivariant library."""
    f, L = _make_2d_field()
    with derivatives.use_backend("spectral"):
        Theta, descs = build_equivariant_library(f, L, L, max_p=2, max_n=3)
    np.testing.assert_allclose(
        _col(Theta, descs, "f").reshape(f.shape), f,
        rtol=1e-12, atol=1e-12,
        err_msg="Equivariant library: the 'f' column should equal f.",
    )
