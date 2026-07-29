"""
Tests for ``library/conversions.py`` — symbolic Leibniz expansion (via
jet_algebra.JetExpr), WSINDy-to-jet mapping, and the multigraph jet
decomposition dictionaries (now generated from jet_algebra formulas rather
than hand-written, for any ndim).

Covers:
- The decomposition dictionaries (completeness, linear independence).
- JetExpr-based Leibniz expansion for simple cases (replaces the old
  hand-rolled _expand_fq_deriv/_terms_to_jet_dict, which no longer exist).
- build_wsindy_to_jet_map correctness.
- map_standard_to_jet, map_multigraph_to_jet, map_mod_multigraph_to_jet.
- Known-equation mapping tests (SH, AC, KS, CH, PME).
"""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.library.conversions import (
    _MULTIGRAPH_TO_JET_DECOMPOSITION,
    _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION,
    build_wsindy_to_jet_map,
    map_multigraph_to_jet,
    map_mod_multigraph_to_jet,
    map_standard_to_jet,
    map_wsindy_to_jet,
)
from equivariant_pde_find.jet_algebra import JetExpr, expr_to_dict
from equivariant_pde_find.library import (
    build_multigraph_library,
    build_mod_multigraph_library,
    build_pdefind_library,
    build_polyjet_library,
    build_wsindy_library,
)
from equivariant_pde_find.equations import (
    true_allen_cahn,
    true_cahn_hilliard,
    true_ks,
    true_pme,
    true_swift_hohenberg,
    true_allen_cahn_mod_multigraph,
    true_cahn_hilliard_mod_multigraph,
    true_ks_mod_multigraph,
    true_swift_hohenberg_mod_multigraph,
    true_pme_mod_multigraph,
    true_pme_wsindy,
)


# ======================================================================
# Decomposition completeness
# ======================================================================


def test_multigraph_decomposition_has_all_20_operators():
    """All 20 multigraph operators must have a jet decomposition."""
    from equivariant_pde_find.multigraph_operators import get_operators

    ops = get_operators(2, 3)
    op_names = {desc for _, desc in ops}
    assert op_names == set(_MULTIGRAPH_TO_JET_DECOMPOSITION.keys()), (
        f"Missing: {op_names - set(_MULTIGRAPH_TO_JET_DECOMPOSITION)}; "
        f"Extra: {set(_MULTIGRAPH_TO_JET_DECOMPOSITION) - op_names}"
    )


def test_mod_multigraph_decomposition_has_all_20_operators():
    """All 20 modified multigraph operators must have a jet decomposition."""
    from equivariant_pde_find.multigraph_operators import get_mod_multigraph_operators

    ops = get_mod_multigraph_operators(2, 3)
    op_names = {desc for _, desc in ops}
    assert op_names == set(_MOD_MULTIGRAPH_TO_JET_DECOMPOSITION.keys()), (
        f"Missing: {op_names - set(_MOD_MULTIGRAPH_TO_JET_DECOMPOSITION)}; "
        f"Extra: {set(_MOD_MULTIGRAPH_TO_JET_DECOMPOSITION) - op_names}"
    )


# ======================================================================
# JetExpr-based Leibniz expansion (replaces the old hand-rolled
# _expand_fq_deriv/_terms_to_jet_dict; same mathematical assertions,
# expressed via .deriv() + expr_to_dict())
# ======================================================================


def _expand(q, *axes):
    """∂_axes(f^q) as a {jet_description: coeff} dict, ndim=2."""
    expr = JetExpr.var(2) ** q
    for axis in axes:
        expr = expr.deriv(axis)
    return expr_to_dict(expr)


class TestJetExprLeibnizExpansion:
    def test_f1_no_derivs(self):
        """∂^0(f^1) = f."""
        assert _expand(1) == {"f": 1.0}

    def test_f1_dx(self):
        """∂_x(f) = f_x."""
        assert _expand(1, 0) == {"f_x": 1.0}

    def test_f2_dx(self):
        """∂_x(f²) = 2f·f_x."""
        result = _expand(2, 0)
        assert result == {"f·f_x": pytest.approx(2.0)}

    def test_f2_dxx(self):
        """∂_xx(f²) = 2f·f_xx + 2(f_x)²."""
        result = _expand(2, 0, 0)
        assert result["f·f_xx"] == pytest.approx(2.0)
        assert result["f_x²"] == pytest.approx(2.0)

    def test_f3_dxx(self):
        """∂_xx(f³) = 6f·(f_x)² + 3f²·f_xx."""
        result = _expand(3, 0, 0)
        assert any(abs(c - 3.0) < 1e-10 for c in result.values())
        assert any(abs(c - 6.0) < 1e-10 for c in result.values())

    def test_f3_dx_dy(self):
        """∂_x ∂_y(f³) = 6f·f_x·f_y + 3f²·f_xy."""
        result = _expand(3, 0, 1)
        assert len(result) == 2

    def test_f2_dx_dy(self):
        """∂_x ∂_y(f²) = 2f·f_xy + 2f_x·f_y."""
        result = _expand(2, 0, 1)
        assert any(abs(c - 2.0) < 1e-10 for c in result.values())


# ======================================================================
# expr_to_dict rendering (replaces the old _terms_to_jet_dict unit tests)
# ======================================================================


class TestExprToDictRendering:
    def test_single_f(self):
        assert expr_to_dict(JetExpr.var(2)) == {"f": 1.0}

    def test_fx(self):
        assert expr_to_dict(JetExpr.var(2).deriv(0)) == {"f_x": 1.0}

    def test_fx_squared(self):
        fx = JetExpr.var(2).deriv(0)
        assert expr_to_dict(2.0 * fx * fx) == {"f_x²": 2.0}

    def test_f_times_fx(self):
        f = JetExpr.var(2)
        fx = f.deriv(0)
        assert expr_to_dict(3.0 * f * fx) == {"f·f_x": 3.0}


# ======================================================================
# WSINDy-to-jet mapping (build_wsindy_to_jet_map)
# ======================================================================


def test_wsindy_to_jet_map_pme():
    """PME in WSINDy: (1/3)(f³)_xx + (1/3)(f³)_yy maps correctly to jet."""
    N, L = 32, 2 * math.pi
    rng = np.random.default_rng(42)
    f = rng.standard_normal((N, N))
    _, desc_ws = build_wsindy_library(f, L, L, max_order=2, max_poly=3)
    _, desc_jet = build_polyjet_library(f, L, L, max_order=4, max_degree=3)

    w_ws = true_pme_wsindy(desc_ws)
    M = build_wsindy_to_jet_map(desc_ws, desc_jet)
    w_jet_from_map = M @ w_ws

    # Manually verify: PME = f²·Δf + 2f·‖∇f‖²
    # = f²·f_xx + f²·f_yy + 2f·f_x² + 2f·f_y²
    expected = {}
    for name in ["f²·f_xx", "f²·f_yy", "f·f_x²", "f·f_y²"]:
        expected[name] = 0.0
    expected["f²·f_xx"] = 1.0
    expected["f²·f_yy"] = 1.0
    expected["f·f_x²"] = 2.0
    expected["f·f_y²"] = 2.0

    for name, val in expected.items():
        j = desc_jet.index(name)
        assert abs(w_jet_from_map[j] - val) < 1e-10, (
            f"Jet {name}: expected {val}, got {w_jet_from_map[j]}"
        )

    # Check no extra nonzero terms
    nz = np.sum(np.abs(w_jet_from_map) > 1e-10)
    assert nz == 4, f"Expected 4 nonzero jet terms, got {nz}"


def test_wsindy_to_jet_map_linear_terms():
    """Linear WSINDy terms f_{x}, f_{xx} map 1:1 to jet."""
    N, L = 32, 2 * math.pi
    rng = np.random.default_rng(99)
    f = rng.standard_normal((N, N))
    _, desc_ws = build_wsindy_library(f, L, L, max_order=2, max_poly=1)
    _, desc_jet = build_polyjet_library(f, L, L, max_order=2, max_degree=2)

    M = build_wsindy_to_jet_map(desc_ws, desc_jet)

    # WSINDy "f_{x}" should map to jet "f_x" with coefficient 1
    idx_ws = desc_ws.index("f_{x}")
    idx_jet = desc_jet.index("f_x")
    assert M[idx_jet, idx_ws] == pytest.approx(1.0)

    # WSINDy "f_{xx}" → jet "f_xx"
    idx_ws = desc_ws.index("f_{xx}")
    idx_jet = desc_jet.index("f_xx")
    assert M[idx_jet, idx_ws] == pytest.approx(1.0)


# ======================================================================
# map_standard_to_jet
# ======================================================================


def test_standard_to_jet_allen_cahn():
    """Standard PDE-FIND library for AC maps to jet with correct coefficients."""
    N, L = 32, 2 * math.pi
    rng = np.random.default_rng(7)
    f = rng.standard_normal((N, N))
    _, desc_std = build_pdefind_library(f, L, L, max_order=2, max_poly=3)
    _, desc_jet = build_polyjet_library(f, L, L, max_order=4, max_degree=3)

    # AC in standard basis: f - f³ + ε²·(f_xx + f_yy)
    eps = 0.15
    w_std = np.zeros(len(desc_std))
    w_std[desc_std.index("f")] = 1.0
    w_std[desc_std.index("f³")] = -1.0
    w_std[desc_std.index("f_xx")] = eps ** 2
    w_std[desc_std.index("f_yy")] = eps ** 2

    w_jet = map_standard_to_jet(w_std, desc_std, desc_jet)
    assert w_jet[desc_jet.index("f")] == pytest.approx(1.0)
    assert w_jet[desc_jet.index("f³")] == pytest.approx(-1.0)
    assert w_jet[desc_jet.index("f_xx")] == pytest.approx(eps ** 2)
    assert w_jet[desc_jet.index("f_yy")] == pytest.approx(eps ** 2)


# ======================================================================
# map_multigraph_to_jet — ground truth consistency
# ======================================================================


@pytest.mark.parametrize("name, true_fn", [
    ("AC", lambda d: true_allen_cahn(d, eps=0.15)),
    ("SH", lambda d: true_swift_hohenberg(d, r=0.3)),
    ("KS", true_ks),
    ("CH", lambda d: true_cahn_hilliard(d, eps=0.2)),
    ("PME", true_pme),
])
def test_multigraph_to_jet_ground_truth(name, true_fn):
    """Ground-truth multigraph coefficients mapped to jet must produce
    the same PDE prediction as the multigraph library."""
    N, L = 32, 2 * math.pi
    rng = np.random.default_rng(42)
    f = rng.standard_normal((N, N))

    T_mg, desc_mg = build_multigraph_library(f, L, L, max_p=2, max_n=3)
    T_jet, desc_jet = build_polyjet_library(f, L, L, max_order=4, max_degree=3)

    w_mg = true_fn(desc_mg)
    w_jet = map_multigraph_to_jet(w_mg, desc_mg, desc_jet)

    pred_mg = T_mg @ w_mg
    pred_jet = T_jet @ w_jet
    scale = max(np.max(np.abs(pred_mg)), 1e-10)
    rel_err = np.max(np.abs(pred_jet - pred_mg)) / scale
    assert rel_err < 1e-10, (
        f"{name}: multigraph→jet mapping error {rel_err:.2e}"
    )


# ======================================================================
# map_mod_multigraph_to_jet — ground truth consistency
# ======================================================================


@pytest.mark.parametrize("name, true_fn", [
    ("AC", lambda d: true_allen_cahn_mod_multigraph(d, eps=0.15)),
    ("SH", lambda d: true_swift_hohenberg_mod_multigraph(d, r=0.3)),
    ("KS", true_ks_mod_multigraph),
    ("CH", lambda d: true_cahn_hilliard_mod_multigraph(d, eps=0.2)),
    ("PME", true_pme_mod_multigraph),
])
def test_mod_multigraph_to_jet_ground_truth(name, true_fn):
    """Ground-truth mod. multigraph coefficients mapped to jet must produce
    the same PDE prediction as the mod. multigraph library."""
    N, L = 32, 2 * math.pi
    rng = np.random.default_rng(99)
    f = rng.standard_normal((N, N))

    T_mod, desc_mod = build_mod_multigraph_library(f, L, L, max_p=2, max_n=3)
    T_jet, desc_jet = build_polyjet_library(f, L, L, max_order=4, max_degree=3)

    w_mod = true_fn(desc_mod)
    w_jet = map_mod_multigraph_to_jet(w_mod, desc_mod, desc_jet)

    pred_mod = T_mod @ w_mod
    pred_jet = T_jet @ w_jet
    scale = max(np.max(np.abs(pred_mod)), 1e-10)
    rel_err = np.max(np.abs(pred_jet - pred_mod)) / scale
    assert rel_err < 1e-10, (
        f"{name}: mod. multigraph→jet mapping error {rel_err:.2e}"
    )


# ======================================================================
# Cross-library equivalence via jet space
# ======================================================================


def test_allen_cahn_all_bases_agree():
    """Allen-Cahn ground truth in all libraries must produce the same jet
    coefficient vector."""
    N, L = 32, 2 * math.pi
    rng = np.random.default_rng(123)
    f = rng.standard_normal((N, N))
    eps = 0.15

    # Build all libraries
    _, desc_mg = build_multigraph_library(f, L, L, max_p=2, max_n=3)
    _, desc_mod = build_mod_multigraph_library(f, L, L, max_p=2, max_n=3)
    _, desc_std = build_pdefind_library(f, L, L, max_order=4, max_poly=3)
    _, desc_ws = build_wsindy_library(f, L, L, max_order=4, max_poly=3)
    _, desc_jet = build_polyjet_library(f, L, L, max_order=4, max_degree=3)

    # Map all to jet
    w_jet_mg = map_multigraph_to_jet(
        true_allen_cahn(desc_mg, eps=eps), desc_mg, desc_jet,
    )
    w_jet_mod = map_mod_multigraph_to_jet(
        true_allen_cahn_mod_multigraph(desc_mod, eps=eps), desc_mod, desc_jet,
    )

    # Standard: f - f³ + ε²·(f_xx + f_yy)
    w_std = np.zeros(len(desc_std))
    w_std[desc_std.index("f")] = 1.0
    w_std[desc_std.index("f³")] = -1.0
    w_std[desc_std.index("f_xx")] = eps ** 2
    w_std[desc_std.index("f_yy")] = eps ** 2
    w_jet_std = map_standard_to_jet(w_std, desc_std, desc_jet)

    # All must agree
    np.testing.assert_allclose(w_jet_mg, w_jet_mod, atol=1e-12)
    np.testing.assert_allclose(w_jet_mg, w_jet_std, atol=1e-12)
