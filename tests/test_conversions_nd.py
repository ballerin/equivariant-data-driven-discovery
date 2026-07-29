"""Tests for the N-D generalization of library/conversions.py.

Covers:
- The new jet_algebra-driven multigraph/mod-multigraph decompositions
  reproduce the historical hand-written 2-D dicts exactly (regression
  guard: if this ever fails, the symbolic formulas have diverged from the
  values the paper's 2-D results were computed with).
- Decompositions are well-defined and internally consistent at ndim=4,5.
- build_wsindy_to_jet_map works at ndim=4.
"""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.library.conversions import (
    _multigraph_to_jet_decomposition,
    _mod_multigraph_to_jet_decomposition,
    _MULTIGRAPH_TO_JET_DECOMPOSITION,
    _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION,
    build_wsindy_to_jet_map,
    map_wsindy_to_jet,
)
from equivariant_pde_find.library import build_wsindy_library, build_polyjet_library
from equivariant_pde_find.multigraph_operators import get_operators, get_mod_multigraph_operators


# Historical hand-written 2-D decomposition (exact values the paper's
# results were computed with), kept here ONLY as a frozen regression
# fixture -- do not "clean up" by importing from conversions.py, since the
# whole point is comparing the new generator's output against these
# independently-recorded values.
_HISTORICAL_COMMON = {
    "1": {"1": 1.0}, "f": {"f": 1.0}, "f²": {"f²": 1.0}, "f³": {"f³": 1.0},
    "Δf": {"f_xx": 1.0, "f_yy": 1.0},
    "‖∇f‖²": {"f_x²": 1.0, "f_y²": 1.0},
    "f·‖∇f‖²": {"f·f_x²": 1.0, "f·f_y²": 1.0},
    "Δ²f": {"f_xxxx": 1.0, "f_xxyy": 2.0, "f_yyyy": 1.0},
    "Δf·Δf": {"f_xx²": 1.0, "f_yy·f_xx": 2.0, "f_yy²": 1.0},
    "‖∇²f‖²": {"f_xx²": 1.0, "f_xy²": 2.0, "f_yy²": 1.0},
    "f·Δf·Δf": {"f·f_xx²": 1.0, "f·f_yy·f_xx": 2.0, "f·f_yy²": 1.0},
    "f·‖∇²f‖²": {"f·f_xx²": 1.0, "f·f_xy²": 2.0, "f·f_yy²": 1.0},
    "∇²f(∇f,∇f)": {"f_x²·f_xx": 1.0, "f_y·f_x·f_xy": 2.0, "f_y²·f_yy": 1.0},
    "‖∇f‖²·Δf": {"f_x²·f_xx": 1.0, "f_x²·f_yy": 1.0, "f_y²·f_xx": 1.0, "f_y²·f_yy": 1.0},
}
_HISTORICAL_MULTIGRAPH_ONLY = {
    "f·Δf": {"f·f_xx": 1.0, "f·f_yy": 1.0},
    "f²·Δf": {"f²·f_xx": 1.0, "f²·f_yy": 1.0},
    "f·Δ²f": {"f·f_xxxx": 1.0, "f·f_xxyy": 2.0, "f·f_yyyy": 1.0},
    "⟨∇Δf,∇f⟩": {"f_x·f_xxx": 1.0, "f_x·f_xyy": 1.0, "f_y·f_xxy": 1.0, "f_y·f_yyy": 1.0},
    "f²·Δ²f": {"f²·f_xxxx": 1.0, "f²·f_xxyy": 2.0, "f²·f_yyyy": 1.0},
    "f·⟨∇Δf,∇f⟩": {"f·f_x·f_xxx": 1.0, "f·f_x·f_xyy": 1.0, "f·f_y·f_xxy": 1.0, "f·f_y·f_yyy": 1.0},
}
_HISTORICAL_MOD_MULTIGRAPH_ONLY = {
    "Δf²": {"f·f_xx": 2.0, "f·f_yy": 2.0, "f_x²": 2.0, "f_y²": 2.0},
    "Δf³": {"f²·f_xx": 3.0, "f²·f_yy": 3.0, "f·f_x²": 6.0, "f·f_y²": 6.0},
    "Δ‖∇f‖²": {"f_x·f_xxx": 2.0, "f_x·f_xyy": 2.0, "f_y·f_xxy": 2.0, "f_y·f_yyy": 2.0,
              "f_xx²": 2.0, "f_xy²": 4.0, "f_yy²": 2.0},
    "Δ(f·‖∇f‖²)": {"f·f_x·f_xxx": 2.0, "f·f_x·f_xyy": 2.0, "f·f_y·f_xxy": 2.0,
                   "f·f_y·f_yyy": 2.0, "f·f_xx²": 2.0, "f·f_xy²": 4.0, "f·f_yy²": 2.0,
                   "f_x²·f_xx": 5.0, "f_x²·f_yy": 1.0, "f_y²·f_xx": 1.0, "f_y²·f_yy": 5.0,
                   "f_y·f_x·f_xy": 8.0},
    "Δ²f²": {"f·f_xxxx": 2.0, "f·f_xxyy": 4.0, "f·f_yyyy": 2.0,
             "f_x·f_xxx": 8.0, "f_x·f_xyy": 8.0, "f_y·f_xxy": 8.0, "f_y·f_yyy": 8.0,
             "f_xx²": 6.0, "f_xy²": 8.0, "f_yy²": 6.0, "f_yy·f_xx": 4.0},
    "Δ²f³": {"f²·f_xxxx": 3.0, "f²·f_xxyy": 6.0, "f²·f_yyyy": 3.0,
             "f·f_x·f_xxx": 24.0, "f·f_x·f_xyy": 24.0, "f·f_y·f_xxy": 24.0,
             "f·f_y·f_yyy": 24.0, "f·f_xx²": 18.0, "f·f_xy²": 24.0, "f·f_yy²": 18.0,
             "f·f_yy·f_xx": 12.0, "f_x²·f_xx": 36.0, "f_y·f_x·f_xy": 48.0,
             "f_y²·f_yy": 36.0, "f_x²·f_yy": 12.0, "f_y²·f_xx": 12.0},
}
_HISTORICAL_MULTIGRAPH = {**_HISTORICAL_COMMON, **_HISTORICAL_MULTIGRAPH_ONLY}
_HISTORICAL_MOD_MULTIGRAPH = {**_HISTORICAL_COMMON, **_HISTORICAL_MOD_MULTIGRAPH_ONLY}


def _as_dict(pairs):
    return dict(pairs)


class TestDecompositionMatchesHistorical2D:
    def test_multigraph(self):
        new = _multigraph_to_jet_decomposition(2)
        assert set(new) == set(_HISTORICAL_MULTIGRAPH)
        for k, expected in _HISTORICAL_MULTIGRAPH.items():
            assert _as_dict(new[k]) == expected, k

    def test_mod_multigraph(self):
        new = _mod_multigraph_to_jet_decomposition(2)
        assert set(new) == set(_HISTORICAL_MOD_MULTIGRAPH)
        for k, expected in _HISTORICAL_MOD_MULTIGRAPH.items():
            assert _as_dict(new[k]) == expected, k

    def test_module_level_dicts_match_function(self):
        """The backward-compatible module-level dicts equal
        _*_decomposition(2)."""
        assert _MULTIGRAPH_TO_JET_DECOMPOSITION == _multigraph_to_jet_decomposition(2)
        assert _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION == _mod_multigraph_to_jet_decomposition(2)


class TestDecompositionAtHighNdim:
    @pytest.mark.parametrize("ndim", [4, 5])
    def test_multigraph_all_20_operators_have_decomposition(self, ndim):
        ops = {desc for _, desc in get_operators(2, 3)}
        decomp = _multigraph_to_jet_decomposition(ndim)
        assert ops == set(decomp)

    @pytest.mark.parametrize("ndim", [4, 5])
    def test_mod_multigraph_all_20_operators_have_decomposition(self, ndim):
        ops = {desc for _, desc in get_mod_multigraph_operators(2, 3)}
        decomp = _mod_multigraph_to_jet_decomposition(ndim)
        assert ops == set(decomp)

    def test_delta_f_decomposition_scales_with_ndim(self):
        """Delta f = sum_i f_ii -- exactly ndim terms, all coefficient 1."""
        for ndim in [2, 3, 4, 5]:
            decomp = dict(_multigraph_to_jet_decomposition(ndim)["Δf"])
            assert len(decomp) == ndim
            assert all(c == 1.0 for c in decomp.values())

    def test_grad_sq_decomposition_scales_with_ndim(self):
        """||grad f||^2 = sum_i f_i^2 -- exactly ndim terms, all coefficient 1."""
        for ndim in [2, 3, 4, 5]:
            decomp = dict(_multigraph_to_jet_decomposition(ndim)["‖∇f‖²"])
            assert len(decomp) == ndim
            assert all(c == 1.0 for c in decomp.values())


class TestWsindyToJetMapNdim4:
    def test_map_matches_direct_library_column(self):
        import math
        ndim, N, L = 4, 8, 2 * math.pi
        rng = np.random.default_rng(0)
        f = 0.3 * rng.standard_normal((N,) * ndim)

        Theta_ws, desc_ws = build_wsindy_library(f, L, L, max_order=2, max_poly=1)
        Theta_jet, desc_jet = build_polyjet_library(f, L, L, max_order=2, max_degree=1)

        M = build_wsindy_to_jet_map(desc_ws, desc_jet, ndim=ndim)
        w_ws = np.zeros(len(desc_ws))
        w_ws[desc_ws.index("f_{xx}")] = 2.0
        w_jet = map_wsindy_to_jet(w_ws, M)

        direct = Theta_ws @ w_ws
        recon = Theta_jet @ w_jet
        np.testing.assert_allclose(recon, direct, atol=1e-8)
