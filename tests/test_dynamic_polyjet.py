"""Tests for ``build_polyjet_from_descriptions`` (dynamic polyjet builder)."""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.library import (
    build_polyjet_library,
    build_polyjet_from_descriptions,
    build_multigraph_library,
)
from equivariant_pde_find.library.conversions import (
    collect_needed_jet_monomials_from_multigraph,
    _MULTIGRAPH_TO_JET_DECOMPOSITION,
)


class TestDynamicPolyjet:
    """Tests for build_polyjet_from_descriptions."""

    def test_dynamic_matches_capped_for_same_descriptions(self):
        """When given the same 105 descriptions as the capped builder,
        the dynamic builder produces identical Theta."""
        N, L = 16, 2 * math.pi
        rng = np.random.default_rng(42)
        f = rng.standard_normal((N, N))

        Theta_capped, descs = build_polyjet_library(
            f, L, L, max_order=4, max_degree=3,
        )
        Theta_dynamic, descs2 = build_polyjet_from_descriptions(
            f, L, L, descs,
        )

        assert descs == descs2
        np.testing.assert_allclose(Theta_capped, Theta_dynamic, atol=1e-12)

    def test_dynamic_subset_matches_capped_subset(self):
        """Dynamic builder with a subset of capped descriptions matches
        the corresponding columns of the capped builder."""
        N, L = 16, 2 * math.pi
        rng = np.random.default_rng(43)
        f = rng.standard_normal((N, N))

        Theta_capped, descs_capped = build_polyjet_library(
            f, L, L, max_order=4, max_degree=3,
        )
        # Pick a subset: some common monomials
        subset = ["1", "f", "f_x", "f_xx", "f³", "f·f_x", "f_x²"]
        Theta_dynamic, descs_dyn = build_polyjet_from_descriptions(
            f, L, L, subset,
        )

        assert descs_dyn == subset
        for j, desc in enumerate(subset):
            idx = descs_capped.index(desc)
            np.testing.assert_allclose(
                Theta_dynamic[:, j], Theta_capped[:, idx], atol=1e-12,
                err_msg=f"Column {desc!r} differs",
            )

    def test_dynamic_includes_degree_four(self):
        """Degree-4 monomials (e.g. f³·f_x) can be built."""
        N, L = 16, 2 * math.pi
        rng = np.random.default_rng(44)
        f = rng.standard_normal((N, N))

        descs = ["f³·f_x", "f³·f_xx"]
        Theta, descs_out = build_polyjet_from_descriptions(f, L, L, descs)

        assert descs_out == descs
        assert Theta.shape == (N * N, 2)
        assert np.all(np.isfinite(Theta))
        # Verify columns are not all zero and not constant
        assert np.std(Theta[:, 0]) > 0, "f³·f_x column should vary"
        assert np.std(Theta[:, 1]) > 0, "f³·f_xx column should vary"
        # f³·f_x should match manual computation: f^3 * f_x
        _, descs_full = build_polyjet_library(f, L, L, max_order=4, max_degree=3)
        # f_x is a jet variable, build via dynamic for the "f_x" column
        T_fx, _ = build_polyjet_from_descriptions(f, L, L, ["f_x"])
        T_fxx, _ = build_polyjet_from_descriptions(f, L, L, ["f_xx"])
        expected_f3_fx = (f ** 3 * T_fx.reshape(N, N)).ravel()
        expected_f3_fxx = (f ** 3 * T_fxx.reshape(N, N)).ravel()
        np.testing.assert_allclose(Theta[:, 0], expected_f3_fx, atol=1e-10)
        np.testing.assert_allclose(Theta[:, 1], expected_f3_fxx, atol=1e-10)

    def test_dynamic_constant_column(self):
        """The '1' column is all ones."""
        N, L = 16, 2 * math.pi
        rng = np.random.default_rng(45)
        f = rng.standard_normal((N, N))

        Theta, descs = build_polyjet_from_descriptions(f, L, L, ["1", "f"])
        np.testing.assert_allclose(Theta[:, 0], 1.0, atol=1e-12)
        np.testing.assert_allclose(Theta[:, 1], f.ravel(), atol=1e-12)

    def test_dynamic_from_collector_equivariant(self):
        """Dynamic builder with multigraph-derived descriptions produces
        the same PDE prediction as the multigraph library."""
        N, L = 32, 2 * math.pi
        rng = np.random.default_rng(99)
        f = rng.standard_normal((N, N))

        T_mg, desc_mg = build_multigraph_library(f, L, L, max_p=2, max_n=3)
        desc_jet = collect_needed_jet_monomials_from_multigraph(desc_mg)
        T_jet, _ = build_polyjet_from_descriptions(f, L, L, desc_jet)

        # Allen-Cahn ground truth in both libraries
        from equivariant_pde_find.equations import true_allen_cahn
        from equivariant_pde_find.library.conversions import map_multigraph_to_jet

        w_mg = true_allen_cahn(desc_mg, eps=0.15)
        w_jet = map_multigraph_to_jet(w_mg, desc_mg, desc_jet)

        pred_mg = T_mg @ w_mg
        pred_jet = T_jet @ w_jet
        scale = max(np.max(np.abs(pred_mg)), 1e-10)
        rel_err = np.max(np.abs(pred_jet - pred_mg)) / scale
        assert rel_err < 1e-10, f"Multigraph→jet prediction error {rel_err:.2e}"
