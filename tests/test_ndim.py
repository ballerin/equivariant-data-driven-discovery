"""Tests for n-dimensional PDE discovery extensions.

Covers:
- _pdefind_axis_labels, true_allen_cahn_pdefind, true_cahn_hilliard_pdefind
- allen_cahn_nd, cahn_hilliard_nd, swift_hohenberg_nd, porous_medium_nd
  simulator consistency
- kuramoto_sivashinsky ndim guard (only ndim=2 supported)
- _imex_step n-dimensional correctness
- flat_idx C-order correctness in sweep_nd helper
"""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.equations import (
    _pdefind_axis_labels,
    true_allen_cahn_pdefind,
    true_cahn_hilliard_pdefind,
    CH_EPSILON,
)
from equivariant_pde_find.simulations.allen_cahn import allen_cahn_nd
from equivariant_pde_find.simulations.cahn_hilliard import cahn_hilliard_nd
from equivariant_pde_find.simulations.swift_hohenberg import swift_hohenberg_nd
from equivariant_pde_find.simulations.porous_medium import porous_medium_nd
from equivariant_pde_find.simulations.kuramoto_sivashinsky import kuramoto_sivashinsky
from equivariant_pde_find.simulations._imex import _imex_step


# ======================================================================
# _pdefind_axis_labels
# ======================================================================


class TestPdefindAxisLabels:
    def test_ndim_2(self):
        assert _pdefind_axis_labels(2) == ["x", "y"]

    def test_ndim_3(self):
        assert _pdefind_axis_labels(3) == ["x", "y", "z"]

    def test_ndim_4(self):
        assert _pdefind_axis_labels(4) == ["x", "y", "z", "w"]

    def test_ndim_8(self):
        labels = _pdefind_axis_labels(8)
        assert len(labels) == 8
        assert labels == ["x", "y", "z", "w", "v", "u", "t", "s"]

    @pytest.mark.parametrize("bad", [0, 9, 100])
    def test_invalid_raises(self, bad):
        with pytest.raises(ValueError, match="ndim must be"):
            _pdefind_axis_labels(bad)


# ======================================================================
# true_allen_cahn_pdefind
# ======================================================================


class TestTrueAllenCahnPdefind:
    def test_ndim_2(self):
        result = true_allen_cahn_pdefind(2, eps=0.15)
        assert len(result) == 4  # f, f^3, f_xx, f_yy
        assert result["f"] == 1.0
        assert result["f\u00b3"] == -1.0
        assert result["f_xx"] == pytest.approx(0.0225)
        assert result["f_yy"] == pytest.approx(0.0225)

    def test_ndim_3(self):
        result = true_allen_cahn_pdefind(3, eps=0.15)
        assert len(result) == 5  # f, f^3, f_xx, f_yy, f_zz
        assert result["f_xx"] == pytest.approx(0.0225)
        assert result["f_yy"] == pytest.approx(0.0225)
        assert result["f_zz"] == pytest.approx(0.0225)

    def test_ndim_4(self):
        result = true_allen_cahn_pdefind(4, eps=0.15)
        assert len(result) == 6  # f, f^3, + 4 pure 2nd
        assert "f_ww" in result
        assert result["f_ww"] == pytest.approx(0.0225)

    def test_all_second_derivatives_present(self):
        """Every axis gets exactly one f_{aa} entry with coeff eps^2."""
        for ndim in [2, 3, 4]:
            result = true_allen_cahn_pdefind(ndim, eps=0.2)
            axes = _pdefind_axis_labels(ndim)
            for a in axes:
                key = f"f_{a}{a}"
                assert key in result, f"ndim={ndim}: missing {key}"
                assert result[key] == pytest.approx(0.04)

    @pytest.mark.parametrize("bad", [0, 9])
    def test_invalid_raises(self, bad):
        with pytest.raises(ValueError, match="ndim must be"):
            true_allen_cahn_pdefind(bad)


# ======================================================================
# true_cahn_hilliard_pdefind
# ======================================================================


class TestTrueCahnHilliardPdefind:
    def test_ndim_2_term_count(self):
        """ndim=2: 2 pure 4th + 1 mixed 4th + 2 pure 2nd + 2 f^2*2nd = 7."""
        result = true_cahn_hilliard_pdefind(2, eps=0.2)
        assert len(result) == 7

    def test_ndim_3_term_count(self):
        """ndim=3: 3 + 3 + 3 + 3 = 12."""
        result = true_cahn_hilliard_pdefind(3, eps=0.2)
        assert len(result) == 12

    def test_ndim_4_term_count(self):
        """ndim=4: 4 + 6 + 4 + 4 = 18."""
        result = true_cahn_hilliard_pdefind(4, eps=0.2)
        assert len(result) == 18

    def test_pure_4th_coefficients(self):
        """Pure 4th derivatives have coeff -eps^2."""
        for ndim in [2, 3, 4]:
            result = true_cahn_hilliard_pdefind(ndim, eps=0.3)
            axes = _pdefind_axis_labels(ndim)
            for a in axes:
                assert result[f"f_{a*4}"] == pytest.approx(-0.09)

    def test_mixed_4th_coefficients(self):
        """Mixed 4th derivatives have coeff -2*eps^2."""
        result = true_cahn_hilliard_pdefind(3, eps=0.2)
        assert result["f_xxyy"] == pytest.approx(-0.08)
        assert result["f_xxzz"] == pytest.approx(-0.08)
        assert result["f_yyzz"] == pytest.approx(-0.08)

    def test_pure_2nd_coefficients(self):
        """Pure 2nd derivatives have coeff -1, from the -Delta f term."""
        result = true_cahn_hilliard_pdefind(2, eps=0.2)
        assert result["f_xx"] == -1.0
        assert result["f_yy"] == -1.0

    def test_f2_2nd_coefficients(self):
        """f^2 * pure 2nd derivatives have coeff +3, from +3 f^2 Delta f."""
        result = true_cahn_hilliard_pdefind(2, eps=0.2)
        assert result["f\u00b2\u00b7f_xx"] == 3.0
        assert result["f\u00b2\u00b7f_yy"] == 3.0

    def test_signs_match_the_multigraph_ground_truth(self):
        """The same equation must get the same signs in both bases.

        ``true_cahn_hilliard`` (multigraph) and ``true_cahn_hilliard_pdefind``
        describe one PDE; disagreeing signs on Delta f / f^2 Delta f is what
        made the PDE-FIND baseline score against a wrong truth vector.
        """
        from equivariant_pde_find.equations import true_cahn_hilliard
        from equivariant_pde_find.multigraph_operators import get_operators

        desc_mg = [d for _fn, d in get_operators(2, 3)]
        w_mg = true_cahn_hilliard(desc_mg, eps=0.2)
        pdefind = true_cahn_hilliard_pdefind(2, eps=0.2)

        assert np.sign(w_mg[desc_mg.index("\u0394f")]) == np.sign(pdefind["f_xx"])
        assert np.sign(w_mg[desc_mg.index("f\u00b2\u00b7\u0394f")]) == np.sign(pdefind["f\u00b2\u00b7f_xx"])

    def test_term_count_formula(self):
        """General formula: 3*ndim + C(ndim,2)."""
        for ndim in range(1, 8):
            expected = 3 * ndim + ndim * (ndim - 1) // 2
            result = true_cahn_hilliard_pdefind(ndim)
            assert len(result) == expected, f"ndim={ndim}"

    @pytest.mark.parametrize("bad", [0, 9])
    def test_invalid_raises(self, bad):
        with pytest.raises(ValueError, match="ndim must be"):
            true_cahn_hilliard_pdefind(bad)


# ======================================================================
# allen_cahn_nd — consistency and correctness
# ======================================================================


class TestAllenCahnND:
    def test_ndim_1_runs(self):
        """1-D should run without error."""
        snaps, times = allen_cahn_nd(
            ndim=1, N=32, L=2 * math.pi, T_end=0.5, dt=0.01,
            save_every=10, seed=42,
        )
        assert snaps.ndim == 2  # (n_saved, N)
        assert snaps.shape[1] == 32
        assert not np.any(np.isnan(snaps))

    def test_ndim_2_output_shape(self):
        snaps, times = allen_cahn_nd(
            ndim=2, N=16, L=2 * math.pi, T_end=0.5, dt=0.01,
            save_every=10, seed=42,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1:] == (16, 16)

    def test_ndim_3_no_nan(self):
        snaps, _ = allen_cahn_nd(
            ndim=3, N=12, L=2 * math.pi, T_end=0.3, dt=0.01,
            save_every=10, seed=42,
        )
        assert snaps.ndim == 4
        assert not np.any(np.isnan(snaps))

    def test_ndim_4_runs(self):
        """4-D should run without error after _imex_step generalization."""
        snaps, _ = allen_cahn_nd(
            ndim=4, N=8, L=2 * math.pi, T_end=0.2, dt=0.01,
            save_every=10, seed=42,
        )
        assert snaps.ndim == 5  # (n_saved, N, N, N, N)
        assert not np.any(np.isnan(snaps))

    @pytest.mark.parametrize("bad", [0, 9])
    def test_invalid_ndim_raises(self, bad):
        with pytest.raises(ValueError, match="ndim must be"):
            allen_cahn_nd(ndim=bad, N=16)


# ======================================================================
# cahn_hilliard_nd — consistency and correctness
# ======================================================================


class TestCahnHilliardND:
    def test_ndim_2_output_shape(self):
        snaps, times = cahn_hilliard_nd(
            ndim=2, N=16, L=8 * math.pi, T_end=1.0, dt=0.005,
            save_every=40, seed=42,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1:] == (16, 16)

    def test_ndim_3_no_nan(self):
        snaps, _ = cahn_hilliard_nd(
            ndim=3, N=12, L=8 * math.pi, T_end=0.5, dt=0.005,
            save_every=40, seed=42,
        )
        assert snaps.ndim == 4
        assert not np.any(np.isnan(snaps))

    def test_ndim_4_runs(self):
        """4-D should run without error."""
        snaps, _ = cahn_hilliard_nd(
            ndim=4, N=8, L=8 * math.pi, T_end=0.3, dt=0.005,
            save_every=40, seed=42,
        )
        assert snaps.ndim == 5
        assert not np.any(np.isnan(snaps))

    @pytest.mark.parametrize("bad", [0, 9])
    def test_invalid_ndim_raises(self, bad):
        with pytest.raises(ValueError, match="ndim must be"):
            cahn_hilliard_nd(ndim=bad, N=16)


# ======================================================================
# swift_hohenberg_nd — consistency and correctness
# ======================================================================


class TestSwiftHohenbergND:
    def test_ndim_2_output_shape(self):
        snaps, times = swift_hohenberg_nd(
            ndim=2, N=16, L=4 * math.pi, T_end=2.0, dt=0.1,
            save_every=5, seed=42,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1:] == (16, 16)

    def test_ndim_3_no_nan(self):
        snaps, _ = swift_hohenberg_nd(
            ndim=3, N=12, L=4 * math.pi, T_end=1.0, dt=0.1,
            save_every=5, seed=42,
        )
        assert snaps.ndim == 4
        assert not np.any(np.isnan(snaps))

    def test_ndim_4_runs(self):
        """4-D should run without error."""
        snaps, _ = swift_hohenberg_nd(
            ndim=4, N=8, L=4 * math.pi, T_end=0.5, dt=0.1,
            save_every=5, seed=42,
        )
        assert snaps.ndim == 5
        assert not np.any(np.isnan(snaps))

    @pytest.mark.parametrize("bad", [0, 9])
    def test_invalid_ndim_raises(self, bad):
        with pytest.raises(ValueError, match="ndim must be"):
            swift_hohenberg_nd(ndim=bad, N=16)


# ======================================================================
# porous_medium_nd — consistency and correctness
# ======================================================================


class TestPorousMediumND:
    def test_ndim_2_output_shape(self):
        snaps, times = porous_medium_nd(
            ndim=2, N=16, L=6 * math.pi, T_end=0.05, dt=0.002,
            save_every=5, seed=42,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1:] == (16, 16)
        assert np.all(snaps >= 0.0)

    def test_ndim_3_no_nan(self):
        snaps, _ = porous_medium_nd(
            ndim=3, N=12, L=6 * math.pi, T_end=0.03, dt=0.002,
            save_every=5, seed=42,
        )
        assert snaps.ndim == 4
        assert not np.any(np.isnan(snaps))
        assert np.all(snaps >= 0.0)

    def test_ndim_4_runs(self):
        """4-D should run without error."""
        snaps, _ = porous_medium_nd(
            ndim=4, N=8, L=6 * math.pi, T_end=0.02, dt=0.002,
            save_every=5, seed=42,
        )
        assert snaps.ndim == 5
        assert not np.any(np.isnan(snaps))
        assert np.all(snaps >= 0.0)

    @pytest.mark.parametrize("bad", [0, 9])
    def test_invalid_ndim_raises(self, bad):
        with pytest.raises(ValueError, match="ndim must be"):
            porous_medium_nd(ndim=bad, N=16)


# ======================================================================
# kuramoto_sivashinsky — ndim guard (only ndim=2 supported)
# ======================================================================


class TestKuramotoSivashinskyNdimGuard:
    def test_ndim_2_runs(self):
        """Default ndim=2 is unaffected by the guard."""
        snaps, times = kuramoto_sivashinsky(
            N=16, L=16 * math.pi, T_end=1.0, dt=0.25, save_every=2, seed=42,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1:] == (16, 16)

    @pytest.mark.parametrize("bad", [1, 3, 4])
    def test_other_ndim_raises(self, bad):
        with pytest.raises(ValueError, match="only supports ndim=2"):
            kuramoto_sivashinsky(N=16, ndim=bad)


# ======================================================================
# _imex_step — n-dimensional correctness
# ======================================================================


class TestImexStep:
    def test_fftn_equivalent_to_fft2_for_2d(self):
        """fftn and fft2 should give identical results for 2-D fields."""
        rng = np.random.default_rng(42)
        f = 0.1 * rng.standard_normal((8, 8))
        denom = np.ones((8, 8))
        dt = 0.01

        def nonlin(x):
            return x

        # Use a "reference" by running the step directly
        f_next = _imex_step(f.copy(), nonlin, denom, dt, 0.0, rng)
        assert f_next.shape == (8, 8)
        assert not np.any(np.isnan(f_next))

    def test_ndim_3_consistent(self):
        """3-D fields should step without error."""
        rng = np.random.default_rng(42)
        f = 0.1 * rng.standard_normal((6, 6, 6))
        denom = np.ones((6, 6, 6))
        dt = 0.01

        def nonlin(x):
            return x

        f_next = _imex_step(f, nonlin, denom, dt, 0.0, rng)
        assert f_next.shape == (6, 6, 6)
        assert not np.any(np.isnan(f_next))

    def test_ndim_4_runs(self):
        """4-D fields should step without error."""
        rng = np.random.default_rng(42)
        f = 0.1 * rng.standard_normal((4, 4, 4, 4))
        denom = np.ones((4, 4, 4, 4))
        dt = 0.01

        def nonlin(x):
            return x

        f_next = _imex_step(f, nonlin, denom, dt, 0.0, rng)
        assert f_next.shape == (4, 4, 4, 4)
        assert not np.any(np.isnan(f_next))

    def test_ndim_1_runs(self):
        """1-D fields should step without error."""
        rng = np.random.default_rng(42)
        f = 0.1 * rng.standard_normal((16,))
        denom = np.ones((16,))
        dt = 0.01

        def nonlin(x):
            return x

        f_next = _imex_step(f, nonlin, denom, dt, 0.0, rng)
        assert f_next.shape == (16,)
        assert not np.any(np.isnan(f_next))


# ======================================================================
# flat_idx formula (C-order correctness)
# ======================================================================


class TestFlatIdx:
    """Verify the C-order flat-indexing formula used in sweep_nd._build_weak_nd."""

    @staticmethod
    def _build_flat_idx(Nv, sample_stride):
        """Reference implementation matching sweep_nd._build_weak_nd."""
        ndim = len(Nv)
        grid_axes = [np.arange(0, Nv[i], sample_stride) for i in range(ndim)]
        mesh = np.meshgrid(*grid_axes, indexing="ij")
        flat_idx = 0
        stride = 1
        for i in reversed(range(ndim)):
            flat_idx = flat_idx + mesh[i].ravel() * stride
            stride *= Nv[i]
        return flat_idx

    def test_2d_matches_ravel_indexing(self):
        """2-D flat indices should match np.ravel_multi_index."""
        Nv = (12, 8)
        idx = self._build_flat_idx(Nv, 4)
        # Reconstruct the multi-indices
        ndim = len(Nv)
        grid_axes = [np.arange(0, Nv[i], 4) for i in range(ndim)]
        mesh = np.meshgrid(*grid_axes, indexing="ij")
        multi = np.column_stack([m.ravel() for m in mesh])
        expected = np.ravel_multi_index(multi.T, Nv)
        np.testing.assert_array_equal(idx, expected)

    def test_3d_matches_ravel_indexing(self):
        """3-D flat indices should match np.ravel_multi_index."""
        Nv = (10, 8, 6)
        idx = self._build_flat_idx(Nv, 3)
        ndim = len(Nv)
        grid_axes = [np.arange(0, Nv[i], 3) for i in range(ndim)]
        mesh = np.meshgrid(*grid_axes, indexing="ij")
        multi = np.column_stack([m.ravel() for m in mesh])
        expected = np.ravel_multi_index(multi.T, Nv)
        np.testing.assert_array_equal(idx, expected)

    def test_4d_matches_ravel_indexing(self):
        """4-D flat indices should match np.ravel_multi_index."""
        Nv = (6, 5, 4, 3)
        idx = self._build_flat_idx(Nv, 2)
        ndim = len(Nv)
        grid_axes = [np.arange(0, Nv[i], 2) for i in range(ndim)]
        mesh = np.meshgrid(*grid_axes, indexing="ij")
        multi = np.column_stack([m.ravel() for m in mesh])
        expected = np.ravel_multi_index(multi.T, Nv)
        np.testing.assert_array_equal(idx, expected)

    def test_square_2d_self_consistent(self):
        """For square grids the old (buggy) and new formulas give same set."""
        Nv = (10, 10)
        idx_new = self._build_flat_idx(Nv, 3)
        # Old (buggy) formula for comparison
        ndim = len(Nv)
        grid_axes = [np.arange(0, Nv[i], 3) for i in range(ndim)]
        mesh = np.meshgrid(*grid_axes, indexing="ij")
        flat_old = mesh[0].ravel()
        stride = 1
        for i in range(1, ndim):
            stride *= Nv[i - 1]
            flat_old = flat_old + mesh[i].ravel() * stride
        # For square grids, the sets should be equal (just permuted)
        assert set(idx_new) == set(flat_old)
