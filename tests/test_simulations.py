"""Smoke tests for simulations.py PDE simulators."""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.simulations import (
    allen_cahn,
    cahn_hilliard,
    kuramoto_sivashinsky,
    make_regression_data,
    porous_medium,
    swift_hohenberg,
)


class TestMakeRegressionData:
    def test_output_shapes(self):
        N, Nt = 32, 10
        snaps = np.random.default_rng(0).standard_normal((Nt, N, N))
        times = np.arange(Nt) * 0.1
        f_snaps, ft_snaps, idx_sets, dt = make_regression_data(
            snaps, times, n_per_snap=100, skip_fraction=0.0, seed=0,
        )
        assert f_snaps.shape == (Nt - 2, N, N)
        assert ft_snaps.shape == (Nt - 2, N, N)
        assert len(idx_sets) == Nt - 2
        assert all(len(idx) == 100 for idx in idx_sets)
        assert dt == pytest.approx(0.1)

    def test_skip_fraction(self):
        N, Nt = 32, 20
        snaps = np.random.default_rng(0).standard_normal((Nt, N, N))
        times = np.arange(Nt) * 0.1
        f_snaps, ft_snaps, idx_sets, dt = make_regression_data(
            snaps, times, n_per_snap=50, skip_fraction=0.5, seed=0,
        )
        # skip_fraction=0.5 → keep last 10 snapshots → 8 interior
        assert f_snaps.shape[0] <= 8


class TestAllenCahn:
    def test_runs_and_produces_output(self):
        N, L = 16, 2 * np.pi
        snaps, times = allen_cahn(
            N=N, L=L, T_end=0.01, eps=0.3, dt=0.005, save_every=2, seed=0,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1] == N
        assert snaps.shape[2] == N
        assert len(times) == snaps.shape[0]
        assert np.all(np.isfinite(snaps))


class TestSwiftHohenberg:
    def test_runs_and_produces_output(self):
        N, L = 16, 8 * np.pi
        snaps, times = swift_hohenberg(
            N=N, L=L, T_end=0.01, r=0.3, dt=0.01, save_every=1, seed=0,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1] == N
        assert len(times) == snaps.shape[0]
        assert np.all(np.isfinite(snaps))


class TestKuramotoSivashinsky:
    def test_runs_and_produces_output(self):
        N, L = 16, 4 * np.pi
        snaps, times = kuramoto_sivashinsky(
            N=N, L=L, T_end=0.05, dt=0.05, save_every=1, seed=0,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1] == N
        assert len(times) == snaps.shape[0]
        assert np.all(np.isfinite(snaps))


class TestCahnHilliard:
    def test_runs_and_produces_output(self):
        N, L = 16, 2 * np.pi
        snaps, times = cahn_hilliard(
            N=N, L=L, T_end=0.01, eps=0.2, dt=0.005, save_every=2, seed=0,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1] == N
        assert len(times) == snaps.shape[0]
        assert np.all(np.isfinite(snaps))


class TestPorousMedium:
    def test_runs_and_produces_output(self):
        N, L = 16, 2 * np.pi
        snaps, times = porous_medium(
            N=N, L=L, T_end=0.01, dt=0.005, save_every=2, seed=0,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1] == N
        assert len(times) == snaps.shape[0]
        assert np.all(np.isfinite(snaps))
        # PME should preserve non-negativity
        assert np.all(snaps >= 0)
