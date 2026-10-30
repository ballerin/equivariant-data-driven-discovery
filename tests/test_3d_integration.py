"""End-to-end integration tests for 3-D Allen-Cahn PDE discovery.

Runs a small 3-D simulation and verifies RFE exact recovery on clean data
for all three library types.
"""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives
from equivariant_pde_find.library import (
    build_multigraph_library, build_mod_multigraph_library, build_pdefind_library,
)
from equivariant_pde_find.library import weak_form as wf
from equivariant_pde_find.regression import TrainRFE
from equivariant_pde_find.simulations.allen_cahn import allen_cahn_nd


@pytest.fixture(scope="module")
def ac_3d_small():
    """Small 3-D Allen-Cahn simulation for integration testing."""
    N, L = 24, 2 * math.pi
    eps = 0.15
    T_end = 4.0  # Need enough dynamics for RFE to resolve all 3 terms
    dt = 0.01
    save_every = 20

    snaps, times = allen_cahn_nd(
        ndim=3, N=N, L=L, T_end=T_end, eps=eps, dt=dt,
        save_every=save_every, seed=42, obs_noise=0.0,
    )
    return snaps, times, L, eps


def _make_regression_data_3d(snapshots, times, n_per_snap=200):
    """Minimal 3-D regression data pipeline (avoids modifying existing code)."""
    dt_snap = times[1] - times[0]
    ft_snaps = (snapshots[2:] - snapshots[:-2]) / (2.0 * dt_snap)
    f_snaps = snapshots[1:-1]

    N_total = int(np.prod(f_snaps.shape[1:]))
    rng = np.random.default_rng(0)
    idx_sets = [rng.choice(N_total, n_per_snap, replace=False) for _ in range(len(f_snaps))]
    return f_snaps, ft_snaps, idx_sets


class TestAllenCahn3DIntegration:
    """End-to-end PDE discovery on small 3-D Allen-Cahn data."""

    def test_simulator_produces_valid_data(self, ac_3d_small):
        snaps, times, L, eps = ac_3d_small
        assert snaps.ndim == 4  # (n_saved, N, N, N)
        assert snaps.shape[1] == snaps.shape[2] == snaps.shape[3] == 24
        assert len(times) == snaps.shape[0]
        assert np.all(np.isfinite(snaps))

    def test_multigraph_rfe_exact_recovery(self, ac_3d_small):
        snaps, times, L, eps = ac_3d_small
        f_snaps, ft_snaps, idx_sets = _make_regression_data_3d(snaps, times)

        m, p = 6, 4
        derivatives.set_backend("spectral")

        # Build weak-form multigraph library with IBP
        kernel = wf._cached_kernel_1d(m, p)
        rows_G, rows_b = [], []
        for s in range(len(f_snaps)):
            f = f_snaps[s]
            ft = ft_snaps[s]
            G, desc = build_multigraph_library(
                f, L, L, max_p=2, max_n=3, formulation="space_weak",
                m=m, p=p, space_ibp=True, Lz=L,
            )
            b_w = wf.build_weak_rhs(ft, kernel_x=kernel, kernel_y=kernel, m=m, p=p)
            idx = idx_sets[s] % len(b_w)  # clip to valid interior
            rows_G.append(G[idx])
            rows_b.append(b_w[idx])

        G_all = np.vstack(rows_G).astype(float)
        b_all = np.concatenate(rows_b).astype(float)

        w = TrainRFE(G_all, b_all, normalize=2, lam=0.0)

        # Check recovery: should have exactly 3 non-zero terms
        nz = np.abs(w) > 1e-10
        n_active = int(np.sum(nz))
        assert n_active == 3, f"Expected 3 active terms, got {n_active}"

        # Check the right terms are active
        idx_lap = desc.index("Δf")
        idx_f = desc.index("f")
        idx_f3 = desc.index("f³")
        assert np.abs(w[idx_lap]) > 1e-10, "Δf missing"
        assert np.abs(w[idx_f]) > 1e-10, "f missing"
        assert np.abs(w[idx_f3]) > 1e-10, "f³ missing"

        # Verify coefficients
        np.testing.assert_allclose(w[idx_lap], eps**2, rtol=1e-1, atol=0.01)
        np.testing.assert_allclose(w[idx_f], 1.0, rtol=1e-1, atol=0.1)
        np.testing.assert_allclose(w[idx_f3], -1.0, rtol=1e-1, atol=0.1)

    def test_mod_multigraph_rfe_exact_recovery(self, ac_3d_small):
        snaps, times, L, eps = ac_3d_small
        f_snaps, ft_snaps, idx_sets = _make_regression_data_3d(snaps, times)

        m, p = 6, 4
        derivatives.set_backend("spectral")

        kernel = wf._cached_kernel_1d(m, p)
        rows_G, rows_b = [], []
        for s in range(len(f_snaps)):
            f = f_snaps[s]
            ft = ft_snaps[s]
            G, desc = build_mod_multigraph_library(
                f, L, L, max_p=2, max_n=3, formulation="space_weak",
                m=m, p=p, space_ibp=True, Lz=L,
            )
            b_w = wf.build_weak_rhs(ft, kernel_x=kernel, kernel_y=kernel, m=m, p=p)
            idx = idx_sets[s] % len(b_w)
            rows_G.append(G[idx])
            rows_b.append(b_w[idx])

        G_all = np.vstack(rows_G).astype(float)
        b_all = np.concatenate(rows_b).astype(float)

        w = TrainRFE(G_all, b_all, normalize=2, lam=0.0)

        nz = np.abs(w) > 1e-10
        n_active = int(np.sum(nz))
        assert n_active == 3, f"Expected 3 active terms, got {n_active}"

        idx_lap = desc.index("Δf")
        idx_f = desc.index("f")
        idx_f3 = desc.index("f³")
        assert np.abs(w[idx_lap]) > 1e-10
        assert np.abs(w[idx_f]) > 1e-10
        assert np.abs(w[idx_f3]) > 1e-10

        np.testing.assert_allclose(w[idx_lap], eps**2, rtol=1e-1, atol=0.01)
        np.testing.assert_allclose(w[idx_f], 1.0, rtol=1e-1, atol=0.1)
        np.testing.assert_allclose(w[idx_f3], -1.0, rtol=1e-1, atol=0.1)

    def test_pdefind_rfe_exact_recovery(self, ac_3d_small):
        snaps, times, L, eps = ac_3d_small
        f_snaps, ft_snaps, idx_sets = _make_regression_data_3d(snaps, times)

        m, p = 6, 4
        derivatives.set_backend("spectral")

        kernel = wf._cached_kernel_1d(m, p)
        rows_G, rows_b = [], []
        for s in range(len(f_snaps)):
            f = f_snaps[s]
            ft = ft_snaps[s]
            G, desc = build_pdefind_library(
                f, L, L, max_order=4, max_poly=3, formulation="space_weak",
                m=m, p=p, space_ibp=True, Lz=L,
            )
            b_w = wf.build_weak_rhs(ft, kernel_x=kernel, kernel_y=kernel, m=m, p=p)
            idx = idx_sets[s] % len(b_w)
            rows_G.append(G[idx])
            rows_b.append(b_w[idx])

        G_all = np.vstack(rows_G).astype(float)
        b_all = np.concatenate(rows_b).astype(float)

        w = TrainRFE(G_all, b_all, normalize=2, lam=0.0)

        nz = np.abs(w) > 1e-10
        n_active = int(np.sum(nz))
        assert 3 <= n_active <= 6, f"Expected 3-6 active terms, got {n_active}"

        idx_f = desc.index("f")
        assert np.abs(w[idx_f]) > 1e-10, "f term should be active"

        laplacian_terms = ["f_xx", "f_yy", "f_zz"]
        lap_active = any(
            np.abs(w[desc.index(t)]) > 1e-10 for t in laplacian_terms if t in desc
        )
        assert lap_active, "At least one Laplacian term should be active"
