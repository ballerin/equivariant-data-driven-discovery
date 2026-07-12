"""
Shared fixtures for the equivariant-pde-find test suite.

Fixtures defined here are automatically available to all test files
without explicit imports, per pytest's conftest.py convention.

Fixtures
--------
test_field_2d  : (f, L) — a smooth band-limited 2-D field on [0, 2π)²
test_field_1d  : (u, L) — a smooth band-limited 1-D field on [0, 2π)
sparse_system  : factory → (G, b, w_true, nz_idx) for synthetic regression
ac_data        : (snaps, times, true_eq, desc) — small Allen-Cahn trajectory
"""

import math

import numpy as np
import pytest

# Ensure the package is importable regardless of how pytest is invoked
import sys
import os

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if CODE_DIR not in sys.path:
    sys.path.insert(0, CODE_DIR)


# ── 2-D field fixtures ────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def test_field_2d():
    """A smooth, band-limited 2-D field on [0, 2π)² with multiple modes."""
    N, L = 64, 2 * math.pi
    x = np.linspace(0.0, L, N, endpoint=False)
    y = np.linspace(0.0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="ij")
    f = np.cos(2 * X + 3 * Y) + 0.5 * np.sin(X - 2 * Y)
    return f, L


@pytest.fixture(scope="session")
def test_field_2d_small():
    """Smaller 2-D field for faster tests."""
    N, L = 32, 2 * math.pi
    x = np.linspace(0.0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    f = np.sin(X) + 0.5 * np.cos(2 * Y)
    return f, L


# ── 1-D field fixtures ────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def test_field_1d():
    """A smooth, band-limited 1-D field on [0, 2π)."""
    N, L = 64, 2 * math.pi
    x = np.linspace(0.0, L, N, endpoint=False)
    u = np.sin(x) + 0.3 * np.cos(2 * x)
    return u, L


# ── Sparse system factory ──────────────────────────────────────────────────


@pytest.fixture
def sparse_system():
    """Return a factory that generates synthetic sparse linear systems.

    Usage in tests::

        G, b, w_true, nz_idx = sparse_system(n=200, d=20, nz=4, noise=0.0)
    """

    def _make(n=200, d=20, nz=4, noise=0.0, seed=42, normalize=True):
        rng = np.random.default_rng(seed)
        G = rng.standard_normal((n, d))
        if normalize:
            G = G / np.sqrt(np.sum(G ** 2, axis=0))[np.newaxis, :]
        w_true = np.zeros(d)
        nz_idx = rng.choice(d, size=nz, replace=False)
        w_true[nz_idx] = rng.uniform(0.5, 2.0, size=nz) * rng.choice(
            [-1, 1], size=nz,
        )
        b = G @ w_true
        if noise > 0:
            b += noise * rng.standard_normal(n)
        return G, b, w_true, nz_idx

    return _make


# ── Allen-Cahn trajectory fixture ─────────────────────────────────────────


@pytest.fixture(scope="session")
def ac_data():
    """Simulate a small Allen-Cahn trajectory for integration tests.

    Returns (snapshots, times, true_eq, descriptions) where
    true_eq is the ground-truth coefficient vector in the equivariant basis.
    """
    from equivariant_pde_find.simulations import allen_cahn
    from equivariant_pde_find.library import build_multigraph_library

    snaps, times = allen_cahn(
        N=32, L=2 * np.pi, T_end=2.0, eps=0.15, dt=0.01,
        save_every=20, seed=42,
    )
    _, desc = build_multigraph_library(snaps[0], 2 * np.pi, 2 * np.pi)

    # Ground-truth Allen-Cahn: f_t = f - f³ + ε²·Δf
    w = np.zeros(len(desc))
    w[desc.index("f")] = 1.0
    w[desc.index("f\u00b3")] = -1.0
    w[desc.index("\u0394f")] = 0.15 ** 2
    return snaps, times, w, desc


# ── Configure to always use spectral backend for consistency ──────────────


@pytest.fixture(autouse=True)
def _reset_derivatives_backend():
    """Ensure each test starts with the 'spectral' backend active."""
    from equivariant_pde_find import derivatives

    derivatives.set_backend("spectral")
    yield
    derivatives.set_backend("spectral")
