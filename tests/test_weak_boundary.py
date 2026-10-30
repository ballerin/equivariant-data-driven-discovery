"""
Tests for the periodic spatial boundary of the weak form
(``weak_form.use_weak_boundary``) and for the 2026-09 renaming: library
aliases (``formulation="weak"``, ``use_ibp``), sweep-CSV normalisation, and
the sweep's method / boundary configurations.
"""

import math
import os
import sys

import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

from equivariant_pde_find import derivatives  # noqa: E402
from equivariant_pde_find.equations import find_op  # noqa: E402
from equivariant_pde_find.library import (  # noqa: E402
    build_mod_multigraph_library, build_multigraph_library,
    build_pdefind_library, build_polyjet_library, build_wsindy_library,
)
from equivariant_pde_find.library import weak_form as wf  # noqa: E402
from equivariant_pde_find.library.spacetime import build_weak_system  # noqa: E402
from equivariant_pde_find.summary_metadata import (  # noqa: E402
    METHOD_MAP, canonical_method, normalize_sweep_row,
)


def small_builders():
    return {
        "multigraph": lambda f, L, **kw: build_multigraph_library(f, L, L, **kw),
        "mod_multigraph": lambda f, L, **kw: build_mod_multigraph_library(
            f, L, L, **kw),
        "pdefind": lambda f, L, **kw: build_pdefind_library(
            f, L, L, max_order=4, max_poly=2, **kw),
        "wsindy": lambda f, L, **kw: build_wsindy_library(
            f, L, L, max_order=4, max_poly=2, **kw),
        "polyjet": lambda f, L, **kw: build_polyjet_library(
            f, L, L, max_order=2, max_degree=2, **kw),
    }


@pytest.fixture(scope="module")
def periodic_field():
    N, L = 40, 2 * math.pi
    x = np.linspace(0.0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    f = 0.8 * np.cos(X + 2 * Y) + 0.5 * np.sin(3 * X) - 0.3 * np.cos(2 * Y)
    return f, L


# ---------------------------------------------------------------------------
# Boundary switch
# ---------------------------------------------------------------------------


def test_default_is_valid():
    assert wf.get_weak_boundary() == "valid"


def test_context_manager_restores_mode_even_on_error():
    with pytest.raises(RuntimeError):
        with wf.use_weak_boundary("periodic"):
            assert wf.get_weak_boundary() == "periodic"
            raise RuntimeError
    assert wf.get_weak_boundary() == "valid"
    with wf.use_weak_boundary(None):
        assert wf.get_weak_boundary() == "valid"


def test_unknown_boundary_rejected():
    with pytest.raises(ValueError):
        with wf.use_weak_boundary("reflect"):
            pass
    assert wf.get_weak_boundary() == "valid"


def test_weak_output_shape():
    assert wf.weak_output_shape((32, 40), 5) == (22, 30)
    assert wf.weak_output_shape((32, 40), 5, "periodic") == (32, 40)
    with wf.use_weak_boundary("periodic"):
        assert wf.weak_output_shape((32, 40), 5) == (32, 40)


# ---------------------------------------------------------------------------
# Periodic convolution
# ---------------------------------------------------------------------------


def test_periodic_convolution_centred_and_wrapped():
    """Entry i is centred at grid point i: the interior equals the 'valid'
    result, and the whole array equals 'valid' on a wrap-padded field."""
    rng = np.random.default_rng(0)
    f = rng.standard_normal((30, 24, 18))
    m = 5
    kernels = [wf._test_function_deriv_1d(m, 9, 1),
               wf._test_function_1d(m, 9),
               wf._test_function_deriv_1d(m, 9, 2)]
    valid = wf.weak_convolution(f, kernels)
    with wf.use_weak_boundary("periodic"):
        per = wf.weak_convolution(f, kernels)
    assert per.shape == f.shape
    np.testing.assert_allclose(per[m:-m, m:-m, m:-m], valid, atol=1e-12)
    padded = np.pad(f, m, mode="wrap")
    np.testing.assert_allclose(wf.weak_convolution(padded, kernels), per,
                               atol=1e-12)


def test_periodic_kernel_longer_than_axis_rejected():
    with wf.use_weak_boundary("periodic"):
        with pytest.raises(ValueError):
            wf.weak_convolution(np.zeros((8, 8)), wf._test_function_1d(5, 6),
                                wf._test_function_1d(3, 6))


@pytest.mark.parametrize("space_ibp", [False, True])
@pytest.mark.parametrize("library", sorted(small_builders()))
def test_every_library_honours_periodic_boundary(periodic_field, library,
                                                 space_ibp):
    """All five builders go through weak_convolution: under 'periodic' they
    return one row per grid point, and the interior rows are the 'valid'
    rows."""
    f, L = periodic_field
    N, m, p = f.shape[0], 6, 8
    build = small_builders()[library]
    with derivatives.use_backend("spectral"):
        G_valid, desc = build(f, L, formulation="space_weak", m=m, p=p,
                              space_ibp=space_ibp)
        with wf.use_weak_boundary("periodic"):
            G_per, desc_per = build(f, L, formulation="space_weak", m=m, p=p,
                                    space_ibp=space_ibp)
    assert desc_per == desc
    K = G_valid.shape[1]
    assert G_per.shape == (N * N, K)
    np.testing.assert_allclose(
        G_per.reshape(N, N, K)[m:N - m, m:N - m],
        G_valid.reshape(N - 2 * m, N - 2 * m, K), rtol=1e-10, atol=1e-10)


def test_periodic_ibp_matches_direct_everywhere(periodic_field):
    """On a torus integration by parts has no boundary terms, so the IBP
    Laplacian column equals the direct one at every grid point, including
    centres whose support wraps around the domain."""
    f, L = periodic_field
    m = 12
    p = wf.select_test_function_degree(m, 4)
    with derivatives.use_backend("spectral"), wf.use_weak_boundary("periodic"):
        G_ibp, desc = build_multigraph_library(f, L, L, formulation="space_weak",
                                               m=m, p=p, space_ibp=True)
        G_dir, _ = build_multigraph_library(f, L, L, formulation="space_weak",
                                            m=m, p=p, space_ibp=False)
    N = f.shape[0]
    interior = np.zeros((N, N), dtype=bool)
    interior[m:N - m, m:N - m] = True
    interior = interior.ravel()
    for name in ("Δf", "Δ²f"):
        j = find_op(desc, name)
        err = np.abs(G_ibp[:, j] - G_dir[:, j]) / np.abs(G_dir[:, j]).max()
        # what remains is the quadrature error of the degree-p bump (~1e-7
        # for Δ²), the same at wrapped centres as in the 'valid' interior
        assert err.max() < 1e-6, (name, err.max())
        assert err[~interior].max() < 10 * err[interior].max() + 1e-12, name


# ---------------------------------------------------------------------------
# build_weak_system with a periodic boundary
# ---------------------------------------------------------------------------


NU, SOURCE = 0.5, 0.3


@pytest.fixture(scope="module")
def heat_data():
    """Exact solution of u_t = ν Δu + c on the periodic square [0, 2π)²."""
    N, L, Nt, dt = 32, 2 * math.pi, 61, 0.02
    x = np.linspace(0.0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    t = np.arange(Nt) * dt
    S = np.zeros((Nt, N, N))
    for (kx, ky), a in [((1, 2), 1.0), ((3, 1), 0.6), ((2, 2), -0.4)]:
        S += (a * np.exp(-NU * (kx ** 2 + ky ** 2) * t)[:, None, None]
              * np.cos(kx * X + ky * Y)[None])
    S += SOURCE * t[:, None, None] + 1.0
    return S, t, L


def test_build_weak_system_periodic_keeps_every_query_point(heat_data):
    """m close to N/2 leaves one valid centre per axis but every grid point
    under 'periodic', and the exact solution still has a tiny residual."""
    S, t, L = heat_data
    N = S.shape[1]
    m = N // 2 - 1
    p = wf.select_test_function_degree(m, 4)

    def build(f, formulation, **kw):
        return build_multigraph_library(f, L, L, formulation=formulation, **kw)

    kw = dict(m=m, p=p, space_ibp=True, m_t=8,
              p_t=wf.select_test_function_degree(8, 1))
    with derivatives.use_backend("spectral"):
        _, _, _, info_v = build_weak_system(S, t, build, "spacetime_weak", **kw)
        G, b, desc, info_p = build_weak_system(S, t, build, "spacetime_weak",
                                               boundary="periodic", **kw)
    assert wf.get_weak_boundary() == "valid"
    assert info_v["n_space_query"] == 4          # (N − 2m)² = 2²
    assert info_p["n_space_query"] == N * N
    assert info_p["boundary"] == "periodic"
    w = np.zeros(len(desc))
    w[find_op(desc, "Δf")], w[find_op(desc, "1")] = NU, SOURCE
    assert np.linalg.norm(G @ w - b) / np.linalg.norm(b) < 1e-6


def test_build_weak_system_boundary_ignored_when_pointwise_in_space(heat_data):
    S, t, L = heat_data

    def build(f, formulation, **kw):
        return build_multigraph_library(f, L, L, formulation=formulation, **kw)

    with derivatives.use_backend("spectral"):
        _, _, _, info = build_weak_system(
            S, t, build, "time_weak", m_t=8,
            p_t=wf.select_test_function_degree(8, 1), boundary="periodic")
    assert info["boundary"] is None
    assert info["n_space_query"] == S.shape[1] ** 2


# ---------------------------------------------------------------------------
# Library aliases
# ---------------------------------------------------------------------------


def test_weak_and_use_ibp_are_aliases(periodic_field):
    f, L = periodic_field
    G_new, d_new = build_wsindy_library(f, L, L, formulation="space_weak",
                                        m=6, p=8, space_ibp=True)
    G_old, d_old = build_wsindy_library(f, L, L, formulation="weak",
                                        m=6, p=8, use_ibp=True)
    assert d_new == d_old
    np.testing.assert_array_equal(G_new, G_old)


def test_alias_and_new_kwarg_together_rejected(periodic_field):
    f, L = periodic_field
    with pytest.raises(TypeError):
        build_multigraph_library(f, L, L, formulation="space_weak", m=6, p=8,
                                 space_ibp=True, use_ibp=True)


def test_unknown_formulation_rejected(periodic_field):
    f, L = periodic_field
    with pytest.raises(ValueError):
        build_pdefind_library(f, L, L, formulation="time_weak")


# ---------------------------------------------------------------------------
# Sweep CSV names
# ---------------------------------------------------------------------------


def test_normalize_sweep_row_translates_old_names():
    row = {"method": "weak_ibp", "weak_m": "8", "fp": "0"}
    assert normalize_sweep_row(row) is row
    assert row == {"method": "space_weak_ibp", "space_m": "8", "fp": "0"}
    nd = normalize_sweep_row({"use_ibp": "direct"})
    assert nd == {"space_ibp": "direct"}
    coll = normalize_sweep_row({"formulation": "weak_direct"})
    assert coll == {"formulation": "space_weak_direct"}
    new = {"method": "spacetime_weak_ibp", "space_m": "4"}
    assert normalize_sweep_row(dict(new)) == new


def test_every_method_has_a_label():
    for old in ("weak_direct", "weak_ibp"):
        assert canonical_method(old) in METHOD_MAP
    sweep = pytest.importorskip("sweep", reason="scripts/sweep.py not importable")
    assert set(sweep.METHODS) == set(METHOD_MAP)


def test_sweep_configs_expand_boundaries():
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    sweep = pytest.importorskip("sweep", reason="scripts/sweep.py not importable")
    both = ("valid", "periodic")
    cfgs = sweep.test_fn_configs("space_weak_ibp", None, [4, 8, 12], [8, 12],
                                 [4], [True], 40, boundary_modes=both)
    assert [(c["m"], c["space_boundary"]) for c in cfgs] == [
        (8, "valid"), (8, "periodic"), (12, "valid"), (12, "periodic")]
    cfgs = sweep.test_fn_configs("time_weak", None, [8], [8], [4, 8],
                                 [True, False], 40, boundary_modes=both)
    assert all("space_boundary" not in c for c in cfgs) and len(cfgs) == 4
    rule = dict(m=15, p=12, m_t=10, p_t=9)
    cfgs = sweep.test_fn_configs("spacetime_weak_direct", rule, [], [], [],
                                 [True], 40, boundary_modes=both)
    assert [c["space_boundary"] for c in cfgs] == ["valid", "periodic"]
    assert sweep.cfg_columns(cfgs[1], "mb")["space_boundary"] == "periodic"
    assert sweep.test_fn_configs("pointwise", rule, [], [], [], [True], 40,
                                 boundary_modes=both) == [{}]
