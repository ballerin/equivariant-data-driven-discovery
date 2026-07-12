"""Tests for the Modified Multigraph (IBP-friendly) operator basis."""

import math
import numpy as np
import pytest

from equivariant_pde_find.multigraph_operators import (
    get_operators,
    get_mod_multigraph_operators,
)
from equivariant_pde_find.library import (
    build_multigraph_library,
    build_mod_multigraph_library,
)
from equivariant_pde_find.basis_transforms import (
    mod_multigraph_to_multigraph_matrix,
    multigraph_to_mod_multigraph_matrix,
    convert_coefficients,
)
from equivariant_pde_find.equation_configs import (
    true_allen_cahn,
    true_swift_hohenberg,
    true_ks,
    true_cahn_hilliard,
    true_pme,
    true_allen_cahn_mod_multigraph,
    true_swift_hohenberg_mod_multigraph,
    true_ks_mod_multigraph,
    true_cahn_hilliard_mod_multigraph,
    true_pme_mod_multigraph,
)


# ── Fixtures ────────────────────────────────────────────────────────────


@pytest.fixture
def test_field():
    N, L = 32, 2 * math.pi
    rng = np.random.default_rng(42)
    f = rng.standard_normal((N, N))
    return f, L


# ── Operator count ──────────────────────────────────────────────────────


def test_operator_count():
    """get_mod_multigraph_operators must return exactly 20 operators."""
    ops = get_mod_multigraph_operators(max_p=2, max_n=3)
    assert len(ops) == 20, f"Expected 20 operators, got {len(ops)}"


def test_subset_counts():
    """Test that max_p/max_n filtering works."""
    assert len(get_mod_multigraph_operators(max_p=2, max_n=2)) == 11
    assert len(get_mod_multigraph_operators(max_p=1, max_n=3)) == 9
    assert len(get_mod_multigraph_operators(max_p=1, max_n=1)) == 3


# ── Basis equivalence ──────────────────────────────────────────────────


def test_basis_equivalence(test_field):
    """Both bases must span the same 20-dimensional column space."""
    f, L = test_field
    T_orig, _ = build_multigraph_library(f, L, L)
    T_mod, _ = build_mod_multigraph_library(f, L, L)
    combined = np.column_stack([T_orig, T_mod])
    U, S, Vt = np.linalg.svd(combined, full_matrices=False)
    rank = np.sum(S > 1e-10 * S[0])
    assert rank == 20, f"Combined rank is {rank}, expected 20"


# ── Change-of-basis ────────────────────────────────────────────────────


def test_change_of_basis_random(test_field):
    """M @ w_mod must equal w_orig for random coefficient vectors."""
    f, L = test_field
    T_orig, _ = build_multigraph_library(f, L, L)
    T_mod, _ = build_mod_multigraph_library(f, L, L)
    M = mod_multigraph_to_multigraph_matrix()
    rng = np.random.default_rng(99)
    w_mod = rng.standard_normal(T_mod.shape[1])
    w_orig = M @ w_mod
    pred_mod = T_mod @ w_mod
    pred_orig = T_orig @ w_orig
    rel_err = np.max(np.abs(pred_mod - pred_orig)) / np.max(np.abs(pred_orig))
    assert rel_err < 1e-10, f"Change-of-basis error {rel_err:.2e}"


def test_change_of_basis_inverse():
    """M and M⁻¹ must be genuine inverses."""
    M = mod_multigraph_to_multigraph_matrix()
    M_inv = multigraph_to_mod_multigraph_matrix()
    assert np.allclose(M @ M_inv, np.eye(20), atol=1e-12)
    assert np.allclose(M_inv @ M, np.eye(20), atol=1e-12)


@pytest.mark.parametrize("name, true_orig, true_mod", [
    ("AC", true_allen_cahn, true_allen_cahn_mod_multigraph),
    ("SH", true_swift_hohenberg, true_swift_hohenberg_mod_multigraph),
    ("KS", true_ks, true_ks_mod_multigraph),
    ("CH", true_cahn_hilliard, true_cahn_hilliard_mod_multigraph),
    ("PME", true_pme, true_pme_mod_multigraph),
])
def test_ground_truth_consistency(test_field, name, true_orig, true_mod):
    """Ground truth in both bases must produce the same PDE prediction."""
    f, L = test_field
    T_orig, d_orig = build_multigraph_library(f, L, L)
    T_mod, d_mod = build_mod_multigraph_library(f, L, L)
    w_orig = true_orig(d_orig)
    w_mod = true_mod(d_mod)
    pred_orig = T_orig @ w_orig
    pred_mod = T_mod @ w_mod
    scale = max(np.max(np.abs(pred_orig)), 1e-10)
    rel_err = np.max(np.abs(pred_mod - pred_orig)) / scale
    assert rel_err < 1e-10, f"{name}: ground truth mismatch ({rel_err:.2e})"


def test_convert_coefficients():
    """convert_coefficients must round-trip correctly."""
    rng = np.random.default_rng(7)
    w_orig = rng.standard_normal(20)
    w_mod = convert_coefficients(w_orig, "multigraph", "mod_multigraph")
    w_back = convert_coefficients(w_mod, "mod_multigraph", "multigraph")
    assert np.allclose(w_back, w_orig, atol=1e-12)


# ── Ordering ────────────────────────────────────────────────────────────


def test_ordering():
    """Operators must be sorted by (total_p, k, description)."""
    ops = get_mod_multigraph_operators()
    descs = [d for _, d in ops]
    from equivariant_pde_find.multigraph_operators import _MOD_MULTIGRAPH_SPEC
    pk = {d: (p, k) for _, d, p, k in _MOD_MULTIGRAPH_SPEC}
    for i in range(len(descs) - 1):
        p0, k0 = pk[descs[i]]
        p1, k1 = pk[descs[i + 1]]
        assert (p0, k0, descs[i]) <= (p1, k1, descs[i + 1]), \
            f"Order violation at {i}: {descs[i]} > {descs[i+1]}"


# ── Equivariance ───────────────────────────────────────────────────────


@pytest.fixture
def smooth_field():
    """A smooth, band-limited periodic field for equivariance tests."""
    N, L = 64, 2 * math.pi
    x = np.linspace(0, L, N, endpoint=False)
    y = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="ij")
    return np.cos(2 * X + 3 * Y) + 0.5 * np.sin(X - 2 * Y), L


def test_mod_multigraph_rotation_equivariance(smooth_field):
    """All mod. multigraph operators must be rotation-equivariant."""
    f, L = smooth_field
    ops = get_mod_multigraph_operators()
    f_rot = np.rot90(f)
    for func, desc in ops:
        v = func(f, L, L)
        v_rot = func(f_rot, L, L)
        np.testing.assert_allclose(
            np.rot90(v), v_rot, rtol=5e-10, atol=5e-10,
            err_msg=f"Rotation equivariance failed for {desc}"
        )


def test_mod_multigraph_reflection_equivariance(smooth_field):
    """All mod. multigraph operators must be reflection-equivariant."""
    f, L = smooth_field
    ops = get_mod_multigraph_operators()
    f_ref = np.fliplr(f)
    for func, desc in ops:
        v = func(f, L, L)
        v_ref = func(f_ref, L, L)
        np.testing.assert_allclose(
            np.fliplr(v), v_ref, rtol=5e-10, atol=5e-10,
            err_msg=f"Reflection equivariance failed for {desc}"
        )


# ── Jet decomposition ─────────────────────────────────────────────────


def test_jet_decomposition():
    """All 20 mod. multigraph operators must have a jet decomposition."""
    from equivariant_pde_find.library_mapping import _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION
    ops = get_mod_multigraph_operators()
    for _, desc in ops:
        assert desc in _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION, \
            f"Missing jet decomposition for {desc}"


def test_jet_decomposition_naming_consistency():
    """Every jet name in the mod. multigraph decomposition must exist
    in the polynomial-jet library produced by build_polyjet_library."""
    from equivariant_pde_find.library_mapping import _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION
    from equivariant_pde_find.library import build_polyjet_library
    N, L = 32, 2 * math.pi
    rng = np.random.default_rng(99)
    f = rng.standard_normal((N, N))
    _, desc_jet = build_polyjet_library(f, L, L, max_order=4, max_degree=3)
    desc_jet_set = set(desc_jet)
    for mg_name, comps in _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION.items():
        for jet_name, _mult in comps:
            assert jet_name in desc_jet_set, \
                f"Jet name {jet_name!r} (from {mg_name}) not in polynomial-jet library"


def test_map_mod_multigraph_to_jet_known_coefficients():
    """Mapping known Swift-Hohenberg coefficients (f, Δf, Δ²f, f³)
    from mod. multigraph to jet basis must give the correct 7 jet terms."""
    from equivariant_pde_find.library_mapping import map_mod_multigraph_to_jet
    from equivariant_pde_find.library import build_polyjet_library

    N, L = 32, 2 * math.pi
    rng = np.random.default_rng(42)
    f = rng.standard_normal((N, N))
    _, desc_mod = build_mod_multigraph_library(f, L, L, max_p=2, max_n=3)
    _, desc_jet = build_polyjet_library(f, L, L, max_order=4, max_degree=3)

    w_mod = true_swift_hohenberg_mod_multigraph(desc_mod)
    w_jet = map_mod_multigraph_to_jet(w_mod, desc_mod, desc_jet)

    # SH truth: f = -0.7, Δf = -2.0, Δ²f = -1.0, f³ = -1.0
    # Jet: f, f_xx, f_yy, f_xxxx, f_xxyy, f_yyyy, f³
    expected = {
        "f": -0.7,
        "f_xx": -2.0,
        "f_yy": -2.0,
        "f_xxxx": -1.0,
        "f_xxyy": -2.0,
        "f_yyyy": -1.0,
        "f³": -1.0,
    }
    for jet_name, expected_val in expected.items():
        j = desc_jet.index(jet_name)
        assert abs(w_jet[j] - expected_val) < 1e-12, \
            f"Jet {jet_name}: expected {expected_val}, got {w_jet[j]}"

    # All other jet terms must be zero
    nz_expected = len(expected)
    nz_actual = int(np.sum(np.abs(w_jet) > 1e-12))
    assert nz_actual == nz_expected, \
        f"Expected {nz_expected} nonzero jet terms, got {nz_actual}"


def test_evaluate_unified_handles_mod_multigraph():
    """evaluate_unified must handle library='mod_multigraph' without error
    and produce consistent fp/fn/E against known ground truth."""
    from equivariant_pde_find.evaluation import evaluate_unified
    from equivariant_pde_find.library_mapping import map_mod_multigraph_to_jet
    from equivariant_pde_find.library import build_polyjet_library

    N, L = 32, 2 * math.pi
    rng = np.random.default_rng(42)
    f = rng.standard_normal((N, N))
    _, desc_mod = build_mod_multigraph_library(f, L, L, max_p=2, max_n=3)
    _, desc_jet = build_polyjet_library(f, L, L, max_order=4, max_degree=3)

    # Perfect recovery: pass the true coefficients
    w_perfect = true_swift_hohenberg_mod_multigraph(desc_mod)
    true_jet = map_mod_multigraph_to_jet(w_perfect, desc_mod, desc_jet)

    # Build a dummy Theta_jet and Ft for evaluation
    n = 100
    Theta_jet = rng.standard_normal((n, len(desc_jet)))
    Ft = Theta_jet @ true_jet  # perfect fit

    result = evaluate_unified(
        w_perfect, "mod_multigraph", desc_mod, desc_jet,
        Theta_jet, Ft, true_jet,
    )
    assert result["fp"] == 0, f"Perfect recovery should have fp=0, got {result['fp']}"
    assert result["fn"] == 0, f"Perfect recovery should have fn=0, got {result['fn']}"
    assert result["exact_recovery"] == 1
    assert result["E_coeff"] < 1e-10, f"E should be ~0, got {result['E']}"

    # Partial recovery: drop Δ²f (index 9 in mod multigraph)
    w_partial = w_perfect.copy()
    w_partial[9] = 0.0  # zero out Δ²f
    result2 = evaluate_unified(
        w_partial, "mod_multigraph", desc_mod, desc_jet,
        Theta_jet, Ft, true_jet,
    )
    assert result2["fn"] >= 1, "Dropping Δ²f should give fn >= 1 in jet space"


def test_mod_vs_orig_condition_number():
    """The modified multigraph library condition number must be finite
    and positive, and the subspace projection residual (QR-based) must
    confirm both bases span the same column space."""
    N, L = 64, 8 * math.pi  # realistic grid
    rng = np.random.default_rng(123)
    f = rng.standard_normal((N, N))

    T_orig, _ = build_multigraph_library(f, L, L, max_p=2, max_n=3)
    T_mod, _ = build_mod_multigraph_library(f, L, L, max_p=2, max_n=3)

    cond_orig = float(np.linalg.cond(T_orig))
    cond_mod = float(np.linalg.cond(T_mod))

    # Both must be finite and positive
    assert 1.0 <= cond_orig < 1e20, f"cond_orig out of range: {cond_orig:.2e}"
    assert 1.0 <= cond_mod < 1e20, f"cond_mod out of range: {cond_mod:.2e}"

    # Subspace identity: QR projection residual must be tiny
    Q_orig, _ = np.linalg.qr(T_orig)
    Q_mod, _ = np.linalg.qr(T_mod)
    proj_residual = np.linalg.norm(Q_mod - Q_orig @ (Q_orig.T @ Q_mod))
    assert proj_residual < 1e-10, \
        f"Subspace projection residual too large: {proj_residual:.2e}"
