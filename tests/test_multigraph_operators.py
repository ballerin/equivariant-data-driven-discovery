"""
Tests for `multigraph_operators.py`.

These tests pin the canonical operator ordering documented in README.md.
Any change to the operator set or ordering must be reflected here.
"""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.multigraph_operators import get_operators  # noqa: E402


# ---------------------------------------------------------------------------
# Canonical operator list for get_operators(max_p=2, max_n=3)
# Order: primary by p (edge count), secondary by k (polynomial degree),
#        tertiary by description string (Unicode lexicographic).
# This table is the source of truth; see README.md for the full index table.
# ---------------------------------------------------------------------------
EXPECTED_DESCRIPTIONS_2_3 = [
    # p=0 (no differential part)
    "1",         # k=0
    "f",         # k=1
    "f²",        # k=2
    "f³",        # k=3
    # p=1 (order-2 operators)
    "Δf",        # k=1
    "f·Δf",      # k=2  ('f' < 'Δ' in Unicode)
    "‖∇f‖²",     # k=2
    "f²·Δf",     # k=3  ('f²·' < 'f·‖' because '²' < '·' in Unicode)
    "f·‖∇f‖²",   # k=3
    # p=2 (order-4 operators)
    "Δ²f",       # k=1
    "f·Δ²f",     # k=2
    "Δf·Δf",     # k=2
    "‖∇²f‖²",    # k=2
    "⟨∇Δf,∇f⟩",  # k=2
    "f²·Δ²f",    # k=3
    "f·Δf·Δf",   # k=3
    "f·‖∇²f‖²",  # k=3
    "f·⟨∇Δf,∇f⟩", # k=3
    "‖∇f‖²·Δf",  # k=3
    "∇²f(∇f,∇f)", # k=3
]


def test_operator_count_2_3():
    """get_operators(2, 3) should return exactly 20 operators."""
    ops = get_operators(2, 3)
    assert len(ops) == 20, f"Expected 20 operators, got {len(ops)}"


def test_operator_descriptions_2_3():
    """Operator descriptions must match the canonical ordered list exactly."""
    ops = get_operators(2, 3)
    actual = [desc for _, desc in ops]
    assert actual == EXPECTED_DESCRIPTIONS_2_3, (
        "Operator ordering has changed.  Diff:\n"
        + "\n".join(
            f"  [{i}] expected {repr(e)}, got {repr(a)}"
            for i, (e, a) in enumerate(zip(EXPECTED_DESCRIPTIONS_2_3, actual))
            if e != a
        )
        + (
            f"\n  (lengths differ: expected {len(EXPECTED_DESCRIPTIONS_2_3)}, "
            f"got {len(actual)})"
            if len(EXPECTED_DESCRIPTIONS_2_3) != len(actual)
            else ""
        )
    )


def test_operator_ordering_rule_2_3():
    """Ordering must satisfy the three-level rule: (p, k, desc) lexicographic."""
    # Re-derive (p, k) from operator structure by calling each op on a
    # test field and comparing with the enumeration internals.
    # Simpler approach: check that the sort order is non-decreasing on (p, k, desc).
    ops = get_operators(2, 3)
    descs = [desc for _, desc in ops]

    # We need (p, k) per operator; derive them from EXPECTED_DESCRIPTIONS_2_3
    # cross-referenced with the README table.
    p_k = [
        (0, 0), (0, 1), (0, 2), (0, 3),
        (1, 1), (1, 2), (1, 2), (1, 3), (1, 3),
        (2, 1), (2, 2), (2, 2), (2, 2), (2, 2),
        (2, 3), (2, 3), (2, 3), (2, 3), (2, 3), (2, 3),
    ]
    keys = [(p, k, d) for (p, k), d in zip(p_k, descs)]
    assert keys == sorted(keys), "Operators are not in (p, k, desc) lexicographic order."


@pytest.mark.parametrize("max_p,max_n", [(0, 3), (1, 2), (1, 3), (2, 3)])
def test_operator_ordering_is_sorted(max_p, max_n):
    """For any (max_p, max_n), the returned list must be sorted by desc strings
    within each (p, k) bucket — verified by re-sorting and comparing."""
    from equivariant_pde_find.multigraph_operators import _all_base_operators

    raw = _all_base_operators(max_p, max_n)
    # raw entries: (total_p, k, func, desc)
    ops = get_operators(max_p, max_n)
    descs = [desc for _, desc in ops]

    # Group by (p, k) in the order they appear and check each group is sorted.
    from itertools import groupby

    # Attach expected (p, k) from raw results (already sorted by (p, k, desc))
    expected_keys = [(p, k, desc) for p, k, _, desc in
                     sorted(raw, key=lambda x: (x[0], x[1], x[3]))]
    actual_keys = [(p, k, d) for (p, k), d in
                   zip([(pk[0], pk[1]) for pk in expected_keys], descs)]
    assert expected_keys == actual_keys, (
        "Within-group ordering is not lexicographic by description."
    )


def test_all_operators_are_callable():
    """Every returned operator must be callable with signature (f, Lx, Ly)."""
    N = 16
    Lx = Ly = 2 * np.pi
    x = np.linspace(0, Lx, N, endpoint=False)
    y = np.linspace(0, Ly, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="ij")
    f = np.sin(X) * np.cos(Y) + 0.5

    ops = get_operators(2, 3)
    for func, desc in ops:
        result = func(f, Lx, Ly)
        assert result.shape == f.shape, (
            f"Operator '{desc}' returned shape {result.shape}, expected {f.shape}"
        )
        assert np.isfinite(result).all(), (
            f"Operator '{desc}' returned non-finite values."
        )


def test_operator_count_1_2():
    """get_operators(1, 2) should return 6 operators."""
    ops = get_operators(1, 2)
    assert len(ops) == 6


def test_operator_descriptions_1_2():
    """Operator descriptions for get_operators(1, 2) in canonical order."""
    expected = ["1", "f", "f²", "Δf", "f·Δf", "‖∇f‖²"]
    ops = get_operators(1, 2)
    actual = [desc for _, desc in ops]
    assert actual == expected


# =============================================================================
# 3-D operator tests
# =============================================================================


class TestMultigraphOperators3D:
    """3-D operator evaluation correctness."""

    def test_all_operators_callable_3d(self):
        """Every operator accepts (f_3d, Lx, Ly, Lz) and returns correct shape."""
        ops = get_operators(2, 3)
        N = 16
        L = 2 * np.pi
        x = np.linspace(0, L, N, endpoint=False)
        X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
        f = np.sin(X) * np.cos(Y) * np.cos(Z) + 0.5

        for func, desc in ops:
            result = func(f, L, L, L)
            assert result.shape == f.shape, (
                f"Operator '{desc}' returned shape {result.shape}, expected {f.shape}"
            )
            assert np.isfinite(result).all(), (
                f"Operator '{desc}' returned non-finite values."
            )

    def test_grad_sq_3d_vs_manual(self):
        """_op_grad_sq in 3D = fx² + fy² + fz²."""
        from equivariant_pde_find.multigraph_operators import _op_grad_sq
        from equivariant_pde_find import derivatives
        N = 16
        L = 2 * np.pi
        x = np.linspace(0, L, N, endpoint=False)
        X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
        f = np.sin(2 * X) * np.cos(3 * Y) * np.cos(Z)

        result = _op_grad_sq(f, L, L, L)
        fx, fy, fz = derivatives.grad(f, L, L, L)
        expected = fx**2 + fy**2 + fz**2
        np.testing.assert_allclose(result, expected, rtol=1e-12, atol=1e-12)

    def test_grad_Delta_grad_3d_vs_manual(self):
        """_op_grad_Delta_grad in 3D = Dfx·fx + Dfy·fy + Dfz·fz."""
        from equivariant_pde_find.multigraph_operators import _op_grad_Delta_grad
        from equivariant_pde_find import derivatives
        N = 16
        L = 2 * np.pi
        x = np.linspace(0, L, N, endpoint=False)
        X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
        f = np.sin(2 * X) * np.cos(3 * Y) * np.cos(Z)

        result = _op_grad_Delta_grad(f, L, L, L)
        Df = derivatives.laplacian(f, L, L, L)
        Dfx, Dfy, Dfz = derivatives.grad(Df, L, L, L)
        fx, fy, fz = derivatives.grad(f, L, L, L)
        expected = Dfx * fx + Dfy * fy + Dfz * fz
        np.testing.assert_allclose(result, expected, rtol=1e-11, atol=1e-11)

    def test_hess_norm_sq_3d_vs_manual(self):
        """_op_hess_norm_sq in 3D = Σ f_{ij}²."""
        from equivariant_pde_find.multigraph_operators import _op_hess_norm_sq
        from equivariant_pde_find import derivatives
        N = 16
        L = 2 * np.pi
        x = np.linspace(0, L, N, endpoint=False)
        X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
        f = np.sin(2 * X) * np.cos(3 * Y) * np.cos(Z)

        result = _op_hess_norm_sq(f, L, L, L)
        fxx, fxy, fxz, fyy, fyz, fzz = derivatives.hessian(f, L, L, L)
        expected = fxx**2 + fyy**2 + fzz**2 + 2 * (fxy**2 + fxz**2 + fyz**2)
        np.testing.assert_allclose(result, expected, rtol=1e-11, atol=1e-11)

    def test_hess_grad_grad_3d_vs_manual(self):
        """_op_hess_grad_grad in 3D = Σ f_{ij} f_i f_j."""
        from equivariant_pde_find.multigraph_operators import _op_hess_grad_grad
        from equivariant_pde_find import derivatives
        N = 16
        L = 2 * np.pi
        x = np.linspace(0, L, N, endpoint=False)
        X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
        f = np.sin(2 * X) * np.cos(3 * Y) * np.cos(Z)

        result = _op_hess_grad_grad(f, L, L, L)
        fx, fy, fz = derivatives.grad(f, L, L, L)
        fxx, fxy, fxz, fyy, fyz, fzz = derivatives.hessian(f, L, L, L)
        expected = (fxx * fx**2 + fyy * fy**2 + fzz * fz**2
                    + 2 * (fxy * fx * fy + fxz * fx * fz + fyz * fy * fz))
        np.testing.assert_allclose(result, expected, rtol=1e-11, atol=1e-11)

    def test_2d_backward_compat(self):
        """Calling operators without Lz on 2D field gives same results as before."""
        from equivariant_pde_find.multigraph_operators import (
            _op_grad_sq, _op_hess_norm_sq, _op_hess_grad_grad, _op_grad_Delta_grad,
        )
        N = 16
        L = 2 * np.pi
        x = np.linspace(0, L, N, endpoint=False)
        X, Y = np.meshgrid(x, x, indexing="ij")
        f = np.sin(2 * X) * np.cos(3 * Y)

        # These should work without Lz (old calling convention)
        r1 = _op_grad_sq(f, L, L)
        r2 = _op_hess_norm_sq(f, L, L)
        r3 = _op_hess_grad_grad(f, L, L)
        r4 = _op_grad_Delta_grad(f, L, L)
        assert r1.shape == f.shape
        assert r2.shape == f.shape
        assert r3.shape == f.shape
        assert r4.shape == f.shape
