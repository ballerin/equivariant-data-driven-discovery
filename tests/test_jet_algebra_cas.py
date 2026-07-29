"""Independent computer-algebra cross-check of the jet decompositions.

Why this file exists
--------------------
``multigraph_operators._sym_*(ndim)`` is the *single* source for both the
numeric operator column (via ``_make_symbolic_op`` → ``_evaluate_jet_expr``)
and the analytic jet decomposition (via ``library/conversions.py``).  That
is deliberate — the two can never silently diverge — but it means every
other N-D test in this suite compares ``jet_algebra`` against itself.  If a
``_sym_*`` formula were wrong at ndim=4, both sides would be wrong in
exactly the same way and every existing test would still pass.  The only
non-circular pin is the frozen ndim=2 fixture in ``test_conversions_nd.py``.

This file supplies genuinely independent ground truth at ndim ≥ 3: every
operator is re-derived from its *mathematical definition* using **sympy**,
whose ``Derivative``/``expand`` machinery implements the product and chain
rules completely independently of ``jet_algebra.JetExpr``.  Nothing here
imports ``jet_algebra``; the monomial-naming helpers below are written from
scratch so a naming divergence shows up as a key mismatch rather than
silently agreeing.

Two independent symbolic engines producing identical decompositions for all
20 multigraph and 20 modified-multigraph operators, at ndim = 2, 3 and 4, is
the strongest correctness evidence available for the N-D generalization.
"""

import os
import sys
from collections import defaultdict
from fractions import Fraction

import pytest

sympy = pytest.importorskip("sympy")

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.library.conversions import (  # noqa: E402
    _multigraph_to_jet_decomposition,
    _mod_multigraph_to_jet_decomposition,
)


# ---------------------------------------------------------------------------
# Local (independent) monomial-naming helpers.
#
# These deliberately duplicate jet_algebra's naming convention instead of
# importing it: the point of this file is to be an outside witness.
# ---------------------------------------------------------------------------

_LABELS = ["x", "y", "z", "w", "v", "u", "t", "s"]
_SUPERSCRIPT = {2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}


def _factor_name(multi_index):
    """(2, 1, 0) -> 'f_xxy';  (0, 0, 0) -> 'f'."""
    if not any(multi_index):
        return "f"
    body = "".join(_LABELS[axis] * order
                   for axis, order in enumerate(multi_index) if order)
    return "f_" + body


def _monomial_name(monomial):
    """Sorted tuple of multi-indices -> 'f_x²·f_xx'.  Empty -> '1'."""
    if not monomial:
        return "1"
    counts = defaultdict(int)
    for factor in monomial:
        counts[factor] += 1
    parts = []
    for factor in sorted(counts, key=lambda mi: (sum(mi), mi)):
        exponent = counts[factor]
        name = _factor_name(factor)
        parts.append(name if exponent == 1 else name + _SUPERSCRIPT[exponent])
    return "·".join(parts)


# ---------------------------------------------------------------------------
# sympy → {monomial name: coefficient}
# ---------------------------------------------------------------------------


def _multi_index_of(term, coords):
    """Multi-index of a single sympy factor that is f or a Derivative of f."""
    ndim = len(coords)
    if isinstance(term, sympy.Derivative):
        mi = [0] * ndim
        for var in term.variables:  # repeated per differentiation order
            mi[coords.index(var)] += 1
        return tuple(mi)
    return (0,) * ndim  # bare f


def _sympy_to_jet_dict(expr, coords):
    """Expand a sympy expression in f and its derivatives into
    ``{monomial name: float coefficient}``.

    The expression must be a polynomial in f and its partial derivatives —
    true for every operator in both bases.
    """
    expr = sympy.expand(sympy.expand(expr).doit())
    out = defaultdict(Fraction)

    for term in sympy.Add.make_args(expr):
        coeff = Fraction(1)
        factors = []
        for factor in sympy.Mul.make_args(term):
            base, exponent = factor.as_base_exp()
            if base.is_number and exponent.is_number:
                coeff *= Fraction(int(sympy.Integer(factor)), 1) if factor.is_Integer \
                    else Fraction(*sympy.Rational(factor).as_numer_denom())
                continue
            if not exponent.is_Integer or exponent < 0:
                raise AssertionError(f"non-polynomial factor {factor!r}")
            factors.extend([_multi_index_of(base, coords)] * int(exponent))
        out[tuple(sorted(factors))] += coeff

    return {_monomial_name(mono): float(c) for mono, c in out.items() if c != 0}


# ---------------------------------------------------------------------------
# Operator definitions, written straight from the mathematics.
# ---------------------------------------------------------------------------


def _setup(ndim):
    coords = sympy.symbols(f"x0:{ndim}", real=True)
    f = sympy.Function("f")(*coords)
    return coords, f


def _laplacian(expr, coords):
    """Δg = Σ_i ∂²g/∂x_i²."""
    return sum(sympy.diff(expr, c, 2) for c in coords)


def _grad_sq(expr, coords):
    """‖∇g‖² = Σ_i (∂g/∂x_i)²."""
    return sum(sympy.diff(expr, c) ** 2 for c in coords)


def _multigraph_operators_sympy(ndim):
    """The 20 multigraph operators (max_p=2, max_n=3) as sympy expressions."""
    coords, f = _setup(ndim)
    lap = _laplacian(f, coords)
    bilap = _laplacian(lap, coords)
    grad_sq = _grad_sq(f, coords)
    # ⟨∇Δf, ∇f⟩ = Σ_i ∂_i(Δf) · ∂_i f
    grad_lap_grad = sum(
        sympy.diff(lap, c) * sympy.diff(f, c) for c in coords
    )
    # ‖∇²f‖²_F = Σ_ij (∂_i∂_j f)²
    hess_norm_sq = sum(
        sympy.diff(f, ci, cj) ** 2 for ci in coords for cj in coords
    )
    # ∇²f(∇f, ∇f) = Σ_ij (∂_i∂_j f) ∂_i f ∂_j f
    hess_grad_grad = sum(
        sympy.diff(f, ci, cj) * sympy.diff(f, ci) * sympy.diff(f, cj)
        for ci in coords for cj in coords
    )
    return coords, {
        "1": sympy.Integer(1),
        "f": f,
        "f²": f ** 2,
        "f³": f ** 3,
        "Δf": lap,
        "‖∇f‖²": grad_sq,
        "f·‖∇f‖²": f * grad_sq,
        "Δ²f": bilap,
        "Δf·Δf": lap * lap,
        "‖∇²f‖²": hess_norm_sq,
        "f·Δf·Δf": f * lap * lap,
        "f·‖∇²f‖²": f * hess_norm_sq,
        "∇²f(∇f,∇f)": hess_grad_grad,
        "‖∇f‖²·Δf": grad_sq * lap,
        "f·Δf": f * lap,
        "f²·Δf": f ** 2 * lap,
        "f·Δ²f": f * bilap,
        "⟨∇Δf,∇f⟩": grad_lap_grad,
        "f²·Δ²f": f ** 2 * bilap,
        "f·⟨∇Δf,∇f⟩": f * grad_lap_grad,
    }


def _mod_multigraph_operators_sympy(ndim):
    """The 20 modified (IBP-friendly) operators as sympy expressions.

    The first 14 coincide with the multigraph basis; the last 6 replace the
    ``f^q ·`` prefactored forms with total-divergence forms Δ(·) / Δ²(·).
    """
    coords, shared = _multigraph_operators_sympy(ndim)
    _, f = _setup(ndim)
    coords = sympy.symbols(f"x0:{ndim}", real=True)
    f = sympy.Function("f")(*coords)
    grad_sq = _grad_sq(f, coords)

    ops = {k: shared[k] for k in [
        "1", "f", "f²", "f³", "Δf", "‖∇f‖²", "f·‖∇f‖²", "Δ²f", "Δf·Δf",
        "‖∇²f‖²", "f·Δf·Δf", "f·‖∇²f‖²", "∇²f(∇f,∇f)", "‖∇f‖²·Δf",
    ]}
    ops["Δf²"] = _laplacian(f ** 2, coords)
    ops["Δf³"] = _laplacian(f ** 3, coords)
    ops["Δ‖∇f‖²"] = _laplacian(grad_sq, coords)
    ops["Δ(f·‖∇f‖²)"] = _laplacian(f * grad_sq, coords)
    ops["Δ²f²"] = _laplacian(_laplacian(f ** 2, coords), coords)
    ops["Δ²f³"] = _laplacian(_laplacian(f ** 3, coords), coords)
    return coords, ops


# ---------------------------------------------------------------------------
# The cross-checks
# ---------------------------------------------------------------------------


def _compare(expected_sympy, actual_pairs, label):
    actual = dict(actual_pairs)
    assert set(actual) == set(expected_sympy), (
        f"{label}: monomial sets differ.\n"
        f"  only in library: {sorted(set(actual) - set(expected_sympy))}\n"
        f"  only in sympy:   {sorted(set(expected_sympy) - set(actual))}"
    )
    for mono, coeff in expected_sympy.items():
        assert actual[mono] == pytest.approx(coeff, rel=0, abs=1e-9), (
            f"{label}: coefficient of {mono} is {actual[mono]}, sympy says {coeff}"
        )


class TestMultigraphDecompositionAgainstSympy:
    @pytest.mark.parametrize("ndim", [2, 3])
    def test_all_operators(self, ndim):
        coords, ops = _multigraph_operators_sympy(ndim)
        library = _multigraph_to_jet_decomposition(ndim)
        assert set(library) == set(ops), "operator name sets differ"
        for desc, expr in ops.items():
            _compare(_sympy_to_jet_dict(expr, coords), library[desc],
                     f"multigraph ndim={ndim} {desc}")

    @pytest.mark.slow
    def test_all_operators_ndim4(self):
        self.test_all_operators(4)


class TestModMultigraphDecompositionAgainstSympy:
    @pytest.mark.parametrize("ndim", [2, 3])
    def test_all_operators(self, ndim):
        coords, ops = _mod_multigraph_operators_sympy(ndim)
        library = _mod_multigraph_to_jet_decomposition(ndim)
        assert set(library) == set(ops), "operator name sets differ"
        for desc, expr in ops.items():
            _compare(_sympy_to_jet_dict(expr, coords), library[desc],
                     f"mod_multigraph ndim={ndim} {desc}")

    @pytest.mark.slow
    def test_all_operators_ndim4(self):
        self.test_all_operators(4)


class TestOracleIsActuallyDiscriminating:
    """Guard against a vacuous oracle: the comparison must reject a wrong
    decomposition, and the sympy side must produce non-trivial expansions."""

    def test_rejects_perturbed_coefficient(self):
        coords, ops = _multigraph_operators_sympy(3)
        truth = _sympy_to_jet_dict(ops["Δ²f"], coords)
        wrong = list(_multigraph_to_jet_decomposition(3)["Δ²f"])
        wrong[0] = (wrong[0][0], wrong[0][1] + 1.0)
        with pytest.raises(AssertionError):
            _compare(truth, wrong, "perturbed")

    def test_rejects_missing_monomial(self):
        coords, ops = _multigraph_operators_sympy(3)
        truth = _sympy_to_jet_dict(ops["‖∇²f‖²"], coords)
        wrong = list(_multigraph_to_jet_decomposition(3)["‖∇²f‖²"])[:-1]
        with pytest.raises(AssertionError):
            _compare(truth, wrong, "truncated")

    def test_expansions_are_nontrivial(self):
        """The hardest operators really do expand to many monomials, so a
        passing comparison is meaningful rather than an empty-set match."""
        coords, ops = _mod_multigraph_operators_sympy(3)
        assert len(_sympy_to_jet_dict(ops["Δ²f³"], coords)) >= 15
        assert len(_sympy_to_jet_dict(ops["Δ(f·‖∇f‖²)"], coords)) >= 8
