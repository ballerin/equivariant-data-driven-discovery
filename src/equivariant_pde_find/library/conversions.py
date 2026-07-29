"""
Change-of-basis conversions between PDE library representations.

This module provides two families of basis conversions:

1. **Multigraph ↔ Modified Multigraph** (matrix-based)
   Linear change-of-basis via precomputed block matrices encoding the IBP
   identities from the reference paper.

2. **Any library → Polynomial Jet** (analytical decomposition)
   Maps coefficient vectors from multigraph, modified multigraph, WSINDy,
   and PDE-FIND libraries into the universal polynomial-jet comparison space.
   The jet basis can represent all terms from every other library, making it
   the natural common space for cross-library evaluation.

We could simplify the code by rewriting the mod_multigraph → multigraph
conversion in terms of the multigraph → jet plus a matrix mapping. But we
keep it this way to keep the code more explicit and to add an equivalence
test to check that the two procedures give the same result.
"""

from functools import lru_cache

import numpy as np

from ..multigraph_operators import (
    get_mod_multigraph_operators,
    get_operators,
    _sym_Delta,
    _sym_grad_sq,
    _sym_bilap,
    _sym_grad_Delta_grad,
    _sym_hess_norm_sq,
    _sym_hess_grad_grad,
)
from ..jet_algebra import JetExpr, AXIS_LABELS, expr_to_dict


# =====================================================================================
# Part 1: Multigraph ↔ Modified Multigraph
# =====================================================================================
# Linear change-of-basis via precomputed block matrices.  For each (p, n)
# subspace where the two bases differ, a block matrix encodes the IBP
# identities.  Subspaces that are identical map through the identity.
# =====================================================================================

# -----------------------------------------------------------------------
# Block matrices mapping mod. multigraph → original (w_orig = M @ w_mod)
# -----------------------------------------------------------------------

_M_12 = np.array([[2.0, 0.0],
                  [2.0, 1.0]])

_M_13 = np.array([[3.0, 0.0],
                  [6.0, 1.0]])

_M_22 = np.array([[2.0, 0.0, 0.0, 0.0],
                  [8.0, 2.0, 0.0, 0.0],
                  [4.0, 2.0, 1.0, 0.0],
                  [2.0, 0.0, 0.0, 1.0]])

_M_23 = np.array([[3.0,  0.0, 0.0, 0.0, 0.0, 0.0],
                  [24.0, 2.0, 0.0, 0.0, 0.0, 0.0],
                  [12.0, 2.0, 1.0, 0.0, 0.0, 0.0],
                  [6.0,  0.0, 0.0, 1.0, 0.0, 0.0],
                  [24.0, 4.0, 0.0, 0.0, 1.0, 0.0],
                  [12.0, 1.0, 0.0, 0.0, 0.0, 1.0]])


# -----------------------------------------------------------------------
# Block specifications
# -----------------------------------------------------------------------

_BLOCKS = [
    (["Δf²", "‖∇f‖²"],
     ["f·Δf", "‖∇f‖²"],
     _M_12),

    (["Δf³", "f·‖∇f‖²"],
     ["f²·Δf", "f·‖∇f‖²"],
     _M_13),

    (["Δ²f²", "Δ‖∇f‖²", "‖∇²f‖²", "Δf·Δf"],
     ["f·Δ²f", "⟨∇Δf,∇f⟩", "‖∇²f‖²", "Δf·Δf"],
     _M_22),

    (["Δ²f³", "Δ(f·‖∇f‖²)", "f·‖∇²f‖²", "f·Δf·Δf",
      "∇²f(∇f,∇f)", "‖∇f‖²·Δf"],
     ["f²·Δ²f", "f·⟨∇Δf,∇f⟩", "f·‖∇²f‖²", "f·Δf·Δf",
      "∇²f(∇f,∇f)", "‖∇f‖²·Δf"],
     _M_23),
]


# -----------------------------------------------------------------------
# Public API: multigraph ↔ mod_multigraph
# -----------------------------------------------------------------------

def mod_multigraph_to_multigraph_matrix(max_p=2, max_n=3):
    """Build the change-of-basis matrix M such that w_orig = M @ w_mod.

    For most (p, n) subspaces the two bases are identical, so the
    corresponding diagonal block of M is an identity matrix.  For the
    subspaces where the modified basis is genuinely different, the blocks
    are filled with the coefficients from the reference paper.

    Parameters
    ----------
    max_p, max_n : int - passed to get_operators / get_mod_multigraph_operators

    Returns
    -------
    M : ndarray (K, K) where K = len(get_operators(max_p, max_n))
    """
    descs_orig = [d for _, d in get_operators(max_p, max_n)]
    descs_mod = [d for _, d in get_mod_multigraph_operators(max_p, max_n)]
    K = len(descs_orig)
    M = np.eye(K)

    for mod_names, orig_names, block in _BLOCKS:
        mod_map = {d: i for i, d in enumerate(descs_mod)}
        orig_map = {d: i for i, d in enumerate(descs_orig)}
        try:
            idx_mod = [mod_map[n] for n in mod_names]
            idx_orig = [orig_map[n] for n in orig_names]
        except KeyError:
            # Block not present in this (max_p, max_n) configuration
            continue
        for i_local, i_full in enumerate(idx_orig):
            for j_local, j_full in enumerate(idx_mod):
                M[i_full, j_full] = block[i_local, j_local]

    return M


def multigraph_to_mod_multigraph_matrix(max_p=2, max_n=3):
    """Build M⁻¹ s.t. w_mod = M⁻¹ @ w_orig."""
    M = mod_multigraph_to_multigraph_matrix(max_p, max_n)
    return np.linalg.inv(M)


def convert_coefficients(w, from_basis, to_basis, max_p=2, max_n=3):
    """Convert a coefficient vector between multigraph bases.

    Parameters
    ----------
    w : ndarray (K,)
        Coefficient vector in ``from_basis``.
    from_basis : str
        One of ``'multigraph'`` or ``'mod_multigraph'``.
    to_basis : str
        One of ``'multigraph'`` or ``'mod_multigraph'``.
    max_p, max_n : int

    Returns
    -------
    w_out : ndarray (K,)
        Coefficient vector in ``to_basis``.
    """
    if from_basis == to_basis:
        return w.copy()
    if from_basis == "mod_multigraph" and to_basis == "multigraph":
        M = mod_multigraph_to_multigraph_matrix(max_p, max_n)
        return M @ w
    if from_basis == "multigraph" and to_basis == "mod_multigraph":
        M_inv = multigraph_to_mod_multigraph_matrix(max_p, max_n)
        return M_inv @ w
    raise ValueError(
        f"Unknown basis conversion: {from_basis} → {to_basis}. "
        f"Expected 'multigraph' or 'mod_multigraph'."
    )


# =====================================================================================
# Part 2: Library → Polynomial Jet mappings
# =====================================================================================
# The polynomial-jet library serves as the universal comparison space.
# Every operator from every other library can be analytically expanded
# as a linear combination of jet monomials.  The mappings below implement
# these expansions.
# =====================================================================================

# -----------------------------------------------------------------------
# Multigraph → Jet decomposition, for any ndim
# -----------------------------------------------------------------------
# Every entry is defined via the SAME symbolic jet_algebra formulas used by
# multigraph_operators.py's numeric evaluators (imported above), so the
# numeric operator and its jet decomposition are derived from one source
# and can never silently disagree, at any ndim.  At ndim=2 these reproduce
# the historical hand-written decomposition dicts exactly (see
# tests/test_conversions_nd.py for the regression check).
# -----------------------------------------------------------------------


@lru_cache(maxsize=None)
def _multigraph_to_jet_decomposition(ndim):
    """{operator_desc: [(jet_desc, coeff), ...]} for the original
    multigraph basis, at the given ndim."""
    F = JetExpr.var(ndim)
    Delta = _sym_Delta(ndim)
    grad_sq = _sym_grad_sq(ndim)
    bilap = _sym_bilap(ndim)
    grad_Delta_grad = _sym_grad_Delta_grad(ndim)
    hess_norm_sq = _sym_hess_norm_sq(ndim)
    hess_grad_grad = _sym_hess_grad_grad(ndim)

    exprs = {
        "1": JetExpr.constant(ndim, 1.0),
        "f": F,
        "f²": F * F,
        "f³": F * F * F,
        "Δf": Delta,
        "‖∇f‖²": grad_sq,
        "f·‖∇f‖²": F * grad_sq,
        "Δ²f": bilap,
        "Δf·Δf": Delta * Delta,
        "‖∇²f‖²": hess_norm_sq,
        "f·Δf·Δf": F * Delta * Delta,
        "f·‖∇²f‖²": F * hess_norm_sq,
        "∇²f(∇f,∇f)": hess_grad_grad,
        "‖∇f‖²·Δf": grad_sq * Delta,
        # multigraph-specific
        "f·Δf": F * Delta,
        "f²·Δf": F * F * Delta,
        "f·Δ²f": F * bilap,
        "⟨∇Δf,∇f⟩": grad_Delta_grad,
        "f²·Δ²f": F * F * bilap,
        "f·⟨∇Δf,∇f⟩": F * grad_Delta_grad,
    }
    return {desc: list(expr_to_dict(e).items()) for desc, e in exprs.items()}


@lru_cache(maxsize=None)
def _mod_multigraph_to_jet_decomposition(ndim):
    """{operator_desc: [(jet_desc, coeff), ...]} for the modified
    (IBP-friendly) multigraph basis, at the given ndim."""
    from ..multigraph_operators import (
        _sym_Delta_f2, _sym_Delta_f3, _sym_Delta_grad_sq,
        _sym_Delta_f_grad_sq, _sym_Delta2_f2, _sym_Delta2_f3,
    )

    F = JetExpr.var(ndim)
    grad_sq = _sym_grad_sq(ndim)
    hess_norm_sq = _sym_hess_norm_sq(ndim)
    hess_grad_grad = _sym_hess_grad_grad(ndim)

    exprs = {
        "1": JetExpr.constant(ndim, 1.0),
        "f": F,
        "f²": F * F,
        "f³": F * F * F,
        "Δf": _sym_Delta(ndim),
        "‖∇f‖²": grad_sq,
        "f·‖∇f‖²": F * grad_sq,
        "Δ²f": _sym_bilap(ndim),
        "Δf·Δf": _sym_Delta(ndim) * _sym_Delta(ndim),
        "‖∇²f‖²": hess_norm_sq,
        "f·Δf·Δf": F * _sym_Delta(ndim) * _sym_Delta(ndim),
        "f·‖∇²f‖²": F * hess_norm_sq,
        "∇²f(∇f,∇f)": hess_grad_grad,
        "‖∇f‖²·Δf": grad_sq * _sym_Delta(ndim),
        # modified-basis-specific replacements
        "Δf²": _sym_Delta_f2(ndim),
        "Δf³": _sym_Delta_f3(ndim),
        "Δ‖∇f‖²": _sym_Delta_grad_sq(ndim),
        "Δ(f·‖∇f‖²)": _sym_Delta_f_grad_sq(ndim),
        "Δ²f²": _sym_Delta2_f2(ndim),
        "Δ²f³": _sym_Delta2_f3(ndim),
    }
    return {desc: list(expr_to_dict(e).items()) for desc, e in exprs.items()}


# Backward-compatible module-level dicts (ndim=2), for existing callers that
# import these names directly as plain dicts (evaluation.py, scripts, tests).
_MULTIGRAPH_TO_JET_DECOMPOSITION = _multigraph_to_jet_decomposition(2)
_MOD_MULTIGRAPH_TO_JET_DECOMPOSITION = _mod_multigraph_to_jet_decomposition(2)


# -----------------------------------------------------------------------
# WSINDy → Jet: via jet_algebra symbolic differentiation (any ndim)
# -----------------------------------------------------------------------
# Each WSINDy column is exactly ∂^multi_index(f^q); jet_algebra's .deriv()
# already implements the full Leibniz product rule symbolically, so we
# reuse it here instead of the hand-rolled combinatorial expansion the
# 2-D-only version used to need.
# -----------------------------------------------------------------------

_SUP_REV = {"²": 2, "³": 3, "⁴": 4, "⁵": 5, "⁶": 6, "⁷": 7, "⁸": 8, "⁹": 9}


def _axis_str_to_multi_index(axis_str, ndim):
    mi = [0] * ndim
    for ch in axis_str:
        mi[AXIS_LABELS.index(ch)] += 1
    return tuple(mi)


def _parse_wsindy_desc(desc, ndim):
    """Parse a WSINDy column description into (multi_index, q).

    '1' -> (zeros, 0);  'f' -> (zeros, 1);  'f²' -> (zeros, 2)
    'f_{xx}' -> (multi_index, 1);  '(f²)_{xy}' -> (multi_index, 2)
    """
    zero_mi = (0,) * ndim
    if desc == "1":
        return zero_mi, 0
    if desc.startswith("(f"):
        close_paren = desc.index(")")
        q = _SUP_REV.get(desc[2:close_paren], 2)
        brace_start = desc.index("{", close_paren)
        axis_str = desc[brace_start + 1: -1] if desc.endswith("}") else desc[brace_start + 1:]
        return _axis_str_to_multi_index(axis_str, ndim), q
    if desc.startswith("f_{"):
        axis_str = desc[3:-1] if desc.endswith("}") else desc[3:]
        return _axis_str_to_multi_index(axis_str, ndim), 1
    if desc.startswith("f"):
        q = _SUP_REV.get(desc[1:], 1) if desc != "f" else 1
        return zero_mi, q
    raise ValueError(f"Cannot parse WSINDy description: {desc!r}")


def build_wsindy_to_jet_map(desc_ws, desc_jet, ndim=2):
    """Build a matrix M such that w_jet = M @ w_wsindy, for any ndim
    (default 2, matching the paper's 2-D domain).

    Each WSINDy column ∂^α(f^q) is analytically expanded to jet monomials
    via jet_algebra's symbolic differentiation (full Leibniz product rule).
    """
    M = np.zeros((len(desc_jet), len(desc_ws)))
    for j, desc in enumerate(desc_ws):
        multi_index, q = _parse_wsindy_desc(desc, ndim)
        if q == 0:
            if "1" in desc_jet:
                M[desc_jet.index("1"), j] = 1.0
            continue

        expr = JetExpr.var(ndim) ** q
        for axis, a in enumerate(multi_index):
            for _ in range(a):
                expr = expr.deriv(axis)

        for jet_desc, coeff in expr_to_dict(expr).items():
            try:
                i = desc_jet.index(jet_desc)
            except ValueError:
                continue
            M[i, j] += coeff

    return M


# -----------------------------------------------------------------------
# Public API: Library → Jet mappers
# -----------------------------------------------------------------------

def map_wsindy_to_jet(w_ws, M_ws_to_jet):
    """Map WSINDy coefficients to jet coefficients."""
    return M_ws_to_jet @ w_ws


def map_standard_to_jet(w_std, desc_std, desc_jet):
    """Map standard-library coefficients to jet coefficients.

    Standard columns are of the form f^p · ∂_x^a ∂_y^b f which are
    individual jet monomials; mapping is 1:1 index lookup.
    Columns beyond the jet library's degree/order limits are silently
    skipped (they cannot contribute to the evaluation).
    """
    w_jet = np.zeros(len(desc_jet))
    for j, desc in enumerate(desc_std):
        if abs(w_std[j]) < 1e-14:
            continue
        if desc in desc_jet:
            w_jet[desc_jet.index(desc)] += w_std[j]
    return w_jet


def map_multigraph_to_jet(w_eq, desc_eq, desc_jet, ndim=2):
    """Map multigraph coefficients to jet coefficients using the
    analytical decomposition, at the given ndim (default 2, matching the
    paper's 2-D domain).

    Raises KeyError if an operator name is not in the decomposition."""
    decomposition = _multigraph_to_jet_decomposition(ndim)
    w_jet = np.zeros(len(desc_jet))
    for j, (name, coeff) in enumerate(zip(desc_eq, w_eq)):
        if abs(coeff) < 1e-14:
            continue
        if name not in decomposition:
            raise KeyError(
                f"Multigraph operator {name!r} has no jet decomposition. "
                f"Available: {list(decomposition)}"
            )
        for jet_name, mult in decomposition[name]:
            w_jet[desc_jet.index(jet_name)] += coeff * mult
    return w_jet


def map_mod_multigraph_to_jet(w_mod, desc_mod, desc_jet, ndim=2):
    """Map modified multigraph coefficients to jet coefficients using the
    analytical decomposition, at the given ndim (default 2, matching the
    paper's 2-D domain).

    Raises KeyError if an operator name is not in the decomposition."""
    decomposition = _mod_multigraph_to_jet_decomposition(ndim)
    w_jet = np.zeros(len(desc_jet))
    for j, (name, coeff) in enumerate(zip(desc_mod, w_mod)):
        if abs(coeff) < 1e-14:
            continue
        if name not in decomposition:
            raise KeyError(
                f"Mod. multigraph operator {name!r} has no jet decomposition. "
                f"Available: {list(decomposition)}"
            )
        for jet_name, mult in decomposition[name]:
            w_jet[desc_jet.index(jet_name)] += coeff * mult
    return w_jet
