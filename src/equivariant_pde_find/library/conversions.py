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

import numpy as np
from collections import defaultdict

from ..multigraph_operators import (
    get_mod_multigraph_operators,
    get_operators,
)


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
# Multigraph → Jet decomposition dictionaries
# -----------------------------------------------------------------------
# The multigraph and modified-multigraph bases share 14 operators with
# identical jet expansions.  Each basis adds 6 specific operators.
# -----------------------------------------------------------------------

_COMMON_JET = {
    "1":                [("1", 1.0)],
    "f":                [("f", 1.0)],
    "f²":               [("f²", 1.0)],
    "f³":               [("f³", 1.0)],
    "Δf":               [("f_xx", 1.0), ("f_yy", 1.0)],
    "‖∇f‖²":           [("f_x²", 1.0), ("f_y²", 1.0)],
    "f·‖∇f‖²":         [("f·f_x²", 1.0), ("f·f_y²", 1.0)],
    "Δ²f":              [("f_xxxx", 1.0), ("f_xxyy", 2.0), ("f_yyyy", 1.0)],
    "Δf·Δf":            [("f_xx²", 1.0), ("f_yy·f_xx", 2.0), ("f_yy²", 1.0)],
    "‖∇²f‖²":          [("f_xx²", 1.0), ("f_xy²", 2.0), ("f_yy²", 1.0)],
    "f·Δf·Δf":          [("f·f_xx²", 1.0), ("f·f_yy·f_xx", 2.0),
                         ("f·f_yy²", 1.0)],
    "f·‖∇²f‖²":        [("f·f_xx²", 1.0), ("f·f_xy²", 2.0),
                         ("f·f_yy²", 1.0)],
    "∇²f(∇f,∇f)":      [("f_x²·f_xx", 1.0), ("f_y·f_x·f_xy", 2.0),
                         ("f_y²·f_yy", 1.0)],
    "‖∇f‖²·Δf":        [("f_x²·f_xx", 1.0), ("f_x²·f_yy", 1.0),
                         ("f_y²·f_xx", 1.0), ("f_y²·f_yy", 1.0)],
}

_MULTIGRAPH_SPECIFIC = {
    "f·Δf":             [("f·f_xx", 1.0), ("f·f_yy", 1.0)],
    "f²·Δf":            [("f²·f_xx", 1.0), ("f²·f_yy", 1.0)],
    "f·Δ²f":            [("f·f_xxxx", 1.0), ("f·f_xxyy", 2.0),
                         ("f·f_yyyy", 1.0)],
    "⟨∇Δf,∇f⟩":        [("f_x·f_xxx", 1.0), ("f_x·f_xyy", 1.0),
                         ("f_y·f_xxy", 1.0), ("f_y·f_yyy", 1.0)],
    "f²·Δ²f":           [("f²·f_xxxx", 1.0), ("f²·f_xxyy", 2.0),
                         ("f²·f_yyyy", 1.0)],
    "f·⟨∇Δf,∇f⟩":      [("f·f_x·f_xxx", 1.0), ("f·f_x·f_xyy", 1.0),
                         ("f·f_y·f_xxy", 1.0), ("f·f_y·f_yyy", 1.0)],
}

_MOD_MULTIGRAPH_SPECIFIC = {
    # Δf² = 2f·Δf + 2‖∇f‖²
    "Δf²":              [("f·f_xx", 2.0), ("f·f_yy", 2.0),
                         ("f_x²", 2.0), ("f_y²", 2.0)],

    # Δf³ = 3f²·Δf + 6f·‖∇f‖²
    "Δf³":              [("f²·f_xx", 3.0), ("f²·f_yy", 3.0),
                         ("f·f_x²", 6.0), ("f·f_y²", 6.0)],

    # Δ‖∇f‖² = 2⟨∇Δf,∇f⟩ + 2‖∇²f‖²
    "Δ‖∇f‖²":          [("f_x·f_xxx", 2.0), ("f_x·f_xyy", 2.0),
                         ("f_y·f_xxy", 2.0), ("f_y·f_yyy", 2.0),
                         ("f_xx²", 2.0), ("f_xy²", 4.0), ("f_yy²", 2.0)],

    # Δ(f·‖∇f‖²) = 2f·⟨∇Δf,∇f⟩ + 2f·‖∇²f‖² + ‖∇f‖²·Δf + 4∇²f(∇f,∇f)
    "Δ(f·‖∇f‖²)":      [("f·f_x·f_xxx", 2.0), ("f·f_x·f_xyy", 2.0),
                         ("f·f_y·f_xxy", 2.0), ("f·f_y·f_yyy", 2.0),
                         ("f·f_xx²", 2.0), ("f·f_xy²", 4.0), ("f·f_yy²", 2.0),
                         ("f_x²·f_xx", 5.0), ("f_x²·f_yy", 1.0),
                         ("f_y²·f_xx", 1.0), ("f_y²·f_yy", 5.0),
                         ("f_y·f_x·f_xy", 8.0)],

    # Δ²f² = 2f·Δ²f + 8⟨∇Δf,∇f⟩ + 4‖∇²f‖² + 2(Δf)²
    "Δ²f²":             [("f·f_xxxx", 2.0), ("f·f_xxyy", 4.0),
                         ("f·f_yyyy", 2.0),
                         ("f_x·f_xxx", 8.0), ("f_x·f_xyy", 8.0),
                         ("f_y·f_xxy", 8.0), ("f_y·f_yyy", 8.0),
                         ("f_xx²", 6.0), ("f_xy²", 8.0), ("f_yy²", 6.0),
                         ("f_yy·f_xx", 4.0)],

    # Δ²f³ = 3f²·Δ²f + 24f·⟨∇Δf,∇f⟩ + 12f·‖∇²f‖² + 6f·(Δf)²
    #        + 24∇²f(∇f,∇f) + 12‖∇f‖²·Δf
    "Δ²f³":             [("f²·f_xxxx", 3.0), ("f²·f_xxyy", 6.0),
                         ("f²·f_yyyy", 3.0),
                         ("f·f_x·f_xxx", 24.0), ("f·f_x·f_xyy", 24.0),
                         ("f·f_y·f_xxy", 24.0), ("f·f_y·f_yyy", 24.0),
                         ("f·f_xx²", 18.0), ("f·f_xy²", 24.0),
                         ("f·f_yy²", 18.0), ("f·f_yy·f_xx", 12.0),
                         ("f_x²·f_xx", 36.0), ("f_y·f_x·f_xy", 48.0),
                         ("f_y²·f_yy", 36.0),
                         ("f_x²·f_yy", 12.0), ("f_y²·f_xx", 12.0)],
}

_MULTIGRAPH_TO_JET_DECOMPOSITION = {**_COMMON_JET, **_MULTIGRAPH_SPECIFIC}
_MOD_MULTIGRAPH_TO_JET_DECOMPOSITION = {**_COMMON_JET, **_MOD_MULTIGRAPH_SPECIFIC}


# -----------------------------------------------------------------------
# WSINDy → Jet: Leibniz-rule expansion
# -----------------------------------------------------------------------

def _apply_deriv_x(terms):
    """Apply ∂_x to a dict of {tuple: (count, coeff)}."""
    result = defaultdict(float)
    for key, (count, coeff) in terms.items():
        key_list = list(key)
        for i, (a, b) in enumerate(key_list):
            new_key = list(key_list)
            new_key[i] = (a + 1, b)
            new_key.sort()
            result[tuple(new_key)] += coeff
    return {k: (len(k), v) for k, v in result.items()}


def _apply_deriv_y(terms):
    """Apply ∂_y."""
    result = defaultdict(float)
    for key, (count, coeff) in terms.items():
        key_list = list(key)
        for i, (a, b) in enumerate(key_list):
            new_key = list(key_list)
            new_key[i] = (a, b + 1)
            new_key.sort()
            result[tuple(new_key)] += coeff
    return {k: (len(k), v) for k, v in result.items()}


def _expand_fq_deriv(a, b, q):
    """Expand ∂_x^a ∂_y^b (f^q) to jet monomials.

    Returns dict {((a₁,b₁), ..., (a_q,b_q)): coefficient}
    where each (a_i,b_i) is a derivative multi-index for one copy of f.
    """
    key = tuple((0, 0) for _ in range(q))
    terms = {key: (q, 1.0)}

    for _ in range(a):
        terms = _apply_deriv_x(terms)
    for _ in range(b):
        terms = _apply_deriv_y(terms)

    return terms


def _terms_to_jet_dict(terms, jet_names):
    """Convert expansion terms to a dict of {jet_description: coefficient}."""
    result = defaultdict(float)
    for key, (_count, coeff) in terms.items():
        exponents = defaultdict(int)
        for a, b in key:
            if a == 0 and b == 0:
                exponents[(0, 0)] += 1
            else:
                exponents[(a, b)] += 1
        parts = []
        for (a, b), exp in sorted(exponents.items(),
                                  key=lambda kv: (kv[0][0] + kv[0][1], kv[0][0])):
            if exp == 0:
                continue
            vname = jet_names.get((a, b), "")
            if not vname:
                if a == 0 and b == 0:
                    vname = "f"
                else:
                    vname = "f_" + "x" * a + "y" * b
            if exp == 1:
                parts.append(vname)
            else:
                sup_map = {2: "²", 3: "³", 4: "⁴", 5: "⁵",
                           6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}
                parts.append(f"{vname}{sup_map.get(exp, f'^{exp}')}")
        desc = "·".join(parts) if parts else "1"
        result[desc] += coeff
    return dict(result)


def build_wsindy_to_jet_map(desc_ws, desc_jet):
    """Build a matrix M such that w_jet = M @ w_wsindy.

    Each WSINDy column ∂_x^a ∂_y^b (f^q) is analytically expanded to
    jet monomials using the general Leibniz rule.
    """
    jet_names = {}
    for total in range(5):  # max_order=4
        for a in range(total + 1):
            b = total - a
            if a == 0 and b == 0:
                jet_names[(a, b)] = "f"
            else:
                jet_names[(a, b)] = "f_" + "x" * a + "y" * b

    M = np.zeros((len(desc_jet), len(desc_ws)))
    for j, desc in enumerate(desc_ws):
        if desc == "1":
            M[desc_jet.index("1"), j] = 1.0
            continue
        sup_map_rev = {"²": 2, "³": 3, "⁴": 4, "⁵": 5,
                       "⁶": 6, "⁷": 7, "⁸": 8, "⁹": 9}
        if desc.startswith("(f"):
            q_str = ""
            idx = 2
            while idx < len(desc) and desc[idx] in sup_map_rev:
                q_str += desc[idx]
                idx += 1
            q = sup_map_rev.get(q_str, 2) if q_str else 2
            rest = desc[desc.index(")_{") + 3: -1] if desc.endswith("}") \
                else desc[desc.index(")_{") + 3:]
            a = rest.count("x")
            b = rest.count("y")
        elif desc.startswith("f_{"):
            q = 1
            rest = desc[3:-1] if desc.endswith("}") else desc[3:]
            a = rest.count("x")
            b = rest.count("y")
        elif desc in ("f", "f²", "f³"):
            M[desc_jet.index(desc), j] = 1.0
            continue
        else:
            continue

        terms = _expand_fq_deriv(a, b, q)
        jet_dict = _terms_to_jet_dict(terms, jet_names)
        for jet_desc, coeff in jet_dict.items():
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


def map_multigraph_to_jet(w_eq, desc_eq, desc_jet):
    """Map multigraph coefficients to jet coefficients using the
    analytical decomposition.

    Raises KeyError if an operator name is not in the decomposition."""
    w_jet = np.zeros(len(desc_jet))
    for j, (name, coeff) in enumerate(zip(desc_eq, w_eq)):
        if abs(coeff) < 1e-14:
            continue
        if name not in _MULTIGRAPH_TO_JET_DECOMPOSITION:
            raise KeyError(
                f"Multigraph operator {name!r} has no jet decomposition. "
                f"Available: {list(_MULTIGRAPH_TO_JET_DECOMPOSITION)}"
            )
        for jet_name, mult in _MULTIGRAPH_TO_JET_DECOMPOSITION[name]:
            w_jet[desc_jet.index(jet_name)] += coeff * mult
    return w_jet


def map_mod_multigraph_to_jet(w_mod, desc_mod, desc_jet):
    """Map modified multigraph coefficients to jet coefficients using the
    analytical decomposition.

    Raises KeyError if an operator name is not in the decomposition."""
    w_jet = np.zeros(len(desc_jet))
    for j, (name, coeff) in enumerate(zip(desc_mod, w_mod)):
        if abs(coeff) < 1e-14:
            continue
        if name not in _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION:
            raise KeyError(
                f"Mod. multigraph operator {name!r} has no jet decomposition. "
                f"Available: {list(_MOD_MULTIGRAPH_TO_JET_DECOMPOSITION)}"
            )
        for jet_name, mult in _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION[name]:
            w_jet[desc_jet.index(jet_name)] += coeff * mult
    return w_jet
