"""
Analytical mappings between PDE term of different libraries.

Maps coefficients from the library bases use din the project
to polynomial-jet terms, which serves as the universal
comparison space.
"""

import numpy as np
from collections import defaultdict

# ===========================================================================
# Multigraph decomposition dictionaries
# ===========================================================================
# Easiest way is to hard-code the decompositions
# ===========================================================================

_MULTIGRAPH_TO_JET_DECOMPOSITION = {
    "1":                [("1", 1.0)],
    "f":                [("f", 1.0)],
    "f²":               [("f²", 1.0)],
    "f³":               [("f³", 1.0)],
    "Δf":               [("f_xx", 1.0), ("f_yy", 1.0)],
    "f·Δf":             [("f·f_xx", 1.0), ("f·f_yy", 1.0)],
    "‖∇f‖²":           [("f_x²", 1.0), ("f_y²", 1.0)],
    "f²·Δf":            [("f²·f_xx", 1.0), ("f²·f_yy", 1.0)],
    "f·‖∇f‖²":         [("f·f_x²", 1.0), ("f·f_y²", 1.0)],
    "Δ²f":              [("f_xxxx", 1.0), ("f_xxyy", 2.0), ("f_yyyy", 1.0)],
    "f·Δ²f":            [("f·f_xxxx", 1.0), ("f·f_xxyy", 2.0),
                         ("f·f_yyyy", 1.0)],
    "Δf·Δf":            [("f_xx²", 1.0), ("f_yy·f_xx", 2.0), ("f_yy²", 1.0)],
    "‖∇²f‖²":          [("f_xx²", 1.0), ("f_xy²", 2.0), ("f_yy²", 1.0)],
    "⟨∇Δf,∇f⟩":        [("f_x·f_xxx", 1.0), ("f_x·f_xyy", 1.0),
                         ("f_y·f_xxy", 1.0), ("f_y·f_yyy", 1.0)],
    "f²·Δ²f":           [("f²·f_xxxx", 1.0), ("f²·f_xxyy", 2.0),
                         ("f²·f_yyyy", 1.0)],
    "f·Δf·Δf":          [("f·f_xx²", 1.0), ("f·f_yy·f_xx", 2.0),
                         ("f·f_yy²", 1.0)],
    "f·‖∇²f‖²":        [("f·f_xx²", 1.0), ("f·f_xy²", 2.0),
                         ("f·f_yy²", 1.0)],
    "f·⟨∇Δf,∇f⟩":      [("f·f_x·f_xxx", 1.0), ("f·f_x·f_xyy", 1.0),
                         ("f·f_y·f_xxy", 1.0), ("f·f_y·f_yyy", 1.0)],
    "‖∇f‖²·Δf":        [("f_x²·f_xx", 1.0), ("f_x²·f_yy", 1.0),
                         ("f_y²·f_xx", 1.0), ("f_y²·f_yy", 1.0)],
    "∇²f(∇f,∇f)":      [("f_x²·f_xx", 1.0), ("f_y·f_x·f_xy", 2.0),
                         ("f_y²·f_yy", 1.0)],
}

_MOD_MULTIGRAPH_TO_JET_DECOMPOSITION = {
    # --- Unchanged operators (identical to _MULTIGRAPH_TO_JET_DECOMPOSITION) ---
    "1":                [("1", 1.0)],
    "f":                [("f", 1.0)],
    "f²":               [("f²", 1.0)],
    "f³":               [("f³", 1.0)],
    "Δf":               [("f_xx", 1.0), ("f_yy", 1.0)],
    "‖∇f‖²":           [("f_x²", 1.0), ("f_y²", 1.0)],
    "f·‖∇f‖²":         [("f·f_x²", 1.0), ("f·f_y²", 1.0)],
    "Δ²f":              [("f_xxxx", 1.0), ("f_xxyy", 2.0), ("f_yyyy", 1.0)],
    "‖∇²f‖²":          [("f_xx²", 1.0), ("f_xy²", 2.0), ("f_yy²", 1.0)],
    "Δf·Δf":            [("f_xx²", 1.0), ("f_yy·f_xx", 2.0), ("f_yy²", 1.0)],
    "f·‖∇²f‖²":        [("f·f_xx²", 1.0), ("f·f_xy²", 2.0),
                         ("f·f_yy²", 1.0)],
    "f·Δf·Δf":          [("f·f_xx²", 1.0), ("f·f_yy·f_xx", 2.0),
                         ("f·f_yy²", 1.0)],
    "∇²f(∇f,∇f)":      [("f_x²·f_xx", 1.0), ("f_y·f_x·f_xy", 2.0),
                         ("f_y²·f_yy", 1.0)],
    "‖∇f‖²·Δf":        [("f_x²·f_xx", 1.0), ("f_x²·f_yy", 1.0),
                         ("f_y²·f_xx", 1.0), ("f_y²·f_yy", 1.0)],

    # --- Modified operators (6 replacements) ---
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

# ===========================================================================
# WSINDy decomposition functions
# ===========================================================================
# WSINDy comes in derivatives of powers of f, e.g. ∂_x^a ∂_y^b (f^q).
# Each such term can be expanded analytically to a linear combination of
# polynomial-jet monomials using the general Leibniz rule.
# ===========================================================================


def _apply_deriv_x(terms):
    """Apply ∂_x to a dict of {frozenset((a,b),...): (count, coeff)}.
    Each key is a frozenset of (a,b) tuples representing copies of f
    with given derivative orders; count is how many copies; coeff is
    the scalar multiplier."""
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
    The jet monomial is ∏_i ∂_x^{a_i} ∂_y^{b_i} f.
    """
    # Start with q copies of undifferentiated f: ((0,0), (0,0), ..., (0,0))
    key = tuple((0, 0) for _ in range(q))
    terms = {key: (q, 1.0)}  # key → (count, coeff)

    for _ in range(a):
        terms = _apply_deriv_x(terms)
    for _ in range(b):
        terms = _apply_deriv_y(terms)

    return terms


def _terms_to_jet_dict(terms, jet_names):
    """Convert expansion terms to a dict of {jet_description: coefficient}.

    Each term key is a tuple of (a,b) pairs.  The jet monomial
    description is determined by the multiset of derivative orders.
    """
    result = defaultdict(float)
    for key, (_count, coeff) in terms.items():
        # Build jet variable exponents
        exponents = defaultdict(int)
        for a, b in key:
            if a == 0 and b == 0:
                exponents[(0, 0)] += 1  # bare f
            else:
                exponents[(a, b)] += 1
        # Build canonical description (same ordering as build_polynomial_jet_library)
        parts = []
        for (a, b), exp in sorted(exponents.items(), key=lambda kv: (kv[0][0] + kv[0][1], kv[0][0])):
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
                sup_map = {2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}
                parts.append(f"{vname}{sup_map.get(exp, f'^{exp}')}")
        desc = "·".join(parts) if parts else "1"
        result[desc] += coeff
    return dict(result)


def build_wsindy_to_jet_map(desc_ws, desc_jet):
    """Build a matrix M such that w_jet = M @ w_wsindy.

    Each WSINDy column ∂_x^a ∂_y^b (f^q) is analytically expanded to
    jet monomials using the general Leibniz rule.
    """
    # Build jet variable names
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
        # Parse WSINDy description: "f_{xx}", "(f²)_{xx}", "f²", etc.
        sup_map_rev = {"²": 2, "³": 3, "⁴": 4, "⁵": 5, "⁶": 6, "⁷": 7, "⁸": 8, "⁹": 9}
        if desc.startswith("(f"):
            # (f^q)_{x...y...}
            q_str = ""
            idx = 2
            while idx < len(desc) and desc[idx] in sup_map_rev:
                q_str += desc[idx]
                idx += 1
            q = sup_map_rev.get(q_str, 2) if q_str else 2
            # Parse derivative: )_{xxy...}
            rest = desc[desc.index(")_{") + 3: -1] if desc.endswith("}") else desc[desc.index(")_{") + 3:]
            a = rest.count("x")
            b = rest.count("y")
        elif desc.startswith("f_{"):
            # f_{xxy...}  (q=1)
            q = 1
            rest = desc[3:-1] if desc.endswith("}") else desc[3:]
            a = rest.count("x")
            b = rest.count("y")
        elif desc in ("f", "f²", "f³"):
            # Pure monomial: look up directly
            M[desc_jet.index(desc), j] = 1.0
            continue
        else:
            continue  # shouldn't happen

        terms = _expand_fq_deriv(a, b, q)
        jet_dict = _terms_to_jet_dict(terms, jet_names)
        for jet_desc, coeff in jet_dict.items():
            try:
                i = desc_jet.index(jet_desc)
            except ValueError:
                continue  # jet monomial not in library (e.g., order > 4)
            M[i, j] += coeff

    return M

# ===========================================================================
# Mapper functions: wrappers to map coefficients to polynomial-jet space
# ===========================================================================


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
    analytical decomposition in _MULTIGRAPH_TO_JET_DECOMPOSITION.

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
    analytical decomposition in _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION.

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
