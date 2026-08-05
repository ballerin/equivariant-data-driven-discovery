"""Equation ground-truth coefficients, utilities, and registry.

This module is the single source of truth for all benchmark PDE
definitions.  It is organised in three dependency-ordered sections:

  1. **Utilities** — general-purpose helpers and constants
     (``find_op``, ``CH_EPSILON``, ``mask_pme``, ``_pdefind_axis_labels``)
  2. **Ground-truth coefficient vectors** — one function per equation × basis
     that returns a length-``n_cols`` NumPy array of true operator coefficients
  3. **Equation registry & naming** — ``EQUATION_META`` maps canonical
     equation keys to their truth functions, and ``normalize_eq_name``
     canonicalises legacy / aliased names
"""

import numpy as np

# ============================================================================
# Section 1 — Utilities
# ============================================================================

# Cahn-Hilliard interface width — fixed at generation time and for ground truth.
CH_EPSILON = 0.2


def find_op(descriptions, name):
    """Return column index of operator *name* in the description list.

    Parameters
    ----------
    descriptions : list of str
        Operator descriptions from the library builder.
    name : str
        Exact description string (e.g. ``'Δf'``, ``'f³'``).

    Returns
    -------
    int

    Raises
    ------
    KeyError
        If *name* is not found.
    """
    try:
        return descriptions.index(name)
    except ValueError:
        raise KeyError(f"Operator {name!r} not in library. Available: {descriptions}")


def mask_pme(field):
    """Column mask for the Porous Medium Equation: discard near-zero values.

    Only spatial points where ``field > 0.05`` are retained.
    """
    return field.ravel() > 0.05


# PDE-FIND coordinate labels (up to 8-D).
_AXIS_LABELS = ['x', 'y', 'z', 'w', 'v', 'u', 't', 's']
_MAX_NDIM = len(_AXIS_LABELS)


def _pdefind_axis_labels(ndim):
    """Return list of axis label strings for the given dimension count."""
    if not 1 <= ndim <= _MAX_NDIM:
        raise ValueError(f"ndim must be 1-{_MAX_NDIM}, got {ndim}")
    return _AXIS_LABELS[:ndim]


# ============================================================================
# Section 2 — Ground-truth coefficient vectors
# ============================================================================
#
# Subsection 2a — Original multigraph basis
# -----------------------------------------

def true_allen_cahn(descriptions, eps=0.15):
    """Ground-truth coefficients for the Allen-Cahn equation.

    .. math::
        f_t = f - f^3 + \\nu\\,\\Delta f, \\qquad \\nu = \\varepsilon^2.
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, 'f')] = 1.0
    w[find_op(descriptions, 'f³')] = -1.0
    w[find_op(descriptions, 'Δf')] = eps ** 2
    return w


def true_swift_hohenberg(descriptions, r=0.3):
    """Ground-truth coefficients for the Swift-Hohenberg equation.

    .. math::
        f_t = -(1+\\Delta)^2 f + \\alpha f - f^3,
        \\qquad \\alpha = r.
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, 'f')] = r - 1.0
    w[find_op(descriptions, 'Δf')] = -2.0
    w[find_op(descriptions, 'Δ²f')] = -1.0
    w[find_op(descriptions, 'f³')] = -1.0
    return w


def true_ks(descriptions):
    """Ground-truth coefficients for the Kuramoto-Sivashinsky equation.

    .. math::
        f_t = -\\Delta f - \\Delta^2 f - 1/2\\|\\nabla f\\|^2.
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, '‖∇f‖²')] = -0.5
    w[find_op(descriptions, 'Δf')] = -1.0
    w[find_op(descriptions, 'Δ²f')] = -1.0
    return w


def true_cahn_hilliard(descriptions, eps=CH_EPSILON):
    """Ground-truth coefficients for the Cahn-Hilliard equation.

    .. math::
        f_t = -\\Delta^2 f + \\Delta(f^3 - f)
            = -\\Delta f + 3f^2\\Delta f + 6f\\|\\nabla f\\|^2
              - \\varepsilon^2 \\Delta^2 f.
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, 'Δf')] = -1.0
    w[find_op(descriptions, 'f²·Δf')] = 3.0
    w[find_op(descriptions, 'f·‖∇f‖²')] = 6.0
    w[find_op(descriptions, 'Δ²f')] = -(eps ** 2)
    return w


def true_pme(descriptions):
    """Ground-truth coefficients for the Porous Medium Equation (m=3).

    .. math::
        f_t = \\nabla\\cdot(f^2\\nabla f) = f^2\\Delta f + 2f\\|\\nabla f\\|^2
            = \\tfrac13\\Delta(f^3).

    In the standard convention :math:`u_t = \\Delta(u^m)` this is m = 3; the
    1/3 prefactor is inessential and absorbed by rescaling time.  The
    coefficients below (1.0, 2.0) belong to the ``div(f^2 grad f)``
    normalisation actually simulated — see ``simulations/porous_medium.py``.
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, 'f²·Δf')] = 1.0
    w[find_op(descriptions, 'f·‖∇f‖²')] = 2.0
    return w


# Subsection 2b — Modified multigraph (IBP-friendly) basis
# --------------------------------------------------------
# For equations whose terms are all in the unchanged set (AC, SH, KS),
# the ground truth is identical.  For CH and PME, the changed operator
# f²·Δf is replaced by Δf³, and coefficients are transformed via the
# change-of-basis matrix.

def true_allen_cahn_mod_multigraph(descriptions, eps=0.15):
    """Allen-Cahn in the modified multigraph basis.  Identical to original."""
    return true_allen_cahn(descriptions, eps=eps)


def true_swift_hohenberg_mod_multigraph(descriptions, r=0.3):
    """Swift-Hohenberg in the modified multigraph basis.  Identical to original."""
    return true_swift_hohenberg(descriptions, r=r)


def true_ks_mod_multigraph(descriptions):
    """Kuramoto-Sivashinsky in the modified multigraph basis.  Identical to original."""
    return true_ks(descriptions)


def true_cahn_hilliard_mod_multigraph(descriptions, eps=CH_EPSILON):
    """Cahn-Hilliard in the modified multigraph basis.

    In the original basis:  -Δf + 3·f²·Δf + 6·f·‖∇f‖² - ε²·Δ²f.
    In the mod. basis: Δf³ = 3·f²·Δf + 6·f·‖∇f‖², so the coefficient
    at Δf³ is 1.0 and f·‖∇f‖² is zero.
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, 'Δf')] = -1.0
    w[find_op(descriptions, 'Δf³')] = 1.0
    w[find_op(descriptions, 'Δ²f')] = -(eps ** 2)
    return w


def true_pme_mod_multigraph(descriptions):
    """Porous Medium Equation (m=3) in the modified multigraph basis.

    Original: f²·Δf + 2·f·‖∇f‖².
    Δf³ = 3·f²·Δf + 6·f·‖∇f‖², so:
      w[Δf³] = 1/3,  w[f·‖∇f‖²] = 0.
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, 'Δf³')] = 1.0 / 3.0
    return w


# Subsection 2c — PDE-FIND (coordinate-basis) ground truth
# --------------------------------------------------------

def true_allen_cahn_pdefind(ndim, eps=0.15):
    """Ground-truth coefficients for Allen-Cahn in the PDE-FIND basis.

    f_t = eps^2 * Delta f + f - f^3

    Delta f splits into *ndim* pure second derivatives (f_xx, f_yy, ...).
    """
    if not 1 <= ndim <= _MAX_NDIM:
        raise ValueError(f"ndim must be 1-{_MAX_NDIM}, got {ndim}")
    axes = _pdefind_axis_labels(ndim)
    result = {"f": 1.0, "f\u00b3": -1.0}
    for a in axes:
        result[f"f_{a}{a}"] = eps ** 2
    return result


def true_cahn_hilliard_pdefind(ndim, eps=CH_EPSILON):
    """Ground-truth coefficients for Cahn-Hilliard in the PDE-FIND basis.

    f_t = -eps^2 * Delta^2 f + Delta(f^3) - Delta f

    Expanding Delta(f^3) = 3 f^2 Delta f + 6 f ||grad f||^2:

    PDE-FIND expansion (partial -- f*(first_derivative)^2 terms omitted,
    since 6*f*||grad f||^2 has no f^p * d^alpha f column):
      - Pure 4th derivatives:        -eps^2     (ndim terms)
      - Mixed 4th derivatives:       -2*eps^2   (C(ndim,2) terms)
      - Pure 2nd derivatives:        -1         (ndim terms)
      - f^2 * pure 2nd derivatives:  +3         (ndim terms)

    The last two groups carry the sign of ``-Delta f`` and ``+3 f^2 Delta f``
    respectively, matching :func:`true_cahn_hilliard` in the multigraph basis.
    """
    if not 1 <= ndim <= _MAX_NDIM:
        raise ValueError(f"ndim must be 1-{_MAX_NDIM}, got {ndim}")
    axes = _pdefind_axis_labels(ndim)
    result = {}
    # Pure 4th: f_xxxx, f_yyyy, ...
    for a in axes:
        result[f"f_{a*4}"] = -(eps ** 2)
    # Mixed 4th: f_xxyy, f_xxzz, ... (lexicographic order)
    for i in range(ndim):
        for j in range(i + 1, ndim):
            a, b = axes[i], axes[j]
            result[f"f_{a*2}{b*2}"] = -2.0 * (eps ** 2)
    # Pure 2nd: f_xx, f_yy, ...  (from -Delta f)
    for a in axes:
        result[f"f_{a*2}"] = -1.0
    # f^2 * pure 2nd: f^2.f_xx, f^2.f_yy, ...  (from +3 f^2 Delta f)
    for a in axes:
        result[f"f\u00b2\u00b7f_{a*2}"] = 3.0
    return result


def true_swift_hohenberg_pdefind(ndim, r=0.3):
    """Ground-truth coefficients for Swift-Hohenberg in the PDE-FIND basis.

    f_t = (r-1) f - 2 Delta f - Delta^2 f - f^3

      - f                              →  r - 1
      - f^3                            → -1
      - Pure 2nd derivatives (Delta f) → -2        (ndim terms)
      - Pure 4th derivatives           → -1        (ndim terms)
      - Mixed 4th derivatives          → -2        (C(ndim,2) terms)

    Total active: 2 + 2*ndim + C(ndim,2).
    """
    if not 1 <= ndim <= _MAX_NDIM:
        raise ValueError(f"ndim must be 1-{_MAX_NDIM}, got {ndim}")
    axes = _pdefind_axis_labels(ndim)
    result = {"f": r - 1.0, "f\u00b3": -1.0}
    # Pure 2nd: f_xx, f_yy, ...  (from -2 Delta f)
    for a in axes:
        result[f"f_{a*2}"] = -2.0
    # Pure 4th: f_xxxx, f_yyyy, ...  (from -Delta^2 f)
    for a in axes:
        result[f"f_{a*4}"] = -1.0
    # Mixed 4th: f_xxyy, f_xxzz, ...
    for i in range(ndim):
        for j in range(i + 1, ndim):
            a, b = axes[i], axes[j]
            result[f"f_{a*2}{b*2}"] = -2.0
    return result


def true_pme_pdefind(ndim):
    """Ground-truth coefficients for Porous Medium in the PDE-FIND basis.

    PME:  f_t = Delta(f^3) / 3 = f^2 Delta f + 2 f ||grad f||^2

    The term f ||grad f||^2 = f * (f_x^2 + f_y^2 + ...) has NO representation
    in the PDE-FIND monomial basis f^p * ∂^α f (products of first derivatives
    are not PDE-FIND columns).  Therefore only the f^2 Delta f portion is
    encoded:

      - f^2 * pure 2nd derivatives → +1   (ndim terms, from f^2 Delta f)

    This is an *incomplete* representation — the truth vector is provided for
    benchmark purposes, but exact sparse recovery in the PDE-FIND basis is
    not expected.  The equivariant library is required for full recovery.
    """
    if not 1 <= ndim <= _MAX_NDIM:
        raise ValueError(f"ndim must be 1-{_MAX_NDIM}, got {ndim}")
    axes = _pdefind_axis_labels(ndim)
    result = {}
    for a in axes:
        result[f"f\u00b2\u00b7f_{a*2}"] = 1.0
    return result


# Subsection 2d — WSINDy basis ground truth
# -----------------------------------------

def true_pme_wsindy(descriptions):
    """Ground-truth coefficients for PME in the WSINDy basis ∂^α(f^q).

    In WSINDy, PME is exactly representable as:

        f_t = (1/3) ∂_x²(f³) + (1/3) ∂_y²(f³)

    because ∂_x²(f³) = 6 f·f_x² + 3 f²·f_{xx} and similarly for ∂_y²,
    so the sum recovers f²·Δf + 2 f·‖∇f‖².
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, '(f³)_{xx}')] = 1.0 / 3.0
    w[find_op(descriptions, '(f³)_{yy}')] = 1.0 / 3.0
    return w


def true_allen_cahn_wsindy(ndim, eps=0.15):
    """Ground-truth coefficients for Allen-Cahn in the WSINDy basis ∂^α(f^q).

    Allen-Cahn in N dimensions:

        f_t = ε² Δf + f − f³

    In WSINDy, Δf = Σ_i ∂_i² f, so the active columns are:

        f          →  1.0
        f³         → −1.0
        f_{aa}     → ε²   for each axis a ∈ {x, y, z, w, v, …}

    Returns a dict {column_name: coefficient} matching the pdefind truth
    convention so the sweep code can use the same ``truth_fn_ndim=True`` path.

    Total active terms: 2 + ndim.
    """
    d = {"f": 1.0, "f³": -1.0}
    for a in _AXIS_LABELS[:ndim]:
        d[f"f_{{{a}{a}}}"] = eps ** 2
    return d


def true_swift_hohenberg_wsindy(ndim, r=0.3):
    """Ground-truth coefficients for Swift-Hohenberg in the WSINDy basis.

    f_t = (r-1) f - 2 Δf - Δ²f - f³

    In WSINDy:
      - f                                → r - 1
      - f³                               → -1
      - f_{aa}                           → -2      (ndim pure 2nd derivs)
      - f_{aaaa}                         → -1      (ndim pure 4th derivs)
      - f_{aabb} (a < b)                 → -2      (C(ndim,2) mixed 4th)

    Total active: 2 + 2*ndim + C(ndim,2).
    """
    axes = _AXIS_LABELS[:ndim]
    d = {"f": r - 1.0, "f³": -1.0}
    # Pure 2nd: f_{xx}, f_{yy}, ...  (from -2 Δf)
    for a in axes:
        d[f"f_{{{a}{a}}}"] = -2.0
    # Pure 4th: f_{xxxx}, f_{yyyy}, ...  (from -Δ²f)
    for a in axes:
        d[f"f_{{{a*4}}}"] = -1.0
    # Mixed 4th: f_{xxyy}, f_{xxzz}, ...
    for i in range(ndim):
        for j in range(i + 1, ndim):
            a, b = axes[i], axes[j]
            d[f"f_{{{a*2}{b*2}}}"] = -2.0
    return d


def true_pme_wsindy_nd(ndim):
    """Ground-truth coefficients for PME in the WSINDy basis (N-D).

    PME:  f_t = (1/3) Δ(f³)

    In WSINDy, Δ(f³) = Σ_i ∂_i²(f³), so:

        (f³)_{aa}  →  1/3    for each axis a ∈ {x, y, z, w, v, …}

    Total active terms: ndim.
    """
    d = {}
    for a in _AXIS_LABELS[:ndim]:
        d[f"(f³)_{{{a}{a}}}"] = 1.0 / 3.0
    return d


def true_cahn_hilliard_wsindy(ndim, eps=CH_EPSILON):
    """Ground-truth coefficients for Cahn-Hilliard in the WSINDy basis.

    CH:  f_t = -ε² Δ²f + Δf - Δ(f³)

    In WSINDy (∂^α(f^q) monomials):
      - f_{aa}           → +1          (ndim terms, from +Δf)
      - f_{aaaa}         → -ε²         (ndim terms, from -ε²Δ²f)
      - f_{aabb} (a<b)   → -2ε²        (C(ndim,2) terms, from -ε²Δ²f)
      - (f³)_{aa}        → -1          (ndim terms, from -Δ(f³))

    Total active: 3*ndim + C(ndim,2).
    """
    axes = _AXIS_LABELS[:ndim]
    d = {}
    # Pure 2nd of f: f_{xx}, f_{yy}, ...  (from +Δf)
    for a in axes:
        d[f"f_{{{a}{a}}}"] = 1.0
    # Pure 4th of f: f_{xxxx}, f_{yyyy}, ...  (from -ε²Δ²f)
    for a in axes:
        d[f"f_{{{a*4}}}"] = -(eps ** 2)
    # Mixed 4th of f: f_{xxyy}, f_{xxzz}, ...
    for i in range(ndim):
        for j in range(i + 1, ndim):
            a, b = axes[i], axes[j]
            d[f"f_{{{a*2}{b*2}}}"] = -2.0 * (eps ** 2)
    # Pure 2nd of f³: (f³)_{xx}, (f³)_{yy}, ...  (from -Δ(f³))
    for a in axes:
        d[f"(f³)_{{{a}{a}}}"] = -1.0
    return d


# ============================================================================
# Section 3 — Equation registry & naming
# ============================================================================

# Maps canonical short name → (true_fn, skip_fraction, mask_fn).
# KS variants (ks_cell, ks_wchaos, ks_chaos) all use true_ks.
#
# Legacy data files use long names (allen_cahn, swift_hohenberg, …);
# ``normalize_eq_name`` maps them to the canonical short keys below.
EQUATION_META = {
    "ac":        (true_allen_cahn,       0.0,  None),
    "sh":        (true_swift_hohenberg,  0.1,  None),
    "ch":        (true_cahn_hilliard,     0.0,  None),
    "pme":       (true_pme,               0.0,  mask_pme),
    "ks_cell":   (true_ks,                0.3,  None),
    "ks_wchaos": (true_ks,                0.3,  None),
    "ks_chaos":  (true_ks,                0.3,  None),
}

# Mapping from any variant (old long name, new short name) → canonical key.
_EQ_NAME_ALIASES = {
    "allen_cahn": "ac",
    "swift_hohenberg": "sh",
    "cahn_hilliard": "ch",
    "porous_medium": "pme",
    "ks_cellular": "ks_cell",
    "ks_weak_chaotic": "ks_wchaos",
    "ks_chaotic": "ks_chaos",
    # canonical → canonical (identity)
    "ac": "ac", "sh": "sh", "ch": "ch", "pme": "pme",
    "ks_cell": "ks_cell", "ks_wchaos": "ks_wchaos", "ks_chaos": "ks_chaos",
}


def normalize_eq_name(name):
    """Map any equation name variant to the canonical short key.

    Accepts both legacy long names (``"allen_cahn"``, ``"swift_hohenberg"``,
    …) and the new short names (``"ac"``, ``"sh"``, …).  Returns the
    canonical short key, or the original *name* unchanged if no mapping
    is known (useful for KS regime sub-names that are already canonical).
    """
    return _EQ_NAME_ALIASES.get(name, name)
