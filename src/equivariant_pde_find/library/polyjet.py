"""
Pointwise polynomial-jet library construction.

Internal module — use ``library.build_polyjet_library()`` (the public
facade) instead of calling ``_build_polyjet_pointwise`` directly.

N-D generalization: both the pointwise and weak-form builders work for any
ndim, spectral or polynomial backend (isotropic domain beyond ndim=3; see
derivatives/spectral.py's module docstring).  Axis names beyond z use
jet_algebra.AXIS_LABELS (w, v, u, t, s), so this supports up to ndim=8.
"""

import numpy as np

from .. import derivatives as _derivatives
from ..jet_algebra import AXIS_LABELS, parse_jet_monomial
from ._nd_common import enumerate_multi_indices, compute_derivs_nd


# Unicode superscript digits for polynomial exponents in descriptions.
_SUP = {2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}


def _sup(n: int) -> str:
    """Return a Unicode superscript string for integer n >= 2."""
    return _SUP.get(n, f"^{n}")


def _jet_var_name(multi_index):
    """'f_xxy' style name for a jet variable, 'f' for the zero multi-index."""
    if not any(multi_index):
        return "f"
    if len(multi_index) > len(AXIS_LABELS):
        raise ValueError(f"ndim={len(multi_index)} exceeds available axis labels")
    parts = "".join(AXIS_LABELS[axis] * a for axis, a in enumerate(multi_index) if a)
    return "f_" + parts


def _build_jet_vars(ndim, max_order):
    """[(multi_index, order, name), ...] for all jet variables up to max_order."""
    return [
        (mi, sum(mi), _jet_var_name(mi))
        for mi in enumerate_multi_indices(ndim, max_order)
    ]


def _describe_exponents(exponents, jet_names):
    """Render a monomial (dict of multi_index -> exponent) as e.g. 'f_x²·f_xx'."""
    parts = []
    for mi, exp in sorted(exponents.items(), key=lambda kv: (sum(kv[0]), kv[0])):
        if exp == 0:
            continue
        name = jet_names[mi]
        parts.append(name if exp == 1 else f"{name}{_sup(exp)}")
    return "1" if not parts else "·".join(parts)


def _build_polyjet_pointwise(f, Lx, Ly, max_order=4, max_degree=3, Lz=None):
    """Build the polynomial jet library: all monomials in jet variables, for
    any ndim.

    The *jet variables* are all partial derivatives ∂^α f for
    0 <= |α| <= max_order (including f itself).  A library column is a
    product (monomial) of these jet variables:

        ∏_α (∂^α f)^{e_α}

    subject to:
      - total degree:  Σ e_α <= max_degree
      - total differential order:  Σ e_α · |α| <= max_order

    For max_order=4 and max_degree=3 this gives 105 columns (2-D).

    Parameters
    ----------
    f          : ndarray, shape (N_1, ..., N_ndim) - data
    Lx, Ly, Lz : float - domain lengths (isotropic beyond ndim=3)
    max_order  : int - maximum total differential order
    max_degree : int - maximum total polynomial degree

    Returns
    -------
    Theta        : ndarray, shape (f.size, n_cols)
    descriptions : list of n_cols strings
    """
    method = _derivatives.get_backend()
    ndim = f.ndim
    Ls = _derivatives.axis_lengths(f, Lx, Ly, Lz)

    jet_vars = _build_jet_vars(ndim, max_order)
    multi_indices = [mi for mi, _order, _name in jet_vars]
    jet_names = {mi: name for mi, _order, name in jet_vars}
    jet_fields = compute_derivs_nd(f, Ls, multi_indices, method, zero_value="f")

    monomials = enumerate_jet_monomials(jet_vars, max_order, max_degree)

    cols, descs = [], []
    for exponents, deg, order in monomials:
        col = np.ones_like(f)
        for mi, exp in exponents.items():
            if exp > 0:
                col = col * (jet_fields[mi] ** exp)
        cols.append(col.ravel())
        descs.append(_describe_exponents(exponents, jet_names))

    return np.column_stack(cols), descs


def _build_polyjet_from_descriptions(f, Lx, Ly, descriptions, Lz=None):
    """Build a polynomial-jet library from an explicit list of monomial
    description strings (rather than enumerating all monomials under caps).

    This is the dynamic counterpart of ``_build_polyjet_pointwise``: instead
    of a recursive enumeration bounded by max_order/max_degree, it builds
    exactly the columns named in *descriptions*.  The union of all needed
    jet variables is computed first so that ``compute_derivs_nd`` is called
    only once.

    Parameters
    ----------
    f : ndarray, shape (N_1, ..., N_ndim)
    Lx, Ly, Lz : float — domain lengths
    descriptions : list of str
        Monomial description strings in the canonical format produced by
        ``_describe_exponents`` / ``describe_monomial``.

    Returns
    -------
    Theta : ndarray, shape (f.size, len(descriptions))
    descriptions : list of str (the same list, returned for API consistency)
    """
    method = _derivatives.get_backend()
    ndim = f.ndim
    Ls = _derivatives.axis_lengths(f, Lx, Ly, Lz)

    # Parse all descriptions and collect the union of needed jet variables
    all_exponents = []
    needed_mi = set()
    for desc in descriptions:
        exponents = parse_jet_monomial(desc, ndim)
        all_exponents.append(exponents)
        needed_mi.update(exponents.keys())

    # Include (0,...,0) for "f" itself, always
    zero_mi = (0,) * ndim
    needed_mi.add(zero_mi)

    # Build jet variable names for description formatting
    jet_names = {}
    jet_names[zero_mi] = "f"
    for mi in needed_mi:
        if mi == zero_mi:
            continue
        parts = "".join(AXIS_LABELS[i] * a for i, a in enumerate(mi) if a)
        jet_names[mi] = "f_" + parts

    multi_indices = list(needed_mi)
    jet_fields = compute_derivs_nd(f, Ls, multi_indices, method, zero_value="f")
    # f itself is available as jet_fields[zero_mi]

    cols = []
    for exponents in all_exponents:
        col = np.ones_like(f)
        for mi, exp in exponents.items():
            if exp > 0:
                col = col * (jet_fields[mi] ** exp)
        cols.append(col.ravel())

    return np.column_stack(cols), list(descriptions)


def enumerate_jet_monomials(jet_vars, max_order, max_degree):
    """Enumerate all monomials in jet variables satisfying constraints.

    Parameters
    ----------
    jet_vars : list of (multi_index, order, name)
        multi_index is any hashable value (historically a 2-D (a, b) tuple;
        now an ndim-tuple) -- used only as an opaque dict key here.
    max_order : int
    max_degree : int

    Returns
    -------
    list of (exponents_dict, total_degree, total_order)
        Sorted by (total_order, total_degree).
    """
    n_vars = len(jet_vars)
    results = []

    def _recurse(idx, exponents, cur_deg, cur_order):
        if idx == n_vars:
            results.append((dict(exponents), cur_deg, cur_order))
            return

        multi_index, order, _name = jet_vars[idx]
        max_by_deg = max_degree - cur_deg
        max_by_order = (
            (max_order - cur_order) // order if order > 0 else max_by_deg
        )
        max_exp = min(max_by_deg, max_by_order)

        for e in range(max_exp + 1):
            exponents[multi_index] = e
            _recurse(idx + 1, exponents, cur_deg + e, cur_order + e * order)

    _recurse(0, {}, 0, 0)

    results.sort(key=lambda x: (x[2], x[1]))
    return results


def _build_polyjet_weak(f, Lx, Ly, max_order=4, max_degree=3,
                        m=8, p=4, use_ibp=False, Lz=None):
    """Build the weak-form polynomial jet library, for any ndim.

    Computes each jet monomial pointwise and convolves with ψ.  When
    ``use_ibp=True``, IBP is applied to single-derivative terms of the
    form ∂^α f (i.e., degree-1 monomials with order > 0), moving the full
    multi-index derivative onto the test function.

    Parameters
    ----------
    f : ndarray, shape (N_1, ..., N_ndim)
    Lx, Ly, Lz : float — domain lengths (isotropic beyond ndim=3)
    max_order  : int — maximum total differential order (default 4)
    max_degree : int — maximum total polynomial degree (default 3)
    m, p : int — test function parameters
    use_ibp : bool — enable IBP for single-derivative columns

    Returns
    -------
    Theta_weak : ndarray of shape (N_valid, K)
    descriptions : list of K str
    """
    from .weak_form import _cached_kernel_1d, _test_function_deriv_1d, weak_convolution

    method = _derivatives.get_backend()
    ndim = f.ndim
    Ls = _derivatives.axis_lengths(f, Lx, Ly, Lz)
    kernel = _cached_kernel_1d(m, p)

    jet_vars = _build_jet_vars(ndim, max_order)
    multi_indices = [mi for mi, _order, _name in jet_vars]
    jet_names = {mi: name for mi, _order, name in jet_vars}
    jet_fields = compute_derivs_nd(f, Ls, multi_indices, method, zero_value="f")

    monomials = enumerate_jet_monomials(jet_vars, max_order, max_degree)

    cols, descs = [], []
    for exponents, deg, order in monomials:
        non_one_terms = [(k, v) for k, v in exponents.items() if v > 0]
        is_single_deriv = (
            use_ibp
            and len(non_one_terms) == 1
            and sum(non_one_terms[0][0]) > 0
            and non_one_terms[0][1] == 1
        )

        if is_single_deriv:
            mi, _exp = non_one_terms[0]
            kernel_list = []
            for axis, a in enumerate(mi):
                if a == 0:
                    kernel_list.append(kernel)
                else:
                    dx_a = (Ls[axis] / f.shape[axis]) ** a
                    kernel_list.append(_test_function_deriv_1d(m, p, a) / dx_a)
            col_weak = weak_convolution(f, kernel_list)
        else:
            col_raw = np.ones_like(f)
            for mi, exp in exponents.items():
                if exp > 0:
                    col_raw = col_raw * (jet_fields[mi] ** exp)
            col_weak = weak_convolution(col_raw, [kernel] * ndim)

        cols.append(col_weak.ravel())
        descs.append(_describe_exponents(exponents, jet_names))

    return np.column_stack(cols), descs
