"""
Pointwise WSINDy library construction.

Internal module — use ``library.build_wsindy_library()`` (the public
facade) instead of calling ``_build_wsindy_pointwise`` directly.

N-D generalization: both the pointwise and weak-form builders work for any
ndim, spectral or polynomial backend (isotropic domain beyond ndim=3; see
derivatives/spectral.py's module docstring).  Axis names beyond z use
jet_algebra.AXIS_LABELS (w, v, u, t, s), so this supports up to ndim=8.
"""

import numpy as np

from .. import derivatives as _derivatives
from ..jet_algebra import AXIS_LABELS
from ._nd_common import enumerate_multi_indices, compute_derivs_nd


# Unicode superscript digits for polynomial exponents in descriptions.
_SUP = {2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}


def _sup(n: int) -> str:
    """Return a Unicode superscript string for integer n >= 2."""
    return _SUP.get(n, f"^{n}")


def _axis_part(multi_index):
    """'xxy' style axis-letter string for a nonzero derivative multi-index."""
    if len(multi_index) > len(AXIS_LABELS):
        raise ValueError(f"ndim={len(multi_index)} exceeds available axis labels")
    return "".join(AXIS_LABELS[axis] * a for axis, a in enumerate(multi_index) if a)


def _deriv_name(multi_index, q):
    """WSINDy column name for ∂^multi_index(f^q), matching the historical
    curly-brace subscript convention: 'f_{xx}' for q=1, '(f²)_{xx}' for q>=2."""
    dx_part = _axis_part(multi_index)
    if q == 1:
        return f"f_{{{dx_part}}}"
    return f"(f{_sup(q)})_{{{dx_part}}}"


def _build_wsindy_pointwise(f, Lx, Ly, max_order=4, max_poly=3, Lz=None):
    """Build the pointwise WSINDy library from data, for any ndim.

    This is the pointwise analogue of the library used in WSINDy, for
    benchmark purposes only.

    For max_order=4 and max_poly=3 this gives 46 columns (2-D).

    Parameters
    ----------
    f         : ndarray, shape (N_1, ..., N_ndim) - data
    Lx, Ly, Lz : float - domain lengths (isotropic beyond ndim=3)
    max_order : int - highest partial derivative order (default 4)
    max_poly  : int - highest monomial power q in ∂^α(f^q) (default 3)

    Returns
    -------
    Theta        : ndarray, shape (f.size, n_cols)
    descriptions : list of n_cols strings
    """
    method = _derivatives.get_backend()
    ndim = f.ndim
    Ls = _derivatives.axis_lengths(f, Lx, Ly, Lz)
    multi_indices = enumerate_multi_indices(ndim, max_order)
    nonzero_indices = [mi for mi in multi_indices if any(mi)]

    cols, descs = [], []

    # Constant term: partial^0(f^0) = 1
    cols.append(np.ones(f.size, dtype=np.float64))
    descs.append("1")

    for q in range(1, max_poly + 1):
        f_pow = f**q
        deriv_map = compute_derivs_nd(f_pow, Ls, nonzero_indices, method, zero_value="f")
        for mi in nonzero_indices:
            cols.append(deriv_map[mi].ravel())
            descs.append(_deriv_name(mi, q))

    # Pure monomials f^q for q >= 1
    for q in range(1, max_poly + 1):
        cols.append((f**q).ravel())
        descs.append("f" if q == 1 else f"f{_sup(q)}")

    return np.column_stack(cols), descs


def _build_wsindy_space_weak(f, Lx, Ly, max_order=4, max_poly=3,
                             m=8, p=4, space_ibp=False, Lz=None):
    """Build the weak-form WSINDy library: ψ ∗ ∂^α(f^q), for any ndim.

    Mirrors ``_build_wsindy_pointwise`` but convolves each column with the
    polynomial-bump test function ψ.  When ``space_ibp=True``, integration by
    parts moves the derivative onto the test function for *every* monomial
    power q, so f^q is never differentiated numerically:

        ∫ ψ · ∂^α(f^q)  =  (-1)^{|α|} ∫ ∂^α ψ · f^q

    This is the formulation of Messenger & Bortz (arXiv:2007.02848), whose
    Gram matrix (their eq. 3.3) is ⟨(-1)^{|α|} ∂^α ψ, f_j(U)⟩ with the
    nonlinearity f_j evaluated pointwise on the raw data.  Applying IBP only
    to q=1 would leave the q >= 2 columns differentiating the data, which
    the weak form exists to avoid.

    Parameters
    ----------
    f : ndarray, shape (N_1, ..., N_ndim)
    Lx, Ly, Lz : float — domain lengths (isotropic beyond ndim=3)
    max_order : int — highest derivative order (default 4)
    max_poly  : int — highest monomial power q (default 3)
    m, p : int — test function parameters
    space_ibp : bool — enable IBP for all ∂^α(f^q) columns

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
    multi_indices = enumerate_multi_indices(ndim, max_order)
    nonzero_indices = [mi for mi in multi_indices if any(mi)]

    cols, descs = [], []

    # Constant term: psi * 1
    col_weak = weak_convolution(np.ones_like(f), [kernel] * ndim)
    cols.append(col_weak.ravel())
    descs.append("1")

    for q in range(1, max_poly + 1):
        f_pow = f**q

        if space_ibp:
            deriv_map = {}
        else:
            deriv_map = compute_derivs_nd(f_pow, Ls, nonzero_indices, method, zero_value="f")

        for mi in nonzero_indices:
            if space_ibp:
                # IBP: move partial^mi entirely onto the test function, so
                # the derivative acts on psi and f^q is left untouched.
                kernel_list = []
                for axis, a in enumerate(mi):
                    if a == 0:
                        kernel_list.append(kernel)
                    else:
                        dx_a = (Ls[axis] / f.shape[axis]) ** a
                        kernel_list.append(_test_function_deriv_1d(m, p, a) / dx_a)
                col_weak = weak_convolution(f_pow, kernel_list)
            else:
                col_weak = weak_convolution(deriv_map[mi], [kernel] * ndim)

            cols.append(col_weak.ravel())
            descs.append(_deriv_name(mi, q))

    # Pure monomials f^q for q >= 1
    for q in range(1, max_poly + 1):
        col_weak = weak_convolution(f**q, [kernel] * ndim)
        cols.append(col_weak.ravel())
        descs.append("f" if q == 1 else f"f{_sup(q)}")

    return np.column_stack(cols), descs
