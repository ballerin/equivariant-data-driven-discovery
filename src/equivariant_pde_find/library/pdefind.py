"""
Pointwise PDE-FIND (coordinate-based) library construction.

Internal module — use ``library.build_pdefind_library()`` (the public
facade) instead of calling ``_build_pdefind_pointwise`` directly.

N-D generalization: both the pointwise and weak-form builders work for any
ndim, spectral or polynomial backend (isotropic domain beyond ndim=3; see
derivatives/spectral.py's module docstring).  Axis names beyond z use
jet_algebra.AXIS_LABELS (w, v, u, t, s), so this supports up to ndim=8.
"""

import numpy as np

from .. import derivatives as _derivatives
from ..derivatives.polynomial import poly_partial_1d
from ..jet_algebra import AXIS_LABELS
from ._nd_common import enumerate_multi_indices, compute_derivs_nd


# Unicode superscript digits for polynomial exponents in descriptions.
_SUP = {2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}


def _sup(n: int) -> str:
    """Return a Unicode superscript string for integer n >= 2."""
    return _SUP.get(n, f"^{n}")


# ---------------------------------------------------------------------------
# N-D multi-index naming (enumeration/derivative computation is shared via
# _nd_common.py; matches the historical 2-D/3-D lexicographic ordering)
# ---------------------------------------------------------------------------


def _multi_index_name(multi_index):
    """'f_xxy' style name for a derivative multi-index, '1' for the zero
    multi-index."""
    if not any(multi_index):
        return "1"
    if len(multi_index) > len(AXIS_LABELS):
        raise ValueError(f"ndim={len(multi_index)} exceeds available axis labels")
    parts = "".join(AXIS_LABELS[axis] * a for axis, a in enumerate(multi_index) if a)
    return "f_" + parts


def _build_pdefind_pointwise(f, Lx, Ly, max_order=4, max_poly=3, Lz=None):
    """Build a coordinate-based library as used in PDE-Find, for any ndim.

    Columns = (all partial derivatives of f up to order max_order) ×
              (monomials f^0, f^1, ..., f^max_poly).
    For max_order=4 and max_poly=3 this gives 15 × 4 = 60 columns (2-D).

    The differentiation method is determined by the active backend in
    ``derivatives.py`` (set via ``derivatives.set_backend()`` or
    ``derivatives.use_backend()``):

    - ``'spectral'`` (default): FFT-based exact differentiation.  Accurate for
      noiseless periodic data; matches the original PDE-FIND Fourier option.
    - ``'polynomial'``: local ndim polynomial least-squares fitting.
      Noise-robust; uses a joint patch fit instead of separable 1-D fits,
      which avoids error compounding for mixed partials.

    Parameters
    ----------
    f         : ndarray, shape (N_1, ..., N_ndim) - data
    Lx, Ly, Lz : float - domain lengths (isotropic beyond ndim=3)
    max_order : int - highest partial derivative order (default 4)
    max_poly  : int - highest polynomial degree factor (default 3)

    Returns
    -------
    Theta        : ndarray, shape (f.size, n_cols)
    descriptions : list of strings
    """
    method = _derivatives.get_backend()
    ndim = f.ndim
    Ls = _derivatives.axis_lengths(f, Lx, Ly, Lz)
    multi_indices = enumerate_multi_indices(ndim, max_order)
    deriv_map = compute_derivs_nd(f, Ls, multi_indices, method, zero_value="ones")

    cols, descs = [], []
    for pw in range(max_poly + 1):
        poly = f**pw if pw > 0 else np.ones_like(f)
        for mi in multi_indices:
            dname = _multi_index_name(mi)
            col = (poly * deriv_map[mi]).ravel()
            if pw == 0:
                desc = dname
            elif pw == 1:
                desc = f"f·{dname}" if dname != "1" else "f"
            else:
                desc = f"f{_sup(pw)}·{dname}" if dname != "1" else f"f{_sup(pw)}"
            cols.append(col)
            descs.append(desc)

    return np.column_stack(cols), descs


def _build_pdefind_pointwise_1d(u, L, max_deriv=3, max_poly=3):
    """Build a coordinate-based library for a 1D periodic scalar field.

    Columns = (derivatives ∂^d u for d = 0, ..., max_deriv) ×
              (monomials u^p for p = 0, ..., max_poly).

    The differentiation method is determined by the active backend in
    ``derivatives.py`` (set via ``derivatives.set_backend()`` or
    ``derivatives.use_backend()``):

    - ``'spectral'`` (default): FFT-based exact differentiation.  Matches
      PDE-FIND's ``space_diff='Fourier'``.
    - ``'polynomial'``: local 1-D polynomial least-squares fitting via
      ``poly_partial_1d``.  Noise-robust, consistent with the ndim
      polynomial backend.

    For max_deriv=3 and max_poly=3 this gives 4 × 4 = 16 columns with the
    same ordering as PDE-FIND's build_linear_system(D=3, P=3).

    Parameters
    ----------
    u         : ndarray, shape (N,) — field values on a uniform periodic grid
    L         : float — domain length
    max_deriv : int — highest derivative order (default 3)
    max_poly  : int — highest polynomial degree factor (default 3)

    Returns
    -------
    Theta        : ndarray, shape (N, n_cols)
    descriptions : list of n_cols strings
    """
    method = _derivatives.get_backend()
    N = len(u)

    derivs = [np.ones(N)]
    deriv_names = ["1"]
    if method == "spectral":
        k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
        uhat = np.fft.fft(u)
        for d in range(1, max_deriv + 1):
            derivs.append(np.real(np.fft.ifft((1j * k) ** d * uhat)))
            deriv_names.append("u_" + "x" * d)
    else:
        for d in range(1, max_deriv + 1):
            derivs.append(poly_partial_1d(u, L, d))
            deriv_names.append("u_" + "x" * d)

    cols, descs = [], []
    for d_arr, d_name in zip(derivs, deriv_names):
        for p in range(max_poly + 1):
            col = u**p * d_arr
            if p == 0:
                desc = d_name
            elif p == 1:
                desc = f"u·{d_name}" if d_name != "1" else "u"
            else:
                desc = f"u{_sup(p)}·{d_name}" if d_name != "1" else f"u{_sup(p)}"
            cols.append(col)
            descs.append(desc)

    return np.column_stack(cols), descs


def _build_pdefind_space_weak(f, Lx, Ly, max_order=4, max_poly=3,
                              m=8, p=4, space_ibp=False, Lz=None):
    """Build the weak-form PDE-FIND library matrix, for any ndim.

    IBP is applied to pure linear terms (poly degree 0): the full
    multi-index derivative is moved onto the test function, one factor
    (plain φ or a derivative of φ) per axis.
    """
    from .weak_form import _cached_kernel_1d, _test_function_deriv_1d, weak_convolution

    method = _derivatives.get_backend()
    ndim = f.ndim
    Ls = _derivatives.axis_lengths(f, Lx, Ly, Lz)
    kernel = _cached_kernel_1d(m, p)

    multi_indices = enumerate_multi_indices(ndim, max_order)
    deriv_map = compute_derivs_nd(f, Ls, multi_indices, method, zero_value="ones")

    cols, descs = [], []
    for poly_deg in range(max_poly + 1):
        poly_field = f**poly_deg if poly_deg > 0 else np.ones_like(f)
        for mi in multi_indices:
            dname = _multi_index_name(mi)
            order = sum(mi)

            if space_ibp and poly_deg == 0 and order > 0:
                # IBP: move ∂^mi entirely onto the test function, one axis
                # at a time (identity axes get plain φ = kernel).
                kernel_list = []
                for axis, a in enumerate(mi):
                    if a == 0:
                        kernel_list.append(kernel)
                    else:
                        dx_a = (Ls[axis] / f.shape[axis]) ** a
                        kernel_list.append(_test_function_deriv_1d(m, p, a) / dx_a)
                col_weak = weak_convolution(f, kernel_list)
            else:
                col_raw = poly_field * deriv_map[mi]
                col_weak = weak_convolution(col_raw, [kernel] * ndim)

            cols.append(col_weak.ravel())

            if poly_deg == 0:
                desc = dname
            elif poly_deg == 1:
                desc = f"f·{dname}" if dname != "1" else "f"
            else:
                desc = f"f{_sup(poly_deg)}·{dname}" if dname != "1" else f"f{_sup(poly_deg)}"
            descs.append(desc)

    return np.column_stack(cols), descs
