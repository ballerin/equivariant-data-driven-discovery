"""
Pointwise polynomial-jet library construction.

Internal module — use ``library.build_polyjet_library()`` (the public
facade) instead of calling ``_build_polyjet_pointwise`` directly.
"""

import numpy as np

from .. import derivatives as _derivatives
from ..derivatives.polynomial import poly_partial


# Unicode superscript digits for polynomial exponents in descriptions.
_SUP = {2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}


def _sup(n: int) -> str:
    """Return a Unicode superscript string for integer n >= 2."""
    return _SUP.get(n, f"^{n}")


def _wavenumbers_2d(Nx, Ny, Lx, Ly):
    """Return 2-D wavenumber grids (KX, KY) for an (Nx, Ny) FFT grid."""
    kx = 2.0 * np.pi / Lx * np.fft.fftfreq(Nx) * Nx
    ky = 2.0 * np.pi / Ly * np.fft.fftfreq(Ny) * Ny
    return np.meshgrid(kx, ky, indexing="ij")


def _build_polyjet_pointwise(f, Lx, Ly, max_order=4, max_degree=3):
    """Build the polynomial jet library: all monomials in jet variables.

    The *jet variables* are all partial derivatives ∂_x^a ∂_y^b f for
    0 ≤ a+b ≤ max_order (including f itself).  A library column is a
    product (monomial) of these jet variables:

        ∏_{a,b} (∂_x^a ∂_y^b f)^{e_{a,b}}

    subject to:
      - total degree:  Σ e_{a,b} ≤ max_degree
      - total differential order:  Σ e_{a,b} · (a+b) ≤ max_order

    For max_order=4 and max_degree=3 this gives 105 columns.

    Parameters
    ----------
    f          : ndarray, shape (Nx, Ny) - data
    Lx, Ly     : float - domain lengths
    max_order  : int - maximum total differential order
    max_degree : int - maximum total polynomial degree

    Returns
    -------
    Theta        : ndarray, shape (Nx*Ny, n_cols)
    descriptions : list of n_cols strings
    """
    method = _derivatives.get_backend()
    Nx, Ny = f.shape

    # --- 1. Enumerate jet variables and compute all needed derivatives ---
    jet_vars = []  # list of ((a, b), order=a+b, name)
    for total in range(max_order + 1):
        for a in range(total + 1):
            b = total - a
            name = (("x" * a) + ("y" * b)) if total > 0 else ""
            name = f"f_{name}" if total > 0 else "f"
            jet_vars.append(((a, b), total, name))

    # Compute each jet variable once
    if method == "spectral":
        KX, KY = _wavenumbers_2d(Nx, Ny, Lx, Ly)
        fh = np.fft.fft2(f)

    jet_fields = {}   # (a,b) → ndarray
    jet_names = {}    # (a,b) → str
    for (a, b), order, name in jet_vars:
        jet_names[(a, b)] = name
        if a + b == 0:
            jet_fields[(a, b)] = f.copy()
        elif method == "spectral":
            d_hat = (1j * KX) ** a * (1j * KY) ** b * fh
            jet_fields[(a, b)] = np.real(np.fft.ifft2(d_hat))
        else:
            jet_fields[(a, b)] = poly_partial(f, Lx, Ly, a, b)

    # --- 2. Enumerate all monomials ---
    monomials = enumerate_jet_monomials(jet_vars, max_order, max_degree)

    # --- 3. Build each column ---
    cols, descs = [], []
    for exponents, deg, order in monomials:
        # Build the field: start with ones, multiply by each jet variable
        # raised to its exponent
        col = np.ones_like(f)
        for (a, b), exp in exponents.items():
            if exp > 0:
                col = col * (jet_fields[(a, b)] ** exp)
        cols.append(col.ravel())

        # Build description string
        parts = []
        # Sort by derivative order for readability
        for (a, b), exp in sorted(exponents.items(), key=lambda kv: (kv[0][0] + kv[0][1], kv[0][0])):
            if exp == 0:
                continue
            var_name = jet_names[(a, b)]
            if exp == 1:
                parts.append(var_name)
            else:
                parts.append(f"{var_name}{_sup(exp)}")
        if not parts:
            descs.append("1")
        else:
            descs.append("·".join(parts))

    return np.column_stack(cols), descs


def enumerate_jet_monomials(jet_vars, max_order, max_degree):
    """Enumerate all monomials in jet variables satisfying constraints.

    Parameters
    ----------
    jet_vars : list of ((a, b), order, name)
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
            # Include the constant (all exponents zero) and all others
            results.append((dict(exponents), cur_deg, cur_order))
            return

        (a, b), order, _name = jet_vars[idx]
        # Maximum exponent for this variable given remaining budget
        max_by_deg = max_degree - cur_deg
        max_by_order = (
            (max_order - cur_order) // order if order > 0 else max_by_deg
        )
        max_exp = min(max_by_deg, max_by_order)

        for e in range(max_exp + 1):
            exponents[(a, b)] = e
            _recurse(idx + 1, exponents, cur_deg + e, cur_order + e * order)

    _recurse(0, {}, 0, 0)

    # Sort by (total_order, total_degree)
    results.sort(key=lambda x: (x[2], x[1]))
    return results


def _build_polyjet_weak(f, Lx, Ly, max_order=4, max_degree=3,
                        m=8, p=4, use_ibp=False):
    """Build the weak-form polynomial jet library.

    Computes each jet monomial pointwise and convolves with ψ.  When
    ``use_ibp=True``, IBP is applied to single-derivative terms of the
    form ∂_x^a ∂_y^b f (i.e., degree-1 monomials with order > 0).

    Parameters
    ----------
    f : ndarray of shape (Nx, Ny)
    Lx, Ly : float
    max_order  : int — maximum total differential order (default 4)
    max_degree : int — maximum total polynomial degree (default 3)
    m, p : int — test function parameters
    use_ibp : bool — enable IBP for single-derivative columns

    Returns
    -------
    Theta_weak : ndarray of shape (N_valid, K)
    descriptions : list of K str
    """
    from .weak_form import (
        _cached_kernel_1d, _test_function_deriv_1d, weak_convolution,
    )

    Nx, Ny = f.shape
    kernel = _cached_kernel_1d(m, p)

    # Enumerate jet variables
    jet_vars = []  # ((a,b), order, name)
    for total in range(max_order + 1):
        for a in range(total + 1):
            b = total - a
            name = (("x" * a) + ("y" * b)) if total > 0 else ""
            name = f"f_{name}" if total > 0 else "f"
            jet_vars.append(((a, b), total, name))

    # Compute jet fields pointwise
    method = _derivatives.get_backend()

    jet_fields = {}
    jet_names = {}
    if method == "spectral":
        kx = 2.0 * np.pi / Lx * np.fft.fftfreq(Nx) * Nx
        ky = 2.0 * np.pi / Ly * np.fft.fftfreq(Ny) * Ny
        KX, KY = np.meshgrid(kx, ky, indexing="ij")
        fh = np.fft.fft2(f)

    for (a, b), order, name in jet_vars:
        jet_names[(a, b)] = name
        if a + b == 0:
            jet_fields[(a, b)] = f.copy()
        elif method == "spectral":
            d_hat = (1j * KX) ** a * (1j * KY) ** b * fh
            jet_fields[(a, b)] = np.real(np.fft.ifft2(d_hat))
        else:
            jet_fields[(a, b)] = poly_partial(f, Lx, Ly, a, b)

    # Enumerate monomials
    monomials = enumerate_jet_monomials(jet_vars, max_order, max_degree)

    # Build weak columns
    cols, descs = [], []
    for exponents, deg, order in monomials:
        # Check if this is a single-derivative term (IBP candidate)
        non_one_terms = [(k, v) for k, v in exponents.items() if v > 0]
        is_single_deriv = (
            use_ibp
            and len(non_one_terms) == 1
            and non_one_terms[0][0][0] + non_one_terms[0][0][1] > 0
            and non_one_terms[0][1] == 1
        )

        if is_single_deriv:
            (a, b), _exp = non_one_terms[0]
            # IBP: ∫ ψ · ∂_x^a ∂_y^b f = (-1)^{a+b} ∫ ∂_x^a ∂_y^b ψ · f
            # The (-1)^{a+b} sign is handled implicitly by np.convolve.
            kx_ibp = _test_function_deriv_1d(m, p, a)
            ky_ibp = _test_function_deriv_1d(m, p, b)
            dx_a = (Lx / Nx) ** a
            dy_b = (Ly / Ny) ** b
            kx_scaled = kx_ibp / dx_a
            ky_scaled = ky_ibp / dy_b
            col_weak = weak_convolution(f, kx_scaled, ky_scaled)
        else:
            # Build pointwise monomial
            col_raw = np.ones_like(f)
            for (a, b), exp in exponents.items():
                if exp > 0:
                    col_raw = col_raw * (jet_fields[(a, b)] ** exp)
            col_weak = weak_convolution(col_raw, kernel, kernel)

        cols.append(col_weak.ravel())

        # Build description
        parts = []
        for (a, b), exp in sorted(exponents.items(),
                                   key=lambda kv: (kv[0][0] + kv[0][1], kv[0][0])):
            if exp == 0:
                continue
            var_name = jet_names[(a, b)]
            if exp == 1:
                parts.append(var_name)
            else:
                parts.append(f"{var_name}{_sup(exp)}")
        descs.append("1" if not parts else "·".join(parts))

    return np.column_stack(cols), descs
