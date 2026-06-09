"""
Library construction for equivariant PDE discovery.

Evaluates the candidate operator bases (equivariant and standard) on a 2-D
scalar field snapshot, producing the design matrix Θ used in sparse regression.

Provides
--------
build_equivariant_library    — O(2)-equivariant multigraph operator basis via
                               get_operators().  Respects the active derivative
                               backend (spectral or polynomial).
build_standard_library       — coordinate-based library (all partial derivatives
                               up to a given order × polynomial prefactors),
                               i.e. the PDE-FIND library f^p · ∂^a f.
build_wsindy_library         — WSINDy library ∂^α(f^q): partial derivatives
                               applied to monomial powers of f.
build_polynomial_jet_library — all monomials in jet variables (f, f_x, f_y,
                               f_{xx}, …) up to given degree and order.
build_standard_library_1d    — same as build_standard_library for 1-D fields.

Sparse regression and utilities live in separate modules:
  regression.py  — STRidge, TrainSTRidge
  utils.py       — print_pde, rotation_90
  derivatives.py — set_backend / use_backend (spectral ↔ polynomial)
"""

import numpy as np

import derivatives as _derivatives


# Unicode superscript digits for polynomial exponents in descriptions.
_SUP = {2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}


def _sup(n: int) -> str:
    """Return a Unicode superscript string for integer n >= 2."""
    return _SUP.get(n, f"^{n}")
from multigraph_operators import get_operators
from polynomial_derivatives import poly_partial, poly_partial_1d


def _wavenumbers_2d(Nx, Ny, Lx, Ly):
    """Return 2-D wavenumber grids (KX, KY) for an (Nx, Ny) FFT grid."""
    kx = 2.0 * np.pi / Lx * np.fft.fftfreq(Nx) * Nx
    ky = 2.0 * np.pi / Ly * np.fft.fftfreq(Ny) * Ny
    return np.meshgrid(kx, ky, indexing="ij")


def spectral_lowpass(f, Lx, Ly, cutoff):
    """Apply a 2-D isotropic spectral low-pass filter.

    Zeroes all Fourier modes with physical wavenumber magnitude |k| > cutoff
    and returns the real-space reconstruction.  The filter is O(2)-symmetric
    (isotropic cutoff sphere), so it commutes with all O(2)-equivariant
    operators in the library.

    Parameters
    ----------
    f      : ndarray, shape (Nx, Ny)
    Lx, Ly : float — domain lengths
    cutoff : float — maximum physical wavenumber to retain (e.g. 3.0)

    Returns
    -------
    f_filtered : ndarray, shape (Nx, Ny), real
    """
    Nx, Ny = f.shape
    fh = np.fft.fft2(f)
    KX, KY = _wavenumbers_2d(Nx, Ny, Lx, Ly)
    fh[KX ** 2 + KY ** 2 > cutoff ** 2] = 0.0
    return np.real(np.fft.ifft2(fh))


# -----------------------------------------------------------------------
# Library builders
# -----------------------------------------------------------------------


def build_equivariant_library(f, Lx, Ly, max_p=2, max_n=3):
    """Build the equivariant library matrix for a 2-D scalar field f.

    Parameters
    ----------
    f     : ndarray, shape (Nx, Ny)
    Lx    : float — domain length in x
    Ly    : float — domain length in y
    max_p : int — maximum multigraph edge count p (differential order = 2*max_p)
    max_n : int — maximum polynomial degree (total number of f-factors per term)

    Returns
    -------
    Theta        : ndarray, shape (Nx*Ny, K)
    descriptions : list of K strings

    Examples
    --------
    >>> Theta, descs = build_equivariant_library(f, Lx, Ly)           # p<=2, n<=3 → 20 cols
    >>> Theta, descs = build_equivariant_library(f, Lx, Ly, max_n=4)  # p<=2, n<=4 → 30 cols
    >>> Theta, descs = build_equivariant_library(f, Lx, Ly, max_p=3)  # p<=3, n<=3 → 40 cols
    >>> Theta, descs = build_equivariant_library(f, Lx, Ly, max_p=3, max_n=4)  # → 68 cols
    """
    ops = get_operators(max_p, max_n)
    cols, descs = [], []
    for op, name in ops:
        cols.append(op(f, Lx, Ly).ravel())
        descs.append(name)
    return np.column_stack(cols), descs


def build_standard_library(f, Lx, Ly, max_order=4, max_poly=3):
    """Build a coordinate-based library as used in PDE-Find.

    Columns = (all partial derivatives of f up to order max_order) ×
              (monomials f^0, f^1, ..., f^max_poly).
    For max_order=4 and max_poly=3 this gives 15 × 4 = 60 columns.

    The differentiation method is determined by the active backend in
    ``derivatives.py`` (set via ``derivatives.set_backend()`` or
    ``derivatives.use_backend()``):

    - ``'spectral'`` (default): FFT-based exact differentiation.  Accurate for
      noiseless periodic data; matches the original PDE-FIND Fourier option.
    - ``'polynomial'``: local 2-D polynomial least-squares fitting via
      ``poly_partial`` from ``polynomial_derivatives.py``.  Noise-robust;
      comparable to PDE-FIND's ``space_diff='poly'`` but uses a joint 2-D
      patch fit instead of separable 1-D Chebyshev fits, which avoids error
      compounding for mixed partials.

    Parameters
    ----------
    f         : ndarray, shape (Nx, Ny)
    Lx, Ly    : float — domain lengths
    max_order : int — highest partial derivative order (default 4)
    max_poly  : int — highest polynomial degree factor (default 3)

    Returns
    -------
    Theta        : ndarray, shape (Nx*Ny, n_cols)
    descriptions : list of strings
    """
    method = _derivatives.get_backend()
    Nx, Ny = f.shape
    if method == "spectral":
        KX, KY = _wavenumbers_2d(Nx, Ny, Lx, Ly)
        fh = np.fft.fft2(f)

    multi_indices = []
    for total in range(max_order + 1):
        for a in range(total + 1):
            multi_indices.append((a, total - a))

    derivs, deriv_names = [], []
    for a, b in multi_indices:
        if a + b == 0:
            deriv = np.ones_like(f)
        elif method == "spectral":
            d_hat = (1j * KX) ** a * (1j * KY) ** b * fh
            deriv = np.real(np.fft.ifft2(d_hat))
        else:
            deriv = poly_partial(f, Lx, Ly, a, b)
        name = (("x" * a) + ("y" * b)) if (a + b) > 0 else "1"
        name = f"f_{name}" if (a + b) > 0 else "1"
        derivs.append(deriv)
        deriv_names.append(name)

    cols, descs = [], []
    for p in range(max_poly + 1):
        for deriv, dname in zip(derivs, deriv_names):
            poly = f**p if p > 0 else np.ones_like(f)
            col = (poly * deriv).ravel()
            if p == 0:
                desc = dname
            elif p == 1:
                desc = f"f·{dname}" if dname != "1" else "f"
            else:
                desc = f"f{_sup(p)}·{dname}" if dname != "1" else f"f{_sup(p)}"
            cols.append(col)
            descs.append(desc)

    return np.column_stack(cols), descs


def build_wsindy_library(f, Lx, Ly, max_order=4, max_poly=3):
    """Build the WSINDy library: ∂_x^a ∂_y^b (f^q).

    Columns = all derivative combinations ∂_x^a ∂_y^b applied to monomial
    powers f^q, for 0 ≤ a+b ≤ max_order and 1 ≤ q ≤ max_poly, plus the
    constant term 1 (a=b=0, q=0).

    This is the pointwise analogue of the library used in Messenger & Bortz
    (2021).  For max_order=4 and max_poly=3 this gives
    1 + (15 × 3) = 46 columns.

    The differentiation method respects the active backend in
    ``derivatives.py``.

    Parameters
    ----------
    f         : ndarray, shape (Nx, Ny)
    Lx, Ly    : float — domain lengths
    max_order : int — highest partial derivative order (default 4)
    max_poly  : int — highest monomial power q in ∂^α(f^q) (default 3)

    Returns
    -------
    Theta        : ndarray, shape (Nx*Ny, n_cols)
    descriptions : list of n_cols strings
    """
    method = _derivatives.get_backend()
    Nx, Ny = f.shape
    if method == "spectral":
        KX, KY = _wavenumbers_2d(Nx, Ny, Lx, Ly)

    # Enumerate derivative multi-indices (same as standard library)
    multi_indices = []
    for total in range(max_order + 1):
        for a in range(total + 1):
            multi_indices.append((a, total - a))

    cols, descs = [], []

    # Constant term: a=b=0, q=0 → ∂^0(f^0) = 1
    cols.append(np.ones(f.size, dtype=np.float64))
    descs.append("1")

    for q in range(1, max_poly + 1):
        # Compute f^q once per polynomial power
        f_pow = f ** q
        if method == "spectral":
            fh_pow = np.fft.fft2(f_pow)

        for a, b in multi_indices:
            if a + b == 0:
                # ∂^0(f^q) = f^q — pure monomial, handled separately below
                # for q≥2 to avoid duplicate with q=1 case
                continue

            if method == "spectral":
                d_hat = (1j * KX) ** a * (1j * KY) ** b * fh_pow
                deriv = np.real(np.fft.ifft2(d_hat))
            else:
                deriv = poly_partial(f_pow, Lx, Ly, a, b)
            cols.append(deriv.ravel())

            # Build description: ∂_x^a ∂_y^b (f^q)
            dx_part = "x" * a + "y" * b
            if q == 1:
                desc = f"f_{{{dx_part}}}"
            else:
                desc = f"(f{_sup(q)})_{{{dx_part}}}"
            descs.append(desc)

    # Add pure monomials f^q for q ≥ 1
    for q in range(1, max_poly + 1):
        cols.append((f ** q).ravel())
        descs.append("f" if q == 1 else f"f{_sup(q)}")

    return np.column_stack(cols), descs


def build_polynomial_jet_library(f, Lx, Ly, max_order=4, max_degree=3):
    """Build the polynomial jet library: all monomials in jet variables.

    The *jet variables* are all partial derivatives ∂_x^a ∂_y^b f for
    0 ≤ a+b ≤ max_order (including f itself).  A library column is a
    product (monomial) of these jet variables:

        ∏_{a,b} (∂_x^a ∂_y^b f)^{e_{a,b}}

    subject to:
      - total degree:  Σ e_{a,b} ≤ max_degree
      - total differential order:  Σ e_{a,b} · (a+b) ≤ max_order

    For max_order=4 and max_degree=3 this gives ~105 columns — the full
    space of polynomial nonlinear differential operators in coordinates,
    compared to 20 O(2)-equivariant operators.

    Parameters
    ----------
    f          : ndarray, shape (Nx, Ny)
    Lx, Ly     : float — domain lengths
    max_order  : int — maximum total differential order (default 4)
    max_degree : int — maximum total polynomial degree (default 3)

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
    monomials = _enumerate_jet_monomials(jet_vars, max_order, max_degree)

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


def _enumerate_jet_monomials(jet_vars, max_order, max_degree):
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


def build_standard_library_1d(u, L, max_deriv=3, max_poly=3):
    """Build a coordinate-based library for a 1D periodic scalar field.

    Columns = (derivatives ∂^d u for d = 0, ..., max_deriv) ×
              (monomials u^p for p = 0, ..., max_poly).

    The differentiation method is determined by the active backend in
    ``derivatives.py`` (set via ``derivatives.set_backend()`` or
    ``derivatives.use_backend()``):

    - ``'spectral'`` (default): FFT-based exact differentiation.  Matches
      PDE-FIND's ``space_diff='Fourier'``.
    - ``'polynomial'``: local 1-D polynomial least-squares fitting via
      ``poly_partial_1d``.  Noise-robust, consistent with the 2-D polynomial
      backend.

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
