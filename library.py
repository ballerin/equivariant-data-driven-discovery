"""
Library construction for equivariant PDE discovery.

Evaluates the candidate operator bases (equivariant and standard) on a 2-D
scalar field snapshot, producing the design matrix Θ used in sparse regression.

Provides
--------
build_equivariant_library  — evaluates the O(2)-equivariant multigraph operator
                             basis via get_operators().  Respects the active
                             derivative backend (spectral or polynomial) set in
                             derivatives.py.
build_standard_library     — coordinate-based library (all partial derivatives up
                             to a given order × polynomial prefactors).  Uses the
                             active derivative backend set in derivatives.py
                             (spectral or polynomial) via derivatives.get_backend().
build_standard_library_1d  — same as above for 1-D periodic fields.

Sparse regression and utilities live in separate modules:
  regression.py  — STRidge, TrainSTRidge
  utils.py       — print_pde, rotation_90
  derivatives.py — set_backend / use_backend (spectral ↔ polynomial)
"""

import numpy as np

import derivatives as _derivatives
from multigraph_operators import get_operators
from polynomial_derivatives import poly_partial, poly_partial_1d


def _wavenumbers_2d(Nx, Ny, Lx, Ly):
    """Return 2-D wavenumber grids (KX, KY) for an (Nx, Ny) FFT grid."""
    kx = 2.0 * np.pi / Lx * np.fft.fftfreq(Nx) * Nx
    ky = 2.0 * np.pi / Ly * np.fft.fftfreq(Ny) * Ny
    return np.meshgrid(kx, ky, indexing='ij')


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
    if method == 'spectral':
        KX, KY = _wavenumbers_2d(Nx, Ny, Lx, Ly)
        fh = np.fft.fft2(f)

    multi_indices = []
    for total in range(max_order + 1):
        for a in range(total + 1):
            multi_indices.append((a, total - a))

    derivs, deriv_names = [], []
    for a, b in multi_indices:
        if method == 'spectral':
            d_hat = (1j * KX)**a * (1j * KY)**b * fh
            deriv = np.real(np.fft.ifft2(d_hat))
        else:
            if a + b == 0:
                deriv = f.copy()
            else:
                deriv = poly_partial(f, Lx, Ly, a, b)
        name = (('x' * a) + ('y' * b)) if (a + b) > 0 else '1'
        name = f'f_{name}' if (a + b) > 0 else '1'
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
                desc = f'f·{dname}' if dname != '1' else 'f'
            else:
                desc = f'f^{p}·{dname}' if dname != '1' else f'f^{p}'
            cols.append(col)
            descs.append(desc)

    return np.column_stack(cols), descs


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
    deriv_names = ['1']
    if method == 'spectral':
        k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
        uhat = np.fft.fft(u)
        for d in range(1, max_deriv + 1):
            derivs.append(np.real(np.fft.ifft((1j * k) ** d * uhat)))
            deriv_names.append('u_' + 'x' * d)
    else:
        for d in range(1, max_deriv + 1):
            derivs.append(poly_partial_1d(u, L, d))
            deriv_names.append('u_' + 'x' * d)

    cols, descs = [], []
    for d_arr, d_name in zip(derivs, deriv_names):
        for p in range(max_poly + 1):
            col = u ** p * d_arr
            if p == 0:
                desc = d_name
            elif p == 1:
                desc = f'u*{d_name}' if d_name != '1' else 'u'
            else:
                desc = f'u^{p}*{d_name}' if d_name != '1' else f'u^{p}'
            cols.append(col)
            descs.append(desc)

    return np.column_stack(cols), descs
