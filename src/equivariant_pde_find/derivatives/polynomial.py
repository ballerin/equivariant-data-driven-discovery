"""
Polynomial (local least-squares) derivative utilities for scalar fields on
periodic N-D domains.

Provides the **same public names** as spectral.py so that the two modules
are interchangeable.  To switch from spectral to polynomial differentiation:

    from equivariant_pde_find import derivatives
    derivatives.set_backend('polynomial')  # switch to polynomial

Method
------
At each grid point an ndim-dimensional polynomial of total degree `deg` is
fitted by least-squares to the ``(2r+1)^ndim`` periodic neighbourhood.  The
required partial derivative is extracted from the polynomial coefficients:

    ∂^α f / ∂x^α  ≈  α! * c_α

where ``c_α`` is the coefficient of the monomial ``∏ x_i^{α_i}`` in the
local fit and ``α! = ∏ α_i!``.

Because the grid is uniform and periodic, the fitting reduces to a **fixed
linear filter** (a ``(2r+1)^ndim`` kernel) applied to the full field via
periodic ndim correlation.  The kernel is computed once via the
pseudoinverse of the local Vandermonde matrix.  Isotropic domain: like
spectral.py, ``ndim<=3`` supports independent per-axis (Lx, Ly, Lz) and grid
sizes; ``ndim>3`` assumes an isotropic domain (single length, reused for
every axis), matching every N-D simulator in this codebase.

Cost note
---------
The stencil has ``(2r+1)^ndim`` points and the Vandermonde pseudoinverse is
``O(n_pts * n_mono^2)`` — both grow quickly with ndim.  The per-order
``(deg, r)`` defaults are unchanged from the 2-D backend; for high ndim this
means a slower (but still correct) fit.  Expect the polynomial backend to
become impractical well before spectral does, as ndim grows.

Default stencil parameters (derivative order → (deg, r)):
  Order 1  (grad)         : deg=9,  r=5
  Order 2  (Hessian, Δ)   : deg=9,  r=5
  Order 3  (3rd partials) : deg=9,  r=5
  Order 4  (Δ²)           : deg=11, r=6
  Order 6  (Δ³)           : deg=13, r=7

Noise robustness vs. spectral
------------------------------
Spectral derivatives are exact for band-limited, noiseless data.  Polynomial
fitting introduces a small approximation bias (O(dx^{deg+1-order})) but
suppresses noise: fitting a low-degree polynomial to many neighbours is
equivalent to implicit smoothing before differentiating.  The bias-variance
trade-off is controlled by deg and r.  For clean pseudospectral data, prefer
spectral.py; add this module when Gaussian noise is injected.

Note: sixth-order polynomial fitting (trilaplacian / Δ³f) is intrinsically
sensitive to noise.  When working with noisy data, consider restricting to
max_p ≤ 2 to avoid the highest-order operators.
"""

import itertools
from functools import lru_cache
from math import factorial

import numpy as np

from ..utils import _apply_prefilter
from ._common import axis_lengths

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Default (deg, r) per derivative order.
# Invariant: (2r+1)^ndim  >=  #monomials of degree <= deg in ndim variables.
# These defaults are chosen to give <0.1% relative error at wavenumber k*dx ~ 0.3
# (e.g. k=3 on N=64, L=2π), which covers typical pseudospectral PDE simulations.
# Accuracy scales as (k·dx)^(deg+1−order); larger deg gives exponentially better
# results at the cost of a wider stencil (still O(N^ndim) via periodic correlation).
_DEFAULTS = {
    1: (9, 5),   # <1e-5 rel. err at k=3 (2-D: 121 pts, 55 mono)
    2: (9, 5),   # <1e-5 rel. err at k=3
    3: (9, 5),   # adequate for 3rd-order terms
    4: (11, 6),  # good for 4th-order (bilaplacian)
    6: (13, 7),  # reasonable for 6th-order (trilaplacian)
}


def _enumerate_monomials(ndim, deg):
    """All ndim-tuples (p_1, ..., p_ndim) of nonneg ints with sum <= deg.

    Ordered so that lower-total-degree monomials come first (matches the
    2-D ``[(p, q) for p in range(deg+1) for q in range(deg+1-p)]`` order
    exactly when ndim=2).
    """
    def _rec(ndim_left, budget):
        if ndim_left == 1:
            for p in range(budget + 1):
                yield (p,)
            return
        for p in range(budget + 1):
            for rest in _rec(ndim_left - 1, budget - p):
                yield (p,) + rest

    return list(_rec(ndim, deg))


@lru_cache(maxsize=256)
def _poly_kernel_nd(multi_index, deg, r, dxs):
    """Build the ``(2r+1)^ndim`` correlation kernel for ``∂^multi_index f``.

    Parameters
    ----------
    multi_index : tuple of int, length ndim — target derivative multi-index
    deg  : int — polynomial degree to fit
    r    : int — stencil half-width (isotropic, same for every axis)
    dxs  : tuple of float, length ndim — grid spacing per axis

    Returns
    -------
    kernel : ndarray, shape (2r+1,)*ndim
        Applying periodic ndim correlation of the field with this kernel
        yields the requested partial derivative at every grid point.
    """
    ndim = len(multi_index)
    offsets = list(range(-r, r + 1))
    pts = list(itertools.product(offsets, repeat=ndim))  # (2r+1)^ndim rows
    monomials = _enumerate_monomials(ndim, deg)

    # Vandermonde matrix: V[pt, mono] = prod_axis (offset_axis * dx_axis)^{mono_axis}
    V = np.empty((len(pts), len(monomials)), dtype=float)
    for row, pt in enumerate(pts):
        for col, mono in enumerate(monomials):
            val = 1.0
            for axis in range(ndim):
                val *= (pt[axis] * dxs[axis]) ** mono[axis]
            V[row, col] = val

    Vp = np.linalg.pinv(V)  # shape (n_mono, n_pts)

    target = monomials.index(multi_index)
    coeff = 1
    for a in multi_index:
        coeff *= factorial(a)
    kernel = (Vp[target] * coeff).reshape((2 * r + 1,) * ndim)
    return kernel


def _apply_nd(f, kernel):
    """Periodic ndim correlation of f with kernel (pure NumPy, no scipy).

    Generalizes the 2-D ``np.roll``-based correlation to any ndim.
    """
    ndim = f.ndim
    r = kernel.shape[0] // 2
    offsets = list(range(-r, r + 1))
    result = np.zeros_like(f, dtype=float)
    for offset_tuple in itertools.product(offsets, repeat=ndim):
        weight = kernel[tuple(o + r for o in offset_tuple)]
        if weight == 0.0:
            continue
        shifted = f
        for axis, o in enumerate(offset_tuple):
            if o:
                shifted = np.roll(shifted, -o, axis=axis)
        result += weight * shifted
    return result


def _deriv_nd(f, Ls, multi_index, deg, r):
    """Evaluate ∂^multi_index f at all grid points via polynomial fitting."""
    dxs = tuple(L / N for L, N in zip(Ls, f.shape))
    return _apply_nd(f, _poly_kernel_nd(tuple(multi_index), deg, r, dxs))


# Default (deg, r) for derivative orders not listed in _DEFAULTS.
# Rule of thumb: deg = order + 1, r = ceil((deg + 1) / 2).
def _default_stencil(order):
    if order in _DEFAULTS:
        return _DEFAULTS[order]
    deg = order + 1
    r = (deg + 2) // 2
    return deg, r


# ------------------ runtime stencil overrides ------------------
# Users can call `set_default_stencil(deg, r, order=None)` to override the
# default (deg, r) used by the named functions (polynomial_grad, ...).
# If `order` is None the provided pair is used for all orders; otherwise it
# sets the pair for that specific total derivative order (a+b).
_STENCIL_OVERRIDE_GLOBAL = None
_STENCIL_OVERRIDE_BY_ORDER = {}


def set_default_stencil(deg=None, r=None, order=None):
    """Override the default polynomial stencil used by the convenience calls.

    Parameters
    ----------
    deg, r : int or None
        If both provided, set the override to (deg, r). If either is None,
        the override for the target `order` (or global) is removed.
    order : int or None
        If provided, set the override only for that total derivative order
        (a+b). If None, set a global override used for all orders.
    """
    global _STENCIL_OVERRIDE_GLOBAL, _STENCIL_OVERRIDE_BY_ORDER
    if order is None:
        if deg is None or r is None:
            _STENCIL_OVERRIDE_GLOBAL = None
        else:
            _STENCIL_OVERRIDE_GLOBAL = (int(deg), int(r))
    else:
        if deg is None or r is None:
            _STENCIL_OVERRIDE_BY_ORDER.pop(order, None)
        else:
            _STENCIL_OVERRIDE_BY_ORDER[order] = (int(deg), int(r))


def get_default_stencil(order):
    """Return the (deg, r) to use for the requested total derivative order.

    Checks per-order overrides first, then a global override, then built-in
    defaults, and finally falls back to the heuristic `_default_stencil`.
    """
    if order in _STENCIL_OVERRIDE_BY_ORDER:
        return _STENCIL_OVERRIDE_BY_ORDER[order]
    if _STENCIL_OVERRIDE_GLOBAL is not None:
        return _STENCIL_OVERRIDE_GLOBAL
    if order in _DEFAULTS:
        return _DEFAULTS[order]
    return _default_stencil(order)


# Canonical per-axis domain lengths; see derivatives/_common.py.
_axis_lengths = axis_lengths


def _unit_multi_index(ndim, axis_tuple):
    mi = [0] * ndim
    for a in axis_tuple:
        mi[a] += 1
    return tuple(mi)


# ---------------------------------------------------------------------------
# Core N-D primitive
# ---------------------------------------------------------------------------


def poly_partial_nd(f, Ls, multi_index, deg=None, r=None):
    """Compute the partial derivative ∂^multi_index f via local ndim
    polynomial fitting.

    Parameters
    ----------
    f        : ndarray, shape (N_1, ..., N_ndim)
    Ls       : sequence of float, length f.ndim — per-axis domain length
    multi_index : sequence of int, length f.ndim
    deg, r   : int, optional — override the default stencil for this order

    Returns
    -------
    ndarray, shape f.shape
    """
    multi_index = tuple(multi_index)
    if len(multi_index) != f.ndim:
        raise ValueError(
            f"multi_index has length {len(multi_index)}, expected {f.ndim} (f.ndim)"
        )
    order = sum(multi_index)
    if order == 0:
        return f.copy()
    if deg is None or r is None:
        _deg, _r = get_default_stencil(order)
        if deg is None:
            deg = _deg
        if r is None:
            r = _r
    return _deriv_nd(f, Ls, multi_index, deg, r)


# ---------------------------------------------------------------------------
# Public APIs (backward-compatible (Lx, Ly, Lz) signature, now N-D-generic)
# ---------------------------------------------------------------------------


def polynomial_grad(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Gradient via local polynomial fit, for any ndim.

    Returns a tuple of f.ndim components (f_x, f_y, ...).
    """
    f = _apply_prefilter(f, Lx, Ly, prefilter_sigma)
    ndim = f.ndim
    Ls = _axis_lengths(f, Lx, Ly, Lz)
    return tuple(
        poly_partial_nd(f, Ls, _unit_multi_index(ndim, (axis,)))
        for axis in range(ndim)
    )


def polynomial_laplacian(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Laplacian via local polynomial fit, for any ndim: Δf = Σ_axis f_{axis,axis}."""
    f = _apply_prefilter(f, Lx, Ly, prefilter_sigma)
    ndim = f.ndim
    Ls = _axis_lengths(f, Lx, Ly, Lz)
    return sum(
        poly_partial_nd(f, Ls, _unit_multi_index(ndim, (axis, axis)))
        for axis in range(ndim)
    )


def polynomial_hessian(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Hessian via local polynomial fit, for any ndim.

    Same canonical unique-component ordering as ``spectral_hessian``.
    """
    f = _apply_prefilter(f, Lx, Ly, prefilter_sigma)
    ndim = f.ndim
    Ls = _axis_lengths(f, Lx, Ly, Lz)
    pairs = list(itertools.combinations_with_replacement(range(ndim), 2))
    return tuple(
        poly_partial_nd(f, Ls, _unit_multi_index(ndim, pair)) for pair in pairs
    )


def polynomial_bilaplacian(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Bilaplacian via local polynomial fit, for any ndim:
    Δ²f = Σ_i Σ_j f_{iijj} (i.e. Σ_i f_iiii + 2 Σ_{i<j} f_iijj)."""
    f = _apply_prefilter(f, Lx, Ly, prefilter_sigma)
    ndim = f.ndim
    Ls = _axis_lengths(f, Lx, Ly, Lz)
    total = np.zeros_like(f, dtype=float)
    for i in range(ndim):
        for j in range(ndim):
            mi = _unit_multi_index(ndim, (i, i, j, j))
            total += poly_partial_nd(f, Ls, mi)
    return total


def polynomial_trilaplacian(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Trilaplacian via local polynomial fit, for any ndim: Δ³f = Σ_i Σ_j Σ_k f_{iijjkk}."""
    f = _apply_prefilter(f, Lx, Ly, prefilter_sigma)
    ndim = f.ndim
    Ls = _axis_lengths(f, Lx, Ly, Lz)
    total = np.zeros_like(f, dtype=float)
    for i in range(ndim):
        for j in range(ndim):
            for k in range(ndim):
                mi = _unit_multi_index(ndim, (i, i, j, j, k, k))
                total += poly_partial_nd(f, Ls, mi)
    return total


def polynomial_third_order(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Third-order partial derivatives via local polynomial fit, for any ndim.

    Same canonical unique-component ordering as ``spectral_third_order``.
    """
    f = _apply_prefilter(f, Lx, Ly, prefilter_sigma)
    ndim = f.ndim
    Ls = _axis_lengths(f, Lx, Ly, Lz)
    triples = list(itertools.combinations_with_replacement(range(ndim), 3))
    return tuple(
        poly_partial_nd(f, Ls, _unit_multi_index(ndim, t)) for t in triples
    )


def poly_partial(f, Lx, Ly, a, b, deg=None, r=None):
    """Compute the mixed partial derivative ∂^(a+b) f / ∂x^a ∂y^b via local
    2-D polynomial fitting.

    This is the 2-D convenience interface (used historically by the 2-D
    library builders); for N-D use ``poly_partial_nd``.

    Parameters
    ----------
    f        : ndarray, shape (Nx, Ny)
    Lx, Ly   : float — domain lengths
    a, b     : int ≥ 0 — derivative orders in x and y
    deg      : int, optional — polynomial degree to fit locally.
               Defaults to the value in _DEFAULTS for order a+b, or
               order+1 if not listed.
    r        : int, optional — stencil half-width ((2r+1)² > #monomials required).
               Defaults to the value in _DEFAULTS for order a+b.

    Returns
    -------
    ndarray, shape (Nx, Ny)

    Examples
    --------
    >>> fx  = poly_partial(f, Lx, Ly, 1, 0)   # ∂f/∂x
    >>> fxy = poly_partial(f, Lx, Ly, 1, 1)   # ∂²f/∂x∂y
    >>> fxxy = poly_partial(f, Lx, Ly, 2, 1)  # ∂³f/∂x²∂y
    """
    return poly_partial_nd(f, (Lx, Ly), (a, b), deg=deg, r=r)


def poly_partial_1d(u, L, d, deg=None, r=None):
    """Compute the d-th derivative of a 1-D periodic field via local polynomial fitting.

    Uses the same ndim kernel approach restricted to 1-D: fits a degree-`deg`
    polynomial on a (2r+1)-point periodic stencil and extracts the d-th
    derivative coefficient.

    Parameters
    ----------
    u   : ndarray, shape (N,)
    L   : float — domain length
    d   : int ≥ 1 — derivative order
    deg : int, optional — polynomial degree (default: d+1 or from _DEFAULTS)
    r   : int, optional — stencil half-width (default: from _DEFAULTS)

    Returns
    -------
    ndarray, shape (N,)
    """
    return poly_partial_nd(u, (L,), (d,), deg=deg, r=r)
