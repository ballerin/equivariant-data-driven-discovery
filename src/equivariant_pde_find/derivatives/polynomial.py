"""
Polynomial (local least-squares) derivative utilities for scalar fields on periodic 2-D domains.

Provides the **same public names** as spectral_derivatives.py so that the two
modules are interchangeable.  To switch from spectral to polynomial
differentiation, change the one import block in multigraph_operators.py::

    import derivatives
    derivatives.set_backend('polynomial')  # switch to polynomial

Method
------
At each grid point (i, j) a 2-D polynomial of total degree `deg` is fitted by
least-squares to the (2r+1)×(2r+1) periodic neighbourhood.  The required
partial derivative is extracted from the polynomial coefficients:

    ∂^(a+b) f / ∂x^a ∂y^b  ≈  a! b! * c_{ab}

where c_{ab} is the coefficient of x^a y^b in the local fit.

Because the grid is uniform and periodic, the fitting reduces to a **fixed
linear filter** (a (2r+1)×(2r+1) kernel) applied to the full field via
periodic 2-D correlation.  The kernel is computed once via the pseudoinverse
of the local Vandermonde matrix.  All pixels are processed simultaneously in
O(N² × (2r+1)²) time — no per-pixel matrix solve.

Default stencil parameters (derivative order → (deg, r, stencil)):
  Order 1  (grad)         : deg=3, r=2  →  5×5  stencil, 10 monomials / 25 points
  Order 2  (Hessian, Δ)   : deg=3, r=2  →  5×5  stencil, 10 monomials / 25 points
  Order 3  (3rd partials) : deg=4, r=2  →  5×5  stencil, 15 monomials / 25 points
  Order 4  (Δ²)           : deg=5, r=3  →  7×7  stencil, 21 monomials / 49 points
  Order 6  (Δ³)           : deg=7, r=4  →  9×9  stencil, 36 monomials / 81 points

Noise robustness vs. spectral
------------------------------
Spectral derivatives are exact for band-limited, noiseless data.  Polynomial
fitting introduces a small approximation bias (O(dx^{deg+1-order})) but
suppresses noise: fitting a low-degree polynomial to many neighbours is
equivalent to implicit smoothing before differentiating.  The bias-variance
trade-off is controlled by deg and r.  For clean pseudospectral data, prefer
spectral_derivatives.py; add this module when Gaussian noise is injected.

Note: sixth-order polynomial fitting (trilaplacian / Δ³f) is intrinsically
sensitive to noise.  When working with noisy data, consider restricting to
max_p ≤ 2 to avoid the highest-order operators.
"""

from functools import lru_cache
from math import factorial

import numpy as np

from ..utils import _apply_prefilter

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Default (deg, r) per derivative order.
# Invariant: (2r+1)^2  >=  (deg+1)(deg+2)/2   (more points than monomials).
# These defaults are chosen to give <0.1% relative error at wavenumber k*dx ~ 0.3
# (e.g. k=3 on N=64, L=2π), which covers typical pseudospectral PDE simulations.
# Accuracy scales as (k·dx)^(deg+1−order); larger deg gives exponentially better
# results at the cost of a wider stencil (still O(N²) via periodic correlation).
_DEFAULTS = {
    1: (9, 5),   # 121 pts,  55 mono — <1e-5 rel. err at k=3
    2: (9, 5),   # 121 pts,  55 mono — <1e-5 rel. err at k=3
    3: (9, 5),   # 121 pts,  55 mono — adequate for 3rd-order terms
    4: (11, 6),  # 169 pts,  78 mono — good for 4th-order (bilaplacian)
    6: (13, 7),  # 225 pts, 105 mono — reasonable for 6th-order (trilaplacian)
}


@lru_cache(maxsize=256)
def _poly_kernel(a, b, deg, r, dx, dy):
    """Build the (2r+1)×(2r+1) correlation kernel for ∂^(a+b)f/∂x^a∂y^b.

    Parameters
    ----------
    a, b : int   — derivative orders in x and y  (a + b <= deg required)
    deg  : int   — polynomial degree to fit
    r    : int   — stencil half-width  ((2r+1)² must exceed #monomials)
    dx   : float — grid spacing in x  (Lx / Nx)
    dy   : float — grid spacing in y  (Ly / Ny)

    Returns
    -------
    kernel : ndarray, shape (2r+1, 2r+1)
        Applying periodic 2-D correlation of the field with this kernel yields
        the requested partial derivative at every grid point.
    """
    offsets = list(range(-r, r + 1))
    pts = [(di, dj) for di in offsets for dj in offsets]  # (2r+1)² rows
    monomials = [(p, q) for p in range(deg + 1) for q in range(deg + 1 - p)]

    # Vandermonde matrix: V[pt, mono] = (di*dx)^p * (dj*dy)^q
    # Note: 0.0**0 == 1.0 in Python, so the constant monomial is handled correctly.
    V = np.array(
        [[(di * dx) ** p * (dj * dy) ** q for p, q in monomials] for di, dj in pts],
        dtype=float,
    )  # shape (n_pts, n_mono)

    Vp = np.linalg.pinv(V)  # shape (n_mono, n_pts)

    # Row of Vp corresponding to monomial x^a y^b gives the LS coefficient c_{ab}.
    # The derivative is a! b! c_{ab}.
    target = next(i for i, (p, q) in enumerate(monomials) if p == a and q == b)
    kernel = (Vp[target] * factorial(a) * factorial(b)).reshape(2 * r + 1, 2 * r + 1)
    return kernel


def _apply(f, kernel):
    """Periodic 2D correlation of f with kernel (pure NumPy, no scipy).

    Computes output[i,j] = Σ_{di,dj} kernel[r+di, r+dj] * f[(i+di)%Nx, (j+dj)%Ny],
    equivalent to scipy.ndimage.correlate(f, kernel, mode='wrap').
    """
    r = kernel.shape[0] // 2
    result = np.zeros_like(f, dtype=float)
    for di in range(-r, r + 1):
        for dj in range(-r, r + 1):
            result += kernel[di + r, dj + r] * np.roll(np.roll(f, -di, axis=0), -dj, axis=1)
    return result


def _deriv(f, Lx, Ly, a, b, deg, r):
    """Evaluate ∂^(a+b)f/∂x^a∂y^b at all grid points via polynomial fitting."""
    Nx, Ny = f.shape
    dx = Lx / Nx
    dy = Ly / Ny
    return _apply(f, _poly_kernel(a, b, deg, r, dx, dy))


# ---------------------------------------------------------------------------
# Public APIs
# ---------------------------------------------------------------------------


def polynomial_grad(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Gradient via local polynomial fit.

    2-D: returns (f_x, f_y).  3-D: not yet supported.
    """
    if f.ndim == 3:
        raise NotImplementedError("Polynomial backend does not support 3-D fields yet")
    f = _apply_prefilter(f, Lx, Ly, prefilter_sigma)
    # allow user overrides via set_default_stencil
    if 'get_default_stencil' in globals():
        deg, r = get_default_stencil(1)
    else:
        deg, r = _DEFAULTS[1]
    fx = _deriv(f, Lx, Ly, 1, 0, deg, r)
    fy = _deriv(f, Lx, Ly, 0, 1, deg, r)
    return fx, fy


def polynomial_laplacian(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Laplacian via local polynomial fit. 3-D not yet supported."""
    if f.ndim == 3:
        raise NotImplementedError("Polynomial backend does not support 3-D fields yet")
    f = _apply_prefilter(f, Lx, Ly, prefilter_sigma)
    if 'get_default_stencil' in globals():
        deg, r = get_default_stencil(2)
    else:
        deg, r = _DEFAULTS[2]
    return _deriv(f, Lx, Ly, 2, 0, deg, r) + _deriv(f, Lx, Ly, 0, 2, deg, r)


def polynomial_hessian(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Hessian via local polynomial fit. 3-D not yet supported."""
    if f.ndim == 3:
        raise NotImplementedError("Polynomial backend does not support 3-D fields yet")
    f = _apply_prefilter(f, Lx, Ly, prefilter_sigma)
    if 'get_default_stencil' in globals():
        deg, r = get_default_stencil(2)
    else:
        deg, r = _DEFAULTS[2]
    fxx = _deriv(f, Lx, Ly, 2, 0, deg, r)
    fxy = _deriv(f, Lx, Ly, 1, 1, deg, r)
    fyy = _deriv(f, Lx, Ly, 0, 2, deg, r)
    return fxx, fxy, fyy


def polynomial_bilaplacian(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Bilaplacian via local polynomial fit. 3-D not yet supported."""
    if f.ndim == 3:
        raise NotImplementedError("Polynomial backend does not support 3-D fields yet")
    f = _apply_prefilter(f, Lx, Ly, prefilter_sigma)
    if 'get_default_stencil' in globals():
        deg, r = get_default_stencil(4)
    else:
        deg, r = _DEFAULTS[4]
    return (
        _deriv(f, Lx, Ly, 4, 0, deg, r)
        + 2 * _deriv(f, Lx, Ly, 2, 2, deg, r)
        + _deriv(f, Lx, Ly, 0, 4, deg, r)
    )


def polynomial_trilaplacian(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Trilaplacian via local polynomial fit. 3-D not yet supported."""
    if f.ndim == 3:
        raise NotImplementedError("Polynomial backend does not support 3-D fields yet")
    f = _apply_prefilter(f, Lx, Ly, prefilter_sigma)
    if 'get_default_stencil' in globals():
        deg, r = get_default_stencil(6)
    else:
        deg, r = _DEFAULTS[6]
    return (
        _deriv(f, Lx, Ly, 6, 0, deg, r)
        + 3 * _deriv(f, Lx, Ly, 4, 2, deg, r)
        + 3 * _deriv(f, Lx, Ly, 2, 4, deg, r)
        + _deriv(f, Lx, Ly, 0, 6, deg, r)
    )


def polynomial_third_order(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Third-order partial derivatives via local polynomial fit. 3-D not yet supported."""
    if f.ndim == 3:
        raise NotImplementedError("Polynomial backend does not support 3-D fields yet")
    f = _apply_prefilter(f, Lx, Ly, prefilter_sigma)
    if 'get_default_stencil' in globals():
        deg, r = get_default_stencil(3)
    else:
        deg, r = _DEFAULTS[3]
    fxxx = _deriv(f, Lx, Ly, 3, 0, deg, r)
    fxxy = _deriv(f, Lx, Ly, 2, 1, deg, r)
    fxyy = _deriv(f, Lx, Ly, 1, 2, deg, r)
    fyyy = _deriv(f, Lx, Ly, 0, 3, deg, r)
    return fxxx, fxxy, fxyy, fyyy


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


def poly_partial(f, Lx, Ly, a, b, deg=None, r=None):
    """Compute the mixed partial derivative ∂^(a+b) f / ∂x^a ∂y^b via local
    2-D polynomial fitting.

    This is the generic public interface to the polynomial differentiation
    backend.  It is used by build_pdefind_library(method='polynomial') to
    compute arbitrary multi-index derivatives, complementing the named
    functions (polynomial_grad, polynomial_laplacian, …) which are optimised
    for specific operators used by the equivariant library.

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
    if deg is None or r is None:
        _deg, _r = get_default_stencil(a + b)  # respects runtime overrides
        if deg is None:
            deg = _deg
        if r is None:
            r = _r
    return _deriv(f, Lx, Ly, a, b, deg, r)


def poly_partial_1d(u, L, d, deg=None, r=None):
    """Compute the d-th derivative of a 1-D periodic field via local polynomial fitting.

    Uses the same 2-D kernel approach restricted to 1-D: fits a degree-`deg`
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
    if deg is None or r is None:
        _deg, _r = _default_stencil(d)
        if deg is None:
            deg = _deg
        if r is None:
            r = _r

    N = len(u)
    dx = L / N
    offsets = list(range(-r, r + 1))

    # 1-D Vandermonde: V[pt, mono] = (di*dx)^p for p = 0, ..., deg
    V = np.array(
        [[(di * dx) ** p for p in range(deg + 1)] for di in offsets], dtype=float
    )  # (2r+1, deg+1)
    Vp = np.linalg.pinv(V)  # (deg+1, 2r+1)

    # Row d of Vp gives the LS coefficient c_d; derivative = d! * c_d.
    kernel = Vp[d] * factorial(d)  # shape (2r+1,)

    # Apply as a periodic 1-D correlation.
    result = np.zeros(N, dtype=float)
    for di, w in zip(offsets, kernel):
        result += w * np.roll(u, -di)
    return result
