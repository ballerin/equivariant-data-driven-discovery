"""
Derivative backend selector for the different differential operators.

Currently supports two backends, both N-D-generic (see spectral.py and
polynomial.py module docstrings for the isotropic-beyond-3D convention):
 - 'spectral' (FFT-based) [DEFAULT]
 - 'polynomial' (local least-squares polynomial fit).

Potential todos:
 - finite-difference
 - stencil-based
 - Tikhonov-regularized polynomial fit

Exposes six backend-neutral named functions (grad, laplacian, hessian,
bilaplacian, trilaplacian, third_order) plus a generic ``partial_nd`` for an
arbitrary derivative multi-index; all dispatch to whichever backend is
currently active.

Usage
-----
Default backend is 'spectral' (FFT, exact for band-limited data).

    import derivatives
    derivatives.set_backend('polynomial')   # switch globally

    # or, temporarily inside a block:
    with derivatives.use_backend('polynomial'):
        ops = get_operators(2, 3)           # evaluated with polynomial kernels

After a backend switch, all subsequent operator evaluations use the new backend;
calls go through the module-level functions below, so the switch is immediate.

Available backends
------------------
'spectral':   
    - derivatives/spectral.py
    - FFT-based, exact for band-limited data
    - [DEFAULT]

'polynomial':
    - derivatives/polynomial.py
    - local least-squares polynomial fit
    - robust to measurement noise at the cost of a small bias
"""

from contextlib import contextmanager

from . import polynomial as _polynomial
from . import spectral as _spectral
from ._common import axis_lengths as _axis_lengths

# Dispatch table: backend name → dict of operation → function.
# Each entry maps the neutral operation name to the backend-specific function.
_DISPATCH = {
    "spectral": {
        "grad": _spectral.spectral_grad,
        "laplacian": _spectral.spectral_laplacian,
        "hessian": _spectral.spectral_hessian,
        "bilaplacian": _spectral.spectral_bilaplacian,
        "trilaplacian": _spectral.spectral_trilaplacian,
        "third_order": _spectral.spectral_third_order,
    },
    "polynomial": {
        "grad": _polynomial.polynomial_grad,
        "laplacian": _polynomial.polynomial_laplacian,
        "hessian": _polynomial.polynomial_hessian,
        "bilaplacian": _polynomial.polynomial_bilaplacian,
        "trilaplacian": _polynomial.polynomial_trilaplacian,
        "third_order": _polynomial.polynomial_third_order,
    },
}

# Neutral multi-index primitive per backend: fn(f, Ls, multi_index) -> ndarray.
_PARTIAL_ND_DISPATCH = {
    "spectral": _spectral.spectral_partial_nd,
    "polynomial": _polynomial.poly_partial_nd,
}

_active_name = "spectral"  # default


def set_backend(name: str) -> None:
    """Switch the active derivative backend globally.

    Parameters
    ----------
    name : 'spectral' or 'polynomial'
    """
    global _active_name
    if name not in _DISPATCH:
        raise ValueError(f"Unknown backend {name!r}. Choose from: {list(_DISPATCH)}")
    _active_name = name


def get_backend() -> str:
    """Return the name of the currently active backend."""
    return _active_name


@contextmanager
def use_backend(name: str):
    """Context manager: temporarily switch backend, then restore the previous one.

    Example
    -------
    with derivatives.use_backend('polynomial'):
        Theta, descs = build_multigraph_library(f, Lx, Ly)
    """
    previous = get_backend()
    set_backend(name)
    try:
        yield
    finally:
        set_backend(previous)


# ---------------------------------------------------------------------------
# Public derivative APIs (backend-neutral names)
#
# Parameters
# ----------
# f : 2D array - function values on a uniform grid
# Lx, Ly : float - physical domain sizes in x and y
# prefilter_sigma - apply Gaussian prefilter by sigma (if provided)
# ---------------------------------------------------------------------------


def grad(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Gradient via the active backend, for any ndim.

    Returns a tuple of f.ndim components (f_x, f_y, ...).
    """
    return _DISPATCH[_active_name]["grad"](f, Lx, Ly, Lz, prefilter_sigma)


# Public re-export of the single canonical implementation (derivatives/
# _common.py), which both backends also use.
axis_lengths = _axis_lengths


def partial_nd(f, Ls, multi_index):
    """Compute ∂^multi_index f via the active backend, for any ndim.

    Parameters
    ----------
    f : ndarray, shape (N_1, ..., N_ndim)
    Ls : sequence of float, length f.ndim — per-axis domain length
    multi_index : sequence of int, length f.ndim

    Returns
    -------
    ndarray, shape f.shape
    """
    return _PARTIAL_ND_DISPATCH[_active_name](f, Ls, multi_index)


def laplacian(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Laplacian via the active backend, for any ndim: Δf = Σ_axis f_{axis,axis}."""
    return _DISPATCH[_active_name]["laplacian"](f, Lx, Ly, Lz, prefilter_sigma)


def hessian(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Unique Hessian entries via the active backend, for any ndim.

    2-D: returns (f_xx, f_xy, f_yy).  3-D: returns 6 components.  N-D:
    returns the C(ndim+1, 2) unique components; see spectral.py's module
    docstring for the canonical ordering.
    """
    return _DISPATCH[_active_name]["hessian"](f, Lx, Ly, Lz, prefilter_sigma)


def bilaplacian(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Bilaplacian Δ²f via the active backend, for any ndim."""
    return _DISPATCH[_active_name]["bilaplacian"](f, Lx, Ly, Lz, prefilter_sigma)


def trilaplacian(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Trilaplacian Δ³f via the active backend, for any ndim."""
    return _DISPATCH[_active_name]["trilaplacian"](f, Lx, Ly, Lz, prefilter_sigma)


def third_order(f, Lx, Ly, Lz=None, prefilter_sigma=None):
    """Unique third-order partial derivatives via the active backend, for any ndim.

    2-D: returns (f_xxx, f_xxy, f_xyy, f_yyy).  3-D: returns 10 components.
    N-D: returns the C(ndim+2, 3) unique components; see spectral.py's
    module docstring for the canonical ordering.
    """
    return _DISPATCH[_active_name]["third_order"](f, Lx, Ly, Lz, prefilter_sigma)


# -------------------- configuration helpers --------------------
def set_polynomial_stencil(deg=None, r=None, order=None):
    """Set the default polynomial stencil used by the polynomial backend.

    Parameters
    ----------
    deg, r : int or None
        If both provided, set the override to (deg, r). If either is None,
        the override for the target `order` (or global) is removed.
    order : int or None
        If provided, set the override only for that total derivative order
        (a+b). If None, set a global override used for all orders.
    """
    # delegate to polynomial backend module
    _polynomial.set_default_stencil(deg=deg, r=r, order=order)


def get_polynomial_stencil(order):
    """Return the (deg, r) pair currently used for ``order`` by the polynomial backend."""
    return _polynomial.get_default_stencil(order)
