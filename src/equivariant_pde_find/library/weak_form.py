"""
Weak-form utilities to evlaluate library terms via convolution with
compactly-supported test functions.

Instead of evaluating the library at individual grid points (pointwise 
formulation), each operator is integrated against compactly-supported
test functions ψ:

    G_{ij,k} = ∫ ψ_{ij}(x,y) · N_k[f](x,y)  dx dy   # Space
    b_{ij}   = ∫ ψ_{ij}(x,y) · f_t(x,y)     dx dy   # Time

where ψ_{ij} is centred at query point (i, j).  The convolution provides a
natural noise averaging due to integrating over a support of width 2m+1 in
both directions.

For the equivariant library, operators like ‖∇f‖² and Δf cannot all be
rewritten as ∂^α(f^β), so we use a **direct** approach:
  - Compute each operator N_k[f] on the full grid (spectral or polynomial
    backend) via the existing ``multigraph_operators`` functions,
  - Convolve with separable polynomial-bump test functions,
  - Subsample at query points to build (G, b).

Regression is then ran by feeding G and b which replace the pointwise Θ and f_t.
"""

from contextlib import contextmanager
from functools import lru_cache
import math

import numpy as np


# ---------------------------------------------------------------------------
# Spatial boundary handling
#
# 'valid' (default, Messenger & Bortz): a test function is only centred where
# its whole support fits inside the grid, so each spatial axis loses 2m
# points.  'periodic': the support wraps around the domain, so every grid
# point is a query point.  Exact for the periodic domains simulated here (the
# integration-by-parts boundary terms cancel on a torus), and the only way to
# keep query points when the test-function rule asks for m close to N/2.
# Applies to every spatial weak convolution (all library builders and
# build_weak_rhs) through weak_convolution; the time axis of a trajectory is
# never periodic and is not affected.
#
# 'periodic' is the setting used for the experiments from 2026-09 on, for every
# equation.  All simulated domains are periodic; the derivative backends
# (spectral, and the polynomial backend's periodic stencils) already assume
# it, so 'valid' would not make the pipeline periodicity-free; KS-cell and
# KS-wchaos need it (the rule asks m ~ N/2); and it keeps one protocol across
# equations and grid sizes.  The consequence is that the method as run assumes
# periodic data, which results must state.  'valid' stays the default so that
# earlier sweeps reproduce unchanged.
# ---------------------------------------------------------------------------

WEAK_BOUNDARIES = ("valid", "periodic")
_weak_boundary = "valid"


def get_weak_boundary():
    """Name of the active spatial boundary mode ('valid' or 'periodic')."""
    return _weak_boundary


def set_weak_boundary(boundary):
    """Switch the spatial boundary mode globally."""
    global _weak_boundary
    if boundary not in WEAK_BOUNDARIES:
        raise ValueError(f"Unknown weak boundary {boundary!r}; expected one "
                         f"of {WEAK_BOUNDARIES}.")
    _weak_boundary = boundary


@contextmanager
def use_weak_boundary(boundary):
    """Context manager: temporarily switch the spatial boundary mode, then
    restore the previous one (``None`` leaves it unchanged).

    Example
    -------
    with weak_form.use_weak_boundary("periodic"):
        G, desc = build_multigraph_library(f, Lx, Ly, formulation="space_weak",
                                           m=24, p=10)   # G.shape[0] == f.size
    """
    previous = get_weak_boundary()
    if boundary is not None:
        set_weak_boundary(boundary)
    try:
        yield
    finally:
        set_weak_boundary(previous)


def weak_output_shape(shape, m, boundary=None):
    """Spatial shape of a weak column for half-support *m*: ``N − 2m`` per
    axis for 'valid', ``N`` for 'periodic' (``boundary=None``: active mode)."""
    boundary = get_weak_boundary() if boundary is None else boundary
    if boundary == "periodic":
        return tuple(shape)
    return tuple(n - 2 * m for n in shape)


# ---------------------------------------------------------------------------
# Test function construction
# ---------------------------------------------------------------------------


def _test_function_1d(m, p):
    """Build a 1-D polynomial bump test function supported on [-m, m] in
    grid-point units (dx=1).

    ψ(x) = (1 - (x / m)²)^p   for |x| < m,  0 otherwise.

    The function is C^{p-1} at the boundary and integrates to a known
    normalisation constant.

    Parameters
    ----------
    m : int - size of half-support (in grid points)
    p : int - polynomial degree

    Returns
    -------
    kernel : ndarray of shape (2*m + 1,) - values of ψ at grid points
    """
    if m < 2 or p < 1:
        raise ValueError("Require m >= 2, p >= 1.")
    n_pts = np.arange(-m, m + 1, dtype=float)
    x = n_pts / m  # scaled to [-1, 1]
    phi = np.maximum(0.0, 1.0 - x ** 2) ** p
    return phi.astype(np.float64)


def _test_function_deriv_1d(m, p, order):
    """Return the *order*-th derivative of φ(t) = (1−t²)^p at grid points.

    Computed analytically via the polynomial expansion of (1−t²)^p and
    scaled by 1/m^order to account for d/dx = (1/m) · d/dt in grid-point
    units (dx=1).  The result is exact.

    For meaningful IBP application ``order < p`` is required, otherwise
    a discontinuity arises at the boundary.  ``order == p`` is borderline
    (the polynomial still exists but has a jump at t = ±1).

    Parameters
    ----------
    m : int - size of half-support (in grid points)
    p : int - polynomial degree (must be > order for a valid IBP kernel)
    order : int (≥ 0; 0 returns φ itself)

    Returns
    -------
    kernel : ndarray of shape (2*m + 1,)

    Raises
    ------
    ValueError
        If ``order > p`` (the polynomial collapses; derivative is
        identically zero or otherwise invalid).
    """
    if order > p:
        raise ValueError(
            f"order={order} must be ≤ p={p}"
        )
    if order == 0:
        return _test_function_1d(m, p)
    # Build polynomial coefficients of (1 − t²)^p = Σ C(p,k)(−1)^k t^{2k}
    coeffs = np.zeros(2 * p + 1)
    for k in range(p + 1):
        coeffs[2 * k] = math.comb(p, k) * ((-1) ** k)
    poly = np.polynomial.Polynomial(coeffs)
    for _ in range(order):
        poly = poly.deriv()
    t = np.arange(-m, m + 1, dtype=float) / m
    result = np.where(np.abs(t) < 1.0, poly(t), 0.0)
    return result.astype(np.float64) / (m ** order)


# ---------------------------------------------------------------------------
# IBP (integration-by-parts) kernels for linear differential operators
# ---------------------------------------------------------------------------


@lru_cache(maxsize=64)
def _ibp_laplacian_kernels(m, p):
    """Return (kernel_xx, kernel_0) for IBP of the Laplacian.

    Uses the **analytic** second derivative of the test function, computed
    from the exact polynomial expansion of φ(t) = (1−t²)^p.  This matches
    the spectral Laplacian more closely than the discrete 3-point stencil
    (relative error ~0.05 % vs ~0.86 % for p=4, m=8 on N=128 grids).

    The kernel is returned in grid-point units (dx=1); the caller scales
    by 1/dx² for physical units.

    Returns
    -------
    kernel_xx : ndarray (2*m+1,) — analytic φ'' (grid-point units)
    kernel_0  : ndarray (2*m+1,) — φ
    """
    phi = _test_function_1d(m, p)
    d2 = _test_function_deriv_1d(m, p, 2)
    return d2, phi


@lru_cache(maxsize=64)
def _ibp_bilaplacian_kernels(m, p):
    """Return (kernel_xxxx, kernel_xx, kernel_0) for IBP of the bilaplacian.

    Uses analytic 4th derivative of φ via ``_test_function_deriv_1d`` when
    p ≥ 5 (φ ∈ C⁴, 4th derivative continuous at boundary).  Falls back to
    discrete D4 = D2∘D2 when p ≤ 4 (φ ∈ C³, analytic φ⁗ has a jump).

    With separable ψ=φ(x)·φ(y):
        Δ²ψ = φ⁗(x)·φ(y)/dx⁴ + 2·φ″(x)·φ″(y)/(dx²·dy²) + φ(x)·φ⁗(y)/dy⁴

    Returns
    -------
    kernel_xxxx : ndarray (2*m+1,) — φ⁗ (analytic, p≥5) or D4_φ (discrete)
    kernel_xx   : ndarray (2*m+1,) — φ″ (analytic) or D2_φ (discrete, p≤4)
    kernel_0    : ndarray (2*m+1,) — φ
    """
    phi = _test_function_1d(m, p)
    if p >= 5:
        # φ ∈ C⁴ — analytic derivatives are accurate
        d2 = _test_function_deriv_1d(m, p, 2)
        d4 = _test_function_deriv_1d(m, p, 4)
    else:
        # φ ∈ C³ — fall back to discrete differences
        d2 = np.roll(phi, 1) + np.roll(phi, -1) - 2.0 * phi
        d4 = np.roll(d2, 1) + np.roll(d2, -1) - 2.0 * d2
    return d4, d2, phi


# ---------------------------------------------------------------------------
# IBP weak-form evaluation for linear operators
# ---------------------------------------------------------------------------


def weak_laplacian_ibp(f, m, p, Lx, Ly, Lz=None):
    """Compute ∫ ψ · Δf ≈ ∫ Δψ · f via integration by parts, for any ndim.

    No derivatives of *f* are computed — the differential operator is
    moved entirely onto the analytically-known test function.

    Δψ = Σ_axis φ″(axis)/dx_axis² · ∏_{other axes} φ(other), so each axis
    contributes one term with that axis's kernel swapped from φ to φ″.

    Parameters
    ----------
    f : ndarray, shape (N_1, ..., N_ndim)
    m, p : int — test function parameters
    Lx, Ly, Lz : float — domain lengths (isotropic beyond ndim=3; see
        derivatives.axis_lengths)

    Returns
    -------
    ndarray, shape (N_1 - 2m, ..., N_ndim - 2m)
    """
    from .. import derivatives

    kxx_grid, k0 = _ibp_laplacian_kernels(m, p)
    ndim = f.ndim
    Ls = derivatives.axis_lengths(f, Lx, Ly, Lz)
    dx2s = [(Ls[axis] / f.shape[axis]) ** 2 for axis in range(ndim)]

    total = None
    for axis in range(ndim):
        kernel_list = [k0] * ndim
        kernel_list[axis] = kxx_grid / dx2s[axis]
        term = weak_convolution(f, kernel_list)
        total = term if total is None else total + term
    return total


def weak_bilaplacian_ibp(f, m, p, Lx, Ly, Lz=None):
    """Compute ∫ ψ · Δ²f ≈ ∫ Δ²ψ · f via integration by parts, for any ndim.

    No derivatives of *f* are computed.

    Δ²ψ = Σ_axis φ⁗(axis)/dx_axis⁴ · ∏_{other}φ(other)
          + 2 Σ_{i<j} φ″(i)φ″(j)/(dx_i²dx_j²) · ∏_{other}φ(other)
    i.e. ndim pure 4th-derivative terms plus C(ndim,2) mixed terms with
    coefficient 2 (2-D: 2 pure + 1 mixed; 3-D: 3 pure + 3 mixed).

    Parameters
    ----------
    f : ndarray, shape (N_1, ..., N_ndim)
    m, p : int — test function parameters
    Lx, Ly, Lz : float — domain lengths (isotropic beyond ndim=3)

    Returns
    -------
    ndarray, shape (N_1 - 2m, ..., N_ndim - 2m)
    """
    from .. import derivatives

    kxxxx_grid, kxx_grid, k0 = _ibp_bilaplacian_kernels(m, p)
    ndim = f.ndim
    Ls = derivatives.axis_lengths(f, Lx, Ly, Lz)
    dx2s = [(Ls[axis] / f.shape[axis]) ** 2 for axis in range(ndim)]

    total = None
    for i in range(ndim):
        kernel_list = [k0] * ndim
        kernel_list[i] = kxxxx_grid / (dx2s[i] ** 2)
        term = weak_convolution(f, kernel_list)
        total = term if total is None else total + term
    for i in range(ndim):
        for j in range(i + 1, ndim):
            kernel_list = [k0] * ndim
            kernel_list[i] = kxx_grid / dx2s[i]
            kernel_list[j] = kxx_grid / dx2s[j]
            total = total + 2.0 * weak_convolution(f, kernel_list)
    return total


# ---------------------------------------------------------------------------
# Weak convolution
# ---------------------------------------------------------------------------


def weak_convolution(f, *kernels):
    """Compute the weak form of field `f` by separable 1-D convolution, for
    any ndim.

    Two calling conventions are supported:

        weak_convolution(f, kernel_0, kernel_1, ...)   # positional, historical
        weak_convolution(f, [kernel_0, kernel_1, ...])  # single list, N-D

    ``None`` entries are dropped (matches the historical
    ``weak_convolution(f, kx, ky, None)`` 2-D calling pattern from callers
    that always pass a 3rd, possibly-None argument).

    Parameters
    ----------
    f : ndarray, shape (N_1, ..., N_ndim)
    kernels : one 1-D ndarray per axis (or a single list of them)

    Returns
    -------
    ndarray, shape (N_1 - 2*m_1, ..., N_ndim - 2*m_ndim) under the default
    'valid' boundary; shape (N_1, ..., N_ndim), entry i centred at grid point
    i, under 'periodic' (see :func:`use_weak_boundary`).
    """
    if len(kernels) == 1 and isinstance(kernels[0], (list, tuple)):
        kernel_list = list(kernels[0])
    else:
        kernel_list = list(kernels)
    kernel_list = [k for k in kernel_list if k is not None]
    if len(kernel_list) != f.ndim:
        raise ValueError(
            f"weak_convolution: got {len(kernel_list)} kernels for a "
            f"{f.ndim}-D field; need exactly f.ndim kernels."
        )

    if get_weak_boundary() == "periodic":
        from scipy.ndimage import convolve1d

        # convolve1d flips the kernel like np.convolve and centres an
        # odd-length kernel, so entry i equals the 'valid' entry centred at i.
        result = np.asarray(f, dtype=float)
        for axis, kernel in enumerate(kernel_list):
            kernel = np.ravel(kernel)
            if kernel.size > result.shape[axis]:
                raise ValueError(
                    f"weak_convolution: a {kernel.size}-point kernel does not "
                    f"fit a periodic axis of {result.shape[axis]} points."
                )
            result = convolve1d(result, kernel, axis=axis, mode="wrap")
        return result

    result = f
    for axis, kernel in enumerate(kernel_list):
        kernel = np.ravel(kernel)
        result = np.apply_along_axis(
            lambda v, _k=kernel: np.convolve(v, _k, mode="valid"), axis, result
        )
    return result


# ---------------------------------------------------------------------------
# Cache kernel-building (same (m, p, dx) → same kernel)
# ---------------------------------------------------------------------------

@lru_cache(maxsize=64)
def _cached_kernel_1d(m, p):
    return _test_function_1d(m, p)


# ---------------------------------------------------------------------------
# Weak library construction (single snapshot)
# ---------------------------------------------------------------------------


def _is_laplacian_of(desc):
    """Determine whether a (modified) multigraph operator is of the form
    Δ(something) or Δ²(something), and if so, return the appropriate
    IBP tag and the callable that computes the *inner* argument.

    Returns
    -------
    (operator_type, inner_func)  or  None
        operator_type is ``'laplacian'`` or ``'bilaplacian'``.
        ``inner_func(f, Lx, Ly)`` computes the argument of Δ / Δ².
    """
    # Laplacian-of operators: Δf, Δf², Δf³, Δ‖∇f‖², Δ(f·‖∇f‖²)
    if desc == "Δf":
        return ("laplacian", lambda f, Lx, Ly, Lz=None: f)
    if desc == "Δf²":
        return ("laplacian", lambda f, Lx, Ly, Lz=None: f**2)
    if desc == "Δf³":
        return ("laplacian", lambda f, Lx, Ly, Lz=None: f**3)
    if desc == "Δ‖∇f‖²":
        from ..multigraph_operators import _op_grad_sq
        return ("laplacian", _op_grad_sq)
    if desc == "Δ(f·‖∇f‖²)":
        from ..multigraph_operators import _op_grad_sq
        return ("laplacian", lambda f, Lx, Ly, Lz=None: f * _op_grad_sq(f, Lx, Ly, Lz))

    # Bilaplacian-of operators: Δ²f, Δ²f², Δ²f³
    if desc == "Δ²f":
        return ("bilaplacian", lambda f, Lx, Ly, Lz=None: f)
    if desc == "Δ²f²":
        return ("bilaplacian", lambda f, Lx, Ly, Lz=None: f**2)
    if desc == "Δ²f³":
        return ("bilaplacian", lambda f, Lx, Ly, Lz=None: f**3)

    return None


def build_weak_rhs(ft, kernel_x=None, kernel_y=None, kernel_z=None, m=8, p=4,
                   kernels=None):
    """Build the weak-form RHS vector b = ψ ∗ f_t, for any ndim.

    Parameters
    ----------
    ft : ndarray, shape (N_1, ..., N_ndim) — time derivative of f.
    kernel_x, kernel_y, kernel_z : ndarray or None
        Per-axis overrides for axes 0, 1, 2 (historical 2-D/3-D calling
        convention).  Any axis without an explicit override (including all
        axes beyond the third, for ndim>3) is built from (m, p).
    kernels : list of ndarray or None
        If given, used directly as the full per-axis kernel list (must have
        length ft.ndim); takes precedence over kernel_x/kernel_y/kernel_z.
    m, p : int
        Used to build the default kernel for axes without an explicit
        override.

    Returns
    -------
    b : ndarray of shape (N_valid,)
    """
    if kernels is not None:
        kernel_list = list(kernels)
        if len(kernel_list) != ft.ndim:
            raise ValueError(
                f"build_weak_rhs: kernels has length {len(kernel_list)}, "
                f"expected {ft.ndim} (ft.ndim)"
            )
        return weak_convolution(ft, kernel_list).ravel()

    overrides = [kernel_x, kernel_y, kernel_z]
    default = _cached_kernel_1d(m, p)
    kernel_list = [
        overrides[axis] if axis < len(overrides) and overrides[axis] is not None
        else default
        for axis in range(ft.ndim)
    ]
    return weak_convolution(ft, kernel_list).ravel()


# ---------------------------------------------------------------------------
# MSTLS dominant-balance scale factors
# ---------------------------------------------------------------------------


def compute_m_scale(G, b):
    """Compute per-column scale factors for the MSTLS dominant-balance bounds.

    Following the WSINDy pipeline (``get_scales.m`` → ``wsindy_pde_RGLS_seq``),
    the scale vector *M* captures the relative expected magnitude of each
    library column compared to the RHS::

        M_j = ‖G_j‖ / ‖b‖

    where *G_j* is the *j*-th column of the weak library and *b* is the weak
    RHS.  This is a first-order approximation to the original ``M_scale``
    (which depends on test function parameters and library structure via
    ``get_scales.m``).  When passed to ``mstls``, it enables per-column
    dominant-balance bounds ``[λ·max(1,bnds_j), (1/λ)·min(1,bnds_j)]``.

    Parameters
    ----------
    G : ndarray, shape (n, d) — weak library matrix
    b : ndarray, shape (n,)   — weak RHS vector

    Returns
    -------
    M : ndarray, shape (d,)
    """
    col_norms = np.linalg.norm(G, axis=0)
    b_norm = np.linalg.norm(b)
    return col_norms / max(b_norm, 1e-30)


# ---------------------------------------------------------------------------
# Automatic test-function parameter selection (changepoint method)
#
# Messenger & Bortz (2021), Appendix A: "Learning Test Functions From Data"
# ---------------------------------------------------------------------------


def estimate_critical_wavenumber(U, axis=-1):
    """Estimate the wavenumber separating signal- from noise-dominated modes.

    Implements the changepoint method from Appendix A of the WSINDy PDE paper:
    the cumulative sum of the time-averaged power spectrum is approximately
    curved over signal-dominated modes and linear over noise-dominated modes.
    The changepoint ``k*`` is found as the best 2-segment piecewise-linear fit.

    Parameters
    ----------
    U : ndarray, shape (..., N_axis)
        Spatiotemporal data.  ``axis`` identifies the dimension along which
        to estimate the critical wavenumber (e.g. ``axis=0`` for x-direction
        of an (Nx, Ny, Nt) array).
    axis : int
        Axis along which to compute the FFT.

    Returns
    -------
    k_star : int
        The critical positive wavenumber index (0 = DC, N//2 = Nyquist).
    """
    N = U.shape[axis]
    # Compute FFT along chosen axis, then average over all other axes.
    fft_U = np.fft.fft(U, axis=axis)
    power = np.abs(fft_U)
    axes_to_avg = tuple(i for i in range(U.ndim) if i != axis)
    avg_power = np.mean(power, axis=axes_to_avg)  # shape (N,)
    # fftshift so that negative freqs are first, DC in middle, positive last
    avg_power_s = np.fft.fftshift(avg_power)
    cumsum = np.cumsum(avg_power_s)

    # Search for best 2-segment piecewise-linear fit in the central region.
    # Avoid the flat DC tail at the extremes.
    lo, hi = max(N // 8, 5), min(7 * N // 8, N - 5)
    best_resid, best_k = np.inf, N // 2
    x_all = np.arange(N, dtype=float)
    for k in range(lo, hi):
        # Left segment [0, k), right segment [k, N)
        x1, y1 = x_all[:k], cumsum[:k]
        x2, y2 = x_all[k:], cumsum[k:]

        # Fit line: y = m*x + c via normal equations
        A1 = np.column_stack([x1, np.ones_like(x1)])
        m1_c1 = np.linalg.lstsq(A1, y1, rcond=None)[0]
        resid1 = np.sum((y1 - (m1_c1[0] * x1 + m1_c1[1])) ** 2)

        A2 = np.column_stack([x2, np.ones_like(x2)])
        m2_c2 = np.linalg.lstsq(A2, y2, rcond=None)[0]
        resid2 = np.sum((y2 - (m2_c2[0] * x2 + m2_c2[1])) ** 2)

        total = resid1 + resid2
        if total < best_resid:
            best_resid, best_k = total, k

    # Convert changepoint index to positive wavenumber (0 = DC)
    k_star = best_k - N // 2
    return max(0, abs(int(k_star)))


def _compute_p_from_m(m, alpha_bar, tau=1e-10):
    """Compute test-function degree *p* from support *m* using the
    real-space decay condition from Algorithm 1 of the paper:

        p = min{ p ≥ ᾱ + 1  :  (1 − (1 − 1/m)²)^p ≤ τ }.

    Parameters
    ----------
    m : int
    alpha_bar : int — max derivative order along this axis
    tau : float — decay tolerance

    Returns
    -------
    p : int
    """
    base = 1.0 - (1.0 - 1.0 / m) ** 2  # = (2m-1)/m²
    if base <= 0:
        return max(alpha_bar + 1, 10)
    # p_min from inequality: p ≥ log(τ) / log(base)
    p_min = int(math.ceil(math.log(tau) / math.log(base)))
    return max(alpha_bar + 1, p_min)


def select_test_function_degree(m, max_deriv=4, tau=1e-10):
    """Public wrapper over :func:`_compute_p_from_m`: the degree *p* that
    Messenger & Bortz (2021), Algorithm 1, prescribe for a support of *m*
    grid points and derivatives up to order *max_deriv*.

    The MATLAB reference computes the same quantity in
    ``WSINDy_PDE/phi_int_weights.m``::

        p = ceil(max(log(tol)/log((2*m-1)/m^2), d+1));

    Using a fixed *p* instead (the sweeps previously hardcoded ``p = 6``)
    leaves φ⁗ large at the edge of the support, and the grid sum standing in
    for ∫ψ⁗f then carries a large quadrature error — 3 % at m=12, 84 % at
    m=8, in the ∂⁴ columns that dominate Cahn-Hilliard, Swift-Hohenberg and
    Kuramoto-Sivashinsky.

    Always check :func:`test_function_is_resolved` before using the result:
    for small *m* the prescribed *p* is large, and the bump becomes too
    narrow for the grid to represent at all.
    """
    return _compute_p_from_m(m, max_deriv, tau)


def test_function_is_resolved(m, p, max_deriv=4):
    """Whether a (m, p) test function can actually be represented on the grid.

    φ(t) = (1 − t²)^p concentrates as p grows: its core has half-width
    ≈ 1/√p in t, i.e. ≈ m/√p grid points.  Resolving the max_deriv-th
    derivative needs several points across that core, so we require

        m >= 2 * sqrt(p)

    (about 4 points across the core).  Measured relative error of the ∂⁴
    weak-IBP column against a refined-grid reference, at the p prescribed by
    :func:`select_test_function_degree`:

        m=4  (p=28)  3.8e+04   <- fails the criterion, unusable
        m=6  (p=20)  9.9e-02   <- fails the criterion
        m=8  (p=16)  4.9e-04   <- passes, marginally
        m=12 (p=13)  5.1e-06   <- passes comfortably

    A support that fails this cannot support IBP of that derivative order at
    *any* degree: at m=4 the fixed p=6 gives 2.4e+02 and the prescribed p=28
    gives 3.8e+04.  Nine grid points cannot carry a fourth derivative.

    Parameters
    ----------
    m : int — test-function half-support, in grid points
    p : int — test-function polynomial degree
    max_deriv : int — highest derivative order moved onto the test function

    Returns
    -------
    bool
    """
    if max_deriv <= 0:
        return True
    return m >= 2.0 * math.sqrt(p)


def auto_test_function_params(
    U, dx, dt=1.0, max_deriv_spatial=4, max_deriv_temporal=1,
    tau=1e-10, tau_hat=2.0,
):
    """Automatically select (m_x, p_x, m_t, p_t) from spatiotemporal data.

    Uses the changepoint method (Appendix A) to estimate critical wavenumbers
    k*_x, k*_t and then solves for (m, p) via the simultaneous decay
    constraints in real and Fourier space.

    Parameters
    ----------
    U : ndarray
        Spatiotemporal data of shape ``(Nx, Ny)`` or ``(Nx, Ny, Nt)``.
    dx : float
        Grid spacing (same for x and y).
    dt : float
        Temporal grid spacing (used only if U is 3-D).  Default 1.0 for
        single-snapshot 2-D data.
    max_deriv_spatial : int
        Maximum spatial derivative order in the library (e.g. 4 for Δ²).
    max_deriv_temporal : int
        Maximum temporal derivative order (usually 1 for ∂_t).
    tau : float
        Real-space decay tolerance (default 1e-10 as in the paper).
    tau_hat : float
        Number of standard deviations into the tail of φ̂ for the corner
        wavenumber (paper uses 1–3).

    Returns
    -------
    m_x : int
    p_x : int
    m_t : int
    p_t : int
    """
    Nx = U.shape[0]
    has_time = (U.ndim == 3)

    # --- spatial: x-axis ---
    k_star_x = estimate_critical_wavenumber(U, axis=0)
    m_x = _solve_m_from_changepoint(k_star_x, Nx, dx, tau, tau_hat)
    p_x = _compute_p_from_m(m_x, max_deriv_spatial, tau)

    # --- temporal (if available) ---
    if has_time:
        Nt = U.shape[2]
        k_star_t = estimate_critical_wavenumber(U, axis=2)
        m_t = _solve_m_from_changepoint(k_star_t, Nt, dt, tau, tau_hat)
    else:
        # Fallback for single snapshot: use spatial m, scaled to temporal
        m_t = max(2, m_x // 2)
    p_t = _compute_p_from_m(m_t, max_deriv_temporal, tau)

    return m_x, p_x, m_t, p_t


def _solve_m_from_changepoint(k_star, N, dx, tau=1e-10, tau_hat=2.0):
    """Solve for support size *m* from the critical wavenumber k*.

    Finds the root of

        F(m) = log((2m−1)/m²)·(4π²k*²m² − 3N²τ̂²) − 2N²τ̂²·log(τ)

    by scanning for a sign change in the feasible interval.

    Parameters
    ----------
    k_star : int
        Critical wavenumber (positive, from ``estimate_critical_wavenumber``).
    N : int
        Number of grid points along this axis.
    dx : float
        Grid spacing (not used in the root — the formula uses grid-point units).
    tau : float
    tau_hat : float

    Returns
    -------
    m : int
    """
    _ = dx  # not used in the root equation (grid-point units)

    def F(m):
        if m <= 1:
            return -np.inf
        term1 = math.log((2.0 * m - 1.0) / (m ** 2))
        k2 = (2.0 * math.pi * k_star) ** 2
        term2 = k2 * m ** 2 - 3.0 * (N ** 2) * (tau_hat ** 2)
        return term1 * term2 - 2.0 * (N ** 2) * (tau_hat ** 2) * math.log(tau)

    # Valid search interval from the paper (after eq. in Appendix A)
    if k_star == 0:
        return max(2, N // 8)

    m_lo = max(2, int(math.sqrt(3) * N * tau_hat / (math.pi * k_star)))
    m_hi = max(
        m_lo + 1,
        int(m_lo * math.sqrt(max(1.0, 1.0 - 8.0 / math.sqrt(3) * math.log(tau)))),
    )
    # Guard against unreasonable values
    m_hi = min(m_hi, N // 2 - 1)

    # Scan for sign change; fall back to m_lo if none found.
    f_prev = F(m_lo)
    for m in range(m_lo + 1, m_hi + 1):
        f_cur = F(m)
        if f_prev * f_cur <= 0:
            return m if abs(f_cur) < abs(f_prev) else m - 1
        f_prev = f_cur

    return max(2, m_lo)


# ---------------------------------------------------------------------------
# Messenger & Bortz (2021) test-function rule, per axis, space and time
#
# arXiv:2007.02848, eq. (4.3) and Appendix A; reference implementation
# WSINDy_PDE/findcorners.m.  The functions above (estimate_critical_wavenumber,
# _solve_m_from_changepoint, auto_test_function_params) predate this section
# and are kept unchanged for reproducibility of earlier sweeps: they search the
# full symmetric spectrum, use twice the paper's lower root bracket, and assume
# a time-last (Nx, Ny, Nt) layout.  The functions below follow the paper and
# take the axis layout explicitly.
# ---------------------------------------------------------------------------


def mb_critical_wavenumber(U, axis):
    """Critical wavenumber k* along ``axis`` by the changepoint method.

    Appendix A, step 1: average |F^d(U)| over all other axes, take the
    cumulative sum H over the negative half-spectrum, and locate the corner of
    the best two-piece linear approximation to H, measured in relative ℓ²
    (weights 1/|H|).  Follows ``findcornerpts`` in ``findcorners.m``: each
    piece is the line through its own range, anchored at the outer end, and
    the returned wavenumber is the corner's distance from DC minus 3 (the
    reference implementation's conservative offset), floored at 1.

    The time axis of a trajectory is not periodic, so the FFT along it
    carries leakage; this is exactly what the reference code does too.

    Parameters
    ----------
    U : ndarray - data with the coordinate of interest on ``axis``
    axis : int

    Returns
    -------
    k_star : int (>= 1)
    """
    N = U.shape[axis]
    spec = np.abs(np.fft.fftshift(np.fft.fft(U, axis=axis), axes=axis))
    other = tuple(i for i in range(U.ndim) if i != axis % U.ndim)
    spec = spec.mean(axis=other) if other else spec
    H = np.cumsum(spec)[: int(math.ceil(N / 2))]
    H = np.maximum(H, np.finfo(float).tiny)
    NN = H.size
    x = np.arange(NN, dtype=float)
    if NN < 3:
        return 1

    errs = np.empty(NN - 2)
    for k in range(1, NN - 1):           # 0-based corner index
        H1, x1 = H[: k + 1], x[: k + 1]
        H2, x2 = H[k:], x[k:]
        s1 = np.ptp(H1) / np.ptp(x1)
        s2 = np.ptp(H2) / np.ptp(x2)
        L1 = H1.min() + s1 * (x1 - x1[0])
        L2 = H2.max() + s2 * (x2 - x2[-1])
        errs[k - 1] = math.sqrt(np.sum(((L1 - H1) / H1) ** 2)
                                + np.sum(((L2 - H2) / H2) ** 2))
    t_star = int(np.argmin(errs)) + 1    # 1-based index of the minimiser
    return max(NN - t_star - 3, 1)


def mb_root_m(k_star, N, tau=1e-10, tau_hat=2.0):
    """Real root m of Appendix A's F(m), capped at N/2 − 1.

        F(m) = log((2m−1)/m²)·(4π²k*²m² − 3N²τ̂²) − 2N²τ̂²·log τ

    F is strictly decreasing for m above m₁ = √3·(N/2)·τ̂/(π k*), where
    F(m₁) = −2N²τ̂² log τ > 0; the paper's upper bracket
    m₁·√(1 − (8/√3) log τ) is used, and doubled until F changes sign in the
    (rare) case it does not.  The root is found by bisection.

    Returns a float; :func:`mb_test_function_params` rounds up after
    averaging over axes, as ``findcorners.m`` does.
    """
    k_star = max(int(k_star), 1)

    def F(m):
        return (math.log((2.0 * m - 1.0) / m ** 2)
                * (4.0 * math.pi ** 2 * k_star ** 2 * m ** 2
                   - 3.0 * N ** 2 * tau_hat ** 2)
                - 2.0 * N ** 2 * tau_hat ** 2 * math.log(tau))

    lo = max(2.0, math.sqrt(3.0) * (N / 2.0) * tau_hat / (math.pi * k_star))
    hi = lo * math.sqrt(1.0 - (8.0 / math.sqrt(3.0)) * math.log(tau))
    while F(hi) > 0.0:
        lo, hi = hi, 2.0 * hi
    if F(lo) < 0.0:          # only possible when lo was clamped up to 2
        return min(lo, N / 2.0 - 1.0)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if F(mid) > 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-9:
            break
    return min(0.5 * (lo + hi), N / 2.0 - 1.0)


def mb_test_function_params(U, spatial_axes, time_axis=None, alpha_x=4,
                            alpha_t=1, tau=1e-10, tau_hat=2.0,
                            tau_hat_t=None):
    """(m, p) for space and time by the Messenger & Bortz rule.

    For each spatial axis: k* from :func:`mb_critical_wavenumber`, then the
    root of F; the spatial support is ``ceil(mean)`` over spatial axes, capped
    at ``floor((min N − 1)/2)`` (Algorithm 4.1's admissibility bound).  The
    same for time on ``time_axis``.  Degrees follow eq. (4.3):

        p = max(ᾱ + 1, ⌈log τ / log(1 − (1 − 1/m)²)⌉)

    Parameters
    ----------
    U : ndarray - the data actually regressed on (noise included)
    spatial_axes : sequence of int
    time_axis : int or None - skip the temporal parameters when None
    alpha_x, alpha_t : int - maximal derivative order moved onto the test
        function along each spatial axis / along time
    tau : float - real-space decay tolerance (paper: 1e-10)
    tau_hat : float - Fourier-tail parameter (paper: 1 to 3)
    tau_hat_t : float or None - separate τ̂ for time (default: ``tau_hat``)

    Returns
    -------
    dict with keys ``m``, ``p``, ``k_x`` (list, per spatial axis) and, when
    ``time_axis`` is given, ``m_t``, ``p_t``, ``k_t``.
    """
    spatial_axes = list(spatial_axes)
    k_x = [mb_critical_wavenumber(U, a) for a in spatial_axes]
    roots = [mb_root_m(k, U.shape[a], tau, tau_hat)
             for k, a in zip(k_x, spatial_axes)]
    n_min = min(U.shape[a] for a in spatial_axes)
    m = int(min((n_min - 1) // 2, math.ceil(np.mean(roots))))
    out = dict(m=m, p=_compute_p_from_m(m, alpha_x, tau), k_x=k_x)

    if time_axis is not None:
        Nt = U.shape[time_axis]
        k_t = mb_critical_wavenumber(U, time_axis)
        root_t = mb_root_m(k_t, Nt, tau,
                           tau_hat if tau_hat_t is None else tau_hat_t)
        m_t = int(min((Nt - 1) // 2, math.ceil(root_t)))
        out.update(m_t=m_t, p_t=_compute_p_from_m(m_t, alpha_t, tau), k_t=k_t)
    return out


# ---------------------------------------------------------------------------
# Rescaling via PDE scale invariance
#
# Messenger & Bortz (2021), Section 3.3.
# ---------------------------------------------------------------------------


def compute_rescaling_scales(
    U, dx, dt, m_x, m_t, p_x, p_t,
    max_deriv_spatial, max_deriv_temporal, max_poly_degree,
):
    r"""Compute coordinate and amplitude scaling factors.

    .. math::

        \gamma_x &= \frac{1}{m_x \Delta x}
            \binom{p_x}{\bar{\alpha}_x/2} \bar{\alpha}_x!^{1/\bar{\alpha}_x} \\
        \gamma_t &= \frac{1}{m_t \Delta t}
            \binom{p_t}{\bar{\alpha}_t/2} \bar{\alpha}_t!^{1/\bar{\alpha}_t} \\
        \gamma_u &= \left\|
            \frac{\mathbf{U}^{\bar{\beta}}}{\|\mathbf{U}\|_F}
            \right\|_F^{-1/\bar{\beta}}

    where :math:`\bar{\alpha}_x` = ``max_deriv_spatial``,
    :math:`\bar{\alpha}_t` = ``max_deriv_temporal``,
    :math:`\bar{\beta}` = ``max_poly_degree``.

    Parameters
    ----------
    U : ndarray
        Data array (any shape; raveled for γ_u computation).
    dx, dt : float
        Grid spacings.
    m_x, m_t : int
        Test-function support half-widths.
    p_x, p_t : int
        Test-function polynomial degrees.
    max_deriv_spatial : int
        Maximum spatial derivative order.
    max_deriv_temporal : int
        Maximum temporal derivative order.
    max_poly_degree : int
        Maximum polynomial degree (β̄).

    Returns
    -------
    gamma_x : float
    gamma_t : float
    gamma_u : float
    """
    # Spatial scale
    if max_deriv_spatial > 0:
        alpha_x = max_deriv_spatial
        gamma_x = (
            1.0
            / (m_x * dx)
            * (math.comb(p_x, alpha_x // 2) * math.factorial(alpha_x))
            ** (1.0 / alpha_x)
        )
    else:
        gamma_x = 1.0

    # Temporal scale
    if max_deriv_temporal > 0:
        alpha_t = max_deriv_temporal
        if alpha_t % 2 == 0 and p_t >= alpha_t // 2:
            gamma_t = (
                1.0
                / (m_t * dt)
                * (math.comb(p_t, alpha_t // 2) * math.factorial(alpha_t))
                ** (1.0 / alpha_t)
            )
        else:
            gamma_t = 1.0 / (m_t * dt)
    else:
        gamma_t = 1.0

    # Amplitude scale
    U_fro = float(np.linalg.norm(U.ravel()))
    if max_poly_degree > 0 and U_fro > 1e-30:
        U_pow = np.abs(U) ** max_poly_degree
        norm_pow = float(np.linalg.norm(U_pow.ravel()))
        gamma_u = (norm_pow / U_fro) ** (-1.0 / max_poly_degree)
    else:
        gamma_u = 1.0

    return float(gamma_x), float(gamma_t), float(gamma_u)


def build_rescale_matrix(gamma_x, gamma_t, gamma_u, max_p, max_n):
    r"""Build the diagonal rescaling matrix :math:`\mathbf{M}` for the
    equivariant multigraph library.

    For an operator with edge count ``p`` (differential order = 2p) and
    polynomial degree ``k``:

    .. math::

        \mu_j = \gamma_u^{-(k-1)} \cdot \gamma_x^{-2p} \cdot \gamma_t

    so that :math:`\hat{w}_j = \mu_j \cdot \tilde{w}_j` recovers the
    coefficients in the original (unscaled) coordinates.

    Parameters
    ----------
    gamma_x, gamma_t, gamma_u : float
        Scaling factors from ``compute_rescaling_scales``.
    max_p, max_n : int
        Library parameters (must match ``get_operators``).

    Returns
    -------
    M : ndarray, shape (K,)
        Diagonal entries of the rescaling matrix.
    """
    from ..multigraph_operators import get_operator_metadata

    meta = get_operator_metadata(max_p, max_n)
    mu = np.empty(len(meta))
    for j, (p, k) in enumerate(meta):
        # μ_j = γ_t · γ_u^{-(k-1)} · γ_x^{-2p}
        mu[j] = gamma_t * (gamma_u ** (-(k - 1))) * (gamma_x ** (-2 * p))
    return mu
