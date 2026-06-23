"""
Weak-form PDE discovery via convolution with compactly-supported test
functions.

Instead of evaluating the library at individual grid points, each operator
is integrated against compactly-supported test functions ψ:

    G_{ij,k} = ∫ ψ_{ij}(x,y) · N_k[f](x,y)  dx dy
    b_{ij}   = ∫ ψ_{ij}(x,y) · f_t(x,y)      dx dy

where ψ_{ij} is centred at query point (i, j).  The convolution provides
natural noise averaging: integrating over a support of width (2m+1)² grid
points reduces the noise standard deviation by a factor of 2m+1.

Supports four library types in weak form:
  - Equivariant (multigraph-operator) library via ``build_weak_library``
  - Standard (PDE-FIND) library via ``build_weak_standard_library``
  - WSINDy library (∂^α(f^q)) via ``build_weak_wsindy_library``
  - Polynomial jet library via ``build_weak_polynomial_jet_library``

For the equivariant library, operators like ‖∇f‖² and Δf cannot all be
rewritten as ∂^α(f^β), so we use a **direct** approach:
  - Compute each operator N_k[f] on the full grid (spectral or polynomial
    backend) via the existing ``multigraph_operators`` functions,
  - Convolve with separable polynomial-bump test functions,
  - Subsample at query points to build (G, b).

Regression is then ran by feeding G and b which replace the pointwise Θ and f_t.

References
----------
.. [1] Messenger, D. A. & Bortz, D. M. (2021). Weak SINDy: Galerkin-based
   data-driven model selection. J. Comput. Phys., 443, 110525.
"""

from functools import lru_cache
import math

import numpy as np

from .multigraph_operators import get_operators


# ---------------------------------------------------------------------------
# Test function construction
# ---------------------------------------------------------------------------


def _test_function_1d(m, p, dx=1.0):
    """Build a 1-D polynomial bump test function supported on [-m·dx, m·dx].

    ψ(x) = (1 - (x / (m·dx))²)^p   for |x| < m·dx,  0 otherwise.

    The function is C^{p-1} at the boundary and integrates to a known
    normalisation constant.

    Parameters
    ----------
    m : int
        Half-support (number of grid points to each side of centre).
    p : int
        Polynomial degree.  Higher p → smoother, broader effective support.
    dx : float
        Grid spacing.  Default 1.0 returns kernel values on integer grid.

    Returns
    -------
    kernel : ndarray of shape (2*m + 1,)
        Values of ψ at grid points −m·dx, …, 0, …, m·dx.
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
    scaled by 1/m^order to account for d/dx = (1/(m·dx)) · d/dt with
    dx = 1 (grid-point units).  The result is exact (no finite-difference
    approximation) and serves as the IBP kernel for linear differential
    operators.

    Parameters
    ----------
    m : int
    p : int
    order : int (≥ 0; 0 returns φ itself)

    Returns
    -------
    kernel : ndarray of shape (2*m + 1,)
    """
    if order == 0:
        return _test_function_1d(m, p, dx=1.0)
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
    phi = _test_function_1d(m, p, dx=1.0)
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
    phi = _test_function_1d(m, p, dx=1.0)
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


def weak_laplacian_ibp(f, m, p, Lx, Ly):
    """Compute ∫ ψ · Δf ≈ ∫ Δψ · f  via integration by parts.

    Uses the analytic second derivative of the test function φ (computed
    from the exact polynomial expansion), scaled to physical units.  This
    gives ~0.05 % relative error against the spectral direct method (vs
    ~0.86 % for the discrete 3-point stencil).

    No derivatives of *f* are computed — the differential operator is
    moved entirely onto the analytically-known test function.

    Parameters
    ----------
    f : ndarray of shape (Nx, Ny)
    m, p : int — test function parameters
    Lx, Ly : float — domain lengths

    Returns
    -------
    ndarray of shape (Nx − 2*m, Ny − 2*m)
    """
    Nx, Ny = f.shape
    dx2 = (Lx / Nx) ** 2
    dy2 = (Ly / Ny) ** 2
    kxx_grid, k0 = _ibp_laplacian_kernels(m, p)
    kxx = kxx_grid / dx2   # analytic φ'' → φ''/dx² (physical units)
    kyy = kxx_grid / dy2
    return weak_convolution(f, kxx, k0) + weak_convolution(f, k0, kyy)


def weak_bilaplacian_ibp(f, m, p, Lx, Ly):
    """Compute ∫ ψ · Δ²f ≈ ∫ Δ²ψ · f  via discrete summation by parts.

    Uses discrete 4th differences of ψ (D4 = D2∘D2).  For p=4 the bump
    function is C³, so its 4th derivative has a jump at the boundary;
    the discrete kernel is preferred over the analytic one in this case.
    Relative error against the spectral direct method is ~5.6 % (O(dx²)).

    No derivatives of *f* are computed.

    Parameters
    ----------
    f : ndarray of shape (Nx, Ny)
    m, p : int — test function parameters
    Lx, Ly : float — domain lengths

    Returns
    -------
    ndarray of shape (Nx − 2*m, Ny − 2*m)
    """
    Nx, Ny = f.shape
    dx2 = (Lx / Nx) ** 2
    dy2 = (Ly / Ny) ** 2
    kxxxx_grid, kxx_grid, k0 = _ibp_bilaplacian_kernels(m, p)
    # Δ²ψ = D4_φ(x)·φ(y)/dx⁴ + 2·D2_φ(x)·D2_φ(y)/(dx²·dy²) + φ(x)·D4_φ(y)/dy⁴
    kxxxx = kxxxx_grid / (dx2 ** 2)
    kxx   = kxx_grid   / dx2
    kyyyy = kxxxx_grid / (dy2 ** 2)
    kyy   = kxx_grid   / dy2
    return (
        weak_convolution(f, kxxxx, k0)
        + 2.0 * weak_convolution(f, kxx, kyy)
        + weak_convolution(f, k0, kyyyy)
    )


# ---------------------------------------------------------------------------
# Weak convolution
# ---------------------------------------------------------------------------


def weak_convolution(f, kernel_x, kernel_y):
    """Compute the 2-D weak form of field `f` by separable 1-D convolution.

    The convolution is

        f_weak[i,j] = Σ_{di} Σ_{dj} ker_x[di] · ker_y[dj] · f[i+di, j+dj]

    with ``mode='valid'``, so the output has shape
    ``(Nx - 2*mx, Ny - 2*my)`` where ``mx = len(kernel_x)//2``.

    Uses ``np.convolve`` applied separably: first along axis=0 (x), then
    axis=1 (y).  No scipy dependency.

    Parameters
    ----------
    f : ndarray of shape (Nx, Ny)
    kernel_x : ndarray of shape (2*mx + 1,)
    kernel_y : ndarray of shape (2*my + 1,)

    Returns
    -------
    f_weak : ndarray of shape (Nx - 2*mx, Ny - 2*my)
    """
    kx = kernel_x.ravel()
    ky = kernel_y.ravel()
    mx, my = len(kx) // 2, len(ky) // 2

    # Convolve along axis=0 (rows) for each column, then axis=1 for each row.
    Nx, Ny = f.shape
    # Step 1: convolve each column (axis=0) with kernel_x
    temp = np.empty((Nx - 2 * mx, Ny), dtype=np.float64)
    for j in range(Ny):
        temp[:, j] = np.convolve(f[:, j], kx, mode="valid")
    # Step 2: convolve each row (axis=1) with kernel_y
    result = np.empty((Nx - 2 * mx, Ny - 2 * my), dtype=np.float64)
    for i in range(temp.shape[0]):
        result[i, :] = np.convolve(temp[i, :], ky, mode="valid")
    return result


# ---------------------------------------------------------------------------
# Cache kernel-building (same (m, p, dx) → same kernel)
# ---------------------------------------------------------------------------

@lru_cache(maxsize=64)
def _cached_kernel_1d(m, p):
    return _test_function_1d(m, p, dx=1.0)


# ---------------------------------------------------------------------------
# Weak library construction (single snapshot)
# ---------------------------------------------------------------------------


# Operators that can be handled via IBP (purely linear differential operators).
_IBP_OPERATORS = {"Δf", "Δ²f"}


def build_weak_library(f, Lx, Ly, max_p=2, max_n=3, m=8, p=4,
                          use_ibp=False):
    """Build the weak-form library matrix for a single 2-D snapshot.

    Parameters
    ----------
    f : ndarray of shape (Nx, Ny)
    Lx, Ly : float
    max_p, max_n : int
        Passed to ``get_operators()``.
    m : int
        Test function half-support in grid points.  The test function kernel
        has length 2*m + 1.
    p : int
        Test function polynomial degree.
    use_ibp : bool
        If True, use integration by parts for linear differential operators
        (Δf, Δ²f).  This moves derivatives off the (possibly noisy) data
        and onto the analytically-known test function, eliminating the
        |k|² and |k|⁴ noise amplification.  Default False (direct weak-form).

    Returns
    -------
    Theta_weak : ndarray of shape (N_valid, K)
        Weak library matrix.  N_valid = (Nx - 2*m) * (Ny - 2*m).
    descriptions : list of K str
    """
    ops = get_operators(max_p, max_n)
    kernel = _cached_kernel_1d(m, p)

    cols = []
    for func, desc in ops:
        if use_ibp and desc in _IBP_OPERATORS:
            # IBP path: convolve raw f with differentiated test function.
            # No derivatives of f are computed — the differential operator
            # is moved onto ψ via Green's identities.
            if desc == "Δf":
                col_weak = weak_laplacian_ibp(f, m, p, Lx, Ly)
            elif desc == "Δ²f":
                col_weak = weak_bilaplacian_ibp(f, m, p, Lx, Ly)
            else:
                # Should not reach here, but fall back to direct.
                col_raw = func(f, Lx, Ly)
                col_weak = weak_convolution(col_raw, kernel, kernel)
        else:
            # Direct path: compute operator on data, then convolve.
            col_raw = func(f, Lx, Ly)
            col_weak = weak_convolution(col_raw, kernel, kernel)
        cols.append(col_weak.ravel())

    descriptions = [desc for _, desc in ops]
    return np.column_stack(cols), descriptions


def build_weak_standard_library(f, Lx, Ly, max_order=4, max_poly=3,
                                   m=8, p=4, use_ibp=False):
    """Build the weak-form **standard** (coordinate-basis) library matrix.

    Mirrors ``build_standard_library`` in ``library.py`` but convolves each
    column with the test function ψ.  IBP is applied to purely linear
    differential operators (∂_x, ∂_y, ∂_xx, ∂_xy, ∂_yy, …) when
    ``use_ibp=True``.

    Parameters
    ----------
    f : ndarray of shape (Nx, Ny)
    Lx, Ly : float
    max_order : int — highest partial derivative order (default 4)
    max_poly  : int — highest polynomial degree factor (default 3)
    m, p : int — test function parameters
    use_ibp : bool — enable IBP for linear partial-derivative columns

    Returns
    -------
    Theta_weak : ndarray of shape (N_valid, K)
    descriptions : list of K str
    """
    from . import derivatives as _derivatives
    from .polynomial_derivatives import poly_partial

    Nx, Ny = f.shape
    kernel = _cached_kernel_1d(m, p)
    method = _derivatives.get_backend()

    # Build list of (a, b) multi-indices up to max_order
    multi_indices = []
    for total in range(max_order + 1):
        for a in range(total + 1):
            multi_indices.append((a, total - a))

    # Compute each pure derivative term on the grid
    derivs, deriv_names = [], []
    for a, b in multi_indices:
        if a + b == 0:
            deriv = np.ones_like(f)
        elif method == "spectral":
            kx = 2.0 * np.pi / Lx * np.fft.fftfreq(Nx) * Nx
            ky = 2.0 * np.pi / Ly * np.fft.fftfreq(Ny) * Ny
            KX, KY = np.meshgrid(kx, ky, indexing="ij")
            fh = np.fft.fft2(f)
            d_hat = (1j * KX) ** a * (1j * KY) ** b * fh
            deriv = np.real(np.fft.ifft2(d_hat))
        else:
            deriv = poly_partial(f, Lx, Ly, a, b)
        name = (("x" * a) + ("y" * b)) if (a + b) > 0 else "1"
        name = f"f_{name}" if (a + b) > 0 else "1"
        derivs.append(deriv)
        deriv_names.append(name)

    # Build columns: polynomial prefactor × derivative, then convolve
    cols, descs = [], []
    N_valid = (Nx - 2 * m) * (Ny - 2 * m)

    for poly_deg in range(max_poly + 1):
        for deriv, dname in zip(derivs, deriv_names):
            poly_field = f ** poly_deg if poly_deg > 0 else np.ones_like(f)
            col_raw = poly_field * deriv

            # Determine if this is a pure linear differential term
            # (poly_deg == 0 and a+b > 0) — candidates for IBP
            a_plus_b = 0
            if dname != "1" and dname.startswith("f_"):
                dx_part = dname[2:]  # e.g. "xxy"
                a_plus_b = len(dx_part)

            if use_ibp and poly_deg == 0 and a_plus_b > 0:
                # IBP: convolve f (not the derivative) with the IBP kernel.
                # For ∂_x^a ∂_y^b: ∫ ψ · ∂_x^a∂_y^b f = (-1)^{a+b} ∫ ∂_x^a∂_y^b ψ · f
                # The (-1)^{a+b} sign is handled implicitly by np.convolve
                # (which reverses the kernel; for odd-order anti-symmetric
                # kernels this introduces the required sign flip).
                # Extract (a, b) from the derivative name (e.g. "f_xxy" → a=2, b=1).
                a_ibp = dx_part.count("x")
                b_ibp = dx_part.count("y")
                kx_ibp = _test_function_deriv_1d(m, p, a_ibp)
                ky_ibp = _test_function_deriv_1d(m, p, b_ibp)
                # Scale by grid spacing to match physical units
                dx_a = (Lx / Nx) ** a_ibp
                dy_b = (Ly / Ny) ** b_ibp
                kx_scaled = kx_ibp / dx_a
                ky_scaled = ky_ibp / dy_b
                col_weak = weak_convolution(f, kx_scaled, ky_scaled)
            else:
                col_weak = weak_convolution(col_raw, kernel, kernel)

            cols.append(col_weak.ravel())

            # Build description matching standard library convention
            if poly_deg == 0:
                desc = dname
            elif poly_deg == 1:
                desc = f"f·{dname}" if dname != "1" else "f"
            else:
                sup = {2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶"}.get(poly_deg, f"^{poly_deg}")
                desc = f"f{sup}·{dname}" if dname != "1" else f"f{sup}"
            descs.append(desc)

    return np.column_stack(cols), descs


def build_weak_wsindy_library(f, Lx, Ly, max_order=4, max_poly=3,
                                  m=8, p=4, use_ibp=False):
    """Build the weak-form WSINDy library: ψ ∗ ∂^α(f^q).

    Mirrors ``build_pointwise_wsindy_library`` in ``library.py`` but convolves each
    column with the polynomial-bump test function ψ.  When ``use_ibp=True``,
    integration by parts is applied to linear terms (q=1), moving the
    derivative onto the test function:

        ∫ ψ · ∂^α f  =  (-1)^{|α|} ∫ ∂^α ψ · f

    For nonlinear terms (q ≥ 2), the pointwise derivative is computed
    using the active backend, then convolved with ψ (direct path).

    Parameters
    ----------
    f : ndarray of shape (Nx, Ny)
    Lx, Ly : float
    max_order : int — highest derivative order (default 4)
    max_poly  : int — highest monomial power q (default 3)
    m, p : int — test function parameters
    use_ibp : bool — enable IBP for q=1 columns

    Returns
    -------
    Theta_weak : ndarray of shape (N_valid, K)
    descriptions : list of K str
    """
    from . import derivatives as _deriv
    from .polynomial_derivatives import poly_partial

    Nx, Ny = f.shape
    kernel = _cached_kernel_1d(m, p)
    method = _deriv.get_backend()

    # Enumerate derivative multi-indices
    multi_indices = []
    for total in range(max_order + 1):
        for a in range(total + 1):
            multi_indices.append((a, total - a))

    # Pre-compute spectral wavenumbers if needed
    if method == "spectral":
        kx = 2.0 * np.pi / Lx * np.fft.fftfreq(Nx) * Nx
        ky = 2.0 * np.pi / Ly * np.fft.fftfreq(Ny) * Ny
        KX, KY = np.meshgrid(kx, ky, indexing="ij")

    cols, descs = [], []

    # Constant term: ψ ∗ 1
    ones_field = np.ones_like(f)
    col_weak = weak_convolution(ones_field, kernel, kernel)
    cols.append(col_weak.ravel())
    descs.append("1")

    for q in range(1, max_poly + 1):
        f_pow = f ** q
        if method == "spectral":
            fh_pow = np.fft.fft2(f_pow)

        for a, b in multi_indices:
            if a + b == 0:
                continue  # pure monomials handled below

            if use_ibp and q == 1:
                # IBP: ∫ ψ · ∂_x^a ∂_y^b f = (-1)^{a+b} ∫ ∂_x^a ∂_y^b ψ · f
                # The (-1)^{a+b} sign is handled implicitly by np.convolve
                # (which reverses the kernel; for odd-order anti-symmetric
                # kernels this introduces the required sign flip).
                kx_ibp = _test_function_deriv_1d(m, p, a)
                ky_ibp = _test_function_deriv_1d(m, p, b)
                dx_a = (Lx / Nx) ** a
                dy_b = (Ly / Ny) ** b
                kx_scaled = kx_ibp / dx_a
                ky_scaled = ky_ibp / dy_b
                col_weak = weak_convolution(f, kx_scaled, ky_scaled)
            else:
                # Direct: compute ∂^a_x ∂^b_y (f^q) then convolve
                if method == "spectral":
                    d_hat = (1j * KX) ** a * (1j * KY) ** b * fh_pow
                    col_raw = np.real(np.fft.ifft2(d_hat))
                else:
                    col_raw = poly_partial(f_pow, Lx, Ly, a, b)
                col_weak = weak_convolution(col_raw, kernel, kernel)

            cols.append(col_weak.ravel())

            dx_part = "x" * a + "y" * b
            if q == 1:
                desc = f"f_{{{dx_part}}}"
            else:
                desc = f"(f{sup(q)})_{{{dx_part}}}"
            descs.append(desc)

    # Pure monomials f^q for q ≥ 1
    for q in range(1, max_poly + 1):
        col_weak = weak_convolution(f ** q, kernel, kernel)
        cols.append(col_weak.ravel())
        descs.append("f" if q == 1 else f"f{sup(q)}")

    return np.column_stack(cols), descs


from .library import sup


def build_weak_polynomial_jet_library(f, Lx, Ly, max_order=4, max_degree=3,
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
    from .library import enumerate_jet_monomials

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
    from . import derivatives as _deriv
    method = _deriv.get_backend()

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
            from .polynomial_derivatives import poly_partial
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
                parts.append(f"{var_name}{sup(exp)}")
        descs.append("1" if not parts else "·".join(parts))

    return np.column_stack(cols), descs


def build_weak_rhs(ft, kernel_x=None, kernel_y=None, m=8, p=4):
    """Build the weak-form RHS vector b = ψ ∗ f_t.

    Parameters
    ----------
    ft : ndarray of shape (Nx, Ny)
        Time derivative of f.
    kernel_x, kernel_y : ndarray or None
        If None, build from (m, p).
    m, p : int
        Used only if kernels are None.

    Returns
    -------
    b : ndarray of shape (N_valid,)
    """
    if kernel_x is None:
        kernel_x = _cached_kernel_1d(m, p)
    if kernel_y is None:
        kernel_y = _cached_kernel_1d(m, p)
    return weak_convolution(ft, kernel_x, kernel_y).ravel()


# ---------------------------------------------------------------------------
# Weak regression data (multiple snapshots + optional subsampling)
# ---------------------------------------------------------------------------


def make_weak_regression_data(
    f_snaps,
    ft_snaps,
    idx_sets,
    Lx,
    Ly,
    max_p=2,
    max_n=3,
    m=8,
    p=4,
    sample_stride=1,
    use_ibp=False,
    rescale=False,
):
    """Build weak-form regression matrices from a stack of snapshots.

    This mirrors ``make_regression_data`` in ``simulations.py``:
    for each snapshot, the weak library and weak RHS are computed, and
    ``n_per_snap`` query points are drawn from the valid interior.

    Parameters
    ----------
    f_snaps : ndarray of shape (n_snaps, Nx, Ny)
    ft_snaps : ndarray of shape (n_snaps, Nx, Ny)
    idx_sets : list of ndarray
        Flat indices into the valid (N_valid,) domain for each snapshot.
    Lx, Ly : float
    max_p, max_n : int
    m, p : int
        Test function parameters.  Ignored when ``rescale=True`` and
        ``auto_test_function_params`` is able to estimate them from data.
    sample_stride : int
        If > 1, subsample the weak grid uniformly with given stride.
        idx_sets is ignored when sample_stride > 1.
    use_ibp : bool
        Forwarded to ``build_weak_library``.
    rescale : bool
        If True, compute scale-invariance factors (γ_x, γ_t, γ_u) from the
        first snapshot and apply column/LHS rescaling.  The returned
        ``Theta_all`` and ``f_t_all`` are scaled so that the solved
        coefficients must be unscaled via ``rescale_coefficients``.

    Returns
    -------
    Theta_all : ndarray of shape (n_total, K)
    f_t_all : ndarray of shape (n_total,)
    descriptions : list of K str
    rescale_info : dict or None
        ``{"col_scale": ..., "lhs_scale": ..., "meta": ...}`` if
        ``rescale=True``, else ``None``.
    """
    kernel = _cached_kernel_1d(m, p)

    Theta_parts, rhs_parts = [], []
    descriptions = None

    for snap_idx in range(len(f_snaps)):
        f = f_snaps[snap_idx]
        Theta_w, descs = build_weak_library(
            f, Lx, Ly, max_p, max_n, m=m, p=p, use_ibp=use_ibp,
        )
        if descriptions is None:
            descriptions = descs

        b_w = build_weak_rhs(ft_snaps[snap_idx], kernel_x=kernel, kernel_y=kernel)

        Nx, Ny = f.shape
        Nv_x, Nv_y = Nx - 2 * m, Ny - 2 * m

        if sample_stride > 1:
            # Strided uniform grid of query points.
            ix = np.arange(0, Nv_x, sample_stride)
            iy = np.arange(0, Nv_y, sample_stride)
            IX, IY = np.meshgrid(ix, iy, indexing="ij")
            flat_idx = (IX * Nv_y + IY).ravel()
            Theta_parts.append(Theta_w[flat_idx, :])
            rhs_parts.append(b_w[flat_idx])
        else:
            # Random indices (provided by caller).
            idx = idx_sets[snap_idx]
            Theta_parts.append(Theta_w[idx, :])
            rhs_parts.append(b_w[idx])

    Theta_all = np.vstack(Theta_parts)
    f_t_all = np.concatenate(rhs_parts)

    rescale_info = None
    if rescale:
        # Compute rescaling factors from first snapshot
        dx = Lx / f_snaps[0].shape[0]
        gamma_x, gamma_t, gamma_u = compute_rescaling_scales(
            f_snaps[0], dx, 1.0, m, m, p, p,
            max_deriv_spatial=2 * max_p,
            max_deriv_temporal=1,
            max_poly_degree=max_n,
        )
        col_scale, lhs_scale = _rescale_factors(
            gamma_x, gamma_t, gamma_u, max_p, max_n,
        )
        # Apply column scaling: G_scaled = G * diag(col_scale)
        Theta_all = Theta_all * col_scale[np.newaxis, :]
        # Apply LHS scaling: b_scaled = lhs_scale * b
        f_t_all = lhs_scale * f_t_all
        rescale_info = {
            "gamma_x": gamma_x,
            "gamma_t": gamma_t,
            "gamma_u": gamma_u,
            "col_scale": col_scale,
            "lhs_scale": lhs_scale,
            "max_p": max_p,
            "max_n": max_n,
        }

    return Theta_all, f_t_all, descriptions, rescale_info


def _rescale_factors(gamma_x, gamma_t, gamma_u, max_p, max_n):
    """Compute per-column and LHS scaling factors for rescaling.

    Returns
    -------
    col_scale : ndarray, shape (K,)
        Column scaling factors: G_scaled = G * diag(col_scale).
    lhs_scale : float
        LHS scaling factor: b_scaled = lhs_scale * b.
    """
    from .multigraph_operators import get_operator_metadata

    meta = get_operator_metadata(max_p, max_n)
    col_scale = np.empty(len(meta))
    for j, (p, k) in enumerate(meta):
        # Column j: N_j[ũ] = γ_u^{-k} γ_x^{-2p} N_j[u]
        col_scale[j] = gamma_u ** (-k) * gamma_x ** (-2 * p)
    # LHS: ∂_t̃ ũ = γ_u^{-1} γ_t^{-1} ∂_t u
    lhs_scale = gamma_u ** (-1) * gamma_t ** (-1)
    return col_scale, lhs_scale


def rescale_coefficients(w_scaled, rescale_info):
    """Recover original-scale coefficients from scaled solution.

    Parameters
    ----------
    w_scaled : ndarray, shape (K,)
        Coefficients solved on rescaled (G_scaled, b_scaled).
    rescale_info : dict
        Dictionary returned by ``make_weak_regression_data`` when
        ``rescale=True``.

    Returns
    -------
    w : ndarray, shape (K,)
        Coefficients in original (unscaled) coordinates.
    """
    col_scale = np.asarray(rescale_info["col_scale"])
    lhs_scale = float(rescale_info["lhs_scale"])
    return w_scaled * col_scale / lhs_scale


# ---------------------------------------------------------------------------
# Convenience: full pipeline from simulation output
# ---------------------------------------------------------------------------


def weak_form_from_snapshots(
    snapshots,
    times,
    Lx,
    Ly,
    max_p=2,
    max_n=3,
    m=8,
    p=4,
    n_per_snap=200,
    skip_fraction=0.0,
    sample_stride=1,
    seed=0,
    use_ibp=False,
    rescale=False,
    auto_params=False,
):
    """One-stop: produce (G, b) for weak-form regression from raw snapshots.

    This wraps ``make_regression_data`` and ``make_weak_regression_data``
    into a single call.

    Parameters
    ----------
    snapshots : ndarray of shape (n_saved, Nx, Ny)
    times : ndarray of shape (n_saved,)
    Lx, Ly : float
    max_p, max_n : int
    m, p : int
        Ignored when ``auto_params=True``.
    n_per_snap : int
        Ignored when sample_stride > 1.
    skip_fraction : float
    sample_stride : int
    seed : int
    use_ibp : bool
        Forwarded to ``build_weak_library``.
    rescale : bool
        If True, apply scale-invariance rescaling to G and b.
    auto_params : bool
        If True, estimate (m, p) from the data via the changepoint method.

    Returns
    -------
    G : ndarray of shape (n_total, K)
    b : ndarray of shape (n_total,)
    descriptions : list of K str
    rescale_info : dict or None
    """
    from .simulations import make_regression_data

    # Auto-select test function params from the data
    if auto_params:
        dx = Lx / snapshots.shape[1]
        dt_est = float(np.mean(np.diff(times))) if len(times) > 1 else 1.0
        m_x, p_x, m_t, p_t = auto_test_function_params(
            snapshots, dx, dt_est, max_deriv_spatial=2 * max_p,
        )
        m, p = m_x, p_x  # use spatial params for the 2D test functions

    f_snaps, ft_snaps, idx_sets, _dt = make_regression_data(
        snapshots, times, n_per_snap=n_per_snap, skip_fraction=skip_fraction, seed=seed
    )
    return make_weak_regression_data(
        f_snaps,
        ft_snaps,
        idx_sets,
        Lx,
        Ly,
        max_p=max_p,
        max_n=max_n,
        m=m,
        p=p,
        sample_stride=sample_stride,
        use_ibp=use_ibp,
        rescale=rescale,
    )


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
    from .multigraph_operators import get_operator_metadata

    meta = get_operator_metadata(max_p, max_n)
    mu = np.empty(len(meta))
    for j, (p, k) in enumerate(meta):
        # μ_j = γ_t · γ_u^{-(k-1)} · γ_x^{-2p}
        mu[j] = gamma_t * (gamma_u ** (-(k - 1))) * (gamma_x ** (-2 * p))
    return mu
