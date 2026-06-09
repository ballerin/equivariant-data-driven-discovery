"""
Weak-form SINDy (WSINDy) for equivariant PDE discovery.

Implements the weak-form / Galerkin discretization of PDE regression from
Messenger & Bortz, J. Comput. Phys. (2021).  Instead of evaluating the
library at individual grid points, each operator is integrated against
compactly supported test functions:

    G_{ij,k} = ∫ ψ_{ij}(x,y) · N_k[f](x,y)  dx dy
    b_{ij}   = ∫ ψ_{ij}(x,y) · f_t(x,y)      dx dy

where ψ_{ij} is a polynomial bump centred at query point (i, j).  The
convolution provides natural noise averaging: integrating over a support
of width (2m+1)² grid points reduces the noise standard deviation by a
factor of 2m+1.

For the equivariant library, operators like ‖∇f‖² and Δf cannot all be
rewritten as ∂^α(f^β), so we use a **direct** approach:
  - Compute each operator N_k[f] on the full grid (spectral or polynomial
    backend) via the existing ``multigraph_operators`` functions,
  - Convolve with separable 2-D polynomial-bump test functions,
  - Subsample at query points to build (G, b).

This is a drop-in replacement for the pointwise Theta matrix: feed G and b
into TrainSTRidge / TrainSR3 exactly as you would Θ and f_t.

References
----------
.. [1] Messenger, D. A. & Bortz, D. M. (2021). Weak SINDy: Galerkin-based
   data-driven model selection. J. Comput. Phys., 443, 110525.
.. [2] Minor, S. (2025). PyWSINDy-for-PDEs.
   https://github.com/SethMinor/PyWSINDy-for-PDEs
"""

from functools import lru_cache
import math

import numpy as np

from multigraph_operators import get_operators


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

    Uses the **discrete** second difference of the test function, which
    makes the discrete summation-by-parts identity exact:

        Σ ψ·ψ · Δ_disc f = Σ (D2_ψ)·ψ · f + Σ ψ·(D2_ψ) · f

    where D2_ψ[j] = ψ[j-1] + ψ[j+1] - 2ψ[j].

    Returns
    -------
    kernel_xx : ndarray (2*m+1,) — discrete D2_φ (grid-point units)
    kernel_0  : ndarray (2*m+1,) — φ
    """
    phi = _test_function_1d(m, p, dx=1.0)
    # Discrete second difference: D2_phi[j] = phi[j-1] + phi[j+1] - 2*phi[j]
    d2 = np.roll(phi, 1) + np.roll(phi, -1) - 2.0 * phi
    return d2, phi


@lru_cache(maxsize=64)
def _ibp_bilaplacian_kernels(m, p):
    """Return (kernel_xxxx, kernel_xx, kernel_0) for IBP of the bilaplacian.

    Uses discrete differences for exact summation-by-parts:

        D2_φ[j]  = φ[j-1] + φ[j+1] - 2φ[j]
        D4_φ[j]  = D2_φ[j-1] + D2_φ[j+1] - 2·D2_φ[j]

    With separable ψ=φ(x)·φ(y):
        Δ²ψ = D4_φ(x)·φ(y) + 2·D2_φ(x)·D2_φ(y) + φ(x)·D4_φ(y)

    Returns
    -------
    kernel_xxxx : ndarray (2*m+1,) — discrete D4_φ
    kernel_xx   : ndarray (2*m+1,) — discrete D2_φ
    kernel_0    : ndarray (2*m+1,) — φ
    """
    phi = _test_function_1d(m, p, dx=1.0)
    d2 = np.roll(phi, 1) + np.roll(phi, -1) - 2.0 * phi
    d4 = np.roll(d2, 1) + np.roll(d2, -1) - 2.0 * d2
    return d4, d2, phi


# ---------------------------------------------------------------------------
# IBP weak-form evaluation for linear operators
# ---------------------------------------------------------------------------


def weak_laplacian_ibp(f, m, p, Lx, Ly):
    """Compute ∫ ψ · Δf ≈ ∫ Δψ · f  via discrete summation by parts.

    Uses the discrete second difference of ψ (D2_ψ[j] = ψ[j-1]+ψ[j+1]-2ψ[j]),
    which makes the discrete Green's identity exact.  The kernel is scaled
    by 1/dx² to match the physical units of the spectral Laplacian.

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
    kxx_grid, k0 = _ibp_laplacian_kernels(m, p)
    kxx = kxx_grid / dx2   # D2_φ → D2_φ/dx²  (physical units)
    kyy = kxx_grid / dy2
    return weak_convolution_2d(f, kxx, k0) + weak_convolution_2d(f, k0, kyy)


def weak_bilaplacian_ibp(f, m, p, Lx, Ly):
    """Compute ∫ ψ · Δ²f ≈ ∫ Δ²ψ · f  via discrete summation by parts.

    Uses discrete 2nd and 4th differences of ψ.  No derivatives of *f*
    are computed.

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
        weak_convolution_2d(f, kxxxx, k0)
        + 2.0 * weak_convolution_2d(f, kxx, kyy)
        + weak_convolution_2d(f, k0, kyyyy)
    )


# ---------------------------------------------------------------------------
# Weak convolution
# ---------------------------------------------------------------------------


def weak_convolution_2d(f, kernel_x, kernel_y):
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


def build_weak_library_2d(f, Lx, Ly, max_p=2, max_n=3, m=8, p=4,
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
        |k|² and |k|⁴ noise amplification.  Default False (pure WSINDy).

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
                col_weak = weak_convolution_2d(col_raw, kernel, kernel)
        else:
            # Direct path: compute operator on data, then convolve.
            col_raw = func(f, Lx, Ly)
            col_weak = weak_convolution_2d(col_raw, kernel, kernel)
        cols.append(col_weak.ravel())

    descriptions = [desc for _, desc in ops]
    return np.column_stack(cols), descriptions


def build_weak_standard_library_2d(f, Lx, Ly, max_order=4, max_poly=3,
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
    import derivatives as _derivatives
    from polynomial_derivatives import poly_partial

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
                # IBP: convolve f (not the derivative) with the IBP kernel
                # For ∂_x^a ∂_y^b: ∫ ψ · ∂_x^a∂_y^b f = (-1)^{a+b} ∫ ∂_x^a∂_y^b ψ · f
                # Build separable IBP kernel using discrete differences
                phi = _test_function_1d(m, p, dx=1.0)
                kx_ibp = _test_function_deriv_1d(m, p, a)
                ky_ibp = _test_function_deriv_1d(m, p, b)
                # Scale by grid spacing to match physical units
                dx_a = (Lx / Nx) ** a
                dy_b = (Ly / Ny) ** b
                sign = (-1) ** (a + b)
                kx_scaled = sign * kx_ibp / dx_a
                ky_scaled = ky_ibp / dy_b
                col_weak = weak_convolution_2d(f, kx_scaled, ky_scaled)
            else:
                col_weak = weak_convolution_2d(col_raw, kernel, kernel)

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


def build_weak_wsindy_library_2d(f, Lx, Ly, max_order=4, max_poly=3,
                                  m=8, p=4, use_ibp=False):
    """Build the weak-form WSINDy library: ψ ∗ ∂^α(f^q).

    Mirrors ``build_wsindy_library`` in ``library.py`` but convolves each
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
    import derivatives as _deriv
    from polynomial_derivatives import poly_partial

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
    col_weak = weak_convolution_2d(ones_field, kernel, kernel)
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
                kx_ibp = _test_function_deriv_1d(m, p, a)
                ky_ibp = _test_function_deriv_1d(m, p, b)
                dx_a = (Lx / Nx) ** a
                dy_b = (Ly / Ny) ** b
                sign = (-1) ** (a + b)
                kx_scaled = sign * kx_ibp / dx_a
                ky_scaled = ky_ibp / dy_b
                col_weak = weak_convolution_2d(f, kx_scaled, ky_scaled)
            else:
                # Direct: compute ∂^a_x ∂^b_y (f^q) then convolve
                if method == "spectral":
                    d_hat = (1j * KX) ** a * (1j * KY) ** b * fh_pow
                    col_raw = np.real(np.fft.ifft2(d_hat))
                else:
                    col_raw = poly_partial(f_pow, Lx, Ly, a, b)
                col_weak = weak_convolution_2d(col_raw, kernel, kernel)

            cols.append(col_weak.ravel())

            dx_part = "x" * a + "y" * b
            if q == 1:
                desc = f"f_{{{dx_part}}}"
            else:
                desc = f"(f{_sup_ws(q)})_{{{dx_part}}}"
            descs.append(desc)

    # Pure monomials f^q for q ≥ 1
    for q in range(1, max_poly + 1):
        col_weak = weak_convolution_2d(f ** q, kernel, kernel)
        cols.append(col_weak.ravel())
        descs.append("f" if q == 1 else f"f{_sup_ws(q)}")

    return np.column_stack(cols), descs


def _sup_ws(n):
    """Unicode superscript for n ≥ 2."""
    _sup = {2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}
    return _sup.get(n, f"^{n}")


def build_weak_polynomial_jet_library_2d(f, Lx, Ly, max_order=4, max_degree=3,
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
    from library import _enumerate_jet_monomials

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
    import derivatives as _deriv
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
            from polynomial_derivatives import poly_partial
            jet_fields[(a, b)] = poly_partial(f, Lx, Ly, a, b)

    # Enumerate monomials
    monomials = _enumerate_jet_monomials(jet_vars, max_order, max_degree)

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
            kx_ibp = _test_function_deriv_1d(m, p, a)
            ky_ibp = _test_function_deriv_1d(m, p, b)
            dx_a = (Lx / Nx) ** a
            dy_b = (Ly / Ny) ** b
            sign = (-1) ** (a + b)
            kx_scaled = sign * kx_ibp / dx_a
            ky_scaled = ky_ibp / dy_b
            col_weak = weak_convolution_2d(f, kx_scaled, ky_scaled)
        else:
            # Build pointwise monomial
            col_raw = np.ones_like(f)
            for (a, b), exp in exponents.items():
                if exp > 0:
                    col_raw = col_raw * (jet_fields[(a, b)] ** exp)
            col_weak = weak_convolution_2d(col_raw, kernel, kernel)

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
                parts.append(f"{var_name}{_sup_ws(exp)}")
        descs.append("1" if not parts else "·".join(parts))

    return np.column_stack(cols), descs


def build_weak_rhs_2d(ft, kernel_x=None, kernel_y=None, m=8, p=4):
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
    return weak_convolution_2d(ft, kernel_x, kernel_y).ravel()


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
        Test function parameters.
    sample_stride : int
        If > 1, subsample the weak grid uniformly with given stride.
        idx_sets is ignored when sample_stride > 1.
    use_ibp : bool
        Forwarded to ``build_weak_library_2d``.

    Returns
    -------
    Theta_all : ndarray of shape (n_total, K)
    f_t_all : ndarray of shape (n_total,)
    descriptions : list of K str

    Notes
    -----
    When ``sample_stride > 1``, a uniform strided grid of query points is
    used instead of random indices.  This is the preferred mode for WSINDy
    because the convolution already decorrelates adjacent points.
    """
    kernel = _cached_kernel_1d(m, p)

    Theta_parts, rhs_parts = [], []
    descriptions = None

    for snap_idx in range(len(f_snaps)):
        f = f_snaps[snap_idx]
        Theta_w, descs = build_weak_library_2d(
            f, Lx, Ly, max_p, max_n, m=m, p=p, use_ibp=use_ibp,
        )
        if descriptions is None:
            descriptions = descs

        b_w = build_weak_rhs_2d(ft_snaps[snap_idx], kernel_x=kernel, kernel_y=kernel)

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
    return Theta_all, f_t_all, descriptions


# ---------------------------------------------------------------------------
# Convenience: full pipeline from simulation output
# ---------------------------------------------------------------------------


def wsindy_from_snapshots(
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
):
    """One-stop: produce (G, b) for WSINDy regression from raw snapshots.

    This wraps ``make_regression_data`` and ``make_weak_regression_data``
    into a single call.

    Parameters
    ----------
    snapshots : ndarray of shape (n_saved, Nx, Ny)
    times : ndarray of shape (n_saved,)
    Lx, Ly : float
    max_p, max_n : int
    m, p : int
    n_per_snap : int
        Ignored when sample_stride > 1.
    skip_fraction : float
    sample_stride : int
    seed : int
    use_ibp : bool
        Forwarded to ``build_weak_library_2d``.

    Returns
    -------
    G : ndarray of shape (n_total, K)
    b : ndarray of shape (n_total,)
    descriptions : list of K str
    """
    from simulations import make_regression_data

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
    )
