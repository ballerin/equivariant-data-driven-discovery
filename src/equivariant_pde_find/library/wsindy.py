"""
Pointwise WSINDy library construction.

Internal module — use ``library.build_wsindy_library()`` (the public
facade) instead of calling ``_build_wsindy_pointwise`` directly.
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


def _build_wsindy_pointwise(f, Lx, Ly, max_order=4, max_poly=3):
    """Build the pointwise WSINDy library from data.

    This is the pointwise analogue of the library used in WSINDy, for
    benchmark purposes only.

    For max_order=4 and max_poly=3 this gives 46 columns.

    Parameters
    ----------
    f         : ndarray shape (Nx, Ny) - data
    Lx, Ly    : float - domain lengths
    max_order : int - highest partial derivative order (default 4)
    max_poly  : int - highest monomial power q in ∂^α(f^q) (default 3)

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


def _build_wsindy_weak(f, Lx, Ly, max_order=4, max_poly=3,
                       m=8, p=4, use_ibp=False):
    """Build the weak-form WSINDy library: ψ ∗ ∂^α(f^q).

    Mirrors ``_build_wsindy_pointwise`` but convolves each column with the
    polynomial-bump test function ψ.  When ``use_ibp=True``, integration by
    parts is applied to linear terms (q=1), moving the derivative onto the
    test function:

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
    from .weak_form import (
        _cached_kernel_1d, _test_function_deriv_1d, weak_convolution,
    )

    Nx, Ny = f.shape
    kernel = _cached_kernel_1d(m, p)
    method = _derivatives.get_backend()

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
                desc = f"(f{_sup(q)})_{{{dx_part}}}"
            descs.append(desc)

    # Pure monomials f^q for q ≥ 1
    for q in range(1, max_poly + 1):
        col_weak = weak_convolution(f ** q, kernel, kernel)
        cols.append(col_weak.ravel())
        descs.append("f" if q == 1 else f"f{_sup(q)}")

    return np.column_stack(cols), descs
