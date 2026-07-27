"""
Pointwise PDE-FIND (coordinate-based) library construction.

Internal module — use ``library.build_pdefind_library()`` (the public
facade) instead of calling ``_build_pdefind_pointwise`` directly.
"""

import numpy as np

from .. import derivatives as _derivatives
from ..derivatives.polynomial import poly_partial, poly_partial_1d


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


def _wavenumbers_3d(Nx, Ny, Nz, Lx, Ly, Lz):
    """Return 3-D wavenumber grids for an (Nx, Ny, Nz) FFT grid."""
    kx = 2.0 * np.pi / Lx * np.fft.fftfreq(Nx) * Nx
    ky = 2.0 * np.pi / Ly * np.fft.fftfreq(Ny) * Ny
    kz = 2.0 * np.pi / Lz * np.fft.fftfreq(Nz) * Nz
    return np.meshgrid(kx, ky, kz, indexing="ij")


def _build_pdefind_pointwise(f, Lx, Ly, max_order=4, max_poly=3, Lz=None):
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
    f         : ndarray, shape (Nx, Ny) - data
    Lx, Ly    : float - domain lengths
    max_order : int - highest partial derivative order (default 4)
    max_poly  : int - highest polynomial degree factor (default 3)

    Returns
    -------
    Theta        : ndarray, shape (Nx*Ny, n_cols)
    descriptions : list of strings
    """
    method = _derivatives.get_backend()
    if f.ndim == 3:
        Nx, Ny, Nz = f.shape
        if method == "spectral":
            KX, KY, KZ = _wavenumbers_3d(Nx, Ny, Nz, Lx, Ly, Lz)
            fh = np.fft.fftn(f)
        # 3-D multi-indices: all (a,b,c) with a+b+c <= max_order
        multi_indices = []
        for total in range(max_order + 1):
            for a in range(total + 1):
                for b in range(total - a + 1):
                    multi_indices.append((a, b, total - a - b))

        derivs, deriv_names = [], []
        for a, b, c in multi_indices:
            if a + b + c == 0:
                deriv = np.ones_like(f)
            elif method == "spectral":
                d_hat = (1j * KX) ** a * (1j * KY) ** b * (1j * KZ) ** c * fh
                deriv = np.real(np.fft.ifftn(d_hat))
            else:
                raise NotImplementedError(
                    "Polynomial backend does not support 3-D fields yet."
                )
            name = ("x" * a) + ("y" * b) + ("z" * c) if (a + b + c) > 0 else "1"
            name = "f_" + name if (a + b + c) > 0 else "1"
            derivs.append(deriv)
            deriv_names.append(name)
    else:
        # 2-D path (original)
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
            name = ("x" * a) + ("y" * b) if (a + b) > 0 else "1"
            name = "f_" + name if (a + b) > 0 else "1"
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


def _build_pdefind_pointwise_1d(u, L, max_deriv=3, max_poly=3):
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


def _build_pdefind_weak(f, Lx, Ly, max_order=4, max_poly=3,
                        m=8, p=4, use_ibp=False, Lz=None):
    """Build the weak-form PDE-FIND library matrix.

    Supports 2-D and 3-D fields.  In 3-D, uses 3-D separable convolution
    and extends IBP to ∂_x^a ∂_y^b ∂_z^c terms.
    """
    from .weak_form import (
        _cached_kernel_1d, _test_function_deriv_1d, weak_convolution,
    )

    kernel = _cached_kernel_1d(m, p)
    method = _derivatives.get_backend()

    if f.ndim == 3:
        Nx, Ny, Nz = f.shape
        if method == "spectral":
            KX, KY, KZ = _wavenumbers_3d(Nx, Ny, Nz, Lx, Ly, Lz)
            fh = np.fft.fftn(f)

        # 3-D multi-indices
        multi_indices = []
        for total in range(max_order + 1):
            for a in range(total + 1):
                for b in range(total - a + 1):
                    multi_indices.append((a, b, total - a - b))

        derivs, deriv_names = [], []
        for a, b, c in multi_indices:
            if a + b + c == 0:
                deriv = np.ones_like(f)
            elif method == "spectral":
                d_hat = (1j * KX) ** a * (1j * KY) ** b * (1j * KZ) ** c * fh
                deriv = np.real(np.fft.ifftn(d_hat))
            else:
                raise NotImplementedError("Polynomial backend not supported for 3-D.")
            name = ("x" * a) + ("y" * b) + ("z" * c) if (a + b + c) > 0 else "1"
            name = "f_" + name if (a + b + c) > 0 else "1"
            derivs.append(deriv)
            deriv_names.append(name)

        # Build columns with 3-D convolution
        cols, descs = [], []
        for poly_deg in range(max_poly + 1):
            for deriv, dname in zip(derivs, deriv_names):
                poly_field = f ** poly_deg if poly_deg > 0 else np.ones_like(f)
                col_raw = poly_field * deriv

                # IBP for pure linear terms: ∫ ψ · ∂_x^a∂_y^b∂_z^c f
                a_plus_bc = 0
                if dname != "1" and dname.startswith("f_"):
                    dx_part = dname[2:]
                    a_plus_bc = len(dx_part)

                if use_ibp and poly_deg == 0 and a_plus_bc > 0:
                    a_i = dx_part.count("x")
                    b_i = dx_part.count("y")
                    c_i = dx_part.count("z")
                    kx_ibp = _test_function_deriv_1d(m, p, a_i)
                    ky_ibp = _test_function_deriv_1d(m, p, b_i)
                    kz_ibp = _test_function_deriv_1d(m, p, c_i)
                    dx_a = (Lx / Nx) ** a_i
                    dy_b = (Ly / Ny) ** b_i
                    dz_c = (Lz / Nz) ** c_i
                    kx_scaled = kx_ibp / dx_a
                    ky_scaled = ky_ibp / dy_b
                    kz_scaled = kz_ibp / dz_c
                    col_weak = weak_convolution(f, kx_scaled, ky_scaled, kz_scaled)
                else:
                    col_weak = weak_convolution(col_raw, kernel, kernel, kernel)

                cols.append(col_weak.ravel())

                if poly_deg == 0:
                    desc = dname
                elif poly_deg == 1:
                    desc = "f·" + dname if dname != "1" else "f"
                else:
                    sup = _sup(poly_deg)
                    desc = f"f{sup}·{dname}" if dname != "1" else f"f{sup}"
                descs.append(desc)

        return np.column_stack(cols), descs

    # ===== 2-D path (original) =====
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
