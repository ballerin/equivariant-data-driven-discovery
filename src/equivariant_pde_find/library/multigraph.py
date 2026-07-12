"""
Pointwise multigraph library construction.

Internal module — use ``library.build_multigraph_library()`` (the public
facade) instead of calling ``_build_multigraph_pointwise`` directly.
"""

import numpy as np

from ..multigraph_operators import get_operators


def _build_multigraph_pointwise(f, Lx, Ly, max_p=2, max_n=3):
    """Build the equivariant library matrix for a 2-D scalar field f.

    Parameters
    ----------
    f     : ndarray (Nx, Ny) - data
    Lx    : float - domain length in x
    Ly    : float - domain length in y
    max_p : int - maximum multigraph edge count p (differential order = 2*max_p)
    max_n : int - maximum polynomial degree

    Returns
    -------
    Theta        : ndarray, shape (Nx*Ny, K)
    descriptions : list of K strings
    """
    ops = get_operators(max_p, max_n)
    cols, descs = [], []
    for op, name in ops:
        cols.append(op(f, Lx, Ly).ravel())
        descs.append(name)
    return np.column_stack(cols), descs


def _build_multigraph_weak(f, Lx, Ly, max_p=2, max_n=3, m=8, p=4,
                           use_ibp=False):
    """Build the weak-form equivariant (multigraph) library matrix.

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
        If True, use integration by parts for operators of the form
        Δ(something) or Δ²(something).  The Δ/Δ² is moved onto the
        analytically-known test function, eliminating |k|² and |k|⁴
        noise amplification.  Default False (direct weak-form).

    Returns
    -------
    Theta_weak : ndarray of shape (N_valid, K)
        Weak library matrix.  N_valid = (Nx - 2*m) * (Ny - 2*m).
    descriptions : list of K str
    """
    from ..multigraph_operators import get_operators
    from .weak_form import (
        _cached_kernel_1d, _is_laplacian_of,
        weak_laplacian_ibp, weak_bilaplacian_ibp, weak_convolution,
    )
    ops = get_operators(max_p, max_n)
    kernel = _cached_kernel_1d(m, p)

    cols = []
    for func, desc in ops:
        ibp_info = _is_laplacian_of(desc) if use_ibp else None
        if ibp_info is not None:
            op_type, inner_func = ibp_info
            inner_field = inner_func(f, Lx, Ly)
            if op_type == "laplacian":
                col_weak = weak_laplacian_ibp(inner_field, m, p, Lx, Ly)
            elif op_type == "bilaplacian":
                col_weak = weak_bilaplacian_ibp(inner_field, m, p, Lx, Ly)
            else:
                # Fallback
                col_raw = func(f, Lx, Ly)
                col_weak = weak_convolution(col_raw, kernel, kernel)
        else:
            # Direct path: compute operator on data, then convolve.
            col_raw = func(f, Lx, Ly)
            col_weak = weak_convolution(col_raw, kernel, kernel)
        cols.append(col_weak.ravel())

    descriptions = [desc for _, desc in ops]
    return np.column_stack(cols), descriptions
