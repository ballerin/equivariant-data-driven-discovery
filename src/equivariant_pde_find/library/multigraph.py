"""
Pointwise multigraph library construction.

Internal module — use ``library.build_multigraph_library()`` (the public
facade) instead of calling ``_build_multigraph_pointwise`` directly.
"""

import numpy as np

from ..multigraph_operators import get_operators


def _build_multigraph_pointwise(f, Lx, Ly, max_p=2, max_n=3, Lz=None):
    """Build the equivariant library matrix for a scalar field f.

    Supports 2-D (Nx, Ny) and 3-D (Nx, Ny, Nz) fields.

    Parameters
    ----------
    f     : ndarray (Nx, Ny) or (Nx, Ny, Nz)
    Lx, Ly, Lz : float — domain lengths (Lz required for 3-D)
    max_p : int
    max_n : int

    Returns
    -------
    Theta        : ndarray, shape (N_total, K)
    descriptions : list of K strings
    """
    ops = get_operators(max_p, max_n)
    cols, descs = [], []
    for op, name in ops:
        cols.append(op(f, Lx, Ly, Lz).ravel())
        descs.append(name)
    return np.column_stack(cols), descs


def _build_multigraph_weak(f, Lx, Ly, max_p=2, max_n=3, m=8, p=4,
                           use_ibp=False, Lz=None):
    """Build the weak-form equivariant (multigraph) library matrix.

    Supports 2-D and 3-D fields.

    Parameters
    ----------
    f : ndarray (Nx, Ny) or (Nx, Ny, Nz)
    Lx, Ly, Lz : float
    max_p, max_n : int
    m : int — test function half-support
    p : int — test function polynomial degree
    use_ibp : bool

    Returns
    -------
    Theta_weak : ndarray, shape (N_valid, K)
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
            inner_field = inner_func(f, Lx, Ly, Lz)
            if op_type == "laplacian":
                col_weak = weak_laplacian_ibp(inner_field, m, p, Lx, Ly, Lz)
            elif op_type == "bilaplacian":
                col_weak = weak_bilaplacian_ibp(inner_field, m, p, Lx, Ly, Lz)
            else:
                col_raw = func(f, Lx, Ly, Lz)
                col_weak = weak_convolution(col_raw, kernel, kernel, kernel if f.ndim == 3 else None)
        else:
            col_raw = func(f, Lx, Ly, Lz)
            col_weak = weak_convolution(col_raw, kernel, kernel, kernel if f.ndim == 3 else None)
        cols.append(col_weak.ravel())

    descriptions = [desc for _, desc in ops]
    return np.column_stack(cols), descriptions
