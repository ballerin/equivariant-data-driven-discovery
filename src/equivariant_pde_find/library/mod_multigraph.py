"""
Pointwise modified-multigraph (IBP-friendly) library construction.

Internal module — use ``library.build_mod_multigraph_library()`` (the public
facade) instead of calling ``_build_mod_multigraph_pointwise`` directly.
"""

import numpy as np


def _build_mod_multigraph_pointwise(f, Lx, Ly, max_p=2, max_n=3):
    """Build the modified multigraph (IBP-friendly) library matrix.

    Uses the operator basis from Lemma 1 / Table 4 of the paper, which
    replaces 6 operators of the form f·(derivative) with Δ(something)
    or Δ²(something) forms.  This basis spans the same 20-dimensional
    space as ``build_multigraph_library`` but is designed for use with
    the weak-form integration-by-parts mechanism.

    Parameters
    ----------
    f     : ndarray, shape (Nx, Ny)
    Lx, Ly : float
    max_p : int - maximum edge count (0, 1, or 2)
    max_n : int - maximum polynomial degree (0..3)

    Returns
    -------
    Theta        : ndarray, shape (Nx*Ny, K)
    descriptions : list of K strings
    """
    from ..multigraph_operators import get_mod_multigraph_operators
    ops = get_mod_multigraph_operators(max_p, max_n)
    cols, descs = [], []
    for op, name in ops:
        cols.append(op(f, Lx, Ly).ravel())
        descs.append(name)
    return np.column_stack(cols), descs


def _build_mod_multigraph_weak(f, Lx, Ly, max_p=2, max_n=3,
                               m=8, p=4, use_ibp=False):
    """Build the weak-form **modified multigraph** library matrix.

    Like ``_build_multigraph_weak`` but uses the IBP-friendly operator basis
    from Lemma 1 / Table 4.  When ``use_ibp=True``, up to 8 operators
    (Δf, Δf², Δf³, Δ‖∇f‖², Δ(f·‖∇f‖²), Δ²f, Δ²f², Δ²f³) are evaluated
    via integration by parts instead of by differentiating noisy data.

    Parameters
    ----------
    f : ndarray of shape (Nx, Ny)
    Lx, Ly : float
    max_p, max_n : int
    m, p : int — test function parameters
    use_ibp : bool

    Returns
    -------
    Theta_weak : ndarray of shape (N_valid, K)
    descriptions : list of K str
    """
    from ..multigraph_operators import get_mod_multigraph_operators
    from .weak_form import (
        _cached_kernel_1d, _is_laplacian_of,
        weak_laplacian_ibp, weak_bilaplacian_ibp, weak_convolution,
    )
    ops = get_mod_multigraph_operators(max_p, max_n)
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
                col_raw = func(f, Lx, Ly)
                col_weak = weak_convolution(col_raw, kernel, kernel)
        else:
            col_raw = func(f, Lx, Ly)
            col_weak = weak_convolution(col_raw, kernel, kernel)
        cols.append(col_weak.ravel())

    descriptions = [desc for _, desc in ops]
    return np.column_stack(cols), descriptions
