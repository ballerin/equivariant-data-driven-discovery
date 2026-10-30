"""
Weak form in space, in time, in both, or in neither.

The per-snapshot library builders (``build_*_library(f, ..., formulation)``)
integrate against a test function in space only; the time derivative on the
left-hand side comes from a central finite difference.  This module adds the
temporal half of the Messenger & Bortz (arXiv:2007.02848) weak form, where
the test function is ψ(x, t) = φ_t(t) · ∏_d φ_x(x_d):

    b_k = ⟨ψ_k, ∂_t u⟩ = −⟨∂_t ψ_k, u⟩          (time IBP: no ∂_t of the data)
    G_kj = ⟨ψ_k, N_j[u]⟩

Every library column N_j[u] is instantaneous in time (spatial derivatives and
powers of u at one t), and ψ is separable, so the space-time convolution
factorises: build each snapshot's spatial column with the existing builder,
stack along t, and convolve along t.  No library builder changes.

The four formulations
---------------------
============== ============ ========== ==========================================
formulation    space        time       b
============== ============ ========== ==========================================
pointwise      pointwise    pointwise  f_t (central FD)
space_weak     weak         pointwise  φ_x ∗ f_t
time_weak      pointwise    weak       φ_t′ ∗ f            (or φ_t ∗ f_t)
spacetime_weak weak         weak       (φ_t′ ⊗ φ_x) ∗ f    (or φ_t ⊗ φ_x ∗ f_t)
============== ============ ========== ==========================================

The parenthesised variants are ``time_ibp=False``: the temporal test function
is still applied, but to the finite-difference f_t.

Unlike the per-snapshot sweep builders, query points share one spatial set
across all frames — required once a test function spans several frames.

Spatial boundary: ``boundary="valid"`` (M&B) keeps only centres whose support
fits in the grid; ``"periodic"`` wraps the spatial support around the domain
so every grid point is a centre (see ``weak_form.use_weak_boundary``).  Time
is always 'valid'.  Experiments from 2026-09 on use ``"periodic"`` (the
simulated domains are periodic; the reasons are in ``weak_form``); ``"valid"``
remains the default for reproducibility.

Sign convention: ``np.convolve`` flips its kernel, which supplies the
(−1)^|α| of integration by parts, so derivative kernels are used as-is (see
``wsindy-test-function-rule.md`` §9 and eq. (3.4) of the paper).
"""

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from .weak_form import (
    _test_function_1d,
    _test_function_deriv_1d,
    get_weak_boundary,
    use_weak_boundary,
    weak_convolution,
    weak_output_shape,
)


FORMULATIONS = {
    # name: (weak in space, weak in time)
    "pointwise": (False, False),
    "space_weak": (True, False),
    "time_weak": (False, True),
    "spacetime_weak": (True, True),
}


def time_kernel(m_t, p_t, order=0, dt=1.0):
    """1-D temporal test-function kernel ∂_t^order φ_t on 2m_t+1 frames,
    in physical units (scaled by dt^−order)."""
    if order == 0:
        return _test_function_1d(m_t, p_t)
    return _test_function_deriv_1d(m_t, p_t, order) / dt ** order


def time_weak_convolve(A, kernel):
    """Valid convolution of ``A`` along axis 0 (time) with a 1-D kernel.

    Matches ``np.convolve(A[:, i], kernel, mode="valid")`` for every trailing
    index i, vectorised over the trailing axes.

    Parameters
    ----------
    A : ndarray, shape (Nt, ...)
    kernel : ndarray, shape (2m_t + 1,)

    Returns
    -------
    ndarray, shape (Nt − 2m_t, ...)
    """
    kernel = np.ravel(kernel)
    if A.shape[0] < kernel.size:
        raise ValueError(
            f"time_weak_convolve: {A.shape[0]} frames cannot hold a "
            f"{kernel.size}-frame test function."
        )
    windows = sliding_window_view(A, kernel.size, axis=0)
    return windows @ kernel[::-1]


def spatial_query_indices(valid_shape, stride=1, point_mask=None):
    """Flat indices of the spatial query grid, shared by every frame.

    Parameters
    ----------
    valid_shape : tuple of int - spatial shape of one column (the full grid
        for pointwise-in-space or a periodic weak boundary, the valid
        interior ``N − 2m`` for a 'valid' weak boundary)
    stride : int - uniform subsampling along every spatial axis
    point_mask : bool ndarray or None, shape ``valid_shape`` - if given,
        grid points where it is False are dropped
    """
    axes = [np.arange(0, n, stride) for n in valid_shape]
    grid = np.meshgrid(*axes, indexing="ij")
    flat = np.ravel_multi_index([g.ravel() for g in grid], valid_shape)
    if point_mask is not None:
        flat = flat[np.asarray(point_mask).ravel()[flat]]
    return flat


def build_weak_system(snapshots, times, builder, formulation, *,
                      m=None, p=None, space_ibp=False,
                      m_t=None, p_t=None, time_ibp=True,
                      spatial_stride=1, time_stride=1, point_mask_fn=None,
                      boundary=None):
    """Assemble (G, b) for any of the four formulations.

    Parameters
    ----------
    snapshots : ndarray, shape (Nt, N_1, ..., N_ndim) - uniformly spaced in
        time, transient already removed
    times : ndarray, shape (Nt,)
    builder : callable ``builder(f, formulation, **kw) -> (Theta, desc)`` for
        a single snapshot ``f``; ``formulation`` is ``"pointwise"`` or
        ``"space_weak"`` and ``kw`` is ``m, p, space_ibp`` for space_weak.  The library
        facade functions fit this after binding their library arguments.
        The derivative backend is whatever is active when this is called.
    formulation : one of :data:`FORMULATIONS`
    m, p : int - spatial test function (weak in space only)
    space_ibp : bool - move spatial derivatives onto φ_x where the library
        supports it (the builders' ``space_ibp``)
    m_t, p_t : int - temporal test function (weak in time only)
    time_ibp : bool - b = φ_t′ ∗ f (True) or φ_t ∗ f_t with central-FD f_t
        (False).  Weak in time only.
    spatial_stride, time_stride : int - query-point subsampling (the paper's
        s_x, s_t)
    point_mask_fn : callable ``f -> bool array (f.size,)`` or None - keep
        only spatial points where the mask holds in *every* frame used.
        Pointwise in space only; ignored for weak-in-space, as in the
        existing space-weak sweep path.
    boundary : {'valid', 'periodic'} or None - spatial boundary of the weak
        convolution (None: the active ``weak_form`` mode, 'valid' by
        default).  Weak in space only.

    Returns
    -------
    G : ndarray, shape (K, n_cols)
    b : ndarray, shape (K,)
    desc : list of str
    info : dict - resolved parameters and query-grid sizes
    """
    if formulation not in FORMULATIONS:
        raise ValueError(f"Unknown formulation {formulation!r}; expected one "
                         f"of {sorted(FORMULATIONS)}.")
    weak_x = FORMULATIONS[formulation][0]
    with use_weak_boundary(boundary if weak_x else None):
        return _build_weak_system(
            snapshots, times, builder, formulation, m=m, p=p,
            space_ibp=space_ibp, m_t=m_t, p_t=p_t, time_ibp=time_ibp,
            spatial_stride=spatial_stride, time_stride=time_stride,
            point_mask_fn=point_mask_fn,
        )


def _build_weak_system(snapshots, times, builder, formulation, *, m, p,
                       space_ibp, m_t, p_t, time_ibp, spatial_stride,
                       time_stride, point_mask_fn):
    """:func:`build_weak_system` under the already-active spatial boundary."""
    if formulation not in FORMULATIONS:
        raise ValueError(f"Unknown formulation {formulation!r}; expected one "
                         f"of {sorted(FORMULATIONS)}.")
    weak_x, weak_t = FORMULATIONS[formulation]

    snapshots = np.asarray(snapshots, dtype=float)
    times = np.asarray(times, dtype=float)
    dts = np.diff(times)
    if not np.allclose(dts, dts[0], rtol=1e-6, atol=0.0):
        raise ValueError("build_weak_system needs uniformly spaced times.")
    dt = float(dts[0])
    Nt = snapshots.shape[0]
    spatial_shape = snapshots.shape[1:]
    ndim = len(spatial_shape)

    if weak_x and (m is None or p is None):
        raise ValueError(f"{formulation} needs spatial m and p.")
    if weak_t and (m_t is None or p_t is None):
        raise ValueError(f"{formulation} needs temporal m_t and p_t.")

    # Which frames enter G, and the per-frame left-hand-side field.
    # Differentiating in time (FD) loses one frame at each end.
    use_fd = (not weak_t) or (not time_ibp)
    if use_fd:
        frames = np.arange(1, Nt - 1)
        lhs_fields = (snapshots[2:] - snapshots[:-2]) / (2.0 * dt)
    else:
        frames = np.arange(Nt)
        lhs_fields = snapshots

    # Spatial query grid, fixed across frames
    if weak_x:
        k_x = _test_function_1d(m, p)
        valid_shape = weak_output_shape(spatial_shape, m)
        if min(valid_shape) < 1 or 2 * m + 1 > min(spatial_shape):
            raise ValueError(f"m={m} does not fit the spatial grid "
                             f"{spatial_shape}.")
        mask = None
    else:
        k_x = None
        valid_shape = spatial_shape
        mask = None
        if point_mask_fn is not None:
            mask = np.ones(int(np.prod(spatial_shape)), dtype=bool)
            for s in frames:
                mask &= np.asarray(point_mask_fn(snapshots[s])).ravel()
            mask = mask.reshape(spatial_shape)
    q_idx = spatial_query_indices(valid_shape, spatial_stride, mask)
    if q_idx.size == 0:
        raise ValueError("No spatial query points left after masking.")

    # Per-frame spatial columns, sliced to the query set before stacking
    G_frames, b_frames, desc = [], [], None
    for i, s in enumerate(frames):
        f = snapshots[s]
        if weak_x:
            Theta, d = builder(f, "space_weak", m=m, p=p, space_ibp=space_ibp)
            lhs = weak_convolution(lhs_fields[i], [k_x] * ndim).ravel()
        else:
            Theta, d = builder(f, "pointwise")
            lhs = lhs_fields[i].ravel()
        if desc is None:
            desc = d
        G_frames.append(Theta[q_idx])
        b_frames.append(lhs[q_idx])
    G_stack = np.stack(G_frames)          # (n_frames, n_q, n_cols)
    b_stack = np.stack(b_frames)          # (n_frames, n_q)

    if weak_t:
        if 2 * m_t + 1 > G_stack.shape[0]:
            raise ValueError(f"m_t={m_t} needs {2 * m_t + 1} frames; only "
                             f"{G_stack.shape[0]} available.")
        phi_t = time_kernel(m_t, p_t, 0, dt)
        G_stack = time_weak_convolve(G_stack, phi_t)
        b_kernel = time_kernel(m_t, p_t, 1, dt) if time_ibp else phi_t
        b_stack = time_weak_convolve(b_stack, b_kernel)

    t_idx = np.arange(0, G_stack.shape[0], time_stride)
    G = G_stack[t_idx].reshape(-1, G_stack.shape[-1])
    b = b_stack[t_idx].reshape(-1)

    info = dict(
        formulation=formulation, dt=dt,
        m=m if weak_x else None, p=p if weak_x else None,
        space_ibp=bool(space_ibp) if weak_x else None,
        boundary=get_weak_boundary() if weak_x else None,
        m_t=m_t if weak_t else None, p_t=p_t if weak_t else None,
        time_ibp=bool(time_ibp) if weak_t else None,
        n_time_query=int(t_idx.size), n_space_query=int(q_idx.size),
    )
    return G, b, desc, info
