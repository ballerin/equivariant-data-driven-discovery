"""
Weak-form regression data pipeline.

Sibling to ``regression_data.py``: takes snapshots and produces (G, b)
matrices ready for sparse regression in weak-form.
"""

import numpy as np

from ..library.weak_form import (
    _cached_kernel_1d,
    build_weak_rhs,
    compute_rescaling_scales,
    auto_test_function_params,
)
from ..multigraph_operators import get_operator_metadata
from .regression_data import make_regression_data


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
    space_ibp=False,
    rescale=False,
):
    """Build weak-form regression matrices from a stack of snapshots.

    This mirrors ``make_regression_data`` in ``simulations.py``:
    for each snapshot, the weak library and weak RHS are computed, and
    query points are drawn from the valid interior.

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
    space_ibp : bool
        Forwarded to ``multigraph._build_multigraph_space_weak``.
    rescale : bool
        If True, compute scale-invariance factors and apply column/LHS
        rescaling.  The returned ``Theta_all`` and ``f_t_all`` are scaled
        so that the solved coefficients must be unscaled via
        ``rescale_coefficients``.

    Returns
    -------
    Theta_all : ndarray of shape (n_total, K)
    f_t_all : ndarray of shape (n_total,)
    descriptions : list of K str
    rescale_info : dict or None
        ``{"col_scale": ..., "lhs_scale": ..., ...}`` if
        ``rescale=True``, else ``None``.
    """
    from ..library.multigraph import _build_multigraph_space_weak

    kernel = _cached_kernel_1d(m, p)

    Theta_parts, rhs_parts = [], []
    descriptions = None

    for snap_idx in range(len(f_snaps)):
        f = f_snaps[snap_idx]
        Theta_w, descs = _build_multigraph_space_weak(
            f, Lx, Ly, max_p, max_n, m=m, p=p, space_ibp=space_ibp,
        )
        if descriptions is None:
            descriptions = descs

        b_w = build_weak_rhs(ft_snaps[snap_idx], kernel_x=kernel, kernel_y=kernel)

        spatial_shape = f.shape
        Nv = tuple(n - 2 * m for n in spatial_shape)

        if sample_stride > 1:
            if f.ndim == 3:
                ix = np.arange(0, Nv[0], sample_stride)
                iy = np.arange(0, Nv[1], sample_stride)
                iz = np.arange(0, Nv[2], sample_stride)
                IX, IY, IZ = np.meshgrid(ix, iy, iz, indexing="ij")
                flat_idx = (IX * Nv[1] * Nv[2] + IY * Nv[2] + IZ).ravel()
            else:
                Nv_x, Nv_y = Nv
                ix = np.arange(0, Nv_x, sample_stride)
                iy = np.arange(0, Nv_y, sample_stride)
                IX, IY = np.meshgrid(ix, iy, indexing="ij")
                flat_idx = (IX * Nv_y + IY).ravel()
            Theta_parts.append(Theta_w[flat_idx, :])
            rhs_parts.append(b_w[flat_idx])
        else:
            idx = idx_sets[snap_idx]
            Theta_parts.append(Theta_w[idx, :])
            rhs_parts.append(b_w[idx])

    Theta_all = np.vstack(Theta_parts)
    f_t_all = np.concatenate(rhs_parts)

    rescale_info = None
    if rescale:
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
        Theta_all = Theta_all * col_scale[np.newaxis, :]
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

    Private internal function

    Returns
    -------
    col_scale : ndarray, shape (K,)
        Column scaling factors: G_scaled = G * diag(col_scale).
    lhs_scale : float
        LHS scaling factor: b_scaled = lhs_scale * b.
    """
    meta = get_operator_metadata(max_p, max_n)
    col_scale = np.empty(len(meta))
    for j, (p, k) in enumerate(meta):
        col_scale[j] = gamma_u ** (-k) * gamma_x ** (-2 * p)
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
# Full pipeline from simulation output, for convenience
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
    space_ibp=False,
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
    space_ibp : bool
        Forwarded to ``multigraph._build_multigraph_space_weak``.
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
    if auto_params:
        dx = Lx / snapshots.shape[1]
        dt_est = float(np.mean(np.diff(times))) if len(times) > 1 else 1.0
        m_x, p_x, m_t, p_t = auto_test_function_params(
            snapshots, dx, dt_est, max_deriv_spatial=2 * max_p,
        )
        m, p = m_x, p_x

    f_snaps, ft_snaps, idx_sets, _dt = make_regression_data(
        snapshots, times, n_per_snap=n_per_snap,
        skip_fraction=skip_fraction, seed=seed,
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
        space_ibp=space_ibp,
        rescale=rescale,
    )
