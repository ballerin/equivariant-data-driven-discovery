"""
Universal evaluation of sparse PDE recovery against ground truth.

Provides:
  build_eval_theta    — pre-build the pointwise Theta from snapshots (slow;
                        do once per dataset).
  evaluate_from_theta — compute fp, fn, E, etc. from pre-built (Theta, Ft)
                        (fast; call for each hyperparameter combination).
  evaluate_recovery   — convenience wrapper: build + evaluate in one call.
  best_of_sweep       — select the best row from a sweep using the
                        lexicographic criterion (exact → fp → fn → E).

All evaluations use the **pointwise** design matrix Theta on all grid
points, regardless of how the coefficients were obtained (pointwise,
WSINDy, IBP, etc.).  This gives a universal metric space.

For sweep scripts, the recommended pattern is::

    Theta_eval, Ft_eval = build_eval_theta(snapshots, times, Lx, Ly)
    for each hyperparameter combination:
        w = regression(...)
        result = evaluate_from_theta(w, true_eq, Theta_eval, Ft_eval)
"""

import math
import numpy as np

from library import build_equivariant_library
from simulations import make_regression_data


def _is_nan(val):
    """Check if a value is NaN (float)."""
    try:
        return math.isnan(float(val))
    except (ValueError, TypeError):
        return False


def build_eval_theta(snapshots, times, Lx, Ly, skip_fraction=0.0, seed=0):
    """Pre-build the pointwise Theta and Ft for fast repeated evaluation.

    This is the expensive step — it constructs Theta on ALL grid points
    across ALL interior snapshots.  Call it **once** per dataset; then
    pass the result to ``evaluate_from_theta`` for each hyperparameter
    combination.

    Parameters
    ----------
    snapshots  : ndarray (n_saved, Nx, Ny)
    times      : ndarray (n_saved,)
    Lx, Ly     : float
    skip_fraction : float — fraction of front snapshots to discard
    seed       : int — unused (compatibility)

    Returns
    -------
    Theta : ndarray (n_total, K) — pointwise equivariant library on ALL grid
                                    points.
    Ft    : ndarray (n_total,)   — time derivative on ALL grid points.
    """
    f_snaps, ft_snaps, _, _ = make_regression_data(
        snapshots, times, n_per_snap=200, skip_fraction=skip_fraction, seed=seed,
    )
    rows_theta, rows_ft = [], []
    for s in range(len(f_snaps)):
        T, _ = build_equivariant_library(f_snaps[s], Lx, Ly)
        rows_theta.append(T)
        rows_ft.append(ft_snaps[s].ravel())
    return np.vstack(rows_theta), np.concatenate(rows_ft)


def evaluate_from_theta(w, true_eq, Theta, Ft):
    """Compute recovery metrics from pre-built (Theta, Ft).

    Parameters
    ----------
    w       : ndarray (K,) — recovered coefficient vector
    true_eq : ndarray (K,) — ground-truth coefficients
    Theta   : ndarray (n, K) — pointwise equivariant library
    Ft      : ndarray (n,)   — time derivative

    Returns
    -------
    dict with keys:
        fp, fn, exact_recovery, n_active, E, max_coeff_err,
        n_samples
    """
    active_mask = np.abs(true_eq) > 1e-10
    nz_mask = np.abs(w) > 1e-10

    fp = int(np.sum(nz_mask & ~active_mask))
    fn = int(np.sum(~nz_mask & active_mask))
    exact = int(fp == 0 and fn == 0)
    n_active = int(np.sum(nz_mask))

    # Coefficient error on true support
    if active_mask.any() and np.any(nz_mask & active_mask):
        max_coeff_err = float(np.max(np.abs(w[active_mask] - true_eq[active_mask])))
    elif active_mask.any():
        max_coeff_err = float(np.max(np.abs(true_eq[active_mask])))
    else:
        max_coeff_err = float("nan")

    # Prediction error E = ||Theta·(w - true)|| / ||Theta·true||
    delta = w - true_eq
    pred_true = Theta @ true_eq
    pred_delta = Theta @ delta
    denom = np.linalg.norm(pred_true) + 1e-30
    E = float(np.linalg.norm(pred_delta) / denom)

    return dict(
        fp=fp,
        fn=fn,
        exact_recovery=exact,
        n_active=n_active,
        E=E,
        max_coeff_err=max_coeff_err,
        n_samples=Theta.shape[0],
    )


def evaluate_standard(w, Theta, Ft):
    """Compute the irreducible fit error for a standard-library recovery.

    Since the standard (coordinate-basis) library cannot express the
    ground-truth equivariant operators (e.g. ‖∇f‖² = f_x² + f_y² is a
    product of two derivatives, not of the form f^p·∂^a f), we cannot
    compute fp/fn or E (distance to truth).  Instead we compute the
    **fit error**:

        E_fit = ||Theta · w − Ft||₂ / (||Ft||₂ + ε)

    which measures how well the recovered PDE reproduces the observed
    time derivative, regardless of interpretability.

    Parameters
    ----------
    w     : ndarray (K_std,) — recovered coefficient vector
    Theta : ndarray (n, K_std) — standard-library design matrix
    Ft    : ndarray (n,) — time-derivative vector

    Returns
    -------
    dict with keys:
        n_active, E_fit, n_samples
    """
    n_active = int(np.sum(np.abs(w) > 1e-10))
    pred = Theta @ w
    denom = np.linalg.norm(Ft) + 1e-30
    E_fit = float(np.linalg.norm(pred - Ft) / denom)
    return dict(
        n_active=n_active,
        E_fit=E_fit,
        n_samples=Theta.shape[0],
    )


# ── Equivariant → standard library decomposition ──────────────────────────

# Maps equivariant operator name → list of (standard_column_name, coefficient)
# for operators that can be exactly expressed in the standard (coordinate-basis)
# library f^p · ∂_x^a ∂_y^b f.  Operators containing products of derivatives
# (‖∇f‖², f·‖∇f‖², etc.) are NOT representable and are omitted here.
_EQUIV_TO_STD_DECOMPOSITION = {
    "f":       [("f", 1.0)],
    "Δf":      [("f_xx", 1.0), ("f_yy", 1.0)],
    "Δ²f":     [("f_xxxx", 1.0), ("f_xxyy", 2.0), ("f_yyyy", 1.0)],
    "f²":      [("f²", 1.0)],
    "f³":      [("f³", 1.0)],
    "f²·Δf":   [("f²·f_xx", 1.0), ("f²·f_yy", 1.0)],
}


def decompose_equivariant_to_standard(true_eq, desc_eq, desc_std):
    """Decompose equivariant ground truth into standard-library coefficients.

    For operators that can be expressed as f^p · ∂_x^a ∂_y^b f, this
    produces a coefficient vector in the standard basis.  Operators
    containing products of derivatives (‖∇f‖², f·‖∇f‖², etc.) are
    silently omitted since they have no representation in the standard
    library.

    Parameters
    ----------
    true_eq  : ndarray (K_eq,) — ground-truth coefficients in equivariant basis
    desc_eq  : list of K_eq str — equivariant operator descriptions
    desc_std : list of K_std str — standard column descriptions

    Returns
    -------
    true_std : ndarray (K_std,) — ground truth in standard basis
    """
    true_std = np.zeros(len(desc_std), dtype=float)
    for i, (name, coeff) in enumerate(zip(desc_eq, true_eq)):
        if abs(coeff) < 1e-12:
            continue
        if name not in _EQUIV_TO_STD_DECOMPOSITION:
            continue  # operator not representable in standard library
        for std_name, mult in _EQUIV_TO_STD_DECOMPOSITION[name]:
            try:
                j = desc_std.index(std_name)
            except ValueError:
                raise KeyError(
                    f"Standard column {std_name!r} not found in library. "
                    f"Available: {desc_std[:10]}..."
                )
            true_std[j] += coeff * mult
    return true_std


def evaluate_standard_with_truth(w, true_eq, desc_eq, Theta_std, desc_std, Ft):
    """Evaluate a standard-library recovery with proper fp/fn via decomposition.

    Decomposes the equivariant ground truth into the standard library
    wherever possible, then computes fp, fn, exact_recovery, E, and
    max_coeff_err in the standard basis.  Also computes E_fit.

    Parameters
    ----------
    w         : ndarray (K_std,)
    true_eq   : ndarray (K_eq,)
    desc_eq   : list of K_eq str
    Theta_std : ndarray (n, K_std)
    desc_std  : list of K_std str
    Ft        : ndarray (n,)

    Returns
    -------
    dict with keys: fp, fn, exact_recovery, n_active, E, max_coeff_err,
    E_fit, n_samples
    """
    true_std = decompose_equivariant_to_standard(true_eq, desc_eq, desc_std)
    return evaluate_from_theta(w, true_std, Theta_std, Ft)


def evaluate_recovery(w, true_eq, snapshots, times, Lx, Ly,
                      skip_fraction=0.0, seed=0):
    """Convenience: build Theta and evaluate in one call.

    For sweep scripts that evaluate many w vectors against the same data,
    prefer ``build_eval_theta`` + ``evaluate_from_theta`` to avoid
    rebuilding Theta for every hyperparameter combination.

    Parameters
    ----------
    w          : ndarray (K,)
    true_eq    : ndarray (K,)
    snapshots  : ndarray (n_saved, Nx, Ny)
    times      : ndarray (n_saved,)
    Lx, Ly     : float
    skip_fraction : float
    seed       : int — unused (compatibility)

    Returns
    -------
    dict with keys: fp, fn, exact_recovery, n_active, E, max_coeff_err,
    n_samples
    """
    Theta, Ft = build_eval_theta(snapshots, times, Lx, Ly, skip_fraction, seed)
    return evaluate_from_theta(w, true_eq, Theta, Ft)


def best_of_sweep(rows, criterion="lexicographic"):
    """Select the best row from a sweep according to a multi-level criterion.

    Parameters
    ----------
    rows : list of dict
        Each dict must have keys: 'exact_recovery', 'fp', 'fn', 'E'.
        (exact_recovery is a bool or int; others are numeric.)
    criterion : str
        "lexicographic" — compare exact_recovery first (True > False),
        then fp (lower better), then fn (lower better), then E (lower
        better).

    Returns
    -------
    best : dict — the row with the best score.
    """
    if criterion != "lexicographic":
        raise ValueError(f"Unknown criterion: {criterion}")

    def key(row):
        # exact_recovery: True=0 (best), False=1
        exact_score = 0 if row["exact_recovery"] else 1
        # NaN fp/fn → treat as infinity (worst possible), so standard-library
        # rows always lose to any row with finite metrics.
        fp_val = row["fp"] if not _is_nan(row["fp"]) else float("inf")
        fn_val = row["fn"] if not _is_nan(row["fn"]) else float("inf")
        return (exact_score, fp_val, fn_val, row["E"])

    return min(rows, key=key)
