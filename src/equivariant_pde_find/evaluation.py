"""
Evaluation of recovered PDEs against data and ground truth.

Provides a universal metric space for comparing PDE recovery across
different libraries (standard, WSINDy, multigraph, modified multigraph,
polynomial-jet) and formulations (pointwise, weak-form).

The idea is to map all library types to a common polynomial-jet basis,
and then evaluate a recovery in such spanning space.

Primary entry points
--------------------
build_eval_theta_multigraphs  - equivariant (multigraph) library on all
                                grid points and all interior snapshots.
build_eval_theta_jet          - polynomial-jet library (universal space).
evaluate_from_theta           - compute fp, fn, E_coeff, E_data from a
                                pre-built (Theta, Ft) pair.
evaluate_unified              - map any library's coefficients to jet space
                                and evaluate there.
best_of_sweep                 - lexicographic best-row selector

Usage pattern
-------------

    # Compute Theta, Ft
    Theta_eval, Ft_eval = build_eval_theta_jet(snapshots, times, Lx, Ly)

    # Recover PDE coefficients w and evaluate
    for each hyperparameter combination:
        w = regression(...)
        result = evaluate_unified(w, library, desc_lib, desc_jet,
                                  Theta_eval, Ft_eval, true_jet)
"""

import math
import numpy as np

from .library import build_multigraph_library, build_polyjet_library
from .library.conversions import (
    build_wsindy_to_jet_map,
    map_multigraph_to_jet,
    map_mod_multigraph_to_jet,
    map_standard_to_jet,
    map_wsindy_to_jet,
)
from .simulations import make_regression_data


def _is_nan(val):
    """Check if a value is NaN (float)."""
    try:
        return math.isnan(float(val))
    except (ValueError, TypeError):
        return False


def build_eval_theta_multigraphs(snapshots, times, Lx, Ly,
                                  skip_fraction=0.0, seed=0, backend="spectral"):
    """Build the pointwise equivariant (multigraph) Theta and Ft.

    Constructs the design matrix in the multigraph library on every grid
    point across all interior snapshots.  By default uses spectral
    (FFT-based) derivatives because they are exact on periodic domains.

    For cross-library comparison, prefer ``build_eval_theta_jet``,
    which uses the polynomial-jet universal space as the evaluation basis.

    Parameters
    ----------
    snapshots       : ndarray (n_saved, Nx, Ny)
    times           : ndarray (n_saved,)
    Lx, Ly          : float - domain lengths
    skip_fraction   : float - fraction of leading snapshots to discard
    seed            : int - RNG seed for spatial subsampling
    backend         : str - derivative backend ("spectral" or
                        "polynomial").  Default is "spectral".

    Returns
    -------
    Theta : ndarray (n_total, K) - multigraph Theta on ALL grid points
    Ft    : ndarray (n_total,)   - time derivative on ALL grid points
    """
    from . import derivatives as _deriv

    with _deriv.use_backend(backend):
        f_snaps, ft_snaps, _, _ = make_regression_data(
            snapshots, times, n_per_snap=200,
            skip_fraction=skip_fraction, seed=seed,
        )
        rows_theta, rows_ft = [], []
        for s in range(len(f_snaps)):
            T, _ = build_multigraph_library(f_snaps[s], Lx, Ly)
            rows_theta.append(T)
            rows_ft.append(ft_snaps[s].ravel())
        return np.vstack(rows_theta), np.concatenate(rows_ft)


def build_eval_theta_jet(snapshots, times, Lx, Ly,
                          skip_fraction=0.0, seed=0, backend="spectral",
                          jet_monomials=None):
    """Build the polynomial-jet Theta and Ft for unified evaluation.

    This is the **recommended** evaluation entry point for cross-library
    comparison as polynomial-jet library can represent all the terms of
    the other libraries (potentially as a linear combination of terms),
    making it a universal comparison space into which ALL library types
    can be analytically mapped.

    The evaluation Theta is always built with **spectral** (FFT-based)
    derivatives by default, because on periodic domains these are exact
    up to Nyquist and give the most faithful representation of the true
    PDE right-hand side.  The evaluation data should be clean (noise-free)
    for the same reason.

    Parameters
    ----------
    snapshots  : ndarray (n_saved, Nx, Ny)
    times      : ndarray (n_saved,)
    Lx, Ly     : float - domain lengths
    skip_fraction : float - fraction of leading snapshots to discard
    seed       : int - RNG seed for spatial subsampling
    backend    : str - derivative backend ("spectral" or "polynomial").
                    Default is "spectral".
    jet_monomials : list of str or None
        If provided, build exactly these monomials via the dynamic builder
        (``build_polyjet_from_descriptions``).  If None (default), use the
        legacy capped builder (``build_polyjet_library`` with
        max_order=4, max_degree=3).

    Returns
    -------
    Theta : ndarray (n_total, K_jet) - jet library on ALL grid points
    Ft    : ndarray (n_total,)       - time derivative on ALL grid points
    """
    from . import derivatives as _deriv
    from .library import build_polyjet_from_descriptions as _bpfj

    with _deriv.use_backend(backend):
        f_snaps, ft_snaps, _, _ = make_regression_data(
            snapshots, times, n_per_snap=200,
            skip_fraction=skip_fraction, seed=seed,
        )
        rows_theta, rows_ft = [], []
        for s in range(len(f_snaps)):
            if jet_monomials is not None:
                T, _ = _bpfj(f_snaps[s], Lx, Ly, jet_monomials)
            else:
                T, _ = build_polyjet_library(
                    f_snaps[s], Lx, Ly, max_order=4, max_degree=3,
                )
            rows_theta.append(T)
            rows_ft.append(ft_snaps[s].ravel())
        return np.vstack(rows_theta), np.concatenate(rows_ft)

# ---------------------------------------------------------------------------
# Evaluation functions
# ---------------------------------------------------------------------------

def evaluate_unified(w, library, desc_lib, desc_jet, Theta_jet, Ft,
                     true_jet, M_ws_to_jet=None, ndim=2):
    """Evaluate a recovery from any library against the polynomial-jet truth.

    Maps the recovered coefficients to the jet basis, then computes
    fp, fn, E in the common jet space.

    Parameters
    ----------
    ndim : int
        Spatial dimension of the data `w` was recovered from (default 2,
        matching the paper's 2-D domain).  Only consulted for the
        multigraph/mod_multigraph jet decomposition, which depends on ndim;
        `M_ws_to_jet` (WSINDy) must itself already be built for the correct
        ndim via ``build_wsindy_to_jet_map(..., ndim=ndim)``.
    """
    if library in ("multigraph",):
        w_jet = map_multigraph_to_jet(w, desc_lib, desc_jet, ndim=ndim)
    elif library == "mod_multigraph":
        w_jet = map_mod_multigraph_to_jet(w, desc_lib, desc_jet, ndim=ndim)
    elif library == "pdefind":
        w_jet = map_standard_to_jet(w, desc_lib, desc_jet)
    elif library == "wsindy" and M_ws_to_jet is not None:
        w_jet = map_wsindy_to_jet(w, M_ws_to_jet)
    elif library == "polyjet":
        # Name-based mapping: the regression polyjet may differ from the
        # evaluation polyjet (e.g., old capped 105 cols vs new dynamic ~65).
        # Map by column name, matching the logic of map_standard_to_jet.
        w_jet = np.zeros(len(desc_jet))
        for j, desc in enumerate(desc_lib):
            if abs(w[j]) < 1e-14:
                continue
            if desc in desc_jet:
                w_jet[desc_jet.index(desc)] += w[j]
    else:
        # Fallback: use library's own Theta
        return evaluate_from_theta(w, true_jet, Theta_jet, Ft)

    return evaluate_from_theta(w_jet, true_jet, Theta_jet, Ft)


def evaluate_from_theta(w, true_eq, Theta, Ft):
    """Compute recovery metrics from pre-built (Theta, Ft).

    Parameters
    ----------
    w       : ndarray (K,) - recovered coefficient vector
    true_eq : ndarray (K,) - ground-truth coefficients
    Theta   : ndarray (n, K) - pointwise equivariant library
    Ft      : ndarray (n,)   - time derivative

    Returns
    -------
    dict with keys:
        fp, fn, exact_recovery, n_active, E, max_coeff_err,
        n_samples
    """
    active_mask = np.abs(true_eq) > 1e-10
    nz_mask = np.abs(w) > 1e-10

    n_active = int(np.sum(nz_mask))

    # Guard: if the ground truth is all zeros, the true PDE is not
    # representable in this library.  fp/fn/E are ill-defined; fall
    # back to E_fit and report NaN for the structural metrics.
    if not active_mask.any():
        pred = Theta @ w
        E_fit = float(
            np.linalg.norm(pred - Ft) / (np.linalg.norm(Ft) + 1e-30)
        )
        return dict(
            fp=float("nan"), fn=float("nan"), exact_recovery=0,
            n_active=n_active, E=float("nan"),
            max_coeff_err=float("nan"),
            n_samples=Theta.shape[0],
            E_fit=E_fit,
        )

    fp = int(np.sum(nz_mask & ~active_mask))
    fn = int(np.sum(~nz_mask & active_mask))
    exact = int(fp == 0 and fn == 0)

    # Coefficient error on true support
    if np.any(nz_mask & active_mask):
        max_coeff_err = float(np.max(np.abs(w[active_mask] - true_eq[active_mask])))
    else:
        max_coeff_err = float(np.max(np.abs(true_eq[active_mask])))

    # Coefficient error E_coeff = ||Theta·(w - true)|| / ||Theta·true||
    delta = w - true_eq
    pred_true = Theta @ true_eq
    pred_delta = Theta @ delta
    denom = np.linalg.norm(pred_true) + 1e-30
    E_coeff = float(np.linalg.norm(pred_delta) / denom)

    # Data-fit error E_data = ||Theta·w - Ft|| / ||Ft||
    pred_w = Theta @ w
    E_data = float(
        np.linalg.norm(pred_w - Ft) / (np.linalg.norm(Ft) + 1e-30)
    )

    return dict(
        fp=fp,
        fn=fn,
        exact_recovery=exact,
        n_active=n_active,
        E_coeff=E_coeff,
        max_coeff_err=max_coeff_err,
        n_samples=Theta.shape[0],
        E_data=E_data,
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
    w     : ndarray (K_std,) - recovered coefficient vector
    Theta : ndarray (n, K_std) - standard-library design matrix
    Ft    : ndarray (n,) - time-derivative vector

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


# ── Evaluate-with-truth functions ────────────────────────────────────────
# These marshal library-specific ground-truth decompositions (now in the
# ``library.conversions`` module) and call ``evaluate_from_theta``.


def evaluate_jet_with_truth(w, true_eq, desc_eq, desc_jet, Theta_jet, Ft):
    """Evaluate a polynomial-jet recovery via decomposition.

    Decomposes the equivariant ground truth into the jet basis, then
    computes fp, fn, exact_recovery, E, and max_coeff_err.

    Parameters
    ----------
    w         : ndarray (K_jet,)
    true_eq   : ndarray (K_eq,)
    desc_eq   : list of K_eq str
    desc_jet  : list of K_jet str
    Theta_jet : ndarray (n, K_jet)
    Ft        : ndarray (n,)

    Returns
    -------
    dict with keys: fp, fn, exact_recovery, n_active, E, max_coeff_err,
    n_samples
    """
    true_jet = map_multigraph_to_jet(true_eq, desc_eq, desc_jet)
    return evaluate_from_theta(w, true_jet, Theta_jet, Ft)


def evaluate_recovery(w, true_eq, snapshots, times, Lx, Ly,
                      skip_fraction=0.0, seed=0, backend="spectral"):
    """Convenience: build Theta and evaluate in one call.

    For sweep scripts that evaluate many w vectors against the same data,
    prefer ``build_eval_theta_multigraphs`` + ``evaluate_from_theta`` to
    avoid rebuilding Theta for every hyperparameter combination.

    Parameters
    ----------
    w          : ndarray (K,)
    true_eq    : ndarray (K,)
    snapshots  : ndarray (n_saved, Nx, Ny)
    times      : ndarray (n_saved,)
    Lx, Ly     : float
    skip_fraction : float
    seed       : int — RNG seed for spatial subsampling
    backend    : str — derivative backend.  Default is "spectral".

    Returns
    -------
    dict with keys: fp, fn, exact_recovery, n_active, E_coeff,
    max_coeff_err, n_samples
    """
    Theta, Ft = build_eval_theta_multigraphs(
        snapshots, times, Lx, Ly,
        skip_fraction=skip_fraction, seed=seed, backend=backend,
    )
    return evaluate_from_theta(w, true_eq, Theta, Ft)


def best_of_sweep(rows, criterion="lexicographic"):
    """Select the best row from a sweep according to a multi-level criterion.

    Parameters
    ----------
    rows : list of dict
        Each dict must have keys: 'exact_recovery', 'fp', 'fn', 'E'.
        (exact_recovery is a bool or int; others are numeric.)
    criterion : str
        "lexicographic" - compare exact_recovery first (True > False),
        then fp (lower better), then fn (lower better), then E (lower
        better).

    Returns
    -------
    best : dict - the row with the best score.
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
        return (exact_score, fp_val, fn_val, row["E_coeff"])

    return min(rows, key=key)
