"""
Sparse regression algorithms for PDE discovery.

Provides:
  STRidge       — sequential threshold ridge regression (inner loop)
  TrainSTRidge  — STRidge with cross-validated sparsity tolerance
"""

import numpy as np


def STRidge(X0, y, lam, maxit, tol, normalize=2):
    """Sequential threshold ridge regression (inner loop).

    Parameters
    ----------
    X0        : ndarray, shape (n, d) — design matrix
    y         : ndarray, shape (n,)   — target
    lam       : float — ridge regularisation parameter
    maxit     : int   — maximum thresholding iterations
    tol       : float — hard-threshold level (entries |w_i| < tol set to zero)
    normalize : int   — norm order for column normalisation (0 = none)

    Returns
    -------
    w : ndarray, shape (d,) — sparse coefficient vector
    """
    n, d = X0.shape
    y = y.ravel()

    if normalize != 0:
        Mreg = np.array([1.0 / (np.linalg.norm(X0[:, i], normalize) + 1e-14) for i in range(d)])
        X = X0 * Mreg  # broadcast (n, d) * (d,)
    else:
        X = X0.astype(float, copy=True)
        Mreg = np.ones(d)

    # Precompute normal equations once — O(n d²) — to avoid repeating in the
    # inner thresholding loop. Inner solves become cheap d×d problems O(d³).
    XtX = X.T @ X   # (d, d)
    Xty = X.T @ y   # (d,)

    # Initial least-squares solve (handles rank-deficient startup correctly)
    w = np.linalg.lstsq(X, y, rcond=None)[0]
    for _ in range(maxit):
        small = np.abs(w) < tol
        w[small] = 0.0
        big = ~small
        if big.sum() == 0:
            break
        # Solve normal equations on the active submatrix: O(d_active³) per step
        idx = np.where(big)[0]
        G = XtX[np.ix_(idx, idx)]
        h = Xty[idx]
        try:
            w[big] = np.linalg.solve(G, h)
        except np.linalg.LinAlgError:
            # Fallback for rank-deficient active submatrix
            w[big] = np.linalg.lstsq(X[:, big], y, rcond=None)[0]
    return Mreg * w


def TrainSTRidge(
    R, ut, lam, d_tol, maxit=25, STR_iters=10, l0_penalty=None, normalize=2, split=0.8
):
    """Train STRidge with cross-validation over the sparsity tolerance.

    Parameters
    ----------
    R          : ndarray, shape (M, K) — library matrix
    ut         : ndarray, shape (M,)   — time-derivative vector
    lam        : float — ridge regularisation parameter
    d_tol      : float — initial tolerance step for cross-validation
    maxit      : int   — number of tolerance search steps
    STR_iters  : int   — inner STRidge iterations per step
    l0_penalty : float or None — L0 penalty weight; defaults to 0.001*cond(R)
    normalize  : int   — column normalisation norm order (passed to STRidge)
    split      : float — fraction of data used for training

    Returns
    -------
    w : ndarray, shape (K,) — best sparse coefficient vector
    """
    np.random.seed(0)
    n = R.shape[0]
    ut = ut.ravel()
    train = np.random.choice(n, int(n * split), replace=False)
    # Use a boolean mask: O(n) instead of O(n * |train|).
    train_mask = np.zeros(n, dtype=bool)
    train_mask[train] = True
    test = np.where(~train_mask)[0]
    Rtr, Rte = R[train], R[test]
    ytr, yte = ut[train], ut[test]

    if l0_penalty is None:
        l0_penalty = 0.001 * np.linalg.cond(R)

    w_best = np.linalg.lstsq(Rtr, ytr, rcond=None)[0]
    err_best = np.linalg.norm(yte - Rte @ w_best) + l0_penalty * np.count_nonzero(w_best)
    tol = float(d_tol)

    for _ in range(maxit):
        w = STRidge(Rtr, ytr, lam, STR_iters, tol, normalize)
        err = np.linalg.norm(yte - Rte @ w) + l0_penalty * np.count_nonzero(w)
        if err <= err_best:
            err_best, w_best = err, w
            tol += d_tol
        else:
            d_tol = 2 * d_tol / (maxit + 1)
            tol = max(0, tol - 2 * d_tol) + d_tol
    return w_best


def SR3(X0, y, lam_l1, kappa=1.0, maxit=100, normalize=2):
    """Sparse Relaxed Regularized Regression (SR3).

    Solves the relaxed problem::

        min_{xi, w}  ½‖X0 @ xi − y‖² + (κ/2)‖xi − w‖² + λ‖w‖₁

    via block-coordinate descent:
        xi ← (X₀ᵀX₀ + κI)⁻¹ (X₀ᵀy + κ w)   [linear solve, O(d³) once]
        w  ← soft_threshold(xi, λ/κ)           [O(d)]

    SR3 is significantly more robust than STRidge for ill-conditioned design
    matrices because κI regularises the normal equations independently of the
    data.

    Parameters
    ----------
    X0     : ndarray, shape (n, d) — design matrix
    y      : ndarray, shape (n,)   — target
    lam_l1 : float — L₁ penalty λ (sparsity)
    kappa  : float — relaxation strength κ (larger → tighter ξ≈w coupling)
    maxit  : int   — coordinate-descent iterations
    normalize : int — column normalisation norm order (0 = none)

    Returns
    -------
    w : ndarray, shape (d,) — sparse coefficient vector
    """
    n, d = X0.shape
    y = y.ravel()

    if normalize != 0:
        Mreg = np.array([1.0 / (np.linalg.norm(X0[:, i], normalize) + 1e-14) for i in range(d)])
        X = X0 * Mreg[np.newaxis, :]
    else:
        X = X0.copy()
        Mreg = np.ones(d)

    XtX = X.T @ X                   # (d, d)
    Xty = X.T @ y                   # (d,)
    A = XtX + kappa * np.eye(d)     # precompute once
    A_inv = np.linalg.inv(A)        # cheap: d is small (20–60)

    threshold = lam_l1 / (kappa + 1e-30)

    w = np.zeros(d)
    xi = A_inv @ Xty                # warm start with w = 0

    for _ in range(maxit):
        xi = A_inv @ (Xty + kappa * w)
        w = np.sign(xi) * np.maximum(np.abs(xi) - threshold, 0.0)

    return Mreg * w


def TrainSR3(
    R, ut, lam, d_tol, maxit=25, SR3_iters=100, normalize=2, split=0.8
):
    """Train SR3 with cross-validation over the L₁ penalty.

    The API mirrors TrainSTRidge for drop-in comparison:
        lam   → κ (relaxation strength)
        d_tol → λ_L1 initial value (cross-validated upward, same logic as
                 TrainSTRidge's tolerance search)

    Parameters
    ----------
    R     : ndarray, shape (M, K) — library matrix
    ut    : ndarray, shape (M,)   — time-derivative vector
    lam   : float — κ, the SR3 relaxation parameter
    d_tol : float — initial L₁ penalty λ_L1
    maxit     : int — number of cross-validation steps
    SR3_iters : int — inner SR3 coordinate-descent iterations per step
    normalize : int — column normalisation (passed to SR3)
    split     : float — train/test split fraction

    Returns
    -------
    w : ndarray, shape (K,) — best sparse coefficient vector
    """
    np.random.seed(0)
    n = R.shape[0]
    ut = ut.ravel()
    train = np.random.choice(n, int(n * split), replace=False)
    train_mask = np.zeros(n, dtype=bool)
    train_mask[train] = True
    test = np.where(~train_mask)[0]
    Rtr, Rte = R[train], R[test]
    ytr, yte = ut[train], ut[test]

    l0_penalty = 0.001 * np.linalg.cond(R)

    w_best = np.linalg.lstsq(Rtr, ytr, rcond=None)[0]
    err_best = np.linalg.norm(yte - Rte @ w_best) + l0_penalty * np.count_nonzero(w_best)
    lam_l1 = float(d_tol)

    for _ in range(maxit):
        w = SR3(Rtr, ytr, lam_l1, kappa=lam, maxit=SR3_iters)
        err = np.linalg.norm(yte - Rte @ w) + l0_penalty * np.count_nonzero(w)
        if err <= err_best:
            err_best, w_best = err, w
            lam_l1 += d_tol
        else:
            d_tol = 2 * d_tol / (maxit + 1)
            lam_l1 = max(0, lam_l1 - 2 * d_tol) + d_tol
    return w_best
