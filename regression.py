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
    X = np.zeros_like(X0, dtype=float)
    y = y.ravel()

    if normalize != 0:
        Mreg = np.array([1.0 / (np.linalg.norm(X0[:, i], normalize) + 1e-14) for i in range(d)])
        for i in range(d):
            X[:, i] = Mreg[i] * X0[:, i]
    else:
        X = X0.copy()
        Mreg = np.ones(d)

    w = np.linalg.lstsq(X, y, rcond=None)[0]
    for _ in range(maxit):
        small = np.abs(w) < tol
        w[small] = 0.0
        big = ~small
        if big.sum() == 0:
            break
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
    test = np.array([i for i in range(n) if i not in train])
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
