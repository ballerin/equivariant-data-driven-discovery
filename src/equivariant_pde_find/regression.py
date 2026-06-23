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
        G = XtX[np.ix_(idx, idx)] + lam * np.eye(len(idx))
        h = Xty[idx]
        try:
            w[big] = np.linalg.solve(G, h)
        except np.linalg.LinAlgError:
            # Fallback for rank-deficient active submatrix
            w[big] = np.linalg.lstsq(X[:, big], y, rcond=None)[0]
    return Mreg * w


def TrainSTRidge(
    R, ut, lam, d_tol, maxit=25, STR_iters=10, l0_penalty=None,
    normalize=2, split=0.8, seed=0,
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
    l0_penalty : float or None — L0 penalty weight; defaults to 0.001*cond(Rtr)
    normalize  : int   — column normalisation norm order (passed to STRidge)
    split      : float — fraction of data used for training
    seed       : int   — RNG seed for reproducible train/test splits

    Returns
    -------
    w : ndarray, shape (K,) — best sparse coefficient vector
    """
    rng = np.random.default_rng(seed)
    n = R.shape[0]
    ut = ut.ravel()
    train = rng.choice(n, int(n * split), replace=False)
    # Use a boolean mask: O(n) instead of O(n * |train|).
    train_mask = np.zeros(n, dtype=bool)
    train_mask[train] = True
    test = np.where(~train_mask)[0]
    Rtr, Rte = R[train], R[test]
    ytr, yte = ut[train], ut[test]

    if l0_penalty is None:
        l0_penalty = 0.001 * np.linalg.cond(Rtr)

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


def SR3(X0, y, lam_l1, kappa=1.0, maxit=100, normalize=2, tol=1e-8):
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
    tol    : float — convergence tolerance on coefficient change

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
    # Use Cholesky-based solve for SPD systems (numerically superior to inv)
    try:
        L = np.linalg.cholesky(A)
        solve_A = lambda rhs: np.linalg.solve(
            L.T, np.linalg.solve(L, rhs)
        )
    except np.linalg.LinAlgError:
        solve_A = lambda rhs: np.linalg.solve(A, rhs)

    threshold = lam_l1 / (kappa + 1e-30)

    w = np.zeros(d)
    xi = solve_A(Xty)                # warm start with w = 0

    for _ in range(maxit):
        xi_new = solve_A(Xty + kappa * w)
        w = np.sign(xi_new) * np.maximum(np.abs(xi_new) - threshold, 0.0)
        if np.max(np.abs(xi_new - xi)) < tol:
            xi = xi_new
            break
        xi = xi_new

    return Mreg * w


def TrainSR3(
    R, ut, lam, d_tol, maxit=25, SR3_iters=100, normalize=2, split=0.8, seed=0,
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
    seed      : int — RNG seed for reproducible train/test splits

    Returns
    -------
    w : ndarray, shape (K,) — best sparse coefficient vector
    """
    rng = np.random.default_rng(seed)
    n = R.shape[0]
    ut = ut.ravel()
    train = rng.choice(n, int(n * split), replace=False)
    train_mask = np.zeros(n, dtype=bool)
    train_mask[train] = True
    test = np.where(~train_mask)[0]
    Rtr, Rte = R[train], R[test]
    ytr, yte = ut[train], ut[test]

    l0_penalty = 0.001 * np.linalg.cond(Rtr)

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


# ======================================================================
# MSTLS — Modified Sequential Thresholding Least Squares
#
# Messenger & Bortz (2021), "Weak SINDy for Partial Differential
# Equations", J. Comput. Phys., Section 3.2.
#
# When column scaling (M_scale) is not provided, the algorithm falls back
# to simple STRidge-style thresholding |w_j| < λ, matching the MATLAB
# ``sparsifyDynamics`` fallback path.
# ======================================================================


def mstls(G, b, lam, maxits=50, normalize=2):
    """Modified Sequential Thresholding Least-Squares (single λ).

    Implements the sequential thresholding from the WSINDy paper
    (Messenger & Bortz, 2021).  When *normalize* is non-zero, columns
    are scaled to unit norm and simple STRidge-style thresholding
    ``|w_j| < λ`` is used (matching the MATLAB ``sparsifyDynamics``
    fallback when ``M_scale`` is not provided).  The dominant-balance
    bounds approach requires per-column scaling (``M_scale``) which is
    not computed here.

    Parameters
    ----------
    G : ndarray, shape (n, d) — design matrix
    b : ndarray, shape (n,)   — target vector
    lam : float — threshold λ > 0; coefficients with |w_j| < λ are zeroed
    maxits : int — maximum thresholding iterations
    normalize : int — column normalisation norm order (0 = none)

    Returns
    -------
    w : ndarray, shape (d,) — sparse coefficient vector
    """
    n, d = G.shape
    y = b.ravel()

    if normalize != 0:
        Mreg = np.array(
            [1.0 / (np.linalg.norm(G[:, i], normalize) + 1e-14) for i in range(d)]
        )
        X = G * Mreg[np.newaxis, :]
    else:
        X = G.copy()
        Mreg = np.ones(d)

    w = np.linalg.lstsq(X, y, rcond=None)[0]
    active_prev = None
    for _ in range(maxits):
        # STRidge-style thresholding (matches MATLAB sparsifyDynamics
        # fallback when M_scale is not provided)
        active = np.abs(w) >= lam
        if active_prev is not None and np.array_equal(active, active_prev):
            break
        active_prev = active.copy()
        w[~active] = 0.0
        if np.sum(active) == 0:
            break
        w[active] = np.linalg.lstsq(X[:, active], y, rcond=None)[0]
    return Mreg * w


def select_lambda_grid(G, b, n_lambdas=50):
    """Auto-generate a log-spaced grid of candidate λ values for MSTLS.

    Based on the MATLAB ``wsindy_pde_RGLS_seq.m``:

    .. math::

        \\lambda_{\\max} = \\max_j |G_j^T b| / \\|G_j\\|^2
        \\lambda_{\\min} = \\|G w_{LS}\\| / (d \\cdot \\max_j \\|G_j\\|)

    Parameters
    ----------
    G : ndarray, shape (n, d)
    b : ndarray, shape (n,)
    n_lambdas : int — number of λ candidates

    Returns
    -------
    lambdas : ndarray, shape (n_lambdas,)
        Log-spaced from λ_min to λ_max.
    """
    _, d = G.shape
    y = b.ravel()
    col_norms = np.sqrt(np.sum(G ** 2, axis=0))
    # λ_max = largest regression coefficient magnitude
    Gty = np.abs(G.T @ y)
    lam_max = np.max(Gty / np.maximum(col_norms ** 2, 1e-30))

    # λ_min from projection magnitude
    w_ls = np.linalg.lstsq(G, y, rcond=None)[0]
    proj_norm = np.linalg.norm(G @ w_ls)
    lam_min = proj_norm / (d * max(np.max(col_norms), 1e-30))

    # Clamp to sensible range
    lam_min = max(lam_min, 1e-10)
    lam_max = min(lam_max, 1e4)
    if lam_min >= lam_max:
        lam_min = lam_max * 1e-4

    return np.logspace(np.log10(lam_min), np.log10(lam_max), n_lambdas)


def mstls_sweep(G, b, lambdas=None, n_lambdas=50, maxits=50, alpha=0.5, normalize=2):
    """MSTLS with automatic λ selection via the loss function.

    For each λ in the candidate grid, runs ``mstls`` and evaluates

    .. math::

        \\mathcal{L}(\\lambda) =
        \\alpha \\frac{\\|G(w^\\lambda - w^{LS})\\|}{\\|G w^{LS}\\|}
        + (1-\\alpha) \\frac{\\#\\{j : w^\\lambda_j \\ne 0\\}}{d}.

    Returns the coefficient vector with the *smallest* λ that minimises
    :math:`\\mathcal{L}`.  Ties are broken by choosing the sparser model.

    Parameters
    ----------
    G : ndarray, shape (n, d)
    b : ndarray, shape (n,)
    lambdas : ndarray or None — candidate λ values; auto-generated if None
    n_lambdas : int — used only when ``lambdas is None``
    maxits : int — inner MSTLS iterations per λ
    alpha : float — projection-cost weight in [0, 1]
    normalize : int — column normalisation (passed to ``mstls``)

    Returns
    -------
    w_best : ndarray, shape (d,)
    lambda_best : float
    loss_vals : ndarray, shape (n_lambdas,)
    """
    _, d = G.shape
    y = b.ravel()

    # Normalize columns (matching what mstls does internally) so that
    # the lambda grid and loss computations use the same scale.
    if normalize != 0:
        Mreg = np.array(
            [1.0 / (np.linalg.norm(G[:, i], normalize) + 1e-14) for i in range(d)]
        )
        G_norm = G * Mreg[np.newaxis, :]
    else:
        G_norm = G.copy()
        Mreg = np.ones(d)

    if lambdas is None:
        lambdas = select_lambda_grid(G_norm, y, n_lambdas)

    w_ls = np.linalg.lstsq(G_norm, y, rcond=None)[0]
    G_w_ls_norm = np.linalg.norm(G_norm @ w_ls)
    if G_w_ls_norm < 1e-30:
        G_w_ls_norm = 1.0

    best_loss = np.inf
    best_w = np.zeros(d)  # unscaled
    best_lambda = lambdas[0]
    loss_vals = np.empty(len(lambdas))

    for i, lam in enumerate(lambdas):
        w = mstls(G, b, lam, maxits=maxits, normalize=normalize)
        # Loss computed on NORMALIZED system for consistency with lambda grid
        w_norm = w / np.maximum(Mreg, 1e-30)
        proj_cost = np.linalg.norm(G_norm @ (w_norm - w_ls)) / G_w_ls_norm
        nnz = np.count_nonzero(w)
        sparsity_cost = nnz / max(d, 1)
        loss = alpha * proj_cost + (1.0 - alpha) * sparsity_cost
        loss_vals[i] = loss

        # Keep the *smallest* λ that achieves the minimum loss.
        # Ties at same loss: prefer sparser model, then smaller λ.
        if (
            loss < best_loss
            or (loss == best_loss and nnz < np.count_nonzero(best_w))
            or (loss == best_loss and nnz == np.count_nonzero(best_w) and lam < best_lambda)
        ):
            best_loss = loss
            best_w = w.copy()
            best_lambda = lam

    return best_w, best_lambda, loss_vals
