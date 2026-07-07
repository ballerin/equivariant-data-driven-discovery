import numpy as np

# ==================================================================================================
# STRidge — Sequential Thresholding Ridge Regression
#
# Brunton, Proctor, Kutz (2016), "Discovering governing equations from data by sparse identification
# of nonlinear dynamical systems", doi: 10.1073/pnas.1517384113.
#
# Rudy, Brunton, Proctor, Kutz (2017), "Data-driven discovery of partial differential equations",
# doi: 10.1126/sciadv.1602614.
#
# Regression-based sparse identification of PDEs (PDE-FIND) uses a sequential thresholding
# ridge regression (STRidge) inner loop to discover sparse coefficients.  The outer loop
# performs cross-validation over the sparsity tolerance.
#
# Code inspired by https://github.com/snagcliffs/PDE-FIND
# ==================================================================================================

def STRidge(X0, y, lam, maxit, tol, normalize=2, early_stop=True, unbias=True):
    """Sequential Threshold Ridge regression.

    Performs ridge regression with sequential hard-thresholding of small coefficients to a set of
    datapoints (X0, y) with ridge regularisation parameter lam.

    Parameters
    ----------
    X0         : ndarray, shape (n, d) - library matrix
    y          : ndarray, shape (n,)   - target
    lam        : float - ridge regularisation parameter
    maxit      : int   - maximum thresholding iterations
    tol        : float - hard-threshold level (entries s.t. |w_i| < tol are set to zero)
    normalize  : int   - norm order for column normalisation (0 = none)
    early_stop : bool  - if True (default), stop when the active set stabilises
                         (PDE-FIND original); if False, run until *maxit*
    unbias     : bool  - if True, final pure least-squares re-fit on the support (PDE-FIND);
                         if False, return the ridge-regularised coefficients

    Returns
    -------
    w : ndarray, shape (d,) - sparse coefficient vector
    """
    n, d = X0.shape
    y = y.ravel()

    # Normalize columns to unit norm so thresholding is scale-invariant
    if normalize != 0:
        Mreg = np.array([1.0 / (np.linalg.norm(X0[:, i], normalize) + 1e-14) for i in range(d)])
        X = X0 * Mreg  # broadcast (n, d) * (d,)
    else:
        X = X0.astype(float, copy=True)
        Mreg = np.ones(d)

    # Precompute normal equations once (avoids repetitions)
    XtX = X.T @ X   # (d, d)
    Xty = X.T @ y   # (d,)

    # Initial ridge solve
    w = np.linalg.lstsq(XtX + lam * np.eye(d), Xty, rcond=None)[0]

    # Sequential thresholding loop to compute the sparsity pattern
    big_prev = None # Previous active set used to evaluate early stopping
    for _ in range(maxit):
        # Threshold small coefficients to zero, and stop if none survives
        small = np.abs(w) < tol
        w[small] = 0.0
        big = ~small
        if big.sum() == 0:
            break

        # Early stopping condition
        if early_stop and big_prev is not None and np.array_equal(big, big_prev):
            break
        big_prev = big.copy()

        # Solve on the active submatrix
        idx = np.where(big)[0]
        G = XtX[np.ix_(idx, idx)] + lam * np.eye(len(idx))
        h = Xty[idx]
        try:
            w[big] = np.linalg.solve(G, h)
        except np.linalg.LinAlgError:
            # Fallback for rank-deficient active submatrix
            w[big] = np.linalg.lstsq(X[:, big], y, rcond=None)[0]

    # Pure least-squares on discovered support (to obtain unbiased coefficients)
    if unbias and big.sum() > 0:
        w[big] = np.linalg.lstsq(X[:, big], y, rcond=None)[0]

    return Mreg * w


def TrainSTRidge(
    Theta, ut, lam, d_tol, maxit=25, STR_iters=10, l0_penalty=None,
    normalize=2, split=0.8, seed=0, early_stop=True, unbias=True,
):
    """Train STRidge with cross-validation over the sparsity tolerance.

    Parameters
    ----------
    Theta      : ndarray, shape (M, K) - library matrix evaluated on data
    ut         : ndarray, shape (M,) - time-derivative vector
    lam        : float - ridge regularisation parameter
    d_tol      : float - initial tolerance step for cross-validation
    maxit      : int   - number of tolerance search steps
    STR_iters  : int   - inner STRidge iterations per step
    l0_penalty : float or None - L0 penalty weight; defaults to 0.001*cond(Theta_tr)
    normalize  : int   - column normalisation norm order (passed to STRidge)
    split      : float - fraction of data used for training
    seed       : int   - RNG seed for reproducible train/test splits
    early_stop : bool  — passed through to STRidge (default True)
    unbias     : bool  - passed through to STRidge (default True)

    Returns
    -------
    w : ndarray, shape (K,) - best sparse coefficient vector
    """
    rng = np.random.default_rng(seed)
    n = Theta.shape[0]
    ut = ut.ravel()
    perm = rng.permutation(n)
    k = int(n * split)
    train, test = perm[:k], perm[k:]
    Theta_tr, Theta_te = Theta[train], Theta[test]
    ytr, yte = ut[train], ut[test]

    if l0_penalty is None:
        l0_penalty = 0.001 * np.linalg.cond(Theta_tr)

    # Initial least-squares fit on the training set to get a baseline error
    w_best = np.linalg.lstsq(Theta_tr, ytr, rcond=None)[0]
    err_best = np.linalg.norm(yte - Theta_te @ w_best) + l0_penalty * np.count_nonzero(w_best)
    tol = float(d_tol)

    # Cross-validation loop: increase tolerance until test error increases
    for _ in range(maxit):
        w = STRidge(Theta_tr, ytr, lam, STR_iters, tol, normalize,
                     early_stop=early_stop, unbias=unbias)
        err = np.linalg.norm(yte - Theta_te @ w) + l0_penalty * np.count_nonzero(w)

        if err <= err_best:
            # Update best solution and increase tolerance for next iteration
            err_best, w_best = err, w
            tol += d_tol
        else:
            # Decrease tolerance for next iteration (but not below zero)
            d_tol = 2 * d_tol / (maxit + 1) # Different from PDE-FIND, which I think has a bug
            tol = max(0, tol - 2 * d_tol) + d_tol
    return w_best


# ==================================================================================================
# SR3 — Sparse Relaxed Regularized Regression
#
# Zheng, Peng, et al. (2018), "A unified framework for sparse relaxed regularized regression:
# SR3", IEEE Access 7 (2018): 1404-1423, doi: 10.1109/ACCESS.2018.2886528.
#
# Champion, K., Zheng, P., Aravkin, A. Y., Brunton, S. L., & Kutz, J. N. (2020), "A unified
# sparse optimization framework to learn parsimonious physics-informed models from data",
# IEEE Access, 8, 169259-169271, doi: 10.1109/ACCESS.2020.3022825.
#
# SR3 solves the sparse regression problem by introducing making w into an auxiliary variable
# relaxing the hard equality constraint ξ = w into a quadratic penalty (κ/2)‖ξ − w‖².
# The result is a block-coordinate scheme:
#
#   1. ξ-step:  (XᵀX + κI)⁻¹ (Xᵀy + κ w)   - a κ-regularised LS solve, always well-posed
#   2. w-step:  soft_threshold(ξ, λ/κ)     - element-wise L1 shrinkage
#
# Code inspired by https://pysindy.readthedocs.io/en/stable/_modules/pysindy/optimizers/sr3.html
# ==================================================================================================

def SR3(X0, y, lam_l1, kappa=1.0, maxit=100, normalize=2, tol=1e-8):
    """Sparse Relaxed Regularized Regression (SR3).

    Solves the relaxed problem::

        min_{xi, w}  ½‖X0 @ xi - y‖² + (κ/2)‖xi - w‖² + λ‖w‖₁

    via block-coordinate descent:
        xi ← (X₀ᵀX₀ + κI)⁻¹ (X₀ᵀy + κ w)
        w  ← soft_threshold(xi, λ/κ)

    The normal-equation matrix is factorised once (Cholesky) and reused
    across all iterations.  SR3 is significantly more robust than STRidge
    for ill-conditioned design matrices because κI regularises the normal
    equations independently of the data.

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

    # Factorise once (Cholesky for SPD, fallback to generic solve)
    try:
        L = np.linalg.cholesky(A)
        solve_A = lambda rhs: np.linalg.solve(
            L.T, np.linalg.solve(L, rhs)
        )
    except np.linalg.LinAlgError:
        solve_A = lambda rhs: np.linalg.solve(A, rhs)

    threshold = lam_l1 / (kappa + 1e-30)

    # Coordinate descent: alternate between LS solve and soft-threshold
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
    Theta, ut, lam, d_tol, maxit=25, SR3_iters=100, normalize=2, split=0.8, seed=0,
):
    """Train SR3 with cross-validation over the L1 penalty.

    The API mirrors TrainSTRidge.

    Parameters
    ----------
    Theta : ndarray, shape (M, K) — library (design) matrix
    ut    : ndarray, shape (M,)   — time-derivative vector
    lam   : float — κ, the SR3 relaxation parameter
    d_tol : float — initial L1 penalty λ_L1
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
    n = Theta.shape[0]
    ut = ut.ravel()
    perm = rng.permutation(n)
    k = int(n * split)
    train, test = perm[:k], perm[k:]
    Theta_tr, Theta_te = Theta[train], Theta[test]
    ytr, yte = ut[train], ut[test]

    l0_penalty = 0.001 * np.linalg.cond(Theta_tr)

    # Initial least-squares fit on the training set to get a baseline error
    w_best = np.linalg.lstsq(Theta_tr, ytr, rcond=None)[0]
    err_best = np.linalg.norm(yte - Theta_te @ w_best) + l0_penalty * np.count_nonzero(w_best)
    lam_l1 = float(d_tol)

    # Cross-validation loop: increase sparsity penalty until test error rises
    for _ in range(maxit):
        w = SR3(Theta_tr, ytr, lam_l1, kappa=lam, maxit=SR3_iters)
        err = np.linalg.norm(yte - Theta_te @ w) + l0_penalty * np.count_nonzero(w)
        if err <= err_best:
            err_best, w_best = err, w
            lam_l1 += d_tol
        else:
            d_tol = 2 * d_tol / (maxit + 1)
            lam_l1 = max(0, lam_l1 - 2 * d_tol) + d_tol
    return w_best


# ==================================================================================================
# MSTLS — Modified Sequential Thresholding Least Squares
#
# Messenger & Bortz (2021), "Weak SINDy for Partial Differential Equations",
# J. Comput. Phys., 443, 110525, doi: 10.1016/j.jcp.2021.110525.
#
# MSTLS is a pure sequential thresholding least-squares algorithm designed for the weak-form
# SINDy (WSINDy) framework.  Instead of ridge regularisation (adding λI to the normal equations) it 
# relies on the weak-form convolution to suppress noise (built-in in the definition).
#
# Code inspired by https://github.com/MathBioCU/WSINDy_PDE
# ==================================================================================================


def mstls(X, y, lam, M=None, maxit=50, normalize=2):
    """Modified Sequential Thresholding Least-Squares (single λ).

    Implements the sequential thresholding from the WSINDy paper
    (Messenger & Bortz, 2021).

    When *M* is None (default), columns are normalised and simple
    STRidge-style thresholding ``|w_j| < λ`` is used (matching the
    MATLAB ``sparsifyDynamics`` fallback when ``M_scale`` is not
    provided).

    When *M* is provided (a vector of length *d*), the dominant-balance
    bounds from the WSINDy paper are used::

        bnds_j = ‖y‖ / ‖X_j‖ * M_j
        LBs_j = λ * max(1, bnds_j)
        UBs_j = (1/λ) * min(1, bnds_j)

    and a coefficient is zeroed if ``|w_j| < LBs_j`` or ``|w_j| > UBs_j``.
    Column normalisation is skipped in this mode since *M* carries the
    scale information.  The coefficients are also scaled by *M* after
    each solve, matching the MATLAB ``sparsifyDynamics`` path A.

    Parameters
    ----------
    X : ndarray, shape (n, d) — design matrix
    y : ndarray, shape (n,)   — target vector
    lam : float — threshold λ > 0
    M : ndarray, shape (d,) or None — per-column scale factors for
        dominant-balance bounds; if None, use simple thresholding
    maxit : int — maximum thresholding iterations
    normalize : int — column normalisation norm order (0 = none);
        ignored when *M* is not None

    Returns
    -------
    w : ndarray, shape (d,) — sparse coefficient vector
    """
    n, d = X.shape
    y = y.ravel()

    if M is not None:
        # Dominant-balance path (sparsifyDynamics path A): per-column
        # bounds based on column norms relative to the target.
        M = np.asarray(M, dtype=float).ravel()
        col_norms = np.linalg.norm(X, axis=0)
        y_norm = np.linalg.norm(y)
        bnds = y_norm / np.maximum(col_norms, 1e-30) * M
        LBs = lam * np.maximum(1.0, bnds)
        UBs = (1.0 / np.maximum(lam, 1e-30)) * np.minimum(1.0, bnds)
        X_work = X
        Mreg = np.ones(d)
    else:
        # Simple thresholding path (sparsifyDynamics fallback): normalise
        # columns and use a single threshold λ for all.
        if normalize != 0:
            Mreg = np.array(
                [1.0 / (np.linalg.norm(X[:, i], normalize) + 1e-14) for i in range(d)]
            )
            X_work = X * Mreg[np.newaxis, :]
        else:
            X_work = X.copy()
            Mreg = np.ones(d)
        LBs = None
        UBs = None

    w = np.linalg.lstsq(X_work, y, rcond=None)[0]

    # Apply initial scaling when using dominant-balance bounds
    if M is not None:
        w = M * w

    active_prev = None
    for _ in range(maxit):
        if M is not None:
            active = (np.abs(w) >= LBs) & (np.abs(w) <= UBs)
        else:
            active = np.abs(w) >= lam

        if active_prev is not None and np.array_equal(active, active_prev):
            break
        active_prev = active.copy()
        w[~active] = 0.0
        if np.sum(active) == 0:
            break
        w[active] = np.linalg.lstsq(X_work[:, active], y, rcond=None)[0]

        # Re-apply scaling after solve (matching MATLAB's M.*(...))
        if M is not None:
            w[active] = M[active] * w[active]

    return Mreg * w


def select_lambda_grid(X, y, n_lambdas=50):
    """Auto-generate a log-spaced grid of candidate λ values for MSTLS.

    .. math::

        \\lambda_{\\max} = \\max_j |X_j^T y| / \\|X_j\\|^2
        \\lambda_{\\min} = \\|X w_{LS}\\| / (d \\cdot \\max_j \\|X_j\\|)

    Parameters
    ----------
    X : ndarray, shape (n, d)
    y : ndarray, shape (n,)
    n_lambdas : int — number of λ candidates

    Returns
    -------
    lambdas : ndarray, shape (n_lambdas,)
        Log-spaced from λ_min to λ_max.
    """
    _, d = X.shape
    y = y.ravel()
    col_norms = np.sqrt(np.sum(X ** 2, axis=0))
    # λ_max = largest regression coefficient magnitude
    Xty = np.abs(X.T @ y)
    lam_max = np.max(Xty / np.maximum(col_norms ** 2, 1e-30))

    # λ_min from projection magnitude
    w_ls = np.linalg.lstsq(X, y, rcond=None)[0]
    proj_norm = np.linalg.norm(X @ w_ls)
    lam_min = proj_norm / (d * max(np.max(col_norms), 1e-30))

    # Clamp to sensible range
    lam_min = max(lam_min, 1e-10)
    lam_max = min(lam_max, 1e4)
    if lam_min >= lam_max:
        lam_min = lam_max * 1e-4

    return np.logspace(np.log10(lam_min), np.log10(lam_max), n_lambdas)


def TrainMSTLS(Theta, ut, lambdas=None, M=None, n_lambdas=50, maxit=50, alpha=0.5, normalize=2):
    """MSTLS with automatic λ selection via the loss function.

    For each λ in the candidate grid, runs ``mstls`` and evaluates

    .. math::

        \\mathcal{L}(\\lambda) =
        \\alpha \\frac{\\|\\Theta(w^\\lambda - w^{LS})\\|}{\\|\\Theta w^{LS}\\|}
        + (1-\\alpha) \\frac{\\#\\{j : w^\\lambda_j \\ne 0\\}}{d}.

    Returns the coefficient vector with the *smallest* λ that minimises
    :math:`\\mathcal{L}`.  Ties are broken by choosing the sparser model.

    Parameters
    ----------
    Theta : ndarray, shape (n, d) — library (design) matrix
    ut    : ndarray, shape (n,)   — time-derivative vector
    lambdas : ndarray or None — candidate λ values; auto-generated if None
    M : ndarray, shape (d,) or None — per-column scale factors for
        dominant-balance bounds in ``mstls``; if None, simple thresholding
    n_lambdas : int — used only when ``lambdas is None``
    maxit : int — inner MSTLS iterations per λ
    alpha : float — projection-cost weight in [0, 1]
    normalize : int — column normalisation (passed to ``mstls``)

    Returns
    -------
    w_best : ndarray, shape (d,)
    lambda_best : float
    loss_vals : ndarray, shape (n_lambdas,)
    """
    _, d = Theta.shape
    y = ut.ravel()

    # Normalize columns (matching what mstls does internally) so that
    # the lambda grid and loss computations use the same scale.
    # When M is provided, skip normalisation here too — M carries the scale.
    if M is not None:
        Theta_norm = Theta
        Mreg = np.ones(d)
    elif normalize != 0:
        Mreg = np.array(
            [1.0 / (np.linalg.norm(Theta[:, i], normalize) + 1e-14) for i in range(d)]
        )
        Theta_norm = Theta * Mreg[np.newaxis, :]
    else:
        Theta_norm = Theta.copy()
        Mreg = np.ones(d)

    if lambdas is None:
        lambdas = select_lambda_grid(Theta_norm, y, n_lambdas)

    w_ls = np.linalg.lstsq(Theta_norm, y, rcond=None)[0]
    Theta_w_ls_norm = np.linalg.norm(Theta_norm @ w_ls)
    if Theta_w_ls_norm < 1e-30:
        Theta_w_ls_norm = 1.0

    best_loss = np.inf
    best_w = np.zeros(d)  # unscaled
    best_lambda = lambdas[0]
    loss_vals = np.empty(len(lambdas))

    for i, lam in enumerate(lambdas):
        w = mstls(Theta, ut, lam, M=M, maxit=maxit, normalize=normalize)
        # Loss computed on NORMALIZED system for consistency with lambda grid
        w_norm = w / np.maximum(Mreg, 1e-30)
        proj_cost = np.linalg.norm(Theta_norm @ (w_norm - w_ls)) / Theta_w_ls_norm
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
