import numpy as np

# ==================================================================================================
# STRidge - Sequential Thresholding Ridge Regression
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
    early_stop : bool  - passed through to STRidge (default True)
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
# SR3 - Sparse Relaxed Regularized Regression
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
    X0     : ndarray, shape (n, d) - design matrix
    y      : ndarray, shape (n,)   - target
    lam_l1 : float - L₁ penalty λ (sparsity)
    kappa  : float - relaxation strength κ (larger → tighter ξ≈w coupling)
    maxit  : int   - coordinate-descent iterations
    normalize : int - column normalisation norm order (0 = none)
    tol    : float - convergence tolerance on coefficient change

    Returns
    -------
    w : ndarray, shape (d,) - sparse coefficient vector
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
    Theta : ndarray, shape (M, K) - library (design) matrix
    ut    : ndarray, shape (M,)   - time-derivative vector
    lam   : float - κ, the SR3 relaxation parameter
    d_tol : float - initial L1 penalty λ_L1
    maxit     : int - number of cross-validation steps
    SR3_iters : int - inner SR3 coordinate-descent iterations per step
    normalize : int - column normalisation (passed to SR3)
    split     : float - train/test split fraction
    seed      : int - RNG seed for reproducible train/test splits

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
# MSTLS - Modified Sequential Thresholding Least Squares
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


def mstls(X, y, lam, M=None, maxit=None):
    """Single-λ MSTLS iteration - dominant-balance thresholding.

    Implements the sequential thresholding from Messenger & Bortz (2021),
    matching the reference ``_mstls_iterate`` from the WSINDy codebase
    (MathBioCU/WSciML-SIAMWorkshop2026,
    ``libs/pywsindy_base/wsindy_obj/_weakform.py``).

    Always uses per-column dominant-balance bounds.  When *M* is None
    it defaults to ``M = ones(d)``, matching the reference behaviour
    (the MATLAB ``sparsifyDynamics`` only uses simple thresholding when
    ``M_scale`` is empty, but the higher-level ``wsindy_pde_RGLS_seq``
    always provides *M* - the simple-thresholding fallback is a
    STRidge-ism, not the MSTLS algorithm).

    The least-squares solves are done in the original ("scaled") space
    spanned by the columns of *X*.  The per-column scale vector *M*
    converts from scaled to physical coefficients: ``w_phys = M ⊙ w_sc``.
    Thresholding is applied to the physical coefficients against
    physical bounds::

        bnds_j = ‖y‖ / ‖X_j‖ · M_j
        LBs_j = λ · max(1, bnds_j)
        UBs_j = (1/λ) · min(1, bnds_j)

    and a coefficient is zeroed if ``|w_j| < LBs_j`` or ``|w_j| > UBs_j``.

    Parameters
    ----------
    X : ndarray, shape (n, d) - design matrix
    y : ndarray, shape (n,)   - target vector
    lam : float - threshold λ > 0
    M : ndarray, shape (d,) or None - per-column scale factors for
        dominant-balance bounds; if None, defaults to ``ones(d)``
    maxit : int or None - maximum thresholding iterations;
        if None, defaults to ``d + 1`` (matching the reference)

    Returns
    -------
    w : ndarray, shape (d,) - sparse coefficient vector (physical)
    """
    n, d = X.shape
    y = y.ravel()

    if M is None:
        M = np.ones(d)
    else:
        M = np.asarray(M, dtype=float).ravel()

    if maxit is None:
        maxit = d + 1

    col_norms = np.linalg.norm(X, axis=0)
    col_norms[col_norms == 0] = np.inf
    bnds = (np.linalg.norm(y) / col_norms) * M
    lower = lam * np.maximum(bnds, 1.0)
    upper = (1.0 / lam) * np.minimum(bnds, 1.0)

    # Initial solve: scaled space → physical coefficients
    w = M * np.linalg.lstsq(X, y, rcond=None)[0]

    prev = None
    for _ in range(maxit):
        small = (np.abs(w) < lower) | (np.abs(w) > upper)
        if prev is not None and np.array_equal(small, prev):
            break
        big = ~small
        w = np.zeros(d)
        if big.any():
            w[big] = M[big] * np.linalg.lstsq(X[:, big], y, rcond=None)[0]
        prev = small

    return w  # physical coefficients


def select_lambda_grid(X, y, n_lambdas=50):
    """Auto-generate a log-spaced grid of candidate λ values for MSTLS.

    .. math::

        \\lambda_{\\max} = \\max_j |X_j^T y| / \\|X_j\\|^2
        \\lambda_{\\min} = \\|X w_{LS}\\| / (d \\cdot \\max_j \\|X_j\\|)

    Parameters
    ----------
    X : ndarray, shape (n, d)
    y : ndarray, shape (n,)
    n_lambdas : int - number of λ candidates

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

    # Safety floor / ceiling to avoid log of zero / overflow;
    # the original MATLAB code has no clamps, so these are deliberately
    # wide (1e-30, 1e30) and should never be hit in practice.
    lam_min = max(lam_min, 1e-30)
    lam_max = min(lam_max, 1e30)
    if lam_min >= lam_max:
        lam_min = lam_max * 1e-4

    return np.logspace(np.log10(lam_min), np.log10(lam_max), n_lambdas)


def TrainMSTLS(Theta, ut, lambdas=None, M=None, n_lambdas=50, maxit=None, alpha=0.5):
    """MSTLS with automatic λ selection via the loss function.

    Matches the reference ``mstls`` + ``_mstls_loss`` from the WSINDy
    codebase (MathBioCU/WSciML-SIAMWorkshop2026).  For each λ in the
    candidate grid, runs :func:`mstls` and evaluates

    .. math::

        \\mathcal{L}(\\lambda) =
        \\alpha \\frac{\\|\\Theta(w_{\\text{sc}}^\\lambda - w^{LS})\\|}{\\|\\Theta w^{LS}\\|}
        + (1-\\alpha) \\frac{\\#\\{j : w^\\lambda_j \\ne 0\\}}{d},

    where :math:`w_{\\text{sc}} = w / M` converts physical coefficients
    back to the scaled space in which :math:`\\Theta` and :math:`w^{LS}`
    live.  The reference uses equal weights (α = 0.5).

    Returns the coefficient vector minimising :math:`\\mathcal{L}`.
    Ties are broken by preferring sparser models, then smaller λ.

    Parameters
    ----------
    Theta : ndarray, shape (n, d) - library (design) matrix
    ut    : ndarray, shape (n,)   - time-derivative vector
    lambdas : ndarray or None - candidate λ values; auto-generated if None
    M : ndarray, shape (d,) or None - per-column scale factors for
        dominant-balance bounds in :func:`mstls`; if None, defaults to
        ``ones(d)`` (matching the reference fallback)
    n_lambdas : int - used only when ``lambdas is None``
    maxit : int or None - inner MSTLS iterations per λ; if None,
        defaults to ``d + 1``
    alpha : float - projection-cost weight in [0, 1]; the reference
        uses α = 0.5 (equal weight)

    Returns
    -------
    w_best : ndarray, shape (d,)
    lambda_best : float
    loss_vals : ndarray, shape (n_lambdas,)
    """
    _, d = Theta.shape
    y = ut.ravel()

    if M is None:
        M = np.ones(d)
    else:
        M = np.asarray(M, dtype=float).ravel()

    if lambdas is None:
        lambdas = select_lambda_grid(Theta, y, n_lambdas)

    # LS solution in original (scaled) space
    w_ls = np.linalg.lstsq(Theta, y, rcond=None)[0]
    denom = np.linalg.norm(Theta @ w_ls)
    if denom < 1e-30:
        denom = 1.0

    best_loss = np.inf
    best_w = np.zeros(d)
    best_lambda = lambdas[0]
    loss_vals = np.empty(len(lambdas))

    for i, lam in enumerate(lambdas):
        w = mstls(Theta, ut, lam, M=M, maxit=maxit)
        # Loss - matches _mstls_loss from the reference:
        #   ls_term   = ‖G @ (w_sc − w_ls)‖ / ‖G @ w_ls‖
        #   sparsity  = nnz(w_phys) / d
        w_sc = w / np.maximum(M, 1e-30)          # physical → scaled
        ls_term = np.linalg.norm(Theta @ (w_sc - w_ls)) / denom
        nnz = np.count_nonzero(w)
        sparsity_term = nnz / max(d, 1)
        loss = alpha * ls_term + (1.0 - alpha) * sparsity_term
        loss_vals[i] = loss

        # Tie-breaking: smaller loss → sparser model → smaller λ
        if (
            loss < best_loss
            or (loss == best_loss and nnz < np.count_nonzero(best_w))
            or (loss == best_loss and nnz == np.count_nonzero(best_w) and lam < best_lambda)
        ):
            best_loss = loss
            best_w = w.copy()
            best_lambda = lam

    return best_w, best_lambda, loss_vals


# ==================================================================================================
# RFE - Recursive Feature Elimination
#
# Stephany & Earls (2022), "PDE-READ: Human-readable Partial Differential Equation
# Discovery using Deep Learning", 10.1016/j.neunet.2022.07.008.
#
# RFE is a parameter-free sparse regression algorithm.  It iteratively eliminates the smallest 
# OLS (Ordinary Least Squares) coefficient, tracks the residual at each step, and selects the model 
# at the step that maximises the relative residual increase:
#
#     j* = argmax_j  (RLS[j+1] - RLS[j]) / RLS[j]
#
# This is equivalent to selecting the step where removing a feature causes the largest
# relative jump in residual, i.e. the "elbow" in the residual curve.
#
# Code inspired by github.com/punkduckable/PDE-READ
# ==================================================================================================


def RFE(X0, y, normalize=2, lam=0.0):
    """Recursive Feature Elimination.

    Sparse regression for PDE discovery.  Iteratively eliminates the smallest
    (optionally ridge-regularised) coefficient, tracks the residual at each
    step, and selects the model at the step that maximises the relative
    residual increase.

    The algorithm:
        1. Normalises columns to unit norm (default: L2).
        2. Fits (ridge-regularised) LS on all features; records residual.
        3. Eliminates the feature with smallest |coefficient|.
        4. Re-fits on remaining features; records new residual.
        5. Repeats until zero features remain (residual = ‖y‖²).
        6. Selects step j* = argmax_j (RLS[j+1] - RLS[j]) / RLS[j].
        7. Returns coefficients from step j*, rescaled to original units.

    Parameters
    ----------
    X0        : ndarray, shape (N, d) - design matrix
    y         : ndarray, shape (N,)   - target vector
    normalize : int - norm order for column normalisation (0 = none, 2 = L2)
    lam       : float - ridge regularisation (0 = pure OLS, >0 = ridge)

    Returns
    -------
    w : ndarray, shape (d,) - sparse coefficient vector (original scale)
    """
    N, d = X0.shape
    y = y.ravel()

    # Column normalisation
    if normalize != 0:
        scale = np.array(
            [np.linalg.norm(X0[:, i], normalize) for i in range(d)]
        )
        # Guard against zero columns
        scale = np.maximum(scale, 1e-30)
        X = X0 / scale[np.newaxis, :]
    else:
        X = X0.copy()
        scale = np.ones(d)

    # Precompute normal equations
    XtX = X.T @ X   # (d, d)
    Xty = X.T @ y   # (d,)

    # Initial (ridge-regularised) LS on all features
    try:
        w = np.linalg.solve(XtX + lam * np.eye(d), Xty)
    except np.linalg.LinAlgError:
        w = np.linalg.lstsq(X, y, rcond=None)[0]

    Xi = np.zeros((d, d))
    Xi[:, 0] = w
    RLS = [np.sum((y - X @ w) ** 2)]

    # Active-set elimination loop
    S = np.ones(d, dtype=bool)      # True = active
    for j in range(d - 1):          # d-1 eliminations. Last one is trivial residual
        # Find smallest |coefficient| among active features
        active_coeffs = np.where(S, np.abs(w), np.inf)
        k_min = int(np.argmin(active_coeffs))
        S[k_min] = False

        # Re-fit (ridge-regularised) LS on remaining active features
        idx = np.where(S)[0]
        w_new = np.zeros(d)
        n_active = len(idx)
        try:
            w_new[idx] = np.linalg.solve(
                XtX[np.ix_(idx, idx)] + lam * np.eye(n_active), Xty[idx],
            )
        except np.linalg.LinAlgError:
            w_new[idx] = np.linalg.lstsq(X[:, idx], y, rcond=None)[0]

        Xi[:, j + 1] = w_new
        RLS.append(np.sum((y - X @ w_new) ** 2))
        w = w_new

    # Zero-feature residual is just the squared norm of y (zero coefficients)
    RLS.append(np.sum(y ** 2))

    # Select best step via relative residual jump
    RLS_arr = np.array(RLS, dtype=float)
    ratios = np.zeros(d)
    for j in range(d):
        denom = max(RLS_arr[j], 1e-30)
        ratios[j] = (RLS_arr[j + 1] - RLS_arr[j]) / denom
    j_max = int(np.argmax(ratios))

    # Rescale to original units
    w_out = Xi[:, j_max] / scale
    return w_out


def TrainRFE(Theta, ut, lam=0.0, **kwargs):
    """Thin wrapper around ``RFE`` for API consistency with other regression
    methods (STRidge, SR3, MSTLS).

    Parameters
    ----------
    Theta : ndarray - design matrix
    ut    : ndarray - target vector
    lam   : float - ridge regularisation (0 = pure OLS, >0 = ridge-RFE)
    """
    return RFE(Theta, ut, lam=lam)


# =====================================================================================
# Exhaustive best-subset selection  (L0-penalised least squares)
#
#   min_{β}   ‖Xβ − y‖² / ‖y‖²   +   λ · ‖β‖₀
#
# Enumerates every subset of columns and returns the Pareto frontier of
# (model size, relative RSS) plus the single best model for a given λ.

# The key precomputation trick is to work entirely with the d×d Gram matrix
# G = ΘᵀΘ and h = ΘᵀFₜ, so each subset costs O(k³) independent of n.
# =====================================================================================

import itertools
import math
import time


def exhaustive_best_subset(X0, y, max_k=None, normalize=2, l0_penalty=0.0,
                           verbose=False):
    """Exhaustive best-subset selection with L₀ sparsity penalty.

    Enumerates every subset of columns of size 1 up to *max_k* (default: all
    d columns), fits OLS on each, and returns the Pareto frontier.  The best
    model minimises::

        relative_RSS  +  λ · k,      relative_RSS = ‖Xβ - y‖² / ‖y‖².

    Parameters
    ----------
    X0 : ndarray, shape (n, d) - design matrix.
    y  : ndarray, shape (n,)  - target vector.
    max_k : int or None - maximum subset size (None = all d columns).
    normalize : int - column normalisation (0 = none, 2 = L2).
    l0_penalty : float - λ ≥ 0, the sparsity penalty per active term.
        λ = 0 picks the full model; larger λ enforces fewer terms.
    verbose : bool - print progress.

    Returns
    -------
    best_by_k : dict
        k → (rss, indices, xi, score) for the best k-term model.
    best_overall : tuple
        (k, indices, xi, rss, score) for the L0-optimal model.
    """
    _, d = X0.shape
    y = y.ravel()

    # column normalisation (same convention as STRidge / SR3)
    if normalize != 0:
        Mreg = np.array([
            1.0 / (np.linalg.norm(X0[:, i], normalize) + 1e-14)
            for i in range(d)
        ])
        X = X0 * Mreg[np.newaxis, :]
    else:
        X = X0.astype(float, copy=True)
        Mreg = np.ones(d)

    # precompute normal equations
    XtX = X.T @ X       # (d, d)
    Xty = X.T @ y       # (d,)
    yty = float(y @ y)  # normalize RSS to relative scale

    max_k = min(max_k if max_k is not None else d, d)
    if max_k == d and d > 25:
        print(
            f"Warning: full enumeration (2^{d} subsets) may be slow. "
            f"Consider setting max_k ≤ 8."
        )

    # enumerate
    best_by_k = {}          # k → (rss, indices, xi, score)
    total_subsets = sum(math.comb(d, k) for k in range(1, max_k + 1))
    count = 0
    t0 = time.perf_counter()

    for k in range(1, max_k + 1):
        for combo in itertools.combinations(range(d), k):
            idx = list(combo)
            G = XtX[np.ix_(idx, idx)]   # (k, k)
            h = Xty[idx]                 # (k,)

            try:
                xi_sub = np.linalg.solve(G, h)
            except np.linalg.LinAlgError:
                xi_sub = np.linalg.lstsq(X[:, idx], y, rcond=None)[0]

            rss = float(yty - xi_sub @ h)
            rss = max(rss, 0.0)          # guard against tiny negatives

            score = rss / yty + l0_penalty * float(k)

            if k not in best_by_k or score < best_by_k[k][3]:
                xi_full = np.zeros(d)
                xi_full[idx] = xi_sub
                best_by_k[k] = (rss, idx, xi_full, score)

            count += 1

        if verbose:
            elapsed = time.perf_counter() - t0
            frac = count / max(total_subsets, 1)
            print(
                f"  k={k:2d} done  |  {count}/{total_subsets} subsets  "
                f"({100*frac:.1f}%)  |  {elapsed:.1f}s elapsed"
            )

    # rescale to original units
    best_by_k_rescaled = {}
    for k, (rss, idx, xi, score) in best_by_k.items():
        best_by_k_rescaled[k] = (rss, idx, Mreg * xi, score)

    # pick best overall (lowest L₀-penalised score)
    best_k = min(best_by_k_rescaled, key=lambda k_: best_by_k_rescaled[k_][3])
    best_rss, best_idx, best_xi, best_score = best_by_k_rescaled[best_k]
    best_overall = (best_k, best_idx, best_xi, best_rss, best_score)

    if verbose:
        elapsed = time.perf_counter() - t0
        score0 = best_by_k_rescaled[best_k][3]
        print(
            f"Done: {count} subsets in {elapsed:.1f}s "
            f"({count / max(elapsed, 1e-6):.0f} subsets/s).  "
            f"Best: k={best_k}, indices={best_idx}, "
            f"relRSS={best_rss / yty:.6f}, score={score0:.6f}"
        )

    return best_by_k_rescaled, best_overall


def TrainExhaustive(Theta, ut, max_k=None, normalize=2, l0_penalty=0.0,
                    verbose=False):
    """Thin wrapper around :func:`exhaustive_best_subset` for API consistency
    with ``TrainSTRidge``, ``TrainSR3``, etc.

    Parameters
    ----------
    Theta : ndarray - design matrix.
    ut    : ndarray - target vector.
    max_k : int or None - maximum subset size.
    normalize : int - column normalisation.
    l0_penalty : float - sparsity penalty per term.
    verbose : bool - print progress.

    Returns
    -------
    w : ndarray, shape (d,) - best sparse coefficient vector.
    """
    _, best = exhaustive_best_subset(
        Theta, ut, max_k=max_k, normalize=normalize,
        l0_penalty=l0_penalty, verbose=verbose,
    )
    return best[2]  # xi vector
