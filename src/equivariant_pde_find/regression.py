"""Sparse regression algorithms for PDE discovery: STRidge, SR3, MSTLS,
RFE, and brute-force exhaustive best-subset search, plus their Train*
wrappers used by the sweep scripts.

Also provides ``centering``, an opt-in wrapper that partials the intercept
out (Frisch-Waugh-Lovell) before delegating to any of the Train* solvers,
so the constant column does not compete in the sparsification."""

import numpy as np


# ==================================================================================================
# SolverCache - reuse data-only work across solver calls on the same matrix
#
# A sweep calls the Train* solvers many times on one (Theta, ut): 20 STRidge (lam, d_tol) settings,
# several RFE lam values, several MSTLS alphas, each with and without centering.  Everything that
# depends only on the data (train/test split, cond, column normalisation, normal equations,
# least-squares refits on a given support, the MSTLS lambda path, the centred system) is identical
# across those calls.  Passing one SolverCache to all of them computes it once.  The arithmetic is
# unchanged, so results are bit-identical to calls without a cache.
# ==================================================================================================

class SolverCache:
    """Memo for the ``cache=`` argument of ``TrainSTRidge``, ``TrainRFE``,
    ``TrainMSTLS`` and ``centering``.

    Entries are keyed by the *identity* of the input arrays (the cache holds a
    reference, so an id cannot be reused while it is alive) plus the options
    that affect the precomputation.  The arrays must therefore not be modified
    in place while the cache is in use.  Create one cache per design matrix
    and drop it when moving to the next, since entries hold copies of the
    matrix (the training split, the normalised or centred columns).
    """

    def __init__(self):
        self._entries = {}

    def entry(self, arrays, options=()):
        """Return the dict stored for ``arrays`` (a tuple) and ``options``,
        creating an empty one on first use."""
        key = tuple(id(a) for a in arrays) + tuple(options)
        e = self._entries.get(key)
        if e is None or any(r is not a for r, a in zip(e["_refs"], arrays)):
            e = {"_refs": tuple(arrays)}
            self._entries[key] = e
        return e

    def clear(self):
        self._entries.clear()


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
    y = y.ravel()
    X, Mreg, XtX, Xty = _stridge_setup(X0, y, normalize)
    return _stridge_iterate(X, y, Mreg, XtX, Xty, lam, maxit, tol,
                            early_stop=early_stop, unbias=unbias)


def _stridge_setup(X0, y, normalize):
    """Column normalisation and normal equations for ``STRidge``.

    These depend only on the data, not on ``lam`` or ``tol``, so
    ``TrainSTRidge`` computes them once and reuses them for every tolerance
    it tries.
    """
    d = X0.shape[1]

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
    return X, Mreg, XtX, Xty


def _support_lstsq(X, y, big, refits):
    """Least squares of ``y`` on the columns ``X[:, big]``.

    ``refits`` (a dict, or None to disable) memoises the result by support.
    The same support recurs many times within a ``TrainSTRidge`` call and
    across the λ path of ``TrainMSTLS``, and with many rows this solve is the
    dominant cost.
    """
    if refits is None:
        return np.linalg.lstsq(X[:, big], y, rcond=None)[0]
    key = big.tobytes()
    if key not in refits:
        refits[key] = np.linalg.lstsq(X[:, big], y, rcond=None)[0]
    return refits[key]


def _stridge_iterate(X, y, Mreg, XtX, Xty, lam, maxit, tol, early_stop=True,
                     unbias=True, refits=None):
    """The thresholding loop of ``STRidge`` on precomputed ``_stridge_setup``
    output.  ``refits`` is passed to ``_support_lstsq``."""
    d = X.shape[1]

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
            w[big] = _support_lstsq(X, y, big, refits)

    # Pure least-squares on discovered support (to obtain unbiased coefficients)
    if unbias and big.sum() > 0:
        w[big] = _support_lstsq(X, y, big, refits)

    return Mreg * w


def TrainSTRidge(
    Theta, ut, lam, d_tol, maxit=125, STR_iters=10, l0_penalty=None,
    normalize=2, split=0.8, seed=0, early_stop=True, unbias=True, cache=None,
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
    cache      : SolverCache or None - reuse the split, cond, normalisation,
                 normal equations and support refits across calls on the
                 same (Theta, ut); results are unchanged

    Returns
    -------
    w : ndarray, shape (K,) - best sparse coefficient vector
    """
    # Everything below up to the loop depends only on (Theta, ut, split, seed,
    # normalize), not on lam or d_tol.  STRidge's own setup is independent of
    # tol too, so it is computed once here rather than in each of up to maxit
    # STRidge calls, and support refits are memoised.  The arithmetic is
    # unchanged, so w is identical to calling STRidge in the loop.
    pre = (cache.entry((Theta, ut), ("stridge", split, seed, normalize))
           if cache is not None else {})
    if "Theta_tr" not in pre:
        rng = np.random.default_rng(seed)
        n = Theta.shape[0]
        y = ut.ravel()
        perm = rng.permutation(n)
        k = int(n * split)
        train, test = perm[:k], perm[k:]
        pre["Theta_tr"], pre["Theta_te"] = Theta[train], Theta[test]
        pre["ytr"], pre["yte"] = y[train], y[test]
        # Initial least-squares fit on the training set to get a baseline error
        pre["w_ls"] = np.linalg.lstsq(pre["Theta_tr"], pre["ytr"], rcond=None)[0]
        pre["res_ls"] = np.linalg.norm(pre["yte"] - pre["Theta_te"] @ pre["w_ls"])
        pre["setup"] = _stridge_setup(pre["Theta_tr"], pre["ytr"], normalize)
        pre["refits"] = {}
    Theta_te, ytr, yte = pre["Theta_te"], pre["ytr"], pre["yte"]

    if l0_penalty is None:
        if "cond" not in pre:
            pre["cond"] = np.linalg.cond(pre["Theta_tr"])
        l0_penalty = 0.001 * pre["cond"]

    w_best = pre["w_ls"].copy()
    err_best = pre["res_ls"] + l0_penalty * np.count_nonzero(w_best)
    tol = float(d_tol)

    X_tr, Mreg, XtX, Xty = pre["setup"]
    refits = pre["refits"]

    # Cross-validation loop: increase tolerance until test error increases
    for _ in range(maxit):
        w = _stridge_iterate(X_tr, ytr, Mreg, XtX, Xty, lam, STR_iters, tol,
                             early_stop=early_stop, unbias=unbias,
                             refits=refits)
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
    Theta, ut, lam, d_tol, maxit=125, SR3_iters=125, normalize=2, split=0.8, seed=0,
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
    y = y.ravel()
    return _mstls_iterate(X, y, lam, M=M, maxit=maxit)


def _mstls_iterate(X, y, lam, M=None, maxit=None, pre=None):
    """``mstls`` body.  ``pre`` (a dict, or None) memoises the λ-independent
    pieces - column norms, ‖y‖, the full least-squares solve - and the
    least-squares refits by support, for reuse across λ on the same data."""
    n, d = X.shape

    if M is None:
        M = np.ones(d)
    else:
        M = np.asarray(M, dtype=float).ravel()

    if maxit is None:
        maxit = d + 1

    if pre is None:
        pre = {}
    if "col_norms" not in pre:
        col_norms = np.linalg.norm(X, axis=0)
        col_norms[col_norms == 0] = np.inf
        pre["col_norms"] = col_norms
        pre["y_norm"] = np.linalg.norm(y)
    if "w_ls" not in pre:
        pre["w_ls"] = np.linalg.lstsq(X, y, rcond=None)[0]
    refits = pre.setdefault("refits", {})

    bnds = (pre["y_norm"] / pre["col_norms"]) * M
    lower = lam * np.maximum(bnds, 1.0)
    upper = (1.0 / lam) * np.minimum(bnds, 1.0)

    # Initial solve: scaled space → physical coefficients
    w = M * pre["w_ls"]

    prev = None
    for _ in range(maxit):
        small = (np.abs(w) < lower) | (np.abs(w) > upper)
        if prev is not None and np.array_equal(small, prev):
            break
        big = ~small
        w = np.zeros(d)
        if big.any():
            w[big] = M[big] * _support_lstsq(X, y, big, refits)
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
    y = y.ravel()
    return _select_lambda_grid(X, y, n_lambdas)


def _select_lambda_grid(X, y, n_lambdas, w_ls=None):
    """``select_lambda_grid`` body; ``w_ls`` is the least-squares solution of
    (X, y) when the caller already has it."""
    _, d = X.shape
    col_norms = np.sqrt(np.sum(X ** 2, axis=0))
    # λ_max = largest regression coefficient magnitude
    Xty = np.abs(X.T @ y)
    lam_max = np.max(Xty / np.maximum(col_norms ** 2, 1e-30))

    # λ_min from projection magnitude
    if w_ls is None:
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


def TrainMSTLS(Theta, ut, lambdas=None, M=None, n_lambdas=50, maxit=None, alpha=0.5,
               cache=None):
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
    cache : SolverCache or None - reuse the λ grid, the least-squares
        solution, the support refits and the per-λ MSTLS solutions (which do
        not depend on α) across calls on the same (Theta, ut); results are
        unchanged

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

    # Data-only work, shared by every call with the same data (and, for the
    # λ path, the same M, grid and maxit).
    data = cache.entry((Theta, ut), ("mstls",)) if cache is not None else {}
    if "w_ls" not in data:
        # LS solution in original (scaled) space
        data["w_ls"] = np.linalg.lstsq(Theta, y, rcond=None)[0]
    w_ls = data["w_ls"]

    if lambdas is None:
        grids = data.setdefault("grids", {})
        if n_lambdas not in grids:
            grids[n_lambdas] = _select_lambda_grid(Theta, y, n_lambdas, w_ls=w_ls)
        lambdas = grids[n_lambdas]
    lambdas = np.asarray(lambdas)

    if "denom" not in data:
        denom = np.linalg.norm(Theta @ w_ls)
        data["denom"] = 1.0 if denom < 1e-30 else denom
    denom = data["denom"]

    paths = data.setdefault("paths", {})
    path_key = (M.tobytes(), lambdas.tobytes(), maxit)
    path = paths.setdefault(path_key, {})

    best_loss = np.inf
    best_w = np.zeros(d)
    best_lambda = lambdas[0]
    loss_vals = np.empty(len(lambdas))

    for i, lam in enumerate(lambdas):
        if i not in path:
            w = _mstls_iterate(Theta, y, lam, M=M, maxit=maxit, pre=data)
            # Loss - matches _mstls_loss from the reference:
            #   ls_term   = ‖G @ (w_sc − w_ls)‖ / ‖G @ w_ls‖
            #   sparsity  = nnz(w_phys) / d
            w_sc = w / np.maximum(M, 1e-30)          # physical → scaled
            path[i] = (w, np.linalg.norm(Theta @ (w_sc - w_ls)) / denom)
        w, ls_term = path[i]
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
    y = y.ravel()
    return _rfe(X0, y, normalize=normalize, lam=lam)


def _rfe_setup(X0, y, normalize):
    """Column normalisation and normal equations for ``RFE`` (independent of
    ``lam``)."""
    d = X0.shape[1]

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
    return X, scale, XtX, Xty


def _rfe(X0, y, normalize=2, lam=0.0, setup=None):
    """``RFE`` body; ``setup`` is a precomputed ``_rfe_setup`` result."""
    d = X0.shape[1]
    X, scale, XtX, Xty = setup if setup is not None else _rfe_setup(X0, y, normalize)

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


def TrainRFE(Theta, ut, lam=0.0, cache=None, **kwargs):
    """Thin wrapper around ``RFE`` for API consistency with other regression
    methods (STRidge, SR3, MSTLS).

    Parameters
    ----------
    Theta : ndarray - design matrix
    ut    : ndarray - target vector
    lam   : float - ridge regularisation (0 = pure OLS, >0 = ridge-RFE)
    cache : SolverCache or None - reuse the column normalisation and normal
            equations across calls (e.g. several ``lam``) on the same
            (Theta, ut); results are unchanged
    """
    if cache is None:
        return RFE(Theta, ut, lam=lam)
    y = ut.ravel()
    pre = cache.entry((Theta, ut), ("rfe", 2))
    if "setup" not in pre:
        pre["setup"] = _rfe_setup(Theta, y, 2)
    return _rfe(Theta, y, normalize=2, lam=lam, setup=pre["setup"])


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


# ==================================================================================================
# Centering - intercept partialling for sparse PDE regression
#
# Frisch, R. & Waugh, F. V. (1933), "Partial time regressions as compared with individual
# trends", Econometrica 1(4), 387-401, doi: 10.2307/1907330.
#
# Lovell, M. C. (1963), "Seasonal adjustment of economic time series and multiple regression
# analysis", J. Amer. Statist. Assoc. 58(304), 993-1010, doi: 10.1080/01621459.1963.10480682.
#
# The libraries in this package carry an explicit constant column "1", and the solvers above
# fit it as an ordinary feature: it is column-normalised alongside the others and competes in
# the thresholding / elimination loop.  The constant column is however exactly the direction
# that mean-centring removes, so it is collinear with the mean component of every other
# column.  This inflates cond(Theta) - and hence TrainSTRidge's default
# l0_penalty = 0.001 * cond(Theta_tr) - and lets the sparsifier trade a dropped constant
# against a corrupted coefficient on a genuine term.
#
# Centering partials the intercept out instead (Frisch-Waugh-Lovell): drop the constant
# column, centre y and the remaining columns, run the sparse solver on the centred system
# where the constant cannot compete, then recover the intercept in closed form
#
#     w0 = mean(y) - w . mean(X00, axis=0)
#
# and zero it when |w0| < tol * mean(|y|).  For plain OLS steps 3-4 reproduce the full-library
# solution exactly (that is the FWL identity); for a sparsifying solver the selection path
# genuinely differs, which is the point of the method.
#
# This is the convention used by scikit-learn's linear models (``fit_intercept=True`` centres
# X and y in ``_preprocess_data`` and leaves the intercept unpenalised).  Note it is the
# *opposite* of the PySINDy convention, whose ``BaseOptimizer`` passes ``fit_intercept=False``
# and documents that "intercept is the responsibility of feature libraries", and of PDE-FIND,
# which carries Theta(U) = [1, U, U^2, ...] and fits the constant as a feature.
# ==================================================================================================

def _resolve_const_index(Theta, descs=None, const_index=None, const_desc="1"):
    """Locate the constant column of a library matrix.

    Resolution order: explicit ``const_index``, then ``descs.index(const_desc)``,
    then a numeric search for a zero-variance column.

    Parameters
    ----------
    Theta      : ndarray, shape (K, m) - library matrix
    descs      : sequence of str or None - column descriptions
    const_index: int or None - explicit column index (may be negative)
    const_desc : str - description marking the constant column

    Returns
    -------
    int or None - index of the constant column, or None if there is none.
    """
    m = Theta.shape[1]

    if const_index is not None:
        j = int(const_index) % m
        return j

    if descs is not None:
        descs = list(descs)
        if len(descs) != m:
            raise ValueError(
                f"descs has length {len(descs)} but Theta has {m} columns."
            )
        if const_desc in descs:
            return descs.index(const_desc)
        return None

    # Numeric fallback: a column is constant if its spread is negligible relative
    # to its own level.  Same convention as diagnostics.vif_factors, which skips
    # columns with std ~ 0.
    col_mean = Theta.mean(axis=0)
    col_std = Theta.std(axis=0)
    scale = np.maximum(1.0, np.abs(col_mean))
    const_cols = np.where(col_std < 1e-12 * scale)[0]
    if const_cols.size == 0:
        return None
    return int(const_cols[0])


def centering(Theta, ut, algorithm, tol=0.0, descs=None, const_index=None,
              const_desc="1", allow_missing=False, return_intercept=False,
              zero_tol=1e-12, cache=None):
    """Fit a sparse model with the intercept partialled out.

    Drops the constant column, centres the target and the remaining columns,
    delegates to ``algorithm`` on the centred system, then recovers the
    intercept in closed form and thresholds it against ``tol * mean(|ut|)``.

    The returned coefficient vector has the *same length and ordering* as
    ``Theta`` has columns, with the recovered intercept written back into the
    constant column's slot.  It is therefore a drop-in replacement for the
    output of ``TrainSTRidge`` / ``TrainRFE`` and needs no changes to
    ``evaluate_native``, ``evaluate_unified``, or the usual
    ``E_data = ||Theta @ w - ut|| / ||ut||``.

    Centring uses the mean over *all* rows, before ``algorithm`` applies any
    internal train/test split (as ``TrainSTRidge`` does).  One scalar mean per
    column therefore crosses from test to train; this is deliberate, matching
    the reference algorithm, and negligible in practice.

    Columns that become degenerate *because* of centring - any other constant
    column, which centres to exactly zero - are held out of the fit and
    returned as exact zeros, since a zero column carries no signal and would
    otherwise be selected arbitrarily.

    Parameters
    ----------
    Theta       : ndarray, shape (K, m) - library matrix including a constant column
    ut          : ndarray, shape (K,) - target vector
    algorithm   : callable ``(Theta, ut) -> w`` - sparse solver, e.g.
                  ``lambda X, y: TrainSTRidge(X, y, lam, d_tol)`` or
                  ``lambda X, y: TrainRFE(X, y, lam=lam)`` (RFERidge).  Solvers
                  returning a tuple (such as ``TrainMSTLS``) are unwrapped by
                  taking element 0.  Note the ``(design, target)`` argument
                  order, matching the rest of this module.
    tol         : float - relative threshold for the intercept; ``w0`` is zeroed
                  when ``|w0| < tol * mean(|ut|)``.  0 (the default) disables the
                  threshold, so that calling this function never silently
                  suppresses a term.  On the weak-form AC/SH/PME/CH data in this
                  repo - none of which carry a genuine constant - the spurious
                  intercept left by STRidge runs at 1e-3 to 2e-2 of
                  ``mean(|ut|)`` across 0-5 % noise, so ``tol=0.05`` removes it
                  with an order of magnitude to spare; that is the recommended
                  starting value.  Note that ``w0`` absorbs the mean of whatever
                  the sparse fit failed to explain, so a large ``|w0|`` is
                  usually a symptom of a bad fit (missing terms) rather than a
                  real constant - raising ``tol`` to hide it is the wrong fix.
    descs       : sequence of str or None - column descriptions, used to locate
                  the constant column by name
    const_index : int or None - explicit index of the constant column
    const_desc  : str - description marking the constant column (default ``"1"``)
    allow_missing : bool - if True, centre even when there is no constant column
                  and return the intercept separately instead of raising
    return_intercept : bool - also return the scalar intercept
    zero_tol    : float - relative norm below which a centred column is treated
                  as degenerate
    cache       : SolverCache or None - reuse the centred system across calls
                  on the same (Theta, ut).  ``algorithm`` then receives the
                  same centred arrays each time, so a cache passed to the
                  solver inside ``algorithm`` is reused as well; results are
                  unchanged

    Returns
    -------
    w : ndarray, shape (m,) - coefficient vector aligned with ``Theta``'s columns.
        The constant column's slot holds ``w0 / c``, where ``c`` is that
        column's value (1 pointwise; the test-function mass ⟨ψ, 1⟩ in the weak
        form), so ``Theta @ w`` is the fitted model in every formulation.
    w0 : float - the intercept in units of ``ut`` (the additive offset), only
        when ``return_intercept`` is True.  Equals the slot when ``c == 1``.  When
        there is no constant column (``allow_missing=True``) the intercept is
        *not* representable in ``w`` and is only available through this.

    Raises
    ------
    ValueError
        If no constant column can be located and ``allow_missing`` is False.
    """
    pre = None
    if cache is not None:
        pre = cache.entry(
            (Theta, ut),
            ("centering", None if descs is None else tuple(descs), const_index,
             const_desc, allow_missing, zero_tol),
        )
    Theta = np.asarray(Theta, dtype=float)
    ut = np.asarray(ut, dtype=float).ravel()
    K, m = Theta.shape

    if pre is not None and "live" in pre:
        j_const, keep_idx, mu_y, mu_X, y1, live, X1_live = (
            pre[k] for k in ("j_const", "keep_idx", "mu_y", "mu_X", "y1", "live",
                             "X1_live"))
    else:
        j_const, keep_idx, mu_y, mu_X, y1, live, X1_live = _centre_system(
            Theta, ut, K, m, descs, const_index, const_desc, allow_missing,
            zero_tol)
        if pre is not None:
            pre.update(j_const=j_const, keep_idx=keep_idx, mu_y=mu_y, mu_X=mu_X,
                       y1=y1, live=live, X1_live=X1_live)

    w_sub = np.zeros(keep_idx.size)
    if np.any(live):
        w_fit = algorithm(X1_live, y1)
        if isinstance(w_fit, tuple):
            w_fit = w_fit[0]
        w_sub[live] = np.asarray(w_fit, dtype=float).ravel()

    # Recover the intercept in closed form and threshold it
    w0 = mu_y - float(w_sub @ mu_X)
    if tol > 0.0:
        y_scale = float(np.mean(np.abs(ut)))
        if abs(w0) < tol * max(y_scale, 1e-300):
            w0 = 0.0

    # Reassemble a full-length vector aligned with Theta's columns.  The
    # constant column need not be all ones: in the weak form it is ⟨ψ, 1⟩, the
    # test function's mass.  Its coefficient is w0 divided by that value, so
    # that Theta @ w reproduces the fit and the slot holds the PDE's constant.
    w = np.zeros(m)
    w[keep_idx] = w_sub
    if j_const is not None:
        c_val = float(Theta[:, j_const].mean())
        w[j_const] = w0 / c_val if c_val != 0.0 else 0.0

    if return_intercept:
        return w, w0
    return w


def _centre_system(Theta, ut, K, m, descs, const_index, const_desc, allow_missing,
                   zero_tol):
    """The data-only part of ``centering``: locate the constant column, centre
    the target and the other columns, and screen out degenerate ones.

    Returns ``(j_const, keep_idx, mu_y, mu_X, y1, live, X1[:, live])``.
    """
    if ut.shape[0] != K:
        raise ValueError(
            f"ut has length {ut.shape[0]} but Theta has {K} rows."
        )

    j_const = _resolve_const_index(Theta, descs, const_index, const_desc)
    if j_const is None and not allow_missing:
        raise ValueError(
            "No constant column found in Theta. Pass descs= or const_index= to "
            "identify it, or allow_missing=True to centre anyway and receive "
            "the intercept separately (it has no slot in the returned vector)."
        )

    # Non-constant columns, in their original order
    keep = np.ones(m, dtype=bool)
    if j_const is not None:
        keep[j_const] = False
    keep_idx = np.where(keep)[0]
    X00 = Theta[:, keep_idx]

    # Centre target and columns
    mu_y = float(ut.mean())
    mu_X = X00.mean(axis=0)
    y1 = ut - mu_y
    X1 = X00 - mu_X[np.newaxis, :]

    # Screen columns that centring made degenerate (any *other* constant column
    # centres to exactly zero).  They are excluded from the fit and returned as
    # exact zeros rather than left to absorb noise.
    col_norms = np.linalg.norm(X1, axis=0)
    ref = np.median(col_norms) if col_norms.size else 0.0
    if ref <= 0.0:
        ref = float(np.max(col_norms)) if col_norms.size else 0.0
    live = col_norms > zero_tol * max(ref, 1e-300)

    return j_const, keep_idx, mu_y, mu_X, y1, live, X1[:, live]
