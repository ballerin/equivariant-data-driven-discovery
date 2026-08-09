"""
Collinearity diagnostics for design matrices.

Provides per-column and global diagnostics to quantify multicollinearity
in PDE regression libraries:

  pairwise_max_correlation  — max |Pearson r| among all column pairs
  vif_factors               — per-column Variance Inflation Factors
  condition_numbers         — raw and column-normalised 2-norm condition numbers
  collinearity_report       — combined summary dict

VIF formula
-----------
Uses the no-intercept regression formula:

    R²_j = 1 - ‖(I - H_{-j}) N_j‖² / ‖N_j - mean(N_j)‖²
    VIF_j = 1 / (1 - R²_j)

where H_{-j} = N_{-j} (N_{-j}^T N_{-j})^{-1} N_{-j}^T is the hat matrix
for the regression of column j on all other columns (no intercept).
Constant columns (std ≈ 0) are skipped with VIF = NaN.
"""

import numpy as np


def pairwise_max_correlation(Theta):
    """Maximum absolute Pearson correlation among all column pairs.

    Parameters
    ----------
    Theta : ndarray of shape (n_rows, n_cols)
        Design matrix.

    Returns
    -------
    float
        max_{i ≠ j} |corr(Theta[:, i], Theta[:, j])|.
        Returns 0.0 if n_cols < 2.
    """
    n_rows, n_cols = Theta.shape
    if n_cols < 2:
        return 0.0

    # Mean-centre columns for correlation
    Theta_c = Theta - Theta.mean(axis=0, keepdims=True)
    # Column norms
    col_norms = np.sqrt(np.sum(Theta_c ** 2, axis=0))
    # Guard against zero-variance columns
    safe_norms = np.where(col_norms < 1e-14, 1.0, col_norms)
    Theta_n = Theta_c / safe_norms

    # Correlation matrix: R = Theta_n.T @ Theta_n
    # R[i,j] = corr(col_i, col_j), but zero-variance columns produce NaNs
    R = Theta_n.T @ Theta_n
    # Mask out diagonal
    n = R.shape[0]
    mask = ~np.eye(n, dtype=bool)
    return float(np.nanmax(np.abs(R[mask])))


def vif_factors(Theta, descriptions=None, tol=1e-14, vif_cap=1e10):
    """Per-column Variance Inflation Factors (no-intercept formula).

    For each column j, computes VIF_j = 1 / (1 - R²_j) where R²_j is
    the R² from an OLS regression of column j on all other columns
    (without intercept).

    Uses the precomputed Gram matrix G = ΘᵀΘ to avoid repeated
    large matrix multiplications, reducing the per-column solve from
    O(n_rows × K²) to O(K³).

    Constant columns are skipped and receive VIF = NaN.  Constancy is
    tested as ``ptp(col) <= tol * max(|col|, 1)``, i.e. relative to the
    column's own magnitude, so the result does not depend on how many
    rows the design matrix happens to have.

    Parameters
    ----------
    Theta : ndarray of shape (n_rows, n_cols)
        Design matrix.
    descriptions : list of str or None
        Column descriptions (used only for skip-mask labels).
    tol : float
        Relative threshold for detecting constant columns (see above).
    vif_cap : float
        Maximum VIF value (capped to avoid inf when R² ≈ 1).

    Returns
    -------
    vif : ndarray of shape (n_cols,)
        VIF_j for each column.  NaN for skipped (constant) columns.
    skipped : ndarray of bool
        True where the column was skipped (constant).
    """
    n_rows, n_cols = Theta.shape
    if n_cols < 2:
        return np.full(n_cols, np.nan), np.zeros(n_cols, dtype=bool)

    # Detect constant columns.  Peak-to-peak is exact (comparisons only),
    # unlike np.std, whose cancellation error on an exactly-constant column
    # grows with n_rows and can exceed an absolute threshold on large data.
    # The tolerance is relative to each column's magnitude so the test is
    # invariant to both column scale and row count.
    col_scale = np.maximum(np.abs(Theta).max(axis=0), 1.0)
    skipped = np.ptp(Theta, axis=0) <= tol * col_scale

    # Total sum of squares (centred) — denominator for R²
    Theta_c = Theta - Theta.mean(axis=0, keepdims=True)
    TSS = np.sum(Theta_c ** 2, axis=0)

    # Precompute Gram matrix: G = ΘᵀΘ  (K×K)
    G = Theta.T @ Theta  # (n_cols, n_cols)

    vif = np.full(n_cols, np.nan)

    for j in range(n_cols):
        if skipped[j]:
            continue

        # Safety net for columns that vary but carry no usable signal.
        # Scaled by n_rows and column magnitude so it, too, is a relative
        # test rather than an absolute one.
        if TSS[j] <= (tol * col_scale[j]) ** 2 * n_rows:
            skipped[j] = True
            continue

        # Build G_{-j,-j} and G_{-j,j}
        keep = np.ones(n_cols, dtype=bool)
        keep[j] = False
        G_sub = G[np.ix_(keep, keep)]     # (K-1)×(K-1)
        g_rhs = G[keep, j]                 # (K-1,)

        # Solve G_{-j,-j} β̂ = G_{-j,j}
        try:
            beta = np.linalg.solve(G_sub, g_rhs)
        except np.linalg.LinAlgError:
            beta, _, _, _ = np.linalg.lstsq(G_sub, g_rhs, rcond=None)

        # RSS = G[j,j] - G_{-j,j}^T β̂
        RSS = G[j, j] - np.dot(g_rhs, beta)

        R2 = 1.0 - RSS / TSS[j]
        R2 = max(0.0, min(R2, 1.0 - 1e-15))

        if R2 >= 1.0 - 1e-15:
            vif[j] = float(vif_cap)
        else:
            vif[j] = min(1.0 / (1.0 - R2), float(vif_cap))

    return vif, skipped


def condition_numbers(Theta):
    """Raw and column-normalised 2-norm condition numbers.

    Parameters
    ----------
    Theta : ndarray of shape (n_rows, n_cols)

    Returns
    -------
    cond_raw : float
        κ(Θ) = σ_max / σ_min of the raw design matrix.
    cond_normalised : float
        κ after scaling each column to unit 2-norm.
    """
    if Theta.shape[1] < 2:
        return 1.0, 1.0

    # Raw condition number
    try:
        cond_raw = float(np.linalg.cond(Theta))
    except np.linalg.LinAlgError:
        cond_raw = np.inf

    # Column-normalised condition number
    col_norms = np.sqrt(np.sum(Theta ** 2, axis=0))
    safe_norms = np.where(col_norms < 1e-14, 1.0, col_norms)
    Theta_norm = Theta / safe_norms
    try:
        cond_norm = float(np.linalg.cond(Theta_norm))
    except np.linalg.LinAlgError:
        cond_norm = np.inf

    return cond_raw, cond_norm


def collinearity_report(Theta, descriptions=None):
    """Combined collinearity summary for a design matrix.

    Parameters
    ----------
    Theta : ndarray of shape (n_rows, n_cols)
    descriptions : list of str or None

    Returns
    -------
    dict with keys:
        n_rows, n_cols              — matrix dimensions
        max_abs_corr                — max |pairwise correlation|
        max_vif                     — max VIF (NaN if all skipped)
        mean_vif                    — mean VIF over non-skipped columns
        cond_raw, cond_normalised   — condition numbers
        skipped_columns             — list of description strings for skipped cols
        vif_per_column              — list of (description, vif) tuples (NaN for skipped)
    """
    n_rows, n_cols = Theta.shape

    max_corr = pairwise_max_correlation(Theta)
    vif, skipped = vif_factors(Theta, descriptions)
    cond_raw, cond_norm = condition_numbers(Theta)

    # Summary VIF stats (exclude NaN/skipped)
    vif_valid = vif[~skipped]
    max_vif = float(np.nanmax(vif_valid)) if len(vif_valid) > 0 else np.nan
    mean_vif = float(np.nanmean(vif_valid)) if len(vif_valid) > 0 else np.nan

    # Skipped column names
    if descriptions is not None:
        skipped_names = [descriptions[i] for i in range(n_cols) if skipped[i]]
        vif_list = [(descriptions[i], float(vif[i])) for i in range(n_cols)]
    else:
        skipped_names = [f"col_{i}" for i in range(n_cols) if skipped[i]]
        vif_list = [(f"col_{i}", float(vif[i])) for i in range(n_cols)]

    return {
        "n_rows": n_rows,
        "n_cols": n_cols,
        "max_abs_corr": max_corr,
        "max_vif": max_vif,
        "mean_vif": mean_vif,
        "cond_raw": cond_raw,
        "cond_normalised": cond_norm,
        "skipped_columns": skipped_names,
        "vif_per_column": vif_list,
    }
