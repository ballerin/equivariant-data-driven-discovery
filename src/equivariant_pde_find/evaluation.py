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
from collections import defaultdict

from .library import build_equivariant_library, build_polynomial_jet_library
from .simulations import make_regression_data


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
    seed       : int — forwarded to make_regression_data for spatial subsampling RNG

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


def build_eval_theta_jet(snapshots, times, Lx, Ly, skip_fraction=0.0, seed=0):
    """Pre-build the polynomial-jet Theta and Ft for unified evaluation.

    All libraries are evaluated in the common polynomial-jet space
    (monomials in jet variables f, f_x, f_y, f_xx, …) so that fp/fn/E
    are comparable across libraries.
    """
    f_snaps, ft_snaps, _, _ = make_regression_data(
        snapshots, times, n_per_snap=200,
        skip_fraction=skip_fraction, seed=seed,
    )
    rows_theta, rows_ft = [], []
    for s in range(len(f_snaps)):
        T, _ = build_polynomial_jet_library(
            f_snaps[s], Lx, Ly, max_order=4, max_degree=3,
        )
        rows_theta.append(T)
        rows_ft.append(ft_snaps[s].ravel())
    return np.vstack(rows_theta), np.concatenate(rows_ft)


# ── WSINDy → jet analytical expansion ───────────────────────────────────

def _apply_deriv_x(terms):
    """Apply ∂_x to a dict of {frozenset((a,b),...): (count, coeff)}.
    Each key is a frozenset of (a,b) tuples representing copies of f
    with given derivative orders; count is how many copies; coeff is
    the scalar multiplier."""
    result = defaultdict(float)
    for key, (count, coeff) in terms.items():
        key_list = list(key)
        for i, (a, b) in enumerate(key_list):
            new_key = list(key_list)
            new_key[i] = (a + 1, b)
            new_key.sort()
            result[tuple(new_key)] += coeff
    return {k: (len(k), v) for k, v in result.items()}


def _apply_deriv_y(terms):
    """Apply ∂_y."""
    result = defaultdict(float)
    for key, (count, coeff) in terms.items():
        key_list = list(key)
        for i, (a, b) in enumerate(key_list):
            new_key = list(key_list)
            new_key[i] = (a, b + 1)
            new_key.sort()
            result[tuple(new_key)] += coeff
    return {k: (len(k), v) for k, v in result.items()}


def _expand_fq_deriv(a, b, q):
    """Expand ∂_x^a ∂_y^b (f^q) to jet monomials.

    Returns dict {((a₁,b₁), ..., (a_q,b_q)): coefficient}
    where each (a_i,b_i) is a derivative multi-index for one copy of f.
    The jet monomial is ∏_i ∂_x^{a_i} ∂_y^{b_i} f.
    """
    # Start with q copies of undifferentiated f: ((0,0), (0,0), ..., (0,0))
    key = tuple((0, 0) for _ in range(q))
    terms = {key: (q, 1.0)}  # key → (count, coeff)

    for _ in range(a):
        terms = _apply_deriv_x(terms)
    for _ in range(b):
        terms = _apply_deriv_y(terms)

    return terms


def _terms_to_jet_dict(terms, jet_names):
    """Convert expansion terms to a dict of {jet_description: coefficient}.

    Each term key is a tuple of (a,b) pairs.  The jet monomial
    description is determined by the multiset of derivative orders.
    """
    result = defaultdict(float)
    for key, (_count, coeff) in terms.items():
        # Build jet variable exponents
        exponents = defaultdict(int)
        for a, b in key:
            if a == 0 and b == 0:
                exponents[(0, 0)] += 1  # bare f
            else:
                exponents[(a, b)] += 1
        # Build canonical description (same ordering as build_polynomial_jet_library)
        parts = []
        for (a, b), exp in sorted(exponents.items(), key=lambda kv: (kv[0][0] + kv[0][1], kv[0][0])):
            if exp == 0:
                continue
            vname = jet_names.get((a, b), "")
            if not vname:
                if a == 0 and b == 0:
                    vname = "f"
                else:
                    vname = "f_" + "x" * a + "y" * b
            if exp == 1:
                parts.append(vname)
            else:
                sup_map = {2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}
                parts.append(f"{vname}{sup_map.get(exp, f'^{exp}')}")
        desc = "·".join(parts) if parts else "1"
        result[desc] += coeff
    return dict(result)


def build_wsindy_to_jet_map(desc_ws, desc_jet):
    """Build a matrix M such that w_jet = M @ w_wsindy.

    Each WSINDy column ∂_x^a ∂_y^b (f^q) is analytically expanded to
    jet monomials using the general Leibniz rule.
    """
    # Build jet variable names
    jet_names = {}
    for total in range(5):  # max_order=4
        for a in range(total + 1):
            b = total - a
            if a == 0 and b == 0:
                jet_names[(a, b)] = "f"
            else:
                jet_names[(a, b)] = "f_" + "x" * a + "y" * b

    M = np.zeros((len(desc_jet), len(desc_ws)))
    for j, desc in enumerate(desc_ws):
        if desc == "1":
            M[desc_jet.index("1"), j] = 1.0
            continue
        # Parse WSINDy description: "f_{xx}", "(f²)_{xx}", "f²", etc.
        sup_map_rev = {"²": 2, "³": 3, "⁴": 4, "⁵": 5, "⁶": 6, "⁷": 7, "⁸": 8, "⁹": 9}
        if desc.startswith("(f"):
            # (f^q)_{x...y...}
            q_str = ""
            idx = 2
            while idx < len(desc) and desc[idx] in sup_map_rev:
                q_str += desc[idx]
                idx += 1
            q = sup_map_rev.get(q_str, 2) if q_str else 2
            # Parse derivative: )_{xxy...}
            rest = desc[desc.index(")_{") + 3: -1] if desc.endswith("}") else desc[desc.index(")_{") + 3:]
            a = rest.count("x")
            b = rest.count("y")
        elif desc.startswith("f_{"):
            # f_{xxy...}  (q=1)
            q = 1
            rest = desc[3:-1] if desc.endswith("}") else desc[3:]
            a = rest.count("x")
            b = rest.count("y")
        elif desc in ("f", "f²", "f³"):
            # Pure monomial: look up directly
            M[desc_jet.index(desc), j] = 1.0
            continue
        else:
            continue  # shouldn't happen

        terms = _expand_fq_deriv(a, b, q)
        jet_dict = _terms_to_jet_dict(terms, jet_names)
        for jet_desc, coeff in jet_dict.items():
            try:
                i = desc_jet.index(jet_desc)
            except ValueError:
                continue  # jet monomial not in library (e.g., order > 4)
            M[i, j] += coeff

    return M


def map_wsindy_to_jet(w_ws, M_ws_to_jet):
    """Map WSINDy coefficients to jet coefficients."""
    return M_ws_to_jet @ w_ws


def map_standard_to_jet(w_std, desc_std, desc_jet):
    """Map standard-library coefficients to jet coefficients.

    Standard columns are of the form f^p · ∂_x^a ∂_y^b f which are
    individual jet monomials; mapping is 1:1 index lookup.
    Columns beyond the jet library's degree/order limits are silently
    skipped (they cannot contribute to the evaluation).
    """
    w_jet = np.zeros(len(desc_jet))
    for j, desc in enumerate(desc_std):
        if abs(w_std[j]) < 1e-14:
            continue
        if desc in desc_jet:
            w_jet[desc_jet.index(desc)] += w_std[j]
    return w_jet


def map_equivariant_to_jet_from_w(w_eq, desc_eq, desc_jet):
    """Map equivariant coefficients to jet coefficients using the
    analytical decomposition in _EQUIV_TO_JET_DECOMPOSITION."""
    w_jet = np.zeros(len(desc_jet))
    for j, (name, coeff) in enumerate(zip(desc_eq, w_eq)):
        if abs(coeff) < 1e-14:
            continue
        if name not in _EQUIV_TO_JET_DECOMPOSITION:
            continue
        for jet_name, mult in _EQUIV_TO_JET_DECOMPOSITION[name]:
            w_jet[desc_jet.index(jet_name)] += coeff * mult
    return w_jet


def map_mod_multigraph_to_jet_from_w(w_mod, desc_mod, desc_jet):
    """Map modified multigraph coefficients to jet coefficients using the
    analytical decomposition in _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION."""
    w_jet = np.zeros(len(desc_jet))
    for j, (name, coeff) in enumerate(zip(desc_mod, w_mod)):
        if abs(coeff) < 1e-14:
            continue
        if name not in _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION:
            continue
        for jet_name, mult in _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION[name]:
            w_jet[desc_jet.index(jet_name)] += coeff * mult
    return w_jet


def evaluate_unified(w, library, desc_lib, desc_jet, Theta_jet, Ft,
                     true_jet, M_ws_to_jet=None):
    """Evaluate a recovery from any library against the polynomial-jet truth.

    Maps the recovered coefficients to the jet basis, then computes
    fp, fn, E in the common jet space.
    """
    if library in ("equivariant", "multigraph"):
        w_jet = map_equivariant_to_jet_from_w(w, desc_lib, desc_jet)
    elif library == "mod_multigraph":
        w_jet = map_mod_multigraph_to_jet_from_w(w, desc_lib, desc_jet)
    elif library == "standard":
        w_jet = map_standard_to_jet(w, desc_lib, desc_jet)
    elif library in ("wsindy", "wsindy_proj", "mwsindy") and M_ws_to_jet is not None:
        w_jet = map_wsindy_to_jet(w, M_ws_to_jet)
    elif library == "polynomial_jet":
        w_jet = w
    else:
        # Fallback: use library's own Theta
        return evaluate_from_theta(w, true_jet, Theta_jet, Ft)

    return evaluate_from_theta(w_jet, true_jet, Theta_jet, Ft)


def evaluate_mwsindy(w_ws, desc_ws, desc_mg, desc_jet, M_ws_to_jet,
                     Theta_jet, Ft, true_jet, G_mg, b_weak):
    """MWSINDy: WSINDy model selection + multigraph fine-tuning.

    Pipeline
    --------
    1. WSINDy regression (already done — w_ws is provided)
    2. Map WSINDy coefficients → jet → identify active multigraph operators
       (all jet components must survive thresholding)
    3. Fine-tune on multigraph submatrix (only active columns)
    4. Evaluate fp/fn/E in multigraph space against ground truth

    Parameters
    ----------
    w_ws       : ndarray (K_ws,) — sparse WSINDy coefficient vector
    desc_ws    : list of str — WSINDy column descriptions
    desc_mg    : list of str — multigraph operator descriptions
    desc_jet   : list of str — polynomial jet descriptions
    M_ws_to_jet : ndarray (K_jet, K_ws) — WSINDy → jet mapping matrix
    Theta_jet  : ndarray (n, K_jet) — pointwise jet evaluation matrix
    Ft         : ndarray (n,) — time derivative
    true_jet   : ndarray (K_jet,) — ground truth in jet space
    G_mg       : ndarray (n_weak, K_mg) — weak-form multigraph matrix
    b_weak     : ndarray (n_weak,) — weak-form RHS

    Returns
    -------
    dict with fp, fn, exact_recovery, n_active, E, max_coeff_err, n_samples
    """
    from .library import build_equivariant_library

    # ── Step 1: Map WSINDy → jet ────────────────────────────────────────
    w_jet = map_wsindy_to_jet(w_ws, M_ws_to_jet)
    nz_jet = np.abs(w_jet) > 1e-10

    # ── Step 2: Identify active multigraph operators ────────────────────
    # An operator is active iff ALL its jet components survive thresholding
    active_mg = []
    for j_mg, mg_name in enumerate(desc_mg):
        if mg_name not in _EQUIV_TO_JET_DECOMPOSITION:
            continue
        # Check all jet components
        all_active = True
        for jet_name, _mult in _EQUIV_TO_JET_DECOMPOSITION[mg_name]:
            j_jet = desc_jet.index(jet_name)
            if not nz_jet[j_jet]:
                all_active = False
                break
        if all_active:
            active_mg.append(j_mg)

    if not active_mg:
        # No multigraph operator fully recovered → hard failure.
        # The data does not support an O(2)-equivariant PDE.
        return dict(
            fp=0, fn=int(np.sum(np.abs(true_jet) > 1e-10)),
            exact_recovery=0, n_active=0, E=1.0,
            max_coeff_err=float(np.max(np.abs(true_jet[np.abs(true_jet) > 1e-10]))),
            n_samples=Theta_jet.shape[0],
        )

    # ── Step 3: Fine-tune on multigraph submatrix ───────────────────────
    G_sub = G_mg[:, active_mg]
    w_mg_sub, _, _, _ = np.linalg.lstsq(G_sub, b_weak, rcond=None)
    w_mg = np.zeros(len(desc_mg))
    for i, idx in enumerate(active_mg):
        w_mg[idx] = w_mg_sub[i]

    # ── Step 4: Map multigraph → jet and evaluate ────────────────────────
    w_jet_refined = map_equivariant_to_jet_from_w(w_mg, desc_mg, desc_jet)
    return evaluate_from_theta(w_jet_refined, true_jet, Theta_jet, Ft)


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
    true_std   : ndarray (K_std,) — ground truth in standard basis
    n_omitted  : int — number of true equivariant operators that could not
                 be decomposed (always-missing; should be added to fn).
    """
    true_std = np.zeros(len(desc_std), dtype=float)
    n_omitted = 0
    for i, (name, coeff) in enumerate(zip(desc_eq, true_eq)):
        if abs(coeff) < 1e-12:
            continue
        if name not in _EQUIV_TO_STD_DECOMPOSITION:
            n_omitted += 1  # operator not representable in standard library
            continue
        for std_name, mult in _EQUIV_TO_STD_DECOMPOSITION[name]:
            try:
                j = desc_std.index(std_name)
            except ValueError:
                raise KeyError(
                    f"Standard column {std_name!r} not found in library. "
                    f"Available: {desc_std[:10]}..."
                )
            true_std[j] += coeff * mult
    return true_std, n_omitted


def evaluate_standard_with_truth(w, true_eq, desc_eq, Theta_std, desc_std, Ft):
    """Evaluate a standard-library recovery with proper fp/fn via decomposition.

    Decomposes the equivariant ground truth into the standard library
    wherever possible, then computes fp, fn, exact_recovery, E, and
    max_coeff_err in the standard basis.  Non-representable equivariant
    operators are counted as always-missing (added to fn).

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
    n_samples
    """
    true_std, n_omitted = decompose_equivariant_to_standard(
        true_eq, desc_eq, desc_std,
    )
    result = evaluate_from_theta(w, true_std, Theta_std, Ft)
    if n_omitted > 0:
        result["fn"] = int(result["fn"]) + n_omitted
        result["exact_recovery"] = int(
            result["fp"] == 0 and result["fn"] == 0
        )
    return result


# ── Equivariant → polynomial jet library decomposition ────────────────────

# Maps equivariant operator name → list of (jet_monomial_name, coefficient).
# The polynomial jet library (monomials in jet variables f, f_x, f_y, f_xx, …)
# is spanning for all O(2)-equivariant operators up to max_order=4, max_degree=3.
# Jet monomial descriptions use canonical sort order: by (total_order, x-count).
_EQUIV_TO_JET_DECOMPOSITION = {
    "1":                [("1", 1.0)],
    "f":                [("f", 1.0)],
    "f²":               [("f²", 1.0)],
    "f³":               [("f³", 1.0)],
    "Δf":               [("f_xx", 1.0), ("f_yy", 1.0)],
    "f·Δf":             [("f·f_xx", 1.0), ("f·f_yy", 1.0)],
    "‖∇f‖²":           [("f_x²", 1.0), ("f_y²", 1.0)],
    "f²·Δf":            [("f²·f_xx", 1.0), ("f²·f_yy", 1.0)],
    "f·‖∇f‖²":         [("f·f_x²", 1.0), ("f·f_y²", 1.0)],
    "Δ²f":              [("f_xxxx", 1.0), ("f_xxyy", 2.0), ("f_yyyy", 1.0)],
    "f·Δ²f":            [("f·f_xxxx", 1.0), ("f·f_xxyy", 2.0),
                         ("f·f_yyyy", 1.0)],
    "Δf·Δf":            [("f_xx²", 1.0), ("f_yy·f_xx", 2.0), ("f_yy²", 1.0)],
    "‖∇²f‖²":          [("f_xx²", 1.0), ("f_xy²", 2.0), ("f_yy²", 1.0)],
    "⟨∇Δf,∇f⟩":        [("f_x·f_xxx", 1.0), ("f_x·f_xyy", 1.0),
                         ("f_y·f_xxy", 1.0), ("f_y·f_yyy", 1.0)],
    "f²·Δ²f":           [("f²·f_xxxx", 1.0), ("f²·f_xxyy", 2.0),
                         ("f²·f_yyyy", 1.0)],
    "f·Δf·Δf":          [("f·f_xx²", 1.0), ("f·f_yy·f_xx", 2.0),
                         ("f·f_yy²", 1.0)],
    "f·‖∇²f‖²":        [("f·f_xx²", 1.0), ("f·f_xy²", 2.0),
                         ("f·f_yy²", 1.0)],
    "f·⟨∇Δf,∇f⟩":      [("f·f_x·f_xxx", 1.0), ("f·f_x·f_xyy", 1.0),
                         ("f·f_y·f_xxy", 1.0), ("f·f_y·f_yyy", 1.0)],
    "‖∇f‖²·Δf":        [("f_x²·f_xx", 1.0), ("f_x²·f_yy", 1.0),
                         ("f_y²·f_xx", 1.0), ("f_y²·f_yy", 1.0)],
    "∇²f(∇f,∇f)":      [("f_x²·f_xx", 1.0), ("f_y·f_x·f_xy", 2.0),
                         ("f_y²·f_yy", 1.0)],
}


# ── Modified Multigraph → polynomial jet library decomposition ────────────

# Maps modified multigraph operator name → list of (jet_monomial, coefficient).
# The 14 unchanged operators reuse entries from _EQUIV_TO_JET_DECOMPOSITION;
# the 6 new operators use corrected formulas from Lemma 1 (see lemma1-errata.md).
_MOD_MULTIGRAPH_TO_JET_DECOMPOSITION = {
    # --- Unchanged operators (identical to _EQUIV_TO_JET_DECOMPOSITION) ---
    "1":                [("1", 1.0)],
    "f":                [("f", 1.0)],
    "f²":               [("f²", 1.0)],
    "f³":               [("f³", 1.0)],
    "Δf":               [("f_xx", 1.0), ("f_yy", 1.0)],
    "‖∇f‖²":           [("f_x²", 1.0), ("f_y²", 1.0)],
    "f·‖∇f‖²":         [("f·f_x²", 1.0), ("f·f_y²", 1.0)],
    "Δ²f":              [("f_xxxx", 1.0), ("f_xxyy", 2.0), ("f_yyyy", 1.0)],
    "‖∇²f‖²":          [("f_xx²", 1.0), ("f_xy²", 2.0), ("f_yy²", 1.0)],
    "Δf·Δf":            [("f_xx²", 1.0), ("f_yy·f_xx", 2.0), ("f_yy²", 1.0)],
    "f·‖∇²f‖²":        [("f·f_xx²", 1.0), ("f·f_xy²", 2.0),
                         ("f·f_yy²", 1.0)],
    "f·Δf·Δf":          [("f·f_xx²", 1.0), ("f·f_yy·f_xx", 2.0),
                         ("f·f_yy²", 1.0)],
    "∇²f(∇f,∇f)":      [("f_x²·f_xx", 1.0), ("f_y·f_x·f_xy", 2.0),
                         ("f_y²·f_yy", 1.0)],
    "‖∇f‖²·Δf":        [("f_x²·f_xx", 1.0), ("f_x²·f_yy", 1.0),
                         ("f_y²·f_xx", 1.0), ("f_y²·f_yy", 1.0)],

    # --- Modified operators (6 replacements, corrected coefficients) ---
    # Δf² = 2f·Δf + 2‖∇f‖²
    "Δf²":              [("f·f_xx", 2.0), ("f·f_yy", 2.0),
                         ("f_x²", 2.0), ("f_y²", 2.0)],

    # Δf³ = 3f²·Δf + 6f·‖∇f‖²
    "Δf³":              [("f²·f_xx", 3.0), ("f²·f_yy", 3.0),
                         ("f·f_x²", 6.0), ("f·f_y²", 6.0)],

    # Δ‖∇f‖² = 2⟨∇Δf,∇f⟩ + 2‖∇²f‖²
    "Δ‖∇f‖²":          [("f_x·f_xxx", 2.0), ("f_x·f_xyy", 2.0),
                         ("f_y·f_xxy", 2.0), ("f_y·f_yyy", 2.0),
                         ("f_xx²", 2.0), ("f_xy²", 4.0), ("f_yy²", 2.0)],

    # Δ(f·‖∇f‖²) = 2f·⟨∇Δf,∇f⟩ + 2f·‖∇²f‖² + ‖∇f‖²·Δf + 4∇²f(∇f,∇f)
    "Δ(f·‖∇f‖²)":      [("f·f_x·f_xxx", 2.0), ("f·f_x·f_xyy", 2.0),
                         ("f·f_y·f_xxy", 2.0), ("f·f_y·f_yyy", 2.0),
                         ("f·f_xx²", 2.0), ("f·f_xy²", 4.0), ("f·f_yy²", 2.0),
                         ("f_x²·f_xx", 5.0), ("f_x²·f_yy", 1.0),
                         ("f_y²·f_xx", 1.0), ("f_y²·f_yy", 5.0),
                         ("f_y·f_x·f_xy", 8.0)],

    # Δ²f² = 2f·Δ²f + 8⟨∇Δf,∇f⟩ + 4‖∇²f‖² + 2(Δf)²  (CORRECTED: 4‖∇²f‖²)
    "Δ²f²":             [("f·f_xxxx", 2.0), ("f·f_xxyy", 4.0),
                         ("f·f_yyyy", 2.0),
                         ("f_x·f_xxx", 8.0), ("f_x·f_xyy", 8.0),
                         ("f_y·f_xxy", 8.0), ("f_y·f_yyy", 8.0),
                         ("f_xx²", 6.0), ("f_xy²", 8.0), ("f_yy²", 6.0),
                         ("f_yy·f_xx", 4.0)],

    # Δ²f³ = 3f²·Δ²f + 24f·⟨∇Δf,∇f⟩ + 12f·‖∇²f‖² + 6f·(Δf)²
    #        + 24∇²f(∇f,∇f) + 12‖∇f‖²·Δf  (CORRECTED: 24f⟨∇Δf,∇f⟩)
    "Δ²f³":             [("f²·f_xxxx", 3.0), ("f²·f_xxyy", 6.0),
                         ("f²·f_yyyy", 3.0),
                         ("f·f_x·f_xxx", 24.0), ("f·f_x·f_xyy", 24.0),
                         ("f·f_y·f_xxy", 24.0), ("f·f_y·f_yyy", 24.0),
                         ("f·f_xx²", 18.0), ("f·f_xy²", 24.0),
                         ("f·f_yy²", 18.0), ("f·f_yy·f_xx", 12.0),
                         ("f_x²·f_xx", 36.0), ("f_y·f_x·f_xy", 48.0),
                         ("f_y²·f_yy", 36.0),
                         ("f_x²·f_yy", 12.0), ("f_y²·f_xx", 12.0)],
}


def decompose_equivariant_to_jet(true_eq, desc_eq, desc_jet):
    """Decompose equivariant ground truth into polynomial-jet coefficients.

    Every O(2)-equivariant operator up to order 4, degree 3 can be
    expressed as a linear combination of jet monomials (unlike the
    standard library, which omits products of derivatives like ‖∇f‖²).

    Parameters
    ----------
    true_eq  : ndarray (K_eq,) — ground-truth coefficients in equivariant basis
    desc_eq  : list of K_eq str — equivariant operator descriptions
    desc_jet : list of K_jet str — jet monomial descriptions

    Returns
    -------
    true_jet : ndarray (K_jet,) — ground truth in jet basis
    """
    true_jet = np.zeros(len(desc_jet), dtype=float)
    for i, (name, coeff) in enumerate(zip(desc_eq, true_eq)):
        if abs(coeff) < 1e-12:
            continue
        if name not in _EQUIV_TO_JET_DECOMPOSITION:
            raise KeyError(
                f"Equivariant operator {name!r} has no jet decomposition. "
                f"Available: {list(_EQUIV_TO_JET_DECOMPOSITION)}"
            )
        for jet_name, mult in _EQUIV_TO_JET_DECOMPOSITION[name]:
            try:
                j = desc_jet.index(jet_name)
            except ValueError:
                raise KeyError(
                    f"Jet monomial {jet_name!r} not found in library. "
                    f"Available (first 10): {desc_jet[:10]}..."
                )
            true_jet[j] += coeff * mult
    return true_jet


def decompose_mod_multigraph_to_jet(true_mod, desc_mod, desc_jet):
    """Decompose mod. multigraph ground truth into polynomial-jet coefficients.

    Parameters
    ----------
    true_mod : ndarray (K_mod,)
    desc_mod : list of K_mod str
    desc_jet : list of K_jet str

    Returns
    -------
    true_jet : ndarray (K_jet,)
    """
    true_jet = np.zeros(len(desc_jet), dtype=float)
    for i, (name, coeff) in enumerate(zip(desc_mod, true_mod)):
        if abs(coeff) < 1e-12:
            continue
        if name not in _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION:
            raise KeyError(
                f"Mod. multigraph operator {name!r} has no jet decomposition. "
                f"Available: {list(_MOD_MULTIGRAPH_TO_JET_DECOMPOSITION)}"
            )
        for jet_name, mult in _MOD_MULTIGRAPH_TO_JET_DECOMPOSITION[name]:
            try:
                j = desc_jet.index(jet_name)
            except ValueError:
                raise KeyError(
                    f"Jet monomial {jet_name!r} not found in library. "
                    f"Available (first 10): {desc_jet[:10]}..."
                )
            true_jet[j] += coeff * mult
    return true_jet


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
    true_jet = decompose_equivariant_to_jet(true_eq, desc_eq, desc_jet)
    return evaluate_from_theta(w, true_jet, Theta_jet, Ft)


# ── Equivariant → WSINDy library decomposition ─────────────────────────

# Maps equivariant operator name → list of (wsindy_column_name, coefficient).
# The WSINDy library ∂^α(f^q) can only represent pure monomials (1, f, f², f³)
# and pure linear derivatives of f (Δf, Δ²f, and their components).
# Operators involving products of derivatives (‖∇f‖², f·Δf, etc.) are NOT in
# the column span — the chain rule for q ≥ 2 inextricably mixes product-of-
# derivative terms with f×derivative terms.
_EQUIV_TO_WSINDY_DECOMPOSITION = {
    "1":    [("1", 1.0)],
    "f":    [("f", 1.0)],
    "f²":   [("f²", 1.0)],
    "f³":   [("f³", 1.0)],
    "Δf":   [("f_{xx}", 1.0), ("f_{yy}", 1.0)],
    "Δ²f":  [("f_{xxxx}", 1.0), ("f_{xxyy}", 2.0), ("f_{yyyy}", 1.0)],
}


def decompose_equivariant_to_wsindy(true_eq, desc_eq, desc_ws):
    """Decompose equivariant ground truth into WSINDy coefficients.

    Only pure monomials and linear derivatives of f are representable
    in the WSINDy basis ∂^α(f^q).  Operators involving products of
    derivatives or products of f with derivatives are silently omitted
    (they are not in the column span).

    Parameters
    ----------
    true_eq  : ndarray (K_eq,) — ground-truth coefficients in equivariant basis
    desc_eq  : list of K_eq str — equivariant operator descriptions
    desc_ws  : list of K_ws str — WSINDy column descriptions

    Returns
    -------
    true_ws : ndarray (K_ws,) — ground truth in WSINDy basis
    """
    true_ws = np.zeros(len(desc_ws), dtype=float)
    for i, (name, coeff) in enumerate(zip(desc_eq, true_eq)):
        if abs(coeff) < 1e-12:
            continue
        if name not in _EQUIV_TO_WSINDY_DECOMPOSITION:
            continue  # operator not representable in WSINDy library
        for ws_name, mult in _EQUIV_TO_WSINDY_DECOMPOSITION[name]:
            try:
                j = desc_ws.index(ws_name)
            except ValueError:
                raise KeyError(
                    f"WSINDy column {ws_name!r} not found in library. "
                    f"Available (first 10): {desc_ws[:10]}..."
                )
            true_ws[j] += coeff * mult
    return true_ws


def evaluate_wsindy_with_truth(w, true_eq, desc_eq, desc_ws, Theta_ws, Ft):
    """Evaluate a WSINDy-library recovery via decomposition.

    Decomposes the equivariant ground truth into the WSINDy basis
    wherever possible, then computes fp, fn, exact_recovery, E, and
    max_coeff_err.  Operators not representable as ∂^α(f^q) are
    silently omitted from the ground truth.

    Parameters
    ----------
    w         : ndarray (K_ws,)
    true_eq   : ndarray (K_eq,)
    desc_eq   : list of K_eq str
    desc_ws   : list of K_ws str
    Theta_ws  : ndarray (n, K_ws)
    Ft        : ndarray (n,)

    Returns
    -------
    dict with keys: fp, fn, exact_recovery, n_active, E, max_coeff_err,
    n_samples
    """
    true_ws = decompose_equivariant_to_wsindy(true_eq, desc_eq, desc_ws)
    return evaluate_from_theta(w, true_ws, Theta_ws, Ft)


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
    seed       : int — forwarded to make_regression_data for spatial subsampling RNG

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
