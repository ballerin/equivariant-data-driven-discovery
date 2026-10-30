"""
Canonical ablation sweep: compare pointwise / space-weak / time-weak /
space-time weak formulations across equations, noise levels, libraries,
backends, and regression algorithms.

Loads pre-generated .npz data files (produced by generate_data.py) and runs
systematic hyperparameter sweeps for every method combination.  Results are
evaluated with the universal pointwise-Theta metric (evaluation.py).

Supports optional post-hoc noise application and spectral low-pass
preprocessing for additional ablation dimensions.

Usage
-----
    # Core equivariant vs standard comparison
    python scripts/sweep.py --methods pointwise --algorithms stridge

    # Noise robustness (space-weak vs pointwise)
    python scripts/sweep.py \\
        --equations swift_hohenberg ks_chaotic \\
        --noise 0.02 0.05 \\
        --methods pointwise space_weak_direct space_weak_ibp \\
        --libraries multigraph

    # With post-hoc noise and spectral low-pass
    python scripts/sweep.py --post-hoc-noise 0.02 0.05 --spectral-cutoff 3.0

    # Pointwise vs space-weak vs time-weak vs space-time weak, all with the
    # Messenger & Bortz (m, p, m_t, p_t) rule, with and without centering
    python scripts/sweep.py \\
        --methods pointwise space_weak_ibp time_weak spacetime_weak_ibp \\
        --test-fn-rule mb --tau-hat 2 --centering off on

    # Large grids: score through the QR factor of the evaluation matrix
    # (E_coeff to round-off, other columns identical, much faster)
    python scripts/sweep.py --data-dir data/large --eval-mode qr ...

    # Current standard (2026-09): space-time weak form with IBP, M&B rule,
    # periodic spatial boundary, 256² datasets
    python scripts/sweep.py --data-dir data/large \\
        --methods spacetime_weak_ibp --test-fn-rule mb --tau-hat 2 \\
        --space-boundary periodic --eval-mode qr

    # valid vs periodic spatial test functions
    python scripts/sweep.py --data-dir data/large \\
        --methods space_weak_ibp spacetime_weak_ibp --test-fn-rule mb \\
        --space-boundary valid periodic

    # Full sweep
    python scripts/sweep.py

Methods (``--methods``)
----------------------
=======================  =============  =============  =====================
method                   test fn space  test fn time   spatial derivatives
=======================  =============  =============  =====================
pointwise                -              -              backend
space_weak_direct        yes            -              backend, then ∗ φ_x
space_weak_ibp           yes            -              moved onto φ_x (IBP)
time_weak                -              yes            backend
spacetime_weak_direct    yes            yes            backend, then ∗ φ_x
spacetime_weak_ibp       yes            yes            moved onto φ_x (IBP)
=======================  =============  =============  =====================

The first three build per-snapshot rows with a finite-difference f_t, exactly
as before 2026-09 (when ``space_weak_direct`` / ``space_weak_ibp`` were called
``weak_direct`` / ``weak_ibp``; the old names are still accepted).  The last
three integrate against a temporal test function too
(``library/spacetime.py``); ``--time-ibp`` selects b = φ_t′ ∗ f (on) or
φ_t ∗ f_t with FD f_t (off), recorded in the ``time_ibp`` column.

Spatial boundary (``--space-boundary``)
---------------------------------------
``valid`` (default, Messenger & Bortz): test functions are centred only where
their support fits in the grid, losing 2m points per axis.  ``periodic``: the
support wraps around the (periodic) domain and every grid point is a centre;
see ``weak_form.use_weak_boundary``.  Applies to every method that is weak in
space.

Experiments from 2026-09 on use ``periodic`` for every equation: the simulated
domains are periodic, the derivative backends already assume it, KS-cell and
KS-wchaos need it (the rule asks m ≈ N/2), and one protocol then covers all
equations and grid sizes.  The method as run therefore assumes periodic data.
``valid`` remains the default only so that earlier sweeps reproduce unchanged.

Test-function parameters (``--test-fn-rule``)
---------------------------------------------
``legacy`` (default): spatial m from ``--space-m`` with p = SPACE_P_FIXED,
temporal m_t from ``--time-m`` with p_t from eq. (4.3) of Messenger & Bortz.
``mb``: every (m, p, m_t, p_t) chosen per dataset from the data by the
Messenger & Bortz rule (arXiv:2007.02848, eq. (4.3) and Appendix A), for all
weak methods alike.  Old invocations are unaffected by any new flag.
"""

import argparse
import csv
import datetime
import itertools
import re
import os
import platform
import subprocess
import sys
import time

import numpy as np

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        print(kwargs.get("desc", "") + f" ({len(list(iterable))} items)...")
        return iterable

# ── path setup ──────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives  # noqa: E402
from equivariant_pde_find.evaluation import (  # noqa: E402
    build_eval_theta_jet, evaluate_unified, FactoredEval,
)
from equivariant_pde_find.eval_setup import prepare_evaluation  # noqa: E402
from equivariant_pde_find.library.conversions import (  # noqa: E402
    build_wsindy_to_jet_map, _MULTIGRAPH_TO_JET_DECOMPOSITION,
    collect_needed_jet_monomials,
)
from equivariant_pde_find.equations import (  # noqa: E402
    find_op, true_allen_cahn, true_swift_hohenberg, true_ks,
    true_cahn_hilliard, true_pme, mask_pme, EQUATION_META,
    normalize_eq_name, CH_EPSILON, REPRESENTABLE_LIBRARIES,
)
from equivariant_pde_find.library import (  # noqa: E402
    build_multigraph_library, build_pdefind_library,
    build_polyjet_library, build_wsindy_library,
    build_mod_multigraph_library,
)
from equivariant_pde_find.utils import spectral_lowpass  # noqa: E402
from equivariant_pde_find.regression import (  # noqa: E402
    TrainSTRidge, TrainMSTLS, TrainRFE, TrainExhaustive, centering,
    SolverCache, _resolve_const_index,
)
from equivariant_pde_find.simulations import make_regression_data  # noqa: E402
from equivariant_pde_find.library import weak_form as wf  # noqa: E402
from equivariant_pde_find.library.spacetime import build_weak_system  # noqa: E402
from equivariant_pde_find.summary_metadata import (  # noqa: E402
    METHOD_ALIASES, canonical_method,
)

RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyperparameter grids ────────────────────────────────────────────────────
LAMS = [1e-6, 1e-4, 1e-2, 1.0]
RFE_LAMS = [0.0]  # ridge regularisation for RFE (0 = pure OLS)
MSTLS_N_LAMBDAS = 50
MSTLS_ALPHAS = [0.25, 0.5, 0.75]
MSTLS_MAXITS = 50
D_TOLS = [1e-3, 1e-2, 0.1, 1.0, 10.0]            # coarse grid (winners at 1.0, 10.0)
BRUTE_L0_PENALTIES = [1e-4, 1e-2, 1.0]  # mild, moderate, strong sparsity
BRUTE_MAX_K = 8  # cap for libraries > 20 cols; None for <= 20
BACKENDS = ["spectral", "polynomial"]
SPACE_M_GRID = [4, 8, 12]
# Test-function degree: we use a fixed p = SPACE_P_FIXED = 6 for both direct
# and IBP space-weak forms.
#
# This used to be `p = 6` for every m, which left phi^(4) large at the
# edge of the support and put a 3 % (m=12) to 84 % (m=8) systematic error
# into the fourth-order weak-IBP columns -- the leading terms of CH, SH and
# KS.
#
# The Messenger & Bortz (2021) Algorithm 1 prescribes p via
#   select_test_function_degree(m, max_deriv=4, tau=1e-10)
# which gives p=13 for m=12.  This satisfies a quadrature-accuracy decay
# condition (phi^(4) < 1e-10 at the support boundary).
#
# However, p=13 produces a narrow test function (effective width ~4 grid
# points) whose columns for different differential operators are nearly
# collinear locally, destroying sparse recovery.  At p=6 the IBP
# discretisation error is 0.27 % — still 20x below the 5 % observation
# noise floor — but the wider effective support (~6 grid points) gives
# each operator a more distinctive spatial signature, which is decisive
# for STRidge/RFE to recover the correct sparse PDE.
#
# We therefore use a fixed p=6 for both the direct and IBP weak forms
# (--test-fn-rule legacy; --test-fn-rule mb follows eq. (4.3) instead).
SPACE_P_FIXED = 6       # spatial test-function degree (direct and IBP)
SPACE_STRIDE = 4        # spatial query-point stride
# Temporal test functions (time_weak / spacetime_weak_* methods).  p_t always
# follows eq. (4.3) of Messenger & Bortz with ᾱ_t = 1.
TIME_M_GRID = [4, 8, 12]
TAU = 1e-10             # real-space decay tolerance (paper value)
TAU_HAT = 2.0           # Fourier-tail parameter (paper: 1 to 3)
ALPHA_X = 4             # max spatial derivative order in every library (Δ²)
ALPHA_T = 1             # ∂_t

# method -> (library.spacetime formulation, spatial IBP); see the module
# docstring.  Per-snapshot methods build rows as before 2026-09; the others go
# through library.spacetime.build_weak_system.
METHODS = {
    "pointwise": ("pointwise", False),
    "space_weak_direct": ("space_weak", False),
    "space_weak_ibp": ("space_weak", True),
    "time_weak": ("time_weak", False),
    "spacetime_weak_direct": ("spacetime_weak", False),
    "spacetime_weak_ibp": ("spacetime_weak", True),
}
PER_SNAPSHOT_METHODS = ("pointwise", "space_weak_direct", "space_weak_ibp")
SPACE_WEAK_METHODS = tuple(k for k, (form, _) in METHODS.items()
                           if form in ("space_weak", "spacetime_weak"))
SPACE_IBP_METHODS = tuple(k for k, (_, ibp) in METHODS.items() if ibp)
EXTRA_FIELDS = ["space_p", "time_m", "time_p", "time_ibp",
                "test_fn_rule", "centering", "space_boundary",
                "const_active", "intercept_rel"]

# (CH epsilon is stored in equations.CH_EPSILON)


def _get_meta(name):
    """Return (true_fn, skip_fraction, mask_fn) for a dataset name."""
    if name in EQUATION_META:
        return EQUATION_META[name]
    for prefix in sorted(EQUATION_META, key=len, reverse=True):
        if name.startswith(prefix):
            return EQUATION_META[prefix]
    raise KeyError(f"No metadata for equation: {name}")


# ── loading ─────────────────────────────────────────────────────────────────


def load_dataset(npz_path):
    """Load a .npz data file.

    Returns
    -------
    snapshots, times, Lx, Ly, N, dt_snap, obs_noise, eq_name,
    profile_type, profile_direction
    """
    data = np.load(npz_path)
    snapshots = data["snapshots"]
    times = data["times"]
    Lx = float(data["Lx"])
    Ly = float(data["Ly"])
    N = int(data["N"])
    dt_snap = float(data["dt_snap"])
    obs_noise = float(data["obs_noise"])

    # eq_name: prefer .npz metadata (added 2026-06-30), fall back to
    # filename parsing for older datasets.
    if "eq_name" in data:
        eq_name = str(data["eq_name"])
    else:
        eq_name = os.path.basename(npz_path).rsplit("_", 1)[0]
        for suffix in ["_clean", "_noise2pct", "_noise5pct"]:
            if eq_name.endswith(suffix):
                eq_name = eq_name[: -len(suffix)]
                break

    # Anisotropic noise metadata (absent for isotropic datasets)
    profile_type = str(data["profile_type"]) if "profile_type" in data else ""
    profile_direction = str(data["profile_direction"]) if "profile_direction" in data else ""

    return (snapshots, times, Lx, Ly, N, dt_snap, obs_noise, eq_name,
            profile_type, profile_direction)


# ── method builders ─────────────────────────────────────────────────────────


def build_pointwise(f_snaps, ft_snaps, idx_sets, Lx, Ly, library, backend):
    """Build pointwise (Theta, Ft) for a given library and backend."""
    with derivatives.use_backend(backend):
        rows_theta, rows_ft = [], []
        for s in range(len(f_snaps)):
            if library == "multigraph":
                T, desc = build_multigraph_library(f_snaps[s], Lx, Ly)
            elif library == "mod_multigraph":
                T, desc = build_mod_multigraph_library(f_snaps[s], Lx, Ly)
            elif library == "polyjet":
                T, desc = build_polyjet_library(
                    f_snaps[s], Lx, Ly, max_order=4, max_degree=3,
                )
            elif library == "wsindy":
                T, desc = build_wsindy_library(
                    f_snaps[s], Lx, Ly, max_order=4, max_poly=3,
                )
            else:
                T, desc = build_pdefind_library(
                    f_snaps[s], Lx, Ly, max_order=4, max_poly=3,
                )
            rows_theta.append(T[idx_sets[s]])
            rows_ft.append(ft_snaps[s].ravel()[idx_sets[s]])
    Theta = np.vstack(rows_theta).astype(float)
    Ft = np.concatenate(rows_ft).astype(float)
    return Theta, Ft, desc


def build_space_weak(f_snaps, ft_snaps, idx_sets, Lx, Ly, library, backend,
                     m, space_ibp, p=SPACE_P_FIXED, boundary="valid"):
    """Build per-snapshot space-weak (G, b) for a given library and backend.

    ``p`` defaults to the legacy fixed degree (see SPACE_P_FIXED)."""
    with derivatives.use_backend(backend), wf.use_weak_boundary(boundary):
        rows_G, rows_b = [], []
        for s in range(len(f_snaps)):
            f = f_snaps[s]
            ft = ft_snaps[s]
            if library == "multigraph":
                G, desc = build_multigraph_library(
                    f, Lx, Ly, max_p=2, max_n=3,
                    formulation="space_weak", m=m, p=p,
                    space_ibp=space_ibp,
                )
            elif library == "mod_multigraph":
                G, desc = build_mod_multigraph_library(
                    f, Lx, Ly, max_p=2, max_n=3,
                    formulation="space_weak", m=m, p=p,
                    space_ibp=space_ibp,
                )
            elif library == "polyjet":
                G, desc = build_polyjet_library(
                    f, Lx, Ly, max_order=4, max_degree=3,
                    formulation="space_weak", m=m, p=p, space_ibp=space_ibp,
                )
            elif library == "wsindy":
                G, desc = build_wsindy_library(
                    f, Lx, Ly, max_order=4, max_poly=3,
                    formulation="space_weak", m=m, p=p, space_ibp=space_ibp,
                )
            else:
                G, desc = build_pdefind_library(
                    f, Lx, Ly, max_order=4, max_poly=3,
                    formulation="space_weak", m=m, p=p, space_ibp=space_ibp,
                )
            b_w = wf.build_weak_rhs(ft, m=m, p=p)

            Nv_x, Nv_y = wf.weak_output_shape(f.shape, m)
            ix = np.arange(0, Nv_x, SPACE_STRIDE)
            iy = np.arange(0, Nv_y, SPACE_STRIDE)
            IX, IY = np.meshgrid(ix, iy, indexing="ij")
            flat_idx = (IX * Nv_y + IY).ravel()

            rng = np.random.default_rng(s)
            n_per = len(idx_sets[s])
            chosen = rng.choice(
                len(flat_idx), size=min(n_per, len(flat_idx)), replace=False,
            )
            q_idx = flat_idx[chosen]

            rows_G.append(G[q_idx])
            rows_b.append(b_w[q_idx])

    G_all = np.vstack(rows_G).astype(float)
    b_all = np.concatenate(rows_b).astype(float)
    return G_all, b_all, desc


def library_builder(library, Lx, Ly):
    """Single-snapshot builder ``(f, formulation, **weak_kw) -> (Theta, desc)``
    with the same library settings as :func:`build_pointwise` /
    :func:`build_space_weak`."""
    def build(f, formulation, **kw):
        if library == "multigraph":
            return build_multigraph_library(f, Lx, Ly, max_p=2, max_n=3,
                                            formulation=formulation, **kw)
        if library == "mod_multigraph":
            return build_mod_multigraph_library(f, Lx, Ly, max_p=2, max_n=3,
                                                formulation=formulation, **kw)
        if library == "polyjet":
            return build_polyjet_library(f, Lx, Ly, max_order=4, max_degree=3,
                                         formulation=formulation, **kw)
        if library == "wsindy":
            return build_wsindy_library(f, Lx, Ly, max_order=4, max_poly=3,
                                        formulation=formulation, **kw)
        return build_pdefind_library(f, Lx, Ly, max_order=4, max_poly=3,
                                     formulation=formulation, **kw)
    return build


def build_weak_in_time(snapshots, times, Lx, Ly, library, backend, method,
                       cfg, space_stride, time_stride, mask_fn=None):
    """Build (G, b) for a time_weak / spacetime_weak_* method from the full
    (transient-removed) trajectory; see ``library/spacetime.py``."""
    formulation, space_ibp = METHODS[method]
    with derivatives.use_backend(backend):
        G, b, desc, _ = build_weak_system(
            snapshots, times, library_builder(library, Lx, Ly), formulation,
            m=cfg.get("m"), p=cfg.get("p"), space_ibp=space_ibp,
            m_t=cfg["m_t"], p_t=cfg["p_t"], time_ibp=cfg["time_ibp"],
            spatial_stride=space_stride, time_stride=time_stride,
            point_mask_fn=mask_fn, boundary=cfg.get("space_boundary"),
        )
    return G.astype(float), b.astype(float), desc


def test_fn_configs(method, rule, space_m_vals, usable_space_ibp_m_vals,
                    time_m_vals, time_ibp_modes, n_frames, tau=TAU,
                    boundary_modes=("valid",)):
    """Test-function configurations to build for *method*.

    Returns a list of dicts with any of ``m, p, space_boundary`` (spatial),
    ``m_t, p_t, time_ibp`` (temporal); ``[{}]`` for pointwise.  ``rule`` is the output of
    ``weak_form.mb_test_function_params`` (``--test-fn-rule mb``) or None
    (legacy grids).  Temporal supports that do not fit in the available
    frames (``n_frames``, two fewer when f_t is finite-differenced) are
    dropped with a notice.
    """
    if method == "pointwise":
        return [{}]
    if rule is not None:
        spatial = [dict(m=rule["m"], p=rule["p"])]
        temporal = [(rule["m_t"], rule["p_t"])]
    else:
        ms = (usable_space_ibp_m_vals if method in SPACE_IBP_METHODS
              else space_m_vals)
        spatial = [dict(m=m, p=SPACE_P_FIXED) for m in ms]
        temporal = [(m_t, wf.select_test_function_degree(m_t, ALPHA_T, tau))
                    for m_t in time_m_vals]
    spatial = [dict(sp, space_boundary=bd)
               for sp in spatial for bd in boundary_modes]
    if method in PER_SNAPSHOT_METHODS:
        return spatial
    if method == "time_weak":
        spatial = [{}]

    cfgs = []
    for sp in spatial:
        for m_t, p_t in temporal:
            for time_ibp in time_ibp_modes:
                frames = n_frames if time_ibp else n_frames - 2
                if 2 * m_t + 1 > frames:
                    print(f"  [{method}] skipping m_t={m_t} "
                          f"(time_ibp={time_ibp}): needs {2 * m_t + 1} "
                          f"frames, {frames} available", flush=True)
                    continue
                cfgs.append(dict(sp, m_t=m_t, p_t=p_t, time_ibp=time_ibp))
    return cfgs


def cfg_columns(cfg, rule_name):
    """CSV columns describing a test-function configuration (blank when the
    method has no such parameter)."""
    return dict(
        space_p=cfg.get("p", ""),
        time_m=cfg.get("m_t", ""),
        time_p=cfg.get("p_t", ""),
        time_ibp=("" if "time_ibp" not in cfg
                  else ("on" if cfg["time_ibp"] else "off")),
        test_fn_rule=rule_name if cfg else "",
        centering="",
        space_boundary=cfg.get("space_boundary", ""),
    )


def const_column_reporter(Theta, Ft, desc):
    """Return ``w -> {const_active, intercept_rel}`` for this design matrix.

    ``regression.centering`` writes the recovered intercept back into the
    constant column's slot, so an intercept that survives the threshold is
    counted in ``fp`` exactly like any other false positive -- yet it is the
    one term centering exists to remove, and a centering comparison that
    cannot separate the two says little.  These columns do separate it:
    whether the constant slot is active, and the intercept's size relative to
    mean|f_t|, which is the quantity ``--centering-tol`` thresholds.

    The weak-form constant column is ⟨ψ, 1⟩ rather than all ones, so the
    intercept is ``w[j] * mean(Theta[:, j])``, matching ``centering``'s own
    back-substitution.  Both columns are blank when the library has no
    constant column.
    """
    j = _resolve_const_index(Theta, desc, None, "1")
    if j is None:
        blank = dict(const_active="", intercept_rel="")
        return lambda w: blank
    col_mean = float(Theta[:, j].mean())
    scale = float(np.mean(np.abs(Ft))) or 1.0

    def report(w):
        coef = float(w[j])
        return dict(const_active=int(abs(coef) > 1e-10),
                    intercept_rel=f"{abs(coef * col_mean) / scale:.6e}")
    return report


# ── sweep one ───────────────────────────────────────────────────────────────


def fit(Theta, Ft, solve, center=False, center_tol=0.0, desc=None, cache=None):
    """``solve(Theta, Ft)`` directly, or through ``regression.centering``
    (constant column dropped, columns and target centred, intercept restored
    into the '1' slot).  ``cache`` is a ``regression.SolverCache`` for the
    centred system."""
    if not center:
        return solve(Theta, Ft)
    return centering(Theta, Ft, solve, tol=center_tol, descs=desc, cache=cache)


def sweep_one(Theta, Ft, true_eq, lam, d_tol, algorithm, mstls_alpha=0.5,
              mstls_M=None, center=False, center_tol=0.0, desc=None,
              cache=None):
    """Run regression and return (n_active, w, lam_out).

    ``cache`` is a ``regression.SolverCache`` shared by every call on this
    (Theta, Ft); it only avoids recomputation, results are unchanged."""
    lam_out = lam
    if algorithm == "mstls":
        # MSTLS does its own internal λ sweep; lam/d_tol are ignored.  Under
        # centering the solver sees a reduced, centred matrix, so the
        # dominant-balance scales are recomputed on it.
        def solve(X, y):
            nonlocal lam_out
            M = mstls_M
            if center and mstls_M is not None:
                M = wf.compute_m_scale(X, y)
            w_, lam_out, _ = TrainMSTLS(
                X, y, lambdas=None, M=M, n_lambdas=MSTLS_N_LAMBDAS,
                maxit=MSTLS_MAXITS, alpha=mstls_alpha, cache=cache,
            )
            return w_
    elif algorithm == "rfe":
        # lam is the ridge regularisation for RFE (0 = pure OLS)
        def solve(X, y):
            return TrainRFE(X, y, lam=lam, cache=cache)
    else:
        def solve(X, y):
            # Sequential Threshold Ridge
            return TrainSTRidge(X, y, lam, d_tol, cache=cache)
    w = fit(Theta, Ft, solve, center, center_tol, desc, cache=cache)
    n_active = int(np.count_nonzero(np.abs(w) > 1e-10))
    return n_active, w, lam_out


def cond_number(Theta):
    """Condition number (2-norm) of Theta."""
    return float(np.linalg.cond(Theta))


# ── main ────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Canonical ablation sweep for equivariant PDE discovery."
    )
    parser.add_argument(
        "--data-dir", default=DATA_DIR,
        help="Directory containing .npz data files.",
    )
    parser.add_argument(
        "--equations", nargs="+", default=None,
        help="Which equation datasets to sweep (default: all .npz in data dir).",
    )
    parser.add_argument(
        "--noise", nargs="+", type=float, default=None,
        help="Noise levels to sweep (default: all found).",
    )
    parser.add_argument(
        "--methods", nargs="+",
        default=list(PER_SNAPSHOT_METHODS),
        choices=list(METHODS) + list(METHOD_ALIASES),
        help="Formulations to build (see module docstring).  The pre-2026-09 "
             "names weak_direct / weak_ibp are accepted for "
             "space_weak_direct / space_weak_ibp.",
    )
    parser.add_argument(
        "--algorithms", nargs="+",
        default=["stridge", "mstls", "rfe"],
        choices=["stridge", "mstls", "rfe", "brute", "oracle"],
    )
    parser.add_argument(
        "--libraries", nargs="+",
        default=["multigraph", "mod_multigraph", "pdefind", "polyjet", "wsindy"],
        choices=["multigraph", "mod_multigraph", "pdefind", "polyjet", "wsindy"],
    )
    parser.add_argument(
        "--representable-only", action="store_true",
        help="Per dataset, sweep only the --libraries that can represent the "
             "equation (equations.REPRESENTABLE_LIBRARIES): PDE-FIND is "
             "dropped for CH, PME and KS, WSINDy for KS.  Off by default so "
             "that earlier invocations reproduce unchanged.",
    )
    parser.add_argument(
        "--backends", nargs="+",
        default=BACKENDS,
        choices=BACKENDS,
    )
    parser.add_argument(
        "--space-m", "--weak-m", dest="space_m", nargs="+", type=int,
        default=SPACE_M_GRID,
        help="Spatial test-function half-supports m for methods that are weak "
             "in space (--test-fn-rule legacy).  --weak-m is the old name.",
    )
    parser.add_argument(
        "--time-m", nargs="+", type=int, default=TIME_M_GRID,
        help="Temporal half-supports m_t for time_weak / spacetime_weak_* "
             "methods (--test-fn-rule legacy); p_t follows the M&B degree "
             "rule.",
    )
    parser.add_argument(
        "--space-boundary", nargs="+", choices=["valid", "periodic"],
        default=["valid"],
        help="Spatial boundary of the weak convolution for methods that are "
             "weak in space: valid (support inside the grid) or periodic "
             "(support wraps; every grid point is a query centre).  Pass both "
             "to compare.  periodic is the setting for current experiments "
             "(periodic domains); valid is the default for reproducibility.",
    )
    parser.add_argument(
        "--test-fn-rule", choices=["legacy", "mb"], default="legacy",
        help="legacy: m from --weak-m with fixed p, m_t from --weak-mt.  "
             "mb: (m, p, m_t, p_t) per dataset by the Messenger & Bortz rule, "
             "for every weak method.",
    )
    parser.add_argument(
        "--tau", type=float, default=TAU,
        help="M&B real-space decay tolerance τ (degree rule, eq. 4.3).",
    )
    parser.add_argument(
        "--tau-hat", type=float, default=TAU_HAT,
        help="M&B Fourier-tail parameter τ̂ (Appendix A; paper uses 1-3). "
             "--test-fn-rule mb only.",
    )
    parser.add_argument(
        "--time-ibp", nargs="+", choices=["on", "off"], default=["on"],
        help="time_weak / spacetime_weak_*: b = φ_t′∗f (on) or φ_t∗f_t with "
             "FD f_t (off).  Pass both to compare.",
    )
    parser.add_argument(
        "--centering", nargs="+", choices=["off", "on"], default=["off"],
        help="Fit through regression.centering (intercept partialled out). "
             "Pass both to compare.",
    )
    parser.add_argument(
        "--centering-tol", type=float, default=0.05,
        help="Relative intercept threshold for centering: w0 is zeroed when "
             "|w0| < tol·mean|b| (0 disables).  0.05 is the value recommended "
             "in regression.centering; with 0 the unthresholded spurious "
             "intercept counts as a false positive in every centred row.",
    )
    parser.add_argument(
        "--space-stride", type=int, default=SPACE_STRIDE,
        help="Spatial query-point stride for time_weak / spacetime_weak_* "
             "methods (the per-snapshot methods always use SPACE_STRIDE).",
    )
    parser.add_argument(
        "--time-stride", type=int, default=1,
        help="Temporal query-point stride for time_weak / spacetime_weak_* "
             "methods.",
    )
    parser.add_argument(
        "--post-hoc-noise", nargs="+", type=float, default=None,
        help="Apply additive Gaussian noise (as fraction of std) to loaded "
             "snapshots. Useful for noise-robustness ablation.",
    )
    parser.add_argument(
        "--spectral-cutoff", type=float, default=None,
        help="Apply spectral low-pass filter to snapshots before building "
             "libraries (cutoff wavenumber).",
    )
    parser.add_argument(
        "--mstls-alpha", nargs="+", type=float, default=None,
        help="MSTLS alpha values (projection-vs-sparsity trade-off). "
             "Default: 0.0 0.25 0.5 0.75 1.0.",
    )
    parser.add_argument(
        "--rfe-lam", nargs="+", type=float, default=None,
        help="Ridge regularisation for RFE (0 = pure OLS). "
             "Default: 0.0.",
    )
    parser.add_argument(
        "--brute-l0", nargs="+", type=float, default=None,
        help="L₀ sparsity penalties for brute-force exhaustive search. "
             f"Default: {BRUTE_L0_PENALTIES}.",
    )
    parser.add_argument(
        "--brute-max-k", type=int, default=None,
        help=f"Max subset size for brute-force (default: None for ≤20 cols, "
             f"{BRUTE_MAX_K} for >20).",
    )
    parser.add_argument(
        "--name", default="",
        help="Optional tag appended to output filename.",
    )
    parser.add_argument(
        "--no-oracle", action="store_true",
        help="Skip oracle (ground-truth residual) computation.",
    )
    parser.add_argument(
        "--aniso", action="store_true",
        help="Include anisotropic-noise datasets (excluded by default).",
    )
    parser.add_argument(
        "--eval-mode", choices=["dense", "qr"], default="dense",
        help="dense: score every row against the full clean evaluation "
             "matrix (legacy).  qr: replace it by the R factor of its QR "
             "decomposition (evaluation.FactoredEval); E_coeff agrees to "
             "round-off, fp/fn/exact_recovery are unchanged, and scoring no "
             "longer scales with the grid size.",
    )
    parser.add_argument(
        "--out", default=None,
        help="Exact output CSV path. When set, overrides the auto-generated "
             "timestamp-based filename. Use when running multiple workers "
             "that each need their own output file.",
    )
    args = parser.parse_args()

    # ── discover datasets ────────────────────────────────────────────────
    data_dir = args.data_dir
    all_npz = sorted([
        f for f in os.listdir(data_dir)
        if f.endswith(".npz") and not f.startswith(".")
        and (args.aniso or "aniso" not in f)  # skip unless --aniso
        # This script is 2-D throughout (it passes no Lz, and evaluates in a
        # 2-D jet basis).  The n-D datasets written by generate_nd_data.py are
        # named <equation>_nd<ndim>_<label>.npz and share the 2-D equation
        # prefix, so an unfiltered `--equations ac` would pull in
        # allen_cahn_nd3_*.npz and fail in axis_lengths().  Use sweep_nd.py
        # for those.
        and not re.search(r"_nd\d+_", f)
    ])

    if args.equations:
        eq_set = set(args.equations)
        # Map short names → filename prefixes (both short and legacy long
        # names; generate_data.py uses short names, old data used long names)
        _eq_to_filename = {
            "ac": {"allen_cahn", "ac"},
            "sh": {"swift_hohenberg", "sh"},
            "ch": {"cahn_hilliard", "ch"},
            "pme": {"porous_medium", "pme"},
            "ks_cell": {"ks_cellular", "ks_cell"},
            "ks_wchaos": {"ks_weak_chaotic", "ks_wchaos"},
            "ks_chaos": {"ks_chaotic", "ks_chaos"},
        }
        eq_prefixes = set()
        for eq in eq_set:
            eq_prefixes |= _eq_to_filename.get(eq, {eq})
        all_npz = [f for f in all_npz
                   if any(f.startswith(eq + "_") for eq in eq_prefixes)]
    if args.noise is not None:
        noise_labels = {0.0: "clean", 0.02: "noise2pct", 0.05: "noise5pct"}
        allowed = {noise_labels[n] for n in args.noise}
        all_npz = [f for f in all_npz
                   if any(f.endswith("_" + nl + ".npz") for nl in allowed)]

    if not all_npz:
        print(f"No .npz files found in {data_dir} matching filters.")
        return

    methods = list(dict.fromkeys(canonical_method(m) for m in args.methods))
    algorithms = args.algorithms
    libraries = args.libraries
    backends = args.backends
    space_m_vals = args.space_m

    # At small m the grid cannot resolve a degree-p test function at all.
    # Drop those m for spatial IBP; the direct weak form only needs φ itself
    # (no derivatives) and tolerates coarser resolution.
    usable_space_ibp_m_vals = [
        m for m in space_m_vals
        if wf.test_function_is_resolved(m, SPACE_P_FIXED, max_deriv=4)
    ]
    dropped = [m for m in space_m_vals if m not in usable_space_ibp_m_vals]
    if (dropped and args.test_fn_rule == "legacy"
            and (set(SPACE_IBP_METHODS) & set(methods))):
        print(
            f"  [space_weak_ibp] skipping m={dropped}: 2m+1 grid points "
            f"cannot resolve the degree-{SPACE_P_FIXED} test function "
            f"for integration by parts "
        )

    post_hoc_noises = args.post_hoc_noise or [None]
    spectral_cutoffs = [args.spectral_cutoff] if args.spectral_cutoff else [None]
    mstls_alphas = args.mstls_alpha or MSTLS_ALPHAS
    rfe_lams = args.rfe_lam or RFE_LAMS
    brute_l0_penalties = args.brute_l0 or BRUTE_L0_PENALTIES
    oracle_only = (algorithms == ["oracle"])
    compute_oracle = not args.no_oracle

    # ── output ───────────────────────────────────────────────────────────
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    if args.out:
        out_path = args.out
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    else:
        tag = f"_{args.name}" if args.name else ""
        out_path = os.path.join(RESULTS_DIR, f"sweep_{timestamp}{tag}.csv")
    print(f"Output → {out_path}", flush=True)

    fieldnames = [
        "equation", "noise", "method", "algorithm", "library", "backend",
        "space_m", "lam", "d_tol", "mstls_alpha",
        "n_active", "fp", "fn", "exact_recovery", "E_coeff", "max_coeff_err",
        "E_data", "cond_theta", "n_samples",
        "post_hoc_noise", "spectral_cutoff",
        "profile_type", "profile_direction",
        "time_ms",
    ] + EXTRA_FIELDS
    time_ibp_modes = [v == "on" for v in args.time_ibp]
    center_modes = [v == "on" for v in args.centering]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        f.flush()

        t_start_total = time.time()
        total_combos = 0
        eq_summaries = []  # per-equation summary rows

        for npz_file in all_npz:
            npz_path = os.path.join(data_dir, npz_file)
            print(f"\n{'='*70}")
            print(f"Dataset: {npz_file}")
            print(f"{'='*70}")

            # Load
            snapshots_orig, times, Lx, Ly, N, dt_snap, obs_noise, eq_name, \
                profile_type, profile_direction = load_dataset(npz_path)
            eq_name = normalize_eq_name(eq_name)
            true_fn, skip_frac, mask_fn = _get_meta(eq_name)

            # Evaluation must always run on clean data (pointwise spectral).
            # For noisy datasets, load the clean counterpart .npz file.
            # For clean datasets, the loaded data is the clean reference.
            if obs_noise > 0:
                # For anisotropic datasets (e.g. ac_aniso_radial_decay_noise2pct),
                # the clean counterpart is {name}_clean.npz, not
                # {name}_aniso_*_clean.npz.
                if "_aniso_" in npz_path:
                    clean_path = os.path.join(
                        os.path.dirname(npz_path),
                        f"{eq_name}_clean.npz",
                    )
                else:
                    clean_path = npz_path.replace(
                        "_noise2pct", "_clean",
                    ).replace("_noise5pct", "_clean")
                if not os.path.exists(clean_path):
                    raise FileNotFoundError(
                        f"Clean evaluation data required but not found: "
                        f"{clean_path}.  Generate it with "
                        f"`python scripts/generate_data.py --equations "
                        f"{eq_name}`."
                    )
                _clean_data = np.load(clean_path)
                snapshots_clean_eval = _clean_data["snapshots"]
                print(f"    clean eval: {os.path.basename(clean_path)}")
            else:
                snapshots_clean_eval = snapshots_orig  # data IS clean

            # CH epsilon is fixed in equations.CH_EPSILON
            is_ch = eq_name == "ch"
            eps_list = [CH_EPSILON] if is_ch else [None]

            eq_exact_count = 0
            eq_total = 0

            for post_hoc_noise in post_hoc_noises:
                # Apply post-hoc noise
                if post_hoc_noise is not None and post_hoc_noise > 0:
                    rng_noise = np.random.default_rng(42)
                    snapshots = snapshots_orig.copy()
                    for s in range(len(snapshots)):
                        std = np.std(snapshots[s])
                        snapshots[s] += post_hoc_noise * std * \
                            rng_noise.standard_normal(snapshots[s].shape)
                    noise_label = obs_noise + post_hoc_noise
                else:
                    snapshots = snapshots_orig
                    noise_label = obs_noise

                # Apply spectral cutoff
                if (args.spectral_cutoff is not None
                        and args.spectral_cutoff > 0):
                    snapshots_proc = snapshots.copy()
                    for s in range(len(snapshots_proc)):
                        snapshots_proc[s] = spectral_lowpass(
                            snapshots_proc[s], Lx, Ly, args.spectral_cutoff,
                        )
                else:
                    snapshots_proc = snapshots

                for ch_eps in eps_list:
                    print(f"  eq={eq_name}  noise={noise_label:.3f}"
                          + (f"  cutoff={args.spectral_cutoff}"
                             if args.spectral_cutoff else ""),
                          flush=True)

                    # Ground truth and the clean evaluation system, shared
                    # with scripts/protocol.py through eval_setup so that a
                    # second driver can score without copying this block.
                    ev = prepare_evaluation(
                        snapshots_proc, snapshots_clean_eval, times, Lx, Ly,
                        eq_name, true_fn, skip_frac, mask_fn=mask_fn,
                        ch_eps=ch_eps, eval_mode=args.eval_mode,
                    )
                    desc_eq_ref, true_eq = ev.desc_eq_ref, ev.true_eq
                    desc_jet, true_jet = ev.desc_jet, ev.true_jet
                    M_ws_to_jet = ev.M_ws_to_jet
                    Theta_eval, Ft_eval = ev.Theta_eval, ev.Ft_eval
                    f_snaps, ft_snaps = ev.f_snaps, ev.ft_snaps
                    idx_sets, n_snaps = ev.idx_sets, ev.n_snaps

                    # ── Test-function configurations ─────────────────────
                    # time_weak / spacetime_weak_* use the whole trajectory (after
                    # the same transient skip as make_regression_data).
                    n_skip = max(0, int(skip_frac * len(snapshots_proc)))
                    snaps_st = snapshots_proc[n_skip:]
                    times_st = times[n_skip:]
                    rule = None
                    if args.test_fn_rule == "mb":
                        rule = wf.mb_test_function_params(
                            snaps_st, spatial_axes=range(1, snaps_st.ndim),
                            time_axis=0, alpha_x=ALPHA_X, alpha_t=ALPHA_T,
                            tau=args.tau, tau_hat=args.tau_hat,
                        )
                        print(f"    M&B rule (tau_hat={args.tau_hat}): "
                              f"k_x={rule['k_x']} m={rule['m']} p={rule['p']}"
                              f"  k_t={rule['k_t']} m_t={rule['m_t']} "
                              f"p_t={rule['p_t']}", flush=True)
                        if not wf.test_function_is_resolved(
                                rule["m"], rule["p"], max_deriv=ALPHA_X):
                            print(f"    [warn] m={rule['m']}, p={rule['p']} "
                                  f"fails test_function_is_resolved; IBP "
                                  f"columns may be inaccurate", flush=True)
                    method_cfgs = {
                        method: test_fn_configs(
                            method, rule, space_m_vals,
                            usable_space_ibp_m_vals, args.time_m,
                            time_ibp_modes, len(snaps_st), tau=args.tau,
                            boundary_modes=args.space_boundary,
                        )
                        for method in methods
                    }

                    # ── Pre-build method matrices ────────────────────────
                    prebuilt = {}

                    # Oracle-only: only need multigraph matrices
                    _libs = ["multigraph"] if oracle_only else libraries
                    if args.representable_only and not oracle_only:
                        _repr = REPRESENTABLE_LIBRARIES.get(eq_name)
                        if _repr is not None:
                            skipped = [lb for lb in _libs if lb not in _repr]
                            _libs = [lb for lb in _libs if lb in _repr]
                            if skipped:
                                print(f"    not representable, skipped: "
                                      f"{' '.join(skipped)}", flush=True)
                    for method, library, backend in itertools.product(
                        methods, _libs, backends
                    ):
                        for cfg in method_cfgs[method]:
                            m_val = cfg.get("m")
                            actual_key = (method, library, backend,
                                          m_val if m_val is not None else 0,
                                          cfg.get("p", 0), cfg.get("m_t", 0),
                                          cfg.get("p_t", 0),
                                          cfg.get("time_ibp"),
                                          cfg.get("space_boundary"))

                            if actual_key in prebuilt:
                                continue

                            t0 = time.time()
                            print(f"    Building: method={method} "
                                  f"lib={library} backend={backend}"
                                  + (f" m={m_val}" if m_val else "")
                                  + (f" p={cfg['p']}"
                                     if "p" in cfg and "m_t" in cfg else "")
                                  + (f" m_t={cfg['m_t']} p_t={cfg['p_t']} "
                                     f"time_ibp={cfg['time_ibp']}"
                                     if "m_t" in cfg else "")
                                  + (f" boundary={cfg['space_boundary']}"
                                     if cfg.get("space_boundary") == "periodic"
                                     else ""),
                                  flush=True)

                            if method == "pointwise":
                                Theta, Ft, desc = build_pointwise(
                                    f_snaps, ft_snaps, idx_sets, Lx, Ly,
                                    library, backend,
                                )
                            elif method in PER_SNAPSHOT_METHODS:
                                Theta, Ft, desc = build_space_weak(
                                    f_snaps, ft_snaps, idx_sets, Lx, Ly,
                                    library, backend, m_val,
                                    METHODS[method][1], p=cfg["p"],
                                    boundary=cfg["space_boundary"],
                                )
                            else:
                                Theta, Ft, desc = build_weak_in_time(
                                    snaps_st, times_st, Lx, Ly, library,
                                    backend, method, cfg,
                                    args.space_stride,
                                    args.time_stride, mask_fn=mask_fn,
                                )

                            cond_val = cond_number(Theta)
                            prebuilt[actual_key] = (Theta, Ft, desc, cond_val,
                                                    cfg)
                            elapsed = time.time() - t0
                            print(f"      shape={Theta.shape}  "
                                  f"cond={cond_val:.2e}  "
                                  f"({elapsed:.1f}s)", flush=True)

                    # ── Oracle: per-configuration ground truth ──────────
                    if compute_oracle:
                        for actual_key, (Theta, Ft, desc, cond_val, cfg) in \
                                prebuilt.items():
                            method, library, backend, m_val = actual_key[:4]
                            if library != "multigraph":
                                continue
                            # Use method-specific (Theta, Ft) directly
                            Ft_pred = Theta @ true_eq
                            E_data = float(
                                np.linalg.norm(Ft_pred - Ft)
                                / max(np.linalg.norm(Ft), 1e-14)
                            )
                            row = dict(
                                equation=eq_name,
                                noise=noise_label,
                                method=method,
                                algorithm="oracle",
                                library=library,
                                backend=backend,
                                space_m=(m_val if m_val != 0 else ""),
                                lam="", d_tol="", mstls_alpha="",
                                n_active=int(np.sum(np.abs(true_eq) > 1e-10)),
                                fp=0, fn=0, exact_recovery=1,
                                E_coeff="0.000000",
                                max_coeff_err="0.000000",
                                E_data=f"{E_data:.6f}",
                                cond_theta=cond_val,
                                n_samples=len(Ft),
                                post_hoc_noise=(
                                    post_hoc_noise if post_hoc_noise else ""),
                                spectral_cutoff=(
                                    args.spectral_cutoff
                                    if args.spectral_cutoff else ""),
                                profile_type=profile_type,
                                profile_direction=profile_direction,
                                time_ms="0.0",
                            )
                            row.update(cfg_columns(cfg, args.test_fn_rule))
                            writer.writerow(row)
                            total_combos += 1

                    if oracle_only:
                        print(f"  [{eq_name}] oracle done", flush=True)
                        continue

                    # ── Sweep hyperparams ──────────────────────────────
                    # prebuilt preserves (method, library, backend, config)
                    # insertion order, so rows come out in the same order as
                    # before the space-time extension.
                    for actual_key, (Theta, Ft, desc, cond_val, cfg) in \
                            prebuilt.items():
                        method, library, backend, m_val = actual_key[:4]
                        cfg_cols = cfg_columns(cfg, args.test_fn_rule)
                        # Data-only solver work, shared by every setting and
                        # both centering modes on this matrix.
                        solver_cache = SolverCache()
                        const_cols = const_column_reporter(Theta, Ft, desc)

                        for center in center_modes:
                            cfg_cols["centering"] = "on" if center else "off"

                            for algorithm in algorithms:
                                # MSTLS: weak-form only (skip pointwise)
                                if algorithm == "mstls":
                                    if method == "pointwise":
                                        continue

                                    # M-scaling only for WSINDy; skip for
                                    # multigraph/standard where column
                                    # norms span orders of magnitude.
                                    if library in ("wsindy", "wsindy_proj"):
                                        mstls_M = wf.compute_m_scale(Theta, Ft)
                                    else:
                                        mstls_M = None

                                    for mstls_alpha in mstls_alphas:
                                        t_row = time.time()
                                        n_active, w, lam_out = sweep_one(
                                            Theta, Ft, true_eq, 0.0, 0.0,
                                            algorithm,
                                            mstls_alpha=mstls_alpha,
                                            mstls_M=mstls_M,
                                            center=center,
                                            center_tol=args.centering_tol,
                                            desc=desc, cache=solver_cache,
                                        )

                                        eval_result = evaluate_unified(
                                            w, library, desc, desc_jet,
                                            Theta_eval, Ft_eval, true_jet,
                                            M_ws_to_jet=M_ws_to_jet,
                                        )
                                        dt_ms = (time.time() - t_row) * 1000

                                        row = dict(
                                            equation=eq_name,
                                            noise=noise_label,
                                            method=method,
                                            algorithm=algorithm,
                                            library=library,
                                            backend=backend,
                                            space_m=(m_val if m_val != 0
                                                      else ""),
                                            lam=lam_out,
                                            d_tol="",
                                            mstls_alpha=mstls_alpha,
                                            n_active=n_active,
                                            fp=eval_result["fp"],
                                            fn=eval_result["fn"],
                                            exact_recovery=eval_result[
                                                "exact_recovery"],
                                            E_coeff=eval_result["E_coeff"],
                                            max_coeff_err=eval_result[
                                                "max_coeff_err"],
                                            E_data=f"{np.linalg.norm(Theta @ w - Ft) / max(np.linalg.norm(Ft), 1e-14):.6f}",
                                            cond_theta=cond_val,
                                            n_samples=eval_result["n_samples"],
                                            post_hoc_noise=(
                                                post_hoc_noise
                                                if post_hoc_noise else ""),
                                            spectral_cutoff=(
                                                args.spectral_cutoff
                                                if args.spectral_cutoff
                                                else ""),
                                            profile_type=profile_type,
                                            profile_direction=profile_direction,
                                            time_ms=f"{dt_ms:.1f}",
                                        )
                                        row.update(cfg_cols)
                                        row.update(const_cols(w))
                                        writer.writerow(row)
                                        total_combos += 1
                                        if eval_result["exact_recovery"]:
                                            eq_exact_count += 1
                                        eq_total += 1
                                    continue

                                # RFE: sweep over ridge regularisation
                                if algorithm == "rfe":
                                    for rfe_lam in rfe_lams:
                                        t_row = time.time()
                                        n_active, w, lam_out = sweep_one(
                                            Theta, Ft, true_eq, rfe_lam, 0.0,
                                            algorithm, center=center,
                                            center_tol=args.centering_tol,
                                            desc=desc, cache=solver_cache,
                                        )

                                        eval_result = evaluate_unified(
                                            w, library, desc, desc_jet,
                                            Theta_eval, Ft_eval, true_jet,
                                            M_ws_to_jet=M_ws_to_jet,
                                        )
                                        dt_ms = (time.time() - t_row) * 1000

                                        row = dict(
                                            equation=eq_name,
                                            noise=noise_label,
                                            method=method,
                                            algorithm=algorithm,
                                            library=library,
                                            backend=backend,
                                            space_m=(m_val if m_val != 0
                                                      else ""),
                                            lam=lam_out,
                                            d_tol="",
                                            mstls_alpha="",
                                            n_active=n_active,
                                            fp=eval_result["fp"],
                                            fn=eval_result["fn"],
                                            exact_recovery=eval_result[
                                                "exact_recovery"],
                                            E_coeff=eval_result["E_coeff"],
                                            max_coeff_err=eval_result[
                                                "max_coeff_err"],
                                            E_data=f"{np.linalg.norm(Theta @ w - Ft) / max(np.linalg.norm(Ft), 1e-14):.6f}",
                                            cond_theta=cond_val,
                                            n_samples=eval_result[
                                                "n_samples"],
                                            post_hoc_noise=(
                                                post_hoc_noise
                                                if post_hoc_noise else ""),
                                            spectral_cutoff=(
                                                args.spectral_cutoff
                                                if args.spectral_cutoff
                                                else ""),
                                            profile_type=profile_type,
                                            profile_direction=
                                                profile_direction,
                                            time_ms=f"{dt_ms:.1f}",
                                        )
                                        row.update(cfg_cols)
                                        row.update(const_cols(w))
                                        writer.writerow(row)
                                        total_combos += 1
                                        eq_total += 1
                                        if eval_result["exact_recovery"]:
                                            eq_exact_count += 1
                                    continue

                                # Brute-force: multigraphs only (same as oracle)
                                if algorithm == "brute":
                                    if library != "multigraph":
                                        continue
                                    max_k = args.brute_max_k
                                    if max_k is None:
                                        n_cols = Theta.shape[1]
                                        max_k = None if n_cols <= 20 else BRUTE_MAX_K

                                    for l0_pen in brute_l0_penalties:
                                        t_row = time.time()
                                        w = fit(
                                            Theta, Ft,
                                            lambda X, y: TrainExhaustive(
                                                X, y, max_k=max_k,
                                                l0_penalty=l0_pen,
                                            ),
                                            center, args.centering_tol, desc,
                                            cache=solver_cache,
                                        )
                                        n_active = int(np.count_nonzero(
                                            np.abs(w) > 1e-10))

                                        eval_result = evaluate_unified(
                                            w, library, desc, desc_jet,
                                            Theta_eval, Ft_eval, true_jet,
                                            M_ws_to_jet=M_ws_to_jet,
                                        )
                                        dt_ms = (time.time() - t_row) * 1000

                                        row = dict(
                                            equation=eq_name,
                                            noise=noise_label,
                                            method=method,
                                            algorithm=algorithm,
                                            library=library,
                                            backend=backend,
                                            space_m=(m_val if m_val != 0
                                                      else ""),
                                            lam="",
                                            d_tol=l0_pen,
                                            mstls_alpha="",
                                            n_active=n_active,
                                            fp=eval_result["fp"],
                                            fn=eval_result["fn"],
                                            exact_recovery=eval_result[
                                                "exact_recovery"],
                                            E_coeff=eval_result["E_coeff"],
                                            max_coeff_err=eval_result[
                                                "max_coeff_err"],
                                            E_data=f"{np.linalg.norm(Theta @ w - Ft) / max(np.linalg.norm(Ft), 1e-14):.6f}",
                                            cond_theta=cond_val,
                                            n_samples=eval_result[
                                                "n_samples"],
                                            post_hoc_noise=(
                                                post_hoc_noise
                                                if post_hoc_noise else ""),
                                            spectral_cutoff=(
                                                args.spectral_cutoff
                                                if args.spectral_cutoff
                                                else ""),
                                            profile_type=profile_type,
                                            profile_direction=
                                                profile_direction,
                                            time_ms=f"{dt_ms:.1f}",
                                        )
                                        row.update(cfg_cols)
                                        row.update(const_cols(w))
                                        writer.writerow(row)
                                        total_combos += 1
                                        eq_total += 1
                                        if eval_result["exact_recovery"]:
                                            eq_exact_count += 1
                                    continue

                                # STRidge: sweep over lam × d_tol
                                for lam in LAMS:
                                    for d_tol in D_TOLS:
                                        t_row = time.time()
                                        n_active, w, lam_out = sweep_one(
                                            Theta, Ft, true_eq, lam, d_tol,
                                            algorithm, center=center,
                                            center_tol=args.centering_tol,
                                            desc=desc, cache=solver_cache,
                                        )

                                        eval_result = evaluate_unified(
                                            w, library, desc, desc_jet,
                                            Theta_eval, Ft_eval, true_jet,
                                            M_ws_to_jet=M_ws_to_jet,
                                        )
                                        dt_ms = (time.time() - t_row) * 1000

                                        row = dict(
                                            equation=eq_name,
                                            noise=noise_label,
                                            method=method,
                                            algorithm=algorithm,
                                            library=library,
                                            backend=backend,
                                            space_m=(m_val if m_val != 0
                                                      else ""),
                                            lam=lam_out,
                                            d_tol=d_tol,
                                            mstls_alpha="",
                                            n_active=n_active,
                                            fp=eval_result["fp"],
                                            fn=eval_result["fn"],
                                            exact_recovery=eval_result[
                                                "exact_recovery"],
                                            E_coeff=eval_result["E_coeff"],
                                            max_coeff_err=eval_result[
                                                "max_coeff_err"],
                                            E_data=f"{np.linalg.norm(Theta @ w - Ft) / max(np.linalg.norm(Ft), 1e-14):.6f}",
                                            cond_theta=cond_val,
                                            n_samples=eval_result["n_samples"],
                                            post_hoc_noise=(
                                                post_hoc_noise
                                                if post_hoc_noise else ""),
                                            spectral_cutoff=(
                                                args.spectral_cutoff
                                                if args.spectral_cutoff
                                                else ""),
                                            profile_type=profile_type,
                                            profile_direction=profile_direction,
                                            time_ms=f"{dt_ms:.1f}",
                                        )
                                        row.update(cfg_cols)
                                        row.update(const_cols(w))
                                        writer.writerow(row)
                                        total_combos += 1
                                        eq_total += 1
                                        if eval_result["exact_recovery"]:
                                            eq_exact_count += 1
                                        if total_combos % 50 == 0:
                                            f.flush()
                                            print(f"    [{eq_name}] "
                                                  f"{total_combos} combos",
                                                  flush=True)
                    f.flush()

            # ── Per-equation summary ─────────────────────────────────────
            eq_summaries.append(
                f"  {eq_name:25s}  exact={eq_exact_count}/{eq_total}  "
                f"({100 * eq_exact_count / max(eq_total, 1):.0f}%)"
            )
            print(f"  [{eq_name}] done: {eq_total} rows, "
                  f"exact={eq_exact_count}", flush=True)

    # ── Final summary ────────────────────────────────────────────────────
    elapsed_total = time.time() - t_start_total
    print(f"\n{'='*70}")
    print("Per-equation exact recovery summary:")
    for s in eq_summaries:
        print(s)
    print(f"\nDone. {total_combos} total rows → {out_path}")
    print(f"Total wall time: {elapsed_total:.1f}s ({elapsed_total/60:.1f} min)")

    # ── Write timing log ────────────────────────────────────────────────
    log_path = out_path.replace(".csv", ".log")
    omp_threads = os.environ.get("OMP_NUM_THREADS", "not set")
    try:
        cpu_info = subprocess.check_output(
            ["lscpu"], stderr=subprocess.DEVNULL, text=True
        )
        cpu_model = [l for l in cpu_info.splitlines() if "Model name" in l]
        cpu_model = cpu_model[0].split(":", 1)[1].strip() if cpu_model else "unknown"
    except Exception:
        cpu_model = "unknown"
    with open(log_path, "w") as lf:
        lf.write(f"# Sweep timing log\n")
        lf.write(f"csv_file: {os.path.basename(out_path)}\n")
        lf.write(f"total_wall_time_s: {elapsed_total:.1f}\n")
        lf.write(f"total_wall_time_min: {elapsed_total/60:.1f}\n")
        lf.write(f"total_rows: {total_combos}\n")
        lf.write(f"omp_num_threads: {omp_threads}\n")
        lf.write(f"hostname: {platform.node()}\n")
        lf.write(f"cpu_model: {cpu_model}\n")
        lf.write(f"python_version: {sys.version.split()[0]}\n")
        lf.write(f"timestamp: {timestamp}\n")
        lf.write(f"equations: {args.equations or 'all'}\n")
        lf.write(f"eval_mode: {args.eval_mode}\n")
    print(f"Timing log → {log_path}")


if __name__ == "__main__":
    main()
