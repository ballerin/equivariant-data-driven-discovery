"""
Canonical ablation sweep: compare pointwise / weak (w/wo IBP) across equations,
noise levels, libraries, backends, and regression algorithms.

Loads pre-generated .npz data files (produced by generate_data.py) and runs
systematic hyperparameter sweeps for every method combination.  Results are
evaluated with the universal pointwise-Theta metric (evaluation.py).

Supports optional post-hoc noise application and spectral low-pass
preprocessing for additional ablation dimensions.

Usage
-----
    # Core equivariant vs standard comparison
    python scripts/sweep.py --methods pointwise --algorithms stridge

    # Noise robustness (weak vs pointwise)
    python scripts/sweep.py \\
        --equations swift_hohenberg ks_chaotic \\
        --noise 0.02 0.05 \\
        --methods pointwise weak_direct weak_ibp \\
        --libraries equivariant

    # With post-hoc noise and spectral low-pass
    python scripts/sweep.py --post-hoc-noise 0.02 0.05 --spectral-cutoff 3.0

    # Full sweep
    python scripts/sweep.py
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
    build_eval_theta_jet, evaluate_unified,
)
from equivariant_pde_find.library.conversions import (  # noqa: E402
    build_wsindy_to_jet_map, _MULTIGRAPH_TO_JET_DECOMPOSITION,
    collect_needed_jet_monomials,
)
from equivariant_pde_find.equations import (  # noqa: E402
    find_op, true_allen_cahn, true_swift_hohenberg, true_ks,
    true_cahn_hilliard, true_pme, mask_pme, EQUATION_META,
    normalize_eq_name, CH_EPSILON,
)
from equivariant_pde_find.library import (  # noqa: E402
    build_multigraph_library, build_pdefind_library,
    build_polyjet_library, build_wsindy_library,
    build_mod_multigraph_library,
)
from equivariant_pde_find.utils import spectral_lowpass  # noqa: E402
from equivariant_pde_find.regression import TrainSTRidge, TrainMSTLS, TrainRFE, TrainExhaustive  # noqa: E402
from equivariant_pde_find.simulations import make_regression_data  # noqa: E402
from equivariant_pde_find.library import weak_form as wf  # noqa: E402

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
WEAK_M = [4, 8, 12]
# Test-function degree: we use a fixed p = WEAK_P = 6 for both direct and
# IBP weak forms.
#
# This used to be `WEAK_P = 6` for every m, which left phi^(4) large at the
# edge of the support and put a 3 % (m=12) to 84 % (m=8) systematic error
# into the fourth-order weak-IBP columns -- the leading terms of CH, SH and
# KS.
#
# The Messenger & Bortz (2021) Algorithm 1 prescribes p via
#   select_test_function_degree(m, max_deriv=4, tau=1e-10)
# which gives p=13 for m=12.  This satisfies a quadrature-accuracy decay
# condition (phi^(4) < 1e-10 at the support boundary), yielding < 0.001 %
# IBP discretisation error for Delta^2 f.
#
# However, p=13 produces a narrow test function (effective width ~4 grid
# points) whose columns for different differential operators are nearly
# collinear locally, destroying sparse recovery.  At p=6 the IBP
# discretisation error is 0.27 % — still 20x below the 5 % observation
# noise floor — but the wider effective support (~6 grid points) gives
# each operator a more distinctive spatial signature, which is decisive
# for STRidge/RFE to recover the correct sparse PDE.
#
# We therefore use a fixed p=6 for both the direct and IBP weak forms.
WEAK_P = 6              # test-function degree (direct and IBP)
SAMPLE_STRIDE = 4

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


def weak_p(m, use_ibp=True):
    """Test-function degree for support *m*.

    Returns a fixed p=6 for both direct and IBP weak forms.  See the
    discussion at the WEAK_P constant above for the quadrature-accuracy
    vs sparse-recovery trade-off.
    """
    return WEAK_P


def build_weak_form(f_snaps, ft_snaps, idx_sets, Lx, Ly, library, backend,
                 m, use_ibp):
    """Build weak-form (G, b) for a given library and backend."""
    WEAK_P = weak_p(m, use_ibp)
    with derivatives.use_backend(backend):
        rows_G, rows_b = [], []
        for s in range(len(f_snaps)):
            f = f_snaps[s]
            ft = ft_snaps[s]
            if library == "multigraph":
                G, desc = build_multigraph_library(
                    f, Lx, Ly, max_p=2, max_n=3,
                    formulation="weak", m=m, p=WEAK_P,
                    use_ibp=use_ibp,
                )
            elif library == "mod_multigraph":
                G, desc = build_mod_multigraph_library(
                    f, Lx, Ly, max_p=2, max_n=3,
                    formulation="weak", m=m, p=WEAK_P,
                    use_ibp=use_ibp,
                )
            elif library == "polyjet":
                G, desc = build_polyjet_library(
                    f, Lx, Ly, max_order=4, max_degree=3,
                    formulation="weak", m=m, p=WEAK_P, use_ibp=use_ibp,
                )
            elif library == "wsindy":
                G, desc = build_wsindy_library(
                    f, Lx, Ly, max_order=4, max_poly=3,
                    formulation="weak", m=m, p=WEAK_P, use_ibp=use_ibp,
                )
            else:
                G, desc = build_pdefind_library(
                    f, Lx, Ly, max_order=4, max_poly=3,
                    formulation="weak", m=m, p=WEAK_P, use_ibp=use_ibp,
                )
            b_w = wf.build_weak_rhs(ft, m=m, p=WEAK_P)

            Nx, Ny = f.shape
            Nv_x, Nv_y = Nx - 2 * m, Ny - 2 * m
            ix = np.arange(0, Nv_x, SAMPLE_STRIDE)
            iy = np.arange(0, Nv_y, SAMPLE_STRIDE)
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


# ── sweep one ───────────────────────────────────────────────────────────────


def sweep_one(Theta, Ft, true_eq, lam, d_tol, algorithm, mstls_alpha=0.5, mstls_M=None):
    """Run regression and return (n_active, w, lam_out)."""
    if algorithm == "mstls":
        # MSTLS does its own internal λ sweep; lam/d_tol are ignored
        w, lam_out, _ = TrainMSTLS(
            Theta, Ft, lambdas=None, M=mstls_M, n_lambdas=MSTLS_N_LAMBDAS,
            maxit=MSTLS_MAXITS, alpha=mstls_alpha,
        )
    elif algorithm == "rfe":
        # lam is the ridge regularisation for RFE (0 = pure OLS)
        w = TrainRFE(Theta, Ft, lam=lam)
        lam_out = lam
    else:
        w = TrainSTRidge(Theta, Ft, lam, d_tol)  # STRidge = Sequential Threshold Ridge
        lam_out = lam
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
        default=["pointwise", "weak_direct", "weak_ibp"],
        choices=["pointwise", "weak_direct", "weak_ibp"],
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
        "--backends", nargs="+",
        default=BACKENDS,
        choices=BACKENDS,
    )
    parser.add_argument(
        "--weak-m", nargs="+", type=int,
        default=WEAK_M,
        help="Test-function half-support values for weak-form methods.",
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

    methods = args.methods
    algorithms = args.algorithms
    libraries = args.libraries
    backends = args.backends
    weak_m_vals = args.weak_m

    # At small m the grid cannot resolve a degree-p test function at all.
    # Drop those m for weak_ibp; the direct weak form only needs φ itself
    # (no derivatives) and tolerates coarser resolution.
    usable_ibp_m_vals = [
        m for m in weak_m_vals
        if wf.test_function_is_resolved(m, WEAK_P, max_deriv=4)
    ]
    dropped = [m for m in weak_m_vals if m not in usable_ibp_m_vals]
    if dropped and "weak_ibp" in methods:
        print(
            f"  [weak_ibp] skipping m={dropped}: 2m+1 grid points cannot "
            f"resolve the degree-{WEAK_P} test function "
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
    tag = f"_{args.name}" if args.name else ""
    out_path = os.path.join(RESULTS_DIR, f"sweep_{timestamp}{tag}.csv")
    print(f"Output → {out_path}", flush=True)

    fieldnames = [
        "equation", "noise", "method", "algorithm", "library", "backend",
        "weak_m", "lam", "d_tol", "mstls_alpha",
        "n_active", "fp", "fn", "exact_recovery", "E_coeff", "max_coeff_err",
        "E_data", "cond_theta", "n_samples",
        "post_hoc_noise", "spectral_cutoff",
        "profile_type", "profile_direction",
        "time_ms",
    ]

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

            # If noisy dataset, load clean counterpart for evaluation
            # (Option B: evaluate against clean data so that E_coeff
            #  is comparable across noise levels)
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
                if os.path.exists(clean_path):
                    _clean_data = np.load(clean_path)
                    snapshots_clean_eval = _clean_data["snapshots"]
                    print(f"    clean eval: {os.path.basename(clean_path)}")
                else:
                    print(f"    WARNING: clean counterpart missing: {clean_path}")
                    snapshots_clean_eval = None
            else:
                snapshots_clean_eval = None  # already clean, no separate copy

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

                    # Build ground truth with appropriate parameters
                    _, desc_eq_ref = build_multigraph_library(
                        snapshots_proc[0], Lx, Ly,
                    )
                    if is_ch:
                        true_eq = true_fn(desc_eq_ref, eps=ch_eps)
                    elif eq_name == "ac":
                        true_eq = true_fn(desc_eq_ref, eps=0.15)
                    elif eq_name == "sh":
                        true_eq = true_fn(desc_eq_ref, r=0.3)
                    else:
                        true_eq = true_fn(desc_eq_ref)

                    # Prepare regression data
                    f_snaps_raw, ft_snaps_raw, idx_sets_raw, _ = \
                        make_regression_data(
                            snapshots_proc, times, n_per_snap=200,
                            skip_fraction=skip_frac, seed=0,
                        )

                    # Apply mask if needed
                    if mask_fn is not None:
                        rng = np.random.default_rng(1)
                        new_idx = []
                        for s in range(len(f_snaps_raw)):
                            valid = np.where(mask_fn(f_snaps_raw[s]))[0]
                            if len(valid) <= 200:
                                new_idx.append(valid)
                            else:
                                new_idx.append(np.sort(
                                    rng.choice(valid, size=200,
                                              replace=False)))
                        f_snaps = f_snaps_raw
                        ft_snaps = ft_snaps_raw
                        idx_sets = new_idx
                    else:
                        f_snaps = f_snaps_raw
                        ft_snaps = ft_snaps_raw
                        idx_sets = idx_sets_raw

                    n_snaps = len(f_snaps)
                    print(f"    interior snaps: {n_snaps}  "
                          f"points: {n_snaps * 200}", flush=True)

                    # ── Jet descriptions & ground truth (backend-independent) ─
                    # Build pdefind descriptions upfront (backend-independent:
                    # column names are the same regardless of spectral/polynomial)
                    _, desc_pdefind = build_pdefind_library(
                        snapshots_proc[0], Lx, Ly, max_order=4, max_poly=3,
                    )
                    # Pre-build WSINDy descriptions
                    from equivariant_pde_find.library import \
                        build_wsindy_library as _bpw
                    _, desc_ws = _bpw(
                        snapshots_proc[0], Lx, Ly, max_order=4, max_poly=3,
                    )
                    # Dynamic jet: union of all library decompositions
                    desc_jet = collect_needed_jet_monomials(
                        desc_mg=desc_eq_ref, desc_pdefind=desc_pdefind,
                        desc_ws=desc_ws,
                    )
                    # Ground truth in jet basis (analytical — no backend)
                    true_jet = np.zeros(len(desc_jet))
                    for j, (name, coeff) in enumerate(zip(desc_eq_ref, true_eq)):
                        if abs(coeff) < 1e-12:
                            continue
                        if name not in _MULTIGRAPH_TO_JET_DECOMPOSITION:
                            continue
                        for jet_name, mult in _MULTIGRAPH_TO_JET_DECOMPOSITION[name]:
                            true_jet[desc_jet.index(jet_name)] += coeff * mult

                    M_ws_to_jet = build_wsindy_to_jet_map(desc_ws, desc_jet)

                    # ── Per-backend evaluation Theta ─────────────────────
                    t_eval = time.time()
                    # Evaluate against clean data when available (Option B)
                    eval_snaps = (
                        snapshots_clean_eval
                        if snapshots_clean_eval is not None
                        else snapshots_proc
                    )
                    eval_thetas = {}
                    for b in backends:
                        Tj, Ftj = build_eval_theta_jet(
                            eval_snaps, times, Lx, Ly,
                            skip_fraction=skip_frac, seed=0, backend=b,
                            jet_monomials=desc_jet,
                        )
                        eval_thetas[b] = (Tj, Ftj)
                    Theta_eval, Ft_eval = eval_thetas[backends[0]]
                    print(f"    eval Theta: {Theta_eval.shape}  "
                          f"active_true="
                          f"{int(np.sum(np.abs(true_jet) > 1e-10))}"
                          f"  ({time.time() - t_eval:.1f}s)", flush=True)

                    # ── Pre-build method matrices ────────────────────────
                    prebuilt = {}

                    # Oracle-only: only need multigraph matrices
                    _libs = ["multigraph"] if oracle_only else libraries
                    for method, library, backend in itertools.product(
                        methods, _libs, backends
                    ):
                        if method == "pointwise":
                            m_vals = [None]
                        elif method == "weak_ibp":
                            m_vals = usable_ibp_m_vals
                        else:
                            m_vals = weak_m_vals

                        for m_val in m_vals:
                            if method == "weak_direct":
                                use_ibp_flag = False
                            elif method == "weak_ibp":
                                use_ibp_flag = True
                            else:
                                use_ibp_flag = False

                            actual_key = (method, library, backend,
                                          m_val if m_val is not None else 0)

                            if actual_key in prebuilt:
                                continue

                            t0 = time.time()
                            print(f"    Building: method={method} "
                                  f"lib={library} backend={backend}"
                                  + (f" m={m_val}" if m_val else ""),
                                  flush=True)

                            if method == "pointwise":
                                Theta, Ft, desc = build_pointwise(
                                    f_snaps, ft_snaps, idx_sets, Lx, Ly,
                                    library, backend,
                                )
                            else:
                                Theta, Ft, desc = build_weak_form(
                                    f_snaps, ft_snaps, idx_sets, Lx, Ly,
                                    library, backend, m_val, use_ibp_flag,
                                )

                            cond_val = cond_number(Theta)
                            prebuilt[actual_key] = (Theta, Ft, desc, cond_val)
                            elapsed = time.time() - t0
                            print(f"      shape={Theta.shape}  "
                                  f"cond={cond_val:.2e}  "
                                  f"({elapsed:.1f}s)", flush=True)

                    # ── Oracle: per-configuration ground truth ──────────
                    if compute_oracle:
                        for actual_key, (Theta, Ft, desc, cond_val) in \
                                prebuilt.items():
                            method, library, backend, m_val = actual_key
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
                                weak_m=(m_val if m_val != 0 else ""),
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
                            writer.writerow(row)
                            total_combos += 1

                    if oracle_only:
                        print(f"  [ac] oracle done", flush=True)
                        continue

                    # ── Sweep hyperparams ──────────────────────────────
                    for method, library, backend in itertools.product(
                        methods, libraries, backends
                    ):
                        if method == "pointwise":
                            m_vals = [0]
                        elif method == "weak_direct":
                            m_vals = weak_m_vals
                        else:
                            m_vals = usable_ibp_m_vals

                        for m_val in m_vals:
                            actual_key = (method, library, backend, m_val)
                            if actual_key not in prebuilt:
                                continue
                            Theta, Ft, desc, cond_val = prebuilt[actual_key]

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
                                        )

                                        eval_result = evaluate_unified(
                                            w, library, desc, desc_jet,
                                            *eval_thetas[backend], true_jet,
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
                                            weak_m=(m_val if m_val != 0
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
                                            algorithm,
                                        )

                                        eval_result = evaluate_unified(
                                            w, library, desc, desc_jet,
                                            *eval_thetas[backend], true_jet,
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
                                            weak_m=(m_val if m_val != 0
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
                                        w = TrainExhaustive(
                                            Theta, Ft, max_k=max_k,
                                            l0_penalty=l0_pen,
                                        )
                                        n_active = int(np.count_nonzero(
                                            np.abs(w) > 1e-10))

                                        eval_result = evaluate_unified(
                                            w, library, desc, desc_jet,
                                            *eval_thetas[backend], true_jet,
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
                                            weak_m=(m_val if m_val != 0
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
                                            algorithm,
                                        )

                                        eval_result = evaluate_unified(
                                            w, library, desc, desc_jet,
                                            *eval_thetas[backend], true_jet,
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
                                            weak_m=(m_val if m_val != 0
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
    print(f"Timing log → {log_path}")


if __name__ == "__main__":
    main()
