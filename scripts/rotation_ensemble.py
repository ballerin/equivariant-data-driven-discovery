#!/usr/bin/env python3
"""
Rotation-ensemble PDE discovery — standalone exploration script.

Runs sparse regression on multiple rotated copies of the same data and
combines results via intersection or voting.  The hypothesis is that
false positives from noise decorrelate under rotation, so ensemble
combination suppresses them.

Usage
-----
    # Quick test on noisy Swift-Hohenberg
    OMP_NUM_THREADS=1 python scripts/rotation_ensemble.py \\
        --data data/swift_hohenberg_noise2pct.npz \\
        --libraries equivariant \\
        --methods pointwise \\
        --algorithms strridge \\
        --backends spectral

    # Full sweep: all libraries, methods, backends
    OMP_NUM_THREADS=1 python scripts/rotation_ensemble.py \\
        --data data/swift_hohenberg_noise2pct.npz \\
        --transforms R_0,R_90,R_180,R_270 \\
        --combine intersection voting3 voting4
"""

import argparse
import csv
import datetime
import os
import sys
import time

import numpy as np

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        desc = kwargs.get("desc", "")
        items = list(iterable)
        print(f"{desc} ({len(items)} items)...")
        return items

# ── path setup ──────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives  # noqa: E402
from equivariant_pde_find.evaluation import (  # noqa: E402
    build_eval_theta_jet,
    evaluate_unified,
)
from equivariant_pde_find.library_mapping import (  # noqa: E402
    _MULTIGRAPH_TO_JET_DECOMPOSITION,
)
from equivariant_pde_find.equation_configs import (  # noqa: E402
    find_op, true_allen_cahn, true_swift_hohenberg, true_ks,
    true_cahn_hilliard, true_pme, mask_pme, EQUATION_META,
)
from equivariant_pde_find.library import (  # noqa: E402
    build_multigraph_library, build_pdefind_library,
    build_polyjet_library, build_wsindy_library,
    build_mod_multigraph_library,
)
from equivariant_pde_find.regression import (  # noqa: E402
    TrainSTRidge, TrainMSTLS, TrainRFE,
)
from equivariant_pde_find.rotation_ensemble import (  # noqa: E402
    apply_transform_snapshots,
    get_rotation_map,
    get_reflection_map,
    combine_by_intersection,
    combine_by_voting,
    TRANSFORM_REGISTRY,
    _rotation_k,
)
from equivariant_pde_find.simulations import make_regression_data  # noqa: E402

RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyperparameter grids ────────────────────────────────────────────────────
LAMS = [1e-6, 1e-4, 1e-2, 1.0]
D_TOLS = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, 0.5, 1.0,
          1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0]
MSTLS_N_LAMBDAS = 50
MSTLS_ALPHAS = [0.0, 0.5, 1.0]
MSTLS_MAXITS = 50


# ── loading ─────────────────────────────────────────────────────────────────


def load_dataset(npz_path):
    """Load a .npz data file. Returns (snapshots, times, Lx, Ly, N, dt_snap,
    obs_noise, eq_name, profile_type, profile_direction)."""
    data = np.load(npz_path)
    snapshots = data["snapshots"]
    times = data["times"]
    Lx = float(data["Lx"])
    Ly = float(data["Ly"])
    N = int(data["N"])
    dt_snap = float(data["dt_snap"])
    obs_noise = float(data["obs_noise"])

    if "eq_name" in data:
        eq_name = str(data["eq_name"])
    else:
        eq_name = os.path.basename(npz_path).rsplit("_", 1)[0]
        for suffix in ["_clean", "_noise2pct", "_noise5pct"]:
            if eq_name.endswith(suffix):
                eq_name = eq_name[: -len(suffix)]
                break

    profile_type = str(data["profile_type"]) if "profile_type" in data else ""
    profile_direction = str(data["profile_direction"]) if "profile_direction" in data else ""

    return (snapshots, times, Lx, Ly, N, dt_snap, obs_noise, eq_name,
            profile_type, profile_direction)


def _get_meta(name):
    """Return (true_fn, skip_fraction, mask_fn) for a dataset name."""
    if name in EQUATION_META:
        return EQUATION_META[name]
    for prefix in sorted(EQUATION_META, key=len, reverse=True):
        if name.startswith(prefix):
            return EQUATION_META[prefix]
    raise KeyError(f"No metadata for equation: {name}")


# ── library building ─────────────────────────────────────────────────────────


def build_library(f, Lx, Ly, library):
    """Build a pointwise library from a single snapshot."""
    if library == "multigraph":
        T, desc = build_multigraph_library(f, Lx, Ly, max_p=2, max_n=3)
    elif library == "mod_multigraph":
        T, desc = build_mod_multigraph_library(f, Lx, Ly, max_p=2, max_n=3)
    elif library == "polyjet":
        T, desc = build_polyjet_library(
            f, Lx, Ly, max_order=4, max_degree=3,
        )
    elif library in ("wsindy",):
        T, desc = build_wsindy_library(
            f, Lx, Ly, max_order=4, max_poly=3,
        )
    else:
        T, desc = build_pdefind_library(f, Lx, Ly, max_order=4, max_poly=3)
    return T, desc


# ── regression ──────────────────────────────────────────────────────────────


def run_regression(Theta, Ft, algorithm, lam, d_tol,
                   mstls_alpha=0.5, mstls_M=None):
    """Run regression and return coefficient vector."""
    if algorithm == "mstls":
        w, _, _ = TrainMSTLS(
            Theta, Ft, lambdas=None, M=mstls_M,
            n_lambdas=MSTLS_N_LAMBDAS,
            maxit=MSTLS_MAXITS, alpha=mstls_alpha,
        )
    elif algorithm == "rfe":
        w = TrainRFE(Theta, Ft)
    else:
        w = TrainSTRidge(Theta, Ft, lam, d_tol)
    return w


def cond_number(Theta):
    """Condition number (2-norm)."""
    return float(np.linalg.cond(Theta))


# ── main ────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Rotation-ensemble PDE discovery.",
    )
    parser.add_argument(
        "--data", required=True,
        help="Path to .npz data file.",
    )
    parser.add_argument(
        "--libraries", nargs="+",
        default=["multigraph"],
        choices=["multigraph", "mod_multigraph", "pdefind",
                 "polyjet", "wsindy"],
    )
    parser.add_argument(
        "--methods", nargs="+",
        default=["pointwise"],
        choices=["pointwise"],
        help="Only pointwise is supported in this script.",
    )
    parser.add_argument(
        "--algorithms", nargs="+",
        default=["strridge"],
        choices=["strridge", "mstls", "rfe"],
    )
    parser.add_argument(
        "--backends", nargs="+",
        default=["spectral", "polynomial"],
        choices=["spectral", "polynomial"],
    )
    parser.add_argument(
        "--transforms",
        default="R_0,R_90,R_180,R_270",
        help="Comma-separated transform names (default: R_0,R_90,R_180,R_270).",
    )
    parser.add_argument(
        "--combine", nargs="+",
        default=["intersection", "voting3", "voting4"],
        help="Combine strategies: intersection, voting3, voting4.",
    )
    parser.add_argument(
        "--name", default="",
        help="Optional tag appended to output filename.",
    )
    parser.add_argument(
        "--n-per-snap", type=int, default=200,
        help="Number of regression points per snapshot.",
    )
    parser.add_argument(
        "--d-tols", nargs="+", type=float, default=None,
        help="Override default D_TOLS grid for STRidge.",
    )
    parser.add_argument(
        "--lams", nargs="+", type=float, default=None,
        help="Override default LAMS grid for STRidge.",
    )
    parser.add_argument(
        "--mstls-alphas", nargs="+", type=float, default=None,
        help="Override default MSTLS alpha values.",
    )
    args = parser.parse_args()

    # Parse transform list
    transform_names = [t.strip() for t in args.transforms.split(",")]
    for t in transform_names:
        if t not in TRANSFORM_REGISTRY:
            print(f"Warning: unknown transform {t!r}, skipping.")
    transform_names = [t for t in transform_names if t in TRANSFORM_REGISTRY]
    if not transform_names:
        raise ValueError("No valid transforms specified.")

    n_ensemble = len(transform_names)
    print(f"Ensemble size: {n_ensemble}  ({', '.join(transform_names)})")

    # Allow override of hyperparameter grids
    d_tols_grid = args.d_tols if args.d_tols is not None else D_TOLS
    lams_grid = args.lams if args.lams is not None else LAMS
    mstls_alphas_grid = args.mstls_alphas if args.mstls_alphas is not None else MSTLS_ALPHAS

    # Parse combine strategies
    combine_strategies = []
    for c in args.combine:
        if c == "intersection":
            combine_strategies.append(("intersection", 0))
        elif c.startswith("voting"):
            threshold = int(c.replace("voting", ""))
            combine_strategies.append(("voting", threshold))
        else:
            print(f"Warning: unknown combine strategy {c!r}, skipping.")
    if not combine_strategies:
        raise ValueError("No valid combine strategies specified.")

    # ── load data ───────────────────────────────────────────────────────
    npz_path = args.data
    print(f"\nLoading: {npz_path}")
    snapshots_orig, times, Lx, Ly, N, dt_snap, obs_noise, eq_name, \
        profile_type, profile_direction = load_dataset(npz_path)
    print(f"  eq={eq_name}  noise={obs_noise:.3f}  "
          f"snaps={len(snapshots_orig)}  grid={N}×{N}")

    true_fn, skip_frac, mask_fn = _get_meta(eq_name)

    # ── pre-build per-backend evaluation Theta (polynomial jet) ─────────
    print("Pre-building evaluation Theta (polynomial jet, all points)...", end=" ", flush=True)
    t0 = time.time()
    eval_thetas = {}
    for b in args.backends:
        Tj, Ftj = build_eval_theta_jet(
            snapshots_orig, times, Lx, Ly,
            skip_fraction=skip_frac, seed=0, backend=b,
        )
        eval_thetas[b] = (Tj, Ftj)
    from equivariant_pde_find.library import build_polyjet_library as _bpl
    _, desc_jet = _bpl(snapshots_orig[0], Lx, Ly, max_order=4, max_degree=3)
    print(f"done ({time.time() - t0:.1f}s)")

    # Build ground truth in jet basis
    _, desc_eq_ref = build_multigraph_library(snapshots_orig[0], Lx, Ly,
                                                max_p=2, max_n=3)
    true_eq = true_fn(desc_eq_ref)
    true_jet = np.zeros(len(desc_jet))
    for j, (name, coeff) in enumerate(zip(desc_eq_ref, true_eq)):
        if abs(coeff) < 1e-12:
            continue
        if name not in _MULTIGRAPH_TO_JET_DECOMPOSITION:
            continue
        for jet_name, mult in _MULTIGRAPH_TO_JET_DECOMPOSITION[name]:
            true_jet[desc_jet.index(jet_name)] += coeff * mult

    # ── prepare regression data from original ───────────────────────────
    f_snaps_orig, ft_snaps_orig, idx_sets_orig, _ = make_regression_data(
        snapshots_orig, times, n_per_snap=args.n_per_snap,
        skip_fraction=skip_frac, seed=0,
    )

    # ── output CSV ──────────────────────────────────────────────────────
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"_{args.name}" if args.name else ""
    out_path = os.path.join(
        RESULTS_DIR, f"rotation_ensemble_{timestamp}{tag}.csv",
    )
    print(f"Output → {out_path}\n")

    fieldnames = [
        "equation", "noise", "library", "backend", "algorithm",
        "lam", "d_tol", "mstls_alpha",
        "n_ensemble", "transforms", "combine_strategy",
        "n_active", "fp", "fn", "exact_recovery", "E", "max_coeff_err",
        "E_fit", "cond_theta_orig",
        "time_ms",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        f.flush()

        t_start = time.time()
        total_combos = 0

        for library_key in args.libraries:
            for backend in args.backends:
                with derivatives.use_backend(backend):
                    # Build Theta_orig once per (library, backend)
                    t_build = time.time()
                    Theta_orig, desc_lib = build_library(
                        snapshots_orig[0], Lx, Ly, library_key,
                    )
                    # Assemble full regression Theta/Ft from original
                    rows_T, rows_F = [], []
                    for s in range(len(f_snaps_orig)):
                        T_s, _ = build_library(
                            f_snaps_orig[s], Lx, Ly, library_key,
                        )
                        rows_T.append(T_s[idx_sets_orig[s]])
                        rows_F.append(
                            ft_snaps_orig[s].ravel()[idx_sets_orig[s]],
                        )
                    Theta_full_orig = np.vstack(rows_T).astype(float)
                    Ft_full_orig = np.concatenate(rows_F).astype(float)
                    cond_orig = cond_number(Theta_full_orig)
                    build_time = time.time() - t_build

                    # ── Pre-build Theta/Ft for each transform ────────
                    Theta_transforms = {}
                    Ft_transforms = {}
                    map_matrices = {}

                    for tname in transform_names:
                        t_field = time.time()
                        if tname == "R_0":
                            Theta_transforms[tname] = Theta_full_orig
                            Ft_transforms[tname] = Ft_full_orig
                            map_matrices[tname] = np.eye(
                                len(desc_lib), dtype=float,
                            )
                        else:
                            # Apply transform to snapshots
                            snaps_t = apply_transform_snapshots(
                                snapshots_orig, tname,
                            )
                            f_s, ft_s, idx_s, _ = make_regression_data(
                                snaps_t, times,
                                n_per_snap=args.n_per_snap,
                                skip_fraction=skip_frac, seed=0,
                            )
                            rows_Tt, rows_Ft = [], []
                            for s_idx in range(len(f_s)):
                                T_s, _ = build_library(
                                    f_s[s_idx], Lx, Ly, library_key,
                                )
                                rows_Tt.append(T_s[idx_s[s_idx]])
                                rows_Ft.append(
                                    ft_s[s_idx].ravel()[idx_s[s_idx]],
                                )
                            Theta_transforms[tname] = np.vstack(
                                rows_Tt,
                            ).astype(float)
                            Ft_transforms[tname] = np.concatenate(
                                rows_Ft,
                            ).astype(float)

                            # Get the rotation/reflection map
                            k = _rotation_k(tname)
                            if k is not None:
                                map_matrices[tname] = get_rotation_map(
                                    library_key, desc_lib, k,
                                )
                            else:
                                map_matrices[tname] = get_reflection_map(
                                    library_key, desc_lib, tname,
                                )
                        t_field_time = time.time() - t_field

                    # Determine hyperparameter grid
                    if "mstls" in args.algorithms:
                        alphas = mstls_alphas_grid
                    else:
                        alphas = [0.5]  # dummy for STRidge/RFE

                    # Build iteration list
                    iterations = []
                    for alg in args.algorithms:
                        if alg == "strridge":
                            for lam in lams_grid:
                                for d_tol in d_tols_grid:
                                    iterations.append(
                                        (alg, lam, d_tol, 0.5),
                                    )
                        elif alg == "mstls":
                            for alpha in alphas:
                                iterations.append(
                                    (alg, 0.0, 0.0, alpha),
                                )
                        elif alg == "rfe":
                            iterations.append(("rfe", 0.0, 0.0, 0.5))

                    print(
                        f"  lib={library_key}  backend={backend}  "
                        f"cond(Θ)={cond_orig:.2e}  "
                        f"combos={len(iterations)}×{n_ensemble}",
                        flush=True,
                    )

                    for alg, lam, d_tol, alpha in tqdm(
                        iterations, desc=f"    {library_key}/{backend}/{alg}",
                    ):
                        t_combo = time.time()

                        # ── Run regression on each transform ──────
                        w_list_mapped = []
                        for tname in transform_names:
                            Theta_t = Theta_transforms[tname]
                            Ft_t = Ft_transforms[tname]
                            P = map_matrices[tname]

                            w_rot = run_regression(
                                Theta_t, Ft_t, alg, lam, d_tol,
                                mstls_alpha=alpha,
                            )
                            w_mapped = P @ w_rot
                            w_list_mapped.append(w_mapped)

                        # ── Combine ───────────────────────────────
                        for combine_name, threshold in combine_strategies:
                            if combine_name == "intersection":
                                w_c, mask_c = combine_by_intersection(
                                    w_list_mapped,
                                    Theta_full_orig,
                                    Ft_full_orig,
                                )
                            else:
                                w_c, mask_c = combine_by_voting(
                                    w_list_mapped,
                                    Theta_full_orig,
                                    Ft_full_orig,
                                    threshold=threshold,
                                )

                            n_active = int(np.sum(mask_c))

                            # Evaluate in jet basis
                            result = evaluate_unified(
                                w_c, library_key, desc_lib, desc_jet,
                                *eval_thetas[backend], true_jet,
                            )
                            fp = result["fp"]
                            fn = result["fn"]
                            E = result["E"]
                            exact = result["exact_recovery"]

                            # E_fit: residual on original regression data
                            if n_active > 0:
                                E_fit = float(
                                    np.linalg.norm(
                                        Theta_full_orig @ w_c
                                        - Ft_full_orig,
                                    )
                                    / np.linalg.norm(Ft_full_orig)
                                )
                                max_coeff_err = float(
                                    np.max(np.abs(w_c[mask_c]))
                                )
                            else:
                                E_fit = 1.0
                                max_coeff_err = 0.0

                            elapsed = (time.time() - t_combo) * 1000

                            writer.writerow({
                                "equation": eq_name,
                                "noise": f"{obs_noise:.3f}",
                                "library": library_key,
                                "backend": backend,
                                "algorithm": alg,
                                "lam": f"{lam:.6f}",
                                "d_tol": f"{d_tol:.6f}",
                                "mstls_alpha": f"{alpha:.6f}",
                                "n_ensemble": n_ensemble,
                                "transforms": ",".join(transform_names),
                                "combine_strategy": (
                                    combine_name if combine_name == "intersection"
                                    else f"voting{threshold}"
                                ),
                                "n_active": n_active,
                                "fp": fp,
                                "fn": fn,
                                "exact_recovery": 1 if exact else 0,
                                "E": f"{E:.6e}",
                                "max_coeff_err": f"{max_coeff_err:.6e}",
                                "E_fit": f"{E_fit:.6e}",
                                "cond_theta_orig": f"{cond_orig:.6e}",
                                "time_ms": f"{elapsed:.1f}",
                            })
                            f.flush()
                            total_combos += 1

        elapsed_total = time.time() - t_start
        print(f"\nDone. {total_combos} combinations in {elapsed_total:.1f}s")


if __name__ == "__main__":
    main()
