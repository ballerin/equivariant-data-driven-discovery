#!/usr/bin/env python
"""n-dimensional PDE discovery sweep.

Compares STRidge and RFE across four library types (multigraph,
mod_multigraph, pdefind, wsindy), equations, dimensions, noise levels,
and IBP on/off.  Uses weak-form with the spectral derivative backend.

Records compute stats (peak RSS, build/fit time, grid size) for
modal.com deployment planning.

Usage
-----
    python scripts/sweep_nd.py
    python scripts/sweep_nd.py --equations ac --ndim 2 3 4 5
    python scripts/sweep_nd.py --libraries multigraph wsindy --algorithms stride
    python scripts/sweep_nd.py --noise 0.0 0.01 --ibp both
"""

import argparse
import csv
import datetime
import os
import resource
import sys
import time

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives
from equivariant_pde_find.library import (
    build_multigraph_library,
    build_mod_multigraph_library,
    build_pdefind_library,
    build_wsindy_library,
)
from equivariant_pde_find.library import weak_form as wf
from equivariant_pde_find.regression import TrainSTRidge, TrainRFE, SolverCache
from equivariant_pde_find.simulations import make_regression_data
from equivariant_pde_find.equations import (
    CH_EPSILON,
    EQUATION_META,
    true_allen_cahn,
    true_cahn_hilliard,
    true_swift_hohenberg,
    true_pme,
    true_allen_cahn_mod_multigraph,
    true_cahn_hilliard_mod_multigraph,
    true_swift_hohenberg_mod_multigraph,
    true_pme_mod_multigraph,
    true_allen_cahn_pdefind,
    true_cahn_hilliard_pdefind,
    true_swift_hohenberg_pdefind,
    true_pme_pdefind,
    true_allen_cahn_wsindy,
    true_swift_hohenberg_wsindy,
    true_pme_wsindy_nd,
    true_cahn_hilliard_wsindy,
)

DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyper-parameters ────────────────────────────────────────────────────────
from equivariant_pde_find.library.weak_form import (  # noqa: E402
    test_function_is_resolved as _wf_is_resolved,
)

# Test-function support m scales with grid resolution N.
# Key constraint: N > 2*m for any valid interior points.
M_VALS_BY_NDIM = {2: [12], 3: [10], 4: [8], 5: [6]}

# Test-function degree is fixed at 6 for both direct and IBP:
# p=6 gives wider effective support (~6 grid points) for better sparse
# recovery, and the IBP discretisation error (~0.27 %) is negligible.
SPACE_P_FIXED = 6


def space_p(m, space_ibp=True):
    """Test-function degree for support *m* (fixed p=6)."""
    return SPACE_P_FIXED


def peak_rss_mb():
    """Peak resident set size in MB (Linux/macOS)."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


N_PER_SNAP = 200       # spatial samples per snapshot
SPACE_STRIDE = 4
LAMS = [1e-6, 1e-4, 1e-2, 1.0]       # STRidge ridge penalty
D_TOLS = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, 0.5,
          1.0, 2.0, 5.0, 10.0]
RFE_LAMS = [0.0]                        # RFE: pure OLS only

# ── library builders ─────────────────────────────────────────────────────────

def _build_space_weak_nd(f_snaps, ft_snaps, idx_sets, Lx, Ly, Lz,
                   library_fn, m, p, space_ibp):
    """Build weak-form (G, b) for n-d data with a given library builder.

    Handles both 2-D (Lz=None) and N-D (Lz=float or None with ndim>3)
    automatically.  Dispatches to the correct builder signature for
    multigraph (max_p/max_n), pdefind/wsindy (max_order/max_poly).
    """
    kernel = wf._cached_kernel_1d(m, p)
    rows_G, rows_b = [], []

    for s in range(len(f_snaps)):
        f = f_snaps[s]
        ft = ft_snaps[s]
        ndim = f.ndim

        if library_fn in (build_pdefind_library, build_wsindy_library):
            G, desc = library_fn(
                f, Lx, Ly, max_order=4, max_poly=3,
                formulation="space_weak", m=m, p=p, space_ibp=space_ibp, Lz=Lz,
            )
        else:
            G, desc = library_fn(
                f, Lx, Ly, max_p=2, max_n=3,
                formulation="space_weak", m=m, p=p, space_ibp=space_ibp, Lz=Lz,
            )

        # Build weak RHS (ravels via build_weak_rhs for flat-index compatibility).
        b_w = wf.build_weak_rhs(ft, kernel_x=kernel, kernel_y=kernel,
                                kernel_z=kernel, m=m, p=p)

        # Strided sampling of valid interior points (C-order / row-major).
        shape = f.shape
        Nv = tuple(s - 2 * m for s in shape)
        grid_axes = [np.arange(0, Nv[i], SPACE_STRIDE) for i in range(ndim)]
        mesh = np.meshgrid(*grid_axes, indexing="ij")
        # Flatten with C-order: last axis varies fastest.
        flat_idx = 0
        stride = 1
        for i in reversed(range(ndim)):
            flat_idx = flat_idx + mesh[i].ravel() * stride
            stride *= Nv[i]

        rng = np.random.default_rng(s)
        n_per = len(idx_sets[s])
        chosen = rng.choice(len(flat_idx),
                            size=min(n_per, len(flat_idx)), replace=False)
        q_idx = flat_idx[chosen]

        rows_G.append(G[q_idx])
        rows_b.append(b_w[q_idx])

    G_all = np.vstack(rows_G).astype(float)
    b_all = np.concatenate(rows_b).astype(float)
    return G_all, b_all, desc


# ── evaluation ───────────────────────────────────────────────────────────────

def evaluate_native(w, desc, true_dict):
    """Evaluate recovered coefficients against ground truth in native basis.

    Returns dict with keys: fp, fn, exact_recovery, E_coeff, n_active.
    """
    true_coeffs = np.zeros(len(desc))
    for i, d in enumerate(desc):
        if d in true_dict:
            true_coeffs[i] = true_dict[d]

    active_true = np.abs(true_coeffs) > 1e-12
    active_rec = np.abs(w) > 1e-10
    n_active = int(np.sum(active_rec))

    fp = int(np.sum(active_rec & ~active_true))
    fn = int(np.sum(active_true & ~active_rec))
    exact = int(fp == 0 and fn == 0)

    if np.any(active_true):
        E_coeff = float(np.linalg.norm(w - true_coeffs) /
                        max(np.linalg.norm(true_coeffs), 1e-14))
    else:
        E_coeff = float(np.linalg.norm(w)) if np.any(active_rec) else 0.0

    return {
        "n_active": n_active, "fp": fp, "fn": fn,
        "exact_recovery": exact, "E_coeff": E_coeff,
    }


def cond_number(Theta):
    """Condition number, capped to avoid overflow."""
    try:
        return float(np.linalg.cond(Theta))
    except np.linalg.LinAlgError:
        return float("nan")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="n-D PDE discovery sweep."
    )
    parser.add_argument(
        "--equations", nargs="+", default=["ac", "sh", "pme", "ch"],
        help="Equations to sweep (ac, sh, pme, ch).",
    )
    parser.add_argument(
        "--ndim", nargs="+", type=int, default=[2, 3, 4, 5],
        help="Spatial dimensions (default: 2 3 4 5).",
    )
    parser.add_argument(
        "--noise", nargs="+", type=float, default=[0.0, 0.02, 0.05],
        help="Noise levels (default: 0.0 0.02 0.05).",
    )
    parser.add_argument(
        "--libraries", nargs="+",
        default=["multigraph", "mod_multigraph", "pdefind", "wsindy"],
        help="Library types to test.",
    )
    parser.add_argument(
        "--algorithms", nargs="+", default=["stride", "rfe"],
        help="Regression algorithms (stride, rfe).",
    )
    parser.add_argument(
        "--space-ibp", "--ibp", dest="space_ibp_mode", default="both",
        choices=["both", "ibp", "direct"],
        help="Spatial IBP mode: both (default), ibp (space_ibp=True), direct "
             "(space_ibp=False).  --ibp is the old name.",
    )
    parser.add_argument(
        "--name", default="", help="Tag appended to output filename.",
    )
    parser.add_argument(
        "--out", default="", help="Exact output CSV path (overrides auto-generated name).",
    )
    parser.add_argument(
        "--data-dir", default="", help="Directory containing .npz datasets (default: ../data).",
    )
    args = parser.parse_args()

    derivatives.set_backend("spectral")

    # ── data directory ────────────────────────────────────────────────────
    data_dir = args.data_dir if args.data_dir else DATA_DIR

    # ── IBP modes ─────────────────────────────────────────────────────────
    if args.space_ibp_mode == "both":
        ibp_modes = [True, False]
    elif args.space_ibp_mode == "ibp":
        ibp_modes = [True]
    else:
        ibp_modes = [False]

    # ── library + truth mapping ──────────────────────────────────────────
    LIB_MAP = {
        "multigraph": {
            "builder": build_multigraph_library,
            "truth_ac": true_allen_cahn,
            "truth_sh": true_swift_hohenberg,
            "truth_pme": true_pme,
            "truth_ch": true_cahn_hilliard,
            "truth_fn_ndim": False,  # truth is dimension-agnostic
        },
        "mod_multigraph": {
            "builder": build_mod_multigraph_library,
            "truth_ac": true_allen_cahn_mod_multigraph,
            "truth_sh": true_swift_hohenberg_mod_multigraph,
            "truth_pme": true_pme_mod_multigraph,
            "truth_ch": true_cahn_hilliard_mod_multigraph,
            "truth_fn_ndim": False,
        },
        "pdefind": {
            "builder": build_pdefind_library,
            "truth_ac": true_allen_cahn_pdefind,
            "truth_sh": true_swift_hohenberg_pdefind,
            "truth_pme": true_pme_pdefind,
            "truth_ch": true_cahn_hilliard_pdefind,
            "truth_fn_ndim": True,  # truth depends on ndim
        },
        "wsindy": {
            "builder": build_wsindy_library,
            "truth_ac": true_allen_cahn_wsindy,
            "truth_sh": true_swift_hohenberg_wsindy,
            "truth_pme": true_pme_wsindy_nd,
            "truth_ch": true_cahn_hilliard_wsindy,
            "truth_fn_ndim": True,
        },
    }

    # ── discover datasets ────────────────────────────────────────────────
    noise_labels = {0.0: "clean", 0.01: "noise1pct", 0.02: "noise2pct", 0.05: "noise5pct"}
    eq_names = {
        "ac": "ac", "sh": "sh",
        "pme": "pme", "ch": "ch",
    }

    # Priority ordering for dataset filename patterns.
    TAG = "ndim_sweep"

    def _find_dataset(eq_name, ndim, noise):
        label = noise_labels[noise]
        # Tagged N-D file (required for all dimensions, including 2D).
        candidates = [
            f"{eq_name}_nd{ndim}_{TAG}_{label}.npz",
            f"{eq_name}_nd{ndim}_{label}.npz",
        ]
        for fname in candidates:
            path = os.path.join(data_dir, fname)
            if os.path.exists(path):
                return fname, path
        return None, None

    datasets = []
    for eq_key in args.equations:
        eq_name = eq_names[eq_key]
        for ndim in args.ndim:
            for noise in args.noise:
                fname, path = _find_dataset(eq_name, ndim, noise)
                if fname:
                    datasets.append((eq_key, ndim, noise, fname, path))
                else:
                    label = noise_labels[noise]
                    print(f"SKIP: {eq_name}_nd{ndim}_{label}.npz not found. "
                          f"Run generate_nd_data.py first.")

    if not datasets:
        print("No datasets found.")
        return

    # ── output ───────────────────────────────────────────────────────────
    if args.out:
        out_path = args.out
    else:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        tag = f"_{args.name}" if args.name else ""
        out_path = os.path.join(RESULTS_DIR, f"sweep_nd_{timestamp}{tag}.csv")

    fieldnames = [
        "equation", "ndim", "noise", "library", "algorithm",
        "space_ibp", "m", "d_tol", "lam",
        "n_active", "fp", "fn", "exact_recovery", "E_coeff",
        "E_data", "cond_theta",
        "n_cols", "n_rows", "n_snapshots", "grid_points",
        "build_time_ms", "fit_time_ms", "peak_rss_mb",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        f.flush()

        total = 0

        for eq_key, ndim, noise, fname, npz_path in datasets:
            print(f"\n{'='*60}")
            print(f"Dataset: {fname}  (eq={eq_key}, ndim={ndim}, noise={noise})")
            print(f"{'='*60}")

            # Load data
            data = np.load(npz_path, allow_pickle=True)
            snapshots = data["snapshots"]
            times = data["times"]
            Lx = float(data["Lx"])
            Ly = float(data["Ly"])
            Lz = float(data["Lz"]) if "Lz" in data else None

            # For legacy 2D data without ndim metadata
            file_ndim = int(data["ndim"]) if "ndim" in data else ndim
            effective_ndim = file_ndim
            grid_points = int(np.prod(snapshots.shape[1:]))
            n_snapshots = len(snapshots)

            # Per-equation skip_fraction (matching sweep.py)
            skip_frac = EQUATION_META[eq_key][1]

            f_snaps, ft_snaps, idx_sets, _ = make_regression_data(
                snapshots, times, n_per_snap=N_PER_SNAP,
                skip_fraction=skip_frac, seed=0,
            )
            print(f"  snaps={len(f_snaps)}, grid={snapshots.shape[1:]}")

            for lib_key in args.libraries:
                # PDE-FIND cannot represent PME (products of derivatives).
                if eq_key == "pme" and lib_key == "pdefind":
                    continue
                lib_info = LIB_MAP[lib_key]
                builder = lib_info["builder"]

                # ── ground truth ──────────────────────────────────────────
                # Per-equation parameters for truth functions.
                _EQ_PARAMS = {
                    "ac": {"eps": 0.15},
                    "sh": {"r": 0.3},
                    "ch": {"eps": CH_EPSILON},
                    "pme": {},
                }
                truth_fn = lib_info.get(f"truth_{eq_key}")
                if truth_fn is None:
                    print(f"  [{lib_key:>16s}] SKIP: no truth for {eq_key}")
                    continue

                if lib_info["truth_fn_ndim"]:
                    # PDE-FIND / WSINDy: truth_fn(ndim, **params) → dict
                    true_dict = truth_fn(effective_ndim, **_EQ_PARAMS.get(eq_key, {}))
                else:
                    # Multigraph / mod_multigraph: truth_fn(descriptions, **params) → array
                    m_tmp = M_VALS_BY_NDIM[effective_ndim][0]
                    p_tmp = space_p(m_tmp, space_ibp=True)
                    if builder in (build_pdefind_library, build_wsindy_library):
                        _, desc_tmp = builder(
                            f_snaps[0], Lx, Ly, max_order=4, max_poly=3,
                            formulation="space_weak", m=m_tmp, p=p_tmp,
                            space_ibp=True, Lz=Lz,
                        )
                    else:
                        _, desc_tmp = builder(
                            f_snaps[0], Lx, Ly, max_p=2, max_n=3,
                            formulation="space_weak", m=m_tmp, p=p_tmp,
                            space_ibp=True, Lz=Lz,
                        )
                    true_w = truth_fn(desc_tmp, **_EQ_PARAMS.get(eq_key, {}))
                    true_dict = {d: true_w[i] for i, d in enumerate(desc_tmp)
                                 if abs(true_w[i]) > 1e-12}

                # ── loop over test-function support ─────────────────────
                m_vals = M_VALS_BY_NDIM.get(effective_ndim, [12])
                for m_val in m_vals:
                    # ── loop over IBP modes ───────────────────────────
                    for space_ibp in ibp_modes:
                        p_val = space_p(m_val, space_ibp)

                        # Build weak-form library once per (m, ibp)
                        t_build_start = time.time()
                        G, b, desc = _build_space_weak_nd(
                            f_snaps, ft_snaps, idx_sets, Lx, Ly, Lz,
                            builder, m_val, p_val, space_ibp,
                        )
                        build_time_ms = (time.time() - t_build_start) * 1000
                        mem_mb = peak_rss_mb()
                        cond_val = cond_number(G)
                        n_cols = len(desc)
                        n_rows = G.shape[0]

                        ibp_label = "ibp" if space_ibp else "direct"
                        print(f"  [{lib_key:>16s}] m={m_val} {ibp_label:6s}  "
                              f"shape=({n_rows},{n_cols})  cond={cond_val:.2e}  "
                              f"build={build_time_ms:.0f}ms  "
                              f"rss={mem_mb:.0f}MB")

                        # Data-only solver work, shared by all settings on
                        # this matrix (results are unchanged).
                        solver_cache = SolverCache()
                        for algorithm in args.algorithms:
                            if algorithm == "stride":
                                for lam in LAMS:
                                    for d_tol in D_TOLS:
                                        t_fit = time.time()
                                        w = TrainSTRidge(G, b, lam, d_tol,
                                                         maxit=25, seed=0,
                                                         cache=solver_cache)
                                        fit_time_ms = (time.time() - t_fit) * 1000
                                        eval_r = evaluate_native(w, desc, true_dict)
                                        E_data = float(
                                            np.linalg.norm(G @ w - b) /
                                            max(np.linalg.norm(b), 1e-14)
                                        )

                                        row = {
                                            "equation": eq_key,
                                            "ndim": effective_ndim,
                                            "noise": noise,
                                            "library": lib_key,
                                            "algorithm": "stride",
                                            "space_ibp": ibp_label,
                                            "m": m_val,
                                            "d_tol": d_tol,
                                            "lam": lam,
                                            "n_active": eval_r["n_active"],
                                            "fp": eval_r["fp"],
                                            "fn": eval_r["fn"],
                                            "exact_recovery": eval_r["exact_recovery"],
                                            "E_coeff": f"{eval_r['E_coeff']:.6f}",
                                            "E_data": f"{E_data:.6f}",
                                            "cond_theta": f"{cond_val:.2e}",
                                            "n_cols": n_cols,
                                            "n_rows": n_rows,
                                            "n_snapshots": n_snapshots,
                                            "grid_points": grid_points,
                                            "build_time_ms": f"{build_time_ms:.0f}",
                                            "fit_time_ms": f"{fit_time_ms:.1f}",
                                            "peak_rss_mb": f"{mem_mb:.0f}",
                                        }
                                        writer.writerow(row)
                                        total += 1

                            elif algorithm == "rfe":
                                for rfe_lam in RFE_LAMS:
                                    t_fit = time.time()
                                    w = TrainRFE(G, b, normalize=2, lam=rfe_lam,
                                                 cache=solver_cache)
                                    fit_time_ms = (time.time() - t_fit) * 1000
                                    eval_r = evaluate_native(w, desc, true_dict)
                                    E_data = float(
                                        np.linalg.norm(G @ w - b) /
                                        max(np.linalg.norm(b), 1e-14)
                                    )

                                    row = {
                                        "equation": eq_key,
                                        "ndim": effective_ndim,
                                        "noise": noise,
                                        "library": lib_key,
                                        "algorithm": "rfe",
                                        "space_ibp": ibp_label,
                                        "m": m_val,
                                        "d_tol": "",
                                        "lam": rfe_lam,
                                        "n_active": eval_r["n_active"],
                                        "fp": eval_r["fp"],
                                        "fn": eval_r["fn"],
                                        "exact_recovery": eval_r["exact_recovery"],
                                        "E_coeff": f"{eval_r['E_coeff']:.6f}",
                                        "E_data": f"{E_data:.6f}",
                                        "cond_theta": f"{cond_val:.2e}",
                                        "n_cols": n_cols,
                                        "n_rows": n_rows,
                                        "n_snapshots": n_snapshots,
                                        "grid_points": grid_points,
                                        "build_time_ms": f"{build_time_ms:.0f}",
                                        "fit_time_ms": f"{fit_time_ms:.1f}",
                                        "peak_rss_mb": f"{mem_mb:.0f}",
                                    }
                                    writer.writerow(row)
                                    total += 1

                        f.flush()
                        if total % 100 == 0:
                            print(f"    [{total} rows written]", flush=True)

    print(f"\nDone. {total} rows -> {out_path}")


if __name__ == "__main__":
    main()
