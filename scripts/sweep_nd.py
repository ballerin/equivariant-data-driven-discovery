#!/usr/bin/env python
"""n-dimensional PDE discovery sweep: AC & CH with multigraph / mod_multigraph / pdefind.

Compares STRidge and RFE across libraries, equations, dimensions, and noise levels.
Uses weak-form IBP with the spectral derivative backend.

Usage
-----
    python scripts/sweep_nd.py
    python scripts/sweep_nd.py --equations ac ch --ndim 2 3
    python scripts/sweep_nd.py --libraries multigraph pdefind --algorithms stride
    python scripts/sweep_nd.py --noise 0.0 0.02
"""

import argparse
import csv
import datetime
import os
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
)
from equivariant_pde_find.library import weak_form as wf
from equivariant_pde_find.regression import TrainSTRidge, TrainRFE
from equivariant_pde_find.simulations import make_regression_data
from equivariant_pde_find.equations import (
    CH_EPSILON,
    true_allen_cahn,
    true_cahn_hilliard,
    true_allen_cahn_mod_multigraph,
    true_cahn_hilliard_mod_multigraph,
    true_allen_cahn_pdefind,
    true_cahn_hilliard_pdefind,
)

DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyper-parameters ────────────────────────────────────────────────────────
from equivariant_pde_find.library.weak_form import (  # noqa: E402
    select_test_function_degree as _wf_select_degree,
    test_function_is_resolved as _wf_is_resolved,
)

M_VALS = [12]          # weak-form half-support
# Test-function degree follows Messenger & Bortz (2021), Algorithm 1, rather
# than a fixed constant: with a fixed p the fourth-order weak-IBP columns pick
# up a 3 % (m=12) to 84 % (m=8) systematic error.
# The direct weak form never differentiates psi, so it keeps the old degree.
WEAK_MAX_DERIV = 4
WEAK_P_DIRECT = 6


def weak_p(m, use_ibp=True):
    """Test-function degree for support *m* (see WEAK_MAX_DERIV above)."""
    if not use_ibp:
        return WEAK_P_DIRECT
    return _wf_select_degree(m, WEAK_MAX_DERIV)


def weak_ibp_supports(m_values):
    """Subset of *m_values* whose grid can represent psi^(WEAK_MAX_DERIV)."""
    return [m for m in m_values
            if _wf_is_resolved(m, weak_p(m), WEAK_MAX_DERIV)]
N_PER_SNAP = 200       # spatial samples per snapshot
SAMPLE_STRIDE = 4
LAM = 1e-6             # STRidge ridge penalty (legacy)
D_TOLS = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, 0.5,
          1.0, 5.0, 15.0, 20.0]
RFE_LAMS = [0.0, 1e-6, 1e-4, 1e-2, 1.0]

# ── library builders ─────────────────────────────────────────────────────────

def _build_weak_nd(f_snaps, ft_snaps, idx_sets, Lx, Ly, Lz,
                   library_fn, m, p, use_ibp):
    """Build weak-form (G, b) for n-d data with a given library builder.

    Handles both 2-D (Lz=None) and 3-D (Lz=float) automatically.
    """
    kernel = wf._cached_kernel_1d(m, p)
    rows_G, rows_b = [], []

    for s in range(len(f_snaps)):
        f = f_snaps[s]
        ft = ft_snaps[s]
        ndim = f.ndim

        if library_fn is build_pdefind_library:
            G, desc = library_fn(
                f, Lx, Ly, max_order=4, max_poly=3,
                formulation="weak", m=m, p=p, use_ibp=use_ibp, Lz=Lz,
            )
        else:
            G, desc = library_fn(
                f, Lx, Ly, max_p=2, max_n=3,
                formulation="weak", m=m, p=p, use_ibp=use_ibp, Lz=Lz,
            )

        # Build weak RHS
        if ndim == 3:
            b_w = wf.build_weak_rhs(ft, kernel_x=kernel, kernel_y=kernel,
                                    kernel_z=kernel, m=m, p=p)
        else:
            b_w = wf.build_weak_rhs(ft, kernel_x=kernel, kernel_y=kernel,
                                    m=m, p=p)

        # Strided sampling of valid interior points (C-order / row-major).
        shape = f.shape
        Nv = tuple(s - 2 * m for s in shape)
        grid_axes = [np.arange(0, Nv[i], SAMPLE_STRIDE) for i in range(ndim)]
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
        "--equations", nargs="+", default=["ac", "ch"],
        help="Equations to sweep (ac, ch).",
    )
    parser.add_argument(
        "--ndim", nargs="+", type=int, default=[2, 3],
        help="Spatial dimensions (default: 2 3).",
    )
    parser.add_argument(
        "--noise", nargs="+", type=float, default=[0.0, 0.02, 0.05],
        help="Noise levels (default: 0.0 0.02 0.05).",
    )
    parser.add_argument(
        "--libraries", nargs="+",
        default=["multigraph", "mod_multigraph", "pdefind"],
        help="Library types to test.",
    )
    parser.add_argument(
        "--algorithms", nargs="+", default=["stride", "rfe"],
        help="Regression algorithms (stride, rfe).",
    )
    parser.add_argument(
        "--name", default="", help="Tag appended to output filename.",
    )
    args = parser.parse_args()

    derivatives.set_backend("spectral")

    # ── library + truth mapping ──────────────────────────────────────────
    LIB_MAP = {
        "multigraph": {
            "builder": build_multigraph_library,
            "truth_ac": true_allen_cahn,
            "truth_ch": true_cahn_hilliard,
            "truth_fn_ndim": False,  # truth is dimension-agnostic
        },
        "mod_multigraph": {
            "builder": build_mod_multigraph_library,
            "truth_ac": true_allen_cahn_mod_multigraph,
            "truth_ch": true_cahn_hilliard_mod_multigraph,
            "truth_fn_ndim": False,
        },
        "pdefind": {
            "builder": build_pdefind_library,
            "truth_ac": true_allen_cahn_pdefind,
            "truth_ch": true_cahn_hilliard_pdefind,
            "truth_fn_ndim": True,  # truth depends on ndim
        },
    }

    # ── discover datasets ────────────────────────────────────────────────
    noise_labels = {0.0: "clean", 0.02: "noise2pct", 0.05: "noise5pct"}
    eq_names = {"ac": "allen_cahn", "ch": "cahn_hilliard"}

    datasets = []
    for eq_key in args.equations:
        eq_name = eq_names[eq_key]
        for ndim in args.ndim:
            for noise in args.noise:
                label = noise_labels[noise]
                fname = f"{eq_name}_nd{ndim}_{label}.npz"
                # Also try legacy 2D names without _nd suffix
                if ndim == 2 and not os.path.exists(
                    os.path.join(DATA_DIR, fname)
                ):
                    legacy = f"{eq_name}_{label}.npz"
                    if os.path.exists(os.path.join(DATA_DIR, legacy)):
                        fname = legacy
                path = os.path.join(DATA_DIR, fname)
                if os.path.exists(path):
                    datasets.append((eq_key, ndim, noise, fname, path))
                else:
                    print(f"SKIP: {fname} not found. Run generate_nd_data.py first.")

    if not datasets:
        print("No datasets found.")
        return

    # ── output ───────────────────────────────────────────────────────────
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"_{args.name}" if args.name else ""
    out_path = os.path.join(RESULTS_DIR, f"sweep_nd_{timestamp}{tag}.csv")

    fieldnames = [
        "equation", "ndim", "noise", "library", "algorithm",
        "m", "d_tol", "lam",
        "n_active", "fp", "fn", "exact_recovery", "E_coeff",
        "E_data", "cond_theta", "time_ms",
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

            # Prepare regression points
            f_snaps, ft_snaps, idx_sets, _ = make_regression_data(
                snapshots, times, n_per_snap=N_PER_SNAP,
                skip_fraction=0.0, seed=0,
            )
            print(f"  snaps={len(f_snaps)}, grid={snapshots.shape[1:]}")

            for lib_key in args.libraries:
                lib_info = LIB_MAP[lib_key]
                builder = lib_info["builder"]

                # Build ground truth
                if lib_info["truth_fn_ndim"]:
                    if eq_key == "ac":
                        truth_fn = lib_info["truth_ac"]
                        true_dict = truth_fn(effective_ndim, eps=0.15)
                    else:
                        truth_fn = lib_info["truth_ch"]
                        true_dict = truth_fn(effective_ndim, eps=CH_EPSILON)
                else:
                    if eq_key == "ac":
                        truth_fn = lib_info["truth_ac"]
                    else:
                        truth_fn = lib_info["truth_ch"]
                    # Build on one snapshot to get descriptions, then build truth
                    if builder is build_pdefind_library:
                        _, desc_tmp = builder(
                            f_snaps[0], Lx, Ly, max_order=4, max_poly=3,
                            formulation="weak", m=M_VALS[0], p=weak_p(M_VALS[0]),
                            use_ibp=True, Lz=Lz,
                        )
                    else:
                        _, desc_tmp = builder(
                            f_snaps[0], Lx, Ly, max_p=2, max_n=3,
                            formulation="weak", m=M_VALS[0], p=weak_p(M_VALS[0]),
                            use_ibp=True, Lz=Lz,
                        )
                    true_w = truth_fn(desc_tmp)
                    true_dict = {d: true_w[i] for i, d in enumerate(desc_tmp)
                                 if abs(true_w[i]) > 1e-12}

                for m_val in M_VALS:
                    # Build weak-form library once per m
                    t_build = time.time()
                    G, b, desc = _build_weak_nd(
                        f_snaps, ft_snaps, idx_sets, Lx, Ly, Lz,
                        builder, m_val, weak_p(m_val), use_ibp=True,
                    )
                    cond_val = cond_number(G)
                    build_ms = (time.time() - t_build) * 1000
                    print(f"  [{lib_key:>16s}] m={m_val}  "
                          f"shape={G.shape}  cond={cond_val:.2e}  "
                          f"build={build_ms:.0f}ms")

                    for algorithm in args.algorithms:
                        if algorithm == "stride":
                            for d_tol in D_TOLS:
                                t_row = time.time()
                                w = TrainSTRidge(G, b, LAM, d_tol,
                                                 maxit=25, seed=0)
                                eval_r = evaluate_native(w, desc, true_dict)
                                E_data = float(
                                    np.linalg.norm(G @ w - b) /
                                    max(np.linalg.norm(b), 1e-14)
                                )
                                dt_ms = (time.time() - t_row) * 1000

                                row = {
                                    "equation": eq_key,
                                    "ndim": effective_ndim,
                                    "noise": noise,
                                    "library": lib_key,
                                    "algorithm": "stride",
                                    "m": m_val,
                                    "d_tol": d_tol,
                                    "lam": LAM,
                                    "n_active": eval_r["n_active"],
                                    "fp": eval_r["fp"],
                                    "fn": eval_r["fn"],
                                    "exact_recovery": eval_r["exact_recovery"],
                                    "E_coeff": f"{eval_r['E_coeff']:.6f}",
                                    "E_data": f"{E_data:.6f}",
                                    "cond_theta": f"{cond_val:.2e}",
                                    "time_ms": f"{dt_ms:.1f}",
                                }
                                writer.writerow(row)
                                total += 1

                        elif algorithm == "rfe":
                            for rfe_lam in RFE_LAMS:
                                t_row = time.time()
                                w = TrainRFE(G, b, normalize=2, lam=rfe_lam)
                                eval_r = evaluate_native(w, desc, true_dict)
                                E_data = float(
                                    np.linalg.norm(G @ w - b) /
                                    max(np.linalg.norm(b), 1e-14)
                                )
                                dt_ms = (time.time() - t_row) * 1000

                                row = {
                                    "equation": eq_key,
                                    "ndim": effective_ndim,
                                    "noise": noise,
                                    "library": lib_key,
                                    "algorithm": "rfe",
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
                                    "time_ms": f"{dt_ms:.1f}",
                                }
                                writer.writerow(row)
                                total += 1

                    f.flush()
                    if total % 50 == 0:
                        print(f"    [{total} rows written]", flush=True)

    print(f"\nDone. {total} rows -> {out_path}")


if __name__ == "__main__":
    main()
