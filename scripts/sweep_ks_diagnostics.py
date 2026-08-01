"""
Diagnostic sweep: KS weak-chaotic — multigraph only, polynomial backend, STRidge.

Compares multigraph weak-form (direct + IBP) using only the polynomial stencil
backend.  Sweeps over a coarse d_tol grid to identify promising
data-hyperparameter combinations.

Usage
-----
    python scripts/sweep_ks_diagnostics.py
    python scripts/sweep_ks_diagnostics.py --data-dirs data/ data/ks_diagnostics/
"""

import argparse
import csv
import datetime
import gc
import itertools
import os
import sys
import time

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives  # noqa: E402
from equivariant_pde_find.evaluation import (  # noqa: E402
    evaluate_unified,
)
from equivariant_pde_find.library.conversions import (  # noqa: E402
    map_multigraph_to_jet, map_mod_multigraph_to_jet,
    collect_needed_jet_monomials,
)
from equivariant_pde_find.equations import (  # noqa: E402
    true_ks, normalize_eq_name,
)
from equivariant_pde_find.library import (  # noqa: E402
    build_multigraph_library,
    build_mod_multigraph_library,
)
from equivariant_pde_find.regression import TrainSTRidge  # noqa: E402
from equivariant_pde_find.simulations import make_regression_data  # noqa: E402
from equivariant_pde_find.library.weak_form import build_weak_rhs  # noqa: E402

RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyperparameter grids ────────────────────────────────────────────────────
LAM = 1e-6
D_TOLS = [1e-3, 1e-2, 0.1, 1.0, 10.0]           # coarse grid
BACKENDS = ["polynomial"]                           # polynomial only
LIBRARIES = ["multigraph", "mod_multigraph"]
from equivariant_pde_find.library.weak_form import (  # noqa: E402
    select_test_function_degree as _wf_select_degree,
    test_function_is_resolved as _wf_is_resolved,
)

WEAK_M = [4, 8, 12]
# Test-function degree follows Messenger & Bortz (2021), Algorithm 1, rather
# than a fixed constant: with a fixed p the fourth-order weak-IBP columns pick
# up a 3 % (m=12) to 84 % (m=8) systematic error.  See audit/FINDINGS.md F-01.
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
SAMPLE_STRIDE = 4

# ── helpers ─────────────────────────────────────────────────────────────────


def cond_number(X):
    s = np.linalg.svd(X, compute_uv=False)
    return float(s[0] / (s[-1] + 1e-30))


def load_dataset(npz_path):
    data = np.load(npz_path, allow_pickle=True)
    snapshots = data["snapshots"]
    times = data["times"]
    Lx = float(data["Lx"])
    Ly = float(data["Ly"])
    N = int(data["N"])
    obs_noise = float(data.get("obs_noise", 0.0))
    eq_name = str(data.get("eq_name", ""))
    if not eq_name and "ks" in os.path.basename(npz_path).lower():
        # Infer from filename for legacy datasets
        fname = os.path.basename(npz_path)
        if "cellular" in fname or "ks_cell" in fname:
            eq_name = "ks_cell"
        elif "weak_chaotic" in fname or "ks_wchaos" in fname:
            eq_name = "ks_wchaos"
        elif "chaotic" in fname or "ks_chaos" in fname:
            eq_name = "ks_chaos"
    sim_eps = data.get("sim_eps", None)
    if sim_eps is not None:
        sim_eps = float(sim_eps)
    # Diagnostic metadata
    sim_N = data.get("sim_N", N)
    sim_L = data.get("sim_L", Lx)
    sim_T_end = str(data.get("sim_T_end", ""))
    sim_dt_snap = str(data.get("sim_dt_snap", ""))
    return (snapshots, times, Lx, Ly, N, obs_noise, eq_name, sim_eps,
            sim_N, sim_L, sim_T_end, sim_dt_snap)


def build_weak_form(f_snaps, ft_snaps, idx_sets, Lx, Ly, m, use_ibp):
    """Build weak-form (Theta, Ft) for multigraph library, polynomial backend."""
    with derivatives.use_backend("polynomial"):
        rows_theta, rows_ft = [], []
        for s in range(len(f_snaps)):
            f = f_snaps[s]
            ft = ft_snaps[s]
            G, desc = build_multigraph_library(
                f, Lx, Ly, max_p=2, max_n=3,
                formulation="weak", m=m, p=weak_p(m, use_ibp), use_ibp=use_ibp,
            )
            b_w = build_weak_rhs(ft, m=m, p=weak_p(m, use_ibp))

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

            rows_theta.append(G[q_idx])
            rows_ft.append(b_w[q_idx])

    Theta = np.vstack(rows_theta).astype(float)
    Ft = np.concatenate(rows_ft).astype(float)
    return Theta, Ft, desc


def _map_eq_key(eq_name):
    """Map dataset eq_name to equations key."""
    eq_name = normalize_eq_name(eq_name)
    mapping = {
        "ks_cell": "ks_cell", "ks_wchaos": "ks_wchaos", "ks_chaos": "ks_chaos",
        "ks_cellular": "ks_cell", "ks_weak_chaotic": "ks_wchaos",
        "ks_chaotic": "ks_chaos",
    }
    return mapping.get(eq_name, eq_name)


# ── main ────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Diagnostic KS sweep: multigraph + polynomial + STRidge."
    )
    parser.add_argument(
        "--data-dirs", nargs="+",
        default=["data/", "data/ks_diagnostics/"],
        help="Directories to scan for .npz files.",
    )
    parser.add_argument(
        "--name", default="",
        help="Optional tag for output filename.",
    )
    args = parser.parse_args()

    # Discover .npz files
    all_npz = []
    for d in args.data_dirs:
        d = os.path.abspath(os.path.join(SCRIPT_DIR, "..", d))
        if os.path.isdir(d):
            for f in sorted(os.listdir(d)):
                if f.endswith(".npz") and not f.startswith("."):
                    all_npz.append((d, f))

    print(f"Found {len(all_npz)} .npz files across {len(args.data_dirs)} dirs.",
          flush=True)

    # Output
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"_{args.name}" if args.name else ""
    out_path = os.path.join(RESULTS_DIR, f"ks_diag_{timestamp}{tag}.csv")

    fieldnames = [
        "equation", "noise", "method", "library", "weak_m", "lam", "d_tol",
        "n_active", "fp", "fn", "exact_recovery",
        "E_coeff", "E_data", "cond_theta",
        "dataset_file", "sim_N", "sim_L", "sim_T_end", "sim_dt_snap",
        "time_ms",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        f.flush()

        for data_dir, npz_file in all_npz:
            npz_path = os.path.join(data_dir, npz_file)
            (snapshots, times, Lx, Ly, N, obs_noise, eq_name, sim_eps,
             sim_N, sim_L, sim_T_end, sim_dt_snap) = load_dataset(npz_path)
            eq_name = normalize_eq_name(eq_name)
            true_fn_key = _map_eq_key(eq_name)

            # Only KS datasets (by filename pattern or stored eq_name)
            is_ks = (npz_file.startswith("ks_") or npz_file.startswith("ks_diag_"))
            if not is_ks:
                print(f"\n{'='*70}\nDataset: {npz_file}\n{'='*70}")
                print(f"  SKIP (not KS)", flush=True)
                continue

            # Normalize eq_name from file or metadata
            if true_fn_key in ("ks_wchaos", "ks_cell", "ks_chaos",
                               "ks_cellular", "ks_weak_chaotic", "ks_chaotic"):
                pass  # already valid
            elif "ks_diag" in npz_file or "ks_wchaos" in str(eq_name):
                true_fn_key = "ks_wchaos"
            elif "ks_cell" in npz_file or "ks_cellular" in npz_file:
                true_fn_key = "ks_cell"
            elif "ks_chaos" in npz_file or "ks_chaotic" in npz_file:
                true_fn_key = "ks_chaos"
            else:
                print(f"  SKIP (unknown KS variant: {true_fn_key})", flush=True)
                continue

            noise_label = f"{obs_noise:.4f}"

            print(f"\n{'='*70}\nDataset: {npz_file}\n{'='*70}")
            print(f"  eq={true_fn_key}  noise={noise_label}"
                  f"  N={sim_N}  L={sim_L}  T_end={sim_T_end}"
                  f"  dt_snap={sim_dt_snap}", flush=True)

            # Prepare regression data
            f_snaps_raw, ft_snaps_raw, idx_sets_raw, _ = make_regression_data(
                snapshots, times, n_per_snap=200,
                skip_fraction=0.3, seed=0,  # KS: skip 30% transient
            )
            f_snaps = f_snaps_raw
            ft_snaps = ft_snaps_raw
            idx_sets = idx_sets_raw

            # Build eval Theta on subsampled points (avoids OOM)
            rows_theta, rows_ft = [], []
            for s in range(len(f_snaps_raw)):
                T_full, _ = build_polyjet_library(
                    f_snaps_raw[s], Lx, Ly, max_order=4, max_degree=3,
                )
                rows_theta.append(T_full[idx_sets_raw[s]])
                rows_ft.append(ft_snaps_raw[s].ravel()[idx_sets_raw[s]])
            Tj = np.vstack(rows_theta)
            Ftj = np.concatenate(rows_ft)
            eval_theta = (Tj, Ftj)
            # Dynamic jet: union of multigraph bases
            _, desc_mg_jet = build_multigraph_library(
                f_snaps[0], Lx, Ly, max_p=2, max_n=3,
            )
            _, desc_mod_jet = build_mod_multigraph_library(
                f_snaps[0], Lx, Ly, max_p=2, max_n=3,
            )
            desc_jet = collect_needed_jet_monomials(
                desc_mg=desc_mg_jet, desc_mod_mg=desc_mod_jet,
            )

            # Build ground truth in jet basis for each library
            true_jet = {}
            for lib in LIBRARIES:
                if lib == "multigraph":
                    _, desc_lib = build_multigraph_library(
                        f_snaps[0], Lx, Ly, max_p=2, max_n=3,
                        formulation="weak", m=4, p=weak_p(4, False), use_ibp=False,
                    )
                else:
                    _, desc_lib = build_mod_multigraph_library(
                        f_snaps[0], Lx, Ly, max_p=2, max_n=3,
                        formulation="weak", m=4, p=weak_p(4, False), use_ibp=False,
                    )
                true_w = true_ks(desc_lib)
                if lib == "multigraph":
                    true_jet[lib] = map_multigraph_to_jet(true_w, desc_lib, desc_jet)
                else:
                    true_jet[lib] = map_mod_multigraph_to_jet(true_w, desc_lib, desc_jet)

            # Build one (m, ibp) at a time, run STRidge, free
            for m_val, use_ibp in itertools.product(WEAK_M, [False, True]):
                # A support too small to represent psi^(WEAK_MAX_DERIV) cannot
                # carry the IBP formulation at any degree; the direct weak form
                # never differentiates psi and is fine.
                if use_ibp and m_val not in weak_ibp_supports([m_val]):
                    print(f"    Skipping: m={m_val} ibp=True (2m+1 points "
                          f"cannot resolve a degree-{weak_p(m_val)} test "
                          f"function; see audit/FINDINGS.md F-01)", flush=True)
                    continue

                method_key = "weak_ibp" if use_ibp else "weak_direct"

                t0 = time.time()
                print(f"    Building: m={m_val} ibp={use_ibp}", flush=True)
                Theta, Ft, desc = build_weak_form(
                    f_snaps, ft_snaps, idx_sets, Lx, Ly, m_val, use_ibp,
                )
                cond_val = cond_number(Theta)
                print(f"      shape={Theta.shape} cond={cond_val:.2e}"
                      f" ({time.time()-t0:.1f}s)", flush=True)

                for d_tol in D_TOLS:
                    t_row = time.time()
                    w = TrainSTRidge(Theta, Ft, LAM, d_tol, maxit=25, seed=0)
                    eval_result = evaluate_unified(
                        w, lib, desc, desc_jet,
                        *eval_theta, true_jet[lib],
                    )
                    dt_ms = (time.time() - t_row) * 1000

                    row = dict(
                        equation=true_fn_key,
                        noise=noise_label,
                        method=method_key,
                        weak_m=m_val,
                        lam=LAM,
                        d_tol=d_tol,
                        n_active=eval_result["n_active"],
                        fp=eval_result["fp"],
                        fn=eval_result["fn"],
                        exact_recovery=eval_result["exact_recovery"],
                        E_coeff=eval_result["E_coeff"],
                        E_data=eval_result.get("E_data", ""),
                        cond_theta=cond_val,
                        dataset_file=npz_file,
                        sim_N=sim_N,
                        sim_L=sim_L,
                        sim_T_end=sim_T_end,
                        sim_dt_snap=sim_dt_snap,
                        time_ms=f"{dt_ms:.1f}",
                    )
                    writer.writerow(row)
                    f.flush()

                del Theta, Ft, desc

            # Free per-dataset memory
            del eval_theta, true_jet, desc_jet, desc_lib
            del f_snaps_raw, ft_snaps_raw, idx_sets_raw, f_snaps, ft_snaps, idx_sets, snapshots
            gc.collect()

    print(f"\nDone. Results → {out_path}", flush=True)


if __name__ == "__main__":
    main()
