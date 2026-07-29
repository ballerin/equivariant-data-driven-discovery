#!/usr/bin/env python3
"""Benchmark combo-level parallelism for STRidge hyperparameter sweeps.

Builds the Theta/Ft matrix once (Allen-Cahn, spectral, multigraph, pointwise),
then runs the 68 (lam × d_tol) STRidge combos with different worker counts
to find the optimal parallelization strategy.

Compares:
  - Sequential (1 worker)
  - ThreadPoolExecutor with 2, 4, 8 workers
  - All with OMP_NUM_THREADS=1 (single-threaded BLAS)

Usage:
  python scripts/benchmark_combo_parallel.py
"""

import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# Force single-threaded BLAS so we isolate combo-level parallelism
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives
from equivariant_pde_find.library import build_multigraph_library, build_polyjet_library
from equivariant_pde_find.library.conversions import _MULTIGRAPH_TO_JET_DECOMPOSITION
from equivariant_pde_find.evaluation import build_eval_theta_jet, evaluate_unified
from equivariant_pde_find.equations import true_allen_cahn
from equivariant_pde_find.regression import TrainSTRidge
from equivariant_pde_find.simulations import make_regression_data

DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))

# ── hyperparameter grid (same as sweep.py STRidge defaults) ─────────────
LAMS = [1e-6, 1e-4, 1e-2, 1.0]
D_TOLS = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, 0.5,
          1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0]
N_RUNS = 3  # repeat each worker count for stability
WORKER_COUNTS = [1, 2, 4, 8]


def load_and_build():
    """Load AC clean data and build Theta, Ft, ground-truth, eval matrices."""
    npz_path = os.path.join(DATA_DIR, "allen_cahn_clean.npz")
    data = np.load(npz_path)
    snapshots = data["snapshots"]
    times = data["times"]
    Lx, Ly = float(data["Lx"]), float(data["Ly"])

    # Regression data
    f_snaps, ft_snaps, idx_sets, _ = make_regression_data(
        snapshots, times, n_per_snap=200, skip_fraction=0.0, seed=0,
    )

    # Build library Theta (once)
    with derivatives.use_backend("spectral"):
        rows_theta, rows_ft = [], []
        for s in range(len(f_snaps)):
            T, desc = build_multigraph_library(f_snaps[s], Lx, Ly)
            rows_theta.append(T[idx_sets[s]])
            rows_ft.append(ft_snaps[s].ravel()[idx_sets[s]])
    Theta = np.vstack(rows_theta).astype(float)
    Ft = np.concatenate(rows_ft).astype(float)

    # Ground truth
    true_eq = true_allen_cahn(desc, eps=0.15)

    # Jet descriptions for evaluation
    _, desc_jet = build_polyjet_library(
        snapshots[0], Lx, Ly, max_order=4, max_degree=3,
    )
    true_jet = np.zeros(len(desc_jet))
    for j, (name, coeff) in enumerate(zip(desc, true_eq)):
        if abs(coeff) < 1e-12:
            continue
        if name not in _MULTIGRAPH_TO_JET_DECOMPOSITION:
            continue
        for jet_name, mult in _MULTIGRAPH_TO_JET_DECOMPOSITION[name]:
            true_jet[desc_jet.index(jet_name)] += coeff * mult

    # Evaluation Theta (for unified metric)
    Theta_eval, Ft_eval = build_eval_theta_jet(
        snapshots, times, Lx, Ly, skip_fraction=0.0, seed=0, backend="spectral",
    )

    return Theta, Ft, desc, true_eq, desc_jet, true_jet, Theta_eval, Ft_eval


def run_one_combo(args):
    """Run STRidge for one (lam, d_tol) pair.  Returns row dict."""
    lam, d_tol, Theta, Ft, desc, desc_jet, Theta_eval, Ft_eval, true_jet = args
    w = TrainSTRidge(Theta, Ft, lam, d_tol)
    eval_result = evaluate_unified(
        w, "multigraph", desc, desc_jet,
        Theta_eval, Ft_eval, true_jet,
        M_ws_to_jet=None,
    )
    n_active = int(np.count_nonzero(np.abs(w) > 1e-10))
    return {
        "lam": lam, "d_tol": d_tol,
        "n_active": n_active,
        "exact": eval_result["exact_recovery"],
        "E_coeff": eval_result["E_coeff"],
    }


def run_benchmark(Theta, Ft, desc, desc_jet, true_jet, Theta_eval, Ft_eval):
    """Time the full STRidge grid with different worker counts."""
    # Build task list (68 combos)
    tasks = [
        (lam, d_tol, Theta, Ft, desc, desc_jet, Theta_eval, Ft_eval, true_jet)
        for lam in LAMS for d_tol in D_TOLS
    ]
    n_tasks = len(tasks)
    print(f"Grid: {len(LAMS)} λ × {len(D_TOLS)} d_tol = {n_tasks} combos")
    print(f"Theta shape: {Theta.shape}, cond={np.linalg.cond(Theta):.2e}")
    print()

    results = {}
    for n_workers in WORKER_COUNTS:
        times = []
        exact_counts = []
        for run_idx in range(N_RUNS):
            t0 = time.perf_counter()
            if n_workers == 1:
                # Sequential
                rows = [run_one_combo(t) for t in tasks]
            else:
                # Parallel with ThreadPoolExecutor
                with ThreadPoolExecutor(max_workers=n_workers) as ex:
                    rows = list(ex.map(run_one_combo, tasks))
            elapsed = time.perf_counter() - t0
            times.append(elapsed)
            exact_counts.append(sum(1 for r in rows if r["exact"]))

        avg_time = np.mean(times)
        std_time = np.std(times)
        speedup = results[1]["avg_time"] / avg_time if 1 in results else 1.0

        results[n_workers] = {
            "avg_time": avg_time,
            "std_time": std_time,
            "speedup": speedup,
            "times": times,
            "exact_counts": exact_counts,
        }

        print(
            f"  workers={n_workers:>2}  "
            f"avg={avg_time:.3f}s  ±{std_time:.3f}s  "
            f"speedup={speedup:.2f}×  "
            f"exact={exact_counts}"
        )

    return results


def main():
    print("=" * 65)
    print("  Combo-level parallelism benchmark (OMP_NUM_THREADS=1)")
    print("  Allen–Cahn, spectral pointwise, STRidge")
    print("=" * 65)
    print()

    t_load = time.perf_counter()
    Theta, Ft, desc, true_eq, desc_jet, true_jet, Theta_eval, Ft_eval = load_and_build()
    print(f"Data + Theta build: {time.perf_counter() - t_load:.1f}s")
    print()

    results = run_benchmark(Theta, Ft, desc, desc_jet, true_jet, Theta_eval, Ft_eval)

    # ── Summary ─────────────────────────────────────────────────────────
    print()
    print("=" * 65)
    print("  Summary")
    print("=" * 65)
    print(f"  {'Workers':<10} {'Avg (s)':<10} {'Speedup':<10} {'Best for'}")
    print(f"  {'-'*10} {'-'*10} {'-'*10} {'-'*30}")
    for n in WORKER_COUNTS:
        r = results[n]
        note = "baseline" if n == 1 else ""
        if n > 1 and r["speedup"] < 1.0:
            note = "⚠ slower than sequential"
        elif n > 1 and r["speedup"] > results.get(n - 1, {}).get("speedup", 0):
            note = "✓ current best"
        print(f"  {n:<10} {r['avg_time']:<10.3f} {r['speedup']:<10.2f} {note}")

    best = min(WORKER_COUNTS, key=lambda n: results[n]["avg_time"])
    print()
    print(f"  ➤ Optimal: {best} workers ({results[best]['avg_time']:.3f}s)")
    print(f"  ➤ Speedup vs sequential: {results[best]['speedup']:.2f}×")


if __name__ == "__main__":
    main()
