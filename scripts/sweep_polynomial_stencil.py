#!/usr/bin/env python3
"""Small sweep of polynomial stencil choices for Allen--Cahn (clean data).

Produces a CSV summarising recovery metrics for each (deg, r) pair.

Usage:
    python scripts/sweep_polynomial_stencil.py [--quick] [--out results.csv]

Notes
-----
This script mirrors the workflow in the notebook `notebooks/allen_cahn.ipynb`:
- simulate Allen–Cahn (spectral IMEX)
- build equivariant & standard libraries using spectral (baseline) and
  polynomial backends (with varying stencil overrides)
- run STRidge and report recovery metrics (false positives/negatives and
  max coefficient error on true terms)
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from typing import List, Tuple

import numpy as np

# Ensure the sibling code directory is on sys.path (works when script is
# executed from the repository or from the `code/` folder).
THIS_DIR = os.path.dirname(__file__)
CODE_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if CODE_DIR not in sys.path:
    sys.path.insert(0, CODE_DIR)

from simulations import allen_cahn, make_regression_data
import derivatives
from library import build_equivariant_library, build_standard_library
from regression import TrainSTRidge


def feasible(deg: int, r: int) -> bool:
    """Return True if (deg, r) satisfies number-of-points >= #monomials."""
    n_pts = (2 * r + 1) ** 2
    n_mono = (deg + 1) * (deg + 2) // 2
    return n_pts >= n_mono


def run_sweep(out_csv: str, quick: bool = False) -> None:
    # Simulation params (match notebook defaults)
    eps = 0.15
    N = 64
    L = 2 * np.pi
    dt = 0.01
    if quick:
        T_end = 10.0
        save_every = 5
    else:
        T_end = 30.0
        save_every = 20

    print(f"Simulating Allen–Cahn (N={N}, L={L:.2f}, T_end={T_end})...")
    t0 = time.time()
    snapshots, times = allen_cahn(
        N=N, L=L, T_end=T_end, eps=eps, dt=dt, save_every=save_every, seed=42
    )
    print(f"  done ({len(snapshots)} snapshots, {time.time()-t0:.1f}s)")

    # Regression subsampling
    n_per_snap = 200
    f_snaps, ft_snaps, idx_sets, dt_snap = make_regression_data(
        snapshots, times, n_per_snap=n_per_snap, seed=0
    )
    Ft_all = np.concatenate([ft_snaps[s].ravel()[idx_sets[s]] for s in range(len(f_snaps))])

    # Baseline: spectral libraries & STRidge
    print("Building spectral baseline libraries...")
    Theta_eq_rows_spec, Theta_std_rows_spec = [], []
    for s in range(len(f_snaps)):
        fs = f_snaps[s]
        idx = idx_sets[s]
        with derivatives.use_backend("spectral"):
            Th_eq, desc_eq = build_equivariant_library(fs, L, L)
            Th_std, desc_std = build_standard_library(fs, L, L, max_order=4, max_poly=3)
        Theta_eq_rows_spec.append(Th_eq[idx])
        Theta_std_rows_spec.append(Th_std[idx])
    Theta_eq_spec = np.vstack(Theta_eq_rows_spec).astype(float)
    Theta_std_spec = np.vstack(Theta_std_rows_spec).astype(float)

    lam = 1e-5
    d_tol = 0.05
    print("Fitting spectral STRidge baseline...")
    w_eq_spec = TrainSTRidge(Theta_eq_spec, Ft_all, lam, d_tol)
    w_std_spec = TrainSTRidge(Theta_std_spec, Ft_all, lam, d_tol)

    # Truth vector (matches notebook): equivariant ordering assumed identical
    true_eq = np.zeros(len(desc_eq))
    true_eq[1] = 1.0
    true_eq[3] = -1.0
    true_eq[4] = eps ** 2
    act = np.abs(true_eq) > 1e-10

    def metrics(w: np.ndarray) -> Tuple[int, int, int, float]:
        nz = int(np.count_nonzero(np.abs(w) > 1e-10))
        fp = int(np.sum((np.abs(w) > 1e-10) & ~act))
        fn = int(np.sum((np.abs(w) < 1e-10) & act))
        max_err = float(np.max(np.abs(w[act] - true_eq[act])))
        return nz, fp, fn, max_err

    base_metrics_eq = metrics(w_eq_spec)
    nz_eq_s, fp_eq_s, fn_eq_s, err_eq_s = base_metrics_eq
    nz_std_s = int(np.count_nonzero(np.abs(w_std_spec) > 1e-10))
    cond_eq_s = float(np.linalg.cond(Theta_eq_spec))
    cond_std_s = float(np.linalg.cond(Theta_std_spec))
    print(
        f"Spectral baseline: eq_nz={nz_eq_s}, fp={fp_eq_s}, fn={fn_eq_s}, "
        f"err={err_eq_s:.3e}, cond_eq={cond_eq_s:.2e}, std_nz={nz_std_s}"
    )

    # Sweep candidates
    deg_candidates = [3, 5, 7, 9, 11]
    r_candidates = [2, 3, 4, 5, 6]

    results: List[dict] = []

    for deg in deg_candidates:
        for r in r_candidates:
            if not feasible(deg, r):
                print(f"Skipping infeasible deg={deg}, r={r}")
                continue

            print(f"\nTesting deg={deg}, r={r}...")
            # Clear any previous overrides (global + per-order), then set
            # per-order overrides only for orders that the chosen `deg` can
            # represent (i.e. deg >= order). This avoids creating invalid
            # small-degree overrides for high-order operators (which would
            # cause _poly_kernel to fail when a+b > deg).
            derivatives.set_polynomial_stencil(None, None)
            for ord_ in (1, 2, 3, 4, 6):
                if deg >= ord_:
                    derivatives.set_polynomial_stencil(deg=deg, r=r, order=ord_)
                else:
                    derivatives.set_polynomial_stencil(None, None, order=ord_)

            Theta_eq_rows_poly, Theta_std_rows_poly = [], []
            for s in range(len(f_snaps)):
                fs = f_snaps[s]
                idx = idx_sets[s]
                with derivatives.use_backend("polynomial"):
                    Th_eq, _ = build_equivariant_library(fs, L, L)
                    Th_std, _ = build_standard_library(fs, L, L, max_order=4, max_poly=3)
                Theta_eq_rows_poly.append(Th_eq[idx])
                Theta_std_rows_poly.append(Th_std[idx])

            Theta_eq_poly = np.vstack(Theta_eq_rows_poly).astype(float)
            Theta_std_poly = np.vstack(Theta_std_rows_poly).astype(float)

            # Fit
            w_eq_poly = TrainSTRidge(Theta_eq_poly, Ft_all, lam, d_tol)
            w_std_poly = TrainSTRidge(Theta_std_poly, Ft_all, lam, d_tol)

            # compute metrics
            nz_eq, fp_eq, fn_eq, err_eq = metrics(w_eq_poly)
            nz_std = int(np.count_nonzero(np.abs(w_std_poly) > 1e-10))

            # library condition numbers
            cond_eq = float(np.linalg.cond(Theta_eq_poly))
            cond_std = float(np.linalg.cond(Theta_std_poly))

            print(
                f"deg={deg} r={r}: eq_nz={nz_eq}, fp={fp_eq}, fn={fn_eq}, err={err_eq:.3e}, cond_eq={cond_eq:.2e}, std_nz={nz_std}"
            )

            results.append(
                {
                    "deg": deg,
                    "r": r,
                    "cond_eq": cond_eq,
                    "cond_std": cond_std,
                    "eq_nz": nz_eq,
                    "eq_fp": fp_eq,
                    "eq_fn": fn_eq,
                    "eq_max_coeff_err": err_eq,
                    "std_nz": nz_std,
                }
            )

    # Write CSV — first row is the spectral baseline for easy comparison
    keys = [
        "deg",
        "r",
        "cond_eq",
        "cond_std",
        "eq_nz",
        "eq_fp",
        "eq_fn",
        "eq_max_coeff_err",
        "std_nz",
    ]
    spectral_row = {
        "deg": "spectral",
        "r": "spectral",
        "cond_eq": cond_eq_s,
        "cond_std": cond_std_s,
        "eq_nz": nz_eq_s,
        "eq_fp": fp_eq_s,
        "eq_fn": fn_eq_s,
        "eq_max_coeff_err": err_eq_s,
        "std_nz": nz_std_s,
    }
    print(f"\nWriting results to {out_csv}")
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerow(spectral_row)
        for row in results:
            writer.writerow(row)

    print("Sweep complete.")


def main(argv: List[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Sweep polynomial stencils (Allen–Cahn, clean data)")
    parser.add_argument("--quick", action="store_true", help="Run a shorter/quick simulation")
    parser.add_argument("--out", default="results/sweep_results.csv", help="Output CSV path")
    args = parser.parse_args(argv)
    run_sweep(out_csv=args.out, quick=args.quick)


if __name__ == "__main__":
    main()
