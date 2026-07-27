#!/usr/bin/env python
"""Focused 3-D Allen-Cahn experiment: test 3 configurations with RFE + weak IBP.

Configurations tested:
  1. multigraph + RFE, weak IBP, m=12, spectral backend
  2. mod_multigraph + RFE, weak IBP, m=12, spectral backend
  3. pdefind + RFE, weak IBP, m=12, spectral backend

Each tested on clean, 2% noise, and 5% noise data.
"""

import os
import sys
import csv
import datetime
import time

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives
from equivariant_pde_find.library import (
    build_multigraph_library, build_mod_multigraph_library, build_pdefind_library,
)
from equivariant_pde_find.library import weak_form as wf
from equivariant_pde_find.regression import TrainRFE
from equivariant_pde_find.simulations import make_regression_data

DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
os.makedirs(RESULTS_DIR, exist_ok=True)

# Experiment parameters
M = 12       # weak form half-support
WEAK_P = 6   # test function polynomial degree
N_PER_SNAP = 200  # spatial samples per snapshot
SAMPLE_STRIDE = 4


def load_3d_dataset(npz_path):
    data = np.load(npz_path)
    return (
        data["snapshots"], data["times"],
        float(data["Lx"]), float(data["Ly"]), float(data["Lz"]),
        int(data["Nx"]), int(data["Ny"]), int(data["Nz"]),
        float(data["dt_snap"]), float(data["obs_noise"]),
        str(data["eq_name"]),
    )


def build_weak_library_3d(f_snaps, ft_snaps, idx_sets, Lx, Ly, Lz, library_fn, m, p):
    """Build weak-form (G, b) for 3-D data with given library builder."""
    kernel = wf._cached_kernel_1d(m, p)
    rows_G, rows_b = [], []

    for s in range(len(f_snaps)):
        f = f_snaps[s]
        ft = ft_snaps[s]

        G, desc = library_fn(
            f, Lx, Ly, max_p=2, max_n=3, formulation="weak",
            m=m, p=p, use_ibp=True, Lz=Lz,
        )
        b_w = wf.build_weak_rhs(ft, kernel_x=kernel, kernel_y=kernel, m=m, p=p)

        Nx, Ny, Nz = f.shape
        Nv_x, Nv_y, Nv_z = Nx - 2 * m, Ny - 2 * m, Nz - 2 * m
        ix = np.arange(0, Nv_x, SAMPLE_STRIDE)
        iy = np.arange(0, Nv_y, SAMPLE_STRIDE)
        iz = np.arange(0, Nv_z, SAMPLE_STRIDE)
        IX, IY, IZ = np.meshgrid(ix, iy, iz, indexing="ij")
        flat_idx = (IX * Nv_y * Nv_z + IY * Nv_z + IZ).ravel()

        rng = np.random.default_rng(s)
        n_per = len(idx_sets[s])
        chosen = rng.choice(len(flat_idx), size=min(n_per, len(flat_idx)), replace=False)
        q_idx = flat_idx[chosen]

        rows_G.append(G[q_idx])
        rows_b.append(b_w[q_idx])

    G_all = np.vstack(rows_G).astype(float)
    b_all = np.concatenate(rows_b).astype(float)
    return G_all, b_all, desc


def build_pdefind_weak_3d(f_snaps, ft_snaps, idx_sets, Lx, Ly, Lz, m, p):
    """Build weak-form PDE-FIND library for 3-D data."""
    kernel = wf._cached_kernel_1d(m, p)
    rows_G, rows_b = [], []

    for s in range(len(f_snaps)):
        f = f_snaps[s]
        ft = ft_snaps[s]

        G, desc = build_pdefind_library(
            f, Lx, Ly, max_order=4, max_poly=3, formulation="weak",
            m=m, p=p, use_ibp=True, Lz=Lz,
        )
        b_w = wf.build_weak_rhs(ft, kernel_x=kernel, kernel_y=kernel, m=m, p=p)

        Nx, Ny, Nz = f.shape
        Nv_x, Nv_y, Nv_z = Nx - 2 * m, Ny - 2 * m, Nz - 2 * m
        ix = np.arange(0, Nv_x, SAMPLE_STRIDE)
        iy = np.arange(0, Nv_y, SAMPLE_STRIDE)
        iz = np.arange(0, Nv_z, SAMPLE_STRIDE)
        IX, IY, IZ = np.meshgrid(ix, iy, iz, indexing="ij")
        flat_idx = (IX * Nv_y * Nv_z + IY * Nv_z + IZ).ravel()

        rng = np.random.default_rng(s)
        n_per = len(idx_sets[s])
        chosen = rng.choice(len(flat_idx), size=min(n_per, len(flat_idx)), replace=False)
        q_idx = flat_idx[chosen]

        rows_G.append(G[q_idx])
        rows_b.append(b_w[q_idx])

    G_all = np.vstack(rows_G).astype(float)
    b_all = np.concatenate(rows_b).astype(float)
    return G_all, b_all, desc


def evaluate_native(w, desc, true_dict):
    """Evaluate recovered coefficients against ground truth in native basis.

    Returns (fp, fn, E_coeff) where fp/fn are computed by checking non-zero
    coefficients against the true dictionary.
    """
    eps = 0.15
    true_coeffs = np.zeros(len(desc))
    for i, d in enumerate(desc):
        if d in true_dict:
            true_coeffs[i] = true_dict[d]

    active_true = np.abs(true_coeffs) > 1e-12
    active_rec = np.abs(w) > 1e-10

    fp = int(np.sum(active_rec & ~active_true))
    fn = int(np.sum(active_true & ~active_rec))
    exact = int(fp == 0 and fn == 0)

    if np.any(active_true):
        E_coeff = np.linalg.norm(w - true_coeffs) / np.linalg.norm(true_coeffs)
    else:
        E_coeff = np.linalg.norm(w) if np.any(active_rec) else 0.0

    return fp, fn, exact, E_coeff


def main():
    # True Allen-Cahn PDE: f_t = eps^2 * Δf + f - f^3
    eps = 0.15
    true_multigraph = {"Δf": eps**2, "f": 1.0, "f³": -1.0}
    true_mod_multigraph = {"Δf": eps**2, "f": 1.0, "f³": -1.0}
    true_pdefind = {"f_xx": eps**2, "f_yy": eps**2, "f_zz": eps**2, "f": 1.0, "f³": -1.0}

    derivatives.set_backend("spectral")

    noise_labels = ["clean", "noise2pct", "noise5pct"]
    configs = [
        ("multigraph", build_multigraph_library, true_multigraph),
        ("mod_multigraph", build_mod_multigraph_library, true_mod_multigraph),
        ("pdefind", None, true_pdefind),  # pdefind handled specially
    ]

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = os.path.join(RESULTS_DIR, f"ac_3d_{timestamp}.csv")

    fieldnames = [
        "library", "noise", "n_active", "fp", "fn", "exact_recovery",
        "E_coeff", "E_data", "time_s",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for noise_label in noise_labels:
            npz_path = os.path.join(DATA_DIR, f"allen_cahn_3d_{noise_label}.npz")
            if not os.path.exists(npz_path):
                print(f"SKIP: {npz_path} not found. Run generate_3d_data.py first.")
                continue

            print(f"\n{'='*60}")
            print(f"Noise: {noise_label}")
            print(f"{'='*60}")

            snapshots, times, Lx, Ly, Lz, Nx, Ny, Nz, dt_snap, obs_noise, eq_name = \
                load_3d_dataset(npz_path)

            f_snaps, ft_snaps, idx_sets, _ = make_regression_data(
                snapshots, times, n_per_snap=N_PER_SNAP, skip_fraction=0.0, seed=0,
            )
            print(f"  Snaps: {len(f_snaps)}, grid: {Nx}x{Ny}x{Nz}, points: {len(f_snaps)*N_PER_SNAP}")

            for lib_name, lib_fn, true_dict in configs:
                t0 = time.time()

                if lib_name == "pdefind":
                    G, b, desc = build_pdefind_weak_3d(
                        f_snaps, ft_snaps, idx_sets, Lx, Ly, Lz, M, WEAK_P,
                    )
                else:
                    G, b, desc = build_weak_library_3d(
                        f_snaps, ft_snaps, idx_sets, Lx, Ly, Lz, lib_fn, M, WEAK_P,
                    )

                # Run RFE
                w = TrainRFE(G, b, normalize=2, lam=0.0)
                n_active = int(np.count_nonzero(np.abs(w) > 1e-10))
                fp, fn, exact, E_coeff = evaluate_native(w, desc, true_dict)

                # Compute E_data
                b_pred = G @ w
                E_data = float(np.linalg.norm(b_pred - b) / max(np.linalg.norm(b), 1e-14))

                dt = time.time() - t0

                row = {
                    "library": lib_name,
                    "noise": noise_label,
                    "n_active": n_active,
                    "fp": fp,
                    "fn": fn,
                    "exact_recovery": exact,
                    "E_coeff": f"{E_coeff:.6f}",
                    "E_data": f"{E_data:.6f}",
                    "time_s": f"{dt:.1f}",
                }
                writer.writerow(row)
                f.flush()

                marker = "✓" if exact else "✗"
                print(f"  [{lib_name:16s}] {marker} active={n_active} fp={fp} fn={fn} "
                      f"E_coeff={E_coeff:.4f} E_data={E_data:.4f} ({dt:.1f}s)")

    print(f"\nResults → {out_path}")


if __name__ == "__main__":
    main()
