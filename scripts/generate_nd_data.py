#!/usr/bin/env python
"""Generate n-dimensional simulation data for PDE discovery experiments.

Produces .npz files for Allen-Cahn and Cahn-Hilliard at specified
dimensions and noise levels.

Usage
-----
    python scripts/generate_nd_data.py
    python scripts/generate_nd_data.py --equations ac ch --ndim 2 3
    python scripts/generate_nd_data.py --noise 0.0 0.02
"""

import argparse
import os
import sys
import time

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.simulations.allen_cahn import allen_cahn_nd
from equivariant_pde_find.simulations.cahn_hilliard import cahn_hilliard_nd

DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
os.makedirs(DATA_DIR, exist_ok=True)

EQUATION_CONFIGS = {
    "ac": {
        "sim_fn": allen_cahn_nd,
        "name": "allen_cahn",
        "L": 2 * np.pi,
    },
    "ch": {
        "sim_fn": cahn_hilliard_nd,
        "name": "cahn_hilliard",
        "L": 8 * np.pi,
    },
}

NOISE_LABELS = {0.0: "clean", 0.02: "noise2pct", 0.05: "noise5pct"}


def main():
    parser = argparse.ArgumentParser(
        description="Generate n-dimensional PDE simulation data."
    )
    parser.add_argument(
        "--equations", nargs="+", default=["ac", "ch"],
        help="Equations to generate (ac, ch).",
    )
    parser.add_argument(
        "--ndim", nargs="+", type=int, default=[2, 3],
        help="Spatial dimensions (default: 2 3).",
    )
    parser.add_argument(
        "--noise", nargs="+", type=float, default=[0.0, 0.02, 0.05],
        help="Observation noise levels (default: 0.0 0.02 0.05).",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="RNG seed for simulation and noise.",
    )
    args = parser.parse_args()

    for eq_key in args.equations:
        if eq_key not in EQUATION_CONFIGS:
            print(f"Unknown equation: {eq_key}.  Choices: {list(EQUATION_CONFIGS)}")
            continue
        cfg = EQUATION_CONFIGS[eq_key]
        sim_fn = cfg["sim_fn"]
        eq_name = cfg["name"]

        for ndim in args.ndim:
            print(f"\n{'='*60}")
            print(f"Equation: {eq_name}  ndim={ndim}")
            print(f"{'='*60}")

            # Generate clean data once
            t0 = time.time()
            snaps_clean, times = sim_fn(ndim=ndim, seed=args.seed, obs_noise=0.0)
            dt = time.time() - t0
            n_snaps = snaps_clean.shape[0]
            shape_str = "x".join(str(s) for s in snaps_clean.shape[1:])
            print(f"  Clean: {n_snaps} snaps, grid {shape_str}, "
                  f"std={snaps_clean.std():.3f}, {dt:.1f}s")

            N = snaps_clean.shape[1]
            L = cfg["L"]

            for noise_level in args.noise:
                label = NOISE_LABELS[noise_level]

                if noise_level > 0:
                    rng = np.random.default_rng(args.seed + 1)
                    snaps = snaps_clean.copy()
                    for s in range(len(snaps)):
                        std = np.std(snaps[s])
                        snaps[s] += noise_level * std * rng.standard_normal(
                            snaps[s].shape
                        )
                else:
                    snaps = snaps_clean

                fname = f"{eq_name}_nd{ndim}_{label}.npz"
                path = os.path.join(DATA_DIR, fname)

                # Build metadata dict with ndim-aware keys
                metadata = {
                    "snapshots": snaps.astype(np.float64),
                    "times": times.astype(np.float64),
                    "dt_snap": float(times[1] - times[0]) if len(times) > 1 else 0.0,
                    "obs_noise": noise_level,
                    "eq_name": eq_name,
                    "ndim": ndim,
                    "N": N,
                    "L": L,
                }
                # Add per-axis params (Lx, Ly, Lz, ...) for compatibility
                axis_labels = ["x", "y", "z", "w", "v", "u", "t", "s"]
                for i in range(ndim):
                    metadata[f"L{axis_labels[i]}"] = L
                    metadata[f"N{axis_labels[i]}"] = N

                np.savez_compressed(path, **metadata)
                size_mb = os.path.getsize(path) / 1e6
                print(f"  [{label:>10s}] -> {fname}  ({size_mb:.1f} MB)")

    print("\nDone.")


if __name__ == "__main__":
    main()
