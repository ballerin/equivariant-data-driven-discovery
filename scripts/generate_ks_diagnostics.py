"""
Generate diagnostic KS (weak chaotic) datasets for hyperparameter investigation.

Varies one parameter at a time (N, L, T_end, dt_snap) across three
noise levels (0%, 2%, 5%) to identify which data hyperparameter most
affects Kuramoto-Sivashinsky recovery.

Usage
-----
    python scripts/generate_ks_diagnostics.py
    python scripts/generate_ks_diagnostics.py --data-dir data/ks_diagnostics/
"""

import argparse
import os
import sys
import time

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.simulations import kuramoto_sivashinsky  # noqa: E402

DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "ks_diagnostics")
os.makedirs(DATA_DIR, exist_ok=True)

NOISE_LEVELS = [0.0, 0.02, 0.05]
NOISE_LABELS = {0.0: "clean", 0.02: "noise2pct", 0.05: "noise5pct"}

# Default KS weak-chaotic parameters
DEFAULT = dict(N=64, L=8 * np.pi, T_end=100.0, dt=0.1,
               save_every=5, seed=42)


def generate_one(name, N, L, T_end, dt, save_every, seed, obs_noise, data_dir):
    filename = f"ks_diag_{name}_{NOISE_LABELS[obs_noise]}.npz"
    path = os.path.join(data_dir, filename)

    if os.path.exists(path):
        print(f"  SKIP {filename} (already exists)", flush=True)
        return

    t0 = time.time()
    snapshots, times = kuramoto_sivashinsky(
        N=N, L=L, T_end=T_end, dt=dt,
        save_every=save_every, seed=seed, obs_noise=obs_noise,
    )
    elapsed = time.time() - t0

    dt_snap = times[1] - times[0]

    np.savez_compressed(
        path,
        snapshots=snapshots.astype(np.float64),
        times=times.astype(np.float64),
        Lx=L, Ly=L, N=N, dt_snap=dt_snap,
        obs_noise=obs_noise,
        eq_name="ks_wchaos",
        sim_N=N, sim_L=L, sim_T_end=T_end,
        sim_dt=dt, sim_save_every=save_every, sim_dt_snap=dt_snap,
    )

    n_snaps = snapshots.shape[0]
    has_nan = bool(np.any(np.isnan(snapshots)))
    print(
        f"  SAVED {filename}  ({n_snaps} snaps, {N}x{N}, "
        f"range=[{snapshots.min():.3f}, {snapshots.max():.3f}], "
        f"nan={has_nan}, {elapsed:.1f}s)",
        flush=True,
    )


# Sweep configurations: vary one parameter at a time
CONFIGS = []

# --- vary N ---
for N in [32, 128]:
    CONFIGS.append(dict(name=f"N{N}", N=N, **{k: v for k, v in DEFAULT.items()
                   if k != "N"}))

# --- vary L ---
for L in [4 * np.pi, 16 * np.pi]:
    CONFIGS.append(dict(name=f"L{4 if L < 20 else 16}pi", L=L,
                   **{k: v for k, v in DEFAULT.items() if k != "L"}))

# --- vary T_end ---
for T_end in [50.0, 200.0]:
    CONFIGS.append(dict(name=f"T{int(T_end)}", T_end=T_end,
                   **{k: v for k, v in DEFAULT.items() if k != "T_end"}))

# --- vary dt_snap via save_every ---
for save_every in [2, 20]:
    dt_snap_label = DEFAULT["dt"] * save_every
    CONFIGS.append(dict(
        name=f"dtsnap{dt_snap_label:.2f}".replace(".", "p"),
        save_every=save_every,
        **{k: v for k, v in DEFAULT.items() if k != "save_every"},
    ))


def main():
    parser = argparse.ArgumentParser(
        description="Generate diagnostic KS weak-chaotic datasets."
    )
    parser.add_argument(
        "--data-dir", default=DATA_DIR,
        help="Output directory for .npz files.",
    )
    args = parser.parse_args()
    data_dir = args.data_dir
    os.makedirs(data_dir, exist_ok=True)

    total = len(CONFIGS) * len(NOISE_LEVELS)
    print(f"Generating {len(CONFIGS)} configs x {len(NOISE_LEVELS)} noise"
          f" levels = {total} files")
    print(f"Output: {data_dir}")
    print()

    for cfg in CONFIGS:
        name = cfg.pop("name")
        for obs_noise in NOISE_LEVELS:
            generate_one(
                name=name, obs_noise=obs_noise,
                data_dir=data_dir, **cfg,
            )

    print(f"\nDone. Generated {total} files → {data_dir}")


if __name__ == "__main__":
    main()
