#!/usr/bin/env python
"""Generate 3-D Allen-Cahn data for PDE discovery experiments.

Produces .npz files with the same structure as the 2-D data.
"""

import os
import sys
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.simulations.allen_cahn import allen_cahn_3d

DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
os.makedirs(DATA_DIR, exist_ok=True)


def main():
    N = 48
    L = 2 * np.pi
    eps = 0.15
    T_end = 8.0
    dt = 0.01
    save_every = 20
    seed = 42

    print(f"Generating 3-D Allen-Cahn data: N={N}³, L={L:.1f}, eps={eps}, T={T_end}")
    print("  (this may take a few minutes)")

    # Generate clean data
    print("\n[1/3] Clean data...")
    snaps_clean, times = allen_cahn_3d(
        N=N, L=L, T_end=T_end, eps=eps, dt=dt,
        save_every=save_every, seed=seed, obs_noise=0.0,
    )
    print(f"  Shape: {snaps_clean.shape}, times: {len(times)} snaps")

    # No noise labels (clean, 2%, 5%)
    noise_configs = [
        ("clean", 0.0, snaps_clean),
        ("noise2pct", 0.02, snaps_clean.copy()),
        ("noise5pct", 0.05, snaps_clean.copy()),
    ]

    for label, noise_level, snaps in noise_configs:
        if noise_level > 0:
            rng = np.random.default_rng(seed + 1)
            for s in range(len(snaps)):
                std = np.std(snaps[s])
                snaps[s] += noise_level * std * rng.standard_normal(snaps[s].shape)

        fname = f"allen_cahn_3d_{label}.npz"
        path = os.path.join(DATA_DIR, fname)
        np.savez_compressed(
            path,
            snapshots=snaps,
            times=times,
            Lx=L, Ly=L, Lz=L,
            Nx=N, Ny=N, Nz=N,
            dt_snap=dt * save_every,
            obs_noise=noise_level,
            eq_name="allen_cahn_3d",
        )
        print(f"[{label}] Saved → {path}  ({os.path.getsize(path)/1e6:.1f} MB)")

    print("\nDone.")


if __name__ == "__main__":
    main()
