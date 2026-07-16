"""
Generate diagnostic CH datasets for noise-robustness investigation.

Varies one parameter at a time (N, L, T_end, dt_snap, eps) across three
noise levels (0%, 0.5%, 1%) to identify which data hyperparameter most
affects Cahn-Hilliard recovery at low noise.

Usage
-----
    python scripts/generate_ch_diagnostics.py
    python scripts/generate_ch_diagnostics.py --data-dir data/ch_diagnostics/
"""

import argparse
import os
import sys
import time

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.simulations import cahn_hilliard  # noqa: E402

DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "ch_diagnostics")
os.makedirs(DATA_DIR, exist_ok=True)

NOISE_LEVELS = [0.0, 0.005, 0.01]
NOISE_LABELS = {0.0: "clean", 0.005: "noise0p5pct", 0.01: "noise1pct"}

# Default CH simulation parameters
DEFAULT = dict(N=128, L=4 * np.pi, T_end=30.0, eps=0.2, dt=0.005,
               save_every=30, seed=42)

# Generate one dataset
# ---------------------------------------------------------------------------


def generate_one(name, N, L, T_end, eps, dt, save_every, seed, obs_noise, data_dir):
    filename = f"ch_diag_{name}_{NOISE_LABELS[obs_noise]}.npz"
    path = os.path.join(data_dir, filename)

    if os.path.exists(path):
        print(f"  SKIP {filename} (already exists)", flush=True)
        return

    t0 = time.time()
    snapshots, times = cahn_hilliard(
        N=N, L=L, T_end=T_end, eps=eps, dt=dt,
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
        eq_name="ch",
        # diagnostic metadata
        sim_N=N, sim_L=L, sim_T_end=T_end, sim_eps=eps,
        sim_dt=dt, sim_save_every=save_every, sim_dt_snap=dt_snap,
    )

    n_snaps = snapshots.shape[0]
    has_nan = bool(np.any(np.isnan(snapshots)))
    print(
        f"  SAVED {filename}  ({n_snaps} snaps, {N}×{N}, "
        f"range=[{snapshots.min():.3f}, {snapshots.max():.3f}], "
        f"nan={has_nan}, {elapsed:.1f}s)",
        flush=True,
    )


# Sweep configurations: one parameter varies, rest at default
# ---------------------------------------------------------------------------

def configs():
    """Yield (name_short, kwargs_dict) for each varied-parameter config.

    The kwarg dict contains the full set of simulation parameters.
    """
    d = DEFAULT.copy()

    # --- Vary N (spatial resolution) ---
    for N in [64, 256]:
        kw = d.copy(); kw["N"] = N
        yield f"N{N}", kw

    # --- Vary L (domain size) ---
    for L in [2 * np.pi, 8 * np.pi]:
        kw = d.copy(); kw["L"] = L
        yield f"L{int(L / np.pi)}pi", kw

    # --- Vary T_end (simulation duration) ---
    for T_end in [15.0, 60.0, 100.0]:
        kw = d.copy(); kw["T_end"] = T_end
        yield f"T{T_end:.0f}", kw

    # --- Vary dt_snap (via save_every, dt=0.005 fixed) ---
    for save_every in [10, 20, 60]:  # dt_snap = 0.05, 0.10, 0.30
        kw = d.copy(); kw["save_every"] = save_every
        dt_snap = kw["dt"] * save_every
        yield f"dtsnap{dt_snap:.2f}", kw

    # --- Vary eps (interface width) ---
    for eps in [0.1, 0.4]:
        kw = d.copy(); kw["eps"] = eps
        yield f"eps{eps}", kw


# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate CH diagnostic datasets (vary one param at a time)."
    )
    parser.add_argument("--data-dir", default=DATA_DIR,
                        help="Output directory for .npz files.")
    args = parser.parse_args()

    data_dir = args.data_dir
    os.makedirs(data_dir, exist_ok=True)

    total = 0
    for name_short, kw in configs():
        for obs_noise in NOISE_LEVELS:
            noise_label = NOISE_LABELS[obs_noise]
            print(f"\n[{name_short}] noise={obs_noise:.3f} ({noise_label})", flush=True)
            generate_one(
                name=name_short,
                N=kw["N"], L=kw["L"], T_end=kw["T_end"],
                eps=kw["eps"], dt=kw["dt"],
                save_every=kw["save_every"], seed=kw["seed"],
                obs_noise=obs_noise,
                data_dir=data_dir,
            )
            total += 1

    print(f"\nDone. {total} datasets requested.", flush=True)


if __name__ == "__main__":
    main()
