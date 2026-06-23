"""
Generate and save simulation data for all benchmark equations.

Run once to pre-generate data for the unified ablation sweep.
Each (equation, noise_level) pair is saved as a .npz file in data/.

Usage
-----
    cd paper-experimental/code
    python scripts/generate_data.py
    python scripts/generate_data.py --equations allen_cahn swift_hohenberg
    python scripts/generate_data.py --noise 0.0 0.02
    python scripts/generate_data.py --data-dir data/
"""

import argparse
import csv
import datetime
import os
import sys
import time

import numpy as np

# ── path setup ──────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from simulations import (  # noqa: E402
    allen_cahn,
    cahn_hilliard,
    kuramoto_sivashinsky,
    porous_medium,
    swift_hohenberg,
)

DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
os.makedirs(DATA_DIR, exist_ok=True)

NOISE_LEVELS = [0.0, 0.02, 0.05]
NOISE_LABELS = {0.0: "clean", 0.02: "noise2pct", 0.05: "noise5pct"}

# ── equation configurations ──────────────────────────────────────────────────
EQUATIONS = [
    dict(
        name="allen_cahn",
        simulate_fn=allen_cahn,
        sim_kwargs=dict(N=64, L=2 * np.pi, T_end=8.0, eps=0.15, dt=0.01,
                        save_every=20, seed=42),
    ),
    dict(
        name="swift_hohenberg",
        simulate_fn=swift_hohenberg,
        sim_kwargs=dict(N=64, L=8 * np.pi, T_end=50.0, r=0.3, dt=0.1,
                        save_every=5, seed=42),
    ),
    dict(
        name="cahn_hilliard",
        simulate_fn=cahn_hilliard,
        sim_kwargs=dict(N=128, L=4 * np.pi, T_end=30.0, eps=0.2, dt=0.005,
                        save_every=30, seed=42),
    ),
    dict(
        name="porous_medium",
        simulate_fn=porous_medium,
        sim_kwargs=dict(N=64, L=6 * np.pi, T_end=2.0, dt=0.002,
                        save_every=25, seed=42),
    ),
    dict(
        name="ks_cellular",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=64, L=4 * np.pi, T_end=50.0, dt=0.1,
                        save_every=5, seed=42),
    ),
    dict(
        name="ks_weak_chaotic",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=64, L=8 * np.pi, T_end=100.0, dt=0.1,
                        save_every=5, seed=42),
    ),
    dict(
        name="ks_chaotic",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=64, L=16 * np.pi, T_end=100.0, dt=0.1,
                        save_every=5, seed=42),
    ),
]


def generate_one(eq_entry, obs_noise, data_dir):
    """Simulate one (equation, noise) pair and save as .npz."""
    name = eq_entry["name"]
    sim_fn = eq_entry["simulate_fn"]
    kwargs = eq_entry["sim_kwargs"].copy()
    kwargs["obs_noise"] = obs_noise

    noise_label = NOISE_LABELS[obs_noise]
    filename = f"{name}_{noise_label}.npz"
    path = os.path.join(data_dir, filename)

    if os.path.exists(path):
        print(f"  SKIP {filename} (already exists)", flush=True)
        return path, True

    t0 = time.time()
    snapshots, times = sim_fn(**kwargs)
    elapsed = time.time() - t0

    dt_snap = times[1] - times[0]
    N = snapshots.shape[1]
    L = kwargs["L"]

    np.savez_compressed(
        path,
        snapshots=snapshots.astype(np.float64),
        times=times.astype(np.float64),
        Lx=L,
        Ly=L,
        N=N,
        dt_snap=dt_snap,
        obs_noise=obs_noise,
    )

    n_snaps = snapshots.shape[0]
    has_nan = bool(np.any(np.isnan(snapshots)))
    print(
        f"  SAVED {filename}  "
        f"({n_snaps} snaps, {N}×{N}, "
        f"range=[{snapshots.min():.3f}, {snapshots.max():.3f}], "
        f"std={snapshots.std():.3f}, "
        f"nan={has_nan}, "
        f"{elapsed:.1f}s)",
        flush=True,
    )
    return path, False


def main():
    parser = argparse.ArgumentParser(
        description="Pre-generate simulation data for the unified sweep."
    )
    parser.add_argument(
        "--equations",
        nargs="+",
        default=None,
        help="Which equations to generate (default: all).",
    )
    parser.add_argument(
        "--noise",
        nargs="+",
        type=float,
        default=NOISE_LEVELS,
        help="Noise levels to generate (default: 0.0 0.02 0.05).",
    )
    parser.add_argument(
        "--data-dir",
        default=DATA_DIR,
        help="Output directory for .npz files.",
    )
    args = parser.parse_args()

    data_dir = args.data_dir
    os.makedirs(data_dir, exist_ok=True)

    eqs = EQUATIONS
    if args.equations is not None:
        eq_names = set(args.equations)
        eqs = [e for e in EQUATIONS if e["name"] in eq_names]
        if not eqs:
            print(f"No matching equations. Available: "
                  f"{[e['name'] for e in EQUATIONS]}")
            return

    noise_levels = args.noise

    print(f"Generating {len(eqs)} equations × {len(noise_levels)} noise levels "
          f"→ {len(eqs) * len(noise_levels)} files")
    print(f"Output directory: {data_dir}")
    print()

    # ── manifest ─────────────────────────────────────────────────────────
    manifest_path = os.path.join(data_dir, "manifest.csv")
    manifest_rows = []

    total_start = time.time()
    n_generated, n_skipped = 0, 0

    for eq_entry in eqs:
        name = eq_entry["name"]
        print(f"[{name}]", flush=True)
        for noise in noise_levels:
            path, skipped = generate_one(eq_entry, noise, data_dir)
            if skipped:
                n_skipped += 1
            else:
                n_generated += 1
            manifest_rows.append(dict(
                equation=name,
                noise=noise,
                noise_label=NOISE_LABELS[noise],
                file=os.path.basename(path),
                path=path,
            ))
        print()

    # Write manifest
    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "equation", "noise", "noise_label", "file", "path",
        ])
        writer.writeheader()
        writer.writerows(manifest_rows)

    total_elapsed = time.time() - total_start
    print(f"Done. {n_generated} generated, {n_skipped} skipped "
          f"in {total_elapsed:.0f}s.")
    print(f"Manifest → {manifest_path}")


if __name__ == "__main__":
    main()
