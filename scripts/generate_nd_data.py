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

from equivariant_pde_find import derivatives
from equivariant_pde_find.simulations.allen_cahn import allen_cahn_nd
from equivariant_pde_find.simulations.cahn_hilliard import cahn_hilliard_nd
from equivariant_pde_find.simulations.swift_hohenberg import swift_hohenberg_nd
from equivariant_pde_find.simulations.porous_medium import porous_medium_nd

DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
os.makedirs(DATA_DIR, exist_ok=True)


def _rhs_ac(f, L, eps):
    return eps**2 * derivatives.laplacian(f, L, L, L) + f - f**3


def _rhs_ch(f, L, eps):
    return (-(eps**2) * derivatives.bilaplacian(f, L, L, L)
            + derivatives.laplacian(f**3 - f, L, L, L))


def _rhs_sh(f, L, r):
    return (r - 1.0) * f - 2.0 * derivatives.laplacian(f, L, L, L) \
           - derivatives.bilaplacian(f, L, L, L) - f**3


EQUATION_CONFIGS = {
    "ac": {
        "sim_fn": allen_cahn_nd,
        "name": "ac",
        "L": 2 * np.pi,
        "eps": 0.15,            # used for PDE residual check
        "sim_kwargs": {"eps": 0.15},
        "rhs": _rhs_ac,
    },
    "ch": {
        "sim_fn": cahn_hilliard_nd,
        "name": "ch",
        "L": 8 * np.pi,
        "eps": 0.2,
        "sim_kwargs": {"eps": 0.2},
        "rhs": _rhs_ch,
    },
    "sh": {
        "sim_fn": swift_hohenberg_nd,
        "name": "sh",
        "L": 8 * np.pi,
        "eps": 0.3,             # used as 'r' for PDE residual check
        "sim_kwargs": {"r": 0.3, "T_end": 50.0, "save_every": 5},
        "rhs": _rhs_sh,
    },
    "pme": {
        "sim_fn": porous_medium_nd,
        "name": "pme",
        "L": 6 * np.pi,
        "eps": None,            # PME has no interface parameter
        "sim_kwargs": {},
        "rhs": None,            # PME residual not checked (nonlinear diffusion)
    },
}

NOISE_LABELS = {0.0: "clean", 0.01: "noise1pct", 0.02: "noise2pct", 0.05: "noise5pct"}


def pde_residual(snaps, times, L, eps, rhs_fn):
    """Relative residual ||f_t - RHS_true(f)|| / ||f_t|| of the saved data.

    Uses the same central time difference the regression pipeline uses
    (``simulations/regression_data.py``) and spectral spatial derivatives.
    A dataset that does not satisfy its own PDE cannot support a recovery
    claim -- the Cahn-Hilliard data misses its by 33 % (audit F-04), which
    is why this is now reported at generation time.
    """
    if len(snaps) < 3:
        return float("nan")
    dt_snap = times[1] - times[0]
    rels = []
    for i in range(1, len(snaps) - 1):
        ft = (snaps[i + 1] - snaps[i - 1]) / (2.0 * dt_snap)
        r = ft - rhs_fn(snaps[i], L, eps)
        rels.append(np.linalg.norm(r) / max(np.linalg.norm(ft), 1e-30))
    return float(np.median(rels))


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
    # Optional simulator overrides.  All default to None, i.e. the simulator's
    # own defaults, so an argument-free run reproduces the previous behaviour.
    parser.add_argument(
        "--N", type=int, default=None,
        help="Grid points per axis (default: simulator's, 64 at ndim=2 else 48).",
    )
    parser.add_argument(
        "--eps", type=float, default=None,
        help="Interface-width parameter (default: simulator's).",
    )
    parser.add_argument(
        "--T-end", type=float, default=None, dest="T_end",
        help="Final time (default: simulator's).",
    )
    parser.add_argument(
        "--save-every", type=int, default=None,
        help="Snapshot interval in steps (default: simulator's).",
    )
    parser.add_argument(
        "--keep-last", type=int, default=None,
        help=(
            "Keep only the final K snapshots.  Snapshot *spacing* sets the "
            "accuracy of the central-difference f_t, so coarsen with this "
            "rather than with --save-every: it keeps the tight spacing and "
            "only bounds file size, which matters at high ndim where one "
            "snapshot is hundreds of MB."
        ),
    )
    parser.add_argument(
        "--tag", type=str, default="",
        help=(
            "Filename tag: <equation>_nd<ndim>_<tag>_<label>.npz.  Use one "
            "whenever non-default parameters are passed, so the new files do "
            "not collide with the canonical datasets."
        ),
    )
    args = parser.parse_args()

    if args.tag == "" and any(
        v is not None for v in (args.N, args.eps, args.T_end, args.save_every)
    ):
        parser.error(
            "non-default simulator parameters were given without --tag; that "
            "would overwrite the canonical datasets. Pass --tag."
        )

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
            sim_kwargs = {"ndim": ndim, "seed": args.seed, "obs_noise": 0.0}
            # Merge equation-specific defaults (e.g. eps for AC/CH, r for SH)
            sim_kwargs.update(cfg.get("sim_kwargs", {}))
            # Always pass L so the domain size matches the RHS checker
            sim_kwargs.setdefault("L", cfg["L"])
            # CLI overrides for common simulator parameters
            for key, value in (("N", args.N), ("T_end", args.T_end),
                               ("save_every", args.save_every)):
                if value is not None:
                    sim_kwargs[key] = value

            t0 = time.time()
            snaps_clean, times = sim_fn(**sim_kwargs)
            dt = time.time() - t0

            if args.keep_last is not None and args.keep_last < len(snaps_clean):
                # .copy() matters: a slice is a *view*, which would keep the
                # whole snapshot array alive.  At ndim=5 that is ~14 GB held
                # for the rest of the run instead of ~0.8 GB.
                snaps_clean = snaps_clean[-args.keep_last:].copy()
                times = times[-args.keep_last:].copy()

            n_snaps = snaps_clean.shape[0]
            shape_str = "x".join(str(s) for s in snaps_clean.shape[1:])
            print(f"  Clean: {n_snaps} snaps, grid {shape_str}, "
                  f"std={snaps_clean.std():.3f}, {dt:.1f}s")

            N = snaps_clean.shape[1]
            L = cfg["L"]
            eps = args.eps if args.eps is not None else cfg["eps"]

            # Does the data satisfy the equation it is labelled with?
            if cfg["rhs"] is not None:
                with derivatives.use_backend("spectral"):
                    resid = pde_residual(snaps_clean, times, L, eps, cfg["rhs"])
                flag = "" if not (resid > 0.01) else "   <-- WARNING: > 1 %"
                print(f"  PDE residual (median, clean): {resid:.3e}{flag}")
            else:
                print(f"  PDE residual: skipped (no RHS checker for {eq_key})")

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

                tag_part = f"{args.tag}_" if args.tag else ""
                fname = f"{eq_name}_nd{ndim}_{tag_part}{label}.npz"
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
                    # eps is a ground-truth coefficient, so it must travel with
                    # the data: a sweep that assumes the default would score
                    # against the wrong truth vector.
                    "eps": eps,
                    "pde_residual": resid,
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
