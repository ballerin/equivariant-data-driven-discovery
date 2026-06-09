"""
WSINDy hyperparameter sweep: weak-form regression with optional IBP.

Sweeps over: equations, noise levels, IBP on/off, library type (equivariant
/ standard), derivative backend, and d_tol.  Saves one CSV per run.

Usage
-----
    cd paper-experimental/code
    python scripts/wsindy_sweep.py --equations allen_cahn swift_hohenberg --noise 0.0 0.02
    python scripts/wsindy_sweep.py --equations kuramoto_sivashinsky --noise 0.02 --ibp
    python scripts/wsindy_sweep.py  # all equations, both IBP modes, clean + 2% + 5%
"""

import argparse
import csv
import datetime
import itertools
import os
import sys
import time

import numpy as np

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        print(kwargs.get("desc", "") + f" ({len(list(iterable))} items)...")
        return iterable

# ── path setup ──────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, CODE_DIR)

import derivatives  # noqa: E402
from library import build_equivariant_library, build_standard_library  # noqa: E402
from regression import TrainSTRidge  # noqa: E402
from simulations import (  # noqa: E402
    allen_cahn,
    cahn_hilliard,
    kuramoto_sivashinsky,
    make_regression_data,
    porous_medium,
    swift_hohenberg,
)
import wsindy as ws  # noqa: E402

RESULTS_DIR = os.path.join(CODE_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyperparameter grid ──────────────────────────────────────────────────────
LAMS = [1e-5]
D_TOLS = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, 5.0, 15.0, 20.0]
BACKENDS = ["spectral", "polynomial"]
WSINDY_M = [8]       # test function half-support
WSINDY_P = [4]       # test function polynomial degree
SAMPLE_STRIDE = 4    # query-point stride in weak grid
IBP_MODES = [False, True]

# ── ground-truth coefficient vectors ────────────────────────────────────────


def _find_op(desc, name):
    try:
        return desc.index(name)
    except ValueError:
        raise KeyError(f"Operator {name!r} not found.")


def _true_allen_cahn(desc_eq, eps=0.15):
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, "f")] = 1.0
    w[_find_op(desc_eq, "f³")] = -1.0
    w[_find_op(desc_eq, "Δf")] = eps ** 2
    return w


def _true_swift_hohenberg(desc_eq, r=0.3):
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, "f")] = r - 1.0
    w[_find_op(desc_eq, "Δf")] = -2.0
    w[_find_op(desc_eq, "Δ²f")] = -1.0
    w[_find_op(desc_eq, "f³")] = -1.0
    return w


def _true_pme(desc_eq):
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, "f²·Δf")] = 1.0
    w[_find_op(desc_eq, "f·‖∇f‖²")] = 2.0
    return w


def _true_cahn_hilliard(desc_eq, eps=0.2):
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, "Δf")] = -1.0
    w[_find_op(desc_eq, "f²·Δf")] = 3.0
    w[_find_op(desc_eq, "f·‖∇f‖²")] = 6.0
    w[_find_op(desc_eq, "Δ²f")] = -(eps ** 2)
    return w


def _true_ks(desc_eq):
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, "‖∇f‖²")] = -0.5
    w[_find_op(desc_eq, "Δf")] = -1.0
    w[_find_op(desc_eq, "Δ²f")] = -1.0
    return w


def _mask_pme(f):
    return f.ravel() > 0.05


# ── equation definitions ─────────────────────────────────────────────────────
EQUATIONS = [
    dict(
        name="allen_cahn",
        simulate_fn=allen_cahn,
        sim_kwargs=dict(N=64, L=2 * np.pi, T_end=8.0, eps=0.15, dt=0.01, save_every=20, seed=42),
        n_per_snap=200,
        true_eq_fn=lambda d: _true_allen_cahn(d, eps=0.15),
        skip_fraction=0.0,
        mask_fn=None,
    ),
    dict(
        name="swift_hohenberg",
        simulate_fn=swift_hohenberg,
        sim_kwargs=dict(N=64, L=8 * np.pi, T_end=50.0, r=0.3, dt=0.1, save_every=5, seed=42),
        n_per_snap=200,
        true_eq_fn=lambda d: _true_swift_hohenberg(d, r=0.3),
        skip_fraction=0.1,
        mask_fn=None,
    ),
    dict(
        name="porous_medium",
        simulate_fn=porous_medium,
        sim_kwargs=dict(N=64, L=6 * np.pi, T_end=2.0, dt=0.002, save_every=25, seed=42),
        n_per_snap=200,
        true_eq_fn=_true_pme,
        skip_fraction=0.0,
        mask_fn=_mask_pme,
    ),
    dict(
        name="cahn_hilliard",
        simulate_fn=cahn_hilliard,
        sim_kwargs=dict(N=128, L=4 * np.pi, T_end=30.0, eps=0.2, dt=0.005, save_every=30, seed=42),
        n_per_snap=200,
        true_eq_fn=lambda d: _true_cahn_hilliard(d, eps=0.2),
        skip_fraction=0.0,
        mask_fn=None,
    ),
    dict(
        name="kuramoto_sivashinsky",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=64, L=16 * np.pi, T_end=200.0, dt=0.1, save_every=5, seed=42),
        n_per_snap=200,
        true_eq_fn=_true_ks,
        skip_fraction=0.3,
        mask_fn=None,
    ),
]

DEFAULT_EQUATIONS = ["allen_cahn", "swift_hohenberg", "porous_medium", "cahn_hilliard"]


# ── build weak-form matrices ────────────────────────────────────────────────


def build_wsindy_matrices(
    snapshots, times, L, n_per_snap, skip_fraction,
    mask_fn, m, p, sample_stride, use_ibp,
):
    """Build weak-form regression matrices for equivariant and standard libraries.

    Returns (mats, Ft, desc_eq, desc_std) where mats is keyed by
    (backend, library_type).
    """
    f_snaps, ft_snaps, idx_sets, _ = make_regression_data(
        snapshots, times, n_per_snap=n_per_snap, seed=0,
    )
    n = len(f_snaps)
    start = int(np.floor(skip_fraction * n))
    f_snaps = f_snaps[start:]
    ft_snaps = ft_snaps[start:]
    idx_sets = idx_sets[start:]

    if mask_fn is not None:
        rng = np.random.default_rng(1)
        new_idx = []
        for s in range(len(f_snaps)):
            valid = np.where(mask_fn(f_snaps[s]))[0]
            if len(valid) <= n_per_snap:
                new_idx.append(valid)
            else:
                new_idx.append(np.sort(rng.choice(valid, size=n_per_snap, replace=False)))
        idx_sets = new_idx

    mats = {}
    desc_eq = desc_std = None

    for backend in BACKENDS:
        derivatives.set_backend(backend)
        eq_rows, std_rows, rhs_rows = [], [], []

        for s in range(len(f_snaps)):
            f = f_snaps[s]
            ft = ft_snaps[s]

            # Equivariant weak library
            Theta_w_eq, d_eq = ws.build_weak_library_2d(
                f, L, L, max_p=2, max_n=3, m=m, p=p, use_ibp=use_ibp,
            )
            # Standard weak library
            Theta_w_std, d_std = ws.build_weak_standard_library_2d(
                f, L, L, max_order=4, max_poly=3, m=m, p=p, use_ibp=use_ibp,
            )
            # Weak RHS
            b_w = ws.build_weak_rhs_2d(ft, m=m, p=p)

            # Subsample: use random indices mapped to valid grid
            Nv_x, Nv_y = f.shape[0] - 2 * m, f.shape[1] - 2 * m
            Nv = Nv_x * Nv_y

            if sample_stride > 1:
                ix = np.arange(0, Nv_x, sample_stride)
                iy = np.arange(0, Nv_y, sample_stride)
                IX, IY = np.meshgrid(ix, iy, indexing="ij")
                flat_idx = (IX * Nv_y + IY).ravel()
                nq = len(flat_idx)
                # Take a random subset of strided points
                rng = np.random.default_rng(s)
                chosen = rng.choice(nq, size=min(n_per_snap, nq), replace=False)
                q_idx = flat_idx[chosen]
            else:
                q_idx = idx_sets[s]

            eq_rows.append(Theta_w_eq[q_idx])
            std_rows.append(Theta_w_std[q_idx])
            rhs_rows.append(b_w[q_idx])

            if desc_eq is None:
                desc_eq = d_eq
                desc_std = d_std

        mats[(backend, "equivariant")] = np.vstack(eq_rows).astype(float)
        mats[(backend, "standard")] = np.vstack(std_rows).astype(float)
        mats[(backend, "_rhs")] = np.concatenate(rhs_rows).astype(float)

    Ft = mats[("spectral", "_rhs")]  # RHS is the same for both backends
    return mats, Ft, desc_eq, desc_std


# ── sweep one configuration ─────────────────────────────────────────────────


def sweep_one(Theta, Ft, true_eq, lam, d_tol):
    w = TrainSTRidge(Theta, Ft, lam, d_tol)
    act = np.abs(true_eq) > 1e-10
    nz = int(np.count_nonzero(np.abs(w) > 1e-10))
    fp = int(np.sum((np.abs(w) > 1e-10) & ~act))
    fn = int(np.sum((np.abs(w) < 1e-10) & act))
    exact = int(fp == 0 and fn == 0)
    found = act & (np.abs(w) > 1e-10)
    err = float(np.max(np.abs(w[act] - true_eq[act]))) if found.any() else float("nan")
    return nz, fp, fn, exact, err


# ── main ─────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(description="WSINDy sweep.")
    parser.add_argument("--equations", nargs="+", default=DEFAULT_EQUATIONS)
    parser.add_argument("--noise", type=float, nargs="+", default=[0.0, 0.02, 0.05],
                        help="Noise levels to sweep (e.g. 0.0 0.02 0.05).")
    parser.add_argument("--ibp", action="store_true", default=False,
                        help="If set, run only use_ibp=True. Otherwise run both False and True.")
    parser.add_argument("--name", default="", help="Tag for output filename.")
    args = parser.parse_args()

    eq_names = set(args.equations)
    equations = [e for e in EQUATIONS if e["name"] in eq_names]
    if not equations:
        print(f"No matching equations. Available: {[e['name'] for e in EQUATIONS]}")
        return

    ibp_modes = [True] if args.ibp else IBP_MODES
    noise_levels = sorted(args.noise)

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"_{args.name}" if args.name else ""
    out_path = os.path.join(RESULTS_DIR, f"wsindy_sweep_{timestamp}{tag}.csv")
    print(f"Output → {out_path}")

    fieldnames = [
        "equation", "backend", "library", "noise", "use_ibp", "m", "p",
        "lam", "d_tol", "n_active", "fp", "fn", "exact_recovery", "max_coeff_err",
    ]

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for eq in equations:
            name = eq["name"]
            L = eq["sim_kwargs"].get("L", 2 * np.pi)

            for noise in noise_levels:
                print(f"\n{'='*60}")
                print(f"  {name}  |  noise={noise}")
                print(f"{'='*60}")

                t0 = time.time()
                print(f"  Simulating...", flush=True)
                snapshots, times = eq["simulate_fn"](
                    **eq["sim_kwargs"], obs_noise=noise,
                )
                print(f"  done in {time.time()-t0:.1f}s  ({len(snapshots)} snaps)", flush=True)

                for ibp in ibp_modes:
                    ibp_label = "IBP" if ibp else "direct"
                    print(f"\n  --- {ibp_label} ---")

                    for m_val in WSINDY_M:
                        p_val = WSINDY_P[0]

                        t0 = time.time()
                        print(f"  Building weak libraries (m={m_val}, p={p_val})...", flush=True)
                        mats, Ft, desc_eq, desc_std = build_wsindy_matrices(
                            snapshots, times, L,
                            eq["n_per_snap"], eq["skip_fraction"],
                            eq.get("mask_fn"), m_val, p_val,
                            SAMPLE_STRIDE, ibp,
                        )
                        print(f"  done in {time.time()-t0:.1f}s  ({len(Ft)} points)", flush=True)

                        true_eq = eq["true_eq_fn"](desc_eq)
                        true_std = np.zeros(len(desc_std))

                        # Condition numbers
                        for (backend, lib), Theta in mats.items():
                            if lib != "_rhs":
                                print(f"  cond({backend:10s} {lib:12s}): {np.linalg.cond(Theta):.2e}")

                        # Sweep
                        for backend in BACKENDS:
                            combos = list(itertools.product(LAMS, D_TOLS))
                            for lam, d_tol in tqdm(combos, desc=f"  {name} {backend} {ibp_label}", unit="co", leave=False):
                                # Equivariant
                                nz, fp, fn, ex, err = sweep_one(
                                    mats[(backend, "equivariant")], Ft, true_eq, lam, d_tol,
                                )
                                writer.writerow(dict(
                                    equation=name, backend=backend, library="equivariant",
                                    noise=noise, use_ibp=int(ibp), m=m_val, p=p_val,
                                    lam=lam, d_tol=d_tol, n_active=nz,
                                    fp=fp, fn=fn, exact_recovery=ex, max_coeff_err=err,
                                ))

                                # Standard
                                nz, fp, fn, ex, err = sweep_one(
                                    mats[(backend, "standard")], Ft, true_std, lam, d_tol,
                                )
                                writer.writerow(dict(
                                    equation=name, backend=backend, library="standard",
                                    noise=noise, use_ibp=int(ibp), m=m_val, p=p_val,
                                    lam=lam, d_tol=d_tol, n_active=nz,
                                    fp=fp if fp != fp else "nan",
                                    fn=fn if fn == fn else "nan",
                                    exact_recovery=ex, max_coeff_err=err,
                                ))
                            fh.flush()

    print(f"\nDone → {out_path}")


if __name__ == "__main__":
    main()
