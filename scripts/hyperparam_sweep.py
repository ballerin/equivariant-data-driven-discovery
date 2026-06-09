"""
Hyperparameter sweep: (lambda, d_tol) grid for STRidge on the equivariant
and standard libraries across all five benchmark equations.

Results are saved to  ../results/hyperparam_sweep.csv
with columns:
    equation, backend, library, lam, d_tol,
    n_active, fp, fn, exact_recovery, max_coeff_err

Usage
-----
    cd paper-experimental/code
    python scripts/hyperparam_sweep.py
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
    # Graceful fallback: tqdm not installed
    def tqdm(iterable, **kwargs):
        desc = kwargs.get("desc", "")
        items = list(iterable)
        print(f"{desc} ({len(items)} items)...")
        return items


# ── path setup ──────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, CODE_DIR)

import derivatives  # noqa: E402
from library import build_equivariant_library, build_standard_library, spectral_lowpass  # noqa: E402
from regression import TrainSTRidge, TrainSR3  # noqa: E402
from simulations import (  # noqa: E402
    allen_cahn,
    cahn_hilliard,
    kuramoto_sivashinsky,
    make_regression_data,
    porous_medium,
    swift_hohenberg,
)

RESULTS_DIR = os.path.join(CODE_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyperparameter grid ──────────────────────────────────────────────────────
# Note: lam is not applied in the STRidge inner loop (historical PDE-FIND
# convention), so sweeping over lam values produces identical results.
# A single representative value is used here.
LAMS = [1e-5]
# STRidge thresholds *normalized* coefficients (column-scaled to unit norm).
# For the KS equation with n≈44800 training rows, the false-positive constant
# column '1' has col_norm≈212 and original-space coefficient≈0.066, giving
# |w_norm|≈14.  The D_TOLS grid must include values that let the search reach
# tol>14 within maxit=25 steps.  Large values (5.0, 15.0, 20.0) achieve this.
D_TOLS = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, 5.0, 15.0, 20.0]
BACKENDS = ["spectral", "polynomial"]

# SR3 relaxation-parameter grid (kappa); replaces LAMS when --algorithm sr3.
# Larger kappa tightens the ξ≈w coupling and improves conditioning of the
# normal equations independently of data scale.
SR3_KAPPAS = [0.01, 0.1, 1.0, 10.0, 100.0]

# ── equation definitions ─────────────────────────────────────────────────────
# Each entry:
#   name          : str
#   simulate_fn   : callable → (snapshots, times)
#   sim_kwargs    : dict
#   n_per_snap    : int
#   true_eq_fn    : callable(desc_eq) → ndarray of true coefficients
#   skip_fraction : float  (fraction of front snapshots to discard as transient)


def _find_op(desc_eq, name):
    """Return the column index of operator `name` in desc_eq.

    Raises KeyError if the name is not found.  Using this instead of
    hardcoded integer indices keeps true-coefficient vectors correct even if
    the library ordering changes (e.g. when max_p or max_n are varied).
    """
    try:
        return desc_eq.index(name)
    except ValueError:
        raise KeyError(
            f"Operator {name!r} not found in library.  "
            f"Available operators: {desc_eq}"
        )


def _true_allen_cahn(desc_eq, eps=0.15):
    """Allen–Cahn: f_t = eps²·Δf + f − f³."""
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, 'f')] = 1.0
    w[_find_op(desc_eq, 'f³')] = -1.0
    w[_find_op(desc_eq, 'Δf')] = eps ** 2
    return w


def _true_swift_hohenberg(desc_eq, r=0.3):
    """Swift–Hohenberg: f_t = (r−1)f − 2Δf − Δ²f − f³."""
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, 'f')] = r - 1.0
    w[_find_op(desc_eq, 'Δf')] = -2.0
    w[_find_op(desc_eq, 'Δ²f')] = -1.0
    w[_find_op(desc_eq, 'f³')] = -1.0
    return w


def _mask_pme(f):
    """Boolean mask: keep only spatial points where f > 0.05.
    This excludes the degenerate support boundary of the porous-medium PDE
    where the diffusivity D=f^2 → 0 and the PDE itself degenerates.
    """
    return f.ravel() > 0.05


def _true_pme(desc_eq):
    """Porous-medium equation: f_t = ∇·(f²∇f) = f²·Δf + 2f·‖∇f‖²."""
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, 'f²·Δf')] = 1.0
    w[_find_op(desc_eq, 'f·‖∇f‖²')] = 2.0
    return w


def _true_cahn_hilliard(desc_eq, eps=0.2):
    """Cahn–Hilliard: f_t = −Δf + 3f²·Δf + 6f·‖∇f‖² − eps²·Δ²f."""
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, 'Δf')] = -1.0
    w[_find_op(desc_eq, 'f²·Δf')] = 3.0
    w[_find_op(desc_eq, 'f·‖∇f‖²')] = 6.0
    w[_find_op(desc_eq, 'Δ²f')] = -(eps ** 2)
    return w


def _true_ks(desc_eq):
    """Kuramoto–Sivashinsky: f_t = −½‖∇f‖² − Δf − Δ²f."""
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, '‖∇f‖²')] = -0.5
    w[_find_op(desc_eq, 'Δf')] = -1.0
    w[_find_op(desc_eq, 'Δ²f')] = -1.0
    return w


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
        # L=8π (4 critical wavelengths per side, λ_c=2π): small enough to simulate
        # cheaply, large enough to contain rich defect texture throughout the sweep.
        # L=4π saturates to a single-stripe state by t≈20, collapsing library
        # columns regardless of T_end.
        # save_every=5 → dt_snap=0.5 ≈ 0.15*τ_c (τ_c=1/r≈3.3): fine enough for
        # accurate central-FD estimates of f_t.
        # skip_fraction=0.1: discards the low-amplitude initial transient where f³≈0.
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
        # Mask out the degenerate support boundary (f → 0) where D=f² vanishes.
        # Without this mask the regression includes near-zero-f points where
        # the PDE degenerates and coefficients are numerically meaningless.
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

# Equations run by default (KS excluded — too slow)
DEFAULT_EQUATIONS = ["allen_cahn", "swift_hohenberg", "porous_medium", "cahn_hilliard"]

# ── helpers ──────────────────────────────────────────────────────────────────


def build_matrices(snapshots, times, n_per_snap, L, skip_fraction, mask_fn=None, spectral_cutoff=None):
    """Simulate → regression data → both library matrices for both backends.

    Parameters
    ----------
    mask_fn : callable(f_snap) → bool array of shape (N*N,), or None
        If provided, applied to each snapshot's flattened field to select
        which spatial points enter the regression.  Replaces the random
        idx_sets sub-sampling for masked equations (masked points + random
        sub-sample of the remaining ones).
    spectral_cutoff : float or None
        If set, each snapshot is passed through an isotropic spectral low-pass
        filter (zeroing modes with |k| > spectral_cutoff) before library
        construction.  Does not affect the f_t target vector (computed from
        unfiltered snapshots via central FD on consecutive snaps).
    """
    f_snaps, ft_snaps, idx_sets, _ = make_regression_data(
        snapshots,
        times,
        n_per_snap=n_per_snap,
        seed=0,
    )
    # optionally skip transient front
    n = len(f_snaps)
    start = int(np.floor(skip_fraction * n))
    f_snaps = f_snaps[start:]
    ft_snaps = ft_snaps[start:]
    idx_sets = idx_sets[start:]

    # Apply per-snapshot spatial mask if provided
    if mask_fn is not None:
        rng = np.random.default_rng(1)
        new_idx_sets = []
        for s in range(len(f_snaps)):
            valid = np.where(mask_fn(f_snaps[s]))[0]  # indices passing mask
            if len(valid) <= n_per_snap:
                new_idx_sets.append(valid)
            else:
                chosen = rng.choice(valid, size=n_per_snap, replace=False)
                new_idx_sets.append(np.sort(chosen))
        idx_sets = new_idx_sets

    Ft = np.concatenate([ft_snaps[s].ravel()[idx_sets[s]] for s in range(len(f_snaps))])

    eq_spec_rows, std_spec_rows = [], []
    eq_poly_rows, std_poly_rows = [], []

    n_snaps = len(f_snaps)
    t_spec_total = t_poly_total = 0.0
    t_loop_start = time.time()
    print(f"  Building libraries for {n_snaps} snapshots...", flush=True)

    for s in tqdm(
        range(n_snaps), desc="  libraries", unit="snap", ncols=90, leave=True
    ):
        fs = f_snaps[s]
        if spectral_cutoff is not None:
            fs = spectral_lowpass(fs, L, L, spectral_cutoff)
        idx = idx_sets[s]

        t_s = time.time()
        with derivatives.use_backend("spectral"):
            Teq, desc_eq = build_equivariant_library(fs, L, L)
            Tstd, desc_std = build_standard_library(fs, L, L, max_order=4, max_poly=3)
        eq_spec_rows.append(Teq[idx])
        std_spec_rows.append(Tstd[idx])
        t_spec_total += time.time() - t_s

        t_p = time.time()
        with derivatives.use_backend("polynomial"):
            Teq, _ = build_equivariant_library(fs, L, L)
            Tstd, _ = build_standard_library(fs, L, L, max_order=4, max_poly=3)
        eq_poly_rows.append(Teq[idx])
        std_poly_rows.append(Tstd[idx])
        t_poly_total += time.time() - t_p

        # Print a timing line every 10 snapshots
        if (s + 1) % 10 == 0 or s == n_snaps - 1:
            elapsed = time.time() - t_loop_start
            rate = (s + 1) / elapsed
            eta = (n_snaps - s - 1) / rate if rate > 0 else float("inf")
            print(
                f"    snap {s+1}/{n_snaps}  "
                f"spectral {t_spec_total/(s+1)*1e3:.0f}ms/snap  "
                f"poly {t_poly_total/(s+1)*1e3:.0f}ms/snap  "
                f"ETA {eta:.0f}s",
                flush=True,
            )

    print(
        f"  Library build done in {time.time()-t_loop_start:.1f}s  "
        f"(spectral {t_spec_total:.1f}s  polynomial {t_poly_total:.1f}s)",
        flush=True,
    )

    mats = {
        ("spectral", "equivariant"): np.vstack(eq_spec_rows).astype(float),
        ("spectral", "standard"): np.vstack(std_spec_rows).astype(float),
        ("polynomial", "equivariant"): np.vstack(eq_poly_rows).astype(float),
        ("polynomial", "standard"): np.vstack(std_poly_rows).astype(float),
    }
    return mats, Ft, desc_eq, desc_std


def sweep_one(Theta, Ft, true_eq, lam, d_tol, algorithm="strridge"):
    """Run TrainSTRidge or TrainSR3 and return (n_active, fp, fn, exact_recovery, max_coeff_err)."""
    act = np.abs(true_eq) > 1e-10
    if algorithm == "sr3":
        w = TrainSR3(Theta, Ft, lam, d_tol)
    else:
        w = TrainSTRidge(Theta, Ft, lam, d_tol)
    nz = int(np.count_nonzero(np.abs(w) > 1e-10))
    fp = int(np.sum((np.abs(w) > 1e-10) & ~act))
    fn = int(np.sum((np.abs(w) < 1e-10) & act))
    exact = int(fp == 0 and fn == 0)
    # coefficient error only on true support (if any true terms were found)
    found = act & (np.abs(w) > 1e-10)
    if found.any():
        err = float(np.max(np.abs(w[act] - true_eq[act])))
    else:
        err = float("nan")
    return nz, fp, fn, exact, err


# ── main ─────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(description="Hyperparameter sweep for STRidge.")
    parser.add_argument(
        "--equations",
        nargs="+",
        default=DEFAULT_EQUATIONS,
        help="Which equations to run.",
    )
    parser.add_argument(
        "--name",
        default="",
        help="Optional distinctive tag appended to the output filename.",
    )
    parser.add_argument(
        "--noise",
        type=float,
        default=0.0,
        help="Observation noise level as a fraction of signal std (e.g. 0.02 = 2%%). "
             "Passed to the simulator as obs_noise.",
    )
    parser.add_argument(
        "--spectral-cutoff",
        type=float,
        default=None,
        dest="spectral_cutoff",
        help="Isotropic spectral low-pass cutoff (physical wavenumber). "
             "If set, each snapshot is filtered before library construction. "
             "Good starting point for Swift-Hohenberg (k_c=1): 3.0.",
    )
    parser.add_argument(
        "--n-per-snap",
        type=int,
        default=None,
        dest="n_per_snap",
        help="Override n_per_snap (spatial sub-samples per snapshot) for all equations.",
    )
    parser.add_argument(
        "--algorithm",
        choices=["strridge", "sr3"],
        default="strridge",
        help="Regression algorithm: strridge (default) or sr3.",
    )
    args = parser.parse_args()

    eq_names_to_run = set(args.equations)
    equations_to_run = [eq for eq in EQUATIONS if eq["name"] in eq_names_to_run]
    if not equations_to_run:
        print(f'No matching equations found. Available: {[e["name"] for e in EQUATIONS]}')
        return

    obs_noise = args.noise
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"_{args.name}" if args.name else ""
    out_path = os.path.join(RESULTS_DIR, f"hyperparam_sweep_{timestamp}{tag}.csv")
    print(f"Output → {out_path}")

    fieldnames = [
        "equation",
        "backend",
        "library",
        "algorithm",
        "lam",
        "d_tol",
        "noise",
        "n_active",
        "fp",
        "fn",
        "exact_recovery",
        "max_coeff_err",
    ]

    # Always write a fresh file (timestamped filename guarantees no collisions).
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for eq in equations_to_run:
            name = eq["name"]
            L = eq["sim_kwargs"].get("L", 2 * np.pi)
            print(f'\n{"="*60}')
            print(f"  Equation: {name}")
            print(f'{"="*60}')

            # ── simulate ────────────────────────────────────────────────
            t0 = time.time()
            print(f"  Simulating... (obs_noise={obs_noise})", flush=True)
            snapshots, times = eq["simulate_fn"](**eq["sim_kwargs"], obs_noise=obs_noise)
            print(f"  done in {time.time()-t0:.1f}s  " f"({len(snapshots)} snapshots)", flush=True)

            # ── build library matrices (once per equation) ───────────────
            n_per_snap = args.n_per_snap if args.n_per_snap is not None else eq["n_per_snap"]
            cutoff_info = f", spectral_cutoff={args.spectral_cutoff}" if args.spectral_cutoff else ""
            t0 = time.time()
            print(f"  Building library matrices (n_per_snap={n_per_snap}{cutoff_info})...", flush=True)
            mats, Ft, desc_eq, desc_std = build_matrices(
                snapshots,
                times,
                n_per_snap,
                L,
                eq["skip_fraction"],
                mask_fn=eq.get("mask_fn"),
                spectral_cutoff=args.spectral_cutoff,
            )
            print(f"  done in {time.time()-t0:.1f}s  " f"({len(Ft)} regression points)", flush=True)

            true_eq = eq["true_eq_fn"](desc_eq)
            true_std = np.zeros(len(desc_std))  # std library: no clean true_eq mapping

            # condition numbers
            for (backend, lib), Theta in mats.items():
                kappa = np.linalg.cond(Theta)
                print(f"  cond({backend:10s} {lib:12s}): {kappa:.2e}")

            # ── hyperparameter sweep ─────────────────────────────────────
            lam_grid = SR3_KAPPAS if args.algorithm == "sr3" else LAMS
            lam_label = "kappa" if args.algorithm == "sr3" else "lam"
            combos = list(itertools.product(BACKENDS, lam_grid, D_TOLS))
            t0 = time.time()
            n_combos = len(combos)
            print(f"  Running {n_combos} (backend, {lam_label}, d_tol) combos "
                  f"[algorithm={args.algorithm}]...", flush=True)

            pbar = tqdm(
                combos,
                desc=f"  {name}",
                unit="combo",
                ncols=90,
                leave=True,
            )
            t_combo_start = time.time()
            for i, (backend, lam, d_tol) in enumerate(pbar):
                pbar.set_postfix(backend=backend[:4], lam=f"{lam:.0e}", d_tol=f"{d_tol:.0e}")
                if i > 0 and i % 20 == 0:
                    avg = (time.time() - t_combo_start) / i
                    eta = avg * (n_combos - i)
                    print(f"    {i}/{n_combos} combos done  avg {avg*1e3:.0f}ms/combo  ETA {eta:.0f}s", flush=True)
                Theta_eq  = mats[(backend, "equivariant")]
                Theta_std = mats[(backend, "standard")]

                # equivariant library
                nz, fp, fn, exact, err = sweep_one(Theta_eq, Ft, true_eq, lam, d_tol, algorithm=args.algorithm)
                writer.writerow(
                    dict(
                        equation=name,
                        backend=backend,
                        library="equivariant",
                        algorithm=args.algorithm,
                        lam=lam,
                        d_tol=d_tol,
                        noise=obs_noise,
                        n_active=nz,
                        fp=fp,
                        fn=fn,
                        exact_recovery=exact,
                        max_coeff_err=err,
                    )
                )

                # standard library — can only count active terms
                nz_s, fp_s, fn_s, exact_s, err_s = sweep_one(Theta_std, Ft, true_std, lam, d_tol, algorithm=args.algorithm)
                writer.writerow(
                    dict(
                        equation=name,
                        backend=backend,
                        library="standard",
                        algorithm=args.algorithm,
                        lam=lam,
                        d_tol=d_tol,
                        noise=obs_noise,
                        n_active=nz_s,
                        fp=float("nan"),
                        fn=float("nan"),
                        exact_recovery=float("nan"),
                        max_coeff_err=float("nan"),
                    )
                )

            fh.flush()
            print(f"  Sweep done in {time.time()-t0:.1f}s", flush=True)

            # ── print best result per (backend, library) ──────────────────────────
            import csv as _csv_mod
            with open(out_path, newline="") as rfh2:
                all_rows = [r for r in _csv_mod.DictReader(rfh2) if r["equation"] == name]
            print(f"  Best equivariant results ({name}):")
            for backend in BACKENDS:
                eq_rows = [
                    r for r in all_rows
                    if r["backend"] == backend and r["library"] == "equivariant"
                    and r["exact_recovery"] not in ("", "nan")
                ]
                if not eq_rows:
                    continue
                exact = [r for r in eq_rows if r["exact_recovery"] == "1"]
                best = min(eq_rows, key=lambda r: (int(r["fn"]), int(r["fp"]), float(r["max_coeff_err"]) if r["max_coeff_err"] not in ("", "nan") else 1e9))
                mark = " ✓ EXACT" if exact else ""
                print(
                    f"    {backend:11s}: best fn={best['fn']}, fp={best['fp']}, "
                    f"err={float(best['max_coeff_err']):.3e}"
                    f" at lam={best['lam']}, d_tol={best['d_tol']}{mark}"
                )

    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
