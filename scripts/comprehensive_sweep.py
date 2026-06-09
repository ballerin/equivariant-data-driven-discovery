"""
Comprehensive hyperparameter sweep for the paper.

Runs all combinations of:
  - 5 equations  (Allen-Cahn, Swift-Hohenberg, Porous Medium, KS,
                  Cahn-Hilliard with eps ∈ {0.1, 0.2, 0.4})
  - Noise levels  : 0%, 2%, 5% observation noise
  - Preprocessing : none, spectral low-pass (cutoff = 3.0)
  - Algorithms    : STRidge, SR3
  - Backends      : spectral, polynomial  (inside matrix build)
  - Libraries     : equivariant, standard (inside matrix build)
  - D_TOLS        : 11-value grid (same as hyperparam_sweep.py)

Efficiency strategy
-------------------
* Simulate once per equation (clean, obs_noise=0).  Noise is applied
  post-hoc by adding sigma*randn to the clean snapshots.  This avoids
  repeating expensive PDE solvers for each noise level.
* Build library matrices once per (equation, noise, spectral_cutoff).
  A second matrix build per noise level handles the lowpass variant.
* Run both STRidge and SR3 regression sweeps on the same matrices.

Output
------
results/comprehensive_sweep_<timestamp>[_<name>].csv
with columns:
    equation, backend, library, algorithm, lam, d_tol,
    noise, spectral_cutoff, n_active, fp, fn, exact_recovery, max_coeff_err

Usage
-----
    cd paper-experimental/code
    python scripts/comprehensive_sweep.py
    python scripts/comprehensive_sweep.py --equations allen_cahn swift_hohenberg
    python scripts/comprehensive_sweep.py --no-ks          # skip slow KS
    python scripts/comprehensive_sweep.py --name test
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
        desc = kwargs.get("desc", "")
        items = list(iterable)
        print(f"{desc} ({len(items)} items)...")
        return items

# ── path setup ───────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, CODE_DIR)

import derivatives                                                     # noqa: E402
from library import (                                                  # noqa: E402
    build_equivariant_library,
    build_standard_library,
    spectral_lowpass,
)
from regression import TrainSTRidge, TrainSR3                         # noqa: E402
from simulations import (                                             # noqa: E402
    allen_cahn,
    cahn_hilliard,
    kuramoto_sivashinsky,
    make_regression_data,
    porous_medium,
    swift_hohenberg,
)

RESULTS_DIR = os.path.join(CODE_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyperparameter grids ─────────────────────────────────────────────────────
NOISE_LEVELS      = [0.0, 0.02, 0.05]
SPECTRAL_CUTOFFS  = [None, 3.0]          # None = no preprocessing
ALGORITHMS        = ["strridge", "sr3"]
BACKENDS          = ["spectral", "polynomial"]
LAMS              = [1e-5]               # STRidge (lam unused in inner loop)
SR3_KAPPAS        = [0.01, 0.1, 1.0, 10.0, 100.0]
D_TOLS            = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, 5.0, 15.0, 20.0]

# Cahn-Hilliard is run at three interface-width values to study how the
# small biharmonic coefficient eps² affects recovery as eps shrinks.
CH_EPSILONS = [0.1, 0.2, 0.4]

# ── operator-name lookup ──────────────────────────────────────────────────────
def _find_op(desc_eq, name):
    try:
        return desc_eq.index(name)
    except ValueError:
        raise KeyError(
            f"Operator {name!r} not found in library. "
            f"Available: {desc_eq}"
        )


# ── true-coefficient vectors ──────────────────────────────────────────────────
def _true_allen_cahn(desc_eq, eps=0.15):
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, 'f')] = 1.0
    w[_find_op(desc_eq, 'f³')] = -1.0
    w[_find_op(desc_eq, 'Δf')] = eps ** 2
    return w


def _true_swift_hohenberg(desc_eq, r=0.3):
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, 'f')] = r - 1.0
    w[_find_op(desc_eq, 'Δf')] = -2.0
    w[_find_op(desc_eq, 'Δ²f')] = -1.0
    w[_find_op(desc_eq, 'f³')] = -1.0
    return w


def _mask_pme(f):
    """Keep spatial points where f > 0.05 (away from the degenerate support boundary)."""
    return f.ravel() > 0.05


def _true_pme(desc_eq):
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, 'f²·Δf')] = 1.0
    w[_find_op(desc_eq, 'f·‖∇f‖²')] = 2.0
    return w


def _true_cahn_hilliard(desc_eq, eps=0.2):
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, 'Δf')] = -1.0
    w[_find_op(desc_eq, 'f²·Δf')] = 3.0
    w[_find_op(desc_eq, 'f·‖∇f‖²')] = 6.0
    w[_find_op(desc_eq, 'Δ²f')] = -(eps ** 2)
    return w


def _true_ks(desc_eq):
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, '‖∇f‖²')] = -0.5
    w[_find_op(desc_eq, 'Δf')] = -1.0
    w[_find_op(desc_eq, 'Δ²f')] = -1.0
    return w


# ── equation specifications ───────────────────────────────────────────────────
# Each dict has:
#   name          : str  (written to CSV; must be unique for CH epsilon variants)
#   simulate_fn   : callable(**sim_kwargs, obs_noise=0) → (snapshots, times)
#   sim_kwargs    : dict  (passed without obs_noise; noise applied post-hoc)
#   n_per_snap    : int
#   true_eq_fn    : callable(desc_eq) → ndarray
#   skip_fraction : float
#   mask_fn       : callable(f_snap) → bool array, or None
#   L             : float  (domain side-length, derived from sim_kwargs)

def _make_equation_specs():
    specs = []

    # Allen-Cahn
    specs.append(dict(
        name="allen_cahn",
        simulate_fn=allen_cahn,
        sim_kwargs=dict(N=64, L=2*np.pi, T_end=8.0, eps=0.15, dt=0.01, save_every=20, seed=42),
        n_per_snap=200,
        true_eq_fn=lambda d: _true_allen_cahn(d, eps=0.15),
        skip_fraction=0.0,
        mask_fn=None,
    ))

    # Swift-Hohenberg
    specs.append(dict(
        name="swift_hohenberg",
        simulate_fn=swift_hohenberg,
        sim_kwargs=dict(N=64, L=8*np.pi, T_end=50.0, r=0.3, dt=0.1, save_every=5, seed=42),
        n_per_snap=200,
        true_eq_fn=lambda d: _true_swift_hohenberg(d, r=0.3),
        skip_fraction=0.1,
        mask_fn=None,
    ))

    # Porous Medium
    specs.append(dict(
        name="porous_medium",
        simulate_fn=porous_medium,
        sim_kwargs=dict(N=64, L=6*np.pi, T_end=2.0, dt=0.002, save_every=25, seed=42),
        n_per_snap=200,
        true_eq_fn=_true_pme,
        skip_fraction=0.0,
        mask_fn=_mask_pme,
    ))

    # Cahn-Hilliard — three epsilon values
    for eps in CH_EPSILONS:
        specs.append(dict(
            name=f"cahn_hilliard_eps{eps:.2g}",
            simulate_fn=cahn_hilliard,
            sim_kwargs=dict(N=128, L=4*np.pi, T_end=30.0, eps=eps, dt=0.005, save_every=30, seed=42),
            n_per_snap=200,
            true_eq_fn=lambda d, e=eps: _true_cahn_hilliard(d, eps=e),
            skip_fraction=0.0,
            mask_fn=None,
        ))

    # Kuramoto-Sivashinsky
    specs.append(dict(
        name="kuramoto_sivashinsky",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=64, L=16*np.pi, T_end=200.0, dt=0.1, save_every=5, seed=42),
        n_per_snap=200,
        true_eq_fn=_true_ks,
        skip_fraction=0.3,
        mask_fn=None,
    ))

    # derive L for each spec
    for s in specs:
        s['L'] = s['sim_kwargs']['L']

    return specs


EQUATION_SPECS = _make_equation_specs()

ALL_NAMES = [s['name'] for s in EQUATION_SPECS]
# Default excludes KS (slow) unless explicitly requested
DEFAULT_NAMES = [n for n in ALL_NAMES if not n.startswith("kuramoto")]


# ── noise application ─────────────────────────────────────────────────────────
def _apply_noise(clean_snapshots, noise_level, seed):
    """Return clean_snapshots + noise_level * std(clean) * randn.

    Uses a fixed seed per (noise_level) so results are reproducible.
    The seed is hashed from the noise level to avoid correlations between
    noise realisations at different levels.
    """
    if noise_level == 0.0:
        return clean_snapshots
    rng = np.random.default_rng(int(abs(noise_level) * 1e6) + seed)
    sigma = noise_level * clean_snapshots.std()
    return clean_snapshots + sigma * rng.standard_normal(clean_snapshots.shape)


# ── library matrix builder ────────────────────────────────────────────────────
def build_matrices(snapshots, times, n_per_snap, L, skip_fraction,
                   mask_fn=None, spectral_cutoff=None):
    """Build all four (backend, library) design matrices.

    Returns
    -------
    mats     : dict {(backend, library): Theta}
    Ft       : target vector (f_t)
    desc_eq  : operator names for equivariant library
    desc_std : operator names for standard library
    """
    f_snaps, ft_snaps, idx_sets, _ = make_regression_data(
        snapshots, times, n_per_snap=n_per_snap, seed=0
    )
    n = len(f_snaps)
    start = int(np.floor(skip_fraction * n))
    f_snaps  = f_snaps[start:]
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
                chosen = rng.choice(valid, size=n_per_snap, replace=False)
                new_idx.append(np.sort(chosen))
        idx_sets = new_idx

    Ft = np.concatenate([ft_snaps[s].ravel()[idx_sets[s]] for s in range(len(f_snaps))])

    eq_spec, std_spec, eq_poly, std_poly = [], [], [], []
    desc_eq = desc_std = None

    n_snaps = len(f_snaps)
    t0 = time.time()
    print(f"  Building libraries ({n_snaps} snaps"
          + (f", cutoff={spectral_cutoff}" if spectral_cutoff else "") + ")...",
          flush=True)

    for s in tqdm(range(n_snaps), desc="  libs", unit="snap", ncols=90, leave=True):
        fs = f_snaps[s]
        if spectral_cutoff is not None:
            fs = spectral_lowpass(fs, L, L, spectral_cutoff)
        idx = idx_sets[s]

        with derivatives.use_backend("spectral"):
            Teq, desc_eq   = build_equivariant_library(fs, L, L)
            Tstd, desc_std = build_standard_library(fs, L, L, max_order=4, max_poly=3)
        eq_spec.append(Teq[idx])
        std_spec.append(Tstd[idx])

        with derivatives.use_backend("polynomial"):
            Teq, _ = build_equivariant_library(fs, L, L)
            Tstd, _ = build_standard_library(fs, L, L, max_order=4, max_poly=3)
        eq_poly.append(Teq[idx])
        std_poly.append(Tstd[idx])

    print(f"  Library build done in {time.time()-t0:.1f}s", flush=True)

    mats = {
        ("spectral",   "equivariant"): np.vstack(eq_spec).astype(float),
        ("spectral",   "standard"):    np.vstack(std_spec).astype(float),
        ("polynomial", "equivariant"): np.vstack(eq_poly).astype(float),
        ("polynomial", "standard"):    np.vstack(std_poly).astype(float),
    }
    return mats, Ft, desc_eq, desc_std


# ── single regression run ─────────────────────────────────────────────────────
def sweep_one(Theta, Ft, true_eq, lam, d_tol, algorithm):
    """Run STRidge or SR3 and return (n_active, fp, fn, exact_recovery, max_coeff_err)."""
    act = np.abs(true_eq) > 1e-10
    if algorithm == "sr3":
        w = TrainSR3(Theta, Ft, lam, d_tol)
    else:
        w = TrainSTRidge(Theta, Ft, lam, d_tol)
    nz = int(np.count_nonzero(np.abs(w) > 1e-10))
    fp = int(np.sum((np.abs(w) > 1e-10) & ~act))
    fn = int(np.sum((np.abs(w) < 1e-10) & act))
    exact = int(fp == 0 and fn == 0)
    found = act & (np.abs(w) > 1e-10)
    err = float(np.max(np.abs(w[act] - true_eq[act]))) if found.any() else float("nan")
    return nz, fp, fn, exact, err


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive sweep: all equations × noise × preprocessing × algorithms."
    )
    parser.add_argument(
        "--equations", nargs="+", default=None,
        help=(
            "Names of equations to run.  Defaults to all except KS.  "
            f"Available: {ALL_NAMES}"
        ),
    )
    parser.add_argument(
        "--no-ks", action="store_true",
        help="Exclude Kuramoto-Sivashinsky (slow; excluded from default).",
    )
    parser.add_argument(
        "--include-ks", action="store_true",
        help="Include Kuramoto-Sivashinsky (excluded by default).",
    )
    parser.add_argument(
        "--name", default="",
        help="Optional tag appended to the output filename.",
    )
    parser.add_argument(
        "--noise-levels", nargs="+", type=float, default=NOISE_LEVELS,
        dest="noise_levels",
        help="Observation noise levels to sweep (default: 0.0 0.02 0.05).",
    )
    parser.add_argument(
        "--algorithms", nargs="+", default=ALGORITHMS,
        choices=["strridge", "sr3"],
        help="Algorithms to run (default: both strridge and sr3).",
    )
    parser.add_argument(
        "--no-lowpass", action="store_true",
        help="Skip the spectral low-pass preprocessing variant.",
    )
    args = parser.parse_args()

    # Resolve equation list
    if args.equations:
        names_to_run = set(args.equations)
    else:
        names_to_run = set(DEFAULT_NAMES)
        if args.include_ks:
            names_to_run.add("kuramoto_sivashinsky")
    if args.no_ks:
        names_to_run.discard("kuramoto_sivashinsky")

    specs_to_run = [s for s in EQUATION_SPECS if s["name"] in names_to_run]
    if not specs_to_run:
        print(f"No matching equations. Available: {ALL_NAMES}")
        return

    cutoffs = [None] if args.no_lowpass else SPECTRAL_CUTOFFS

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"_{args.name}" if args.name else ""
    out_path = os.path.join(RESULTS_DIR, f"comprehensive_sweep_{timestamp}{tag}.csv")
    print(f"Output → {out_path}")
    print(f"Equations : {[s['name'] for s in specs_to_run]}")
    print(f"Noise     : {args.noise_levels}")
    print(f"Cutoffs   : {cutoffs}")
    print(f"Algorithms: {args.algorithms}")

    fieldnames = [
        "equation", "backend", "library", "algorithm",
        "lam", "d_tol", "noise", "spectral_cutoff",
        "n_active", "fp", "fn", "exact_recovery", "max_coeff_err",
    ]

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for spec in specs_to_run:
            name  = spec["name"]
            L     = spec["L"]
            print(f'\n{"="*70}')
            print(f"  Equation: {name}")
            print(f'{"="*70}')

            # ── simulate once (clean) ────────────────────────────────────
            t0 = time.time()
            print("  Simulating (clean)...", flush=True)
            clean_snaps, times = spec["simulate_fn"](**spec["sim_kwargs"], obs_noise=0.0)
            print(f"  done in {time.time()-t0:.1f}s  ({len(clean_snaps)} snapshots)", flush=True)

            # ── loop over noise levels ───────────────────────────────────
            for noise in args.noise_levels:
                print(f'\n  {"─"*60}')
                print(f"  noise = {noise:.0%}", flush=True)

                noisy_snaps = _apply_noise(clean_snaps, noise, seed=7919)

                # ── loop over preprocessing ──────────────────────────────
                for cutoff in cutoffs:
                    cutoff_str = f"{cutoff}" if cutoff is not None else "none"
                    print(f'\n  spectral_cutoff = {cutoff_str}', flush=True)

                    mats, Ft, desc_eq, desc_std = build_matrices(
                        noisy_snaps, times,
                        n_per_snap=spec["n_per_snap"],
                        L=L,
                        skip_fraction=spec["skip_fraction"],
                        mask_fn=spec.get("mask_fn"),
                        spectral_cutoff=cutoff,
                    )

                    # Condition numbers
                    for (backend, lib), Theta in mats.items():
                        kappa = np.linalg.cond(Theta)
                        print(f"    cond({backend:10s} {lib:12s}): {kappa:.2e}")

                    true_eq  = spec["true_eq_fn"](desc_eq)
                    true_std = np.zeros(len(desc_std))

                    # ── loop over algorithms ─────────────────────────────
                    for algorithm in args.algorithms:
                        lam_grid = SR3_KAPPAS if algorithm == "sr3" else LAMS
                        combos = list(itertools.product(BACKENDS, lam_grid, D_TOLS))
                        print(f"\n    algorithm={algorithm}  ({len(combos)} combos)...", flush=True)

                        pbar = tqdm(combos, desc=f"    {name[:20]}", unit="combo",
                                    ncols=90, leave=True)
                        t_reg = time.time()
                        for backend, lam, d_tol in pbar:
                            Theta_eq  = mats[(backend, "equivariant")]
                            Theta_std = mats[(backend, "standard")]

                            # equivariant
                            nz, fp, fn, exact, err = sweep_one(
                                Theta_eq, Ft, true_eq, lam, d_tol, algorithm
                            )
                            writer.writerow(dict(
                                equation=name, backend=backend,
                                library="equivariant", algorithm=algorithm,
                                lam=lam, d_tol=d_tol,
                                noise=noise, spectral_cutoff=cutoff_str,
                                n_active=nz, fp=fp, fn=fn,
                                exact_recovery=exact, max_coeff_err=err,
                            ))

                            # standard
                            nz_s, fp_s, fn_s, exact_s, err_s = sweep_one(
                                Theta_std, Ft, true_std, lam, d_tol, algorithm
                            )
                            writer.writerow(dict(
                                equation=name, backend=backend,
                                library="standard", algorithm=algorithm,
                                lam=lam, d_tol=d_tol,
                                noise=noise, spectral_cutoff=cutoff_str,
                                n_active=nz_s, fp=float("nan"), fn=float("nan"),
                                exact_recovery=float("nan"),
                                max_coeff_err=float("nan"),
                            ))

                        print(f"    algorithm={algorithm} done in {time.time()-t_reg:.1f}s",
                              flush=True)

                    fh.flush()

            # ── print summary for this equation ──────────────────────────
            print(f'\n  Summary for {name}:')
            import csv as _csv_mod
            with open(out_path, newline="") as rfh:
                all_rows = [r for r in _csv_mod.DictReader(rfh) if r["equation"] == name]
            exact_rows = [r for r in all_rows if r.get("exact_recovery") == "1"]
            if exact_rows:
                for r in exact_rows:
                    print(
                        f"    ✓ EXACT: backend={r['backend']}, lib={r['library']}, "
                        f"algo={r['algorithm']}, noise={r['noise']}, "
                        f"cutoff={r['spectral_cutoff']}, "
                        f"d_tol={r['d_tol']}, err={float(r['max_coeff_err']):.3e}"
                    )
            else:
                eq_rows = [r for r in all_rows if r["library"] == "equivariant"
                           and r["exact_recovery"] not in ("", "nan")]
                if eq_rows:
                    best = min(eq_rows, key=lambda r: (
                        int(r["fn"]), int(r["fp"]),
                        float(r["max_coeff_err"]) if r["max_coeff_err"] not in ("", "nan") else 1e9
                    ))
                    print(
                        f"    Best (no exact): fn={best['fn']}, fp={best['fp']}, "
                        f"err={best['max_coeff_err']}, "
                        f"backend={best['backend']}, algo={best['algorithm']}, "
                        f"noise={best['noise']}, cutoff={best['spectral_cutoff']}"
                    )
                else:
                    print("    No equivariant results.")

    print(f"\nDone. Results saved to {out_path}")


if __name__ == "__main__":
    main()
