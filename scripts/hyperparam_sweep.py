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

import sys
import os
import argparse
import itertools
import time
import csv

try:
    from tqdm import tqdm
except ImportError:
    # Graceful fallback: tqdm not installed
    def tqdm(iterable, **kwargs):
        desc = kwargs.get('desc', '')
        items = list(iterable)
        print(f'{desc} ({len(items)} items)...')
        return items

import numpy as np

# ── path setup ──────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR   = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, CODE_DIR)

import derivatives
from simulations import (
    allen_cahn, swift_hohenberg, cahn_hilliard,
    kuramoto_sivashinsky, porous_medium, make_regression_data,
)
from library import build_equivariant_library, build_standard_library
from regression import TrainSTRidge

RESULTS_DIR = os.path.join(CODE_DIR, 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyperparameter grid ──────────────────────────────────────────────────────
LAMS   = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
D_TOLS = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2]
BACKENDS = ['spectral', 'polynomial']

# ── equation definitions ─────────────────────────────────────────────────────
# Each entry:
#   name          : str
#   simulate_fn   : callable → (snapshots, times)
#   sim_kwargs    : dict
#   n_per_snap    : int
#   true_eq_fn    : callable(desc_eq) → ndarray of true coefficients
#   skip_fraction : float  (fraction of front snapshots to discard as transient)

def _true_allen_cahn(desc_eq):
    w = np.zeros(20)
    w[4]  =  0.15**2   # eps^2 * Delta f
    w[1]  =  1.0       # f
    w[3]  = -1.0       # -f^3
    return w

def _true_swift_hohenberg(desc_eq):
    w = np.zeros(20)
    w[1]  = 0.3 - 1.0  # (r-1)*f,  r=0.3  → -0.7
    w[4]  = -2.0        # -2 Δf
    w[10] = -1.0        # -Δ²f
    w[3]  = -1.0        # -f³
    return w


def _mask_pme(f):
    """Boolean mask: keep only spatial points where f > 0.05.
    This excludes the degenerate support boundary of the porous-medium PDE
    where the diffusivity D=f^2 → 0 and the PDE itself degenerates.
    """
    return f.ravel() > 0.05

def _true_pme(desc_eq):
    w = np.zeros(20)
    w[6] =  1.0   # f^2 Delta f
    w[8] =  2.0   # f |grad f|^2   (note: N8 = f * ||grad f||^2)
    return w

def _true_cahn_hilliard(desc_eq, eps=0.2):
    w = np.zeros(20)
    w[4]  = -1.0        # -Delta f
    w[6]  =  3.0        # 3 f^2 Delta f
    w[8]  =  6.0        # 6 f ||grad f||^2
    w[10] = -(eps**2)   # -eps^2 Delta^2 f
    return w

def _true_ks(desc_eq):
    w = np.zeros(20)
    w[7]  = -0.5   # -1/2 ||grad f||^2
    w[4]  = -1.0   # -Delta f
    w[10] = -1.0   # -Delta^2 f
    return w

EQUATIONS = [
    dict(
        name         = 'allen_cahn',
        simulate_fn  = allen_cahn,
        sim_kwargs   = dict(N=64, L=2*np.pi, T_end=8.0, eps=0.15, dt=0.01,
                            save_every=20, seed=42),
        n_per_snap   = 200,
        true_eq_fn   = _true_allen_cahn,
        skip_fraction= 0.0,
        mask_fn      = None,
    ),
    dict(
        name         = 'swift_hohenberg',
        simulate_fn  = swift_hohenberg,
        # T_end=50 keeps us in the active coarsening/transient phase.
        # At T_end=200 the solution is a near-pure single-wavenumber stripe
        # whose library columns (f, Δf, Δ²f, f³) become collinear, making
        # the regression ill-posed regardless of hyperparameters.
        sim_kwargs   = dict(N=64, L=4*np.pi, T_end=50.0, r=0.3, dt=0.1,
                            save_every=10, seed=42),
        n_per_snap   = 200,
        true_eq_fn   = _true_swift_hohenberg,
        skip_fraction= 0.0,
        mask_fn      = None,
    ),
    dict(
        name         = 'porous_medium',
        simulate_fn  = porous_medium,
        sim_kwargs   = dict(N=64, L=6*np.pi, T_end=2.0, dt=0.002,
                            save_every=25, seed=42),
        n_per_snap   = 200,
        true_eq_fn   = _true_pme,
        skip_fraction= 0.0,
        # Mask out the degenerate support boundary (f → 0) where D=f² vanishes.
        # Without this mask the regression includes near-zero-f points where
        # the PDE degenerates and coefficients are numerically meaningless.
        mask_fn      = _mask_pme,
    ),
    dict(
        name         = 'cahn_hilliard',
        simulate_fn  = cahn_hilliard,
        sim_kwargs   = dict(N=128, L=4*np.pi, T_end=30.0, eps=0.2, dt=0.005,
                            save_every=30, seed=42),
        n_per_snap   = 200,
        true_eq_fn   = lambda d: _true_cahn_hilliard(d, eps=0.2),
        skip_fraction= 0.0,
        mask_fn      = None,
    ),
    dict(
        name         = 'kuramoto_sivashinsky',
        simulate_fn  = kuramoto_sivashinsky,
        sim_kwargs   = dict(N=64, L=16*np.pi, T_end=400.0, dt=0.05,
                            save_every=5, seed=42),
        n_per_snap   = 200,
        true_eq_fn   = _true_ks,
        skip_fraction= 0.4,
        mask_fn      = None,
    ),
]

# Equations run by default (KS excluded — too slow)
DEFAULT_EQUATIONS = ['allen_cahn', 'swift_hohenberg', 'porous_medium', 'cahn_hilliard']

# ── helpers ──────────────────────────────────────────────────────────────────

def build_matrices(snapshots, times, n_per_snap, L, skip_fraction, mask_fn=None):
    """Simulate → regression data → both library matrices for both backends.

    Parameters
    ----------
    mask_fn : callable(f_snap) → bool array of shape (N*N,), or None
        If provided, applied to each snapshot's flattened field to select
        which spatial points enter the regression.  Replaces the random
        idx_sets sub-sampling for masked equations (masked points + random
        sub-sample of the remaining ones).
    """
    f_snaps, ft_snaps, idx_sets, _ = make_regression_data(
        snapshots, times, n_per_snap=n_per_snap, seed=0,
    )
    # optionally skip transient front
    n = len(f_snaps)
    start = int(np.floor(skip_fraction * n))
    f_snaps   = f_snaps[start:]
    ft_snaps  = ft_snaps[start:]
    idx_sets  = idx_sets[start:]

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

    Ft = np.concatenate([ft_snaps[s].ravel()[idx_sets[s]]
                         for s in range(len(f_snaps))])

    eq_spec_rows, std_spec_rows = [], []
    eq_poly_rows, std_poly_rows = [], []

    for s in tqdm(range(len(f_snaps)), desc='  building libraries', unit='snap',
                  ncols=90, leave=False):
        fs  = f_snaps[s]
        idx = idx_sets[s]
        with derivatives.use_backend('spectral'):
            Teq,  desc_eq  = build_equivariant_library(fs, L, L)
            Tstd, desc_std = build_standard_library(fs, L, L,
                                                    max_order=4, max_poly=3)
        eq_spec_rows.append(Teq[idx]);  std_spec_rows.append(Tstd[idx])
        with derivatives.use_backend('polynomial'):
            Teq,  _ = build_equivariant_library(fs, L, L)
            Tstd, _ = build_standard_library(fs, L, L, max_order=4, max_poly=3)
        eq_poly_rows.append(Teq[idx]);  std_poly_rows.append(Tstd[idx])

    mats = {
        ('spectral',    'equivariant'): np.vstack(eq_spec_rows).astype(float),
        ('spectral',    'standard'):    np.vstack(std_spec_rows).astype(float),
        ('polynomial',  'equivariant'): np.vstack(eq_poly_rows).astype(float),
        ('polynomial',  'standard'):    np.vstack(std_poly_rows).astype(float),
    }
    return mats, Ft, desc_eq, desc_std


def sweep_one(Theta, Ft, true_eq, lam, d_tol):
    """Run TrainSTRidge and return (n_active, fp, fn, exact_recovery, max_coeff_err)."""
    act = np.abs(true_eq) > 1e-10
    w   = TrainSTRidge(Theta, Ft, lam, d_tol)
    nz  = int(np.count_nonzero(np.abs(w) > 1e-10))
    fp  = int(np.sum((np.abs(w) > 1e-10) & ~act))
    fn  = int(np.sum((np.abs(w) < 1e-10) &  act))
    exact = int(fp == 0 and fn == 0)
    # coefficient error only on true support (if any true terms were found)
    found = act & (np.abs(w) > 1e-10)
    if found.any():
        err = float(np.max(np.abs(w[act] - true_eq[act])))
    else:
        err = float('nan')
    return nz, fp, fn, exact, err


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Hyperparameter sweep for STRidge.')
    parser.add_argument(
        '--equations', nargs='+', default=DEFAULT_EQUATIONS,
        help='Which equations to run. Appends to existing CSV if it exists.'
    )
    args = parser.parse_args()

    eq_names_to_run = set(args.equations)
    equations_to_run = [eq for eq in EQUATIONS if eq['name'] in eq_names_to_run]
    if not equations_to_run:
        print(f'No matching equations found. Available: {[e["name"] for e in EQUATIONS]}')
        return

    out_path = os.path.join(RESULTS_DIR, 'hyperparam_sweep.csv')
    fieldnames = [
        'equation', 'backend', 'library', 'lam', 'd_tol',
        'n_active', 'fp', 'fn', 'exact_recovery', 'max_coeff_err',
    ]

    # Append mode if file exists (don't rewrite already-computed equations)
    file_exists = os.path.isfile(out_path)
    mode = 'a' if file_exists else 'w'
    if file_exists:
        # Remove rows for equations we're about to re-run so we don't duplicate
        import csv as _csv
        rows = []
        with open(out_path, newline='') as rfh:
            reader = _csv.DictReader(rfh)
            rows = [r for r in reader if r['equation'] not in eq_names_to_run]
        with open(out_path, 'w', newline='') as wfh:
            writer_tmp = _csv.DictWriter(wfh, fieldnames=fieldnames)
            writer_tmp.writeheader()
            writer_tmp.writerows(rows)
        print(f'Kept {len(rows)} existing rows; will append {[e["name"] for e in equations_to_run]}.')
        mode = 'a'

    with open(out_path, mode, newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        if mode == 'w':
            writer.writeheader()

        for eq in equations_to_run:
            name = eq['name']
            L    = eq['sim_kwargs'].get('L', 2*np.pi)
            print(f'\n{"="*60}')
            print(f'  Equation: {name}')
            print(f'{"="*60}')

            # ── simulate ────────────────────────────────────────────────
            t0 = time.time()
            print(f'  Simulating...', flush=True)
            snapshots, times = eq['simulate_fn'](**eq['sim_kwargs'])
            print(f'  done in {time.time()-t0:.1f}s  '
                  f'({len(snapshots)} snapshots)', flush=True)

            # ── build library matrices (once per equation) ───────────────
            t0 = time.time()
            print(f'  Building library matrices...', flush=True)
            mats, Ft, desc_eq, desc_std = build_matrices(
                snapshots, times, eq['n_per_snap'], L, eq['skip_fraction'],
                mask_fn=eq.get('mask_fn'),
            )
            print(f'  done in {time.time()-t0:.1f}s  '
                  f'({len(Ft)} regression points)', flush=True)

            true_eq  = eq['true_eq_fn'](desc_eq)
            true_std = np.zeros(len(desc_std))   # std library: no clean true_eq mapping

            # condition numbers
            for (backend, lib), Theta in mats.items():
                kappa = np.linalg.cond(Theta)
                print(f'  cond({backend:10s} {lib:12s}): {kappa:.2e}')

            # ── hyperparameter sweep ─────────────────────────────────────
            combos = list(itertools.product(BACKENDS, LAMS, D_TOLS))
            t0 = time.time()

            pbar = tqdm(
                combos,
                desc=f'  {name}',
                unit='combo',
                ncols=90,
                leave=True,
            )
            for backend, lam, d_tol in pbar:
                pbar.set_postfix(backend=backend, lam=f'{lam:.0e}',
                                 d_tol=f'{d_tol:.0e}')
                Theta_eq  = mats[(backend, 'equivariant')]
                Theta_std = mats[(backend, 'standard')]

                # equivariant library
                nz, fp, fn, exact, err = sweep_one(
                    Theta_eq, Ft, true_eq, lam, d_tol)
                writer.writerow(dict(
                    equation=name, backend=backend, library='equivariant',
                    lam=lam, d_tol=d_tol, n_active=nz, fp=fp, fn=fn,
                    exact_recovery=exact, max_coeff_err=err,
                ))

                # standard library — can only count active terms
                nz_s, fp_s, fn_s, exact_s, err_s = sweep_one(
                    Theta_std, Ft, true_std, lam, d_tol)
                writer.writerow(dict(
                    equation=name, backend=backend, library='standard',
                    lam=lam, d_tol=d_tol, n_active=nz_s, fp=float('nan'),
                    fn=float('nan'), exact_recovery=float('nan'),
                    max_coeff_err=float('nan'),
                ))

            fh.flush()
            print(f'  Sweep done in {time.time()-t0:.1f}s', flush=True)

    print(f'\nResults saved to {out_path}')


if __name__ == '__main__':
    main()
