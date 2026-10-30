#!/usr/bin/env python
"""Pointwise-library cost scaling across dimensions 2–5 for Allen-Cahn, clean data.

Measures build time, STRidge fit time, and peak RAM for all five libraries
(multigraph, mod_multigraph, pdefind, wsindy, polyjet) at dimensions 2 through 5
on a fixed N=24 grid.  Generates fresh clean data for each dimension on first run.

Output
------
CSV with columns: ndim, library, n_cols, n_grid_points, build_time_s, fit_time_s,
peak_rss_gb.  One row per (ndim, library).

Usage
-----
    OMP_NUM_THREADS=1 python scripts/cost_scaling.py
    OMP_NUM_THREADS=1 python scripts/cost_scaling.py --name run1
    OMP_NUM_THREADS=1 python scripts/cost_scaling.py --ndim 2 3 4 5 --libraries all
"""

import argparse
import csv
import datetime
import gc
import os
import resource
import sys
import time

import numpy as np

# ── Enforce single-threaded BLAS *before* numpy re-imports anything ─────
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives
from equivariant_pde_find.library import (
    build_mod_multigraph_library,
    build_multigraph_library,
    build_pdefind_library,
    build_polyjet_library,
    build_wsindy_library,
)
from equivariant_pde_find.regression import STRidge
from equivariant_pde_find.simulations import make_regression_data
from equivariant_pde_find.simulations.allen_cahn import allen_cahn_nd

DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── Fixed parameters ───────────────────────────────────────────────────
N_GRID = 24                     # grid points per direction (user-specified)
L_DOMAIN = 2 * np.pi            # domain side length
EPS = 0.15                      # Allen-Cahn interface width
DT = 0.1
T_END = 20.0
SAVE_EVERY = 5                  # → 41 snapshots
SEED = 42
N_ROWS_REGRESSION = 10_000      # rows for STRidge timing
STRIDGE_LAM = 1e-4
STRIDGE_TOL = 1.0
STRIDGE_MAXIT = 25
STRIDGE_INNER = 10

# Per-dimension snapshot cap to keep build times manageable.
# At ndim=5, even 1 snapshot is 7.96M grid points — enough for a fair
# per-snapshot timing comparison.  The full Θ build time scales linearly
# with snapshots, so per-snapshot times are the meaningful metric.
MAX_SNAPS_BY_NDIM = {2: 39, 3: 39, 4: 20, 5: 5}

# ── Library builders (pointwise, all accept Lx, Ly, and Lz for 3-D fields) ─

def _lz_for_field(f):
    """Return Lz=L_DOMAIN for 3-D fields, None otherwise.

    ndim=3 requires explicit Lz; ndim>3 reuses Lx (isotropic); ndim<=2 optional.
    """
    return L_DOMAIN if f.ndim == 3 else None


LIBRARY_BUILDERS = {
    "multigraph":     lambda f: build_multigraph_library(
        f, L_DOMAIN, L_DOMAIN, max_p=2, max_n=3,
        formulation="pointwise", Lz=_lz_for_field(f)),
    "mod_multigraph": lambda f: build_mod_multigraph_library(
        f, L_DOMAIN, L_DOMAIN, max_p=2, max_n=3,
        formulation="pointwise", Lz=_lz_for_field(f)),
    "pdefind":        lambda f: build_pdefind_library(
        f, L_DOMAIN, L_DOMAIN, max_order=4, max_poly=3,
        formulation="pointwise", Lz=_lz_for_field(f)),
    "wsindy":         lambda f: build_wsindy_library(
        f, L_DOMAIN, L_DOMAIN, max_order=4, max_poly=3,
        formulation="pointwise", Lz=_lz_for_field(f)),
    "polyjet":        lambda f: build_polyjet_library(
        f, L_DOMAIN, L_DOMAIN, max_order=4, max_degree=3,
        formulation="pointwise", Lz=_lz_for_field(f)),
}

# ── Helpers ────────────────────────────────────────────────────────────


def _read_status_field(field):
    """Read a memory field from /proc/self/status (Linux only)."""
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith(field + ":"):
                return float(line.split()[1])  # value in kB
    return 0.0


def current_rss_gb():
    """Current resident set size (VmRSS) in GiB."""
    return _read_status_field("VmRSS") / (1024.0 * 1024.0)


def peak_rss_gb():
    """Peak resident set size (VmHWM, the high-water mark) in GiB.

    This is the maximum RSS the process has used since birth.  We record
    the delta between successive calls to isolate the peak attributable to
    a single library build.
    """
    return _read_status_field("VmHWM") / (1024.0 * 1024.0)


def dataset_path(ndim):
    return os.path.join(DATA_DIR, f"ac_nd{ndim}_coststudy_clean.npz")


def generate_data(ndim):
    """Simulate clean Allen-Cahn at dimension *ndim*, N=24, and save as .npz."""
    path = dataset_path(ndim)
    if os.path.exists(path):
        print(f"  Data already exists: {os.path.basename(path)} — skipping")
        return

    print(f"  Generating Allen-Cahn ndim={ndim}, N={N_GRID} ... ", end="", flush=True)
    t0 = time.perf_counter()
    snaps, times = allen_cahn_nd(
        ndim=ndim, N=N_GRID, L=L_DOMAIN, eps=EPS,
        dt=DT, T_end=T_END, save_every=SAVE_EVERY, seed=SEED,
    )
    elapsed = time.perf_counter() - t0
    np.savez_compressed(
        path,
        snapshots=snaps,
        times=times,
        Lx=L_DOMAIN,
        Ly=L_DOMAIN,
        Lz=L_DOMAIN,
        ndim=ndim,
        N=N_GRID,
        eps=EPS,
    )
    print(f"done ({elapsed:.1f}s, shape={snaps.shape})")


# ── Memory guard ─────────────────────────────────────────────────────
MAX_RAM_GB = 96               # hard ceiling (system limit)
SAFE_RAM_GB = 80              # soft ceiling — skip if estimated peak exceeds this
BYTES_PER_FLOAT64 = 8


def _estimate_peak_gb(n_snaps, n_grid, n_cols):
    """Estimate peak RAM for building Theta on *n_snaps* snapshots.

    Accounts for: Theta array, derivative cache (~ n_derivs × n_grid × 8),
    FFT workspace (~ 2 × n_grid × 8), and a 50 % safety margin.
    """
    theta_gb = n_snaps * n_grid * n_cols * BYTES_PER_FLOAT64 / (1024**3)
    # derivative cache: roughly n_derivs copies of f (one per multi-index)
    # plus 2× for FFT workspace
    overhead_gb = (n_cols // 2) * n_grid * BYTES_PER_FLOAT64 / (1024**3)
    return 1.5 * (theta_gb + overhead_gb)


# ── Main measurement logic ─────────────────────────────────────────────


def measure_one(ndim, lib_key, snaps, times):
    """Build pointwise Theta, time it, run STRidge, return a result dict."""
    result = {
        "ndim": ndim,
        "library": lib_key,
        "n_cols": "",
        "n_grid_points": "",
        "build_time_s": "",
        "fit_time_s": "",
        "peak_rss_gb": "",
        "status": "ok",
    }

    n_grid = int(np.prod(snaps.shape[1:]))  # N^ndim
    result["n_grid_points"] = n_grid

    # ── Build regression data ────────────────────────────────────────
    f_snaps, ft_snaps, _idx, _dt = make_regression_data(
        snaps, times, n_per_snap=1, seed=SEED,
    )

    # ── Memory guard: probe first snapshot to learn n_cols ───────────
    gc.collect()
    try:
        T_probe, desc = LIBRARY_BUILDERS[lib_key](f_snaps[0])
        del T_probe
    except Exception as exc:
        result["status"] = f"error (probe): {exc}"
        print(f"    {lib_key:>16s}: PROBE ERROR — {exc}")
        return result

    n_cols = len(desc)
    result["n_cols"] = n_cols
    full_snaps = len(f_snaps)
    max_snaps = MAX_SNAPS_BY_NDIM.get(ndim, full_snaps)
    full_snaps = min(full_snaps, max_snaps)
    f_snaps = f_snaps[:full_snaps]
    ft_snaps = ft_snaps[:full_snaps]

    # Determine how many snapshots we can afford
    est_full = _estimate_peak_gb(full_snaps, n_grid, n_cols)
    n_use = full_snaps
    while n_use > 1 and _estimate_peak_gb(n_use, n_grid, n_cols) > SAFE_RAM_GB:
        n_use = max(1, n_use // 2)

    if _estimate_peak_gb(n_use, n_grid, n_cols) > SAFE_RAM_GB:
        result["status"] = "skipped (est > 80 GB)"
        print(f"    {lib_key:>16s}: {n_cols:>4d} cols  SKIP — "
              f"est peak {est_full:.1f} GB > {SAFE_RAM_GB} GB guard")
        return result

    n_rows_total = n_use * n_grid
    desc_str = f"{n_cols} cols × {n_use}/{full_snaps} snaps = {n_rows_total:,} rows"

    # ── Build Theta (pre-allocated, no list accumulation) ─────────────
    try:
        gc.collect()
        hwm_before = peak_rss_gb()
        current_rss = current_rss_gb()
        if current_rss > MAX_RAM_GB - 2:
            result["status"] = "skipped (RSS > 94 GB before build)"
            result["peak_rss_gb"] = f"{current_rss:.2f}"
            print(f"    {lib_key:>16s}: SKIP — RSS already {current_rss:.1f} GB")
            return result

        t0 = time.perf_counter()

        # Pre-allocate and fill row-by-row to avoid 2× memory spike
        Theta = np.empty((n_rows_total, n_cols), dtype=np.float64)
        Ft = np.empty(n_rows_total, dtype=np.float64)
        row = 0
        for s in range(n_use):
            T_s, _ = LIBRARY_BUILDERS[lib_key](f_snaps[s])
            n_rows_s = T_s.shape[0]
            Theta[row:row + n_rows_s] = T_s
            Ft[row:row + n_rows_s] = ft_snaps[s].ravel()
            row += n_rows_s
            # Monitor RSS during build; abort if we exceed hard ceiling
            if s % 5 == 0 and current_rss_gb() > MAX_RAM_GB - 1:
                del Theta, Ft, T_s
                gc.collect()
                result["status"] = "oom_guard"
                result["peak_rss_gb"] = f"{peak_rss_gb():.2f}"
                print(f"    {lib_key:>16s}: ABORT — RSS exceeded {MAX_RAM_GB - 1:.0f} GB "
                      f"at snapshot {s}/{n_use}")
                return result

        build_time = time.perf_counter() - t0
        hwm_after = peak_rss_gb()
        build_peak = max(hwm_after - hwm_before, 0.001)
        result["build_time_s"] = f"{build_time:.3f}"
        result["peak_rss_gb"] = f"{hwm_after:.3f}"

        print(f"    {lib_key:>16s}: {desc_str:>40s}  "
              f"build={build_time:.1f}s  HWM={hwm_after:.1f} GiB  "
              f"Δ={build_peak:.1f} GiB", end="", flush=True)

        # ── Subsample for regression ─────────────────────────────────
        rng = np.random.default_rng(SEED)
        if n_rows_total > N_ROWS_REGRESSION:
            idx_sub = rng.choice(n_rows_total, size=N_ROWS_REGRESSION, replace=False)
            Theta_sub = Theta[idx_sub]
            Ft_sub = Ft[idx_sub]
        else:
            Theta_sub = Theta
            Ft_sub = Ft

        del Theta, Ft
        gc.collect()

        # ── STRidge timing ───────────────────────────────────────────
        t0 = time.perf_counter()
        w = STRidge(
            Theta_sub, Ft_sub,
            lam=STRIDGE_LAM, maxit=STRIDGE_MAXIT, tol=STRIDGE_TOL,
            normalize=2, early_stop=True, unbias=True,
        )
        fit_time = time.perf_counter() - t0

        result["fit_time_s"] = f"{fit_time:.4f}"
        n_active = int(np.sum(np.abs(w) > 1e-10))
        print(f"  fit={fit_time:.3f}s  active={n_active}")

    except MemoryError:
        result["status"] = "oom"
        result["peak_rss_gb"] = f"{peak_rss_gb():.2f}"
        print(f"    {lib_key:>16s}: OOM  HWM={peak_rss_gb():.2f} GiB")
    except Exception as exc:
        result["status"] = f"error: {exc}"
        print(f"    {lib_key:>16s}: ERROR — {exc}")

    gc.collect()
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ndim", nargs="+", type=int, default=[2, 3, 4, 5],
                    help="Dimensions to test (default: 2 3 4 5).")
    ap.add_argument("--libraries", nargs="+",
                    default=["multigraph", "mod_multigraph", "pdefind", "wsindy", "polyjet"],
                    help="Libraries to test (default: all five).")
    ap.add_argument("--name", default="", help="Tag for output filename.")
    ap.add_argument("--skip-generate", action="store_true",
                    help="Skip data generation (assume .npz already exist).")
    ap.add_argument("--no-fit", action="store_true",
                    help="Skip regression fit; only measure build time/RAM.")
    args = ap.parse_args()

    derivatives.set_backend("spectral")

    # ── Phase 1: Generate data ───────────────────────────────────────
    if not args.skip_generate:
        print("=== Phase 1: Data generation ===")
        for ndim in args.ndim:
            generate_data(ndim)
        print()

    # ── Phase 2: Measure build + fit for each (ndim, library) ────────
    print("=== Phase 2: Build & regression timing ===")
    results = []

    for ndim in args.ndim:
        path = dataset_path(ndim)
        if not os.path.exists(path):
            print(f"SKIP ndim={ndim}: {os.path.basename(path)} not found. "
                  f"Run without --skip-generate first.")
            continue

        data = np.load(path, allow_pickle=True)
        snaps = data["snapshots"]
        times = data["times"]
        n_grid = int(np.prod(snaps.shape[1:]))
        print(f"\nndim={ndim}  N={N_GRID}  grid_points={n_grid:,}  "
              f"snapshots={len(snaps)}")

        for lib_key in args.libraries:
            row = measure_one(ndim, lib_key, snaps, times)
            results.append(row)

    # ── Phase 3: Write CSV ───────────────────────────────────────────
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    suffix = f"_{args.name}" if args.name else ""
    out_path = os.path.join(RESULTS_DIR, f"cost_scaling_{stamp}{suffix}.csv")

    fieldnames = ["ndim", "library", "n_cols", "n_grid_points",
                  "build_time_s", "fit_time_s", "peak_rss_gb", "status"]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"\nResults saved to: {out_path}")

    # ── Quick summary ─────────────────────────────────────────────────
    print(f"\n{'ndim':>5s}  {'library':>16s}  {'cols':>5s}  "
          f"{'grid_pts':>12s}  {'build(s)':>9s}  {'fit(s)':>9s}  {'RSS(GiB)':>9s}  {'status'}")
    print("-" * 90)
    for row in results:
        print(f"{row['ndim']:>5d}  {row['library']:>16s}  {str(row['n_cols']):>5s}  "
              f"{str(row['n_grid_points']):>12s}  {str(row['build_time_s']):>9s}  "
              f"{str(row['fit_time_s']):>9s}  {str(row['peak_rss_gb']):>9s}  "
              f"{row['status']}")


if __name__ == "__main__":
    main()
