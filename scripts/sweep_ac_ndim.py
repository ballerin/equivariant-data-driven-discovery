#!/usr/bin/env python
"""Allen-Cahn dimension-scaling study: how recovery and cost change with ndim.

Fixed configuration (everything except the dimension is held constant, so
ndim is the only variable):

    equation    Allen-Cahn,  N=32, eps=0.3, L=2*pi at every ndim
    formulation weak form with integration by parts, m=12
    backend     spectral only
    algorithm   STRidge, one (lam, d_tol) pair calibrated at ndim=2
    libraries   multigraph / mod_multigraph / pdefind

The libraries scale very differently with dimension, which is the point of
the study: the equivariant bases stay at 20 columns in every dimension while
PDE-FIND's coordinate library grows 60 -> 140 -> 280 -> 504.

Usage
-----
    # 1. calibrate (lam, d_tol) on the 2-D data
    python scripts/sweep_ac_ndim.py --calibrate

    # 2. run the study with the winning pair
    python scripts/sweep_ac_ndim.py --lam 1e-4 --d-tol 1.0 --name study

    # single cell, e.g. to probe the cost of the top dimension first
    python scripts/sweep_ac_ndim.py --ndim 5 --noise 0.0 \\
        --libraries multigraph --lam 1e-4 --d-tol 1.0 --name probe
"""

import argparse
import csv
import datetime
import os
import resource
import sys
import time

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src")))
sys.path.insert(0, SCRIPT_DIR)

from equivariant_pde_find import derivatives  # noqa: E402
from equivariant_pde_find.equations import (  # noqa: E402
    true_allen_cahn,
    true_allen_cahn_mod_multigraph,
    true_allen_cahn_pdefind,
)
from equivariant_pde_find.library import (  # noqa: E402
    build_mod_multigraph_library,
    build_multigraph_library,
    build_pdefind_library,
)
from equivariant_pde_find.library import weak_form as wf  # noqa: E402
from equivariant_pde_find.regression import TrainSTRidge  # noqa: E402
from equivariant_pde_find.simulations import make_regression_data  # noqa: E402

# evaluate_native defines fp/fn/exact/E_coeff for a native-basis recovery;
# imported rather than copied so this study and sweep_nd.py cannot drift.
from sweep_nd import evaluate_native  # noqa: E402

DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
os.makedirs(RESULTS_DIR, exist_ok=True)

TAG = "dimstudy"
WEAK_M = 12
WEAK_MAX_DERIV = 4          # highest derivative order moved onto the test function
MAX_ROWS = 2000             # same row budget at every ndim, so cost differences
                            # come from n_cols and N^ndim, not from sample size
NOISE_LABELS = {0.0: "clean", 0.02: "noise2pct", 0.05: "noise5pct"}

# Calibration grids, matching scripts/sweep.py
CAL_LAMS = [1e-6, 1e-4, 1e-2, 1.0]
CAL_D_TOLS = [3e-3, 1e-2, 3e-2, 0.1, 0.2, 0.3, 0.5, 1.0]

LIBRARIES = {
    "multigraph": (build_multigraph_library, true_allen_cahn),
    "mod_multigraph": (build_mod_multigraph_library, true_allen_cahn_mod_multigraph),
    "pdefind": (build_pdefind_library, None),  # truth is ndim-parameterised
}


def peak_rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def dataset_path(ndim, noise):
    return os.path.join(DATA_DIR, f"allen_cahn_nd{ndim}_{TAG}_{NOISE_LABELS[noise]}.npz")


def load(ndim, noise):
    d = np.load(dataset_path(ndim, noise))
    eps = float(d["eps"]) if "eps" in d else 0.15
    return d["snapshots"], d["times"], float(d["L"]), eps


def truth_vector(library, desc, ndim, eps):
    """Ground-truth coefficients in the library's own basis.

    Allen-Cahn is f_t = eps^2 Delta f + f - f^3, so the support is
    {f, f^3, Delta f} for the equivariant bases and {f, f^3} plus one pure
    second derivative per axis for PDE-FIND.
    """
    if library == "pdefind":
        truth = true_allen_cahn_pdefind(ndim, eps=eps)
        w = np.zeros(len(desc))
        for name, coeff in truth.items():
            w[desc.index(name)] = coeff
        expected_nnz = 2 + ndim
    else:
        w = LIBRARIES[library][1](desc, eps=eps)
        expected_nnz = 3

    nnz = int(np.count_nonzero(w))
    if nnz != expected_nnz:
        raise RuntimeError(
            f"{library} at ndim={ndim}: truth vector has {nnz} active terms, "
            f"expected {expected_nnz}. A silently-empty or mis-wired truth "
            f"vector would make fp/fn/E_coeff meaningless."
        )
    return w


def build_weak(library, snaps, times, L, ndim, max_rows, seed=0):
    """Weak-form (G, b) with IBP at m=WEAK_M, spectral backend.

    Uses only as many interior snapshots as the row budget needs, so the
    build cost is not inflated by data the regression will not use -- and
    takes them from the **end** of the trajectory.

    That last part is not cosmetic.  Allen-Cahn starts from a small random
    field, so early snapshots sit in the linear-growth regime where f^3 is
    invisible and the library is near-degenerate.  Taking the first
    snapshots instead of the last gave cond(G) = 3.2e7 at ndim=3 (versus
    3.0e3 at ndim=2) purely because the higher dimension needed fewer
    snapshots and therefore stopped earlier in time -- an artefact that
    looks exactly like a dimension effect.  Sampling the tail makes every
    dimension see the same developed state.
    """
    p = wf.select_test_function_degree(WEAK_M, WEAK_MAX_DERIV)
    if not wf.test_function_is_resolved(WEAK_M, p, WEAK_MAX_DERIV):
        raise RuntimeError(f"m={WEAK_M} cannot resolve a degree-{p} test function")

    f_snaps, ft_snaps, _idx, _dt = make_regression_data(
        snaps, times, n_per_snap=1, seed=seed
    )
    interior_times = times[1:len(f_snaps) + 1]

    n_valid = (snaps.shape[1] - 2 * WEAK_M) ** ndim
    n_needed = max(1, int(np.ceil(max_rows / max(n_valid, 1))))
    n_use = min(len(f_snaps), n_needed)

    builder = LIBRARIES[library][0]
    rng = np.random.default_rng(seed)
    rows_G, rows_b, desc = [], [], None
    per_snap = int(np.ceil(max_rows / n_use))

    for s in range(len(f_snaps) - n_use, len(f_snaps)):
        kw = dict(formulation="weak", m=WEAK_M, p=p, use_ibp=True, Lz=L)
        if library == "pdefind":
            G, desc = builder(f_snaps[s], L, L, max_order=4, max_poly=3, **kw)
        else:
            G, desc = builder(f_snaps[s], L, L, max_p=2, max_n=3, **kw)
        b = wf.build_weak_rhs(ft_snaps[s], m=WEAK_M, p=p)

        take = min(per_snap, G.shape[0])
        sel = rng.choice(G.shape[0], size=take, replace=False)
        rows_G.append(G[sel])
        rows_b.append(b[sel])
        del G

    G_all = np.vstack(rows_G)[:max_rows]
    b_all = np.concatenate(rows_b)[:max_rows]
    used = interior_times[len(f_snaps) - n_use:]
    window = (float(used[0]), float(used[-1])) if len(used) else (float("nan"),) * 2
    return G_all, b_all, desc, n_use, window


def normalize_system(G, b):
    """Scale (G, b) by 1/||b||.

    The recovered coefficients are unchanged -- both sides are divided by the
    same scalar -- but the quantity STRidge thresholds, |w_j| * ||G_j||,
    becomes the dimensionless ratio |w_j| * ||G_j|| / ||b||.

    This is what makes a single d_tol transferable across dimensions.  The
    weak form integrates against a test function supported on (2m+1)^ndim
    grid points, so ||G_j|| grows by more than an order of magnitude per
    dimension: measured on the equivariant library, |w_j| * ||G_j|| for the
    *true* terms runs 7.7e2 (ndim=2) -> 6.9e3 (ndim=3) -> 2.2e4 (ndim=4),
    while the spurious terms sit at 5.9 / 8.7 / 24.  The separation is
    excellent everywhere, but a fixed threshold calibrated at ndim=2 lands
    above the true terms at ndim=2 and far below them at ndim=4, so it keeps
    every column.  After this scaling the same threshold works at all three.
    """
    scale = np.linalg.norm(b)
    if scale < 1e-30:
        return G, b
    return G / scale, b / scale


def run_cell(library, ndim, noise, lam, d_tol, max_rows, timeout_s,
             normalize=True):
    """One (library, ndim, noise) cell.  Returns a result row."""
    row = dict(ndim=ndim, noise=noise, library=library, lam=lam, d_tol=d_tol,
               status="ok", n_cols="", n_rows="", n_snapshots="", n_active="",
               fp="", fn="", exact_recovery="", E_coeff="", E_coeff_ls="",
               build_time_s="", fit_time_s="", total_time_s="",
               peak_rss_mb="", t_first="", t_last="")
    t_start = time.time()
    try:
        snaps, times, L, eps = load(ndim, noise)

        t0 = time.time()
        G, b, desc, n_use, window = build_weak(
            library, snaps, times, L, ndim, max_rows)
        build_time = time.time() - t0

        if timeout_s and build_time > timeout_s:
            row.update(status="timeout", build_time_s=f"{build_time:.2f}",
                       total_time_s=f"{time.time() - t_start:.2f}")
            return row

        w_true = truth_vector(library, desc, ndim, eps)

        if normalize:
            G, b = normalize_system(G, b)

        t0 = time.time()
        w = TrainSTRidge(G, b, lam=lam, d_tol=d_tol)
        fit_time = time.time() - t0

        true_dict = {d: w_true[i] for i, d in enumerate(desc) if w_true[i] != 0.0}
        ev = evaluate_native(w, desc, true_dict)

        # Unthresholded least squares on the *true* support: if this is
        # accurate then the library and the weak form are fine at this
        # dimension, and any failure above is the sparsity threshold.
        sup = np.flatnonzero(w_true)
        w_ls_sup = np.linalg.lstsq(G[:, sup], b, rcond=None)[0]
        E_ls = float(np.linalg.norm(w_ls_sup - w_true[sup])
                     / max(np.linalg.norm(w_true[sup]), 1e-30))

        row.update(
            n_cols=len(desc), n_rows=G.shape[0], n_snapshots=n_use,
            t_first=f"{window[0]:.2f}", t_last=f"{window[1]:.2f}",
            n_active=ev["n_active"], fp=ev["fp"], fn=ev["fn"],
            exact_recovery=ev["exact_recovery"], E_coeff=f"{ev['E_coeff']:.6f}",
            E_coeff_ls=f"{E_ls:.6f}",
            build_time_s=f"{build_time:.2f}", fit_time_s=f"{fit_time:.3f}",
            total_time_s=f"{time.time() - t_start:.2f}",
            peak_rss_mb=f"{peak_rss_mb():.0f}",
        )
    except MemoryError:
        row.update(status="oom", total_time_s=f"{time.time() - t_start:.2f}",
                   peak_rss_mb=f"{peak_rss_mb():.0f}")
    except FileNotFoundError:
        row.update(status="missing_data")
    return row


def calibrate(ndim, noises, libraries, max_rows, normalize=True):
    """Rank (lam, d_tol) on the study's exact configuration at one dimension."""
    print(f"Calibrating (lam, d_tol) at ndim={ndim} on {len(libraries)} libraries, "
          f"{len(noises)} noise levels\n")
    # Build each (library, noise) once; the grid only changes the regression.
    built = {}
    for library in libraries:
        for noise in noises:
            snaps, times, L, eps = load(ndim, noise)
            G, b, desc, _, _ = build_weak(library, snaps, times, L, ndim, max_rows)
            if normalize:
                G, b = normalize_system(G, b)
            built[(library, noise)] = (G, b, desc,
                                       truth_vector(library, desc, ndim, eps))

    scores = []
    for lam in CAL_LAMS:
        for d_tol in CAL_D_TOLS:
            n_exact, E_sum, n = 0, 0.0, 0
            for (library, noise), (G, b, desc, w_true) in built.items():
                w = TrainSTRidge(G, b, lam=lam, d_tol=d_tol)
                true_dict = {d: w_true[i] for i, d in enumerate(desc)
                             if w_true[i] != 0.0}
                ev = evaluate_native(w, desc, true_dict)
                n_exact += ev["exact_recovery"]
                E_sum += ev["E_coeff"]
                n += 1
            scores.append((n_exact, E_sum / n, lam, d_tol, n))

    scores.sort(key=lambda s: (-s[0], s[1]))
    print(f"{'lam':>10} {'d_tol':>8} {'exact':>10} {'mean E_coeff':>13}")
    for n_exact, E_mean, lam, d_tol, n in scores:
        print(f"{lam:>10.0e} {d_tol:>8g} {n_exact:>4}/{n:<5} {E_mean:>13.4f}")
    best = scores[0]
    print(f"\nBest: lam={best[2]:g}  d_tol={best[3]:g}  "
          f"({best[0]}/{best[4]} exact, mean E_coeff {best[1]:.4f})")
    return best[2], best[3]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ndim", nargs="+", type=int, default=[2, 3, 4, 5])
    ap.add_argument("--noise", nargs="+", type=float, default=[0.0, 0.02])
    ap.add_argument("--libraries", nargs="+", default=list(LIBRARIES),
                    choices=list(LIBRARIES))
    ap.add_argument("--lam", type=float, default=None,
                    help="STRidge ridge penalty (required unless --calibrate).")
    ap.add_argument("--d-tol", type=float, default=None, dest="d_tol",
                    help="STRidge initial tolerance step.")
    ap.add_argument("--max-rows", type=int, default=MAX_ROWS)
    ap.add_argument("--timeout", type=float, default=0.0,
                    help="Per-cell build budget in seconds (0 = none).")
    ap.add_argument("--calibrate", action="store_true",
                    help="Rank (lam, d_tol) at the lowest requested ndim and exit.")
    ap.add_argument("--no-normalize", action="store_true",
                    help="Skip the 1/||b|| scaling (see normalize_system).")
    ap.add_argument("--name", default="")
    args = ap.parse_args()

    derivatives.set_backend("spectral")

    if args.calibrate:
        calibrate(min(args.ndim), args.noise, args.libraries, args.max_rows,
                  normalize=not args.no_normalize)
        return

    if args.lam is None or args.d_tol is None:
        ap.error("--lam and --d-tol are required (or use --calibrate first)")

    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    suffix = f"_{args.name}" if args.name else ""
    out_path = os.path.join(RESULTS_DIR, f"ac_ndim_{stamp}{suffix}.csv")

    fields = ["ndim", "noise", "library", "lam", "d_tol", "status", "n_cols",
              "n_rows", "n_snapshots", "n_active", "fp", "fn",
              "exact_recovery", "E_coeff", "E_coeff_ls", "build_time_s",
              "fit_time_s",
              "total_time_s", "peak_rss_mb", "t_first", "t_last"]

    print(f"Allen-Cahn dimension study  ->  {out_path}")
    print(f"  ndim={args.ndim}  noise={args.noise}  libraries={args.libraries}")
    print(f"  lam={args.lam:g}  d_tol={args.d_tol:g}  max_rows={args.max_rows}")
    print(f"  weak IBP, m={WEAK_M}, p={wf.select_test_function_degree(WEAK_M, WEAK_MAX_DERIV)},"
          f" spectral backend\n")

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        # ascending ndim, so partial results are useful if the top is too slow
        for ndim in sorted(args.ndim):
            for noise in args.noise:
                for library in args.libraries:
                    print(f"  ndim={ndim} noise={noise:g} {library:15s} ...",
                          end="", flush=True)
                    row = run_cell(library, ndim, noise, args.lam, args.d_tol,
                                   args.max_rows, args.timeout,
                                   normalize=not args.no_normalize)
                    writer.writerow(row)
                    fh.flush()
                    if row["status"] == "ok":
                        print(f" cols={row['n_cols']:>4} fp={row['fp']} fn={row['fn']}"
                              f" exact={row['exact_recovery']}"
                              f" E={row['E_coeff']}"
                              f" build={row['build_time_s']}s"
                              f" fit={row['fit_time_s']}s")
                    else:
                        print(f" {row['status'].upper()}")

    print(f"\nDone -> {out_path}")


if __name__ == "__main__":
    main()
