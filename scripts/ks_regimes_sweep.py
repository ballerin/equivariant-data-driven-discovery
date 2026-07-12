"""
KS domain-size regime sweep: run the hyperparameter sweep on four KS configurations
spanning ordered-cellular → weakly-chaotic → chaotic → strongly-chaotic regimes.

The Kuramoto–Sivashinsky linear growth rate is σ(K) = K² − K⁴, with peak at
K* = 1/√2 and unstable band 0 < K < 1.  The number of unstable modes scales
as N_unstable ≈ L²/(4π):

    L =  4π → N_unstable ≈   8  (ordered cellular, quasi-periodic)
    L =  8π → N_unstable ≈  44  (weakly chaotic)
    L = 16π → N_unstable ≈ 192  (chaotic — confirmed exact recovery)
    L = 24π → N_unstable ≈ 436  (strongly chaotic)

NOTE: L=4π gives a quasi-periodic field (autocorr sign-change at Δt≈10) that
collapses the design matrix (cond~200k vs ~155 for L=8π).  Regression is
expected to fail at L=4π regardless of library choice — this serves as a
controlled negative example showing that spatiotemporal diversity is required.

Usage
-----
    cd paper-experimental/code
    python scripts/ks_regimes_sweep.py [--algorithm strridge|sr3] [--name TAG]
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


# ── path setup ──────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives  # noqa: E402
from equivariant_pde_find.library import build_multigraph_library, build_pdefind_library  # noqa: E402
from equivariant_pde_find.regression import TrainSTRidge, TrainSR3  # noqa: E402
from equivariant_pde_find.simulations import kuramoto_sivashinsky, make_regression_data  # noqa: E402
from equivariant_pde_find.equation_configs import true_ks  # noqa: E402

RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyperparameter grid ──────────────────────────────────────────────────────
# STRidge thresholds (normalised coefficients).  The KS false-positive '1'
# column has |w_norm| ≈ sqrt(n_train) × |coeff_raw| ≈ 14–16; D_TOLS must
# include values that let STRidge exceed that within maxit=25 steps.
LAMS = [1e-5]
D_TOLS = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, 5.0, 15.0, 20.0]
BACKENDS = ["spectral", "polynomial"]
SR3_KAPPAS = [0.01, 0.1, 1.0, 10.0, 100.0]

# ── KS regime configurations ─────────────────────────────────────────────────
# All regimes use the same time-discretisation (dt=0.1, save_every=5) and
# spatial resolution (N=64) so that the regression infrastructure is identical.
# Only L varies.  skip_fraction=0.3 discards the linear-growth transient.
#
# Diagnostic statistics (T_end=200, skip 30%, seed=42):
#   L=4π:  std≈3.04, drift≈-1.57/t, autocorr(Δt=10)≈-0.57, cond(samp)≈200k
#   L=8π:  std≈3.47, drift≈-2.60/t, autocorr(Δt=10)≈+0.05, cond(samp)≈155
#   L=16π: std≈4.00, drift≈-2.59/t, autocorr(Δt=10)≈+0.14, cond(samp)≈76
#   L=24π: std≈4.49, drift≈-2.88/t, autocorr(Δt=10)≈+0.00, cond(samp)≈522
KS_REGIMES = [
    dict(
        # Ordered cellular: only 8 unstable modes, quasi-periodic field.
        # Design matrix condition number ~200k; regression expected to fail.
        name="ks_4pi",
        label="L=4π  (cellular, N_u≈8)",
        L=4 * np.pi,
        N_unstable=8,
        note="quasi-periodic / degenerate for regression",
        sim_kwargs=dict(N=64, L=4 * np.pi, T_end=200.0, dt=0.1, save_every=5, seed=42),
        n_per_snap=200,
        skip_fraction=0.3,
    ),
    dict(
        # Weakly chaotic: 44 unstable modes, decaying autocorrelation.
        # Well-conditioned design matrix (~155); exact recovery expected.
        name="ks_8pi",
        label="L=8π  (weakly chaotic, N_u≈44)",
        L=8 * np.pi,
        N_unstable=44,
        note="weakly chaotic",
        sim_kwargs=dict(N=64, L=8 * np.pi, T_end=200.0, dt=0.1, save_every=5, seed=42),
        n_per_snap=200,
        skip_fraction=0.3,
    ),
    dict(
        # Chaotic: 192 unstable modes.  Already confirmed: exact recovery at
        # d_tol=5.0 and 15.0 with spectral equivariant library.
        name="ks_16pi",
        label="L=16π (chaotic, N_u≈192)",
        L=16 * np.pi,
        N_unstable=192,
        note="chaotic — confirmed exact recovery",
        sim_kwargs=dict(N=64, L=16 * np.pi, T_end=200.0, dt=0.1, save_every=5, seed=42),
        n_per_snap=200,
        skip_fraction=0.3,
    ),
    dict(
        # Strongly chaotic: 436 unstable modes, fastest temporal decorrelation.
        name="ks_24pi",
        label="L=24π (strongly chaotic, N_u≈436)",
        L=24 * np.pi,
        N_unstable=436,
        note="strongly chaotic",
        sim_kwargs=dict(N=64, L=24 * np.pi, T_end=200.0, dt=0.1, save_every=5, seed=42),
        n_per_snap=200,
        skip_fraction=0.3,
    ),
]

# Ground truth: KS = -Δf - Δ²f - ½‖∇f‖² (imported from equation_configs)


# ── helpers ──────────────────────────────────────────────────────────────────

def build_matrices(snapshots, times, n_per_snap, L, skip_fraction):
    """Build equivariant and standard library matrices for both backends."""
    f_snaps, ft_snaps, idx_sets, _ = make_regression_data(
        snapshots, times, n_per_snap=n_per_snap, seed=0,
    )
    n = len(f_snaps)
    start = int(np.floor(skip_fraction * n))
    f_snaps  = f_snaps[start:]
    ft_snaps = ft_snaps[start:]
    idx_sets = idx_sets[start:]

    Ft = np.concatenate([ft_snaps[s].ravel()[idx_sets[s]] for s in range(len(f_snaps))])

    eq_spec_rows, std_spec_rows = [], []
    eq_poly_rows, std_poly_rows = [], []
    t_spec_total = t_poly_total = 0.0
    t_start = time.time()
    n_snaps = len(f_snaps)
    print(f"  Building libraries for {n_snaps} snapshots...", flush=True)

    for s in tqdm(range(n_snaps), desc="  libraries", unit="snap", ncols=90, leave=True):
        fs  = f_snaps[s]
        idx = idx_sets[s]

        t0 = time.time()
        with derivatives.use_backend("spectral"):
            Teq,  desc_eq  = build_multigraph_library(fs, L, L)
            Tstd, desc_std = build_pdefind_library(fs, L, L, max_order=4, max_poly=3)
        eq_spec_rows.append(Teq[idx])
        std_spec_rows.append(Tstd[idx])
        t_spec_total += time.time() - t0

        t0 = time.time()
        with derivatives.use_backend("polynomial"):
            Teq,  _ = build_multigraph_library(fs, L, L)
            Tstd, _ = build_pdefind_library(fs, L, L, max_order=4, max_poly=3)
        eq_poly_rows.append(Teq[idx])
        std_poly_rows.append(Tstd[idx])
        t_poly_total += time.time() - t0

        if (s + 1) % 10 == 0 or s == n_snaps - 1:
            elapsed = time.time() - t_start
            rate = (s + 1) / elapsed
            eta  = (n_snaps - s - 1) / rate if rate > 0 else float("inf")
            print(
                f"    snap {s+1}/{n_snaps}  "
                f"spectral {t_spec_total/(s+1)*1e3:.0f}ms/snap  "
                f"poly {t_poly_total/(s+1)*1e3:.0f}ms/snap  "
                f"ETA {eta:.0f}s",
                flush=True,
            )

    print(
        f"  Library build done in {time.time()-t_start:.1f}s  "
        f"(spectral {t_spec_total:.1f}s  polynomial {t_poly_total:.1f}s)",
        flush=True,
    )
    mats = {
        ("spectral",    "multigraph"): np.vstack(eq_spec_rows).astype(float),
        ("spectral",    "pdefind"):    np.vstack(std_spec_rows).astype(float),
        ("polynomial",  "multigraph"): np.vstack(eq_poly_rows).astype(float),
        ("polynomial",  "pdefind"):    np.vstack(std_poly_rows).astype(float),
    }
    return mats, Ft, desc_eq, desc_std


def sweep_one(Theta, Ft, true_eq, lam, d_tol, algorithm="strridge"):
    """Run regression and return (n_active, fp, fn, exact_recovery, max_coeff_err)."""
    act = np.abs(true_eq) > 1e-10
    w   = TrainSR3(Theta, Ft, lam, d_tol) if algorithm == "sr3" else TrainSTRidge(Theta, Ft, lam, d_tol)
    nz  = int(np.count_nonzero(np.abs(w) > 1e-10))
    fp  = int(np.sum((np.abs(w) > 1e-10) & ~act))
    fn  = int(np.sum((np.abs(w) < 1e-10) &  act))
    exact = int(fp == 0 and fn == 0)
    found = act & (np.abs(w) > 1e-10)
    err = float(np.max(np.abs(w[act] - true_eq[act]))) if found.any() else float("nan")
    return nz, fp, fn, exact, err


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="KS domain-size regime sweep (four L values)."
    )
    parser.add_argument(
        "--algorithm",
        choices=["strridge", "sr3"],
        default="strridge",
        help="Regression algorithm (default: strridge).",
    )
    parser.add_argument(
        "--name",
        default="",
        help="Optional tag appended to the output filename.",
    )
    args = parser.parse_args()

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"_{args.name}" if args.name else ""
    out_path = os.path.join(RESULTS_DIR, f"ks_regimes_{timestamp}{tag}.csv")
    print(f"Output → {out_path}")
    print(f"Algorithm: {args.algorithm}")
    print()

    fieldnames = [
        "regime",          # ks_4pi / ks_8pi / ks_16pi / ks_24pi
        "L_over_pi",       # 4 / 8 / 16 / 24  (integer, for easy filtering)
        "N_unstable",      # approximate number of unstable modes
        "backend",
        "library",
        "algorithm",
        "lam",
        "d_tol",
        "n_active",
        "fp",
        "fn",
        "exact_recovery",
        "max_coeff_err",
    ]

    lam_grid   = SR3_KAPPAS if args.algorithm == "sr3" else LAMS
    lam_label  = "kappa"    if args.algorithm == "sr3" else "lam"

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for regime in KS_REGIMES:
            name        = regime["name"]
            label       = regime["label"]
            L           = regime["L"]
            L_over_pi   = int(round(L / np.pi))
            N_u         = regime["N_unstable"]
            note        = regime["note"]
            sim_kwargs  = regime["sim_kwargs"]
            n_per_snap  = regime["n_per_snap"]
            skip_frac   = regime["skip_fraction"]

            print(f'\n{"="*60}')
            print(f"  Regime: {label}")
            print(f"  Note: {note}")
            print(f'{"="*60}')

            # ── simulate ──────────────────────────────────────────────────
            t0 = time.time()
            print("  Simulating...", flush=True)
            snapshots, times = kuramoto_sivashinsky(**sim_kwargs)
            print(f"  done in {time.time()-t0:.1f}s  ({len(snapshots)} snapshots)", flush=True)

            # ── brief temporal diagnostics ────────────────────────────────
            skip = int(0.3 * len(snapshots))
            late = snapshots[skip:]
            f_ms = late - late.mean(axis=(1,2), keepdims=True)
            spatial_stds = np.std(f_ms, axis=(1,2))
            mean_vals    = late.mean(axis=(1,2))
            drift_rate   = np.polyfit(times[skip:], mean_vals, 1)[0]
            f0_flat = f_ms[0].ravel()
            f0_norm = np.dot(f0_flat, f0_flat)
            lag10_idx = min(20, len(f_ms) - 1)  # lag=20 snapshots → Δt=10
            autocorr10 = np.dot(f_ms[lag10_idx].ravel(), f0_flat) / f0_norm if f0_norm > 0 else float("nan")
            print(f"  Diagnostics (late {len(late)} snaps): "
                  f"std={np.mean(spatial_stds):.2f}, "
                  f"drift={drift_rate:.2f}/t, "
                  f"autocorr(Δt=10)={autocorr10:+.3f}",
                  flush=True)

            # ── build library matrices ────────────────────────────────────
            t0 = time.time()
            print(f"  Building library matrices (n_per_snap={n_per_snap})...", flush=True)
            mats, Ft, desc_eq, desc_std = build_matrices(
                snapshots, times, n_per_snap, L, skip_frac
            )
            print(f"  done in {time.time()-t0:.1f}s  ({len(Ft)} regression points)", flush=True)

            true_eq  = true_ks(desc_eq)
            true_std = np.zeros(len(desc_std))

            # condition numbers
            for (backend, lib), Theta in mats.items():
                kappa = np.linalg.cond(Theta)
                print(f"  cond({backend:11s} {lib:12s}): {kappa:.2e}")

            # ── hyperparameter sweep ──────────────────────────────────────
            combos = list(itertools.product(BACKENDS, lam_grid, D_TOLS))
            t0     = time.time()
            n_combos = len(combos)
            print(f"  Running {n_combos} (backend, {lam_label}, d_tol) combos "
                  f"[algorithm={args.algorithm}]...", flush=True)

            pbar = tqdm(combos, desc=f"  {name}", unit="combo", ncols=90, leave=True)
            t_combo_start = time.time()
            for i, (backend, lam, d_tol) in enumerate(pbar):
                pbar.set_postfix(backend=backend[:4], lam=f"{lam:.0e}", d_tol=f"{d_tol:.0e}")
                if i > 0 and i % 20 == 0:
                    avg = (time.time() - t_combo_start) / i
                    eta = avg * (n_combos - i)
                    print(f"    {i}/{n_combos}  avg {avg*1e3:.0f}ms/combo  ETA {eta:.0f}s", flush=True)

                Theta_eq  = mats[(backend, "multigraph")]
                Theta_std = mats[(backend, "pdefind")]

                nz, fp, fn, exact, err = sweep_one(
                    Theta_eq, Ft, true_eq, lam, d_tol, algorithm=args.algorithm
                )
                writer.writerow(dict(
                    regime=name, L_over_pi=L_over_pi, N_unstable=N_u,
                    backend=backend, library="multigraph", algorithm=args.algorithm,
                    lam=lam, d_tol=d_tol,
                    n_active=nz, fp=fp, fn=fn, exact_recovery=exact, max_coeff_err=err,
                ))

                # standard library — no ground truth mapping; fp/fn are nan
                nz_s, _, _, _, _ = sweep_one(
                    Theta_std, Ft, true_std, lam, d_tol, algorithm=args.algorithm
                )
                writer.writerow(dict(
                    regime=name, L_over_pi=L_over_pi, N_unstable=N_u,
                    backend=backend, library="pdefind", algorithm=args.algorithm,
                    lam=lam, d_tol=d_tol,
                    n_active=nz_s, fp=float("nan"), fn=float("nan"),
                    exact_recovery=float("nan"), max_coeff_err=float("nan"),
                ))

            fh.flush()
            print(f"  Sweep done in {time.time()-t0:.1f}s", flush=True)

            # ── best-result summary ───────────────────────────────────────
            import csv as _csv
            with open(out_path, newline="") as rfh:
                all_rows = [r for r in _csv.DictReader(rfh) if r["regime"] == name]
            print(f"  Best equivariant results ({label}):")
            for backend in BACKENDS:
                eq_rows = [
                    r for r in all_rows
                    if r["backend"] == backend
                    and r["library"] == "multigraph"
                    and r["exact_recovery"] not in ("", "nan")
                ]
                if not eq_rows:
                    print(f"    {backend:11s}: no valid rows")
                    continue
                exact_rows = [r for r in eq_rows if r["exact_recovery"] == "1"]
                best = min(
                    eq_rows,
                    key=lambda r: (int(r["fn"]), int(r["fp"]),
                                   float(r["max_coeff_err"]) if r["max_coeff_err"] not in ("", "nan") else 1e9),
                )
                mark = " ✓ EXACT" if exact_rows else ""
                print(
                    f"    {backend:11s}: best fn={best['fn']}, fp={best['fp']}, "
                    f"err={float(best['max_coeff_err']):.3e}"
                    f" at d_tol={best['d_tol']}{mark}"
                )

    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
