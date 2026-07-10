"""Quick test: exhaustive best-subset on Cahn-Hilliard, equivariant only."""
import numpy as np
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from equivariant_pde_find.library import build_equivariant_library
from equivariant_pde_find.regression import exhaustive_best_subset
from equivariant_pde_find.equation_configs import true_cahn_hilliard
from equivariant_pde_find.simulations import make_regression_data
from equivariant_pde_find.derivatives import set_backend

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')

CONFIGS = [
    ("cahn_hilliard_clean.npz",     "clean", None),
    ("cahn_hilliard_noise2pct.npz", "2%",    None),
]
BACKENDS = ["spectral", "polynomial"]

for fname, noise_label, _ in CONFIGS:
    data = np.load(os.path.join(DATA_DIR, fname))
    snapshots = data["snapshots"]
    times = data["times"]
    Lx, Ly = float(data["Lx"]), float(data["Ly"])
    eps = float(data["eps"]) if "eps" in data else 0.2

    print(f"\n{'='*70}")
    print(f"CH {noise_label}  |  eps={eps}  |  snapshots={snapshots.shape}  |  L={Lx:.1f}")
    print(f"{'='*70}")

    for backend in BACKENDS:
        set_backend(backend)

        # Build regression data with spatial subsampling (matching sweep.py)
        f_snaps, ft_snaps, idx_sets, _ = make_regression_data(
            snapshots, times, n_per_snap=200, skip_fraction=0.0, seed=0,
        )

        # Build pointwise library, subsample spatially
        rows_T, rows_Ft = [], []
        for s in range(len(f_snaps)):
            T, desc = build_equivariant_library(f_snaps[s], Lx, Ly)
            rows_T.append(T[idx_sets[s]])
            rows_Ft.append(ft_snaps[s].ravel()[idx_sets[s]])
        Theta = np.vstack(rows_T).astype(float)
        Ft = np.concatenate(rows_Ft).astype(float)

        # Ground truth
        true_w = true_cahn_hilliard(desc, eps=eps)
        true_idx = list(np.where(np.abs(true_w) > 1e-10)[0])
        true_k = len(true_idx)

        print(f"\n  [{backend}]  Theta={Theta.shape}  cond={np.linalg.cond(Theta):.2e}")
        print(f"  True model: k={true_k}, indices={true_idx}")
        for i in true_idx:
            print(f"    [{i:2d}] {desc[i]:30s} = {true_w[i]:+.6f}")

        # Exhaustive search with L0 penalty
        t0 = time.perf_counter()
        best_by_k, best = exhaustive_best_subset(
            Theta, Ft, max_k=6, l0_penalty=0.0, verbose=False,
        )
        elapsed = time.perf_counter() - t0

        k_opt, idx_opt, xi_opt, rss_opt, score_opt = best
        idx_sorted = sorted(idx_opt)

        # Metrics
        active_set = set(idx_sorted)
        true_set = set(true_idx)
        fp = len(active_set - true_set)
        fn = len(true_set - active_set)
        exact = fp == 0 and fn == 0
        max_err = max(abs(xi_opt[i] - true_w[i]) for i in true_idx) if fn == 0 else float("nan")

        # Also try with l0_penalty tuned to prefer the true model size
        best_tuned = min(best_by_k, key=lambda k: best_by_k[k][3])

        print(f"  Exhaustive (max_k=6): {elapsed:.1f}s")
        print(f"  Best (λ=0): k={k_opt}, fp={fp}, fn={fn}, exact={exact}, "
              f"relRSS={rss_opt / (Ft @ Ft):.6f}, max_err={max_err}")
        for i in idx_sorted:
            marker = " *" if i in true_set else "  "
            print(f"    [{i:2d}] {desc[i]:30s} = {xi_opt[i]:+.6f}  (true={true_w[i]:+.6f}){marker}")

        # Show score by k for a range of λ values
        yty = float(Ft @ Ft)
        print(f"  Pareto (k → relRSS, scores for λ values):")
        for k in sorted(best_by_k):
            rss_k, idx_k, _, _ = best_by_k[k]
            rel = rss_k / yty
            scores = [f"λ={lm:.0e}:{rel + lm*k:.6f}" for lm in [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]]
            idx_str = str(sorted(idx_k))
            true_overlap = len(set(idx_k) & true_set)
            print(f"    k={k:2d}: relRSS={rel:.6f}  overlap={true_overlap}/{true_k}  {', '.join(scores[:4])}")
            if len(scores) > 4:
                print(f"           {' ' * 32}{', '.join(scores[4:])}")

print("\nDone.")
