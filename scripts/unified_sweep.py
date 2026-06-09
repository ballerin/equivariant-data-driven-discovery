"""
Unified ablation sweep: compare pointwise / WSINDy / IBP across equations,
noise levels, libraries, backends, and regression algorithms.

Loads pre-generated .npz data files and runs systematic hyperparameter
sweeps for every method combination.  Results are evaluated with the
universal pointwise-Theta metric (evaluation.py).

Usage
-----
    cd paper-experimental/code

    # Tier 1: core equivariant vs standard comparison
    python scripts/unified_sweep.py --methods pointwise --algorithms strridge

    # Tier 2: noise robustness (WSINDy vs pointwise on SH + KS)
    python scripts/unified_sweep.py \\
        --equations swift_hohenberg ks_chaotic ks_weak_chaotic \\
        --noise 0.02 0.05 \\
        --methods pointwise wsindy_direct wsindy_ibp \\
        --libraries equivariant

    # Tier 3: algorithm sanity check
    python scripts/unified_sweep.py --noise 0.0 --algorithms strridge sr3 \\
        --libraries equivariant --backends spectral

    # Full sweep (everything)
    python scripts/unified_sweep.py
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
from evaluation import build_eval_theta, evaluate_from_theta, evaluate_standard, evaluate_standard_with_truth  # noqa: E402
from library import build_equivariant_library, build_standard_library  # noqa: E402
from regression import TrainSTRidge, TrainSR3  # noqa: E402
from simulations import make_regression_data  # noqa: E402
import wsindy as ws  # noqa: E402

RESULTS_DIR = os.path.join(CODE_DIR, "results")
DATA_DIR = os.path.join(CODE_DIR, "data")
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyperparameter grids ────────────────────────────────────────────────────
LAMS = [1e-5]                       # STRidge lam (historical: not applied)
SR3_KAPPAS = [0.01, 0.1, 1.0, 10.0, 100.0]
D_TOLS = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, 5.0, 15.0, 20.0]
BACKENDS = ["spectral", "polynomial"]
WSINDY_M_VALUES = [4, 8, 12]
WSINDY_P = 4
SAMPLE_STRIDE = 4

# ── equation metadata ───────────────────────────────────────────────────────
# Maps equation names to their skip_fraction, mask_fn, true_eq_fn.
# Loaded from the .npz file's equation name.

def _find_op(desc, name):
    """Return column index of operator `name` in description list."""
    try:
        return desc.index(name)
    except ValueError:
        raise KeyError(f"Operator {name!r} not in library. Available: {desc}")


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


def _true_ks(desc_eq):
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, '‖∇f‖²')] = -0.5
    w[_find_op(desc_eq, 'Δf')] = -1.0
    w[_find_op(desc_eq, 'Δ²f')] = -1.0
    return w


def _true_cahn_hilliard(desc_eq, eps=0.2):
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, 'Δf')] = -1.0
    w[_find_op(desc_eq, 'f²·Δf')] = 3.0
    w[_find_op(desc_eq, 'f·‖∇f‖²')] = 6.0
    w[_find_op(desc_eq, 'Δ²f')] = -(eps ** 2)
    return w


def _true_pme(desc_eq):
    w = np.zeros(len(desc_eq))
    w[_find_op(desc_eq, 'f²·Δf')] = 1.0
    w[_find_op(desc_eq, 'f·‖∇f‖²')] = 2.0
    return w


def _mask_pme(f):
    return f.ravel() > 0.05


# Equation registry: maps name prefix → (true_fn, skip_fraction, mask_fn)
# KS variants (ks_cellular, ks_weak_chaotic, ks_chaotic) all use _true_ks.
EQUATION_META = {
    "allen_cahn":       (_true_allen_cahn,       0.0,  None),
    "swift_hohenberg":  (_true_swift_hohenberg,  0.1,  None),
    "cahn_hilliard":    (_true_cahn_hilliard,     0.0,  None),
    "porous_medium":    (_true_pme,               0.0,  _mask_pme),
    "ks_cellular":      (_true_ks,                0.3,  None),
    "ks_weak_chaotic":  (_true_ks,                0.3,  None),
    "ks_chaotic":       (_true_ks,                0.3,  None),
}


def _get_meta(name):
    """Return (true_fn, skip_fraction, mask_fn) for a dataset name."""
    # Try exact match first, then prefix match
    if name in EQUATION_META:
        return EQUATION_META[name]
    for prefix in sorted(EQUATION_META, key=len, reverse=True):
        if name.startswith(prefix):
            return EQUATION_META[prefix]
    raise KeyError(f"No metadata for equation: {name}")


# ── loading ─────────────────────────────────────────────────────────────────


def load_dataset(npz_path):
    """Load a .npz data file and return (snapshots, times, Lx, Ly, N, dt_snap,
    obs_noise, eq_name)."""
    data = np.load(npz_path)
    snapshots = data["snapshots"]
    times = data["times"]
    Lx = float(data["Lx"])
    Ly = float(data["Ly"])
    N = int(data["N"])
    dt_snap = float(data["dt_snap"])
    obs_noise = float(data["obs_noise"])
    eq_name = os.path.basename(npz_path).rsplit("_", 1)[0]
    # Strip noise label suffix
    for suffix in ["_clean", "_noise2pct", "_noise5pct"]:
        if eq_name.endswith(suffix):
            eq_name = eq_name[: -len(suffix)]
            break
    return snapshots, times, Lx, Ly, N, dt_snap, obs_noise, eq_name


# ── method builders ─────────────────────────────────────────────────────────


def build_pointwise(f_snaps, ft_snaps, idx_sets, Lx, Ly, library, backend):
    """Build pointwise (Theta, Ft) for a given library and backend."""
    with derivatives.use_backend(backend):
        rows_theta, rows_ft = [], []
        for s in range(len(f_snaps)):
            if library == "equivariant":
                T, desc = build_equivariant_library(f_snaps[s], Lx, Ly)
            else:
                T, desc = build_standard_library(f_snaps[s], Lx, Ly,
                                                  max_order=4, max_poly=3)
            rows_theta.append(T[idx_sets[s]])
            rows_ft.append(ft_snaps[s].ravel()[idx_sets[s]])
    Theta = np.vstack(rows_theta).astype(float)
    Ft = np.concatenate(rows_ft).astype(float)
    return Theta, Ft, desc


def build_wsindy(f_snaps, ft_snaps, idx_sets, Lx, Ly, library, backend,
                 m, use_ibp):
    """Build weak-form (G, b) for a given library and backend."""
    with derivatives.use_backend(backend):
        rows_G, rows_b = [], []
        for s in range(len(f_snaps)):
            f = f_snaps[s]
            ft = ft_snaps[s]
            if library == "equivariant":
                G, desc = ws.build_weak_library_2d(
                    f, Lx, Ly, max_p=2, max_n=3, m=m, p=WSINDY_P,
                    use_ibp=use_ibp,
                )
            else:
                G, desc = ws.build_weak_standard_library_2d(
                    f, Lx, Ly, max_order=4, max_poly=3, m=m, p=WSINDY_P,
                    use_ibp=use_ibp,
                )
            b_w = ws.build_weak_rhs_2d(ft, m=m, p=WSINDY_P)

            Nx, Ny = f.shape
            Nv_x, Nv_y = Nx - 2 * m, Ny - 2 * m

            # Strided subsample
            ix = np.arange(0, Nv_x, SAMPLE_STRIDE)
            iy = np.arange(0, Nv_y, SAMPLE_STRIDE)
            IX, IY = np.meshgrid(ix, iy, indexing="ij")
            flat_idx = (IX * Nv_y + IY).ravel()

            # Take n_per_snap random from strided points
            rng = np.random.default_rng(s)
            n_per = len(idx_sets[s])  # same number as pointwise
            chosen = rng.choice(len(flat_idx), size=min(n_per, len(flat_idx)),
                                replace=False)
            q_idx = flat_idx[chosen]

            rows_G.append(G[q_idx])
            rows_b.append(b_w[q_idx])

    G_all = np.vstack(rows_G).astype(float)
    b_all = np.concatenate(rows_b).astype(float)
    return G_all, b_all, desc


# ── sweep one ───────────────────────────────────────────────────────────────


def sweep_one(Theta, Ft, true_eq, lam, d_tol, algorithm):
    """Run regression and return (n_active, w)."""
    if algorithm == "sr3":
        w = TrainSR3(Theta, Ft, lam, d_tol)
    else:
        w = TrainSTRidge(Theta, Ft, lam, d_tol)
    n_active = int(np.count_nonzero(np.abs(w) > 1e-10))
    return n_active, w


# ── compute condition number ────────────────────────────────────────────────


def cond_number(Theta):
    """Condition number (2-norm) of Theta."""
    return float(np.linalg.cond(Theta))


# ── main ────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Unified ablation sweep for equivariant PDE discovery."
    )
    parser.add_argument(
        "--data-dir", default=DATA_DIR,
        help="Directory containing .npz data files.",
    )
    parser.add_argument(
        "--equations", nargs="+", default=None,
        help="Which equation datasets to sweep (default: all .npz in data dir).",
    )
    parser.add_argument(
        "--noise", nargs="+", type=float, default=None,
        help="Noise levels to sweep (default: all found).",
    )
    parser.add_argument(
        "--methods", nargs="+",
        default=["pointwise", "wsindy_direct", "wsindy_ibp"],
        choices=["pointwise", "wsindy_direct", "wsindy_ibp"],
    )
    parser.add_argument(
        "--algorithms", nargs="+",
        default=["strridge", "sr3"],
        choices=["strridge", "sr3"],
    )
    parser.add_argument(
        "--libraries", nargs="+",
        default=["equivariant", "standard"],
        choices=["equivariant", "standard"],
    )
    parser.add_argument(
        "--backends", nargs="+",
        default=BACKENDS,
        choices=BACKENDS,
    )
    parser.add_argument(
        "--wsindy-m", nargs="+", type=int,
        default=WSINDY_M_VALUES,
        help="Test-function half-support values for WSINDy methods.",
    )
    parser.add_argument(
        "--name", default="",
        help="Optional tag appended to output filename.",
    )
    args = parser.parse_args()

    # ── discover datasets ────────────────────────────────────────────────
    data_dir = args.data_dir
    all_npz = sorted([
        f for f in os.listdir(data_dir)
        if f.endswith(".npz") and not f.startswith(".")
    ])

    if args.equations:
        eq_set = set(args.equations)
        all_npz = [f for f in all_npz
                   if any(f.startswith(eq + "_") for eq in eq_set)]
    if args.noise is not None:
        noise_labels = {0.0: "clean", 0.02: "noise2pct", 0.05: "noise5pct"}
        allowed = {noise_labels[n] for n in args.noise}
        all_npz = [f for f in all_npz
                   if any(f.endswith("_" + nl + ".npz") for nl in allowed)]

    if not all_npz:
        print(f"No .npz files found in {data_dir} matching filters.")
        return

    methods = args.methods
    algorithms = args.algorithms
    libraries = args.libraries
    backends = args.backends
    wsindy_m_vals = args.wsindy_m

    # ── output ───────────────────────────────────────────────────────────
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"_{args.name}" if args.name else ""
    out_path = os.path.join(RESULTS_DIR, f"unified_sweep_{timestamp}{tag}.csv")
    print(f"Output → {out_path}", flush=True)

    fieldnames = [
        "equation", "noise", "method", "algorithm", "library", "backend",
        "wsindy_m", "lam", "d_tol",
        "n_active", "fp", "fn", "exact_recovery", "E", "max_coeff_err",
        "cond_theta", "n_samples",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        f.flush()

        for npz_file in all_npz:
            npz_path = os.path.join(data_dir, npz_file)
            print(f"\n{'='*70}")
            print(f"Dataset: {npz_file}")
            print(f"{'='*70}")

            # Load
            snapshots, times, Lx, Ly, N, dt_snap, obs_noise, eq_name = \
                load_dataset(npz_path)
            true_fn, skip_frac, mask_fn = _get_meta(eq_name)

            print(f"  eq={eq_name}  noise={obs_noise}  "
                  f"snaps={snapshots.shape[0]}  N={N}  L={Lx:.1f}",
                  flush=True)

            # Prepare regression data (shared for all pointwise methods)
            f_snaps_raw, ft_snaps_raw, idx_sets_raw, _ = make_regression_data(
                snapshots, times, n_per_snap=200, skip_fraction=skip_frac,
                seed=0,
            )

            # Apply mask if needed
            if mask_fn is not None:
                rng = np.random.default_rng(1)
                new_idx = []
                for s in range(len(f_snaps_raw)):
                    valid = np.where(mask_fn(f_snaps_raw[s]))[0]
                    if len(valid) <= 200:
                        new_idx.append(valid)
                    else:
                        new_idx.append(np.sort(
                            rng.choice(valid, size=200, replace=False)))
                f_snaps = f_snaps_raw
                ft_snaps = ft_snaps_raw
                idx_sets = new_idx
            else:
                f_snaps = f_snaps_raw
                ft_snaps = ft_snaps_raw
                idx_sets = idx_sets_raw

            n_snaps = len(f_snaps)
            print(f"  interior snaps: {n_snaps}  "
                  f"points: {n_snaps * 200}", flush=True)

            # ── Pre-build evaluation Theta + ground truth (once per dataset)
            t_eval = time.time()
            Theta_eval, Ft_eval = build_eval_theta(
                snapshots, times, Lx, Ly,
                skip_fraction=skip_frac, seed=0,
            )
            # Ground truth is always in the equivariant basis.
            _, desc_eq_ref = build_equivariant_library(f_snaps[0], Lx, Ly)
            true_eq = true_fn(desc_eq_ref)
            print(f"  eval Theta: {Theta_eval.shape}  "
                  f"active_true={int(np.sum(np.abs(true_eq)>1e-10))}  "
                  f"({time.time() - t_eval:.1f}s)", flush=True)

            # ── Pre-build method matrices ────────────────────────────────
            # We build once per (method, library, backend, wsindy_m) and
            # sweep over (algorithm, lam, d_tol).
            # Store as: (Theta, Ft, desc, cond_val)

            prebuilt = {}  # key → (Theta, Ft, desc, cond_val)

            total_combos = 0
            for method, library, backend in itertools.product(
                methods, libraries, backends
            ):
                if method == "pointwise":
                    m_vals = [None]
                else:
                    m_vals = wsindy_m_vals

                for m_val in m_vals:
                    key = (method, library, backend, m_val, False)
                    if method == "wsindy_ibp":
                        key = (method, library, backend, m_val, True)

                    # Also key without use_ibp for wsindy_direct
                    if method == "wsindy_direct":
                        use_ibp_flag = False
                    elif method == "wsindy_ibp":
                        use_ibp_flag = True
                    else:
                        use_ibp_flag = False

                    actual_key = (method, library, backend,
                                  m_val if m_val is not None else 0)

                    if actual_key in prebuilt:
                        continue

                    t0 = time.time()
                    print(f"  Building: method={method} lib={library} "
                          f"backend={backend}"
                          + (f" m={m_val}" if m_val is not None else ""),
                          flush=True)

                    if method == "pointwise":
                        Theta, Ft, desc = build_pointwise(
                            f_snaps, ft_snaps, idx_sets, Lx, Ly,
                            library, backend,
                        )
                    else:
                        Theta, Ft, desc = build_wsindy(
                            f_snaps, ft_snaps, idx_sets, Lx, Ly,
                            library, backend, m_val, use_ibp_flag,
                        )

                    cond_val = cond_number(Theta)
                    prebuilt[actual_key] = (Theta, Ft, desc, cond_val)
                    elapsed = time.time() - t0

                    print(f"    shape={Theta.shape}  cond={cond_val:.2e}  "
                          f"({elapsed:.1f}s)", flush=True)

            # ── Sweep hyperparams ──────────────────────────────────────
            for method, library, backend in itertools.product(
                methods, libraries, backends
            ):
                if method == "pointwise":
                    m_vals = [0]  # placeholder
                    use_ibp_flag = False
                elif method == "wsindy_direct":
                    m_vals = wsindy_m_vals
                    use_ibp_flag = False
                else:  # wsindy_ibp
                    m_vals = wsindy_m_vals
                    use_ibp_flag = True

                for m_val in m_vals:
                    actual_key = (method, library, backend, m_val)
                    if actual_key not in prebuilt:
                        continue
                    Theta, Ft, desc, cond_val = prebuilt[actual_key]

                    for algorithm in algorithms:
                        if algorithm == "sr3":
                            lam_grid = SR3_KAPPAS
                        else:
                            lam_grid = LAMS

                        for lam in lam_grid:
                            for d_tol in D_TOLS:
                                n_active, w = sweep_one(
                                    Theta, Ft, true_eq, lam, d_tol, algorithm,
                                )

                                # Evaluate: for standard library, decompose
                                # equivariant truth into standard columns and
                                # compute proper fp/fn where representable.
                                if library == "standard":
                                    result = evaluate_standard_with_truth(
                                        w, true_eq, desc_eq_ref,
                                        Theta, desc, Ft,
                                    )
                                    eval_result = dict(
                                        fp=result["fp"],
                                        fn=result["fn"],
                                        exact_recovery=result["exact_recovery"],
                                        n_active=result["n_active"],
                                        E=result["E"],
                                        max_coeff_err=result["max_coeff_err"],
                                        n_samples=result["n_samples"],
                                    )
                                else:
                                    eval_result = evaluate_from_theta(
                                        w, true_eq, Theta_eval, Ft_eval)

                                row = dict(
                                    equation=eq_name,
                                    noise=obs_noise,
                                    method=method,
                                    algorithm=algorithm,
                                    library=library,
                                    backend=backend,
                                    wsindy_m=(m_val if m_val != 0 else ""),
                                    lam=lam,
                                    d_tol=d_tol,
                                    n_active=n_active,
                                    fp=eval_result["fp"],
                                    fn=eval_result["fn"],
                                    exact_recovery=eval_result["exact_recovery"],
                                    E=eval_result["E"],
                                    max_coeff_err=eval_result["max_coeff_err"],
                                    cond_theta=cond_val,
                                    n_samples=eval_result["n_samples"],
                                )
                                writer.writerow(row)
                                total_combos += 1
                                if total_combos % 50 == 0:
                                    f.flush()
                                    print(f"  [{eq_name} noise={obs_noise}] "
                                          f"{total_combos} combos, last: "
                                          f"{method} {library} {backend} "
                                          f"d_tol={d_tol} n_active={n_active}",
                                          flush=True)

            f.flush()
            print(f"  [{eq_name} noise={obs_noise}] done: {total_combos} rows",
                  flush=True)

    print(f"\nDone. {total_combos} rows → {out_path}")


if __name__ == "__main__":
    main()
