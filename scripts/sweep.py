"""
Canonical ablation sweep: compare pointwise / weak (w/wo IBP) across equations,
noise levels, libraries, backends, and regression algorithms.

Loads pre-generated .npz data files (produced by generate_data.py) and runs
systematic hyperparameter sweeps for every method combination.  Results are
evaluated with the universal pointwise-Theta metric (evaluation.py).

Supports optional post-hoc noise application and spectral low-pass
preprocessing for additional ablation dimensions.

Usage
-----
    # Core equivariant vs standard comparison
    python scripts/sweep.py --methods pointwise --algorithms strridge

    # Noise robustness (weak vs pointwise)
    python scripts/sweep.py \\
        --equations swift_hohenberg ks_chaotic \\
        --noise 0.02 0.05 \\
        --methods pointwise weak_direct weak_ibp \\
        --libraries equivariant

    # With post-hoc noise and spectral low-pass
    python scripts/sweep.py --post-hoc-noise 0.02 0.05 --spectral-cutoff 3.0

    # Full sweep
    python scripts/sweep.py
"""

import argparse
import csv
import datetime
import itertools
import os
import platform
import subprocess
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
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives  # noqa: E402
from equivariant_pde_find.evaluation import (  # noqa: E402
    build_eval_theta_jet, evaluate_from_theta, evaluate_unified,
    evaluate_mwsindy,
    build_wsindy_to_jet_map, _EQUIV_TO_JET_DECOMPOSITION,
)
from equivariant_pde_find.equation_configs import (  # noqa: E402
    find_op, true_allen_cahn, true_swift_hohenberg, true_ks,
    true_cahn_hilliard, true_pme, mask_pme, EQUATION_META,
)
from equivariant_pde_find.library import (  # noqa: E402
    build_equivariant_library, build_standard_library,
    build_polynomial_jet_library, build_pointwise_wsindy_library,
    build_mod_multigraph_library,
    spectral_lowpass,
)
from equivariant_pde_find.regression import TrainSTRidge, TrainSR3, mstls_sweep  # noqa: E402
from equivariant_pde_find.simulations import make_regression_data  # noqa: E402
from equivariant_pde_find import weak_form as wf  # noqa: E402

RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyperparameter grids ────────────────────────────────────────────────────
LAMS = [1e-6, 1e-4, 1e-2, 1.0]
SR3_KAPPAS = [0.01, 0.1, 1.0, 10.0, 100.0]
MSTLS_N_LAMBDAS = 50
MSTLS_ALPHAS = [0.0, 0.25, 0.5, 0.75, 1.0]
MSTLS_MAXITS = 50
D_TOLS = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0]
BACKENDS = ["spectral", "polynomial"]
WEAK_M = [4, 8, 12]
WEAK_P = 6
SAMPLE_STRIDE = 4

# CH epsilon values for Cahn-Hilliard ablation
CH_EPSILONS = [0.1, 0.2, 0.4]


def _get_meta(name):
    """Return (true_fn, skip_fraction, mask_fn) for a dataset name."""
    if name in EQUATION_META:
        return EQUATION_META[name]
    for prefix in sorted(EQUATION_META, key=len, reverse=True):
        if name.startswith(prefix):
            return EQUATION_META[prefix]
    raise KeyError(f"No metadata for equation: {name}")


# ── loading ─────────────────────────────────────────────────────────────────


def load_dataset(npz_path):
    """Load a .npz data file."""
    data = np.load(npz_path)
    snapshots = data["snapshots"]
    times = data["times"]
    Lx = float(data["Lx"])
    Ly = float(data["Ly"])
    N = int(data["N"])
    dt_snap = float(data["dt_snap"])
    obs_noise = float(data["obs_noise"])
    eq_name = os.path.basename(npz_path).rsplit("_", 1)[0]
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
            elif library == "mod_multigraph":
                T, desc = build_mod_multigraph_library(f_snaps[s], Lx, Ly)
            elif library == "polynomial_jet":
                T, desc = build_polynomial_jet_library(
                    f_snaps[s], Lx, Ly, max_order=4, max_degree=3,
                )
            elif library in ("wsindy", "wsindy_proj", "mwsindy"):
                T, desc = build_pointwise_wsindy_library(
                    f_snaps[s], Lx, Ly, max_order=4, max_poly=3,
                )
            else:
                T, desc = build_standard_library(
                    f_snaps[s], Lx, Ly, max_order=4, max_poly=3,
                )
            rows_theta.append(T[idx_sets[s]])
            rows_ft.append(ft_snaps[s].ravel()[idx_sets[s]])
    Theta = np.vstack(rows_theta).astype(float)
    Ft = np.concatenate(rows_ft).astype(float)
    return Theta, Ft, desc


def build_weak_form(f_snaps, ft_snaps, idx_sets, Lx, Ly, library, backend,
                 m, use_ibp):
    """Build weak-form (G, b) for a given library and backend."""
    with derivatives.use_backend(backend):
        rows_G, rows_b = [], []
        for s in range(len(f_snaps)):
            f = f_snaps[s]
            ft = ft_snaps[s]
            if library == "equivariant":
                G, desc = wf.build_weak_library(
                    f, Lx, Ly, max_p=2, max_n=3, m=m, p=WEAK_P,
                    use_ibp=use_ibp,
                )
            elif library == "mod_multigraph":
                G, desc = wf.build_weak_mod_multigraph_library(
                    f, Lx, Ly, max_p=2, max_n=3, m=m, p=WEAK_P,
                    use_ibp=use_ibp,
                )
            elif library == "polynomial_jet":
                G, desc = wf.build_weak_polynomial_jet_library(
                    f, Lx, Ly, max_order=4, max_degree=3,
                    m=m, p=WEAK_P, use_ibp=use_ibp,
                )
            elif library in ("wsindy", "wsindy_proj", "mwsindy"):
                G, desc = wf.build_weak_wsindy_library(
                    f, Lx, Ly, max_order=4, max_poly=3,
                    m=m, p=WEAK_P, use_ibp=use_ibp,
                )
            else:
                G, desc = wf.build_weak_standard_library(
                    f, Lx, Ly, max_order=4, max_poly=3, m=m, p=WEAK_P,
                    use_ibp=use_ibp,
                )
            b_w = wf.build_weak_rhs(ft, m=m, p=WEAK_P)

            Nx, Ny = f.shape
            Nv_x, Nv_y = Nx - 2 * m, Ny - 2 * m
            ix = np.arange(0, Nv_x, SAMPLE_STRIDE)
            iy = np.arange(0, Nv_y, SAMPLE_STRIDE)
            IX, IY = np.meshgrid(ix, iy, indexing="ij")
            flat_idx = (IX * Nv_y + IY).ravel()

            rng = np.random.default_rng(s)
            n_per = len(idx_sets[s])
            chosen = rng.choice(
                len(flat_idx), size=min(n_per, len(flat_idx)), replace=False,
            )
            q_idx = flat_idx[chosen]

            rows_G.append(G[q_idx])
            rows_b.append(b_w[q_idx])

    G_all = np.vstack(rows_G).astype(float)
    b_all = np.concatenate(rows_b).astype(float)
    return G_all, b_all, desc


# ── sweep one ───────────────────────────────────────────────────────────────


def sweep_one(Theta, Ft, true_eq, lam, d_tol, algorithm, mstls_alpha=0.5):
    """Run regression and return (n_active, w, lam_out)."""
    if algorithm == "sr3":
        w = TrainSR3(Theta, Ft, lam, d_tol)
        lam_out = lam
    elif algorithm == "mstls":
        # MSTLS does its own internal λ sweep; lam/d_tol are ignored
        w, lam_out, _ = mstls_sweep(
            Theta, Ft, lambdas=None, n_lambdas=MSTLS_N_LAMBDAS,
            maxits=MSTLS_MAXITS, alpha=mstls_alpha,
        )
    else:
        w = TrainSTRidge(Theta, Ft, lam, d_tol)
        lam_out = lam
    n_active = int(np.count_nonzero(np.abs(w) > 1e-10))
    return n_active, w, lam_out


def cond_number(Theta):
    """Condition number (2-norm) of Theta."""
    return float(np.linalg.cond(Theta))


# ── main ────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Canonical ablation sweep for equivariant PDE discovery."
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
        default=["pointwise", "weak_direct", "weak_ibp"],
        choices=["pointwise", "weak_direct", "weak_ibp"],
    )
    parser.add_argument(
        "--algorithms", nargs="+",
        default=["strridge", "sr3", "mstls"],
        choices=["strridge", "sr3", "mstls"],
    )
    parser.add_argument(
        "--libraries", nargs="+",
        default=["equivariant", "standard", "polynomial_jet", "wsindy", "wsindy_proj", "mwsindy"],
        choices=["equivariant", "mod_multigraph", "standard", "polynomial_jet", "wsindy", "wsindy_proj", "mwsindy"],
    )
    parser.add_argument(
        "--backends", nargs="+",
        default=BACKENDS,
        choices=BACKENDS,
    )
    parser.add_argument(
        "--weak-m", nargs="+", type=int,
        default=WEAK_M,
        help="Test-function half-support values for weak-form methods.",
    )
    parser.add_argument(
        "--post-hoc-noise", nargs="+", type=float, default=None,
        help="Apply additive Gaussian noise (as fraction of std) to loaded "
             "snapshots. Useful for noise-robustness ablation.",
    )
    parser.add_argument(
        "--spectral-cutoff", type=float, default=None,
        help="Apply spectral low-pass filter to snapshots before building "
             "libraries (cutoff wavenumber).",
    )
    parser.add_argument(
        "--ch-epsilons", nargs="+", type=float, default=None,
        help="Cahn-Hilliard epsilon values to sweep (default: 0.1 0.2 0.4). "
             "Only applies when an equation contains 'cahn_hilliard'.",
    )
    parser.add_argument(
        "--mstls-alpha", nargs="+", type=float, default=None,
        help="MSTLS alpha values (projection-vs-sparsity trade-off). "
             "Default: 0.0 0.25 0.5 0.75 1.0.",
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
    weak_m_vals = args.weak_m
    post_hoc_noises = args.post_hoc_noise or [None]
    spectral_cutoffs = [args.spectral_cutoff] if args.spectral_cutoff else [None]
    ch_epsilons = args.ch_epsilons or CH_EPSILONS
    mstls_alphas = args.mstls_alpha or MSTLS_ALPHAS

    # ── output ───────────────────────────────────────────────────────────
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"_{args.name}" if args.name else ""
    out_path = os.path.join(RESULTS_DIR, f"sweep_{timestamp}{tag}.csv")
    print(f"Output → {out_path}", flush=True)

    fieldnames = [
        "equation", "noise", "method", "algorithm", "library", "backend",
        "weak_m", "lam", "d_tol", "mstls_alpha",
        "n_active", "fp", "fn", "exact_recovery", "E", "max_coeff_err",
        "E_fit", "cond_theta", "n_samples",
        "post_hoc_noise", "spectral_cutoff", "ch_eps",
        "time_ms",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        f.flush()

        t_start_total = time.time()
        total_combos = 0
        eq_summaries = []  # per-equation summary rows

        for npz_file in all_npz:
            npz_path = os.path.join(data_dir, npz_file)
            print(f"\n{'='*70}")
            print(f"Dataset: {npz_file}")
            print(f"{'='*70}")

            # Load
            snapshots_orig, times, Lx, Ly, N, dt_snap, obs_noise, eq_name = \
                load_dataset(npz_path)
            true_fn, skip_frac, mask_fn = _get_meta(eq_name)

            # Determine CH epsilon sweep if applicable
            is_ch = "cahn_hilliard" in eq_name
            eps_list = ch_epsilons if is_ch else [None]

            eq_exact_count = 0
            eq_total = 0

            for post_hoc_noise in post_hoc_noises:
                # Apply post-hoc noise
                if post_hoc_noise is not None and post_hoc_noise > 0:
                    rng_noise = np.random.default_rng(42)
                    snapshots = snapshots_orig.copy()
                    for s in range(len(snapshots)):
                        std = np.std(snapshots[s])
                        snapshots[s] += post_hoc_noise * std * \
                            rng_noise.standard_normal(snapshots[s].shape)
                    noise_label = obs_noise + post_hoc_noise
                else:
                    snapshots = snapshots_orig
                    noise_label = obs_noise

                # Apply spectral cutoff
                if (args.spectral_cutoff is not None
                        and args.spectral_cutoff > 0):
                    snapshots_proc = snapshots.copy()
                    for s in range(len(snapshots_proc)):
                        snapshots_proc[s] = spectral_lowpass(
                            snapshots_proc[s], Lx, Ly, args.spectral_cutoff,
                        )
                else:
                    snapshots_proc = snapshots

                for ch_eps in eps_list:
                    print(f"  eq={eq_name}  noise={noise_label:.3f}"
                          + (f"  cutoff={args.spectral_cutoff}"
                             if args.spectral_cutoff else "")
                          + (f"  ch_eps={ch_eps}" if ch_eps else ""),
                          flush=True)

                    # Build ground truth with appropriate parameters
                    _, desc_eq_ref = build_equivariant_library(
                        snapshots_proc[0], Lx, Ly,
                    )
                    if is_ch:
                        true_eq = true_fn(desc_eq_ref, eps=ch_eps)
                    elif eq_name.startswith("allen_cahn"):
                        true_eq = true_fn(desc_eq_ref)
                    elif eq_name.startswith("swift_hohenberg"):
                        true_eq = true_fn(desc_eq_ref)
                    else:
                        true_eq = true_fn(desc_eq_ref)

                    # Prepare regression data
                    f_snaps_raw, ft_snaps_raw, idx_sets_raw, _ = \
                        make_regression_data(
                            snapshots_proc, times, n_per_snap=200,
                            skip_fraction=skip_frac, seed=0,
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
                                    rng.choice(valid, size=200,
                                              replace=False)))
                        f_snaps = f_snaps_raw
                        ft_snaps = ft_snaps_raw
                        idx_sets = new_idx
                    else:
                        f_snaps = f_snaps_raw
                        ft_snaps = ft_snaps_raw
                        idx_sets = idx_sets_raw

                    n_snaps = len(f_snaps)
                    print(f"    interior snaps: {n_snaps}  "
                          f"points: {n_snaps * 200}", flush=True)

                    # ── Pre-build evaluation Theta (polynomial jet) ──────
                    t_eval = time.time()
                    Theta_eval, Ft_eval = build_eval_theta_jet(
                        snapshots_proc, times, Lx, Ly,
                        skip_fraction=skip_frac, seed=0,
                    )
                    # Build jet descriptions for mapping
                    from equivariant_pde_find.library import \
                        build_polynomial_jet_library as _bpl
                    _, desc_jet = _bpl(
                        snapshots_proc[0], Lx, Ly, max_order=4, max_degree=3,
                    )
                    # Ground truth in jet basis
                    true_jet = np.zeros(len(desc_jet))
                    for j, (name, coeff) in enumerate(zip(desc_eq_ref, true_eq)):
                        if abs(coeff) < 1e-12:
                            continue
                        if name not in _EQUIV_TO_JET_DECOMPOSITION:
                            continue
                        for jet_name, mult in _EQUIV_TO_JET_DECOMPOSITION[name]:
                            true_jet[desc_jet.index(jet_name)] += coeff * mult

                    # Pre-build WSINDy → jet mapping matrix
                    from equivariant_pde_find.library import \
                        build_pointwise_wsindy_library as _bpw
                    _, desc_ws = _bpw(
                        snapshots_proc[0], Lx, Ly, max_order=4, max_poly=3,
                    )
                    M_ws_to_jet = build_wsindy_to_jet_map(desc_ws, desc_jet)

                    print(f"    eval Theta: {Theta_eval.shape}  "
                          f"active_true="
                          f"{int(np.sum(np.abs(true_jet) > 1e-10))}"
                          f"  ({time.time() - t_eval:.1f}s)", flush=True)

                    # ── Pre-build method matrices ────────────────────────
                    prebuilt = {}

                    # MWSINDy needs equivariant G for fine-tuning — ensure
                    # it's built even if not explicitly requested.
                    build_libs = list(libraries)
                    if "mwsindy" in libraries and "equivariant" not in build_libs:
                        build_libs = ["equivariant"] + build_libs

                    for method, library, backend in itertools.product(
                        methods, build_libs, backends
                    ):
                        if method == "pointwise":
                            m_vals = [None]
                        else:
                            m_vals = weak_m_vals

                        for m_val in m_vals:
                            if method == "weak_direct":
                                use_ibp_flag = False
                            elif method == "weak_ibp":
                                use_ibp_flag = True
                            else:
                                use_ibp_flag = False

                            actual_key = (method, library, backend,
                                          m_val if m_val is not None else 0)

                            if actual_key in prebuilt:
                                continue

                            t0 = time.time()
                            print(f"    Building: method={method} "
                                  f"lib={library} backend={backend}"
                                  + (f" m={m_val}" if m_val else ""),
                                  flush=True)

                            if method == "pointwise":
                                Theta, Ft, desc = build_pointwise(
                                    f_snaps, ft_snaps, idx_sets, Lx, Ly,
                                    library, backend,
                                )
                            else:
                                Theta, Ft, desc = build_weak_form(
                                    f_snaps, ft_snaps, idx_sets, Lx, Ly,
                                    library, backend, m_val, use_ibp_flag,
                                )

                            cond_val = cond_number(Theta)
                            prebuilt[actual_key] = (Theta, Ft, desc, cond_val)
                            elapsed = time.time() - t0
                            print(f"      shape={Theta.shape}  "
                                  f"cond={cond_val:.2e}  "
                                  f"({elapsed:.1f}s)", flush=True)

                    # ── Sweep hyperparams ──────────────────────────────
                    for method, library, backend in itertools.product(
                        methods, libraries, backends
                    ):
                        if method == "pointwise":
                            m_vals = [0]
                        elif method == "weak_direct":
                            m_vals = weak_m_vals
                        else:
                            m_vals = weak_m_vals

                        for m_val in m_vals:
                            actual_key = (method, library, backend, m_val)
                            if actual_key not in prebuilt:
                                continue
                            Theta, Ft, desc, cond_val = prebuilt[actual_key]

                            for algorithm in algorithms:
                                if algorithm == "mstls":
                                    # MSTLS does its own internal λ sweep;
                                    # iterate over α values
                                    for mstls_alpha in mstls_alphas:
                                        t_row = time.time()
                                        n_active, w, lam_out = sweep_one(
                                            Theta, Ft, true_eq, 0.0, 0.0,
                                            algorithm,
                                            mstls_alpha=mstls_alpha,
                                        )

                                        eval_result = evaluate_unified(
                                            w, library, desc, desc_jet,
                                            Theta_eval, Ft_eval, true_jet,
                                            M_ws_to_jet=M_ws_to_jet,
                                        )
                                        # MWSINDy: fine-tune on multigraph submatrix
                                        if library == "mwsindy":
                                            mg_key = (method, "equivariant", backend, m_val)
                                            if mg_key in prebuilt:
                                                G_mg, b_mg, desc_mg, _ = prebuilt[mg_key]
                                                eval_result = evaluate_mwsindy(
                                                    w, desc, desc_mg, desc_jet,
                                                    M_ws_to_jet, Theta_eval,
                                                    Ft_eval, true_jet,
                                                    G_mg, b_mg,
                                                )
                                        dt_ms = (time.time() - t_row) * 1000

                                        row = dict(
                                            equation=eq_name,
                                            noise=noise_label,
                                            method=method,
                                            algorithm=algorithm,
                                            library=library,
                                            backend=backend,
                                            weak_m=(m_val if m_val != 0
                                                      else ""),
                                            lam=lam_out,
                                            d_tol="",
                                            mstls_alpha=mstls_alpha,
                                            n_active=n_active,
                                            fp=eval_result["fp"],
                                            fn=eval_result["fn"],
                                            exact_recovery=eval_result[
                                                "exact_recovery"],
                                            E=eval_result["E"],
                                            max_coeff_err=eval_result[
                                                "max_coeff_err"],
                                            E_fit=eval_result.get("E_fit", ""),
                                            cond_theta=cond_val,
                                            n_samples=eval_result["n_samples"],
                                            post_hoc_noise=(
                                                post_hoc_noise
                                                if post_hoc_noise else ""),
                                            spectral_cutoff=(
                                                args.spectral_cutoff
                                                if args.spectral_cutoff
                                                else ""),
                                            ch_eps=(ch_eps if ch_eps else ""),
                                            time_ms=f"{dt_ms:.1f}",
                                        )
                                        writer.writerow(row)
                                        total_combos += 1
                                        if eval_result["exact_recovery"]:
                                            eq_exact_count += 1
                                        eq_total += 1
                                    continue

                                lam_grid = (SR3_KAPPAS if algorithm == "sr3"
                                            else LAMS)

                                for lam in lam_grid:
                                    for d_tol in D_TOLS:
                                        t_row = time.time()
                                        n_active, w, lam_out = sweep_one(
                                            Theta, Ft, true_eq, lam, d_tol,
                                            algorithm,
                                        )

                                        eval_result = evaluate_unified(
                                            w, library, desc, desc_jet,
                                            Theta_eval, Ft_eval, true_jet,
                                            M_ws_to_jet=M_ws_to_jet,
                                        )
                                        # MWSINDy: fine-tune on multigraph submatrix
                                        if library == "mwsindy":
                                            mg_key = (method, "equivariant", backend, m_val)
                                            if mg_key in prebuilt:
                                                G_mg, b_mg, desc_mg, _ = prebuilt[mg_key]
                                                eval_result = evaluate_mwsindy(
                                                    w, desc, desc_mg, desc_jet,
                                                    M_ws_to_jet, Theta_eval,
                                                    Ft_eval, true_jet,
                                                    G_mg, b_mg,
                                                )
                                        dt_ms = (time.time() - t_row) * 1000

                                        row = dict(
                                            equation=eq_name,
                                            noise=noise_label,
                                            method=method,
                                            algorithm=algorithm,
                                            library=library,
                                            backend=backend,
                                            weak_m=(m_val if m_val != 0
                                                      else ""),
                                            lam=lam_out,
                                            d_tol=d_tol,
                                            mstls_alpha="",
                                            n_active=n_active,
                                            fp=eval_result["fp"],
                                            fn=eval_result["fn"],
                                            exact_recovery=eval_result[
                                                "exact_recovery"],
                                            E=eval_result["E"],
                                            max_coeff_err=eval_result[
                                                "max_coeff_err"],
                                            E_fit=eval_result.get("E_fit", ""),
                                            cond_theta=cond_val,
                                            n_samples=eval_result["n_samples"],
                                            post_hoc_noise=(
                                                post_hoc_noise
                                                if post_hoc_noise else ""),
                                            spectral_cutoff=(
                                                args.spectral_cutoff
                                                if args.spectral_cutoff
                                                else ""),
                                            ch_eps=(ch_eps if ch_eps else ""),
                                            time_ms=f"{dt_ms:.1f}",
                                        )
                                        writer.writerow(row)
                                        total_combos += 1
                                        eq_total += 1
                                        if eval_result["exact_recovery"]:
                                            eq_exact_count += 1
                                        if total_combos % 50 == 0:
                                            f.flush()
                                            print(f"    [{eq_name}] "
                                                  f"{total_combos} combos",
                                                  flush=True)
                    f.flush()

            # ── Per-equation summary ─────────────────────────────────────
            eq_summaries.append(
                f"  {eq_name:25s}  exact={eq_exact_count}/{eq_total}  "
                f"({100 * eq_exact_count / max(eq_total, 1):.0f}%)"
            )
            print(f"  [{eq_name}] done: {eq_total} rows, "
                  f"exact={eq_exact_count}", flush=True)

    # ── Final summary ────────────────────────────────────────────────────
    elapsed_total = time.time() - t_start_total
    print(f"\n{'='*70}")
    print("Per-equation exact recovery summary:")
    for s in eq_summaries:
        print(s)
    print(f"\nDone. {total_combos} total rows → {out_path}")
    print(f"Total wall time: {elapsed_total:.1f}s ({elapsed_total/60:.1f} min)")

    # ── Write timing log ────────────────────────────────────────────────
    log_path = out_path.replace(".csv", ".log")
    omp_threads = os.environ.get("OMP_NUM_THREADS", "not set")
    try:
        cpu_info = subprocess.check_output(
            ["lscpu"], stderr=subprocess.DEVNULL, text=True
        )
        cpu_model = [l for l in cpu_info.splitlines() if "Model name" in l]
        cpu_model = cpu_model[0].split(":", 1)[1].strip() if cpu_model else "unknown"
    except Exception:
        cpu_model = "unknown"
    with open(log_path, "w") as lf:
        lf.write(f"# Sweep timing log\n")
        lf.write(f"csv_file: {os.path.basename(out_path)}\n")
        lf.write(f"total_wall_time_s: {elapsed_total:.1f}\n")
        lf.write(f"total_wall_time_min: {elapsed_total/60:.1f}\n")
        lf.write(f"total_rows: {total_combos}\n")
        lf.write(f"omp_num_threads: {omp_threads}\n")
        lf.write(f"hostname: {platform.node()}\n")
        lf.write(f"cpu_model: {cpu_model}\n")
        lf.write(f"python_version: {sys.version.split()[0]}\n")
        lf.write(f"timestamp: {timestamp}\n")
        lf.write(f"equations: {args.equations or 'all'}\n")
    print(f"Timing log → {log_path}")


if __name__ == "__main__":
    main()
