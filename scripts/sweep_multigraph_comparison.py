"""
Focused sweep: compare original vs modified multigraph libraries.

Uses ONLY weak-form methods (weak_direct and weak_ibp), runs STRidge for
both libraries on all benchmark equations, and writes a CSV for later
summary generation.

Usage
-----
    python scripts/sweep_multigraph_comparison.py
    python scripts/sweep_multigraph_comparison.py --equations swift_hohenberg cahn_hilliard
    python scripts/sweep_multigraph_comparison.py --noise 0.02 --name quick
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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives  # noqa: E402
from equivariant_pde_find.evaluation import (  # noqa: E402
    build_eval_theta_jet, evaluate_unified,
)
from equivariant_pde_find.library.conversions import (  # noqa: E402
    map_multigraph_to_jet, map_mod_multigraph_to_jet,
)
from equivariant_pde_find.equations import (  # noqa: E402
    true_allen_cahn, true_swift_hohenberg, true_ks,
    true_cahn_hilliard, true_pme, mask_pme, EQUATION_META,
    true_allen_cahn_mod_multigraph, true_swift_hohenberg_mod_multigraph,
    true_ks_mod_multigraph, true_cahn_hilliard_mod_multigraph,
    true_pme_mod_multigraph, normalize_eq_name, CH_EPSILON,
)
from equivariant_pde_find.library import (  # noqa: E402
    build_polyjet_library, build_multigraph_library,
    build_mod_multigraph_library,
)
from equivariant_pde_find.regression import TrainSTRidge  # noqa: E402
from equivariant_pde_find.simulations import make_regression_data  # noqa: E402
from equivariant_pde_find.library.weak_form import build_weak_rhs  # noqa: E402

RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyperparameter grids ────────────────────────────────────────────────────
LAMS = [1e-6]
D_TOLS = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
BACKENDS = ["spectral", "polynomial"]
WEAK_M = [4, 8, 12]
WEAK_P = 6
SAMPLE_STRIDE = 4

# ── truth functions per library and equation ─────────────────────────────────
TRUTH_FNS = {
    "multigraph": {
        "ac": true_allen_cahn,
        "sh": true_swift_hohenberg,
        "ks": true_ks,
        "ch": true_cahn_hilliard,
        "pme": true_pme,
    },
    "mod_multigraph": {
        "ac": true_allen_cahn_mod_multigraph,
        "sh": true_swift_hohenberg_mod_multigraph,
        "ks": true_ks_mod_multigraph,
        "ch": true_cahn_hilliard_mod_multigraph,
        "pme": true_pme_mod_multigraph,
    },
}


def _get_meta(name):
    if name in EQUATION_META:
        return EQUATION_META[name]
    for prefix in sorted(EQUATION_META, key=len, reverse=True):
        if name.startswith(prefix):
            return EQUATION_META[prefix]
    raise KeyError(f"No metadata for equation: {name}")


def load_dataset(npz_path):
    data = np.load(npz_path)
    snapshots = data["snapshots"]
    times = data["times"]
    Lx = float(data["Lx"])
    Ly = float(data["Ly"])
    N = int(data["N"])
    obs_noise = float(data["obs_noise"])
    eq_name = os.path.basename(npz_path).rsplit("_", 1)[0]
    for suffix in ["_clean", "_noise2pct", "_noise5pct", "_noise1pct", "_noise0p5pct", "_noise0p1pct"]:
        if eq_name.endswith(suffix):
            eq_name = eq_name[: -len(suffix)]
            break
    return snapshots, times, Lx, Ly, N, obs_noise, eq_name


def build_weak_form(f_snaps, ft_snaps, idx_sets, Lx, Ly, library, backend,
                    m, use_ibp):
    with derivatives.use_backend(backend):
        rows_G, rows_b = [], []
        for s in range(len(f_snaps)):
            f = f_snaps[s]
            ft = ft_snaps[s]
            if library == "multigraph":
                G, desc = build_multigraph_library(
                    f, Lx, Ly, max_p=2, max_n=3,
                    formulation="weak", m=m, p=WEAK_P, use_ibp=use_ibp,
                )
            else:  # mod_multigraph
                G, desc = build_mod_multigraph_library(
                    f, Lx, Ly, max_p=2, max_n=3,
                    formulation="weak", m=m, p=WEAK_P, use_ibp=use_ibp,
                )
            b_w = build_weak_rhs(ft, m=m, p=WEAK_P)

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


def cond_number(Theta):
    return float(np.linalg.cond(Theta))


def _map_equivariant_key(eq_name):
    """Map canonical equation name to truth-function key."""
    name = normalize_eq_name(eq_name)
    if name in ("ac", "sh", "ch", "pme"):
        return name
    if name in ("ks_cell", "ks_wchaos", "ks_chaos"):
        return "ks"
    raise KeyError(f"Unknown equation: {eq_name}")


def main():
    parser = argparse.ArgumentParser(
        description="Sweep: compare original vs modified multigraph (weak-form only)."
    )
    parser.add_argument("--data-dir", default=DATA_DIR)
    parser.add_argument("--equations", nargs="+", default=None)
    parser.add_argument("--noise", nargs="+", type=float, default=None)
    parser.add_argument("--backends", nargs="+", default=BACKENDS, choices=BACKENDS)
    parser.add_argument("--weak-m", nargs="+", type=int, default=WEAK_M)
    parser.add_argument("--name", default="", help="Tag appended to output filename.")
    args = parser.parse_args()

    # Discover datasets
    all_npz = sorted([
        f for f in os.listdir(args.data_dir)
        if f.endswith(".npz") and not f.startswith(".")
    ])
    if args.equations:
        eq_set = set(args.equations)
        all_npz = [f for f in all_npz
                   if any(f.startswith(eq + "_") or f.startswith(eq) for eq in eq_set)]
    if args.noise is not None:
        noise_labels = {0.0: "clean", 0.02: "noise2pct", 0.05: "noise5pct"}
        allowed = {noise_labels[n] for n in args.noise}
        all_npz = [f for f in all_npz
                   if any(f.endswith("_" + nl + ".npz") for nl in allowed)]

    if not all_npz:
        print(f"No .npz files found in {args.data_dir} matching filters.")
        return

    libraries = ["multigraph", "mod_multigraph"]
    weak_m_vals = args.weak_m

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"_{args.name}" if args.name else ""
    out_path = os.path.join(RESULTS_DIR, f"mg_comparison_{timestamp}{tag}.csv")
    print(f"Output → {out_path}", flush=True)

    fieldnames = [
        "equation", "noise", "method", "library", "backend",
        "weak_m", "lam", "d_tol",
        "n_active", "fp", "fn", "exact_recovery", "E", "cond_theta",
        "time_ms",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        f.flush()

        total_combos = 0

        for npz_file in all_npz:
            npz_path = os.path.join(args.data_dir, npz_file)
            print(f"\n{'='*70}")
            print(f"Dataset: {npz_file}")
            print(f"{'='*70}")

            snapshots, times, Lx, Ly, N, obs_noise, eq_name = load_dataset(npz_path)
            eq_name = normalize_eq_name(eq_name)
            true_fn_key = _map_equivariant_key(eq_name)
            true_fn, skip_frac, mask_fn = _get_meta(eq_name)

            is_ch = eq_name == "ch"
            eps_list = [CH_EPSILON] if is_ch else [None]

            for ch_eps in eps_list:
                noise_label = obs_noise
                print(f"  eq={eq_name}  noise={noise_label:.4f}", flush=True)

                # Prepare regression data
                f_snaps_raw, ft_snaps_raw, idx_sets_raw, _ = make_regression_data(
                    snapshots, times, n_per_snap=200,
                    skip_fraction=skip_frac, seed=0,
                )
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

                # Build per-backend evaluation Theta (jet basis)
                eval_thetas = {}
                for b in args.backends:
                    Tj, Ftj = build_eval_theta_jet(
                        snapshots, times, Lx, Ly,
                        skip_fraction=skip_frac, seed=0, backend=b,
                    )
                    eval_thetas[b] = (Tj, Ftj)
                _, desc_jet = build_polyjet_library(
                    snapshots[0], Lx, Ly, max_order=4, max_degree=3,
                )

                # Build ground truth in jet basis for each library
                true_jet = {}
                for lib in libraries:
                    # Get the correct desc for this lib (cheap: just build on one snapshot)
                    if lib == "multigraph":
                        _, desc_lib = build_multigraph_library(
                            f_snaps[0], Lx, Ly, max_p=2, max_n=3,
                            formulation="weak", m=4, p=WEAK_P, use_ibp=False)
                    else:
                        _, desc_lib = build_mod_multigraph_library(
                            f_snaps[0], Lx, Ly, max_p=2, max_n=3,
                            formulation="weak", m=4, p=WEAK_P, use_ibp=False)
                    if is_ch:
                        true_w = TRUTH_FNS[lib][true_fn_key](desc_lib, eps=ch_eps)
                    elif true_fn_key == "swift_hohenberg":
                        true_w = TRUTH_FNS[lib][true_fn_key](desc_lib, r=0.3)
                    elif true_fn_key == "allen_cahn":
                        true_w = TRUTH_FNS[lib][true_fn_key](desc_lib, eps=0.15)
                    else:
                        true_w = TRUTH_FNS[lib][true_fn_key](desc_lib)
                    if lib == "multigraph":
                        true_jet[lib] = map_multigraph_to_jet(true_w, desc_lib, desc_jet)
                    else:
                        true_jet[lib] = map_mod_multigraph_to_jet(true_w, desc_lib, desc_jet)

                # Pre-build all (library, backend, m, ibp) matrices
                prebuilt = {}
                for lib, backend, m_val, use_ibp in itertools.product(
                    libraries, args.backends, weak_m_vals, [False, True],
                ):
                    method_key = "weak_ibp" if use_ibp else "weak_direct"
                    actual_key = (method_key, lib, backend, m_val)

                    t0 = time.time()
                    print(f"    Building: lib={lib} backend={backend} m={m_val}"
                          f" ibp={use_ibp}", flush=True)
                    Theta, Ft, desc = build_weak_form(
                        f_snaps, ft_snaps, idx_sets, Lx, Ly,
                        lib, backend, m_val, use_ibp,
                    )
                    cond_val = cond_number(Theta)
                    prebuilt[actual_key] = (Theta, Ft, desc, cond_val)
                    print(f"      shape={Theta.shape} cond={cond_val:.2e}"
                          f" ({time.time()-t0:.1f}s)", flush=True)

                # Sweep STRidge over all d_tol values
                for (method_key, lib, backend, m_val), \
                    (Theta, Ft, desc, cond_val) in prebuilt.items():
                    for d_tol in tqdm(D_TOLS, desc=f"  {lib} m={m_val} {method_key}", leave=False):
                        t_row = time.time()
                        w = TrainSTRidge(Theta, Ft, 1e-6, d_tol, maxit=25, seed=0)
                        eval_result = evaluate_unified(
                            w, lib, desc, desc_jet,
                            *eval_thetas[backend], true_jet[lib],
                        )
                        dt_ms = (time.time() - t_row) * 1000

                        # Debug: if result is missing expected keys
                        if "E" not in eval_result:
                            print(f"WARNING: eval_result missing 'E' key. "
                                  f"Keys: {sorted(eval_result.keys())}. "
                                  f"lib={lib} eq={eq_name} noise={noise_label}",
                                  flush=True)

                        row = dict(
                            equation=eq_name,
                            noise=noise_label,
                            method=method_key,
                            library=lib,
                            backend=backend,
                            weak_m=m_val,
                            lam=1e-6,
                            d_tol=d_tol,
                            n_active=eval_result.get("n_active", 0),
                            fp=eval_result.get("fp", float("nan")),
                            fn=eval_result.get("fn", float("nan")),
                            exact_recovery=eval_result.get("exact_recovery", 0),
                            E=eval_result.get("E_coeff", eval_result.get("E", float("nan"))),
                            cond_theta=cond_val,
                            time_ms=f"{dt_ms:.1f}",
                        )
                        writer.writerow(row)
                        total_combos += 1

    print(f"\nDone. {total_combos} total rows → {out_path}")


if __name__ == "__main__":
    main()
