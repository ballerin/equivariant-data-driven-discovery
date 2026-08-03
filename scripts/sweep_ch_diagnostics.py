"""
Diagnostic sweep: CH noise robustness — polynomial backend, coarse d_tol.

Compares original vs modified multigraph libraries for weak-form (direct + IBP)
using only the polynomial stencil backend.  Sweeps over a coarse d_tol grid
to identify promising data-hyperparameter combinations quickly.

Usage
-----
    # Sweep only the new diagnostic datasets
    python scripts/sweep_ch_diagnostics.py

    # Include existing baseline CH datasets from data/
    python scripts/sweep_ch_diagnostics.py --data-dirs data/ data/ch_diagnostics/

    # Run a quick subset
    python scripts/sweep_ch_diagnostics.py --data-dirs data/ch_diagnostics/ --name quick
"""

import argparse
import csv
import datetime
import gc
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
    collect_needed_jet_monomials,
)
from equivariant_pde_find.equations import (  # noqa: E402
    true_allen_cahn, true_swift_hohenberg, true_ks,
    true_cahn_hilliard, true_pme, mask_pme, EQUATION_META,
    true_allen_cahn_mod_multigraph, true_swift_hohenberg_mod_multigraph,
    true_ks_mod_multigraph, true_cahn_hilliard_mod_multigraph,
    true_pme_mod_multigraph, normalize_eq_name,
)
from equivariant_pde_find.library import (  # noqa: E402
    build_multigraph_library,
    build_mod_multigraph_library,
)
from equivariant_pde_find.regression import TrainSTRidge  # noqa: E402
from equivariant_pde_find.simulations import make_regression_data  # noqa: E402
from equivariant_pde_find.library.weak_form import build_weak_rhs  # noqa: E402

RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyperparameter grids (polynomial-only, coarse d_tol) ────────────────────
LAM = 1e-6
D_TOLS = [1e-3, 1e-2, 0.1, 1.0, 10.0]           # coarse grid
BACKENDS = ["polynomial"]                           # polynomial only
from equivariant_pde_find.library.weak_form import (  # noqa: E402
    select_test_function_degree as _wf_select_degree,
    test_function_is_resolved as _wf_is_resolved,
)

WEAK_M = [4, 8, 12]
# Test-function degree follows Messenger & Bortz (2021), Algorithm 1, rather
# than a fixed constant: with a fixed p the fourth-order weak-IBP columns pick
# up a 3 % (m=12) to 84 % (m=8) systematic error.
# The direct weak form never differentiates psi, so it keeps the old degree.
WEAK_MAX_DERIV = 4
WEAK_P_DIRECT = 6


def weak_p(m, use_ibp=True):
    """Test-function degree for support *m* (see WEAK_MAX_DERIV above)."""
    if not use_ibp:
        return WEAK_P_DIRECT
    return _wf_select_degree(m, WEAK_MAX_DERIV)


def weak_ibp_supports(m_values):
    """Subset of *m_values* whose grid can represent psi^(WEAK_MAX_DERIV)."""
    return [m for m in m_values
            if _wf_is_resolved(m, weak_p(m), WEAK_MAX_DERIV)]
SAMPLE_STRIDE = 4
LIBRARIES = ["multigraph", "mod_multigraph"]

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
    """Load a .npz data file.  Extracts eq_name and sim_eps if present."""
    data = np.load(npz_path, allow_pickle=True)
    snapshots = data["snapshots"]
    times = data["times"]
    Lx = float(data["Lx"])
    Ly = float(data["Ly"])
    N = int(data["N"])
    obs_noise = float(data["obs_noise"])

    # eq_name: prefer .npz metadata, fall back to filename parsing
    if "eq_name" in data:
        eq_name = str(data["eq_name"])
    else:
        eq_name = os.path.basename(npz_path).rsplit("_", 1)[0]
        for suffix in ["_clean", "_noise2pct", "_noise5pct",
                        "_noise0p5pct", "_noise1pct", "_noise0pct"]:
            if eq_name.endswith(suffix):
                eq_name = eq_name[: -len(suffix)]
                break

    # sim_eps: for ch_diag datasets, the generating epsilon is stored
    sim_eps = None
    if "sim_eps" in data:
        sim_eps = float(data["sim_eps"])

    return snapshots, times, Lx, Ly, N, obs_noise, eq_name, sim_eps


def build_weak_form(f_snaps, ft_snaps, idx_sets, Lx, Ly, library, backend,
                    m, use_ibp):
    """Build weak-form (G, b) for a given library and backend."""
    with derivatives.use_backend(backend):
        rows_G, rows_b = [], []
        for s in range(len(f_snaps)):
            f = f_snaps[s]
            ft = ft_snaps[s]
            if library == "multigraph":
                G, desc = build_multigraph_library(
                    f, Lx, Ly, max_p=2, max_n=3,
                    formulation="weak", m=m, p=weak_p(m, use_ibp),
                    use_ibp=use_ibp,
                )
            else:  # mod_multigraph
                G, desc = build_mod_multigraph_library(
                    f, Lx, Ly, max_p=2, max_n=3,
                    formulation="weak", m=m, p=weak_p(m, use_ibp),
                    use_ibp=use_ibp,
                )
            b_w = build_weak_rhs(ft, m=m, p=weak_p(m, use_ibp))

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
        description="CH diagnostic sweep: polynomial backend, coarse d_tol."
    )
    parser.add_argument("--data-dirs", nargs="+",
                        default=[os.path.join(os.path.dirname(SCRIPT_DIR),
                                              "data", "ch_diagnostics")],
                        help="Directories to scan for .npz files.")
    parser.add_argument("--name", default="", help="Tag appended to output filename.")
    args = parser.parse_args()

    # Collect all .npz files across data dirs
    all_npz = []
    for d in args.data_dirs:
        if not os.path.isdir(d):
            print(f"WARNING: data dir not found: {d}", flush=True)
            continue
        for f in sorted(os.listdir(d)):
            if f.endswith(".npz") and not f.startswith("."):
                all_npz.append((d, f))

    if not all_npz:
        print("No .npz files found.", flush=True)
        return

    print(f"Found {len(all_npz)} .npz files across {len(args.data_dirs)} dirs.",
          flush=True)

    # ── output ───────────────────────────────────────────────────────────
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"_{args.name}" if args.name else ""
    out_path = os.path.join(RESULTS_DIR, f"ch_diag_{timestamp}{tag}.csv")
    print(f"Output → {out_path}", flush=True)

    fieldnames = [
        "equation", "noise", "method", "library", "backend",
        "weak_m", "lam", "d_tol",
        "n_active", "fp", "fn", "exact_recovery", "E_coeff", "cond_theta",
        "ch_eps", "dataset_file", "sim_N", "sim_L", "sim_T_end",
        "sim_dt_snap", "time_ms",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        f.flush()

        for data_dir, npz_file in all_npz:
            npz_path = os.path.join(data_dir, npz_file)
            print(f"\n{'='*70}")
            print(f"Dataset: {npz_file}")
            print(f"{'='*70}")

            snapshots, times, Lx, Ly, N, obs_noise, eq_name, sim_eps = \
                load_dataset(npz_path)
            eq_name = normalize_eq_name(eq_name)
            true_fn_key = _map_equivariant_key(eq_name)
            true_fn, skip_frac, mask_fn = _get_meta(eq_name)

            is_ch = eq_name == "ch"
            if not is_ch:
                print(f"  SKIP (not CH: {eq_name})", flush=True)
                continue

            # Determine CH epsilon: use sim_eps if stored, else default 0.2
            if sim_eps is not None:
                eps_list = [sim_eps]
            else:
                eps_list = [0.2]  # default for legacy baseline datasets

            # Extract metadata for CSV
            sim_N = N
            sim_L = Lx
            sim_T_end = str(times[-1]) if len(times) > 0 else ""
            sim_dt_snap = str(times[1] - times[0]) if len(times) > 1 else ""

            for ch_eps in eps_list:
                noise_label = obs_noise
                print(f"  eq={eq_name}  noise={noise_label:.4f}"
                      f"  ch_eps={ch_eps}"
                      f"  N={sim_N}  L={sim_L:.1f}  T_end={sim_T_end}"
                      f"  dt_snap={sim_dt_snap}", flush=True)

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

                # Build evaluation Theta (jet basis) on SUBSAMPLED points only.
                # build_polyjet_library returns (N², 105) per snap — for
                # N=256 that's 55 MB/snap × 199 snaps = 11 GB.  We index into
                # each snap's library with idx_sets_raw to keep only 200
                # points/snap, reducing this to ~33 MB total.
                rows_theta, rows_ft = [], []
                for s in range(len(f_snaps_raw)):
                    T_full, _ = build_polyjet_library(
                        f_snaps_raw[s], Lx, Ly, max_order=4, max_degree=3,
                    )
                    rows_theta.append(T_full[idx_sets_raw[s]])
                    rows_ft.append(ft_snaps_raw[s].ravel()[idx_sets_raw[s]])
                Tj = np.vstack(rows_theta)
                Ftj = np.concatenate(rows_ft)
                eval_thetas = {b: (Tj, Ftj) for b in BACKENDS}
                del rows_theta, rows_ft, Tj, Ftj
                # Dynamic jet descriptions from multigraph bases
                _, desc_mg_jet = build_multigraph_library(
                    f_snaps[0], Lx, Ly, max_p=2, max_n=3,
                )
                _, desc_mod_jet = build_mod_multigraph_library(
                    f_snaps[0], Lx, Ly, max_p=2, max_n=3,
                )
                desc_jet = collect_needed_jet_monomials(
                    desc_mg=desc_mg_jet, desc_mod_mg=desc_mod_jet,
                )

                # Build ground truth in jet basis for each library
                true_jet = {}
                for lib in LIBRARIES:
                    if lib == "multigraph":
                        _, desc_lib = build_multigraph_library(
                            f_snaps[0], Lx, Ly, max_p=2, max_n=3,
                            formulation="weak", m=4, p=weak_p(4, False), use_ibp=False,
                        )
                    else:
                        _, desc_lib = build_mod_multigraph_library(
                            f_snaps[0], Lx, Ly, max_p=2, max_n=3,
                            formulation="weak", m=4, p=weak_p(4, False), use_ibp=False,
                        )
                    if true_fn_key == "ch":
                        true_w = TRUTH_FNS[lib][true_fn_key](desc_lib, eps=ch_eps)
                    elif true_fn_key == "ac":
                        true_w = TRUTH_FNS[lib][true_fn_key](desc_lib, eps=0.15)
                    elif true_fn_key == "sh":
                        true_w = TRUTH_FNS[lib][true_fn_key](desc_lib, r=0.3)
                    else:
                        true_w = TRUTH_FNS[lib][true_fn_key](desc_lib)

                    if lib == "multigraph":
                        true_jet[lib] = map_multigraph_to_jet(
                            true_w, desc_lib, desc_jet)
                    else:
                        true_jet[lib] = map_mod_multigraph_to_jet(
                            true_w, desc_lib, desc_jet)

                # Build one (lib, m, ibp) at a time — run STRidge, then free memory.
                # This avoids storing all 12 prebuilt matrices simultaneously,
                # preventing OOM on large grids (e.g. N=256).
                for lib, backend, m_val, use_ibp in itertools.product(
                    LIBRARIES, BACKENDS, WEAK_M, [False, True],
                ):
                    # A support too small to represent psi^(WEAK_MAX_DERIV)
                    # cannot carry the IBP formulation at any degree; the
                    # direct weak form never differentiates psi and is fine.
                    if use_ibp and m_val not in weak_ibp_supports([m_val]):
                        print(f"    Skipping: lib={{lib}} m={{m_val}} ibp=True "
                              f"(2m+1 points cannot resolve a degree-"
                              f"{{weak_p(m_val)}} test function; ", flush=True)
                        continue

                    method_key = "weak_ibp" if use_ibp else "weak_direct"

                    t0 = time.time()
                    print(f"    Building: lib={lib} backend={backend} m={m_val}"
                          f" ibp={use_ibp}", flush=True)
                    Theta, Ft, desc = build_weak_form(
                        f_snaps, ft_snaps, idx_sets, Lx, Ly,
                        lib, backend, m_val, use_ibp,
                    )
                    cond_val = cond_number(Theta)
                    print(f"      shape={Theta.shape} cond={cond_val:.2e}"
                          f" ({time.time()-t0:.1f}s)", flush=True)

                    # Run STRidge for all d_tol values immediately
                    for d_tol in D_TOLS:
                        t_row = time.time()
                        w = TrainSTRidge(Theta, Ft, LAM, d_tol, maxit=25, seed=0)
                        eval_result = evaluate_unified(
                            w, lib, desc, desc_jet,
                            *eval_thetas[backend], true_jet[lib],
                        )
                        dt_ms = (time.time() - t_row) * 1000

                        row = dict(
                            equation=eq_name,
                            noise=noise_label,
                            method=method_key,
                            library=lib,
                            backend=backend,
                            weak_m=m_val,
                            lam=LAM,
                            d_tol=d_tol,
                            n_active=eval_result["n_active"],
                            fp=eval_result["fp"],
                            fn=eval_result["fn"],
                            exact_recovery=eval_result["exact_recovery"],
                            E_coeff=eval_result["E_coeff"],
                            cond_theta=cond_val,
                            ch_eps=ch_eps,
                            dataset_file=npz_file,
                            sim_N=sim_N,
                            sim_L=sim_L,
                            sim_T_end=sim_T_end,
                            sim_dt_snap=sim_dt_snap,
                            time_ms=f"{dt_ms:.1f}",
                        )
                        writer.writerow(row)
                        f.flush()

                    # Free memory before next combination
                    del Theta, Ft, desc

                # Free all per-dataset large arrays to prevent OOM accumulation
                del (eval_thetas, true_jet, desc_jet,
                     f_snaps_raw, ft_snaps_raw, idx_sets_raw,
                     f_snaps, ft_snaps, idx_sets, snapshots)
                gc.collect()

    print(f"\nDone. Results → {out_path}", flush=True)


if __name__ == "__main__":
    main()
