"""
Anisotropic-noise ablation sweep: weak-form only, STRidge + RFE only.

Compares equivariant vs non-equivariant libraries under space-dependent
(anisotropic) observation noise at 2 % and 5 % levels.  Clean data is
loaded for evaluation only — no clean sweep rows.

Restrictions (vs the canonical ``sweep.py``):
    - Methods: ``weak_direct`` and ``weak_ibp`` only (no pointwise)
    - Algorithms: ``stridge`` and ``rfe`` only (no MSTLS, brute, oracle)
    - Noise: 2 % and 5 % anisotropic datasets only
    - Libraries: only those whose column space contains the ground truth
      for each equation (``REPRESENTABLE_LIBRARIES``)

Usage
-----
    # Quick test: Allen-Cahn, multigraph only, spectral backend
    python scripts/sweep_aniso.py --equations ac --libraries multigraph --backends spectral

    # Full sweep (all equations, all representable libraries, both backends)
    python scripts/sweep_aniso.py

    # Restrict to specific libraries and weak-m values
    python scripts/sweep_aniso.py --libraries multigraph,pdefind --weak-m 8,12

    # With custom RFE ridge regularisation
    python scripts/sweep_aniso.py --rfe-lam 0.0,0.01,0.1
"""

import argparse
import csv
import datetime
import itertools
import os
import platform
import re
import subprocess
import sys
import time

import numpy as np

# ── path setup ──────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives  # noqa: E402
from equivariant_pde_find.evaluation import (  # noqa: E402
    build_eval_theta_jet, evaluate_unified,
)
from equivariant_pde_find.library.conversions import (  # noqa: E402
    build_wsindy_to_jet_map, _MULTIGRAPH_TO_JET_DECOMPOSITION,
    collect_needed_jet_monomials,
)
from equivariant_pde_find.equations import (  # noqa: E402
    normalize_eq_name, EQUATION_META, CH_EPSILON,
)
from equivariant_pde_find.library import (  # noqa: E402
    build_multigraph_library, build_pdefind_library,
    build_polyjet_library, build_wsindy_library,
    build_mod_multigraph_library,
)
from equivariant_pde_find.library import weak_form as wf  # noqa: E402
from equivariant_pde_find.regression import TrainSTRidge, TrainRFE  # noqa: E402
from equivariant_pde_find.simulations import make_regression_data  # noqa: E402
from equivariant_pde_find.utils import spectral_lowpass  # noqa: E402

RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
os.makedirs(RESULTS_DIR, exist_ok=True)

# ── hyperparameter grids ────────────────────────────────────────────────────
LAMS = [1e-6, 1e-4, 1e-2, 1.0]            # STRidge ridge penalty
RFE_LAMS_DEFAULT = [0.0]                    # RFE ridge (0 = pure OLS)
D_TOLS = [1e-3, 1e-2, 0.1, 1.0, 10.0]     # STRidge threshold sweep
WEAK_M = [4, 8, 12]                         # test-function half-supports
WEAK_P = 6                                   # test-function degree (direct & IBP)
SAMPLE_STRIDE = 4                            # weak-form query-point stride

# Both multigraph bases are always swept together.  Keys are the canonical
# short equation name as returned by ``normalize_eq_name``.
REPRESENTABLE_LIBRARIES = {
    "ac":  {"multigraph", "mod_multigraph", "pdefind", "wsindy", "polyjet"},
    "sh":  {"multigraph", "mod_multigraph", "pdefind", "wsindy", "polyjet"},
    "ch":  {"multigraph", "mod_multigraph", "wsindy", "polyjet"},
    "pme": {"multigraph", "mod_multigraph", "wsindy", "polyjet"},
    "ks_cell":   {"multigraph", "mod_multigraph", "polyjet"},
    "ks_wchaos": {"multigraph", "mod_multigraph", "polyjet"},
}
# KS regimes have no anisotropic datasets.
ALL_LIBRARIES = ["multigraph", "mod_multigraph", "pdefind", "polyjet", "wsindy"]

# ── equation-name→prefix mapping (for dataset discovery) ────────────────────
_EQ_TO_FILENAME = {
    "ac":       {"allen_cahn", "ac"},
    "sh":       {"swift_hohenberg", "sh"},
    "ch":       {"cahn_hilliard", "ch"},
    "pme":      {"porous_medium", "pme"},
    "ks_cell":   {"ks_cellular", "ks_cell"},
    "ks_wchaos": {"ks_weak_chaotic", "ks_wchaos"},
}


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers (adapted from sweep.py to keep this script self-contained)
# ═══════════════════════════════════════════════════════════════════════════════

def load_dataset(npz_path):
    """Load a .npz data file.

    Returns
    -------
    snapshots, times, Lx, Ly, N, dt_snap, obs_noise, eq_name,
    profile_type, profile_direction
    """
    data = np.load(npz_path)
    snapshots = data["snapshots"]
    times = data["times"]
    Lx = float(data["Lx"])
    Ly = float(data["Ly"])
    N = int(data["N"])
    dt_snap = float(data["dt_snap"])
    obs_noise = float(data["obs_noise"])

    # eq_name from metadata (added 2026-06-30), fallback to filename parsing
    if "eq_name" in data:
        eq_name = str(data["eq_name"])
    else:
        eq_name = os.path.basename(npz_path).rsplit("_", 1)[0]
        for suffix in ["_clean", "_noise2pct", "_noise5pct"]:
            if eq_name.endswith(suffix):
                eq_name = eq_name[: -len(suffix)]
                break

    profile_type = str(data["profile_type"]) if "profile_type" in data else ""
    profile_direction = str(data["profile_direction"]) if "profile_direction" in data else ""

    return (snapshots, times, Lx, Ly, N, dt_snap, obs_noise, eq_name,
            profile_type, profile_direction)


def _get_meta(name):
    """Return (true_fn, skip_fraction, mask_fn) for a dataset name."""
    if name in EQUATION_META:
        return EQUATION_META[name]
    for prefix in sorted(EQUATION_META, key=len, reverse=True):
        if name.startswith(prefix):
            return EQUATION_META[prefix]
    raise KeyError(f"No metadata for equation: {name}")


def weak_p(m, use_ibp=True):
    """Test-function degree for support *m*.

    Returns a fixed p=6 for both direct and IBP weak forms.  See the
    discussion in sweep.py for the quadrature-accuracy vs sparse-recovery
    trade-off.
    """
    return WEAK_P


def build_weak_form(f_snaps, ft_snaps, idx_sets, Lx, Ly, library, backend,
                    m, use_ibp):
    """Build weak-form (G, b) for a given library and backend."""
    p = weak_p(m, use_ibp)
    with derivatives.use_backend(backend):
        rows_G, rows_b = [], []
        for s in range(len(f_snaps)):
            f = f_snaps[s]
            ft = ft_snaps[s]
            if library == "multigraph":
                G, desc = build_multigraph_library(
                    f, Lx, Ly, max_p=2, max_n=3,
                    formulation="weak", m=m, p=p,
                    use_ibp=use_ibp,
                )
            elif library == "mod_multigraph":
                G, desc = build_mod_multigraph_library(
                    f, Lx, Ly, max_p=2, max_n=3,
                    formulation="weak", m=m, p=p,
                    use_ibp=use_ibp,
                )
            elif library == "polyjet":
                G, desc = build_polyjet_library(
                    f, Lx, Ly, max_order=4, max_degree=3,
                    formulation="weak", m=m, p=p, use_ibp=use_ibp,
                )
            elif library == "wsindy":
                G, desc = build_wsindy_library(
                    f, Lx, Ly, max_order=4, max_poly=3,
                    formulation="weak", m=m, p=p, use_ibp=use_ibp,
                )
            else:
                G, desc = build_pdefind_library(
                    f, Lx, Ly, max_order=4, max_poly=3,
                    formulation="weak", m=m, p=p, use_ibp=use_ibp,
                )
            b_w = wf.build_weak_rhs(ft, m=m, p=p)

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


def sweep_one(Theta, Ft, lam, d_tol, algorithm):
    """Run regression and return (n_active, w, lam_out).

    Only ``"stridge"`` and ``"rfe"`` are supported (no MSTLS).
    """
    if algorithm == "rfe":
        w = TrainRFE(Theta, Ft, lam=lam)
        lam_out = lam
    else:
        w = TrainSTRidge(Theta, Ft, lam, d_tol)
        lam_out = lam
    n_active = int(np.count_nonzero(np.abs(w) > 1e-10))
    return n_active, w, lam_out


def cond_number(Theta):
    """Condition number (2-norm) of Theta."""
    return float(np.linalg.cond(Theta))


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Anisotropic-noise ablation: weak-form only, STRidge + RFE."
    )
    parser.add_argument(
        "--data-dir", default=DATA_DIR,
        help="Directory containing .npz data files.",
    )
    parser.add_argument(
        "--equations", nargs="+", default=None,
        help="Which equations to sweep (default: all with aniso data).",
    )
    parser.add_argument(
        "--libraries", nargs="+",
        default=ALL_LIBRARIES,
        choices=ALL_LIBRARIES,
        help="Libraries to sweep (intersected with representable set per eq).",
    )
    parser.add_argument(
        "--noise", nargs="+", type=float,
        default=[0.02, 0.05],
        help="Noise levels to sweep (default: 0.02 0.05).",
    )
    parser.add_argument(
        "--backends", nargs="+",
        default=["spectral", "polynomial"],
        choices=["spectral", "polynomial"],
    )
    parser.add_argument(
        "--weak-m", nargs="+", type=int,
        default=WEAK_M,
        help="Test-function half-support values.",
    )
    parser.add_argument(
        "--rfe-lam", nargs="+", type=float,
        default=RFE_LAMS_DEFAULT,
        help="Ridge regularisation for RFE (0 = pure OLS).",
    )
    parser.add_argument(
        "--name", default="",
        help="Optional tag appended to output filename.",
    )
    parser.add_argument(
        "--out", default=None,
        help="Exact output CSV path (overrides auto-generated name).",
    )
    args = parser.parse_args()

    # ── discover anisotropic datasets ────────────────────────────────────
    data_dir = args.data_dir
    all_npz = sorted([
        f for f in os.listdir(data_dir)
        if f.endswith(".npz") and "aniso" in f and not f.startswith(".")
    ])

    if args.equations:
        eq_set = set(args.equations)
        eq_prefixes = set()
        for eq in eq_set:
            eq_prefixes |= _EQ_TO_FILENAME.get(eq, {eq})
        all_npz = [f for f in all_npz
                   if any(f.startswith(eq + "_") for eq in eq_prefixes)]

    # Filter by noise level
    if args.noise is not None:
        noise_labels = {0.0: "clean", 0.02: "noise2pct", 0.05: "noise5pct"}
        allowed = {noise_labels[n] for n in args.noise}
        all_npz = [f for f in all_npz
                   if any(f.endswith("_" + nl + ".npz") for nl in allowed)]

    if not all_npz:
        print(f"No anisotropic .npz files found in {data_dir} matching filters.")
        return

    backends = args.backends
    rfe_lams = args.rfe_lam
    weak_m_vals = args.weak_m

    # Filter out m values that can't resolve the degree-p test function for IBP
    usable_ibp_m_vals = [
        m for m in weak_m_vals
        if wf.test_function_is_resolved(m, WEAK_P, max_deriv=4)
    ]
    dropped = [m for m in weak_m_vals if m not in usable_ibp_m_vals]
    if dropped:
        print(
            f"  [weak_ibp] skipping m={dropped}: 2m+1 grid points cannot "
            f"resolve the degree-{WEAK_P} test function for IBP"
        )

    # ── output ───────────────────────────────────────────────────────────
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    if args.out:
        out_path = args.out
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    else:
        tag = f"_{args.name}" if args.name else ""
        out_path = os.path.join(RESULTS_DIR, f"sweep_aniso_{timestamp}{tag}.csv")
    print(f"Output → {out_path}", flush=True)

    fieldnames = [
        "equation", "noise", "method", "algorithm", "library", "backend",
        "weak_m", "lam", "d_tol", "mstls_alpha",
        "n_active", "fp", "fn", "exact_recovery", "E_coeff", "max_coeff_err",
        "E_data", "cond_theta", "n_samples",
        "post_hoc_noise", "spectral_cutoff",
        "profile_type", "profile_direction",
        "time_ms",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        f.flush()

        t_start_total = time.time()
        total_combos = 0
        eq_summaries = []

        for npz_file in all_npz:
            npz_path = os.path.join(data_dir, npz_file)
            print(f"\n{'='*70}")
            print(f"Dataset: {npz_file}")
            print(f"{'='*70}")

            # ── Load ─────────────────────────────────────────────────────
            snapshots_orig, times, Lx, Ly, N, dt_snap, obs_noise, eq_name, \
                profile_type, profile_direction = load_dataset(npz_path)
            eq_name = normalize_eq_name(eq_name)
            true_fn, skip_frac, mask_fn = _get_meta(eq_name)

            # ── Filter libraries by representability ─────────────────────
            representable = REPRESENTABLE_LIBRARIES.get(eq_name, set(ALL_LIBRARIES))
            lib_set = sorted(set(args.libraries) & representable,
                             key=lambda x: ALL_LIBRARIES.index(x)
                             if x in ALL_LIBRARIES else 99)
            if not lib_set:
                print(f"  No representable libraries for {eq_name}. "
                      f"Skipping.", flush=True)
                continue

            # ── Load clean counterpart for evaluation (mandatory) ───────
            # Evaluation always uses pointwise spectral derivatives on clean
            # data, regardless of the regression method or dataset noise level.
            clean_path = os.path.join(
                os.path.dirname(npz_path),
                f"{eq_name}_clean.npz",
            )
            if not os.path.exists(clean_path):
                raise FileNotFoundError(
                    f"Clean evaluation data required but not found: "
                    f"{clean_path}.  Generate it with "
                    f"`python scripts/generate_data.py --equations "
                    f"{eq_name}`."
                )
            snapshots_clean_eval = np.load(clean_path)["snapshots"]
            print(f"    clean eval: {os.path.basename(clean_path)}")

            noise_label = obs_noise
            print(f"  eq={eq_name}  noise={noise_label:.3f}  "
                  f"profile={profile_type}/{profile_direction}  "
                  f"libraries={lib_set}", flush=True)

            # ── Build ground truth (backend-independent) ─────────────────
            _, desc_eq_ref = build_multigraph_library(
                snapshots_orig[0], Lx, Ly,
            )
            # CH and SH take extra parameters
            is_ch = eq_name == "ch"
            if is_ch:
                true_eq = true_fn(desc_eq_ref, eps=CH_EPSILON)
            elif eq_name == "sh":
                true_eq = true_fn(desc_eq_ref, r=0.3)
            elif eq_name == "ac":
                true_eq = true_fn(desc_eq_ref, eps=0.15)
            else:
                true_eq = true_fn(desc_eq_ref)

            # ── Prepare regression data ──────────────────────────────────
            f_snaps_raw, ft_snaps_raw, idx_sets_raw, _ = \
                make_regression_data(
                    snapshots_orig, times, n_per_snap=200,
                    skip_fraction=skip_frac, seed=0,
                )

            # Apply mask if needed (PME)
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
            print(f"    interior snaps: {n_snaps}  "
                  f"points: {n_snaps * 200}", flush=True)

            # ── Jet descriptions & ground truth (backend-independent) ────
            _, desc_pdefind = build_pdefind_library(
                snapshots_orig[0], Lx, Ly, max_order=4, max_poly=3,
            )
            _, desc_ws = build_wsindy_library(
                snapshots_orig[0], Lx, Ly, max_order=4, max_poly=3,
            )
            desc_jet = collect_needed_jet_monomials(
                desc_mg=desc_eq_ref, desc_pdefind=desc_pdefind,
                desc_ws=desc_ws,
            )
            # Ground truth in jet basis
            true_jet = np.zeros(len(desc_jet))
            for j, (name, coeff) in enumerate(zip(desc_eq_ref, true_eq)):
                if abs(coeff) < 1e-12:
                    continue
                if name not in _MULTIGRAPH_TO_JET_DECOMPOSITION:
                    continue
                for jet_name, mult in _MULTIGRAPH_TO_JET_DECOMPOSITION[name]:
                    true_jet[desc_jet.index(jet_name)] += coeff * mult

            M_ws_to_jet = build_wsindy_to_jet_map(desc_ws, desc_jet)

            # ── Build evaluation Theta (pointwise, spectral) ────────────
            # Always use pointwise formulation + spectral derivatives on
            # clean data.  The polynomial jet basis is the universal
            # comparison space.
            Theta_eval, Ft_eval = build_eval_theta_jet(
                snapshots_clean_eval, times, Lx, Ly,
                skip_fraction=skip_frac, seed=0, backend="spectral",
                jet_monomials=desc_jet,
            )
            n_active_true = int(np.sum(np.abs(true_jet) > 1e-10))
            print(f"    eval Theta: {Theta_eval.shape}  "
                  f"active_true={n_active_true}", flush=True)

            # ── Pre-build weak-form matrices ─────────────────────────────
            prebuilt = {}

            for library, backend in itertools.product(lib_set, backends):
                for m_val in weak_m_vals:
                    # weak_ibp: only usable m values
                    if m_val not in usable_ibp_m_vals:
                        continue
                    for use_ibp_flag in (False, True):
                        method = "weak_ibp" if use_ibp_flag else "weak_direct"
                        actual_key = (method, library, backend, m_val)
                        if actual_key in prebuilt:
                            continue

                        t0 = time.time()
                        print(f"    Building: method={method} "
                              f"lib={library} backend={backend} m={m_val}",
                              flush=True)

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

            # ── Sweep hyperparams ────────────────────────────────────────
            eq_exact_count = 0
            eq_total = 0

            for library, backend in itertools.product(lib_set, backends):
                for m_val in weak_m_vals:
                    if m_val not in usable_ibp_m_vals:
                        continue
                    for use_ibp_flag in (False, True):
                        method = "weak_ibp" if use_ibp_flag else "weak_direct"
                        actual_key = (method, library, backend, m_val)
                        if actual_key not in prebuilt:
                            continue
                        Theta, Ft, desc, cond_val = prebuilt[actual_key]

                        # ── STRidge: lam × d_tol grid ────────────────────
                        for lam in LAMS:
                            for d_tol in D_TOLS:
                                t_row = time.time()
                                n_active, w, lam_out = sweep_one(
                                    Theta, Ft, lam, d_tol, "stridge",
                                )

                                eval_result = evaluate_unified(
                                    w, library, desc, desc_jet,
                                    Theta_eval, Ft_eval, true_jet,
                                    M_ws_to_jet=M_ws_to_jet,
                                )
                                dt_ms = (time.time() - t_row) * 1000

                                row = dict(
                                    equation=eq_name,
                                    noise=noise_label,
                                    method=method,
                                    algorithm="stridge",
                                    library=library,
                                    backend=backend,
                                    weak_m=m_val,
                                    lam=lam_out,
                                    d_tol=d_tol,
                                    mstls_alpha="",
                                    n_active=n_active,
                                    fp=eval_result["fp"],
                                    fn=eval_result["fn"],
                                    exact_recovery=eval_result["exact_recovery"],
                                    E_coeff=eval_result["E_coeff"],
                                    max_coeff_err=eval_result["max_coeff_err"],
                                    E_data=f"{np.linalg.norm(Theta @ w - Ft) / max(np.linalg.norm(Ft), 1e-14):.6f}",
                                    cond_theta=cond_val,
                                    n_samples=eval_result["n_samples"],
                                    post_hoc_noise="",
                                    spectral_cutoff="",
                                    profile_type=profile_type,
                                    profile_direction=profile_direction,
                                    time_ms=f"{dt_ms:.1f}",
                                )
                                writer.writerow(row)
                                total_combos += 1
                                eq_total += 1
                                if eval_result["exact_recovery"]:
                                    eq_exact_count += 1

                        # ── RFE: lam grid ────────────────────────────────
                        for rfe_lam in rfe_lams:
                            t_row = time.time()
                            n_active, w, lam_out = sweep_one(
                                Theta, Ft, rfe_lam, 0.0, "rfe",
                            )

                            eval_result = evaluate_unified(
                                w, library, desc, desc_jet,
                                Theta_eval, Ft_eval, true_jet,
                                M_ws_to_jet=M_ws_to_jet,
                            )
                            dt_ms = (time.time() - t_row) * 1000

                            row = dict(
                                equation=eq_name,
                                noise=noise_label,
                                method=method,
                                algorithm="rfe",
                                library=library,
                                backend=backend,
                                weak_m=m_val,
                                lam=lam_out,
                                d_tol="",
                                mstls_alpha="",
                                n_active=n_active,
                                fp=eval_result["fp"],
                                fn=eval_result["fn"],
                                exact_recovery=eval_result["exact_recovery"],
                                E_coeff=eval_result["E_coeff"],
                                max_coeff_err=eval_result["max_coeff_err"],
                                E_data=f"{np.linalg.norm(Theta @ w - Ft) / max(np.linalg.norm(Ft), 1e-14):.6f}",
                                cond_theta=cond_val,
                                n_samples=eval_result["n_samples"],
                                post_hoc_noise="",
                                spectral_cutoff="",
                                profile_type=profile_type,
                                profile_direction=profile_direction,
                                time_ms=f"{dt_ms:.1f}",
                            )
                            writer.writerow(row)
                            total_combos += 1
                            eq_total += 1

                    f.flush()

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
        lf.write(f"# Anisotropic sweep timing log\n")
        lf.write(f"csv_file: {os.path.basename(out_path)}\n")
        lf.write(f"total_wall_time_s: {elapsed_total:.1f}\n")
        lf.write(f"total_wall_time_min: {elapsed_total/60:.1f}\n")
        lf.write(f"total_rows: {total_combos}\n")
        lf.write(f"omp_num_threads: {omp_threads}\n")
        lf.write(f"hostname: {platform.node()}\n")
        lf.write(f"cpu_model: {cpu_model}\n")
        lf.write(f"python_version: {sys.version.split()[0]}\n")
        lf.write(f"timestamp: {timestamp}\n")
        lf.write(f"equations: {args.equations or 'all aniso'}\n")
    print(f"Timing log → {log_path}")


if __name__ == "__main__":
    main()
