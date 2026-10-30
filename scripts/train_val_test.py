#!/usr/bin/env python3
"""Hyperparameter selection on held-out trajectories, then a frozen test.

``sweep.py`` answers "for this dataset, fit every setting and score each
against the truth".  That makes every number an oracle: the setting is chosen
knowing the answer.  This script answers the honest question instead — "choose
a setting without looking at the truth, then see how often it works on data it
has never seen" — by running the three tiers of new-plan.md §2 and §6:

    train       for each setting, fit on each train trajectory -> w
    validation  score every w on held-out trajectories -> E_val
    select      pick a setting from E_val alone (no ground truth)
    test        refit at the selected setting, once per test trajectory,
                and score against the truth

**The unit of independence is the trajectory** (§4.1).  Rows are never split:
neighbouring frames are near-duplicates, a weak row centred at t integrates
t−m_t … t+m_t, and spectral derivatives make every row in a frame depend on the
whole frame.  Trajectories differ by initial-condition seed and carry
independent noise draws, so a seed split leaks nothing.

**What transfers from train to test is the setting, not the coefficients.**
Each test trajectory is refit from scratch at the frozen setting, which is
what a practitioner does with a fresh dataset, and gives one independent
recovery outcome per test trajectory — the error bars a single-seed sweep
cannot produce.

**The M&B test-function rule is not a hyperparameter.**  It is a deterministic
function of the data, so each trajectory legitimately gets its own
(m, p, m_t, p_t) and nothing needs freezing.  Coefficients are in PDE units,
so a ``w`` fitted where the rule gave m_t=20 scores correctly on a trajectory
where it gave 17 (checked: the true coefficient vector gives residuals
0.0025 / 0.0024 / 0.0031 across three ks_cell trajectories).

**This script applies no selection rule.**  It writes every (setting, train
seed, validation seed) residual and lets the summariser apply the sparsest-
within-(1+δ) rule, the 1-SE rule, or anything else, to the same numbers.  The
one thing it must decide is which setting to carry into the test tier, so it
does that with ``--rule`` and records the choice.

Centering is a **compared factor**, not a hyperparameter (§7.2): selecting over
it would be a second oracle.  Each arm therefore gets its own selection.

Usage
-----
    # all equations and noise levels in a split dataset directory
    python scripts/train_val_test.py --data-dir data/paper_split \\
        --out-dir results/tvt_20260919

    # one cell, as the runner does it
    python scripts/train_val_test.py --data-dir data/paper_split \\
        --equations ks_cell --noise 0.02 --out-dir results/tvt_20260919
"""

import argparse
import csv
import itertools
import os
import sys
import time
from collections import defaultdict

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src")))
sys.path.insert(0, SCRIPT_DIR)

import sweep as S                                                # noqa: E402
from equivariant_pde_find.eval_setup import prepare_evaluation   # noqa: E402
from equivariant_pde_find.equations import (                     # noqa: E402
    REPRESENTABLE_LIBRARIES, normalize_eq_name,
)
from equivariant_pde_find.evaluation import evaluate_unified     # noqa: E402
from equivariant_pde_find.library import weak_form as wf         # noqa: E402
from equivariant_pde_find.regression import SolverCache          # noqa: E402

NOISE_LABELS = {0.0: "clean", 0.02: "noise2pct", 0.05: "noise5pct"}
SPLITS = ("train", "val", "test")


# ── hyperparameter grid ─────────────────────────────────────────────────────


def settings_grid(algorithms, lams, d_tols, rfe_lams, mstls_alphas):
    """Every hyperparameter setting, as (algorithm, lam, d_tol, alpha) tuples.

    The unused slots are None rather than "" so that the tuple can be used as
    a dict key and still round-trip through the CSV unambiguously.
    """
    out = []
    if "stridge" in algorithms:
        out += [("stridge", lam, d_tol, None)
                for lam, d_tol in itertools.product(lams, d_tols)]
    if "rfe" in algorithms:
        out += [("rfe", lam, None, None) for lam in rfe_lams]
    if "mstls" in algorithms:
        out += [("mstls", None, None, a) for a in mstls_alphas]
    return out


def fit_setting(Theta, Ft, desc, setting, center, center_tol, cache,
                mstls_M=None):
    """Fit one setting.  Returns (n_active, w, lam_out).

    ``sweep_one`` takes a ``true_eq`` argument it never reads; None is passed
    explicitly so that it is obvious no ground truth reaches the train or
    validation tiers.
    """
    algorithm, lam, d_tol, alpha = setting
    return S.sweep_one(
        Theta, Ft, None,
        0.0 if lam is None else lam,
        0.0 if d_tol is None else d_tol,
        algorithm,
        mstls_alpha=0.5 if alpha is None else alpha,
        mstls_M=mstls_M,
        center=center, center_tol=center_tol, desc=desc, cache=cache,
    )


# ── one trajectory ──────────────────────────────────────────────────────────


def load_trajectory(path, args):
    """Load a dataset and return everything the tiers need from it."""
    (snapshots, times, Lx, Ly, N, dt_snap, obs_noise, eq_name,
     _ptype, _pdir) = S.load_dataset(path)
    eq_name = normalize_eq_name(eq_name)
    true_fn, skip_frac, mask_fn = S._get_meta(eq_name)
    n_skip = max(0, int(skip_frac * len(snapshots)))
    rule = wf.mb_test_function_params(
        snapshots[n_skip:], spatial_axes=range(1, snapshots.ndim),
        time_axis=0, alpha_x=S.ALPHA_X, alpha_t=S.ALPHA_T,
        tau=args.tau, tau_hat=args.tau_hat,
    )
    return dict(
        path=path, snapshots=snapshots, times=times, Lx=Lx, Ly=Ly,
        eq_name=eq_name, true_fn=true_fn, skip_frac=skip_frac,
        mask_fn=mask_fn, n_skip=n_skip, rule=rule, obs_noise=obs_noise,
    )


def build_design(traj, library, args):
    """Build the space-time weak system for one trajectory and library."""
    r = traj["rule"]
    cfg = dict(m=r["m"], p=r["p"], m_t=r["m_t"], p_t=r["p_t"],
               time_ibp=(args.time_ibp == "on"),
               space_boundary=args.space_boundary)
    with wf.use_weak_boundary(args.space_boundary):
        Theta, Ft, desc = S.build_weak_in_time(
            traj["snapshots"][traj["n_skip"]:], traj["times"][traj["n_skip"]:],
            traj["Lx"], traj["Ly"], library, args.backend, args.method,
            cfg, args.space_stride, args.time_stride, mask_fn=traj["mask_fn"],
        )
    return Theta, Ft, desc, cfg


def relative_residual(Theta, Ft, w):
    """``‖Θw − b‖ / ‖b‖`` — the held-out residual the selection rule reads."""
    return float(np.linalg.norm(Theta @ w - Ft)
                 / max(np.linalg.norm(Ft), 1e-14))


# ── the three tiers ─────────────────────────────────────────────────────────


def run_cell(eq, noise_label, files, args, sel_writer, test_writer):
    """Train, validate, select and test one (equation, noise) cell."""
    libraries = sorted(set(args.libraries) & REPRESENTABLE_LIBRARIES[eq]) \
        if args.representable_only else list(args.libraries)
    skipped = sorted(set(args.libraries) - set(libraries))
    print(f"\n{'=' * 70}\n{eq}  {noise_label}   libraries: "
          f"{' '.join(libraries)}"
          + (f"   (not representable: {' '.join(skipped)})" if skipped else "")
          + f"\n{'=' * 70}", flush=True)

    grid = settings_grid(args.algorithms, args.lams, args.d_tols,
                         args.rfe_lams, args.mstls_alphas)
    center_modes = [(m, m == "on") for m in args.centering]
    print(f"  {len(grid)} settings x {len(center_modes)} centering arm(s) "
          f"x {len(files['train'])} train trajectories", flush=True)

    for library in libraries:
        # ── TRAIN: fit every setting on every train trajectory ───────────
        # Only the coefficient vectors are kept, so one design matrix is
        # held at a time -- at 256 x 256 each is several hundred MB.
        fits = {}          # (arm, train_seed, setting) -> w
        for seed, path in files["train"]:
            t0 = time.time()
            traj = load_trajectory(path, args)
            Theta, Ft, desc, cfg = build_design(traj, library, args)
            cache = SolverCache()
            mstls_M = (wf.compute_m_scale(Theta, Ft)
                       if library in ("wsindy", "wsindy_proj") else None)
            for arm, center in center_modes:
                for setting in grid:
                    n_active, w, _ = fit_setting(
                        Theta, Ft, desc, setting, center,
                        args.centering_tol, cache, mstls_M)
                    fits[(arm, seed, setting)] = (w, n_active)
            print(f"  [{library}] train s{seed}: {Theta.shape} "
                  f"m={cfg['m']} m_t={cfg['m_t']} "
                  f"({time.time() - t0:.0f}s)", flush=True)
            del Theta, Ft, cache

        # ── VALIDATION: score every fit on every validation trajectory ───
        # No ground truth is read here; E_val is a pure held-out residual.
        residuals = defaultdict(list)     # (arm, setting) -> [E_val, ...]
        for seed, path in files["val"]:
            t0 = time.time()
            traj = load_trajectory(path, args)
            Theta, Ft, desc, cfg = build_design(traj, library, args)
            for (arm, train_seed, setting), (w, n_active) in fits.items():
                e_val = relative_residual(Theta, Ft, w)
                residuals[(arm, setting)].append(e_val)
                sel_writer.writerow(dict(
                    equation=eq, noise=noise_label, library=library,
                    centering=arm, algorithm=setting[0],
                    lam=_fmt(setting[1]), d_tol=_fmt(setting[2]),
                    mstls_alpha=_fmt(setting[3]),
                    train_seed=train_seed, val_seed=seed,
                    n_active=n_active, E_val=f"{e_val:.8f}",
                    space_m=cfg["m"], time_m=cfg["m_t"],
                ))
            print(f"  [{library}] val   s{seed}: {Theta.shape} "
                  f"m={cfg['m']} m_t={cfg['m_t']} "
                  f"({time.time() - t0:.0f}s)", flush=True)
            del Theta, Ft

        # ── SELECT: one setting per (arm, solver), truth-free ────────────
        chosen = {}
        for arm, _ in center_modes:
            for algorithm in args.algorithms:
                cand = [(s, residuals[(arm, s)],
                         fits[(arm, files["train"][0][0], s)][1])
                        for s in grid if s[0] == algorithm
                        and residuals.get((arm, s))]
                if not cand:
                    continue
                chosen[(arm, algorithm)] = select(cand, args.rule, args.delta)
                s, mean_e, _ = chosen[(arm, algorithm)]
                print(f"  [{library}] select {algorithm:8s} centering={arm}: "
                      f"lam={_fmt(s[1])} d_tol={_fmt(s[2])} "
                      f"alpha={_fmt(s[3])}  E_val={mean_e:.5f}", flush=True)

        # ── TEST: refit at the frozen setting, score against truth ───────
        for seed, path in files["test"]:
            t0 = time.time()
            traj = load_trajectory(path, args)
            clean_path = _clean_counterpart(path)
            snaps_clean = (np.load(clean_path)["snapshots"]
                           if traj["obs_noise"] > 0 else traj["snapshots"])
            ev = prepare_evaluation(
                traj["snapshots"], snaps_clean, traj["times"],
                traj["Lx"], traj["Ly"], eq, traj["true_fn"],
                traj["skip_frac"], mask_fn=traj["mask_fn"],
                eval_mode=args.eval_mode, verbose=False,
            )
            Theta, Ft, desc, cfg = build_design(traj, library, args)
            cache = SolverCache()
            mstls_M = (wf.compute_m_scale(Theta, Ft)
                       if library in ("wsindy", "wsindy_proj") else None)
            const_cols = S.const_column_reporter(Theta, Ft, desc)
            for (arm, algorithm), (setting, mean_e, se_e) in chosen.items():
                center = arm == "on"
                n_active, w, lam_out = fit_setting(
                    Theta, Ft, desc, setting, center, args.centering_tol,
                    cache, mstls_M)
                res = evaluate_unified(
                    w, library, desc, ev.desc_jet, ev.Theta_eval, ev.Ft_eval,
                    ev.true_jet, M_ws_to_jet=ev.M_ws_to_jet,
                )
                row = dict(
                    equation=eq, noise=noise_label, library=library,
                    centering=arm, algorithm=algorithm, test_seed=seed,
                    sel_rule=args.rule,
                    lam=_fmt(setting[1]), d_tol=_fmt(setting[2]),
                    mstls_alpha=_fmt(setting[3]),
                    E_val_selected=f"{mean_e:.8f}",
                    E_val_se=("" if se_e is None else f"{se_e:.8f}"),
                    n_active=n_active, fp=res["fp"], fn=res["fn"],
                    exact_recovery=res["exact_recovery"],
                    E_coeff=res["E_coeff"],
                    max_coeff_err=res["max_coeff_err"],
                    E_data=f"{relative_residual(Theta, Ft, w):.6f}",
                    space_m=cfg["m"], time_m=cfg["m_t"],
                    n_samples=res["n_samples"],
                )
                row.update(const_cols(w))
                test_writer.writerow(row)
            print(f"  [{library}] test  s{seed}: {Theta.shape} "
                  f"({time.time() - t0:.0f}s)", flush=True)
            del Theta, Ft, cache, ev


def select(candidates, rule, delta):
    """Pick one setting from validation residuals alone.

    ``candidates`` is [(setting, [E_val, ...], n_active), ...].  Returns
    ``(setting, mean_E_val, se_E_val)``.  Both rules take the *sparsest* model
    whose residual clears a floor: sparsity is the objective and the residual
    is the constraint, never the other way round -- minimising the residual
    picks the densest candidate every time (§5.2).
    """
    stats = []
    for setting, evs, n_active in candidates:
        evs = [e for e in evs if np.isfinite(e)]
        if not evs:
            continue
        mean = float(np.mean(evs))
        se = (float(np.std(evs, ddof=1) / np.sqrt(len(evs)))
              if len(evs) > 1 else None)
        stats.append((setting, mean, se, n_active))
    if not stats:
        return None
    best = min(s[1] for s in stats)
    if rule == "one-se":
        # Threshold is one standard error above the best mean, using that
        # setting's own SE.  With a single train x val pair there is no SE,
        # so this degenerates to argmin; the caller is warned at startup.
        best_se = next((s[2] for s in stats if s[1] == best), None) or 0.0
        cutoff = best + best_se
    else:
        cutoff = best * (1.0 + delta)
    keep = [s for s in stats if s[1] <= cutoff] or stats
    setting, mean, se, _ = min(
        keep, key=lambda s: (s[3], s[1],
                             float("inf") if s[0][1] is None else s[0][1]))
    return setting, mean, se


# ── plumbing ────────────────────────────────────────────────────────────────


def _fmt(v):
    return "" if v is None else v


def _clean_counterpart(path):
    return (path.replace("_noise2pct", "_clean")
                .replace("_noise5pct", "_clean"))


def discover(data_dir, equations, noises):
    """Group the split manifest into {(eq, noise_label): {split: [(seed, path)]}}."""
    manifest = os.path.join(data_dir, "manifest.csv")
    if not os.path.exists(manifest):
        raise SystemExit(
            f"No manifest at {manifest}.  Generate a split dataset with "
            f"`generate_data.py --profile paper --splits N_TRAIN N_VAL N_TEST`."
        )
    cells = defaultdict(lambda: {s: [] for s in SPLITS})
    for r in csv.DictReader(open(manifest)):
        if r["split"] not in SPLITS:
            continue
        eq = normalize_eq_name(r["equation"])
        if equations and eq not in equations:
            continue
        if noises and r["noise_label"] not in noises:
            continue
        cells[(eq, r["noise_label"])][r["split"]].append(
            (int(r["ic_seed"]), r["path"]))
    for cell in cells.values():
        for s in SPLITS:
            cell[s].sort()
    missing = {k: [s for s in SPLITS if not v[s]] for k, v in cells.items()}
    for k, m in missing.items():
        if m:
            raise SystemExit(f"{k}: no {'/'.join(m)} trajectories in manifest.")
    return cells


SEL_FIELDS = ["equation", "noise", "library", "centering", "algorithm",
              "lam", "d_tol", "mstls_alpha", "train_seed", "val_seed",
              "n_active", "E_val", "space_m", "time_m"]
TEST_FIELDS = ["equation", "noise", "library", "centering", "algorithm",
               "test_seed", "sel_rule", "lam", "d_tol", "mstls_alpha",
               "E_val_selected", "E_val_se", "n_active", "fp", "fn",
               "exact_recovery", "E_coeff", "max_coeff_err", "E_data",
               "space_m", "time_m", "n_samples", "const_active",
               "intercept_rel"]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True,
                    help="Split dataset directory (must carry a manifest.csv "
                         "with train/val/test labels).")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--equations", nargs="+", default=None)
    ap.add_argument("--noise", nargs="+", type=float, default=None)
    ap.add_argument("--libraries", nargs="+",
                    default=["multigraph", "mod_multigraph", "pdefind",
                             "polyjet", "wsindy"])
    ap.add_argument("--representable-only", action="store_true", default=True)
    ap.add_argument("--all-libraries", dest="representable_only",
                    action="store_false",
                    help="Keep libraries that cannot represent the equation.")
    ap.add_argument("--method", default="spacetime_weak_ibp")
    ap.add_argument("--backend", default="polynomial")
    ap.add_argument("--space-boundary", default="periodic",
                    choices=["valid", "periodic"])
    ap.add_argument("--time-ibp", default="on", choices=["on", "off"])
    ap.add_argument("--centering", nargs="+", default=["off", "on"],
                    choices=["off", "on"])
    ap.add_argument("--centering-tol", type=float, default=0.05)
    ap.add_argument("--algorithms", nargs="+",
                    default=["stridge", "rfe", "mstls"],
                    choices=["stridge", "rfe", "mstls"])
    ap.add_argument("--lams", nargs="+", type=float, default=S.LAMS)
    ap.add_argument("--d-tols", nargs="+", type=float, default=S.D_TOLS)
    ap.add_argument("--rfe-lams", nargs="+", type=float,
                    default=[0.0, 1e-6, 1e-4, 1e-2])
    ap.add_argument("--mstls-alphas", nargs="+", type=float,
                    default=S.MSTLS_ALPHAS)
    ap.add_argument("--rule", default="delta", choices=["delta", "one-se"],
                    help="Which rule picks the setting carried into the test "
                         "tier.  Every residual is written out either way, so "
                         "the summariser can apply both.")
    ap.add_argument("--delta", type=float, default=0.5)
    ap.add_argument("--space-stride", type=int, default=S.SPACE_STRIDE)
    ap.add_argument("--time-stride", type=int, default=1)
    ap.add_argument("--tau", type=float, default=S.TAU)
    ap.add_argument("--tau-hat", type=float, default=S.TAU_HAT)
    ap.add_argument("--eval-mode", default="qr", choices=["dense", "qr"])
    args = ap.parse_args()

    noises = ({NOISE_LABELS[n] for n in args.noise} if args.noise else None)
    cells = discover(args.data_dir, set(args.equations or []), noises)
    if not cells:
        raise SystemExit("No (equation, noise) cells matched.")
    os.makedirs(args.out_dir, exist_ok=True)

    n_train = len(next(iter(cells.values()))["train"])
    if args.rule == "one-se" and n_train < 2:
        print("[warn] --rule one-se needs at least 2 train trajectories for a "
              "standard error; with 1 it degenerates to argmin E_val. "
              "Falling back to the delta rule for the test tier.", flush=True)
        args.rule = "delta"

    t_start = time.time()
    for (eq, noise_label), files in sorted(cells.items()):
        tag = f"{eq}_{noise_label}"
        sel_path = os.path.join(args.out_dir, f"selection_{tag}.csv")
        test_path = os.path.join(args.out_dir, f"test_{tag}.csv")
        with open(sel_path, "w", newline="") as fs, \
                open(test_path, "w", newline="") as ft:
            sel_writer = csv.DictWriter(fs, fieldnames=SEL_FIELDS)
            test_writer = csv.DictWriter(ft, fieldnames=TEST_FIELDS)
            sel_writer.writeheader()
            test_writer.writeheader()
            run_cell(eq, noise_label, files, args, sel_writer, test_writer)
        print(f"  -> {os.path.basename(sel_path)}, "
              f"{os.path.basename(test_path)}", flush=True)
    print(f"\nDone in {time.time() - t_start:.0f}s -> {args.out_dir}")


if __name__ == "__main__":
    main()
