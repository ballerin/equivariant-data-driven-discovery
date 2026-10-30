"""
Shared machinery for the paper notebooks.

Every notebook follows the same seven-section script: generate the clean and
noisy datasets with the *exact* parameters of ``scripts/generate_data.py``,
draw the motion, then run the sweep-optimal configuration for the best
equivariant library and the best non-equivariant one.

The point of this module is that the notebooks must reproduce the numbers in
``results/sweep_paper_v2_*.csv`` rather than approximate them, so the design
matrices are built by calling ``scripts/sweep.py`` directly (same ``SPACE_P_FIXED``,
same ``SPACE_STRIDE``, same query-point RNG) and the metrics come from the
same polynomial-jet evaluation space.

This lives in ``notebooks/`` rather than in the package because it imports
matplotlib, and the package depends only on numpy.
"""

import os
import sys

# ── path setup (must precede the package imports below) ─────────────────────

NB_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(NB_DIR, ".."))
FIG_DIR = os.path.join(CODE_DIR, "results", "figures")
DATA_DIR = os.path.join(CODE_DIR, "data")

for _rel in ("src", "scripts"):
    _p = os.path.join(CODE_DIR, _rel)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np                                            # noqa: E402
import matplotlib.pyplot as plt                               # noqa: E402

import sweep                                                  # noqa: E402
from equivariant_pde_find.library import weak_form as wf      # noqa: E402
from equivariant_pde_find import equations                    # noqa: E402
from equivariant_pde_find.equations import (                  # noqa: E402
    EQUATION_META, normalize_eq_name,
)
from equivariant_pde_find.evaluation import (                 # noqa: E402
    build_eval_theta_jet, evaluate_unified,
)
from equivariant_pde_find.library import (                    # noqa: E402
    build_multigraph_library, build_mod_multigraph_library,
    build_pdefind_library, build_wsindy_library,
)
from equivariant_pde_find.library.conversions import (        # noqa: E402
    build_wsindy_to_jet_map, collect_needed_jet_monomials,
    map_multigraph_to_jet, map_mod_multigraph_to_jet,
    map_standard_to_jet, map_wsindy_to_jet,
    _MULTIGRAPH_TO_JET_DECOMPOSITION,
)
from equivariant_pde_find.regression import (                 # noqa: E402
    TrainSTRidge, TrainRFE, TrainMSTLS,
)
from equivariant_pde_find.simulations import (                # noqa: E402
    make_regression_data,
)
from equivariant_pde_find.utils import print_pde              # noqa: E402

N_PER_SNAP = 200          # matches sweep.py and generate_data.py
EQUIVARIANT = ("multigraph", "mod_multigraph")

LIBRARY_LABEL = {
    "multigraph": "Multigraph",
    "mod_multigraph": "Modified multigraph",
    "pdefind": "PDE-FIND",
    "wsindy": "WSINDy",
    "polyjet": "PolyJet",
}
METHOD_LABEL = {
    "pointwise": "pointwise",
    "space_weak_direct": "space-weak (direct)",
    "space_weak_ibp": "space-weak (IBP)",
}
ALGORITHM_LABEL = {"stridge": "STRidge", "rfe": "RFE", "mstls": "MSTLS"}


# ═══════════════════════════════════════════════════════════════════════════
# Figures
# ═══════════════════════════════════════════════════════════════════════════

def _color_limits(snapshots, symmetric=True, percentile=99.5):
    """Robust colour limits over the whole run.

    A hard-coded [-1, 1] clips Swift-Hohenberg and Kuramoto-Sivashinsky
    badly, and a full min/max is hostage to a single outlier pixel, so take
    a high percentile of |f| (symmetric fields) or of f (positive fields).
    """
    if symmetric:
        v = float(np.percentile(np.abs(snapshots), percentile))
        return -v, v
    lo = float(np.percentile(snapshots, 100.0 - percentile))
    hi = float(np.percentile(snapshots, percentile))
    return lo, hi


def _save(fig, save_as):
    """Write a figure as both PDF (for the paper) and PNG (for previews)."""
    if save_as is None:
        return []
    os.makedirs(FIG_DIR, exist_ok=True)
    stem = os.path.join(FIG_DIR, save_as)
    paths = [f"{stem}.pdf", f"{stem}.png"]
    fig.savefig(paths[0], bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight", dpi=200)
    return paths


def plot_motion(snapshots, times, L, *, n_frames=5, cmap="RdBu_r",
                symmetric=True, percentile=99.5, save_as=None, title=None,
                frames=None, detrend=False, time_fmt=".3g"):
    """Draw the evolution as a row of evenly spaced snapshots.

    Parameters
    ----------
    snapshots : ndarray (n_saved, N, N)
    times     : ndarray (n_saved,)
    L         : float — domain side length, used for the axis extent
    n_frames  : int — how many frames to show (ignored if *frames* given)
    frames    : sequence of int or None — explicit frame indices
    cmap      : str — diverging for signed fields, sequential for positive
    symmetric : bool — centre the colour scale on zero
    detrend   : bool — subtract each frame's spatial mean before plotting.
        Needed for Kuramoto-Sivashinsky, whose mean drifts far outside the
        pattern amplitude; the right-hand side is unchanged by ``f → f + c``,
        so the drift carries no information about the dynamics and would
        otherwise saturate the colour scale.
    time_fmt  : str — format spec for the frame labels.  The integrators save
        the state *after* each step, so the first snapshot sits at ``t = dt``
        rather than ``t = 0``; ``".0f"`` rounds those one-step offsets away
        (0.01 → 0, 2.01 → 2) where the run is long enough for that to be the
        honest reading.
    save_as   : str or None — basename under ``results/figures/``

    Returns
    -------
    fig, list of written paths
    """
    if detrend:
        snapshots = snapshots - snapshots.mean(
            axis=tuple(range(1, snapshots.ndim)), keepdims=True,
        )
    if frames is None:
        frames = np.linspace(0, len(snapshots) - 1, n_frames).astype(int)
    frames = list(frames)
    vmin, vmax = _color_limits(snapshots, symmetric, percentile)

    fig, axes = plt.subplots(
        1, len(frames), figsize=(2.15 * len(frames), 2.5),
    )
    axes = np.atleast_1d(axes)
    for ax, idx in zip(axes, frames):
        im = ax.imshow(snapshots[idx], cmap=cmap, vmin=vmin, vmax=vmax,
                       origin="lower", extent=[0, L, 0, L],
                       interpolation="bilinear")
        ax.set_title(f"$t = {times[idx]:{time_fmt}}$", fontsize=11, pad=6)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.6)
            spine.set_color("0.5")

    fig.subplots_adjust(bottom=0.20, wspace=0.08, top=0.86)
    cbar_ax = fig.add_axes([0.33, 0.10, 0.34, 0.045])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.outline.set_linewidth(0.6)
    cbar.ax.tick_params(labelsize=9, length=3)
    if title:
        fig.suptitle(title, fontsize=12, y=1.02)

    return fig, _save(fig, save_as)


def plot_radial_profiles(snapshots, times, L, *, frames=None, n_frames=5,
                         save_as=None, center=None, time_fmt=".3g"):
    """Azimuthally averaged profile f(r) at several times.

    For a radially symmetric spreading solution the heatmap shows the shape
    but hides the motion — the front advances by a few pixels per frame.  The
    profile makes the advancing front and the decaying peak explicit.
    """
    if frames is None:
        frames = np.linspace(0, len(snapshots) - 1, n_frames).astype(int)
    frames = list(frames)

    N = snapshots.shape[1]
    h = L / N
    if center is None:
        center = ((N - 1) / 2.0, (N - 1) / 2.0)
    yy, xx = np.indices((N, snapshots.shape[2]))
    r = np.hypot(xx - center[1], yy - center[0]) * h
    nbins = N // 2
    edges = np.linspace(0, r.max(), nbins + 1)
    which = np.clip(np.digitize(r.ravel(), edges) - 1, 0, nbins - 1)
    counts = np.bincount(which, minlength=nbins)
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    colors = plt.cm.viridis(np.linspace(0.05, 0.85, len(frames)))
    for color, idx in zip(colors, frames):
        prof = np.bincount(which, weights=snapshots[idx].ravel(),
                           minlength=nbins) / np.maximum(counts, 1)
        ax.plot(centers, prof, color=color, lw=1.6,
                label=f"$t = {times[idx]:{time_fmt}}$")
    ax.set_xlabel("$r$")
    ax.set_ylabel("$f$")
    ax.set_xlim(0, centers[-1])
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig, _save(fig, save_as)


def plot_clean_vs_noisy(clean, noisy, times, L, *, frame=-1, cmap="RdBu_r",
                        symmetric=True, percentile=99.5, noise_level=0.02,
                        save_as=None, detrend=False):
    """Show one frame clean, noisy, and their difference.

    This is what the regression at 2 % noise is actually up against, and it
    is far more informative than quoting the noise level alone.  *detrend*
    removes each frame's spatial mean, as in :func:`plot_motion`; the
    difference panel is unaffected by it.
    """
    if detrend:
        axes_sp = tuple(range(1, clean.ndim))
        clean = clean - clean.mean(axis=axes_sp, keepdims=True)
        noisy = noisy - noisy.mean(axis=axes_sp, keepdims=True)
    vmin, vmax = _color_limits(clean, symmetric, percentile)
    diff = noisy[frame] - clean[frame]
    dv = float(np.max(np.abs(diff))) or 1.0

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6))
    panels = [
        (clean[frame], cmap, vmin, vmax, "clean", False),
        (noisy[frame], cmap, vmin, vmax, f"+ {noise_level * 100:.0f}% noise",
         True),
        (diff, "PuOr_r", -dv, dv, "difference", True),
    ]
    for ax, (field, cm, lo, hi, label, bar) in zip(axes, panels):
        im = ax.imshow(field, cmap=cm, vmin=lo, vmax=hi, origin="lower",
                       extent=[0, L, 0, L], interpolation="bilinear")
        ax.set_title(label, fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.6)
            spine.set_color("0.5")
        if bar:  # first two panels share a scale — one colourbar is enough
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03).ax.tick_params(
                labelsize=8, length=3,
            )
    fig.suptitle(f"$t = {times[frame]:.3g}$", fontsize=11)
    fig.tight_layout()
    return fig, _save(fig, save_as)


# ═══════════════════════════════════════════════════════════════════════════
# Data and design matrices
# ═══════════════════════════════════════════════════════════════════════════

def check_matches_dataset(snapshots, eq, noise_label):
    """Assert the freshly simulated data equals the stored sweep dataset.

    The notebooks re-simulate rather than load, so this is the guarantee that
    they are regressing on exactly the arrays the sweep consumed.  Returns
    None (and says so) when the .npz is absent.
    """
    path = os.path.join(DATA_DIR, f"{eq}_{noise_label}.npz")
    if not os.path.exists(path):
        print(f"  (no stored dataset at data/{eq}_{noise_label}.npz — skipped)")
        return None
    stored = np.load(path)["snapshots"]
    err = float(np.max(np.abs(snapshots - stored)))
    status = "identical" if err < 1e-12 else "DIFFERS"
    print(f"  vs data/{eq}_{noise_label}.npz: max|Δ| = {err:.2e}  → {status}")
    return err


def regression_data(snapshots, times, eq, n_per_snap=N_PER_SNAP, seed=0):
    """Sub-sampled (f, f_t) pairs, with the sweep's per-equation conventions.

    Applies the equation's ``skip_fraction`` (transient discard) and mask
    from ``EQUATION_META`` — for the porous-medium equation the mask keeps
    the sample points inside the support, mirroring ``sweep.py``.
    """
    eq = normalize_eq_name(eq)
    _, skip_frac, mask_fn = EQUATION_META[eq]

    f_snaps, ft_snaps, idx_sets, dt_snap = make_regression_data(
        snapshots, times, n_per_snap=n_per_snap,
        skip_fraction=skip_frac, seed=seed,
    )

    if mask_fn is not None:
        # Same masking block as sweep.py, including its RNG seed, so the
        # query points match row for row.
        rng = np.random.default_rng(1)
        masked = []
        for s in range(len(f_snaps)):
            valid = np.where(mask_fn(f_snaps[s]))[0]
            if len(valid) <= n_per_snap:
                masked.append(valid)
            else:
                masked.append(np.sort(
                    rng.choice(valid, size=n_per_snap, replace=False)))
        idx_sets = masked

    return f_snaps, ft_snaps, idx_sets, dt_snap


def build_design(f_snaps, ft_snaps, idx_sets, Lx, Ly, cfg):
    """Design matrix (Theta, Ft) for one frozen configuration.

    Delegates to ``sweep.build_pointwise`` / ``sweep.build_space_weak`` so the
    result is identical to the sweep's, including the weak test-function
    degree ``sweep.SPACE_P_FIXED`` and query-point stride ``sweep.SPACE_STRIDE``.
    """
    method = cfg["method"]
    if method == "pointwise":
        return sweep.build_pointwise(
            f_snaps, ft_snaps, idx_sets, Lx, Ly, cfg["library"], cfg["backend"],
        )
    return sweep.build_space_weak(
        f_snaps, ft_snaps, idx_sets, Lx, Ly, cfg["library"], cfg["backend"],
        cfg["m"], space_ibp=(method == "space_weak_ibp"),
    )


def run_fit(Theta, Ft, cfg):
    """Run the configuration's sparse regression. Returns the coefficients."""
    algorithm = cfg["algorithm"]
    if algorithm == "stridge":
        return TrainSTRidge(Theta, Ft, cfg["lam"], cfg["d_tol"])
    if algorithm == "rfe":
        return TrainRFE(Theta, Ft, lam=cfg["lam"])
    if algorithm == "mstls":
        # sweep.py applies M-scaling to WSINDy only; the other libraries have
        # column norms spanning orders of magnitude and are left unscaled.
        M = (wf.compute_m_scale(Theta, Ft)
             if cfg["library"] == "wsindy" else None)
        w, _, _ = TrainMSTLS(
            Theta, Ft, lambdas=None, M=M, n_lambdas=sweep.MSTLS_N_LAMBDAS,
            maxit=sweep.MSTLS_MAXITS, alpha=cfg["alpha"],
        )
        return w
    raise ValueError(f"Unknown algorithm: {algorithm}")


def describe(cfg):
    """One-line human-readable summary of a frozen configuration."""
    bits = [
        LIBRARY_LABEL[cfg["library"]],
        METHOD_LABEL[cfg["method"]],
        ALGORITHM_LABEL[cfg["algorithm"]],
        f"{cfg['backend']} derivatives",
    ]
    if cfg["method"] != "pointwise":
        bits.append(f"m={cfg['m']}, p={sweep.SPACE_P_FIXED}")
    if cfg["algorithm"] == "stridge":
        bits.append(f"λ={cfg['lam']:g}, d_tol={cfg['d_tol']:g}")
    elif cfg["algorithm"] == "rfe":
        bits.append(f"λ={cfg['lam']:g}")
    elif cfg["algorithm"] == "mstls":
        bits.append(f"α={cfg['alpha']:g}")
    return " · ".join(bits)


# ═══════════════════════════════════════════════════════════════════════════
# Evaluation in the common polynomial-jet space
# ═══════════════════════════════════════════════════════════════════════════

class JetEvaluator:
    """Scores any library's recovery against the truth in one common basis.

    The libraries do not share a basis — the multigraph library has 20
    rotation-invariant columns, PDE-FIND has 60 coordinate ones — so fp/fn
    counts are not comparable term by term.  Following ``sweep.py``, every
    recovery is mapped analytically into the polynomial-jet basis, where the
    ground truth has a single unambiguous representation, and scored there.

    The evaluation library is always built with **spectral** derivatives on
    **clean** data: it defines truth, not the method under test.
    """

    def __init__(self, snapshots_clean, times, Lx, Ly, eq, snapshot_stride=1,
                 truth_kwargs=None):
        eq = normalize_eq_name(eq)
        true_fn, skip_frac, _ = EQUATION_META[eq]
        self.eq = eq
        self.snapshot_stride = snapshot_stride

        # Ground truth in the multigraph basis (the reference form)
        _, self.desc_mg = build_multigraph_library(snapshots_clean[0], Lx, Ly)
        self.true_mg = true_fn(self.desc_mg, **(truth_kwargs or {}))

        # Column names of the other libraries (backend-independent)
        _, self.desc_pdefind = build_pdefind_library(
            snapshots_clean[0], Lx, Ly, max_order=4, max_poly=3,
        )
        _, self.desc_wsindy = build_wsindy_library(
            snapshots_clean[0], Lx, Ly, max_order=4, max_poly=3,
        )
        self.desc_jet = collect_needed_jet_monomials(
            desc_mg=self.desc_mg, desc_pdefind=self.desc_pdefind,
            desc_ws=self.desc_wsindy,
        )

        # Ground truth in jet space — analytic, no derivatives involved
        self.true_jet = np.zeros(len(self.desc_jet))
        for name, coeff in zip(self.desc_mg, self.true_mg):
            if abs(coeff) < 1e-12:
                continue
            for jet_name, mult in _MULTIGRAPH_TO_JET_DECOMPOSITION[name]:
                self.true_jet[self.desc_jet.index(jet_name)] += coeff * mult

        self.M_ws_to_jet = build_wsindy_to_jet_map(
            self.desc_wsindy, self.desc_jet,
        )

        # Evaluation library on every grid point of every (strided) snapshot
        self.Theta, self.Ft = build_eval_theta_jet(
            snapshots_clean[::snapshot_stride], times[::snapshot_stride],
            Lx, Ly, skip_fraction=skip_frac, seed=0, backend="spectral",
            jet_monomials=self.desc_jet,
        )

    def to_jet(self, w, library, desc_lib):
        """Map coefficients from a library's own basis into jet space."""
        if library == "multigraph":
            return map_multigraph_to_jet(w, desc_lib, self.desc_jet)
        if library == "mod_multigraph":
            return map_mod_multigraph_to_jet(w, desc_lib, self.desc_jet)
        if library == "pdefind":
            return map_standard_to_jet(w, desc_lib, self.desc_jet)
        if library == "wsindy":
            return map_wsindy_to_jet(w, self.M_ws_to_jet)
        if library == "polyjet":
            w_jet = np.zeros(len(self.desc_jet))
            for j, desc in enumerate(desc_lib):
                if abs(w[j]) > 1e-14 and desc in self.desc_jet:
                    w_jet[self.desc_jet.index(desc)] += w[j]
            return w_jet
        raise ValueError(f"Unknown library: {library}")

    def evaluate(self, w, library, desc_lib):
        """fp, fn, exact_recovery, E_coeff, max_coeff_err, E_data.

        ``n_active`` is the count in the library's *own* basis, matching the
        sweep CSV; ``evaluate_unified`` counts in jet space, where a single
        equivariant term such as ``Δf`` expands into ``f_xx + f_yy``.
        """
        res = evaluate_unified(
            w, library, desc_lib, self.desc_jet, self.Theta, self.Ft,
            self.true_jet, M_ws_to_jet=self.M_ws_to_jet,
        )
        res["n_active_jet"] = res["n_active"]
        res["n_active"] = int(np.count_nonzero(np.abs(w) > 1e-10))
        return res

    def summary(self):
        n_true = int(np.sum(np.abs(self.true_jet) > 1e-10))
        print(f"Evaluation space: {len(self.desc_jet)} jet monomials, "
              f"{self.Theta.shape[0]:,} rows, {n_true} true active terms")


# ═══════════════════════════════════════════════════════════════════════════
# Reporting
# ═══════════════════════════════════════════════════════════════════════════

def show_recovery(w, desc, cfg, ev, label=None):
    """Print one recovery: the PDE, then its metrics. Returns the metrics."""
    res = ev.evaluate(w, cfg["library"], desc)
    head = label or LIBRARY_LABEL[cfg["library"]]
    print(f"── {head} ──")
    print(f"   {describe(cfg)}")
    print("   ", end="")
    print_pde(w, desc)
    verdict = "EXACT RECOVERY" if res["exact_recovery"] else "not exact"
    print(f"   {res['n_active']} active terms · fp={res['fp']} fn={res['fn']}"
          f" · {verdict}")
    print(f"   E_coeff = {res['E_coeff']:.4g}   "
          f"max|Δc| = {res['max_coeff_err']:.3g}")
    return res


def coefficient_table(w, desc, cfg, ev, tol=1e-10):
    """Recovered vs true coefficients, term by term, in the jet basis.

    The jet basis is used because it is the only space in which an
    equivariant and a non-equivariant recovery can be compared coefficient by
    coefficient: ``Δf`` and ``f_xx + f_yy`` are the same object there.
    """
    w_jet = ev.to_jet(w, cfg["library"], desc)
    rows = [
        (name, float(w_jet[i]), float(ev.true_jet[i]))
        for i, name in enumerate(ev.desc_jet)
        if abs(w_jet[i]) > tol or abs(ev.true_jet[i]) > tol
    ]
    rows.sort(key=lambda r: -abs(r[2]))

    width = max(len(r[0]) for r in rows) + 2
    print(f"{'term':<{width}}{'recovered':>12}{'true':>12}{'error':>12}")
    print("─" * (width + 36))
    for name, got, want in rows:
        flag = ""
        if abs(want) <= tol:
            flag = "  ← false positive"
        elif abs(got) <= tol:
            flag = "  ← MISSED"
        print(f"{name:<{width}}{got:>12.5f}{want:>12.5f}"
              f"{got - want:>12.2e}{flag}")


def compare_noise_levels(results, title="Clean vs noisy recovery"):
    """Table over (label, noise, metrics) triples, as collected in §6."""
    print(title)
    print(f"{'library':<22}{'noise':>7}{'n_act':>7}{'fp':>5}{'fn':>5}"
          f"{'exact':>7}{'E_coeff':>11}{'max|Δc|':>11}")
    print("─" * 75)
    for label, noise, res in results:
        print(f"{label:<22}{noise * 100:>6.0f}%{res['n_active']:>7}"
              f"{res['fp']:>5}{res['fn']:>5}"
              f"{'yes' if res['exact_recovery'] else 'no':>7}"
              f"{res['E_coeff']:>11.4g}{res['max_coeff_err']:>11.3g}")
