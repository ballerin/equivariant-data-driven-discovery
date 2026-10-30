"""
Generate and save simulation data for all benchmark equations.

Run once to pre-generate data for the unified ablation sweep.
Each (equation, noise_level) pair is saved as a .npz file in data/.

Usage
-----
    cd paper-experimental/code
    python scripts/generate_data.py
    python scripts/generate_data.py --equations ac sh
    python scripts/generate_data.py --noise 0.0 0.02
    python scripts/generate_data.py --data-dir data/
    python scripts/generate_data.py --profile large --anisotropic   # → data/large/
    python scripts/generate_data.py --profile gridtest   # → data/gridtest/{fine,coarse}/

Profiles
--------
``legacy`` (default, → ``data/``): the configurations every sweep in
``results/`` up to 2026-09 was run on.  64² grids (CH 128²) with 41–201
frames.

``large`` (→ ``data/large/``): sized for the Messenger & Bortz (2021)
test-function rule and the time-weak formulations, following the
dimensions of their 2-D examples (reaction-diffusion 256²×201, Navier-Stokes
324×149×201).  On the legacy grids the rule, at 2–5 % noise, finds only
k* ≈ 1–3 signal modes and asks for a test-function half-width m ≈ N/2.

* AC, SH, CH, PME: 256², at the legacy grid spacing Δx (so the domain side
  grows: AC 2π→8π, SH 8π→32π, CH 8π→16π, PME 6π→24π), same PDE parameters
  and time step, and ~200 frames after the transient skip of
  ``equations.EQUATION_META``.  The larger domain holds more interfaces /
  rolls / blobs rather than finer versions of the same ones.  PME starts from
  a jittered 4×4 lattice of the legacy Gaussian blob (same width, amplitudes
  in [0.5, 1]) instead of one blob, since a single blob stays a few-mode
  signal on any domain.
* KS (cell / wchaos / chaos): the regime is set by L (number of linearly
  unstable modes ∝ L², see ``ks_regimes_sweep.py``), so L and the legacy
  64² grid (already ~10 dealiased modes per unstable wavelength) are kept.
  T_end is kept too: observation noise is a fraction of the *global* std,
  which for KS is dominated by the linear drift of the mean, so a longer run
  would make the same "5 %" a larger fraction of the spatial fluctuations
  (T_end=150 roughly doubles it).  More frames come from a finer snapshot
  spacing instead (ETDRK4 dt 0.1→0.05, Δt_snap 0.5→0.25): 141 (cell) and 281
  (wchaos, chaos) frames after the 30 % skip.  At 2–5 % noise the rule still
  asks m ≈ N/2 in space for cell and wchaos, which requires
  ``sweep.py --space-boundary periodic``; experiments from 2026-09 use that
  boundary for every equation.

``gridtest`` (→ ``data/gridtest/fine/`` and ``data/gridtest/coarse/``): does a
finer grid help the test-function rule?  AC, SH, CH and PME on the *legacy*
domain L with the ``large`` time settings (same T_end and snapshot spacing),
simulated once on a 256² grid (Δx = legacy/4; CH legacy/2).  ``fine`` is that
data; ``coarse`` keeps every 4th grid point (CH every 2nd) of the *same* noisy
trajectory, i.e. the legacy 64² (CH 128²) spacing.  Noise is added once at the
fine resolution (σ = level · std of the clean fine field) and subsampled, so
the coarse data is literally a subset of the fine measurements and the pair
differs only in spatial resolution.  The simulators' random initial conditions
are drawn on the grid, so simulating each resolution separately would give
different systems.  PME uses dt = 2.5e-4 (explicit Euler at Δx/4; converged to
1e-5 against dt/2) and the legacy single-blob initial condition.

``paper`` (→ ``data/paper/``): the configuration the 2026-09 space-time weak
sweep is run on.  It is ``large`` with three changes, each answering a way in
which ``large`` still starved the M&B rule:

* **Noise is a fraction of each snapshot's own spatial std**, not of the
  global std (``--noise-scaling``, default ``per-snapshot`` for this profile).
  The two differ by ≤1.4× for AC/SH/CH/PME but by 7.6× (ks_cell) and
  21-24× (ks_wchaos, ks_chaos), because the KS global std is dominated by
  the drift of the spatial mean.  Under the old rule KS "5 %" was of order
  the spatial fluctuation itself, which is what drove k* to 1-2 there.
* **PME runs to T_end = 50** instead of 2 (dt unchanged; save_every 5 → 125
  keeps 201 frames, Δt_snap 0.01 → 0.25).  On every earlier profile PME's
  temporal spectrum had k_t = 1, so m_t sat at the (N_t−1)/2 cap and the
  space-time weak form had a *single* time window.  A longer trajectory is
  the only fix: the blob lattice spreads and merges over T ≈ 50, giving
  k_t = 8/6/4 and 131/111/71 windows at 0/2/5 % noise.  T_end = 20 is not
  enough (still capped at 5 %).
* **The KS regimes are put on a common Δx = 4π/64 = 0.196**: ks_cell 4π at
  64² (unchanged), ks_wchaos 8π at 128², ks_chaos 16π at 256².  The three
  KS datasets previously had Δx = 0.196 / 0.393 / 0.785, so the regime family
  confounded domain size with resolution.  At fixed Δx, m/N falls
  monotonically with L at every noise level (2 %: 0.48 / 0.34 / 0.14), which
  isolates the domain-size effect.

Note that KS cannot be given a larger domain at a fixed regime.  The regime
depends only on L/λ with λ = 2π/k*, k* = √(μ/2ν) for
f_t + νΔ²f + μΔf + ½|∇f|² = 0, and rescaling (ν, μ) to hold L/λ fixed while
L → αL maps the solution to the old one stretched by α (time by α²).  The
rule's k* is a mode index, so it is invariant under that stretch, and at
fixed Δx so is m/N.  ks_cell holds ≈ 1.4 cells and is therefore genuinely
single-mode at any noise level; that is a property of the cellular regime,
not of the pipeline.

M&B rule on the ``large`` data (τ̂ = 2; noise 0 / 2 / 5 %), m/N and m_t/Nt:
AC 0.05/0.07/0.08 and 0.09–0.10; SH 0.03/0.09/0.09 and 0.13–0.15;
CH 0.04/0.04/0.05 and 0.11–0.12; PME 0.05/0.14/0.15 and 0.5 (a monotone
decay: a single time window); KS chaos 0.12/0.27/0.31 and 0.07–0.08, KS cell
0.20/0.48/0.48 and 0.12–0.13.
"""

import argparse
import csv
import datetime
import os
import sys
import time

import numpy as np

# ── path setup ──────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.utils import (  # noqa: E402
    make_spatial_grid, build_noise_profile, normalize_profile,
    apply_anisotropic_noise,
)
from equivariant_pde_find.simulations import (  # noqa: E402
    allen_cahn,
    cahn_hilliard,
    kuramoto_sivashinsky,
    porous_medium,
    swift_hohenberg,
)

DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
os.makedirs(DATA_DIR, exist_ok=True)

DATA_DIRS = {"legacy": DATA_DIR,
             "large": os.path.join(DATA_DIR, "large"),
             "gridtest": os.path.join(DATA_DIR, "gridtest"),
             "paper": os.path.join(DATA_DIR, "paper")}

NOISE_LEVELS = [0.0, 0.02, 0.05]
NOISE_LABELS = {0.0: "clean", 0.02: "noise2pct", 0.05: "noise5pct"}

# ── equation configurations ──────────────────────────────────────────────────


def pme_blob_lattice(N, L, n_side=4, seed=42):
    """Jittered n_side × n_side lattice of Gaussian blobs on the periodic
    [0, L]², each of width (L / n_side) / 8 — the legacy single blob's width
    relative to its cell — with amplitudes in [0.5, 1], normalised to max 1."""
    rng = np.random.default_rng(seed)
    x = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    cell = L / n_side
    sigma = cell / 8
    f0 = np.zeros((N, N))
    for i in range(n_side):
        for j in range(n_side):
            cx = (i + 0.5 + rng.uniform(-0.2, 0.2)) * cell
            cy = (j + 0.5 + rng.uniform(-0.2, 0.2)) * cell
            amp = rng.uniform(0.5, 1.0)
            dxp = (X - cx + L / 2) % L - L / 2      # periodic distance
            dyp = (Y - cy + L / 2) % L - L / 2
            f0 += amp * np.exp(-(dxp**2 + dyp**2) / (2 * sigma**2))
    return f0 / f0.max()


LEGACY_EQUATIONS = [
    dict(
        name="ac",
        simulate_fn=allen_cahn,
        sim_kwargs=dict(N=64, L=2 * np.pi, T_end=8.0, eps=0.15, dt=0.01,
                        save_every=20, seed=42),
    ),
    dict(
        name="sh",
        simulate_fn=swift_hohenberg,
        sim_kwargs=dict(N=64, L=8 * np.pi, T_end=50.0, r=0.3, dt=0.1,
                        save_every=5, seed=42),
    ),
    dict(
        name="ch",
        simulate_fn=cahn_hilliard,
        sim_kwargs=dict(N=128, L=8 * np.pi, T_end=30.0, eps=0.2, dt=0.005,
                        save_every=30, seed=42),
    ),
    dict(
        name="pme",
        simulate_fn=porous_medium,
        sim_kwargs=dict(N=64, L=6 * np.pi, T_end=2.0, dt=0.002,
                        save_every=25, seed=42),
    ),
    dict(
        name="ks_cell",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=64, L=4 * np.pi, T_end=50.0, dt=0.1,
                        save_every=5, seed=42),
    ),
    dict(
        name="ks_wchaos",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=64, L=8 * np.pi, T_end=100.0, dt=0.1,
                        save_every=5, seed=42),
    ),
    dict(
        name="ks_chaos",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=64, L=16 * np.pi, T_end=100.0, dt=0.1,
                        save_every=5, seed=42),
    ),
]

LARGE_EQUATIONS = [
    dict(
        name="ac",
        simulate_fn=allen_cahn,
        sim_kwargs=dict(N=256, L=8 * np.pi, T_end=8.0, eps=0.15, dt=0.01,
                        save_every=4, seed=42),
    ),
    dict(
        name="sh",
        simulate_fn=swift_hohenberg,
        sim_kwargs=dict(N=256, L=32 * np.pi, T_end=50.0, r=0.3, dt=0.05,
                        save_every=5, seed=42),
    ),
    dict(
        name="ch",
        simulate_fn=cahn_hilliard,
        sim_kwargs=dict(N=256, L=16 * np.pi, T_end=30.0, eps=0.2, dt=0.005,
                        save_every=30, seed=42),
    ),
    dict(
        name="pme",
        simulate_fn=porous_medium,
        sim_kwargs=dict(N=256, L=24 * np.pi, T_end=2.0, dt=0.002,
                        save_every=5, seed=42),
        initial_condition=pme_blob_lattice,
    ),
    dict(
        name="ks_cell",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=64, L=4 * np.pi, T_end=50.0, dt=0.05,
                        save_every=5, seed=42),
    ),
    dict(
        name="ks_wchaos",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=64, L=8 * np.pi, T_end=100.0, dt=0.05,
                        save_every=5, seed=42),
    ),
    dict(
        name="ks_chaos",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=64, L=16 * np.pi, T_end=100.0, dt=0.05,
                        save_every=5, seed=42),
    ),
]

PAPER_EQUATIONS = [
    # AC, SH, CH: identical to ``large``.
    dict(
        name="ac",
        simulate_fn=allen_cahn,
        sim_kwargs=dict(N=256, L=8 * np.pi, T_end=8.0, eps=0.15, dt=0.01,
                        save_every=4, seed=42),
    ),
    dict(
        name="sh",
        simulate_fn=swift_hohenberg,
        sim_kwargs=dict(N=256, L=32 * np.pi, T_end=50.0, r=0.3, dt=0.05,
                        save_every=5, seed=42),
    ),
    dict(
        name="ch",
        simulate_fn=cahn_hilliard,
        sim_kwargs=dict(N=256, L=16 * np.pi, T_end=30.0, eps=0.2, dt=0.005,
                        save_every=30, seed=42),
    ),
    # PME: ``large`` in space, but T_end 2 -> 50 (save_every 5 -> 125 keeps
    # 201 frames, dt_snap 0.01 -> 0.25).  See the module docstring.
    dict(
        name="pme",
        simulate_fn=porous_medium,
        sim_kwargs=dict(N=256, L=24 * np.pi, T_end=50.0, dt=0.002,
                        save_every=125, seed=42),
        initial_condition=pme_blob_lattice,
    ),
    # KS ladder at the common spacing dx = 4*pi/64 = 0.196: L doubles, N
    # doubles with it.  ks_cell is the ``large`` configuration unchanged.
    dict(
        name="ks_cell",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=64, L=4 * np.pi, T_end=50.0, dt=0.05,
                        save_every=5, seed=42),
    ),
    dict(
        name="ks_wchaos",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=128, L=8 * np.pi, T_end=100.0, dt=0.05,
                        save_every=5, seed=42),
    ),
    dict(
        name="ks_chaos",
        simulate_fn=kuramoto_sivashinsky,
        sim_kwargs=dict(N=256, L=16 * np.pi, T_end=100.0, dt=0.05,
                        save_every=5, seed=42),
    ),
]

# Grid test: legacy L, ``large`` time settings, fine grid; ``coarsen`` is the
# subsampling factor to the legacy spacing (see the module docstring).
GRIDTEST_EQUATIONS = [
    dict(
        name="ac",
        simulate_fn=allen_cahn,
        sim_kwargs=dict(N=256, L=2 * np.pi, T_end=8.0, eps=0.15, dt=0.01,
                        save_every=4, seed=42),
        coarsen=4,
    ),
    dict(
        name="sh",
        simulate_fn=swift_hohenberg,
        sim_kwargs=dict(N=256, L=8 * np.pi, T_end=50.0, r=0.3, dt=0.05,
                        save_every=5, seed=42),
        coarsen=4,
    ),
    dict(
        name="ch",
        simulate_fn=cahn_hilliard,
        sim_kwargs=dict(N=256, L=8 * np.pi, T_end=30.0, eps=0.2, dt=0.005,
                        save_every=30, seed=42),
        coarsen=2,
    ),
    dict(
        name="pme",
        simulate_fn=porous_medium,
        sim_kwargs=dict(N=256, L=6 * np.pi, T_end=2.0, dt=2.5e-4,
                        save_every=40, seed=42),
        coarsen=4,
    ),
]

PROFILES = {"legacy": LEGACY_EQUATIONS, "large": LARGE_EQUATIONS,
            "gridtest": GRIDTEST_EQUATIONS, "paper": PAPER_EQUATIONS}

# How a noise level is turned into a standard deviation, per profile.
# "global": fraction of the std over every snapshot at once, applied
# inside the simulator (the pre-2026-09 behaviour, kept byte-identical).
# "per-snapshot": fraction of each snapshot's own spatial std, applied
# after simulation -- the definition generate_nd_data.py and
# sweep.py --post-hoc-noise already use.  Overridden by --noise-scaling.
NOISE_SCALING_DEFAULTS = {"legacy": "global", "large": "global",
                         "gridtest": "global", "paper": "per-snapshot"}

# ── seeded trajectories (train / validation / test) ─────────────────────────
# The unit of independence is the trajectory, so a split is a set of initial
# condition seeds, never a set of rows (see new-plan.md §4.1).  Pools are fixed
# once and never re-drawn; the ranges are disjoint so that a seed's split is
# readable from the seed itself.
SPLIT_SEED_BASE = {"train": 100, "val": 200, "test": 300}
DEFAULT_IC_SEED = 42          # the unseeded legacy/large/paper trajectories


def split_seeds(n_train, n_val, n_test):
    """Seed pools for a split, e.g. (3, 3, 5) -> train 101-103, val 201-203,
    test 301-305.  Returns {split: [ic_seed, ...]}."""
    counts = dict(train=n_train, val=n_val, test=n_test)
    return {s: [SPLIT_SEED_BASE[s] + i + 1 for i in range(counts[s])]
            for s in ("train", "val", "test")}


def noise_seed_for(ic_seed, noise):
    """Deterministic noise seed: ``ic_seed * 1000 + round(noise * 100)``.

    The same initial condition can therefore be re-used at several noise
    levels with a different realisation at each, and the mapping is
    reproducible from the filename alone.
    """
    return int(ic_seed) * 1000 + int(round(float(noise) * 100))


def seed_suffix(ic_seed):
    """``_s{ic_seed}`` for a seeded trajectory, empty for the default seed.

    Keeping the default unsuffixed is what lets every existing file, glob and
    results CSV stay valid.
    """
    return "" if int(ic_seed) == DEFAULT_IC_SEED else f"_s{int(ic_seed)}"

# ── anisotropic noise configurations ────────────────────────────────────────

ANISOTROPIC_NOISE_LEVELS = [0.02, 0.05]        # clean data is the same regardless
ANISOTROPIC_NOISE_LABELS = {0.02: "noise2pct", 0.05: "noise5pct"}
ANISOTROPIC_EQUATIONS = ["ac", "sh", "pme", "ch", "ks_cell", "ks_wchaos"]

ANISOTROPIC_PROFILES = [
    dict(profile_type="gaussian_decay", direction="radial", label="radial_decay"),
    dict(profile_type="gaussian_growth", direction="radial", label="radial_growth"),
    dict(profile_type="gaussian_decay", direction="axial", label="axial_decay"),
    dict(profile_type="gaussian_growth", direction="axial", label="axial_growth"),
]


def generate_one(eq_entry, obs_noise, data_dir, profile="legacy"):
    """Simulate one (equation, noise) pair and save as .npz."""
    name = eq_entry["name"]
    sim_fn = eq_entry["simulate_fn"]
    kwargs = eq_entry["sim_kwargs"].copy()
    kwargs["obs_noise"] = obs_noise
    if "initial_condition" in eq_entry:
        kwargs["f0"] = eq_entry["initial_condition"](
            kwargs["N"], kwargs["L"], seed=kwargs["seed"])

    noise_label = NOISE_LABELS[obs_noise]
    filename = f"{name}_{noise_label}.npz"
    path = os.path.join(data_dir, filename)

    if os.path.exists(path):
        print(f"  SKIP {filename} (already exists)", flush=True)
        return path, True

    t0 = time.time()
    snapshots, times = sim_fn(**kwargs)
    elapsed = time.time() - t0

    dt_snap = times[1] - times[0]
    N = snapshots.shape[1]
    L = kwargs["L"]

    np.savez_compressed(
        path,
        snapshots=snapshots.astype(np.float64),
        times=times.astype(np.float64),
        Lx=L,
        Ly=L,
        N=N,
        dt_snap=dt_snap,
        obs_noise=obs_noise,
        eq_name=name,
        profile=profile,
    )

    n_snaps = snapshots.shape[0]
    has_nan = bool(np.any(np.isnan(snapshots)))
    print(
        f"  SAVED {filename}  "
        f"({n_snaps} snaps, {N}×{N}, "
        f"range=[{snapshots.min():.3f}, {snapshots.max():.3f}], "
        f"std={snapshots.std():.3f}, "
        f"nan={has_nan}, "
        f"{elapsed:.1f}s)",
        flush=True,
    )
    return path, False


def _save_npz(path, snapshots, times, L, obs_noise, name, profile,
              ic_seed=None, noise_seed=None):
    extra = {}
    if ic_seed is not None:
        extra["ic_seed"] = int(ic_seed)
    if noise_seed is not None:
        extra["noise_seed"] = int(noise_seed)
    np.savez_compressed(
        path,
        snapshots=snapshots.astype(np.float64),
        times=times.astype(np.float64),
        Lx=L,
        Ly=L,
        N=snapshots.shape[1],
        dt_snap=times[1] - times[0],
        obs_noise=obs_noise,
        eq_name=name,
        profile=profile,
        **extra,
    )


def per_snapshot_sigma(clean):
    """Mean over snapshots of each snapshot's own spatial standard deviation.

    A single scale for the whole trajectory, so the noise is homoscedastic in
    time as it is under the ``global`` rule; what changes is *which* std sets
    it.  For KS the global std is dominated by the drift of the spatial mean
    (it is 7.6x the per-snapshot std for ks_cell and 21-24x for ks_wchaos and
    ks_chaos), so "5 % noise" there is of order the spatial fluctuation
    itself.  For AC, SH, CH and PME the two differ by at most 1.4x.
    """
    return float(np.mean([snap.std() for snap in clean]))


def generate_per_snapshot(eq_entry, noise_levels, data_dir, profile,
                          ic_seed=None):
    """Simulate once, then add each noise level scaled by the per-snapshot std.

    One clean trajectory serves every noise level (as in the gridtest path),
    so an equation is simulated once instead of once per level.  Returns a
    list of (noise, path, skipped).

    ``ic_seed`` overrides the configuration's seed, which drives the random
    initial condition in every simulator (and, for PME, the blob lattice).
    The noise realisation is drawn from :func:`noise_seed_for` instead, so
    the same trajectory gets an independent noise draw at each level and two
    trajectories never share one.  With ``ic_seed=None`` (or the default 42)
    the filenames and contents are exactly as before.
    """
    name = eq_entry["name"]
    kwargs = eq_entry["sim_kwargs"].copy()
    kwargs["obs_noise"] = 0.0
    if ic_seed is not None:
        kwargs["seed"] = int(ic_seed)
    if "initial_condition" in eq_entry:
        kwargs["f0"] = eq_entry["initial_condition"](
            kwargs["N"], kwargs["L"], seed=kwargs["seed"])
    suffix = seed_suffix(kwargs["seed"])

    paths = {noise: os.path.join(
                 data_dir, f"{name}{suffix}_{NOISE_LABELS[noise]}.npz")
             for noise in noise_levels}
    if all(os.path.exists(q) for q in paths.values()):
        for noise in noise_levels:
            print(f"  SKIP {os.path.basename(paths[noise])} (already exists)",
                  flush=True)
        return [(noise, paths[noise], True) for noise in noise_levels]

    t0 = time.time()
    clean, times = eq_entry["simulate_fn"](**kwargs)
    print(f"  simulated {name} {clean.shape} in {time.time() - t0:.1f}s",
          flush=True)
    L = kwargs["L"]
    sigma_unit = per_snapshot_sigma(clean)
    print(f"  sigma: per-snapshot={sigma_unit:.4f}  global={clean.std():.4f}  "
          f"ratio={clean.std() / sigma_unit:.2f}", flush=True)

    out = []
    for noise in noise_levels:
        path = paths[noise]
        if os.path.exists(path):
            print(f"  SKIP {os.path.basename(path)} (already exists)",
                  flush=True)
            out.append((noise, path, True))
            continue
        if noise > 0.0:
            rng = np.random.default_rng(noise_seed_for(kwargs["seed"], noise))
            data = clean + noise * sigma_unit * rng.standard_normal(clean.shape)
        else:
            data = clean
        _save_npz(path, data, times, L, noise, name, profile,
                  ic_seed=kwargs["seed"],
                  noise_seed=noise_seed_for(kwargs["seed"], noise))
        print(f"  SAVED {os.path.basename(path)}  "
              f"({data.shape[0]} snaps, {data.shape[1]}x{data.shape[2]}, "
              f"range=[{data.min():.3f}, {data.max():.3f}], "
              f"std={data.std():.3f}, "
              f"nan={bool(np.isnan(data).any())})", flush=True)
        out.append((noise, path, False))
    return out


def generate_gridtest_one(eq_entry, noise_levels, root_dir):
    """Grid-test profile: simulate once on the fine grid, add each noise level
    once at fine resolution, and save the fine data and its every-``coarsen``-th
    grid point subsample (see the module docstring).

    Returns a list of (subdir, noise, path) for the manifest.
    """
    name = eq_entry["name"]
    c = eq_entry["coarsen"]
    kwargs = eq_entry["sim_kwargs"].copy()
    kwargs["obs_noise"] = 0.0
    dirs = {sub: os.path.join(root_dir, sub) for sub in ("fine", "coarse")}
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)

    paths = {(sub, noise): os.path.join(d, f"{name}_{NOISE_LABELS[noise]}.npz")
             for sub, d in dirs.items() for noise in noise_levels}
    out = [(sub, noise, p) for (sub, noise), p in paths.items()]
    if all(os.path.exists(p) for p in paths.values()):
        print(f"  SKIP {name} (all files exist)", flush=True)
        return out

    t0 = time.time()
    clean, times = eq_entry["simulate_fn"](**kwargs)
    print(f"  simulated {name} {clean.shape} in {time.time() - t0:.1f}s",
          flush=True)
    L = kwargs["L"]
    sigma_unit = clean.std()
    for noise in noise_levels:
        if noise > 0.0:
            # One noise realisation per level, drawn at fine resolution.
            rng = np.random.default_rng([kwargs["seed"], int(round(noise * 1000))])
            data = clean + noise * sigma_unit * rng.standard_normal(clean.shape)
        else:
            data = clean
        for sub, arr in (("fine", data), ("coarse", data[:, ::c, ::c])):
            _save_npz(paths[(sub, noise)], arr, times, L, noise, name,
                      f"gridtest_{sub}")
            print(f"  SAVED {sub}/{os.path.basename(paths[(sub, noise)])}  "
                  f"({arr.shape[0]} snaps, {arr.shape[1]}², "
                  f"std={arr.std():.3f}, nan={bool(np.isnan(arr).any())})",
                  flush=True)
    return out


def generate_anisotropic_one(eq_entry, obs_noise, profile_spec, data_dir,
                             profile="legacy"):
    """Simulate clean, apply anisotropic noise, and save as .npz.

    Parameters
    ----------
    eq_entry : dict
        Equation config with keys ``name``, ``simulate_fn``, ``sim_kwargs``.
    obs_noise : float
        Overall noise fraction (e.g. 0.02).
    profile_spec : dict
        Keys ``profile_type``, ``direction``, ``label``.
    data_dir : str
        Output directory.
    """
    name = eq_entry["name"]

    noise_label = ANISOTROPIC_NOISE_LABELS[obs_noise]
    label = profile_spec["label"]
    filename = f"{name}_aniso_{label}_{noise_label}.npz"
    path = os.path.join(data_dir, filename)

    if os.path.exists(path):
        print(f"  SKIP {filename} (already exists)", flush=True)
        return path, True

    t0 = time.time()

    # Load existing clean data (re-use pre-generated clean .npz)
    clean_path = os.path.join(data_dir, f"{name}_clean.npz")
    if not os.path.exists(clean_path):
        print(f"  ERROR: clean data not found at {clean_path}", flush=True)
        return None, False
    clean_data = np.load(clean_path)
    clean = clean_data["snapshots"]
    times = clean_data["times"]
    L = float(clean_data["Lx"])
    N = int(clean_data["N"])

    # Build and normalise spatial profile
    X, Y = make_spatial_grid(N, L, L)
    P_raw = build_noise_profile(
        X, Y,
        profile_type=profile_spec["profile_type"],
        direction=profile_spec["direction"],
    )
    P = normalize_profile(P_raw)

    # Apply anisotropic noise
    rng = np.random.default_rng(42)  # fixed seed for reproducibility
    snapshots = apply_anisotropic_noise(clean, obs_noise, P, rng)

    elapsed = time.time() - t0
    dt_snap = times[1] - times[0]

    np.savez_compressed(
        path,
        snapshots=snapshots.astype(np.float64),
        times=times.astype(np.float64),
        Lx=L,
        Ly=L,
        N=N,
        dt_snap=dt_snap,
        obs_noise=obs_noise,
        eq_name=name,
        profile_type=profile_spec["profile_type"],
        profile_direction=profile_spec["direction"],
        profile_label=label,
        profile=profile,
    )

    n_snaps = snapshots.shape[0]
    has_nan = bool(np.any(np.isnan(snapshots)))
    print(
        f"  SAVED {filename}  "
        f"({n_snaps} snaps, {N}×{N}, "
        f"range=[{snapshots.min():.3f}, {snapshots.max():.3f}], "
        f"std={snapshots.std():.3f}, "
        f"nan={has_nan}, "
        f"{elapsed:.1f}s)",
        flush=True,
    )
    return path, False


def main():
    parser = argparse.ArgumentParser(
        description="Pre-generate simulation data for the unified sweep."
    )
    parser.add_argument(
        "--profile", choices=sorted(PROFILES), default="legacy",
        help="Dataset configuration set (see module docstring): legacy "
             "(64² grids, data/), large (256² AC/SH/CH/PME and longer KS, "
             "data/large/), gridtest (AC/SH/CH/PME at legacy L, fine 256² "
             "and its coarse subsample, data/gridtest/{fine,coarse}/) or "
             "paper (large, with PME run to T=50 and the KS ladder at a "
             "common Δx; per-snapshot noise; data/paper/).",
    )
    parser.add_argument(
        "--equations",
        nargs="+",
        default=None,
        help="Which equations to generate (default: all).",
    )
    parser.add_argument(
        "--noise",
        nargs="+",
        type=float,
        default=NOISE_LEVELS,
        help="Noise levels to generate (default: 0.0 0.02 0.05).",
    )
    parser.add_argument(
        "--anisotropic",
        action="store_true",
        default=False,
        help="Also generate anisotropic (space-dependent) noise datasets "
             f"for {', '.join(ANISOTROPIC_EQUATIONS)}.",
    )
    parser.add_argument(
        "--data-dir",
        default=None,
        help="Output directory for .npz files (default: data/ for the legacy "
             "profile, data/large/ for large).",
    )
    parser.add_argument(
        "--splits", nargs=3, type=int, metavar=("N_TRAIN", "N_VAL", "N_TEST"),
        default=None,
        help="Generate independent trajectories for a train / validation / "
             "test split, e.g. --splits 3 3 5.  Seeds come from disjoint "
             "pools (train 101.., val 201.., test 301..), filenames gain "
             "_s{ic_seed}, and the manifest records the split so a sweep "
             "cannot mistake one tier for another.  Mutually exclusive with "
             "--ic-seeds.",
    )
    parser.add_argument(
        "--ic-seeds", nargs="+", type=int, default=None,
        help="Explicit initial-condition seeds, for a one-off set of "
             "trajectories with no split labels.  The default seed 42 keeps "
             "its unsuffixed filename.",
    )
    parser.add_argument(
        "--noise-scaling", choices=["global", "per-snapshot"], default=None,
        help="What a noise level is a fraction of: the std over every "
             "snapshot at once (global, applied inside the simulator) or "
             "each snapshot's own spatial std (per-snapshot, applied after "
             "simulation).  Default: per-snapshot for the paper profile, "
             "global for every other profile, which reproduces the "
             "pre-2026-09 files byte for byte.",
    )
    args = parser.parse_args()

    if args.splits and args.ic_seeds:
        raise SystemExit("--splits and --ic-seeds are mutually exclusive.")
    if args.splits:
        pools = split_seeds(*args.splits)
        # (ic_seed, split) in a fixed order: train, then val, then test.
        trajectories = [(s, split) for split in ("train", "val", "test")
                        for s in pools[split]]
    elif args.ic_seeds:
        trajectories = [(s, "") for s in args.ic_seeds]
    else:
        trajectories = [(None, "")]       # the configuration's own seed

    noise_scaling = (args.noise_scaling
                     or NOISE_SCALING_DEFAULTS[args.profile])
    if (args.splits or args.ic_seeds) and noise_scaling != "per-snapshot":
        raise SystemExit(
            "--splits/--ic-seeds require --noise-scaling per-snapshot: the "
            "global path draws its noise inside the simulator from the same "
            "seed as the initial condition, so two trajectories could not "
            "have independent noise."
        )
    data_dir = args.data_dir or DATA_DIRS[args.profile]
    os.makedirs(data_dir, exist_ok=True)

    all_eqs = PROFILES[args.profile]
    eqs = all_eqs
    if args.equations is not None:
        eq_names = set(args.equations)
        eqs = [e for e in all_eqs if e["name"] in eq_names]
        if not eqs:
            print(f"No matching equations. Available: "
                  f"{[e['name'] for e in all_eqs]}")
            return

    noise_levels = args.noise

    print(f"Generating {len(eqs)} equations × {len(noise_levels)} noise levels "
          f"→ {len(eqs) * len(noise_levels)} files")
    print(f"Profile: {args.profile}   Output directory: {data_dir}")
    print(f"Noise scaling: {noise_scaling}")
    if args.splits or args.ic_seeds:
        shown = ", ".join(f"{s}{'(' + sp + ')' if sp else ''}"
                          for s, sp in trajectories)
        print(f"Trajectories ({len(trajectories)}): ic_seed {shown}")
    print()

    if args.profile == "gridtest":
        rows = {"fine": [], "coarse": []}
        for eq_entry in eqs:
            print(f"[{eq_entry['name']}]", flush=True)
            for sub, noise, path in generate_gridtest_one(
                    eq_entry, noise_levels, data_dir):
                rows[sub].append(dict(
                    equation=eq_entry["name"], noise=noise,
                    noise_label=NOISE_LABELS[noise],
                    file=os.path.basename(path), path=path,
                    profile_type="", profile_direction="", profile_label="",
                ))
        for sub, sub_rows in rows.items():
            with open(os.path.join(data_dir, sub, "manifest.csv"), "w",
                      newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(sub_rows[0]))
                writer.writeheader()
                writer.writerows(sub_rows)
        print(f"Done → {data_dir}/{{fine,coarse}}")
        return

    # ── manifest ─────────────────────────────────────────────────────────
    manifest_path = os.path.join(data_dir, "manifest.csv")
    manifest_rows = []

    total_start = time.time()
    n_generated, n_skipped = 0, 0

    for eq_entry in eqs:
        name = eq_entry["name"]
        for ic_seed, split in trajectories:
            seed_shown = (eq_entry["sim_kwargs"]["seed"] if ic_seed is None
                          else ic_seed)
            print(f"[{name}  ic_seed={seed_shown}"
                  + (f"  {split}]" if split else "]"), flush=True)
            if noise_scaling == "per-snapshot":
                emitted = generate_per_snapshot(
                    eq_entry, noise_levels, data_dir, args.profile,
                    ic_seed=ic_seed)
            else:
                emitted = [(noise,) + generate_one(eq_entry, noise, data_dir,
                                                   args.profile)
                           for noise in noise_levels]
            for noise, path, skipped in emitted:
                if skipped:
                    n_skipped += 1
                else:
                    n_generated += 1
                manifest_rows.append(dict(
                    equation=name,
                    noise=noise,
                    noise_label=NOISE_LABELS[noise],
                    file=os.path.basename(path),
                    path=path,
                    profile_type="",
                    profile_direction="",
                    profile_label="",
                    ic_seed=seed_shown,
                    noise_seed=(noise_seed_for(seed_shown, noise)
                                if noise > 0 else ""),
                    split=split,
                ))
            print()

    # ── anisotropic datasets ──────────────────────────────────────────────
    if args.anisotropic and noise_scaling != "global":
        raise SystemExit(
            "--anisotropic is only implemented for --noise-scaling global: "
            "apply_anisotropic_noise scales by the global std, so mixing it "
            "with per-snapshot isotropic data would put two different "
            "definitions of the same noise level in one manifest."
        )
    if args.anisotropic:
        aniso_eqs = [e for e in all_eqs
                     if e["name"] in ANISOTROPIC_EQUATIONS]
        # If --equations was specified, only generate for those that also
        # appear in ANISOTROPIC_EQUATIONS.
        if args.equations is not None:
            req = set(args.equations)
            aniso_eqs = [e for e in aniso_eqs if e["name"] in req]

        if aniso_eqs:
            n_profiles = len(ANISOTROPIC_PROFILES)
            n_aniso = (len(aniso_eqs)
                       * len(ANISOTROPIC_NOISE_LEVELS)
                       * n_profiles)
            print(f"Anisotropic: {len(aniso_eqs)} equations "
                  f"× {len(ANISOTROPIC_NOISE_LEVELS)} noise levels "
                  f"× {n_profiles} profiles → {n_aniso} files")
            print()

            for eq_entry in aniso_eqs:
                name = eq_entry["name"]
                print(f"[{name} anisotropic]", flush=True)
                for noise in ANISOTROPIC_NOISE_LEVELS:
                    for prof in ANISOTROPIC_PROFILES:
                        path, skipped = generate_anisotropic_one(
                            eq_entry, noise, prof, data_dir, args.profile,
                        )
                        if skipped:
                            n_skipped += 1
                        else:
                            n_generated += 1
                        manifest_rows.append(dict(
                            equation=name,
                            noise=noise,
                            noise_label=ANISOTROPIC_NOISE_LABELS[noise],
                            file=os.path.basename(path),
                            path=path,
                            profile_type=prof["profile_type"],
                            profile_direction=prof["direction"],
                            profile_label=prof["label"],
                            ic_seed=eq_entry["sim_kwargs"]["seed"],
                            noise_seed="",
                            split="",
                        ))
                print()

    # Write manifest
    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "equation", "noise", "noise_label", "file", "path",
            "profile_type", "profile_direction", "profile_label",
            "ic_seed", "noise_seed", "split",
        ])
        writer.writeheader()
        writer.writerows(manifest_rows)

    total_elapsed = time.time() - total_start
    print(f"Done. {n_generated} generated, {n_skipped} skipped "
          f"in {total_elapsed:.0f}s.")
    print(f"Manifest → {manifest_path}")


if __name__ == "__main__":
    main()
