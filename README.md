# Equivariant PDE discovery

This repository is intended as reference for the paper "Equivariant PDE discovery" by Ballerin & Grong [1] for discovery of partial differential equations in a 2D flat domain with periodic boundary conditions. It contains the core Python code used to:
- simulate benchmark PDEs, optionally with noise;
- construct an O(2)-equivariant library of nonlinear partial differential operators;
- build the standard coordinate libraries (non-equivariant) for comparison;
- apply sparse regression to recover PDEs from data;


The implementation provides several differentiation backends, library building strategies and regression algorithms, so experiments can trade accuracy for noise-robustness and an ablation study can be performed. Details on the implemented methods can be found below.

## Motivation
The code exploits the complete classification of the multigraph-operators from Ballerin & Grong [2]: a minimal, complete, $O(2)$-equivariant basis of polynomial nonlinear partial differential operators. The general aim of the research is to show that whenever the underlying dynamics has the additional assumption of $O(2)$-equivariance, introducing bias in the form of a complete and non-redundant library of multigraph operators yields better recovery when deducing the generating PDE.

The code is organised to make it straightforward to (re)generate datasets, assemble library matrices, and run hyperparameter sweeps that compare the equivariant basis from the paper against standard non-equivariant libraries such as those used in SINDy/PDE-FIND  and WSINDy.

## What is implemented

The repository implements the following strategies, of which the first one is our novel contribution, and the other two are benchmarks from prior published work:
- $O(2)$-equivariant library based on multigraph operators;
- PDE-FIND library (based on SINDy algorithm) [3];
- WSINDy library and algorithm [4].

Two differentiation backends are available:
- **Spectral** (FFT-based): exact for band-limited data, but highly susceptible to noise.
- **Polynomial stencils** (local least-squares): fit a 2-D polynomial on a $(2r+1)\times(2r+1)$ stencil at each grid point. More noise-robust than spectral at the cost of a small approximation bias.

Several regression algorithms:
- STRidge, originally from PDE-FIND [3];
- MSTLS, as described in [4];
- RFE (recursive feature elimination);
- Brute-force exhaustive search;
- Oracle (ground-truth support).

## Pipeline

The recommended workflow has four stages: generate data, run all sweeps,
extract the best results.

### Setup

**Option A: Conda environment (recommended)**

A conda environment file ``environment.yml`` is provided for reproducible setup:

```bash
conda env create -f environment.yml
conda activate PDEDisc
pip install -e .
```

This creates an environment named ``PDEDisc`` with Python 3.11 and all required dependencies (numpy, scipy, pytest, tqdm, sympy).

**Option B: pip only**

Install in editable mode:

    pip install -e .

This makes the ``equivariant_pde_find`` package importable from anywhere.
Scripts in the ``scripts/`` directory automatically add ``src/`` to ``sys.path`` so they also work without a pip install.

Requires Python >= 3.10. `pyproject.toml` declares loose lower bounds; `requirements-lock.txt` pins the exact versions this codebase has been tested against, if you want an exact reproduction of the development environment (`pip install -r requirements-lock.txt`).

#### 1. Generate simulation data

Generate all datasets with the default parameters (all equations, all noise levels):

```bash
python scripts/generate_data.py
```

To generate only specific equations or noise levels:

```bash
# Generate only Allen-Cahn and Swift-Hohenberg equations
python scripts/generate_data.py --equations ac sh

# Generate only noise levels 0.0 and 0.02
python scripts/generate_data.py --noise 0.0 0.02
```

Output: `.npz` files in `data/` plus a `data/manifest.csv` index.

##### Dataset profiles

`--profile` selects a configuration set (details and literature
references in the `generate_data.py` docstring):

| Profile | Output | Grids and frames |
|---|---|---|
| `legacy` (default) | `data/` | 64² (CH 128²), 41–201 frames: every sweep in `results/` up to 2026-09 |
| `large` | `data/large/` | AC, SH, CH, PME at 256² with the legacy grid spacing (so a 2–4× larger domain) and ~200 frames; KS at the legacy 64² and L, with twice as many frames |
| `gridtest` | `data/gridtest/{fine,coarse}/` | AC, SH, CH, PME at the legacy L, simulated at 256² and subsampled to the legacy spacing: resolution varied with nothing else |
| `paper` | `data/paper/` | `large`, with PME run to T=50 and the KS regimes put on a common Δx; per-snapshot noise. The 2026-09 space-time weak sweep |

```bash
python scripts/generate_data.py --profile large --anisotropic
python scripts/sweep.py --data-dir data/large ...

# The 2026-09 configuration
python scripts/generate_data.py --profile paper
python scripts/probe_test_fn_rule.py --data-dir data/paper   # rule pre-flight
```

Use `large` for the time-weak formulations and the Messenger & Bortz test-function
rule (`sweep.py --test-fn-rule mb`). On the legacy grids the rule finds only
1–3 signal-dominated Fourier modes at 2–5 % noise and asks for a test function as
wide as half the domain. The larger domains hold more structures (interfaces,
rolls, blobs) at the same resolution, and the rule then picks half-widths of
0.04–0.15 N in space and 0.07–0.15 Nt in time. The exception is PME, whose
monotone spreading always gets one wide time window.

In `large`, KS keeps its legacy L and 64² grid, because L sets the regime
(cellular, weakly chaotic, chaotic), and its legacy T_end: the noise there is a
fraction of the *global* std, which for KS is dominated by the drifting mean, so
a longer run would turn the same "5 %" into a harsher noise level. At 2–5 %
noise KS-cell and KS-wchaos then still need `--space-boundary periodic`, the
boundary used for every equation (see *Spatial boundary* below). The `paper`
profile addresses both points directly.

The `paper` profile differs from `large` in three ways, each closing a way in
which `large` still starved the test-function rule:

- **PME runs to `T_end = 50`** instead of 2 (`save_every` 5 → 125 keeps 201
  frames). On every earlier profile PME's temporal spectrum had k\* = 1, so
  m_t sat at its (N_t−1)/2 cap and the space-time weak form had a *single*
  time window. At T=50 the rule gives m_t = 35/45/65 and 131/111/71 windows at
  0/2/5 % noise. T=20 is not enough — still capped at 5 %.
- **The KS regimes share Δx = 4π/64 = 0.196**: ks_cell 4π at 64², ks_wchaos 8π
  at 128², ks_chaos 16π at 256². Previously the three had Δx = 0.196 / 0.393 /
  0.785, so the regime family confounded domain size with resolution. At fixed
  Δx, m/N falls monotonically with L at every noise level (at 2 %: 0.48 / 0.34
  / 0.14).
- **Noise scales with each snapshot's own std** (see below).

KS cannot be given a larger domain at a fixed regime. The regime depends only
on L/λ with λ = 2π/k\*, k\* = √(μ/2ν) for f_t + ν∆²f + µ∆f + ½|∇f|² = 0, and
rescaling (ν, µ) to hold L/λ fixed while L → αL maps the solution to the old
one stretched by α (time by α²). The rule's k\* is a mode index, so it is
invariant under that stretch, and at fixed Δx so is m/N. ks_cell holds ≈ 1.4
cells and is therefore genuinely single-mode at any noise level — a property
of the cellular regime, not of the pipeline.

##### How the noise level is scaled

A noise level is a fraction of the field's standard deviation — but *which* standard
deviation differs by code path, so the same number does not mean the same thing everywhere:

| Path | Noise standard deviation |
|---|---|
| `generate_data.py --noise-scaling global` (`legacy`, `large`, `gridtest`), via the simulators' `obs_noise` | fraction of the **global** std, taken over every snapshot at once |
| `generate_data.py --noise-scaling per-snapshot` (the `paper` default) | fraction of the mean **per-snapshot** spatial std, added after simulation |
| `generate_nd_data.py` (all `*_nd*` datasets) | fraction of **each snapshot's own** std |
| `sweep.py --post-hoc-noise` | fraction of **each snapshot's own** std |

The two definitions differ by at most 1.4× for AC, SH, CH and PME, but by 7.6×
(ks_cell) and 21–24× (ks_wchaos, ks_chaos), because the KS global std is
dominated by the drift of the spatial mean. Under the global rule a KS "5 %"
was of order the spatial fluctuation itself, which is what drove the rule's
k\* to 1–2 there; rescaling lifts ks_chaos from k\* = 5–6 to 8–9 and takes
ks_wchaos off its cap at 2 %.

`--noise-scaling` defaults to `global` for every profile except `paper`, so
the pre-2026-09 files regenerate byte for byte. `--anisotropic` is implemented
for `global` only, since `apply_anisotropic_noise` scales by the global std;
combining them is refused rather than silently mixing two definitions of the
same noise level in one manifest.

#### 2. Run the main ablation sweep

Full sweep (all libraries, formulations, backends, and regression algorithms).

```bash
python scripts/sweep.py
```

Targeted sweeps for faster iteration:

```bash
# Multigraph library vs PDE-FIND library comparison (pointwise only, STRidge only)
python scripts/sweep.py --methods pointwise --algorithms stridge --libraries multigraph pdefind

# Noise-robustness comparison (space-weak vs pointwise, multigraphs only)
python scripts/sweep.py \
    --equations sh ks_chaos \
    --noise 0.02 0.05 \
    --methods pointwise space_weak_direct space_weak_ibp \
    --libraries multigraph

# With post-hoc noise injection and spectral low-pass prefiltering
python scripts/sweep.py --post-hoc-noise 0.02 0.05 --spectral-cutoff 3.0

# MSTLS-only with custom alpha sweep
python scripts/sweep.py --algorithms mstls --mstls-alpha 0.0 0.25 0.5

# Polynomial jet library only
python scripts/sweep.py --libraries polyjet

# WSINDy library only
python scripts/sweep.py --libraries wsindy

# Current standard (2026-09): space-time weak form with IBP, test functions by
# the Messenger & Bortz rule, periodic spatial boundary, 256² data
python scripts/sweep.py --data-dir data/large \
    --methods spacetime_weak_ibp --test-fn-rule mb --tau-hat 2 \
    --space-boundary periodic --centering off on --eval-mode qr

# The 2026-09 main sweep, one (equation, noise) at a time.  One mode, every
# library that can represent the equation, both centering arms.  The driver
# results/paper_20260918/run.sh runs this over all 21 datasets and skips any
# chunk whose CSV already exists, so it is resumable.
python scripts/sweep.py --data-dir data/paper --equations ac --noise 0.02 \
    --methods spacetime_weak_ibp --backends polynomial \
    --libraries multigraph mod_multigraph pdefind polyjet wsindy \
    --representable-only \
    --test-fn-rule mb --tau-hat 2 --space-boundary periodic \
    --centering off on --centering-tol 0.05 \
    --algorithms stridge rfe mstls \
    --rfe-lam 0 1e-6 1e-4 1e-2 --mstls-alpha 0.25 0.5 0.75 \
    --eval-mode qr --out results/paper_20260918/sweep_ac_0.02.csv

# Pointwise vs space-weak vs time-weak vs space-time weak, valid vs periodic boundary
python scripts/sweep.py --data-dir data/large \
    --methods pointwise space_weak_ibp time_weak spacetime_weak_ibp \
    --test-fn-rule mb --space-boundary valid periodic --centering off on
```

Methods (`--methods`). The name says where the test function lives. For methods
that are weak in space, the suffix says how spatial derivatives are handled:
`direct` computes them on the data with the derivative backend and then integrates
the column; `ibp` moves them onto the test function. Whether the *time* derivative
is integrated by parts is a separate choice, `--time-ibp`, recorded in the
`time_ibp` CSV column.

| Method | Test function in space | Test function in time | Spatial derivatives |
|---|---|---|---|
| `pointwise` | – | – | backend |
| `space_weak_direct` | ✓ | – | backend, then integrated |
| `space_weak_ibp` | ✓ | – | on the test function |
| `time_weak` | – | ✓ | backend |
| `spacetime_weak_direct` | ✓ | ✓ | backend, then integrated |
| `spacetime_weak_ibp` | ✓ | ✓ | on the test function |

Before 2026-09, `space_weak_direct` and `space_weak_ibp` were called `weak_direct`
and `weak_ibp`, and the `space_m` column was `weak_m` (`use_ibp` in `sweep_nd.py`).
The old method names and `--weak-m` are still accepted on the command line. All
summary and table builders translate old CSVs through
`summary_metadata.normalize_sweep_row`. In the library API, `formulation="weak"` and
`use_ibp=` are likewise aliases of `formulation="space_weak"` and `space_ibp=`.

Spatial boundary (`--space-boundary`). With `valid` (the default, as in Messenger &
Bortz), a test function is centred only where its whole support fits inside the
grid, which loses 2m points per axis. With `periodic`, the support wraps around the
domain and every grid point is a centre. On the periodic domains simulated here
this is exact: the integration-by-parts boundary terms cancel on a torus. It is
the only way to keep query points when the rule asks for m ≈ N/2 (KS-cell,
KS-wchaos). The switch lives in `library/weak_form.py`
(`use_weak_boundary("periodic")`) and affects every spatial weak convolution in
every library; the time axis always uses `valid`. `sweep.py` is currently the only
sweep exposing it.

**Experiments from 2026-09 on use `periodic` for every equation.** The reasons:

- All simulated domains are periodic, so the wrapped support is exact.
- The derivative backends already assume periodicity: the spectral backend by
  construction, and the polynomial backend fits its stencils on a periodic
  neighbourhood. IBP moves only the outermost derivative onto the test function,
  so `valid` would not make the pipeline periodicity-free.
- KS-cell and KS-wchaos need it, because the rule asks for m ≈ N/2 there. So do
  small grids in general: on the legacy 64² AC data a `valid` space-time system
  has a single row.
- One boundary then covers every equation and grid size.

The consequence is that the method, as run, assumes periodic data, and results
must say so. `valid` remains the CLI default only so that earlier sweeps reproduce
unchanged.

Key CLI options for `sweep.py`:

| Flag | Default | Description |
|---|---|---|
| `--methods` | `pointwise space_weak_direct space_weak_ibp` | Formulations (table above); also `time_weak`, `spacetime_weak_direct`, `spacetime_weak_ibp` (see `library/spacetime.py`) |
| `--test-fn-rule` | `legacy` | `legacy`: `--space-m` with fixed p, `--time-m`; `mb`: (m, p, m_t, p_t) per dataset by the Messenger & Bortz rule |
| `--space-boundary` | `valid` | `valid` or `periodic` spatial test-function support (pass both to compare). Current experiments use `periodic`; the default is kept for reproducibility |
| `--time-m` | `4 8 12` | Temporal half-supports (legacy rule) |
| `--tau`, `--tau-hat` | `1e-10`, `2` | M&B decay tolerance and Fourier-tail parameter |
| `--time-ibp` | `on` | `on`: b = φ_t′∗f; `off`: φ_t∗f_t with FD f_t (pass both to compare) |
| `--centering` | `off` | Fit through `regression.centering` (pass `off on` to compare) |
| `--eval-mode` | `dense` | `dense`: score against the full clean evaluation matrix. `qr`: use the R factor of its QR decomposition (`evaluation.FactoredEval`); E_coeff agrees to round-off, all other columns are identical, and scoring cost no longer grows with the grid (also in `sweep_aniso.py`, `sweep_ac_ndim.py`, `sweep_multigraph_comparison.py`) |
| `--space-stride`, `--time-stride` | `4`, `1` | Query-point strides for the methods that are weak in time |
| `--algorithms` | `stridge mstls rfe` | Regression algorithms (`stridge`, `mstls`, `rfe`, `brute`, `oracle`) |
| `--libraries` | `multigraph pdefind polyjet wsindy` | Library types (`multigraph`, `mod_multigraph`, `pdefind`, `polyjet`, `wsindy`) |
| `--representable-only` | off | Per dataset, keep only the libraries that can represent the equation (`equations.REPRESENTABLE_LIBRARIES`): PDE-FIND is dropped for CH, PME and KS, WSINDy for KS |
| `--backends` | `spectral polynomial` | Differentiation backends |
| `--equations` | all | Filter by equation name prefix |
| `--noise` | all | Filter by noise level (e.g. `0.0 0.02 0.05`) |
| `--space-m` | `4 8 12` | Spatial test-function half-supports (legacy rule; `--weak-m` is the old name) |
| `--post-hoc-noise` | none | Add noise to clean data as fraction of std |
| `--spectral-cutoff` | none | Low-pass filter cutoff wavenumber |
| `--name` | `""` | Tag appended to output filename |

Output: `results/sweep_<YYYYMMDD>_<HHMMSS>_<name>.csv`

Two CSV columns isolate a spurious constant term. `regression.centering` writes
the recovered intercept back into the constant column's slot, so an intercept
that survives `--centering-tol` is counted in `fp` exactly like any other false
positive — yet it is the one term centering exists to remove. `const_active`
says whether that slot is active, and `intercept_rel` gives |w₀| / mean|f_t|,
the quantity `--centering-tol` thresholds. Both are blank for a library with no
constant column. (On `ks_cell`, all 20 STRidge settings activate a constant with
centering off and none do with it on.)

#### 3. Run auxiliary sweeps

**KS regime sweep** (domain-size exploration across cellular → chaotic):

```bash
python scripts/ks_regimes_sweep.py
python scripts/ks_regimes_sweep.py --algorithm sr3 --name sr3_only
```

**Polynomial stencil sweep** (degree & radius trade-off on clean Allen-Cahn):

```bash
python scripts/sweep_polynomial_stencil.py
python scripts/sweep_polynomial_stencil.py --quick          # faster, lower resolution
python scripts/sweep_polynomial_stencil.py --out my_stencil_sweep.csv
```

#### 4. Train / validation / test with seeded trajectories

`sweep.py` fits every setting and scores each against the truth, so the setting
is chosen knowing the answer — an oracle. `train_val_test.py` runs the honest
protocol instead: choose a setting from held-out trajectories without reading
the truth, then see how often it works on data it has never seen.

```bash
# independent trajectories, one per split
python scripts/generate_data.py --profile paper --splits 3 3 5 \
    --data-dir data/paper_split

python scripts/train_val_test.py \
    --data-dir data/paper_split --out-dir results/tvt_<date>
```

| tier | what happens |
|---|---|
| train | for each setting, fit on each train trajectory → `w` |
| validation | score every `w` on held-out trajectories → `E_val` (no ground truth) |
| select | pick a setting from `E_val` alone |
| test | refit at the frozen setting, once per test trajectory, score against truth |

**The unit of independence is the trajectory.** Rows are never split:
neighbouring frames are near-duplicates, a weak row centred at *t* integrates
*t*−m_t … *t*+m_t, and spectral derivatives make every row in a frame depend on
the whole frame. Trajectories differ by initial-condition seed and carry
independent noise draws (`noise_seed = ic_seed·1000 + noise_pct`).

**What transfers is the setting, not the coefficients.** Each test trajectory is
refit from scratch, which is what a practitioner does with a fresh dataset and
gives one independent recovery outcome per test trajectory — the error bars a
single-seed sweep cannot produce.

**The M&B rule is not a hyperparameter.** It is a deterministic function of the
data, so each trajectory gets its own (m, p, m_t, p_t) and nothing is frozen.
Coefficients are in PDE units, so a `w` fitted where the rule gave m_t = 20
scores correctly where it gave 17 — checked directly: the true coefficient
vector gives residuals 0.0025 / 0.0024 / 0.0031 across three `ks_cell`
trajectories.

Two outputs per cell, in `--out-dir`: `selection_<eq>_<noise>.csv` (every
setting, with `E_val` for each train × validation pair) and
`test_<eq>_<noise>.csv` (one row per test trajectory). The script applies no
rule to the residuals beyond choosing what to carry into the test tier
(`--rule delta|one-se`), so the summariser can apply both to the same numbers.
`--rule one-se` needs at least 2 train trajectories for a standard error and
falls back to the δ rule with a warning below that.

Centering stays a **compared factor**, not a hyperparameter (§7.2 of
`new-plan.md`): each arm gets its own selection, because selecting over it
would be a second oracle.

#### 5. Summarise the space-time weak sweep (centering as a compared factor)

```bash
python scripts/probe_test_fn_rule.py --data-dir data/paper

python scripts/summarize_paper.py results/paper_20260918/sweep_*.csv \
    --delta 0.5 \
    --out results/paper_20260918/summary.txt \
    --csv results/paper_20260918/summary.csv
```

`probe_test_fn_rule.py` is the cheap pre-flight: it runs no regression and
prints, per (equation, noise), the rule's (m, p, m_t, p_t), the critical
wavenumbers behind them, whether either support sits at its cap, and the row
count the sweep will have to fit.

`summarize_paper.py` reports, per equation × noise × library × solver and
separately for each centering arm: the truth-free "sparsest within (1+δ) of the
best residual" selection, `best_of_sweep`'s oracle selection for continuity, the
fraction of settings that recover, E_coeff over recovering rows only, and the
spurious-constant columns. It also pairs centering on against off at identical
settings and reports fixed / broken counts with an exact McNemar test.

It replaces the former `summarize_best.py` and `build_summary.py`, removed in
the 2026-09-19 cleanup: neither knew the `centering` column, so both would have
grouped the two arms together and taken the best of them. Centering is a
compared factor, not a hyperparameter — selecting over it is a second oracle,
subtler than the λ one and just as wrong.

Two caveats belong in any caption built from its output. The residual driving
the δ rule is `E_data`, measured **in sample**, because the protocol's held-out
block is not implemented yet; and there is one realisation per dataset, so
there are no error bars over seeds.

### Results (`results/`)

Sweep outputs are CSV files named `sweep_<YYYYMMDD>_<HHMMSS>_<config>.csv`, or,
for a chunked run, one file per `(equation, noise)` under a dated directory.
What is kept here is deliberately narrow: the provenance for the published
tables, plus current work.

| Path | What it backs |
|---|---|
| `sweep_paper_v2_*_20260804_223951.*` | the 2D results and ablation tables |
| `sweep_aniso_default_*_20260805_215410.*` | the space-varying-noise table |
| `sweep_nd_ndim_2d4d_v2_20260806_145657.*` | the higher-dimension table |
| `cost_scaling_20260809_*.csv` | the computational-cost table |
| `collinearity_diagnostics_spectral.*`, `collinearity_pct_tables.tex` | the conditioning appendix |
| `ablation_tables.tex` | the ablation appendix |
| `figures/` | the equation-motion figures (built from `notebooks/`) |
| `gridtest_20260917/` | the grid and domain study — see `grid-and-domain-study.md` |
| `paper_20260918/` | the space-time weak sweep on `data/paper` |

Superseded development runs were removed in the 2026-09-19 cleanup; the
`grid-and-domain-study.md` note records the conclusions of the studies whose
data no longer exists, so they need not be repeated.

## N-dimensional extension (experimental)

The paper's results are all 2D. In parallel, the whole discovery pipeline — simulators, derivative backends (spectral and polynomial), all 5 libraries (multigraph, mod_multigraph, pdefind, polyjet, wsindy) in both pointwise and weak/IBP formulations, and cross-library jet-space scoring — generalizes to arbitrary spatial dimension via an explicit `ndim` parameter. Domains are isotropic beyond 3D (a single length reused for every axis, matching every N-D simulator); axis names beyond `z` use `w, v, u, t, s` (up to ndim=8).

- **Equations**: Allen-Cahn and Cahn-Hilliard have `ndim`-generic simulators (`allen_cahn_nd`, `cahn_hilliard_nd`, 1 ≤ ndim ≤ 8) and ground truth. Swift-Hohenberg and Porous Medium also have `ndim`-generic simulators (`swift_hohenberg_nd`, `porous_medium_nd`), but are **not yet wired into `generate_nd_data.py`/`sweep_nd.py`** below. Kuramoto-Sivashinsky explicitly does **not** generalize: `kuramoto_sivashinsky(ndim=...)` raises `ValueError` for any value other than 2, because its ETDRK4/IMEX integrators (`simulations/_etdrk4.py`) are hardcoded around 2-D FFTs and have not been validated at higher dimension.
- **Libraries**: all 5 library types (`multigraph`, `mod_multigraph`, `pdefind`, `polyjet`, `wsindy`) are dimension-generic for any ndim, pointwise and weak (with IBP). The equivariant (multigraph/mod_multigraph) operator basis is generalized for `max_p<=2` (the 20-operator basis this repo's results actually use) — `max_p=3` raises `NotImplementedError` outside ndim in `{2, 3}`, since 3 of its 11 operators already required a full 3rd-order derivative tensor contraction that was never implemented beyond 3D in the first place. Verified end-to-end at ndim=4: simulate → build weak-form multigraph library with IBP → RFE regression → exact recovery (see `tests/test_library_nd.py`).
- **Polynomial derivative backend**: also generalized to N-D (local least-squares fit over an ndim stencil), alongside the spectral backend — but cost grows as `(2r+1)^ndim`, so expect it to become impractical well before the spectral backend does as ndim grows.

Generate n-D data (Allen-Cahn and Cahn-Hilliard only):

```bash
python scripts/generate_nd_data.py
python scripts/generate_nd_data.py --equations ac ch --ndim 2 3
python scripts/generate_nd_data.py --noise 0.0 0.02
```

Output: `data/<equation>_nd<ndim>_<label>.npz`.

Run the n-D sweep (multigraph / mod_multigraph / pdefind, STRidge / RFE, weak-form IBP):

```bash
python scripts/sweep_nd.py
python scripts/sweep_nd.py --equations ac ch --ndim 2 3
python scripts/sweep_nd.py --libraries multigraph pdefind --algorithms stride
```

Output: `results/sweep_nd_<YYYYMMDD>_<HHMMSS>_<name>.csv`.

## Libraries

We compare classes of libraries with maximum order `max_p=2` and maximum polynomial degree `max_n=3` for the 2D flat settings.

Notice that in general the terms in the library are not linearly independent.

For the equivariant library, whether the multigraph operators are linearly independent **depends on the dimension of the underlying space** — this is one of the main points of [2], whose abstract describes the classification as revealing "non-trivial linear dependence relations between nonlinear differential operators relative to the dimension of the manifold". Concretely:

- **Up to total order 4** (`max_p=2`, the setting used throughout this repository), the count is the same in every dimension $d \ge 2$ and all **20 operators are linearly independent**, in 2D as well as in higher dimensions. This is the complete, non-redundant basis of [1, Table 2]. `tests/test_multigraph_operators.py` pins the set and its ordering; `audit/repro/A3_basis_rank.py` measures the rank as 20 for $d \in \{2, 3, 4\}$, both from the symbolic jet decomposition and numerically from operator values.
- **At total order 6** (`max_p=3`) the dimension starts to matter. Of the 11 multigraph operators with 3 edges, exactly **two become linearly dependent when $d = 2$** — they are identically zero there, and non-zero for every $d \ge 3$. This is the theorem in [2] on the operators $P^{\beta_0}, Q^{\beta_0}$, which generate the kernel of the multigraph-to-operator map restricted to 3-edge multigraphs precisely when $d = 2$:

  $$P f = \|\nabla f\|^2 (\Delta f)^2 - \|\nabla f\|^2\|\nabla^2 f\|^2 + 2\|\nabla^2 f(\nabla f,\cdot)\|^2 - 2\Delta f\cdot\nabla^2 f(\nabla f,\nabla f)$$
  $$Q f = (\Delta f)^3 - 3\Delta f\,\|\nabla^2 f\|^2 + 2\,\mathrm{tr}\,(\nabla^2 f)^3$$

  ($Q$ vanishing in 2D is the Cayley–Hamilton identity for the $2\times2$ Hessian.) So a $d = 2$ order-6 basis has 9 independent operators where a $d \ge 3$ basis has 11. `multigraph_operators.py` encodes this with the `excluded_2d` flag, and `audit/repro/A2_A6_equivariance_identities.py` verifies both identities hold exactly at $d=2$ and fail at $d=3$.

In short: the two dependent operators are a **$d=2$, order-6** phenomenon, and they do not touch the 20-operator order-4 library used for the results below.

| Library | Form | # terms | $O(2)$-equivariant | Linear | Spanning |
|---|---|---|---|---|---|
| PDE-FIND | $f^p \cdot \partial^a f$ | 60 | ❌ | ❌ | ❌ |
| WSINDy | $\partial^\alpha f^\beta$ | 46 | ❌ | ✅ | ❌ |
| Polynomial jet | monomials in jet variables | 105 (capped) | ❌ | ❌ | ✅ |
| Multigraph | multigraph operators | 20 | ✅ | ❌ | ✅ |
| Modified Multigraph | multigraph operators | 20 | ✅ | ~ | ✅ |

The polynomial-jet count of 105 is for the **capped** enumeration
(``build_polyjet_library`` with ``max_order=4, max_degree=3``).  For
cross-library evaluation the code uses a **dynamic** polynomial-jet library
(``build_polyjet_from_descriptions``) that contains exactly the monomials
needed to represent the operators in all libraries being compared, typically
60–70 columns instead of 105.  See ``library/conversions.py``
(``collect_needed_jet_monomials``) and ``evaluation.py``
(``build_eval_theta_jet``).

Here **Linear** means the library's terms are linear differential operators applied to monomials in $f$ — so integration by parts can move the derivative onto the test function in the weak formulation.

## Multigraph operator library — O(2)-equivariant

`get_operators(max_p, max_n)` returns operators in a **fixed three-level
lexicographic order**:

| Level | Key | Meaning |
|-------|-----|---------|
| 1 (primary) | `p` | total edge count in the multigraph (= half the differential order, i.e. order = 2p) |
| 2 (secondary) | `k` | polynomial degree (total vertex count, non-isolated + isolated) |
| 3 (tertiary) | description string | Unicode lexicographic order — deterministic tie-breaker within each (p, k) group |

This ordering is hardcoded, and is **enforced by `tests/test_multigraph_operators.py`**.  Any change to the operator set or ordering must be accompanied by an update to that test.

For `get_operators(max_p=2, max_n=3)` the 20 operators are, in order:

| Index | p | k | Description |
|-------|---|---|-------------|
| 0  | 0 | 0 | `1` |
| 1  | 0 | 1 | `f` |
| 2  | 0 | 2 | `f²` |
| 3  | 0 | 3 | `f³` |
| 4  | 1 | 1 | `Δf` |
| 5  | 1 | 2 | `f·Δf` |
| 6  | 1 | 2 | `‖∇f‖²` |
| 7  | 1 | 3 | `f²·Δf` |
| 8  | 1 | 3 | `f·‖∇f‖²` |
| 9  | 2 | 1 | `Δ²f` |
| 10 | 2 | 2 | `f·Δ²f` |
| 11 | 2 | 2 | `Δf·Δf` |
| 12 | 2 | 2 | `‖∇²f‖²` |
| 13 | 2 | 2 | `⟨∇Δf,∇f⟩` |
| 14 | 2 | 3 | `f²·Δ²f` |
| 15 | 2 | 3 | `f·Δf·Δf` |
| 16 | 2 | 3 | `f·‖∇²f‖²` |
| 17 | 2 | 3 | `f·⟨∇Δf,∇f⟩` |
| 18 | 2 | 3 | `‖∇f‖²·Δf` |
| 19 | 2 | 3 | `∇²f(∇f,∇f)` |

### Modified multigraph library

We have also introduced a modified version of the multigraph library, which replaces 6 basis elements with elements that behave better under IBP in the weak formulation. The modified library is also O(2)-equivariant, non-redundant, and linearly independent.

For `get_mod_multigraph_operators(max_p=2, max_n=3)` the 20 operators are, in order:

| Index | p | k | Description | Replaces |
|-------|---|---|-------------|----------|
| 0  | 0 | 0 | `1` | — |
| 1  | 0 | 1 | `f` | — |
| 2  | 0 | 2 | `f²` | — |
| 3  | 0 | 3 | `f³` | — |
| 4  | 1 | 1 | `Δf` | — |
| 5  | 1 | 2 | `Δf²` | `f·Δf` |
| 6  | 1 | 2 | `‖∇f‖²` | — |
| 7  | 1 | 3 | `f·‖∇f‖²` | — |
| 8  | 1 | 3 | `Δf³` | `f²·Δf` |
| 9  | 2 | 1 | `Δ²f` | — |
| 10 | 2 | 2 | `Δf·Δf` | — |
| 11 | 2 | 2 | `Δ²f²` | `f·Δ²f` |
| 12 | 2 | 2 | `Δ‖∇f‖²` | `⟨∇Δf,∇f⟩` |
| 13 | 2 | 2 | `‖∇²f‖²` | — |
| 14 | 2 | 3 | `f·Δf·Δf` | — |
| 15 | 2 | 3 | `f·‖∇²f‖²` | — |
| 16 | 2 | 3 | `Δ(f·‖∇f‖²)` | `f·⟨∇Δf,∇f⟩` |
| 17 | 2 | 3 | `Δ²f³` | `f²·Δ²f` |
| 18 | 2 | 3 | `‖∇f‖²·Δf` | — |
| 19 | 2 | 3 | `∇²f(∇f,∇f)` | — |

The 6 replaced operators are those of the form f·(derivative) which are replaced by Δ(something) or Δ²(something) forms that are better suited for integration by parts in the weak formulation.



## Error metric

Recovery is measured in terms of how many false positives (fp) and false negatives (fn) are present in the recovered PDE, as well as the relative L2 error between the recovered PDE and the true PDE, evaluated on the same dataset:

$$E = \frac{\|\Theta (\xi - \xi_{\mathrm{true}})\|_2}{\|\Theta \xi_{\mathrm{true}}\|_2}$$

where $\Theta$ is the library matrix, $\xi$ is the recovered coefficient vector, and $\xi_{\mathrm{true}}$ is the ground-truth coefficient vector. The denominator uses $\Theta \xi_{\mathrm{true}}$ (the true dynamics), *not* the observed $f_t$, so the metric is independent of observation noise in the time derivative.

This is a good metric if the PDE to recover is representable in the library, but it can be misleading if the library is missing key operators. For this reason we convert all terms to a **common polynomial-jet basis** and compute fp/fn counts and $E$ in that basis, so that all libraries are compared on equal footing.

The jet basis is built **dynamically**: ``collect_needed_jet_monomials`` (in ``library/conversions.py``) unions the analytic jet-space decompositions of every operator in every library being compared, and ``build_polyjet_from_descriptions`` (in ``library/polyjet.py``) builds exactly those columns.  This guarantees that every library operator is representable in the evaluation space without silent truncation.  The top-level entry point is ``evaluation.build_eval_theta_jet(jet_monomials=...)``.


# References

[1] Ballerin, F. & Grong, E. (2026). Equivariant PDE discovery. (This paper; currently unpublished.)
<!-- TODO: replace with the arXiv/DOI link and full citation once the PDE-discovery paper is published. -->

[2] Ballerin, F. & Grong, E. (2026). Equivariant nonlinear partial differential operators on constant curvature spaces. *arXiv preprint arXiv:2605.16847*. https://arxiv.org/abs/2605.16847

[3] Samuel H. Rudy et al. ,Data-driven discovery of partial differential equations.Sci. Adv.3,e1602614(2017). DOI:10.1126/sciadv.1602614 (PDE-FIND)

[4] Messenger, D. A. & Bortz, D. M. (2021). Weak SINDy for Partial Differential Equations. *J. Comput. Phys.*, 443, 110525. https://arxiv.org/abs/2007.02848 (WSINDy)

[5] Patra, M. & Karttunen, M. (2006). Stencils with isotropic discretization error for differential operators. *Numer. Methods Partial Differential Eq.*, 22(4), 936–953.