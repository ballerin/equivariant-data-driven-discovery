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

#### 2. Run the main ablation sweep

Full sweep (all libraries, formulations, backends, and regression algorithms).

```bash
python scripts/sweep.py
```

Targeted sweeps for faster iteration:

```bash
# Multigraph library vs PDE-FIND library comparison (pointwise only, STRidge only)
python scripts/sweep.py --methods pointwise --algorithms stride --libraries multigraph pdefind

# Noise-robustness comparison (weak vs pointwise, multigraphs only)
python scripts/sweep.py \
    --equations sh ks_chaos \
    --noise 0.02 0.05 \
    --methods pointwise weak_direct weak_ibp \
    --libraries multigraph

# With post-hoc noise injection and spectral low-pass prefiltering
python scripts/sweep.py --post-hoc-noise 0.02 0.05 --spectral-cutoff 3.0

# MSTLS-only with custom alpha sweep
python scripts/sweep.py --algorithms mstls --mstls-alpha 0.0 0.25 0.5

# Polynomial jet library only
python scripts/sweep.py --libraries polyjet

# WSINDy library only
python scripts/sweep.py --libraries wsindy
```

Key CLI options for `sweep.py`:

| Flag | Default | Description |
|---|---|---|
| `--methods` | `pointwise weak_direct weak_ibp` | Library evaluation methods |
| `--algorithms` | `stride mstls rfe` | Regression algorithms (`stride`, `mstls`, `rfe`, `brute`, `oracle`) |
| `--libraries` | `multigraph pdefind polyjet wsindy` | Library types (`multigraph`, `mod_multigraph`, `pdefind`, `polyjet`, `wsindy`) |
| `--backends` | `spectral polynomial` | Differentiation backends |
| `--equations` | all | Filter by equation name prefix |
| `--noise` | all | Filter by noise level (e.g. `0.0 0.02 0.05`) |
| `--weak-m` | `4 8 12` | Test-function half-support sizes |
| `--post-hoc-noise` | none | Add noise to clean data as fraction of std |
| `--spectral-cutoff` | none | Low-pass filter cutoff wavenumber |
| `--name` | `""` | Tag appended to output filename |

Output: `results/sweep_<YYYYMMDD>_<HHMMSS>_<name>.csv`

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

#### 4. Extract best results

**Overall best** (one row per equation + noise + method + library + backend group):

```bash
python scripts/summarize_best.py results/sweep_20260612_120000.csv
```

This writes `results/sweep_20260612_120000_best.csv` and prints a formatted table.

**Per-equation best** (grouped by equation only — the default `--group-by`):

```bash
python scripts/summarize_best.py results/sweep_20260612_120000.csv \
    --group-by equation
```

Custom grouping (e.g. by equation + noise level):

```bash
python scripts/summarize_best.py results/sweep_20260612_120000.csv \
    --group-by equation noise
```

**Process multiple sweep files at once:**

```bash
python scripts/summarize_best.py results/sweep_*.csv
```

The lexicographic criterion (exact match → fp → fn → $E$) is used to select the
single best row within each group.

#### 5. Build per-equation markdown summaries

Generate a human-readable markdown summary for each equation from a sweep CSV:

```bash
# Build summaries for all equations in the CSV
python scripts/build_summary.py results/sweep_20260619_170245.csv

# Restrict to a single equation
python scripts/build_summary.py results/sweep_20260619_170245.csv --equation pme

# Specify output path
python scripts/build_summary.py results/sweep_20260619_170245.csv --out summaries/my_summary.md
```

Output: `summaries/<equation>.md` by default.

### Results (`results/`)

Sweep outputs are CSV files with naming convention `sweep_<YYYYMMDD>_<HHMMSS>_<config>.csv`. A companion `_best.csv` file is generated by `summarize_best.py` containing only the best row per (PDE, noise) group.

### Summaries (`summaries/`)

Per-PDE markdown summaries of best recovery results, one per equation. Also includes a `README.md` summary index.

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

Notice that in general the terms in the library are not linearly independent. In particular, the equivariant library has 20 operators but only 18 are linearly independent as shown in [1].

| Library | Form | # terms | $O(2)$-equivariant | Linear | Spanning |
|---|---|---|---|---|---|
| PDE-FIND | $f^p \cdot \partial^a f$ | 60 | ❌ | ❌ | ❌ |
| WSINDy | $\partial^\alpha f^\beta$ | 46 | ❌ | ✅ | ❌ |
| Polynomial jet | monomials in jet variables | 105 | ❌ | ❌ | ✅ |
| Multigraph | multigraph operators | 20 | ✅ | ❌ | ✅ |
| Modified Multigraph | multigraph operators | 20 | ✅ | ~ | ✅ |

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

This is a good metric if the PDE to recover is representable in the library, but it can be misleading if the library is missing key operators. For this reason we convert all terms to polynomial jet form and compute fp/fn counts and $E$ in that basis, so that all libraries are compared on equal footing.


# References

[1] Ballerin, F. & Grong, E. (2026). Equivariant PDE discovery. (This paper; currently unpublished.)
<!-- TODO: replace with the arXiv/DOI link and full citation once the PDE-discovery paper is published. -->

[2] Ballerin, F. & Grong, E. (2026). Equivariant nonlinear partial differential operators on constant curvature spaces. *arXiv preprint arXiv:2605.16847*. https://arxiv.org/abs/2605.16847

[3] Samuel H. Rudy et al. ,Data-driven discovery of partial differential equations.Sci. Adv.3,e1602614(2017). DOI:10.1126/sciadv.1602614 (PDE-FIND)

[4] Messenger, D. A. & Bortz, D. M. (2021). Weak SINDy for Partial Differential Equations. *J. Comput. Phys.*, 443, 110525. https://arxiv.org/abs/2007.02848 (WSINDy)

[5] Patra, M. & Karttunen, M. (2006). Stencils with isotropic discretization error for differential operators. *Numer. Methods Partial Differential Eq.*, 22(4), 936–953.