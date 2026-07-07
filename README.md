# Equivariant PDE discovery

This repository is intended as reference for the paper "Equivariant PDE discovery" by Ballerin & Grong [1] for discovery of partial differential equations in a 2D flat domain with periodic boundary conditions. It contains the core Python code used to:
- simulate benchmark PDEs, optionally with noise;
- construct an O(2)-equivariant library of nonlinear partial differential operators;
- build the standard coordinate libraries (non-equivariant) for comparison;
- apply sparse regression to recover PDEs from data;


The implementation provides several differentiation backends, library building strategies and regression algorithms, so experiments can trade accuracy for noise-robustness and an ablation study can be performed. Details on the implemented methods can be found below.

## Setup

Install in editable mode:

    pip install -e .

This makes the ``equivariant_pde_find`` package importable from anywhere.
Scripts add ``src/`` to ``sys.path`` so they also work without a pip install.

## Motivation
The code exploits the complete classification of the multigraph-operators from Ballerin & Grong [2]: a minimal, complete, $O(2)$-equivariant basis of polynomial nonlinear partial differential operators. The general aim of the research is to show that whenever the underlying dynamics has the additional assumption of $O(2)$-equivariance, introducing bias in the form of a complete and non-redundant library of multigraph operators yields better recovery when deducing the generating PDE.

The code is organised to make it straightforward to (re)generate datasets, assemble library matrices, and run hyperparameter sweeps that compare the equivariant basis from the paper against standard non-equivariant libraries such as those used in SINDy/PDE-FIND  and WSINDy.

## What is implemented

The repository implements the following strategies, of which the first one is our novel contribution, and the other two are benchmarks from prior published work:
- $O(2)$-equivariant library based on multigraph operators;
- PDE-FIND library (based on SINDy algorithm) [3];
- WSINDy library and algorithm [4].

Two differentiation backends are available:
- **Spectral** (FFT-based): exact for band-limited data, but highly susceptible to measurement noise.
- **Polynomial stencils** (local least-squares): fit a 2-D polynomial on a $(2r+1)\times(2r+1)$ stencil at each grid point. More noise-robust than spectral at the cost of a small approximation bias.

Several regression algorithms:
- STridge, originally from PDE-FIND [3];
- SR3, TODO.
- MSTLS, TODO.

## Repository structure

```
pyproject.toml                  # Package metadata, dependencies, tool config
README.md
src/equivariant_pde_find/       # Core library (installable package)
    __init__.py
    multigraph_operators.py     # O(2)-equivariant operator basis
    derivatives.py              # Backend selector (spectral ↔ polynomial)
    spectral_derivatives.py     # FFT-based differentiation
    polynomial_derivatives.py   # Local least-squares stencil differentiation
    library.py                  # Design-matrix builders (pointwise)
    weak_form.py                # Weak-form library builders & test functions
    regression.py               # Sparse regression (STRidge, SR3, MSTLS)
    evaluation.py               # Recovery metrics & sweep selection
    simulations.py              # PDE simulators (IMEX / ETDRK4)
    equation_configs.py         # Ground-truth coefficients & equation registry
    utils.py                    # Small utilities (print_pde, rotation_90, filter)
data/                           # Pre-generated .npz datasets + manifest.csv
scripts/
    generate_data.py            # Generate & save benchmark datasets
    sweep.py                    # Canonical ablation sweep (pointwise / weak / IBP)
    sweep_polynomial_stencil.py # Stencil parameter sweep (degree, radius)
    ks_regimes_sweep.py         # KS regime exploration (domain size sweep)
    summarize_best.py           # Parse sweep CSVs → best-per-config tables
notebooks/                      # Demo notebooks, one per benchmark PDE
    allen_cahn.ipynb
    cahn_hilliard.ipynb
    kuramoto_sivashinsky.ipynb
    nonlinear_diffusion.ipynb
    swift_hohenberg.ipynb
results/                        # Sweep output CSVs (auto-generated)
summaries/                      # Per-PDE markdown summaries of best results
tests/                          # pytest test suite
```

### Core library (`src/equivariant_pde_find/`)

| Module | Purpose |
|---|---|
| `multigraph_operators.py` | O(2)-equivariant multigraph operator basis from [2]. Implements all connected multigraph operators for $p \le 3$ edges and $n \le 4$ nodes. Exports `get_operators(max_p, max_n) → list[(func, desc)]`. |
| `derivatives.py` | Backend selector. Exposes six neutral operations (`grad`, `laplacian`, `hessian`, `bilaplacian`, `trilaplacian`, `third_order`) that dispatch to the active backend. Use `set_backend('spectral'|'polynomial')` or the `use_backend(...)` context manager. |
| `spectral_derivatives.py` | FFT-based differentiation. Exact for band-limited data but highly susceptible to noise. |
| `polynomial_derivatives.py` | Local least-squares polynomial fit differentiation. Fits a 2-D polynomial of degree `deg` on a $(2r+1)\times(2r+1)$ stencil and extracts derivatives from the coefficients. More noise-robust than spectral at the cost of a small bias. Default stencil parameters are tuned per derivative order (e.g. order-2 uses deg=9, r=5 for <1e-5 relative error at typical wavenumbers). Based on the isotropic stencil construction method of Patra & Karttunen [5]. |
| `library.py` | Pointwise library construction. Builds the design matrix $\Theta$ from field snapshots. Five builders: `build_equivariant_library` (multigraph basis), `build_standard_library` (PDE-FIND: $f^p \cdot \partial^a f$), `build_pointwise_wsindy_library` ($\partial^\alpha(f^\beta)$), `build_polynomial_jet_library` (monomials in jet variables), and `build_standard_library_1d` (1-D variant). Also provides `spectral_lowpass` for isotropic Fourier filtering. |
| `weak_form.py` | Weak-form PDE discovery via convolution with compactly-supported polynomial-bump test functions. Five builders: `build_weak_library` (equivariant), `build_weak_standard_library` (PDE-FIND), `build_weak_wsindy_library` ($\partial^\alpha(f^q)$), `build_weak_polynomial_jet_library` (monomials in jet variables), and `build_weak_rhs` (weak RHS only). Includes integration by parts for linear operators, automatic test-function parameter selection via the changepoint method [4], and PDE scale-invariance rescaling ($\gamma_x, \gamma_t, \gamma_u$) for improved conditioning. |
| `regression.py` | Sparse regression algorithms: `STRidge` (sequential threshold ridge regression), `TrainSTRidge` (STRidge with CV over tolerance), `SR3` (sparse relaxed regularized regression), `TrainSR3` (SR3 with CV over L₁ penalty), `mstls` (modified sequential thresholding least squares), `TrainMSTLS` (MSTLS with λ sweep). |
| `evaluation.py` | Recovery evaluation. `build_eval_theta` pre-builds a universal pointwise Θ for fast repeated evaluation. `evaluate_from_theta` computes fp, fn, and relative $L_2$ error $E$. `best_of_sweep` selects the best row using the lexicographic criterion (exact match → fp → fn → $E$). Also supports standard-library evaluation via `evaluate_standard` / `evaluate_standard_with_truth` and operator decomposition via `decompose_equivariant_to_standard`. |
| `simulations.py` | PDE simulators for benchmark datasets. Uses IMEX (implicit-explicit) and ETDRK4 time-steppers. Supports observation noise (additive Gaussian post-simulation) and dynamical noise (Euler-Maruyama stochastic forcing). Generates $(t, x, y)$ snapshots for five PDEs: `allen_cahn`, `swift_hohenberg`, `kuramoto_sivashinsky` (three regimes: cellular / weak-chaotic / chaotic), `cahn_hilliard`, and `porous_medium`. Also provides `make_regression_data` for spatial subsampling. |
| `equation_configs.py` | Shared equation metadata. Provides ground-truth coefficient vectors (`true_allen_cahn`, `true_swift_hohenberg`, `true_ks`, `true_cahn_hilliard`, `true_pme`), a `find_op` helper for looking up operator indices by name, spatial masks (e.g. `mask_pme` for porous medium), and the `EQUATION_META` registry mapping equation names to (truth_fn, skip_fraction, mask_fn) tuples. |
| `utils.py` | Small utilities: `print_pde` (pretty-print recovered PDE from coefficient vector), `rotation_90` (rotate field 90° CCW for equivariance tests), `_apply_prefilter` (periodic Gaussian prefilter shared by both differentiation backends). |

### Scripts (`scripts/`)

| Script | Purpose |
|---|---|
| `generate_data.py` | Generate benchmark datasets (snapshots + times) for all PDEs and noise levels (0%, 2%, 5%). Saves `.npz` files to `data/` with a `manifest.csv` index. Supports `--equations` and `--noise` filters. |
| `sweep.py` | Canonical ablation sweep. Compares equivariant vs standard libraries across PDEs, noise levels, differentiation backends, methods (pointwise / weak-direct / weak-IBP), and regression algorithms (STRidge / SR3 / MSTLS). Supports optional post-hoc noise injection and spectral low-pass preprocessing. Outputs timestamped CSVs to `results/`. |
| `sweep_polynomial_stencil.py` | Sweep over polynomial-derivative stencil parameters (degree, radius) to study accuracy-noise tradeoffs. |
| `ks_regimes_sweep.py` | Kuramoto-Sivashinsky regime exploration across domain sizes (cellular → weak-chaotic → chaotic). |
| `summarize_best.py` | Parse one or more sweep CSV files and print the best row per (PDE, noise) configuration using the lexicographic criterion. |
| `build_summary.py` | Build a per-equation markdown summary file from a sweep CSV, with representability tables, equivariant-vs-standard comparison tables, and detailed pointwise/weak-form results. |

### Data (`data/`)

Pre-generated simulation snapshots stored as `.npz` files, one per (equation, noise level) pair. Each `.npz` contains `snapshots` (ndarray, shape `(n_saved, N, N)`) and `times` (ndarray, shape `(n_saved,)`). Noise levels: `clean` (0%), `noise2pct` (2%), `noise5pct` (5%). The `manifest.csv` file indexes all available datasets.

| Equation | Simulator | # Datasets |
|---|---|---|
| Allen-Cahn | `allen_cahn` | 3 (clean, 2%, 5%) |
| Swift-Hohenberg | `swift_hohenberg` | 3 |
| Cahn-Hilliard | `cahn_hilliard` | 3 |
| Porous medium | `porous_medium` | 3 |
| KS cellular | `kuramoto_sivashinsky` | 3 |
| KS weak-chaotic | `kuramoto_sivashinsky` | 3 |
| KS chaotic | `kuramoto_sivashinsky` | 3 |

### Notebooks (`notebooks/`)

Demonstration notebooks with sample recovery, one per benchmark PDE:

| Notebook | PDE |
|---|---|
| `allen_cahn.ipynb` | Allen-Cahn: $f_t = f - f^3 + \nu \Delta f$, $\nu = \varepsilon^2$ |
| `cahn_hilliard.ipynb` | Cahn-Hilliard: $f_t = -\Delta^2 f + \Delta(f^3 - f)$ |
| `kuramoto_sivashinsky.ipynb` | Kuramoto-Sivashinsky: $f_t = -\Delta f - \Delta^2 f - \\|\nabla f\\|^2$ |
| `nonlinear_diffusion.ipynb` | Nonlinear diffusion / porous medium: $f_t = \nabla\cdot(f^2\nabla f)$ |
| `swift_hohenberg.ipynb` | Swift-Hohenberg: $f_t = -(1+\Delta)^2 f + \alpha f - f^3$ |

### Pipeline

The recommended workflow has four stages: generate data, run all sweeps,
extract the best results (overall), then extract best results per equation.

#### 1. Generate simulation data

Generate all 21 datasets (7 equations × 3 noise levels):

```bash
python scripts/generate_data.py
```

To generate only specific equations or noise levels:

```bash
python scripts/generate_data.py --equations allen_cahn swift_hohenberg
python scripts/generate_data.py --noise 0.0 0.02
```

Output: `.npz` files in `data/` plus a `data/manifest.csv` index.

#### 2. Run the main ablation sweep

Full sweep (all four libraries, pointwise + weak-direct + weak-IBP,
STRidge + SR3 + MSTLS, spectral + polynomial backends, all equations,
all noise levels):

```bash
python scripts/sweep.py
```

Targeted sweeps for faster iteration:

```bash
# Core equivariant vs standard comparison (pointwise only, STRidge only)
python scripts/sweep.py --methods pointwise --algorithms strridge --libraries equivariant standard

# Noise-robustness comparison (weak vs pointwise, equivariant only)
python scripts/sweep.py \
    --equations swift_hohenberg ks_chaotic \
    --noise 0.02 0.05 \
    --methods pointwise weak_direct weak_ibp \
    --libraries equivariant

# With post-hoc noise injection and spectral low-pass prefiltering
python scripts/sweep.py --post-hoc-noise 0.02 0.05 --spectral-cutoff 3.0

# MSTLS-only with custom alpha sweep
python scripts/sweep.py --algorithms mstls --mstls-alpha 0.0 0.25 0.5

# Polynomial jet library (monomials in jet variables, 105 columns)
python scripts/sweep.py --libraries polynomial_jet

# WSINDy library (∂^α(f^q), 46 columns)
python scripts/sweep.py --libraries wsindy

# Cahn-Hilliard epsilon ablation
python scripts/sweep.py --equations cahn_hilliard --ch-epsilons 0.1 0.2 0.4
```

Key CLI options for `sweep.py`:

| Flag | Default | Description |
|---|---|---|
| `--methods` | `pointwise weak_direct weak_ibp` | Library evaluation methods |
| `--algorithms` | `strridge sr3 mstls` | Regression algorithms (`strridge`, `sr3`, `mstls`) |
| `--libraries` | `equivariant standard polynomial_jet wsindy` | Library types (`equivariant`, `standard`, `polynomial_jet`, `wsindy`) |
| `--backends` | `spectral polynomial` | Differentiation backends |
| `--equations` | all | Filter by equation name prefix |
| `--noise` | all | Filter by noise level (e.g. `0.0 0.02 0.05`) |
| `--weak-m` | `4 8 12` | Test-function half-support sizes |
| `--post-hoc-noise` | none | Add noise to clean data as fraction of std |
| `--spectral-cutoff` | none | Low-pass filter cutoff wavenumber |
| `--ch-epsilons` | `0.1 0.2 0.4` | Cahn-Hilliard ε values to sweep |
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
python scripts/build_summary.py results/sweep_20260619_170245.csv --equation porous_medium

# Specify output path
python scripts/build_summary.py results/sweep_20260619_170245.csv --out summaries/my_summary.md
```

Output: `summaries/<equation>.md` by default.  Each file contains:
- A representability table (which terms each library can express)
- A summary table comparing equivariant vs best non-equivariant method with a
  5-level colour indicator (🟢 clear equivariant win → 🔴 clear non-equivariant win)
- Detailed pointwise and weak-form tables for both backends

### Results (`results/`)

Sweep outputs are CSV files with naming convention `sweep_<YYYYMMDD>_<HHMMSS>_<config>.csv`. A companion `_best.csv` file is generated by `summarize_best.py` containing only the best row per (PDE, noise) group.

### Summaries (`summaries/`)

Per-PDE markdown summaries of best recovery results, one per equation. Also includes a `README.md` summary index.

## Error metric

Recovery is measured in terms of how many false positives (fp) and false negatives (fn) are present in the recovered PDE, as well as the relative L2 error between the recovered PDE and the true PDE, evaluated on the same dataset:

$$E = \frac{\|\Theta (\xi - \xi_{\mathrm{true}})\|_2}{\|\Theta \xi_{\mathrm{true}}\|_2}$$

where $\Theta$ is the library matrix, $\xi$ is the recovered coefficient vector, and $\xi_{\mathrm{true}}$ is the ground-truth coefficient vector. The denominator uses $\Theta \xi_{\mathrm{true}}$ (the true dynamics), *not* the observed $f_t$, so the metric is independent of observation noise in the time derivative.

This is a good metric if the PDE to recover is representable in the library, but it can be misleading if the library is missing key operators. In that case, a low E_fit (fit error) can be achieved by using many spurious terms to approximate the missing dynamics, which is what happens sometimes.

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

### Overlap between PDE-FIND and WSINDy

Both libraries share the pure monomials ($1$, $f$, $f^2$, $f^3$) and the pure
derivatives ($\partial_x f$, $\partial_{xx} f$, …, $\partial_{xxxx} f$).
The structural difference is in how they handle nonlinear derivative terms:

* **PDE-FIND unique** (42 columns): $f^p \cdot \partial^a f$ for $p \ge 1$.
  Multiplies a monomial by a **single** derivative factor.  Simple and
  interpretable, but cannot represent products of derivatives like
  $\|\nabla f\|^2 = f_x^2 + f_y^2$.

* **WSINDy unique** (28 columns): $\partial^a(f^q)$ for $q \ge 2$.
  Differentiates a monomial power.  By the chain rule these expand to
  **mixtures** of terms — e.g.\ $\partial_x^2(f^2) = 2f_x^2 + 2f f_{xx}$
  bundles the product-of-derivatives term $f_x^2$ with the PDE-FIND-style
  term $f f_{xx}$.  WSINDy can *reach* $f_x^2$ but cannot isolate it as a
  sparse column.

Consequently, **neither** standard library can represent $\|\nabla f\|^2$ as an independent sparse term. For this reason the **polynomial jet** library has been introduced.

## Equivariant operator library — canonical ordering

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


# References

[1] Ballerin, F. & Grong, E. (2026). Equivariant PDE discovery. (This paper.)

[2] Ballerin, F. & Grong, E. (2026). Equivariant nonlinear partial differential operators on constant curvature spaces. *arXiv preprint arXiv:2605.16847*. https://arxiv.org/abs/2605.16847

[3] Samuel H. Rudy et al. ,Data-driven discovery of partial differential equations.Sci. Adv.3,e1602614(2017). DOI:10.1126/sciadv.1602614 (PDE-FIND)

[4] Messenger, D. A. & Bortz, D. M. (2021). Weak SINDy for Partial Differential Equations. *J. Comput. Phys.*, 443, 110525. https://arxiv.org/abs/2007.02848 (WSINDy)

[5] Patra, M. & Karttunen, M. (2006). Stencils with isotropic discretization error for differential operators. *Numer. Methods Partial Differential Eq.*, 22(4), 936–953.