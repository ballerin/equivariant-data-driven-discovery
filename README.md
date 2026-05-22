# Equivariant PDE discovery

## Summary
This directory contains the core Python code used to
- construct an O(2)-equivariant library of nonlinear partial differential operators;
- build the standard coordinate library for comparison;
- simulate benchmark PDEs;
- run sparse regression(STRidge) experiments to recover PDEs from data.

The implementation provides two derivative backends (spectral and local-polynomial) so experiments can trade accuracy for noise-robustness.

## Motivation
The code implements the multigraph-operator approach from the paper: a minimal, complete, rotation-equivariant basis of polynomial nonlinear differential operators. It is organised to make it straightforward to (re)generate datasets, assemble library matrices, and run hyperparameter sweeps that compare the equivariant basis against the standard PDE-FIND style coordinate library.

## Repository structure (brief)
- `derivatives.py`        : Backend selector and context manager for derivative methods.
- `spectral_derivatives.py`: FFT-based derivative routines (accurate for periodic, clean data).
- `polynomial_derivatives.py`: Local least-squares polynomial derivative kernels (noise-robust).
- `multigraph_operators.py`: Construction of the O(2)-equivariant multigraph operators and
	the public `get_operators(max_p, max_n)` API.
- `library.py`            : Builders for the equivariant (`build_equivariant_library`) and
	standard (`build_standard_library`) design matrices used for regression.
- `regression.py`         : STRidge implementation and `TrainSTRidge` cross-validation helper.
- `simulations.py`        : Spectral/IMEX simulation drivers for benchmark PDEs (Allen–Cahn, Swift–Hohenberg, Cahn–Hilliard, Kuramoto–Sivashinsky, Porous Medium), plus `make_regression_data` for dataset preparation.
- `utils.py`              : Small helpers (`print_pde`, `rotation_90`, etc.).
- `scripts/hyperparam_sweep.py`: Script to sweep `(lambda, d_tol)` for STRidge and save results to `results/hyperparam_sweep.csv`.
- `PDE_FIND_DATA/`        : Example/legacy datasets used for comparison (MAT files, helper scripts).
- `notebooks/`            : Jupyter notebooks with interactive experiments and visualisations.
- `results/`              : Output from scripts (CSV, figures).

## Assumptions
- 2D domain with periodic boundary conditions