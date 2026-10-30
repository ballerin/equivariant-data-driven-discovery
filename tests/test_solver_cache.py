"""
Tests for ``regression.SolverCache`` and the data-only precomputation it
enables in ``TrainSTRidge``, ``TrainRFE``, ``TrainMSTLS`` and ``centering``.

The cache must never change a result: every solution computed with a cache
(on first use and on later cache hits) has to be bitwise identical to the
solution computed without one.
"""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.regression import (  # noqa: E402
    SolverCache, STRidge, TrainSTRidge, TrainRFE, TrainMSTLS, centering,
    mstls, select_lambda_grid,
)


def _system(n=1500, d=24, k=4, noise=0.05, seed=0, collinear=True):
    """Design matrix with a constant column, badly scaled and (optionally)
    near-duplicate columns, as in weak-form libraries."""
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, d))
    if collinear:
        X[:, d // 2:] = X[:, : d - d // 2] + 1e-3 * rng.standard_normal((n, d - d // 2))
    X *= np.logspace(0, 4, d)[rng.permutation(d)]
    X += 2.0 * rng.standard_normal(d)
    X[:, 0] = 1.0
    w = np.zeros(d)
    w[1 + rng.choice(d - 1, k, replace=False)] = rng.uniform(0.5, 2.0, k)
    y = X @ w
    y = y + noise * np.linalg.norm(y) / np.sqrt(n) * rng.standard_normal(n)
    descs = ["1"] + [f"c{j}" for j in range(1, d)]
    return X, y, descs


LAMS = [1e-6, 1e-2, 1.0]
D_TOLS = [1e-3, 0.1, 10.0]


def _solvers(cache):
    """(label, solve) pairs covering the grids the sweeps use."""
    kw = {} if cache is None else {"cache": cache}
    out = []
    for lam in LAMS:
        for d_tol in D_TOLS:
            out.append((f"stridge lam={lam} d_tol={d_tol}",
                        lambda A, b, lam=lam, d_tol=d_tol:
                        TrainSTRidge(A, b, lam, d_tol, **kw)))
    for lam in [0.0, 1e-4]:
        out.append((f"rfe lam={lam}",
                    lambda A, b, lam=lam: TrainRFE(A, b, lam=lam, **kw)))
    for alpha in [0.25, 0.75]:
        out.append((f"mstls alpha={alpha}",
                    lambda A, b, alpha=alpha: TrainMSTLS(A, b, alpha=alpha, **kw)[0]))
        out.append((f"mstls M alpha={alpha}",
                    lambda A, b, alpha=alpha: TrainMSTLS(
                        A, b, M=np.linalg.norm(A, axis=0) / np.linalg.norm(A),
                        alpha=alpha, **kw)[0]))
    return out


@pytest.mark.parametrize("center", [False, True])
@pytest.mark.parametrize("seed,collinear", [(0, True), (1, False)])
def test_cached_solvers_bit_identical(center, seed, collinear):
    X, y, descs = _system(seed=seed, collinear=collinear)
    cache = SolverCache()

    def run(solve, c):
        if center:
            return centering(X, y, solve, tol=0.05, descs=descs, cache=c)
        return solve(X, y)

    reference = [run(solve, None) for _, solve in _solvers(None)]
    for pass_ in range(2):          # first use, then all cache hits
        for (label, solve), ref in zip(_solvers(cache), reference):
            got = run(solve, cache)
            assert np.array_equal(got, ref), f"{label} (pass {pass_})"


def test_mstls_lambda_and_loss_bit_identical():
    X, y, _ = _system(seed=2)
    cache = SolverCache()
    ref = TrainMSTLS(X, y, alpha=0.5)
    for alpha in (0.25, 0.5, 0.75, 0.5):
        got = TrainMSTLS(X, y, alpha=alpha, cache=cache)
        if alpha == 0.5:
            assert np.array_equal(got[0], ref[0])
            assert got[1] == ref[1]
            assert np.array_equal(got[2], ref[2])


def test_public_single_solves_unchanged():
    """STRidge / mstls / select_lambda_grid keep their signatures and results
    across repeated calls (no hidden state)."""
    X, y, _ = _system(seed=3)
    for f in (lambda: STRidge(X, y, 1e-5, 10, 0.1),
              lambda: mstls(X, y, 0.05),
              lambda: select_lambda_grid(X, y)):
        assert np.array_equal(f(), f())


def test_cache_distinguishes_arrays_and_options():
    X, y, _ = _system(seed=4)
    X2, y2, _ = _system(seed=5)
    cache = SolverCache()
    a = TrainSTRidge(X, y, 1e-5, 0.1, cache=cache)
    b = TrainSTRidge(X2, y2, 1e-5, 0.1, cache=cache)
    assert np.array_equal(b, TrainSTRidge(X2, y2, 1e-5, 0.1))
    # different split / seed must not reuse the first entry
    c = TrainSTRidge(X, y, 1e-5, 0.1, seed=1, cache=cache)
    assert np.array_equal(c, TrainSTRidge(X, y, 1e-5, 0.1, seed=1))
    assert np.array_equal(a, TrainSTRidge(X, y, 1e-5, 0.1, cache=cache))


def test_cache_does_not_alias_returned_vectors():
    """Mutating a returned solution must not corrupt later cached results."""
    X, y, _ = _system(seed=6)
    cache = SolverCache()
    ref = TrainSTRidge(X, y, 1e-6, 1e6)       # huge d_tol: baseline LS wins
    w = TrainSTRidge(X, y, 1e-6, 1e6, cache=cache)
    w[:] = 0.0
    assert np.array_equal(TrainSTRidge(X, y, 1e-6, 1e6, cache=cache), ref)
