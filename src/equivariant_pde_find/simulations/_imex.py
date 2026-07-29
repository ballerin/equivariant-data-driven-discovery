"""Generic IMEX (implicit-linear / explicit-nonlinear) spectral time-stepper.

Shared by allen_cahn, allen_cahn_nd, cahn_hilliard, cahn_hilliard_nd,
swift_hohenberg, swift_hohenberg_nd, and kuramoto_sivashinsky's legacy
``scheme="imex"`` path. Works for any number of spatial dimensions via
``fftn``/``ifftn``, unlike ``_etdrk4.py`` which is 2-D only.
"""

import numpy as np


def _imex_step(f, nonlin_fn, denom, dt, dyn_noise, rng):
    """Single IMEX step: implicit linear + explicit nonlinear + optional noise.

    Supports any ndim >= 1 via fftn/ifftn.
    """
    fhat = np.fft.fftn(f) + dt * np.fft.fftn(nonlin_fn(f))
    if dyn_noise > 0.0:
        xi = rng.standard_normal(f.shape)
        fhat += np.sqrt(dt) * dyn_noise * np.fft.fftn(xi)
    return np.real(np.fft.ifftn(fhat / denom))


def _run(f0, nonlin_fn, denom, dt, n_steps, save_every, dyn_noise, rng):
    """Run the time loop, returning (snapshots, times).

    A snapshot saved when ``step % save_every == 0`` is the state *after*
    that iteration's update, i.e. after ``step + 1`` total steps -- so it
    is labeled with time ``(step + 1) * dt``, not ``step * dt``.
    """
    f = f0.copy()
    snapshots, steps = [], []
    for step in range(n_steps):
        f = _imex_step(f, nonlin_fn, denom, dt, dyn_noise, rng)
        if step % save_every == 0:
            snapshots.append(f.copy())
            steps.append(step + 1)
    return np.array(snapshots), np.array(steps, dtype=float) * dt
