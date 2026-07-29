"""Swift-Hohenberg equation simulator (spectral IMEX), 2-D and n-D variants."""

import numpy as np

from ._imex import _run


def swift_hohenberg(
    N=64,
    L=4 * np.pi,
    T_end=200.0,
    r=0.3,
    dt=0.1,
    save_every=50,
    f0=None,
    seed=42,
    obs_noise=0.0,
    dyn_noise=0.0,
):
    """
    Simulate the Swift-Hohenberg equation on [0, L]^2 via spectral IMEX.

    PDE:
        f_t = (r-1) f - 2 Delta f - Delta^2 f - f^3

    IMEX scheme (linear-implicit / nonlinear-explicit):
        f̂^{n+1} = (f̂^n + dt * N̂(f^n)) / (1 - dt * L_k)

    Parameters are analogous to :func:`allen_cahn`.
    """
    rng = np.random.default_rng(seed)

    k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    K2 = KX**2 + KY**2
    K4 = K2**2

    # Linear Fourier symbol: (r-1) + 2 K2 - K4
    Lk = (r - 1.0) + 2.0 * K2 - K4
    denom = 1.0 - dt * Lk

    if f0 is None:
        f0 = 0.01 * rng.standard_normal((N, N))

    def nonlin(f):
        # cubic reaction term
        return -(f**3)

    n_steps = int(round(T_end / dt)) + 1
    snapshots, times = _run(f0, nonlin, denom, dt, n_steps, save_every, dyn_noise, rng)

    if obs_noise > 0.0:
        sigma = obs_noise * snapshots.std()
        snapshots = snapshots + sigma * rng.standard_normal(snapshots.shape)

    return snapshots, times


_MAX_NDIM = 8


def swift_hohenberg_nd(
    ndim=2,
    N=None,
    L=4 * np.pi,
    T_end=200.0,
    r=0.3,
    dt=0.1,
    save_every=50,
    f0=None,
    seed=42,
    obs_noise=0.0,
    dyn_noise=0.0,
):
    """Simulate the Swift-Hohenberg equation on [0, L]^ndim via spectral IMEX.

    PDE:
        f_t = (r-1) f - 2 Delta f - Delta^2 f - f^3

    Parameters
    ----------
    ndim : int
        Number of spatial dimensions (1 <= ndim <= 8).
    N : int or None
        Grid points per direction.  Default: 64 for ndim=2, 48 for ndim>=3.
    L, T_end, r, dt, save_every, f0, seed, obs_noise, dyn_noise :
        Same as ``swift_hohenberg``.
    """
    if not 1 <= ndim <= _MAX_NDIM:
        raise ValueError(f"ndim must be 1-{_MAX_NDIM}, got {ndim}")

    if N is None:
        N = 64 if ndim == 2 else 48

    rng = np.random.default_rng(seed)

    k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
    K_grids = np.meshgrid(*([k] * ndim), indexing="ij")
    K2 = sum(Ki**2 for Ki in K_grids)
    K4 = K2**2

    # Linear Fourier symbol: (r-1) + 2 K2 - K4
    Lk = (r - 1.0) + 2.0 * K2 - K4
    denom = 1.0 - dt * Lk

    if f0 is None:
        f0 = 0.01 * rng.standard_normal((N,) * ndim)

    def nonlin(f):
        # cubic reaction term
        return -(f**3)

    n_steps = int(round(T_end / dt)) + 1
    snapshots, times = _run(f0, nonlin, denom, dt, n_steps, save_every, dyn_noise, rng)

    if obs_noise > 0.0:
        sigma = obs_noise * snapshots.std()
        snapshots = snapshots + sigma * rng.standard_normal(snapshots.shape)

    return snapshots, times
