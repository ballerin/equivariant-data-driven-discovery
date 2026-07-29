"""Cahn-Hilliard equation simulator (spectral IMEX), 2-D and n-D variants."""

import numpy as np

from ._imex import _run


def cahn_hilliard(
    N=128,
    L=8 * np.pi,
    T_end=30.0,
    eps=0.2,
    dt=0.005,
    save_every=30,
    f0=None,
    seed=42,
    obs_noise=0.0,
    dyn_noise=0.0,
):
    """
    Simulate the Cahn-Hilliard equation on [0, L]^2 via spectral IMEX.

    PDE:
        f_t = -eps^2 * Delta^2 f + Delta(f^3 - f)

    IMEX scheme (stabilised: linear bilaplacian implicit):
        f̂^{n+1} = (f̂^n + dt * (-|k|^2) * F̂(f^3 - f)) / (1 + eps^2 |k|^4 dt)

    Parameters
    ----------
    (same as allen_cahn, with eps = interface width)

    Returns
    -------
    snapshots : ndarray of shape (n_saved, N, N)
    times : ndarray of shape (n_saved,)
    """
    rng = np.random.default_rng(seed)

    k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    K2 = KX**2 + KY**2
    K4 = K2**2

    denom = 1.0 + eps**2 * K4 * dt

    if f0 is None:
        f0 = 0.5 * rng.standard_normal((N, N))

    # Explicit part: Delta(f^3 - f)  -->  in Fourier: -|k|^2 * FFT(f^3 - f)
    def nonlin(f):
        # return spatial-domain representation of -K2 * fft(f^3-f) → go via iFFT
        # but _run calls fft2(nonlin(f)), so nonlin must return a *real-space* field
        # We compute Delta(f^3 - f) spectrally and return the real-space result.
        return np.real(np.fft.ifft2(-K2 * np.fft.fft2(f**3 - f)))

    n_steps = int(round(T_end / dt)) + 1  # +1 to ensure we compute also time T_end
    snapshots, times = _run(f0, nonlin, denom, dt, n_steps, save_every, dyn_noise, rng)

    if obs_noise > 0.0:
        sigma = obs_noise * snapshots.std()
        snapshots = snapshots + sigma * rng.standard_normal(snapshots.shape)

    return snapshots, times


_MAX_NDIM = 8


def cahn_hilliard_nd(
    ndim=2,
    N=None,
    L=8 * np.pi,
    T_end=30.0,
    eps=0.2,
    dt=0.005,
    save_every=30,
    f0=None,
    seed=42,
    obs_noise=0.0,
    dyn_noise=0.0,
):
    """Simulate the Cahn-Hilliard equation on [0, L]^ndim via spectral IMEX.

    PDE:
        f_t = -eps^2 * Delta^2 f + Delta(f^3 - f)

    IMEX scheme (stabilised: linear bilaplacian implicit):
        f^{n+1} = (f^n + dt * (-|k|^2) * F(f^3 - f)) / (1 + eps^2 |k|^4 dt)

    Parameters
    ----------
    ndim : int
        Number of spatial dimensions (1 <= ndim <= 8).
    N : int or None
        Grid points per direction.  Default: 128 for ndim=2, 48 for ndim>=3.
    L, T_end, eps, dt, save_every, f0, seed, obs_noise, dyn_noise :
        Same as ``cahn_hilliard``.
    """
    if not 1 <= ndim <= _MAX_NDIM:
        raise ValueError(f"ndim must be 1-{_MAX_NDIM}, got {ndim}")

    if N is None:
        N = 128 if ndim == 2 else 48

    rng = np.random.default_rng(seed)

    k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
    K_grids = np.meshgrid(*([k] * ndim), indexing="ij")
    K2 = sum(Ki**2 for Ki in K_grids)
    K4 = K2**2

    denom = 1.0 + eps**2 * K4 * dt

    if f0 is None:
        f0 = 0.5 * rng.standard_normal((N,) * ndim)

    # Explicit: Delta(f^3 - f) = -K2 * FFT(f^3 - f) in Fourier, returned in real space
    def nonlin(f):
        return np.real(np.fft.ifftn(-K2 * np.fft.fftn(f**3 - f)))

    n_steps = int(round(T_end / dt)) + 1
    snapshots, times = _run(f0, nonlin, denom, dt, n_steps, save_every, dyn_noise, rng)

    if obs_noise > 0.0:
        sigma = obs_noise * snapshots.std()
        snapshots = snapshots + sigma * rng.standard_normal(snapshots.shape)

    return snapshots, times
