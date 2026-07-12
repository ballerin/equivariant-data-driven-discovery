import numpy as np

from ._imex import _run


def allen_cahn(
    N=64,
    L=2 * np.pi,
    T_end=8.0,
    eps=0.15,
    dt=0.01,
    save_every=20,
    f0=None,
    seed=42,
    obs_noise=0.0,
    dyn_noise=0.0,
):
    """
    Simulate the Allen-Cahn equation on [0, L]^2 via spectral IMEX.

    PDE:
        f_t = eps^2 * Delta f + f - f^3

    IMEX scheme (linear-implicit / nonlinear-explicit):
        f̂^{n+1} = (f̂^n + dt * F̂(f^n)) / (1 + eps^2 |k|^2 dt)

    Parameters
    ----------
    N : int
        Grid points per direction (power of 2 recommended).
    L : float
        Domain side length (periodic, [0, L]^2).
    T_end : float
        End time.
    eps : float
        Interface width parameter.
    dt : float
        Time step.
    save_every : int
        Save a snapshot every this many steps.
    f0 : ndarray of shape (N, N) or None
        Initial condition.  If None, a small-amplitude random field is used.
    seed : int
        RNG seed (controls both initial condition and all noise).
    obs_noise : float
        Observation noise level as a fraction of ``clean_snapshots.std()``.
        Gaussian noise with std = obs_noise * std(clean) is added to every
        saved snapshot independently.  0 means no observation noise.
    dyn_noise : float
        Absolute standard deviation sigma of the per-step Wiener increment
        in physical space.  0 means deterministic (no dynamical noise).

    Returns
    -------
    snapshots : ndarray of shape (n_saved, N, N)
    times : ndarray of shape (n_saved,)
    """
    rng = np.random.default_rng(seed)

    k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    K2 = KX**2 + KY**2

    denom = 1.0 + eps**2 * K2 * dt

    if f0 is None:
        f0 = 0.05 * rng.standard_normal((N, N))

    def nonlin(f):
        return f - f**3

    n_steps = int(round(T_end / dt)) + 1  # +1 to ensure we compute also time T_end
    snapshots, times = _run(f0, nonlin, denom, dt, n_steps, save_every, dyn_noise, rng)

    if obs_noise > 0.0:
        sigma = obs_noise * snapshots.std()
        snapshots = snapshots + sigma * rng.standard_normal(snapshots.shape)

    return snapshots, times
