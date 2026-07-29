"""Porous medium equation simulator (explicit Euler + spectral derivatives),
2-D and n-D variants."""

import numpy as np


def porous_medium(
    N=64,
    L=6 * np.pi,
    T_end=2.0,
    dt=0.002,
    save_every=25,
    f0=None,
    seed=42,
    obs_noise=0.0,
):
    """
    Simulate the porous medium equation (m=2) via explicit Euler + spectral derivatives.

    PDE:
        f_t = div(f^2 grad f) = f^2 Delta f + 2 f |grad f|^2

    Scheme:
        f^{n+1} = clip(f^n + dt * (f^2 * Delta f + 2 f |grad f|^2), 0)

    Spectral derivatives are used for accuracy on the periodic domain.
    The non-negativity constraint is enforced by clipping after each step.

    Parameters
    ----------
    N : int
        Grid points per direction.
    L : float
        Domain side length (periodic, [0, L]^2).
    T_end : float
        End time.
    dt : float
        Time step (explicit Euler; must be small for stability).
    save_every : int
        Save a snapshot every this many steps.
    f0 : ndarray of shape (N, N) or None
        Initial condition.  If None, a centred Gaussian of width L/8 is used.
    seed : int
        RNG seed.
    obs_noise : float
        Observation noise level as a fraction of ``clean_snapshots.std()``.

    Returns
    -------
    snapshots : ndarray of shape (n_saved, N, N)
    times : ndarray of shape (n_saved,)
    """
    rng = np.random.default_rng(seed)

    x = np.linspace(0, L, N, endpoint=False)
    k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    K2 = KX**2 + KY**2

    if f0 is None:
        X, Y = np.meshgrid(x, x, indexing="ij")
        cx, cy, sigma = L / 2, L / 2, L / 8
        f0 = np.exp(-((X - cx) ** 2 + (Y - cy) ** 2) / (2 * sigma**2))
        f0 /= f0.max()

    f = f0.copy()
    snapshots_list, times_list = [], []
    n_steps = int(round(T_end / dt)) + 1

    for step in range(n_steps):
        fhat = np.fft.fft2(f)
        Df = np.real(np.fft.ifft2(-K2 * fhat))
        fx = np.real(np.fft.ifft2(1j * KX * fhat))
        fy = np.real(np.fft.ifft2(1j * KY * fhat))
        rhs = f**2 * Df + 2.0 * f * (fx**2 + fy**2)
        f = np.clip(f + dt * rhs, 0.0, None)
        if step % save_every == 0:
            snapshots_list.append(f.copy())
            times_list.append((step + 1) * dt)

    snapshots = np.array(snapshots_list)
    times = np.array(times_list)

    if obs_noise > 0.0:
        sigma_n = obs_noise * snapshots.std()
        snapshots = np.clip(snapshots + sigma_n * rng.standard_normal(snapshots.shape), 0.0, None)

    return snapshots, times


_MAX_NDIM = 8


def porous_medium_nd(
    ndim=2,
    N=None,
    L=6 * np.pi,
    T_end=2.0,
    dt=0.002,
    save_every=25,
    f0=None,
    seed=42,
    obs_noise=0.0,
):
    """
    Simulate the porous medium equation (m=2) on [0, L]^ndim via explicit
    Euler + spectral derivatives.

    PDE:
        f_t = div(f^2 grad f) = f^2 Delta f + 2 f |grad f|^2

    Scheme:
        f^{n+1} = clip(f^n + dt * (f^2 * Delta f + 2 f |grad f|^2), 0)

    Parameters
    ----------
    ndim : int
        Number of spatial dimensions (1 <= ndim <= 8).
    N : int or None
        Grid points per direction.  Default: 64 for ndim=2, 48 for ndim>=3.
    L, T_end, dt, save_every, f0, seed, obs_noise :
        Same as ``porous_medium``.
    """
    if not 1 <= ndim <= _MAX_NDIM:
        raise ValueError(f"ndim must be 1-{_MAX_NDIM}, got {ndim}")

    if N is None:
        N = 64 if ndim == 2 else 48

    rng = np.random.default_rng(seed)

    x = np.linspace(0, L, N, endpoint=False)
    k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
    K_grids = np.meshgrid(*([k] * ndim), indexing="ij")
    K2 = sum(Ki**2 for Ki in K_grids)

    if f0 is None:
        X_grids = np.meshgrid(*([x] * ndim), indexing="ij")
        centre, sigma = L / 2, L / 8
        r2 = sum((Xi - centre) ** 2 for Xi in X_grids)
        f0 = np.exp(-r2 / (2 * sigma**2))
        f0 /= f0.max()

    f = f0.copy()
    snapshots_list, times_list = [], []
    n_steps = int(round(T_end / dt)) + 1

    for step in range(n_steps):
        fhat = np.fft.fftn(f)
        Df = np.real(np.fft.ifftn(-K2 * fhat))
        grad_sq = sum(
            np.real(np.fft.ifftn(1j * Ki * fhat)) ** 2 for Ki in K_grids
        )
        rhs = f**2 * Df + 2.0 * f * grad_sq
        f = np.clip(f + dt * rhs, 0.0, None)
        if step % save_every == 0:
            snapshots_list.append(f.copy())
            times_list.append((step + 1) * dt)

    snapshots = np.array(snapshots_list)
    times = np.array(times_list)

    if obs_noise > 0.0:
        sigma_n = obs_noise * snapshots.std()
        snapshots = np.clip(snapshots + sigma_n * rng.standard_normal(snapshots.shape), 0.0, None)

    return snapshots, times
