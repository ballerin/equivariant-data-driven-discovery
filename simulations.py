"""
Simulation utilities for generating the ground truth for the experiments.

Each simulator returns ``(snapshots, times)`` where ``snapshots`` has shape
``(n_saved, N, N)`` and ``times`` has shape ``(n_saved,)``.

Noise
-----
Two noise models are supported and can be combined:

observation noise (``obs_noise > 0``)
    Independent Gaussian noise added to every saved snapshot *after* the
    simulation:  ``snapshot += obs_noise * std(clean_snapshots) * randn``.
    Models imperfect sensors.  The underlying trajectory is deterministic;
    the PDE being discovered is unchanged.  Note that finite-difference
    estimates of f_t inherit O(obs_noise / dt_snap) amplification, so
    large ``obs_noise`` quickly destroys the regression target.

dynamical noise (``dyn_noise > 0``)
    Spatially white Wiener forcing injected at every time step (Euler-
    Maruyama).  The SPDE is  df = L(f)dt + dyn_noise * dW_t  where W_t
    is a cylindrical Wiener process on the grid.  The IMEX scheme becomes

        f̂^{n+1} = (f̂^n + dt * N̂(f^n) + sqrt(dt) * dyn_noise * ξ̂^n)
                   / denom(k)

    with ξ^n ~ N(0, I_{N²}).  The implicit denominator naturally damps
    high-frequency noise modes.  The *drift* (Allen-Cahn) is unchanged;
    regression is harder because f_t estimates absorb the stochastic term.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Internal IMEX kernel
# ---------------------------------------------------------------------------


def _imex_step(f, nonlin_fn, denom, dt, dyn_noise, rng):
    """Single IMEX step: implicit linear + explicit nonlinear + optional noise."""
    fhat = np.fft.fft2(f) + dt * np.fft.fft2(nonlin_fn(f))
    if dyn_noise > 0.0:
        xi = rng.standard_normal(f.shape)
        fhat += np.sqrt(dt) * dyn_noise * np.fft.fft2(xi)
    return np.real(np.fft.ifft2(fhat / denom))


def _run(f0, nonlin_fn, denom, dt, n_steps, save_every, dyn_noise, rng):
    """Run the time loop, returning (snapshots, step_indices)."""
    f = f0.copy()
    snapshots, steps = [], []
    for step in range(n_steps):
        f = _imex_step(f, nonlin_fn, denom, dt, dyn_noise, rng)
        if step % save_every == 0:
            snapshots.append(f.copy())
            steps.append(step)
    return np.array(snapshots), np.array(steps, dtype=float) * dt


# ---------------------------------------------------------------------------
# Allen-Cahn
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Swift-Hohenberg
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Kuramoto–Sivashinsky
# ---------------------------------------------------------------------------


def kuramoto_sivashinsky(
    N=64,
    L=16 * np.pi,
    T_end=400.0,
    dt=0.05,
    save_every=5,
    f0=None,
    seed=42,
    obs_noise=0.0,
    dyn_noise=0.0,
):
    """Simulate the 2D Kuramoto–Sivashinsky-type dynamics via spectral IMEX.

    We integrate the form used in the notebooks:
        f_t = -1/2 |∇f|^2 - Δf - Δ^2 f

    The linear Fourier symbol is L_k = K^2 - K^4 and the nonlinear term is
    the gradient-square contribution computed spectrally.
    """
    rng = np.random.default_rng(seed)

    k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    K2 = KX**2 + KY**2
    K4 = K2**2

    Lk = K2 - K4
    denom = 1.0 - dt * Lk

    if f0 is None:
        f0 = 0.1 * rng.standard_normal((N, N))

    def nonlin(f):
        # compute spectral gradients then form -1/2*(fx^2 + fy^2)
        fhat = np.fft.fft2(f)
        fx = np.real(np.fft.ifft2(1j * KX * fhat))
        fy = np.real(np.fft.ifft2(1j * KY * fhat))
        return -0.5 * (fx * fx + fy * fy)

    n_steps = int(round(T_end / dt)) + 1
    snapshots, times = _run(f0, nonlin, denom, dt, n_steps, save_every, dyn_noise, rng)

    if obs_noise > 0.0:
        sigma = obs_noise * snapshots.std()
        snapshots = snapshots + sigma * rng.standard_normal(snapshots.shape)

    return snapshots, times


# ---------------------------------------------------------------------------
# Cahn-Hilliard
# ---------------------------------------------------------------------------


def cahn_hilliard(
    N=64,
    L=2 * np.pi,
    T_end=20.0,
    eps=0.3,
    dt=0.01,
    save_every=5,
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


# ---------------------------------------------------------------------------
# Porous Medium Equation (PME, m=2)
# ---------------------------------------------------------------------------


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
            times_list.append(step * dt)

    snapshots = np.array(snapshots_list)
    times = np.array(times_list)

    if obs_noise > 0.0:
        sigma_n = obs_noise * snapshots.std()
        snapshots = np.clip(snapshots + sigma_n * rng.standard_normal(snapshots.shape), 0.0, None)

    return snapshots, times


# ---------------------------------------------------------------------------
# Helpers for regression dataset preparation
# ---------------------------------------------------------------------------


def make_regression_data(snapshots, times, n_per_snap=200, skip_fraction=0.0, seed=0):
    """
    Extract (f_snaps, ft_snaps, sample_indices) from a snapshot array.

    Uses central finite differences in time to estimate f_t.  Optionally
    skips an initial transient fraction of the trajectory.

    Parameters
    ----------
    snapshots : ndarray of shape (n_saved, N, N)
    times : ndarray of shape (n_saved,)
    n_per_snap : int
        Number of spatial points sampled per snapshot.
    skip_fraction : float in [0, 1)
        Skip this fraction of snapshots from the start (transient removal).
    seed : int
        RNG seed for spatial subsampling.

    Returns
    -------
    f_snaps : ndarray of shape (n_interior, N, N)
        Field snapshots at interior time points.
    ft_snaps : ndarray of shape (n_interior, N, N)
        Time derivatives estimated by central finite difference.
    idx_sets : list of ndarray
        Spatial indices sampled for each snapshot.
    dt_snap : float
        Time between consecutive saved snapshots.
    """
    skip = max(0, int(skip_fraction * len(snapshots)))
    snapshots = snapshots[skip:]
    times = times[skip:]

    dt_snap = times[1] - times[0]
    ft_snaps = (snapshots[2:] - snapshots[:-2]) / (2.0 * dt_snap)
    f_snaps = snapshots[1:-1]

    N2 = f_snaps.shape[1] * f_snaps.shape[2]
    rng = np.random.default_rng(seed)
    idx_sets = [rng.choice(N2, n_per_snap, replace=False) for _ in range(len(f_snaps))]

    return f_snaps, ft_snaps, idx_sets, dt_snap
