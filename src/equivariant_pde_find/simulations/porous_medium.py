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
            times_list.append(step * dt)

    snapshots = np.array(snapshots_list)
    times = np.array(times_list)

    if obs_noise > 0.0:
        sigma_n = obs_noise * snapshots.std()
        snapshots = np.clip(snapshots + sigma_n * rng.standard_normal(snapshots.shape), 0.0, None)

    return snapshots, times
